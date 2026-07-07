argument <- commandArgs(T)

MM_GS_file <- trimws(argument[1])
Net_file <- trimws(argument[2])
ref_10x_dir <- trimws(argument[3])
sample_cellType_file <- trimws(argument[4])
ID_type <- trimws(argument[5])                #entrez,symbol,cpg,ensembl
MM <- as.numeric(trimws(argument[6]))
MM_pval <- as.numeric(trimws(argument[7]))
GS <- as.numeric(trimws(argument[8]))
GS_pval <- as.numeric(trimws(argument[9]))
out_prefix <- trimws(argument[10])

rep_ <- 1000
n_cores <- parallel::detectCores()
if(Sys.info()["sysname"] == "Windows"){
  n_cores <- 1 #in Windows EWCE does not support parallel computing
}
set.seed(123456)
dir.create(dirname(out_prefix))
############################################################################
message("Input arguments:")
message("      Module Membesrhip and Gene Significance file: " , MM_GS_file)
message("      WGCNA network file: " , Net_file)
message("      Reference single cell data directory: " , ref_10x_dir)
message("      Cell type information data file " , sample_cellType_file)
message("      Genes/Probes ID type (entrez, symbol, ensembl or cpg): " , ID_type)
message("      Module Membesrhip threshold for selecting hub genes/probes: " , MM)
message("      Module Membesrhip P-value threshold for selecting hub genes/probes: " , MM_pval)
message("      Gene significance threshold for selecting hub genes/probes: " , GS)
message("      Gene significance P-value threshold threshold for selecting hub genes/probes: " , GS_pval)
message("      Output files prefix: " , out_prefix)
#############################################################################
message("")

if(!(ID_type %in% c("entrz","symbol","cpg","enseble"))){
  stop("ID type must be one of the following values:\n entrez,symbol,cpg,ensembl")
}

message("Loading required packages...")


methylation <- ifelse(tolower(ID_type) == "cpg" , T , F)

suppressPackageStartupMessages({
  library(EWCE)
  library(Seurat)
  library(SummarizedExperiment)
  library(ggplot2)
  library(cowplot)
  library(parallel)
  
})
if(methylation){
  suppressPackageStartupMessages({
    library(minfi)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    })
}

if(!methylation){
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(clusterProfiler)
  })
}

dir.create(dirname(out_prefix) , recursive = T)

message("Reading inputs...")
MM_GS <- read.csv(MM_GS_file , stringsAsFactors = F)
net <- readRDS(Net_file)

gene_list <- MM_GS$ID[(MM_GS$MM > MM) & 
                        (MM_GS$GS > GS) & 
                        (MM_GS$MM.Pval < MM_pval) & 
                        (MM_GS$GS.Pval < GS_pval)]

if(length(gene_list) == 0){
  stop("No gene/probe passed the specified cutoffs.")
}else{
  message(length(gene_list)," genes/probes passed the specified cutoffs.")
}

univers_list <- names(net$colors)

if(methylation){
  message("Mapping CpGs to genes using Illumina annotation data (hg19)...")
  
  annot_all <- data.frame( minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  gene_symbol <- annot_all$UCSC_RefGene_Name[annot_all$Name %in% gene_list]
  gene_symbol <- gene_symbol[!is.na(gene_symbol)]
  gene_symbol <- gene_symbol[gene_symbol != ""]
  gene_symbol <- unique(unlist(strsplit(gene_symbol,";")))
  
  univers_symbol <- annot_all$UCSC_RefGene_Name[annot_all$Name %in% univers_list]
  univers_symbol <- univers_symbol[!is.na(univers_symbol)]
  univers_symbol <- univers_symbol[univers_symbol != ""]
  univers_symbol <- unique(unlist(strsplit(univers_symbol,";")))
  
}

if(!methylation){
  
  if(ID_type == "entrez"){
    
    message("Converting gene ENTREZ ids to gene symbol...")
    gene_symbol <- clusterProfiler::bitr(
      gene_list, 
      fromType = "ENTREZID", 
      toType   = "SYMBOL", 
      OrgDb    = org.Hs.eg.db)
    univers_symbol <- clusterProfiler::bitr(
      univers_list, 
      fromType = "ENTREZID", 
      toType   = "SYMBOL", 
      OrgDb    = org.Hs.eg.db)
    
    }else if(ID_type == "ensembl"){
      
      message("Converting gene ensembl ids to gene symbol...")
      gene_symbol <- clusterProfiler::bitr(
        gene_list, 
        fromType = "ENSEMBL",    # Current format
        toType   = "SYMBOL",   # Desired format
        OrgDb    = org.Hs.eg.db
      )
      univers_symbol <- clusterProfiler::bitr(
        univers_list, 
        fromType = "ENSEMBL",    # Current format
        toType   = "SYMBOL",   # Desired format
        OrgDb    = org.Hs.eg.db
      )
    }
    
}

message("Reading the single cell reference data...")

if(!file.exists(paste0(out_prefix,"_EWCERef.rda"))){
  
  sample_cellType <- read.csv(sample_cellType_file , stringsAsFactors = F)
  if(!all(c("broad.cell.type","Subcluster") %in% colnames(sample_cellType))){
    stop("sample_cellType_file must contains the following columns:
         'broad.cell.type'
         'Subcluster'")
  }
  sce <- Seurat::Read10X(ref_10x_dir, gene.column = 1)
  
  message("      [1]: Removing bad HGNC symbols in reference data...")
  sce1 <- fix_bad_hgnc_symbols(sce)
  
  message("      [2]: Droping uninformative genes...")
  sce2 <- drop_uninformative_genes(sce1, drop_nonhuman_genes=T , input_species = "human" ,output_species = "human" , 
                                   level2annot = sample_cellType$Subcluster, no_cores = n_cores)
  
  message("      [3]: Generating cell type dataset in approperiate format...\n")
  
  
  annotLevels = list(level1class=sample_cellType$broad.cell.type,
                     level2class= sample_cellType$Subcluster)
  ctd_file <- generate_celltype_data(exp=sce2,
                                     annotLevels=annotLevels,
                                     groupName="EWCERef",
                                     savePath=dirname(out_prefix), 
                                     file_prefix = basename(out_prefix),
                                     no_cores = n_cores) 
}

ctd <- EWCE::load_rdata(ctd_file)

message("Cell type enrichment analysis...")
result<- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = gene_symbol,
  bg = univers_symbol,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = rep_,
  annotLevel = 1, no_cores = n_cores)

p_ctd <- EWCE::plot_ctd(ctd = ctd , genes = gene_symbol[1:5],level = 1,metric = "specificity",show_plot = F)+
  ggtitle("EWCE Reference Dataset")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

p_result <- EWCE::ewce_plot(result$results)[[1]] +
  ggtitle("EWCE Cell Type Enrichment Result")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

message("Saving results...")
pdf(paste0(out_prefix , ".EWCE.pdf"))
print(p_ctd)
print(p_result)
graphics.off()

save(result , file = paste0(out_prefix , "EWCE.rdat"))
write.csv(result$results , file = paste(out_prefix , ".EWCE.csv") , row.names = F)


