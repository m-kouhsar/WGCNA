args <- commandArgs(TRUE)

counts.file <- trimws(args[1])
pheno.file <- trimws(args[2])
var.trait <- trimws(args[3])
normalize.method <- trimws(args[4])
OutPrefix <- trimws(args[5])


message("Input arguments:")
message("        Count matrix: ", counts.file)
message("        Phenotype file: ", pheno.file)
message("        Trait variable: ", var.trait)
message("        Normalization method (cpm or vst): ", normalize.method)
message("        Output files prefix: ", OutPrefix)
cat("\n")

##############################################################################

cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(edgeR))

###############################################################################

dir.create(path = dirname(OutPrefix) , recursive = T , showWarnings = F)
cat("Reading expression data...\n")
counts <- read.table(file = counts.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

if(!identical(colnames(counts) , rownames(pheno))){
  message("Warning message:\nColnames in the count matrix are not equal to the rownames in the phenotype file!\nShared names will be considered.")
  index <- intersect(colnames(counts), rownames(pheno))
  message("Number of shared names: " , length(index))
  if(length(index) < 1){
    stop("There is no shared ID between Phneotype and Count data!")
  }
  counts <- counts[,index]
  pheno <- pheno[index , ]
}
pheno[,var.trait] <- as.factor(pheno[,var.trait])
message("Total number of probes:",dim(counts)[1],"\nTotal number of samples:",dim(counts)[2])

######################################################################
message("filtering low count genes...")
keep <- edgeR::filterByExpr(counts,group = pheno[,var.trait])
message(sum(!keep),"/",nrow(counts)," genes removed. Remaining genes:", sum(keep))
counts <- counts[keep,]

if(normalize.method == "cpm"){
  message("Normalizing gene counts uisng cpm function in edgeR...")
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)
  counts <- edgeR::cpm(dge , log = TRUE)
  
}else{
  if(normalize.method == "vst"){
    message("Normalizing gene counts uisng vst method in DESeq2...")
    dds <- DESeqDataSetFromMatrix(countData = counts , 
                                  colData = data.frame(row.names = colnames(counts), group = rep(1, ncol(counts))),
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    dds_vst <- DESeq2::vst(dds , blind = T)
    counts <- SummarizedExperiment::assay(dds_vst)
  }
}

cat("Saving filtered data...\n")
write.table(counts , file = paste0(OutPrefix,".Filtered.tsv") , row.names = T , col.names = T , sep = "\t" , quote = F)

cat("All done!")
cat("\n")
