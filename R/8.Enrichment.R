argument <- commandArgs(T)

MM_GS_file <- trimws(argument[1])
Net_file <- trimws(argument[2])
ID_type <- trimws(argument[3])                #entrez,symbol,cpg,ensembl
MM <- as.numeric(trimws(argument[4]))
MM_pval <- as.numeric(trimws(argument[5]))
GS <- as.numeric(trimws(argument[6]))
GS_pval <- as.numeric(trimws(argument[7]))
out_prefix <- trimws(argument[8])

############################################################################
message("Input arguments:")
message("      Module Membesrhip and Gene Significance file: " , MM_GS_file)
message("      WGCNA network file: " , Net_file)
message("      Genes/Probes ID type (entrez, symbol, ensembl or cpg): " , ID_type)
message("      Module Membesrhip threshold for selecting hub genes/probes: " , MM)
message("      Module Membesrhip P-value threshold for selecting hub genes/probes: " , MM_pval)
message("      Gene significance threshold for selecting hub genes/probes: " , GS)
message("      Gene significance P-value threshold threshold for selecting hub genes/probes: " , GS_pval)
message("      Output files prefix: " , out_prefix)
#############################################################################
message("")

message("Loading required packages...")

if(!(ID_type %in% c("entrz","symbol","cpg","enseble"))){
  stop("ID type must be one of the following values:\n entrez,symbol,cpg,ensembl")
}
methylation <- ifelse(tolower(ID_type) == "cpg" , T , F)

if(methylation){
  suppressPackageStartupMessages(library(missMethyl))
}

if(!methylation){
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
  })
}
suppressPackageStartupMessages({
  library(funr)
  library(ggplot2)  
  library(dplyr)
})

source(paste0(dirname(sys.script()),"/Enrichment.Function.R"))
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
  message("CpG set enrichment analysis using missMethyl::gometh...")
  
  message("      GO enrichment...")
  go_results <- gometh(
    sig.cpg = gene_list,          # Vector of significant CpG IDs
    all.cpg = univers_list,          # Vector of background CpGs (the universe)
    collection = "GO",           # Specify "GO" or "KEGG"
    array.type = "EPIC",         # Change to "450K" if using older arrays
    prior.prob = TRUE  ,          # Adjusts for the varying number of CpGs per gene
    sig.genes = T
  )
  
  message("      KEGG enrichment...")
  kegg_results <- gometh(
    sig.cpg = gene_list,
    all.cpg = univers_list,
    collection = "KEGG",
    array.type = "EPIC",
    prior.prob = TRUE,
    sig.genes = T
  )
  
  go_results$GeneRatio <- go_results$DE/length(gene_list)
  go_results$BgRatio <- go_results$N/length(univers_list)
  go_results_BP <- go_results[(!is.na(go_results$ONTOLOGY)) & (go_results$ONTOLOGY == "BP"),]
  go_results_MF <- go_results[(!is.na(go_results$ONTOLOGY)) & (go_results$ONTOLOGY == "MF"),]
  go_results_CC <- go_results[(!is.na(go_results$ONTOLOGY)) & (go_results$ONTOLOGY == "CC"),]
  
  message("      Generating plots...")
  plot.BP <- gometh_dotplot(gometh_res = go_results_BP , showCategory = 20,
                            plot.title = "Top 20 GO Biological Process enrichment results" )
  plot.MF <- gometh_dotplot(gometh_res = go_results_MF , showCategory = 20 ,
                            plot.title = "Top 20 GO Molecular Function enrichment results")
  plot.CC <- gometh_dotplot(gometh_res = go_results_CC , showCategory = 20,
                            plot.title = "Top 20 GO Cellular Component enrichment results" )
  
  kegg_results$GeneRatio <- kegg_results$DE/length(gene_list)
  kegg_results$BgRatio <- kegg_results$N/length(univers_list)
  
  plot.KEGG <- gometh_dotplot(gometh_res = kegg_results , showCategory = 20,
                              plot.title = "Top 20 KEGG enrichment results")
  
}

if(!methylation){
  message("Gene set enrichment analysis using clusterprofiler...")
  
  if(ID_type == "symbol"){
    gene_list_entrez <- clusterProfiler::bitr(
      gene_list, 
      fromType = "SYMBOL", 
      toType   = "ENTREZID", 
      OrgDb    = org.Hs.eg.db)
    }else if(ID_type == "ensembl"){
      gene_list_entrez <- clusterProfiler::bitr(
        gene_list, 
        fromType = "ENSEMBL",    # Current format
        toType   = "ENTREZID",   # Desired format
        OrgDb    = org.Hs.eg.db
      )
    }else{
      gene_list_entrez <- gene_list
    }
    
    message("      GO enrichment...")
    go_results <- enrichGO(
      gene          = gene_list_entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "ALL",       # "BP" (Biological Process), "MF", "CC", or "ALL"
      pAdjustMethod = "BH", 
      universe = univers_list,     # Benjamini-Hochberg FDR correction
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable      = TRUE,         # Converts Entrez IDs back to symbols in the output table
      pool = F
    )
    
    message("      KEGG enrichment...")
    kegg_results <- enrichKEGG(
      gene          = gene_list_entrez,
      organism      = "hsa",       # "hsa" stands for Homo sapiens
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = univers_list,
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = F   #FALSE: use latest online db. TRUE: use KEGG.db package
    )
    
    go_results_BP <- go_results |> filter(ONTOLOGY == "BP")
    go_results_MF <- go_results |> filter(ONTOLOGY == "MF")
    go_results_CC <- go_results |> filter(ONTOLOGY == "CC")
    
    message("Generating plots...")
    plot.BP <- dotplot(go_results_BP , 
                       x = "GeneRatio",
                       color = "p.adjust",
                       showCategory = 20, 
                       title="Top 20 GO Biological Process enrichment results"
                      )
    plot.MF <- dotplot(go_results_MF , 
                       x = "GeneRatio",
                       color = "p.adjust",
                       showCategory = 20, 
                       title="Top 20 GO Molecular Function enrichment results"
                      )
    plot.CC <- dotplot(go_results_CC , 
                       x = "GeneRatio",
                       color = "p.adjust",
                       showCategory = 20, 
                       title="Top 20 GO Cellular Component enrichment results"
                      )
    
    go_results_BP <- as.data.frame(go_results_BP@result)
    go_results_MF <- as.data.frame(go_results_MF@result)
    go_results_CC <- as.data.frame(go_results_CC@result)
    
    plot.KEGG <- dotplot(kegg_results , 
                         x = "GeneRatio",
                         color = "p.adjust",
                         showCategory = 20, 
                         title="Top 20 KEGG enrichment results"
                         )
    
    kegg_results <- as.data.frame(kegg_results@result)
    
}

message("Saving results...")
pdf(paste0(out_prefix , ".Enrichment.pdf"))
print(plot.KEGG)
print(plot.BP)
print(plot.MF)
print(plot.CC)
graphics.off()

write.csv(kegg_results , file = paste0(out_prefix , "Enrichment.KEGG.csv"))
write.csv(go_results_BP , file = paste0(out_prefix , "Enrichment.BP.csv"))
write.csv(go_results_MF , file = paste0(out_prefix , "Enrichment.MF.csv"))
write.csv(go_results_CC , file = paste0(out_prefix , "Enrichment.CC.csv"))



