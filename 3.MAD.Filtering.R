args <- commandArgs(TRUE)

beta_file <- args[1] 
mad_thr <- as.numeric(args[2])
out_pref <- args[3]

sink(paste0(out_pref,".MAD.",mad_thr,".log.txt"))

cat("Input arguments:\n")
cat("    Beta value file:",beta_file,"\n")
cat("    MAD Threshold:",mad_thr,"\n")
cat("    Output prefix:",out_pref,"\n")
cat("\n")


cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dendextend))

cat("Reading beta values...\n")
betas <- readRDS(file=beta_file)

cat("Keeping top ",(mad_thr*100),"% of CpGs with high Median absolute deviation\n")
## calculate gene Median Absolute Deviation across all samples
gene_mad <- apply(betas,1,mad)

gene_mad <- sort(gene_mad,decreasing = T)
gene_mad1 <- gene_mad[1:round(length(gene_mad)*mad_thr,digits = 0)]

cat("Total number of CpGs: ",length(gene_mad),"\n")
cat("Number of CpGs after filtering: ",length(gene_mad1),"\n")
## removing genes with low MAD
betas <- betas[rownames(betas) %in% rownames(as.data.frame(gene_mad1)),]

## checking genes and samples with too many missing value
passed_gene_samples <- goodSamplesGenes(betas,verbose = 3)
cat("goodSamplesGenes function: All samples and genes passed? ",ifelse(passed_gene_samples$allOK,"Yes","NO"))
cat("\n")

cat("Saving filtered data...\n")
saveRDS(betas, file = paste0(out_pref,".MAD.",mad_thr,".rds"))

cat("All done!")
cat("\n")

sink()

