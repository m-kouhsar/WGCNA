args <- commandArgs(TRUE)

beta_file <- trimws(args[1] )
meth_filter <- trimws(args[2]) # set it to y/yes/T/true/1 if you need to filter the SNP,Sex and Cross Hybridising probes
convert2M <- trimws(args[3])   # set it to y/yes/T/true/1 if you need to convert Beta values to M values
mad_thr <- as.numeric(trimws(args[4]))
out_prefix <- trimws(args[5])
ScriptDir <- trimws(args[6])


cat("Input arguments:\n")
cat("    Methylation data file (rds format):",beta_file,"\n")
cat("    Filtering methylation data (Sex, SNP and Cross Hybridising probes)?",meth_filter,"\n")
cat("    Convert Beta values to M values?",convert2M,"\n")
cat("    MAD Threshold:",mad_thr,"\n")
cat("    Output prefix:",out_prefix,"\n")
cat("    Scripts directory:",ScriptDir,"\n")
cat("\n")

##############################################################################

cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(lumi))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
cat("\n")
###############################################################################

dir.create(path = dirname(out_prefix) , recursive = T , showWarnings = F)
cat("Reading beta values...\n")
betas <- readRDS(file=beta_file)

message("Total number of probes:",dim(betas)[1],"\nTotal number of samples:",dim(betas)[2])
meth_filter <- tolower(meth_filter)
if((meth_filter == "y")|(meth_filter == "yes")|(meth_filter == "t")|(meth_filter == "true")|(meth_filter == "1")){
  message("Reading Cross Hydridising Probes data..")
  crosslist<-read.table(paste0(ScriptDir , "/References/CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)[,1]
  
  message("Reading SNP Probes data..")
  snpProbes<-read.table(paste0(ScriptDir , "/References/SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
  SNPlist<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),1]
  
  message("Reading manifest file...")
  manifest<-fread(paste0(ScriptDir , "/References/EPIC.V1.Manifest.tsv"), data.table = F)
  
  message("**********************************************************************\n")
  
  message("Step 1: Removing Cross Hybridising Probes ...")
  index <- rownames(betas) %in% crosslist
  betas <-betas[!index,]
  message("        ",sum(index)," probes are Cross Hybridising.")
  
  message("Step 2: Removing SNP probes...")
  index <- rownames(betas) %in% SNPlist
  betas <- betas[!index,]
  message("        ",sum(index)," probes are SNPs.")
  
  message("Step 3: Removing probes in the sex chromosomes...")
  sexchrX<-manifest[which(manifest$CHR=="X"),]
  sexchrY<- manifest[which(manifest$CHR== "Y"),]
  index<- (rownames(betas) %in% sexchrX$IlmnID)|(rownames(betas) %in% sexchrY$IlmnID)
  betas <- betas[!index,]
  message("        ",sum(index)," are in the sex chromosomes.")
  
  message("Step 4: Removing probes with rs name...")
  index <- substr(rownames(betas),1,2)=="rs"
  betas <- betas[!index,]
  message("        ",sum(index)," probes start with 'rs'")
  
  message("Final number of probes after filtering: ",nrow(betas))
  
}

######################################################################
message("Keeping top ",(mad_thr*100),"% of CpGs with high Median Absolute Deviation")
## calculate gene Median Absolute Deviation across all samples
gene_mad <- apply(betas,1,mad)

message("MAD Summary:")
print(summary(gene_mad))

pdf(file = paste0(out_prefix , ".MAD.Hist.pdf") , height = 8 , width = 8)
ggplot(data.frame(MAD = gene_mad), aes(x = MAD)) + 
  # Histogram: Note the y = after_stat(density) argument
  geom_histogram(aes(y = after_stat(density)),
                 fill = "lightblue", 
                 color = "black", bins = 30) +
  # Density line: Overlay the line
  geom_density(color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Histogram of Median Absolote Deviation", x = "MAD", y = "Density")
graphics.off()

gene_mad <- sort(gene_mad,decreasing = T)
gene_mad1 <- gene_mad[1:round(length(gene_mad)*mad_thr,digits = 0)]

message("Total number of CpGs: ",length(gene_mad))
message("Number of CpGs after filtering: ",length(gene_mad1))
## removing genes with low MAD
betas <- betas[rownames(betas) %in% rownames(as.data.frame(gene_mad1)),]

convert2M <- tolower(convert2M)
if((convert2M == "y")|(convert2M == "yes")|(convert2M == "t")|(convert2M == "true")|(convert2M == "1")){
  message("Converting Beta values to M values...")
  betas <- lumi::beta2m(beta = betas)
  out_prefix <- paste0(out_prefix,".Mval")
}
message("Saving filtered data...")
saveRDS(betas, file = paste0(out_prefix,".Filtered.MAD.",(mad_thr*100),".rds"))

message("All done!")
cat("\n")
