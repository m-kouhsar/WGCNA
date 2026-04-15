args <- commandArgs(TRUE)

beta_file <- trimws(args[1] )
convert2M <- trimws(args[2])   # set it to y/yes/T/true/1 if you need to convert Beta values to M values
out_prefix <- trimws(args[3])
ScriptDir <- trimws(args[4])


cat("Input arguments:\n")
cat("    Methylation data file (rds format):",beta_file,"\n")
cat("    Convert Beta values to M values?",convert2M,"\n")
# Decide between M value and Beta value: https://pubmed.ncbi.nlm.nih.gov/21118553/
cat("    Output prefix:",out_prefix,"\n")
cat("    Scripts directory:",ScriptDir,"\n")
cat("\n")

##############################################################################

cat("Loading libraries...\n")
suppressMessages(library(lumi))
suppressMessages(library(data.table))
cat("\n")
###############################################################################

dir.create(path = dirname(out_prefix) , recursive = T , showWarnings = F)
cat("Reading beta values...\n")
betas <- readRDS(file=beta_file)

message("Total number of probes:",dim(betas)[1],"\nTotal number of samples:",dim(betas)[2])

message("Reading Cross Hybridising Probes data..")
crosslist<-read.table(paste0(ScriptDir , "/References/CrossHybridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)[,1]

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

######################################################################

convert2M <- tolower(convert2M)
if((convert2M == "y")|(convert2M == "yes")|(convert2M == "t")|(convert2M == "true")|(convert2M == "1")){
  message("Converting Beta values to M values...")
  betas <- lumi::beta2m(beta = betas)
  out_prefix <- paste0(out_prefix,".Mval")
}
message("Saving filtered data...")
saveRDS(betas, file = paste0(out_prefix,".Filtered.rds"))

message("All done!")
cat("\n")
