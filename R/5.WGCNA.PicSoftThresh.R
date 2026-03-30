args<-commandArgs(TRUE)

Data.File <- trimws(args[1])
block_size <- as.numeric(trimws(args[2]))
out_prefix <- trimws(args[3])

cat("Input arguments:\n")
cat("    Methylation/Expression data file: ",Data.File,"\n")
cat("    Block size: ",block_size,"\n")
cat("    Output files prefix: ",out_prefix,"\n")
cat("\n")

cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

cat("Reading input data...")
if(str_ends(string = Data.File , pattern = ".rds")){
  Data <- readRDS(Data.File)
}else{
  Data <- read.table(file = Data.File, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

wgcna_input = t(Data)
gsg <- goodSamplesGenes(wgcna_input, verbose = 3)

message("goodSamplesGenes function in WGCNA:")
if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) 
    printFlush(paste("    Removing genes:", paste(colnames(wgcna_input)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) 
    printFlush(paste("    Removing samples:", paste(rownames(wgcna_input)[!gsg$goodSamples], collapse = ", ")))
  
  wgcna_input <- wgcna_input[gsg$goodSamples, gsg$goodGenes]
  print("Bad samples/genes removed. Data is now ready for WGCNA.")
} else {
  print("    All OK! No samples or genes need to be removed.")
}

cat("Calculating soft powers...\n")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(wgcna_input, powerVector = powers, verbose = 5, dataIsExpr = TRUE,blockSize = block_size)
saveRDS(sft,file = paste0(out_prefix,".WGCNA.SoftPower.rds"))

cat("Generating plots...\n")

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

pdf(file = paste0(out_prefix,".WGCNA.SoftPower.pdf"))

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

cat("All done!")

