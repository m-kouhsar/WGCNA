args<-commandArgs(TRUE)

Data.File <- trimws(args[1])
block_size <- as.numeric(trimws(args[2]))
out_prefix <- trimws(args[3])

message("Input arguments:")
message("    Methylation/Expression data file: ",Data.File)
message("    Block size: ",block_size)
message("    Output files prefix: ",out_prefix)
message("\n")

message("Loading libraries...")
suppressMessages(library(stringr))
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

message("Reading input data...")
if(str_ends(string = Data.File , pattern = ".rds")){
  Data <- readRDS(Data.File)
}else{
  Data <- read.table(file = Data.File, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

wgcna_input = t(Data)
message("goodSamplesGenes function in WGCNA:")
gsg <- goodSamplesGenes(wgcna_input, verbose = 3)
if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) 
    message(paste("    Removing genes:", paste(colnames(wgcna_input)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) 
    message(paste("    Removing samples:", paste(rownames(wgcna_input)[!gsg$goodSamples], collapse = ", ")))
  
  wgcna_input <- wgcna_input[gsg$goodSamples, gsg$goodGenes]
  
} else {
  message("    All OK! No samples or genes need to be removed.")
}

message("Calculating soft powers...")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(wgcna_input, powerVector = powers, verbose = 5, dataIsExpr = TRUE,blockSize = block_size)
saveRDS(sft,file = paste0(out_prefix,".WGCNA.SoftPower.rds"))

message("Generating plots...")


pdf(file = paste0(out_prefix,".WGCNA.SoftPower.pdf"))

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

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

graphics.off()

message("All done!")

