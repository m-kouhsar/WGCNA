args<-commandArgs(TRUE)

Data.File <- args[1]
block_size <- as.numeric(args[2])
out_pref <- args[3]

sink(paste0(out_pref,".WGCNA.SoftPower.log.txt"))

cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

Data <- readRDS(Data.File)

cat("Calculating soft powers...\n")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(t(Data), powerVector = powers, verbose = 5, dataIsExpr = TRUE,blockSize = block_size)
saveRDS(sft,file = paste0(out_pref,".WGCNA.SoftPower.rds"))

cat("Generating plots...\n")

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

pdf(file = paste0(out_pref,".WGCNA.SoftPower.pdf"))

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

sink()
