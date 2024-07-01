args<-commandArgs(TRUE)

Data.File <- args[1]
SoftPow <- as.numeric(args[2])
Block.Size <- as.numeric(args[3])
min.Module.Size <- as.numeric(args[4])
Save.TOM <- args[5]
Plot.Dendro <- args[6]
OutPrefix <- args[7]

sink(paste0(OutPrefix,".WGCNA.BlocwiseNet.log.txt"))

cat("Input arguments:\n")
cat("    Input data file:",Data.File,"\n")
cat("    Soft Power threshold:",SoftPow,"\n")
cat("    Maximum block size:",Block.Size,"\n")
cat("    Minimum module size:",min.Module.Size,"\n")
cat("    Do you want to save TOM data?",Save.TOM,"\n")
cat("    Do you want to generate dendrogram for each block?",Plot.Dendro,"\n")
cat("    Output prefix:",OutPrefix,"\n")
cat("\n")

cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

#########################################

Save.TOM <- ifelse(trimws(tolower(Save.TOM))=="yes",T,F)
Plot.Dendro <- ifelse(trimws(tolower(Plot.Dendro))=="yes",T,F)

Data <- readRDS(Data.File)

cat("Generating network...\n")

nthr = max(1, parallel::detectCores(), na.rm = TRUE)
net = blockwiseModules(datExpr = t(Data), maxBlockSize = Block.Size, power = SoftPow, TOMType = "unsigned", 
                       minModuleSize = min.Module.Size, reassignThreshold = 0, mergeCutHeight = 0.25, 
                       numericLabels = F, saveTOMs= Save.TOM, verbose = 3, nThreads = nthr)

cat("Saving the network object...\n")
saveRDS(net,file = paste0(OutPrefix,".WGCNA.Net.rds"))

if(Plot.Dendro){
  cat("Generating plots...\n")
  sizeGrWindow(6,6)
  pdf(paste0(OutPrefix,".WGCNA.Net.Dendrogram.pdf"))
  for (i in 1:length(net$dendrograms)) {
    plotDendroAndColors(net$dendrograms[[i]], net$blockGenes[[i]],
                        "Module colors", main = "Gene dendrogram and module colors in block 1", 
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }
  dev.off()
}

print("All done!")


sink()
