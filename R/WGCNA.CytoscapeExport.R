args<-commandArgs(TRUE)

net_file <- trimws(args[1])
module_name <- trimws(args[2])
weighted_ <- trimws(ifelse(args[3]=="T",T,F))
threshold_ <- trimws(as.numeric(args[4]))
out_pref <- trimws(args[5])
TOM.Dir <- trimws(args[6])

print("Input arguments:")
print(paste("    Network file:",net_file))
print(paste("    Module name:",module_name))
print(paste("    Is weighted?",weighted_))
print(paste("    Adjacency threshold:",threshold_))
print(paste("    Output prefix:",out_pref))
cat("\n")
print("Loading libraries...")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

#########################################
if(is.na(TOM.Dir) | (TOM.Dir=="")){
  TOM.Dir <- dirname(net_file)
}
message("Openning network file...")
net = readRDS(net_file)

message("Finding Module nodes and blocks...")
module.nodes.index = which(net$colors==module_name)
module.blocks = net$blocks[module.nodes.index]

message("Module are in the following blocks: ",unique(module.blocks))

for (b in unique(module.blocks)) {
  
  message("Loading block ",b,"...")
  load(paste0(dirname(net_file),"/",net$TOMFiles[b]))
  
  message("Convert block ",b," To matrix ...")
  TOM.mat = as.matrix(TOM)
  module.nodes.TOM.index = net$blockGenes[[b]] %in% module.nodes.index
  module.TOM.mat=TOM.mat[module.nodes.TOM.index,module.nodes.TOM.index]
  
  message("Export Cytoscape data for block ",b,"...")
  cytoscape = exportNetworkToCytoscape(adjMat = module.TOM.mat,weighted = weighted_ , threshold = threshold_,
                                       edgeFile = paste0(out_pref,".",module_name,".Block.",b,".Thresh.",threshold_,".egdeData.txt"),
                                       nodeFile = paste0(out_pref,".",module_name,".Block.",b,".Thresh.",threshold_,".nodeData.txt"),
                                       nodeNames = names(module.nodes.index[module.blocks==b]))
}

