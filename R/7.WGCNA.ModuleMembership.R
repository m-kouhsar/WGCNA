arguments <- commandArgs(T)

net_file <- trimws(arguments[1])
expr_file <- trimws(arguments[2])
pheno_file <- trimws(arguments[3])
trait <- trimws(arguments[4])
categorical_trait <- trimws(arguments[5])
cofounders_num <- trimws(arguments[6])
cofounders_cat <- trimws(arguments[7])
GS_legend_pvalue <- as.numeric(trimws(arguments[8]))
GS_label_pvalue <- as.numeric(trimws(arguments[9]))
adjusted_pvalue_method <- trimws(arguments[10])
modules <- trimws(arguments[11])
softPow <- as.numeric(trimws(arguments[12]))
out_prefix <- trimws(arguments[13])

message("Input arguments")
message("      Network file: ",net_file)
message("      Expression file: ",expr_file)
message("      Metadata file: ",pheno_file)
message("      Trait variable: ",trait)
message("      Is trait categorical? ",categorical_trait)
message("      Categorical cofounders: ",cofounders_cat)
message("      Numeric cofounders: ",cofounders_num)
message("      Gene Significance P-value threshold used for coloring genes in the scatter plot: ",GS_legend_pvalue)
message("      Gene Significance P-value threshold used for labeling genes in the scatter plot: ",GS_label_pvalue)
message("      Adjusted P-value method: ",adjusted_pvalue_method)
message("      Analysis modules: ",modules)
message("      Soft threshold power: ",softPow)
message("      Output files prefix: ",out_prefix)
message("#####################################################################################################")
message("")
message("Loading required packages...")
suppressPackageStartupMessages({
  library(WGCNA)
  library(stringr)
  library(funr)
})
message("")
source(paste0(dirname(sys.script()),"/WGCNA.ModuleMembership.Function.R"))
##################################################################################################

message("Reading inputs...")

net <- readRDS(net_file)
pheno = read.csv(pheno_file,row.names = 1,stringsAsFactors = F)
if(str_ends(expr_file , pattern = ".rds")){
  expr = readRDS(expr_file)
}else{
  expr = read.table(file = expr_file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

if(!identical(names(net$colors),rownames(expr))){
  warning("CpG/gene IDs in network and beta files are not match. It will be matched using shared CpGs...")
  index = intersect(names(net$colors),rownames(expr))
  expr = expr[index , ]
  net$colors = net$colors[index]
}
if(!identical(colnames(expr) , rownames(pheno))){
  warning("Sample IDs in phenotype and beta files are not match. It will be matched using shared samples...")
  index = intersect(colnames(expr) , rownames(pheno))
  expr = expr[,index]
  pheno = pheno[index ,]
}

if(!(trait %in% colnames(pheno))){
  stop("Cannot find trait variable (",trait , ") in metadata!")
}

categorical_trait <- ifelse(tolower(categorical_trait) == "yes" , T , F)

cofounders_num <- trimws(str_split_1(cofounders_num , pattern = ","))
cofounders_cat <- trimws(str_split_1(cofounders_cat , pattern = ","))



if(tolower(modules=="all")){
  modules = unique(net$colors)
}else{
  modules <- trimws(str_split_1(modules , pattern = ","))
}
message("")
#################################################################
message("Running the analysis...")
message("")
out_dir <- dirname(out_prefix)
if(!dir.exists(out_dir)){
  dir.create(out_dir , recursive = T)
}

for(v in cofounders_cat){
  pheno[,v] <- as.factor(pheno[,v])
}

for(v in cofounders_num){
  pheno[,v] <- as.numeric(pheno[,v])
}

MM <- Module.Membership.Plot(net.colors = net$colors , expr.mat = expr , metadata = pheno , trait = trait , 
                               categorical.trait = categorical_trait, cofounders = c(cofounders_num , cofounders_cat),
                             modules = modules,legend_GS_pval = GS_legend_pvalue , label_GS_pval = GS_label_pvalue , 
                             p.adjust.method = adjusted_pvalue_method , soft.power = softPow)
message("")
message("Saving the results...")
for (m in modules) {
  tiff(filename = paste0(out_prefix,".",m , ".MM_vs_GS.Plot.tif") , width = 8 , height = 8 , units = "in" , res = 300)
  print(MM$Plots[[m]])
  graphics.off()
  write.csv(MM$Data[[m]] , file = paste0(out_prefix,".",m,".MM.GS.csv"),row.names = F)
}
write.csv(MM$Cor.test , file =paste0(out_prefix,".","MM_vs_GS.csv") )




