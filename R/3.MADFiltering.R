
args <- commandArgs(TRUE)

expr.file <- trimws(args[1] )
mad_thr <- as.numeric(trimws(args[2] ))
out_prefix <- trimws(args[3] )

message("Input arguments:")
message("    Methylation/Expression data file: ",expr.file)
message("    Median Absolute Deviation (MAD) threshold: ",mad_thr)
message("    Output prefix: ",out_prefix)
cat("\n")

########################################################################################################################################
message("Loading libraries...")
library(stringr)
suppressWarnings(suppressMessages(library(ggplot2)))
########################################################################################################################################
message("Reading input data...")
if(str_ends(string = expr.file , pattern = ".rds")){
  expr_mat <- readRDS(expr.file)
}else{
  expr_mat <- read.table(file = expr.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

message("Keeping top ",(mad_thr*100),"% of CpGs/genes with highest Median Absolute Deviation")
## calculate gene Median Absolute Deviation across all samples
gene_mad <- apply(expr_mat,1,mad)

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

message("Total number of CpGs/genes: ",length(gene_mad))
message("Number of CpGs/genes after filtering: ",length(gene_mad1))
## removing genes with low MAD
expr_mat <- expr_mat[rownames(expr_mat) %in% rownames(as.data.frame(gene_mad1)),]

message("Saving filtered data...")
if(str_ends(string = expr.file , pattern = ".rds")){
  saveRDS(expr_mat, file = paste0(out_prefix,".MAD.",(mad_thr*100),".rds"))
}else{
  write.table(expr_mat , file = paste0(out_prefix,".MAD.",(mad_thr*100),".tsv") , row.names = T , col.names = T , sep = "\t" , quote = F)
}
message("All done!")

