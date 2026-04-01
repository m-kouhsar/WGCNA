
args <- commandArgs(TRUE)

expr.file <- trimws(args[1] )
mad_thr <- as.numeric(trimws(args[2] ))
out_prefix <- trimws(args[3] )

cat("Input arguments:\n")
cat("    Methylation/Expression data file:",expr.file,"\n")
cat("    Median Absolute Deviation (MAD) threshold:",mad_thr,"\n")
cat("    Output prefix:",out_prefix,"\n")
cat("\n")

cat("Reading input data...\n")
if(str_ends(string = expr.file , pattern = ".rds")){
  expr_mat <- readRDS(expr.file)
}else{
  expr_mat <- read.table(file = expr.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

message("Keeping top ",(mad_thr*100),"% of CpGs with high Median Absolute Deviation")
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

message("Total number of CpGs: ",length(gene_mad))
message("Number of CpGs after filtering: ",length(gene_mad1))
## removing genes with low MAD
expr_mat <- expr_mat[rownames(expr_mat) %in% rownames(as.data.frame(gene_mad1)),]

cat("Saving regressed data...\n")
if(str_ends(string = expr.file , pattern = ".rds")){
  saveRDS(expr_mat, file = paste0(out_prefix,".MAD.",mad_thr,".rds"))
}else{
  write.table(expr_mat , file = paste0(out_prefix,".MAD.",mad_thr,".tsv") , row.names = T , col.names = T , sep = "\t" , quote = F)
}
cat("All done!\n")

