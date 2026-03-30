
args = commandArgs(T)

expr_file <- trimws(args[1])  #Methylation matrix in rds format or expression matrix in tsv format
pheno_file <- trimws(args[2])
outliers <- trimws(str_split_1(trimws(args[3])))
out_prefix <- trimws(args[4])

message("Input arguments:")
message("        Methylatio/Expression data file: ", expr_file)
message("        Phenotype file: ", pheno_file)
message("        Outlier samples: ", outliers)
message("        Output files prefix: ", out_prefix)
cat("\n")

##################################################################################################
dir.create(dirname(OutPrefix) , recursive = T , showWarnings = F)

message("Reading the data ...")
if(str_ends(string = expr.file , pattern = ".rds")){
  expr_mat <- readRDS(expr.file)
}else{
  expr_mat <- read.table(file = counts.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

if(!identical(colnames(expr_mat) , rownames(pheno))){
  message("Warning message:\nColnames in the count matrix are not equal to the rownames in the phenotype file!\nShared names will be considered.")
  index <- intersect(colnames(expr_mat), rownames(pheno))
  message("Number of shared names: " , length(index))
  if(length(index) < 1){
    stop("There is no shared ID between Phneotype and Count data!")
  }
  expr_mat <- expr_mat[,index]
  pheno <- pheno[index , ]
}

####################################################################################

if(all(outliers != "")){
  message("Removing outliers...")
  if(! all (outliers %in% colnames(expr_mat))){
    warning("The followin outlier IDs are not exist in the data: ")
    paste(outliers[!(outliers %in% colnames(expr_mat))] , collapse = ";")
  }
  expr_mat <- expr_mat[,!(colnames(expr_mat) %in% outliers)]
  pheno <- pheno[!(rownames(pheno) %in% outliers),]
}

####################################################################################
message("Saving the results...")
if(str_ends(string = expr.file , pattern = ".rds")){
  saveRDS(expr_mat, file = paste0(out_prefix,".noOutlr.rds"))
}else{
  write.table(expr_mat , file = paste0(out_prefix , ".noOutlr.tsv") , row.names = T , col.names = T , sep = "\t" , quote = F)
}



