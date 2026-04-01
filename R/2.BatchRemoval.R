args <- commandArgs(TRUE)

expr.file <- trimws(args[1])
pheno.file <- trimws(args[2])
model_protect <- trimws(args[3]) 
model_remove <- trimws(args[4])
out_prefix <- trimws(args[5])

cat("Input arguments:\n")
cat("    Input expression file: ",expr.file,"\n")
cat("    Input phenotype file: ",pheno.file,"\n")
cat("    Regression model for interesting variables (e.g ~ AD + Braak): ",model_protect,"\n")
cat("    Regression model for batch variables (e.g. ~ Plate + RIN + PMI): ",model_remove,"\n")
cat("    Output prefix: ",out_prefix,"\n") 
cat("\n")

cat("Loading libraries...\n")
suppressMessages(library(stringr))
suppressMessages(library(limma))

############################################################
cat("Reading input data...\n")
if(str_ends(string = expr.file , pattern = ".rds")){
  expr_mat <- readRDS(expr.file)
}else{
  expr_mat <- read.table(file = expr.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

pheno <- read.csv(file = pheno.file , stringsAsFactors = F , row.names = 1)

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

if(model_protect==""){
  all_vars <- all.vars(as.formula(model_remove))
}else{
  all_vars <- c(all.vars(as.formula(model_protect)),all.vars(as.formula(model_remove)))
}

for (v in all_vars) {
  if(any(is.na(pheno[,v]) | trimws(pheno[,v])=="")){
    stop("Missing value found in ",v)
  }
}

##############################################################

cat("Removing batch effects using removeBatchEffect function in limma package...\n")

if(model_protect==""){
  protect_design <- matrix(1,ncol(expr_mat),1)
}else{
  protect_design <- model.matrix(as.formula(model_protect), data = pheno)
}

remove_covars <- model.matrix(as.formula(model_remove), data = pheno)[, -1]

clean_expr_mat <- limma::removeBatchEffect(
  x = expr_mat,
  covariates = remove_covars,
  design = protect_design
)

######################################################
cat("Saving regressed data...\n")
if(str_ends(string = expr.file , pattern = ".rds")){
  saveRDS(clean_expr_mat, file = paste0(out_prefix,".noBatch.rds"))
}else{
  write.table(clean_expr_mat , file = paste0(out_prefix , ".noBatch.tsv") , row.names = T , col.names = T , sep = "\t" , quote = F)
}
cat("All done!\n")
