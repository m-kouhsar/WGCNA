args <- commandArgs(TRUE)

input_expression <- args[1]
input_phenotype <- args[2]
trait <- args[3]
variables_fact <- args[4]
variables_num <- args[5]
model_lm <- args[6]
out_pref <- args[7]

sink(paste0(out_pref,".RegressOut.log.txt"))

cat("Input arguments:\n")
cat("    Input expression file: ",input_expression,"\n")
cat("    Input phenotype file: ",input_phenotype,"\n")
cat("    Trait variable: ",trait,"\n")
cat("    Factor variables: ",variables_fact,"\n")
cat("    Numeric variables: ",variables_num,"\n")
cat("    Regression model: ",model_lm,"\n")
cat("    Output prefix: ",out_pref,"\n") 
cat("\n")
cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dendextend))
suppressMessages(library(parallel))


expr <- readRDS(input_expression)
pheno <- read.csv(file = input_phenotype , stringsAsFactors = F , row.names = 1)

cat("Are phenotype file and expr file matched? ",ifelse(identical(rownames(pheno),colnames(expr)),"Yes","No"))
cat("\n")

variables_num <- str_split(variables_num,pattern = ',',simplify = T)[1,]
variables_fact <- str_split(variables_fact,pattern = ',',simplify = T)[1,]

pheno2 <- as.numeric(as.factor(pheno[,trait]))
for (i in 1:length(variables_fact)) {
  pheno2 <- cbind.data.frame(pheno2,as.numeric(as.factor(pheno[,variables_fact[i]])))
}
for (i in 1:length(variables_num)) {
  pheno2 <- cbind.data.frame(pheno2,as.numeric(pheno[,variables_num[i]]))
}
names(pheno2) <- c(trait,variables_fact,variables_num)

##############################################################

cat("Regressing out covariates using linear regression...\n")
fit.lm <- function(x,model=model_lm,data=pheno2){
  x <- as.numeric(x)
  lm_vars <- all.vars(as.formula(model))
  for(j in 1:length(lm_vars)){
      
      assign(lm_vars[j], data[,lm_vars[j]])
  }
  model <- paste0("x",model)
  fit<-try(lm(formula = as.formula(model),na.action =na.omit))
  if(inherits(fit,'try-error')){
    return(rep(NA,nrow(data)))
  }else{
    return(fit$residuals) 
  }
}

nThreads <- max(1, parallel::detectCores(), na.rm = TRUE)-1

cl<- makeCluster(nThreads)
expr.norm<-t(parApply(cl, expr , 1, fit.lm, model_lm,pheno2))

stopCluster(cl) 

expr.norm <- as.data.frame(expr.norm)
colnames(expr.norm) = colnames(expr)
expr = expr.norm

######################################################
cat("Saving regressed data...\n")
saveRDS(expr, file = paste0(out_pref,".Regressed.rds"))
cat("\n")
cat("All done!\n")
sink()
