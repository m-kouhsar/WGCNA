#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


################################# Argument description ##################################################

# modules = Names of modules separated with comma, "all" for running analysis on all modules
# analysis.type= one of the following options (or more than one separated by comma)
#        "p" -> Pearson correlation
#        "s" -> Spearman correlation
#        "l" -> Linear regression
#        "ps" -> Pearson correlation for numeric variables and Spearman correlation for categorical variables
#        "tps" -> Pearson correlation for two numeric variables,
#                   Spearman correlation for two categorical variables or one categorical (more than one groups) 
#                   and one numeric variables, t-test for a binary and a numeric variable

# calc_ME= "yes" or "no", Calculate Module Eigengene from expression matrix
# corr.plot= "yes" or "no" , correlation plot between MEs and phenotype of interest
# scatter.plot= "yes", "no" , Scatter plot for MEs and Trait variable
# save_csv = "yes" , or "no", Save the results in csv format

#########################################################################################################


net_file="./WGCNA.Net.rds"
expr_file="./betas.rds"
pheno_file="./pheno.csv"
trait=Psychosis
covars_fact=Sex,BraakStage,Plate
covars_num=Age,CellProportion
modules=all
analysis_type=tsp,s
calc_ME=no
SoftPow=3
corr_plot=no
scatter_plot=no
save_csv=yes
out_pref="./results_moduleTrait"

ScriptDir=./R

Rscript ${ScriptDir}/WGCNA.ModuleTrait.R $net_file $expr_file $pheno_file $trait $covars_fact $covars_num $modules \
									   $analysis_type $calc_ME $SoftPow $corr_plot $scatter_plot $save_csv $out_pref 



