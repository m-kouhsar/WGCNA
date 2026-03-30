#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

expr_file="./Results/WGCNA/Immunisation.Mval.Filtered.MAD.25.rds"  #Expression data in tsv format or methylation data in rds format
pheno_file="./Raw/Immunisation.pheno.csv"
model_protect="~AD+Braak+Immunised"                                #Regression model contains all Variable of interest (optional)
model_remove="~Sex+Age+Plate"                                      #Regression model contains all batch variables (DO NOT add your variable of interest)
out_prefix=./Results/WGCNA/Immunisation.Mval.Filtered.MAD.25


ScriptDir=./WGCNA
#################################################################################################################

Rscript "${ScriptDir}"/R/4.BatchRemoval.R "$expr_file" "$pheno_file" "$model_protect" "$model_remove" "$out_prefix"
