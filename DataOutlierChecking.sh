#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=01:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#########################################################################

expr_file=./methylation.rds   #Methylation matrix in rds format or expression matrix in tsv format
pheno_file=./phenotype.csv
var_cat=Sex,Plate
var_num=Age,PMI
OutPrefix=./Cohort1

ScriptDir=./WGCNA
##########################################################################

Rscript "$ScriptDir"/R/DataOutlierChecking.R "$ScriptDir" "$expr_file" "$pheno_file" "$var_cat" "$var_num" "$OutPrefix"
