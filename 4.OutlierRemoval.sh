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

expr_file=./Methylation.rds  #Methylation matrix in rds format or expression matrix in tsv format
pheno_file=./phenotype.csv
outliers=sample1,sample2,sample3
out_prefix=./Results/Methylation.Cohort1

ScriptDir=./WGCNS
##########################################################################

Rscript "${ScriptDir}"/R/4.OutlierRemoval.R "$expr_file" "$pheno_file" "$outliers" "$out_prefix"
