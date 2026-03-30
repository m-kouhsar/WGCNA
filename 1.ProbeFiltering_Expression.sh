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

counts_file=./Expression.counts.tsv
pheno_file=./Phenotype.csv
var_trait=Braak
normalize_method=vst
OutPrefix=./Results/Cohort1

ScriptDir=./WGCNA
##########################################################################

Rscript "${ScriptDir}"/R/1.ProbeFiltering_Expression.R "$counts_file" "$pheno_file" "$var_trait" "$normalize_method" "$OutPrefix"

