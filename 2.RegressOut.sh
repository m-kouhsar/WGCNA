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

input_expression=./betas.rds
input_phenotype=./pheno.csv
variables_fact=Sex
variables_num=Age
model_lm=~Age+Sex
out_pref=./data


ScriptDir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Revision.June2024/Scripts


Rscript ${ScriptDir}/2.RegressOut.R $input_expression $input_phenotype $variables_fact $variables_num $model_lm $out_pref
