#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=24:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


data_file=./Methylation.rds #Methylation data in rds format or Expression data in tsv format
block_size=30000
out_prefix=./Results/Cohort1


script_dir=./WGCNA
#################################################################################

Rscript "$script_dir"/R/4.WGCNA.PicSoftPower.R "$data_file" "$block_size" "$out_prefix"


