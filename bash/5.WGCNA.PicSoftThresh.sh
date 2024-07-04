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


data_file=./betas.Regressed.MAD.0.5.rds
block_size=30000
out_pref=./betas.Regressed.MAD.0.5


script_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Revision.June2024/Scripts


Rscript $script_dir/5.WGCNA.PicSoftThresh.R $data_file $block_size $out_pref


