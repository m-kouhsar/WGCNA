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

beta_file=./betas.Regressed.rds
mad_thr=0.5
out_pref=./betas.Regressed

ScriptDir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Revision.June2024/Scripts


Rscript ${ScriptDir}/3.MAD.Filtering.R $beta_file $mad_thr $out_pref
