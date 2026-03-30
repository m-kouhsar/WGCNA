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

beta_file=./Methylation.rds
meth_filter=yes     # set it to y/yes/T/true/1 if you need to filter the SNP,Sex and Cross Hybridising probes
convert2M=yes       # set it to y/yes/T/true/1 if you need to convert Beta values to M values
mad_thr=0.25
out_prefix=./Results/Methylation.Cohort1

ScriptDir=./WGCNS
##########################################################################

Rscript "${ScriptDir}"/R/1.ProbeFiltering_Methylation.R "$beta_file" "$meth_filter" "$convert2M" "$mad_thr" "$out_prefix" "$ScriptDir"
