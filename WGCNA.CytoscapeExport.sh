#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=4:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err


net_file=./WGCNA.Net.rds
module_name=darkgreen
weighted=T
threshold=0.03
out_pref=./betas
TOM_Dir=

script_dir=./R


Rscript $script_dir/WGCNA.CytoscapeExport.R $net_file $module_name $weighted $threshold $out_pref $TOM_dir


