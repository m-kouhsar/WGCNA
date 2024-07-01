#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=100:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

exp_file_1=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/mad.0.5.lm.v2.pow3/Pitts.All.Psycho.0.2.v2.lmCorrected.betas.rds
exp_file_2=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/BDR_data/BDR.All.Psycho.0.2.AD.braak.lmCorrected.betas.rds
net_file_1=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/mad.0.5.lm.v2.pow3/Pitts.All.Psycho.0.2.v2.lmCorrected.wgcna.network.rds
modules=seashell4,mediumpurple3,thistle2,green4,firebrick3,sienna2,lightcoral,antiquewhite,magenta1,lightcyan,burlywood,linen,thistle3,lavenderblush,darkgoldenrod1,tan2,greenyellow,whitesmoke,orange,midnightblue,mistyrose,darkred,lightblue1
max_size=1000
max_gold=1000
n_permut=1000
out_pref=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/mad.0.5.lm.v2.pow3/Pitts.BDR.All.Braak

script_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/Scripts

module load R

Rscript ${script_dir}/WGCNA.ModulePreservation.R  $exp_file_1  $exp_file_2  $net_file_1  $modules  $max_size $max_gold $n_permut $out_pref
