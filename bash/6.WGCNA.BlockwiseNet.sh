#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


Data_File=./betas.Regressed.MAD.0.5.rds
SoftPow=6
Block_Size=30000
min_Module_Size=100
Save_TOM=No
Plot_Dendro=No
OutPrefix=./betas.Regressed.MAD.0.5

ScriptDir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Revision.June2024/Scripts


Rscript ${ScriptDir}/6.WGCNA.BlockwiseNet.R $Data_File $SoftPow $Block_Size $min_Module_Size $Save_TOM $Plot_Dendro $OutPrefix

