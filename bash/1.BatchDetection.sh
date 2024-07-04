#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


Data_File=./betas.rds
Pheno_File=./pheno.csv
Covars_Factor=Sex,Plate
Covars_Num=Age,CellProportion
Out_Prefix=./betas

ScriptDir=./R

Rscript ${ScriptDir}/1.BatchDetection.R $Data_File $Pheno_File $Covars_Factor $Covars_Num $Out_Prefix 

