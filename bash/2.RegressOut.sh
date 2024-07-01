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

input_expression=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Raw/Pitts.NoControl.Quantile.0.2.prepared.rds
input_phenotype=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Raw/Pitts.Pheno.NoControl.0.2.csv
trait=Psychosis
variables_fact=Sex,BraakStage,Plate,TissueType
variables_num=Age,CellProportion
model_lm=~Age+CellProportion+Sex+Plate+BraakStage+TissueType
out_pref=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Revision.June2024/PITT


ScriptDir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/wgcna/Revision.June2024/Scripts


Rscript ${ScriptDir}/2.RegressOut.R $input_expression $input_phenotype $trait $variables_fact $variables_num $model_lm $out_pref
