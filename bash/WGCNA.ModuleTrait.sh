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

net_file="./betas.Regressed.MAD.0.5.WGCNA.Net.rds"
expr_file="./betas.Regressed.MAD.0.5.rds"
pheno_file="./betas.csv"
trait=Psychosis
covars_fact=Sex,BraakStage,Plate
covars_num=Age,CellProportion
modules=all
analysis_type=tsp,s
calc_ME=no
SoftPow=3
corr_plot=no
scatter_plot=no
save_csv=yes
out_pref="./betas.Regressed.MAD.0.5.WGCNA"

ScriptDir=./R

Rscript ${ScriptDir}/WGCNA.ModuleTrait.R $net_file $expr_file $pheno_file $trait $covars_fact $covars_num $modules \
									   $analysis_type $calc_ME $SoftPow $corr_plot $scatter_plot $save_csv $out_pref 



