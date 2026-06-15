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


################################# Argument description ##################################################

# Net_file:                     Network file (result of blockwiseModules function) in rds format.  
# Expr_file:                    Expression/Methylation matrix in rds or tsv format. Genes/Probes must be in rows and samples must be in columns
# Pheno_file:                   Metadata in csv format
# Trait:                        Trait variable (a column name in metadata)
# Categorical_trait:            Set it to 'yes' if your trait is categorical (it will be used to calculate Gene Significance (GS) based on a linear regression analysis)
# Cofounders_num:               Numerical cofounders you need to add to the lm model in GS calculation
# Cofounders_cat:               Categorical cofounders you need to add to the lm model in GS calculation
# GS_legend_pvalue:             All genes/probes with GS P-value smaller than this threshold will be shown in deifferent color in module membership scatter plots.
# GS_label_pvalue:              The genes/probes with GS P-value smaller than this threshold will be labeled by gene/probe IDs in module membership scatter plots.
# Adjusted_pvalue_method:       The method to adjust Module Membership (MM) and GS P-values based on number of genes (all genes for GS and genes inside each module for MM).
# modules:                      The selected module (separated by camma) for this analysis. Set it to 'all' if you need to run the analysis on all modules in the Net_file.
# SoftPow:                      The soft power threshold (will be used for calculating Module Eigenegenes (MEs)
# out_prefix:                   Output files prefix (can be included path. The directory will be created)
# ScriptDir:                    The WGCNA repository path in your computer (eg. ./WGCNA/)

#########################################################################################################


Net_file="./WGCNA.Net.rds"
Expr_file="./betas.rds"
Pheno_file="./pheno.csv"
Trait="Group"
Categorical_trait="yes"
Cofounders_num="Age,CellProportion"
Cofounders_cat="Sex,Plate"
GS_legend_pvalue=0.05
GS_label_pvalue=1e-4
Adjusted_pvalue_method="BH"
modules="module28,module124"
SoftPow=10
out_prefix="./WGCNA/Results/ModuleMembership"

ScriptDir=./WGCNA
#############################################################################################################################

Rscript "${ScriptDir}"/R/6.WGCNA.ModuleTrait.R "$net_file" "$expr_file" "$pheno_file" "$covars_fact" "$covars_num" "$modules" \
									   "$analysis_type" "$calc_ME" "$SoftPow" "$corr_plot" "$box_plot" "$save_csv" "$out_prefix" 



