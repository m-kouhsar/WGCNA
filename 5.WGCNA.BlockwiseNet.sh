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
#SBATCH --output=WGCNA.BlockwiseNet.%j.out

###################################################  Input Agruments  #########################################################

# Data_File:                     Methylation/Expression matrix is rds or tsv format. Genes/Probes must be in rows and samples must be in columns
# SoftPow:                       Soft threshold power calculated by the previous step
# Block_Size:                    maximum number of gene/probes in each block (maxBlockSize in blockwiseModules function)
# min_Module_Size:               minimum number of genes/probes in each module
# Save_TOM:                      Set it ti 'yes' if you want to save TOM matrices
# Plot_Dendro:                   Set it to 'yes' if you want to save denrograms
# Numeric_Labels:                set it to 'yes' if you want to have numeric module names (eg module1, module2, ...) instead of colors 
# OutPrefix:                     Output files prefix (can be included path. The directory will be created)
# ScriptDir:                     The WGCNA repository path in your computer (eg. ./WGCNA/)

###############################################################################################################################
Data_File=./Methylation.rds     
SoftPow=6
Block_Size=30000
min_Module_Size=100
Save_TOM=No
Plot_Dendro=No
Numeric_Labels=Yes
OutPrefix=./Results/Cohort1

ScriptDir=./WGCNA
###############################################################################################################################

Rscript "${ScriptDir}"/R/5.WGCNA.BlockwiseNet.R "$Data_File" "$SoftPow" "$Block_Size" "$min_Module_Size" "$Save_TOM" "$Plot_Dendro" "$Numeric_Labels" "$OutPrefix"

