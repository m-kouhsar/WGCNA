#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out


################################# Argument description ##################################################

# MM_GS_file: Module Membership and Gene Significance csv file (for an specific module) obtained from step 7. This is a csv file contains the following columns:
#				ID:         Gene or probe ID
#           	MM:         Module Membership value
#               MM.Pval:    Module Membership P-value
#           	GS:         Gene Significance
#           	GS.Pval:    Gene Significance P-value
#
# Net_file:   WGCNA network opbject file in rds format. This is the results of step 5
# id_type:    gene/probe id type. can be 'cpg' for methylation data and CpG IDs and 'entrez','symbol' or 'ensembl' for expression data
# MM:         Module Membership threshold for selecting hub genes
# MM_pval:    Module Membership P-value threshold for selecting hub genes
# GS: Gene    Significance threshold for selecting hub genes
# GS_pval:    Gene Significance P-value for selecting hub genes. 
# out_prefix: Output files prefix.
# ScriptDir:  repository directory (./WGCNA)

#########################################################################################################

MM_GS_file="./module28.csv"
Net_file="./net.rds"
id_type="cpg"                 
MM=0.8
MM_pval=1
GS=0
GS_pval=0.05
out_prefix="./Results/Enrichment/module28"

ScriptDir="./WGCNA"
#########################################################################

Rscript "$ScriptDir"/R/Enrichment.R "$MM_GS_file" "$Net_file" "$methylation" "$MM" "$MM_pval" "$GS" "$GS_pval" "$out_prefix"

