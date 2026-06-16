# WGCNA

This repository contains some scripts to run Weighted Gene Co-expression Network Analysis (WGCNA) on DMA Methylation or RNA Expression data. 
The expression (or methylation) matrix, in `tsv` or `rds` format,  with matched metadata in `csv` format are  requiered to run the analysis. You can start with count matrix and normalized beta/M values in expression and methylation data, respectively. 
All the input arguments need to be set in `sh` scripts and after that you can run each script simply by `bash` command. All results and plots will be save in the output directory you set. 
The analysis pipeline is designed based on the following steps:

## Step1: Gene or CpG filtering

Some basic filtering, for example removing low count genes in expression data and removing SNP probes in methylation data, can be done using `1.ProbeFiltering_Expression.sh` and `1.ProbeFiltering_Methylation.sh`scripts. 

## Step2: Removing batch effects

To remove the effect of any cofounder or batch in your data, you can use `2.BatchRemoval.sh`. It will use `limma::removeBatchEffect()` function to regress out the effect of cofounders.

## Step3: Filtering low variance genes or CpGs using Median Absolute Deviation (MAD)

The third optional but highly recommended step is removing low variance features from the data by filtering out those features with low MAD value in the data. 
You can run this step using `3.MADFiltering.sh`. 

## Step4: Finding approperiate soft-thresholding power for network construction

Before the network reconstruction and finding modules, you need to check the scale free topology index to find the best soft-thresholding power. 
You can use `4.WGCNA.PicSoftPower.sh`script to do this.

## Step5: Network reconstruction and module finding

In this step, you can find the co-expressed gene or co-methylated CpG modules in your data using `WGCNA::blockwiseModules()` function. 
The main output of this step would be a network object saved in `rds` format.

## Step6: Module-Trait relationship analysis

To run this step you need the network object in `rds` format and the metadata in `csv` format. 
The `6.ModuleTrait.sh` script can help you to check the relationship between the modules and your variables of interest using linear regression (`analysis_type = "lm"`), 
t-test or ANOVA and Tukey test (`analysis_type = "test"`) or even a simple Pearson or Spearman correlation test (`analysis_type = "cor"`). 

The recommended approach is using `test` method on all modules (set module argument to `all`) to find those modules significantly related to the trait and then useing `lm` 
method to confirm your finding by adjusting the Module-Trait test based on the cofounders (e.g. Sex and Age).

This script can also visualize Module-Trait relationships using a heatmap (`heatmap = yes`), 
box plot based on Module Eigengene (ME) and categorical variables and a scatter plot to show the correlation between ME and numeric variables (`ME_plot = yes`).

## Step7: Module Membership VS Gene Significance analysis

In this step, you can use `7.WGCAN.ModuleMembership.sh` script to calculate Module Membership (MM) and Gene Significance (GS) for all genes in each module based on the variable of interest (Trait) and visualized it in a scatter plot. 
The plot can show those genes/CpGs that are significantly related to the trait based on an input P-value threshold (`GS_legend_pvalue`). 
It can also tag the most significant genes/CpGs in the plot by adding their ID as label based on another p-value threshold (`GS_label_pvalue`).

For categorical trait, the script will use `limma::lmFit()` and a linear regression model to calculate GS. 
So, you can optionally add any cofounders to the model by setting `Cofounders_num` for numeric cofounders and `Cofounders_cat`for categorical cofounders. 

To have a general view of the data and find any possible outlier samples using principal componenet analysis and Mahalanobis distance, the `DataOutlierChecking.sh` is provided. 
Outlier samples can be removed using `OutlierRemoval.sh`. 

Some other `R` and `bash` scripts are also available. For example:

* `WGCNA.CytoscapeExport.sh` for exporting modules in [Cytoscape](https://cytoscape.org/) format.
* `WGCNA.ModulePreservation.sh` for module preservation analysis and visualization.
* `ME.box.plot.R` for creating a box plot based on Module Eigengene and trait variable in one or two datasets (see [this figure part C](https://alz-journals.onlinelibrary.wiley.com/cms/asset/62e080f4-55ea-4b2d-a90b-af5e16b4198f/alz14501-fig-0002-m.jpg)).
* `Merged.ModuleMembershipPlot.R` for visualizing the MM vs GS in a module in two different datasets in one scatter plot (see [this figure part D](https://alz-journals.onlinelibrary.wiley.com/cms/asset/62e080f4-55ea-4b2d-a90b-af5e16b4198f/alz14501-fig-0002-m.jpg))
