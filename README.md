# WGCNA

This repository contains scripts for running Weighted Gene Co-expression Network Analysis (WGCNA) on DNA methylation or RNA expression data.

An expression (or methylation) matrix in `tsv` or `rds` format, together with matched metadata in `csv` format, is required to run the analysis. You can start with a count matrix for expression data or normalized beta/M values for methylation data.

All input arguments should be set in the corresponding `sh` scripts. After that, each script can be executed using the `bash` command. All results and plots will be saved to the output directory you specify.

## Analysis Steps

The analysis is based on the following steps:

### Step 1: Gene or CpG Filtering

Basic filtering, such as removing low-count genes in expression data or removing SNP probes in methylation data, can be performed using `1.ProbeFiltering_Expression.sh` and `1.ProbeFiltering_Methylation.sh`.

### Step 2: Removing Batch Effects

To remove the effects of confounders or batch variables in your data, you can use `2.BatchRemoval.sh`. This script uses the `limma::removeBatchEffect()` function to regress out the effects of confounding variables.

### Step 3: Filtering Low-Variance Genes or CpGs Using Median Absolute Deviation (MAD)

The third step, which is optional but highly recommended, removes low-variance features from the data by filtering features with low MAD values.

You can run this step using `3.MADFiltering.sh`.

### Step 4: Finding an Appropriate Soft-Thresholding Power for Network Construction

Before network construction and module detection, you should evaluate the scale-free topology fit index to determine the optimal soft-thresholding power.

You can use `4.WGCNA.PicSoftPower.sh` to perform this step.

### Step 5: Network Construction and Module Detection

In this step, co-expressed gene modules or co-methylated CpG modules are identified using the `WGCNA::blockwiseModules()` function.

The main output of this step is a network object saved in `rds` format.

### Step 6: Module–Trait Relationship Analysis

To run this step, you need the network object in `rds` format and the metadata in `csv` format.

The `6.ModuleTrait.sh` script can be used to assess relationships between modules and variables of interest using linear regression (`analysis_type = "lm"`), t-tests or ANOVA with Tukey post-hoc tests (`analysis_type = "test"`), or Pearson/Spearman correlation analysis (`analysis_type = "cor"`).

The recommended approach is to use the `test` method on all modules (`module = "all"`) to identify modules significantly associated with the trait and then use the `lm` method to confirm these associations while adjusting for confounders (e.g., sex and age).

This script can also visualize module–trait relationships using a heatmap (`heatmap = yes`), box plots of Module Eigengenes (ME) against categorical variables, and scatter plots showing correlations between ME values and numeric variables (`ME_plot = yes`).

### Step 7: Module Membership vs. Gene Significance Analysis

In this step, `7.WGCNA.ModuleMembership.sh` can be used to calculate Module Membership (MM) and Gene Significance (GS) for all genes within each module with respect to a trait of interest and to visualize the results in a scatter plot.

The plot can highlight genes/CpGs significantly associated with the trait based on a user-defined p-value threshold (`GS_legend_pvalue`). It can also label the most significant genes/CpGs using a separate threshold (`GS_label_pvalue`).

For categorical traits, the script uses `limma::lmFit()` and a linear model to calculate GS values.

You can optionally adjust for confounders by specifying `Cofounders_num` for numeric confounders and `Cofounders_cat` for categorical confounders.

A `csv` file will be saved for each module, containing MM and GS values along with their corresponding p-values. These files can be used to identify hub genes/CpGs within the modules for downstream analyses.

## Other Available Scripts

To obtain an overview of the data and identify potential outlier samples using principal component analysis (PCA) and Mahalanobis distance, `DataOutlierChecking.sh` is provided. Outlier samples can then be removed using `OutlierRemoval.sh`.

Some additional `R` and `bash` scripts are also available. For example:

* `WGCNA.CytoscapeExport.sh` for exporting modules in [Cytoscape](https://cytoscape.org/) format.
* `WGCNA.ModulePreservation.sh` for module preservation analysis and visualization.
* `ME.box.plot.R` for creating box plots based on Module Eigengenes and trait variables in one or two datasets (see [Figure 2C](https://alz-journals.onlinelibrary.wiley.com/cms/asset/62e080f4-55ea-4b2d-a90b-af5e16b4198f/alz14501-fig-0002-m.jpg)).
* `Merged.ModuleMembershipPlot.R` for visualizing MM vs. GS relationships within a module across two different datasets in a single scatter plot (see [Figure 2D](https://alz-journals.onlinelibrary.wiley.com/cms/asset/62e080f4-55ea-4b2d-a90b-af5e16b4198f/alz14501-fig-0002-m.jpg)).
