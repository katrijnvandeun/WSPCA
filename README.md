# WSPCA
Weighted sparse PCA

This github repository contains all matlab and R code and functions that support the weighted sparse PCA manuscript. Also included are (derived) data and results. The repository contains a MATLAB branch for all matlab code, a R branch for all R code, a DATA branch with data that were derived from public data, and a RESULTS branch with the results obtained on simulated and empirical data.

## 1. CREATE DATA
  
-Gene expression data: Run the R scripts *R/SCRIPT_ProcessingGSE29614_2007.R* and *R/SCRIPT_ProcessingGSE29617_2008.R*. This retrieves the expression data from the online database and creates RMA pre-processed expression matrices for both seasons. The resulting pre-processed gene expression data are stored in *DATA/GSE29614_rma.txt* (2008 season) and *DATA/GSE29617_rma.txt* (2007 season). Besides, phenotype data with subject ids and the number of days after vaccination are created (*DATA/GSE29614_pdata.txt* and *DATA/GSE29617_pdata.txt*) along with a table containing the variance parameters in the two-component model resulting from the GESTr package (*DATA/GSE29614__RLparsD0.txt* and *DATA/GSE29617__RLparsD0.txt*). 
-Antibody titers:  These can be created also from the R scripts *R/SCRIPT_ProcessingGSE29614_2007.R* and *R/SCRIPT_ProcessingGSE29617_2008.R*. In these scripts the antibody titers are extracted from the CEL files and written to *DATA/GSE29614_titers.txt* (2008 season) and *DATA/GSE29617_titers.txt* (2007 season). 
  
These files contain the antibody titers for three influenza strains and measured just before and 28 days after vaccination

Because we use the baseline-corrected gene expression data obtained three days after vaccination, the day 0 and day 3 data have to be matched on their subject identifiers and the difference scores have to be calculated. Note that the analyses are on the raw data (no centering nor scaling to unit variance). To reproduce these steps, run the following MATLAB scripts.  

 a. 2007 sample  
  * *MATLAB/Script_Preprocessing_2007.m* creates *DATA/DATA2007.mat* and *TIVD3_2007_rev.txt*
  * *MATLAB/ScriptHAI_TIVD28vsD0_2007.m* creates *TIVtiter2007.m* and *TIVtiter2007.txt*  
  
 b. 2008 sample  
  * *MATLAB/Script_Preprocessing_2008.m* creates *DATA/DATA2008.mat* and *DATA/TIVD3_rev.txt*
  * *MATLAB/ScriptHAI_TIVD28vsD0_2008.m* creates *DATA/TIVtiter.m* and *DATA/TIVtiter.txt*

## 2. ANALYZE DATA + POST-PROCESS RESULTS

1. PMA (for comparison) using R and the package PMA  
  *R/Script_sgcca_spls.R*
2. Ordinary PCovR analysis using MATLAB and creation of the Figures 1 & 2 in the paper  
  *MATLAB/plot_PCovR.m*: this script requires two external matlab functions: fig.m and exportfig.m
	Available from: http://www.mathworks.com/matlabcentral/fileexchange/30736 and
	https://nl.mathworks.com/matlabcentral/fileexchange/727-exportfig
3. Sparse PCovR using MATLAB  
  *MATLAB/Script_SPCovRanalysis.m*: Calls different function that implement Algorithm 1 and Algorithm 2.


## 3. ANNOTATION OF SELECTED PROBE SETS

Annotate the probe-sets with non-zero weights resulting from SPCovR and SGCCA, SPLS
  * For SGCCA and SPLS: use *DATA/GENEIDS_spls.txt* and *DATA/GENEIDS_sgcca.txt* as input to AmiGO: http://amigo.geneontology.org/amigo
  * For SPCovR, retrieve the AFFYIDS within *MATLAB/Script_SPCovRanalysis.m* (see the bottom of this file) and convert these to the official gene symbols using DAVID (as the annotation file included in the data folder may be outdated): https://david.ncifcrf.gov/. 
Next submit the official gene symbols to AmiGO.

## 4. SIMULATION STUDY
