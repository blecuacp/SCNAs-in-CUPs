# Project Title

Here one can find data and code to reproduce the K_c constants of the paper: "A novel workflow for computational identification of somatic copy number alterations using DNA methylation microarrays in cancers of unknown primary", as well as key elements of Fig.1 (i.e., Figures 1c and 1d).

The necessary input files (or extra files needed to run the code) can be found under Data/. All the necessary (comented) code is under Code/.

## Description

We suggest the interested party to run the code as follows:

* 1- Code to download the TCGA necessary data can be found in Code/Biolinks.R:
  * 1.a Run this code once for the Training Set (Data/Additional.File.2_TableS1.csv). 
  * 1.b Run this code once for the Validation Set (Data/Additional.File.4_TableS3.csv).

2- Code to process and QC the methyaltion arrays (.idat files for both cohorts) can be found under Code/Run_Minfi.R:
	* 2.a Run this code once for the Training Set.
	* 2.b Run this code once for the Validation Set.

3- Code to download the whole blood controls from GEO is provided in Code/GEO.R.

4- Run Code/Run_Minfi.R to process the raw .idat files for the whole blood control group above in 3).

5- To determine the K_c constants mentioned in the text (Fig. 1a), we provide code (for c="Amp10" as an illustrative example) under Code/Kc-fit.R:
	* 5.a Run it once for the Training Cohort.
	* 5.b To reproduce the other K_c constants is straight forward.

6- Run Code/Run_CONUMEE.R for the Validation Cohort only as Query and whole blood cohort as Control (thhis step may take several hours):
	* 6.a Use as input data the output from Run_Minfi.R for the Validation cohort (Query) and 
	* 6.b Use as input data the output from Run_Minfi.R for the GEO data set (Control).
	* 6.c Segemented files and full genome (excluding sex chromosomes) plots will be generated.

7- To Calibrate and determine the diferent Copy Number Calls from the previous step, run Code/Calibration_CONUMEE_ValidationSet.R:
	* 7a Produce Segmented files "withCalls".

8- Reproduce Figures 1d) and 1c):
	* 8.a True Positive Rate: Run Code/Comparison_CONUMEE_ASCAT_TruePositiveRate.R for the Validation Cohort.
	* 8.b Fase Positive Rate: Run Code/Comparison_CONUMEE_ASCAT_FalsePositiveRate.R for the Validation Cohort.
	* 8.c Run Code/Figure1c.R
	* 8.d Run Code/Figure1d.R 

## Authors

Contributors names and contact info

Pedro Blecua  
pblecua@carrerasresearch.org

