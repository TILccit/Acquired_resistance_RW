# Acquired_resistance_RW
Analysis scripts associated with the manuscript "Depth of Response Predicts Acquired Resistance after PD-1/PD-L1 blockade Across Melanoma, RCC, and NSCLC: A Nationwide Danish Cohort of 2 127 Responders"

## **Analysis Workflow** ##

The analysis workflow consists of four main scripts:
1) Survival_analysis_melanoma.Rmd
2) Survival_analysis_renal.Rmd
3) Survival_analysis_lung.Rmd
4) PropensityScore_joined_cohorts_mel_rcc_lung.Rmd

and two helper scripts defining functions:
* Survival_comparisons_function.R: used to facilitate the display of pairwise comparisons between the clinical endpoints of the cohorts (OS, PFS); 
* PSM_survival_function_updated.R: used to perform pairwise propensity score matching (based on the MatchIt package) and comparison of the matched cohorts;

**Scripts 1, 2, and 3** perform the following for the three cohorts:
* Data cleaning (i.e. removing duplicates and excluded cases - see methods at the manuscript) and naming harmonization;
* Multiple univariate analysis on all covariates;
* Multivariate analysis on all the covariates with p<0.2 in the previous analysis;
* Calculate 5-year survival estimates in complete and partial responders;
* Perform a 18-month landmark analyses to rule out immortal-time bias.
Those analyses are performed for both the OS and PFS clinical endpoints. 
Of note, in the script for the Melanoma cohort, the previous analyses have been performed also for the "Melanoma-specific survival" clinical endpoint

**Script 4** instead performs the following:
* Cohort joining;
* Create of descriptive statistics tables (Table 1, and Supplementary table 1 in the manuscript);
* Remove of cases with missing data (complete*case analysis);
* Compute pairwise comparisons in the clinical endpoints of the cohorts using helper script A;
* Perform propensity score matching and compare clinical endpoints in the matched cohorts using helper script B;
* Compute median/max follow up;
* Calculate the percentage of censoring;
* Study the differential distribution of complete vs partial responses among the cohorts;
* Calculate the median time to progression in each cohort.

## **Instructions** ##

To reproduce the analysis conducted in the paper, please follow these steps:
1.	Clone or download this repository to your local machine.
2.	Navigate to the "data" folder and ensure that all the required data files are present. These files serve as the input for the analysis scripts.
3.	Open the R Markdown scripts located in the main folder. These scripts outline the analysis workflow and provide instructions for reproducing the results.
4.	Ensure that the necessary R packages are installed and up to date. Refer to the SessionInfo.txt files within the respective "session_info" subfolders for the specific package versions used during each analysis.
5.	Set the working directory within each R Markdown script to the root folder of this project. This ensures that the relative paths to the data files and other resources are correctly resolved. Make sure that the default directory of all the chunks is set to the one also used in the console (change via global options). For knitting the .RMD file, change the root directory in the knit setup chunk at the top of the document.
6.	Execute the R Markdown scripts 1,2, and 3 in the specified order, following the instructions and guidelines provided within each script. Subsequently, run script 4.
7.	The outputs, including tables, figures and models, will be generated and stored as markdown reports in your work directory. The expected output has been copied to the "Reports" folder.
Refer to the generated outputs for the corresponding results mentioned in the paper.

Please note that reproducing the analysis may require sufficient computational resources and dependencies. In order to be able to run the scripts, several R packages need to be installed and loaded, details on the packages (incl. versions) are listed in the session_info files provided for each of the main scripts.

It is recommended to review the documentation and understand the workflow before running the scripts.

- - - -
For any further inquiries or clarifications, please contact **Marco Donia, M.D., Ph.D** at _marco.donia@regionh.dk_.
- - - -

## **Citation**  ##

Please acknowledge the original authors of the paper and cite the appropriate references when using or building upon this work.
If you use the code or findings from this project in your research, please consider citing:

XXXXXXXX