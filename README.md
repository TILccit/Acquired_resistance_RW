# Acquired_resistance_RW
Analysis scripts associated with the manuscript "Depth of Response Predicts Acquired Resistance after PD-1/PD-L1 blockade Across Melanoma, RCC, and NSCLC: A Nationwide Danish Cohort of 2 127 Responders"

The analysis workflow consists of four main scripts:
1) Survival_analysis_melanoma.Rmd
2) Survival_analysis_renal.Rmd
3) Survival_analysis_lung.Rmd
4) PropensityScore_joined_cohorts_mel_rcc_lung.Rmd

and two helper scripts defining functions:
A) Survival_comparisons_function.R: used to facilitate the display of pairwise comparisons between the clinical endpoints of the cohorts (OS, PFS); 
B) PSM_survival_function_updated.R: used to perform pairwise propensity score matching (based on the MatchIt package) and comparison of the matched cohorts;

Scripts 1, 2, and 3 perform the following for the three cohorts:
- Data cleaning (i.e. removing duplicates and excluded cases - see methods at the manuscript) and naming harmonization;
- Multiple univariate analysis on all covariates;
- Multivariate analysis on all the covariates with p<0.2 in the previous analysis;
- Calculate 5-year survival estimates in complete and partial responders;
- Perform a 18-month landmark analyses to rule out immortal-time bias.
Those analyses are performed for both the OS and PFS clinical endpoints. 
Of note, in the script for the Melanoma cohort, the previous analyses have been performed also for the "Melanoma-specific survival" clinical endpoint

Script 4 instead performs the following:
- Cohort joining;
- Create of descriptive statistics tables (Table 1, and Supplementary table 1 in the manuscript);
- Remove of cases with missing data (complete-case analysis);
- Compute pairwise comparisons in the clinical endpoints of the cohorts using helper script A;
- Perform propensity score matching and compare clinical endpoints in the matched cohorts;
- Compute median/max follow up;
- Calculate the percentage of censoring;
- Study the differential distribution of complete vs partial responses among the cohorts;
- Calculate the median time to progression in each cohort.