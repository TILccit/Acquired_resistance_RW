Survival analysis renal
================
Mario Presti
First created on Feb 2025. Updated on 06 August 2025

- [Introduction](#introduction)
- [Loading Libraries](#loading-libraries)
  - [Pre-processing of data](#pre-processing-of-data)
- [OS analysis](#os-analysis)
  - [Multivariate Survival Analysis - only complete
    cases](#multivariate-survival-analysis---only-complete-cases)
- [Combined Univariable vs Multivariable Forest Plots -
  OS](#combined-univariable-vs-multivariable-forest-plots---os)
- [PFS analysis](#pfs-analysis)
  - [Multiple Univariate analysis - Incomplete complete
    cases](#multiple-univariate-analysis---incomplete-complete-cases)
- [Multivariate Survival Analysis](#multivariate-survival-analysis)
  - [All responders - PFS](#all-responders---pfs)
- [Combined Univariable vs Multivariable Forest Plots -
  PFS](#combined-univariable-vs-multivariable-forest-plots---pfs)
- [Five-Year Survival Estimates by
  Response](#five-year-survival-estimates-by-response)

# Introduction

Multivariate Survival Analysis on the Renal database.

# Loading Libraries

``` r
# Load required libraries
library(pacman)
pacman::p_load(data.table,stringr,tidyverse,ggplot2,MatchIt,survival,survminer,survey,compareGroups, forestmodel, forestplot, openxlsx,kableExtra,ggsurvfit,extrafont,DT,tibble,strex,hrbrthemes,ggstatsplot,reshape,pander,ggrepel,scales,dataMaid,gridExtra,tidytidbits,survivalAnalysis,gtsummary, Cairo,Amelia,officer,mice,naniar)
```

## Pre-processing of data

- Check that all columns have no spaces etc, especially character
  columns,
- transform character columns to factors.
- Check if all factors are ok.
- Check numerical columns.

We first need to define which columns we want as character
(categorical), and which need to be numerical. Sometimes things get
messed up with excel files.

``` r
setwd("E:/PhD_projects/Realworld/Data/")
dataDF <- read.xlsx("renal.xlsx")

# Check for duplicates
any(duplicated(dataDF$patient_id))
```

    ## [1] FALSE

``` r
# Remove duplicates
dataDF <- dataDF[!duplicated(dataDF), ]
dim(dataDF)
```

    ## [1] 261  49

``` r
# Rename some columns
names(dataDF)[names(dataDF) == 'progression_Marco'] <- 'Progressed'
names(dataDF)[names(dataDF) == 'age'] <- 'age_1st_treat'
names(dataDF)[names(dataDF) == 'Treatment_line'] <- 'line_correct'
names(dataDF)[names(dataDF) == 'Treatment'] <- 'regime_correct'
names(dataDF)[names(dataDF) == 'CNS_mets'] <- 'Brain_metastases'
names(dataDF)[names(dataDF) == 'BOR'] <- 'bor'
names(dataDF)[names(dataDF) == 'ecog_ps.1'] <- 'PS'

# Define the features of interest
features <- c("sex", "age_1st_treat","line_correct", "regime_correct", 
              "Brain_metastases","PS", "bor")

features_cancer_spec <- c("clear_cell","sarcomatoid", "smoking", "IMDC")

## Define which columns include survival data
OS_time <- c("OS_days")
OS_status <- c("Dead")
PFS_time <- c("PFS_days")
PFS_status <- c("Progressed")

## DEFINE WHICH COLUMNS ARE CATEGORICAL
categ_feats <-c("sex","regime_correct","PS", "bor", "line_correct",features_cancer_spec)

## DEFINE WHICH COLUMNS ARE NUMERICAL
numeric_feats <- c("OS_days", "Dead","PFS_days","Progressed","age_1st_treat")

# Remove spaces from column names and features.
colnames(dataDF) <- gsub(" ","_",colnames(dataDF))
features <- gsub(" ","_",features)
categ_feats <- gsub(" ","_",categ_feats)
numeric_feats <- gsub(" ","_",numeric_feats)
OS_time <- gsub(" ","_",OS_time)
OS_status <- gsub(" ","_",OS_status)
PFS_time <- gsub(" ","_",PFS_time)
PFS_status <- gsub(" ","_",PFS_status)

## MAKE SURE THERE ARE NO BLANKS ANYWHERE IN THE DATA
## Find non-Numeric columns
num_cols <- sapply(dataDF, is.numeric)

dataDF <- cbind(apply(dataDF[,!num_cols], 2, str_remove_all, "NA"),dataDF[,num_cols])   # Remove blanks

### MAKE SURE THAT THE COLUMNS ARE IN THE CORRECT FORMAT NEEDED FOR ANALYSIS.
## Numeric columns
dataDF <- dataDF %>% mutate(across(all_of(numeric_feats), as.numeric))

## Categorical
dataDF <- dataDF %>% mutate(across(all_of(categ_feats), as.character))

# Transform character columns to categorical, to have the different levels
dataDF <- dataDF %>% dplyr::mutate(across(all_of(categ_feats), as.factor))
dataDF <- dataDF %>% dplyr::select(all_of(c(colnames(dataDF)[1], OS_time,OS_status, PFS_time,PFS_status, features, features_cancer_spec))) %>% distinct()


print("Numbers per treatment")
```

    ## [1] "Numbers per treatment"

``` r
table(dataDF$regime_correct)
```

    ## 
    ##  ipi/nivo nivolumab 
    ##       200        61

``` r
print("Numbers per response")
```

    ## [1] "Numbers per response"

``` r
table(dataDF$bor)
```

    ## 
    ##  CR  PR 
    ##  79 182

``` r
## SUMMARY OF THE DATA
summary(dataDF)
```

    ##   record_id            OS_days          Dead           PFS_days     
    ##  Length:261         Min.   : 164   Min.   :0.0000   Min.   :  60.0  
    ##  Class :character   1st Qu.:1021   1st Qu.:0.0000   1st Qu.: 384.0  
    ##  Mode  :character   Median :1370   Median :0.0000   Median : 918.0  
    ##                     Mean   :1372   Mean   :0.3333   Mean   : 953.8  
    ##                     3rd Qu.:1770   3rd Qu.:1.0000   3rd Qu.:1412.0  
    ##                     Max.   :3313   Max.   :1.0000   Max.   :2756.0  
    ##                                                                     
    ##    Progressed       sex      age_1st_treat   line_correct   regime_correct
    ##  Min.   :0.0000   1   :200   Min.   :27.00   First:189    ipi/nivo :200   
    ##  1st Qu.:0.0000   2   : 59   1st Qu.:55.00   Other: 72    nivolumab: 61   
    ##  Median :1.0000   NA's:  2   Median :63.00                                
    ##  Mean   :0.6054              Mean   :62.19                                
    ##  3rd Qu.:1.0000              3rd Qu.:69.50                                
    ##  Max.   :1.0000              Max.   :85.00                                
    ##                              NA's   :2                                    
    ##  Brain_metastases      PS      bor      clear_cell sarcomatoid smoking   
    ##  Length:261         0   :128   CR: 79   ja :225    ja : 79     0   : 74  
    ##  Class :character   1   :111   PR:182   nej: 36    nej:182     1   : 39  
    ##  Mode  :character   2   : 18                                   2   :125  
    ##                     3   :  1                                   99  :  5  
    ##                     NA's:  3                                   NA's: 18  
    ##                                                                          
    ##                                                                          
    ##                 IMDC    
    ##  Good/Intermediate:187  
    ##  Poor             : 74  
    ##                         
    ##                         
    ##                         
    ##                         
    ## 

``` r
# remove some white spaces in the PS column
dataDF$PS <- trimws(dataDF$PS)
dataDF$PS <- as.factor(dataDF$PS)

# Rename some factors for plots
dataDF <- dataDF %>%
  dplyr::mutate(
    Sex = case_when(
      sex == "2"   ~ "Female",
      sex == "1"   ~ "Male",
      sex == "99"  ~ NA_character_,
      TRUE         ~ NA_character_
    ),
    Sex = factor(Sex, levels = c("Male", "Female")),

    regime_correct = case_when(
      regime_correct == "nivolumab" ~ "Anti-PD1",
      regime_correct == "ipi/nivo"  ~ "Ipi+Nivo",
      regime_correct == "99"        ~ NA_character_,
      TRUE                           ~ NA_character_
    ),
    regime_correct = factor(
      regime_correct,
      levels = c("Anti-PD1", "Ipi+Nivo")
    ),

    PS = case_when(
      PS %in% c("1", "2", "3") ~ "PS≥1",
      PS == "0"               ~ "PS=0",
      PS == "99"              ~ NA_character_,
      TRUE                    ~ NA_character_
    ),
    PS = factor(PS, levels = c("PS=0", "PS≥1")),
    
    regime_correct = case_when(
      regime_correct == "Ipi+Nivo" ~ "Anti-PD1+Anti-CTLA4",
      regime_correct == "Anti-PD1" ~ "Anti-PD1",
      regime_correct == "99"       ~ NA_character_,
      TRUE                    ~ NA_character_
    ),
    regime_correct = factor(regime_correct, levels = c("Anti-PD1", "Anti-PD1+Anti-CTLA4")),

    smoking = case_when(
      smoking == "0"           ~ "Non smoker",
      smoking == "1" ~ "Current smoker",
      smoking == "2" ~ "Former smoker",
      smoking == "99"          ~ NA_character_,
      TRUE                     ~ NA_character_
    ),
    smoking = factor(
      smoking,
      levels = c("Non smoker", "Current smoker","Former smoker")
    ),
    
    Brain_metastases = case_when(
     Brain_metastases == "No"           ~ "No",
     Brain_metastases == "Yes" ~ "Yes",
      TRUE                     ~ NA_character_
    ),
    Brain_metastases = factor(
      Brain_metastases,
      levels = c("No", "Yes")
    ),

    clear_cell = case_when(
      clear_cell == "ja"   ~ "Clear cell",
      clear_cell == "nej"  ~ "Non clear cell",
      clear_cell == "99"   ~ NA_character_,
      TRUE                  ~ NA_character_
    ),
    clear_cell = factor(
      clear_cell,
      levels = c("Clear cell", "Non clear cell")
    ),
    sarcomatoid = case_when(
      sarcomatoid == "ja"   ~ "Sarcomatoid",
      sarcomatoid == "nej"  ~ "Non sarcomatoid",
      sarcomatoid == "99"   ~ NA_character_,
      TRUE                  ~ NA_character_
    ),
    sarcomatoid = factor(
      sarcomatoid,
      levels = c("Sarcomatoid", "Non sarcomatoid")
    ),
    IMDC = case_when(
      IMDC == "Poor"   ~ "Poor",
      IMDC == "Good/Intermediate"  ~ "Good/Intermediate",
      IMDC == "99"   ~ NA_character_,
      TRUE                  ~ NA_character_
    ),
    IMDC = factor(
      IMDC,
      levels = c("Poor", "Good/Intermediate")
    )
  )

#remove some useless columns
dataDF$sex <- NULL

# rename some of the columns for plotting
names(dataDF)[names(dataDF) == 'age_1st_treat'] <- 'Age'
names(dataDF)[names(dataDF) == 'regime_correct'] <- 'CPI Regimen'
names(dataDF)[names(dataDF) == 'PS'] <- 'ECOG Performance Status'
# names(dataDF)[names(dataDF) == 'stage_2'] <- 'AJCC 8th stage'
names(dataDF)[names(dataDF) == 'Brain_metastases'] <- 'Brain metastases'
# names(dataDF)[names(dataDF) == 'braf_correct'] <- 'BRAF mutation'
names(dataDF)[names(dataDF) == 'bor'] <- 'Objective response'
names(dataDF)[names(dataDF) == 'line_correct'] <- 'Treatment line'
names(dataDF)[names(dataDF) == 'clear_cell'] <- 'RCC subtype'
names(dataDF)[names(dataDF) == 'sarcomatoid'] <- 'Sarcomatoid subtype'
names(dataDF)[names(dataDF) == 'smoking'] <- 'Smoking status'
features <- c("Sex", "Age", "CPI Regimen", "Treatment line", "RCC subtype", "Sarcomatoid subtype","IMDC","ECOG Performance Status", "Objective response")

write.csv(dataDF, file = "RCC_data_polished.csv", sep = ",", append = F, row.names = F)

#Administrative censoring at 60 months
dataDF$censor_time_OS   <- pmin(dataDF$OS_days/30, 60)
dataDF$censor_status_OS <- ifelse(dataDF$OS_days/30 > 60, 0, dataDF$Dead)
#Administrative censoring at 60 months
dataDF$censor_time_PFS   <- pmin(dataDF$PFS_days/30, 60)
dataDF$censor_status_PFS <- ifelse(dataDF$PFS_days/30 > 60, 0, dataDF$Progressed)
```

``` r
descr_features <- paste0("`", c(features, "CPI Regimen"), "`", collapse = " + ")
descr_formula <- as.formula(paste0("~", descr_features))

Baseline_Characteristics_table_Partial_patients <- descrTable(descr_formula , data = dataDF)
```

# OS analysis

\##Multiple Univariate analysis - Including missing data

``` r
surv_time <- "censor_time_OS"
surv_status <- "censor_status_OS"
surv_time_label <- "OS"
# Find the best reference for each feature
best_refs <- character(length(features))
names(best_refs) <- features

for (i in seq_along(features)) {
  var <- features[i]

  # Skip this iteration if this column isn't a factor
  if (! is.factor(dataDF[[var]])) {
    next
  }

  vals <- na.omit(unique(dataDF[[var]]))
  # ensure factor
  vals <- as.character(vals)

  # storage for max HR per candidate ref
  hr_max <- numeric(length(vals))

  for (j in seq_along(vals)) {
    ref_level <- vals[j]
    df2 <- dataDF
    if(!"ordered" %in% class(df2[[var]])){
      df2[[var]] <- relevel(as.factor(df2[[var]]), ref = ref_level)
    }

    uni <- analyse_survival(
      df2,
      time_status = c(surv_time, surv_status),
      by          = .data[[var]]
    )

    fit   <- uni$coxph    # the coxph object
    coefs <- coef(fit)    # log HRs
    hrs   <- exp(coefs)   # HRs
    hr_max[j] <- if (length(hrs) > 0) max(hrs, na.rm = TRUE) else NA
  }

  # pick the ref that yielded the largest HR
  best_refs[var] <- vals[which.max(hr_max)]
}

# Re‐level dataDF based on those best refs
for (var in features) {
  if (! is.factor(dataDF[[var]])) {
    next
  }
  if(!"ordered" %in% class(dataDF[[var]])){
    dataDF[[var]] <- relevel(
    as.factor(dataDF[[var]]),
    ref = best_refs[var]
  )
  }
}

multiple_uni <- map(features, function(by){analyse_multivariate(dataDF,
                       c(surv_time,surv_status),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = NULL,
                       covariate_label_dict = NULL)})
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `factor.name = map_chr(coefficient_labels, symbol_substring,
    ##   call_symbols)`.
    ## Caused by warning:
    ## ! `as_logical()` is deprecated as of rlang 0.4.0
    ## Please use `vctrs::vec_cast()` instead.
    ## This warning is displayed once every 8 hours.

``` r
max_cis <- max(sapply(
  multiple_uni,
  function(x) max(x$summary$conf.int, na.rm = TRUE))
)
min_cis <- min(sapply(
  multiple_uni,
  function(x) min(x$summary$conf.int, na.rm = TRUE))
)

HR_x_limits <- c(min_cis, max_cis) * c(1/1.5, 1.5)
HR_x_breaks <- round(10^pretty(log10(HR_x_limits), n = 7),digits = 1)

forest_plot(multiple_uni,
            #factor_labeller = covariate_names,
            endpoint_labeller = c(time=surv_time_label),
            # orderer = ~order(HR),
            labels_displayed = c("factor"),#"endpoint",
            ggtheme = ggplot2::theme_bw(base_size = 10),
            #values_displayed = c("HR","CI","p"),
            HR_x_limits = c(min_cis*1.5,max_cis*1.5),
            HR_x_breaks = c(1,2,ceiling(10^pretty(log10(HR_x_limits), n = 7)*2/2))
            #p_lessthan_cutoff = 0.05
            )
```

    ## Warning: `as_list()` is deprecated as of rlang 0.4.0
    ## Please use `vctrs::vec_cast()` instead.
    ## This warning is displayed once every 8 hours.

    ## Warning: `switch_type()` is soft-deprecated as of rlang 0.4.0.
    ## Please use `switch(typeof())` or `switch(my_typeof())` instead.
    ## This warning is displayed once every 8 hours.

![](Survival_analysis_renal_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data-1.png)<!-- -->

``` r
#save this for the comparison at the end with the 18 month landmark

hr_OS_resp <-multiple_uni[[length(multiple_uni)]]$coxph
```

## Multivariate Survival Analysis - only complete cases

``` r
# # To remove both (NAs and empty):
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 237  20

``` r
multiv_feats_sign <- c()
for (i in seq_along(features)){
  coef <- as.data.frame(multiple_uni[[i]]$summary$coefficients)
  if(nrow(coef) >1){
    pval <- min(coef$`Pr(>|z|)`)
  } else {
  pval <- coef$`Pr(>|z|)`
  }
  if(any(pval < 0.2)){
    multiv_feats_sign <- c(multiv_feats_sign, multiple_uni[[i]]$overall$covariates)
  }
}

multiv_feats_sign_formula <- paste0("`", multiv_feats_sign, "`", collapse = " + ")
fmla_sign <- as.formula(paste0(
  "Surv(dataDF_complete[[surv_time]], dataDF_complete[[surv_status]]) ~ ",
  multiv_feats_sign_formula
))

multiv_feats_all<- paste0("`", features, "`", collapse = " + ")
fmla_all <- as.formula(paste0(
  "Surv(dataDF_complete[[surv_time]], dataDF_complete[[surv_status]]) ~ ",
  multiv_feats_all
))

cox_model_all <- coxph(fmla_all, data = dataDF_complete)
cox_model_sign <- coxph(fmla_sign, data = dataDF_complete)

summary(cox_model_all)
```

    ## Call:
    ## coxph(formula = fmla_all, data = dataDF_complete)
    ## 
    ##   n= 237, number of events= 70 
    ## 
    ##                                          coef exp(coef) se(coef)      z
    ## SexFemale                             0.53739   1.71153  0.27460  1.957
    ## Age                                   0.01960   1.01980  0.01376  1.425
    ## `CPI Regimen`Anti-PD1                 0.69547   2.00465  0.43681  1.592
    ## `Treatment line`Other                -0.25881   0.77197  0.44297 -0.584
    ## `RCC subtype`Non clear cell           0.73011   2.07531  0.33712  2.166
    ## `Sarcomatoid subtype`Non sarcomatoid  0.10535   1.11110  0.30251  0.348
    ## IMDCPoor                              0.78535   2.19316  0.25169  3.120
    ## `ECOG Performance Status`PS≥1         0.17059   1.18600  0.26304  0.649
    ## `Objective response`PR                2.61710  13.69600  0.59802  4.376
    ##                                      Pr(>|z|)    
    ## SexFemale                             0.05035 .  
    ## Age                                   0.15415    
    ## `CPI Regimen`Anti-PD1                 0.11135    
    ## `Treatment line`Other                 0.55904    
    ## `RCC subtype`Non clear cell           0.03033 *  
    ## `Sarcomatoid subtype`Non sarcomatoid  0.72766    
    ## IMDCPoor                              0.00181 ** 
    ## `ECOG Performance Status`PS≥1         0.51664    
    ## `Objective response`PR               1.21e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                      exp(coef) exp(-coef) lower .95 upper .95
    ## SexFemale                                1.712    0.58427    0.9992     2.932
    ## Age                                      1.020    0.98059    0.9927     1.048
    ## `CPI Regimen`Anti-PD1                    2.005    0.49884    0.8516     4.719
    ## `Treatment line`Other                    0.772    1.29539    0.3240     1.839
    ## `RCC subtype`Non clear cell              2.075    0.48185    1.0718     4.018
    ## `Sarcomatoid subtype`Non sarcomatoid     1.111    0.90001    0.6141     2.010
    ## IMDCPoor                                 2.193    0.45596    1.3392     3.592
    ## `ECOG Performance Status`PS≥1            1.186    0.84317    0.7083     1.986
    ## `Objective response`PR                  13.696    0.07301    4.2419    44.221
    ## 
    ## Concordance= 0.781  (se = 0.027 )
    ## Likelihood ratio test= 72.87  on 9 df,   p=4e-12
    ## Wald test            = 49.86  on 9 df,   p=1e-07
    ## Score (logrank) test = 65.85  on 9 df,   p=1e-10

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 237, number of events= 70 
    ## 
    ##                                          coef exp(coef) se(coef)      z
    ## SexFemale                             0.53739   1.71153  0.27460  1.957
    ## Age                                   0.01960   1.01980  0.01376  1.425
    ## `CPI Regimen`Anti-PD1                 0.69547   2.00465  0.43681  1.592
    ## `Treatment line`Other                -0.25881   0.77197  0.44297 -0.584
    ## `RCC subtype`Non clear cell           0.73011   2.07531  0.33712  2.166
    ## `Sarcomatoid subtype`Non sarcomatoid  0.10535   1.11110  0.30251  0.348
    ## IMDCPoor                              0.78535   2.19316  0.25169  3.120
    ## `ECOG Performance Status`PS≥1         0.17059   1.18600  0.26304  0.649
    ## `Objective response`PR                2.61710  13.69600  0.59802  4.376
    ##                                      Pr(>|z|)    
    ## SexFemale                             0.05035 .  
    ## Age                                   0.15415    
    ## `CPI Regimen`Anti-PD1                 0.11135    
    ## `Treatment line`Other                 0.55904    
    ## `RCC subtype`Non clear cell           0.03033 *  
    ## `Sarcomatoid subtype`Non sarcomatoid  0.72766    
    ## IMDCPoor                              0.00181 ** 
    ## `ECOG Performance Status`PS≥1         0.51664    
    ## `Objective response`PR               1.21e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                      exp(coef) exp(-coef) lower .95 upper .95
    ## SexFemale                                1.712    0.58427    0.9992     2.932
    ## Age                                      1.020    0.98059    0.9927     1.048
    ## `CPI Regimen`Anti-PD1                    2.005    0.49884    0.8516     4.719
    ## `Treatment line`Other                    0.772    1.29539    0.3240     1.839
    ## `RCC subtype`Non clear cell              2.075    0.48185    1.0718     4.018
    ## `Sarcomatoid subtype`Non sarcomatoid     1.111    0.90001    0.6141     2.010
    ## IMDCPoor                                 2.193    0.45596    1.3392     3.592
    ## `ECOG Performance Status`PS≥1            1.186    0.84317    0.7083     1.986
    ## `Objective response`PR                  13.696    0.07301    4.2419    44.221
    ## 
    ## Concordance= 0.781  (se = 0.027 )
    ## Likelihood ratio test= 72.87  on 9 df,   p=4e-12
    ## Wald test            = 49.86  on 9 df,   p=1e-07
    ## Score (logrank) test = 65.85  on 9 df,   p=1e-10

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_renal_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_renal_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases-2.png)<!-- -->

# Combined Univariable vs Multivariable Forest Plots - OS

``` r
# Prepare empty list to collect rows
rows <- vector("list", 0)

# Summaries of the multivariable model
s_multi    <- summary(cox_model_sign)
coef_multi <- s_multi$coefficients      # matrix with exp(coef) and p
ci_multi   <- s_multi$conf.int          # matrix with lower/upper .95

# Loop over each univariable result
for (i in seq_along(multiple_uni)) {
  u    <- multiple_uni[[i]]
  cov  <- u$overall$covariates           # the covariate name
  fit  <- u$coxph
  s_u  <- summary(fit)
  coef_u <- s_u$coefficients             # matrix with exp(coef) etc.
  ci_u   <- s_u$conf.int
  if(u$summaryAsFrame$factor.value[1] == "<continuous>"){
    reference <- "NA"
    levs <- "Fitted as continuous"
  } else {
    reference <- levels(as.factor(u$data[,ncol(u$data)]))[1]
    levs <- setdiff(unique(u$data[,ncol(u$data)]), reference)
  }
  # for each row (level) in the univariable fit:
  for (j in seq_len(nrow(coef_u))) {
    HRu  <- coef_u[j, "exp(coef)"]
    LUu  <- ci_u[j,   "lower .95"]
    UCu  <- ci_u[j,   "upper .95"]
    pu   <- coef_u[j, "Pr(>|z|)"]
    
    if (cov %in% multiv_feats_sign) {
      HRm <- coef_multi[grep(cov, rownames(coef_multi))[j], "exp(coef)"]
      LUm <- ci_multi[ grep(cov, rownames(coef_multi))[j], "lower .95"]
      UCm <- ci_multi[ grep(cov, rownames(coef_multi))[j], "upper .95"]
      pm  <- coef_multi[grep(cov, rownames(coef_multi))[j], "Pr(>|z|)"]
    } else {
      HRm <- NA; LUm <- NA; UCm <- NA; pm <- NA
    }
    
    rows[[length(rows) + 1]] <- data.frame(
      Covariate       = cov,
      Level           = levs[j],
      Reference       = reference,
      HR_uni          = HRu,
      CI_lo_uni       = LUu,
      CI_hi_uni       = UCu,
      p_uni           = pu,
      HR_multi        = HRm,
      CI_lo_multi     = LUm,
      CI_hi_multi     = UCm,
      p_multi         = pm,
      stringsAsFactors = FALSE
    )
  }
}

# Bind into one data.frame
results_df <- do.call(rbind, rows)

# Format up to three digits
results_df <- results_df %>%
  mutate(across(starts_with("HR_") | starts_with("CI_") | matches("^p_"),
                ~ signif(.x, 3)))

results_df$CI_lo_uni <- paste0(results_df$CI_lo_uni, "-", results_df$CI_hi_uni)
results_df$CI_hi_uni <- NULL
results_df$CI_lo_multi <- paste0(results_df$CI_lo_multi, "-", results_df$CI_hi_multi)
results_df$CI_hi_multi <- NULL
colnames(results_df) <- c("Covariate", "Level", "Reference" , "Univariable HR", "Univariable CI", "Univariable P-value", "Multivariable HR", "Multivariable CI", "Multivariable P-value")

last_col <- ncol(results_df)
results_df <- results_df %>%
  mutate(across(
    .cols = -last_col,
    .fns  = ~ ifelse(is.na(.) | grepl("NA", ., fixed = TRUE), "", .)
  )) %>%
  mutate(across(
    .cols = last_col,
    .fns  = ~ ifelse(is.na(.), "Excluded from multivariable analysis", as.character(.))
  ))

results_df <- results_df %>%
  group_by(Covariate) %>%
  mutate(Covariate = ifelse(row_number() == 1, Covariate, "")) %>%
  ungroup()

n_cols <- ncol(results_df)

# align: left for first two, center for the rest
aligns <- c("l", "l", rep("c", n_cols - 2))

results_df %>%
  kable(
    format = "html",
    caption = "Univariable vs Multivariable Cox Results - OS",
    align   = aligns,
    escape  = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width        = FALSE
  ) %>%
  # make header row bold and centered
  row_spec(
    0,
    bold  = TRUE,
    align = "center"
  ) %>%
  # add vertical lines right of col 3 and col 6
  column_spec(3, border_right = TRUE) %>%
  column_spec(6, border_right = TRUE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Univariable vs Multivariable Cox Results - OS
</caption>

<thead>

<tr>

<th style="text-align:left;font-weight: bold;text-align: center;">

Covariate
</th>

<th style="text-align:left;font-weight: bold;text-align: center;">

Level
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Reference
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Univariable HR
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Univariable CI
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Univariable P-value
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Multivariable HR
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Multivariable CI
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Multivariable P-value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Sex
</td>

<td style="text-align:left;">

Female
</td>

<td style="text-align:center;border-right:1px solid;">

Male
</td>

<td style="text-align:center;">

1.86
</td>

<td style="text-align:center;">

1.16-2.97
</td>

<td style="text-align:center;border-right:1px solid;">

1.01e-02
</td>

<td style="text-align:center;">

1.710
</td>

<td style="text-align:center;">

0.999-2.93
</td>

<td style="text-align:center;">

0.0504
</td>

</tr>

<tr>

<td style="text-align:left;">

Age
</td>

<td style="text-align:left;">

Fitted as continuous
</td>

<td style="text-align:center;border-right:1px solid;">

</td>

<td style="text-align:center;">

1.03
</td>

<td style="text-align:center;">

1-1.05
</td>

<td style="text-align:center;border-right:1px solid;">

2.62e-02
</td>

<td style="text-align:center;">

1.020
</td>

<td style="text-align:center;">

0.993-1.05
</td>

<td style="text-align:center;">

0.154
</td>

</tr>

<tr>

<td style="text-align:left;">

CPI Regimen
</td>

<td style="text-align:left;">

Anti-PD1
</td>

<td style="text-align:center;border-right:1px solid;">

Anti-PD1+Anti-CTLA4
</td>

<td style="text-align:center;">

1.65
</td>

<td style="text-align:center;">

1.02-2.65
</td>

<td style="text-align:center;border-right:1px solid;">

4.04e-02
</td>

<td style="text-align:center;">

2.000
</td>

<td style="text-align:center;">

0.852-4.72
</td>

<td style="text-align:center;">

0.111
</td>

</tr>

<tr>

<td style="text-align:left;">

Treatment line
</td>

<td style="text-align:left;">

Other
</td>

<td style="text-align:center;border-right:1px solid;">

First
</td>

<td style="text-align:center;">

1.47
</td>

<td style="text-align:center;">

0.919-2.34
</td>

<td style="text-align:center;border-right:1px solid;">

1.09e-01
</td>

<td style="text-align:center;">

0.772
</td>

<td style="text-align:center;">

0.324-1.84
</td>

<td style="text-align:center;">

0.559
</td>

</tr>

<tr>

<td style="text-align:left;">

RCC subtype
</td>

<td style="text-align:left;">

Non clear cell
</td>

<td style="text-align:center;border-right:1px solid;">

Clear cell
</td>

<td style="text-align:center;">

1.74
</td>

<td style="text-align:center;">

0.975-3.1
</td>

<td style="text-align:center;border-right:1px solid;">

6.09e-02
</td>

<td style="text-align:center;">

2.080
</td>

<td style="text-align:center;">

1.07-4.02
</td>

<td style="text-align:center;">

0.0303
</td>

</tr>

<tr>

<td style="text-align:left;">

Sarcomatoid subtype
</td>

<td style="text-align:left;">

Non sarcomatoid
</td>

<td style="text-align:center;border-right:1px solid;">

Sarcomatoid
</td>

<td style="text-align:center;">

1.46
</td>

<td style="text-align:center;">

0.873-2.46
</td>

<td style="text-align:center;border-right:1px solid;">

1.48e-01
</td>

<td style="text-align:center;">

1.110
</td>

<td style="text-align:center;">

0.614-2.01
</td>

<td style="text-align:center;">

0.728
</td>

</tr>

<tr>

<td style="text-align:left;">

IMDC
</td>

<td style="text-align:left;">

Poor
</td>

<td style="text-align:center;border-right:1px solid;">

Good/Intermediate
</td>

<td style="text-align:center;">

2.28
</td>

<td style="text-align:center;">

1.46-3.57
</td>

<td style="text-align:center;border-right:1px solid;">

3.05e-04
</td>

<td style="text-align:center;">

2.190
</td>

<td style="text-align:center;">

1.34-3.59
</td>

<td style="text-align:center;">

0.00181
</td>

</tr>

<tr>

<td style="text-align:left;">

ECOG Performance Status
</td>

<td style="text-align:left;">

PS≥1
</td>

<td style="text-align:center;border-right:1px solid;">

PS=0
</td>

<td style="text-align:center;">

2.08
</td>

<td style="text-align:center;">

1.31-3.31
</td>

<td style="text-align:center;border-right:1px solid;">

1.95e-03
</td>

<td style="text-align:center;">

1.190
</td>

<td style="text-align:center;">

0.708-1.99
</td>

<td style="text-align:center;">

0.517
</td>

</tr>

<tr>

<td style="text-align:left;">

Objective response
</td>

<td style="text-align:left;">

PR
</td>

<td style="text-align:center;border-right:1px solid;">

CR
</td>

<td style="text-align:center;">

14.40
</td>

<td style="text-align:center;">

4.53-45.6
</td>

<td style="text-align:center;border-right:1px solid;">

6.10e-06
</td>

<td style="text-align:center;">

13.700
</td>

<td style="text-align:center;">

4.24-44.2
</td>

<td style="text-align:center;">

1.21e-05
</td>

</tr>

</tbody>

</table>

# PFS analysis

## Multiple Univariate analysis - Incomplete complete cases

``` r
surv_time <- "censor_time_PFS"
surv_status <- "censor_status_PFS"
surv_time_label <- "PFS"
# Find the best reference for each feature
best_refs <- character(length(features))
names(best_refs) <- features

for (i in seq_along(features)) {
  var <- features[i]

  # Skip this iteration if this column isn't a factor
  if (! is.factor(dataDF[[var]])) {
    next
  }

  vals <- na.omit(unique(dataDF[[var]]))
  # ensure factor
  vals <- as.character(vals)

  # storage for max HR per candidate ref
  hr_max <- numeric(length(vals))

  for (j in seq_along(vals)) {
    ref_level <- vals[j]
    df2 <- dataDF
    if(!"ordered" %in% class(df2[[var]])){
      df2[[var]] <- relevel(as.factor(df2[[var]]), ref = ref_level)
    }

    uni <- analyse_survival(
      df2,
      time_status = c(surv_time, surv_status),
      by          = .data[[var]]
    )

    fit   <- uni$coxph    # the coxph object
    coefs <- coef(fit)    # log HRs
    hrs   <- exp(coefs)   # HRs
    hr_max[j] <- if (length(hrs) > 0) max(hrs, na.rm = TRUE) else NA
  }

  # pick the ref that yielded the largest HR
  best_refs[var] <- vals[which.max(hr_max)]
}

# Re‐level dataDF based on those best refs
for (var in features) {
  if (! is.factor(dataDF[[var]])) {
    next
  }
  if(!"ordered" %in% class(dataDF[[var]])){
    dataDF[[var]] <- relevel(
    as.factor(dataDF[[var]]),
    ref = best_refs[var]
  )
  }
}

multiple_uni <- map(features, function(by){analyse_multivariate(dataDF,
                       c(surv_time,surv_status),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = NULL,
                       covariate_label_dict = NULL)})
max_cis <- max(sapply(
  multiple_uni,
  function(x) max(x$summary$conf.int, na.rm = TRUE))
)
min_cis <- min(sapply(
  multiple_uni,
  function(x) min(x$summary$conf.int, na.rm = TRUE))
)

HR_x_limits <- c(min_cis, max_cis) * c(1/1.5, 1.5)
HR_x_breaks <- round(10^pretty(log10(HR_x_limits), n = 7),digits = 1)

forest_plot(multiple_uni,
            #factor_labeller = covariate_names,
            endpoint_labeller = c(time=surv_time_label),
            # orderer = ~order(HR),
            labels_displayed = c("factor"),#"endpoint",
            ggtheme = ggplot2::theme_bw(base_size = 10),
            #values_displayed = c("HR","CI","p"),
            HR_x_limits = c(min_cis*1.5,max_cis*1.5),
            HR_x_breaks = c(1,2,ceiling(10^pretty(log10(HR_x_limits), n = 7)*2/2))
            #p_lessthan_cutoff = 0.05
            )
```

![](Survival_analysis_renal_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20PFS-1.png)<!-- -->

# Multivariate Survival Analysis

## All responders - PFS

``` r
multiv_feats_sign <- c()
for (i in seq_along(features)){
  coef <- as.data.frame(multiple_uni[[i]]$summary$coefficients)
  if(nrow(coef) >1){
    pval <- min(coef$`Pr(>|z|)`)
  } else {
  pval <- coef$`Pr(>|z|)`
  }
  if(any(pval < 0.2)){
    multiv_feats_sign <- c(multiv_feats_sign, multiple_uni[[i]]$overall$covariates)
  }
}

multiv_feats_sign_formula <- paste0("`", multiv_feats_sign, "`", collapse = " + ")
fmla_sign <- as.formula(paste0(
  "Surv(dataDF_complete[[surv_time]], dataDF_complete[[surv_status]]) ~ ",
  multiv_feats_sign_formula
))

multiv_feats_all<- paste0("`", features, "`", collapse = " + ")
fmla_all <- as.formula(paste0(
  "Surv(dataDF_complete[[surv_time]], dataDF_complete[[surv_status]]) ~ ",
  multiv_feats_all
))

cox_model_all <- coxph(fmla_all, data = dataDF_complete)
cox_model_sign <- coxph(fmla_sign, data = dataDF_complete)

summary(cox_model_all)
```

    ## Call:
    ## coxph(formula = fmla_all, data = dataDF_complete)
    ## 
    ##   n= 237, number of events= 141 
    ## 
    ##                                           coef exp(coef)  se(coef)      z
    ## SexFemale                             0.098200  1.103183  0.205726  0.477
    ## Age                                   0.001019  1.001020  0.009178  0.111
    ## `CPI Regimen`Anti-PD1                -0.025460  0.974861  0.314368 -0.081
    ## `Treatment line`Other                 0.232608  1.261886  0.303809  0.766
    ## `RCC subtype`Non clear cell           0.399935  1.491728  0.243357  1.643
    ## `Sarcomatoid subtype`Non sarcomatoid  0.120939  1.128556  0.204963  0.590
    ## IMDCPoor                              0.152932  1.165246  0.187496  0.816
    ## `ECOG Performance Status`PS≥1         0.019422  1.019612  0.179454  0.108
    ## `Objective response`PR                1.551614  4.719080  0.249898  6.209
    ##                                      Pr(>|z|)    
    ## SexFemale                               0.633    
    ## Age                                     0.912    
    ## `CPI Regimen`Anti-PD1                   0.935    
    ## `Treatment line`Other                   0.444    
    ## `RCC subtype`Non clear cell             0.100    
    ## `Sarcomatoid subtype`Non sarcomatoid    0.555    
    ## IMDCPoor                                0.415    
    ## `ECOG Performance Status`PS≥1           0.914    
    ## `Objective response`PR               5.33e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                      exp(coef) exp(-coef) lower .95 upper .95
    ## SexFemale                               1.1032     0.9065    0.7371     1.651
    ## Age                                     1.0010     0.9990    0.9832     1.019
    ## `CPI Regimen`Anti-PD1                   0.9749     1.0258    0.5264     1.805
    ## `Treatment line`Other                   1.2619     0.7925    0.6957     2.289
    ## `RCC subtype`Non clear cell             1.4917     0.6704    0.9259     2.403
    ## `Sarcomatoid subtype`Non sarcomatoid    1.1286     0.8861    0.7552     1.687
    ## IMDCPoor                                1.1652     0.8582    0.8069     1.683
    ## `ECOG Performance Status`PS≥1           1.0196     0.9808    0.7173     1.449
    ## `Objective response`PR                  4.7191     0.2119    2.8916     7.701
    ## 
    ## Concordance= 0.688  (se = 0.021 )
    ## Likelihood ratio test= 64.97  on 9 df,   p=1e-10
    ## Wald test            = 49.53  on 9 df,   p=1e-07
    ## Score (logrank) test = 58.7  on 9 df,   p=2e-09

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 237, number of events= 141 
    ## 
    ##                                         coef exp(coef) se(coef)     z Pr(>|z|)
    ## `Treatment line`Other                0.21117   1.23513  0.20237 1.044   0.2967
    ## `RCC subtype`Non clear cell          0.40555   1.50013  0.23898 1.697   0.0897
    ## `Sarcomatoid subtype`Non sarcomatoid 0.10677   1.11268  0.20379 0.524   0.6003
    ## `ECOG Performance Status`PS≥1        0.05978   1.06161  0.17369 0.344   0.7307
    ## `Objective response`PR               1.54601   4.69272  0.24764 6.243 4.29e-10
    ##                                         
    ## `Treatment line`Other                   
    ## `RCC subtype`Non clear cell          .  
    ## `Sarcomatoid subtype`Non sarcomatoid    
    ## `ECOG Performance Status`PS≥1           
    ## `Objective response`PR               ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                      exp(coef) exp(-coef) lower .95 upper .95
    ## `Treatment line`Other                    1.235     0.8096    0.8307     1.836
    ## `RCC subtype`Non clear cell              1.500     0.6666    0.9391     2.396
    ## `Sarcomatoid subtype`Non sarcomatoid     1.113     0.8987    0.7463     1.659
    ## `ECOG Performance Status`PS≥1            1.062     0.9420    0.7553     1.492
    ## `Objective response`PR                   4.693     0.2131    2.8883     7.625
    ## 
    ## Concordance= 0.69  (se = 0.021 )
    ## Likelihood ratio test= 64.01  on 5 df,   p=2e-12
    ## Wald test            = 48.71  on 5 df,   p=3e-09
    ## Score (logrank) test = 57.88  on 5 df,   p=3e-11

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_renal_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_renal_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

# Combined Univariable vs Multivariable Forest Plots - PFS

``` r
# Prepare empty list to collect rows
rows <- vector("list", 0)

# Summaries of the multivariable model
s_multi    <- summary(cox_model_sign)
coef_multi <- s_multi$coefficients      # matrix with exp(coef) and p
ci_multi   <- s_multi$conf.int          # matrix with lower/upper .95

# Loop over each univariable result
for (i in seq_along(multiple_uni)) {
  u    <- multiple_uni[[i]]
  cov  <- u$overall$covariates           # the covariate name
  fit  <- u$coxph
  s_u  <- summary(fit)
  coef_u <- s_u$coefficients             # matrix with exp(coef) etc.
  ci_u   <- s_u$conf.int
  if(u$summaryAsFrame$factor.value[1] == "<continuous>"){
    reference <- "NA"
    levs <- "Fitted as continuous"
  } else {
    reference <- levels(as.factor(u$data[,ncol(u$data)]))[1]
    levs <- setdiff(unique(u$data[,ncol(u$data)]), reference)
  }
  # for each row (level) in the univariable fit:
  for (j in seq_len(nrow(coef_u))) {
    HRu  <- coef_u[j, "exp(coef)"]
    LUu  <- ci_u[j,   "lower .95"]
    UCu  <- ci_u[j,   "upper .95"]
    pu   <- coef_u[j, "Pr(>|z|)"]
    
    if (cov %in% multiv_feats_sign) {
      HRm <- coef_multi[grep(cov, rownames(coef_multi))[j], "exp(coef)"]
      LUm <- ci_multi[ grep(cov, rownames(coef_multi))[j], "lower .95"]
      UCm <- ci_multi[ grep(cov, rownames(coef_multi))[j], "upper .95"]
      pm  <- coef_multi[grep(cov, rownames(coef_multi))[j], "Pr(>|z|)"]
    } else {
      HRm <- NA; LUm <- NA; UCm <- NA; pm <- NA
    }
    
    rows[[length(rows) + 1]] <- data.frame(
      Covariate       = cov,
      Level           = levs[j],
      Reference       = reference,
      HR_uni          = HRu,
      CI_lo_uni       = LUu,
      CI_hi_uni       = UCu,
      p_uni           = pu,
      HR_multi        = HRm,
      CI_lo_multi     = LUm,
      CI_hi_multi     = UCm,
      p_multi         = pm,
      stringsAsFactors = FALSE
    )
  }
}

# Bind into one data.frame
results_df <- do.call(rbind, rows)

# Format up to three digits 
results_df <- results_df %>%
  mutate(across(starts_with("HR_") | starts_with("CI_") | matches("^p_"),
                ~ signif(.x, 3)))

results_df$CI_lo_uni <- paste0(results_df$CI_lo_uni, "-", results_df$CI_hi_uni)
results_df$CI_hi_uni <- NULL
results_df$CI_lo_multi <- paste0(results_df$CI_lo_multi, "-", results_df$CI_hi_multi)
results_df$CI_hi_multi <- NULL
colnames(results_df) <- c("Covariate", "Level", "Reference" , "Univariable HR", "Univariable CI", "Univariable P-value", "Multivariable HR", "Multivariable CI", "Multivariable P-value")

last_col <- ncol(results_df)
results_df <- results_df %>%
  mutate(across(
    .cols = -last_col,
    .fns  = ~ ifelse(is.na(.) | grepl("NA", ., fixed = TRUE), "", .)
  )) %>%
  mutate(across(
    .cols = last_col,
    .fns  = ~ ifelse(is.na(.), "Excluded from multivariable analysis", as.character(.))
  ))

results_df$`Univariable P-value` <- round(results_df$`Univariable P-value`, digits = 3)
# Now print
knitr::kable(
  results_df,
  caption = "Univariable vs Multivariable Cox Results",
  row.names = FALSE
)
```

| Covariate | Level | Reference | Univariable HR | Univariable CI | Univariable P-value | Multivariable HR | Multivariable CI | Multivariable P-value |
|:---|:---|:---|---:|:---|---:|:---|:---|:---|
| Sex | Female | Male | 1.12 | 0.773-1.63 | 0.544 |  |  | Excluded from multivariable analysis |
| Age | Fitted as continuous |  | 1.01 | 0.992-1.03 | 0.318 |  |  | Excluded from multivariable analysis |
| CPI Regimen | Anti-PD1 | Anti-PD1+Anti-CTLA4 | 1.14 | 0.795-1.64 | 0.471 |  |  | Excluded from multivariable analysis |
| Treatment line | Other | First | 1.33 | 0.947-1.87 | 0.100 | 1.24 | 0.831-1.84 | 0.297 |
| RCC subtype | Non clear cell | Clear cell | 1.49 | 0.972-2.29 | 0.067 | 1.5 | 0.939-2.4 | 0.0897 |
| Sarcomatoid subtype | Non sarcomatoid | Sarcomatoid | 1.49 | 1.03-2.15 | 0.034 | 1.11 | 0.746-1.66 | 0.6 |
| IMDC | Poor | Good/Intermediate | 1.21 | 0.861-1.7 | 0.270 |  |  | Excluded from multivariable analysis |
| ECOG Performance Status | PS≥1 | PS=0 | 1.32 | 0.957-1.81 | 0.091 | 1.06 | 0.755-1.49 | 0.731 |
| Objective response | PR | CR | 4.70 | 2.96-7.47 | 0.000 | 4.69 | 2.89-7.62 | 4.29e-10 |

Univariable vs Multivariable Cox Results

``` r
n_cols <- ncol(results_df)

# align: left for first two, center for the rest
aligns <- c("l", "l", rep("c", n_cols - 2))

results_df %>%
  kable(
    format = "html",
    caption = "Univariable vs Multivariable Cox Results - PFS",
    align   = aligns,
    escape  = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width        = FALSE
  ) %>%
  # make header row bold and centered
  row_spec(
    0,
    bold  = TRUE,
    align = "center"
  ) %>%
  # add vertical lines right of col 3 and col 6
  column_spec(3, border_right = TRUE) %>%
  column_spec(6, border_right = TRUE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Univariable vs Multivariable Cox Results - PFS
</caption>

<thead>

<tr>

<th style="text-align:left;font-weight: bold;text-align: center;">

Covariate
</th>

<th style="text-align:left;font-weight: bold;text-align: center;">

Level
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Reference
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Univariable HR
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Univariable CI
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Univariable P-value
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Multivariable HR
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Multivariable CI
</th>

<th style="text-align:center;font-weight: bold;text-align: center;">

Multivariable P-value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Sex
</td>

<td style="text-align:left;">

Female
</td>

<td style="text-align:center;border-right:1px solid;">

Male
</td>

<td style="text-align:center;">

1.12
</td>

<td style="text-align:center;">

0.773-1.63
</td>

<td style="text-align:center;border-right:1px solid;">

0.544
</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

Excluded from multivariable analysis
</td>

</tr>

<tr>

<td style="text-align:left;">

Age
</td>

<td style="text-align:left;">

Fitted as continuous
</td>

<td style="text-align:center;border-right:1px solid;">

</td>

<td style="text-align:center;">

1.01
</td>

<td style="text-align:center;">

0.992-1.03
</td>

<td style="text-align:center;border-right:1px solid;">

0.318
</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

Excluded from multivariable analysis
</td>

</tr>

<tr>

<td style="text-align:left;">

CPI Regimen
</td>

<td style="text-align:left;">

Anti-PD1
</td>

<td style="text-align:center;border-right:1px solid;">

Anti-PD1+Anti-CTLA4
</td>

<td style="text-align:center;">

1.14
</td>

<td style="text-align:center;">

0.795-1.64
</td>

<td style="text-align:center;border-right:1px solid;">

0.471
</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

Excluded from multivariable analysis
</td>

</tr>

<tr>

<td style="text-align:left;">

Treatment line
</td>

<td style="text-align:left;">

Other
</td>

<td style="text-align:center;border-right:1px solid;">

First
</td>

<td style="text-align:center;">

1.33
</td>

<td style="text-align:center;">

0.947-1.87
</td>

<td style="text-align:center;border-right:1px solid;">

0.100
</td>

<td style="text-align:center;">

1.24
</td>

<td style="text-align:center;">

0.831-1.84
</td>

<td style="text-align:center;">

0.297
</td>

</tr>

<tr>

<td style="text-align:left;">

RCC subtype
</td>

<td style="text-align:left;">

Non clear cell
</td>

<td style="text-align:center;border-right:1px solid;">

Clear cell
</td>

<td style="text-align:center;">

1.49
</td>

<td style="text-align:center;">

0.972-2.29
</td>

<td style="text-align:center;border-right:1px solid;">

0.067
</td>

<td style="text-align:center;">

1.5
</td>

<td style="text-align:center;">

0.939-2.4
</td>

<td style="text-align:center;">

0.0897
</td>

</tr>

<tr>

<td style="text-align:left;">

Sarcomatoid subtype
</td>

<td style="text-align:left;">

Non sarcomatoid
</td>

<td style="text-align:center;border-right:1px solid;">

Sarcomatoid
</td>

<td style="text-align:center;">

1.49
</td>

<td style="text-align:center;">

1.03-2.15
</td>

<td style="text-align:center;border-right:1px solid;">

0.034
</td>

<td style="text-align:center;">

1.11
</td>

<td style="text-align:center;">

0.746-1.66
</td>

<td style="text-align:center;">

0.6
</td>

</tr>

<tr>

<td style="text-align:left;">

IMDC
</td>

<td style="text-align:left;">

Poor
</td>

<td style="text-align:center;border-right:1px solid;">

Good/Intermediate
</td>

<td style="text-align:center;">

1.21
</td>

<td style="text-align:center;">

0.861-1.7
</td>

<td style="text-align:center;border-right:1px solid;">

0.270
</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

Excluded from multivariable analysis
</td>

</tr>

<tr>

<td style="text-align:left;">

ECOG Performance Status
</td>

<td style="text-align:left;">

PS≥1
</td>

<td style="text-align:center;border-right:1px solid;">

PS=0
</td>

<td style="text-align:center;">

1.32
</td>

<td style="text-align:center;">

0.957-1.81
</td>

<td style="text-align:center;border-right:1px solid;">

0.091
</td>

<td style="text-align:center;">

1.06
</td>

<td style="text-align:center;">

0.755-1.49
</td>

<td style="text-align:center;">

0.731
</td>

</tr>

<tr>

<td style="text-align:left;">

Objective response
</td>

<td style="text-align:left;">

PR
</td>

<td style="text-align:center;border-right:1px solid;">

CR
</td>

<td style="text-align:center;">

4.70
</td>

<td style="text-align:center;">

2.96-7.47
</td>

<td style="text-align:center;border-right:1px solid;">

0.000
</td>

<td style="text-align:center;">

4.69
</td>

<td style="text-align:center;">

2.89-7.62
</td>

<td style="text-align:center;">

4.29e-10
</td>

</tr>

</tbody>

</table>

``` r
wb <- createWorkbook()
addWorksheet(wb, "Supp_Table_Cox")

# Write data
writeData(wb, sheet = "Supp_Table_Cox", x = results_df, startRow = 1, startCol = 1)

# Compute dimensions
n_rows <- nrow(results_df) + 1  # +1 for header
n_cols <- ncol(results_df)

# Define styles
headerStyle <- createStyle(
  fontSize        = 12,
  textDecoration  = "bold",
  halign          = "center",
  border          = "bottom",
  borderStyle     = "medium"
)
evenRowStyle <- createStyle(fgFill = "#F2F2F2")
sepBorderStyle <- createStyle(border = "right", borderStyle = "thin")

# Apply alternating shading
dataRows <- 2:n_rows
evenRows <- dataRows[seq(1, length(dataRows), by = 2)]

# Add vertical separators after col 3 and 6
addStyle(
  wb, "Supp_Table_Cox",
  style = sepBorderStyle,
  rows  = 1:n_rows,
  cols  = 3,
  gridExpand = TRUE
)
addStyle(
  wb, "Supp_Table_Cox",
  style = sepBorderStyle,
  rows  = 1:n_rows,
  cols  = 6,
  gridExpand = TRUE
)

#alternate shading
addStyle(
  wb, "Supp_Table_Cox",
  style = evenRowStyle,
  rows  = evenRows,
  cols  = 1:n_cols,
  gridExpand = TRUE
)

# Apply header style
addStyle(
  wb, "Supp_Table_Cox",
  style = headerStyle,
  rows  = 1, cols = 1:n_cols,
  gridExpand = TRUE
)

# 7) Autofit column widths & freeze header
setColWidths(wb, "Supp_Table_Cox", cols = 1:n_cols, widths = "auto")
freezePane(wb, "Supp_Table_Cox", firstRow = TRUE)

# 8) Save file
saveWorkbook(
  wb,
  file = "E:/PhD_projects/Realworld/Scripts/Acquired_resistance_RW/Tables/Supplementary_Table_3.xlsx",
  overwrite = TRUE
)
```

# Five-Year Survival Estimates by Response

``` r
library(survival)
library(dplyr)

# the two response groups
responses <- c("CR", "PR")

# helper to get survival ± CI at 60 months
get5yr_ci <- function(fit) {
  s <- summary(fit, times = 60)
  if (length(s$surv) == 0) {
    s_all   <- summary(fit)
    surv    <- tail(s_all$surv, 1)
    lower   <- tail(s_all$lower, 1)
    upper   <- tail(s_all$upper, 1)
  } else {
    surv  <- s$surv
    lower <- s$lower
    upper <- s$upper
  }
  data.frame(
    Estimate = surv,
    Lower95  = lower,
    Upper95  = upper
  )
}

# loop over response levels
results <- lapply(responses, function(resp) {
  d <- filter(dataDF, `Objective response` == resp)
  fit_OS  <- survfit(Surv(censor_time_OS,  censor_status_OS)  ~ 1, data = d)
  fit_PFS <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1, data = d)

  os_ci  <- get5yr_ci(fit_OS)
  pfs_ci <- get5yr_ci(fit_PFS)

  data.frame(
    Response       = resp,
    OS_5yr         = os_ci$Estimate,
    OS_5yr_Lower   = os_ci$Lower95,
    OS_5yr_Upper   = os_ci$Upper95,
    PFS_5yr        = pfs_ci$Estimate,
    PFS_5yr_Lower  = pfs_ci$Lower95,
    PFS_5yr_Upper  = pfs_ci$Upper95
  )
})

five_year_table <- bind_rows(results) %>%
  mutate(across(matches("_5yr"), ~ . * 100, .names = "{.col}_pct")) %>%
  select(Response,
         OS_5yr_pct, OS_5yr_Lower_pct, OS_5yr_Upper_pct,
         PFS_5yr_pct, PFS_5yr_Lower_pct, PFS_5yr_Upper_pct)
print(five_year_table)
```

    ##   Response OS_5yr_pct OS_5yr_Lower_pct OS_5yr_Upper_pct PFS_5yr_pct
    ## 1       CR   95.46956         90.55291         100.0000    68.38697
    ## 2       PR   52.46009         44.60644          61.6965    23.70498
    ##   PFS_5yr_Lower_pct PFS_5yr_Upper_pct
    ## 1          57.54669          81.26927
    ## 2          17.83307          31.51035

``` r
# pull out labels from five_year_table
x_max=60
rownames(five_year_table) <- five_year_table$Response
os_lbls <- five_year_table %>%
  transmute(
    Response,
    os_lbls = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                       Response, OS_5yr_pct, OS_5yr_Lower_pct, OS_5yr_Upper_pct))

pfs_lbls <- five_year_table %>%
  transmute(
    Response,
    pfs_lbls = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                        Response, PFS_5yr_pct, PFS_5yr_Lower_pct, PFS_5yr_Upper_pct))

#adding a syntactically correct name for OR, as it was failing in the Survfit step
dataDF$Objective_response = dataDF$`Objective response`

ypos <- c(CR = five_year_table$OS_5yr_pct[1], PR = five_year_table$OS_5yr_pct[2])

fit_OS <- survfit(
  Surv(censor_time_OS, censor_status_OS) ~ Objective_response,
  data = dataDF
)
names(fit_OS$strata) <- gsub("Objective_response=", "", names(fit_OS$strata))
# OS plot
p_os <- ggsurvfit(fit_OS, size = 1.5) +
    add_censor_mark() +
    add_confidence_interval()+
    scale_ggsurvfit() +
    add_risktable(
      size            = 7,
      theme           = theme_risktable_default(axis.text.y.size = 10,
                                                 plot.title.size  = 20),
      risktable_stats = "{n.risk}",  # ({cum.event}) removed because of space
      stats_label = "Number at risk"
    ) +
    add_risktable_strata_symbol(symbol = "•", size = 20)+
    labs(
      x        = "Months after treatment initiation",
      y        = "PFS (%)",
      title = "RCC - Overall Survival"
    ) +
    scale_x_continuous(breaks = seq(0, x_max, by = 12), limits = c(0, x_max)) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 25),
      axis.title.x    = element_text(size = 20),
      axis.title.y    = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text.x     = element_text(size = 15),
      axis.text.y     = element_text(size = 15),
      legend.position = "bottom",
      legend.direction= "horizontal",
      legend.text     = element_text(size = 18),
      legend.key.size = unit(10, "bigpts"),
      legend.title    = element_blank(),
      plot.margin = unit(c(0,0.2,0,1), 'lines')
    ) +
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#619CFF", "PR"="#F564E3")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#619CFF", "PR"="#F564E3")) +
   annotate(
    "text",
    x     = x_max*0.25,
    y     = ypos[names(ypos)]*0.9/100,
    label = paste0("5 year OS for ",os_lbls$os_lbls),
    size  = 10,
    hjust = 0
    )
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

``` r
# PFS plot
ypos <- c(CR = five_year_table$PFS_5yr_pct[1], PR = five_year_table$PFS_5yr_pct[2])

fit_PFS <- survfit(
  Surv(censor_time_PFS, censor_status_PFS) ~ Objective_response,
  data = dataDF
)
names(fit_PFS$strata) <- gsub("Objective_response=", "", names(fit_PFS$strata))

p_pfs <- ggsurvfit(fit_PFS, size = 1.5) +
    add_censor_mark() +
    add_confidence_interval()+
    scale_ggsurvfit() +
    add_risktable(
      size            = 7,
      theme           = theme_risktable_default(axis.text.y.size = 10,
                                                 plot.title.size  = 20),
      risktable_stats = "{n.risk}",  # ({cum.event}) removed because of space
      stats_label = "Number at risk"
    ) +
    add_risktable_strata_symbol(symbol = "•", size = 20)+
    labs(
      x        = "Months after treatment initiation",
      y        = "PFS (%)",
      title = "RCC - Progression-free Survival"
    ) +
    scale_x_continuous(breaks = seq(0, x_max, by = 12), limits = c(0, x_max)) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 25),
      axis.title.x    = element_text(size = 20),
      axis.title.y    = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text.x     = element_text(size = 15),
      axis.text.y     = element_text(size = 15),
      legend.position = "bottom",
      legend.direction= "horizontal",
      legend.text     = element_text(size = 18),
      legend.key.size = unit(10, "bigpts"),
      legend.title    = element_blank(),
      plot.margin = unit(c(0,0.2,0,1), 'lines')
    ) +
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#619CFF", "PR"="#F564E3")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#619CFF", "PR"="#F564E3")) +
   annotate(
    "text",
    x     = x_max*0.25,
    y     = ypos[names(ypos)]*0.9/100,
    label = paste0("5 year PFS for ",pfs_lbls$pfs_lbls),
    size  = 10,
    hjust = 0
    )
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

``` r
p_pfs
```

![](Survival_analysis_renal_files/figure-gfm/-%205%20year%20KM%20curves-1.png)<!-- -->

``` r
p_os
```

![](Survival_analysis_renal_files/figure-gfm/-%205%20year%20KM%20curves-2.png)<!-- -->
