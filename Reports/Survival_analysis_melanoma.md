Survival analysis melanoma
================
Mario Presti
First created on Feb 2025. Updated on 28 November 2025

- [Introduction](#introduction)
- [Loading Libraries](#loading-libraries)
- [Pre-processing of data](#pre-processing-of-data)
- [OS analysis](#os-analysis)
- [Combined Univariable vs Multivariable Forest Plots -
  OS](#combined-univariable-vs-multivariable-forest-plots---os)
- [PFS analysis](#pfs-analysis)
- [Combined Univariable vs Multivariable Forest Plots -
  PFS](#combined-univariable-vs-multivariable-forest-plots---pfs)
- [Melanoma specific survival](#melanoma-specific-survival)
- [Combined Univariable vs Multivariable Forest Plots -
  DSS](#combined-univariable-vs-multivariable-forest-plots---dss)
  - [Numeric table with 5-year OS, PFS,
    DSS](#numeric-table-with-5-year-os-pfs-dss)
- [Five-Year Survival Estimates by
  Response](#five-year-survival-estimates-by-response)
  - [18 Month landmark analysis](#18-month-landmark-analysis)

# Introduction

Multivariate Survival Analysis on Melanoma database

# Loading Libraries

``` r
# Load required libraries
library(pacman)
pacman::p_load(data.table,stringr,tidyverse,ggplot2,MatchIt,survival,survminer,survey,compareGroups, forestmodel, forestplot, openxlsx,kableExtra,ggsurvfit,extrafont,DT,tibble,strex,hrbrthemes,ggstatsplot,reshape,pander,ggrepel,scales,dataMaid,gridExtra,tidytidbits,survivalAnalysis,gtsummary, Cairo,Amelia,officer,mice,naniar)
```

# Pre-processing of data

``` r
setwd("E:/PhD_projects/Realworld/Data/")
dataDF <- read.xlsx("melanoma.xlsx")

# remove patients with ipilimumab as treatment is considered outdated per analysis protocol
dataDF <- subset(dataDF, regime_correct != "Ipilimumab")

dim(dataDF)
```

    ## [1] 1199  154

``` r
# Check for duplicates
any(duplicated(dataDF$patient_id))
```

    ## [1] FALSE

``` r
# Remove duplicates
dataDF <- dataDF[!duplicated(dataDF), ]
dim(dataDF)
```

    ## [1] 1199  154

``` r
features <- c("sex", "age_1st_treat","line_correct", "regime_correct", 
              "Brain_metastases","performance_status", "bor")

features_cancer_spec <- c("recieved_adjuvant_therapy", "melanoma_diagnosis", "braf_correct", "stage_2", "LDH")

## Define which columns include survival data
OS_time <- c("OS_days")
OS_status <- c("Dead")
PFS_time <- c("PFS_days")
PFS_status <- c("Progressed")
MEL_status <- "Dead_Mel"

## Rename some column names
names(dataDF)[names(dataDF) == 'Line_correct'] <- 'line_correct'
names(dataDF)[names(dataDF) == 'bedste_respons'] <- 'bor'
names(dataDF)[names(dataDF) == 'Days_OS'] <- 'OS_days'
names(dataDF)[names(dataDF) == 'Days_PFS'] <- 'PFS_days'

categ_feats <-c("sex","regime_correct","line_correct","Brain_metastases", "performance_status", "bor",features_cancer_spec)

## DEFINE WHICH COLUMNS ARE NUMERICAL
# we also add the status and time columns
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
dataDF <- dataDF %>% dplyr::select(all_of(c(colnames(dataDF)[1], OS_time,OS_status, PFS_time,PFS_status, MEL_status, features, features_cancer_spec))) %>% distinct()

print("Numbers per treatment")
```

    ## [1] "Numbers per treatment"

``` r
table(dataDF$regime_correct)
```

    ## 
    ## Anti-PD1 ipi+nivo 
    ##      789      410

``` r
print("Numbers per response")
```

    ## [1] "Numbers per response"

``` r
table(dataDF$bor)
```

    ## 
    ##  CR  PR 
    ## 576 623

``` r
## SUMMARY OF THE DATA
summary(dataDF)
```

    ##   patient_id           OS_days            Dead           PFS_days     
    ##  Length:1199        Min.   :  42.0   Min.   :0.0000   Min.   :  42.0  
    ##  Class :character   1st Qu.: 629.5   1st Qu.:0.0000   1st Qu.: 330.5  
    ##  Mode  :character   Median :1167.0   Median :0.0000   Median : 729.0  
    ##                     Mean   :1417.0   Mean   :0.2994   Mean   : 986.7  
    ##                     3rd Qu.:2114.0   3rd Qu.:1.0000   3rd Qu.:1602.5  
    ##                     Max.   :4438.0   Max.   :1.0000   Max.   :3692.0  
    ##    Progressed        Dead_Mel          sex      age_1st_treat   line_correct
    ##  Min.   :0.0000   Min.   :0.0000   Female:487   Min.   :24.58   First:1031  
    ##  1st Qu.:0.0000   1st Qu.:0.0000   Male  :712   1st Qu.:59.13   Other: 168  
    ##  Median :0.0000   Median :0.0000                Median :70.02               
    ##  Mean   :0.4637   Mean   :0.2127                Mean   :67.45               
    ##  3rd Qu.:1.0000   3rd Qu.:0.0000                3rd Qu.:76.85               
    ##  Max.   :1.0000   Max.   :1.0000                Max.   :94.10               
    ##   regime_correct Brain_metastases performance_status bor     
    ##  Anti-PD1:789    No :1011         0   :787           CR:576  
    ##  ipi+nivo:410    Yes: 188         1   :332           PR:623  
    ##                                   2   : 74                   
    ##                                   3   :  2                   
    ##                                   4   :  1                   
    ##                                   NA's:  3                   
    ##  recieved_adjuvant_therapy                melanoma_diagnosis      braf_correct
    ##  0:1098                    Cutaneous               :949      BRAF_mutant:504  
    ##  1: 101                    Melanoma unknown primary:212      BRAF_normal:667  
    ##                            Mucosal                 : 38      NA's       : 28  
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##         stage_2        LDH     
    ##  III/M1a/M1b:503   High  :339  
    ##  M1c/M1d    :696   Normal:821  
    ##                    NA's  : 39  
    ##                                
    ##                                
    ## 

``` r
# remove some white spaces in the PS column
dataDF$performance_status <- trimws(dataDF$performance_status)
dataDF$performance_status <- as.factor(dataDF$performance_status)

dim(dataDF)
```

    ## [1] 1199   18

``` r
# 1053   17

# Rename some factors for plots
dataDF <- dataDF %>%
  dplyr::mutate(
    regime_correct = case_when(
      regime_correct == "Anti-PD1" ~ "Anti-PD1",
      regime_correct == "ipi+nivo"  ~ "Ipi+Nivo",
      TRUE                           ~ NA_character_
    ),
    regime_correct = factor(
      regime_correct,
      levels = c("Anti-PD1", "Ipi+Nivo")
    ),

    performance_status = case_when(
      performance_status %in% c("1", "2", "3","4") ~ "PS≥1",
      performance_status == "0"               ~ "PS=0",
      performance_status == "99"              ~ NA_character_,
      TRUE                    ~ NA_character_
    ),
    PS = factor(performance_status, levels = c("PS=0", "PS≥1")
                ),
    
    regime_correct = case_when(
      regime_correct == "Ipi+Nivo" ~ "Anti-PD1+Anti-CTLA4",
      regime_correct == "Anti-PD1" ~ "Anti-PD1",
      regime_correct == "99"       ~ NA_character_,
      TRUE                    ~ NA_character_
    ),
    regime_correct = factor(regime_correct, levels = c("Anti-PD1", "Anti-PD1+Anti-CTLA4")
    ),
    
    LDH = case_when(
      LDH == "High" ~ "Elevated",
      LDH == "Normal" ~ "Normal",
      LDH == "99" ~ NA_character_,
      TRUE ~ NA_character_
    ),
    LDH = factor(LDH, levels = c("Normal", "Elevated")
    ),
    
    recieved_adjuvant_therapy = case_when(
      recieved_adjuvant_therapy == 0 ~ "No",
      recieved_adjuvant_therapy == "1" ~ "Yes",
      recieved_adjuvant_therapy == "99"              ~ NA_character_,
      TRUE                    ~ NA_character_
    ),
    recieved_adjuvant_therapy = factor(recieved_adjuvant_therapy, levels = c("No", "Yes")
    ),

    melanoma_diagnosis = case_when(
      melanoma_diagnosis == "Cutaneous" ~ "Cutaneous",
      melanoma_diagnosis == "Melanoma unknown primary" ~ "Unk. primary",
      melanoma_diagnosis == "Mucosal" ~ "Mucosal",
      melanoma_diagnosis == "99"          ~ NA_character_,
      TRUE                     ~ NA_character_
    ),
    melanoma_diagnosis = factor(
      melanoma_diagnosis,
      levels = c("Cutaneous", "Unk. primary","Mucosal")
    ),
    
    `braf_correct` = case_when(
      `braf_correct` == "BRAF_normal" ~ "Wild type",
      `braf_correct` == "BRAF_mutant" ~ "Mutant",
      `braf_correct` == "99"          ~ NA_character_,
      TRUE                     ~ NA_character_
    ),
    `braf_correct` = factor(
      `braf_correct`,
      levels = c("Wild type", "Mutant")
    ),

    Brain_metastases = case_when(
     Brain_metastases == "No"  ~ "No",
     Brain_metastases == "Yes" ~ "Yes",
      TRUE                     ~ NA_character_
    ),
    Brain_metastases = factor(
      Brain_metastases,
      levels = c("No", "Yes")
    )
)

# rename some of the columns for plotting
names(dataDF)[names(dataDF) == 'recieved_adjuvant_therapy'] <- 'Previous adjuvant therapy'
names(dataDF)[names(dataDF) == 'sex'] <- 'Sex'
names(dataDF)[names(dataDF) == 'age_1st_treat'] <- 'Age'
names(dataDF)[names(dataDF) == 'regime_correct'] <- 'CPI Regimen'
names(dataDF)[names(dataDF) == 'performance_status'] <- 'ECOG Performance Status'
names(dataDF)[names(dataDF) == 'stage_2'] <- 'AJCC 8th stage'
names(dataDF)[names(dataDF) == 'Brain_metastases'] <- 'Brain metastases'
names(dataDF)[names(dataDF) == 'braf_correct'] <- 'BRAF mutation'
names(dataDF)[names(dataDF) == 'bor'] <- 'Objective response'
names(dataDF)[names(dataDF) == 'line_correct'] <- 'Treatment line'
names(dataDF)[names(dataDF) == 'melanoma_diagnosis'] <- 'Melanoma subtype'
names(dataDF)[names(dataDF) == 'patient_id'] <- 'record_id'

# redefine feature names accordingly
features <- c("Sex", "Age", "CPI Regimen", "Treatment line", "Brain metastases", "Melanoma subtype", "ECOG Performance Status","AJCC 8th stage","BRAF mutation", "Previous adjuvant therapy","LDH", "Objective response")

write.csv(dataDF, file = "Melanoma_data_polished.csv", sep = ",", append = F,row.names = F)

#Administrative censoring at 60 months
dataDF$censor_time_OS   <- pmin(dataDF$OS_days/30, 60)
dataDF$censor_status_OS <- ifelse(dataDF$OS_days/30 > 60, 0, dataDF$Dead)
#Administrative censoring at 60 months
dataDF$censor_time_PFS   <- pmin(dataDF$PFS_days/30, 60)
dataDF$censor_status_PFS <- ifelse(dataDF$PFS_days/30 > 60, 0, dataDF$Progressed)
#Administrative censoring at 60 months
dataDF$censor_status_DSS <- ifelse(dataDF$OS_days/30 > 60, 0, dataDF$Dead_Mel)
```

``` r
descr_features <- paste0("`", c(features, "CPI Regimen"), "`", collapse = " + ")
descr_formula <- as.formula(paste0("~", descr_features))

Baseline_Characteristics_table_Partial_patients <- descrTable(descr_formula , data = dataDF)
```

# OS analysis

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

    ## `height` was translated to `width`.

![](Survival_analysis_melanoma_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20OS-1.png)<!-- -->

``` r
# # To remove both (NAs and empty):
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 1131   24

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
    ##   n= 1131, number of events= 291 
    ## 
    ##                                     coef exp(coef)  se(coef)      z Pr(>|z|)
    ## SexMale                         0.212432  1.236682  0.124384  1.708  0.08766
    ## Age                             0.016529  1.016666  0.006035  2.739  0.00616
    ## `CPI Regimen`Anti-PD1           0.024973  1.025287  0.154168  0.162  0.87132
    ## `Treatment line`First           0.177051  1.193692  0.176780  1.002  0.31657
    ## `Brain metastases`No            0.136106  1.145804  0.183535  0.742  0.45834
    ## `Melanoma subtype`Unk. primary -0.080602  0.922560  0.155163 -0.519  0.60343
    ## `Melanoma subtype`Mucosal       0.304822  1.356384  0.272662  1.118  0.26359
    ## `ECOG Performance Status`PS≥1   0.201578  1.223332  0.124626  1.617  0.10578
    ## `AJCC 8th stage`M1c/M1d        -0.096231  0.908254  0.135925 -0.708  0.47896
    ## `BRAF mutation`Wild type        0.085708  1.089488  0.130380  0.657  0.51094
    ## `Previous adjuvant therapy`Yes  0.457529  1.580164  0.231876  1.973  0.04848
    ## LDHNormal                       0.062775  1.064788  0.139304  0.451  0.65225
    ## `Objective response`PR          1.965169  7.136121  0.160055 12.278  < 2e-16
    ##                                   
    ## SexMale                        .  
    ## Age                            ** 
    ## `CPI Regimen`Anti-PD1             
    ## `Treatment line`First             
    ## `Brain metastases`No              
    ## `Melanoma subtype`Unk. primary    
    ## `Melanoma subtype`Mucosal         
    ## `ECOG Performance Status`PS≥1     
    ## `AJCC 8th stage`M1c/M1d           
    ## `BRAF mutation`Wild type          
    ## `Previous adjuvant therapy`Yes *  
    ## LDHNormal                         
    ## `Objective response`PR         ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                           1.2367     0.8086    0.9691     1.578
    ## Age                               1.0167     0.9836    1.0047     1.029
    ## `CPI Regimen`Anti-PD1             1.0253     0.9753    0.7579     1.387
    ## `Treatment line`First             1.1937     0.8377    0.8441     1.688
    ## `Brain metastases`No              1.1458     0.8727    0.7996     1.642
    ## `Melanoma subtype`Unk. primary    0.9226     1.0839    0.6806     1.250
    ## `Melanoma subtype`Mucosal         1.3564     0.7373    0.7949     2.315
    ## `ECOG Performance Status`PS≥1     1.2233     0.8174    0.9582     1.562
    ## `AJCC 8th stage`M1c/M1d           0.9083     1.1010    0.6958     1.186
    ## `BRAF mutation`Wild type          1.0895     0.9179    0.8438     1.407
    ## `Previous adjuvant therapy`Yes    1.5802     0.6328    1.0031     2.489
    ## LDHNormal                         1.0648     0.9392    0.8104     1.399
    ## `Objective response`PR            7.1361     0.1401    5.2146     9.766
    ## 
    ## Concordance= 0.765  (se = 0.012 )
    ## Likelihood ratio test= 260.1  on 13 df,   p=<2e-16
    ## Wald test            = 197.1  on 13 df,   p=<2e-16
    ## Score (logrank) test = 257.1  on 13 df,   p=<2e-16

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 1131, number of events= 291 
    ## 
    ##                                     coef exp(coef)  se(coef)      z Pr(>|z|)
    ## SexMale                         0.223797  1.250817  0.124112  1.803  0.07136
    ## Age                             0.018408  1.018578  0.005981  3.078  0.00209
    ## `CPI Regimen`Anti-PD1           0.061135  1.063043  0.145829  0.419  0.67505
    ## `Melanoma subtype`Unk. primary -0.101162  0.903787  0.154482 -0.655  0.51257
    ## `Melanoma subtype`Mucosal       0.297440  1.346408  0.270650  1.099  0.27178
    ## `ECOG Performance Status`PS≥1   0.193829  1.213888  0.123770  1.566  0.11734
    ## `BRAF mutation`Wild type        0.080815  1.084170  0.128018  0.631  0.52786
    ## `Previous adjuvant therapy`Yes  0.458082  1.581039  0.230898  1.984  0.04726
    ## `Objective response`PR          1.925358  6.857604  0.158383 12.156  < 2e-16
    ##                                   
    ## SexMale                        .  
    ## Age                            ** 
    ## `CPI Regimen`Anti-PD1             
    ## `Melanoma subtype`Unk. primary    
    ## `Melanoma subtype`Mucosal         
    ## `ECOG Performance Status`PS≥1     
    ## `BRAF mutation`Wild type          
    ## `Previous adjuvant therapy`Yes *  
    ## `Objective response`PR         ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                           1.2508     0.7995    0.9807     1.595
    ## Age                               1.0186     0.9818    1.0067     1.031
    ## `CPI Regimen`Anti-PD1             1.0630     0.9407    0.7988     1.415
    ## `Melanoma subtype`Unk. primary    0.9038     1.1065    0.6677     1.223
    ## `Melanoma subtype`Mucosal         1.3464     0.7427    0.7921     2.289
    ## `ECOG Performance Status`PS≥1     1.2139     0.8238    0.9524     1.547
    ## `BRAF mutation`Wild type          1.0842     0.9224    0.8436     1.393
    ## `Previous adjuvant therapy`Yes    1.5810     0.6325    1.0055     2.486
    ## `Objective response`PR            6.8576     0.1458    5.0276     9.354
    ## 
    ## Concordance= 0.764  (se = 0.012 )
    ## Likelihood ratio test= 256.5  on 9 df,   p=<2e-16
    ## Wald test            = 193.8  on 9 df,   p=<2e-16
    ## Score (logrank) test = 253.1  on 9 df,   p=<2e-16

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_melanoma_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20OS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_melanoma_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20OS-2.png)<!-- -->

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

Male
</td>

<td style="text-align:center;border-right:1px solid;">

Female
</td>

<td style="text-align:center;">

1.31
</td>

<td style="text-align:center;">

1.04-1.66
</td>

<td style="text-align:center;border-right:1px solid;">

0.02260
</td>

<td style="text-align:center;">

1.25
</td>

<td style="text-align:center;">

0.981-1.6
</td>

<td style="text-align:center;">

0.0714
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

1.02-1.04
</td>

<td style="text-align:center;border-right:1px solid;">

0.00000
</td>

<td style="text-align:center;">

1.02
</td>

<td style="text-align:center;">

1.01-1.03
</td>

<td style="text-align:center;">

0.00209
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

1.24
</td>

<td style="text-align:center;">

0.963-1.59
</td>

<td style="text-align:center;border-right:1px solid;">

0.09620
</td>

<td style="text-align:center;">

1.06
</td>

<td style="text-align:center;">

0.799-1.41
</td>

<td style="text-align:center;">

0.675
</td>

</tr>

<tr>

<td style="text-align:left;">

Treatment line
</td>

<td style="text-align:left;">

First
</td>

<td style="text-align:center;border-right:1px solid;">

Other
</td>

<td style="text-align:center;">

1.16
</td>

<td style="text-align:center;">

0.834-1.61
</td>

<td style="text-align:center;border-right:1px solid;">

0.38200
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

Brain metastases
</td>

<td style="text-align:left;">

No
</td>

<td style="text-align:center;border-right:1px solid;">

Yes
</td>

<td style="text-align:center;">

1.08
</td>

<td style="text-align:center;">

0.79-1.48
</td>

<td style="text-align:center;border-right:1px solid;">

0.62900
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

Melanoma subtype
</td>

<td style="text-align:left;">

Unk. primary
</td>

<td style="text-align:center;border-right:1px solid;">

Cutaneous
</td>

<td style="text-align:center;">

1.07
</td>

<td style="text-align:center;">

0.799-1.43
</td>

<td style="text-align:center;border-right:1px solid;">

0.65500
</td>

<td style="text-align:center;">

0.904
</td>

<td style="text-align:center;">

0.668-1.22
</td>

<td style="text-align:center;">

0.513
</td>

</tr>

<tr>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

Mucosal
</td>

<td style="text-align:center;border-right:1px solid;">

Cutaneous
</td>

<td style="text-align:center;">

2.28
</td>

<td style="text-align:center;">

1.39-3.73
</td>

<td style="text-align:center;border-right:1px solid;">

0.00103
</td>

<td style="text-align:center;">

1.35
</td>

<td style="text-align:center;">

0.792-2.29
</td>

<td style="text-align:center;">

0.272
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

1.89
</td>

<td style="text-align:center;">

1.51-2.36
</td>

<td style="text-align:center;border-right:1px solid;">

0.00000
</td>

<td style="text-align:center;">

1.21
</td>

<td style="text-align:center;">

0.952-1.55
</td>

<td style="text-align:center;">

0.117
</td>

</tr>

<tr>

<td style="text-align:left;">

AJCC 8th stage
</td>

<td style="text-align:left;">

M1c/M1d
</td>

<td style="text-align:center;border-right:1px solid;">

III/M1a/M1b
</td>

<td style="text-align:center;">

1.13
</td>

<td style="text-align:center;">

0.897-1.41
</td>

<td style="text-align:center;border-right:1px solid;">

0.30600
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

BRAF mutation
</td>

<td style="text-align:left;">

Wild type
</td>

<td style="text-align:center;border-right:1px solid;">

Mutant
</td>

<td style="text-align:center;">

1.33
</td>

<td style="text-align:center;">

1.06-1.68
</td>

<td style="text-align:center;border-right:1px solid;">

0.01580
</td>

<td style="text-align:center;">

1.08
</td>

<td style="text-align:center;">

0.844-1.39
</td>

<td style="text-align:center;">

0.528
</td>

</tr>

<tr>

<td style="text-align:left;">

Previous adjuvant therapy
</td>

<td style="text-align:left;">

Yes
</td>

<td style="text-align:center;border-right:1px solid;">

No
</td>

<td style="text-align:center;">

1.40
</td>

<td style="text-align:center;">

0.907-2.17
</td>

<td style="text-align:center;border-right:1px solid;">

0.12800
</td>

<td style="text-align:center;">

1.58
</td>

<td style="text-align:center;">

1.01-2.49
</td>

<td style="text-align:center;">

0.0473
</td>

</tr>

<tr>

<td style="text-align:left;">

LDH
</td>

<td style="text-align:left;">

Normal
</td>

<td style="text-align:center;border-right:1px solid;">

Elevated
</td>

<td style="text-align:center;">

1.01
</td>

<td style="text-align:center;">

0.785-1.3
</td>

<td style="text-align:center;border-right:1px solid;">

0.93200
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

Objective response
</td>

<td style="text-align:left;">

PR
</td>

<td style="text-align:center;border-right:1px solid;">

CR
</td>

<td style="text-align:center;">

6.87
</td>

<td style="text-align:center;">

5.14-9.17
</td>

<td style="text-align:center;border-right:1px solid;">

0.00000
</td>

<td style="text-align:center;">

6.86
</td>

<td style="text-align:center;">

5.03-9.35
</td>

<td style="text-align:center;">

5.31e-34
</td>

</tr>

</tbody>

</table>

# PFS analysis

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

    ## `height` was translated to `width`.

![](Survival_analysis_melanoma_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20PFS-1.png)<!-- -->

``` r
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 1131   24

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
    ##   n= 1131, number of events= 501 
    ## 
    ##                                       coef exp(coef)  se(coef)      z Pr(>|z|)
    ## SexMale                           0.028112  1.028511  0.093546  0.301  0.76379
    ## Age                               0.003448  1.003454  0.004346  0.793  0.42765
    ## `CPI Regimen`Anti-PD1+Anti-CTLA4 -0.061582  0.940276  0.116009 -0.531  0.59553
    ## `Treatment line`First             0.187675  1.206441  0.134278  1.398  0.16222
    ## `Brain metastases`Yes            -0.122763  0.884473  0.133773 -0.918  0.35878
    ## `Melanoma subtype`Cutaneous       0.188078  1.206928  0.120091  1.566  0.11732
    ## `Melanoma subtype`Mucosal         0.589339  1.802797  0.241757  2.438  0.01478
    ## `ECOG Performance Status`PS≥1    -0.050933  0.950343  0.097621 -0.522  0.60185
    ## `AJCC 8th stage`M1c/M1d          -0.078724  0.924295  0.104504 -0.753  0.45127
    ## `BRAF mutation`Mutant             0.273437  1.314474  0.098347  2.780  0.00543
    ## `Previous adjuvant therapy`Yes    0.485732  1.625364  0.164846  2.947  0.00321
    ## LDHNormal                         0.111471  1.117921  0.106144  1.050  0.29363
    ## `Objective response`PR            1.930409  6.892327  0.113510 17.007  < 2e-16
    ##                                     
    ## SexMale                             
    ## Age                                 
    ## `CPI Regimen`Anti-PD1+Anti-CTLA4    
    ## `Treatment line`First               
    ## `Brain metastases`Yes               
    ## `Melanoma subtype`Cutaneous         
    ## `Melanoma subtype`Mucosal        *  
    ## `ECOG Performance Status`PS≥1       
    ## `AJCC 8th stage`M1c/M1d             
    ## `BRAF mutation`Mutant            ** 
    ## `Previous adjuvant therapy`Yes   ** 
    ## LDHNormal                           
    ## `Objective response`PR           ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                  exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                             1.0285     0.9723    0.8562     1.235
    ## Age                                 1.0035     0.9966    0.9949     1.012
    ## `CPI Regimen`Anti-PD1+Anti-CTLA4    0.9403     1.0635    0.7490     1.180
    ## `Treatment line`First               1.2064     0.8289    0.9273     1.570
    ## `Brain metastases`Yes               0.8845     1.1306    0.6805     1.150
    ## `Melanoma subtype`Cutaneous         1.2069     0.8285    0.9538     1.527
    ## `Melanoma subtype`Mucosal           1.8028     0.5547    1.1224     2.896
    ## `ECOG Performance Status`PS≥1       0.9503     1.0523    0.7848     1.151
    ## `AJCC 8th stage`M1c/M1d             0.9243     1.0819    0.7531     1.134
    ## `BRAF mutation`Mutant               1.3145     0.7608    1.0840     1.594
    ## `Previous adjuvant therapy`Yes      1.6254     0.6152    1.1766     2.245
    ## LDHNormal                           1.1179     0.8945    0.9079     1.376
    ## `Objective response`PR              6.8923     0.1451    5.5176     8.610
    ## 
    ## Concordance= 0.746  (se = 0.01 )
    ## Likelihood ratio test= 392.7  on 13 df,   p=<2e-16
    ## Wald test            = 325.5  on 13 df,   p=<2e-16
    ## Score (logrank) test = 405.3  on 13 df,   p=<2e-16

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 1131, number of events= 501 
    ## 
    ##                                     coef exp(coef)  se(coef)      z Pr(>|z|)
    ## Age                             0.002615  1.002619  0.003908  0.669  0.50338
    ## `Melanoma subtype`Cutaneous     0.217209  1.242604  0.119336  1.820  0.06874
    ## `Melanoma subtype`Mucosal       0.555000  1.741940  0.239095  2.321  0.02027
    ## `ECOG Performance Status`PS≥1  -0.044754  0.956233  0.096565 -0.463  0.64303
    ## `AJCC 8th stage`M1c/M1d        -0.185723  0.830503  0.094351 -1.968  0.04902
    ## `Previous adjuvant therapy`Yes  0.473308  1.605296  0.160215  2.954  0.00313
    ## `Objective response`PR          1.910522  6.756615  0.113038 16.902  < 2e-16
    ##                                   
    ## Age                               
    ## `Melanoma subtype`Cutaneous    .  
    ## `Melanoma subtype`Mucosal      *  
    ## `ECOG Performance Status`PS≥1     
    ## `AJCC 8th stage`M1c/M1d        *  
    ## `Previous adjuvant therapy`Yes ** 
    ## `Objective response`PR         ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                exp(coef) exp(-coef) lower .95 upper .95
    ## Age                               1.0026     0.9974    0.9950    1.0103
    ## `Melanoma subtype`Cutaneous       1.2426     0.8048    0.9835    1.5700
    ## `Melanoma subtype`Mucosal         1.7419     0.5741    1.0902    2.7832
    ## `ECOG Performance Status`PS≥1     0.9562     1.0458    0.7913    1.1555
    ## `AJCC 8th stage`M1c/M1d           0.8305     1.2041    0.6903    0.9992
    ## `Previous adjuvant therapy`Yes    1.6053     0.6229    1.1727    2.1975
    ## `Objective response`PR            6.7566     0.1480    5.4139    8.4323
    ## 
    ## Concordance= 0.737  (se = 0.01 )
    ## Likelihood ratio test= 381.4  on 7 df,   p=<2e-16
    ## Wald test            = 315.4  on 7 df,   p=<2e-16
    ## Score (logrank) test = 395.1  on 7 df,   p=<2e-16

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_melanoma_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20PFS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_melanoma_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20PFS-2.png)<!-- -->

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
| Sex | Male | Female | 1.12 | 0.937-1.33 | 0.218 |  |  | Excluded from multivariable analysis |
| Age | Fitted as continuous |  | 1.01 | 1-1.02 | 0.002 | 1 | 0.995-1.01 | 0.503 |
| CPI Regimen | Anti-PD1+Anti-CTLA4 | Anti-PD1 | 1.00 | 0.832-1.2 | 0.999 |  |  | Excluded from multivariable analysis |
| Treatment line | First | Other | 1.05 | 0.825-1.34 | 0.677 |  |  | Excluded from multivariable analysis |
| Brain metastases | Yes | No | 1.07 | 0.854-1.35 | 0.538 |  |  | Excluded from multivariable analysis |
| Melanoma subtype | Cutaneous | Unk. primary | 1.03 | 0.824-1.3 | 0.777 | 1.24 | 0.983-1.57 | 0.0687 |
| Melanoma subtype | Mucosal | Unk. primary | 2.00 | 1.28-3.14 | 0.003 | 1.74 | 1.09-2.78 | 0.0203 |
| ECOG Performance Status | PS≥1 | PS=0 | 1.34 | 1.12-1.59 | 0.001 | 0.956 | 0.791-1.16 | 0.643 |
| AJCC 8th stage | M1c/M1d | III/M1a/M1b | 1.18 | 0.987-1.4 | 0.069 | 0.831 | 0.69-0.999 | 0.049 |
| BRAF mutation | Mutant | Wild type | 1.09 | 0.917-1.3 | 0.325 |  |  | Excluded from multivariable analysis |
| Previous adjuvant therapy | Yes | No | 1.63 | 1.2-2.22 | 0.002 | 1.61 | 1.17-2.2 | 0.00313 |
| LDH | Normal | Elevated | 1.01 | 0.837-1.23 | 0.881 |  |  | Excluded from multivariable analysis |
| Objective response | PR | CR | 6.44 | 5.24-7.91 | 0.000 | 6.76 | 5.41-8.43 | 4.38e-64 |

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

Male
</td>

<td style="text-align:center;border-right:1px solid;">

Female
</td>

<td style="text-align:center;">

1.12
</td>

<td style="text-align:center;">

0.937-1.33
</td>

<td style="text-align:center;border-right:1px solid;">

0.218
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

1-1.02
</td>

<td style="text-align:center;border-right:1px solid;">

0.002
</td>

<td style="text-align:center;">

1
</td>

<td style="text-align:center;">

0.995-1.01
</td>

<td style="text-align:center;">

0.503
</td>

</tr>

<tr>

<td style="text-align:left;">

CPI Regimen
</td>

<td style="text-align:left;">

Anti-PD1+Anti-CTLA4
</td>

<td style="text-align:center;border-right:1px solid;">

Anti-PD1
</td>

<td style="text-align:center;">

1.00
</td>

<td style="text-align:center;">

0.832-1.2
</td>

<td style="text-align:center;border-right:1px solid;">

0.999
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

First
</td>

<td style="text-align:center;border-right:1px solid;">

Other
</td>

<td style="text-align:center;">

1.05
</td>

<td style="text-align:center;">

0.825-1.34
</td>

<td style="text-align:center;border-right:1px solid;">

0.677
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

Brain metastases
</td>

<td style="text-align:left;">

Yes
</td>

<td style="text-align:center;border-right:1px solid;">

No
</td>

<td style="text-align:center;">

1.07
</td>

<td style="text-align:center;">

0.854-1.35
</td>

<td style="text-align:center;border-right:1px solid;">

0.538
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

Melanoma subtype
</td>

<td style="text-align:left;">

Cutaneous
</td>

<td style="text-align:center;border-right:1px solid;">

Unk. primary
</td>

<td style="text-align:center;">

1.03
</td>

<td style="text-align:center;">

0.824-1.3
</td>

<td style="text-align:center;border-right:1px solid;">

0.777
</td>

<td style="text-align:center;">

1.24
</td>

<td style="text-align:center;">

0.983-1.57
</td>

<td style="text-align:center;">

0.0687
</td>

</tr>

<tr>

<td style="text-align:left;">

Melanoma subtype
</td>

<td style="text-align:left;">

Mucosal
</td>

<td style="text-align:center;border-right:1px solid;">

Unk. primary
</td>

<td style="text-align:center;">

2.00
</td>

<td style="text-align:center;">

1.28-3.14
</td>

<td style="text-align:center;border-right:1px solid;">

0.003
</td>

<td style="text-align:center;">

1.74
</td>

<td style="text-align:center;">

1.09-2.78
</td>

<td style="text-align:center;">

0.0203
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

1.34
</td>

<td style="text-align:center;">

1.12-1.59
</td>

<td style="text-align:center;border-right:1px solid;">

0.001
</td>

<td style="text-align:center;">

0.956
</td>

<td style="text-align:center;">

0.791-1.16
</td>

<td style="text-align:center;">

0.643
</td>

</tr>

<tr>

<td style="text-align:left;">

AJCC 8th stage
</td>

<td style="text-align:left;">

M1c/M1d
</td>

<td style="text-align:center;border-right:1px solid;">

III/M1a/M1b
</td>

<td style="text-align:center;">

1.18
</td>

<td style="text-align:center;">

0.987-1.4
</td>

<td style="text-align:center;border-right:1px solid;">

0.069
</td>

<td style="text-align:center;">

0.831
</td>

<td style="text-align:center;">

0.69-0.999
</td>

<td style="text-align:center;">

0.049
</td>

</tr>

<tr>

<td style="text-align:left;">

BRAF mutation
</td>

<td style="text-align:left;">

Mutant
</td>

<td style="text-align:center;border-right:1px solid;">

Wild type
</td>

<td style="text-align:center;">

1.09
</td>

<td style="text-align:center;">

0.917-1.3
</td>

<td style="text-align:center;border-right:1px solid;">

0.325
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

Previous adjuvant therapy
</td>

<td style="text-align:left;">

Yes
</td>

<td style="text-align:center;border-right:1px solid;">

No
</td>

<td style="text-align:center;">

1.63
</td>

<td style="text-align:center;">

1.2-2.22
</td>

<td style="text-align:center;border-right:1px solid;">

0.002
</td>

<td style="text-align:center;">

1.61
</td>

<td style="text-align:center;">

1.17-2.2
</td>

<td style="text-align:center;">

0.00313
</td>

</tr>

<tr>

<td style="text-align:left;">

LDH
</td>

<td style="text-align:left;">

Normal
</td>

<td style="text-align:center;border-right:1px solid;">

Elevated
</td>

<td style="text-align:center;">

1.01
</td>

<td style="text-align:center;">

0.837-1.23
</td>

<td style="text-align:center;border-right:1px solid;">

0.881
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

Objective response
</td>

<td style="text-align:left;">

PR
</td>

<td style="text-align:center;border-right:1px solid;">

CR
</td>

<td style="text-align:center;">

6.44
</td>

<td style="text-align:center;">

5.24-7.91
</td>

<td style="text-align:center;border-right:1px solid;">

0.000
</td>

<td style="text-align:center;">

6.76
</td>

<td style="text-align:center;">

5.41-8.43
</td>

<td style="text-align:center;">

4.38e-64
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

# 5) Apply alternating shading
dataRows <- 2:n_rows
evenRows <- dataRows[seq(1, length(dataRows), by = 2)]

# 6) Add vertical separators after col 3 and 6
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
  file = "E:/PhD_projects/Realworld/Scripts/Acquired_resistance_RW/Tables/Supplementary_Table_2.xlsx",
  overwrite = TRUE
)
```

# Melanoma specific survival

``` r
surv_time <- "censor_time_OS"
surv_status <- "censor_status_DSS"
surv_time_label <- "DSS"
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

    ## `height` was translated to `width`.

![](Survival_analysis_melanoma_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20DSS-1.png)<!-- -->

``` r
surv_time <- "censor_time_OS"
surv_status <- "censor_status_DSS"
surv_time_label <- "OS"
# # To remove both (NAs and empty):
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 1131   24

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
    ##   n= 1131, number of events= 214 
    ## 
    ##                                     coef exp(coef)  se(coef)      z Pr(>|z|)
    ## SexMale                         0.067613  1.069951  0.142018  0.476   0.6340
    ## Age                             0.005350  1.005365  0.006660  0.803   0.4218
    ## `CPI Regimen`Anti-PD1           0.057056  1.058715  0.175269  0.326   0.7448
    ## `Treatment line`First           0.112718  1.119316  0.196848  0.573   0.5669
    ## `Brain metastases`Yes          -0.149651  0.861008  0.206285 -0.725   0.4682
    ## `Melanoma subtype`Unk. primary -0.049307  0.951889  0.181150 -0.272   0.7855
    ## `Melanoma subtype`Mucosal       0.491474  1.634724  0.296576  1.657   0.0975
    ## `ECOG Performance Status`PS≥1   0.127085  1.135513  0.147256  0.863   0.3881
    ## `AJCC 8th stage`M1c/M1d        -0.008496  0.991540  0.160236 -0.053   0.9577
    ## `BRAF mutation`Wild type       -0.018962  0.981217  0.151195 -0.125   0.9002
    ## `Previous adjuvant therapy`Yes  0.559134  1.749157  0.252932  2.211   0.0271
    ## LDHNormal                       0.188333  1.207236  0.163472  1.152   0.2493
    ## `Objective response`PR          2.376063 10.762452  0.213753 11.116   <2e-16
    ##                                   
    ## SexMale                           
    ## Age                               
    ## `CPI Regimen`Anti-PD1             
    ## `Treatment line`First             
    ## `Brain metastases`Yes             
    ## `Melanoma subtype`Unk. primary    
    ## `Melanoma subtype`Mucosal      .  
    ## `ECOG Performance Status`PS≥1     
    ## `AJCC 8th stage`M1c/M1d           
    ## `BRAF mutation`Wild type          
    ## `Previous adjuvant therapy`Yes *  
    ## LDHNormal                         
    ## `Objective response`PR         ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                           1.0700    0.93462    0.8100     1.413
    ## Age                               1.0054    0.99466    0.9923     1.019
    ## `CPI Regimen`Anti-PD1             1.0587    0.94454    0.7509     1.493
    ## `Treatment line`First             1.1193    0.89340    0.7610     1.646
    ## `Brain metastases`Yes             0.8610    1.16143    0.5747     1.290
    ## `Melanoma subtype`Unk. primary    0.9519    1.05054    0.6674     1.358
    ## `Melanoma subtype`Mucosal         1.6347    0.61172    0.9141     2.923
    ## `ECOG Performance Status`PS≥1     1.1355    0.88066    0.8508     1.515
    ## `AJCC 8th stage`M1c/M1d           0.9915    1.00853    0.7243     1.357
    ## `BRAF mutation`Wild type          0.9812    1.01914    0.7296     1.320
    ## `Previous adjuvant therapy`Yes    1.7492    0.57170    1.0654     2.872
    ## LDHNormal                         1.2072    0.82834    0.8763     1.663
    ## `Objective response`PR           10.7625    0.09292    7.0789    16.363
    ## 
    ## Concordance= 0.772  (se = 0.013 )
    ## Likelihood ratio test= 223.7  on 13 df,   p=<2e-16
    ## Wald test            = 144.7  on 13 df,   p=<2e-16
    ## Score (logrank) test = 216.9  on 13 df,   p=<2e-16

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 1131, number of events= 214 
    ## 
    ##                                     coef exp(coef)  se(coef)      z Pr(>|z|)
    ## Age                             0.007759  1.007789  0.006099  1.272   0.2033
    ## `Melanoma subtype`Unk. primary -0.068259  0.934018  0.179484 -0.380   0.7037
    ## `Melanoma subtype`Mucosal       0.507761  1.661567  0.290555  1.748   0.0805
    ## `ECOG Performance Status`PS≥1   0.117351  1.124514  0.145647  0.806   0.4204
    ## `AJCC 8th stage`M1c/M1d        -0.118831  0.887958  0.144277 -0.824   0.4101
    ## `Previous adjuvant therapy`Yes  0.560253  1.751115  0.246771  2.270   0.0232
    ## `Objective response`PR          2.363094 10.623776  0.213594 11.063   <2e-16
    ##                                   
    ## Age                               
    ## `Melanoma subtype`Unk. primary    
    ## `Melanoma subtype`Mucosal      .  
    ## `ECOG Performance Status`PS≥1     
    ## `AJCC 8th stage`M1c/M1d           
    ## `Previous adjuvant therapy`Yes *  
    ## `Objective response`PR         ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                exp(coef) exp(-coef) lower .95 upper .95
    ## Age                                1.008    0.99227    0.9958     1.020
    ## `Melanoma subtype`Unk. primary     0.934    1.07064    0.6570     1.328
    ## `Melanoma subtype`Mucosal          1.662    0.60184    0.9401     2.937
    ## `ECOG Performance Status`PS≥1      1.125    0.88927    0.8453     1.496
    ## `AJCC 8th stage`M1c/M1d            0.888    1.12618    0.6692     1.178
    ## `Previous adjuvant therapy`Yes     1.751    0.57106    1.0796     2.840
    ## `Objective response`PR            10.624    0.09413    6.9898    16.147
    ## 
    ## Concordance= 0.77  (se = 0.013 )
    ## Likelihood ratio test= 220.8  on 7 df,   p=<2e-16
    ## Wald test            = 141.8  on 7 df,   p=<2e-16
    ## Score (logrank) test = 213.8  on 7 df,   p=<2e-16

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_melanoma_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20DSS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_melanoma_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20DSS-2.png)<!-- -->

# Combined Univariable vs Multivariable Forest Plots - DSS

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

n_cols <- ncol(results_df)

# align: left for first two, center for the rest
aligns <- c("l", "l", rep("c", n_cols - 2))

results_df %>%
  kable(
    format = "html",
    caption = "Univariable vs Multivariable Cox Results - DSS",
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

Univariable vs Multivariable Cox Results - DSS
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

Male
</td>

<td style="text-align:center;border-right:1px solid;">

Female
</td>

<td style="text-align:center;">

1.13
</td>

<td style="text-align:center;">

0.864-1.47
</td>

<td style="text-align:center;border-right:1px solid;">

0.376000
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

1.02
</td>

<td style="text-align:center;">

1.01-1.03
</td>

<td style="text-align:center;border-right:1px solid;">

0.004240
</td>

<td style="text-align:center;">

1.01
</td>

<td style="text-align:center;">

0.996-1.02
</td>

<td style="text-align:center;">

0.203
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

1.11
</td>

<td style="text-align:center;">

0.833-1.47
</td>

<td style="text-align:center;border-right:1px solid;">

0.479000
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

First
</td>

<td style="text-align:center;border-right:1px solid;">

Other
</td>

<td style="text-align:center;">

1.04
</td>

<td style="text-align:center;">

0.719-1.5
</td>

<td style="text-align:center;border-right:1px solid;">

0.836000
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

Brain metastases
</td>

<td style="text-align:left;">

Yes
</td>

<td style="text-align:center;border-right:1px solid;">

No
</td>

<td style="text-align:center;">

1.06
</td>

<td style="text-align:center;">

0.749-1.5
</td>

<td style="text-align:center;border-right:1px solid;">

0.738000
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

Melanoma subtype
</td>

<td style="text-align:left;">

Unk. primary
</td>

<td style="text-align:center;border-right:1px solid;">

Cutaneous
</td>

<td style="text-align:center;">

1.08
</td>

<td style="text-align:center;">

0.767-1.52
</td>

<td style="text-align:center;border-right:1px solid;">

0.666000
</td>

<td style="text-align:center;">

0.934
</td>

<td style="text-align:center;">

0.657-1.33
</td>

<td style="text-align:center;">

0.704
</td>

</tr>

<tr>

<td style="text-align:left;">

Melanoma subtype
</td>

<td style="text-align:left;">

Mucosal
</td>

<td style="text-align:center;border-right:1px solid;">

Cutaneous
</td>

<td style="text-align:center;">

2.77
</td>

<td style="text-align:center;">

1.63-4.7
</td>

<td style="text-align:center;border-right:1px solid;">

0.000155
</td>

<td style="text-align:center;">

1.66
</td>

<td style="text-align:center;">

0.94-2.94
</td>

<td style="text-align:center;">

0.0805
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

1.68
</td>

<td style="text-align:center;">

1.29-2.18
</td>

<td style="text-align:center;border-right:1px solid;">

0.000120
</td>

<td style="text-align:center;">

1.12
</td>

<td style="text-align:center;">

0.845-1.5
</td>

<td style="text-align:center;">

0.42
</td>

</tr>

<tr>

<td style="text-align:left;">

AJCC 8th stage
</td>

<td style="text-align:left;">

M1c/M1d
</td>

<td style="text-align:center;border-right:1px solid;">

III/M1a/M1b
</td>

<td style="text-align:center;">

1.32
</td>

<td style="text-align:center;">

1.01-1.73
</td>

<td style="text-align:center;border-right:1px solid;">

0.043700
</td>

<td style="text-align:center;">

0.888
</td>

<td style="text-align:center;">

0.669-1.18
</td>

<td style="text-align:center;">

0.41
</td>

</tr>

<tr>

<td style="text-align:left;">

BRAF mutation
</td>

<td style="text-align:left;">

Wild type
</td>

<td style="text-align:center;border-right:1px solid;">

Mutant
</td>

<td style="text-align:center;">

1.13
</td>

<td style="text-align:center;">

0.865-1.47
</td>

<td style="text-align:center;border-right:1px solid;">

0.372000
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

Previous adjuvant therapy
</td>

<td style="text-align:left;">

Yes
</td>

<td style="text-align:center;border-right:1px solid;">

No
</td>

<td style="text-align:center;">

1.65
</td>

<td style="text-align:center;">

1.03-2.65
</td>

<td style="text-align:center;border-right:1px solid;">

0.037400
</td>

<td style="text-align:center;">

1.75
</td>

<td style="text-align:center;">

1.08-2.84
</td>

<td style="text-align:center;">

0.0232
</td>

</tr>

<tr>

<td style="text-align:left;">

LDH
</td>

<td style="text-align:left;">

Normal
</td>

<td style="text-align:center;border-right:1px solid;">

Elevated
</td>

<td style="text-align:center;">

1.03
</td>

<td style="text-align:center;">

0.766-1.39
</td>

<td style="text-align:center;border-right:1px solid;">

0.841000
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

Objective response
</td>

<td style="text-align:left;">

PR
</td>

<td style="text-align:center;border-right:1px solid;">

CR
</td>

<td style="text-align:center;">

10.50
</td>

<td style="text-align:center;">

7.12-15.6
</td>

<td style="text-align:center;border-right:1px solid;">

0.000000
</td>

<td style="text-align:center;">

10.6
</td>

<td style="text-align:center;">

6.99-16.1
</td>

<td style="text-align:center;">

1.89e-28
</td>

</tr>

</tbody>

</table>

## Numeric table with 5-year OS, PFS, DSS

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
  fit_DSS <- survfit(Surv(censor_time_OS,  censor_status_DSS) ~ 1, data = d)

  os_ci  <- get5yr_ci(fit_OS)
  pfs_ci <- get5yr_ci(fit_PFS)
  DSS_ci <- get5yr_ci(fit_DSS)

  data.frame(
    Response       = resp,
    OS_5yr         = os_ci$Estimate,
    OS_5yr_Lower   = os_ci$Lower95,
    OS_5yr_Upper   = os_ci$Upper95,
    PFS_5yr        = pfs_ci$Estimate,
    PFS_5yr_Lower  = pfs_ci$Lower95,
    PFS_5yr_Upper  = pfs_ci$Upper95,
    DSS_5yr        = DSS_ci$Estimate,
    DSS_5yr_Lower  = DSS_ci$Lower95,
    DSS_5yr_Upper  = DSS_ci$Upper95
  )
})

five_year_table <- bind_rows(results) %>%
  mutate(across(matches("_5yr"), ~ . * 100, .names = "{.col}_pct")) %>%
  select(Response,
         OS_5yr_pct, OS_5yr_Lower_pct, OS_5yr_Upper_pct,
         PFS_5yr_pct, PFS_5yr_Lower_pct, PFS_5yr_Upper_pct,
         DSS_5yr_pct, DSS_5yr_Lower_pct, DSS_5yr_Upper_pct)

print(five_year_table)
```

    ##   Response OS_5yr_pct OS_5yr_Lower_pct OS_5yr_Upper_pct PFS_5yr_pct
    ## 1       CR   87.48618         84.45029         90.63121    73.74845
    ## 2       PR   44.01466         39.30324         49.29087    18.77321
    ##   PFS_5yr_Lower_pct PFS_5yr_Upper_pct DSS_5yr_pct DSS_5yr_Lower_pct
    ## 1          69.79731          77.92325    93.44162          91.13900
    ## 2          15.07256          23.38245    52.69265          47.70295
    ##   DSS_5yr_Upper_pct
    ## 1          95.80241
    ## 2          58.20427

``` r
# pull out labels from five_year_table
x_max=60
rownames(five_year_table) <- five_year_table$Response
os_lbls <- five_year_table %>%
  transmute(
    Response,
    os_label = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                       Response, OS_5yr_pct, OS_5yr_Lower_pct, OS_5yr_Upper_pct))

pfs_lbls <- five_year_table %>%
  transmute(
    Response,
    pfs_label = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                        Response, PFS_5yr_pct, PFS_5yr_Lower_pct, PFS_5yr_Upper_pct))

DSS_lbls <- five_year_table %>%
  transmute(
    Response,
    DSS_label = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                        Response, DSS_5yr_pct, DSS_5yr_Lower_pct, DSS_5yr_Upper_pct))

#adding a syntactically correct name for OR, as it was failing in the Survfit step
dataDF$Objective_response = dataDF$`Objective response`

ypos <- c(CR = five_year_table$OS_5yr_pct[1], PR = five_year_table$OS_5yr_pct[2])

# OS plot
fit_OS <- survfit(
  Surv(censor_time_OS, censor_status_OS) ~ Objective_response,
  data = dataDF
)
names(fit_OS$strata) <- gsub("Objective_response=", "", names(fit_OS$strata))

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
      y        = "OS (%)",
      title = "Melanoma - Overall Survival"
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
   annotate(
    "text",
    x     = x_max*0.25,
    y     = (ypos[names(ypos)]*0.9)/100,
    label = paste0("5 year OS for ",os_lbls$os_label),
    size  = 10,
    hjust = 0
    )
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the ggsurvfit package.
    ##   Please report the issue at <https://github.com/pharmaverse/ggsurvfit/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

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
      title = "Melanoma - Progression-free Survival"
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
   annotate(
    "text",
    x     = x_max*0.25,
    y     = (ypos[names(ypos)]*0.9)/100,
    label = paste0("5 year PFS for ",pfs_lbls$pfs_label),
    size  = 10,
    hjust = 0
    )
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

``` r
# DSS plot
ypos <- c(CR = five_year_table$DSS_5yr_pct[1], PR = five_year_table$DSS_5yr_pct[2])

fit_DSS <- survfit(
  Surv(censor_time_OS, censor_status_DSS) ~ Objective_response,
  data = dataDF
)
names(fit_DSS$strata) <- gsub("Objective_response=", "", names(fit_DSS$strata))

p_DSS <- ggsurvfit(fit_DSS, size = 1.5) +
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
      y        = "DSS (%)",
      title = "Melanoma - Disease-specific Survival"
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
   annotate(
    "text",
    x     = x_max*0.25,
    y     = ypos[names(ypos)]*0.9/100,
    label = paste0("5 year DSS for ",DSS_lbls$DSS_label),
    size  = 10,
    hjust = 0
    )
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

``` r
p_pfs
```

![](Survival_analysis_melanoma_files/figure-gfm/-%205%20year%20KM%20curves-1.png)<!-- -->

``` r
p_os
```

![](Survival_analysis_melanoma_files/figure-gfm/-%205%20year%20KM%20curves-2.png)<!-- -->

``` r
p_DSS
```

![](Survival_analysis_melanoma_files/figure-gfm/-%205%20year%20KM%20curves-3.png)<!-- -->

## 18 Month landmark analysis

``` r
#remove all patients with an event before 18 months
dataDF_18mo <- subset(dataDF, OS_days > 540 & PFS_days > 540)

#reset time
dataDF_18mo$censor_time_OS <- dataDF_18mo$censor_time_OS-18
dataDF_18mo$censor_time_PFS <- dataDF_18mo$censor_time_PFS-18

surv_time <- "censor_time_PFS"
surv_status <- "censor_status_PFS"
surv_time_label <- "PFS"
# PFS plot
results <- lapply(responses, function(resp) {
  d <- filter(dataDF_18mo, `Objective response` == resp)
  fit_OS  <- survfit(Surv(censor_time_OS,  censor_status_OS)  ~ 1, data = d)
  fit_PFS <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1, data = d)
  fit_DSS <- survfit(Surv(censor_time_OS,  censor_status_DSS) ~ 1, data = d)

  os_ci  <- get5yr_ci(fit_OS)
  pfs_ci <- get5yr_ci(fit_PFS)
  DSS_ci <- get5yr_ci(fit_DSS)

  data.frame(
    Response       = resp,
    OS_5yr         = os_ci$Estimate,
    OS_5yr_Lower   = os_ci$Lower95,
    OS_5yr_Upper   = os_ci$Upper95,
    PFS_5yr        = pfs_ci$Estimate,
    PFS_5yr_Lower  = pfs_ci$Lower95,
    PFS_5yr_Upper  = pfs_ci$Upper95,
    DSS_5yr        = DSS_ci$Estimate,
    DSS_5yr_Lower  = DSS_ci$Lower95,
    DSS_5yr_Upper  = DSS_ci$Upper95
  )
})

five_year_table <- bind_rows(results) %>%
  mutate(across(matches("_5yr"), ~ . * 100, .names = "{.col}_pct")) %>%
  select(Response,
         OS_5yr_pct, OS_5yr_Lower_pct, OS_5yr_Upper_pct,
         PFS_5yr_pct, PFS_5yr_Lower_pct, PFS_5yr_Upper_pct,
         DSS_5yr_pct, DSS_5yr_Lower_pct, DSS_5yr_Upper_pct)


pfs_lbls <- five_year_table %>%
  transmute(
    Response,
    pfs_label = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                        Response, PFS_5yr_pct, PFS_5yr_Lower_pct, PFS_5yr_Upper_pct))

ypos <- c(CR = five_year_table$PFS_5yr_pct[1], PR = five_year_table$PFS_5yr_pct[2])

fit_PFS <- survfit(
  Surv(censor_time_PFS, censor_status_PFS) ~ Objective_response,
  data = dataDF_18mo
)
names(fit_PFS$strata) <- gsub("Objective_response=", "", names(fit_PFS$strata))

x_max <- 42
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
      title = "Melanoma - 18-month landmark analysis (PFS)"
    ) +
    scale_x_continuous(breaks = seq(0, x_max, by = 6), limits = c(0, x_max)) +
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) 

os_lbls <- five_year_table %>%
  transmute(
    Response,
    os_label = sprintf("%s: %.0f%% (%.0f–%.0f%%)",
                        Response, OS_5yr_pct, OS_5yr_Lower_pct, OS_5yr_Upper_pct))

ypos <- c(CR = five_year_table$OS_5yr_pct[1], PR = five_year_table$OS_5yr_pct[2])

fit_OS <- survfit(
  Surv(censor_time_OS, censor_status_OS) ~ Objective_response,
  data = dataDF_18mo
)
names(fit_OS$strata) <- gsub("Objective_response=", "", names(fit_OS$strata))

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
      y        = "OS (%)",
      title = "Melanoma - 18-month landmark analysis (OS)"
    ) +
    scale_x_continuous(breaks = seq(0, x_max, by = 6), limits = c(0, x_max)) +
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#00BA38", "PR"="#00BFC4"))


hr_OS_resp_18  <- coxph(Surv(censor_time_OS, censor_status_OS) ~ Objective_response, data = dataDF_18mo)
s_cox_18   <- summary(hr_OS_resp_18)
# extract HR, CI and p
hr_18     <- s_cox_18$coefficients[, "exp(coef)"]
ci_lo_18  <- s_cox_18$conf.int[, "lower .95"]
ci_hi_18  <- s_cox_18$conf.int[, "upper .95"]
pval_18   <- s_cox_18$coefficients[, "Pr(>|z|)"]
hr_lab_18 <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)", hr_18, ci_lo_18, ci_hi_18)
p_lab_18  <- ifelse(pval_18 < 0.001, "p < 0.001", 
                    paste0("p = ", sub("^0\\.", ".", sprintf("%.4f", pval_18))))

print(paste0("18-month landmark OS ",
             hr_lab_18, " - ", p_lab_18))
```

    ## [1] "18-month landmark OS HR=3.97 (95% CI 2.65–5.96) - p < 0.001"

``` r
p_os+
   annotate(
    "text",
    x     = 0.25,
    y     = 0.25,
    label = paste0(hr_lab_18, " \n ", p_lab_18),
    size  = 10,
    hjust = 0
    )
```

![](Survival_analysis_melanoma_files/figure-gfm/landmark%20analysis-1.png)<!-- -->

``` r
hr_PFS_resp_18  <- coxph(Surv(censor_time_PFS, censor_status_PFS) ~ Objective_response, data = dataDF_18mo)
s_cox_18   <- summary(hr_PFS_resp_18)
# extract HR, CI and p
hr_18     <- s_cox_18$coefficients[, "exp(coef)"]
ci_lo_18  <- s_cox_18$conf.int[, "lower .95"]
ci_hi_18  <- s_cox_18$conf.int[, "upper .95"]
pval_18   <- s_cox_18$coefficients[, "Pr(>|z|)"]
hr_lab_18 <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)", hr_18, ci_lo_18, ci_hi_18)
p_lab_18  <- ifelse(pval_18 < 0.001, "p < 0.001", 
                    paste0("p = ", sub("^0\\.", ".", sprintf("%.4f", pval_18))))

print(paste0("18-month landmark PFS ",
             hr_lab_18, " - ", p_lab_18))
```

    ## [1] "18-month landmark PFS HR=4.38 (95% CI 3.28–5.86) - p < 0.001"

``` r
p_pfs+
   annotate(
    "text",
    x     = 0.25,
    y     = 0.25,
    label = paste0(hr_lab_18, " \n ", p_lab_18),
    size  = 10,
    hjust = 0
    )
```

![](Survival_analysis_melanoma_files/figure-gfm/landmark%20analysis-2.png)<!-- -->

\`\`\`
