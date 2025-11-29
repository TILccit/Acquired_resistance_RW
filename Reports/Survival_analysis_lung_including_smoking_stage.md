Survival analysis Lung with smoking and stage
================
Mario Presti
First created on Feb 2025. Updated on 29 November 2025

- [Introduction](#introduction)
- [Loading Libraries](#loading-libraries)
  - [Pre-processing of data](#pre-processing-of-data)
- [OS analysis](#os-analysis)
- [Combined Univariable vs Multivariable Forest Plots -
  OS](#combined-univariable-vs-multivariable-forest-plots---os)
- [PFS analysis](#pfs-analysis)
- [Combined Univariable vs Multivariable Forest Plots -
  PFS](#combined-univariable-vs-multivariable-forest-plots---pfs)
- [Five-Year Survival Estimates by
  Response](#five-year-survival-estimates-by-response)

# Introduction

Multivariate Survival Analysis on Lung database

# Loading Libraries

``` r
# Load required libraries
library(pacman)
pacman::p_load(data.table,stringr,tidyverse,ggplot2,MatchIt,survival,survminer,survey,compareGroups, forestmodel, forestplot, openxlsx,kableExtra,ggsurvfit,extrafont,DT,tibble,strex,hrbrthemes,ggstatsplot,reshape,pander,ggrepel,scales,dataMaid,gridExtra,tidytidbits,survivalAnalysis,gtsummary, Cairo,Amelia,officer,mice,naniar)
```

## Pre-processing of data

``` r
setwd("E:/PhD_projects/Realworld/Data/")
dataDF <- read.xlsx("lung.xlsx")

# Check for duplicates
any(duplicated(dataDF$patient_id))
```

    ## [1] FALSE

``` r
# Remove duplicates
dataDF <- dataDF[!duplicated(dataDF), ]
dim(dataDF)
```

    ## [1] 667  56

``` r
#filter out only first line patients
dataDF <- subset(dataDF, Line_correct != "First")

names(dataDF)[names(dataDF) == 'bedste_respons'] <- 'bor'
#correct PDL1 column
names(dataDF)[names(dataDF) == 'PD-L1_correct'] <- 'PDL1_correct'
dataDF$PDL1_correct <- gsub("&lt;", "<", dataDF$PDL1_correct)
features <- c("sex", "age_1st_treat","Line_correct", 
              "Brain_metastases","PS", "bor")

features_cancer_spec <- c("Histology", "PDL1_correct","Smoking_", "RA_Stage")

## Define which columns include survival data
OS_time <- c("Days_OS")
OS_status <- c("Dead")
PFS_time <- c("PFS_days")
PFS_status <- c("Progressed")

## DEFINE WHICH COLUMNS ARE CATEGORICAL
categ_feats <-c("sex","Line_correct","Brain_metastases", "PS", "bor",features_cancer_spec)

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
#remove white spaces at the end of features (typos)
dataDF <- as.data.frame(apply(dataDF, MARGIN = 2, FUN = trimws))

# Transform character columns to categorical, to have the different levels
dataDF <- dataDF %>% mutate(across(all_of(categ_feats), as.factor))
dataDF <- dataDF %>% dplyr::select(all_of(c(colnames(dataDF)[1], OS_time,OS_status, PFS_time,PFS_status, features, features_cancer_spec))) %>% distinct()

dim(dataDF)
```

    ## [1] 190  15

``` r
# 264   16

print("Numbers per response")
```

    ## [1] "Numbers per response"

``` r
table(dataDF$bor)
```

    ## 
    ##  CR  PR 
    ##  21 169

``` r
dataDF$PS <- as.factor(dataDF$PS)


#identify some categories in the histology column
squam_hist <- unique(grep("\\bsquamous\\b", dataDF$Histology, value = T, ignore.case = T))
nonsquam_hist <- setdiff(unique(dataDF$Histology), squam_hist)
stage_I_III <- unique(grep("1|I-III|3|2B|2b|2|<4", dataDF$RA_Stage, value = T, ignore.case = T))
dataDF$RA_Stage[which(dataDF$RA_Stage=="Unspecified")] <- NA
# Rename some factors for plots
dataDF <- dataDF %>%
  dplyr::mutate(
    sex = case_when(
      sex %in% c("female", "Female") ~ "Female",
      sex == "Male" ~ "Male",
      TRUE                    ~ NA_character_
    ),
    sex = factor(sex, levels = c("Male", "Female")
    ),
    Brain_metastases = case_when(
     Brain_metastases == "No"  ~ "No",
     Brain_metastases == "Yes" ~ "Yes",
      TRUE                     ~ NA_character_
    ),
    Brain_metastases = factor(
      Brain_metastases,
      levels = c("No", "Yes")
    ),
    PS = case_when(
      PS == "PS=0" ~ "PS=0",
      PS == "PSover0" ~ "PS≥1",
      TRUE                     ~ NA_character_
    ),
    PS = factor(PS, levels = c("PS=0", "PS≥1")
    ),
    
    Histology_correct = case_when(
      Histology %in% squam_hist  ~ "Squamous",
      Histology %in% nonsquam_hist ~ "Nonsquamous",
      TRUE                     ~ NA_character_
    ),
    Histology_correct = factor(Histology_correct, levels = c("Nonsquamous", "Squamous")
    ),
    
    PDL1_correct = case_when(
      PDL1_correct == "over49" ~ "PDL1≥50%",
      PDL1_correct == "under50" ~ "PDL1<50%",
      TRUE                     ~ NA_character_
    ),
    PDL1_correct = factor(PDL1_correct, levels = c("PDL1≥50%", "PDL1<50%")
    ),
    RA_Stage = case_when(
      RA_Stage %in% stage_I_III ~ "< IV",
      RA_Stage %in% c("4", "IV") ~ "IV",
      TRUE                     ~ NA_character_
    ),
    RA_Stage = factor(RA_Stage, levels = c("< IV", "IV")
    ),
    Smoking_ = case_when(
      Smoking_ %in% c("Current smoker", "Active", "active") ~ "Current",
      Smoking_ %in% c("Previous smoker", "Former", "former") ~ "Former",
      Smoking_ %in% c("Never") ~ "Never",
      TRUE                     ~ NA_character_
    ),
    Smoking_ = factor(Smoking_, levels = c("Current", "Former", "Never"))
)

dataDF[,"CPI Regimen"] <- "Anti-PD1"

## SUMMARY OF THE DATA
summary(dataDF)
```

    ##  RA_Patient-ID        Days_OS              Dead             PFS_days        
    ##  Length:190         Length:190         Length:190         Length:190        
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##   Progressed            sex      age_1st_treat      Line_correct
    ##  Length:190         Male  : 81   Length:190         Other:190   
    ##  Class :character   Female:109   Class :character               
    ##  Mode  :character                Mode  :character               
    ##                                                                 
    ##                                                                 
    ##                                                                 
    ##                                                                 
    ##  Brain_metastases    PS      bor                        Histology  
    ##  No :159          PS=0: 42   CR: 21   Adeno                  :100  
    ##  Yes: 31          PS≥1:147   PR:169   Squamous               : 27  
    ##                   NA's:  1            adeno                  : 22  
    ##                                       squamous               : 12  
    ##                                       Adenocarcinoma         : 10  
    ##                                       Squamous cell carcinoma:  7  
    ##                                       (Other)                : 12  
    ##    PDL1_correct    Smoking_   RA_Stage     Histology_correct CPI Regimen       
    ##  PDL1≥50%:93    Current: 58   < IV: 26   Nonsquamous:144     Length:190        
    ##  PDL1<50%:70    Former :126   IV  :161   Squamous   : 46     Class :character  
    ##  NA's    :27    Never  :  5   NA's:  3                       Mode  :character  
    ##                 NA's   :  1                                                    
    ##                                                                                
    ##                                                                                
    ## 

``` r
names(dataDF)[names(dataDF) == 'Days_OS'] <- 'OS_days'
names(dataDF)[names(dataDF) == 'Days_PFS'] <- 'PFS_days'

### NOW MAKE SURE THAT THE COLUMNS ARE IN THE CORRECT FORMAT NEEDED FOR ANALYSIS.
## Numeric columns
dataDF <- dataDF %>% dplyr::mutate(across(all_of(numeric_feats), as.numeric))

print("Summary after releveling some variables: ")
```

    ## [1] "Summary after releveling some variables: "

``` r
summary(dataDF)
```

    ##  RA_Patient-ID         OS_days            Dead           PFS_days     
    ##  Length:190         Min.   :  82.0   Min.   :0.0000   Min.   :  49.0  
    ##  Class :character   1st Qu.: 638.2   1st Qu.:0.0000   1st Qu.: 335.5  
    ##  Mode  :character   Median :1300.5   Median :1.0000   Median : 731.0  
    ##                     Mean   :1444.5   Mean   :0.6211   Mean   :1061.9  
    ##                     3rd Qu.:2235.2   3rd Qu.:1.0000   3rd Qu.:1733.0  
    ##                     Max.   :3187.0   Max.   :1.0000   Max.   :3037.0  
    ##                                                                       
    ##    Progressed         sex      age_1st_treat   Line_correct Brain_metastases
    ##  Min.   :0.0000   Male  : 81   Min.   :40.00   Other:190    No :159         
    ##  1st Qu.:1.0000   Female:109   1st Qu.:62.00                Yes: 31         
    ##  Median :1.0000                Median :67.27                                
    ##  Mean   :0.7684                Mean   :66.57                                
    ##  3rd Qu.:1.0000                3rd Qu.:72.00                                
    ##  Max.   :1.0000                Max.   :80.00                                
    ##                                                                             
    ##     PS      bor                        Histology     PDL1_correct    Smoking_  
    ##  PS=0: 42   CR: 21   Adeno                  :100   PDL1≥50%:93    Current: 58  
    ##  PS≥1:147   PR:169   Squamous               : 27   PDL1<50%:70    Former :126  
    ##  NA's:  1            adeno                  : 22   NA's    :27    Never  :  5  
    ##                      squamous               : 12                  NA's   :  1  
    ##                      Adenocarcinoma         : 10                               
    ##                      Squamous cell carcinoma:  7                               
    ##                      (Other)                : 12                               
    ##  RA_Stage     Histology_correct CPI Regimen       
    ##  < IV: 26   Nonsquamous:144     Length:190        
    ##  IV  :161   Squamous   : 46     Class :character  
    ##  NA's:  3                       Mode  :character  
    ##                                                   
    ##                                                   
    ##                                                   
    ## 

``` r
# rename some of the columns for plotting
names(dataDF)[names(dataDF) == 'sex'] <- 'Sex'
names(dataDF)[names(dataDF) == 'age_1st_treat'] <- 'Age'
names(dataDF)[names(dataDF) == 'PS'] <- 'ECOG Performance Status'
names(dataDF)[names(dataDF) == 'Brain_metastases'] <- 'Brain metastases'
names(dataDF)[names(dataDF) == 'bor'] <- 'Objective response'
names(dataDF)[names(dataDF) == 'Line_correct'] <- 'Treatment line'
names(dataDF)[names(dataDF) == 'RA_Patient-ID'] <- 'record_id'
names(dataDF)[names(dataDF) == 'Histology_correct'] <- 'Histology subtype'
names(dataDF)[names(dataDF) == 'PDL1_correct'] <- 'PD-L1 status'
names(dataDF)[names(dataDF) == 'RA_Stage'] <- 'Stage'
names(dataDF)[names(dataDF) == 'Smoking_'] <- 'Smoking status'

# redefine feature names accordingly
features <- c("Sex", "Age", "ECOG Performance Status", "Brain metastases", "Histology subtype", "PD-L1 status","Stage", "Smoking status","Objective response")

#write.csv(dataDF, file = "Lung_data_polished.csv", sep = ",", append = F,row.names = F)

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

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20OS-1.png)<!-- -->

``` r
#save this for the comparison at the end with the 18 month landmark

hr_OS_resp <-multiple_uni[[length(multiple_uni)]]$coxph
```

``` r
# # To remove both (NAs and empty):
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 159  21

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
    ##   n= 159, number of events= 83 
    ## 
    ##                                   coef exp(coef) se(coef)      z Pr(>|z|)    
    ## SexMale                        0.29572   1.34409  0.23851  1.240 0.215032    
    ## Age                            0.05435   1.05585  0.01683  3.229 0.001244 ** 
    ## `ECOG Performance Status`PS≥1  0.05525   1.05680  0.29168  0.189 0.849762    
    ## `Brain metastases`Yes          1.02057   2.77478  0.29027  3.516 0.000438 ***
    ## `Histology subtype`Squamous    0.12216   1.12993  0.31135  0.392 0.694809    
    ## `PD-L1 status`PDL1<50%         0.09833   1.10332  0.24012  0.409 0.682175    
    ## Stage< IV                      0.21263   1.23692  0.31473  0.676 0.499307    
    ## `Smoking status`Current        0.10284   1.10832  0.24523  0.419 0.674953    
    ## `Smoking status`Never         -0.04215   0.95873  0.62188 -0.068 0.945964    
    ## `Objective response`PR         1.60440   4.97490  0.59830  2.682 0.007327 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                          1.3441     0.7440    0.8422     2.145
    ## Age                              1.0558     0.9471    1.0216     1.091
    ## `ECOG Performance Status`PS≥1    1.0568     0.9462    0.5967     1.872
    ## `Brain metastases`Yes            2.7748     0.3604    1.5709     4.901
    ## `Histology subtype`Squamous      1.1299     0.8850    0.6138     2.080
    ## `PD-L1 status`PDL1<50%           1.1033     0.9064    0.6892     1.766
    ## Stage< IV                        1.2369     0.8085    0.6675     2.292
    ## `Smoking status`Current          1.1083     0.9023    0.6854     1.792
    ## `Smoking status`Never            0.9587     1.0430    0.2834     3.244
    ## `Objective response`PR           4.9749     0.2010    1.5400    16.071
    ## 
    ## Concordance= 0.686  (se = 0.029 )
    ## Likelihood ratio test= 37  on 10 df,   p=6e-05
    ## Wald test            = 30.75  on 10 df,   p=6e-04
    ## Score (logrank) test = 33.03  on 10 df,   p=3e-04

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 159, number of events= 83 
    ## 
    ##                                   coef exp(coef) se(coef)      z Pr(>|z|)    
    ## SexMale                        0.30418   1.35551  0.23723  1.282 0.199771    
    ## Age                            0.05613   1.05773  0.01658  3.385 0.000713 ***
    ## `ECOG Performance Status`PS≥1  0.03334   1.03390  0.28875  0.115 0.908073    
    ## `Brain metastases`Yes          0.97785   2.65872  0.28339  3.451 0.000559 ***
    ## `Histology subtype`Squamous    0.15826   1.17147  0.29751  0.532 0.594766    
    ## `Smoking status`Current        0.10869   1.11482  0.24264  0.448 0.654178    
    ## `Smoking status`Never         -0.13297   0.87549  0.60746 -0.219 0.826731    
    ## `Objective response`PR         1.58596   4.88400  0.59646  2.659 0.007838 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                          1.3555     0.7377    0.8515     2.158
    ## Age                              1.0577     0.9454    1.0239     1.093
    ## `ECOG Performance Status`PS≥1    1.0339     0.9672    0.5871     1.821
    ## `Brain metastases`Yes            2.6587     0.3761    1.5256     4.633
    ## `Histology subtype`Squamous      1.1715     0.8536    0.6539     2.099
    ## `Smoking status`Current          1.1148     0.8970    0.6929     1.794
    ## `Smoking status`Never            0.8755     1.1422    0.2662     2.880
    ## `Objective response`PR           4.8840     0.2048    1.5173    15.721
    ## 
    ## Concordance= 0.687  (se = 0.03 )
    ## Likelihood ratio test= 36.4  on 8 df,   p=1e-05
    ## Wald test            = 30.58  on 8 df,   p=2e-04
    ## Score (logrank) test = 32.74  on 8 df,   p=7e-05

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20OS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20OS-2.png)<!-- -->

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

1.64
</td>

<td style="text-align:center;">

1.12-2.39
</td>

<td style="text-align:center;border-right:1px solid;">

0.011100
</td>

<td style="text-align:center;">

1.36
</td>

<td style="text-align:center;">

0.851-2.16
</td>

<td style="text-align:center;">

0.2
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

1.05
</td>

<td style="text-align:center;">

1.03-1.08
</td>

<td style="text-align:center;border-right:1px solid;">

0.000184
</td>

<td style="text-align:center;">

1.06
</td>

<td style="text-align:center;">

1.02-1.09
</td>

<td style="text-align:center;">

0.000713
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

1.42
</td>

<td style="text-align:center;">

0.862-2.33
</td>

<td style="text-align:center;border-right:1px solid;">

0.170000
</td>

<td style="text-align:center;">

1.03
</td>

<td style="text-align:center;">

0.587-1.82
</td>

<td style="text-align:center;">

0.908
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

1.50
</td>

<td style="text-align:center;">

0.928-2.41
</td>

<td style="text-align:center;border-right:1px solid;">

0.098500
</td>

<td style="text-align:center;">

2.66
</td>

<td style="text-align:center;">

1.53-4.63
</td>

<td style="text-align:center;">

0.000559
</td>

</tr>

<tr>

<td style="text-align:left;">

Histology subtype
</td>

<td style="text-align:left;">

Squamous
</td>

<td style="text-align:center;border-right:1px solid;">

Nonsquamous
</td>

<td style="text-align:center;">

1.49
</td>

<td style="text-align:center;">

0.981-2.26
</td>

<td style="text-align:center;border-right:1px solid;">

0.061800
</td>

<td style="text-align:center;">

1.17
</td>

<td style="text-align:center;">

0.654-2.1
</td>

<td style="text-align:center;">

0.595
</td>

</tr>

<tr>

<td style="text-align:left;">

PD-L1 status
</td>

<td style="text-align:left;">

PDL1\<50%
</td>

<td style="text-align:center;border-right:1px solid;">

PDL1≥50%
</td>

<td style="text-align:center;">

1.25
</td>

<td style="text-align:center;">

0.81-1.92
</td>

<td style="text-align:center;border-right:1px solid;">

0.317000
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

Stage
</td>

<td style="text-align:left;">

\< IV
</td>

<td style="text-align:center;border-right:1px solid;">

IV
</td>

<td style="text-align:center;">

1.17
</td>

<td style="text-align:center;">

0.687-1.99
</td>

<td style="text-align:center;border-right:1px solid;">

0.565000
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

Smoking status
</td>

<td style="text-align:left;">

Current
</td>

<td style="text-align:center;border-right:1px solid;">

Former
</td>

<td style="text-align:center;">

1.39
</td>

<td style="text-align:center;">

0.932-2.09
</td>

<td style="text-align:center;border-right:1px solid;">

0.106000
</td>

<td style="text-align:center;">

1.11
</td>

<td style="text-align:center;">

0.693-1.79
</td>

<td style="text-align:center;">

0.654
</td>

</tr>

<tr>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

</td>

<td style="text-align:center;border-right:1px solid;">

Former
</td>

<td style="text-align:center;">

1.01
</td>

<td style="text-align:center;">

0.318-3.22
</td>

<td style="text-align:center;border-right:1px solid;">

0.984000
</td>

<td style="text-align:center;">

0.875
</td>

<td style="text-align:center;">

0.266-2.88
</td>

<td style="text-align:center;">

0.827
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

4.37
</td>

<td style="text-align:center;">

1.61-11.9
</td>

<td style="text-align:center;border-right:1px solid;">

0.003820
</td>

<td style="text-align:center;">

4.88
</td>

<td style="text-align:center;">

1.52-15.7
</td>

<td style="text-align:center;">

0.00784
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

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20PFS-1.png)<!-- -->

``` r
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 159  21

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
    ##   n= 159, number of events= 108 
    ## 
    ##                                   coef exp(coef) se(coef)      z Pr(>|z|)    
    ## SexMale                        0.10163   1.10697  0.21618  0.470 0.638272    
    ## Age                            0.03656   1.03724  0.01452  2.519 0.011775 *  
    ## `ECOG Performance Status`PS≥1  0.06343   1.06548  0.26488  0.239 0.810749    
    ## `Brain metastases`Yes          0.99097   2.69386  0.26032  3.807 0.000141 ***
    ## `Histology subtype`Squamous    0.27814   1.32067  0.26762  1.039 0.298661    
    ## `PD-L1 status`PDL1<50%         0.27659   1.31862  0.21336  1.296 0.194865    
    ## Stage< IV                      0.11301   1.11964  0.29102  0.388 0.697777    
    ## `Smoking status`Former        -0.16211   0.85035  0.53747 -0.302 0.762945    
    ## `Smoking status`Current       -0.13724   0.87176  0.55907 -0.245 0.806088    
    ## `Objective response`PR         1.80144   6.05837  0.52383  3.439 0.000584 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                          1.1070     0.9034    0.7246     1.691
    ## Age                              1.0372     0.9641    1.0081     1.067
    ## `ECOG Performance Status`PS≥1    1.0655     0.9385    0.6340     1.791
    ## `Brain metastases`Yes            2.6939     0.3712    1.6173     4.487
    ## `Histology subtype`Squamous      1.3207     0.7572    0.7816     2.231
    ## `PD-L1 status`PDL1<50%           1.3186     0.7584    0.8680     2.003
    ## Stage< IV                        1.1196     0.8931    0.6329     1.981
    ## `Smoking status`Former           0.8503     1.1760    0.2966     2.438
    ## `Smoking status`Current          0.8718     1.1471    0.2914     2.608
    ## `Objective response`PR           6.0584     0.1651    2.1701    16.914
    ## 
    ## Concordance= 0.665  (se = 0.027 )
    ## Likelihood ratio test= 44.72  on 10 df,   p=2e-06
    ## Wald test            = 34.87  on 10 df,   p=1e-04
    ## Score (logrank) test = 39.22  on 10 df,   p=2e-05

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 159, number of events= 108 
    ## 
    ##                                  coef exp(coef) se(coef)     z Pr(>|z|)    
    ## SexMale                       0.11765   1.12485  0.21197 0.555 0.578865    
    ## Age                           0.03726   1.03796  0.01441 2.585 0.009741 ** 
    ## `ECOG Performance Status`PS≥1 0.05394   1.05542  0.25537 0.211 0.832709    
    ## `Brain metastases`Yes         0.97994   2.66430  0.25674 3.817 0.000135 ***
    ## `Histology subtype`Squamous   0.27925   1.32214  0.26328 1.061 0.288843    
    ## `PD-L1 status`PDL1<50%        0.27303   1.31395  0.20878 1.308 0.190960    
    ## `Objective response`PR        1.77422   5.89568  0.51622 3.437 0.000588 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                           1.125     0.8890    0.7425     1.704
    ## Age                               1.038     0.9634    1.0090     1.068
    ## `ECOG Performance Status`PS≥1     1.055     0.9475    0.6398     1.741
    ## `Brain metastases`Yes             2.664     0.3753    1.6108     4.407
    ## `Histology subtype`Squamous       1.322     0.7563    0.7892     2.215
    ## `PD-L1 status`PDL1<50%            1.314     0.7611    0.8727     1.978
    ## `Objective response`PR            5.896     0.1696    2.1435    16.216
    ## 
    ## Concordance= 0.664  (se = 0.027 )
    ## Likelihood ratio test= 44.5  on 7 df,   p=2e-07
    ## Wald test            = 35.02  on 7 df,   p=1e-05
    ## Score (logrank) test = 38.95  on 7 df,   p=2e-06

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20PFS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20PFS-2.png)<!-- -->

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
| Sex | Male | Female | 1.44 | 1.03-2.02 | 0.035 | 1.12 | 0.742-1.7 | 0.579 |
| Age | Fitted as continuous |  | 1.04 | 1.01-1.06 | 0.002 | 1.04 | 1.01-1.07 | 0.00974 |
| ECOG Performance Status | PS≥1 | PS=0 | 1.35 | 0.873-2.09 | 0.177 | 1.06 | 0.64-1.74 | 0.833 |
| Brain metastases | Yes | No | 1.53 | 0.995-2.35 | 0.053 | 2.66 | 1.61-4.41 | 0.000135 |
| Histology subtype | Squamous | Nonsquamous | 1.72 | 1.19-2.5 | 0.004 | 1.32 | 0.789-2.22 | 0.289 |
| PD-L1 status | PDL1\<50% | PDL1≥50% | 1.45 | 0.997-2.11 | 0.052 | 1.31 | 0.873-1.98 | 0.191 |
| Stage | \< IV | IV | 1.03 | 0.627-1.69 | 0.908 |  |  | Excluded from multivariable analysis |
| Smoking status | Current | Never | 1.00 | 0.368-2.73 | 0.996 |  |  | Excluded from multivariable analysis |
| Smoking status | Former | Never | 1.23 | 0.441-3.42 | 0.694 |  |  | Excluded from multivariable analysis |
| Objective response | PR | CR | 5.28 | 2.16-12.9 | 0.000 | 5.9 | 2.14-16.2 | 0.000588 |

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

1.44
</td>

<td style="text-align:center;">

1.03-2.02
</td>

<td style="text-align:center;border-right:1px solid;">

0.035
</td>

<td style="text-align:center;">

1.12
</td>

<td style="text-align:center;">

0.742-1.7
</td>

<td style="text-align:center;">

0.579
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

1.04
</td>

<td style="text-align:center;">

1.01-1.06
</td>

<td style="text-align:center;border-right:1px solid;">

0.002
</td>

<td style="text-align:center;">

1.04
</td>

<td style="text-align:center;">

1.01-1.07
</td>

<td style="text-align:center;">

0.00974
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

1.35
</td>

<td style="text-align:center;">

0.873-2.09
</td>

<td style="text-align:center;border-right:1px solid;">

0.177
</td>

<td style="text-align:center;">

1.06
</td>

<td style="text-align:center;">

0.64-1.74
</td>

<td style="text-align:center;">

0.833
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

1.53
</td>

<td style="text-align:center;">

0.995-2.35
</td>

<td style="text-align:center;border-right:1px solid;">

0.053
</td>

<td style="text-align:center;">

2.66
</td>

<td style="text-align:center;">

1.61-4.41
</td>

<td style="text-align:center;">

0.000135
</td>

</tr>

<tr>

<td style="text-align:left;">

Histology subtype
</td>

<td style="text-align:left;">

Squamous
</td>

<td style="text-align:center;border-right:1px solid;">

Nonsquamous
</td>

<td style="text-align:center;">

1.72
</td>

<td style="text-align:center;">

1.19-2.5
</td>

<td style="text-align:center;border-right:1px solid;">

0.004
</td>

<td style="text-align:center;">

1.32
</td>

<td style="text-align:center;">

0.789-2.22
</td>

<td style="text-align:center;">

0.289
</td>

</tr>

<tr>

<td style="text-align:left;">

PD-L1 status
</td>

<td style="text-align:left;">

PDL1\<50%
</td>

<td style="text-align:center;border-right:1px solid;">

PDL1≥50%
</td>

<td style="text-align:center;">

1.45
</td>

<td style="text-align:center;">

0.997-2.11
</td>

<td style="text-align:center;border-right:1px solid;">

0.052
</td>

<td style="text-align:center;">

1.31
</td>

<td style="text-align:center;">

0.873-1.98
</td>

<td style="text-align:center;">

0.191
</td>

</tr>

<tr>

<td style="text-align:left;">

Stage
</td>

<td style="text-align:left;">

\< IV
</td>

<td style="text-align:center;border-right:1px solid;">

IV
</td>

<td style="text-align:center;">

1.03
</td>

<td style="text-align:center;">

0.627-1.69
</td>

<td style="text-align:center;border-right:1px solid;">

0.908
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

Smoking status
</td>

<td style="text-align:left;">

Current
</td>

<td style="text-align:center;border-right:1px solid;">

Never
</td>

<td style="text-align:center;">

1.00
</td>

<td style="text-align:center;">

0.368-2.73
</td>

<td style="text-align:center;border-right:1px solid;">

0.996
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

Smoking status
</td>

<td style="text-align:left;">

Former
</td>

<td style="text-align:center;border-right:1px solid;">

Never
</td>

<td style="text-align:center;">

1.23
</td>

<td style="text-align:center;">

0.441-3.42
</td>

<td style="text-align:center;border-right:1px solid;">

0.694
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

5.28
</td>

<td style="text-align:center;">

2.16-12.9
</td>

<td style="text-align:center;border-right:1px solid;">

0.000
</td>

<td style="text-align:center;">

5.9
</td>

<td style="text-align:center;">

2.14-16.2
</td>

<td style="text-align:center;">

0.000588
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
  file = "E:/PhD_projects/Realworld/Scripts/Acquired_resistance_RW/Tables/Supplementary_Table_4.xlsx",
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
    ## 1       CR   80.95238         65.78531         99.61629    76.19048
    ## 2       PR   38.27328         31.53960         46.44461    21.80056
    ##   PFS_5yr_Lower_pct PFS_5yr_Upper_pct
    ## 1          59.98805          96.76909
    ## 2          16.31199          29.13590

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
    add_risktable_strata_symbol(symbol = "•", size = 20, family = "Arial Unicode MS")+
    labs(
      x        = "Months after treatment initiation",
      y        = "PFS (%)",
      title = "NSCLC - Overall Survival"
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#F8766D", "PR"="#B79F00")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#F8766D", "PR"="#B79F00")) +
   annotate(
    "text",
    x     = x_max*0.25,
    y     = ypos[names(ypos)]*0.9/100,
    label = paste0("5 year OS for ",os_lbls$os_lbls),
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
      title = "NSCLC - Progression-free Survival"
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
    scale_color_manual(name=c("CR", "PR"),values = c("CR"="#F8766D", "PR"="#B79F00")) +
    scale_fill_manual(name=c("CR", "PR"),values = c("CR"="#F8766D", "PR"="#B79F00")) +
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

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/-%205%20year%20KM%20curves-1.png)<!-- -->

``` r
p_os
```

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

![](Survival_analysis_lung_including_smoking_stage_files/figure-gfm/-%205%20year%20KM%20curves-2.png)<!-- -->
