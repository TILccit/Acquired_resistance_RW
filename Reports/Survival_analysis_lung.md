Survival analysis Lung
================
Mario Presti
First created on Feb 2025. Updated on 06 August 2025

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
# 264  49

names(dataDF)[names(dataDF) == 'bedste_respons'] <- 'bor'
#correct PDL1 column
names(dataDF)[names(dataDF) == 'PD-L1_correct'] <- 'PDL1_correct'
dataDF$PDL1_correct <- gsub("&lt;", "<", dataDF$PDL1_correct)
features <- c("sex", "age_1st_treat","Line_correct", 
              "Brain_metastases","PS", "bor")

features_cancer_spec <- c("Histology", "PDL1_correct")

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

    ## [1] 667  13

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
    ##  73 594

``` r
dataDF$PS <- as.factor(dataDF$PS)


#identify some categories in the histology column
squam_hist <- unique(grep("\\bsquamous\\b", dataDF$Histology, value = T, ignore.case = T))
nonsquam_hist <- setdiff(unique(dataDF$Histology), squam_hist)

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
    PDL1_correct = factor(PDL1_correct, levels = c("PDL1≥50%", "PDL1<50%"))
)

dataDF[,"CPI Regimen"] <- "Anti-PD1"

## SUMMARY OF THE DATA
summary(dataDF)
```

    ##  RA_Patient-ID        Days_OS              Dead             PFS_days        
    ##  Length:667         Length:667         Length:667         Length:667        
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##   Progressed            sex      age_1st_treat      Line_correct
    ##  Length:667         Male  :281   Length:667         First:477   
    ##  Class :character   Female:386   Class :character   Other:190   
    ##  Mode  :character                Mode  :character               
    ##                                                                 
    ##                                                                 
    ##                                                                 
    ##                                                                 
    ##  Brain_metastases    PS      bor                        Histology  
    ##  No :593          PS=0:186   CR: 73   Adeno                  :261  
    ##  Yes: 74          PS≥1:478   PR:594   adeno                  :129  
    ##                   NA's:  3            Adenocarcinoma         :101  
    ##                                       Squamous               : 60  
    ##                                       squamous               : 32  
    ##                                       Squamous cell carcinoma: 29  
    ##                                       (Other)                : 55  
    ##    PDL1_correct   Histology_correct CPI Regimen       
    ##  PDL1≥50%:546   Nonsquamous:546     Length:667        
    ##  PDL1<50%: 87   Squamous   :121     Class :character  
    ##  NA's    : 34                       Mode  :character  
    ##                                                       
    ##                                                       
    ##                                                       
    ## 

``` r
names(dataDF)[names(dataDF) == 'Days_OS'] <- 'OS_days'
names(dataDF)[names(dataDF) == 'Days_PFS'] <- 'PFS_days'

### NOW MAKE SURE THAT THE COLUMNS ARE IN THE CORRECT FORMAT NEEDED FOR ANALYSIS.
## Numeric columns
dataDF <- dataDF %>% dplyr::mutate(across(all_of(numeric_feats), as.numeric))

## Categorical
dataDF <- dataDF %>%  dplyr::mutate(across(all_of(categ_feats), as.character))

print("Summary after releveling some variables: ")
```

    ## [1] "Summary after releveling some variables: "

``` r
summary(dataDF)
```

    ##  RA_Patient-ID         OS_days          Dead           PFS_days     
    ##  Length:667         Min.   :  63   Min.   :0.0000   Min.   :  49.0  
    ##  Class :character   1st Qu.: 670   1st Qu.:0.0000   1st Qu.: 324.0  
    ##  Mode  :character   Median :1176   Median :1.0000   Median : 661.0  
    ##                     Mean   :1335   Mean   :0.5832   Mean   : 943.1  
    ##                     3rd Qu.:1952   3rd Qu.:1.0000   3rd Qu.:1341.5  
    ##                     Max.   :4163   Max.   :1.0000   Max.   :4163.0  
    ##    Progressed         sex            age_1st_treat   Line_correct      
    ##  Min.   :0.0000   Length:667         Min.   :40.00   Length:667        
    ##  1st Qu.:1.0000   Class :character   1st Qu.:63.00   Class :character  
    ##  Median :1.0000   Mode  :character   Median :69.00   Mode  :character  
    ##  Mean   :0.7826                      Mean   :68.45                     
    ##  3rd Qu.:1.0000                      3rd Qu.:75.00                     
    ##  Max.   :1.0000                      Max.   :92.00                     
    ##  Brain_metastases        PS                bor             Histology        
    ##  Length:667         Length:667         Length:667         Length:667        
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  PDL1_correct         Histology_correct CPI Regimen       
    ##  Length:667         Nonsquamous:546     Length:667        
    ##  Class :character   Squamous   :121     Class :character  
    ##  Mode  :character                       Mode  :character  
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

# redefine feature names accordingly
features <- c("Sex", "Age", "ECOG Performance Status", "Treatment line", "Brain metastases", "Histology subtype", "PD-L1 status","Objective response")

write.csv(dataDF, file = "Lung_data_polished.csv", sep = ",", append = F,row.names = F)

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

![](Survival_analysis_lung_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20OS-1.png)<!-- -->

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

    ## [1] 631  19

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
    ##   n= 631, number of events= 319 
    ## 
    ##                                    coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## SexMale                        0.319730  1.376756  0.114853  2.784 0.005372 ** 
    ## Age                            0.033223  1.033781  0.007351  4.519 6.20e-06 ***
    ## `ECOG Performance Status`PS≥1  0.008374  1.008409  0.127097  0.066 0.947470    
    ## `Treatment line`Other         -0.001681  0.998320  0.153060 -0.011 0.991236    
    ## `Brain metastases`Yes          0.622959  1.864437  0.167623  3.716 0.000202 ***
    ## `Histology subtype`Squamous   -0.048106  0.953033  0.158073 -0.304 0.760878    
    ## `PD-L1 status`PDL1≥50%        -0.120286  0.886667  0.190356 -0.632 0.527454    
    ## `Objective response`PR         1.274949  3.578518  0.265421  4.803 1.56e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                          1.3768     0.7263    1.0992     1.724
    ## Age                              1.0338     0.9673    1.0190     1.049
    ## `ECOG Performance Status`PS≥1    1.0084     0.9917    0.7861     1.294
    ## `Treatment line`Other            0.9983     1.0017    0.7396     1.348
    ## `Brain metastases`Yes            1.8644     0.5364    1.3424     2.590
    ## `Histology subtype`Squamous      0.9530     1.0493    0.6991     1.299
    ## `PD-L1 status`PDL1≥50%           0.8867     1.1278    0.6106     1.288
    ## `Objective response`PR           3.5785     0.2794    2.1270     6.020
    ## 
    ## Concordance= 0.636  (se = 0.016 )
    ## Likelihood ratio test= 75.13  on 8 df,   p=5e-13
    ## Wald test            = 62.77  on 8 df,   p=1e-10
    ## Score (logrank) test = 66.45  on 8 df,   p=2e-11

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 631, number of events= 319 
    ## 
    ##                                    coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## SexMale                        0.321615  1.379354  0.114742  2.803 0.005064 ** 
    ## Age                            0.032742  1.033284  0.007203  4.545 5.48e-06 ***
    ## `ECOG Performance Status`PS≥1  0.006900  1.006924  0.126285  0.055 0.956424    
    ## `Brain metastases`Yes          0.618586  1.856302  0.167284  3.698 0.000217 ***
    ## `Histology subtype`Squamous   -0.023884  0.976399  0.153788 -0.155 0.876582    
    ## `Objective response`PR         1.278672  3.591868  0.265309  4.820 1.44e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                          1.3794     0.7250    1.1016     1.727
    ## Age                              1.0333     0.9678    1.0188     1.048
    ## `ECOG Performance Status`PS≥1    1.0069     0.9931    0.7861     1.290
    ## `Brain metastases`Yes            1.8563     0.5387    1.3374     2.577
    ## `Histology subtype`Squamous      0.9764     1.0242    0.7223     1.320
    ## `Objective response`PR           3.5919     0.2784    2.1354     6.042
    ## 
    ## Concordance= 0.636  (se = 0.016 )
    ## Likelihood ratio test= 74.61  on 6 df,   p=5e-14
    ## Wald test            = 62.45  on 6 df,   p=1e-11
    ## Score (logrank) test = 66.09  on 6 df,   p=3e-12

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_lung_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20OS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_lung_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20OS-2.png)<!-- -->

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

1.450
</td>

<td style="text-align:center;">

1.17-1.79
</td>

<td style="text-align:center;border-right:1px solid;">

5.40e-04
</td>

<td style="text-align:center;">

1.38
</td>

<td style="text-align:center;">

1.1-1.73
</td>

<td style="text-align:center;">

0.00506
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

1.030
</td>

<td style="text-align:center;">

1.02-1.05
</td>

<td style="text-align:center;border-right:1px solid;">

3.20e-06
</td>

<td style="text-align:center;">

1.03
</td>

<td style="text-align:center;">

1.02-1.05
</td>

<td style="text-align:center;">

5.48e-06
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

1.170
</td>

<td style="text-align:center;">

0.926-1.49
</td>

<td style="text-align:center;border-right:1px solid;">

1.85e-01
</td>

<td style="text-align:center;">

1.01
</td>

<td style="text-align:center;">

0.786-1.29
</td>

<td style="text-align:center;">

0.956
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

1.080
</td>

<td style="text-align:center;">

0.856-1.35
</td>

<td style="text-align:center;border-right:1px solid;">

5.32e-01
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

1.390
</td>

<td style="text-align:center;">

1.02-1.89
</td>

<td style="text-align:center;border-right:1px solid;">

3.74e-02
</td>

<td style="text-align:center;">

1.86
</td>

<td style="text-align:center;">

1.34-2.58
</td>

<td style="text-align:center;">

0.000217
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

1.260
</td>

<td style="text-align:center;">

0.974-1.64
</td>

<td style="text-align:center;border-right:1px solid;">

7.77e-02
</td>

<td style="text-align:center;">

0.976
</td>

<td style="text-align:center;">

0.722-1.32
</td>

<td style="text-align:center;">

0.877
</td>

</tr>

<tr>

<td style="text-align:left;">

PD-L1 status
</td>

<td style="text-align:left;">

PDL1≥50%
</td>

<td style="text-align:center;border-right:1px solid;">

PDL1\<50%
</td>

<td style="text-align:center;">

0.884
</td>

<td style="text-align:center;">

0.649-1.21
</td>

<td style="text-align:center;border-right:1px solid;">

4.36e-01
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

3.390
</td>

<td style="text-align:center;">

2.08-5.52
</td>

<td style="text-align:center;border-right:1px solid;">

1.00e-06
</td>

<td style="text-align:center;">

3.59
</td>

<td style="text-align:center;">

2.14-6.04
</td>

<td style="text-align:center;">

1.44e-06
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

![](Survival_analysis_lung_files/figure-gfm/Multiple%20Univariate%20analysis%20-%20Including%20missing%20data%20-%20PFS-1.png)<!-- -->

``` r
dataDF_complete <- dataDF %>%
  drop_na() %>%  # Remove NA values
  filter_all(all_vars(. != "")) %>%  # Remove empty strings
  filter_all(all_vars(trimws(.) != ""))  # Remove strings with only spaces

dim(dataDF_complete)
```

    ## [1] 631  19

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
    ##   n= 631, number of events= 450 
    ## 
    ##                                    coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## SexMale                        0.172866  1.188707  0.097758  1.768  0.07701 .  
    ## Age                            0.017655  1.017812  0.005966  2.959  0.00308 ** 
    ## `ECOG Performance Status`PS≥1 -0.086098  0.917504  0.106797 -0.806  0.42013    
    ## `Treatment line`Other         -0.302345  0.739083  0.134381 -2.250  0.02445 *  
    ## `Brain metastases`Yes          0.526161  1.692422  0.145188  3.624  0.00029 ***
    ## `Histology subtype`Squamous   -0.020813  0.979402  0.133937 -0.155  0.87651    
    ## `PD-L1 status`PDL1≥50%        -0.351837  0.703395  0.163228 -2.155  0.03112 *  
    ## `Objective response`PR         1.544903  4.687516  0.225251  6.859 6.96e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                               exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                          1.1887     0.8413    0.9814    1.4397
    ## Age                              1.0178     0.9825    1.0060    1.0298
    ## `ECOG Performance Status`PS≥1    0.9175     1.0899    0.7442    1.1311
    ## `Treatment line`Other            0.7391     1.3530    0.5679    0.9618
    ## `Brain metastases`Yes            1.6924     0.5909    1.2733    2.2495
    ## `Histology subtype`Squamous      0.9794     1.0210    0.7533    1.2734
    ## `PD-L1 status`PDL1≥50%           0.7034     1.4217    0.5108    0.9686
    ## `Objective response`PR           4.6875     0.2133    3.0145    7.2892
    ## 
    ## Concordance= 0.613  (se = 0.014 )
    ## Likelihood ratio test= 106.8  on 8 df,   p=<2e-16
    ## Wald test            = 77.26  on 8 df,   p=2e-13
    ## Score (logrank) test = 87.25  on 8 df,   p=2e-15

``` r
summary(cox_model_sign)
```

    ## Call:
    ## coxph(formula = fmla_sign, data = dataDF_complete)
    ## 
    ##   n= 631, number of events= 450 
    ## 
    ##                                  coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## SexMale                      0.164343  1.178618  0.097284  1.689 0.091160 .  
    ## Age                          0.019448  1.019638  0.005913  3.289 0.001005 ** 
    ## `Brain metastases`Yes        0.496353  1.642720  0.144735  3.429 0.000605 ***
    ## `Histology subtype`Squamous  0.003778  1.003785  0.132839  0.028 0.977310    
    ## `PD-L1 status`PDL1≥50%      -0.146702  0.863552  0.135422 -1.083 0.278680    
    ## `Objective response`PR       1.540276  4.665878  0.225008  6.845 7.63e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                             exp(coef) exp(-coef) lower .95 upper .95
    ## SexMale                        1.1786     0.8485    0.9740     1.426
    ## Age                            1.0196     0.9807    1.0079     1.032
    ## `Brain metastases`Yes          1.6427     0.6087    1.2370     2.182
    ## `Histology subtype`Squamous    1.0038     0.9962    0.7737     1.302
    ## `PD-L1 status`PDL1≥50%         0.8636     1.1580    0.6622     1.126
    ## `Objective response`PR         4.6659     0.2143    3.0020     7.252
    ## 
    ## Concordance= 0.612  (se = 0.014 )
    ## Likelihood ratio test= 100.5  on 6 df,   p=<2e-16
    ## Wald test            = 71.48  on 6 df,   p=2e-13
    ## Score (logrank) test = 81.41  on 6 df,   p=2e-15

``` r
forest_model_all <- forest_model(cox_model_all,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)
forest_model_sign <- forest_model(cox_model_sign,  # Set x-axis limits
                             exponentiate = TRUE)  # Display hazard ratios on log scale)

# Display the plot
print(forest_model_all)
```

![](Survival_analysis_lung_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20PFS-1.png)<!-- -->

``` r
print(forest_model_sign)
```

![](Survival_analysis_lung_files/figure-gfm/Multivariate%20Survival%20Analysis%20-%20only%20complete%20cases%20-%20PFS-2.png)<!-- -->

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
| Sex | Male | Female | 1.240 | 1.04-1.49 | 0.017 | 1.18 | 0.974-1.43 | 0.0912 |
| Age | Fitted as continuous |  | 1.020 | 1.01-1.03 | 0.000 | 1.02 | 1.01-1.03 | 0.00101 |
| ECOG Performance Status | PS≥1 | PS=0 | 1.020 | 0.834-1.24 | 0.858 |  |  | Excluded from multivariable analysis |
| Treatment line | Other | First | 0.900 | 0.738-1.1 | 0.296 |  |  | Excluded from multivariable analysis |
| Brain metastases | Yes | No | 1.300 | 0.992-1.7 | 0.057 | 1.64 | 1.24-2.18 | 0.000605 |
| Histology subtype | Squamous | Nonsquamous | 1.220 | 0.976-1.53 | 0.080 | 1 | 0.774-1.3 | 0.977 |
| PD-L1 status | PDL1≥50% | PDL1\<50% | 0.836 | 0.645-1.08 | 0.176 | 0.864 | 0.662-1.13 | 0.279 |
| Objective response | PR | CR | 4.450 | 2.92-6.77 | 0.000 | 4.67 | 3-7.25 | 7.63e-12 |

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

1.240
</td>

<td style="text-align:center;">

1.04-1.49
</td>

<td style="text-align:center;border-right:1px solid;">

0.017
</td>

<td style="text-align:center;">

1.18
</td>

<td style="text-align:center;">

0.974-1.43
</td>

<td style="text-align:center;">

0.0912
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

1.020
</td>

<td style="text-align:center;">

1.01-1.03
</td>

<td style="text-align:center;border-right:1px solid;">

0.000
</td>

<td style="text-align:center;">

1.02
</td>

<td style="text-align:center;">

1.01-1.03
</td>

<td style="text-align:center;">

0.00101
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

1.020
</td>

<td style="text-align:center;">

0.834-1.24
</td>

<td style="text-align:center;border-right:1px solid;">

0.858
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

0.900
</td>

<td style="text-align:center;">

0.738-1.1
</td>

<td style="text-align:center;border-right:1px solid;">

0.296
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

1.300
</td>

<td style="text-align:center;">

0.992-1.7
</td>

<td style="text-align:center;border-right:1px solid;">

0.057
</td>

<td style="text-align:center;">

1.64
</td>

<td style="text-align:center;">

1.24-2.18
</td>

<td style="text-align:center;">

0.000605
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

1.220
</td>

<td style="text-align:center;">

0.976-1.53
</td>

<td style="text-align:center;border-right:1px solid;">

0.080
</td>

<td style="text-align:center;">

1
</td>

<td style="text-align:center;">

0.774-1.3
</td>

<td style="text-align:center;">

0.977
</td>

</tr>

<tr>

<td style="text-align:left;">

PD-L1 status
</td>

<td style="text-align:left;">

PDL1≥50%
</td>

<td style="text-align:center;border-right:1px solid;">

PDL1\<50%
</td>

<td style="text-align:center;">

0.836
</td>

<td style="text-align:center;">

0.645-1.08
</td>

<td style="text-align:center;border-right:1px solid;">

0.176
</td>

<td style="text-align:center;">

0.864
</td>

<td style="text-align:center;">

0.662-1.13
</td>

<td style="text-align:center;">

0.279
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

4.450
</td>

<td style="text-align:center;">

2.92-6.77
</td>

<td style="text-align:center;border-right:1px solid;">

0.000
</td>

<td style="text-align:center;">

4.67
</td>

<td style="text-align:center;">

3-7.25
</td>

<td style="text-align:center;">

7.63e-12
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
    ## 1       CR   74.90839         65.14914         86.12957    66.54732
    ## 2       PR   39.94039         35.94781         44.37642    18.62724
    ##   PFS_5yr_Lower_pct PFS_5yr_Upper_pct
    ## 1          56.17685          78.83222
    ## 2          15.55600          22.30483

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

![](Survival_analysis_lung_files/figure-gfm/-%205%20year%20KM%20curves-1.png)<!-- -->

``` r
p_os
```

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

![](Survival_analysis_lung_files/figure-gfm/-%205%20year%20KM%20curves-2.png)<!-- -->
