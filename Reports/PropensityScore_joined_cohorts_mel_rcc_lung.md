Propensity Scores Calculation and Survival Analysis
================
Mario Presti
First created on May 2025 Updated on 25 November 2025

- [Introduction](#introduction)
- [Loading Libraries](#loading-libraries)
  - [Pre-processing of data](#pre-processing-of-data)
- [Baseline Characteristics of Patients Before Propensity Score
  Matching](#baseline-characteristics-of-patients-before-propensity-score-matching)
- [Survival analysis using unmatched patient
  data](#survival-analysis-using-unmatched-patient-data)
  - [Overall survival](#overall-survival)
    - [Administrative censoring](#administrative-censoring)
    - [Combined survival analysis OS](#combined-survival-analysis-os)
      - [OS - paired analysis stratified by response before
        PSM](#os---paired-analysis-stratified-by-response-before-psm)
  - [Progression free survival](#progression-free-survival)
    - [Combined survival analysis PFS](#combined-survival-analysis-pfs)
    - [PFS - paired analysis stratified by response before
      PSM](#pfs---paired-analysis-stratified-by-response-before-psm)
- [Propensity Score Matching for Treatment Groups Based on Clinical
  Features - matching cohorts in pairs of
  2](#propensity-score-matching-for-treatment-groups-based-on-clinical-features---matching-cohorts-in-pairs-of-2)
  - [Matching with standardized caliper - using custom
    function](#matching-with-standardized-caliper---using-custom-function)
  - [PSM only on CR](#psm-only-on-cr)
  - [PSM only on PR](#psm-only-on-pr)
- [Sensitivity analysis - PSM on all using different
  calipers](#sensitivity-analysis---psm-on-all-using-different-calipers)
- [Median / Max follow up](#median--max-follow-up)
- [Percentages of administrative censored
  patients](#percentages-of-administrative-censored-patients)
  - [Administrative censoring rates per
    tumour](#administrative-censoring-rates-per-tumour)
- [Differential distribution of CRs](#differential-distribution-of-crs)
- [Post-hoc power analysis](#post-hoc-power-analysis)

# Introduction

This analysis investigates propensity score matching in patients based
on PD-L1 expression and other clinical parameters.

# Loading Libraries

``` r
# Load required libraries
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MatchIt)
library(survival)
library(survminer)
library(survey)
library(compareGroups)
library(cobalt)
library(ggsurvfit)
library(survivalAnalysis)
library(patchwork)
library(broom)
library(ggrepel)
library(openxlsx)
library(kableExtra)

source("PSM_survival_function_updated.R")
source("Survival_comparisons_function.R")

#this is to prevent problems with previous r versions
Tailcall <- function(f, x, char, symmetrical, recursive) {
  f(x, char, symmetrical, recursive)
}
```

## Pre-processing of data

- Check that all columns have no spaces etc, especially character
  columns,
- transform character columns to factors.
- Check if all factors are ok.
- Check numerical columns.

\#Joining melanoma and RCC datasets, plus data cleaning

``` r
setwd("E:/PhD_projects/Realworld/Data")
dataMel <- read.csv("Melanoma_data_polished.csv",check.names = F)
dataLung <- read.csv("Lung_data_polished.csv",check.names = F)
dataMel$Tumor <- "Melanoma"
dataRCC <- read.csv("RCC_data_polished.csv",check.names = F)
dataRCC$Tumor <- "RCC"
dataLung$Tumor <- "NSCLC"
common_cols <- intersect(names(dataMel), names(dataRCC))
common_cols <- intersect(names(dataLung),common_cols)
dataDF <- dataMel %>%
  bind_rows(
    dataRCC, dataLung)
dataDF <- dataDF[,common_cols]

colnames(dataDF) <- gsub("\\.", " ", colnames(dataDF))

names(dataDF)[names(dataDF) == 'Objective response'] <- 'Objective_response'

xpos <-  60 
y_cr <- 0.85  # height for the CR p‐value
y_pr <- 0.30  # height for the PR p‐value
```

We first need to define which columns we want as character
(categorical), and which need to be numerical. Sometimes things get
messed up with excel files.

``` r
# Define caliper 
caliper = 0.2

cat(paste0("Your will use a caliper of ", print(caliper)))
```

    ## [1] 0.2
    ## Your will use a caliper of 0.2

``` r
# Check for duplicates
any(duplicated(dataDF$record_id))
```

    ## [1] FALSE

``` r
# Remove duplicates
dataDF <- dataDF[!duplicated(dataDF), ]
dim(dataDF)
```

    ## [1] 2127   13

``` r
features <- c("Sex", "Age","Treatment line","Brain metastases", "ECOG Performance Status", 
              "Objective_response","Tumor")

## Define which columns include survival data
OS_time <- c("OS_days")
OS_status <- c("Dead")

PFS_time <- c("PFS_days")
PFS_status <- c("Progressed")

## DEFINE WHICH COLUMNS ARE CATEGORICAL
categ_feats <-c("Sex", "Treatment line", "Brain metastases", "ECOG Performance Status", "Objective_response","Tumor")

## DEFINE WHICH COLUMNS ARE NUMERICAL
# we also add the status and time columns
numeric_feats <- c("OS_days", "Dead","PFS_days","Progressed","Age")

OS_time <- gsub(" ","_",OS_time)
OS_status <- gsub(" ","_",OS_status)
PFS_time <- gsub(" ","_",PFS_time)
PFS_status <- gsub(" ","_",PFS_status)


## MAKE SURE THERE ARE NO BLANKS ANYWHERE IN THE DATA
## Find non-Numeric columns
num_cols <- unlist(lapply(dataDF, is.numeric))

### MAKE SURE THAT THE COLUMNS ARE IN THE CORRECT FORMAT NEEDED FOR ANALYSIS.
## Numeric columns
dataDF <- dataDF %>% mutate(across(all_of(numeric_feats), as.numeric))

## Categorical
dataDF <- dataDF %>% mutate(across(all_of(categ_feats), as.character))

#Transform character columns to categorical, to have the different levels
dataDF <- dataDF %>% mutate(across(all_of(categ_feats), as.factor))
dataDF <- dataDF %>% select(all_of(c(colnames(dataDF)[1], OS_time,OS_status, PFS_time,PFS_status, features, "CPI Regimen"))) %>% distinct()

print("Numbers per response")
```

    ## [1] "Numbers per response"

``` r
table(dataDF$bor_text)
```

    ## < table of extent 0 >

``` r
print("Numbers per treatment")
```

    ## [1] "Numbers per treatment"

``` r
table(dataDF$regime_correct)
```

    ## < table of extent 0 >

``` r
# Days to Months
dataDF$OS_months <- dataDF$OS_days /30.44
dataDF$PFS_months <- dataDF$PFS_days /30.44

## SUMMARY OF THE DATA
summary(dataDF)
```

    ##   record_id            OS_days            Dead           PFS_days     
    ##  Length:2127        Min.   :  42.0   Min.   :0.0000   Min.   :  42.0  
    ##  Class :character   1st Qu.: 668.5   1st Qu.:0.0000   1st Qu.: 330.5  
    ##  Mode  :character   Median :1211.0   Median :0.0000   Median : 726.0  
    ##                     Mean   :1385.8   Mean   :0.3926   Mean   : 969.0  
    ##                     3rd Qu.:2009.5   3rd Qu.:1.0000   3rd Qu.:1485.5  
    ##                     Max.   :4438.0   Max.   :1.0000   Max.   :4163.0  
    ##    Progressed         Sex            Age        Treatment line Brain metastases
    ##  Min.   :0.0000   Female: 932   Min.   :24.58   First:1697     No :1854        
    ##  1st Qu.:0.0000   Male  :1193   1st Qu.:60.13   Other: 430     Yes: 273        
    ##  Median :1.0000   NA's  :   2   Median :69.00                                  
    ##  Mean   :0.5811                 Mean   :67.44                                  
    ##  3rd Qu.:1.0000                 3rd Qu.:75.68                                  
    ##  Max.   :1.0000                 Max.   :94.10                                  
    ##  ECOG Performance Status Objective_response      Tumor      CPI Regimen       
    ##  PS=0:1101               CR: 728            Melanoma:1199   Length:2127       
    ##  PS≥1:1017               PR:1399            NSCLC   : 667   Class :character  
    ##  NA's:   9                                  RCC     : 261   Mode  :character  
    ##                                                                               
    ##                                                                               
    ##                                                                               
    ##    OS_months        PFS_months    
    ##  Min.   :  1.38   Min.   :  1.38  
    ##  1st Qu.: 21.96   1st Qu.: 10.86  
    ##  Median : 39.78   Median : 23.85  
    ##  Mean   : 45.53   Mean   : 31.83  
    ##  3rd Qu.: 66.02   3rd Qu.: 48.80  
    ##  Max.   :145.80   Max.   :136.76

``` r
names(dataDF)[names(dataDF) == 'Objective_response'] <- 'DOR'
features <- c("Sex", "Age","Treatment line","Brain metastases", "ECOG Performance Status", 
              "DOR","Tumor")

dataDF_unfilter<-dataDF

dataDF <- dataDF %>%
  mutate(across(
    where(is.character),
    ~ str_replace_all(., "\\bNA\\b", "")  
    )) %>%
  filter(if_all(everything(), ~ !is.na(.)))
```

# Baseline Characteristics of Patients Before Propensity Score Matching

``` r
common_cols <- intersect(names(dataMel), names(dataRCC))
lung_specific_cols <- setdiff(names(dataLung),common_cols)
melanoma_specific_cols <- setdiff(names(dataMel),c(common_cols, "Dead_Mel", "PS"))
renal_specific_cols <- setdiff(names(dataRCC),c(common_cols, "Smoking status")) #removing smoking status as it's only present in renal and not the others

descr_features <- paste0("`", c(features[ features != "Tumor" ], "CPI Regimen"), "`", collapse = " + ")
descr_features <- gsub("DOR","Objective Response",descr_features)
descr_formula <- as.formula(paste0("Tumor ~", descr_features))

dataDF_unfilter$Tumor <- factor(dataDF_unfilter$Tumor, levels = c("Melanoma", "NSCLC", "RCC"))

dataDF_extended <- full_join(dataDF_unfilter,dataMel[,c("record_id",melanoma_specific_cols)],
                             by="record_id")
dataDF_extended <- full_join(dataDF_extended,dataLung[,c("record_id",lung_specific_cols)],
                              by="record_id")
dataDF_extended <- full_join(dataDF_extended,dataRCC[,c("record_id",renal_specific_cols)],,
                              by="record_id")

#some little adjustments here
names(dataDF_extended)[names(dataDF_extended) == 'Histology subtype'] <- 'NSCLC subtype'
names(dataLung)[names(dataLung) == 'Histology subtype'] <- 'NSCLC subtype'
lung_specific_cols[which(lung_specific_cols == 'Histology subtype')] <- 'NSCLC subtype'
names(dataDF_extended)[names(dataDF_extended) == 'DOR'] <- 'Objective Response'

descr_features_ext <- paste0("`", c(features[ features != "Tumor" ],"CPI Regimen",
                                    melanoma_specific_cols,lung_specific_cols,renal_specific_cols), "`", collapse = " + ")
descr_features_ext <- gsub("DOR","Objective Response",descr_features_ext)
descr_formula_ext <- as.formula(paste0("Tumor ~", descr_features_ext))

#transform all categorical variables into factors

#now a bit of clean up of the table to prevent "missing = 100%" where is not needed, i.e. have the % of missingness only in tumor-specific columns actually missing data
na_counts_all <- colSums(is.na(dataDF_unfilter))
# names of columns with at least one NA

na_counts_mel <- colSums(is.na(dataMel[,melanoma_specific_cols]))
# names of columns with at least one NA
cols_with_na_mel <- names(na_counts_mel)[na_counts_mel > 0]

na_counts_rcc <- colSums(is.na(dataRCC[,renal_specific_cols]))
# names of columns with at least one NA
cols_with_na_rcc <- names(na_counts_rcc)[na_counts_rcc > 0]

na_counts_lung <- colSums(is.na(dataLung[,lung_specific_cols]))
# names of columns with at least one NA
cols_with_na_lung <- names(na_counts_lung)[na_counts_lung > 0]

cols_with_na <- c(cols_with_na_lung, cols_with_na_mel, cols_with_na_rcc)

dataDF_extended <- dataDF_extended %>%
  mutate(across(all_of(cols_with_na), as.factor))

#reorder for consistency 
dataDF_extended$Tumor <- factor(dataDF_extended$Tumor, levels = c("Melanoma","RCC","NSCLC"))

Baseline_Characteristics_table_Partial_patients <- descrTable(descr_formula_ext , data = dataDF_extended,include.miss=T, method=2, lab.missing = "Missing data", show.p.mul=T, show.p.overall=F)

Baseline_Characteristics_pvals <- descrTable(descr_formula , data = dataDF_extended,include.miss=T, method=2, lab.missing = "Missing data", show.p.mul=T, show.p.overall=F)
Baseline_Characteristics_pvals$descr <-  Baseline_Characteristics_table_Partial_patients$descr[,4:6]

Baseline_Characteristics_table_Partial_patients$descr <- Baseline_Characteristics_table_Partial_patients$descr[,1:3]
#Baseline_Characteristics_table_Partial_patients, file = file.path(path.expand(result), "Baseline_Characteristics_table_all_patients.docx"))
index_na <- which(sub(":.*$", "", rownames(Baseline_Characteristics_table_Partial_patients$descr)) %in% cols_with_na)
descr <- Baseline_Characteristics_table_Partial_patients$descr

# compute mask on subset
sub <- descr[index_na, , drop = FALSE]
mask <- str_detect(
  sub,
  regex("^\\s*(?:0 \\(\\s*0(?:\\.0+)?%\\s*\\)|\\d+ \\(\\s*100(?:\\.0+)?%\\s*\\))$")
)

# embed into full-size mask
full_mask <- matrix(FALSE, nrow = nrow(descr), ncol = ncol(descr))
full_mask[index_na, ] <- mask

# assign NA where TRUE
descr[full_mask] <- "NA"
Baseline_Characteristics_table_Partial_patients$descr <- descr


Baseline_Characteristics_table_Partial_patients$descr <- gsub("0 *\\(\\.?%\\)", "NA", Baseline_Characteristics_table_Partial_patients$descr)

truncate_pct <- function(x) {
  str_replace_all(x, "\\(([0-9]+\\.?[0-9]*)%\\)", function(m) {
    num <- as.numeric(sub("^\\((.*)%\\)$", "\\1", m))
    truncated <- floor(num * 10) / 10
    sprintf("(%0.1f%%)", truncated)
  })
}

Baseline_Characteristics_table_Partial_patients$descr <-
  apply(Baseline_Characteristics_table_Partial_patients$descr,MARGIN = c(1,2), truncate_pct)

export2word(Baseline_Characteristics_table_Partial_patients, caption = "Table 1. Summary of the patient characteristics", file = "E:/PhD_projects/Realworld/Scripts/Acquired_resistance_RW/Tables/Table_1.docx") 
export2word(Baseline_Characteristics_pvals, caption = "Supplementary table 1. Pairwise comparisons among baseline features in the cohorts", file = "E:/PhD_projects/Realworld/Scripts/Acquired_resistance_RW/Tables/Supplementary_table_1.docx") 

Baseline_Characteristics_table_Partial_patients
```

    ## 
    ## --------Summary descriptives table by 'Tumor'---------
    ## 
    ## _____________________________________________________________________________ 
    ##                                Melanoma           RCC             NSCLC       
    ##                                 N=1199           N=261            N=667       
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
    ## Sex:                                                                          
    ##     Female                   487 (40.6%)       59 (22.6%)      386 (57.9%)    
    ##     Male                     712 (59.4%)      200 (76.6%)      281 (42.1%)    
    ##     Missing data               0 (0.0%)         2 (0.7%)         0 (0.0%)     
    ## Age                        70.0 [59.1;76.9] 65.4 [58.4;72.3] 69.0 [63.0;75.0] 
    ## Treatment line:                                                               
    ##     First                    1031 (86.0%)     189 (72.4%)      477 (71.5%)    
    ##     Other                    168 (14.0%)       72 (27.6%)      190 (28.5%)    
    ## Brain metastases:                                                             
    ##     No                       1011 (84.3%)     250 (95.8%)      593 (88.9%)    
    ##     Yes                      188 (15.7%)       11 (4.2%)        74 (11.1%)    
    ## ECOG Performance Status:                                                      
    ##     PS=0                     787 (65.6%)      128 (49.0%)      186 (27.9%)    
    ##     PS≥1                     409 (34.1%)      130 (49.8%)      478 (71.7%)    
    ##     Missing data               3 (0.2%)         3 (1.1%)         3 (0.4%)     
    ## Objective Response:                                                           
    ##     CR                       576 (48.0%)       79 (30.3%)       73 (10.9%)    
    ##     PR                       623 (52.0%)      182 (69.7%)      594 (89.1%)    
    ## CPI Regimen:                                                                  
    ##     Anti-PD1                 789 (65.8%)       61 (23.4%)      667 (100.0%)   
    ##     Anti-PD1+Anti-CTLA4      410 (34.2%)      200 (76.6%)        0 (0.0%)     
    ## Previous adjuvant therapy:                                                    
    ##     No                       1098 (91.6%)          NA               NA        
    ##     Yes                       101 (8.4%)           NA               NA        
    ## Melanoma subtype:                                                             
    ##     Cutaneous                949 (79.1%)           NA               NA        
    ##     Mucosal                   38 (3.1%)            NA               NA        
    ##     Unk. primary             212 (17.7%)           NA               NA        
    ## BRAF mutation:                                                                
    ##     Mutant                   504 (42.0%)           NA               NA        
    ##     Wild type                667 (55.6%)           NA               NA        
    ##     Missing data              28 (2.3%)            NA               NA        
    ## AJCC 8th stage:                                                               
    ##     III/M1a/M1b              503 (42.0%)           NA               NA        
    ##     M1c/M1d                  696 (58.0%)           NA               NA        
    ## LDH:                                                                          
    ##     Elevated                 339 (28.3%)           NA               NA        
    ##     Normal                   821 (68.5%)           NA               NA        
    ##     Missing data              39 (3.2%)            NA               NA        
    ## PD-L1 status:                                                                 
    ##     PDL1<50%                      NA               NA           87 (13.0%)    
    ##     PDL1≥50%                      NA               NA          546 (81.9%)    
    ##     Missing data                  NA               NA           34 (5.1%)     
    ## NSCLC subtype:                                                                
    ##     Nonsquamous                   NA               NA          546 (81.9%)    
    ##     Squamous                      NA               NA          121 (18.1%)    
    ## Dead_RCC                       . [.;.]      0.00 [0.00;1.00]     . [.;.]      
    ## RCC subtype:                                                                  
    ##     Clear cell                    NA          225 (86.2%)           NA        
    ##     Non clear cell                NA           36 (13.8%)           NA        
    ## Sarcomatoid subtype:                                                          
    ##     Non sarcomatoid               NA          182 (69.7%)           NA        
    ##     Sarcomatoid                   NA           79 (30.3%)           NA        
    ## IMDC:                                                                         
    ##     Good/Intermediate             NA          187 (71.6%)           NA        
    ##     Poor                          NA           74 (28.4%)           NA        
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

# Survival analysis using unmatched patient data

## Overall survival

### Administrative censoring

``` r
# Administrative censoring at 60 months
dataDF$censor_time_OS   <- pmin(dataDF$OS_months, 60)
dataDF$censor_status_OS <- ifelse(dataDF$OS_months > 60, 0, dataDF$Dead)
```

### Combined survival analysis OS

``` r
# Compute the Kaplan-Meier model
fit_overall = survfit(Surv(dataDF$censor_time_OS, dataDF$censor_status_OS) ~ Tumor+DOR, data=dataDF)
print(fit_overall)
```

    ## Call: survfit(formula = Surv(dataDF$censor_time_OS, dataDF$censor_status_OS) ~ 
    ##     Tumor + DOR, data = dataDF)
    ## 
    ##                          n events median 0.95LCL 0.95UCL
    ## Tumor=Melanoma, DOR=CR 576     57     NA      NA      NA
    ## Tumor=Melanoma, DOR=PR 620    251   46.0    38.7    55.4
    ## Tumor=NSCLC, DOR=CR     73     17     NA      NA      NA
    ## Tumor=NSCLC, DOR=PR    591    334   42.6    38.2    48.4
    ## Tumor=RCC, DOR=CR       79      3     NA      NA      NA
    ## Tumor=RCC, DOR=PR      178     75     NA    50.6      NA

``` r
# Modify the names of the strata to remove "regime_correct="
names(fit_overall$strata) <- sub("Tumor=", "", names(fit_overall$strata))
names(fit_overall$strata) <- sub("DOR=", "", names(fit_overall$strata))

# Cox model for HR & CI
dataDF$Tumor <- relevel(dataDF$Tumor, ref = "Melanoma")
dataDF$DOR <- relevel(dataDF$DOR, ref = "CR")
fit_cox <- coxph(Surv(dataDF$censor_time_OS, dataDF$censor_status_OS) ~ Tumor+DOR, data=dataDF)
s_cox   <- summary(fit_cox)
hr      <- s_cox$coefficients[, "exp(coef)"]
ci_lo   <- s_cox$conf.int[, "lower .95"]
ci_hi   <- s_cox$conf.int[, "upper .95"]
hr_label <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)",
                        hr, ci_lo, ci_hi)

# Pull out the table of HRs+CIs
tidy_cox <- broom::tidy(
  fit_cox,
  exponentiate = TRUE,  # get HRs
  conf.int     = TRUE   # get 95% CIs
)

# Prepare input for forest_plot.df()
tbl <- tidy_cox %>%
  mutate(
    # identify which factor each term belongs to:
    factor.name  = case_when(
      startsWith(term, "Tumor")              ~ "Tumor",
      startsWith(term, "DOR") ~ "DOR",
      TRUE                                   ~ NA_character_
    ),
    # pull out the level name by dropping the var name prefix
    factor.value = mapply(function(t, f) substring(t, nchar(f) + 1),
                          term, factor.name,
                          USE.NAMES = FALSE),
    # define the rest of the required columns:
    survivalResult = "OS",
    endpoint       = "OS_months",
    factor.id      = paste0(factor.name, ":", factor.value),
    HR             = estimate,
    Lower_CI       = conf.low,
    Upper_CI       = conf.high,
    p              = p.value,
    n              = NA_integer_,
    subgroup_n     = NA_integer_
  ) %>%
  select(survivalResult, endpoint,
         factor.name, factor.value, factor.id,
         HR, Lower_CI, Upper_CI, p,
         n, subgroup_n)

# Build the reference‐level rows dynamically
refs <- tbl %>%
  distinct(factor.name) %>%
  rowwise() %>%
  dplyr::mutate(
    # all levels of this factor:
    all_lvls = list(levels(dataDF[[factor.name]])),
    # those present in the Cox output:
    used_lvls = list(tbl$factor.value[tbl$factor.name == factor.name]),
    # the one leftover is the reference:
    factor.value = setdiff(all_lvls[[1]], used_lvls[[1]]),
    survivalResult = tbl$survivalResult[1],
    endpoint       = tbl$endpoint[1],
    factor.id      = paste0(factor.name, ":", factor.value),
    HR             = 1,
    Lower_CI       = 1,
    Upper_CI       = 1,
    p              = NA_real_,
    n              = sum(dataDF[[factor.name]] == factor.value),
    subgroup_n     = n
  ) %>%
  ungroup() %>%
  select(names(tbl))

# Combine REF + non‐REF, ordering by factor.name
full_tbl <- bind_rows(refs, tbl) %>%
  arrange(factor.name,
          factor.value)


OVS_all <- ggsurvfit(fit_overall, size = 1.5)
p_km <- OVS_all +
  add_censor_mark() +
  add_risktable(
    size = 5,
    theme = theme_risktable_default(
      axis.text.y.size = 15,
      plot.title.size  = 15
    ),
    risktable_stats = "{n.risk} ({cum.event})"
  ) +
  labs(
    x     = "Months",
    y     = "Survival probability (%)",
    title = "Overall Survival before propensity matching"
  ) +
  xlim(0, 60) +
  scale_x_continuous(breaks = seq(0, 60, by = 12)) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1.0), limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 70), expand = FALSE) +
  theme_classic() +
  theme(
    plot.title        = element_text(hjust = 0.5, size = 18),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    axis.text.x       = element_text(size = 15),
    axis.text.y       = element_text(size = 15),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 18),
    legend.key.size   = unit(10, "bigpts"),
    legend.title      = element_blank()
  )

HR_x_limits <- c(min(full_tbl$HR), max(full_tbl$HR)) * c(1/1.5, 1.5)
HR_x_breaks <- round(10^pretty(log10(HR_x_limits), n = 7),digits = 1)

full_tbl$factor.id <- sub("Tumor:", "", full_tbl$factor.id)
full_tbl$factor.id <- sub("DOR:", "", full_tbl$factor.id)
## Forest plot panel with REF included
p_forest2 <- survivalAnalysis:::forest_plot.df(
  full_tbl,
  labels_displayed = c("factor"),
  values_displayed = c("HR", "CI"),
  HR_x_limits = HR_x_limits,
  HR_x_breaks = HR_x_breaks
) +
  labs(x = "Hazard Ratio") +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 12)
  )

# Patchwork them together
final_plot <- p_km + p_forest2 +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "OS curves with corresponding HRs (including reference)"
  )

print(final_plot)
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Survival%20analysis%20OS-1.png)<!-- -->

#### OS - paired analysis stratified by response before PSM

``` r
comps <- list(
  c("Melanoma",   "RCC"),
  c("Melanoma",   "NSCLC"),
  c("RCC","NSCLC")
)

res_OS_CR <- run_survival_comparisons(
  subset(dataDF,DOR=="CR"),
  comparisons       = comps,
  covariate_to_split= "Tumor",
  descr_formula     = NULL,
  endpoint          = "OS",
  truncate_month    = 60
)

OS_CR <- res_OS_CR[[3]]
OS_CR +
  plot_annotation(
    title    = "Overall Survival in Complete Responders",
    subtitle = "Before PSM",
    theme    = theme(
      plot.title    = element_text(size = 40, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 30,, face = "italic", hjust = 0.5)
    )
  )
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Pairwise%20OS%20analysis-1.png)<!-- -->

``` r
res_OS_PR <- run_survival_comparisons(
  dataDF = subset(dataDF,DOR=="PR"),
  comparisons       = comps,
  covariate_to_split= "Tumor",
  descr_formula     = NULL,
  endpoint          = "OS",
  truncate_month    = 60
)

OS_PR <- res_OS_PR[[3]]
OS_PR +
  plot_annotation(
    title    = "Overall Survival in Partial Responders",
    subtitle = "Before PSM",
    theme    = theme(
      plot.title    = element_text(size = 40, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 30,, face = "italic", hjust = 0.5)
    )
  )
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Pairwise%20OS%20analysis-2.png)<!-- -->

## Progression free survival

### Combined survival analysis PFS

``` r
# Administrative censoring at 60 months
dataDF$censor_time_PFS   <- pmin(dataDF$PFS_months, 60)
dataDF$censor_status_PFS <- ifelse(dataDF$PFS_months > 60, 0, dataDF$Progressed)

dataDF$Tumor <- relevel(dataDF$Tumor, ref = "Melanoma")
dataDF$DOR <- relevel(dataDF$DOR, ref = "CR")
# Kaplan–Meier on the censored data
fit_PFS <- survfit(
  Surv(censor_time_PFS, censor_status_PFS) ~ Tumor + DOR,
  data = dataDF
)
# clean up strata names
names(fit_PFS$strata) <- sub("Tumor=",             "", names(fit_PFS$strata))
names(fit_PFS$strata) <- sub("DOR=", "", names(fit_PFS$strata))

# Plot
# Cox model for HR & CI
dataDF$Tumor <- relevel(dataDF$Tumor, ref = "Melanoma")
dataDF$DOR <- relevel(dataDF$DOR, ref = "CR")
fit_cox <- coxph(Surv(dataDF$censor_time_PFS, dataDF$censor_status_PFS) ~ Tumor+DOR, data=dataDF)
s_cox   <- summary(fit_cox)
hr      <- s_cox$coefficients[, "exp(coef)"]
ci_lo   <- s_cox$conf.int[, "lower .95"]
ci_hi   <- s_cox$conf.int[, "upper .95"]
hr_label <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)",
                        hr, ci_lo, ci_hi)

# Pull out the table of HRs+CIs
tidy_cox <- broom::tidy(
  fit_cox,
  exponentiate = TRUE,  # get HRs
  conf.int     = TRUE   # get 95% CIs
)

# Prepare input for forest_plot.df()
tbl <- tidy_cox %>%
  mutate(
    # identify which factor each term belongs to:
    factor.name  = case_when(
      startsWith(term, "Tumor")              ~ "Tumor",
      startsWith(term, "DOR") ~ "DOR",
      TRUE                                   ~ NA_character_
    ),
    # pull out the level name by dropping the var name prefix
    factor.value = mapply(function(t, f) substring(t, nchar(f) + 1),
                          term, factor.name,
                          USE.NAMES = FALSE),
    # define the rest of the required columns:
    survivalResult = "PFS",
    endpoint       = "PFS_months",
    factor.id      = paste0(factor.name, ":", factor.value),
    HR             = estimate,
    Lower_CI       = conf.low,
    Upper_CI       = conf.high,
    p              = p.value,
    # fill n & subgroup_n in the next step
    n              = NA_integer_,
    subgroup_n     = NA_integer_
  ) %>%
  select(survivalResult, endpoint,
         factor.name, factor.value, factor.id,
         HR, Lower_CI, Upper_CI, p,
         n, subgroup_n)

# Build the reference‐level rows dynamically
refs <- tbl %>%
  distinct(factor.name) %>%
  rowwise() %>%
  dplyr::mutate(
    # all levels of this factor:
    all_lvls = list(levels(dataDF[[factor.name]])),
    # those present in the Cox output:
    used_lvls = list(tbl$factor.value[tbl$factor.name == factor.name]),
    # the one leftover is the reference:
    factor.value = setdiff(all_lvls[[1]], used_lvls[[1]]),
    survivalResult = tbl$survivalResult[1],
    endpoint       = tbl$endpoint[1],
    factor.id      = paste0(factor.name, ":", factor.value),
    HR             = 1,
    Lower_CI       = 1,
    Upper_CI       = 1,
    p              = NA_real_,
    n              = sum(dataDF[[factor.name]] == factor.value),
    subgroup_n     = n
  ) %>%
  ungroup() %>%
  select(names(tbl))

# Combine REF + non‐REF, ordering by factor.name
full_tbl <- bind_rows(refs, tbl) %>%
  arrange(factor.name,
          factor.value)


PFS_all <- ggsurvfit(fit_PFS, size = 1.5)
p_km <- PFS_all +
  add_censor_mark() +
  add_risktable(
    size = 5,
    theme = theme_risktable_default(
      axis.text.y.size = 15,
      plot.title.size  = 15
    ),
    risktable_stats = "{n.risk} ({cum.event})"
  ) +
  labs(
    x     = "Months",
    y     = "Survival probability (%)",
    title = "Progression free survival before propensity matching"
  ) +
  xlim(0, 60) +
  scale_x_continuous(breaks = seq(0, 60, by = 12)) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1.0), limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 70), expand = FALSE) +
  theme_classic() +
  theme(
    plot.title        = element_text(hjust = 0.5, size = 18),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    axis.text.x       = element_text(size = 15),
    axis.text.y       = element_text(size = 15),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 18),
    legend.key.size   = unit(10, "bigpts"),
    legend.title      = element_blank()
  )

HR_x_limits <- c(min(full_tbl$HR), max(full_tbl$HR)) * c(1/1.5, 1.5)
HR_x_breaks <- round(10^pretty(log10(HR_x_limits), n = 7),digits = 1)

full_tbl$factor.id <- sub("Tumor:", "", full_tbl$factor.id)
full_tbl$factor.id <- sub("DOR:", "", full_tbl$factor.id)
## Forest plot panel with REF included
p_forest2 <- survivalAnalysis:::forest_plot.df(
  full_tbl,
  labels_displayed = c("factor"),
  values_displayed = c("HR", "CI", "p"),
  HR_x_limits = HR_x_limits,
  HR_x_breaks = HR_x_breaks
) +
  labs(x = "Hazard Ratio") +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 12)
  )

# Patchwork them together
final_plot <- p_km + p_forest2 +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "PFS curves with corresponding HRs (including reference)"
  )

print(p_km)
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Survival%20analysis%20-%20PFS-1.png)<!-- -->

``` r
print(final_plot)
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Survival%20analysis%20-%20PFS-2.png)<!-- -->

``` r
#perform pairwise comparisons among the PFS curves
pairs <- combn(unique(dataDF$Tumor), 2, simplify = FALSE)
PR_results <- lapply(pairs, function(pair) {
  # subset to PR & these two tumours
  subdf <- dataDF %>%
    filter(DOR == "PR", Tumor %in% pair) %>%
    mutate(Tumor = factor(Tumor, levels = pair))
  
  # fit Cox on PFS
  fit <- coxph(Surv(censor_time_PFS, censor_status_PFS) ~ Tumor, data = subdf)
  s   <- summary(fit)
  
  # pull out the one coefficient
  hr      <- s$coefficients[1, "exp(coef)"]
  ci_lo   <- s$conf.int[1, "lower .95"]
  ci_hi   <- s$conf.int[1, "upper .95"]
  pval    <- s$coefficients[1, "Pr(>|z|)"]
  
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    HR         = hr,
    Lower_95   = ci_lo,
    Upper_95   = ci_hi,
    p.value    = pval,
    stringsAsFactors = FALSE
  )
})

# Combine and adjust
PR_pvals <- bind_rows(PR_results) %>%
  mutate(p.adj = p.adjust(p.value, method = "bonferroni"))
PR_pvals
```

    ##          comparison        HR  Lower_95  Upper_95   p.value     p.adj
    ## 1   Melanoma vs RCC 0.8646867 0.7092980 1.0541171 0.1502914 0.4508741
    ## 2 Melanoma vs NSCLC 0.8674949 0.7580637 0.9927233 0.0388172 0.1164516
    ## 3      RCC vs NSCLC 1.0253783 0.8444288 1.2451029 0.8002748 1.0000000

``` r
##Now for CR patients
pairs <- combn(unique(dataDF$Tumor), 2, simplify = FALSE)
CR_results <- lapply(pairs, function(pair) {
  # subset to PR & these two tumours
  subdf <- dataDF %>%
    filter(DOR == "CR", Tumor %in% pair) %>%
    mutate(Tumor = factor(Tumor, levels = pair))
  
  # fit Cox on PFS
  fit <- coxph(Surv(censor_time_PFS, censor_status_PFS) ~ Tumor, data = subdf)
  s   <- summary(fit)
  
  # pull out the one coefficient
  hr      <- s$coefficients[1, "exp(coef)"]
  ci_lo   <- s$conf.int[1, "lower .95"]
  ci_hi   <- s$conf.int[1, "upper .95"]
  pval    <- s$coefficients[1, "Pr(>|z|)"]
  
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    HR         = hr,
    Lower_95   = ci_lo,
    Upper_95   = ci_hi,
    p.value    = pval,
    stringsAsFactors = FALSE
  )
})

# Combine and adjust
CR_pvals <- bind_rows(CR_results) %>%
  mutate(p.adj = p.adjust(p.value, method = "bonferroni"))
CR_pvals
```

    ##          comparison       HR  Lower_95 Upper_95   p.value     p.adj
    ## 1   Melanoma vs RCC 1.123322 0.7073093 1.784018 0.6222048 1.0000000
    ## 2 Melanoma vs NSCLC 1.298940 0.8391503 2.010659 0.2406806 0.7220419
    ## 3      RCC vs NSCLC 1.067709 0.5908494 1.929430 0.8281998 1.0000000

``` r
### Median time to progression in PR vs CR
fit_PFS <- survfit(
  Surv(censor_time_PFS, censor_status_PFS) ~ DOR,
  data = dataDF
)

print(final_plot)
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Survival%20analysis%20-%20PFS-3.png)<!-- -->

``` r
PFS_all +
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
    add_risktable_strata_symbol(symbol = "•", size = 20) +
    labs(
        x        = "Months after treatment initiation",
        y        = paste0("PFS", "(%)")
    ) +
    scale_x_continuous(breaks = seq(0, 60, by = 12), limits = c(0, 60)) +
    theme_classic() +
    theme(
        plot.title      = element_text(hjust = 0.5, size = 18),
        plot.subtitle   = element_text(hjust = 0.5, size = 30),
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
    labs(
    title = "Progression free survival before propensity matching"
  ) 
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Survival%20analysis%20-%20PFS-4.png)<!-- -->

### PFS - paired analysis stratified by response before PSM

``` r
comps <- list(
  c("Melanoma",   "RCC"),
  c("Melanoma",   "NSCLC"),
  c("RCC","NSCLC")
)

res_PFS_CR <- run_survival_comparisons(
  subset(dataDF,DOR=="CR"),
  comparisons       = comps,
  covariate_to_split= "Tumor",
  descr_formula     = NULL,
  endpoint          = "PFS",
  truncate_month    = 60
)


res_PFS_CR[[3]] +
  plot_annotation(
    title    = "Progression-Free Survival in Complete Responders",
    subtitle = "Before PSM",
    theme    = theme(
      plot.title    = element_text(size = 40, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 30,, face = "italic", hjust = 0.5)
    )
  )
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Pairwise%20PFS-1.png)<!-- -->

``` r
res_PFS_PR <- run_survival_comparisons(dataDF = subset(dataDF,DOR=="PR"),
  comparisons       = comps,
  covariate_to_split= "Tumor",
  descr_formula     = NULL,
  endpoint          = "PFS",
  truncate_month    = 60
)

PFS_PR <- res_PFS_PR[[3]]
PFS_PR +
  plot_annotation(
    title    = "Progression-Free Survival in Partial Responders",
    subtitle = "Before PSM",
    theme    = theme(
      plot.title    = element_text(size = 40, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 30,, face = "italic", hjust = 0.5)
    )
  )
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/Pairwise%20PFS-2.png)<!-- -->

# Propensity Score Matching for Treatment Groups Based on Clinical Features - matching cohorts in pairs of 2

## Matching with standardized caliper - using custom function

``` r
matching_features <- paste0("`", setdiff(features, c("CPI Regimen", "Brain metastases", "Tumor")), "`", collapse = " + ")
matching_formula <- as.formula(paste0("Tumor ~", matching_features))

descr_features <- paste0("`", setdiff(features, c("CPI Regimen", "DOR", "Brain metastases")), "`", collapse = " + ")
descr_formula <- as.formula(paste0("DOR ~", descr_features))
comparisons <- list(c("Melanoma", "RCC"), c("Melanoma", "NSCLC"), c("RCC", "NSCLC"))

OS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons,matching_formula =  matching_formula, covariate_to_split = "Tumor",descr_formula = NULL, endpoint = "OS", truncate_month = 60, caliper=0.2)
PFS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons, matching_formula, covariate_to_split = "Tumor",descr_formula = NULL, endpoint = "PFS", truncate_month = 60, caliper = 0.2)

print(paste0("PS matching with the following formula:", paste0("Tumor ~", matching_features)))
```

    ## [1] "PS matching with the following formula:Tumor ~`Sex` + `Age` + `Treatment line` + `ECOG Performance Status` + `DOR`"

``` r
OS_psm_all$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20all-1.png)<!-- -->

``` r
OS_psm_all$love_plots
```

    ## $`Melanoma vs RCC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20all-2.png)<!-- -->

    ## 
    ## $`Melanoma vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20all-3.png)<!-- -->

    ## 
    ## $`RCC vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20all-4.png)<!-- -->

``` r
PFS_psm_all$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20all-5.png)<!-- -->

``` r
ratio = 2
OS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons, matching_formula, covariate_to_split = "Tumor",descr_formula = NULL, endpoint = "OS", truncate_month = 60, caliper=0.2,ratio = ratio)
PFS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons, matching_formula, covariate_to_split = "Tumor",descr_formula = NULL, endpoint = "PFS", truncate_month = 60, caliper = 0.2,ratio = ratio)

print(paste0("PS matching with the following formula:", paste0("Tumor ~", matching_features), "and ratio: ", ratio))
```

    ## [1] "PS matching with the following formula:Tumor ~`Sex` + `Age` + `Treatment line` + `ECOG Performance Status` + `DOR`and ratio: 2"

``` r
OS_psm_all$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-1.png)<!-- -->

``` r
OS_psm_all$love_plots
```

    ## $`Melanoma vs RCC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-2.png)<!-- -->

    ## 
    ## $`Melanoma vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-3.png)<!-- -->

    ## 
    ## $`RCC vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-4.png)<!-- -->

``` r
PFS_psm_all$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-5.png)<!-- -->

``` r
ratio = 3
OS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons, matching_formula, covariate_to_split = "Tumor",descr_formula = NULL, endpoint = "OS", truncate_month = 60, caliper=0.2,ratio = ratio)
PFS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons, matching_formula, covariate_to_split = "Tumor",descr_formula = NULL, endpoint = "PFS", truncate_month = 60, caliper = 0.2,ratio = ratio)

print(paste0("PS matching with the following formula:", paste0("Tumor ~", matching_features), "and ratio: ", ratio))
```

    ## [1] "PS matching with the following formula:Tumor ~`Sex` + `Age` + `Treatment line` + `ECOG Performance Status` + `DOR`and ratio: 3"

``` r
OS_psm_all$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-6.png)<!-- -->

``` r
OS_psm_all$love_plots
```

    ## $`Melanoma vs RCC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-7.png)<!-- -->

    ## 
    ## $`Melanoma vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-8.png)<!-- -->

    ## 
    ## $`RCC vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-9.png)<!-- -->

``` r
PFS_psm_all$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/PSM%20using%202%20to%201%20or%203%20to%201-10.png)<!-- -->

## PSM only on CR

``` r
matching_features_no_DOR <- paste0("`", setdiff(features, c("CPI Regimen", "Brain metastases", "DOR", "Tumor")), "`", collapse = " + ")
matching_formula_no_DOR <- as.formula(paste0("Tumor ~", matching_features_no_DOR))

comparisons <- list(c("Melanoma", "RCC"), c("Melanoma", "NSCLC"), c("RCC", "NSCLC"))

OS_psm_CR <- run_psm_and_survival(dataDF = subset(dataDF, DOR == "CR"), matching_formula = matching_formula_no_DOR, comparisons = comparisons,covariate_to_split = "Tumor", endpoint = "OS", truncate_month = 60)
PFS_psm_CR <- run_psm_and_survival(dataDF = subset(dataDF, DOR == "CR"),matching_formula= matching_formula_no_DOR, comparisons = comparisons,covariate_to_split = "Tumor", endpoint = "PFS", truncate_month = 60)

print(paste0("PS matching with the following formula:", paste0("Tumor ~", matching_features)))
```

    ## [1] "PS matching with the following formula:Tumor ~`Sex` + `Age` + `Treatment line` + `ECOG Performance Status` + `DOR`"

``` r
OS_psm_CR$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-1.png)<!-- -->

``` r
PFS_psm_CR$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-2.png)<!-- -->

``` r
OS_psm_CR$love_plots
```

    ## $`Melanoma vs RCC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-3.png)<!-- -->

    ## 
    ## $`Melanoma vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-4.png)<!-- -->

    ## 
    ## $`RCC vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-5.png)<!-- -->

``` r
PFS_psm_CR$love_plots
```

    ## $`Melanoma vs RCC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-6.png)<!-- -->

    ## 
    ## $`Melanoma vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-7.png)<!-- -->

    ## 
    ## $`RCC vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20CR-8.png)<!-- -->

## PSM only on PR

``` r
matching_features <- paste0("`", setdiff(features, c("CPI Regimen", "Brain metastases", "DOR", "Tumor")), "`", collapse = " + ")
matching_formula <- as.formula(paste0("Tumor ~", matching_features))

comparisons <- list(c("Melanoma", "RCC"), c("Melanoma", "NSCLC"), c("RCC", "NSCLC"))
matching_formula
```

    ## Tumor ~ Sex + Age + `Treatment line` + `ECOG Performance Status`

``` r
OS_psm_PR <- run_psm_and_survival(dataDF = subset(dataDF, DOR == "PR"), comparisons = comparisons, matching_formula = matching_formula_no_DOR, covariate_to_split = "Tumor", endpoint = "OS", truncate_month = 60)
PFS_psm_PR <- run_psm_and_survival(dataDF = subset(dataDF, DOR == "PR"), comparisons = comparisons, matching_formula = matching_formula_no_DOR,covariate_to_split = "Tumor", endpoint = "PFS", truncate_month = 60)

print(paste0("PS matching with the following formula:", paste0("Tumor ~", matching_features)))
```

    ## [1] "PS matching with the following formula:Tumor ~`Sex` + `Age` + `Treatment line` + `ECOG Performance Status`"

``` r
OS_psm_PR$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20PR-1.png)<!-- -->

``` r
PFS_psm_PR$combined_surv_plot
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20PR-2.png)<!-- -->

# Sensitivity analysis - PSM on all using different calipers

``` r
matching_features <- paste0("`", setdiff(features, c("CPI Regimen", "Brain metastases", "Tumor")), "`", collapse = " + ")
matching_formula <- as.formula(paste0("Tumor ~", matching_features))

descr_features <- paste0("`", setdiff(features, c("CPI Regimen", "DOR", "Brain metastases")), "`", collapse = " + ")
descr_formula <- as.formula(paste0("DOR ~", descr_features))
comparisons <- list(c("Melanoma", "RCC"), c("Melanoma", "NSCLC"), c("RCC", "NSCLC"))

calipers <- c(0.05,0.1,0.2,0.4,0.8)
for (caliper in calipers){
  PFS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons,
                                      matching_formula, covariate_to_split = "Tumor",
                                      descr_formula = NULL, endpoint = "PFS",
                                      truncate_month = 60, 
                                      caliper = caliper)

  plot_title <- paste0("PS matching with the full formula and caliper:", paste0(caliper))
  print(PFS_psm_all$combined_surv_plot+
          plot_annotation(title=plot_title,
                          theme=theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5),
                                      plot.subtitle = element_text(size = 15, hjust = 0.5)
    )
  )
  )
}
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-1.png)<!-- -->![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-2.png)<!-- -->![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-3.png)<!-- -->![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-4.png)<!-- -->![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-5.png)<!-- -->

``` r
PFS_psm_all <- run_psm_and_survival(dataDF = dataDF, comparisons = comparisons,
                                      matching_formula, covariate_to_split = "Tumor",
                                      descr_formula = NULL, endpoint = "PFS",
                                      truncate_month = 60, 
                                      caliper = 0.2, 
                                    mahvars = "Age")
PFS_psm_all$combined_surv_plot+plot_annotation(title="PSM using mahvars for Age",
                          theme=theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5),
                                      plot.subtitle = element_text(size = 15, hjust = 0.5)))
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-6.png)<!-- -->

``` r
PFS_psm_all$love_plots
```

    ## $`Melanoma vs RCC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-7.png)<!-- -->

    ## 
    ## $`Melanoma vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-8.png)<!-- -->

    ## 
    ## $`RCC vs NSCLC`

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/propensity%20matching%20different%20calipers-9.png)<!-- -->

# Median / Max follow up

``` r
for(i in unique(dataDF$Tumor)) {
  rev_surv_obj <- with(
    subset(dataDF, Tumor == i),
    Surv(censor_time_OS, censor_status_OS == 0)
  )

  # Fit the reversed Kaplan-Meier curve
  rev_fit <- survfit(rev_surv_obj ~ 1)

  # Grab the summary table
  tab <- summary(rev_fit)$table
  # It has named elements "median", "0.95LCL", "0.95UCL"
  med   <- tab["median"]
  lower <- tab["0.95LCL"]
  upper <- tab["0.95UCL"]

  # Convert from months to years (divide by 12)
  med_yr   <- med   / 12
  lower_yr <- lower / 12
  upper_yr <- upper / 12

  cat(sprintf(
    "Tumour: %s\n  Median follow-up = %.2f years (95%% CI %.2f–%.2f)\n",
    i, med_yr, lower_yr, upper_yr
  ))

  # And the maximum follow-up (no CI available for max)
  max_fu <- max(rev_fit$time, na.rm = TRUE) / 12
  cat(sprintf("  Maximum follow-up = %.2f years\n\n", max_fu))
}
```

    ## Tumour: Melanoma
    ##   Median follow-up = 4.92 years (95% CI 4.66–NA)
    ##   Maximum follow-up = 5.00 years
    ## 
    ## Tumour: RCC
    ##   Median follow-up = 4.50 years (95% CI 4.22–4.74)
    ##   Maximum follow-up = 5.00 years
    ## 
    ## Tumour: NSCLC
    ##   Median follow-up = 5.00 years (95% CI NA–NA)
    ##   Maximum follow-up = 5.00 years

``` r
# Fit the reverse‐KM: event = still under follow‐up (i.e. censor_status_OS == 0)
rev_fit_all <- survfit(
  Surv(censor_time_OS, censor_status_OS == 0) ~ Tumor,
  data = dataDF
)

# Plot as percentage survival (i.e. % still being followed)
p_rev <- ggsurvplot(
  rev_fit_all,
  fun             = "pct",        # multiply S(t) by 100
  conf.int        = TRUE,         # 95% CI bands
  risk.table      = TRUE,         # show # at risk = # still under follow-up
  risk.table.title= "Still under follow-up",
  xlab            = "Months since start",
  ylab            = "Follow-up probability (%)",
  title           = "Reverse Kaplan–Meier: Percent under follow-up",
  legend.title    = "Tumour",
  palette         = scales::hue_pal()(length(levels(dataDF$Tumor))),
  ggtheme         = theme_classic(base_size = 14)
)

# Tweak the y-axis so it shows percentages
p_rev$plot <- p_rev$plot +
  scale_y_continuous(
    limits = c(0,100),
    breaks = seq(0,100, by = 20),
    labels = function(x) paste0(x, "%")
  )

p_rev<- p_rev$plot+geom_hline(yintercept = 50, linetype = "dashed")

print(p_rev)
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/follow%20up%20calculation-1.png)<!-- -->

# Percentages of administrative censored patients

## Administrative censoring rates per tumour

``` r
library(dplyr)

admin_censoring_OS <- dataDF %>%
  dplyr::summarise(
    total          = n(),
    admin_censored = sum(OS_months > 60),
    pct_censored   = admin_censored / total * 100
  )

admin_censoring_OS %>%
  mutate(
    summary = sprintf(
      "%d/%d (%.1f%%) administratively censored at 60 months (OS)",
      admin_censored, total, pct_censored
    )
  ) %>%
  pull(summary) %>%
  cat(sep = "\n")
```

    ## 645/2117 (30.5%) administratively censored at 60 months (OS)

``` r
#PFS
admin_censoring_PFS <- dataDF %>%
  dplyr::summarise(
    total          = n(),
    admin_censored = sum(PFS_months > 60),
    pct_censored   = admin_censored / total * 100
  )

admin_censoring_PFS %>%
  mutate(
    summary = sprintf(
      "%d/%d (%.1f%%) administratively censored at 60 months (PFS)",
      admin_censored, total, pct_censored
    )
  ) %>%
  pull(summary) %>%
  cat(sep = "\n")
```

    ## 374/2117 (17.7%) administratively censored at 60 months (PFS)

# Differential distribution of CRs

``` r
# Contingency table: Tumour × Response (CR vs PR)
resp_tab <- table(dataDF$Tumor, dataDF$DOR)
print(resp_tab)
```

    ##           
    ##             CR  PR
    ##   Melanoma 576 620
    ##   NSCLC     73 591
    ##   RCC       79 178

``` r
# χ² test of independence
chi2_res <- chisq.test(resp_tab)
print(chi2_res)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  resp_tab
    ## X-squared = 263.12, df = 2, p-value < 2.2e-16

``` r
# Extract and display CR proportions by tumour
cr_props <- prop.table(resp_tab, 1)[, "CR"] * 100
cat("CR proportion by tumour (%):\n")
```

    ## CR proportion by tumour (%):

``` r
print(round(cr_props, 1))
```

    ## Melanoma    NSCLC      RCC 
    ##     48.2     11.0     30.7

``` r
# Build a matrix of counts: rows = Tumour, cols = Response
resp_tab <- table(dataDF$Tumor, dataDF$DOR)

###checking here with fisher's exact test as well
# Pairwise tests of proportions, Bonferroni correction
library(rstatix)

tumours <- rownames(resp_tab)
pairs   <- combn(tumours, 2, simplify = FALSE)

fisher_results <- list()
for (i in 1:length(pairs)){
  pair <- pairs[[i]]
  sub_tab <- resp_tab[pair, , drop = FALSE]
  pval    <- fisher.test(sub_tab)$p.value
  fisher_results[[i]] <- data.frame(
    comparison = paste(unlist(pair), collapse = " vs "),
    p.value    = pval,
    stringsAsFactors = FALSE
    )
}


fisher_df <- do.call(rbind, fisher_results) %>%
  mutate(p.adj = p.adjust(p.value, method = "bonferroni"))

cat("\nPairwise Fisher’s Exact (with Bonferroni adj):\n")
```

    ## 
    ## Pairwise Fisher’s Exact (with Bonferroni adj):

``` r
print(fisher_df)
```

    ##          comparison      p.value        p.adj
    ## 1 Melanoma vs NSCLC 1.487667e-64 4.463001e-64
    ## 2   Melanoma vs RCC 2.824761e-07 8.474284e-07
    ## 3      NSCLC vs RCC 5.161822e-12 1.548547e-11

\#Median time to event

``` r
responses <- c("PR","CR")
tumors   <- c("Melanoma","RCC","NSCLC")
abbr     <- c("Melanoma"="MM", "NSCLC"="NSCLC", "RCC"="RCC")

dataDF$censor_time_PFS   <- pmin(dataDF$PFS_months, 60)
dataDF$censor_status_PFS <- ifelse(dataDF$PFS_months > 60, 0, dataDF$Progressed)

# Loop over responses (PR and CR)
for (d in responses) {
  medians <- c()
  for (t in tumors) {
    df_sub <- dplyr::filter(dataDF, DOR == d, Tumor == t, censor_status_PFS == 1)
    if (nrow(df_sub) > 0) {
      fit <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1, data = df_sub)
      median_time <- summary(fit)$table["median"]
      medians <- c(medians, median_time)  # Store median time for each tumor type
    } else {
      # If no progression, add NA
      medians <- c(medians, NA)
    }
  }
  # Format each median with the tumor abbreviation
  parts <- sprintf("%.2f months (%s)", medians, abbr[tumors])
  
  # Combine the parts into a single string with commas and "and"
  txt <- paste0(
    paste(parts[1:2], collapse = ", "),
    " and ",
    parts[3]
  )
  
  # Print the final result
  cat(sprintf(
    "Median time-to-progression among %s was %s.\n",
    d, txt
  ))
}
```

    ## Median time-to-progression among PR was 10.87 months (MM), 13.70 months (RCC) and 15.18 months (NSCLC).
    ## Median time-to-progression among CR was 23.95 months (MM), 31.27 months (RCC) and 31.95 months (NSCLC).

# Post-hoc power analysis

``` r
library(survival)
library(stats)

# Schoenfeld-based functions
schoenfeld_power <- function(D, p = 0.5, HR, alpha = 0.05) {
  if (HR <= 0) stop("HR must be > 0")
  z_alpha <- qnorm(1 - alpha / 2)
  delta <- sqrt(D * p * (1 - p)) * abs(log(HR))
  power <- pnorm(delta - z_alpha)
  return(power)
}

schoenfeld_required_events <- function(target_power = 0.8, p = 0.5, HR, alpha = 0.05) {
  if (HR <= 0) stop("HR must be > 0")
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta  <- qnorm(target_power)                # 1 - beta quantile
  D <- ((z_alpha + z_beta)^2) / (p * (1 - p) * (log(HR))^2)
  return(ceiling(D))
}

# Helper to extract info from a matched_df (one comparison)
posthoc_power_from_matched <- function(matched_df, group_var = "Tumor",
                                       timevar = "PFS_months", statusvar = "Progressed",
                                       alpha = 0.05) {
  # make sure group is a factor and compute allocation p
  matched_df[[group_var]] <- droplevels(as.factor(matched_df[[group_var]]))
  tbl <- table(matched_df[[group_var]])
  # pick first level as "reference" for p
  levels_group <- levels(matched_df[[group_var]])
  if (length(levels_group) != 2) stop("Need exactly 2 groups for pairwise power.")
  n1 <- as.integer(tbl[1]); n2 <- as.integer(tbl[2])
  p <- n1 / (n1 + n2)

  # events
  D <- sum(matched_df[[statusvar]] == 1, na.rm = TRUE)

  # fit Cox (naive)
  f <- as.formula(paste0("Surv(", timevar, ", ", statusvar, ") ~ ", group_var))
  fit <- coxph(f, data = matched_df)
  s <- summary(fit)
  logHR <- s$coefficients[1, "coef"]
  se_logHR <- s$coefficients[1, "se(coef)"]
  HR <- exp(logHR)

  # observed power via Schoenfeld:
  power_obs <- schoenfeld_power(D = D, p = p, HR = HR, alpha = alpha)
  # required events for 80% and 90%:
  req80 <- schoenfeld_required_events(0.8, p = p, HR = HR, alpha = alpha)
  req90 <- schoenfeld_required_events(0.9, p = p, HR = HR, alpha = alpha)

  list(
    n1 = n1, n2 = n2, D = D, p = p,
    HR = HR, logHR = logHR, se_logHR = se_logHR,
    observed_power = power_obs,
    required_events_80 = req80,
    required_events_90 = req90,
    cox_summary = s
  )
}
```

``` r
design_info_twoarm <- function(df,
                               group_var  = "Tumor",
                               timevar    = "censor_time_PFS",
                               statusvar  = "censor_status_PFS") {
  vars_needed <- c(group_var, timevar, statusvar)
  df <- df[complete.cases(df[, vars_needed]), ]

  g <- droplevels(as.factor(df[[group_var]]))
  if (nlevels(g) != 2L)
    stop("Need exactly 2 levels in ", group_var, " (got: ",
         paste(levels(g), collapse = ", "), ")")

  tab <- table(g)
  n1  <- as.integer(tab[1])
  n2  <- as.integer(tab[2])
  p   <- n1 / (n1 + n2)

  D <- sum(df[[statusvar]] == 1L, na.rm = TRUE)

  list(n1 = n1, n2 = n2, p = p, D = D)
}

power_for_pairs_target_HR <- function(df,
                                      comparisons,
                                      response_level,
                                      HR_target,
                                      timevar   = "censor_time_PFS",
                                      statusvar = "censor_status_PFS",
                                      alpha     = 0.05) {

  out_list <- list()

  # Restrict to chosen response level
  if(!is.na(response_level)){
    df <- subset(df, DOR == response_level)
  } else {
    df=df
  }

  for (pair in comparisons) {
    pair <- as.character(pair)

    sub_df <- df[df$Tumor %in% pair, , drop = FALSE]
    # fix ordering (important for n1 / n2 labelling, but not for power)
    sub_df$Tumor <- droplevels(factor(sub_df$Tumor, levels = pair))

    info <- design_info_twoarm(sub_df,
                               group_var = "Tumor",
                               timevar   = timevar,
                               statusvar = statusvar)

    power <- schoenfeld_power(D = info$D,
                              p = info$p,
                              HR = HR_target,
                              alpha = alpha)

    out_list[[paste(pair, collapse = " vs ")]] <- data.frame(
      comparison = paste(pair, collapse = " vs "),
      response   = response_level,
      HR_target  = HR_target,
      n1         = info$n1,
      n2         = info$n2,
      events_D   = info$D,
      power      = power,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

power_pairs_CR_ls <- list()
for (i in seq_along(comparisons)){
  data <- subset(dataDF, Tumor %in% unlist(comparisons[i]))
  fit_CR <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1,
                      data = data)
  m_ref_CR   <- as.numeric(summary(fit_CR)$table["median"])
  delta_CR   <- 12
  HR_target_CR <- m_ref_CR / (m_ref_CR + delta_CR)

  power_pairs_CR <- power_for_pairs_target_HR(
  df           = data,
  comparisons  = comparisons[i],
  response_level = "CR",
  HR_target    = HR_target_CR,
  timevar      = "censor_time_PFS",
  statusvar    = "censor_status_PFS",
  alpha        = 0.05)
  
  power_pairs_CR_ls[[i]] <- power_pairs_CR 
}

power_CR <- do.call(rbind, power_pairs_CR_ls)
power_CR$group <- "Unmatched CR"

power_pairs_PR_ls <- list()
for (i in seq_along(comparisons)){
  data <- subset(dataDF, Tumor %in% unlist(comparisons[i]))
  fit_PR <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1,
                      data = data)
  m_ref_PR   <- as.numeric(summary(fit_PR)$table["median"])
  delta_PR   <- 12
  HR_target_PR <- m_ref_PR / (m_ref_PR + delta_PR)

  power_pairs_PR <- power_for_pairs_target_HR(
  df           = data,
  comparisons  = comparisons[i],
  response_level = "PR",
  HR_target    = HR_target_PR,
  timevar      = "censor_time_PFS",
  statusvar    = "censor_status_PFS",
  alpha        = 0.05)
  
  power_pairs_PR_ls[[i]] <- power_pairs_PR 
}

power_PR <- do.call(rbind, power_pairs_PR_ls)
power_PR$group <- "Unmatched PR"

power_pairs_psm_ls <- list()
for (i in seq_along(comparisons)){
  matched_data <- PFS_psm_all$matched_dfs[[i]]
  fit_psm <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1,
                      data = matched_data)
  m_ref_psm   <- as.numeric(summary(fit_psm)$table["median"])
  delta_psm   <- 12
  HR_target_psm <- m_ref_psm / (m_ref_psm + delta_psm)

  power_pairs_psm <- power_for_pairs_target_HR(
  df           = matched_data,
  comparisons  = comparisons[i],
  response_level = NA,
  HR_target    = HR_target_psm,
  timevar      = "censor_time_PFS",
  statusvar    = "censor_status_PFS",
  alpha        = 0.05)
  
  power_pairs_psm_ls[[i]] <- power_pairs_psm 
}

power_psm <- do.call(rbind, power_pairs_psm_ls)
power_psm$group <- "PSM"

power_table_final <- rbind(power_CR,power_PR, power_psm)
power_table_final
```

    ##                           comparison response HR_target  n1  n2 events_D
    ## Melanoma vs RCC      Melanoma vs RCC       CR 0.7771705 576  79      147
    ## Melanoma vs NSCLC  Melanoma vs NSCLC       CR 0.7201520 576  73      150
    ## RCC vs NSCLC            RCC vs NSCLC       CR 0.6736473  79  73       45
    ## Melanoma vs RCC1     Melanoma vs RCC       PR 0.7771705 620 178      533
    ## Melanoma vs NSCLC1 Melanoma vs NSCLC       PR 0.7201520 620 591      862
    ## RCC vs NSCLC1           RCC vs NSCLC       PR 0.6736473 178 591      591
    ## Melanoma vs RCC2     Melanoma vs RCC     <NA> 0.7222797 256 256      276
    ## Melanoma vs NSCLC2 Melanoma vs NSCLC     <NA> 0.6661915 472 472      606
    ## RCC vs NSCLC2           RCC vs NSCLC     <NA> 0.6813344 219 219      303
    ##                        power        group
    ## Melanoma vs RCC    0.1673872 Unmatched CR
    ## Melanoma vs NSCLC  0.2452294 Unmatched CR
    ## RCC vs NSCLC       0.2624001 Unmatched CR
    ## Melanoma vs RCC1   0.6782874 Unmatched PR
    ## Melanoma vs NSCLC1 0.9978682 Unmatched PR
    ## RCC vs NSCLC1      0.9817206 Unmatched PR
    ## Melanoma vs RCC2   0.7711189          PSM
    ## Melanoma vs NSCLC2 0.9988151          PSM
    ## RCC vs NSCLC2      0.9161399          PSM

``` r
# kable(
#   power_table_final,
#   format    = "html",
#   digits    = 2,              # round numeric columns to 2 digits
#   caption   = "Power Calculations",row.names = F
# ) %>%
#   kable_styling(
#     bootstrap_options = c("striped", "hover", "condensed"), 
#     full_width        = FALSE
#   )
```

``` r
# Choose the deltas you want to explore (months)
deltas <- seq(3,42, by=3)

power_all_delta_ls <- list()

for (delta in deltas) {

  ## ---------- CR: unmatched ----------
  power_pairs_CR_ls <- list()
  for (i in seq_along(comparisons)) {
    data <- subset(dataDF, Tumor %in% unlist(comparisons[i]))
    
    fit_CR <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1,
                      data = data)
    m_ref_CR     <- as.numeric(summary(fit_CR)$table["median"])
    delta_CR     <- delta
    HR_target_CR <- m_ref_CR / (m_ref_CR + delta_CR)

    power_pairs_CR <- power_for_pairs_target_HR(
      df            = data,
      comparisons   = comparisons[i],
      response_level = "CR",
      HR_target     = HR_target_CR,
      timevar       = "censor_time_PFS",
      statusvar     = "censor_status_PFS",
      alpha         = 0.05
    )
    
    power_pairs_CR_ls[[i]] <- power_pairs_CR 
  }
  power_CR <- do.call(rbind, power_pairs_CR_ls)
  power_CR$group <- "Unmatched CR"
  power_CR$delta <- delta
  
  ## ---------- PR: unmatched ----------
  power_pairs_PR_ls <- list()
  for (i in seq_along(comparisons)) {
    data <- subset(dataDF, Tumor %in% unlist(comparisons[i]))
    
    fit_PR <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1,
                      data = data)
    m_ref_PR     <- as.numeric(summary(fit_PR)$table["median"])
    delta_PR     <- delta
    HR_target_PR <- m_ref_PR / (m_ref_PR + delta_PR)

    power_pairs_PR <- power_for_pairs_target_HR(
      df            = data,
      comparisons   = comparisons[i],
      response_level = "PR",
      HR_target     = HR_target_PR,
      timevar       = "censor_time_PFS",
      statusvar     = "censor_status_PFS",
      alpha         = 0.05
    )
    
    power_pairs_PR_ls[[i]] <- power_pairs_PR 
  }
  power_PR <- do.call(rbind, power_pairs_PR_ls)
  power_PR$group <- "Unmatched PR"
  power_PR$delta <- delta
  
  ## ---------- PSM: all patients (no DOR restriction) ----------
  power_pairs_psm_ls <- list()
  for (i in seq_along(comparisons)) {
    matched_data <- PFS_psm_all$matched_dfs[[i]]
    
    fit_psm <- survfit(Surv(censor_time_PFS, censor_status_PFS) ~ 1,
                       data = matched_data)
    m_ref_psm     <- as.numeric(summary(fit_psm)$table["median"])
    delta_psm     <- delta
    HR_target_psm <- m_ref_psm / (m_ref_psm + delta_psm)

    power_pairs_psm <- power_for_pairs_target_HR(
      df            = matched_data,
      comparisons   = comparisons[i],
      response_level = NA,   # as in your original code
      HR_target     = HR_target_psm,
      timevar       = "censor_time_PFS",
      statusvar     = "censor_status_PFS",
      alpha         = 0.05
    )
    
    power_pairs_psm_ls[[i]] <- power_pairs_psm 
  }
  power_psm <- do.call(rbind, power_pairs_psm_ls)
  power_psm$group <- "PSM"
  power_psm$delta <- delta
  
  ## ---------- bind for this delta ----------
  power_all_delta_ls[[as.character(delta)]] <- 
    rbind(power_CR, power_PR, power_psm)
}

# Final long table with all deltas
power_table_heat <- do.call(rbind, power_all_delta_ls)

# Make sure things are nicely formatted for plotting
power_table_heat <- power_table_heat %>%
  mutate(
    delta      = as.factor(delta),
    group      = factor(group, levels = c("Unmatched PR", "Unmatched CR", "PSM")),
    comparison = as.factor(comparison)
  )

power_table_heat_plot <- power_table_heat %>%
  mutate(
    power_trunc = pmax(power, 0.7) 
  )

ggplot(power_table_heat_plot,
       aes(x = delta, y = comparison, fill = power_trunc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(
    name   = "Power",
    breaks = c(0.7, 0.8, 0.9, 1.0),
    labels = c("<0.7", "0.8", "0.9", "1.0"),
    option = "D"
  ) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(
    title = "Post-hoc power to detect Δ median PFS",
    x     = "Δ PFS (months)",
    y     = "Tumour comparison"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    strip.text      = element_text(size = 14, face = "bold"),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    legend.position = "right"
  )
```

![](PropensityScore_joined_cohorts_mel_rcc_lung_files/figure-gfm/power%20heatmap-1.png)<!-- -->
