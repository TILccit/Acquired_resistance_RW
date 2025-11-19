#’ Run propensity score matching & plot survival for pairwise multiple cohort comparisons
#’
#’ @description
#’ Takes a data frame with a cohort‐defining column and runs
#’ nearest‐neighbor propensity score matching (caliper = 0.2, standardized)
#’ for each pair you supply.  It returns:
#’ 1) balance tables, 2) love plots, 3) Kaplan–Meier curves (w/ HR & 95% CI),
#’ and—if you supply one—4) descriptive tables via `descrTable()`.  
#’ Finally, it assembles all the PSM survival plots into one combined grid
#’ (with the third plot centered if there are exactly three).
#’
#’ @param dataDF A `data.frame` containing at least:
#’   * your cohort column (e.g. `Cohort` or `Tumor`)  
#’   * `OS_months`, `Dead`            — for overall survival  
#’   * `PFS_months`, `Progressed`     — for progression‐free survival  
#’   * plus whatever covariates your `matching_formula` uses  
#’
#’ @param comparisons A `list` of character vectors (length 2).  Each element
#’   names the two cohorts to compare, e.g.
#’   `list(c("Melanoma","Renal cancer"), c("Lung cancer","Melanoma"))`.  
#’
#’ @param covariate_to_split A single variable name (string) or simple RHS
#’   formula (e.g. `"Tumor"`, `"Treatment"`) used to define the strata in
#’   the KM curves and Cox model.  
#’
#’ @param descr_formula Optional survival‐independent formula (e.g.  
#’   `~ Age + Sex + Stage`).  If non‐NULL, for each matched cohort the
#’   function will call  
#’   `descrTable(descr_formula, data = matched_data_i)` and return those
#’   tables in `descriptive_tables`.  
#’
#’ @param endpoint Character, either `"OS"` (overall survival) or `"PFS"`
#’   (progression‐free survival).  Default: `"OS"`.  
#’
#’ @param truncate_month Optional numeric.  If non‐NULL, all times > this
#’   value are set to `truncate_month` and status to 0 (administrative
#’   censoring), and both curves & HRs are computed on the truncated data.
#’   Default: `NULL` (no truncation).  
#’
#’ @return A named list with five elements:
#’   * `balance_tables`     : named list of before/after matching tables  
#’   * `love_plots`         : named list of `ggplot` love plots  
#’   * `surv_plots`         : named list of individual `ggsurvfit` survival plots (w/ HR & CI)  
#’   * `descriptive_tables` : named list of `descrTable()` outputs (or `NULL`)  
#’   * `combined_surv_plot` : a single `ggpubr::ggarrange` object assembling all survival plots  
#’
#’ @examples
#’ \dontrun{
#’ matching_formula <- Cohort ~ Age + Sex + Stage
#’ comps <- list(c("Melanoma","Renal cancer"), c("Lung cancer","Melanoma"))
#’ descr_f  <- ~ Age + Sex + Stage
#’ out1 <- run_psm_and_survival(
#’   myData, comps, "Cohort",
#’   descr_formula  = descr_f,
#’   endpoint       = "OS",
#’   truncate_month = 24
#’ )
#’ # names:
#’ names(out1$balance_tables)
#’ # individual plot:
#’ print(out1$surv_plots[[1]])
#’ # arranged:
#’ print(out1$combined_surv_plot)
#’ }
#’ @export
run_psm_and_survival <- function(dataDF,
                                 comparisons,
                                 covariate_to_split,
                                 matching_formula = NULL,
                                 descr_formula    = NULL,
                                 mahvars=NULL,
                                 endpoint         = c("OS", "PFS"),
                                 truncate_month   = NULL,
                                 caliper = 0.2,
                                 ratio = 1) {
  endpoint <- match.arg(endpoint)
  stopifnot(require(MatchIt), require(cobalt),
            require(survival), require(ggsurvfit),
            require(ggplot2), require(dplyr),
            require(ggpubr))   # for ggarrange()
  
  ncmp <- length(comparisons)
  comp_names <- vapply(
    comparisons,
    function(p) paste(p[1], "vs", p[2]),
    character(1)
  )
  
  #generate a color palette for the covariate to split
  all_levels <- levels(factor(dataDF[[covariate_to_split]]))
  pal <- scales::hue_pal()(length(all_levels))
  names(pal) <- all_levels
  
  # prepare storage
  balance_tbls       <- vector("list", ncmp)
  love_plots         <- vector("list", ncmp)
  surv_plots         <- vector("list", ncmp)
  descriptive_tables <- if (!is.null(descr_formula)) vector("list", ncmp) else NULL
  matched_dfs        <- vector("list", ncmp)
  
  # PSM & love plots (and descriptives)
  for (i in seq_along(comparisons)) {
    pair   <- comparisons[[i]]
    df_sub <- subset(dataDF, dataDF[[covariate_to_split]] %in% pair)
    df_sub[[covariate_to_split]] <- droplevels(df_sub[[covariate_to_split]])
    
    fit_i <- matchit(
      matching_formula,
      data        = df_sub,
      method      = "nearest",
      distance    = "glm",
      link="linear.logit",
      ratio       = ratio,
      caliper     = 0.2,
      std.caliper = TRUE,
      mahvars = mahvars
    )
    
    balance_tbls[[i]] <- bal.tab(fit_i, un = TRUE)
    orig <- rownames(balance_tbls[[i]]$Balance)
    #clean labels
    new <- sub("_.*",    "", orig)      # drop “:Level”
    if("DOR" %in% new){ #this bit is slightly hardcoded but it was the quickest fix
      new <- sub("DOR", "Objective Response", new)
    }
    # stick them into a named vector
    lab_map <- setNames(new, orig)
    
    love_plots[[i]]   <- love.plot(fit_i, abs = TRUE, var.names = lab_map) +
      ggtitle(sprintf("SMD Before & After Matching: %s", comp_names[i])) +
      geom_vline(xintercept = 0.2, linetype = "dashed", size = 1)
    
    matched_dfs[[i]] <- match.data(fit_i)
    
    if (!is.null(descr_formula)) {
      descriptive_tables[[i]] <-
        descrTable(descr_formula, data = matched_dfs[[i]])
    }
  }
  
  # Survival plots w/ HR & optional truncation
  for (i in seq_along(matched_dfs)) {
    md   <- matched_dfs[[i]]
    
    if (endpoint == "OS") {
      timevar    <- "OS_months"; statusvar <- "Dead"
      main_title <- "Propensity-Matched Overall Survival"
    } else {
      timevar    <- "PFS_months"; statusvar <- "Progressed"
      main_title <- "Propensity-Matched PFS"
    }
    
    # administrative censoring?
    if (!is.null(truncate_month)) {
      orig_t <- md[["OS_months"]]
      md[["OS_months"]]   <- pmin(orig_t, truncate_month)
      md[["Dead"]] <- ifelse(orig_t > truncate_month, 0, md[["Dead"]])
    }
    
    surv_formula <- as.formula(
      paste0("Surv(", timevar, ", ", statusvar, ") ~ ", covariate_to_split)
    )
    
    fit_km  <- survfit(surv_formula, data = md)
    fit_cox <- coxph(surv_formula, data = md)
    s_cox   <- summary(fit_cox)
    hr      <- s_cox$coefficients[, "exp(coef)"]
    ci_lo   <- s_cox$conf.int[,    "lower .95"]
    ci_hi   <- s_cox$conf.int[,    "upper .95"]
    hr_lab  <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)",
                       hr[1], ci_lo[1], ci_hi[1])
    pval   <- s_cox$coefficients[, "Pr(>|z|)"]
    p_lab  <- paste0("p = ",sub("^0\\.", ".", sprintf("%.4f", pval)))
    
    names(fit_km$strata) <- sub(
      paste0(covariate_to_split, "="), "", names(fit_km$strata)
    )
    
    x_max <- if (is.null(truncate_month))
      max(md[[timevar]], na.rm = TRUE) else truncate_month
    
    ################
    surv_plots[[i]] <- ggsurvfit(fit_km, size = 1.5) +
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
        y        = paste0(endpoint, "(%)"),
        subtitle = comp_names[i]
      ) +
      scale_x_continuous(breaks = seq(0, x_max, by = 12), limits = c(0, x_max)) +
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
      scale_color_manual(name   = covariate_to_split, values = pal) +
      scale_fill_manual(name    = covariate_to_split, values = pal) +
      annotate(
        "text",
        x     = x_max * 0.1,
        y     = 0.2,
        label = hr_lab,
        size  = 10,
        hjust = 0
      ) +
      annotate(
        "text",
        x     = x_max * 0.1,
        y     = 0.1,
        label = p_lab,
        size  = 10,
        hjust = 0
      )
    surv_plots[[i]] <- ggsurvfit_build(surv_plots[[i]], combine_plots = TRUE)
  }
  
  
  # name each list for easy extraction
  names(balance_tbls)       <- comp_names
  names(love_plots)         <- comp_names
  names(surv_plots)         <- comp_names
  if (!is.null(descriptive_tables))
    names(descriptive_tables) <- comp_names
  
  # assemble into one combined plot
  combined_surv <- ggarrange(
    plotlist       = surv_plots,
    ncol           = length(surv_plots),
    nrow           = 1,
    common.legend  = F,
    legend         = "top",
    align = "v"
  ) +theme(
    plot.margin   = margin(0,0,0,0, "pt"),   
    panel.spacing = unit(0, "pt") 
  )
  
  endpoint_long <- ifelse(endpoint == "PFS", "Progression-Free Survival",
                          ifelse(endpoint == "OS", "Overall Survival", "NA"))
  
  combined_surv <- combined_surv + 
    plot_annotation(
      title    = paste0(endpoint_long," after Propensity Score Matching"),
      theme    = theme(
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold")))
        
  
  list(
    balance_tables      = balance_tbls,
    love_plots          = love_plots,
    surv_plots          = surv_plots,
    descriptive_tables  = descriptive_tables,
    combined_surv_plot  = combined_surv,
    matched_dfs         = matched_dfs
  )
}
