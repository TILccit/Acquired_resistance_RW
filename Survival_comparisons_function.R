run_survival_comparisons <- function(dataDF,
                                     comparisons,
                                     covariate_to_split,
                                     descr_formula   = NULL,
                                     endpoint        = c("OS", "PFS"),
                                     truncate_month  = NULL) {
  endpoint <- match.arg(endpoint)
  stopifnot(require(survival), require(ggsurvfit),
            require(ggplot2),    require(dplyr),
            require(ggpubr))
  if (!is.null(descr_formula)) stopifnot(require(compareGroups))
  
  # prepare names & palette
  comp_names <- vapply(
    comparisons,
    function(p) paste(p[1], "vs", p[2]),
    character(1)
  )
  all_lvls <- levels(factor(dataDF[[covariate_to_split]]))
  pal      <- scales::hue_pal()(length(all_lvls))
  names(pal) <- all_lvls
  
  # storage
  surv_plots         <- vector("list", length(comparisons))
  descriptive_tables <- if (!is.null(descr_formula)) vector("list", length(comparisons)) else NULL
  
  # loop over each pair
  for (i in seq_along(comparisons)) {
    pair   <- comparisons[[i]]
    df_sub <- subset(dataDF, dataDF[[covariate_to_split]] %in% pair)
    df_sub[[covariate_to_split]] <- droplevels(df_sub[[covariate_to_split]])
    
    # optional descriptives on the unmatched subset
    if (!is.null(descr_formula)) {
      descriptive_tables[[i]] <- 
        compareGroups::descrTable(descr_formula, data = df_sub)
    }
    
    # pick OS vs PFS
    if (endpoint == "OS") {
      timevar    <- "OS_months"; statusvar <- "Dead"
      main_title <- "Overall Survival"
    } else {
      timevar    <- "PFS_months"; statusvar <- "Progressed"
      main_title <- "Progression-Free Survival"
    }
    
    # apply administrative censoring if requested
    if (!is.null(truncate_month)) {
      orig_t <- df_sub[["OS_months"]]
      df_sub[["OS_months"]]   <- pmin(orig_t, truncate_month)
      df_sub[["Dead"]] <- ifelse(orig_t > truncate_month, 0, df_sub[["Dead"]])
    }

    # fit KM + Cox
    surv_formula <- as.formula(paste0("Surv(", timevar, ", ", statusvar, ") ~ ",
                                      covariate_to_split))
    fit_km  <- surv_fit(surv_formula, data = df_sub)
    fit_cox <- coxph(surv_formula, data = df_sub)
    s_cox   <- summary(fit_cox)
    
    # extract HR, CI and p
    hr     <- s_cox$coefficients[, "exp(coef)"]
    ci_lo  <- s_cox$conf.int[, "lower .95"]
    ci_hi  <- s_cox$conf.int[, "upper .95"]
    pval   <- s_cox$coefficients[, "Pr(>|z|)"]
    hr_lab <- sprintf("HR=%.2f (95%% CI %.2f–%.2f)", hr, ci_lo, ci_hi)
    p_lab  <- paste0("p = ", sub("^0\\.", ".", sprintf("%.4f", pval)))
    
    # clean up strata labels
    names(fit_km$strata) <- sub(
      paste0(covariate_to_split, "="), "", names(fit_km$strata)
    )
    
    # determine x-axis maximum
    x_max <- if (is.null(truncate_month)) {
      max(df_sub[[timevar]], na.rm = TRUE)
    } else truncate_month
    
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
      add_risktable_strata_symbol(symbol = "•", size = 20) +
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
  names(surv_plots) <- comp_names
  if (!is.null(descriptive_tables)) names(descriptive_tables) <- comp_names
  
  # assemble combined plot (2 per row)
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
  
  # return
  list(
    surv_plots         = surv_plots,
    descriptive_tables = descriptive_tables,
    combined_surv_plot = combined_surv
  )
}
