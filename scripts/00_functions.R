#############################################################################
## Functions for the ASIS project:
#############################################################################

#############################################################################
# 1. Functions for nice tables in the reports:
#############################################################################

## A) Table 1
bold <- function(x){paste0('{\\textbf{', x, '}}')}

## table1 as data frame
tab_df <- function(vars, vars_cat, vars_nonnorm, strat, dat,
                   caption, label, filename,
                   sideways = "table", fontsize = "footnotesize", SMD=TRUE,
                   overall = TRUE, p = FALSE, Levels=TRUE, Missing = TRUE) {
  tab <- CreateTableOne(
    vars = vars,
    factorVars = vars_cat,
    data = dat,
    strata = strat,
    addOverall = overall,
    includeNA = TRUE
  )
  tab <- print(
    tab,
    nonnormal = vars_nonnorm,
    smd = SMD, 
    contDigits = 1,
    showAllLevels = Levels,
    missing = Missing
  )
  
  ## sanitation of row names
  # get quad-command in it for indentation of factor levels
  rownames(tab) <- str_replace_all(string = rownames(tab), pattern = "   ",
                                   replacement = "\\\\quad ")
  # handle %-symbol in tex
  rownames(tab) <- str_replace_all(string = rownames(tab), pattern = "%", 
                                   replacement = "\\\\%")
  # handle _-symbol in tex
  rownames(tab) <- str_replace_all(string = rownames(tab), pattern = "_", 
                                   replacement = "\\\\_")
  # replace NA in tex
  rownames(tab) <- str_replace_all(string = rownames(tab), pattern = "NA", 
                                   replacement = "Missing")
  tab[,1] <- str_replace_na(string = tab[,1], replacement = "Missing")
  
  colnames(tab)[which(colnames(tab) == "Missing")] <- "Missing (\\%)"
  colnames(tab)[which(colnames(tab) == "p")] <- "p-val"
  colnames(tab)[which(colnames(tab) == "level")] <- "Levels"
  if(p == FALSE) {ind <- which(colnames(tab) %in% c("p", "test"))}
  if(p == TRUE) {ind <- which(colnames(tab) == "test")}
  tab.df <- cbind(Variable = rownames(tab), tab[, -ind])
  
  ## table1 as latex table
  xtab.tab <- xtable(tab.df, caption = caption, label = label, 
                     align = paste0("ll|", strrep("r", ncol(tab.df)-1)))
  hlines <- c(-1, 0, 1, nrow(xtab.tab))
  
  return(
    print(xtab.tab, 
          include.rownames = FALSE, 
          sanitize.text.function = identity,
          sanitize.colnames.function = bold,
          hline.after = hlines,
          # table.placement = "H",
          caption.placement = "top",
          floating = TRUE, 
          type = "latex",
          floating.environment = sideways,
          size = paste0("\\", fontsize),
          file = filename))}


### B) Summary tables for the regression models:
SummaryTable <- function(fit, n, digits=2){
  if ("glmmTMB" %in% class(fit)) {
    coefs <- summary(fit)$coefficients$cond
    logodds <- coefs[-c(-1:-n), 1]
    estimates <- exp(logodds)
    se <- coefs[-c(-1:-n), 2]
    df <- df.residual(fit)
    z_crit <- qnorm(0.975)
    CI_low <- round(exp(logodds - z_crit * se), digits)
    CI_upper <- round(exp(logodds + z_crit * se), digits)
    pval <- coefs[-c(-1:-n), "Pr(>|z|)"]
  } else {
    coefs <- summary(fit)$coefficients
    df_res <- fit$df.residual
    logodds <- coefs[-c(-1:-n), 1]
    estimates <- exp(logodds)
    se <- coefs[-c(-1:-n), 2]
    CI_low <- round(exp(logodds - qnorm(0.975) * se), digits)
    CI_upper <- round(exp(logodds + qnorm(0.975) * se), digits)
    tval <- coefs[-c(-1:-n), 3]
    pval <- 2 * pt(abs(tval), df = df_res, lower.tail = FALSE)
  }
  
  pval <- biostatUZH::formatPval(pval)
  
  mat <- cbind("log(OR)" = round(logodds, digits),
               "OR" = round(estimates, digits),
               "95%-Confidence Interval" = paste(CI_low, "to", CI_upper),
               "p-value" = pval)
  return(mat)}



SummaryTableclmm <- function(fit,n, t=6){
  logodds <- summary(fit)$coefficients[c(seq(t+1,t+n,1)),1]
  estimates <- exp(logodds)
  se <- summary(fit)$coefficients[c(seq(t+1,t+n,1)),2]
  CI_low <- round(exp(logodds-qnorm(0.975)*se),2)
  CI_upper <- round(exp(logodds+qnorm(0.975)*se),2)
  pval <- summary(fit)$coefficients[c(seq(t+1,t+n,1)),4]
  pval <- formatPval(pval)
  mat= cbind("log(OR)" = round(logodds,2),
             "OR" = round(estimates,2),
             "95%-Confidence Interval" =
               paste(CI_low, "to ", CI_upper),
             "p-value" = pval)
  return(mat)}


#############################################################################
# 2. Functions for the forest plots of the regression models:
#############################################################################
### Plot the MICE Forest plots for the regression models:
ForestPlot_MICE <- function(model, title=NULL, col_def="darkred"){
  sum_df <- summary(model, conf.int = TRUE, exponentiate = TRUE,
                    conf.level = 0.95, digits = 3) %>% 
    select(-c(`2.5 %`, `97.5 %`, statistic, df))
  
  fmi_df <- model$pooled[, c("term", "fmi")]
  
  res_df <- left_join(sum_df, fmi_df, by = "term")
  
  res_for_plot <- res_df %>%
    filter(term != "(Intercept)") %>%  
    mutate(term = factor(term, levels = rev(term)))  %>% 
    mutate(term= recode(term,"genderMale"= "Male", 
                        "Positive_familial_historyNo"="No Family History" ,
                        "Positive_familial_historyYes" = "Family History",
                        "prf_hypertension_awarenessYes"="Aware of Hypertension" ,
                        "prf_smokedtobaccoFormer" = "Former Smoker",
                        "prf_smokedtobaccoCurrent" = "Current Smoker",
                        "aneuriskMedium" = "Medium Risk Location",
                        "aneuriskHigh" = "High Risk Location",
                        "scale(maxdiam)" = "IA Size, scaled, log",
                        "typeSaccular - Side-Wall"= "Side-Wall",
                        "sideMidline" = "IA at Mideline",
                        "sideRight" = "IA Right-sided" ,
                        "timepointDischarge"= "Discharge" ,
                        "timepoint1-year FU" = "1-year Follow-up",
                        "modalityEndovascular" = "Endovascular",
                        "modalityMicroneurosurgery" = "Microneurosurgery",
                        "scale(age_Consultation)" = "Age at Time Point, scaled" ,
                        "modalityEndovascular:scale(age_Consultation)" = "Age x Endovascular" ,
                        "modalityMicroneurosurgery:scale(age_Consultation)" = "Age x Microneurosurgery"))
  
  x_min <- min(res_for_plot$conf.low, na.rm = TRUE)
  x_max <- max(res_for_plot$conf.high, na.rm = TRUE)
  x_range <- x_max - x_min
  
  # Textpositionen
  x_left <- x_min - 0.05 * x_range
  x_right <- x_max/1.12 + 0.05 * x_range
  
  # Y-Position unterhalb der Achse
  y_min <- 0.5  # kleiner als die kleinste Term-Position
  y_annot <- 0  # eine Zeile darunter
  
  p <- ggplot(res_for_plot, aes(x = estimate, y = term, color=""))+
    geom_point(show.legend = F, size=3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3,
                   show.legend = F) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    theme_classic() +
    labs(
      x = "IRR",
      y = NULL)+
    ggtitle(title)+
    coord_cartesian(clip = "off") +  
    annotate("text", x = x_left, y = y_annot, 
             label = "Better \nCondition", hjust = 1, vjust = 2, size = 4) +
    annotate("text", x = x_right, y = y_annot, 
             label = "Worse \nCondition", hjust = 0, vjust = 2, size = 4) +
    theme(plot.margin = margin(20, 20, 40, 20),
          text = element_text(size=14),
          axis.text = element_text(size=14))+
    scale_color_manual(values = col_def)
  return(list(data = res_df, plot = p))
}


ForestPlot_MICE_multi <- function(models, model_names, title=NULL, col_defs=colorPalette[c(1,2,3)]){
  
  if(length(models) != length(model_names)){
    stop("Please hand over the same number of models and model names.")
  }
  
  extract_model_df <- function(model, model_name){
    sum_df <- summary(model, conf.int = TRUE, exponentiate = TRUE,
                      conf.level = 0.95, digits = 3) %>% 
      select(-c(`2.5 %`, `97.5 %`, statistic, df))
    
    fmi_df <- model$pooled[, c("term", "fmi")]
    
    left_join(sum_df, fmi_df, by = "term") %>%
      filter(term != "(Intercept)") %>%
      mutate(model = model_name)
  }
  
  all_models_df <- purrr::map2_df(models, model_names, extract_model_df)
  
  # Recode terms
  all_models_df <- all_models_df %>%
    mutate(term_clean = recode(term,
                               "genderMale"= "Male", 
                               "Positive_familial_historyNo"="No Family History" ,
                               "Positive_familial_historyYes" = "Family History",
                               "prf_hypertension_awarenessYes"="Aware of Hypertension",
                               "prf_smokedtobaccoFormer" = "Former Smoker",
                               "prf_smokedtobaccoCurrent" = "Current Smoker",
                               "aneuriskMedium" = "Medium Risk Location",
                               "aneuriskHigh" = "High Risk Location",
                               "scale(maxdiam)" = "IA Size",
                               "typeSaccular - Side-Wall"= "Side-Wall",
                               "sideMidline" = "IA at Mideline",
                               "sideRight" = "IA Right-sided",
                               "modalityEndovascular" = "Endovascular",
                               "modalityMicroneurosurgery" = "Microneurosurgery",
                               "scale(age_tp)" = "Age at Time Point",
                               "modalityEndovascular:scale(age_tp)" = "Age x Endovascular",
                               "modalityMicroneurosurgery:scale(age_tp)" = "Age x Microneurosurgery",
                               "timepointDischarge"= "Discharge",
                               "timepoint1-year FU" = "1-year Follow-up"))
  
  # Wide format: estimates side-by-side
  est_wide <- all_models_df %>%
    select(term_clean, model, estimate) %>%
    tidyr::pivot_wider(names_from = model, 
                       values_from = estimate, 
                       names_prefix = "estimate_")
  
  # Merge back
  all_models_df <- left_join(all_models_df, est_wide, by = "term_clean")
  
  # Create final label: (IRR: x / y)
  if(length(model_names) == 2){
    m1 <- model_names[1]
    m2 <- model_names[2]
    
    all_models_df <- all_models_df %>%
      mutate(term_label = paste0(
        term_clean,
        " (",
        sprintf("%.2f", .data[[paste0("estimate_", m2)]]), " / ",
        sprintf("%.2f", .data[[paste0("estimate_", m1)]]),
        ")"
      ))
  } else {
    stop("Label format currently implemented for exactly 2 models.")
  }
  
  # Fix ordering
  all_models_df$term_label <- factor(all_models_df$term_label,
                                     levels = rev(unique(all_models_df$term_label)))
  
  all_models_df$model <- factor(all_models_df$model, levels = model_names)
  
  ## Plotting ##
  x_min <- min(all_models_df$conf.low)
  x_max <- max(all_models_df$conf.high)
  x_range <- x_max - x_min
  
  x_left  <- x_min - 0.05*x_range
  x_right <- x_max/1.12 + 0.05*x_range
  
  y_annot <- 0
  
  p <- ggplot(all_models_df, aes(x = estimate, y = term_label,
                                 color = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.6), size=2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   position = position_dodge(width = 0.6), height = 0.4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    theme_classic() +
    labs(x = "IRR", y = NULL, color = "Model") +
    ggtitle(title) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = col_defs, name="", limits=rev(model_names)) +
    scale_shape_manual(values = c(17,19), name="", limits=rev(model_names)) +
    coord_cartesian(clip = "off") +
    annotate("text", x = x_left, y = y_annot,
             label = "Better \nCondition", hjust = 1, vjust = 2, size = 4) +
    annotate("text", x = x_right, y = y_annot,
             label = "Worse \nCondition", hjust = 0, vjust = 2, size = 4) +
    theme(plot.margin = margin(20,20,40,20),
          text = element_text(size=14),
          axis.text = element_text(size=14))
  
  return(list(data = all_models_df, plot = p))
}



ForestPlot_GLM <- function(model, title = NULL, col_def = "darkred", exponentiate = TRUE,
                           term1="left", term2="right", measure="OR"){
  # Extract coefficients and standard errors
  sum_df <- broom::tidy(model, conf.int = TRUE, exponentiate = exponentiate)
  
  # Remove intercept & keep order as-is
  res_for_plot <- sum_df %>%
    filter(term != "(Intercept)") %>%
    
    # Recode terms to readable labels
    mutate(term_label = recode(term,
                               "genderMale" = "Male",
                               "Positive_familial_historyNo" = "No Family History",
                               "Positive_familial_historyYes" = "Family History",
                               "prf_hypertension_awarenessYes" = "Aware of Hypertension",
                               "prf_smokedtobaccoFormer" = "Former Smoker",
                               "prf_smokedtobaccoCurrent" = "Current Smoker",
                               "aneuriskMedium" = "Medium Risk Location",
                               "aneuriskHigh" = "High Risk Location",
                               "scale(maxdiam)" = "IA Size",
                               "typeSaccular - Side-Wall" = "Side-Wall",
                               "sideMidline" = "IA at Midline",
                               "sideRight" = "IA Right-sided",
                               "timepointDischarge" = "Discharge",
                               "timepoint1-year FU" = "1-year Follow-up",
                               "modalityEndovascular" = "Endovascular",
                               "modalityMicroneurosurgery" = "Microneurosurgery",
                               "scale(age_tp)" = "Age at Time Point",
                               "scale(age_Consultation)" = "Age at Consultation",
                               "modalityEndovascular:scale(age_tp)" = "Age x Endovascular",
                               "modalityMicroneurosurgery:scale(age_tp)" = "Age x Microneurosurgery",
                               "as.numeric(mR_Discharge)" = "mRS at Discharge",
                               .default = as.character(term))) %>%
    
    # Add (OR: x.xx)
    mutate(term = paste0(term_label,
                         " (", measure, ": ",
                         sprintf("%.2f", estimate), ")"))
  
  # Determine plotting range
  x_min <- min(res_for_plot$conf.low, na.rm = TRUE)
  x_max <- max(res_for_plot$conf.high, na.rm = TRUE)
  x_range <- x_max - x_min
  
  x_left <- x_min - 0.05 * x_range
  x_right <- x_max + 0.05 * x_range
  y_annot <- 0
  
  # Keep order exactly as provided
  res_for_plot$term <- factor(res_for_plot$term, levels = rev(res_for_plot$term))
  
  # Forest plot
  p <- ggplot(res_for_plot, aes(x = estimate, y = term, color = "")) +
    geom_point(show.legend = FALSE, size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.3, show.legend = FALSE) +
    geom_vline(xintercept = ifelse(exponentiate, 1, 0),
               linetype = "dashed", color = "grey") +
    theme_classic() +
    labs(x = ifelse(exponentiate, measure, "Estimate"), y = NULL) +
    ggtitle(title) +
    coord_cartesian(clip = "off") +
    annotate("text", x = x_left, y = y_annot,
             label = term1, hjust = 1, vjust = 2.5, size = 4.25) +
    annotate("text", x = x_right - 0.2 * x_range, y = y_annot,
             label = term2, hjust = 0.1, vjust = 2.5, size = 4.25) +
    theme(plot.margin = margin(20,20,40,20),
          text = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = col_def)
  
  return(list(data = res_for_plot, plot = p))
}


ForestPlot_MICE_multi_n <- function(models, model_names, title = NULL,
                                    col_defs = colorPalette[c(1,2,3,4,5,6)]) {
  
  if (length(models) != length(model_names)) {
    stop("Please hand over the same number of models and model names.")
  }
  
  extract_model_df <- function(model, model_name) {
    sum_df <- summary(
      model,
      conf.int = TRUE,
      exponentiate = TRUE,
      conf.level = 0.95,
      digits = 3
    ) %>%
      select(-c(`2.5 %`, `97.5 %`, statistic, df))
    
    fmi_df <- model$pooled[, c("term", "fmi")]
    
    left_join(sum_df, fmi_df, by = "term") %>%
      filter(term != "(Intercept)") %>%
      mutate(model = model_name)
  }
  
  all_models_df <- purrr::map2_df(models, model_names, extract_model_df)
  
  # Recode terms
  all_models_df <- all_models_df %>%
    mutate(
      term_clean = recode(
        term,
        "genderMale" = "Male",
        "Positive_familial_historyNo" = "No Family History",
        "Positive_familial_historyYes" = "Family History",
        "prf_hypertension_awarenessYes" = "Aware of Hypertension",
        "prf_smokedtobaccoFormer" = "Former Smoker",
        "prf_smokedtobaccoCurrent" = "Current Smoker",
        "aneuriskMedium" = "Medium Risk Location",
        "aneuriskHigh" = "High Risk Location",
        "scale(maxdiam)" = "IA Size",
        "typeSaccular - Side-Wall" = "Side-Wall",
        "sideMidline" = "IA at Mideline",
        "sideRight" = "IA Right-sided",
        "modalityEndovascular" = "Endovascular",
        "modalityMicroneurosurgery" = "Microneurosurgery",
        "scale(age_tp)" = "Age at Time Point",
        "modalityEndovascular:scale(age_tp)" = "Age x Endovascular",
        "modalityMicroneurosurgery:scale(age_tp)" = "Age x Microneurosurgery",
        "timepointDischarge" = "Discharge",
        "timepoint1-year FU" = "1-year Follow-up"
      )
    )
  
  # Wide format: estimates side-by-side
  est_wide <- all_models_df %>%
    select(term_clean, model, estimate) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = estimate,
      names_prefix = "estimate_"
    )
  
  # Merge back
  all_models_df <- left_join(all_models_df, est_wide, by = "term_clean")
  
  all_models_df <- all_models_df %>%
    rowwise() %>%
    mutate(
      term_label = {
        row_vals <- cur_data()
        paste0(
          term_clean,
          " (",
          paste(
            sapply(
              model_names,
              function(m) {
                sprintf(
                  "%s: %.2f",
                  m,
                  row_vals[[paste0("estimate_", m)]]
                )
              }
            ),
            collapse = " | "
          ),
          ")"
        )
      }
    ) %>%
    ungroup()
  
  # Fix ordering
  all_models_df$term_label <- factor(
    all_models_df$term_label,
    levels = rev(unique(all_models_df$term_label))
  )
  
  all_models_df$model <- factor(
    all_models_df$model,
    levels = rev(model_names)
  )
  
  
  ## Plotting ##
  x_min <- min(all_models_df$conf.low)
  x_max <- max(all_models_df$conf.high)
  x_range <- x_max - x_min
  
  x_left  <- x_min - 0.05 * x_range
  x_right <- x_max / 1.12 + 0.05 * x_range
  
  y_annot <- 0
  
  p <- ggplot(
    all_models_df,
    aes(x = estimate, y = term_label, color = model, shape = model)
  ) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      position = position_dodge(width = 0.6),
      height = 0.4
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    theme_classic() +
    labs(x = "IRR", y = NULL, color = "Model") +
    ggtitle(title) +
    theme(legend.position = "bottom") +
    scale_color_manual(
      values = col_defs[seq_along(model_names)],
      name = "",
      limits = (model_names)
    ) +
    scale_shape_manual(
      values = c(17, 19, 15, 18, 8, 7)[seq_along(model_names)],
      name = "",
      limits = (model_names)
    ) +
    coord_cartesian(clip = "off") +
    annotate(
      "text",
      x = x_left,
      y = y_annot,
      label = "Better \nCondition",
      hjust = 1,
      vjust = 2,
      size = 4
    ) +
    annotate(
      "text",
      x = x_right,
      y = y_annot,
      label = "Worse \nCondition",
      hjust = 0,
      vjust = 2,
      size = 4
    ) +
    theme(
      plot.margin = margin(20, 20, 40, 20),
      text = element_text(size = 14),
      axis.text = element_text(size = 14)
    )
  
  return(list(data = all_models_df, plot = p))
}


# get individual and pooled estimates per function:
# function to extract one parameter from all fits and add pooled estimate
forestplot_mice_param <- function(fits_list, param, exponentiate = TRUE, title = NULL, xlim = c(-1, 1)) {
  
  # 1. Extract from single fits
  extract_param <- function(fit, param){
    tidy_fit <- broom.mixed::tidy(fit, effects = "fixed")
    df_resid <- df.residual(fit)
    res <- tidy_fit %>%
      dplyr::filter(term == param) %>%
      dplyr::select(term, estimate, std.error) %>%
      dplyr::mutate(
        conf.low = estimate - qnorm(0.975) * std.error,
        conf.high = estimate + qnorm(0.975) * std.error
      )
    return(res)
  }
  
  all_params <- lapply(fits_list, extract_param, param = param)
  all_params_df <- dplyr::bind_rows(all_params, .id = "run_id") %>%
    dplyr::mutate(run = paste0("Imputation_", run_id),
                  source = "single") %>%
    dplyr::select(-run_id)
  
  # 2. Add pooled estimate
  fit_pool <- pool(fits_list)
  sum_pool <- summary(fit_pool, conf.int = TRUE, exponentiate = FALSE)
  
  pool_df <- sum_pool %>%
    dplyr::filter(term == param) %>%
    dplyr::transmute(
      term,
      estimate,
      std.error = std.error,
      conf.low = `2.5 %`,
      conf.high = `97.5 %`,
      run = "Pooled",
      source = "pooled"
    )
  
  # 3. Combine and ensure factor order
  plot_df <- dplyr::bind_rows(all_params_df, pool_df)
  
  # 4. Exponentiate if requested
  if (exponentiate) {
    plot_df <- plot_df %>%
      dplyr::mutate(
        estimate = exp(estimate),
        conf.low = exp(conf.low),
        conf.high = exp(conf.high)
      )
    xlab <- "IRR"
    vline <- 1
  } else {
    xlab <- "Log(IRR)"
    vline <- 0
  }
  
  # 5. Add label text
  plot_df <- plot_df %>%
    dplyr::mutate(
      label = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
    )
  
  # 6. Factor order: pooled at bottom
  levels_order <- c(unique(all_params_df$run), "Pooled")
  plot_df <- plot_df %>%
    dplyr::mutate(run = factor(run, levels = levels_order)) %>%
    mutate(source = factor(source, levels = c("single", "pooled")))
  
  # 7. Define position for text (to the right of all CIs)
  max_x <- max(plot_df$conf.high, na.rm = TRUE)
  text_x <- max_x + 0.05 * diff(range(plot_df$conf.high, plot_df$conf.low))  # spacing
  
  # 8. Forest plot
  p <- ggplot(plot_df, aes(x = estimate, y = run, color = source, shape = source)) +
    geom_point(aes(size = source)) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3) +
    geom_vline(xintercept = vline, linetype = "dashed") +
    geom_text(aes(x = text_x, label = label), hjust = 0, color = "black", size = 3.2) +  # right-aligned text
    labs(
      x = xlab, y = NULL,
      title = ifelse(is.null(title), paste("Forest Plot for", param), title),
      color = "", shape = "", size = ""
    ) +
    scale_color_manual(values = c("single" = "grey50", "pooled" = "black")) +
    scale_shape_manual(values = c("single" = 16, "pooled" = 18)) +
    scale_size_manual(values = c("single" = 2, "pooled" = 4)) +
    scale_y_discrete(limits = rev) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      element_text()
    ) +
    xlim(xlim[1], text_x * 1.2)  # extend x-axis to fit text
  
  return(list(data = plot_df, plot = p))
}

#################################################################################
## 3. ATE function
#################################################################################
# Function to compute pooled pairwise contrasts (RR) from glmmTMB fits
pool_contrasts <- function(fits_list, factor_name) {
  # 1. Number of imputations
  m <- length(fits_list)
  
  # 2. Parameter names (fixed effects)
  param_names <- names(fixef(fits_list[[1]])$cond)
  
  # 3. Variance-covariance matrices per fit
  vcovs_list <- lapply(fits_list, function(fit) {
    V <- vcov(fit)$cond   # conditional part only
    rownames(V) <- colnames(V) <- names(fixef(fit)$cond)
    return(V)
  })
  
  # 4. Coefficient matrix (rows = parameters, cols = imputations)
  Q_mat <- sapply(fits_list, function(fit) fixef(fit)$cond[param_names])
  
  # 5. Array of within-imputation variances
  U_array <- array(NA, dim = c(length(param_names), length(param_names), m),
                   dimnames = list(param_names, param_names, NULL))
  for (i in 1:m) {
    U_array[,,i] <- vcovs_list[[i]][param_names, param_names]
  }
  
  # 6. Rubin’s rules
  Qbar <- rowMeans(Q_mat)                   # pooled coefficients
  B <- cov(t(Q_mat))                        # between-imputation variance
  Ubar <- apply(U_array, c(1,2), mean)      # within-imputation variance
  Tmat <- Ubar + (1 + 1/m) * B              # total variance-covariance
  
  # 7. Factor levels
  factor_levels <- grep(paste0("^", factor_name), param_names, value = TRUE)
  f_levels <- levels(fits_list[[1]]$frame[[factor_name]])  # includes reference
  ref_level <- f_levels[1]                                 # reference level
  other_levels <- f_levels[-1]
  
  # 8. All pairwise comparisons, but always put reference on the right
  comparisons <- c(
    lapply(other_levels, function(lvl) c(lvl, ref_level)),  # each vs reference
    combn(other_levels, 2, simplify = FALSE)                # non-reference vs non-reference
  )
  
  # 9. Compute contrasts
  result_list <- lapply(comparisons, function(pair) {
    L <- rep(0, length(Qbar))
    for (lvl in pair) {
      if (lvl == ref_level) next
      coef_name <- paste0(factor_name, lvl)
      L[param_names == coef_name] <- ifelse(pair[1] == lvl, 1, -1)
    }
    
    est <- sum(L * Qbar)
    se  <- sqrt(t(L) %*% Tmat %*% L)
    CI_low <- exp(est - 1.96*se)
    CI_high <- exp(est + 1.96*se)
    RR <- exp(est)
    
    data.frame(
      Comparison = paste(pair[1], "vs", pair[2]),
      Estimate = est,
      SE = se,
      CI_Lower = CI_low,
      CI_Upper = CI_high,
      RR = RR
    )
  })
  
  # 10. Bind results into a single table
  ATE_table <- dplyr::bind_rows(result_list)
  return(ATE_table)
}

#################################################################################
# 4. Functions for model fit statistics:
###############################################################################

## Mean BIC for imputed data sets:
bic_summary <- function(fitlist) {
  bics <- sapply(fitlist, function(fit) {
    # BIC() sollte für glmmTMB funktionieren
    BIC(fit)
  })
  list(
    bic_values = bics,
    mean_bic = mean(bics),
    median_bic = median(bics),
    sd_bic = sd(bics),
    m = length(bics)
  )
}
