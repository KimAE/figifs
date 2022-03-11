#=============================================================================#
# Meta Analyses
#
# Functions to organize data, conduct meta-analyses, and create plots
#=============================================================================#


#' remove_low_cell_count
#'
#' low cell counts of outcome, exposure, and study_gxe cause regression errors. If exposure is categorical, 3-way tabulate exposure/outcome/study_gxe, remove studies with cell counts <= 1. If numeric, 2-way tabulate outcome/study_gxe, remove studies with cell counts == 0. variables 'outcome' and 'study_gxe' are hardcoded
#'
#' @param dat dataset
#' @param exposure string containing exposure variable
#' @param min_cell_size minimum cell size allowed (numeric)
#'
#' @return returns processed input dataframe
#' @export
#'
remove_low_cell_count <- function(dat, exposure, min_cell_size = 1) {
  
  # only binary variables are modeled categorically
  # all else are modeled as continuous (including Q4 variables)
  categorical <- ifelse(length(table(dat[, exposure])) <= 2, T, F)
  
  if(categorical) {
    drops <- data.frame(table(dat$outcome, dat[, exposure], dat$study_gxe)) %>% 
      filter(Freq <= min_cell_size)
    tmp <- dat %>% 
      filter(!study_gxe %in% drops$Var3) %>% 
      mutate(study_gxe = fct_drop(study_gxe))
  } else {
    drops <- data.frame(table(dat$outcome, dat$study_gxe)) %>% 
      filter(Freq <= min_cell_size)
    tmp <- dat %>% 
      filter(!study_gxe %in% drops$Var2) %>% 
      mutate(study_gxe = fct_drop(study_gxe))
  }
  return(tmp)
}


#' get_counts_outcome_by_group
#'
#' Tabulate outcome by study/platform to inclusion in forest plots. Input data should not include studies with zero counts when tabulating exposure/outcome/study.
#'
#' @param dat Input data
#' @param outcome Outcome variable
#' @param group Group variable e.g. study_gxe
#'
#' @return Tibble with case/control counts by group.
#' @export
#'
get_counts_outcome_by_group <- function(dat, outcome, group) {
  
  # get counts for rmeta function
  dat <- dat %>%
    dplyr::group_by( .data[[outcome]], .data[[group]] ) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    tidyr::spread(key = {{ outcome }}, value = count) %>%
    dplyr::rename(Control = `0`, Case = `1`) %>%
    dplyr::mutate(N = Control + Case)
  return(dat)
  
}



#' get_estimates_e_by_group_pwalk
#'
#' Fit regression for exposure main effect by group. Currently group should ALWAYS be "study_gxe". Ensure data doesn't contain studies with empty cells. Ensure that model_formula is appropriate (shouldn't include study_gxe, for example)
#'
#' @param dataset Dataset
#' @param exposure Exposure variable
#' @param model_formula string character of glm formula
#'
#' @return Tibble with exposure GLM estimates/se/stats/pval/95%CI, by group
#' @export
#'
get_estimates_e_by_group_pwalk <- function(dataset, exposure, model_formula) {
  
  # group by study_gxe
  dataset <- dataset %>%
    group_by(study_gxe)
  
  # run glm with model_formula
  results_beta <- dplyr::do(dataset, broom::tidy(glm(model_formula, data = . , family = 'binomial'), conf.int = T))
  
  # clean up output, return exposure main effect
  results <- results_beta %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl(exposure, term)) %>% 
    dplyr::arrange(study_gxe)
  return(results)
}






#' meta_analysis_wrapper_pwrap
#'
#' Wrapper function to run meta-analysis, create forest and funnel plots. Uses package 'meta'.
#'
#' @param dataset Data frame containing counts and regression estimates BY GROUP
#' @param exposure string containing exposure variable
#' @param model_formula string glm formula
#' @param subset_var string analysis subset ('all' means everyone)
#' @param output_dir string output directory
#' @param forest_plot_title plot title
#' @param filename_suffix filename suffix following exposure (usually covariates)
#' @param forest_height png height
#' @param forest_width png width
#' @param funnel_height png height
#' @param funnel_width png width
#' 
#' @return Outputs png files for forest and funnel plots
#' @export
#'
meta_analysis_wrapper_pwrap <- function(dataset, exposure, model_formula, subset_var, output_dir, forest_plot_title, filename_suffix, forest_height, forest_width, funnel_height, funnel_width) {
  
  # analysis subset (matches GxEScanR input)
  exposure_subset <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_glm.rds"))[, 'vcfid']
  
  # create study_design data.frame
  study_design <- dataset %>%
    dplyr::select(study_gxe, study_design) %>%
    filter(!duplicated(.)) %>%
    arrange(study_gxe)
  
  # create temporary analysis dataset subset
  tmp <- dataset %>% 
    filter(vcfid %in% exposure_subset) %>% 
    filter(case_when(subset_var == "female" ~ sex == 0, 
                     subset_var == "male" ~ sex == 1,
                     subset_var == "proximal" ~ cancer_site_sum2 == "proximal" | outcome == 0,
                     subset_var == "distal" ~ cancer_site_sum2 == "distal" | outcome == 0,
                     subset_var == "rectal" ~ cancer_site_sum2 == "rectal" | outcome == 0,
                     TRUE ~ TRUE))
  dataset_tmp <- remove_low_cell_count(tmp, exposure)
  
  # get study sample sizes
  count_gxe_subset <- get_counts_outcome_by_group(dataset_tmp, 'outcome', 'study_gxe') %>%
    mutate(study_gxe = as.character(study_gxe))
  count_gxe_all <- full_join(study_design, count_gxe_subset, 'study_gxe')
  
  # run GLM by study_gxe, merge with counts table for meta-analysis
  gxe_glm <- get_estimates_e_by_group_pwalk(dataset_tmp, exposure, model_formula) %>%
    mutate(study_gxe = as.character(study_gxe))
  gxe_meta <- dplyr::full_join(count_gxe_all, gxe_glm, 'study_gxe')
  
  # don't display warnings
  oldw <- getOption("warn")
  options(warn = -1)
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=gxe_meta,
                                studlab=paste(study_gxe),
                                comb.fixed = FALSE,
                                comb.random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                byvar = study_design)
  options(warn = oldw)
  
  png(paste0(output_dir, "meta_analysis_", exposure,  "_", filename_suffix, ".png"), height = forest_height, width = forest_width, units = 'in', res = 150)
  meta::forest(results_meta,
               layout = "JAMA",
               # text.predict = "95% CI",
               # col.predict = "black",
               leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
               digits.addcols=0,
               study.results=T,
               prediction = F,
               col.random = 'red')
  grid.text(forest_plot_title, 0.5, .98, gp=gpar(cex=1))
  dev.off()
  
  png(paste0(output_dir, "funnel_plot_", exposure,  "_", filename_suffix, ".png"), height = funnel_height, width = funnel_width, units = 'in', res = 150)
  oldw <- getOption("warn")
  options(warn = -1)
  meta::funnel(results_meta, sm="OR", studlab = T, pos = 4, col.random = 'red')
  options(warn = oldw)
  dev.off()
}




#' get_estimates_gxe_by_group
#'
#' Function to get GLM estimates by group (e.g. study_gxe). Make sure data is properly formatted yeah (maybe add info on what that means exactly..). This is an alternate version that returns the gxe term in the summary table
#'
#' @param df Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable (typically study_g)
#' @param dosage SNP for interaction with E
#' @param ... Adjustment Covariates
#'
#' @return Tibble with GLM estimates/se/stats/pval/95%CI by group
#' @export
#'
get_estimates_gxe_by_group <- function(df, outcome, exposure, group, dosage, ...) {

  outcome <- enquo(outcome)
  exposure <- enquo(exposure)
  group <- enquo(group)
  dosage <- enquo(dosage)
  covariates <- enquos(...)

  # data prep
  df <- format_data_meta_analysis(df, !! outcome, !! exposure, !! group) %>%
    dplyr::group_by(!! group)

  glm_formula <- reformulate(termlabels = c(paste0(quo_name(exposure),"*", quo_name(dosage)),
                                            as.vector(unlist(lapply(covariates, quo_name)))),
                             response = quo_name(outcome))

  results_beta <- dplyr::do(df, broom::tidy(glm(glm_formula, data = . , family = 'binomial')))
  results_ci   <- dplyr::do(df, broom::confint_tidy(glm(glm_formula, data = . , family = 'binomial')))

  results <- dplyr::bind_cols(results_beta, results_ci) %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl(paste0(quo_name(exposure), ":"), term))
  return(results)

}









#' meta_analysis_execute
#' 
#' call meta_analysis_wrapper_pwrap. makes it easier to run this for all the exposures
#' also makes it easier to run on multiple adjustment covariate sets
#'
#' @param dataset 
#' @param exposure 
#' @param hrc_version 
#' @param covariates 
#' @param output_dir 
#' @param filename_suffix2 
#' @param nosex 
#'
#' @return
#' @export
# meta_analysis_execute <- function(dataset, exposure, hrc_version, covariates, path, filename_suffix2 = "", nosex=F) {
#   
#   # assemble data.frame of arguments for meta_analysis_wrapper function
#   covariates_sorted <- sort(covariates)
#   
#   ## stratifying variables - ALWAYS these. change in future if necessary
#   stratifying_vars <- c("all", "proximal", "distal", "rectal", "female", "male")
#   
#   ## model formulas
#   model_formula_all   <- Reduce(paste, deparse(reformulate(c(exposure, covariates_sorted), response = 'outcome')))
#   model_formula_nosex <- Reduce(paste, deparse(reformulate(c(exposure, covariates_sorted[!covariates_sorted %in% c("sex")]), response = 'outcome')))
#   model_formula <- c(rep(model_formula_all, 4), rep(model_formula_nosex, 2))
#   
#   ## figure titles
#   forest_plot_title <- purrr::map2_chr(.x = stratifying_vars, .y = model_formula, ~ paste0(.x, "\n", .y))
#   
#   ## filename
#   
#   filename <- paste0(glue("meta_analysis_{exposure}_{hrc_version}_{glue_collapse(covariates, sep = '_')}_strata_"), stratifying_vars, ".png")
#   
#   ## argument data.frame
#   meta_analysis_args <- data.frame(
#     exposure = rep(exposure, length(model_formula)),
#     model_formula = model_formula,
#     subset_var = stratifying_vars,
#     output_dir = rep(path, length(model_formula)),
#     forest_plot_title = forest_plot_title,
#     filename_suffix = filename,
#     forest_height = 17, 
#     forest_width = 8.5, 
#     funnel_height = 8, 
#     funnel_width = 8.5, 
#     stringsAsFactors = F)
#   
#   # execute
#   purrr::pwalk(meta_analysis_args, meta_analysis_wrapper_pwrap, dataset = dataset)
# }
# 
# 






#' create_forest_plot
#'
#' @param data_epi 
#' @param exposure 
#' @param covariates 
#' @param hrc_version 
#' @param path 
#' @param forest_height 
#' @param forest_width 
#' @param funnel_height 
#' @param funnel_width 
#' @param strata 
#' @param categorical 
#'
#' @return
#' @export
#'
#' @examples
create_forest_plot <- function(data_epi, exposure, covariates, hrc_version, path, forest_height = 17, forest_width = 8.5, funnel_height = 8, funnel_width = 8.5, strata = 'all', categorical = T) {
  
  # note - data_epi should be the exposure subset you submitted to gxescanR
  # make sure study_gxe is not in the covariates vector
  covariates <- covariates[which(covariates != 'study_gxe')]
  
  if (strata == 'female') {
    data_epi <- dplyr::filter(data_epi, sex == 0)
    covariates <- covariates[which(!covariates %in% c('sex'))]
  } else if (strata == 'male') {
    data_epi <- dplyr::filter(data_epi, sex == 1)
    covariates <- covariates[which(!covariates %in% c('sex'))]
  } else if (strata == 'proximal') {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 'proximal' | data_epi$outcome == 0)
  } else if (strata == 'distal') {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 'distal' | data_epi$outcome == 0)
  } else if (strata == 'rectal') { 
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 'rectal' | data_epi$outcome == 0)
  }
  
  # some subsets generate empty cells, need to remove them
  if(categorical == T) {
    drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome, data_epi[, exposure])) %>% 
      dplyr::filter(Freq == 0) %>% 
      dplyr::pull(Var1) %>% 
      unique(.)
    
    data_epi <- data_epi %>% 
      dplyr::filter(!study_gxe %in% drops)
  } else if(categorical == F) {
    drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome)) %>% 
      dplyr::filter(Freq <= 5) %>% 
      dplyr::pull(Var1) %>% 
      unique(.)
    
    data_epi <- data_epi %>% 
      dplyr::filter(!study_gxe %in% drops)
  }
  
 
  
  # create model term
  # model_formula <- Reduce(paste, deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome')))
  # model_formula <- deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome'))
    model_formula <- glue("outcome ~ {exposure} + {glue_collapse(sort(covariates), sep = '+')}")


  # create study_design data.frame
  study_design <- data_epi %>%
    dplyr::select(study_gxe, study_design) %>%
    filter(!duplicated(.)) %>%
    arrange(study_gxe)
  
  glm_out <- data_epi %>%
    tidyr::nest(data = -study_gxe) %>% 
    dplyr::mutate(
      fit = purrr::map(data, ~ glm(model_formula, data = .x, family = 'binomial')),
      tidied = purrr::map(fit, ~ tidy(.x)), 
      quality = purrr::map(fit, ~ glance(.x))
    ) %>% 
    dplyr::select(-data, -fit) %>% 
    tidyr::unnest(tidied) %>% tidyr::unnest(quality) %>% 
    dplyr::filter(grepl(exposure, term)) %>%
                  #null.deviance > quantile(null.deviance, 0.01)) %>% 
    dplyr::arrange(study_gxe)
  
  # include study sample sizes and study design information
  meta_input <- get_counts_outcome_by_group(data_epi, 'outcome', 'study_gxe') %>%
    mutate(study_gxe = as.character(study_gxe)) %>% 
    inner_join(study_design, 'study_gxe') %>% 
    full_join(glm_out, 'study_gxe')
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=meta_input,
                                studlab=paste(study_gxe),
                                fixed = TRUE,
                                random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                subgroup = study_design)
  
  # plot_title = glue(strata, "\n{model_formula}", .na = "All")
  plot_title = glue(strata, " (N=", sum(glm_out$nobs), ")", "\n{model_formula}")
  
  png(glue("{path}/forest_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"), height = forest_height, width = forest_width, units = 'in', res = 150)
  meta::forest(results_meta,
               layout = "JAMA",
               # text.predict = "95% CI",
               # col.predict = "black",
               leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
               digits.addcols=0,
               study.results=T,
               prediction = F,
               col.random = 'red', 
               xlim = c(0.25, 4))
  grid.text(plot_title, 0.5, .98, gp=gpar(cex=1))
  dev.off()
  
  png(glue("{path}/funnel_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"), height = funnel_height, width = funnel_width, units = 'in', res = 150)
  meta::funnel(results_meta, sm="OR", studlab = T, pos = 4, col.random = 'red')
  dev.off()
  
  #return(glue("{path}/forest_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"))
}







#' create_glm_stratified_plot
#'
#' Uses package 'effects' to create a stratified plot of GxE interactions. Incidentally, whenever you specify interactions in glm, ALWAYS USE G*E (otherwise, specify on 'gxe' flag)
#'
#' @param model GLM output with higher order (interaction) term in the format G*E
#' @param G Genotype dosage variable (string)
#' @param E Exposure variable (string)
#' @param gxe T/F, just to tell function if GxE or ExG (sorry)
#'
#' @return A stratified GxE plot (NOT exported as *.png)
#' @export
#'
#' @examples create_glm_stratified_plot(model1, "X5.12345678", "exposure", gxe = T)
create_glm_stratified_plot <- function(model, G, E, gxe) {
  
  # how you define list names right off the bat? grr
  xlevel_list <- list(c(0,1,2))
  if(gxe == T) {
    names(xlevel_list) <- G
  } else {
    names(xlevel_list) <- E # bad but let's get on with it
  }
  
  model_eff <- effect(paste(G,E, sep = "*"),
                      model,
                      xlevels=xlevel_list, # might have to add more here, depending on coding of E
                      se=TRUE,
                      confidence.level=.95,
                      typical=mean)
  
  model_eff <- as.data.frame(model_eff) # factors should be defined before fitting model
  
  ggplot(data=model_eff, aes(x = !! sym(G), y = fit, group = !! sym(E))) +
    # coord_cartesian(ylim = c(0.25,.5)) +
    geom_line(size=2, aes(color=!! sym(E))) +
    ylab("Model Fit")+
    xlab(paste0(G, " Dosage")) +
    ggtitle(paste0("GxE Stratified Analysis - ", G, "*", E)) +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    scale_fill_grey()
  
}




#' pooled_analysis_glm
#'
#' conduct pooled analysis of exposure main effect. Assess interaction by sex and study_design
#'
#' @param ds dataset
#' @param exposure string containing exposure variable
#' @param covariates vector of adjustment covariates
#' @param strata string containing stratifying variable (either sex or study_design)
#' @param filename_suffix string containing filename suffix
#' @param output_dir string output directory
#'
#' @return HTML file
#' @export
#'
#' @examples pooled_analysis_glm(figi, 'asp_ref', c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'), strata = 'sex', 'testing')
pooled_analysis_glm <- function(ds, exposure, hrc_version, covariates, strata = c('sex', 'study_design', 'bmic3'), filename_suffix, output_dir) {

  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  strata <- match.arg(strata)
  number_of_levels <- nlevels(factor(ds[, strata]))
  ds[, 'strata_num'] <- as.numeric(factor(ds[, strata]))-1
  
  # formulas
  covariates_nostrata <- paste0(covariates[! covariates %in% strata], collapse = " + ")
  formula_all               <- paste0("outcome ~ ", exposure, " + ", paste0(covariates, collapse = "+"))
  formula_original          <- paste0("outcome ~ ", exposure, " + ", covariates_nostrata)
  formula_interaction       <- paste0("outcome ~ ", exposure, " * ", strata, " + ", covariates_nostrata)
  
  # fit models (original, stratified, interaction)
  out <- list()
  
  m_original    <- glm(formula_all, data = ds, family = 'binomial')
  out[['original']] <- m_original
  
  for(level in seq(number_of_levels) - 1) {
    # fit models
    if(strata == "cancer_site_sum2") {
      index_vector <- which(ds[,'strata_num'] == level | ds[,'outcome'] == 0)
    } else {
      index_vector <- which(ds[,'strata_num'] == level)
    }
    
    m_subset    <- glm(formula_original,          data = ds[index_vector,], family = 'binomial')
    out[[paste0(strata, "_", as.character(level))]] <- m_subset
  }
  
  m_interaction <- glm(formula_interaction, data = ds, family = 'binomial')
  out[['interaction']] <- m_interaction
  coefs <- lapply(out, function(x) (exp(coef(x))))
  list_of_samplesizes <- lapply(out, function(x) paste0(c("Co=", "Ca="), as.character(table(x$model$outcome)), collapse = ','))
  
  # output to stargazer html file
  if(strata == 'sex') {
    col_label = paste0(c("All", "Female", "Male", "Interaction"), " (", list_of_samplesizes, ")")
  } else if(strata == 'study_design') {
    col_label = paste0(c("All", "Cohort", "Case-Control", "Interaction"), " (", list_of_samplesizes, ")")
  } else if(strata == 'bmic3') {
    col_label = paste0(c("All", "Normal", "Overweight", "Obese"), " (", list_of_samplesizes, ")")
  } else {
    col_label = c("All", seq(0, number_of_levels - 1))
  }
  
  # save output html from stargazer
  out_html <- stargazer_helper(out,
                               title=paste0(gsub('\\_', '\\\\_', strata), " stratified ", gsub('\\_', '\\\\_', exposure), " main effects"), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=c("(Study specific and PC estimates omitted from table)"), single.row = T)
  
  # write object to html file
  cat(paste(out_html, collapse = "\n"), "\n",
      file = paste0(output_dir, "main_effects_pooled_analysis_", exposure, "_", hrc_version, "_", filename_suffix, ".html"), append = F)
}


# pooled_analysis_glm(figi, 'asp_ref', c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'), strata = 'sex', 'testing')




#' pooled_analysis_multinom
#' 
#' conduct pooled analysis of exposure main effect. Fits multinomial models to test differences in tumor subsite
#'
#' @param ds dataset
#' @param exposure string containing exposure variable
#' @param covariates vector of adjustment covariates
#' @param filename_suffix string containing filename suffix
#' @param output_dir string output directory
#'
#' @return HTML file
#' @export
#'
#' @examples pooled_analysis_multinom(ds, exposure, covariates, 'test')
pooled_analysis_multinom <- function(ds, exposure, hrc_version, covariates, filename_suffix, output_dir) {
  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  # no output in messages from multinom function
  # quiet <- function(x) {
  #   sink(tempfile())
  #   on.exit(sink())
  #   invisible(force(x))
  # }
  
  tmp <- ds %>% 
    mutate(outcome_multinomial = factor(ifelse(outcome == 0 & is.na(cancer_site_sum2), "control",
                                               ifelse(cancer_site_sum2 == "proximal", "proximal",
                                                      ifelse(cancer_site_sum2 == "distal", "distal",
                                                             ifelse(cancer_site_sum2 == "rectal", "rectal", NA)))), 
                                        levels = c("control", "proximal", "distal", "rectal")))
  
  
  stratified_formula <- as.formula(paste0("outcome_multinomial ~ ", exposure, " + ", paste0(covariates, collapse = " + ")))
  
  # quiet(x1 <- multinom(stratified_formula, data = tmp))
  x1 <- multinom(stratified_formula,  data = tmp)
  
  # summary(x1)
  coefficients <- coef(x1)
  coefname_e <- x1$coefnames[2]
  pval_prox_dist <- linearHypothesis(x1,paste0("proximal:", coefname_e," = distal:", coefname_e),test="Chisq")[2,3]
  pval_prox_rect <- linearHypothesis(x1,paste0("proximal:", coefname_e," = rectal:", coefname_e),test="Chisq")[2,3]
  pval_dist_rect <- linearHypothesis(x1,paste0("distal:", coefname_e," = rectal:", coefname_e),test="Chisq")[2,3]
  
  sample_sizes = as.character(table(tmp$outcome_multinomial))
  
  out_html <- stargazer_helper(x1,
                               title=paste0("Tumor site multinomial logistic regression (Control N = ", sample_sizes[1], ")"), 
                               # column.labels=rep("vs. controls", 3),
                               column.labels=paste0("(N=", sample_sizes[2:4], ")"),
                               coef=list(exp(coefficients)), 
                               notes=c("(Study specific and PC estimates omitted from table)"), single.row = T, 
                               add.lines = list(c("Proximal vs. Distal pval", "", formatC(pval_prox_dist, format = "e", digits = 2)),
                                                c("Proximal vs. Rectal pval", "", formatC(pval_prox_rect, format = "e", digits = 2)),
                                                c("Distal vs. Rectal pval",   "", formatC(pval_dist_rect, format = "e", digits = 2))))
  cat(paste(out_html, collapse = "\n"), "\n",
      file = paste0(output_dir, "main_effects_pooled_analysis_", exposure, "_", hrc_version, "_", filename_suffix, ".html"), append = F)
  
}



reload_figifs <- function() {
  devtools::reload(pkgload::inst("figifs"))
}








