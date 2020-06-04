#=============================================================================#
# Meta Analyses
#
# Functions to organize data, conduct meta-analyses, and create plots
#=============================================================================#


#' remove_low_cell_count
#'
#' low cell counts of outcome, exposure, and study_gxe cause regression errors. If exposure is categorical, 3-way tabulate exposure/outcome/study_gxe, remove studies with cell counts <= 1. If numeric, 2-way tabulate outcome/study_gxe, remove studies with cell counts == 0. variables 'outcome' and 'study_gxe' are hardcoded
#'
#' @param ds dataset
#' @param exposure string containing exposure variable
#' @param min_cell_size minimum cell size allowed (numeric)
#'
#' @return returns processed input dataframe
#' @export
#'
#' @examples remove_low_cell_count(ds, exposure, min_cell_size = 1)
remove_low_cell_count <- function(ds, min_cell_size = 1) {
  
  # only binary variables are modeled categorically
  # all else are modeled as continuous (including Q4 variables)
  categorical <- ifelse(length(table(ds[, exposure])) <= 2, T, F)
  
  if(categorical) {
    drops <- data.frame(table(ds$outcome, ds[, exposure], ds$study_gxe)) %>% 
      filter(Freq <= min_cell_size)
    tmp <- ds %>% 
      filter(!study_gxe %in% drops$Var3) %>% 
      mutate(study_gxe = fct_drop(study_gxe))
  } else {
    drops <- data.frame(table(ds$outcome, ds$study_gxe)) %>% 
      filter(Freq <= min_cell_size)
    tmp <- ds %>% 
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
#' @examples get_counts_by_outcome(gxe, outcome, aspirin, study_gxe)
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


#' get_estimates_e_by_group
#'
#' Fit regression for exposure main effect by group. Input data should not include studies with zero counts when tabulating exposure/outcome/study.
#'
#' @param dat Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable e.g. study_gxe
#' @param ... Adjustment covariates
#'
#' @return Tibble with exposure GLM estimates/se/stats/pval/95%CI by group
#' @export
#'
#' @examples run_glm_e_by_group(gxe, outcome, asp_ref, study_gxe, age_ref_imp, sex, PC1, PC2, PC3)
get_estimates_e_by_group <- function(dat, outcome, exposure, group, ...) {

  # model covariates as vector
  covariates <- c(...)

  # group data, perform grouped glm
  dat <- dat %>%
    group_by( .data[[group]] )

  glm_formula <- reformulate(termlabels = c( exposure, covariates) , response = outcome )
  results_beta <- dplyr::do(dat, broom::tidy(glm(glm_formula, data = . , family = 'binomial')))
  results_ci   <- dplyr::do(dat, broom::confint_tidy(glm(glm_formula, data = . , family = 'binomial')))

  results <- dplyr::bind_cols(results_beta, results_ci[,c(2,3)]) %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl(exposure, term))
  return(results)

}



#' meta_analysis_wrapper
#'
#' Wrapper function to run meta-analysis, create forest and funnel plots. Uses package 'meta'.
#'
#' @param dat Data frame containing counts and regression estimates BY GROUP
#' @param exposure string containing exposure variable
#' @param covariates vector of adjustment covariates
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
#' @examples meta_analysis_wrapper(gxe_meta, 'asp_ref', c('sex'), "title", "age_ref_imp_sex_study_gxe", 13, 8.5, 8, 8.5)
meta_analysis_wrapper <- function(dat, exposure, covariates, forest_plot_title, filename_suffix, forest_height, forest_width, funnel_height, funnel_width) {
  
  output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  
  tmp <- remove_low_cell_count(dat)
  count_gxe_subset <- get_counts_outcome_by_group(tmp, 'outcome', 'study_gxe') %>% 
    mutate(study_gxe = as.character(study_gxe))
  count_gxe_all <- full_join(study_design, count_gxe_subset, 'study_gxe')
  
  gxe_glm <- get_estimates_e_by_group(tmp, "outcome", exposure, group = "study_gxe", covariates) %>% 
    mutate(study_gxe = as.character(study_gxe))
  # forest_plot_title <- paste0("All subjects\nOutcome ~ ", exposure, " + ", paste0(covariates, collapse = " + "))
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
  
  # png(paste0("/media/work/gwis/posthoc/", exposure, "/meta_analysis_", exposure,  "_", filename_suffix, ".png"), height = forest_height, width = forest_width, units = 'in', res = 150)
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
#' @examples run_glm_gxe_by_group(gxe, outcome, asp_ref, study_gxe, X5.40252294, age_ref_imp, sex, PC1, PC2, PC3)
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


