#=============================================================================#
# Functions for posthoc analysis
#=============================================================================#

#' stargazer_helper 
#' 
#' convenience function to call stargazer with commonly used arguments
#'
#' @param ... exposure specific arguments
#' @param title string containing title of plot
#' @param column.labels vector of column labels
#' @param coef list of GLM coefficients (I usually exponentiate these values on the wrapper function that calls stargazer)
#' @param notes a vector of comments to be added at the bottom of the stargazer table
#'
#' @return raw HTML from stargazer package. FYI output is saved to file which is then included in an Rmarkdown document
#' @export
#'
#' @examples stargazer_helper(list_of_glms, title = 'test', column.labels = c("a", "b", "c"), )
stargazer_helper <- function(..., title, column.labels, coef, notes) {
  capture.output(stargazer(..., 
         title = title,
         align = T, 
         type = 'html', 
         ci=TRUE, 
         ci.level=0.95, 
         # not sure if this generates error if those coefficients are absent from model
         omit = c("pc", "study_gxe"), 
         keep.stat = "n", 
         column.labels=column.labels,
         star.cutoffs = c(0.05, 0.01, 0.001), 
         column.sep.width = '10pt', 
         coef=coef, 
         p.auto = F, 
         notes = notes))
}


#' fit_gxe
#' 
#' @description
#' fits GxE GLM for a specific exposure, SNP, and covariate set. 
#' Calculates likelihood ratio test chisq values for GxE, 2DF, and E|G associations
#' 
#' @section Warning:
#' FIGI CRC variable should always be called 'outcome' (0/1)
#'
#' @param ds dataset
#' @param exposure string containing name of exposure
#' @param snp string containing data variable name of SNP (should match names in dataset)
#' @param covariates vector of adjustment covariates
#'
#' @return a list with 2 elements - GLM model object, and a vector of chisq values (GxE, 2DF, E|G)
#' @export
#'
#' @examples fit_gxe(df = figi, exposure = 'asp_ref', snp = '6:12577203:T:C', covariates = c('age_ref_imp', 'sex', 'study_gxe'))
fit_gxe <- function(ds, exposure, snp, covariates) {
  
  # linear/numeric version of exposure to fit E|G model
  ds[, 'exposure_num'] = as.numeric(ds[, exposure])
  
  # formulas
  gxe_formula          <- paste0("outcome ~ ", snp, " * ", exposure, " + ", paste0(covariates, collapse = " + "))
  gxe_formula_base_1df <- paste0("outcome ~ ", snp, " + ", exposure, " + ", paste0(covariates, collapse = " + "))
  gxe_formula_base_2df <- paste0("outcome ~ ",             exposure, " + ", paste0(covariates, collapse = " + "))
  formula_eg           <- paste0("exposure_num ~ ", snp, " + ", paste0(covariates, collapse = " + "))
  formula_eg_base      <- paste0("exposure_num ~ ",             paste0(covariates, collapse = " + "))
  
  # fit models
  m_original <- glm(gxe_formula,          data = ds, family = 'binomial')
  m_base_1df <- glm(gxe_formula_base_1df, data = ds, family = 'binomial')
  m_base_2df <- glm(gxe_formula_base_2df, data = ds, family = 'binomial')
  m_eg       <- lm(formula_eg,            data = ds)
  m_eg_base  <- lm(formula_eg_base,       data = ds)
  
  # calculate lrtest chisq for various baseline models
  gxe_lrtest_chisq_1df <- lrtest(m_original, m_base_1df)$'Chisq'[2]
  gxe_lrtest_chisq_2df <- lrtest(m_original, m_base_2df)$'Chisq'[2]
  eg_lrtest_chisq      <- lrtest(m_eg, m_eg_base)$'Chisq'[2]
  
  # assemble results into a list
  return(list(m_original, c(gxe_lrtest_chisq_1df, gxe_lrtest_chisq_2df, eg_lrtest_chisq)))
}


#' fit_gxe_stratified
#' 
#' @description
#' generates GLM summaries of GxE interaction models, overall and stratified by group e.g. sex or tumor subsite. Outputs an HTML file that is meant to be inserted in Rmarkdown documents
#'
#' @section Warning:
#' files paths are hardcoded -- /media/work/gwis/posthoc/exposure folder
#' 
#' @param strata string describing stratifying variable. Possible choices include sex, study_design, and cancer_site_sum2
#' @param method string describing GxE methods used. Only determines which LR test p-value to report in notes section of stargazer table. Possible choices include chiSqGxE, two-step, chiSqCase, chiSq2df, chiSq3df
#' @inheritParams fit_gxe
#'
#' @return saves a raw HTML file to be used in Rmarkdown documents. Output file naming pattern is 'gxe_method_snp_exposure_stratified_strata'. 
#' @export
#'
#' @examples fit_gxe_stratified(ds = figi, exposure = 'asp_ref', snp = 'chr1_8559660_G_A', covariates = c('age_ref_imp', 'study_gxe'), strata = 'sex', method = 'chiSqGxE')
fit_gxe_stratified <- function(ds, 
                               exposure, 
                               snp, 
                               covariates, 
                               strata = c('sex', 'study_design', 'cancer_site_sum2'), 
                               method = c('chiSqGxE', 'two-step', 'chiSqCase', 'chiSq2df', 'chiSq3df')) {
  
  # limit possible argument choices
  strata <- match.arg(strata)
  method <- match.arg(method)
  covariates_nostrata <- paste0(covariates[! covariates %in% strata], collapse = " + ")
  
  # flip dosage coding if linear model parameter is negative (protective)
  # dg_model <- lm(as.formula(paste0("outcome ~ ", snp)), data = ds)
  dg_model <- lm(paste0("outcome ~ ", snp, "*", exposure, "+", paste0(covariates, collapse = "+")), data = ds)
  if(dg_model$coefficients[2] < 0) {
    # flip dosages
    ds[, paste0(snp)] <- abs(ds[, paste0(snp)] - 2)
  }
  
  # create numeric exposure and strata variables
  ds[, 'strata_num'] <- as.numeric(factor(ds[, strata]))-1
  ds[, 'exposure_num'] = as.numeric(ds[, exposure])
  
  # compile results as list  
  out <- list()
  
  ## overall GLM
  out_all <- fit_gxe(ds, exposure, snp, covariates_nostrata)
  out[['all']] <- out_all
  
  ## stratified GLM
  number_of_levels <- nlevels(factor(ds[, strata]))
  for(level in seq(number_of_levels) - 1) {
    # analysis subsets for each strata
    # need specific case for cancer_site_sum2 to capture controls
    if(strata == "cancer_site_sum2") {
      index_vector <- which(ds[, 'strata_num'] == level | ds[,'outcome'] == 0)
    } else {
      index_vector <- which(ds[, 'strata_num'] == level)
    }
    out_level <- fit_gxe(ds[index_vector,], exposure, snp, covariates_nostrata)
    out[[paste0(strata, "_", as.character(level))]] <- out_level
  }
  
  # ----------------------------------- #
  # process output, create stargazer HTML
  # ----------------------------------- #
  output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  
  # exponentiated coefficients
  list_of_glms <- lapply(out, function(x) x[[1]])
  list_of_samplesizes <- lapply(list_of_glms, function(x) paste0(c("Ca=", "Co="), rev(as.character(table(x$model$outcome))), collapse = ','))
  coefs <- lapply(list_of_glms, function(x) (exp(coef(x))))
  
  # column names for stargazer summary
  if(strata == 'sex') {
    col_label = paste0(c("All", "Female", "Male"), " (", list_of_samplesizes, ")")
  } else if(strata == 'study_design') {
    col_label = paste0(c("All", "Cohort", "Case-Control"), " (", list_of_samplesizes, ")")
  } else if(strata == 'cancer_site_sum2') {
    col_label = paste0(c("All", "Proximal", "Distal", "Rectal"), " (", list_of_samplesizes, ")")
  }
  
  # using chisq, calculate p values
  # create 'notes' vector to insert as comment in stargazer table 
  list_of_chisq <- lapply(out, function(x) x[[2]])
  
  if(method %in% c('chiSqGxE', 'two-step', 'chiSqCase')) {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[1]], df = 1, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0(col_label, ", LRtest GxE p = ", gxe_pvalues))
  } else if(method == "chiSq2df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[2]], df = 2, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0(col_label, ", LRtest 2DF p = ", gxe_pvalues))
  } else if(method == "chiSq3df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq((x[[2]] + x[[3]]), df = 3, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0(col_label, ", LRtest 3DF p = ", gxe_pvalues))
  }
  
  # call stargazer
  out_html <- stargazer_helper(list_of_glms,
                               title=paste0(gsub('\\_', '\\\\_', strata), " stratified ", gsub("\\_", "\\\\_", snp), " x ", gsub('\\_', '\\\\_', exposure)), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=notes, single.row = T)
  
  # output to file
  cat(paste(out_html, collapse = "\n"), "\n",
      file = paste0(output_dir, "gxe_", method, "_", snp, "_", exposure, "_", paste0(covariates, collapse = '_'), "_stratified_", strata, ".html"), append = F)
}


# fit_gxe_stratified(figi, 'asp_ref', 'chr1_8559660_G_A', covariates = c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'), strata = 'sex', method = 'chiSqGxE')





#' fit_gxe_covars
#' 
#' @description
#' generates GLM summaries of GxE interaction models for multiple covariate sets. (This is much easier - just apply function calls to the list of covariates). Outputs an HTML file that is meant to be inserted in Rmarkdown documents
#'
#' @section Warning:
#' files paths are hardcoded -- /media/work/gwis/posthoc/exposure folder
#' 
#' @param ds dataset
#' @param exposure string containing name of exposure
#' @param snp string containing data variable name of SNP (should match names in dataset)
#' @param covariates_list vector of adjustment covariates
#' @param method string describing GxE methods used. Only determines which LR test p-value to report in notes section of stargazer table. Possible choices include chiSqGxE, two-step, chiSqCase, chiSq2df, chiSq3df
#' @param output_dir string output directory
#' 
#' @return saves a raw HTML file to be used in Rmarkdown documents. Output file naming pattern is 'gxe_method_snp_exposure_covariate_sets'. 
#' @export
#'
#' @examples fit_gxe_covars(ds = figi, exposure = 'asp_ref', snp = 'chr1_8559660_G_A', covariates = list(c('age_ref_imp', 'study_gxe'), c('age_ref_imp', 'study_gxe', 'bmi5')), strata = 'sex', method = 'chiSqGxE')
fit_gxe_covars <- function(ds, 
                           exposure, 
                           snp, 
                           covariates_list, 
                           method = c('chiSqGxE', 'two-step', 'chiSqCase', 'chiSq2df', 'chiSq3df'),
                           output_dir) {
  
  method <- match.arg(method)
  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = "_"))
  
  # check allele frequency
  total_dosage <- sum(ds[,snp])
  total_alleles <- nrow(ds) * 2
  aaf <- total_dosage / total_alleles
  
  # ---- recode SNPs so that major allele is the reference group
  if(aaf >= 0.5) {
    # flip dosages
    ds[, snp] <- abs(ds[, paste0(snp)] - 2)
    ds[, paste0(snp, '_dose')] <- abs(ds[, paste0(snp, '_dose')] - 2)
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }

  
  # apply 'fit_gxe' over covariate_list
  out <- lapply(covariates_list, function(x) fit_gxe(ds, exposure, snp, covariates = x))
  
  # combine them to call stargazer
  list_of_glms <- lapply(out, function(x) x[[1]])
  list_of_samplesizes <- lapply(list_of_glms, function(x) paste0(c("Ca=", "Co="), rev(as.character(table(x$model$outcome))), collapse = ','))
  coefs <- lapply(list_of_glms, function(x) (exp(coef(x))))
  
  # need to calculate p values
  list_of_chisq <- lapply(out, function(x) x[[2]])
  col_label = paste0(paste0("Covariate Set ", seq(1, length(out))), " (", list_of_samplesizes, ")")
  
  if(method == "chiSqGxE" | method == "twostep" | method == "chiSqCase") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[1]], df = 1, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("Reference allele = ", ref_allele),
               paste0(col_label, ", LRtest GxE p = ", gxe_pvalues))
  } else if(method == "chiSq2df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[2]], df = 2, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("Reference allele = ", ref_allele),
               paste0(col_label, ", LRtest 2DF p = ", gxe_pvalues))
  } else if(method == "chiSq3df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq((x[[2]] + x[[3]]), df = 3, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("Reference allele = ", ref_allele),
               paste0(col_label, ", LRtest 3DF p = ", gxe_pvalues))
  }
  
  # save output html from stargazer
  out_html <- stargazer_helper(list_of_glms,
                               title=paste0(gsub("\\_", "\\\\_", snp), " x ", gsub('\\_', '\\\\_', exposure)), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=notes, single.row = T)
  
  # write object to html file
  # cat(paste(out_html, collapse = "\n"), "\n", 
  #     file = paste0(output_dir, "model_stargazer", "_", method, "_", snp, "_", exposure, "_", "covariate_sets.html"), append = F)
  cat(paste(out_html, collapse = "\n"), "\n", 
      file = paste0(output_dir, "gxe_", method, "_", snp, "_", exposure, "_covariate_sets.html"), append = F)
}


# fit_gxe_covars(figi, 'asp_ref', 'chr1_8559660_G_A', covariates_list = list(c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'), c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe', 'bmi5')), method = 'chiSqGxE')



# ---------------------------------------------------------------------------- #
# functions to generate table of p-values
# stratified by sex, study_design, and cancer_site_sum2
# for multiple covariate sets, if requested
# ----------------------------------------------------------------------------


#' calc_pval
#' 
#' calculates likelihood ratio test p value for a given exposure and SNP, for a specific method (see params for possible choices).
#'
#' @inheritParams fit_gxe
#' @param method string describing GxE methods used. Possible choices include chiSqG, chiSqGxE, chiSqGE, chiSqEDGE, chiSqCase, chiSq2df, chiSq3df
#'
#' @return a single p-value (string, scientific notation)
#' @export
#'
#' @examples calc_pval(ds = figi, exposure = 'asp_ref', snp = 'chr1_8559660_G_A', covariates = c('age_ref_imp', 'sex'), method = 'chiSqG')
calc_pval <- function(ds, exposure, snp, covariates, 
                        method = c('chiSqG', 'chiSqGxE', 'chiSqGE', 'chiSqEDGE', 'chiSqCase', 'chiSq2df', 'chiSq3df')) {
  
  method <- match.arg(method)
  covariates_formula <- sort(paste0(covariates, collapse = '+'))
  exposure_levels <- function(x) nlevels(factor(x))
  
  #-----------------#
  ## model formula
  formula_g            <- paste0("outcome ~ ", snp, " + ", exposure, " + ", covariates_formula)
  formula_g_base       <- paste0("outcome ~ "            , exposure, " + ", covariates_formula)
  formula_gxe          <- paste0("outcome ~ ", snp, " * ", exposure, " + ", covariates_formula)
  formula_gxe_base     <- paste0("outcome ~ ", snp, " + ", exposure, " + ", covariates_formula)
  formula_gxe_base_2df <- paste0("outcome ~ ",             exposure, " + ", covariates_formula)
  
  ## for E|G, categorical exposures needs to be converted to numeric
  if(exposure_levels(ds[, exposure]) <= 4) {
    ds[, "exposure_num"] <- as.numeric(factor(ds[, exposure]))-1
    formula_eg           <- paste0("exposure_num ~ ", snp, " + ", covariates_formula)
    formula_eg_base      <- paste0("exposure_num ~ ",             covariates_formula)
  } else {
    ds[, "exposure"]     <- ds[, exposure]
    formula_eg           <- paste0("exposure ~ ", snp, " + ", covariates_formula)
    formula_eg_base      <- paste0("exposure ~ ",             covariates_formula)
  }
  
  #----------------------------#
  # fit models, calculate lrtest
  if(method == 'chiSqG') {
    model_g         <- glm(formula_g,      data = ds, family = 'binomial')
    model_g_base    <- glm(formula_g_base, data = ds, family = 'binomial')
    out_pvalue      <- lrtest(model_g, model_g_base)$'Pr(>Chisq)'[2]
  } else if(method == "chiSqGxE") {
    model_gxe       <- glm(formula_gxe,      data = ds, family = 'binomial')
    model_gxe_base  <- glm(formula_gxe_base, data = ds, family = 'binomial')
    out_pvalue      <- lrtest(model_gxe, model_gxe_base)$'Pr(>Chisq)'[2]
  } else if(method == "chiSqCase") {
    model_case      <- lm(formula_eg,            data = ds[which(ds[, 'outcome'] == "1"), ])
    model_case_base <- lm(formula_eg_base,       data = ds[which(ds[, 'outcome'] == "1"), ])
    out_pvalue      <- lrtest(model_case, model_case_base)$'Pr(>Chisq)'[2]
  } else if(method == "chiSqGE") {
    model_eg            <- lm(formula_eg,            data = ds)
    model_eg_base       <- lm(formula_eg_base,       data = ds)
    out_pvalue          <- lrtest(model_eg, model_eg_base)$'Pr(>Chisq)'[2]
  } else if(method == 'chiSqEDGE') {
    model_g         <- glm(formula_g,      data = ds, family = 'binomial')
    model_g_base    <- glm(formula_g_base, data = ds, family = 'binomial')
    chisq_g         <- lrtest(model_g, model_g_base)$'Chisq'[2]
    model_eg        <- lm(formula_eg,            data = ds)
    model_eg_base   <- lm(formula_eg_base,       data = ds)
    chisq_ge        <- lrtest(model_eg, model_eg_base)$'Chisq'[2]
    out_pvalue      <- pchisq(chisq_g + chisq_ge, df = 2, lower.tail = F)
  } else if(method == "chiSq2df") {
    model_gxe           <- glm(formula_gxe,          data = ds, family = 'binomial')
    model_gxe_base_2df  <- glm(formula_gxe_base_2df, data = ds, family = 'binomial')
    out_pvalue          <- lrtest(model_gxe, model_gxe_base_2df)$'Pr(>Chisq)'[2]
  } else if(method == "chiSq3df") {
    model_gxe           <- glm(formula_gxe,          data = ds, family = 'binomial')
    model_gxe_base_2df  <- glm(formula_gxe_base_2df, data = ds, family = 'binomial')
    chisq_2df           <- lrtest(model_gxe, model_gxe_base_2df)$'Chisq'[2]
    model_eg            <- lm(formula_eg,            data = ds)
    model_eg_base       <- lm(formula_eg_base,       data = ds)
    chisq_eg            <- lrtest(model_eg, model_eg_base)$'Chisq'[2]
    out_pvalue          <- pchisq(chisq_2df + chisq_eg, df = 3, lower.tail = F)
  } else {
    out_pvalue <- NA
  }
  # return(out_pvalue)
  if(!is.na(out_pvalue)){
    return(formatC(out_pvalue, format = "e", digits = 3))
  } else {
    return(out_pvalue)
  }
}







#' pval_summary
#' 
#' @description 
#' create a table of GxE p-values, stratified by sex, study_design, and cancer_site_sum2 and using multiple covariate sets if applicable. Provides broad overview of potentially significant interactions for further follow up
#' 
#' @section Warning: 
#' only stratified by the variables listed above. For exposures such as HRT, I explicitly exclude sex stratification 
#'
#' @param ds dataset
#' @param exposure string containing name of exposure
#' @param snp string containing data variable name of SNP (should match names in dataset)
#' @param covariates_list list of vectors of adjustment covariate sets
#' @param method string describing GxE methods used. Possible choices include chiSqG, chiSqGxE, chiSqGE, chiSqEDGE, chiSqCase, chiSq2df, chiSq3df
#' @param output_dir string output directory
#'
#' @return a dataframe of counts and pvalues by strata and covariate set
#' @export
#'
#' @examples pval_summary(figi, exposure = 'asp_ref', snp = 'chr1_8559660_G_A', covariates_list = list(c('age', 'sex'), c('age', 'sex', 'bmi')), method = 'chiSqGxE')
pval_summary <- function(ds, 
                         exposure, 
                         snp, 
                         covariates_list, 
                         method = c('chiSqG', 'chiSqGxE', 'chiSqGE', 'chiSqEDGE', 'chiSqCase', 'chiSq2df', 'chiSq3df'),
                         output_dir) {

  method <- match.arg(method)
  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  
  # apply helper function to list of covariate set vectors
  helper <- function(covars, analysis = c("counts", "pvalues")) {
    analysis <- match.arg(analysis)
    
    # generate temporary complete case data based on covariate set
    ds_covar_subset <- ds %>% 
      filter(complete.cases(.[, covars]))
    
    # data subsets (vectors of booleans)
    subset_all      <- rep(T, nrow(ds_covar_subset))
    subset_female   <- ds_covar_subset$sex == 0
    subset_male     <- ds_covar_subset$sex == 1
    subset_cohort   <- ds_covar_subset$study_design == "Cohort"
    subset_cc       <- ds_covar_subset$study_design == "Case-Control"
    subset_proximal <- ds_covar_subset$cancer_site_sum2 == 'proximal' | ds_covar_subset$outcome == 0
    subset_distal   <- ds_covar_subset$cancer_site_sum2 == 'distal'   | ds_covar_subset$outcome == 0
    subset_rectal   <- ds_covar_subset$cancer_site_sum2 == 'rectal'   | ds_covar_subset$outcome == 0
    
    # remove sex stratification if variables are HRT related
    if(exposure %in% c("hrt_ref_pm2", "eo_ref_pm_gxe", "ep_ref_pm_gxe")) {
      subset_list <- list(subset_all, subset_cohort, subset_cc, subset_proximal, subset_distal, subset_rectal)
      covars <- paste0(covars[! covars %in% 'sex'], collapse = " + ")
    } else {
      subset_list <- list(subset_all, subset_female, subset_male, subset_cohort, subset_cc, subset_proximal, subset_distal, subset_rectal)
    }
    
    # get counts or pvalues
    if(analysis == "counts") {
      tally_wrapper <- function(x) { 
        counts <- ds_covar_subset[x, ] %>% 
          group_by(outcome) %>% 
          tally()
        paste0(counts[2,2], "/", counts[1,2])
      }
      return(lapply(subset_list, tally_wrapper))
    } else if(analysis == "pvalues") {
      return(lapply(subset_list, function(x) calc_pval(ds_covar_subset[x,], exposure, snp, covars, method)))
    }
  }
  
  # call helper function
  counts <- sapply(covariates_list, helper, analysis = 'counts')
  pvalues <- sapply(covariates_list, helper, analysis = 'pvalues')
  
  # column labels for table
  if(exposure %in% c("hrt_ref_pm2", "eo_ref_pm_gxe", "ep_ref_pm_gxe")) {
    col_labels <- c('All', 'Cohort', 'Case-Control', 'Proximal', 'Distal', 'Rectal')
  } else {
    col_labels <- c('All', 'Females', 'Males', 'Cohort', 'Case-Control', 'Proximal', 'Distal', 'Rectal')
  }

  # combine outputs in the right order
  # (useful when there are multiple covariate sets)
  xx <- (2*(seq(ncol(counts)) - 1) + 1) + 1
  yy <- (2*seq(ncol(pvalues))) + 1
  neworder <- order(c(1, xx, yy))
  out <- data.frame(cbind(col_labels, counts, pvalues)[, neworder])
  
  # fix names
  z1 <- paste0("Set", seq(1, length(covariates_list)))
  z2 <- c('Ca/Co', 'p-value')
  z3 <- expand.grid(z1, z2)
  z4 <- z3[order(z3$Var1),]
  labels2 <- apply(z4, 1 , paste, collapse=" ")
  names(out) <- c("Groups", labels2)
  
  # output
  # return(out)
  saveRDS(out, file = paste0(output_dir, "pvalues_dataframe_", method, "_", snp, "_", exposure, ".rds"))
}



#' posthoc_input
#' 
#' put together input data for posthoc analysis. a note about SNPs - makes it so that data contains chrX_BP_REF_ALT (dose), in addition to chrX_BP_REF_ALT_dose/p0/p1/p2. This function is for very specific use case, edit as needed! Also, paths are hardcoded, should be ok because I'm not changing it for a while now
#'
#' @param exposure 
#' @param hrc_version 
#'
#' @return
#' @export
#'
#' @examples
posthoc_input <- function(exposure, hrc_version, dosage_filename) {
  # exposure subset as determined by analysis plan
  exposure_subset <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version , "_gxeset_", exposure, "_basic_covars_glm.rds"))[,'vcfid']
  
  # binary dosage output
  # geno <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_output_", exposure, ".rds")) 
  
  geno <- readRDS(glue("/media/work/gwis/posthoc/{exposure}/{dosage_filename}"))
  
  geno_dose <- geno %>%
    dplyr::select(-contains(c('p0', 'p1', 'p2')))
  
  # clean up SNP names from binarydosage output
  # removes '_dose' from SNP name
  snpname_clean <- function(x) {
    tmp <- gsub("\\.", "\\_", x)
    tmp <- gsub("X", "chr", tmp)
    tmp <- gsub("\\_dose", "", tmp)
    return(tmp)
  }
  
  # cleans SNPs but keeps original suffixes ("_dose", "p0/p1/p2")
  snpname_clean2 <- function(x) {
    tmp <- gsub("\\.", "\\_", x)
    tmp <- gsub("X", "chr", tmp)
    return(tmp)
  }
  
  names(geno_dose) <- snpname_clean(names(geno_dose))
  names(geno) <- snpname_clean2(names(geno))
  
  # dataset that contains all gxe samples and variables (not specific to any exposure)
  out <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_analysis_data_glm.rds")) %>%
    filter(vcfid %in% exposure_subset) %>% 
    inner_join(geno_dose, 'vcfid') %>%
    inner_join(geno, 'vcfid') %>% 
    mutate(cancer_site_sum2 = factor(cancer_site_sum2, levels = c("proximal", "distal", "rectal")))
  
  return(out)
  
}


