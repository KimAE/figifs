#=============================================================================#
# functions to perform some posthoc analyses
# stratified models (stargazer)
# stratify by covariate sets
# can stratify by sex, study design, tumor site
#=============================================================================#





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
fit_gxe_covars <- function(data, 
                           exposure, 
                           snp, 
                           covariates_list, 
                           method = c('chiSqGxE', 'two-step', 'chiSqCase', 'chiSq2df', 'chiSq3df'),
                           path) {
  
  method <- match.arg(method)
  
  
  # convert 12:50610976:C:T to chr12_50610976_C_T_dose
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  
  # SNP - need to change the SNP name to match the file (because of binarydosage output)
  tmp <- unlist(strsplit(snpfix_short, "_"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]
  
  # we probably don't want to flip SNPs for each covariates set, so let's base direction on the first set (always the simpler model), in the interaction model. 
  model_check <- glm(glue("outcome ~ {exposure}*{snpfix} + {glue_collapse(covariates_list[[1]], sep = '+')}"), family = 'binomial', data = data)
  
  if (model_check[[1]][snpfix] < 0) {
    snp_new <- glue("{chr}_{bp}_{alt}_{ref}_dose_flipped")
    data[[snp_new]] <- abs(2-data[, snpfix])
    ref_allele = alt
  } else {
    snp_new <- snpfix
    ref_allele = ref
  }
  
  # apply 'fit_gxe' over covariate_list
  out <- lapply(covariates_list, function(x) fit_gxe(data, exposure, snp_new, covariates = x))
  
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
                               title=paste0(gsub("\\_", "\\\\_", snp_new), " x ", gsub('\\_', '\\\\_', exposure)), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=notes, single.row = T)
  
  # write object to html file
  cat(paste(out_html, collapse = "\n"), "\n", 
      file = glue("{path}/gxe_{method}_{snp}_{exposure}_covariate_sets.html"), append = F)
}