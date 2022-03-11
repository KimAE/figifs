#=============================================================================#
# functions to perform some posthoc analyses
#=============================================================================#

#' posthoc_sig_wrapper
#'
#' @param filename 
#' @param output_dir 
#'
#' @return
#' @export
#'
#' @examples
posthoc_sig_wrapper <- function(filename, output_dir) {
  
  # check if data.frame is empty (no rows)
  tmp <- readRDS(glue(output_dir, filename)) %>% 
    mutate(method = gsub(".rds", "", filename))
  
  if (grepl("manhattan", filename)) {
    out <- tmp %>% 
      dplyr::rename(Pval = P) %>% 
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval, method) %>% 
      dplyr::arrange(Chromosome, Location)
  } else {
    out <- tmp %>% 
      dplyr::rename(Pval = step2p) %>% 
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval, method) %>% 
      dplyr::arrange(Chromosome, Location)
  }
}



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# some posthoc functions
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# i'm gonna read in the epi file (gxeset) in the targets file
# then read a SNPs at the time you run some posthoc stuff.. maybe that's easier? 

# so you can create the dosage data.frames, save it, and then _targets should know not to run it again?

# figi <- qread("../data/FIGI_v3.0_gxeset_analysis_data_glm.qs")
# exposure_subset <- readRDS(glue("/scratch/andreeki/gwis/gxescan/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))[, 'vcfid']




# function to get SNP, then merge it with the working dataset
#' get_dosage
#'
#' @param exposure_subset_vcfids 
#' @param exposure 
#' @param hrc_version 
#' @param snp 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
get_dosage <- function(exposure_subset_vcfids, exposure, hrc_version, snp, path) {
  
  # data_subset <- data %>% 
  #   dplyr::filter(vcfid %in% exposure_subset)
  
  # wdir = "output/posthoc"
  wdir = glue("{path}/posthoc")
  dir.create(file.path(glue("{wdir}")), showWarnings = F)
  
  tmp <- unlist(strsplit(snp, ":"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]
  
  snpid          <- glue("chr{chr}:{bp}")
  snpid_filename <- glue("chr{chr}_{bp}")
  
  # get dosages
  bdose <- readRDS(paste0("/project/dconti_250/HRC_BDose/FIGI_snpid_fix_chr", chr, ".rds"))
  dosages <- GetSNPValues(bdose, snp, geneProb = T)
 
  # Don't filter by exposure subset vcfids (just get all) 
  dosages_df <- data.frame(dosages) %>% 
    dplyr::mutate(vcfid = rownames(.)) #%>% 
    #dplyr::filter(vcfid %in% exposure_subset_vcfids) 
    # inner_join(data_subset, 'vcfid')
  
  snpname_clean <- function(x) {
    tmp <- gsub("\\.", "\\_", x)
    tmp <- gsub("X", "chr", tmp)
    return(tmp)
  }
  
  names(dosages_df) <- snpname_clean(names(dosages_df))
  
  # save the file, return the filename for targets
  qsave(dosages_df, file = glue("{wdir}/dosage_{snpid_filename}.qs"))
  return(glue("{wdir}/dosage_{snpid_filename}.qs"))
  
  
}







#------------------------------------------------------------------------------#
# stratified gxe analysis (pooled)
#------------------------------------------------------------------------------#

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
	 report = 'vcs*p',
         notes = notes))
}




#' fit_gxe_covars
#'
#' @param data_epi 
#' @param exposure 
#' @param snp 
#' @param covariates_list 
#' @param method 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
fit_gxe_covars <- function(data_epi, 
                           exposure, 
                           snp, 
                           covariates_list, 
                           method = c('chiSqGxE', 'two-step', 'chiSqCase', 'chiSq2df', 'chiSq3df'),
                           path, 
			   flip_allele = F) {
  
  method <- match.arg(method)
  
  wdir = glue("{path}/posthoc")
  
  snp_info <- unlist(strsplit(snp, split = ":"))
  
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  # check if SNP has to be recoded 
  #model_check <- glm(glue("outcome ~ {exposure}*{snpfix} + {glue_collapse(covariates_list[[1]], sep = '+')}"), family = 'binomial', data = data)
  #
  #if (model_check[[1]][snpfix] < 0) {
  #  snp_old <- snpfix
  #  snp_tmp <- unlist(strsplit(snpfix, split = "_"))
  #  chr <- snp_tmp[1]
  #  bp <- snp_tmp[2]
  #  a1 <- snp_tmp[3]
  #  a2 <- snp_tmp[4]
  #  snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
  #  data[[snp_new]] <- abs(2-data[, snp_old])
  #} else {
  #  snp_new <- snpfix
  #}
  
  if (flip_allele == T) {
    snp_old <- snpfix
    snp_tmp <- unlist(strsplit(snpfix, split = "_"))
    chr <- snp_tmp[1]
    bp <- snp_tmp[2]
    a1 <- snp_tmp[3]
    a2 <- snp_tmp[4]
    snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
    data[[snp_new]] <- abs(2-data[, snp_old])
  } else {
    snp_new <- snpfix
  }

  # apply 'fit_gxe' over covariate_list
  out <- lapply(covariates_list, function(x) fit_gxe(data, exposure, snp_new, covariates = x, method = method))
  
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
               paste0("GxE term LRtest p = ", gxe_pvalues))
  } else if(method == "chiSq2df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[2]], df = 2, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("2DF LRtest p = ", gxe_pvalues))
  } else if(method == "chiSq3df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq((x[[2]] + x[[3]]), df = 3, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("3DF LRtest p = ", gxe_pvalues))
  }
  
  # save output html from stargazer
  out_html <- stargazer_helper(list_of_glms,
                               title=paste0(gsub("\\_", "\\\\_", snp_new), " x ", gsub('\\_', '\\\\_', exposure)), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=notes, single.row = T)
  
  # write object to html file
  cat(paste(out_html, collapse = "\n"), "\n", 
      file = glue("{wdir}/gxe_models_{exposure}_{hrc_version}_{snpfix}_covariate_sets.html"), append = F)
  
  # return(glue("{wdir}/gxe_models_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}_covariate_sets.html"))

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
#' @param method string statistical method (e.g. chiSqGxE)
#'
#' @return a list with 2 elements - GLM model object, and a vector of chisq values (GxE, 2DF, E|G)
#' @export
#'
#' @examples fit_gxe(df = figi, exposure = 'asp_ref', snp = '6:12577203:T:C', covariates = c('age_ref_imp', 'sex', 'study_gxe'), method = 'chiSqGxE')
fit_gxe <- function(ds, exposure, snp, covariates, method) {
  
  # create linear/numeric version of exposure to fit E|G model
  ds[, 'exposure_num'] = as.numeric(ds[, exposure])
  
  # formulas
  gxe_formula          <- paste0("outcome ~ ", snp, " * ", exposure, " + ", paste0(covariates, collapse = " + "))
  gxe_formula_base_1df <- paste0("outcome ~ ", snp, " + ", exposure, " + ", paste0(covariates, collapse = " + "))
  gxe_formula_base_2df <- paste0("outcome ~ ",             exposure, " + ", paste0(covariates, collapse = " + "))
  formula_eg           <- paste0("exposure_num ~ ", snp, " + ", paste0(covariates, collapse = " + "))
  formula_eg_base      <- paste0("exposure_num ~ ",             paste0(covariates, collapse = " + "))
  
  # fit models
  if (method == "chiSqGxE") {
    glm_model      <- glm(gxe_formula,          data = ds, family = 'binomial')
    glm_model_base <- glm(gxe_formula_base_1df, data = ds, family = 'binomial')
    lrtest_pvalue  <- lrtest(glm_model, glm_model_base)$'Chisq'[2]

  } else if (method == 'chiSq2df') {
    # report GxE 1df for 2df models
    glm_model          <- glm(gxe_formula,          data = ds, family = 'binomial')
    glm_model_2df      <- glm(gxe_formula_base_2df, data = ds, family = 'binomial')
    lrtest_pvalue      <- lrtest(glm_model, glm_model_2df)$'Chisq'[2]

  } else if (method == 'chiSq3df') {
    # report GxE 1df for 3df models
    glm_model          <- glm(gxe_formula,          data = ds, family = 'binomial')
    glm_model_2df      <- glm(gxe_formula_base_2df, data = ds, family = 'binomial')
    m_eg               <- lm(formula_eg,            data = ds)
    m_eg_base          <- lm(formula_eg_base,       data = ds)
    pval_2df           <- lrtest(glm_model, glm_model_2df)$'Chisq'[2]
    pval_eg            <- lrtest(m_eg, m_eg_base)$'Chisq'[2]
    lrtest_pvalue      <- pval_2df + pval_eg
    
  } else if (method == "chiSqCase") {
    # care - outcome variable 'hard-coded'
    glm_model          <- lm(formula_eg,            data = ds %>% filter(outcome == 1))
    m_eg_base          <- lm(formula_eg_base,       data = ds %>% filter(outcome == 1))
    lrtest_pvalue      <- lrtest(glm_model, m_eg_base)$'Chisq'[2]
    
  } else if (method == "chiSqGE") {
    glm_model          <- lm(formula_eg,            data = ds)
    m_eg_base          <- lm(formula_eg_base,       data = ds)
    lrtest_pvalue      <- lrtest(glm_model, m_eg_base)$'Chisq'[2]
    
  }
  
  # assemble results into a list
  return(list(glm_model, lrtest_pvalue))
}


#' fit_gxe_stratified
#' 
#' @description
#' generates GLM summaries of GxE interaction models, overall and stratified by group e.g. sex or tumor subsite. Outputs an HTML file that is meant to be inserted in Rmarkdown documents
#'
#' @param data_epi dataset
#' @param exposure exposure
#' @param snp snp
#' @param covariates covars
#' @param strata string describing stratifying variable. Possible choices include sex, study_design, cancer_site_sum2, bmic3 (add as needed)
#' @param method. determines which GLM model is reported + LR test p-value. For GxE, Two-step, 2/3DF, I simply output GxE 1DF results
#' @param flip_allele flip or no
#' @param path path
#'
#' @return saves a raw HTML file to be used in Rmarkdown documents. Output file naming pattern is 'gxe_method_snp_exposure_stratified_strata'. 
#' @export
#'
#' @examples fit_gxe_stratified(ds = figi, exposure = 'asp_ref', snp = 'chr1_8559660_G_A', covariates = c('age_ref_imp', 'study_gxe'), strata = 'sex', method = 'chiSqGxE')
fit_gxe_stratified <- function(data_epi,
                               exposure,
                               snp,
                               covariates,
                               strata = c('sex', 'study_design', 'cancer_site_sum2', 'bmic3'),
                               method = c('chiSqGxE', 'chiSqCase', 'chiSq2df', 'chiSq3df', 'chiSqGE'),
                               flip_allele = F,
                               path) {
  
  # make path the exposure directory e.g. /media/work/gwis_test/exposure
  wdir = glue("{path}/output/posthoc")
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = ":"))
  
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  
  data_dose <-
    qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  if (flip_allele == T) {
    snp_old <- snpfix
    snp_tmp <- unlist(strsplit(snpfix, split = "_"))
    chr <- snp_tmp[1]
    bp <- snp_tmp[2]
    a1 <- snp_tmp[3]
    a2 <- snp_tmp[4]
    snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
    data[[snp_new]] <- abs(2 - data[, snp_old])
  } else {
    snp_new <- snpfix
  }
  
  # limit possible argument choices
  strata <- match.arg(strata)
  method <- match.arg(method)
  covariates_nostrata <-
    paste0(covariates[!covariates %in% strata], collapse = " + ")
  
  # create numeric exposure and strata variables
  data[, 'strata_num'] <- as.numeric(factor(data[, strata])) - 1
  data[, 'exposure_num'] = as.numeric(data[, exposure])
  
  # compile results as list
  out <- list()
  
  ## overall GLM
  ## fit_gxe defined elsewhere
  out_all <- fit_gxe(data, exposure, snp_new, covariates, method)
  out[['all']] <- out_all
  
  ## stratified GLM
  number_of_levels <- nlevels(factor(data[, strata]))
  
  for (level in seq(number_of_levels) - 1) {
    # analysis subsets for each strata
    # cancer_site_sum2 - get specific tumor site + controls 
    if (strata == "cancer_site_sum2") {
      index_vector <-
        which(data[, 'strata_num'] == level |
                data[, 'outcome'] == 0)
    } else {
      index_vector <- which(data[, 'strata_num'] == level)
    }
    out_level <-
      fit_gxe(data[index_vector,], exposure, snp_new, covariates_nostrata, method)
    out[[paste0(strata, "_", as.character(level))]] <- out_level
  }
  
  
  # ----------------------------------- #
  # process output, create stargazer HTML
  # ----------------------------------- #
  
  # exponentiated coefficients
  list_of_glms <- lapply(out, function(x)
    x[[1]])
  list_of_samplesizes <-
    lapply(list_of_glms, function(x)
      paste0(c("Ca=", "Co="), rev(as.character(
        table(x$model$outcome)
      )), collapse = ','))
  coefs <- lapply(list_of_glms, function(x)
    (exp(coef(x)))) # need this for stargazer call
  
  
  # column names for stargazer summary
  if (strata == 'sex') {
    col_label = paste0(c("All", "Female", "Male"), " (", list_of_samplesizes, ")")
  } else if (strata == 'study_design') {
    col_label = paste0(c("All", "Cohort", "Case-Control"),
                       " (",
                       list_of_samplesizes,
                       ")")
  } else if (strata == 'cancer_site_sum2') {
    col_label = paste0(c("All", "Proximal", "Distal", "Rectal"),
                       " (",
                       list_of_samplesizes,
                       ")")
  } else if (strata == 'bmic3') {
    col_label = paste0(c("All", "Normal", "Overweight", "Obese"), 
                       " (", list_of_samplesizes, ")")
  }
  
  
  
  # using chisq, calculate p values
  # create 'notes' vector to insert as comment in stargazer table
  list_of_chisq <- lapply(out, function(x)
    x[[2]])
  
  if (method %in% c('chiSqGxE', 'chiSqCase', 'chiSqGE')) {
    gxe_pvalues <-
      do.call(c, lapply(list_of_chisq, function(x)
        formatC(
          pchisq(x, df = 1, lower.tail = F),
          format = "e",
          digits = 3
        )))
    notes <- c(
      "(PC and Study estimates omitted from table)",
      paste0("LRtest p = ", gxe_pvalues)
    )
  } else if (method == "chiSq2df") {
    gxe_pvalues <-
      do.call(c, lapply(list_of_chisq, function(x)
        formatC(
          pchisq(x, df = 2, lower.tail = F),
          format = "e",
          digits = 3
        )))
    notes <- c("(PC and Study estimates omitted from table)",
               paste0("2DF LRtest p = ", gxe_pvalues))
  } else if (method == "chiSq3df") {
    gxe_pvalues <-
      do.call(c, lapply(list_of_chisq, function(x)
        formatC(
          pchisq(x, df = 3, lower.tail = F),
          format = "e",
          digits = 3
        )))
    notes <- c("(PC and Study estimates omitted from table)",
               paste0("3DF LRtest p = ", gxe_pvalues))
  }
  
  
  # call stargazer
  out_html <- stargazer_helper(
    list_of_glms,
    title = paste0(
      gsub('\\_', '\\\\_', strata),
      " stratified ",
      gsub("\\_", "\\\\_", snp_new),
      " x ",
      gsub('\\_', '\\\\_', exposure)
    ),
    column.labels = col_label,
    coef = coefs,
    notes = notes,
    single.row = T
  )
  
  # output to file

  filename = glue("{wdir}/gxe_models_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}_{deparse(substitute(data_epi))}_stratified_by_{strata}.html")

  cat(
    paste(out_html, collapse = "\n"),
    "\n",
    file = filename,
    append = F)
  
  return(filename)
}








# create function that takes the bed files from functional annotation plots and creates a table summarizing overlaps between scacheri and finucane 

# the function should simply define input and output files and then call the perl script in ../data/Annotation etc 
# just SNP should suffice (chr5_12345)
# snp = "5:40252294:C:T"
#

# you kno what, i'm going to put this at the end of the functional annotation plot scirpt.. 
# NO ... since you never know if you'll expand upon these annotations.. 
#' output_chromatin_mark_overlap
#'
#' @param snp 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
output_chromatin_mark_overlap <- function(snp, path) {
  
  idir <- glue("{path}/functional_plot")
  odir <- glue("{path}/functional_annotation")
  
  dir.create(file.path(glue("{odir}")), showWarnings = F)
  
  tmp <- unlist(strsplit(snp, ":"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]
  
  file_input <- (glue("{idir}/functional_annotation_chr{chr}_{bp}.bed"))
  file_output <- (glue("{odir}/functional_annotation_chr{chr}_{bp}_chromatin_marks.tsv"))
  
  system(glue("perl ../data/Annotation_Workflow/chromatin_marks/summarize-overlap-regulatory-region.pl {file_input} {file_output}"))
  return(file_output)
}











#' interactionR_table2
#'
#' @param obj 
#' @param file_path 
#'
#' @return
#' @export
#'
interactionR_table2 <- function (obj, pvalue, file_path = NA) 
{
  if (class(obj) != "interactionR") {
    stop("Argument 'obj' must be an object of class 'interactionR',\n             use the interactionR() function to generate such object ")
  }
  beta1 <- obj$exp_names[1]
  beta2 <- obj$exp_names[2]
  em <- obj$analysis
  d <- obj$dframe
  d$Estimates <- as.character(round(d$Estimates, 2))
  d$CI.ll <- as.character(round(d$CI.ll, 2))
  d$CI.ul <- as.character(round(d$CI.ul, 2))
  E1.absent <- paste(beta1, "absent", sep = " ")
  E1.present <- paste(beta1, "present", sep = " ")
  E2.absent <- paste(beta2, "absent", sep = " ")
  E2.present <- paste(beta2, "present", sep = " ")
  WithinStrataEffect1 <- paste("Effect of", beta2, "within the strata of", 
                               beta1, sep = " ")
  WithinStrataEffect2 <- paste("Effect of", beta1, "within the strata of", 
                               beta2, sep = " ")
  if (grepl("\\blog\\b", obj$call[3]) || grepl("poisson", obj$call[3])) {
    effect_measure <- "RR [95% CI]"
  }
  else {
    effect_measure <- "OR [95% CI]"
  }
  if (em) {
    t <- data.frame(c(NA, NA, E1.absent, E1.present, "Multiplicative scale", 
                      "RERI"), c(NA, effect_measure, NA, NA, NA, NA), c(NA, 
                                                                        effect_measure, NA, NA, NA, NA), c(NA, effect_measure, 
                                                                                                           NA, NA, NA, NA), stringsAsFactors = FALSE)
    names(t) <- c("*", E2.absent, E2.present, WithinStrataEffect1)
    t[3, 2] <- paste("1", "[Reference]", sep = " ")
    t[3, 3] <- paste(d[2, 2], " [", d[2, 3], ", ", d[2, 4], 
                     "]", sep = "")
    t[3, 4] <- paste(d[5, 2], " [", d[5, 3], ", ", d[5, 4], 
                     "]", sep = "")
    t[4, 2] <- paste(d[3, 2], " [", d[3, 3], ", ", d[3, 4], 
                     "]", sep = "")
    t[4, 3] <- paste(d[4, 2], " [", d[4, 3], ", ", d[4, 4], 
                     "]", sep = "")
    t[4, 4] <- paste(d[6, 2], " [", d[6, 3], ", ", d[6, 4], 
                     "]", sep = "")
    t[5, 2] <- paste(d[7, 2], " [", d[7, 3], ", ", d[7, 4], 
                     "]", sep = "")
    t[6, 2] <- paste(d[8, 2], " [", d[8, 3], ", ", d[8, 4], 
                     "]", sep = "")
    t2 <- flextable(t)
    t2 <- set_caption(t2, paste("Modification of the effect of", 
                                beta1, "and", beta2, sep = " "))
  }
  else {
    t <- data.frame(c(NA, NA, E1.absent, E1.present, WithinStrataEffect2, 
                      "Multiplicative scale", "RERI", "AP", "SI"), c(NA, 
                                                                     effect_measure, NA, NA, NA, NA, NA, NA, NA), c(NA, 
                                                                                                                    effect_measure, NA, NA, NA, NA, NA, NA, NA), c(NA, 
                                                                                                                                                                   effect_measure, NA, NA, NA, NA, NA, NA, NA), stringsAsFactors = FALSE)
    names(t) <- c("*", E2.absent, E2.present, WithinStrataEffect1)
    t[3, 2] <- paste("1", "[Reference]", sep = " ")
    t[3, 3] <- paste(d[2, 2], " [", d[2, 3], ", ", d[2, 4], 
                     "]", sep = "")
    t[3, 4] <- paste(d[5, 2], " [", d[5, 3], ", ", d[5, 4], 
                     "]", sep = "")
    t[4, 2] <- paste(d[3, 2], " [", d[3, 3], ", ", d[3, 4], 
                     "]", sep = "")
    t[4, 3] <- paste(d[4, 2], " [", d[4, 3], ", ", d[4, 4], 
                     "]", sep = "")
    t[4, 4] <- paste(d[6, 2], " [", d[6, 3], ", ", d[6, 4], 
                     "]", sep = "")
    t[5, 2] <- paste(d[7, 2], " [", d[7, 3], ", ", d[7, 4], 
                     "]", sep = "")
    t[5, 3] <- paste(d[8, 2], " [", d[8, 3], ", ", d[8, 4], 
                     "]", sep = "")
    t[6, 2] <- paste(d[9, 2], " [", d[9, 3], ", ", d[9, 4], 
                     "]", sep = "")
    t[7, 2] <- paste(d[10, 2], " [", d[10, 3], ", ", d[10, 
                                                       4], "]", sep = "")
    t[8, 2] <- paste(d[11, 2], " [", d[11, 3], ", ", d[11, 
                                                       4], "]", sep = "")
    t[9, 2] <- paste(d[12, 2], " [", d[12, 3], ", ", d[12, 
                                                       4], "]", sep = "")
    t[7, 3] <- paste("p=", pvalue, sep = "")
    t2 <- flextable(t)
    t2 <- fontsize(t2, size = 15)
    t2 <- fontsize(t2, size = 15, part = 'header')
    t2 <- set_caption(t2, paste("Interaction of", beta1, 
                                "and", beta2, sep = " "))
  }
  
  print(t2)
  invisible(t2)
  
}












#' reri_wrapper
#'
#' @param data_epi 
#' @param exposure 
#' @param snp 
#' @param covariates 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
reri_wrapper <- function(data_epi, exposure, snp, covariates, path){
  
  wdir <- glue("{path}/posthoc")
  
  snp_info <- unlist(strsplit(snp, split = ":"))
  
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  # check if SNP has to be recoded
  model_check <- glm(glue("outcome ~ {exposure}*{snpfix} + {glue_collapse(covariates, sep = '+')}"), family = 'binomial', data = data)
  
  if (model_check[[1]][snpfix] < 0) {
    snp_old <- snpfix
    snp_tmp <- unlist(strsplit(snpfix, split = "_"))
    chr <- snp_tmp[1]
    bp <- snp_tmp[2]
    a1 <- snp_tmp[3]
    a2 <- snp_tmp[4]
    snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
    data[[snp_new]] <- abs(2-data[, snp_old])
  } else {
    snp_new <- snpfix
  }
  
  model <- glm(glue("outcome ~ {exposure}*{snp_new} + {glue_collapse(covariates, sep = '+')}"), family = binomial(link = "logit"), data = data)
  # summary(model)
  
  
  # calculate p value
  
  ## (get coef positions..) 
  coef_names <- names(coef(model))
  coef_exposure <- grep(exposure, coef_names)[1]
  coef_snp <- grep(snp_new, coef_names)[1]
  coef_interaction <- grep(exposure, coef_names)[2]
  
  ## calculation
  reri_est = epi.interaction(model = model, coef = c(coef_exposure,coef_snp,coef_interaction), param = 'product', type  = 'RERI', conf.level = 0.95)
  coef_keep <- coef_names[c(coef_exposure, coef_snp, coef_interaction)]
  cov.mat <- vcov(model)
  V2 = cov.mat[coef_keep, coef_keep]
  
  reri_se = deltamethod( ~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1, 
                         mean = c( coef(model)[coef_exposure], coef(model)[coef_snp], coef(model)[coef_interaction]), 
                         cov = V2)
  
  reri_pval = format.pval(2*pnorm(-abs(reri_est[1, 1] / reri_se)), digits = 4)
  
  # output 
  value = interactionR(model, exposure_names = c(exposure, snp_new), ci.type = "delta", ci.level = 0.95, em = F, recode = F)
  out <- interactionR_table2(value, pvalue = reri_pval) # just save as RDS and use flextable to print in rmarkdown docs.. 
  
  saveRDS(out, file = glue("{wdir}/reri_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  
  
}




#------------------------------------------------------------------------------#
# create AAF plot by study
#------------------------------------------------------------------------------#



# snp = chr:bp:ref:alt
# path = '/media/work/gwis_test/exposure/
#' create_aaf_study_plot
#'
#' @param data 
#' @param exposure 
#' @param hrc_version 
#' @param snp 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
create_aaf_study_plot <- function(data, exposure, hrc_version, snp, path) {

  # use dosage to calculate AAF
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



  tmp <- qread(glue("{path}/output/posthoc/dosage_{chr}_{bp}.qs")) %>%
    inner_join(data, 'vcfid')

  aaf <- function(x) {
    sum(x) / nrow(x)
  }

  out <- tmp %>%
    group_by(study_gxe) %>%
    summarise(total = n(),
              study_aaf = sum(.data[[snpfix]]) / (total*2)) %>%
    arrange(study_aaf) %>%
    mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

  ggplot(aes(x = study_gxe, y = study_aaf), data = out) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270)) +
    xlab("Study") +
    ylab("Alternate Allele Frequency")
  ggsave(glue("{path}/output/posthoc/aaf_by_studygxe_{exposure}_{hrc_version}_{snpfix}.png"), height = 4, width = 6)

}

