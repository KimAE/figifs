#-----------------------------------------------------------------------------#
# Regression functions ------
#
#-----------------------------------------------------------------------------#

# plot stratified analysis Gxe

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

# model1 <- glm(outcome ~ X5.40252294*aspirin+age_ref_imp+sex+study_gxe+PC1+PC2+PC3, data = posthoc_df, family = 'binomial')
# create_glm_stratified_plot(model1, "X5.40252294", "aspirin", gxe = T)
#
# model1 <- glm(outcome ~ aspirin*X5.40252294+age_ref_imp+sex+study_gxe+PC1+PC2+PC3, data = posthoc_df, family = 'binomial')
# create_glm_stratified_plot(model1, "aspirin", "X5.40252294", gxe = F)
#




# strat_wrapper <- function(data, snp, exposure) {
#
#   # get counts
#   ctrl <- data %>%  filter(type == 0)
#   ta <- data.frame(table(ctrl[, snp], ctrl[, exposure]))
#
#   case <- data %>% filter(type == 1)
#   tb <- data.frame(table(case[, snp], case[, exposure]))
#
#   tf <- bind_cols(ta, tb[, 3, drop = F]) %>%
#     mutate(counts = paste(Freq, Freq1, sep = "/")) %>%
#     dplyr::select(-starts_with('Freq')) %>%
#     spread(Var2, counts)
#
#   # get E estimates stratified by genotype (0,1,2)
#   covarsE <- c(covars_full, exposure)
#   form <- as.formula(paste0('type ~ ', paste(covarsE, collapse = "+")))
#
#   # just make sure ORs are correct..
#   b <-  data[, c(covarsE, 'type', snp)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(snp) %>%
#     do(confint_tidy(glm(form, data = ., family = 'binomial')))
#
#   a <- data[, c(covarsE, 'type', snp)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(snp) %>%
#     do(tidy(glm(form, data = ., family = 'binomial'))) %>%
#     bind_cols(b[, c("conf.low", "conf.high")]) %>%
#     mutate(OR = paste0(round(exp(estimate), 2), " (", round(exp(conf.low), 2), "-", round(exp(conf.high), 2), ")")) %>%
#     filter(grepl(exposure, term)) %>% dplyr::select(!!snp, OR, term) %>%
#     spread_(snp, value = 'OR')
#
#   af <- rbind(c("baseline", "1.0 (Ref)", "1.0 (Ref)", "1.0 (Ref)"), a)
#
#   # get G estimate (additive) for each E level
#   covarsG <- c(covars_full, snp)
#   form <- as.formula(paste0('type ~ ', paste(covarsG, collapse = "+")))
#
#   b <-  data[, c(covarsG, 'type', exposure)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(exposure) %>%
#     do(confint_tidy(glm(form, data = ., family = 'binomial')))
#
#   aa <- data[, c(covarsG, 'type', exposure)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(exposure) %>%
#     do(tidy(glm(form, data = ., family = 'binomial'))) %>%
#     bind_cols(b[, c("conf.low", "conf.high")]) %>%
#     mutate(OR = paste0(round(exp(estimate), 2), " (", round(exp(conf.low), 2), "-", round(exp(conf.high), 2), ")")) %>%
#     filter(grepl(snp, term)) %>% dplyr::select(!!exposure, OR, term)
#
#   final <- cbind(tf, af, aa[,2])
#
#
# }

# logit2prob <- function(logit){
#   odds <- exp(logit)
#   prob <- odds / (1 + odds)
#   return(prob)
# }




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
pooled_analysis_glm <- function(ds, exposure, covariates, strata = c('sex', 'study_design'), filename_suffix, output_dir) {
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
  } else {
    col_label = c("All", seq(0, number_of_levels - 1))
  }
  
  # save output html from stargazer
  out_html <- stargazer_helper(out,
                               title=paste0(gsub('\\_', '\\\\_', strata), " stratified ", gsub('\\_', '\\\\_', exposure), " main effects"), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=c("(Study specific estimates omitted from table)"), single.row = T)
  
  # write object to html file
  cat(paste(out_html, collapse = "\n"), "\n",
      file = paste0(output_dir, "main_effects_pooled_analysis_", exposure, "_", filename_suffix, ".html"), append = F)
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
pooled_analysis_multinom <- function(ds, exposure, covariates, filename_suffix, output_dir) {
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
                               notes=c("(Study specific estimates omitted from table)"), single.row = T, 
                               add.lines = list(c("Proximal vs. Distal pval", "", formatC(pval_prox_dist, format = "e", digits = 2)),
                                                c("Proximal vs. Rectal pval", "", formatC(pval_prox_rect, format = "e", digits = 2)),
                                                c("Distal vs. Rectal pval",   "", formatC(pval_dist_rect, format = "e", digits = 2))))
  cat(paste(out_html, collapse = "\n"), "\n",
      file = paste0(output_dir, "main_effects_pooled_analysis_", exposure, "_", filename_suffix, ".html"), append = F)
  
}


# pooled_analysis_multinom(ds, exposure, covariates, 'test')




# reload package (convenience)

reload_figifs <- function() {
  devtools::reload(pkgload::inst("figifs"))
}
