#=============================================================================#
# functions to perform some posthoc analyses
# interaction plot
#=============================================================================#



#' iplot_wrapper
#'
#' @param data 
#' @param exposure 
#' @param hrc_version 
#' @param snp 
#' @param covariates 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
iplot_wrapper <- function(data_epi, exposure, hrc_version, snp, covariates, path) {
  
  wdir = glue("{path}/posthoc")
  
  # SNP - need to change the SNP name to match the file (because of binarydosage output)
  tmp <- unlist(strsplit(snp, ":"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]
  
  # convert 12:50610976:C:T to chr12_50610976_C_T_dose
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  
  #  need to merge with dosage information
  data_dose <- qread(glue("{wdir}/dosage_chr{chr}_{bp}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  # check if SNP has to be recoded (for consistency with RERI model)
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
  
  png(glue("{wdir}/interaction_plot_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.png"), height = 720, width = 1280)
  if (is.factor(data[,exposure])) {
    print(interact_plot(model, modx = !! exposure , pred = !! snp_new, plot.points = F, interval = T, outcome.scale = 'link', y.label = 'predicted log odds') + theme(text = element_text(size = 26)))
    # johnson_neyman(model, pred = folate_totqc2, modx = chr2_55255870_C_T, alpha = 0.05)
  } else {
    print(interact_plot(model, modx = !! exposure , pred = !! snp_new, plot.points = F, interval = T, modx.values = c(0,1,2,3), outcome.scale = 'link', y.label = 'predicted log odds') + theme(text = element_text(size = 26)))
  }
  dev.off()
  
  return(glue("{wdir}/interaction_plot_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.png"))
  # saveRDS(out, file = glue("{output_dir}reri_{exposure}_{snp}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
}
