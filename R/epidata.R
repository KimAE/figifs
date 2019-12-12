#=============================================================================#
# Convenience scripts to clean up epidata for gxescan / glm etc
#=============================================================================#



#' format_data_glm
#'
#' @description
#' Format data for analysis.
#' Performs the following steps:
#' * subset gxe == 1
#' * drop unused outcome factor level
#' * selects vcfid, outcome, exposure, and basic covariates (edit as needed)
#' * drops samples with missing information
#' * removes case only / control only studies. Can also drop studies with low cells sizes.
#'
#' @param d Input data - 'master' set I create
#' @param exposure Exposure variable for GxEScanR
#' @param is_e_categorical Is exposure categorical T/F
#' @param min_cell_size Minimum cell size for study removal - tabulation of outcome, exposure, study if categorical. if numeric, tabulation of outcome and study
#'
#' @return Cleaned dataset
#' @export
#'
#' @examples format_data_glm(figi_gwas, 'asp_ref', T, 0)
format_data_glm <- function(d, exposure, is_e_categorical, min_cell_size = 0) {

  tmp <- d %>%
    dplyr::filter(gxe == 1) %>%
    dplyr::mutate(outcome = fct_drop(outcome)) %>%
    dplyr::select(vcfid, outcome, exposure, age_ref_imp, sex, study_gxe, paste0(rep("PC", 3), seq(1,3))) %>%
    filter(complete.cases(.))

  # note - you need to be very careful with factor order when relabeling levels

  if (is_e_categorical == T) {
    drops <- data.frame(table(tmp$outcome, tmp[, exposure], tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    outcome = factor(outcome, labels = c(1,0)),
                    sex = factor(sex, labels = seq(from = 0, length(levels(sex)) - 1)),
                    {{exposure}} := factor(get(exposure), labels = seq(from = 0, length(levels(get(exposure))) - 1)))

    # return(tmp)
  }
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    outcome = factor(outcome, labels = c(1,0)),
                    sex = factor(sex, labels = seq(from = 0, length(levels(sex)) - 1)))
    # return(tmp)
  }

  return(tmp)

}




#' format_data_glm_nosex
#'
#' @description
#' Same function as format_data_glm, but no sex recoding. Meant to be used for HRT or other sex specific analyses..
#' Incidentally - will need an entire set of functions for stratified analyses.
#'
#' @param d Input data - 'master' set I create
#' @param exposure Exposure variable for GxEScanR
#' @param is_e_categorical Is exposure categorical T/F
#' @param min_cell_size Minimum cell size for study removal - tabulation of outcome, exposure, study if categorical. if numeric, tabulation of outcome and study
#'
#' @return Cleaned dataset
#' @export
#'
#' @examples format_data_glm_nosex(figi_gwas, 'asp_ref', T, 0)
format_data_glm_nosex <- function(d, exposure, is_e_categorical, min_cell_size = 0) {

  tmp <- d %>%
    dplyr::filter(gxe == 1) %>%
    dplyr::mutate(outcome = fct_drop(outcome)) %>%
    dplyr::select(vcfid, outcome, exposure, age_ref_imp, sex, study_gxe, paste0(rep("PC", 3), seq(1,3))) %>%
    filter(complete.cases(.))

  if (is_e_categorical == T) {
    drops <- data.frame(table(tmp$outcome, tmp[, exposure], tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    outcome = factor(outcome, labels = c(1,0)),
                    {{exposure}} := factor(get(exposure), labels = seq(from = 0, length(levels(get(exposure))) - 1)))

    # return(tmp)
  }
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    outcome = factor(outcome, labels = c(1,0)))
    # return(tmp)
  }

  return(tmp)

}




#' format_data_glm_energy
#'
#' @description
#' Same as format_data_glm, but adding energytot as part of the dataset
#'
#' using energytot_imp
#'
#' @param d Input data - 'master' set I create
#' @param exposure Exposure variable for GxEScanR
#' @param is_e_categorical Is exposure categorical T/F
#' @param min_cell_size Minimum cell size for study removal - tabulation of outcome, exposure, study if categorical. if numeric, tabulation of outcome and study
#'
#' @return Cleaned dataset
#' @export
#'
#' @examples format_data_glm_energy(figi_gwas, 'asp_ref', T, 0)
format_data_glm_energy <- function(d, exposure, is_e_categorical, min_cell_size = 0) {

  tmp <- d %>%
    dplyr::filter(gxe == 1) %>%
    dplyr::mutate(outcome = fct_drop(outcome)) %>%
    dplyr::select(vcfid, outcome, exposure, age_ref_imp, sex, energytot_imp, study_gxe, paste0(rep("PC", 3), seq(1,3))) %>%
    filter(complete.cases(.))

  if (is_e_categorical == T) {
    drops <- data.frame(table(tmp$outcome, tmp[, exposure], tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    outcome = factor(outcome, labels = c(1,0)),
                    sex = factor(sex, labels = seq(from = 0, length(levels(sex)) - 1)),
                    {{exposure}} := factor(get(exposure), labels = seq(from = 0, length(levels(get(exposure))) - 1)))

    # return(tmp)
  }
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    outcome = factor(outcome, labels = c(1,0)),
                    sex = factor(sex, labels = seq(from = 0, length(levels(sex)) - 1)))
    # return(tmp)
  }

  return(tmp)

}





#' format_data_gxescan
#'
#' @description
#' Create study indicator variables. Set reference study to avoid any issues with gxescanR (just in case). Study subsets vary according to exposure, set the first study in the list as reference
#'
#' Set the interaction exposure as the last variable in the phenotype file
#'
#' Make sure all factors are numeric
#'
#' @param d Output data from the format_data_glm function.
#'
#' @return Phenotype file for GxEScanR
#' @export
#'
#' @examples format_data_gxescan(asp_ref_glm, 'asp_ref')
format_data_gxescan <- function(d, exposure) {

  tmp <- d
  ref_study <- as.character(unique(d[, 'study_gxe'])[1])

  for(t in unique(tmp$study_gxe)) {
    tmp[paste0(t)] <- ifelse(tmp$study_gxe==t,1,0)
  }

  # # tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  if(length(unique(d$sex)) == 1) {
    tmp <- dplyr::select(tmp, -ref_study, -sex, -study_gxe, -exposure, exposure)
  } else {
    tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  }

}
