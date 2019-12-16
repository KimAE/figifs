#=============================================================================#
# Convenience scripts to clean up epidata for gxescan / glm etc
#=============================================================================#

#' format_data_glm
#'
#' @description
#' Format data for analysis.
#'
#' Performs the following steps:
#' * subset gxe == 1
#' * subset !is.na(exposure)
#' * change outcome and sex into numeric variables
#' * only keep variables necessary for scans
#' * removes case only / control only studies. Alternatively, drop studies with low cells counts
#'
#' Note that all main adjustment variables do not have missing values. Sufficient to subset by exposure to create phenotype file
#'
#' @param d Input data (FIGI_EpiData02_master.R)
#' @param exposure Exposure variable for GxEScanR
#' @param is_e_categorical Is exposure categorical - T/F
#' @param min_cell_size Minimum cell size for study removal. Tabulate outcome+exposure+study for categorical variable, tabulate outcome+study for continuous variable
#' @param vars_to_exclude Vector of variables to remove from dataset. Useful for stratified analyses e.g. females only
#'
#' @return Cleaned dataset
#' @export
#'
#' @examples format_data_glm(figi_gwas, 'asp_ref', T, 0, c("energytot"))
format_data_glm <- function(d, exposure, is_e_categorical, min_cell_size = 0, vars_to_exclude = c('energytot_imp'), eur_only=T) {

  vars_to_keep <- c("vcfid", "outcome", exposure, "age_ref_imp", "sex", "energytot_imp", "study_gxe", "PC1", "PC2", "PC3")
  vars_to_keep <- vars_to_keep[!vars_to_keep %in% vars_to_exclude]

  # note that in gxe set, outcome+age_ref_imp+sex+study_gxe+energytot_imp do NOT have missing values
  # OK to subset simply by using is.na(exposure)
  if(eur_only == T) {
    tmp <- d %>%
      dplyr::filter(gxe == 1,
                    EUR == 1,
                    !is.na(get(exposure))) %>%
      dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                    sex = ifelse(sex == "Female", 0, 1))
  } else {
    tmp <- d %>%
      dplyr::filter(gxe == 1,
                    !is.na(get(exposure))) %>%
      dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                    sex = ifelse(sex == "Female", 0, 1))
  }

  # drop zero cells, keep vars_to_keep
  if (is_e_categorical == T) {
    tmp <- mutate(tmp, {{exposure}} := as.numeric(get(exposure)) - 1)

    drops <- data.frame(table(tmp$outcome, tmp[, exposure], tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)

    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe)) %>%
      dplyr::select(vars_to_keep)  }
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var2)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe)) %>%
      dplyr::select(vars_to_keep)  }

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
#' @param exposure Exposure variable for GxEScanR
#'
#' @return Phenotype file for GxEScanR
#' @export
#'
#' @examples format_data_gxescan(gxe, 'asp_ref')
format_data_gxescan <- function(d, exposure) {

  tmp <- d
  ref_study <- as.character(unique(d[, 'study_gxe'])[1])

  for(t in unique(tmp$study_gxe)) {
    tmp[paste0(t)] <- ifelse(tmp$study_gxe==t,1,0)
  }

  # # tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)

}
