#=============================================================================#
# Convenience scripts to clean up epidata for gxescan / glm etc
#=============================================================================#



#' format_data_gxescan
#'
#' @param d Input data - 'master' set I create
#' @param exposure Exposure variable for GxEScanR
#' @param is_e_categorical Is exposure categorical T/F
#' @param min_cell_size Minimum cell size for study removal - tabulation of outcome, exposure, study if categorical. if numeric, tabulation of outcome and study
#'
#' @return Cleaned dataset
#' @export
#'
#' @examples format_data_gxescan(figi_gwas, 'asp_ref', T, 0)
format_data_gxescan <- function(d, exposure, is_e_categorical, min_cell_size = 0) {

  tmp <- d %>%
    dplyr::filter(gxe == 1) %>%
    dplyr::mutate(outcome = fct_drop(outcome)) %>%
    dplyr::select(vcfid, outcome, exposure, age_ref_imp, sex, study_gxe, paste0(rep("PC", 3), seq(1,3))) %>%
    filter(complete.cases(.))

  if (is_e_categorical == T) {
    drops <- data.frame(table(tmp$outcome, tmp[, exposure], tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3))
    return(tmp)
  }
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var2))
    return(tmp)
  }

}
