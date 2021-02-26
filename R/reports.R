# convenience functions to output markdown reports
#' gwis_report
#'
#' @param exposure 
#' @param hrc_version 
#' @param covariates 
#'
#' @return
#' @export
#'
#' @examples
gwis_report <- function(exposure, hrc_version, covariates) {
  rmarkdown::render('~/git/figi/results.Rmd',
                    params = list(exposure = exposure, 
                                  hrc_version = hrc_version,
                                  covariates = covariates), 
                    output_file = glue('~/Dropbox/FIGI/Results/{exposure}_gwis.html'))
}


#' posthoc_report
#'
#' @param exposure 
#'
#' @return
#' @export
#'
#' @examples
posthoc_report <- function(exposure) {
  rmarkdown::render(glue('/home/rak/git/figi/{exposure}_posthoc.Rmd'),
                    params = list(exposure = exposure), 
                    output_file = glue('~/Dropbox/FIGI/Results/{exposure}_posthoc.html'))
}
