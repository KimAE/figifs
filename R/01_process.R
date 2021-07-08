#=============================================================================#
# functions to process GxEScanR output
#=============================================================================#

#' process_gxe_output
#'
#' @param exposure 
#' @param hrc_version 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
process_gxe_output <- function(exposure, hrc_version, path) {
  
  rsq_filter <- readRDS("/scratch2/andreeki/gwis_test/data/rsq_0.8_filter_snplist.rds")
  
  out <- do.call(rbind, lapply(list.files(path = path, full.names = T, pattern = glue("FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_chr")), fread, stringsAsFactors = F, data.table = F)) %>%
    dplyr::filter(SNP %in% rsq_filter) %>%
    dplyr::mutate(SNP2 = paste0(Chromosome, ":", Location),
                  chiSqEDGE = chiSqG + chiSqGE,
                  chiSq3df = chiSqG + chiSqGxE + chiSqGE,
                  chiSqEDGE_p = pchisq(chiSqEDGE, df = 2, lower.tail = F),
                  chiSqG_p = pchisq(chiSqG, df = 1, lower.tail = F),
                  chiSqGxE_p = pchisq(chiSqGxE, df = 1, lower.tail = F),
                  chiSq2df_p = pchisq(chiSq2df, df = 2, lower.tail = F),
                  chiSq3df_p = pchisq(chiSq3df, df = 3, lower.tail = F),
                  chiSqGE_p = pchisq(chiSqGE, df = 1, lower.tail = F),
                  chiSqCase_p = pchisq(chiSqCase, df = 1, lower.tail = F),
                  chiSqControl_p = pchisq(chiSqControl, df = 1, lower.tail = F))
  
  file <- glue("{path}/processed/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_results.qs")
  dir.create(file.path(glue("{path}/processed")))
  
  qsave(out, file = file)
  # return(file)
}


# some older results have chr:bp as SNP id.. need to change that
process_gxe_output <- function(exposure, hrc_version, path) {
  
  rsq_filter <- readRDS("/scratch2/andreeki/gwis_test/data/rsq_0.8_filter_snplist.rds")
  
  out <- do.call(rbind, lapply(list.files(path = path, full.names = T, pattern = glue("FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_chr")), fread, stringsAsFactors = F, data.table = F)) %>%
    dplyr::filter(SNP %in% rsq_filter) %>%
    dplyr::mutate(SNP2 = paste0(Chromosome, ":", Location),
                  chiSqEDGE = chiSqG + chiSqGE,
                  chiSq3df = chiSqG + chiSqGxE + chiSqGE,
                  chiSqEDGE_p = pchisq(chiSqEDGE, df = 2, lower.tail = F),
                  chiSqG_p = pchisq(chiSqG, df = 1, lower.tail = F),
                  chiSqGxE_p = pchisq(chiSqGxE, df = 1, lower.tail = F),
                  chiSq2df_p = pchisq(chiSq2df, df = 2, lower.tail = F),
                  chiSq3df_p = pchisq(chiSq3df, df = 3, lower.tail = F),
                  chiSqGE_p = pchisq(chiSqGE, df = 1, lower.tail = F),
                  chiSqCase_p = pchisq(chiSqCase, df = 1, lower.tail = F),
                  chiSqControl_p = pchisq(chiSqControl, df = 1, lower.tail = F))
  
  file <- glue("{path}/processed/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_results.qs")
  dir.create(file.path(glue("{path}/processed")))
  
  qsave(out, file = file)
  # return(file)
}






#' output_clump_pvals
#' 
#' output .txt files that contain SNP and test P values for clumping with plink. 
#' Set lower p-value threshold to 5e-4 to make clumping faster
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param statistic string
#' @param hrc_version string
#' @param path string
#' @param gwas_snps vector of GWAS SNPs
#'
#' @return
#' @export
#'
#' @examples
output_clump_pvals <- function(data, exposure, statistic, hrc_version, path, gwas_snps=NULL) {

    # create output data.frame
    out <- data %>%
	{if (!is.null(gwas_snps)) dplyr::filter(., !SNP2 %in% gwas_snps) else . } %>% 
	dplyr::rename(P = glue("{statistic}_p")) %>%
	dplyr::filter(P < 5e-4) %>%
	dplyr::select(SNP, P)

    # name file differently depending on whether you excluded GWAS snps or not
    if (!is.null(gwas_snps)) {
	file = glue("{path}/FIGI_{hrc_version}_gxeset_{exposure}_{statistic}_no_gwas_ldclump.txt")
    } else {
	file = glue("{path}/FIGI_{hrc_version}_gxeset_{exposure}_{statistic}_ldclump.txt")
    }


    write.table(out, file = file, quote = F, row.names = F, sep = '\t')
    return(file) # for targets workflow (tracks file name)
}




#' output_locuszoom_pvals
#' 
#' output p values for creating regional plots using LocusZoom
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param statistic string
#' @param hrc_version string
#' @param path string
#'
#' @return
#' @export
#'
#' @examples
output_locuszoom_pvals <- function(data, exposure, statistic, hrc_version, path) {

    # create output data.frame
    locuszoom <- data %>%
	dplyr::rename(`P-value` = paste0(statistic, "_p")) %>%
	dplyr::mutate(MarkerName = paste0("chr", Chromosome, ":", Location)) %>%
	dplyr::select(MarkerName, `P-value`)

    file <- glue("{path}/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_{statistic}_locuszoom.txt")
    write.table(locuszoom, file = file,	quote = F, row.names = F, sep = "\t")
    return(file)
}




#' clump_plink
#' 
#' system call to clump p values. uses output from output_clump_pvals function. 
#' returns the log text file so R targets can track it - if the log files changes, then presumably targets will run this component again
#' also remember that LD is estimated using a sample of 1000 FIGI controls
#'
#' @param file string
#' @param clump_p1 numeric
#'
#' @return
#' @export
#'
#' @examples
clump_plink <- function(file, clump_p1 = 5e-8) {
    system(glue("plink --bfile /scratch/andreeki/clump/figi_controls_1000 --memory 8000 --clump {file} --clump-p1 {clump_p1} --clump-r2 0.15 --clump-p2 1 --out {gsub('.txt', '', {file})}"), ignore.stdout = T, ignore.stderr = T)
    # uses the figi_controls_1000 file (random sample of FIGI controls) for LD estimation

    output_file <- glue("{gsub('.txt', '.clumped', {file})}")
    if(file.exists(output_file)) {
	return(output_file)
    } else {
	return(glue("{gsub('.txt', '.log', {file})}"))
    }
}





