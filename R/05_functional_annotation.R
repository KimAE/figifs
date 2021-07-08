#======================================================================#
# these functions are based on Jeroen's functional annotation workflow
#======================================================================#

# run chromatin accesibility marker (scacheri and finucane) perl script
# for calling these functions, need to supply file name with full path (for convenience)
#' output_chromatin
#'
#' @param exposure 
#' @param hrc_version 
#' @param snp 
#' @param path_in 
#' @param path_out 
#'
#' @return
#' @export
#'
#' @examples
output_chromatin <- function(exposure, hrc_version, snp, path_in, path_out) {

    # i expect SNP in format chr:bp:ref:alt
    snp_sep <- unlist(strsplit(snp, split = ":"))
    input_file <- glue("{path_in}/functional_annotation_chr{snp_sep[1]}_{snp_sep[2]}")
    output_file <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{snp_sep[1]}_{snp_sep[2]}")
    file_full <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{snp_sep[1]}_{snp_sep[2]}_chromatin")

    # chromatin marker overlap 
    system(glue("bash ../functional_annotation/colon-crypt-dhs.sh {input_file} {path_out}"))
    system(glue("bash ../functional_annotation/colon-crypt-h3k27ac.sh {input_file} {path_out}"))
    system(glue("bash ../functional_annotation/colon-mucosa-atacseq.sh {input_file} {path_out}"))
    system(glue("bash ../functional_annotation/crc-dhs.sh {input_file} {path_out}"))
    system(glue("bash ../functional_annotation/crc-h3k27ac.sh {input_file} {path_out}"))
    system(glue("bash ../functional_annotation/finucane.sh {input_file} {path_out}"))
    system(glue("perl ../functional_annotation/summarize-overlap-regulatory-region.pl {input_file} {path_out} {file_full}"))

    # return file for targets workflow
    return(glue("{path_out}/figi_{exposure}_{hrc_version}_chr{snp_sep[1]}_{snp_sep[2]}_chromatin.tsv"))
}


# this function outputs two files
# 1) another bed file with ref/alt information (used in liftOver later))
# 2) bed file edited to use in VEP query (tsv)
#' format_vep
#'
#' @param exposure 
#' @param hrc_version 
#' @param snp 
#' @param path_in 
#' @param path_out 
#'
#' @return
#' @export
#'
#' @examples
format_vep <- function(exposure, hrc_version, snp, path_in, path_out) {

    snp_sep <- unlist(strsplit(snp, split = ":")) 
    chr <- snp_sep[1]
    bp <- snp_sep[2]

    # hrc ref/alt information
    hrc <- readRDS(glue("../functional_annotation/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr{chr}.rds"))

    # create another bedfile
    # r2 > 0.2 (from plink call)
    tmp1 <- fread(glue("{path_in}/EUR_chr{chr}_{bp}_ld.ld"))
    tmp2 <- hrc %>%
	filter(POS %in% tmp1$BP_B) %>%
	mutate(V1 = paste0("chr", `#CHROM`),
	    V2 = POS - 1,
	    V3 = POS,
	    V4 = paste(`#CHROM`, POS, REF, ALT, sep = "-")) %>%
    dplyr::select(V1, V2, V3, V4)

    bed_new <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}.bed")
    write.table(tmp2, file = bed_new, quote = F, row.names = F, col.names = F, sep = '\t')

    # output tsv formatted for VEP
    tsv_new <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}.tsv")
    system(glue("perl -lane '@F=split(/\t/); $snp=$F[3]; ($chr, $pos, $ref, $alt) = split(/-/,$snp); print \"$chr\t$pos\t$pos\t$ref/$alt\t+\t$snp\";' {bed_new} | sort -gk1 -gk2 - > {tsv_new}"))

    # return file for targets workflow
    return(glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}.tsv"))
}



# eQTL (GTEx v8)
#' output_eqtl_gtex
#'
#' @param exposure 
#' @param hrc_version 
#' @param snp 
#' @param path_out 
#'
#' @return
#' @export
#'
#' @examples
output_eqtl_gtex <- function(exposure, hrc_version, snp,  path_out) {

    snp_sep <- unlist(strsplit(snp, split = ":"))
    chr <- snp_sep[1]
    bp <- snp_sep[2]
    bed_new <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}.bed")
    bed_liftover <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}_liftover.bed")
    bed_unlifted <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}_unlifted.bed")

    # liftover (v8 uses build38)
    system(glue("../functional_annotation/liftOver {bed_new}  ../functional_annotation/hg19ToHg38.over.chain.gz {bed_liftover} {bed_unlifted}"))

    # query GTEx eQTL
    eqtl <- glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}")
    system(glue("perl ../functional_annotation/01-eqtl-overlap.pl {eqtl} {path_out}")) 
    system(glue("perl ../functional_annotation/02-parse-eqtl-overlap-results.pl {eqtl}"))
    system(glue("perl ../functional_annotation/03-summarize-eqtl-overlap.pl {eqtl}"))

    # return file for targets workflow
    return(glue("{path_out}/figi_{exposure}_{hrc_version}_chr{chr}_{bp}_eQTL_overlap_summary.tsv"))
}


