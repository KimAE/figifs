#=============================================================================#
# functions to create results plots
#=============================================================================#


#' qq_ramwas_edit
#' 
#' I modified the ramwas qq function to include lambda1000 estimates in the output
#'
#' @param x 
#' @param ntests 
#' @param ismlog10 
#' @param ci.level 
#' @param ylim 
#' @param newplot 
#' @param col 
#' @param cex 
#' @param yaxmax 
#' @param lwd 
#' @param axistep 
#' @param col.band 
#' @param makelegend 
#' @param lambda1000 
#' @param xlab 
#' @param ylab 
#'
#' @return
#' @export
#'
#' @examples
qq_ramwas_edit <- function(x, ntests = NULL, ismlog10 = FALSE, ci.level = 0.05,
    ylim = NULL, newplot = TRUE, col = "#D94D4C", cex = 0.5,
    yaxmax = NULL, lwd = 3, axistep = 2, col.band = "#ECA538",
    makelegend = TRUE, lambda1000,
    xlab = expression(paste("–", " log"[10] * "(", italic("P"), "), null")), 
    ylab = expression(paste("–", " log"[10] * "(", italic("P"), "), observed"))) {

    if (methods::is(x, "qqPlotInfo")) {
	qq = x
    }
    else {
	qq = qqPlotPrepare(pvalues = x, ntests = ntests, ismlog10 = ismlog10)
    }
    mx = head(qq$xpvs, 1) * 1.05
    if (is.null(ylim)) {
	my = max(mx, head(qq$ypvs, 1) * 1.05)
	ylim = c(0, my)
    }
    else {
	my = ylim[2]
    }
    if (is.null(yaxmax))
	yaxmax = floor(my)
    if (newplot) {
	plot(x = NA, y = NA, ylim = ylim, xlim = c(0, mx), xaxs = "i",
	    yaxs = "i", xlab = xlab, ylab = ylab, axes = FALSE)
	axis(1, seq(0, mx + axistep, axistep), lwd = lwd)
	axis(2, seq(0, yaxmax, axistep), lwd = lwd)
    }
    abline(a = 0, b = 1, col = "grey", lwd = lwd)
    points(qq$xpvs, qq$ypvs, col = col, cex = cex, pch = 19)
    if (!is.null(ci.level)) {
	if ((ci.level > 0) & (ci.level < 1)) {
	    quantiles = qbeta(p = rep(c(ci.level/2, 1 - ci.level/2),
		    each = length(qq$xpvs)), shape1 = qq$keep, shape2 = qq$ntests -
		qq$keep + 1)
	    quantiles = matrix(quantiles, ncol = 2)
	    lines(qq$xpvs, -log10(quantiles[, 1]), col = col.band,
		lwd = lwd)
	    lines(qq$xpvs, -log10(quantiles[, 2]), col = col.band,
		lwd = lwd)
	}
    }
    if (makelegend) {
	if (!is.null(ci.level)) {
	    legend("topleft", legend = c(expression(paste(italic("P"),
			    " value")), sprintf("%.0f%% CI", 100 - ci.level *
			100)), lwd = c(0, lwd), pch = c(19, NA_integer_),
		lty = c(0, 1), col = c(col, col.band), box.col = "transparent",
		bg = "transparent")
	}
	else {
	    legend("topleft", legend = expression(paste(italic("P"),
			" value")), lwd = 0, pch = 19, lty = 0, col = col,
		box.col = "transparent", bg = "transparent")
	}
    }
    if (!is.null(qq$lambda)) {
	lastr = sprintf("%.3f", qq$lambda)
	# legend("bottom", legend = bquote(lambda == .(lastr)), bty = "n")
	legend("bottom", legend = bquote(lambda == .(lastr) ~~ lambda[1000] == .(lambda1000)), bty = "n")
    }
    return(invisible(qq))
}





#' create_qqplot_ramwas
#' 
#' create qq plots using package ramwas
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param statistic string
#' @param hrc_version string
#' @param path string
#' @param filename_suffix string - not used but there just in case I modify the code
#'
#' @return
#' @export
#'
#' @examples
create_qqplot_ramwas <- function(data, exposure, statistic, hrc_version, path, filename_suffix = "") {

    # df
    if (statistic == 'chiSq2df') {
	df = 2
    } else if(statistic == 'chiSq3df') {
	df = 3
    } else {
	df = 1
    }

    ## get pvalue
    pvals <- data %>%
	dplyr::filter(!is.na(.data[[glue(statistic, "_p")]])) %>%
	dplyr::pull(.data[[glue(statistic, "_p")]])

    ## calculate lambda
    lambda <- round( (median(qchisq(1-pvals, df)) / qchisq(0.5, df)), 4)

    ## calculate lambda1000
    if (statistic %in% c("chiSqControl", "chiSqCase")) {
	lambda1000 <- "NA"
    } else {
	cases <- unique(data[, 'Cases'])
	controls <- unique(data[, 'Subjects']) - unique(data[, 'Cases'])
	total <- cases + controls
	cases1000 <- (cases/total) * 1000
	controls1000 <- (controls/total) * 1000
	lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000)))
	lambda1000 <- sprintf("%.3f", lambda1000)
    }

    # call modified qq_ramwas function
    png(glue("{path}/qq_{exposure}_{hrc_version}_{statistic}.png"), height = 720, width = 1280)
    qq_ramwas_edit(pvals, lambda1000 = lambda1000)
    dev.off()

    # return image file to targets workflow
    return(glue("{path}/qq_{exposure}_{hrc_version}_{statistic}.png"))
}







#' create_manhattanplot_ramwas
#' 
#' much simpler to use than EasyStrata, especially when recreating plots 
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param statistic string
#' @param hrc_version string
#' @param path string
#' @param filename_suffix string 
#' @param sig_line numeric
#' @param gwas_snps vector of GWAS SNPs
#'
#' @return
#' @export
#'
#' @examples
create_manhattanplot_ramwas <- function(data, exposure, statistic, hrc_version, path, filename_suffix = "", sig_line = 5e-8, gwas_snps = NULL) {

    # output data.frame of significant results
  
    if(statistic == 'chiSqGxE') {
      sig_line = 5e-8
    }  else {
      sig_line = sig_line
    }
  
    data_df <- data %>%
	{if (!is.null(gwas_snps)) dplyr::filter(., !SNP2 %in% gwas_snps) else . } %>%
	dplyr::filter(.data[[glue(statistic, "_p")]] <= sig_line)

    if (!is.null(gwas_snps)) {
	saveRDS(data_df , file = glue("{path}/manhattan_{exposure}_{hrc_version}_{statistic}_no_gwas_df.rds"))
    } else {
	saveRDS(data_df , file = glue("{path}/manhattan_{exposure}_{hrc_version}_{statistic}_df.rds"))
    }

    # remove missing stats (causes problems otherwise)
    data_man <- data %>%
	{if (!is.null(gwas_snps)) dplyr::filter(., !SNP2 %in% gwas_snps) else . } %>%
	dplyr::filter(!is.na(.data[[glue(statistic, "_p")]]))

    # run ramwas 
    pvals <- as.numeric(data_man[, glue(statistic, "_p")])
    chr <- as.factor(data_man$Chromosome)
    pos <- as.numeric(data_man$Location)

    man_prep  <- manPlotPrepare(pvalues = pvals, 
	chr = chr, 
	pos = pos)

    if (!is.null(gwas_snps)) {
	filename = glue("{path}/manhattan_{exposure}_{hrc_version}_{statistic}_no_gwas.png")
	png(filename, height = 720, width = 1280)
    } else {
	filename = glue("{path}/manhattan_{exposure}_{hrc_version}_{statistic}.png")
	png(filename, height = 720, width = 1280)
    }
    manPlotFast(man_prep)
    abline(h = -log10(sig_line), col = 'red')
    dev.off()

    return(filename)
}








#' output_expectation_bin_snps
#' 
#' for creating two-step plots that use expectation based binning. outputs vector of SNPs in each bin, for each step 1 filtering statistic we used (chiSqG, chiSqGE, chiSqEDGE)
#' important to note that expectation based binning assumes a total of 1000000 tests to calculate step1 p-value thresholds
#' 
#' @param data Processed GxESranR output
#' @param exposure string
#' @param hrc_version string
#' @param stats_step1 string
#' @param step1_source string - figi or gecco 
#' @param path string
#' @param bin numeric. targets workflow runs this with bins 1:7
#' @param sizeBin0 numeric
#'
#' @return
#' @export
#'
#' @examples
output_expectation_bin_snps <- function(data, exposure, hrc_version, stats_step1, step1_source = "figi", path, bin, sizeBin0 = 5, m = 1000000) {

    # ---- create bin numbers based on step1 statistic ---- #
    # expectation assumption (1 million tests)
    # m = 1000000
    # (fyi number of bins equal 18 with one million tests)
    nbins = floor(log2(m/sizeBin0 + 1))
    nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
    sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) )
    sizeBin_endpt = cumsum(sizeBin)

    # step 1 bin p value cutoffs (see expectation based slides)
    alphaBin_step1 = sizeBin_endpt/m
    alphaBin_step1_cut <- c(-Inf, alphaBin_step1, Inf)

    # create data.frame with expectation based partition assignment
    out <- data %>% 
	dplyr::mutate(bin_number = as.numeric(cut(.data[[glue(stats_step1, "_p")]], breaks = alphaBin_step1_cut))) %>% 
	dplyr::filter(bin_number == bin) %>% 
	dplyr::arrange(.data[[glue(stats_step1, "_p")]]) %>% 
	pull(SNP)

    # create output directory.. (specific to workflow)
    dir.create(file.path(glue("{path}/expectation_bin_dosages")), showWarnings = F)

    if(step1_source == "gecco") {
	filename = glue("{path}/expectation_bin_dosages/expectation_based_snplist_{exposure}_{hrc_version}_{stats_step1}_gecco_bin{bin}.rds")
	saveRDS(out, file = filename)
    } else {
	filename = glue("{path}/expectation_bin_dosages/expectation_based_snplist_{exposure}_{hrc_version}_{stats_step1}_bin{bin}.rds")
	saveRDS(out, file = filename)
    }

    return(filename)
}





#' create_twostep_plot
#' 
#' original function to create two-step plots using WHT
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param covars vector
#' @param binsToPlot numeric - typically set to 7 or 8
#' @param stats_step1 string
#' @param sizeBin0 numeric - typically 5
#' @param alpha numeric - overall alpha, typically 0.05
#' @param path string
#' @param filename_suffix string - currently not in use
#'
#' @return
#' @export
#'
#' @examples
create_twostep_plot <- function(data, exposure, covars, binsToPlot, stats_step1, sizeBin0, alpha, path, exclude_gwas = F, filename_suffix = "") {

    # ------- Functions ------- #
    # plot title and file name
    write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
	gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
	    paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
	    paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
	names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
	return(gxescan_tests[statistic])
    }

    ## add mapinfo so that points in plot reflect chromosome/location on x-axis
    create_mapinfo <- function(x) {
	x %>% 
	    arrange(Chromosome, Location) %>% 
	    mutate(mapinfo = seq(unique(bin_number) - 1 + 0.1, unique(bin_number) - 1 + 0.9, length.out = nrow(.)))
    }

    # ------- create working dataset ------- #

    if(exclude_gwas == T) {
        gwas_snps <- readRDS("../data/conditioning_snps_v20200930_filter_GWAS_SNPS_1000k.rds")
        # remove from data
        data2 <- data %>%
            dplyr::filter(!SNP2 %in% gwas_snps)
    } else {
       data2 <- data
    }

    # assign SNPs to bins
    m = nrow(data2) 
    nbins = floor(log2(m/sizeBin0 + 1)) 
    nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
    sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) 
    endpointsBin = cumsum(sizeBin)

    rk.pv <- c(1:m)
    grp <- ceiling(log(rk.pv/sizeBin0+1,base=2))

    rep_helper <- c(table(grp))
    alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin 
    alphaBin_dat <- rep(alphaBin, rep_helper)

    # prep data
    tmp <- data2 %>% 
	dplyr::mutate(step1p = .data[[paste0(stats_step1, "_p")]],
	    step2p = chiSqGxE_p) %>% 
    dplyr::arrange(step1p) %>% 
    dplyr::mutate(bin_number = as.numeric(grp), 
	step2p_sig = as.numeric(alphaBin_dat), 
	log_step2p_sig = -log10(step2p_sig), 
	log_step2p = -log10(step2p))

    significant_hits <- filter(tmp, step2p < step2p_sig)

    # output list of bins for plotting
    tmp_plot <- tmp %>% 
	arrange(Chromosome, Location) %>% 
	group_by(bin_number) %>% 
	group_split()

    # add mapinfo (see functions)
    tmp_plot <- map(tmp_plot, create_mapinfo)

    cases <- unique(tmp[, 'Cases'])
    controls <- unique(tmp[, 'Subjects']) - unique(tmp[, 'Cases'])
    total <- cases + controls
    logp_plot_limit = 12
    last.sig = alphaBin[binsToPlot]

    # plots
    if(exclude_gwas == T) {
        png_filename_prefix = glue("{path}/twostep_{exposure}_{hrc_version}_{stats_step1}_nogwas")
    } else {
        png_filename_prefix = glue("{path}/twostep_{exposure}_{hrc_version}_{stats_step1}")
    }

    png(glue("{png_filename_prefix}.png"), height = 720, width = 1280)
    color <- rep(c("#377EB8","#4DAF4A"),100)
    par(mar=c(6, 7, 6, 3))

    bin_to_plot = tmp_plot[[1]]
    plot(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p),
	col = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], '#E41A1C','#377EB8'),
	pch = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], 19, 20),
	cex = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], 1.3, 1.7),
	xlab="Bin number for step1 p value",
	ylab="-log10(step2 chiSqGxE p value)",
	xlim=c(0, binsToPlot),
	ylim=c(0, logp_plot_limit),
	axes=F,
	cex.main = 1.7,
	cex.axis = 1.7,
	cex.lab = 1.7,
	cex.sub = 1.7)
    lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig), col = "black", lwd=1)

    # remaining bins
    for(i in 2:binsToPlot) {
	bin_to_plot = tmp_plot[[i]]
	points(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p),
	    col = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, '#E41A1C', color[i]),
	    pch = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, 19, 20),
	    cex = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, 1.3, 1.7),
	    cex.main = 1.7,
	    cex.axis = 1.7,
	    cex.lab = 1.7,
	    cex.sub = 1.7)
	lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig),
	    col = "black",lwd = 1)
    }

    axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.2, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
    axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
    title(main = write_twostep_weightedHT_plot_title(stats_step1, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)

    dev.off()

    saveRDS(significant_hits, file = glue("{png_filename_prefix}_df.rds"))
    return(glue("{png_filename_prefix}.png"))

}







#' create_twostep_plot_eh
#' 
#' modified twostep plots with expectation based bins + step 2 GxE test adjustment for effective number of tests in each bin
#' I don't call this function directly, I use another wrapper defined below..  
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param covars vector
#' @param binsToPlot numeric
#' @param stats_step1 string
#' @param sizeBin0 numeric
#' @param alpha numeric
#' @param path string
#' @param meff_method string - gao or liji
#' @param number_of_snps vector - number of snps for each bin
#' @param number_of_tests vector - effective number of tests for each bin
#'
#' @return
#' @export
#'
#' @examples
create_twostep_plot_eh <- function(data, exposure, covars, binsToPlot, stats_step1, sizeBin0, alpha, path, m = 1000000, meff_method = "", number_of_snps, number_of_tests) { 

    # ------- Some Functions ------- #
    # plot title and file name
    write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
	gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
	    paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
	    paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
	names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
	return(gxescan_tests[statistic])
    }

    # add mapinfo so that points in plot reflect chromosome/location on x-axis
    create_mapinfo <- function(x) {
	x %>% 
	    arrange(Chromosome, Location) %>% 
	    mutate(mapinfo = seq(unique(bin_number) - 1 + 0.1, unique(bin_number) - 1 + 0.9, length.out = nrow(.)))
    }


    # ------- create working dataset ------- #
    # assign SNPs to bins
    # expectation assumption
    # m = 1000000 
    # (number of bins should always equal 18 with one million tests)
    nbins = floor(log2(m/sizeBin0 + 1))
    nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
    # bin sizes
    sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) )
    sizeBin_endpt = cumsum(sizeBin) 
    # step 1 bin p value cutoffs (see expectation based slides)
    alphaBin_step1 = sizeBin_endpt/m
    alphaBin_step1_cut <- c(-Inf, alphaBin_step1, Inf)

    tmp <- data %>%
	dplyr::mutate(step1p = .data[[glue(stats_step1, "_p")]],
	    step2p = chiSqGxE_p,
	    bin_number = as.numeric(cut(step1p, breaks = alphaBin_step1_cut))) %>%
    dplyr::arrange(step1p)


    # add step2 significance threshold, adjusted for effective number of tests
    # make sure SNPs are sorted by step1p!
    meff <- c(number_of_tests[1:binsToPlot])
    alphaBin_step2_simpleM = alpha * 2 ^ -(1:binsToPlot) / meff
    rep_helper <- c(table(tmp[, 'bin_number']))[as.character(1:binsToPlot)]
    rep_helper <- replace(rep_helper, is.na(rep_helper), 0)
    
    # index_helper <- as.numeric(names(rep_helper[!is.na(names(rep_helper))]))
    index_helper <- names(rep_helper[!is.na(names(rep_helper))])
    
    step2p_sig_simpleM <- rep(alphaBin_step2_simpleM, rep_helper)
    
    # ------- Plot ------- #
    
    tmp_plot <- tmp %>%
    dplyr::filter(bin_number <= binsToPlot) %>%
    dplyr::mutate(step2p_sig = step2p_sig_simpleM,
	log_step2p_sig = -log10(step2p_sig), 
	log_step2p = -log10(step2p))

    # output data.frame of significant findings if any
    significant_hits <- dplyr::filter(tmp_plot, step2p < step2p_sig_simpleM)

    # index to name the list (for convenience when plotting)
    list_names <- unique(tmp_plot$bin_number)

    # output list of bins for plotting
    tmp_plot <- tmp_plot %>% 
	arrange(Chromosome, Location) %>% 
	group_by(bin_number) %>% 
	group_split()

    names(tmp_plot) <- list_names

    # subset the label vectors too
    number_of_snps <- number_of_snps[list_names]
    number_of_tests <- number_of_tests[list_names]

    # add mapinfo (see functions)
    tmp_plot <- map(tmp_plot, create_mapinfo)

    cases <- unique(tmp[, 'Cases'])
    controls <- unique(tmp[, 'Subjects']) - unique(tmp[, 'Cases'])
    total <- cases + controls
    logp_plot_limit = 12

    # plots
    png(glue("{path}/twostep_{meff_method}_{exposure}_{hrc_version}_{stats_step1}.png"), height = 720, width = 1280)
    color <- rep(c("#377EB8","#4DAF4A"),100)
    par(mar=c(6, 7, 6, 3))
    bin_to_plot = tmp_plot[[1]]

    plot(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p),
	col = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], '#E41A1C','#377EB8'),
	pch = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], 19, 20),
	cex = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], 1.3, 1.7),
	xlab="Bin number for step1 p value",
	ylab="-log10(step2 chiSqGxE p value)",
	xlim=c(0, binsToPlot),
	ylim=c(0, logp_plot_limit),
	axes=F,
	cex.main = 1.7,
	cex.axis = 1.7,
	cex.lab = 1.7,
	cex.sub = 1.7)
    lines(c(unique(pull(bin_to_plot, bin_number)) - 1, 
	    unique(pull(bin_to_plot, bin_number))), 
	rep(unique(pull(bin_to_plot, log_step2p_sig)), 2), 
	col = "black", lwd = 1)
    text(unique(pull(bin_to_plot, bin_number)) - 1 + 0.5, pull(bin_to_plot, log_step2p_sig)[1] + 2, paste0("SNPs: ", number_of_snps[1]))
    text(unique(pull(bin_to_plot, bin_number)) - 1 + 0.5, pull(bin_to_plot, log_step2p_sig)[1] + 1, paste0("Meff: ", number_of_tests[1]))

    # remaining bins
    for(i in 2:length(tmp_plot)) {
	bin_to_plot = tmp_plot[[i]]

	points(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p),
	    col = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, '#E41A1C', color[i]),
	    pch = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, 19, 20),
	    cex = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, 1.3, 1.7),
	    cex.main = 1.7,
	    cex.axis = 1.7,
	    cex.lab = 1.7,
	    cex.sub = 1.7)
	lines(c(unique(pull(bin_to_plot, bin_number)) - 1,
		unique(pull(bin_to_plot, bin_number))), 
	    rep(unique(pull(bin_to_plot, log_step2p_sig)), 2),
	    col = "black", lwd = 1)
	text(unique(pull(bin_to_plot, bin_number)) - 1 + 0.5, unique(pull(bin_to_plot, log_step2p_sig)) + 2, paste0("SNPs: ", number_of_snps[i]))
	text(unique(pull(bin_to_plot, bin_number)) - 1 + 0.5, unique(pull(bin_to_plot, log_step2p_sig)) + 1, paste0("Meff: ", number_of_tests[i]))
    }

    axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.2, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
    axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
    title(main = write_twostep_weightedHT_plot_title(stats_step1, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)

    dev.off()

    saveRDS(significant_hits, file = glue("{path}/twostep_{meff_method}_{exposure}_{hrc_version}_{stats_step1}_df.rds"))
    return(glue("{path}/twostep_{meff_method}_{exposure}_{hrc_version}_{stats_step1}.png"))
}



#' simplem_wrap
#' 
#' wrapper to create twostep plots with expectation based binning + adjustment of step2 tests by bin specific effective number of tests
#' the name might be confusing since simplem is the gao method but don't change it at this point. 
#' 
#' I modified the code to rely on the poolr package. I think for Gao method, it uses default window size (which is ok)
#' Also, instead of calculating r^2 per chromosome, now perform combined for all SNPs. 
#' 
#' Also modified code to accommodate different assumed M for expectation based hybrid plots
#'
#' @param data Processed GxEScanR output
#' @param exposure string
#' @param covariates vector
#' @param simplem_step1_statistic string
#' @param path string
#' @param meff_method string - gao or liji
#' @param gwas_snps vector - to remove GWAS SNPs if desired
#' @param gao_pca_cutoff numeric - default is 0.995
#'
#' @return
#' @export
#'
#' @examples
simplem_wrap <- function(data, exposure, hrc_version, covariates, simplem_step1_statistic, path, meff_method, exclude_gwas = F, gao_pca_cutoff = 0.995) {

    # create list of bin SNP dosage data.frames
    # FYI - 'glue' doesn't seem regex friendly
    files_input <- mixedsort(list.files(path, pattern = paste0(paste0("expectation_based_snplist_", exposure, "_", hrc_version, "_", simplem_step1_statistic, "_bin"), "(?:.+)", "output.rds"), full.names = T)) 
    files_list <- map(files_input, ~ readRDS(.x)) 

    # ------------------------------------------- #
    # function to filter GWAS SNPs if desired
    remove_snps <- function(zz, gwas_snps_vector) {
        zznames <- substr(names(zz), 1, nchar(names(zz)) - 4)
        zz_index <- !zznames %in% gwas_snps_vector
        zz_out <- zz[, zz_index]
        return(zz_out)
    }

    # filter out GWAS SNPs if provided
    # be careful - make sure gwas_snps is provided as chr:bp - so i can change it to Xchr.bp
    if(exclude_gwas == T) {
	# GWAS SNPs to remove (vector of chr:bp)
	gwas_snps <- readRDS("../data/conditioning_snps_v20200930_filter_GWAS_SNPS_1000k.rds")
	gwas_snps_v2 <- paste0("X", gsub(":", ".", gwas_snps))

	# remove from bins
	files_list <- map(files_list, ~ remove_snps(.x, gwas_snps_v2))

	# remove from data
	data2 <- data %>%
	    dplyr::filter(!SNP2 %in% gwas_snps)
    } else {
       data2 <- data
    }    

    # output correlation matrix for each bin 
    # overall, no longer doing it by chromosome
    snps_cor <- map(files_list, ~ cor(.x[,-1])) # '-1' to remove vcfid column
    
    # ------------------------------------------ #
    number_of_snps <- map_int(files_list, ~ ncol(.x[,-1]))

    # calculate effective number of tests, output results in vector
    if(meff_method == "gao") {
        number_of_tests <- map_dbl(snps_cor, ~ poolr::meff(R = .x, method = meff_method, C = gao_pca_cutoff))
    } else if(meff_method == "lea") {
	number_of_tests <- map_dbl(snps_cor, ~ meff_lea(s = .x))
    } else {
	number_of_tests <- map_dbl(snps_cor, ~ poolr::meff(R = .x, method = meff_method))
    }

    # helper for writing filename
    if(meff_method == "gao") {
	meff_method_out = glue("gao_{gao_pca_cutoff}")
    } else {
	meff_method_out = meff_method
    }

    # helper for writing 'gwas' if you filtered gwas 
    if(exclude_gwas) {
	meff_method_out = glue(meff_method_out, "_excludeGWAS")
    }

    # call plotting function
    output_filename <- create_twostep_plot_eh(
            data = data2,
            exposure = exposure,
            covars = covariates,
            binsToPlot = 7,
            stats_step1 = simplem_step1_statistic,
            sizeBin0 = 5,
            alpha = 0.05,
            path = "output",
            meff_method = meff_method_out, # for the file name text
            number_of_snps = number_of_snps,
            number_of_tests = number_of_tests)

    return(output_filename)
}



# function to apply lea method (the other li)
meff_lea <- function(s) {
    # s is cor(snps)
    y <- eigen(s)
    out <- floor(ncol(s) - sum((y$values > 1) * (y$values - 1)))
    return(out)
}







#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# create functional annotation plots 
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' output_functional_plot
#'
#' @param snp 
#' @param ld_window 
#' @param plot_window 
#' @param panel 
#' @param r2 
#'
#' @return
#' @export
#'
#' @examples
output_functional_plot <- function(snp, ld_window = 500000, plot_window = 500000, panel = "kg", r2 = 0.5) {

  # path for plink files (LD calculation) and writing out bed, ini, pdf, and png files (folder contains all functional annotation stuff)
  # note that the tmp1.ini and tmp2.ini are hardcoded (just find a permanent place for them)
  if(panel == "kg") {
      eur_path = "/project/dconti_250/locuszoom/data/1000G/genotypes/2014-10-14/EUR"
  } else if(panel == "hrc") {
      eur_path = "/project/dconti_250/locuszoom/data/1000G/genotypes/2014-10-14/ASN"
  }

  wdir = "output/functional_plot"
  tdir = "../data"

  # create the wdir
  dir.create(file.path(glue("{wdir}")), showWarnings = F)

  tmp <- unlist(strsplit(snp, ":"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]

  snpid          <- glue("chr{chr}:{bp}")
  snpid_filename <- glue("chr{chr}_{bp}")

  # get LD SNPs for the locus you're plotting
  system(glue("plink --allow-extra-chr --bfile {eur_path}/chr{chr} --ld-snp {snpid} --out {wdir}/EUR_{snpid_filename}_ld --r2 --ld-window {ld_window}"))

  # create bed file for functional annotation plot
  # filter by r2 < 0.5
  bed_out <- read.table(glue("{wdir}/EUR_{snpid_filename}_ld.ld"), header = T) %>%
    dplyr::filter(R2 >= r2) %>%
    mutate(v1 = paste0('chr', CHR_B),
           v2 = BP_B - 1,
           v3 = BP_B,
           v4 = paste0(CHR_B, "_", BP_B),
           v5 = round(R2*1000),
           v6 = '.',
           v7 = BP_B - 1,
           v8 = BP_B - 1,
           v9 = case_when(R2 < 0.2             ~ '0,0,128',
	                  R2 >= 0.2 & R2 < 0.4 ~ '135,206,250',
			  R2 >= 0.4 & R2 < 0.6 ~ '0,255,0', 
			  R2 >= 0.6 & R2 < 0.8 ~ '255,165,0', 
			  R2 >= 0.8            ~ '255,0,0')) %>%
    filter(!duplicated(.)) %>%
    dplyr::select(v1,v2,v3,v4,v5,v6,v7,v8,v9)

  write.table(bed_out,
              file = glue("{wdir}/functional_annotation_{snpid_filename}.bed"),
              quote = F, row.names = F, col.names = F, sep = '\t')

  # create bed file for annotation main hit with vertical line
  bed_chr_vlines <- bed_out %>%
    filter(v3 == bp)

  write.table(bed_chr_vlines,
              file = glue("{wdir}/functional_annotation_{snpid_filename}_vlines.bed"),
              quote = F, row.names = F, col.names = F, sep = '\t')

  # create .ini files for pygenometracks
  cat(paste("[snps]",
            glue("file={wdir}/functional_annotation_{snpid_filename}.bed"),
            glue("title = r^2 > {r2}"),
            "height = 1",
            "color = bed_rgb",
            "border_color=none",
            "labels=false",
            "display=collapsed",
            "fontsize=11",
            "file_type=bed", sep = '\n'),
      file = glue("{wdir}/functional_annotation_{snpid_filename}_track.ini"), append = F)


  cat(paste("\n",
            "[vlines]",
            glue("file={wdir}/functional_annotation_{snpid_filename}_vlines.bed"),
            "type=vlines",
            "labels=true",
            "\n", sep = '\n'), # this is only available in updated version of pygenometracks..
      file = glue("{wdir}/functional_annotation_{snpid_filename}_track.ini"), append = T) # appending to the previous file

  # assemble .ini file using system call
  system(glue("cat {tdir}/tmp1.ini {wdir}/functional_annotation_{snpid_filename}_track.ini {tdir}/tmp2.ini > {wdir}/functional_annotation_{snpid_filename}.ini"))

  # create a bash script to call under R
  # system calls don't seem to work
  # create shell script
  cat(glue("#!/bin/bash",
	   #"conda init bash", 
	   #"conda activate pygenometracks", 
           "pyGenomeTracks --tracks {wdir}/functional_annotation_{snpid_filename}.ini --fontSize 8 --dpi 60 --region chr{chr}:{bp-plot_window}-{bp+plot_window} --outFileName {wdir}/functional_annotation_{snpid_filename}.pdf",
           "gs -o {wdir}/functional_annotation_{snpid_filename}.png -sDEVICE=png16m -dTextAlphaBits=4 -r300 -dLastPage=1 {wdir}/functional_annotation_{snpid_filename}.pdf",
           .sep = "\n"), file = glue("{wdir}/functional_annotation_{snpid_filename}.sh"), append = F)

  # call shell script
  system(glue("bash  {wdir}/functional_annotation_{snpid_filename}.sh"))

  # return the png file
  return(glue("{wdir}/functional_annotation_{snpid_filename}.png"))

}




#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# create locuszoom plots
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' output_locuszoom_plot
#'
#' @param snp 
#' @param exposure 
#' @param hrc_version 
#' @param statistic 
#' @param plot_window 
#' @param path 
#' @param panel 
#'
#' @return
#' @export
#'
#' @examples
output_locuszoom_plot <- function(snp, exposure, hrc_version, statistic, plot_window = 500, path, panel = 'kg') {
  
  # path for plink files (LD calculation) and writing out bed, ini, pdf, and png files (folder contains all functional annotation stuff)
  # note that the tmp1.ini and tmp2.ini are hardcoded (just find a permanent place for them)
  annotation_file = "../data/locuszoom_gecco_uk_asia_gwas_annotation_203"
  wdir = glue("{path}/locuszoom_plot")
  # wdir = "output/locuszoom_plot"
  pval_file = glue("data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_{statistic}_locuszoom.txt")
  
  # create the wdir 
  dir.create(file.path(glue("{wdir}")), showWarnings = F)
  
  tmp <- unlist(strsplit(snp, ":"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]
  
  snpid          <- glue("chr{chr}:{bp}")
  snpid_filename <- glue("chr{chr}_{bp}")
  
  # gwas results
  gwas_loci_tmp <- read.table(glue("{annotation_file}.txt"), header = T, stringsAsFactors = F)
  
  # need to create blank denote-markers-file. you get errors if you specify a file and doesn't overlap with your hit
  tmp1 <- c(snpid, "", "white")
  gwas_loci <- rbind(gwas_loci_tmp, tmp1)
  write.table(gwas_loci, file = glue("{wdir}/locuszoom_plot_{exposure}_{hrc_version}_{statistic}_{snpid_filename}.txt"), quote = F, row.names = F, sep = "\t")
 
  if(panel == 'kg') {
      pop = "EUR"
  } else if (panel == 'hrc') {
      pop = "ASN"
  }

  system(glue("/project/dconti_250/locuszoom/bin/locuszoom --snpset NULL --metal {pval_file} --pop {pop} --build hg19 --source 1000G_Nov2014 --refsnp {snpid} --flank {plot_window}kb --prefix {wdir}/locuszoom_plot_{exposure}_{hrc_version}_{statistic}  --no-date title='{exposure} x {snpid} - {statistic}' signifLine=7.30103 signifLineColor='blue' --denote-markers-file {wdir}/locuszoom_plot_{exposure}_{hrc_version}_{statistic}_{snpid_filename}.txt"))
  
  system(glue("gs -o {wdir}/locuszoom_plot_{exposure}_{hrc_version}_{statistic}_{snpid_filename}.png -sDEVICE=png16m -dTextAlphaBits=4 -r300 -dLastPage=1 {wdir}/locuszoom_plot_{exposure}_{hrc_version}_{statistic}_{snpid_filename}/chr{chr}_{bp-(plot_window*1000)}-{bp+(plot_window*1000)}.pdf"))
  
  # return png file
  return(glue("{wdir}/locuszoom_plot_{exposure}_{hrc_version}_{statistic}_{snpid_filename}.png"))
}




