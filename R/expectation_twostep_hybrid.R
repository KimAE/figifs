#' Meff_PCA
#'
#' From Gao et al 2008 (simpleM)
#'
#' @param eigenValues 
#' @param percentCut 
#'
#' @return
#' @export
#'
#' @examples
Meff_PCA <- function(eigenValues, percentCut){
  totalEigenValues <- sum(eigenValues)
  myCut <- percentCut*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if(myEigenSum <= myCut){
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    }
    else{
      break
    }
  }	
  return(index_Eigen)
}




#' inferCutoff
#' 
#' From Gao et al 2008 (simpleM)
#'
#' @param dt_My 
#'
#' @return
#' @export
#'
#' @examples
inferCutoff <- function(dt_My, PCA_cutoff){
  CLD <- cor(dt_My)
  eigen_My <- eigen(CLD)
  
  # PCA approach
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}



#' meff_r
#' 
#' wrapper function to run simpleM
#'
#' @param dat output from binarydosage getSNPs function
#' @param PCA_cutoff Gao et al cutoff (0.995)
#' @param fixLength Default at 150
#'
#' @return integer. number of effective tests for the input dosages
#' @export
#'
#' @examples

meff_r <- function(dat, PCA_cutoff = 0.995, fixLength = 150) {
  
  # separate data.frame into chromosome specific lists of SNP vectors 
  # assume that columns are in chr:bp order (they should be provided input vector for binarydosage was ordered)
  dat_colnames <- colnames(dat)[!colnames(dat) %in% "vcfid"]
  tmp1 <- sapply(dat_colnames, function(x) strsplit(x, split = "\\."))
  tmp2 <- sapply(tmp1, function(x) as.numeric(gsub("X", "", x[[1]][1])))
  snps_list <- split(dat_colnames, tmp2) # names of SNPs by chromosome
  
  # nested function to apply Meff method to chromosome specific lists of SNP vectors  (`snps_list`)
  run_meff <- function(snps_vector) {
    
    workdata <- dplyr::select(dat, vcfid, all_of(snps_vector)) %>% 
      pivot_longer(-vcfid, ) %>% 
      pivot_wider(names_from = vcfid, values_from = value) %>% 
      separate(name, into = c("chr", "bp", "ref", "alt"), remove = F) %>% 
      mutate(chr = as.numeric(gsub("X", "", chr)), 
             bp = as.numeric(bp)) %>% 
      arrange(chr, bp) %>% 
      dplyr::select(-chr, -bp, -ref, -alt, -name)
    
    numLoci <- length(pull(workdata, 1))
    
    simpleMeff <- NULL
    
    fixLength <- fixLength
    i <- 1
    myStart <- 1
    myStop <- 1
    iteration <- 0
    
    while(myStop < numLoci){
      myDiff <- numLoci - myStop 
      if(myDiff <= fixLength) break
      
      myStop <- myStart + i*fixLength - 1
      snpInBlk <- t(workdata[myStart:myStop, ])
      MeffBlk <- inferCutoff(snpInBlk, PCA_cutoff)
      simpleMeff <- c(simpleMeff, MeffBlk)
      myStart <- myStop+1
      iteration <- iteration+1
      print(iteration)
    }
    
    snpInBlk <- t(workdata[myStart:numLoci, ])
    MeffBlk <- inferCutoff(snpInBlk, PCA_cutoff)
    simpleMeff <- c(simpleMeff, MeffBlk)
    # cat("Total number of SNPs is: ", numLoci, "\n")
    # cat("Inferred Meff is: ", sum(simpleMeff), "\n")
    
    return(sum(simpleMeff))
  }
  
  # apply nested function here
  # return integer
  eff_tests_list <- map(snps_list, run_meff)
  return(do.call(sum, eff_tests_list))
  
}
















#' format_data_twostep_expectation_hybrid
#' 
#' expectation based twostep method requires calculation of effective number of tests for multiple testing adjustment
#' use this function to assign SNPs into bins IN EXPECTATION (uniform distribution of 1 million tests), then you can use output to write snp vectors to file and extract dosages using binaryDosage package
#'
#' @param dat 
#' @param stats_step1 
#' @param sizeBin0 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
format_data_twostep_expectation_hybrid <- function(dat, stats_step1, sizeBin0, alpha) {
  
  # function to calculate p values from chisq stats
  create_pval_info <- function(dat, stats_step1, df=1) {
    tmp <- data.table(dat)
    tmp[, step1p := pchisq(tmp[, get(stats_step1)], df = df, lower.tail = F)
        ][
          , step2p := pchisq(tmp[, get('chiSqGxE')],  df = 1, lower.tail = F)
          ][
            , y := -log10(step2p)
            ][
              order(step1p)
              ][
                , MapInfo := Location
                ]
  }
  
  if(stats_step1 == 'chiSqEDGE') {
    pv <- create_pval_info(dat, stats_step1, df = 2)
  } else {
    pv <- create_pval_info(dat, stats_step1, df = 1)
  }
  
  # assign SNPs to bins
  m = 1000000 # assuming 1 million tests in expectation
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} # number of bins should always equal 18 with one million tests
  
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin_step1 = endpointsBin/1000000 # step 1 bin p value cutoffs (see Jim's slides)
  alphaBin = alpha * 2 ^ -(1:nbins) # this is how it's laid out in jim's proposal
  
  alphaBinCut <- c(-Inf, alphaBin_step1, Inf)
  pv[ , grp:=as.numeric(cut(pv[,step1p], breaks = alphaBinCut )) ] # actually assign groups based on step 1 pvalues
  # rep_helper <- c(table(pv[,grp]))
  # test <- alphaBin / rep_helper
  # pv[ , step2p_sig:=rep(test, rep_helper)]
  
  # return the data.table
  return(pv)
}







#' create_twostep_plot_expectation_hybrid
#' 
#' create plot with significant threshold based on effective number of tests (simpleM)
#'
#' @param dat 
#' @param exposure 
#' @param covars 
#' @param binsToPlot 
#' @param statistic 
#' @param filename_suffix 
#'
#' @return
#' @export
#'
#' @examples
create_twostep_plot_expectation_hybrid <- function(dat, exposure, covars, binsToPlot, statistic, output_dir, filename_suffix = "") {
  
  cases <- unique(data.frame(dat[, 'Cases']))
  controls <- unique(data.frame(dat[, 'Subjects'])) - unique(data.frame(dat[, 'Cases']))
  total <- cases + controls
  
  # plot title and file name
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }
  
  # create list of bin specific data.tables
  logp_plot_limit = 12
  
  glist<-list()
  setkey(dat,grp)
  for(i in 1:binsToPlot){
    z <- dat[J(i)][order(Chromosome, Location)]
    z[, log_binalpha := -1*log10(min(z[,step2p_sig_simpleM]))] # -log10 of bin specific alpha
    z[, logstep2p := -log10(step2p)]
    z[, mapinfo_tmp := seq(Chromosome, Chromosome+1, length.out = .N), by = Chromosome]
    z[, mapinfo := 0.9*( (z[,mapinfo_tmp]-min(z[,mapinfo_tmp])) / (max(z[,mapinfo_tmp])-min(z[,mapinfo_tmp])) ) + 0.05 + i - 1]
    glist[[i]] <- z[order(Chromosome, Location)]
    rm(z)
  }
  
  significant_hits <- dat[step2p <= step2p_sig_simpleM]
  
  # plots
  png(paste0(output_dir, "twostep_wht_", statistic, "_", exposure,  filename_suffix, ".png"), height = 720, width = 1280)
  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  
  # first bin
  plot(glist[[1]][,mapinfo], glist[[1]][,logstep2p],
       col = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], '#E41A1C','#377EB8'),
       pch = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 19, 20),
       cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0,binsToPlot),
       ylim=c(0,logp_plot_limit),
       axes=F,
       cex.main = 1.7,
       cex.axis = 1.7,
       cex.lab = 1.7,
       cex.sub = 1.7)
  lines(glist[[1]][,mapinfo], glist[[1]][,log_binalpha], col = "black", lwd=1)
  text(glist[[1]][,mapinfo][1]+0.3, glist[[1]][,log_binalpha][1]+2, paste0("SNPs: ", number_of_snps[1]))
  text(glist[[1]][,mapinfo][1]+0.3, glist[[1]][,log_binalpha][1]+1, paste0("Meff: ", number_of_tests[1]))
  
  # remaining bins
  for(i in 2:binsToPlot){
    points(glist[[i]][,mapinfo], glist[[i]][,logstep2p],
           col = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, '#E41A1C', color[i]),
           pch = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, 19, 20),
           cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
           cex.main = 1.7,
           cex.axis = 1.7,
           cex.lab = 1.7,
           cex.sub = 1.7)
    lines(glist[[i]][,mapinfo], glist[[i]][,log_binalpha],
          col = "black",lwd = 1)
    text(glist[[i]][,mapinfo][i]+0.3, glist[[i]][,log_binalpha][i]+2, paste0("SNPs: ", number_of_snps[i]))
    text(glist[[i]][,mapinfo][i]+0.3, glist[[i]][,log_binalpha][i]+1, paste0("Meff: ", number_of_tests[i]))
  }
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.2, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
  
}










#' plot_twostep_eh
#' 
#' Create two-step weighted hypothesis testing plot using expectation based hybrid method. Must supply vectors of number of SNPs and number of effective tests. 
#'
#' @param dat 
#' @param exposure 
#' @param covars 
#' @param binsToPlot 
#' @param stats_step1 
#' @param sizeBin0 
#' @param alpha 
#' @param output_dir 
#' @param filename_suffix 
#' @param number_of_snps 
#' @param number_of_tests 
#'
#' @return
#' @export
#'
#' @examples

plot_twostep_eh <- function(x, exposure, covars, binsToPlot, stats_step1, sizeBin0, alpha, output_dir, filename_suffix = "", number_of_snps, number_of_tests) { 
  
  # ------- Some Functions ------- #
  ## plot title and file name
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }
  
  ## add mapinfo so that points in plot reflect chromosome/location on x-axis
  # create_mapinfo <- function(x) {
  #   
  #   if (nrow(x) == 1) {
  #     x %>% 
  #       mutate(mapinfo_tmp = 0.05, 
  #              mapinfo = 0.05)
  #   } else {
  #     x %>% 
  #       group_by(Chromosome) %>% 
  #       mutate(mapinfo_tmp = seq(unique(Chromosome), unique(Chromosome) + 1, length.out = n())) %>% 
  #       ungroup() %>% 
  #       mutate(mapinfo = 0.9*( (mapinfo_tmp-min(mapinfo_tmp)) / (max(mapinfo_tmp)-min(mapinfo_tmp))) + 0.05 + bin_number - 1)
  #   }
  # }
  
  create_mapinfo <- function(x) {
    x %>% 
      arrange(Chromosome, Location) %>% 
      mutate(mapinfo = seq(unique(bin_number) - 1 + 0.1, unique(bin_number) - 1 + 0.9, length.out = nrow(.)))
  }
  
  
  
  
  # ------- create working dataset ------- #
  # assign SNPs to bins
  ## expectation assumption
  m = 1000000 
  ## (number of bins should always equal 18 with one million tests)
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
  ## bin sizes
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) )
  sizeBin_endpt = cumsum(sizeBin) 
  ## step 1 bin p value cutoffs (see expectation based slides)
  alphaBin_step1 = sizeBin_endpt/1000000
  alphaBin_step1_cut <- c(-Inf, alphaBin_step1, Inf)
  
  # working data.frame
  # number of effective tests is set to bonferroni for bins 11:Inf
  # (will never plot those bins, and calculating effective number of tests is troublesome)
  tmp <- x %>% 
    mutate(step1p = case_when(stats_step1 == "chiSqEDGE" ~ pchisq(.data[[stats_step1]], df = 2, lower.tail = F),
                              TRUE ~ pchisq(.data[[stats_step1]], df = 1, lower.tail = F)), 
           step2p = pchisq(.data[['chiSqGxE']], df = 1, lower.tail = F), 
           bin_number = as.numeric(cut(step1p, breaks = alphaBin_step1_cut))) %>% 
    arrange(step1p)
  
  
  ## add step2 significance threshold, adjusted for effective number of tests
  ## make sure SNPs are sorted by step1p!
  meff <- c(number_of_tests[1:binsToPlot])
  alphaBin_step2_simpleM = alpha * 2 ^ -(1:binsToPlot) / meff
  rep_helper <- c(table(tmp[, 'bin_number']))[as.character(1:binsToPlot)]
  rep_helper <- replace(rep_helper, is.na(rep_helper), 0)
  
  # index_helper <- as.numeric(names(rep_helper[!is.na(names(rep_helper))]))
  index_helper <- names(rep_helper[!is.na(names(rep_helper))])
  
  step2p_sig_simpleM <- rep(alphaBin_step2_simpleM, rep_helper)
  
  # ------- Plot ------- #
  tmp_plot <- tmp %>%
    filter(bin_number <= binsToPlot) %>% 
    mutate(step2p_sig = step2p_sig_simpleM,
           log_step2p_sig = -log10(step2p_sig), 
           log_step2p = -log10(step2p)) 
  
  ## output data.frame of significant findings if any
  significant_hits <- filter(tmp_plot, step2p < step2p_sig_simpleM)
  
  ## index to name the list (for convenience when plotting)
  list_names <- unique(tmp_plot$bin_number)
  
  
  ## output list of bins for plotting
  tmp_plot <- tmp_plot %>% 
    arrange(Chromosome, Location) %>% 
    group_by(bin_number) %>% 
    group_split()
  
  names(tmp_plot) <- list_names
  
  # subset the label vectors too
  number_of_snps <- number_of_snps[list_names]
  number_of_tests <- number_of_tests[list_names]
  
  
  ## add mapinfo (see functions)
  tmp_plot <- map(tmp_plot, create_mapinfo)
  
  cases <- unique(tmp[, 'Cases'])
  controls <- unique(tmp[, 'Subjects']) - unique(tmp[, 'Cases'])
  total <- cases + controls
  logp_plot_limit = 12
  
  # plots
  png(glue(output_dir, "twostep_wht_{stats_step1}_{exposure}{filename_suffix}.png"), height = 720, width = 1280)
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
  # lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig), col = "black", lwd=1)
  text(unique(pull(bin_to_plot, bin_number)) - 1 + 0.5, pull(bin_to_plot, log_step2p_sig)[1] + 2, paste0("SNPs: ", number_of_snps[1]))
  text(unique(pull(bin_to_plot, bin_number)) - 1 + 0.5, pull(bin_to_plot, log_step2p_sig)[1] + 1, paste0("Meff: ", number_of_tests[1]))
  # text(pull(bin_to_plot, mapinfo)[1]+0.3, pull(bin_to_plot, log_step2p_sig)[1]+2, paste0("SNPs: ", number_of_snps[1]))
  # text(pull(bin_to_plot, mapinfo)[1]+0.3, pull(bin_to_plot, log_step2p_sig)[1]+1, paste0("Meff: ", number_of_tests[1]))
  
  # remaining bins
  for(i in 2:length(tmp_plot)){
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
  
  # return data.frame of significant results!
  saveRDS(significant_hits, file = glue(output_dir, "twostep_wht_{stats_step1}_{exposure}{filename_suffix}_df.rds"))
  return(significant_hits)
  
}




#' data_twostep_eh
#' 
#' output gxescan results with bin numbers based on step1p statistic
#'
#' @param dat dataframe - gxescan output
#' @param stats_step1 string - step1 filtering statistic 
#'
#' @return
#' @export
#'
#' @examples
data_twostep_eh <- function(dat, stats_step1, sizeBin0 = 5) { 

  # ------- create working dataset ------- #
  # assign SNPs to bins
  ## expectation assumption
  m = 1000000 
  ## (number of bins should always equal 18 with one million tests)
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
  ## bin sizes
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) )
  sizeBin_endpt = cumsum(sizeBin) 
  ## step 1 bin p value cutoffs (see expectation based slides)
  alphaBin_step1 = sizeBin_endpt/1000000
  alphaBin_step1_cut <- c(-Inf, alphaBin_step1, Inf)
  
  # working data.frame
  # number of effective tests is set to bonferroni for bins 11:Inf
  # (will never plot those bins, and calculating effective number of tests is troublesome)
  tmp <- dat %>% 
    mutate(step1p = case_when(stats_step1 == "chiSqEDGE" ~ pchisq(.data[[stats_step1]], df = 2, lower.tail = F),
                              TRUE ~ pchisq(.data[[stats_step1]], df = 1, lower.tail = F)), 
           step2p = pchisq(.data[['chiSqGxE']], df = 1, lower.tail = F), 
           bin_number = as.numeric(cut(step1p, breaks = alphaBin_step1_cut))) %>% 
    arrange(step1p)
  
  return(tmp)
}














