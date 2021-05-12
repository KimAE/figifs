
# data <- output
# sizeBin0 = 5
# stats_step1 = 'wtg_gp'
# alpha = 0.05
# binsToPlot = 10


#' create_mapinfo
#'
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
create_mapinfo <- function(y) {
  y %>% 
    arrange(CHR, BP) %>% 
    mutate(mapinfo = seq(unique(bin_number) - 1 + 0.1, unique(bin_number) - 1 + 0.9, length.out = nrow(.)))
}



#' create_twostep_plot
#'
#' @param data 
#' @param binsToPlot 
#' @param stats_step1 
#' @param sizeBin0 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
create_twostep_plot <- function(data, binsToPlot, stats_step1, sizeBin0, alpha) {
  
  # assign SNPs to bins
  m = nrow(data) 
  nbins = floor(log2(m/sizeBin0 + 1)) 
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) 
  sizeBin = c(sizeBin0 * 2^(0:(nbins-1))) 
  endpointsBin = cumsum(sizeBin)
  
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2))
  
  rep_helper <- c(table(grp))
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin 
  alphaBin_dat <- rep(alphaBin, rep_helper)
  
  # prep data
  tmp <- data %>% 
    dplyr::mutate(step1p = .data[[stats_step1]],
                  step2p = wtgxep) %>% 
    dplyr::arrange(step1p) %>% 
    dplyr::mutate(bin_number = as.numeric(grp), 
                  step2p_sig = as.numeric(alphaBin_dat), 
                  log_step2p_sig = -log10(step2p_sig), 
                  log_step2p = -log10(step2p))
  
  significant_hits <- filter(tmp, step2p < step2p_sig)
  
  # output list of bins for plotting
  tmp_plot <- tmp %>% 
    arrange(CHR, BP) %>% 
    group_by(bin_number) %>% 
    group_split()
  
  # add mapinfo
  tmp_plot <- lapply(tmp_plot, create_mapinfo)
  logp_plot_limit = 12
  last.sig = alphaBin[binsToPlot]

  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  bin_to_plot = tmp_plot[[1]]
  plot(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p),
       col = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], '#E41A1C','#377EB8'),
       pch = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], 19, 20),
       cex = ifelse(pull(bin_to_plot, SNP) %in% significant_hits[, 'SNP'], 0.8, 1),
       xlab="Bin number (based on Step 1 p-value)",
       ylab="-log10 GxE p-value)",
       xlim=c(0, binsToPlot),
       ylim=c(0, logp_plot_limit),
       axes=F,
       cex.main = 1,
       cex.axis = 1,
       cex.lab = 1,
       cex.sub = 0.9)
  lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig), col = "black", lwd=1)
  
  # remaining bins
  for(i in 2:binsToPlot) {
    bin_to_plot = tmp_plot[[i]]
    points(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p),
           col = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, '#E41A1C', color[i]),
           pch = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, 19, 20),
           cex = ifelse(pull(bin_to_plot, SNP) %in% significant_hits$SNP, 0.8, 1),
           cex.main = 1,
           cex.axis = 1,
           cex.lab = 1,
           cex.sub = 0.9)
    lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig),
          col = "black",lwd = 1)
  }
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.2, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1)
  axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1)
  title(main = "Two-step plot", sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.2, cex.sub = 1)
  
}
