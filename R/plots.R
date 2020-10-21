#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots
# Two-Step Weighted Hypothesis Test Plots
#-----------------------------------------------------------------------------#

#' calculate_pval
#'
#' Takes a vector of ChiSq statistics, output a vector of p values
#'
#' @param chiSq Vector of chiSq statistics - this is most applicable to GxEScanR results (e.g. a single column)
#' @param df ChiSq distribution degrees of freedom
#' @return Returns a vector of chi-square p values
#' @export
#'
#' @examples calculate_pval(3.6, 1)
#' ex <- rchisq(10, df=1)
#' calculate_pval(ex, df=1)
calculate_pval <- function(chiSq, df) {
  pchisq(chiSq, df = df, lower.tail = F)
}


#' write_plot_title
#'
#' Function to create figure titles based on the results statistic plotted
#'
#' @param stat GxEScan chi-square statistic (look up how to set/list parameter choices)
#' @param exposure Exposure variable (string)
#' @param covars Covariate string vector
#' @param N Sample size
#'
#' @return Plot title string
#' @export
#'
#' @examples write_plot_title('chiSqG', 'aspirin', c('sex', 'age'), 100000)
write_plot_title <- function(stat, exposure, covars, N) {
  gxescan_tests <- c(paste0("G Main Effects Results (N = ", N, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                     paste0("GxE Results (N = ", N, ")\noutc ~ G*", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("2DF Results (N = ", N, ")\noutc ~ G+G*", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("G|E Results (N = ", N, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("Case-Only Results (N = ", N, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("Control-Only Results (N = ", N, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("3DF Results (N = ", N, ")\nchiSqG+chiSqGxE+chiSqGE"))
  names(gxescan_tests) <- c("chiSqG", "chiSqGxE", "chiSq2df", "chiSqGE", "chiSqCase", "chiSqControl", "chiSq3df")
  return(gxescan_tests[stat])
}


#' create_qqplot
#'
#' Function creates a QQ plot and outputs a PNG file
#'
#' @param dat Input data
#' @param exposure Exposure variable
#' @param stat GxEScan results chi-square statistic to plot
#' @param df Degrees of Freedom
#' @param filename_suffix For convenience when you're creating test plots
#'
#' @return Writes a PNG file into /media/work/gwis/posthoc/exposure
#' @export
#'
#' @examples create_qqplot(x, 'aspirin', c('sex', 'age'), 'chiSqGxE', df = 1, "_test")
create_qqplot <- function(dat, exposure, covars, stat, df, filename_suffix = "", output_dir) {

  # calculate p value + lambda
  pvals <- pchisq(dat[, stat], df = df, lower.tail = F)
  lambda <- round( (median(qchisq(1-pvals, df)) / qchisq(0.5, df)), 4)

  # calculate lambda1000
  cases <- unique(dat[, 'Cases'])
  controls <- unique(dat[, 'Subjects']) - unique(dat[, 'Cases'])
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000)))

  # plotting function
  png(paste0(output_dir, "qq_", stat, "_", exposure, filename_suffix, ".png"), height = 720, width = 1280)
  qqman::qq(pvals,
            xlab = "Expected -log10(p)",
            ylab = "Observed -log10(p)",
            main = write_plot_title(stat, exposure, covars, total),
            cex.main = 1.6,
            cex.axis = 1.3,
            cex.lab = 1.3,
            cex.sub = 1.3,
            cex = 1.4,
            pch = 1,
            col = 'blue4')
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3) # FYI ~~ adds spaces when using signif
  dev.off()
}



#' create_qqplot_ge
#'
#' Function creates a QQ plot and outputs a PNG file. Using a separate function for GE, GE among Controls, and GE among Cases (Case-only analysis) just for cosmetic reasons and to omit lambda1000
#'
#' @param dat Input data
#' @param exposure Exposure variable
#' @param stat GxEScan results chi-square statistic to plot
#' @param df Degrees of Freedom
#' @param filename_suffix For convenience when you're creating test plots
#'
#' @return Writes a .png file into location /media/work/gwis/posthoc/exposure
#' @export
#'
#' @examples create_qqplot(x, 'aspirin', c('sex', 'age'), 'chiSqGxE', df = 1, "_test")
create_qqplot_ge <- function(dat, exposure, covars, stat, df, filename_suffix = "") {

  # calculate p value + lambda
  # for control/case only, just report lambda
  pvals <- calculate_pval(dat[, get(stat)], df)
  lambda <- round( (median(qchisq(1-pvals, df)) / qchisq(0.5, df)), 4)

  # get totals for plot title
  cases <- unique(dat[, 'Cases'])
  controls <- unique(dat[, 'Subjects']) - unique(dat[, 'Cases'])

  if (stat == 'chiSqControl') {
    total = controls
    } else if (stat == "chiSqCase") {
    total = cases
    } else {
    total = cases + controls
    }

   # plotting function
  png(paste0("/media/work/gwis/posthoc/", exposure, "/qq_", stat, "_", exposure, filename_suffix, ".png"), height = 720, width = 1280)
  qqman::qq(pvals,
            xlab = "Expected -log10(p)",
            ylab = "Observed -log10(p)",
            main = write_plot_title(stat, exposure, covars, total),
            cex.main = 1.7,
            cex.axis = 1.7,
            cex.lab = 1.7,
            cex.sub = 1.7,
            cex = 1.4,
            pch = 1,
            col = 'blue4')
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4))), cex.sub = 1.7) # FYI ~~ adds spaces when using signif
  dev.off()
}



#' create_manhattanplot
#'
#' Function to create manhattanplot. Uses package EasyStrata - need to write 2 files: results flat text (with appropriate columns), and .ecf file which contains EasyStrata parameters.
#'
#' NOTE: using hardcoded paths to store temporary results text files, and to read annotation file. Be careful!
#'
#' @param dat Input data
#' @param exposure Exposure variable (string)
#' @param stat GxEScan results chi-square statistic
#' @param df Degrees of freedom
#' @param annotation_file File to be used for annotation. Should located in /media/work/tmp
#' @param filename_suffix For convenience when you're creating test plots
#'
#' @return Writes a .png file into location /media/work/gwis/posthoc/exposure
#' @export
#'
#' @examples create_manhattanplot(x, 'aspirin', c('sex', 'age'), 'chiSqGxE', df = 1, "_test")

create_manhattanplot <- function(x, exposure, stat, df, annotation_file, output_dir, filename_suffix = "") {
  
  # calculate p value and filter
  if (stat == "chiSqGxE") {
    tmpdf <- x %>% 
      mutate(P = pchisq(x[, stat], df = df, lower.tail = F)) %>% 
      filter(!(P > 5e-6))
  } else {
    tmpdf <- x %>% 
      mutate(P = pchisq(x[, stat], df = df, lower.tail = F)) %>% 
      filter(!(P > 5e-8))
  }
  
  # output ata.frame of significant results
  saveRDS(tmpdf , file = glue("{output_dir}manhattan_{stat}_{exposure}{filename_suffix}_df.rds"))
  
  # create temporary data set for plotting
  x <- x %>%
    mutate(P = calculate_pval(.data[[stat]], df)) %>%
    filter(!(P > 0.01)) %>%
    dplyr::rename(CHR = Chromosome,
                  BP = Location) %>%
    dplyr::select(SNP, CHR, BP, P)
  
  write.table(x, file = paste0("/media/work/tmp/manhattan_", stat, "_", exposure, filename_suffix), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0(output_dir)
  ecf2 <- "SNP;CHR;BP;P"
  ecf3 <- "character;numeric;numeric;numeric"
  ecf4 <- paste0("/media/work/tmp/manhattan_", stat, "_", exposure, filename_suffix)
  ecf_file_name <- paste0(output_dir,"EasyStrata_", stat, "_", exposure, filename_suffix, ".ecf")
  
  cat(paste0("DEFINE	--pathOut ", ecf1, "
      --acolIn ", ecf2, "
      --acolInClasses ", ecf3, "
      --strMissing NA
      --strSeparator TAB

      EASYIN	--fileIn ", ecf4, "

      START EASYX

      ################
      ## MHplot
      ################

      MHPLOT
      --colMHPlot P
      --colInChr CHR
      --colInPos BP
      --numWidth 1280
      --numHeight 720
      #--astrDefaultColourChr gray51;gray66
      --astrDefaultColourChr gray70;gray80
      --blnYAxisBreak 1
      --numYAxisBreak 22
      #--numPvalOffset 0.01
      # Annotation
      --fileAnnot /home/rak/data/Annotations/", annotation_file, "
      --numAnnotPosLim 1
      # Horizontal lines
      --anumAddPvalLine 5e-6;5e-8
      --anumAddPvalLineLty 6;6
      --astrAddPvalLineCol blue;red
      # Other Graphical Params
      --anumParMar 7;5;7;3
      --numDefaultSymbol 1
      --numDefaultCex 0.6
      --numCexAxis 1.7
      --numCexLab 1.7
      --arcdSymbolCrit  (P>5e-8 & P<5e-6);P<5e-8
      --anumSymbol 20;19
      --arcdColourCrit (P>5e-8 & P<5e-6);P<5e-8
      --astrColour gray55;black
      --arcdCexCrit (P>5e-8 & P<5e-6);P<5e-8
      --anumCex 0.8;0.6

      STOP EASYX"), file = ecf_file_name, append = F)
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}





#-----------------------------------------------------------------------------#
# Weighted Hypothesis Testing for 2-step methods ------
# kooperberg, murcray, edge
#-----------------------------------------------------------------------------#

#' format_twostep_data
#'
#' Function to format step 1 (arrange by p value) --- assume data is the gxe object from GxEScanR
#'
#' @param dat Input data
#' @param stats_step1 Step1 chiSq statistic
#' @param sizeBin0 Size of initial bin (this value needs to be power optimized)
#' @param alpha Overall alpha level
#'
#' @return A DATA.TABLE (!) - mainly because the code I copied pasted used a dt. Contains bin number, bin logp threshold, normalized X axis information for plotting.
#' @export
#' @import data.table
#'
#' @examples format_data_twostep_data(gxe, 'chiSqG', 5, 0.05)
format_twostep_data <- function(dat, stats_step1, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
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
  
  
  # format output for plotting..
  m = nrow(pv)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  # add grp, wt (p value threshold for each bin), rename p value to 'y'
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
  pv[,grp:=grp] # assigning group to the p values..
  setkey(pv,grp)
  for(i in 1:max(grp))
  {
    pv[J(i),wt:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
  }
  
  # return the data.table
  return(pv)
}


#' create_twostep_weighted_plot
#'
#' Create weighted hypothesis testing plot. First step is creating a list object based on the data.table object output from "format_data_twostep_data" function. Second step is actual plotting. Remember that it only plots the first 10 bins for legibility.
#'
#' @param dat Input data
#' @param sizeBin0 Initial bin size (seeds size of remaining bins)
#' @param alpha Overall alpha value (default should be 0.05)
#' @param binsToPlot Number of bins to include in plot
#' @param statistic Step 1 filtering chi-square statistic (depends on method desired)
#'
#' @return A weighted hypothesis plot (png file)
#' @export
#'
#' @examples create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')
create_twostep_weighted_plot <- function(dat, exposure, covars, sizeBin0, alpha, binsToPlot, statistic, filename_suffix = "", output_dir) {
  
  cases <- unique(data.frame(dat[, 'Cases']))
  controls <- unique(data.frame(dat[, 'Subjects'])) - unique(data.frame(dat[, 'Cases']))
  total <- cases + controls
  
  # some plot title and file names
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }
  
  # write_twostep_weightedHT_plot_title("chiSqG", exposure, covars, total)
  
  # bin information
  m = nrow(dat)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} # add +1 bin if condition met
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha level for each bbin 1,2, ... N bin tested
  
  # create list and vars for plotting
  min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
  last.sig = alphaBin[binsToPlot]
  
  
  # create list where each component contains a bin
  # log transform bin alpha value
  # create 'x', normalized position information for each bin
  
  dat <- data.table(dat)
  
  glist<-list()
  for(i in 1:binsToPlot){
    tmp <- dat[J(i)][order(Chromosome, Location)]
    tmp[, ref := -1*log10(min(tmp[,wt]))] # -log10 of bin specific alpha
    tmp[, mapinfo_tmp := seq(Chromosome, Chromosome+1, length.out = .N), by = Chromosome]
    tmp[, x := 0.8*((tmp[,mapinfo_tmp]-min(tmp[,mapinfo_tmp])) / (max(tmp[,mapinfo_tmp])-min(tmp[,mapinfo_tmp]))) + 0.1 + i - 1]
    glist[[i]]<-tmp[order(Chromosome, Location)]
    rm(tmp)
  }
  
  significant_hits <- dat[step2p <= wt]
  saveRDS(significant_hits, file = paste0(output_dir, "/twostep_wht_", statistic, "_", exposure, filename_suffix, "_df.rds"))
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(paste0(output_dir, "/twostep_wht_", statistic, "_", exposure, filename_suffix, ".png"), height = 720, width = 1280)
  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  plot(glist[[1]][,x], glist[[1]][,y],
       col = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], '#E41A1C','#377EB8'),
       pch = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 19, 20),
       cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0,binsToPlot),
       ylim=c(0,min.p),
       axes=F,
       cex.main = 1.7,
       cex.axis = 1.7,
       cex.lab = 1.7,
       cex.sub = 1.7)
  # cex = 1.2)
  lines(glist[[1]][,x], glist[[1]][,ref], col = "black", lwd=1)
  
  # the rest of the points
  for(i in 2:binsToPlot){
    points(glist[[i]][,x], glist[[i]][,y],
           col = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, '#E41A1C', color[i]),
           pch = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, 19, 20),
           cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
           cex.main = 1.7,
           cex.axis = 1.7,
           cex.lab = 1.7,
           cex.sub = 1.7)
    # cex = 1.2)
    lines(glist[[i]][,x], glist[[i]][,ref], col = "black",lwd = 1)
  }
  
  ## the last bin..
  ## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
  ## let's only plot the first 15 bins for now, so change code a bit above.
  # points(glist[[num]][,x], glist[[num]][,y],
  #        col= color[num], pch = 19, cex = 0.5)
  # lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
}










#' format_twostep_data_expectation
#'
#' Function to format step 1 (arrange by p value) --- assume data is the gxe object from GxEScanR
#'
#' @param dat Input data
#' @param stats_step1 Step1 chiSq statistic
#' @param sizeBin0 Size of initial bin (this value needs to be power optimized)
#' @param alpha Overall alpha level
#'
#' @return A DATA.TABLE (!) - mainly because the code I copied pasted used a dt. Contains bin number, bin logp threshold, normalized X axis information for plotting.
#' @export
#' @import data.table
#'
#' @examples format_twostep_data_expectation(gxe, 'chiSqG', 5, 0.05)
format_twostep_data_expectation <- function(dat, stats_step1, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
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
  
  # format output for plotting..
  # m = nrow(pv)
  m = 1000000
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  
  ## sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin_step1 = endpointsBin/1000000
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  
  # expectation based grouping of step 1 pvalues
  alphaBinCut <- c(-Inf, alphaBin_step1, Inf)
  pv[ , grp:=as.numeric(cut(pv[,step1p], breaks = alphaBinCut )) ]
  rep_helper <- c(table(pv[,grp]))
  pv[ , wt:=rep(alphaBin, rep_helper)]
  
  # return the data.table
  return(pv)
}

#' create_twostep_weighted_plot_expectation
#'
#' Create weighted hypothesis testing plot. First step is creating a list object based on the data.table object output from "format_data_twostep_data" function. Second step is actual plotting. Remember that it only plots the first 10 bins for legibility.
#'
#' @param dat Input data
#' @param sizeBin0 Initial bin size (seeds size of remaining bins)
#' @param alpha Overall alpha value (default should be 0.05)
#' @param binsToPlot Number of bins to include in plot
#' @param statistic Step 1 filtering chi-square statistic (depends on method desired)
#'
#' @return A weighted hypothesis plot (png file)
#' @export
#'
#' @examples create_twostep_weighted_plot_expectation(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')
create_twostep_weighted_plot_expectation <- function(dat, exposure, covars, sizeBin0, alpha, binsToPlot, statistic, filename_suffix = "") {
  
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
  
  # get number of bins and bin sizes based on expectation
  m = 1000000
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin_step1 = endpointsBin/1000000
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  # create list of bin specific data.tables
  logp_plot_limit = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
  last.sig = alphaBin[binsToPlot]
  
  glist<-list()
  setkey(dat,grp)
  for(i in 1:binsToPlot){
    z <- dat[J(i)][order(Chromosome, Location)]
    z[, log_binalpha := -1*log10(min(z[,wt]))] # -log10 of bin specific alpha
    z[, logstep2p := -log10(step2p)]
    z[, mapinfo_tmp := seq(Chromosome, Chromosome+1, length.out = .N), by = Chromosome]
    z[, mapinfo := 0.9*( (z[,mapinfo_tmp]-min(z[,mapinfo_tmp])) / (max(z[,mapinfo_tmp])-min(z[,mapinfo_tmp])) ) + 0.05 + i - 1]
    glist[[i]] <- z[order(Chromosome, Location)]
    rm(z)
  }
  
  significant_hits <- dat[step2p <= wt]
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(paste0("/media/work/gwis/posthoc/", exposure, "/twostep_wht_", statistic, "_", exposure,  filename_suffix, ".png"), height = 720, width = 1280)
  
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
  }
  
  ## the last bin..
  ## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
  ## let's only plot the first 15 bins for now, so change code a bit above.
  # points(glist[[num]][,x], glist[[num]][,y],
  #        col= color[num], pch = 19, cex = 0.5)
  # lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
  
}



















#' plot_twostep
#'
#' @param x 
#' @param exposure 
#' @param covars 
#' @param binsToPlot 
#' @param stats_step1 
#' @param sizeBin0 
#' @param alpha 
#' @param output_dir 
#' @param filename_suffix 
#'
#' @return
#' @export
#'
#' @examples
plot_twostep <- function(x, exposure, covars, binsToPlot, stats_step1, sizeBin0, alpha, output_dir, filename_suffix = "") {
  
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
  create_mapinfo <- function(x) {
    x %>% 
      group_by(Chromosome) %>% 
      mutate(mapinfo_tmp = seq(unique(Chromosome), unique(Chromosome) + 1, length.out = n())) %>% 
      ungroup() %>% 
      mutate(mapinfo = 0.9*( (mapinfo_tmp-min(mapinfo_tmp)) / (max(mapinfo_tmp)-min(mapinfo_tmp))) + 0.05 + bin_number - 1)
  }
  
  
  # ------- create working dataset ------- #
  # assign SNPs to bins
  m = nrow(x) 
  nbins = floor(log2(m/sizeBin0 + 1)) 
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) 
  endpointsBin = cumsum(sizeBin)
  
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2))
  
  rep_helper <- c(table(grp))
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin 
  alphaBin_dat <- rep(alphaBin, rep_helper)
  
  
  # working data.frame
  tmp <- x %>% 
    mutate(step1p = case_when(stats_step1 == "chiSqEDGE" ~ pchisq(.data[[stats_step1]], df = 2, lower.tail = F),
                              TRUE ~ pchisq(.data[[stats_step1]], df = 1, lower.tail = F)), 
           step2p = pchisq(.data[['chiSqGxE']], df = 1, lower.tail = F)) %>% 
    arrange(step1p) %>% 
    mutate(bin_number = as.numeric(grp), 
           step2p_sig = as.numeric(alphaBin_dat), 
           log_step2p_sig = -log10(step2p_sig), 
           log_step2p = -log10(step2p))
  
  significant_hits <- filter(tmp, step2p < step2p_sig)
  
  ## output list of bins for plotting
  tmp_plot <- tmp %>% 
    arrange(Chromosome, Location) %>% 
    group_by(bin_number) %>% 
    group_split()
  
  
  ## add mapinfo (see functions)
  tmp_plot <- map(tmp_plot, create_mapinfo)
  
  cases <- unique(tmp[, 'Cases'])
  controls <- unique(tmp[, 'Subjects']) - unique(tmp[, 'Cases'])
  total <- cases + controls
  logp_plot_limit = 12
  last.sig = alphaBin[binsToPlot]
  
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
  lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig), col = "black", lwd=1)
  
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
    lines(pull(bin_to_plot, mapinfo), pull(bin_to_plot, log_step2p_sig),
          col = "black",lwd = 1)
  }
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.2, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(stats_step1, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
  
  # return data.frame of significant results!
  saveRDS(significant_hits, file = glue(output_dir, "twostep_wht_{stats_step1}_{exposure}{filename_suffix}_df.rds"))
  return(significant_hits)
  
}











#' create_qqplot_ramwas
#' 
#' wrapper function to create qqplot using package ramwas (much faster)
#'
#' @param x 
#' @param stat 
#' @param df 
#' @param filename_suffix 
#' @param output_dir 
#'
#' @return
#' @export
#'
#' @examples
create_qqplot_ramwas <- function(x, stat, df, filename_suffix = "", output_dir) {
  
  # calculate p value
  pvals <- pchisq(x[, stat], df = df, lower.tail = F)
  # remove NAs.. (usually for controls and case only)
  pvals <- pvals[!is.na(pvals)]
  
  lambda <- round( (median(qchisq(1-pvals, df)) / qchisq(0.5, df)), 4)
  
  # calculate lambda1000
  if (stat %in% c("chiSqControl", "chiSqCase")) {
    lambda1000 <- "NA"
  } else {
    cases <- unique(x[, 'Cases'])
    controls <- unique(x[, 'Subjects']) - unique(x[, 'Cases'])
    total <- cases + controls
    cases1000 <- (cases/total) * 1000
    controls1000 <- (controls/total) * 1000
    lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000)))
    lambda1000 <- sprintf("%.3f", lambda1000)
  }

  # plotting function
  png(paste0(output_dir, "qq_", stat, "_", exposure, filename_suffix, ".png"), height = 720, width = 1280)
  qqPlotFast_figifs(pvals, lambda1000 = lambda1000)
  dev.off()
}





#' qqPlotFast_figifs
#' 
#' From ramwas package, changed so it also reports lambda1000
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
qqPlotFast_figifs <- function (x, ntests = NULL, ismlog10 = FALSE, ci.level = 0.05, 
                        ylim = NULL, newplot = TRUE, col = "#D94D4C", cex = 0.5, 
                        yaxmax = NULL, lwd = 3, axistep = 2, col.band = "#ECA538", 
                        makelegend = TRUE, lambda1000, 
                        xlab = expression(paste("–", " log"[10] * "(", italic("P"), "), null")), ylab = expression(paste("–", " log"[10] * "(", italic("P"), "), observed"))) 
{
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




