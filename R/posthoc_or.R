#=============================================================================#
# stratified odds ratios from yi's code
#=============================================================================#


#' Model.all.new
#'
#' create stratified odds ratios. Written by yi lin
#'
#' @param data dataset
#' @param mod string of model to be fitted, exclude dependent var e.g. "age_ref_imp+sex+study_gxe+pc1+pc2+pc3"
#' @param env string of exposure variable e.g. 'asp_ref'
#'
#' @return a list containing model estimates and standard errors for stratified table cells
#' @export
#'
#' @examples Model.all.new(figi, "age_ref_imp+sex+study_gxe+pc1+pc2+pc3", "asp_ref")
Model.all.new <- function(data,mod,env){
  if(length(table(data$p1))<=1 | length(table(data$p2))<=1 | 
     var(2*data$p2+data$p1,na.rm=T)==0 | 
     sum((data$p2+data$p1)>1.1,na.rm=T)>0) {
    res = NA
  } else {
    data[,env]=factor(data[,env])
    mModel <- paste0('outcome~',env,'*p1+',env,'*p2+',mod)
    tmp <- summary(glm(mModel, family='binomial', data=data))
    COV <- tmp$cov.unscaled
    env.idx = env
    if(is.factor(data[,env])) env.idx = paste0(env,levels(data[,env])[-1])
    
    # stratified by GE  
    eg0 <- tmp$coef[env.idx,c(1,2)]
    beta.eg1 <- tmp$coef[env.idx,1]+tmp$coef['p1',1]+tmp$coef[paste0(env.idx,':p1'),1]
    beta.eg2 <- tmp$coef[env.idx,1]+tmp$coef['p2',1]+tmp$coef[paste0(env.idx,':p2'),1]
    e0g <- tmp$coef[c('p1','p2'),c(1,2)]
    I3 =  rep(1,3)
    se.eg1 = se.eg2 = NULL
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs = COV[c(env.idx[i],'p1',paste0(env.idx[i],':p1')),c(env.idx[i],'p1',paste0(env.idx[i],':p1'))]		
      se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
      covs = COV[c(env.idx[i],'p2',paste0(env.idx[i],':p2')),c(env.idx[i],'p2',paste0(env.idx[i],':p2'))]		
      se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
    }
    GE<- c(eg0,e0g,beta.eg1,se.eg1,beta.eg2,se.eg2)
    
    # stratified by G 
    g0 <-tmp$coef[env.idx,1:2]
    if(length(env.idx)==1){
      idx.g1 = c(env.idx,paste0(env.idx,':p1'))
      idx.g2 = c(env.idx,paste0(env.idx,':p2'))
      g1 = sum(tmp$coef[idx.g1,1],na.rm=T)
      g2 = sum(tmp$coef[idx.g2,1],na.rm=T)
      Is1 =  rep(1,length(idx.g1))
      Is2 =  rep(1,length(idx.g2))
    }else{
      est.g1 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p1'),1])
      est.g2 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p2'),1])
      g1 = rowSums(est.g1,na.rm=T)
      g2 = rowSums(est.g2,na.rm=T)
      Is1 =  rep(1,ncol(est.g1))
      Is2 =  rep(1,ncol(est.g2))
    }
    
    se.g1 = se.g2 = NULL
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs1 = COV[c(env.idx[i],paste0(env.idx[i],':p1')),c(env.idx[i],paste0(env.idx[i],':p1'))]  
      covs2 = COV[c(env.idx[i],paste0(env.idx[i],':p2')),c(env.idx[i],paste0(env.idx[i],':p2'))]  
      se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
      se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
    }
    
    G<- c(g0,g1,se.g1,g2,se.g2)
    
    # stratified by E
    
    e0 <-tmp$coef[c('p1','p2'),1:2]
    if(length(env.idx)==1){
      idx.e11 = c('p1',paste0(env.idx,':p1'))
      idx.e12 = c('p2',paste0(env.idx,':p2'))
      e11 = sum(tmp$coef[idx.e11,1],na.rm=T)
      e12 = sum(tmp$coef[idx.e12,1],na.rm=T)
      Is1 =  rep(1,length(idx.e11))
      Is2 =  rep(1,length(idx.e12))
    }else{
      est.e11 <- data.frame(tmp$coef['p1',1],tmp$coef[paste0(env.idx,':p1'),1])
      est.e12 <- data.frame(tmp$coef['p2',1],tmp$coef[paste0(env.idx,':p2'),1])
      e11 = rowSums(est.e11,na.rm=T)
      e12 = rowSums(est.e12,na.rm=T)
      Is1 =  rep(1,ncol(est.e11))
      Is2 =  rep(1,ncol(est.e12))
    }
    
    se.e11 = se.e12 = NULL
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs1 = COV[c('p1',paste0(env.idx[i],':p1')),c('p1',paste0(env.idx[i],':p1'))]  
      covs2 = COV[c('p2',paste0(env.idx[i],':p2')),c('p2',paste0(env.idx[i],':p2'))]  
      se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
      se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
    }
    E <- c(e0[1,1],e11,e0[1,2],se.e11,e0[2,1],e12,e0[2,2],se.e12)
  }
  res<-list(GE=GE,G=G,E=E)
}



#' ORtab
#' 
#' clean up outputs from Model.all.new
#'
#' @param x 
#' @param elvl 
#' @param glvl 
#' @param res 
#'
#' @return
#' @export
#'
#' @examples
ORtab = function(x,elvl,glvl,res) {
  ORs = data.frame(matrix(NA,nrow=length(elvl),ncol=length(glvl)))
  colnames(ORs) = paste0('p',glvl)
  rownames(ORs) = paste0('OR',elvl)
  ORs[1,1] = 1
  for(c in colnames(ORs))
    for(r in rownames(ORs))
      if(!(c=='p0' & r==paste0('OR',elvl[1]))) 
        ORs[r,c] =  as.character(res[res$snp %in% x,paste0(r,c)])
  ORs      
} 


#' ptab
#' 
#' clean up output from Model.all.new
#'
#' @param x 
#' @param elvl 
#' @param glvl 
#' @param res 
#'
#' @return
#' @export
#'
#' @examples
ptab = function(x,elvl,glvl,res) {
  pval = data.frame(matrix(NA,nrow=length(elvl),ncol=length(glvl)))
  colnames(pval) = paste0('p',glvl)
  rownames(pval) = elvl
  pval[1,1] = 1
  for(c in colnames(pval))
    for(r in rownames(pval))
      if(!(c=='p0' & r==elvl[1])) 
        pval[r,c] =  as.numeric(as.vector(res[res$snp %in% x,paste0('pval',r,c)]))
  pval = apply(pval,2,function(y) paste0('P= ',formatC(y,format='g',digit=2)))    
}


#' format_res
#' 
#' clean up output from Model.all.new
#'
#' @param res 
#'
#' @return
#' @export
#'
#' @examples
format_res <- function(res) {
  betas = colnames(res)[grep('beta',colnames(res),fixed=T)]
  ses = colnames(res)[grep('se',colnames(res),fixed=T)]
  ORs = sapply(seq(length(betas)),function(x,betas,ses,res) 
  {
    or = paste0(rnd2(exp(res[,betas[x]])),' (',
                rnd2(exp(res[,betas[x]]-qnorm(0.975)*res[,ses[x]])),'-',
                rnd2(exp(res[,betas[x]]+qnorm(0.975)*res[,ses[x]])),')')
  },betas=betas,ses=ses,res=res)
  if(nrow(res)==1)  ORs = data.frame(t(ORs)) 
  colnames(ORs) = sub('beta','OR',betas)
  pvals = sapply(seq(length(betas)),function(x,betas,ses,res) 
  {pval = 2*pnorm(-abs(res[,betas[x]]/res[,ses[x]]))},betas=betas,ses=ses,res=res)
  if(nrow(res)==1) pvals = data.frame(t(pvals)) 
  colnames(pvals) = sub('beta','pval',betas)
  pvals.p = apply(pvals,2, function(y) paste0('P= ',formatC(y,format='g',digit=2)))
  if(nrow(res)==1) pvals.p = data.frame(t(pvals.p)) 
  colnames(pvals.p) = paste0('P',colnames(pvals.p))
  res = data.frame(res,ORs,pvals,pvals.p,stringsAsFactors=F)
}


#' rnd2
#' 
#' rounds numbers
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
rnd2 <- function(x) {
  round(x, 2)
}










#' Model.all.new.dosage
#' 
#' similar to original function, but creating outputs for dosage odds ratios. this is so you can create stratified odds ratios but modeling E specific dosage associations (which is modeled continuously, meaning it's a per allele Odds Ratios)
#'
#' @param data 
#' @param mod 
#' @param env 
#'
#' @return
#' @export
#'
#' @examples
Model.all.new.dosage <- function(data, mod, env) {
  if(length(table(data$p1))<=1 | length(table(data$p2))<=1 | var(2*data$p2+data$p1,na.rm=T)==0 | sum((data$p2+data$p1)>1.1,na.rm=T)>0) {
    res = NA
  } else {
    data[,env]=factor(data[,env])
    mModel <- paste0('outcome~',env,'*dosage+',mod)
    tmp <- summary(glm(mModel, family='binomial', data=data))
    COV <- tmp$cov.unscaled
    env.idx = env
    if(is.factor(data[,env])) env.idx = paste0(env,levels(data[,env])[-1])
    
    # stratified by GE  
    eg0 <- tmp$coef[env.idx,c(1,2)]
    beta.eg1 <- tmp$coef[env.idx,1]+tmp$coef['dosage',1]+tmp$coef[paste0(env.idx,':dosage'),1]
    # beta.eg2 <- tmp$coef[env.idx,1]+tmp$coef['p2',1]+tmp$coef[paste0(env.idx,':p2'),1]
    e0g <- tmp$coef[c('dosage'),c(1,2)]
    I3 =  rep(1,3)
    se.eg1 = se.eg2 = NULL
    for(i in 1:(nlevels(data[,env])-1)) {
      covs = COV[c(env.idx[i],'dosage',paste0(env.idx[i],':dosage')), c(env.idx[i],'dosage',paste0(env.idx[i],':dosage'))]		
      se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
      # covs = COV[c(env.idx[i],'p2',paste0(env.idx[i],':p2')),c(env.idx[i],'p2',paste0(env.idx[i],':p2'))]		
      # se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
    }
    
    GE <- c(eg0,e0g,beta.eg1,se.eg1)
    
    # stratified by G 
    g0 <-tmp$coef[env.idx,1:2]
    if(length(env.idx)==1) {
      idx.g1 = c(env.idx,paste0(env.idx,':dosage'))
      g1 = sum(tmp$coef[idx.g1,1],na.rm=T)
      Is1 =  rep(1,length(idx.g1))
    } else {
      est.g1 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p1'),1])
      est.g2 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p2'),1])
      g1 = rowSums(est.g1,na.rm=T)
      g2 = rowSums(est.g2,na.rm=T)
      Is1 =  rep(1,ncol(est.g1))
      Is2 =  rep(1,ncol(est.g2))
    }
    
    se.g1 = se.g2 = NULL
    for(i in 1:(nlevels(data[,env])-1)) {
      covs1 = COV[c(env.idx[i],paste0(env.idx[i],':dosage')),c(env.idx[i],paste0(env.idx[i],':dosage'))]  
      # covs2 = COV[c(env.idx[i],paste0(env.idx[i],':p2')),c(env.idx[i],paste0(env.idx[i],':p2'))]  
      se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
      # se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
    }
    
    G <- c(g0,g1,se.g1)
    
    # stratified by E
    
    e0 <-tmp$coef[c('dosage'),1:2]
    if(length(env.idx)==1) {
      idx.e11 = c('dosage',paste0(env.idx,':dosage'))
      # idx.e12 = c('p2',paste0(env.idx,':p2'))
      e11 = sum(tmp$coef[idx.e11,1],na.rm=T)
      # e12 = sum(tmp$coef[idx.e12,1],na.rm=T)
      Is1 =  rep(1,length(idx.e11))
      # Is2 =  rep(1,length(idx.e12))
    } else {
      est.e11 <- data.frame(tmp$coef['p1',1],tmp$coef[paste0(env.idx,':p1'),1])
      est.e12 <- data.frame(tmp$coef['p2',1],tmp$coef[paste0(env.idx,':p2'),1])
      e11 = rowSums(est.e11,na.rm=T)
      e12 = rowSums(est.e12,na.rm=T)
      Is1 =  rep(1,ncol(est.e11))
      Is2 =  rep(1,ncol(est.e12))
    }
    
    se.e11 = se.e12 = NULL
    for(i in 1:(nlevels(data[,env])-1)) {
      covs1 = COV[c('dosage',paste0(env.idx[i],':dosage')),c('dosage',paste0(env.idx[i],':dosage'))]  
      # covs2 = COV[c('p2',paste0(env.idx[i],':p2')),c('p2',paste0(env.idx[i],':p2'))]  
      se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
      # se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
    }
    E <- c(e0[1],e11,e0[2],se.e11)
  }
  res<-list(GE=GE,G=G,E=E)
}




#-----------------------------------------------------------------------------#
# wrapper functions to generate two types of stratified odds ratios tables
# 1) original that Yi wrote
# 2) similar to 1, but includes dosage per allele effect in the G column
#-----------------------------------------------------------------------------#


#' fit_stratified_or
#' 
#' create stratified odds ratios table (code written by yi)
#'
#' @param ds dataset
#' @param exposure string; exposure variable
#' @param snp string; snp name (chr1_bp_ref_alt)
#' @param hrc_version string; e.g. v2.4
#' @param covariates vector; adjustment covarites
#' @param mod string; model of exposure + covariates only. can specify study_gxe interactions if needed (that's how Yi's original code modeled these interaction models)
#' @param dosage whether you want to report dosage or genotype probabilities
#' @param output_dir string output directory
#'
#' @return a matrix of odds ratios, p-values, and strata counts
#' @export
#'
#' @examples fit_stratified_or(figi, 'asp_ref', 'chr5_1234_A_C', 'v2.4', c("age_ref_imp", "sex"), "age_ref_imp+sex")
fit_stratified_or <- function(ds, exposure, snp, hrc_version, covariates, mod, dosage = F, output_dir) {
  
  # output directory
  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = "_"))
  
  # check allele frequency
  total_dosage <- sum(ds[,snp])
  total_alleles <- nrow(ds) * 2
  aaf <- total_dosage / total_alleles
  
  # ---- recode SNPs so that major allele is the reference group
  if(aaf >= 0.5) {
    # flip dosages
    ds[, snp] <- abs(ds[, paste0(snp)] - 2)
    ds[, paste0(snp, '_dose')] <- abs(ds[, paste0(snp, '_dose')] - 2)
    # flip genotype probabilities
    pp <- ds[,paste0(snp, "_p2")]
    ds[,paste0(snp, "_p2")] <- ds[, paste0(snp, "_p0")]
    ds[,paste0(snp, "_p0")] <- pp
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }
  
  
  # create data subset
  tmp1 = ds[, c('outcome', exposure, covariates)]
  # tmp2 = ds[, grepl(snp, names(ds))]
  tmp2 = ds[, grepl(paste(paste0(snp, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(ds))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
  
  ds_tmp = cbind(tmp1, tmp2) %>% 
    na.omit(.[, c(exposure,'outcome', covariates, 'p1','p2','dosage')])
  
  res.pool.un = res.pool.g = res.pool.e = NULL
  Ncaco = data.frame(snps=snp,matrix(NA,ncol=6*2))
  colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                      paste0('Co1',1:2),paste0('Ca1',1:2),
                      paste0('Co2',1:2),paste0('Ca2',1:2))
  rownames(Ncaco) = snp
  
  
  #---- Calculate counts for each cell ----------------
  Ncaco[snp,c(paste0('Co',1:2),paste0('Ca',1:2))] = t(table(ds_tmp[,c('outcome',exposure)]))
  Ncaco[snp,c(paste0('Co1',1:2),paste0('Ca1',1:2))] = c(t(tapply(ds_tmp$p1, ds_tmp[,c('outcome',exposure)],sum,na.rm=T)))
  Ncaco[snp,c(paste0('Co2',1:2),paste0('Ca2',1:2))] = c(t(tapply(ds_tmp$p2, ds_tmp[,c('outcome',exposure)],sum,na.rm=T)))
  
  #-- Fit unrestricted model --------
  ds_tmp[,exposure] = factor(ds_tmp[,exposure])
  tmp = as.vector(Model.all.new(ds_tmp,mod,exposure))
  res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))
  
  ## organize results ######
  elvl=c(0:1) ; glvl = c(0,1,2)
  
  colnames(res.pool.un) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                            'beta0p1','beta0p2','se0p1','se0p2',
                            paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'),
                            paste0('beta',elvl[-1],'p2'),paste0('se',elvl[-1],'p2'))
  res.pool.un = format_res(res.pool.un)
  
  
  ##== stratified by G results #######
  colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1','beta2','se2'),each=1),
                                              rep(elvl[-1],6)))
  res.pool.g = format_res(res.pool.g)
  
  ##== stratified by E results #######
  colnames(res.pool.e) = c('snp','env',paste0('beta1',elvl),paste0('se1',elvl),
                           paste0('beta2',elvl),paste0('se2',elvl))
  res.pool.e = format_res(res.pool.e)
  
  ##== Put into table 
  OR.tab = ORtab(snp,elvl=elvl,glvl=glvl,res=res.pool.un)
  pval.tab = ptab(snp,elvl=elvl,glvl=glvl,res=res.pool.un)
  
  ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('OR0',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR1',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR2',elvl[-1])])))), ncol=3)
  pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('Ppval0',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval1',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval2',elvl[-1])])))),ncol=3)
  
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))
  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)],
                                         res.pool.e[,paste0('OR2',elvl)]))),ncol=2)
  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)],
                                         res.pool.e[,paste0('Ppval2',elvl)]))),ncol=2,byrow=T)
  
  #== calculate counts for G=0 and put counts into format ca/co
  for(i in 1:2){
    eval(parse(text=paste0('Ncaco$caco0',i,"=paste0(round(Ncaco$Ca",i,'-Ncaco$Ca1',i,'-Ncaco$Ca2',i,
                           ",1),'/',round(Ncaco$Co",i,'-Ncaco$Co1',i,'-Ncaco$Co2',i,',1))')))
    eval(parse(text=paste0('Ncaco$caco1',i,'=paste0(round(Ncaco$Ca1',i,",1),'/',round(Ncaco$Co1",i,",1))")))
    eval(parse(text=paste0('Ncaco$caco2',i,'=paste0(round(Ncaco$Ca2',i,",1),'/',round(Ncaco$Co2",i,",1))")))
  }
  
  #=== write into table
  est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
                    ORg.tab[1,],pg.tab[1,]),
              rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)),
              N0 = c(Ncaco[snp,'caco01'],NA,Ncaco[snp,'caco02'],rep(NA,3)),
              N1 = c(Ncaco[snp,'caco11'],NA,Ncaco[snp,'caco12'],rep(NA,3)),
              N2 = c(Ncaco[snp,'caco21'],NA,Ncaco[snp,'caco22'],rep(NA,3)))
  
  # slightly modify output in case anyone wants per allele effects rather than p1/p2
  if(dosage == T) {
    
    tmp = as.vector(Model.all.new.dosage(ds_tmp,mod,exposure))
    
    res.pool.un = res.pool.g = res.pool.e = NULL
    res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
    # res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
    # res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))
    
    elvl=c(0:1) ; glvl = c(0,1)
    tmp2 = NULL
    tmp2 = rbind(tmp2,data.frame(snp,exposure,t(tmp$GE)))
    
    colnames(tmp2) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                       'beta0p1','beta0p2','se0p1','se0p2')
    tmp2 = format_res(tmp2)
    
    res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
    colnames(res.pool.e) = c('snp',
                             'env',
                             paste0('beta1',elvl),paste0('se1',elvl))
    tmp2 <- format_res(res.pool.e)
    
    
    ORe.tab = matrix(as.character(unlist(c(tmp2[,paste0('OR1',elvl)]))),ncol=2)
    pe.tab  = matrix(as.character(unlist(c(tmp2[,paste0('Ppval1',elvl)]))),ncol=2,byrow=T)
    est2 <- rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2], NA, NA)
    
    final_out <- est %>% 
      dplyr::select(-c(4,5)) %>% 
      add_column(est2, .after = 3)
    
    colnames(final_out) <- c(paste0(ref_allele, ref_allele), 
                             paste0(ref_allele, alt_allele), 
                             paste0(alt_allele, alt_allele), 
                             paste0(alt_allele, " Allelic Dosage"), 
                             "N0", "N1", "N2")
    # this is only for FACTOR variables (might have to modify when you run Q4 variables)
    exposure_level <- levels(ds[,exposure])
    rownames(final_out) <- c(paste0(exposure,"=",exposure_level[1]), 
                             "p0",
                             paste0(exposure,"=",exposure_level[2]), 
                             "p1", 
                             exposure, 
                             "p")
    saveRDS(final_out, file = paste0(output_dir, "stratified_oddsratio_", snp, "_", exposure, "_dosage.rds"))
  } else {
    colnames(est) <- c(paste0(ref_allele, ref_allele), 
                       paste0(ref_allele, alt_allele), 
                       paste0(alt_allele, alt_allele), 
                       paste0(ref_allele, alt_allele), 
                       paste0(alt_allele, alt_allele),
                       "N0", "N1", "N2")
    exposure_level <- levels(ds[,exposure])
    rownames(est) <- c(paste0(exposure,"=",exposure_level[1]), 
                       "p0",
                       paste0(exposure,"=",exposure_level[2]), 
                       "p1", 
                       exposure, 
                       "p")
    
    saveRDS(est, file = paste0(output_dir, "stratified_oddsratio_", snp, "_", exposure, ".rds"))
  }
  
}





#' fit_stratified_or_q4
#'
#' use this one for 
#'
#' @param ds 
#' @param exposure 
#' @param snp 
#' @param hrc_version 
#' @param covariates 
#' @param mod 
#' @param dosage 
#' @param output_dir 
#'
#' @return
#' @export
#'
#' @examples
fit_stratified_or_q4 <- function(ds, exposure, snp, hrc_version, covariates, mod, dosage = F, output_dir) {
  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  # snp_info <- unlist(strsplit(snp, split = "_"))
  # 
  # # ---- check D|G association, recode SNPs as necessary ---- #
  # dg_model <- lm(paste0("outcome ~ ", paste0(snp, "_dose"), "*", exposure, "+", paste0(covariates, collapse = "+")), data = ds)
  # if(dg_model$coefficients[2] < 0) {
  #   # flip dosages
  #   ds[, paste0(snp)] <- abs(ds[, paste0(snp)] - 2)
  #   ds[, paste0(snp, '_dose')] <- abs(ds[, paste0(snp, '_dose')] - 2)
  #   # flip genotype probabilities
  #   pp <- ds[,paste0(snp, "_p2")]
  #   ds[,paste0(snp, "_p2")] <- ds[, paste0(snp, "_p0")]
  #   ds[,paste0(snp, "_p0")] <- pp
  # }
  # 
  # # risk allele
  # if(dg_model$coefficients[2] < 0) {
  #   risk_allele = snp_info[3]
  # } else {
  #   risk_allele = snp_info[4]
  # }
  
  # output directory
  # output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = "_"))
  
  # check allele frequency
  total_dosage <- sum(ds[,snp])
  total_alleles <- nrow(ds) * 2
  aaf <- total_dosage / total_alleles
  
  # ---- recode SNPs so that major allele is the reference group
  if(aaf >= 0.5) {
    # flip dosages
    ds[, snp] <- abs(ds[, paste0(snp)] - 2)
    ds[, paste0(snp, '_dose')] <- abs(ds[, paste0(snp, '_dose')] - 2)
    # flip genotype probabilities
    pp <- ds[,paste0(snp, "_p2")]
    ds[,paste0(snp, "_p2")] <- ds[, paste0(snp, "_p0")]
    ds[,paste0(snp, "_p0")] <- pp
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }
  
  # create data subset
  tmp1 = ds[, c('outcome', exposure, covariates)]
  # tmp2 = ds[, grepl(snp, names(ds))]
  tmp2 = ds[, grepl(paste(paste0(snp, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(ds))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
  
  ds_tmp = cbind(tmp1, tmp2) %>% 
    na.omit(.[, c(exposure,'outcome', covariates, 'p1','p2','dosage')])
  
  
  
  
  
  
  # start here.... #
  
  exposure_levels <- seq(1,4)
  res.pool.un = res.pool.g = res.pool.e = NULL
  Ncaco = data.frame(snps=snp,matrix(NA,ncol=6*4))
  colnames(Ncaco) = c('snps',paste0('Co',exposure_levels),paste0('Ca',exposure_levels),
                      paste0('Co1',exposure_levels),paste0('Ca1',exposure_levels),
                      paste0('Co2',exposure_levels),paste0('Ca2',exposure_levels))
  rownames(Ncaco) = snp
  
  # res.pool.un = res.pool.g = res.pool.e = NULL
  # Ncaco = data.frame(snps=snp,matrix(NA,ncol=6*2))
  # colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
  #                     paste0('Co1',1:2),paste0('Ca1',1:2),
  #                     paste0('Co2',1:2),paste0('Ca2',1:2))
  # rownames(Ncaco) = snp
  
  
  #---- Calculate counts for each cell ----------------
  Ncaco[snp,c(paste0('Co',1:4),paste0('Ca',1:4))] = t(table(ds_tmp[,c('outcome',exposure)]))
  Ncaco[snp,c(paste0('Co1',1:4),paste0('Ca1',1:4))] = c(t(tapply(ds_tmp$p1, ds_tmp[,c('outcome',exposure)],sum,na.rm=T)))
  Ncaco[snp,c(paste0('Co2',1:4),paste0('Ca2',1:4))] = c(t(tapply(ds_tmp$p2, ds_tmp[,c('outcome',exposure)],sum,na.rm=T)))
  
  
  
  #-- Fit unrestricted model --------
  ds_tmp[,exposure] = factor(ds_tmp[,exposure])
  tmp = as.vector(Model.all.new(ds_tmp,mod,exposure))
  res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))
  
  ## organize results ######
  elvl=c(0:1) ; glvl = c(0,1,2)
  elvl=c(0:3) ; glvl = c(0,1,2)
  
  # elvl=c(0:3) ; glvl = c(0,1,2)
  

  colnames(res.pool.un) = c('snp','env',
                            paste0('beta',elvl[-1],'p0'),
                            paste0('se',elvl[-1],'p0'),
                            'beta0p1','beta0p2','se0p1','se0p2',
                            paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'),
                            paste0('beta',elvl[-1],'p2'),paste0('se',elvl[-1],'p2'))
  res.pool.un = format_res(res.pool.un)
  
  
  ##== stratified by G results #######
  # colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1','beta2','se2'),each=1),
  #                                             rep(elvl[-1],6)))
  colnames(res.pool.g) = c('snp','env',
                           paste0('beta',elvl[-1], 'p0'), 
                           paste0('se',elvl[-1], 'p0'),
                           paste0('beta',elvl[-1],'p1'),
                           paste0('se',elvl[-1],'p1'),
                           paste0('beta',elvl[-1],'p2'),
                           paste0('se',elvl[-1],'p2'))
  res.pool.g = format_res(res.pool.g)
  
  ##== stratified by E results #######
  colnames(res.pool.e) = c('snp','env',
                           paste0('beta1',elvl),
                           paste0('se1',elvl),
                           paste0('beta2',elvl),
                           paste0('se2',elvl))
  res.pool.e = format_res(res.pool.e)
  
  ##== Put into table 
  OR.tab = ORtab(snp,elvl=elvl,glvl=glvl,res=res.pool.un)
  pval.tab = ptab(snp,elvl=elvl,glvl=glvl,res=res.pool.un)
  
  ORg.tab = matrix(as.character(unlist(c(as.vector(res.pool.g[,paste0('OR',elvl[-1],'p0')]),
                                         as.vector(res.pool.g[,paste0('OR',elvl[-1],'p1')]),
                                         as.vector(res.pool.g[,paste0('OR',elvl[-1],'p2')])))), ncol=3)
  pg.tab = matrix(as.character(unlist(c(as.vector(res.pool.g[,paste0('Ppval',elvl[-1],'p0')]),
                                        as.vector(res.pool.g[,paste0('Ppval',elvl[-1],'p1')]),
                                        as.vector(res.pool.g[,paste0('Ppval',elvl[-1],'p2')])))),ncol=3)
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))
  
  
  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)],
                                         res.pool.e[,paste0('OR2',elvl)]))),ncol=2)
  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)],
                                         res.pool.e[,paste0('Ppval2',elvl)]))),ncol=2,byrow=T)
  
  #== calculate counts for G=0 and put counts into format ca/co
  for(i in 1:4){
    eval(parse(text=paste0('Ncaco$caco0',i,"=paste0(round(Ncaco$Ca",i,'-Ncaco$Ca1',i,'-Ncaco$Ca2',i,
                           ",1),'/',round(Ncaco$Co",i,'-Ncaco$Co1',i,'-Ncaco$Co2',i,',1))')))
    eval(parse(text=paste0('Ncaco$caco1',i,'=paste0(round(Ncaco$Ca1',i,",1),'/',round(Ncaco$Co1",i,",1))")))
    eval(parse(text=paste0('Ncaco$caco2',i,'=paste0(round(Ncaco$Ca2',i,",1),'/',round(Ncaco$Co2",i,",1))")))
  }
  
  #=== write into table
  # est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
  #                   ORg.tab[1,],pg.tab[1,]),
  #             rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)),
  #             N0 = c(Ncaco[snp,'caco01'],NA,Ncaco[snp,'caco02'],rep(NA,3)),
  #             N1 = c(Ncaco[snp,'caco11'],NA,Ncaco[snp,'caco12'],rep(NA,3)),
  #             N2 = c(Ncaco[snp,'caco21'],NA,Ncaco[snp,'caco22'],rep(NA,3)))
  # 
  
  est = cbind(rbind(OR.tab[1,],pval.tab[1,],
                    OR.tab[2,],pval.tab[2,],
                    OR.tab[3,],pval.tab[3,],
                    OR.tab[4,],pval.tab[4,],
                    ORg.tab[1,],pg.tab[1,],
                    ORg.tab[2,],pg.tab[2,],
                    ORg.tab[3,],pg.tab[3,]),
              rbind(ORe.tab[1,],pe.tab[1,],
                    ORe.tab[2,],pe.tab[2,],
                    ORe.tab[3,],pe.tab[3,],
                    ORe.tab[4,],pe.tab[4,],
                    rep(NA,2),rep(NA,2),rep(NA,2),rep(NA,2),rep(NA,2),rep(NA,2)),
              N0 = c(Ncaco[snp,'caco01'],NA,Ncaco[snp,'caco02'],NA,Ncaco[snp,'caco03'],NA,Ncaco[snp,'caco04'],rep(NA,7)),
              N1 = c(Ncaco[snp,'caco11'],NA,Ncaco[snp,'caco12'],NA,Ncaco[snp,'caco13'],NA,Ncaco[snp,'caco14'],rep(NA,7)),
              N2 = c(Ncaco[snp,'caco21'],NA,Ncaco[snp,'caco22'],NA,Ncaco[snp,'caco23'],NA,Ncaco[snp,'caco24'],rep(NA,7)))
  
  # colnames(est) <- c("G=0", "G=1", "G=2", "G=1", "G=2", "N0", "N1", "N2")
  # rownames(est) <- c("E=0", "p0", "E=1", "p1", "E=2", "p2", "E=3", "p3", "E1", "P01", "E2", "P02", "E3", "P03")
  
  colnames(est) <- c(paste0(ref_allele, ref_allele), 
                     paste0(ref_allele, alt_allele), 
                     paste0(alt_allele, alt_allele), 
                     paste0(ref_allele, alt_allele), 
                     paste0(alt_allele, alt_allele),
                     "N0", "N1", "N2")
  exposure_level <- levels(as.factor(ds[,exposure]))
  rownames(est) <- c(paste0(exposure,"=",exposure_level[1]), 
                     "p0",
                     paste0(exposure,"=",exposure_level[2]), 
                     "p1", 
                     paste0(exposure,"=",exposure_level[3]), 
                     "p2",
                     paste0(exposure,"=",exposure_level[4]), 
                     "p3",
                     paste0(exposure,"=",exposure_level[2], "(G)"), 
                     "p1 (G)", 
                     paste0(exposure,"=",exposure_level[3], "(G)"), 
                     "p2 (G)",
                     paste0(exposure,"=",exposure_level[4], "(G)"), 
                     "p3 (G)")
  
  
  
  # options(knitr.kable.NA = '')
  # 
  # kable(est) %>%
  #   kable_styling('bordered', bootstrap_options = c("striped", "hover", "condensed"), full_width = F, position = 'left') %>% 
  #   add_header_above(c(" " = 4, "G param by E" = 2, "Counts (Ca/Co)" = 3)) %>% 
  #   pack_rows("E param by G", 5, 6, indent = F) %>% 
  #   footnote(general = paste0("Risk allele - ", risk_allele)) %>% 
  #   save_kable(file = paste0(output_dir, "stratified_oddsratio_", snp, "_", exposure, ".html"), self_contained = T)
  
  saveRDS(est, file = paste0(output_dir, "stratified_oddsratio_", snp, "_", exposure, ".rds"))
  # return(est)
  
}















