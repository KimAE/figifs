#=============================================================================#
# functions to perform some posthoc analyses
# stratified odds ratios table
#=============================================================================#

#' format_res
#' 
#' clean up output from Model.all.new
#'
#' @param res 
#'
#' @return
#' @export
format_res <- function(res) {
  betas = colnames(res)[grep('beta',colnames(res),fixed=T)]
  ses = colnames(res)[grep('se',colnames(res),fixed=T)]
  ORs = sapply(seq(length(betas)),
               function(x, betas, ses, res) {
                 or = paste0(round(exp(res[,betas[x]]), 2),' (',
                             round(exp(res[,betas[x]]-qnorm(0.975)*res[,ses[x]]), 2), '-',
                             round(exp(res[,betas[x]]+qnorm(0.975)*res[,ses[x]]), 2), ')')
               }, 
               betas=betas,ses=ses,res=res)
  
  if(nrow(res)==1)  ORs = data.frame(t(ORs)) 
  
  colnames(ORs) = sub('beta','OR',betas)
  pvals = sapply(seq(length(betas)),
                 function(x,betas,ses,res) {
                   pval = 2*pnorm(-abs(res[,betas[x]]/res[,ses[x]]))
                   },
                 betas=betas,ses=ses,res=res)
  
  if(nrow(res)==1) pvals = data.frame(t(pvals))
  
  colnames(pvals) = sub('beta','pval',betas)
  
  pvals.p = apply(pvals,2, function(y) paste0('P= ',formatC(y,format='g',digit=2)))
  
  if(nrow(res)==1) pvals.p = data.frame(t(pvals.p)) 
  colnames(pvals.p) = paste0('P',colnames(pvals.p))
  res = data.frame(res,ORs,pvals,pvals.p,stringsAsFactors=F)
}



#' Model.all.new
#'
#' Estimate stratified odds ratios + CIs. Written by Yi, with a few modifications to accomodate continuous E. 
#'
#' @param data dataset
#' @param mod string of model to be fitted, exclude dependent var e.g. "age_ref_imp+sex+study_gxe+pc1+pc2+pc3"
#' @param env string of exposure variable e.g. 'asp_ref'
#'
#' @return a list containing model estimates and standard errors for stratified table cells
#' @export
#'
#' @examples Model.all.new(figi, "age_ref_imp+sex+study_gxe+pc1+pc2+pc3", "asp_ref")
Model.all.new <- function(data, mod, env) {
  
  mModel <- paste0('outcome~', env, '*p1+', env, '*p2+', mod)
  tmp <- summary(glm(mModel, family='binomial', data=data))
  COV <- tmp$cov.unscaled
  
  # if factor, need to get correct glm model parameter
  env.idx = env
  if(is.factor(data[,env])) env.idx = paste0(env, levels(data[,env])[-1])
  
  # stratified by GE  
  eg0 <- tmp$coef[env.idx,c(1,2)]
  e0g <- tmp$coef[c('p1','p2'), c(1,2)]
  beta.eg1 <- tmp$coef[env.idx,1] + tmp$coef['p1',1] + tmp$coef[paste0(env.idx,':p1'),1]
  beta.eg2 <- tmp$coef[env.idx,1] + tmp$coef['p2',1] + tmp$coef[paste0(env.idx,':p2'),1]
  
  I3 =  rep(1,3)
  se.eg1 = se.eg2 = NULL
  if(is.factor(data[,env])) {
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs = COV[c(env.idx[i], 'p1', paste0(env.idx[i], ':p1')), 
                 c(env.idx[i], 'p1', paste0(env.idx[i], ':p1'))]		
      se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
      covs = COV[c(env.idx[i], 'p2', paste0(env.idx[i], ':p2')), 
                 c(env.idx[i], 'p2', paste0(env.idx[i], ':p2'))]		
      se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
    }
  } else {
    covs = COV[c(env.idx, 'p1', paste0(env.idx,':p1')), 
               c(env.idx, 'p1', paste0(env.idx,':p1'))]		
    se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
    covs = COV[c(env.idx, 'p2', paste0(env.idx,':p2')),
               c(env.idx, 'p2', paste0(env.idx,':p2'))]		
    se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
  }
  
  GE<- c(eg0,e0g,beta.eg1,se.eg1,beta.eg2,se.eg2)
  
  
  # stratified by G 
  # (remember there are variables that can be categorical with more than 2 levels)
  g0 <-tmp$coef[env.idx, 1:2]
  if(length(env.idx) == 1) { # two level factor or continuous E... 
    idx.g1 = c(env.idx, paste0(env.idx,':p1'))
    idx.g2 = c(env.idx, paste0(env.idx,':p2'))
    g1 = sum(tmp$coef[idx.g1, 1],na.rm=T)
    g2 = sum(tmp$coef[idx.g2, 1],na.rm=T)
    Is1 =  rep(1,length(idx.g1))
    Is2 =  rep(1,length(idx.g2))
  } else { # categorical with 2+ levels
    est.g1 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p1'),1])
    est.g2 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p2'),1])
    g1 = rowSums(est.g1,na.rm=T)
    g2 = rowSums(est.g2,na.rm=T)
    Is1 =  rep(1,ncol(est.g1))
    Is2 =  rep(1,ncol(est.g2))
  }
  
  se.g1 = se.g2 = NULL
  if(is.factor(data[,env])) {
    for(i in 1:(nlevels(data[,env])-1)) {
      covs1 <- COV[c(env.idx[i], paste0(env.idx[i], ':p1')),
                   c(env.idx[i], paste0(env.idx[i], ':p1'))]  
      covs2 <- COV[c(env.idx[i], paste0(env.idx[i], ':p2')),
                   c(env.idx[i], paste0(env.idx[i], ':p2'))]  
      se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
      se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
    }
  } else {
    covs1 <- COV[c(env.idx, paste0(env.idx, ':p1')),
                 c(env.idx, paste0(env.idx, ':p1'))]  
    covs2 <- COV[c(env.idx, paste0(env.idx, ':p2')),
                 c(env.idx, paste0(env.idx, ':p2'))]  
    se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
    se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
  }
  G <- c(g0,g1,se.g1,g2,se.g2)
  
  # stratified by E
  e0 <-tmp$coef[c('p1','p2'),1:2]
  if(length(env.idx) == 1) {
    idx.e11 = c('p1',paste0(env.idx,':p1'))
    idx.e12 = c('p2',paste0(env.idx,':p2'))
    e11 = sum(tmp$coef[idx.e11,1],na.rm=T)
    e12 = sum(tmp$coef[idx.e12,1],na.rm=T)
    Is1 =  rep(1,length(idx.e11))
    Is2 =  rep(1,length(idx.e12))
  } else {
    est.e11 <- data.frame(tmp$coef['p1',1],tmp$coef[paste0(env.idx,':p1'),1])
    est.e12 <- data.frame(tmp$coef['p2',1],tmp$coef[paste0(env.idx,':p2'),1])
    e11 = rowSums(est.e11,na.rm=T)
    e12 = rowSums(est.e12,na.rm=T)
    Is1 =  rep(1,ncol(est.e11))
    Is2 =  rep(1,ncol(est.e12))
  }
  
  se.e11 = se.e12 = NULL
  if(is.factor(data[,env])) {
    for(i in 1:(nlevels(data[,env])-1)) {
      covs1 = COV[c('p1', paste0(env.idx[i], ':p1')),
                  c('p1', paste0(env.idx[i], ':p1'))]  
      covs2 = COV[c('p2', paste0(env.idx[i], ':p2')),
                  c('p2', paste0(env.idx[i], ':p2'))]  
      se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
      se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
    }
  } else {
    covs1 = COV[c('p1', paste0(env.idx, ':p1')),
                c('p1', paste0(env.idx, ':p1'))]  
    covs2 = COV[c('p2', paste0(env.idx, ':p2')),
                c('p2', paste0(env.idx, ':p2'))]  
    se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
    se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
  }
  E <- c(e0[1,1],e11,e0[1,2],se.e11,e0[2,1],e12,e0[2,2],se.e12)
  
  res <- list(GE=GE,G=G,E=E)
  
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
ptab = function(x, elvl, glvl, res) {
  pval = data.frame(matrix(NA,
                           nrow=length(elvl),
                           ncol=length(glvl)))
  
  colnames(pval) = paste0('p',glvl)
  rownames(pval) = elvl
  
  pval[1,1] = 1
  
  for(c in colnames(pval))
    for(r in rownames(pval))
      if(!(c=='p0' & r==elvl[1])) 
        pval[r,c] =  as.numeric(as.vector(res[res$snp %in% x, paste0('pval', r, c)]))
  pval = apply(pval, 2, function(y) paste0('P= ', formatC(y, format='g', digit=2)))    
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
ORtab = function(x, elvl, glvl, res) {
  ORs = data.frame(matrix(NA, nrow=length(elvl), ncol=length(glvl)))
  rownames(ORs) = paste0('OR',elvl)
  colnames(ORs) = paste0('p',glvl)
  
  # reference group
  ORs[1,1] = 1
  
  # loop
  for(c in colnames(ORs))
    for(r in rownames(ORs))
      if(!(c=='p0' & r==paste0('OR', elvl[1]))) 
        ORs[r,c] = as.character(res[res$snp %in% x, paste0(r,c)])
  ORs      
} 




#' fit_stratified_or
#'
#' @param data_epi 
#' @param exposure 
#' @param snp 
#' @param hrc_version 
#' @param covariates 
#' @param dosage 
#' @param path 
#' @param flip_allele 
#'
#' @return
#' @export
#'
#' @examples
fit_stratified_or <- function(data_epi, exposure, snp, hrc_version, covariates, dosage = F, path, flip_allele = F) {
  
  mod = glue_collapse(covariates, sep = "+")
  wdir = glue("{path}/posthoc")
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = ":"))
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  
  # read dosage file
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  if (flip_allele == T) {
    # flip dosages
    # data[, snp] <- abs(data[, paste0(snp)] - 2)
    data[, snpfix] <- abs(data[, snpfix] - 2)
    # flip genotype probabilities
    pp <- data[,paste0(snpfix_short, "_p2")]
    data[,paste0(snpfix_short, "_p2")] <- data[, paste0(snpfix_short, "_p0")]
    data[,paste0(snpfix_short, "_p0")] <- pp
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }
  
  # create data subset
  tmp1 = data[, c('outcome', exposure, covariates)]
  tmp2 = data[, grepl(paste(paste0(snpfix_short, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(data))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
  
  ds_tmp = cbind(tmp1, tmp2) %>%
    na.omit(.[, c(exposure, 'outcome', covariates, 'p1', 'p2', 'dosage')]) # from data.table. remove rows with NA in these columns)
  

  
  #---- Calculate counts for each cell ----------------
  Ncaco = data.frame(snps=snpfix_short,matrix(NA,ncol=6*2))
  colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                      paste0('Co1',1:2),paste0('Ca1',1:2),
                      paste0('Co2',1:2),paste0('Ca2',1:2))
  rownames(Ncaco) = snpfix_short
  
  Ncaco[snpfix_short,c(paste0('Co',1:2),paste0('Ca',1:2))] = t(table(ds_tmp[,c('outcome',exposure)]))
  Ncaco[snpfix_short,c(paste0('Co1',1:2),paste0('Ca1',1:2))] = c(t(tapply(ds_tmp$p1, ds_tmp[,c('outcome',exposure)], sum, na.rm=T)))
  Ncaco[snpfix_short,c(paste0('Co2',1:2),paste0('Ca2',1:2))] = c(t(tapply(ds_tmp$p2, ds_tmp[,c('outcome',exposure)], sum, na.rm=T)))
  
  #== calculate counts for G=0 and put counts into format ca/co
  for(i in 1:2){
    eval(parse(text=paste0('Ncaco$caco0',i,"=paste0(round(Ncaco$Ca",i,'-Ncaco$Ca1',i,'-Ncaco$Ca2',i,
                           ",1),'/',round(Ncaco$Co",i,'-Ncaco$Co1',i,'-Ncaco$Co2',i,',1))')))
    eval(parse(text=paste0('Ncaco$caco1',i,'=paste0(round(Ncaco$Ca1',i,",1),'/',round(Ncaco$Co1",i,",1))")))
    eval(parse(text=paste0('Ncaco$caco2',i,'=paste0(round(Ncaco$Ca2',i,",1),'/',round(Ncaco$Co2",i,",1))")))
  }
  
  
  
  #-- Fit unrestricted model --------
  # ds_tmp[,exposure] = factor(ds_tmp[,exposure])
  tmp = as.vector(Model.all.new(ds_tmp, mod, exposure))
  
  res.pool.un = res.pool.g = res.pool.e = NULL
  res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(snpfix_short,exposure,t(tmp$G)))
  
  ## organize results ######
  ## need to automate elvl eventually.. 
  if(is.factor(ds_tmp[,exposure])) {
    elvl <- c(0:(nlevels(ds_tmp[, exposure])-1))
  } else {
    elvl <- c(0,1)
  }
  
  # always set to genotype (have different function for per allelic dosage.. )
  glvl <- c(0:2)
  
  ##== GE stratified results ########
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
  OR.tab = ORtab(snpfix_short,
                 elvl=elvl,
                 glvl=glvl,
                 res=res.pool.un)
  pval.tab = ptab(snpfix_short,elvl=elvl,glvl=glvl,res=res.pool.un)
  
  ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('OR0',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR1',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR2',elvl[-1])])))), ncol=3)
  pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('Ppval0',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval1',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval2',elvl[-1])])))),ncol=3)
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))
  
  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)],
                                         res.pool.e[,paste0('OR2',elvl)]))), ncol=2)
  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)],
                                         res.pool.e[,paste0('Ppval2',elvl)]))),ncol=2,byrow=T)
  


  #=== write into table
  if(is.factor(ds_tmp[, exposure])) {
    est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
                      ORg.tab[1,],pg.tab[1,]),
                rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)),
                N0 = c(Ncaco[snpfix_short,'caco01'],NA,Ncaco[snpfix_short,'caco02'],rep(NA,3)),
                N1 = c(Ncaco[snpfix_short,'caco11'],NA,Ncaco[snpfix_short,'caco12'],rep(NA,3)),
                N2 = c(Ncaco[snpfix_short,'caco21'],NA,Ncaco[snpfix_short,'caco22'],rep(NA,3)))
    colnames(est) <- c(paste0(ref_allele, ref_allele),
                       paste0(ref_allele, alt_allele),
                       paste0(alt_allele, alt_allele),
                       paste0(ref_allele, alt_allele),
                       paste0(alt_allele, alt_allele),
                       "N0", "N1", "N2")
    exposure_level <- levels(ds_tmp[,exposure])
    rownames(est) <- c(paste0(exposure,"=",exposure_level[1]),
                       "p0",
                       paste0(exposure,"=",exposure_level[2]),
                       "p1",
                       exposure,
                       "p")
    
  } else {
    est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
                      ORg.tab[1,],pg.tab[1,]),
                rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)))
    colnames(est) <- c(paste0(ref_allele, ref_allele),
                       paste0(ref_allele, alt_allele),
                       paste0(alt_allele, alt_allele),
                       paste0(ref_allele, alt_allele),
                       paste0(alt_allele, alt_allele))
    exposure_level <- levels(ds_tmp[,exposure])
    rownames(est) <- c(paste0(exposure,"=Ref"),
                       "p0",
                       paste0(exposure,"=UnitInc"),
                       "p1",
                       exposure,
                       "p")
  }
  
  saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))

  
  # # slightly modify output in case anyone wants per allele effects rather than p1/p2
  # if(dosage == T) {
  #   
  #   tmp = as.vector(Model.all.new.dosage(ds_tmp,mod,exposure))
  #   res.pool.un = res.pool.g = res.pool.e = NULL
  #   res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
  #   
  #   elvl=c(0:1) ; glvl = c(0,1)
  #   tmp2 = NULL
  #   tmp2 = rbind(tmp2,data.frame(snpfix_short,exposure,t(tmp$GE)))
  #   
  #   colnames(tmp2) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
  #                      'beta0p1','beta0p2','se0p1','se0p2')
  #   tmp2 = format_res(tmp2)
  #   
  #   res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
  #   colnames(res.pool.e) = c('snp',
  #                            'env',
  #                            paste0('beta1',elvl),paste0('se1',elvl))
  #   tmp2 <- format_res(res.pool.e)
  #   
  #   
  #   ORe.tab = matrix(as.character(unlist(c(tmp2[,paste0('OR1',elvl)]))),ncol=2)
  #   pe.tab  = matrix(as.character(unlist(c(tmp2[,paste0('Ppval1',elvl)]))),ncol=2,byrow=T)
  #   est2 <- rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2], NA, NA)
  #   
  #   final_out <- est %>%
  #     dplyr::select(-c(4,5)) %>%
  #     add_column(est2, .after = 3)
  #   
  #   colnames(final_out) <- c(paste0(ref_allele, ref_allele),
  #                            paste0(ref_allele, alt_allele),
  #                            paste0(alt_allele, alt_allele),
  #                            paste0(alt_allele, " Allelic Dosage"),
  #                            "N0", "N1", "N2")
  #   # this is only for FACTOR variables (might have to modify when you run Q4 variables)
  #   exposure_level <- levels(data[,exposure])
  #   rownames(final_out) <- c(paste0(exposure,"=",exposure_level[1]),
  #                            "p0",
  #                            paste0(exposure,"=",exposure_level[2]),
  #                            "p1",
  #                            exposure,
  #                            "p")
  #   saveRDS(final_out, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}_dosage.rds"))
  #   return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}_dosage.rds"))
  # } else {
  #   colnames(est) <- c(paste0(ref_allele, ref_allele),
  #                      paste0(ref_allele, alt_allele),
  #                      paste0(alt_allele, alt_allele),
  #                      paste0(ref_allele, alt_allele),
  #                      paste0(alt_allele, alt_allele),
  #                      "N0", "N1", "N2")
  #   exposure_level <- levels(data[,exposure])
  #   rownames(est) <- c(paste0(exposure,"=",exposure_level[1]),
  #                      "p0",
  #                      paste0(exposure,"=",exposure_level[2]),
  #                      "p1",
  #                      exposure,
  #                      "p")
  #   
  #   saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  #   return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  #   # saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}.rds"))
  # }
  
}

# ---------------------------------------------------------- #
# create stratified odds ratios for a continuous variable per unit increase
# by genotype
# this is the bit of code i had written before and  couldn't find so i ended up writing it again anyway
# ---------------------------------------------------------- #


##' Model.all.new.continuous
##' 
##' for red/proc meat analysis, Dr. stern wanted stratified odds ratios, and i think it's worth doing it for a unit increase in exposure
##' for each genotype (informative). I think this is the correct way of fitting the model. 
##'
##' @param data 
##' @param mod 
##' @param env 
##'
##' @return
##' @export
##'
##' @examples
#Model.all.new.continuous <- function(data,mod,env){
#  
#  # error check
#  if(length(table(data$p1))<=1 | length(table(data$p2))<=1 | 
#     var(2*data$p2+data$p1,na.rm=T)==0 | 
#     sum((data$p2+data$p1)>1.1,na.rm=T)>0) {
#    res = NA
#  } else {
#    mModel <- paste0('outcome~',env,'*p1+',env,'*p2+',mod)
#    tmp <- summary(glm(mModel, family='binomial', data=data))
#    COV <- tmp$cov.unscaled
#    env.idx = env
#    
#    # stratified by GE  
#    eg0 <- tmp$coef[env.idx,c(1,2)] # E
#    beta.eg1 <- tmp$coef[env.idx,1]+tmp$coef['p1',1]+tmp$coef[paste0(env.idx,':p1'),1] # E+p1+E:p1
#    beta.eg2 <- tmp$coef[env.idx,1]+tmp$coef['p2',1]+tmp$coef[paste0(env.idx,':p2'),1] # E+p2+E:p2
#    e0g <- tmp$coef[c('p1','p2'),c(1,2)] # p1, p2
#    I3 =  rep(1,3)
#    
#    se.eg1 = se.eg2 = NULL
#    covs = COV[c(env.idx,'p1',paste0(env.idx,':p1')),c(env.idx,'p1',paste0(env.idx,':p1'))]		
#    se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3)) # calculating SE of the combined E+G+E:G estimate 
#    covs = COV[c(env.idx,'p2',paste0(env.idx,':p2')),c(env.idx,'p2',paste0(env.idx,':p2'))]		
#    se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3)) # calculating SE of the combined E+G+E:G estimate 
#    GE<- c(eg0,e0g,beta.eg1,se.eg1,beta.eg2,se.eg2) # E, G, E:p1, se.E:p1, E:p2, se.E:p2
#    
#    
#    
#    
#    # stratified by G 
#    g0 <-tmp$coef[env.idx,1:2]
#    if(length(env.idx)==1){ # confused about this check.. 
#      idx.g1 = c(env.idx,paste0(env.idx,':p1'))
#      idx.g2 = c(env.idx,paste0(env.idx,':p2'))
#      g1 = sum(tmp$coef[idx.g1,1],na.rm=T)
#      g2 = sum(tmp$coef[idx.g2,1],na.rm=T)
#      Is1 =  rep(1,length(idx.g1))
#      Is2 =  rep(1,length(idx.g2))
#    }else{
#      est.g1 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p1'),1])
#      est.g2 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p2'),1])
#      g1 = rowSums(est.g1,na.rm=T)
#      g2 = rowSums(est.g2,na.rm=T)
#      Is1 =  rep(1,ncol(est.g1))
#      Is2 =  rep(1,ncol(est.g2))
#    }
#    
#
#    se.g1 = se.g2 = NULL
#    covs1 = COV[c(env.idx,paste0(env.idx,':p1')),c(env.idx,paste0(env.idx,':p1'))]  
#    covs2 = COV[c(env.idx,paste0(env.idx,':p2')),c(env.idx,paste0(env.idx,':p2'))]  
#    se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
#    se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
#    G<- c(g0,g1,se.g1,g2,se.g2) # E, E+E:p1, se.E+E:p1, E+E:p2, se.E+E:p2
#    
#    
#    # stratified by E
#    
#    e0 <-tmp$coef[c('p1','p2'),1:2]
#    if(length(env.idx)==1){
#      idx.e11 = c('p1',paste0(env.idx,':p1'))
#      idx.e12 = c('p2',paste0(env.idx,':p2'))
#      e11 = sum(tmp$coef[idx.e11,1],na.rm=T) # p1 + E:p1
#      e12 = sum(tmp$coef[idx.e12,1],na.rm=T)
#      Is1 =  rep(1,length(idx.e11))
#      Is2 =  rep(1,length(idx.e12))
#    }else{
#      est.e11 <- data.frame(tmp$coef['p1',1],tmp$coef[paste0(env.idx,':p1'),1])
#      est.e12 <- data.frame(tmp$coef['p2',1],tmp$coef[paste0(env.idx,':p2'),1])
#      e11 = rowSums(est.e11,na.rm=T)
#      e12 = rowSums(est.e12,na.rm=T)
#      Is1 =  rep(1,ncol(est.e11))
#      Is2 =  rep(1,ncol(est.e12))
#    }
#    
#    
#    se.e11 = se.e12 = NULL
#      covs1 = COV[c('p1',paste0(env.idx,':p1')),c('p1',paste0(env.idx,':p1'))]  
#      covs2 = COV[c('p2',paste0(env.idx,':p2')),c('p2',paste0(env.idx,':p2'))]  
#      se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
#      se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
#    E <- c(e0[1,1],e11,e0[1,2],se.e11,e0[2,1],e12,e0[2,2],se.e12) # p1, p1+E:p1, se.p1, se.p1+E:p1, p2, p2+E:p2, se.p2, se.p2+E:p2
#    
#    
#  }
#  res<-list(GE=GE,G=G,E=E)
#}
#
#
#
##' fit_stratified_or_continuous
##'
##' @param data_epi
##' @param exposure
##' @param snp
##' @param hrc_version
##' @param covariates
##' @param dosage
##' @param path
##' @param flip_allele
##'
##' @return
##' @export
##'
##' @examples
#fit_stratified_or_continuous <- function(data_epi, exposure, snp, hrc_version, covariates, dosage = F, path, flip_allele = F) {
#  
#  mod = glue_collapse(covariates, sep = "+")
#  
#  # wdir = glue("{path}/output/posthoc")
#  wdir = glue("{path}/posthoc")
#  
#  
#  # SNP information
#  snp_info <- unlist(strsplit(snp, split = ":"))
#  
#  snpname_clean <- function(x) {
#    tmp <- gsub("\\:", "\\_", x)
#    # tmp <- gsub("X", "chr", tmp)
#    tmp <- glue("chr{tmp}_dose")
#    return(tmp)
#  }
#  
#  snpfix <- snpname_clean(snp)
#  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
#  
#  #data_dose <- qread(glue("{wdir}/dosage_chr{chr}_{bp}.qs"))
#  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
#  
#  data <- inner_join(data_epi, data_dose, 'vcfid')
#
#  
#  if (flip_allele == T) {
#    snp_old <- snpfix_short
#    snp_tmp <- unlist(strsplit(snp, split = ":"))
#    chr <- snp_tmp[1]; bp <- snp_tmp[2]; a1 <- snp_tmp[3]; a2 <- snp_tmp[4]
#    snp_new <- glue("chr{snp_info[1]}_{snp_info[2]}_{snp_info[4]}_{snp_info[3]}")
#    
#    data[, paste0(snp_new, "_dose")] <- abs(data[, paste0(snp_old, "_dose")] - 2)
#    data[, paste0(snp_new, "_p2")] <- data[, paste0(snp_old, "_p0")]
#    data[, paste0(snp_new, "_p1")] <- data[, paste0(snp_old, "_p1")]
#    data[, paste0(snp_new, "_p0")] <- data[, paste0(snp_old, "_p2")]
#
#    ref_allele <- snp_info[4]
#    alt_allele <- snp_info[3]
#    
#  } else {
#    snp_old <- snpfix_short
#    snp_tmp <- unlist(strsplit(snp, split = ":"))
#    chr <- snp_tmp[1]; bp <- snp_tmp[2]; a1 <- snp_tmp[3]; a2 <- snp_tmp[4]
#    snp_new <- glue("chr{snp_info[1]}_{snp_info[2]}_{snp_info[3]}_{snp_info[4]}")
#    
#    ref_allele <- snp_info[3]
#    alt_allele <- snp_info[4]
#    
#  }
#  
#  
#  
#  
#  # create data subset
#  tmp1 = data[, c('outcome', exposure, covariates)]
#  # tmp2 = data[, grepl(paste(paste0(snp_new, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(data))]
#  tmp2 = data[, paste(paste0(snp_new, c("_dose", "_p0", "_p1", "_p2")))]
#  names(tmp2) <- c("dosage", "p0", "p1", "p2")
#  
#  ds_tmp = cbind(tmp1, tmp2) %>%
#    na.omit(.[, c(exposure,'outcome', covariates, 'p1','p2','dosage')])
#  
#  
#  # ---------------------------------- #
#
#  res.pool.un = res.pool.g = res.pool.e = NULL
#  
#  #-- Fit unrestricted model --------
#  tmp = as.vector(Model.all.new.continuous(ds_tmp,mod,exposure))
#  res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
#  res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
#  res.pool.g = rbind(res.pool.g,data.frame(snpfix_short,exposure,t(tmp$G)))
#  
#  ## organize results ######
#  elvl=c(0:1) ; glvl = c(0,1,2)
#  
#  colnames(res.pool.un) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
#                            'beta0p1','beta0p2','se0p1','se0p2',
#                            paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'),
#                            paste0('beta',elvl[-1],'p2'),paste0('se',elvl[-1],'p2'))
#  res.pool.un = format_res(res.pool.un)
#  
#  
#  ##== stratified by G results #######
#  colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1','beta2','se2'),each=1),
#                                              rep(elvl[-1],6)))
#  res.pool.g = format_res(res.pool.g)
#  
#  
#  
#  ##== stratified by E results #######
#  colnames(res.pool.e) = c('snp','env',paste0('beta1',elvl),paste0('se1',elvl),
#                           paste0('beta2',elvl),paste0('se2',elvl))
#  res.pool.e = format_res(res.pool.e)
#  
#  
#  ##== Put into table
#  OR.tab = ORtab(snpfix_short,elvl=elvl,glvl=glvl,res=res.pool.un)
#  pval.tab = ptab(snpfix_short,elvl=elvl,glvl=glvl,res=res.pool.un)
#  
#  ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('OR0',elvl[-1])]),
#                                         as.character(res.pool.g[,paste0('OR1',elvl[-1])]),
#                                         as.character(res.pool.g[,paste0('OR2',elvl[-1])])))), ncol=3)
#  pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('Ppval0',elvl[-1])]),
#                                        as.character(res.pool.g[,paste0('Ppval1',elvl[-1])]),
#                                        as.character(res.pool.g[,paste0('Ppval2',elvl[-1])])))),ncol=3)
#  
#  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))
#  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)],
#                                         res.pool.e[,paste0('OR2',elvl)]))),ncol=2)
#  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)],
#                                         res.pool.e[,paste0('Ppval2',elvl)]))),ncol=2,byrow=T)
#  
#  
#  #=== write into table
#  # est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
#  #                   ORg.tab[1,],pg.tab[1,]),
#  #             rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)),
#  #             N0 = c(Ncaco[snpfix_short,'caco01'],NA,Ncaco[snpfix_short,'caco02'],rep(NA,3)),
#  #             N1 = c(Ncaco[snpfix_short,'caco11'],NA,Ncaco[snpfix_short,'caco12'],rep(NA,3)),
#  #             N2 = c(Ncaco[snpfix_short,'caco21'],NA,Ncaco[snpfix_short,'caco22'],rep(NA,3)))
#  
#  est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
#                    ORg.tab[1,],pg.tab[1,]),
#              rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)))
#  
#  
#  
#  
#  # slightly modify output in case anyone wants per allele effects rather than p1/p2
#  if(dosage == T) {
#    
#    tmp = as.vector(Model.all.new.dosage(ds_tmp,mod,exposure))
#    res.pool.un = res.pool.g = res.pool.e = NULL
#    res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
#    
#    elvl=c(0:1) ; glvl = c(0,1)
#    tmp2 = NULL
#    tmp2 = rbind(tmp2,data.frame(snpfix_short,exposure,t(tmp$GE)))
#    
#    colnames(tmp2) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
#                       'beta0p1','beta0p2','se0p1','se0p2')
#    tmp2 = format_res(tmp2)
#    
#    res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
#    colnames(res.pool.e) = c('snp',
#                             'env',
#                             paste0('beta1',elvl),paste0('se1',elvl))
#    tmp2 <- format_res(res.pool.e)
#    
#    
#    ORe.tab = matrix(as.character(unlist(c(tmp2[,paste0('OR1',elvl)]))),ncol=2)
#    pe.tab  = matrix(as.character(unlist(c(tmp2[,paste0('Ppval1',elvl)]))),ncol=2,byrow=T)
#    est2 <- rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2], NA, NA)
#    
#    final_out <- est %>%
#      dplyr::select(-c(4,5)) %>%
#      add_column(est2, .after = 3)
#    
#    colnames(final_out) <- c(paste0(ref_allele, ref_allele),
#                             paste0(ref_allele, alt_allele),
#                             paste0(alt_allele, alt_allele),
#                             paste0(alt_allele, " Allelic Dosage"))
#    # this is only for FACTOR variables (might have to modify when you run Q4 variables)
#    # exposure_level <- levels(data[,exposure])
#    exposure_level <- c(0,1)
#    rownames(final_out) <- c(paste0(exposure,"=",exposure_level[1]),
#                             "p0",
#                             paste0(exposure,"=",exposure_level[2]),
#                             "p1",
#                             exposure,
#                             "p")
#    saveRDS(final_out, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}_dosage.rds"))
#    return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}_dosage.rds"))
#  } else {
#    colnames(est) <- c(paste0(ref_allele, ref_allele),
#                       paste0(ref_allele, alt_allele),
#                       paste0(alt_allele, alt_allele),
#                       paste0(ref_allele, alt_allele),
#                       paste0(alt_allele, alt_allele))
#    # exposure_level <- levels(data[,exposure])
#    exposure_level <- c(0,1)
#    rownames(est) <- c(paste0(exposure,"=",exposure_level[1]),
#                       "p0",
#                       paste0(exposure,"=",exposure_level[2]),
#                       "p1",
#                       exposure,
#                       "p")
#    
#    #return(est)
#    saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
#    return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
#    # saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}.rds"))
#  }
#  
#}





#' fit_stratified_or_q3
#'
#' @param data_epi
#' @param exposure
#' @param snp
#' @param hrc_version
#' @param covariates
#' @param dosage
#' @param path
#' @param flip_allele
#'
#' @return
#' @export
#'
#' @examples
fit_stratified_or_q3 <- function(data_epi, exposure, snp, hrc_version, covariates, dosage = F, path, flip_allele = F) {

  mod = glue_collapse(covariates, sep = "+")
  wdir = glue("{path}/posthoc")

  # SNP information
  snp_info <- unlist(strsplit(snp, split = ":"))
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }

  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))

  # read dosage file
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')

  if (flip_allele == T) {
    # flip dosages
    # data[, snp] <- abs(data[, paste0(snp)] - 2)
    data[, snpfix] <- abs(data[, snpfix] - 2)
    # flip genotype probabilities
    pp <- data[,paste0(snpfix_short, "_p2")]
    data[,paste0(snpfix_short, "_p2")] <- data[, paste0(snpfix_short, "_p0")]
    data[,paste0(snpfix_short, "_p0")] <- pp
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }

  # create data subset
  tmp1 = data[, c('outcome', exposure, covariates)]
  tmp2 = data[, grepl(paste(paste0(snpfix_short, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(data))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")

  ds_tmp = cbind(tmp1, tmp2) %>%
    na.omit(.[, c(exposure, 'outcome', covariates, 'p1', 'p2', 'dosage')]) # from data.table. remove rows with NA in these columns)



  #---- Calculate counts for each cell ----------------
  n_levels = 3
  Ncaco = data.frame(snps=snpfix_short,matrix(NA,ncol=6*n_levels))
  colnames(Ncaco) = c('snps',paste0('Co',1:n_levels),paste0('Ca',1:n_levels),
                      paste0('Co1',1:n_levels),paste0('Ca1',1:n_levels),
                      paste0('Co2',1:n_levels),paste0('Ca2',1:n_levels))
  rownames(Ncaco) = snpfix_short

  Ncaco[snpfix_short,c(paste0('Co',1:n_levels),paste0('Ca',1:n_levels))] = t(table(ds_tmp[,c('outcome',exposure)]))
  Ncaco[snpfix_short,c(paste0('Co1',1:n_levels),paste0('Ca1',1:n_levels))] = c(t(tapply(ds_tmp$p1, ds_tmp[,c('outcome',exposure)], sum, na.rm=T)))
  Ncaco[snpfix_short,c(paste0('Co2',1:n_levels),paste0('Ca2',1:n_levels))] = c(t(tapply(ds_tmp$p2, ds_tmp[,c('outcome',exposure)], sum, na.rm=T)))

  #== calculate counts for G=0 and put counts into format ca/co
  for(i in 1:n_levels){
    eval(parse(text=paste0('Ncaco$caco0',i,"=paste0(round(Ncaco$Ca",i,'-Ncaco$Ca1',i,'-Ncaco$Ca2',i,
                           ",1),'/',round(Ncaco$Co",i,'-Ncaco$Co1',i,'-Ncaco$Co2',i,',1))')))
    eval(parse(text=paste0('Ncaco$caco1',i,'=paste0(round(Ncaco$Ca1',i,",1),'/',round(Ncaco$Co1",i,",1))")))
    eval(parse(text=paste0('Ncaco$caco2',i,'=paste0(round(Ncaco$Ca2',i,",1),'/',round(Ncaco$Co2",i,",1))")))
  }



  #-- Fit unrestricted model --------
  # ds_tmp[,exposure] = factor(ds_tmp[,exposure])
  tmp = as.vector(Model.all.new(ds_tmp, mod, exposure))

  res.pool.un = res.pool.g = res.pool.e = NULL
  res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(snpfix_short,exposure,t(tmp$G)))

  ## organize results ######
  ## need to automate elvl eventually..
  if(is.factor(ds_tmp[,exposure])) {
    elvl <- c(0:(nlevels(ds_tmp[, exposure])-1))
  } else {
    elvl <- c(0,1)
  }

  # always set to genotype (have different function for per allelic dosage.. )
  glvl <- c(0:2)

  ##== GE stratified results ########
  colnames(res.pool.un) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                            'beta0p1','beta0p2','se0p1','se0p2',
                            paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'),
                            paste0('beta',elvl[-1],'p2'),paste0('se',elvl[-1],'p2'))
  res.pool.un = format_res(res.pool.un)

  ##== stratified by G results #######
  colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1','beta2','se2'),each=2),
                                              rep(elvl[-1],6)))
  res.pool.g = format_res(res.pool.g)

  ##== stratified by E results #######
  colnames(res.pool.e) = c('snp','env',paste0('beta1',elvl),paste0('se1',elvl),
                           paste0('beta2',elvl),paste0('se2',elvl))
  res.pool.e = format_res(res.pool.e)


  ##== Put into table
  OR.tab = ORtab(snpfix_short,
                 elvl=elvl,
                 glvl=glvl,
                 res=res.pool.un)
  pval.tab = ptab(snpfix_short,elvl=elvl,glvl=glvl,res=res.pool.un)

  ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('OR0',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR1',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR2',elvl[-1])])))), ncol=3)
  pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('Ppval0',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval1',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval2',elvl[-1])])))),ncol=3)
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))

  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)],
                                         res.pool.e[,paste0('OR2',elvl)]))), ncol=2)
  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)],
                                         res.pool.e[,paste0('Ppval2',elvl)]))),ncol=2,byrow=T)




  est = cbind(rbind(OR.tab[1,],pval.tab[1,],
                    OR.tab[2,],pval.tab[2,],
                    OR.tab[3,],pval.tab[3,],
                    ORg.tab[1,],pg.tab[1,],
                    ORg.tab[2,],pg.tab[2,]),
              rbind(ORe.tab[1,],pe.tab[1,],
                    ORe.tab[2,],pe.tab[2,],
                    ORe.tab[3,],pe.tab[3,],
                    matrix(data = NA, nrow = (n_levels-1)*2, ncol = 2)),
              N0 = c(Ncaco[snpfix_short,'caco01'],NA,Ncaco[snpfix_short,'caco02'],NA,Ncaco[snpfix_short,'caco03'],rep(NA,5)),
              N1 = c(Ncaco[snpfix_short,'caco11'],NA,Ncaco[snpfix_short,'caco12'],NA,Ncaco[snpfix_short,'caco13'],rep(NA,5)),
              N2 = c(Ncaco[snpfix_short,'caco21'],NA,Ncaco[snpfix_short,'caco22'],NA,Ncaco[snpfix_short,'caco23'],rep(NA,5)))

  colnames(est) <- c(paste0(ref_allele, ref_allele),
                     paste0(ref_allele, alt_allele),
                     paste0(alt_allele, alt_allele),
                     paste0(ref_allele, alt_allele),
                     paste0(alt_allele, alt_allele),
                     "N0", "N1", "N2")
  exposure_level <- levels(as.factor(ds_tmp[,exposure]))
  rownames(est) <- c(paste0(exposure,"=",exposure_level[1]),
                     "p0",
                     paste0(exposure,"=",exposure_level[2]),
                     "p1",
                     paste0(exposure,"=",exposure_level[3]),
                     "p2",
                     paste0(exposure,"=",exposure_level[2], "(G)"),
                     "p1 (G)",
                     paste0(exposure,"=",exposure_level[3], "(G)"),
                     "p2 (G)")

}

