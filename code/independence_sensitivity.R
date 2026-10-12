options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(TwoSampleMR);library(MVMR);library(MendelianRandomization)})
args<-commandArgs(trailingOnly=TRUE)
src<-if(length(args)>0) args[1] else 'bmc_reproducibility'
out<-if(length(args)>1) args[2] else 'audit'
flags<-read.csv(file.path(out,'instrument_flags.csv'),fileEncoding='UTF-8-BOM')
stopifnot('SNP' %in% names(flags),sum(as.logical(flags$exclude_primary))>0)
keys<-c('cad','pad_finngen_r5','ais','laa','svs','ces')
res<-list();removed<-list()
for(scheme in c('all','primary','expanded')) {
 exc<-if(scheme=='all') character() else flags$SNP[as.logical(flags[[paste0('exclude_',scheme)]])]
 for(k in c(keys,'x_to_m','reverse_bmi_to_lifetime')) {
  d0<-read.csv(file.path(src,'data_archive',paste0(k,'_harmonised.csv')))
  if('mr_keep' %in% names(d0)) d0<-d0[d0$mr_keep,]
  d<-d0[!d0$SNP %in% exc,]
  fit<-mr_ivw_mre(d$beta.exposure,d$beta.outcome,d$se.exposure,d$se.outcome)
  res[[paste(scheme,k,'ivw')]]<-data.frame(scheme=scheme,key=k,model='univariable',exposure=if(k=='reverse_bmi_to_lifetime')'bmi' else 'smoking',n=nrow(d),removed=nrow(d0)-nrow(d),b=fit$b,se=fit$se,p=fit$pval,conditional_F=NA)
  removed[[paste(scheme,k,'ivw')]]<-data.frame(scheme=scheme,key=k,model='univariable',removed_snps=paste(setdiff(d0$SNP,d$SNP),collapse=';'))
  if(!k %in% keys)next
  x0<-readRDS(file.path(src,'data_archive',paste0(k,'_mvmr_input.rds')))
  x<-x0[!x0$SNP %in% exc,];class(x)<-class(x0)
  cf<-NULL;invisible(capture.output(cf<-strength_mvmr(x,gencov=0)))
  # Exact original MVMR implementation: weighted no-intercept regression; P uses residual-df t test.
  invisible(capture.output(f<-ivw_mvmr(x,gencov=0)))
  for(j in 1:2)res[[paste(scheme,k,j)]]<-data.frame(scheme=scheme,key=k,model='multivariable',exposure=c('smoking','bmi')[j],n=nrow(x),removed=nrow(x0)-nrow(x),b=f[j,1],se=f[j,2],p=f[j,4],conditional_F=cf[1,j])
  removed[[paste(scheme,k,'mv')]]<-data.frame(scheme=scheme,key=k,model='multivariable',removed_snps=paste(setdiff(x0$SNP,x$SNP),collapse=';'))
 }
}
r<-do.call(rbind,res);r$fdr<-NA_real_
for(s in unique(r$scheme))for(m in unique(r$model))for(e in unique(r$exposure)) {
 ix<-which(r$scheme==s & r$model==m & r$exposure==e & r$key %in% keys)
 if(length(ix)>0)r$fdr[ix]<-p.adjust(r$p[ix],method='BH')
}
r$lo<-r$b-qnorm(.975)*r$se;r$hi<-r$b+qnorm(.975)*r$se
r$OR<-ifelse(r$key %in% keys,exp(r$b),NA);r$OR_lo<-ifelse(r$key %in% keys,exp(r$lo),NA);r$OR_hi<-ifelse(r$key %in% keys,exp(r$hi),NA)
core<-read.csv(file.path(src,'results/core_results_refitted.csv'))
for(k in keys)for(e in c('smoking','bmi')) {
 a<-r[r$scheme=='all' & r$key==k & r$model=='multivariable' & r$exposure==e,]
  b<-core[core$key==k,];stopifnot(abs(a$b-b[[paste0(e,'_b')]])<1e-9,abs(a$se-b[[paste0(e,'_se')]])<1e-9,abs(a$p-b[[paste0(e,'_p')]])<1e-9)
}
write.csv(r,file.path(out,'independence_sensitivity_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,removed),file.path(out,'excluded_snps_by_model.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo_independence.txt'))
print(r[r$scheme!='all' & r$model=='multivariable',c('scheme','key','exposure','n','OR','fdr','conditional_F')],row.names=FALSE)
