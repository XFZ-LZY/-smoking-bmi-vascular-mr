options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(TwoSampleMR);library(MVMR);library(MendelianRandomization)})
argv <- commandArgs(trailingOnly=TRUE)
src <- if(length(argv)>=1) argv[1] else 'output/bmc_reproducibility'
dest <- if(length(argv)>=2) argv[2] else 'output/strict_native_results'
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
keys <- c('cad','pad_finngen_r5','ais','laa','svs','ces')
core <- read.csv(file.path(src,'results/core_results_refitted.csv'))
checks <- list(); covsens <- list(); raw_outliers <- list()
for(k in keys) {
  d <- read.csv(file.path(src,'data_archive',paste0(k,'_harmonised.csv')))
  fit <- mr_ivw_mre(d$beta.exposure,d$beta.outcome,d$se.exposure,d$se.outcome)
  r <- core[core$key==k,]
  stopifnot(abs(fit$b-r$total_b)<1e-10,abs(fit$se-r$total_se)<1e-10)
  x <- readRDS(file.path(src,'data_archive',paste0(k,'_mvmr_input.rds')))
  se <- as.matrix(x[,c('sebetaX1','sebetaX2')])
  for(rho in seq(-0.9,0.9,by=0.3)) {
    rho <- round(rho,1)
    # Explicit hypothetical correlation of sampling errors, not measured overlap.
    gl <- lapply(seq_len(nrow(x)),function(i) {
      matrix(c(se[i,1]^2,rho*prod(se[i,]),rho*prod(se[i,]),se[i,2]^2),2,2)
    })
    invisible(capture.output(cf <- suppressWarnings(strength_mvmr(x,gencov=gl))))
    covsens[[paste(k,rho)]] <- data.frame(key=k,error_correlation=rho,F_smoking=cf[1,1],F_bmi=cf[1,2])
  }
  p <- readRDS(file.path(src,'data_archive/diagnostic_fits',paste0(k,'_mr_presso_3000.rds')))
  oo <- as.data.frame(p[['MR-PRESSO results']][['Outlier Test']])
  if(nrow(oo)>0) {
    oo$key <- k; oo$row_in_original_input <- rownames(oo)
    oo$SNP <- d$SNP[as.integer(rownames(oo))]
    raw_outliers[[k]] <- oo
  }
  checks[[k]] <- data.frame(key=k,abs_beta_difference=abs(fit$b-r$total_b),abs_SE_difference=abs(fit$se-r$total_se))
}
write.csv(do.call(rbind,covsens),file.path(dest,'S17_covariance_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(dest,'native_core_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,raw_outliers),file.path(dest,'PRESSO_SNP_outlier_tests.csv'),row.names=FALSE)

sd <- list()
for(pair in list(c('smoking_to_bmi','x_to_m_harmonised.csv'),c('bmi_to_smoking','reverse_bmi_to_lifetime_harmonised.csv'))) {
  d <- read.csv(file.path(src,'data_archive',pair[2]))
  rx <- get_r_from_bsen(d$beta.exposure,d$se.exposure,d$samplesize.exposure)
  ry <- get_r_from_bsen(d$beta.outcome,d$se.outcome,d$samplesize.outcome)
  tt <- mr_steiger(d$pval.exposure,d$pval.outcome,d$samplesize.exposure,d$samplesize.outcome,rx,ry)
  sd[[pair[1]]] <- data.frame(direction=pair[1],R2_exposure=tt$r2_exp,R2_outcome=tt$r2_out,correct_direction=tt$correct_causal_direction,p=tt$steiger_test,sensitivity_ratio=tt$sensitivity_ratio)
}
write.csv(do.call(rbind,sd),file.path(dest,'S16_native_Steiger.csv'),row.names=FALSE)
# Check normal-distribution reconstructed sensitivity CIs against retained S4 fits.
sc <- read.csv(file.path(src,'results/table_s4_mvmr_bmi_sensitivity.csv'))
sc$lci <- exp(sc$b-qnorm(.975)*sc$se)
sc$uci <- exp(sc$b+qnorm(.975)*sc$se)
write.csv(sc,file.path(dest,'S4_exact_normal_CI.csv'),row.names=FALSE)
sc_checks <- list()
for(k in c('laa','svs','ces')) {
  x <- readRDS(file.path(src,'data_archive/diagnostic_fits',paste0(k,'_robust_mvmr_fits.rds')))
  for(nm in c('ivw_random_robust','mv_egger','mv_median')) {
    m <- c(ivw_random_robust='Robust IVW',mv_egger='MVMR-Egger',mv_median='Multivariable median')[[nm]]
    f <- x[[nm]]
    if(nm=='mv_egger'){b<-f@Estimate[2];s<-f@StdError.Est[2];lo<-f@CILower.Est[2];hi<-f@CIUpper.Est[2]}
    else {b<-f@Estimate[2];s<-f@StdError[2];lo<-f@CILower[2];hi<-f@CIUpper[2]}
    r <- sc[sc$key==k & sc$method==m,]
    sc_checks[[paste(k,nm)]] <- data.frame(key=k,method=m,SE_original=s,SE_reconstructed=r$se,abs_SE_difference=abs(s-r$se),CI_low_original=exp(lo),CI_high_original=exp(hi),abs_low_difference=abs(exp(lo)-r$lci),abs_high_difference=abs(exp(hi)-r$uci))
  }
}
ss<-do.call(rbind,sc_checks);write.csv(ss,file.path(dest,'S4_original_fit_CI_checks.csv'),row.names=FALSE)
stopifnot(all(ss$abs_SE_difference<1e-9),all(ss$abs_low_difference<1e-9),all(ss$abs_high_difference<1e-9))
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo_native.txt'))
cat('Native verification passed; conditional F across hypothetical sampling-error correlations:\n')
cf<-do.call(rbind,covsens);print(c(min_F_smoking=min(cf$F_smoking),max_F_smoking=max(cf$F_smoking),min_F_BMI=min(cf$F_bmi),max_F_BMI=max(cf$F_bmi)))
print(do.call(rbind,sd))
