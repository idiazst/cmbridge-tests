# Read-only calculations on saved fits; no nuisance or estimator is fitted.
study_source<-Sys.getenv('STUDY_SOURCE','simulations/observed_history_corrected')
source(file.path(study_source,'study.R'))
out<-Sys.getenv('SIM_OUTPUT','results/observed-history-corrected/missing-both-n10000')
parts<-list()
for(name in c('binary_longitudinal','discrete_dose')) {
  x<-readRDS(file.path(out,'jobs',paste0(name,'.rds')))
  g<-make_mechanism(name,visits='missing_both',dose_max=3L)
  p<-x$population;truth<-population_truth_functions(p,g)
  for(j in seq_along(x$folds)) {
    # SDR/TMLE share these bridge/adjoint functions; use either representation.
    nu<-x$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',j,sep='/')]]
    measured<-p$R3==1L
    term<-numeric(nrow(p))
    term[measured]<--p$probability[measured]*(nu$lambda[measured,2L]-truth$lambda[measured,2L])*
      (nu$beta[measured,2L]-truth$beta[measured,2L])
    z<-data.table(R2=p$R2,term=term)[,.(bridge_remainder=sum(term)),by=R2]
    z[,`:=`(mechanism=name,fold=j,fold_weight=length(x$folds[[j]]$validation_set)/x$design$n)]
    parts[[length(parts)+1L]]<-z
  }
}
parts<-rbindlist(parts)
summary<-parts[,.(bridge_remainder=sum(fold_weight*bridge_remainder)),by=.(mechanism,R2)]
fwrite(parts,file.path(out,'bridge-remainder-by-measurement-fold.csv'))
fwrite(summary,file.path(out,'bridge-remainder-by-measurement.csv'))
w<-fread(file.path(out,'ensemble_weights.csv'))
ws<-w[kind %in% c('beta','adjoint') & horizon_or_depth==2L,
  .(mean_weight=mean(weight)),by=.(mechanism,kind,candidate)]
fwrite(ws,file.path(out,'final-outcome-cmbridge-weights.csv'))
print(summary);print(ws)
cat('Saved-fit diagnostics completed without fitting.\n')
