# Component validation of the unchanged numerical-dose mechanism. Exact
# functions are evaluation targets only, never learner inputs.
study_source <- Sys.getenv('STUDY_SOURCE')
source(file.path(study_source, 'study.R'))
library(ggplot2)
out <- Sys.getenv('DEBUG_OUTPUT', 'results/observed-history-corrected/dose-cmbridge-debug')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, 'figures'), showWarnings = FALSE)
n_large <- as.integer(Sys.getenv('DEBUG_N', '1000000'))
seed_large <- 5204001L
g <- make_mechanism('discrete_dose', dose_max = 3L)
require_valid_mechanism(g)
p <- enumerate_data(g)
pt <- study_task(p, g)$task
truth <- population_truth_functions(p, g)
theta <- exact_truth(g)

design_at <- function(time) {
 H <- pt$vars$history('A', time); A <- pt$vars$A[[time]]
 list(B = as.matrix(pt$natural[, c(H, A), drop = FALSE]),
  V = as.matrix(cbind(pt$natural[, H, drop = FALSE],
   p[, c(paste0('C', time+1L, '_covariate'), paste0('Y', time+1L)), drop = FALSE])))
}
diagnose <- function(beta, lambda, time, n, seed, fold, source) {
 z <- design_at(time); m <- p[[paste0('R', time+1L)]] == 1
 w <- p$probability; y <- ifelse(m, p[[paste0('Y', time+1L)]], 0)
 omega <- apply(truth$ratios[, seq_len(time), drop = FALSE], 1, prod)
 beta_error <- sum(w * omega * y * beta) - theta[time]
 residual <- rep(-1, nrow(p)); residual[m] <- beta[m] - 1
 correction <- -sum(w * lambda * residual)
 paired <- beta_error + correction
 product <- -sum(w[m] * (lambda[m] - truth$lambda[m, time]) *
                         (beta[m] - truth$beta[m, time]))
 stopifnot(abs(paired - product) < 1e-8)
 bkey <- cell_key(z$B); vkey <- cell_key(z$V[m, , drop = FALSE])
 phi <- conditional_mean(omega[m]*y[m], w[m], vkey)
 adj <- conditional_mean(lambda[m], w[m], vkey)
 br <- conditional_mean(residual, w, bkey)
 bv <- data.table(cell=vkey, true=truth$beta[m,time], estimate=beta[m],
                  probability=w[m])[, .(true=unique(true),estimate=unique(estimate),
                  probability=sum(probability)), by=cell]
 stopifnot(nrow(bv)==uniqueN(vkey))
 bv[, `:=`(n=n, seed=seed, fold=fold, time=time, source=source)]
 av <- data.table(cell=vkey, probability=w[m])[,.(probability=sum(probability)),by=cell]
 av[, `:=`(true=unname(phi[cell]), estimate=unname(adj[cell]),
           n=n,seed=seed,fold=fold,time=time,source=source)]
 summary <- data.table(n=n,seed=seed,fold=fold,time=time,outcome_time=time+1L,
  source=source, beta_representation_error=beta_error,
  lambda_representation_error=if(time==1L) NA_real_ else sum(w*lambda)-theta[time],
  bridge_correction=correction,bridge_remainder=paired,product_identity_error=paired-product,
  beta_mean_error=sum(w[m]*(beta[m]-truth$beta[m,time]))/sum(w[m]),
  beta_rmse=sqrt(sum(w[m]*(beta[m]-truth$beta[m,time])^2)/sum(w[m])),
  bridge_equation_rmse=sqrt(sum(w*unname(br[bkey])^2)),
  bridge_equation_max_error=max(abs(br)),
  adjoint_equation_rmse=if(time==1L)NA_real_ else
   sqrt(sum(w[m]*(unname(adj[vkey])-unname(phi[vkey]))^2)/sum(w[m])),
  adjoint_equation_max_error=if(time==1L)NA_real_ else max(abs(adj-phi)))
 list(summary=summary,beta=bv,adjoint=av)
}

paths <- 'results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds'
reference <- readRDS(paths)
summaries <- beta_points <- adjoint_points <- list()
for(path in paths) {
 x <- readRDS(path)
 for(j in seq_along(x$folds)) for(time in seq_len(g$tau)) {
  nu <- x$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',j,sep='/')]]
  stopifnot(identical(x$population,p))
  got <- diagnose(nu$beta[,time],nu$lambda[,time],time,x$design$n,x$design$seed,j,'saved check')
  summaries[[length(summaries)+1L]] <- got$summary
  if(time==2L) {
   beta_points[[length(beta_points)+1L]] <- got$beta
   adjoint_points[[length(adjoint_points)+1L]] <- got$adjoint
  }
 }
}
old_summary <- rbindlist(summaries)
fwrite(old_summary[order(-abs(bridge_remainder))],file.path(out,'saved-fit-diagnostics.csv'))
fit_rank <- rbindlist(list(
 old_summary[,.(n,seed,fold,time,outcome_time,function_name='beta',
  representation_error=beta_representation_error)],
 old_summary[time==2L,.(n,seed,fold,time,outcome_time,function_name='lambda',
  representation_error=lambda_representation_error)]))
fwrite(fit_rank[order(-abs(representation_error))],file.path(out,'individual-fit-ranking.csv'))
worst <- fit_rank[function_name=='beta'][which.max(abs(representation_error))]
print(fit_rank[order(-abs(representation_error))][1:8])
cat('Largest bridge interaction error:\n')
print(old_summary[which.max(abs(bridge_remainder))])
cat('Large check: n=',n_large,' seed=',seed_large,' outer sample=',worst$fold,'\n',sep='')

# Reproduce the selected n=4,000 fit using its actual saved assignments.
reference_path <- file.path(out,'reference-reproduction.rds')
if(!file.exists(reference_path)) {
 set.seed(reference$design$seed); d <- draw_data(reference$design$n,g)
 prepared <- study_task(d,g,folds=reference$folds,learner_folds=reference$learner_folds,
  learner_groups=3L,balance_measurement=TRUE)
 task <- prepared$task; j <- worst$fold; rows <- task$folds[[j]]$training_set
 args <- study_arguments(g); ctrl <- study_control(3L)
 ctrl$.nested_split_function <- lmtp:::make_bridge_nested_folds(task$id[rows],
  d[rows,args$measurement,drop=FALSE],task$learner_folds[[j]])
 original <- d; original$..bridge_row <- seq_len(nrow(d))
 fitted <- lmtp:::bridge_nuisance_fit(original,prepared$encoded,task,j,args$measurement,
  args$health,cm_library(FALSE,TRUE),cm_library(FALSE,FALSE),trt_library,TRUE,ctrl)
 predicted <- predict_bridge_functions(fitted,pt,p,g)
 saved <- reference$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',j,sep='/')]]
 errors <- c(beta=max(abs(predicted$beta-saved$beta)),lambda=max(abs(predicted$lambda-saved$lambda)))
 print(errors); stopifnot(all(errors<1e-8))
 saveRDS(list(n=reference$design$n,seed=reference$design$seed,fold=j,
  prediction_difference=errors,beta=predicted$beta,lambda=predicted$lambda),reference_path)
 rm(d,prepared,task,original,fitted,predicted);gc()
}

checkpoint <- file.path(out,'large-fit.rds')
if(!file.exists(checkpoint)) {
 set.seed(seed_large); d <- draw_data(n_large,g)
 stopifnot(all(d$R1==1L),all(d$R2==1L),all(d$A1_policy %in% 0:3),all(d$A2_policy %in% 0:3))
 prepared <- study_task(d,g,learner_groups=3L,balance_measurement=TRUE)
 task <- prepared$task; j <- worst$fold; rows <- task$folds[[j]]$training_set
 nested <- lmtp:::make_bridge_nested_folds(task$id[rows],d[rows,study_arguments(g)$measurement,drop=FALSE],
                                          task$learner_folds[[j]])
 original <- d; original$..bridge_row <- seq_len(n_large)
 ctrl <- study_control(3L); ctrl$.nested_split_function <- nested
 args <- study_arguments(g)
 trace('fit_bridge_ensemble',where=asNamespace('cmbridge'),print=FALSE,
  tracer=quote(cat('CM bridge ensemble, training people:',length(M),'\n')))
 trace('fit_adjoint_ensemble',where=asNamespace('cmbridge'),print=FALSE,
  tracer=quote(cat('CM adjoint ensemble, training people:',length(M),'\n')))
 started <- proc.time()[3]; warnings <- character()
 fitted <- withCallingHandlers(lmtp:::bridge_nuisance_fit(original,prepared$encoded,task,j,
  args$measurement,args$health,cm_library(FALSE,TRUE),cm_library(FALSE,FALSE),
  trt_library,TRUE,ctrl), warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 untrace('fit_bridge_ensemble',where=asNamespace('cmbridge'))
 untrace('fit_adjoint_ensemble',where=asNamespace('cmbridge'))
 tuning <- rbindlist(lapply(c('beta','adjoint','ratio'),function(kind)
  fit_tuning(fitted[[paste0(kind,'_fits')]],kind,j)),fill=TRUE)
 predicted <- predict_bridge_functions(fitted,pt,p,g)
 compact <- function(model) {
  if(is.null(model))return(NULL)
  list(weights=model$weights,candidate_cv_loss=model$candidate_cv_loss,
   candidates=lapply(model$candidates,function(a)a[c('coefficients','target_levels',
    'instrument_levels','penalty_cv','tuning','moment_loss')]))
 }
 counts <- list()
 for(time in 1:2) {
  H<-task$vars$history('A',time); A<-task$vars$A[[time]]
  key<-cell_key(task$natural[rows,c(H,A),drop=FALSE])
  counts[[time]]<-data.table(cell=key,measured=d[[args$measurement[time]]][rows])[
   ,.(training=.N,measured=sum(measured)),by=cell][,time:=time]
 }
 large <- list(n=n_large,seed=seed_large,fold=j,training_n=length(rows),
  beta=predicted$beta,lambda=predicted$lambda,tuning=tuning,counts=rbindlist(counts),
  beta_fit=compact(fitted$beta_fits[[2]]),adjoint_fit=compact(fitted$adjoint_fits[[2]]),
  folds=task$folds,learner_folds=task$learner_folds,
  warnings=unique(warnings),seconds=proc.time()[3]-started,
  cmbridge_version=as.character(packageVersion('cmbridge')),
  lmtp_version=as.character(packageVersion('lmtp')))
 tmp<-paste0(checkpoint,'.tmp');saveRDS(large,tmp);stopifnot(file.rename(tmp,checkpoint))
} else large <- readRDS(checkpoint)
stopifnot(large$n==n_large,large$seed==seed_large)
got <- diagnose(large$beta[,2],large$lambda[,2],2L,n_large,seed_large,large$fold,'large check')
all_summary <- rbindlist(list(old_summary,got$summary))
fwrite(all_summary,file.path(out,'fit-diagnostics.csv'))
fwrite(large$tuning,file.path(out,'large-penalty-selection.csv'))
fwrite(large$counts,file.path(out,'large-training-counts.csv'))
print(got$summary)
cat('Large check completed in',large$seconds,'seconds.\n')
cat('Warnings:',paste(large$warnings,collapse=' | '),'\n')

bp <- rbindlist(c(beta_points,list(got$beta)))
ap <- rbindlist(c(adjoint_points,list(got$adjoint)))
fwrite(bp,file.path(out,'bridge-predictions.csv'));fwrite(ap,file.path(out,'adjoint-equation-predictions.csv'))
save_plot <- function(plot,name,width=7,height=6) {
 ggsave(file.path(out,'figures',paste0(name,'.png')),plot,width=width,height=height,dpi=180)
 ggsave(file.path(out,'figures',paste0(name,'.pdf')),plot,width=width,height=height)
}
base_plot <- function(points,title,subtitle,xlab='True value from the DGP',ylab='Estimated value')
 ggplot(points,aes(true,estimate))+geom_abline(slope=1,intercept=0,color='#3265a8',linewidth=.65)+
 geom_point(aes(size=probability),alpha=.55,color='#2a423e')+
 scale_size_area(max_size=5,guide='none')+coord_equal()+
 labs(title=title,subtitle=subtitle,x=xlab,y=ylab)+theme_minimal(base_size=12)+
 theme(panel.grid.minor=element_blank())
save_plot(base_plot(got$beta,'Bridge at time 2: estimate versus truth',
 sprintf('n = %s; training n = %s; cross-validated positive L1 penalty',
  format(n_large,big.mark=','),format(large$training_n,big.mark=','))), 'bridge-y-equals-x-large')
selected <- old_summary[time==2L, .SD[which.max(abs(beta_representation_error))],by=.(n,seed)]
select_keys <- paste(selected$n,selected$seed,selected$fold)
comparison <- bp[paste(n,seed,fold) %in% select_keys | source=='large check']
comparison[,label:=paste0('n = ',format(n,big.mark=','),'\nseed ',seed,'; training sample ',fold)]
comparison[,label:=factor(label,levels=unique(comparison[order(n,seed)]$label))]
save_plot(base_plot(comparison,'Bridge at time 2: estimate versus truth',
 'Reference: seed 5203002, n = 4,000, training sample 2; same learner for the large check.')+
 facet_wrap(~label,ncol=2),'bridge-y-equals-x-comparison',11,6)
save_plot(base_plot(got$adjoint,'Adjoint at time 2: conditional equation',
 'Compare conditional means: multiple adjoint solutions can be valid.',
 'True conditional mean from the DGP','Conditional mean of fitted adjoint'),
 'adjoint-equation-y-equals-x-large')
writeLines(c(capture.output(sessionInfo()),paste('n',n_large),paste('seed',seed_large)),file.path(out,'session.txt'))
cat('Graphs and numerical diagnostics written to',out,'\n')
