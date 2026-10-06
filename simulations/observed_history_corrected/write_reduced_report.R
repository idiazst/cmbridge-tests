library(data.table)
outdir<-Sys.getenv('SIM_OUTPUT','results/observed-history-corrected/missing-both-study-v2')
s<-fread(file.path(outdir,'summary.csv'));design<-fread(file.path(outdir,'design.csv'))
log<-fread(file.path(outdir,'job_log.csv'));details<-readRDS(file.path(outdir,'design-details.rds'))
full<-nrow(design)==1200L && all(design[, .N,by=.(mechanism,n)]$N==200L) &&
 setequal(design$n,c(500L,1000L,4000L)) && setequal(design$mechanism,c('binary_longitudinal','discrete_dose'))
complete<-full && nrow(log)==nrow(design)
stopped<-file.exists(file.path(outdir,'run-stopped-at.txt'))
status<-if(stopped)'**The study was stopped at the user\'s request. Saved results are preliminary.**' else
 if(complete)'The planned datasets have finished. Failed fits and warnings are retained and counted.' else
 if(full)'**The study is running. These are preliminary results.**' else
 '**Implementation and timing check only; this is not the 200-replication study.**'
lines<-c('# Observed-history simulation results','',status,'',
 sprintf('Completed %s of %s planned datasets.',format(nrow(log),big.mark=','),format(nrow(design),big.mark=',')),'',
 sprintf('Both two-time mechanisms use n = %s, with %d replications per sample size and %d parallel R workers. The numerical dose takes values 0, 1, 2, and 3, with the policy increasing it by one, capped at 3.',
 paste(format(sort(details$sizes),big.mark=','),collapse=', '),details$replications,details$workers),'',
 'Bridge candidates are sieve minimum distance with fixed intercept and all main effects plus joint-category deviations, preconditioned Landweber with fixed intercept and all main effects, and PMMR with Gaussian features. The first two inverse-expit classes contain a valid bridge solution for both mechanisms; all three bridge candidates use the inverse-expit parameterization. Adjoint candidates are a saturated joint-category sieve, Landweber with splines, and PMMR, all with an unrestricted identity link. Saturated L1 is excluded from the cmbridge libraries. Outcome regressions combine saturated L1, MARS, and a mean; treatment-ratio classification combines saturated L1 and a mean. MARS is used only for outcome regressions. SDR and TMLE share generated data, sample assignments, and bridge and adjoint fits. Both receive the one-step bridge correction. Positive sieve ridge and L1 penalties are chosen by cross-validation inside the fitting samples. This study does not assess robustness under deliberate misspecification.','',
 'Bridge sieve ridge penalizes link coefficients after main-effect column scaling, with category deviations receiving 100 times the main-effect penalty; the intercept is unpenalized. Adjoint sieve ridge penalizes centered function values. Gaussian kernel U-statistics select ensemble weights and sieve penalties, retaining every history and treatment predictor. The ensemble Gram matrix is projected to positive semidefinite before weight optimization. Failed or nonfinite tuning trials cannot be selected, failed candidates receive zero weight, and these events are recorded. A boundary penalty minimum must become interior after extending the grid, or the candidate fails. cmbridge fits are not clipped after fitting; the estimator applies the previously approved wider bridge interval [-100,100].','',
 if(identical(details$visits,'missing_both')) 'Both follow-up health vectors can be missing. Bridge and adjoint functions are estimated for both outcome times. At a missed intermediate visit, treatment follows the deterministic rule and the final outcome equals the observed baseline outcome.' else 'The intermediate outcome is completely observed in this historical design; only the final outcome requires bridge estimation.','',
 'The [frozen simulation plan](source/reports/simulation-plan.md) gives the conditional distributions, true parameters, learners, and sample-splitting construction.','',
 'Coverage uncertainty describes variation due to the finite number of replications; it is not a confidence interval for the causal parameter.','')
table<-function(x){x<-as.data.frame(x);if(!nrow(x))return(character());c(
 paste0('| ',paste(names(x),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
 apply(x,1,function(row)paste0('| ',paste(row,collapse=' | '),' |')),'')}
labels<-c(binary_longitudinal='Two-time binary treatment',discrete_dose='Two-time numerical dose')
component_file<-file.path(outdir,'eif-component-means.csv')
if(file.exists(component_file)) {
 components<-fread(component_file)
 lines<-c(lines,'## Mean of each EIF component','',
  'For each replication, average the two EIF components over the outer validation observations, using the known true parameter in the sequential component. Then average these sample means over successful replications. Their sum is the mean estimation error.','')
 for(h in c(2L,1L))lines<-c(lines,paste0('### Outcome at time ',h+1L),'',
   table(components[horizon==h][order(mechanism,n,estimator),.(
     Mechanism=labels[mechanism],n,Estimator=toupper(estimator),Successful=successful,
     `Mean sequential EIF`=sprintf('%.6f',mean_sequential_eif),
     `Mean bridge/adjoint EIF`=sprintf('%.6f',mean_bridge_adjoint_eif))]))
 lines<-c(lines,'[Component means](eif-component-means.csv) and [means for each replication](eif-component-means-by-replication.csv) retain the numerical values.','')
}
for(h in c(2L,1L)){
 lines<-c(lines,paste0('## Outcome at time ',h+1L),'')
 for(name in names(labels)){
  x<-s[mechanism==name & horizon==h][order(n,estimator)]
  lines<-c(lines,paste0('### ',labels[[name]]),'',table(x[,.(n,Estimator=toupper(estimator),
   Successful=successful,Failed=failed,Pending=pending,Truth=sprintf('%.5f',truth),
   Mean=sprintf('%.5f',mean_estimate),Bias=sprintf('%.5f',bias),SD=sprintf('%.5f',empirical_sd),
   RMSE=sprintf('%.5f',rmse),`Mean SE`=sprintf('%.5f',mean_se),
   `Interval length`=sprintf('%.5f',mean_interval_length),Coverage=sprintf('%.3f',coverage),
   `Coverage uncertainty`=sprintf('[%.3f, %.3f]',coverage_mc_lower,coverage_mc_upper),
   Warnings=replications_with_warnings)]))
 }
}
lines<-c(lines,'## Counterfactual mean curves','',
 'The fraction of simultaneous bands containing both true means is reported separately from pointwise coverage.','',
 table(s[horizon==2L][order(mechanism,n,estimator),.(Mechanism=labels[mechanism],n,
 Estimator=toupper(estimator),Successful=successful,`Joint coverage`=sprintf('%.3f',joint_curve_coverage))]),
 '![Estimated curves and true means](figures/correctly-specified-curves.png)','',
 '## Estimation error and interval performance','')
for(name in names(labels))for(kind in c('bias','root-n-bias','coverage','standard-errors'))
 lines<-c(lines,sprintf('![%s: %s](figures/%s-%s.png)',labels[[name]],kind,name,kind),'')
lines<-c(lines,'## Results requiring examination','')
flag<-s[horizon==2L & n==max(design$n) & successful>=100L &
 (coverage_mc_upper<.95 | abs(bias)>1.96*bias_mcse)]
if(nrow(flag))lines<-c(lines,'These results require examination. A flag alone does not demonstrate an implementation error. The fitted conditional equations and separate population remainder contributions must be checked.','',
 table(flag[,.(mechanism,estimator,n,successful,bias,coverage,coverage_mc_upper)])) else
 lines<-c(lines,if(complete)'No largest-sample result met this automated screening criterion. This does not establish that the asymptotic conditions hold.'else
 'The planned replications are incomplete. Coverage and convergence conclusions remain preliminary.','')
lines<-c(lines,'## Computation and validation','',
 sprintf('Recorded dataset processing times sum to %.2f hours. This is the sum across workers, not elapsed wall time.',sum(log$seconds)/3600),'',
 'For every successful fit, the direct prediction-based one-step calculation agrees with the public package output. Exact population calculations check the separate remainder contributions. Unknown nuisance functions are fitted from generated data; true functions are used only for diagnostics.','',
 '[Numerical summary](summary.csv), [replication estimates](replicates.csv), [population diagnostics](population_diagnostics.csv), [penalties and weights](penalty_and_weight_selections.csv), [errors](errors.csv), and [job log](job_log.csv) retain the numerical records. Checkpoints also store predictor-combination counts, fitted functions on the finite support, seeds, and sample assignments.','',
 '[Failed nested candidates](nested-candidate-failures.csv) and [failed penalty trials](nested-penalty-trial-failures.csv) retain recorded bridge fitting failures even when the remaining candidates produce a complete estimator. These events are distinct from a failed final estimator.','',
 'The bridge and adjoint may have multiple solutions; conditional-equation errors and the actual remainder determine their accuracy. Difference from one selected exact solution alone is not treated as misspecification.','')
timing_path<-file.path(outdir,'execution-time-estimate.csv')
if(file.exists(timing_path)){
 timing<-fread(timing_path)
 if(!complete && !anyNA(timing$mean_seconds)) lines<-c(lines,
  sprintf('Based on observed processing times in each mechanism and sample size, the remaining fits are estimated to require approximately %.1f hours with %d workers. This estimate includes contention measured during this run but excludes final analysis and is revised as more datasets finish.',
   sum(timing$estimated_remaining_worker_hours)/details$workers,details$workers),'')
}
if(file.exists(file.path(outdir,'ensemble_weight_summary.csv'))) lines<-c(lines,
 '[All candidate weights](ensemble_weights.csv) and [weight summaries](ensemble_weight_summary.csv) retain the expanded library, including candidates selected with zero weight.','')
diagnostic_file<-file.path(outdir,'population_diagnostic_summary.csv')
if(file.exists(diagnostic_file)){
 dg<-fread(diagnostic_file)
 ds<-dcast(dg[horizon==2L & variable %in% c('population_bias','sequential_remainder','bridge_remainder')],
  mechanism+n+estimator~variable,value.var='mean')
 lines<-c(lines,'## Conditional equations and estimation error','',
 'For each outer training sample, evaluate its fitted functions over the complete generating distribution. The population error below is the expectation of its estimator contribution minus the true parameter. Average these expectations using the outer validation sample sizes, then average across completed replications. The sequential-regression and bridge contributions add to the population error. They describe fitted-function estimation error separately from variation in the final evaluation sample.','',
 table(ds[order(mechanism,n,estimator),.(Mechanism=labels[mechanism],n,Estimator=toupper(estimator),
  `Population error`=sprintf('%.5f',population_bias),
  `Sequential-regression contribution`=sprintf('%.5f',sequential_remainder),
  `Bridge contribution`=sprintf('%.5f',bridge_remainder))]),
  'The libraries contain the required functions. Coverage additionally depends on how accurately those functions are estimated. Sparse full-history combinations and regularization may affect finite-sample performance. The study retains this behavior and does not alter the generating distributions or learner settings in response to coverage.','')
}
writeLines(lines,file.path(outdir,'REPORT.md'))
if(requireNamespace('markdown',quietly=TRUE))markdown::mark_html(file.path(outdir,'REPORT.md'),output=file.path(outdir,'REPORT.html'))
cat('Report written to',file.path(outdir,'REPORT.md'),'\n')
