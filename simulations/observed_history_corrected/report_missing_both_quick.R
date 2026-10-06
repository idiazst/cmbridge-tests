library(data.table)
out <- Sys.getenv('SIM_OUTPUT','results/observed-history-corrected/missing-both-n10000')
e <- fread(file.path(out,'estimates.csv'))
s <- fread(file.path(out,'eif-summary.csv'))
p <- fread(file.path(out,'population-remainders.csv'))
missing <- fread(file.path(out,'sample-missingness.csv'))
stopifnot(nrow(e)==8L,nrow(s)==8L,nrow(p)==8L)
v <- merge(e,s,by=c('mechanism','estimator','horizon'))
setorder(v,mechanism,horizon,estimator)
v[,example:=ifelse(mechanism=='binary_longitudinal','Binary','Dose')]
v[,outcome_time:=horizon+1L]
v[,`:=`(sequential_eif_mean_at_truth=sequential_estimate-truth,
  total_eif_mean_at_truth=final_estimate-truth)]
stopifnot(max(abs(v$estimate-v$final_estimate))<1e-10,
  max(abs(v$se-v$total_se))<1e-10,
  max(abs(v$sequential_estimate+v$bridge_adjoint_correction-v$estimate))<1e-10)
fwrite(v,file.path(out,'combined-results.csv'))
means<-v[,.(example,outcome_time,estimator=toupper(estimator),
  mean_sequential_eif=sequential_eif_mean_at_truth,
  mean_bridge_adjoint_eif=bridge_adjoint_eif_mean)]
fwrite(means,file.path(out,'eif-component-means.csv'))
components<-rbindlist(lapply(file.path(out,'jobs',paste0(c('binary_longitudinal','discrete_dose'),'.rds')),
  function(path)readRDS(path)$eif_components))
if('truth' %in% names(components))components[,truth:=NULL]
components<-merge(components,unique(v[,.(mechanism,estimator,horizon,truth)]),
  by=c('mechanism','estimator','horizon'),all.x=TRUE,sort=FALSE)
components[,`:=`(sequential_eif_at_truth=sequential_score-truth,
  total_eif_at_truth=sequential_score+bridge_adjoint-truth)]
stopifnot(nrow(components)==80000L,!anyNA(components$truth),
  max(abs(components$sequential_eif_at_truth+components$bridge_adjoint-
    components$total_eif_at_truth))<1e-12)
fwrite(components,file.path(out,'person-eif-contributions.csv.gz'))
if(Sys.getenv('REPORT_PLOTS','1')=='1' && requireNamespace('ggplot2',quietly=TRUE)) {
  library(ggplot2)
  dir.create(file.path(out,'figures'),showWarnings=FALSE)
  plotdata<-copy(v);plotdata[,estimator:=toupper(estimator)]
  truth<-unique(plotdata[,.(example,outcome_time,truth)])
  fig<-ggplot(plotdata,aes(factor(outcome_time),estimate,colour=estimator))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=.12,position=position_dodge(.4))+
    geom_point(size=2.6,position=position_dodge(.4))+
    geom_point(data=truth,aes(factor(outcome_time),truth),inherit.aes=FALSE,shape=18,size=3)+
    facet_wrap(~example)+labs(x='Outcome time',y='Estimate and 95% interval',colour='Estimator',
      caption='Black diamonds: exact true values. One dataset per example; these intervals do not assess coverage.')+
    theme_bw(base_size=12)
  ggsave(file.path(out,'figures/estimates.png'),fig,width=8,height=4.5,dpi=160)
  sddata<-melt(plotdata,id.vars=c('example','outcome_time','estimator'),
    measure.vars=c('sequential_eif_sd','bridge_adjoint_eif_sd','total_eif_sd'),
    variable.name='component',value.name='sd')
  sddata[,component:=factor(component,levels=c('sequential_eif_sd','bridge_adjoint_eif_sd','total_eif_sd'),
    labels=c('Sequential regression','Bridge/adjoint','Total'))]
  fig<-ggplot(sddata,aes(factor(outcome_time),sd,fill=component))+
    geom_col(position='dodge')+facet_grid(estimator~example)+
    labs(x='Outcome time',y='Empirical EIF standard deviation',fill='Contribution',
      caption='The total includes covariance between contributions; their standard deviations do not add.')+
    theme_bw(base_size=12)+theme(legend.position='bottom')
  ggsave(file.path(out,'figures/eif-contributions.png'),fig,width=8,height=6,dpi=160)
}
table <- function(x,cols,labels,digits=6L) {
  y <- x[,..cols]
  for(j in seq_along(y))if(is.numeric(y[[j]]) && !is.integer(y[[j]]))
    set(y,j=j,value=formatC(ifelse(abs(y[[j]])<.5*10^(-digits),0,y[[j]]),digits=digits,format='f'))
  for(j in seq_along(y))set(y,j=j,value=gsub('[\r\n]+',' ',gsub('|',' / ',as.character(y[[j]]),fixed=TRUE)))
  c(paste0('| ',paste(labels,collapse=' | '),' |'),
    paste0('| ',paste(rep(':--',length(cols)),collapse=' | '),' |'),
    apply(as.data.frame(y),1L,function(z)paste0('| ',paste(z,collapse=' | '),' |')))
}
mean_lines<-c('# Mean of each EIF component','',
  'Sample means for the n = 10,000 datasets. The sequential component uses the known true parameter in the paper’s expression. All fitted functions are evaluated on their outer validation observations.','',
  table(means,c('example','outcome_time','estimator','mean_sequential_eif','mean_bridge_adjoint_eif'),
    c('Example','Outcome time','Estimator','Mean sequential EIF','Mean bridge/adjoint EIF')),'',
  '[Numerical values](eif-component-means.csv).')
writeLines(mean_lines,file.path(out,'EIF-COMPONENT-MEANS.md'))
text <- c('# One-dataset check with missingness at both follow-ups','',
  'One independently generated n = 10,000 dataset was used for each example: binary seed 6100001 and numerical-dose seed 6200001. SDR and TMLE used the same dataset, outer assignments, learner assignments, bridge fits, and adjoint fits within each example. The repeated-sample simulation remains stopped.','',
  '## Mean of each EIF component','',
  'These are sample means over the outer validation observations, using the known true parameter in the sequential EIF component.','',
  table(means,c('example','outcome_time','estimator','mean_sequential_eif','mean_bridge_adjoint_eif'),
    c('Example','Outcome time','Estimator','Mean sequential EIF','Mean bridge/adjoint EIF')),'',
  '[Mean table](EIF-COMPONENT-MEANS.md) and [numerical values](eif-component-means.csv).','',
  '## Generating distributions and estimators','',
  'All four function classes are correctly specified: saturated joint-category sieve and cross-validated saturated L1 are included in the bridge/adjoint ensemble, alongside Landweber and PMMR. Sequential regressions combine saturated cross-validated L1, MARS, and the mean; treatment ratios combine saturated cross-validated L1 and the mean. No true nuisance function is supplied to fitting. Positive L1 penalties are selected within training observations.','',
  'At a missed intermediate visit, treatment is (0, 0), and the final outcome equals the observed baseline outcome. The final covariate remains random. This is the revised distribution audited against the paper, not the earlier completely observed intermediate design.','',
  '## Estimates and contributions to the estimating equation','',
  'The sequential column is the sample mean of the uncentered sequential-regression part of equation (7). The bridge/adjoint column is the sample mean of the correction below. Their sum is the final estimate.','',
  '$$',
  'D_t^{\\mathrm{br}}(O)=-\\lambda_t(A_t,H_t)',
  '\\{R_{t+1}\\beta_t(H_t,C_{t+1})-1\\}.',
  '$$','',
  table(v,c('example','outcome_time','estimator','truth','sequential_estimate',
    'bridge_adjoint_correction','estimate','se','lower','upper'),
    c('Example','Outcome time','Estimator','True value','Sequential','Bridge/adjoint','Final estimate','SE','CI lower','CI upper')),'',
  '![Estimates and exact true values](figures/estimates.png)','',
  '## Contributions to the estimated EIF','',
  'These are the two components defined in Section 3.1 of the paper. The paper establishes their sum as a gradient; it does not establish efficiency in the observed-data model. The EIF labels here follow the requested terminology.','',
  '$$',
  '\\begin{aligned}',
  'D_t^{\\mathrm{seq}}(O)',
  '&=q_{t+1,1}(A_1,H_1)-\\theta_{t+1}(\\bar d_t)',
  '+\\sum_{s=1}^{t}\\Omega_s',
  '\\{\\eta_{t+1,s+1}(O)-m_{t+1,s}(A_s,H_s)\\},\\\\',
  'D_t(O)&=D_t^{\\mathrm{seq}}(O)+D_t^{\\mathrm{br}}(O).',
  '\\end{aligned}',
  '$$','',
  'For this diagnostic, the sequential component subtracts the known true parameter. The bridge/adjoint component retains the expression above. The means of these two components add to the observed estimation error.','',
  'The mean of each component is reported at the beginning of this document.','',
  'For confidence intervals, the sequential component instead subtracts the final corrected estimate. Its mean then cancels the bridge/adjoint mean, giving total mean zero by construction. That cancellation is not evidence that either conditional equation is accurately estimated. Centering by either value gives the same component variances and covariance:','',
  table(v,c('example','outcome_time','estimator','sequential_eif_sd',
    'bridge_adjoint_eif_sd','eif_covariance','total_eif_sd'),
    c('Example','Outcome time','Estimator','Sequential SD','Bridge/adjoint SD','Covariance','Total SD')),'',
  'The standard error uses the variance of the sum, including twice the covariance. It does not add the two component standard errors. The [person-level file](person-eif-contributions.csv.gz) retains all 80,000 contributions (two examples × two estimators × two outcomes × 10,000 people).','',
  '![Variability of EIF contributions](figures/eif-contributions.png)','',
  '## Exact population estimation errors for the fitted functions','',
  'These are the separate second-order remainder contributions evaluated by exact summation over the generating distribution, averaged over the outer training fits with validation-sample weights. They are different quantities from the empirical EIF means above. The sequential-regression and bridge/adjoint columns sum to the population estimation error.','')
p[,example:=ifelse(mechanism=='binary_longitudinal','Binary','Dose')]
p[,outcome_time:=horizon+1L];setorder(p,mechanism,horizon,estimator)
details<-character()
if(file.exists(file.path(out,'bridge-remainder-by-measurement.csv'))) {
  b<-fread(file.path(out,'bridge-remainder-by-measurement.csv'))
  dose<-b[mechanism=='discrete_dose']
  weights<-fread(file.path(out,'final-outcome-cmbridge-weights.csv'))
  l1<-weights[mechanism=='discrete_dose' & kind=='beta' & candidate=='saturated_cv',mean_weight]
  details<-c('','### Final numerical-dose outcome','',
    sprintf('The bridge/adjoint remainder is %.6f. Its contribution from histories with a measured intermediate visit is %.6f (%.1f%% of the total); the contribution from missed intermediate visits is %.6f. Thus most of this error occurs in the measured-visit histories, where the outcome remains random.',
      sum(dose$bridge_remainder),dose[R2==1L,bridge_remainder],
      100*dose[R2==1L,bridge_remainder]/sum(dose$bridge_remainder),dose[R2==0L,bridge_remainder]),'',
    sprintf('The final-outcome bridge assigns %.1f%% average ensemble weight to saturated L1. Every outer training sample selected penalty scale 0.1, the largest candidate scale. These are descriptions of the saved fits; they do not establish that the penalty caused the remainder. Full conditional-equation errors and weights remain available in the numerical files.',100*l1),'',
    'The generated-sample estimates are close to the true mean, but evaluation on those observations is about +0.0244 (SDR) and +0.0243 (TMLE) above the corresponding fitted-function population expectations. This offsets the negative population remainder in this dataset. The point estimates alone therefore do not establish negligible second-order error.','',
    '[Remainders by intermediate measurement](bridge-remainder-by-measurement.csv) and [final-outcome candidate weights](final-outcome-cmbridge-weights.csv) were computed from the saved fits, without refitting.')
}
text <- c(text,table(p,c('example','outcome_time','estimator','sequential_remainder',
  'bridge_remainder','population_bias'),c('Example','Outcome time','Estimator',
    'Sequential regression','Bridge/adjoint','Sum')),details,'',
  '## Missingness and checks','',table(missing,c('mechanism','time','sample_missing_percent'),
    c('Example','Time','Sample missing (%)'),digits=2L),'',
  'The estimator output agrees with direct person-level evaluation of equation (7). Component means reconstruct every estimate, and component variances and covariance reconstruct every standard error. The exact bridge remainder agrees with its product identity. The generating-distribution audits are preserved in [validation/dgp](validation/dgp/REPORT.md).','',
  'One dataset per example assesses implementation and displays contributions. It does not assess confidence-interval coverage or convergence across replications.','',
  'Sources, seeds, datasets, person assignments, population predictions, candidate weights, selected penalties, warnings, and checkpoints are preserved. See [REPRODUCE.md](REPRODUCE.md).')
warning_rows <- unique(e[nzchar(warnings),.(mechanism,estimator,warnings)])
if(nrow(warning_rows))text<-c(text,'','## Recorded warnings','',
  table(warning_rows,c('mechanism','estimator','warnings'),c('Example','Estimator','Warnings')))
writeLines(text,file.path(out,'REPORT.md'))
if(requireNamespace('markdown',quietly=TRUE))markdown::mark_html(
  file.path(out,'REPORT.md'),output=file.path(out,'REPORT.html'))
cat('Quick-test report written.\n')
