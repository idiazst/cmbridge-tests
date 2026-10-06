# Audit and graph already saved estimated-model comparisons; no refitting.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
library(ggplot2)
input<-Sys.getenv('BRIDGE_CELL_INPUT');root<-Sys.getenv('BRIDGE_CELL_ROOT')
original<-Sys.getenv('BRIDGE_ORIGINAL_ROOT')
stopifnot(nzchar(input),nzchar(root),nzchar(original))
x<-readRDS(input);config<-x$design
set.seed(config$seed);g<-make_mechanism(config$mechanism,dose_max=3L,visits='missing_both')
d<-draw_data(config$n,g);prepared<-study_task(d,g,folds=3L,learner_groups=3L,balance_measurement=TRUE)
stopifnot(identical(prepared$task$folds,x$folds),identical(prepared$task$learner_folds,x$learner_folds))
H<-prepared$task$vars$history('A',2L);A<-prepared$task$vars$A[[2L]]
rows<-x$folds[[1L]]$training_set
B<-as.matrix(prepared$encoded[rows,c(H,A),drop=FALSE])
configurations<-data.table(scorer=c('Gaussian','Gaussian','Cell','Cell'),
  basis=c('Additive','Joint categories','Additive','Joint categories'),
  path=c(file.path(original,'polynomial'),file.path(original,'joint-categories'),
    file.path(root,'polynomial'),file.path(root,'joint_categories')))
metrics<-weights<-points<-penalties<-selected<-audits<-scores<-pairs<-list()
for(i in seq_len(nrow(configurations))) {
  setting<-configurations[i];m<-readRDS(file.path(setting$path,'bridge-fit.rds'))
  label<-paste(setting$scorer,setting$basis,sep='; ')
  stopifnot(length(m$candidate_failures)==0L,identical(m$fold_id,prepared$task$learner_folds[[1L]]),
    all(is.finite(m$weights)),all(m$weights>=0),abs(sum(m$weights)-1)<1e-10)
  for(name in names(m$candidates)) {
    s<-m$candidates[[name]]$solver
    stopifnot(if(name=='landweber')isTRUE(s$stopped_by_tolerance) else isTRUE(s$converged))
  }
  cv_list<-c(list(m$candidates$sieve_md$penalty_cv),
    unlist(m$fold_penalty_cv,recursive=FALSE))
  checked<-0L
  for(cv in cv_list) {
    if(is.null(cv))next
    chosen<-cv[cv$selected,,drop=FALSE];valid<-cv[is.finite(cv$loss)&!cv$failed,,drop=FALSE]
    stopifnot(nrow(chosen)==1L,chosen$scale>0,!chosen$failed,
      chosen$loss==min(valid$loss),chosen$scale>min(valid$scale),chosen$scale<max(valid$scale))
    checked<-checked+1L
    selected[[length(selected)+1L]]<-data.table(configuration=label,
      fit=if(checked==1L)'final candidate' else paste('learner training sample',checked-1L),
      scale=chosen$scale,loss=chosen$loss,failed_trials=sum(cv$failed),interior=TRUE)
  }
  stopifnot(checked==4L)
  gram_u<-gram_v<-self<-matrix(0,3L,3L,dimnames=list(names(m$weights),names(m$weights)))
  for(f in sort(unique(m$fold_id))) {
    validation<-which(m$fold_id==f);r<-m$cv_residuals[validation,,drop=FALSE]
    z<-B[validation,,drop=FALSE];n<-nrow(r)
    if(setting$scorer=='Gaussian') {
      spec<-cmbridge:::.ensemble_kernel(B[m$fold_id!=f,,drop=FALSE],
        'rbf',m$scoring_kernel$control,1L+f)
      features<-cmbridge:::.ensemble_features(z,spec)
      moment<-crossprod(features,r)
      diagonal<-crossprod(r*sqrt(rowSums(features^2)))
      u<-(crossprod(moment)-diagonal)/(n*(n-1))
      v<-crossprod(moment)/n^2;diagonal<-diagonal/n^2
    } else {
      # Independent loop formula verifies that all self-products were removed.
      index<-split(seq_len(n),cell_key(z));u<-v<-diagonal<-matrix(0,3L,3L)
      for(group in index) {
        residual<-r[group,,drop=FALSE];k<-length(group)
        total<-colSums(residual);product<-tcrossprod(total)
        within<-crossprod(residual)
        v<-v+product/(n*k);diagonal<-diagonal+within/(n*k)
        if(k>1L)u<-u+(product-within)/(n*(k-1L))
      }
      pairs[[length(pairs)+1L]]<-data.table(configuration=label,learner_sample=f,
        rows=n,cells=length(index),singleton_cells=sum(lengths(index)==1L),
        paired_rows=sum(lengths(index)[lengths(index)>1L]))
    }
    stopifnot(max(abs(u-m$fold_gram[[f]]))<1e-10)
    fraction<-n/sum(m$fold_n_scored)
    gram_u<-gram_u+fraction*u;gram_v<-gram_v+fraction*v;self<-self+fraction*diagonal
  }
  difference<-max(abs(gram_u-m$raw_gram));stopifnot(difference<1e-10)
  audits[[length(audits)+1L]]<-data.table(configuration=label,
    raw_gram_reproduction_max_difference=difference,
    candidate_failure_records=length(m$candidate_failures),selected_penalties_checked=checked,
    removed_negative_eigenvalues=m$psd_projection$removed_negative_eigenvalues,
    projection_adjustment_norm=m$psd_projection$adjustment_norm)
  scores[[length(scores)+1L]]<-data.table(configuration=label,candidate=names(m$weights),
    weight=as.numeric(m$weights),raw_u_cv_loss=diag(gram_u),v_cv_loss=diag(gram_v),
    v_self_product_term=diag(self),projected_cv_loss=diag(m$gram))
  table<-fread(file.path(setting$path,'population-bridge-checks.csv'))
  table[,`:=`(configuration=label,scorer=setting$scorer,conditioning_basis=setting$basis)]
  metrics[[i]]<-table
  table<-fread(file.path(setting$path,'ensemble-weights.csv'));table[,configuration:=label];weights[[i]]<-table
  table<-fread(file.path(setting$path,'bridge-plot-points.csv'));table[,configuration:=label];points[[i]]<-table
  table<-as.data.table(m$candidates$sieve_md$penalty_cv);table[,configuration:=label];penalties[[i]]<-table
}
for(pair in list(list(metrics,'population-bridge-checks.csv'),list(weights,'ensemble-weights.csv'),
  list(selected,'selected-penalty-audit.csv'),list(audits,'selection-audit.csv'),
  list(scores,'candidate-selection-scores.csv'),list(pairs,'cell-pairs.csv'),list(penalties,'penalty-curves.csv')))
  fwrite(rbindlist(pair[[1L]]),file.path(root,pair[[2L]]))
points<-rbindlist(points)
points[,candidate:=factor(candidate,levels=c('ensemble','sieve_md','landweber','pmmr'),
  labels=c('Ensemble','Sieve minimum distance','Landweber','PMMR'))]
points[,configuration:=factor(configuration,levels=paste(configurations$scorer,configurations$basis,sep='; '))]
dir.create(file.path(root,'figures'),showWarnings=FALSE)
limits<-range(points$reference,points$fitted)
figure<-ggplot(points,aes(reference,fitted,size=probability))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(alpha=.55,colour='#215d80')+facet_grid(configuration~candidate)+
  coord_equal(xlim=limits,ylim=limits)+scale_size_continuous(range=c(.6,3),guide='none')+
  labs(x='Population bridge solution',y='Estimated bridge',
    title='Binary treatment, n = 4,000, replication 2, first training sample',
    subtitle='Outcome at time 3; R2 = 1 and R3 = 1',
    caption='All combinations are shown on equal axes, without clipping. Point size reflects probability.')+
  theme_bw(base_size=10)+theme(strip.text.y=element_text(size=9))
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('all-candidates-and-scores.',ext)),
  figure,width=14,height=14,dpi=180)
ensemble<-points[candidate=='Ensemble'];limits<-range(ensemble$reference,ensemble$fitted)
figure<-ggplot(ensemble,aes(reference,fitted,size=probability))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(alpha=.65,colour='#215d80')+facet_wrap(~configuration,ncol=2)+
  coord_equal(xlim=limits,ylim=limits)+scale_size_continuous(range=c(.8,3.5),guide='none')+
  labs(x='Population bridge solution',y='Estimated ensemble bridge',
    title='Effect of fitting equations and validation scores in the same training sample',
    caption='All ensemble predictions are shown on equal axes. This is one diagnostic, not a coverage study.')+
  theme_bw(base_size=11)
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('ensemble-comparison.',ext)),
  figure,width=9,height=8.5,dpi=180)
curves<-rbindlist(penalties)
figure<-ggplot(curves,aes(scale,loss))+
  geom_line()+geom_point(data=curves[selected==TRUE],colour='#b63a3a',size=2.5)+
  scale_x_log10()+facet_wrap(~configuration,scales='free_y',ncol=2)+
  labs(x='Positive penalty',y='Raw U-statistic validation loss',
    title='Sieve penalty selection within the same training sample',
    caption='Selected penalties are red. Loss scales differ across scoring criteria; vertical values are not directly comparable.')+
  theme_bw(base_size=11)
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('penalty-curves.',ext)),
  figure,width=10,height=7,dpi=180)
cat('Verified all four saved score matrices and all sixteen positive interior sieve selections; no refits.\n')
