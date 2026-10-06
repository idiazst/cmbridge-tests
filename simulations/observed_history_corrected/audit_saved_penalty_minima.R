# Audit exact saved double-precision losses; do not refit or select a new penalty.
library(data.table)
input<-Sys.getenv('DIAGNOSTIC_INPUT');output<-Sys.getenv('DIAGNOSTIC_OUTPUT')
stopifnot(nzchar(input),nzchar(output))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
records<-list()
paths<-list.files(file.path(input,'jobs'),'\\.rds$',full.names=TRUE)
stopifnot(length(paths)>0L)
for(path in paths) {
  x<-readRDS(path);tuning<-x$tuning
  if(is.null(tuning)||!nrow(tuning))next
  keys<-intersect(c('kind','fold','horizon_or_depth','candidate','constant',
    'beta_correct','ratio_correct','m_correct','estimator'),names(tuning))
  audit<-tuning[,{
    cv<-.SD[order(scale)]
    selected<-which(cv$selected)
    stopifnot(length(selected)==1L)
    chosen<-selected[1L];successful<-which(!cv$failed & is.finite(cv$loss))
    stopifnot(length(successful)>0L,chosen%in%successful)
    minimum<-min(cv$loss[successful]);minima<-which(cv$loss==minimum)
    list(scale=cv$scale[chosen],selected_loss=cv$loss[chosen],
      selected_is_exact_minimum=cv$loss[chosen]==minimum,
      selected_is_interior=chosen>1L && chosen<nrow(cv),
      selected_at_successful_grid_edge=chosen%in%range(successful),
      left_neighbour_failed=if(chosen>1L)cv$failed[chosen-1L]else TRUE,
      right_neighbour_failed=if(chosen<nrow(cv))cv$failed[chosen+1L]else TRUE,
      exact_minimum_count=length(minima),
      minimum_also_at_lower_boundary=1L%in%minima,
      minimum_also_at_upper_boundary=nrow(cv)%in%minima,
      grid_points=nrow(cv),failed_trials=sum(cv$failed),
      penalty_status=cv$penalty_status[chosen])
  },by=keys]
  audit[,`:=`(mechanism=x$design$mechanism,n=x$design$n,
    replicate=x$design$replicate,seed=x$design$seed)]
  stopifnot(all(audit$scale>0),all(is.finite(audit$selected_loss)),
    all(audit$selected_is_exact_minimum),all(audit$selected_is_interior),
    !any(audit$selected_at_successful_grid_edge),
    !any(audit$left_neighbour_failed),!any(audit$right_neighbour_failed),
    all(startsWith(audit$penalty_status,'interior')))
  records[[length(records)+1L]]<-audit
}
records<-rbindlist(records,use.names=TRUE,fill=TRUE)
stopifnot(nrow(records)>0L)
fwrite(records,file.path(output,'saved-penalty-minima.csv'))
print(records[,.(selected_records=.N,
  exact_ties=sum(exact_minimum_count>1L),
  ties_also_at_boundary=sum(minimum_also_at_lower_boundary|minimum_also_at_upper_boundary),
  successful_grid_edges=sum(selected_at_successful_grid_edge),
  adjacent_failed_trials=sum(left_neighbour_failed|right_neighbour_failed)),by=kind])
cat('Saved penalty audit complete; no trial refits or source changes.\n')
