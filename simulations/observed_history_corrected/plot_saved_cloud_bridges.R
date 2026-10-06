# Plot stored ensemble predictions from the completed cloud n=4000 datasets.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
library(ggplot2)
root<-Sys.getenv('SIM_OUTPUT');stopifnot(nzchar(root))
points<-metrics<-list()
for(path in list.files(file.path(root,'jobs'),'\\.rds$',full.names=TRUE)) {
  x<-readRDS(path);if(x$design$n!=4000L||!length(x$fitted_functions))next
  g<-make_mechanism(x$design$mechanism,dose_max=3L,visits='missing_both')
  p<-x$population;truth<-population_truth_functions(p,g)
  for(fold in seq_along(x$folds)) {
    nu<-x$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',fold,sep='/')]]
    if(is.null(nu))next
    selected<-p$R2==1 & p$R3==1
    points[[length(points)+1L]]<-data.table(mechanism=x$design$mechanism,
      replicate=x$design$replicate,fold=factor(fold),probability=p$probability[selected],
      reference=truth$beta[selected,2],fitted=nu$beta[selected,2])
    metrics[[length(metrics)+1L]]<-data.table(mechanism=x$design$mechanism,
      replicate=x$design$replicate,seed=x$design$seed,fold=fold,
      rmse=sqrt(sum(p$probability[selected]*(nu$beta[selected,2]-truth$beta[selected,2])^2)/
        sum(p$probability[selected])),
      boundary_fraction=mean(nu$beta[selected,2]<=-100 | nu$beta[selected,2]>=100))
  }
}
if(!length(points))stop('No n=4000 fitted functions saved yet.')
points<-rbindlist(points);metrics<-rbindlist(metrics)
fwrite(metrics,file.path(root,'saved-cloud-bridge-rmse.csv'))
dir.create(file.path(root,'figures'),showWarnings=FALSE)
for(wanted_mechanism in unique(points$mechanism)) {
  selected<-points[mechanism==wanted_mechanism]
  selected[,replication_label:=factor(paste('Replication',replicate),
    levels=paste('Replication',sort(unique(replicate))))]
  limits<-range(selected$reference,selected$fitted)
  plot<-ggplot(selected,aes(reference,fitted,colour=fold))+
    geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
    geom_point(alpha=.4,size=.8)+facet_wrap(~replication_label,ncol=3)+
    coord_equal(xlim=limits,ylim=limits)+
    labs(x='Valid bridge solution',y='Stored fitted ensemble bridge',colour='Training sample',
      title=paste(if(wanted_mechanism=='binary_longitudinal')'Binary treatment'else'Numerical dose',
        ', n = 4,000, outcome at time 3',sep=''),
      caption='All saved predictions with R2 = 1 and R3 = 1 are shown. Equal axes; no cropped outliers.')+
    theme_bw(base_size=11)+theme(legend.position='bottom')
  for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0(wanted_mechanism,'-saved-n4000-bridges.',ext)),
    plot,width=10,height=3.5*ceiling(uniqueN(selected$replicate)/3)+1,dpi=180)
}
print(metrics)
