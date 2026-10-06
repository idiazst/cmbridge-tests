# Scientific graphs from saved fitted-function evaluations, without refitting.
library(data.table)
library(ggplot2)
root<-Sys.getenv('DIAGNOSTIC_OUTPUT');stopifnot(nzchar(root))
a<-fread(file.path(root,'adjoint-conditional-equations.csv'))
dir.create(file.path(root,'figures'),showWarnings=FALSE)
selected<-a[mechanism=='binary_longitudinal' & n==1000 & replicate==1 & horizon==2]
stopifnot(nrow(selected)>0L)
selected[,fold_label:=paste('Training sample',fold)]
selected[,training_rows:=factor(ifelse(measured_training_n==0L,'0',
  ifelse(measured_training_n<=2L,'1-2','3 or more')),levels=c('0','1-2','3 or more'))]
limits<-range(selected$conditional_true_loading,selected$conditional_fitted_adjoint)
p<-ggplot(selected,aes(conditional_true_loading,conditional_fitted_adjoint,
  colour=training_rows,size=probability))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(alpha=.7)+facet_wrap(~fold_label,nrow=1)+
  coord_equal(xlim=limits,ylim=limits)+
  scale_size_continuous(range=c(.8,4),guide='none')+
  labs(x='Conditional mean of the true adjoint loading',y='Conditional mean of the fitted adjoint',
    colour='Measured training rows\nsharing conditioning values',
    title='Binary treatment, n = 1,000, replication 1, outcome at time 3',
    caption='Exact finite-population equations. All conditioning-value combinations are shown; point size reflects probability.')+
  theme_bw(base_size=11)+theme(legend.position='bottom')
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('binary-n1000-r001-adjoint-equations.',ext)),
  p,width=12,height=5.5,dpi=180)
d<-fread(file.path(root,'equation-diagnostics-by-dataset.csv'))[horizon==2]
d[,mechanism_label:=c(binary_longitudinal='Binary treatment',discrete_dose='Numerical dose')[mechanism]]
q<-ggplot(d,aes(factor(n,levels=c(500,1000,4000)),
  observed_dictionary_adjoint_class_minimum_equation_rmse))+
  geom_point(size=2,position=position_jitter(width=.07,height=0,seed=1))+
  facet_wrap(~mechanism_label)+
  labs(x='Sample size',y='Smallest attainable adjoint equation RMSE',
    title='Adjoint class defined by each measured training dictionary',
    caption='Exact population projection, allowing any valid solution. Each point averages the three training samples of a saved dataset; no nuisance fits are run.')+
  theme_bw(base_size=11)
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('observed-dictionary-adjoint-class.',ext)),
  q,width=10,height=5,dpi=180)
cat('Saved defining-equation and class-audit figures.\n')
