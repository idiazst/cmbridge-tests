# Scientific figures from actual saved cloud fits; no nuisance models are fitted.
library(data.table)
library(ggplot2)
root<-Sys.getenv('DIAGNOSTIC_OUTPUT');input<-Sys.getenv('DIAGNOSTIC_INPUT')
stopifnot(nzchar(root),nzchar(input))
dir.create(file.path(root,'figures'),showWarnings=FALSE)
d<-fread(file.path(root,'equation-diagnostics-by-dataset.csv'))[horizon==2]
flags<-d[order(-abs(bridge_remainder))][1:2]
fwrite(flags,file.path(root,'largest-remainder-datasets.csv'))
a<-fread(file.path(root,'adjoint-conditional-equations.csv'))
for(i in seq_len(nrow(flags))) {
  chosen<-flags[i]
  selected<-a[mechanism==chosen$mechanism & n==chosen$n & replicate==chosen$replicate & horizon==2]
  stopifnot(nrow(selected)>0L)
  selected[,fold_label:=paste('Training sample',fold)]
  selected[,training_rows:=factor(ifelse(measured_training_n==0L,'0',
    ifelse(measured_training_n<=2L,'1-2','3 or more')),levels=c('0','1-2','3 or more'))]
  limits<-range(selected$conditional_true_loading,selected$conditional_fitted_adjoint)
  figure<-ggplot(selected,aes(conditional_true_loading,conditional_fitted_adjoint,
    colour=training_rows,size=probability))+
    geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
    geom_point(alpha=.65)+facet_wrap(~fold_label,nrow=1)+
    coord_equal(xlim=limits,ylim=limits)+
    scale_size_continuous(range=c(.65,3.5),guide='none')+
    labs(x='Conditional mean of the true adjoint loading',y='Conditional mean of the fitted adjoint',
      colour='Measured training rows\nsharing conditioning values',
      title=sprintf('Binary treatment, n = %s, replication %s, outcome at time 3',chosen$n,chosen$replicate),
      caption='Exact population conditional equations. All combinations are shown; point size reflects probability.')+
    theme_bw(base_size=11)+theme(legend.position='bottom')
  stem<-sprintf('binary-n%s-r%03d-adjoint-equations',chosen$n,chosen$replicate)
  for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0(stem,'.',ext)),
    figure,width=12,height=5.5,dpi=180)
}
population<-fread(file.path(input,'population_diagnostics.csv'))[horizon==2]
population<-population[,.(sequential_component=sum(fold_weight*sequential_remainder),
  bridge_component=sum(fold_weight*bridge_remainder),population_error=sum(fold_weight*population_bias)),
  by=.(mechanism,n,replicate,seed,estimator)]
stopifnot(max(abs(population$sequential_component+population$bridge_component-population$population_error))<1e-9)
fwrite(population,file.path(root,'population-components-by-dataset.csv'))
population[,estimator_label:=toupper(estimator)]
figure<-ggplot(population,aes(sequential_component,bridge_component,colour=factor(n)))+
  geom_hline(yintercept=0,colour='grey75')+geom_vline(xintercept=0,colour='grey75')+
  geom_abline(slope=-1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(size=2)+facet_wrap(~estimator_label)+coord_equal()+
  labs(x='Sequential component: fixed-function population mean',
    y='Bridge/adjoint component: fixed-function population mean',colour='Sample size',
    title='Outcome at time 3 in the 24 saved cloud datasets',
    caption='The two coordinates sum to population error. The dashed line has zero sum. These are diagnostics, not sample EIF means.')+
  theme_bw(base_size=11)+theme(legend.position='bottom')
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('population-components.',ext)),
  figure,width=10,height=5.5,dpi=180)
class_values<-melt(d,id.vars=c('n','replicate'),
  measure.vars=c('adjoint_equation_rmse','observed_dictionary_adjoint_class_minimum_equation_rmse'),
  variable.name='quantity',value.name='rmse')
class_values[,quantity:=factor(quantity,levels=c('adjoint_equation_rmse',
  'observed_dictionary_adjoint_class_minimum_equation_rmse'),
  labels=c('Actual ensemble equation RMSE','Minimum within observed sieve dictionary'))]
figure<-ggplot(class_values,aes(factor(n,levels=c(500,1000,4000)),rmse,colour=quantity))+
  geom_point(size=2,position=position_jitter(width=.07,height=0,seed=1))+
  labs(x='Sample size',y='Adjoint conditional-equation RMSE',colour=NULL,
    title='Actual fits and observed adjoint dictionaries at time 3',
    caption='Each point averages three training samples. The dictionary minimum is not a lower bound for the whole ensemble.')+
  theme_bw(base_size=11)+theme(legend.position='bottom')
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('adjoint-equations-and-dictionaries.',ext)),
  figure,width=10,height=5.5,dpi=180)
cat('Saved cloud equation figures without refitting.\n')
