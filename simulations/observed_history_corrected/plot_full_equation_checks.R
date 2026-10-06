# Graph captured functions/equations without refitting or changing frozen inputs.
library(data.table)
library(ggplot2)
root<-Sys.getenv('SIM_OUTPUT');stopifnot(nzchar(root))
folders<-list.dirs(root,recursive=FALSE,full.names=TRUE)
folders<-folders[grepl('-n4000-seed',basename(folders))]
read_captures<-function(stem)rbindlist(lapply(folders,function(folder) {
  paths<-list.files(folder,paste0('^',stem,'-fold[1-3][.]csv$'),full.names=TRUE)
  if(length(paths)!=3L)stop('All three outer captures are required for these figures.')
  table<-rbindlist(lapply(paths,fread),fill=TRUE)
  table[,`:=`(mechanism=sub('-n.*','',basename(folder)),n=4000L)]
  table
}),fill=TRUE)
points<-read_captures('predictions');equations<-read_captures('equations')
metrics<-points[,.(rmse=sqrt(sum(probability*(estimate-true)^2)/sum(probability)),
  signed_error=sum(probability*(estimate-true))/sum(probability),
  maximum_absolute_error=max(abs(estimate-true)),minimum=min(estimate),maximum=max(estimate),
  combinations=.N,unique_solution=all(unique_solution)),
  by=.(mechanism,n,kind,horizon,current_R,fold,candidate)]
eq<-equations[,.(equation_rmse=sqrt(sum(probability*(estimate-true)^2)/sum(probability)),
  equation_max_error=max(abs(estimate-true))),
  by=.(mechanism,n,kind,horizon,current_R,fold,candidate)]
metrics<-merge(metrics,eq,by=c('mechanism','n','kind','horizon','current_R','fold','candidate'))
fwrite(metrics,file.path(root,'function-and-equation-errors.csv'))
fwrite(points,file.path(root,'all-function-predictions.csv'))
fwrite(equations,file.path(root,'all-equation-predictions.csv'))
labels<-c(ensemble='Ensemble',sieve_md='Sieve minimum distance',landweber='Landweber',pmmr='PMMR')
for(mechanism_name in c('binary_longitudinal','discrete_dose')) {
  name<-if(mechanism_name=='binary_longitudinal')'Binary treatment' else 'Numerical dose'
  selected<-points[mechanism==mechanism_name & kind=='beta' & horizon==2 & current_R==1]
  stopifnot(nrow(selected)>0L,all(selected$unique_solution),
    all(is.finite(selected$estimate)),all(selected$estimate>=1))
  selected[,candidate_label:=factor(labels[candidate],levels=labels)]
  selected[,training_sample:=factor(fold)]
  for(ensemble_only in c(FALSE,TRUE)) {
    table<-if(ensemble_only)selected[candidate=='ensemble'] else selected
    limits<-range(table$true,table$estimate)
    figure<-ggplot(table,aes(true,estimate,colour=training_sample))+
      geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey35')+
      geom_point(alpha=.55,size=1.3)+coord_equal(xlim=limits,ylim=limits)+
      scale_colour_manual(values=c('#2166ac','#b35806','#5e3c99'))+
      labs(x='Population bridge solution',y='Estimated bridge',colour='Training sample',
        title=paste0(name,', n = 4,000: ',if(ensemble_only)'ensemble bridge' else 'all bridge candidates'),
        subtitle='Outcome at time 3; R2 = 1 and R3 = 1',
        caption='All combinations and all three training samples are shown on equal axes, without clipping.')+
      theme_bw(base_size=11)+theme(legend.position='bottom')
    if(!ensemble_only)figure<-figure+facet_wrap(~candidate_label,ncol=2)
    stem<-paste0(mechanism_name,if(ensemble_only)'-bridge-ensemble' else '-all-bridge-learners')
    for(ext in c('png','pdf'))ggsave(file.path(root,paste0(stem,'.',ext)),figure,
      width=if(ensemble_only)7 else 10,height=if(ensemble_only)6.5 else 9,dpi=180)
  }
  selected<-equations[mechanism==mechanism_name & kind=='adjoint' & horizon==2 & candidate=='ensemble']
  selected[,training_sample:=paste('Training sample',fold)]
  selected[,visit:=factor(current_R,levels=c(0,1),labels=c('Intermediate visit missing','Intermediate visit observed'))]
  limits<-range(selected$true,selected$estimate)
  figure<-ggplot(selected,aes(true,estimate,colour=visit))+
    geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey35')+
    geom_point(alpha=.65,size=1.3)+facet_wrap(~training_sample,nrow=1)+
    coord_equal(xlim=limits,ylim=limits)+
    labs(x='True right-hand side of the adjoint equation',
      y='Conditional mean of the fitted adjoint',colour=NULL,
      title=paste0(name,', n = 4,000: adjoint conditional equations at time 3'),
      caption='All conditioning combinations are shown. Agreement with one chosen adjoint solution is not required.')+
    theme_bw(base_size=11)+theme(legend.position='bottom')
  for(ext in c('png','pdf'))ggsave(file.path(root,paste0(mechanism_name,'-adjoint-equations.',ext)),
    figure,width=12,height=5.5,dpi=180)
}
cat('Saved full-range bridge and defining-equation plots; no estimator-completion claim.\n')
