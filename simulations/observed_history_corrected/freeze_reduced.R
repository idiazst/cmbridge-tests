out<-normalizePath('results/observed-history-corrected',mustWork=TRUE)
out<-file.path(out,'reduced-study')
if(dir.exists(file.path(out,'jobs'))&&length(list.files(file.path(out,'jobs'),'\\.rds$')))
 stop('This run already has results. Resume its frozen sources; do not replace them.')
for(name in c('packages','source','library','validation'))dir.create(file.path(out,name),recursive=TRUE,showWarnings=FALSE)
old<-normalizePath('results/observed-history-corrected/study',mustWork=TRUE)
archives<-file.path(old,'packages',c('cmbridge_0.3.0.9002.tar.gz','lmtp_1.6.0.9004.tar.gz'))
destinations<-file.path(out,'packages',basename(archives))
stopifnot(all(file.copy(archives,destinations,overwrite=TRUE)))
for(archive in destinations)untar(archive,exdir=file.path(out,'source'))
scripts<-file.path(out,'source','simulations','observed_history_corrected')
dir.create(scripts,recursive=TRUE,showWarnings=FALSE)
stopifnot(all(file.copy(list.files('simulations/observed_history_corrected',full.names=TRUE),scripts,overwrite=TRUE)))
reports<-file.path(out,'source','reports');dir.create(reports,showWarnings=FALSE)
file.copy('reports/observed-history-simulation.md',file.path(reports,'simulation-plan.md'),overwrite=TRUE)
stopifnot(all(file.copy(list.files(file.path(old,'library'),full.names=TRUE),file.path(out,'library'),recursive=TRUE)))
file.copy(list.files(file.path(old,'validation'),full.names=TRUE),file.path(out,'validation'),overwrite=TRUE)
file.copy(list.files('results/observed-history-corrected/revised-design',full.names=TRUE),file.path(out,'validation'),overwrite=TRUE)
for(name in c('cmbridge','lmtp'))stopifnot(any(grepl('Status: OK',readLines(file.path(out,'validation',paste0(name,'-check.log'))))))
files<-c(destinations,list.files(scripts,full.names=TRUE),file.path(reports,'simulation-plan.md'))
hashes<-vapply(files,function(path)digest::digest(file=path,algo='sha256'),'')
writeLines(paste(hashes,substring(files,nchar(out)+2L),sep='  '),file.path(out,'SHA256SUMS'))
quote<-shQuote
launcher<-c('#!/bin/sh','set -eu',paste('cd',quote(normalizePath('.'))),
 paste0('export R_LIBS=',quote(file.path(out,'library'))),
 paste0('export STUDY_SOURCE=',quote(scripts)),
 paste0('export LMTP_SOURCE=',quote(file.path(out,'source','lmtp'))),
 paste0('export CMBRIDGE_SOURCE=',quote(file.path(out,'source','cmbridge'))),
 paste0('export SIM_OUTPUT=',quote(out)),
 'export SIM_WORKERS=8','export SIM_POSTPROCESS=1',
 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1',
 'unset SIM_REPS SIM_SIZES',
 paste(quote(file.path(R.home('bin'),'Rscript')),quote(file.path(scripts,'run_reduced.R'))))
writeLines(launcher,file.path(out,'run-study.sh'));Sys.chmod(file.path(out,'run-study.sh'),'0755')
cat('Revised study frozen at',out,'\n')
