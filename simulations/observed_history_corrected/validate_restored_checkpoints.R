# Read-only validation before the unchanged dispatcher resumes unfinished fits.
library(data.table)
outdir <- Sys.getenv('SIM_OUTPUT')
shard <- as.integer(Sys.getenv('SIM_SHARD'))
design <- fread(file.path(outdir,'design.csv'))
fingerprint <- unname(readLines(file.path(outdir,'source-fingerprint.txt'))[1L])
stopifnot(nrow(design)==1200L,uniqueN(design$job)==1200L,shard%in%1:40)
paths <- list.files(file.path(outdir,'jobs'),'\\.rds$',full.names=TRUE)
for(path in paths) {
  wanted_job <- sub('\\.rds$','',basename(path))
  config <- design[job==wanted_job]
  x <- readRDS(path)
  stopifnot(nrow(config)==1L,config$shard==shard,
    identical(x$fingerprint,fingerprint),x$design$seed==config$seed,
    x$design$n==config$n,x$design$mechanism==config$mechanism,
    x$design$replicate==config$replicate)
  if(is.null(x$errors)||!nrow(x$errors))stopifnot(!is.null(x$results),nrow(x$results)==4L,
    setequal(x$results$estimator,c('sdr','tmle')),setequal(x$results$horizon,1:2))
  if(!is.null(x$results)&&nrow(x$results))stopifnot(!anyDuplicated(x$results[,.(estimator,horizon)]))
}
cat(length(paths),'restored checkpoints validated; saved failures retained.\n')
