# Reproduce only data and sample assignments, never fit known nuisance functions.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
proof <- rbindlist(lapply(c('binary_longitudinal','discrete_dose'), function(mechanism) {
  seed <- 5000000L + match(mechanism,c('binary_longitudinal','discrete_dose')) * 100000L + 3000L + 1L
  set.seed(seed)
  g <- make_mechanism(mechanism,dose_max=3L,visits='missing_both')
  require_valid_mechanism(g)
  stopifnot(max(abs(exact_truth(g)-backward_truth(g)))<1e-12)
  d <- draw_data(4000L,g)
  task <- study_task(d,g,folds=3L,learner_groups=3L,balance_measurement=TRUE)$task
  hash <- function(x) digest::digest(x,serializeVersion=2L)
  data.table(mechanism=mechanism,n=4000L,seed=seed,
    data_hash=hash(as.data.frame(d)),outer_hash=hash(task$folds),learner_hash=hash(task$learner_folds))
}))
path <- Sys.getenv('SPLIT_PROOF_OUTPUT')
stopifnot(nzchar(path))
fwrite(proof,path)
expected <- Sys.getenv('SPLIT_PROOF_EXPECTED')
if(nzchar(expected)) stopifnot(identical(as.data.frame(proof),as.data.frame(fread(expected))))
print(proof)
