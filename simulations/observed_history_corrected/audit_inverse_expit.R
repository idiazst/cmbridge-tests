library(data.table)
library(cmbridge)
root <- Sys.getenv('SIM_OUTPUT')
reference <- Sys.getenv('SIM_REFERENCE')
rows <- list()
for (mechanism in c('binary_longitudinal','discrete_dose')) for (sample_n in c(4000L,20000L)) {
  seed <- if (mechanism == 'binary_longitudinal') 5103006L else 5203006L
  folder <- file.path(root, paste0(mechanism,'-n',sample_n,'-seed',seed))
  for (fold in 1:3) {
    fit <- readRDS(file.path(folder,paste0('nuisance-fits-fold',fold,'.rds')))
    for (s in 1:2) for (kind in c('beta','adjoint')) {
      model <- fit[[paste0(kind,'_fits')]][[s]]
      stopifnot(identical(names(model$weights),c('sieve_md','landweber','pmmr')),
                abs(sum(model$weights)-1)<1e-10)
      sieve <- model$candidates$sieve_md
      if (kind == 'beta') {
        solver <- sieve$solver
        stopifnot(sieve$tuning$link == 'inverse_logit', all(is.finite(sieve$link_coefficients)),
          solver$converged, max(solver$coefficient_gradient,solver$function_gradient)<=solver$tolerance,
          all(fitted(sieve)[!is.na(fitted(sieve))]>=1))
        values <- 1+exp(-sieve$link_coefficients)
        stopifnot(all(is.finite(values)),all(values>=1))
        penalty <- sieve$tuning$lambda * sum(sieve$coefficients[-1L]^2)
        stopifnot(abs(solver$penalized_loss - (sieve$moment_loss + penalty)) < 1e-7)
        rows[[length(rows)+1L]] <- data.table(mechanism,n=sample_n,fold,horizon=s,
          minimum_value=min(values),maximum_value=max(values),cells=length(values),
          near_one_cells=sum(values-1<1e-6),minimum_coefficient=min(sieve$link_coefficients),
          maximum_coefficient=max(sieve$link_coefficients),
          attempts=nrow(solver$trace),restarts=sum(grepl('restart',solver$trace$method)),
          optimizer_stop_code=solver$optimizer$convergence,
          coefficient_gradient=solver$coefficient_gradient,function_gradient=solver$function_gradient,
          accepted_threshold=solver$tolerance,penalized_loss=solver$penalized_loss)
      } else stopifnot(sieve$tuning$link == 'identity',is.null(sieve$solver))
    }
  }
}
fwrite(rbindlist(rows),file.path(root,'inverse-expit-optimizer-audit.csv'))
new <- fread(file.path(root,'all-function-predictions.csv'))
old <- fread(file.path(reference,'all-function-predictions.csv'))
keys <- c('mechanism','n','kind','horizon','current_R','fold','candidate','cell')
# Sieve bridge and bridge ensemble change. All other functions should agree.
a <- new[!(kind=='beta' & candidate %in% c('sieve_md','ensemble'))]
b <- old[!(kind=='beta' & candidate %in% c('sieve_md','ensemble'))]
comparison <- merge(a[,c(keys,'estimate'),with=FALSE],b[,c(keys,'estimate'),with=FALSE],by=keys,suffixes=c('_new','_old'))
stopifnot(nrow(comparison)==nrow(a),max(abs(comparison$estimate_new-comparison$estimate_old))<1e-10)
positive <- new[kind=='beta' & candidate=='sieve_md']
stopifnot(all(is.finite(positive$estimate)),all(positive$estimate>=1))
writeLines(c(paste('Compared unchanged function rows:',nrow(comparison)),
 paste('Maximum difference:',max(abs(comparison$estimate_new-comparison$estimate_old))),
 paste('Minimum inverse-expit sieve prediction:',min(positive$estimate))),file.path(root,'inverse-expit-audit.txt'))
cat('Verified unrestricted finite coefficients, inverse-expit predictions, objective consistency, and unchanged other learners.\n')
