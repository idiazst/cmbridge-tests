study_source <- Sys.getenv("STUDY_SOURCE", "simulations/observed_history_corrected")
source(file.path(study_source, "study.R"))
outdir <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/study")
dir.create(file.path(outdir, "jobs"), recursive = TRUE, showWarnings = FALSE)
workers <- as.integer(Sys.getenv("SIM_WORKERS", "4"))
replications <- as.integer(Sys.getenv("SIM_REPS", "500"))
sizes <- as.integer(strsplit(Sys.getenv("SIM_SIZES", "1000,2000,4000,8000,16000,32000,64000"), ",")[[1]])
mechanisms <- strsplit(Sys.getenv("SIM_MECHANISMS", "binary_point,binary_longitudinal,discrete_dose"), ",")[[1]]
grid <- CJ(mechanism = mechanisms, n = sizes, replicate = seq_len(replications))
# Seeds are invariant to worker count, pilot grids, batching, or resumption.
grid[, seed := 2000000L + match(mechanism, c("binary_point", "binary_longitudinal", "discrete_dose")) * 100000L +
       match(n, c(1000L, 2000L, 4000L, 8000L, 16000L, 32000L, 64000L)) * 1000L + replicate]
grid[, job := sprintf("%s-n%d-r%03d", mechanism, n, replicate)]
fwrite(grid, file.path(outdir, "design.csv"))
lmtp_source <- Sys.getenv("LMTP_SOURCE", "../lmtp")
cmbridge_source <- Sys.getenv("CMBRIDGE_SOURCE", "../cmbridge")
sources <- c(file.path(study_source, c("dgp.R", "study.R", "run.R")),
  list.files(file.path(lmtp_source, "R"), "\\.R$", full.names = TRUE),
  list.files(file.path(cmbridge_source, "R"), "\\.R$", full.names = TRUE),
  file.path(lmtp_source, "DESCRIPTION"), file.path(cmbridge_source, "DESCRIPTION"))
fingerprint <- digest::digest(c(tools::md5sum(sort(sources)), packageVersion("lmtp"), packageVersion("cmbridge")))
manifest <- file.path(outdir, "source-fingerprint.txt")
if (file.exists(manifest) && readLines(manifest)[1] != fingerprint)
  stop("Source differs from this run. Archive the run and rerun the original seeds; do not mix versions.")
writeLines(fingerprint, manifest)
writeLines(c(capture.output(sessionInfo()), paste("workers", workers), paste("fingerprint", fingerprint)),
           file.path(outdir, "session.txt"))
selected <- order(grid$replicate, grid$n, grid$mechanism)
if (nzchar(Sys.getenv("SIM_JOBS"))) selected <- as.integer(strsplit(Sys.getenv("SIM_JOBS"), ",")[[1]])
one <- function(i) {
  config <- grid[i]; path <- file.path(outdir, "jobs", paste0(config$job, ".rds"))
  if (file.exists(path)) {
    previous <- readRDS(path)
    stopifnot(identical(previous$fingerprint, fingerprint), identical(previous$design$seed, config$seed))
    return(invisible(TRUE))
  }
  started <- proc.time()[3]
  initial_warnings <- character()
  result <- tryCatch(withCallingHandlers(run_dataset(config$mechanism, config$n, config$replicate, config$seed),
    warning = function(w) {initial_warnings <<- c(initial_warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) list(design = as.list(config[, .(mechanism, n, replicate, seed)]),
      errors = data.table(error = conditionMessage(e))))
  result$fingerprint <- fingerprint; result$seconds <- proc.time()[3] - started
  result$job_warnings <- unique(initial_warnings)
  tmp <- paste0(path, ".tmp")
  saveRDS(result, tmp, compress = "gzip"); stopifnot(file.rename(tmp, path))
  status <- if (nrow(result$errors)) paste(unique(result$errors$error), collapse = " | ") else "OK"
  cat(config$job, sprintf("%.1fs", result$seconds), status, "\n")
  invisible(TRUE)
}
batches <- split(selected, grid$replicate[selected])
for (batch in batches) {
  answers <- parallel::mclapply(batch, one, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  saved <- file.exists(file.path(outdir, "jobs", paste0(grid$job[batch], ".rds")))
  if (length(answers) != length(batch) || !all(vapply(answers, isTRUE, FALSE)) || !all(saved))
    stop("A worker failed to save its checkpoint; inspect and resume the unfinished jobs.")
  replicate <- grid$replicate[batch[1L]]
  if (Sys.getenv("SIM_POSTPROCESS", "0") == "1" && replicate %in% c(1L, 5L, 10L, 25L, 50L, 100L, 250L, 500L)) {
    for (script in c("summarize.R", "write_report.R")) {
      status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(study_source, script)))
      if (status != 0) stop("Postprocessing failed; saved simulation checkpoints remain available.")
    }
  }
}
cat("Requested jobs completed.\n")
