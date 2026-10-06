# Prepare immutable source/package copies before starting the long run.
outdir <- normalizePath("results/observed-history-corrected", mustWork = TRUE)
outdir <- file.path(outdir, "study")
if (dir.exists(file.path(outdir, "jobs")) &&
    length(list.files(file.path(outdir, "jobs"), "\\.rds$")))
  stop("This study already has checkpoints. Resume its frozen run; do not replace its sources.")
for (name in c("packages", "source", "library", "validation"))
  dir.create(file.path(outdir, name), recursive = TRUE, showWarnings = FALSE)
archives <- c("../cmbridge/cmbridge_0.3.0.9002.tar.gz", "../lmtp/lmtp_1.6.0.9004.tar.gz")
stopifnot(all(file.exists(archives)))
destinations <- file.path(outdir, "packages", basename(archives))
stopifnot(all(file.copy(archives, destinations, overwrite = TRUE)))
for (archive in destinations) untar(archive, exdir = file.path(outdir, "source"))
scripts <- file.path(outdir, "source", "simulations", "observed_history_corrected")
dir.create(scripts, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.copy(list.files("simulations/observed_history_corrected", full.names = TRUE),
  scripts, overwrite = TRUE)))
# OrdMonReg is an existing imported dependency, although these curves do not
# use isotonic projection. Keep that dependency with the frozen installations.
ord <- normalizePath("../lmtp/.library/OrdMonReg", mustWork = TRUE)
if (!dir.exists(file.path(outdir, "library", "OrdMonReg")))
  stopifnot(file.copy(ord, file.path(outdir, "library"), recursive = TRUE))
Sys.setenv(R_LIBS = file.path(outdir, "library"))
for (archive in destinations) {
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
    shQuote(paste0("--library=", file.path(outdir, "library"))), shQuote(archive)))
  if (status != 0) stop("Frozen installation failed.")
}
for (name in c("cmbridge", "lmtp")) {
  check <- file.path("..", name, paste0(name, ".Rcheck"), "00check.log")
  stopifnot(file.exists(check), any(grepl("Status: OK", readLines(check))))
  file.copy(check, file.path(outdir, "validation", paste0(name, "-check.log")), overwrite = TRUE)
}
files <- c(destinations, list.files(scripts, full.names = TRUE))
hashes <- vapply(files, function(path) digest::digest(file = path, algo = "sha256"), "")
writeLines(paste(hashes, sub(paste0(outdir, "/"), "", files, fixed = TRUE), sep = "  "),
  file.path(outdir, "SHA256SUMS"))
cwd <- normalizePath(".")
quote <- shQuote
launcher <- c("#!/bin/sh", "set -eu", paste("cd", quote(cwd)),
  paste0("export R_LIBS=", quote(file.path(outdir, "library"))),
  paste0("export STUDY_SOURCE=", quote(scripts)),
  paste0("export LMTP_SOURCE=", quote(file.path(outdir, "source", "lmtp"))),
  paste0("export CMBRIDGE_SOURCE=", quote(file.path(outdir, "source", "cmbridge"))),
  paste0("export SIM_OUTPUT=", quote(outdir)),
  "export SIM_WORKERS=4", "export SIM_POSTPROCESS=1",
  "unset SIM_REPS SIM_SIZES SIM_MECHANISMS SIM_JOBS",
  paste(quote(file.path(R.home("bin"), "Rscript")), quote(file.path(scripts, "run.R"))))
writeLines(launcher, file.path(outdir, "run-study.sh"))
Sys.chmod(file.path(outdir, "run-study.sh"), "0755")
cat("Frozen study prepared at", outdir, "\n")
