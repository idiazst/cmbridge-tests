library(data.table)
outdir <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/study")
s <- fread(file.path(outdir, "summary.csv"))
design <- fread(file.path(outdir, "design.csv")); log <- fread(file.path(outdir, "job_log.csv"))
full_design <- nrow(design) == 10500L && all(design[, .N, by = .(mechanism, n)]$N == 500L)
complete <- full_design && nrow(log) == nrow(design)
lines <- c("# Observed-history simulation results", "", if (complete)
  "The planned datasets have finished. Numerical failures are retained and counted below." else
  if (full_design) "**The study is running. These are preliminary outputs; they do not establish bias or coverage.**" else
  "**Implementation check only. These are not results of the planned 500-replication study.**", "",
  sprintf("Completed %s of %s planned datasets. The study uses 500 replications for each of three mechanisms and seven sample sizes, with both estimators and all 16 specifications evaluated on each dataset.",
    format(nrow(log), big.mark = ","), format(nrow(design), big.mark = ",")), "",
  "The data-generating distributions and estimator specifications are in [the approved setup](../../../reports/observed-history-simulation.md). All penalties are selected by cross-validation. No known-function or zero-penalty simulation is included.", "",
  "## All four function classes correctly specified", "",
  "Tables distinguish successful fits, failed fits, and datasets still pending. Intervals for coverage describe simulation uncertainty. Correct specification describes the function class; the population diagnostics quantify the error of the fitted functions.", "")
table <- function(x) {
  if (!nrow(x)) return(character())
  x <- as.data.frame(x)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")), "")
}
primary <- s[horizon == ifelse(mechanism == "binary_point", 1, 2)]
correct <- primary[beta_correct & lambda_correct & ratio_correct & m_correct]
for (name in c("binary_point", "binary_longitudinal", "discrete_dose")) {
  x <- correct[mechanism == name][order(n, estimator)]
  lines <- c(lines, paste0("### ", gsub("_", " ", name)), "", table(x[, .(
    n, Estimator = toupper(estimator), Successful = successful, Failed = failed, Pending = pending,
    Truth = sprintf("%.5f", truth), Mean = sprintf("%.5f", mean_estimate),
    Bias = sprintf("%.5f", bias), SD = sprintf("%.5f", empirical_sd),
    RMSE = sprintf("%.5f", rmse), `Mean SE` = sprintf("%.5f", mean_se),
    `Interval length` = sprintf("%.5f", mean_interval_length), Coverage = sprintf("%.3f", coverage),
    `Coverage uncertainty` = sprintf("[%.3f, %.3f]", coverage_mc_lower, coverage_mc_upper))]))
}
lines <- c(lines, "## All 16 specifications", "",
  "C means a correctly specified function class and M means the stated misspecified model. The columns identify β, λ, ω, and m individually. The nine sufficient consistency conditions in the setup do not by themselves guarantee valid confidence intervals under persistent misspecification.", "")
for (name in c("binary_point", "binary_longitudinal", "discrete_dose")) for (method in c("sdr", "tmle")) {
  x <- primary[mechanism == name & estimator == method]
  x[, specification := paste(ifelse(beta_correct, "C", "M"), ifelse(lambda_correct, "C", "M"),
    ifelse(ratio_correct, "C", "M"), ifelse(m_correct, "C", "M"), sep = "/")]
  x[, entry := ifelse(pending > 0, sprintf("%.4f; %.3f (%d)", bias, coverage, successful),
    sprintf("%.4f; %.3f", bias, coverage))]
  wide <- dcast(x, specification ~ n, value.var = "entry")
  setnames(wide, "specification", "β / λ / ω / m")
  lines <- c(lines, paste0("### ", gsub("_", " ", name), ": ", toupper(method)), "",
    "Each entry gives mean estimation error; coverage. While the study runs, parentheses give completed successful replications.", "", table(wide))
}
lines <- c(lines, "## Graphs", "",
  "Estimation error and its √n multiple show whether the remaining error decreases sufficiently quickly. Coverage is compared with 0.95; standard-error plots compare estimated uncertainty with variation across replications.", "")
for (name in c("binary_point", "binary_longitudinal", "discrete_dose")) for (kind in c("bias", "root-n-bias", "coverage", "standard-errors"))
  lines <- c(lines, sprintf("![%s: %s](figures/%s-%s.png)", name, kind, name, kind), "")
lines <- c(lines, "![Estimated curves compared with truth](figures/correctly-specified-curves.png)", "",
  "## Results requiring examination", "")
largest <- max(design$n)
flag <- correct[n == largest & successful >= 100L &
  (coverage_mc_upper < .95 | abs(bias) > 1.96 * bias_mcse)]
if (!nrow(flag)) lines <- c(lines, if (complete)
  "The automated check found no largest-sample, fully specified result with coverage clearly below 0.95 or mean error clearly different from zero. This screening is not a proof of the paper's asymptotic conditions." else
  "The run is incomplete. Coverage and convergence conclusions are deferred until the planned replications finish.", "") else {
  lines <- c(lines, "The following fully specified results require examination. A flag is not proof of an implementation error. The separate population remainder contributions, conditional equations, selected penalties, and training counts must be considered before attributing the result to a package defect.", "",
    table(flag[, .(mechanism, estimator, n, successful, bias, coverage, coverage_mc_upper)]))
}
lines <- c(lines, "## Reproduction and checks", "",
  "Every successful comparison verifies the public estimator's prediction-based one-step calculation and the exact population identity for the bridge contribution. All conditional-moment and regression fits evaluated on estimator validation observations exclude those observations. Candidate penalty selection excludes the candidate's validation observations.", "",
  "[Full numerical summary](summary.csv), [replication results](replicates.csv), [population diagnostics](population_diagnostics.csv), [penalties and weights](penalty_and_weight_selections.csv), and [job log](job_log.csv) retain the information behind the tables. Checkpoints additionally retain the fitted functions on the entire finite support, predictor-combination counts, and sample assignments.", "",
  "The adjoint equation may have multiple solutions. Difference from the particular exact solution used in one diagnostic is not by itself treated as misspecification; conditional-equation errors and the actual remainder are also recorded.", "")
writeLines(lines, file.path(outdir, "REPORT.md"))
if (requireNamespace("markdown", quietly = TRUE))
  markdown::mark_html(file.path(outdir, "REPORT.md"), output = file.path(outdir, "REPORT.html"))
cat("Report written to", file.path(outdir, "REPORT.md"), "\n")
