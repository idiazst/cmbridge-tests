"""Summarize the completed binary million-person adjoint diagnostic."""
from pathlib import Path
import csv
import hashlib
import shutil

repo = Path(__file__).resolve().parents[2]
root = repo / "results/observed-history-corrected/missing-both-study-v2"
out = root / "diagnostics/binary-dose-n1000000-seed5103006-adjoint"
if not (out / "large-adjoint-fits.rds").exists():
    raise SystemExit("The adjoint fits have not finished; do not report completion.")
with (out / "fit-diagnostics.csv").open() as stream:
    metrics = list(csv.DictReader(stream))
assert len(metrics) == 6
assert all(int(row["represented_combinations"]) == int(row["target_combinations"]) for row in metrics)
with (out / "population-uniqueness.csv").open() as stream:
    uniqueness = list(csv.DictReader(stream))
assert len(uniqueness) == 40 and all(row["unique_solution"] == "TRUE" for row in uniqueness)
with (out / "adjoint-predictions.csv").open() as stream:
    points = list(csv.DictReader(stream))
assert len(points) == 408
observed_metrics = [row for row in metrics if row["R2"] == "1"]
missing_metrics = [row for row in metrics if row["R2"] == "0"]
observed_points = [row for row in points if row["R2"] == "1"]
worst = max(observed_points, key=lambda row: abs(float(row["estimate"]) - float(row["true"])))
largest_error = abs(float(worst["estimate"]) - float(worst["true"]))
lines = [
    "# Binary-treatment saturated sieve adjoint at n = 1,000,000", "",
    "The sample is exactly the preceding bridge diagnostic: n = 1,000,000, seed 5103006, missingness at both follow-ups. Regeneration reproduces the saved post-draw random-number-generator state, and all three saved outer splits are reused. The full observed history and original frozen package versions are retained: cmbridge 0.3.0.9005 and lmtp 1.6.0.9005.", "",
    "Only the saturated joint-category sieve_md adjoint for the outcome at time 3 is fitted. Its target and conditioning bases are unchanged, both ridge values are 1e−8, and its predictions are unrestricted. The adjoint function can have negative values. No other cmbridge candidate, sequential regression, final SDR/TMLE estimator, or confidence interval is fitted.", "",
    "The adjoint training response uses the product of the two estimated treatment ratios times the measured outcome. For every outer training sample, three inner splits supply out-of-sample treatment-ratio predictions. Each ratio fit uses the original SuperLearner library, SL.saturated_l1_cv and SL.mean, with shared person assignments and positive L1 scales 0.001, 0.01, and 0.1 chosen entirely within that ratio fit's training data. Both stacked copies of a person stay in the same group. No outer validation observation participates in these fits. True nuisance functions are used only for plots and diagnostics.", "",
    "## The adjoint solution is unique here", "",
    "An exact population check covers all 40 reachable time-2 histories. When R₂ = 1, four reachable treatment combinations have an adjoint equation matrix of rank four. When R₂ = 0, the one deterministic treatment combination has rank one. Thus the adjoint is unique on the reachable support in both cases, and direct comparison to its true value is meaningful. This differs from the nonunique bridge at R₂ = 0 and from adjoints in some numerical-dose settings.", "",
    "All 136 adjoint-predictor combinations are represented in measured training observations in every fit, so the fitted saturated classes contain this unique true adjoint.", "",
    "![True versus estimated adjoint](adjoint-true-versus-estimated.png)", "",
    "| Training fit | R₂ | Probability-weighted function error | Largest function error | Probability-weighted equation error | Largest equation error |",
    "| --- | ---: | ---: | ---: | ---: | ---: |",
]
for row in metrics:
    lines.append(
        f"| {row['fold']} | {row['R2']} | {float(row['adjoint_rmse']):.6f} | "
        f"{float(row['adjoint_max_absolute_error']):.6f} | {float(row['equation_rmse']):.6f} | "
        f"{float(row['equation_max_absolute_error']):.6f} |"
    )
lines += ["",
    "Function error is the square root of mean squared difference from the true adjoint, weighted by the population history/treatment probabilities conditional on the indicated R₂ value. Equation error compares the conditional mean of the fitted adjoint given the full bridge predictors and measured outcome with the true right-hand side of the adjoint equation. It is weighted by measured-outcome population probabilities conditional on the indicated R₂ value. Their weights and units differ from the bridge-function error, so these numbers should not be compared as identical accuracy measures.", "",
    "![Adjoint equation](adjoint-equation-true-versus-estimated.png)", "",
    f"Most probability is near y = x, but rare large errors remain. The largest function error is {largest_error:.6f} in fit {worst['fold']}: true value {float(worst['true']):.6f}, estimated value {float(worst['estimate']):.6f}. In fit 1 the largest departures occur in the same history as the preceding bridge's worst point, with 292 training observations and 129 measured final outcomes. Counts are relevant to interpreting recovery; this check does not establish that they explain every error.", "",
    "Combining each new adjoint with its saved sieve bridge reproduces the paper's product-of-errors identity to within 1e−8. The separately saved bridge remainder diagnostic uses true treatment ratios to evaluate this identity, not to fit an estimator. It is neither a reported mean EIF component nor a final estimator's bias.", "",
    "The estimated-ratio-target column in the diagnostics is a separate population evaluation of the mixture of the three fitted ratio functions. It helps assess loading error but does not reproduce the exact finite-sample mean of the training responses or isolate the entire cause of adjoint estimation error.", "",
    "The full simulation remains stopped at 118 completed datasets. This diagnostic is stored outside its job checkpoints.", "",
    "[Function PDF](adjoint-true-versus-estimated.pdf), [equation PDF](adjoint-equation-true-versus-estimated.pdf), [predictions](adjoint-predictions.csv), [equation values](adjoint-equation-predictions.csv), [fit diagnostics](fit-diagnostics.csv), [population uniqueness](population-uniqueness.csv), [ratio diagnostics](ratio-fit-diagnostics.csv), [ratio weights](ratio-learner-weights.csv), and [ratio penalty selection](ratio-penalty-selection.csv).",
]
(out / "README.md").write_text("\n".join(lines) + "\n")
shutil.copy2(repo / "simulations/observed_history_corrected/check_large_missing_both_adjoint.R", out)
shutil.copy2(Path(__file__), out)
shutil.copy2(root / "diagnostics/binary-large-adjoint-run.log", out / "run.log")
report = root / "REPORT.md"
text = report.read_text()
section = "\n\n## Binary saturated sieve adjoint at n = 1,000,000\n\n"
section += "The adjoint diagnostic uses exactly the large binary bridge dataset (seed 5103006) and its three outer splits. Only the saturated sieve adjoint is fitted, using training-out-of-sample treatment ratios from the original SuperLearner library. All 136 adjoint predictor combinations are represented, and an exact population audit confirms a unique adjoint on the reachable support for both R₂ values.\n\n"
section += "![Large binary sieve adjoint](diagnostics/binary-dose-n1000000-seed5103006-adjoint/adjoint-true-versus-estimated.png)\n\n"
section += "Probability-weighted function errors where R₂ = 1 are " + ", ".join(f"{float(row['adjoint_rmse']):.4f}" for row in observed_metrics) + ". Where R₂ = 0 they are " + ", ".join(f"{float(row['adjoint_rmse']):.4f}" for row in missing_metrics) + f". Most probability is near y = x, with rare errors up to {largest_error:.2f}. [The diagnostic](diagnostics/binary-dose-n1000000-seed5103006-adjoint/README.md) also reports the adjoint equations, treatment-ratio fits, and product identity. The full simulation remains stopped.\n"
if "## Binary saturated sieve adjoint at n = 1,000,000" not in text:
    report.write_text(text + section)
for folder in (out,):
    checksums = []
    for file in sorted(folder.iterdir()):
        if file.is_file() and file.name != "SHA256SUMS":
            checksums.append(f"{hashlib.sha256(file.read_bytes()).hexdigest()}  {file.name}\n")
    (folder / "SHA256SUMS").write_text("".join(checksums))
manifest = root / "REPORT-SHA256SUMS"
checksums = []
for line in manifest.read_text().splitlines():
    _, name = line.split("  ", 1)
    checksums.append(f"{hashlib.sha256((root / name).read_bytes()).hexdigest()}  {name}\n")
manifest.write_text("".join(checksums))
print("Adjoint report and checksums saved. Full-study dataset count:", len(list((root / "jobs").glob("*.rds"))))
