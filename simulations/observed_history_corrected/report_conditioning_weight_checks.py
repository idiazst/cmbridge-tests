"""Report actual complete checks without assuming favorable outcomes."""
import csv
import math
import os
from pathlib import Path

root = Path(os.environ["SIM_OUTPUT"])
reference = Path(os.environ["SIM_REFERENCE"])
names = {"binary_longitudinal": "Binary treatment", "discrete_dose": "Numerical dose"}
candidates = ["ensemble", "sieve_md", "landweber", "pmmr"]


def read(path):
    return list(csv.DictReader(path.open())) if path.stat().st_size else []


def table(headers, rows):
    return ("| " + " | ".join(headers) + " |\n| " + " | ".join(["---"] * len(headers)) +
            " |\n" + "".join("| " + " | ".join(map(str, row)) + " |\n" for row in rows) + "\n")


def error(rows, mechanism, candidate):
    values = [float(row["rmse"]) ** 2 for row in rows if
              row["mechanism"] == mechanism and row["candidate"] == candidate and
              row["kind"] == "beta" and row["horizon"] == "2" and row["current_R"] == "1"]
    assert len(values) == 3
    return math.sqrt(sum(values) / 3)


estimates, weights, failures, trials, timings = [], [], [], [], []
for mechanism in names:
    folder = next(root.glob(mechanism + "-n4000-seed*"))
    status = (folder / "STATUS.txt").read_text()
    assert "Estimator error rows 0" in status
    timings.append([names[mechanism], float(next(line.split()[1] for line in
                                               status.splitlines() if line.startswith("Seconds ")))])
    current = read(folder / "estimates.csv")
    assert len(current) == 4
    assert all(abs(float(row["mean_sequential_eif"]) + float(row["mean_bridge_adjoint_eif"]) -
                   float(row["estimate"]) + float(row["truth"])) < 1e-10 for row in current)
    estimates += current
    weights += [dict(row, mechanism=mechanism) for row in read(folder / "all-weights.csv")]
    failures += read(folder / "nested-candidate-failures.csv")
    trials += read(folder / "nested-penalty-trial-failures.csv")

metrics = read(root / "function-and-equation-errors.csv")
old = read(reference / "function-and-equation-errors.csv")
solvers = read(root / "solver-tolerance-audit.csv")
solver_failures = [row for row in solvers if row["tolerance_reached"] == "FALSE"]
scores = read(root / "scorer-audit.csv")
clips = read(root / "applied-bridge-truncation.csv")
assert len(scores) == 72 and len(clips) == 12
assert all(int(row["validation_clipped_n"]) == 0 for row in clips)

text = """# Training-only conditioning-weight CV: saved n = 4,000 checks

**This revision is not adopted.** The binary bridge improves, but the numerical-dose bridge remains less accurate than the earlier saved configuration. The separate original cloud study continues unchanged; these diagnostic datasets are never pooled with it.

Both saved datasets finished with four estimates and zero final estimator errors. Seeds are 5103006 and 5203006. The read-only audit verifies unchanged data, outer assignments, learner assignments and EIF component sums. The stricter final-solver audit **fails**: the table below retains every recorded tolerance miss. Successful estimator returns do not imply that this stricter audit passed.

This revision starts from full joint-category conditioning and cell U-statistic selection. It adds training-only CV for the Landweber bridge conditioning-weight ridge, starting with 13 positive values and extending the grid until an interior minimum or a recorded failure. Sieve penalties retain their 41-value initial positive grid and extension rules. All raw selection scores exclude self-products; only the ensemble Gram is projected to positive semidefinite. No generating function enters fitting or selection.

Packages are cmbridge 0.3.0.9021 and modified lmtp 1.6.0.9021. Target classes, predictors, inverse-expit bridges, unrestricted adjoints, PMMR, regression/treatment libraries, strict nested response rebuilding, shared training-only splits, SDR, logistic TMLE and one-step bridge correction are unchanged. No old nuisance caches were copied. Saturated L1 is absent from cmbridge and remains in the SuperLearner regression/classification libraries with CV; MARS appears only in regression.

"""
text += table(["Example", "Measured seconds", "Measured minutes"],
              [[name, f"{seconds:.3f}", f"{seconds/60:.2f}"] for name, seconds in timings])
text += "Both checks ran concurrently on the Mac; these are elapsed times, not a cloud completion forecast.\n\n"
text += "## True versus estimated bridge at time 3\n\n"
text += table(["Example", "Candidate", "Earlier RMSE", "Conditioning-weight CV RMSE"],
              [[names[m], c, f"{error(old,m,c):.6f}", f"{error(metrics,m,c):.6f}"]
               for m in names for c in candidates])
text += "RMSE uses population probabilities within the group with both follow-up visits observed, then averages squared errors across the three training fits. It measures function error, not interval coverage. All combinations appear on equal axes without clipping.\n\n"
for mechanism in names:
    text += f"![{names[mechanism]} ensemble]({mechanism}-bridge-ensemble.png)\n\n"
    text += f"![{names[mechanism]} candidates]({mechanism}-all-bridge-learners.png)\n\n"

for kind, label in [("beta", "Bridge"), ("adjoint", "Adjoint")]:
    rows = []
    for mechanism in names:
        for fold in range(1, 4):
            selected = {row["candidate"]: float(row["weight"]) for row in weights if
                        row["mechanism"] == mechanism and row["kind"] == kind and
                        row["horizon_or_depth"] == "2" and row["fold"] == str(fold)}
            assert set(selected) == set(candidates[1:])
            rows.append([names[mechanism], fold] + [f"{selected[c]:.6f}" for c in candidates[1:]])
    text += f"## {label} ensemble weights at time 3\n\n" + table(
        ["Example", "Training fit", "Sieve", "Landweber", "PMMR"], rows)

text += "## Estimates and sample mean EIF components\n\n" + table(
    ["Example", "Outcome time", "Estimator", "Truth", "Estimate", "SE", "Sequential EIF mean", "Bridge/adjoint EIF mean", "Pointwise CI contains truth"],
    [[names[row["mechanism"]], int(row["horizon"]) + 1, row["estimator"].upper()] +
     [f"{float(row[column]):.6f}" for column in
      ["truth", "estimate", "se", "mean_sequential_eif", "mean_bridge_adjoint_eif"]] +
     ["Yes" if float(row["lower"]) <= float(row["truth"]) <= float(row["upper"]) else "No"]
     for row in estimates])
misses = [row for row in estimates if not float(row["lower"]) <= float(row["truth"]) <= float(row["upper"])]
text += f"{len(misses)} of eight pointwise intervals exclude the truth. One dataset per mechanism cannot estimate coverage or establish a convergence rate. These are actual sample EIF means; their sum equals estimate minus truth.\n\n"
text += "## Recorded final solver tolerance misses\n\n" + table(
    ["Example", "Training fit", "Outcome time", "Function", "Candidate", "Recorded coefficient gradient"],
    [[names[row["mechanism"]], row["fold"], int(row["horizon"]) + 1, row["kind"], row["candidate"], row["gradient"]]
     for row in solver_failures])
text += "No tolerance was relaxed. Separate same-positive-ridge iteration checks reached the original tolerance without removing the poor dose fit; original models and the failed audit are preserved. Increasing an iteration ceiling is distinct from choosing a penalty using truth.\n\n"
text += "## Defining equations and retained failures\n\n"
for mechanism in names:
    text += f"![{names[mechanism]} adjoint equations]({mechanism}-adjoint-equations.png)\n\n"
text += f"Independent ordered-pair cell formulas reproduce all 24 raw ensemble Grams and verify PSD projection. Selected penalties are positive, finite, successful exact interior minima with shared training-only person IDs. All final/learner-training CV paths are checked; the solver gate still fails as stated above. There is no clipping in cmbridge and no changes from the approved [-100,100] application interval in these outer checks. The complete estimators retain {len(failures)} nested candidate-failure records and {len(trials)} failed penalty-trial records. These are records, not independent statistical events.\n\n"
text += "Full-support bridge class containment does not guarantee accurate fitted coefficients. Adjoint diagnostics check the defining conditional equation, allowing multiple solutions; the sample-dependent sieve dictionary can omit combinations absent from training. The observed dose scatter is still under investigation. The exact-kernel diagnostic is preserved separately and cannot establish estimator performance on its own.\n\n"
text += "[All function/equation errors](function-and-equation-errors.csv) · [Score audit](scorer-audit.csv) · [Positive sieve penalties](sieve-selected-penalties.csv) · [Solver flags](solver-tolerance-audit.csv) · [Applied bounds](applied-bridge-truncation.csv).\n"
(root / "REPORT.md").write_text(text)
(root / "README.md").write_text("# Saved n = 4,000 conditioning-weight CV checks\n\nBoth complete-estimator checks have finished. This configuration is **not adopted**: binary bridge error improves, numerical-dose bridge error remains worse, and three final Landweber fits miss their recorded tolerance. All actual estimates, sample EIF means, weights, raw plots and failures are in [the report](REPORT.md). The separate original cloud study remains unchanged.\n")
print(f"Reported eight completed estimates, {len(misses)} CI misses, and {len(solver_failures)} final tolerance misses.")
