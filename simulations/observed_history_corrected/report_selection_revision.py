"""Write the report from the four saved selection-revision diagnostics."""
from pathlib import Path
import csv
import os
from datetime import datetime

root = Path(os.environ.get("SIM_OUTPUT", "results/observed-history-corrected/selection-revision"))


def read(path):
    if not path.exists() or not path.stat().st_size:
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream))


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


labels = {"binary_longitudinal": "Binary treatment", "discrete_dose": "Numerical dose"}
seeds = {"binary_longitudinal": 5103006, "discrete_dose": 5203006}
metrics = read(root / "function-and-equation-errors.csv")
main = [r for r in metrics if r["horizon"] == "2" and r["current_R"] == "1"]


def value(m, n, kind, fold, candidate, field):
    return float(next(r[field] for r in main if r["mechanism"] == m and
                      int(r["n"]) == n and r["kind"] == kind and
                      int(r["fold"]) == fold and r["candidate"] == candidate))


bridge_rows = []
adjoint_rows = []
for m in labels:
    for n in (4000, 20000):
        bridge_rows.append([labels[m], f"{n:,}"] +
                           [f"{value(m, n, 'beta', j, 'ensemble', 'rmse'):.3f}" for j in (1, 2, 3)])
        adjoint_rows.append([labels[m], f"{n:,}"] +
                            [f"{value(m, n, 'adjoint', j, 'ensemble', 'equation_rmse'):.3f}" for j in (1, 2, 3)])

weight_rows = []
boundary_rows = []
estimator_rows = []
eif_rows = []
failures = []
for m, seed in seeds.items():
    for n in (4000, 20000):
        folder = root / f"{m}-n{n}-seed{seed}"
        weights = sum((read(folder / f"weights-fold{j}.csv") for j in (1, 2, 3)), [])
        for kind, displayed in (("beta", "Bridge"), ("adjoint", "Adjoint")):
            for j in (1, 2, 3):
                selected = {r["candidate"]: float(r["weight"]) for r in weights if
                            r["kind"] == kind and r["horizon"] == "2" and int(r["fold"]) == j}
                weight_rows.append([labels[m], f"{n:,}", displayed, j] +
                                   [f"{100 * selected.get(c, 0):.1f}%" if c in selected else "—"
                                    for c in ("sieve_md", "landweber", "pmmr", "saturated_cv")])
        for j in (1, 2, 3):
            for row in read(folder / f"boundary-attempts-fold{j}.csv"):
                boundary_rows.append([labels[m], f"{n:,}", row["kind"], row["horizon"], j,
                                      row["learner_fold"], row["candidate"], row["lower"], row["upper"]])
        status = folder / "STATUS.txt"
        if not status.exists():
            estimator_rows.append([labels[m], f"{n:,}", "SDR/TMLE", "Still running when report was written"])
        else:
            errors = read(folder / "estimator-errors.csv")
            failures.extend(errors)
            for method in ("sdr", "tmle"):
                error = next((r for r in errors if r["estimator"] == method), None)
                results = [r for r in read(folder / "estimates.csv") if r["estimator"] == method]
                if error:
                    estimator_rows.append([labels[m], f"{n:,}", method.upper(), "Failed; see saved error and log"])
                elif results:
                    for r in results:
                        estimator_rows.append([labels[m], f"{n:,}", method.upper(),
                                               f"Time {int(r['horizon']) + 1}: estimate {float(r['estimate']):.6f}; truth {float(r['truth']):.6f}; 95% interval [{float(r['lower']):.6f}, {float(r['upper']):.6f}]"])
                        eif_rows.append([labels[m], f"{n:,}", method.upper(), int(r['horizon']) + 1,
                                         f"{float(r['mean_sequential_eif']):+.6f}",
                                         f"{float(r['mean_bridge_adjoint_eif']):+.6f}"])
                else:
                    estimator_rows.append([labels[m], f"{n:,}", method.upper(), "No successful result recorded"])

support = [r for r in read(root / "saturated-sieve-support.csv") if
           r["horizon"] == "2" and r["current_R"] == "1"]
support_rows = []
for m in labels:
    for n in (4000, 20000):
        for kind, name in (("beta", "Bridge"), ("adjoint", "Adjoint")):
            rows = [r for r in support if r["mechanism"] == m and int(r["n"]) == n and r["kind"] == kind]
            count = [int(r["represented_combinations"]) for r in rows]
            missing = [100 * float(r["missing_probability"]) for r in rows]
            support_rows.append([labels[m], f"{n:,}", name,
                                 f"{min(count)}–{max(count)} / {rows[0]['population_combinations']}",
                                 f"{min(missing):.2f}%–{max(missing):.2f}%"])

old_path = root.parent / "missing-both-study-v2/diagnostics/binary-numerical-bridge-comparison/bridge-error-comparison.csv"
old = read(old_path)
comparison = []
for m, label in labels.items():
    for j in (1, 2, 3):
        before = next(float(r["rmse"]) for r in old if r["mechanism"] == label and
                      r["method"] == "Ensemble" and r["R2"] == "1" and int(r["fold"]) == j)
        comparison.append([label, j, f"{before:.3f}", f"{value(m, 4000, 'beta', j, 'ensemble', 'rmse'):.3f}"])

text = f"""# Cell-score and penalty revision: four saved datasets

The requested bridge and adjoint fits are complete for both examples at n = 4,000 and n = 20,000, with three outer training fits per dataset. The numerical-dose bridge improves in every n = 4,000 fit relative to the previous implementation. Binary recovery is mixed, including a substantially worse first training fit. Adjoint equation errors remain visible. These results do not establish corrected confidence-interval coverage or convergence.

The repeated-sample simulation remains stopped at 118 completed datasets. These four datasets and the earlier million-person checks are separate diagnostics. The [binary million-person adjoint check](../missing-both-study-v2/diagnostics/binary-dose-n1000000-seed5103006-adjoint/README.md) retains its original package versions and scores.

## What changed

- All cell validation scores use distinct observations within each complete conditioning combination. Self-products are excluded in the ensemble score, saturated bridge penalty selection, and `lmtp:::bridge_response_engine`. Treatment A and all observed-history columns remain in the conditioning matrix.
- Singleton conditioning combinations supply no distinct pair. If every combination in a validation sample is a singleton, the scorer fails explicitly. No V-statistic fallback is used.
- The averaged raw U-statistic matrix can have negative eigenvalues. Before choosing ensemble weights, the package replaces those eigenvalues by zero. This is the requested positive semidefinite projection. Raw scores, eigenvalues and the size of the adjustment are saved. Penalty selection uses the raw scalar U-statistic, which can be negative.
- Saturated L1 penalty selection starts with 21 positive scales from `1e-5` through `1`, spaced by quarter decades and divided by the square root of the training person count. A minimum at a boundary is a failed attempt; that side is extended by one decade with four new values. Only an interior minimum is accepted. After 16 unresolved extensions the fit fails and carries its evaluated grid and boundary history. No zero-penalty comparison was run.
- The adjoint library contains saturated joint-category `sieve_md`, `landweber`, and `pmmr`; `saturated_l1` has been removed. The bridge library contains those three plus saturated L1. The non-L1 candidates keep their previously documented fixed regularization defaults. MARS remains only in outcome regression, alongside saturated L1 and the mean. Treatment-ratio regression uses saturated L1 and the mean.
- Raw cmbridge functions are unbounded. Bridge predictions supplied to the longitudinal estimator use the wider interval [-100, 100], as authorized. The one-step correction uses the same applied predictions. The logistic TMLE update is retained; no Gaussian mean update was introduced.

## Datasets, splits and truth

Binary treatment uses seed 5103006; numerical dose uses seed 5203006 and values 0, 1, 2, 3. Missingness occurs at both follow-ups. The n = 4,000 datasets reuse replication 6 and its saved outer and learner splits. The n = 20,000 datasets use the same respective seeds and new three-way outer and learner assignments. Using the same seed at another sample size does not make these nested datasets.

All nuisance functions are estimated. The exact generating distribution is used only to check assumptions, compute truth and evaluate the saved fitted functions. Its valid bridge and adjoint solutions were rechecked. The complete generating distributions and true parameter calculations are unchanged; θ₃ is 0.649747477399 for binary treatment and 0.633624046382 for numerical dose.

The primary figures concern the final outcome, where β₂ and λ₂ are estimated from the time-2 history. R₂ = 1 means the intermediate health vector is observed; R₂ = 0 means it is unobserved and A₂ follows its deterministic assignment. The plotted bridge is unique on the reachable support where R₂ = 1. The binary adjoint is also unique there. For numerical dose, multiple adjoints may solve the equation, so its primary adjoint plot evaluates the defining conditional equation rather than distance from one reference solution.

Each figure has n = 4,000 in the top row and n = 20,000 in the bottom row. Columns show all three outer training fits. Each point is one reachable predictor combination, with equal point sizes; error summaries weight combinations by their exact population probability. Full-range figures retain every point.

## Bridge fitted versus true

The following values are root mean squared errors among measured final outcomes conditional on R₂ = 1. They compare the raw ensemble bridge with the unique true bridge.

{table(['Mechanism', 'n', 'Training fit 1', 'Training fit 2', 'Training fit 3'], bridge_rows)}

![Binary ensemble bridge](binary_longitudinal-bridge-ensemble.png)

[Binary PDF](binary_longitudinal-bridge-ensemble.pdf), [restricted view](binary_longitudinal-bridge-ensemble-central.png), and [saturated sieve alone](binary_longitudinal-bridge-sieve_md.png).

![Numerical-dose ensemble bridge](discrete_dose-bridge-ensemble.png)

[Numerical-dose PDF](discrete_dose-bridge-ensemble.pdf), [restricted view](discrete_dose-bridge-ensemble-central.png), and [saturated sieve alone](discrete_dose-bridge-sieve_md.png).

The same n = 4,000 datasets before and after this combined revision give:

{table(['Mechanism', 'Training fit', 'Previous bridge error', 'Revised bridge error'], comparison)}

Several changes were made together, so this comparison cannot attribute the improvement or deterioration to any one change. A larger sample reduces the numerical-dose bridge errors in all three fits, but the binary fits still include rare large errors. The maximum binary bridge error at n = 20,000 is 9.005 in training fit 3.

## Adjoint and its conditional equation

The equation checked is the paper's condition:

$$
E\{{\lambda_t(A_t,H_t)\mid R_{{t+1}}=1,H_t,C_{{t+1}}\}}
=\phi_t(H_t,C_{{t+1}}).
$$

The fitted left side is integrated over the exact conditional distribution and compared with the true right side. This checks the fitted adjoint, including error inherited from its estimated treatment-ratio loading. It is not its empirical training residual. The reported equation errors are root mean squared errors among measured final outcomes conditional on R₂ = 1.

{table(['Mechanism', 'n', 'Training fit 1', 'Training fit 2', 'Training fit 3'], adjoint_rows)}

![Binary adjoint fitted versus true](binary_longitudinal-adjoint-ensemble.png)

[Binary adjoint PDF](binary_longitudinal-adjoint-ensemble.pdf), [saturated sieve alone](binary_longitudinal-adjoint-sieve_md.png), and [binary adjoint equation](binary_longitudinal-adjoint-equation-ensemble.png).

![Numerical-dose adjoint equation](discrete_dose-adjoint-equation-ensemble.png)

[Numerical-dose equation PDF](discrete_dose-adjoint-equation-ensemble.pdf) and [saturated sieve equation](discrete_dose-adjoint-equation-sieve_md.png). Differences from the selected reference numerical adjoint are saved in the CSV, but are not used to judge whether it solves the equation.

## Learner weights for the final outcome

The weights below were selected with the U-statistic and the requested matrix projection. A dash means the candidate is absent from that library.

{table(['Mechanism', 'n', 'Function', 'Training fit', 'Sieve', 'Landweber', 'PMMR', 'Saturated L1'], weight_rows)}

The projection materially affects some selections. For example, the first binary n = 4,000 bridge matrix has two negative eigenvalues; its projection adjustment is 0.106. At binary n = 20,000, the second adjoint fit selects sieve alone; its raw sieve score is -3.905, compared with 6.724 using the historical squared cell-mean score. These validation scores remain variable despite removing self-products. [All raw and historical scores](scorer-audit.csv) and [scoring counts](scoring-cell-counts.csv) document the calculation. Historical scores were computed only for this audit and did not select any revised fit.

## Penalty boundary checks and truncation

All completed primary nuisance penalty selections have interior minima. The recorded failed boundary attempts are:

{table(['Mechanism', 'n', 'Function', 't', 'Outer fit', 'Learner fit (0 = full refit)', 'Candidate', 'Lower boundary', 'Upper boundary'], boundary_rows) if boundary_rows else 'No initial boundary minimum occurred in the completed primary nuisance fits.'}

For the binary n = 4,000 final bridge refit in training fit 1, the initial minimum was `1e-5`. The extended grid reached `1e-6`; `1e-5` then became an interior minimum because all newly added smaller scales had larger losses. The failed initial attempt and the accepted extended grid are both saved.

There was no [-100, 100] truncation in any primary outer-validation bridge prediction or population ensemble prediction, and no inner training bridge value was recorded at the interval limits. Thus widening the interval did not change the bridge curves shown here. [Truncation counts](applied-bridge-truncation.csv) retain the checks. Deeper response-construction fits are distinct from these primary nuisance fits.

## Actual fitted saturated dictionaries

The complete joint-category family can represent the required solutions. The fitted sieve dictionary contains only combinations represented in its training data, with a common constant used for unseen combinations. The following counts and missing probabilities concern R₂ = 1. Bridge probabilities are conditional on measured final outcomes; adjoint probabilities refer to the natural treatment/history distribution at that intermediate visit.

{table(['Mechanism', 'n', 'Function', 'Represented / reachable combinations', 'Probability of unseen combinations'], support_rows)}

Most of the bridge squared error occurs in represented combinations. Unseen combinations therefore cannot alone account for the poor sieve or ensemble recovery. These diagnostics distinguish representability of the complete family from the accuracy and support of each finite fitted function. [Complete support audit](saturated-sieve-support.csv).

## Additional full-estimator checks

The nuisance-function plots above completed independently of these additional SDR/TMLE attempts.

{table(['Mechanism', 'n', 'Estimator', 'Outcome'], estimator_rows)}

Both n = 4,000 full-estimator attempts failed for each example: deeply nested bridge-validation samples had no repeated full conditioning combination, so the distinct-pair score could not be computed. SuperLearner reported a failed saturated candidate, and the study rejected the incomplete library fit. The original U-statistic error is retained in the dataset log. This is a limitation of the current recursive validation construction at this sample size. Neither changing A, silently using a V-statistic, nor discarding failed datasets was used to conceal it. The broader study should not resume before that validation construction is addressed.

For successful fits, the two component means below center the sequential component at the true parameter. Their sum is the estimate minus truth. Individual contributions and the separately centered estimated gradient are retained in each dataset's `eif-components.csv`.

{table(['Mechanism', 'n', 'Estimator', 'Outcome time', 'Mean sequential component', 'Mean bridge/adjoint component'], eif_rows) if eif_rows else 'No successful full-estimator component means are available yet.'}

## Verification and saved files

The cmbridge suite and the full lmtp suite passed. Regression checks cover the exact distinct-pair formula, treatment retained in conditioning, removal of self-product variance, matrix projection, both penalty-extension directions, unresolved-boundary errors, and isolation of training responses from validation observations.

An independent audit recomputed every primary raw ensemble matrix from its saved residuals and full conditioning inputs, agreeing within `1e-10`; projected matrices had no eigenvalue below `-1e-10`. The audit also checked the adjoint libraries, applied prediction values and generating assumptions.

Package versions are cmbridge `0.3.0.9006` and lmtp `1.6.0.9006`. [Frozen sources](source/), private library, test/install logs and per-dataset checkpoints are retained here. Installable source archives are [cmbridge](cmbridge_0.3.0.9006.tar.gz) and [lmtp](lmtp_1.6.0.9006.tar.gz). Each dataset directory contains fitted objects, predictions, conditional equations, penalty grids, boundary attempts, weights and estimator status. The two package repositories are the sibling `cmbridge` and `lmtp` directories; the analysis code lives in `cmbridge-tests/simulations/observed_history_corrected`.

[All function values](all-function-predictions.csv), [all conditional-equation values](all-equation-predictions.csv), [errors by candidate and fit](function-and-equation-errors.csv), [cmbridge tests](cmbridge-tests.log), [lmtp tests](lmtp-all-tests.log), and [audit log](scorer-audit.log).

The plot and report can be regenerated from `cmbridge-tests` without fitting another dataset:

```sh
Rscript simulations/observed_history_corrected/plot_selection_revision.R
python3 simulations/observed_history_corrected/report_selection_revision.py
```

The independent fit audit additionally requires the private revised library:

```sh
R_LIBS=results/observed-history-corrected/selection-revision/library:results/observed-history-corrected/missing-both-study-v2/library Rscript simulations/observed_history_corrected/audit_selection_revision.R
```

Report generated {datetime.now().astimezone().isoformat(timespec='seconds')} from saved diagnostics. Initial diagnostic export/setup errors are preserved with `initial-*` filenames; those were corrected without altering the fitted nuisance objects.
"""
(root / "REPORT.md").write_text(text)
(root / "README.md").write_text("# Selection revision diagnostics\n\n[Results, plots and limitations](REPORT.md).\n\nFour isolated datasets; the repeated-sample study remains stopped.\n")
print(root / "REPORT.md")
