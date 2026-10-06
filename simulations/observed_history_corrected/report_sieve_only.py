"""Sieve-only figures and saved conditional-moment ensemble weights."""
from pathlib import Path
import csv
import os

root = Path(os.environ.get("SIM_OUTPUT", "results/observed-history-corrected/selection-revision"))
labels = {"binary_longitudinal": "Binary treatment", "discrete_dose": "Numerical dose"}
seeds = {"binary_longitudinal": 5103006, "discrete_dose": 5203006}


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


weights = []
rows = []
for mechanism, seed in seeds.items():
    for n in (4000, 20000):
        folder = root / f"{mechanism}-n{n}-seed{seed}"
        for fold in (1, 2, 3):
            for row in read(folder / f"weights-fold{fold}.csv"):
                weights.append({"mechanism": mechanism, "n": n, **row})
        for kind, name in (("beta", "Bridge"), ("adjoint", "Adjoint")):
            values = []
            for candidate in ("sieve_md", "landweber", "pmmr", "saturated_cv"):
                entries = [next((float(r["weight"]) for r in weights if
                                  r["mechanism"] == mechanism and int(r["n"]) == n and
                                  r["kind"] == kind and r["horizon"] == "2" and
                                  int(r["fold"]) == fold and r["candidate"] == candidate), None)
                           for fold in (1, 2, 3)]
                values.append("—" if all(v is None for v in entries) else
                              " / ".join(f"{100 * v:.1f}" for v in entries))
            rows.append([labels[mechanism], f"{n:,}", name, *values])

with (root / "conditional-moment-ensemble-weights.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(weights[0]))
    writer.writeheader()
    writer.writerows(weights)

metrics = read(root / "function-and-equation-errors.csv")
errors = []
for mechanism in labels:
    for n in (4000, 20000):
        for kind, name, field in (("beta", "Bridge function", "rmse"),
                                  ("adjoint", "Adjoint equation", "equation_rmse")):
            selected = [next(r for r in metrics if r["mechanism"] == mechanism and
                             int(r["n"]) == n and r["horizon"] == "2" and
                             r["current_R"] == "1" and r["kind"] == kind and
                             r["candidate"] == "sieve_md" and int(r["fold"]) == fold)
                        for fold in (1, 2, 3)]
            errors.append([labels[mechanism], f"{n:,}", name,
                           *[f"{float(r[field]):.3f}" for r in selected]])

output = f"""# Saturated sieve plots and ensemble weights

These figures use only the saved saturated joint-category `sieve_md` candidate.
It includes interactions of every order. They show its raw predictions, with no
ensemble averaging or clipping, from the same n = 4,000 and n = 20,000 datasets.
No new dataset or nuisance fit was run to produce this presentation.

The complete saturated function family contains the true bridge and adjoint
solutions. The finite fitted dictionary includes only training combinations;
unseen combinations use its fitted common constant. Its quadratic regularization
remains `lambda = 1e-8`, with `weight_ridge = 1e-8`. The expanded positive L1
penalty grid applies to the separate saturated L1 candidates. The sieve's fitted
dictionary support is recorded in [the support audit](saturated-sieve-support.csv).

## Ensemble weights

The entries below are percentages for training fits **1 / 2 / 3**, for the final
outcome at time 3. A dash means the candidate is absent. These are the ensemble
weights; the following plots use sieve alone irrespective of its weight.

{table(['Example', 'n', 'Function', 'Sieve', 'Landweber', 'PMMR', 'Saturated L1'], rows)}

[Exact weights for both outcome times, raw scores and projection diagnostics](conditional-moment-ensemble-weights.csv).

## Sieve bridge fitted versus true

The top row is n = 4,000 and the bottom row is n = 20,000. Columns show the three
outer training fits. Each point is one reachable predictor combination. The
blue line is y = x. These figures concern R₂ = 1, where the time-2 health vector
is observed and the bridge is unique. Full-range figures retain every point.

![Binary treatment: saturated sieve bridge](binary_longitudinal-bridge-sieve_md.png)

[Binary PDF](binary_longitudinal-bridge-sieve_md.pdf) and
[binary restricted view](binary_longitudinal-bridge-sieve_md-central.png).

![Numerical dose: saturated sieve bridge](discrete_dose-bridge-sieve_md.png)

[Numerical-dose PDF](discrete_dose-bridge-sieve_md.pdf) and
[numerical-dose restricted view](discrete_dose-bridge-sieve_md-central.png).

## Sieve adjoint

The binary adjoint is unique on this reachable support, so fitted versus true
function values are shown directly.

![Binary treatment: saturated sieve adjoint](binary_longitudinal-adjoint-sieve_md.png)

[Binary adjoint PDF](binary_longitudinal-adjoint-sieve_md.pdf) and
[binary defining-equation plot](binary_longitudinal-adjoint-equation-sieve_md.png).

The numerical-dose adjoint can have multiple solutions. The y = x plot below
therefore compares the conditional mean of its fitted adjoint with the true
right-hand side of the paper's defining equation. It does not require recovery
of one arbitrarily selected reference solution.

![Numerical dose: saturated sieve adjoint equation](discrete_dose-adjoint-equation-sieve_md.png)

[Numerical-dose adjoint equation PDF](discrete_dose-adjoint-equation-sieve_md.pdf).

## Sieve errors

Root mean squared errors are weighted by exact population probabilities among
measured final outcomes conditional on R₂ = 1. The bridge rows measure function
error; the adjoint rows measure error in its defining equation.

{table(['Example', 'n', 'Quantity', 'Training fit 1', 'Training fit 2', 'Training fit 3'], errors)}

The sieve bridge still has substantial errors at both sizes. Representability
of the complete saturated family does not establish accurate recovery by these
finite fitted functions. [All saved candidate errors](function-and-equation-errors.csv)
and [the full diagnostic record](REPORT.md) retain the underlying calculations.
"""
(root / "SIEVE-REPORT.md").write_text(output)
print(root / "SIEVE-REPORT.md")
print(table(['Example', 'n', 'Function', 'Sieve', 'Landweber', 'PMMR', 'Saturated L1'], rows))
