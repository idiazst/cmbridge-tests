# cmbridge tests

Independent large-sample validation of the bridge and adjoint learners in
[`idiazst/cmbridge`](https://github.com/idiazst/cmbridge).

Each GitHub Actions run:

1. installs a fresh R release on Ubuntu;
2. clones and installs the current `cmbridge` package;
3. records the exact package commit SHA;
4. runs known-truth large-sample tests for the bridge equation;
5. runs known-truth large-sample tests for the adjoint equation;
6. validates cross-validated ensemble selection with exactly one correctly specified candidate;
7. stores CSV results, seed-1 truth-vs-estimate curves, PDFs, and `sessionInfo()`;
8. commits the results under `results/run-<github-run-id>/`.

The validation workflow runs manually (`workflow_dispatch`) or when validation
scripts or its workflow change. Commits containing only results do not trigger
new runs.

## Ensemble selection validation

`scripts/validate_ensemble.R` evaluates both bridge and adjoint ensembles over
sample sizes 2,000, 20,000, and 200,000 and five seeds, with three inner folds.
Two scenarios change which candidate is correctly specified:

| Truth | Correct candidate | Misspecified candidates |
|---|---|---|
| Cubic polynomial | Cubic sieve minimum distance | Linear Landweber; one-center Gaussian PMMR |
| Sum of five Gaussian kernels | PMMR in the exact five-center Gaussian span | Linear sieve; quadratic Landweber |

Correct specification is about the candidate's function class, rather than the
method name. The truth lies outside the other two classes. Candidate controls,
the scoring kernel, and acceptance criteria are fixed in the script; the
ensemble never sees the true function or true candidate label during fitting.
The same Gaussian scoring kernel is used across all candidates, with a
training-derived 81-center Nyström approximation. Bridge V is missing when M=0.
For the adjoint, the observation probability depends on B beyond V, and the
known common loading makes the true conditional-moment equation hold exactly.

At the largest size, every seed and both scenarios must:

- give the correct candidate at least 90% weight and select it by largest weight;
- recover the true function with grid RMSE below 0.06;
- have ensemble RMSE below half that of the best misspecified candidate;
- satisfy the simplex optimizer's scaled KKT-gap tolerance of 1e-8.

Smaller sizes describe finite-sample behavior and need not select a simplex
vertex. Results include all candidate weights, candidate RMSE, cross-validated
moment losses, summary tables, curves, a PDF, configuration, and R session
information. Passing these deliberately specified diagnostics does not establish
performance when all learners are misspecified or identification is weak.

```bash
Rscript scripts/validate_ensemble.R results/ensemble
```

For a shorter diagnostic, set `CMBRIDGE_ENSEMBLE_SIZES` and
`CMBRIDGE_ENSEMBLE_SEEDS` to comma-separated values. Acceptance criteria are
applied to the largest requested size.

## Run again

From a terminal:

```bash
gh workflow run validate.yml --repo idiazst/cmbridge-tests
```

Then inspect the newest run:

```bash
gh run list --repo idiazst/cmbridge-tests --workflow validate.yml --limit 5
```
