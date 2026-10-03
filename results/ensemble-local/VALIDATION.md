The bridge and adjoint ensembles were validated in actual R 4.3.3 using three inner folds, sample sizes 2,000, 20,000, and 200,000, and seeds 1 through 5. There are 60 ensemble fits in total and 20 largest-sample cases. All largest-sample cases passed the prespecified criteria in scripts/validate_ensemble.R.

The cubic scenario has exactly one correctly specified cubic sieve candidate; the alternative candidates are linear Landweber and one-center Gaussian PMMR. The Gaussian scenario has exactly one correctly specified five-center PMMR candidate; the alternatives are linear sieve and quadratic Landweber. The candidate function class, rather than the method name, determines correct specification. Truth and candidate labels are used only for simulation and evaluation, never for selecting weights.

At n=200,000:

| Truth | Problem | Mean correct weight | Minimum correct weight | Mean ensemble RMSE | Maximum ensemble RMSE | Selection rate |
|---|---|---:|---:|---:|---:|---:|
| cubic_truth | adjoint | 99.9583% | 99.7917% | 0.001389 | 0.002509 | 100% |
| cubic_truth | bridge | 97.8858% | 96.0172% | 0.027618 | 0.039089 | 100% |
| rbf_truth | adjoint | 99.8137% | 99.6060% | 0.006094 | 0.008732 | 100% |
| rbf_truth | bridge | 99.4244% | 99.2867% | 0.036727 | 0.055539 | 100% |

The correct candidate was selected by largest weight in all 20 large-sample cases. Every such case gave it at least 90% weight, had ensemble grid RMSE below 0.06 and below half the best misspecified candidate's RMSE, and met the scaled KKT-gap tolerance of 1e-8. Actual minimum correct weight was 96.0172%.

The smaller samples are reported without enforcing selection. At n=2,000, the mean correct weight in the bridge scenarios was 18.75% for the cubic truth and 39.96% for the Gaussian truth. This diagnostic therefore demonstrates the requested large-sample behavior, not uniform selection at all sample sizes.

Scoring uses the common Gaussian RKHS maximum-moment criterion with a training-derived 81-center Nystrom approximation, independent of the candidate fitting bases. Bridge scoring includes all rows without evaluating missing V when M=0. Adjoint scoring uses only complete cases and a known common loading; the observation probability depends on B beyond V. An estimated-loading callback is covered separately by package unit tests.

The package provides exact blocked Gaussian scoring as an option. Equation (22), simplex interior/boundary/singular cases, fold isolation, loading callback calls, complete-case handling, prediction combination, and random-state preservation were checked by 38 passing unit assertions. R CMD check completed with zero errors, zero warnings, and one note about locally unavailable optional suggested packages np, npiv, and gmm.

cmbridge_commit.txt records the tested local package commit. The updated GitHub Actions workflow is included, but these new commits could not be pushed because the connector rejected repository writes and the local Git remote had no authentication. No GitHub Actions run of this version is claimed.

The simulation is deliberately specified and strongly informative. It does not establish general convergence rates, performance with all candidates misspecified, behavior under weak identification, or correctness of a complete longitudinal causal estimator.
