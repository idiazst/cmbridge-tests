"""Report the completed saved n=4000 checks and all fitting failures."""
import collections
import csv
import json
import math
import os
from pathlib import Path

root = Path(os.environ['SIM_OUTPUT'])
reference = Path(os.environ['SIM_REFERENCE'])
def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))
def table(headers, rows):
    return ('| ' + ' | '.join(headers) + ' |\n| ' + ' | '.join(['---'] * len(headers)) +
            ' |\n' + ''.join('| ' + ' | '.join(map(str, row)) + ' |\n' for row in rows) + '\n')
names = {'binary_longitudinal': 'Binary treatment', 'discrete_dose': 'Numerical dose'}
estimates, weights, candidate_failures, trial_failures = [], [], [], []
for mechanism in names:
    seed = 5103006 if mechanism == 'binary_longitudinal' else 5203006
    folder = root / f'{mechanism}-n4000-seed{seed}'
    assert 'Estimator error rows 0' in (folder / 'STATUS.txt').read_text()
    rows = read(folder / 'estimates.csv')
    assert len(rows) == 4
    for row in rows:
        assert abs(float(row['mean_sequential_eif']) + float(row['mean_bridge_adjoint_eif']) -
                   float(row['estimate']) + float(row['truth'])) < 1e-10
    estimates.extend(rows)
    weights.extend(dict(row, mechanism=mechanism) for row in read(folder / 'all-weights.csv')
                   if row['kind'] in ['beta', 'adjoint'])
    candidate_failures.extend(read(folder / 'nested-candidate-failures.csv'))
    trial_failures.extend(read(folder / 'nested-penalty-trial-failures.csv'))
metrics = read(root / 'function-and-equation-errors.csv')
before = read(reference / 'function-and-equation-errors.csv')
def rmse(values):
    grouped = collections.defaultdict(list)
    for row in values:
        if row['kind'] == 'beta' and row['horizon'] == '2' and row['current_R'] == '1' and row['n'] == '4000':
            grouped[(row['mechanism'], row['candidate'])].append(float(row['rmse']) ** 2)
    return {key: math.sqrt(sum(value) / len(value)) for key, value in grouped.items()}
now, previous = rmse(metrics), rmse(before)
text = '''# Ensemble stability: saved n = 4,000 datasets

Both binary-treatment and numerical-dose datasets completed for SDR and logistic TMLE at both outcome times. The data, seeds (5103006 and 5203006), and shared person assignments match the saved checks. Both follow-up health vectors can be missing. This comparison uses one dataset per mechanism and does not assess repeated-sample coverage.

Packages: **cmbridge 0.3.0.9020** and **lmtp 1.6.0.9021**. All three bridge candidates use inverse-expit with unrestricted real coefficients. All three adjoint candidates remain unrestricted with identity link. Saturated L1 remains only in SuperLearner classification and regression; MARS remains only in regression. Neither the DGP nor the true parameters changed. No generating probabilities or true coefficients were supplied to fitted learners.

## Repairs and final learner configuration

The bridge sieve now includes a fixed intercept and every main effect, followed by joint-category deviations. New combinations retain the main-effect prediction. Ridge penalizes scaled link coefficients, with 100 times the main-effect penalty on category deviations and no penalty on the intercept. Training-only cross-validation chooses a positive penalty from the 41-point initial grid; boundary minima extend the grid until interior or a recorded failure. The adjoint sieve retains saturated joint-category bases and centered function-value ridge.

The bridge Landweber includes a fixed intercept and every main effect. Its full-rank coordinate transformation improves conditioning without dropping any column. The previous quadratic basis remains available and its diagnostic attempts are preserved separately. The adjoint Landweber remains spline based. PMMR retains its Gaussian-feature approximation and fixed penalty. The supplied histories and current treatment are retained throughout.

An algebraic audit checks the actual sieve and Landweber bases on every reachable history, including missed visits. Both represent a valid bridge solution to numerical precision, even when initialized from one training combination. Class containment does not assert that estimated coefficients are accurate.

Gaussian kernel U-statistics, with training-only scaling and centers, are used for both ensemble weights and sieve penalty selection. All self-products are excluded. Scalar penalty scores remain raw, including negative values; only the ensemble Gram is projected to positive semidefinite before simplex optimization. The cell U-statistic remains available and also excludes self-products.

Large ridge values previously obscured the unpenalized common direction. The centered-cell solve and linked main-effect sieve now use separate common/intercept and penalized directions. Optimizer refinement retains the original convergence checks. A very small nested training fit also produced overflowing validation scores. Such penalty trials are now recorded as failed and cannot be chosen. Failed candidates are excluded consistently across folds and the final refit, with zero weight and retained failure details. All-candidate failures remain estimator errors.

Successful penalty batches retain vectorized evaluation and glmnet warm starts. Individual trials are retried only when an error in the batch prevents identifying failures. This repairs a runtime regression introduced by evaluating every penalty separately. The CV criterion, split assignments, positive grid, and interior minimum rule are unchanged.

Only beta and density-ratio fits with identical training people, training data, learner configuration, and split callback are shared between SDR and TMLE. Sequential regressions and responses remain estimator specific. Responses and upstream nuisance fits continue to be rebuilt inside candidate and penalty training samples. Cached outer fits predate the final error-handling repair; independent fresh refits reproduce their predictions, weights, selected penalties, and split assignments exactly. Previous attempts are retained separately and are not counted as completed estimator results.

## Fitted versus true bridges

Plots show raw predictions for every reachable predictor combination at time 3 with R₂ = 1, where the bridge solution is unique. Equal axes and the diagonal show y = x. The ensemble clusters much more closely around the diagonal, with visible differences between the three training fits.

![Binary treatment ensemble](binary_longitudinal-bridge-ensemble.png)

![Numerical-dose ensemble](discrete_dose-bridge-ensemble.png)

![Binary treatment: all learners](binary_longitudinal-all-bridge-learners.png)

![Numerical dose: all learners](discrete_dose-all-bridge-learners.png)

Population-probability-weighted RMSE below averages squared errors over the three outer fits. The earlier inverse-expit comparison and the current revision use the same datasets and splits. Several repairs were applied together, so their separate effects cannot be inferred from this table.

'''
candidates = ['sieve_md', 'landweber', 'pmmr', 'ensemble']
text += table(['Example', 'Revision', 'Sieve', 'Landweber', 'PMMR', 'Ensemble'],
              [[label, stage] + [f'{values[(mechanism,c)]:.4f}' for c in candidates]
               for mechanism, label in names.items() for stage, values in [('Previous', previous), ('Current', now)]])
text += '## Ensemble weights at time 3\n\n'
for kind, label in [('beta', 'Bridge'), ('adjoint', 'Adjoint')]:
    rows = []
    for mechanism in names:
        for fold in [1, 2, 3]:
            selected = {r['candidate']: float(r['weight']) for r in weights if r['mechanism'] == mechanism and
                        r['kind'] == kind and int(r['horizon_or_depth']) == 2 and int(r['fold']) == fold}
            rows.append([names[mechanism], fold] + [f'{selected[c]:.4f}' for c in candidates[:3]])
    text += f'### {label}\n\n' + table(['Example', 'Training fit', 'Sieve', 'Landweber', 'PMMR'], rows)
text += '## Estimates and mean EIF components\n\nThe sequential component is centered at the true parameter. Its mean plus the bridge/adjoint mean equals estimate minus truth.\n\n'
text += table(['Example', 'Outcome time', 'Estimator', 'Truth', 'Estimate', 'SE', 'Sequential mean', 'Bridge/adjoint mean'],
              [[names[r['mechanism']], int(r['horizon']) + 1, r['estimator'].upper()] +
               [f'{float(r[k]):.6f}' for k in ['truth','estimate','se','mean_sequential_eif','mean_bridge_adjoint_eif']]
               for r in estimates])
misses = sum(not float(r['lower']) <= float(r['truth']) <= float(r['upper']) for r in estimates)
text += f'{misses} of eight pointwise intervals exclude the truth in these two datasets.\n\n'
diagnostic_path = root / 'complete-population-diagnostics.csv'
if diagnostic_path.exists():
    grouped = collections.defaultdict(lambda: collections.defaultdict(float))
    for row in read(diagnostic_path):
        key = (row['mechanism'], row['estimator'], row['horizon'])
        for column in ['population_bias', 'sequential_remainder', 'bridge_remainder']:
            grouped[key][column] += float(row[column]) * float(row['fold_weight'])
    text += '## Examination of the excluded intervals\n\n'
    text += 'Both excluded intervals concern the numerical-dose outcome at time 2 in the same dataset. To examine this, evaluate each fitted estimator contribution over every point of the known generating distribution and weight the three training fits by their validation sample sizes. The resulting expected estimation errors below are much smaller than the observed error of about -0.028. This points toward variation in the evaluation observations in this dataset; it does not establish repeated-sample coverage or exclude other finite-sample effects. The full study will examine their frequency and standard-error calibration.\n\n'
    text += table(['Example', 'Outcome time', 'Estimator', 'Expected estimation error'],
                  [[names[m], int(h) + 1, e.upper(), f'{value["population_bias"]:.6f}']
                   for (m, e, h), value in grouped.items()])
    text += 'The [complete population diagnostics](complete-population-diagnostics.csv) also retain the separate second-order contributions and conditional-equation errors. Those second-order contributions are different from the sample mean EIF components in the preceding table. The sequential second-order contribution compares the sequential estimate with the parameter defined using its fitted bridge; the bridge/adjoint contribution then accounts for the bridge error and correction. Their sum equals the expected estimation error. True functions are used only for these diagnostics.\n\n'
text += '## Recorded fitting failures\n\n'
text += f'The final checks have zero estimator error rows. Their shared nested bridge caches record {len(candidate_failures)} candidate-failure events and {len(trial_failures)} failed penalty-trial records. Records are not counts of independent people or distinct statistical failures: the same training fit may have multiple trials or scoring folds. All records are retained in the case folders.\n\n'
text += '''## Audits and full study

Audits verify unchanged data and outer/learner assignments, training-only penalty splits, positive interior selected penalties, direct reconstruction of public estimates, EIF component sums, kernel U-statistics, and PSD ensemble matrices. cmbridge predictions are not clipped after fitting. The sequential estimator retains the approved wider applied interval [-100,100]; changed values are recorded in [applied-bridge-truncation.csv](applied-bridge-truncation.csv). Plots retain raw values.

Bridge solutions at R₂ = 0, and some numerical-dose adjoints, need not be unique. Their defining conditional equations remain the relevant diagnostic. Difference from a selected reference solution alone is not treated as a fitting error.

[Full-support class audit](bridge-function-class.csv) · [Fresh outer-fit equality](outer-fit-reuse-audit.csv) · [Scoring audit](scorer-audit.csv) · [Penalties](sieve-selected-penalties.csv) · [Function and equation errors](function-and-equation-errors.csv) · [All predictions](all-function-predictions.csv).

The approved full study has its own frozen sources and package archives. It includes both mechanisms at n = 500, 1,000, and 4,000, with 200 replications each (1,200 datasets). It compares SDR and logistic TMLE with the one-step bridge correction, using ten workers. Misspecification comparisons remain deferred. Full-study failures, bias, interval coverage, component means, and convergence plots will be reported from the actual results.

'''
github = json.loads((root / 'validation/github-validation.json').read_text())
text += f"## GitHub validation\n\n[Run {github['run_id']}](https://github.com/idiazst/cmbridge-tests/actions/runs/{github['run_id']}) checks package commit **{github['package_commit']}**. Recorded outcomes:\n\n```text\n{github['outcomes'].strip()}\n```\n\n"
text += '![GitHub ensemble recovery](github-run/ensemble/truth_vs_estimate.png)\n\n![GitHub inverse-expit and adjoint checks](github-run/linked-bridge/truth_vs_estimate.png)\n'
(root / 'REPORT.md').write_text(text)
(root / 'README.md').write_text('# Ensemble stability checks\n\nSee [REPORT.md](REPORT.md) for both n = 4,000 examples, all learner plots, weights, EIF means, and recorded failures.\n')
print(f'Reported {len(estimates)} estimates and all recorded nested failures.')
