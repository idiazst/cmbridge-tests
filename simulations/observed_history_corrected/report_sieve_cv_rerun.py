"""Report the two requested n=4000 checks without running a coverage study."""
from pathlib import Path
import collections
import csv
import math
import os
import json

root = Path(os.environ['SIM_OUTPUT'])
reference = Path(os.environ['SIM_REFERENCE'])
github = json.loads((root/'validation/github-validation.json').read_text())
def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))
def label(mechanism):
    return 'Binary' if mechanism == 'binary_longitudinal' else 'Numerical dose'
def table(headers, rows):
    return ('| ' + ' | '.join(headers) + ' |\n| ' + ' | '.join(['---']*len(headers)) +
            ' |\n' + ''.join('| ' + ' | '.join(map(str, row)) + ' |\n' for row in rows) + '\n')
def number(value):
    return f'{float(value):.6f}'

mechanisms = ['binary_longitudinal', 'discrete_dose']
candidates = ['sieve_md', 'landweber', 'pmmr', 'ensemble']
estimates, weights, previous = [], [], []
for mechanism in mechanisms:
    seed = 5103006 if mechanism == mechanisms[0] else 5203006
    folder = root/f'{mechanism}-n4000-seed{seed}'
    assert 'Estimator error rows 0' in (folder/'STATUS.txt').read_text()
    rows = read(folder/'estimates.csv')
    assert len(rows) == 4
    estimates.extend(rows)
    previous.extend(read(reference/folder.name/'estimates.csv'))
    weights.extend(dict(row, mechanism=mechanism, n='4000') for row in read(folder/'all-weights.csv')
                   if row['kind'] in ['beta', 'adjoint'])
with (root/'conditional-moment-ensemble-weights.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(weights[0]))
    writer.writeheader()
    writer.writerows(weights)

metrics = read(root/'function-and-equation-errors.csv')
oldmetrics = read(reference/'function-and-equation-errors.csv')
def rmse_rows(values):
    squared = collections.defaultdict(list)
    for row in values:
        if row['kind'] == 'beta' and row['horizon'] == '2' and row['current_R'] == '1' and row['n'] == '4000':
            squared[(row['mechanism'], row['candidate'])].append(float(row['rmse'])**2)
    return {key: math.sqrt(sum(value)/len(value)) for key, value in squared.items()}
now, before = rmse_rows(metrics), rmse_rows(oldmetrics)
errors = [[label(m), stage] + [f'{values[(m,c)]:.4f}' for c in candidates]
          for m in mechanisms for stage, values in [('Previous', before), ('Current', now)]]

penalties = read(root/'sieve-selected-penalties.csv')
penaltyrows = [[label(r['mechanism']), int(r['horizon'])+1, r['kind'], r['fold'],
                f'{float(r["lambda"]):.6g}', r['grid_values'], r['boundary_failures']]
               for r in penalties]
eifrows = []
for r in estimates:
    assert abs(float(r['mean_sequential_eif'])+float(r['mean_bridge_adjoint_eif'])-
               (float(r['estimate'])-float(r['truth']))) < 1e-10
    eifrows.append([label(r['mechanism']), int(r['horizon'])+1, r['estimator'].upper()] +
                  [number(r[k]) for k in ['truth', 'estimate', 'se', 'mean_sequential_eif', 'mean_bridge_adjoint_eif']])
misses = sum(not float(r['lower']) <= float(r['truth']) <= float(r['upper']) for r in estimates)
text = '''# Sieve penalty CV and quadratic Landweber: n = 4,000

The binary and numerical-dose examples were rerun at **n = 4,000 only**, using the same saved datasets, seeds 5103006 and 5203006, and the same three outer and learner-validation groups. Both follow-up outcomes can be missing. SDR and TMLE were fitted at both outcome times. The data-generating distributions and parameter truths are unchanged. The repeated-sample simulation and n = 20,000 checks remain stopped. Successful outer nuisance fits from the second attempt are retained because the final restart repair only affects fits that previously failed the unchanged convergence criteria. Sequential response construction, SDR and TMLE are refitted from scratch. This cache reuse and its audits are recorded in `validation/outer-fit-reuse.txt`.

Packages: **cmbridge 0.3.0.9017** and **lmtp 1.6.0.9018**. All three bridge candidates use inverse-expit with unrestricted coefficients; their convex ensemble is at least one. Adjoint candidates use identity and remain unrestricted. Saturated L1 is absent from the conditional-moment libraries and remains a regression/classification candidate in SuperLearner. MARS remains a regression candidate only. The logistic TMLE update is unchanged.

## Changes

The sieve function-value ridge penalty is selected by cross-validation, starting with 41 positive values from 1e-10 to 1 in quarter-decade steps. Selection uses the cell U-statistic with all conditioning variables, including treatment. Penalty validation occurs entirely within the corresponding training sample, using the common person-level split function. The grid is extended if its minimum is at a boundary; the attempted boundary choice is recorded as a failure. An unresolved boundary selection raises an error. The final refit and each ensemble training fit select their own penalties. The adjoint loading is also rebuilt within penalty-training splits. The table below contains actual selected penalties rather than fixed defaults.

The Landweber bridge now uses a fixed quadratic target basis: intercept, all main effects, all squares, and every pairwise interaction among the supplied predictors. Its conditioning basis uses joint-category indicators. Early stopping still regularizes the fitted coefficients. The adjoint Landweber basis is unchanged.

CV exposed a small nested inverse-expit sieve fit whose coefficients appeared stationary while the function-scale gradient was too large. The optimizer now refines such fits with coefficient-step scaling based on the link derivative. This changes optimization step sizes only. Coefficients remain unrestricted; the objective and convergence checks are unchanged. After an optimizer restart revives a flat coefficient, step scaling is refined again at the new location. The failing examples are retained as regression fixtures. Local package and lmtp integration tests pass. The first n = 4,000 attempt also exposed validation subsets containing only singleton conditioning cells. Their U-statistic is the zero matrix, since no distinct within-cell pairs exist. The scorer now returns that matrix and records paired-row counts and an informative flag rather than throwing. An entirely uninformative grid uses the existing deterministic interior tie rule; it does not identify an optimal penalty. The first two attempts are retained in the adjacent `sieve-cv-quadratic-landweber` and `sieve-cv-quadratic-landweber-v2` folders.

## Does the Landweber bridge class contain the truth?

An algebraic check evaluated the actual package basis on every reachable combination at both time points, including histories with R₂ = 0 and R₂ = 1. It represents the true bridge throughout the support, with maximum absolute error below 1.4e-14. The basis specification is identical whether initialized from one arbitrary row or the full support. These checks verify function-class membership; they do not replace the fitted coefficients with known functions or measure an oracle estimator.

'''
classrows = [[label(r['mechanism']), int(r['horizon'])+1, r['current_R'],
              r['reachable_combinations'], r['basis_columns'],
              f'{float(r["maximum_bridge_representation_error"]):.3g}']
             for r in read(root/'landweber-function-class.csv')]
text += table(['Example', 'Outcome time', 'Current R', 'Combinations', 'Basis columns', 'Maximum representation error'], classrows)
text += '''The sieve still builds its joint-category dictionary from its own training sample. It contains the truth on represented combinations, but can omit population combinations. CV does not change this dictionary. Missing combinations may also affect estimates at represented combinations through the conditional-moment equations, so the fraction of error on represented combinations cannot by itself rule out this explanation.

## Fitted versus true bridges

These figures display all three candidates and their ensemble at outcome time 3, on histories with R₂ = 1, where β₂ is unique. They show every reachable combination and all three outer training fits. The plotted estimates are untruncated. Each panel uses equal horizontal and vertical scales and the dashed line y = x.

![Binary: true versus fitted bridge](binary_longitudinal-all-bridge-learners.png)

![Numerical dose: true versus fitted bridge](discrete_dose-all-bridge-learners.png)

[Binary PDF](binary_longitudinal-all-bridge-learners.pdf) · [Numerical-dose PDF](discrete_dose-all-bridge-learners.pdf).

The table reports population-probability-weighted root mean squared error, averaging squared errors across the three fits. Both changes were applied together; their separate effects cannot be identified from this comparison.

'''
text += table(['Example', 'Revision', 'Sieve', 'Landweber', 'PMMR', 'Ensemble'], errors)
text += '''Sieve error decreased in both datasets but substantial disagreement with y = x remains. The current Landweber fit is less accurate than the previous spline fit, despite the new class containing the true bridge. The ensemble improves very slightly in the numerical-dose example and becomes less accurate in the binary example. These results do not establish that CV has resolved the bridge estimation problem.

Every plotted Landweber fit reached its limit of 2,000 iterations without satisfying its coefficient-change tolerance. Its training moment loss decreased, but the fit is still regularized by iteration stopping. This provides a concrete reason to investigate its remaining shrinkage; it does not show that iteration stopping explains all its error. The iteration limit was retained for this requested comparison. [Stopping diagnostics](landweber-stopping.csv).

'''
text += '## Ensemble weights at outcome time 3\n\n'
for kind, title in [('beta', 'Bridge'), ('adjoint', 'Adjoint')]:
    rows = []
    for m in mechanisms:
        for fold in [1,2,3]:
            z = {r['candidate']: float(r['weight']) for r in weights if r['mechanism'] == m and
                 r['kind'] == kind and int(r['horizon_or_depth']) == 2 and int(r['fold']) == fold}
            assert set(z) == set(candidates[:3])
            rows.append([label(m), fold] + [f'{z[c]:.4f}' for c in candidates[:3]])
    text += f'### {title}\n\n' + table(['Example', 'Training fit', 'Sieve', 'Landweber', 'PMMR'], rows)
text += '[Complete weights at both outcome times](conditional-moment-ensemble-weights.csv).\n\n'
text += '## Estimates and mean EIF components\n\nThe sequential mean is centered at the true parameter. Adding the signed bridge/adjoint mean gives estimate minus truth.\n\n'
text += table(['Example', 'Outcome time', 'Estimator', 'Truth', 'Estimate', 'SE', 'Mean sequential component', 'Mean bridge/adjoint component'], eifrows)
text += f'{misses} of the eight pointwise intervals exclude the truth. These two datasets do not estimate repeated-sample coverage. All fits and results are retained.\n\n'
text += '## Selected sieve penalties\n\nValues below are the final candidate refits within each outer training sample. All selected values are positive and interior. CV tables for the smaller ensemble-training fits are also retained in each case folder.\n\n'
text += table(['Example', 'Outcome time', 'Function', 'Training fit', 'Lambda', 'Grid values', 'Boundary failures'], penaltyrows)
text += '''## Checks and reproduction

Audits verify exact reuse of the saved data and split assignments, absence of outer validation observations from penalty training, the shared nested split callback, positive interior penalty selections, finite predictions and inverse-expit bridge values, raw cell U-statistics, and positive-semidefinite matrices before ensemble weight selection. The longitudinal applied bridge interval remains [-100, 100]; its use is recorded in [applied-bridge-truncation.csv](applied-bridge-truncation.csv). The figures show the raw predictions.

[Function and conditional-equation errors](function-and-equation-errors.csv) · [Selected penalties](sieve-selected-penalties.csv) · [Landweber class audit](landweber-function-class.csv) · [U-statistic audit](scorer-audit.csv) · [All predictions](all-function-predictions.csv) · [Conditional-equation values](all-equation-predictions.csv).

The full local record contains package archives, source snapshots, private installed libraries, scripts, fitted objects and validation logs. The public copy includes the package archives, scripts, numerical summaries, plots and audit records. A snapshot manifest records file checksums. The earlier population calculation remains applicable: β₁ is unique; β₂ is unique at R₂ = 1 and nonunique at R₂ = 0. The numerical-dose adjoint may be nonunique, so its defining conditional equation is also checked.

## GitHub validation

[Run {GITHUB_RUN}](https://github.com/idiazst/cmbridge-tests/actions/runs/{GITHUB_RUN}) passed every check using cmbridge commit **{PACKAGE_COMMIT}**. This includes all package unit tests, the original bridge and adjoint recovery examples, the ensemble selection examples, and the additional inverse-expit bridge and signed-adjoint recovery checks. The original specifications and acceptance tolerances are retained.

[All GitHub validation results](https://github.com/idiazst/cmbridge-tests/tree/main/results/run-{GITHUB_RUN}).

![Inverse-expit bridge and signed-adjoint recovery](github-run/linked-bridge/truth_vs_estimate.png)

![Original ensemble recovery checks](github-run/ensemble/truth_vs_estimate.png)

Local test logs retain all warnings. The earlier failed attempts retain their learner failures and convergence fixtures. Final n = 4,000 results, warnings and error records are reported in the case folders; no failed attempt is silently substituted with a successful result.
'''
text = text.replace('{GITHUB_RUN}', str(github['run_id'])).replace('{PACKAGE_COMMIT}', github['package_commit'])
(root/'REPORT.md').write_text(text)
(root/'README.md').write_text('# Sieve CV and quadratic Landweber\n\nSee [REPORT.md](REPORT.md) for the two n = 4,000 checks, plots, penalties, ensemble weights, EIF means and audits.\n')
(root/'validation/report-integrity.txt').write_text(f'Two datasets; {len(estimates)} estimates; zero fit errors.\nAll EIF sums agree with estimate minus truth to 1e-10.\n')
print(f'Wrote report with {len(estimates)} estimates.')
