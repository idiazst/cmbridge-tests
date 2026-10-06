"""Summarize the saved exact-score diagnostic, including unfavorable fits."""
import csv
import os
from pathlib import Path

root = Path(os.environ['EXACT_BRIDGE_OUTPUT'])


def read(path):
    return list(csv.DictReader(path.open())) if path.stat().st_size else []


def table(headers, rows):
    return '| ' + ' | '.join(headers) + ' |\n| ' + ' | '.join(['---'] * len(headers)) + ' |\n' + ''.join(
        '| ' + ' | '.join(map(str, row)) + ' |\n' for row in rows) + '\n'


names = {'binary_longitudinal': 'Binary treatment', 'discrete_dose': 'Numerical dose'}
summary = read(root / 'summary.csv')
comparisons = read(root / 'fixed-candidate-score-comparison.csv')
penalties, solvers, scores = [], [], []
for mechanism in names:
    folder = root / mechanism
    assert (folder / 'FINISHED.txt').exists()
    penalties += read(folder / 'penalty-audit.csv')
    solvers += read(folder / 'solver-audit.csv')
    scores += read(folder / 'score-audit.csv')
assert len(penalties) == 24 and all(row['interior'] == 'TRUE' and float(row['selected']) > 0 for row in penalties)
assert len(solvers) == 18 and all(row['tolerance_reached'] == 'TRUE' for row in solvers)
assert len(scores) == 6 and all(int(row['candidate_failures']) == 0 for row in scores)

text = '''# Exact Gaussian U-statistic selection: saved bridge checks

**This configuration is not adopted.** Exact kernel scoring does not resolve the bridge discrepancy: binary ensemble RMSE is 0.420815 and numerical-dose RMSE is 0.246350. These checks fit only the time-3 bridge; they do not produce new SDR/TMLE estimates or confidence-interval coverage results. The running original cloud study remains unchanged.

Both examples reuse the saved n = 4,000 datasets, seeds 5103006/5203006, with all three original outer training sets and learner labels. The sieve retains its full joint-category conditioning basis. Landweber uses its original degree-one conditioning and target classes, with the original positive conditioning-weight ridge. PMMR is unchanged. Sieve penalties use positive training-only CV and the existing boundary extension/failure rules. The existing exact Gaussian kernel option is used for penalty selection and ensemble selection, with training-only scaling and bandwidth; self-products are excluded and only the ensemble Gram is projected to positive semidefinite.

This configuration changes both sieve conditioning and score approximation relative to the baseline. Therefore a baseline comparison cannot isolate the effect of approximation. The separate fixed-candidate score comparison below changes only the scoring approximation, without refitting or changing any selected penalty.

cmbridge 0.3.0.9021 and modified lmtp 1.6.0.9021 are unchanged. All bridges retain inverse-expit. No generating functions enter fitting, no predictors are removed, no zero-penalty or known-function performance simulation is run, and no Gaussian mean update is used. Gaussian here describes the kernel used to score moments.

## Fitted bridge and complete conditional equations

'''
text += table(['Example', 'Candidate', 'Bridge RMSE', 'Conditional-equation RMSE'],
              [[names[row['mechanism']], row['candidate'], f"{float(row['bridge_rmse']):.6f}",
                f"{float(row['equation_rmse']):.6f}"] for row in summary])
text += 'Bridge RMSE is probability weighted in the group with intermediate and final visits observed, then averaged in squared units across the three training fits. Equation RMSE uses the full population conditional equation, including the missing-visit group. All plotted points appear on equal axes without cropping or clipping.\n\n'
for mechanism in names:
    text += f"![{names[mechanism]} ensemble]({mechanism}-ensemble.png)\n\n![{names[mechanism]} candidates]({mechanism}-all-candidates.png)\n\n"
text += '## Actual ensemble weights\n\n' + table(
    ['Example', 'Training fit', 'Sieve', 'Landweber', 'PMMR'],
    [[names[m], fold] + [f"{float(next(row['weight'] for row in read(root/m/'ensemble-weights.csv') if row['fold']==str(fold) and row['candidate']==c)):.6f}"
                        for c in ['sieve_md','landweber','pmmr']] for m in names for fold in range(1,4)])
text += '## Changing only the score approximation for fixed fitted candidates\n\n' + table(
    ['Example', 'Training fit', 'Maximum raw Gram change', 'Maximum ensemble weight change'],
    [[names[row['mechanism']], row['fold'], f"{float(row['maximum_raw_gram_change']):.8g}",
      f"{float(row['maximum_weight_change']):.8g}"] for row in comparisons])
text += 'The largest weight change is ' + f"{max(float(row['maximum_weight_change']) for row in comparisons):.6f}" + '. Thus the approximation has little effect on ensemble weights for these fixed candidate predictions. This does not measure changes in penalty selection that could occur if the candidates were refitted using the other score. It does not establish that the Gaussian score is sufficient to recover every fitted bridge accurately.\n\n'
text += '## Checks and measured runtime\n\n'
text += f'All 24 selected sieve CV records are positive, successful, finite interior minima; every failed trial remains recorded. All 18 final candidates reach their recorded solver tolerance. Direct full Gaussian ordered-pair formulas reproduce all six raw Grams; independent approximate-kernel formulas also reproduce the package U score. Projected Grams are PSD, final penalty person IDs exclude the outer validation people, and the shared learner labels are unchanged. No candidate fails and no fitted population value reaches the approved application bounds. Source files and model objects are preserved separately.\n\n'
text += table(['Example', 'Bridge fitting seconds across three training samples'],
              [[names[m], f"{sum(float(row['elapsed_seconds']) for row in scores if row['mechanism']==m):.3f}"] for m in names])
text += 'These are bridge-only fitting times, excluding strict nested sequential-response estimation; they cannot forecast complete-dataset runtimes. The complete-estimator conditioning-weight checks remain in [their separate report](../ensemble-stability-v8-weight-cv/REPORT.md).\n\n[Actual errors](all-bridge-errors.csv) · [All weights](all-ensemble-weights.csv) · [Fixed-candidate score comparison](fixed-candidate-score-comparison.csv) · [Configuration](EXPERIMENT.md).\n'
(root/'REPORT.md').write_text(text)
print('Reported all six saved bridge fits and the unfavorable exact-score comparison.')
