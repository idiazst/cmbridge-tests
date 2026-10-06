"""Reports from saved diagnostic tables. This does not fit statistical models."""
import csv
import os
from pathlib import Path


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def num(row, name):
    return float(row[name])


def write(root, text):
    (root / 'REPORT.md').write_text(text)
    print(root / 'REPORT.md')


base = Path(os.environ['DIAGNOSTIC_BASE'])
root = base / 'bridge-equation-identification-v1'
rows = read(root / 'bridge-equation-identification.csv')
table = ['| Treatment | Outcome time | Independent fitted equations | Independent full conditional equations | Full equation RMSE of the alternative | Largest fitted-moment error |',
         '| --- | --- | --- | --- | --- | --- |']
for r in rows:
    treatment = 'Binary' if r['mechanism'] == 'binary_longitudinal' else 'Numerical dose'
    table.append(f"| {treatment} | {r['outcome_time']} | {r['polynomial_operator_rank']} | "
                 f"{r['full_conditional_operator_rank']} | {num(r,'alternative_full_equation_rmse'):.6f} | "
                 f"{num(r,'alternative_polynomial_moment_max_error'):.2e} |")
write(root, f'''# Which bridge equations are enforced by the current basis?

The saturated sieve **function class contains a valid bridge**, but its current additive conditioning basis does not enforce all the paper's conditional equations. These are separate requirements. This audit demonstrates the distinction directly on the full discrete population in both approved mechanisms.

The paper requires the conditional mean of the measured bridge, given the observed history and treatment, to equal one. The present degree-one conditioning basis checks averages against the intercept and individual predictors. With a flexible joint-category bridge, satisfying those averages can leave errors within particular combinations of predictors.

For each outcome time, the audit constructed a different bridge in the complete joint-category inverse-expit class. Every value remains strictly above one. It satisfies every fitted additive moment to numerical precision, yet violates the full conditional equation. No datasets or nuisance estimators were fitted; this is an algebraic calculation, not a known-function or zero-penalty simulation.

{chr(10).join(table)}

The reference bridge satisfies the full equations with maximum error below {max(num(r,'reference_full_equation_max_error') for r in rows):.2e}. The alternative functions have values between {min(num(r,'alternative_minimum') for r in rows):.3f} and {max(num(r,'alternative_maximum') for r in rows):.3f}, and their representation errors in the sieve class are zero to recorded precision. Thus the discrepancy is not caused by clipping or lack of a valid bridge in the DGP.

This does **not** show that an actual fit selects the constructed alternative, or that the whole ensemble is inconsistent. The fixed main-effect Landweber class is smaller: the number of independent linearized fitted equations equals the number of independent target-basis coefficients at the valid bridge in all four cases. That is a local check around the valid solution, not a proof of global uniqueness. PMMR was not given this identification audit.

Adding a full joint-category conditioning basis removes this particular population limitation. However, the [separate estimated-model comparison](../bridge-conditioning-check-v1/REPORT.md) made finite-sample predictions worse in an existing n=4000 dataset. Therefore that change has not been substituted into the running study. More conditional equations also require estimating more conditional averages from the same observations; population correctness alone does not establish better finite-sample performance.

The frozen study and its statistical sources are preserved. Any future fitting change requires a separate version and new checks. [Numerical audit](bridge-equation-identification.csv), [all conditional errors of the constructed alternatives](alternative-conditional-equations.csv), and [reproducible algebraic script](../../../simulations/observed_history_corrected/audit_bridge_defining_equations.R) retain the evidence.
''')

root = base / 'bridge-conditioning-check-v1'
rows = read(root / 'population-bridge-checks.csv')
weights = read(root / 'ensemble-weights.csv')
audit = read(root / 'saved-fit-audit.csv')
penalties = read(root / 'selected-penalty-audit.csv')
names = {'ensemble': 'Ensemble', 'sieve_md': 'Sieve minimum distance', 'landweber': 'Landweber', 'pmmr': 'PMMR'}
table = ['| Conditioning basis | Candidate | Bridge RMSE, R2 = 1 | Conditional-equation RMSE | Largest conditional-equation error |',
         '| --- | --- | --- | --- | --- |']
for r in rows:
    label = 'Original additive' if r['basis'] == 'polynomial' else 'Joint categories'
    table.append(f"| {label} | {names[r['candidate']]} | {num(r,'bridge_rmse_R2_1'):.6f} | "
                 f"{num(r,'defining_equation_rmse'):.6f} | {num(r,'defining_equation_max_error'):.6f} |")
wt = ['| Conditioning basis | Sieve minimum distance | Landweber | PMMR |',
      '| --- | --- | --- | --- |']
for basis in ['polynomial', 'joint-categories']:
    values = {r['candidate']: num(r, 'weight') for r in weights if r['basis'] == basis}
    wt.append(f"| {'Original additive' if basis == 'polynomial' else 'Joint categories'} | "
              f"{values['sieve_md']:.6f} | {values['landweber']:.6f} | {values['pmmr']:.6f} |")
diff = next(num(r, 'baseline_prediction_max_difference') for r in audit if r['basis'] == 'polynomial')
write(root, f'''# Conditioning-basis comparison in one saved n=4000 dataset

**The joint-category conditioning basis made the ensemble and sieve bridge predictions worse in this check. It has not been adopted as a repair.**

This uses the binary-treatment dataset from the completed Mac runtime check: sample size 4,000, replication 2, seed 5103002. Only the first outer training sample and the bridge for the outcome at time 3 are fitted. It is a diagnostic of bridge estimation, not another replication of the complete estimator or a coverage study.

Both fits use identical data, learner sample assignments, positive penalty cross-validation, Gaussian U-statistic scoring, and ensemble optimization. The comparison changes the conditioning basis to joint categories for the sieve and Landweber candidates; PMMR is unchanged. The response and predictors contain only estimated-model inputs. Population true functions are used after fitting to evaluate predictions and conditional equations.

The original fit reproduces the already saved ensemble predictions with maximum difference {diff:.2e}. All three final candidate fits satisfy their recorded convergence criteria. There are zero final candidate failures. All eight checked sieve penalty selections are positive, finite, successful and interior, including selections inside the learner training samples. No application-bound hits occurred.

{chr(10).join(table)}

The bridge RMSE column uses the population solution where the intermediate and final visits are observed. The conditional-equation columns use the full population, including missing intermediate visits. They measure different quantities.

![All bridge candidates](figures/all-bridge-candidates.png)

Every point is shown with equal axes, using the same range across all panels. Point size reflects population probability. The wider range is required by the joint-category sieve estimates; points have not been cropped to make the comparison look better.

{chr(10).join(wt)}

The original fit took 9.934 seconds and the joint-category fit 16.852 seconds. These are bridge-only timings, not full-dataset simulation times. The joint-category sieve receives greater ensemble weight despite its larger population error. Ensemble selection uses estimated out-of-sample scores, not the population truth, so the weights do not establish which candidate is most accurate in this dataset.

The [population equation audit](../bridge-equation-identification-v1/REPORT.md) establishes a limitation of the original sieve conditioning basis. This estimated-model comparison shows that simply saturating that basis does not solve the finite-sample problem. It does not identify a proven package repair, and the current cloud study remains unchanged.

Both original diagnostic invocations saved successful fits and all metrics, then failed at the final reporting expression because of a missing closing parenthesis. The attempt logs are preserved. A corrected retry was prevented from overwriting the saved models by the existing guard. A separate read-only report script completed the audit and graphs without refitting. The reporting error is distinct from a statistical fit failure.

[Prediction errors](population-bridge-checks.csv), [weights](ensemble-weights.csv), [selected penalties](selected-penalty-audit.csv), [solver records](solver-audit.csv), and [saved-fit verification](saved-fit-audit.csv) retain the calculations. [Fit script](../../../simulations/observed_history_corrected/check_bridge_conditioning_basis.R) and [read-only report script](../../../simulations/observed_history_corrected/report_bridge_conditioning_comparison.R) reproduce the procedure. The fitted objects remain local and are excluded from the published results.
''')

root = base / 'full-ensemble-cloud-v1' / 'early-equation-diagnostics'
rows = read(root / 'equation-diagnostics-by-dataset.csv')
time3 = [r for r in rows if r['horizon'] == '2']
flags = sorted(time3, key=lambda r: -abs(num(r, 'bridge_remainder')))[:2]
pop = read(root / 'population-components-by-dataset.csv')
subgroups = read(root / 'time3-R2-contributions-by-dataset.csv')
penalties = read(root / 'saved-penalty-minima.csv')
weights = read(base / 'full-ensemble-cloud-v1' / 'early-summary' / 'ensemble_weights.csv')
assert len(time3) == 24 and len(penalties) == 720
assert all(r['selected_is_exact_minimum'] == 'TRUE' and r['selected_is_interior'] == 'TRUE'
           and r['selected_at_successful_grid_edge'] == 'FALSE' for r in penalties)
flag_tables = []
for flag in flags:
    match = lambda r: (r['n'], r['replicate']) == (flag['n'], flag['replicate'])
    components = [r for r in pop if match(r)]
    table = ['| Estimator | Sequential population component | Bridge/adjoint population component | Sum |',
             '| --- | --- | --- | --- |']
    for r in components:
        table.append(f"| {r['estimator'].upper()} | {num(r,'sequential_component'):+.6f} | "
                     f"{num(r,'bridge_component'):+.6f} | {num(r,'population_error'):+.6f} |")
    table += ['', '| Intermediate visit | Bridge/adjoint contribution | Bridge equation RMSE | Adjoint equation RMSE |',
              '| --- | --- | --- | --- |']
    for r in subgroups:
        if match(r):
            table.append(f"| {'Missing (R2 = 0)' if r['R2'] == '0' else 'Observed (R2 = 1)'} | "
                         f"{num(r,'bridge_remainder_contribution'):+.6f} | {num(r,'bridge_equation_rmse'):.6f} | "
                         f"{num(r,'adjoint_equation_rmse'):.6f} |")
    table += ['', '| Training sample | Sieve adjoint weight | Landweber adjoint weight | PMMR adjoint weight |',
              '| --- | --- | --- | --- |']
    for fold in ['1', '2', '3']:
        selected = [r for r in weights if match(r) and r['kind'] == 'adjoint'
                    and r['horizon_or_depth'] == '2' and r['fold'] == fold]
        values = {r['candidate']: num(r, 'weight') for r in selected}
        assert set(values) == {'sieve_md', 'landweber', 'pmmr'}
        table.append(f"| {fold} | {values['sieve_md']:.6f} | {values['landweber']:.6f} | {values['pmmr']:.6f} |")
    figure = f"figures/binary-n{flag['n']}-r{int(flag['replicate']):03d}-adjoint-equations.png"
    total = num(flag, 'bridge_remainder')
    unseen = num(flag, 'remainder_from_unseen_adjoint_target_combinations')
    few = num(flag, 'remainder_from_adjoint_target_combinations_with_at_most_two_rows')
    flag_tables.append(f'''### Binary treatment, n = {flag['n']}, replication {flag['replicate']}

{chr(10).join(table)}

The missing and observed contributions add to {total:+.6f}; the large discrepancy is mostly in observations with the intermediate visit observed. Measured target combinations absent from the adjoint training sample contribute {unseen:+.6f}. Combinations with at most two measured training rows, including absent combinations, contribute {few:+.6f}. Other combinations contribute {total-few:+.6f}. Signed contributions can offset one another; both absent and present combinations affect the remainder. Counts alone do not determine the source of their fitted errors.

The actual adjoint conditional-equation RMSE is {num(flag,'adjoint_equation_rmse'):.6f}. The minimum attainable within the observed sieve dictionary is {num(flag,'observed_dictionary_adjoint_class_minimum_equation_rmse'):.6f}. That comparison concerns the sieve dictionary; it is not a lower bound for Landweber, PMMR, or the whole ensemble.

![Adjoint conditional equations]({figure})
''')
ties = sum(int(r['exact_minimum_count']) > 1 for r in penalties)
boundary_ties = sum(r['minimum_also_at_lower_boundary'] == 'TRUE'
                    or r['minimum_also_at_upper_boundary'] == 'TRUE' for r in penalties)
write(root, f'''# Equation checks for the 24 saved cloud datasets

**These are partial results from a running study, not the final 200-replication comparisons.** They contain nine binary-treatment datasets at sample size 500, nine at 1,000, and six at 4,000. No numerical-dose datasets are saved in this snapshot. Completion by the early time limit may depend on data and fitting difficulty, so this subset should not be treated as a random sample of planned replications.

All diagnostics evaluate the actual saved estimated functions. The original data, outer sample assignments and learner assignments were regenerated and matched exactly. No nuisance models were refitted. The full discrete population permits exact evaluation of conditional equations and fixed-function population contributions, without a Monte Carlo approximation to the true parameter.

The two datasets below were selected by the largest absolute bridge/adjoint population contribution at outcome time 3 among all 24 saved datasets. The population components are **different from the sample means of EIF components** in the simulation report. They hold estimated functions fixed and average over the DGP. TMLE targeting uses validation outcomes, so this calculation is a diagnostic and not a literal conditional expectation given training data alone.

{chr(10).join(flag_tables)}

## Other saved datasets

![Population components](figures/population-components.png)

The two coordinates sum to the population estimation error. A point on the dashed line has zero sum. The picture neither estimates final coverage nor establishes a convergence rate.

![Adjoint equation and class errors](figures/adjoint-equations-and-dictionaries.png)

The complete joint-category adjoint class contains a solution of the paper's equations over the whole support. The implemented sieve dictionary uses only combinations in measured training rows, with mean continuation for an unseen combination. That smaller sample-dependent class can exclude every exact solution in a finite sample. The class audit permits any solution; it does not require agreement with a single chosen reference. It is population algebra, not a prohibited known-function or zero-penalty performance study.

Large actual errors remain beyond this dictionary limitation. Treatment-ratio estimation, finite-sample equation estimation and selection also require consideration. The fitted-ratio columns use final outer ratio predictions; the actual adjoint training labels used nested out-of-sample predictions. Those two sets of predictions are not assumed equal. Signed conditional errors add, whereas RMSEs do not.

## Penalty checks

All 720 selected records are positive, finite, successful, exact minima and interior to their successful grids. None has a failed immediate neighbour. The audit therefore did not find selection of an invalid trial or an actual boundary penalty choice.

There are {ties} records with more than one exactly minimizing penalty, including {boundary_ties} whose set of tied minima also contains a boundary. An interior selected minimizer is not a uniquely isolated optimum. This distinction is retained rather than silently claiming every CV curve has a unique minimum. The time-3 tied adjoint choice in binary sample size 500, replication 9, has zero ensemble weight; it does not explain the two flagged datasets above. No boundary rule or fit source was changed for this audit.

The separate [bridge equation audit](../../bridge-equation-identification-v1/REPORT.md) found that the saturated sieve's additive conditioning moments do not by themselves enforce the full paper equation. The [estimated conditioning-basis comparison](../../bridge-conditioning-check-v1/REPORT.md) tested a straightforward change on an existing n=4000 training sample and found worse errors. Both findings are preserved; neither has been used to alter or discard the running study's results.

[All dataset equations](equation-diagnostics-by-dataset.csv), [fold equations](equation-diagnostics-by-fold.csv), [population components](population-components-by-dataset.csv), [intermediate-visit contributions](time3-R2-contributions-by-dataset.csv), [penalty audit](saved-penalty-minima.csv), and [all conditional adjoint values](adjoint-conditional-equations.csv) retain the numerical evidence. Actual sample EIF means and provisional coverage remain in the [separate running-study report](../early-summary/REPORT.md).
''')
