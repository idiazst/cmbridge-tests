"""Report exact equation audits of saved fits; no nuisance models are fitted."""
import csv
import os
from pathlib import Path

root = Path(os.environ['DIAGNOSTIC_OUTPUT'])
rows = list(csv.DictReader((root / 'equation-diagnostics-by-dataset.csv').open()))
time3 = [r for r in rows if r['horizon'] == '2']
flag = next(r for r in time3 if r['mechanism'] == 'binary_longitudinal'
            and r['n'] == '1000' and r['replicate'] == '1')
def f(row, column):
    return float(row[column])

table = ['| Treatment | n | Replication | Bridge/adjoint remainder | Adjoint equation RMSE | Smallest equation RMSE using the observed dictionary | Measured probability with unseen adjoint target values |',
         '| --- | --- | --- | --- | --- | --- | --- |']
for r in time3:
    name = 'Binary' if r['mechanism'] == 'binary_longitudinal' else 'Numerical dose'
    table.append(f"| {name} | {r['n']} | {r['replicate']} | {f(r,'bridge_remainder'):.6f} | "
                 f"{f(r,'adjoint_equation_rmse'):.4f} | "
                 f"{f(r,'observed_dictionary_adjoint_class_minimum_equation_rmse'):.4f} | "
                 f"{100*f(r,'measured_probability_in_unseen_adjoint_target_combinations'):.2f}% |")
full_max = max(f(r,'full_support_adjoint_class_maximum_equation_error') for r in rows)
unseen = f(flag, 'remainder_from_unseen_adjoint_target_combinations')
few = f(flag, 'remainder_from_adjoint_target_combinations_with_at_most_two_rows')
total = f(flag, 'bridge_remainder')
report = f'''# Defining-equation checks on the ten saved Mac datasets

These are diagnostics of the initial local timing run, not the cloud study's final results. No nuisance models were refitted. The original data were regenerated from each saved seed, and outer and learner sample assignments were checked for exact equality. The saved fitted functions were evaluated over the complete discrete population. The true functions enter only these diagnostic calculations.

## Function classes

The fixed main-effect inverse-expit bridge classes contain a valid solution in both mechanisms, as established by the separate full-support bridge audit. This remains true when a joint category is absent from a training sample.

The saturated adjoint class on the **whole support** also contains a valid solution: the largest defining-equation error in its algebraic representation audit is {full_max:.2e}. However, the fitted sieve learns its category dictionary from **measured training rows** and assigns the mean of the fitted category values to an unseen combination. This produces a smaller, sample-dependent class. In the saved datasets, that class does not contain an exact population solution at outcome time 3.

This qualifies the earlier report's statement that every fitted class contains the required function. The complete saturated class does; the realized observed-category dictionary need not. This is a finite-sample restriction, and the fraction of missing categories declines with sample size. The audit permits any adjoint solution, so its conclusion does not rely on uniqueness or equality to one chosen reference solution. It solves population linear equations algebraically; it is not a zero-penalty simulation or a performance experiment using known nuisance functions.

The table reports the smallest probability-weighted conditional-equation RMSE attainable by that observed dictionary, including its mean continuation. It also reports the error of the actual ensemble. Each dataset row averages its three outer training samples. These are different from sample means of EIF components.

For the flagged binary n=1000 dataset, allowing the unseen-value intercept to vary freely still leaves minimum equation RMSE {f(flag,'observed_dictionary_with_free_intercept_minimum_equation_rmse'):.4f}. The failure is therefore not solely the rule that fixes that intercept at the category mean. The retained operator condition numbers are far below the inverse singular-value cutoff used by the audit, so discarding numerical near-zero directions does not explain the residual.

{chr(10).join(table)}

![Observed dictionary class audit](figures/observed-dictionary-adjoint-class.png)

## The binary n=1000 dataset with a large remainder

For replication 1 at outcome time 3, the bridge/adjoint remainder is {total:+.6f}. Combinations absent from the measured adjoint training sample contribute {unseen:+.6f}; combinations with at most two measured training observations contribute {few:+.6f}. The latter includes the absent combinations. The remaining contribution is {total-few:+.6f}. Thus absent or infrequent combinations explain part of the error, and are not its sole cause.

The actual adjoint equation RMSE is {f(flag,'adjoint_equation_rmse'):.4f}, compared with {f(flag,'observed_dictionary_adjoint_class_minimum_equation_rmse'):.4f} for the best function in the observed dictionary. The ensemble places all adjoint weight on Landweber in two outer training samples and all weight on the sieve in the third. The gap between the attainable and fitted errors requires examination alongside selection and the estimated treatment ratios; expanding a dictionary alone is not evidence that the fitted estimator would improve.

![Actual adjoint equations](figures/binary-n1000-r001-adjoint-equations.png)

The adjoint loading also depends on estimated treatment ratios. Replacing the true loading with predictions from the final outer ratio fits gives conditional-equation RMSE {f(flag,'adjoint_error_using_fitted_ratios_rmse'):.4f}; the conditional difference attributable to those ratio predictions has RMSE {f(flag,'adjoint_difference_from_ratio_estimation_rmse'):.4f}. The signed conditional errors add exactly, but their RMSEs do not. Actual adjoint training labels use nested out-of-sample ratio predictions, so this diagnostic does not assert that final outer predictions equal the training labels.

These diagnostics establish neither asymptotic failure nor confidence-interval coverage. The frozen cloud study retains the approved fitting settings and will measure their actual finite-sample behavior. No source changes have been made during its run. Any later repair must be evaluated in a separate version, preserving these results.

[Fold-level values](equation-diagnostics-by-fold.csv), [dataset values](equation-diagnostics-by-dataset.csv), [adjoint conditional equations](adjoint-conditional-equations.csv), and [bridge conditional equations](bridge-conditional-equations.csv) retain all numerical calculations.
'''
(root / 'REPORT.md').write_text(report)
print(root / 'REPORT.md')
