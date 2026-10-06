"""Write the component-validation report from completed numerical records."""
import csv
import hashlib
from pathlib import Path

out = Path("results/observed-history-corrected/dose-cmbridge-debug")


def read(name):
    with (out / name).open() as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    return float(row[key])


summary = read("fit-diagnostics.csv")
reference = next(row for row in summary if row["source"] == "saved check" and
                 row["time"] == "2" and row["fold"] == "2")
large = next(row for row in summary if row["source"] == "large check")
ranking = read("individual-fit-ranking.csv")
tuning = read("large-penalty-selection.csv")
counts = [row for row in read("large-training-counts.csv") if row["time"] == "2"]
training_n = sum(int(row["training"]) for row in counts)
selected = [row for row in tuning if row["selected"] == "TRUE" and
            row["kind"] in ("beta", "adjoint")]
equations = read("true-equation-and-parameter-audit.csv")
operators = read("full-history-operator-audit.csv")
alternative = read("alternative-adjoint-check.csv")[0]
constant = read("constant-term-criterion-diagnostic.csv")[0]
sample_equations = read("sample-equation-summary.csv")
overflow = read("integer-overflow-audit.csv")[0]
large_bridge_errors = read("bridge-large-errors.csv")
large_adjoint_errors = read("adjoint-equation-large-errors.csv")
classes = read("function-class-membership.csv")
point_summary = read("large-bridge-point-summary.csv")[0]
shrinkage = read("positive-penalty-shrinkage-comparison.csv")
regression_classes = read("sequential-regression-class-membership.csv")
classifier_classes = read("treatment-classifier-class-membership.csv")

lines = [
    "# Debugging cmbridge in the numerical dose example", "",
    "This check uses the same two-time, four-level numerical dose mechanism as "
    "the earlier checks. It uses the second n = 4,000 dataset, seed 5203002, "
    "and its saved estimator and learner assignments. The true mean at outcome "
    "time 3 is 0.871138516760. R₁ and R₂ are always 1; R₃ is informative. "
    "Treatment is random conditional on Rₜ = 1, and its coded no-visit branch "
    "is the deterministic continuation rule (0, 0). That branch is not reached "
    "at either treatment time in this mechanism.", "",
    "## Time point and fit with the largest error", "",
    "The cmbridge contribution is at treatment time 2, for outcome time 3. "
    "At time 1 the intermediate outcome is completely observed: β₁ is exactly "
    "1, its correction is zero, and no adjoint is fitted.", "",
    "The largest error in the selected dataset occurs in outer training sample "
    "2. Both fitted functions use the saturated L1 conditional-moment learner, "
    "with a positive penalty selected by cross-validation. Its single-candidate "
    "ensemble weight is one.", "",
    "To compare individual fits on the same parameter scale, evaluate the two "
    "representations from the paper using the known generating probabilities "
    "for evaluation only:", "",
    "$$", r"E\{R_3\Omega_2Y_3\widehat\beta_2(H_2,C_3)\},", "$$", "",
    "$$", r"E\{\widehat\lambda_2(A_2,H_2)\}.", "$$", "",
    "The errors below are these implied means minus the true parameter, "
    "calculated over the entire generating distribution. They are not biases "
    "averaged over repeated datasets.", "",
    "| Function | Training sample | Error in implied parameter mean |",
    "| --- | ---: | ---: |",
]
for row in ranking:
    if row["time"] == "2":
        name = "β₂" if row["function_name"] == "beta" else "λ₂"
        lines.append(f"| {name} | {row['fold']} | {number(row, 'representation_error'):.6f} |")
lines += [
    "", "λ₂ in sample 2 has the largest absolute error by this criterion; "
    "β₂ in the same sample has a very similar error. Their errors interact in "
    "the one-step correction, so that contribution cannot be uniquely "
    "attributed to just one of the two functions.", "",
    f"For sample 2, the bridge representation error is "
    f"{number(reference, 'beta_representation_error'):.6f}. The bridge correction "
    f"adds {number(reference, 'bridge_correction'):.6f}, leaving a bridge "
    f"contribution to the remainder of {number(reference, 'bridge_remainder'):.6f}. "
    "The separate sequential-regression contribution is not included in this number.", "",
    "The selected n = 4,000 nuisance fit was regenerated from its original "
    "seed and saved assignments. Its β and λ predictions matched the saved "
    "predictions exactly. The direct correction also agrees with the product "
    "identity in the second-order remainder.", "",
    "## Large-dataset check", "",
    f"A separate sample of {int(large['n']):,} people was generated with seed "
    f"{large['seed']}, and training sample 2 was fitted using {training_n:,} "
    "people. The generating mechanism, "
    "complete observed histories, learner, bounds, penalty candidates, and "
    "splitting construction remain the same. The bridge and adjoint are fitted "
    "with the public cmbridge ensemble routines used by the original nuisance "
    "fitting procedure. The treatment ratios used to construct the adjoint "
    "response are estimated with SuperLearner; generating probabilities are "
    "not supplied as learner responses or predictors.", "",
    "The positive penalty scales remain 0.001, 0.01, and 0.1, divided by the "
    "square root of the number of training people in each fit. Cross-validation "
    "selects the scale within the corresponding training sample.", "",
    "Bridge prediction and adjoint equation RMSEs use generating probabilities "
    "among measured outcomes. Bridge equation RMSE uses the natural distribution "
    "of treatment and history. The bridge contribution to the remainder uses "
    "the true treatment ratios for evaluation, isolating this component; it "
    "is not the bias of a final SDR or TMLE fit.", "",
    "| Quantity | n = 4,000, sample 2 | Large dataset, sample 2 |",
    "| --- | ---: | ---: |",
]
metrics = [
    ("Mean bridge prediction error", "beta_mean_error"),
    ("Bridge prediction RMSE", "beta_rmse"),
    ("Bridge conditional-equation RMSE", "bridge_equation_rmse"),
    ("Adjoint conditional-equation RMSE", "adjoint_equation_rmse"),
    ("Error in mean implied by bridge", "beta_representation_error"),
    ("Error in mean implied by adjoint", "lambda_representation_error"),
    ("Bridge contribution to remainder", "bridge_remainder"),
]
for name, key in metrics:
    lines.append(f"| {name} | {number(reference, key):.6f} | {number(large, key):.6f} |")
lines += [
    "", f"The bridge prediction RMSE falls from {number(reference, 'beta_rmse'):.3f} "
    f"to {number(large, 'beta_rmse'):.3f}, and the adjoint equation RMSE falls from "
    f"{number(reference, 'adjoint_equation_rmse'):.3f} to "
    f"{number(large, 'adjoint_equation_rmse'):.3f}. The bridge contribution to the "
    f"one-step remainder falls from {number(reference, 'bridge_remainder'):.6f} "
    f"to {number(large, 'bridge_remainder'):.6f}. This is substantial improvement "
    "using the same learners, with the positive penalties "
    "chosen by cross-validation.", "",
    "The largest errors have not disappeared. The largest bridge prediction "
    f"error is {abs(number(large_bridge_errors[0], 'error')):.3f}, in a combination "
    f"with {large_bridge_errors[0]['training_measured']} measured training "
    "observations. The largest adjoint equation error is "
    f"{abs(number(large_adjoint_errors[0], 'error')):.3f}, in a combination with "
    f"{large_adjoint_errors[0]['training_measured']} measured observations and "
    f"generating probability {number(large_adjoint_errors[0], 'probability'):.2e} "
    "among measured outcomes. Thus even this large overall sample contains "
    "very small groups for some complete histories. The large-data adjoint "
    "error includes error from its estimated treatment-ratio responses; this "
    "check does not isolate the contribution of that error.", "",
]
lines += ["", "Selected large-data penalty scales:", ""]
for row in selected:
    name = "bridge" if row["kind"] == "beta" else "adjoint"
    lines.append(f"- {name}: {row['scale']}.")
lines += [
    "", f"The large training sample contains {len(counts)} current-treatment "
    f"and history combinations. Its smallest combination has "
    f"{min(int(row['training']) for row in counts):,} "
    f"{'person' if min(int(row['training']) for row in counts) == 1 else 'people'}, and "
    f"{sum(int(row['measured']) == 0 for row in counts)} combinations have no "
    "measured final outcome.", "",
    "## Estimate against true value", "",
    "For β₂, the exact true bridge is the inverse of the generating final "
    "measurement probability. The horizontal axis is that true value and the "
    "vertical axis is the fitted value. The blue line is y = x. Each point "
    "represents a complete combination of H₂ and C₃; its area reflects the "
    "generating probability among measured outcomes. The plots evaluate all "
    "such combinations, including ones absent from a smaller training sample.", "",
    "![Bridge estimate against truth, large dataset](figures/bridge-y-equals-x-large.png)", "",
    "![Bridge comparison with the selected n = 4,000 fit](figures/bridge-y-equals-x-comparison.png)", "",
    "For λ₂, the DGP admits multiple valid solutions. A pointwise comparison "
    "with one arbitrarily chosen solution can therefore mislabel a valid fit. "
    "The next plot compares its conditional equation, which has a unique "
    "target from the generating distribution:", "",
    "$$", r"E\{\widehat\lambda_2(A_2,H_2)\mid R_3=1,H_2,C_3\}",
    r"\quad\text{against}\quad\phi_2(H_2,C_3).", "$$", "",
    "![Adjoint conditional equation against truth](figures/adjoint-equation-y-equals-x-large.png)", "",
    "![Adjoint conditional-equation comparison](figures/adjoint-equation-y-equals-x-comparison.png)", "",
    "These are component checks, not a new repeated-sample coverage study. "
    "The million-person run does not fit the final SDR or TMLE estimator or "
    "establish confidence-interval coverage.", "",
    "## Reproduction and numerical records", "",
    "[Individual-fit ranking](individual-fit-ranking.csv), "
    "[all fit diagnostics](fit-diagnostics.csv), "
    "[bridge predictions](bridge-predictions.csv), "
    "[adjoint conditional-equation predictions](adjoint-equation-predictions.csv), "
    "[large-data penalty selection](large-penalty-selection.csv), "
    "[large training counts](large-training-counts.csv), "
    "[saved-fit reproduction](reference-reproduction.rds), "
    "[large fitted functions and coefficients](large-fit.rds), "
    "[execution script](run-debug-v2.R), and [execution log](run-v2.log).", "",
    "The earlier run.log records a diagnostic-script setup error before the "
    "large learner fit; run-debug-v2.R corrected that setup error. The fits "
    "use the frozen package versions identified by the input checksums and "
    "archives. A scoring-arithmetic correction made after the large fit is "
    "documented below; it leaves these fitted functions unchanged.", "",
]
lines += [
    "## Empirical equation values in the n = 4,000 sample", "",
    "The paper's equations can also be assessed using the observed sample. "
    "These calculations use the same n = 4,000 dataset, seed 5203002. Each "
    "person's bridge, adjoint, and treatment-ratio predictions come from the "
    "saved fit whose estimator training sample excludes that person. No new "
    "learner was fitted for this calculation.", "",
    "For the bridge, within each observed treatment-and-history combination, "
    "calculate the empirical left side and compare it with 1:", "",
    "$$", r"\widehat E\{R_3\widehat\beta_2(H_2,C_3)\mid A_2,H_2\}=1.", "$$", "",
    "For the adjoint, within each observed complete-health-and-history "
    "combination among measured outcomes, compare the empirical left side "
    "with the estimated right side:", "",
    "$$", r"\widehat E\{\widehat\lambda_2(A_2,H_2)\mid R_3=1,H_2,C_3\}",
    r"\quad\text{against}\quad", 
    r"\widehat E\{\widehat\Omega_2Y_3\mid R_3=1,H_2,C_3\}.", "$$", "",
    "The adjoint right side uses estimated treatment ratios from outside "
    "each person's estimator fold. It does not use generating probabilities. "
    "Its unknown right side must be estimated when checking this equation "
    "in actual data.", "",
    "The following averages weight each observed conditioning group by its "
    "number of evaluation observations. The residual is left side minus "
    "right side. The RMSE first forms a residual within each conditioning "
    "group and then averages its square using these weights.", "",
    "| Evaluation observations | Equation | Mean left side | Mean right side | Mean residual | Conditional-equation RMSE |",
    "| --- | --- | ---: | ---: | ---: | ---: |",
]
for row in sample_equations:
    if row["fold"] in ("0", "2"):
        scope = "All 4,000, cross-fitted" if row["fold"] == "0" else "Sample 2 fit: held-out observations"
        name = "Bridge" if row["equation"] == "bridge" else "Adjoint"
        lines.append(f"| {scope} | {name} | {number(row, 'mean_lhs'):.6f} | "
                     f"{number(row, 'mean_rhs'):.6f} | {number(row, 'mean_residual'):.6f} | "
                     f"{number(row, 'conditional_equation_rmse'):.6f} |")
lines += [
    "", "The combined sample has 409 observed treatment-and-history "
    "groups for the bridge and 203 observed complete-health-and-history "
    "groups among 1,730 measured outcomes for the adjoint. Their median "
    "evaluation group sizes are 6 and 5; 68 bridge groups and 33 adjoint "
    "groups contain just one observation. The selected sample-2 fit is "
    "evaluated on 1,333 held-out people, including 563 measured outcomes.", "",
    "These are empirical conditional means, so finite-sample variation is "
    "part of the reported errors, particularly in small groups. They are "
    "distinct from the exact generating-distribution equation errors above. "
    "The combined bridge mean below 1 shows a systematic calibration "
    "shortfall in this sample. This table alone is not a formal statistical "
    "test or a decomposition of the adjoint error into treatment-ratio "
    "and adjoint-fitting errors.", "",
    "At time 1, R₂ and β₁ both equal 1, so every empirical bridge equation "
    "has left side exactly 1. No adjoint is fitted at time 1 in this estimator.", "",
    "Every observed conditioning-group equation value is saved in "
    "[bridge equation values](sample-bridge-equation-values.csv) and "
    "[adjoint equation values](sample-adjoint-equation-values.csv). "
    "[Summary by estimator fold](sample-equation-summary.csv), "
    "[observation-level predictions](sample-equation-observations.csv), "
    "[calculation script](evaluate_dose_sample_equations.R), and "
    "[calculation log](sample-equations.log) are also included.", "",
]
lines += [
    "## Paper equations and existence of solutions", "",
    "Section 3 requires the bridge equation (1) and adjoint equation (4):", "",
    "$$", r"E\{R_{t+1}\beta_t(H_t,C_{t+1})\mid A_t,H_t\}=1,", "$$", "",
    "$$", r"E\{\lambda_t(A_t,H_t)\mid R_{t+1}=1,H_t,C_{t+1}\}",
    r"=\phi_t(H_t,C_{t+1}).", "$$", "",
    "Theorem 2 states that the parameter is invariant across square-integrable "
    "solutions to these equations. Equality to one chosen function is therefore "
    "not a general validity requirement. The conditional-equation checks and "
    "parameter representations are the relevant assessments. Appendix A.6, "
    "Assumption 16 additionally selects a minimum-norm bridge path and requires "
    "its differentiability for the gradient argument.", "",
    "An exact audit used the full observed histories from this DGP, including "
    "every history with positive natural or policy probability. There are two "
    "baseline histories at time 1 and 64 histories at time 2. Within each, "
    "there are eight treatment values and four complete-health values. The "
    "bridge equation has an 8-by-4 matrix of rank 4, so its solution is unique "
    "for this DGP. The adjoint equation has a 4-by-8 matrix of rank 4, so it "
    "admits solutions with four dimensions of nonuniqueness.", "",
    "The bridge is the inverse of the positive generating measurement "
    "probability. The adjoint target is in the range of its operator at every "
    "reachable history. The exact calculations give:", "",
    "| Time | True bridge equation maximum error | True adjoint equation maximum error | Error in bridge parameter representation | Error in adjoint parameter representation |",
    "| --- | ---: | ---: | ---: | ---: |",
]
for row in equations:
    lines.append(f"| {row['time']} | {number(row, 'true_bridge_equation_max_error'):.2e} | "
                 f"{number(row, 'true_adjoint_equation_max_error'):.2e} | "
                 f"{number(row, 'bridge_representation_error'):.2e} | "
                 f"{number(row, 'adjoint_representation_error'):.2e} |")
lines += [
    "", "A second exact adjoint was constructed with a maximum pointwise "
    f"difference of {number(alternative, 'maximum_pointwise_adjoint_difference'):.1f} "
    "from the first. It still satisfies the conditional equation, changes the "
    f"parameter mean by only {number(alternative, 'parameter_mean_difference'):.2e}, "
    "and changes the population correction using the selected fitted bridge "
    f"by only {number(alternative, 'correction_mean_difference'):.2e}. "
    "This directly demonstrates why pointwise adjoint disagreement is insufficient "
    "to diagnose a bad fit.", "",
    "The selected n = 4,000 functions instead have appreciable equation errors "
    "(reported above). Their errors cannot be explained by the existence of "
    "multiple exact solutions.", "",
    "| Paper condition | DGP audit |",
    "| --- | --- |",
    "| Assumption 1 | Treatment is generated using observed history when Rₜ = 1; the coded Rₜ = 0 branch is deterministic. R₁ = R₂ = 1 under both the natural process and policy. |",
    "| Assumption 2 | Independent structural draws generate treatment and subsequent health and measurement; treatment uses only observed history. This is established by construction. |",
    "| Assumption 3 | Every treatment assigned by the policy has positive natural conditional probability, and every reachable policy history has positive natural probability. |",
    "| Assumption 4 | The measurement equation depends on Hₜ and Cₜ₊₁, with no direct Aₜ term, and uses an independent measurement draw. |",
    "| Assumption 5 | The smallest final measurement probability is 0.299432857526; intermediate measurement is certain. |",
    "| Assumption 6 | Its continuous-treatment sufficient condition is unnecessary here: the discrete policy probability mass and supported ratios exist directly. |",
    "| Assumption 7 | The adjoint equation is solved at every reachable full history; all its four-dimensional complete-health targets are in the operator range. |",
    "| Assumption 16 | Full column rank makes the bridge unique and hence the selected minimum-norm bridge. Positive finite support and full rank persist locally, giving a smooth bridge path within the specified model. |",
    "", "These checks establish existence of solutions for this generating "
    "distribution. They do not establish the product-rate or gradient-convergence "
    "conditions of the fitted learners. The coded no-visit histories have zero "
    "probability here; this audit makes no claim for an alternative mechanism "
    "with missed intermediate visits.", "",
]
lines += [
    "## Do the function classes contain the required functions?", "",
    "The full saturated bridge class assigns a separate coefficient to each "
    "combination of H₂ and C₃. All 256 true bridge values are representable, "
    "and their range, approximately 1.607 to 3.340, lies within the learner's "
    "bounds of 1 to 6. The full saturated adjoint class assigns separate "
    "coefficients to the 512 combinations of A₂ and H₂. It has no sign or "
    "finite magnitude bound; it contains valid solutions to the adjoint equation.", "",
    "The actual million-person bridge fit includes all 256 combinations. "
    "Its adjoint includes 509 of 512, with one common prediction for the "
    "remaining three. An exact linear-system calculation finds an adjoint "
    "solution even within that restricted space: its maximum equation "
    "error is about 1.8 × 10⁻¹⁴. Thus both actual large-fit spaces can "
    "represent the required functions. Points away from y = x in the large "
    "bridge plot are not explained by exclusion of the true bridge from "
    "its function class.", "",
    "At n = 4,000, some combinations are absent from measured training "
    "observations. The code creates independent coefficients only for "
    "combinations observed in training; other combinations receive a common "
    "prediction. This is a finite-data approximation and extrapolation "
    "restriction, also encountered with empty groups in regression. It "
    "should not be called structural misspecification of the population "
    "model. The following table distinguishes the full population class "
    "from the finite coefficient space constructed for each fit.", "",
    "| Space | Function | Separate combination coefficients | Can represent the true bridge or any valid adjoint? |",
    "| --- | --- | ---: | --- |",
]
for row in classes:
    if row["n"] == "0" or row["fold"] == "2":
        scope = "Full population class" if row["n"] == "0" else f"n = {int(row['n']):,}, sample 2"
        name = "Bridge" if row["function_name"] == "beta" else "Adjoint"
        passed = "Yes" if row["exact_member_or_solution"] == "TRUE" else "No, with the common unseen-combination prediction"
        lines.append(f"| {scope} | {name} | {row['represented_targets']} / {row['population_targets']} | {passed} |")
lines += [
    "", "For the selected n = 4,000 bridge fit, 66 combinations are absent, "
    "accounting for 2.61% of the measured generating population. Even the "
    "best representation with a common prediction for these combinations "
    "has RMSE 0.0908; the actual fitted RMSE is 0.6694. The empty-combination "
    "restriction therefore explains only part of the error. Representability "
    "and accuracy of the fitted coefficients are different questions.", "",
    "The sequential regressions were checked for both SDR and TMLE using "
    "the actual pooled predictor construction, including time and the full "
    "observed history. On complete support, their classes contain the true "
    "functions for m₂,₁, m₃,₂, and m₃,₁. The true scaled regression values "
    "range from about 0.0248 to 0.1647, within the TMLE prediction bounds. "
    "Their small-sample combination tables have the same empty-group "
    "limitation; the earlier regression m₃,₁ has all 16 combinations represented.", "",
    "For the treatment-ratio classifier, the categorical predictors contain "
    "the required treatment and history information. The dose policy creates "
    "zero shifted probabilities at dose 0. A saturated logistic model "
    "represents these zero classifier probabilities as a boundary limit, "
    "rather than with finite logit coefficients. This qualification is "
    "recorded explicitly. All positive true classifier probabilities are "
    "below 0.708, so the ratio calculation's upper cap at 0.999 excludes "
    "none of the true ratios. No omitted history term or incorrect pooling "
    "restriction was found in these class checks.", "",
    "[Bridge and adjoint representability](function-class-membership.csv), "
    "[sequential-regression class checks](sequential-regression-class-membership.csv), "
    "[treatment-classifier checks](treatment-classifier-class-membership.csv), "
    "[class audit script](audit_dose_function_classes.R), and "
    "[sequential-class audit script](audit_dose_sequential_classes.R) "
    "record the calculations. They are algebraic evaluations, not known-function "
    "performance simulations.", "",
    "## Why the small-data bridge is nearly constant", "",
    "The saved plot shows most fitted β₂ values near "
    f"{number(constant, 'most_common_fitted_bridge_value'):.6f}, "
    "although the true bridge varies from about 1.607 to 3.340. The selected "
    "L1 penalty suppresses differences across the many complete-history "
    "combinations. There is also an effect from estimating conditional "
    "measurement probabilities in small treatment-and-history groups: squaring "
    "their noisy empirical means changes the fitting criterion.", "",
    "As an algebraic check of the constant term in that criterion, its "
    "preferred constant is "
    f"{number(constant, 'empirical_constant_criterion_minimizer'):.6f} "
    "using the actual training counts, compared with "
    f"{number(constant, 'population_constant_criterion_minimizer'):.6f} "
    "under the exact population criterion. This helps explain why the "
    "regularized fit concentrates near 1.893. This calculation is a criterion "
    "diagnostic; no additional constant learner or zero-penalty simulation was run.", "",
    "The full-history population operators are well conditioned here: their "
    "largest condition number at time 2 is about 2.54. The checks therefore "
    "do not support a failure of population rank or nonexistence of the "
    "solutions as the explanation. The evidence points to finite-sample "
    "estimation with small history groups and regularization. The full "
    "population class contains the required functions; the fitted functions "
    "in this sample "
    "do not satisfy their equations accurately.", "",
    "[Full-history operator audit](full-history-operator-audit.csv), "
    "[true equation and parameter checks](true-equation-and-parameter-audit.csv), "
    "[alternative valid adjoint](alternative-adjoint-check.csv), "
    "[assumption checks](assumption-audit.csv), "
    "[audit log](assumption-audit.log), and "
    "[constant-term criterion calculation](constant-term-criterion-diagnostic.csv). "
    "Paper references: Sections 3 and 5.3.2, Theorem 2, Proposition 1, "
    "and Appendix A.6 of the supplied long_missing-10.pdf.", "",
]
lines += [
    "## Shrinkage in the million-person bridge fit", "",
    "The fitted bridge has one unpenalized common constant plus a separate "
    "penalized deviation for each combination. A deviation shrunk to zero "
    "therefore gives the common bridge value, not a bridge value of zero.", "",
    f"Of its 256 deviation coefficients, {int(number(point_summary, 'zero_deviation_coefficients'))} "
    f"are zero. The fitted common constant is {number(point_summary, 'fitted_common_constant'):.6f}. "
    f"There are {int(number(point_summary, 'points_with_absolute_error_over_half'))} "
    "combinations with absolute bridge error greater than 0.5; "
    f"{int(number(point_summary, 'their_zero_deviation_coefficients'))} of these "
    "have zero deviations. These combinations have between "
    f"{int(number(point_summary, 'their_minimum_measured_training_n'))} and "
    f"{int(number(point_summary, 'their_maximum_measured_training_n'))} measured "
    "training observations and together account for "
    f"{100 * number(point_summary, 'their_measured_population_probability'):.3f}% "
    "of the measured generating population.", "",
    "The two largest errors have true bridge 3.339647 but estimate "
    "2.205937. Their deviations are zero, and they have 24 and 3 measured "
    "training observations. Reconstructing the fitting objective shows "
    "that their coefficient gradients lie below the selected L1 threshold. "
    "The penalty suppresses their deviations. All 64 empirical history "
    "matrices have full bridge rank 4; no empirical rank loss explains these points.", "",
    "To check whether less shrinkage helps, refit only this bridge on the "
    "same 666,667 training people with the three original positive penalty "
    "candidates. The CV-selected fit reproduces the saved predictions to "
    "numerical precision. No generating probabilities or true bridge "
    "values enter any fit; the true values are used afterward for evaluation.", "",
    "| Positive penalty scale | Original CV loss | Bridge prediction RMSE | Maximum absolute error | Zero deviation coefficients |",
    "| ---: | ---: | ---: | ---: | ---: |",
]
for row in shrinkage:
    if number(row, "tolerance") == 1e-8:
        lines.append(f"| {row['scale']} | {number(row, 'original_cv_loss'):.6f} | "
                     f"{number(row, 'beta_rmse'):.6f} | {number(row, 'maximum_absolute_error'):.6f} | "
                     f"{int(number(row, 'zero_deviation_coefficients'))} |")
lines += [
    "", "Reducing the scale from 0.01 to 0.001 improves 12 of the 16 "
    "originally distant points, confirming that shrinkage contributes to "
    "their errors. It also increases overall weighted RMSE from 0.0833 "
    "to 0.0875 and the maximum error from 1.13 to 2.42. The original "
    "CV selection has the lowest generating-distribution weighted RMSE "
    "among these three candidates in this training sample. Thus there is "
    "strong shrinkage in some rare combinations, but the comparison does "
    "not support simply reducing the penalty as an overall improvement.", "",
    "A further numerical check keeps the selected positive penalty and "
    "tightens solver tolerance from 10⁻⁸ to 10⁻¹². Its RMSE is 0.0835, "
    "and its largest prediction change is 0.00573. Numerical stopping "
    "tolerance is therefore too small an effect to explain the large "
    "distances from y = x.", "",
    "![Positive penalty comparison](figures/positive-penalty-comparison.png)", "",
    "These diagnostics identify a tradeoff between shrinkage and noisy "
    "estimation within rare histories. They do not prove that the present "
    "learner or penalty grid is best possible, or establish coverage at "
    "the planned smaller sample sizes. The full repeated-sample study "
    "remains unstarted.", "",
    "[Point-by-point diagnosis](large-bridge-point-diagnosis.csv), "
    "[empirical operator ranks](large-empirical-bridge-operators.csv), "
    "[positive penalty comparison](positive-penalty-shrinkage-comparison.csv), "
    "[distant points under the smaller penalty](shrinkage-far-point-comparison.csv), "
    "[saved candidate coefficients](positive-penalty-bridge-fits.rds), "
    "[comparison execution log](shrinkage-comparison.log), and "
    "[comparison script](check_large_bridge_shrinkage.R) are included.", "",
    "## Integer-overflow warning and package correction", "",
    "The large fit reported an integer-overflow warning. The conditional-moment "
    "scoring code multiplied two integer counts before taking their square "
    "root. In the final bridge training score, 77 count products exceed R's "
    "32-bit integer maximum, giving a missing diagnostic score.", "",
    "The package now converts the sample count to double precision before "
    "multiplication. A regression test uses 50,000 observations in one group, "
    "for which the correct Gram matrix is known exactly. All cmbridge tests "
    "pass with the correction.", "",
    "The original split counts also show that every bridge validation "
    f"count product is at most {int(number(overflow, 'maximum_bridge_validation_count_product_bound')):,}, "
    "and every adjoint count product is at most "
    f"{int(number(overflow, 'maximum_adjoint_count_product_bound')):,}. "
    "Both are below the integer maximum of 2,147,483,647. Hence the selected "
    "penalties and ensemble weights are unaffected in this check. The final "
    "bridge training score is computed after coefficients, predictions, "
    "penalty selection, and weights have been determined.", "",
    "Recomputing that score from the saved fitted coefficients gives "
    f"{number(overflow, 'corrected_final_bridge_training_score'):.9f}, which "
    "matches the direct weighted mean of squared empirical conditional "
    "residuals. No learner refit is needed for this arithmetic correction. "
    "At n = 4,000 the integer products are too small to overflow, so this "
    "bug does not explain the earlier estimation error.", "",
    "[Overflow audit](integer-overflow-audit.csv), "
    "[audit log](integer-overflow-audit.log), "
    "[package test log](package-test.log), "
    "[audit script](audit_dose_scoring.R), and "
    "[corrected package source](cmbridge-ensemble-corrected.R) document the fix. "
    "[Remaining bridge errors](bridge-large-errors.csv) and "
    "[remaining adjoint equation errors](adjoint-equation-large-errors.csv) "
    "include each combination's measured training count and generating probability.", "",
]
(out / "REPORT.md").write_text("\n".join(lines))
files = sorted(path for path in out.rglob("*") if path.is_file() and path.name != "SHA256SUMS")
(out / "SHA256SUMS").write_text("".join(
    hashlib.sha256(path.read_bytes()).hexdigest() + "  " + str(path.relative_to(out)) + "\n"
    for path in files))
print(out / "REPORT.md")
