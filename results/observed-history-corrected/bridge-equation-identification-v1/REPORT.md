# Which bridge equations are enforced by the current basis?

The saturated sieve **function class contains a valid bridge**, but its current additive conditioning basis does not enforce all the paper's conditional equations. These are separate requirements. This audit demonstrates the distinction directly on the full discrete population in both approved mechanisms.

The paper requires the conditional mean of the measured bridge, given the observed history and treatment, to equal one. The present degree-one conditioning basis checks averages against the intercept and individual predictors. With a flexible joint-category bridge, satisfying those averages can leave errors within particular combinations of predictors.

For each outcome time, the audit constructed a different bridge in the complete joint-category inverse-expit class. Every value remains strictly above one. It satisfies every fitted additive moment to numerical precision, yet violates the full conditional equation. No datasets or nuisance estimators were fitted; this is an algebraic calculation, not a known-function or zero-penalty simulation.

| Treatment | Outcome time | Independent fitted equations | Independent full conditional equations | Full equation RMSE of the alternative | Largest fitted-moment error |
| --- | --- | --- | --- | --- | --- |
| Binary | 2 | 4 | 8 | 0.062738 | 5.55e-16 |
| Binary | 3 | 9 | 136 | 0.056425 | 4.44e-16 |
| Numerical dose | 2 | 4 | 8 | 0.044378 | 2.22e-16 |
| Numerical dose | 3 | 9 | 272 | 0.036217 | 8.88e-16 |

The reference bridge satisfies the full equations with maximum error below 1.33e-15. The alternative functions have values between 1.383 and 3.528, and their representation errors in the sieve class are zero to recorded precision. Thus the discrepancy is not caused by clipping or lack of a valid bridge in the DGP.

This does **not** show that an actual fit selects the constructed alternative, or that the whole ensemble is inconsistent. The fixed main-effect Landweber class is smaller: the number of independent linearized fitted equations equals the number of independent target-basis coefficients at the valid bridge in all four cases. That is a local check around the valid solution, not a proof of global uniqueness. PMMR was not given this identification audit.

Adding a full joint-category conditioning basis removes this particular population limitation. However, the [separate estimated-model comparison](../bridge-conditioning-check-v1/REPORT.md) made finite-sample predictions worse in an existing n=4000 dataset. Therefore that change has not been substituted into the running study. More conditional equations also require estimating more conditional averages from the same observations; population correctness alone does not establish better finite-sample performance.

The frozen study and its statistical sources are preserved. Any future fitting change requires a separate version and new checks. [Numerical audit](bridge-equation-identification.csv), [all conditional errors of the constructed alternatives](alternative-conditional-equations.csv), and [reproducible algebraic script](../../../simulations/observed_history_corrected/audit_bridge_defining_equations.R) retain the evidence.
