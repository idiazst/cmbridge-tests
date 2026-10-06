# Observed-history simulation study

## 1. Revised study with missingness at both follow-ups

**The revised full study is authorized after the ensemble and complete estimator checks pass.** The current checks use the saved binary and numerical-dose datasets at n = 4,000 with cmbridge 0.3.0.9020 and lmtp 1.6.0.9021. The new full run has its own sources, library, and result directory, `full-ensemble-study-v2`. The historical study remains stopped with 118 completed datasets in `missing-both-study-v2`; an earlier stopped design saved 176 datasets. Earlier n = 20,000 checks and the n = 10,000 checks remain preserved separately. This study concerns the estimator under observed-history adjustment in Section 3 of the paper. The corrected cross-fitting and cross-validation construction is retained. The revised distributions have passed exact probability checks.

The main study will compare SDR and TMLE, both with the one-step bridge correction, when all four nuisance-function classes are correctly specified. It will assess confidence-interval coverage and estimation error as sample size increases. The two-time binary treatment and two-time numerical dose mechanisms will each use n = 500, 1,000, and 4,000, with 200 replications at every sample size: **1,200 generated datasets in total**. The numerical treatment takes values 0, 1, 2, and 3.

Misspecified models and the empirical assessment of double robustness are deferred. The sample sizes, replication counts, and limits of this study are stated in Section 1.7. The full run will use ten parallel R workers. The earlier 8–9-hour estimate applied to the stopped design with a completely observed intermediate outcome. The revised design requires bridge and adjoint estimation at both times; the initial full-study batches will provide an updated runtime estimate.

Baseline is completely observed. Both follow-up health vectors can be missing, and measurement depends on their outcomes and covariates. Treatment at a missed intermediate visit follows the deterministic rule specified below. Earlier results, including the 176 datasets from the most recent stopped run, remain separate from this revision.

### 1.1 Observed data, treatment policies, and parameters

We use C₁ for baseline information, Cₜ for the complete health vector, Aₜ for the treatment vector, and Rₜ for the measurement indicator. Baseline is always observed, so R₁ is 1 and L₁ is C₁. At subsequent times,

$$
L_t=
\begin{cases}
C_t,&R_t=1,\\
\zeta,&R_t=0,
\end{cases}
\qquad
H_t=(L_1,A_1,\ldots,L_{t-1},A_{t-1},L_t).
$$

The marker ζ has the meaning given in the paper and is distinct from every complete health value. For a vector Xₜ, use the first subscript for its component and the second for time:

$$
X_{j,t}=\text{component }j\text{ of }X_t.
$$

The components used in the study are

$$
\begin{aligned}
A_{1,t}&:\ \text{treatment of interest, binary or numerical dose},\\
A_{2,t}&:\ \text{additional binary treatment used as a shadow variable},\\
C_{1,t}&:\ \text{time-varying covariate},\\
C_{2,t}&:\ \text{outcome component},\\
C_{3,1}&:\ \text{binary baseline characteristic},\\
L_{1,t}&:\ \text{recorded time-varying covariate},\\
L_{2,t}&:\ \text{recorded outcome}.
\end{aligned}
$$

Thus the vectors and outcome projection are

$$
\begin{gathered}
A_t=(A_{1,t},A_{2,t}),\qquad
C_1=(C_{1,1},C_{2,1},C_{3,1}),\\
C_t=(C_{1,t},C_{2,t}),\quad t\geq2,\\
Y_{t+1}=h(C_{t+1})=C_{2,t+1}.
\end{gathered}
$$

At a measured visit, L₁,ₜ and L₂,ₜ equal C₁,ₜ and C₂,ₜ, respectively. Baseline also includes L₃,₁ = C₃,₁. At an unmeasured time, the entire Lₜ is ζ. The policy modifies A₁,ₜ and retains A₂,ₜ. The additional treatment A₂,ₜ predicts subsequent health and provides variation needed for the dual range condition in Assumption 7; the generating measurement equation also respects the shadow variable exclusion in Assumption 4. These conditions will be checked as specified in Section 1.3.

The baseline covariate C₁,₁ is zero. The baseline outcome C₂,₁ equals the binary baseline characteristic C₃,₁. All generated variables have discrete, finite support, including numerical dose.

| Data-generating mechanism | Treatment times | A₁,ₜ | Treatment policy at a visit | Main parameter |
|:--|:--|:--|:--|:--|
| Two-time binary treatment | τ = 2 | 0 or 1 | At each time, set A₁,ₜ to 1; retain A₂,ₜ | θ₃(d̄₂) |
| Two-time numerical dose | τ = 2 | 0, 1, 2, or 3 | At each time, increase A₁,ₜ by one, capped at 3; retain A₂,ₜ | θ₃(d̄₂) |

The parameters are complete population means under the modified treatment policies:

$$
\theta_{t+1}(\bar d_t)=E\{Y_{t+1}(\bar d_t)\}.
$$

Both mechanisms will also estimate θ₂(d̄₁), so that both methods return the entire counterfactual mean curve. R₂ and R₃ are both informative, so both outcomes require estimated bridge and adjoint functions. Neither β₁ nor β₂ is supplied as a known function to an estimator.

Treatment is generated randomly only when Rₜ = 1. When R₂ = 0, the code sets the entire second treatment vector to the deterministic rule α₂(H₂) = (0, 0), and the policy leaves it unchanged. The final outcome then equals the observed baseline outcome. The final covariate continues to be drawn randomly. Section 1.3 explains why this outcome distribution is needed for the adjoint condition in this design.

### 1.2 Data-generating mechanisms

In the equations below, logit⁻¹ converts log-odds to a probability:

$$
\operatorname{logit}^{-1}(x)=\frac{\exp(x)}{1+\exp(x)}.
$$

**Both mechanisms.** The random treatment distributions below apply conditional on Rₜ = 1. At such times, A₁,ₜ and A₂,ₜ are conditionally independent given Hₜ. At times without a visit, both treatment components are assigned deterministically:

$$
P\{A_t=\alpha_t(H_t)\mid H_t,R_t=0\}=1,
\qquad
\alpha_t(H_t)=(0,0).
$$

Independent draws are used for each random conditional distribution. R₁ is always one; both values of R₂ occur.

$$
C_{3,1}\sim\operatorname{Bernoulli}(0.5),\qquad
C_{1,1}=0,\qquad C_{2,1}=C_{3,1},\qquad
R_1\sim\operatorname{Bernoulli}(1).
$$

$$
A_{2,t}\mid H_t,R_t=1
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(-0.15+0.30L_{2,t}-0.20C_{3,1}\right)
\right].
$$

$$
\begin{gathered}
C_{1,t+1}\mid A_t,H_t,R_t=1\\
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(
-1.8+3.6A_{2,t}+0.35A_{1,t}+0.30L_{1,t}+0.25C_{3,1}
\right)\right].
\end{gathered}
$$

At a missed intermediate visit:

$$
\begin{gathered}
C_{1,3}\mid H_2,R_2=0
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(-1.8+0.25C_{3,1}\right)\right],\\[6pt]
P\{Y_3=L_{2,1}\mid H_2,R_2=0\}=1.
\end{gathered}
$$

At each measured treatment visit, t = 1, 2:

$$
\begin{gathered}
R_{t+1}\mid H_t,C_{t+1},R_t=1\\
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(
0.25-0.65Y_{t+1}-0.45C_{1,t+1}+0.15L_{2,t}+0.10C_{3,1}
\right)\right].
\end{gathered}
$$

At a missed intermediate visit:

$$
\begin{gathered}
R_3\mid H_2,C_3,R_2=0\\
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(
0.25-0.65Y_3-0.45C_{1,3}+0.10C_{3,1}
\right)\right].
\end{gathered}
$$

**Two-time binary treatment:** τ = 2, t = 1, 2. It uses

$$
\begin{gathered}
A_{1,t}\mid H_t,R_t=1
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(-0.20+0.45L_{1,t}-0.40L_{2,t}+0.15C_{3,1}\right)
\right],\\[6pt]
Y_{t+1}\mid A_t,H_t,C_{1,t+1},R_t=1\\
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(
-1.8+3.6A_{1,t}+0.35C_{1,t+1}+0.30L_{2,t}+0.20C_{3,1}
\right)\right].
\end{gathered}
$$

**Two-time numerical dose:** τ = 2, t = 1, 2. The conditional probability mass function for A₁,ₜ is normalized to sum to 1:

$$
p(A_{1,t}=a\mid H_t,R_t=1)\propto
\begin{cases}
1,&a=0,\\[3pt]
\exp\left(0.15+0.30L_{1,t}-0.20L_{2,t}\right),&a=1,\\[3pt]
\exp\left(-0.10+0.40L_{1,t}-0.30L_{2,t}+0.10C_{3,1}\right),&a=2,\\[3pt]
\exp\left(-0.35+0.50L_{1,t}-0.40L_{2,t}+0.20C_{3,1}\right),&a=3.
\end{cases}
$$

$$
\begin{gathered}
Y_{t+1}\mid A_t,H_t,C_{1,t+1},R_t=1\\
\sim\operatorname{Bernoulli}\!\left[
\operatorname{logit}^{-1}\!\left(
-1.8+1.8A_{1,t}+0.35C_{1,t+1}+0.30L_{2,t}+0.20C_{3,1}
\right)\right].
\end{gathered}
$$

### 1.3 Assumptions and calculation of the true parameters

The probability distributions will be checked before any estimator simulations. These checks apply to every history with positive probability under the natural treatment mechanism or the policy. A failed condition will prevent simulation under that mechanism.

| Condition in the paper | Verification in the proposed study |
|:--|:--|
| Assumption 1: observed-history treatment mechanism | For Rₜ = 1, treatment is drawn from the specified conditional distributions using Hₜ. For R₂ = 0, the code sets A₂ = α₂(H₂) = (0, 0), and the policy retains that value. Both branches occur with positive probability. |
| Assumption 2: sequential randomization conditional on observed history | The independent treatment, health, and measurement draws establish the required independence by construction. |
| Assumption 3: supported intervention | Compute natural and policy probabilities exactly and check that every treatment assigned by the policy has positive conditional natural probability on the relevant histories. |
| Assumption 4: shadow variable exclusion | Measurement depends on Hₜ and Cₜ₊₁, with no direct dependence on Aₜ once these are fixed. |
| Assumption 5: positivity of missingness mechanism | Every relevant measurement probability is strictly positive; the smallest is approximately 0.2994. |
| Assumption 7: dual range condition for the representer | Calculate Tₜ⋆ and φₜ exactly and verify the equation Tₜ⋆λₜ = φₜ at every relevant history. Treatment–health association alone will not be used as a substitute for this check. |
| Assumption 16: regularity of the selected visit-dependent bridge path | On the finite support, the conditional bridge matrix has full column rank at measured treatment visits and rank one at missed visits. Those ranks persist locally, and probabilities are positive on the reachable support. The weighted minimum-norm solution is therefore differentiable along regular submodels on that support. |

Because treatment is discrete, the existence and support of gₜᵈ are checked directly. We do not require the continuous-treatment sufficient condition in Assumption 6. The exact bridge and adjoint equations have errors below 4 × 10⁻¹⁵ on every reachable full observed history. These checks establish properties of the specified distributions; they do not establish convergence rates of fitted nuisance functions.

At R₂ = 0, A₂ is fixed by H₂, so the adjoint can only return a function of H₂. Keeping a random final outcome at that fixed history would generally violate Assumption 7 whenever Ω₁ is positive. The revised distribution sets the final outcome equal to the observed baseline outcome at these histories. Hence

$$
\omega_2(A_2,H_2)=1,\qquad
\phi_2(H_2,C_3)=\Omega_1 L_{2,1},\qquad
\lambda_2(A_2,H_2)=\Omega_1 L_{2,1},
\qquad R_2=0.
$$

This gives a solution to the adjoint equation. A valid bridge at the same histories is

$$
\beta_2(H_2,C_3)=\frac{1}{P(R_3=1\mid H_2)},\qquad R_2=0.
$$

It satisfies the bridge equation but need not be the only solution. The audit also evaluates a different solution and verifies that both produce the same parameter. At R₂ = 1, the outcome remains random according to the stated conditional distribution.

The expected proportions missing in naturally generated data, calculated by exact summation, are:

| Time | Two-time binary treatment | Two-time numerical dose |
|:--|--:|--:|
| t = 1 | 0% | 0% |
| t = 2 | 54.32% | 56.82% |
| t = 3 | 53.54% | 54.33% |

For the true parameter, we will use the known generating probabilities, without fitting models and without estimating the bridge:

1. Begin with the exact distribution of C₁.
2. At each treatment time, calculate gₛᵈ by applying dₛ to every possible natural treatment value.
3. Sum over each subsequent complete health vector and measurement indicator, updating the distribution of the observed history under the policy.
4. At the outcome time, average h(Cₜ₊₁) over this complete distribution, including unmeasured outcomes.

A separate calculation will work backward through the known conditional distributions, beginning with h(Cₜ₊₁). The two calculations must agree within 10⁻¹². As further checks, the bridge equation, adjoint equation, and both representations in Theorem 2 will be evaluated by exact summation:

$$
\theta_{t+1}(\bar d_t)
=E\{R_{t+1}\Omega_tY_{t+1}\beta_t(H_t,C_{t+1})\}
=E\{\lambda_t(A_t,H_t)\}.
$$

The forward and backward calculations already agree for these proposed distributions. The values below are known parameter values, not estimates from replacement simulations.

| Data-generating mechanism | θ₂(d̄₁) | θ₃(d̄₂) |
|:--|--:|--:|
| Two-time binary treatment | 0.900308733205 | 0.649747477399 |
| Two-time numerical dose | 0.851702239961 | 0.633624046382 |

### 1.4 Estimators and nuisance-function estimation

The comparison will use the cross-fitted sequentially doubly robust estimator and longitudinal TMLE supplied by the modified lmtp package. Both will receive the same one-step bridge correction in equation (7):

$$
\widehat\theta_{t+1}
=\widehat\theta^{\mathrm{seq}}_{t+1}
+\frac{1}{n}\sum_{i=1}^{n}
\widehat D^{\mathrm{br}}_{t,j(i)}(O_i).
$$

The terminal pseudo-outcome and bridge contribution are

$$
\eta_{t+1,t+1}(O)
=R_{t+1}Y_{t+1}\beta_t(H_t,C_{t+1}),
$$

$$
D_t^{\mathrm{br}}(O)
=-\lambda_t(A_t,H_t)
\{R_{t+1}\beta_t(H_t,C_{t+1})-1\}.
$$

When Rₜ₊₁ is 0, the terminal pseudo-outcome is zero and the expression in braces is −1. Both complete health components will be withheld from the fitted estimators at that time. Computing either expression does not require evaluating unobserved Cₜ₊₁ or Yₜ₊₁.

The bridge and adjoint will be estimated using the conditional moment restrictions in the paper:

$$
E\{R_{t+1}\beta_t(H_t,C_{t+1})-1\mid A_t,H_t\}=0,
$$

$$
E\{\lambda_t(A_t,H_t)\mid R_{t+1}=1,H_t,C_{t+1}\}
=\phi_t(H_t,C_{t+1}).
$$

The policy-induced probability mass function gₛᵈ is obtained by applying dₛ to the natural treatment distribution gₛ and summing the probabilities of all treatment values mapped to the same destination. The policy density ratios and cumulative products are

$$
\omega_s(a_s,h_s)=\frac{g_s^d(a_s\mid h_s)}{g_s(a_s\mid h_s)},
\qquad
\Omega_s=\prod_{j=1}^{s}\omega_j(A_j,H_j).
$$

Policy density ratios will be estimated by the probabilistic-classification construction in Section 5.2.1. The sequential regressions mₜ₊₁,ₛ will be estimated backward from the terminal pseudo-outcome. Following Section 3.2, regressions on the same diagonal will be pooled and time included as a predictor. Saturation in time permits separate regression functions at each time and is a special case of the paper's pool-and-smooth construction. This study will assess the estimator with these saturated learners; it will not assume that pooling improves performance.

SuperLearner will combine treatment and regression learners, and cmbridge will combine bridge and adjoint learners. The candidate libraries are:

| Function estimated | Candidates |
|:--|:--|
| βₜ | `sieve_md` with fixed main effects and saturated joint-category deviations, with CV-selected ridge penalty; `landweber` with fixed main effects in every predictor; `pmmr`. All use inverse-expit. |
| λₜ | Saturated joint-category `sieve_md`, `landweber`, and `pmmr` |
| ωₛ, through probabilistic classification | `SL.saturated_l1_cv` and `SL.mean` |
| mₜ₊₁,ₛ | `SL.saturated_l1_cv`, MARS (`SL.earth`), and `SL.mean` |

The joint-category basis gives every distinct combination of the supplied discrete inputs a separate coefficient. It includes interactions of every order. The revised bridge sieve also includes an intercept and a linear term for every supplied predictor. These terms are defined before observing the training combinations, and an unseen combination therefore retains its main-effect prediction. Joint-category deviations allow additional interactions. The fixed main-effect terms in both the sieve and Landweber contain a valid true bridge solution on the entire reachable support of both mechanisms. An algebraic audit verifies this using a basis constructed from a single training combination; the coefficients from that audit are never supplied to fitted learners. At missed visits, the represented solution may differ from the particular solution used in earlier tables, because the bridge need not be unique there. Comparisons against a single true function therefore focus on the histories where it is unique, with conditional-equation checks for every history.

The preceding complete quadratic Landweber basis remains available in cmbridge and in the preserved diagnostic runs. Its class also contained a true solution, but its fitted values were variable in the n = 4,000 checks. The current Landweber uses a fixed degree-one polynomial basis in every input. No predictor is omitted using the generating equations. Saturated L1 is absent from both conditional-moment libraries and remains in regression and classification.

| cmbridge candidate | Basis and regularization |
|:--|:--|
| `sieve_md` | Bridge: an intercept, every fixed main-effect term, and joint-category indicators (`cell_linear`). Ridge penalizes the scaled link coefficients, with a factor of 100 on joint-category deviation coefficients; the intercept is unpenalized. The conditioning basis has an intercept and every main-effect term. Adjoint: unrestricted joint-category target and conditioning bases with the preceding function-value ridge. In both fits, training-only CV chooses lambda from the initial 41 positive values from 1e-10 to 1, extending boundary minima until interior or failure. |
| `landweber` | Bridge: fixed degree-one polynomial target and conditioning bases in every predictor. An invertible change in coefficient coordinates improves numerical conditioning and preserves the entire function class. Adjoint: the preceding additive cubic B-spline bases, with 6 target and 9 conditioning degrees of freedom per input. Fitting uses at most 2,000 iterations and tolerance 1e-10. |
| `pmmr` | Gaussian-kernel approximation with 40 target and 80 conditioning centers chosen within training data, with lambda = 1e-4. This finite approximation need not individually contain the true functions. |

Joint-category conditioning indicators include a redundant intercept. The revised implementation removes that redundancy algebraically while preserving the ridge-weighted fitting loss. Large adjoint ridge values are solved with the unpenalized common direction separate from the penalized directions, so the penalty cannot erase its curvature numerically. The original function-penalized bridge sieve is likewise evaluated with separate moment and penalty terms at large ridge values.

The two SuperLearner libraries also retain saturated L1 models. `SL.mean` estimates a common mean. MARS fits piecewise-linear spline terms and permits pairwise interactions with the standard `SL.earth` setting `degree = 2`; it is an additional regression candidate, not a saturated model. It uses its standard training-data pruning rule and no separate internal cross-validation (`nfold = 0`). Its ensemble validation uses the supplied person-based learner splits. MARS is used only for the sequential outcome regressions, not in cmbridge or treatment-ratio classification. The [SuperLearner implementation of SL.earth](https://github.com/ecpolley/SuperLearner/blob/master/R/SL.earth.R) records these settings.

The L1 penalty is inside `SL.saturated_l1_cv` for ωₛ and mₜ₊₁,ₛ. Neither conditional-moment library contains this learner. The penalty adds a cost for the sum of absolute deviation coefficients to the regression or classification fitting loss. The constant is unpenalized.

Each saturated L1 regression or classification learner's penalty will be chosen by cross-validation within its training sample using validation prediction loss. The initial grid contains 21 positive penalty scales equally spaced on the logarithmic scale from 10⁻⁵ to 1, divided by the square root of the number of distinct people in the corresponding training fit. A minimum at either grid boundary is recorded as a failed selection attempt. Extend that side by one decade with four logarithmically spaced values and repeat on the same training-only splits. Accept only an interior minimum. After 16 unsuccessful extensions, record an unresolved boundary error rather than fitting a selected boundary penalty. Exact ties may select an interior point only when it also attains the minimum. Refit the selected learner on its training sample. SuperLearner and cmbridge will choose convex-combination weights across their respective candidate libraries using shared person-based validation assignments. Every candidate's weight and validation loss will be recorded. The sieve selects its ridge penalty using raw Gaussian-kernel U-statistics on shared person-based training-only splits. It tunes separately within every ensemble training fold and again for its final refit; nested bridges in the sequential-response construction use the same selection rule. Landweber retains iteration stopping and PMMR retains the fixed penalty listed above. When a candidate learner is evaluated on a validation sample, its own penalty selection will use only the corresponding training observations. These choices will use no outer validation observations and no information about the true parameter or final coverage.

All three bridge candidates optimize unrestricted real coefficients through the inverse-expit parameterization. The sieve and PMMR penalize the linear predictor inside that link; Landweber uses iteration stopping. Predictions are mathematically greater than one, with one as a limiting value, and their convex ensemble has the same range. All three adjoint candidates retain unrestricted identity parameterizations. No conditional-moment learner clips predictions after fitting.

For the sequential estimator, the wider interval [−100, 100] may be applied to bridge predictions. Record the number and magnitude of changed values and their moment-equation errors, and retain raw predictions for plots. Outcomes in these mechanisms are binary, so the same interval contains the transformed terminal outcomes. Keep the logistic TMLE update. Both bridges are estimated, and both outcome means receive the one-step bridge correction.

All fitted nuisance functions will use the arguments given in the paper: βₜ uses Hₜ and Cₜ₊₁; λₜ uses Aₜ and Hₜ; ωₛ and mₜ₊₁,ₛ use Aₛ and Hₛ. No correctly specified fit will omit parts of these histories because the generating equations happen to make them unnecessary. Correct specification means that the function class contains the relevant target. For the sequential regressions, that target is defined using the supplied bridge, as in Proposition 1. For diagnostics, the exact sequential-regression target will be recalculated using each estimated bridge. Class containment will not be treated as proof of the finite-sample accuracy or product-rate conditions of a fitted learner.

Ensemble weights and bridge/adjoint ridge choices in this revision use a Gaussian-kernel U-statistic, removing all self-products. Kernel scales and centers are chosen inside the corresponding training sample and retain every supplied predictor, including treatment. The kernel uses the standard training-data Gaussian bandwidth and at most 100 centers. A penalty trial with a fitting error or nonfinite score is retained as failed and cannot be selected. A conditional-moment candidate with no valid tuning/fitting/validation result receives weight zero across every fold and its final refit; record its failure and the remaining library. If every candidate fails, record an estimator failure. This does not accept a boundary minimum. The raw fold matrices are averaged with their observation-count weights, retained, and projected to PSD before choosing ensemble weights. Penalty selection uses the raw scalar U-statistic, which may be negative. Conditional-cell U-statistics remain implemented and verified, including the rule that singleton cells contribute zero. Their cell counts will still be reported as diagnostics. No candidate uses outer validation data for its fitting, tuning, or weight selection.

### 1.5 Correctly specified functions and the second-order remainder

Every fitted replication in the main study will use correctly specified classes for all four nuisance functions. No misspecified-model comparison is included in this reduced design. All nuisance functions will be estimated from the generated data.

| Nuisance function | Model used in every replication |
|:--|:--|
| βₜ at the informatively measured outcome | Three-candidate cmbridge ensemble using Hₜ and Cₜ₊₁, with inverse-expit for all three candidates, including the fixed main-effect and joint-category sieve; no saturated L1 candidate |
| λₜ | Three-candidate cmbridge ensemble using Aₜ and Hₜ, including the saturated sieve class; no saturated L1 candidate |
| ωₛ | SuperLearner classification ensemble using treatment and Hₛ, including a saturated L1 class |
| mₜ₊₁,ₛ | SuperLearner regression ensemble for the corresponding pseudo-outcome, including saturated L1, MARS, and a mean |

“Correctly specified” means that the function class can represent a required function or a valid solution. The fitted functions still have estimation error. Positive L1 penalties will be selected by cross-validation within training samples, as described in Section 1.4; zero-penalty and known-function performance simulations are excluded. All nuisance functions are estimated.

The adjoint target uses the policy density ratios through

$$
\phi_t(H_t,C_{t+1})
=\Omega_{t-1}
E\{Y_{t+1}\omega_t(A_t,H_t)
\mid R_{t+1}=1,H_t,C_{t+1}\}.
$$

Those ratios will also be estimated with correctly specified classes and the shared sample-splitting restrictions. All fits retain the complete observed histories specified in the paper; none omit predictors using the generating restrictions.

Proposition 1 gives the bound

$$
\begin{aligned}
\left|R_{t,2}(P,\widehat P)\right|
\lesssim{}&
\sum_{s=1}^{t}
\left\|\omega_{s,P}-\omega_{s,\widehat P}\right\|_{2,P}
\left\|m^{\beta_{\widehat P}}_{t+1,s,P}
      -m^{\beta_{\widehat P}}_{t+1,s,\widehat P}\right\|_{2,P}\\
&+
\left\|\lambda_{t,P}-\lambda_{t,\widehat P}\right\|_{2,P}
\left\|\beta_{t,P}-\beta_{t,\widehat P}\right\|_{2,P}.
\end{aligned}
$$

The reduced study will measure the estimation errors and the separate contributions to the remainder when all classes are correct. It will examine whether the tested sample sizes are sufficient for the remainder to be negligible relative to the sampling uncertainty. It will not test double robustness under misspecification, because that question requires comparisons in which some nuisance functions are misspecified.

If the adjoint equation has more than one solution, difference from one chosen exact solution alone will not be called misspecification. Its fitted conditional equation and the remainder will be evaluated directly.

### 1.6 Shared cross-fitting and cross-validation

For each generated dataset, divide people into J = 3 validation samples 𝒱ⱼ and corresponding training samples

$$
\mathcal T_j=\{1,\ldots,n\}\setminus\mathcal V_j.
$$

The outer program will generate these assignments once and pass them to both lmtp and cmbridge. Within each 𝒯ⱼ, three validation samples will be used for cross-validation of candidate learners and convex-combination weights. Penalty selection inside a candidate fit will use further training and validation assignments restricted to that fit's training observations; these assignments will also be shared by the two packages where applicable. Within each outer training sample, learner-validation assignments will balance the measurement patterns using only that training sample. Three learner-validation groups retain more observations in each training fit than two groups. Further groups inside candidate training samples remain shared between packages. Observations from the same person will retain the same assignment across time and across numerator and denominator rows used in classification.

Preliminary nuisance functions evaluated on 𝒱ⱼ will be fitted without 𝒱ⱼ. As in lmtp, the TMLE intercept update uses 𝒱ⱼ after those preliminary fits; that update does not enter the learner training samples. Within 𝒯ⱼ, estimated pseudo-outcomes and adjoint responses used for cross-validation will also use nuisance predictions obtained without their own validation observations. Additional training and validation samples will be formed inside 𝒯ⱼ where this is needed to construct those responses. For each candidate or penalty validation trial, the upstream functions used to construct its training responses will be fitted only within that trial's training people. Their predictions for those training responses will themselves be out of sample. The upstream functions used to construct its validation responses will likewise use only that trial's training people. Candidate fits, tuning choices, and ensemble weights will all exclude the outer validation sample.

The same assignments will be propagated through the entire backward recursion. The sequentially doubly robust estimator and TMLE will use the same datasets, assignments, and bridge and adjoint fits.

Additional validation inside training samples normally uses two groups. If a current training sample has 16 or fewer measured people at any incompletely measured outcome time, it uses the configured three groups to retain more measured observations in subsequent training fits. The rule uses only measurement indicators in that training sample, and both packages receive the same assignments. The primary estimator groups and learner-validation groups remain three. A bridge validation group with no measured outcomes still has a defined bridge residual and can be scored; its training group must have measured outcomes. An adjoint validation group must have measured outcomes. Failed fits remain recorded.

Before the expanded study, cmbridge's spline basis was corrected for training-constant inputs and for the case where all internal knots coincided with a boundary. The latter knots are moved inward by one eighth of the training range, matching the spacing used by the standard spline routine's usual boundary adjustment. This prevents nonfinite extrapolation when a new discrete category is outside the observed training range. Those historical runs applied bridge bounds after fitting; the current revision removes that operation from the learners. Earlier attempts, their errors, and their source snapshots are preserved.

The first launch with missingness at both follow-ups exposed a further spline prediction defect when a validation group had no measured outcomes. cmbridge 0.3.0.9005 returns an empty design with the fitted column count for these zero-row inputs. The learner definitions and predictions on nonempty inputs are unchanged. Its 115 package tests pass. The interrupted attempt is preserved in `missing-both-study`; the corrected full study uses `missing-both-study-v2`.

### 1.7 Sample sizes, replications, and expected runtime

Sample size n is the number of independent people in a generated dataset. A replication is a separately generated dataset with its recorded random seed. Every replication will fit both estimators with all four nuisance-function classes correctly specified.

| Data-generating mechanism | Main parameter | Sample sizes n | Replications at each n | Total datasets |
|:--|:--|:--|--:|--:|
| Two-time binary treatment | θ₃(d̄₂), also reporting θ₂(d̄₁) | 500; 1,000; 4,000 | 200 | 600 |
| Two-time numerical dose | θ₃(d̄₂), also reporting θ₂(d̄₁) | 500; 1,000; 4,000 | 200 | 600 |
| **Total** | | | | **1,200** |

The smallest sample size is 500 for every mechanism. Each of the six mechanism–sample-size combinations has 200 independent datasets. Both estimators use each dataset, giving 2,400 estimator fits from 1,200 independent datasets.

This is the entire proposed grid; there is no separate larger-sample extension. The computer has 10 CPU cores: eight performance cores and two efficiency cores, with 32 GB of RAM. The proposed run will use ten parallel R workers. Each worker will process a separate dataset, and both estimators will continue to share the data and sample assignments within that dataset.

The earlier 8–9-hour estimate used the stopped design with no missingness at t = 2. It does not apply directly to this revision, which estimates bridge and adjoint functions at both times. The initial full-study batches will provide timings for all six mechanism/sample-size combinations, without adding separate pilot datasets. The separate million-person component check is outside the planned grid and is not used to extrapolate this duration.

With 200 replications, an observed coverage fraction near 0.95 has simulation standard error approximately 0.015, or 1.5 percentage points, for every mechanism and sample size. Tables and graphs will show this uncertainty. Consequently, the study can identify substantial undercoverage more reliably than small differences from 0.95.

Each mechanism has three sample sizes. These permit comparison of estimation error and √n times average estimation error as sample size increases. Because the largest n is 4,000, the study will describe behavior over this range; it cannot establish that a limiting convergence regime has been reached.

The revised counts reduce computation; they have not been selected to obtain a desired coverage fraction. Failures and warnings will be retained and reported, and replication counts will not be silently reduced to meet the deadline.

### 1.8 Confidence intervals, tables, graphs, and investigation of unexpected results

Standard errors will use the combined gradient from the paper,

$$
D_t(O)=D_t^{\mathrm{seq}}(O)+D_t^{\mathrm{br}}(O).
$$

The empirical variance of the cross-fitted estimated gradient, divided by n, will estimate the variance of θ̂ₜ₊₁. This includes covariance between its sequential-regression and bridge contributions. Pointwise confidence intervals will use the estimate plus or minus 1.96 times its estimated standard error. They will not be clipped to the outcome range. The estimates entering the comparison will likewise retain the value produced by equation (7).

For the two-time counterfactual mean curves, the report will also evaluate simultaneous 95% confidence bands using the joint multiplier-bootstrap distribution described in Section 3.2, with 2,000 multiplier draws per fitted replication. It will distinguish the fraction of intervals containing an individual true parameter from the fraction of bands containing both true curve values.

Tables will give the true value, average estimate, average estimation error, standard deviation across replications, square root of the average squared estimation error, average estimated standard error, average interval length, and fraction of intervals containing the truth. They will state the number of planned, successful, and failed replications. The number of warnings will also be reported. Separate tables will report the mean of each EIF component over outer validation observations, evaluated at the known true parameter, both for each replication and averaged across successful replications. Bias and coverage conclusions will be shown separately for every mechanism, outcome, estimator, and sample size. All rows will use the same correctly specified function classes.

Graphs will show average estimation error against n, √n times the average estimation error against n, coverage against n with a reference at 0.95, and the standard deviation of estimates alongside the average estimated standard error. For the two-time mechanisms, graphs will also show the estimated counterfactual mean curves against their true values.

For correctly specified models, the paper's confidence-interval expectation will be assessed together with its requirements: convergence of the estimated gradient and a second-order remainder that is negligible relative to n⁻¹ᐟ².

For each fitted replication, store the nuisance fits and calculate conditional-moment errors, policy ratios, and the separate sequential-regression and bridge contributions to the population estimation error by summing over the known generating distribution. Summaries of these quantities will use all replications, rather than one selected dataset. Also record the numbers of training observations and measured observations at each predictor combination, the selected penalties and ensemble weights, and any numerical failure.

If a result is unexpected, first recheck the true parameter and assumption calculations, then compare the public package output with a direct calculation of equation (7) using the same predictions. Check sample assignments and pseudo-outcome construction, conditional equations, penalization, observations available at each predictor combination, and gradient-based standard errors. A demonstrated implementation error will be fixed and the affected comparisons repeated using the original datasets and assignments. Results before and after the fix will be retained. A result caused by fitted-function estimation error will be reported as such; increasing n or changing a learner will be a separately documented comparison. Simulation results will not be discarded or modified to make them agree with an expected conclusion.

The distributions with missingness at both follow-ups have been checked by exact calculation. For both mechanisms, independent forward and backward calculations and both representations in Theorem 2 agree within 10⁻¹². Bridge and adjoint equations were checked on every reachable full observed history. Full joint-category classes contain bridge and adjoint solutions at both times. The audit itself generates no sample and fits no estimator. The stopped run, its 176 checkpoints, and its recorded fit failure are preserved separately; those results do not validate the revised distributions.

The separate [n = 10,000 check](../results/observed-history-corrected/missing-both-n10000/REPORT.md) completed successfully for both mechanisms and both estimators in approximately seven minutes with two workers. It saves each person's sequential-regression and bridge/adjoint contributions, their means, standard deviations and covariance, and the separate second-order population remainders. Direct evaluation reconstructs the package estimates and standard errors. All eight pointwise intervals contain their true values in this check; one dataset per example cannot assess coverage or convergence. The dose final-outcome bridge/adjoint population remainder is approximately −0.0202 and remains visible in the report despite the final estimates being close to the true mean.

The earlier dose timing checks, under the scheduled intermediate-visit design, had substantial downward estimation error for θ₃(d̄₂), including at n = 4,000. The saved population calculations show nonnegligible bridge and sequential-regression remainder contributions. At n = 4,000, many full-history combinations still have few training observations, and all three outer training samples selected the largest candidate bridge and adjoint penalties. These historical checks do not estimate repeated-sample coverage or validate the current revision. They show why correctly specified classes alone are insufficient to predict good finite-sample coverage. The revised distributions were changed to add the requested missingness while satisfying the paper's assumptions, rather than in response to coverage results.

The preceding larger study was stopped after 19 completed datasets; those preliminary results remain separately identified. The generating probabilities will be used for assumption checks, true-value calculations, and evaluation of fitted-function errors; the simulations will not supply those probabilities or exact nuisance functions to the fitted estimators.

Source: Sections 2, 3.1–3.2, and 5, and Appendix A.2 and A.6–A.7 of the supplied paper. The [assumption checks](../results/observed-history-corrected/missing-both-design/assumption_audit.csv), [full-history conditional equations](../results/observed-history-corrected/missing-both-design/full-history-equations.csv), [function classes](../results/observed-history-corrected/missing-both-design/function-classes.csv), and [independent true-value calculations](../results/observed-history-corrected/missing-both-design/truth.csv) are available separately.
