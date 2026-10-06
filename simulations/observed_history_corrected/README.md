# New diagnostic configurations (October6)

cmbridge0.3.0.9021 now supports optional positive training-only CV for Landweber's conditioning-weight ridge. Modified lmtp remains1.6.0.9021. All GitHub validation checks passed for package847e187026ab5a68c456db0f9d49496ff7a869fc. The already running full cloud study retains its immutable0.3.0.9020 archive and original fingerprint.

The [v7 saved n4000 checks](../../results/observed-history-corrected/ensemble-stability-v7-full-equations/REPORT.md) use full joint-category bridge conditioning and cell U-statistic selection. They completed, but worsened dose bridge RMSE; this is not an adopted repair. [Conditioning-weight diagnostics](../../results/observed-history-corrected/bridge-conditioning-weight-check-v1/REPORT.md) include all comparisons, independent score checks, positive interior penalties, solver limits and analytic moment calculations. [v8 complete checks](../../results/observed-history-corrected/ensemble-stability-v8-weight-cv/README.md) add CV for Landweber's conditioning-weight ridge, with the same strict shared training-only splits, data and other libraries. They finished with eight estimates and no final estimator errors, but three final Landweber tolerance misses and worse dose bridge error; this configuration is not adopted. [Exact Gaussian score checks](../../results/observed-history-corrected/bridge-exact-score-check-v1/REPORT.md) separately preserve all bridge fits, raw score audits, weights and unfavorable results. No revised full study is launched.

Frozen sources and old audit helpers are preserved. Separate current postprocessing scripts implement the actual cell-scored checks; no failed criterion is loosened. No known-function or zero-penalty performance study, generating restrictions or Gaussian mean update is introduced. True population values enter only diagnostics after fitting.

# Running study configuration and archived earlier revisions

# Current ensemble stability revision

Current packages: cmbridge 0.3.0.9020 and lmtp 1.6.0.9021. The two saved
n = 4,000 checks are in `ensemble-stability-v6`. Fixed main effects in the
sieve and Landweber classes contain a valid bridge solution on the full
support; the sieve additionally includes joint-category deviations.
All bridge candidates use inverse-expit; all adjoints remain unrestricted.
Gaussian-kernel U-statistics are used consistently for penalties and ensemble
weights. Failed numerical trials and candidates are recorded explicitly.
The approved repeated study uses both mechanisms, n = 500, 1,000, and 4,000,
and 200 replications per combination. It starts only after the complete
SDR/TMLE checks pass. See the frozen source and status in its own result
directory; older stopped studies below remain historical.

# Previous sieve penalty and Landweber class revision

The current packages are cmbridge 0.3.0.9017 and lmtp 1.6.0.9018.
The current rerun uses only the saved n = 4,000 binary and numerical-dose
datasets, with SDR and TMLE. [Report and plots](../../results/observed-history-corrected/sieve-cv-quadratic-landweber-v3/REPORT.md).
Bridge and adjoint sieves select positive ridge penalties by training-only CV,
using shared person splits and raw cell U-statistics. The initial 41-point grid
spans 1e-10 to 1 and boundary minima extend the grid or fail. Singleton-only
subsets return a zero score and are marked uninformative. Optimizer precision
was repaired for a saved nested sieve fit; convergence criteria are unchanged.
The simulation Landweber bridge uses a fixed complete quadratic target basis
(all main effects, squares and pairwise interactions) with joint-category
conditioning indicators. The full-support membership audit is
`audit_landweber_quadratic_class.R`. Previous results and plots are frozen
with their original packages. The full repeated-sample study remains stopped.

# Previous bridge parameterization

The previous [all-candidate report](../../results/observed-history-corrected/selection-inverse-expit-all-bridges-v4/REPORT.md) uses cmbridge 0.3.0.9014 and lmtp 1.6.0.9015.
All three bridge candidates fit inverse-expit with unrestricted coefficients.
All three adjoint candidates use identity. The four isolated datasets retain
the original data, seeds, outer splits and learner-validation splits.
The full repeated-sample study remains stopped. Historical descriptions below
are retained with their original package versions and designs.

# Observed-history simulation study

## Current inverse-expit sieve revision (October 5, 2026)

The bridge sieve optimizes unrestricted real coefficients through `1 / expit`.
The adjoint sieve uses its identity parameterization; Landweber and PMMR retain
their original fits. Therefore the bridge ensemble is not guaranteed to be at
least one. Both conditional-moment libraries exclude saturated L1. The current
packages are cmbridge 0.3.0.9010 and lmtp 1.6.0.9011.

All four isolated datasets completed successfully for SDR and logistic TMLE
with the one-step correction. The saved n = 4,000 datasets, corresponding
n = 20,000 datasets, and all sample assignments match the preceding revision.
[Results, plots, weights and audits](../../results/observed-history-corrected/selection-inverse-expit-final/REPORT.md)
are stored separately. The full repeated-sample study remains stopped.

## Previous three-candidate revision (October 5, 2026)

Both bridge and adjoint libraries now contain saturated joint-category sieve_md,
Landweber and PMMR, with no saturated_l1 candidate. Saturated L1 remains in
regression and classification. lmtp is 1.6.0.9007; cmbridge remains 0.3.0.9006.
Four isolated checks reuse the exact datasets and assignments from
`selection-revision` and save their outputs in `selection-no-saturated-l1`.
The full repeated-sample study remains stopped. The earlier four-candidate
bridge results and their plots are preserved separately.

## Previous selection revision (October 5, 2026)

The repeated study remains stopped at 118 datasets in `missing-both-study-v2`.
Its frozen packages and results are preserved. That revision uses
cmbridge 0.3.0.9006 and lmtp 1.6.0.9006. Every cell validation score now excludes
self-products; ensemble weights use a PSD projection, and raw scores are retained.
Treatment remains in the conditioning cells and saturated regression designs.
The adjoint library has sieve_md, Landweber and PMMR, without saturated L1.
The initial positive penalty grid has 21 values from 1e-5 to 1; boundary choices
are recorded and extend the grid until the minimum is interior, or fail after
16 extensions. Learners fit unbounded. The sequential estimator retains the
logistic TMLE and uses the wider [-100, 100] interval if truncation is needed.

`check_selection_revision.R` runs an isolated dataset specified by mechanism
and sample size. It reuses the original replication-6 seed and, at n=4000, the
saved outer and learner splits. Its four requested checks use n=4000 and
n=20000 for both mechanisms. Results, failures, penalty grids, and raw function
plots are saved separately under `selection-revision`. The full simulation
is not resumed.

The [selection revision report](../../results/observed-history-corrected/selection-revision/REPORT.md)
includes all three fitted-versus-true panels at each size, weights, penalty
extensions, conditional equations, fitted-support checks and full-estimator
failures. `plot_selection_revision.R`, `audit_selection_revision.R` and
`report_selection_revision.py` rebuild these outputs from the saved fits.

## Previous revision: missingness at both follow-ups

The earlier repeated-sample study is stopped, with 176 datasets and one recorded
SDR failure preserved in `ensemble-study`. The full revised study is authorized
for launch in a separate `missing-both-study-v2` directory. It makes measurement
informative at both t = 2 and t = 3. Its constructor is
`make_mechanism(name, visits = "missing_both", dose_max = 3L)`.
Baseline is observed, with the baseline outcome equal to the baseline binary
characteristic. If the intermediate visit is missed, both treatment components
are zero and the final outcome equals that observed baseline outcome; the final
covariate remains random. This makes the target loading satisfy the paper's
adjoint equation at the deterministic-treatment histories. Changing only the
measurement indicator would fail that condition.

Expected missingness at t = 2 and t = 3 is 54.3% and 53.5% for binary treatment,
and 56.8% and 54.3% for numerical dose. Both bridge and adjoint functions are
estimated at both times. The exact checks verify both conditional equations,
the two identifying representations, supported policies, and containment of
valid solutions in the saturated joint-category classes. At missed visits,
bridge solutions need not be unique. The independent true-value calculations
agree within 1e-12. See the current [simulation plan](../../reports/observed-history-simulation.md).

Run the exact audit from the cmbridge-tests directory:

```sh
R_LIBS=results/observed-history-corrected/ensemble-study/library \
  Rscript simulations/observed_history_corrected/audit_missing_both.R
```

The audit generates no sample and fits no estimator. Its outputs are in
`results/observed-history-corrected/missing-both-design`. The current reduced
runner explicitly uses this revision and a separate `missing-both-study-v2`
output directory, retaining the approved sample sizes and replication counts.
The old timing estimate does not apply to this revision. `freeze_ensemble.py`
defaults to that separate directory; freezing sources does not start a run.
The initial `missing-both-study` attempt was interrupted for a zero-row spline
prediction repair. cmbridge 0.3.0.9005 handles empty measured validation groups;
all 115 package tests pass. The DGP, learners and split rules are unchanged.
The scheduled constructor default and earlier results remain available for
historical reproduction. The sections below describe those historical checks.

The separate [n = 10,000 check](../../results/observed-history-corrected/missing-both-n10000/REPORT.md)
has now completed for both mechanisms and both SDR and TMLE, with all four
function classes correctly specified. It retains person-level contributions,
component means and variances, covariance, and exact population remainders.
Its frozen `run_missing_both_quick.R` launches only those two datasets.
It does not resume the repeated-sample study.

## Historical redesign and checks

The study is being redesigned to impose the paper's observed-history treatment
restriction (Assumption 1), as well as the assumptions used by its identification
and estimation results. The previous point-treatment and scheduled-intermediate
designs did not have the identified no-visit treatment violation. The previous
main two-time designs did: they allowed random treatment at a missed visit.

The manuscript explicitly says Assumption 1 is not strictly required for the
technical results in Section 3; Theorem 2 lists Assumptions 2–5 and 7. Thus an
Assumption 1 violation alone does **not** establish a failure of that theorem or
explain the previous dose bias. The replacement study will nonetheless respect
the requested visit-dependent treatment setting.

## Candidate mechanisms and audit

`dgp.R` implements three finite-support candidate mechanisms using components
of C and A rather than separate report notation. At a missed treatment visit,
A is the known deterministic zero prescription. The policy leaves it unchanged.
At a measured treatment visit, the policy sets the first component to one for
binary treatment or increases the dose by one, capped at two. The second
component is retained.

The candidate family with scheduled treatment visits sets intermediate
measurement to one; final-outcome measurement remains informative. The
one-time mechanism has informative measurement of its single follow-up outcome.
The two-time candidates have treatment-confounder feedback through the measured
intermediate covariate and outcome, but do not study intermittent missing health
before later treatment decisions. This is a limitation of that candidate family.

| Paper condition | Verification for scheduled-treatment candidates |
|:--|:--|
| Assumption 1: observed-history treatment rule | Treatment is generated from recorded history; the no-visit branch is deterministic and the policy leaves it unchanged. All treatment-time health is observed in these candidate samples. |
| Assumption 2: sequential randomization | Treatment, subsequent health, and measurement use independent random draws; treatment depends only on recorded history. This is established by construction, not a test of simulated correlations. |
| Assumption 3: supported intervention | Exact natural and policy probabilities are checked on reachable states. The policy does not assign a treatment outside conditional observed support. |
| Assumption 4: shadow-variable exclusion | The measurement probability is a function only of history and next complete health, without a direct current-treatment term. Independent measurement draws establish the conditional exclusion. |
| Assumption 5: positive measurement probability | Every reachable complete-health configuration has strictly positive measurement probability; intermediate scheduled measurement may equal one. |
| Assumption 7: adjoint range | On every reachable baseline/current-state block, the conditional adjoint equation is solved and its residual checked. The health transitions depend only on that state, so the same block applies to each full recorded history with that state; the preceding cumulative ratio is a history-measurable scalar multiplier. |
| Existence of treatment ratios | All variables are discrete. Policy probability mass is calculated by summing over the treatment values mapped to each destination, and absolute continuity is checked directly. The continuous-policy sufficient condition in Assumption 6 is not used. |

The proposed scheduled family passes these finite-support checks. The exact
adjoint-equation residual is below 1e-10 in every audited block. An initial
1,000-person generation pilot was run for each mechanism to check that
measurement before treatment is indeed scheduled. Independent forward and
backward calculations of the true parameter also agree within 1e-12.

Changing only the treatment rule while retaining the previous intermittent
health mechanism is **not** sufficient. At a no-visit history with deterministic
treatment, the adjoint can only produce a conditional constant. The candidate
intermittent mechanism retains a variable next outcome, and its target loading
does not lie in that range. Its exact Assumption 7 check fails. It is rejected
before even a generation pilot or estimator simulation can run.

More specifically, at a no-visit history with deterministic treatment and a
supported policy leaving it unchanged, the current treatment ratio is one.
The adjoint range contains only functions constant in next health at that fixed
history. The paper's loading is the preceding cumulative ratio times the next
outcome. At a history with positive preceding policy weight and a variable next
outcome, it cannot be represented by that adjoint. An intermittent replacement
therefore needs additional structure satisfying this target-specific condition;
it cannot be obtained by simply imposing deterministic treatment on the old
health equations.

## Current outputs

- `results/observed-history-corrected/assumption_audit.csv`: reachable-block audit
  for the scheduled-treatment candidates.
- `results/observed-history-corrected/truth.csv`: independent forward target
  calculations for these candidates.
- `results/observed-history-corrected/rejected_intermittent_design.csv`: the failed
  range checks for the intermittent candidate with deterministic treatment.

Run from the cmbridge-tests directory:

```sh
Rscript simulations/observed_history_corrected/audit.R
```

## Approved reduced study

The current report specifies only the two-time binary and numerical dose
mechanisms. Each uses n = 500, 1,000, and 4,000, with 200 replications per size,
both SDR and TMLE, and all nuisance-function classes correctly specified.
There are 1,200 datasets. The numerical dose now takes values 0, 1, 2, and 3;
the policy increases it by one, capped at 3. The outcome effect per dose unit
is retained. Its new fourth probability weight continues the coefficient trend
from the preceding dose categories.

`make_mechanism("discrete_dose", dose_max = 3)` constructs this revised
mechanism. The default `dose_max = 2` preserves the historical design for
reproduction. `revised_audit.R` verifies both revised mechanisms by deterministic
summation, without generating samples or fitting estimators. Its outputs are
in `results/observed-history-corrected/revised-design`:

```sh
Rscript simulations/observed_history_corrected/revised_audit.R
```

The earlier 10,500-dataset, all-16-specification study is stopped after 19
completed datasets. Its frozen sources and results are preserved separately.
The existing full-study driver describes that earlier design; it is not the
runner for the revised plan. `run_reduced.R` implements the approved reduced study; `freeze_reduced.R` saves its source, validated package archives, and private library before execution.
No known-function or zero-penalty performance simulation is proposed, and
correct learners retain the entire observed history.

## Estimator study

`study.R` supports shared estimator and learner splits. The reduced study uses
three outer groups and three measurement-balanced learner-validation groups.
Additional assignments within candidate training samples are shared between
packages. The historical driver defaults to two learner groups. Saturated L1 candidates select a positive penalty scale by further
cross-validation inside their training samples. For sequential outcome
regressions, SuperLearner combines `SL.saturated_l1_cv`, MARS (`SL.earth`),
and `SL.mean`. Treatment-ratio classification retains `SL.saturated_l1_cv`
and `SL.mean`; MARS is not used for treatment ratios or conditional moments.
The correctly specified conditional-moment library now combines `sieve_md`,
`landweber`, `pmmr`, and `saturated_l1`, with weights chosen by shared-fold
conditional-moment validation. `sieve_md` uses saturated joint-category target
and conditioning bases. Landweber and PMMR retain their default bases;
all three retain their default regularization. Saturated sieve and L1
candidates supply the full discrete target class. The previous single-candidate results are preserved
and do not validate this expanded library. The wrong
bridge and adjoint restrict the target to a constant, while retaining the
full conditioning variables in the moment criterion. Wrong ratio models
omit history; wrong regressions use a separate mean at each regression and
outcome time. In the historical scheduled design, the intermediate outcome
is fully measured, beta is one, and its bridge correction is zero. This
exception does not apply to the current revision.

The adjoint always uses additional correctly specified ratio fits, including
comparisons with wrong ratios in the sequential estimator. Candidate loading
responses are regenerated inside their training/validation samples. The
same person-based penalty assignments are used in both packages when their
training samples coincide. The original learner assignments are preserved.
Further groups inside a training sample alternate person IDs within measurement
profiles, keeping repeated observations together. These assignments use no
health outcomes or observations outside that training sample. Sequential
regression responses are rebuilt separately for each learner and penalty
validation trial: upstream fits for its training responses use only its
training people, with further out-of-sample predictions within that sample.
Responses for its validation people use upstream fits restricted to its
training people. The same construction applies through the backward recursion.

Each checkpoint includes estimates, standard errors, pointwise and joint
intervals, separate population remainder contributions, moment errors,
selected penalties and weights, cell counts, and assignments. Fitted functions
are stored as their values on the entire finite support: these tables exactly
represent their predictions for this study's population. Exact distribution
calculations are diagnostics, not inputs to fitting. Joint intervals use
2,000 Gaussian multiplier draws; drawing from their exact conditional
covariance avoids allocating a person-by-draw matrix.

The package output is checked against a direct evaluation of equation (7) in
every successful comparison. The population bridge remainder is also checked
against its exact product identity. A failed SuperLearner candidate is a
failure, rather than a successful fit with the candidate silently removed.
Completed jobs are saved atomically and source fingerprints prevent mixing
estimator versions on resumption. Partial results and errors are retained.

From the cmbridge-tests directory, after installing the package archives into
`../lmtp/.library-corrected`:

```sh
R_LIBS=../lmtp/.library-corrected:../lmtp/.library \
  Rscript simulations/observed_history_corrected/run.R
Rscript simulations/observed_history_corrected/summarize.R
```

The `run.R` defaults implement the stopped larger design, not the reduced study. `SIM_WORKERS` controls
parallel workers; `SIM_OUTPUT` supplies a separate output directory for pilots
and benchmarks. `SIM_REPS`, `SIM_SIZES`, and `SIM_MECHANISMS` are available for
implementation checks; pilot results are kept separate from the main study.
Seeds depend on mechanism, sample size, and replication, not worker count or
run order. Checkpoints are processed across mechanisms and sizes before moving
to the next replication. The archived previous numerical study is preserved
and is not relabeled as results from this design.

## Execution of the reduced study

The six final timing checks are saved in
`results/observed-history-corrected/reduced-pilot-v2`. All fits completed and
the package output agreed with direct prediction-based calculations. The
earlier failed small-sample check is preserved in `reduced-pilot-v1`.

The previous single-candidate study uses the frozen `reduced-study/run-study.sh`, with eight R
workers and one numerical-library thread per worker. Do not change its frozen
sources or private library while it is running. Completed jobs are atomic
checkpoints; rerunning that launcher resumes unfinished jobs without replacing
completed fits. Its `design.csv` must contain 1,200 rows.

That frozen launcher retains the previous library. It must not be used to
launch the expanded four-candidate cmbridge and MARS configuration.
`freeze_ensemble.py --name ensemble-study` created that separate source/library
freeze and launcher. Its current default prepares `missing-both-study` for
the revised distributions. Earlier attempts are preserved separately.

Tables, graphs, numerical diagnostics, and status are refreshed after the first
batch, periodically thereafter, and at completion. `REPORT.md` labels interim
results as preliminary. The dose timing checks had substantial estimation
error; all prescribed replications will be retained so that the report can
evaluate this behavior without changing the generating distributions or
learners to obtain a desired answer.

## Numerical dose component check

`debug_dose_cmbridge.R` reproduces the saved n = 4,000 dataset's nuisance
functions (seed 5203002), ranks their implied parameter errors, and fits
the same learners in a separate 1,000,000-person sample (seed 5204001).
Only training sample 2 is fitted in that large check, using 666,667 people.
It preserves positive cross-validated penalties and estimated treatment ratios.
The final SDR and TMLE estimators are not fitted in the large sample.

The completed [report](../../results/observed-history-corrected/dose-cmbridge-debug/REPORT.md)
contains the bridge y = x plots and adjoint conditional-equation plots.
`audit_dose_debug.R` verifies both equations, parameter representations,
operator ranks, and the assumptions on every reachable full observed history.
`evaluate_dose_sample_equations.R` evaluates the equations on the original
n = 4,000 observations using predictions from their saved estimator folds.
`plot_dose_cmbridge_debug.R` regenerates figures from saved predictions.
`audit_dose_scoring.R` documents a corrected integer-count overflow in the
large bridge's training-score diagnostic; it does not affect these fitted
functions or the n = 4,000 results.

`audit_dose_function_classes.R` and `audit_dose_sequential_classes.R` check
representability of true functions and valid solutions. They distinguish the
full saturated classes from finite training tables with empty combinations.
`diagnose_large_bridge_points.R` traces distant points to coefficients,
training counts, and conditional-moment matrices. The positive-penalty
comparison in `check_large_bridge_shrinkage.R` uses the original three
candidate scales and the same large training data; no zero penalty or known
function is supplied to a learner. Twelve of sixteen originally distant
points improve with the smaller penalty, while total weighted RMSE increases.

The output folder preserves the execution scripts, numerical records, package
test log, and checksums. Its `REPRODUCE.md` gives commands and source locations.
The expanded repeated-sample study was authorized after this component check
and the expanded-library checks below.

## Expanded-library checks

The expanded [n = 4,000 report](../../results/observed-history-corrected/dose-ensemble-r002/REPORT.md)
reuses seed 5203002 and the original assignments. SDR and TMLE use all four
cmbridge candidates, with joint-category saturation in sieve, and MARS only
for sequential outcome regressions. `rerun_dose_ensemble_r002.R` runs the fits;
`report_dose_ensemble_r002.R` writes the comparisons and plots.

The [million-person report](../../results/observed-history-corrected/dose-ensemble-large/REPORT.md)
reuses seed 5204001 and every original assignment, fitting training sample 2
(666,667 people). `check_large_dose_ensemble.R` estimates all four bridge and
adjoint candidates with estimated out-of-sample treatment ratios;
`audit_large_dose_ensemble.R` checks function-class containment and identifies
the training counts behind large point errors; `report_large_dose_ensemble.R`
produces per-candidate and ensemble y = x graphs. The adjoint graphs compare
the two sides of its conditional equation because its solution is not unique.
This large check does not fit the final SDR/TMLE estimators or outcome
regressions. Each output directory preserves source, diagnostics, plots,
logs, and reproduction instructions.

Before the expanded repeated-sample study, the n = 500 checks revealed a
spline extrapolation defect and overly small training samples in deeper
validation. cmbridge 0.3.0.9004 corrects the spline case and allows bridge
validation groups containing only unmeasured outcomes: their residuals are
defined. lmtp 1.6.0.9005 uses three additional groups when a current training
sample has at most 16 measured people at an incompletely measured outcome
time, and two additional groups otherwise. These assignments use only
training-sample measurement indicators and remain common to both packages.
The primary estimator and learner-validation groups remain three.
`ensemble-pilot-v4` checks both mechanisms and all three sample sizes using
two replications each and eight workers. The earlier configuration failures
and timing trial are retained in separate result directories.

All 12 pilot datasets completed successfully for both SDR and TMLE. Concurrent
processing times averaged 161–213 seconds by mechanism/sample-size combination,
supporting an estimate of eight to nine hours for fitting and analysis. The
1,200-dataset study launched on October 5, 2026 with eight workers from the
separate frozen `ensemble-study/run-study.sh`. Its sources, package archives,
and installed private-library checksums were verified before launch.
It was stopped on October 5, 2026 after 176 saved datasets, at the user's
request. One SDR failure is recorded. These results use the completely
observed intermediate outcome and do not validate the revised distributions.
See its [reproduction instructions](../../results/observed-history-corrected/ensemble-study/REPRODUCE.md).
