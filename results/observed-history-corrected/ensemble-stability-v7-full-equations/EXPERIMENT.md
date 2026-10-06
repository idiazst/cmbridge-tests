# Separate full-conditional-equation checks

These two checks use the already saved binary and numerical-dose n4000 datasets, seeds 5103006 and 5203006, with exactly the previous outer and learner assignments. They are not additional simulation replications and are not pooled with the running cloud study.

The original frozen study remains unchanged. Packages remain cmbridge 0.3.0.9020 and modified lmtp 1.6.0.9021. This separate configuration uses joint-category conditioning for the bridge sieve and Landweber, and the existing cell U-statistic for conditional-moment penalty CV and ensemble selection. The ensemble Gram is projected to PSD before weight optimization. The aim is to check the full defining conditional equations, following the algebraic identification audit. The one-training-sample diagnostic improved with this combination; that result alone is insufficient to adopt it.

Bridge target classes, inverse-expit parameterization, unrestricted adjoints, all history/treatment predictors, positive training-only penalty CV with boundary extension/failure, and convergence criteria are unchanged. No saturated L1 in cmbridge. Regression retains saturated L1, MARS, mean; treatment classification retains saturated L1, mean. Both complete SDR and logistic TMLE have the one-step bridge correction. No known-function fitting, zero penalties, generating restrictions or Gaussian mean update. Truth is diagnostic only.

Do not launch a revised full study before these complete-estimator checks and raw y=x/conditional-equation plots have been examined. Preserve estimator failures and unfavorable comparisons. Do not copy nuisance caches from an older configuration. The check script may resume only its own saved nuisance fits in this separately frozen directory.
