# Ensemble checks: current status

Both outer n = 4,000 bridge ensembles now cluster around y = x. Complete SDR/TMLE checks are still running; no final estimates are reported here yet. The full 1,200-dataset study starts after these complete checks pass.

![Binary treatment](binary_longitudinal-bridge-ensemble.png)

![Numerical dose](discrete_dose-bridge-ensemble.png)

The untruncated plots show every reachable combination at time 3 with R₂ = 1, where the bridge is unique. Probability-weighted ensemble RMSE is 0.20639 for binary treatment and 0.24602 for numerical dose. All three bridge candidates use inverse-expit; adjoints remain unrestricted. Sieve and Landweber fixed main-effect classes contain a valid bridge solution on the full support. All histories and treatment predictors are retained.

Packages: cmbridge 0.3.0.9020 and lmtp 1.6.0.9021. Positive sieve and saturated L1 penalties use training-only CV; failed trials and candidates are recorded, and unresolved boundary minima fail. Gaussian kernel U-statistics remove all self-products, with PSD projection only for ensemble weights. Vectorized penalty paths are retained. The accepted DGP and strict response rebuilding are unchanged.

[GitHub validation run 37412467587](https://github.com/idiazst/cmbridge-tests/actions/runs/37412467587) passed every recorded check for package commit 4ad2f90870b4c3282dbaddc67dc2a0d61a7ec9f5. These unit and recovery tests do not replace complete simulation validation.

Earlier unsuccessful and interrupted attempts are preserved separately. The final report will include actual estimates, mean EIF components, all learner weights, errors, and nested fitting failures. Follow the run status for measured timing and progress.
