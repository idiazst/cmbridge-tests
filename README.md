# cmbridge-tests: migrated

Package algorithms, unit tests and numerical recovery checks now live in
[cmbridge](https://github.com/idiazst/cmbridge), under `tests/testthat` and
`tests/validation`. The package workflow publishes recovery figures as artifacts.

The self-censoring implementation is the
[self-censoring branch of lmtp](https://github.com/idiazst/lmtp/tree/self-censoring).

The latest longitudinal simulation and multiple-imputation comparison are in
[lmtp_selfcensor_sim](https://github.com/idiazst/lmtp_selfcensor_sim).
Earlier simulation code and results have been removed from this repository's
current tree. Historical commits remain accessible in Git history.
