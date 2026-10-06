# Frozen study on additional workers

This run uses the exact package archives and fitting source in [the frozen study](../full-ensemble-study-v2/README.md). The execution dispatcher partitions its unchanged 1,200 seeds into 40 disjoint groups of 30 datasets. At most 20 standard Ubuntu 24.04 jobs run simultaneously, with four R workers each. No statistical fitting settings are changed.

The study includes binary and numerical dose mechanisms, n=500,1000,4000, and 200 replications per mechanism and sample size. Both SDR and logistic TMLE receive the one-step bridge correction and the same training-only sample assignments. Known-function performance simulations, zero-penalty simulations and Gaussian mean updates are excluded.

R is fixed at 4.5.2. A preparation job verifies archive/source checksums, installs the exact package versions, tests the packages, and checks that generated n=4000 data and outer/learner splits match the Mac hashes for seeds 5103001 and 5203001. It builds one shared Linux dependency library used by every shard. The relative-path source fingerprint also includes dispatcher code and package versions. Linux and Mac results are kept separately; repeated seeds from the Mac timing attempt are not counted twice.

Each shard uploads an early four-dataset batch for numerical and runtime diagnostics, then completes its 30 datasets. Errors, warnings, failed candidates and failed penalty trials remain in the saved checkpoints. A timeout preserves partial output. The aggregation verifies unique seeds and one source fingerprint before making tables and graphs. Incomplete replications must be described as preliminary.

The [workflow](../../../.github/workflows/ensemble-study.yml) uses standard public-repository runners, with one-day artifact retention. [GitHub documents standard public runner compute as free](https://docs.github.com/en/billing/concepts/product-billing/github-actions); no larger paid runners are used. Checkpoints and the dependency archive are approximately 150 MB based on the first local files, below the free account's 500 MB storage allowance; actual artifact sizes must be monitored. The preparation library and checkpoints are downloaded and verified before their own artifacts are removed.

Status: workflow prepared; no completed repeated-study result is claimed here.
