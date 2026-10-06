# Frozen ensemble study on additional workers

The approved 1,200-dataset study is running on GitHub. Its fitting sources and cmbridge/lmtp archives match the frozen study; the prepared Linux dependency library is archived and checksummed. Linux glmnet and SuperLearner are newer than their Mac counterparts, and all version differences are recorded. Linux and Mac checkpoints are never pooled.

[Actual progress, timings and qualifications](early-summary/CLOUD-PROGRESS.md) accompany [preliminary tables and scientific graphs](early-summary/REPORT.md). The current early snapshot has 24 complete datasets, all binary treatment, and zero final estimator errors. It is not the final study. Fourteen first-wave batches reached their 55-minute early-stage limit without saving a dataset; their logs are preserved and the longer stages are running. These unfinished attempts make full completion within twelve hours doubtful.

[Execution-state snapshot](early-stage-outcomes.csv), [retained early logs](execution-logs), [matching-seed platform comparisons](platform-comparison), and [passing Linux preparation](preparation-metadata) retain the evidence. GitHub normalizes early step conclusions under continue-on-error; exact execution outcomes are recorded separately in the final shard artifacts.

The [adjoint class and defining-equation audit](../full-ensemble-study-v2/diagnostic-evaluation/REPORT.md) qualifies the previous claim that every realized finite-sample class contains a valid solution. The complete saturated class does, while the dictionary learned from measured training combinations can fail this condition.

The original pinned workflow has a reporting-only fingerprint-comparison error. The corrected standalone report workflow will validate all 1,200 final checkpoints; successful reporting must not be confused with estimator performance or with the original workflow's outcome.

[Live computation](https://github.com/idiazst/cmbridge-tests/actions/runs/37417223244). No incomplete replications are presented as final results.
