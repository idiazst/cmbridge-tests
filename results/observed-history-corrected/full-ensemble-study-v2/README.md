# Full observed-history simulation

The approved study started at 00:28 EDT on October 6, 2026. It uses both mechanisms at n = 500, 1,000, and 4,000, with 200 replications each (1,200 datasets), ten workers, SDR and logistic TMLE, the one-step bridge correction, and shared training-only splits. The two saved n = 4,000 checks and all audits passed before launch.

See [progress](PROGRESS.txt), [design](design.csv), and [the completed check report](../ensemble-stability-v6/REPORT.md). Repeated-sample results are not yet available. The fitting sources and package versions are frozen separately from earlier studies.

The cached-outer n = 4,000 checks took 17.5 and 18.7 minutes with two concurrent processes. These are not full-study timings: the full study refits its outer nuisance functions and uses ten concurrent workers. The first completed datasets will provide a measured runtime estimate. Completion within 12 hours is not established.
