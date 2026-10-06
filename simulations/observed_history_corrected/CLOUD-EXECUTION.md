# Resume and report the approved frozen study

The original study is [GitHub run 37417223244](https://github.com/idiazst/cmbridge-tests/actions/runs/37417223244), pinned to `e96efbc9ec26a0f5bff036c4b4a8c9025173b8b5`. The 1,200 planned datasets, seeds, package archives, fitting code and common Linux dependency library remain unchanged. Keep Mac diagnostic fits separate from Linux study checkpoints.

The initial jobs have explicit time limits and upload saved checkpoints even when their computation stops. Download these artifacts before their one-day retention expires. A saved final estimator failure counts as a completed attempt; retain it for investigation. An interrupted dataset without a final checkpoint is unfinished.

`ensemble-study-resume.yml` restores prior final/partial artifacts for the specified shards, verifies byte equality for duplicate checkpoints, validates their R metadata and seeds, and runs the **unchanged original dispatcher**. It skips saved datasets and computes only unfinished ones. Original execution logs remain under `prior-execution`. All restored source/dependency fingerprints must match. Missing artifacts, expired artifacts and conflicting copies cause an error rather than an unrecorded restart.

Dispatch this only after the selected shards have uploaded final or partial `cloud-shard` artifacts, and after inspecting which shards are unfinished. Standard runners share the existing concurrency limit; this does not create additional paid capacity.

```sh
gh workflow run ensemble-study-resume.yml --repo idiazst/cmbridge-tests --ref main \
  --field study_run=37417223244 \
  --field checkpoint_runs=37417223244 \
  --field shards='[1,2]'
```

The two shard numbers above are illustrative, not a claim that those shards are currently unfinished. A subsequent resume should list the original run and all prior resume run IDs in `checkpoint_runs`, separated by commas. Result artifact names are `cloud-resume-<resume run ID>-<shard>`. Download them into separate directories; do not overwrite original artifacts.

The corrected reporting workflow accepts original and resumed checkpoints. It merges byte-identical duplicates once per shard before aggregation and verifies unique seeds and a single statistical/dependency fingerprint. A final report requires all 200 datasets in each of the six mechanism/sample-size cells. A failed fit yields no interval in the additional coverage table; pending datasets are excluded from completed-dataset denominators.

```sh
gh workflow run ensemble-study-summary.yml --repo idiazst/cmbridge-tests --ref main \
  --field study_run=37417223244 \
  --field checkpoint_runs='<actual resume run IDs, separated by commas>' \
  --field early_checkpoints=false
```

Omit `checkpoint_runs` when no resume has run. `early_checkpoints=true` explicitly creates a preliminary report from the early artifacts. It is not a final 200-replication result. The original pinned workflow has a fingerprint-name reporting bug; use the corrected reporting workflow instead of refitting datasets to change that workflow's status.

Restoration was checked with five automated integrity tests and the 24 actual saved Linux checkpoints. All 24 RDS objects were preserved byte-for-byte. The restored aggregation reproduced `summary.csv`, `replicates.csv`, `coverage-including-failures.csv` and `population-error-by-outcome-time.csv` exactly. This verification ran no nuisance fits.
