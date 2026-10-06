"""Describe actual completed cloud work; do not turn timing projections into results."""
import csv
import math
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

root = Path(os.environ['SIM_OUTPUT'])
jobs = list(csv.DictReader((root / 'job_log.csv').open()))
groups = defaultdict(list)
for job in jobs:
    groups[job['mechanism'], int(job['n'])].append(job)
seconds = [float(job['seconds']) for job in jobs]
assert seconds and len(jobs) <= 1200
lines = ['# Cloud study progress', '',
         f"Updated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}. **Partial results: {len(jobs)} of 1,200 planned datasets.**", '',
         '[Live computation](https://github.com/idiazst/cmbridge-tests/actions/runs/37417223244) uses the approved frozen fitting source and exact archived package versions. There are forty shards of thirty datasets, up to twenty simultaneous GitHub jobs and four R workers per job. Saved seeds, sample assignments, source fingerprints and four estimates per successful dataset are validated before aggregation. Mac timing checkpoints are kept separate.', '',
         '| Treatment | n | Completed datasets | Dataset processing time, minutes: mean [minimum, maximum] | Datasets with final errors |',
         '| --- | --- | --- | --- | --- |']
for (mechanism, n), group in sorted(groups.items()):
    times = [float(job['seconds'])/60 for job in group]
    errors = sum(int(job['errors']) > 0 for job in group)
    name = 'Binary' if mechanism == 'binary_longitudinal' else 'Numerical dose'
    lines.append(f'| {name} | {n} | {len(group)} | {sum(times)/len(times):.1f} [{min(times):.1f}, {max(times):.1f}] | {errors} |')
ideal_hours = 1200 * sum(seconds)/len(seconds) / 80 / 3600
lines += ['',
          f'Mean processing time across saved datasets is {sum(seconds)/len(seconds)/60:.1f} minutes. Applying that mean to all 1,200 datasets at eighty fully occupied workers gives {ideal_hours:.1f} hours of idealized computation. This omits setup, scheduling, uneven shard loads and restarts, and is not a measured completion time. Different mechanisms or slow pending datasets can change this estimate. Early completed datasets are not a random timing sample.', '',
          'Bias, confidence-interval coverage and convergence conclusions require the planned 200 replications per cell. Pending datasets are excluded from preliminary calculations; recorded final and nested fitting failures remain in the output. An execution timeout is kept distinct from a saved statistical fitting error.', '',
          '[Actual preliminary tables and scientific graphs](REPORT.md), [job timings](job_log.csv), [sample EIF component means](eif-component-means-by-replication.csv), and [fixed-function population contributions by outcome time](population-error-by-outcome-time.csv) retain the completed observations.', '',
          'The full-support saturated adjoint class contains a valid solution, while a dictionary learned only from measured training combinations can fail to contain an exact solution across the whole population. The [saved Mac class and equation audit](https://github.com/idiazst/cmbridge-tests/blob/main/results/observed-history-corrected/full-ensemble-study-v2/diagnostic-evaluation/REPORT.md) documents this distinction. The current cloud run is preserved so that any later change can be assessed against the original version.', '',
          'The original pinned workflow has a known reporting-only fingerprint comparison error. A corrected standalone reporting workflow will verify and summarize all final checkpoints. This does not affect the running statistical fits. Package tests and source/split checks passed in preparation.']
early_file = root.parent / 'early-stage-outcomes.csv'
if early_file.exists():
    early = list(csv.DictReader(early_file.open()))
    censored = []
    for row in early:
        if row['saved_early_datasets'] != '0' or not row['early_finished']:
            continue
        elapsed = (datetime.fromisoformat(row['early_finished'].replace('Z','+00:00'))
                   - datetime.fromisoformat(row['early_started'].replace('Z','+00:00'))).total_seconds()
        if elapsed >= 3299:
            censored.append(row)
    if censored:
        attempted = len(jobs) + 4 * len(censored)
        lower_mean = (sum(seconds) + 4*len(censored)*3300) / attempted / 60
        lines += ['', '## Early batches that reached their limit', '',
                  f'{len(censored)} of twenty first-wave batches reached the 55-minute limit with none of their four datasets saved. All twenty jobs have moved to the longer computation stage. The logs and [{early_file.name}](../{early_file.name}) preserve these incomplete execution attempts; they are not counted as completed statistical replications.', '',
                  f'Using 55 minutes as a lower bound for each unfinished dataset, the mean time among all {attempted} initially started datasets is at least {lower_mean:.1f} minutes. If that mean represented the full study, even this bound would imply {1200*lower_mean/80/60:.1f} hours at eighty fully occupied workers, before setup, restarts or load imbalance. The actual duration can be longer. The completed-only timing projection above therefore understates the current evidence.', '',
                  'GitHub normalizes the early step conclusion under continue-on-error. A displayed successful step is not proof that its four datasets finished. Final shard artifacts record the actual early and compute outcomes separately. No incomplete replication is presented as final.']
comparison_file = root.parent / 'platform-comparison' / 'dependency-versions.csv'
if comparison_file.exists():
    versions = list(csv.DictReader(comparison_file.open()))
    different = [r for r in versions if r['same_version']=='FALSE']
    lines += ['', '## Cloud and Mac dependency versions', '',
              'The cmbridge and modified lmtp archives match exactly. The prepared Linux environment has newer regression dependencies, so numerical equality across platforms is not claimed. Linux checkpoints share one archived dependency library and one source/version fingerprint.', '',
              '| Package | Mac | Linux |', '| --- | --- | --- |']
    lines += [f"| {r['package']} | {r['local_version']} | {r['cloud_version']} |" for r in different]
    lines += ['', '[Matching-seed comparisons](../platform-comparison/estimate-differences.csv) record the actual differences and exact agreement of sample assignments. These results are never pooled across platforms.']
bridge_file = root / 'saved-cloud-bridge-rmse.csv'
if bridge_file.exists():
    bridge_rows = list(csv.DictReader(bridge_file.open()))
    mechanisms = sorted({r['mechanism'] for r in bridge_rows})
    lines += ['', '## Completed cloud n=4000 bridge plots', '',
              'The plots show every saved ensemble prediction with R2 = 1 and R3 = 1, on equal axes. Some training samples have visibly larger deviations; they are retained. A favorable bridge plot alone does not establish coverage or accurate adjoint estimation.']
    for mechanism in mechanisms:
        lines += ['', f'![{mechanism}: completed n=4000 bridge predictions](figures/{mechanism}-saved-n4000-bridges.png)']
penalties = list(csv.DictReader((root / 'penalty_and_weight_selections.csv').open()))
selected = [r for r in penalties if r['selected']=='TRUE']
assert all(float(r['scale'])>0 and math.isfinite(float(r['loss'])) and r['failed']=='FALSE'
           and r['penalty_status'].startswith('interior') for r in selected)
candidate_records = len(list(csv.DictReader((root / 'nested-candidate-failures.csv').open())))
trial_records = len(list(csv.DictReader((root / 'nested-penalty-trial-failures.csv').open())))
lines += ['', '## Recorded fitting checks', '',
          f'All {len(selected)} selected penalty records have positive penalties, finite successful trials and an interior status, including minima that became interior after grid extension. The output retains {candidate_records} nested candidate-failure records and {trial_records} failed penalty-trial records. These are repeated nested fitting events, not counts of failed datasets; they must not be interpreted as independent failures or omitted from the report.']
(root / 'CLOUD-PROGRESS.md').write_text('\n'.join(lines)+'\n')
print(len(jobs), 'actual cloud checkpoints; idealized computation hours', round(ideal_hours, 2))
