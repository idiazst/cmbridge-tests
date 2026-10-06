"""Wait for completed checks, audit them, then run the approved full study."""
import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--checks', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--check-pids', type=int, nargs=2, required=True)
args = parser.parse_args()
checks = args.checks.resolve(); output = args.output.resolve()
project = Path(__file__).resolve().parents[3]
tests = project / 'cmbridge-tests'; scripts = Path(__file__).resolve().parent
lock = checks / 'continuation-lock'
try:
    lock.mkdir()
except FileExistsError:
    raise SystemExit('Continuation already has a lock; inspect its process and status before restarting.')
status = checks / 'CONTINUATION.txt'
status.write_text('Waiting for complete SDR/TMLE checks and lmtp integration tests. Full study not started.\n')
def alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))
env = os.environ.copy()
env.update(R_LIBS=':'.join(map(str, [checks / 'library', project / 'lmtp/.library-ensemble',
    tests / 'results/observed-history-corrected/missing-both-study-v2/library'])),
    SIM_OUTPUT=str(checks), SIM_REFERENCE=str(tests / 'results/observed-history-corrected/selection-inverse-expit-all-bridges-v4'),
    STUDY_SOURCE=str(checks / 'source/simulation'), SIM_SIZES='4000',
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
cases = [checks / 'binary_longitudinal-n4000-seed5103006', checks / 'discrete_dose-n4000-seed5203006']
try:
    while True:
        ready = []
        for folder, pid in zip(cases, args.check_pids):
            path = folder / 'STATUS.txt'
            if path.exists():
                if 'Estimator error rows 0' not in path.read_text():
                    raise RuntimeError('Estimator check failed: ' + str(folder))
                if len(rows(folder / 'estimates.csv')) != 4:
                    raise RuntimeError('Incomplete estimates: ' + str(folder))
                ready.append(True)
            elif not alive(pid):
                raise RuntimeError('Check process exited without completed outputs: ' + str(folder))
            else:
                ready.append(False)
        log = checks / 'validation/lmtp-tests.log'
        text = log.read_text() if log.exists() else ''
        if 'Error:' in text or '══ Failed' in text:
            raise RuntimeError('lmtp integration tests failed; inspect their log.')
        if all(ready) and '══ DONE' in text:
            break
        time.sleep(15)
    status.write_text('Checks complete. Running all audits before the full study.\n')
    for script in ['plot_selection_revision.R', 'plot_all_bridge_learners.R',
                   'audit_selection_revision.R', 'audit_sieve_cv_rerun.R', 'audit_bridge_library_class.R']:
        with (checks / 'validation' / (script[:-2] + '.log')).open('w') as log:
            subprocess.run(['Rscript', str(scripts / script)], cwd=tests, env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run([sys.executable, str(scripts / 'report_ensemble_stability.py')],
                   cwd=tests, env=env, check=True)
    subprocess.run([sys.executable, str(scripts / 'freeze_stability_study.py'),
                    '--checks', str(checks), '--output', str(output), '--workers', '10'],
                   cwd=tests, env=env, check=True)
    status.write_text('Audits passed. Full study launched in ' + str(output) + '\n')
    with (output / 'run.log').open('w') as log:
        subprocess.run([str(output / 'run-study.sh')], cwd=tests, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    status.write_text('Full study finished. Reports and unexpected-result investigation require final review and publication.\n')
except Exception as error:
    status.write_text('Continuation stopped: ' + str(error) + '\nInspect the preserved logs and outputs.\n')
    raise
