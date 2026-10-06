"""Freeze the approved repeated study from the validated package archives."""
import argparse
import hashlib
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tarfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--checks', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
parser.add_argument('--workers', type=int, default=10)
args = parser.parse_args()
project = Path(__file__).resolve().parents[3]
tests = project / 'cmbridge-tests'
out = args.output.resolve()
checks = args.checks.resolve()
if out.exists():
    raise SystemExit('Preserve existing runs; choose a new output directory.')
if not 1 <= args.workers <= 10:
    raise SystemExit('This host has ten logical CPUs.')
for name in ('source', 'packages', 'library', 'validation'):
    (out / name).mkdir(parents=True)
for name, version in (('cmbridge', '0.3.0.9020'), ('lmtp', '1.6.0.9021')):
    archive = checks / 'packages' / f'{name}_{version}.tar.gz'
    shutil.copy2(archive, out / 'packages' / archive.name)
    with tarfile.open(archive) as stream:
        # Support the host's Python 3.9 while rejecting traversal and links.
        for member in stream.getmembers():
            destination = (out / 'source' / member.name).resolve()
            if not str(destination).startswith(str(out / 'source') + os.sep) or member.issym() or member.islnk():
                raise ValueError('Unsafe archive member: ' + member.name)
        stream.extractall(out / 'source')
scripts = out / 'source/simulations/observed_history_corrected'
shutil.copytree(tests / 'simulations/observed_history_corrected', scripts)
(out / 'source/reports').mkdir()
shutil.copy2(tests / 'reports/observed-history-simulation.md', out / 'source/reports/simulation-plan.md')
for name in ('cmbridge-tests.log', 'lmtp-tests.log', 'function-class.log'):
    shutil.copy2(checks / 'validation' / name, out / 'validation' / name)
shutil.copytree(tests / 'results/observed-history-corrected/missing-both-design', out / 'validation/dgp')
fallbacks = [project / 'lmtp/.library-ensemble', tests / 'results/observed-history-corrected/missing-both-study-v2/library']
env = os.environ.copy()
env['R_LIBS'] = ':'.join(map(str, [out / 'library', *fallbacks]))
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    env[key] = '1'
for archive in sorted((out / 'packages').glob('*.tar.gz')):
    with (out / 'validation' / (archive.name.split('_')[0] + '-install.log')).open('w') as log:
        subprocess.run(['R', 'CMD', 'INSTALL', '--library=' + str(out / 'library'), str(archive)],
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
with (out / 'validation/library-versions.log').open('w') as log:
    subprocess.run(['Rscript', '-e', 'stopifnot(packageVersion("cmbridge")=="0.3.0.9020", '
                    'packageVersion("lmtp")=="1.6.0.9021", requireNamespace("earth",quietly=TRUE)); '
                    'print(.libPaths()); print(sessionInfo())'], env=env, stdout=log,
                   stderr=subprocess.STDOUT, check=True)
for directory, filename in [('source', 'SOURCE-SHA256SUMS'), ('packages', 'PACKAGE-SHA256SUMS')]:
    with (out / filename).open('w') as log:
        for path in sorted((out / directory).rglob('*')):
            if path.is_file():
                log.write(hashlib.sha256(path.read_bytes()).hexdigest() + '  ' + str(path.relative_to(out)) + '\n')
q = shlex.quote
lines = ['#!/bin/sh', 'set -eu', 'cd ' + q(str(tests)),
         'export R_LIBS=' + q(env['R_LIBS']), 'export STUDY_SOURCE=' + q(str(scripts)),
         'export LMTP_SOURCE=' + q(str(out / 'source/lmtp')),
         'export CMBRIDGE_SOURCE=' + q(str(out / 'source/cmbridge')),
         'export SIM_OUTPUT=' + q(str(out)), f'export SIM_WORKERS={args.workers}',
         'export SIM_REPS=200 SIM_SIZES=500,1000,4000 SIM_POSTPROCESS=1',
         'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1',
         'exec caffeinate -i Rscript ' + q(str(scripts / 'run_reduced.R'))]
(out / 'run-study.sh').write_text('\n'.join(lines) + '\n')
(out / 'run-study.sh').chmod(0o755)
(out / 'STATUS.txt').write_text('Frozen approved design: 1,200 datasets; simulation has not started.\n')
print(out)
