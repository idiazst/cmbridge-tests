"""Freeze the expanded study without replacing earlier sources or results."""
from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import argparse
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--name', default='missing-both-selection-study', help='New result directory name.')
args = parser.parse_args()
if not re.fullmatch(r'[a-z][a-z0-9-]*', args.name):
    raise SystemExit('Use a lowercase result directory name without path separators.')

project = Path(__file__).resolve().parents[3]
tests = project / "cmbridge-tests"
base = tests / "results/observed-history-corrected"
out = base / args.name
if out.exists():
    raise SystemExit("Result directory already exists. Preserve it; use a new name for a changed design.")
if not (base / 'missing-both-design/assumption_audit.csv').exists():
    raise SystemExit('Run audit_missing_both.R before freezing the revised design.')
for name in ("source", "packages", "library", "validation"):
    (out / name).mkdir(parents=True, exist_ok=True)
shutil.copytree(project / "cmbridge", out / "source/cmbridge",
                ignore=shutil.ignore_patterns(".git", ".DS_Store", "*.Rcheck"))
shutil.copytree(project / "lmtp", out / "source/lmtp",
                ignore=shutil.ignore_patterns(".git", ".library*", ".DS_Store", "*.Rcheck", "results"))
shutil.copytree(tests / "simulations/observed_history_corrected",
                out / "source/simulations/observed_history_corrected")
(out / "source/reports").mkdir()
shutil.copy2(tests / "reports/observed-history-simulation.md", out / "source/reports/simulation-plan.md")
shutil.copytree(base / "reduced-study/library", out / "library", dirs_exist_ok=True)
shutil.copytree(base / "dose-ensemble-large/library", out / "library", dirs_exist_ok=True)
shutil.copytree(base / "missing-both-design", out / "validation/dgp")
shutil.copy2(base / "cmbridge-small-sample-tests.log", out / "validation/cmbridge-small-sample-tests.log")
shutil.copy2(base / "cmbridge-empty-predictions-tests.log", out / "validation/cmbridge-tests.log")
shutil.copy2(base / "lmtp-adaptive-folds-tests.log", out / "validation/lmtp-tests.log")
env = os.environ.copy()
env["R_LIBS"] = str(out / "library")
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    env[key] = "1"
r = shutil.which("R")
rscript = shutil.which("Rscript")
for package, version in (("cmbridge", "0.3.0.9010"), ("lmtp", "1.6.0.9011")):
    with (out / f"validation/{package}-build.log").open("w") as log:
        subprocess.run([r, "CMD", "build", str(out / "source" / package)],
                       cwd=out / "packages", env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    archive = out / "packages" / f"{package}_{version}.tar.gz"
    if not archive.exists():
        raise SystemExit(f"Expected archive missing: {archive}")
    with (out / f"validation/{package}-install.log").open("w") as log:
        subprocess.run([r, "CMD", "INSTALL", "--library=" + str(out / "library"), str(archive)],
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
with (out / "validation/library-versions.log").open("w") as log:
    subprocess.run([rscript, "-e", 'stopifnot(packageVersion("cmbridge")=="0.3.0.9010", '
                    'packageVersion("lmtp")=="1.6.0.9011", requireNamespace("earth",quietly=TRUE)); '
                    'print(.libPaths()); print(sessionInfo())'], env=env, stdout=log,
                   stderr=subprocess.STDOUT, check=True)

def manifest(directory, filename):
    with (out / filename).open("w") as stream:
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                stream.write(hashlib.sha256(path.read_bytes()).hexdigest() + "  " +
                             str(path.relative_to(out)) + "\n")

manifest(out / "source", "SOURCE-SHA256SUMS")
manifest(out / "packages", "PACKAGE-SHA256SUMS")
manifest(out / "library", "LIBRARY-SHA256SUMS")
import shlex
q = shlex.quote
scripts = out / "source/simulations/observed_history_corrected"
launcher = ["#!/bin/sh", "set -eu", "cd " + q(str(tests)),
    "export R_LIBS=" + q(str(out / "library")),
    "export STUDY_SOURCE=" + q(str(scripts)),
    "export LMTP_SOURCE=" + q(str(out / "source/lmtp")),
    "export CMBRIDGE_SOURCE=" + q(str(out / "source/cmbridge")),
    "export SIM_OUTPUT=" + q(str(out)), "export SIM_WORKERS=8", "export SIM_POSTPROCESS=1",
    "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1",
    "unset SIM_REPS SIM_SIZES", q(rscript) + " " + q(str(scripts / "run_reduced.R"))]
(out / "run-study.sh").write_text("\n".join(launcher) + "\n")
(out / "run-study.sh").chmod(0o755)
(out / "STATUS.txt").write_text("Packages, learner settings, and sources frozen for the design with missingness at both follow-ups. No simulation has started.\n")
print("Expanded study frozen at", out)
