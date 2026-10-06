"""Restore byte-identical checkpoints for one shard without changing fit code.

The input can contain both early and final artifacts. A dataset is copied once;
any conflicting copy is an error. Original logs and execution metadata remain
under prior-execution. R validates seeds, designs, and fingerprints before fits
resume. This script never reads or modifies the serialized R objects.
"""
import argparse
import csv
import hashlib
import re
import shutil
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_copy(source, destination):
    if destination.exists():
        if sha256(source) != sha256(destination):
            raise ValueError(f"Conflicting existing file: {destination}")
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)


def restore(incoming, output, shard):
    if not 1 <= shard <= 40:
        raise ValueError("Shard must be in 1..40")
    artifact_pattern = re.compile(rf"cloud-(?:early|shard|resume)-[0-9]+-{shard}$")
    directories = sorted(path for path in incoming.iterdir()
                         if path.is_dir() and artifact_pattern.fullmatch(path.name))
    if not directories:
        raise ValueError(f"No saved artifact for shard {shard}")
    metadata = ("design.csv", "specifications.csv", "design-details.rds",
                "source-fingerprint.txt")
    fingerprints = {path.joinpath("source-fingerprint.txt").read_text().strip()
                    for path in directories}
    if len(fingerprints) != 1:
        raise ValueError("Conflicting statistical/dependency fingerprints")
    reference = directories[0]
    with reference.joinpath("design.csv").open(newline="") as handle:
        design = list(csv.DictReader(handle))
    if len(design) != 1200 or len({row["job"] for row in design}) != 1200:
        raise ValueError("Invalid 1,200-dataset design")
    assigned = {row["job"] for row in design if int(row["shard"]) == shard}
    if len(assigned) != 30:
        raise ValueError("Invalid 30-dataset shard assignment")
    # Check all conflicts before writing any restored checkpoints.
    jobs = {}
    origins = []
    for directory in directories:
        for name in ("design.csv", "specifications.csv"):
            if sha256(directory / name) != sha256(reference / name):
                raise ValueError(f"Conflicting {name} in {directory.name}")
        for path in sorted(directory.joinpath("jobs").glob("*.rds")):
            if path.stem not in assigned:
                raise ValueError(f"Unexpected checkpoint for shard {shard}: {path.name}")
            digest = sha256(path)
            if path.name in jobs and jobs[path.name][1] != digest:
                raise ValueError(f"Conflicting copies of checkpoint: {path.name}")
            jobs.setdefault(path.name, (path, digest))
            destination = output / "jobs" / path.name
            if destination.exists() and sha256(destination) != digest:
                raise ValueError(f"Conflicting existing checkpoint: {path.name}")
            origins.append({"job": path.stem, "artifact": directory.name,
                            "sha256": digest, "bytes": path.stat().st_size})
    for name in metadata:
        checked_copy(reference / name, output / name)
    for path, _ in jobs.values():
        checked_copy(path, output / "jobs" / path.name)
    for directory in directories:
        # Retain earlier timeouts, warnings, and attempts rather than replacing
        # their logs with a later successful execution.
        for path in sorted(directory.rglob("*")):
            if path.is_file() and "jobs" not in path.relative_to(directory).parts:
                checked_copy(path, output / "prior-execution" / directory.name /
                             path.relative_to(directory))
        for path in directory.glob("started-shard-*.txt"):
            destination = output / path.name
            if not destination.exists():
                checked_copy(path, destination)
    output.mkdir(parents=True, exist_ok=True)
    with output.joinpath("restored-checkpoints.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["job", "artifact", "sha256", "bytes"])
        writer.writeheader()
        writer.writerows(origins)
    output.joinpath("RESTORED.txt").write_text(
        f"Restored {len(jobs)} unique checkpoints of 30 assigned datasets.\n"
        "Duplicate copies were required to be byte-identical.\n"
        "Original execution logs retained under prior-execution.\n"
        "Final fit failures remain saved and are not silently retried.\n")
    return len(jobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--incoming", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--shard", required=True, type=int)
    args = parser.parse_args()
    count = restore(args.incoming, args.output, args.shard)
    print(f"Restored {count} unique checkpoints; R will validate before resuming.")
