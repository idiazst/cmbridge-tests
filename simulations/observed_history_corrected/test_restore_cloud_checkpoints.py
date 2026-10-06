"""Integrity tests only: no simulated estimates or nuisance fits are run."""
import csv
import tempfile
import unittest
from pathlib import Path

from restore_cloud_checkpoints import restore


class RestoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.incoming = self.root / "incoming"
        self.incoming.mkdir()
        self.output = self.root / "restored"

    def artifact(self, name, jobs, fingerprint="test-fingerprint"):
        path = self.incoming / name
        path.joinpath("jobs").mkdir(parents=True)
        with path.joinpath("design.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["job", "shard"])
            writer.writeheader()
            writer.writerows({"job": f"job-{i}", "shard": i % 40 + 1}
                             for i in range(1200))
        path.joinpath("specifications.csv").write_text("all four classes\n")
        path.joinpath("design-details.rds").write_bytes(b"unchanged-design")
        path.joinpath("source-fingerprint.txt").write_text(fingerprint + "\n")
        path.joinpath("run.log").write_text(f"Original execution {name}\n")
        for job, contents in jobs.items():
            path.joinpath("jobs", f"{job}.rds").write_bytes(contents)
        return path

    def test_identical_duplicates_and_idempotence(self):
        self.artifact("cloud-shard-111-1", {"job-0": b"old"})
        self.artifact("cloud-resume-222-1", {"job-0": b"old", "job-40": b"new"})
        self.assertEqual(restore(self.incoming, self.output, 1), 2)
        self.assertEqual(restore(self.incoming, self.output, 1), 2)
        self.assertEqual(len(list(self.output.joinpath("jobs").glob("*.rds"))), 2)
        self.assertEqual(self.output.joinpath("jobs/job-0.rds").read_bytes(), b"old")
        self.assertTrue(self.output.joinpath("prior-execution/cloud-shard-111-1/run.log").exists())

    def test_conflicting_checkpoint_is_rejected_before_copying(self):
        self.artifact("cloud-shard-111-1", {"job-0": b"old"})
        self.artifact("cloud-resume-222-1", {"job-0": b"changed"})
        with self.assertRaisesRegex(ValueError, "Conflicting copies"):
            restore(self.incoming, self.output, 1)
        self.assertFalse(self.output.exists())

    def test_changed_fingerprint_is_rejected(self):
        self.artifact("cloud-shard-111-1", {"job-0": b"old"})
        self.artifact("cloud-resume-222-1", {"job-40": b"new"}, fingerprint="changed")
        with self.assertRaisesRegex(ValueError, "Conflicting statistical"):
            restore(self.incoming, self.output, 1)

    def test_cross_shard_checkpoint_is_rejected(self):
        self.artifact("cloud-shard-111-1", {"job-1": b"wrong-shard"})
        with self.assertRaisesRegex(ValueError, "Unexpected checkpoint"):
            restore(self.incoming, self.output, 1)

    def test_existing_different_result_is_preserved(self):
        self.artifact("cloud-shard-111-1", {"job-0": b"old"})
        self.output.joinpath("jobs").mkdir(parents=True)
        destination = self.output / "jobs/job-0.rds"
        destination.write_bytes(b"do-not-overwrite")
        with self.assertRaisesRegex(ValueError, "Conflicting existing checkpoint"):
            restore(self.incoming, self.output, 1)
        self.assertEqual(destination.read_bytes(), b"do-not-overwrite")


if __name__ == "__main__":
    unittest.main()
