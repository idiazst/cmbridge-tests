# cmbridge tests

Independent large-sample validation of the bridge and adjoint learners in
[`idiazst/cmbridge-validation`](https://github.com/idiazst/cmbridge-validation).

Each GitHub Actions run:

1. installs a fresh R release on Ubuntu;
2. clones and installs the current `cmbridge-validation` package;
3. records the exact package commit SHA;
4. runs known-truth large-sample tests for the bridge equation;
5. runs known-truth large-sample tests for the adjoint equation;
6. stores CSV results, seed-1 truth-vs-estimate curves, PDFs, and `sessionInfo()`;
7. commits the results under `results/run-<github-run-id>/`.

The validation workflow is manual (`workflow_dispatch`) so commits containing
results do not recursively launch new runs.

## Run again

From a terminal:

```bash
gh workflow run validate.yml --repo idiazst/cmbridge-tests
```

Then inspect the newest run:

```bash
gh run list --repo idiazst/cmbridge-tests --workflow validate.yml --limit 5
```
