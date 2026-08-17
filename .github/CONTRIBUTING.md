# Contributing to STAPLE

Thank you for your interest in contributing! This document covers the conventions and workflow used in this repository. It is a living draft — please open an issue or PR if something is unclear or missing.

---

## Table of Contents

1. [Getting started](#getting-started)
2. [Branching and pull requests](#branching-and-pull-requests)
3. [Development environment](#development-environment)
4. [Testing](#testing)
5. [Lint and formatting](#lint-and-formatting)
6. [Adding a new tool or subworkflow](#adding-a-new-tool-or-subworkflow)
7. [Reporting outputs](#reporting-outputs)
8. [Resource labels](#resource-labels)
9. [Commit style](#commit-style)

---

## Getting started

1. Fork the repository and clone your fork.
2. Create a feature branch off `dev` (not `main`).
3. Open a draft PR early so others can follow your work.
4. Target `dev` with your PR; `main` is used for releases only.

---

## Branching and pull requests

- Branch naming: `feature/<short-description>`, `fix/<short-description>`, `docs/<short-description>`.
- Keep PRs focused. A single logical change per PR is easier to review and bisect.
- Fill in the PR template, including a summary, motivation, and test evidence.
- PRs must pass the CI test suite (`nf-test`) before merge.

---

## Development environment

| Requirement | Version |
|---|---|
| Nextflow | ≥ 25.x |
| nf-test | as pinned in `nf-test.config` |
| Docker | any recent release |
| Python | 3.10+ (samplesheet validator only) |

Run the pipeline locally against the bundled test data with Docker:

```bash
cd tests/data/samplesheets
nextflow run ../../../main.nf \
  --input multisample-test.csv \
  --max_memory 8GB \
  --max_cpus 4 \
  --outdir outs \
  -profile docker
```

---

## Testing

Tests are split into two scopes. Run what is relevant to your change:

```bash
# All process-level tests
nf-test test --tag process

# All pipeline-level tests
nf-test test --tag pipeline

# A single module test file
nf-test test tests/modules/local/rctd/rctd.nf.test

# A single pipeline test file
nf-test test tests/main.nf.test

# A narrow tag slice inside a file
nf-test test tests/main.nf.test --tag visiumhd
nf-test test tests/modules/local/squidpy/main.nf.test --tag ligrec
```

`nf-test.config` points tests at `tests/nextflow.config` and uses the `docker` profile by default.

When adding a new process, add a corresponding test file under `tests/modules/local/<tool>/`. When changing workflow wiring, update or add a pipeline-level test in `tests/main.nf.test`.

---

## Lint and formatting

Nextflow files do not have an auto-formatter at this time; follow the existing indentation and style (I know, it's not perfect). Use `nf-core lint` to check for common issues:

```bash
nf-core lint .
```

## Adding a new tool or subworkflow

1. **Process definition** goes in `modules/local/<tool>/main.nf`. Follow the tuple shape `tuple val(meta), path(...)` used by all existing processes.

2. **Params** live under a nested key such as `params.deconvolve.<tool>` or `params.analyze.<tool>`. Do not add top-level params.

3. **Wire the process** into the appropriate subworkflow under `subworkflows/local/` (e.g. `deconvolve/main.nf` or `analyze/main.nf`), then thread the result channel through to `workflows/spatial.nf`.

4. **Propagate outputs** — new analysis outputs must be attached back to `adata` and/or funnelled into `ch_report` / `ch_multiqc_files`; otherwise they will not appear in MultiQC. See the _Reporting outputs_ section below.

5. **Resource label** — annotate the process with an appropriate label (see _Resource labels_).

6. **Dockerfile** — if the tool requires a new image, add a `Dockerfile` alongside `main.nf` in the module directory.

7. **Test** — add a test file under `tests/modules/local/<tool>/` tagged with `process`, and extend `tests/main.nf.test` for any new end-to-end path.

### Compatibility matrix

When a new tool does not support all input formats, update the compatibility table in [README.md](../README.md).

---

## Reporting outputs

The reporting path is as important as the primary output path. Two things need to happen:

- The output must be attached to the `adata` object (via `STAPLE_ATTACH_*` or equivalent utility step in `modules/local/util/main.nf`).
- The output must be published to `ch_report` / `ch_multiqc_files` so it appears in the MultiQC HTML report.

If only the first step is done, the output exists on disk but is invisible to the report. If only the second step is done, it will not be included in the final `adata`.

---

## Resource labels

Use nf-core-style labels defined in `conf/base.config`:

| Label | Typical use |
|---|---|
| `process_low` | fast steps, small memory |
| `process_medium` | most Python/R tools |
| `process_high` | multi-core steps |
| `process_high_memory` | large atlas loading, BayesTME |

Do not hard-code CPU/memory in a process; always use a label.

---

## Commit style

- Use the imperative mood in the subject line: _"Add CoGAPS cross-sample integration"_ not _"Added…"_.
- Keep the subject line under 72 characters.
- Reference relevant issues or PRs in the body: `Closes #123`.
- Prefer small, atomic commits over large "fix everything" commits.

---

Questions? Open a [GitHub Discussion](https://github.com/break-through-cancer/staple/discussions) or reach out via the issue tracker.
