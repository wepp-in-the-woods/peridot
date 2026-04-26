# Peridot Benchmarks and Evidence Discipline
> Methodology and claim boundaries for Peridot performance, scalability, topology, and parallelization statements.

## Purpose

Peridot documentation may describe architectural advantages, but measured performance claims must be reproducible. This document separates confirmed implementation facts, architecture inferences, and benchmark hypotheses so Peridot is presented as an explicit graph abstraction rather than as an unsupported speedup story.

## Claim Discipline

Use these labels when citing evidence:

- `confirmed`: directly evidenced by Peridot source, tests, generated artifacts, or a recorded benchmark run.
- `inference`: reasoned from source structure or observed behavior, but not independently measured as a benchmark result.
- `hypothesis`: plausible and directionally useful, but requiring a dedicated benchmark or publication-quality validation before being presented as a result.

Do not convert an `inference` or `hypothesis` into a measured claim without adding benchmark artifacts that include dataset, command, hardware, version, and output.

## Communication Kit

Clean claim statement:

> Peridot replaces the implicit raster legacy model with an explicit watershed graph abstraction, then emits WEPP slope files and auditable tabular contracts from that graph.

Figure specification:

- Use one watershed and draw two panels.
- Left panel: legacy raster discretization with channel pixels, hillslope raster zones, D8 arrows, and TOPAZ/TOP2WEPP-style file labels.
- Right panel: Peridot graph abstraction with channel edges, hillslope nodes, optional sub-flowpath bundles, and generated output tables.
- Caption: "The terrain can be the same; the abstraction surface changes from implicit raster convention to explicit graph topology."

Core metric definitions:

| Metric | Definition | Minimum evidence |
| --- | --- | --- |
| Scalability | Wall time, peak RSS, and output bytes as functions of raster cells, channel count, hillslope count, and optional flowpath row count. | Dataset dimensions, command, hardware, binary version, `/usr/bin/time -v` output, and generated manifest. |
| Topology flexibility/correctness | Preservation of explicit channel/hillslope relationships and geometry/provenance fields across TOPAZ, WBT, representative-flowpath, and sub-field modes. | Invariant checks for ID classes, network reachability, row counts, length/width sanity, and schema presence. |
| Parallelization potential | Amount of abstraction and output work that can be partitioned by graph element or writer task without changing outputs. | Thread-scaling runs with fixed inputs plus checksum/schema comparisons across `--ncpu` values. |

## Current Evidence Register

| Claim | Status | Evidence |
| --- | --- | --- |
| Peridot watershed CLIs emit Parquet tables and a generated manifest. | `confirmed` | `src/watershed_abstraction/watershed_abstraction.rs`, `src/wbt/wbt_watershed_abstraction.rs`, `src/watershed_abstraction/watershed_manifest.rs`, and `tests/watershed_parquet_manifest.rs`. |
| `--representative-flowpath` forces full flowpath export off. | `confirmed` | `src/bin/wbt_abstract_watershed.rs` and `src/wbt/wbt_watershed_abstraction.rs`. |
| Peridot exposes a graph-oriented abstraction boundary rather than only a TOPAZ/TOP2WEPP file conversion. | `inference` | Flowpath collections, channel network parsing, explicit output schemas, and WBT/TOPAZ input modes share downstream abstraction routines. |
| Representative-flowpath mode can reduce output volume for large runs by omitting per-pixel flowpath outputs. | `confirmed` for output omission; `hypothesis` for run-specific speedup until measured on the target dataset. | Flag contract and output writers confirm omission. Wall-time gains require benchmark artifacts. |
| Historical statements such as "3x to 10x faster" or "10x to 100x" are generally true for all workloads. | `hypothesis` unless tied to a specific recorded benchmark. | Existing narrative claims require dataset/hardware/version evidence before reuse as measured results. |

## Reproducible Benchmark Template

Record one benchmark file per dataset and mode. Include all fields below.

```text
Dataset label:
Run directory or fixture:
Data date/source:
Raster dimensions:
Channel count:
Hillslope count:
Flowpath rows if exported:
Peridot commit:
Binary path and --version output:
Host CPU model:
Host RAM:
Storage type:
GDAL version:
Command:
/usr/bin/time -v output:
Generated watershed/README.md excerpt:
Notes on cold/warm cache:
```

Recommended command shape:

```bash
RUST_LOG=info /usr/bin/time -v \
  ./target/release/wbt_abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --bieger2015-widths

RUST_LOG=info /usr/bin/time -v \
  ./target/release/wbt_abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --representative-flowpath
```

For topology and schema checks, preserve the generated manifest and run a table inspection command appropriate to the environment, for example DuckDB or a small Python/PyArrow script. The artifact must include the exact command used.

## Benchmark Matrix

A useful benchmark set should cover:

| Mode | Required comparison |
| --- | --- |
| TOPAZ input, full flowpaths | Legacy compatibility and output contract baseline. |
| WBT input, full flowpaths | Current production-oriented full-output mode. |
| WBT input, `--skip-flowpaths` | Output-volume and memory-pressure reduction from omitting full flowpaths. |
| WBT input, `--representative-flowpath` | Intentional abstraction-mode shift with representative hillslope paths. |
| Sub-fields | Field/subcatchment intersection output volume and metadata correctness. |

For each mode, record `--ncpu 1` and at least one production-like `--ncpu` value before making parallel speedup claims.

## Reporting Rules

- Report speedups only as `baseline_time / candidate_time` for the same dataset, host, storage, and binary version unless the difference is explicitly labeled.
- Report memory as peak RSS from the same measurement tool across runs.
- Report output volume as total bytes under the generated output directory plus row counts from the manifest or table reader.
- State whether flowpaths were exported, skipped, or replaced by representative-flowpath mode.
- State whether the result is `confirmed`, `inference`, or `hypothesis`.
