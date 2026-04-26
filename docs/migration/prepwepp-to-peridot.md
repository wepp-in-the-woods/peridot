# prepwepp/TOPAZ to Peridot Migration Guide
> Migration framing for moving from implicit legacy watershed abstraction to Peridot's explicit graph abstraction.

## Audience

Use this guide when replacing a prepwepp, TOPAZ/TOP2WEPP, or older WEPPcloud abstraction path with Peridot outputs. It is written for maintainers who need to understand which behavior is expected to match, which behavior intentionally differs, and where downstream consumers should point for canonical contracts.

## Mental Model Shift

Legacy abstraction workflows commonly treat the watershed as a collection of raster classes and generated text files. Relationships between channels, hillslopes, and flowpaths are implied by TOPAZ IDs, D8 flow directions, and naming conventions.

Peridot makes the relationship layer explicit. It reads TOPAZ or WBT delineation artifacts, constructs channel and hillslope flowpath collections, and writes a documented set of WEPP slope files plus tabular outputs. The migration is therefore not just "same workflow in Rust"; it is a move from implicit raster convention to explicit graph abstraction and output contracts.

## Replacement Scope

Peridot replaces:

- channel and hillslope abstraction after delineation,
- generation of WEPP channel and hillslope `.slp` files,
- optional full flowpath table and `.slps` bundle generation,
- watershed metadata table generation,
- run-local output manifest generation,
- sub-field slope/profile generation for AgFields when `sub_fields_abstraction` is used.

Peridot does not replace:

- DEM acquisition or conditioning,
- WhiteboxTools or TOPAZ terrain delineation,
- stream-network pruning policy,
- WEPP model execution,
- WEPPpy NoDb state, RQ orchestration, reports, or query-engine activation,
- scientific validation of a new delineation rule.

## Legacy-to-Peridot Mapping

| Legacy expectation | Peridot behavior | Migration note |
| --- | --- | --- |
| TOPAZ `.ARC` inputs drive abstraction. | `abstract_watershed` still consumes TOPAZ `.ARC` inputs. | Use for compatibility paths that still produce `dem/topaz/*`. |
| WBT-derived TOPAZ-compatible IDs drive abstraction. | `wbt_abstract_watershed` consumes `dem/wbt/*` rasters and `netw.tsv`. | Preferred for current WEPPpy WBT workflows. |
| Slope files are the primary output. | Slope files remain required, but Parquet tables and `watershed/README.md` are also first-class outputs. | Downstream consumers should validate both files and tables. |
| CSV watershed summaries may be expected. | Current Peridot watershed CLIs write Parquet tables directly, not watershed CSV summaries. | Update readers to use `channels.parquet`, `hillslopes.parquet`, and optional `flowpaths.parquet`; WEPPpy compatibility code may handle historical CSV-only runs separately. |
| Flowpaths are always materialized. | Flowpaths are optional and are omitted by `--skip-flowpaths` and by WBT representative-flowpath mode. | Consumers must tolerate absent `flowpaths.parquet` when the manifest records flowpaths as skipped. |
| Length/width derivation is implicit. | Hillslope rows include length provenance fields. | Use `length_estimate_mode`, `length_area_over_channel`, and `length_edge_median` to audit side-hillslope behavior. |
| Field boundaries are external to watershed abstraction. | `sub_fields_abstraction` can intersect field boundaries with WBT hillslopes and emit sub-field profiles. | This is an extension of the explicit graph model, not a legacy parity requirement. |

## Migration Checklist

1. Confirm the delineation stage produces one supported input layout: `dem/topaz/*` for `abstract_watershed` or `dem/wbt/*` for `wbt_abstract_watershed`.
2. Choose the abstraction mode: full flowpaths, `--skip-flowpaths`, or WBT `--representative-flowpath`.
3. Run the Peridot binary with explicit `--ncpu`, `--max-points`, and width/clipping flags.
4. Validate required outputs against [the output contract](../contracts/watershed-output-contract.md).
5. Update downstream readers to use Parquet-first watershed tables.
6. If migrating historical runs, keep compatibility code separate from the current Peridot output contract.
7. Record benchmark claims using [the benchmark discipline](../benchmarks.md) before describing performance improvements as measured results.

## Intentional Differences

Flowpath export can be intentionally absent. This is not a failed abstraction when `--skip-flowpaths` or `--representative-flowpath` is recorded in `watershed/README.md`.

Representative-flowpath mode is an abstraction-mode choice. It uses one deterministic source-cell flowpath per hillslope and records representative length provenance. It should not be described as bit-for-bit parity with per-pixel flowpath expansion.

Peridot's current watershed contract is Parquet-first. Documentation or code that still requires `watershed/channels.csv`, `watershed/hillslopes.csv`, or `watershed/flowpaths.csv` is depending on historical or compatibility behavior, not the current Peridot CLI output.

Sub-field output currently uses CSV metadata and has a known duplicated `topaz_id` header in `field_flowpaths.csv`. Treat that as current behavior until a schema cleanup package is approved.

## WEPPpy Integration Notes

WEPPpy commonly calls Peridot as part of watershed abstraction and then performs additional post-processing. That post-processing may add canonical WEPPpy columns, refresh catalogs, or normalize historical outputs. Do not infer those WEPPpy additions as direct Peridot CLI outputs.

When updating WEPPpy docs, point readers to these canonical Peridot docs for direct Peridot behavior:

- [Output contract](../contracts/watershed-output-contract.md)
- [Operations runbook](../operations.md)
- [Benchmarks and evidence discipline](../benchmarks.md)

## Evidence Notes

- `confirmed`: Peridot has separate TOPAZ and WBT CLI entrypoints and shared watershed abstraction routines.
- `confirmed`: Current watershed CLIs write Parquet outputs and generated manifests.
- `confirmed`: WBT representative-flowpath mode disables full flowpath export.
- `inference`: The explicit graph abstraction improves migration legibility because topology and provenance are exposed as data structures and tables.
- `hypothesis`: Workload-specific speedups require benchmark artifacts before being reused as measured migration claims.
