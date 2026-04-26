# peridot
> Programmable Environmental Rust Interface for Drainage & Operational Topography.

Peridot is the watershed abstraction layer used by WEPPpy to turn hydrologic delineation outputs into WEPP-ready slope profiles and auditable tabular summaries. Its role is not only to make an older TOPAZ/TOP2WEPP path faster. Peridot changes the abstraction surface: watershed topology becomes an explicit graph of channels, hillslopes, flowpaths, and field intersections instead of remaining implicit in raster IDs and flow-direction side effects.

## Why This Matters

WEPP workflows depend on the abstraction step between terrain delineation and model execution. If that step is opaque, downstream users only see generated slope files and have little evidence about how channels, hillslopes, lengths, widths, and optional flowpaths were derived. Peridot makes that boundary inspectable by producing both WEPP slope files and machine-readable watershed tables with a generated manifest.

The practical result is a watershed abstraction layer that is easier to validate, migrate, and operate at WEPPpy scale. Speed and reduced output volume are downstream effects; the primary shift is from implicit raster segmentation to explicit topology plus documented outputs.

## Category Shift, Not Just Modernization

Clean claim statement:

> Legacy TOPAZ/TOP2WEPP workflows treat watershed structure as an implicit raster-and-file convention. Peridot treats watershed structure as an explicit graph abstraction, then emits WEPP slope files and tabular contracts from that graph.

| Legacy mental model | Peridot mental model |
| --- | --- |
| Raster classes, D8 flow directions, and sidecar text files imply the watershed relationships. | Channels, hillslopes, sub-flowpaths, and field intersections are represented as explicit graph elements in Rust data structures. |
| The abstraction is tightly coupled to a TOPAZ/TOP2WEPP-style file workflow. | The abstraction layer can consume TOPAZ or WhiteboxTools-derived inputs while keeping the downstream contract stable. |
| Validation often means checking whether expected slope files exist. | Validation includes slope files, Parquet tables, generated run flags, schema summaries, and diagnostics. |
| Scaling behavior is dominated by per-pixel flowpath expansion unless outputs are manually pruned. | Full flowpath export can be disabled, and WBT representative-flowpath mode can intentionally avoid per-pixel flowpath output. |

## Replacement Scope and Boundaries

Peridot replaces the watershed abstraction stage: it consumes existing delineation artifacts, walks the topology, builds representative hillslope/channel profiles, and writes WEPP-facing outputs.

Peridot does not replace:

- DEM preparation, stream pruning, or basin delineation. WEPPpy commonly uses `weppcloud-wbt` for those terrain-processing steps.
- WEPP physics or WEPP input/output execution.
- WEPPpy post-processing, reports, query activation, or run-state controllers.
- Every historical TOPAZ/TOP2WEPP behavior bit-for-bit. Some differences, such as representative-flowpath mode and side-hillslope length provenance fields, are intentional and documented.

The legacy `abstract_watershed` binary remains for TOPAZ `.ARC` inputs. The WBT path, `wbt_abstract_watershed`, is the current production-oriented path for WhiteboxTools-derived WEPPcloud runs. `sub_fields_abstraction` extends the abstraction model to agricultural field/subcatchment intersections for AgFields workflows.

## Canonical Documentation

Use these documents as the source of truth before relying on README summaries:

- [Watershed output contract](docs/contracts/watershed-output-contract.md) - required outputs, optional outputs, schema fields, flags, and current error-boundary notes.
- [Benchmarks and evidence discipline](docs/benchmarks.md) - how to measure scalability, topology correctness/flexibility, and parallelization potential without overstating unmeasured speedups.
- [prepwepp/TOPAZ to Peridot migration guide](docs/migration/prepwepp-to-peridot.md) - replacement scope, parity expectations, intentional differences, and migration checklist.
- [Operations runbook](docs/operations.md) - commands, preflight checks, failure signatures, and output validation.
- [Sub-fields design spec](sub_fields_abstraction.spec.md) - implementation context for AgFields sub-field abstraction.
- [Development notes](dev-notes.md) - implementation notes that are useful after reading the canonical contract docs.

## Communication Kit

Claim statement:

- Legacy model: raster IDs and flow-direction files implicitly encode watershed structure.
- Peridot model: an explicit graph abstraction encodes topology first, then writes WEPP slope files and tabular outputs.

Figure specification:

- Show the same watershed twice. On the left, label the legacy view as raster discretization: colored raster zones, D8 arrows, channel pixels, and file names such as `SUBWTA.ARC` or `subwta.tif`. On the right, label the Peridot view as graph abstraction: channel nodes/edges, hillslope nodes linked to channels, optional sub-flowpath bundles, and output tables (`channels.parquet`, `hillslopes.parquet`, optional `flowpaths.parquet`). The caption should state that both views can start from the same delineated terrain, but Peridot exposes the topology as the primary abstraction surface.

Metric definitions:

- Scalability: wall time, peak memory, and output volume as functions of raster cells, channel count, hillslope count, and optional flowpath row count.
- Topology flexibility/correctness: whether explicit channel/hillslope relationships, IDs, lengths, widths, and provenance fields remain internally consistent across TOPAZ, WBT, representative-flowpath, and sub-field modes.
- Parallelization potential: how much work can be decomposed by independent graph elements or output writers, measured by thread-count scaling and by the amount of state that can be partitioned without changing results.

## Binaries

Peridot ships these command-line entrypoints:

- `abstract_watershed` consumes TOPAZ outputs stored as `.ARC` rasters.
- `wbt_abstract_watershed` consumes WhiteboxTools outputs generated by `weppcloud-wbt` conventions.
- `sub_fields_abstraction` consumes WhiteboxTools outputs plus a rasterized field-boundary map to split hillslopes into sub-fields for AgFields runs.
- `trace_downslope_flowpath` traces roads/downstream flowpaths for integrations that need a focused point-source trace.

`abstract_watershed` and `wbt_abstract_watershed` share the core flowpath routines in `peridot::watershed_abstraction`; WBT-specific entrypoints live under `peridot::wbt`.

## Expected Data Layout

The working directory passed to any binary is mutated in place. The following structure is assumed to exist before execution:

```text
<run>/
  dem/
    topaz/        # for abstract_watershed
      SUBWTA.ARC
      RELIEF.ARC
      FLOVEC.ARC
      FVSLOP.ARC
      TASPEC.ARC
      NETW.TAB
    wbt/          # for WBT-based tools
      subwta.tif  # or subwta.vrt
      relief.tif  # or relief.vrt
      flovec.tif  # or flovec.vrt
      fvslop.tif  # or fvslop.vrt
      taspec.tif  # or taspec.vrt
      netw.tsv    # for wbt_abstract_watershed
      discha.tif  # or discha.vrt; required for representative-flowpath mode
  ag_fields/      # for sub_fields_abstraction
    field_boundaries.tif
```

WBT-based tools remap WhiteboxTools D8 flow directions to match TOPAZ conventions before processing, so downstream abstractions can preserve TOPAZ-style hillslope/channel IDs.

## Output Contract Summary

Running `abstract_watershed` or `wbt_abstract_watershed` creates or refreshes `<run>/watershed/`. Required current watershed outputs are:

- `slope_files/channels.slp`
- `slope_files/hillslopes/hill_<topaz_id>.slp`
- `channels.parquet`
- `hillslopes.parquet`
- `channels.geojson`
- `network.txt`
- `README.md`, generated as the per-run output manifest

When full flowpath export is enabled, Peridot also writes:

- `flowpaths.parquet`
- `slope_files/flowpaths/fps_<topaz_id>.slps`

`--skip-flowpaths` omits those flowpath outputs. `--representative-flowpath` is WBT-only and forces flowpath export off while using one deterministic representative flowpath per hillslope.

Running `sub_fields_abstraction` creates or refreshes `<run>/ag_fields/sub_fields/` by default:

- `sub_field_id_map.tif`
- `slope_files/field_<field_id>_<topaz_id>.slp`
- `slope_files/flowpaths/field_<field_id>_<topaz_id>.slps`
- `fields.csv`
- `field_flowpaths.csv`

See the [watershed output contract](docs/contracts/watershed-output-contract.md) for schemas, conditional behavior, and current error-boundary notes.

## Representative Flowpath Mode (WBT Only)

`wbt_abstract_watershed --representative-flowpath` replaces the per-pixel hillslope flowpath expansion with one deterministic flowpath per hillslope chosen from source-cell candidates using `dem/wbt/discha.tif` and tie-breakers. This intentionally changes the abstraction strategy for large batch workflows:

- Full flowpath table and `.slps` export is disabled.
- Hillslope profile generation still emits `hillslopes.parquet` and slope files.
- Side hillslopes still apply the side-length selection rules and record provenance fields.
- Top hillslopes use `length_estimate_mode = top_representative_flowpath`.

Treat this as a documented abstraction mode, not a hidden optimization of the legacy per-pixel model.

## Building

```bash
cargo build --release
```

GDAL bindings are enabled through the `gdal` dependency's `bindgen` feature in `Cargo.toml`.

To build inside the `wepppy310` Conda environment, activate it first:

```bash
source set_wepppy310_env.sh
cargo build --release
```

## Running the Binaries

WEPP supports up to 100 points per slope profile; WEPPpy usually keeps `--max-points` at or below 30.

```bash
./target/release/abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --clip-hillslopes --clip-hillslope-length 300

./target/release/wbt_abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --bieger2015-widths

./target/release/wbt_abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --representative-flowpath

./target/release/sub_fields_abstraction /path/to/run \
  --ncpu 8 --max-points 30 \
  --sub-field-min-area-threshold-m2 500 \
  --field-raster ag_fields/field_boundaries.tif \
  --output-dir ag_fields/sub_fields
```

Key options:

- `--ncpu` sets the Rayon thread pool size.
- `--max-points` caps generated WEPP slope profile point counts.
- `--clip-hillslopes` and `--clip-hillslope-length` clip hillslope slope profiles to a maximum physical length.
- `--bieger2015-widths` uses the Bieger 2015 regressions to infer channel widths from drainage area.
- `--skip-flowpaths` omits full flowpath table and `.slps` exports for watershed abstraction.
- `--representative-flowpath` uses a single WBT representative path per hillslope and forces `--skip-flowpaths` behavior.
- `--sub-field-min-area-threshold-m2` drops field/subcatchment intersections smaller than the specified area.
- `--field-raster` and `--output-dir` override the default AgFields raster location and output directory for the sub-field tool.

## Integrating With WEPPpy

Copy the compiled binaries used by WEPPpy into the location expected by the Python package:

```bash
cp /workdir/peridot/target/release/abstract_watershed \
   /workdir/wepppy/wepppy/topo/peridot/bin/
cp /workdir/peridot/target/release/wbt_abstract_watershed \
   /workdir/wepppy/wepppy/topo/peridot/bin/
```

If you need the sub-fields workflow, copy `sub_fields_abstraction` into the same `bin/` directory.

For the Conda-based deployment:

```bash
cp /workdir/peridot/target/release/abstract_watershed \
   /workdir/wepppy/wepppy/topo/peridot/bin/abstract_watershed.conda310.ub2404
cp /workdir/peridot/target/release/wbt_abstract_watershed \
   /workdir/wepppy/wepppy/topo/peridot/bin/wbt_abstract_watershed.conda310.ub2404
```

Tests live under `tests/` and can be executed with `cargo test`.

## Source Tree

```text
peridot/src
.
├── bin
│   ├── abstract_watershed.rs
│   ├── sub_fields_abstraction.rs
│   ├── trace_downslope_flowpath.rs
│   └── wbt_abstract_watershed.rs
├── lib.rs
├── logging.rs
├── rasters
├── roads_trace
├── support
├── topaz
├── watershed_abstraction
└── wbt
```
