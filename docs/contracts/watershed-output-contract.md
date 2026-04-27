# Peridot Watershed Output Contract
> Canonical contract for Peridot watershed abstraction outputs, conditional output modes, and schema fields.

## Normative Status

This document defines the current Peridot output contract for repository documentation and WEPPpy integration. If this document and the implementation diverge, either the implementation must be corrected in a runtime package or this contract must be updated in the same change set.

This contract is documentation-only. It does not change runtime behavior.

## Scope

Covered commands:

- `abstract_watershed`: TOPAZ `.ARC` input mode.
- `wbt_abstract_watershed`: WhiteboxTools-derived WBT input mode.
- `sub_fields_abstraction`: AgFields sub-field output mode built from WBT watershed inputs plus a field-boundary raster.

Out of scope:

- DEM preprocessing, watershed delineation, and stream pruning.
- WEPP model execution and WEPP output reports.
- WEPPpy post-processing that may add derived columns such as `wepp_id`, `chn_enum`, or query-engine catalog metadata after Peridot completes.
- Historical TOPAZ/TOP2WEPP output formats not emitted by current Peridot binaries.

## Terms

- Run directory: the working directory passed as the first positional argument to a Peridot binary.
- `topaz_id`: TOPAZ-style watershed element ID preserved across TOPAZ and WBT input modes.
- Channel ID: a `topaz_id` ending in `4`.
- Hillslope ID: a positive `topaz_id` not ending in `4`; common suffixes are `1`, `2`, and `3`.
- Flowpath export: optional per-pixel/sub-flowpath output under `watershed/flowpaths.parquet` and `watershed/slope_files/flowpaths/`.
- Representative-flowpath mode: WBT-only mode that chooses one deterministic source-cell flowpath per hillslope and disables full flowpath export.

## Input Layout

`abstract_watershed` expects TOPAZ files under `dem/topaz/`:

- `SUBWTA.ARC`
- `RELIEF.ARC`
- `FLOVEC.ARC`
- `FVSLOP.ARC`
- `TASPEC.ARC`
- `NETW.TAB`

`wbt_abstract_watershed` expects WBT files under `dem/wbt/`:

- `subwta.tif` or `subwta.vrt`
- `relief.tif` or `relief.vrt`
- `flovec.tif` or `flovec.vrt`
- `fvslop.tif` or `fvslop.vrt`
- `taspec.tif` or `taspec.vrt`
- `netw.tsv`
- `discha.tif` or `discha.vrt` when `--representative-flowpath` is enabled

`sub_fields_abstraction` expects the WBT raster stack used to build sub-field slope profiles, but it does not read `netw.tsv` or `discha`:

- `subwta.tif` or `subwta.vrt`
- `relief.tif` or `relief.vrt`
- `flovec.tif` or `flovec.vrt`
- `fvslop.tif` or `fvslop.vrt`
- `taspec.tif` or `taspec.vrt`
- field raster from `--field-raster`, defaulting to `ag_fields/field_boundaries.tif` relative to the run directory

## Mutation Contract

`abstract_watershed` and `wbt_abstract_watershed` delete any existing `<run>/watershed/` directory before regenerating outputs.

`sub_fields_abstraction` deletes its configured output directory before regenerating outputs. The default output directory is `<run>/ag_fields/sub_fields/`.

Callers must not store unrelated durable files inside those output directories unless they are prepared for them to be removed on the next abstraction run.

## Watershed Outputs

Current `abstract_watershed` and `wbt_abstract_watershed` runs write Parquet outputs directly. They do not currently write watershed CSV compatibility files from the CLI path.

Required outputs when the command succeeds and full flowpath export is disabled:

| Path | Status | Notes |
| --- | --- | --- |
| `watershed/slope_files/channels.slp` | required | WEPP channel slope profile bundle. |
| `watershed/slope_files/hillslopes/hill_<topaz_id>.slp` | required | One WEPP hillslope slope file per hillslope. |
| `watershed/channels.parquet` | required | Channel metadata table. |
| `watershed/hillslopes.parquet` | required | Hillslope metadata table. |
| `watershed/channels.geojson` | required | Lightweight channel diagnostic geometry. |
| `watershed/network.txt` | required | Connectivity diagnostic emitted from the parsed network table. |
| `watershed/README.md` | required | Generated manifest containing flags, file listing, and schema summaries. |

Additional outputs when full flowpath export is enabled:

| Path | Status | Notes |
| --- | --- | --- |
| `watershed/flowpaths.parquet` | conditional | Omitted when `--skip-flowpaths` is set or representative-flowpath mode is enabled. |
| `watershed/slope_files/flowpaths/fps_<topaz_id>.slps` | conditional | Omitted under the same conditions as `flowpaths.parquet`. |

The `watershed/slope_files/flowpaths/` directory may exist even when flowpath export is disabled; consumers should check for expected files, not only the directory.

## Watershed Schema: `channels.parquet`

| Column | Type | Meaning |
| --- | --- | --- |
| `topaz_id` | `Int32` | Channel element ID. |
| `slope_scalar` | `Float64` | Channel slope summary. |
| `length` | `Float64` | Channel length. |
| `width` | `Float64` | Channel width. |
| `direction` | `Float64` | Flowpath direction summary. |
| `order` | `Int32` | Channel order from the parsed network. |
| `aspect` | `Float64` | Aspect summary. |
| `area` | `Float64` | Channel area in square meters. |
| `elevation` | `Float64` | Elevation summary. |
| `centroid_px` | `Int32` | Centroid pixel x coordinate. |
| `centroid_py` | `Int32` | Centroid pixel y coordinate. |
| `centroid_lon` | `Float64` | Centroid longitude. |
| `centroid_lat` | `Float64` | Centroid latitude. |

## Watershed Schema: `hillslopes.parquet`

| Column | Type | Meaning |
| --- | --- | --- |
| `topaz_id` | `Int32` | Hillslope element ID. |
| `slope_scalar` | `Float64` | Hillslope slope summary. |
| `length` | `Float64` | Representative hillslope length. |
| `width` | `Float64` | Representative hillslope width. |
| `direction` | `Float64` | Flowpath direction summary. |
| `aspect` | `Float64` | Aspect summary. |
| `length_estimate_mode` | `Utf8` | Provenance label for the selected length rule. |
| `length_area_over_channel` | `Float64` | Side-hillslope candidate length from area divided by channel length; `NaN` when not applicable. |
| `length_edge_median` | `Float64` | Edge/source or representative length candidate; `NaN` when unavailable. |
| `area` | `Int32` | Hillslope area in square meters, currently stored as an integer cast. |
| `elevation` | `Float64` | Elevation summary. |
| `centroid_px` | `Int32` | Centroid pixel x coordinate. |
| `centroid_py` | `Int32` | Centroid pixel y coordinate. |
| `centroid_lon` | `Float64` | Centroid longitude. |
| `centroid_lat` | `Float64` | Centroid latitude. |

Current `length_estimate_mode` values include:

- `top_edge_median`
- `top_representative_flowpath`
- `side_edge_median_capped`
- `side_area_over_channel`
- `side_area_over_channel_no_edge`

## Watershed Schema: `flowpaths.parquet`

`flowpaths.parquet` is conditional. It is emitted only when full flowpath export is enabled.

| Column | Type | Meaning |
| --- | --- | --- |
| `topaz_id` | `Int32` | Parent hillslope ID for the sub-flowpath. |
| `fp_id` | `Int32` | Flowpath ID within the parent hillslope. |
| `slope_scalar` | `Float64` | Flowpath slope summary. |
| `length` | `Float64` | Flowpath length. |
| `width` | `Float64` | Flowpath width. |
| `direction` | `Float64` | Flow direction summary. |
| `aspect` | `Float64` | Aspect summary. |
| `area` | `Float64` | Flowpath area in square meters. |
| `elevation` | `Float64` | Elevation summary. |
| `order` | `Int32` | Flowpath order metadata. |
| `centroid_px` | `Int32` | Centroid pixel x coordinate. |
| `centroid_py` | `Int32` | Centroid pixel y coordinate. |
| `centroid_lon` | `Float64` | Centroid longitude. |
| `centroid_lat` | `Float64` | Centroid latitude. |

## Flag Behavior

| Flag | Commands | Contract |
| --- | --- | --- |
| `--ncpu` | all commands | Configures Rayon thread-pool size for the process. |
| `--max-points` | all commands | Caps generated WEPP slope profile point counts; it does not cap table row counts. |
| `--clip-hillslopes` | all commands | Clips generated hillslope slope profiles when enabled. |
| `--clip-hillslope-length` | all commands | Sets the clipping length used when clipping is enabled. |
| `--bieger2015-widths` | watershed commands | Uses Bieger 2015 channel-width regressions. |
| `--skip-flowpaths` | watershed commands | Omits `flowpaths.parquet` and `slope_files/flowpaths/*.slps`. |
| `--representative-flowpath` | `wbt_abstract_watershed` only | Forces flowpath export off and builds hillslope summaries from one deterministic source-cell flowpath per hillslope. |
| `--sub-field-min-area-threshold-m2` | `sub_fields_abstraction` only | Drops field/subcatchment intersections smaller than the threshold. |
| `--field-raster` | `sub_fields_abstraction` only | Selects the field-boundary raster relative to the run directory. |
| `--output-dir` | `sub_fields_abstraction` only | Selects the sub-field output directory relative to the run directory. |

## Generated Manifest Contract

`watershed/README.md` is generated by Peridot for each watershed abstraction run. It records:

- command name and runtime flags,
- file manifest with formats, sizes, and row counts where known,
- tabular schema summaries for generated Parquet tables,
- conditional notes for flowpath export, representative-flowpath mode, clipping, and channel-width mode.

Consumers should treat the generated manifest as run-local evidence, not as a replacement for this repository-level contract.

## Sub-Field Outputs

`sub_fields_abstraction` writes CSV and raster/slope artifacts under its configured output directory. The default is `ag_fields/sub_fields/`.

| Path | Status | Notes |
| --- | --- | --- |
| `sub_field_id_map.tif` | required | Raster map of retained field/subcatchment intersections. |
| `slope_files/field_<field_id>_<topaz_id>.slp` | required | Representative slope file per retained sub-field. |
| `slope_files/flowpaths/field_<field_id>_<topaz_id>.slps` | required | Sub-field flowpath bundle. |
| `fields.csv` | required | Sub-field metadata table. |
| `field_flowpaths.csv` | required | Sub-field flowpath metadata table. |

`fields.csv` columns:

| Column | Meaning |
| --- | --- |
| `field_id` | Source field ID from the field-boundary raster. |
| `topaz_id` | Parent hillslope ID. |
| `sub_field_id` | Peridot-assigned sub-field ID. |
| `slope_scalar` | Slope summary. |
| `length` | Representative length. |
| `width` | Representative width. |
| `direction` | Flow direction summary. |
| `aspect` | Aspect summary. |
| `area` | Area in square meters. |
| `elevation` | Elevation summary. |
| `centroid_px` | Centroid pixel x coordinate. |
| `centroid_py` | Centroid pixel y coordinate. |
| `centroid_lon` | Centroid longitude. |
| `centroid_lat` | Centroid latitude. |

`field_flowpaths.csv` columns:

| Column | Meaning |
| --- | --- |
| `field_id` | Source field ID from the field-boundary raster. |
| `topaz_id` | Parent hillslope ID. |
| `sub_field_id` | Peridot-assigned sub-field ID. |
| `flowpath_topaz_id` | Topaz ID stored on the sub-field flowpath record. |
| `fp_id` | Flowpath ID within the sub-field flowpath bundle. |
| `slope_scalar` | Slope summary. |
| `length` | Flowpath length. |
| `width` | Flowpath width. |
| `direction` | Flow direction summary. |
| `aspect` | Aspect summary. |
| `area` | Area in square meters. |
| `elevation` | Elevation summary. |
| `order` | Flowpath order. |
| `centroid_px` | Centroid pixel x coordinate. |
| `centroid_py` | Centroid pixel y coordinate. |
| `centroid_lon` | Centroid longitude. |
| `centroid_lat` | Centroid latitude. |

## Current Error Boundary

Many missing-input cases fail by panic because raster and network reads use explicit `unwrap()` calls. Those failures normally produce a non-zero process exit.

The `abstract_watershed`, `wbt_abstract_watershed`, and `sub_fields_abstraction` CLI entrypoints return the underlying abstraction `io::Result<()>`. Propagated write-stage errors therefore make the process exit non-zero.

Operational callers should still validate required outputs and generated manifests after each run because a zero exit status only confirms that the command returned success, not that downstream post-processing or deployment-specific expectations were satisfied.
