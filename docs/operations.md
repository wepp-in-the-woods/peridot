# Peridot Operations Runbook
> Operational checks for running Peridot watershed and sub-field abstraction safely.

## Scope

This runbook covers local or WEPPpy-hosted execution of Peridot binaries. It is documentation-only and does not change runtime behavior.

## Preflight Checklist

Before running watershed abstraction:

- Confirm the run directory is writable.
- Confirm the expected input layout exists under `dem/topaz/` or `dem/wbt/`.
- Confirm `dem/wbt/discha.tif` or `dem/wbt/discha.vrt` exists before using `--representative-flowpath`.
- Confirm the caller is prepared for Peridot to delete and recreate `watershed/`.
- Decide whether full flowpath export is required; use `--skip-flowpaths` or `--representative-flowpath` for workflows that do not consume per-pixel flowpaths.

Before running sub-field abstraction:

- Confirm WBT watershed inputs exist.
- Confirm the field-boundary raster exists, defaulting to `ag_fields/field_boundaries.tif`.
- Confirm the caller is prepared for Peridot to delete and recreate the configured output directory, defaulting to `ag_fields/sub_fields/`.

## Standard Commands

TOPAZ input mode:

```bash
RUST_LOG=info ./target/release/abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --clip-hillslopes --clip-hillslope-length 300
```

WBT input mode:

```bash
RUST_LOG=info ./target/release/wbt_abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --bieger2015-widths
```

WBT representative-flowpath mode:

```bash
RUST_LOG=info ./target/release/wbt_abstract_watershed /path/to/run \
  --ncpu 8 --max-points 30 --representative-flowpath
```

Sub-field mode:

```bash
RUST_LOG=info ./target/release/sub_fields_abstraction /path/to/run \
  --ncpu 8 --max-points 30 \
  --sub-field-min-area-threshold-m2 500 \
  --field-raster ag_fields/field_boundaries.tif \
  --output-dir ag_fields/sub_fields
```

## Post-Run Validation

Treat non-zero process status as a failed abstraction. Still validate required outputs after a zero exit status so deployment-specific expectations and downstream post-processing assumptions are checked explicitly.

Minimum watershed validation:

```bash
run=/path/to/run

test -s "$run/watershed/README.md"
test -s "$run/watershed/channels.parquet"
test -s "$run/watershed/hillslopes.parquet"
test -s "$run/watershed/slope_files/channels.slp"
test -s "$run/watershed/channels.geojson"
test -s "$run/watershed/network.txt"
find "$run/watershed/slope_files/hillslopes" -type f -name 'hill_*.slp' -print -quit | grep -q .
```

If full flowpaths were requested, also validate:

```bash
test -s "$run/watershed/flowpaths.parquet"
find "$run/watershed/slope_files/flowpaths" -type f -name 'fps_*.slps' -print -quit | grep -q .
```

If `--skip-flowpaths` or `--representative-flowpath` was used, `flowpaths.parquet` and `fps_*.slps` are expected to be absent.

Minimum sub-field validation:

```bash
run=/path/to/run
out="$run/ag_fields/sub_fields"

test -s "$out/sub_field_id_map.tif"
test -s "$out/fields.csv"
test -s "$out/field_flowpaths.csv"
find "$out/slope_files" -maxdepth 1 -type f -name 'field_*_*.slp' -print -quit | grep -q .
find "$out/slope_files/flowpaths" -type f -name 'field_*_*.slps' -print -quit | grep -q .
```

## Failure Signatures

| Symptom | Likely cause | Operator action |
| --- | --- | --- |
| Panic mentioning missing raster or `unwrap()` | Required input is missing or unreadable. | Check input layout and file permissions. |
| Panic or failure when representative mode starts | `dem/wbt/discha` is missing or incompatible. | Regenerate WBT distance-to-channel output or run without `--representative-flowpath`. |
| `watershed/README.md` missing after command exits | Abstraction did not complete successfully or output validation found an incomplete run. | Treat the run as failed; inspect logs and required inputs. |
| `flowpaths.parquet` missing | Flowpaths were skipped or representative mode was used, or full export failed. | Check generated manifest flags; require flowpath output only when `skip_flowpaths=false` and `representative_flowpath=false`. |
| Process killed by the OS or container runtime | Peak memory exceeded available limits, often during full flowpath expansion/export. | Re-run with `--skip-flowpaths` if consumers do not require flowpaths, or evaluate `--representative-flowpath` for WBT workflows. |
| Empty or implausible channels/hillslopes | Upstream delineation or network table mismatch. | Revalidate `subwta`, `flovec`, and network inputs before rerunning Peridot. |
| Sub-field CSVs missing | Field raster missing, no retained intersections, or output write failure. | Check `--field-raster`, threshold, and output directory permissions. |

## Current Error Boundary

The watershed CLI entrypoints return the underlying abstraction `io::Result<()>`. Propagated write-stage errors therefore return non-zero. Many missing-input cases still fail by panic because raster and network reads use explicit `unwrap()` calls; those failures also normally return non-zero.

Operational automation should combine process status with the post-run validation checks above. Missing required outputs remain authoritative failure evidence even if a caller did not capture process status.

## Deployment Notes for WEPPpy

After building release binaries, copy them into the WEPPpy Peridot binary directory expected by the deployment:

```bash
cp /workdir/peridot/target/release/abstract_watershed \
   /workdir/wepppy/wepppy/topo/peridot/bin/
cp /workdir/peridot/target/release/wbt_abstract_watershed \
   /workdir/wepppy/wepppy/topo/peridot/bin/
cp /workdir/peridot/target/release/sub_fields_abstraction \
   /workdir/wepppy/wepppy/topo/peridot/bin/
```

For Conda-specific binaries, use the deployment names expected by WEPPpy, for example `abstract_watershed.conda310.ub2404` and `wbt_abstract_watershed.conda310.ub2404`.

## Escalation Boundaries

Escalate to a runtime package, not a docs-only patch, when any of these are required:

- changing CLI exit-code behavior,
- renaming or removing output columns,
- changing representative-flowpath selection,
- adding new required outputs,
- changing deletion/mutation behavior for output directories,
- altering WEPPpy orchestration or queue behavior.
