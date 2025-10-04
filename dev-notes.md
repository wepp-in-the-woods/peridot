# Developer Notes

## Code layout
- `src/bin/abstract_watershed.rs` and `src/bin/wbt_abstract_watershed.rs` parse CLI flags with `clap`, configure the Rayon thread pool, then call into `peridot::watershed_abstraction`.
- `src/watershed_abstraction.rs` holds the core domain logic: raster IO, flow path walking, slope aggregation, and report writers.
- Supporting modules live alongside it:
  - `src/raster.rs` wraps GDAL datasets behind a typed `Raster<T>` abstraction.
  - `src/netw.rs` and `src/wbt_netw.rs` parse the channel connectivity tables produced by TOPAZ and Whitebox Tools respectively.
  - `src/whiteboxtools_wrappers.rs` reconciles D8 encoding differences between Whitebox Tools and TOPAZ.
  - `src/support.rs` and `src/douglas_peucker.rs` provide geometric helpers used during slope profile synthesis.

## Execution flow
1. The binary switches `cwd` to the run directory supplied on the command line. **All paths are relative from there**, so tests typically stage fixtures under `tests/data/<scenario>/`.
2. Existing `watershed/` output is purged to guarantee clean regeneration.
3. Rasters are loaded from either `dem/topaz/` (`*.ARC`) or `dem/wbt/` (`*.tif`) into `Raster<T>` instances. Metadata (cell size, transforms) travels with the data.
4. Channel connectivity is read via `read_netw_tab` (TOPAZ) or `read_wbt_netw_tab` (Whitebox). The result is a vector of network nodes plus an adjacency list dumped to `watershed/network.txt` for debugging.
5. `walk_channels` traverses the stream network, recording the polyline geometry, slope, and width metadata for each channel. Widths honour the Bieger (2015) regression when requested.
6. `abstract_subcatchments` iterates over each hillslope ID, hydrologically samples internal flow paths (`walk_flowpaths`), then collapses them into a representative `FlowPath` using weighted slope averaging (`FlowpathCollection::weighted_slope_average_from_fps`).
7. The final `FlowpathCollection` objects write slope files (`*.slp`), CSV tables, and GeoJSON in parallel using Rayon tasks.

## Data expectations
- TOPAZ inputs must include `SUBWTA.ARC`, `RELIEF.ARC`, `FLOVEC.ARC`, `FVSLOP.ARC`, `TASPEC.ARC`, and `NETW.TAB` under `dem/topaz/`.
- Whitebox inputs must include `subwta.tif`, `relief.tif`, `flovec.tif`, `fvslop.tif`, `taspec.tif`, and `netw.tsv` under `dem/wbt/`.
- Rasters share a common grid; mismatched dimensions will panic early when `Raster::read` performs shape checks.
- `FlowPath.topaz_id` mirrors TOPAZ's numbering (channels end in `4`, contributing hillslopes end in `1`, `2`, or `3`). Whitebox runs adopt the same convention after the D8 remap step.

## Developing locally
- Use `cargo test` to execute the regression suite; tests live under `tests/` and stage watershed fixtures for both TOPAZ and Whitebox scenarios.
- Enable logging with `RUST_LOG=debug` when running binaries to trace flowpath decisions (the core module uses `log::debug!` extensively).
- Long-running operations are parallelised with Rayon; prefer extending the existing task vectors instead of spawning threads manually.
- When adjusting raster IO, keep the ASCII-only requirement in mind: the slope file writers assume US-ASCII output compatible with WEPP importers.
- The repo ships with `set_wepppy310_env.sh` to activate the Conda toolchain used by WEPPcloud deployments; source it before building if you rely on that environment.

## Extending the abstraction tools
- Share new per-hillslope metrics by extending `FlowPath` in `src/watershed_abstraction.rs`, then thread the data through the CSV and GeoJSON writers.
- Add backend-specific behaviour using feature gates inside `wbt_abstract_watershed` or `abstract_watershed`. The shared code path makes it easy to keep the outputs aligned—avoid duplicating logic between the two entry points.
- If you need additional rasters, follow the pattern in `Raster::<T>::read` to guarantee GDAL closes datasets deterministically and honours the `wgs_transform` metadata.
