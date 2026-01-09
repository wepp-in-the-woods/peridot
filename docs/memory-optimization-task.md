# Memory Optimization Task Plan

## Phase 0: Baseline and profiling
- [ ] Capture peak RSS and wall time for representative 1.0 m runs.
- [ ] Record raster sizes, counts, and output volume for comparison.

## Phase 1: Raster footprint reductions
- [x] Read relief, fvslop, and taspec as f32.
- [x] Read flovec as u8 and remap in place.
- [x] Apply the above changes in WBT sub-fields abstraction.
- [ ] Evaluate tiled or windowed GDAL reads for very large rasters.

## Phase 2: Index precomputation
- [x] Precompute topaz_id -> indices map for subwta.
- [x] Reuse indices map for WBT sub-fields intersection_subwta.
- [ ] Consider replacing HashSet-heavy paths with Vec-based indexing where safe.

## Phase 3: Flowpaths optional outputs
- [x] Add a CLI flag to skip flowpaths outputs.
- [x] Avoid storing subflow collections when flowpaths are disabled.
- [ ] Add a low-memory mode preset that combines multiple flags.

## Phase 4: Flowpath memory trimming
- [ ] Make slopes_r and distance_to_chn_r lazy or on-demand.
- [ ] Avoid cloning large index vectors where possible (centroid/aspect helpers).

## Phase 5: Output streaming and parallelism
- [ ] Stream channel SLP and GeoJSON writes instead of buffering full strings.
- [ ] Gate parallel output tasks for low-memory runs.
