# Bug: `walk_channel` panics on zero-elevation-drop channels

## Completion (2026-02-16)

Implemented in this repo:

1. Replaced the `assert!(total_elev > 0.0)` branch in `walk_channel()` with a graceful fallback.
2. Added `MIN_CHANNEL_SLOPE_SCALAR` (`1e-4`) and assign it when `total_elev <= 0.0`.
3. Added a warning log when fallback is used, including `topaz_id` and `total_elev`.
4. Added integration test `tests/walk_channels_sooke03.rs` that:
   - loads the Sooke03 WBT rasters,
   - runs `walk_channels()` without panic,
   - asserts 27 channels are returned,
   - verifies every channel has positive, finite `slope_scalar`,
   - verifies at least one channel uses `MIN_CHANNEL_SLOPE_SCALAR`.

Verification run:

- `cargo test --test walk_channels_sooke03` ✅
- `cargo test --test edge_flowpaths` ✅
- `cargo test` ⚠️ fails due pre-existing unrelated tests/environment issues (GDAL fixture opens and two panic-expectation tests).

## Problem

`wbt_abstract_watershed` crashes when a channel subcatchment has zero elevation
difference between its highest and lowest cells. The panic originates in
`walk_channel()` at `src/watershed_abstraction/watershed_abstraction.rs:336`:

```rust
let total_elev: f64 = elevs[0] - elevs[n - 1]; // in meters
assert!(total_elev > 0.0, "total_elev is 0.0 for topaz_id: {}", topaz_id);
slope_scalar = total_elev / distances[n - 1];
```

When every cell in a channel shares the same elevation (flat valley floor,
lakeshore, or a very small subcatchment whose cells round to the same value),
`total_elev` is 0.0 and the assert kills the process. Because the panic
happens inside a Rayon worker thread during `walk_channels()`, it aborts the
entire binary before any output files are written.

## Impact

The binary exits without producing `hillslopes.csv`, `channels.csv`, slope
files, or any other output. The Python caller (`peridot_runner.py`) does not
check the exit code, so the failure surfaces later as a misleading
`FileNotFoundError: hillslopes.csv not found`.

Three watersheds in the `victoria-2026-sbs` production batch failed this way:

| Watershed | Panicking topaz_ids   | Hillslopes | Channels |
|-----------|-----------------------|------------|----------|
| Leech     | 8154, 3704, 8834      | 2002       | 913      |
| Sooke20   | 294, 184              | 98         | 41       |
| Sooke03   | 214                   | 51         | 27       |

The `_peridot.log` for each shows the panic message followed by normal raster
loading (the panic text is flushed to the log out of order because it comes
from stderr on a worker thread).

## Reproduction

A fixture for Sooke03 has been added at:

```
tests/fixtures/wbt/sooke03/dem/wbt/
```

Before this fix, running `wbt_abstract_watershed` against this fixture would
reproduce the panic:

```
cargo run --bin wbt_abstract_watershed -- \
    tests/fixtures/wbt/sooke03 \
    --bieger2015-widths \
    --ncpu 1
```

After this fix, the same command is a verification run and should complete
without panic, producing watershed outputs.

## Requested fix

Replace the `assert!` with a graceful fallback that assigns a minimum
`slope_scalar` when `total_elev` is zero or negative. A flat channel is valid
geography and should not crash the binary.

Suggested approach in `walk_channel()` (line 334-338):

```rust
} else {
    let total_elev: f64 = elevs[0] - elevs[n - 1];
    if total_elev <= 0.0 {
        log::warn!(
            "total_elev is {:.6} for channel topaz_id: {}, assigning minimum slope",
            total_elev, topaz_id
        );
        slope_scalar = 1e-4; // minimum slope for flat channels
    } else {
        slope_scalar = total_elev / distances[n - 1];
    }
}
```

The `n == 1` branch (line 329-333) already handles single-cell channels by
falling back to `slopes[0]`. This fix extends the same spirit to multi-cell
channels that happen to be flat.

## Requested tests

Add an integration test using the `sooke03` fixture that:

1. Loads the Sooke03 rasters from `tests/fixtures/wbt/sooke03/dem/wbt/`
2. Calls `walk_channels()` and verifies it completes without panic
3. Asserts that the returned `FlowpathCollection` contains the expected number
   of channel flowpaths (27 channels per the `netw.tsv`)
4. Verifies that channels with `slope_scalar` set to the fallback value have
   a positive, finite `slope_scalar`

Example test skeleton:

```rust
#[test]
fn sooke03_flat_channels_do_not_panic() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/wbt/sooke03/dem/wbt");

    let subwta = Raster::<i32>::read(fixture.join("subwta.tif").to_str().unwrap()).unwrap();
    let relief = Raster::<f32>::read(fixture.join("relief.tif").to_str().unwrap()).unwrap();
    let mut flovec = Raster::<u8>::read(fixture.join("flovec.tif").to_str().unwrap()).unwrap();
    let fvslop = Raster::<f32>::read(fixture.join("fvslop.tif").to_str().unwrap()).unwrap();
    let taspec = Raster::<f32>::read(fixture.join("taspec.tif").to_str().unwrap()).unwrap();

    remap_whitebox_d8_to_topaz_in_place(&mut flovec);

    let (netw, _network) = read_wbt_netw_tab(
        fixture.parent().unwrap()  // sooke03/dem/wbt -> sooke03/dem
            .parent().unwrap()     // sooke03/dem -> sooke03
            .join("dem/wbt/netw.tsv")
            .to_str().unwrap()
    ).unwrap();

    let subwta_indices = subwta.indices_map();

    // This must not panic
    let channels = walk_channels(
        &subwta, &subwta_indices,
        &relief, &flovec, &fvslop, &taspec,
        &netw, true,
    );

    // Sooke03 has 27 channels
    assert_eq!(channels.flowpaths.len(), 27);

    // Every channel must have a positive, finite slope_scalar
    for fp in &channels.flowpaths {
        assert!(
            fp.slope_scalar > 0.0 && fp.slope_scalar.is_finite(),
            "channel topaz_id {} has invalid slope_scalar: {}",
            fp.topaz_id, fp.slope_scalar
        );
    }
}
```

A full end-to-end test that calls `wbt_abstract_watershed()` on the Sooke03
fixture and verifies that `watershed/hillslopes.csv` and
`watershed/channels.csv` are written would also be valuable but requires the
function to run in a temp directory (it does `env::set_current_dir`).
