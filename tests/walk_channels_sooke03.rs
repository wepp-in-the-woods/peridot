use std::path::PathBuf;

use peridot::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use peridot::raster::Raster;
use peridot::watershed_abstraction::{walk_channels, MIN_CHANNEL_SLOPE_SCALAR};
use peridot::wbt_netw::read_wbt_netw_tab;

fn fixture_wbt_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/wbt/sooke03/dem/wbt")
}

#[test]
fn sooke03_flat_channels_do_not_panic() {
    let fixture = fixture_wbt_path();

    let subwta = Raster::<i32>::read(fixture.join("subwta.tif").to_str().unwrap()).unwrap();
    let relief = Raster::<f32>::read(fixture.join("relief.tif").to_str().unwrap()).unwrap();
    let mut flovec = Raster::<u8>::read(fixture.join("flovec.tif").to_str().unwrap()).unwrap();
    let fvslop = Raster::<f32>::read(fixture.join("fvslop.tif").to_str().unwrap()).unwrap();
    let taspec = Raster::<f32>::read(fixture.join("taspec.tif").to_str().unwrap()).unwrap();

    remap_whitebox_d8_to_topaz_in_place(&mut flovec);

    let (netw, _network) = read_wbt_netw_tab(fixture.join("netw.tsv")).unwrap();
    let subwta_indices = subwta.indices_map();

    let channels = walk_channels(
        &subwta,
        &subwta_indices,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &netw,
        true,
    );

    assert_eq!(channels.flowpaths.len(), 27);

    let mut fallback_channels = 0;
    let mut non_fallback_channels = 0;
    for flowpath in &channels.flowpaths {
        assert!(
            flowpath.slope_scalar.is_finite() && flowpath.slope_scalar > 0.0,
            "channel topaz_id {} has invalid slope_scalar: {}",
            flowpath.topaz_id,
            flowpath.slope_scalar
        );

        // Single-cell channels take the dedicated n == 1 branch in walk_channel().
        if flowpath.indices.len() == 1 {
            continue;
        }

        let total_elev = flowpath.elevs[0] - flowpath.elevs[flowpath.elevs.len() - 1];
        if total_elev <= 0.0 {
            assert!(
                (flowpath.slope_scalar - MIN_CHANNEL_SLOPE_SCALAR).abs() < 1e-12,
                "channel topaz_id {} expected fallback slope {}, got {}",
                flowpath.topaz_id,
                MIN_CHANNEL_SLOPE_SCALAR,
                flowpath.slope_scalar
            );
            fallback_channels += 1;
        } else {
            let expected = total_elev / flowpath.length;
            assert!(
                (flowpath.slope_scalar - expected).abs() < 1e-12,
                "channel topaz_id {} expected computed slope {}, got {}",
                flowpath.topaz_id,
                expected,
                flowpath.slope_scalar
            );
            non_fallback_channels += 1;
        }
    }

    assert!(
        fallback_channels > 0,
        "expected at least one channel to use MIN_CHANNEL_SLOPE_SCALAR"
    );
    assert!(
        non_fallback_channels > 0,
        "expected at least one multi-cell channel to use non-fallback computed slope"
    );
}
