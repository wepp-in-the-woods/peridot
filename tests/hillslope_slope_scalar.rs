use std::path::PathBuf;

use peridot::netw::read_netw_tab;
use peridot::raster::Raster;
use peridot::watershed_abstraction::{abstract_subcatchments, walk_channels};

fn fixture_topaz_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/watershed_abstraction/mdobre-scarce-belch/dem/topaz")
}

fn zonal_median_slope(fvslop: &Raster<f32>, indices: &Vec<usize>) -> f64 {
    let no_data = fvslop.no_data;
    let mut values: Vec<f64> = indices
        .iter()
        .filter_map(|&index| {
            let slope = fvslop.data[index];
            if !slope.is_finite() {
                return None;
            }
            if let Some(no_data_val) = no_data {
                if slope == no_data_val {
                    return None;
                }
            }
            Some(slope as f64)
        })
        .collect();

    assert!(
        !values.is_empty(),
        "no finite fvslop values for indices len={}",
        indices.len()
    );

    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n % 2 == 0 {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    } else {
        values[n / 2]
    }
}

fn is_side_mode(mode: &str) -> bool {
    matches!(
        mode,
        "side_edge_median_capped" | "side_area_over_channel" | "side_area_over_channel_no_edge"
    )
}

#[test]
fn nonrepresentative_hillslope_slope_scalar_matches_zonal_median() {
    let fixture = fixture_topaz_path();

    let subwta = Raster::<i32>::read(fixture.join("SUBWTA.ARC").to_str().unwrap()).unwrap();
    let relief = Raster::<f32>::read(fixture.join("RELIEF.ARC").to_str().unwrap()).unwrap();
    let flovec = Raster::<u8>::read(fixture.join("FLOVEC.ARC").to_str().unwrap()).unwrap();
    let fvslop = Raster::<f32>::read(fixture.join("FVSLOP.ARC").to_str().unwrap()).unwrap();
    let taspec = Raster::<f32>::read(fixture.join("TASPEC.ARC").to_str().unwrap()).unwrap();

    let (netw, _) = read_netw_tab(fixture.join("NETW.TAB").to_str().unwrap(), &subwta).unwrap();
    let indices_map = subwta.indices_map();

    let channels = walk_channels(
        &subwta,
        &indices_map,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &netw,
        false,
    );

    let hillslopes = abstract_subcatchments(
        &subwta,
        &indices_map,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &channels,
        false,
    );

    for fp in &hillslopes.flowpaths {
        let indices = indices_map
            .get(&fp.topaz_id)
            .expect("missing hillslope indices");
        let expected = zonal_median_slope(&fvslop, indices);
        assert!(
            (fp.slope_scalar - expected).abs() < 1e-12,
            "topaz_id={} expected={} got={}",
            fp.topaz_id,
            expected,
            fp.slope_scalar
        );
    }
}

#[test]
fn nonrepresentative_side_hillslope_length_modes_follow_selection_contract() {
    let fixture = fixture_topaz_path();

    let subwta = Raster::<i32>::read(fixture.join("SUBWTA.ARC").to_str().unwrap()).unwrap();
    let relief = Raster::<f32>::read(fixture.join("RELIEF.ARC").to_str().unwrap()).unwrap();
    let flovec = Raster::<u8>::read(fixture.join("FLOVEC.ARC").to_str().unwrap()).unwrap();
    let fvslop = Raster::<f32>::read(fixture.join("FVSLOP.ARC").to_str().unwrap()).unwrap();
    let taspec = Raster::<f32>::read(fixture.join("TASPEC.ARC").to_str().unwrap()).unwrap();

    let (netw, _) = read_netw_tab(fixture.join("NETW.TAB").to_str().unwrap(), &subwta).unwrap();
    let indices_map = subwta.indices_map();

    let channels = walk_channels(
        &subwta,
        &indices_map,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &netw,
        false,
    );

    let hillslopes = abstract_subcatchments(
        &subwta,
        &indices_map,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &channels,
        false,
    );

    let cellsize2 = subwta.cellsize * subwta.cellsize;
    for fp in &hillslopes.flowpaths {
        let indices = indices_map
            .get(&fp.topaz_id)
            .expect("missing hillslope indices");
        let expected_area = indices.len() as f64 * cellsize2;
        assert!(
            (fp.area_m2() - expected_area).abs() < 1e-6,
            "area mismatch topaz_id={} expected={} got={}",
            fp.topaz_id,
            expected_area,
            fp.area_m2()
        );

        match fp.topaz_id % 10 {
            1 => {
                assert_eq!(
                    fp.length_estimate_mode, "top_edge_median",
                    "top hillslope must keep edge-median mode"
                );
                assert!(
                    fp.length_edge_median.is_finite() && fp.length_edge_median > 0.0,
                    "top hillslope should carry finite edge-median candidate"
                );
                assert!(
                    fp.length_area_over_channel.is_nan(),
                    "top hillslope should not populate side area-over-channel candidate"
                );
            }
            2 | 3 => {
                assert!(
                    is_side_mode(&fp.length_estimate_mode),
                    "invalid side mode {} for topaz_id {}",
                    fp.length_estimate_mode,
                    fp.topaz_id
                );
                assert!(
                    fp.length_area_over_channel.is_finite() && fp.length_area_over_channel > 0.0,
                    "side hillslope must record area-over-channel candidate"
                );

                if fp.length_estimate_mode == "side_edge_median_capped" {
                    assert!(
                        fp.length_edge_median.is_finite() && fp.length_edge_median > 0.0,
                        "capped side hillslope must have finite edge-median candidate"
                    );
                    assert!((fp.length - fp.length_edge_median).abs() < 1e-6);
                    assert!(fp.length_edge_median < fp.length_area_over_channel);
                } else {
                    assert!((fp.length - fp.length_area_over_channel).abs() < 1e-6);
                }
            }
            _ => panic!("unexpected hillslope class for topaz_id={}", fp.topaz_id),
        }
    }
}
