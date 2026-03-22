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
