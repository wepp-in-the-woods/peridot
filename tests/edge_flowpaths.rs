use std::collections::HashSet;
use std::path::PathBuf;

use peridot::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use peridot::raster::Raster;
use peridot::watershed_abstraction::walk_flowpaths;

fn fixture_wbt_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/wbt/amiss-polyhedron/dem/wbt")
}

#[test]
fn edge_flowpaths_match_source_cell_detection() {
    let fixture = fixture_wbt_path();

    let subwta = Raster::<i32>::read(fixture.join("subwta.tif").to_str().unwrap()).unwrap();
    let relief = Raster::<f32>::read(fixture.join("relief.tif").to_str().unwrap()).unwrap();
    let mut flovec = Raster::<u8>::read(fixture.join("flovec.tif").to_str().unwrap()).unwrap();
    let fvslop = Raster::<f32>::read(fixture.join("fvslop.tif").to_str().unwrap()).unwrap();
    let taspec = Raster::<f32>::read(fixture.join("taspec.tif").to_str().unwrap()).unwrap();

    remap_whitebox_d8_to_topaz_in_place(&mut flovec);

    let indices_map = subwta.indices_map();
    for (topaz_id, indices) in indices_map {
        if topaz_id <= 0 || topaz_id % 10 == 4 {
            continue;
        }

        let flowpaths = walk_flowpaths(
            topaz_id,
            &indices,
            &subwta,
            &relief,
            &flovec,
            &fvslop,
            &taspec,
        );

        #[allow(deprecated)]
        let edge_a = flowpaths.get_edge_flowpaths();
        let edge_b = flowpaths.get_edge_flowpaths2(&subwta, &flovec);

        let set_a: HashSet<usize> = edge_a
            .iter()
            .filter_map(|fp| fp.indices.first().copied())
            .collect();
        let set_b: HashSet<usize> = edge_b
            .iter()
            .filter_map(|fp| fp.indices.first().copied())
            .collect();

        assert_eq!(set_a, set_b, "topaz_id {}", topaz_id);
    }
}
