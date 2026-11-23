// Validate the generated geometry against an independent flow-accumulation raster.
use std::{collections::HashMap, fs, path::Path, process::Command};

use peridot::raster::Raster;
use serde_json::Value;

#[test]
fn upstream_node_area_reaches_bieger_geometry_and_slope_file() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/wbt/amiss-polyhedron");
    let run = std::env::temp_dir().join(format!("peridot-bieger-{}", std::process::id()));
    fs::create_dir(&run).unwrap();
    let wbt = run.join("dem/wbt");
    fs::create_dir_all(&wbt).unwrap();
    for name in [
        "subwta.tif",
        "relief.tif",
        "flovec.tif",
        "fvslop.tif",
        "taspec.tif",
        "netw.tsv",
    ] {
        fs::copy(fixture.join("dem/wbt").join(name), wbt.join(name)).unwrap();
    }
    let result = Command::new(env!("CARGO_BIN_EXE_wbt_abstract_watershed"))
        .arg(&run)
        .args(["--bieger2015-widths", "--skip-flowpaths"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let accumulation =
        Raster::<f32>::read(fixture.join("dem/wbt/floaccum.tif").to_str().unwrap()).unwrap();
    let table = fs::read_to_string(wbt.join("netw.tsv")).unwrap();
    let mut expected = HashMap::new();
    for line in table.lines().skip(1) {
        let cols: Vec<_> = line.split('\t').collect();
        let row = cols[4].parse::<usize>().unwrap();
        let col = cols[5].parse::<usize>().unwrap();
        let upstream = accumulation.data[row * accumulation.width + col] as f64 - 1.0;
        let area = upstream * accumulation.cellsize.powi(2);
        let width = (1.24 * (area * 1e-6).powf(0.435)).max(0.324017);
        expected.insert(cols[1].parse::<i64>().unwrap(), (area, width));
    }
    let geojson: Value =
        serde_json::from_slice(&fs::read(run.join("watershed/channels.geojson")).unwrap()).unwrap();
    let reaches = geojson["features"].as_array().unwrap();
    assert_eq!(reaches.len(), 7);
    for feature in reaches {
        let properties = &feature["properties"];
        let id = properties["topaz_id"].as_i64().unwrap();
        let (area, width) = expected.remove(&id).unwrap();
        assert_eq!(properties["areaup"].as_f64().unwrap(), area);
        assert!((properties["width"].as_f64().unwrap() - width).abs() < 1e-10);
    }
    assert!(expected.is_empty());
    let slopes = fs::read_to_string(run.join("watershed/slope_files/channels.slp")).unwrap();
    // The outlet, TOPAZ 24, is the final record. Its 2,316 upstream cells give 1.7 m.
    let outlet = slopes.lines().rev().nth(2).unwrap();
    assert_eq!(outlet.split_whitespace().nth(1), Some("1.7"));
    fs::remove_dir_all(run).unwrap();
}
