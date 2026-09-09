use std::collections::HashMap;
use std::fs::{self, File};
use std::time::{SystemTime, UNIX_EPOCH};

use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use peridot::raster::{MapType, PixelToWgs84, Raster};
use peridot::watershed_abstraction::{Flowpath, FlowpathCollection};

fn raster() -> Raster<i32> {
    Raster::new(
        2044,
        1984,
        10.0,
        vec![],
        None,
        [572803.0, 10.0, 0.0, 5038749.0, 0.0, -10.0],
        Some("EPSG:32610".into()),
        String::new(),
        String::new(),
        MapType::SUBWTA,
    )
}

fn assert_location(actual: (f64, f64), expected: (f64, f64)) {
    assert!(
        (actual.0 - expected.0).abs() < 1e-9,
        "longitude: {actual:?}"
    );
    assert!((actual.1 - expected.1).abs() < 1e-9, "latitude: {actual:?}");
}

#[test]
fn projects_utm_centroids_without_two_corner_approximation() {
    let raster = raster();
    // Independent pyproj reference using the existing pixel-corner convention.
    let expected = (-122.01533619935442, 45.474715610246655);
    let projector = raster.centroid_projector().unwrap();
    assert_location(projector.convert(416, 259).unwrap(), expected);
    // The incident's old approximation differs by over 100 m on each axis.
    let approximate = peridot::raster::px_to_wgs(&raster.wgs_transform, 416, 259);
    assert!((approximate.0 - expected.0).abs() > 0.001);
    assert!((approximate.1 - expected.1).abs() > 0.001);
    assert_location(projector.convert(416, 259).unwrap(), expected);
}

#[test]
fn projects_southern_hemisphere_and_rotated_affine() {
    let south =
        PixelToWgs84::new([650000.0, 30.0, 0.0, 6200000.0, 0.0, -30.0], "EPSG:32755").unwrap();
    assert_location(
        south.convert(913, 715).unwrap(),
        (148.93260146054337, -34.51945229497312),
    );
    let rotated =
        PixelToWgs84::new([572803.0, 10.0, 2.0, 5038749.0, 1.0, -10.0], "EPSG:32610").unwrap();
    assert_location(
        rotated.convert(416, 259).unwrap(),
        (-122.00864425551065, 45.47840221955142),
    );
}

#[test]
fn projection_errors_are_explicit() {
    let mut source = raster();
    source.proj4 = None;
    assert!(source.centroid_projector().is_err());
    assert!(PixelToWgs84::new(source.geo_transform, "invalid-crs").is_err());
    assert!(PixelToWgs84::new(source.geo_transform, "").is_err());
    let mut nonfinite = source.geo_transform;
    nonfinite[0] = f64::NAN;
    assert!(PixelToWgs84::new(nonfinite, "EPSG:32610").is_err());
    // Use an invalid UTM coordinate to exercise PROJ's conversion error path.
    let invalid_utm = PixelToWgs84::new([1e100, 1.0, 0.0, 1e100, 0.0, 1.0], "EPSG:32610").unwrap();
    assert!(invalid_utm.convert(0, 0).is_err());
}

fn collection() -> FlowpathCollection {
    let fp = Flowpath::new(
        vec![0, 1, 2],
        (0, 0),
        (1, 1),
        (416, 259),
        vec![0.0, 0.5, 1.0],
        vec![0.1, 0.08, 0.05],
        vec![500.0, 495.0, 490.0],
        202,
        1,
        100.0,
        2.5,
        45.0,
        30.0,
        0.05,
        1.0,
        500.0,
        1,
        10.0,
    );
    FlowpathCollection {
        flowpaths: vec![fp.clone()],
        subflows: Some(HashMap::from([(
            202,
            FlowpathCollection {
                flowpaths: vec![fp],
                subflows: None,
            },
        )])),
    }
}

#[test]
fn every_metadata_export_uses_projected_centroid() {
    let dir = std::env::temp_dir().join(format!(
        "peridot_projection_{}_{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&dir).unwrap();
    let r = raster();
    let c = collection();
    let lookup = HashMap::from([((7, 11), 202)]);
    let path = |name: &str| dir.join(name).to_str().unwrap().to_owned();
    c.write_chn_metadata_to_parquet(&path("channels.parquet"), &r)
        .unwrap();
    c.write_metadata_to_parquet(&path("hillslopes.parquet"), &r)
        .unwrap();
    c.write_subflows_metadata_to_parquet(&path("flowpaths.parquet"), &r)
        .unwrap();
    c.write_chn_metadata_to_csv(&path("channels.csv"), &r)
        .unwrap();
    c.write_metadata_to_csv(&path("hillslopes.csv"), &r)
        .unwrap();
    c.write_subflows_metadata_to_csv(&path("flowpaths.csv"), &r)
        .unwrap();
    c.write_field_metadata_to_csv(&path("fields.csv"), &r, &lookup)
        .unwrap();
    c.write_field_subflows_metadata_to_csv(&path("field_flowpaths.csv"), &r, &lookup)
        .unwrap();
    let expected = (-122.01533619935442, 45.474715610246655);
    for name in ["channels", "hillslopes", "flowpaths"] {
        let reader =
            SerializedFileReader::new(File::open(path(&format!("{name}.parquet"))).unwrap())
                .unwrap();
        let row = reader.get_row_iter(None).unwrap().next().unwrap().unwrap();
        let columns = reader.metadata().file_metadata().schema_descr().columns();
        let lon = columns
            .iter()
            .position(|col| col.name() == "centroid_lon")
            .unwrap();
        let lat = columns
            .iter()
            .position(|col| col.name() == "centroid_lat")
            .unwrap();
        assert_location(
            (row.get_double(lon).unwrap(), row.get_double(lat).unwrap()),
            expected,
        );
    }
    for name in [
        "channels",
        "hillslopes",
        "flowpaths",
        "fields",
        "field_flowpaths",
    ] {
        let mut reader = csv::Reader::from_path(path(&format!("{name}.csv"))).unwrap();
        let headers = reader.headers().unwrap();
        let lon = headers.iter().position(|s| s == "centroid_lon").unwrap();
        let lat = headers.iter().position(|s| s == "centroid_lat").unwrap();
        let row = reader.records().next().unwrap().unwrap();
        assert_location(
            (row[lon].parse().unwrap(), row[lat].parse().unwrap()),
            expected,
        );
    }
    let mut missing_crs = r;
    missing_crs.proj4 = None;
    assert!(c
        .write_metadata_to_parquet(&path("invalid.parquet"), &missing_crs)
        .is_err());
    assert!(!dir.join("invalid.parquet").exists());
    fs::remove_dir_all(dir).unwrap();
}
