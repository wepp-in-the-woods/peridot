use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use peridot::watershed_abstraction::{Flowpath, FlowpathCollection};

fn unique_temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock before epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "peridot_{}_{}_{}",
        label,
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("failed to create temp dir");
    dir
}

fn sample_flowpath(topaz_id: i32, fp_id: i32) -> Flowpath {
    Flowpath::new(
        vec![0, 1, 2],
        (0, 0),
        (1, 1),
        (10, 11),
        vec![0.0, 0.5, 1.0],
        vec![0.1, 0.08, 0.05],
        vec![500.0, 495.0, 490.0],
        topaz_id,
        fp_id,
        100.0,
        2.5,
        45.0,
        30.0,
        0.05,
        1.0,
        500.0,
        1,
        10.0,
    )
}

#[test]
fn field_flowpaths_csv_headers_are_unique_and_explicit() {
    let wd = unique_temp_dir("field_flowpaths_schema");
    let csv_path = wd.join("field_flowpaths.csv");
    let fake_topaz_id = 9001;

    let mut subflows = HashMap::new();
    subflows.insert(
        fake_topaz_id,
        FlowpathCollection {
            flowpaths: vec![sample_flowpath(11, 1)],
            subflows: None,
        },
    );
    let collection = FlowpathCollection {
        flowpaths: Vec::new(),
        subflows: Some(subflows),
    };
    let lookup = HashMap::from([((7, 11), fake_topaz_id)]);

    collection
        .write_field_subflows_metadata_to_csv(
            csv_path.to_str().unwrap(),
            &[-117.0, 46.0, 0.1, 0.1],
            &lookup,
        )
        .expect("failed to write field flowpaths csv");

    let mut reader = csv::Reader::from_path(&csv_path).expect("failed to read csv");
    let headers = reader.headers().expect("missing headers").clone();
    let header_values: Vec<&str> = headers.iter().collect();
    let unique: HashSet<&str> = header_values.iter().copied().collect();

    assert_eq!(header_values.len(), unique.len(), "headers must be unique");
    assert_eq!(
        &header_values[..5],
        [
            "field_id",
            "topaz_id",
            "sub_field_id",
            "flowpath_topaz_id",
            "fp_id"
        ]
    );
    assert_eq!(
        header_values
            .iter()
            .filter(|&&name| name == "topaz_id")
            .count(),
        1
    );
    assert_eq!(
        header_values
            .iter()
            .filter(|&&name| name == "flowpath_topaz_id")
            .count(),
        1
    );

    let _ = fs::remove_dir_all(wd);
}
