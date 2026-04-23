use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use parquet::file::reader::{FileReader, SerializedFileReader};

use peridot::watershed_abstraction::{
    channels_summary, flowpaths_summary, hillslopes_summary, write_watershed_readme, Flowpath,
    FlowpathCollection, ManifestRunFlags,
};

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

fn sample_flowpath(topaz_id: i32, fp_id: i32, order: i32, width: f64) -> Flowpath {
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
        width,
        45.0,
        30.0,
        0.05,
        1.0,
        500.0,
        order,
        10.0,
    )
}

fn build_collections() -> (FlowpathCollection, FlowpathCollection) {
    let channels = FlowpathCollection {
        flowpaths: vec![sample_flowpath(14, -1, 1, 5.0)],
        subflows: None,
    };

    let mut subflows: HashMap<i32, FlowpathCollection> = HashMap::new();
    subflows.insert(
        11,
        FlowpathCollection {
            flowpaths: vec![
                sample_flowpath(11, 1, 1, 2.5),
                sample_flowpath(11, 2, 1, 2.0),
            ],
            subflows: None,
        },
    );

    let hillslopes = FlowpathCollection {
        flowpaths: vec![
            sample_flowpath(11, -1, 1, 4.0),
            sample_flowpath(12, -1, 1, 4.5),
        ],
        subflows: Some(subflows),
    };

    (channels, hillslopes)
}

fn parquet_rows(path: &Path) -> i64 {
    let file = fs::File::open(path).expect("failed to open parquet");
    let reader = SerializedFileReader::new(file).expect("failed to read parquet");
    reader.metadata().file_metadata().num_rows()
}

fn parquet_columns(path: &Path) -> Vec<String> {
    let file = fs::File::open(path).expect("failed to open parquet");
    let reader = SerializedFileReader::new(file).expect("failed to read parquet");
    reader
        .metadata()
        .file_metadata()
        .schema_descr()
        .columns()
        .iter()
        .map(|column| column.name().to_string())
        .collect()
}

#[test]
fn writes_watershed_parquet_outputs() {
    let wd = unique_temp_dir("parquet_outputs");
    let watershed_dir = wd.join("watershed");
    fs::create_dir_all(&watershed_dir).expect("failed to create watershed dir");

    let (channels, hillslopes) = build_collections();
    let wgs_transform = [-117.0, 46.0, 0.0001, 0.0001];

    channels
        .write_chn_metadata_to_parquet(
            watershed_dir.join("channels.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing channels parquet");
    hillslopes
        .write_metadata_to_parquet(
            watershed_dir.join("hillslopes.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing hillslopes parquet");
    hillslopes
        .write_subflows_metadata_to_parquet(
            watershed_dir.join("flowpaths.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing flowpaths parquet");

    assert_eq!(parquet_rows(&watershed_dir.join("channels.parquet")), 1);
    assert_eq!(parquet_rows(&watershed_dir.join("hillslopes.parquet")), 2);
    assert_eq!(parquet_rows(&watershed_dir.join("flowpaths.parquet")), 2);

    let channels_columns = parquet_columns(&watershed_dir.join("channels.parquet"));
    assert!(channels_columns.contains(&"topaz_id".to_string()));
    assert!(channels_columns.contains(&"order".to_string()));

    let flowpath_columns = parquet_columns(&watershed_dir.join("flowpaths.parquet"));
    assert!(flowpath_columns.contains(&"fp_id".to_string()));

    let hillslope_columns = parquet_columns(&watershed_dir.join("hillslopes.parquet"));
    assert!(hillslope_columns.contains(&"length_estimate_mode".to_string()));
    assert!(hillslope_columns.contains(&"length_area_over_channel".to_string()));
    assert!(hillslope_columns.contains(&"length_edge_median".to_string()));

    let _ = fs::remove_dir_all(wd);
}

#[test]
fn writes_watershed_readme_manifest_with_flags_and_schema() {
    let wd = unique_temp_dir("manifest");
    let watershed_dir = wd.join("watershed");
    fs::create_dir_all(watershed_dir.join("slope_files/hillslopes"))
        .expect("failed to create slope dir");

    let (channels, hillslopes) = build_collections();
    let wgs_transform = [-117.0, 46.0, 0.0001, 0.0001];

    channels
        .write_chn_metadata_to_parquet(
            watershed_dir.join("channels.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing channels parquet");
    hillslopes
        .write_metadata_to_parquet(
            watershed_dir.join("hillslopes.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing hillslopes parquet");
    hillslopes
        .write_subflows_metadata_to_parquet(
            watershed_dir.join("flowpaths.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing flowpaths parquet");

    fs::write(watershed_dir.join("network.txt"), "14|0\n").expect("failed to write network");
    fs::write(watershed_dir.join("channels.geojson"), "{}").expect("failed to write geojson");
    fs::write(
        watershed_dir.join("slope_files/hillslopes/hill_11.slp"),
        "demo slope file",
    )
    .expect("failed to write hill_11 slope file");
    fs::write(
        watershed_dir.join("slope_files/hillslopes/hill_12.slp"),
        "demo slope file",
    )
    .expect("failed to write hill_12 slope file");

    let tabular_outputs = vec![
        channels_summary(channels.flowpaths.len(), "parquet"),
        hillslopes_summary(hillslopes.flowpaths.len(), "parquet"),
        flowpaths_summary(hillslopes.subflow_row_count(), "parquet"),
    ];
    let flags = ManifestRunFlags {
        command: "wbt_abstract_watershed",
        max_points: 80,
        clip_hillslopes: true,
        clip_hillslope_length: 250.0,
        bieger2015_widths: true,
        skip_flowpaths: false,
        representative_flowpath: false,
    };

    write_watershed_readme(&watershed_dir, &flags, &tabular_outputs)
        .expect("failed writing watershed README");

    let readme =
        fs::read_to_string(watershed_dir.join("README.md")).expect("failed to read README");

    assert!(readme.contains("clip_hillslopes"));
    assert!(readme.contains("clip_hillslope_length"));
    assert!(readme.contains("bieger2015_widths"));
    assert!(readme.contains("skip_flowpaths"));
    assert!(readme.contains("representative_flowpath"));

    assert!(readme.contains("watershed/hillslopes.parquet"));
    assert!(readme.contains("watershed/channels.parquet"));
    assert!(readme.contains("watershed/flowpaths.parquet"));
    assert!(readme.contains("watershed/slope_files/hillslopes/*"));
    assert!(readme.contains("slp bundle"));
    assert!(readme.contains("2 files"));
    assert!(!readme.contains("watershed/slope_files/hillslopes/hill_11.slp"));
    assert!(readme.contains("| topaz_id | int32 |"));
    assert!(readme.contains("| fp_id | int32 |"));
    assert!(readme.contains("| length_estimate_mode | utf8 |"));
    assert!(readme.contains("| length_area_over_channel | float64 |"));
    assert!(readme.contains("| length_edge_median | float64 |"));
    assert!(readme.contains(
        "Flowpath tabular and slope outputs are expected because `skip_flowpaths=false`."
    ));

    let _ = fs::remove_dir_all(wd);
}

#[test]
fn writes_watershed_readme_conditional_notes_for_skip_and_representative_modes() {
    let wd = unique_temp_dir("manifest_conditional_notes");
    let watershed_dir = wd.join("watershed");
    fs::create_dir_all(watershed_dir.join("slope_files/hillslopes"))
        .expect("failed to create slope dir");

    let (channels, hillslopes) = build_collections();
    let wgs_transform = [-117.0, 46.0, 0.0001, 0.0001];

    channels
        .write_chn_metadata_to_parquet(
            watershed_dir.join("channels.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing channels parquet");
    hillslopes
        .write_metadata_to_parquet(
            watershed_dir.join("hillslopes.parquet").to_str().unwrap(),
            &wgs_transform,
        )
        .expect("failed writing hillslopes parquet");

    let tabular_outputs = vec![
        channels_summary(channels.flowpaths.len(), "parquet"),
        hillslopes_summary(hillslopes.flowpaths.len(), "parquet"),
    ];
    let flags = ManifestRunFlags {
        command: "wbt_abstract_watershed",
        max_points: 64,
        clip_hillslopes: false,
        clip_hillslope_length: 300.0,
        bieger2015_widths: false,
        skip_flowpaths: true,
        representative_flowpath: true,
    };

    write_watershed_readme(&watershed_dir, &flags, &tabular_outputs)
        .expect("failed writing watershed README");

    let readme =
        fs::read_to_string(watershed_dir.join("README.md")).expect("failed to read README");

    assert!(
        readme.contains("Flowpath tabular/slope outputs are skipped because `skip_flowpaths=true`")
    );
    assert!(readme.contains(
        "`representative_flowpath=true` forces a single representative hillslope flowpath summary"
    ));
    assert!(readme.contains("Hillslope slope profiles are not clipped (`clip_hillslopes=false`)."));
    assert!(readme.contains("Channel widths use the legacy width model."));
    assert!(!readme.contains("watershed/flowpaths.parquet"));

    let _ = fs::remove_dir_all(wd);
}
