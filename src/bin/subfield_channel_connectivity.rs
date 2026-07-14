extern crate clap;

use clap::Parser;
use serde::Serialize;
use std::error::Error;
use std::fs;
use std::path::Path;

use peridot::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use peridot::raster::Raster;
use peridot::subfield_channel_connectivity::{
    analyze_subfield_channel_connectivity, SubfieldChannelConnectivityDetail,
    SubfieldChannelConnectivitySummary,
};

const DEFINITION: &str = "A retained subfield has direct channel drainage when the first cell outside at least one generated per-cell flowpath is a channel cell. Peridot starts one flowpath at every retained subfield cell and stops after appending its first cell outside the subfield.";

#[derive(Parser)]
#[command(
    version = env!("PERIDOT_VERSION_STRING"),
    long_version = env!("PERIDOT_VERSION_STRING"),
    disable_version_flag = true
)]
struct Opts {
    /// Path to the retained subfield ID raster produced by sub_fields_abstraction
    #[arg(long)]
    sub_field_map: String,

    /// Path to the WBT SUBWTA raster; channel IDs ending in 4 are used by default
    #[arg(long)]
    subwta: String,

    /// Path to the raw WhiteboxTools D8 flow-direction raster
    #[arg(long)]
    wbt_flovec: String,

    /// Optional channel mask raster; positive cells override the SUBWTA suffix rule
    #[arg(long)]
    channel_mask: Option<String>,

    /// Optional output JSON path. If omitted, JSON is printed to stdout.
    #[arg(long)]
    out_json: Option<String>,

    /// Optional versioned per-subfield routing detail JSON path
    #[arg(long)]
    out_subfields_json: Option<String>,

    /// Show version information and exit
    #[clap(short = 'v', long = "version", action = clap::ArgAction::Version, short_alias = 'V')]
    _version: Option<bool>,
}

#[derive(Serialize)]
struct InputResources<'a> {
    sub_field_map: &'a str,
    subwta: &'a str,
    wbt_flovec: &'a str,
    channel_mask: Option<&'a str>,
}

#[derive(Serialize)]
struct ConnectivityReport<'a> {
    schema_version: u8,
    definition: &'static str,
    channel_detection: &'static str,
    inputs: InputResources<'a>,
    metrics: SubfieldChannelConnectivitySummary,
}

#[derive(Serialize)]
struct SubfieldConnectivityReport<'a> {
    schema_version: u8,
    peridot_version: &'static str,
    definition: &'static str,
    channel_detection: &'static str,
    inputs: InputResources<'a>,
    metrics: SubfieldChannelConnectivitySummary,
    subfields: Vec<SubfieldChannelConnectivityDetail>,
}

fn write_atomic(path: &Path, contents: &str) -> Result<(), Box<dyn Error>> {
    let file_name = path.file_name().ok_or("output path has no file name")?;
    let temporary = path.with_file_name(format!(
        ".{}.{}.tmp",
        file_name.to_string_lossy(),
        std::process::id()
    ));
    fs::write(&temporary, contents)?;
    match fs::rename(&temporary, path) {
        Ok(()) => Ok(()),
        Err(error) => {
            let _ = fs::remove_file(&temporary);
            Err(error.into())
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let opts = Opts::parse();

    let sub_field_map = Raster::<i32>::read(&opts.sub_field_map)?;
    let subwta = Raster::<i32>::read(&opts.subwta)?;
    let mut flovec = Raster::<u8>::read(&opts.wbt_flovec)?;
    remap_whitebox_d8_to_topaz_in_place(&mut flovec);
    let channel_mask = opts
        .channel_mask
        .as_deref()
        .map(Raster::<i32>::read)
        .transpose()?;

    let analysis = analyze_subfield_channel_connectivity(
        &sub_field_map,
        &subwta,
        &flovec,
        channel_mask.as_ref(),
    )?;
    let channel_detection = if channel_mask.is_some() {
        "positive channel_mask cells"
    } else {
        "SUBWTA IDs whose final decimal digit is 4"
    };
    let report = ConnectivityReport {
        schema_version: 1,
        definition: DEFINITION,
        channel_detection,
        inputs: InputResources {
            sub_field_map: &opts.sub_field_map,
            subwta: &opts.subwta,
            wbt_flovec: &opts.wbt_flovec,
            channel_mask: opts.channel_mask.as_deref(),
        },
        metrics: analysis.summary.clone(),
    };

    let json = serde_json::to_string_pretty(&report)?;
    if let Some(path) = opts.out_json {
        write_atomic(Path::new(&path), &format!("{}\n", json))?;
    } else {
        println!("{}", json);
    }

    if let Some(path) = opts.out_subfields_json {
        let detail_report = SubfieldConnectivityReport {
            schema_version: 1,
            peridot_version: env!("PERIDOT_VERSION_STRING"),
            definition: DEFINITION,
            channel_detection,
            inputs: InputResources {
                sub_field_map: &opts.sub_field_map,
                subwta: &opts.subwta,
                wbt_flovec: &opts.wbt_flovec,
                channel_mask: opts.channel_mask.as_deref(),
            },
            metrics: analysis.summary,
            subfields: analysis.subfields,
        };
        let json = serde_json::to_string_pretty(&detail_report)?;
        write_atomic(Path::new(&path), &format!("{}\n", json))?;
    }

    Ok(())
}
