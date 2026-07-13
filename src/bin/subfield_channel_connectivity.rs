extern crate clap;

use clap::Parser;
use serde::Serialize;
use std::error::Error;
use std::fs;

use peridot::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use peridot::raster::Raster;
use peridot::subfield_channel_connectivity::{
    summarize_subfield_channel_connectivity, SubfieldChannelConnectivitySummary,
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

    let metrics = summarize_subfield_channel_connectivity(
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
        metrics,
    };

    let json = serde_json::to_string_pretty(&report)?;
    if let Some(path) = opts.out_json {
        fs::write(path, format!("{}\n", json))?;
    } else {
        println!("{}", json);
    }

    Ok(())
}
