extern crate clap;

use clap::Parser;
use log::info;
use rayon::ThreadPoolBuilder;
use std::io;
use std::path::Path;

use peridot::logging::init_logging;
use peridot::wbt_sub_fields_abstraction::wbt_sub_fields_abstraction;

#[derive(Parser)]
#[command(
    version = env!("PERIDOT_VERSION_STRING"),
    long_version = env!("PERIDOT_VERSION_STRING"),
    disable_version_flag = true
)]
struct Opts {
    /// Path to the WEPPcloud run directory
    path_to_wd: String,

    /// Show version information and exit
    #[clap(short = 'v', long = "version", action = clap::ArgAction::Version, short_alias = 'V')]
    _version: Option<bool>,

    /// Number of CPU threads
    #[clap(short, long, default_value = "4")]
    ncpu: usize,

    /// Maximum number of points
    #[clap(short, long, default_value = "99")]
    max_points: usize,

    /// Whether to clip hillslopes or not
    #[clap(short, long, default_value = "false")]
    clip_hillslopes: bool,

    /// Clip hillslope length
    #[clap(long, default_value = "300.0")]
    clip_hillslope_length: f64,

    /// Minimum area threshold (m^2) for retaining sub-fields
    #[clap(long = "sub-field-min-area-threshold-m2", default_value = "0.0")]
    sub_field_min_area_threshold_m2: f64,

    /// Path to the rasterised field boundaries relative to the run directory
    #[clap(long = "field-raster", default_value = "ag_fields/field_boundaries.tif")]
    field_raster: String,

    /// Output directory for sub-field artifacts relative to the run directory
    #[clap(long = "output-dir", default_value = "ag_fields/sub_fields")]
    output_dir: String,
}

fn main() -> io::Result<()> {
    let opts: Opts = Opts::parse();
    init_logging(Path::new(&opts.path_to_wd));
    info!(
        "sub_fields_abstraction start: wd={}, ncpu={}, max_points={}, clip_hillslopes={}, clip_hillslope_length={}, sub_field_min_area_threshold_m2={}, field_raster={}, output_dir={}",
        opts.path_to_wd,
        opts.ncpu,
        opts.max_points,
        opts.clip_hillslopes,
        opts.clip_hillslope_length,
        opts.sub_field_min_area_threshold_m2,
        opts.field_raster,
        opts.output_dir
    );

    ThreadPoolBuilder::new()
        .num_threads(opts.ncpu)
        .build_global()
        .expect("Failed to initialise Rayon thread pool");

    wbt_sub_fields_abstraction(
        &opts.path_to_wd,
        opts.max_points,
        opts.clip_hillslopes,
        opts.clip_hillslope_length,
        opts.sub_field_min_area_threshold_m2,
        &opts.field_raster,
        &opts.output_dir,
    )
}
