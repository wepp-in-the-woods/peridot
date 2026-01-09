extern crate clap;

use clap::Parser;
use gdal::errors::GdalError;
use log::info;
use rayon::ThreadPoolBuilder;
use std::path::Path;

use peridot::logging::init_logging;
use peridot::wbt_watershed_abstraction::wbt_abstract_watershed;

#[derive(Parser)]
#[command(
    version = env!("PERIDOT_VERSION_STRING"),
    long_version = env!("PERIDOT_VERSION_STRING"),
    disable_version_flag = true
)]
struct Opts {
    /// Path to the watershed directory
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

    /// Whether to calculate channel widths using Bieger 2015 fits based on uparea
    #[clap(short, long, default_value = "false")]
    bieger2015_widths: bool,

    /// Skip writing flowpaths outputs (flowpaths.csv and slope_files/flowpaths)
    #[clap(long, default_value = "false")]
    skip_flowpaths: bool,
}

fn main() -> Result<(), GdalError> {
    let opts: Opts = Opts::parse();
    init_logging(Path::new(&opts.path_to_wd));
    info!(
        "wbt_abstract_watershed start: wd={}, ncpu={}, max_points={}, clip_hillslopes={}, clip_hillslope_length={}, bieger2015_widths={}, skip_flowpaths={}",
        opts.path_to_wd,
        opts.ncpu,
        opts.max_points,
        opts.clip_hillslopes,
        opts.clip_hillslope_length,
        opts.bieger2015_widths,
        opts.skip_flowpaths
    );

    ThreadPoolBuilder::new()
        .num_threads(opts.ncpu)
        .build_global()
        .unwrap();

    let _ = wbt_abstract_watershed(
        &opts.path_to_wd,
        opts.max_points,
        opts.clip_hillslopes,
        opts.clip_hillslope_length,
        opts.bieger2015_widths,
        !opts.skip_flowpaths);

    Ok(())
}

// sudo -u www-data BACK_TRACE=1 ./abstract_watershed /geodata/weppcloud_runs/falling-validity/
