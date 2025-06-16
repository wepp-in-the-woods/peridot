extern crate clap;

use gdal::errors::GdalError;
use rayon::ThreadPoolBuilder;
use clap::Parser;

use peridot::watershed_abstraction::wbt_abstract_watershed;

#[derive(Parser)]
struct Opts {
    /// Path to the watershed directory
    path_to_wd: String,

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
}

fn main() -> Result<(), GdalError> {
    let opts: Opts = Opts::parse();

    ThreadPoolBuilder::new()
        .num_threads(opts.ncpu)
        .build_global()
        .unwrap();

    let _ = wbt_abstract_watershed(
        &opts.path_to_wd,
        opts.max_points,
        opts.clip_hillslopes,
        opts.clip_hillslope_length,
        opts.bieger2015_widths);

    Ok(())
}

// sudo -u www-data BACK_TRACE=1 ./abstract_watershed /geodata/weppcloud_runs/falling-validity/