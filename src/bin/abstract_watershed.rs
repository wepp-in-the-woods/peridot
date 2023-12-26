extern crate clap;

use gdal::errors::GdalError;
use rayon::ThreadPoolBuilder;
use clap::Parser;

use peridot::watershed_abstraction::abstract_watershed;

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
}

fn main() -> Result<(), GdalError> {
    let opts: Opts = Opts::parse();

    ThreadPoolBuilder::new()
        .num_threads(opts.ncpu)
        .build_global()
        .unwrap();

    let _ = abstract_watershed(
        &opts.path_to_wd, opts.max_points, opts.clip_hillslopes, opts.clip_hillslope_length);

    Ok(())
}
