extern crate clap;

use clap::Parser;
use std::error::Error;
use std::fs;

use peridot::raster::Raster;
use peridot::roads_trace::trace_downslope_flowpath;

#[derive(Parser)]
#[command(
    version = env!("PERIDOT_VERSION_STRING"),
    long_version = env!("PERIDOT_VERSION_STRING"),
    disable_version_flag = true
)]
struct Opts {
    /// Path to SUBWTA raster
    #[arg(long)]
    subwta: String,

    /// Path to FLOVEC raster
    #[arg(long)]
    flovec: String,

    /// Path to RELIEF raster
    #[arg(long)]
    relief: String,

    /// Optional path to channel mask raster
    #[arg(long)]
    channel: Option<String>,

    /// Seed row (0-based)
    #[arg(long)]
    seed_row: usize,

    /// Seed col (0-based)
    #[arg(long)]
    seed_col: usize,

    /// Maximum D8 steps before termination
    #[arg(long, default_value_t = 20000)]
    max_steps: usize,

    /// Optional output JSON path. If omitted, JSON is printed to stdout.
    #[arg(long)]
    out_json: Option<String>,

    /// Show version information and exit
    #[clap(short = 'v', long = "version", action = clap::ArgAction::Version, short_alias = 'V')]
    _version: Option<bool>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let opts = Opts::parse();

    let subwta = Raster::<i32>::read(&opts.subwta)?;
    let flovec = Raster::<u8>::read(&opts.flovec)?;
    let relief = Raster::<f32>::read(&opts.relief)?;

    let channel_mask = match opts.channel {
        Some(path) => Some(Raster::<i32>::read(&path)?),
        None => None,
    };

    let result = trace_downslope_flowpath(
        &subwta,
        &flovec,
        &relief,
        opts.seed_row,
        opts.seed_col,
        channel_mask.as_ref(),
        opts.max_steps,
    )?;

    let json = serde_json::to_string_pretty(&result)?;

    if let Some(path) = opts.out_json {
        fs::write(path, format!("{}\n", json))?;
    } else {
        println!("{}", json);
    }

    Ok(())
}
