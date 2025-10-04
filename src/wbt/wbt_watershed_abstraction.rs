use std::env;
use std::fs;
use std::io::Result;
use std::path::Path;

use rayon::prelude::*;

use crate::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz;
use crate::flowpaths::FlowpathCollection;
use crate::netw::write_network;
use crate::raster::Raster;
use crate::watershed_abstraction::{abstract_subcatchments, walk_channels};
use crate::wbt_netw::read_wbt_netw_tab;
#[allow(dead_code)]
pub fn wbt_abstract_watershed(
    wd: &str,
    max_points: usize,
    clip_hillslopes: bool,
    clip_hillslope_length: f64,
    bieger2015_widths: bool
) -> std::io::Result<()> {
    env::set_current_dir(&wd).unwrap();

    let watershed_path = Path::new("watershed");

    if watershed_path.exists() {
        let _ = fs::remove_dir_all(&watershed_path).unwrap();
    }

    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/hillslopes"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/flowpaths"));

    let subwta: Raster<i32> = Raster::<i32>::read("dem/wbt/subwta.tif").unwrap(); // hillslope with topaz_ids file, channels end with 4 (e.g 24, 34, 44), subcatchments end with 1, 2, 3. It starts at 22
    let relief: Raster<f64> = Raster::<f64>::read("dem/wbt/relief.tif").unwrap(); // dem
    let flovec_wbt: Raster<i32> = Raster::<i32>::read("dem/wbt/flovec.tif").unwrap(); // d8 flowvec
    let flovec = remap_whitebox_d8_to_topaz(&flovec_wbt);
    let fvslop: Raster<f64> = Raster::<f64>::read("dem/wbt/fvslop.tif").unwrap(); // slope
    let taspec: Raster<f64> = Raster::<f64>::read("dem/wbt/taspec.tif").unwrap(); // aspect

    // this is a ASCII tabular report with channel node connection informaiton
    let (netw, network) = read_wbt_netw_tab("dem/wbt/netw.tsv").unwrap();
    let _ = write_network("watershed/network.txt", &network);

    let channels: FlowpathCollection = walk_channels(&subwta, &relief, &flovec, &fvslop, &taspec, &netw, bieger2015_widths);
    let hillslopes: FlowpathCollection = abstract_subcatchments(&subwta, &relief, &flovec, &fvslop, &taspec, &channels);

    let tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        Box::new(|| channels.write_channel_slp("watershed/slope_files/channels.slp", max_points)),
        Box::new(|| channels.write_chn_metadata_to_csv("watershed/channels.csv", &subwta.wgs_transform)),
        Box::new(|| hillslopes.write_slps("watershed/slope_files/hillslopes/", max_points, clip_hillslopes, clip_hillslope_length)),
        Box::new(|| hillslopes.write_metadata_to_csv("watershed/hillslopes.csv", &subwta.wgs_transform)),
        Box::new(|| hillslopes.write_subflows_metadata_to_csv("watershed/flowpaths.csv", &subwta.wgs_transform)),
        Box::new(|| hillslopes.write_subflow_slps("watershed/slope_files/flowpaths/", max_points)),
        Box::new(|| channels.write_geojson(&subwta, "watershed/channels.geojson")),
    ];

    // Execute tasks in parallel
    tasks.into_par_iter().map(|f| f()).collect::<Result<Vec<_>>>()?;

    Ok(())
}
