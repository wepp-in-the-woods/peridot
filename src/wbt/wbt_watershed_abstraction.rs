use std::env;
use std::fs;
use std::io::Result;
use std::mem::size_of;
use std::path::Path;
use std::time::Instant;

use log::info;
use rayon::prelude::*;

use crate::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use crate::flowpath_collection::FlowpathCollection;
use crate::netw::write_network;
use crate::raster::Raster;
use crate::watershed_abstraction::{abstract_subcatchments, walk_channels};
use crate::wbt_netw::read_wbt_netw_tab;

/// Find raster file, checking .tif first then .vrt as fallback
fn find_raster_path(base_path: &str) -> String {
    let tif_path = format!("{}.tif", base_path);
    if Path::new(&tif_path).exists() {
        return tif_path;
    }
    let vrt_path = format!("{}.vrt", base_path);
    if Path::new(&vrt_path).exists() {
        return vrt_path;
    }
    // Return .tif path to get original error message if neither exists
    tif_path
}

fn log_flowpath_collection_stats(label: &str, collection: &FlowpathCollection) {
    let flowpaths = &collection.flowpaths;
    let indices_total: usize = flowpaths.iter().map(|fp| fp.indices.len()).sum();
    let points_total: usize = flowpaths.iter().map(|fp| fp.distances_norm.len()).sum();
    let bytes = (indices_total * size_of::<usize>()) + (points_total * size_of::<f64>() * 5);
    info!(
        "{}: flowpaths={} indices_total={} points_total={} approx_bytes={} (~{:.2} MiB)",
        label,
        flowpaths.len(),
        indices_total,
        points_total,
        bytes,
        bytes as f64 / (1024.0 * 1024.0)
    );

    if let Some(subflows) = &collection.subflows {
        let mut subflow_paths = 0usize;
        let mut subflow_indices = 0usize;
        let mut subflow_points = 0usize;
        for sub in subflows.values() {
            subflow_paths += sub.flowpaths.len();
            subflow_indices += sub.flowpaths.iter().map(|fp| fp.indices.len()).sum::<usize>();
            subflow_points += sub.flowpaths.iter().map(|fp| fp.distances_norm.len()).sum::<usize>();
        }
        let subflow_bytes =
            (subflow_indices * size_of::<usize>()) + (subflow_points * size_of::<f64>() * 5);
        info!(
            "{} subflows: collections={} flowpaths={} indices_total={} points_total={} approx_bytes={} (~{:.2} MiB)",
            label,
            subflows.len(),
            subflow_paths,
            subflow_indices,
            subflow_points,
            subflow_bytes,
            subflow_bytes as f64 / (1024.0 * 1024.0)
        );
    }
}
#[allow(dead_code)]
pub fn wbt_abstract_watershed(
    wd: &str,
    max_points: usize,
    clip_hillslopes: bool,
    clip_hillslope_length: f64,
    bieger2015_widths: bool,
    write_flowpaths: bool
) -> std::io::Result<()> {
    env::set_current_dir(&wd).unwrap();
    info!(
        "wbt_abstract_watershed: wd={}, max_points={}, clip_hillslopes={}, clip_hillslope_length={}, bieger2015_widths={}, write_flowpaths={}",
        wd,
        max_points,
        clip_hillslopes,
        clip_hillslope_length,
        bieger2015_widths,
        write_flowpaths
    );

    let watershed_path = Path::new("watershed");

    if watershed_path.exists() {
        let _ = fs::remove_dir_all(&watershed_path).unwrap();
    }

    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/hillslopes"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/flowpaths"));

    let t0 = Instant::now();
    let subwta: Raster<i32> = Raster::<i32>::read(&find_raster_path("dem/wbt/subwta")).unwrap(); // hillslope with topaz_ids file, channels end with 4 (e.g 24, 34, 44), subcatchments end with 1, 2, 3. It starts at 22
    info!(
        "loaded subwta: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        subwta.width,
        subwta.height,
        subwta.data.len(),
        subwta.data.len() * size_of::<i32>(),
        (subwta.data.len() * size_of::<i32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let relief: Raster<f32> = Raster::<f32>::read(&find_raster_path("dem/wbt/relief")).unwrap(); // dem
    info!(
        "loaded relief: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        relief.width,
        relief.height,
        relief.data.len(),
        relief.data.len() * size_of::<f32>(),
        (relief.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let mut flovec: Raster<u8> = Raster::<u8>::read(&find_raster_path("dem/wbt/flovec")).unwrap(); // d8 flowvec
    remap_whitebox_d8_to_topaz_in_place(&mut flovec);
    info!(
        "loaded flovec: {}x{} cells={} bytes={} (~{:.2} MiB) and remapped in {:.2}s",
        flovec.width,
        flovec.height,
        flovec.data.len(),
        flovec.data.len() * size_of::<u8>(),
        (flovec.data.len() * size_of::<u8>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let fvslop: Raster<f32> = Raster::<f32>::read(&find_raster_path("dem/wbt/fvslop")).unwrap(); // slope
    info!(
        "loaded fvslop: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        fvslop.width,
        fvslop.height,
        fvslop.data.len(),
        fvslop.data.len() * size_of::<f32>(),
        (fvslop.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let taspec: Raster<f32> = Raster::<f32>::read(&find_raster_path("dem/wbt/taspec")).unwrap(); // aspect
    info!(
        "loaded taspec: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        taspec.width,
        taspec.height,
        taspec.data.len(),
        taspec.data.len() * size_of::<f32>(),
        (taspec.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let subwta_indices = subwta.indices_map();
    let channel_ids = subwta_indices.keys().filter(|&&id| id % 10 == 4).count();
    let hillslope_ids = subwta_indices.keys().filter(|&&id| id % 10 != 4).count();
    let indices_total: usize = subwta_indices.values().map(|v| v.len()).sum();
    info!(
        "built subwta indices map: ids={} (channels={}, hillslopes={}) total_cells={} (~{:.2} MiB) in {:.2}s",
        subwta_indices.len(),
        channel_ids,
        hillslope_ids,
        indices_total,
        (indices_total * size_of::<usize>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    // this is a ASCII tabular report with channel node connection informaiton
    let t0 = Instant::now();
    let (netw, network) = read_wbt_netw_tab("dem/wbt/netw.tsv").unwrap();
    let _ = write_network("watershed/network.txt", &network);
    info!(
        "read netw.tsv: links={} nodes={} in {:.2}s",
        netw.len(),
        network.len(),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let channels: FlowpathCollection = walk_channels(
        &subwta,
        &subwta_indices,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &netw,
        bieger2015_widths,
    );
    info!(
        "walked channels: flowpaths={} in {:.2}s",
        channels.flowpaths.len(),
        t0.elapsed().as_secs_f64()
    );
    log_flowpath_collection_stats("channels", &channels);

    let t0 = Instant::now();
    let hillslopes: FlowpathCollection = abstract_subcatchments(
        &subwta,
        &subwta_indices,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
        &channels,
        write_flowpaths,
    );
    let subflow_count = hillslopes
        .subflows
        .as_ref()
        .map(|m| m.len())
        .unwrap_or(0);
    info!(
        "abstracted subcatchments: hillslopes={} subflows={} in {:.2}s",
        hillslopes.flowpaths.len(),
        subflow_count,
        t0.elapsed().as_secs_f64()
    );
    log_flowpath_collection_stats("hillslopes", &hillslopes);

    let mut tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        Box::new(|| {
            info!("writing channels.slp");
            let result = channels.write_channel_slp("watershed/slope_files/channels.slp", max_points);
            info!("wrote channels.slp");
            result
        }),
        Box::new(|| {
            info!("writing channels.csv");
            let result = channels.write_chn_metadata_to_csv("watershed/channels.csv", &subwta.wgs_transform);
            info!("wrote channels.csv");
            result
        }),
        Box::new(|| {
            info!("writing hillslope slps");
            let result = hillslopes.write_slps("watershed/slope_files/hillslopes/", max_points, clip_hillslopes, clip_hillslope_length);
            info!("wrote hillslope slps");
            result
        }),
        Box::new(|| {
            info!("writing hillslopes.csv");
            let result = hillslopes.write_metadata_to_csv("watershed/hillslopes.csv", &subwta.wgs_transform);
            info!("wrote hillslopes.csv");
            result
        }),
        Box::new(|| {
            info!("writing channels.geojson");
            let result = channels.write_geojson(&subwta, "watershed/channels.geojson");
            info!("wrote channels.geojson");
            result
        }),
    ];

    if write_flowpaths {
        tasks.push(Box::new(|| {
            info!("writing flowpaths.csv");
            let result = hillslopes.write_subflows_metadata_to_csv("watershed/flowpaths.csv", &subwta.wgs_transform);
            info!("wrote flowpaths.csv");
            result
        }));
        tasks.push(Box::new(|| {
            info!("writing flowpath slps");
            let result = hillslopes.write_subflow_slps("watershed/slope_files/flowpaths/", max_points);
            info!("wrote flowpath slps");
            result
        }));
    }

    // Execute tasks in parallel
    info!("starting output tasks: {}", tasks.len());
    tasks.into_par_iter().map(|f| f()).collect::<Result<Vec<_>>>()?;
    info!("completed output tasks");

    Ok(())
}
