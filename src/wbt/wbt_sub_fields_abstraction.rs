use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::io::Result;
use std::mem::size_of;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use log::info;
use rayon::prelude::*;

use crate::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use crate::flowpath::Flowpath;
use crate::flowpath_collection::FlowpathCollection;
use crate::raster::Raster;
use crate::watershed_abstraction::walk_flowpaths;

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
pub fn wbt_sub_fields_abstraction(
    wd: &str,
    max_points: usize,
    clip_hillslopes: bool,
    clip_hillslope_length: f64,
    sub_field_min_area_threshold_m2: f64,
    field_raster: &str, // ag_fields/field_boundaries.tif
    output_dir: &str,   // ag_fields/sub_fields
) -> std::io::Result<()> {
    env::set_current_dir(&wd).unwrap();
    info!(
        "wbt_sub_fields_abstraction: wd={}, max_points={}, clip_hillslopes={}, clip_hillslope_length={}, sub_field_min_area_threshold_m2={}, field_raster={}, output_dir={}",
        wd,
        max_points,
        clip_hillslopes,
        clip_hillslope_length,
        sub_field_min_area_threshold_m2,
        field_raster,
        output_dir
    );
    let sub_fields_path = Path::new(output_dir);
    if sub_fields_path.exists() {
        let _ = fs::remove_dir_all(&sub_fields_path).unwrap();
    }

    let _ = std::fs::create_dir_all(Path::new(output_dir).join("slope_files"));
    let _ = std::fs::create_dir_all(Path::new(output_dir).join("slope_files/flowpaths"));

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
    let field: Raster<i32> = Raster::<i32>::read(field_raster).unwrap(); // rasterized field boundaries
    info!(
        "loaded field raster: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        field.width,
        field.height,
        field.data.len(),
        field.data.len() * size_of::<i32>(),
        (field.data.len() * size_of::<i32>()) as f64 / (1024.0 * 1024.0),
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

    let mut intersection_subwta = subwta.clone();

    // iterate over field and relief rasters
    let t0 = Instant::now();
    let mut fake_id: i32 = 0; // start sub_field_id starts at 1
    let mut fake_topaz_id_lookup: HashMap<(i32, i32), i32> = HashMap::new();
    let mut fake_topaz_areas_px: HashMap<i32, i32> = HashMap::new(); // for area thresholding

    for i in 0..subwta.data.len() {
        let field_id = field.data[i];
        let topaz_id = subwta.data[i];

        // continue if no field or no topaz id or channel pixel
        if field_id <= 0 || topaz_id <= 0 || topaz_id % 10 == 4 {
            intersection_subwta.data[i] = 0;
            continue;
        }

        let key = (field_id, topaz_id);
        if fake_topaz_id_lookup.contains_key(&key) {
            intersection_subwta.data[i] = *fake_topaz_id_lookup.get(&key).unwrap();
        } else {
            fake_id += 1;
            fake_topaz_id_lookup.insert(key, fake_id);
            intersection_subwta.data[i] = fake_id;
        }

        fake_topaz_areas_px
            .entry(intersection_subwta.data[i])
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }
    info!(
        "built intersection subwta: fake_topaz_ids={} in {:.2}s",
        fake_topaz_id_lookup.len(),
        t0.elapsed().as_secs_f64()
    );

    // remove small sub-fields based on area threshold
    let cellsize2 = subwta.cellsize * subwta.cellsize;
    let min_area_px = (sub_field_min_area_threshold_m2 / cellsize2).ceil() as i32;
    let valid_fake_ids: HashSet<i32> = fake_topaz_areas_px.iter()
        .filter(|(_k, &v)| v >= min_area_px)
        .map(|(&k, _v)| k)
        .collect();

    let t0 = Instant::now();
    for i in 0..intersection_subwta.data.len() {
        let fake_topaz_id = intersection_subwta.data[i];
        if !valid_fake_ids.contains(&fake_topaz_id) {
            intersection_subwta.data[i] = 0;
        }
    }
    info!(
        "filtered sub-fields: valid_fake_ids={} min_area_px={} in {:.2}s",
        valid_fake_ids.len(),
        min_area_px,
        t0.elapsed().as_secs_f64()
    );

    // write the intersection_subwta raster
    let t0 = Instant::now();
    let _ = intersection_subwta.write(&format!("{}/sub_field_id_map.tif", output_dir));
    info!(
        "wrote sub_field_id_map.tif in {:.2}s",
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let hillslopes = Arc::new(abstract_subfieldcatchments(
        &intersection_subwta,
        &relief,
        &flovec,
        &fvslop,
        &taspec,
    ));
    let subflow_count = hillslopes
        .subflows
        .as_ref()
        .map(|m| m.len())
        .unwrap_or(0);
    info!(
        "abstracted sub-fields: hillslopes={} subflows={} in {:.2}s",
        hillslopes.flowpaths.len(),
        subflow_count,
        t0.elapsed().as_secs_f64()
    );
    log_flowpath_collection_stats("sub-fields", &hillslopes);
    let fake_topaz_id_lookup = Arc::new(fake_topaz_id_lookup);

    let tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let out_dir = format!("{}/slope_files", output_dir);
            Box::new(move || {
                info!("writing field slps to {}", out_dir);
                let result = hillslopes.write_field_slps(
                    &out_dir,
                    max_points,
                    clip_hillslopes,
                    clip_hillslope_length,
                    lookup.as_ref(),
                );
                info!("wrote field slps to {}", out_dir);
                result
            })
        },
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let csv_path = format!("{}/fields.csv", output_dir);
            Box::new(move || {
                info!("writing fields.csv to {}", csv_path);
                let result = hillslopes.write_field_metadata_to_csv(
                    &csv_path,
                    &subwta.wgs_transform,
                    lookup.as_ref(),
                );
                info!("wrote fields.csv to {}", csv_path);
                result
            })
        },
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let csv_path = format!("{}/field_flowpaths.csv", output_dir);
            Box::new(move || {
                info!("writing field_flowpaths.csv to {}", csv_path);
                let result = hillslopes.write_field_subflows_metadata_to_csv(
                    &csv_path,
                    &subwta.wgs_transform,
                    lookup.as_ref(),
                );
                info!("wrote field_flowpaths.csv to {}", csv_path);
                result
            })
        },
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let out_dir = format!("{}/slope_files/flowpaths/", output_dir);
            Box::new(move || {
                info!("writing flowpath slps to {}", out_dir);
                let result = hillslopes.write_field_subflow_slps(&out_dir, max_points, lookup.as_ref());
                info!("wrote flowpath slps to {}", out_dir);
                result
            })
        },
    ];

    // Execute tasks in parallel
    info!("starting output tasks: {}", tasks.len());
    tasks
        .into_par_iter()
        .map(|f| f())
        .collect::<Result<Vec<_>>>()?;
    info!("completed output tasks");

    Ok(())
}

#[allow(dead_code)]
pub fn abstract_subfieldcatchments(
    intersection_subwta: &Raster<i32>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
) -> FlowpathCollection {
    let indices_map = intersection_subwta.indices_map();
    let indices_total: usize = indices_map.values().map(|v| v.len()).sum();
    info!(
        "built sub-field indices map: ids={} total_cells={} (~{:.2} MiB)",
        indices_map.len(),
        indices_total,
        (indices_total * size_of::<usize>()) as f64 / (1024.0 * 1024.0)
    );
    let fake_topaz_ids: Vec<i32> = indices_map
        .keys()
        .filter(|&&fake_topaz_id| fake_topaz_id != 0)
        .cloned()
        .collect();
    let mut hillslope_abstractions: FlowpathCollection = FlowpathCollection {
        flowpaths: Vec::new(),
        subflows: Some(HashMap::<i32, FlowpathCollection>::new()),
    };

    let results: Vec<(Flowpath, i32, FlowpathCollection)> = fake_topaz_ids
        .into_par_iter()
        .map(|fake_topaz_id| {
            let indices = indices_map
                .get(&fake_topaz_id)
                .expect("missing indices for fake_topaz_id");
            let flowpaths: FlowpathCollection = walk_flowpaths(
                fake_topaz_id,
                indices,
                &intersection_subwta,
                &relief,
                &flovec,
                &fvslop,
                &taspec,
            );
            let subcatchment: Flowpath =
                flowpaths.abstract_subfieldcatchment(&intersection_subwta, &taspec, indices);
            (subcatchment, fake_topaz_id, flowpaths)
        })
        .collect();

    for (subcatchment, fake_topaz_id, flowpaths) in results {
        hillslope_abstractions.flowpaths.push(subcatchment);
        hillslope_abstractions
            .subflows
            .as_mut()
            .unwrap()
            .insert(fake_topaz_id, flowpaths);
    }

    hillslope_abstractions
}
