use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::io::Result;
use std::path::Path;
use std::sync::Arc;

use rayon::prelude::*;

use crate::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz;
use crate::flowpath::FlowPath;
use crate::flowpath_collection::FlowpathCollection;
use crate::raster::Raster;
use crate::watershed_abstraction::walk_flowpaths;

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
    let sub_fields_path = Path::new(output_dir);
    if sub_fields_path.exists() {
        let _ = fs::remove_dir_all(&sub_fields_path).unwrap();
    }

    let _ = std::fs::create_dir_all(Path::new(output_dir).join("slope_files"));
    let _ = std::fs::create_dir_all(Path::new(output_dir).join("slope_files/flowpaths"));

    let subwta: Raster<i32> = Raster::<i32>::read("dem/wbt/subwta.tif").unwrap(); // hillslope with topaz_ids file, channels end with 4 (e.g 24, 34, 44), subcatchments end with 1, 2, 3. It starts at 22
    let field: Raster<i32> = Raster::<i32>::read(field_raster).unwrap(); // rasterized field boundaries
    let relief: Raster<f64> = Raster::<f64>::read("dem/wbt/relief.tif").unwrap(); // dem
    let flovec_wbt: Raster<i32> = Raster::<i32>::read("dem/wbt/flovec.tif").unwrap(); // d8 flowvec
    let flovec = remap_whitebox_d8_to_topaz(&flovec_wbt);
    let fvslop: Raster<f64> = Raster::<f64>::read("dem/wbt/fvslop.tif").unwrap(); // slope
    let taspec: Raster<f64> = Raster::<f64>::read("dem/wbt/taspec.tif").unwrap(); // aspect

    let mut intersection_subwta = subwta.clone();

    // iterate over field and relief rasters
    let mut fake_id: i32 = 0;  // start sub_field_id starts at 1
    let mut fake_topaz_id_lookup: HashMap<(i32, i32), i32> = HashMap::new();
    let mut fake_topaz_areas_px: HashMap<i32, i32> = HashMap::new();  // for area thresholding

    for i in 0..subwta.data.len() {

        let field_id = field.data[i];
        let topaz_id = subwta.data[i];

        // continue if no field or no topaz id or channel pixel
        if field_id == 0 || topaz_id == 0 || topaz_id % 10 == 4 {
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

        fake_topaz_areas_px.entry(intersection_subwta.data[i])
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }

    // remove small sub-fields based on area threshold
    let cellsize2 = subwta.cellsize * subwta.cellsize;
    let min_area_px = (sub_field_min_area_threshold_m2 / cellsize2).ceil() as i32;
    let valid_fake_ids: HashSet<i32> = fake_topaz_areas_px.iter()
        .filter(|(_k, &v)| v >= min_area_px)
        .map(|(&k, _v)| k)
        .collect();

    for i in 0..intersection_subwta.data.len() {
        let fake_topaz_id = intersection_subwta.data[i];
        if fake_topaz_id == 0 {
            continue;
        }
        if !valid_fake_ids.contains(&fake_topaz_id) {
            intersection_subwta.data[i] = 0;
        }
    }

    // write the intersection_subwta raster
    let _ = intersection_subwta.write(&format!("{}/sub_field_id_map.tif", output_dir));

    let hillslopes = Arc::new(abstract_subfieldcatchments(&intersection_subwta, &relief, &flovec, &fvslop, &taspec));
    let fake_topaz_id_lookup = Arc::new(fake_topaz_id_lookup);

    let tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let out_dir = format!("{}/slope_files", output_dir);
            Box::new(move || {
                hillslopes.write_field_slps(
                    &out_dir,
                    max_points,
                    clip_hillslopes,
                    clip_hillslope_length,
                    lookup.as_ref(),
                )
            })
        },
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let csv_path = format!("{}/fields.csv", output_dir);
            Box::new(move || {
                hillslopes.write_field_metadata_to_csv(
                    &csv_path,
                    &subwta.wgs_transform,
                    lookup.as_ref(),
                )
            })
        },
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let csv_path = format!("{}/field_flowpaths.csv", output_dir);
            Box::new(move || {
                hillslopes.write_field_subflows_metadata_to_csv(
                    &csv_path,
                    &subwta.wgs_transform,
                    lookup.as_ref(),
                )
            })
        },
        {
            let lookup = Arc::clone(&fake_topaz_id_lookup);
            let hillslopes = Arc::clone(&hillslopes);
            let out_dir = format!("{}/slope_files/flowpaths/", output_dir);
            Box::new(move || {
                hillslopes.write_field_subflow_slps(
                    &out_dir,
                    max_points,
                    lookup.as_ref(),
                )
            })
        },
    ];

    // Execute tasks in parallel
    tasks.into_par_iter().map(|f| f()).collect::<Result<Vec<_>>>()?;

    Ok(())
}

#[allow(dead_code)]
pub fn abstract_subfieldcatchments(
    intersection_subwta: &Raster<i32>,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>
) -> FlowpathCollection {

    let mut unique: HashSet<i32> = intersection_subwta.unique_values();
    unique.remove(&0);

    let fake_topaz_ids: Vec<i32> = unique.into_iter().collect();
    let mut hillslope_abstractions: FlowpathCollection = FlowpathCollection {
        flowpaths: Vec::new(),
        subflows: Some(HashMap::<i32, FlowpathCollection>::new())
    };

    let results: Vec<(FlowPath, i32, FlowpathCollection)> = fake_topaz_ids.into_par_iter()
        .map(|fake_topaz_id| {
            let flowpaths: FlowpathCollection = walk_flowpaths(fake_topaz_id, &intersection_subwta, &relief, &flovec, &fvslop, &taspec);
            let subcatchment: FlowPath = flowpaths.abstract_subfieldcatchment(
                &intersection_subwta,
                &taspec);
            (subcatchment, fake_topaz_id, flowpaths)
        })
        .collect();

    for (subcatchment, fake_topaz_id, flowpaths) in results {
        hillslope_abstractions.flowpaths.push(subcatchment);
        hillslope_abstractions.subflows.as_mut().unwrap().insert(fake_topaz_id, flowpaths);
    }

    hillslope_abstractions
}
