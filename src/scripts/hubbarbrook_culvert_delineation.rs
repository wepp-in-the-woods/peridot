extern crate maplit;

use std::collections::{HashSet, HashMap};
use maplit::hashset;
use peridot::raster::{utm_to_px, Raster};
use peridot::watershed_abstraction::{FlowpathCollection};
use peridot::catchment_trace::{
    trace_catchment_lnglat, 
    dinf_trace_catchment_utm,
    trace_catchment_utm, 
    walk_flowpath_to_indx,
    flowpath_from_indices,
};

use peridot::whiteboxtools_wrappers::{
    combine_geojson_files,
    polygonize_raster,
    rescale_raster,
    d8_flow_direction_raster,
    d8_flowaccum_raster,
    dinf_flowaccum_raster,
    stream_raster,
    valley_raster,
    remap_whitebox_d8_to_topaz,
    dinf_raster,
    slope_raster,
    aspect_raster,
    smooth_raster,
    fill_depressions_raster,
    breach_depressions_raster
};

use geojson::{GeoJson, Geometry, Value};
use std::fs::File;
use std::io::Read;
use serde_json::Value as OtherValue;
use std::process::Command;

use std::fs;
use serde_json::json;
use std::path::Path;


fn main() {

    let max_points = 1000;
    let clip_hillslopes = false;
    let clip_hillslope_length = 1000.0;

    let lidar_dem = "/geodata/locales/hubbar_brook/dem/HBEF_1m_LiDAR_DEM.tif";
    let epsg = "32619";

    // resolutions to test
    // 10 5, 3m

    let resolution = 5.0;

    let wd = format!("/workdir/peridot/tests/fixtures/catchment_trace/hubbar_brook/{}", resolution as i32);

    if Path::new(&wd).exists() {
        fs::remove_dir_all(&wd).unwrap();
    }
    fs::create_dir_all(&wd).unwrap();


    let scaled_dem = format!("{}/scaled_dem.tif", wd);
    if let Err(e) = rescale_raster(&lidar_dem, &scaled_dem, resolution) {
        eprintln!("Error: {}", e);
    }


    let breached_fn = format!("{}/breached.tif", wd);
    if let Err(e) = breach_depressions_raster(&scaled_dem, &breached_fn) {
        eprintln!("Error: {}", e);
    }

    let relief_fn = format!("{}/relief.tif", wd);
    if let Err(e) = smooth_raster(&breached_fn, &relief_fn, 2.0) {
        eprintln!("Error: {}", e);
    }

//        let filled_raster = format!("{}/filled.tif", wd);
//
//        if let Err(e) = fill_depressions_raster(&relief_fn, &filled_raster) {
//            eprintln!("Error: {}", e);
//        }
//        // scale 1m DEM to 10m using whitebox tools
//
    let dinf_fn = format!("{}/dinf.tif", wd);
    if let Err(e) = dinf_raster(&relief_fn, &dinf_fn) {
        eprintln!("Error: {}", e);
    }

    let d8_fn = format!("{}/d8.tif", wd);
    if let Err(e) = d8_flow_direction_raster(&relief_fn, &d8_fn) {
        eprintln!("Error: {}", e);
    }

    /*
    let d8_accum_fn = format!("{}/d8_accum.tif", wd);
    if let Err(e) = d8_flowaccum_raster(&d8_fn, &d8_accum_fn) {
        eprintln!("Error: {}", e);
    }

    let stream_fn = format!("{}/streams.tif", wd);
    if let Err(e) = stream_raster(&d8_accum_fn, &stream_fn, 5.0) { // threshold in pixels
        eprintln!("Error: {}", e);
    }
    */

    let network_fn = format!("{}/network.tif", wd);
    if let Err(e) = valley_raster(&relief_fn, &network_fn) {
        eprintln!("Error: {}", e);
    }


    let slope_fn = format!("{}/slope.tif", wd);
    if let Err(e) = slope_raster(&relief_fn, &slope_fn) {
        eprintln!("Error: {}", e);
    }

    let aspect_fn = format!("{}/aspect.tif", wd);
    if let Err(e) = aspect_raster(&relief_fn, &aspect_fn) {
        eprintln!("Error: {}", e);
    }

    let dinf_accum_fn = format!("{}/dinf_accum.tif", wd);
    if let Err(e) = dinf_flowaccum_raster(&dinf_fn, &dinf_accum_fn) {
        eprintln!("Error: {}", e);
    }
    // as whiteboxtools command

    // relief = fill depressions in scaled DEM
    // flovec = produce d-infinity flow map (equivalent of flow vector map)
    // fvslop = produce slope
    // taspec = aspect map

    // should have all the products that would be generated from TOPAZ

    //let flovec = Raster::<i32>::read("/geodata/weppcloud_runs/offending-ebb/dem/taudem/d8_flow.tif").unwrap();
    //let relief = Raster::<f64>::read("/geodata/weppcloud_runs/offending-ebb/dem/taudem/fel.tif").unwrap();

    // rlew-rank-folliculitis is a weppcloud project created on the interface
    // this code is ran on dev.wepp.cloud

    // flovec is DINF flow direction
    let d8_flovec = Raster::<i32>::read(&d8_fn).unwrap();

    let flovec = Raster::<f64>::read(&dinf_fn).unwrap();

    // slope
    let fvslop = Raster::<f64>::read(&slope_fn).unwrap();

    // aspect
    let taspec = Raster::<f64>::read(&aspect_fn).unwrap();

    // DEM
    let relief = Raster::<f64>::read(&relief_fn).unwrap();

    // Read the .geojson file in UTM
    let mut file = File::open("/geodata/locales/hubbar_brook/Hubbar_Brook_EF.Culvert_Database_HBEF_adjusted.geojson").expect("File not found");
    let mut data = String::new();
    file.read_to_string(&mut data).expect("Unable to read file");


    // Parse the GeoJSON

    let geojson = data.parse::<GeoJson>().expect("Unable to parse GeoJSON");

    let mut i = 1;

    let mut overlap_key = 1000;
    let mut lookup_table: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut reverse_lookup: HashMap<Vec<i32>, i32> = HashMap::new();

    let mut uparea_geojson_fns: Vec<String> = Vec::new();

    match geojson {
        GeoJson::FeatureCollection(collection) => {
            for feature in collection.features {
                if let Some(properties) = feature.properties {
                    let id = properties.get("ID").and_then(|v| v.as_str()).unwrap_or_default();

                    if let Some(geometry) = feature.geometry {
                        match geometry.value {
                            Value::Point(coordinates) => {

                                // get the culvert location from the geometry, not from the attribute table
                                let easting = coordinates[0];
                                let northing = coordinates[1];

                                let (indices, flowpaths) = dinf_trace_catchment_utm(easting, northing, &flovec, &d8_flovec.empty_clone());

                                let (px, py) = utm_to_px(&flovec.geo_transform, easting, northing);
                                let tail_index = flovec.xy_to_index(px as usize, py as usize);

                                let mut _culvert_upareas = d8_flovec.empty_clone();

                                for &indx in &indices {
                                    _culvert_upareas.data[indx] = 1;
                                }

                                let culvert_upareas_fn = format!("{}/culvert{}_upareas_i32.tif", wd, i);
                                let _ = _culvert_upareas.write(&culvert_upareas_fn);

                                let culvert_upareas_opt_fn = format!("{}/culvert{}_upareas.tif", wd, i);

                                let output = Command::new("gdal_translate")
                                .arg("-ot")
                                .arg("Byte")
                                .arg("-co")
                                .arg("COMPRESS=LZW")
                                .arg(&culvert_upareas_fn)
                                .arg(&culvert_upareas_opt_fn)
                                .output()
                                .expect("Failed to execute gdal_translate");

                                if output.status.success() {
                                    println!("gdal_translate executed successfully");
                                } else {
                                    let error_message = String::from_utf8_lossy(&output.stderr);
                                    println!("Error executing gdal_translate: {}", error_message);
                                }

                                if let Err(err) = fs::remove_file(&culvert_upareas_fn) {
                                    println!("Failed to delete file {}: {}", &culvert_upareas_fn, err);
                                } else {
                                    println!("File {} deleted successfully", &culvert_upareas_fn);
                                }

                                let vec_indices: Vec<usize> = indices.iter().cloned().collect();

                                let mut flowpath_collection = FlowpathCollection  {
                                    flowpaths: Vec::new(),
                                    subflows: Some(HashMap::<i32, FlowpathCollection>::new())
                                };

                                let mut max_n = 0;
                                for flowpath_vec in flowpaths {
                                    if flowpath_vec.len() > max_n {
                                        max_n = flowpath_vec.len();
                                    }
                                    flowpath_collection.flowpaths.push(flowpath_from_indices(flowpath_vec, &relief, &d8_flovec, &fvslop, &taspec, i));
                                }

                                let hillslope = flowpath_collection.abstract_hillslope(
                                    &d8_flovec, &taspec, &vec_indices);

                                let _ = hillslope.write_slp(
                                    &format!("{}/culvert_hillslope_{}.slp", wd, i), 
                                    max_points, clip_hillslopes, clip_hillslope_length);

                                    // build shapefile with upareas and attribute table of this data
                                let culvert_elevation = relief.data[tail_index];

                                let properties = json!({
                                    "ID": id.to_string(),
                                    "enum": i,
                                    "easting": easting,
                                    "northing": northing,
                                    "n": indices.len(),
                                    "length": hillslope.length,
                                    "width": hillslope.width,
                                    "culvert_elevation": culvert_elevation,
                                    "flowpath_longest_n": max_n,
                                    "terminal_flowpaths": flowpath_collection.flowpaths.len(),
                                    "hillslope_pts": hillslope.elevs.len(),
                                });

                                // add properties to the geojson
                                let culvert_geojson_fn = format!("{}/culvert{}.geojson", wd, i);
                                let _ polygonize_raster(&culvert_upareas_opt_fn, &culvert_geojson_fn, &properties);
                                uparea_geojson_fns.push(culvert_geojson_fn);

                                println!("{}", serde_json::to_string(&properties).unwrap());

                                i += 1;
                            },
                            _ => panic!("Expected a Point geometry, got {:?}", geometry),
                        }
                    }
                }
            }
        },
        _ => println!("Expected a FeatureCollection"),
    }
    let sort_key = Some("culvert_elevation");
    let _ = combine_geojson_files(&uparea_geojson_fns, &format!("{}/culvert_upareas.geojson", wd), &epsg, sort_key);
}
// cargo run --bin hubbarbrook_culvert_delineation
