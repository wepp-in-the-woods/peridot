use lazy_static::lazy_static;

use geojson::{GeoJson, Geometry, Value, Feature};
use std::fs::File;
use std::io::Write;
use std::io::Read;
use serde_json::Value as OtherValue;
use std::collections::{HashMap, HashSet, VecDeque};

use crate::raster::{Raster, precise_wgs_to_px, utm_to_px};

use crate::support::{
    compute_direction,
    circmean,
    interpolate_slp};
use crate::watershed_abstraction::{PATHS, FlowPath};


fn get_upstream_neighbors(current_index: usize, flovec: &Raster<i32>) -> Vec<usize> {
    // Return a list of neighboring cells that flow into the given cell
    let mut neighbors = Vec::new();

    for neighbor in flovec.get_neighbors(current_index) {    
        
        let flow_dir = flovec.data[neighbor];
        if flow_dir == 5 {
            continue;
        }
        let (dx, dy) = PATHS.get(&flow_dir).unwrap_or_else(|| {
            panic!("flow_dir {} not found in paths!", flow_dir)
        });
        let (x, y) = flovec.index_to_xy(neighbor);  // Convert the neighbor index to (x, y)
        let down_x: isize = x as isize + dx;
        let down_y: isize = y as isize + dy;
        let down_indx: usize = flovec.xy_to_index(down_x as usize, down_y as usize);
        
        if down_indx == current_index {
            neighbors.push(neighbor);
        }
    }
    neighbors
}

pub fn trace_catchment(
    start: usize, 
    flovec: &Raster<i32>,
    uppoundments: &Raster<i32>
) -> HashSet<usize> {
    let mut to_visit = VecDeque::new();
    let mut visited = HashSet::new();
//    visited.insert(start);

    to_visit.push_back(start);

    while let Some(indx) = to_visit.pop_front() {
        for neighbor in get_upstream_neighbors(indx, flovec) {
            if !visited.contains(&neighbor) && uppoundments.data[neighbor] == 0 {
                visited.insert(neighbor);
                to_visit.push_back(neighbor);
            }
        }
    }

    visited
}

pub fn trace_catchment_segment(
    segment_indices: HashSet<usize>, 
    flovec: &Raster<i32>,
    uppoundments: &Raster<i32>
) -> HashSet<usize> {   
    let mut all_ups = HashSet::new();

    for start in segment_indices {
        all_ups.extend(trace_catchment(start, &flovec, &uppoundments));
    }
    all_ups
}

pub fn trace_catchment_lnglat(
    lon: f64, lat: f64, 
    flovec: &Raster<i32>,
    uppoundments: &Raster<i32>
) -> HashSet<usize> {
    let (px, py) = precise_wgs_to_px(flovec.proj4.as_ref(), flovec.geo_transform.as_ref(), lon, lat);

    println!("{}, {}, {:?}", px, py, flovec.proj4.as_ref());
    let start = py as usize * flovec.width + px as usize;
    trace_catchment(start, flovec, uppoundments)
}

pub fn trace_catchment_utm(
    easting: f64, northing: f64, 
    flovec: &Raster<i32>,
    uppoundments: &Raster<i32>
) -> HashSet<usize> {

    let (px, py) = utm_to_px(&flovec.geo_transform, easting, northing);
    let start = py as usize * flovec.width + px as usize;
    trace_catchment(start, flovec, uppoundments)
}

#[allow(dead_code)]
pub fn walk_flowpath_to_indx(
    head_index: usize,
    tail_index: usize,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
    fp_id: i32
) -> FlowPath {
    let cellsize: f64 = flovec.cellsize;
    
    let mut current_index: usize = head_index;
    let mut sorted_indices: Vec<usize> = vec![head_index];
    let mut indices_hash: HashSet<usize> = HashSet::new();
    indices_hash.insert(head_index);
    let mut distances: Vec<f64> = vec![0.0];
    let mut i = 0;
    loop {
        let flow_dir: i32 = flovec.data[current_index];
        let (dx, dy) = *PATHS.get(&flow_dir).expect(&format!("flow_dir {} not found in paths!", flow_dir));

        let (x, y) = flovec.index_to_xy(current_index);
        let next_x: isize = x as isize + dx;
        let next_y: isize = y as isize + dy;
        let next_indx: usize = flovec.xy_to_index(next_x as usize, next_y as usize);

        // check if we walked in a circle
        if indices_hash.contains(&next_indx) {
            break;
        }

        sorted_indices.push(next_indx);
        indices_hash.insert(next_indx);
        distances.push(distances[distances.len() - 1] + cellsize * ((dx as f64).abs() + (dy as f64).abs()).sqrt());

        if next_indx != tail_index {
            break;
        }
        current_index = next_indx;

        i += 1;

        if i > 10000 {
            break;
        }
    }

    let n: usize = sorted_indices.len();

    let mut slopes: Vec<f64> = Vec::new();
    let mut rad_aspects: Vec<f64> = Vec::new();
    let mut elevs: Vec<f64> = Vec::new();
    for i in 0..n {
        let index: usize = sorted_indices[i];
        let slope: f64 = fvslop.data[index];
        slopes.push(slope);

        let deg_aspect = taspec.data[index];
        rad_aspects.push(deg_aspect.to_radians());

        let elev: f64 = relief.data[index];
        elevs.push(elev);
    }
    let elevation: f64 = elevs[0];

    let aspect: f64 = circmean(&rad_aspects).to_degrees();

    let (centroid_x, centroid_y) = flovec.centroid_of(&sorted_indices);

    let (x_usize, y_usize) = flovec.index_to_xy(sorted_indices[0]);
    let head = (x_usize as i32, y_usize as i32);

    let (x_usize, y_usize) = flovec.index_to_xy(sorted_indices[n - 1]);
    let tail = (x_usize as i32, y_usize as i32);

    let direction: f64 = compute_direction(head, tail);

    let width: f64 = cellsize;

    let length: f64 = distances[n - 1];  // in meters
    let total_elev: f64 = elevs[0] - elevs[n - 1]; // in meters
    let slope_scalar: f64 = total_elev / length;

    // need to normalize distances to 0-1
    let mut distances_norm: Vec<f64> = Vec::new();
    for d in &distances {
        distances_norm.push(d / length);
    }

    FlowPath::new(
        sorted_indices,
        head,
        tail,
        (centroid_x as i32, centroid_y as i32),
        distances_norm,
        slopes,
        elevs,
        -1,
        fp_id,
        length,
        width,
        aspect,
        direction,
        slope_scalar,
        cellsize,
        elevation,
        -1
    )
}


fn write_lookup_table_to_file(lookup_table: &HashMap<i32, Vec<i32>>, lookup_table_fn: &str) {
    let mut sorted_keys: Vec<i32> = lookup_table.keys().cloned().collect();
    sorted_keys.sort();

    let mut sorted_values: Vec<Vec<i32>> = Vec::new();
    for key in &sorted_keys {
        let mut values: Vec<i32> = lookup_table.get(key).unwrap().iter().cloned().collect();
        values.sort();
        sorted_values.push(values);
    }

    let mut sorted_lookup_table: HashMap<i32, HashSet<i32>> = HashMap::new();
    for (i, key) in sorted_keys.iter().enumerate() {
        sorted_lookup_table.insert(*key, sorted_values[i].iter().cloned().collect());
    }

    // Serialize the sorted lookup table to a JSON string
    let json = serde_json::to_string_pretty(&sorted_lookup_table).expect("Error serializing lookup table");

    // Write the JSON string to a file
    let mut file = File::create(lookup_table_fn).expect("Error creating file");
    file.write_all(json.as_bytes()).expect("Error writing to file");
}


#[cfg(test)]
mod tests {
    extern crate maplit;

    use super::Raster;  // Assuming Raster is in the parent module
    use std::collections::{HashSet, HashMap};
    use maplit::hashset;
    use crate::raster::utm_to_px;
    use crate::watershed_abstraction::{FlowpathCollection};
    use crate::catchment_trace::{trace_catchment_lnglat, trace_catchment_utm, write_lookup_table_to_file, walk_flowpath_to_indx};
    use geojson::{GeoJson, Geometry, Value};
    use std::fs::File;
    use std::io::Read;
    use serde_json::Value as OtherValue;
use std::process::Command;

    #[test]
    fn test_trace_catchment() {
        let path = "tests/fixtures/watershed_abstraction/litigious-sagacity/dem/topaz/FLOVEC.ARC";
        let flovec = Raster::<i32>::read(&path).unwrap();

        let lon = -115.78378912415252;
        let lat = 44.465605265553585;
        let indices = trace_catchment_lnglat(lon, lat, &flovec, &flovec.empty_clone());

        println!("indices: {:?} ({})", indices, indices.len());

        let path2 = "tests/fixtures/watershed_abstraction/litigious-sagacity/dem/topaz/BOUND.ARC";
        let bound = Raster::<i32>::read(&path2).unwrap();

        println!("bound: {:?}", bound);

        for indx in 0..bound.data.len() - 1 {
            let contains =indices.contains(&indx);

            if indices.contains(&indx) {
                assert_eq!(bound.data[indx], 1);
            }
            else {
                assert_eq!(bound.data[indx],  0);
            }
        }


    }

    #[test]
    fn test_walk_flowpath_to_indx() {
        let flovec = Raster::<i32>::read("tests/fixtures/watershed_abstraction/litigious-sagacity/dem/topaz/FLOVEC.ARC").unwrap();
        let relief = Raster::<i32>::read("tests/fixtures/watershed_abstraction/litigious-sagacity/dem/topaz/RELIEF.ARC").unwrap();


    }
    
    use std::fs;
    use serde_json::json;

    #[test]
    fn test_hubbard_brook_culvert() {

        let max_points = 1000;
        let clip_hillslopes = false;
        let clip_hillslope_length = 1000.0;

        // resolutions to test
        // 10 5, 3m

        // scale 1m DEM to other scale (start with 10m)

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

        // flovec is D8 flow direction
        let flovec = Raster::<i32>::read("/geodata/weppcloud_runs/rlew-rank-folliculitis/dem/topaz/FLOVEC.ARC").unwrap();

        // slope
        let fvslop = Raster::<f64>::read("/geodata/weppcloud_runs/rlew-rank-folliculitis/dem/topaz/FVSLOP.ARC").unwrap();

        // aspect
        let taspec = Raster::<f64>::read("/geodata/weppcloud_runs/rlew-rank-folliculitis/dem/topaz/TASPEC.ARC").unwrap();

        // DEM
        let relief = Raster::<f64>::read("/geodata/weppcloud_runs/rlew-rank-folliculitis/dem/topaz/RELIEF.ARC").unwrap();

        // Read the .geojson file in UTM
        let mut file = File::open("/geodata/locales/hubbar_brook/Hubbar_Brook_EF.Culvert_Database_HBEF_adjusted.geojson").expect("File not found");
        let mut data = String::new();
        file.read_to_string(&mut data).expect("Unable to read file");

        
        let mut culvert_upareas = flovec.empty_clone();

        // Parse the GeoJSON

        let geojson = data.parse::<GeoJson>().expect("Unable to parse GeoJSON");

        let mut i = 1;

        let mut overlap_key = 1000;
        let mut lookup_table: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut reverse_lookup: HashMap<Vec<i32>, i32> = HashMap::new();
    

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
    
                                    let indices = trace_catchment_utm(easting, northing, &flovec, &flovec.empty_clone());
                                    for &indx in &indices {
                                        let existing_val = culvert_upareas.data[indx];
                                        if existing_val != 0 {
                                            let mut vec = Vec::new();
                                            if existing_val < 300 {
                                                // First overlap at this location
                                                vec.push(existing_val);
                                            } else {
                                                // Existing overlap area, get the existing set
                                                if let Some(existing_vec) = lookup_table.get(&existing_val) {
                                                    vec = existing_vec.clone();
                                                }
                                            }
                                            vec.push(i);
                                            vec.sort();
                                            vec.dedup();
                                
                                            // Check if this set of overlaps already exists
                                            if let Some(&existing_key) = reverse_lookup.get(&vec) {
                                                // Use the existing key
                                                culvert_upareas.data[indx] = existing_key;
                                            } else {
                                                // Create a new key
                                                lookup_table.insert(overlap_key, vec.clone());
                                                reverse_lookup.insert(vec, overlap_key);
                                                culvert_upareas.data[indx] = overlap_key;
                                                overlap_key += 1;
                                            }
                                        } else {
                                            // No overlap
                                            culvert_upareas.data[indx] = i;
                                        }
                                    }

                                    let json_data = json!({
                                        "ID": id.to_string(),
                                        "enum": i,
                                        "easting": easting,
                                        "northing": northing,
                                        "n": indices.len()
                                    });

                                    println!("{}", serde_json::to_string(&json_data).unwrap());

                                    let (px, py) = utm_to_px(&flovec.geo_transform, easting, northing);
                                    let tail_index = flovec.xy_to_index(px as usize, py as usize);
                                                                        
                                    let mut flowpath_collection: FlowpathCollection = FlowpathCollection {
                                        flowpaths: Vec::new(),
                                        subflows: None
                                    };
                                
                                    for (fp_id, &head_index) in indices.iter().enumerate() {
                                        let flowpath = walk_flowpath_to_indx(
                                            head_index, tail_index, &relief,  &flovec, &fvslop, &taspec, fp_id as i32);
                                        flowpath_collection.flowpaths.push(flowpath);
                                    }

                                    let vec_indices: Vec<usize> = indices.iter().cloned().collect();

                                    if vec_indices.len() > 0 {
                                        let hillslope = flowpath_collection.abstract_hillslope(
                                            &flovec, &taspec, &vec_indices);
    
                                        hillslope.write_slp(
                                            &format!("/geodata/weppcloud_runs/rlew-rank-folliculitis/culvert_hillslope_{}.slp", i), 
                                            max_points, clip_hillslopes, clip_hillslope_length);
                                    }

    
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

        // Assuming lookup_table and lookup_table_fn are already defined
        let lookup_table_fn = "/geodata/weppcloud_runs/rlew-rank-folliculitis/culver_lookup.json";
        write_lookup_table_to_file(&lookup_table, lookup_table_fn);

        let culvert_upareas_fn = "/geodata/weppcloud_runs/rlew-rank-folliculitis/culvert_upareas_10m.i32.tif";
        culvert_upareas.write(culvert_upareas_fn);

        let culvert_upareas_opt_fn = "/geodata/weppcloud_runs/rlew-rank-folliculitis/culvert_upareas_10m.tif";

        let output = Command::new("gdal_translate")
        .arg("-ot")
        .arg("Int16")
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
        
    }

}
