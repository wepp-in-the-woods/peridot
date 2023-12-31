use lazy_static::lazy_static;

use std::io;

use geojson::{GeoJson, Geometry, Value, Feature};
use std::fs::File;
use std::io::Write;
use std::io::Read;
use serde_json::Value as OtherValue;
use std::collections::{HashMap, HashSet, VecDeque};

use std::process::Command;

use std::path::{Path, PathBuf};

use crate::raster::{Raster, precise_wgs_to_px, utm_to_px};

use crate::support::{
    compute_direction,
    circmean,
    interpolate_slp};
use crate::watershed_abstraction::{PATHS, FlowPath, FlowpathCollection};

fn remap_whitebox_d8_to_topaz(flovec: &Raster<i32>) -> Raster<i32> {
    let mut remapped_flovec = flovec.empty_clone();

    for i in 0..flovec.data.len() {
        let flow_dir = flovec.data[i];
        let new_flow_dir = match flow_dir {
            1 => 6,    // East
            2 => 3,    // Northeast
            4 => 2,    // North
            8 => 1,    // Northwest
            16 => 4,   // West
            32 => 7,   // Southwest
            64 => 8,   // South
            128 => 9,  // Southeast
            _ => 0,    // No flow or undefined
        };
        remapped_flovec.data[i] = new_flow_dir;
    }

    remapped_flovec
}

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

/*
 * dinf flovec
 *         180
 *    270       90
 *        0/360
 *     
 */

 fn unwrap_degrees(angle: f64) -> f64 {
    let mut normalized_angle = angle % 360.0;
    if normalized_angle < 0.0 {
        normalized_angle += 360.0;
    }
    normalized_angle
}

fn get_dinf_upstream_neighbors(current_index: usize, flovec: &Raster<f64>) -> Vec<usize> {
    let mut neighbors = Vec::new();
    let height = flovec.height as isize;
    let width = flovec.width as isize;

    // Convert the current index to (x, y)
    let (current_x, current_y) = flovec.index_to_xy(current_index);

    // Iterate over the 8 neighboring cells
    for dy in -1..=1 {
        for dx in -1..=1 {
            if dx == 0 && dy == 0 {
                continue; // Skip the current cell
            }

            let nx = current_x as isize + dx;
            let ny = current_y as isize + dy;

            // Check if the neighbor is within the raster bounds
            if nx >= 0 && nx < width && ny >= 0 && ny < height {
                let neighbor_index = flovec.xy_to_index(nx as usize, ny as usize);
                let flow_dir = 360.0 - flovec.data[neighbor_index];

                // Calculate the angle from the neighbor to the current cell
                let angle_to_current = unwrap_degrees((-dy as f64).atan2(dx as f64).to_degrees() + 90.0);
                // (dy-> 1, dx-> 0) = 0
                // (dy-> 0, dx-> 1) = 90
                // (dy-> -1, dx-> 0) = 180
                // (dy-> 0, dx-> -1) = 270

                // Check if the flow direction from the neighbor points to the current cell
                if is_flow_into_current_cell(flow_dir, angle_to_current) {
                    neighbors.push(neighbor_index);
                }
                
            }
        }
    }

    neighbors
}

fn is_flow_into_current_cell(flow_dir: f64, angle_to_current: f64) -> bool {
    // This threshold determines how directly the flow must be aimed at the current cell
    // to consider it an upstream neighbor.
    let threshold = 45.0; // You may adjust this threshold as needed

    let lower_bound = (angle_to_current - threshold).rem_euclid(360.0);
    let upper_bound = (angle_to_current + threshold).rem_euclid(360.0);

    if lower_bound < upper_bound {
        flow_dir >= lower_bound && flow_dir <= upper_bound
    } else {
        flow_dir >= lower_bound || flow_dir <= upper_bound
    }
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

pub fn dinf_trace_catchment(
    start: usize, 
    flovec: &Raster<f64>,
    uppoundments: &Raster<i32>
) -> (HashSet<usize>, Vec<Vec<usize>>) {
    let mut to_visit = VecDeque::new();
    let mut visited = HashSet::new();
//    visited.insert(start);

    let mut flow_links: HashMap<usize, Vec<usize>> = HashMap::new();

    to_visit.push_back(start);

    while let Some(indx) = to_visit.pop_front() {
        for neighbor in get_dinf_upstream_neighbors(indx, flovec) {
            flow_links.entry(neighbor).or_insert(Vec::new()).push(indx);
            if !visited.contains(&neighbor) && uppoundments.data[neighbor] == 0 {
                visited.insert(neighbor);
                to_visit.push_back(neighbor);
            }
        }
    }

    let mut terminal_cells: HashSet<usize> = flow_links.keys().cloned().collect();
    for values in flow_links.values() {
        for &value in values {
            terminal_cells.remove(&value);
        }
    }

    // now build list of flowpaths by walking down from the terminal cells
    // flowpaths should always converge but not branch when walking from terminal cells?

    let mut flowpaths: Vec<Vec<usize>> = Vec::new();
    for terminal_cell in terminal_cells {
        let mut flowpath: Vec<usize> = Vec::new();
        let mut current_cell = terminal_cell;
        loop {
            flowpath.push(current_cell);

            if let Some(&next_cell) = flow_links.get(&current_cell).and_then(|v| v.first()) {
                current_cell = next_cell;
            } else {
                break;
            }
        }

        if current_cell != start {
            continue;
        }

        assert!(flowpath.len() > 1);
        flowpaths.push(flowpath);
    }

    (visited, flowpaths)
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

pub fn dinf_trace_catchment_utm(
    easting: f64, northing: f64, 
    flovec: &Raster<f64>,
    uppoundments: &Raster<i32>
) -> (HashSet<usize>, Vec<Vec<usize>>) {

    let (px, py) = utm_to_px(&flovec.geo_transform, easting, northing);
    let start = py as usize * flovec.width + px as usize;
    dinf_trace_catchment(start, flovec, uppoundments)
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
    let mut i = 0;
    loop {
        let flow_dir: i32 = flovec.data[current_index];
        
        // Check if flow_dir is in PATHS and assign dx and dy
        let (dx, dy) = if let Some(&direction) = PATHS.get(&flow_dir) {
            direction
        } else {
            // If flow_dir is not found in PATHS, break from the loop
            break;
        };

        let (x, y) = flovec.index_to_xy(current_index);
        let next_x: isize = x as isize + dx;
        let next_y: isize = y as isize + dy;

        if next_x < 0 || next_x >= flovec.width as isize || next_y < 0 || next_y >= flovec.height as isize {
            break;
        }
        
        let next_indx: usize = flovec.xy_to_index(next_x as usize, next_y as usize);

        // check if we walked in a circle
        if indices_hash.contains(&next_indx) {
            break;
        }

        sorted_indices.push(next_indx);
        indices_hash.insert(next_indx);

        if next_indx != tail_index {
            break;
        }
        current_index = next_indx;

        i += 1;

        if i > 10000 {
            break;
        }
    }

    flowpath_from_indices(sorted_indices, &relief, &flovec, &fvslop, &taspec, fp_id)
}

fn flowpath_from_indices(
    indices: Vec<usize>,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
    fp_id: i32
) -> FlowPath {

    // sort indices by elevation
    let mut sorted_indices: Vec<usize> = indices.clone();
    sorted_indices.sort_by(|a, b| {
        let elev_a = relief.data[*a];
        let elev_b = relief.data[*b];
        elev_a.partial_cmp(&elev_b).unwrap()
    });

    let n: usize = sorted_indices.len();

    let cellsize = flovec.cellsize;

    let mut slopes: Vec<f64> = Vec::new();
    let mut rad_aspects: Vec<f64> = Vec::new();
    let mut elevs: Vec<f64> = Vec::new();
    let mut distances: Vec<f64> = vec![0.0];
    for i in 0..n {
        let index: usize = sorted_indices[i];
        let slope: f64 = fvslop.data[index];
        slopes.push(slope);

        let deg_aspect = taspec.data[index];
        rad_aspects.push(deg_aspect.to_radians());

        let elev: f64 = relief.data[index];
        elevs.push(elev);

        if i > 0 {
            let (x0, y0) = flovec.index_to_xy(sorted_indices[i - 1]);
            let (x1, y1) = flovec.index_to_xy(sorted_indices[i]);
            let dx: f64 = (x1 as f64 - x0 as f64) * cellsize;
            let dy: f64 = (y1 as f64 - y0 as f64) * cellsize;
            distances.push((dx.powi(2) + dy.powi(2)).sqrt());
        }
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

fn rescale_raster(src_fn: &str, dst_fn: &str, resolution: f64) -> io::Result<()> {
    let output = Command::new("gdal_translate")
        .arg("-tr")
        .arg(resolution.to_string())
        .arg(resolution.to_string())
        .arg("-r")
        .arg("cubic")
        .arg(src_fn)
        .arg(dst_fn)
        .output()?;

    if output.status.success() {
        println!("gdal_translate executed successfully.");
        Ok(())
    } else {
        eprintln!("gdal_translate failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "gdal_translate command failed"))
    }
}

fn split_path(src_fn: &str) -> (PathBuf, PathBuf) {
    let path = Path::new(src_fn);

    let parent = path.parent().unwrap_or_else(|| Path::new("")).to_path_buf();
    let file_name = path.file_name()
                        .map(Path::new)  // Convert OsStr to Path
                        .unwrap_or_else(|| Path::new(""))
                        .to_path_buf(); // Convert Path to PathBuf

    (parent, file_name)
}

fn smooth_raster(src_fn: &str, dst_fn: &str, sigma: f64) -> io::Result<()> {

    let (wd, _src_fn) = split_path(src_fn);
    let (wd2, _dst_fn) = split_path(dst_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=GaussianFilter")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _src_fn.to_string_lossy()))
        .arg(format!("-o={}", _dst_fn.to_string_lossy()))
        .arg(format!("--sigma={}", sigma.to_string()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

fn fill_depressions_raster(src_fn: &str, dst_fn: &str) -> io::Result<()> {
    let (wd, _src_fn) = split_path(src_fn);
    let (wd2, _dst_fn) = split_path(dst_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=FillDepressions")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _src_fn.to_string_lossy()))
        .arg(format!("-o={}", _dst_fn.to_string_lossy()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}


fn breach_depressions_raster(src_fn: &str, dst_fn: &str) -> io::Result<()> {
    let (wd, _src_fn) = split_path(src_fn);
    let (wd2, _dst_fn) = split_path(dst_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=BreachDepressions")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _src_fn.to_string_lossy()))
        .arg(format!("-o={}", _dst_fn.to_string_lossy()))
        .arg("fill_pits=True")
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}


fn slope_raster(dem_fn: &str, slope_fn: &str)-> io::Result<()> {
    let (wd, _dem_fn) = split_path(dem_fn);
    let (wd2, _slope_fn) = split_path(slope_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=Slope")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _dem_fn.to_string_lossy()))
        .arg(format!("-o={}", _slope_fn.to_string_lossy()))
        .output()?;

    // scale the slope to match TOPAZ by multiplying by 0.01

    let mut slope = Raster::<f64>::read(&slope_fn).unwrap();
    slope.data.iter_mut().for_each(|v| *v *= 0.01);
    slope.write(&slope_fn);

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

fn aspect_raster(dem_fn: &str, aspect_fn: &str)-> io::Result<()> {
    let (wd, _dem_fn) = split_path(dem_fn);
    let (wd2, _aspect_fn) = split_path(aspect_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=Aspect")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _dem_fn.to_string_lossy()))
        .arg(format!("-o={}", _aspect_fn.to_string_lossy()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

fn d8_flow_direction_raster(dem_fn: &str, flow_dir_fn: &str)-> io::Result<()> {
    let (wd, _dem_fn) = split_path(dem_fn);
    let (wd2, _flow_dir_fn) = split_path(flow_dir_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=D8Pointer")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _dem_fn.to_string_lossy()))
        .arg(format!("-o={}", _flow_dir_fn.to_string_lossy()))
        .output()?;

    // remap the flow direction to match TOPAZ
    let flovec = Raster::<i32>::read(&flow_dir_fn).unwrap();
    let remapped_flovec = remap_whitebox_d8_to_topaz(&flovec);
    remapped_flovec.write(&flow_dir_fn);

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

fn dinf_raster(dem_fn: &str,  dinf_fn: &str)-> io::Result<()> {
    let (wd, _dem_fn) = split_path(dem_fn);
    let (wd2, _dinf_fn) = split_path(dinf_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=DInfPointer")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _dem_fn.to_string_lossy()))
        .arg(format!("-o={}", _dinf_fn.to_string_lossy()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

fn d8_flowaccum_raster(d8_fn: &str,  d8_accum_fn: &str)-> io::Result<()> {
    let (wd2, _d8_fn) = split_path(d8_fn);
    let (wd, _d8_accum_fn) = split_path(d8_accum_fn);

    // default output type is 'cells'
    let output = Command::new("whitebox_tools")
        .arg("-r=D8FlowAccumulation")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _d8_fn.to_string_lossy()))
        .arg(format!("-o={}", _d8_accum_fn.to_string_lossy()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}


fn dinf_flowaccum_raster(dinf_fn: &str,  dinf_accum_fn: &str)-> io::Result<()> {
    let (wd2, _dinf_fn) = split_path(dinf_fn);
    let (wd, _dinf_accum_fn) = split_path(dinf_accum_fn);

    // default output type is 'sca'
    // specific catchment area (SCA), which is the upslope contributing area 
    // divided by the contour length (taken as the grid resolution)
    let output = Command::new("whitebox_tools")
        .arg("-r=DInfFlowAccumulation")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-i={}", _dinf_fn.to_string_lossy()))
        .arg(format!("-o={}", _dinf_accum_fn.to_string_lossy()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

// threshold represents the minimum area (area is used here as a surrogate for discharge) 
// required to initiate and maintain a channel. Smaller threshold values result in more 
// extensive stream networks
fn stream_raster(dinf_accum_fn: &str, stream_fn: &str, threshold: f64) -> io::Result<()> {
    let (wd, _dinf_accum_fn) = split_path(dinf_accum_fn);
    let (wd2, _stream_fn) = split_path(stream_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=ExtractStreams")
        .arg("-v")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("--flow_accum={}", _dinf_accum_fn.to_string_lossy()))
        .arg(format!("-o={}", _stream_fn.to_string_lossy()))
        .arg(format!("--threshold={}", threshold.to_string()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

fn valley_raster(dem_fn: &str, valley_fn: &str)-> io::Result<()> {
    let (wd, _dem_fn) = split_path(dem_fn);
    let (wd2, _valley_fn) = split_path(valley_fn);

    let output = Command::new("whitebox_tools")
        .arg("-r=ExtractValleys")
        .arg("-v")
        .arg("--line_thin")
        .arg("--filter=10")
        .arg(format!("--wd={}", wd.to_string_lossy()))
        .arg(format!("-dem={}", _dem_fn.to_string_lossy()))
        .arg(format!("-o={}", _valley_fn.to_string_lossy()))
        .output()?;

    if output.status.success() {
        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        //Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    extern crate maplit;

    use super::Raster;  // Assuming Raster is in the parent module
    use std::collections::{HashSet, HashMap};
    use maplit::hashset;
    use crate::raster::utm_to_px;
    use crate::watershed_abstraction::{FlowpathCollection};
    use crate::catchment_trace::{
        trace_catchment_lnglat, 
        dinf_trace_catchment_utm,
        trace_catchment_utm, 
        write_lookup_table_to_file, 
        walk_flowpath_to_indx,
        rescale_raster,
        d8_flow_direction_raster,
        d8_flowaccum_raster,
        dinf_flowaccum_raster,
        flowpath_from_indices,
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
    use std::path::Path;

    #[test]
    fn test_hubbard_brook_culvert() {

        let max_points = 1000;
        let clip_hillslopes = false;
        let clip_hillslope_length = 1000.0;

        let lidar_dem = "/geodata/locales/hubbar_brook/dem/HBEF_1m_LiDAR_DEM.tif";

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
    
        let mut culvert_upareas = d8_flovec.empty_clone();

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

                                    let (px, py) = utm_to_px(&flovec.geo_transform, easting, northing);
                                    let tail_index = flovec.xy_to_index(px as usize, py as usize);

                                    let mut _culvert_upareas = d8_flovec.empty_clone();

                                    for &indx in &indices {
                                        _culvert_upareas.data[indx] = 1;
                                    }

                                    let culvert_upareas_fn = format!("{}/culvert{}_upareas_i32.tif", wd, i); 
                                    _culvert_upareas.write(&culvert_upareas_fn);

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

                                    hillslope.write_slp(
                                        &format!("{}/culvert_hillslope_{}.slp", wd, i), 
                                        max_points, clip_hillslopes, clip_hillslope_length);

                                        let json_data = json!({
                                            "ID": id.to_string(),
                                            "enum": i,
                                            "easting": easting,
                                            "northing": northing,
                                            "n": indices.len(),
                                            "flowpath_longest_n": max_n,
                                            "terminal_flowpaths": flowpath_collection.flowpaths.len(),
                                            "hillslope_pts": hillslope.elevs.len(),
                                        });
    
                                        println!("{}", serde_json::to_string(&json_data).unwrap());

                                        
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
      
        let culvert_upareas_fn = format!("{}/culvert_upareas_i32.tif", wd); 
        culvert_upareas.write(&culvert_upareas_fn);

        let culvert_upareas_opt_fn = format!("{}/culvert_upareas.tif", wd); 

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

        // Assuming lookup_table and lookup_table_fn are already defined
        let lookup_table_fn = format!("{}/culver_lookup.json", wd);
        write_lookup_table_to_file(&lookup_table, &lookup_table_fn);

        
    }

}
