extern crate gdal;

use std::fs::File;
use std::io::prelude::*;
use std::collections::{HashSet, HashMap, VecDeque};
use std::path::Path;
use std::process::Command;

use crate::catchment_trace::trace_catchment;
use crate::raster::{Raster, wgs_to_px};
use crate::watershed_abstraction::PATHS;

fn rasterize_skid_trails(skid_geojson: &str, template_fn: &str, dst_fn: &str) ->
    i32
{

    let dataset = gdal::Dataset::open(template_fn).unwrap();
    let geo_transform = dataset.geo_transform().unwrap();

    let (width, height) = dataset.raster_size();
    let min_x = geo_transform[0];
    let max_y = geo_transform[3];
    let max_x = min_x + (geo_transform[1] * width as f64);
    let min_y = max_y + (geo_transform[5] * height as f64); // because geo_transform[5] is typically negative
    
    let cellsize = geo_transform[1];

    let wkt = dataset.projection();


    Command::new("gdal_rasterize")
        .arg("-burn").arg("1")
        .arg("-te").arg(min_x.to_string()).arg(min_y.to_string()).arg(max_x.to_string()).arg(max_y.to_string())
        .arg("-tr").arg(cellsize.to_string()).arg(cellsize.to_string()) 
        .arg("-a_srs").arg(&wkt)
        .arg(&skid_geojson)
        .arg(&dst_fn)
        .output()
        .expect("Failed to run gdal_rasterize");

    1
}

fn segment_skid_trails_raster(
    relief: &Raster<f64>,
    skid: &Raster<i32>,
    segmented_skid_raster_fn: &str) {

    let mut skid_copy = skid.clone();
    let mut segmented_skid = skid.empty_clone();
    let mut skid_segment_id = 2;

    let mut j = 0;
    while skid_copy.data.iter().any(|&x| x != 0) && j < 10_000 {

        let indices: HashSet<usize> = skid_copy.indices_of(1);

        let mut single_neighbor_indices = HashSet::new();

        for indx in &indices {
            let neighbors = skid_copy.get_neighbors(*indx);
            let mut i = 0;

            for neighbor in neighbors {
                i += skid_copy.data[neighbor];
                if i > 1 {
                    break;
                }
            }

            if i == 1 {
                single_neighbor_indices.insert(*indx);
            }
        }

        for indx in &single_neighbor_indices {
            let trail = walk_skid_trail(&relief, &skid_copy, *indx);
            for indx in trail {
                skid_copy.data[indx] = 0;
                segmented_skid.data[indx] = skid_segment_id * 100 + (j > 0) as i32;
            }
            skid_segment_id += 1;
        }

        if single_neighbor_indices.len() == 0 {
            let mut start: i32 = -1;
            let mut max_z = 0.0;
            let mut min_z: f64 = 1e38;

            for indx in &indices {
                let z = relief.data[*indx];
                if z > max_z {
                    max_z = z;
                    start = *indx as i32;
                }
            }

            if start >= 0 {
                let trail = walk_skid_trail(&relief, &skid_copy, start as usize);
                for indx in trail {
                    skid_copy.data[indx] = 0;
                    segmented_skid.data[indx] = skid_segment_id * 100 + 1;
                }
                skid_segment_id += 1;
            }
        }

        j += 1;

        println!("j {} {}", j, &skid_copy.data.iter().sum::<i32>().to_string());
    }

    segmented_skid.write(&segmented_skid_raster_fn);
}

fn is_neighbors(
    indx0: usize, 
    indx1: usize, 
    skid: &Raster<i32>,
) -> bool {
    let neighbors = skid.get_neighbors(indx0);
    return neighbors.contains(&indx1);
}

fn is_junction(
    neighbors: &Vec<usize>, 
    skid: &Raster<i32>,
) -> bool {
    if neighbors.len() == 1 {
        return false;  // Only one neighbor, so it's an end point, not a junction
    }

    for i in 0..neighbors.len() {
        for j in i+1..neighbors.len() {
            if !is_neighbors(neighbors[i], neighbors[j], &skid) {
                return true;  // The cell has neighbors that aren't neighboring each other
            }
        }
    }
    false
}

fn walk_skid_trail(
    relief: &Raster<f64>,
    skid: &Raster<i32>,
    start: usize,
) -> Vec<usize> {

    let mut elevs: Vec<f64> = Vec::new();
    elevs.push(relief.data[start]);
    
    let mut trail: Vec<usize> = Vec::new();
    trail.push(start);

    let mut stop = false;
    let mut current_indx = start;
    let mut direction = 0.0; // -1 going down, 1 going up
    while !stop {

        // find neighbors where skid == 1
        let all_neighbors = relief.get_neighbors(current_indx);
        let mut neighbors = Vec::new();
        for neighbor in all_neighbors {
            if  skid.data[neighbor] == 1 && !trail.contains(&neighbor) {
                neighbors.push(neighbor);
            }
        }
        let neighbors = neighbors;

        if neighbors.len() == 0 {
            stop = true;
        } else {
            let z = relief.data[current_indx];
            let mut k: i32 = -1;

            if direction == 0.0 {
                // pick direction of largest gradient
                let mut dz = 0.0;
                for neighbor in neighbors {
                    let _z = (z - relief.data[neighbor]).abs();
                    if _z > dz {
                        dz = _z;
                        k = neighbor as i32;

                        if relief.data[neighbor] - z > 0.0 {
                            direction = 1.0;
                        } else {
                            direction = -1.0;
                        }
                    }
                } 
            } else {
                let mut max_dz = -1e38;
                // pick the greates gradient in the cooresponding direction
                for neighbor in neighbors {
                    let dz = (relief.data[neighbor] - z) * direction as f64;
                    if dz >= 0.0 {
                        if dz > max_dz {
                            max_dz = dz;
                            k = neighbor as i32;
                        }
                    }
                }
            }

            if k >= 0 {
                current_indx = k as usize;
                let z = relief.data[current_indx];
                trail.push(current_indx);
                elevs.push(z);
            } else {
                stop = true;
            }
        }
    }

    trail
}


pub fn find_hillslope_id(
    seg_map: &HashMap<i32, HashSet<usize>>,
    indices: &HashSet<usize>,
    flovec: &Raster<i32>,
) -> i32 {
    for &index in indices.iter() {
        let neighbors = flovec.get_4d_neighbors(index);
        for (&key, seg_indices) in seg_map.iter() {
            if neighbors.iter().any(|&neighbor| seg_indices.contains(&neighbor)) {
                return key;
            }
        }
    }
    0
}

pub fn unmerge_4d_neighbors(input_set: HashSet<usize>, flovec: &Raster<i32>) -> Vec<HashSet<usize>> {
    let mut result: Vec<HashSet<usize>> = Vec::new();
    let mut visited: HashSet<usize> = HashSet::new();

    for &index in &input_set {
        if visited.contains(&index) {
            continue;
        }

        let mut current_group: HashSet<usize> = HashSet::new();
        let mut queue: VecDeque<usize> = VecDeque::new();

        queue.push_back(index);

        while let Some(current_index) = queue.pop_front() {
            if visited.insert(current_index) {
                current_group.insert(current_index);
                let neighbors = flovec.get_4d_neighbors(current_index);

                for &neighbor in &neighbors {
                    if input_set.contains(&neighbor) {
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        result.push(current_group);
    }

    result
}

pub fn merge_catchments(
    catchment_indices:  Vec<HashSet<usize>>,
    flovec: &Raster<i32>,
) -> Vec<HashSet<usize>> {
    let mut catchment_vec: Vec<HashSet<usize>> = catchment_indices.into_iter().collect();

    let mut merged = true;

    while merged {
        merged = false;

        for i in 0..catchment_vec.len() {
            for j in (i+1)..catchment_vec.len() {
                if i == j {
                    continue;
                }
                // Check if catchment_vec[i] and catchment_vec[j] have any 4D neighboring points
                if any_4d_neighbors(&catchment_vec[i], &catchment_vec[j], flovec) {
                    // Merge catchment_vec[j] into catchment_vec[i] and remove catchment_vec[j]
                    catchment_vec[i] = catchment_vec[i].union(&catchment_vec[j]).cloned().collect();
                    catchment_vec.remove(j);
                    merged = true;
                    break;
                }
            }

            if merged {
                break;
            }
        }
    }

    catchment_vec
}

// Function to check if two sets have any 4D neighboring points
fn any_4d_neighbors(set1: &HashSet<usize>, set2: &HashSet<usize>, flovec: &Raster<i32>) -> bool {
    for &index1 in set1 {
        let neighbors = flovec.get_4d_neighbors(index1);
        for &neighbor in &neighbors {
            if set2.contains(&neighbor) {
                return true;
            }
        }
    }
    false
}


pub fn catchment_skid_trails_raster(
    skid: &Raster<i32>,
    seg_skids: &Raster<i32>,
    flovec: &Raster<i32>,
    relief: &Raster<f64>,
    catchment_skid_raster_fn: &str,
) {

    // identify skid ids
    let skid_ids = seg_skids.unique_values();

    // create an empty raster to populate with skid catchment upareas
    let mut up_seg = seg_skids.empty_clone();

    // loop over the skid_ids
    for skid_id in skid_ids {
        if skid_id == 0 {
            continue;
        }

        // hillslope enumerator
        let mut hillslope_id = skid_id + 11;

        // find the indices of the skid segment and sort by reverse elevation
        let mut indices_vec: Vec<usize> = seg_skids.indices_of(skid_id).into_iter().collect();
        indices_vec.sort_by(|&a, &b| relief.data[a].partial_cmp(&relief.data[b]).unwrap().reverse());

        let mut fragments: Vec<HashSet<usize>> = Vec::new();

        // iterate over each index
        for seg_indx in indices_vec {

            // walk up to find the catchment
            let indices = trace_catchment(seg_indx, &flovec, &seg_skids);

            // if the pixel is in a channel it can have catchments from different directions so we
            // need to segregate the catchments
            for unmerged_indices in unmerge_4d_neighbors(indices, &flovec) {

                // check if the flowpath drains all the way to the segment without over-writing values in the upsegment map
                let mut has_pourpoint = false;
                for neighbor in seg_skids.get_neighbors(seg_indx) {
                    for &indx in &unmerged_indices {
                        if neighbor == indx && up_seg.data[indx] == 0 && seg_skids.data[indx] != skid_id {
                            has_pourpoint = true;
                        }
                    }
                }

                if has_pourpoint {
                    for &indx in &unmerged_indices {
                        if up_seg.data[indx] == 0 {
                            up_seg.data[indx] = hillslope_id;
                        }
                    }
                    hillslope_id += 1;
                }
                else {
                    fragments.push(unmerged_indices);
                }
            }
        }    
        
        // try to place the fragments
        for fragment in fragments {
            // figure out if the fragment drains to a hillslope_id and add pixels as this hillslope
            hillslope_id = -1;
            for indx in &fragment {
                let flow_dir: i32 = flovec.data[*indx];
                let (dx, dy) = *PATHS.get(&flow_dir).expect(&format!("flow_dir {} not found in paths!", flow_dir));
    
                let (x, y) = flovec.index_to_xy(*indx);
                let next_x: isize = x as isize + dx;
                let next_y: isize = y as isize + dy;
                let next_indx: usize = flovec.xy_to_index(next_x as usize, next_y as usize);

                if up_seg.data[next_indx] != 0 {
                    hillslope_id = up_seg.data[next_indx];
                    println!("routed fragment -> {}", hillslope_id);
                    break;
                }
            }

            if hillslope_id > 0 {
                for indx in &fragment {
                    up_seg.data[*indx] = hillslope_id;

                }
            }
        }    
    }

    // Mask skid values
    let mut inv_skid_mask = skid.clone();
    inv_skid_mask.data = inv_skid_mask.data.into_iter().map(|x| 1 - x).collect();
    let up_seg = up_seg.from_mask(&inv_skid_mask, 0);

    // Write file
    up_seg.write(catchment_skid_raster_fn);
}


pub fn catchment_skid_trails_raster_hillslopes(
    skid: &Raster<i32>,
    seg_skids: &Raster<i32>,
    flovec: &Raster<i32>,
    catchment_skid_raster_fn: &str,
) {
    let skid_ids = seg_skids.unique_values();
    let mut up_seg = seg_skids.empty_clone();

    for skid_id in skid_ids {
        if skid_id == 0 {
            continue;
        }

        // First pass: Trace the catchment and store indices
        let mut catchment_indices: Vec<HashSet<usize>> = Vec::new();
        let mut segment_indices = seg_skids.indices_of(skid_id);

        // want to trace single neighbor index first if it exists (equivalent of top slope)
        let mut endpoint_indx: usize = std::usize::MAX;
        for indx in &segment_indices {
            let neighbors = skid.get_neighbors(*indx);
            let mut i = 0;

            for neighbor in neighbors {
                i += skid.data[neighbor];
                if i > 1 {
                    break;
                }
            }

            if i == 1 {
                let indices = trace_catchment(*indx, &flovec, &seg_skids);
                catchment_indices.push(indices);
                endpoint_indx = *indx;
            }
        }

        // trace non-single neighbor indices
        for indx in segment_indices {
            if indx != endpoint_indx {
                let indices = trace_catchment(indx, &flovec, &seg_skids);
                for unmerged_indices in unmerge_4d_neighbors(indices, &flovec) {
                    catchment_indices.push(unmerged_indices);
                }
            }
        }

        // if skid is a single neighbor segment then create a separate hillslope for the top
        let mut hillslope_id = skid_id + 10;
        if endpoint_indx != std::usize::MAX {

            for indx in &catchment_indices[0] {
                up_seg.data[*indx] = hillslope_id;
            }
            catchment_indices.remove(0);
        }

        // Second pass: Merge catchments and write values to raster
       let catchment_vec = merge_catchments(catchment_indices, &flovec);
//        println!("{} {}", skid_id, catchment_vec.len());

        hillslope_id = skid_id + 11;
        for indices in &catchment_vec {
            if indices.len() < 3 {
                continue;
            }
            for &indx in indices {
                up_seg.data[indx] = hillslope_id;
            }

            hillslope_id += 1;
        }
        
    }

    // Mask skid values
    let mut inv_skid_mask = skid.clone();
    inv_skid_mask.data = inv_skid_mask.data.into_iter().map(|x| 1 - x).collect();
    let up_seg = up_seg.from_mask(&inv_skid_mask, 0);

    // Write file
    up_seg.write(catchment_skid_raster_fn);
}

pub fn identify_skid_channels(
    seg_skids: &Raster<i32>,
    flovec: &Raster<i32>,
    relief: &Raster<f64>,
    bound: &Raster<i32>,
    skid_channel_raster_fn: &str,
    segment_links_fn: &str,
) {
    let skid_ids = seg_skids.unique_values();

    // we are going to build a raster of channels
    let mut skid_chns = seg_skids.clone();

    // a structure map specifying the routing of skid segments to segments/ channels
    let mut segment_links: HashMap<i32, i32> = HashMap::new();

    let mut skid_pourpoints: Vec<(i32, usize, f64)> = Vec::new();
    
    for skid_id in skid_ids {
        if skid_id == 0 {
            continue;
        }

        let mut segment_indices = seg_skids.indices_of(skid_id);

        let mut min_z = 1e38;
        let mut pour_indx: usize = std::usize::MAX;

        for indx in &segment_indices {
            let z = relief.data[*indx];
            if z < min_z {
                pour_indx = *indx;
                min_z = z;
            }
        }

        skid_pourpoints.push((skid_id, pour_indx, min_z));
    }

    // Sort Vec by elevation
    skid_pourpoints.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

    
    // Iterate over sorted Vec
    for (skid_id, pour_indx, elevation) in skid_pourpoints {
//        println!("Skid ID: {}, Pour Index: {}, Elevation: {}", skid_id, pour_indx, elevation);
        
        let mut success = false;
        let mut routes_to: i32 = 0;


        // check if there is a junction to a neighbor skid segment
        let neighbors = relief.get_neighbors(pour_indx);
        let mut max_dz = 0.0;
        for neighbor in &neighbors {
            let neighbor_id = seg_skids.data[*neighbor];
            if neighbor_id == 0 || neighbor_id == skid_id {
                continue;
            }

            let dz = relief.data[*neighbor] - elevation;
            if dz > max_dz{
                success = true;
                max_dz = dz;
                routes_to = neighbor_id
            }

            if max_dz > 0.0 {
                println!("skid junction {} => {}", skid_id, routes_to);
                segment_links.insert(skid_id, routes_to);
            }
        }

        if !success {
            let mut current_indx = pour_indx;
            let mut fp: Vec<usize> = Vec::new();
            fp.push(current_indx);
            while !success{
                let flow_dir: i32 = flovec.data[current_indx];
                let (dx, dy) = *PATHS.get(&flow_dir).expect(&format!("flow_dir {} not found in paths!", flow_dir));
    
                let (x, y) = flovec.index_to_xy(current_indx);
                let next_x: isize = x as isize + dx;
                let next_y: isize = y as isize + dy;

                current_indx = flovec.xy_to_index(next_x as usize, next_y as usize);
                fp.push(current_indx);
                let current_id = seg_skids.data[current_indx];

                // have we reached a segment?
                if current_id != skid_id && current_id != 0 {
                    success = true;
                    routes_to = current_id
                } 

                // have we reached a channel?
                let current_id = skid_chns.data[current_indx];
                if current_id != 0 {
                    success = true;
                    routes_to = current_id
                } 

                // have we reached the edge of the watershed?
                if bound.data[current_indx] == 0 {
                    success = true;
                    routes_to = 24;
                }
            }

            if success {
                println!("skid channel to skid segment {} => {}", skid_id, routes_to);
                segment_links.insert(skid_id, routes_to);

                let channel_id = (skid_id as f64 / 100.0).round() as i32 * 100 + 4;
                for indx in fp {
                    if skid_chns.data[indx] == 0 {
                        skid_chns.data[indx] = channel_id;
                    }
                }
            }
        }
    }

    // Write file
    skid_chns.write(skid_channel_raster_fn);

    // process structure
    let mut structure: HashMap<i32, Vec<i32>> = HashMap::new();
    for (k, v) in segment_links.iter() {
        structure.entry(*v).or_insert_with(Vec::new).push(*k);
    }

    let mut file = File::create(segment_links_fn).unwrap();

    // Iterate over the HashMap and write each key-value pair to the file
    for (key, value) in structure {
        writeln!(file, "{}:{:?}", key, value);
    }
    
}

#[cfg(test)]
mod tests {
    extern crate maplit;

    use std::collections::HashSet;
    use maplit::hashset;

    use crate::raster::Raster;
    use crate::skid_trail_segmentation::{
        rasterize_skid_trails,
        segment_skid_trails_raster,
        catchment_skid_trails_raster,
        identify_skid_channels
    };


    #[test]
    fn test_rasterize() {
        let geojson_fn = "/geodata/salvage_logging/north_star/Skid_segments.utm.geojson";
        let template_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/dem/dem.tif";
        let dst_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/Skid.tif";

        rasterize_skid_trails(geojson_fn, template_fn, dst_fn);

        let skid = Raster::<i32>::read(&dst_fn).unwrap();

        let relief_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/dem/topaz/RELIEF.ARC";
        let relief = Raster::<f64>::read(&relief_fn).unwrap();

        let bound_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/dem/topaz/BOUND.ARC";
        let bound = Raster::<i32>::read(&bound_fn).unwrap();

        let bounded_skid = skid.from_mask(&bound, 0);

        let dst_fn2 = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/Segmented_Skid.tif";
        segment_skid_trails_raster(&relief, &bounded_skid, &dst_fn2);
        let seg_skids = Raster::<i32>::read(&dst_fn2).unwrap();

        let flovec_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/dem/topaz/FLOVEC.ARC";
        let flovec = Raster::<i32>::read(&flovec_fn).unwrap();
        let dst_fn3 = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/Skid_Catchments.tif";
        catchment_skid_trails_raster(&skid, &seg_skids, &flovec, &relief, &dst_fn3);

        let dst_fn4 = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/Skid_Channels.tif";
        let seg_links_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/segment_links.txt";
        identify_skid_channels(&seg_skids, &flovec, &relief, &bound, &dst_fn4, &seg_links_fn);
        
    }
}

