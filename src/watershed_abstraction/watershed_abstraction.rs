use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::io::Result;
use std::path::Path;

use lazy_static::lazy_static;
use rayon::prelude::*;

use crate::flowpath::Flowpath;
use crate::flowpath_collection::FlowpathCollection;
use crate::netw::{read_netw_tab, write_network, ChannelNode};
use crate::raster::Raster;
use crate::support::{circmean, compute_direction};
use crate::watershed_abstraction::watershed_manifest::{
    channels_summary, flowpaths_summary, hillslopes_summary, write_watershed_readme,
    ManifestRunFlags, TabularOutputSummary,
};
use crate::wbt_netw::Link;

pub trait LinkAttributes {
    fn get_areaup(&self) -> f64;
    fn get_order(&self) -> i32;
}

impl LinkAttributes for ChannelNode {
    fn get_areaup(&self) -> f64 {
        self.areaup
    }
    fn get_order(&self) -> i32 {
        self.order
    }
}

impl LinkAttributes for Link {
    fn get_areaup(&self) -> f64 {
        self.areaup
    }
    fn get_order(&self) -> i32 {
        self.order
    }
}

// 1     2     3
//   \   |   /
//    \  | /
// 4 <-- o --> 6
//     / | \
//   /   |   \
// 7     8     9

lazy_static! {
    pub static ref PATHS: HashMap<i32, (isize, isize)> = {
        let mut m = HashMap::new();
        m.insert(1, (-1, -1));
        m.insert(2, (0, -1));
        m.insert(3, (1, -1));
        m.insert(4, (-1, 0));
        m.insert(6, (1, 0));
        m.insert(7, (-1, 1));
        m.insert(8, (0, 1));
        m.insert(9, (1, 1));
        m
    };
}

pub const MIN_CHANNEL_SLOPE_SCALAR: f64 = 1e-4;

fn compute_channel_slope_scalar(total_elev: f64, length: f64, topaz_id: i32) -> f64 {
    if total_elev <= 0.0 {
        log::warn!(
            "total_elev is {:.6} for channel topaz_id: {}, assigning minimum slope",
            total_elev,
            topaz_id
        );
        MIN_CHANNEL_SLOPE_SCALAR
    } else {
        total_elev / length
    }
}

#[allow(dead_code)]
pub fn abstract_watershed(
    wd: &str,
    max_points: usize,
    clip_hillslopes: bool,
    clip_hillslope_length: f64,
    bieger2015_widths: bool,
    write_flowpaths: bool,
) -> std::io::Result<()> {
    env::set_current_dir(&wd).unwrap();

    let watershed_path = Path::new("watershed");

    if watershed_path.exists() {
        let _ = fs::remove_dir_all(&watershed_path).unwrap();
    }

    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/hillslopes"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/flowpaths"));

    let subwta: Raster<i32> = Raster::<i32>::read("dem/topaz/SUBWTA.ARC").unwrap(); // hillslope with topaz_ids file, channels end with 4 (e.g 24, 34, 44), subcatchments end with 1, 2, 3. It starts at 22
    let relief: Raster<f32> = Raster::<f32>::read("dem/topaz/RELIEF.ARC").unwrap(); // dem
    let flovec: Raster<u8> = Raster::<u8>::read("dem/topaz/FLOVEC.ARC").unwrap(); // d8 flowvec
    let fvslop: Raster<f32> = Raster::<f32>::read("dem/topaz/FVSLOP.ARC").unwrap(); // slope
    let taspec: Raster<f32> = Raster::<f32>::read("dem/topaz/TASPEC.ARC").unwrap(); // aspect

    let subwta_indices = subwta.indices_map();

    // this is a ASCII tabular report with channel node connection information
    let (netw, network) = read_netw_tab("dem/topaz/NETW.TAB", &subwta).unwrap();
    let _ = write_network("watershed/network.txt", &network);

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

    let mut tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        Box::new(|| channels.write_channel_slp("watershed/slope_files/channels.slp", max_points)),
        Box::new(|| {
            channels
                .write_chn_metadata_to_parquet("watershed/channels.parquet", &subwta.wgs_transform)
        }),
        Box::new(|| {
            hillslopes.write_slps(
                "watershed/slope_files/hillslopes/",
                max_points,
                clip_hillslopes,
                clip_hillslope_length,
            )
        }),
        Box::new(|| {
            hillslopes
                .write_metadata_to_parquet("watershed/hillslopes.parquet", &subwta.wgs_transform)
        }),
        Box::new(|| channels.write_geojson(&subwta, "watershed/channels.geojson")),
    ];

    if write_flowpaths {
        tasks.push(Box::new(|| {
            hillslopes.write_subflow_slps("watershed/slope_files/flowpaths/", max_points)
        }));
    }

    // Execute tasks in parallel
    tasks
        .into_par_iter()
        .map(|f| f())
        .collect::<Result<Vec<_>>>()?;

    // Flowpaths parquet export is memory-heavy; run it after the parallel pool to
    // reduce peak memory contention with other writers.
    if write_flowpaths {
        hillslopes.write_subflows_metadata_to_parquet(
            "watershed/flowpaths.parquet",
            &subwta.wgs_transform,
        )?;
    }

    let mut tabular_outputs: Vec<TabularOutputSummary> = vec![
        channels_summary(channels.flowpaths.len(), "parquet"),
        hillslopes_summary(hillslopes.flowpaths.len(), "parquet"),
    ];
    if write_flowpaths {
        let flowpath_rows = hillslopes.subflow_row_count();
        tabular_outputs.push(flowpaths_summary(flowpath_rows, "parquet"));
    }

    let run_flags = ManifestRunFlags {
        command: "abstract_watershed",
        max_points,
        clip_hillslopes,
        clip_hillslope_length,
        bieger2015_widths,
        skip_flowpaths: !write_flowpaths,
        representative_flowpath: false,
    };
    write_watershed_readme(Path::new("watershed"), &run_flags, &tabular_outputs)?;

    Ok(())
}

#[allow(dead_code)]
pub fn abstract_subcatchments(
    subwta: &Raster<i32>,
    subwta_indices: &HashMap<i32, Vec<usize>>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
    channels: &FlowpathCollection,
    write_flowpaths: bool,
) -> FlowpathCollection {
    let topaz_ids: Vec<i32> = subwta_indices
        .keys()
        .filter(|&&topaz_id| topaz_id % 10 != 4)
        .cloned()
        .collect();
    let mut hillslope_abstractions: FlowpathCollection = FlowpathCollection {
        flowpaths: Vec::new(),
        subflows: if write_flowpaths {
            Some(HashMap::<i32, FlowpathCollection>::new())
        } else {
            None
        },
    };

    if write_flowpaths {
        let results: Vec<(Flowpath, i32, FlowpathCollection)> = topaz_ids
            .into_par_iter()
            .map(|topaz_id| {
                let indices = subwta_indices
                    .get(&topaz_id)
                    .expect("missing indices for topaz_id");
                let flowpaths: FlowpathCollection = walk_flowpaths(
                    topaz_id, indices, &subwta, &relief, &flovec, &fvslop, &taspec,
                );
                let subcatchment: Flowpath = flowpaths
                    .abstract_subcatchment(&subwta, &taspec, &fvslop, &flovec, &channels, indices);
                (subcatchment, topaz_id, flowpaths)
            })
            .collect();

        for (subcatchment, topaz_id, flowpaths) in results {
            hillslope_abstractions.flowpaths.push(subcatchment);
            hillslope_abstractions
                .subflows
                .as_mut()
                .unwrap()
                .insert(topaz_id, flowpaths);
        }
    } else {
        let results: Vec<(Flowpath, i32)> = topaz_ids
            .into_par_iter()
            .map(|topaz_id| {
                let indices = subwta_indices
                    .get(&topaz_id)
                    .expect("missing indices for topaz_id");
                let flowpaths: FlowpathCollection = walk_flowpaths(
                    topaz_id, indices, &subwta, &relief, &flovec, &fvslop, &taspec,
                );
                let subcatchment: Flowpath = flowpaths
                    .abstract_subcatchment(&subwta, &taspec, &fvslop, &flovec, &channels, indices);
                (subcatchment, topaz_id)
            })
            .collect();

        for (subcatchment, _topaz_id) in results {
            hillslope_abstractions.flowpaths.push(subcatchment);
        }
    }

    hillslope_abstractions
}

#[allow(dead_code)]
pub fn walk_channels<T: LinkAttributes>(
    subwta: &Raster<i32>,
    subwta_indices: &HashMap<i32, Vec<usize>>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
    netw: &HashMap<i32, T>,
    bieger2015_widths: bool,
) -> FlowpathCollection
where
    T: LinkAttributes + Sync,
{
    let mut topaz_ids: Vec<i32> = subwta_indices
        .keys()
        .filter(|&&topaz_id| topaz_id % 10 == 4)
        .cloned()
        .collect();
    topaz_ids.sort();

    let flowpaths: Vec<Flowpath> = topaz_ids
        .into_par_iter()
        .map(|topaz_id| {
            walk_channel(
                topaz_id,
                subwta_indices
                    .get(&topaz_id)
                    .expect("missing indices for topaz_id"),
                &subwta,
                &relief,
                &flovec,
                &fvslop,
                &taspec,
                netw[&topaz_id].get_areaup(),
                netw[&topaz_id].get_order(),
                bieger2015_widths,
            )
        })
        .collect();

    FlowpathCollection {
        flowpaths: flowpaths,
        subflows: None,
    }
}

#[allow(dead_code)]
pub fn walk_channel(
    topaz_id: i32,
    indices: &Vec<usize>,
    subwta: &Raster<i32>,
    relief: &Raster<f32>,
    _flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
    areaup: f64,
    order: i32,
    bieger2015_widths: bool,
) -> Flowpath {
    // get hashset of indices of topaz_id
    assert!(
        !indices.is_empty(),
        "indices.len() == 0 for topaz_id: {}",
        topaz_id
    );

    let cellsize = subwta.cellsize;

    // build hashmap of indices and elevations
    let mut elevs_d: HashMap<usize, f64> = HashMap::new();
    for i in indices {
        elevs_d.insert(*i, relief.data[*i] as f64);
    }

    // build Vec<i32> of indices in descending order of elevation
    let mut indices_elevs: Vec<_> = elevs_d.iter().collect();
    indices_elevs.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));

    let sorted_indices: Vec<usize> = indices_elevs
        .clone()
        .into_iter()
        .map(|(index, _)| *index)
        .collect();
    let elevs: Vec<f64> = indices_elevs
        .clone()
        .into_iter()
        .map(|(_, elev)| *elev)
        .collect();
    let elevation: f64 = elevs[0];

    let mut n: usize = sorted_indices.len();

    // calculate distances between the indices
    let mut distances: Vec<f64> = vec![0.0];
    for i in 1..n {
        let index1 = sorted_indices[i - 1];
        let index2 = sorted_indices[i];
        let distance = subwta.distance_between(index1, index2);
        assert!(distance > 0.0, "distance is 0.0 for topaz_id: {}", topaz_id);
        distances.push(distances[i - 1] + distance);
    }

    let mut slopes: Vec<f64> = Vec::new();
    let mut rad_aspects: Vec<f64> = Vec::new();
    for i in 0..n {
        let index = sorted_indices[i];
        let slope = fvslop.data[index] as f64;
        slopes.push(slope);

        let deg_aspect = taspec.data[index] as f64;
        rad_aspects.push(deg_aspect.to_radians());
    }

    let mut aspect: f64 = circmean(rad_aspects.as_slice()).to_degrees();
    if aspect < 0.0 {
        aspect += 360.0; // ensure aspect is in [0, 360) range
    }

    let (centroid_x, centroid_y) = subwta.centroid_of(indices);

    let (x_usize, y_usize) = subwta.index_to_xy(sorted_indices[0]);
    let head = (x_usize as i32, y_usize as i32);

    let (x_usize, y_usize) = subwta.index_to_xy(sorted_indices[n - 1]);
    let tail: (i32, i32) = (x_usize as i32, y_usize as i32);

    let direction: f64 = compute_direction(head, tail);

    let slope_scalar: f64;
    if n == 1 {
        distances.push(cellsize);
        slopes.push(slopes[0]);
        slope_scalar = slopes[0];
        n += 1;
    } else {
        let total_elev: f64 = elevs[0] - elevs[n - 1]; // in meters
        slope_scalar = compute_channel_slope_scalar(total_elev, distances[n - 1], topaz_id);
    }

    let length: f64 = distances[n - 1];
    assert!(
        length > 0.0,
        "length is 0.0 for topaz_id: {}, {:?} {:?}",
        topaz_id,
        distances,
        slopes
    );

    // need to normalize distances to 0-1
    let mut distances_norm: Vec<f64> = Vec::new();
    for d in &distances {
        distances_norm.push(d / length);
    }

    //    let areaup: f64 = netw.get(&topaz_id).map(|node| node.areaup as f64).unwrap() * cellsize * cellsize;  // area in m^2
    //    let order = netw.get(&topaz_id).map_or(8, |node| std::cmp::min(node.order, 8));

    //     1    2    3    4    5    6    7    8
    let mut width: f64 = [-1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0][order as usize];

    if bieger2015_widths {
        let da: f64 = areaup * 1e-6; // area in km^2
                                     // width = 2.70 * da.powf(0.352); // USA model https://github.com/rogerlew/wepppy/issues/268
        width = 1.24 * da.powf(0.435); // Rocky Mountain System

        if width < 0.324017 {
            width = 0.324017;
        }
    }

    Flowpath::new(
        sorted_indices,
        head,
        tail,
        (centroid_x as i32, centroid_y as i32),
        distances_norm,
        slopes,
        elevs,
        topaz_id,
        -1,
        length,
        width,
        aspect,
        direction,
        slope_scalar,
        cellsize,
        elevation,
        order,
        areaup,
    )
}

#[allow(dead_code)]
pub fn walk_flowpaths(
    topaz_id: i32,
    indices: &Vec<usize>,
    subwta: &Raster<i32>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
) -> FlowpathCollection {
    let flowpaths: Vec<Flowpath> = indices
        .iter()
        .enumerate()
        .map(|(i, &head_index)| {
            walk_flowpath(
                topaz_id,
                head_index,
                &subwta,
                &relief,
                &flovec,
                &fvslop,
                &taspec,
                i.try_into().unwrap(),
            )
        })
        .collect();

    FlowpathCollection {
        flowpaths: flowpaths,
        subflows: None,
    }
}

#[allow(dead_code)]
pub fn walk_flowpath(
    topaz_id: i32,
    head_index: usize,
    subwta: &Raster<i32>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
    fp_id: i32,
) -> Flowpath {
    assert_eq!(subwta.data[head_index], topaz_id);

    let cellsize: f64 = subwta.cellsize;

    let mut current_index: usize = head_index;
    let mut sorted_indices: Vec<usize> = vec![head_index];
    let mut indices_hash: HashSet<usize> = HashSet::new();
    indices_hash.insert(head_index);
    let mut distances: Vec<f64> = vec![0.0];
    let mut i = 0;
    loop {
        let flow_dir: i32 = flovec.data[current_index] as i32;

        if flow_dir == 0 {
            println!(
                "Warning: flow_dir is 0 at index {} for topaz_id {}",
                current_index, topaz_id
            );
            break;
        }

        let (dx, dy) = *PATHS
            .get(&flow_dir)
            .expect(&format!("flow_dir {} not found in paths!", flow_dir));

        let (x, y) = subwta.index_to_xy(current_index);
        let next_x: isize = x as isize + dx;
        let next_y: isize = y as isize + dy;
        let next_indx: usize = subwta.xy_to_index(next_x as usize, next_y as usize);

        // check if we walked in a circle
        if indices_hash.contains(&next_indx) {
            break;
        }

        sorted_indices.push(next_indx);
        indices_hash.insert(next_indx);
        distances.push(
            distances[distances.len() - 1]
                + cellsize * ((dx as f64).abs() + (dy as f64).abs()).sqrt(),
        );

        let next_topaz_id: i32 = subwta.data[next_indx];
        if next_topaz_id != topaz_id {
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
        let slope: f64 = fvslop.data[index] as f64;
        slopes.push(slope);

        let deg_aspect = taspec.data[index] as f64;
        rad_aspects.push(deg_aspect.to_radians());

        let elev: f64 = relief.data[index] as f64;
        elevs.push(elev);
    }
    let elevation: f64 = elevs[0];

    let aspect: f64 = circmean(&rad_aspects).to_degrees();

    let (centroid_x, centroid_y) = subwta.centroid_of(&sorted_indices);

    let (x_usize, y_usize) = subwta.index_to_xy(sorted_indices[0]);
    let head = (x_usize as i32, y_usize as i32);

    let (x_usize, y_usize) = subwta.index_to_xy(sorted_indices[n - 1]);
    let tail = (x_usize as i32, y_usize as i32);

    let direction: f64 = compute_direction(head, tail);

    let width: f64 = cellsize;

    let length: f64 = distances[n - 1]; // in meters
    let total_elev: f64 = elevs[0] - elevs[n - 1]; // in meters
    let slope_scalar: f64 = total_elev / length;

    // need to normalize distances to 0-1
    let mut distances_norm: Vec<f64> = Vec::new();
    for d in &distances {
        distances_norm.push(d / length);
    }

    Flowpath::new(
        sorted_indices,
        head,
        tail,
        (centroid_x as i32, centroid_y as i32),
        distances_norm,
        slopes,
        elevs,
        topaz_id,
        fp_id,
        length,
        width,
        aspect,
        direction,
        slope_scalar,
        cellsize,
        elevation,
        -1,
        -1.0,
    )
}

#[cfg(test)]
mod tests {
    use super::{compute_channel_slope_scalar, MIN_CHANNEL_SLOPE_SCALAR};

    #[test]
    fn compute_channel_slope_scalar_returns_ratio_for_positive_drop() {
        let slope = compute_channel_slope_scalar(5.0, 20.0, 1234);
        assert!((slope - 0.25).abs() < 1e-12);
    }

    #[test]
    fn compute_channel_slope_scalar_returns_minimum_for_zero_drop() {
        let slope = compute_channel_slope_scalar(0.0, 20.0, 1234);
        assert!((slope - MIN_CHANNEL_SLOPE_SCALAR).abs() < 1e-12);
    }

    #[test]
    fn compute_channel_slope_scalar_returns_minimum_for_negative_drop() {
        let slope = compute_channel_slope_scalar(-1.0, 20.0, 1234);
        assert!((slope - MIN_CHANNEL_SLOPE_SCALAR).abs() < 1e-12);
    }
}
