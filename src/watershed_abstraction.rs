
use crate::raster::{Raster, px_to_wgs};
use crate::support::{
    compute_direction,
    circmean,
    interpolate_slp};
use crate::netw::{read_netw_tab, ChannelNode, write_network};
use crate::wbt_netw::{read_wbt_netw_tab, Link};
use crate::douglas_peucker::douglas_peucker;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;
use std::fs;
use std::io::{Write, Result};
use std::fs::OpenOptions;

use std::env;

use rayon::prelude::*;

use geojson::{Feature, Geometry, FeatureCollection};
use geojson::Value::LineString;

use interp::interp;

extern crate serde_json;
use serde_json::to_string_pretty;

use maplit::hashmap;

use lazy_static::lazy_static;

trait LinkAttributes {
    fn get_areaup(&self) -> f64;
    fn get_order(&self) -> i32;
}

impl LinkAttributes for ChannelNode {
    fn get_areaup(&self) -> f64 { self.areaup }
    fn get_order(&self) -> i32 {self.order }
}

impl LinkAttributes for Link {
    fn get_areaup(&self) -> f64 { self.areaup }
    fn get_order(&self) -> i32 {self.order }
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


#[allow(dead_code)]
pub fn abstract_watershed(
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

    let subwta: Raster<i32> = Raster::<i32>::read("dem/topaz/SUBWTA.ARC").unwrap(); // hillslope with topaz_ids file, channels end with 4 (e.g 24, 34, 44), subcatchments end with 1, 2, 3. It starts at 22
    let relief: Raster<f64> = Raster::<f64>::read("dem/topaz/RELIEF.ARC").unwrap(); // dem
    let flovec: Raster<i32> = Raster::<i32>::read("dem/topaz/FLOVEC.ARC").unwrap(); // d8 flowvec
    let fvslop: Raster<f64> = Raster::<f64>::read("dem/topaz/FVSLOP.ARC").unwrap(); // slope
    let taspec: Raster<f64> = Raster::<f64>::read("dem/topaz/TASPEC.ARC").unwrap(); // aspect

    // this is a ASCII tabular report with channel node connection informaiton
    let (netw, network) = read_netw_tab("dem/topaz/NETW.TAB", &subwta).unwrap();
    let _ = write_network("watershed/network.txt", &network);

    let channels: FlowpathCollection = walk_channels(&subwta, &relief, &flovec, &fvslop, &taspec, &netw, bieger2015_widths);
    let hillslopes: FlowpathCollection = abstract_subcatchments(&subwta, &relief, &flovec, &fvslop, &taspec, &channels);

    let tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        Box::new(|| channels.write_channel_slp("watershed/slope_files/channels.slp", max_points)),
        Box::new(|| channels.write_chn_metadata_to_csv("watershed/channels.csv", &subwta.wgs_transform)),
        Box::new(|| hillslopes.write_slps("watershed/slope_files/hillslopes/", max_points, clip_hillslopes, clip_hillslope_length)),
        Box::new(|| hillslopes.write_metadata_to_csv("watershed/hillslopes.csv", &subwta.wgs_transform)),
        Box::new(|| hillslopes.write_subflows_metadata_to_csv("watershed/flowpaths.csv", &subwta.wgs_transform)),
        Box::new(|| hillslopes.write_subflow_slps("watershed/slope_files/flowpaths/", max_points, clip_hillslopes, clip_hillslope_length)),
        Box::new(|| channels.write_geojson(&subwta, "watershed/channels.geojson")),
    ];

    // Execute tasks in parallel
    tasks.into_par_iter().map(|f| f()).collect::<Result<Vec<_>>>()?;

    Ok(())

}

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
    let flovec: Raster<i32> = Raster::<i32>::read("dem/wbt/flovec.tif").unwrap(); // d8 flowvec
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
        Box::new(|| hillslopes.write_subflow_slps("watershed/slope_files/flowpaths/", max_points, clip_hillslopes, clip_hillslope_length)),
        Box::new(|| channels.write_geojson(&subwta, "watershed/channels.geojson")),
    ];

    // Execute tasks in parallel
    tasks.into_par_iter().map(|f| f()).collect::<Result<Vec<_>>>()?;

    Ok(())

}

#[allow(dead_code)]
pub fn abstract_subcatchments(
    subwta: &Raster<i32>,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
    channels: &FlowpathCollection
) -> FlowpathCollection {

    let unique: HashSet<i32> = subwta.unique_values();
    let topaz_ids: Vec<i32> = unique.into_iter().filter(|&topaz_id| topaz_id % 10 != 4).collect();
    let mut hillslope_abstractions: FlowpathCollection = FlowpathCollection {
        flowpaths: Vec::new(),
        subflows: Some(HashMap::<i32, FlowpathCollection>::new())
    };

    let results: Vec<(FlowPath, i32, FlowpathCollection)> = topaz_ids.into_par_iter()
        .map(|topaz_id| {
            let flowpaths: FlowpathCollection = walk_flowpaths(topaz_id, &subwta, &relief, &flovec, &fvslop, &taspec);
            let subcatchment: FlowPath = flowpaths.abstract_subcatchment(
                &subwta,
                &taspec,
                &channels);
            (subcatchment, topaz_id, flowpaths)
        })
        .collect();

    for (subcatchment, topaz_id, flowpaths) in results {
        hillslope_abstractions.flowpaths.push(subcatchment);
        hillslope_abstractions.subflows.as_mut().unwrap().insert(topaz_id, flowpaths);
    }

    hillslope_abstractions
}

#[derive(Debug, Clone)]
pub struct FlowpathCollection {
    pub flowpaths: Vec<FlowPath>,
    pub subflows: Option<HashMap<i32, FlowpathCollection>>
}

impl FlowpathCollection {

    #[allow(dead_code)]
    pub fn get_fp_by_topaz_id(&self, topaz_id: i32) -> Option<&FlowPath> {
        self.flowpaths.iter().find(|fp| fp.topaz_id == topaz_id)
    }

    #[allow(dead_code)]
    pub fn get_longest_fp(&self) -> &FlowPath {
        let mut max_length: f64 = 0.0;
        let mut max_index: usize = 0;
        for (i, fp) in self.flowpaths.iter().enumerate() {
            if fp.length > max_length {
                max_length = fp.length;
                max_index = i;
            }
        }
        &self.flowpaths[max_index]
    }

    /// calculates the length of a subcatchment from the flowpaths
    /// contained within the subcatchment. The length is weighted_flowpaths
    /// by the flowpaths length relative to its area.
    ///
    /// distances should be an array of distances between cells along the
    /// flowpath (not cumulative distance)
    ///
    /// eq. 3.4 in Thomas Cochrane's Dissertation
    #[allow(dead_code)]
    pub fn garbrecht_length(&self) -> f64 {

        let mut sum_xa: f64 = 0.0;
        let mut sum_a: f64 = 0.0;
        let mut n: f64 = 0.0;

        for fp in &self.flowpaths {
            let length = fp.length;
            let area = fp.area_m2();

            sum_xa += length * area;
            sum_a += area;
            n += 1.0;
        }

        sum_xa / (sum_a * n)
    }

    ///calculates weighted slopes based on the flowpaths contained on the hillslope
    #[allow(dead_code)]
    pub fn weighted_slope_average_from_fps(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>)  {

        let longest_fp: &FlowPath = self.get_longest_fp();
        let mut num_points: usize = longest_fp.distances_norm.len();
        let longest_length: f64 = longest_fp.length;

        if num_points == 1 {
            let slope: f64 = longest_fp.slopes[0];
            let eps: Vec<f64> = vec![slope, slope];
            let distances_norm: Vec<f64> = vec![0.0, 1.0];
            let distance_p: Vec<f64> = vec![0.0, longest_length];
            return (eps, distance_p, distances_norm);
        }

        // for each flowpath determine the distance from channel
        // this requires reversing the elements in the distance
        // array and calculating the cumulative sum

        num_points += 1;
        let mut distance_p: Vec<f64> = Vec::new();
        let step = 1.0 / (num_points - 1) as f64;
        for i in 0..(num_points) {
            let distance = longest_length * step * i as f64;
            distance_p.push(distance);
        }
        assert_eq!(distance_p.len(), num_points);
        assert_eq!(distance_p[0], 0.0);
        assert!((distance_p[num_points - 1] - longest_length).abs() < 1e-3);

        let mut eps: Vec<f64> = Vec::new();

        // we will weight the slope at each distance away from the channel
        for d_p in &distance_p {
            let mut num: f64 = 0.0;   // to hold numerator value
            let mut kpsum: f64 = 0.0; // to hold k_p sum

            for fp in &self.flowpaths {
                if d_p > &(fp.length + 1e-6)  {
                    continue;
                }

                let mut slp_p: f64 = fp.interp_slp_at_distance_to_channel(*d_p);
                if slp_p <= 0.0 {
                    slp_p = 0.001;
                }
                assert!(slp_p.is_finite());
                assert!(slp_p >= 0.0);

                num += slp_p * fp.kp;
                kpsum += fp.kp;
            }

            assert!(kpsum > 0.0, "kpsum is 0.0, d_p: {}, num_points: {}, longest_length: {}", d_p, num_points, longest_length);
            let weighted_slp: f64 = num / kpsum;


            // store the weighted slope estimate
            eps.push(weighted_slp);
        }

        // normalize distance_p array
        let distances_norm: Vec<f64> = distance_p.iter().map(|&x| x / longest_length).collect();

        // reverse weighted slopes
        eps.reverse();

        (eps, distance_p, distances_norm)
    }

    #[allow(dead_code)]
    pub fn abstract_subcatchment(&self,
        subwta: &Raster<i32>,
        taspec: &Raster<f64>,
        channels: &FlowpathCollection) -> FlowPath {

        let cellsize: f64 = subwta.cellsize;
        let cellsize2: f64 = cellsize * cellsize;
        let topaz_id: i32 = self.flowpaths[0].topaz_id;

        // get indices of subcatchment
        let indices: HashSet<usize> = subwta.indices_of(topaz_id);
        let area = indices.len() as f64 * cellsize2;

        // find corresponding chn_id
        let chn_id: i32 = ((topaz_id as f64 / 10.0).floor() * 10.0) as i32 + 4;
        let chn_summary: &FlowPath = &channels.get_fp_by_topaz_id(chn_id).unwrap();

        let longest_fp = self.get_longest_fp();

        // If subcatchment is a source type then we calculate the distance
        // by taking a weighted average based on the length of the flowpaths
        // contained in the subcatchment
        let length: f64;
        let width: f64;
        if topaz_id % 10 == 1 {
            // find length by finding each pixels and then taking the median length of the flowpaths originating from those edge pixels.
            let edge_flowpaths = self.get_edge_flowpaths();
            length = flowpaths_median_length(&edge_flowpaths).unwrap_or(0.0);
            width = area / length;
        } else {
            // Otherwise the  width of the subcatchment is determined by the
            // channel that the subcatchment drains into. The length is
            // then determined by the area / width

            width = chn_summary.length;
            length = area / width;
        }

        let mut direction: f64 = chn_summary.direction;
        match topaz_id % 10 {
            2 => direction += 90.0,
            3 => direction -= 90.0,
            _ => (),
        }

        // determine aspect
        let aspect: f64 = taspec.determine_aspect(&indices);

        // calculate weighted slope from flowpaths
        let (w_slopes, distances, distances_norm): (Vec<f64>, Vec<f64>, Vec<f64>) =
            self.weighted_slope_average_from_fps();

        let longest_fp: &FlowPath = self.get_longest_fp();
        let centroid_px = subwta.centroid_of(&indices);

        assert!(distances.len() > 1, "distances {:?}", distances);

        // iterate over distances and slopes and calculate elevations
        // for each point
        let mut elevs: Vec<f64> = vec![longest_fp.elevation];
        for i in 0..distances.len() - 1 {
            let dx: f64 = distances[i+1] - distances[i];
            let dy: f64 = w_slopes[i];
            let elevation: f64 = elevs[i] - (dx * dy);
            elevs.push(elevation);
        }

        let slope_scalar: f64 = (elevs[0] - elevs[elevs.len() - 1]) / length;
        let elevation: f64 = elevs[0];

        let vec_indices: Vec<usize> = indices.into_iter().collect();
        FlowPath::new(
            vec_indices,
            longest_fp.head,
            longest_fp.tail,
            (centroid_px.0 as i32, centroid_px.1 as i32),
            distances_norm,
            w_slopes,
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
            -1,
            -1.0
        )
    }

    #[allow(dead_code)]
    pub fn abstract_hillslope(&self,
        flovec: &Raster<i32>,
        taspec: &Raster<f64>,
        vec_indices: &Vec<usize>
    ) -> FlowPath {

        let cellsize: f64 = flovec.cellsize;
        let cellsize2: f64 = cellsize * cellsize;
        let topaz_id: i32 = self.flowpaths[0].topaz_id;

        let indices: HashSet<usize> = vec_indices.iter().cloned().collect();

        // get indices of subcatchment
        let area = indices.len() as f64 * cellsize2;

        let longest_fp: &FlowPath = self.get_longest_fp();

        // If subcatchment is a source type then we calculate the distance
        // by taking a weighted average based on the length of the flowpaths
        // contained in the subcatchment
        //let length: f64 = longest_fp.length;
        let edge_flowpaths = self.get_edge_flowpaths();
        let length = flowpaths_median_length(&edge_flowpaths).unwrap_or(0.0);
        let width: f64 = area / length;

        let mut direction: f64 = 0.0;

        // determine aspect
        let aspect: f64 = taspec.determine_aspect(vec_indices);

        // calculate weighted slope from flowpaths
        let (w_slopes, distances, distances_norm): (Vec<f64>, Vec<f64>, Vec<f64>) =
            self.weighted_slope_average_from_fps();

        let centroid_px = flovec.centroid_of(vec_indices);

        assert!(distances.len() > 1, "distances {:?}", distances);

        // iterate over distances and slopes and calculate elevations
        // for each point
        let mut elevs: Vec<f64> = vec![longest_fp.elevation];
        for i in 0..distances.len() - 1 {
            let dx: f64 = distances[i+1] - distances[i];
            let dy: f64 = w_slopes[i];
            let elevation: f64 = elevs[i] - (dx * dy);
            elevs.push(elevation);
        }

        let slope_scalar: f64 = (elevs[0] - elevs[elevs.len() - 1]) / length;
        let elevation: f64 = elevs[0];

        let vec_indices: Vec<usize> = indices.into_iter().collect();
        FlowPath::new(
            vec_indices,
            longest_fp.head,
            longest_fp.tail,
            (centroid_px.0 as i32, centroid_px.1 as i32),
            distances_norm,
            w_slopes,
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
            -1,
            -1.0
        )
    }

    #[allow(dead_code)]
    pub fn to_geojson_feature_collection(&self, raster: &Raster<i32>) -> FeatureCollection {
        let features: Vec<Feature> = self.flowpaths.iter()
            .map(|fp| fp.to_geojson_feature(raster))
            .collect();

        FeatureCollection {
            bbox: None,
            features: features,
            foreign_members: None,
        }
    }

    #[allow(dead_code)]
    pub fn to_geojson(&self, raster: &Raster<i32>) -> String {
        let feature_collection = self.to_geojson_feature_collection(raster);

        let mut geojson = serde_json::Map::new();
        geojson.insert(String::from("type"),
        serde_json::Value::String(String::from("FeatureCollection")));

        geojson.insert(
            String::from("features"),
            serde_json::Value::Array(
                feature_collection.features.into_iter().map(|f| serde_json::to_value(f).unwrap()).collect()
            )
        );

        // Optionally add CRS if it exists in the raster
        if let Some(proj) = &raster.proj4 {
            let mut crs = serde_json::Map::new();
            crs.insert(String::from("type"), serde_json::Value::String(String::from("name")));
            crs.insert(String::from("properties"), serde_json::Value::Object(serde_json::Map::from_iter(
                std::iter::once((String::from("name"), serde_json::Value::String(proj.clone())))
            )));
            geojson.insert(String::from("crs"), serde_json::Value::Object(crs));
        }

        to_string_pretty(&serde_json::Value::Object(geojson)).unwrap()
    }

    #[allow(dead_code)]
    pub fn write_geojson(&self, raster: &Raster<i32>, path: &str) -> std::io::Result<()> {
        let geojson_string = self.to_geojson(raster);

        // Open a file in write mode
        let mut file = File::create(path)?;

        // Write the GeoJSON string to the file
        file.write_all(geojson_string.as_bytes())?;

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_subflow_slps(&self, out_dir: &str, max_points: usize, clip_hillslopes: bool, clip_hillslope_length: f64) -> std::io::Result<()> {
        if let Some(subflows_map) = &self.subflows {
            subflows_map.par_iter().for_each(|(topaz_id, subflowpath_collection)| {
                subflowpath_collection.write_fp_slps(&out_dir, max_points).unwrap();
            });
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_channel_slp(&self, path: &str, max_points: usize) -> std::io::Result<()> {

        let mut all_strings = Vec::new();
        all_strings.push(format!("2023.1\n{}\n", &self.flowpaths.len()));

        for fp in self.flowpaths.iter().rev() {

            let simplified: (Vec<f64>, Vec<f64>);
            let (d0, s0) = if fp.distances_norm.len() > 3 {
                simplified = douglas_peucker(&fp.distances_norm, &fp.slopes, 0.01).unwrap();
                println!("douglas_peuker {} -> {}", fp.distances_norm.len(), simplified.0.len());
                (&simplified.0, &simplified.1)
            } else {
                (&fp.distances_norm, &fp.slopes)
            };

            let interpolated: (Vec<f64>, Vec<f64>);
            let (d, s) = if d0.len() > max_points {
                interpolated = interpolate_slp(&d0, &s0, max_points).unwrap();
                (&interpolated.0, &interpolated.1)
            } else {
                (d0, s0)
            };

            let npts: usize = d.len();

            // Build the defs string
            let defs: Vec<String> = d.iter()
                .zip(s.iter())
                .map(|(&dist, &slope)| format!("{:.4}, {:.4}", dist, slope))
                .collect();

            let slp = format!(
                "{:.4} {:.1} {:.1} {}\n{} {:.1}\n{} \n",
                fp.aspect, fp.width, fp.elevation, fp.order, npts, fp.length, defs.join(" ")
            );
            all_strings.push(slp);
        }
        let contents = all_strings.join("");

        // Write to file
        let mut file = std::fs::File::create(path)?;
        file.write_all(contents.as_bytes())?;

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_slps(&self, out_dir: &str, max_points: usize, 
        clip_hillslopes: bool, clip_hillslope_length: f64) -> std::io::Result<()> {

        let results: Vec<std::io::Result<()>> = self.flowpaths.par_iter()
            .map(|fp| {
                let fname;
                if fp.fp_id == -1 {
                    if fp.topaz_id % 10 == 4 {
                        fname = format!("chn_{}.slp", fp.topaz_id);
                    } else {
                        fname = format!("hill_{}.slp", fp.topaz_id);
                    }
                } else {
                    fname = format!("fp_{}_{}.slp", fp.topaz_id, fp.fp_id);
                }

                let path = format!("{}/{}", out_dir, fname);
                fp.write_slp(&path, max_points, clip_hillslopes, clip_hillslope_length)
            })
            .collect();

        // Check for any errors
        for result in results {
            result?;
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_fp_slps(&self, out_dir: &str, max_points: usize) -> Result<()> {

        let fname = format!("fps_{}.slps", self.flowpaths[0].topaz_id);
        let path = format!("{}/{}", out_dir, fname);

        // Create the file before writing
        let mut file: File = File::create(&path)?;

        for fp in &self.flowpaths {
            // Ensure the directory exists
            let path_dir = Path::new(&out_dir);
            if !path_dir.exists() {
                std::fs::create_dir_all(path_dir)?;
            }
    
            let mut slp = String::new();
            slp.push_str(&format!("# fp_{}_{}.slp\n", fp.topaz_id, fp.fp_id));
            slp.push_str(&format!("# {:?}\n", fp.indices));
            slp.push_str(&format!("# {:?}\n", fp.distances_norm));
            let _ = file.write(slp.as_bytes());
            
            // Write the SLP file
            fp._write_slp(&file, max_points, false, 300.0)?;

            // Add a newline after each flowpath
            file.write(b"\n")?;
        }
    
        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_chn_metadata_to_csv(&self, path: &str, wgs_transform: &[f64; 4]) -> std::io::Result<()> {
        let file = File::create(path).unwrap();
        let mut writer = csv::Writer::from_writer(file);

        let headers: Vec<String> = vec![
            String::from("topaz_id"),
            String::from("slope_scalar"),
            String::from("length"),
            String::from("width"),
            String::from("direction"),
            String::from("order"),
            String::from("aspect"),
            String::from("area"),
            String::from("elevation"),
            String::from("centroid_px"),
            String::from("centroid_py"),
            String::from("centroid_lon"),
            String::from("centroid_lat"),
        ];

        writer.write_record(headers).unwrap();

        for fp in &self.flowpaths {
            let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);

            let record: Vec<String> = vec![
                fp.topaz_id.to_string(),
                fp.slope_scalar.to_string(),
                fp.length.to_string(),
                fp.width.to_string(),
                fp.direction.to_string(),
                fp.order.to_string(),
                fp.aspect.to_string(),
                (fp.area_m2() as i32).to_string(),
                fp.elevation.to_string(),
                fp.centroid_px.0.to_string(),
                fp.centroid_px.1.to_string(),
                lon.to_string(),
                lat.to_string(),
            ];

            writer.write_record(record).unwrap();
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_metadata_to_csv(&self, path: &str, wgs_transform: &[f64; 4]) -> std::io::Result<()> {
        let file = File::create(path).unwrap();
        let mut writer = csv::Writer::from_writer(file);

        let headers: Vec<String> = vec![
            String::from("topaz_id"),
            String::from("slope_scalar"),
            String::from("length"),
            String::from("width"),
            String::from("direction"),
            String::from("aspect"),
            String::from("area"),
            String::from("elevation"),
            String::from("centroid_px"),
            String::from("centroid_py"),
            String::from("centroid_lon"),
            String::from("centroid_lat"),
        ];

        writer.write_record(headers).unwrap();

        for fp in &self.flowpaths {
            let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);

            let record: Vec<String> = vec![
                fp.topaz_id.to_string(),
                fp.slope_scalar.to_string(),
                fp.length.to_string(),
                fp.width.to_string(),
                fp.direction.to_string(),
                fp.aspect.to_string(),
                (fp.area_m2() as i32).to_string(),
                fp.elevation.to_string(),
                fp.centroid_px.0.to_string(),
                fp.centroid_px.1.to_string(),
                lon.to_string(),
                lat.to_string(),
            ];

            writer.write_record(record).unwrap();
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_subflows_metadata_to_csv(&self, path: &str, wgs_transform: &[f64; 4]) -> std::io::Result<()> {
        let file = File::create(path).unwrap();
        let mut writer = csv::Writer::from_writer(file);

        let headers: Vec<String> = vec![
            String::from("topaz_id"),
            String::from("fp_id"),
            String::from("slope_scalar"),
            String::from("length"),
            String::from("width"),
            String::from("direction"),
            String::from("aspect"),
            String::from("area"),
            String::from("elevation"),
            String::from("order"),
            String::from("centroid_px"),
            String::from("centroid_py"),
            String::from("centroid_lon"),
            String::from("centroid_lat"),
        ];

        writer.write_record(headers).unwrap();

        if let Some(subflows_map) = &self.subflows {
            for (topaz_id, subflow_collection) in subflows_map {
                for fp in subflow_collection.flowpaths.iter() {
                    let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);

                    let record: Vec<String> = vec![
                        topaz_id.to_string(),
                        fp.fp_id.to_string(),
                        fp.slope_scalar.to_string(),
                        fp.length.to_string(),
                        fp.width.to_string(),
                        fp.direction.to_string(),
                        fp.aspect.to_string(),
                        fp.area_m2().to_string(),
                        fp.elevation.to_string(),
                        fp.order.to_string(),
                        fp.centroid_px.0.to_string(),
                        fp.centroid_px.1.to_string(),
                        lon.to_string(),
                        lat.to_string(),
                    ];

                    writer.write_record(record).unwrap();
                }
            }
        } else {
            // This will panic if the None case occurs
            panic!("Unexpected None in subflows!");
        }

        Ok(())
    }

    // Method to find edge flowpaths
    pub fn get_edge_flowpaths(&self) -> Vec<FlowPath> {
        let mut edge_flowpaths = Vec::new();

        // Iterate over flowpaths with their index for comparison
        for (i, flowpath) in self.flowpaths.iter().enumerate() {
            // Assume current flowpath is an edge until proven otherwise
            let mut is_edge = true;
            let first_index = match flowpath.indices.first() {
                Some(&index) => index,
                None => continue, // Skip if flowpath has no indices
            };

            // Check against all other flowpaths
            for (j, other_flowpath) in self.flowpaths.iter().enumerate() {
                if i == j { continue; } // Skip self-comparison

                // If first_index is found in any other flowpath, it's not an edge
                if other_flowpath.indices.contains(&first_index) {
                    is_edge = false;
                    break; // Early exit if we find a match
                }
            }

            // If is_edge remains true, add to result vector
            if is_edge {
                edge_flowpaths.push(self.flowpaths[i].clone());
            }
        }

        edge_flowpaths
    }
}


// Method to calculate the median length of flowpaths
fn flowpaths_median_length(flowpaths: &Vec<FlowPath>) -> Option<f64> {
    let num_flowpaths = flowpaths.len();

    // Return None if the vector is empty
    if num_flowpaths == 0 {
        return None;
    }

    // Extract the lengths and sort them
    let mut lengths: Vec<f64> = flowpaths.iter().map(|fp| fp.length).collect();
    lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Calculate median based on odd or even number of flowpaths
    if num_flowpaths % 2 == 0 {
        // Even number of flowpaths: average the two middle elements
        let mid_right = num_flowpaths / 2;
        let mid_left = mid_right - 1;
        Some((lengths[mid_left] + lengths[mid_right]) / 2.0)
    } else {
        // Odd number of flowpaths: middle element
        Some(lengths[num_flowpaths / 2])
    }
}


#[derive(Debug, Clone)]
pub struct FlowPath {
    pub indices: Vec<usize>,
    pub head: (i32, i32),
    pub tail: (i32, i32),
    pub centroid_px: (i32, i32),
    pub distances_norm: Vec<f64>,
    pub slopes: Vec<f64>,
    pub elevs: Vec<f64>,
    pub topaz_id: i32,
    pub fp_id: i32,
    pub length: f64,
    pub width: f64,
    pub aspect: f64,
    pub direction: f64,
    pub slope_scalar: f64,
    pub cellsize: f64,
    pub slopes_r: Vec<f64>,
    pub distance_to_chn_r: Vec<f64>,
    pub kp: f64,
    pub elevation: f64,
    pub order: i32,
    pub areaup: f64,
}

impl Default for FlowPath {
    fn default() -> Self {
        FlowPath {
            indices: Vec::new(),
            centroid_px: (0, 0),
            head: (0, 0),
            tail: (0, 0),
            distances_norm: Vec::new(),
            slopes: Vec::new(),
            elevs: Vec::new(),
            topaz_id: 0,
            fp_id: -1,
            length: -1.0,
            width: -1.0,
            aspect: -1.0,
            direction: -1.0,
            slope_scalar: -1.0,
            cellsize: -1.0,
            slopes_r: Vec::new(),
            distance_to_chn_r: Vec::new(),
            kp: -1.0,
            elevation: -1.0,
            order: -1,
            areaup: -1.0,
        }
    }
}

impl FlowPath {

    #[allow(dead_code)]
    pub fn new(
        indices: Vec<usize>,
        head: (i32, i32),
        tail: (i32, i32),
        centroid_px: (i32, i32),
        distances_norm: Vec<f64>,
        slopes: Vec<f64>,
        elevs: Vec<f64>,
        topaz_id: i32,
        fp_id: i32,
        length: f64,
        width: f64,
        aspect: f64,
        direction: f64,
        slope_scalar: f64,
        cellsize: f64,
        elevation: f64,
        order: i32,
        areaup: f64
    ) -> Self {
        let mut slopes_r = slopes.clone();
        slopes_r.reverse();

        let rcd: Vec<f64> = distances_norm.iter().map(|&x| 1.0 - x).collect();
        let mut distance_to_chn_r: Vec<f64> = rcd.iter().map(|&x| length * x).collect();
        distance_to_chn_r.reverse();

        let kp: f64 = length * indices.len() as f64 * cellsize * cellsize;

        FlowPath {
            indices,
            head,
            tail,
            centroid_px,
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
            slopes_r,
            distance_to_chn_r,
            kp,
            elevation,
            order,
            areaup
        }
    }

    #[allow(dead_code)]
    pub fn interp_slp_at_distance_to_channel(&self, distance_to_chn: f64) -> f64 {
        let x = interp(&self.distance_to_chn_r, &self.slopes_r, distance_to_chn);

        if x.is_finite() {
            return x as f64;
        }
        // return slopes_r value at index where distance_to_chn is closest to self.distance_to_chn_r
        let mut min_dist: f64 = 1e10;
        let mut min_index: usize = 0;
        for (i, &d) in self.distance_to_chn_r.iter().enumerate() {
            let dist = (d - distance_to_chn).abs();
            if dist < min_dist {
                min_dist = dist;
                min_index = i;
            }
        }

        self.slopes_r[min_index]
        
    }

    #[allow(dead_code)]
    pub fn distances(&self) -> Vec<f64> {
        self.distances_norm.iter().map(|&d| d * self.length).collect()
    }

    #[allow(dead_code)]
    pub fn area_m2(&self) -> f64 {
        self.width * self.length
    }

    #[allow(dead_code)]
    pub fn write_slp(&self,
        path: &str,
        max_points: usize,
        clip_hillslopes: bool,
        clip_hillslope_length: f64
    ) -> std::io::Result<()> {
        
        // Open a file in write mode
        let mut file = File::create(path)?;

        self._write_slp(&file, max_points, clip_hillslopes, clip_hillslope_length)?;

        Ok(())
    }

    #[allow(dead_code)]
    pub fn _write_slp(&self,
        mut file: &File,
        max_points: usize,
        clip_hillslopes: bool,
        clip_hillslope_length: f64
    ) -> std::io::Result<()> {

        let simplified: (Vec<f64>, Vec<f64>);
        let (d0, s0) = if self.distances_norm.len() > 3 {
            simplified = douglas_peucker(&self.distances_norm, &self.slopes, 0.01).unwrap();
            (&simplified.0, &simplified.1)
        } else {
            (&self.distances_norm, &self.slopes)
        };

        let interpolated: (Vec<f64>, Vec<f64>);
        let (d, s) = if d0.len() > max_points {
            interpolated = interpolate_slp(&d0, &s0, max_points).unwrap();
            (&interpolated.0, &interpolated.1)
        } else {
            (d0, s0)
        };

        let nofes: i32 = 1;
        let npts: usize = d.len();

        // Build the defs string
        let defs: Vec<String> = d.iter()
            .zip(s.iter())
            .map(|(&dist, &slope)| format!("{:.4}, {:.4}", dist, slope))
            .collect();

        let mut width = self.width;
        let mut length = self.length;

        if clip_hillslopes {
            if length > clip_hillslope_length {
                length = clip_hillslope_length;
                width = self.area_m2() / length;
            }
        }

        let slp = format!(
            "2023.3\n{}\n{:.4} {:.1} {:.1}\n{} {:.1}\n{} ",
            nofes, self.aspect, width, self.elevation, npts, length, defs.join(" ")
        );
        
        file.write_all(slp.as_bytes())?;
    
        Ok(())
    }

    #[allow(dead_code)]
    pub fn to_geojson_feature(&self, raster: &Raster<i32>) -> Feature {
        let coords = raster.coordinates_of(&self.indices);
        let line_string = Geometry::new(LineString(coords));
        let mut properties = serde_json::Map::new();
        properties.insert(String::from("topaz_id"),
        serde_json::Value::Number(
            serde_json::Number::from(
                self.topaz_id as i64)));

        if let Some(number) = serde_json::Number::from_f64(self.width) {
            properties.insert(String::from("width"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.length) {
            properties.insert(String::from("length"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.aspect) {
            properties.insert(String::from("aspect"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.direction) {
            properties.insert(String::from("direction"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.slope_scalar) {
            properties.insert(String::from("slope_scalar"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.cellsize) {
            properties.insert(String::from("cellsize"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.order as f64) {
            properties.insert(String::from("order"), serde_json::Value::Number(number));
        }

        if let Some(number) = serde_json::Number::from_f64(self.areaup) {
            properties.insert(String::from("areaup"), serde_json::Value::Number(number));
        }

        Feature {
            bbox: None,
            geometry: Some(line_string),
            id: None,
            properties: Some(properties),
            foreign_members: None,
        }
    }
}

#[allow(dead_code)]
pub fn walk_channels<T: LinkAttributes>(
    subwta: &Raster<i32>,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
    netw:   &HashMap<i32, T>,
    bieger2015_widths: bool
) -> FlowpathCollection
where
    T: LinkAttributes + Sync, {

    let unique_vals = subwta.unique_values();
    let mut topaz_ids: Vec<i32> = unique_vals.iter().filter(|&&topaz_id| topaz_id % 10 == 4).cloned().collect();
    topaz_ids.sort();

    let flowpaths: Vec<FlowPath> = topaz_ids.into_par_iter()
        .map(|topaz_id| walk_channel(
            topaz_id, 
            &subwta, 
            &relief, 
            &flovec, 
            &fvslop, 
            &taspec, 
            netw[&topaz_id].get_areaup(), 
            netw[&topaz_id].get_order(), 
            bieger2015_widths))
        .collect();

    FlowpathCollection {
        flowpaths: flowpaths,
        subflows: None
    }
}

#[allow(dead_code)]
pub fn walk_channel(
    topaz_id: i32,
    subwta: &Raster<i32>,
    relief: &Raster<f64>,
    _flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
    areaup: f64,
    order: i32,
    bieger2015_widths: bool
) -> FlowPath {
    // get hashset of indices of topaz_id
    let indices: HashSet<usize> = subwta.indices_of(topaz_id);
    assert!(indices.len() > 0, "indices.len() == 0 for topaz_id: {}", topaz_id);

    let cellsize = subwta.cellsize;

    // build hashmap of indices and elevations
    let mut elevs_d: HashMap<usize, f64> = HashMap::new();
    for i in &indices {
        elevs_d.insert(*i, relief.data[*i]);
    }

    // build Vec<i32> of indices in descending order of elevation
    let mut indices_elevs: Vec<_> = elevs_d.iter().collect();
    indices_elevs.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));

    let sorted_indices: Vec<usize> = indices_elevs.clone().into_iter().map(|(index, _)| *index).collect();
    let elevs: Vec<f64> = indices_elevs.clone().into_iter().map(|(_, elev)| *elev).collect();
    let elevation: f64 = elevs[0];

    let mut n: usize = sorted_indices.len();

    // calculate distances between the indices
    let mut distances: Vec<f64> = vec![0.0];
    for i in 1..n {
        let index1 = sorted_indices[i - 1];
        let index2 = sorted_indices[i];
        let distance = subwta.distance_between(index1, index2);
        assert!(distance > 0.0, "distance is 0.0 for topaz_id: {}", topaz_id);
        distances.push(distances[i-1] + distance);
    }

    let mut slopes: Vec<f64> = Vec::new();
    let mut rad_aspects: Vec<f64> = Vec::new();
    for i in 0..n {
        let index = sorted_indices[i];
        let slope = fvslop.data[index];
        slopes.push(slope);

        let deg_aspect = taspec.data[index];
        rad_aspects.push(deg_aspect.to_radians());
    }

    let aspect: f64 = circmean(rad_aspects.as_slice()).to_degrees();

    let (centroid_x, centroid_y) = subwta.centroid_of(&indices);

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
        assert!(total_elev > 0.0, "total_elev is 0.0 for topaz_id: {}", topaz_id);
        slope_scalar = total_elev / distances[n - 1];
    }

    let length: f64 = distances[n - 1];
    assert!(length > 0.0, "length is 0.0 for topaz_id: {}, {:?} {:?}", topaz_id, distances, slopes);

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
        let da: f64 = areaup * 1e-6;  // area in km^2
        // width = 2.70 * da.powf(0.352); // USA model https://github.com/rogerlew/wepppy/issues/268
        width = 1.24 * da.powf(0.435); // Rocky Mountain System
    }

    FlowPath::new(
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
        areaup
    )
}

#[allow(dead_code)]
pub fn walk_flowpaths(
    topaz_id: i32,
    subwta: &Raster<i32>,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
) -> FlowpathCollection {

    let indices = subwta.indices_of(topaz_id);
    let indices_vec: Vec<_> = indices.iter().collect();
    let flowpaths: Vec<FlowPath> = indices_vec.iter()
        .enumerate()
        .map(|(i, &head_index)| {
            walk_flowpath(
                topaz_id,
                *head_index,
                &subwta,
                &relief,
                &flovec,
                &fvslop,
                &taspec,
                i.try_into().unwrap()
            )
        })
        .collect();

    FlowpathCollection {
        flowpaths: flowpaths,
        subflows: None
    }
}

#[allow(dead_code)]
pub fn walk_flowpath(
    topaz_id: i32,
    head_index: usize,
    subwta: &Raster<i32>,
    relief: &Raster<f64>,
    flovec: &Raster<i32>,
    fvslop: &Raster<f64>,
    taspec: &Raster<f64>,
    fp_id: i32
) -> FlowPath {

    assert_eq!(subwta.data[head_index], topaz_id);

    let cellsize: f64 = subwta.cellsize;
    
    let mut current_index: usize = head_index;
    let mut sorted_indices: Vec<usize> = vec![head_index];
    let mut indices_hash: HashSet<usize> = HashSet::new();
    indices_hash.insert(head_index);
    let mut distances: Vec<f64> = vec![0.0];
    let mut i = 0;
    loop {
        let flow_dir: i32 = flovec.data[current_index];
        let (dx, dy) = *PATHS.get(&flow_dir).expect(&format!("flow_dir {} not found in paths!", flow_dir));

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
        distances.push(distances[distances.len() - 1] + cellsize * ((dx as f64).abs() + (dy as f64).abs()).sqrt());

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
        let slope: f64 = fvslop.data[index];
        slopes.push(slope);

        let deg_aspect = taspec.data[index];
        rad_aspects.push(deg_aspect.to_radians());

        let elev: f64 = relief.data[index];
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
        -1.0
    )
}

