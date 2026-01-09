use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Result, Write};
use std::path::Path;

use geojson::{Feature, FeatureCollection};
use rayon::prelude::*;
use serde_json::to_string_pretty;

use crate::douglas_peucker::douglas_peucker;
use crate::raster::{px_to_wgs, Raster};
use crate::support::interpolate_slp;
use crate::flowpath::Flowpath;
use crate::watershed_abstraction::PATHS;

#[derive(Debug, Clone)]
pub struct FlowpathCollection {
    pub flowpaths: Vec<Flowpath>,
    pub subflows: Option<HashMap<i32, FlowpathCollection>>
}

impl FlowpathCollection {

    #[allow(dead_code)]
    pub fn get_fp_by_topaz_id(&self, topaz_id: i32) -> Option<&Flowpath> {
        self.flowpaths.iter().find(|fp| fp.topaz_id == topaz_id)
    }

    #[allow(dead_code)]
    pub fn get_longest_fp(&self) -> &Flowpath {
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

    fn resolve_field_lookup(
        fake_topaz_id: i32,
        fake_topaz_id_lookup: &HashMap<(i32, i32), i32>,
    ) -> Option<(i32, i32)> {
        fake_topaz_id_lookup
            .iter()
            .find_map(|(&(field_id, topaz_id), &v)| {
                if v == fake_topaz_id {
                    Some((field_id, topaz_id))
                } else {
                    None
                }
            })
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

        let longest_fp: &Flowpath = self.get_longest_fp();
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
        taspec: &Raster<f32>,
        flovec: &Raster<u8>,
        channels: &FlowpathCollection,
        indices: &Vec<usize>) -> Flowpath {

        let cellsize: f64 = subwta.cellsize;
        let cellsize2: f64 = cellsize * cellsize;
        let topaz_id: i32 = self.flowpaths[0].topaz_id;
        let area = indices.len() as f64 * cellsize2;

        // find corresponding chn_id
        let chn_id: i32 = ((topaz_id as f64 / 10.0).floor() * 10.0) as i32 + 4;
        let chn_summary: &Flowpath = &channels.get_fp_by_topaz_id(chn_id).unwrap();

        // If subcatchment is a source type then we calculate the distance
        // by taking a weighted average based on the length of the flowpaths
        // contained in the subcatchment
        let length: f64;
        let width: f64;
        if topaz_id % 10 == 1 {
            // find length by finding each pixels and then taking the median length of the flowpaths originating from those edge pixels.
            let edge_flowpaths = self.get_edge_flowpaths2(subwta, flovec);
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
        let aspect: f64 = taspec.determine_aspect(indices);

        // calculate weighted slope from flowpaths
        let (w_slopes, distances, distances_norm): (Vec<f64>, Vec<f64>, Vec<f64>) =
            self.weighted_slope_average_from_fps();

        let longest_fp: &Flowpath = self.get_longest_fp();
        let centroid_px = subwta.centroid_of(indices);

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

        let vec_indices: Vec<usize> = indices.clone();
        Flowpath::new(
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
    pub fn abstract_subfieldcatchment(&self,
        intersection_subwta: &Raster<i32>,
        taspec: &Raster<f32>,
        flovec: &Raster<u8>,
        indices: &Vec<usize>
    ) -> Flowpath {


        let cellsize: f64 = intersection_subwta.cellsize;
        let cellsize2: f64 = cellsize * cellsize;
        let fake_topaz_id: i32 = self.flowpaths[0].topaz_id;
        let area = indices.len() as f64 * cellsize2;

        let longest_fp = self.get_longest_fp();

        let length: f64;
        let width: f64;
        
        // find length by finding each pixels and then taking the median length of the flowpaths originating from those edge pixels.
        let edge_flowpaths = self.get_edge_flowpaths2(intersection_subwta, flovec);
        length = flowpaths_median_length(&edge_flowpaths).unwrap_or(0.0);
        width = area / length;

        let direction : f64 = longest_fp.direction;

        // determine aspect
        let aspect: f64 = taspec.determine_aspect(indices);

        // calculate weighted slope from flowpaths
        let (w_slopes, distances, distances_norm): (Vec<f64>, Vec<f64>, Vec<f64>) =
            self.weighted_slope_average_from_fps();

        let longest_fp: &Flowpath = self.get_longest_fp();
        let centroid_px = intersection_subwta.centroid_of(indices);

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

        let vec_indices: Vec<usize> = indices.clone();
        Flowpath::new(
            vec_indices,
            longest_fp.head,
            longest_fp.tail,
            (centroid_px.0 as i32, centroid_px.1 as i32),
            distances_norm,
            w_slopes,
            elevs,
            fake_topaz_id,
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
    pub fn write_subflow_slps(&self, out_dir: &str, max_points: usize) -> std::io::Result<()> {
        if let Some(subflows_map) = &self.subflows {
            subflows_map.par_iter().for_each(|(topaz_id, subflowpath_collection)| {
                let path = format!("{}/fps_{}.slps", out_dir, topaz_id);
                subflowpath_collection.write_fp_slps(&path, max_points).unwrap();
            });
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_channel_slp(&self, path: &str, max_points: usize) -> std::io::Result<()> {

        let mut all_strings = Vec::new();
        all_strings.push(format!("2025.8\n{}\n", &self.flowpaths.len()));

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
    pub fn write_field_slps(
        &self,
        out_dir: &str,
        max_points: usize,
        clip_hillslopes: bool,
        clip_hillslope_length: f64,
        fake_topaz_id_lookup: &HashMap<(i32, i32), i32>,
    ) -> std::io::Result<()> {
        for fp in &self.flowpaths {
            if let Some((field_id, topaz_id)) = Self::resolve_field_lookup(fp.topaz_id, fake_topaz_id_lookup)
            {
                let fname = if fp.fp_id == -1 {
                    format!("field_{}_{}.slp", field_id, topaz_id)
                } else {
                    format!("field_{}_{}_fp_{}.slp", field_id, topaz_id, fp.fp_id)
                };
                let path = format!("{}/{}", out_dir, fname);
                fp.write_slp(&path, max_points, clip_hillslopes, clip_hillslope_length)?;
            }
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_fp_slps(&self, path: &str, max_points: usize) -> Result<()> {
        if let Some(parent) = Path::new(path).parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                std::fs::create_dir_all(parent)?;
            }
        }

        let mut file: File = File::create(path)?;

        for fp in &self.flowpaths {
            let mut slp = String::new();
            slp.push_str(&format!("# fp_{}_{}.slp\n", fp.topaz_id, fp.fp_id));
            slp.push_str(&format!("# {:?}\n", fp.indices));
            slp.push_str(&format!("# {:?}\n", fp.distances_norm));
            let _ = file.write(slp.as_bytes());

            fp._write_slp(&file, max_points, false, 300.0)?;

            file.write(b"\n")?;
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub fn write_field_subflow_slps(
        &self,
        out_dir: &str,
        max_points: usize,
        fake_topaz_id_lookup: &HashMap<(i32, i32), i32>,
    ) -> std::io::Result<()> {
        if let Some(subflows_map) = &self.subflows {
            for (fake_topaz_id, subflow_collection) in subflows_map {
                if let Some((field_id, topaz_id)) =
                    Self::resolve_field_lookup(*fake_topaz_id, fake_topaz_id_lookup)
                {
                    let path = format!("{}/field_{}_{}.slps", out_dir, field_id, topaz_id);
                    subflow_collection.write_fp_slps(&path, max_points)?;
                }
            }
        } else {
            panic!("Unexpected None in subflows!");
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
                fp.area_m2().to_string(),
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
    pub fn write_field_metadata_to_csv(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
        fake_topaz_id_lookup: &HashMap<(i32, i32), i32>,
    ) -> std::io::Result<()> {
        let file = File::create(path).unwrap();
        let mut writer = csv::Writer::from_writer(file);

        writer
            .write_record([
                "field_id",
                "topaz_id",
                "sub_field_id",
                "slope_scalar",
                "length",
                "width",
                "direction",
                "aspect",
                "area",
                "elevation",
                "centroid_px",
                "centroid_py",
                "centroid_lon",
                "centroid_lat",
            ])
            .unwrap();

        for fp in &self.flowpaths {
            if let Some((field_id, topaz_id)) = Self::resolve_field_lookup(fp.topaz_id, fake_topaz_id_lookup)
            {
                let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);

                let record = [
                    field_id.to_string(),
                    topaz_id.to_string(),
                    fp.topaz_id.to_string(),
                    fp.slope_scalar.to_string(),
                    fp.length.to_string(),
                    fp.width.to_string(),
                    fp.direction.to_string(),
                    fp.aspect.to_string(),
                    fp.area_m2().to_string(),
                    fp.elevation.to_string(),
                    fp.centroid_px.0.to_string(),
                    fp.centroid_px.1.to_string(),
                    lon.to_string(),
                    lat.to_string(),
                ];

                writer.write_record(record).unwrap();
            }
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

    #[allow(dead_code)]
    pub fn write_field_subflows_metadata_to_csv(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
        fake_topaz_id_lookup: &HashMap<(i32, i32), i32>,
    ) -> std::io::Result<()> {
        let file = File::create(path).unwrap();
        let mut writer = csv::Writer::from_writer(file);

        writer
            .write_record([
                "field_id",
                "topaz_id",
                "sub_field_id",
                "topaz_id",
                "fp_id",
                "slope_scalar",
                "length",
                "width",
                "direction",
                "aspect",
                "area",
                "elevation",
                "order",
                "centroid_px",
                "centroid_py",
                "centroid_lon",
                "centroid_lat",
            ])
            .unwrap();

        if let Some(subflows_map) = &self.subflows {
            for (fake_topaz_id, subflow_collection) in subflows_map {
                if let Some((field_id, topaz_id)) =
                    Self::resolve_field_lookup(*fake_topaz_id, fake_topaz_id_lookup)
                {
                    for fp in &subflow_collection.flowpaths {
                        let (lon, lat) =
                            px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);

                        let record = [
                            field_id.to_string(),
                            topaz_id.to_string(),
                            fake_topaz_id.to_string(),
                            fp.topaz_id.to_string(),
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
            }
        } else {
            panic!("Unexpected None in subflows!");
        }

        Ok(())
    }

    // Method to find edge flowpaths
    #[deprecated(note = "Use get_edge_flowpaths2 with subwta+flovec for faster, source-cell detection.")]
    pub fn get_edge_flowpaths(&self) -> Vec<Flowpath> {
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

    #[allow(dead_code)]
    pub fn get_edge_flowpaths2(
        &self,
        subwta: &Raster<i32>,
        flovec: &Raster<u8>,
    ) -> Vec<Flowpath> {
        let mut mask_indices: HashSet<usize> = HashSet::new();
        for fp in &self.flowpaths {
            for &index in &fp.indices {
                mask_indices.insert(index);
            }
        }

        let mut has_upstream: HashSet<usize> = HashSet::new();
        for &index in &mask_indices {
            let flow_dir = flovec.data[index] as i32;
            if flow_dir == 0 {
                continue;
            }

            let Some((dx, dy)) = PATHS.get(&flow_dir) else {
                continue;
            };
            let (x, y) = subwta.index_to_xy(index);
            let next_x = x as isize + dx;
            let next_y = y as isize + dy;
            if next_x < 0
                || next_y < 0
                || next_x >= subwta.width as isize
                || next_y >= subwta.height as isize
            {
                continue;
            }
            let next_index = subwta.xy_to_index(next_x as usize, next_y as usize);
            if mask_indices.contains(&next_index) {
                has_upstream.insert(next_index);
            }
        }

        let mut edge_flowpaths = Vec::new();
        for fp in &self.flowpaths {
            let first_index = match fp.indices.first() {
                Some(&index) => index,
                None => continue,
            };
            if !has_upstream.contains(&first_index) {
                edge_flowpaths.push(fp.clone());
            }
        }

        edge_flowpaths
    }
}


// Method to calculate the median length of flowpaths
fn flowpaths_median_length(flowpaths: &Vec<Flowpath>) -> Option<f64> {
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
