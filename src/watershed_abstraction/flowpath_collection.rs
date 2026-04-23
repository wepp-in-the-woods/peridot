use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Result, Write};
use std::path::Path;
use std::sync::Arc;

use arrow_array::{ArrayRef, Float64Array, Int32Array, RecordBatch, StringArray};
use arrow_schema::{DataType, Field, Schema};
use geojson::{Feature, FeatureCollection};
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use serde_json::to_string_pretty;

use crate::douglas_peucker::douglas_peucker;
use crate::flowpath::Flowpath;
use crate::raster::{px_to_wgs, Raster};
use crate::support::interpolate_slp;
use crate::watershed_abstraction::PATHS;

#[derive(Debug, Clone)]
pub struct FlowpathCollection {
    pub flowpaths: Vec<Flowpath>,
    pub subflows: Option<HashMap<i32, FlowpathCollection>>,
}

fn parquet_error(context: &str, err: impl std::fmt::Display) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Other, format!("{}: {}", context, err))
}

fn write_record_batch_parquet(
    path: &str,
    schema: Schema,
    columns: Vec<ArrayRef>,
) -> std::io::Result<()> {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            std::fs::create_dir_all(parent)?;
        }
    }

    let schema_ref = Arc::new(schema);
    let batch = RecordBatch::try_new(schema_ref.clone(), columns)
        .map_err(|err| parquet_error("failed to build record batch", err))?;

    let file = File::create(path)?;
    let props = WriterProperties::builder().build();
    let mut writer = ArrowWriter::try_new(file, schema_ref, Some(props))
        .map_err(|err| parquet_error("failed to create parquet writer", err))?;
    writer
        .write(&batch)
        .map_err(|err| parquet_error("failed to write parquet batch", err))?;
    writer
        .close()
        .map_err(|err| parquet_error("failed to close parquet writer", err))?;

    Ok(())
}

pub(crate) fn zonal_median_slope(fvslop: &Raster<f32>, indices: &Vec<usize>) -> f64 {
    let no_data = fvslop.no_data;
    let mut slopes: Vec<f64> = indices
        .iter()
        .filter_map(|&index| {
            let slope = fvslop.data[index];
            if !slope.is_finite() {
                return None;
            }
            if let Some(no_data_val) = no_data {
                if slope == no_data_val {
                    return None;
                }
            }
            Some(slope as f64)
        })
        .collect();

    assert!(
        !slopes.is_empty(),
        "zonal_median_slope found no finite values for indices len={}",
        indices.len()
    );

    slopes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = slopes.len();
    if n % 2 == 0 {
        (slopes[n / 2 - 1] + slopes[n / 2]) / 2.0
    } else {
        slopes[n / 2]
    }
}

fn median_length(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n % 2 == 0 {
        Some((values[n / 2 - 1] + values[n / 2]) / 2.0)
    } else {
        Some(values[n / 2])
    }
}

fn source_cells_from_indices(
    indices: &[usize],
    subwta: &Raster<i32>,
    flovec: &Raster<u8>,
) -> Vec<usize> {
    let Some(&first_index) = indices.first() else {
        return Vec::new();
    };

    let topaz_id = subwta.data[first_index];
    debug_assert!(
        indices.iter().all(|&index| subwta.data[index] == topaz_id),
        "source-cell detection expects one topaz_id per indices slice"
    );

    let mut has_upstream: HashSet<usize> = HashSet::with_capacity(indices.len() / 4 + 1);
    for &index in indices {
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
        if subwta.data[next_index] == topaz_id {
            has_upstream.insert(next_index);
        }
    }

    indices
        .iter()
        .copied()
        .filter(|index| !has_upstream.contains(index))
        .collect()
}

fn source_flowpath_median_length(flowpaths: &[Flowpath], source_cells: &[usize]) -> Option<f64> {
    if source_cells.is_empty() {
        return None;
    }
    let source_set: HashSet<usize> = source_cells.iter().copied().collect();

    let mut lengths: Vec<f64> = flowpaths
        .iter()
        .filter_map(|fp| {
            let head_index = fp.indices.first().copied()?;
            if source_set.contains(&head_index) {
                Some(fp.length)
            } else {
                None
            }
        })
        .collect();

    median_length(lengths.as_mut_slice())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SideLengthEstimateMode {
    AreaOverChannel,
    EdgeMedianCapped,
    AreaOverChannelNoEdge,
}

impl SideLengthEstimateMode {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::AreaOverChannel => "side_area_over_channel",
            Self::EdgeMedianCapped => "side_edge_median_capped",
            Self::AreaOverChannelNoEdge => "side_area_over_channel_no_edge",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct SideLengthSelection {
    pub length: f64,
    pub width: f64,
    pub length_area_over_channel: f64,
    pub edge_median: Option<f64>,
    pub mode: SideLengthEstimateMode,
}

pub(crate) fn select_side_hillslope_geometry(
    area: f64,
    channel_length: f64,
    edge_median: Option<f64>,
) -> SideLengthSelection {
    assert!(
        channel_length.is_finite() && channel_length > 0.0,
        "channel_length must be finite and > 0.0, got {}",
        channel_length
    );
    let length_area_over_channel = area / channel_length;

    let edge_median = edge_median.filter(|v| v.is_finite() && *v > 0.0);
    let (length, mode) = match edge_median {
        Some(edge) if edge < length_area_over_channel => {
            (edge, SideLengthEstimateMode::EdgeMedianCapped)
        }
        Some(_) => (
            length_area_over_channel,
            SideLengthEstimateMode::AreaOverChannel,
        ),
        None => (
            length_area_over_channel,
            SideLengthEstimateMode::AreaOverChannelNoEdge,
        ),
    };
    assert!(
        length.is_finite() && length > 0.0,
        "invalid selected side length {}",
        length
    );

    SideLengthSelection {
        length,
        width: area / length,
        length_area_over_channel,
        edge_median,
        mode,
    }
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
    pub fn weighted_slope_average_from_fps(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
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
            let mut num: f64 = 0.0; // to hold numerator value
            let mut kpsum: f64 = 0.0; // to hold k_p sum

            for fp in &self.flowpaths {
                if d_p > &(fp.length + 1e-6) {
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

            assert!(
                kpsum > 0.0,
                "kpsum is 0.0, d_p: {}, num_points: {}, longest_length: {}",
                d_p,
                num_points,
                longest_length
            );
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
    pub fn abstract_subcatchment(
        &self,
        subwta: &Raster<i32>,
        taspec: &Raster<f32>,
        fvslop: &Raster<f32>,
        flovec: &Raster<u8>,
        channels: &FlowpathCollection,
        indices: &Vec<usize>,
    ) -> Flowpath {
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
        let source_cells = source_cells_from_indices(indices, subwta, flovec);
        let edge_median = source_flowpath_median_length(&self.flowpaths, &source_cells);
        let (length, width, length_estimate_mode, length_area_over_channel, length_edge_median) =
            if topaz_id % 10 == 1 {
                // find length by taking the median length of flowpaths
                // that originate from source cells of the subcatchment.
                let length = edge_median.unwrap_or(0.0);
                (
                    length,
                    area / length,
                    "top_edge_median",
                    None,
                    if length.is_finite() && length > 0.0 {
                        Some(length)
                    } else {
                        None
                    },
                )
            } else {
                // Otherwise the width is inferred from area consistency after
                // selecting the side-length candidate.
                let selection =
                    select_side_hillslope_geometry(area, chn_summary.length, edge_median);
                (
                    selection.length,
                    selection.width,
                    selection.mode.as_str(),
                    Some(selection.length_area_over_channel),
                    selection.edge_median,
                )
            };

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
            let dx: f64 = distances[i + 1] - distances[i];
            let dy: f64 = w_slopes[i];
            let elevation: f64 = elevs[i] - (dx * dy);
            elevs.push(elevation);
        }

        let slope_scalar = zonal_median_slope(fvslop, indices);
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
            -1.0,
        )
        .with_length_estimate_metadata(
            length_estimate_mode,
            length_area_over_channel,
            length_edge_median,
        )
    }

    #[allow(dead_code)]
    pub fn abstract_subfieldcatchment(
        &self,
        intersection_subwta: &Raster<i32>,
        taspec: &Raster<f32>,
        flovec: &Raster<u8>,
        indices: &Vec<usize>,
    ) -> Flowpath {
        let cellsize: f64 = intersection_subwta.cellsize;
        let cellsize2: f64 = cellsize * cellsize;
        let fake_topaz_id: i32 = self.flowpaths[0].topaz_id;
        let area = indices.len() as f64 * cellsize2;

        let longest_fp = self.get_longest_fp();

        let length: f64;
        let width: f64;

        // find length by taking median of flowpaths that originate from source cells.
        let source_cells = source_cells_from_indices(indices, intersection_subwta, flovec);
        length = source_flowpath_median_length(&self.flowpaths, &source_cells).unwrap_or(0.0);
        width = area / length;

        let direction: f64 = longest_fp.direction;

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
            let dx: f64 = distances[i + 1] - distances[i];
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
            -1.0,
        )
    }

    #[allow(dead_code)]
    pub fn to_geojson_feature_collection(&self, raster: &Raster<i32>) -> FeatureCollection {
        let features: Vec<Feature> = self
            .flowpaths
            .iter()
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
        geojson.insert(
            String::from("type"),
            serde_json::Value::String(String::from("FeatureCollection")),
        );

        geojson.insert(
            String::from("features"),
            serde_json::Value::Array(
                feature_collection
                    .features
                    .into_iter()
                    .map(|f| serde_json::to_value(f).unwrap())
                    .collect(),
            ),
        );

        // Optionally add CRS if it exists in the raster
        if let Some(proj) = &raster.proj4 {
            let mut crs = serde_json::Map::new();
            crs.insert(
                String::from("type"),
                serde_json::Value::String(String::from("name")),
            );
            crs.insert(
                String::from("properties"),
                serde_json::Value::Object(serde_json::Map::from_iter(std::iter::once((
                    String::from("name"),
                    serde_json::Value::String(proj.clone()),
                )))),
            );
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
            subflows_map
                .par_iter()
                .for_each(|(topaz_id, subflowpath_collection)| {
                    let path = format!("{}/fps_{}.slps", out_dir, topaz_id);
                    subflowpath_collection
                        .write_fp_slps(&path, max_points)
                        .unwrap();
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
                println!(
                    "douglas_peuker {} -> {}",
                    fp.distances_norm.len(),
                    simplified.0.len()
                );
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
            let defs: Vec<String> = d
                .iter()
                .zip(s.iter())
                .map(|(&dist, &slope)| format!("{:.4}, {:.4}", dist, slope))
                .collect();

            let slp = format!(
                "{:.4} {:.1} {:.1} {}\n{} {:.1}\n{} \n",
                fp.aspect,
                fp.width,
                fp.elevation,
                fp.order,
                npts,
                fp.length,
                defs.join(" ")
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
    pub fn write_slps(
        &self,
        out_dir: &str,
        max_points: usize,
        clip_hillslopes: bool,
        clip_hillslope_length: f64,
    ) -> std::io::Result<()> {
        let results: Vec<std::io::Result<()>> = self
            .flowpaths
            .par_iter()
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
            if let Some((field_id, topaz_id)) =
                Self::resolve_field_lookup(fp.topaz_id, fake_topaz_id_lookup)
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
    pub fn subflow_row_count(&self) -> usize {
        self.subflows
            .as_ref()
            .map(|subflows_map| {
                subflows_map
                    .values()
                    .map(|collection| collection.flowpaths.len())
                    .sum()
            })
            .unwrap_or(0)
    }

    #[allow(dead_code)]
    pub fn write_chn_metadata_to_parquet(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
    ) -> std::io::Result<()> {
        let mut topaz_ids: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut slope_scalars: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut lengths: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut widths: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut directions: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut orders: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut aspects: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut areas: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut elevations: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_px: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_py: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_lon: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_lat: Vec<f64> = Vec::with_capacity(self.flowpaths.len());

        for fp in &self.flowpaths {
            let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);
            topaz_ids.push(fp.topaz_id);
            slope_scalars.push(fp.slope_scalar);
            lengths.push(fp.length);
            widths.push(fp.width);
            directions.push(fp.direction);
            orders.push(fp.order);
            aspects.push(fp.aspect);
            areas.push(fp.area_m2());
            elevations.push(fp.elevation);
            centroid_px.push(fp.centroid_px.0);
            centroid_py.push(fp.centroid_px.1);
            centroid_lon.push(lon);
            centroid_lat.push(lat);
        }

        let schema = Schema::new(vec![
            Field::new("topaz_id", DataType::Int32, false),
            Field::new("slope_scalar", DataType::Float64, false),
            Field::new("length", DataType::Float64, false),
            Field::new("width", DataType::Float64, false),
            Field::new("direction", DataType::Float64, false),
            Field::new("order", DataType::Int32, false),
            Field::new("aspect", DataType::Float64, false),
            Field::new("area", DataType::Float64, false),
            Field::new("elevation", DataType::Float64, false),
            Field::new("centroid_px", DataType::Int32, false),
            Field::new("centroid_py", DataType::Int32, false),
            Field::new("centroid_lon", DataType::Float64, false),
            Field::new("centroid_lat", DataType::Float64, false),
        ]);

        write_record_batch_parquet(
            path,
            schema,
            vec![
                Arc::new(Int32Array::from(topaz_ids)) as ArrayRef,
                Arc::new(Float64Array::from(slope_scalars)),
                Arc::new(Float64Array::from(lengths)),
                Arc::new(Float64Array::from(widths)),
                Arc::new(Float64Array::from(directions)),
                Arc::new(Int32Array::from(orders)),
                Arc::new(Float64Array::from(aspects)),
                Arc::new(Float64Array::from(areas)),
                Arc::new(Float64Array::from(elevations)),
                Arc::new(Int32Array::from(centroid_px)),
                Arc::new(Int32Array::from(centroid_py)),
                Arc::new(Float64Array::from(centroid_lon)),
                Arc::new(Float64Array::from(centroid_lat)),
            ],
        )
    }

    #[allow(dead_code)]
    pub fn write_metadata_to_parquet(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
    ) -> std::io::Result<()> {
        let mut topaz_ids: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut slope_scalars: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut lengths: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut widths: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut directions: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut aspects: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut length_estimate_modes: Vec<String> = Vec::with_capacity(self.flowpaths.len());
        let mut lengths_area_over_channel: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut lengths_edge_median: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut areas: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut elevations: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_px: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_py: Vec<i32> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_lon: Vec<f64> = Vec::with_capacity(self.flowpaths.len());
        let mut centroid_lat: Vec<f64> = Vec::with_capacity(self.flowpaths.len());

        for fp in &self.flowpaths {
            let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);
            topaz_ids.push(fp.topaz_id);
            slope_scalars.push(fp.slope_scalar);
            lengths.push(fp.length);
            widths.push(fp.width);
            directions.push(fp.direction);
            aspects.push(fp.aspect);
            length_estimate_modes.push(fp.length_estimate_mode.clone());
            lengths_area_over_channel.push(fp.length_area_over_channel);
            lengths_edge_median.push(fp.length_edge_median);
            areas.push(fp.area_m2() as i32);
            elevations.push(fp.elevation);
            centroid_px.push(fp.centroid_px.0);
            centroid_py.push(fp.centroid_px.1);
            centroid_lon.push(lon);
            centroid_lat.push(lat);
        }

        let schema = Schema::new(vec![
            Field::new("topaz_id", DataType::Int32, false),
            Field::new("slope_scalar", DataType::Float64, false),
            Field::new("length", DataType::Float64, false),
            Field::new("width", DataType::Float64, false),
            Field::new("direction", DataType::Float64, false),
            Field::new("aspect", DataType::Float64, false),
            Field::new("length_estimate_mode", DataType::Utf8, false),
            Field::new("length_area_over_channel", DataType::Float64, false),
            Field::new("length_edge_median", DataType::Float64, false),
            Field::new("area", DataType::Int32, false),
            Field::new("elevation", DataType::Float64, false),
            Field::new("centroid_px", DataType::Int32, false),
            Field::new("centroid_py", DataType::Int32, false),
            Field::new("centroid_lon", DataType::Float64, false),
            Field::new("centroid_lat", DataType::Float64, false),
        ]);

        write_record_batch_parquet(
            path,
            schema,
            vec![
                Arc::new(Int32Array::from(topaz_ids)) as ArrayRef,
                Arc::new(Float64Array::from(slope_scalars)),
                Arc::new(Float64Array::from(lengths)),
                Arc::new(Float64Array::from(widths)),
                Arc::new(Float64Array::from(directions)),
                Arc::new(Float64Array::from(aspects)),
                Arc::new(StringArray::from(length_estimate_modes)),
                Arc::new(Float64Array::from(lengths_area_over_channel)),
                Arc::new(Float64Array::from(lengths_edge_median)),
                Arc::new(Int32Array::from(areas)),
                Arc::new(Float64Array::from(elevations)),
                Arc::new(Int32Array::from(centroid_px)),
                Arc::new(Int32Array::from(centroid_py)),
                Arc::new(Float64Array::from(centroid_lon)),
                Arc::new(Float64Array::from(centroid_lat)),
            ],
        )
    }

    #[allow(dead_code)]
    pub fn write_subflows_metadata_to_parquet(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
    ) -> std::io::Result<()> {
        let mut topaz_ids: Vec<i32> = Vec::new();
        let mut fp_ids: Vec<i32> = Vec::new();
        let mut slope_scalars: Vec<f64> = Vec::new();
        let mut lengths: Vec<f64> = Vec::new();
        let mut widths: Vec<f64> = Vec::new();
        let mut directions: Vec<f64> = Vec::new();
        let mut aspects: Vec<f64> = Vec::new();
        let mut areas: Vec<f64> = Vec::new();
        let mut elevations: Vec<f64> = Vec::new();
        let mut orders: Vec<i32> = Vec::new();
        let mut centroid_px: Vec<i32> = Vec::new();
        let mut centroid_py: Vec<i32> = Vec::new();
        let mut centroid_lon: Vec<f64> = Vec::new();
        let mut centroid_lat: Vec<f64> = Vec::new();

        if let Some(subflows_map) = &self.subflows {
            for (topaz_id, subflow_collection) in subflows_map {
                for fp in &subflow_collection.flowpaths {
                    let (lon, lat) = px_to_wgs(wgs_transform, fp.centroid_px.0, fp.centroid_px.1);
                    topaz_ids.push(*topaz_id);
                    fp_ids.push(fp.fp_id);
                    slope_scalars.push(fp.slope_scalar);
                    lengths.push(fp.length);
                    widths.push(fp.width);
                    directions.push(fp.direction);
                    aspects.push(fp.aspect);
                    areas.push(fp.area_m2());
                    elevations.push(fp.elevation);
                    orders.push(fp.order);
                    centroid_px.push(fp.centroid_px.0);
                    centroid_py.push(fp.centroid_px.1);
                    centroid_lon.push(lon);
                    centroid_lat.push(lat);
                }
            }
        } else {
            panic!("Unexpected None in subflows!");
        }

        let schema = Schema::new(vec![
            Field::new("topaz_id", DataType::Int32, false),
            Field::new("fp_id", DataType::Int32, false),
            Field::new("slope_scalar", DataType::Float64, false),
            Field::new("length", DataType::Float64, false),
            Field::new("width", DataType::Float64, false),
            Field::new("direction", DataType::Float64, false),
            Field::new("aspect", DataType::Float64, false),
            Field::new("area", DataType::Float64, false),
            Field::new("elevation", DataType::Float64, false),
            Field::new("order", DataType::Int32, false),
            Field::new("centroid_px", DataType::Int32, false),
            Field::new("centroid_py", DataType::Int32, false),
            Field::new("centroid_lon", DataType::Float64, false),
            Field::new("centroid_lat", DataType::Float64, false),
        ]);

        write_record_batch_parquet(
            path,
            schema,
            vec![
                Arc::new(Int32Array::from(topaz_ids)) as ArrayRef,
                Arc::new(Int32Array::from(fp_ids)),
                Arc::new(Float64Array::from(slope_scalars)),
                Arc::new(Float64Array::from(lengths)),
                Arc::new(Float64Array::from(widths)),
                Arc::new(Float64Array::from(directions)),
                Arc::new(Float64Array::from(aspects)),
                Arc::new(Float64Array::from(areas)),
                Arc::new(Float64Array::from(elevations)),
                Arc::new(Int32Array::from(orders)),
                Arc::new(Int32Array::from(centroid_px)),
                Arc::new(Int32Array::from(centroid_py)),
                Arc::new(Float64Array::from(centroid_lon)),
                Arc::new(Float64Array::from(centroid_lat)),
            ],
        )
    }

    #[allow(dead_code)]
    pub fn write_chn_metadata_to_csv(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
    ) -> std::io::Result<()> {
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
    pub fn write_metadata_to_csv(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
    ) -> std::io::Result<()> {
        let file = File::create(path).unwrap();
        let mut writer = csv::Writer::from_writer(file);

        let headers: Vec<String> = vec![
            String::from("topaz_id"),
            String::from("slope_scalar"),
            String::from("length"),
            String::from("width"),
            String::from("direction"),
            String::from("aspect"),
            String::from("length_estimate_mode"),
            String::from("length_area_over_channel"),
            String::from("length_edge_median"),
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
                fp.length_estimate_mode.to_string(),
                fp.length_area_over_channel.to_string(),
                fp.length_edge_median.to_string(),
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
            if let Some((field_id, topaz_id)) =
                Self::resolve_field_lookup(fp.topaz_id, fake_topaz_id_lookup)
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
    pub fn write_subflows_metadata_to_csv(
        &self,
        path: &str,
        wgs_transform: &[f64; 4],
    ) -> std::io::Result<()> {
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
    #[deprecated(
        note = "Use get_edge_flowpaths2 with subwta+flovec for faster, source-cell detection."
    )]
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
                if i == j {
                    continue;
                } // Skip self-comparison

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
    pub fn get_edge_flowpaths2(&self, subwta: &Raster<i32>, flovec: &Raster<u8>) -> Vec<Flowpath> {
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

#[cfg(test)]
mod tests {
    use crate::raster::{MapType, Raster};

    use super::{select_side_hillslope_geometry, zonal_median_slope, SideLengthEstimateMode};

    fn mock_fvslop(values: Vec<f32>, no_data: Option<f32>) -> Raster<f32> {
        Raster::new(
            values.len(),
            1,
            1.0,
            values,
            no_data,
            [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
            None,
            String::from(""),
            String::from("fvslop"),
            MapType::FVSLOP,
        )
    }

    #[test]
    fn zonal_median_slope_even_count() {
        let fvslop = mock_fvslop(vec![0.5, 0.1, 0.3, 0.9], None);
        let indices = vec![0, 1, 2, 3];
        let median = zonal_median_slope(&fvslop, &indices);
        assert!((median - 0.4).abs() < 1e-6);
    }

    #[test]
    fn zonal_median_slope_ignores_no_data_and_nan() {
        let fvslop = mock_fvslop(vec![0.2, -9999.0, f32::NAN, 0.6, 0.4], Some(-9999.0));
        let indices = vec![0, 1, 2, 3, 4];
        let median = zonal_median_slope(&fvslop, &indices);
        assert!((median - 0.4).abs() < 1e-6);
    }

    #[test]
    fn side_length_selection_caps_to_edge_median_and_preserves_area() {
        let area = 10_000.0;
        let selection = select_side_hillslope_geometry(area, 50.0, Some(120.0));
        assert_eq!(selection.mode, SideLengthEstimateMode::EdgeMedianCapped);
        assert!((selection.length_area_over_channel - 200.0).abs() < 1e-9);
        assert!((selection.length - 120.0).abs() < 1e-9);
        assert!((selection.width * selection.length - area).abs() < 1e-6);
    }

    #[test]
    fn side_length_selection_uses_area_over_channel_when_edge_is_longer() {
        let area = 10_000.0;
        let selection = select_side_hillslope_geometry(area, 50.0, Some(260.0));
        assert_eq!(selection.mode, SideLengthEstimateMode::AreaOverChannel);
        assert!((selection.length_area_over_channel - 200.0).abs() < 1e-9);
        assert!((selection.length - 200.0).abs() < 1e-9);
        assert!((selection.width - 50.0).abs() < 1e-9);
        assert!((selection.width * selection.length - area).abs() < 1e-6);
    }

    #[test]
    fn side_length_selection_falls_back_when_edge_is_missing() {
        let area = 10_000.0;
        let selection = select_side_hillslope_geometry(area, 50.0, None);
        assert_eq!(
            selection.mode,
            SideLengthEstimateMode::AreaOverChannelNoEdge
        );
        assert!((selection.length - 200.0).abs() < 1e-9);
        assert!((selection.width - 50.0).abs() < 1e-9);
    }
}
