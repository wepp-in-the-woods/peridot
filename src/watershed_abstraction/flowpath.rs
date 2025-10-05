use std::fs::File;
use std::io::Write;

use interp::interp;

use geojson::{Feature, Geometry};
use geojson::Value::LineString;

use crate::support::{douglas_peucker, interpolate_slp};
use crate::raster::Raster;

#[derive(Debug, Clone)]
pub struct Flowpath {
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

impl Default for Flowpath {
    fn default() -> Self {
        Flowpath {
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

impl Flowpath {

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

        Flowpath {
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
        let file = File::create(path)?;

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
        let mut properties: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();
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
