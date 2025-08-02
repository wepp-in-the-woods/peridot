
use std::process::Command;

use std::io;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::raster::Raster;

use serde_json::{Value, json};
use std::fs;

pub fn combine_geojson_files(filenames: &Vec<String>, output_filename: &str, epsg: &str, sort_key: Option<&str>)
    -> Result<(), io::Error> {
    let mut combined_features: Vec<Value> = Vec::new();

    // Read each GeoJSON file and combine features
    for filename in filenames {
        let file_contents = fs::read_to_string(filename)?;
        let geojson: Value = serde_json::from_str(&file_contents)?;

        if let Some(features) = geojson["features"].as_array() {
            combined_features.extend(features.clone());
        }
    }

    // Sort the features by the specified sort_key (if provided) in ascending order
    if let Some(key) = sort_key {
        combined_features.sort_by(|a, b| {
            let a_value = a["properties"][key].as_f64().unwrap_or(f64::MAX);
            let b_value = b["properties"][key].as_f64().unwrap_or(f64::MAX);
            a_value.partial_cmp(&b_value).unwrap()
        });
    }

    // Create the CRS JSON object
    let crs = json!({
        "type": "name",
        "properties": {
            "name": format!("urn:ogc:def:crs:EPSG::{}", epsg)
        }
    });


    // Create a new GeoJSON FeatureCollection object
    let combined_geojson = json!({
        "crs": crs,
        "type": "FeatureCollection",
        "features": combined_features
    });

    // Serialize and write to output file
    let serialized = serde_json::to_string_pretty(&combined_geojson)?;
    fs::write(output_filename, serialized)?;

    Ok(())
}

pub fn polygonize_raster(raster_src_fn: &str, geojson_dst_fn: &str, properties: &Value) -> io::Result<()> {
    let output = Command::new("python3")
        .arg("/usr/bin/gdal_polygonize.py")
        .arg(raster_src_fn)
        .arg("-mask")
        .arg(raster_src_fn)
        .arg("-f")
        .arg("GeoJSON")
        .arg(geojson_dst_fn)
        .output()?;

    if !output.status.success() {
        eprintln!("gdal_translate failed: {}", String::from_utf8_lossy(&output.stderr));
        return Err(io::Error::new(io::ErrorKind::Other, "gdal_translate command failed"));
    }

     // Read the existing GeoJSON file
     let file_contents = fs::read_to_string(geojson_dst_fn)?;
     let mut geojson: Value = serde_json::from_str(&file_contents)?;
 
     // Modify the GeoJSON
     if let Some(features) = geojson["features"].as_array_mut() {
         for feature in features {
             if let Some(feature_properties) = feature["properties"].as_object_mut() {                
                // Remove the DN field
                feature_properties.remove("DN");

                 // Merge the new properties
                 for (key, value) in properties.as_object().unwrap() {
                     feature_properties.insert(key.clone(), value.clone());
                 }
             }
         }
     }
 
     // Serialize the modified GeoJSON object back to a string
     let modified_geojson = serde_json::to_string_pretty(&geojson)?;
 
     // Write back to the file
     let mut file = fs::File::create(geojson_dst_fn)?;
     file.write_all(modified_geojson.as_bytes())?;
 
     Ok(())
}

pub fn remap_whitebox_d8_to_topaz(flovec: &Raster<i32>) -> Raster<i32> {
    let mut remapped_flovec = flovec.empty_clone();

    for i in 0..flovec.data.len() {
        let flow_dir = flovec.data[i];
        let new_flow_dir = match flow_dir {
            1 => 3,    // East
            2 => 6,    // Northeast
            4 => 9,    // North
            8 => 8,    // Northwest
            16 => 7,   // West
            32 => 4,   // Southwest
            64 => 1,   // South
            128 => 2,  // Southeast
            _ => 0,    // No flow or undefined
        };
        remapped_flovec.data[i] = new_flow_dir;
    }

    remapped_flovec
}

// 1     2     3
//   \   |   /
//    \  | /
// 4 <-- o --> 6
//     / | \
//   /   |   \
// 7     8     9

//64	128	1
//32	0	2
//16	8	4


pub fn rescale_raster(src_fn: &str, dst_fn: &str, resolution: f64) -> io::Result<()> {
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
//        println!("gdal_translate executed successfully.");
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

pub fn smooth_raster(src_fn: &str, dst_fn: &str, sigma: f64) -> io::Result<()> {

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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

pub fn fill_depressions_raster(src_fn: &str, dst_fn: &str) -> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}


pub fn breach_depressions_raster(src_fn: &str, dst_fn: &str) -> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}


pub fn slope_raster(dem_fn: &str, slope_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

pub fn aspect_raster(dem_fn: &str, aspect_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

pub fn d8_flow_direction_raster(dem_fn: &str, flow_dir_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

pub fn dinf_raster(dem_fn: &str,  dinf_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

pub fn d8_flowaccum_raster(d8_fn: &str,  d8_accum_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}


pub fn dinf_flowaccum_raster(dinf_fn: &str,  dinf_accum_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

// threshold represents the minimum area (area is used here as a surrogate for discharge) 
// required to initiate and maintain a channel. Smaller threshold values result in more 
// extensive stream networks
pub fn stream_raster(dinf_accum_fn: &str, stream_fn: &str, threshold: f64) -> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
    }
}

pub fn valley_raster(dem_fn: &str, valley_fn: &str)-> io::Result<()> {
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
//        println!("whitebox_tools executed successfully.");
        Ok(())
    } else {
        eprintln!("whitebox_tools failed: {}", String::from_utf8_lossy(&output.stderr));
        //Err(io::Error::new(io::ErrorKind::Other, "whitebox_tools command failed"))
        Ok(())
    }
}
