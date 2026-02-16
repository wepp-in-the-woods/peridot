use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::io::Result;
use std::mem::size_of;
use std::path::Path;
use std::time::Instant;

use log::{info, warn};
use rayon::prelude::*;

use crate::d8_wbt_to_topaz::remap_whitebox_d8_to_topaz_in_place;
use crate::flowpath::Flowpath;
use crate::flowpath_collection::FlowpathCollection;
use crate::netw::write_network;
use crate::raster::Raster;
use crate::watershed_abstraction::{abstract_subcatchments, walk_channels, walk_flowpath, PATHS};
use crate::wbt_netw::read_wbt_netw_tab;

/// Find raster file, checking .tif first then .vrt as fallback
fn find_raster_path(base_path: &str) -> String {
    let tif_path = format!("{}.tif", base_path);
    if Path::new(&tif_path).exists() {
        return tif_path;
    }
    let vrt_path = format!("{}.vrt", base_path);
    if Path::new(&vrt_path).exists() {
        return vrt_path;
    }
    // Return .tif path to get original error message if neither exists
    tif_path
}

fn log_flowpath_collection_stats(label: &str, collection: &FlowpathCollection) {
    let flowpaths = &collection.flowpaths;
    let indices_total: usize = flowpaths.iter().map(|fp| fp.indices.len()).sum();
    let points_total: usize = flowpaths.iter().map(|fp| fp.distances_norm.len()).sum();
    let bytes = (indices_total * size_of::<usize>()) + (points_total * size_of::<f64>() * 5);
    info!(
        "{}: flowpaths={} indices_total={} points_total={} approx_bytes={} (~{:.2} MiB)",
        label,
        flowpaths.len(),
        indices_total,
        points_total,
        bytes,
        bytes as f64 / (1024.0 * 1024.0)
    );

    if let Some(subflows) = &collection.subflows {
        let mut subflow_paths = 0usize;
        let mut subflow_indices = 0usize;
        let mut subflow_points = 0usize;
        for sub in subflows.values() {
            subflow_paths += sub.flowpaths.len();
            subflow_indices += sub
                .flowpaths
                .iter()
                .map(|fp| fp.indices.len())
                .sum::<usize>();
            subflow_points += sub
                .flowpaths
                .iter()
                .map(|fp| fp.distances_norm.len())
                .sum::<usize>();
        }
        let subflow_bytes =
            (subflow_indices * size_of::<usize>()) + (subflow_points * size_of::<f64>() * 5);
        info!(
            "{} subflows: collections={} flowpaths={} indices_total={} points_total={} approx_bytes={} (~{:.2} MiB)",
            label,
            subflows.len(),
            subflow_paths,
            subflow_indices,
            subflow_points,
            subflow_bytes,
            subflow_bytes as f64 / (1024.0 * 1024.0)
        );
    }
}
#[allow(dead_code)]
pub fn wbt_abstract_watershed(
    wd: &str,
    max_points: usize,
    clip_hillslopes: bool,
    clip_hillslope_length: f64,
    bieger2015_widths: bool,
    write_flowpaths: bool,
    representative_flowpath: bool,
) -> std::io::Result<()> {
    env::set_current_dir(&wd).unwrap();
    let write_flowpaths = if representative_flowpath {
        false
    } else {
        write_flowpaths
    };
    info!(
        "wbt_abstract_watershed: wd={}, max_points={}, clip_hillslopes={}, clip_hillslope_length={}, bieger2015_widths={}, write_flowpaths={}, representative_flowpath={}",
        wd,
        max_points,
        clip_hillslopes,
        clip_hillslope_length,
        bieger2015_widths,
        write_flowpaths,
        representative_flowpath
    );

    let watershed_path = Path::new("watershed");

    if watershed_path.exists() {
        let _ = fs::remove_dir_all(&watershed_path).unwrap();
    }

    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/hillslopes"));
    let _ = std::fs::create_dir_all(Path::new("watershed").join("slope_files/flowpaths"));

    let t0 = Instant::now();
    let subwta: Raster<i32> = Raster::<i32>::read(&find_raster_path("dem/wbt/subwta")).unwrap(); // hillslope with topaz_ids file, channels end with 4 (e.g 24, 34, 44), subcatchments end with 1, 2, 3. It starts at 22
    info!(
        "loaded subwta: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        subwta.width,
        subwta.height,
        subwta.data.len(),
        subwta.data.len() * size_of::<i32>(),
        (subwta.data.len() * size_of::<i32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let relief: Raster<f32> = Raster::<f32>::read(&find_raster_path("dem/wbt/relief")).unwrap(); // dem
    info!(
        "loaded relief: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        relief.width,
        relief.height,
        relief.data.len(),
        relief.data.len() * size_of::<f32>(),
        (relief.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let mut flovec: Raster<u8> = Raster::<u8>::read(&find_raster_path("dem/wbt/flovec")).unwrap(); // d8 flowvec
    remap_whitebox_d8_to_topaz_in_place(&mut flovec);
    info!(
        "loaded flovec: {}x{} cells={} bytes={} (~{:.2} MiB) and remapped in {:.2}s",
        flovec.width,
        flovec.height,
        flovec.data.len(),
        flovec.data.len() * size_of::<u8>(),
        (flovec.data.len() * size_of::<u8>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let fvslop: Raster<f32> = Raster::<f32>::read(&find_raster_path("dem/wbt/fvslop")).unwrap(); // slope
    info!(
        "loaded fvslop: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        fvslop.width,
        fvslop.height,
        fvslop.data.len(),
        fvslop.data.len() * size_of::<f32>(),
        (fvslop.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
    let taspec: Raster<f32> = Raster::<f32>::read(&find_raster_path("dem/wbt/taspec")).unwrap(); // aspect
    info!(
        "loaded taspec: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
        taspec.width,
        taspec.height,
        taspec.data.len(),
        taspec.data.len() * size_of::<f32>(),
        (taspec.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    let discha: Option<Raster<f32>> = if representative_flowpath {
        let t0 = Instant::now();
        let discha_raster: Raster<f32> =
            Raster::<f32>::read(&find_raster_path("dem/wbt/discha")).unwrap();
        info!(
            "loaded discha: {}x{} cells={} bytes={} (~{:.2} MiB) in {:.2}s",
            discha_raster.width,
            discha_raster.height,
            discha_raster.data.len(),
            discha_raster.data.len() * size_of::<f32>(),
            (discha_raster.data.len() * size_of::<f32>()) as f64 / (1024.0 * 1024.0),
            t0.elapsed().as_secs_f64()
        );
        Some(discha_raster)
    } else {
        None
    };

    let t0 = Instant::now();
    let subwta_indices = subwta.indices_map();
    let channel_ids = subwta_indices.keys().filter(|&&id| id % 10 == 4).count();
    let hillslope_ids = subwta_indices.keys().filter(|&&id| id % 10 != 4).count();
    let indices_total: usize = subwta_indices.values().map(|v| v.len()).sum();
    info!(
        "built subwta indices map: ids={} (channels={}, hillslopes={}) total_cells={} (~{:.2} MiB) in {:.2}s",
        subwta_indices.len(),
        channel_ids,
        hillslope_ids,
        indices_total,
        (indices_total * size_of::<usize>()) as f64 / (1024.0 * 1024.0),
        t0.elapsed().as_secs_f64()
    );

    // this is a ASCII tabular report with channel node connection informaiton
    let t0 = Instant::now();
    let (netw, network) = read_wbt_netw_tab("dem/wbt/netw.tsv").unwrap();
    let _ = write_network("watershed/network.txt", &network);
    info!(
        "read netw.tsv: links={} nodes={} in {:.2}s",
        netw.len(),
        network.len(),
        t0.elapsed().as_secs_f64()
    );

    let t0 = Instant::now();
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
    info!(
        "walked channels: flowpaths={} in {:.2}s",
        channels.flowpaths.len(),
        t0.elapsed().as_secs_f64()
    );
    log_flowpath_collection_stats("channels", &channels);

    let t0 = Instant::now();
    let hillslopes: FlowpathCollection = if representative_flowpath {
        let discha = discha
            .as_ref()
            .expect("discha raster required for representative flowpath mode");
        abstract_subcatchments_representative(
            &subwta,
            &subwta_indices,
            &relief,
            &flovec,
            &fvslop,
            &taspec,
            discha,
            &channels,
        )
    } else {
        abstract_subcatchments(
            &subwta,
            &subwta_indices,
            &relief,
            &flovec,
            &fvslop,
            &taspec,
            &channels,
            write_flowpaths,
        )
    };
    let subflow_count = hillslopes.subflows.as_ref().map(|m| m.len()).unwrap_or(0);
    info!(
        "abstracted subcatchments: hillslopes={} subflows={} in {:.2}s",
        hillslopes.flowpaths.len(),
        subflow_count,
        t0.elapsed().as_secs_f64()
    );
    log_flowpath_collection_stats("hillslopes", &hillslopes);

    let mut tasks: Vec<Box<dyn FnOnce() -> Result<()> + Send>> = vec![
        Box::new(|| {
            info!("writing channels.slp");
            let result =
                channels.write_channel_slp("watershed/slope_files/channels.slp", max_points);
            info!("wrote channels.slp");
            result
        }),
        Box::new(|| {
            info!("writing channels.csv");
            let result =
                channels.write_chn_metadata_to_csv("watershed/channels.csv", &subwta.wgs_transform);
            info!("wrote channels.csv");
            result
        }),
        Box::new(|| {
            info!("writing hillslope slps");
            let result = hillslopes.write_slps(
                "watershed/slope_files/hillslopes/",
                max_points,
                clip_hillslopes,
                clip_hillslope_length,
            );
            info!("wrote hillslope slps");
            result
        }),
        Box::new(|| {
            info!("writing hillslopes.csv");
            let result =
                hillslopes.write_metadata_to_csv("watershed/hillslopes.csv", &subwta.wgs_transform);
            info!("wrote hillslopes.csv");
            result
        }),
        Box::new(|| {
            info!("writing channels.geojson");
            let result = channels.write_geojson(&subwta, "watershed/channels.geojson");
            info!("wrote channels.geojson");
            result
        }),
    ];

    if write_flowpaths {
        tasks.push(Box::new(|| {
            info!("writing flowpaths.csv");
            let result = hillslopes
                .write_subflows_metadata_to_csv("watershed/flowpaths.csv", &subwta.wgs_transform);
            info!("wrote flowpaths.csv");
            result
        }));
        tasks.push(Box::new(|| {
            info!("writing flowpath slps");
            let result =
                hillslopes.write_subflow_slps("watershed/slope_files/flowpaths/", max_points);
            info!("wrote flowpath slps");
            result
        }));
    }

    // Execute tasks in parallel
    info!("starting output tasks: {}", tasks.len());
    tasks
        .into_par_iter()
        .map(|f| f())
        .collect::<Result<Vec<_>>>()?;
    info!("completed output tasks");

    Ok(())
}

#[derive(Clone, Copy)]
struct SeedCandidate {
    index: usize,
    discha: f32,
    relief: f32,
}

fn build_source_cells(
    indices: &Vec<usize>,
    subwta: &Raster<i32>,
    flovec: &Raster<u8>,
) -> Vec<usize> {
    let mut mask_indices: HashSet<usize> = HashSet::with_capacity(indices.len());
    for &index in indices {
        mask_indices.insert(index);
    }

    let mut has_upstream: HashSet<usize> = HashSet::new();
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
        if mask_indices.contains(&next_index) {
            has_upstream.insert(next_index);
        }
    }

    indices
        .iter()
        .copied()
        .filter(|index| !has_upstream.contains(index))
        .collect()
}

fn build_seed_candidates(
    source_cells: &[usize],
    discha: &Raster<f32>,
    relief: &Raster<f32>,
) -> Vec<SeedCandidate> {
    let mut candidates = Vec::new();
    let no_data = discha.no_data;
    for &index in source_cells {
        let d = discha.data[index];
        if !d.is_finite() || d <= 0.0 {
            continue;
        }
        if let Some(no_data_val) = no_data {
            if d == no_data_val {
                continue;
            }
        }
        let r = relief.data[index];
        candidates.push(SeedCandidate {
            index,
            discha: d,
            relief: r,
        });
    }
    candidates
}

fn flowpath_reaches_channel(
    flowpath: &Flowpath,
    subwta: &Raster<i32>,
    flovec: &Raster<u8>,
    topaz_id: i32,
) -> bool {
    let tail_index = match flowpath.indices.last() {
        Some(&index) => index,
        None => return false,
    };
    let flow_dir = flovec.data[tail_index] as i32;
    if flow_dir == 0 {
        return false;
    }
    let Some((dx, dy)) = PATHS.get(&flow_dir) else {
        return false;
    };
    let (x, y) = subwta.index_to_xy(tail_index);
    let next_x = x as isize + dx;
    let next_y = y as isize + dy;
    if next_x < 0
        || next_y < 0
        || next_x >= subwta.width as isize
        || next_y >= subwta.height as isize
    {
        return false;
    }
    let next_index = subwta.xy_to_index(next_x as usize, next_y as usize);
    let next_topaz_id = subwta.data[next_index];
    next_topaz_id % 10 == 4 && next_topaz_id != topaz_id
}

fn select_representative_flowpath(
    topaz_id: i32,
    indices: &Vec<usize>,
    subwta: &Raster<i32>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
    discha: &Raster<f32>,
) -> Flowpath {
    let source_cells = build_source_cells(indices, subwta, flovec);
    let mut candidates = build_seed_candidates(&source_cells, discha, relief);

    if candidates.is_empty() {
        warn!(
            "representative_flowpath: no discha candidates for topaz_id={}, falling back to first index",
            topaz_id
        );
        let seed_index = *indices.first().expect("missing hillslope indices");
        return walk_flowpath(
            topaz_id, seed_index, subwta, relief, flovec, fvslop, taspec, -1,
        );
    }

    candidates.sort_by(|a, b| {
        a.discha
            .partial_cmp(&b.discha)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.relief
                    .partial_cmp(&a.relief)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.index.cmp(&b.index))
    });

    let mid = candidates.len() / 2;
    let mut fallback_fp: Option<Flowpath> = None;
    let mut order: Vec<usize> = Vec::with_capacity(candidates.len());
    order.push(mid);
    for offset in 1..candidates.len() {
        if mid >= offset {
            order.push(mid - offset);
        }
        let right = mid + offset;
        if right < candidates.len() {
            order.push(right);
        }
    }

    for candidate_index in order {
        let seed_index = candidates[candidate_index].index;
        let fp = walk_flowpath(
            topaz_id, seed_index, subwta, relief, flovec, fvslop, taspec, -1,
        );
        if fallback_fp.is_none() {
            fallback_fp = Some(fp.clone());
        }
        if flowpath_reaches_channel(&fp, subwta, flovec, topaz_id) {
            return fp;
        }
    }

    warn!(
        "representative_flowpath: no candidate reached channel for topaz_id={}, using median fallback",
        topaz_id
    );
    fallback_fp.expect("missing fallback representative flowpath")
}

fn build_representative_hillslope(
    representative_fp: &Flowpath,
    subwta: &Raster<i32>,
    taspec: &Raster<f32>,
    channels: &FlowpathCollection,
    indices: &Vec<usize>,
) -> Flowpath {
    let cellsize: f64 = subwta.cellsize;
    let cellsize2: f64 = cellsize * cellsize;
    let topaz_id: i32 = representative_fp.topaz_id;
    let area = indices.len() as f64 * cellsize2;

    let chn_id: i32 = ((topaz_id as f64 / 10.0).floor() * 10.0) as i32 + 4;
    let chn_summary: &Flowpath = &channels.get_fp_by_topaz_id(chn_id).unwrap();

    let length: f64;
    let width: f64;
    if topaz_id % 10 == 1 {
        length = representative_fp.length;
        width = area / length;
    } else {
        width = chn_summary.length;
        length = area / width;
    }

    let mut direction: f64 = chn_summary.direction;
    match topaz_id % 10 {
        2 => direction += 90.0,
        3 => direction -= 90.0,
        _ => (),
    }

    let aspect: f64 = taspec.determine_aspect(indices);
    let centroid_px = subwta.centroid_of(indices);
    let elevation = representative_fp.elevs[0];
    let elev_drop =
        representative_fp.elevs[0] - representative_fp.elevs[representative_fp.elevs.len() - 1];
    let slope_scalar = if length > 0.0 {
        elev_drop / length
    } else {
        0.0
    };

    Flowpath::new(
        representative_fp.indices.clone(),
        representative_fp.head,
        representative_fp.tail,
        (centroid_px.0 as i32, centroid_px.1 as i32),
        representative_fp.distances_norm.clone(),
        representative_fp.slopes.clone(),
        representative_fp.elevs.clone(),
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
}

fn abstract_subcatchments_representative(
    subwta: &Raster<i32>,
    subwta_indices: &HashMap<i32, Vec<usize>>,
    relief: &Raster<f32>,
    flovec: &Raster<u8>,
    fvslop: &Raster<f32>,
    taspec: &Raster<f32>,
    discha: &Raster<f32>,
    channels: &FlowpathCollection,
) -> FlowpathCollection {
    let topaz_ids: Vec<i32> = subwta_indices
        .keys()
        .filter(|&&topaz_id| topaz_id % 10 != 4)
        .cloned()
        .collect();

    let results: Vec<Flowpath> = topaz_ids
        .into_par_iter()
        .map(|topaz_id| {
            let indices = subwta_indices
                .get(&topaz_id)
                .expect("missing indices for topaz_id");
            let representative_fp = select_representative_flowpath(
                topaz_id, indices, subwta, relief, flovec, fvslop, taspec, discha,
            );
            build_representative_hillslope(&representative_fp, subwta, taspec, channels, indices)
        })
        .collect();

    FlowpathCollection {
        flowpaths: results,
        subflows: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::watershed_abstraction::walk_flowpaths;
    use std::path::PathBuf;

    fn fixture_wbt_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/wbt/amiss-polyhedron/dem/wbt")
    }

    fn median(values: &mut [f64]) -> f64 {
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = values.len();
        if n % 2 == 0 {
            (values[n / 2 - 1] + values[n / 2]) / 2.0
        } else {
            values[n / 2]
        }
    }

    #[test]
    fn representative_flowpath_uses_source_cells() {
        let fixture = fixture_wbt_path();

        let subwta = Raster::<i32>::read(fixture.join("subwta.tif").to_str().unwrap()).unwrap();
        let relief = Raster::<f32>::read(fixture.join("relief.tif").to_str().unwrap()).unwrap();
        let mut flovec = Raster::<u8>::read(fixture.join("flovec.tif").to_str().unwrap()).unwrap();
        let fvslop = Raster::<f32>::read(fixture.join("fvslop.tif").to_str().unwrap()).unwrap();
        let taspec = Raster::<f32>::read(fixture.join("taspec.tif").to_str().unwrap()).unwrap();
        let discha = Raster::<f32>::read(fixture.join("discha.tif").to_str().unwrap()).unwrap();

        remap_whitebox_d8_to_topaz_in_place(&mut flovec);

        let indices_map = subwta.indices_map();
        for (topaz_id, indices) in indices_map {
            if topaz_id <= 0 || topaz_id % 10 == 4 {
                continue;
            }

            let flowpaths = walk_flowpaths(
                topaz_id, &indices, &subwta, &relief, &flovec, &fvslop, &taspec,
            );
            let edge_flowpaths = flowpaths.get_edge_flowpaths2(&subwta, &flovec);
            assert!(
                !edge_flowpaths.is_empty(),
                "missing edge flowpaths for topaz_id {}",
                topaz_id
            );

            let edge_heads: HashSet<usize> = edge_flowpaths
                .iter()
                .filter_map(|fp| fp.indices.first().copied())
                .collect();

            let rep_fp = select_representative_flowpath(
                topaz_id, &indices, &subwta, &relief, &flovec, &fvslop, &taspec, &discha,
            );
            let rep_head = *rep_fp
                .indices
                .first()
                .expect("representative flowpath missing head");
            assert!(
                edge_heads.contains(&rep_head),
                "representative head not a source cell for topaz_id {}",
                topaz_id
            );

            let mut edge_lengths: Vec<f64> = edge_flowpaths.iter().map(|fp| fp.length).collect();
            let min_len = edge_lengths.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_len = edge_lengths
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            let rep_len = rep_fp.length;
            assert!(
                rep_len >= min_len - 1e-6 && rep_len <= max_len + 1e-6,
                "representative length out of edge range for topaz_id {}",
                topaz_id
            );

            if topaz_id % 10 == 1 {
                let median_len = median(&mut edge_lengths);
                println!(
                    "representative_flowpath topaz_id={} edge_count={} median_len={:.2} rep_len={:.2} ratio={:.3}",
                    topaz_id,
                    edge_lengths.len(),
                    median_len,
                    rep_len,
                    rep_len / median_len
                );
            }
        }
    }
}
