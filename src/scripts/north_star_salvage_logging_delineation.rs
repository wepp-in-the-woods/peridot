
use std::collections::{HashSet, HashMap};

use peridot::watershed_abstraction::{FlowpathCollection};

use peridot::skid_trail_segmentation::{
    rasterize_skid_trails, 
    segment_skid_trails_raster, 
    catchment_skid_trails_raster, 
    identify_skid_channels
};


use peridot::catchment_trace::{
    walk_skid_flowpath
};

use peridot::raster::{utm_to_px, Raster};

fn main() {

    let max_points = 1000;
    let clip_hillslopes = false;
    let clip_hillslope_length = 1000.0;

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

    let fvslop_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/dem/topaz/FVSLOP.ARC";
    let fvslop = Raster::<f64>::read(&fvslop_fn).unwrap();
    
    let taspec_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/dem/topaz/TASPEC.ARC";
    let taspec = Raster::<f64>::read(&taspec_fn).unwrap();

    
    let wd = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/";

    let dst_fn3 = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/Skid_Catchments.tif";
    catchment_skid_trails_raster(&skid, &seg_skids, &flovec, &relief, &dst_fn3);


    let skid_catchments = Raster::<i32>::read(&dst_fn3).unwrap();

    let hillslope_ids = skid_catchments.unique_values();

    // iterate over skid ids and build flowpath collection for each by walking dwon from each pixel in the catchment
    // and recording the flowpath

    for hill_id in hillslope_ids {
        println!("hill_id: {:?}", hill_id);

        if hill_id == 0 {
            continue;
        }

        let mut flowpaths = vec![];

        // need upseg map for each hill_id
        let starting_indices = skid_catchments.indices_of(hill_id);

        let mut _indices: HashSet<usize> = HashSet::new();

        for start in &starting_indices {
            println!("hill_id: {:?}, start: {:?}", hill_id, start);
            let mut fp = walk_skid_flowpath(
                *start,
                &skid_catchments,
                &relief,
                &flovec,
                &fvslop,
                &taspec,
                hill_id);
            // add fp.indices to vec_indices
            _indices.extend(fp.indices.iter().cloned());

            flowpaths.push(fp);
        }

        let mut vec_indices: Vec<usize> = _indices.into_iter().collect();

        println!("flowpaths: {:?}", flowpaths);

        if flowpaths.len() == 0 {
            continue;
        }

        let mut flowpath_collection = FlowpathCollection  {
            flowpaths: flowpaths,
            subflows: Some(HashMap::<i32, FlowpathCollection>::new())
        };

        let hillslope = flowpath_collection.abstract_hillslope(
            &flovec, &taspec, &vec_indices);

        hillslope.write_slp(
            &format!("{}/skid_hillslope_{}.slp", wd, hill_id), 
            max_points, clip_hillslopes, clip_hillslope_length);
    }

    let dst_fn4 = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/Skid_Channels.tif";
    let seg_links_fn = "tests/fixtures/watershed_abstraction/bearded-flyover/skid_trails/segment_links.txt";
    identify_skid_channels(&seg_skids, &flovec, &relief, &bound, &dst_fn4, &seg_links_fn);

}

//> cd /workdir/peridot
//> RUST_BACKTRACE=1 cargo run --bin north_star_salvage
