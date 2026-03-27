use peridot::raster::{MapType, Raster};
use peridot::roads_trace::{
    trace_downslope_flowpath, TraceDownslopeResult, TraceTerminationReason,
};

fn raster_i32(width: usize, height: usize, data: Vec<i32>) -> Raster<i32> {
    Raster::new(
        width,
        height,
        1.0,
        data,
        None,
        [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
        None,
        String::from("in_memory"),
        String::from("subwta"),
        MapType::OTHER,
    )
}

fn raster_u8(width: usize, height: usize, data: Vec<u8>) -> Raster<u8> {
    Raster::new(
        width,
        height,
        1.0,
        data,
        None,
        [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
        None,
        String::from("in_memory"),
        String::from("flovec"),
        MapType::OTHER,
    )
}

fn raster_f32(width: usize, height: usize, data: Vec<f32>) -> Raster<f32> {
    Raster::new(
        width,
        height,
        1.0,
        data,
        None,
        [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
        None,
        String::from("in_memory"),
        String::from("relief"),
        MapType::OTHER,
    )
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-12,
        "expected {}, got {}",
        expected,
        actual
    );
}

fn assert_vec_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "vector length mismatch");
    for (a, e) in actual.iter().zip(expected.iter()) {
        assert_close(*a, *e);
    }
}

fn assert_lengths_consistent(trace: &TraceDownslopeResult) {
    assert_eq!(trace.rows.len(), trace.cols.len());
    assert_eq!(trace.rows.len(), trace.indices.len());
    assert_eq!(trace.rows.len(), trace.distance_m.len());
    assert_eq!(trace.rows.len(), trace.elevation_m.len());
    assert_eq!(trace.segment_slope.len() + 1, trace.rows.len());
}

#[test]
fn channel_hit_with_channel_mask_has_deterministic_profile_vectors() {
    let subwta = raster_i32(4, 1, vec![11, 11, 11, 11]);
    let flovec = raster_u8(4, 1, vec![6, 6, 6, 0]);
    let relief = raster_f32(4, 1, vec![100.0, 99.0, 98.0, 97.0]);
    let channel = raster_i32(4, 1, vec![0, 0, 1, 1]);

    let trace = trace_downslope_flowpath(&subwta, &flovec, &relief, 0, 0, Some(&channel), 100)
        .expect("trace should succeed");

    assert_eq!(trace.termination_reason, TraceTerminationReason::HitChannel);
    assert!(trace.reaches_channel);
    assert_eq!(trace.channel_row, Some(0));
    assert_eq!(trace.channel_col, Some(2));
    assert_eq!(trace.channel_topaz_id, Some(11));
    assert_eq!(trace.rows, vec![0, 0, 0]);
    assert_eq!(trace.cols, vec![0, 1, 2]);
    assert_eq!(trace.indices, vec![0, 1, 2]);
    assert_vec_close(&trace.distance_m, &[0.0, 1.0, 2.0]);
    assert_vec_close(&trace.elevation_m, &[100.0, 99.0, 98.0]);
    assert_vec_close(&trace.segment_slope, &[1.0, 1.0]);
    assert_close(trace.path_length_m, 2.0);
    assert_close(trace.drop_m, 2.0);
    assert_close(trace.mean_slope, 1.0);
    assert_close(trace.max_slope, 1.0);

    assert_lengths_consistent(&trace);
}

#[test]
fn channel_hit_without_channel_mask_uses_topaz_suffix_rule() {
    let subwta = raster_i32(2, 1, vec![11, 14]);
    let flovec = raster_u8(2, 1, vec![6, 0]);
    let relief = raster_f32(2, 1, vec![5.0, 4.0]);

    let trace = trace_downslope_flowpath(&subwta, &flovec, &relief, 0, 0, None, 100)
        .expect("trace should succeed");

    assert_eq!(trace.termination_reason, TraceTerminationReason::HitChannel);
    assert!(trace.reaches_channel);
    assert_eq!(trace.channel_col, Some(1));
    assert_eq!(trace.channel_topaz_id, Some(14));
    assert_lengths_consistent(&trace);
}

#[test]
fn invalid_flow_direction_terminates_explicitly() {
    let subwta = raster_i32(2, 1, vec![11, 11]);
    let flovec = raster_u8(2, 1, vec![5, 0]);
    let relief = raster_f32(2, 1, vec![10.0, 9.0]);

    let trace = trace_downslope_flowpath(&subwta, &flovec, &relief, 0, 0, None, 100)
        .expect("trace should succeed");

    assert_eq!(
        trace.termination_reason,
        TraceTerminationReason::InvalidFlowDirection
    );
    assert!(!trace.reaches_channel);
    assert_eq!(trace.rows, vec![0]);
    assert_eq!(trace.cols, vec![0]);
    assert_eq!(trace.indices, vec![0]);
    assert_close(trace.path_length_m, 0.0);
    assert_lengths_consistent(&trace);
}

#[test]
fn loop_detection_terminates_before_revisiting_cell() {
    let subwta = raster_i32(2, 1, vec![11, 11]);
    let flovec = raster_u8(2, 1, vec![6, 4]);
    let relief = raster_f32(2, 1, vec![10.0, 9.0]);

    let trace = trace_downslope_flowpath(&subwta, &flovec, &relief, 0, 0, None, 100)
        .expect("trace should succeed");

    assert_eq!(
        trace.termination_reason,
        TraceTerminationReason::LoopDetected
    );
    assert!(!trace.reaches_channel);
    assert_eq!(trace.rows, vec![0, 0]);
    assert_eq!(trace.cols, vec![0, 1]);
    assert_eq!(trace.indices, vec![0, 1]);
    assert_close(trace.path_length_m, 1.0);
    assert_lengths_consistent(&trace);
}

#[test]
fn raster_edge_termination_is_explicit() {
    let subwta = raster_i32(1, 1, vec![11]);
    let flovec = raster_u8(1, 1, vec![6]);
    let relief = raster_f32(1, 1, vec![10.0]);

    let trace = trace_downslope_flowpath(&subwta, &flovec, &relief, 0, 0, None, 100)
        .expect("trace should succeed");

    assert_eq!(trace.termination_reason, TraceTerminationReason::RasterEdge);
    assert!(!trace.reaches_channel);
    assert_eq!(trace.rows, vec![0]);
    assert_eq!(trace.cols, vec![0]);
    assert_eq!(trace.indices, vec![0]);
    assert_close(trace.path_length_m, 0.0);
    assert_lengths_consistent(&trace);
}

#[test]
fn max_step_limit_terminates_trace() {
    let subwta = raster_i32(3, 1, vec![11, 11, 11]);
    let flovec = raster_u8(3, 1, vec![6, 6, 0]);
    let relief = raster_f32(3, 1, vec![10.0, 9.0, 8.0]);

    let trace = trace_downslope_flowpath(&subwta, &flovec, &relief, 0, 0, None, 1)
        .expect("trace should succeed");

    assert_eq!(
        trace.termination_reason,
        TraceTerminationReason::MaxStepsExceeded
    );
    assert!(!trace.reaches_channel);
    assert_eq!(trace.rows, vec![0, 0]);
    assert_eq!(trace.cols, vec![0, 1]);
    assert_eq!(trace.indices, vec![0, 1]);
    assert_close(trace.path_length_m, 1.0);
    assert_lengths_consistent(&trace);
}
