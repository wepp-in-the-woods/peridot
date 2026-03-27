use std::collections::HashSet;
use std::fmt;

use serde::Serialize;

use crate::raster::Raster;

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TraceTerminationReason {
    HitChannel,
    InvalidFlowDirection,
    LoopDetected,
    RasterEdge,
    MaxStepsExceeded,
}

#[derive(Debug, Clone, Serialize)]
pub struct TraceDownslopeResult {
    pub seed_row: i32,
    pub seed_col: i32,
    pub seed_topaz_id: i32,
    pub reaches_channel: bool,
    pub channel_row: Option<i32>,
    pub channel_col: Option<i32>,
    pub channel_topaz_id: Option<i32>,
    pub termination_reason: TraceTerminationReason,
    pub rows: Vec<i32>,
    pub cols: Vec<i32>,
    pub indices: Vec<usize>,
    pub distance_m: Vec<f64>,
    pub elevation_m: Vec<f64>,
    pub segment_slope: Vec<f64>,
    pub path_length_m: f64,
    pub drop_m: f64,
    pub mean_slope: f64,
    pub max_slope: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TraceError {
    SeedOutOfBounds {
        seed_row: usize,
        seed_col: usize,
        width: usize,
        height: usize,
    },
    RasterShapeMismatch {
        raster_name: &'static str,
        expected_width: usize,
        expected_height: usize,
        actual_width: usize,
        actual_height: usize,
    },
    InvalidMaxSteps {
        max_steps: usize,
    },
}

impl fmt::Display for TraceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TraceError::SeedOutOfBounds {
                seed_row,
                seed_col,
                width,
                height,
            } => write!(
                f,
                "seed row/col ({}, {}) is outside raster bounds width={} height={}",
                seed_row, seed_col, width, height
            ),
            TraceError::RasterShapeMismatch {
                raster_name,
                expected_width,
                expected_height,
                actual_width,
                actual_height,
            } => write!(
                f,
                "{} raster shape mismatch: expected {}x{}, got {}x{}",
                raster_name, expected_width, expected_height, actual_width, actual_height
            ),
            TraceError::InvalidMaxSteps { max_steps } => {
                write!(f, "max_steps must be >= 1, got {}", max_steps)
            }
        }
    }
}

impl std::error::Error for TraceError {}

#[derive(Debug)]
struct TraceState {
    rows: Vec<i32>,
    cols: Vec<i32>,
    indices: Vec<usize>,
    distance_m: Vec<f64>,
    elevation_m: Vec<f64>,
    segment_slope: Vec<f64>,
    channel_row: Option<i32>,
    channel_col: Option<i32>,
    channel_topaz_id: Option<i32>,
    termination_reason: TraceTerminationReason,
}

pub fn trace_downslope_flowpath(
    subwta: &Raster<i32>,
    flovec: &Raster<u8>,
    relief: &Raster<f32>,
    seed_row: usize,
    seed_col: usize,
    channel_mask: Option<&Raster<i32>>,
    max_steps: usize,
) -> Result<TraceDownslopeResult, TraceError> {
    if max_steps == 0 {
        return Err(TraceError::InvalidMaxSteps { max_steps });
    }

    validate_raster_shapes(subwta, flovec, relief, channel_mask)?;

    if seed_row >= subwta.height || seed_col >= subwta.width {
        return Err(TraceError::SeedOutOfBounds {
            seed_row,
            seed_col,
            width: subwta.width,
            height: subwta.height,
        });
    }

    let seed_index = subwta.xy_to_index(seed_col, seed_row);
    let seed_topaz_id = subwta.data[seed_index];

    let mut state = TraceState {
        rows: vec![seed_row as i32],
        cols: vec![seed_col as i32],
        indices: vec![seed_index],
        distance_m: vec![0.0],
        elevation_m: vec![relief.data[seed_index] as f64],
        segment_slope: Vec::new(),
        channel_row: None,
        channel_col: None,
        channel_topaz_id: None,
        termination_reason: TraceTerminationReason::MaxStepsExceeded,
    };

    if is_channel_cell(seed_index, subwta, channel_mask) {
        state.channel_row = Some(seed_row as i32);
        state.channel_col = Some(seed_col as i32);
        state.channel_topaz_id = Some(seed_topaz_id);
        state.termination_reason = TraceTerminationReason::HitChannel;
        return Ok(finalize_result(seed_row, seed_col, seed_topaz_id, state));
    }

    let mut visited: HashSet<usize> = HashSet::new();
    visited.insert(seed_index);

    let mut current_index = seed_index;
    let mut steps = 0usize;

    loop {
        if steps >= max_steps {
            state.termination_reason = TraceTerminationReason::MaxStepsExceeded;
            break;
        }

        let flow_dir = flovec.data[current_index];
        let (dx, dy) = match d8_offset(flow_dir) {
            Some(offset) => offset,
            None => {
                state.termination_reason = TraceTerminationReason::InvalidFlowDirection;
                break;
            }
        };

        let (x, y) = subwta.index_to_xy(current_index);
        let next_x = x as isize + dx;
        let next_y = y as isize + dy;

        if next_x < 0
            || next_y < 0
            || next_x >= subwta.width as isize
            || next_y >= subwta.height as isize
        {
            state.termination_reason = TraceTerminationReason::RasterEdge;
            break;
        }

        let next_index = subwta.xy_to_index(next_x as usize, next_y as usize);
        if visited.contains(&next_index) {
            state.termination_reason = TraceTerminationReason::LoopDetected;
            break;
        }

        let segment_distance = subwta.distance_between(current_index, next_index);
        let previous_distance = *state.distance_m.last().unwrap_or(&0.0);
        state.distance_m.push(previous_distance + segment_distance);

        let previous_elev = *state.elevation_m.last().unwrap_or(&0.0);
        let next_elev = relief.data[next_index] as f64;
        state.elevation_m.push(next_elev);

        let slope = if segment_distance > 0.0 {
            (previous_elev - next_elev) / segment_distance
        } else {
            0.0
        };
        state.segment_slope.push(slope);

        state.rows.push(next_y as i32);
        state.cols.push(next_x as i32);
        state.indices.push(next_index);

        if is_channel_cell(next_index, subwta, channel_mask) {
            state.channel_row = Some(next_y as i32);
            state.channel_col = Some(next_x as i32);
            state.channel_topaz_id = Some(subwta.data[next_index]);
            state.termination_reason = TraceTerminationReason::HitChannel;
            break;
        }

        visited.insert(next_index);
        current_index = next_index;
        steps += 1;
    }

    Ok(finalize_result(seed_row, seed_col, seed_topaz_id, state))
}

fn validate_raster_shapes(
    subwta: &Raster<i32>,
    flovec: &Raster<u8>,
    relief: &Raster<f32>,
    channel_mask: Option<&Raster<i32>>,
) -> Result<(), TraceError> {
    if flovec.width != subwta.width || flovec.height != subwta.height {
        return Err(TraceError::RasterShapeMismatch {
            raster_name: "flovec",
            expected_width: subwta.width,
            expected_height: subwta.height,
            actual_width: flovec.width,
            actual_height: flovec.height,
        });
    }

    if relief.width != subwta.width || relief.height != subwta.height {
        return Err(TraceError::RasterShapeMismatch {
            raster_name: "relief",
            expected_width: subwta.width,
            expected_height: subwta.height,
            actual_width: relief.width,
            actual_height: relief.height,
        });
    }

    if let Some(mask) = channel_mask {
        if mask.width != subwta.width || mask.height != subwta.height {
            return Err(TraceError::RasterShapeMismatch {
                raster_name: "channel_mask",
                expected_width: subwta.width,
                expected_height: subwta.height,
                actual_width: mask.width,
                actual_height: mask.height,
            });
        }
    }

    Ok(())
}

fn d8_offset(flow_dir: u8) -> Option<(isize, isize)> {
    match flow_dir {
        1 => Some((-1, -1)),
        2 => Some((0, -1)),
        3 => Some((1, -1)),
        4 => Some((-1, 0)),
        6 => Some((1, 0)),
        7 => Some((-1, 1)),
        8 => Some((0, 1)),
        9 => Some((1, 1)),
        _ => None,
    }
}

fn is_channel_cell(index: usize, subwta: &Raster<i32>, channel_mask: Option<&Raster<i32>>) -> bool {
    match channel_mask {
        Some(mask) => mask.data[index] > 0,
        None => subwta.data[index] % 10 == 4,
    }
}

fn finalize_result(
    seed_row: usize,
    seed_col: usize,
    seed_topaz_id: i32,
    state: TraceState,
) -> TraceDownslopeResult {
    let path_length_m = *state.distance_m.last().unwrap_or(&0.0);

    let first_elev = *state.elevation_m.first().unwrap_or(&0.0);
    let last_elev = *state.elevation_m.last().unwrap_or(&first_elev);
    let drop_m = first_elev - last_elev;

    let mean_slope = if path_length_m > 0.0 {
        drop_m / path_length_m
    } else {
        0.0
    };

    let max_slope = state
        .segment_slope
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .reduce(f64::max)
        .unwrap_or(0.0);

    let reaches_channel = state.termination_reason == TraceTerminationReason::HitChannel;

    TraceDownslopeResult {
        seed_row: seed_row as i32,
        seed_col: seed_col as i32,
        seed_topaz_id,
        reaches_channel,
        channel_row: state.channel_row,
        channel_col: state.channel_col,
        channel_topaz_id: state.channel_topaz_id,
        termination_reason: state.termination_reason,
        rows: state.rows,
        cols: state.cols,
        indices: state.indices,
        distance_m: state.distance_m,
        elevation_m: state.elevation_m,
        segment_slope: state.segment_slope,
        path_length_m,
        drop_m,
        mean_slope,
        max_slope,
    }
}
