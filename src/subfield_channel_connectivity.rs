use std::collections::HashSet;
use std::fmt;

use serde::Serialize;

use crate::raster::Raster;
use crate::watershed_abstraction::PATHS;

/// Simple, raster-derived metrics for direct subfield-to-channel connectivity.
#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct SubfieldChannelConnectivitySummary {
    pub subfields_total: usize,
    pub subfields_with_direct_channel_drainage: usize,
    pub subfields_without_direct_channel_drainage: usize,
    pub direct_channel_outlet_cells: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubfieldChannelConnectivityError {
    RasterShapeMismatch {
        raster_name: &'static str,
        expected_width: usize,
        expected_height: usize,
        actual_width: usize,
        actual_height: usize,
    },
    RasterGeoTransformMismatch {
        raster_name: &'static str,
    },
    RasterProjectionMismatch {
        raster_name: &'static str,
    },
}

impl fmt::Display for SubfieldChannelConnectivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RasterShapeMismatch {
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
            Self::RasterGeoTransformMismatch { raster_name } => write!(
                f,
                "{} raster geotransform does not align with sub_field_map",
                raster_name
            ),
            Self::RasterProjectionMismatch { raster_name } => write!(
                f,
                "{} raster projection does not match sub_field_map",
                raster_name
            ),
        }
    }
}

impl std::error::Error for SubfieldChannelConnectivityError {}

/// Count retained subfields with at least one per-pixel flowpath that enters a channel directly.
///
/// Peridot generates one subfield flowpath from every positive cell in `sub_field_map` and
/// stops each path after appending its first cell outside that subfield. A subfield therefore
/// has a direct channel-draining flowpath exactly when at least one of its cells has a valid
/// D8 successor that is a channel cell. When `channel_mask` is supplied, positive mask values
/// identify channels; otherwise the established TOPAZ convention (`subwta_id % 10 == 4`) is
/// used.
///
/// `flovec` must already use Peridot/TOPAZ D8 codes. WBT callers should remap the raster
/// with `remap_whitebox_d8_to_topaz_in_place` before calling this function.
pub fn summarize_subfield_channel_connectivity(
    sub_field_map: &Raster<i32>,
    subwta: &Raster<i32>,
    flovec: &Raster<u8>,
    channel_mask: Option<&Raster<i32>>,
) -> Result<SubfieldChannelConnectivitySummary, SubfieldChannelConnectivityError> {
    validate_raster_grid("subwta", sub_field_map, subwta)?;
    validate_raster_grid("flovec", sub_field_map, flovec)?;
    if let Some(mask) = channel_mask {
        validate_raster_grid("channel_mask", sub_field_map, mask)?;
    }

    let mut subfield_ids = HashSet::new();
    let mut directly_connected_subfield_ids = HashSet::new();
    let mut direct_channel_outlet_cells = 0usize;

    for (index, &subfield_id) in sub_field_map.data.iter().enumerate() {
        if subfield_id <= 0 {
            continue;
        }
        subfield_ids.insert(subfield_id);

        let flow_direction = flovec.data[index] as i32;
        let Some(&(dx, dy)) = PATHS.get(&flow_direction) else {
            continue;
        };

        let (x, y) = sub_field_map.index_to_xy(index);
        let next_x = x as isize + dx;
        let next_y = y as isize + dy;
        if next_x < 0
            || next_y < 0
            || next_x >= sub_field_map.width as isize
            || next_y >= sub_field_map.height as isize
        {
            continue;
        }

        let next_index = sub_field_map.xy_to_index(next_x as usize, next_y as usize);
        if sub_field_map.data[next_index] == subfield_id {
            continue;
        }
        let is_channel = match channel_mask {
            Some(mask) => mask.data[next_index] > 0,
            None => subwta.data[next_index] % 10 == 4,
        };
        if is_channel {
            directly_connected_subfield_ids.insert(subfield_id);
            direct_channel_outlet_cells += 1;
        }
    }

    let subfields_total = subfield_ids.len();
    let subfields_with_direct_channel_drainage = directly_connected_subfield_ids.len();
    Ok(SubfieldChannelConnectivitySummary {
        subfields_total,
        subfields_with_direct_channel_drainage,
        subfields_without_direct_channel_drainage: subfields_total
            - subfields_with_direct_channel_drainage,
        direct_channel_outlet_cells,
    })
}

fn validate_raster_grid<T, U>(
    raster_name: &'static str,
    expected: &Raster<T>,
    actual: &Raster<U>,
) -> Result<(), SubfieldChannelConnectivityError> {
    if actual.width != expected.width || actual.height != expected.height {
        return Err(SubfieldChannelConnectivityError::RasterShapeMismatch {
            raster_name,
            expected_width: expected.width,
            expected_height: expected.height,
            actual_width: actual.width,
            actual_height: actual.height,
        });
    }
    if !geotransforms_match(&expected.geo_transform, &actual.geo_transform) {
        return Err(SubfieldChannelConnectivityError::RasterGeoTransformMismatch { raster_name });
    }
    if !projections_match(expected.proj4.as_deref(), actual.proj4.as_deref()) {
        return Err(SubfieldChannelConnectivityError::RasterProjectionMismatch { raster_name });
    }
    Ok(())
}

fn geotransforms_match(expected: &[f64; 6], actual: &[f64; 6]) -> bool {
    expected.iter().zip(actual.iter()).all(|(&left, &right)| {
        let scale = left.abs().max(right.abs()).max(1.0);
        (left - right).abs() <= 1e-9 * scale
    })
}

fn projections_match(expected: Option<&str>, actual: Option<&str>) -> bool {
    match (expected, actual) {
        (Some(left), Some(right)) => left.split_whitespace().eq(right.split_whitespace()),
        (None, None) => true,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raster::MapType;

    fn raster_i32(width: usize, height: usize, data: Vec<i32>) -> Raster<i32> {
        Raster::new(
            width,
            height,
            1.0,
            data,
            Some(-9999),
            [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
            None,
            String::new(),
            String::new(),
            MapType::OTHER,
        )
    }

    fn raster_u8(width: usize, height: usize, data: Vec<u8>) -> Raster<u8> {
        Raster::new(
            width,
            height,
            1.0,
            data,
            Some(0),
            [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
            None,
            String::new(),
            String::new(),
            MapType::OTHER,
        )
    }

    #[test]
    fn counts_subfields_with_a_direct_channel_successor() {
        // Two outlet cells in subfield 10 enter channel 14, but the subfield is counted once.
        // Subfield 20 flows into non-channel cell 31 and is not directly connected.
        let sub_field_map = raster_i32(3, 3, vec![10, 10, 0, 10, 10, 0, 20, 20, 0]);
        let subwta = raster_i32(3, 3, vec![11, 11, 14, 11, 11, 14, 31, 31, 31]);
        let flovec = raster_u8(3, 3, vec![6, 6, 0, 6, 6, 0, 6, 6, 0]);

        let summary =
            summarize_subfield_channel_connectivity(&sub_field_map, &subwta, &flovec, None)
                .expect("connectivity summary failed");

        assert_eq!(summary.subfields_total, 2);
        assert_eq!(summary.subfields_with_direct_channel_drainage, 1);
        assert_eq!(summary.subfields_without_direct_channel_drainage, 1);
        assert_eq!(summary.direct_channel_outlet_cells, 2);
    }

    #[test]
    fn positive_channel_mask_overrides_topaz_suffix_rule() {
        let sub_field_map = raster_i32(3, 1, vec![10, 10, 0]);
        let subwta = raster_i32(3, 1, vec![11, 11, 31]);
        let flovec = raster_u8(3, 1, vec![6, 6, 0]);
        let channel_mask = raster_i32(3, 1, vec![0, 0, 1]);

        let summary = summarize_subfield_channel_connectivity(
            &sub_field_map,
            &subwta,
            &flovec,
            Some(&channel_mask),
        )
        .expect("connectivity summary failed");

        assert_eq!(summary.subfields_with_direct_channel_drainage, 1);
        assert_eq!(summary.direct_channel_outlet_cells, 1);
    }

    #[test]
    fn channel_mask_only_counts_the_first_cell_outside_the_subfield() {
        let sub_field_map = raster_i32(3, 1, vec![10, 10, 0]);
        let subwta = raster_i32(3, 1, vec![11, 11, 31]);
        let flovec = raster_u8(3, 1, vec![6, 6, 0]);
        let channel_mask = raster_i32(3, 1, vec![0, 1, 0]);

        let summary = summarize_subfield_channel_connectivity(
            &sub_field_map,
            &subwta,
            &flovec,
            Some(&channel_mask),
        )
        .expect("connectivity summary failed");

        assert_eq!(summary.subfields_with_direct_channel_drainage, 0);
        assert_eq!(summary.direct_channel_outlet_cells, 0);
    }

    #[test]
    fn rejects_shape_mismatch() {
        let sub_field_map = raster_i32(2, 1, vec![10, 0]);
        let subwta = raster_i32(2, 1, vec![11, 14]);
        let flovec = raster_u8(1, 1, vec![0]);

        let error = summarize_subfield_channel_connectivity(&sub_field_map, &subwta, &flovec, None)
            .expect_err("shape mismatch should fail");

        assert_eq!(
            error,
            SubfieldChannelConnectivityError::RasterShapeMismatch {
                raster_name: "flovec",
                expected_width: 2,
                expected_height: 1,
                actual_width: 1,
                actual_height: 1,
            }
        );
    }

    #[test]
    fn rejects_geotransform_mismatch() {
        let sub_field_map = raster_i32(2, 1, vec![10, 0]);
        let mut subwta = raster_i32(2, 1, vec![11, 14]);
        let flovec = raster_u8(2, 1, vec![6, 0]);
        subwta.geo_transform[0] = 10.0;

        let error = summarize_subfield_channel_connectivity(&sub_field_map, &subwta, &flovec, None)
            .expect_err("shifted raster should fail");

        assert_eq!(
            error,
            SubfieldChannelConnectivityError::RasterGeoTransformMismatch {
                raster_name: "subwta",
            }
        );
    }

    #[test]
    fn rejects_projection_mismatch() {
        let mut sub_field_map = raster_i32(2, 1, vec![10, 0]);
        let mut subwta = raster_i32(2, 1, vec![11, 14]);
        let mut flovec = raster_u8(2, 1, vec![6, 0]);
        sub_field_map.proj4 = Some("+proj=utm +zone=11".to_string());
        subwta.proj4 = Some("+proj=utm +zone=12".to_string());
        flovec.proj4 = sub_field_map.proj4.clone();

        let error = summarize_subfield_channel_connectivity(&sub_field_map, &subwta, &flovec, None)
            .expect_err("different projection should fail");

        assert_eq!(
            error,
            SubfieldChannelConnectivityError::RasterProjectionMismatch {
                raster_name: "subwta",
            }
        );
    }
}
