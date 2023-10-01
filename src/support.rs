use std::f64::consts::PI;
use interp::interp_slice;

#[allow(dead_code)]
pub fn compute_direction(head: (i32, i32), tail: (i32, i32)) -> f64 {
    let a = (tail.1 as f64 - head.1 as f64).atan2(head.0 as f64 - tail.0 as f64) 
           * (180.0 / PI) - 180.0;

    if a < 0.0 {
        return a + 360.0;
    }
    a
}

/// Computes the circular mean of a slice of angles in radians.
///
/// The circular mean is calculated using the trigonometric
/// representation of the angles. The function expects angles
/// to be in radians and returns the mean angle in radians.
///
/// # Arguments
///
/// * `angles` - A slice of angles in radians.
///
/// # Returns
///
/// Returns the circular mean of the given angles in radians.
#[allow(dead_code)]
pub fn circmean(angles: &[f64]) -> f64 {
    let mut sum_sin = 0.0;
    let mut sum_cos = 0.0;

    for &angle in angles {
        sum_sin += angle.sin();
        sum_cos += angle.cos();
    }

    sum_sin /= angles.len() as f64;
    sum_cos /= angles.len() as f64;

    sum_sin.atan2(sum_cos)
}

#[allow(dead_code)]
pub fn interpolate_slp(
    distances: &Vec<f64>, 
    slopes: &Vec<f64>, 
    max_points: usize) -> Result<(Vec<f64>, Vec<f64>), &'static str> {

    // Check for the special case where slope is a single cell
    if slopes.len() == 1{
        return Ok((vec![0.0, 1.0], vec![slopes[0], slopes[0]]));
    }

    assert!(distances.len() == slopes.len());
    assert!(max_points < distances.len(), "{} < {}", max_points, distances.len());

    let new_distances: Vec<f64> = (0..max_points).map(|i| i as f64 / (max_points - 1) as f64).collect();
    let new_slopes: Vec<f64> = interp_slice(distances, slopes, &new_distances);

    Ok((new_distances, new_slopes))
}



#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{self, BufRead, BufReader};

    use super::*;

    #[test]
    fn test_interpolate_single_slope() {
        let mut distances = vec![0.0, 1.0];
        let mut slopes = vec![0.2];
        let result = interpolate_slp(&mut distances, &mut slopes, 10).unwrap();
        assert_eq!(result.0, vec![0.0, 1.0]);
        assert_eq!(result.1, vec![0.2, 0.2]);
    }

    #[test]
    fn test_interpolate_many_points() {
        let mut distances = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let mut slopes = vec![0.2, 0.3, 0.4, 0.5, 0.6];
        let result = interpolate_slp(&mut distances, &mut slopes, 3).unwrap();

        // Since we're interpolating to only 3 points, the expected values might be approximate.
        assert_eq!(result.0, vec![0.0, 0.5, 1.0]);
        println!("{:?}", result);
        assert!((result.1[0] - 0.2).abs() < 1e-10);
        assert!((result.1[1] - 0.4).abs() < 1e-10);
        assert!((result.1[2] - 0.6).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "distances[0] should be 0")]
    fn test_invalid_start_distance() {
        let mut distances = vec![0.1, 1.0];
        let mut slopes = vec![0.2];
        let _ = interpolate_slp(&mut distances, &mut slopes, 10).unwrap();
    }

    #[test]
    #[should_panic(expected = "distances[-1] should be 1")]
    fn test_invalid_end_distance() {
        let mut distances = vec![0.0, 0.9];
        let mut slopes = vec![0.2];
        let _ = interpolate_slp(&mut distances, &mut slopes, 10).unwrap();
    }

}
