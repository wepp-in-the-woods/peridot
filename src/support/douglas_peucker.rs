fn point_line_distance(x0: f64, y0: f64, x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let numerator = ((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1).abs();
    let denominator = ((y2 - y1).powi(2) + (x2 - x1).powi(2)).sqrt();
    numerator / denominator
}

fn douglas_peucker_recursive(
    distances: &[f64], 
    slopes: &[f64], 
    start: usize, 
    end: usize, 
    tolerance: f64, 
    result_distances: &mut Vec<f64>, 
    result_slopes: &mut Vec<f64>
) {
    if start >= end {
        return;
    }

    let mut max_distance = 0.0;
    let mut max_index = 0;

    for i in start + 1..end {
        let distance = point_line_distance(
            distances[i], slopes[i],
            distances[start], slopes[start],
            distances[end], slopes[end]
        );
        if distance > max_distance {
            max_distance = distance;
            max_index = i;
        }
    }

    if max_distance > tolerance {
        douglas_peucker_recursive(distances, slopes, start, max_index, tolerance, result_distances, result_slopes);
        douglas_peucker_recursive(distances, slopes, max_index, end, tolerance, result_distances, result_slopes);
    } else {
        result_distances.push(distances[start]);
        result_slopes.push(slopes[start]);
        if end == distances.len() - 1 {
            result_distances.push(distances[end]);
            result_slopes.push(slopes[end]);
        }
    }
}

pub fn douglas_peucker(
    distances_norm: &[f64], 
    slopes: &[f64], 
    tolerance: f64
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    let mut result_distances = Vec::new();
    let mut result_slopes = Vec::new();
    douglas_peucker_recursive(distances_norm, slopes, 0, distances_norm.len() - 1, tolerance, &mut result_distances, &mut result_slopes);
    Ok((result_distances, result_slopes))
}
