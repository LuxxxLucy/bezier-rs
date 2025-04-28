//! Heuristics for t-parameter estimation used in the least square curve fitting

use crate::data::Point;

/// Heuristics for t-parameter estimation used in the least square curve fitting
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum THeuristic {
    /// Chord length - t values based on linear distance between points
    ///
    /// This assigns t values proportionally to the distance traveled along the polyline
    /// formed by the input points, producing good parameterization for most curves.
    #[default]
    ChordLength,

    /// Uniform spacing - t values are evenly distributed
    ///
    /// This assigns t values uniformly in \[0,1\], which works well for evenly spaced points
    /// but may perform poorly for unevenly distributed points.
    Uniform,

    /// Centripetal - t values based on square root of chord length
    ///
    /// This assigns t values based on the square root of the chord length, which can
    /// help prevent overshooting in curves with sharp turns.
    Centripetal,
}

/// Estimate parameter t values using chord length parameterization
///
/// This implementation directly calculates the chord length parameterization as described
/// in the Bezier primer's Curve Fitting chapter. It assigns t values proportionally to
/// the distance traveled along the polyline formed by the input points.
pub fn estimate_t_values_chord_length(points: &[Point]) -> Vec<f64> {
    if points.is_empty() {
        return Vec::new();
    }

    if points.len() == 1 {
        return vec![0.0];
    }

    // Calculate the path length to parameterize the points
    let mut path_lengths = vec![0.0];
    let mut total_length = 0.0;

    for i in 1..points.len() {
        let segment_length = points[i].distance(&points[i - 1]);
        total_length += segment_length;
        path_lengths.push(total_length);
    }

    // Normalize path lengths to get parameter t values
    path_lengths
        .iter()
        .map(|&length| {
            if total_length > 0.0 {
                length / total_length
            } else {
                0.0
            }
        })
        .collect()
}

/// Estimate parameter t values using uniform spacing
pub fn estimate_t_values_uniform(points: &[Point]) -> Vec<f64> {
    if points.is_empty() {
        return Vec::new();
    }

    if points.len() == 1 {
        return vec![0.0];
    }

    (0..points.len())
        .map(|i| i as f64 / (points.len() - 1) as f64)
        .collect()
}

/// Estimate parameter t values using centripetal parameterization
pub fn estimate_t_values_centripetal(points: &[Point]) -> Vec<f64> {
    if points.is_empty() {
        return Vec::new();
    }

    if points.len() == 1 {
        return vec![0.0];
    }

    // Calculate the path length to parameterize the points
    let mut path_lengths = vec![0.0];
    let mut total_length = 0.0;

    for i in 1..points.len() {
        let segment_length = points[i].distance(&points[i - 1]).sqrt();
        total_length += segment_length;
        path_lengths.push(total_length);
    }

    // Normalize path lengths to get parameter t values
    path_lengths
        .iter()
        .map(|&length| {
            if total_length > 0.0 {
                length / total_length
            } else {
                0.0
            }
        })
        .collect()
}

/// Estimate parameter t values using specified heuristic
pub fn estimate_t_values_with_heuristic(points: &[Point], heuristic: THeuristic) -> Vec<f64> {
    match heuristic {
        THeuristic::ChordLength => estimate_t_values_chord_length(points),
        THeuristic::Uniform => estimate_t_values_uniform(points),
        THeuristic::Centripetal => estimate_t_values_centripetal(points),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::Point;

    // Test data from main.rs
    const TEST_T_VALUES: [f64; 7] = [0.0, 0.1, 0.15, 0.3, 0.7, 0.85, 1.0];

    fn create_test_curve() -> Vec<Point> {
        let segment = crate::cubic!([(50, 200), (100, 50), (200, 50), (250, 200)]);
        let (points, _) = segment.sample_at_t_values(&TEST_T_VALUES);
        points
    }

    #[test]
    fn test_chord_length_heuristic() {
        let points = create_test_curve();
        let t_values = estimate_t_values_chord_length(&points);

        // Verify basic properties
        assert_eq!(t_values.len(), points.len());
        assert_eq!(t_values[0], 0.0);
        assert_eq!(t_values[t_values.len() - 1], 1.0);
        assert!(t_values.windows(2).all(|w| w[0] <= w[1])); // Monotonic increasing

        // Verify against expected values (with some tolerance)
        let expected_t_values = vec![
            0.0,
            0.14100483007137046,
            0.20293241410471896,
            0.3574976594382045,
            0.6427698285208382,
            0.7973350738543237,
            1.0,
        ];

        for (actual, expected) in t_values.iter().zip(expected_t_values.iter()) {
            assert!((actual - expected).abs() < 0.01);
        }
    }

    #[test]
    fn test_uniform_heuristic() {
        let points = create_test_curve();
        let t_values = estimate_t_values_uniform(&points);

        // Verify basic properties
        assert_eq!(t_values.len(), points.len());
        assert_eq!(t_values[0], 0.0);
        assert_eq!(t_values[t_values.len() - 1], 1.0);
        assert!(t_values.windows(2).all(|w| w[0] <= w[1])); // Monotonic increasing

        // Verify uniform spacing
        let expected_t_values = vec![
            0.0,
            0.16666666666666666,
            0.3333333333333333,
            0.5,
            0.6666666666666666,
            0.8333333333333334,
            1.0,
        ];

        for (actual, expected) in t_values.iter().zip(expected_t_values.iter()) {
            assert!((actual - expected).abs() < 0.0001);
        }
    }

    #[test]
    fn test_centripetal_heuristic() {
        let points = create_test_curve();
        let t_values = estimate_t_values_centripetal(&points);

        // Verify basic properties
        assert_eq!(t_values.len(), points.len());
        assert_eq!(t_values[0], 0.0);
        assert_eq!(t_values[t_values.len() - 1], 1.0);
        assert!(t_values.windows(2).all(|w| w[0] <= w[1])); // Monotonic increasing

        // Verify against expected values (with some tolerance)
        let expected_t_values = vec![
            0.0,
            0.1567910277889129,
            0.26069838064703205,
            0.4248556563192105,
            0.6478705739969236,
            0.812027849669102,
            1.0,
        ];

        for (actual, expected) in t_values.iter().zip(expected_t_values.iter()) {
            assert!((actual - expected).abs() < 0.01);
        }
    }

    #[test]
    fn test_heuristic_selection() {
        let points = create_test_curve();

        // Test each heuristic
        let heuristics = [
            THeuristic::ChordLength,
            THeuristic::Uniform,
            THeuristic::Centripetal,
        ];

        for heuristic in heuristics.iter() {
            let t_values = estimate_t_values_with_heuristic(&points, *heuristic);

            // Verify basic properties for all heuristics
            assert_eq!(t_values.len(), points.len());
            assert_eq!(t_values[0], 0.0);
            assert_eq!(t_values[t_values.len() - 1], 1.0);
            assert!(t_values.windows(2).all(|w| w[0] <= w[1])); // Monotonic increasing
        }
    }
}
