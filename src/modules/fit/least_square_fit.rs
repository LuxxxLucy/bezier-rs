//! Fitting bezier curves to a set of points using the least squares method
//!
//! Based on "Least Squares Bezier Fit" by Jim Herold
//! <https://web.archive.org/web/20180403213813/http://jimherold.com/2012/04/20/least-squares-bezier-fit/>
//!
//! Altenatively see bezier primer chap 35 Curve fitting <https://pomax.github.io/bezierinfo/#curvefitting>
//!
//! This approach would rely on a heuristic of guessing the `t` parameter and then apply a least
//! square solving procedure.
//!
//! # Example
//!
//! ```rust
//! use bezier_rs::{modules::fit::least_square_fit, Point, BezierResult};
//!
//! // Some sample points
//! let points = vec![
//!     Point::new(0.0, 0.0),
//!     Point::new(1.0, 1.5),
//!     Point::new(2.0, 1.8),
//!     Point::new(3.0, 0.0),
//! ];
//!
//! // Fit a cubic Bezier to the points
//! let result: BezierResult<_> = least_square_fit::fit_cubic_bezier_default(&points);
//! let fitted = result.unwrap();
//! ```

use super::t_heuristic::{estimate_t_values_with_heuristic, THeuristic};
use crate::data::{BezierSegment, Point};
use crate::error::{BezierError, BezierResult};
use crate::modules::fit::least_square_fit_common::least_square_solve_p_given_t;

/// Fit a cubic bezier curve to a set of points using least squares
///
/// This implementation uses the chord length parameterization for t-value estimation
/// as described in Jim Herold's blog post and the Bezier primer's Curve Fitting chapter.
pub fn fit_cubic_bezier_default(points: &[Point]) -> BezierResult<BezierSegment> {
    fit_cubic_bezier_with_heuristic(points, THeuristic::default())
}

/// Fit a cubic bezier curve to a set of points using least squares with specified heuristic
pub fn fit_cubic_bezier_with_heuristic(
    points: &[Point],
    heuristic: THeuristic,
) -> BezierResult<BezierSegment> {
    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
    }

    // Estimate t values using the specified heuristic
    let t_values = estimate_t_values_with_heuristic(points, heuristic);

    // solve P given T
    least_square_solve_p_given_t(points, &t_values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic; // Import macros

    #[test]
    fn test_fitting() {
        // Create a bezier curve
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);

        // Sample points from the curve
        let samples = original.sample_points(20);

        // Fit a curve to the sampled points
        let fitted = fit_cubic_bezier_default(&samples).unwrap();

        // Get points from both segments
        let original_points = original.points();
        let fitted_points = fitted.points();

        // The fitted curve should have the same number of points
        assert_eq!(fitted_points.len(), original_points.len());

        // The first and last points should be close to the original
        let first_original = original_points[0];
        let first_fitted = fitted_points[0];
        let last_original = original_points[3];
        let last_fitted = fitted_points[3];

        // Check if the fitted curve has similar endpoints
        assert!(
            (first_original.distance(&first_fitted) < 0.5)
                || (first_original.distance(&last_fitted) < 0.5)
        );
        assert!(
            (last_original.distance(&last_fitted) < 0.5)
                || (last_original.distance(&first_fitted) < 0.5)
        );
    }

    #[test]
    fn test_demo_curve_fitting() {
        // Create a sample cubic bezier curve similar to the demo
        let segment = cubic!([(50.0, 200.0), (100.0, 50.0), (200.0, 50.0), (250.0, 200.0)]);

        // Sample points along the curve
        let samples = segment.sample_points(4);

        // Using the default chord length parameterization
        let fitted_segment1 = fit_cubic_bezier_default(&samples).unwrap();

        // Manually calculate t values
        let t_values = estimate_t_values_with_heuristic(&samples, THeuristic::default());

        // Fit with explicit t values
        let fitted_segment2 = least_square_solve_p_given_t(&samples, &t_values).unwrap();

        // Test that both methods produce the same result
        let points1 = fitted_segment1.points();
        let points2 = fitted_segment2.points();

        for i in 0..4 {
            assert!((points1[i].x - points2[i].x).abs() < 0.001);
            assert!((points1[i].y - points2[i].y).abs() < 0.001);
        }

        // Compare with original segment (endpoints should be preserved)
        let original_points = segment.points();
        let fitted_points = fitted_segment1.points();

        // The first and last control points should be very close to the original
        assert!(original_points[0].distance(&fitted_points[0]) < 0.1);
        assert!(original_points[3].distance(&fitted_points[3]) < 0.1);
    }

    #[test]
    fn test_lock_demo_fitted_points() {
        // Create the same cubic bezier curve as in the demo
        let original = cubic!([(50.0, 200.0), (100.0, 50.0), (200.0, 50.0), (250.0, 200.0)]);

        // Sample the same number of points as the demo (4)
        let samples = original.sample_points(4);

        // Fit using the default chord length parameterization
        let fitted = fit_cubic_bezier_default(&samples).unwrap();

        // Get the fitted control points
        let fitted_points = fitted.points();

        // Expected fitted control points from the demo output
        let expected_p1 = Point::new(50.0, 200.0);
        let expected_p2 = Point::new(38.61, 58.62);
        let expected_p3 = Point::new(261.39, 58.62);
        let expected_p4 = Point::new(250.0, 200.0);

        // Assert that the fitted points match the expected values with a small epsilon
        assert!((fitted_points[0].x - expected_p1.x).abs() < 0.01);
        assert!((fitted_points[0].y - expected_p1.y).abs() < 0.01);

        assert!((fitted_points[1].x - expected_p2.x).abs() < 0.01);
        assert!((fitted_points[1].y - expected_p2.y).abs() < 0.01);

        assert!((fitted_points[2].x - expected_p3.x).abs() < 0.01);
        assert!((fitted_points[2].y - expected_p3.y).abs() < 0.01);

        assert!((fitted_points[3].x - expected_p4.x).abs() < 0.01);
        assert!((fitted_points[3].y - expected_p4.y).abs() < 0.01);
    }

    #[test]
    fn test_different_parameterization() {
        // Create a bezier curve
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);

        // Sample points from the curve
        let samples = original.sample_points(10);

        // Standard chord length parameterization
        let t_values_chord = estimate_t_values_with_heuristic(&samples, THeuristic::default());
        let fitted_chord = least_square_solve_p_given_t(&samples, &t_values_chord).unwrap();

        // Try a different parameterization - uniform spacing
        let t_values_uniform: Vec<f64> = (0..samples.len())
            .map(|i| i as f64 / (samples.len() - 1) as f64)
            .collect();
        let fitted_uniform = least_square_solve_p_given_t(&samples, &t_values_uniform).unwrap();

        // Both should produce valid curves, but they may differ
        // At the very least, endpoints should be the same
        let chord_points = fitted_chord.points();
        let uniform_points = fitted_uniform.points();

        // Endpoints should be preserved in both parameterizations
        assert!(samples[0].distance(&chord_points[0]) < 0.1);
        assert!(samples[0].distance(&uniform_points[0]) < 0.1);
        assert!(samples.last().unwrap().distance(&chord_points[3]) < 0.1);
        assert!(samples.last().unwrap().distance(&uniform_points[3]) < 0.1);
    }
}
