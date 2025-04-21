//! Fitting bezier curves to a set of points using the least squares method
//!
//! Based on "Least Squares Bezier Fit" by Jim Herold
//! https://web.archive.org/web/20180403213813/http://jimherold.com/2012/04/20/least-squares-bezier-fit/
//!
//! Altenatively see bezier primer chap 35 Curve fitting https://pomax.github.io/bezierinfo/#curvefitting
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
use crate::{cubic, pt}; // Import macros
use nalgebra::{DMatrix, DVector};

/// Fit a cubic bezier curve to a set of points using least squares with given t values
fn least_square_fit_routine(points: &[Point], t_values: &[f64]) -> BezierResult<BezierSegment> {
    if points.len() < 4 || points.len() != t_values.len() {
        return Err(BezierError::FitError(
            "At least 4 points are required and number of points must match number of t values"
                .to_string(),
        ));
    }

    // Create the polynomial matrix CT (n x 4)
    // Each row is [1, t, t², t³]
    let mut ct = DMatrix::zeros(points.len(), 4);
    for i in 0..points.len() {
        let t = t_values[i];
        ct[(i, 0)] = 1.0;
        ct[(i, 1)] = t;
        ct[(i, 2)] = t.powi(2);
        ct[(i, 3)] = t.powi(3);
    }

    // Create the Bernstein matrix B (4 x 4)
    let b = DMatrix::from_row_slice(
        4,
        4,
        &[
            1.0, 0.0, 0.0, 0.0, //  1  0  0  0
            -3.0, 3.0, 0.0, 0.0, // -3  3  0  0
            3.0, -6.0, 3.0, 0.0, //  3 -6  3  0
            -1.0, 3.0, -3.0, 1.0, // -1  3 -3  1
        ],
    );

    // Compose A = CT * B
    let a = &ct * &b;

    // Create vectors for x and y coordinates
    let b_x = DVector::from_iterator(points.len(), points.iter().map(|p| p.x));
    let b_y = DVector::from_iterator(points.len(), points.iter().map(|p| p.y));

    // Calculate the control points using least squares (A^T * A) * x = A^T * b
    let a_t = a.transpose();
    let a_ta = &a_t * &a;

    // Compute the inverse or pseudo-inverse of a_ta
    let a_ta_inv = match a_ta.try_inverse() {
        Some(inv) => inv,
        None => {
            return Err(BezierError::FitError(
                "Could not compute matrix inverse for least squares solution".to_string(),
            ))
        }
    };

    // Compute the solution: x = (A^T * A)^-1 * A^T * b
    let a_tb_x = &a_t * &b_x;
    let a_tb_y = &a_t * &b_y;

    let cx = a_ta_inv.clone() * a_tb_x;
    let cy = a_ta_inv * a_tb_y;

    // Create bezier segment from control points
    // Note: The ordering here needs to match the expected ordering for the curve
    // Let's check if we need to flip the control point ordering to match the original curve
    let p1 = pt!(cx[3], cy[3]); // First control point
    let p2 = pt!(cx[2], cy[2]); // Second control point
    let p3 = pt!(cx[1], cy[1]); // Third control point
    let p4 = pt!(cx[0], cy[0]); // Fourth control point

    // Create the default segment
    let segment = cubic!([(p1.x, p1.y), (p2.x, p2.y), (p3.x, p3.y), (p4.x, p4.y)]);

    // Compare with the original endpoints
    let start_point = points.first().unwrap();

    // Check if the orientation of the fitted curve matches the original points
    let start_dist_to_p1 = start_point.distance(&p1);
    let start_dist_to_p4 = start_point.distance(&p4);

    // If the start point is closer to p4 than to p1, flip the control points
    if start_dist_to_p4 < start_dist_to_p1 {
        return Ok(cubic!([
            (p4.x, p4.y),
            (p3.x, p3.y),
            (p2.x, p2.y),
            (p1.x, p1.y)
        ]));
    }

    Ok(segment)
}

/// Fit a cubic bezier curve to a set of points using least squares
///
/// This implementation uses the chord length parameterization for t-value estimation
/// as described in Jim Herold's blog post and the Bezier primer's Curve Fitting chapter.
pub fn fit_cubic_bezier_default(points: &[Point]) -> BezierResult<BezierSegment> {
    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
    }

    // Estimate t values using the chord length parameterization
    let t_value = estimate_t_values_with_heuristic(points, TParameterizationHeuristic::default());

    // Perform the least squares fitting with the estimated t values
    least_square_fit_routine(points, &t_value)
}

#[cfg(test)]
mod tests {
    use super::*;

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
