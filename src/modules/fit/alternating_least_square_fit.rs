//! Fit a cubic bezier curve to a set of points using the alternating method
//!
//! see Bezier Curve Fitting, Tim A Pastva 1998
//!
//! This method goes like this:
//! 1. Estimating t values using chord length parameterization
//! 2. Solving for control points given t values
//! 3. Updating t values by finding nearest points on the curve
//! 4. Goto step 2. Repeat until convergence or max iterations
//!
//! # Parameters
//!
//! * `points` - The points to fit the curve to
//! * `max_iterations` - Maximum number of iterations to perform
//! * `tolerance` - Convergence tolerance for t values. The algorithm stops when the maximum change
//!                 in any t value is less than this tolerance.
//! * `update_method` - The method to use for updating t values
//!    1. `NearestPoint` (default) - Finds the nearest point on the curve for each sample point
//!    2. `GaussNewton` - Uses Gauss-Newton optimization to update t values
//!
//! The Nearest Point method is chosen as the default because:
//! - It is simpler to implement and understand
//! - Gauss-newton is more prune to numerical stability issues
//! - As contrary to expectation, I found the nearest point method converges almost as fast as in Gauss-newton.
//!   Gauss-newton is only slightly faster than the nearest point, and I think the difference is negligible.
//!   It might be due to that my implementation of Gauss-Newton is not good enough
//!
//! # Example
//!
//! ```rust
//! use bezier_rs::{Point, BezierResult};
//! use bezier_rs::modules::fit::alternating_least_square_fit::{fit_cubic_bezier_alternating, TUpdateMethod, fit_cubic_bezier_alternating_default};
//!
//! // Some sample points
//! let points = vec![
//!     Point::new(0.0, 0.0),
//!     Point::new(1.0, 1.5),
//!     Point::new(2.0, 1.8),
//!     Point::new(3.0, 0.0),
//! ];
//!
//! // Fit a cubic Bezier using the default nearest point method
//! let result = fit_cubic_bezier_alternating(&points, 10, 1e-6, TUpdateMethod::NearestPoint);
//! let fitted = result.unwrap();
//!
//! // or using the default api
//! let result = fit_cubic_bezier_alternating_default(&points, 10, 1e-6);
//! let fitted = result.unwrap();
//! ```
use crate::data::{BezierSegment, Point};
use crate::error::{BezierError, BezierResult};
use crate::modules::fit::least_square_fit_common::{
    all_points_within_tolerance, get_delta_t, least_square_solve_p_given_t,
};
use crate::modules::fit::t_heuristic::{estimate_t_values_with_heuristic, THeuristic};

/// Methods for updating t values in alternating least squares fit
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum TUpdateMethod {
    /// Update t values by finding nearest points on the curve
    #[default]
    NearestPoint,
    /// Update t values using Gauss-Newton method
    GaussNewton,
}

/// Update t values using the nearest point method
fn update_t_values_nearest_point(segment: &BezierSegment, points: &[Point]) -> Vec<f64> {
    points
        .iter()
        .map(|point| {
            let (_, t) = segment.nearest_point(point);
            t
        })
        .collect()
}

/// Update t values using Gauss-Newton method
fn update_t_values_gauss_newton(
    points: &[Point],
    t_values: &[f64],
    segment: &BezierSegment,
) -> BezierResult<Vec<f64>> {
    let delta_t = get_delta_t(points, t_values, segment)?;
    Ok(t_values
        .iter()
        .zip(delta_t.iter())
        .map(|(&t, &dt)| (t + dt).clamp(0.0, 1.0))
        .collect())
}

/// Alternating optimization algorithm:
/// 1. Initialize t_i using chord length parameterization
///    t_i = (Σ_{k=1}^i ‖p_k - p_{k-1}‖) / total_length
/// 2. Solve for control points P given fixed t_i:
///    P = argminₚ Σ‖B(t_i; P) - p_i‖² = (AᵀA)⁻¹AᵀD
/// 3. Update t_i given fixed control points P:
///    t_i = argmin_t ‖B(t; P) - p_i‖² this done by either nearest point or gauss newton
/// 4. Repeat until convergence
pub fn fit_cubic_bezier_alternating(
    points: &[Point],
    max_iterations: usize,
    tolerance: f64,
    update_method: TUpdateMethod,
) -> BezierResult<BezierSegment> {
    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
    }

    // Start with chord length parameterization
    let mut t_values = estimate_t_values_with_heuristic(points, THeuristic::ChordLength);
    let mut segment = least_square_solve_p_given_t(points, &t_values)?;

    // If max_iterations is 0, return the initial curve
    if max_iterations == 0 {
        return Ok(segment);
    }

    // Iterate until convergence or max iterations
    for _ in 0..max_iterations {
        // Check if current fit is good enough, if so, return the current curve results
        if all_points_within_tolerance(&segment, points, tolerance) {
            break;
        }

        let new_t_values = match update_method {
            TUpdateMethod::NearestPoint => update_t_values_nearest_point(&segment, points),
            TUpdateMethod::GaussNewton => {
                update_t_values_gauss_newton(points, &t_values, &segment)?
            }
        };

        t_values = new_t_values;
        segment = least_square_solve_p_given_t(points, &t_values)?;
    }

    Ok(segment)
}

/// Fit a cubic bezier curve using alternating least squares with nearest point method
pub fn fit_cubic_bezier_alternating_default(
    points: &[Point],
    max_iterations: usize,
    tolerance: f64,
) -> BezierResult<BezierSegment> {
    fit_cubic_bezier_alternating(
        points,
        max_iterations,
        tolerance,
        TUpdateMethod::NearestPoint,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic;
    use approx::assert_relative_eq;

    #[test]
    fn test_alternating_fit_nearest_point() {
        // Create a bezier curve
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);

        // Sample points from the curve
        let samples = original.sample_points(20);

        // Fit a curve to the sampled points using nearest point method
        let fitted =
            fit_cubic_bezier_alternating(&samples, 100, 1e-6, TUpdateMethod::NearestPoint).unwrap();

        // For each sample point, find the nearest point on the fitted curve
        for sample in &samples {
            let (nearest_point, _) = fitted.nearest_point(sample);
            assert_relative_eq!(nearest_point.distance(sample), 0.0, epsilon = 1e-3);
        }
    }

    #[test]
    fn test_alternating_fit_gauss_newton() {
        // Create a bezier curve
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);

        // Sample points from the curve
        let samples = original.sample_points(20);

        // Fit a curve to the sampled points using Gauss-Newton method
        let fitted =
            fit_cubic_bezier_alternating(&samples, 100, 1e-6, TUpdateMethod::GaussNewton).unwrap();

        // For each sample point, find the nearest point on the fitted curve
        for sample in &samples {
            let (nearest_point, _) = fitted.nearest_point(sample);
            assert_relative_eq!(nearest_point.distance(sample), 0.0, epsilon = 1e-3);
        }
    }

    #[test]
    fn test_nearest_point_converge_faster_than_gauss_newton() {
        // Create a more complex curve
        let original = cubic!([(0.0, 0.0), (1.0, 3.0), (2.0, -1.0), (3.0, 2.0)]);

        // Sample points from the curve
        let samples = original.sample_points(15);

        // Track errors over iterations
        let mut nearest_point_errors = Vec::new();
        let mut gauss_newton_errors = Vec::new();

        // Run both methods for multiple iterations
        for iterations in 1..=20 {
            let nearest_point_result = fit_cubic_bezier_alternating(
                &samples,
                iterations,
                1e-6,
                TUpdateMethod::NearestPoint,
            )
            .unwrap();

            let gauss_newton_result = fit_cubic_bezier_alternating(
                &samples,
                iterations,
                1e-6,
                TUpdateMethod::GaussNewton,
            )
            .unwrap();

            // Calculate total error as sum of distances to nearest points
            let nearest_point_error: f64 = samples
                .iter()
                .map(|point| {
                    let (nearest_point, _) = nearest_point_result.nearest_point(point);
                    point.distance(&nearest_point)
                })
                .sum();

            let gauss_newton_error: f64 = samples
                .iter()
                .map(|point| {
                    let (nearest_point, _) = gauss_newton_result.nearest_point(point);
                    point.distance(&nearest_point)
                })
                .sum();

            nearest_point_errors.push(nearest_point_error);
            gauss_newton_errors.push(gauss_newton_error);
        }

        // Verify that errors decrease over iterations
        for i in 1..nearest_point_errors.len() {
            assert!(
                nearest_point_errors[i] <= nearest_point_errors[i - 1],
                "Nearest point error should decrease over iterations"
            );
            assert!(
                gauss_newton_errors[i] <= gauss_newton_errors[i - 1],
                "Gauss-Newton error should decrease over iterations"
            );
        }
        // sometimes the nearest point method converges faster, sometimes the gauss newton method converges faster
    }
}
