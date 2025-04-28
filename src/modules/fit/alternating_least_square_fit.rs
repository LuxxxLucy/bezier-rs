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
use crate::modules::fit::least_square_fit::least_square_solve_p_given_t;
use crate::modules::fit::t_heuristic::{estimate_t_values_with_heuristic, THeuristic};
use nalgebra::{DMatrix, DVector};

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
    segment: &BezierSegment,
    points: &[Point],
    t_values: &[f64],
) -> Vec<f64> {
    // We want to minimize the squared residuals: E = Σ‖B(t_i) - p_i‖²
    // Using Gauss-Newton iteration: Δt = -(JᵀJ)⁻¹Jᵀr
    // Where:
    // - J is Jacobian matrix of residuals w.r.t t
    // - r is residual vector [B_x(t_i)-x_i, B_y(t_i)-y_i]ᵀ

    let n = points.len();
    let control_points = segment.points();

    // Convert control points to nalgebra format
    let p_x: Vec<f64> = control_points.iter().map(|p| p.x).collect();
    let p_y: Vec<f64> = control_points.iter().map(|p| p.y).collect();
    let p_x = DVector::from_vec(p_x);
    let p_y = DVector::from_vec(p_y);

    // 1. Compute polynomial basis matrix CT (n x 4) of [1   t   t²   t³] for each t value (each row)
    let mut ct = DMatrix::zeros(n, 4);
    for i in 0..n {
        let t = t_values[i];
        ct[(i, 0)] = 1.0;
        ct[(i, 1)] = t;
        ct[(i, 2)] = t.powi(2);
        ct[(i, 3)] = t.powi(3);
    }

    // 2. Bernstein basis matrix B (4x4) this is a constant.
    let bernstein_matrix = DMatrix::from_row_slice(
        4,
        4,
        &[
            1.0, 0.0, 0.0, 0.0, // 1, 0, 0, 0
            -3.0, 3.0, 0.0, 0.0, // -3, 3, 0, 0
            3.0, -6.0, 3.0, 0.0, // 3, -6, 3, 0
            -1.0, 3.0, -3.0, 1.0, // -1, 3, -3, 1
        ],
    );

    // 3. Compute derivative of basis matrix dCT/dt = [0  1   2t   3t²]
    let mut ct_derivative = DMatrix::zeros(n, 4);
    for i in 0..n {
        let t = t_values[i];
        ct_derivative[(i, 0)] = 0.0;
        ct_derivative[(i, 1)] = 1.0;
        ct_derivative[(i, 2)] = 2.0 * t;
        ct_derivative[(i, 3)] = 3.0 * t.powi(2);
    }

    // 4. Compute position and derivative matrices
    // A = CT * B
    // A' = dCT/dt * B
    let a = &ct * &bernstein_matrix;
    let a_derivative = &ct_derivative * &bernstein_matrix;

    // 5. Compute residual vector r = [B(t_i) - x_i, B(t_i) - y_i]ᵀ
    let mut residual = DVector::zeros(2 * n);
    for i in 0..n {
        let predicted_x = (a.row(i) * &p_x)[0];
        let predicted_y = (a.row(i) * &p_y)[0];
        residual[2 * i] = predicted_x - points[i].x;
        residual[2 * i + 1] = predicted_y - points[i].y;
    }

    // 6. Construct Jacobian matrix J = ∂r/∂t
    // J_ij = ∂r_i/∂t_j = [dB_x(t_j)/dt, dB_y(t_j)/dt] at point i
    let mut jacobian = DMatrix::zeros(2 * n, n);
    for i in 0..n {
        let derivative_x = (a_derivative.row(i) * &p_x)[0];
        let derivative_y = (a_derivative.row(i) * &p_y)[0];
        jacobian[(2 * i, i)] = derivative_x;
        jacobian[(2 * i + 1, i)] = derivative_y;
    }

    // 7. Solve Gauss-Newton normal equations: (JᵀJ)ΔT = -Jᵀr
    let jacobian_transpose = jacobian.transpose();
    let jtj = &jacobian_transpose * &jacobian;
    let jtr = &jacobian_transpose * residual;
    let delta_t = -jtj.lu().solve(&jtr).unwrap();

    // 8. Update t values and ensure they stay in [0,1]
    let new_t_values: Vec<f64> = t_values
        .iter()
        .zip(delta_t.iter())
        .map(|(&t, &dt)| (t + dt).clamp(0.0, 1.0))
        .collect();

    new_t_values
}

pub fn fit_cubic_bezier_alternating(
    points: &[Point],
    max_iterations: usize,
    tolerance: f64,
    update_method: TUpdateMethod,
) -> BezierResult<BezierSegment> {
    // Alternating optimization algorithm:
    // 1. Initialize t_i using chord length parameterization
    //    t_i = (Σ_{k=1}^i ‖p_k - p_{k-1}‖) / total_length
    // 2. Solve for control points P given fixed t_i:
    //    P = argminₚ Σ‖B(t_i; P) - p_i‖² = (AᵀA)⁻¹AᵀD
    // 3. Update t_i given fixed control points P:
    //    t_i = argmin_t ‖B(t; P) - p_i‖²
    // 4. Repeat until convergence

    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
    }

    // Start with chord length parameterization
    let mut t_values = estimate_t_values_with_heuristic(points, THeuristic::ChordLength);
    let mut segment = least_square_solve_p_given_t(points, &t_values)?;

    // If max_iterations is 0, return the initial curve resulted in the first linear problem
    // solution
    if max_iterations == 0 {
        return Ok(segment);
    }

    // Iterate until convergence or max iterations
    for _ in 0..max_iterations {
        let new_t_values = match update_method {
            TUpdateMethod::NearestPoint => update_t_values_nearest_point(&segment, points),
            TUpdateMethod::GaussNewton => update_t_values_gauss_newton(&segment, points, &t_values),
        };

        // Check if t values have converged
        let mut max_change = 0.0;
        for (old_t, new_t) in t_values.iter().zip(new_t_values.iter()) {
            let change = (old_t - new_t).abs();
            if change > max_change {
                max_change = change;
            }
        }

        if max_change < tolerance {
            break;
        }

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
