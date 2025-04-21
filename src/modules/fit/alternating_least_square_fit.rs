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
//! # Example
//!
//! ```rust
//! use bezier_rs::{modules::fit::alternating_least_square_fit, Point, BezierResult};
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
//! let result: BezierResult<_> = alternating_least_square_fit::fit_cubic_bezier_alternating(&points, 10, 1e-6);
//! let fitted = result.unwrap();
//! ```
use crate::data::{BezierSegment, Point};
use crate::error::{BezierError, BezierResult};
use crate::modules::fit::least_square_fit::least_square_solve_p_given_t;
use crate::modules::fit::t_heuristic::{estimate_t_values_with_heuristic, THeuristic};

pub fn fit_cubic_bezier_alternating(
    points: &[Point],
    max_iterations: usize,
    tolerance: f64,
) -> BezierResult<BezierSegment> {
    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
    }

    // Start with chord length parameterization
    let mut t_values = estimate_t_values_with_heuristic(points, THeuristic::ChordLength);
    let mut segment = least_square_solve_p_given_t(points, &t_values)?;

    // Iterate until convergence or max iterations
    for _ in 0..max_iterations {
        // Update t values by finding nearest points
        let mut new_t_values = Vec::with_capacity(points.len());
        for point in points {
            let (_, t) = segment.nearest_point(point);
            new_t_values.push(t);
        }

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic;

    #[test]
    fn test_alternating_fit() {
        // Create a bezier curve
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);

        // Sample points from the curve
        let samples = original.sample_points(20);

        // Fit a curve to the sampled points
        let fitted = fit_cubic_bezier_alternating(&samples, 10, 1e-6).unwrap();

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
}
