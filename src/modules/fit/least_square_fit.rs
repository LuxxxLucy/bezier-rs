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
//! let result: BezierResult<_> = least_square_fit::fit_cubic_bezier(&points);
//! let fitted = result.unwrap();
//! ```

use crate::data::{BezierSegment, Point};
use crate::error::{BezierError, BezierResult};
use crate::{cubic, pt}; // Import macros
use nalgebra::{DMatrix, DVector};

/// Fit a cubic bezier curve to a set of points using least squares
pub fn fit_cubic_bezier(points: &[Point]) -> BezierResult<BezierSegment> {
    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
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
    let t_values: Vec<f64> = path_lengths
        .iter()
        .map(|&length| length / total_length)
        .collect();

    // For each point, calculate the basis function values
    let mut a = DMatrix::zeros(points.len(), 4);
    for i in 0..points.len() {
        let t = t_values[i];
        a[(i, 0)] = (1.0 - t).powi(3);
        a[(i, 1)] = 3.0 * t * (1.0 - t).powi(2);
        a[(i, 2)] = 3.0 * t.powi(2) * (1.0 - t);
        a[(i, 3)] = t.powi(3);
    }

    // Create vectors for x and y coordinates
    let b_x = DVector::from_iterator(points.len(), points.iter().map(|p| p.x));
    let b_y = DVector::from_iterator(points.len(), points.iter().map(|p| p.y));

    // Calculate the control points using least squares (A^T * A) * x = A^T * b
    let a_t = a.transpose();
    let a_ta = &a_t * &a;

    // Compute the inverse or pseudo-inverse of a_ta
    // For simplicity, we'll use a direct matrix inverse here
    // In a production system, you might want to use SVD or another method
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
        let fitted = fit_cubic_bezier(&samples).unwrap();

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
