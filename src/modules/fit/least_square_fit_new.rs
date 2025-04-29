use crate::data::Point;
use crate::modules::fit::least_square_fit::least_square_solve_p_given_t;
use crate::modules::fit::t_heuristic::{estimate_t_values_with_heuristic, THeuristic};
use crate::{BezierError, BezierResult, BezierSegment};

/// Check if all points are within the given tolerance of the curve
fn all_points_within_tolerance(segment: &BezierSegment, points: &[Point], tolerance: f64) -> bool {
    points.iter().all(|p| {
        let (nearest, _) = segment.nearest_point(p);
        nearest.distance(p) <= tolerance
    })
}

/// Update t values by finding the nearest point on the curve for each input point
fn update_t_values_nearest_point(segment: &BezierSegment, points: &[Point]) -> Vec<f64> {
    points
        .iter()
        .map(|p| {
            let (_, t) = segment.nearest_point(p);
            t
        })
        .collect()
}

pub fn fit_cubic_bezier_new(
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

        let new_t_values = update_t_values_nearest_point(&segment, points);
        t_values = new_t_values;
        segment = least_square_solve_p_given_t(points, &t_values)?;
    }

    Ok(segment)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cubic;
    use approx::assert_relative_eq;

    #[test]
    fn test_argmin_gauss_newton() {
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 3.0)]);
        let samples = original.sample_points(20);

        let fitted = fit_cubic_bezier_new(&samples, 10, 0.001).unwrap();

        samples.iter().for_each(|p| {
            let (nearest, _) = fitted.nearest_point(p);
            assert_relative_eq!(nearest.distance(p), 0.0, epsilon = 0.01);
        });
    }
}
