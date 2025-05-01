use crate::data::Point;
use crate::modules::fit::least_square_fit_common::{
    adjust_t_values, compute_residual, get_delta_t, least_square_solve_p_given_t,
};
use crate::modules::fit::t_heuristic::{estimate_t_values_with_heuristic, THeuristic};
use crate::{BezierError, BezierResult, BezierSegment};
use rand::prelude::*;
use rand_distr::Normal;

/// Parameters for gradient-based optimization
#[derive(Debug, Clone)]
pub struct GradientParams {
    pub min_step_size: f64,
    pub max_step_size: f64,
    pub num_random_samples: usize,
    pub random_scale: f64,
}

impl Default for GradientParams {
    fn default() -> Self {
        Self {
            min_step_size: 1e-6,
            max_step_size: 10.0,
            num_random_samples: 10,
            random_scale: 0.01,
        }
    }
}

/// Find the best step size using golden ratio line search
fn find_best_step_size_linesearch(
    points: &[Point],
    t_values: &[f64],
    delta_t: &[f64],
    min_step: f64,
    max_step: f64,
    num_steps: usize,
) -> f64 {
    const GOLDEN_RATIO: f64 = 0.618033988749895; // (sqrt(5) - 1) / 2

    let mut a = min_step;
    let mut b = max_step;
    let mut c = b - GOLDEN_RATIO * (b - a);
    let mut d = a + GOLDEN_RATIO * (b - a);

    // Compute initial function values
    let mut fc = compute_step_loss(points, t_values, delta_t, c);
    let mut fd = compute_step_loss(points, t_values, delta_t, d);

    // Golden ratio search
    for _ in 0..num_steps {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - GOLDEN_RATIO * (b - a);
            fc = compute_step_loss(points, t_values, delta_t, c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GOLDEN_RATIO * (b - a);
            fd = compute_step_loss(points, t_values, delta_t, d);
        }
    }

    // Return the midpoint of the final interval
    (a + b) / 2.0
}

/// Helper function to compute loss for a given step size
fn compute_step_loss(points: &[Point], t_values: &[f64], delta_t: &[f64], step_size: f64) -> f64 {
    let new_t_values: Vec<f64> = t_values
        .iter()
        .zip(delta_t.iter())
        .map(|(&t, &dt)| (t + step_size * dt).clamp(0.0, 1.0))
        .collect();

    compute_residual(
        points,
        &new_t_values,
        &least_square_solve_p_given_t(points, &new_t_values).unwrap(),
    )
    .norm()
}

/// Generate a list of t-value variations around the base point
fn generate_t_variations(
    base_t_values: &[f64],
    delta_t: &[f64],
    params: &GradientParams,
) -> Vec<Vec<f64>> {
    let normal = Normal::new(0.0, params.random_scale).unwrap();
    let mut rng = thread_rng();
    let mut variations = Vec::with_capacity(params.num_random_samples + 1);

    // Add the base point as the first element (index 0)
    variations.push(base_t_values.to_vec());

    // Generate random variations
    for _ in 0..params.num_random_samples {
        let random_t_values: Vec<f64> = base_t_values
            .iter()
            .zip(delta_t.iter())
            .map(|(&t, &dt)| {
                let random_factor = normal.sample(&mut rng);
                (t + random_factor * dt).clamp(0.0, 1.0)
            })
            .collect();
        // adjust the t values to ensure the first t is 0 and the last t is 1
        let adjusted_t_values = adjust_t_values(&random_t_values);
        variations.push(adjusted_t_values);
    }

    variations
}

/// Find the best t-values from a list of variations
fn find_best_t_values(points: &[Point], variations: &[Vec<f64>]) -> (Vec<f64>, f64) {
    let mut best_loss = f64::INFINITY;
    let mut best_t_values = variations[0].clone();

    for t_values in variations {
        let segment = least_square_solve_p_given_t(points, t_values).unwrap();
        let loss = compute_residual(points, t_values, &segment).norm();
        if loss < best_loss {
            best_loss = loss;
            best_t_values = t_values.clone();
        }
    }

    (best_t_values, best_loss)
}

/// Update t values using weak variable projection with line search and variations
fn update_t_values_weak_varpro(
    points: &[Point],
    t_values: &[f64],
    segment: &BezierSegment,
    params: &GradientParams,
) -> BezierResult<(Vec<f64>, f64)> {
    // Get step direction using Gauss-Newton on the temporary loss
    let delta_t = get_delta_t(points, t_values, segment)?;

    // Find the best step size using golden ratio line search
    let best_step_size = find_best_step_size_linesearch(
        points,
        t_values,
        &delta_t,
        params.min_step_size,
        params.max_step_size,
        10, // Try 10 different step sizes
    );

    // Use the best step size found as the base point
    let base_t_values: Vec<f64> = t_values
        .iter()
        .zip(delta_t.iter())
        .map(|(&t, &dt)| (t + best_step_size * dt).clamp(0.0, 1.0))
        .collect();
    let base_t_values = adjust_t_values(&base_t_values);

    // Generate variations and append the base point
    let mut variations = generate_t_variations(&base_t_values, &delta_t, params);
    variations.push(base_t_values);

    // Find the best t-values from all variations
    let (best_t_values, best_loss) = find_best_t_values(points, &variations);

    // Compute original loss
    let original_loss = compute_residual(points, t_values, segment).norm();

    // Only update if the new loss is better
    if best_loss < original_loss {
        Ok((best_t_values, best_loss))
    } else {
        Ok((t_values.to_vec(), original_loss))
    }
}

/// Check if all points are within the given tolerance of the curve
fn all_points_within_tolerance(segment: &BezierSegment, points: &[Point], tolerance: f64) -> bool {
    points.iter().all(|p| {
        let (nearest, _) = segment.nearest_point(p);
        nearest.distance(p) <= tolerance
    })
}

pub fn fit_cubic_bezier_weak_varpro(
    points: &[Point],
    max_iterations: usize,
    tolerance: f64,
    gradient_params: Option<GradientParams>,
) -> BezierResult<BezierSegment> {
    if points.len() < 4 {
        return Err(BezierError::FitError(
            "At least 4 points are required for cubic bezier fitting".to_string(),
        ));
    }

    let params = gradient_params.unwrap_or_default();

    // Start with chord length parameterization
    let mut t_values = estimate_t_values_with_heuristic(points, THeuristic::ChordLength);
    let mut segment = least_square_solve_p_given_t(points, &t_values)?;
    let mut prev_loss = compute_residual(points, &t_values, &segment).norm();

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

        // Update t-values using weak variable projection with line search and variations
        let (new_t_values, new_loss) =
            update_t_values_weak_varpro(points, &t_values, &segment, &params)?;

        // Check if loss improvement is too small
        if prev_loss < new_loss {
            break;
        }

        t_values = new_t_values;
        prev_loss = new_loss;
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
    fn test_gradient_fit() {
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 3.0)]);
        let samples = original.sample_points(20);

        let fitted = fit_cubic_bezier_weak_varpro(&samples, 10, 0.001, None).unwrap();

        samples.iter().for_each(|p| {
            let (nearest, _) = fitted.nearest_point(p);
            assert_relative_eq!(nearest.distance(p), 0.0, epsilon = 0.02);
        });
    }

    #[test]
    fn test_gradient_fit_with_params() {
        let original = cubic!([(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 3.0)]);
        let samples = original.sample_points(20);

        let params = GradientParams {
            min_step_size: 1e-8,
            max_step_size: 10.0,
            num_random_samples: 10,
            random_scale: 0.1,
        };

        let fitted = fit_cubic_bezier_weak_varpro(&samples, 10, 0.001, Some(params)).unwrap();

        samples.iter().for_each(|p| {
            let (nearest, _) = fitted.nearest_point(p);
            assert_relative_eq!(nearest.distance(p), 0.0, epsilon = 0.02);
        });
    }
}
