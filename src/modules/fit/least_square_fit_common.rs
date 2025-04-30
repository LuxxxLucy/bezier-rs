use crate::data::{BezierSegment, Point};
use crate::error::{BezierError, BezierResult};
use crate::{cubic, pt}; // Import macros
use nalgebra::{DMatrix, DVector};

/// Compute the polynomial basis matrix CT (n x 4) of [1   t   t²   t³] for each t value
pub fn compute_polynomial_basis(t_values: &[f64]) -> DMatrix<f64> {
    let n = t_values.len();
    let mut ct = DMatrix::zeros(n, 4);
    for i in 0..n {
        let t = t_values[i];
        ct[(i, 0)] = 1.0;
        ct[(i, 1)] = t;
        ct[(i, 2)] = t.powi(2);
        ct[(i, 3)] = t.powi(3);
    }
    ct
}

/// Compute the derivative of polynomial basis matrix dCT/dt = [0  1   2t   3t²]
pub fn compute_polynomial_basis_derivative(t_values: &[f64]) -> DMatrix<f64> {
    let n = t_values.len();
    let mut ct_derivative = DMatrix::zeros(n, 4);
    for i in 0..n {
        let t = t_values[i];
        ct_derivative[(i, 0)] = 0.0;
        ct_derivative[(i, 1)] = 1.0;
        ct_derivative[(i, 2)] = 2.0 * t;
        ct_derivative[(i, 3)] = 3.0 * t.powi(2);
    }
    ct_derivative
}

/// Compute the Bernstein basis matrix B (4x4)
pub fn compute_bernstein_basis() -> DMatrix<f64> {
    DMatrix::from_row_slice(
        4,
        4,
        &[
            1.0, 0.0, 0.0, 0.0, // 1, 0, 0, 0
            -3.0, 3.0, 0.0, 0.0, // -3, 3, 0, 0
            3.0, -6.0, 3.0, 0.0, // 3, -6, 3, 0
            -1.0, 3.0, -3.0, 1.0, // -1, 3, -3, 1
        ],
    )
}

/// Compute the residual vector r = [B(t_i) - x_i, B(t_i) - y_i]ᵀ
pub fn compute_residual(
    points: &[Point],
    t_values: &[f64],
    segment: &BezierSegment,
) -> DVector<f64> {
    let n = points.len();
    let control_points = segment.points();
    let p_x: Vec<f64> = control_points.iter().map(|p| p.x).collect();
    let p_y: Vec<f64> = control_points.iter().map(|p| p.y).collect();
    let p_x = DVector::from_vec(p_x);
    let p_y = DVector::from_vec(p_y);

    // Compute polynomial basis and Bernstein matrix
    let ct = compute_polynomial_basis(t_values);
    let bernstein_matrix = compute_bernstein_basis();
    let a = ct * bernstein_matrix;

    // Compute residual vector
    let mut residual = DVector::zeros(2 * n);
    for i in 0..n {
        let predicted_x = (a.row(i) * &p_x)[0];
        let predicted_y = (a.row(i) * &p_y)[0];
        residual[2 * i] = predicted_x - points[i].x;
        residual[2 * i + 1] = predicted_y - points[i].y;
    }

    residual
}

/// Compute the Jacobian matrix J
pub fn compute_jacobian(
    points: &[Point],
    t_values: &[f64],
    segment: &BezierSegment,
) -> (DMatrix<f64>, DVector<f64>) {
    let n = points.len();
    let control_points = segment.points();
    let p_x: Vec<f64> = control_points.iter().map(|p| p.x).collect();
    let p_y: Vec<f64> = control_points.iter().map(|p| p.y).collect();
    let p_x = DVector::from_vec(p_x);
    let p_y = DVector::from_vec(p_y);

    // Compute polynomial basis derivative and Bernstein matrix
    let ct_derivative = compute_polynomial_basis_derivative(t_values);
    let bernstein_matrix = compute_bernstein_basis();

    // Compute derivative matrix
    let a_derivative = &ct_derivative * &bernstein_matrix;

    // Compute JᵀJ and Jᵀr directly without forming the full Jacobian
    let mut jtj = DMatrix::zeros(n, n);
    let mut jtr = DVector::zeros(n);

    for i in 0..n {
        let derivative_x = (a_derivative.row(i) * &p_x)[0];
        let derivative_y = (a_derivative.row(i) * &p_y)[0];

        // JᵀJ is diagonal since each t_i only affects its own residual
        jtj[(i, i)] = derivative_x.powi(2) + derivative_y.powi(2);

        // Jᵀr
        let residual = compute_residual(points, t_values, segment);
        jtr[i] = derivative_x * residual[2 * i] + derivative_y * residual[2 * i + 1];
    }

    (jtj, jtr)
}

/// Compute the step direction ΔT for updating t-values
pub fn get_delta_t(
    points: &[Point],
    t_values: &[f64],
    segment: &BezierSegment,
) -> BezierResult<Vec<f64>> {
    let (jtj, jtr) = compute_jacobian(points, t_values, segment);

    // Solve (JᵀJ)ΔT = -Jᵀr
    let delta_t = -jtj.lu().solve(&jtr).ok_or_else(|| {
        BezierError::FitError(
            "Failed to solve linear system for getting the Gauss-Newton Delta t".to_string(),
        )
    })?;

    Ok(delta_t.data.into())
}

/// Adjust and ensure the first t is always 0 and the last t is always 1
pub fn adjust_t_values(t_values: &[f64]) -> Vec<f64> {
    // Force first t-value to 0 and last t-value to 1
    let n = t_values.len();
    let mut adjusted_t_values = t_values.to_vec();
    adjusted_t_values[0] = 0.0;
    adjusted_t_values[n - 1] = 1.0;
    adjusted_t_values
}

/// Fit a cubic bezier curve to a set of points using least squares with given t values
pub fn least_square_solve_p_given_t(
    points: &[Point],
    t_values: &[f64],
) -> BezierResult<BezierSegment> {
    let n = points.len();

    if n < 4 || n != t_values.len() {
        return Err(BezierError::FitError(
            "At least 4 points are required and number of points must match number of t values"
                .to_string(),
        ));
    }

    // Create the polynomial matrix ct (n x 4) from `t_values` where each row is [1, t, t², t³]
    let ct = compute_polynomial_basis(t_values);

    // Create the Bernstein matrix B (4 x 4)
    let cubic_bernstein_matrix = compute_bernstein_basis();

    // Now we are going to convert the formulation into the simple minimize of ||Ax - b ||^2
    // where A = ct * cubic_bernstein_matrix
    // b is the input points: &[Point]
    let a = &ct * &cubic_bernstein_matrix;

    // Create vectors for x and y coordinates
    let b_x = DVector::from_iterator(n, points.iter().map(|p| p.x));
    let b_y = DVector::from_iterator(n, points.iter().map(|p| p.y));

    // Calculate the control points using least squares (A^T * A) * x = A^T * b
    let a_t = a.transpose();
    let a_ta = &a_t * &a;

    // Compute the inverse or pseudo-inverse of a_ta
    let a_ta_inv = match a_ta.try_inverse() {
        Some(inv) => inv,
        None => {
            return Err(BezierError::FitError(
                "Could not compute matrix inverse for least squares solve p given t".to_string(),
            ))
        }
    };

    // Compute the solution: x = (A^T * A)^-1 * A^T * b
    let a_tb_x = &a_t * &b_x;
    let a_tb_y = &a_t * &b_y;

    let cx = a_ta_inv.clone() * a_tb_x;
    let cy = a_ta_inv * a_tb_y;

    // Create bezier segment from control points
    // this is because we might have reversed the order, so we call `reorder_control_points`
    let (p1, p2, p3, p4) = reorder_control_points(
        pt!(cx[0], cy[0]),
        pt!(cx[1], cy[1]),
        pt!(cx[2], cy[2]),
        pt!(cx[3], cy[3]),
        points.first().unwrap(),
    );

    Ok(cubic!([
        (p1.x, p1.y),
        (p2.x, p2.y),
        (p3.x, p3.y),
        (p4.x, p4.y)
    ]))
}

/// Reorders control points to match the start point, used in `least_square_solve_p_given_t`
fn reorder_control_points(
    p1: Point,
    p2: Point,
    p3: Point,
    p4: Point,
    start_point: &Point,
) -> (Point, Point, Point, Point) {
    if start_point.distance(&p4) < start_point.distance(&p1) {
        return (p4, p3, p2, p1);
    }
    (p1, p2, p3, p4)
}

/// Check if all sample points are within tolerance distance of the fitted curve
pub fn all_points_within_tolerance(
    segment: &BezierSegment,
    points: &[Point],
    tolerance: f64,
) -> bool {
    points.iter().all(|point| {
        let (nearest_point, _) = segment.nearest_point(point);
        point.distance(&nearest_point) <= tolerance
    })
}
