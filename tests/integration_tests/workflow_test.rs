use bezier_rs::modules::export::svg_path::ToSvgPath;
use bezier_rs::modules::fit::least_square_fit::fit_cubic_bezier_default;
use bezier_rs::{cubic, curve_from};

#[test]
fn test_complete_workflow() {
    // Create a simple cubic Bezier curve
    let original_segment = cubic!([
        (50.0, 200.0),  // Start point
        (100.0, 50.0),  // Control point 1
        (200.0, 50.0),  // Control point 2
        (250.0, 200.0)  // End point
    ]);
    let original_curve = curve_from!(original_segment.clone());

    // Sample points along the curve
    let points = original_segment.sample_points(10);

    // Fit a new curve to the sampled points
    let fitted_segment = fit_cubic_bezier_default(&points).unwrap();
    let fitted_curve = curve_from!(fitted_segment);

    // Convert both curves to SVG path data
    let original_path = original_curve.to_svg_path();
    let fitted_path = fitted_curve.to_svg_path();

    // Verify the exact path data
    assert_eq!(original_path, "M 50,200 C 100,50 200,50 250,200");

    // For the fitted curve, we can't predict the exact path data since it depends on the fitting algorithm
    // But we can verify it has the correct structure
    assert!(fitted_path.starts_with("M"));
    assert!(fitted_path.contains("C"));
}
