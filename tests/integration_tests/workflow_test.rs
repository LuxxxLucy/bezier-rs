use bezier_rs::{
    modules::export, modules::fit::least_square_fit, BezierCurve, BezierSegment, Point,
};

#[test]
fn test_complete_workflow() {
    // 1. Create a bezier curve
    let original = BezierSegment::cubic(
        Point::new(10.0, 20.0),
        Point::new(30.0, 40.0),
        Point::new(50.0, 60.0),
        Point::new(70.0, 80.0),
    );

    // 2. Sample points from the curve
    let samples = original.sample_points(20);
    println!("Sampled {} points from the curve", samples.len());

    // 3. Fit a curve to the sampled points
    let fitted = least_square_fit::fit_cubic_bezier_default(&samples).unwrap();

    // 4. Export both to SVG
    let width = 300;
    let height = 300;
    let original_curve = BezierCurve::from_segment(&original);
    let fitted_curve = BezierCurve::from_segment(&fitted);

    let original_svg = export::svg::to_svg(&original_curve, width, height);
    let fitted_svg = export::svg::to_svg(&fitted_curve, width, height);
    let points_svg = export::svg::points_to_svg(&samples, width, height);

    // 5. Verify the results
    assert!(original_svg.contains("M10,20"));
    assert!(original_svg.contains("C30,40,50,60,70,80"));

    // Verify fitted curve creates valid SVG
    assert!(fitted_svg.contains("<path"));
    assert!(fitted_svg.contains("fill=\"none\""));

    // Verify points SVG contains circles
    assert!(points_svg.contains("<circle"));
    assert!(points_svg.contains("fill=\"red\""));

    // The workflow test is primarily to demonstrate the complete
    // capabilities of the library in one integrated test
    println!("Workflow test completed successfully!");
}
