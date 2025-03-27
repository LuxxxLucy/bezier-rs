//! SVG export utilities for Bezier curves
//!
//! This module provides functionality to export Bezier curves and points
//! to SVG format for visualization and web display.
//!
//! # Features
//!
//! - Export a Bezier curve to an SVG path
//! - Export a collection of points to SVG circles
//! - Configure the output SVG dimensions
//!
//! # Examples
//!
//! ## Exporting a simple Bezier curve to SVG
//!
//! ```rust
//! use bezier_rs::{cubic, curve_from, modules::export::svg};
//!
//! // Create a simple cubic Bezier curve
//! let cubic_bezier = curve_from!(cubic!([
//!     (50.0, 200.0),    // Start point
//!     (100.0, 50.0),    // Control point 1
//!     (200.0, 50.0),    // Control point 2
//!     (250.0, 200.0)    // End point
//! ]));
//!
//! // Convert the curve to SVG with a 300x300 canvas
//! let svg_string = svg::to_svg(&cubic_bezier, 300, 300);
//!
//! assert!(svg_string.contains("C100,50,200,50,250,200"));
//! ```
//!
//! ## Exporting multiple curve segments
//!
//! ```rust
//! use bezier_rs::{cubic, curve, quad, modules::export::svg};
//!
//! // Create a curve with multiple segments (a cubic followed by a quadratic)
//! let multi_segment = curve!([
//!     cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),  // Cubic segment
//!     quad!([(40, 50), (50, 60), (60, 70)])            // Quadratic segment
//! ]);
//!
//! // Convert to SVG
//! let svg_string = svg::to_svg(&multi_segment, 100, 100);
//!
//! // The SVG now contains a path with both cubic and quadratic segments
//! assert!(svg_string.contains("C20,30,30,40,40,50"));
//! assert!(svg_string.contains("Q50,60,60,70"));
//! ```
//!
//! ## Exporting points for visualization
//!
//! ```rust
//! use bezier_rs::{cubic, modules::export::svg, pt};
//!
//! // Create a bezier curve
//! let curve = cubic!([(10, 20), (30, 40), (50, 60), (70, 80)]);
//!
//! // Sample points along the curve
//! let points = curve.sample_points(10);
//!
//! // Export the points to SVG with red circles
//! let svg_string = svg::points_to_svg(&points, 100, 100);
//!
//! // The SVG now contains 10 circle elements representing the points
//! assert!(svg_string.contains("<circle"));
//! assert!(svg_string.contains("fill=\"red\""));
//! ```

use crate::data::{BezierCurve, BezierSegment, Point};
use svg::node::element::{path::Data, Path};
use svg::Document;

/// Export a bezier curve to an SVG string
///
/// This function converts a `BezierCurve` to an SVG string representation,
/// containing a path element that follows the curve segments.
///
/// # Parameters
///
/// * `curve` - The bezier curve to export
/// * `width` - The width of the SVG canvas
/// * `height` - The height of the SVG canvas
///
/// # Returns
///
/// A string containing the SVG representation of the curve
///
/// # Example
///
/// ```
/// use bezier_rs::{cubic, curve_from, modules::export::svg};
///
/// // Create a simple cubic Bezier curve
/// let cubic_bezier = curve_from!(cubic!([
///     (50.0, 200.0),    // Start point
///     (100.0, 50.0),    // Control point 1
///     (200.0, 50.0),    // Control point 2
///     (250.0, 200.0)    // End point
/// ]));
///
/// // Export to SVG
/// let svg_string = svg::to_svg(&cubic_bezier, 300, 300);
///
/// // The SVG string contains a path element with the curve
/// assert!(svg_string.contains("<path"));
/// ```
pub fn to_svg(curve: &BezierCurve, width: u32, height: u32) -> String {
    let mut document = Document::new()
        .set("width", width)
        .set("height", height)
        .set("viewBox", (0, 0, width, height));

    if curve.segments.is_empty() {
        return document.to_string();
    }

    // Get the first point from the first segment
    let first_points = curve.segments[0].points();
    let mut data = Data::new().move_to((first_points[0].x, first_points[0].y));

    for segment in &curve.segments {
        match segment {
            BezierSegment::Cubic { points } => {
                data = data.cubic_curve_to((
                    points[1].x,
                    points[1].y,
                    points[2].x,
                    points[2].y,
                    points[3].x,
                    points[3].y,
                ));
            }
            BezierSegment::Quadratic { points } => {
                data =
                    data.quadratic_curve_to((points[1].x, points[1].y, points[2].x, points[2].y));
            }
        }
    }

    let path = Path::new()
        .set("fill", "none")
        .set("stroke", "black")
        .set("stroke-width", 1)
        .set("d", data);

    document = document.add(path);

    document.to_string()
}

/// Export sampled points from a curve to an SVG string
///
/// This function converts a collection of points to an SVG string representation,
/// with each point rendered as a small red circle.
///
/// # Parameters
///
/// * `points` - The collection of points to export
/// * `width` - The width of the SVG canvas
/// * `height` - The height of the SVG canvas
///
/// # Returns
///
/// A string containing the SVG representation of the points
///
/// # Example
///
/// ```
/// use bezier_rs::{cubic, modules::export::svg, pt};
///
/// // Create some sample points
/// let points = vec![
///     pt!(10.0, 20.0),
///     pt!(30.0, 50.0),
///     pt!(50.0, 30.0),
///     pt!(70.0, 60.0)
/// ];
///
/// // Export to SVG
/// let svg_string = svg::points_to_svg(&points, 100, 100);
///
/// // The resulting SVG contains circle elements for each point
/// assert!(svg_string.contains("<circle cx=\"10\""));
/// assert!(svg_string.contains("<circle cx=\"30\""));
/// ```
pub fn points_to_svg(points: &[Point], width: u32, height: u32) -> String {
    let mut document = Document::new()
        .set("width", width)
        .set("height", height)
        .set("viewBox", (0, 0, width, height));

    // Create a path with circles for each point
    for point in points {
        let circle = svg::node::element::Circle::new()
            .set("cx", point.x)
            .set("cy", point.y)
            .set("r", 2)
            .set("fill", "red");

        document = document.add(circle);
    }

    document.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{cubic, curve, curve_from, quad};

    #[test]
    fn test_points_export_to_svg() {
        // Sample points from a cubic segment
        let samples = cubic!([(10, 20), (30, 40), (50, 60), (70, 80)]).sample_points(10);

        // Expected SVG string with exact points
        let expected_svg = "<svg height=\"100\" viewBox=\"0 0 100 100\" width=\"100\" xmlns=\"http://www.w3.org/2000/svg\">\n\
            <circle cx=\"10\" cy=\"20\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"16.666666666666668\" cy=\"26.666666666666664\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"23.333333333333332\" cy=\"33.333333333333336\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"30\" cy=\"40\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"36.666666666666664\" cy=\"46.666666666666664\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"43.33333333333334\" cy=\"53.333333333333336\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"50\" cy=\"59.99999999999999\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"56.66666666666667\" cy=\"66.66666666666667\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"63.333333333333336\" cy=\"73.33333333333334\" fill=\"red\" r=\"2\"/>\n\
            <circle cx=\"70\" cy=\"80\" fill=\"red\" r=\"2\"/>\n\
            </svg>";

        // Export sampled points to SVG and compare
        let svg_string = points_to_svg(&samples, 100, 100);
        assert_eq!(svg_string, expected_svg);
    }

    #[test]
    fn test_curve_export_to_svg() {
        struct SvgExportTestCase<'a> {
            name: &'a str,
            curve: BezierCurve,
            expected_svg: &'a str,
        }

        fn run_svg_export_test(test_case: SvgExportTestCase) {
            let svg_string = to_svg(&test_case.curve, 100, 100);
            assert_eq!(
                svg_string, test_case.expected_svg,
                "Test case: {}",
                test_case.name
            );
        }

        let test_cases = [
            SvgExportTestCase {
                name: "cubic_segment",
                curve: curve_from!(cubic!([(10, 20), (20, 30), (30, 40), (40, 50)])),
                expected_svg: "<svg height=\"100\" viewBox=\"0 0 100 100\" width=\"100\" xmlns=\"http://www.w3.org/2000/svg\">\n\
                    <path d=\"M10,20 C20,30,30,40,40,50\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n\
                    </svg>",
            },
            SvgExportTestCase {
                name: "quadratic_segment",
                curve: curve_from!(quad!([(10, 20), (40, 50), (70, 80)])),
                expected_svg: "<svg height=\"100\" viewBox=\"0 0 100 100\" width=\"100\" xmlns=\"http://www.w3.org/2000/svg\">\n\
                    <path d=\"M10,20 Q40,50,70,80\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n\
                    </svg>",
            },
            SvgExportTestCase {
                name: "multi_segment_curve",
                curve: curve!([
                    cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),
                    quad!([(40, 50), (50, 60), (60, 70)])
                ]),
                expected_svg: "<svg height=\"100\" viewBox=\"0 0 100 100\" width=\"100\" xmlns=\"http://www.w3.org/2000/svg\">\n\
                    <path d=\"M10,20 C20,30,30,40,40,50 Q50,60,60,70\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n\
                    </svg>",
            },
            SvgExportTestCase {
                name: "closed_shape",
                curve: curve!([
                    cubic!([(20, 20), (40, 20), (60, 40), (60, 60)]),
                    cubic!([(60, 60), (60, 80), (40, 80), (20, 60)]),
                    cubic!([(20, 60), (0, 40), (0, 20), (20, 20)])
                ]),
                expected_svg: "<svg height=\"100\" viewBox=\"0 0 100 100\" width=\"100\" xmlns=\"http://www.w3.org/2000/svg\">\n\
                    <path d=\"M20,20 C40,20,60,40,60,60 C60,80,40,80,20,60 C0,40,0,20,20,20\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n\
                    </svg>",
            },
        ];

        // Run all test cases
        for test_case in test_cases {
            run_svg_export_test(test_case);
        }
    }

    #[test]
    fn test_svg_consistency() {
        // Two different ways to create the same curve
        // Both approaches should produce identical SVG

        let curve1 = &curve!([
            cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),
            quad!([(40, 50), (50, 60), (60, 70)])
        ]);

        let segments = vec![
            cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),
            quad!([(40, 50), (50, 60), (60, 70)]),
        ];
        let curve2 = &curve!(segments);

        let svg1 = to_svg(curve1, 100, 100);
        let svg2 = to_svg(curve2, 100, 100);
        assert_eq!(svg1, svg2);
    }
}
