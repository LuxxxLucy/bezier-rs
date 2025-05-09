//! SVG export utilities for Bezier curves
//!
//! This module provides functionality to export Bezier curves and points
//! to SVG format for visualization and web display.
//!
//! # Features
//!
//! - Export a Bezier curve to an SVG path data string
//! - Export a collection of points to SVG path data
//!
//! # Examples
//!
//! ## Exporting a simple Bezier curve to SVG path data
//!
//! ```rust
//! use bezier_rs::{cubic, curve_from, modules::export::svg_path::ToSvgPath};
//!
//! // Create a simple cubic Bezier curve
//! let cubic_bezier = curve_from!(cubic!([
//!     (50.0, 200.0),    // Start point
//!     (100.0, 50.0),    // Control point 1
//!     (200.0, 50.0),    // Control point 2
//!     (250.0, 200.0)    // End point
//! ]));
//!
//! // Convert the curve to SVG path data
//! let path_data = cubic_bezier.to_svg_path();
//!
//! assert_eq!(path_data, "M50,200 C100,50,200,50,250,200");
//! ```
//!
//! ## Exporting multiple curve segments
//!
//! ```rust
//! use bezier_rs::{cubic, curve, quad, modules::export::svg_path::ToSvgPath};
//!
//! // Create a curve with multiple segments (a cubic followed by a quadratic)
//! let multi_segment = curve!([
//!     cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),  // Cubic segment
//!     quad!([(40, 50), (50, 60), (60, 70)])            // Quadratic segment
//! ]);
//!
//! // Convert to SVG path data
//! let path_data = multi_segment.to_svg_path();
//!
//! assert_eq!(path_data, "M10,20 C20,30,30,40,40,50 Q50,60,60,70");
//! ```

use crate::data::{BezierCurve, BezierSegment};

/// Trait for types that can be converted to SVG path data
pub trait ToSvgPath {
    /// Convert to SVG path data string
    fn to_svg_path(&self) -> String;
}

impl ToSvgPath for BezierCurve {
    fn to_svg_path(&self) -> String {
        if self.segments.is_empty() {
            return String::new();
        }

        let mut result = String::new();
        let mut first = true;

        for segment in &self.segments {
            match segment {
                BezierSegment::Cubic { points } => {
                    if first {
                        result.push_str(&format!("M{},{}", points[0].x, points[0].y));
                        first = false;
                    }
                    result.push_str(&format!(
                        " C{},{},{},{},{},{}",
                        points[1].x,
                        points[1].y,
                        points[2].x,
                        points[2].y,
                        points[3].x,
                        points[3].y
                    ));
                }
                BezierSegment::Quadratic { points } => {
                    if first {
                        result.push_str(&format!("M{},{}", points[0].x, points[0].y));
                        first = false;
                    }
                    result.push_str(&format!(
                        " Q{},{},{},{}",
                        points[1].x, points[1].y, points[2].x, points[2].y
                    ));
                }
            }
        }

        // Add closing command for closed curves
        if self.is_closed() {
            result.push('Z');
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{cubic, curve, curve_from, quad};

    #[test]
    fn test_curve_export_to_svg_path() {
        struct SvgPathExportTestCase<'a> {
            name: &'a str,
            curve: BezierCurve,
            expected_path: &'a str,
        }

        fn run_svg_path_export_test(test_case: SvgPathExportTestCase) {
            let path_data = test_case.curve.to_svg_path();
            assert_eq!(
                path_data, test_case.expected_path,
                "Test case: {}",
                test_case.name
            );
        }

        let test_cases = [
            SvgPathExportTestCase {
                name: "cubic_segment",
                curve: curve_from!(cubic!([(10, 20), (20, 30), (30, 40), (40, 50)])),
                expected_path: "M10,20 C20,30,30,40,40,50",
            },
            SvgPathExportTestCase {
                name: "quadratic_segment",
                curve: curve_from!(quad!([(10, 20), (40, 50), (70, 80)])),
                expected_path: "M10,20 Q40,50,70,80",
            },
            SvgPathExportTestCase {
                name: "multi_segment_curve",
                curve: curve!([
                    cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),
                    quad!([(40, 50), (50, 60), (60, 70)])
                ]),
                expected_path: "M10,20 C20,30,30,40,40,50 Q50,60,60,70",
            },
        ];

        // Run all test cases
        for test_case in test_cases {
            run_svg_path_export_test(test_case);
        }
    }

    #[test]
    fn test_closed_curve_export() {
        let segments = vec![
            cubic!([(10, 10), (20, 20), (40, 20), (50, 10)]).into(),
            quad!([(50, 10), (60, 0), (70, 10)]).into(),
            quad!([(70, 10), (60, 20), (10, 10)]).into(),
        ];
        let closed = BezierCurve::new_closed(segments).unwrap();
        let path_data = closed.to_svg_path();

        // The path data should contain both curves and a closing command
        assert!(path_data.contains("M10,10 C20,20,40,20,50,10"));
        assert!(path_data.contains("Q60,0,70,10"));
        assert!(path_data.contains("Q60,20,10,10"));
        assert!(path_data.ends_with("Z"));
    }

    #[test]
    fn test_round_trip_export_and_then_parse() {
        use crate::modules::parse::svg_path::FromSvgPath;

        let test_cases = [
            curve_from!(cubic!([(10, 20), (20, 30), (30, 40), (40, 50)])),
            curve_from!(quad!([(10, 20), (40, 50), (70, 80)])),
            curve!([
                cubic!([(10, 20), (20, 30), (30, 40), (40, 50)]),
                quad!([(40, 50), (50, 60), (60, 70)])
            ]),
        ];

        for (i, original_curve) in test_cases.iter().enumerate() {
            // Export to SVG path
            let path_data = original_curve.to_svg_path();

            // Parse back to curve
            let parsed_curve = BezierCurve::from_svg_path(&path_data)
                .unwrap_or_else(|e| panic!("Failed to parse path data for test case {}: {}", i, e));

            // Compare segments
            assert_eq!(
                original_curve.segments.len(),
                parsed_curve.segments.len(),
                "Segment count mismatch in test case {}",
                i
            );

            for (j, (original, parsed)) in original_curve
                .segments
                .iter()
                .zip(parsed_curve.segments.iter())
                .enumerate()
            {
                assert_eq!(
                    original, parsed,
                    "Segment {} mismatch in test case {}",
                    j, i
                );
            }

            // Compare closed state
            assert_eq!(
                original_curve.is_closed(),
                parsed_curve.is_closed(),
                "Closed state mismatch in test case {}",
                i
            );
        }
    }
}
