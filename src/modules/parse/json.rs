//! Parsing JSON into Bezier Curves
//!
//! The expected JSON format is an array of point objects:
//!
//! ```json
//! [
//!   {"x": 0.0, "y": 0.0, "on": true},   // Starting point (on-curve)
//!   {"x": 1.0, "y": 1.0, "on": false},  // Control point (off-curve)
//!   {"x": 2.0, "y": 1.0, "on": false},  // Control point (off-curve)
//!   {"x": 3.0, "y": 0.0, "on": true}    // End point (on-curve)
//! ]
//! ```
//!
//! Where:
//! - `x` and `y` are the coordinates of the point (required)
//! - `on` is a boolean indicating if the point is on the curve (optional, defaults to `true`)
//!
//! # Curve Interpretation
//!
//! The algorithm creates Bezier segments based on the following rules:
//!
//! - A sequence of "on-curve point → off-curve point → on-curve point" creates a quadratic Bezier
//! - A sequence of "on-curve point → off-curve point → off-curve point → on-curve point" creates a cubic Bezier
//! - A sequence of "on-curve point → on-curve point" creates a straight line (implemented as a quadratic with control point at midpoint)
//!
//! # Example
//!
//! ```rust
//! use bezier_rs::modules::parse::json;
//!
//! let json_str = r#"[
//!     {"x": 0.0, "y": 0.0, "on": true},
//!     {"x": 1.0, "y": 1.0, "on": false},
//!     {"x": 2.0, "y": 0.0, "on": true}
//! ]"#;
//!
//! let curve = json::parse(json_str).unwrap();
//! println!("Parsed a curve with {} segments", curve.segments.len());
//! ```

use crate::data::{BezierCurve, BezierSegment, Point};
use crate::error::{BezierError, BezierResult};
use crate::{cubic, curve, quad};
use serde::{Deserialize, Serialize};

/// Information about a point on the curve - on or off the curve
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
struct PointInfo {
    pub point: Point,
    pub on_curve: bool,
}

impl PointInfo {
    fn new(x: f64, y: f64, on_curve: bool) -> Self {
        Self {
            point: Point::new(x, y),
            on_curve,
        }
    }
}

/// Parsed JSON format for a point with on/off curve information
#[derive(Debug, Clone, Deserialize, Serialize)]
struct JsonPointInfo {
    x: f64,
    y: f64,
    #[serde(rename = "on", default = "default_on_curve")]
    on_curve: bool,
}

fn default_on_curve() -> bool {
    true
}

/// Parse a JSON string into a BezierCurve
/// Format expected: [{"x": x1, "y": y1, "on": true/false}, {...}, ...]
///
/// # Examples
///
/// ```
/// use bezier_rs::modules::parse::json;
///
/// // A simple quadratic Bezier curve
/// let json_str = r#"[
///     {"x": 0.0, "y": 0.0, "on": true},
///     {"x": 1.0, "y": 1.0, "on": false},
///     {"x": 2.0, "y": 0.0, "on": true}
/// ]"#;
///
/// let curve = json::parse(json_str).unwrap();
/// assert_eq!(curve.segments.len(), 1); // One quadratic segment
///
/// // A cubic Bezier curve
/// let cubic_json = r#"[
///     {"x": 0.0, "y": 0.0, "on": true},
///     {"x": 1.0, "y": 1.0, "on": false},
///     {"x": 2.0, "y": 1.0, "on": false},
///     {"x": 3.0, "y": 0.0, "on": true}
/// ]"#;
///
/// let cubic_curve = json::parse(cubic_json).unwrap();
/// assert_eq!(cubic_curve.segments.len(), 1); // One cubic segment
/// ```
pub fn parse(json_str: &str) -> BezierResult<BezierCurve> {
    let points: Vec<JsonPointInfo> = serde_json::from_str(json_str)
        .map_err(|e| BezierError::ParseError(format!("JSON parse error: {}", e)))?;

    if points.is_empty() {
        return Err(BezierError::ParseError("Empty points array".to_string()));
    }

    // Convert to PointInfo structs
    let point_infos: Vec<PointInfo> = points
        .iter()
        .map(|p| PointInfo::new(p.x, p.y, p.on_curve))
        .collect();

    // Now convert the points to segments
    let segments = create_segments_from_points(&point_infos)?;

    // Create and return the BezierCurve
    Ok(curve!(segments))
}

/// Create bezier segments from a list of points with on/off curve information
fn create_segments_from_points(points: &[PointInfo]) -> BezierResult<Vec<BezierSegment>> {
    let mut segments = Vec::new();
    let mut i = 0;

    while i < points.len() {
        // Each segment starts with an on-curve point
        if !points[i].on_curve {
            return Err(BezierError::ParseError(format!(
                "Expected on-curve point at index {}",
                i
            )));
        }

        let start_point = points[i].point;
        i += 1;

        // If we're at the end, we're done
        if i >= points.len() {
            break;
        }

        // Next point is off-curve (control point)
        if !points[i].on_curve {
            let control1 = points[i].point;
            i += 1;

            if i >= points.len() {
                return Err(BezierError::ParseError(
                    "Curve cannot end with an off-curve point".to_string(),
                ));
            }

            // Next point can be on or off curve
            if points[i].on_curve {
                // Quadratic bezier: start(on) -> control(off) -> end(on)
                let end_point = points[i].point;
                segments.push(quad!([
                    (start_point.x, start_point.y),
                    (control1.x, control1.y),
                    (end_point.x, end_point.y)
                ]));

                // Stay at current point for next segment
                continue;
            } else {
                // Two off-curve points in a row means cubic bezier
                let control2 = points[i].point;
                i += 1;

                if i >= points.len() || !points[i].on_curve {
                    return Err(BezierError::ParseError(
                        "Expected an on-curve point after two off-curve points".to_string(),
                    ));
                }

                // Cubic bezier: start(on) -> control1(off) -> control2(off) -> end(on)
                let end_point = points[i].point;
                segments.push(cubic!([
                    (start_point.x, start_point.y),
                    (control1.x, control1.y),
                    (control2.x, control2.y),
                    (end_point.x, end_point.y)
                ]));

                // Stay at current point for next segment
                continue;
            }
        } else {
            // Two on-curve points in a row means a straight line
            // We'll represent it as a quadratic bezier with control point at midpoint
            let end_point = points[i].point;
            let mid_x = (start_point.x + end_point.x) / 2.0;
            let mid_y = (start_point.y + end_point.y) / 2.0;

            segments.push(quad!([
                (start_point.x, start_point.y),
                (mid_x, mid_y),
                (end_point.x, end_point.y)
            ]));

            // Stay at current point for next segment
            continue;
        }
    }

    Ok(segments)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pt;

    #[test]
    fn test_json_parsing() {
        let json = r#"[
            {"x": 0.0, "y": 0.0, "on": true},
            {"x": 1.0, "y": 1.0, "on": false},
            {"x": 2.0, "y": 1.0, "on": false},
            {"x": 3.0, "y": 0.0, "on": true}
        ]"#;

        let curve = parse(json).unwrap();
        assert_eq!(curve.segments.len(), 1);

        let segment = &curve.segments[0];

        // Check that it's a cubic segment by pattern matching
        match segment {
            BezierSegment::Cubic { .. } => (),
            _ => panic!("Expected a cubic segment"),
        }

        // Check the points
        let points = segment.points();
        assert_eq!(points[0], pt!(0.0, 0.0));
        assert_eq!(points[3], pt!(3.0, 0.0));
    }
}
