//! A Bezier curve: a collection of Bezier segments.

use crate::data::point::Point;
use crate::data::segment::BezierSegment;
use std::fmt;

/// A Bezier curve consisting of one or more Bezier segments
#[derive(Clone, PartialEq)] // we deliberately don't derive Debug
pub struct BezierCurve {
    /// The segments that make up this curve
    pub segments: Vec<BezierSegment>,
    /// Whether this curve is closed (end point connects to start point)
    is_closed: bool,
}

fn get_first_point(segments: &[BezierSegment]) -> Point {
    if segments.is_empty() {
        panic!("calling `get_first_point` on a bezier  empty list of segments");
    }
    segments[0].points()[0]
}

fn get_last_point(segments: &[BezierSegment]) -> Point {
    if segments.is_empty() {
        panic!("calling `get_last_point` on a bezier  empty list of segments");
    }
    if let Some(last_segment) = segments.last() {
        if let Some(end_point) = last_segment.points().last() {
            return *end_point;
        }
    }
    panic!("calling `get_last_point` on a bezier curve with no segments");
}

// Private helper to check if segments form a closed curve
fn is_segments_closed(segments: &[BezierSegment]) -> bool {
    if segments.is_empty() {
        panic!("calling `is_segments_closed` on an empty list of segments");
    }
    let start_point = get_first_point(segments);
    let end_point = get_last_point(segments);
    start_point == end_point
}

impl BezierCurve {
    /// Create a new curve from segments, automatically detecting if it's closed.
    /// Returns None if the segments list is empty.
    pub fn new(segments: Vec<BezierSegment>) -> Self {
        if segments.is_empty() {
            return Self {
                segments,
                is_closed: false,
            };
        }
        let is_closed = is_segments_closed(&segments);
        Self {
            segments,
            is_closed,
        }
    }

    /// Create a new closed curve from segments, returns None if:
    /// - The end point doesn't match the start point (for non-empty segments)
    pub fn new_closed(segments: Vec<BezierSegment>) -> Option<Self> {
        if segments.is_empty() {
            return None;
        }

        let mut segments = segments;
        if !is_segments_closed(&segments) {
            // add a new segment if line to from last point to the initial point
            let first_point = get_first_point(&segments);
            let last_point = get_last_point(&segments);
            segments.push(BezierSegment::Line {
                points: [last_point, first_point],
            });
        }

        Some(Self {
            segments,
            is_closed: true,
        })
    }

    /// Check if this curve is closed
    pub fn is_closed(&self) -> bool {
        self.is_closed
    }
}

impl fmt::Display for BezierCurve {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "BezierCurve [closed: {}]", self.is_closed)?;
        for (i, seg) in self.segments.iter().enumerate() {
            writeln!(f, "  {}: {}", i, seg)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quad;

    #[test]
    fn test_new_closed() {
        // Single segment with same start/end point can be closed
        let segment = quad!([(0, 0), (1, 1), (0, 0)]);
        let curve = BezierCurve::new(vec![segment]);
        assert!(!curve.segments.is_empty());
        assert!(curve.is_closed());

        // Segments that don't form a loop cannot be closed
        let segment1 = quad!([(0, 0), (1, 1), (2, 2)]);
        let segment2 = quad!([(2, 2), (3, 3), (4, 4)]);
        let curve = BezierCurve::new(vec![segment1, segment2]);
        assert!(!curve.segments.is_empty());
        assert!(!curve.is_closed());
    }

    #[test]
    fn test_new_auto_detect_closed() {
        // Single segment with same start/end point is detected as closed
        let segment = quad!([(0, 0), (1, 1), (0, 0)]);
        let curve = BezierCurve::new(vec![segment]);
        assert!(!curve.segments.is_empty());
        assert!(curve.is_closed());

        // Open curve is detected as open
        let segment = quad!([(0, 0), (1, 1), (2, 2)]);
        let curve = BezierCurve::new(vec![segment]);
        assert!(!curve.segments.is_empty());
        assert!(!curve.is_closed());

        // Multiple segments forming a loop are detected as closed
        let segments = vec![
            quad!([(0, 0), (1, 1), (2, 2)]),
            quad!([(2, 2), (1, 1), (0, 0)]),
        ];
        let curve = BezierCurve::new(segments);
        assert!(!curve.segments.is_empty());
        assert!(curve.is_closed());
    }
}
