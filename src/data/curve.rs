//! A Bezier curve: a collection of Bezier segments.

use crate::data::segment::BezierSegment;

/// A Bezier curve consisting of one or more Bezier segments
#[derive(Debug, Clone)]
pub struct BezierCurve {
    /// The segments that make up this curve
    pub segments: Vec<BezierSegment>,
    /// Whether this curve is closed (end point connects to start point)
    is_closed: bool,
}

impl BezierCurve {
    // Private helper to check if segments form a closed curve
    fn is_segments_closed(segments: &[BezierSegment]) -> bool {
        if segments.is_empty() {
            panic!("calling `is_segments_closed` on an empty list of segments");
        }
        let start_point = segments[0].points()[0];
        if let Some(last_segment) = segments.last() {
            if let Some(end_point) = last_segment.points().last() {
                return start_point == *end_point;
            }
        }
        false
    }

    /// Create a new curve from segments, automatically detecting if it's closed.
    /// Returns None if the segments list is empty.
    pub fn new(segments: Vec<BezierSegment>) -> Self {
        let is_closed = Self::is_segments_closed(&segments);
        Self {
            segments,
            is_closed,
        }
    }

    /// Create a new closed curve from segments, returns None if:
    /// - The segments list is empty
    /// - The end point doesn't match the start point
    pub fn new_closed(segments: Vec<BezierSegment>) -> Option<Self> {
        if segments.is_empty() {
            return None;
        }
        if !Self::is_segments_closed(&segments) {
            return None;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quad;

    #[test]
    #[should_panic]
    fn test_new_empty_segments() {
        let empty_curve = BezierCurve::new(vec![]);
        assert!(empty_curve.segments.is_empty());
        assert!(!empty_curve.is_closed());
    }

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
