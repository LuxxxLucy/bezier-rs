//! A Bezier curve: a collection of Bezier segments.

use crate::data::segment::BezierSegment;

/// A Bezier curve consisting of one or more Bezier segments
#[derive(Debug, Clone)]
pub struct BezierCurve {
    pub segments: Vec<BezierSegment>,
}

impl BezierCurve {
    /// Create a new bezier curve from a list of segments
    pub fn new(segments: Vec<BezierSegment>) -> Self {
        Self { segments }
    }

    /// Create a bezier curve from a single segment
    pub fn from_segment(segment: &BezierSegment) -> Self {
        Self {
            segments: vec![segment.clone()],
        }
    }
}
