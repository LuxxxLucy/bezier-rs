//! Data structure definitions for Bezier curves
//!
//! This module contains the fundamental data structures:
//!
//! 1. `Point`: representing a point in 2D space
//! 2. `BezierSegment`: a bezier segment, either Cubic or Quadratic
//! 3. `BezierCurve`: a bezier curve consisting of multiple segments
//! 4. `macros`: convenient macros for creating these structures
//!     - `pt!(x, y)`: Creates a Point
//!     - `cubic!([(x1, y1), (x2, y2), (x3, y3), (x4, y4)])`: Creates a cubic BezierSegment
//!     - `quad!([(x1, y1), (x2, y2), (x3, y3)])`: Creates a quadratic BezierSegment
//!     - `curve_from!(segment)`: Creates a BezierCurve from a single segment
//!     - `curve!(segments)`: Creates a BezierCurve from an existing vector of segments
//!     - `curve!([segment1, segment2, ...])`: Creates a BezierCurve from a list of segments

pub mod curve;
pub mod format;
pub mod macros;
pub mod point;
pub mod segment;

pub use curve::BezierCurve;
pub use point::Point;
pub use segment::BezierSegment;

#[doc(inline)]
pub use crate::{cubic, curve, curve_from, line, pt, quad};
