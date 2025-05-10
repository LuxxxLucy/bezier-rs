//! This module provides convenient macros for creating points, segments, and curves.

/// Macro for creating a Point
#[macro_export]
macro_rules! pt {
    ($x:expr, $y:expr) => {
        $crate::data::Point::new($x as f64, $y as f64)
    };
}

/// Macro for creating a cubic bezier segment
#[macro_export]
macro_rules! cubic {
    ([$($point:expr),*]) => {{
        let points = [$($point),*];
        assert_eq!(points.len(), 4, "Cubic bezier requires exactly 4 points");
        $crate::data::BezierSegment::cubic(
            $crate::pt!(points[0].0, points[0].1),
            $crate::pt!(points[1].0, points[1].1),
            $crate::pt!(points[2].0, points[2].1),
            $crate::pt!(points[3].0, points[3].1),
        )
    }};
}

/// Macro for creating a quadratic bezier segment
#[macro_export]
macro_rules! quad {
    ([$($point:expr),*]) => {{
        let points = [$($point),*];
        assert_eq!(points.len(), 3, "Quadratic bezier requires exactly 3 points");
        $crate::data::BezierSegment::quadratic(
            $crate::pt!(points[0].0, points[0].1),
            $crate::pt!(points[1].0, points[1].1),
            $crate::pt!(points[2].0, points[2].1),
        )
    }};
}

/// Macro for creating a line segment
#[macro_export]
macro_rules! line {
    ([$($point:expr),*]) => {{
        let points = [$($point),*];
        assert_eq!(points.len(), 2, "Line segment requires exactly 2 points");
        $crate::data::BezierSegment::line(
            $crate::pt!(points[0].0, points[0].1),
            $crate::pt!(points[1].0, points[1].1),
        )
    }};
}

/// Macro for creating a Bezier curve from segments
#[macro_export]
macro_rules! curve {
    // Create from a list of segments
    ([$($segment:expr),*]) => {{
        let segments = vec![$($segment),*];
        $crate::data::BezierCurve::new(segments)
    }};

    // Create from an existing vector of segments
    ($segments:expr) => {
        $crate::data::BezierCurve::new($segments)
    };
}

/// Macro for creating a Bezier curve from a single segment
#[macro_export]
macro_rules! curve_from {
    ($segment:expr) => {
        $crate::data::BezierCurve::new(vec![$segment])
    };
}
