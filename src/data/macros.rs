//! This module provides convenient macros for creating points, segments, and curves.

/// Macro for creating a Point
#[macro_export]
macro_rules! pt {
    ($x:expr, $y:expr) => {
        $crate::data::Point::new($x as f64, $y as f64)
    };
}

/// Macro for creating a cubic bezier segment that accepts either four points or an array of coordinate tuples
#[macro_export]
macro_rules! cubic {
    // Form 1: Direct points (handles both Point and &mut Point)
    ($p1:expr, $p2:expr, $p3:expr, $p4:expr) => {{
        let p1: &$crate::data::Point = &$p1;
        let p2: &$crate::data::Point = &$p2;
        let p3: &$crate::data::Point = &$p3;
        let p4: &$crate::data::Point = &$p4;
        $crate::data::BezierSegment::cubic(*p1, *p2, *p3, *p4)
    }};
    // Form 2: Array of coordinate tuples
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

/// Macro for creating a quadratic bezier segment that accepts either three points or an array of coordinate tuples
#[macro_export]
macro_rules! quad {
    // Form 1: Direct points (handles both Point and &mut Point)
    ($p1:expr, $p2:expr, $p3:expr) => {{
        let p1: &$crate::data::Point = &$p1;
        let p2: &$crate::data::Point = &$p2;
        let p3: &$crate::data::Point = &$p3;
        $crate::data::BezierSegment::quadratic(*p1, *p2, *p3)
    }};
    // Form 2: Array of coordinate tuples
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

/// Macro for creating a line segment that accepts either two points or an array of coordinate tuples
#[macro_export]
macro_rules! line {
    // Form 1: Direct points (handles both Point and &mut Point)
    ($p1:expr, $p2:expr) => {{
        let p1: &$crate::data::Point = &$p1;
        let p2: &$crate::data::Point = &$p2;
        $crate::data::BezierSegment::line(*p1, *p2)
    }};
    // Form 2: Array of coordinate tuples
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

/// Macro for creating an arc segment that accepts either direct parameters or a tuple form
#[macro_export]
macro_rules! arc {
    // Form 1: Direct parameters
    ($start:expr, $end:expr, $rx:expr, $ry:expr, $angle:expr, $large_arc:expr, $sweep:expr) => {{
        let start: &$crate::data::Point = &$start;
        let end: &$crate::data::Point = &$end;
        $crate::data::BezierSegment::arc(
            *start,
            *end,
            $rx as f64,
            $ry as f64,
            $angle as f64,
            $large_arc,
            $sweep,
        )
    }};
    // Form 2: Tuple form with all parameters
    ([$start:expr, $end:expr, $rx:expr, $ry:expr, $angle:expr, $large_arc:expr, $sweep:expr]) => {{
        let (start_x, start_y) = $start;
        let (end_x, end_y) = $end;
        $crate::data::BezierSegment::arc(
            $crate::pt!(start_x, start_y),
            $crate::pt!(end_x, end_y),
            $rx as f64,
            $ry as f64,
            $angle as f64,
            $large_arc,
            $sweep,
        )
    }};
}
