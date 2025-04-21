//! Bezier segment: quadratic or cubic bezier curve segment

use crate::data::point::Point;

/// A bezier segment, either cubic or quadratic
#[derive(Debug, Clone, PartialEq)]
pub enum BezierSegment {
    /// Cubic bezier with 4 control points
    Cubic {
        /// Control points: start point, control1, control2, end point
        points: [Point; 4],
    },
    /// Quadratic bezier with 3 control points
    Quadratic {
        /// Control points: start point, control point, end point
        points: [Point; 3],
    },
}

impl BezierSegment {
    /// Create a cubic segment with 4 control points
    pub fn cubic(p1: Point, p2: Point, p3: Point, p4: Point) -> Self {
        Self::Cubic {
            points: [p1, p2, p3, p4],
        }
    }

    /// Create a quadratic segment with 3 control points
    pub fn quadratic(p1: Point, p2: Point, p3: Point) -> Self {
        Self::Quadratic {
            points: [p1, p2, p3],
        }
    }

    /// Get all control points for this segment
    pub fn points(&self) -> Vec<Point> {
        match self {
            Self::Cubic { points } => points.to_vec(),
            Self::Quadratic { points } => points.to_vec(),
        }
    }

    /// Get a point on the bezier curve at parameter t (0 <= t <= 1)
    pub fn point_at(&self, t: f64) -> Point {
        match self {
            Self::Cubic { points } => {
                let p1 = points[0];
                let p2 = points[1];
                let p3 = points[2];
                let p4 = points[3];

                let t1 = 1.0 - t;

                // B(t) = (1-t)^3 * p1 + 3(1-t)^2 * t * p2 + 3(1-t) * t^2 * p3 + t^3 * p4
                let x = t1.powi(3) * p1.x
                    + 3.0 * t1.powi(2) * t * p2.x
                    + 3.0 * t1 * t.powi(2) * p3.x
                    + t.powi(3) * p4.x;

                let y = t1.powi(3) * p1.y
                    + 3.0 * t1.powi(2) * t * p2.y
                    + 3.0 * t1 * t.powi(2) * p3.y
                    + t.powi(3) * p4.y;

                Point::new(x, y)
            }
            Self::Quadratic { points } => {
                let p1 = points[0];
                let p2 = points[1];
                let p3 = points[2];

                let t1 = 1.0 - t;

                // B(t) = (1-t)^2 * p1 + 2(1-t) * t * p2 + t^2 * p3
                let x = t1.powi(2) * p1.x + 2.0 * t1 * t * p2.x + t.powi(2) * p3.x;
                let y = t1.powi(2) * p1.y + 2.0 * t1 * t * p2.y + t.powi(2) * p3.y;

                Point::new(x, y)
            }
        }
    }

    /// Sample a point at parameter t (0 <= t <= 1)
    /// This is an alias for point_at to maintain consistency with sample_points
    pub fn sample_point_at(&self, t: f64) -> Point {
        self.point_at(t)
    }

    /// Generate a series of points along the bezier curve
    pub fn sample_points(&self, num_points: usize) -> Vec<Point> {
        (0..num_points)
            .map(|i| {
                let t = i as f64 / (num_points - 1) as f64;
                self.sample_point_at(t)
            })
            .collect()
    }

    /// Sample points at specific t values
    pub fn sample_at_t_values(&self, t_values: &[f64]) -> (Vec<Point>, Vec<f64>) {
        let points = t_values.iter().map(|&t| self.point_at(t)).collect();
        (points, t_values.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_bezier_endpoint() {
        let segment = BezierSegment::cubic(
            Point::new(0.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(2.0, 1.0),
            Point::new(3.0, 0.0),
        );

        // At t=0, should be at first control point
        let start = segment.point_at(0.0);
        assert_eq!(start, Point::new(0.0, 0.0));

        // At t=1, should be at last control point
        let end = segment.point_at(1.0);
        assert_eq!(end, Point::new(3.0, 0.0));
    }

    #[test]
    fn test_segment_sampling() {
        // Create a quadratic bezier segment
        let segment = BezierSegment::quadratic(
            Point::new(0.0, 0.0),    // Start point
            Point::new(50.0, 100.0), // Control point
            Point::new(100.0, 0.0),  // End point
        );

        let mid_point = segment.sample_point_at(0.5);
        assert_eq!(mid_point, Point::new(50.0, 50.0));

        // Sample 5 points from the segment
        let samples = segment.sample_points(5);

        // Check we got the right number of points
        assert_eq!(samples.len(), 5);

        // start point
        assert_eq!(samples[0], Point::new(0.0, 0.0));

        // Check sample 1 (t=0.25)
        // For a quadratic curve: B(0.25) = 0.5625*P0 + 0.375*P1 + 0.0625*P2
        // = 0.5625*(0,0) + 0.375*(50,100) + 0.0625*(100,0)
        // = (0,0) + (18.75,37.5) + (6.25,0) = (25,37.5)
        assert_eq!(samples[1], Point::new(25.0, 37.5));

        // The middle point (t=0.5) should be at (x=50, y=50) for this specific curve
        // For a quadratic curve: B(0.5) = 0.25*P0 + 0.5*P1 + 0.25*P2
        // = 0.25*(0,0) + 0.5*(50,100) + 0.25*(100,0)
        // = (0,0) + (25,50) + (25,0) = (50,50)
        assert_eq!(samples[2], Point::new(50.0, 50.0));

        // Check sample 3 (t=0.75)
        // For a quadratic curve: B(0.75) = 0.0625*P0 + 0.375*P1 + 0.5625*P2
        // = 0.0625*(0,0) + 0.375*(50,100) + 0.5625*(100,0)
        // = (0,0) + (18.75,37.5) + (56.25,0) = (75,37.5)
        assert_eq!(samples[3], Point::new(75.0, 37.5));

        // end point
        assert_eq!(samples[4], Point::new(100.0, 0.0));
    }
}
