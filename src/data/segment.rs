//! Bezier segment: quadratic or cubic bezier curve segment

use crate::data::point::Point;
use std::fmt;

/// A bezier segment, either cubic or quadratic
#[derive(Clone, PartialEq)]
pub enum BezierSegment {
    /// Line segment with 2 points
    Line {
        /// Points: start point, end point
        points: [Point; 2],
    },
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
    /// Elliptical arc segment
    Arc {
        /// Start point
        start: Point,
        /// End point
        end: Point,
        /// Radii of the ellipse
        rx: f64,
        ry: f64,
        /// Rotation angle in degrees
        angle: f64,
        /// Whether to use the large arc
        large_arc: bool,
        /// Whether to sweep clockwise
        sweep: bool,
    },
}

impl BezierSegment {
    /// Create a line segment with 2 points
    pub fn line(p1: Point, p2: Point) -> Self {
        Self::Line { points: [p1, p2] }
    }

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

    /// Create an elliptical arc segment
    pub fn arc(
        start: Point,
        end: Point,
        rx: f64,
        ry: f64,
        angle: f64,
        large_arc: bool,
        sweep: bool,
    ) -> Self {
        Self::Arc {
            start,
            end,
            rx,
            ry,
            angle,
            large_arc,
            sweep,
        }
    }

    /// Get all control points for this segment
    pub fn points(&self) -> Vec<Point> {
        match self {
            Self::Line { points } => points.to_vec(),
            Self::Cubic { points } => points.to_vec(),
            Self::Quadratic { points } => points.to_vec(),
            Self::Arc { start, end, .. } => vec![*start, *end],
        }
    }

    /// Get a point on the bezier curve at parameter t (0 <= t <= 1)
    pub fn point_at(&self, t: f64) -> Point {
        match self {
            Self::Line { points } => {
                let p1 = points[0];
                let p2 = points[1];

                let x = p1.x + t * (p2.x - p1.x);
                let y = p1.y + t * (p2.y - p1.y);

                Point::new(x, y)
            }
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
            Self::Arc {
                start,
                end,
                rx,
                ry,
                angle,
                large_arc,
                sweep,
            } => {
                // Convert angle to radians
                let angle_rad = angle.to_radians();

                // Calculate center point
                let dx = end.x - start.x;
                let dy = end.y - start.y;

                // Calculate center point
                let cx = (start.x + end.x) / 2.0;
                let cy = (start.y + end.y) / 2.0;

                // Calculate angle between start and end points
                let theta = dy.atan2(dx);

                // Calculate angle for the arc
                let phi = if *large_arc {
                    if *sweep {
                        theta + std::f64::consts::PI
                    } else {
                        theta - std::f64::consts::PI
                    }
                } else if *sweep {
                    theta
                } else {
                    theta + std::f64::consts::PI
                };

                // Calculate point on arc
                let x = cx + rx * phi.cos() * angle_rad.cos() - ry * phi.sin() * angle_rad.sin();
                let y = cy + rx * phi.cos() * angle_rad.sin() + ry * phi.sin() * angle_rad.cos();

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

    /// Find the nearest point on the curve to a given point using a two-step approach:
    ///     1. Linear sampling to get a good initial guess
    ///     2. Binary search refinement around the initial guess.
    /// this is probably not the best way to do this.
    pub fn nearest_point(&self, point: &Point) -> (Point, f64) {
        // Step 1: Linear sampling to get initial guess
        const LUT_SIZE: usize = 100;
        let mut best_t = 0.0;
        let mut best_point = self.point_at(best_t);
        let mut best_distance = point.distance(&best_point);

        // Create a lookup table of points on the curve
        for i in 0..=LUT_SIZE {
            let t = i as f64 / LUT_SIZE as f64;
            let curve_point = self.point_at(t);
            let distance = point.distance(&curve_point);

            if distance < best_distance {
                best_distance = distance;
                best_t = t;
                best_point = curve_point;
            }
        }

        // Step 2: Binary search refinement around the initial guess
        let mut left = (best_t - 1.0 / LUT_SIZE as f64).max(0.0);
        let mut right = (best_t + 1.0 / LUT_SIZE as f64).min(1.0);
        let tolerance = 0.001;

        while right - left > tolerance {
            let mid1 = left + (right - left) / 3.0;
            let mid2 = right - (right - left) / 3.0;

            let point1 = self.point_at(mid1);
            let point2 = self.point_at(mid2);

            let dist1 = point.distance(&point1);
            let dist2 = point.distance(&point2);

            if dist1 < best_distance {
                best_distance = dist1;
                best_t = mid1;
                best_point = point1;
                right = mid2;
            } else if dist2 < best_distance {
                best_distance = dist2;
                best_t = mid2;
                best_point = point2;
                left = mid1;
            } else {
                left = mid1;
                right = mid2;
            }
        }

        (best_point, best_t)
    }
}

impl fmt::Display for BezierSegment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BezierSegment::Line { points } => {
                write!(f, "Line[{} -> {}]", points[0], points[1])
            }
            BezierSegment::Cubic { points } => {
                write!(
                    f,
                    "Cubic[{} -> {} -> {} -> {}]",
                    points[0], points[1], points[2], points[3]
                )
            }
            BezierSegment::Quadratic { points } => {
                write!(
                    f,
                    "Quadratic[{} -> {} -> {}]",
                    points[0], points[1], points[2]
                )
            }
            BezierSegment::Arc {
                start,
                end,
                rx,
                ry,
                angle,
                large_arc,
                sweep,
            } => {
                write!(
                    f,
                    "Arc[{} -> {}, rx: {:.2}, ry: {:.2}, angle: {:.2}, large_arc: {}, sweep: {}]",
                    start, end, rx, ry, angle, large_arc, sweep
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nearest_point_boundary_cases() {
        // Create a simple cubic bezier segment
        let segment = BezierSegment::cubic(
            Point::new(0.0, 0.0), // Start point
            Point::new(1.0, 1.0), // Control point 1
            Point::new(2.0, 1.0), // Control point 2
            Point::new(3.0, 0.0), // End point
        );

        // Test 1: Point exactly at start (t=0)
        let test_point = Point::new(0.0, 0.0);
        let expected_point = Point::new(0.0, 0.0);
        let (point, t) = segment.nearest_point(&test_point);
        assert_eq!(point, expected_point);
        assert!(t.abs() < 1e-6);

        // Test 2: Point very close to start
        let test_point = Point::new(0.01, 0.01);
        let expected_point = Point::new(0.0, 0.0);
        let (point, t) = segment.nearest_point(&test_point);
        assert!(t < 0.1, "t value was {}", t);
        assert!(point.distance(&expected_point) < 0.1);

        // Test 3: Point exactly at end (t=1)
        let test_point = Point::new(3.0, 0.0);
        let expected_point = Point::new(3.0, 0.0);
        let (point, t) = segment.nearest_point(&test_point);
        assert_eq!(point, expected_point);
        assert!((t - 1.0).abs() < 1e-6);

        // Test 4: Point very close to end
        let test_point = Point::new(2.99, 0.01);
        let expected_point = Point::new(3.0, 0.0);
        let (point, t) = segment.nearest_point(&test_point);
        assert!(t > 0.9, "t value was {}", t);
        assert!(point.distance(&expected_point) < 0.1);
    }

    #[test]
    fn test_nearest_point_middle_cases() {
        // Create a simple cubic bezier segment
        let segment = BezierSegment::cubic(
            Point::new(0.0, 0.0), // Start point
            Point::new(1.0, 1.0), // Control point 1
            Point::new(2.0, 1.0), // Control point 2
            Point::new(3.0, 0.0), // End point
        );

        // Test 1: Point at middle of curve (t=0.5)
        let middle_point = segment.point_at(0.5);
        let test_point = middle_point;
        let expected_point = middle_point;
        let (point, t) = segment.nearest_point(&test_point);
        assert!((t - 0.5).abs() < 0.1, "t value was {}", t);
        assert!(point.distance(&expected_point) < 1e-6);

        // Test 2: Point above middle of curve
        let test_point = Point::new(1.5, 2.0);
        let (_point, t) = segment.nearest_point(&test_point);
        assert!(t > 0.4 && t < 0.6, "t value was {}", t);

        // Test 3: Point below middle of curve
        let test_point = Point::new(1.5, 0.5);
        let (_point, t) = segment.nearest_point(&test_point);
        assert!(t > 0.4 && t < 0.6, "t value was {}", t);
    }

    #[test]
    fn test_nearest_point_complex_curve() {
        // Create a more complex cubic bezier segment
        let segment = BezierSegment::cubic(
            Point::new(0.0, 0.0),  // Start point
            Point::new(0.5, 1.0),  // Control point 1
            Point::new(1.5, -1.0), // Control point 2
            Point::new(2.0, 0.0),  // End point
        );

        // Test 1: Point near a local minimum
        let test_point = Point::new(1.0, 0.2);
        let (_point, t) = segment.nearest_point(&test_point);
        assert!(t > 0.45 && t < 0.55, "t value was {}", t);

        // Test 2: Point near a local maximum
        let test_point = Point::new(0.5, 0.5);
        let (_point, t) = segment.nearest_point(&test_point);
        assert!(t > 0.24 && t < 0.26, "t value was {}", t);

        // Test 3: Point outside the curve's bounding box
        let test_point = Point::new(3.0, 2.0);
        let (_point, t) = segment.nearest_point(&test_point);
        assert!(t > 0.9, "t value was {}", t);
    }
}
