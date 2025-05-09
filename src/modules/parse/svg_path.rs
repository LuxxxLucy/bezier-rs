use crate::curve;
use crate::data::{BezierCurve, BezierSegment, Point};
use std::error::Error;

/// Parse SVG path data into a BezierCurve
pub trait FromSvgPath: Sized {
    /// Parse from SVG path data string
    fn from_svg_path(data: &str) -> Result<Self, Box<dyn Error>>;
}

impl FromSvgPath for BezierCurve {
    fn from_svg_path(data: &str) -> Result<Self, Box<dyn Error>> {
        let mut segments = vec![];
        let mut current_point = Point::new(0.0, 0.0);
        let mut start_point = current_point;
        let mut is_closed = false;
        let mut current_command = ' ';

        let mut numbers = vec![];
        let mut current_number = String::new();

        // Process each character
        for c in data.chars() {
            match c {
                'M' | 'm' | 'C' | 'c' | 'Q' | 'q' | 'L' | 'l' | 'Z' | 'z' => {
                    // Process any pending number
                    if !current_number.is_empty() {
                        numbers.push(current_number.parse::<f64>()?);
                        current_number.clear();
                    }

                    // Process previous command's numbers
                    if !numbers.is_empty() {
                        process_command(
                            current_command,
                            &numbers,
                            &mut current_point,
                            &mut start_point,
                            &mut segments,
                        )?;
                        numbers.clear();
                    }

                    current_command = c;
                }
                '0'..='9' | '.' | '-' | '+' | 'e' | 'E' => {
                    current_number.push(c);
                }
                ',' | ' ' => {
                    if !current_number.is_empty() {
                        numbers.push(current_number.parse::<f64>()?);
                        current_number.clear();
                    }
                }
                _ => {}
            }
        }

        // Process any remaining number
        if !current_number.is_empty() {
            numbers.push(current_number.parse::<f64>()?);
        }

        // Process final command
        if !numbers.is_empty() {
            process_command(
                current_command,
                &numbers,
                &mut current_point,
                &mut start_point,
                &mut segments,
            )?;
        }

        // Handle Z command if it was the last one
        if current_command == 'Z' || current_command == 'z' {
            if current_point != start_point {
                segments.push(BezierSegment::quadratic(
                    current_point,
                    Point::new(
                        (current_point.x + start_point.x) / 2.0,
                        (current_point.y + start_point.y) / 2.0,
                    ),
                    start_point,
                ));
            }
            is_closed = true;
        }

        if segments.is_empty() {
            return Err("Cannot create curve from empty path".into());
        }

        if is_closed {
            BezierCurve::new_closed(segments)
                .ok_or_else(|| "Failed to close path: end point does not match start point".into())
        } else {
            Ok(curve!(segments))
        }
    }
}

fn process_command(
    command: char,
    numbers: &[f64],
    current_point: &mut Point,
    start_point: &mut Point,
    segments: &mut Vec<BezierSegment>,
) -> Result<(), Box<dyn Error>> {
    match command {
        'M' | 'm' => {
            if numbers.len() >= 2 {
                *current_point = Point::new(numbers[0], numbers[1]);
                *start_point = *current_point;
            }
        }
        'C' | 'c' => {
            if numbers.len() >= 6 {
                let p1 = Point::new(numbers[0], numbers[1]);
                let p2 = Point::new(numbers[2], numbers[3]);
                let end = Point::new(numbers[4], numbers[5]);
                segments.push(BezierSegment::cubic(*current_point, p1, p2, end));
                *current_point = end;
            }
        }
        'Q' | 'q' => {
            if numbers.len() >= 4 {
                let control = Point::new(numbers[0], numbers[1]);
                let end = Point::new(numbers[2], numbers[3]);
                segments.push(BezierSegment::quadratic(*current_point, control, end));
                *current_point = end;
            }
        }
        'L' | 'l' => {
            if numbers.len() >= 2 {
                let end = Point::new(numbers[0], numbers[1]);
                let control = Point::new(
                    (current_point.x + end.x) / 2.0,
                    (current_point.y + end.y) / 2.0,
                );
                segments.push(BezierSegment::quadratic(*current_point, control, end));
                *current_point = end;
            }
        }
        _ => {}
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_empty_path() {
        assert!(BezierCurve::from_svg_path("").is_err());
    }

    #[test]
    fn test_parse_single_path() {
        let curve = BezierCurve::from_svg_path("M10,10 C20,20 40,20 50,10").unwrap();

        match &curve.segments[0] {
            BezierSegment::Cubic { points } => {
                assert_eq!(points[0], Point::new(10.0, 10.0));
                assert_eq!(points[1], Point::new(20.0, 20.0));
                assert_eq!(points[2], Point::new(40.0, 20.0));
                assert_eq!(points[3], Point::new(50.0, 10.0));
            }
            _ => panic!("Expected cubic segment"),
        }
    }

    #[test]
    fn test_parse_multiple_subpaths() {
        let curve = BezierCurve::from_svg_path("M10,10 C20,20 40,20 50,10 Q60,0 70,10").unwrap();

        match &curve.segments[0] {
            BezierSegment::Cubic { points } => {
                assert_eq!(points[0], Point::new(10.0, 10.0));
                assert_eq!(points[1], Point::new(20.0, 20.0));
                assert_eq!(points[2], Point::new(40.0, 20.0));
                assert_eq!(points[3], Point::new(50.0, 10.0));
            }
            _ => panic!("Expected cubic segment"),
        }

        match &curve.segments[1] {
            BezierSegment::Quadratic { points } => {
                assert_eq!(points[0], Point::new(50.0, 10.0));
                assert_eq!(points[1], Point::new(60.0, 0.0));
                assert_eq!(points[2], Point::new(70.0, 10.0));
            }
            _ => panic!("Expected quadratic segment"),
        }
    }

    #[test]
    fn test_parse_closed_path() {
        let curve = BezierCurve::from_svg_path("M10,10 C20,20 40,20 50,10 Q60,0 70,10 Z").unwrap();
        assert!(curve.is_closed());
    }

    #[test]
    fn test_parse_line_segments() {
        let curve = BezierCurve::from_svg_path("M10,10 L20,20 L30,30").unwrap();
        assert_eq!(curve.segments.len(), 2);
        match &curve.segments[0] {
            BezierSegment::Quadratic { points } => {
                assert_eq!(points[0], Point::new(10.0, 10.0));
                assert_eq!(points[1], Point::new(15.0, 15.0));
                assert_eq!(points[2], Point::new(20.0, 20.0));
            }
            _ => panic!("Expected quadratic segment"),
        }

        match &curve.segments[1] {
            BezierSegment::Quadratic { points } => {
                assert_eq!(points[0], Point::new(20.0, 20.0));
                assert_eq!(points[1], Point::new(25.0, 25.0));
                assert_eq!(points[2], Point::new(30.0, 30.0));
            }
            _ => panic!("Expected quadratic segment"),
        }
    }

    #[test]
    fn test_round_trip_parse_and_then_export() {
        use crate::modules::export::svg_path::ToSvgPath;

        let test_cases = [
            "M10,20 C20,30,30,40,40,50",
            "M10,20 Q40,50,70,80",
            "M10,20 C20,30,30,40,40,50 Q50,60,60,70",
            "M10,10 C20,20,40,20,50,10 Q60,0,70,10 Q40,10,10,10Z",
        ];

        for (i, original_path) in test_cases.iter().enumerate() {
            // Parse SVG path
            let curve = BezierCurve::from_svg_path(original_path)
                .unwrap_or_else(|e| panic!("Failed to parse path data for test case {}: {}", i, e));

            // Export back to SVG path
            let exported_path = curve.to_svg_path();

            // Compare the paths
            assert_eq!(
                *original_path, exported_path,
                "Path mismatch in test case {}",
                i
            );
        }
    }
}
