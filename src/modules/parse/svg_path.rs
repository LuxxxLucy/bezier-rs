use crate::curve;
use crate::data::{BezierCurve, BezierSegment, Point};
use crate::{cubic, line, pt, quad};
use std::error::Error; // Import macros
use std::fmt;

/// A generic parsing entity that tracks position and length in a string
struct ParsingEntity {
    start: usize,
    len: usize,
}

impl ParsingEntity {
    fn new() -> Self {
        ParsingEntity { start: 0, len: 0 }
    }

    fn reset(&mut self) {
        self.start = 0;
        self.len = 0;
    }

    fn parse<T: std::str::FromStr>(&self, data: &str) -> Option<T> {
        if self.len == 0 {
            return None;
        }
        data[self.start..self.start + self.len].parse::<T>().ok()
    }

    fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// Parse SVG path data into a BezierCurve
pub trait FromSvgPath: Sized {
    /// Parse from SVG path data string
    fn from_svg_path(data: &str) -> Result<Self, Box<dyn Error>>;
}

impl FromSvgPath for BezierCurve {
    fn from_svg_path(data: &str) -> Result<Self, Box<dyn Error>> {
        let (result, _) = Self::parse_one_svg_path(data)?;
        Ok(result)
    }
}

/// Custom error for SVG path parsing
#[derive(Debug)]
pub enum SvgPathParseError {
    MultiplePaths { bytes_consumed: usize },
    Other(Box<dyn Error>),
}

impl fmt::Display for SvgPathParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SvgPathParseError::MultiplePaths { bytes_consumed } => {
                write!(
                    f,
                    "it seems it has multiple path, bytes_consumed={}",
                    bytes_consumed
                )
            }
            SvgPathParseError::Other(e) => write!(f, "{}", e),
        }
    }
}

impl Error for SvgPathParseError {}

impl BezierCurve {
    /// Parse a single SVG path and return both the result and bytes consumed
    pub fn parse_one_svg_path(data: &str) -> Result<(Self, usize), Box<dyn Error>> {
        let mut segments = vec![];
        let mut current_point = pt!(0.0, 0.0);
        let mut start_point = current_point;
        let mut current_command = ' ';
        let mut last_was_move = false; // Track if the last command was a move
        let mut bytes_consumed = 0;

        let mut numbers = vec![];
        let mut current_number = ParsingEntity::new();

        let chars = data.char_indices();
        for (i, c) in chars {
            bytes_consumed = i + c.len_utf8();
            match c {
                'M' | 'm' | 'C' | 'c' | 'Q' | 'q' | 'L' | 'l' | 'H' | 'h' | 'V' | 'v' | 'S'
                | 's' | 'Z' | 'z' | 'A' | 'a' | 'T' | 't' => {
                    // Process any pending number
                    if !current_number.is_empty() {
                        if let Some(num) = current_number.parse::<f64>(data) {
                            numbers.push(num);
                        }
                        current_number.reset();
                    }
                    // Process previous command's numbers
                    if !numbers.is_empty() {
                        process_command(
                            current_command,
                            &numbers,
                            &mut current_point,
                            &mut start_point,
                            &mut segments,
                            &mut last_was_move,
                        )?;
                        numbers.clear();
                    }

                    // Special handling for Z/z
                    if c == 'Z' || c == 'z' {
                        // Handle Z command
                        if current_point != start_point {
                            segments.push(line!(current_point, start_point));
                        }
                        // Return immediately after Z/z
                        // bytes_consumed should include this Z/z
                        bytes_consumed = i + c.len_utf8();
                        // If there are no segments, return an empty curve
                        if segments.is_empty() {
                            if bytes_consumed != data.len() {
                                return Err(Box::new(SvgPathParseError::MultiplePaths {
                                    bytes_consumed,
                                }));
                            }
                            return Ok((curve!([]), bytes_consumed));
                        }
                        let curve = BezierCurve::new_closed(segments).ok_or_else(|| {
                            Box::new(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                "Failed to create closed curve",
                            ))
                        })?;
                        if bytes_consumed != data.len() {
                            return Err(Box::new(SvgPathParseError::MultiplePaths {
                                bytes_consumed,
                            }));
                        }
                        return Ok((curve, bytes_consumed));
                    }

                    current_command = c;
                }
                '0'..='9' | '.' | '+' | 'e' | 'E' | '-' => {
                    if c == '-' {
                        if !current_number.is_empty()
                            && !data
                                [current_number.start..current_number.start + current_number.len]
                                .ends_with(['e', 'E'])
                        {
                            if let Some(num) = current_number.parse::<f64>(data) {
                                numbers.push(num);
                            }
                            current_number.reset();
                            current_number.start = i;
                            current_number.len = 1;
                        } else {
                            if current_number.is_empty() {
                                current_number.start = i;
                            }
                            current_number.len += 1;
                        }
                    } else if c == '+' {
                        // Handle '+' as separator if it's not part of scientific notation
                        if !current_number.is_empty()
                            && !data
                                [current_number.start..current_number.start + current_number.len]
                                .ends_with(['e', 'E'])
                        {
                            if let Some(num) = current_number.parse::<f64>(data) {
                                numbers.push(num);
                            }
                            current_number.reset();
                        } else {
                            if current_number.is_empty() {
                                current_number.start = i;
                            }
                            current_number.len += 1;
                        }
                    } else if c == '.' {
                        if !current_number.is_empty()
                            && data[current_number.start..current_number.start + current_number.len]
                                .contains('.')
                        {
                            if let Some(num) = current_number.parse::<f64>(data) {
                                numbers.push(num);
                            }
                            current_number.reset();
                            current_number.start = i;
                            current_number.len = 1;
                        } else {
                            if current_number.is_empty() {
                                current_number.start = i;
                            }
                            current_number.len += 1;
                        }
                    } else {
                        if current_number.is_empty() {
                            current_number.start = i;
                        }
                        current_number.len += 1;
                    }
                }
                ',' | ' ' | '\n' | '\r' | '\t' => {
                    if !current_number.is_empty() {
                        if let Some(num) = current_number.parse::<f64>(data) {
                            numbers.push(num);
                        }
                        current_number.reset();
                    }
                }
                _ => {}
            }
        }

        // Process any remaining number
        if !current_number.is_empty() {
            if let Some(num) = current_number.parse::<f64>(data) {
                numbers.push(num);
            }
        }

        // Process final command
        if !numbers.is_empty() {
            process_command(
                current_command,
                &numbers,
                &mut current_point,
                &mut start_point,
                &mut segments,
                &mut last_was_move,
            )?;
        }

        // If there are no segments, return an empty curve
        if segments.is_empty() {
            if bytes_consumed != data.len() {
                return Err(Box::new(SvgPathParseError::MultiplePaths {
                    bytes_consumed,
                }));
            }
            Ok((curve!([]), bytes_consumed))
        } else {
            if bytes_consumed != data.len() {
                return Err(Box::new(SvgPathParseError::MultiplePaths {
                    bytes_consumed,
                }));
            }
            Ok((curve!(segments), bytes_consumed))
        }
    }

    // Helper function to parse one SVG path, ignoring the 'multiple paths' error
    fn parse_one_svg_path_ignore_multiple(data: &str) -> Result<(Self, usize), Box<dyn Error>> {
        match Self::parse_one_svg_path(data) {
            Ok(result) => Ok(result),
            Err(e) => {
                if let Some(mp) = e.downcast_ref::<SvgPathParseError>() {
                    match mp {
                        SvgPathParseError::MultiplePaths { bytes_consumed } => {
                            let (curve, _) = Self::parse_one_svg_path(&data[..*bytes_consumed])?;
                            Ok((curve, *bytes_consumed))
                        }
                        SvgPathParseError::Other(_) => Err(e),
                    }
                } else {
                    Err(e)
                }
            }
        }
    }

    /// Parse potentially multiple SVG paths separated by Z commands
    ///
    /// This function is designed to handle SVG path data that may contain multiple paths
    /// separated by Z/z (closepath) commands. It will parse each path individually and
    /// return a vector of curves.
    ///
    /// # Examples
    ///
    /// ```
    /// use bezier_rs::BezierCurve;
    ///
    /// // Parse a single path with cubic Bézier curve
    /// let single_path = "M 10,10 C 20,20 40,20 50,10 Z";
    /// let curves = BezierCurve::parse_maybe_multiple(single_path).unwrap();
    /// assert_eq!(curves.len(), 1);
    /// assert_eq!(curves[0].segments.len(), 2);
    ///
    /// // Parse multiple paths with cubic Bézier curves
    /// let multiple_paths = "M 10,10 C 20,20 40,20 50,10 Z M 30,30 C 40,40 50,50 60,60 Z";
    /// let curves = BezierCurve::parse_maybe_multiple(multiple_paths).unwrap();
    /// assert_eq!(curves.len(), 2);
    /// assert_eq!(curves[0].segments.len(), 2);
    /// assert_eq!(curves[1].segments.len(), 2);
    /// ```
    ///
    /// # Notes
    ///
    /// - Each path must end with a Z/z command to be properly separated
    /// - If a path doesn't end with Z/z, it will be treated as a single path
    /// - Use this function when you expect multiple paths in the input
    /// - For single paths, consider using `parse_one_svg_path` instead
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The path data is invalid
    /// - A path command is malformed
    /// - The path cannot be parsed into valid curves
    pub fn parse_maybe_multiple(data: &str) -> Result<Vec<Self>, Box<dyn Error>> {
        let mut curves = Vec::new();
        let mut remaining = data;

        while !remaining.is_empty() {
            let (curve, consumed) = Self::parse_one_svg_path_ignore_multiple(remaining)?;
            curves.push(curve);

            // Skip any whitespace after the consumed bytes
            let mut next_start = consumed;
            while next_start < remaining.len()
                && remaining[next_start..]
                    .chars()
                    .next()
                    .map_or(false, |c| c.is_whitespace())
            {
                next_start += remaining[next_start..].chars().next().unwrap().len_utf8();
            }

            remaining = &remaining[next_start..];
        }

        Ok(curves)
    }
}

fn process_command(
    command: char,
    numbers: &[f64],
    current_point: &mut Point,
    start_point: &mut Point,
    segments: &mut Vec<BezierSegment>,
    last_was_move: &mut bool,
) -> Result<(), Box<dyn Error>> {
    let mut curr = 0;
    while curr < numbers.len() {
        match command {
            'M' | 'm' => {
                if curr + 1 >= numbers.len() {
                    return Err(format!(
                        "M/m command requires at least 2 numbers, but only {} remaining at position {}",
                        numbers.len() - curr,
                        curr
                    )
                    .into());
                }
                let is_relative = command == 'm';
                // If last command was also a move, add a line segment
                if *last_was_move {
                    let end_point = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };
                    segments.push(line!(current_point, end_point));
                }
                // First point is a move
                *current_point = if is_relative {
                    pt!(
                        current_point.x + numbers[curr],
                        current_point.y + numbers[curr + 1]
                    )
                } else {
                    pt!(numbers[curr], numbers[curr + 1])
                };
                *start_point = *current_point;
                *last_was_move = true;
                curr += 2;

                // Process remaining points as implicit L commands
                while curr + 1 < numbers.len() {
                    let end_point = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };
                    segments.push(line!(current_point, end_point));
                    *current_point = end_point;
                    curr += 2;
                }
            }
            'C' | 'c' => {
                let is_relative = command == 'c';
                // Process all coordinate triplets as implicit cubic Bézier commands
                while curr + 5 < numbers.len() {
                    let p1 = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };
                    let p2 = if is_relative {
                        pt!(
                            current_point.x + numbers[curr + 2],
                            current_point.y + numbers[curr + 3]
                        )
                    } else {
                        pt!(numbers[curr + 2], numbers[curr + 3])
                    };
                    let end = if is_relative {
                        pt!(
                            current_point.x + numbers[curr + 4],
                            current_point.y + numbers[curr + 5]
                        )
                    } else {
                        pt!(numbers[curr + 4], numbers[curr + 5])
                    };
                    segments.push(cubic!(*current_point, p1, p2, end));
                    *current_point = end;
                    curr += 6;
                }
            }
            'Q' | 'q' => {
                let is_relative = command == 'q';
                // Process all coordinate pairs as implicit quadratic Bézier commands
                while curr + 3 < numbers.len() {
                    let p1 = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };
                    let end = if is_relative {
                        pt!(
                            current_point.x + numbers[curr + 2],
                            current_point.y + numbers[curr + 3]
                        )
                    } else {
                        pt!(numbers[curr + 2], numbers[curr + 3])
                    };
                    segments.push(quad!(*current_point, p1, end));
                    *current_point = end;
                    curr += 4;
                }
            }
            'L' | 'l' => {
                if curr + 1 >= numbers.len() {
                    return Err(format!(
                        "L/l command requires at least 2 numbers, but only {} remaining at position {}",
                        numbers.len() - curr,
                        curr
                    )
                    .into());
                }
                let is_relative = command == 'l';
                // Process all coordinate pairs as implicit line commands
                while curr + 1 < numbers.len() {
                    let end_point = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };
                    segments.push(line!(current_point, end_point));
                    *current_point = end_point;
                    curr += 2;
                }
            }
            'H' | 'h' => {
                if curr >= numbers.len() {
                    return Err(format!(
                        "H/h command requires at least 1 number, but none remaining at position {}",
                        curr
                    )
                    .into());
                }
                let is_relative = command == 'h';
                while curr < numbers.len() {
                    let end_point = if is_relative {
                        pt!(current_point.x + numbers[curr], current_point.y)
                    } else {
                        pt!(numbers[curr], current_point.y)
                    };
                    segments.push(line!(current_point, end_point));
                    *current_point = end_point;
                    curr += 1;
                }
            }
            'V' | 'v' => {
                if curr >= numbers.len() {
                    return Err(format!(
                        "V/v command requires at least 1 number, but none remaining at position {}",
                        curr
                    )
                    .into());
                }
                let is_relative = command == 'v';
                while curr < numbers.len() {
                    let end_point = if is_relative {
                        pt!(current_point.x, current_point.y + numbers[curr])
                    } else {
                        pt!(current_point.x, numbers[curr])
                    };
                    segments.push(line!(current_point, end_point));
                    *current_point = end_point;
                    curr += 1;
                }
            }
            'S' | 's' => {
                let is_relative = command == 's';
                while curr + 3 < numbers.len() {
                    let p2 = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };
                    let end = if is_relative {
                        pt!(
                            current_point.x + numbers[curr + 2],
                            current_point.y + numbers[curr + 3]
                        )
                    } else {
                        pt!(numbers[curr + 2], numbers[curr + 3])
                    };

                    // Calculate first control point based on previous segment
                    let p1 = if let Some(prev) = segments.last() {
                        match prev {
                            BezierSegment::Cubic { points } => {
                                // Reflect the second control point of the previous curve
                                let prev_ctrl = points[2];
                                pt!(
                                    2.0 * current_point.x - prev_ctrl.x,
                                    2.0 * current_point.y - prev_ctrl.y
                                )
                            }
                            _ => *current_point, // If previous segment is not cubic, use current point
                        }
                    } else {
                        *current_point // If no previous segment, use current point
                    };

                    segments.push(cubic!(*current_point, p1, p2, end));
                    *current_point = end;
                    curr += 4;
                }
            }
            'T' | 't' => {
                let is_relative = command == 't';
                while curr + 1 < numbers.len() {
                    let end = if is_relative {
                        pt!(
                            current_point.x + numbers[curr],
                            current_point.y + numbers[curr + 1]
                        )
                    } else {
                        pt!(numbers[curr], numbers[curr + 1])
                    };

                    // Calculate control point based on previous segment
                    let p1 = if let Some(prev) = segments.last() {
                        match prev {
                            BezierSegment::Quadratic { points } => {
                                // Reflect the control point of the previous curve
                                let prev_ctrl = points[1];
                                pt!(
                                    2.0 * current_point.x - prev_ctrl.x,
                                    2.0 * current_point.y - prev_ctrl.y
                                )
                            }
                            _ => *current_point, // If previous segment is not quadratic, use current point
                        }
                    } else {
                        *current_point // If no previous segment, use current point
                    };

                    segments.push(quad!(*current_point, p1, end));
                    *current_point = end;
                    curr += 2;
                }
            }
            'A' | 'a' => {
                let is_relative = command == 'a';
                while curr + 6 < numbers.len() {
                    let rx = numbers[curr];
                    let ry = numbers[curr + 1];
                    let angle = numbers[curr + 2];
                    let large_arc = numbers[curr + 3] != 0.0;
                    let sweep = numbers[curr + 4] != 0.0;
                    let end = if is_relative {
                        pt!(
                            current_point.x + numbers[curr + 5],
                            current_point.y + numbers[curr + 6]
                        )
                    } else {
                        pt!(numbers[curr + 5], numbers[curr + 6])
                    };

                    segments.push(BezierSegment::arc(
                        *current_point,
                        end,
                        rx,
                        ry,
                        angle,
                        large_arc,
                        sweep,
                    ));
                    *current_point = end;
                    curr += 7;
                }
            }
            'Z' | 'z' => {
                // Close the path by adding a line to the start point
                if *current_point != *start_point {
                    segments.push(line!(current_point, *start_point));
                    *current_point = *start_point;
                }
                curr += 1;
            }
            _ => {
                return Err(format!("Unknown command '{}' at position {}", command, curr).into());
            }
        }
    }

    // Ensure all numbers were processed
    if curr != numbers.len() {
        return Err(format!(
            "Command '{}' has {} unprocessed numbers at position {}",
            command,
            numbers.len() - curr,
            curr
        )
        .into());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{arc, cubic, line, quad};

    // Helper function to run test cases
    fn run_test_cases(test_cases: &[(&str, &str, Vec<BezierSegment>)]) {
        for (test_name, path, expected_segments) in test_cases.iter() {
            let (curve, bytes_consumed) = BezierCurve::parse_one_svg_path(path)
                .unwrap_or_else(|e| panic!("Failed to parse path in test '{}': {}", test_name, e));

            let expected_curve = curve!(expected_segments.to_vec());
            assert!(
                curve == expected_curve,
                "Wrong segment count in test '{}', path string is {} get is {}, expected is {}",
                test_name,
                path,
                curve,
                expected_curve
            );

            assert!(
                bytes_consumed == path.len(),
                "Wrong bytes consumed in test '{}', expected {} but got {}",
                test_name,
                path.len(),
                bytes_consumed
            );
        }
    }

    #[test]
    fn test_parse_move_close_path() {
        let test_cases = [(
            "Empty path with just M and Z",
            "M 138.42576,121.18355 Z",
            vec![],
        )];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_single_path() {
        let test_cases = [(
            "Single cubic bezier path",
            "M 10,10 C 20,20 40,20 50,10",
            vec![cubic!([
                (10.0, 10.0),
                (20.0, 20.0),
                (40.0, 20.0),
                (50.0, 10.0)
            ])],
        )];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_multiple_subpaths() {
        let test_cases = [(
            "Multiple subpaths with cubic and quadratic",
            "M 10,10 C 20,20 40,20 50,10 Q 60,0 70,10",
            vec![
                cubic!([(10.0, 10.0), (20.0, 20.0), (40.0, 20.0), (50.0, 10.0)]),
                quad!([(50.0, 10.0), (60.0, 0.0), (70.0, 10.0)]),
            ],
        )];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_line_segments() {
        let test_cases = [
            (
                "Multiple line segments",
                "M 10,10 L 20,20 L 30,30",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                ],
            ),
            (
                "Multiple coordinates in L command",
                "M 10,10 L 20,20 30,30 40,40",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Multiple coordinates in l command",
                "M 10,10 l 10,10 10,10 10,10",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Mixed absolute and relative line commands",
                "M 10,10 L 20,20 l 10,10 L 40,40 l 10,10",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                    line!([(40.0, 40.0), (50.0, 50.0)]),
                ],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_almost_all_commands() {
        let test_cases = [
            (
                "Absolute line",
                "M 10,10 L 20,20",
                vec![line!([(10.0, 10.0), (20.0, 20.0)])],
            ),
            (
                "Relative line",
                "M 10,10 l 10,10",
                vec![line!([(10.0, 10.0), (20.0, 20.0)])],
            ),
            (
                "Horizontal line",
                "M 10,10 H 20",
                vec![line!([(10.0, 10.0), (20.0, 10.0)])],
            ),
            (
                "Relative horizontal line",
                "M 10,10 h 10",
                vec![line!([(10.0, 10.0), (20.0, 10.0)])],
            ),
            (
                "Vertical line",
                "M 10,10 V 20",
                vec![line!([(10.0, 10.0), (10.0, 20.0)])],
            ),
            (
                "Relative vertical line",
                "M 10,10 v 10",
                vec![line!([(10.0, 10.0), (10.0, 20.0)])],
            ),
            (
                "Multiple absolute lines",
                "M 10,10 L 20,20 30,30 40,40",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Multiple relative lines",
                "M 10,10 l 10,10 10,10 10,10",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Multiple points in M command",
                "M 10,10 20,20 30,30 40,40",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Multiple points in m command",
                "M 10,10 m 10,10 10,10 10,10",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Consecutive move commands",
                "M 10,10 M 20,20 M 30,30",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                ],
            ),
            (
                "Consecutive relative move commands",
                "M 10,10 m 10,10 m 10,10",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                ],
            ),
            (
                "Mixed absolute and relative move commands",
                "M 10,10 m 10,10 M 30,30 m 10,10",
                vec![
                    line!([(10.0, 10.0), (20.0, 20.0)]),
                    line!([(20.0, 20.0), (30.0, 30.0)]),
                    line!([(30.0, 30.0), (40.0, 40.0)]),
                ],
            ),
            (
                "Simple smooth cubic with previous cubic",
                "M 10,90 C 30,90 25,10 50,10 S 70,90 90,90",
                vec![
                    cubic!([(10.0, 90.0), (30.0, 90.0), (25.0, 10.0), (50.0, 10.0)]),
                    cubic!([(50.0, 10.0), (75.0, 10.0), (70.0, 90.0), (90.0, 90.0)]),
                ],
            ),
            (
                "Smooth cubic without previous cubic",
                "M10,10 S20,20 30,30",
                vec![cubic!([
                    (10.0, 10.0),
                    (10.0, 10.0),
                    (20.0, 20.0),
                    (30.0, 30.0)
                ])],
            ),
            (
                "Smooth quadratic with previous quadratic",
                "M10,10 Q20,20 30,30 T50,50",
                vec![
                    quad!([(10.0, 10.0), (20.0, 20.0), (30.0, 30.0)]),
                    quad!([(30.0, 30.0), (40.0, 40.0), (50.0, 50.0)]),
                ],
            ),
            (
                "Smooth quadratic without previous quadratic",
                "M10,10 T30,30",
                vec![quad!([(10.0, 10.0), (10.0, 10.0), (30.0, 30.0)])],
            ),
            (
                "Relative smooth quadratic",
                "M10,10 t20,20",
                vec![quad!([(10.0, 10.0), (10.0, 10.0), (30.0, 30.0)])],
            ),
            (
                "Simple arc",
                "M 10,10 A 5,5 0,0,1 20,20",
                vec![arc!([
                    (10.0, 10.0),
                    (20.0, 20.0),
                    5.0,
                    5.0,
                    0.0,
                    false,
                    true
                ])],
            ),
            (
                "Relative arc",
                "M10,10 a5,5 0 0 1 10,10",
                vec![arc!([
                    (10.0, 10.0),
                    (20.0, 20.0),
                    5.0,
                    5.0,
                    0.0,
                    false,
                    true
                ])],
            ),
            (
                "Arc with rotation and flags",
                "M10,10 A5,5 45 1 0 20,20",
                vec![arc!([
                    (10.0, 10.0),
                    (20.0, 20.0),
                    5.0,
                    5.0,
                    45.0,
                    true,
                    false
                ])],
            ),
            (
                "Implicit arc commands",
                "M10,10 A5,5 0 0 1 20,20 5,5 0 0 1 30,30",
                vec![
                    arc!([(10.0, 10.0), (20.0, 20.0), 5.0, 5.0, 0.0, false, true]),
                    arc!([(20.0, 20.0), (30.0, 30.0), 5.0, 5.0, 0.0, false, true]),
                ],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_plus_as_separator() {
        let test_cases = [
            (
                "Simple path with plus separator",
                "M10+20L30+40",
                vec![line!([(10.0, 20.0), (30.0, 40.0)])],
            ),
            (
                "Complex path with plus separator",
                "M334.763+412.827C334.899+412.844+332.531+416.181+332.303+416.574",
                vec![cubic!([
                    (334.763, 412.827),
                    (334.899, 412.844),
                    (332.531, 416.181),
                    (332.303, 416.574)
                ])],
            ),
            (
                "Mixed separators",
                "M10+20,30+40",
                vec![line!([(10.0, 20.0), (30.0, 40.0)])],
            ),
            (
                "Scientific notation with plus separator",
                "M10+20,1e+4+2e-4",
                vec![line!([(10.0, 20.0), (10000.0, 0.0002)])],
            ),
            (
                "Multiple segments with plus separator",
                "M10+20L30+40L50+60",
                vec![
                    line!([(10.0, 20.0), (30.0, 40.0)]),
                    line!([(30.0, 40.0), (50.0, 60.0)]),
                ],
            ),
            (
                "Cubic bezier with plus separator",
                "M10+20C30+40+50+60+70+80",
                vec![cubic!([
                    (10.0, 20.0),
                    (30.0, 40.0),
                    (50.0, 60.0),
                    (70.0, 80.0)
                ])],
            ),
            (
                "Quadratic bezier with plus separator",
                "M10+20Q30+40+50+60",
                vec![quad!([(10.0, 20.0), (30.0, 40.0), (50.0, 60.0)])],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_scientific_notation() {
        let test_cases = [
            (
                "Simple scientific notation",
                "M 1e-4,2e+4 C 3e-6,4e+6 5e-8,6e+8 7e-10,8e+10",
                vec![cubic!([
                    (0.0001, 20000.0),
                    (0.000003, 4000000.0),
                    (0.00000005, 600000000.0),
                    (0.0000000007, 80000000000.0)
                ])],
            ),
            (
                "Decimal scientific notation",
                "M 1.23e-4,2.34e+4 C 3.45e-6,4.56e+6 5.67e-8,6.78e+8 7.89e-10,8.90e+10",
                vec![cubic!([
                    (0.000123, 23400.0),
                    (0.00000345, 4560000.0),
                    (0.0000000567, 678000000.0),
                    (0.000000000789, 89000000000.0)
                ])],
            ),
            (
                "Uppercase scientific notation",
                "M 1E-4,2E+4 C 3E-6,4E+6 5E-8,6E+8 7E-10,8E+10",
                vec![cubic!([
                    (0.0001, 20000.0),
                    (0.000003, 4000000.0),
                    (0.00000005, 600000000.0),
                    (0.0000000007, 80000000000.0)
                ])],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_numbers_starting_with_dot() {
        let test_cases = [
            (
                "Numbers starting with decimal point",
                "M.133 51.647l241.4.534-10.67 46.39-167.95-30.92L.136 64.31z",
                vec![
                    line!([(0.133, 51.647), (241.533, 52.181)]),
                    line!([(241.533, 52.181), (230.863, 98.571)]),
                    line!([(230.863, 98.571), (62.913, 67.651)]),
                    line!([(62.913, 67.651), (0.136, 64.31)]),
                    line!([(0.136, 64.31), (0.133, 51.647)]),
                ],
            ),
            (
                "Simple decimal point numbers",
                "M.1.2 L.3.4",
                vec![line!([(0.1, 0.2), (0.3, 0.4)])],
            ),
            (
                "Decimal point numbers in cubic bezier",
                "M.123.456 C.789.012 .345.678 .901.234",
                vec![cubic!([
                    (0.123, 0.456),
                    (0.789, 0.012),
                    (0.345, 0.678),
                    (0.901, 0.234)
                ])],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_adjacent_numbers() {
        let test_cases = [
            (
                "Numbers immediately followed by negative numbers",
                "M 0,0 M 0-10.254 20-30.5",
                vec![
                    line!([(0.0, 0.0), (0.0, -10.254)]),
                    line!([(0.0, -10.254), (20.0, -30.5)]),
                ],
            ),
            (
                "Scientific notation with adjacent numbers",
                "M 0,0 M 1e-4-2e+4",
                vec![line!([(0.0, 0.0), (0.0001, -20000.0)])],
            ),
            (
                "Decimal points with adjacent numbers",
                "M 0,0 M .1-.2 .3-.4",
                vec![
                    line!([(0.0, 0.0), (0.1, -0.2)]),
                    line!([(0.1, -0.2), (0.3, -0.4)]),
                ],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_long_cubic() {
        let test_cases = [
            (
                "Simple long cubic",
                "M 10,10 c 25 3 45 1 49 -6 13 -20 24 -11 42 34",
                vec![
                    cubic!([(10.0, 10.0), (35.0, 13.0), (55.0, 11.0), (59.0, 4.0)]),
                    cubic!([(59.0, 4.0), (72.0, -16.0), (83.0, -7.0), (101.0, 38.0)]),
                ],
            ),
            (
                "Simple long cubic with newlines",
                "M 10,10 c 25 3 45 1 49 -6\n13 -20 24 -11 42 34",
                vec![
                    cubic!([(10.0, 10.0), (35.0, 13.0), (55.0, 11.0), (59.0, 4.0)]),
                    cubic!([(59.0, 4.0), (72.0, -16.0), (83.0, -7.0), (101.0, 38.0)]),
                ],
            ),
            (
                "Long cubic with inconsistent spacing",
                "M 10,10 c 25    3       45      1       49      -6\n13    -20     24      -11     42      34",
                vec![
                    cubic!([(10.0, 10.0), (35.0, 13.0), (55.0, 11.0), (59.0, 4.0)]),
                    cubic!([(59.0, 4.0), (72.0, -16.0), (83.0, -7.0), (101.0, 38.0)]),
                ],
            ),
        ];
        run_test_cases(&test_cases);
    }

    #[test]
    fn test_parse_multiple_paths() {
        let input = "M 10,10 C 20,20 40,20 50,10 Z M 30,30 C 40,40 50,50 60,60 Z";
        let curves = BezierCurve::parse_maybe_multiple(input).unwrap();

        // Verify we got two curves
        assert_eq!(curves.len(), 2, "Expected 2 curves, got {}", curves.len());

        // First path should be a closed line
        assert_eq!(
            curves[0].segments.len(),
            2,
            "First curve should have 2 segment"
        );
        assert!(
            curves[0].segments[0]
                == cubic!([(10.0, 10.0), (20.0, 20.0), (40.0, 20.0), (50.0, 10.0)])
        );
        assert!(curves[0].is_closed(), "First curve should be closed");

        // Second path should be a closed line
        assert_eq!(
            curves[1].segments.len(),
            2,
            "Second curve should have 2 segment"
        );
        assert!(
            curves[1].segments[0]
                == cubic!([(30.0, 30.0), (40.0, 40.0), (50.0, 50.0), (60.0, 60.0)])
        );
        assert!(curves[1].is_closed(), "Second curve should be closed");

        // Test case for relative coordinates with multiple subpaths
        let complex_input = "m 842.88566,3568.9387 2.92314,-5.5685 6.92066,0.8095 5.31975,-10.1339 -4.66022,-5.1155 2.85051,-5.4301 20.1684,23.183 -2.84143,5.4128 z m 14.97766,-4.1596 11.24133,1.4452 -7.6101,-8.3625 z";
        let complex_curves = BezierCurve::parse_maybe_multiple(complex_input).unwrap();

        // Verify we got two curves
        assert_eq!(
            complex_curves.len(),
            2,
            "Expected 2 curves for complex input, got {}",
            complex_curves.len()
        );

        // First path should be closed
        assert!(
            complex_curves[0].is_closed(),
            "First complex curve should be closed"
        );

        // Second path should be closed
        assert!(
            complex_curves[1].is_closed(),
            "Second complex curve should be closed"
        );
    }
}
