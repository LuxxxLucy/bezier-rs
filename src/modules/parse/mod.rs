//! Parsing module for bezier curves
//!
//! Now supported format:
//! - JSON:
//!     in the form of `[{"x": 0.0, "y": 0.0, "on": true}, {"x": 1.0, "y": 1.0, "on": false}, {"x": 2.0, "y": 0.0, "on": true}]`.
//!     See the `json` module for more detailed information on the JSON format.
//! - SVG:
//!     Parse SVG paths and convert them to bezier curves.

pub mod json;
pub mod svg_path;

use crate::data::format::Format;
use crate::data::BezierCurve;
use crate::error::{BezierError, BezierResult};
use crate::modules::parse::svg_path::FromSvgPath;

/// Parse input data into a BezierCurve
///
/// # Arguments
///
/// * `input` - The input string to parse
/// * `format` - Optional format specification. If None, will attempt to auto-detect
///
/// # Returns
///
/// A Result containing either the parsed BezierCurve or an error
///
/// # Examples
///
/// ```
/// use bezier_rs::modules::parse::parse;
/// use bezier_rs::data::format::Format;
///
/// // Parse with explicit format
/// let svg_input = "M10 10 L20 20";
/// let curve = parse(svg_input, Some(Format::SvgPath)).unwrap();
///
/// // Parse with auto-detection
/// let json_input = r#"[{"x":0,"y":0,"on":true},{"x":10,"y":10,"on":false},{"x":20,"y":20,"on":true}]"#;
/// let curve = parse(json_input, None).unwrap();
/// ```
pub fn parse(input: &str, format: Option<Format>) -> BezierResult<BezierCurve> {
    let format = match format {
        Some(f) => f,
        None => Format::detect(input).ok_or_else(|| {
            BezierError::ParseError(
                "Could not detect input format. Please specify format explicitly.".to_string(),
            )
        })?,
    };

    match format {
        Format::SvgPath => {
            BezierCurve::from_svg_path(input).map_err(|e| BezierError::ParseError(e.to_string()))
        }
        Format::Json => json::parse(input),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_detection() {
        // Test JSON detection
        let json_input =
            r#"[{"x":0,"y":0,"on":true},{"x":10,"y":10,"on":false},{"x":20,"y":20,"on":true}]"#;
        assert_eq!(Format::detect(json_input), Some(Format::Json));

        // Test SVG detection
        let svg_input = "M10 10 L20 20";
        assert_eq!(Format::detect(svg_input), Some(Format::SvgPath));

        // Test unknown format
        let unknown_input = "not a valid format";
        assert_eq!(Format::detect(unknown_input), None);
    }

    #[test]
    fn test_parse_with_detection() {
        // Test JSON parsing (ends with on-curve point)
        let json_input =
            r#"[{"x":0,"y":0,"on":true},{"x":10,"y":10,"on":false},{"x":20,"y":20,"on":true}]"#;
        let curve = parse(json_input, None).unwrap();
        assert_eq!(curve.segments.len(), 1);

        // Test SVG parsing
        let svg_input = "M10 10 L20 20";
        let curve = parse(svg_input, None).unwrap();
        assert_eq!(curve.segments.len(), 1);
    }

    #[test]
    fn test_parse_with_explicit_format() {
        // Test JSON parsing with explicit format (ends with on-curve point)
        let json_input =
            r#"[{"x":0,"y":0,"on":true},{"x":10,"y":10,"on":false},{"x":20,"y":20,"on":true}]"#;
        let curve = parse(json_input, Some(Format::Json)).unwrap();
        assert_eq!(curve.segments.len(), 1);

        // Test SVG parsing with explicit format
        let svg_input = "M10 10 L20 20";
        let curve = parse(svg_input, Some(Format::SvgPath)).unwrap();
        assert_eq!(curve.segments.len(), 1);
    }
}
