//! Export Bezier curves to various formats
//!
//! This module provides functionality to export bezier curves to different
//! formats for visualization, sharing, or further processing.
//!
//! # Available Export Formats
//! Supported formats:
//! - SVG: Export as SVG path string
//! - [SVG](svg/index.html) - Export curves and points to SVG format
//! - [Json](json/index.html) - Export curves and points to json format

pub mod json;
pub mod svg_path;

use crate::data::format::Format;
use crate::data::BezierCurve;
use crate::error::BezierResult;
use crate::modules::export::json::ToJson;
use crate::modules::export::svg_path::ToSvgPath;

/// Export a BezierCurve to the specified format
///
/// # Arguments
///
/// * `curve` - The BezierCurve to export
/// * `format` - The format to export to
///
/// # Returns
///
/// A Result containing either the exported string or an error
///
/// # Examples
///
/// ```
/// use bezier_rs::modules::parse::parse;
/// use bezier_rs::modules::export::export;
/// use bezier_rs::data::format::Format;
///
/// let json_input = r#"[{"x":0,"y":0,"on":true},{"x":10,"y":10,"on":false},{"x":20,"y":20,"on":true}]"#;
/// let curve = parse(json_input, Some(Format::Json)).unwrap();
/// let svg = export(&curve, Format::SvgPath).unwrap();
/// let json = export(&curve, Format::Json).unwrap();
/// ```
pub fn export(curve: &BezierCurve, format: Format) -> BezierResult<String> {
    match format {
        Format::SvgPath => Ok(curve.to_svg_path()),
        Format::Json => Ok(curve.to_json()?),
    }
}
