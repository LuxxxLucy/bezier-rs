//! String format for all input/output parsing and exporting
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Format {
    SvgPath,
    Json,
}

impl Format {
    /// Try to detect the format, quite naive.
    pub fn detect(input: &str) -> Option<Self> {
        if input.trim().starts_with('[') && serde_json::from_str::<serde_json::Value>(input).is_ok()
        {
            return Some(Format::Json);
        }
        if input.trim().chars().next().map_or(false, |c| {
            matches!(c, 'M' | 'L' | 'C' | 'Q' | 'H' | 'V' | 'Z')
        }) {
            return Some(Format::SvgPath);
        }
        None
    }
}

impl FromStr for Format {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "svg" | "svgpath" | "path" => Ok(Format::SvgPath),
            "json" => Ok(Format::Json),
            _ => Err(format!("Unknown format: {}. Expected 'svg' or 'json'", s)),
        }
    }
}
