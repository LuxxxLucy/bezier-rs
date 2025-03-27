//! Error types for the bezier-rs crate

use std::error::Error;
use std::fmt;

/// Common error type for bezier-rs crate
#[derive(Debug)]
pub enum BezierError {
    /// Error occurred while parsing data
    ParseError(String),
    /// Error occurred during fit operations
    FitError(String),
    /// Generic error
    Other(String),
}

impl fmt::Display for BezierError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BezierError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            BezierError::FitError(msg) => write!(f, "Fit error: {}", msg),
            BezierError::Other(msg) => write!(f, "Error: {}", msg),
        }
    }
}

impl Error for BezierError {}

/// Result type that uses BezierError as the error type
pub type BezierResult<T> = Result<T, BezierError>;
