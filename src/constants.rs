//! Constants used throughout the library

/// Tolerance used for floating point comparisons
///
/// Used in:
/// - Point equality comparisons (`PartialEq` implementation for `Point`)
pub const FLOAT_TOLERANCE: f64 = 1e-10;
