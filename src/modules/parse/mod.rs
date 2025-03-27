//! Parsing module for bezier curves
//!
//! Now supported format:
//! - JSON:
//!     in the form of `[{"x": 0.0, "y": 0.0, "on": true}, {"x": 1.0, "y": 1.0, "on": false}, {"x": 2.0, "y": 0.0, "on": true}]`.
//!     See the `json` module for more detailed information on the JSON format.

pub mod json;
