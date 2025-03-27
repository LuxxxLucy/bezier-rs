// Module definitions
pub mod constants;
pub mod data;
pub mod error;
pub mod modules;

// export the core data structure at crate level
pub use data::curve::BezierCurve;
pub use data::point::Point;
pub use data::segment::BezierSegment;
pub use error::{BezierError, BezierResult};
