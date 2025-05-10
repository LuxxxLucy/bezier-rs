//! Export Bezier curves to JSON format

use crate::data::BezierCurve;

pub trait ToJson {
    /// Export the curve as a JSON string (array of points)
    fn to_json(&self) -> serde_json::Result<String>;
}

impl ToJson for BezierCurve {
    fn to_json(&self) -> serde_json::Result<String> {
        let points: Vec<_> = self
            .segments
            .iter()
            .flat_map(|segment| {
                let points = segment.points();
                points.into_iter().map(|point| {
                    serde_json::json!({
                        "x": point.x,
                        "y": point.y,
                        "on": true
                    })
                })
            })
            .collect();
        serde_json::to_string(&points)
    }
}
