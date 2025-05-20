use crate::data::curve::BezierCurve;

/// Compare two Bezier curves for equality
pub fn smt_eq(curve1: &BezierCurve, curve2: &BezierCurve) -> bool {
    curve1 == curve2
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quad;

    #[test]
    fn test_smt_eq() {
        let segment = quad!([(0, 0), (1, 1), (2, 2)]);
        let curve1 = curve!([segment.clone()]);
        let curve2 = curve1.clone();
        
        assert!(smt_eq(&curve1, &curve2));
    }
}
