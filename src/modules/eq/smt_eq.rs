//! SMT-based equality checking for Bezier curves.
//!
//! given two curves, we export a smtlib-2 format and then run z3 to check if the curves are equal.

use crate::data::curve::BezierCurve;
use crate::data::segment::BezierSegment;
use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Compare two Bezier curves for equality
///
/// given two curves, we export a smtlib-2 format and then run z3 to check if the curves are equal.
///
/// # Arguments
///
/// * `curve1` - The first curve to compare
/// * `curve2` - The second curve to compare
pub fn smt_eq(curve1: &BezierCurve, curve2: &BezierCurve) -> bool {
    let smt = create_smt_lib(curve1, curve2);

    verify(smt)
}

fn create_smt_lib(curve1: &BezierCurve, curve2: &BezierCurve) -> String {
    // check all segments are cubic
    pre_check(curve1, curve2);

    // p_values: control points from curve1
    let mut p_values = Vec::new();
    for segment in &curve1.segments {
        if let BezierSegment::Cubic { points } = segment {
            for &point in points.iter() {
                p_values.push((point.x, point.y));
            }
        }
    }
    // q_values: control points from curve2
    let mut q_values = Vec::new();
    for segment in &curve2.segments {
        if let BezierSegment::Cubic { points } = segment {
            for &point in points.iter() {
                q_values.push((point.x, point.y));
            }
        }
    }

    let mut smt = String::new();
    smt.push_str(&smt_bezier_cubic_def());
    smt.push('\n');
    smt.push_str(&smt_curve_point_decls(curve1, "p"));
    smt.push('\n');
    smt.push_str(&smt_curve_point_asserts(curve1, "p", &p_values));
    smt.push('\n');
    smt.push_str(&smt_curve_point_decls(curve2, "q"));
    smt.push('\n');
    smt.push_str(&smt_curve_point_asserts(curve2, "q", &q_values));
    smt.push('\n');
    smt.push_str(&add_smt_assert_with_prefix(
        "p",
        "q",
        0,
        curve2.segments.len(),
    ));
    smt.push('\n');
    smt.push_str(&add_smt_assert_with_prefix(
        "p",
        "q",
        1,
        curve2.segments.len(),
    ));
    smt.push('\n');
    smt.push_str(&add_smt_assert_with_prefix(
        "q",
        "p",
        0,
        curve1.segments.len(),
    ));
    smt.push('\n');
    smt.push_str(&add_smt_assert_with_prefix(
        "q",
        "p",
        1,
        curve1.segments.len(),
    ));
    smt.push('\n');
    smt.push_str("(check-sat)");
    smt
}

fn verify(smt: String) -> bool {
    // Create a temporary file
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    file.write_all(smt.as_bytes())
        .expect("Failed to write SMT to file");
    let path = file.path();

    // Run z3 on the file
    let output = Command::new("z3").arg("-smt2").arg(path).output();

    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout);
            stdout.contains("sat")
        }
        Err(_) => false,
    }
}

// check all segments are cubic
fn pre_check(curve1: &BezierCurve, curve2: &BezierCurve) {
    for (i, segment) in curve1.segments.iter().enumerate() {
        if !matches!(segment, BezierSegment::Cubic { .. }) {
            panic!(
                "curve1: segment {} is not cubic. Here only cubic segments are supported.",
                i
            );
        }
    }
    for (i, segment) in curve2.segments.iter().enumerate() {
        if !matches!(segment, BezierSegment::Cubic { .. }) {
            panic!(
                "curve2: segment {} is not cubic. Here only cubic segments are supported.",
                i
            );
        }
    }
}

fn smt_bezier_cubic_def() -> String {
    let mut s = String::new();
    s.push_str(
        "(define-fun bezier_cubic ((t Real) (p0 Real) (p1 Real) (p2 Real) (p3 Real)) Real\n",
    );
    s.push_str("    (+\n");
    s.push_str("        (* (* (* (- 1 t) (- 1 t)) (- 1 t)) p0)\n");
    s.push_str("        (* (* (* (* 3 t) (- 1 t)) (- 1 t)) p1)\n");
    s.push_str("        (* (* (* (* 3 t) t) (- 1 t)) p2)\n");
    s.push_str("        (* (* (* t t) t) p3)))\n");
    s
}

fn smt_curve_point_decls(curve: &BezierCurve, prefix: &str) -> String {
    let mut s = String::new();
    for (i, segment) in curve.segments.iter().enumerate() {
        if let BezierSegment::Cubic { points } = segment {
            for (j, _) in points.iter().enumerate() {
                s.push_str(&format!(
                    "(declare-const {prefix}{}_x Real) (declare-const {prefix}{}_y Real)\n",
                    i * 4 + j,
                    i * 4 + j,
                    prefix = prefix
                ));
            }
        }
    }
    s
}

fn smt_curve_point_asserts(_curve: &BezierCurve, prefix: &str, values: &[(f64, f64)]) -> String {
    let mut s = String::new();
    s.push_str("(assert (and\n");
    for (i, &(x, y)) in values.iter().enumerate() {
        s.push_str(&format!(
            "    (= {prefix}{}_x {:.1}) (= {prefix}{}_y {:.1})",
            i,
            x,
            i,
            y,
            prefix = prefix
        ));
        s.push('\n');
    }
    s.push_str("))\n");
    s
}

fn add_smt_assert_with_prefix(
    prefix: &str,
    other_prefix: &str,
    segment_idx: usize,
    num_segments: usize,
) -> String {
    let mut s = String::new();
    let p_offset = segment_idx * 4;
    s.push_str("(assert (forall ((t1 Real))\n");
    s.push_str("    (=>\n");
    s.push_str("        (and (>= t1 0.0) (<= t1 1.0))\n");
    s.push_str("        (exists ((t2 Real))\n");
    s.push_str("            (and\n");
    s.push_str("                (>= t2 0.0) (<= t2 1.0)\n");
    s.push_str(&format!(
        "                (let (({}_x (bezier_cubic t1 {}{}_x {}{}_x {}{}_x {}{}_x))\n",
        prefix,
        prefix,
        p_offset,
        prefix,
        p_offset + 1,
        prefix,
        p_offset + 2,
        prefix,
        p_offset + 3
    ));
    s.push_str(&format!(
        "                      ({}_y (bezier_cubic t1 {}{}_y {}{}_y {}{}_y {}{}_y)))\n",
        prefix,
        prefix,
        p_offset,
        prefix,
        p_offset + 1,
        prefix,
        p_offset + 2,
        prefix,
        p_offset + 3
    ));
    s.push_str("                    (or\n");

    // Generate assertions for each segment of the other curve
    for i in 0..num_segments {
        let q_offset = i * 4;
        s.push_str("                        (and\n");
        s.push_str(&format!(
            "                            (= {}_x (bezier_cubic t2 {}{}_x {}{}_x {}{}_x {}{}_x))\n",
            prefix,
            other_prefix,
            q_offset,
            other_prefix,
            q_offset + 1,
            other_prefix,
            q_offset + 2,
            other_prefix,
            q_offset + 3
        ));
        s.push_str(&format!(
            "                            (= {}_y (bezier_cubic t2 {}{}_y {}{}_y {}{}_y {}{}_y))\n",
            prefix,
            other_prefix,
            q_offset,
            other_prefix,
            q_offset + 1,
            other_prefix,
            q_offset + 2,
            other_prefix,
            q_offset + 3
        ));
        s.push_str("                        )");
        if i < num_segments - 1 {
            s.push('\n');
        }
    }

    s.push_str("\n                    )\n");
    s.push_str("                )\n");
    s.push_str("            )\n");
    s.push_str("        )\n");
    s.push_str("    )\n");
    s.push_str("))\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{cubic, curve};

    #[test]
    fn test_smt_eq_curve_point_decls() {
        let segment1 = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);
        let segment2 = cubic!([(3.0, 0.0), (2.0, 2.0), (1.0, 2.0), (0.0, 0.0)]);
        let curve = curve!([segment1, segment2]);
        let expected = "(declare-const p0_x Real) (declare-const p0_y Real)\n(declare-const p1_x Real) (declare-const p1_y Real)\n(declare-const p2_x Real) (declare-const p2_y Real)\n(declare-const p3_x Real) (declare-const p3_y Real)\n(declare-const p4_x Real) (declare-const p4_y Real)\n(declare-const p5_x Real) (declare-const p5_y Real)\n(declare-const p6_x Real) (declare-const p6_y Real)\n(declare-const p7_x Real) (declare-const p7_y Real)\n";
        assert_eq!(smt_curve_point_decls(&curve, "p"), expected);
    }

    #[test]
    fn test_smt_eq_add_smt_assert() {
        let expected = r#"(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((p_x (bezier_cubic t1 p4_x p5_x p6_x p7_x))
                      (p_y (bezier_cubic t1 p4_y p5_y p6_y p7_y)))
                    (or
                        (and
                            (= p_x (bezier_cubic t2 q0_x q1_x q2_x q3_x))
                            (= p_y (bezier_cubic t2 q0_y q1_y q2_y q3_y))
                        )
                        (and
                            (= p_x (bezier_cubic t2 q4_x q5_x q6_x q7_x))
                            (= p_y (bezier_cubic t2 q4_y q5_y q6_y q7_y))
                        )
                    )
                )
            )
        )
    )
))
"#;

        let actual = add_smt_assert_with_prefix("p", "q", 1, 2);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_smt_eq_bezier_curves_output() {
        let segment1 = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);
        let segment2 = cubic!([(3.0, 0.0), (2.0, 2.0), (1.0, 2.0), (0.0, 0.0)]);
        let curve1 = curve!([segment1.clone(), segment2.clone()]);
        let curve2 = curve!([segment1, segment2]);

        let expected = r#"(define-fun bezier_cubic ((t Real) (p0 Real) (p1 Real) (p2 Real) (p3 Real)) Real
    (+
        (* (* (* (- 1 t) (- 1 t)) (- 1 t)) p0)
        (* (* (* (* 3 t) (- 1 t)) (- 1 t)) p1)
        (* (* (* (* 3 t) t) (- 1 t)) p2)
        (* (* (* t t) t) p3)))

(declare-const p0_x Real) (declare-const p0_y Real)
(declare-const p1_x Real) (declare-const p1_y Real)
(declare-const p2_x Real) (declare-const p2_y Real)
(declare-const p3_x Real) (declare-const p3_y Real)
(declare-const p4_x Real) (declare-const p4_y Real)
(declare-const p5_x Real) (declare-const p5_y Real)
(declare-const p6_x Real) (declare-const p6_y Real)
(declare-const p7_x Real) (declare-const p7_y Real)

(assert (and
    (= p0_x 0.0) (= p0_y 0.0)
    (= p1_x 1.0) (= p1_y 2.0)
    (= p2_x 2.0) (= p2_y 2.0)
    (= p3_x 3.0) (= p3_y 0.0)
    (= p4_x 3.0) (= p4_y 0.0)
    (= p5_x 2.0) (= p5_y 2.0)
    (= p6_x 1.0) (= p6_y 2.0)
    (= p7_x 0.0) (= p7_y 0.0)
))

(declare-const q0_x Real) (declare-const q0_y Real)
(declare-const q1_x Real) (declare-const q1_y Real)
(declare-const q2_x Real) (declare-const q2_y Real)
(declare-const q3_x Real) (declare-const q3_y Real)
(declare-const q4_x Real) (declare-const q4_y Real)
(declare-const q5_x Real) (declare-const q5_y Real)
(declare-const q6_x Real) (declare-const q6_y Real)
(declare-const q7_x Real) (declare-const q7_y Real)

(assert (and
    (= q0_x 0.0) (= q0_y 0.0)
    (= q1_x 1.0) (= q1_y 2.0)
    (= q2_x 2.0) (= q2_y 2.0)
    (= q3_x 3.0) (= q3_y 0.0)
    (= q4_x 3.0) (= q4_y 0.0)
    (= q5_x 2.0) (= q5_y 2.0)
    (= q6_x 1.0) (= q6_y 2.0)
    (= q7_x 0.0) (= q7_y 0.0)
))

(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((p_x (bezier_cubic t1 p0_x p1_x p2_x p3_x))
                      (p_y (bezier_cubic t1 p0_y p1_y p2_y p3_y)))
                    (or
                        (and
                            (= p_x (bezier_cubic t2 q0_x q1_x q2_x q3_x))
                            (= p_y (bezier_cubic t2 q0_y q1_y q2_y q3_y))
                        )
                        (and
                            (= p_x (bezier_cubic t2 q4_x q5_x q6_x q7_x))
                            (= p_y (bezier_cubic t2 q4_y q5_y q6_y q7_y))
                        )
                    )
                )
            )
        )
    )
))

(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((p_x (bezier_cubic t1 p4_x p5_x p6_x p7_x))
                      (p_y (bezier_cubic t1 p4_y p5_y p6_y p7_y)))
                    (or
                        (and
                            (= p_x (bezier_cubic t2 q0_x q1_x q2_x q3_x))
                            (= p_y (bezier_cubic t2 q0_y q1_y q2_y q3_y))
                        )
                        (and
                            (= p_x (bezier_cubic t2 q4_x q5_x q6_x q7_x))
                            (= p_y (bezier_cubic t2 q4_y q5_y q6_y q7_y))
                        )
                    )
                )
            )
        )
    )
))

(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((q_x (bezier_cubic t1 q0_x q1_x q2_x q3_x))
                      (q_y (bezier_cubic t1 q0_y q1_y q2_y q3_y)))
                    (or
                        (and
                            (= q_x (bezier_cubic t2 p0_x p1_x p2_x p3_x))
                            (= q_y (bezier_cubic t2 p0_y p1_y p2_y p3_y))
                        )
                        (and
                            (= q_x (bezier_cubic t2 p4_x p5_x p6_x p7_x))
                            (= q_y (bezier_cubic t2 p4_y p5_y p6_y p7_y))
                        )
                    )
                )
            )
        )
    )
))

(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((q_x (bezier_cubic t1 q4_x q5_x q6_x q7_x))
                      (q_y (bezier_cubic t1 q4_y q5_y q6_y q7_y)))
                    (or
                        (and
                            (= q_x (bezier_cubic t2 p0_x p1_x p2_x p3_x))
                            (= q_y (bezier_cubic t2 p0_y p1_y p2_y p3_y))
                        )
                        (and
                            (= q_x (bezier_cubic t2 p4_x p5_x p6_x p7_x))
                            (= q_y (bezier_cubic t2 p4_y p5_y p6_y p7_y))
                        )
                    )
                )
            )
        )
    )
))

(check-sat)"#;

        let actual = create_smt_lib(&curve1, &curve2);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_smt_eq_verify() {
        let segment1 = cubic!([(0.0, 0.0), (1.0, 2.0), (2.0, 2.0), (3.0, 0.0)]);
        let segment2 = cubic!([(3.0, 0.0), (2.0, 2.0), (1.0, 2.0), (0.0, 0.0)]);
        let curve1 = curve!([segment1.clone(), segment2.clone()]);
        let curve2 = curve!([segment1, segment2]);

        let smt = create_smt_lib(&curve1, &curve2);
        assert!(
            verify(smt),
            "Verification should return true for identical curves"
        );
    }
}
