;; Define function for cubic Bezier curve coordinates
(define-fun bezier_cubic ((t Real) (p0 Real) (p1 Real) (p2 Real) (p3 Real)) Real
    (+
        (* (* (* (- 1 t) (- 1 t)) (- 1 t)) p0)
        (* (* (* (* 3 t) (- 1 t)) (- 1 t)) p1)
        (* (* (* (* 3 t) t) (- 1 t)) p2)
        (* (* (* t t) t) p3)))

;; Define control points for first curve (p) with two segments
;; First segment: points (0,0) -> (1,2) -> (2,2) -> (3,0)
(declare-const p0_x Real) (declare-const p0_y Real)
(declare-const p1_x Real) (declare-const p1_y Real)
(declare-const p2_x Real) (declare-const p2_y Real)
(declare-const p3_x Real) (declare-const p3_y Real)

;; Second segment: points (3,0) -> (2,2) -> (1,2) -> (0,0)
(declare-const p4_x Real) (declare-const p4_y Real)
(declare-const p5_x Real) (declare-const p5_y Real)
(declare-const p6_x Real) (declare-const p6_y Real)
(declare-const p7_x Real) (declare-const p7_y Real)

;; Set control points for first curve
(assert (and
    (= p0_x 0.0) (= p0_y 0.0)
    (= p1_x 1.0) (= p1_y 2.0)
    (= p2_x 2.0) (= p2_y 2.0)
    (= p3_x 3.0) (= p3_y 0.0)
    (= p4_x 3.0) (= p4_y 0.0)
    (= p5_x 2.0) (= p5_y 2.0)
    (= p6_x 1.0) (= p6_y 2.0)
    (= p7_x 0.0) (= p7_y 0.0)))

;; Define control points for second curve (q) with two segments
;; First segment: points (0,0) -> (1,2) -> (2,2) -> (3,0)
(declare-const q0_x Real) (declare-const q0_y Real)
(declare-const q1_x Real) (declare-const q1_y Real)
(declare-const q2_x Real) (declare-const q2_y Real)
(declare-const q3_x Real) (declare-const q3_y Real)

;; Second segment: points (3,0) -> (2,2) -> (1,2) -> (0,0)
(declare-const q4_x Real) (declare-const q4_y Real)
(declare-const q5_x Real) (declare-const q5_y Real)
(declare-const q6_x Real) (declare-const q6_y Real)
(declare-const q7_x Real) (declare-const q7_y Real)

;; Set control points for second curve
(assert (and
    (= q0_x 0.0) (= q0_y 0.0)
    (= q1_x 1.0) (= q1_y 2.0)
    (= q2_x 2.0) (= q2_y 2.0)
    (= q3_x 3.0) (= q3_y 0.0)
    (= q4_x 3.0) (= q4_y 0.0)
    (= q5_x 2.0) (= q5_y 2.0)
    (= q6_x 1.0) (= q6_y 2.0)
    (= q7_x 0.0) (= q7_y 0.0)))

;; For first segment of p, check if it matches with either segment of q
(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((p_x (bezier_cubic t1 p0_x p1_x p2_x p3_x))
                      (p_y (bezier_cubic t1 p0_y p1_y p2_y p3_y)))
                    (or
                        ;; Check against first segment of q
                        (and
                            (= p_x (bezier_cubic t2 q0_x q1_x q2_x q3_x))
                            (= p_y (bezier_cubic t2 q0_y q1_y q2_y q3_y)))
                        ;; Check against second segment of q
                        (and
                            (= p_x (bezier_cubic t2 q4_x q5_x q6_x q7_x))
                            (= p_y (bezier_cubic t2 q4_y q5_y q6_y q7_y))))))))))

;; For second segment of p, check if it matches with either segment of q
(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((p_x (bezier_cubic t1 p4_x p5_x p6_x p7_x))
                      (p_y (bezier_cubic t1 p4_y p5_y p6_y p7_y)))
                    (or
                        ;; Check against first segment of q
                        (and
                            (= p_x (bezier_cubic t2 q0_x q1_x q2_x q3_x))
                            (= p_y (bezier_cubic t2 q0_y q1_y q2_y q3_y)))
                        ;; Check against second segment of q
                        (and
                            (= p_x (bezier_cubic t2 q4_x q5_x q6_x q7_x))
                            (= p_y (bezier_cubic t2 q4_y q5_y q6_y q7_y))))))))))

;; For first segment of q, check if it matches with either segment of p
(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((q_x (bezier_cubic t1 q0_x q1_x q2_x q3_x))
                      (q_y (bezier_cubic t1 q0_y q1_y q2_y q3_y)))
                    (or
                        ;; Check against first segment of p
                        (and
                            (= q_x (bezier_cubic t2 p0_x p1_x p2_x p3_x))
                            (= q_y (bezier_cubic t2 p0_y p1_y p2_y p3_y)))
                        ;; Check against second segment of p
                        (and
                            (= q_x (bezier_cubic t2 p4_x p5_x p6_x p7_x))
                            (= q_y (bezier_cubic t2 p4_y p5_y p6_y p7_y))))))))))

;; For second segment of q, check if it matches with either segment of p
(assert (forall ((t1 Real))
    (=>
        (and (>= t1 0.0) (<= t1 1.0))
        (exists ((t2 Real))
            (and
                (>= t2 0.0) (<= t2 1.0)
                (let ((q_x (bezier_cubic t1 q4_x q5_x q6_x q7_x))
                      (q_y (bezier_cubic t1 q4_y q5_y q6_y q7_y)))
                    (or
                        ;; Check against first segment of p
                        (and
                            (= q_x (bezier_cubic t2 p0_x p1_x p2_x p3_x))
                            (= q_y (bezier_cubic t2 p0_y p1_y p2_y p3_y)))
                        ;; Check against second segment of p
                        (and
                            (= q_x (bezier_cubic t2 p4_x p5_x p6_x p7_x))
                            (= q_y (bezier_cubic t2 p4_y p5_y p6_y p7_y))))))))))

(check-sat) 