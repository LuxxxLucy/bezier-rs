from z3 import *
from typing import List, Tuple

def bezier_x(t: ArithRef, points: List[Tuple[float, float]]) -> ArithRef:
    """Calculate x coordinate on cubic Bezier curve at parameter t."""
    p0, p1, p2, p3 = points
    return (1 - t)**3 * p0[0] + 3 * t * (1 - t)**2 * p1[0] + 3 * t**2 * (1 - t) * p2[0] + t**3 * p3[0]

def bezier_y(t: ArithRef, points: List[Tuple[float, float]]) -> ArithRef:
    """Calculate y coordinate on cubic Bezier curve at parameter t."""
    p0, p1, p2, p3 = points
    return (1 - t)**3 * p0[1] + 3 * t * (1 - t)**2 * p1[1] + 3 * t**2 * (1 - t) * p2[1] + t**3 * p3[1]

def main():
    s = Solver()

    # Define control points for two curves, each with two segments
    # First curve (p) with two segments
    p_segments = [
        # First segment: points (0,0) -> (1,2) -> (2,2) -> (3,0)
        [
            (0.0, 0.0),  # p0
            (1.0, 2.0),  # p1
            (2.0, 2.0),  # p2
            (3.0, 0.0)   # p3
        ],
        # Second segment: points (3,0) -> (2,2) -> (1,2) -> (0,0)
        [
            (3.0, 0.0),  # p0
            (2.0, 2.0),  # p1
            (1.0, 2.0),  # p2
            (0.0, 0.0)   # p3
        ]
    ]

    # Second curve (q) with two segments
    q_segments = [
        # First segment: points (0,0) -> (1,2) -> (2,2) -> (3,0)
        [
            (0.0, 0.0),  # q0
            (1.0, 2.0),  # q1
            (2.0, 2.0),  # q2
            (3.0, 0.0)   # q3
        ],
        # Second segment: points (3,0) -> (2,2) -> (1,2) -> (0,0)
        [
            (3.0, 0.0),  # q0
            (2.0, 2.0),  # q1
            (1.0, 2.0),  # q2
            (0.0, 0.0)   # q3
        ]
    ]

    # Create parameters
    t1, t2 = Real('t1'), Real('t2')

    # For first segment of p, check if it matches with either segment of q
    s.add(ForAll([t1],
        Implies(
            And(t1 >= 0.0, t1 <= 1.0),
            Exists([t2],
                And(
                    t2 >= 0.0, t2 <= 1.0,
                    Or(
                        # Check against first segment of q
                        And(
                            bezier_x(t1, p_segments[0]) == bezier_x(t2, q_segments[0]),
                            bezier_y(t1, p_segments[0]) == bezier_y(t2, q_segments[0])
                        ),
                        # Check against second segment of q
                        And(
                            bezier_x(t1, p_segments[0]) == bezier_x(t2, q_segments[1]),
                            bezier_y(t1, p_segments[0]) == bezier_y(t2, q_segments[1])
                        )
                    )
                )
            )
        )
    ))

    # For second segment of p, check if it matches with either segment of q
    s.add(ForAll([t1],
        Implies(
            And(t1 >= 0.0, t1 <= 1.0),
            Exists([t2],
                And(
                    t2 >= 0.0, t2 <= 1.0,
                    Or(
                        # Check against first segment of q
                        And(
                            bezier_x(t1, p_segments[1]) == bezier_x(t2, q_segments[0]),
                            bezier_y(t1, p_segments[1]) == bezier_y(t2, q_segments[0])
                        ),
                        # Check against second segment of q
                        And(
                            bezier_x(t1, p_segments[1]) == bezier_x(t2, q_segments[1]),
                            bezier_y(t1, p_segments[1]) == bezier_y(t2, q_segments[1])
                        )
                    )
                )
            )
        )
    ))

    # For first segment of q, check if it matches with either segment of p
    s.add(ForAll([t1],
        Implies(
            And(t1 >= 0.0, t1 <= 1.0),
            Exists([t2],
                And(
                    t2 >= 0.0, t2 <= 1.0,
                    Or(
                        # Check against first segment of p
                        And(
                            bezier_x(t1, q_segments[0]) == bezier_x(t2, p_segments[0]),
                            bezier_y(t1, q_segments[0]) == bezier_y(t2, p_segments[0])
                        ),
                        # Check against second segment of p
                        And(
                            bezier_x(t1, q_segments[0]) == bezier_x(t2, p_segments[1]),
                            bezier_y(t1, q_segments[0]) == bezier_y(t2, p_segments[1])
                        )
                    )
                )
            )
        )
    ))

    # For second segment of q, check if it matches with either segment of p
    s.add(ForAll([t1],
        Implies(
            And(t1 >= 0.0, t1 <= 1.0),
            Exists([t2],
                And(
                    t2 >= 0.0, t2 <= 1.0,
                    Or(
                        # Check against first segment of p
                        And(
                            bezier_x(t1, q_segments[1]) == bezier_x(t2, p_segments[0]),
                            bezier_y(t1, q_segments[1]) == bezier_y(t2, p_segments[0])
                        ),
                        # Check against second segment of p
                        And(
                            bezier_x(t1, q_segments[1]) == bezier_x(t2, p_segments[1]),
                            bezier_y(t1, q_segments[1]) == bezier_y(t2, p_segments[1])
                        )
                    )
                )
            )
        )
    ))

    # Check and print results
    result = s.check()
    print(f"Result: {result}")

if __name__ == "__main__":
    main()
