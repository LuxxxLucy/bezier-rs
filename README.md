# bezier-rs

A Rust library for working with Bezier curves.

## Features

- Data structures for representing Bezier curves (cubic and quadratic segments)
- Fitting curves
    - least square method (see [Cubic Bezier Fitting with Least Squares](https://luxxxlucy.github.io/projects/2025_bezielogue_1_cubic_fitting/index.html))
      - simple least square
      - Alternating method (nearest point and gauss-newton)
      - an improvised gauss-newton with lineseach
