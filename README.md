# bezier-rs

A Rust library for working with Bezier curves.

## Features

- Data structures for representing Bezier curves (cubic and quadratic segments)
- parse/export:
  - SVG path `M 10,20 C 20,30 30,40 40,50 Z`.
  - json format ` {"x": 0.0, "y": 0.0, "on": true}, {"x": 1.0, "y": 1.0, "on": false}, {"x": 2.0, "y": 0.0, "on": true}`
- Fitting curves
    - least square method (see [Cubic Bezier Fitting with Least Squares](https://luxxxlucy.github.io/projects/2025_bezielogue_1_cubic_fitting/index.html))
      - simple least square
      - Alternating method (nearest point and gauss-newton)
      - an improvised gauss-newton with lineseach
