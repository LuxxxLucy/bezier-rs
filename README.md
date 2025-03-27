# bezier-rs

A Rust library for working with Bezier curves.

## Features

- Data structures for representing Bezier curves (cubic and quadratic segments)
## Usage

### Basic usage

```rust
use bezier_rs::{BezierSegment, Point};

// Create a cubic Bezier segment
let segment = BezierSegment::new_cubic(
    Point::new(0.0, 0.0),   // Start point
    Point::new(1.0, 2.0),   // Control point 1
    Point::new(2.0, 2.0),   // Control point 2
    Point::new(3.0, 0.0),   // End point
);

// Get a point on the curve at parameter t (0.0 <= t <= 1.0)
let point = segment.point_at(0.5);
println!("Point at t=0.5: ({}, {})", point.x, point.y);

// Sample multiple points along the curve
let points = segment.sample_points(10);
```

or using macro shorthand

```rust
use bezier_rs::{cubic, quad, curve, curve_from};

// Create a cubic Bezier segment with tuples
let cubic_segment = cubic!([
    (0.0, 0.0),    // Start point
    (1.0, 2.0),    // Control point 1
    (2.0, 2.0),    // Control point 2
    (3.0, 0.0)     // End point
]);

// Create a quadratic Bezier segment
let quad_segment = quad!([
    (0.0, 0.0),    // Start point
    (1.0, 2.0),    // Control point
    (2.0, 0.0)     // End point
]);

// Create a curve from a single segment
let curve1 = curve_from!(cubic_segment);

// Create a curve with multiple segments
let curve2 = curve!([cubic_segment, quad_segment]);
```

### parsing from json

You can parse Bezier curves from JSON input:

```rust
use bezier_rs::modules::parse::json;

let json = r#"[
    {"x": 0.0, "y": 0.0, "on": true},
    {"x": 1.0, "y": 1.0, "on": false},
    {"x": 2.0, "y": 1.0, "on": false},
    {"x": 3.0, "y": 0.0, "on": true}
]"#;

let curve = json::parse(json).unwrap();
```

### SVG export

```rust
use bezier_rs::{cubic, curve_from, modules::export};

// Create a curve
let curve = curve_from!(cubic!([
    (50, 200),
    (100, 50),
    (200, 50),
    (250, 200)
]));

// Export to SVG
let svg = export::svg::to_svg(&curve, 300, 200);
std::fs::write("curve.svg", &svg).unwrap();
```

