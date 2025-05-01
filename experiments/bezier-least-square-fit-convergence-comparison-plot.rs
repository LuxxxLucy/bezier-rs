use bezier_rs::{
    cubic,
    modules::fit::{
        alternating_least_square_fit::{fit_cubic_bezier_alternating, TUpdateMethod},
        least_square_fit_weak_varpro::fit_cubic_bezier_weak_varpro,
    },
};
use plotters::prelude::*;
use rand::prelude::*;
use rand::rngs::StdRng;
use std::fs;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create output directory if it doesn't exist
    let output_dir = Path::new("demo_temporary");
    if !output_dir.exists() {
        fs::create_dir_all(output_dir)?;
    }

    // Create a complex curve for testing
    let original = cubic!([(0.0, 0.0), (1.0, 3.0), (2.0, -1.0), (3.0, 2.0)]);

    // Define some interesting t values for sampling
    let t_values = vec![0.0, 0.1, 0.15, 0.3, 0.7, 0.85, 1.0];

    // Sample points at the specified t values
    let (samples, _) = original.sample_at_t_values(&t_values);
    println!("Sampled {} points from the curve", samples.len());

    // Store errors for plotting
    let mut iterations = Vec::new();
    let mut nearest_errors = Vec::new();
    let mut gauss_errors = Vec::new();
    let mut new_errors = Vec::new();

    // Set random seed for reproducibility
    let _rng = StdRng::seed_from_u64(42);
    // Run all methods for multiple iterations
    for iteration in 1..=20 {
        // Run nearest point method
        let nearest_point_result =
            fit_cubic_bezier_alternating(&samples, iteration, 1e-6, TUpdateMethod::NearestPoint)?;

        // Run Gauss-Newton method
        let gauss_newton_result =
            fit_cubic_bezier_alternating(&samples, iteration, 1e-6, TUpdateMethod::GaussNewton)?;

        // Run new method
        let new_result = fit_cubic_bezier_weak_varpro(&samples, iteration, 1e-6, None)?;

        // Calculate total error for nearest point method
        let nearest_error: f64 = samples
            .iter()
            .map(|point| {
                let (nearest_point, _) = nearest_point_result.nearest_point(point);
                point.distance(&nearest_point)
            })
            .sum();

        // Calculate total error for Gauss-Newton method
        let gauss_error: f64 = samples
            .iter()
            .map(|point| {
                let (nearest_point, _) = gauss_newton_result.nearest_point(point);
                point.distance(&nearest_point)
            })
            .sum();

        // Calculate total error for new method
        let new_error: f64 = samples
            .iter()
            .map(|point| {
                let (nearest_point, _) = new_result.nearest_point(point);
                point.distance(&nearest_point)
            })
            .sum();

        // Store data for plotting
        iterations.push(iteration as f64);
        nearest_errors.push(nearest_error);
        gauss_errors.push(gauss_error);
        new_errors.push(new_error);

        // Output in CSV format for easy plotting
        println!(
            "{},{:.6},{:.6},{:.6}",
            iteration, nearest_error, gauss_error, new_error
        );
    }

    // Create the plot
    let plot_path = output_dir.join("convergence_plot.png");
    let root = BitMapBackend::new(&plot_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Convergence Comparison", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..20f64, 0f64..nearest_errors[0])?;

    chart
        .configure_mesh()
        .x_desc("Iterations")
        .y_desc("Error")
        .draw()?;

    // Plot nearest point method
    chart
        .draw_series(LineSeries::new(
            iterations
                .iter()
                .zip(nearest_errors.iter())
                .map(|(&x, &y)| (x, y)),
            &RED,
        ))?
        .label("Pastva Variant 1 Nearest Point")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

    // Plot Gauss-Newton method
    chart
        .draw_series(LineSeries::new(
            iterations
                .iter()
                .zip(gauss_errors.iter())
                .map(|(&x, &y)| (x, y)),
            &BLUE,
        ))?
        .label("Pastva Variant 2 Gauss-Newton")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    // Plot new method
    chart
        .draw_series(LineSeries::new(
            iterations
                .iter()
                .zip(new_errors.iter())
                .map(|(&x, &y)| (x, y)),
            &BLACK,
        ))?
        .label("Weak Variable Projection")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    println!("Plot saved as {}", plot_path.display());

    Ok(())
}
