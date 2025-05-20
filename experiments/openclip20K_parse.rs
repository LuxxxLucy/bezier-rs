// all all files in the "dataset/subset" directory and parse them into bezier curves.

// parse SVG and get all the holes

// get the path data and then call parsing to convert to bezier curve
// assert correctness and show success rate

use bezier_rs::app::svg_simplifier::{Segment, SvgDoc};
use bezier_rs::data::BezierCurve;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use log::{error, info, warn};
use std::fs;
use std::path::Path;
use std::time::{Duration, Instant};

// const DATASET_DIR: &str = "dataset/subset";
const DATASET_DIR: &str = "dataset/OpenClipArts20K";
const PATH_TIMEOUT: Duration = Duration::from_secs(5);
const FILE_TIMEOUT: Duration = Duration::from_secs(30);
const MAX_FILE_SIZE: u64 = 10 * 1024 * 1024; // 10MB

// Helper function to get context around an error position
fn get_error_context(path_data: &str, error_pos: usize, context_size: usize) -> String {
    let start = error_pos.saturating_sub(context_size);
    let end = (error_pos + context_size).min(path_data.len());
    let mut context = String::new();

    if start > 0 {
        context.push_str("...");
    }
    context.push_str(&path_data[start..end]);
    if end < path_data.len() {
        context.push_str("...");
    }

    let indicator_pos = if start > 0 {
        context_size + 3
    } else {
        error_pos - start
    };
    let indicator = format!("\n{}^", " ".repeat(indicator_pos));
    context.push_str(&indicator);

    context
}

fn format_file_size(size: u64) -> String {
    const KB: u64 = 1024;
    const MB: u64 = KB * 1024;

    if size >= MB {
        format!("{:.1}MB", size as f64 / MB as f64)
    } else if size >= KB {
        format!("{:.1}KB", size as f64 / KB as f64)
    } else {
        format!("{}B", size)
    }
}

fn validate_svg_content(content: &str) -> bool {
    // Basic SVG validation
    content.contains("<svg") && content.contains("</svg>")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    info!("Starting SVG parsing process");

    let start_time = Instant::now();
    let dataset_dir = Path::new(DATASET_DIR);

    let svg_files: Vec<_> = fs::read_dir(dataset_dir)?
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.extension().and_then(|ext| ext.to_str()) == Some("svg") {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    let total_files = svg_files.len();
    info!("Found {} SVG files to process", total_files);

    let multi_progress = MultiProgress::new();

    // Create progress bar for files
    let file_pb = multi_progress.add(ProgressBar::new(total_files as u64));
    file_pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} Files [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));

    // Create progress bar for paths
    let path_pb = multi_progress.add(ProgressBar::new(0));
    path_pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.yellow} Paths [{elapsed_precise}] [{bar:40.yellow/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));

    // Create status line for current file
    let status_pb = multi_progress.add(ProgressBar::new(0));
    status_pb.set_style(ProgressStyle::default_bar().template("{msg}").unwrap());

    let mut total_holes = 0;
    let mut processed_holes = 0;
    let mut failed_paths = 0;
    let mut timeout_paths = 0;
    let mut skipped_files = 0;

    // Process files one by one
    for path in svg_files {
        let file_start_time = Instant::now();
        let file_name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        let file_size = fs::metadata(&path)?.len();

        // Check file size
        if file_size > MAX_FILE_SIZE {
            warn!(
                "Skipping {}: File too large ({} > {})",
                file_name,
                format_file_size(file_size),
                format_file_size(MAX_FILE_SIZE)
            );
            skipped_files += 1;
            file_pb.inc(1);
            continue;
        }

        let svg_content = match fs::read_to_string(&path) {
            Ok(content) => content,
            Err(e) => {
                error!("Failed to read file {}: {}", file_name, e);
                skipped_files += 1;
                file_pb.inc(1);
                continue;
            }
        };

        // Validate SVG content
        if !validate_svg_content(&svg_content) {
            warn!("Skipping {}: Invalid SVG content", file_name);
            skipped_files += 1;
            file_pb.inc(1);
            continue;
        }

        // Parse SVG with timeout and panic protection
        let parse_result = std::panic::catch_unwind(|| {
            let mut svg_doc = SvgDoc::new(svg_content);
            svg_doc.parse();
            svg_doc
        });

        let svg_doc = match parse_result {
            Ok(doc) => doc,
            Err(_) => {
                error!("Failed to parse SVG structure in {}", file_name);
                skipped_files += 1;
                file_pb.inc(1);
                continue;
            }
        };

        // Update path progress bar total
        let file_holes = svg_doc
            .segments
            .iter()
            .filter(|s| matches!(s, Segment::Hole(_)))
            .count();
        total_holes += file_holes;
        path_pb.set_length(total_holes as u64);

        let mut current_file_processed = 0;

        // Process holes immediately
        for (hole_idx, segment) in svg_doc.segments.iter().enumerate() {
            if let Segment::Hole(hole) = segment {
                // Check file processing timeout
                if file_start_time.elapsed() > FILE_TIMEOUT {
                    error!(
                        "File processing timeout for {} after {:.2?}",
                        file_name,
                        file_start_time.elapsed()
                    );
                    skipped_files += 1;
                    break;
                }

                let path_data =
                    svg_doc.content[hole.start_idx..hole.start_idx + hole.len].to_string();

                // Update status line
                status_pb.set_message(format!(
                    "Current file: {} ({}) - Paths: {}/{} - Time: {:.1}s",
                    file_name,
                    format_file_size(file_size),
                    current_file_processed + 1,
                    file_holes,
                    file_start_time.elapsed().as_secs_f64()
                ));

                // Set up timeout
                let parse_start = Instant::now();
                match BezierCurve::parse_maybe_multiple(&path_data) {
                    Ok(_) => {
                        processed_holes += 1;
                        current_file_processed += 1;
                        path_pb.inc(1);
                    }
                    Err(e) => {
                        if parse_start.elapsed() > PATH_TIMEOUT {
                            error!(
                                "Timeout processing path {} in {} ({}): took {:.2?}",
                                hole_idx + 1,
                                file_name,
                                path.display(),
                                parse_start.elapsed()
                            );
                            timeout_paths += 1;
                        } else {
                            let error_context = if let Some(_pos) = e.to_string().find("position") {
                                if let Some(pos_str) = e
                                    .to_string()
                                    .split_whitespace()
                                    .find(|s| s.parse::<usize>().is_ok())
                                {
                                    if let Ok(pos) = pos_str.parse::<usize>() {
                                        get_error_context(&path_data, pos, 20)
                                    } else {
                                        path_data
                                    }
                                } else {
                                    path_data
                                }
                            } else {
                                path_data
                            };

                            error!(
                                "Failed to parse path {} in {} ({}): {}\nError context:\n{}",
                                hole_idx + 1,
                                file_name,
                                path.display(),
                                e,
                                error_context
                            );
                        }
                        failed_paths += 1;
                        current_file_processed += 1;
                        path_pb.inc(1);
                    }
                }
            }
        }
        file_pb.inc(1);
    }

    let duration = start_time.elapsed();
    info!("\nProcessing completed in {:.2?}", duration);
    info!("Total files processed: {}", total_files);
    info!("Files skipped: {}", skipped_files);
    info!("Total paths found: {}", total_holes);
    info!("Successfully processed: {}", processed_holes);
    info!("Failed paths: {}", failed_paths);
    info!("Timeout paths: {}", timeout_paths);
    info!(
        "Success rate: {:.1}%",
        (processed_holes as f64 / total_holes as f64) * 100.0
    );

    Ok(())
}
