//! n2core/bamrep/src/main.rs
//! 

use clap::Parser;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
use plotters::prelude::*;
use rayon::prelude::*;

use n2core::bam::{ BamReader, BamHeader, BamRecord, BamStats };

// ============================================================================
// Args
// ============================================================================

#[derive(Parser, Debug)]
#[command(author, version, about = "BAM Alignment Stats Extractor")]
struct Args {
    /// Input name-sorted BAM file
    #[arg(short = 'b', long)]
    bam: PathBuf,

    /// Output JSON report file
    #[arg(short = 'r', long)]
    report: PathBuf,

    /// Generate histogram plots
    #[arg(short = 'p', long)]
    plot: bool,

    /// Minimum MAPQ score for insert size calculation
    #[arg(short = 'q', long, default_value_t = 40)]
    mapq: usize,

    /// Max insert size calculation
    #[arg(short = 'i', long, default_value_t = 1000)]
    max_ins: usize,
}

// ============================================================================
// Stats
// ============================================================================

/// JSON output format for a single statistic
#[derive(Serialize, Default, Debug)]
struct StatSummary {
    count: f64,
    mean: f64,
    median: f64,
    stdev: f64,
    min: f64,
    max: f64,
}

impl StatSummary {
    fn calculate(data: &mut [f64]) -> Self {
        if data.is_empty() {
            return Self::default();
        }
        // Sort data for median, min, and max calculations
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let count: f64 = data.len() as f64;
        let sum:   f64 = data.iter().sum();
        let mean:  f64 = sum / count;
        let min:   f64 = data[0];
        let max:   f64 = data[data.len() - 1];
        
        let median: f64 = if data.len() % 2 == 0 {
            let mid: usize = data.len() / 2;
            (data[mid - 1] + data[mid]) / 2.0
        } else {
            data[data.len() / 2]
        };

        let variance: f64 = data.iter().map(|&value| {
            let diff: f64 = mean - value;
            diff * diff
        }).sum::<f64>() / count;
        
        let stdev: f64 = variance.sqrt();

        Self { count, mean, median, stdev, min, max }
    }
}

/// Holds raw extracted values to compute summaries and plot histograms later
#[derive(Default)]
struct StatsAccumulator {
    pe_insert_size: Vec<f64>,

    r1_mapq: Vec<f64>,
    r1_align_score: Vec<f64>,
    r1_align_length: Vec<f64>,
    r1_as_al: Vec<f64>,
    r1_align_proportion: Vec<f64>,
    r1_align_accuracy: Vec<f64>,

    r2_mapq: Vec<f64>,
    r2_align_score: Vec<f64>,
    r2_align_length: Vec<f64>,
    r2_as_al: Vec<f64>,
    r2_align_proportion: Vec<f64>,
    r2_align_accuracy: Vec<f64>,
}

// ============================================================================
// Plot - histogram
// ============================================================================

struct PlotConfig {
    min: f64,
    max: f64,
    bin_size: f64,
}

impl PlotConfig {
    /// Determines the best plotting bounds based on the stat name and summary data
    fn from_stat(stat_name: &str, summary: &StatSummary) -> Self {
        // Strip the read identifier prefix so R1 and R2 use the same scale
        let base_name = if stat_name.starts_with("r1_") || stat_name.starts_with("r2_") {
            &stat_name[3..]
        } else {
            stat_name
        };

        match base_name {
            // MAPQ is canonically 0 to 60 (sometimes up to 255, but usually 60). Bin by 1.
            "mapq" => PlotConfig { min: 0.0, max: 62.0, bin_size: 2.0 },
            // Align score is max 300 for 150 base reads **change with arg if needed**
            "align_score" => PlotConfig { min: 0.0, max: 306.0, bin_size: 6.0 },
            // Align len is max 150 for 150 base reads **change with arg if needed**
            "align_length" => PlotConfig { min: 0.0, max: 152.0, bin_size: 2.0 },
            // Align score/Align len is max 2 for 150 base reads **change with arg if needed**
            "as_al" => PlotConfig { min: 0.0, max: 2.1, bin_size: 0.05 },
            // Proportion is exactly 0.0 to 1.0. Use 50 bins of 0.02.
            "align_proportion" => PlotConfig { min: 0.0, max: 1.05, bin_size: 0.05 },
            // Accuracy is a percentage 0.0 to 100.0. Bin by 1.0.
            "align_accuracy" => PlotConfig { min: 0.0, max: 102.0, bin_size: 2.0 },
            // pe_insert_size range 0 to args.max_ins
            //"pe_insert_size" => PlotConfig { min: 0.0, max: max_insert, bin_size: 10.0 },
            // Default dynamic scaling for insert size, align score, as_al, and align length
            _ => {
                let mut min: f64 = summary.min;
                let mut max: f64 = summary.max;
                
                // Fix the "missing bars" issue if all reads have the exact same value
                if (max - min).abs() < f64::EPSILON {
                    min -= 1.0;
                    max += 1.0;
                }

                // Add a 5% buffer to the max so the highest bar isn't cut off by the right axis
                let padding: f64 = (max - min) * 0.05;
                max += padding;

                // Calculate a dynamic bin size targeting ~50 bins
                let bin_size: f64 = (max - min) / 20.0;
                
                PlotConfig { min, max, bin_size }
            }
        }
    }
}

fn plot_histogram(
    data: &[f64], 
    stat_name: &str, 
    summary: &StatSummary,
    report_path: &PathBuf
) -> Result<(), Box<dyn std::error::Error>> {
    if data.is_empty() { return Ok(()); }

    let config: PlotConfig = PlotConfig::from_stat(stat_name, summary);

    // Define a custom light blue-grey color
    let light_blue_grey: RGBColor = RGBColor(160, 184, 206);
    
    // Determine plot filename
    let mut plot_path: PathBuf = report_path.clone();
    plot_path.set_extension(format!("{}.svg", stat_name));

    let root = SVGBackend::new(&plot_path, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    // Calculate how many bins we actually need based on the config
    let num_bins: usize = ((config.max - config.min) / config.bin_size).ceil() as usize;
    let mut bins: Vec<u32> = vec![0u32; num_bins.max(1)]; // Ensure at least 1 bin

    // Safely assign data to bins
    for &val in data {
        // Ignore severe outliers that fall completely outside our expected hardcoded ranges
        if val < config.min || val >= config.max { continue; }
        
        let mut idx: usize = ((val - config.min) / config.bin_size) as usize;
        if idx >= bins.len() { idx = bins.len() - 1; }
        bins[idx] += 1;
    }

    let max_freq: u32 = *bins.iter().max().unwrap_or(&1);
    // Add 1% padding to the Y-axis so the tallest bar doesn't touch the top
    let y_max: u32 = (max_freq as f64 * 1.01) as u32;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Histogram of {}", stat_name), ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(config.min..config.max, 0..y_max)?; // Plain f64 range

    let total_bins: usize = ((config.max - config.min) / config.bin_size).round() as usize;

    chart.configure_mesh()
        .x_labels(total_bins) 
        .x_label_formatter(&|v| format!("{:.1}   ", v))
        .x_label_style(("sans-serif", 14).into_font().transform(FontTransform::Rotate270))
        .y_desc("Count")
        .draw()?;

    // Draw bins manually using Rectangles
    chart.draw_series(
        bins.into_iter().enumerate().map(|(i, count)| {
            let bin_left: f64 = config.min + (i as f64 * config.bin_size);
            let bin_right: f64 = bin_left + config.bin_size; 
            let rect: Rectangle<(f64, u32)> = Rectangle::new(
                [(bin_left, 0), (bin_right, count)], 
                light_blue_grey.filled()
            );
            rect
        })
    )?;

    root.present()?;
    Ok(())
}

// ============================================================================
// Main
// ============================================================================

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Args = Args::parse();

    let mut bam_reader: BamReader = BamReader::open(args.bam.to_str().unwrap())?;
    let _header: BamHeader = bam_reader.read_header()?;
    
    let mut current_record: BamRecord = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8> = Vec::new();
    
    let mut total_pairs: i32 = 0;
    let mut stats: StatsAccumulator = StatsAccumulator::default();

    // Closure to extract alignment stats from individual reads
    let mut extract_read_stats = |r: &BamRecord, is_r1: bool| {
        // Route the data to the correct vectors
        let (
            mapq_vec, 
            align_score_vec,
            align_len_vec, 
            as_al_vec, 
            align_prop_vec, 
            align_acc_vec,
        ) = if is_r1 {
            (
                &mut stats.r1_mapq,
                &mut stats.r1_align_score, 
                &mut stats.r1_align_length, 
                &mut stats.r1_as_al, 
                &mut stats.r1_align_proportion, 
                &mut stats.r1_align_accuracy
            )
        } else {
            (
                &mut stats.r2_mapq,
                &mut stats.r2_align_score, 
                &mut stats.r2_align_length, 
                &mut stats.r2_as_al, 
                &mut stats.r2_align_proportion, 
                &mut stats.r2_align_accuracy
            )
        };  

        // Only push alignment metrics if the read is mapped
        if r.mapq > 0 {
            mapq_vec.push(r.mapq as f64);
            if let Some(val) = r.get_int_tag(b"AS") { align_score_vec.push(val as f64); }
            if let Some(val) = r.calculate_alignment_length() { align_len_vec.push(val as f64); }
            if let Some(val) = r.calculate_as_al() { as_al_vec.push(val as f64); }
            if let Some(val) = r.calculate_alignment_proportion() { align_prop_vec.push(val as f64); }
            if let Some(val) = r.calculate_alignment_accuracy() { align_acc_vec.push(val as f64); }
        }
    };

    // Closure to process the full pair
    let mut process_pair = |r1: &BamRecord, r2: &BamRecord| {
        extract_read_stats(r1, true);
        extract_read_stats(r2, false);

        // Mapq filtering and insert size logic
        // Only calculate insert size for highly confident, uniquely mapped pairs
        if r1.mapq as usize >= args.mapq && r2.mapq as usize >= args.mapq {
            if r1.ref_id == r2.ref_id && r1.ref_id != -1 {
                let (fwd, rev) = if r1.pos <= r2.pos { (r1, r2) } else { (r2, r1) };
                let ref_span: i32 = rev.calculate_ref_span().unwrap_or(0) as i32;
                let insert_size: i32 = (rev.pos + ref_span) - fwd.pos;
                
                if insert_size > 0 && insert_size <= args.max_ins as i32 {
                    stats.pe_insert_size.push(insert_size as f64);
                }
            }
        }
    };

    // Read loop
    while bam_reader.read_record(&mut current_record)? {
        if current_record.read_name != prev_qname {
            if !prev_qname.is_empty() {
                total_pairs += 1;
                if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
                    process_pair(r1, r2);
                }
            }
            r1_record = None;
            r2_record = None;
            prev_qname = current_record.read_name.clone();
        }

        // FLAG 0x40 = first in pair (R1), FLAG 0x80 = second in pair (R2)
        if current_record.flag & 0x40 != 0 {
            r1_record = Some(current_record.clone());
        } else if current_record.flag & 0x80 != 0 {
            r2_record = Some(current_record.clone());
        }
    }

    // Process final record pair
    if !prev_qname.is_empty() {
        total_pairs += 1;
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            process_pair(r1, r2);
        }
    }

    if total_pairs == 0 {
        // Return your custom error here
        eprintln!("No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?");
        std::process::exit(1);
    }

   println!("BAM reading complete. Processed {} pairs. Generating summaries and plots...", total_pairs);

    // Group all 11 vectors into a list so Rayon can process them concurrently
    let mut stats_to_process = vec![
        ("pe_insert_size", &mut stats.pe_insert_size),
        ("r1_mapq", &mut stats.r1_mapq),
        ("r2_mapq", &mut stats.r2_mapq),
        ("r1_align_score", &mut stats.r1_align_score),
        ("r2_align_score", &mut stats.r2_align_score),
        ("r1_align_length", &mut stats.r1_align_length),
        ("r2_align_length", &mut stats.r2_align_length),
        ("r1_as_al", &mut stats.r1_as_al),
        ("r2_as_al", &mut stats.r2_as_al),
        ("r1_align_proportion", &mut stats.r1_align_proportion),
        ("r2_align_proportion", &mut stats.r2_align_proportion),
        ("r1_align_accuracy", &mut stats.r1_align_accuracy),
        ("r2_align_accuracy", &mut stats.r2_align_accuracy),
    ];

    // Process all 11 stats concurrently using Rayon
    let results_vec: Vec<(String, StatSummary)> = stats_to_process
        .par_iter_mut()
        .map(|(name, data)| {
            let summary: StatSummary = StatSummary::calculate(data);

            if args.plot {
                // Pass `&summary` here!
                if let Err(e) = plot_histogram(data, name, &summary, &args.report) {
                    eprintln!("Warning: Failed to plot {}: {}", name, e);
                }
            }
            
            (name.to_string(), summary)
        })
        .collect();

    // Reconstruct the HashMap for JSON serialization
    let mut results: HashMap<String, StatSummary> = HashMap::new();
    for (name, summary) in results_vec {
        results.insert(name, summary);
    }

    // Write JSON Report
    let report_file: File = File::create(&args.report)?;
    serde_json::to_writer_pretty(report_file, &results)?;

    println!("Analysis complete.");
    println!("Processed {} pairs.", total_pairs);
    println!("Results saved to {:?}", args.report);
    
    Ok(())
}

