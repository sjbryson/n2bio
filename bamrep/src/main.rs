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
#[command(author, version, about = "BAM Alignment Stats")]
struct Args {
    /// Input name-sorted BAM file
    #[arg(short = 'b', long)]
    bam: PathBuf,

    /// Output JSON report file
    #[arg(short = 'r', long)]
    report: PathBuf,

    /// Generate html plots
    #[arg(long)]
    html: bool,

    /// Generate svg plots
    #[arg(long)]
    plot: bool,

    /// Minimum MAPQ score for insert size calculation
    #[arg(short = 'q', long, default_value_t = 40)]
    min_mapq: usize,

    /// Max insert size to use for summary stats calculation
    #[arg(short = 'i', long, default_value_t = 1000)]
    max_ins: usize,

    /// Max read length to use
    #[arg(short = 'l', long, default_value_t = 150)]
    max_len: usize,
}

// ============================================================================
// Stats
// ============================================================================

#[derive(Serialize, Default)]
struct HistogramData {
    bin_min: f64,
    bin_size: f64,
    counts: Vec<u32>,
}

#[derive(Serialize, Default)]
struct StatSummary {
    count: f64,
    mean: f64,
    median: f64,
    stdev: f64,
    min: f64,
    max: f64,
    histogram: HistogramData,
}

impl StatSummary {
    fn calculate(name: &str, data: &mut [f64], max_ins: f64, max_len: f64) -> Self {
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

        // population vaiance
        let variance: f64 = data.iter().map(|&value| {
            let diff: f64 = mean - value;
            diff * diff
        }).sum::<f64>() / count;

        // sample variance
        //let variance_denom = if count > 1.0 { count - 1.0 } else { 1.0 };
        //let variance: f64 = data.iter().map(|&value| {
        //    let diff: f64 = mean - value;
        //    diff * diff
        //}).sum::<f64>() / variance_denom;
        
        let stdev: f64 = variance.sqrt();

        let config: ReportConfig = ReportConfig::from_stat(name, min, max, max_ins, max_len);
        let num_bins: usize = ((config.max - config.min) / config.bin_size).ceil() as usize;
        let mut bins: Vec<u32> = vec![0u32; num_bins.max(1)];

        for &val in data.iter() {
            if val < config.min || val >= config.max { continue; }
            let mut idx: usize = ((val - config.min) / config.bin_size) as usize;
            if idx >= bins.len() { idx = bins.len() - 1; }
            bins[idx] += 1;
        }

        StatSummary {
            count, mean, median, stdev, min, max,
            histogram: HistogramData {
                bin_min: config.min,
                bin_size: config.bin_size,
                counts: bins,
            }
        }
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
// Report configuration
// ============================================================================

struct ReportConfig {
    min: f64,
    max: f64,
    bin_size: f64,
}

impl ReportConfig {
    /// Takes the calculated min and max to determine the binning rules
    fn from_stat(stat_name: &str, mut min: f64, mut max: f64, max_insert: f64, max_len: f64) -> Self {
        let base_name = if stat_name.starts_with("r1_") || stat_name.starts_with("r2_") {
            &stat_name[3..]
        } else {
            stat_name
        };

        let max_as: f64 = 2.0 * max_len + 1.0;
        let max_al: f64 = max_len + 1.0;
       
        match base_name {
            // MAPQ is canonically 0 to 60 (sometimes up to 255, but usually 60). Bin by 1.
            "mapq" => ReportConfig { min: 0.0, max: 61.0, bin_size: 1.0 },
            // Align score is max 300 for 150 base reads **change with arg if needed**
            "align_score" => ReportConfig { min: 0.0, max: max_as, bin_size: 1.0 },
            // Align len is max 150 for 150 base reads **change with arg if needed**
            "align_length" => ReportConfig { min: 0.0, max: max_al, bin_size: 1.0 },
            // Align score/Align len is max 2 for 150 base reads **change with arg if needed**
            "as_al" => ReportConfig { min: 0.0, max: 2.01, bin_size: 0.01 },
            // Proportion is exactly 0.0 to 1.0. Use 50 bins of 0.02.
            "align_proportion" => ReportConfig { min: 0.0, max: 1.01, bin_size: 0.01 },
            // Accuracy is a percentage 0.0 to 100.0. Bin by 1.0.
            "align_accuracy" => ReportConfig { min: 0.0, max: 101.0, bin_size: 1.0 },
            // pe_insert_size range 0 to args.max_ins
            "pe_insert_size" => ReportConfig { min: 0.0, max: max_insert, bin_size: 10.0 },
            // Default dynamic scaling
            _ => {
                if (max - min).abs() < f64::EPSILON {
                    min -= 1.0;
                    max += 1.0;
                }
                let padding: f64 = (max - min) * 0.05;
                max += padding;
                // Target ~100 bins for dynamic stats
                let bin_size: f64 = (max - min) / 100.0;
                ReportConfig { min, max, bin_size }
            }
        }
    }
}

// ============================================================================
// HTML output
// ============================================================================

fn generate_html_report(results: &HashMap<String, StatSummary>, report_path: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let mut html_path = report_path.clone();
    html_path.set_extension("html");

    let json_data: String = serde_json::to_string(results)?;

    // Helper to format the top single-column table (Insert Size)
    let get_single_table = |name: &str| {
        let s: &StatSummary = results.get(name).unwrap();
        format!(
            r#"<table>
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Count</td><td>{:.0}</td></tr>
                <tr><td>Mean</td><td>{:.2}</td></tr>
                <tr><td>Median</td><td>{:.2}</td></tr>
                <tr><td>StdDev</td><td>{:.2}</td></tr>
                <tr><td>Min</td><td>{:.2}</td></tr>
                <tr><td>Max</td><td>{:.2}</td></tr>
            </table>"#,
            s.count, s.mean, s.median, s.stdev, s.min, s.max
        )
    };

    // Helper to format the combined R1/R2 tables
    let get_combined_table = |base_name: &str, r1_name: &str, r2_name: &str| {
        let r1: &StatSummary = results.get(r1_name).unwrap();
        let r2: &StatSummary = results.get(r2_name).unwrap();
        format!(
            r#"<table>
                <tr><th>Metric</th><th>R1</th><th>R2</th></tr>
                <tr><td>Count</td><td>{:.0}</td><td>{:.0}</td></tr>
                <tr><td>Mean</td><td>{:.2}</td><td>{:.2}</td></tr>
                <tr><td>Median</td><td>{:.2}</td><td>{:.2}</td></tr>
                <tr><td>StdDev</td><td>{:.2}</td><td>{:.2}</td></tr>
                <tr><td>Min</td><td>{:.2}</td><td>{:.2}</td></tr>
                <tr><td>Max</td><td>{:.2}</td><td>{:.2}</td></tr>
                <tr style="background-color: #fff3e0; font-weight: bold;">
                    <td>Count &ge; <span id="{}_thresh_val">0</span></td>
                    <td id="{}_thresh_count">-</td>
                    <td id="{}_thresh_count">-</td>
                </tr>
            </table>"#,
            r1.count, r2.count, r1.mean, r2.mean, r1.median, r2.median, r1.stdev, r2.stdev, r1.min, r2.min, r1.max, r2.max,
            base_name, r1_name, r2_name
        )
    };

    let html_content: String = format!(r#"
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>BAM Alignment Report</title>
    <script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
    <style>
        body {{ font-family: sans-serif; background-color: #f9f9f9; margin: 0; padding: 20px; }}
        .report-container {{ max-width: 1400px; margin: 0 auto; background: #fff; padding: 30px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
        
        h1 {{ text-align: center; margin-bottom: 40px; }}
        
        /* Centered horizontal divider with text */
        .divider {{ display: flex; align-items: center; text-align: center; margin: 50px 0 20px 0; font-size: 1.2em; font-weight: bold; color: #444; }}
        .divider::before, .divider::after {{ content: ''; flex: 1; border-bottom: 2px solid #eee; }}
        .divider:not(:empty)::before {{ margin-right: 15px; }}
        .divider:not(:empty)::after {{ margin-left: 15px; }}

        /* CSS Grid Layouts */
        .grid-top {{ display: grid; grid-template-columns: 80% 20%; gap: 10px; align-items: center; margin-bottom: 20px; }}
        /* 2fr 1fr 2fr creates a 40% 20% 40% distribution */
        .grid-row {{ display: grid; grid-template-columns: 2fr 1fr 2fr; gap: 15px; align-items: center; margin-bottom: 20px; }}
        
        /* Table Styling */
        table {{ width: 100%; border-collapse: collapse; font-size: 13px; text-align: right; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; }}
        th {{ background-color: #f4f4f4; text-align: center; font-weight: bold; }}
        tr:nth-child(even) {{ background-color: #fafafa; }}
        td:first-child {{ text-align: left; font-weight: bold; }}

        /* Text Styling */
        p.section-desc {{ 
            text-align: center; 
            color: #666; 
            font-size: 0.95em; 
            margin: -10px auto 20px auto;
            max-width: 800px; 
            line-height: 1.5;
        }}

        /* Slider styles */
        .slider-row {{ display: flex; justify-content: center; align-items: center; margin-bottom: 40px; gap: 15px; font-size: 14px; font-weight: bold; color: #444; }}
        input[type=range] {{ width: 300px; }}
    </style>
</head>
<body>
    <div class="report-container">
        <h1>BAM Alignment Report</h1>

        <div class="divider">PE Insert Sizes</div>
        <p class="section-desc">
            Distribution of paired-end insert sizes for uniquely mapped, high-quality pairs. 
            Inserts are only calculated for concordant alignments where both had a MAPQ 
            value above the "--min_mapq" (default=40) threshold. Reads exceeding the 
            max insert size threshold, "--max_ins" (default=1000), are excluded.
        </p>
        <div class="grid-top">
            <div id="pe_insert_size_plot"></div>
            {table_insert}
        </div>

        <div class="divider">MAPQ Distribution</div>
        <p class="section-desc">
            Mapping quality (MAPQ) scores representing the aligner's confidence in the read's origin. 
            Higher scores indicate greater probability of correct placement. Unaligned reads are not included.
        </p>
        <div class="grid-row">
            <div id="r1_mapq"></div> {table_mapq} <div id="r2_mapq"></div>
        </div>
        <div class="slider-row">
            <label>Threshold:</label>
            <input type="range" id="mapq_slider">
        </div>

        <div class="divider">Alignment Scores (AS)</div>
        <p class="section-desc">
            Raw alignment scores indicating how well each read matches the reference genome, 
            accounting for matches, mismatches, and gaps. This is the value of the "AS" tag 
            in the sam/bam record. Unaligned reads are not included. Only scores >= 0 are plotted.
        </p>
        <div class="grid-row">
            <div id="r1_align_score"></div> {table_as} <div id="r2_align_score"></div>
        </div>
        <div class="slider-row">
            <label>Threshold:</label>
            <input type="range" id="align_score_slider">
        </div>

        <div class="divider">Alignment Lengths (AL)</div>
        <p class="section-desc">
            Alignment lengths are calculated from the CIGAR string. 
            Matches, mismatches, and indels are counted; clipped regions are not. 
            Unaligned reads are not included. AL is capped at the max read length.
        </p>
        <div class="grid-row">
            <div id="r1_align_length"></div> {table_al} <div id="r2_align_length"></div>
        </div>
        <div class="slider-row">
            <label>Threshold:</label>
            <input type="range" id="align_length_slider">
        </div>

        <div class="divider">AS per Base</div>
        <p class="section-desc">
            This is the record's Alignment Score divided by the Alignment Length (AS/AL).
            Unaligned reads are not included. Only scores >= 0 are plotted.
        </p>
        <div class="grid-row">
            <div id="r1_as_al"></div> {table_asal} <div id="r2_as_al"></div>
        </div>
        <div class="slider-row">
            <label>Threshold:</label>
            <input type="range" id="as_al_slider">
        </div>

        <div class="divider">Alignment Proportions (AP)</div>
        <p class="section-desc">
            This is the record's Alignment Length divided by the Read Length (AL/RL).
            Unaligned reads are not included.
        </p>
        <div class="grid-row">
            <div id="r1_align_proportion"></div> {table_ap} <div id="r2_align_proportion"></div>
        </div>
        <div class="slider-row">
            <label>Threshold:</label>
            <input type="range" id="align_proportion_slider">
        </div>

        <div class="divider">Alignment Percent Identity (PI)</div>
        <p class="section-desc">
            This is the record's number of matches (sam/bam tag "NM") divided by the Alignment Length (100 * NM/AL).
            Unaligned reads are not included.
        </p>
        <div class="grid-row">
            <div id="r1_align_accuracy"></div> {table_acc} <div id="r2_align_accuracy"></div>
        </div>
        <div class="slider-row">
            <label>Threshold:</label>
            <input type="range" id="align_accuracy_slider">
        </div>
    </div>

    <script>
        const data = {json_data};
        
        function draw(id, data_key, plot_title) {{
            const s = data[data_key];
            const trace = {{
                x: s.histogram.counts.map((_, i) => s.histogram.bin_min + (i * s.histogram.bin_size)),
                y: s.histogram.counts, 
                type: 'bar', 
                name: plot_title,
                marker: {{ 
                    color: 'rgb(218, 235, 254)', 
                    line: {{ color: 'rgb(160, 184, 206)', width: 1 }} // Added a slight border for crispness
                }}
            }};
            Plotly.newPlot(id, [trace], {{ 
                title: plot_title, 
                margin: {{t:40, b:40, l:50, r:20}} 
            }});
        }}

        // Draw Insert Size
        draw('pe_insert_size_plot', 'pe_insert_size', 'PE Insert Sizes');

        // Draw MAPQ
        draw('r1_mapq', 'r1_mapq', 'R1 MAPQ Distribution');
        draw('r2_mapq', 'r2_mapq', 'R2 MAPQ Distribution');

        // Draw Align Score
        draw('r1_align_score', 'r1_align_score', 'R1 Alignment Scores (AS)');
        draw('r2_align_score', 'r2_align_score', 'R2 Alignment Scores (AS)');

        // Draw Align Length
        draw('r1_align_length', 'r1_align_length', 'R1 Alignment Lengths (AL)');
        draw('r2_align_length', 'r2_align_length', 'R2 Alignment Lengths (AL)');

        // Draw AS per Base
        draw('r1_as_al', 'r1_as_al', 'R1 AS per Base');
        draw('r2_as_al', 'r2_as_al', 'R2 AS per Base');

        // Draw Align Proportion
        draw('r1_align_proportion', 'r1_align_proportion', 'R1 Alignment Proportions (AP)');
        draw('r2_align_proportion', 'r2_align_proportion', 'R2 Alignment Proportions (AP)');

        // Draw Align Accuracy (PI)
        draw('r1_align_accuracy', 'r1_align_accuracy', 'R1 Alignment Percent Identity (PI)');
        draw('r2_align_accuracy', 'r2_align_accuracy', 'R2 Alignment Percent Identity (PI)');

        // Interactive Threshold Logic
        function setupSlider(base_name, r1_name, r2_name) {{
            const slider = document.getElementById(base_name + '_slider');
            if (!slider) return;

            const r1_data = data[r1_name];
            const r2_data = data[r2_name];
            
            const min = r1_data.histogram.bin_min;
            const bin_size = r1_data.histogram.bin_size;
            const max = min + (r1_data.histogram.counts.length * bin_size);

            // Configure slider dynamically based on the bin configuration
            slider.min = min;
            slider.max = max;
            slider.step = bin_size;
            slider.value = min;

            slider.addEventListener('input', function(e) {{
                const thresh = parseFloat(e.target.value);
                
                // 1. Update the table labels
                document.getElementById(base_name + '_thresh_val').innerText = thresh.toFixed(2);

                // 2. Calculate R1 Count >= Threshold
                let r1_sum = 0;
                r1_data.histogram.counts.forEach((count, i) => {{
                    if (min + (i * bin_size) >= thresh) r1_sum += count;
                }});
                document.getElementById(r1_name + '_thresh_count').innerText = r1_sum;

                // 3. Calculate R2 Count >= Threshold
                let r2_sum = 0;
                r2_data.histogram.counts.forEach((count, i) => {{
                    if (min + (i * bin_size) >= thresh) r2_sum += count;
                }});
                document.getElementById(r2_name + '_thresh_count').innerText = r2_sum;

                // 4. Draw a vertical line on both plots to show the threshold
                const line_update = {{
                    shapes: [{{
                        type: 'line',
                        x0: thresh, x1: thresh,
                        y0: 0, y1: 1, yref: 'paper',
                        line: {{ color: 'red', width: 2, dash: 'dot' }}
                    }}]
                }};
                Plotly.relayout(r1_name, line_update);
                Plotly.relayout(r2_name, line_update);
            }});

            // Trigger once to initialize table values
            slider.dispatchEvent(new Event('input'));
        }}

        // Setup interactive sliders
        setupSlider('mapq', 'r1_mapq', 'r2_mapq');
        setupSlider('align_score', 'r1_align_score', 'r2_align_score');
        setupSlider('align_length', 'r1_align_length', 'r2_align_length');
        setupSlider('as_al', 'r1_as_al', 'r2_as_al');
        setupSlider('align_proportion', 'r1_align_proportion', 'r2_align_proportion');
        setupSlider('align_accuracy', 'r1_align_accuracy', 'r2_align_accuracy');
    </script>
</body>
</html>
"#, 
    json_data = json_data,
    table_insert = get_single_table("pe_insert_size"),
    table_mapq = get_combined_table("mapq","r1_mapq", "r2_mapq"),
    table_as = get_combined_table("align_score","r1_align_score", "r2_align_score"),
    table_al = get_combined_table("align_length","r1_align_length", "r2_align_length"),
    table_asal = get_combined_table("as_al","r1_as_al", "r2_as_al"),
    table_ap = get_combined_table("align_proportion","r1_align_proportion", "r2_align_proportion"),
    table_acc = get_combined_table("align_accuracy","r1_align_accuracy", "r2_align_accuracy")
    );

    std::fs::write(html_path, html_content)?;
    Ok(())
}

// ============================================================================
// Plot histogram
// ============================================================================

fn plot_histogram(
    stat_name: &str, 
    summary: &StatSummary,
    report_path: &PathBuf
) -> Result<(), Box<dyn std::error::Error>> {
    let hist: &HistogramData = &summary.histogram;
    if hist.counts.is_empty() { return Ok(()); }

    let ucla_lightest_blue: RGBColor = RGBColor(218, 235, 254);
    let mut plot_path: PathBuf = report_path.clone();
    plot_path.set_extension(format!("{}.svg", stat_name));

    let root: DrawingArea<SVGBackend<'_>, plotters::coord::Shift> = SVGBackend::new(&plot_path, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_freq: u32 = *hist.counts.iter().max().unwrap_or(&1);
    let y_max: u32 = (max_freq as f64 * 1.05).ceil() as u32;
    let x_max: f64 = hist.bin_min + (hist.counts.len() as f64 * hist.bin_size);

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Histogram of {}", stat_name), ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(60)
        .y_label_area_size(60)
        .build_cartesian_2d(hist.bin_min..x_max, 0..y_max)?;

    chart.configure_mesh()
        .x_labels(hist.counts.len())
        .x_label_formatter(&|v| format!("{:.2}  ", v))
        .x_label_style(("sans-serif", 12).into_font().transform(FontTransform::Rotate270))
        .y_desc("Count")
        .draw()?;

    chart.draw_series(
        hist.counts.iter().enumerate().map(|(i, &count)| {
            let bin_left: f64 = hist.bin_min + (i as f64 * hist.bin_size);
            let bin_right: f64 = bin_left + hist.bin_size; 
            
            Rectangle::new(
                [(bin_left, 0), (bin_right, count)], 
                ucla_lightest_blue.filled()
            )
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

    // Set some max parapeter values based on args.max_len
    let max_mapq: f64 = 60.0;
    let max_as: f64 = 2.0 * args.max_len as f64;
    let max_al: f64 = args.max_len as f64;
    let max_as_al: f64 = 2.0;
    let max_align_prop: f64 = 1.0;
    let max_align_pi: f64 = 100.0;

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
        // If read is unmapped push 0 to each vec
        //if r.mapq == 0 {
        //    mapq_vec.push(0 as f64);
        //    align_score_vec.push(0 as f64);
        //    align_len_vec.push(0 as f64);
        //    as_al_vec.push(0 as f64);
        //    align_prop_vec.push(0 as f64);
        //    align_acc_vec.push(0 as f64);
        //}

        // Only push alignment metrics if the read is mapped
        if r.mapq > 0 {
            mapq_vec.push((r.mapq as f64).min(max_mapq));
            
            if let Some(val) = r.get_int_tag(b"AS") { 
                align_score_vec.push((val as f64).min(max_as)); 
            }
            if let Some(val) = r.calculate_alignment_length() { 
                align_len_vec.push((val as f64).min(max_al)); 
            }
            if let Some(val) = r.calculate_as_al() { 
                as_al_vec.push((val as f64).min(max_as_al)); 
            }
            if let Some(val) = r.calculate_alignment_proportion() { 
                align_prop_vec.push((val as f64).min(max_align_prop)); 
            }
            if let Some(val) = r.calculate_alignment_accuracy() { 
                align_acc_vec.push((val as f64).min(max_align_pi)); 
            }
        }
    };

    // Closure to process the full pair
    let mut process_pair = |r1: &BamRecord, r2: &BamRecord| {
        extract_read_stats(r1, true);
        extract_read_stats(r2, false);

        // Mapq filtering and insert size logic
        // Only calculate insert size for highly confident, uniquely mapped pairs
        if r1.mapq as usize >= args.min_mapq && r2.mapq as usize >= args.min_mapq {
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
        eprintln!("No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?");
        std::process::exit(1);
    }

   println!("BAM reading complete. Processed {} pairs. Generating summaries and plots...", total_pairs);

    // Group all 11 vectors into a list so Rayon can process them concurrently
    let mut stats_to_process: Vec<(&str, &mut Vec<f64>)> = vec![
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

    let max_ins_val: f64 = args.max_ins as f64;
    let max_len_val: f64 = args.max_len as f64;
    // Process all 11 stats concurrently using Rayon
    let results_vec: Vec<(String, StatSummary)> = stats_to_process
        .par_iter_mut()
        .map(|(name, data)| {
            // Pass the name to get the correct config
            let summary: StatSummary = StatSummary::calculate(name, data, max_ins_val, max_len_val);
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

    // Handle Optional Output Formats
    if args.plot {
        for (name, summary) in &results {
            if let Err(e) = plot_histogram(name, summary, &args.report) {
                eprintln!("Warning: Failed to plot {}: {}", name, e);
            }
        }
    }

    if args.html {
        if let Err(e) = generate_html_report(&results, &args.report) {
            eprintln!("Warning: Failed to generate HTML report: {}", e);
        }
    }

    println!("Analysis complete.");
    println!("Processed {} pairs.", total_pairs);
    println!("Results saved to {:?}", args.report);
    
    Ok(())
}

