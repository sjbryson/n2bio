//! n2core/bamrep/src/report.rs
//! 

use std::path::PathBuf;
use crate::stats::{ ReportData, StatSummary, ReadClassStats };

// ============================================================================
// Report configuration
// ============================================================================

pub struct ReportConfig {
    pub min: f64,
    pub max: f64,
    pub bin_size: f64,
}

impl ReportConfig {
    /// Takes the calculated min and max to determine the binning rules
    pub fn from_stat(stat_name: &str, mut min: f64, mut max: f64, max_insert: f64, max_len: f64) -> Self {
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

pub fn generate_html_report(report_data: &ReportData, html_path: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    
    let json_data: String = serde_json::to_string(report_data)?;

    // Helper to format the top Global Stats table
    let r1: &ReadClassStats = &report_data.global_stats.r1;
    let r2: &ReadClassStats = &report_data.global_stats.r2;
    let table_global: String = format!(
        r#"<table>
            <tr><th>Metric</th><th>R1</th><th>R2</th></tr>
            <tr><td>Total Reads</td><td>{}</td><td>{}</td></tr>
            <tr><td>Primary Mapped</td><td>{}</td><td>{}</td></tr>
            <tr><td>Primary MAPQ > 0</td><td>{}</td><td>{}</td></tr>
            <tr><td>Primary Concordant</td><td>{}</td><td>{}</td></tr>
            <tr><td>Primary Discordant</td><td>{}</td><td>{}</td></tr>
            <tr><td>Primary Singletons</td><td>{}</td><td>{}</td></tr>
            
            <tr><td colspan="3" style="background-color: #f8f9fa; font-weight: bold; text-align: center;">Multi-Mapping & Ambiguity</td></tr>
            <tr><td>Secondary/Supplementary Mapped</td><td>{}</td><td>{}</td></tr>
            <tr><td>MAPQ = 0</td><td>{}</td><td>{}</td></tr>
        </table>"#,
        r1.total_reads, r2.total_reads,
        r1.primary_mapped, r2.primary_mapped,
        r1.primary_mapq, r2.primary_mapq,
        r1.concordant_mapped, r2.concordant_mapped,
        r1.discordant_mapped, r2.discordant_mapped,
        r1.singletons, r2.singletons,
        r1.secondary_mapped, r2.secondary_mapped,
        r1.mapq_0, r2.mapq_0
    );


    // Helper to format the insert size table
    let get_single_table = |name: &str| {
        let s: &StatSummary = report_data.alignment_stats.get(name).unwrap();
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
        let r1: &StatSummary = report_data.alignment_stats.get(r1_name).unwrap();
        let r2: &StatSummary = report_data.alignment_stats.get(r2_name).unwrap();
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
        .report-container {{ max-width: 1200px; margin: 0 auto; background: #fff; padding: 20px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
        
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
        table {{ width: 100%; border-collapse: collapse; font-size: 12px; text-align: right; }}
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
            max-width: 900px; 
            line-height: 1.5;
        }}

        /* Slider styles */
        .slider-row {{ display: flex; justify-content: center; align-items: center; margin-bottom: 40px; gap: 15px; font-size: 12px; font-weight: bold; color: #444; }}
        input[type=range] {{ width: 300px; }}
    </style>
</head>
<body>
    <div class="report-container">
        <h1>BAM Alignment Report</h1>

        <div class="divider">Global Metrics</div>
        <p class="section-desc">Overall sequence counts and mapping classifications for R1 and R2 reads. 
        Only primary alignments with MAPQ>0 are used for alignment stat distributions.</p>
        <div style="max-width: 800px; margin: 0 auto 40px auto;">
            {table_global}
        </div>

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
            Higher scores indicate greater probability of correct placement. Unaligned reads (MAPQ = 0) are not included.
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
            in the sam/bam record. Only scores >= 0 are plotted. Scores are capped at 2 x max read length.
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
            AL is capped at the max read length.
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
            Only scores >= 0 are plotted. The maximum value is 2.
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
        const data = {json_data}.alignment_stats;
        
        function draw(id, data_key, plot_title) {{
            const s = data[data_key];
            const trace = {{
                x: s.histogram.counts.map((_, i) => s.histogram.bin_min + (i * s.histogram.bin_size)),
                y: s.histogram.counts, 
                type: 'bar', 
                name: plot_title,
                marker: {{ 
                    color: 'rgb(218, 235, 254)', 
                    line: {{ color: 'rgb(160, 184, 206)', width: 1 }} 
                }}
            }};
            const layout = {{ 
                title: plot_title, 
                margin: {{t:40, b:40, l:50, r:20}} 
            }};
            const config = {{
                toImageButtonOptions: {{
                    format: 'svg', 
                    filename: plot_title.replace(/\s+/g, '_') 
                }},
                displaylogo: false 
            }};
            Plotly.newPlot(id, [trace], layout, config);
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