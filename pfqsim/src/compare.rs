//! n2bio/pfqsim/src/compare.rs
//! 

use std::io::{self, Error, ErrorKind};
use std::fs::File;
use std::io::Write;
use serde::Serialize;
use std::collections::BTreeMap;
use std::collections::HashMap;

use n2core::hist::Histogram;

use crate::cli::CompareArgs;
use crate::config::CompareConfig;
use crate::analyze_report::AnalyzeReportData;

// ============================================================================
// Compare::run()
// ============================================================================

pub(crate) fn run(args: CompareArgs) -> io::Result<()> {
    // Read config
    println!("Loading comparison configuration: {}", args.config);
    let config: CompareConfig = CompareConfig::from_tsv(&args.config)?;
    
    let mut payloads: Vec<ComparePayload> = Vec::with_capacity(config.rows.len());

    for row in config.rows {
        
        let file: File = File::open(&row.report)?;
        let report: AnalyzeReportData = serde_json::from_reader(file)?; // Uses your struct!

        let mut curves: HashMap<String, CurveData> = std::collections::HashMap::new();
        let p: usize = report.total_expected_positives;
        let n: usize = report.total_expected_negatives;

        // Generate curves for each alignment metric
        curves.insert("mapq".to_string(), compute_curves(&report.true_positives.mapq, &report.false_positives_control.mapq, &report.false_positives_cross.mapq, p, n));
        curves.insert("align_score".to_string(), compute_curves(&report.true_positives.align_score, &report.false_positives_control.align_score, &report.false_positives_cross.align_score, p, n));
        curves.insert("align_length".to_string(), compute_curves(&report.true_positives.align_length, &report.false_positives_control.align_length, &report.false_positives_cross.align_length, p, n));
        curves.insert("as_al".to_string(), compute_curves(&report.true_positives.as_al, &report.false_positives_control.as_al, &report.false_positives_cross.as_al, p, n));
        curves.insert("align_proportion".to_string(), compute_curves(&report.true_positives.align_proportion, &report.false_positives_control.align_proportion, &report.false_positives_cross.align_proportion, p, n));
        curves.insert("percent_identity".to_string(), compute_curves(&report.true_positives.percent_identity, &report.false_positives_control.percent_identity, &report.false_positives_cross.percent_identity, p, n));

        payloads.push(ComparePayload {
            id: row.id,
            curves,
        });
    }

    // Generate output file paths
    let html_path: String = format!("{}.html", args.output);
    let json_path: String = format!("{}.json", args.output);

    // Save the aggregated JSON payload for programmatic use
    let json_file: File = File::create(&json_path)?;
    serde_json::to_writer_pretty(json_file, &payloads)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;
    println!("Aggregated JSON written to: {}", json_path);

    // Generate the HTML Report
    generate_html_report(&payloads, &html_path)?;
    println!("Comparison HTML report written to: {}", html_path);

    Ok(())
}

// ============================================================================
// Curve/AUC Calculation
// ============================================================================

#[derive(Serialize, Debug, Clone)]
pub(crate) struct CurveData {
    pub roc_curve: Vec<[f64; 2]>,
    pub pr_curve: Vec<[f64; 2]>,
    pub roc_auc: f64,
    pub pr_auc: f64,
}

#[derive(Serialize, Debug, Clone)]
struct ComparePayload {
    id: String,
    // Maps a metric name (e.g., "mapq") to its computed curves
    curves: std::collections::HashMap<String, CurveData>,
}

/// Sweeps across the histograms to build ROC and PR curves
fn compute_curves(
    tp_hist: &Histogram,
    fp_ctrl_hist: &Histogram,
    fp_cross_hist: &Histogram,
    total_p: usize,
    total_n: usize,
) -> CurveData {
    // 1. Map values to a unified grid to handle trimmed histograms seamlessly
    let mut grid: BTreeMap<i64, (usize, usize)> = BTreeMap::new();

    let mut insert_hist = |hist: &Histogram, is_tp: bool| {
        for (i, &count) in hist.counts.iter().enumerate() {
            if count > 0 {
                let val: f64 = hist.min_val + (i as f64) * hist.bin_width;
                let grid_idx: i64 = (val / hist.bin_width).round() as i64;
                let entry: &mut (usize, usize) = grid.entry(grid_idx).or_insert((0, 0));
                if is_tp { entry.0 += count; } else { entry.1 += count; }
            }
        }
    };

    insert_hist(tp_hist, true);
    insert_hist(fp_ctrl_hist, false);
    insert_hist(fp_cross_hist, false);

    // FIX: Apply the exact domain logic from the JS implementation
    // Multiply by 2 (e.g., for paired-end reads) and add total cross FPs to the negative pool
    let cross_fp_total: usize = fp_cross_hist.counts.iter().sum();
    let total_p_f64: f64 = (total_p * 2).max(1) as f64;
    let total_n_f64: f64 = ((total_n * 2) + cross_fp_total).max(1) as f64;

    let mut roc_curve: Vec<[f64; 2]> = vec![[0.0, 0.0]];
    let mut pr_curve: Vec<[f64; 2]> = vec![]; 

    let mut tp_sum: usize = 0;
    let mut fp_sum: usize = 0;
    let mut roc_auc: f64 = 0.0;
    let mut pr_auc: f64 = 0.0;

    let mut prev_fpr: f64 = 0.0;
    let mut prev_tpr: f64 = 0.0;
    let mut prev_recall: f64 = 0.0;
    
    // PR precision starts at 1.0 when nothing is selected
    let mut prev_precision: f64 = 1.0; 
    pr_curve.push([0.0, 1.0]);

    // 2. Sweep from highest score downwards
    for (_, (tp_count, fp_count)) in grid.iter().rev() {
        tp_sum += tp_count;
        fp_sum += fp_count;

        let fpr: f64 = fp_sum as f64 / total_n_f64;
        let tpr: f64 = tp_sum as f64 / total_p_f64; // Recall
        let precision: f64 = if (tp_sum + fp_sum) > 0 {
            tp_sum as f64 / (tp_sum + fp_sum) as f64
        } else {
            1.0
        };

        // Trapezoidal integration for AUC
        roc_auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2.0;
        pr_auc += (tpr - prev_recall) * (precision + prev_precision) / 2.0;

        roc_curve.push([fpr, tpr]);
        pr_curve.push([tpr, precision]);

        prev_fpr = fpr;
        prev_tpr = tpr;
        prev_recall = tpr;
        prev_precision = precision;
    }

    // 3. Anchor the curves to negative infinity (unmapped reads)
    // Ensures lines extend to the right and AUC is fully calculated on a [0, 1] domain
    if prev_fpr < 1.0 || prev_tpr < 1.0 {
        let final_fpr: f64 = 1.0;
        let final_tpr: f64 = 1.0;
        let final_precision: f64 = total_p_f64 / (total_p_f64 + total_n_f64);

        roc_auc += (final_fpr - prev_fpr) * (final_tpr + prev_tpr) / 2.0;
        pr_auc += (final_tpr - prev_recall) * (final_precision + prev_precision) / 2.0;

        roc_curve.push([final_fpr, final_tpr]);
        pr_curve.push([final_tpr, final_precision]);
    }

    CurveData { roc_curve, pr_curve, roc_auc, pr_auc }
}

// ============================================================================
// HTML Generation
// ============================================================================

fn generate_html_report(payloads: &[ComparePayload], output_path: &str) -> io::Result<()> {
    let json_data: String = serde_json::to_string(payloads)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let html_content: String = format!(
r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Analysis Comparison Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{
            font-family: system-ui, -apple-system, sans-serif;
            background-color: #f4f4f9;
            color: #333;
            margin: 0;
            padding: 20px;
        }}
        .header-container {{ text-align: center; margin-bottom: 30px; border-bottom: 2px solid #ddd; padding-bottom: 20px; }}
        select {{ padding: 10px; font-size: 16px; border-radius: 5px; margin-top: 10px; cursor: pointer; }}
        
        /* Flexbox layout for Rows */
        .dashboard-row {{
            display: flex;
            justify-content: space-between;
            align-items: stretch;
            background: #fff;
            margin-bottom: 20px;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        
        /* Plot gets 70% width, Table gets 28% */
        .plot-container {{ flex: 0 0 70%; }}
        .table-container {{ flex: 0 0 28%; overflow-x: auto; }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 40px;
        }}
        th, td {{
            text-align: left;
            padding: 10px;
            border-bottom: 1px solid #ddd;
        }}
        th {{ background-color: #f8f9fa; font-weight: 600; }}
        tr:hover {{ background-color: #f1f1f1; }}
    </style>
</head>
<body>
    <div class="header-container">
        <h1>Analysis Comparison Report</h1>
        <label for="metricSelect"><strong>Evaluation Metric: </strong></label>
        <select id="metricSelect">
            <option value="mapq">Mapping Quality (MAPQ)</option>
            <option value="align_score">Alignment Score</option>
            <option value="align_length">Alignment Length</option>
            <option value="as_al">Score per Aligned Base (AS/AL)</option>
            <option value="align_proportion">Alignment Proportion</option>
            <option value="percent_identity">Alignment Accuracy</option>
        </select>
    </div>

    <div class="dashboard-row">
        <div class="plot-container" id="rocPlot"></div>
        <div class="table-container">
            <table>
                <thead><tr><th>ID</th><th>ROC AUC</th></tr></thead>
                <tbody id="rocTableBody"></tbody>
            </table>
        </div>
    </div>

    <div class="dashboard-row">
        <div class="plot-container" id="prPlot"></div>
        <div class="table-container">
            <table>
                <thead><tr><th>ID</th><th>PR AUC</th></tr></thead>
                <tbody id="prTableBody"></tbody>
            </table>
        </div>
    </div>

    <script>
        const data = {json_data}; // Injected from Rust

        function renderDashboard(metricKey) {{
            const rocTraces = [];
            const prTraces = [];
            const rocTable = document.getElementById('rocTableBody');
            const prTable = document.getElementById('prTableBody');
            
            rocTable.innerHTML = '';
            prTable.innerHTML = '';

            data.forEach(dataset => {{
                const curves = dataset.curves[metricKey];
                
                // ROC
                rocTraces.push({{
                    x: curves.roc_curve.map(pt => pt[0]),
                    y: curves.roc_curve.map(pt => pt[1]),
                    mode: 'lines', name: dataset.id, line: {{ width: 2 }}
                }});
                rocTable.innerHTML += `<tr><td>${{dataset.id}}</td><td>${{curves.roc_auc.toFixed(4)}}</td></tr>`;

                // PR
                prTraces.push({{
                    x: curves.pr_curve.map(pt => pt[0]),
                    y: curves.pr_curve.map(pt => pt[1]),
                    mode: 'lines', name: dataset.id, line: {{ width: 2 }}
                }});
                prTable.innerHTML += `<tr><td>${{dataset.id}}</td><td>${{curves.pr_auc.toFixed(4)}}</td></tr>`;
            }});

            const getPlotConfig = (prefix) => ({{
                toImageButtonOptions: {{
                    format: 'svg', 
                    filename: prefix + '_' + metricKey 
                }},
                displaylogo: false 
            }});

            // 2. ROC Plot with its specific config
            Plotly.react('rocPlot', rocTraces, {{
                title: 'Receiver Operating Characteristic (ROC)',
                xaxis: {{ title: 'False Positive Rate (FPR)', range: [-0.02, 1.02] }},
                yaxis: {{ title: 'True Positive Rate (TPR)', range: [-0.02, 1.02] }},
                margin: {{ l: 60, r: 20, t: 40, b: 50 }},
                shapes: [
                    {{
                        type: 'line',
                        x0: 0,
                        y0: 0,
                        x1: 1,
                        y1: 1,
                        line: {{
                            color: 'rgba(0,0,0,0.3)',
                            width: 2,
                            dash: 'dash'
                        }}
                    }}
                ]
            }}, getPlotConfig('ROC_Curve'));

            // 3. PR Plot with its specific config
            Plotly.react('prPlot', prTraces, {{
                title: 'Precision-Recall (PR) Curve',
                xaxis: {{ title: 'Recall', range: [-0.02, 1.02] }},
                yaxis: {{ title: 'Precision', range: [-0.02, 1.02] }},
                margin: {{ l: 60, r: 20, t: 40, b: 50 }}
            }}, getPlotConfig('PR_Curve'));
        }}

        // Listen for dropdown changes
        document.getElementById('metricSelect').addEventListener('change', (e) => {{
            renderDashboard(e.target.value);
        }});

        // Initial render
        renderDashboard('mapq');
    </script>
</body>
</html>"#
    );

    let mut file: File = File::create(output_path)?;
    file.write_all(html_content.as_bytes())?;
    Ok(())
}