//! n2bio/pfqsim/src/report.rs
//! 

use std::fs;
use std::path::PathBuf;
use std::io;
use serde::{Serialize, Deserialize};


// ============================================================================
// Report data
// ============================================================================

/// Serialized data packaged for HTML report
#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct AnalyzeReportData {
    pub(crate) report_name: String,
    pub(crate) total_expected_pairs: usize,
    pub(crate) total_observed_pairs: usize,
    
    // Arrays of score points extracted from alignment pairs
    pub(crate) positives: MetricPayload,
    pub(crate) negatives: MetricPayload,
}

/// Pre-sorted vectors of stats for binary thresholding
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub(crate) struct MetricPayload {
    pub(crate) mapq: Vec<f64>,
    pub(crate) align_score: Vec<f64>,
    pub(crate) align_length: Vec<f64>,
    pub(crate) as_al: Vec<f64>,
    pub(crate) align_proportion: Vec<f64>,
    pub(crate) align_accuracy: Vec<f64>,
}

impl MetricPayload {
    /// Sort vectors in-place so JavaScript can execute binary searches
    pub(crate) fn sort_all(&mut self) {
        self.mapq.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.align_score.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.align_length.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.as_al.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.align_proportion.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.align_accuracy.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }
}

// ============================================================================
// HTML report
// ============================================================================

pub(crate) fn generate_evaluation_report(
    data: &mut AnalyzeReportData, 
    html_path: &PathBuf
) -> io::Result<()> {
    
    // Sort elements for client-side binary search logic
    data.positives.sort_all();
    data.negatives.sort_all();

    let json_raw: String = serde_json::to_string(data).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("JSON serialization error: {}", e))
    })?;

    let html_content: String = format!(r#"
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>pfqsim Benchmark Report: {name}</title>
    <script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background-color: #f4f6f9; margin: 0; padding: 20px; color: #333; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: #fff; padding: 30px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); border-radius: 8px; }}
        
        h1 {{ text-align: center; color: #1e293b; margin-top: 0; margin-bottom: 5px; font-size: 24px; }}
        .subtitle {{ text-align: center; color: #64748b; font-size: 14px; margin-bottom: 30px; }}

        .divider {{ display: flex; align-items: center; text-align: center; margin: 40px 0 20px 0; font-size: 16px; font-weight: bold; color: #475569; }}
        .divider::before, .divider::after {{ content: ''; flex: 1; border-bottom: 2px solid #e2e8f0; }}
        .divider:not(:empty)::before {{ margin-right: 15px; }}
        .divider:not(:empty)::after {{ margin-left: 15px; }}

        /* Layout Grids */
        .dashboard-grid {{ display: grid; grid-template-columns: 350px 1fr 1fr; gap: 25px; margin-bottom: 30px; }}
        .curves-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 30px; }}

        /* Panel Cards */
        .panel {{ background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; }}
        .panel-title {{ font-size: 14px; font-weight: bold; color: #334155; margin-bottom: 15px; text-transform: uppercase; letter-spacing: 0.5px; border-bottom: 1px solid #e2e8f0; padding-bottom: 5px; }}

        /* Interactive Widgets */
        .control-group {{ margin-bottom: 20px; }}
        label {{ display: block; font-size: 12px; font-weight: 600; color: #475569; margin-bottom: 5px; }}
        select, input[type=range] {{ width: 100%; box-sizing: border-box; }}
        .slider-meta {{ display: flex; justify-content: space-between; font-size: 11px; color: #64748b; margin-top: 4px; }}
        .active-val {{ font-weight: bold; color: #2563eb; }}

        /* Matrix Table Styling */
        .matrix-table {{ width: 100%; border-collapse: collapse; text-align: center; font-size: 12px; }}
        .matrix-table th, .matrix-table td {{ padding: 10px; border: 1px solid #cbd5e1; }}
        .matrix-table th {{ background: #f1f5f9; color: #475569; font-weight: 600; }}
        .matrix-label {{ text-align: left; font-weight: bold; background: #f1f5f9; color: #475569; }}
        .cell-tp {{ background: #dcfce7; font-weight: bold; color: #166534; }}
        .cell-fp {{ background: #fee2e2; font-weight: bold; color: #991b1b; }}
        .cell-fn {{ background: #fef9c3; font-weight: bold; color: #854d0e; }}

        /* Summary Stats List */
        .stat-list {{ list-style: none; padding: 0; margin: 0; font-size: 13px; }}
        .stat-list li {{ display: flex; justify-content: space-between; padding: 8px 0; border-bottom: 1px dashed #e2e8f0; }}
        .stat-list li:last-child {{ border: none; }}
        .stat-label {{ color: #475569; font-weight: 500; }}
        .stat-value {{ font-weight: 600; font-family: monospace; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>pfqsim Classifier Assessment Benchmark</h1>
        <div class="subtitle">Dataset Baseline Target: {name}</div>

        <div class="divider">Dynamic Operating Characteristics Threshold Engine</div>

        <div class="dashboard-grid">
            <div class="panel">
                <div class="panel-title">Threshold Controllers</div>
                
                <div class="control-group">
                    <label for="metric_selector">Target Diagnostic Feature</label>
                    <select id="metric_selector">
                        <option value="mapq">Mapping Quality (MAPQ)</option>
                        <option value="align_score">Alignment Score (AS)</option>
                        <option value="align_length">Alignment Length (AL)</option>
                        <option value="as_al">AS / AL (Score per Base)</option>
                        <option value="align_proportion">Alignment Proportion (AL/RL)</option>
                        <option value="align_accuracy">Percent Identity (PI)</option>
                    </select>
                </div>

                <div class="control-group">
                    <label for="threshold_slider">Cutoff Threshold (&ge;)</label>
                    <input type="range" id="threshold_slider">
                    <div class="slider-meta">
                        <span id="slider_min">0</span>
                        <span>Current: <span id="current_threshold_view" class="active-val">0</span></span>
                        <span id="slider_max">100</span>
                    </div>
                </div>

                <div class="panel-title" style="margin-top:25px;">Dynamic Metrics</div>
                <ul class="stat-list">
                    <li><span class="stat-label">Precision (PPV)</span><span id="stat_precision" class="stat-value">-</span></li>
                    <li><span class="stat-label">Recall / TPR (Sensitivity)</span><span id="stat_recall" class="stat-value">-</span></li>
                    <li><span class="stat-label">False Positive Rate (FPR)</span><span id="stat_fpr" class="stat-value">-</span></li>
                    <li><span class="stat-label">F1-Score</span><span id="stat_f1" class="stat-value">-</span></li>
                </ul>
            </div>

            <div class="panel">
                <div class="panel-title">Dynamic Confusion Matrix</div>
                <table class="matrix-table">
                    <tr>
                        <th colspan="2" rowspan="2"></th>
                        <th colspan="2">True Ground Truth Condition</th>
                    </tr>
                    <tr>
                        <th>Positive (Aligned)</th>
                        <th>Negative (Cross-Contam / Misplaced)</th>
                    </tr>
                    <tr>
                        <th rowspan="2" style="writing-mode: vertical-lr; transform: rotate(180deg); font-size:11px;">Decision Outcome</th>
                        <td class="matrix-label">Pass Filter (&ge;)</td>
                        <td id="cell_tp" class="cell-tp">-</td>
                        <td id="cell_fp" class="cell-fp">-</td>
                    </tr>
                    <tr>
                        <td class="matrix-label">Fail Filter (&lt;)</td>
                        <td id="cell_fn" class="cell-fn">-</td>
                        <td style="color:#64748b; font-style:italic; background:#f8fafc;">Excluded</td>
                    </tr>
                </table>
                <p style="font-size:11px; color:#64748b; line-height:1.4; margin-top:15px;">
                    * <strong>False Negatives (FN)</strong> represents simulated read allocations that were either rejected by the filter or dropped entirely from the BAM sequence pipeline during upstream alignment.
                </p>
            </div>

            <div class="panel">
                <div class="panel-title">Simulation Summary Profile</div>
                <ul class="stat-list">
                    <li><span class="stat-label">Expected Simulated Pairs</span><span class="stat-value">{expected}</span></li>
                    <li><span class="stat-label">Observed Primary Pairs</span><span class="stat-value">{observed}</span></li>
                    <li><span class="stat-label">Total TP Alignments (Bases)</span><span class="stat-value">{tp_count}</span></li>
                    <li><span class="stat-label">Total FP Alignments (Bases)</span><span class="stat-value">{fp_count}</span></li>
                </ul>
            </div>
        </div>

        <div class="curves-row">
            <div id="roc_plot" class="panel" style="height: 450px;"></div>
            <div id="pr_plot" class="panel" style="height: 450px;"></div>
        </div>
    </div>

    <script>
        const payload = {json_raw};

        // Binary Search lookup helper: counts values >= threshold inside a sorted array
        fn countGreaterOrEqual(arr, threshold) {{
            let low = 0;
            let high = arr.length;
            while (low < high) {{
                let mid = (low + high) >> 1;
                if (arr[mid] >= threshold) {{
                    high = mid;
                }} else {{
                    low = mid + 1;
                }}
            }}
            return arr.length - low;
        }}

        const metricSelector = document.getElementById('metric_selector');
        const thresholdSlider = document.getElementById('threshold_slider');

        // Dynamic Configuration Rules matching the alignment characteristics
        const metricConfig = {{
            mapq: {{ min: 0, max: 60, step: 1 }},
            align_score: {{ min: 0, max: 300, step: 1 }},
            align_length: {{ min: 0, max: 250, step: 1 }},
            as_al: {{ min: 0, max: 2.0, step: 0.01 }},
            align_proportion: {{ min: 0, max: 1.0, step: 0.01 }},
            align_accuracy: {{ min: 0, max: 100, step: 1 }}
        }};

        // Initialize curves data arrays
        let rocTraceIdx = null;
        let prTraceIdx = null;

        fn updateMetricLayout() {{
            const key = metricSelector.value;
            const config = metricConfig[key];

            thresholdSlider.min = config.min;
            thresholdSlider.max = config.max;
            thresholdSlider.step = config.step;
            thresholdSlider.value = config.min;

            document.getElementById('slider_min').innerText = config.min;
            document.getElementById('slider_max').innerText = config.max;

            generateStaticPlots(key);
            recalculateMetrics();
        }}

        fn generateStaticPlots(key) {{
            const config = metricConfig[key];
            const steps = 100;
            
            let roc_x = [], roc_y = [];
            let pr_x = [], pr_y = [];

            // Compute ROC and PR curve space step coordinates
            for(let i = 0; i <= steps; i++) {{
                let t = config.min + ((config.max - config.min) * (i / steps));
                
                let tp = countGreaterOrEqual(payload.positives[key], t);
                let fp = countGreaterOrEqual(payload.negatives[key], t);
                
                // Ground truths
                let total_simulated_positives = payload.total_expected_simulated * 2; // R1 + R2
                let total_negatives = payload.negatives[key].length;

                let fn = total_simulated_positives - tp;

                let tpr = total_simulated_positives > 0 ? (tp / total_simulated_positives) : 0;
                let fpr = total_negatives > 0 ? (fp / total_negatives) : 0;
                let precision = (tp + fp) > 0 ? (tp / (tp + fp)) : 1.0;

                roc_x.push(fpr);
                roc_y.push(tpr);
                
                pr_x.push(tpr); // Recall on X axis for PR curve
                pr_y.push(precision);
            }}

            const rocTrace = {{ x: roc_x, y: roc_y, type: 'scatter', mode: 'lines', name: 'ROC Curve', line: {{ color: '#2563eb', width: 2.5 }} }};
            const rocBaseline = {{ x: [0, 1], y: [0, 1], type: 'scatter', mode: 'lines', name: 'Random Guess', line: {{ color: '#94a3b8', dash: 'dash' }} }};
            const rocCurrent = {{ x: [0], y: [0], type: 'scatter', mode: 'markers', name: 'Current Cutoff', marker: {{ color: 'red', size: 10, symbol: 'cross' }} }};

            const prTrace = {{ x: pr_x, y: pr_y, type: 'scatter', mode: 'lines', name: 'PR Curve', line: {{ color: '#16a34a', width: 2.5 }} }};
            const prCurrent = {{ x: [0], y: [0], type: 'scatter', mode: 'markers', name: 'Current Cutoff', marker: {{ color: 'red', size: 10, symbol: 'cross' }} }};

            Plotly.newPlot('roc_plot', [rocTrace, rocBaseline, rocCurrent], {{
                title: 'Receiver Operating Characteristic (ROC)',
                xaxis: {{ title: 'False Positive Rate (FPR)', range: [-0.02, 1.02] }},
                yaxis: {{ title: 'True Positive Rate (TPR)', range: [-0.02, 1.02] }},
                margin: {{ t:50, b:50, l:50, r:20 }},
                showlegend: false
            }}, {{ displaylogo: false }});

            Plotly.newPlot('pr_plot', [prTrace, prCurrent], {{
                title: 'Precision-Recall Curve (PR)',
                xaxis: {{ title: 'Recall (TPR)', range: [-0.02, 1.02] }},
                yaxis: {{ title: 'Precision (PPV)', range: [-0.02, 1.02] }},
                margin: {{ t:50, b:50, l:50, r:20 }},
                showlegend: false
            }}, {{ displaylogo: false }});
        }}

        fn recalculateMetrics() {{
            const key = metricSelector.value;
            const threshold = parseFloat(thresholdSlider.value);

            document.getElementById('current_threshold_view').innerText = threshold.toFixed(metricConfig[key].step % 1 === 0 ? 0 : 2);

            // Compute active counts using fast binary searches
            let tp = countGreaterOrEqual(payload.positives[key], threshold);
            let fp = countGreaterOrEqual(payload.negatives[key], threshold);
            
            // Expected ground truth positives totals across both channels
            let total_sim_positives = payload.total_expected_simulated * 2;
            let total_negatives = payload.negatives[key].length;

            let fn = total_sim_positives - tp;

            // Equation evaluations
            let tpr = total_sim_positives > 0 ? (tp / total_sim_positives) : 0;
            let fpr = total_negatives > 0 ? (fp / total_negatives) : 0;
            let precision = (tp + fp) > 0 ? (tp / (tp + fp)) : 1.0;
            let f1 = (precision + tpr) > 0 ? (2 * (precision * tpr) / (precision + tpr)) : 0;

            // Update DOM text labels
            document.getElementById('cell_tp').innerText = tp.toLocaleString();
            document.getElementById('cell_fp').innerText = fp.toLocaleString();
            document.getElementById('cell_fn').innerText = fn.toLocaleString();

            document.getElementById('stat_precision').innerText = (precision * 100).toFixed(2) + '%';
            document.getElementById('stat_recall').innerText = (tpr * 100).toFixed(2) + '%';
            document.getElementById('stat_fpr').innerText = (fpr * 100).toFixed(2) + '%';
            document.getElementById('stat_f1').innerText = f1.toFixed(4);

            // Update marker icons on Plotly traces
            Plotly.animate('roc_plot', {{
                data: [{{}}, {{}}, {{ x: [fpr], y: [tpr] }}],
                traces: [0, 1, 2]
            }}, {{ transition: {{ duration: 0 }}, frame: {{ duration: 0, redraw: false }} }});

            Plotly.animate('pr_plot', {{
                data: [{{}}, {{ x: [tpr], y: [precision] }}],
                traces: [0, 1]
            }}, {{ transition: {{ duration: 0 }}, frame: {{ duration: 0, redraw: false }} }});
        }}

        metricSelector.addEventListener('change', updateMetricLayout);
        thresholdSlider.addEventListener('input', recalculateMetrics);

        // Bootstrap application setup
        updateMetricLayout();
    </script>
</body>
</html>
"#,
    name = data.report_name,
    expected = data.total_expected_pairs,
    observed = data.total_observed_pairs,
    tp_count = data.positives.mapq.len(),
    fp_count = data.negatives.mapq.len(),
    json_raw = json_raw
    );

    fs::write(html_path, html_content)?;
    Ok(())
}