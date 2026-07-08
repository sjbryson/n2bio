//! n2bio/pfqsim/src/report.rs
//! 

use std::fs::{self, File};
use std::path::PathBuf;
use std::io;
use serde::{Serialize, Deserialize};

use crate::hist::Histogram;

// ============================================================================
// Report data
// ============================================================================

/// Serialized data packaged for HTML report
#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct AnalyzeReportData {
    pub(crate) report_name: String,
    
    // Ground Truth
    pub(crate) total_expected_positives: usize, 
    pub(crate) total_expected_negatives: usize, 
    pub(crate) total_observed: usize,
    
    // Unmapped Counts
    pub(crate) unmapped_true_negatives: usize, 
    pub(crate) unmapped_false_negatives: usize,

    // Mapped Histograms
    // Reads that mapped correctly (Subject to becoming FN if threshold raised)
    pub(crate) true_positives: MetricPayload, 
    
    // Control reads that incorrectly mapped (Subject to becoming TN if threshold raised)
    pub(crate) false_positives_control: MetricPayload, 
    
    // Target reads that incorrectly mapped (FP)
    pub(crate) false_positives_cross: MetricPayload, 
}

impl AnalyzeReportData {
    /// Cascades the trim operation to all target pools
    pub(crate) fn trim(&mut self) {
        self.true_positives.trim();
        self.false_positives_control.trim();
        self.false_positives_cross.trim();
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct MetricPayload {
    pub(crate) mapq: Histogram,
    pub(crate) align_score: Histogram,
    pub(crate) align_length: Histogram,
    pub(crate) as_al: Histogram,
    pub(crate) align_proportion: Histogram,
    pub(crate) align_accuracy: Histogram,
}

impl MetricPayload {
    /// Initializes histograms with sensible boundaries for typical short-read alignments.
    /// Adjust max limits (like align_score 1000) if you are working with long reads!
    pub(crate) fn new() -> Self {
        Self {
            mapq: Histogram::new(0.0, 60.0, 1.0),
            align_score: Histogram::new(0.0, 1000.0, 1.0),
            align_length: Histogram::new(0.0, 1000.0, 1.0),
            as_al: Histogram::new(0.0, 1000.0, 0.1),
            align_proportion: Histogram::new(0.0, 1.0, 0.01),
            align_accuracy: Histogram::new(0.0, 1.0, 0.01),
        }
    }

    /// Cascades the trim operation to all histograms in the payload
    pub(crate) fn trim(&mut self) {
        self.mapq.trim();
        self.align_score.trim();
        self.align_length.trim();
        self.as_al.trim();
        self.align_proportion.trim();
        self.align_accuracy.trim();
    }
}

// ============================================================================
// HTML report
// ============================================================================

/// Generates both the HTML dashboard and the raw JSON data file
pub(crate) fn generate_evaluation_reports(
    data: &AnalyzeReportData, 
    html_path: &PathBuf,
    json_path: &PathBuf,
) -> io::Result<()> {
    
    // 1. Write the JSON file report
    let json_file = File::create(json_path)?;
    serde_json::to_writer_pretty(json_file, data).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("JSON file write error: {}", e))
    })?;

    // 2. Serialize raw JSON string for HTML embedding
    let json_raw: String = serde_json::to_string(data).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("JSON serialization error: {}", e))
    })?;

    // 3. Generate HTML
    let html_content: String = format!(r#"
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>PFQSIM Report: {name}</title>
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
        .dashboard-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-bottom: 15px; }}
        .curves-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-bottom: 15px; }}

        /* Panel Cards */
        .panel {{ background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; }}
        .panel-title {{ text-align: center; font-size: 12px; font-weight: bold; color: #334155; margin-bottom: 15px; text-transform: uppercase; letter-spacing: 0.5px; border-bottom: 1px solid #e2e8f0; padding-bottom: 5px; }}

        /* Interactive Widgets */
        .control-group {{ margin-bottom: 15px; }}
        label {{ display: block; font-size: 12px; font-weight: 600; color: #475569; margin-bottom: 5px; }}
        select, input[type=range] {{ width: 100%; box-sizing: border-box; }}
        .slider-meta {{ display: flex; justify-content: space-between; font-size: 10px; color: #64748b; margin-top: 4px; }}
        .active-val {{ font-weight: bold; color: #2563eb; }}

        /* Matrix Table Styling */
        .matrix-table {{ width: 100%; border-collapse: collapse; text-align: center; font-size: 11px; }}
        .matrix-table th, .matrix-table td {{ padding: 8px; border: 1px solid #cbd5e1; }}
        .matrix-table th {{ background: #f1f5f9; color: #475569; font-weight: 600; }}
        .matrix-label {{ text-align: left; font-weight: bold; background: #f1f5f9; color: #475569; }}
        .cell-tp {{ background: #dcfce7; font-weight: bold; color: #166534; }}
        .cell-fp {{ background: #fee2e2; font-weight: bold; color: #991b1b; }}
        .cell-fn {{ background: #fef9c3; font-weight: bold; color: #854d0e; }}
        .cell-tn {{ background: #f3f4f6; font-weight: bold; color: #374151; }}

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
        <h1>PfqSim Classifier Assessment</h1>
        <div class="subtitle">{name}</div>

        <div class="divider">Dynamic Threshold Report</div>

        <div class="dashboard-grid">
            <div class="panel">
                <div class="panel-title">Summary</div>
                <ul class="stat-list">
                    <li><span class="stat-label">Expected Positive Pairs</span><span class="stat-value">{expected}</span></li>
                    <li><span class="stat-label">Expected Negative Pairs</span><span class="stat-value">{expected_neg}</span></li>
                    <li><span class="stat-label">Observed Primary Alignments</span><span class="stat-value">{observed}</span></li>
                    <li><span class="stat-label">Total Valid Placements</span><span class="stat-value">{tp_count}</span></li>
                    <li><span class="stat-label">Total Error Placements</span><span class="stat-value">{fp_count}</span></li>
                </ul>

                <div class="panel-title" style="margin-top:25px;">Dynamic Metrics</div>
                <ul class="stat-list">
                    <li><span class="stat-label">Precision (PPV)</span><span id="stat_precision" class="stat-value">-</span></li>
                    <li><span class="stat-label">Recall - TPR (Sensitivity)</span><span id="stat_recall" class="stat-value">-</span></li>
                    <li><span class="stat-label">False Positive Rate (FPR)</span><span id="stat_fpr" class="stat-value">-</span></li>
                    <li><span class="stat-label">Specificity (TNR)</span><span id="stat_specificity" class="stat-value">-</span></li>
                    <li><span class="stat-label">Accuracy</span><span id="stat_accuracy" class="stat-value">-</span></li>
                    <li><span class="stat-label">Negative Pred Value (NPV)</span><span id="stat_npv" class="stat-value">-</span></li>
                    <li><span class="stat-label">Prevalence</span><span id="stat_prevalence" class="stat-value">-</span></li>
                    <li><span class="stat-label">F1-Score</span><span id="stat_f1" class="stat-value">-</span></li>
                </ul>
            </div>

            <div class="panel">
                <div class="panel-title">Dynamic Confusion Matrix</div>
                <table class="matrix-table">
                    <tr>
                        <th colspan="2" rowspan="2"></th>
                        <th colspan="2">Predicted</th>
                    </tr>
                    <tr>
                        <th>Positive (&ge; Threshold)</th>
                        <th>Negative (&lt; Threshold)</th>
                    </tr>
                    <tr>
                        <th rowspan="2" style="writing-mode: vertical-lr; transform: rotate(180deg); font-size:10px;">Actual</th>
                        <td class="matrix-label">Positive</td>
                        <td class="cell-tp">
                            <div style="font-size:8px; font-weight:normal; margin-bottom:4px;">True Positive (TP)</div>
                            <div id="cell_tp" style="font-size:16px;">-</div>
                        </td>
                        <td class="cell-fn">
                            <div style="font-size:8px; font-weight:normal; margin-bottom:4px;">False Negative (FN)</div>
                            <div id="cell_fn" style="font-size:16px;">-</div>
                        </td>
                    </tr>
                    <tr>
                        <td class="matrix-label">Negative</td>
                        <td class="cell-fp">
                            <div style="font-size:8px; font-weight:normal; margin-bottom:4px;">False Positive (FP)</div>
                            <div id="cell_fp" style="font-size:16px;">-</div>
                        </td>
                        <td class="cell-tn">
                            <div style="font-size:8px; font-weight:normal; margin-bottom:4px;">True Negative (TN)</div>
                            <div id="cell_tn" style="font-size:16px;">-</div>
                        </td>
                    </tr>
                </table>
                <p style="font-size:10px; color:#64748b; line-height:1.4; margin-top:15px; margin-bottom: 25px;">
                    * <strong>False Negatives (FN)</strong> represent simulated expected alignments that were rejected by the threshold filter or dropped entirely as unmapped.
                </p>

                <div class="panel-title">Threshold Controller</div>
                <div class="control-group">
                    <label for="metric_selector">Target Diagnostic Feature</label>
                    <select id="metric_selector">
                        <option value="align_score">Alignment Score (AS)</option>
                        <option value="align_length">Alignment Length (AL)</option>
                        <option value="as_al">AS / AL (Score per Base)</option>
                        <option value="align_proportion">Alignment Proportion (AL/RL)</option>
                        <option value="align_accuracy">Percent Identity (PI)</option>
                        <option value="mapq">Mapping Quality (MAPQ)</option>
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
            </div>
        </div>

        <div class="curves-row">
            <div id="roc_plot" class="panel" style="height: 450px;"></div>
            <div id="pr_plot" class="panel" style="height: 450px;"></div>
        </div>
    </div>

    <script>
        const payload = {json_raw};
        const metrics = ['align_score', 'align_length', 'as_al', 'align_proportion', 'align_accuracy', 'mapq'];

        // 1. O(N) Pre-computation: Build Suffix Sums on page load
        function initHistograms() {{
            for (let key of metrics) {{
                [payload.true_positives[key], payload.false_positives_control[key], payload.false_positives_cross[key]].forEach(hist => {{
                    let total_bins = hist.counts.length;
                    hist.suffix_sums = new Array(total_bins).fill(0);
                    
                    let running = 0;
                    for (let i = total_bins - 1; i >= 0; i--) {{
                        running += hist.counts[i];
                        hist.suffix_sums[i] = running;
                    }}
                }});
            }}
        }}
        initHistograms();

        // 2. O(1) Threshold Lookup
        function countGreaterOrEqual(histogram, threshold) {{
            let bin_index = Math.round((threshold - histogram.min_val) / histogram.bin_width);
            
            if (bin_index < 0) return histogram.suffix_sums[0] || 0;
            if (bin_index >= histogram.counts.length) return 0;
            
            return histogram.suffix_sums[bin_index];
        }}

        // 3. Get bounds directly from Histogram structs
        function getDynamicConfig(key) {{
            const hist = payload.true_positives[key];
            return {{
                min: hist.min_val,
                max: hist.max_val,
                step: hist.bin_width
            }};
        }}

        // 4. Approximate AUC using trapezoidal rule
        function calculateAUC(x_arr, y_arr) {{
            let area = 0;
            for (let i = 1; i < x_arr.length; i++) {{
                let dx = Math.abs(x_arr[i] - x_arr[i - 1]);
                let dy = (y_arr[i] + y_arr[i - 1]) / 2.0;
                area += dx * dy;
            }}
            return area;
        }}

        const metricSelector = document.getElementById('metric_selector');
        const thresholdSlider = document.getElementById('threshold_slider');

        function updateMetricLayout() {{
            const key = metricSelector.value;
            const config = getDynamicConfig(key);

            thresholdSlider.min = config.min;
            thresholdSlider.max = config.max;
            thresholdSlider.step = config.step;
            thresholdSlider.value = config.min;

            document.getElementById('slider_min').innerText = config.min.toFixed(config.step % 1 === 0 ? 0 : 2);
            document.getElementById('slider_max').innerText = config.max.toFixed(config.step % 1 === 0 ? 0 : 2);

            generateStaticPlots(key, config);
            recalculateMetrics();
        }}

        function generateStaticPlots(key, config) {{
            const steps = 100;
            let roc_x = [], roc_y = [];
            let pr_x = [], pr_y = [];

            let total_sim_positives = payload.total_expected_positives * 2; 
            let total_sim_negatives = payload.total_expected_negatives * 2;
            let total_negatives_pool = total_sim_negatives + countGreaterOrEqual(payload.false_positives_cross[key], config.min); // Using counts for baseline

            for(let i = 0; i <= steps; i++) {{
                let t = config.min + ((config.max - config.min) * (i / steps));
                let tp = countGreaterOrEqual(payload.true_positives[key], t);
                let fp_control = countGreaterOrEqual(payload.false_positives_control[key], t);
                let fp_cross = countGreaterOrEqual(payload.false_positives_cross[key], t);
                let fp = fp_control + fp_cross;
                let tpr = total_sim_positives > 0 ? (tp / total_sim_positives) : 0;
                let fpr = total_negatives_pool > 0 ? (fp / total_negatives_pool) : 0;
                let precision = (tp + fp) > 0 ? (tp / (tp + fp)) : 1.0;

                roc_x.push(fpr);
                roc_y.push(tpr);
                
                pr_x.push(tpr); 
                pr_y.push(precision);
            }}

            // Enforce ROC anchors
            if (roc_x[0] < 1.0 || roc_y[0] < 1.0) {{
                roc_x.unshift(1.0);
                roc_y.unshift(1.0);
            }}
            if (roc_x[roc_x.length - 1] > 0.0 || roc_y[roc_y.length - 1] > 0.0) {{
                roc_x.push(0.0);
                roc_y.push(0.0);
            }}

            // Enforce PR anchor
            if (pr_x[pr_x.length - 1] > 0.0) {{
                pr_x.push(0.0);
                pr_y.push(1.0); 
            }}

            const roc_auc = calculateAUC(roc_x, roc_y);
            const pr_auc = calculateAUC(pr_x, pr_y);

            const rocTrace = {{ x: roc_x, y: roc_y, type: 'scatter', mode: 'lines', name: 'ROC Curve', line: {{ color: '#2563eb', width: 2.5 }} }};
            const rocBaseline = {{ x: [0, 1], y: [0, 1], type: 'scatter', mode: 'lines', name: 'Random Guess', line: {{ color: '#94a3b8', dash: 'dash' }} }};
            const rocCurrent = {{ x: [0], y: [0], type: 'scatter', mode: 'markers', name: 'Current Cutoff', marker: {{ color: 'red', size: 10, symbol: 'cross' }} }};

            const prTrace = {{ x: pr_x, y: pr_y, type: 'scatter', mode: 'lines', name: 'PR Curve', line: {{ color: '#16a34a', width: 2.5 }} }};
            const prCurrent = {{ x: [0], y: [0], type: 'scatter', mode: 'markers', name: 'Current Cutoff', marker: {{ color: 'red', size: 10, symbol: 'cross' }} }};

            Plotly.newPlot('roc_plot', [rocTrace, rocBaseline, rocCurrent], {{
                title: 'Receiver Operating Characteristic (ROC)<br><span style="font-size:14px;color:#64748b;">AUC: ' + roc_auc.toFixed(4) + '</span>',
                xaxis: {{ title: 'False Positive Rate (FPR)', range: [-0.02, 1.02] }},
                yaxis: {{ title: 'True Positive Rate (TPR)', range: [-0.02, 1.02] }},
                margin: {{ t:60, b:50, l:50, r:20 }},
                showlegend: false
            }}, {{ displaylogo: false }});

            Plotly.newPlot('pr_plot', [prTrace, prCurrent], {{
                title: 'Precision-Recall Curve (PR)<br><span style="font-size:14px;color:#64748b;">AUC: ' + pr_auc.toFixed(4) + '</span>',
                xaxis: {{ title: 'Recall (TPR)', range: [-0.02, 1.02] }},
                yaxis: {{ title: 'Precision (PPV)', range: [-0.02, 1.02] }},
                margin: {{ t:60, b:50, l:50, r:20 }},
                showlegend: false
            }}, {{ displaylogo: false }});
        }}

        function recalculateMetrics() {{
            const key = metricSelector.value;
            const config = getDynamicConfig(key);
            const threshold = parseFloat(thresholdSlider.value);

            document.getElementById('current_threshold_view').innerText = threshold.toFixed(config.step % 1 === 0 ? 0 : 2);

            let tp = countGreaterOrEqual(payload.true_positives[key], threshold);
            let fp_control = countGreaterOrEqual(payload.false_positives_control[key], threshold);
            let fp_cross = countGreaterOrEqual(payload.false_positives_cross[key], threshold);
            
            let fp = fp_control + fp_cross;
            
            let total_sim_positives = payload.total_expected_positives * 2;
            let total_sim_negatives = payload.total_expected_negatives * 2;
            let total_negatives_pool = total_sim_negatives + countGreaterOrEqual(payload.false_positives_cross[key], config.min);

            let fn = total_sim_positives - tp;
            let tn = total_negatives_pool - fp; // Safely calculated relative to the unified error pool

            let tpr = total_sim_positives > 0 ? (tp / total_sim_positives) : 0;
            let fpr = total_negatives_pool > 0 ? (fp / total_negatives_pool) : 0;
            let precision = (tp + fp) > 0 ? (tp / (tp + fp)) : 1.0;
            let f1 = (precision + tpr) > 0 ? (2 * (precision * tpr) / (precision + tpr)) : 0;
            
            let specificity = total_negatives_pool > 0 ? (tn / total_negatives_pool) : 0;
            let accuracy = (total_sim_positives + total_negatives_pool) > 0 ? ((tp + tn) / (total_sim_positives + total_negatives_pool)) : 0;
            let npv = (tn + fn) > 0 ? (tn / (tn + fn)) : 1.0;
            let prevalence = (total_sim_positives + total_negatives_pool) > 0 ? (total_sim_positives / (total_sim_positives + total_negatives_pool)) : 0;

            document.getElementById('cell_tp').innerText = tp.toLocaleString();
            document.getElementById('cell_fp').innerText = fp.toLocaleString();
            document.getElementById('cell_fn').innerText = fn.toLocaleString();
            document.getElementById('cell_tn').innerText = tn.toLocaleString();

            document.getElementById('stat_precision').innerText = (precision * 100).toFixed(2) + '%';
            document.getElementById('stat_recall').innerText = (tpr * 100).toFixed(2) + '%';
            document.getElementById('stat_fpr').innerText = (fpr * 100).toFixed(2) + '%';
            document.getElementById('stat_f1').innerText = f1.toFixed(4);
            
            document.getElementById('stat_specificity').innerText = (specificity * 100).toFixed(2) + '%';
            document.getElementById('stat_accuracy').innerText = (accuracy * 100).toFixed(2) + '%';
            document.getElementById('stat_npv').innerText = (npv * 100).toFixed(2) + '%';
            document.getElementById('stat_prevalence').innerText = (prevalence * 100).toFixed(2) + '%';

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
    expected = data.total_expected_positives,
    expected_neg = data.total_expected_negatives,
    observed = data.total_observed,
    tp_count = data.true_positives.align_score.total_count(),
    fp_count = data.false_positives_control.align_score.total_count() + data.false_positives_cross.align_score.total_count(),
    json_raw = json_raw
    );

    fs::write(html_path, html_content)?;
    Ok(())
}