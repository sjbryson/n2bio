//! n2bio/fastcov/src/main.rs

use clap::Parser;
use crossbeam::channel::bounded;
use std::io::{self, BufRead};
use std::thread;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;
use std::io::Write;
use std::collections::HashMap;
use rusqlite::Connection;
use n2core::sam::{SamReader, SamStr, SamFields, SamFlags, SamTags, AlignmentStats};


#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "High-performance SAM filter", long_about = None)]
struct Args {

    /// Number of worker threads for parsing
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,

    /// Name of the run/sample for the JSON report
    #[arg(short = 'r', long, required = true)]
    run_name: String,

    /// Optional: Min Alignment Proportion - sam.calculate_alignment_proportion()
    #[arg(long)]
    min_ap: Option<f32>,
    
    /// Optional: Min Percent Identity - sam.calculate_alignment_accuracy()
    #[arg(long)]
    min_pi: Option<f32>,
    
    /// Optional: Min Alignment Score - sam.get_int_tag("AS")
    #[arg(long)]
    min_as: Option<i32>,

    /// Optional: Min Alignment Lenth - sam.calculate_alignment_length()
    #[arg(long)]
    min_al: Option<u32>,
    
    /// Optional: Min AS/AL score - sam.calculate_as_al()
    #[arg(long)]
    min_sl: Option<f32>,

    /// Optional: Min MAPQ score - sam.calculate_as_al()
    #[arg(long)]
    min_mq: Option<u32>,

    /// Optional path to an SQLite taxonomy database (see vref2db)
    #[arg(long)]
    db: Option<String>,
}

/// Filter logic for whether an alignment passes - aligned well in this use case.
fn sam_filter(sam: &SamStr, args: &Args) -> bool {
   
    if sam.is_mapped() {
    
        // Evaluate optional filters
        if args.min_ap.is_some_and(|min: f32| sam.calculate_alignment_proportion().is_some_and(|val: f32| val >= min)) {
            return true;
        }
        if args.min_pi.is_some_and(|min: f32| sam.calculate_alignment_accuracy().is_some_and(|val: f32| val >= min)) {
            return true;
        }
        if args.min_as.is_some_and(|min: i32| sam.get_int_tag("AS").is_some_and(|val: i32| val >= min)) {
            return true;
        }
        if args.min_al.is_some_and(|min: u32| sam.calculate_alignment_length().is_some_and(|val: u32| val >= min)) {
            return true;
        }
        if args.min_sl.is_some_and(|min: f32| sam.calculate_as_al().is_some_and(|val: f32| val >= min)) {
            return true;
        }
        if args.min_mq.is_some_and(|min: u32| sam.mapq() >= min) {
            return true;
        }
    }
    // If it is mapped but didn't pass any of the specified min thresholds or is unmapped
    false
}

#[derive(Debug, Clone)]
struct RefStats {
    ref_name: String,
    ref_length: usize,
    num_primary: u32,
    num_secondary: u32,
    num_primary_concordant: u32,
    num_secondary_concordant: u32,
    primary_coverage: Vec<u32>,
    secondary_coverage: Vec<u32>,
    mismatch_coverage: Vec<u32>,
}

impl RefStats {
    fn new(ref_name: String, ref_length: usize) -> Self {
        Self {
            ref_name,
            ref_length,
            num_primary: 0,
            num_secondary: 0,
            num_primary_concordant: 0,
            num_secondary_concordant: 0,
            primary_coverage: vec![0; ref_length],
            secondary_coverage: vec![0; ref_length],
            mismatch_coverage: vec![0; ref_length],
        }
    }
}

enum PipelineMsg {
    RefLength(String, usize),
    Hit {
        rname: String,
        is_primary: bool,
        is_concordant: bool,
        start_pos: usize,
        cigar: String,
    },
}
fn main() -> io::Result<()> {
    let start_time: Instant = Instant::now();
    let args: Args          = Args::parse();
    let out_json: String    = format!("{}.json", args.run_name);

    // Global Counters
    let total_alignments: Arc<AtomicU64> = Arc::new(AtomicU64::new(0));
    let passed_primary: Arc<AtomicU64>   = Arc::new(AtomicU64::new(0));
    let passed_secondary: Arc<AtomicU64> = Arc::new(AtomicU64::new(0));

    let (line_tx, line_rx) = bounded::<String>(10_000);
    let (hit_tx, hit_rx) = bounded::<PipelineMsg>(10_000);

    // ---------------------------------------------------------
    // AGGREGATOR THREAD: Collects Stats lock-free
    // ---------------------------------------------------------
    let aggregator_handle: thread::JoinHandle<HashMap<String, RefStats>> = thread::spawn(move || -> HashMap<String, RefStats> {
        let mut ref_map: HashMap<String, RefStats> = HashMap::new();

        for msg in hit_rx {
            match msg {
                PipelineMsg::RefLength(rname, length) => {
                    ref_map.entry(rname.clone()).or_insert_with(|| RefStats::new(rname, length));
                }
                PipelineMsg::Hit { rname, is_primary, is_concordant, start_pos, cigar } => {
                    if let Some(stats) = ref_map.get_mut(&rname) {
                        
                        if is_primary {
                            stats.num_primary += 1;
                            if is_concordant { stats.num_primary_concordant += 1; }
                        } else {
                            stats.num_secondary += 1;
                            if is_concordant { stats.num_secondary_concordant += 1; }
                        }

                        let mut current_ref_pos: usize = start_pos.saturating_sub(1);
                        let mut current_num: usize = 0;
                        let ref_len: usize = stats.ref_length;

                        for c in cigar.chars() {
                            if c.is_ascii_digit() {
                                current_num = current_num * 10 + c.to_digit(10).unwrap() as usize;
                            } else {
                                match c {
                                    'M' | '=' => {
                                        let end: usize = std::cmp::min(current_ref_pos + current_num, ref_len);
                                        for i in current_ref_pos..end {
                                            if is_primary { stats.primary_coverage[i] += 1; } 
                                            else { stats.secondary_coverage[i] += 1; }
                                        }
                                        current_ref_pos += current_num;
                                    }
                                    'X' => {
                                        // X is a sequence mismatch (requires minimap2 --eqx)
                                        let end: usize = std::cmp::min(current_ref_pos + current_num, ref_len);
                                        for i in current_ref_pos..end {
                                            if is_primary { stats.primary_coverage[i] += 1; } 
                                            else { stats.secondary_coverage[i] += 1; }
                                            stats.mismatch_coverage[i] += 1;
                                        }
                                        current_ref_pos += current_num;
                                    }
                                    'D' | 'N' => current_ref_pos += current_num,
                                    'I' | 'S' | 'H' | 'P' => {}
                                    _ => {}
                                }
                                current_num = 0;
                            }
                        }
                    }
                }
            }
        }
        ref_map
    });

    // ---------------------------------------------------------
    // WORKER THREADS: Parse & Filter
    // ---------------------------------------------------------
    let mut worker_handles: Vec<thread::JoinHandle<()>> = Vec::with_capacity(args.threads);
    for _ in 0..args.threads {
        let rx: crossbeam::channel::Receiver<String>      = line_rx.clone();
        let h_tx: crossbeam::channel::Sender<PipelineMsg> = hit_tx.clone();
        
        let t_align: Arc<AtomicU64> = Arc::clone(&total_alignments);
        let p_prim: Arc<AtomicU64>  = Arc::clone(&passed_primary);
        let p_sec: Arc<AtomicU64>   = Arc::clone(&passed_secondary);
        let worker_args: Args       = args.clone();

        let handle: thread::JoinHandle<()> = thread::spawn(move || {
            for line in rx {
                let sam: SamStr<'_> = SamStr::new(&line);
                t_align.fetch_add(1, Ordering::Relaxed);
                
                if sam_filter(&sam, &worker_args) {
                    let is_primary: bool = sam.is_primary();
                    if is_primary {
                        p_prim.fetch_add(1, Ordering::Relaxed);
                    } else {
                        p_sec.fetch_add(1, Ordering::Relaxed);
                    }

                    // Flag 0x2 checks if reads are properly aligned/concordant
                    let is_concordant: bool = sam.is_proper(); 

                    let _ = h_tx.send(PipelineMsg::Hit {
                        rname: sam.rname().to_string(),
                        is_primary,
                        is_concordant,
                        start_pos: sam.pos() as usize,
                        cigar: sam.cigar().to_string(),
                    });
                }
            }
        });
        worker_handles.push(handle);
    }

    // ---------------------------------------------------------
    // MAIN THREAD: Tee off of input stream (Stdin -> Stdout + Workers)
    // ---------------------------------------------------------
    let mut sam_reader: SamReader = SamReader::from_stdin();
    let mut stdout: io::StdoutLock<'_> = io::stdout().lock();
    let mut line_buffer: String = String::new();

    while let Ok(bytes) = sam_reader.reader.read_line(&mut line_buffer) {
        if bytes == 0 { break; }
        
        // 1. Pass the stream to stdout
        stdout.write_all(line_buffer.as_bytes())?;

        // 2. Process for JSON report
        let clean_line: String = line_buffer.trim_end().to_string();
        
        if clean_line.starts_with('@') {
            if clean_line.starts_with("@SQ") {
                let mut sn: String = String::new();
                let mut ln: usize = 0;
                for part in clean_line.split('\t') {
                    if let Some(name) = part.strip_prefix("SN:") { sn = name.to_string(); }
                    if let Some(len)  = part.strip_prefix("LN:") { ln = len.parse().unwrap_or(0); }
                }
                if !sn.is_empty() && ln > 0 {
                    let _ = hit_tx.send(PipelineMsg::RefLength(sn, ln));
                }
            }
        } else {
            if line_tx.send(clean_line).is_err() { break; }
        }
        line_buffer.clear();
    }

    drop(line_tx);
    drop(hit_tx);

    for handle in worker_handles { handle.join().unwrap(); }
    let ref_map: HashMap<String, RefStats> = aggregator_handle.join().unwrap();

    // ---------------------------------------------------------
    // JSON EXPORT
    // ---------------------------------------------------------
    let mut num_refs_primary: u32 = 0;
    let mut num_refs_secondary: u32 = 0;
    let mut coverage_stats: Vec<serde_json::Value> = Vec::new();


    // Setup SQLite Connection & Prepare Statement (if --db is provided)
    let conn: Option<Connection> = args.db.as_ref().map(|db_path| {
        rusqlite::Connection::open(db_path).expect("Failed to open SQLite database")
    });

    let mut stmt: Option<rusqlite::Statement<'_>> = conn.as_ref().map(|c| {
        c.prepare("SELECT * FROM viral_taxonomy WHERE accession = ?")
        .expect("Failed to prepare SQL statement")
    });

    // Array of ranks
    let ranks: [&str; 8] = [
        "realm", "kingdom", "phylum", "class", "order", "family","genus", "species"
    ];
    

    for (_, stats) in ref_map {
        if stats.num_primary   > 0 { num_refs_primary += 1; }
        if stats.num_secondary > 0 { num_refs_secondary += 1; }
    
        // Calculate sums
        let total_reads: u32      = stats.num_primary + stats.num_secondary;
        let primary_sum: u64      = stats.primary_coverage.iter().map(|&x| x as u64).sum();
        let secondary_sum: u64    = stats.secondary_coverage.iter().map(|&x| x as u64).sum();
        let total_mismatches: u64 = stats.mismatch_coverage.iter().map(|&x| x as u64).sum();
        let total_bases: u64      = primary_sum + secondary_sum;
        let ref_len_f64: f64      = stats.ref_length as f64;

        // Calculate depths
        let average_coverage: f64     = if ref_len_f64 > 0.0 { total_bases as f64 / ref_len_f64 } else { 0.0 };
        let primary_avg_coverage: f64 = if ref_len_f64 > 0.0 { primary_sum as f64 / ref_len_f64 } else { 0.0 };

        // Calculate breadths
        let primary_bases_cov: u32 = stats.primary_coverage.iter().filter(|&&p| p > 0).count() as u32;
        let bases_covered: u32     = stats.primary_coverage.iter()
            .zip(stats.secondary_coverage.iter())
            .filter(|&(p, s)| *p > 0 || *s > 0)
            .count() as u32;
        let bases_mismatched: u32  = stats.mismatch_coverage.iter().filter(|&&m| m > 0).count() as u32;
        let aligned_coverage: f64  = bases_covered as f64 / ref_len_f64;

        // Calculate Identity
        let aligned_identity: f64  = if total_bases > 0 {
            (total_bases.saturating_sub(total_mismatches)) as f64 / total_bases as f64
        } else {
            0.0
        };

        // Calculate RPK
        let k: f64   = ref_len_f64 / 1000.0;
        let rpk: f64 = total_reads as f64 / k;

        // Query the SQLite Database for the Accession
        let lineage_json: serde_json::Value = if let Some(ref mut statement) = stmt {
            // Query the DB.
            let query_result: Result<serde_json::Value, rusqlite::Error> = statement.query_row([&stats.ref_name], |row| {
                let mut map: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();

                // Get the base taxonomy ID and name first
                let base_tax_id: Option<i64> = row.get("tax_id").ok();
                let base_name: Option<String> = row.get("name").ok();
                if let (Some(id), Some(name)) = (base_tax_id, base_name) {
                    map.insert("organism".to_string(), serde_json::json!({ "tax_id": id, "name": name }));
                }

                // Loop through the ranks
                for rank in ranks.iter() {
                    let id_col: String = format!("{}_id", rank);
                    let name_col: String = format!("{}_name", rank);
                    let id: Option<i64> = row.get(id_col.as_str()).ok();
                    let name: Option<String> = row.get(name_col.as_str()).ok();

                    // Only add the rank to JSON if both ID and Name exist in the row
                    if let (Some(rank_id), Some(rank_name)) = (id, name) {
                        map.insert(rank.to_string(), serde_json::json!({ "tax_id": rank_id, "name": rank_name }));
                    }
                }
                
                Ok(serde_json::Value::Object(map))
            });

        // If the query was successful, return it. Otherwise return null.
            query_result.unwrap_or(serde_json::json!(null))
        
        } else {
            // If --db was not provided (stmt is None), return null.
            serde_json::json!(null)
        };

        // Helper to format the arrays as comma-separated Strings
        let vec_to_string = |v: &Vec<u32>| v.iter().map(|n| n.to_string()).collect::<Vec<String>>().join(",");

        if stats.num_primary > 0 {
            coverage_stats.push(serde_json::json!({
                stats.ref_name: {
                    "aligned_bases"                      : total_bases,
                    "aligned_cov_ratio"                  : f64::trunc(aligned_coverage * 10000.0) / 10000.0,
                    "aligned_cov_identity"               : f64::trunc(aligned_identity * 10000.0) / 10000.0,
                    "aligned_pos_avg_cov"                : f64::trunc(total_bases as f64 / bases_covered as f64 * 10000.0) / 10000.0,
                    "aligned_pos_covered"                : bases_covered,
                    "aligned_pos_mismatched"             : bases_mismatched,
                    "aligned_pos_primary"                : primary_bases_cov,
                    "aligned_reads"                      : total_reads,
                    "aligned_reads_primary"              : stats.num_primary,
                    "aligned_reads_secondary"            : stats.num_secondary,
                    "aligned_bases_mismatches"           : total_mismatches,
                    "aligned_reads_primary_concordant"   : stats.num_primary_concordant,
                    "aligned_reads_secondary_concordant" : stats.num_secondary_concordant,
                    "average_coverage"                   : f64::trunc(average_coverage * 10000.0) / 10000.0,
                    "average_coverage_primary"           : f64::trunc(primary_avg_coverage * 10000.0) / 10000.0,
                    "ref_length"                         : stats.ref_length,
                    "rpk"                                : f64::trunc(rpk * 10000.0) / 10000.0,
                    "virus_lineage"                      : lineage_json,
                    "x_coverage" : {
                        "primary_coverage"   : vec_to_string(&stats.primary_coverage),
                        "secondary_coverage" : vec_to_string(&stats.secondary_coverage),
                        "mismatch_coverage"  : vec_to_string(&stats.mismatch_coverage),
                    }
                    
                }
            }));
        }
    }

    let summary: serde_json::Value = serde_json::json!({
        "1-run_stats": {
            "run_name"                    : args.run_name,
            "total_run_time_seconds"      : start_time.elapsed().as_secs_f64(),
            "total_alignments"            : total_alignments.load(Ordering::Relaxed),
            "passed_primary_alignments"   : passed_primary.load(Ordering::Relaxed),
            "passed_secondary_alignments" : passed_secondary.load(Ordering::Relaxed),
            "num_refs_primary"            : num_refs_primary,
            "num_refs_secondary"          : num_refs_secondary
        },
        "2-coverage_stats"                : coverage_stats
    });

    // Write JSON to file
    let mut json_file: std::fs::File = std::fs::File::create(&out_json)?;
    json_file.write_all(serde_json::to_string_pretty(&summary).unwrap().as_bytes())?;

    Ok(())
}