//! n2core/src/align.rs
//! 


// ============================================================================
// Needleman-Wunch Alignment
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub enum EditOp {
    Match,
    Substitution(u8),
    Insertion(u8),
    Deletion,
}

/// Holds scoring parameters for alignment.
#[derive(Debug, Clone)]
pub struct AlignmentScores {
    pub match_score:      i32,
    pub mismatch_penalty: i32, // Passed as a negative value (e.g., -1)
    pub gap_penalty:      i32, // Passed as a negative value (e.g., -3)
}

impl Default for AlignmentScores {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_penalty: -1,
            gap_penalty: -3,
        }
    }
}

/// 1D flattened Needleman-Wunsch implementation 
pub fn align_seqs(
    ref_seq: &[u8], 
    qry_seq: &[u8], 
    scores: &AlignmentScores
) -> Vec<EditOp> {
    let r_len: usize = ref_seq.len();
    let q_len: usize = qry_seq.len();

    let cols: usize = q_len + 1;
    let rows: usize = r_len + 1;

    let mut dp: Vec<i32> = vec![0i32; rows * cols];
    let mut tb: Vec<u8>  = vec![0u8; rows * cols]; // 0=Diag, 1=Up, 2=Left

    for i in 1..rows {
        dp[i * cols] = i as i32 * scores.gap_penalty;
        tb[i * cols] = 1; 
    }
    for j in 1..cols {
        dp[j] = j as i32 * scores.gap_penalty;
        tb[j] = 2; 
    }

    for i in 1..rows {
        for j in 1..cols {
            let diag_score: i32 = dp[(i - 1) * cols + (j - 1)] + 
                if ref_seq[i - 1] == qry_seq[j - 1] { 
                    scores.match_score 
                } else { 
                    scores.mismatch_penalty 
                };
                
            let up_score: i32 = dp[(i - 1) * cols + j] + scores.gap_penalty;
            let left_score: i32 = dp[i * cols + (j - 1)] + scores.gap_penalty;

            let max_score: i32 = diag_score.max(up_score).max(left_score);
            dp[i * cols + j] = max_score;

            // Tie-breaking preference: Match/Mismatch > Deletion > Insertion
            if max_score == diag_score {
                tb[i * cols + j] = 0;
            } else if max_score == up_score {
                tb[i * cols + j] = 1;
            } else {
                tb[i * cols + j] = 2;
            }
        }
    }

    let mut i: usize = r_len;
    let mut j: usize = q_len;
    let mut ops: Vec<EditOp> = Vec::with_capacity(r_len.max(q_len));

    while i > 0 || j > 0 {
        match tb[i * cols + j] {
            0 => {
                if ref_seq[i - 1] == qry_seq[j - 1] {
                    ops.push(EditOp::Match);
                } else {
                    ops.push(EditOp::Substitution(qry_seq[j - 1]));
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(EditOp::Deletion);
                i -= 1;
            }
            2 => {
                ops.push(EditOp::Insertion(qry_seq[j - 1]));
                j -= 1;
            }
            _ => unreachable!(),
        }
    }

    ops.reverse();
    ops
}

// ============================================================================
// Partial Order Alignment (Graph-based Needleman-Wunsch)
// ============================================================================

/// Represents a single cell in the Partial Order Alignment DP matrix.
#[derive(Clone, Debug)]
pub struct PoaCell {
    pub score: i32,
    pub op: EditOp,
    // Tracks the column index (j) of the predecessor node that yielded the best score.
    // Essential for tracing back through graph branches.
    pub pred_j: Option<usize>, 
}

/// Aligns a query sequence against a DAG (Directed Acyclic Graph) subgraph.
/// 
/// * `query` - The linear sequence of bases to align.
/// * `node_bases` - The nucleotides corresponding to the topologically sorted graph nodes.
/// * `predecessors` - A list where `predecessors[j]` contains the indices `p` (where p < j) 
///                    of the incoming parent nodes in the sorted array.
/// * `scores` - Standard substitution/gap penalties.
/// 
/// Returns a tuple of: 
/// 1. The optimal `EditOp` sequence
/// 2. The sequence of graph array indices (`j`) traversed (the specific path chosen through the bubble)
pub fn align_poa(
    query: &[u8],
    node_bases: &[u8],
    predecessors: &[Vec<usize>],
    scores: &AlignmentScores,
) -> (Vec<EditOp>, Vec<usize>) {
    let q_len: usize = query.len();
    let g_len: usize = node_bases.len();

    if q_len == 0 || g_len == 0 {
        return (Vec::new(), Vec::new());
    }

    // matrix[i][j]: i = query index, j = topologically sorted graph node index
    let mut matrix: Vec<Vec<PoaCell>> = vec![vec![PoaCell { score: 0, op: EditOp::Match, pred_j: None }; g_len]; q_len + 1];

    // 1. Initialize top row (i=0): Query is empty, so these are Deletions tracing through the graph
    for j in 1..g_len {
        let mut best_score: i32 = std::i32::MIN;
        let mut best_p: usize = 0;
        
        for &p in &predecessors[j] {
            // Gap penalty applies for moving through graph nodes without a query base
            let score: i32 = matrix[0][p].score + scores.gap_penalty; 
            if score > best_score {
                best_score = score;
                best_p = p;
            }
        }
        matrix[0][j] = PoaCell { score: best_score, op: EditOp::Deletion, pred_j: Some(best_p) };
    }

    // 2. Initialize first column (j=0): Graph is empty, so these are Insertions of query bases
    for i in 1..=q_len {
        matrix[i][0] = PoaCell {
            score: matrix[i-1][0].score + scores.gap_penalty,
            op: EditOp::Insertion(query[i-1]),
            pred_j: None, // Insertions step down in 'i'; 'j' stays the same
        };
    }

    // 3. Populate DP Matrix
    for i in 1..=q_len {
        for j in 1..g_len {
            let mut best_score: i32 = std::i32::MIN;
            let mut best_op: EditOp = EditOp::Match;
            let mut best_pred_j: Option<usize> = None;

            // Option A: Insertion (Query advances, Graph stays on the same node j)
            let ins_score: i32 = matrix[i-1][j].score + scores.gap_penalty;
            if ins_score > best_score {
                best_score = ins_score;
                best_op = EditOp::Insertion(query[i-1]);
                best_pred_j = Some(j); 
            }

            // Check all valid predecessors in the graph for Match/Mismatch and Deletion
            for &p in &predecessors[j] {
                // Option B: Deletion (Graph advances from p->j, Query stays on i)
                let del_score: i32 = matrix[i][p].score + scores.gap_penalty;
                if del_score > best_score {
                    best_score = del_score;
                    best_op = EditOp::Deletion;
                    best_pred_j = Some(p);
                }

                // Option C: Match/Mismatch (Both advance: Query i-1->i, Graph p->j)
                let is_match: bool = query[i-1] == node_bases[j];
                let match_score: i32 = if is_match { scores.match_score } else { scores.mismatch_penalty };
                let diag_score: i32 = matrix[i-1][p].score + match_score;
                
                if diag_score > best_score {
                    best_score = diag_score;
                    best_op = if is_match { EditOp::Match } else { EditOp::Substitution(query[i-1]) };
                    best_pred_j = Some(p);
                }
            }

            matrix[i][j] = PoaCell { score: best_score, op: best_op, pred_j: best_pred_j };
        }
    }

    // 4. Traceback starting from the bottom-right cell
    let mut ops: Vec<EditOp> = Vec::new();
    let mut path: Vec<usize> = Vec::new(); // Tracks which graph indices were traversed

    let mut curr_i: usize = q_len;
    let mut curr_j: usize = g_len - 1;

    while curr_i > 0 || curr_j > 0 {
        if curr_i == 0 && curr_j == 0 { break; }

        let cell: &PoaCell = &matrix[curr_i][curr_j];
        ops.push(cell.op.clone());

        match cell.op {
            EditOp::Insertion(_) => {
                // Query was consumed, Graph node stays the same
                curr_i -= 1;
            }
            EditOp::Deletion => {
                // Graph node was consumed, Query stays the same
                path.push(curr_j);
                curr_j = cell.pred_j.expect("Deletion must have a predecessor");
            }
            EditOp::Match | EditOp::Substitution(_) => {
                // Both were consumed
                path.push(curr_j);
                curr_i -= 1;
                curr_j = cell.pred_j.expect("Match/Sub must have a predecessor");
            }
        }
    }

    ops.reverse();
    path.reverse();

    (ops, path)
}