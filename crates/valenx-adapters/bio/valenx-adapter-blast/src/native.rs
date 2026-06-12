//! Native Rust BLAST-class sequence search for the BLAST adapter.
//!
//! Uses `valenx_align::search`: k-mer seeding + seed-and-extend with
//! Karlin-Altschul E-values, plus `smith_waterman` on each HSP for
//! exact alignment statistics (pident, mismatches, gap opens).
//!
//! ## Limitations vs real BLAST+
//!
//! - Database input must be a FASTA file (not a preformatted BLAST
//!   database). In `prepare()`, the adapter looks for a companion FASTA
//!   next to the database prefix (same stem + `.fa`, `.fasta`, `.fna`,
//!   `.faa`, or the prefix itself if it is a file). Formatted databases
//!   (`.nhr/.nin/.nsq` etc.) require the real BLAST+ binary.
//! - Only tabular format 6 output is produced; pairwise/XML output
//!   requires the real binary.
//! - Translated BLAST programs (blastx, tblastn, tblastx) run the
//!   search using the raw sequence alphabet (no 6-frame translation);
//!   sensitivity for translated comparisons is lower than real BLAST.
//!
//! ## Reference validation
//!
//! A self-search of a 20-aa protein sequence (ACDEFGHIKLMNPQRSTVWY) vs
//! a 2-sequence database containing that exact sequence + a random one
//! produces exactly one HSP at 100% identity with E-value near zero.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Instant;

use valenx_align::{
    search::{KarlinAltschul, KmerIndex, SeedParams, SeedSearch},
    smith_waterman, ScoringScheme,
};
use valenx_core::{error::RunPhase, AdapterError, RunContext, RunReport};

/// Parameters stored in `native_params.toml` by `prepare()`.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct NativeBlastParams {
    /// Absolute path to the resolved query multi-FASTA file.
    pub query_path: String,
    /// Absolute path to the resolved database FASTA file.
    pub db_fasta_path: String,
    /// E-value cutoff (hits with e_value > this are discarded).
    pub evalue: f64,
    /// Maximum number of hits to report per query (0 = unlimited).
    pub max_hits: usize,
    /// Output filename (relative to workdir).
    pub output_name: String,
}

/// File name used to persist the native params in the workdir.
pub const PARAMS_FILE: &str = "native_params.toml";

/// Command sentinel: run native search, don't spawn a child process.
pub const NATIVE_SENTINEL: &str = "valenx:native:blast";

/// Writes `native_params.toml` to `workdir`.
pub fn write_params(workdir: &Path, params: &NativeBlastParams) -> Result<(), AdapterError> {
    let s = toml::to_string(params)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("serialize native_params: {e}")))?;
    let dest = workdir.join(PARAMS_FILE);
    valenx_core::io_caps::atomic_write_str(&dest, &s)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", dest.display())))
}

/// Reads `native_params.toml` from `workdir`.
pub fn read_params(workdir: &Path) -> Result<NativeBlastParams, AdapterError> {
    let path = workdir.join(PARAMS_FILE);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", path.display())))?;
    toml::from_str(&text)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("parse {}: {e}", path.display())))
}

/// Tries to find a FASTA file for `db_prefix`. Checks, in order:
/// 1. `db_prefix` itself if it is a file
/// 2. `db_prefix` + `.fa`, `.fasta`, `.fna`, `.faa`
///
/// Returns the first path that exists as a file, or `None`.
pub fn find_db_fasta(db_prefix: &Path) -> Option<PathBuf> {
    if db_prefix.is_file() {
        return Some(db_prefix.to_path_buf());
    }
    for ext in &["fa", "fasta", "fna", "faa"] {
        let candidate = db_prefix.with_extension(ext);
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

/// Runs native BLAST-class seed-and-extend search.
///
/// Produces tabular format 6 output at `workdir/<params.output_name>`.
pub fn run_native(workdir: &Path, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
    let params = read_params(workdir)?;
    let start = Instant::now();

    ctx.report_progress(5.0, "native BLAST — loading query");
    let query_text = read_text_file(&params.query_path)?;
    let query_records = parse_fasta_simple(&query_text)?;

    ctx.report_progress(10.0, "native BLAST — loading database");
    let db_text = read_text_file(&params.db_fasta_path)?;
    let db_records = parse_fasta_simple(&db_text)?;

    if db_records.is_empty() {
        return Err(AdapterError::InvalidCase {
            case_path: PathBuf::from(&params.db_fasta_path),
            reason: "database FASTA contains no sequences".to_string(),
        });
    }

    ctx.report_progress(15.0, "native BLAST — detecting alphabet");

    // Use combined query + db sequences to determine scoring scheme.
    let mut all_for_detection = query_records.clone();
    all_for_detection.extend(db_records.iter().take(5).cloned());
    let scheme = detect_scheme(&all_for_detection);

    let is_protein = scheme.matrix.name().contains("BLOSUM");

    // GATE: native nucleotide E-value statistics are NOT yet correctly
    // calibrated. `dna_ungapped()` Karlin-Altschul lambda=1.374 is derived
    // for a +1/-3 scoring system, but `dna_default()` scores +5/-4 — so DNA
    // E-values/bit-scores would be reported on a wildly wrong scale (off by
    // tens of orders of magnitude). Rather than emit confidently-wrong
    // significance values silently, refuse nucleotide search in native mode.
    // The real BLAST+ `blastn` binary still handles nucleotide search.
    // TODO(bio): add a KA parameter set calibrated for the nucleotide
    // scoring system, validate against NCBI blastn, then lift this gate.
    if !is_protein {
        return Err(AdapterError::Run {
            exit_code: 2,
            stderr: "native BLAST supports PROTEIN search only: nucleotide \
                     Karlin-Altschul E-value statistics are not yet calibrated and \
                     would be reported on a wrong scale. Install the BLAST+ binary \
                     (blastn) for nucleotide search, or submit a protein query."
                .to_string(),
            phase: RunPhase::Startup,
        });
    }

    let k = if is_protein { 5 } else { 11 };
    let ka = if is_protein {
        KarlinAltschul::blosum62_ungapped()
    } else {
        KarlinAltschul::dna_ungapped()
    };

    ctx.report_progress(20.0, "native BLAST — building database index");

    let db_seqs_owned: Vec<Vec<u8>> = db_records.iter().map(|(_, s)| s.clone()).collect();
    let db_seqs: Vec<&[u8]> = db_seqs_owned.iter().map(|s| s.as_slice()).collect();

    let index = KmerIndex::build_many(&db_seqs, k).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("k-mer index build failed: {e}"),
        phase: RunPhase::Startup,
    })?;

    let searcher =
        SeedSearch::new(&index, db_seqs.clone(), &scheme, SeedParams::default()).with_stats(ka);

    ctx.report_progress(30.0, "native BLAST — searching");

    let mut output_lines: Vec<String> = Vec::new();
    let total = query_records.len();

    for (qi, (qname, qseq)) in query_records.iter().enumerate() {
        if qi % 10 == 0 {
            let pct = 30.0 + (qi as f32 / total as f32) * 60.0;
            ctx.report_progress(pct, &format!("native BLAST — query {}/{}", qi + 1, total));
        }

        let qid = qname
            .as_deref()
            .unwrap_or("query")
            .split_whitespace()
            .next()
            .unwrap_or("query");

        let mut hsps = searcher.search(qseq);
        // Filter by E-value.
        hsps.retain(|h| h.e_value.is_none_or(|e| e <= params.evalue));
        // Sort by E-value ascending (best first).
        hsps.sort_by(|a, b| {
            a.e_value
                .unwrap_or(f64::MAX)
                .partial_cmp(&b.e_value.unwrap_or(f64::MAX))
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let limit = if params.max_hits == 0 {
            hsps.len()
        } else {
            hsps.len().min(params.max_hits)
        };

        for hsp in &hsps[..limit] {
            let (sid, sseq) = &db_records[hsp.seq_id];
            let sname = sid
                .as_deref()
                .unwrap_or("subject")
                .split_whitespace()
                .next()
                .unwrap_or("subject");

            // Extract the segments and run Smith-Waterman for exact stats.
            let (qs, qe) = hsp.query_span;
            let (ss, se) = hsp.target_span;
            let qs_clamped = qs.min(qseq.len());
            let qe_clamped = qe.min(qseq.len());
            let ss_clamped = ss.min(sseq.len());
            let se_clamped = se.min(sseq.len());

            let qslice = &qseq[qs_clamped..qe_clamped];
            let sslice = &sseq[ss_clamped..se_clamped];

            let (pident, aln_len, mismatches, gap_opens) = if qslice.is_empty() || sslice.is_empty()
            {
                (0.0_f64, 0usize, 0usize, 0usize)
            } else {
                match smith_waterman(qslice, sslice, &scheme) {
                    Ok(aln) => {
                        let stats = aln.stats(&scheme.matrix);
                        let pct = if stats.columns == 0 {
                            0.0
                        } else {
                            stats.identities as f64 / stats.columns as f64 * 100.0
                        };
                        (
                            pct,
                            stats.columns,
                            stats.columns.saturating_sub(stats.identities + stats.gaps),
                            stats.gap_opens,
                        )
                    }
                    Err(_) => {
                        // SW failed on this segment; use rough approximation.
                        let len = qslice.len().min(sslice.len());
                        let matches = qslice
                            .iter()
                            .zip(sslice.iter())
                            .filter(|(a, b)| a == b)
                            .count();
                        let pct = if len == 0 {
                            0.0
                        } else {
                            matches as f64 / len as f64 * 100.0
                        };
                        (pct, len, len - matches, 0)
                    }
                }
            };

            let evalue = hsp.e_value.unwrap_or(f64::MAX);
            let bitscore = hsp.bit_score.unwrap_or(0.0);

            // BLAST tabular format 6: 1-based, inclusive coordinates.
            output_lines.push(format!(
                "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
                qid,
                sname,
                pident,
                aln_len,
                mismatches,
                gap_opens,
                qs_clamped + 1, // qstart (1-based)
                qe_clamped,     // qend (1-based inclusive = 0-based exclusive)
                ss_clamped + 1, // sstart
                se_clamped,     // send
                evalue,
                bitscore,
            ));
        }
    }

    ctx.report_progress(95.0, "native BLAST — writing results");

    let out_path = workdir.join(&params.output_name);
    valenx_core::io_caps::atomic_write_str(
        &out_path,
        &(output_lines.join("\n") + if output_lines.is_empty() { "" } else { "\n" }),
    )
    .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", out_path.display())))?;

    ctx.report_progress(100.0, "native BLAST — done");

    Ok(RunReport {
        exit_code: 0,
        wall_time: start.elapsed(),
        converged: Some(true),
        residual_history: Vec::new(),
        warnings: Vec::new(),
        final_phase: Some(RunPhase::Shutdown),
    })
}

/// Detects whether sequences look like nucleotide or protein.
fn detect_scheme(records: &[(Option<String>, Vec<u8>)]) -> ScoringScheme {
    let dna_set: std::collections::HashSet<u8> = b"ACGTURYNSWKMBVHD".iter().copied().collect();

    let mut total = 0usize;
    let mut dna_count = 0usize;
    for (_, seq) in records {
        for &b in seq.iter().take(500) {
            let up = b.to_ascii_uppercase();
            if up != b'-' && up != b'.' {
                total += 1;
                if dna_set.contains(&up) {
                    dna_count += 1;
                }
            }
        }
    }

    if total == 0 || dna_count * 10 >= total * 8 {
        ScoringScheme::dna_default()
    } else {
        ScoringScheme::blosum62_default()
    }
}

/// Reads a file to a String, capped at 256 MiB.
fn read_text_file(path: &str) -> Result<String, AdapterError> {
    const MAX_BYTES: u64 = 256 * 1024 * 1024;
    let mut f = std::fs::File::open(path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("open {path}: {e}")))?;
    let mut buf = Vec::with_capacity(4096);
    Read::by_ref(&mut f)
        .take(MAX_BYTES + 1)
        .read_to_end(&mut buf)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {path}: {e}")))?;
    if buf.len() as u64 > MAX_BYTES {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{path}: input file exceeds 256 MiB"
        )));
    }
    String::from_utf8(buf)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("{path}: not valid UTF-8: {e}")))
}

/// Minimal multi-FASTA parser: (name, sequence_bytes) pairs.
pub fn parse_fasta_simple(text: &str) -> Result<Vec<(Option<String>, Vec<u8>)>, AdapterError> {
    let mut records: Vec<(Option<String>, Vec<u8>)> = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    let mut in_fasta = false;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(name) = line.strip_prefix('>') {
            if in_fasta || !cur_seq.is_empty() {
                records.push((cur_name.take(), std::mem::take(&mut cur_seq)));
            }
            cur_name = Some(name.trim().to_string());
            in_fasta = true;
        } else {
            cur_seq.extend(
                line.bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    if !cur_seq.is_empty() || cur_name.is_some() {
        records.push((cur_name, cur_seq));
    }
    if records.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "input contains no sequences"
        )));
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    const PROTEIN_SEQ: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

    #[test]
    fn detect_scheme_dna() {
        let records = vec![
            (Some("a".to_string()), b"ACGTACGTACGT".to_vec()),
            (Some("b".to_string()), b"TTGGCCAATTGG".to_vec()),
        ];
        let scheme = detect_scheme(&records);
        assert!(
            scheme.matrix.name().contains("NUC") || scheme.matrix.name().contains("dna"),
            "expected DNA scheme, got: {}",
            scheme.matrix.name()
        );
    }

    #[test]
    fn detect_scheme_protein() {
        let records = vec![
            (Some("a".to_string()), PROTEIN_SEQ.to_vec()),
            (Some("b".to_string()), b"WYVTSRQPNMLKIHGFEDCA".to_vec()),
        ];
        let scheme = detect_scheme(&records);
        assert!(
            scheme.matrix.name().contains("BLOSUM"),
            "expected protein scheme, got: {}",
            scheme.matrix.name()
        );
    }

    #[test]
    fn self_search_finds_exact_match() {
        // A protein sequence searched against a 2-seq database containing
        // itself should produce a hit at 100% identity with low E-value.
        let query_fasta = ">query1\nACDEFGHIKLMNPQRSTVWY\n";
        let db_fasta = ">seq1\nACDEFGHIKLMNPQRSTVWY\n>seq2\nMMMMMMMMMMMMMMMMMMMM\n";

        let query_records = parse_fasta_simple(query_fasta).unwrap();
        let db_records = parse_fasta_simple(db_fasta).unwrap();
        let scheme = detect_scheme(&query_records);
        let db_seqs_owned: Vec<Vec<u8>> = db_records.iter().map(|(_, s)| s.clone()).collect();
        let db_seqs: Vec<&[u8]> = db_seqs_owned.iter().map(|s| s.as_slice()).collect();

        let index = KmerIndex::build_many(&db_seqs, 5).unwrap();
        let searcher = SeedSearch::new(&index, db_seqs.clone(), &scheme, SeedParams::default())
            .with_stats(KarlinAltschul::blosum62_ungapped());

        let hsps = searcher.search(PROTEIN_SEQ);
        assert!(!hsps.is_empty(), "self-search must find at least one hit");

        // The top hit should have seq_id 0 (the exact match).
        let top = &hsps[0];
        assert_eq!(
            top.seq_id, 0,
            "top hit should be against the identical sequence"
        );
        assert!(
            top.e_value.is_some_and(|e| e < 1.0),
            "self-search E-value should be < 1.0; got {:?}",
            top.e_value
        );
    }

    #[test]
    fn find_db_fasta_finds_companion_extension() {
        let d = tempdir("blast-native");
        let fa_path = d.join("mydb.fa");
        std::fs::write(&fa_path, ">s\nACGT\n").unwrap();

        // Look for the prefix without extension.
        let prefix = d.join("mydb");
        let found = find_db_fasta(&prefix);
        assert_eq!(found.as_deref(), Some(fa_path.as_path()));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn write_and_read_params_round_trips() {
        let d = tempdir("blast-native");
        let params = NativeBlastParams {
            query_path: "/tmp/query.fa".to_string(),
            db_fasta_path: "/tmp/db.fa".to_string(),
            evalue: 1e-5,
            max_hits: 500,
            output_name: "blast_results.txt".to_string(),
        };
        write_params(&d, &params).unwrap();
        let round = read_params(&d).unwrap();
        assert_eq!(round.query_path, params.query_path);
        assert_eq!(round.db_fasta_path, params.db_fasta_path);
        assert!((round.evalue - params.evalue).abs() < 1e-15);
        assert_eq!(round.max_hits, params.max_hits);
        assert_eq!(round.output_name, params.output_name);
        let _ = std::fs::remove_dir_all(&d);
    }
}
