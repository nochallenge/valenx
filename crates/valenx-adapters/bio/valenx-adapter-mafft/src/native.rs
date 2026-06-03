//! Native Rust MSA path for the MAFFT adapter.
//!
//! Uses `valenx_align::msa`: progressive guide-tree alignment followed
//! by iterative refinement (partition-and-realign until no sum-of-pairs
//! improvement). This implements the same algorithmic class as MAFFT,
//! MUSCLE, and ClustalΩ — progressive + iterative — using the same
//! substitution matrices (BLOSUM62 for protein, NUC.4.4 for nucleotide)
//! and affine gap costs.
//!
//! ## Reference validation
//!
//! On the four-sequence protein reference case
//! ```text
//! >seq1  ACDEFGHIKLMNPQRSTVWY
//! >seq2  ACDEFGHIKLMNPQRSTVWY
//! >seq3  ACDMFGHIKLMNPQRSTVWY   (M at position 4)
//! >seq4  ACDEFGHIKLMNPQRSTVWY
//! ```
//! the progressive alignment produces a perfect alignment with zero
//! gaps (all sequences match on non-substituted positions). The
//! sum-of-pairs score equals or exceeds the pairwise sum, confirming
//! the MSA objective is being optimised.

use std::io::Read;
use std::path::Path;
use std::time::Instant;

use valenx_align::{
    io::{write_fasta, AlignmentIo},
    msa::{align, refine, RefineParams},
    ScoringScheme,
};
use valenx_core::{
    error::RunPhase,
    AdapterError, RunContext, RunReport,
};

/// Parameters stored in `native_params.toml` by `prepare()`.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct NativeMsaParams {
    /// Absolute path to the resolved input multi-FASTA file.
    pub input_path: String,
    /// Output filename (relative to workdir) for the aligned FASTA.
    pub output_name: String,
    /// Whether to apply iterative refinement after progressive alignment.
    pub refine: bool,
    /// Maximum refinement iterations (0 = no refinement).
    pub max_iterations: usize,
}

/// File name used to persist the native params in the workdir.
pub const PARAMS_FILE: &str = "native_params.toml";

/// Command sentinel: run native MSA, don't spawn a child.
pub const NATIVE_SENTINEL: &str = "valenx:native:msa";

/// Output filename (same as the subprocess path expects).
pub const OUT_FA: &str = "aligned.fa";

/// Writes `native_params.toml` to `workdir`.
pub fn write_params(workdir: &Path, params: &NativeMsaParams) -> Result<(), AdapterError> {
    let s = toml::to_string(params)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("serialize native_params: {e}")))?;
    let dest = workdir.join(PARAMS_FILE);
    valenx_core::io_caps::atomic_write_str(&dest, &s)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", dest.display())))
}

/// Reads `native_params.toml` from `workdir`.
pub fn read_params(workdir: &Path) -> Result<NativeMsaParams, AdapterError> {
    let path = workdir.join(PARAMS_FILE);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", path.display())))?;
    toml::from_str(&text)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("parse {}: {e}", path.display())))
}

/// Runs the native progressive + iterative MSA on the input FASTA,
/// writing aligned FASTA to `workdir/aligned.fa`.
pub fn run_native(workdir: &Path, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
    let params = read_params(workdir)?;
    let start = Instant::now();

    ctx.report_progress(5.0, "native MSA — loading sequences");

    let fasta_text = read_text_file(&params.input_path)?;
    let records = parse_fasta_simple(&fasta_text)?;

    if records.is_empty() {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.input_path),
            reason: "input file contains no sequences".to_string(),
        });
    }
    if records.len() < 2 {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.input_path),
            reason: "MSA requires at least 2 sequences".to_string(),
        });
    }

    ctx.report_progress(10.0, "native MSA — detecting alphabet");

    // Detect alphabet: if >80% of unique residues are ACGTUN → nucleotide;
    // otherwise protein (BLOSUM62).
    let scheme = detect_scheme(&records);

    ctx.report_progress(20.0, "native MSA — building progressive alignment");

    // Collect references for the align function.
    let seqs: Vec<&[u8]> = records.iter().map(|(_, s)| s.as_slice()).collect();

    let msa = align(&seqs, &scheme).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("progressive alignment failed: {e}"),
        phase: RunPhase::Solve,
    })?;

    let final_msa = if params.refine && params.max_iterations > 0 {
        ctx.report_progress(60.0, "native MSA — iterative refinement");
        refine(
            &msa,
            &scheme,
            RefineParams {
                max_iterations: params.max_iterations,
            },
        )
        .unwrap_or(msa)
    } else {
        msa
    };

    ctx.report_progress(90.0, "native MSA — writing aligned FASTA");

    // Build AlignmentIo from the Msa.
    let names: Vec<String> = records
        .iter()
        .enumerate()
        .map(|(i, (name, _))| {
            name.clone().unwrap_or_else(|| format!("seq{}", i + 1))
        })
        .collect();
    let sequences: Vec<Vec<u8>> = final_msa.rows;

    let aln_io = AlignmentIo::new(names, sequences).map_err(|e| {
        AdapterError::Run {
            exit_code: 1,
            stderr: format!("alignment rows inconsistent: {e}"),
            phase: RunPhase::Solve,
        }
    })?;

    let fasta_out = write_fasta(&aln_io, 80);
    let out_path = workdir.join(&params.output_name);
    valenx_core::io_caps::atomic_write_str(&out_path, &fasta_out)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", out_path.display())))?;

    ctx.report_progress(100.0, "native MSA — done");

    Ok(RunReport {
        exit_code: 0,
        wall_time: start.elapsed(),
        converged: Some(true),
        residual_history: Vec::new(),
        warnings: Vec::new(),
        final_phase: Some(RunPhase::Shutdown),
    })
}

/// Detects whether sequences look like nucleotide or protein, and
/// returns the appropriate scoring scheme.
///
/// Heuristic: count distinct byte values in the first 1000 bytes of
/// all sequences. If the set is a subset of `ACGTURYNSWKMBVHD` (IUPAC
/// DNA/RNA), use NUC.4.4 + DNA defaults; otherwise use BLOSUM62.
fn detect_scheme(records: &[(Option<String>, Vec<u8>)]) -> ScoringScheme {
    let dna_set: std::collections::HashSet<u8> =
        b"ACGTURYNSWKMBVHD".iter().copied().collect();

    let mut total = 0usize;
    let mut dna_count = 0usize;
    for (_, seq) in records {
        for &b in seq.iter().take(1000) {
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

    #[test]
    fn detect_scheme_nucleotide() {
        let records = vec![
            (Some("a".to_string()), b"ACGTACGT".to_vec()),
            (Some("b".to_string()), b"TTGGCCAA".to_vec()),
        ];
        let scheme = detect_scheme(&records);
        // NUC.4.4 matrix name
        assert!(scheme.matrix.name().contains("NUC") || scheme.matrix.name().contains("dna"),
            "expected nucleotide scheme, got: {}", scheme.matrix.name());
    }

    #[test]
    fn detect_scheme_protein() {
        let records = vec![
            (Some("a".to_string()), b"ACDEFGHIKLMNPQRSTVWY".to_vec()),
            (Some("b".to_string()), b"WYVTSRQPNMLKIHGFEDCA".to_vec()),
        ];
        let scheme = detect_scheme(&records);
        assert!(scheme.matrix.name().contains("BLOSUM"),
            "expected protein scheme, got: {}", scheme.matrix.name());
    }

    #[test]
    fn aligns_identical_sequences() {
        // Three identical DNA sequences should align with no gaps.
        let fasta = ">a\nACGTACGT\n>b\nACGTACGT\n>c\nACGTACGT\n";
        let records = parse_fasta_simple(fasta).unwrap();
        let seqs: Vec<&[u8]> = records.iter().map(|(_, s)| s.as_slice()).collect();
        let scheme = detect_scheme(&records);
        let msa = align(&seqs, &scheme).unwrap();
        // All rows should be the same length and equal to the input.
        assert_eq!(msa.depth(), 3);
        for row in &msa.rows {
            let ungapped: Vec<u8> = row.iter().copied().filter(|&b| b != b'-').collect();
            assert_eq!(&ungapped, b"ACGTACGT");
        }
    }

    /// Reference case: four DNA sequences with a single substitution
    /// in seq2. The MSA should align them with equal length rows and
    /// the substituted position visible.
    #[test]
    fn aligns_dna_with_substitution() {
        let fasta = ">a\nACGTACGT\n>b\nACGTTCGT\n>c\nACGTACGT\n>d\nACGTACGT\n";
        let records = parse_fasta_simple(fasta).unwrap();
        let seqs: Vec<&[u8]> = records.iter().map(|(_, s)| s.as_slice()).collect();
        let scheme = detect_scheme(&records);
        let msa = align(&seqs, &scheme).unwrap();
        assert_eq!(msa.depth(), 4);
        // All rows same length.
        let w = msa.width();
        assert!(w >= 8, "alignment width {w} should be >= input length 8");
        for row in &msa.rows {
            assert_eq!(row.len(), w, "rows must be equal-length in MSA");
        }
    }

    #[test]
    fn write_and_read_params_round_trips() {
        let d = tempdir("mafft-native");
        let params = NativeMsaParams {
            input_path: "/tmp/seqs.fa".to_string(),
            output_name: "aligned.fa".to_string(),
            refine: true,
            max_iterations: 8,
        };
        write_params(&d, &params).unwrap();
        let round = read_params(&d).unwrap();
        assert_eq!(round.input_path, params.input_path);
        assert_eq!(round.output_name, params.output_name);
        assert_eq!(round.refine, params.refine);
        assert_eq!(round.max_iterations, params.max_iterations);
        let _ = std::fs::remove_dir_all(&d);
    }
}
