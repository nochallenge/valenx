//! Native Rust profile-HMM search for the HMMER adapter.
//!
//! Uses `valenx_align::hmm::ProfileHmm` built from a FASTA alignment,
//! then scores each target sequence with Viterbi log-probability.
//!
//! ## Limitations vs real HMMER
//!
//! - The `profile` field in native mode must point to a **FASTA multiple
//!   alignment** file (`.fa`, `.fasta`, `.aln`, `.sto`), NOT a pre-built
//!   `.hmm` database. HMMER's binary `.hmm` format is not parsed here.
//! - The score reported is the Viterbi log-probability (ln P), not a
//!   calibrated E-value. Real HMMER computes E-values by fitting an
//!   EVD/Gumbel to per-database random scores; this v1 uses a fixed
//!   score threshold instead.
//! - Only `hmmsearch`-style operation (profile-vs-sequences) is
//!   implemented. `hmmscan` (sequences-vs-profile-database) is
//!   equivalent for a single profile and is handled identically.
//! - The profile HMM uses fixed Plan7 transition probabilities and
//!   emission probabilities estimated from the alignment columns
//!   (Laplace-smoothed observed frequencies — NOT a BLOSUM matrix).
//! - **Protein profiles only.** The emission alphabet is the 20 amino
//!   acids; nucleotide symbols (A/C/G/T) would collide on that index, so
//!   native mode rejects DNA/RNA alignments. Use the HMMER binary
//!   `nhmmer` for nucleotide profiles.
//!
//! ## Reference validation
//!
//! A 20-aa protein profile built from a 3-sequence ACDEFGHIKLMNPQRSTVWY
//! alignment scores the member sequence well above an unrelated all-M
//! sequence, so the profile discriminates family membership. Scores are
//! raw Viterbi log-probabilities (ln P), not length-normalized log-odds,
//! so they are not directly comparable across sequences of different
//! lengths — a calibrated E-value is future work.

use std::io::Read;
use std::path::Path;
use std::time::Instant;

use valenx_align::{hmm::ProfileHmm, msa::Msa};
use valenx_core::{error::RunPhase, AdapterError, RunContext, RunReport};

/// Parameters stored in `native_params.toml` by `prepare()`.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct NativeHmmerParams {
    /// Absolute path to the profile FASTA alignment (native-mode input).
    pub profile_fasta_path: String,
    /// Absolute path to the target sequences FASTA file.
    pub sequences_path: String,
    /// Minimum Viterbi log-probability to report (0.0 = only exact
    /// matches; -200.0 = report most homologs; default -100.0).
    pub min_score: f64,
    /// Output filenames (relative to workdir).
    pub tblout_name: String,
    pub report_name: String,
}

/// File name used to persist the native params in the workdir.
pub const PARAMS_FILE: &str = "native_params.toml";

/// Command sentinel: run native profile search, don't spawn a child.
pub const NATIVE_SENTINEL: &str = "valenx:native:hmmer";

/// Writes `native_params.toml` to `workdir`.
pub fn write_params(workdir: &Path, params: &NativeHmmerParams) -> Result<(), AdapterError> {
    let s = toml::to_string(params)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("serialize native_params: {e}")))?;
    let dest = workdir.join(PARAMS_FILE);
    valenx_core::io_caps::atomic_write_str(&dest, &s)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", dest.display())))
}

/// Reads `native_params.toml` from `workdir`.
pub fn read_params(workdir: &Path) -> Result<NativeHmmerParams, AdapterError> {
    let path = workdir.join(PARAMS_FILE);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", path.display())))?;
    toml::from_str(&text)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("parse {}: {e}", path.display())))
}

/// Returns true if `profile_path` looks like a FASTA alignment file
/// (by extension) rather than a prebuilt `.hmm` database.
pub fn profile_is_fasta(path: &Path) -> bool {
    match path.extension().and_then(|s| s.to_str()) {
        Some(e) => matches!(
            e.to_ascii_lowercase().as_str(),
            "fa" | "fasta" | "aln" | "sto" | "mfa"
        ),
        None => false,
    }
}

/// Runs native profile-HMM search.
///
/// Writes a simplified tblout to `workdir/<params.tblout_name>` and a
/// summary report to `workdir/<params.report_name>`.
pub fn run_native(workdir: &Path, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
    let params = read_params(workdir)?;
    let start = Instant::now();

    ctx.report_progress(5.0, "native HMMER — loading profile alignment");
    let profile_text = read_text_file(&params.profile_fasta_path)?;
    let profile_records = parse_fasta_simple(&profile_text)?;

    if profile_records.len() < 2 {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.profile_fasta_path),
            reason: "profile FASTA must contain at least 2 aligned sequences".to_string(),
        });
    }

    // GATE: the native profile HMM uses a 20-letter amino-acid emission
    // alphabet, so nucleotide symbols (A/C/G/T) collide on the index (T
    // aliases A) and a DNA/RNA profile would score meaningless values.
    // Refuse nucleotide profiles in native mode rather than emit garbage;
    // the real HMMER `nhmmer` binary handles nucleotide profiles.
    // TODO(bio): add a 4-symbol nucleotide emission alphabet, then lift.
    if looks_like_nucleotide(&profile_records) {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.profile_fasta_path),
            reason: "native HMMER supports PROTEIN profiles only: its profile-HMM uses \
                     a 20-letter amino-acid emission alphabet, so nucleotide profiles \
                     score incorrectly. Install the HMMER binary (nhmmer) for nucleotide \
                     profiles, or supply a protein alignment."
                .to_string(),
        });
    }

    ctx.report_progress(15.0, "native HMMER — building profile HMM");

    // Build an Msa from the aligned FASTA rows.
    let rows: Vec<Vec<u8>> = profile_records.iter().map(|(_, s)| s.clone()).collect();
    let max_len = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    // Pad rows to equal length (gap-pad to the longest row).
    let rows_padded: Vec<Vec<u8>> = rows
        .into_iter()
        .map(|mut r| {
            r.resize(max_len, b'-');
            r
        })
        .collect();
    let msa = Msa::new(rows_padded).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("MSA construction failed: {e}"),
        phase: RunPhase::Startup,
    })?;

    let profile = ProfileHmm::from_msa(&msa, 0.5).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("profile-HMM build failed: {e}"),
        phase: RunPhase::Startup,
    })?;

    ctx.report_progress(20.0, "native HMMER — loading target sequences");
    let seqs_text = read_text_file(&params.sequences_path)?;
    let seq_records = parse_fasta_simple(&seqs_text)?;

    ctx.report_progress(30.0, "native HMMER — scoring sequences");

    let profile_name = std::path::Path::new(&params.profile_fasta_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("profile")
        .to_string();

    let mut hits: Vec<(String, f64)> = Vec::new();
    let total = seq_records.len();

    for (qi, (name, seq)) in seq_records.iter().enumerate() {
        if qi % 50 == 0 {
            let pct = 30.0 + (qi as f32 / total as f32) * 60.0;
            ctx.report_progress(
                pct,
                &format!("native HMMER — {}/{} sequences scored", qi + 1, total),
            );
        }

        let seqid = name
            .as_deref()
            .unwrap_or("seq")
            .split_whitespace()
            .next()
            .unwrap_or("seq")
            .to_string();

        // Strip gaps from the target sequence before scoring.
        let ungapped: Vec<u8> = seq.iter().copied().filter(|&b| b != b'-').collect();
        if ungapped.is_empty() {
            continue;
        }

        match profile.viterbi(&ungapped) {
            Ok(score) if score >= params.min_score => {
                hits.push((seqid, score));
            }
            Ok(_) => {}
            Err(_) => {}
        }
    }

    // Sort by score descending (best hit first).
    hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    ctx.report_progress(95.0, "native HMMER — writing results");

    // Write simplified tblout: target_name, query_name, score, description
    let tblout_lines: Vec<String> = {
        let mut lines = vec![
            "# valenx native HMMER (profile-HMM Viterbi log-probabilities)".to_string(),
            "# Fields: target_name  profile_name  viterbi_logprob  description".to_string(),
        ];
        for (name, score) in &hits {
            lines.push(format!("{name}\t{profile_name}\t{score:.4}\t-"));
        }
        lines
    };
    let tblout_path = workdir.join(&params.tblout_name);
    valenx_core::io_caps::atomic_write_str(&tblout_path, &(tblout_lines.join("\n") + "\n"))
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("write {}: {e}", tblout_path.display()))
        })?;

    // Write summary report.
    let report = format!(
        "# valenx native HMMER report\n\
         # Profile: {profile_name} ({} match columns)\n\
         # Target:  {} sequences scanned\n\
         # Hits:    {} sequences above min_score {}\n\n",
        profile.length(),
        total,
        hits.len(),
        params.min_score,
    );
    let report_path = workdir.join(&params.report_name);
    valenx_core::io_caps::atomic_write_str(&report_path, &report).map_err(|e| {
        AdapterError::Other(anyhow::anyhow!("write {}: {e}", report_path.display()))
    })?;

    ctx.report_progress(100.0, "native HMMER — done");

    Ok(RunReport {
        exit_code: 0,
        wall_time: start.elapsed(),
        converged: Some(true),
        residual_history: Vec::new(),
        warnings: Vec::new(),
        final_phase: Some(RunPhase::Shutdown),
    })
}

/// Heuristic: does this alignment look like nucleotide (DNA/RNA) rather
/// than protein? Samples up to 500 residues per row and returns true when
/// >=80% of non-gap residues are in the IUPAC nucleotide set. Mirrors the
/// > BLAST adapter's alphabet detection.
fn looks_like_nucleotide(records: &[(Option<String>, Vec<u8>)]) -> bool {
    let nuc: std::collections::HashSet<u8> = b"ACGTURYNSWKMBVHD".iter().copied().collect();
    let mut total = 0usize;
    let mut nuc_count = 0usize;
    for (_, seq) in records {
        for &b in seq.iter().take(500) {
            let up = b.to_ascii_uppercase();
            if up != b'-' && up != b'.' {
                total += 1;
                if nuc.contains(&up) {
                    nuc_count += 1;
                }
            }
        }
    }
    total > 0 && nuc_count * 10 >= total * 8
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
    fn profile_built_from_alignment_scores_member_higher_than_random() {
        // Build a profile from 3 copies of a 20-aa protein + minor variant.
        let profile_fasta = ">s1\nACDEFGHIKLMNPQRSTVWY\n\
                             >s2\nACDEFGHIKLMNPQRSTVWY\n\
                             >s3\nACDEFGHIKLMNPQRSTVWY\n";
        let profile_records = parse_fasta_simple(profile_fasta).unwrap();
        let rows: Vec<Vec<u8>> = profile_records.iter().map(|(_, s)| s.clone()).collect();
        let msa = Msa::new(rows).unwrap();
        let profile = ProfileHmm::from_msa(&msa, 0.5).unwrap();

        // Score the member sequence (should score high / less negative).
        let member_score = profile.viterbi(PROTEIN_SEQ).unwrap();

        // Score a completely different sequence (all Methionine).
        let random_seq: Vec<u8> = b"MMMMMMMMMMMMMMMMMMMM".to_vec();
        let random_score = profile.viterbi(&random_seq).unwrap();

        assert!(
            member_score > random_score,
            "member score ({member_score:.2}) should exceed random score ({random_score:.2})"
        );
    }

    #[test]
    fn nucleotide_profile_is_detected_and_protein_is_not() {
        let dna = parse_fasta_simple(">a\nACGTACGTACGT\n>b\nACGTACGTACGT\n").unwrap();
        assert!(
            looks_like_nucleotide(&dna),
            "DNA alignment should be detected"
        );
        let protein = parse_fasta_simple(">a\nEFILPQEFILPQ\n>b\nEFILPQEFILPQ\n").unwrap();
        assert!(
            !looks_like_nucleotide(&protein),
            "protein alignment should not be detected as nucleotide"
        );
    }

    #[test]
    fn profile_is_fasta_detects_extensions() {
        assert!(profile_is_fasta(Path::new("pf.fa")));
        assert!(profile_is_fasta(Path::new("pf.fasta")));
        assert!(profile_is_fasta(Path::new("pf.aln")));
        assert!(!profile_is_fasta(Path::new("pf.hmm")));
        assert!(!profile_is_fasta(Path::new("pf")));
    }

    #[test]
    fn write_and_read_params_round_trips() {
        let d = tempdir("hmmer-native");
        let params = NativeHmmerParams {
            profile_fasta_path: "/tmp/profile.fa".to_string(),
            sequences_path: "/tmp/seqs.fa".to_string(),
            min_score: -100.0,
            tblout_name: "tblout.txt".to_string(),
            report_name: "hmmer.out".to_string(),
        };
        write_params(&d, &params).unwrap();
        let round = read_params(&d).unwrap();
        assert_eq!(round.profile_fasta_path, params.profile_fasta_path);
        assert_eq!(round.sequences_path, params.sequences_path);
        assert!((round.min_score - params.min_score).abs() < 1e-10);
        assert_eq!(round.tblout_name, params.tblout_name);
        assert_eq!(round.report_name, params.report_name);
        let _ = std::fs::remove_dir_all(&d);
    }
}
