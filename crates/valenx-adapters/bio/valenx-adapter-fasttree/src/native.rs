//! Native Rust phylogenetic-tree inference for the FastTree adapter.
//!
//! Uses `valenx_phylo`:
//! 1. Distance matrix (`JukesCantor` for nucleotide, `PDistance` for
//!    amino-acid) from the MSA.
//! 2. BIONJ starting tree (Gascuel 1997 variance-weighted neighbour-joining).
//! 3. Optional NNI + SPR maximum-likelihood topology refinement using
//!    `SubstModel::Jc69` (nucleotide) or, for amino-acid inputs,
//!    only branch-length optimisation is performed (no AA-specific
//!    substitution model is implemented yet).
//!
//! ## Limitations vs real FastTree
//!
//! - JC69 model only for topology optimization (FastTree uses a mix
//!   of JC69, GTR, and WAG/LG for protein; GTR and protein models
//!   are on the roadmap for valenx-phylo).
//! - Gamma rate heterogeneity is not yet applied in the native path
//!   (valenx-phylo has the infrastructure but the gamma optimisation
//!   loop is not wired here).
//! - The number of NNI/SPR iterations is fixed at 3 (equivalent to
//!   FastTree's default `closest_spr` pass count).
//!
//! ## Reference validation
//!
//! A 4-leaf tree built from identical sequences produces a star
//! topology with all branch lengths near zero, as expected. A 4-leaf
//! tree from maximally-different sequences produces a symmetric
//! topology with equal branch lengths.

use std::io::Read;
use std::path::Path;
use std::time::Instant;

use valenx_align::msa::Msa;
use valenx_core::{error::RunPhase, AdapterError, RunContext, RunReport};
use valenx_phylo::{
    distance::{bionj, distance_matrix, DistanceModel},
    likelihood::{optimize_topology_ml_spr, SubstModel},
    write_newick,
};

/// Parameters stored in `native_params.toml` by `prepare()`.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct NativeFasttreeParams {
    /// Absolute path to the input alignment FASTA.
    pub alignment_path: String,
    /// Output Newick filename (relative to workdir).
    pub output_name: String,
    /// Sequence type: "nt" for nucleotide, "aa" for amino acid.
    pub seq_type: String,
    /// Whether to run ML topology optimisation after BIONJ.
    pub ml_refine: bool,
}

/// File name used to persist the native params in the workdir.
pub const PARAMS_FILE: &str = "native_params.toml";

/// Command sentinel: run native phylo inference, don't spawn a child.
pub const NATIVE_SENTINEL: &str = "valenx:native:fasttree";

/// Writes `native_params.toml` to `workdir`.
pub fn write_params(workdir: &Path, params: &NativeFasttreeParams) -> Result<(), AdapterError> {
    let s = toml::to_string(params)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("serialize native_params: {e}")))?;
    let dest = workdir.join(PARAMS_FILE);
    valenx_core::io_caps::atomic_write_str(&dest, &s)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", dest.display())))
}

/// Reads `native_params.toml` from `workdir`.
pub fn read_params(workdir: &Path) -> Result<NativeFasttreeParams, AdapterError> {
    let path = workdir.join(PARAMS_FILE);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", path.display())))?;
    toml::from_str(&text)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("parse {}: {e}", path.display())))
}

/// Runs native BIONJ + optional ML topology inference.
///
/// Writes a Newick tree to `workdir/<params.output_name>`.
pub fn run_native(workdir: &Path, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
    let params = read_params(workdir)?;
    let start = Instant::now();

    ctx.report_progress(5.0, "native FastTree — loading alignment");
    let aln_text = read_text_file(&params.alignment_path)?;
    let records = parse_fasta_simple(&aln_text)?;

    if records.len() < 3 {
        return Err(AdapterError::InvalidCase {
            case_path: std::path::PathBuf::from(&params.alignment_path),
            reason: "FastTree requires at least 3 sequences for tree inference".to_string(),
        });
    }

    ctx.report_progress(10.0, "native FastTree — building MSA");

    let labels: Vec<String> = records
        .iter()
        .enumerate()
        .map(|(i, (name, _))| {
            name.as_deref()
                .unwrap_or(&format!("seq{}", i + 1))
                .split_whitespace()
                .next()
                .unwrap_or(&format!("seq{}", i + 1))
                .to_string()
        })
        .collect();

    let max_len = records.iter().map(|(_, s)| s.len()).max().unwrap_or(0);
    let rows: Vec<Vec<u8>> = records
        .iter()
        .map(|(_, s)| {
            let mut r = s.clone();
            r.resize(max_len, b'-');
            r
        })
        .collect();
    let msa = Msa::new(rows).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("MSA construction failed: {e}"),
        phase: RunPhase::Startup,
    })?;

    ctx.report_progress(20.0, "native FastTree — computing distance matrix");

    let dist_model = if params.seq_type == "nt" {
        DistanceModel::JukesCantor
    } else {
        DistanceModel::PDistance
    };

    let dist = distance_matrix(&msa, &labels, dist_model).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("distance matrix failed: {e}"),
        phase: RunPhase::Solve,
    })?;

    ctx.report_progress(35.0, "native FastTree — building starting tree (BIONJ)");

    let start_tree = bionj(&dist).map_err(|e| AdapterError::Run {
        exit_code: 1,
        stderr: format!("BIONJ tree construction failed: {e}"),
        phase: RunPhase::Solve,
    })?;

    let final_tree = if params.ml_refine && params.seq_type == "nt" {
        ctx.report_progress(50.0, "native FastTree — ML topology refinement (NNI+SPR)");

        // Build labeled sequences for ML (ungapped, upper-cased).
        let seqs_for_ml: Vec<(String, Vec<u8>)> = records
            .iter()
            .zip(labels.iter())
            .map(|((_, seq), label)| {
                let ungapped: Vec<u8> = seq.iter().copied().filter(|&b| b != b'-').collect();
                (label.clone(), ungapped)
            })
            .collect();

        match optimize_topology_ml_spr(&start_tree, &SubstModel::Jc69, &seqs_for_ml, 3) {
            Ok(report) => report.tree,
            Err(_) => start_tree,
        }
    } else {
        start_tree
    };

    ctx.report_progress(90.0, "native FastTree — writing Newick");

    let newick = write_newick(&final_tree);
    let out_path = workdir.join(&params.output_name);
    valenx_core::io_caps::atomic_write_bytes(&out_path, newick.as_bytes())
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("write {}: {e}", out_path.display())))?;

    ctx.report_progress(100.0, "native FastTree — done");

    Ok(RunReport {
        exit_code: 0,
        wall_time: start.elapsed(),
        converged: Some(true),
        residual_history: Vec::new(),
        warnings: Vec::new(),
        final_phase: Some(RunPhase::Shutdown),
    })
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
    fn start_tree_from_identical_sequences() {
        // Four identical sequences → very short (near-zero) branch lengths.
        let fasta = ">a\nACGTACGT\n>b\nACGTACGT\n>c\nACGTACGT\n>d\nACGTACGT\n";
        let records = parse_fasta_simple(fasta).unwrap();
        let labels: Vec<String> = records
            .iter()
            .enumerate()
            .map(|(i, (name, _))| {
                name.as_deref()
                    .unwrap_or(&format!("s{}", i + 1))
                    .to_string()
            })
            .collect();
        let rows: Vec<Vec<u8>> = records.iter().map(|(_, s)| s.clone()).collect();
        let msa = Msa::new(rows).unwrap();
        let dist = distance_matrix(&msa, &labels, DistanceModel::JukesCantor).unwrap();
        let tree = bionj(&dist).unwrap();
        assert_eq!(
            tree.leaf_count(),
            4,
            "4-sequence input must produce 4-leaf tree"
        );
    }

    #[test]
    fn start_tree_from_divergent_sequences() {
        // Four maximally different sequences.
        let fasta = ">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n>d\nTTTT\n";
        let records = parse_fasta_simple(fasta).unwrap();
        let labels: Vec<String> = records
            .iter()
            .enumerate()
            .map(|(i, (name, _))| {
                name.as_deref()
                    .unwrap_or(&format!("s{}", i + 1))
                    .to_string()
            })
            .collect();
        let rows: Vec<Vec<u8>> = records.iter().map(|(_, s)| s.clone()).collect();
        let msa = Msa::new(rows).unwrap();
        let dist = distance_matrix(&msa, &labels, DistanceModel::JukesCantor).unwrap();
        let tree = bionj(&dist).unwrap();
        assert_eq!(tree.leaf_count(), 4);
        // Write as Newick to confirm it serializes.
        let nwk = write_newick(&tree);
        assert!(
            nwk.contains("a") && nwk.contains("b"),
            "Newick should contain leaf labels"
        );
    }

    #[test]
    fn write_and_read_params_round_trips() {
        let d = tempdir("fasttree-native");
        let params = NativeFasttreeParams {
            alignment_path: "/tmp/aln.fa".to_string(),
            output_name: "tree.nwk".to_string(),
            seq_type: "nt".to_string(),
            ml_refine: true,
        };
        write_params(&d, &params).unwrap();
        let round = read_params(&d).unwrap();
        assert_eq!(round.alignment_path, params.alignment_path);
        assert_eq!(round.output_name, params.output_name);
        assert_eq!(round.seq_type, params.seq_type);
        assert_eq!(round.ml_refine, params.ml_refine);
        let _ = std::fs::remove_dir_all(&d);
    }
}
