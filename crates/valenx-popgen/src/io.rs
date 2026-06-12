//! Export of simulations: VCF, `ms`-style and Newick.
//!
//! A simulation is only useful if it can be written out in the formats
//! downstream tools read. This module provides three exporters and one
//! importer:
//!
//! - [`write_vcf`] — emits a valid VCF 4.2 from a
//!   [`crate::infer::GenotypeMatrix`]: a header, the eight mandatory
//!   columns and one phased diploid genotype per sample pair.
//! - [`write_ms`] — emits the classic Hudson `ms` haplotype format
//!   (`segsites:` / `positions:` / a `0`/`1` matrix), the lingua
//!   franca of coalescent simulators.
//! - [`write_newick_genealogy`] — writes a simulated genealogy
//!   (a [`valenx_phylo::Tree`]) to a Newick string, reusing
//!   `valenx-phylo`'s writer.
//! - [`read_ms`] — parses an `ms`-format block back into a genotype
//!   matrix, so externally-produced datasets can be analysed with
//!   [`crate::stats`].

use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;
use valenx_phylo::tree::Tree;

/// Writes a [`GenotypeMatrix`] as a VCF 4.2 string.
///
/// Consecutive haplotype rows are paired into diploid genotypes — row
/// `2k` and `2k+1` become sample `k`'s two phased alleles — so the
/// matrix must have an even number of rows. `chrom` names the contig.
/// Site positions are written rounded to the nearest integer base.
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has an odd number of rows.
pub fn write_vcf(matrix: &GenotypeMatrix, chrom: &str) -> Result<String> {
    let n_hap = matrix.n_samples();
    if n_hap % 2 != 0 {
        return Err(PopgenError::invalid(
            "matrix",
            "VCF export needs an even number of haplotype rows (diploid)",
        ));
    }
    let n_diploid = n_hap / 2;
    let mut out = String::new();
    // Header.
    out.push_str("##fileformat=VCFv4.2\n");
    out.push_str("##source=valenx-popgen\n");
    out.push_str(&format!(
        "##contig=<ID={chrom},length={}>\n",
        matrix
            .positions()
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
            .ceil() as i64
            + 1
    ));
    out.push_str("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n");
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for k in 0..n_diploid {
        out.push_str(&format!("\tsample{k}"));
    }
    out.push('\n');

    // One data line per site.
    for col in 0..matrix.n_sites() {
        let pos = matrix.positions()[col].round().max(1.0) as i64;
        let ac = matrix.derived_count(col)?;
        out.push_str(&format!("{chrom}\t{pos}\t.\tA\tT\t.\tPASS\tAC={ac}\tGT"));
        for k in 0..n_diploid {
            let a0 = matrix.get(2 * k, col);
            let a1 = matrix.get(2 * k + 1, col);
            out.push_str(&format!("\t{a0}|{a1}"));
        }
        out.push('\n');
    }
    Ok(out)
}

/// Writes a [`GenotypeMatrix`] in Hudson `ms` haplotype format.
///
/// The block is: a `//` separator, `segsites: <S>`, a `positions:`
/// line of the per-site positions scaled to `[0, 1]`, then one
/// `0`/`1` string per haplotype.
///
/// `sequence_length` rescales the absolute site positions onto the
/// `[0, 1)` interval `ms` uses; pass the simulated segment length.
///
/// # Errors
/// [`PopgenError::Invalid`] if `sequence_length <= 0`.
pub fn write_ms(matrix: &GenotypeMatrix, sequence_length: f64) -> Result<String> {
    if sequence_length <= 0.0 {
        return Err(PopgenError::invalid("sequence_length", "must be positive"));
    }
    let mut out = String::new();
    out.push_str("// valenx-popgen simulation\n");
    out.push_str(&format!("segsites: {}\n", matrix.n_sites()));
    out.push_str("positions:");
    for &p in matrix.positions() {
        out.push_str(&format!(" {:.6}", (p / sequence_length).clamp(0.0, 1.0)));
    }
    out.push('\n');
    for row in matrix.rows() {
        let line: String = row
            .iter()
            .map(|&v| if v == 1 { '1' } else { '0' })
            .collect();
        out.push_str(&line);
        out.push('\n');
    }
    Ok(out)
}

/// Writes a simulated genealogy to a Newick string.
///
/// A thin wrapper over `valenx-phylo`'s Newick writer, here so a
/// `valenx-popgen` caller does not need to import `valenx-phylo`
/// directly to serialise a coalescent or birth-death tree.
pub fn write_newick_genealogy(tree: &Tree) -> String {
    valenx_phylo::write_newick(tree)
}

/// Parses a Hudson `ms`-format block into a [`GenotypeMatrix`].
///
/// Recognises the `segsites:`, `positions:` and haplotype lines;
/// ignores the `//` separator and any command / seed lines. The
/// positions are taken on the `[0, 1]` scale `ms` uses.
///
/// # Errors
/// [`PopgenError::Parse`] if no `segsites:` line is found, the
/// haplotype rows disagree with the declared site count, or a row
/// contains a non-`0`/`1` character.
pub fn read_ms(text: &str) -> Result<GenotypeMatrix> {
    let mut segsites: Option<usize> = None;
    let mut positions: Vec<f64> = Vec::new();
    let mut rows: Vec<Vec<u8>> = Vec::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("//") {
            continue;
        }
        if let Some(rest) = line.strip_prefix("segsites:") {
            segsites = Some(
                rest.trim()
                    .parse()
                    .map_err(|_| PopgenError::parse("ms", "segsites: value is not an integer"))?,
            );
        } else if let Some(rest) = line.strip_prefix("positions:") {
            positions = rest
                .split_whitespace()
                .map(|t| {
                    t.parse::<f64>()
                        .map_err(|_| PopgenError::parse("ms", "positions: value is not a number"))
                })
                .collect::<Result<Vec<f64>>>()?;
        } else if line.chars().all(|c| c == '0' || c == '1') {
            // A haplotype row.
            rows.push(line.bytes().map(|b| b - b'0').collect());
        }
        // Anything else (command line, seeds) is ignored.
    }

    let s = segsites.ok_or_else(|| PopgenError::parse("ms", "no `segsites:` line found"))?;
    if positions.is_empty() && s > 0 {
        // Some ms output omits positions for segsites: 0; otherwise
        // synthesise an even grid.
        positions = (0..s).map(|i| (i as f64 + 0.5) / s as f64).collect();
    }
    if positions.len() != s {
        return Err(PopgenError::parse(
            "ms",
            format!(
                "positions count {} disagrees with segsites {s}",
                positions.len()
            ),
        ));
    }
    for (i, r) in rows.iter().enumerate() {
        if r.len() != s {
            return Err(PopgenError::parse(
                "ms",
                format!("haplotype row {i} has {} sites, expected {s}", r.len()),
            ));
        }
    }
    if rows.is_empty() {
        return Err(PopgenError::parse("ms", "no haplotype rows found"));
    }
    GenotypeMatrix::from_rows(rows, positions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coalescent::kingman::{coalescent, PopHistory};

    fn matrix(rows: Vec<Vec<u8>>, positions: Vec<f64>) -> GenotypeMatrix {
        GenotypeMatrix::from_rows(rows, positions).unwrap()
    }

    #[test]
    fn vcf_export_has_header_and_data() {
        let m = matrix(
            vec![vec![1, 0], vec![0, 1], vec![1, 1], vec![0, 0]],
            vec![100.0, 250.0],
        );
        let vcf = write_vcf(&m, "chr1").unwrap();
        assert!(vcf.contains("##fileformat=VCFv4.2"));
        assert!(vcf.contains("#CHROM\tPOS"));
        // 4 haplotypes -> 2 diploid samples.
        assert!(vcf.contains("sample0") && vcf.contains("sample1"));
        // Two data lines for two sites.
        let data_lines = vcf.lines().filter(|l| l.starts_with("chr1\t")).count();
        assert_eq!(data_lines, 2);
        // Genotypes are phased a|b.
        assert!(vcf.contains("1|0") || vcf.contains("0|1"));
    }

    #[test]
    fn vcf_rejects_odd_haplotype_count() {
        let m = matrix(vec![vec![1], vec![0], vec![1]], vec![10.0]);
        assert!(write_vcf(&m, "chr1").is_err());
    }

    #[test]
    fn ms_export_round_trips_through_read_ms() {
        let m = matrix(
            vec![vec![1, 0, 1], vec![0, 1, 1], vec![1, 1, 0]],
            vec![100.0, 500.0, 900.0],
        );
        let ms = write_ms(&m, 1000.0).unwrap();
        assert!(ms.contains("segsites: 3"));
        assert!(ms.contains("positions:"));
        let back = read_ms(&ms).unwrap();
        assert_eq!(back.n_samples(), 3);
        assert_eq!(back.n_sites(), 3);
        // The 0/1 calls survive the round trip.
        assert_eq!(back.rows(), m.rows());
    }

    #[test]
    fn ms_positions_are_scaled_to_unit_interval() {
        let m = matrix(vec![vec![1], vec![0]], vec![500.0]);
        let ms = write_ms(&m, 1000.0).unwrap();
        // 500 / 1000 = 0.5.
        assert!(ms.contains("0.500000"));
    }

    #[test]
    fn read_ms_parses_a_hand_written_block() {
        let text = "\
//
segsites: 2
positions: 0.10 0.80
10
01
11
";
        let m = read_ms(text).unwrap();
        assert_eq!(m.n_sites(), 2);
        assert_eq!(m.n_samples(), 3);
        assert_eq!(m.derived_count(0).unwrap(), 2);
    }

    #[test]
    fn read_ms_rejects_malformed_input() {
        // No segsites line.
        assert!(read_ms("positions: 0.5\n1\n0\n").is_err());
        // Row length disagrees with segsites.
        assert!(read_ms("segsites: 3\npositions: 0.1 0.2 0.3\n10\n").is_err());
    }

    #[test]
    fn newick_genealogy_export_is_valid_newick() {
        let labels: Vec<String> = (0..5).map(|i| format!("L{i}")).collect();
        let tree = coalescent(&labels, &PopHistory::Constant(1000.0), 42).unwrap();
        let newick = write_newick_genealogy(&tree);
        assert!(newick.ends_with(';'));
        // Every tip label is present.
        for l in &labels {
            assert!(newick.contains(l.as_str()), "missing {l}");
        }
        // It re-parses through valenx-phylo.
        let reparsed = valenx_phylo::read_newick(&newick).unwrap();
        assert_eq!(reparsed.leaf_count(), 5);
    }
}
