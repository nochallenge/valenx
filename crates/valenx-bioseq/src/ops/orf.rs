//! Open reading frame (ORF) finding across all six frames.
//!
//! An [`Orf`] is a run of codons from a start codon to the next
//! in-frame stop codon. [`find_orfs`] scans all three forward frames
//! plus all three reverse frames; [`OrfOptions`] controls whether only
//! `ATG` counts as a start (vs. any table start codon) and the minimum
//! protein length.

use crate::error::{BioseqError, Result};
use crate::ops::revcomp::reverse_complement;
use crate::ops::translate::{translate_default, GeneticCode};
use crate::record::{Span, Strand};
use crate::seq::{Seq, SeqKind};

/// A single open reading frame found on a nucleotide sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Orf {
    /// Frame index `0..6` (`0..3` forward, `3..6` reverse).
    pub frame: u8,
    /// Strand the ORF lies on.
    pub strand: Strand,
    /// Location on the *original* (forward) sequence, as a [`Span`].
    /// `start..end` always covers start codon through stop codon
    /// inclusive, in forward-strand coordinates.
    pub span: Span,
    /// The translated protein (includes the trailing `*` if a stop was
    /// reached; no `*` if the ORF runs off the sequence end).
    pub protein: Seq,
    /// `true` if the ORF terminates in a stop codon (vs. running off
    /// the end of the sequence).
    pub has_stop: bool,
}

impl Orf {
    /// Length of the ORF in nucleotides (start codon → stop codon
    /// inclusive).
    pub fn nt_len(&self) -> usize {
        self.span.len()
    }

    /// Length of the translated protein, *excluding* a trailing stop.
    pub fn protein_len(&self) -> usize {
        let p = self.protein.as_bytes();
        if p.last() == Some(&b'*') {
            p.len() - 1
        } else {
            p.len()
        }
    }
}

/// Options for [`find_orfs`].
#[derive(Copy, Clone, Debug)]
pub struct OrfOptions {
    /// If `true`, only `ATG` starts an ORF. If `false`, any start
    /// codon of the supplied [`GeneticCode`] does (e.g. `GTG`, `TTG`
    /// for the bacterial table).
    pub atg_only: bool,
    /// Minimum translated-protein length (excluding the stop) for an
    /// ORF to be reported.
    pub min_protein_len: usize,
    /// If `true`, report ORFs that run off the sequence end without a
    /// stop codon. If `false`, only stop-terminated ORFs are kept.
    pub allow_no_stop: bool,
}

impl Default for OrfOptions {
    fn default() -> Self {
        OrfOptions {
            atg_only: true,
            min_protein_len: 30,
            allow_no_stop: false,
        }
    }
}

/// Finds ORFs across all six reading frames.
///
/// Returns ORFs sorted by descending protein length (longest first).
/// Returns [`BioseqError::Invalid`] for a protein input.
pub fn find_orfs(
    seq: &Seq,
    code: &GeneticCode,
    opts: OrfOptions,
) -> Result<Vec<Orf>> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "ORF finding needs a nucleotide sequence",
        ));
    }
    let n = seq.len();
    let mut orfs = Vec::new();

    // Forward frames.
    for offset in 0..3usize {
        scan_frame(
            seq.as_bytes(),
            offset,
            n,
            Strand::Forward,
            offset as u8,
            code,
            &opts,
            &mut orfs,
        );
    }

    // Reverse frames — scan the reverse complement, then map
    // coordinates back to forward-strand space.
    let rc = reverse_complement(seq)?;
    for offset in 0..3usize {
        let mut rc_orfs = Vec::new();
        scan_frame(
            rc.as_bytes(),
            offset,
            n,
            Strand::Reverse,
            3 + offset as u8,
            code,
            &opts,
            &mut rc_orfs,
        );
        // rc coordinate r maps to forward coordinate n - r.
        for mut o in rc_orfs {
            let fwd_start = n - o.span.end;
            let fwd_end = n - o.span.start;
            o.span = Span::with_strand(fwd_start, fwd_end, Strand::Reverse);
            orfs.push(o);
        }
    }

    orfs.sort_by_key(|o| std::cmp::Reverse(o.protein_len()));
    Ok(orfs)
}

/// Scans one reading frame (a fixed offset on `bytes`) for ORFs.
#[allow(clippy::too_many_arguments)]
fn scan_frame(
    bytes: &[u8],
    offset: usize,
    seq_len: usize,
    strand: Strand,
    frame: u8,
    code: &GeneticCode,
    opts: &OrfOptions,
    out: &mut Vec<Orf>,
) {
    let n = bytes.len();
    if offset >= n {
        return;
    }
    let mut i = offset;
    while i + 3 <= n {
        let codon = &bytes[i..i + 3];
        let is_start = if opts.atg_only {
            codon.eq_ignore_ascii_case(b"ATG") || codon.eq_ignore_ascii_case(b"AUG")
        } else {
            code.is_start_codon(codon)
        };
        if is_start {
            // Walk to the next in-frame stop.
            let mut j = i;
            let mut has_stop = false;
            while j + 3 <= n {
                if code.is_stop_codon(&bytes[j..j + 3]) {
                    has_stop = true;
                    j += 3;
                    break;
                }
                j += 3;
            }
            // j is now one-past the stop codon (or off the end).
            if has_stop || opts.allow_no_stop {
                let sub = Seq::new_unchecked(
                    if bytes_look_rna(bytes) {
                        SeqKind::Rna
                    } else {
                        SeqKind::Dna
                    },
                    bytes[i..j].to_vec(),
                    crate::seq::Topology::Linear,
                );
                if let Ok(protein) = translate_default(&sub, code) {
                    let prot_len = if protein.as_bytes().last() == Some(&b'*') {
                        protein.len() - 1
                    } else {
                        protein.len()
                    };
                    if prot_len >= opts.min_protein_len {
                        // Coordinates here are in `bytes` space; the
                        // caller remaps reverse-strand ones.
                        let span = Span::with_strand(i, j, strand);
                        debug_assert!(j <= seq_len || strand == Strand::Reverse);
                        out.push(Orf {
                            frame,
                            strand,
                            span,
                            protein,
                            has_stop,
                        });
                    }
                }
            }
            // Continue scanning *after* this ORF's stop so nested
            // downstream-start ORFs in the same frame are still found,
            // but the same start is not re-reported.
            i = j;
        } else {
            i += 3;
        }
    }
}

/// Heuristic: does this byte slice contain `U`/`u` but no `T`/`t`?
fn bytes_look_rna(bytes: &[u8]) -> bool {
    let mut has_u = false;
    let mut has_t = false;
    for &b in bytes {
        match b.to_ascii_uppercase() {
            b'U' => has_u = true,
            b'T' => has_t = true,
            _ => {}
        }
    }
    has_u && !has_t
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(min: usize) -> OrfOptions {
        OrfOptions {
            atg_only: true,
            min_protein_len: min,
            allow_no_stop: false,
        }
    }

    #[test]
    fn finds_a_simple_forward_orf() {
        let code = GeneticCode::standard();
        // 5 prefix bases, then ATG AAA GGG TAA, then trailing.
        let dna = Seq::new(SeqKind::Dna, "CCCCCATGAAAGGGTAACCCCC").unwrap();
        let orfs = find_orfs(&dna, &code, opts(2)).unwrap();
        assert_eq!(orfs.len(), 1);
        let o = &orfs[0];
        assert_eq!(o.strand, Strand::Forward);
        assert_eq!(o.span.start, 5);
        assert_eq!(o.span.end, 17); // 12 nt: ATGAAAGGGTAA
        assert_eq!(o.protein.as_str(), "MKG*");
        assert_eq!(o.protein_len(), 3);
        assert!(o.has_stop);
    }

    #[test]
    fn finds_reverse_strand_orf() {
        let code = GeneticCode::standard();
        // Build an ORF on the reverse strand: reverse complement of
        // "ATGAAATAA" is "TTATTTCAT"; embed that on the forward strand.
        let dna = Seq::new(SeqKind::Dna, "GGGTTATTTCATGGG").unwrap();
        let orfs = find_orfs(&dna, &code, opts(1)).unwrap();
        let rev: Vec<&Orf> = orfs.iter().filter(|o| o.strand == Strand::Reverse).collect();
        assert!(!rev.is_empty(), "expected a reverse-strand ORF");
        let o = rev[0];
        assert_eq!(o.protein.as_str(), "MK*");
        // Forward-coordinate span should cover the embedded 9 nt.
        assert_eq!(o.span.len(), 9);
        assert_eq!(o.span.start, 3);
    }

    #[test]
    fn min_length_filters() {
        let code = GeneticCode::standard();
        let dna = Seq::new(SeqKind::Dna, "ATGAAATAA").unwrap(); // protein MK*
        // protein_len = 2; require 5 -> filtered out.
        assert!(find_orfs(&dna, &code, opts(5)).unwrap().is_empty());
        assert_eq!(find_orfs(&dna, &code, opts(2)).unwrap().len(), 1);
    }

    #[test]
    fn any_start_mode_uses_table_starts() {
        let code = GeneticCode::by_id(11).unwrap();
        // GTG is a bacterial start codon.
        let dna = Seq::new(SeqKind::Dna, "GTGAAATAA").unwrap();
        let atg_only = OrfOptions {
            atg_only: true,
            min_protein_len: 1,
            allow_no_stop: false,
        };
        assert!(find_orfs(&dna, &code, atg_only).unwrap().is_empty());
        let any = OrfOptions {
            atg_only: false,
            ..atg_only
        };
        assert_eq!(find_orfs(&dna, &code, any).unwrap().len(), 1);
    }

    #[test]
    fn no_stop_is_dropped_unless_allowed() {
        let code = GeneticCode::standard();
        let dna = Seq::new(SeqKind::Dna, "ATGAAAGGG").unwrap(); // runs off end
        let strict = opts(1);
        assert!(find_orfs(&dna, &code, strict).unwrap().is_empty());
        let lenient = OrfOptions {
            allow_no_stop: true,
            ..strict
        };
        let orfs = find_orfs(&dna, &code, lenient).unwrap();
        assert_eq!(orfs.len(), 1);
        assert!(!orfs[0].has_stop);
    }

    #[test]
    fn protein_input_rejected() {
        let code = GeneticCode::standard();
        let p = Seq::new(SeqKind::Protein, "MKVL").unwrap();
        assert!(find_orfs(&p, &code, opts(1)).is_err());
    }
}
