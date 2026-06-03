//! In-silico PCR — predict the amplicon(s) a primer pair produces.
//!
//! Given a template (linear or circular) and a forward + reverse
//! primer, [`simulate_pcr`] finds every place each primer anneals
//! (allowing a configurable number of 5′ mismatches but requiring an
//! exact 3′ anchor), pairs compatible forward/reverse sites, and
//! reports the predicted amplicons with their sizes.

use crate::error::{BioseqError, Result};
use crate::ops::revcomp::reverse_complement;
use crate::seq::{Seq, SeqKind, Topology};

/// A predicted PCR amplicon.
#[derive(Clone, Debug, PartialEq)]
pub struct Amplicon {
    /// 0-based start on the template top strand (the 5′ end of the
    /// forward primer's binding site).
    pub start: usize,
    /// 0-based end (exclusive) — one past the forward extent of the
    /// reverse primer's binding site.
    pub end: usize,
    /// Amplicon length in bp.
    pub size: usize,
    /// `true` if the amplicon spans the origin of a circular template.
    pub wraps_origin: bool,
    /// The amplicon sequence (with the primer sequences as the
    /// termini).
    pub product: Seq,
}

/// A binding site of a primer on the template top strand.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct BindingSite {
    /// 0-based start of the footprint on the top strand.
    start: usize,
    /// 0-based end (exclusive).
    end: usize,
    /// Number of mismatches in the 5′ region.
    mismatches: usize,
}

/// Options for [`simulate_pcr`].
#[derive(Copy, Clone, Debug)]
pub struct PcrOptions {
    /// Maximum mismatches allowed in the 5′ portion of a primer. The
    /// 3′ anchor (see `anchor_len`) must always match exactly.
    pub max_mismatches: usize,
    /// Length of the 3′ anchor that must match perfectly. Mismatches
    /// here abolish extension, so PCR realism requires an exact 3′ end.
    pub anchor_len: usize,
    /// Maximum amplicon size to report, bp. Amplicons longer than this
    /// are discarded (models the processivity limit of a polymerase).
    pub max_amplicon: usize,
}

impl Default for PcrOptions {
    fn default() -> Self {
        PcrOptions {
            max_mismatches: 0,
            anchor_len: 12,
            max_amplicon: 20_000,
        }
    }
}

/// Simulates PCR with a forward + reverse primer on a template.
///
/// The forward primer must match the top strand; the reverse primer
/// matches the bottom strand (i.e. its reverse complement matches the
/// top strand). Amplicons are formed by every forward site paired with
/// a downstream reverse site. For a circular template the search wraps
/// the origin.
///
/// Returns [`BioseqError::Invalid`] for a non-DNA input or empty
/// primer.
pub fn simulate_pcr(
    template: &Seq,
    forward: &Seq,
    reverse: &Seq,
    opts: PcrOptions,
) -> Result<Vec<Amplicon>> {
    if template.kind() != SeqKind::Dna
        || forward.kind() != SeqKind::Dna
        || reverse.kind() != SeqKind::Dna
    {
        return Err(BioseqError::invalid("kind", "in-silico PCR needs DNA inputs"));
    }
    if forward.is_empty() || reverse.is_empty() {
        return Err(BioseqError::invalid("primer", "primers must be non-empty"));
    }
    let n = template.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Build the search string. For a circular template, search a
    // doubled copy so origin-spanning binding sites are found.
    let top: Vec<u8> = if template.is_circular() {
        let mut v = template.as_bytes().to_vec();
        let dup = v.clone();
        v.extend_from_slice(&dup);
        v
    } else {
        template.as_bytes().to_vec()
    };

    // Forward primer binds the top strand directly. Each forward site
    // is unique modulo n, so wrap-copies in the doubled string are
    // discarded.
    let fwd_sites =
        find_binding_sites(&top, forward.as_bytes(), &opts, n, template.is_circular(), false);
    // Reverse primer binds the bottom strand: its reverse complement
    // matches the top strand. The reverse primer's 3' end corresponds
    // to the *low-coordinate* end of the matched top-strand footprint.
    // For a circular template the reverse-site search KEEPS the
    // wrap-copies in the doubled string [n, 2n): that is how a forward
    // site near the 3' end pairs with a reverse site that lies just past
    // the origin, producing an origin-spanning amplicon.
    let rev_rc = reverse_complement(reverse)?;
    let rev_sites_raw =
        find_binding_sites(&top, rev_rc.as_bytes(), &opts, n, template.is_circular(), true);

    // The 3' anchor of the reverse primer is at the START of the
    // top-strand footprint (because rev_rc is reverse-complemented),
    // so re-screen with the anchor at the left end.
    let rev_sites: Vec<BindingSite> = rev_sites_raw
        .into_iter()
        .filter(|site| {
            // Anchor check at the 5' (left) end of the rc footprint.
            anchor_matches_left(&top, site.start, rev_rc.as_bytes(), opts.anchor_len)
        })
        .collect();

    let mut amplicons = Vec::new();
    for f in &fwd_sites {
        for r in &rev_sites {
            // A valid amplicon: forward site, then a reverse site whose
            // top-strand footprint is downstream of (or overlapping the
            // 3' side of) the forward site.
            let amp_start = f.start;
            let amp_end = r.end;
            if amp_end <= amp_start {
                continue;
            }
            let size = amp_end - amp_start;
            if size > opts.max_amplicon || size == 0 {
                continue;
            }
            // On a circular template a reverse site may only pair with a
            // forward site that lies within one full turn upstream — a
            // reverse wrap-copy more than n away is the same physical
            // site reached the "long way round" and is not a product.
            if template.is_circular() && size > n {
                continue;
            }
            // Skip the degenerate case where the two binding sites are
            // the same span (a primer pairing with itself).
            if f.start == r.start && f.end == r.end {
                continue;
            }
            let wraps = template.is_circular() && amp_end > n;
            let product_bytes: Vec<u8> = top[amp_start..amp_end].to_vec();
            let product = Seq::new_unchecked(SeqKind::Dna, product_bytes, Topology::Linear);
            amplicons.push(Amplicon {
                start: amp_start % n,
                end: if wraps { amp_end % n } else { amp_end },
                size,
                wraps_origin: wraps,
                product,
            });
        }
    }
    // Deterministic order: by start, then size.
    amplicons.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.size.cmp(&b.size)));
    amplicons.dedup_by(|a, b| a.start == b.start && a.size == b.size);
    Ok(amplicons)
}

/// Finds every binding site of `primer` on `search` (the top strand,
/// possibly doubled for a circular template). A site requires the 3′
/// anchor to match exactly and at most `max_mismatches` mismatches in
/// the rest.
///
/// `keep_wrap_copies` controls the circular doubled-string copies in
/// `[template_len, 2·template_len)`: pass `false` for a primer whose
/// sites are unique modulo the length (the forward primer), `true` for
/// the downstream (reverse) primer so an origin-spanning amplicon — a
/// forward site near the 3′ end paired with a reverse site just past
/// the origin — can be formed.
fn find_binding_sites(
    search: &[u8],
    primer: &[u8],
    opts: &PcrOptions,
    template_len: usize,
    circular: bool,
    keep_wrap_copies: bool,
) -> Vec<BindingSite> {
    let m = primer.len();
    let mut sites = Vec::new();
    if m == 0 || m > search.len() {
        return sites;
    }
    let anchor = opts.anchor_len.min(m);
    for start in 0..=search.len() - m {
        // A site beginning at or past the original length is a wrap
        // duplicate (circular). Skip it to avoid double-reporting unless
        // the caller wants the wrap-copies (the reverse primer).
        if circular && start >= template_len && !keep_wrap_copies {
            continue;
        }
        let window = &search[start..start + m];
        // 3' anchor must match exactly (the last `anchor` bases).
        let anchor_ok = window[m - anchor..]
            .iter()
            .zip(&primer[m - anchor..])
            .all(|(a, b)| a.eq_ignore_ascii_case(b));
        if !anchor_ok {
            continue;
        }
        // Count mismatches over the whole primer.
        let mismatches = window
            .iter()
            .zip(primer)
            .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
            .count();
        if mismatches <= opts.max_mismatches {
            sites.push(BindingSite {
                start,
                end: start + m,
                mismatches,
            });
        }
    }
    sites
}

/// Anchor check at the LEFT end of a footprint — for the reverse
/// primer's reverse complement, where the 3′ anchor sits at the start.
fn anchor_matches_left(search: &[u8], start: usize, primer_rc: &[u8], anchor_len: usize) -> bool {
    let anchor = anchor_len.min(primer_rc.len());
    if start + anchor > search.len() {
        return false;
    }
    search[start..start + anchor]
        .iter()
        .zip(&primer_rc[..anchor])
        .all(|(a, b)| a.eq_ignore_ascii_case(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn amplifies_a_simple_region() {
        // Template: [fwd primer][insert][reverse-comp of rev primer]
        let fwd = "ATGCATGCATGC";
        let insert = "AAAATTTTGGGGCCCC";
        let rev = "TTAATTAATTAA";
        let rev_rc: Vec<u8> = crate::ops::revcomp::reverse_complement_dna_bytes(rev.as_bytes());
        let mut tmpl = fwd.as_bytes().to_vec();
        tmpl.extend_from_slice(insert.as_bytes());
        tmpl.extend_from_slice(&rev_rc);
        let template = Seq::new(SeqKind::Dna, tmpl).unwrap();
        let f = Seq::new(SeqKind::Dna, fwd).unwrap();
        let r = Seq::new(SeqKind::Dna, rev).unwrap();

        let amps = simulate_pcr(&template, &f, &r, PcrOptions::default()).unwrap();
        assert_eq!(amps.len(), 1);
        let a = &amps[0];
        assert_eq!(a.start, 0);
        // amplicon = fwd + insert + rev_rc = 12 + 16 + 12 = 40.
        assert_eq!(a.size, 40);
        assert!(a.product.as_str().starts_with("ATGCATGC"));
        assert!(!a.wraps_origin);
    }

    #[test]
    fn no_product_when_primers_absent() {
        let template = Seq::new(SeqKind::Dna, "AAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        let f = Seq::new(SeqKind::Dna, "GGGGGGGGGGGG").unwrap();
        let r = Seq::new(SeqKind::Dna, "CCCCCCCCCCCC").unwrap();
        let amps = simulate_pcr(&template, &f, &r, PcrOptions::default()).unwrap();
        assert!(amps.is_empty());
    }

    #[test]
    fn amplicon_size_limit_filters() {
        let fwd = "ATGCATGCATGC";
        let rev = "TTAATTAATTAA";
        let rev_rc = crate::ops::revcomp::reverse_complement_dna_bytes(rev.as_bytes());
        let mut tmpl = fwd.as_bytes().to_vec();
        tmpl.extend(std::iter::repeat_n(b'A', 100));
        tmpl.extend_from_slice(&rev_rc);
        let template = Seq::new(SeqKind::Dna, tmpl).unwrap();
        let f = Seq::new(SeqKind::Dna, fwd).unwrap();
        let r = Seq::new(SeqKind::Dna, rev).unwrap();
        let opts = PcrOptions {
            max_amplicon: 50, // amplicon would be 124 bp -> filtered
            ..Default::default()
        };
        assert!(simulate_pcr(&template, &f, &r, opts).unwrap().is_empty());
        let opts_ok = PcrOptions {
            max_amplicon: 500,
            ..Default::default()
        };
        assert_eq!(simulate_pcr(&template, &f, &r, opts_ok).unwrap().len(), 1);
    }

    #[test]
    fn three_prime_mismatch_abolishes_binding() {
        let fwd = "ATGCATGCATGC";
        let rev = "TTAATTAATTAA";
        let rev_rc = crate::ops::revcomp::reverse_complement_dna_bytes(rev.as_bytes());
        let mut tmpl = fwd.as_bytes().to_vec();
        tmpl.extend_from_slice(b"AAAATTTT");
        tmpl.extend_from_slice(&rev_rc);
        let template = Seq::new(SeqKind::Dna, tmpl).unwrap();
        // Mutate the forward primer's 3' end.
        let bad_fwd = Seq::new(SeqKind::Dna, "ATGCATGCATGA").unwrap();
        let r = Seq::new(SeqKind::Dna, rev).unwrap();
        let opts = PcrOptions {
            max_mismatches: 3,
            anchor_len: 4,
            ..Default::default()
        };
        // The 3' anchor mismatch should abolish the product.
        assert!(simulate_pcr(&template, &bad_fwd, &r, opts).unwrap().is_empty());
    }

    #[test]
    fn five_prime_mismatch_tolerated_when_allowed() {
        let fwd = "ATGCATGCATGC";
        let rev = "TTAATTAATTAA";
        let rev_rc = crate::ops::revcomp::reverse_complement_dna_bytes(rev.as_bytes());
        let mut tmpl = fwd.as_bytes().to_vec();
        tmpl.extend_from_slice(b"AAAATTTT");
        tmpl.extend_from_slice(&rev_rc);
        let template = Seq::new(SeqKind::Dna, tmpl).unwrap();
        // Mutate the forward primer's 5' end (position 0).
        let tagged_fwd = Seq::new(SeqKind::Dna, "GTGCATGCATGC").unwrap();
        let r = Seq::new(SeqKind::Dna, rev).unwrap();
        let opts = PcrOptions {
            max_mismatches: 1,
            anchor_len: 6,
            ..Default::default()
        };
        // A single 5' mismatch within the budget still amplifies.
        assert_eq!(simulate_pcr(&template, &tagged_fwd, &r, opts).unwrap().len(), 1);
    }

    #[test]
    fn circular_template_amplicon_can_wrap_origin() {
        // Build a circular template where the product spans the origin:
        // [insert tail][rev_rc][filler][fwd][insert head]
        let fwd = "ATGCATGCATGC";
        let rev = "TTAATTAATTAA";
        let rev_rc = crate::ops::revcomp::reverse_complement_dna_bytes(rev.as_bytes());
        let mut tmpl = Vec::new();
        tmpl.extend_from_slice(&rev_rc); // 0..12
        tmpl.extend_from_slice(b"GGGGGGGG"); // filler
        tmpl.extend_from_slice(fwd.as_bytes()); // forward site near the end
        let template =
            Seq::with_topology(SeqKind::Dna, tmpl, Topology::Circular).unwrap();
        let f = Seq::new(SeqKind::Dna, fwd).unwrap();
        let r = Seq::new(SeqKind::Dna, rev).unwrap();
        let amps = simulate_pcr(&template, &f, &r, PcrOptions::default()).unwrap();
        // The forward primer is near the 3' end; the reverse site is at
        // the start -> the amplicon wraps the origin.
        assert!(!amps.is_empty(), "circular PCR should yield a product");
        assert!(amps.iter().any(|a| a.wraps_origin));
    }

    #[test]
    fn non_dna_rejected() {
        let t = Seq::new(SeqKind::Rna, "ACGUACGU").unwrap();
        let f = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let r = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        assert!(simulate_pcr(&t, &f, &r, PcrOptions::default()).is_err());
    }
}
