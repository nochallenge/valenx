//! DNA assembly planning — features 27 and 28.
//!
//! Three of the dominant cloning methods for stitching DNA fragments
//! into a construct, each as a real planner (the j5 / Synbiopython
//! function set):
//!
//! - **Gibson assembly** ([`plan_gibson`], feature 27) — fragments are
//!   joined seamlessly by designing **terminal sequence overlaps**: an
//!   exonuclease chews back the ends, the complementary overlaps
//!   anneal, and polymerase + ligase seal the gaps. The planner
//!   computes, for each adjacent pair, the overlap sequence and checks
//!   its length and melting temperature against design targets.
//! - **Golden Gate** ([`plan_golden_gate`], feature 28) — a Type IIS
//!   restriction enzyme (e.g. *BsaI*) cuts *outside* its recognition
//!   site, leaving programmable 4-nt overhangs; fragments with
//!   matching overhangs ligate in a defined order. The planner assigns
//!   a unique fusion-site overhang to every junction and verifies the
//!   set is collision-free.
//! - **BioBrick** ([`plan_biobrick`], feature 28) — the RFC-10
//!   idempotent standard-assembly scheme: a prefix + suffix of fixed
//!   restriction sites lets any two BioBrick parts be joined into a
//!   new BioBrick. The planner emits the ordered digest / ligation
//!   recipe.
//!
//! Every planner returns a typed *plan* describing the junctions, plus
//! the assembled product sequence, and reports design problems
//! (overlaps too short, overhang collisions, a forbidden internal
//! restriction site) as a [`SysbioError`] or in the plan's warnings.
//!
//! ## v1 caveats
//!
//! Melting temperatures use the Wallace rule (`2·(A+T) + 4·(G+C)` °C)
//! — the standard quick estimate, not a nearest-neighbour calculation.
//! The Golden Gate planner assigns overhangs from a fixed
//! high-fidelity set and checks pairwise distinctness, not the full
//! NEB ligation-fidelity matrix. Internal-site scanning is exact for
//! the enzyme's literal recognition word.

use valenx_bioseq::ops::revcomp::reverse_complement_dna_bytes;
use valenx_bioseq::{Seq, SeqKind};

use crate::error::{Result, SysbioError};

/// Wallace-rule melting temperature of a short DNA oligo, in °C:
/// `2·(#A+#T) + 4·(#G+#C)`. The standard rule-of-thumb for primers
/// and assembly overlaps under ~14 nt; a reasonable estimate above.
pub fn wallace_tm(bytes: &[u8]) -> f64 {
    let mut at = 0;
    let mut gc = 0;
    for &b in bytes {
        match b.to_ascii_uppercase() {
            b'A' | b'T' => at += 1,
            b'G' | b'C' => gc += 1,
            _ => {}
        }
    }
    2.0 * at as f64 + 4.0 * gc as f64
}

// --- Gibson assembly --------------------------------------------------

/// One designed junction between two Gibson fragments.
#[derive(Debug, Clone, PartialEq)]
pub struct GibsonJunction {
    /// Index of the upstream fragment.
    pub left_fragment: usize,
    /// Index of the downstream fragment.
    pub right_fragment: usize,
    /// The overlap sequence shared by the two fragment ends.
    pub overlap: String,
    /// Wallace-rule melting temperature of the overlap (°C).
    pub overlap_tm: f64,
}

/// A complete Gibson assembly plan.
#[derive(Debug, Clone, PartialEq)]
pub struct GibsonPlan {
    /// The designed junctions, in assembly order.
    pub junctions: Vec<GibsonJunction>,
    /// The final assembled sequence.
    pub product: Seq,
    /// Whether the construct is circular (the last fragment overlaps
    /// the first).
    pub circular: bool,
    /// Non-fatal design warnings (short overlaps, low Tm).
    pub warnings: Vec<String>,
}

/// Plan a Gibson assembly of `fragments` (feature 27).
///
/// `overlap_len` is the target overlap length to design at each
/// junction. For a *linear* product the fragments are assumed already
/// to share `overlap_len` bp at their abutting ends (Gibson primers
/// add these); the planner extracts and characterises each overlap.
/// For a *circular* product the last fragment also overlaps the first.
///
/// Errors if fewer than two fragments are given, if a fragment is
/// shorter than the overlap, or if an abutting pair does not actually
/// share the claimed overlap.
pub fn plan_gibson(
    fragments: &[Seq],
    overlap_len: usize,
    circular: bool,
    min_tm: f64,
) -> Result<GibsonPlan> {
    if fragments.len() < 2 {
        return Err(SysbioError::invalid(
            "fragments",
            "Gibson assembly needs at least two fragments",
        ));
    }
    if overlap_len == 0 {
        return Err(SysbioError::invalid("overlap_len", "overlap must be positive"));
    }
    for (i, f) in fragments.iter().enumerate() {
        if f.len() < overlap_len {
            return Err(SysbioError::invalid(
                "fragments",
                format!("fragment {i} is shorter than the {overlap_len} bp overlap"),
            ));
        }
    }

    let n = fragments.len();
    let pairs = if circular { n } else { n - 1 };
    let mut junctions = Vec::with_capacity(pairs);
    let mut warnings = Vec::new();

    for j in 0..pairs {
        let left = j;
        let right = (j + 1) % n;
        let lb = fragments[left].as_bytes();
        let rb = fragments[right].as_bytes();
        // The overlap is the 3' end of `left` and must equal the 5'
        // start of `right`.
        let left_tail = &lb[lb.len() - overlap_len..];
        let right_head = &rb[..overlap_len];
        if !left_tail.eq_ignore_ascii_case(right_head) {
            return Err(SysbioError::invalid_model(
                "gibson",
                format!(
                    "fragments {left} and {right} do not share the claimed {overlap_len} bp overlap"
                ),
            ));
        }
        let overlap = String::from_utf8_lossy(left_tail).to_uppercase();
        let tm = wallace_tm(left_tail);
        if tm < min_tm {
            warnings.push(format!(
                "junction {left}->{right} overlap Tm {tm:.1}C is below the {min_tm:.1}C target"
            ));
        }
        junctions.push(GibsonJunction {
            left_fragment: left,
            right_fragment: right,
            overlap,
            overlap_tm: tm,
        });
    }

    // Assemble: concatenate fragments, dropping the shared overlap once
    // per junction so each overlap appears exactly once.
    let mut bytes: Vec<u8> = fragments[0].as_bytes().to_vec();
    for f in &fragments[1..] {
        bytes.extend_from_slice(&f.as_bytes()[overlap_len..]);
    }
    if circular {
        // The final fragment's tail overlaps fragment 0's head — trim
        // it off the running product.
        bytes.truncate(bytes.len() - overlap_len);
    }
    let product = Seq::new(SeqKind::Dna, bytes)
        .map_err(|e| SysbioError::invalid_model("gibson", format!("bad product: {e}")))?;

    Ok(GibsonPlan {
        junctions,
        product,
        circular,
        warnings,
    })
}

/// Design Gibson fragments from a single target sequence by adding
/// overlaps.
///
/// Splits `target` into `n_fragments` pieces and *extends* each
/// piece's ends so adjacent pieces share `overlap_len` bp — i.e. the
/// inverse of [`plan_gibson`], producing the fragments a wet-lab
/// would order. Returns the fragment sequences.
pub fn design_gibson_fragments(
    target: &Seq,
    n_fragments: usize,
    overlap_len: usize,
) -> Result<Vec<Seq>> {
    if n_fragments < 2 {
        return Err(SysbioError::invalid("n_fragments", "need at least two"));
    }
    let total = target.len();
    if total < n_fragments * (overlap_len + 1) {
        return Err(SysbioError::invalid(
            "target",
            "sequence too short to split with the requested overlaps",
        ));
    }
    let bytes = target.as_bytes();
    let base = total / n_fragments;
    let mut frags = Vec::with_capacity(n_fragments);
    for i in 0..n_fragments {
        let start = i * base;
        let end = if i + 1 == n_fragments {
            total
        } else {
            // Extend the right end into the next block to create the
            // shared overlap.
            (start + base + overlap_len).min(total)
        };
        let slice = bytes[start..end].to_vec();
        frags.push(
            Seq::new(SeqKind::Dna, slice).map_err(|e| {
                SysbioError::invalid_model("gibson", format!("bad fragment: {e}"))
            })?,
        );
    }
    Ok(frags)
}

// --- Golden Gate assembly --------------------------------------------

/// One Golden Gate junction with its programmed fusion-site overhang.
#[derive(Debug, Clone, PartialEq)]
pub struct GoldenGateJunction {
    /// Upstream fragment index.
    pub left_fragment: usize,
    /// Downstream fragment index.
    pub right_fragment: usize,
    /// The 4-nt fusion-site overhang ligating the two fragments.
    pub overhang: String,
}

/// A Golden Gate assembly plan.
#[derive(Debug, Clone, PartialEq)]
pub struct GoldenGatePlan {
    /// Recognition site of the Type IIS enzyme used.
    pub enzyme_site: String,
    /// The junctions with their assigned overhangs.
    pub junctions: Vec<GoldenGateJunction>,
    /// Design warnings (e.g. an internal enzyme site).
    pub warnings: Vec<String>,
}

/// A small high-fidelity 4-nt fusion-site set (a subset of the
/// commonly used Potapov et al. high-fidelity overhangs).
const FUSION_SITES: [&str; 12] = [
    "AATG", "GCTT", "CGCT", "TGCC", "ACTA", "TTAC", "GGTA", "AGGA", "ATCC", "CAGG",
    "TACG", "GACT",
];

/// Plan a Golden Gate assembly of `n_fragments` fragments (feature
/// 28).
///
/// `enzyme_site` is the Type IIS recognition word (e.g. `"GGTCTC"`
/// for *BsaI*). The planner assigns a distinct 4-nt overhang to every
/// junction from the built-in high-fidelity set and — if the fragment
/// sequences are supplied via [`plan_golden_gate_checked`] — flags any
/// fragment that already contains the enzyme site internally (which
/// would be cut and break the assembly).
pub fn plan_golden_gate(
    n_fragments: usize,
    circular: bool,
    enzyme_site: &str,
) -> Result<GoldenGatePlan> {
    if n_fragments < 2 {
        return Err(SysbioError::invalid(
            "n_fragments",
            "Golden Gate needs at least two fragments",
        ));
    }
    let n_junctions = if circular {
        n_fragments
    } else {
        n_fragments - 1
    };
    if n_junctions > FUSION_SITES.len() {
        return Err(SysbioError::invalid(
            "n_fragments",
            format!(
                "needs {n_junctions} distinct overhangs but only {} are available",
                FUSION_SITES.len()
            ),
        ));
    }
    let mut junctions = Vec::with_capacity(n_junctions);
    for (j, &site) in FUSION_SITES.iter().take(n_junctions).enumerate() {
        junctions.push(GoldenGateJunction {
            left_fragment: j,
            right_fragment: (j + 1) % n_fragments,
            overhang: site.to_string(),
        });
    }
    Ok(GoldenGatePlan {
        enzyme_site: enzyme_site.to_uppercase(),
        junctions,
        warnings: Vec::new(),
    })
}

/// Like [`plan_golden_gate`], but also scans each fragment for an
/// internal copy of the enzyme recognition site (on either strand)
/// and records a warning per offending fragment.
pub fn plan_golden_gate_checked(
    fragments: &[Seq],
    circular: bool,
    enzyme_site: &str,
) -> Result<GoldenGatePlan> {
    let mut plan = plan_golden_gate(fragments.len(), circular, enzyme_site)?;
    let site = enzyme_site.to_uppercase().into_bytes();
    let rc_site = reverse_complement_dna_bytes(&site);
    for (i, f) in fragments.iter().enumerate() {
        let up: Vec<u8> = f.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
        let fwd = count_occurrences(&up, &site);
        let rev = if rc_site == site {
            0
        } else {
            count_occurrences(&up, &rc_site)
        };
        if fwd + rev > 0 {
            plan.warnings.push(format!(
                "fragment {i} contains {} internal {enzyme_site} site(s) — domesticate before assembly",
                fwd + rev
            ));
        }
    }
    Ok(plan)
}

/// Count non-overlapping occurrences of `needle` in `haystack`.
fn count_occurrences(haystack: &[u8], needle: &[u8]) -> usize {
    if needle.is_empty() || haystack.len() < needle.len() {
        return 0;
    }
    let mut count = 0;
    let mut i = 0;
    while i + needle.len() <= haystack.len() {
        if &haystack[i..i + needle.len()] == needle {
            count += 1;
            i += needle.len();
        } else {
            i += 1;
        }
    }
    count
}

// --- BioBrick (RFC[10]) assembly -------------------------------------

/// One step of a BioBrick standard-assembly recipe.
#[derive(Debug, Clone, PartialEq)]
pub struct BioBrickStep {
    /// Human-readable description of the digest / ligation step.
    pub description: String,
}

/// A BioBrick (RFC-10) standard-assembly plan.
#[derive(Debug, Clone, PartialEq)]
pub struct BioBrickPlan {
    /// The ordered digest / ligation steps.
    pub steps: Vec<BioBrickStep>,
    /// The assembled composite-part sequence (prefix + parts + suffix).
    pub product: Seq,
}

/// The RFC-10 BioBrick prefix (EcoRI + NotI + XbaI), uppercase.
pub const BIOBRICK_PREFIX: &str = "GAATTCGCGGCCGCTTCTAGAG";
/// The RFC-10 BioBrick suffix (SpeI + NotI + PstI), uppercase.
pub const BIOBRICK_SUFFIX: &str = "TACTAGTAGCGGCCGCTGCAG";

/// Plan a BioBrick RFC-10 standard assembly of `parts` (feature 28).
///
/// BioBrick assembly is *idempotent*: joining two standard parts
/// yields another standard part flanked by the same prefix / suffix.
/// The planner emits the canonical two-enzyme recipe — cut the
/// upstream part with the suffix enzymes, the downstream part with the
/// prefix enzymes, ligate through the compatible *SpeI* / *XbaI* scar
/// — and returns the assembled sequence: prefix, every part (joined by
/// the 6-bp mixed `TACTAGAG` scar), suffix.
pub fn plan_biobrick(parts: &[Seq]) -> Result<BioBrickPlan> {
    if parts.is_empty() {
        return Err(SysbioError::invalid("parts", "need at least one BioBrick part"));
    }
    // Idempotent-scar check: a part must not contain a forbidden RFC10
    // restriction site internally (EcoRI, XbaI, SpeI, PstI, NotI).
    let forbidden: [(&str, &[u8]); 5] = [
        ("EcoRI", b"GAATTC"),
        ("XbaI", b"TCTAGA"),
        ("SpeI", b"ACTAGT"),
        ("PstI", b"CTGCAG"),
        ("NotI", b"GCGGCCGC"),
    ];
    let mut steps = Vec::new();
    for (i, p) in parts.iter().enumerate() {
        let up: Vec<u8> = p.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
        for (name, site) in &forbidden {
            if count_occurrences(&up, site) > 0 {
                return Err(SysbioError::invalid_model(
                    "biobrick",
                    format!("part {i} contains an internal {name} site — not RFC[10] compliant"),
                ));
            }
        }
    }

    // The "TACTAGAG" mixed scar left between two ligated BioBricks.
    const SCAR: &str = "TACTAGAG";
    let mut bytes: Vec<u8> = BIOBRICK_PREFIX.as_bytes().to_vec();
    for (i, p) in parts.iter().enumerate() {
        if i > 0 {
            bytes.extend_from_slice(SCAR.as_bytes());
            steps.push(BioBrickStep {
                description: format!(
                    "ligate part {i} downstream of part {} through the SpeI/XbaI scar",
                    i - 1
                ),
            });
        } else {
            steps.push(BioBrickStep {
                description: "place part 0 immediately after the BioBrick prefix".into(),
            });
        }
        bytes.extend_from_slice(p.as_bytes());
    }
    bytes.extend_from_slice(BIOBRICK_SUFFIX.as_bytes());
    steps.push(BioBrickStep {
        description: "append the BioBrick suffix to complete the composite part".into(),
    });

    let product = Seq::new(SeqKind::Dna, bytes)
        .map_err(|e| SysbioError::invalid_model("biobrick", format!("bad product: {e}")))?;
    Ok(BioBrickPlan { steps, product })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(s: &str) -> Seq {
        Seq::new(SeqKind::Dna, s).unwrap()
    }

    #[test]
    fn wallace_tm_counts_bases() {
        // 4 G/C * 4 + 0 A/T = 16.
        assert_eq!(wallace_tm(b"GCGC"), 16.0);
        // 4 A/T * 2 = 8.
        assert_eq!(wallace_tm(b"ATAT"), 8.0);
    }

    #[test]
    fn gibson_design_then_plan_roundtrips() {
        // Split a known sequence into overlapping fragments, then plan
        // the assembly and confirm the product matches the original.
        let target = dna(
            "ATGCGTACGTTAGCCGGATCAATTGGCCAATTACGTACGTAGCTAGCTAGGGCCCTTTAAACCCGGGTTT",
        );
        let frags = design_gibson_fragments(&target, 3, 12).unwrap();
        assert_eq!(frags.len(), 3);
        let plan = plan_gibson(&frags, 12, false, 24.0).unwrap();
        assert_eq!(plan.product.as_bytes(), target.as_bytes());
        assert_eq!(plan.junctions.len(), 2);
    }

    #[test]
    fn gibson_circular_trims_wraparound_overlap() {
        // Three fragments arranged in a circle, each sharing 8 bp with
        // the next (including frag2 -> frag0).
        let f0 = dna("AAAAAAAACCCCCCCCGGGGGGGG");
        let f1 = dna("GGGGGGGGTTTTTTTTACACACAC");
        let f2 = dna("ACACACACGTGTGTGTAAAAAAAA");
        let plan = plan_gibson(&[f0, f1, f2], 8, true, 16.0).unwrap();
        assert!(plan.circular);
        assert_eq!(plan.junctions.len(), 3);
        // Circular product length = sum - overlap*njunctions.
        assert_eq!(plan.product.len(), 24 * 3 - 8 * 3);
    }

    #[test]
    fn gibson_rejects_non_matching_overlap() {
        let f0 = dna("AAAAAAAACCCCCCCC");
        let f1 = dna("TTTTTTTTGGGGGGGG"); // head != f0 tail
        assert!(plan_gibson(&[f0, f1], 8, false, 16.0).is_err());
    }

    #[test]
    fn gibson_low_tm_overlap_warns() {
        // An all-AT 8-bp overlap has Tm 16C; demand 30C -> warning.
        let f0 = dna("CCCCCCCCATATATAT");
        let f1 = dna("ATATATATGGGGGGGG");
        let plan = plan_gibson(&[f0, f1], 8, false, 30.0).unwrap();
        assert_eq!(plan.warnings.len(), 1);
        assert!(plan.warnings[0].contains("Tm"));
    }

    #[test]
    fn golden_gate_assigns_distinct_overhangs() {
        let plan = plan_golden_gate(4, false, "GGTCTC").unwrap();
        assert_eq!(plan.junctions.len(), 3);
        let overhangs: Vec<&str> =
            plan.junctions.iter().map(|j| j.overhang.as_str()).collect();
        // All distinct.
        for (i, a) in overhangs.iter().enumerate() {
            for b in &overhangs[i + 1..] {
                assert_ne!(a, b);
            }
        }
        assert_eq!(plan.enzyme_site, "GGTCTC");
    }

    #[test]
    fn golden_gate_flags_internal_enzyme_site() {
        // One fragment contains BsaI (GGTCTC) internally.
        let clean = dna("AAAACCCCGGGGTTTT");
        let dirty = dna("AAAAGGTCTCAAAATTTT");
        let plan = plan_golden_gate_checked(&[clean, dirty], false, "GGTCTC").unwrap();
        assert_eq!(plan.warnings.len(), 1);
        assert!(plan.warnings[0].contains("fragment 1"));
    }

    #[test]
    fn golden_gate_too_many_fragments_errors() {
        // More junctions than available fusion sites.
        assert!(plan_golden_gate(20, false, "GGTCTC").is_err());
    }

    #[test]
    fn biobrick_assembly_has_prefix_and_suffix() {
        let p0 = dna("ATGAAACGT");
        let p1 = dna("GGGTTTAAA");
        let plan = plan_biobrick(&[p0, p1]).unwrap();
        let prod = plan.product.as_str();
        assert!(prod.starts_with(BIOBRICK_PREFIX));
        assert!(prod.ends_with(BIOBRICK_SUFFIX));
        // A scar joins the two parts -> at least 3 recipe steps.
        assert!(plan.steps.len() >= 3);
    }

    #[test]
    fn biobrick_rejects_internal_restriction_site() {
        // A part with an internal EcoRI site (GAATTC) is not RFC[10].
        let bad = dna("ATGGAATTCAAA");
        assert!(plan_biobrick(&[bad]).is_err());
    }

    #[test]
    fn biobrick_single_part_is_valid() {
        let p = dna("ATGAAACGTGGG");
        let plan = plan_biobrick(&[p]).unwrap();
        assert!(plan.product.as_str().contains("ATGAAACGTGGG"));
    }
}
