//! **Feature 7 — fragment library (PDB-mined-style realistic geometry).**
//!
//! A fragment library is the heart of Rosetta-class ab-initio
//! prediction. For each position in the target sequence it holds a
//! set of short backbone *fragments* — strings of (φ, ψ, ω) dihedral
//! angles taken from real protein geometry — that are plausible for
//! that sequence window. The Monte-Carlo assembler builds a fold by
//! repeatedly splicing fragments in.
//!
//! A production Rosetta pipeline *mines* fragments from a curated
//! ~10⁴-structure PDB set: for each window it finds the
//! structurally-best matches by sequence + predicted-secondary-
//! structure similarity. This module builds a **PDB-derived-style
//! curated fragment library**: rather than the previous canonical-
//! Ramachandran-basin generator, fragments are drawn from a small
//! curated table of **published (φ, ψ, ω) means for the standard
//! protein backbone fragment classes** — α-helix interior / N-cap /
//! C-cap, β-strand interior / edge, the four canonical β-turn types
//! (I, II, I', II'), γ-turn (classic + inverse), π-helix, 3₁₀-helix
//! interior, polyproline-II — keyed by sequence window + predicted SS
//! class.
//!
//! The geometry comes from **the published per-residue (φ, ψ) means
//! and spreads** for each class (Lovell *et al.* 2003, Hovmöller *et
//! al.* 2002, the PROCHECK / Rost-Sander reviews) rather than from
//! single idealised Ramachandran basin centres. Multiple geometric
//! realisations per class are produced by perturbing the published
//! mean by **the published one-σ spread** (typically ±10°-20°
//! depending on class — α is tight, turns are looser), so the library
//! returns the realistic spread real PDB-mined fragments show.
//!
//! **Honest scope.** This is a *curated* library, not millions of
//! mined chains: it covers the canonical backbone motifs Rosetta uses
//! for its core fragment classes, but a production Rosetta library
//! adds the much larger long tail of position-specific PDB matches.
//! Treat this as the working classical fragment-assembly v1.

use serde::{Deserialize, Serialize};

use crate::abinitio::ss::{predict_secondary_structure, SecondaryStructure};
use crate::error::{Result, StructPredictError};

/// One PDB-curated backbone fragment class. Each class holds the
/// published `(φ, ψ, ω)` mean dihedrals (degrees) for *every residue
/// position* in the class plus the published one-σ spread used to
/// sample neighbours.
#[derive(Copy, Clone, Debug, PartialEq)]
struct FragmentClass {
    /// Which secondary-structure state this fragment class produces
    /// for its central residue.
    state: SecondaryStructure,
    /// The class label (used for diagnostics / documentation).
    #[allow(dead_code)]
    label: &'static str,
    /// Per-position `(φ_mean, ψ_mean, ω_mean)` in degrees. `len()` is
    /// the fragment length (one tuple per fragment residue).
    motifs: &'static [(f64, f64, f64)],
    /// Per-position `(σ_φ, σ_ψ)` one-σ spreads in degrees — the
    /// published PDB Ramachandran-cluster spreads (Lovell 2003).
    spreads: &'static [(f64, f64)],
}

/// The PDB-derived-style curated fragment classes. Each class's
/// `motifs` are the published per-residue (φ, ψ, ω) means; the
/// `spreads` are the published one-σ Ramachandran-cluster widths
/// (Lovell *et al.* 2003 — the canonical analysis), used to generate
/// the realistic per-class realisations the assembler needs.
///
/// **Classes encoded:** α-helix interior, α-helix N-cap, α-helix
/// C-cap, β-strand interior, β-strand edge, β-turn Type I, β-turn
/// Type II, β-turn Type I', β-turn Type II', γ-turn (classic),
/// γ-turn (inverse), polyproline-II, 3₁₀-helix interior, π-helix.
const FRAGMENT_CLASSES: &[FragmentClass] = &[
    // ── Helical motifs ────────────────────────────────────────────
    // Right-handed α-helix interior — Lovell 2003 mean (-63.0, -42.5).
    // Trans ω = 180°. Tight one-σ spread of ~7°, the published
    // "general-case" α-helix cluster width.
    FragmentClass {
        state: SecondaryStructure::Helix,
        label: "alpha-helix interior",
        motifs: &[
            (-63.0, -42.5, 180.0),
            (-63.0, -42.5, 180.0),
            (-63.0, -42.5, 180.0),
        ],
        spreads: &[(7.0, 8.0), (7.0, 8.0), (7.0, 8.0)],
    },
    // α-helix N-cap — Aurora & Rose 1998: the *first* residue often
    // sits in a more open (ψ ≈ +130-150°) "Ncap" geometry, then the
    // helix proper begins. The last residue here is a normal helix
    // residue.
    FragmentClass {
        state: SecondaryStructure::Helix,
        label: "alpha-helix N-cap",
        motifs: &[
            (-95.0, 140.0, 180.0), // cap (loop-like)
            (-63.0, -42.5, 180.0), // helix begins
            (-63.0, -42.5, 180.0),
        ],
        spreads: &[(15.0, 20.0), (7.0, 8.0), (7.0, 8.0)],
    },
    // α-helix C-cap — Aurora & Rose 1998: the C-cap residue (often
    // Gly) flips into a "Schellman" left-handed (φ ≈ +95°) geometry
    // at the C-terminus, terminating the helix.
    FragmentClass {
        state: SecondaryStructure::Helix,
        label: "alpha-helix C-cap (Schellman)",
        motifs: &[
            (-63.0, -42.5, 180.0),
            (-63.0, -42.5, 180.0),
            (95.0, 0.0, 180.0), // C-cap left-handed (Schellman)
        ],
        spreads: &[(7.0, 8.0), (7.0, 8.0), (15.0, 20.0)],
    },
    // 3₁₀-helix interior — Lovell 2003: tighter helix, (φ, ψ) ≈ (-57, -30).
    FragmentClass {
        state: SecondaryStructure::Helix,
        label: "3-10 helix interior",
        motifs: &[
            (-57.0, -30.0, 180.0),
            (-57.0, -30.0, 180.0),
            (-57.0, -30.0, 180.0),
        ],
        spreads: &[(10.0, 12.0), (10.0, 12.0), (10.0, 12.0)],
    },
    // π-helix interior — rarer, wider helix, (φ, ψ) ≈ (-76, -41).
    // The PDB cluster spread is broader than α.
    FragmentClass {
        state: SecondaryStructure::Helix,
        label: "pi-helix interior",
        motifs: &[
            (-76.0, -41.0, 180.0),
            (-76.0, -41.0, 180.0),
            (-76.0, -41.0, 180.0),
        ],
        spreads: &[(12.0, 14.0), (12.0, 14.0), (12.0, 14.0)],
    },
    // ── β-sheet motifs ─────────────────────────────────────────────
    // β-strand interior (parallel/antiparallel-general) — Lovell 2003:
    // mean ≈ (-120, +130), spread ~15° in φ, ~20° in ψ.
    FragmentClass {
        state: SecondaryStructure::Strand,
        label: "beta-strand interior",
        motifs: &[
            (-120.0, 130.0, 180.0),
            (-120.0, 130.0, 180.0),
            (-120.0, 130.0, 180.0),
        ],
        spreads: &[(15.0, 20.0), (15.0, 20.0), (15.0, 20.0)],
    },
    // β-strand edge — Hovmöller 2002 antiparallel-edge geometry,
    // slightly more negative ψ at the edge residue.
    FragmentClass {
        state: SecondaryStructure::Strand,
        label: "beta-strand edge (antiparallel)",
        motifs: &[
            (-140.0, 155.0, 180.0),
            (-130.0, 145.0, 180.0),
            (-120.0, 130.0, 180.0),
        ],
        spreads: &[(18.0, 22.0), (15.0, 20.0), (15.0, 20.0)],
    },
    // ── β-turn motifs (Lewis-Momany-Scheraga / Hutchinson-Thornton) ──
    // Type I β-turn — the most common turn class (~40 % of all turns).
    // Residue i+1: φ ≈ -60, ψ ≈ -30. Residue i+2: φ ≈ -90, ψ ≈ 0.
    // Padded by one residue on each side (a turn fragment is 4 long
    // by convention; we pad to length 3 for the standard 3-mer slot
    // by emitting the central two turn residues + 1 flank).
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "beta-turn Type I (residues 2-4)",
        motifs: &[
            (-60.0, -30.0, 180.0),
            (-90.0, 0.0, 180.0),
            (-90.0, 0.0, 180.0),
        ],
        spreads: &[(15.0, 15.0), (15.0, 20.0), (15.0, 20.0)],
    },
    // Type II β-turn — ~15 % of turns. Residue i+1: φ ≈ -60, ψ ≈
    // +120. Residue i+2: φ ≈ +80, ψ ≈ 0 (commonly glycine).
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "beta-turn Type II (residues 2-4)",
        motifs: &[
            (-60.0, 120.0, 180.0),
            (80.0, 0.0, 180.0),
            (-90.0, 0.0, 180.0),
        ],
        spreads: &[(15.0, 20.0), (15.0, 20.0), (15.0, 20.0)],
    },
    // Type I' β-turn (mirror of Type I) — common in β-hairpins.
    // Residue i+1: φ ≈ +60, ψ ≈ +30. Residue i+2: φ ≈ +90, ψ ≈ 0.
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "beta-turn Type I' (residues 2-4)",
        motifs: &[
            (60.0, 30.0, 180.0),
            (90.0, 0.0, 180.0),
            (-90.0, 0.0, 180.0),
        ],
        spreads: &[(15.0, 15.0), (15.0, 20.0), (15.0, 20.0)],
    },
    // Type II' β-turn (mirror of Type II) — common in β-hairpins
    // (especially with Gly at i+1).
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "beta-turn Type II' (residues 2-4)",
        motifs: &[
            (60.0, -120.0, 180.0),
            (-80.0, 0.0, 180.0),
            (-90.0, 0.0, 180.0),
        ],
        spreads: &[(15.0, 20.0), (15.0, 20.0), (15.0, 20.0)],
    },
    // ── γ-turn motifs ─────────────────────────────────────────────
    // Classic γ-turn (Némethy-Printz 1972) — residue i+1: φ ≈ +75, ψ ≈ -65.
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "gamma-turn classic",
        motifs: &[
            (-120.0, 130.0, 180.0),
            (75.0, -65.0, 180.0),
            (-120.0, 130.0, 180.0),
        ],
        spreads: &[(15.0, 20.0), (15.0, 15.0), (15.0, 20.0)],
    },
    // Inverse γ-turn (more common in proteins) — residue i+1: φ ≈ -75, ψ ≈ +65.
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "gamma-turn inverse",
        motifs: &[
            (-90.0, 0.0, 180.0),
            (-75.0, 65.0, 180.0),
            (-90.0, 0.0, 180.0),
        ],
        spreads: &[(15.0, 20.0), (15.0, 15.0), (15.0, 20.0)],
    },
    // ── Polyproline-II (left-handed extended) ────────────────────
    // Adzhubei-Sternberg PPII region — mean (-75, +145), narrow spread.
    FragmentClass {
        state: SecondaryStructure::Coil,
        label: "polyproline-II",
        motifs: &[
            (-75.0, 145.0, 180.0),
            (-75.0, 145.0, 180.0),
            (-75.0, 145.0, 180.0),
        ],
        spreads: &[(12.0, 15.0), (12.0, 15.0), (12.0, 15.0)],
    },
];

/// One backbone fragment: a run of `(φ, ψ, ω)` dihedral tuples plus
/// the class label it was drawn from.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Fragment {
    /// Sequence-window start index this fragment is keyed to.
    pub start: usize,
    /// `(phi, psi)` dihedral pairs in degrees, one per fragment
    /// residue. `len()` is the fragment length.
    pub torsions: Vec<(f64, f64)>,
    /// Per-residue ω peptide-bond dihedrals (degrees). Mostly 180°
    /// (trans), but some published motifs and *cis*-proline insertions
    /// drift; we record the per-residue value the class published.
    pub omega: Vec<f64>,
    /// Class label of the PDB-curated source motif. Useful for
    /// diagnostics and for the assembler's class-biased move logic.
    pub class_label: String,
}

impl Fragment {
    /// Fragment length (number of residues).
    pub fn len(&self) -> usize {
        self.torsions.len()
    }

    /// `true` when the fragment carries no residues.
    pub fn is_empty(&self) -> bool {
        self.torsions.is_empty()
    }
}

/// A fragment library: per sequence-window position, a set of
/// candidate backbone fragments.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FragmentLibrary {
    /// Fragment length (residues) — typically 3 or 9.
    pub fragment_length: usize,
    /// Target sequence length the library was built for.
    pub sequence_length: usize,
    /// `fragments[w]` are the candidate fragments at window start `w`,
    /// for `w` in `0..=sequence_length - fragment_length`.
    pub fragments: Vec<Vec<Fragment>>,
}

impl FragmentLibrary {
    /// Total fragment count across all windows.
    pub fn total_fragments(&self) -> usize {
        self.fragments.iter().map(|f| f.len()).sum()
    }

    /// The candidate fragments at a given window start, or `None` if
    /// the start index is out of range.
    pub fn at(&self, window: usize) -> Option<&[Fragment]> {
        self.fragments.get(window).map(|v| v.as_slice())
    }
}

/// Builds a PDB-curated-style fragment library for a sequence.
///
/// For each sliding window, the predicted secondary structure of the
/// window's central residue selects fragment *classes* (with a bias
/// toward the predicted state). Each class produces multiple
/// realisations by perturbing its published mean (φ, ψ) by the
/// published one-σ spread. The result is a per-window pool of
/// realistic, class-labelled fragments.
///
/// `fragment_length` is the fragment size (3 and 9 are the classical
/// choices). `per_window` is how many candidate fragments to emit per
/// window. Secondary structure is predicted internally and used to
/// bias the class mix.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty sequence, a
/// non-positive fragment length, a `fragment_length` exceeding the
/// sequence, or `per_window == 0`.
pub fn build_fragment_library(
    sequence: &str,
    fragment_length: usize,
    per_window: usize,
) -> Result<FragmentLibrary> {
    let sequence = sequence.trim();
    if sequence.is_empty() {
        return Err(StructPredictError::invalid("sequence", "empty"));
    }
    if fragment_length == 0 {
        return Err(StructPredictError::invalid(
            "fragment_length",
            "must be at least 1",
        ));
    }
    if fragment_length > sequence.len() {
        return Err(StructPredictError::invalid(
            "fragment_length",
            format!(
                "{} exceeds sequence length {}",
                fragment_length,
                sequence.len()
            ),
        ));
    }
    if per_window == 0 {
        return Err(StructPredictError::invalid(
            "per_window",
            "must be at least 1",
        ));
    }

    let ss = predict_secondary_structure(sequence)?;
    let n = sequence.len();
    let n_windows = n - fragment_length + 1;
    let mut fragments = Vec::with_capacity(n_windows);

    // A tiny deterministic LCG so the library is reproducible without
    // pulling in the MD RNG here (the assembler owns the seeded RNG;
    // the library only needs deterministic per-class sampling).
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;

    for w in 0..n_windows {
        // The window's *central* residue's predicted SS — used to
        // bias the class mix. A helical centre prefers helical
        // classes (interior, N-cap, C-cap, 3₁₀, π); a strand centre
        // prefers strand classes; a coil centre prefers turns / PPII.
        let centre = w + fragment_length / 2;
        let want = ss.states[centre];
        let mut window_frags = Vec::with_capacity(per_window);
        for _ in 0..per_window {
            // 80 % of the time honour the predicted state; 20 % draw
            // any class (the conformational diversity the assembler
            // needs to escape minima).
            let honour = (next_u32(&mut state) % 100) < 80;
            let class = pick_class(want, honour, &mut state);
            // Build the fragment by tiling the class's motif to the
            // requested fragment_length (most published classes are
            // 3 residues long; a 9-mer library tiles 3 of them).
            let mut torsions = Vec::with_capacity(fragment_length);
            let mut omega = Vec::with_capacity(fragment_length);
            for k in 0..fragment_length {
                let m = k % class.motifs.len();
                let (phi_mean, psi_mean, om_mean) = class.motifs[m];
                let (sphi, spsi) = class.spreads[m];
                torsions.push((
                    phi_mean + uniform(&mut state, sphi),
                    psi_mean + uniform(&mut state, spsi),
                ));
                omega.push(om_mean);
            }
            window_frags.push(Fragment {
                start: w,
                torsions,
                omega,
                class_label: class.label.to_owned(),
            });
        }
        fragments.push(window_frags);
    }

    Ok(FragmentLibrary {
        fragment_length,
        sequence_length: n,
        fragments,
    })
}

/// One step of the tiny deterministic LCG: advances `state` and
/// returns the next 32-bit word. Used to drive the per-class spread
/// sampling without dragging in the MD RNG (the assembler owns that).
fn next_u32(state: &mut u64) -> u32 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*state >> 33) as u32
}

/// A uniform draw in `[-scale, +scale]` from the LCG. Used to sample
/// each per-residue (φ, ψ) inside the published one-σ Ramachandran
/// cluster spread of its class — a Gaussian draw would be slightly
/// more faithful, but uniform-in-(-σ,+σ) is the standard cheap
/// fragment-sampling convention.
fn uniform(state: &mut u64, scale: f64) -> f64 {
    let r = (next_u32(state) as f64) / (u32::MAX as f64);
    (2.0 * r - 1.0) * scale
}

/// Picks a fragment class: when `honour` is set, a class whose central
/// residue produces `want`; otherwise any class. The class probability
/// is uniform within each pool (the published Rosetta frequencies are
/// position-specific PDB statistics — not encoded here; a curated
/// uniform-within-class draw is the honest v1).
fn pick_class(
    want: SecondaryStructure,
    honour: bool,
    state: &mut u64,
) -> &'static FragmentClass {
    if honour {
        let matching: Vec<&FragmentClass> =
            FRAGMENT_CLASSES.iter().filter(|c| c.state == want).collect();
        if !matching.is_empty() {
            let idx = (next_u32(state) as usize) % matching.len();
            return matching[idx];
        }
    }
    &FRAGMENT_CLASSES[(next_u32(state) as usize) % FRAGMENT_CLASSES.len()]
}

/// Returns the count of fragment classes in the curated library —
/// useful for diagnostics. Currently 14: α-interior / N-cap / C-cap,
/// 3₁₀, π, β-interior / edge, β-turn I / II / I' / II', γ-turn
/// classic / inverse, PPII.
pub fn fragment_class_count() -> usize {
    FRAGMENT_CLASSES.len()
}

/// Returns the labels of every curated fragment class — useful for
/// telemetry and documentation that wants to enumerate coverage.
pub fn fragment_class_labels() -> Vec<&'static str> {
    FRAGMENT_CLASSES.iter().map(|c| c.label).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_has_a_window_per_position() {
        let seq = "ACDEFGHIKLMNPQRST"; // 17 residues
        let lib = build_fragment_library(seq, 3, 25).expect("lib");
        assert_eq!(lib.fragments.len(), 17 - 3 + 1);
        assert_eq!(lib.fragment_length, 3);
        for window in &lib.fragments {
            assert_eq!(window.len(), 25);
            for frag in window {
                assert_eq!(frag.len(), 3);
                assert_eq!(frag.omega.len(), 3);
                // every fragment carries its source-class label
                assert!(!frag.class_label.is_empty());
            }
        }
    }

    #[test]
    fn helical_sequence_biases_helical_fragments() {
        // Strong helix formers → most fragments helical.
        let seq = "EEEEAAAALLLLEEEEAAAALLLL";
        let lib = build_fragment_library(seq, 9, 50).expect("lib");
        // Count fragment residues in the helical (φ,ψ) basin (α + 3₁₀ + π).
        let mut helical = 0usize;
        let mut total = 0usize;
        for window in &lib.fragments {
            for frag in window {
                for &(phi, psi) in &frag.torsions {
                    total += 1;
                    if (-95.0..-40.0).contains(&phi) && (-70.0..0.0).contains(&psi) {
                        helical += 1;
                    }
                }
            }
        }
        let frac = helical as f64 / total as f64;
        assert!(frac > 0.5, "helical fraction {frac}");
    }

    #[test]
    fn fragments_are_sane_per_ss_class_geometry() {
        // The curated library must produce realistic (φ,ψ) per class:
        // α-helix fragments cluster near (-60, -45); β-strand near
        // (-120, +120). Build a long all-Ala (helix-loving) library
        // and a long all-Val (strand-loving) library, count their
        // dominant clusters.
        let helix_lib = build_fragment_library(&"A".repeat(20), 3, 30).expect("a");
        let strand_lib = build_fragment_library(&"V".repeat(20), 3, 30).expect("v");

        let dominant_basin = |lib: &FragmentLibrary| -> (f64, f64, f64) {
            // Sample window 5 (well inside the homopolymer) and count
            // how many residues fall near α vs β.
            let frags = lib.at(5).expect("window 5");
            let mut alpha = 0;
            let mut beta = 0;
            let mut total = 0;
            for f in frags {
                for &(phi, psi) in &f.torsions {
                    total += 1;
                    // α basin: φ in [-95, -40], ψ in [-70, 0].
                    if (-95.0..-40.0).contains(&phi) && (-70.0..0.0).contains(&psi) {
                        alpha += 1;
                    }
                    // β basin: φ in [-160, -85], ψ in [60, 180].
                    if (-160.0..-85.0).contains(&phi) && (60.0..180.0).contains(&psi) {
                        beta += 1;
                    }
                }
            }
            (
                alpha as f64 / total as f64,
                beta as f64 / total as f64,
                total as f64,
            )
        };
        let (a_alpha, a_beta, _) = dominant_basin(&helix_lib);
        let (v_alpha, v_beta, _) = dominant_basin(&strand_lib);
        // The Ala library is dominated by α and barely populates β.
        assert!(a_alpha > a_beta, "Ala α={a_alpha} should beat β={a_beta}");
        // The Val library leans the other way.
        assert!(v_beta > v_alpha, "Val β={v_beta} should beat α={v_alpha}");
    }

    #[test]
    fn fragments_carry_a_realistic_spread_per_class() {
        // A class's per-residue published spread must produce
        // non-identical realisations: two fragments from the same
        // class are NOT bit-identical — the perturbation samples the
        // PDB spread band the library encodes.
        let lib = build_fragment_library(&"A".repeat(20), 3, 100).expect("lib");
        let frags = lib.at(5).expect("frags");
        // Find the dominant class label (most helical fragments will
        // be from the α-helix-interior class); confirm at least two
        // realisations from that class differ in (φ, ψ).
        let mut by_class: std::collections::HashMap<&str, Vec<&Fragment>> =
            std::collections::HashMap::new();
        for f in frags {
            by_class
                .entry(f.class_label.as_str())
                .or_default()
                .push(f);
        }
        let mut found_distinct = false;
        for v in by_class.values() {
            if v.len() < 2 {
                continue;
            }
            for k in 1..v.len() {
                if v[0].torsions != v[k].torsions {
                    found_distinct = true;
                    break;
                }
            }
            if found_distinct {
                break;
            }
        }
        assert!(
            found_distinct,
            "fragments within a class must NOT all be identical — \
             the published spread must drive realistic variation"
        );
    }

    #[test]
    fn fragments_omega_is_trans_for_standard_classes() {
        // Every published class encoded here uses ω = 180° (trans);
        // cis-Pro is left as a documented v1 simplification, but
        // every fragment must report its ω so the assembler / scorer
        // can consume it.
        let lib = build_fragment_library(&"ALAGLY".repeat(3), 3, 5).expect("lib");
        for w in &lib.fragments {
            for f in w {
                for &om in &f.omega {
                    assert!((om - 180.0).abs() < 1e-6, "trans ω = 180, got {om}");
                }
            }
        }
    }

    #[test]
    fn curated_class_count_covers_canonical_classes() {
        // The curated library encodes the standard canonical PDB-mined
        // backbone fragment classes: α (interior/Ncap/Ccap), 3₁₀,
        // π, β (interior/edge), β-turns (I, II, I', II'), γ-turns
        // (classic, inverse), and PPII — 14 classes.
        assert!(
            fragment_class_count() >= 14,
            "expected ≥14 curated classes, got {}",
            fragment_class_count()
        );
        let labels = fragment_class_labels();
        // The label set must include every key turn class — a
        // production-realistic library cannot omit Type I'.
        assert!(labels.iter().any(|l| l.contains("Type I")));
        assert!(labels.iter().any(|l| l.contains("Type II")));
        assert!(labels.iter().any(|l| l.contains("Type I'")));
        assert!(labels.iter().any(|l| l.contains("Type II'")));
        assert!(labels.iter().any(|l| l.contains("gamma-turn")));
        assert!(labels.iter().any(|l| l.contains("polyproline")));
    }

    #[test]
    fn deterministic() {
        let a = build_fragment_library("ACDEFGHIKL", 3, 10).expect("lib");
        let b = build_fragment_library("ACDEFGHIKL", 3, 10).expect("lib");
        assert_eq!(a, b);
    }

    #[test]
    fn bad_arguments_rejected() {
        assert!(build_fragment_library("", 3, 5).is_err());
        assert!(build_fragment_library("ACDEF", 0, 5).is_err());
        assert!(build_fragment_library("ACDEF", 99, 5).is_err());
        assert!(build_fragment_library("ACDEF", 3, 0).is_err());
    }
}
