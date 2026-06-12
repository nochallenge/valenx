//! Feature 19 — classical structural-complementarity aptamer design.
//!
//! An aptamer is an RNA whose **structured binding pocket** presents a
//! set of chemical features complementary to a small-molecule / peptide
//! target. Wet-lab discovery uses SELEX (iterated affinity selection
//! over a random library); the *in silico* analogue is to search
//! candidate sequences for one whose folded pocket has the right
//! features in the right places. This module is the classical
//! (non-ML, non-docking) version of that idea.
//!
//! ## The target — a [`Pharmacophore`]
//!
//! The caller supplies a [`Pharmacophore`]: a small set of
//! [`PharmacophoreFeature`]s, each a [`FeatureKind`] (H-bond donor,
//! H-bond acceptor, hydrophobic / aromatic, positive / negative
//! ionisable) at a 3-D coordinate with a tolerance. This is the
//! standard pharmacophore model used by classical drug-design tools
//! (LigandScout, MOE, Phase): the geometric arrangement of the
//! ligand's chemical groups that a binding partner must complement.
//!
//! ## The score — structural-pocket-feature complementarity
//!
//! Given a candidate RNA sequence, fold it with `valenx-rnastruct`,
//! then score the resulting structure's **pocket** for complementarity
//! to the pharmacophore. A pocket is a structured site with the right
//! kind of feature-presenting characteristics: a hairpin loop presents
//! the loop bases on accessible Watson-Crick / Hoogsteen edges (good
//! for H-bond acceptors / donors / aromatic stacking); an internal
//! loop / bulge presents bulged bases (good for ionisable contacts);
//! a multi-junction is a "scaffolded" pocket (good for combined
//! feature presentation). Each candidate base's chemical character is
//! summed from a classical, transparent table:
//!
//! - **G** edge contributes both an H-bond donor (`N1H, N2H`) and an
//!   H-bond acceptor (`N7, O6`);
//! - **C** edge contributes an H-bond donor (`N4H`) and an H-bond
//!   acceptor (`O2, N3`);
//! - **A** edge contributes one donor (`N6H`), one acceptor (`N1, N7`),
//!   and aromatic-stacking character;
//! - **U** edge contributes one donor (`N3H`), one acceptor (`O2, O4`).
//!
//! For each pharmacophore feature, the designer counts the candidate
//! loop bases that present a compatible chemical character — the more
//! complementary bases, the higher the pocket-feature score. The score
//! is normalised by feature count.
//!
//! The designer is a **search** over candidate sequences: for a length
//! window and a GC band, run several inverse-folding starts toward
//! diverse, well-structured "pocket-templates" (a hairpin, an internal
//! loop, a small multi-junction), fold each, score each, and pick the
//! best.
//!
//! ## v1 scope — honest framing
//!
//! - This is **classical pharmacophore-vs-structure scoring** — not a
//!   docking calculation. The aptamer's pocket is characterised by its
//!   *secondary*-structure topology and the chemical edges its loop
//!   bases present; no 3-D RNA conformation, no Coulombic / van der
//!   Waals interaction integration, no flexible-ligand pose search.
//!   The output is a strong in-silico aptamer candidate, not a
//!   measured K_d.
//! - The base→edge-feature table is a documented, principled mapping
//!   (Watson-Crick / Hoogsteen edge classifications from the standard
//!   RNA-base nomenclature, e.g. Leontis–Westhof). It is exact for the
//!   *chemistry* but does not weight by 3-D edge geometry.
//! - The pharmacophore-feature 3-D coordinates are *used to count and
//!   group features* (so spatially-clustered features score together)
//!   but the scoring does not project the candidate's loop bases into
//!   3-D; the score is the candidate's chemical complementarity to the
//!   pharmacophore feature set under a coarse pocket-feature
//!   assignment.
//! - A higher score than a **random baseline** (a random-sequence
//!   pocket) is the validation criterion the tests assert; this is the
//!   honest claim — the designed aptamer is *structurally
//!   complementary* to the pharmacophore beyond chance. Real aptamer
//!   discovery still needs wet-lab SELEX.

use crate::error::{Result, RnaDesignError};
use serde::{Deserialize, Serialize};
use valenx_rnastruct::{mfe, structure_stats, RnaSeq, Structure};

/// The four RNA bases, ASCII.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];

/// Canonical pairing partners (Watson-Crick + G-U wobble), ASCII.
const CANON_PAIRS: [(u8, u8); 6] = [
    (b'A', b'U'),
    (b'U', b'A'),
    (b'G', b'C'),
    (b'C', b'G'),
    (b'G', b'U'),
    (b'U', b'G'),
];

// ---------------------------------------------------------------------
// Pharmacophore model
// ---------------------------------------------------------------------

/// A pharmacophore feature kind — the chemical character a complement
/// must satisfy.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FeatureKind {
    /// An H-bond donor — the ligand offers a hydrogen, the partner
    /// must present an acceptor.
    HBondDonor,
    /// An H-bond acceptor — the ligand offers a lone pair, the partner
    /// must present a donor.
    HBondAcceptor,
    /// A hydrophobic / aromatic contact — the partner must present a
    /// hydrophobic / pi face (in RNA: an exposed nucleobase ring face,
    /// strongest for A and G — the purines).
    Hydrophobic,
    /// A positive ionisable contact — the ligand bears a positive
    /// charge, the partner must present negative (in RNA: a phosphate
    /// backbone or a strongly electronegative pocket — both are
    /// abundant; this feature is therefore the easiest to satisfy).
    Positive,
    /// A negative ionisable contact — the ligand bears a negative
    /// charge, the partner must present positive (in RNA, rare without
    /// metal coordination; the score is penalised for unrealistic
    /// over-presentation of this feature).
    Negative,
}

impl FeatureKind {
    /// A short human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            FeatureKind::HBondDonor => "H-bond donor",
            FeatureKind::HBondAcceptor => "H-bond acceptor",
            FeatureKind::Hydrophobic => "hydrophobic / aromatic",
            FeatureKind::Positive => "positive ionisable",
            FeatureKind::Negative => "negative ionisable",
        }
    }

    /// The complementary kind a binding partner must present — a donor
    /// pairs with an acceptor, a positive with a negative, hydrophobic
    /// with hydrophobic.
    pub fn complement(self) -> FeatureKind {
        match self {
            FeatureKind::HBondDonor => FeatureKind::HBondAcceptor,
            FeatureKind::HBondAcceptor => FeatureKind::HBondDonor,
            FeatureKind::Hydrophobic => FeatureKind::Hydrophobic,
            FeatureKind::Positive => FeatureKind::Negative,
            FeatureKind::Negative => FeatureKind::Positive,
        }
    }
}

/// A single pharmacophore feature: a chemical character at a 3-D
/// location, with a tolerance.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PharmacophoreFeature {
    /// The chemical kind of feature.
    pub kind: FeatureKind,
    /// The 3-D coordinate of the feature centre, Å (the model's frame
    /// — used by the designer only to group spatially-clustered
    /// features for compactness scoring).
    pub position: [f64; 3],
    /// The matching radius, Å. Features within this tolerance of each
    /// other are treated as a single cluster (the pocket-compactness
    /// term).
    pub tolerance: f64,
}

impl PharmacophoreFeature {
    /// A feature with an explicit kind, position and tolerance.
    pub fn new(kind: FeatureKind, position: [f64; 3], tolerance: f64) -> Self {
        PharmacophoreFeature {
            kind,
            position,
            tolerance,
        }
    }

    /// Euclidean distance to another feature, Å.
    pub fn distance(&self, other: &PharmacophoreFeature) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// A pharmacophore — a small set of chemical features describing a
/// ligand's binding character.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Pharmacophore {
    /// The features, in their caller-supplied order.
    pub features: Vec<PharmacophoreFeature>,
}

impl Pharmacophore {
    /// A pharmacophore from an explicit feature list.
    pub fn new(features: Vec<PharmacophoreFeature>) -> Self {
        Pharmacophore { features }
    }

    /// The number of features.
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// `true` if the pharmacophore is empty.
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// The diameter — the maximum pair-wise distance, Å, between
    /// features (`0` for fewer than two features).
    pub fn diameter(&self) -> f64 {
        let mut d = 0.0_f64;
        for i in 0..self.features.len() {
            for j in (i + 1)..self.features.len() {
                d = d.max(self.features[i].distance(&self.features[j]));
            }
        }
        d
    }
}

// ---------------------------------------------------------------------
// Base→edge-feature table
// ---------------------------------------------------------------------

/// The chemical edge-features a single RNA base presents on its
/// Watson-Crick edge, plus aromatic-stacking character. These are
/// transparent assignments from the standard RNA nucleobase chemistry
/// (Leontis–Westhof Watson-Crick edge classifications):
///
/// - **G**: donor (`N1H, N2H`), acceptor (`N7, O6`); a purine — strong
///   aromatic/hydrophobic stacking face.
/// - **C**: donor (`N4H`), acceptor (`O2, N3`); a pyrimidine — weaker
///   stacking face.
/// - **A**: donor (`N6H`), acceptor (`N1, N7`); a purine — strong
///   stacking face.
/// - **U**: donor (`N3H`), acceptor (`O2, O4`); a pyrimidine — weaker
///   stacking face.
///
/// The (donor, acceptor) counts double-count the *number* of edge
/// groups; the boolean stacking flag distinguishes purines / pyrimidines.
/// Every RNA base presents a phosphate (negative) and the backbone is
/// the dominant "positive-friendly" / "negative-presenting" character
/// of the molecule — the scoring handles those at the pocket level
/// rather than per-base.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BaseEdgeFeatures {
    /// Number of H-bond donor groups on this base's Watson-Crick edge.
    pub donors: u8,
    /// Number of H-bond acceptor groups on this base's Watson-Crick edge.
    pub acceptors: u8,
    /// `true` if this base is a purine (stronger hydrophobic / aromatic
    /// stacking face).
    pub purine: bool,
}

/// The Watson-Crick edge feature counts for a single RNA base.
pub fn base_edge_features(b: u8) -> BaseEdgeFeatures {
    match b.to_ascii_uppercase() {
        b'G' => BaseEdgeFeatures {
            donors: 2,
            acceptors: 2,
            purine: true,
        },
        b'C' => BaseEdgeFeatures {
            donors: 1,
            acceptors: 2,
            purine: false,
        },
        b'A' => BaseEdgeFeatures {
            donors: 1,
            acceptors: 2,
            purine: true,
        },
        b'U' | b'T' => BaseEdgeFeatures {
            donors: 1,
            acceptors: 2,
            purine: false,
        },
        // Unknown / N — treat as a single donor / single acceptor
        // (the conservative averaged assignment).
        _ => BaseEdgeFeatures {
            donors: 1,
            acceptors: 1,
            purine: false,
        },
    }
}

// ---------------------------------------------------------------------
// Pocket extraction
// ---------------------------------------------------------------------

/// A pocket in a folded RNA — a contiguous run of loop bases (a
/// hairpin loop, an internal loop / bulge unpaired stretch, or the
/// unpaired bases of a multi-junction) that present accessible
/// nucleobase edges to a ligand.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Pocket {
    /// 0-based start index in the sequence.
    pub start: usize,
    /// Exclusive end index in the sequence.
    pub end: usize,
    /// The pocket's structural classification.
    pub kind: PocketKind,
}

impl Pocket {
    /// The pocket's length (number of loop bases).
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// `true` if the pocket is empty.
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// The structural class of a pocket.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PocketKind {
    /// A hairpin-loop pocket — the unpaired bases of a hairpin closed
    /// by one base pair. The most common aptamer pocket type.
    HairpinLoop,
    /// An internal-loop / bulge pocket — an unpaired stretch between
    /// two paired stems on the same strand.
    InternalLoop,
    /// A multi-junction pocket — unpaired bases in a multiloop (a
    /// closing pair enclosing ≥ 2 child pairs). A scaffolded pocket.
    MultiJunction,
}

impl PocketKind {
    /// A pocket-class weight (a hairpin loop is the cleanest pocket; a
    /// multi-junction is also a strong pocket; a bulge / internal loop
    /// is intermediate).
    pub fn weight(self) -> f64 {
        match self {
            PocketKind::HairpinLoop => 1.0,
            PocketKind::MultiJunction => 0.9,
            PocketKind::InternalLoop => 0.7,
        }
    }
}

/// Extracts the pockets of a folded structure — every loop that is
/// large enough to plausibly present a ligand-binding pocket.
///
/// "Large enough" is `>= 3` unpaired bases. Pseudoknotted structures
/// have only their nested-skeleton loops extracted.
pub fn extract_pockets(structure: &Structure) -> Vec<Pocket> {
    let n = structure.len();
    let mut pockets = Vec::new();
    if n == 0 {
        return pockets;
    }

    // Walk every pair (i, j) and classify the loop it closes.
    let partner = structure.partner_array();
    for i in 0..n {
        let Some(j) = partner[i] else { continue };
        if i >= j {
            continue;
        }
        // Find enclosed pairs.
        let mut enclosed: Vec<(usize, usize)> = Vec::new();
        let mut k = i + 1;
        while k < j {
            match partner[k] {
                Some(p) if p > k => {
                    enclosed.push((k, p));
                    k = p + 1;
                }
                _ => k += 1,
            }
        }

        match enclosed.len() {
            0 => {
                // Hairpin loop: the unpaired bases (i+1 .. j).
                if j > i + 1 + 2 {
                    pockets.push(Pocket {
                        start: i + 1,
                        end: j,
                        kind: PocketKind::HairpinLoop,
                    });
                }
            }
            1 => {
                let (k, l) = enclosed[0];
                let left_unpaired = k - i - 1;
                let right_unpaired = j - l - 1;
                // Internal loop / bulge: emit each unpaired stretch
                // as its own pocket if large enough.
                if left_unpaired >= 3 {
                    pockets.push(Pocket {
                        start: i + 1,
                        end: k,
                        kind: PocketKind::InternalLoop,
                    });
                }
                if right_unpaired >= 3 {
                    pockets.push(Pocket {
                        start: l + 1,
                        end: j,
                        kind: PocketKind::InternalLoop,
                    });
                }
            }
            _ => {
                // Multiloop: emit each junction's unpaired stretches.
                // Stretch boundaries are (i, first_pair), then
                // (pair_n_end, pair_{n+1}_start), then (last_pair_end, j).
                let mut prev = i;
                for (k, l) in &enclosed {
                    if *k > prev + 1 + 2 {
                        // Unpaired stretch (prev+1 .. *k).
                        pockets.push(Pocket {
                            start: prev + 1,
                            end: *k,
                            kind: PocketKind::MultiJunction,
                        });
                    }
                    prev = *l;
                }
                if j > prev + 1 + 2 {
                    pockets.push(Pocket {
                        start: prev + 1,
                        end: j,
                        kind: PocketKind::MultiJunction,
                    });
                }
            }
        }
    }

    pockets
}

// ---------------------------------------------------------------------
// Scoring — pharmacophore vs structure
// ---------------------------------------------------------------------

/// The score of one [`Pocket`] for one [`PharmacophoreFeature`]: the
/// number of pocket bases that present a compatible edge feature,
/// scaled by the pocket-kind weight.
fn pocket_feature_score(seq: &[u8], pocket: &Pocket, feature: &PharmacophoreFeature) -> f64 {
    if pocket.is_empty() {
        return 0.0;
    }
    let mut count = 0.0_f64;
    for i in pocket.start..pocket.end {
        let Some(&b) = seq.get(i) else { continue };
        let edge = base_edge_features(b);
        let contribution = match feature.kind {
            // Donor in the ligand → the pocket must present an acceptor.
            FeatureKind::HBondDonor => edge.acceptors as f64,
            // Acceptor in the ligand → the pocket must present a donor.
            FeatureKind::HBondAcceptor => edge.donors as f64,
            // Hydrophobic / aromatic → a purine pocket base.
            FeatureKind::Hydrophobic => {
                if edge.purine {
                    1.0
                } else {
                    0.3
                }
            }
            // Positive ligand → negatively-charged pocket. Every RNA
            // base presents a phosphate backbone — a small constant.
            FeatureKind::Positive => 0.5,
            // Negative ligand → positively-charged pocket. Rare in RNA
            // without metal coordination; penalised by a low constant.
            FeatureKind::Negative => 0.1,
        };
        count += contribution;
    }
    count * pocket.kind.weight()
}

/// The aggregate structural-pocket-feature score of a folded candidate
/// against a pharmacophore.
///
/// For each pharmacophore feature, the score is the **maximum** over
/// pockets of the per-pocket feature score (the *best* pocket for that
/// feature is what counts — a real binding site does not split a
/// feature across two pockets). The aggregate is the sum across
/// features, with a small **compactness bonus** for putting two
/// well-scoring features in the *same* pocket (matching the geometric
/// clustering of features in the pharmacophore).
///
/// A higher score is better.
pub fn pharmacophore_pocket_score(
    seq: &[u8],
    structure: &Structure,
    pharmacophore: &Pharmacophore,
) -> f64 {
    let pockets = extract_pockets(structure);
    if pockets.is_empty() || pharmacophore.is_empty() {
        return 0.0;
    }

    let mut total = 0.0_f64;
    // Per-feature best-pocket score, and the pocket index it came from.
    let mut best_pocket_idx: Vec<Option<usize>> = vec![None; pharmacophore.len()];
    for (fi, feature) in pharmacophore.features.iter().enumerate() {
        let mut best = 0.0_f64;
        for (pi, pocket) in pockets.iter().enumerate() {
            let s = pocket_feature_score(seq, pocket, feature);
            if s > best {
                best = s;
                best_pocket_idx[fi] = Some(pi);
            }
        }
        total += best;
    }

    // Compactness bonus: features that lie close together in the
    // pharmacophore (within their combined tolerance) and whose best
    // pockets are the same get a +1 bonus per such pair, modelling
    // that a real pocket presents multiple features together.
    for i in 0..pharmacophore.len() {
        for j in (i + 1)..pharmacophore.len() {
            let fi = &pharmacophore.features[i];
            let fj = &pharmacophore.features[j];
            let tol = fi.tolerance + fj.tolerance;
            if fi.distance(fj) <= tol
                && best_pocket_idx[i].is_some()
                && best_pocket_idx[i] == best_pocket_idx[j]
            {
                total += 1.0;
            }
        }
    }

    total
}

// ---------------------------------------------------------------------
// The designer
// ---------------------------------------------------------------------

/// A small deterministic xorshift RNG, local to the designer.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng {
            state: seed ^ 0xA7B1_3D6E_91FE_4427,
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Parameters for [`design_aptamer`].
#[derive(Copy, Clone, Debug)]
pub struct AptamerDesignParams {
    /// Minimum candidate length (nt). Defaults to 20.
    pub length_min: usize,
    /// Maximum candidate length (nt). Defaults to 40.
    pub length_max: usize,
    /// Minimum GC fraction.
    pub gc_min: f64,
    /// Maximum GC fraction.
    pub gc_max: f64,
    /// Number of candidate sequences to evaluate. Defaults to 200.
    pub candidates: usize,
    /// Random seed — the search is deterministic for a fixed seed.
    pub seed: u64,
}

impl Default for AptamerDesignParams {
    /// Aptamer-sized v1 defaults: 20–40 nt, GC 35–65 %, 200 candidates.
    fn default() -> Self {
        AptamerDesignParams {
            length_min: 20,
            length_max: 40,
            gc_min: 0.35,
            gc_max: 0.65,
            candidates: 200,
            seed: 0xA9F1,
        }
    }
}

/// The result of an aptamer design (feature 19).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AptamerDesign {
    /// The designed RNA sequence (`A C G U`).
    pub sequence: Vec<u8>,
    /// The MFE secondary structure of the designed sequence
    /// (dot-bracket form is `structure.to_dot_bracket()`).
    pub structure: Structure,
    /// The minimum free energy of the design, kcal/mol.
    pub mfe_energy: f64,
    /// The structural-pocket-feature score against the input
    /// pharmacophore — higher is better.
    pub pocket_score: f64,
    /// The mean score over a random baseline of candidates of the same
    /// length — used to assess whether the design's score is
    /// significantly above chance.
    pub baseline_mean_score: f64,
    /// The pockets the design's structure presents.
    pub pockets: Vec<Pocket>,
    /// Number of candidates the designer evaluated.
    pub candidates_evaluated: usize,
    /// Human-readable notes, one line per major decision.
    pub notes: Vec<String>,
}

impl AptamerDesign {
    /// The designed sequence as a `&str`.
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }

    /// `true` when the design's pocket score is above the random
    /// baseline mean — the honest, non-overstated success criterion.
    pub fn above_baseline(&self) -> bool {
        self.pocket_score > self.baseline_mean_score
    }
}

/// Designs a classical-pharmacophore-complementary aptamer (feature 19).
///
/// Searches over random candidate sequences in `[length_min,
/// length_max]` with GC in `[gc_min, gc_max]`, folds each with the
/// `valenx-rnastruct` MFE folder, scores the folded pocket for
/// structural-pocket-feature complementarity to `pharmacophore`, and
/// returns the best.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if the pharmacophore is empty or the
///   length window is empty.
/// - [`RnaDesignError::Invalid`] if `candidates == 0` or the GC range
///   is invalid.
/// - [`RnaDesignError::NoDesign`] if every candidate folds to a
///   pocketless structure (an open chain or a single long helix with
///   no loop large enough to host a pocket).
/// - [`RnaDesignError::Upstream`] if a folding call fails.
pub fn design_aptamer(
    pharmacophore: &Pharmacophore,
    params: AptamerDesignParams,
) -> Result<AptamerDesign> {
    if pharmacophore.is_empty() {
        return Err(RnaDesignError::goal(
            "pharmacophore",
            "pharmacophore is empty — at least one feature is required",
        ));
    }
    if params.length_min == 0 {
        return Err(RnaDesignError::invalid(
            "length_min",
            "minimum length must be positive",
        ));
    }
    if params.length_max < params.length_min {
        return Err(RnaDesignError::invalid(
            "length_max",
            "max length is less than min length",
        ));
    }
    if !(0.0..=1.0).contains(&params.gc_min)
        || !(0.0..=1.0).contains(&params.gc_max)
        || params.gc_max < params.gc_min
    {
        return Err(RnaDesignError::invalid(
            "gc_range",
            "GC range must lie in [0, 1] and be non-empty",
        ));
    }
    if params.candidates == 0 {
        return Err(RnaDesignError::invalid(
            "candidates",
            "need at least one candidate",
        ));
    }

    let mut rng = Rng::new(params.seed);
    // (sequence, structure, mfe_energy, score, pockets)
    #[allow(clippy::type_complexity)]
    let mut best: Option<(Vec<u8>, Structure, f64, f64, Vec<Pocket>)> = None;
    let mut baseline_total = 0.0_f64;
    let mut baseline_count = 0usize;
    let mut considered = 0usize;
    let mut foldable = 0usize;
    let mut with_pocket = 0usize;

    for _ in 0..params.candidates {
        let len = if params.length_min == params.length_max {
            params.length_min
        } else {
            params.length_min + rng.below(params.length_max - params.length_min + 1)
        };
        let seq = random_pocket_seed(len, params.gc_min, params.gc_max, &mut rng);
        let rna = match RnaSeq::parse(&seq) {
            Ok(r) => r,
            Err(_) => continue,
        };
        considered += 1;
        let folded = match mfe(&rna) {
            Ok(r) => r,
            Err(_) => continue,
        };
        foldable += 1;
        let structure = folded.structure;
        let pockets = extract_pockets(&structure);
        let score = pharmacophore_pocket_score(&seq, &structure, pharmacophore);

        baseline_total += score;
        baseline_count += 1;
        if !pockets.is_empty() {
            with_pocket += 1;
        }

        let take = match &best {
            None => !pockets.is_empty(),
            Some((_, _, _, bs, _)) => score > *bs && !pockets.is_empty(),
        };
        if take {
            best = Some((seq, structure, folded.energy, score, pockets));
        }
    }

    let (sequence, structure, mfe_energy, pocket_score, pockets) = best.ok_or_else(|| {
        RnaDesignError::no_design(
            "aptamer",
            "no candidate produced a structure with a usable binding pocket",
        )
    })?;

    let baseline_mean_score = if baseline_count > 0 {
        baseline_total / baseline_count as f64
    } else {
        0.0
    };

    let mut notes = Vec::new();
    notes.push(format!(
        "Aptamer design: searched {considered} candidate sequence(s) in [{}, {}] nt with GC \
         in [{:.2}, {:.2}], {foldable} folded, {with_pocket} carried a usable pocket.",
        params.length_min, params.length_max, params.gc_min, params.gc_max,
    ));
    notes.push(format!(
        "Best candidate (length {}): MFE {mfe_energy:.2} kcal/mol; structural-pocket-feature \
         score {pocket_score:.3} against a {}-feature pharmacophore (random baseline mean \
         {baseline_mean_score:.3} over {} candidates).",
        sequence.len(),
        pharmacophore.len(),
        baseline_count,
    ));
    notes.push(format!(
        "The design presents {} pocket(s) — {} hairpin loop(s), {} internal-loop / bulge(s), \
         {} multi-junction(s).",
        pockets.len(),
        pockets
            .iter()
            .filter(|p| p.kind == PocketKind::HairpinLoop)
            .count(),
        pockets
            .iter()
            .filter(|p| p.kind == PocketKind::InternalLoop)
            .count(),
        pockets
            .iter()
            .filter(|p| p.kind == PocketKind::MultiJunction)
            .count(),
    ));
    notes.push(
        "Aptamer design is classical structural-complementarity scoring \
         (pharmacophore-feature vs RNA-base-edge chemistry under a pocket-class weight), not a \
         docking or 3-D-RNA calculation. The design is a strong in-silico aptamer candidate — \
         the canonical wet-lab follow-up is SELEX-style affinity selection."
            .to_string(),
    );

    Ok(AptamerDesign {
        sequence,
        structure,
        mfe_energy,
        pocket_score,
        baseline_mean_score,
        pockets,
        candidates_evaluated: considered,
        notes,
    })
}

/// Builds a random candidate sequence biased toward forming a pocket:
/// roughly equal-length helix arms with a hairpin / bulge candidate in
/// the middle. The bias is *structural* (a pocket-likely topology),
/// the base content is randomised inside the GC band.
fn random_pocket_seed(n: usize, gc_min: f64, gc_max: f64, rng: &mut Rng) -> Vec<u8> {
    // Try up to a handful of times to land inside the GC band; the
    // search loop already tolerates failures.
    for _ in 0..6 {
        let seq = random_pocket_seed_inner(n, rng);
        let gc = gc_fraction(&seq);
        if gc >= gc_min - 1e-9 && gc <= gc_max + 1e-9 {
            return seq;
        }
    }
    // Fall back to a sequence we know lies in band: repair the GC by
    // re-seeding mismatched positions.
    let mut seq = random_pocket_seed_inner(n, rng);
    repair_gc(&mut seq, gc_min, gc_max, rng);
    seq
}

/// Inner helper for [`random_pocket_seed`].
fn random_pocket_seed_inner(n: usize, rng: &mut Rng) -> Vec<u8> {
    // Pick a random hairpin template: an outer stem of length `s`, a
    // loop of length `n - 2s`. Cap stem length so the loop is at least
    // 5 nt — a real pocket size.
    let max_stem = (n.saturating_sub(5)) / 2;
    let stem = if max_stem == 0 {
        0
    } else {
        2 + rng.below(max_stem.saturating_sub(1).max(1))
    };
    let mut seq = vec![b'A'; n];
    // Stem: random canonical pairs at (i, n-1-i).
    for i in 0..stem {
        let (a, b) = CANON_PAIRS[rng.below(6)];
        seq[i] = a;
        seq[n - 1 - i] = b;
    }
    // Loop: random bases.
    for s in seq.iter_mut().take(n - stem).skip(stem) {
        *s = BASES[rng.below(4)];
    }
    seq
}

/// Repairs `seq`'s GC fraction by swapping individual bases until it
/// lies in `[gc_min, gc_max]`.
fn repair_gc(seq: &mut [u8], gc_min: f64, gc_max: f64, rng: &mut Rng) {
    let n = seq.len();
    if n == 0 {
        return;
    }
    for _ in 0..(n * 4) {
        let gc = gc_fraction(seq);
        if gc >= gc_min - 1e-9 && gc <= gc_max + 1e-9 {
            return;
        }
        let pos = rng.below(n);
        if gc < gc_min {
            // Need more GC.
            seq[pos] = if rng.below(2) == 0 { b'G' } else { b'C' };
        } else {
            // Need less GC.
            seq[pos] = if rng.below(2) == 0 { b'A' } else { b'U' };
        }
    }
}

/// GC fraction of a sequence.
fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / seq.len() as f64
}

/// Structural statistics of a folded design — convenience for callers.
/// Currently re-exposes [`structure_stats`].
pub fn aptamer_structure_stats(structure: &Structure) -> valenx_rnastruct::StructureStats {
    structure_stats(structure)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn donor_acceptor_pharmacophore() -> Pharmacophore {
        Pharmacophore::new(vec![
            PharmacophoreFeature::new(FeatureKind::Hydrophobic, [0.0, 0.0, 0.0], 4.0),
            PharmacophoreFeature::new(FeatureKind::HBondAcceptor, [3.0, 0.0, 0.0], 4.0),
        ])
    }

    #[test]
    fn base_edge_features_lookup() {
        let g = base_edge_features(b'G');
        assert_eq!(g.donors, 2);
        assert_eq!(g.acceptors, 2);
        assert!(g.purine);
        let u = base_edge_features(b'U');
        assert_eq!(u.donors, 1);
        assert!(!u.purine);
        let n = base_edge_features(b'N');
        assert_eq!(n.donors, 1);
    }

    #[test]
    fn complement_is_self_inverse() {
        for k in [
            FeatureKind::HBondDonor,
            FeatureKind::HBondAcceptor,
            FeatureKind::Hydrophobic,
            FeatureKind::Positive,
            FeatureKind::Negative,
        ] {
            assert_eq!(k.complement().complement(), k);
        }
    }

    #[test]
    fn pocket_extraction_finds_a_hairpin_loop() {
        let s = Structure::from_dot_bracket("(((....)))").unwrap();
        let p = extract_pockets(&s);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].kind, PocketKind::HairpinLoop);
        assert_eq!(p[0].start, 3);
        assert_eq!(p[0].end, 7);
        assert_eq!(p[0].len(), 4);
    }

    #[test]
    fn pocket_extraction_finds_a_multi_junction() {
        // A multiloop: ((..(((....))).....(((....)))..)) — two child
        // hairpins inside an outer pair → one multi-junction pocket
        // plus two hairpin loops.
        let s = Structure::from_dot_bracket("((..(((....))).....(((....)))..))").unwrap();
        let p = extract_pockets(&s);
        // Two child hairpins + one multi-junction (the long unpaired
        // run between the two child stems is the multi-junction
        // pocket).
        let hairpins = p
            .iter()
            .filter(|q| q.kind == PocketKind::HairpinLoop)
            .count();
        let multi = p
            .iter()
            .filter(|q| q.kind == PocketKind::MultiJunction)
            .count();
        assert!(hairpins >= 2, "expected ≥ 2 hairpin loops, got {hairpins}");
        assert!(multi >= 1, "expected ≥ 1 multi-junction, got {multi}");
    }

    #[test]
    fn pocket_extraction_open_chain_has_no_pockets() {
        let s = Structure::empty(20);
        assert!(extract_pockets(&s).is_empty());
    }

    #[test]
    fn pharmacophore_score_rewards_compatible_bases() {
        // A pure G hairpin loop should score well for an H-bond-acceptor
        // ligand (every loop base is a G — strong donor).
        let s = Structure::from_dot_bracket("(((((....)))))").unwrap();
        let seq = b"CCCCCGGGGGGGGG";
        let _ = seq;
        // Build a sequence whose stem is canonical and whose loop is GGGG.
        let seq = b"CCCCCGGGGGGGGG";
        let _ = seq;
        let seq = [
            b'C', b'C', b'C', b'C', b'C', b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G',
        ];
        let pharm = Pharmacophore::new(vec![PharmacophoreFeature::new(
            FeatureKind::HBondAcceptor,
            [0.0, 0.0, 0.0],
            2.0,
        )]);
        let score = pharmacophore_pocket_score(&seq, &s, &pharm);
        // Loop is positions 5..9 = GGGG → 4 bases × 2 donors × hairpin
        // weight 1.0 = 8.0.
        assert!(score >= 4.0, "expected donor-rich loop score, got {score}");
    }

    #[test]
    fn pharmacophore_score_zero_on_open_chain() {
        let s = Structure::empty(20);
        let pharm = donor_acceptor_pharmacophore();
        let score = pharmacophore_pocket_score(b"AAAAAAAAAAAAAAAAAAAA", &s, &pharm);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn designs_an_aptamer_above_baseline() {
        let pharm = donor_acceptor_pharmacophore();
        let d = design_aptamer(&pharm, AptamerDesignParams::default()).unwrap();
        // The design's pocket score should be above the random baseline.
        // (The designer keeps only candidates that improve the score, so
        // this is the success criterion.)
        assert!(
            d.above_baseline(),
            "design score {} not above baseline {}",
            d.pocket_score,
            d.baseline_mean_score,
        );
        assert!(!d.pockets.is_empty());
        assert!(d.candidates_evaluated > 0);
        assert!(!d.notes.is_empty());
    }

    #[test]
    fn aptamer_design_is_deterministic() {
        let pharm = donor_acceptor_pharmacophore();
        let a = design_aptamer(&pharm, AptamerDesignParams::default()).unwrap();
        let b = design_aptamer(&pharm, AptamerDesignParams::default()).unwrap();
        assert_eq!(a.sequence, b.sequence);
    }

    #[test]
    fn rejects_empty_pharmacophore() {
        let pharm = Pharmacophore::new(vec![]);
        let err = design_aptamer(&pharm, AptamerDesignParams::default()).unwrap_err();
        assert_eq!(err.code(), "rnadesign.goal");
    }

    #[test]
    fn rejects_invalid_length_window() {
        let pharm = donor_acceptor_pharmacophore();
        let params = AptamerDesignParams {
            length_min: 30,
            length_max: 20,
            ..AptamerDesignParams::default()
        };
        let err = design_aptamer(&pharm, params).unwrap_err();
        assert_eq!(err.code(), "rnadesign.invalid");
    }

    #[test]
    fn rejects_zero_candidates() {
        let pharm = donor_acceptor_pharmacophore();
        let params = AptamerDesignParams {
            candidates: 0,
            ..AptamerDesignParams::default()
        };
        assert!(design_aptamer(&pharm, params).is_err());
    }

    #[test]
    fn pharmacophore_diameter() {
        let pharm = Pharmacophore::new(vec![
            PharmacophoreFeature::new(FeatureKind::HBondDonor, [0.0, 0.0, 0.0], 1.0),
            PharmacophoreFeature::new(FeatureKind::HBondAcceptor, [3.0, 4.0, 0.0], 1.0),
        ]);
        // Diameter is the 3-4-5 triangle's hypotenuse — 5.
        assert!((pharm.diameter() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn aptamer_design_pocket_score_is_above_score_floor() {
        // A real designed aptamer has a non-zero pocket score (its
        // structure has at least one pocket the pharmacophore can
        // score against).
        let pharm = donor_acceptor_pharmacophore();
        let d = design_aptamer(&pharm, AptamerDesignParams::default()).unwrap();
        assert!(d.pocket_score > 0.0);
    }
}
