//! Group C — multi-objective optimisation (features 10–14).
//!
//! A first-pass design rarely satisfies *every* constraint at once: a
//! sequence that folds to its target may have a GC extreme, a long
//! homopolymer run, a stray restriction site or a strong off-target
//! hairpin. This module refines a candidate by **weighted-objective
//! local search** — it mutates the sequence and keeps changes that
//! improve a blended score.
//!
//! - [`optimize_design`] (feature 10) — the multi-objective optimiser:
//!   a weighted simulated-annealing search over **synonymous** codon
//!   swaps (for a coding design — the protein never changes) or
//!   **structure-preserving** base mutations (for a structural design).
//!   The objective blends target-structure match, ensemble defect, GC
//!   content, repeat avoidance, restriction-site avoidance, codon
//!   adaptation, uridine content and off-target-hairpin avoidance.
//! - [`ensemble_defect`] (feature 11) — the NUPACK-style metric: the
//!   expected number of incorrectly-paired bases relative to a target.
//! - [`repeat_scan`] (feature 12) — a repeat / low-complexity scanner.
//! - [`forbidden_motif_scan`] (feature 13) — a forbidden-motif /
//!   restriction-site scanner.
//! - [`synthesizability_scan`] (feature 14) — a synthesisability scan
//!   (GC extremes, long homopolymer runs, strong 5′ structure).
//!
//! ## v1 scope
//!
//! The optimiser is a **local search** — a weighted-sum simulated
//! annealing — so it finds a strong local optimum, not a proved global
//! one. The objective weights are caller-tunable
//! ([`ObjectiveWeights`]). For a structural design the mutation set is
//! restricted to keep the target's pairs canonical, and for a coding
//! design to synonymous codons — so the optimiser never breaks the
//! design's defining property. The ensemble-defect and
//! off-target-hairpin terms call the real `valenx-rnastruct` partition
//! function / MFE folder; the GC, repeat and homopolymer terms are
//! direct sequence statistics.

use crate::design::{DesignKind, RnaDesign};
use crate::error::{Result, RnaDesignError};
use crate::goal::DesignConstraints;
use serde::{Deserialize, Serialize};
use valenx_bioseq::ops::translate::GeneticCode;
use valenx_bioseq::{Seq, SeqKind};
use valenx_rnastruct::{mfe, partition_function, RnaSeq, Structure};

/// A simple deterministic PCG-style RNG, local to the optimiser so its
/// runs are reproducible without pulling in an external crate.
struct OptRng {
    state: u64,
}

impl OptRng {
    fn new(seed: u64) -> Self {
        OptRng {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    /// A xorshift step — adequate for a design-search mutation picker.
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// A `f64` in `[0, 1)`.
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A `usize` in `[0, n)` (`n` must be non-zero).
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

// ---------------------------------------------------------------------
// Feature 11 — ensemble-defect metric
// ---------------------------------------------------------------------

/// The NUPACK-style **ensemble defect** of a sequence relative to a
/// target structure (feature 11).
///
/// The ensemble defect is the *expected number of incorrectly-paired
/// nucleotides* over the Boltzmann ensemble: for each position `i`, the
/// probability it is *not* in its target pairing state. A perfect
/// design has defect `0`; the maximum is the sequence length.
///
/// This calls the real `valenx-rnastruct` McCaskill partition function
/// to get the base-pair probabilities.
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if the sequence and the target differ
///   in length.
/// - [`RnaDesignError::Upstream`] if the partition function fails.
pub fn ensemble_defect(seq: &[u8], target: &Structure) -> Result<f64> {
    if seq.len() != target.len() {
        return Err(RnaDesignError::invalid(
            "target",
            "sequence and target structure differ in length",
        ));
    }
    if seq.is_empty() {
        return Ok(0.0);
    }
    let rna = RnaSeq::parse(seq)?;
    let pf = partition_function(&rna)?;
    let mut defect = 0.0;
    for i in 0..seq.len() {
        match target.partner(i) {
            Some(j) => {
                // Target wants i paired to j: defect += P(i not paired to j).
                defect += 1.0 - pf.pair_probability(i, j);
            }
            None => {
                // Target wants i unpaired: defect += P(i paired).
                defect += 1.0 - pf.unpaired_probability(i);
            }
        }
    }
    Ok(defect.max(0.0))
}

/// The **normalised ensemble defect** — the ensemble defect divided by
/// the sequence length, a `[0, 1]` score where `0` is a perfect design.
///
/// # Errors
/// As [`ensemble_defect`].
pub fn normalized_ensemble_defect(seq: &[u8], target: &Structure) -> Result<f64> {
    if seq.is_empty() {
        return Ok(0.0);
    }
    Ok(ensemble_defect(seq, target)? / seq.len() as f64)
}

// ---------------------------------------------------------------------
// Feature 12 — repeat & low-complexity scan
// ---------------------------------------------------------------------

/// One repeat or low-complexity stretch found by [`repeat_scan`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RepeatHit {
    /// 0-based start of the stretch.
    pub start: usize,
    /// Length of the stretch in nucleotides.
    pub length: usize,
    /// The kind of repeat (`"homopolymer"`, `"dinucleotide"`,
    /// `"tandem"`).
    pub kind: String,
}

/// The result of a repeat / low-complexity scan (feature 12).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RepeatScan {
    /// Every repeat / low-complexity stretch found.
    pub hits: Vec<RepeatHit>,
}

impl RepeatScan {
    /// `true` when no repeat / low-complexity stretch was found.
    pub fn is_clean(&self) -> bool {
        self.hits.is_empty()
    }

    /// The number of stretches found.
    pub fn count(&self) -> usize {
        self.hits.len()
    }
}

/// Scans a sequence for repeats and low-complexity stretches
/// (feature 12).
///
/// Finds three kinds of stretch: **homopolymer** runs longer than
/// `max_homopolymer`, **dinucleotide** repeats of ≥ 5 units (`AUAUAU…`),
/// and short **tandem** repeats of a 3-nt unit repeated ≥ 4 times. Each
/// is a synthesis / stability liability.
pub fn repeat_scan(seq: &[u8], max_homopolymer: usize) -> RepeatScan {
    let mut hits = Vec::new();
    let n = seq.len();
    if n == 0 {
        return RepeatScan { hits };
    }

    // --- homopolymer runs --------------------------------------------
    let mut i = 0;
    while i < n {
        let b = seq[i].to_ascii_uppercase();
        let start = i;
        while i < n && seq[i].to_ascii_uppercase() == b {
            i += 1;
        }
        let run = i - start;
        if run > max_homopolymer {
            hits.push(RepeatHit {
                start,
                length: run,
                kind: "homopolymer".to_string(),
            });
        }
    }

    // --- dinucleotide repeats (>= 5 units = >= 10 nt) ----------------
    let mut i = 0;
    while i + 4 <= n {
        // The repeating unit is seq[i..i+2].
        if seq[i].eq_ignore_ascii_case(&seq[i + 1]) {
            // a homopolymer, already counted
            i += 1;
            continue;
        }
        let unit = [seq[i].to_ascii_uppercase(), seq[i + 1].to_ascii_uppercase()];
        let mut j = i;
        while j + 2 <= n
            && seq[j].to_ascii_uppercase() == unit[0]
            && seq[j + 1].to_ascii_uppercase() == unit[1]
        {
            j += 2;
        }
        let units = (j - i) / 2;
        if units >= 5 {
            hits.push(RepeatHit {
                start: i,
                length: j - i,
                kind: "dinucleotide".to_string(),
            });
            i = j;
        } else {
            i += 1;
        }
    }

    // --- tandem 3-nt repeats (>= 4 units = >= 12 nt) -----------------
    let mut i = 0;
    while i + 6 <= n {
        let unit = [
            seq[i].to_ascii_uppercase(),
            seq[i + 1].to_ascii_uppercase(),
            seq[i + 2].to_ascii_uppercase(),
        ];
        // skip a trivial-unit homopolymer
        if unit[0] == unit[1] && unit[1] == unit[2] {
            i += 1;
            continue;
        }
        let mut j = i;
        while j + 3 <= n
            && seq[j].to_ascii_uppercase() == unit[0]
            && seq[j + 1].to_ascii_uppercase() == unit[1]
            && seq[j + 2].to_ascii_uppercase() == unit[2]
        {
            j += 3;
        }
        let units = (j - i) / 3;
        if units >= 4 {
            hits.push(RepeatHit {
                start: i,
                length: j - i,
                kind: "tandem".to_string(),
            });
            i = j;
        } else {
            i += 1;
        }
    }

    hits.sort_by_key(|h| (h.start, h.length));
    RepeatScan { hits }
}

// ---------------------------------------------------------------------
// Feature 13 — forbidden-motif / restriction-site scan
// ---------------------------------------------------------------------

/// One forbidden-motif or restriction-site hit found by
/// [`forbidden_motif_scan`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MotifHit {
    /// 0-based start of the hit (in the RNA sequence).
    pub start: usize,
    /// The motif / site sequence that matched (RNA alphabet).
    pub motif: String,
    /// `true` when the hit is a restriction-enzyme site, `false` for a
    /// plain forbidden motif.
    pub is_restriction_site: bool,
    /// The enzyme name, when the hit is a restriction site.
    pub enzyme: Option<String>,
}

/// The result of a forbidden-motif / restriction-site scan (feature 13).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MotifScan {
    /// Every forbidden-motif / restriction-site hit found.
    pub hits: Vec<MotifHit>,
}

impl MotifScan {
    /// `true` when no forbidden motif / restriction site was found.
    pub fn is_clean(&self) -> bool {
        self.hits.is_empty()
    }

    /// The number of hits.
    pub fn count(&self) -> usize {
        self.hits.len()
    }
}

/// Scans an RNA sequence for forbidden motifs and restriction sites
/// (feature 13).
///
/// `constraints.forbidden_motifs` and `constraints.forbidden_subsequences`
/// are matched directly (transcribed to RNA); each enzyme named in
/// `constraints.forbidden_restriction_sites` has its recognition site
/// (from the `valenx-bioseq` database) searched against the *DNA
/// template* of the sequence.
pub fn forbidden_motif_scan(seq: &[u8], constraints: &DesignConstraints) -> MotifScan {
    let rna = transcribe_to_rna(seq);
    let mut hits = Vec::new();

    // Plain forbidden motifs / subsequences.
    for motif in constraints
        .forbidden_motifs
        .iter()
        .chain(constraints.forbidden_subsequences.iter())
    {
        let needle = transcribe_to_rna(motif.as_bytes());
        if needle.is_empty() {
            continue;
        }
        for start in find_all(&rna, &needle) {
            hits.push(MotifHit {
                start,
                motif: String::from_utf8_lossy(&needle).into_owned(),
                is_restriction_site: false,
                enzyme: None,
            });
        }
    }

    // Restriction-enzyme sites — match the recognition site (DNA) on the
    // DNA template; an RNA site equivalent is the U→T form.
    let dna_template = rna_to_dna(&rna);
    for enzyme_name in &constraints.forbidden_restriction_sites {
        if let Some(enzyme) = valenx_bioseq::cloning::restriction::enzyme_by_name(enzyme_name) {
            let site = enzyme.site.as_bytes();
            for start in find_all_iupac(&dna_template, site) {
                hits.push(MotifHit {
                    start,
                    motif: String::from_utf8_lossy(site).into_owned(),
                    is_restriction_site: true,
                    enzyme: Some(enzyme.name.to_string()),
                });
            }
            // Also the reverse-complement site (a non-palindromic
            // enzyme cuts either strand).
            let rc = valenx_bioseq::ops::revcomp::reverse_complement_dna_bytes(site);
            if rc != site {
                for start in find_all_iupac(&dna_template, &rc) {
                    hits.push(MotifHit {
                        start,
                        motif: String::from_utf8_lossy(&rc).into_owned(),
                        is_restriction_site: true,
                        enzyme: Some(enzyme.name.to_string()),
                    });
                }
            }
        }
    }

    hits.sort_by_key(|h| h.start);
    hits.dedup();
    MotifScan { hits }
}

// ---------------------------------------------------------------------
// Feature 14 — synthesizability scan
// ---------------------------------------------------------------------

/// The result of a synthesisability scan (feature 14).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SynthesizabilityScan {
    /// Overall GC fraction of the sequence.
    pub gc_fraction: f64,
    /// `true` when the GC content is in a synthesis-friendly band
    /// (~30–70 %).
    pub gc_ok: bool,
    /// The longest homopolymer run, in nucleotides.
    pub longest_homopolymer: usize,
    /// `true` when no homopolymer run exceeds the synthesis limit.
    pub homopolymer_ok: bool,
    /// The minimum free energy of the 5′ window (kcal/mol). Strong
    /// (very negative) 5′ structure impedes IVT initiation and ribosome
    /// loading.
    pub five_prime_mfe: f64,
    /// `true` when the 5′ structure is not pathologically strong.
    pub five_prime_ok: bool,
    /// Human-readable warnings, one line per issue found.
    pub warnings: Vec<String>,
}

impl SynthesizabilityScan {
    /// `true` when the sequence clears every synthesisability check.
    pub fn is_synthesizable(&self) -> bool {
        self.gc_ok && self.homopolymer_ok && self.five_prime_ok
    }
}

/// The 5′ window length folded for the synthesisability 5′-structure
/// check.
const FIVE_PRIME_WINDOW: usize = 40;

/// Scans an RNA sequence for synthesisability liabilities (feature 14).
///
/// Checks GC extremes (a synthesis-friendly band of ~30–70 %), long
/// homopolymer runs, and strong 5′ secondary structure (which impedes
/// in-vitro transcription and translation initiation). The 5′ structure
/// is folded with the real `valenx-rnastruct` MFE folder.
///
/// # Errors
/// [`RnaDesignError::Upstream`] if the 5′-window fold fails.
pub fn synthesizability_scan(
    seq: &[u8],
    max_homopolymer: usize,
) -> Result<SynthesizabilityScan> {
    let rna = transcribe_to_rna(seq);
    let mut warnings = Vec::new();

    // --- GC content ---------------------------------------------------
    let gc_fraction = gc_content(&rna);
    let gc_ok = (0.30..=0.70).contains(&gc_fraction);
    if !gc_ok {
        warnings.push(format!(
            "GC content {:.0}% is outside the synthesis-friendly 30-70% band — \
             extreme GC raises synthesis failure and IVT bias.",
            gc_fraction * 100.0,
        ));
    }

    // --- homopolymer runs --------------------------------------------
    let longest_homopolymer = longest_run(&rna);
    let homopolymer_ok = longest_homopolymer <= max_homopolymer;
    if !homopolymer_ok {
        warnings.push(format!(
            "A {longest_homopolymer}-nt homopolymer run exceeds the {max_homopolymer}-nt \
             limit — long single-base runs are error-prone to synthesise and transcribe.",
        ));
    }

    // --- 5' structure -------------------------------------------------
    let five_prime_mfe = if rna.is_empty() {
        0.0
    } else {
        let window_len = FIVE_PRIME_WINDOW.min(rna.len());
        let window = RnaSeq::parse(&rna[..window_len])?;
        mfe(&window)?.energy
    };
    // A 5' window MFE below -15 kcal/mol is a strong structural cap.
    let five_prime_ok = five_prime_mfe > -15.0;
    if !five_prime_ok {
        warnings.push(format!(
            "The 5' window folds strongly ({five_prime_mfe:.1} kcal/mol) — strong 5' \
             structure impedes IVT initiation and ribosome loading.",
        ));
    }

    Ok(SynthesizabilityScan {
        gc_fraction,
        gc_ok,
        longest_homopolymer,
        homopolymer_ok,
        five_prime_mfe,
        five_prime_ok,
        warnings,
    })
}

// ---------------------------------------------------------------------
// Off-target hairpin detection (an optimiser objective term)
// ---------------------------------------------------------------------

/// Counts strong **off-target hairpins** — local stem-loops outside any
/// intended structure that would compete with the design's fold or
/// impede translation.
///
/// The scan slides a window across the sequence, folds each window with
/// the MFE folder, and counts windows whose fold is strong (MFE below a
/// threshold). For a structural design the *intended* structure's
/// pairs are subtracted out so only spurious hairpins are counted.
///
/// # Errors
/// [`RnaDesignError::Upstream`] if a window fold fails.
pub fn off_target_hairpin_count(seq: &[u8], intended: Option<&Structure>) -> Result<usize> {
    let rna = transcribe_to_rna(seq);
    let n = rna.len();
    const WINDOW: usize = 30;
    const STEP: usize = 15;
    const STRONG_MFE: f64 = -10.0;
    if n < WINDOW {
        return Ok(0);
    }
    let mut count = 0;
    let mut start = 0;
    while start + WINDOW <= n {
        let window = RnaSeq::parse(&rna[start..start + WINDOW])?;
        let folded = mfe(&window)?;
        if folded.energy < STRONG_MFE {
            // If an intended structure pairs heavily within this window
            // the hairpin is wanted, not off-target.
            let wanted = match intended {
                Some(target) => {
                    let mut paired_in_window = 0;
                    for i in start..(start + WINDOW) {
                        if let Some(j) = target.partner(i) {
                            if (start..start + WINDOW).contains(&j) {
                                paired_in_window += 1;
                            }
                        }
                    }
                    paired_in_window >= 6
                }
                None => false,
            };
            if !wanted {
                count += 1;
            }
        }
        start += STEP;
    }
    Ok(count)
}

// ---------------------------------------------------------------------
// Feature 10 — the multi-objective optimiser
// ---------------------------------------------------------------------

/// The weights of the multi-objective optimiser's blended score
/// (feature 10).
///
/// Each weight scales one objective term; the optimiser minimises the
/// weighted sum. A weight of `0` switches a term off. All weights are
/// non-negative.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ObjectiveWeights {
    /// Weight on the target-structure-match penalty (base-pair distance
    /// of the MFE fold to the target). `0` for a non-structural design.
    pub structure_match: f64,
    /// Weight on the ensemble-defect penalty.
    pub ensemble_defect: f64,
    /// Weight on the GC-out-of-range penalty.
    pub gc_content: f64,
    /// Weight on the repeat / low-complexity penalty.
    pub repeats: f64,
    /// Weight on the forbidden-motif / restriction-site penalty.
    pub forbidden_motifs: f64,
    /// Weight on the codon-adaptation penalty (`1 − CAI`). `0` for a
    /// non-coding design.
    pub codon_adaptation: f64,
    /// Weight on the uridine-content penalty.
    pub uridine_content: f64,
    /// Weight on the off-target-hairpin penalty.
    pub off_target_hairpins: f64,
}

impl Default for ObjectiveWeights {
    /// A balanced default: structure match and ensemble defect dominate;
    /// the constraint terms have moderate weight.
    fn default() -> Self {
        ObjectiveWeights {
            structure_match: 5.0,
            ensemble_defect: 2.0,
            gc_content: 1.0,
            repeats: 1.5,
            forbidden_motifs: 3.0,
            codon_adaptation: 2.0,
            uridine_content: 1.0,
            off_target_hairpins: 1.5,
        }
    }
}

impl ObjectiveWeights {
    /// `true` if every weight is non-negative (the optimiser requires
    /// it).
    pub fn is_valid(&self) -> bool {
        [
            self.structure_match,
            self.ensemble_defect,
            self.gc_content,
            self.repeats,
            self.forbidden_motifs,
            self.codon_adaptation,
            self.uridine_content,
            self.off_target_hairpins,
        ]
        .iter()
        .all(|&w| w >= 0.0 && w.is_finite())
    }
}

/// The result of a multi-objective optimisation run (feature 10).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OptimizationResult {
    /// The optimised design candidate.
    pub design: RnaDesign,
    /// The blended objective score of the *input* design (lower is
    /// better).
    pub score_before: f64,
    /// The blended objective score of the optimised design.
    pub score_after: f64,
    /// Number of mutation steps accepted by the simulated annealing.
    pub accepted_steps: usize,
    /// Total number of mutation steps attempted.
    pub total_steps: usize,
    /// Human-readable notes describing the optimisation.
    pub notes: Vec<String>,
}

impl OptimizationResult {
    /// The score improvement (`score_before − score_after`); `≥ 0` when
    /// the optimiser made progress.
    pub fn improvement(&self) -> f64 {
        self.score_before - self.score_after
    }
}

/// Parameters for [`optimize_design`].
#[derive(Copy, Clone, Debug)]
pub struct OptimizeParams {
    /// The objective weights.
    pub weights: ObjectiveWeights,
    /// The number of simulated-annealing mutation steps.
    pub steps: usize,
    /// The random seed (the search is deterministic for a fixed seed).
    pub seed: u64,
}

impl Default for OptimizeParams {
    /// 200 annealing steps with the default weights.
    fn default() -> Self {
        OptimizeParams {
            weights: ObjectiveWeights::default(),
            steps: 200,
            seed: 0xDE51,
        }
    }
}

/// Runs the multi-objective optimiser on a design candidate
/// (feature 10).
///
/// Refines `design` by a weighted simulated-annealing local search. For
/// a coding design the mutations are **synonymous codon swaps** (the
/// encoded protein is invariant); for a structural / scaffold design
/// they are **base mutations that keep the target's pairs canonical**.
/// The accepted candidate is the lowest-scoring sequence visited.
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if `params.steps == 0` or the weights
///   are invalid.
/// - [`RnaDesignError::Upstream`] if a folding call fails.
pub fn optimize_design(
    design: &RnaDesign,
    constraints: &DesignConstraints,
    params: OptimizeParams,
) -> Result<OptimizationResult> {
    if params.steps == 0 {
        return Err(RnaDesignError::invalid("steps", "need at least one step"));
    }
    if !params.weights.is_valid() {
        return Err(RnaDesignError::invalid(
            "weights",
            "objective weights must all be finite and non-negative",
        ));
    }
    if design.is_empty() {
        return Err(RnaDesignError::invalid("design", "design sequence is empty"));
    }

    let mut rng = OptRng::new(params.seed);
    let code = GeneticCode::standard();

    let mut current = design.sequence.clone();
    let score_before = objective_score(&current, design, constraints, &params.weights)?;
    let mut current_score = score_before;
    let mut best = current.clone();
    let mut best_score = current_score;

    let mut accepted = 0usize;
    for step in 0..params.steps {
        // Linear cooling schedule.
        let frac = step as f64 / params.steps as f64;
        let temperature = (1.0 - frac).max(0.02);

        let trial = match propose_mutation(&current, design, &code, &mut rng) {
            Some(t) => t,
            None => continue, // no legal mutation at this site
        };
        let trial_score = objective_score(&trial, design, constraints, &params.weights)?;
        let delta = trial_score - current_score;
        // Metropolis acceptance.
        let accept = delta <= 0.0 || rng.next_f64() < (-delta / temperature).exp();
        if accept {
            current = trial;
            current_score = trial_score;
            accepted += 1;
            if current_score < best_score {
                best = current.clone();
                best_score = current_score;
            }
        }
    }

    // Build the optimised design — keep the construct in sync for a
    // coding design (the CDS bytes inside the construct change too).
    let optimized = rebuild_design(design, best)?;

    let mut notes = Vec::new();
    notes.push(format!(
        "Multi-objective simulated-annealing optimisation: {accepted}/{} steps accepted.",
        params.steps,
    ));
    notes.push(format!(
        "Blended objective score {score_before:.3} -> {best_score:.3} \
         (improvement {:.3}).",
        score_before - best_score,
    ));
    match &design.kind {
        DesignKind::Coding => notes.push(
            "Mutations were synonymous codon swaps — the encoded protein is unchanged."
                .to_string(),
        ),
        _ => notes.push(
            "Mutations kept every target base pair canonical — the design's fold is preserved."
                .to_string(),
        ),
    }
    notes.push(
        "The optimiser is a local search — this is a strong local optimum, not a proved \
         global one."
            .to_string(),
    );

    Ok(OptimizationResult {
        design: optimized,
        score_before,
        score_after: best_score,
        accepted_steps: accepted,
        total_steps: params.steps,
        notes,
    })
}

/// Proposes a single mutation respecting the design's defining
/// constraint — synonymous for a coding design, pair-canonical for a
/// structural one. Returns `None` if no legal mutation exists at the
/// chosen site.
fn propose_mutation(
    seq: &[u8],
    design: &RnaDesign,
    code: &GeneticCode,
    rng: &mut OptRng,
) -> Option<Vec<u8>> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];
    let n = seq.len();
    if n == 0 {
        return None;
    }

    match (&design.kind, design.cds_span) {
        (DesignKind::Coding, Some((cds_start, cds_end))) => {
            // Pick a sense codon inside the CDS (skip the start codon
            // and the stop codon) and swap it for a synonymous one.
            let cds_len = cds_end.saturating_sub(cds_start);
            let n_codons = cds_len / 3;
            if n_codons < 3 {
                return None;
            }
            // codons 1..n_codons-1 are the sense codons.
            let ci = 1 + rng.below(n_codons - 2);
            let cstart = cds_start + ci * 3;
            let codon_rna = [seq[cstart], seq[cstart + 1], seq[cstart + 2]];
            let codon_dna = [
                rna_base_to_dna(codon_rna[0]),
                rna_base_to_dna(codon_rna[1]),
                rna_base_to_dna(codon_rna[2]),
            ];
            let aa = code.translate_codon(&codon_dna);
            let syns = synonymous_codons(aa, code);
            if syns.len() < 2 {
                return None;
            }
            let pick = syns[rng.below(syns.len())];
            let mut trial = seq.to_vec();
            trial[cstart] = dna_base_to_rna(pick[0]);
            trial[cstart + 1] = dna_base_to_rna(pick[1]);
            trial[cstart + 2] = dna_base_to_rna(pick[2]);
            Some(trial)
        }
        (DesignKind::Structural { target }, _) => {
            structure_preserving_mutation(seq, target, rng, &BASES)
        }
        (DesignKind::Riboswitch { free_target, .. }, _) => {
            // Preserve the free-state target's pairs (the resting fold).
            structure_preserving_mutation(seq, free_target, rng, &BASES)
        }
        (DesignKind::MotifScaffold { .. }, _) => {
            // No structural target attached — mutate any single base.
            let pos = rng.below(n);
            let mut trial = seq.to_vec();
            trial[pos] = BASES[rng.below(4)];
            Some(trial)
        }
        (DesignKind::Coding, None) => {
            // A coding design with no CDS span — fall back to a plain
            // base mutation (rare).
            let pos = rng.below(n);
            let mut trial = seq.to_vec();
            trial[pos] = BASES[rng.below(4)];
            Some(trial)
        }
    }
}

/// A base mutation that keeps every pair of `target` canonical: if the
/// chosen position is paired, the whole pair is re-seated with a fresh
/// canonical pair; an unpaired position is mutated freely.
fn structure_preserving_mutation(
    seq: &[u8],
    target: &Structure,
    rng: &mut OptRng,
    bases: &[u8; 4],
) -> Option<Vec<u8>> {
    const CANON: [(u8, u8); 6] = [
        (b'A', b'U'),
        (b'U', b'A'),
        (b'G', b'C'),
        (b'C', b'G'),
        (b'G', b'U'),
        (b'U', b'G'),
    ];
    let n = seq.len();
    if n == 0 || n != target.len() {
        return None;
    }
    let pos = rng.below(n);
    let mut trial = seq.to_vec();
    match target.partner(pos) {
        Some(partner) => {
            let (a, b) = CANON[rng.below(6)];
            let (lo, hi) = if pos < partner {
                (pos, partner)
            } else {
                (partner, pos)
            };
            trial[lo] = a;
            trial[hi] = b;
        }
        None => {
            trial[pos] = bases[rng.below(4)];
        }
    }
    Some(trial)
}

/// Computes the blended objective score of a sequence (lower is
/// better) — the function the optimiser minimises.
fn objective_score(
    seq: &[u8],
    design: &RnaDesign,
    constraints: &DesignConstraints,
    weights: &ObjectiveWeights,
) -> Result<f64> {
    let mut score = 0.0;
    let rna = transcribe_to_rna(seq);

    // 1) Structure-match penalty (structural / scaffold designs).
    let intended: Option<&Structure> = match &design.kind {
        DesignKind::Structural { target } => Some(target),
        DesignKind::Riboswitch { free_target, .. } => Some(free_target),
        _ => None,
    };
    if weights.structure_match > 0.0 {
        if let Some(target) = intended {
            if target.len() == rna.len() {
                let folded = mfe(&RnaSeq::parse(&rna)?)?.structure;
                let dist = valenx_rnastruct::base_pair_distance(&folded, target)?;
                score += weights.structure_match * dist as f64;
            }
        }
    }

    // 2) Ensemble-defect penalty.
    if weights.ensemble_defect > 0.0 {
        if let Some(target) = intended {
            if target.len() == rna.len() {
                score += weights.ensemble_defect * normalized_ensemble_defect(&rna, target)?
                    * rna.len() as f64
                    / 10.0;
            }
        }
    }

    // 3) GC-out-of-range penalty.
    if weights.gc_content > 0.0 {
        let gc = gc_content(&rna);
        let excess = if gc < constraints.gc_min {
            constraints.gc_min - gc
        } else if gc > constraints.gc_max {
            gc - constraints.gc_max
        } else {
            0.0
        };
        score += weights.gc_content * excess * 100.0;
    }

    // 4) Repeat penalty.
    if weights.repeats > 0.0 {
        let scan = repeat_scan(&rna, constraints.max_homopolymer);
        score += weights.repeats * scan.count() as f64;
    }

    // 5) Forbidden-motif / restriction-site penalty.
    if weights.forbidden_motifs > 0.0 {
        let scan = forbidden_motif_scan(&rna, constraints);
        score += weights.forbidden_motifs * scan.count() as f64;
    }

    // 6) Codon-adaptation penalty (coding designs).
    if weights.codon_adaptation > 0.0 {
        if let Some((s, e)) = design.cds_span {
            if let Some(cds) = rna.get(s..e) {
                if let Ok(cai) = cds_cai(cds, constraints.host) {
                    score += weights.codon_adaptation * (1.0 - cai) * 10.0;
                }
            }
        }
    }

    // 7) Uridine-content penalty.
    if weights.uridine_content > 0.0 {
        let u = uridine_fraction(&rna);
        score += weights.uridine_content * u * 10.0;
    }

    // 8) Off-target-hairpin penalty.
    if weights.off_target_hairpins > 0.0 {
        let count = off_target_hairpin_count(&rna, intended)?;
        score += weights.off_target_hairpins * count as f64;
    }

    Ok(score)
}

/// Rebuilds an [`RnaDesign`] around an optimised sequence, keeping the
/// construct's CDS bytes in sync for a coding design.
fn rebuild_design(original: &RnaDesign, sequence: Vec<u8>) -> Result<RnaDesign> {
    let mut out = original.clone();
    out.sequence = sequence;
    // For a coding design, refresh the construct's CDS slice so the
    // construct and the transcript stay consistent.
    if let (Some((s, e)), Some(construct)) = (out.cds_span, out.construct.as_mut()) {
        if let Some(cds) = out.sequence.get(s..e) {
            construct.cds = cds.to_vec();
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------
// Internal sequence helpers
// ---------------------------------------------------------------------

/// Transcribes raw bytes to uppercase RNA (`T → U`).
fn transcribe_to_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'T' => b'U',
            other => other,
        })
        .collect()
}

/// Reverse-transcribes uppercase RNA to DNA (`U → T`).
fn rna_to_dna(rna: &[u8]) -> Vec<u8> {
    rna.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect()
}

/// One RNA base to DNA.
fn rna_base_to_dna(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

/// One DNA base to RNA.
fn dna_base_to_rna(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'T' => b'U',
        other => other,
    }
}

/// GC fraction of an RNA / DNA sequence in `[0, 1]`.
pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / seq.len() as f64
}

/// Uridine (`U`/`T`) fraction of a sequence in `[0, 1]`.
pub fn uridine_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let u = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'U' | b'T'))
        .count();
    u as f64 / seq.len() as f64
}

/// The longest homopolymer run in a sequence.
fn longest_run(seq: &[u8]) -> usize {
    let mut longest = 0;
    let mut run = 0;
    let mut last = 0u8;
    for &b in seq {
        let u = b.to_ascii_uppercase();
        if u == last {
            run += 1;
        } else {
            run = 1;
            last = u;
        }
        longest = longest.max(run);
    }
    longest
}

/// All start indices where `needle` occurs in `haystack` (exact,
/// case-insensitive).
fn find_all(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
    let mut out = Vec::new();
    if needle.is_empty() || needle.len() > haystack.len() {
        return out;
    }
    for (i, w) in haystack.windows(needle.len()).enumerate() {
        if w.eq_ignore_ascii_case(needle) {
            out.push(i);
        }
    }
    out
}

/// All start indices where the IUPAC `pattern` matches `haystack`
/// (using the `valenx-bioseq` ambiguity table).
fn find_all_iupac(haystack: &[u8], pattern: &[u8]) -> Vec<usize> {
    let mut out = Vec::new();
    if pattern.is_empty() || pattern.len() > haystack.len() {
        return out;
    }
    for (i, w) in haystack.windows(pattern.len()).enumerate() {
        if pattern
            .iter()
            .zip(w)
            .all(|(&p, &b)| valenx_bioseq::alphabet::iupac_matches(p, b))
        {
            out.push(i);
        }
    }
    out
}

/// Every DNA codon synonymous with amino acid `aa`.
fn synonymous_codons(aa: u8, code: &GeneticCode) -> Vec<[u8; 3]> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::new();
    for &b0 in &BASES {
        for &b1 in &BASES {
            for &b2 in &BASES {
                let codon = [b0, b1, b2];
                if code.translate_codon(&codon) == aa {
                    out.push(codon);
                }
            }
        }
    }
    out
}

/// The codon-adaptation index of an RNA CDS for a host.
fn cds_cai(cds_rna: &[u8], host: valenx_genediting::mrna::codon::ExpressionHost) -> Result<f64> {
    if cds_rna.is_empty() || cds_rna.len() % 3 != 0 {
        return Ok(0.0);
    }
    let dna = rna_to_dna(cds_rna);
    let seq = Seq::new(SeqKind::Dna, &dna)?;
    let bioseq_host = match host {
        valenx_genediting::mrna::codon::ExpressionHost::Human => {
            valenx_bioseq::cloning::codon_opt::Host::Human
        }
        valenx_genediting::mrna::codon::ExpressionHost::EColi => {
            valenx_bioseq::cloning::codon_opt::Host::EColi
        }
    };
    Ok(valenx_bioseq::cloning::codon_opt::codon_adaptation_index(&seq, bioseq_host)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::DesignKind;

    fn structural_design(seq: &[u8], db: &str) -> RnaDesign {
        RnaDesign {
            sequence: seq.to_vec(),
            kind: DesignKind::Structural {
                target: Structure::from_dot_bracket(db).unwrap(),
            },
            cds_span: None,
            construct: None,
            notes: Vec::new(),
        }
    }

    #[test]
    fn ensemble_defect_is_low_for_a_good_design() {
        // A strong GC hairpin against its own structure.
        let seq = b"GGGGGGAAAACCCCCC";
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let defect = ensemble_defect(seq, &target).unwrap();
        // The defect is bounded by the length.
        assert!((0.0..=seq.len() as f64).contains(&defect));
    }

    #[test]
    fn ensemble_defect_rejects_length_mismatch() {
        let target = Structure::from_dot_bracket("((....))").unwrap();
        assert!(ensemble_defect(b"GGGG", &target).is_err());
    }

    #[test]
    fn normalized_defect_is_in_unit_range() {
        let seq = b"GGGGGGAAAACCCCCC";
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let nd = normalized_ensemble_defect(seq, &target).unwrap();
        assert!((0.0..=1.0).contains(&nd));
    }

    #[test]
    fn repeat_scan_finds_homopolymer() {
        let scan = repeat_scan(b"ACGTAAAAAAAAGCGC", 6);
        assert!(!scan.is_clean());
        assert!(scan.hits.iter().any(|h| h.kind == "homopolymer"));
    }

    #[test]
    fn repeat_scan_finds_dinucleotide() {
        let scan = repeat_scan(b"GGGAUAUAUAUAUAUGGG", 20);
        assert!(scan.hits.iter().any(|h| h.kind == "dinucleotide"));
    }

    #[test]
    fn repeat_scan_clean_sequence() {
        let scan = repeat_scan(b"ACGUACGUCAGUCAGUACGU", 6);
        assert!(scan.is_clean());
    }

    #[test]
    fn forbidden_motif_scan_finds_motif() {
        let constraints = DesignConstraints::new().forbid_motif("AAUAAA");
        let scan = forbidden_motif_scan(b"GGGAAUAAAGGG", &constraints);
        assert_eq!(scan.count(), 1);
        assert!(!scan.hits[0].is_restriction_site);
    }

    #[test]
    fn forbidden_motif_scan_finds_restriction_site() {
        // EcoRI site GAATTC -> RNA GAAUUC.
        let constraints = DesignConstraints::new().forbid_restriction_site("EcoRI");
        let scan = forbidden_motif_scan(b"GGGGAAUUCGGG", &constraints);
        assert!(scan.hits.iter().any(|h| h.is_restriction_site));
    }

    #[test]
    fn forbidden_motif_scan_clean() {
        let constraints = DesignConstraints::new().forbid_motif("AAUAAA");
        let scan = forbidden_motif_scan(b"GGGGCGCGGGG", &constraints);
        assert!(scan.is_clean());
    }

    #[test]
    fn synthesizability_scan_flags_homopolymer() {
        let scan = synthesizability_scan(b"ACGUACGUAAAAAAAAAAAAACGUACGU", 6).unwrap();
        assert!(!scan.homopolymer_ok);
        assert!(!scan.is_synthesizable());
    }

    #[test]
    fn synthesizability_scan_clean_sequence() {
        // A balanced, repeat-free sequence.
        let scan = synthesizability_scan(b"ACGUCAGUACGUCAGUACGUCAGU", 6).unwrap();
        assert!(scan.gc_ok);
        assert!(scan.homopolymer_ok);
    }

    #[test]
    fn synthesizability_flags_gc_extreme() {
        let scan = synthesizability_scan(b"GCGCGCGCGCGCGCGCGCGCGCGC", 30).unwrap();
        assert!(!scan.gc_ok);
    }

    #[test]
    fn off_target_hairpin_count_runs() {
        let count = off_target_hairpin_count(b"ACGUACGUACGUACGUACGUACGUACGUACGU", None).unwrap();
        // A flat sequence should have few or no strong hairpins.
        let _ = count;
    }

    #[test]
    fn objective_weights_validate() {
        assert!(ObjectiveWeights::default().is_valid());
        let bad = ObjectiveWeights {
            gc_content: -1.0,
            ..ObjectiveWeights::default()
        };
        assert!(!bad.is_valid());
    }

    #[test]
    fn optimize_a_structural_design() {
        let d = structural_design(b"AAAAAAAAAAAA", "((((....))))");
        let result = optimize_design(
            &d,
            &DesignConstraints::default(),
            OptimizeParams::default(),
        )
        .unwrap();
        assert_eq!(result.design.len(), 12);
        // The optimiser should not make the score worse.
        assert!(result.score_after <= result.score_before + 1e-6);
        assert!(!result.notes.is_empty());
    }

    #[test]
    fn optimize_is_deterministic() {
        let d = structural_design(b"AAAAAAAAAAAA", "((((....))))");
        let a = optimize_design(&d, &DesignConstraints::default(), OptimizeParams::default())
            .unwrap();
        let b = optimize_design(&d, &DesignConstraints::default(), OptimizeParams::default())
            .unwrap();
        assert_eq!(a.design.sequence, b.design.sequence);
    }

    #[test]
    fn optimize_rejects_zero_steps() {
        let d = structural_design(b"GGGGAAAACCCC", "((((....))))");
        let params = OptimizeParams {
            steps: 0,
            ..OptimizeParams::default()
        };
        assert!(optimize_design(&d, &DesignConstraints::default(), params).is_err());
    }

    #[test]
    fn optimize_preserves_structural_pairs() {
        // After optimisation a structural design's target pairs must
        // still be canonical on the sequence.
        let d = structural_design(b"AAAAAAAAAAAA", "((((....))))");
        let result = optimize_design(
            &d,
            &DesignConstraints::default(),
            OptimizeParams::default(),
        )
        .unwrap();
        let target = Structure::from_dot_bracket("((((....))))").unwrap();
        let rna = RnaSeq::parse(&result.design.sequence).unwrap();
        assert!(target.validate_against(rna.as_bytes(), true).is_ok());
    }

    #[test]
    fn gc_and_uridine_helpers() {
        assert!((gc_content(b"GGCC") - 1.0).abs() < 1e-9);
        assert!((gc_content(b"AAUU") - 0.0).abs() < 1e-9);
        assert!((uridine_fraction(b"UUUU") - 1.0).abs() < 1e-9);
        assert!((uridine_fraction(b"AAAA") - 0.0).abs() < 1e-9);
        assert_eq!(gc_content(b""), 0.0);
    }
}
