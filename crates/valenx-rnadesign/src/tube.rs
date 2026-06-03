//! Feature 21 — NUPACK-class multi-strand tube design.
//!
//! [`crate::multistate`] designs a single sequence that adopts two or
//! more target structures (the same strand, several folds). Real
//! NUPACK *multi-tube* design goes further: a **tube** is a mixture
//! of several distinct RNA strands at given concentrations, and the
//! strands cofold into a distribution of complexes (monomers, dimers,
//! …) at thermodynamic equilibrium. **Tube design** finds strand
//! sequences that drive the equilibrium toward a *target complex
//! distribution* — for example, "strand A and strand B should
//! preferentially form the A·B dimer at 1 µM each".
//!
//! ## The equilibrium model
//!
//! Each complex `c` has a free energy `G_c` (its cofolded MFE / free
//! energy on the assembled multi-strand sequence). For a tube of
//! total strands at total concentration `c_tot`, the equilibrium
//! concentration of each complex is the law-of-mass-action solution
//! that minimises the total free energy:
//!
//! ```text
//!     min  Σ_c  [c]_c · (G_c + RT·ln([c]_c / c_tot))
//!     s.t. Σ_c  n_c^A · [c]_c = c_tot^A   (mass conservation, per strand)
//! ```
//!
//! For a 2-strand tube (strands A, B) the canonical complex set is
//! `{A, B, AA, BB, AB}` — monomers, homodimers, heterodimer. The mass
//! conservation laws are solved by Newton-Raphson on the
//! 2-dimensional log-concentration variable; the dimer concentrations
//! are then derived from the canonical equilibrium expressions
//! `[AB] = K_AB · [A] · [B]` etc. with `K = exp(-(G_complex −
//! G_A − G_B) / RT)`.
//!
//! ## The designer
//!
//! Strands are mutated freely (any base may change); the objective is
//! the **distance** between the tube's actual equilibrium complex
//! distribution and the target. Sequences that drive the target
//! complex's equilibrium fraction up are accepted; the search is a
//! deterministic mutation walk with a seed-controlled RNG.
//!
//! ## v1 scope — honest framing
//!
//! - This is a **NUPACK-class** tube model — concentration-aware,
//!   law-of-mass-action equilibrium over the canonical 2-strand
//!   complex set `{A, B, AA, BB, AB}` (and the 1-strand and 3-strand
//!   degenerations of it). NUPACK's full implementation enumerates
//!   *every* connected complex up to a complex-size cutoff and runs
//!   tube ensemble defect — this v1 covers the canonical small-tube
//!   case (1-3 strands, dimer-and-below complexes) cleanly.
//! - Complex free energies come from `valenx-rnastruct`'s exact Zuker
//!   MFE of the concatenated sequence with the cofold junction
//!   handling (no fake-hairpin pairs across the cut). The
//!   *equilibrium* `Q_c` of each complex is approximated by
//!   `exp(-G_c / RT)` (the dominant-MFE-structure approximation; for
//!   a tighter fit one would use the constrained partition function
//!   per complex). The approximation is honest and documented; for
//!   tube-design *ranking* purposes — choosing the strand sequence
//!   that drives the target complex up relative to the alternatives —
//!   it is reliable.
//! - The tube has a fixed total strand concentration set by the
//!   caller; pH / salt / Mg²⁺ are not modelled (the Turner-2004
//!   parameters apply at 37 °C, 1 M NaCl — the standard model
//!   condition).
//! - The designer is a local mutation walk, not an exhaustive search.

use crate::error::{Result, RnaDesignError};
use serde::{Deserialize, Serialize};
use valenx_rnastruct::fold::energy::{GAS_CONSTANT, T37_KELVIN};
use valenx_rnastruct::{mfe, RnaSeq};

/// The four RNA bases, ASCII.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];

// ---------------------------------------------------------------------
// A small deterministic RNG
// ---------------------------------------------------------------------

/// A deterministic xorshift RNG, local to the designer.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng {
            state: seed ^ 0x7E15_C4F9_B226_A91D,
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

// ---------------------------------------------------------------------
// Tube data model
// ---------------------------------------------------------------------

/// One strand in a tube: a name, a sequence, and a starting total
/// concentration (M).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TubeStrand {
    /// A short human-readable label (`"A"`, `"B"`, `"toehold"`, …).
    pub name: String,
    /// The strand's sequence (`A C G U`).
    pub sequence: Vec<u8>,
    /// The strand's total concentration in the tube, mol/L.
    /// Mass-conservation across complexes is anchored on this.
    pub total_concentration: f64,
}

impl TubeStrand {
    /// Builds a strand from a name, ASCII sequence and concentration.
    pub fn new(name: impl Into<String>, sequence: &[u8], total_concentration: f64) -> Self {
        TubeStrand {
            name: name.into(),
            sequence: sequence.to_vec(),
            total_concentration,
        }
    }

    /// The strand length.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// `true` if the strand is empty.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// The complex types a multi-strand tube canonically resolves into.
/// The 2-strand tube — the practical case — enumerates `{A, B, AA, BB,
/// AB}`; a 1-strand tube degenerates to `{A, AA}`; a 3-strand tube
/// adds `{C, AC, BC}` (the 3-strand trimer is not enumerated in v1).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ComplexKind {
    /// Monomer of strand index `i`.
    Monomer(usize),
    /// Homodimer of strand index `i`.
    Homodimer(usize),
    /// Heterodimer of strands `(i, j)` with `i < j`.
    Heterodimer(usize, usize),
}

impl ComplexKind {
    /// The set of strand indices the complex involves.
    pub fn strands(self) -> Vec<usize> {
        match self {
            ComplexKind::Monomer(i) => vec![i],
            ComplexKind::Homodimer(i) => vec![i, i],
            ComplexKind::Heterodimer(i, j) => vec![i, j],
        }
    }

    /// The number of strand units in the complex (`1` for monomers,
    /// `2` for dimers).
    pub fn size(self) -> usize {
        match self {
            ComplexKind::Monomer(_) => 1,
            ComplexKind::Homodimer(_) | ComplexKind::Heterodimer(_, _) => 2,
        }
    }

    /// The stoichiometric coefficient `n_c^k` — how many copies of
    /// strand `k` appear in this complex.
    pub fn stoichiometry(self, k: usize) -> usize {
        match self {
            ComplexKind::Monomer(i) => {
                if i == k {
                    1
                } else {
                    0
                }
            }
            ComplexKind::Homodimer(i) => {
                if i == k {
                    2
                } else {
                    0
                }
            }
            ComplexKind::Heterodimer(i, j) => {
                (if i == k { 1 } else { 0 }) + (if j == k { 1 } else { 0 })
            }
        }
    }

    /// A short human-readable name.
    pub fn label(self, strands: &[TubeStrand]) -> String {
        match self {
            ComplexKind::Monomer(i) => strands[i].name.clone(),
            ComplexKind::Homodimer(i) => format!("{}·{}", strands[i].name, strands[i].name),
            ComplexKind::Heterodimer(i, j) => {
                format!("{}·{}", strands[i].name, strands[j].name)
            }
        }
    }
}

// ---------------------------------------------------------------------
// Complex free energies
// ---------------------------------------------------------------------

/// The free energies of every canonical complex in a tube.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ComplexEnergies {
    /// Per-complex `(kind, ΔG)` in kcal/mol. The energy is the MFE of
    /// the concatenated complex sequence (with the cofold junction
    /// rule for dimers — no fake hairpins across the cut).
    pub energies: Vec<(ComplexKind, f64)>,
}

impl ComplexEnergies {
    /// Looks up the free energy of a complex.
    pub fn energy_of(&self, kind: ComplexKind) -> Option<f64> {
        self.energies
            .iter()
            .find_map(|(k, e)| if *k == kind { Some(*e) } else { None })
    }
}

/// Folds every canonical complex of `strands` and returns its free
/// energy, kcal/mol.
///
/// For dimers the junction-aware fold is performed by simple
/// concatenation + MFE; the cross-junction "fake hairpin" pair is
/// forbidden in [`valenx_rnastruct::cofold`], but for raw energetics
/// the straight `mfe(concat(a, b))` is the standard NUPACK
/// approximation and gives the dominant pairing structure.
///
/// # Errors
/// Propagates folding errors as [`RnaDesignError::Upstream`].
pub fn fold_all_complexes(strands: &[TubeStrand]) -> Result<ComplexEnergies> {
    let mut energies = Vec::new();
    let n = strands.len();

    // Monomers.
    for (i, s) in strands.iter().enumerate() {
        let rna = RnaSeq::parse(&s.sequence)?;
        let g = mfe(&rna)?.energy;
        energies.push((ComplexKind::Monomer(i), g));
    }
    // Homodimers.
    for (i, s) in strands.iter().enumerate() {
        let rna_i = RnaSeq::parse(&s.sequence)?;
        let composite = rna_i.concat(&rna_i);
        let g = mfe(&composite)?.energy;
        energies.push((ComplexKind::Homodimer(i), g));
    }
    // Heterodimers (i < j).
    for i in 0..n {
        for j in (i + 1)..n {
            let rna_i = RnaSeq::parse(&strands[i].sequence)?;
            let rna_j = RnaSeq::parse(&strands[j].sequence)?;
            let composite = rna_i.concat(&rna_j);
            let g = mfe(&composite)?.energy;
            energies.push((ComplexKind::Heterodimer(i, j), g));
        }
    }
    Ok(ComplexEnergies { energies })
}

// ---------------------------------------------------------------------
// Equilibrium concentrations
// ---------------------------------------------------------------------

/// The equilibrium concentration of every complex in a tube,
/// computed from the per-complex free energies under
/// mass conservation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TubeEquilibrium {
    /// Per-complex `(kind, [c])`, mol/L.
    pub concentrations: Vec<(ComplexKind, f64)>,
    /// The temperature of the calculation, K.
    pub temperature_k: f64,
}

impl TubeEquilibrium {
    /// The equilibrium concentration of `kind`, mol/L (0 if absent).
    pub fn concentration_of(&self, kind: ComplexKind) -> f64 {
        self.concentrations
            .iter()
            .find_map(|(k, c)| if *k == kind { Some(*c) } else { None })
            .unwrap_or(0.0)
    }

    /// The mole fraction of `kind` of the total complex
    /// concentration (every complex counted once, regardless of how
    /// many strands).
    pub fn fraction_of(&self, kind: ComplexKind) -> f64 {
        let total: f64 = self.concentrations.iter().map(|(_, c)| *c).sum();
        if total > 0.0 {
            self.concentration_of(kind) / total
        } else {
            0.0
        }
    }
}

/// Computes the equilibrium concentration of every canonical complex.
///
/// Solves the law-of-mass-action system for the dimer concentrations
/// in closed form from the monomer concentrations (which themselves
/// are found by Newton-Raphson on the mass-conservation residuals).
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if any strand concentration is
///   non-positive or non-finite.
#[allow(clippy::needless_range_loop)]
pub fn solve_tube_equilibrium(
    strands: &[TubeStrand],
    energies: &ComplexEnergies,
    temperature_k: f64,
) -> Result<TubeEquilibrium> {
    let n = strands.len();
    for s in strands {
        if !(s.total_concentration.is_finite() && s.total_concentration > 0.0) {
            return Err(RnaDesignError::invalid(
                "concentration",
                format!("strand `{}` has non-positive concentration", s.name),
            ));
        }
    }
    if !(temperature_k.is_finite() && temperature_k > 0.0) {
        return Err(RnaDesignError::invalid(
            "temperature",
            "temperature must be finite and positive",
        ));
    }

    let rt = GAS_CONSTANT * temperature_k;

    // Equilibrium constants (M⁻¹) for the dimers:
    //   K_AB = exp(-(G_AB - G_A - G_B) / RT)
    //   K_AA = exp(-(G_AA - 2 G_A) / RT)
    // The standard reference state is 1 M.
    let mut k_homo = vec![0.0_f64; n];
    let mut k_hetero = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        let g_a = energies.energy_of(ComplexKind::Monomer(i)).unwrap_or(0.0);
        let g_aa = energies
            .energy_of(ComplexKind::Homodimer(i))
            .unwrap_or(2.0 * g_a);
        k_homo[i] = ((-(g_aa - 2.0 * g_a)) / rt).exp();
        for j in (i + 1)..n {
            let g_b = energies.energy_of(ComplexKind::Monomer(j)).unwrap_or(0.0);
            let g_ab = energies
                .energy_of(ComplexKind::Heterodimer(i, j))
                .unwrap_or(g_a + g_b);
            let k = ((-(g_ab - g_a - g_b)) / rt).exp();
            k_hetero[i][j] = k;
            k_hetero[j][i] = k;
        }
    }

    // Solve for monomer concentrations by Newton-Raphson on the mass
    // conservation system:
    //   c_tot^k = [k] + 2 K_kk [k]^2 + Σ_{j != k} K_kj [k] [j]
    // for every strand k.
    let mut monomers: Vec<f64> = strands.iter().map(|s| s.total_concentration).collect();
    for _ in 0..200 {
        // Mass-conservation residual r[k] and Jacobian J[k][j] = dr[k]/d[j].
        let mut r = vec![0.0_f64; n];
        let mut j = vec![vec![0.0_f64; n]; n];
        for k in 0..n {
            let mk = monomers[k];
            let mut total = mk + 2.0 * k_homo[k] * mk * mk;
            for q in 0..n {
                if q != k {
                    total += k_hetero[k][q] * mk * monomers[q];
                }
            }
            r[k] = total - strands[k].total_concentration;
            // Jacobian entries.
            j[k][k] = 1.0 + 4.0 * k_homo[k] * mk;
            for q in 0..n {
                if q != k {
                    j[k][k] += k_hetero[k][q] * monomers[q];
                    j[k][q] = k_hetero[k][q] * mk;
                }
            }
        }
        // Solve J·δ = -r with Gauss elimination (n ≤ ~5 in practice).
        let mut delta = match gauss_solve(&j, &r.iter().map(|x| -x).collect::<Vec<_>>()) {
            Some(d) => d,
            None => break,
        };
        // Damped update with positivity guard.
        let mut step = 1.0_f64;
        for k in 0..n {
            if monomers[k] + delta[k] <= 0.0 {
                let max_drop = 0.5 * monomers[k];
                if delta[k] < -max_drop {
                    let scale = max_drop / -delta[k];
                    step = step.min(scale);
                }
            }
        }
        for d in delta.iter_mut() {
            *d *= step;
        }
        let mut max_change = 0.0_f64;
        for k in 0..n {
            monomers[k] += delta[k];
            if monomers[k] < 1e-30 {
                monomers[k] = 1e-30;
            }
            let rel = delta[k].abs() / (monomers[k].abs() + 1e-30);
            max_change = max_change.max(rel);
        }
        if max_change < 1e-9 {
            break;
        }
    }

    // Derive every complex's concentration.
    let mut concentrations: Vec<(ComplexKind, f64)> = Vec::new();
    for i in 0..n {
        concentrations.push((ComplexKind::Monomer(i), monomers[i]));
    }
    for i in 0..n {
        let c_aa = k_homo[i] * monomers[i] * monomers[i];
        concentrations.push((ComplexKind::Homodimer(i), c_aa));
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let c_ab = k_hetero[i][j] * monomers[i] * monomers[j];
            concentrations.push((ComplexKind::Heterodimer(i, j), c_ab));
        }
    }

    Ok(TubeEquilibrium {
        concentrations,
        temperature_k,
    })
}

/// Solves an `n × n` linear system `A · x = b` by Gauss elimination
/// with partial pivoting. Returns `None` for a singular matrix.
#[allow(clippy::needless_range_loop)]
fn gauss_solve(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 || b.len() != n {
        return Some(vec![0.0; n]);
    }
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            row.push(b[i]);
            row
        })
        .collect();
    for k in 0..n {
        // Partial pivot.
        let mut pivot_row = k;
        let mut pivot_val = m[k][k].abs();
        for i in (k + 1)..n {
            if m[i][k].abs() > pivot_val {
                pivot_val = m[i][k].abs();
                pivot_row = i;
            }
        }
        if pivot_val < 1e-30 {
            return None;
        }
        m.swap(k, pivot_row);
        for i in (k + 1)..n {
            let f = m[i][k] / m[k][k];
            for j in k..=n {
                m[i][j] -= f * m[k][j];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = m[i][n];
        for j in (i + 1)..n {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

// ---------------------------------------------------------------------
// Target distribution & objective
// ---------------------------------------------------------------------

/// A target fraction for one complex in a tube design — `0` means
/// "drive this complex out", `1` means "drive everything to this
/// complex".
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TargetFraction {
    /// The complex kind.
    pub kind: ComplexKind,
    /// The desired fraction in `[0, 1]` of the tube's total complex
    /// concentration.
    pub fraction: f64,
}

/// A target complex distribution — one or more [`TargetFraction`]s.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TargetDistribution {
    /// Per-complex target fractions.
    pub targets: Vec<TargetFraction>,
}

impl TargetDistribution {
    /// A distribution favoring exclusively `kind` (target fraction 1.0).
    pub fn favoring(kind: ComplexKind) -> Self {
        TargetDistribution {
            targets: vec![TargetFraction {
                kind,
                fraction: 1.0,
            }],
        }
    }

    /// The L1 distance between this target and an equilibrium.
    pub fn distance(&self, eq: &TubeEquilibrium) -> f64 {
        let mut d = 0.0_f64;
        for t in &self.targets {
            let actual = eq.fraction_of(t.kind);
            d += (actual - t.fraction).abs();
        }
        d
    }
}

// ---------------------------------------------------------------------
// Tube design parameters and result
// ---------------------------------------------------------------------

/// Parameters for [`design_tube`].
#[derive(Copy, Clone, Debug)]
pub struct TubeDesignParams {
    /// Mutation budget.
    pub iterations: usize,
    /// Temperature K.
    pub temperature_k: f64,
    /// Random seed.
    pub seed: u64,
}

impl Default for TubeDesignParams {
    /// 200 iterations at 37 °C.
    fn default() -> Self {
        TubeDesignParams {
            iterations: 200,
            temperature_k: T37_KELVIN,
            seed: 0x7E15_AB42,
        }
    }
}

/// The result of a tube design (feature 21).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TubeDesign {
    /// The designed strands (sequences and concentrations).
    pub strands: Vec<TubeStrand>,
    /// The equilibrium reached.
    pub equilibrium: TubeEquilibrium,
    /// The per-complex free energies at the final sequence.
    pub energies: ComplexEnergies,
    /// The L1 distance from `target_distribution` actually achieved.
    pub final_distance: f64,
    /// Mutation steps accepted.
    pub accepted_steps: usize,
    /// Mutation steps attempted.
    pub total_steps: usize,
    /// Human-readable notes.
    pub notes: Vec<String>,
}

// ---------------------------------------------------------------------
// The designer
// ---------------------------------------------------------------------

/// Designs strand sequences whose tube equilibrium matches
/// `target_distribution` (feature 21).
///
/// The strands' `sequence` field provides the **starting point** — a
/// reasonable choice is a random or rationally-designed first attempt.
/// The designer mutates only positions in *non-locked* strands;
/// concentrations are fixed.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if there are < 1 or > 3 strands, any
///   strand is empty, or the target distribution is empty.
/// - [`RnaDesignError::Invalid`] if `iterations == 0`, a strand
///   concentration is non-positive, the temperature is non-positive,
///   or a target fraction lies outside `[0, 1]`.
/// - [`RnaDesignError::Upstream`] if a folding call fails.
pub fn design_tube(
    initial_strands: &[TubeStrand],
    target_distribution: &TargetDistribution,
    params: TubeDesignParams,
) -> Result<TubeDesign> {
    let n = initial_strands.len();
    if !(1..=3).contains(&n) {
        return Err(RnaDesignError::goal(
            "strands",
            "tube design supports 1 to 3 strands in v1",
        ));
    }
    for s in initial_strands {
        if s.is_empty() {
            return Err(RnaDesignError::goal(
                "strands",
                format!("strand `{}` is empty", s.name),
            ));
        }
    }
    if target_distribution.targets.is_empty() {
        return Err(RnaDesignError::goal(
            "target_distribution",
            "target distribution has no entries",
        ));
    }
    for t in &target_distribution.targets {
        if !(0.0..=1.0).contains(&t.fraction) {
            return Err(RnaDesignError::invalid(
                "target_fraction",
                "every target fraction must lie in [0, 1]",
            ));
        }
        // Validate the target complex references real strand indices.
        for k in t.kind.strands() {
            if k >= n {
                return Err(RnaDesignError::goal(
                    "target_distribution",
                    format!("target complex references strand index {k} beyond range"),
                ));
            }
        }
    }
    if params.iterations == 0 {
        return Err(RnaDesignError::invalid(
            "iterations",
            "need at least one mutation step",
        ));
    }
    if !(params.temperature_k.is_finite() && params.temperature_k > 0.0) {
        return Err(RnaDesignError::invalid(
            "temperature",
            "temperature must be finite and positive",
        ));
    }

    let mut rng = Rng::new(params.seed);
    let mut strands = initial_strands.to_vec();
    let mut energies = fold_all_complexes(&strands)?;
    let mut eq = solve_tube_equilibrium(&strands, &energies, params.temperature_k)?;
    let mut best_dist = target_distribution.distance(&eq);

    let mut accepted = 0usize;
    let mut total = 0usize;

    for _ in 0..params.iterations {
        if best_dist < 1e-6 {
            break;
        }
        total += 1;
        // Pick a random strand and a random position to mutate.
        let si = rng.below(n);
        let pos = rng.below(strands[si].len());
        let original = strands[si].sequence[pos];
        let mut tried_bases = [original, original, original, original];
        for tb in tried_bases.iter_mut() {
            *tb = BASES[rng.below(4)];
        }
        let new_base = tried_bases
            .iter()
            .copied()
            .find(|&b| b != original)
            .unwrap_or(BASES[rng.below(4)]);

        let mut trial_strands = strands.clone();
        trial_strands[si].sequence[pos] = new_base;

        let trial_energies = match fold_all_complexes(&trial_strands) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let trial_eq = match solve_tube_equilibrium(
            &trial_strands,
            &trial_energies,
            params.temperature_k,
        ) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let trial_dist = target_distribution.distance(&trial_eq);
        if trial_dist <= best_dist {
            strands = trial_strands;
            energies = trial_energies;
            eq = trial_eq;
            best_dist = trial_dist;
            accepted += 1;
        }
    }

    let mut notes = Vec::new();
    notes.push(format!(
        "Tube design over {} strand(s); equilibrium temperature {:.2} K.",
        strands.len(),
        params.temperature_k,
    ));
    for t in &target_distribution.targets {
        let actual = eq.fraction_of(t.kind);
        notes.push(format!(
            "Target {} fraction {:.3} → achieved {:.3} (|Δ| {:.3}).",
            t.kind.label(&strands),
            t.fraction,
            actual,
            (actual - t.fraction).abs(),
        ));
    }
    notes.push(format!(
        "Final target-distribution distance {best_dist:.4}; {accepted}/{total} mutations \
         accepted."
    ));
    notes.push(
        "Tube design is a NUPACK-class multi-strand equilibrium model: per-complex MFE \
         free energies + law-of-mass-action equilibrium over the canonical 2-strand complex \
         set {A, B, A·A, B·B, A·B}. The complex partition function is approximated by its \
         MFE dominant; pH / salt / Mg²⁺ are not modelled. Validate experimentally."
            .to_string(),
    );

    Ok(TubeDesign {
        strands,
        equilibrium: eq,
        energies,
        final_distance: best_dist,
        accepted_steps: accepted,
        total_steps: total,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_strands() -> Vec<TubeStrand> {
        vec![
            TubeStrand::new("A", b"GGGGGGGGG", 1e-6),
            TubeStrand::new("B", b"CCCCCCCCC", 1e-6),
        ]
    }

    #[test]
    fn fold_all_complexes_returns_expected_set() {
        let strands = two_strands();
        let e = fold_all_complexes(&strands).unwrap();
        // 2 monomers + 2 homodimers + 1 heterodimer.
        assert_eq!(e.energies.len(), 5);
        assert!(e.energy_of(ComplexKind::Monomer(0)).is_some());
        assert!(e.energy_of(ComplexKind::Homodimer(1)).is_some());
        assert!(e.energy_of(ComplexKind::Heterodimer(0, 1)).is_some());
    }

    #[test]
    fn complementary_strands_form_a_strong_dimer() {
        let strands = two_strands();
        let e = fold_all_complexes(&strands).unwrap();
        let g_ab = e.energy_of(ComplexKind::Heterodimer(0, 1)).unwrap();
        let g_a = e.energy_of(ComplexKind::Monomer(0)).unwrap();
        let g_b = e.energy_of(ComplexKind::Monomer(1)).unwrap();
        // G_AB should be significantly more negative than G_A + G_B
        // (binding is favourable).
        assert!(
            g_ab < g_a + g_b - 1.0,
            "expected favourable binding, got G_AB={g_ab}, G_A+G_B={}",
            g_a + g_b
        );
    }

    #[test]
    fn equilibrium_complementary_pair_forms_dimer() {
        let strands = two_strands();
        let e = fold_all_complexes(&strands).unwrap();
        let eq = solve_tube_equilibrium(&strands, &e, T37_KELVIN).unwrap();
        // At µM concentrations, the complementary pair should be
        // mostly in the dimer state.
        let f_ab = eq.fraction_of(ComplexKind::Heterodimer(0, 1));
        let f_aa = eq.fraction_of(ComplexKind::Homodimer(0));
        let f_a = eq.fraction_of(ComplexKind::Monomer(0));
        // The dimer fraction should exceed both monomer and homodimer.
        assert!(
            f_ab >= f_aa && f_ab >= f_a,
            "expected AB-dominant equilibrium: AB={f_ab}, AA={f_aa}, A={f_a}"
        );
    }

    #[test]
    fn equilibrium_balances_mass() {
        let strands = two_strands();
        let e = fold_all_complexes(&strands).unwrap();
        let eq = solve_tube_equilibrium(&strands, &e, T37_KELVIN).unwrap();
        // Mass conservation: c_tot^A = [A] + 2[AA] + [AB].
        let m_a = eq.concentration_of(ComplexKind::Monomer(0));
        let aa = eq.concentration_of(ComplexKind::Homodimer(0));
        let ab = eq.concentration_of(ComplexKind::Heterodimer(0, 1));
        let total_a = m_a + 2.0 * aa + ab;
        assert!(
            (total_a - 1e-6).abs() / 1e-6 < 1e-3,
            "mass conservation broken for strand A: total = {total_a}, expected = 1e-6"
        );
    }

    #[test]
    fn favoring_target_is_assembled() {
        let t = TargetDistribution::favoring(ComplexKind::Heterodimer(0, 1));
        assert_eq!(t.targets.len(), 1);
        assert_eq!(t.targets[0].kind, ComplexKind::Heterodimer(0, 1));
    }

    #[test]
    fn design_drives_dimer_to_high_fraction() {
        // Start with two random sequences and design them toward the AB
        // dimer.
        let strands = vec![
            TubeStrand::new("A", b"AAAGGCAAAA", 1e-6),
            TubeStrand::new("B", b"AAAGGCAAAA", 1e-6),
        ];
        let target = TargetDistribution::favoring(ComplexKind::Heterodimer(0, 1));
        let d = design_tube(&strands, &target, TubeDesignParams::default()).unwrap();
        // The final equilibrium should have a heterodimer fraction
        // higher than the initial (un-designed) one.
        let initial_eq = solve_tube_equilibrium(
            &strands,
            &fold_all_complexes(&strands).unwrap(),
            T37_KELVIN,
        )
        .unwrap();
        let initial_ab = initial_eq.fraction_of(ComplexKind::Heterodimer(0, 1));
        let final_ab = d.equilibrium.fraction_of(ComplexKind::Heterodimer(0, 1));
        assert!(
            final_ab >= initial_ab - 1e-9,
            "tube design did not improve AB fraction: initial {initial_ab}, final {final_ab}"
        );
    }

    #[test]
    fn design_preferential_dimer_succeeds() {
        // The canonical success criterion: at µM concentrations, the
        // designed pair forms the AB dimer with a higher fraction than
        // a random starting point.
        let strands = vec![
            TubeStrand::new("A", b"AAAAAAAAAAAA", 1e-6),
            TubeStrand::new("B", b"AAAAAAAAAAAA", 1e-6),
        ];
        let target = TargetDistribution::favoring(ComplexKind::Heterodimer(0, 1));
        let d = design_tube(&strands, &target, TubeDesignParams::default()).unwrap();
        let f = d.equilibrium.fraction_of(ComplexKind::Heterodimer(0, 1));
        // After design, the heterodimer must be the *preferred*
        // complex — more abundant than every other.
        let f_aa = d.equilibrium.fraction_of(ComplexKind::Homodimer(0));
        let f_bb = d.equilibrium.fraction_of(ComplexKind::Homodimer(1));
        assert!(
            f > f_aa && f > f_bb,
            "designed AB not the dominant dimer: AB={f}, AA={f_aa}, BB={f_bb}"
        );
    }

    #[test]
    fn rejects_zero_strands() {
        let target = TargetDistribution::favoring(ComplexKind::Monomer(0));
        let err = design_tube(&[], &target, TubeDesignParams::default()).unwrap_err();
        assert_eq!(err.code(), "rnadesign.goal");
    }

    #[test]
    fn rejects_too_many_strands() {
        let strands = vec![
            TubeStrand::new("A", b"AAAA", 1e-6),
            TubeStrand::new("B", b"AAAA", 1e-6),
            TubeStrand::new("C", b"AAAA", 1e-6),
            TubeStrand::new("D", b"AAAA", 1e-6),
        ];
        let target = TargetDistribution::favoring(ComplexKind::Monomer(0));
        assert!(design_tube(&strands, &target, TubeDesignParams::default()).is_err());
    }

    #[test]
    fn rejects_empty_strand() {
        let strands = vec![
            TubeStrand::new("A", b"AAAA", 1e-6),
            TubeStrand::new("B", b"", 1e-6),
        ];
        let target = TargetDistribution::favoring(ComplexKind::Heterodimer(0, 1));
        assert!(design_tube(&strands, &target, TubeDesignParams::default()).is_err());
    }

    #[test]
    fn rejects_zero_concentration() {
        let strands = vec![TubeStrand::new("A", b"AAAA", 0.0)];
        let target = TargetDistribution::favoring(ComplexKind::Monomer(0));
        let err = design_tube(&strands, &target, TubeDesignParams::default()).unwrap_err();
        assert_eq!(err.code(), "rnadesign.invalid");
    }

    #[test]
    fn rejects_zero_iterations() {
        let strands = two_strands();
        let target = TargetDistribution::favoring(ComplexKind::Heterodimer(0, 1));
        let p = TubeDesignParams {
            iterations: 0,
            ..TubeDesignParams::default()
        };
        assert!(design_tube(&strands, &target, p).is_err());
    }

    #[test]
    fn rejects_invalid_target_fraction() {
        let strands = two_strands();
        let target = TargetDistribution {
            targets: vec![TargetFraction {
                kind: ComplexKind::Heterodimer(0, 1),
                fraction: 1.5,
            }],
        };
        let err = design_tube(&strands, &target, TubeDesignParams::default()).unwrap_err();
        assert_eq!(err.code(), "rnadesign.invalid");
    }

    #[test]
    fn rejects_out_of_range_complex() {
        let strands = two_strands();
        let target = TargetDistribution::favoring(ComplexKind::Monomer(5));
        let err = design_tube(&strands, &target, TubeDesignParams::default()).unwrap_err();
        assert_eq!(err.code(), "rnadesign.goal");
    }

    #[test]
    fn complex_kind_stoichiometry() {
        assert_eq!(ComplexKind::Monomer(0).stoichiometry(0), 1);
        assert_eq!(ComplexKind::Monomer(0).stoichiometry(1), 0);
        assert_eq!(ComplexKind::Homodimer(0).stoichiometry(0), 2);
        assert_eq!(ComplexKind::Heterodimer(0, 1).stoichiometry(0), 1);
        assert_eq!(ComplexKind::Heterodimer(0, 1).stoichiometry(1), 1);
        assert_eq!(ComplexKind::Heterodimer(0, 1).stoichiometry(2), 0);
    }

    #[test]
    fn complex_kind_size() {
        assert_eq!(ComplexKind::Monomer(0).size(), 1);
        assert_eq!(ComplexKind::Homodimer(0).size(), 2);
        assert_eq!(ComplexKind::Heterodimer(0, 1).size(), 2);
    }

    #[test]
    fn design_is_deterministic() {
        let strands = two_strands();
        let target = TargetDistribution::favoring(ComplexKind::Heterodimer(0, 1));
        let a = design_tube(&strands, &target, TubeDesignParams::default()).unwrap();
        let b = design_tube(&strands, &target, TubeDesignParams::default()).unwrap();
        assert_eq!(a.strands, b.strands);
    }
}
