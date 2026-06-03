//! Coarse-grained DNA — oxDNA-class model (**roadmap feature 30a**).
//!
//! All-atom MD cannot reach the microsecond–millisecond scales of DNA
//! hybridisation and nanostructure self-assembly. **Coarse-grained**
//! models replace each nucleotide with a few interaction sites and a
//! purpose-built force field. The **oxDNA** model is the best-known of
//! these.
//!
//! This module implements a *simplified oxDNA-class* model: each
//! nucleotide is a single rigid bead carrying a position and a base
//! identity, joined by four interaction terms:
//!
//! 1. **Backbone (FENE)** — consecutive beads on a strand are linked
//!    by a finitely-extensible (FENE) spring, so the backbone cannot
//!    be stretched past a maximum length.
//! 2. **Excluded volume** — a short-range repulsion (a truncated
//!    Lennard-Jones-like core) keeps non-bonded beads from
//!    overlapping.
//! 3. **Hydrogen bonding** — a distance-gated attractive well between
//!    Watson-Crick-complementary bases (A·T, G·C) on different
//!    strands. This is the term that pairs two strands into a duplex.
//! 4. **Stacking** — a short-range attraction between consecutive
//!    bases *along* a strand, the dominant stabiliser of the helix.
//!
//! ## v1 caveat — honest scope
//!
//! This reproduces the *qualitative physics* — a complementary strand
//! pair lowers its energy by forming hydrogen bonds, the backbone
//! resists overstretching, beads do not interpenetrate. It is **not**
//! the published oxDNA / oxDNA2 parameterisation: oxDNA uses rigid
//! bodies with orientation vectors and an anisotropic, mutually-
//! gated HB+stacking potential fitted to DNA thermodynamics. Here the
//! beads are point particles and the base-pairing term is isotropic
//! and distance-gated only. It is a working coarse-grained DNA force
//! field for the engine, with the orientational refinement documented
//! as the next step.

use nalgebra::Vector3;

use crate::error::{MdError, Result};

/// A DNA base identity.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Base {
    /// Adenine.
    A,
    /// Thymine.
    T,
    /// Guanine.
    G,
    /// Cytosine.
    C,
}

impl Base {
    /// Parses a single-letter base code (case-insensitive).
    ///
    /// # Errors
    /// [`MdError::Parse`] on any character other than A/T/G/C.
    pub fn from_char(c: char) -> Result<Self> {
        match c.to_ascii_uppercase() {
            'A' => Ok(Base::A),
            'T' => Ok(Base::T),
            'G' => Ok(Base::G),
            'C' => Ok(Base::C),
            other => Err(MdError::parse(
                "dna-sequence",
                format!("`{other}` is not a DNA base (A/T/G/C)"),
            )),
        }
    }

    /// The Watson-Crick complement.
    pub fn complement(self) -> Base {
        match self {
            Base::A => Base::T,
            Base::T => Base::A,
            Base::G => Base::C,
            Base::C => Base::G,
        }
    }

    /// Whether two bases form a Watson-Crick pair.
    pub fn pairs_with(self, other: Base) -> bool {
        self.complement() == other
    }
}

/// One coarse-grained nucleotide bead.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DnaBead {
    /// Bead position (nm).
    pub position: Vector3<f64>,
    /// Base identity.
    pub base: Base,
    /// Index of the strand this bead belongs to.
    pub strand: usize,
}

/// Parameters of the simplified oxDNA-class force field. The defaults
/// are a self-consistent set in the crate's nm / kJ-mol units; they
/// are *representative*, not the fitted oxDNA constants.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct OxdnaParams {
    /// FENE backbone equilibrium length r₀ (nm).
    pub backbone_r0: f64,
    /// FENE backbone stiffness (kJ/mol/nm²).
    pub backbone_k: f64,
    /// FENE maximum extension Δ from r₀ (nm).
    pub backbone_delta: f64,
    /// Excluded-volume radius σ (nm).
    pub excluded_sigma: f64,
    /// Excluded-volume strength ε (kJ/mol).
    pub excluded_epsilon: f64,
    /// Hydrogen-bond well depth (kJ/mol).
    pub hbond_strength: f64,
    /// Hydrogen-bond equilibrium distance (nm).
    pub hbond_r0: f64,
    /// Hydrogen-bond range (nm) — the Gaussian well width.
    pub hbond_range: f64,
    /// Stacking well depth (kJ/mol).
    pub stacking_strength: f64,
    /// Stacking equilibrium distance (nm).
    pub stacking_r0: f64,
    /// Stacking range (nm).
    pub stacking_range: f64,
}

impl Default for OxdnaParams {
    fn default() -> Self {
        OxdnaParams {
            backbone_r0: 0.34,
            backbone_k: 2000.0,
            backbone_delta: 0.15,
            excluded_sigma: 0.30,
            excluded_epsilon: 5.0,
            hbond_strength: 12.0,
            hbond_r0: 0.34,
            hbond_range: 0.12,
            stacking_strength: 8.0,
            stacking_r0: 0.36,
            stacking_range: 0.10,
        }
    }
}

/// A coarse-grained DNA system: a set of beads and the force-field
/// parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct CoarseDna {
    /// The nucleotide beads. Beads on the same strand are assumed
    /// consecutive (5'→3' order) within the slice.
    pub beads: Vec<DnaBead>,
    /// Force-field parameters.
    pub params: OxdnaParams,
}

impl CoarseDna {
    /// Builds a coarse-grained strand from a sequence and a per-base
    /// position list, on strand index `strand`.
    ///
    /// # Errors
    /// [`MdError::Parse`] on a bad base character;
    /// [`MdError::DimensionMismatch`] if the sequence length and the
    /// position count differ.
    pub fn from_strand(
        sequence: &str,
        positions: &[Vector3<f64>],
        strand: usize,
    ) -> Result<Self> {
        let bases: Vec<Base> = sequence
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(Base::from_char)
            .collect::<Result<_>>()?;
        if bases.len() != positions.len() {
            return Err(MdError::dimension(format!(
                "{} bases but {} positions",
                bases.len(),
                positions.len()
            )));
        }
        let beads = bases
            .into_iter()
            .zip(positions)
            .map(|(base, &position)| DnaBead {
                position,
                base,
                strand,
            })
            .collect();
        Ok(CoarseDna {
            beads,
            params: OxdnaParams::default(),
        })
    }

    /// Appends another strand's beads to this system.
    pub fn add_strand(&mut self, other: &CoarseDna) {
        self.beads.extend_from_slice(&other.beads);
    }

    /// Number of beads.
    pub fn len(&self) -> usize {
        self.beads.len()
    }

    /// Whether the system has no beads.
    pub fn is_empty(&self) -> bool {
        self.beads.is_empty()
    }

    /// The total coarse-grained potential energy (kJ/mol) — the sum of
    /// the backbone, excluded-volume, hydrogen-bond and stacking
    /// terms.
    pub fn energy(&self) -> f64 {
        self.energy_breakdown().total()
    }

    /// The energy split into its four contributions.
    pub fn energy_breakdown(&self) -> OxdnaEnergy {
        let p = &self.params;
        let n = self.beads.len();
        let mut backbone = 0.0;
        let mut excluded = 0.0;
        let mut hbond = 0.0;
        let mut stacking = 0.0;

        for i in 0..n {
            for j in (i + 1)..n {
                let bi = self.beads[i];
                let bj = self.beads[j];
                let r = (bi.position - bj.position).norm();
                let same_strand = bi.strand == bj.strand;
                let consecutive = same_strand && (j == i + 1);

                if consecutive {
                    // FENE backbone spring.
                    backbone += fene_energy(r, p.backbone_r0, p.backbone_k, p.backbone_delta);
                    // Consecutive bases also stack.
                    stacking -= gaussian_well(
                        r,
                        p.stacking_strength,
                        p.stacking_r0,
                        p.stacking_range,
                    );
                } else {
                    // Excluded volume for every non-bonded pair.
                    excluded += excluded_volume(r, p.excluded_sigma, p.excluded_epsilon);
                    // Hydrogen bonding only between complementary
                    // bases on *different* strands.
                    if !same_strand && bi.base.pairs_with(bj.base) {
                        hbond -= gaussian_well(
                            r,
                            p.hbond_strength,
                            p.hbond_r0,
                            p.hbond_range,
                        );
                    }
                }
            }
        }
        OxdnaEnergy {
            backbone,
            excluded,
            hbond,
            stacking,
        }
    }

    /// Counts the formed Watson-Crick base pairs — complementary
    /// inter-strand beads within `1.5·hbond_r0` of each other.
    pub fn count_base_pairs(&self) -> usize {
        let p = &self.params;
        let cutoff = 1.5 * p.hbond_r0;
        let mut pairs = 0;
        let n = self.beads.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let bi = self.beads[i];
                let bj = self.beads[j];
                if bi.strand != bj.strand
                    && bi.base.pairs_with(bj.base)
                    && (bi.position - bj.position).norm() < cutoff
                {
                    pairs += 1;
                }
            }
        }
        pairs
    }
}

/// The four-way split of the coarse-grained DNA energy.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct OxdnaEnergy {
    /// FENE backbone energy (kJ/mol).
    pub backbone: f64,
    /// Excluded-volume energy (kJ/mol).
    pub excluded: f64,
    /// Hydrogen-bond energy (kJ/mol, negative when paired).
    pub hbond: f64,
    /// Stacking energy (kJ/mol, negative when stacked).
    pub stacking: f64,
}

impl OxdnaEnergy {
    /// The total energy.
    pub fn total(&self) -> f64 {
        self.backbone + self.excluded + self.hbond + self.stacking
    }
}

/// FENE finitely-extensible spring energy.
///
/// `V = −½·k·Δ²·ln(1 − ((r−r₀)/Δ)²)`. Diverges as `r` approaches
/// `r₀ ± Δ`; outside that window the bond is "broken" and a large
/// finite penalty is returned so the optimisation stays well-behaved.
fn fene_energy(r: f64, r0: f64, k: f64, delta: f64) -> f64 {
    let x = (r - r0) / delta;
    if x.abs() >= 0.999 {
        // Past the extensibility limit: large but finite.
        return 0.5 * k * delta * delta * 10.0;
    }
    -0.5 * k * delta * delta * (1.0 - x * x).ln()
}

/// A purely-repulsive excluded-volume core: a truncated-shifted
/// Lennard-Jones repulsion, zero beyond `σ`.
fn excluded_volume(r: f64, sigma: f64, epsilon: f64) -> f64 {
    if r >= sigma || r < 1e-9 {
        return 0.0;
    }
    let sr6 = (sigma / r).powi(6);
    // WCA-style: 4ε(sr12 − sr6) + ε, shifted to zero at r = σ.
    4.0 * epsilon * (sr6 * sr6 - sr6) + epsilon
}

/// A Gaussian attractive well of depth `strength` centred at `r0` with
/// width `range`. Returned positive; callers subtract it for an
/// attractive contribution.
fn gaussian_well(r: f64, strength: f64, r0: f64, range: f64) -> f64 {
    let x = (r - r0) / range;
    strength * (-0.5 * x * x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_complement_and_pairing() {
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::G.complement(), Base::C);
        assert!(Base::A.pairs_with(Base::T));
        assert!(Base::C.pairs_with(Base::G));
        assert!(!Base::A.pairs_with(Base::G));
    }

    #[test]
    fn base_parsing() {
        assert_eq!(Base::from_char('a').unwrap(), Base::A);
        assert_eq!(Base::from_char('G').unwrap(), Base::G);
        assert!(Base::from_char('X').unwrap_err().category() == "parse");
    }

    /// Builds a straight strand along x with the backbone spacing.
    fn straight_strand(seq: &str, strand: usize, offset: Vector3<f64>) -> CoarseDna {
        let spacing = OxdnaParams::default().backbone_r0;
        let pos: Vec<Vector3<f64>> = (0..seq.chars().count())
            .map(|i| offset + Vector3::new(i as f64 * spacing, 0.0, 0.0))
            .collect();
        CoarseDna::from_strand(seq, &pos, strand).unwrap()
    }

    #[test]
    fn from_strand_validates_lengths() {
        let pos = vec![Vector3::zeros(); 2];
        assert!(CoarseDna::from_strand("ATG", &pos, 0).is_err());
        assert!(CoarseDna::from_strand("AT", &pos, 0).is_ok());
    }

    #[test]
    fn backbone_at_equilibrium_costs_little() {
        let dna = straight_strand("ATGC", 0, Vector3::zeros());
        let e = dna.energy_breakdown();
        // FENE at r0 is exactly zero.
        assert!(e.backbone.abs() < 1e-6, "backbone = {}", e.backbone);
    }

    #[test]
    fn complementary_strands_form_base_pairs() {
        // Two antiparallel-ish complementary strands placed a HB
        // distance apart pair up.
        let hb = OxdnaParams::default().hbond_r0;
        let strand_a = straight_strand("ATGC", 0, Vector3::zeros());
        // Complement of ATGC is TACG; lay it directly across at the HB
        // distance in y. (Orientation is ignored in this v1.)
        let strand_b = straight_strand("TACG", 1, Vector3::new(0.0, hb, 0.0));
        let mut duplex = strand_a;
        duplex.add_strand(&strand_b);
        // Every base across should pair.
        assert_eq!(duplex.count_base_pairs(), 4);
        // Hydrogen bonding lowers the energy.
        assert!(duplex.energy_breakdown().hbond < 0.0);
    }

    #[test]
    fn non_complementary_strands_do_not_pair() {
        let hb = OxdnaParams::default().hbond_r0;
        let strand_a = straight_strand("AAAA", 0, Vector3::zeros());
        // Place a non-complementary strand across.
        let strand_b = straight_strand("AAAA", 1, Vector3::new(0.0, hb, 0.0));
        let mut system = strand_a;
        system.add_strand(&strand_b);
        assert_eq!(system.count_base_pairs(), 0);
        assert!(system.energy_breakdown().hbond.abs() < 1e-9);
    }

    #[test]
    fn consecutive_bases_stack() {
        let dna = straight_strand("GGGG", 0, Vector3::zeros());
        // Stacking is an attractive (negative) contribution.
        assert!(dna.energy_breakdown().stacking < 0.0);
    }

    #[test]
    fn overlapping_beads_have_excluded_volume_penalty() {
        // Two beads from different strands placed almost on top of
        // each other.
        let strand_a = CoarseDna::from_strand("A", &[Vector3::zeros()], 0).unwrap();
        let strand_b =
            CoarseDna::from_strand("A", &[Vector3::new(0.05, 0.0, 0.0)], 1).unwrap();
        let mut system = strand_a;
        system.add_strand(&strand_b);
        assert!(system.energy_breakdown().excluded > 0.0);
    }

    #[test]
    fn stretched_backbone_is_penalised() {
        let p = OxdnaParams::default();
        // A bond stretched near its FENE limit costs a lot.
        let e_eq = fene_energy(p.backbone_r0, p.backbone_r0, p.backbone_k, p.backbone_delta);
        let e_stretched = fene_energy(
            p.backbone_r0 + 0.99 * p.backbone_delta,
            p.backbone_r0,
            p.backbone_k,
            p.backbone_delta,
        );
        assert!(e_stretched > e_eq + 10.0);
        assert!(e_stretched.is_finite());
    }
}
