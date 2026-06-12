//! Ewald / particle-mesh Ewald — **roadmap feature 9**.
//!
//! A cut-off Coulomb sum (see [`super::coulomb`]) is cheap but
//! approximate. **Ewald summation** computes the *exact* electrostatic
//! energy of a periodic system by splitting the conditionally-
//! convergent `Σ 1/r` lattice sum into two absolutely-convergent
//! pieces, each cheap in its own space:
//!
//! ```text
//! E = E_real + E_recip − E_self  ( + E_excluded correction )
//! ```
//!
//! - **Real space** — each point charge is screened by a Gaussian of
//!   width `1/β`. The screened interaction
//!   `qᵢqⱼ·erfc(β·r)/r` decays fast, so it is summed only within a
//!   cutoff. `β` is the Ewald splitting parameter.
//! - **Reciprocal space** — the smooth compensating Gaussian charge
//!   cloud is summed over reciprocal-lattice vectors **k** as a
//!   structure factor `S(k) = Σ qⱼ·exp(i·k·rⱼ)`:
//!
//!   ```text
//!   E_recip = (1/2V·ε₀') Σ_{k≠0} (4π/k²)·exp(−k²/4β²)·|S(k)|²
//!   ```
//! - **Self term** — removes each charge's spurious interaction with
//!   its own screening Gaussian: `E_self = (β/√π)·Σ qⱼ²`.
//! - **Excluded correction** — subtracts the reciprocal-space
//!   contribution of 1-2 / 1-3 pairs whose direct interaction is
//!   carried by the bonded terms.
//!
//! ## v1 caveat — this is *Ewald*, not a mesh
//!
//! A true particle-mesh-Ewald (smooth-PME) interpolates the charges
//! onto a grid and evaluates `E_recip` with an FFT, giving
//! `O(N log N)`. This v1 evaluates the reciprocal sum **directly** as
//! a structure-factor double loop over `N` charges and `K` k-vectors —
//! `O(N·K)`. The *result is identical* to smooth-PME within the
//! k-space cutoff (it is the exact Ewald reciprocal sum); it is simply
//! slower. The public type is still called [`Pme`] because it occupies
//! the PME slot in a force field and is interchangeable with one; the
//! mesh/FFT acceleration is the documented future step.
//!
//! Both real- and reciprocal-space forces are analytic and are
//! cross-checked against an energy finite difference in the tests.

use std::f64::consts::PI;

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::nonbonded::ExclusionSet;
use crate::system::System;
use crate::units::COULOMB;

/// `erf` / `erfc` via the Abramowitz-&-Stegun 7.1.26 rational
/// approximation (max abs error ~1.5e-7) — enough for an MD force
/// field and dependency-free.
fn erfc(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * ax);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-ax * ax).exp();
    // y is erf(ax); erf is odd.
    1.0 - sign * y
}

/// The Ewald / PME electrostatics force term.
#[derive(Clone, Debug, PartialEq)]
pub struct Pme {
    /// Per-atom partial charges (e).
    charges: Vec<f64>,
    /// Real-space cutoff (nm).
    real_cutoff: f64,
    /// Ewald splitting parameter β (1/nm).
    beta: f64,
    /// Largest reciprocal-lattice index summed in each dimension.
    k_max: i32,
    /// Excluded 1-2 / 1-3 pairs (their reciprocal contribution is
    /// corrected away).
    exclusions: ExclusionSet,
}

impl Pme {
    /// Builds an Ewald term, choosing β and `k_max` automatically from
    /// the real-space cutoff and a target relative accuracy.
    ///
    /// The standard recipe: pick β so the real-space sum is converged
    /// at `real_cutoff` to `accuracy`, i.e. `β = √(−ln accuracy)/r_c`;
    /// then pick `k_max` so the reciprocal sum matches.
    ///
    /// # Errors
    /// [`MdError::Invalid`] for a non-positive cutoff or an accuracy
    /// outside `(0, 1)`; [`MdError::Invalid`] if the system is not
    /// periodic (Ewald is only defined under PBC).
    pub fn from_system(system: &System, real_cutoff: f64, accuracy: f64) -> Result<Self> {
        if !(real_cutoff.is_finite() && real_cutoff > 0.0) {
            return Err(MdError::invalid(
                "real_cutoff",
                "must be finite and positive",
            ));
        }
        if !(accuracy.is_finite() && accuracy > 0.0 && accuracy < 1.0) {
            return Err(MdError::invalid("accuracy", "must lie in (0, 1)"));
        }
        if !system.cell.is_periodic() {
            return Err(MdError::invalid(
                "cell",
                "Ewald summation requires a periodic box",
            ));
        }
        let beta = (-accuracy.ln()).sqrt() / real_cutoff;
        // Reciprocal-space convergence: choose k_max so the largest
        // omitted Gaussian factor exp(-k^2/4β^2) is below `accuracy`.
        // k for index n is ~2π n / L; solve for n.
        let l_min = system
            .cell
            .edge_lengths()
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let k_needed = 2.0 * beta * (-accuracy.ln()).sqrt();
        let k_max = ((k_needed * l_min / (2.0 * PI)).ceil() as i32).clamp(1, 12);
        Ok(Pme {
            charges: system.topology.atoms.iter().map(|a| a.charge).collect(),
            real_cutoff,
            beta,
            k_max,
            exclusions: ExclusionSet::from_topology(&system.topology),
        })
    }

    /// Builds an Ewald term with explicit β and `k_max` (for tests and
    /// for callers who want full control).
    ///
    /// # Errors
    /// [`MdError::Invalid`] on non-positive parameters or a
    /// non-periodic box.
    pub fn with_parameters(
        system: &System,
        real_cutoff: f64,
        beta: f64,
        k_max: i32,
    ) -> Result<Self> {
        if !(real_cutoff.is_finite() && real_cutoff > 0.0) {
            return Err(MdError::invalid(
                "real_cutoff",
                "must be finite and positive",
            ));
        }
        if !(beta.is_finite() && beta > 0.0) {
            return Err(MdError::invalid("beta", "must be finite and positive"));
        }
        if k_max < 1 {
            return Err(MdError::invalid("k_max", "must be at least 1"));
        }
        if !system.cell.is_periodic() {
            return Err(MdError::invalid(
                "cell",
                "Ewald summation requires a periodic box",
            ));
        }
        Ok(Pme {
            charges: system.topology.atoms.iter().map(|a| a.charge).collect(),
            real_cutoff,
            beta,
            k_max,
            exclusions: ExclusionSet::from_topology(&system.topology),
        })
    }

    /// Replaces the exclusion set.
    pub fn with_exclusions(mut self, exclusions: ExclusionSet) -> Self {
        self.exclusions = exclusions;
        self
    }

    /// The Ewald splitting parameter β (1/nm).
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// The reciprocal-space index range.
    pub fn k_max(&self) -> i32 {
        self.k_max
    }

    /// The screened **real-space** energy + forces.
    ///
    /// `Σ qᵢqⱼ·erfc(β·r)/r` over pairs within `real_cutoff`, with the
    /// analytic gradient. Excluded pairs are skipped here *and* their
    /// full-Coulomb contribution is removed in [`reciprocal`] so the
    /// bonded terms own them.
    fn real_space(
        &self,
        system: &System,
        pairs: &[(usize, usize)],
        out: &mut EnergyForce,
    ) -> Result<()> {
        let n = system.len();
        let rc2 = self.real_cutoff * self.real_cutoff;
        let two_beta_over_sqrtpi = 2.0 * self.beta / PI.sqrt();
        for &(i, j) in pairs {
            if i >= n || j >= n {
                return Err(MdError::dimension("pair index out of range"));
            }
            if self.exclusions.contains(i, j) {
                continue;
            }
            let qq = COULOMB * self.charges[i] * self.charges[j];
            if qq == 0.0 {
                continue;
            }
            let d = system
                .cell
                .min_image(system.positions[i] - system.positions[j]);
            let r2 = d.norm_squared();
            if r2 >= rc2 || r2 < 1e-24 {
                continue;
            }
            let r = r2.sqrt();
            let inv_r = 1.0 / r;
            let br = self.beta * r;
            let erfc_br = erfc(br);
            out.energy += qq * erfc_br * inv_r;
            // dV/dr = -qq[ erfc(βr)/r² + 2β/√π·exp(-β²r²)/r ].
            let dv_dr =
                -qq * (erfc_br * inv_r * inv_r + two_beta_over_sqrtpi * (-br * br).exp() * inv_r);
            let fij = -(dv_dr * inv_r) * d;
            out.forces[i] += fij;
            out.forces[j] -= fij;
            out.virial += d.dot(&fij);
        }
        Ok(())
    }

    /// The **reciprocal-space** energy + forces plus the self and
    /// excluded-pair corrections.
    fn reciprocal(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        let cell = &system.cell;
        let volume = cell.volume();
        if !volume.is_finite() || volume <= 0.0 {
            return Err(MdError::invalid("cell", "needs a finite positive volume"));
        }
        // Reciprocal lattice: rows of 2π·(h⁻¹)ᵀ. For a general cell
        // the reciprocal vectors are 2π·(h⁻¹)ᵀ columns; build them
        // from the inverse-transpose.
        let h = cell.matrix();
        let h_inv = h
            .try_inverse()
            .ok_or_else(|| MdError::invalid("cell", "lattice matrix not invertible"))?;
        let recip = h_inv.transpose() * (2.0 * PI);
        let ra = recip.column(0).into_owned();
        let rb = recip.column(1).into_owned();
        let rc = recip.column(2).into_owned();

        let four_beta_sq = 4.0 * self.beta * self.beta;

        // Loop k-vectors. Use the half-space {kx>0} ∪ {kx=0,ky>0} ∪
        // {kx=0,ky=0,kz>0} and double, to avoid summing k and -k.
        //
        // Per-k energy weight (Gaussian units, before the COULOMB
        // prefactor): the full-sphere Ewald term is
        // (1/2V)·(4π/k²)·exp(−k²/4β²)·|S(k)|². Summing only the
        // half-space and doubling folds the leading ½ away, leaving
        //   ak = (4π / (V·k²)) · exp(−k²/4β²).
        for kx in 0..=self.k_max {
            let ymin = if kx == 0 { 0 } else { -self.k_max };
            for ky in ymin..=self.k_max {
                let zmin = if kx == 0 && ky == 0 { 1 } else { -self.k_max };
                for kz in zmin..=self.k_max {
                    let kvec = kx as f64 * ra + ky as f64 * rb + kz as f64 * rc;
                    let k2 = kvec.norm_squared();
                    if k2 < 1e-12 {
                        continue;
                    }
                    let gauss = (-k2 / four_beta_sq).exp();
                    if gauss < 1e-12 {
                        continue;
                    }
                    let ak = 4.0 * PI / (volume * k2) * gauss;
                    // Compute the structure factor S(k) = Σ q exp(i k·r).
                    let mut s_re = 0.0;
                    let mut s_im = 0.0;
                    for a in 0..n {
                        let phase = kvec.dot(&system.positions[a]);
                        s_re += self.charges[a] * phase.cos();
                        s_im += self.charges[a] * phase.sin();
                    }
                    let s2 = s_re * s_re + s_im * s_im;
                    out.energy += ak * s2 * COULOMB;
                    // Force on atom a:
                    // f_a = q_a · 2·ak·COULOMB · k · [sin(k·r_a)·S_re − cos(k·r_a)·S_im].
                    for a in 0..n {
                        let phase = kvec.dot(&system.positions[a]);
                        let (sp, cp) = phase.sin_cos();
                        let coeff = self.charges[a] * 2.0 * ak * COULOMB * (sp * s_re - cp * s_im);
                        out.forces[a] += coeff * kvec;
                    }
                }
            }
        }

        // Self-energy correction: -β/√π · Σ q² · COULOMB.
        let self_e: f64 = self.charges.iter().map(|q| q * q).sum::<f64>();
        out.energy -= self.beta / PI.sqrt() * self_e * COULOMB;

        // Net-charge ("background plasma") correction. For a non-
        // neutral cell Ewald adds a uniform neutralising background;
        // its energy is −π/(2 V β²)·(Σq)²·COULOMB.
        let net: f64 = self.charges.iter().sum();
        if net.abs() > 1e-12 {
            out.energy -= PI / (2.0 * volume * self.beta * self.beta) * net * net * COULOMB;
        }

        // Excluded-pair correction: the reciprocal sum implicitly
        // includes the *full* 1/r interaction of every pair, including
        // the bonded 1-2/1-3 pairs. Subtract the part of it those
        // pairs should not have — the smooth erf(βr)/r piece. The
        // correction energy is E = −qq·erf(βr)/r.
        let two_beta_over_sqrtpi = 2.0 * self.beta / PI.sqrt();
        for &(i, j) in self.exclusion_pairs().iter() {
            let qq = COULOMB * self.charges[i] * self.charges[j];
            if qq == 0.0 {
                continue;
            }
            let d = cell.min_image(system.positions[i] - system.positions[j]);
            let r = d.norm();
            if r < 1e-12 {
                continue;
            }
            let inv_r = 1.0 / r;
            let br = self.beta * r;
            let erf_br = 1.0 - erfc(br);
            out.energy -= qq * erf_br * inv_r;
            // E = −qq·erf(βr)/r, so
            // dE/dr = −qq·[ (2β/√π)·exp(−β²r²)/r − erf(βr)/r² ].
            let de_dr =
                -qq * (two_beta_over_sqrtpi * (-br * br).exp() * inv_r - erf_br * inv_r * inv_r);
            let fij = -(de_dr * inv_r) * d;
            out.forces[i] += fij;
            out.forces[j] -= fij;
        }
        Ok(())
    }

    /// The excluded pairs as a flat vector (the [`ExclusionSet`] keeps
    /// them in a hash set; reciprocal-space correction needs to
    /// iterate them).
    fn exclusion_pairs(&self) -> Vec<(usize, usize)> {
        let n = self.charges.len();
        let mut v = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if self.exclusions.contains(i, j) {
                    v.push((i, j));
                }
            }
        }
        v
    }

    /// Evaluates the complete Ewald energy + forces over an explicit
    /// real-space pair list.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] on a size mismatch.
    pub fn accumulate_pairs(
        &self,
        system: &System,
        real_pairs: &[(usize, usize)],
        out: &mut EnergyForce,
    ) -> Result<()> {
        if out.forces.len() != system.len() {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        self.real_space(system, real_pairs, out)?;
        self.reciprocal(system, out)?;
        Ok(())
    }
}

impl ForceTerm for Pme {
    fn name(&self) -> &str {
        "ewald-pme"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        let mut pairs = Vec::with_capacity(n * n.saturating_sub(1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                pairs.push((i, j));
            }
        }
        self.accumulate_pairs(system, &pairs, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    #[test]
    fn erfc_endpoints() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        assert!(erfc(3.0) < 1e-3);
        assert!((erfc(-1.0) - (2.0 - erfc(1.0))).abs() < 1e-6);
    }

    fn nacl_pair() -> System {
        // A +1 / -1 pair in a small cubic box.
        let mut top = Topology::new();
        top.push_atom(Atom::new("Na", 23.0, 1.0).unwrap());
        top.push_atom(Atom::new("Cl", 35.0, -1.0).unwrap());
        System::new(
            top,
            vec![Vector3::new(0.5, 0.5, 0.5), Vector3::new(1.5, 0.5, 0.5)],
        )
        .unwrap()
        .with_cell(SimBox::cubic(3.0).unwrap())
    }

    #[test]
    fn rejects_nonperiodic_and_bad_params() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 1.0).unwrap());
        let open = System::new(top, vec![Vector3::zeros()]).unwrap();
        assert!(Pme::from_system(&open, 1.0, 1e-5).is_err());

        let sys = nacl_pair();
        assert!(Pme::from_system(&sys, -1.0, 1e-5).is_err());
        assert!(Pme::from_system(&sys, 1.0, 2.0).is_err());
        assert!(Pme::with_parameters(&sys, 1.0, 0.0, 4).is_err());
        assert!(Pme::with_parameters(&sys, 1.0, 3.0, 0).is_err());
    }

    #[test]
    fn ewald_energy_is_negative_for_opposite_charges() {
        let sys = nacl_pair();
        let pme = Pme::with_parameters(&sys, 1.4, 3.0, 6).unwrap();
        let mut ef = EnergyForce::zeros(2);
        pme.accumulate(&sys, &mut ef).unwrap();
        // Opposite charges -> the Madelung-like energy is negative.
        assert!(ef.energy < 0.0, "energy = {}", ef.energy);
        assert!(ef.energy.is_finite());
    }

    #[test]
    fn ewald_is_independent_of_beta() {
        // The Ewald total energy must not depend on the (arbitrary)
        // splitting parameter — the real and reciprocal parts trade
        // off. Check two β values agree.
        let sys = nacl_pair();
        let mut e1 = EnergyForce::zeros(2);
        Pme::with_parameters(&sys, 1.4, 2.5, 8)
            .unwrap()
            .accumulate(&sys, &mut e1)
            .unwrap();
        let mut e2 = EnergyForce::zeros(2);
        Pme::with_parameters(&sys, 1.4, 4.0, 10)
            .unwrap()
            .accumulate(&sys, &mut e2)
            .unwrap();
        assert!(
            (e1.energy - e2.energy).abs() < 0.05 * e1.energy.abs().max(1.0),
            "beta dependence: {} vs {}",
            e1.energy,
            e2.energy
        );
    }

    #[test]
    fn ewald_force_matches_finite_difference() {
        let base = nacl_pair();
        let pme = Pme::with_parameters(&base, 1.4, 3.0, 8).unwrap();
        let mut ef = EnergyForce::zeros(2);
        pme.accumulate(&base, &mut ef).unwrap();

        let h = 1e-6;
        for comp in 0..3 {
            let energy_at = |delta: f64| {
                let mut s = base.clone();
                s.positions[0][comp] += delta;
                let mut e = EnergyForce::zeros(2);
                pme.accumulate(&s, &mut e).unwrap();
                e.energy
            };
            let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
            assert!(
                (ef.forces[0][comp] - fd).abs() < 5e-2,
                "comp {comp}: {} vs {}",
                ef.forces[0][comp],
                fd
            );
        }
    }

    #[test]
    fn auto_parameters_are_sane() {
        let sys = nacl_pair();
        let pme = Pme::from_system(&sys, 1.2, 1e-5).unwrap();
        assert!(pme.beta() > 0.0);
        assert!((1..=12).contains(&pme.k_max()));
    }
}
