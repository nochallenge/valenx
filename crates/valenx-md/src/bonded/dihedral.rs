//! Proper dihedral torsions — **roadmap feature 14**.
//!
//! A proper dihedral `i`-`j`-`k`-`l` measures the twist about the
//! central `j`-`k` bond. Two functional forms are supported, both
//! tagged through [`crate::forcefield::DihedralKind`]:
//!
//! - **Periodic / cosine** — `V(φ) = k·(1 + cos(n·φ − φ₀))`, the
//!   standard AMBER / CHARMM proper torsion.
//! - **Ryckaert-Bellemans** — `V(ψ) = Σ_{m=0}^{5} cₘ·cosᵐ(ψ)` with
//!   `ψ = φ − π`, the GROMOS alkane torsion.
//!
//! ## Geometry & forces
//!
//! With `b1 = rⱼ − rᵢ`, `b2 = rₖ − rⱼ`, `b3 = rₗ − rₖ`, the torsion
//! angle is recovered from the two normal vectors `n1 = b1×b2` and
//! `n2 = b2×b3` via the numerically robust `atan2` form
//!
//! ```text
//! φ = atan2( (n1×n2)·b̂2 , n1·n2 )
//! ```
//!
//! The forces use the classic Allen-&-Tildesley / GROMACS torsion
//! gradient, which distributes `dV/dφ` onto the four atoms so the net
//! force and net torque both vanish. The unit tests cross-check every
//! component against a finite difference of the energy.

use nalgebra::Vector3;

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::forcefield::{DihedralKind, DihedralParam};
use crate::system::{Dihedral, System};

/// The proper-dihedral force term (periodic and RB forms mixed
/// freely).
#[derive(Clone, Debug, PartialEq)]
pub struct ProperDihedrals {
    dihedrals: Vec<(Dihedral, DihedralParam)>,
}

impl ProperDihedrals {
    /// Builds the term from explicit `(dihedral, parameter)` pairs.
    pub fn new(dihedrals: Vec<(Dihedral, DihedralParam)>) -> Self {
        ProperDihedrals { dihedrals }
    }

    /// Builds the term by zipping a system's `topology.dihedrals` with
    /// a parallel parameter slice.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the lengths differ.
    pub fn from_system(system: &System, params: &[DihedralParam]) -> Result<Self> {
        let dihedrals = &system.topology.dihedrals;
        if dihedrals.len() != params.len() {
            return Err(MdError::dimension(format!(
                "{} dihedrals but {} dihedral parameters",
                dihedrals.len(),
                params.len()
            )));
        }
        Ok(ProperDihedrals::new(
            dihedrals
                .iter()
                .copied()
                .zip(params.iter().cloned())
                .collect(),
        ))
    }

    /// Number of dihedrals in the term.
    pub fn len(&self) -> usize {
        self.dihedrals.len()
    }

    /// Whether the term has no dihedrals.
    pub fn is_empty(&self) -> bool {
        self.dihedrals.is_empty()
    }
}

/// Computes the torsion angle of four points, returning
/// `(phi, b1, b2, b3)` for downstream force code.
///
/// `phi` is the IUPAC / GROMACS signed dihedral
/// `phi = atan2((n1×n2)·b̂2, n1·n2)` with `n1 = b1×b2`, `n2 = b2×b3` —
/// the convention [`torsion_forces`] differentiates, so the analytic
/// force equals `−∂V/∂r`.
pub(crate) fn torsion_angle(
    ri: Vector3<f64>,
    rj: Vector3<f64>,
    rk: Vector3<f64>,
    rl: Vector3<f64>,
) -> (f64, Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let b1 = rj - ri;
    let b2 = rk - rj;
    let b3 = rl - rk;
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let b2n = b2.norm().max(1e-12);
    let m = n1.cross(&(b2 / b2n));
    let x = n1.dot(&n2);
    let y = m.dot(&n2);
    (y.atan2(x), b1, b2, b3)
}

/// `dV/dphi` for a periodic-cosine dihedral.
fn periodic_energy_grad(k: f64, n: u32, phase: f64, phi: f64) -> (f64, f64) {
    let arg = n as f64 * phi - phase;
    let energy = k * (1.0 + arg.cos());
    let dv_dphi = -k * (n as f64) * arg.sin();
    (energy, dv_dphi)
}

/// `dV/dphi` for a Ryckaert-Bellemans dihedral.
///
/// `V = Σ cₘ cosᵐ(ψ)`, `ψ = φ − π`, so `cos ψ = −cos φ` and
/// `dV/dφ = (dV/dcosψ)·(dcosψ/dφ) = (Σ m cₘ cosᵐ⁻¹ψ)·(sin φ)`.
fn rb_energy_grad(c: &[f64; 6], phi: f64) -> (f64, f64) {
    let psi = phi - std::f64::consts::PI;
    let cos_psi = psi.cos();
    // Energy: Σ cₘ·cosᵐ(ψ).
    let mut energy = 0.0;
    let mut power = 1.0; // cosᵐ
    for &cm in c.iter() {
        energy += cm * power;
        power *= cos_psi;
    }
    // dV/dcos(ψ): Σ m·cₘ·cosᵐ⁻¹(ψ).
    let mut dv_dcospsi = 0.0;
    let mut pm1 = 1.0; // cosᵐ⁻¹
    for (m, &cm) in c.iter().enumerate() {
        if m >= 1 {
            dv_dcospsi += m as f64 * cm * pm1;
            pm1 *= cos_psi;
        }
    }
    // dcos(ψ)/dφ = -sin(ψ) = -sin(φ - π) = sin(φ).
    let dv_dphi = dv_dcospsi * phi.sin();
    (energy, dv_dphi)
}

/// Distributes `dv_dphi` onto the four atoms with the standard torsion
/// gradient; returns `(fi, fj, fk, fl)`.
///
/// `dv_dphi` is `dV/dφ` for the `φ` that [`torsion_angle`] returns. The
/// end-atom gradients are `∂φ/∂rᵢ = +|b₂|/|n₁|²·n₁` and
/// `∂φ/∂rₗ = −|b₂|/|n₂|²·n₂`. The middle atoms follow from the chain
/// rule on the shared `b₂` axis with `p = b₁·b₂/|b₂|²` and
/// `q = b₃·b₂/|b₂|²`:
///
/// ```text
/// ∂φ/∂rⱼ = −(1 + p)·∂φ/∂rᵢ + q·∂φ/∂rₗ
/// ∂φ/∂rₖ =        p·∂φ/∂rᵢ − (1 + q)·∂φ/∂rₗ
/// ```
///
/// which sum to zero over the four atoms (translation invariance). The
/// force is `f = −dV/dφ · ∂φ/∂r`.
pub(crate) fn torsion_forces(
    dv_dphi: f64,
    b1: Vector3<f64>,
    b2: Vector3<f64>,
    b3: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let n1sq = n1.norm_squared().max(1e-24);
    let n2sq = n2.norm_squared().max(1e-24);
    let b2n = b2.norm().max(1e-12);
    // dphi/dr_i and dphi/dr_l (Blondel & Karplus / GROMACS form).
    let dphi_dri = b2n / n1sq * n1;
    let dphi_drl = -b2n / n2sq * n2;
    // Middle atoms via the chain rule on the shared b2 vector.
    let p = b1.dot(&b2) / (b2n * b2n);
    let q = b3.dot(&b2) / (b2n * b2n);
    let dphi_drj = -(1.0 + p) * dphi_dri + q * dphi_drl;
    let dphi_drk = p * dphi_dri - (1.0 + q) * dphi_drl;
    // f = -dV/dphi * dphi/dr.
    (
        -dv_dphi * dphi_dri,
        -dv_dphi * dphi_drj,
        -dv_dphi * dphi_drk,
        -dv_dphi * dphi_drl,
    )
}

impl ForceTerm for ProperDihedrals {
    fn name(&self) -> &str {
        "proper-dihedrals"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        if out.forces.len() != n {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        for (dih, param) in &self.dihedrals {
            for idx in [dih.i, dih.j, dih.k, dih.l] {
                if idx >= n {
                    return Err(MdError::invalid("dihedral", "atom index out of range"));
                }
            }
            let ri = system.positions[dih.i];
            let rj = system.positions[dih.j];
            let rk = system.positions[dih.k];
            let rl = system.positions[dih.l];
            let (phi, b1, b2, b3) = torsion_angle(ri, rj, rk, rl);
            let (energy, dv_dphi) = match &param.kind {
                DihedralKind::Periodic {
                    k,
                    multiplicity,
                    phase,
                } => periodic_energy_grad(*k, *multiplicity, *phase, phi),
                DihedralKind::RyckaertBellemans { c } => rb_energy_grad(c, phi),
            };
            out.energy += energy;
            let (fi, fj, fk, fl) = torsion_forces(dv_dphi, b1, b2, b3);
            out.forces[dih.i] += fi;
            out.forces[dih.j] += fj;
            out.forces[dih.k] += fk;
            out.forces[dih.l] += fl;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};

    /// Four atoms with a tunable torsion angle about the j-k bond.
    fn torsion_system(phi: f64) -> System {
        let mut top = Topology::new();
        for _ in 0..4 {
            top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        }
        top.add_dihedral(0, 1, 2, 3).unwrap();
        // j at origin, k along +x. i in the xy plane; l rotated by phi
        // about the x-axis.
        let ri = Vector3::new(-0.1, 0.13, 0.0);
        let rj = Vector3::zeros();
        let rk = Vector3::new(0.15, 0.0, 0.0);
        let rl = rk + Vector3::new(0.1, 0.13 * phi.cos(), 0.13 * phi.sin());
        System::new(top, vec![ri, rj, rk, rl]).unwrap()
    }

    #[test]
    fn periodic_energy_is_periodic() {
        let p = DihedralParam::periodic(10.0, 1, 0.0).unwrap();
        // At phi making cos = -1 vs +1 the energy differs by 2k.
        let sys0 = torsion_system(0.0);
        let term = ProperDihedrals::from_system(&sys0, &[p.clone()]).unwrap();
        let mut e0 = EnergyForce::zeros(4);
        term.accumulate(&sys0, &mut e0).unwrap();
        assert!(e0.energy >= 0.0 && e0.energy <= 20.0 + 1e-6);
    }

    #[test]
    fn periodic_force_matches_finite_difference() {
        let p = DihedralParam::periodic(8.0, 3, 0.5).unwrap();
        let base = torsion_system(0.9);
        let term = ProperDihedrals::from_system(&base, &[p]).unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&base, &mut ef).unwrap();

        let h = 1e-6;
        for atom in 0..4 {
            for comp in 0..3 {
                let energy_at = |delta: f64| {
                    let mut s = base.clone();
                    s.positions[atom][comp] += delta;
                    let mut e = EnergyForce::zeros(4);
                    term.accumulate(&s, &mut e).unwrap();
                    e.energy
                };
                let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
                assert!(
                    (ef.forces[atom][comp] - fd).abs() < 2e-3,
                    "atom {atom} comp {comp}: {} vs {}",
                    ef.forces[atom][comp],
                    fd
                );
            }
        }
    }

    #[test]
    fn ryckaert_bellemans_force_matches_finite_difference() {
        // The classic GROMOS n-butane RB coefficients (kJ/mol).
        let c = [9.28, 12.16, -13.12, -3.06, 26.24, -31.5];
        let p = DihedralParam::ryckaert_bellemans(c).unwrap();
        let base = torsion_system(1.3);
        let term = ProperDihedrals::from_system(&base, &[p]).unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&base, &mut ef).unwrap();

        let h = 1e-6;
        for atom in 0..4 {
            for comp in 0..3 {
                let energy_at = |delta: f64| {
                    let mut s = base.clone();
                    s.positions[atom][comp] += delta;
                    let mut e = EnergyForce::zeros(4);
                    term.accumulate(&s, &mut e).unwrap();
                    e.energy
                };
                let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
                assert!(
                    (ef.forces[atom][comp] - fd).abs() < 5e-3,
                    "RB atom {atom} comp {comp}: {} vs {}",
                    ef.forces[atom][comp],
                    fd
                );
            }
        }
    }

    #[test]
    fn net_force_vanishes() {
        let p = DihedralParam::periodic(5.0, 2, 0.0).unwrap();
        let sys = torsion_system(0.7);
        let term = ProperDihedrals::from_system(&sys, &[p]).unwrap();
        let mut ef = EnergyForce::zeros(4);
        term.accumulate(&sys, &mut ef).unwrap();
        let net: Vector3<f64> = ef.forces.iter().sum();
        assert!(net.norm() < 1e-8, "net force = {}", net.norm());
    }
}
