//! Harmonic angle bending — **roadmap feature 13**.
//!
//! Each angle `i`-`j`-`k` (vertex at `j`) carries
//!
//! ```text
//! V(θ) = ½ · k · (θ − θ₀)²
//! ```
//!
//! where `θ` is the angle between the vectors **j→i** and **j→k**.
//!
//! The forces come from the standard analytic gradient. With
//! `rᵢⱼ = rᵢ − rⱼ`, `rₖⱼ = rₖ − rⱼ` and `cos θ = (rᵢⱼ·rₖⱼ)/(|rᵢⱼ||rₖⱼ|)`,
//! the force is `fᵢ = −∂V/∂rᵢ`. Since `θ = arccos(cos θ)` has
//! `∂θ/∂cos θ = −1/sin θ`, the two minus signs cancel and
//!
//! ```text
//! fᵢ = +(k·Δθ / sin θ) · ∂cos θ/∂rᵢ
//! ```
//!
//! and likewise for `k`; the vertex force is `fⱼ = −(fᵢ + fₖ)` so the
//! total force on the angle is zero. The `1/sin θ` is guarded near the
//! 0° / 180° collinear singularity.

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::forcefield::AngleParam;
use crate::system::{Angle, System};

/// The harmonic-angle force term.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicAngles {
    angles: Vec<(Angle, AngleParam)>,
}

impl HarmonicAngles {
    /// Builds the term from explicit `(angle, parameter)` pairs.
    pub fn new(angles: Vec<(Angle, AngleParam)>) -> Self {
        HarmonicAngles { angles }
    }

    /// Builds the term by zipping a system's `topology.angles` with a
    /// parallel parameter slice.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the lengths differ.
    pub fn from_system(system: &System, params: &[AngleParam]) -> Result<Self> {
        let angles = &system.topology.angles;
        if angles.len() != params.len() {
            return Err(MdError::dimension(format!(
                "{} angles but {} angle parameters",
                angles.len(),
                params.len()
            )));
        }
        Ok(HarmonicAngles::new(
            angles.iter().copied().zip(params.iter().copied()).collect(),
        ))
    }

    /// Number of angles in the term.
    pub fn len(&self) -> usize {
        self.angles.len()
    }

    /// Whether the term has no angles.
    pub fn is_empty(&self) -> bool {
        self.angles.is_empty()
    }
}

impl ForceTerm for HarmonicAngles {
    fn name(&self) -> &str {
        "harmonic-angles"
    }

    fn accumulate(&self, system: &System, out: &mut EnergyForce) -> Result<()> {
        let n = system.len();
        if out.forces.len() != n {
            return Err(MdError::dimension(
                "force accumulator size does not match the system",
            ));
        }
        for (angle, param) in &self.angles {
            if angle.i >= n || angle.j >= n || angle.k >= n {
                return Err(MdError::invalid("angle", "atom index out of range"));
            }
            let rij = system
                .cell
                .min_image(system.positions[angle.i] - system.positions[angle.j]);
            let rkj = system
                .cell
                .min_image(system.positions[angle.k] - system.positions[angle.j]);
            let lij = rij.norm();
            let lkj = rkj.norm();
            if lij < 1e-12 || lkj < 1e-12 {
                continue;
            }
            let cos_theta = (rij.dot(&rkj) / (lij * lkj)).clamp(-1.0, 1.0);
            let theta = cos_theta.acos();
            let dtheta = theta - param.theta0;
            out.energy += 0.5 * param.k * dtheta * dtheta;

            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-8);
            // dV/dtheta = k*dtheta.
            let dv_dtheta = param.k * dtheta;
            // f_i = -dV/dr_i = -(dV/dtheta)(dtheta/dcos)(dcos/dr_i).
            // dtheta/dcos = -1/sin(theta), so the two signs cancel:
            // f_i = +(dV/dtheta / sin theta) * dcos/dr_i.
            let prefac = dv_dtheta / sin_theta;
            let fi = prefac * (rkj / (lij * lkj) - cos_theta * rij / (lij * lij));
            let fk = prefac * (rij / (lij * lkj) - cos_theta * rkj / (lkj * lkj));
            let fj = -(fi + fk);
            out.forces[angle.i] += fi;
            out.forces[angle.j] += fj;
            out.forces[angle.k] += fk;
            // Virial: sum over the two bond vectors from the vertex.
            out.virial += rij.dot(&fi) + rkj.dot(&fk);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// Builds an i-j-k system with a right angle at j.
    fn bent_system(angle_rad: f64) -> System {
        let mut top = Topology::new();
        for _ in 0..3 {
            top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        }
        top.add_angle(0, 1, 2).unwrap();
        // j at origin, i along +x, k at the requested angle in xy.
        let pos = vec![
            Vector3::new(0.15, 0.0, 0.0),
            Vector3::zeros(),
            Vector3::new(0.15 * angle_rad.cos(), 0.15 * angle_rad.sin(), 0.0),
        ];
        System::new(top, pos).unwrap()
    }

    #[test]
    fn zero_energy_at_equilibrium() {
        let theta0 = 1.911; // ~109.5 degrees
        let sys = bent_system(theta0);
        let term =
            HarmonicAngles::from_system(&sys, &[AngleParam::new(theta0, 400.0).unwrap()]).unwrap();
        let mut ef = EnergyForce::zeros(3);
        term.accumulate(&sys, &mut ef).unwrap();
        assert!(ef.energy.abs() < 1e-6, "energy = {}", ef.energy);
    }

    #[test]
    fn distorted_angle_has_positive_energy_and_zero_net_force() {
        let sys = bent_system(std::f64::consts::FRAC_PI_2); // 90 deg
        let term = HarmonicAngles::from_system(
            &sys,
            &[AngleParam::new(1.911, 400.0).unwrap()], // wants 109.5
        )
        .unwrap();
        let mut ef = EnergyForce::zeros(3);
        term.accumulate(&sys, &mut ef).unwrap();
        assert!(ef.energy > 0.0);
        let net: Vector3<f64> = ef.forces.iter().sum();
        assert!(net.norm() < 1e-9, "net force = {}", net.norm());
    }

    #[test]
    fn analytic_force_matches_finite_difference() {
        let param = AngleParam::new(1.7, 350.0).unwrap();
        let base = bent_system(2.0);
        let term = HarmonicAngles::from_system(&base, &[param]).unwrap();
        let mut ef = EnergyForce::zeros(3);
        term.accumulate(&base, &mut ef).unwrap();

        let h = 1e-6;
        for atom in 0..3 {
            for comp in 0..3 {
                let energy_at = |delta: f64| {
                    let mut s = base.clone();
                    s.positions[atom][comp] += delta;
                    let mut e = EnergyForce::zeros(3);
                    term.accumulate(&s, &mut e).unwrap();
                    e.energy
                };
                let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
                assert!(
                    (ef.forces[atom][comp] - fd).abs() < 1e-3,
                    "atom {atom} comp {comp}: {} vs {}",
                    ef.forces[atom][comp],
                    fd
                );
            }
        }
    }

    #[test]
    fn len_and_is_empty_report_the_angle_count() {
        let empty = HarmonicAngles::new(Vec::new());
        assert_eq!(empty.len(), 0);
        assert!(empty.is_empty());

        let sys = bent_system(1.911);
        let term =
            HarmonicAngles::from_system(&sys, &[AngleParam::new(1.911, 400.0).unwrap()]).unwrap();
        assert_eq!(term.len(), 1);
        assert!(!term.is_empty());
    }

    #[test]
    fn from_system_rejects_a_parameter_count_mismatch() {
        // The `DimensionMismatch` error branch of `from_system`.
        let sys = bent_system(1.911); // exactly 1 angle
        let two_params = [
            AngleParam::new(1.911, 400.0).unwrap(),
            AngleParam::new(1.911, 400.0).unwrap(),
        ];
        let err = HarmonicAngles::from_system(&sys, &two_params).unwrap_err();
        assert!(
            err.to_string().contains("angle"),
            "error should mention angles: {err}",
        );
    }

    #[test]
    fn accumulate_rejects_a_mis_sized_force_buffer() {
        // The accumulator-size guard in `accumulate`.
        let sys = bent_system(1.911);
        let term =
            HarmonicAngles::from_system(&sys, &[AngleParam::new(1.911, 400.0).unwrap()]).unwrap();
        let mut wrong = EnergyForce::zeros(99); // not 3
        assert!(term.accumulate(&sys, &mut wrong).is_err());
    }

    #[test]
    fn collinear_atoms_are_skipped_without_panicking() {
        // When i and j coincide (zero-length bond vector) the term
        // hits the `lij < 1e-12` guard and skips that angle — no NaN,
        // no panic, zero contribution.
        let mut top = Topology::new();
        for _ in 0..3 {
            top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        }
        top.add_angle(0, 1, 2).unwrap();
        // atom 0 sits exactly on atom 1 (the vertex) -> r_ij = 0.
        let pos = vec![
            Vector3::zeros(),
            Vector3::zeros(),
            Vector3::new(0.15, 0.0, 0.0),
        ];
        let sys = System::new(top, pos).unwrap();
        let term =
            HarmonicAngles::from_system(&sys, &[AngleParam::new(1.911, 400.0).unwrap()]).unwrap();
        let mut ef = EnergyForce::zeros(3);
        term.accumulate(&sys, &mut ef).unwrap();
        assert_eq!(ef.energy, 0.0, "degenerate angle contributes nothing");
        for f in &ef.forces {
            assert!(
                f.norm().is_finite(),
                "no NaN forces from a degenerate angle"
            );
        }
    }
}
