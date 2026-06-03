//! The velocity-Verlet integrator — the standard time-reversible,
//! energy-conserving molecular-dynamics integrator. Pure mechanics: it
//! knows nothing about quantum chemistry, so it is unit-tested in
//! isolation against an analytic harmonic oscillator.
//!
//! Everything is in atomic units (see [`crate::units`]): positions in
//! bohr, velocities in bohr / atomic-time, masses in electron-masses,
//! forces in hartree / bohr, `dt` in atomic-time. Then `a = F/m` is
//! bohr / atomic-time² and the kinetic energy comes out in hartree.

/// Advance one velocity-Verlet step, in place.
///
/// `forces` are the forces at the *current* positions (the previous
/// step's freshly-computed forces); `forces_at(positions) -> forces` is
/// called exactly once, for the new positions, and may **fail** (the
/// AIMD backend's force evaluation is an SCF that can not converge), so
/// it returns a `Result` whose error is propagated. The new forces are
/// returned so the caller reuses them next step — one force evaluation
/// per step. The scheme is half-kick → drift → recompute force →
/// half-kick.
pub fn velocity_verlet_step<E>(
    pos: &mut [[f64; 3]],
    vel: &mut [[f64; 3]],
    forces: &[[f64; 3]],
    masses: &[f64],
    dt: f64,
    forces_at: impl Fn(&[[f64; 3]]) -> core::result::Result<Vec<[f64; 3]>, E>,
) -> core::result::Result<Vec<[f64; 3]>, E> {
    let n = pos.len();
    // First half-kick, then drift.
    for i in 0..n {
        let inv_m = 1.0 / masses[i];
        for d in 0..3 {
            vel[i][d] += 0.5 * forces[i][d] * inv_m * dt;
            pos[i][d] += vel[i][d] * dt;
        }
    }
    // Recompute forces at the new positions (fallible).
    let new_forces = forces_at(pos)?;
    // Second half-kick with the new forces.
    for i in 0..n {
        let inv_m = 1.0 / masses[i];
        for d in 0..3 {
            vel[i][d] += 0.5 * new_forces[i][d] * inv_m * dt;
        }
    }
    Ok(new_forces)
}

/// Total kinetic energy `Σ ½ mᵢ|vᵢ|²` in hartree (atomic units).
pub fn kinetic_energy(vel: &[[f64; 3]], masses: &[f64]) -> f64 {
    vel.iter()
        .zip(masses)
        .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 1-D harmonic oscillator `F(x) = -k·x` embedded along the x axis.
    /// Analytic: period `T = 2π√(m/k)`, total energy `½mv² + ½kx²`
    /// constant. velocity-Verlet must (a) conserve energy to a tight
    /// tolerance and (b) return to the start after exactly one period.
    /// The force closure is infallible here (error type `()`).
    #[test]
    fn harmonic_oscillator_conserves_energy_over_a_period() {
        let (k, m) = (1.0_f64, 1.0_f64);
        let mut pos = [[1.0, 0.0, 0.0]]; // amplitude A = 1, v = 0
        let mut vel = [[0.0, 0.0, 0.0]];
        let masses = [m];
        let force =
            move |p: &[[f64; 3]]| -> core::result::Result<Vec<[f64; 3]>, ()> {
                Ok(vec![[-k * p[0][0], 0.0, 0.0]])
            };
        let mut forces = force(&pos).unwrap();

        let energy = |pos: &[[f64; 3]], vel: &[[f64; 3]]| {
            0.5 * k * pos[0][0] * pos[0][0] + kinetic_energy(vel, &masses)
        };
        let e0 = energy(&pos, &vel);

        let period = 2.0 * std::f64::consts::PI * (m / k).sqrt();
        let dt = period / 2000.0;
        for _ in 0..2000 {
            // exactly one period
            forces = velocity_verlet_step(&mut pos, &mut vel, &forces, &masses, dt, force).unwrap();
            let e = energy(&pos, &vel);
            assert!((e - e0).abs() / e0 < 1e-3, "energy drift: {e} vs {e0}");
        }
        // After one full period, back to (x≈A, v≈0) — confirms the period.
        assert!((pos[0][0] - 1.0).abs() < 5e-3, "x after one period: {}", pos[0][0]);
        assert!(vel[0][0].abs() < 5e-3, "v after one period: {}", vel[0][0]);
    }

    #[test]
    fn kinetic_energy_is_half_m_v_squared() {
        // One atom, m = 2, v = (3,0,0) → KE = ½·2·9 = 9.
        let ke = kinetic_energy(&[[3.0, 0.0, 0.0]], &[2.0]);
        assert!((ke - 9.0).abs() < 1e-12);
    }
}
