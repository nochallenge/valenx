//! First-order combustion thermochemistry — compute chamber conditions from
//! the propellants instead of assuming them.
//!
//! Given a propellant pair, mixture ratio and chamber pressure, this solves
//! the **chemical equilibrium** of the combustion products (including the
//! dissociation that caps the real flame temperature) and returns the
//! adiabatic flame temperature, the equilibrium mean molar mass and ratio of
//! specific heats, and the characteristic velocity `c*`. Those are exactly
//! the numbers [`crate::propulsion::EngineDesign`] otherwise takes as fixed
//! hand-entered inputs — so an engine can now be predicted from *what it
//! burns*.
//!
//! Method: a Gibbs-energy equilibrium of the products at a trial temperature
//! (minor species expressed through reaction equilibrium constants, the two
//! majors solved from the element + pressure balance by Newton), wrapped in a
//! bisection on temperature that enforces the adiabatic enthalpy balance.
//!
//! Honest scope: **first-order** thermochemistry. Species thermodynamics use
//! standard heats of formation + entropies with a constant (high-temperature
//! averaged) heat capacity, the reactants enter at the 298 K reference (no
//! cryogenic sensible-heat debit), and only the H/O product set is modelled
//! so far (H₂/LOX). NASA-polynomial thermo, a hydrocarbon (C-bearing) product
//! set, and finite-rate kinetics are the documented refinements — this is not
//! NASA CEA, and it is not combustion CFD.

/// Universal gas constant (J/(mol·K)).
const R: f64 = 8.314_462_618;
/// Thermo reference temperature (K).
const T_REF: f64 = 298.15;

/// A combustion-product species: standard heat of formation + entropy at
/// 298 K and a constant (high-temperature-averaged) molar heat capacity.
#[derive(Clone, Copy)]
struct Species {
    /// Molar mass (g/mol).
    m: f64,
    /// Standard enthalpy of formation at 298 K (J/mol).
    dhf: f64,
    /// Standard entropy at 298 K (J/(mol·K)).
    s298: f64,
    /// Constant high-T-averaged heat capacity (J/(mol·K)).
    cp: f64,
}

// H/O product species. ΔHf°/S°298 are standard JANAF/NIST values; Cp is a
// high-temperature average (monatomics are exactly 5R/2).
const H2O: Species = Species {
    m: 18.015,
    dhf: -241_826.0,
    s298: 188.84,
    cp: 49.0,
};
const H2: Species = Species {
    m: 2.016,
    dhf: 0.0,
    s298: 130.68,
    cp: 33.0,
};
const O2: Species = Species {
    m: 31.998,
    dhf: 0.0,
    s298: 205.15,
    cp: 36.0,
};
const OH: Species = Species {
    m: 17.007,
    dhf: 38_990.0,
    s298: 183.74,
    cp: 34.0,
};
const H: Species = Species {
    m: 1.008,
    dhf: 217_999.0,
    s298: 114.72,
    cp: 20.786,
};
const O: Species = Species {
    m: 15.999,
    dhf: 249_180.0,
    s298: 161.06,
    cp: 21.0,
};

impl Species {
    /// Molar enthalpy at `t` (J/mol).
    fn h(self, t: f64) -> f64 {
        self.dhf + self.cp * (t - T_REF)
    }
    /// Molar entropy at `t` (J/(mol·K)).
    fn s(self, t: f64) -> f64 {
        self.s298 + self.cp * (t / T_REF).ln()
    }
    /// Standard molar Gibbs energy at `t` (J/mol).
    fn g(self, t: f64) -> f64 {
        self.h(t) - t * self.s(t)
    }
}

/// Equilibrium constant of a reaction `Σ νᵢ·speciesᵢ` (products positive),
/// with partial pressures referenced to 1 bar.
fn kp(t: f64, terms: &[(f64, Species)]) -> f64 {
    let dg: f64 = terms.iter().map(|&(nu, sp)| nu * sp.g(t)).sum();
    (-dg / (R * t)).exp()
}

/// Equilibrium product partial pressures (bar) of the H/O system at
/// temperature `t` (K), total pressure `p` (bar), and oxygen-to-hydrogen
/// atom ratio `o_over_h`. Order: `[H2O, H2, O2, OH, H, O]`.
fn equilibrium_ho(t: f64, p: f64, o_over_h: f64) -> [f64; 6] {
    let k1 = kp(t, &[(1.0, H2), (0.5, O2), (-1.0, H2O)]); // H2O ⇌ H2 + ½O2
    let k2 = kp(t, &[(1.0, OH), (0.5, H2), (-1.0, H2O)]); // H2O ⇌ OH + ½H2
    let k3 = kp(t, &[(1.0, H), (-0.5, H2)]); // ½H2 ⇌ H
    let k4 = kp(t, &[(1.0, O), (-0.5, O2)]); // ½O2 ⇌ O

    // Minor species expressed through the equilibria, in terms of the two
    // majors p(H2O) and p(H2).
    let parts = |ph2o: f64, ph2: f64| -> [f64; 6] {
        let ph2 = ph2.max(1e-30);
        let ph2o = ph2o.max(1e-30);
        let po2 = (k1 * ph2o / ph2).powi(2);
        let poh = k2 * ph2o / ph2.sqrt();
        let ph = k3 * ph2.sqrt();
        let po = k4 * po2.sqrt();
        [ph2o, ph2, po2, poh, ph, po]
    };
    // Residuals: total pressure and the O/H atom ratio.
    let resid = |x: [f64; 2]| -> [f64; 2] {
        let [ph2o, ph2, po2, poh, ph, po] = parts(x[0], x[1]);
        let sum = ph2o + ph2 + po2 + poh + ph + po;
        let o_at = ph2o + 2.0 * po2 + poh + po;
        let h_at = 2.0 * ph2o + 2.0 * ph2 + poh + ph;
        [sum - p, o_at / h_at - o_over_h]
    };
    // 2-D Newton with a numerical Jacobian, kept positive + bounded.
    let mut x = [0.7 * p, 0.2 * p];
    for _ in 0..80 {
        let f = resid(x);
        if f[0].abs() < 1e-7 * p && f[1].abs() < 1e-10 {
            break;
        }
        let e0 = (x[0].abs() * 1e-6).max(1e-9);
        let e1 = (x[1].abs() * 1e-6).max(1e-9);
        let f0 = resid([x[0] + e0, x[1]]);
        let f1 = resid([x[0], x[1] + e1]);
        let j = [
            [(f0[0] - f[0]) / e0, (f1[0] - f[0]) / e1],
            [(f0[1] - f[1]) / e0, (f1[1] - f[1]) / e1],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() < 1e-30 {
            break;
        }
        let dx0 = (j[1][1] * f[0] - j[0][1] * f[1]) / det;
        let dx1 = (-j[1][0] * f[0] + j[0][0] * f[1]) / det;
        // Damp to keep the step from leaving the positive orthant.
        let damp = 0.7;
        x = [
            (x[0] - damp * dx0).clamp(1e-12 * p, p),
            (x[1] - damp * dx1).clamp(1e-12 * p, p),
        ];
    }
    parts(x[0], x[1])
}

/// Equilibrium chamber conditions from a combustion calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CombustionResult {
    /// Adiabatic flame (chamber) temperature (K).
    pub chamber_temperature: f64,
    /// Equilibrium ratio of specific heats (frozen at the chamber state).
    pub gamma: f64,
    /// Equilibrium mean molar mass of the products (g/mol).
    pub molar_mass: f64,
    /// Characteristic velocity `c*` (m/s).
    pub c_star: f64,
}

const HO_SPECIES: [Species; 6] = [H2O, H2, O2, OH, H, O];

/// Solve the equilibrium H₂/LOX combustion at the given oxidizer/fuel mass
/// ratio and chamber pressure (bar), returning the adiabatic chamber
/// conditions.
///
/// Stoichiometric H₂/LOX is `MR = 8`; flight engines run fuel-rich
/// (`MR ≈ 5–6`) for a lower molar mass and a higher `c*`.
pub fn combust_h2_lox(mixture_ratio: f64, chamber_pressure_bar: f64) -> CombustionResult {
    // Basis: 1 mol H₂ (2.016 g). Oxidizer mass = MR · 2.016 g → mol O₂.
    let n_o2 = mixture_ratio.max(1e-6) * 2.016 / 31.998;
    let h_atoms = 2.0;
    let o_atoms = 2.0 * n_o2;
    let o_over_h = o_atoms / h_atoms;
    let p = chamber_pressure_bar.max(0.1);

    // Product enthalpy (J) at temperature `t`, scaled to the 2-H-atom basis.
    let h_products = |t: f64| -> f64 {
        let pp = equilibrium_ho(t, p, o_over_h);
        let h_at_p = 2.0 * pp[0] + 2.0 * pp[1] + pp[3] + pp[4];
        let scale = h_atoms / h_at_p;
        (0..6).map(|i| pp[i] * scale * HO_SPECIES[i].h(t)).sum()
    };

    // Adiabatic balance: product enthalpy at T_c equals the reactant enthalpy
    // (zero at the 298 K formation reference). Bisect on T_c.
    let (mut lo, mut hi) = (1_000.0_f64, 5_000.0_f64);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if h_products(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let tc = 0.5 * (lo + hi);

    // Equilibrium mixture properties at the converged chamber temperature.
    let pp = equilibrium_ho(tc, p, o_over_h);
    let n_tot: f64 = pp.iter().sum();
    let mass: f64 = (0..6).map(|i| pp[i] * HO_SPECIES[i].m).sum();
    let molar_mass = mass / n_tot;
    let cp_mix: f64 = (0..6).map(|i| pp[i] / n_tot * HO_SPECIES[i].cp).sum();
    let gamma = cp_mix / (cp_mix - R);
    let r_spec = 8_314.462_618 / molar_mass; // J/(kg·K)
    let gam_f = gamma.sqrt() * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let c_star = (r_spec * tc).sqrt() / gam_f;

    CombustionResult {
        chamber_temperature: tc,
        gamma,
        molar_mass,
        c_star,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dissociation_caps_the_flame_temperature() {
        // H₂/LOX at MR 6, 100 bar. With equilibrium dissociation the flame
        // temperature sits near ~3300–3700 K — far below the ~4500 K a
        // no-dissociation calc would give, which is the whole point of an
        // equilibrium solve.
        let r = combust_h2_lox(6.0, 100.0);
        assert!(
            (3_100.0..3_900.0).contains(&r.chamber_temperature),
            "T_c = {} K",
            r.chamber_temperature
        );
        // c* for H₂/LOX is ~2300–2500 m/s (high, from the low molar mass).
        assert!(
            (2_100.0..2_700.0).contains(&r.c_star),
            "c* = {} m/s",
            r.c_star
        );
        // Light products ⇒ low mean molar mass and a low γ.
        assert!(
            (9.0..18.0).contains(&r.molar_mass),
            "M = {} g/mol",
            r.molar_mass
        );
        assert!((1.10..1.30).contains(&r.gamma), "γ = {}", r.gamma);
    }

    #[test]
    fn richer_mixture_lowers_molar_mass() {
        // More fuel ⇒ more leftover light H₂ ⇒ lower mean molar mass.
        let lean = combust_h2_lox(7.5, 100.0);
        let rich = combust_h2_lox(4.5, 100.0);
        assert!(
            rich.molar_mass < lean.molar_mass,
            "rich M {} should be < lean M {}",
            rich.molar_mass,
            lean.molar_mass
        );
    }

    #[test]
    fn higher_pressure_raises_the_flame_temperature() {
        // Higher chamber pressure suppresses dissociation (Le Chatelier),
        // raising the flame temperature.
        let lo = combust_h2_lox(6.0, 30.0);
        let hi = combust_h2_lox(6.0, 200.0);
        assert!(
            hi.chamber_temperature > lo.chamber_temperature,
            "T_c should rise with pressure: {} -> {}",
            lo.chamber_temperature,
            hi.chamber_temperature
        );
    }

    #[test]
    fn is_deterministic_and_finite() {
        let a = combust_h2_lox(6.0, 100.0);
        let b = combust_h2_lox(6.0, 100.0);
        assert_eq!(a, b);
        assert!(a.chamber_temperature.is_finite() && a.c_star.is_finite());
    }
}
