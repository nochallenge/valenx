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
//! Two product systems are modelled: H/O (H₂/LOX, 6 species) and C/H/O
//! (RP-1/LOX and CH₄/LOX, 8 species, adding CO₂/CO). Minor species follow from
//! reaction equilibrium constants and the majors solved from the element +
//! pressure balance by Newton, wrapped in a temperature bisection that
//! enforces the adiabatic enthalpy balance.
//!
//! Honest scope: **first-order** thermochemistry. Species thermodynamics use
//! standard heats of formation + entropies with a constant (high-temperature
//! averaged) heat capacity, the reactants enter at the 298 K reference (no
//! cryogenic sensible-heat debit), and the product set is curated (the major
//! combustion + dissociation species). Accuracy is honest about its limits:
//! the H/O system (H₂/LOX) lands on the validated ~3500 K flame temperature
//! and ~2400 m/s c*, but the C/H/O system (RP-1/LOX) reproduces only the right
//! *trends* — the constant-Cp model overpredicts the absolute kerolox flame
//! temperature by ~10 % versus NASA CEA (~3670 K). NASA-polynomial thermo and
//! finite-rate kinetics are the documented refinements — this is not NASA CEA,
//! and it is not combustion CFD.

/// Universal gas constant (J/(mol·K)).
const R: f64 = 8.314_462_618;
/// Thermo reference temperature (K).
const T_REF: f64 = 298.15;

/// A combustion-product species: standard heat of formation + entropy at
/// 298 K, a constant (high-temperature-averaged) molar heat capacity, and the
/// hydrogen-atom count (for the enthalpy-balance basis).
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
    /// Number of H atoms in the molecule.
    nh: f64,
}

// Product species. ΔHf°/S°298 are standard JANAF/NIST values; Cp is a
// high-temperature average (monatomics are exactly 5R/2).
const CO2: Species = Species {
    m: 44.010,
    dhf: -393_522.0,
    s298: 213.79,
    cp: 62.0,
    nh: 0.0,
};
const CO: Species = Species {
    m: 28.010,
    dhf: -110_527.0,
    s298: 197.66,
    cp: 37.0,
    nh: 0.0,
};
const H2O: Species = Species {
    m: 18.015,
    dhf: -241_826.0,
    s298: 188.84,
    cp: 49.0,
    nh: 2.0,
};
const H2: Species = Species {
    m: 2.016,
    dhf: 0.0,
    s298: 130.68,
    cp: 33.0,
    nh: 2.0,
};
const O2: Species = Species {
    m: 31.998,
    dhf: 0.0,
    s298: 205.15,
    cp: 36.0,
    nh: 0.0,
};
const OH: Species = Species {
    m: 17.007,
    dhf: 38_990.0,
    s298: 183.74,
    cp: 34.0,
    nh: 1.0,
};
const H: Species = Species {
    m: 1.008,
    dhf: 217_999.0,
    s298: 114.72,
    cp: 20.786,
    nh: 1.0,
};
const O: Species = Species {
    m: 15.999,
    dhf: 249_180.0,
    s298: 161.06,
    cp: 21.0,
    nh: 0.0,
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
/// partial pressures referenced to 1 bar.
fn kp(t: f64, terms: &[(f64, Species)]) -> f64 {
    let dg: f64 = terms.iter().map(|&(nu, sp)| nu * sp.g(t)).sum();
    (-dg / (R * t)).exp()
}

/// Solve a 3×3 linear system `a·x = b` by Cramer's rule.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let d = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let det = d(a);
    if det.abs() < 1e-300 {
        return [0.0; 3];
    }
    let mut out = [0.0; 3];
    for (c, slot) in out.iter_mut().enumerate() {
        let mut m = a;
        for r in 0..3 {
            m[r][c] = b[r];
        }
        *slot = d(m) / det;
    }
    out
}

/// Equilibrium partial pressures (bar) of the **H/O** system at temperature
/// `t` (K), total pressure `p` (bar) and oxygen/hydrogen atom ratio.
/// Order: `[H2O, H2, O2, OH, H, O]`.
fn equilibrium_ho(t: f64, p: f64, o_over_h: f64) -> [f64; 6] {
    let k1 = kp(t, &[(1.0, H2), (0.5, O2), (-1.0, H2O)]);
    let k2 = kp(t, &[(1.0, OH), (0.5, H2), (-1.0, H2O)]);
    let k3 = kp(t, &[(1.0, H), (-0.5, H2)]);
    let k4 = kp(t, &[(1.0, O), (-0.5, O2)]);

    let parts = |ph2o: f64, ph2: f64| -> [f64; 6] {
        let ph2 = ph2.max(1e-30);
        let ph2o = ph2o.max(1e-30);
        let po2 = (k1 * ph2o / ph2).powi(2);
        let poh = k2 * ph2o / ph2.sqrt();
        let ph = k3 * ph2.sqrt();
        let po = k4 * po2.sqrt();
        [ph2o, ph2, po2, poh, ph, po]
    };
    let resid = |x: [f64; 2]| -> [f64; 2] {
        let [ph2o, ph2, po2, poh, ph, po] = parts(x[0], x[1]);
        let sum = ph2o + ph2 + po2 + poh + ph + po;
        let o_at = ph2o + 2.0 * po2 + poh + po;
        let h_at = 2.0 * ph2o + 2.0 * ph2 + poh + ph;
        [sum - p, o_at / h_at - o_over_h]
    };
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
        x = [
            (x[0] - 0.7 * dx0).clamp(1e-12 * p, p),
            (x[1] - 0.7 * dx1).clamp(1e-12 * p, p),
        ];
    }
    parts(x[0], x[1])
}

/// Equilibrium partial pressures (bar) of the **C/H/O** system at temperature
/// `t`, total pressure `p`, and carbon/hydrogen + oxygen/hydrogen atom ratios.
/// Order: `[CO2, CO, H2O, H2, O2, OH, H, O]`. The majors are `(H2, H2O, CO)`;
/// O2/CO2/OH/H/O follow from the reaction equilibria.
fn equilibrium_cho(t: f64, p: f64, c_over_h: f64, o_over_h: f64) -> [f64; 8] {
    // Formation equilibria from the majors.
    let k_h2o = kp(t, &[(-1.0, H2), (-0.5, O2), (1.0, H2O)]); // H2+½O2 → H2O
    let k_co2 = kp(t, &[(-1.0, CO), (-0.5, O2), (1.0, CO2)]); // CO+½O2 → CO2
    let k_oh = kp(t, &[(-0.5, H2), (-0.5, O2), (1.0, OH)]); // ½H2+½O2 → OH
    let k_h = kp(t, &[(-0.5, H2), (1.0, H)]); // ½H2 → H
    let k_o = kp(t, &[(-0.5, O2), (1.0, O)]); // ½O2 → O

    let parts = |ph2: f64, ph2o: f64, pco: f64| -> [f64; 8] {
        let ph2 = ph2.max(1e-30);
        let ph2o = ph2o.max(1e-30);
        let pco = pco.max(1e-30);
        let so2 = ph2o / (k_h2o * ph2); // = p(O2)^½
        let po2 = so2 * so2;
        let pco2 = k_co2 * pco * so2;
        let poh = k_oh * ph2.sqrt() * so2;
        let ph = k_h * ph2.sqrt();
        let po = k_o * so2;
        [pco2, pco, ph2o, ph2, po2, poh, ph, po]
    };
    let resid = |x: [f64; 3]| -> [f64; 3] {
        let [pco2, pco, ph2o, ph2, po2, poh, ph, po] = parts(x[0], x[1], x[2]);
        let sum = pco2 + pco + ph2o + ph2 + po2 + poh + ph + po;
        let c_at = pco2 + pco;
        let h_at = 2.0 * ph2o + 2.0 * ph2 + poh + ph;
        let o_at = 2.0 * pco2 + pco + ph2o + 2.0 * po2 + poh + po;
        [sum - p, c_at / h_at - c_over_h, o_at / h_at - o_over_h]
    };
    // 3-D damped Newton with a numerical Jacobian.
    let mut x = [0.2 * p, 0.4 * p, (c_over_h * 0.5 * p).max(1e-6 * p)];
    for _ in 0..150 {
        let f = resid(x);
        if f[0].abs() < 1e-7 * p && f[1].abs() < 1e-9 && f[2].abs() < 1e-9 {
            break;
        }
        let mut j = [[0.0; 3]; 3];
        for k in 0..3 {
            let mut xe = x;
            let e = (x[k].abs() * 1e-6).max(1e-9);
            xe[k] += e;
            let fe = resid(xe);
            for r in 0..3 {
                j[r][k] = (fe[r] - f[r]) / e;
            }
        }
        let dx = solve3(j, f);
        for k in 0..3 {
            x[k] = (x[k] - 0.6 * dx[k]).clamp(1e-14 * p, p);
        }
    }
    parts(x[0], x[1], x[2])
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

/// A modelled propellant combination.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Propellant {
    /// Liquid hydrogen / liquid oxygen.
    H2Lox,
    /// RP-1 (kerosene, ≈ CH₁.₉₅) / liquid oxygen.
    Rp1Lox,
    /// Liquid methane (CH₄) / liquid oxygen — the SpaceX Raptor propellant.
    Ch4Lox,
}

const HO_SPECIES: [Species; 6] = [H2O, H2, O2, OH, H, O];
const CHO_SPECIES: [Species; 8] = [CO2, CO, H2O, H2, O2, OH, H, O];

/// Shared finish: bisect on temperature for the adiabatic enthalpy balance,
/// then read off the equilibrium mixture properties. `equil(t)` returns the
/// product partial pressures (same order as `species`); `h_react` is the
/// reactant enthalpy and `reactant_h` its hydrogen-atom basis.
fn finish_combustion(
    species: &[Species],
    h_react: f64,
    reactant_h: f64,
    equil: impl Fn(f64) -> Vec<f64>,
) -> CombustionResult {
    let h_balance = |t: f64| -> f64 {
        let pp = equil(t);
        let h_at: f64 = species.iter().zip(&pp).map(|(s, &n)| s.nh * n).sum();
        let scale = reactant_h / h_at.max(1e-30);
        species
            .iter()
            .zip(&pp)
            .map(|(s, &n)| n * scale * s.h(t))
            .sum::<f64>()
            - h_react
    };
    let (mut lo, mut hi) = (1_000.0_f64, 5_000.0_f64);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if h_balance(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let tc = 0.5 * (lo + hi);
    let pp = equil(tc);
    let n_tot: f64 = pp.iter().sum();
    let mass: f64 = species.iter().zip(&pp).map(|(s, &n)| s.m * n).sum();
    let molar_mass = mass / n_tot;
    let cp_mix: f64 = species
        .iter()
        .zip(&pp)
        .map(|(s, &n)| n / n_tot * s.cp)
        .sum();
    let gamma = cp_mix / (cp_mix - R);
    let r_spec = 8_314.462_618 / molar_mass;
    let gam_f = gamma.sqrt() * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let c_star = (r_spec * tc).sqrt() / gam_f;
    CombustionResult {
        chamber_temperature: tc,
        gamma,
        molar_mass,
        c_star,
    }
}

/// Solve the equilibrium combustion of `propellant` at the oxidizer/fuel mass
/// ratio and chamber pressure (bar), returning the adiabatic chamber
/// conditions.
pub fn combust(
    propellant: Propellant,
    mixture_ratio: f64,
    chamber_pressure_bar: f64,
) -> CombustionResult {
    let p = chamber_pressure_bar.max(0.1);
    let mr = mixture_ratio.max(1e-6);
    match propellant {
        Propellant::H2Lox => {
            // Fuel H₂ (2.016 g/mol). Per mol fuel: 2 H atoms, MR·2.016 g O₂.
            let n_o2 = mr * 2.016 / 31.998;
            let o_over_h = 2.0 * n_o2 / 2.0;
            finish_combustion(&HO_SPECIES, 0.0, 2.0, move |t| {
                equilibrium_ho(t, p, o_over_h).to_vec()
            })
        }
        Propellant::Rp1Lox => {
            // Fuel RP-1 ≈ CH₁.₉₅ (13.977 g/mol, ΔHf ≈ −24.7 kJ/mol per unit).
            let (c_at, h_at, fuel_m, fuel_dhf) = (1.0, 1.95, 13.977, -24_700.0);
            let n_o2 = mr * fuel_m / 31.998;
            let c_over_h = c_at / h_at;
            let o_over_h = 2.0 * n_o2 / h_at;
            finish_combustion(&CHO_SPECIES, fuel_dhf, h_at, move |t| {
                equilibrium_cho(t, p, c_over_h, o_over_h).to_vec()
            })
        }
        Propellant::Ch4Lox => {
            // Fuel CH₄ (16.043 g/mol, ΔHf ≈ −74.6 kJ/mol). The Raptor runs
            // methalox at MR ≈ 3.6, ~330 bar.
            let (c_at, h_at, fuel_m, fuel_dhf) = (1.0, 4.0, 16.043, -74_600.0);
            let n_o2 = mr * fuel_m / 31.998;
            let c_over_h = c_at / h_at;
            let o_over_h = 2.0 * n_o2 / h_at;
            finish_combustion(&CHO_SPECIES, fuel_dhf, h_at, move |t| {
                equilibrium_cho(t, p, c_over_h, o_over_h).to_vec()
            })
        }
    }
}

/// Equilibrium H₂/LOX combustion (a thin wrapper over [`combust`]).
///
/// Stoichiometric H₂/LOX is `MR = 8`; flight engines run fuel-rich
/// (`MR ≈ 5–6`) for a lower molar mass and a higher `c*`.
pub fn combust_h2_lox(mixture_ratio: f64, chamber_pressure_bar: f64) -> CombustionResult {
    combust(Propellant::H2Lox, mixture_ratio, chamber_pressure_bar)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h2lox_dissociation_caps_the_flame_temperature() {
        // H₂/LOX at MR 6, 100 bar: equilibrium dissociation holds the flame
        // near ~3300–3700 K, far below the ~4500 K a no-dissociation calc
        // gives — the whole point of the equilibrium solve.
        let r = combust_h2_lox(6.0, 100.0);
        assert!(
            (3_100.0..3_900.0).contains(&r.chamber_temperature),
            "T_c = {} K",
            r.chamber_temperature
        );
        assert!(
            (2_100.0..2_700.0).contains(&r.c_star),
            "c* = {} m/s",
            r.c_star
        );
        assert!(
            (9.0..18.0).contains(&r.molar_mass),
            "M = {} g/mol",
            r.molar_mass
        );
        assert!((1.10..1.30).contains(&r.gamma), "γ = {}", r.gamma);
    }

    #[test]
    fn rp1lox_carbon_system_is_first_order() {
        // RP-1/LOX at MR 2.7, 100 bar. The carbon (C/H/O) system is FIRST
        // ORDER and honestly overpredicts: with constant-Cp thermo the
        // absolute flame temperature comes out ~10% above CEA's ~3670 K (the
        // dissociation entropy is mis-weighted near ~4000 K). The mean molar
        // mass (~22–24 g/mol) and c* (~1800–1900 m/s) sit in the kerolox
        // range and the TRENDS are correct (see rp1_is_heavier_and_slower);
        // NASA-polynomial thermo is the documented fix for absolute accuracy.
        let r = combust(Propellant::Rp1Lox, 2.7, 100.0);
        assert!(
            (3_500.0..4_200.0).contains(&r.chamber_temperature),
            "T_c = {} K (first-order, ~10% high vs CEA's ~3670 K)",
            r.chamber_temperature
        );
        assert!(
            (1_600.0..2_050.0).contains(&r.c_star),
            "c* = {} m/s",
            r.c_star
        );
        assert!(
            (18.0..27.0).contains(&r.molar_mass),
            "M = {} g/mol",
            r.molar_mass
        );
    }

    #[test]
    fn rp1_is_heavier_and_slower_than_h2() {
        // Carbon products ⇒ a much higher molar mass and a lower c* than
        // hydrogen — the reason H₂/LOX upper stages exist.
        let h2 = combust(Propellant::H2Lox, 6.0, 100.0);
        let rp1 = combust(Propellant::Rp1Lox, 2.7, 100.0);
        assert!(rp1.molar_mass > h2.molar_mass + 5.0);
        assert!(rp1.c_star < h2.c_star);
    }

    #[test]
    fn methalox_raptor_propellant() {
        // CH₄/LOX at MR 3.6, 330 bar — the Raptor's operating point. The
        // first-order C/H/O model overpredicts the absolute flame temperature
        // (~10 %, like RP-1), but the comparison physics is right: methalox is
        // lighter (more H) than kerolox, so its c* is higher, and it sits below
        // hydrolox.
        let ch4 = combust(Propellant::Ch4Lox, 3.6, 330.0);
        assert!(
            (3_500.0..4_400.0).contains(&ch4.chamber_temperature),
            "T_c = {} K",
            ch4.chamber_temperature
        );
        assert!(
            (1_650.0..2_150.0).contains(&ch4.c_star),
            "c* = {} m/s",
            ch4.c_star
        );
        assert!(
            (16.0..26.0).contains(&ch4.molar_mass),
            "M = {} g/mol",
            ch4.molar_mass
        );
        let rp1 = combust(Propellant::Rp1Lox, 2.7, 330.0);
        let h2 = combust(Propellant::H2Lox, 6.0, 330.0);
        assert!(
            ch4.c_star > rp1.c_star,
            "methalox c* {} > kerolox {}",
            ch4.c_star,
            rp1.c_star
        );
        assert!(
            ch4.c_star < h2.c_star,
            "methalox c* {} < hydrolox {}",
            ch4.c_star,
            h2.c_star
        );
    }

    #[test]
    fn richer_mixture_lowers_molar_mass() {
        let lean = combust_h2_lox(7.5, 100.0);
        let rich = combust_h2_lox(4.5, 100.0);
        assert!(rich.molar_mass < lean.molar_mass);
    }

    #[test]
    fn higher_pressure_raises_the_flame_temperature() {
        let lo = combust_h2_lox(6.0, 30.0);
        let hi = combust_h2_lox(6.0, 200.0);
        assert!(hi.chamber_temperature > lo.chamber_temperature);
    }

    #[test]
    fn is_deterministic_and_finite() {
        let a = combust(Propellant::Rp1Lox, 2.7, 100.0);
        let b = combust(Propellant::Rp1Lox, 2.7, 100.0);
        assert_eq!(a, b);
        assert!(a.chamber_temperature.is_finite() && a.c_star.is_finite());
    }
}
