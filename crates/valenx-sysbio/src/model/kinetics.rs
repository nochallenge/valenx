//! Kinetic rate laws — features 4, 5 and 6.
//!
//! A [`RateLaw`] turns the current species amounts (and the model
//! parameters) into a scalar reaction *rate* (flux). Three families
//! are supported as real v1 implementations:
//!
//! - **Mass action** (feature 4): rate = `k · ∏ [reactant]^stoich`.
//!   This is the COPASI / PySB default and the only law that is
//!   thermodynamically exact for elementary reactions.
//! - **Michaelis-Menten** (feature 5): rate = `Vmax · [S] / (Km + [S])`
//!   — the standard single-substrate enzyme law, valid under the
//!   quasi-steady-state approximation.
//! - **Hill** (feature 6): rate = `Vmax · [S]^n / (Kd^n + [S]^n)` for
//!   activation, or its `1 - …` complement for repression — the
//!   cooperative-binding law used throughout the gene-circuit modules.
//!
//! A fourth variant, [`RateLaw::Constant`], models a zeroth-order
//! influx (a constitutive promoter, a buffered source). Every law
//! reads species by **index** into the model's species vector so the
//! ODE / SSA engines never do a string lookup in their inner loop.

use serde::{Deserialize, Serialize};

/// A kinetic rate law evaluated against species amounts.
///
/// Indices (`reactant`, `substrate`, …) are positions in the owning
/// model's species vector. `k_*` constants are the *resolved* numeric
/// values — the [`crate::model::Model`] builder substitutes any
/// parameter references before simulation, so the hot loop is
/// allocation- and lookup-free.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum RateLaw {
    /// Zeroth-order constant flux (e.g. a constitutive source term).
    Constant {
        /// Flux value (amount · time⁻¹).
        rate: f64,
    },

    /// Mass-action kinetics: `k · ∏ amount[i]^order[i]`.
    ///
    /// `reactants` pairs a species index with its kinetic order
    /// (usually the stoichiometric coefficient, but kept separate so a
    /// modified mass-action law can use fractional orders).
    MassAction {
        /// Rate constant `k`.
        k: f64,
        /// `(species index, kinetic order)` pairs.
        reactants: Vec<(usize, f64)>,
    },

    /// Irreversible single-substrate Michaelis-Menten:
    /// `vmax · S / (km + S)`.
    MichaelisMenten {
        /// Maximum velocity `Vmax`.
        vmax: f64,
        /// Michaelis constant `Km`.
        km: f64,
        /// Substrate species index.
        substrate: usize,
    },

    /// Hill kinetics — cooperative binding.
    ///
    /// With `repress = false`: `vmax · S^n / (kd^n + S^n)`.
    /// With `repress = true`:  `vmax · kd^n / (kd^n + S^n)`.
    Hill {
        /// Maximum velocity `Vmax`.
        vmax: f64,
        /// Dissociation constant `Kd` (the half-saturation amount).
        kd: f64,
        /// Hill coefficient `n` (cooperativity).
        n: f64,
        /// Regulator species index.
        regulator: usize,
        /// `true` for repression, `false` for activation.
        repress: bool,
    },
}

impl RateLaw {
    /// Evaluate the rate against the `amounts` vector.
    ///
    /// Out-of-range indices contribute a zero amount (defensive — a
    /// validated [`crate::model::Model`] never produces them). Negative
    /// amounts are clamped to zero so a transient negative excursion in
    /// a stiff ODE step cannot produce a `NaN` from a fractional power.
    pub fn rate(&self, amounts: &[f64]) -> f64 {
        let amt = |i: usize| amounts.get(i).copied().unwrap_or(0.0).max(0.0);
        match self {
            RateLaw::Constant { rate } => *rate,
            RateLaw::MassAction { k, reactants } => {
                let mut v = *k;
                for &(idx, order) in reactants {
                    let a = amt(idx);
                    v *= if order == 1.0 { a } else { a.powf(order) };
                }
                v
            }
            RateLaw::MichaelisMenten {
                vmax,
                km,
                substrate,
            } => {
                let s = amt(*substrate);
                if km + s <= 0.0 {
                    0.0
                } else {
                    vmax * s / (km + s)
                }
            }
            RateLaw::Hill {
                vmax,
                kd,
                n,
                regulator,
                repress,
            } => {
                let s = amt(*regulator);
                let sn = s.powf(*n);
                let kdn = kd.powf(*n);
                let denom = kdn + sn;
                if denom <= 0.0 {
                    return if *repress { *vmax } else { 0.0 };
                }
                if *repress {
                    vmax * kdn / denom
                } else {
                    vmax * sn / denom
                }
            }
        }
    }

    /// Indices of every species this law reads. Used by the
    /// next-reaction method to build its reaction-dependency graph.
    pub fn dependencies(&self) -> Vec<usize> {
        match self {
            RateLaw::Constant { .. } => Vec::new(),
            RateLaw::MassAction { reactants, .. } => reactants.iter().map(|&(i, _)| i).collect(),
            RateLaw::MichaelisMenten { substrate, .. } => vec![*substrate],
            RateLaw::Hill { regulator, .. } => vec![*regulator],
        }
    }

    /// A short human-readable tag (`"mass_action"`, `"hill"`, …) for
    /// the SBML writer and diagnostics.
    pub fn kind_tag(&self) -> &'static str {
        match self {
            RateLaw::Constant { .. } => "constant",
            RateLaw::MassAction { .. } => "mass_action",
            RateLaw::MichaelisMenten { .. } => "michaelis_menten",
            RateLaw::Hill { .. } => "hill",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_action_first_order() {
        let law = RateLaw::MassAction {
            k: 0.5,
            reactants: vec![(0, 1.0)],
        };
        assert!((law.rate(&[4.0]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn mass_action_second_order_bimolecular() {
        // A + B -> C  with two distinct reactants.
        let law = RateLaw::MassAction {
            k: 2.0,
            reactants: vec![(0, 1.0), (1, 1.0)],
        };
        assert!((law.rate(&[3.0, 5.0]) - 30.0).abs() < 1e-12);
    }

    #[test]
    fn michaelis_menten_half_saturation() {
        let law = RateLaw::MichaelisMenten {
            vmax: 10.0,
            km: 2.0,
            substrate: 0,
        };
        // At S == Km the rate is exactly Vmax/2.
        assert!((law.rate(&[2.0]) - 5.0).abs() < 1e-12);
        // Saturating limit approaches Vmax.
        assert!(law.rate(&[1e6]) > 9.99);
    }

    #[test]
    fn michaelis_menten_limits_and_quarter_points() {
        // GROUND TRUTH for v = Vmax·S/(Km + S). The half-saturation
        // point v(Km)=Vmax/2 is checked in the test above; this pins the
        // two defining ASYMPTOTIC limits plus an exact interior value the
        // closed form must reproduce:
        //   • S → 0   ⇒ v → 0   (exact 0 at S=0: Vmax·0/(Km+0))
        //   • S → ∞   ⇒ v → Vmax
        //   • S = 3·Km ⇒ v = Vmax·3Km/(4Km) = ¾·Vmax  (exact)
        let vmax = 10.0_f64;
        let km = 2.0_f64;
        let law = RateLaw::MichaelisMenten {
            vmax,
            km,
            substrate: 0,
        };
        // Zero-substrate limit is algebraically exact ⇒ tol 1e-12.
        assert!(
            law.rate(&[0.0]).abs() < 1e-12,
            "v(0) = {} ≠ 0",
            law.rate(&[0.0])
        );
        // Saturating limit: at S = 1e12, v = Vmax·(1 − Km/(Km+S)) differs
        // from Vmax by Vmax·Km/(Km+S) ≈ 2e-11 — comfortably under 1e-9.
        let v_inf = law.rate(&[1.0e12]);
        assert!(
            (v_inf - vmax).abs() < 1e-9,
            "v(1e12) = {v_inf} ≠ Vmax = {vmax} within 1e-9"
        );
        // Three-quarter-saturation exact value at S = 3·Km.
        let v_three_km = law.rate(&[3.0 * km]);
        assert!(
            (v_three_km - 0.75 * vmax).abs() < 1e-12,
            "v(3·Km) = {v_three_km} ≠ ¾·Vmax = {}",
            0.75 * vmax
        );
    }

    #[test]
    fn hill_activation_and_repression_are_complementary() {
        let act = RateLaw::Hill {
            vmax: 1.0,
            kd: 1.0,
            n: 2.0,
            regulator: 0,
            repress: false,
        };
        let rep = RateLaw::Hill {
            vmax: 1.0,
            kd: 1.0,
            n: 2.0,
            regulator: 0,
            repress: true,
        };
        for &s in &[0.0, 0.5, 1.0, 2.0, 8.0] {
            let sum = act.rate(&[s]) + rep.rate(&[s]);
            assert!((sum - 1.0).abs() < 1e-12, "s={s} sum={sum}");
        }
    }

    #[test]
    fn negative_amount_is_clamped() {
        let law = RateLaw::MassAction {
            k: 1.0,
            reactants: vec![(0, 0.5)],
        };
        // sqrt of a clamped-to-zero amount, not NaN.
        assert_eq!(law.rate(&[-4.0]), 0.0);
    }

    #[test]
    fn dependencies_listed() {
        let law = RateLaw::MassAction {
            k: 1.0,
            reactants: vec![(2, 1.0), (5, 1.0)],
        };
        assert_eq!(law.dependencies(), vec![2, 5]);
    }
}
