//! Parameter targeting — the shared machinery behind the scan,
//! sensitivity and bifurcation analyses.
//!
//! The analysis modules all need the same primitive: *take a model,
//! change one named knob to a value, hand back the modified model*.
//! Because [`RateLaw`]s store resolved numeric
//! constants (so the simulation hot loop never does a string lookup),
//! a "knob" is addressed structurally by a [`ParamTarget`] — a
//! reaction index plus which constant of that reaction's rate law to
//! overwrite.
//!
//! [`ParamTarget::apply`] returns a fresh model with the one value
//! replaced; the original is untouched, so a scan can fan out over a
//! grid without cloning bugs. A global [`crate::model::Parameter`] can
//! also be targeted by id via [`ParamTarget::Global`].

use crate::error::{Result, SysbioError};
use crate::model::{Model, RateLaw};

/// Which numeric constant of a model a sweep should vary.
#[derive(Debug, Clone, PartialEq)]
pub enum ParamTarget {
    /// The rate constant `k` of reaction `reaction`'s mass-action law.
    MassActionK {
        /// Reaction index.
        reaction: usize,
    },
    /// The `Vmax` of a Michaelis-Menten or Hill law.
    Vmax {
        /// Reaction index.
        reaction: usize,
    },
    /// The `Km` of a Michaelis-Menten law.
    Km {
        /// Reaction index.
        reaction: usize,
    },
    /// The `Kd` of a Hill law.
    Kd {
        /// Reaction index.
        reaction: usize,
    },
    /// The Hill coefficient `n` of a Hill law.
    HillN {
        /// Reaction index.
        reaction: usize,
    },
    /// The flux of a constant-rate law.
    ConstantRate {
        /// Reaction index.
        reaction: usize,
    },
    /// A named global [`crate::model::Parameter`].
    Global {
        /// Parameter id.
        id: String,
    },
    /// The initial amount of a species.
    InitialAmount {
        /// Species index.
        species: usize,
    },
}

impl ParamTarget {
    /// A short human-readable label for axis titles / reports.
    pub fn label(&self) -> String {
        match self {
            ParamTarget::MassActionK { reaction } => format!("k[r{reaction}]"),
            ParamTarget::Vmax { reaction } => format!("Vmax[r{reaction}]"),
            ParamTarget::Km { reaction } => format!("Km[r{reaction}]"),
            ParamTarget::Kd { reaction } => format!("Kd[r{reaction}]"),
            ParamTarget::HillN { reaction } => format!("n[r{reaction}]"),
            ParamTarget::ConstantRate { reaction } => format!("rate[r{reaction}]"),
            ParamTarget::Global { id } => id.clone(),
            ParamTarget::InitialAmount { species } => format!("y0[s{species}]"),
        }
    }

    /// Read the current value of this target from `model`.
    pub fn read(&self, model: &Model) -> Result<f64> {
        let rxn = |i: usize| -> Result<&RateLaw> {
            model
                .reactions
                .get(i)
                .map(|r| &r.rate_law)
                .ok_or_else(|| SysbioError::invalid("reaction", "reaction index out of range"))
        };
        match self {
            ParamTarget::MassActionK { reaction } => match rxn(*reaction)? {
                RateLaw::MassAction { k, .. } => Ok(*k),
                _ => Err(mismatch("mass-action")),
            },
            ParamTarget::Vmax { reaction } => match rxn(*reaction)? {
                RateLaw::MichaelisMenten { vmax, .. } | RateLaw::Hill { vmax, .. } => Ok(*vmax),
                _ => Err(mismatch("MM/Hill")),
            },
            ParamTarget::Km { reaction } => match rxn(*reaction)? {
                RateLaw::MichaelisMenten { km, .. } => Ok(*km),
                _ => Err(mismatch("Michaelis-Menten")),
            },
            ParamTarget::Kd { reaction } => match rxn(*reaction)? {
                RateLaw::Hill { kd, .. } => Ok(*kd),
                _ => Err(mismatch("Hill")),
            },
            ParamTarget::HillN { reaction } => match rxn(*reaction)? {
                RateLaw::Hill { n, .. } => Ok(*n),
                _ => Err(mismatch("Hill")),
            },
            ParamTarget::ConstantRate { reaction } => match rxn(*reaction)? {
                RateLaw::Constant { rate } => Ok(*rate),
                _ => Err(mismatch("constant")),
            },
            ParamTarget::Global { id } => model
                .parameter_index(id)
                .map(|i| model.parameters[i].value)
                .ok_or_else(|| SysbioError::invalid("parameter", "unknown global parameter")),
            ParamTarget::InitialAmount { species } => model
                .species
                .get(*species)
                .map(|s| s.initial)
                .ok_or_else(|| SysbioError::invalid("species", "species index out of range")),
        }
    }

    /// Return a clone of `model` with this target set to `value`.
    pub fn apply(&self, model: &Model, value: f64) -> Result<Model> {
        let mut m = model.clone();
        match self {
            ParamTarget::MassActionK { reaction } => {
                if let RateLaw::MassAction { k, .. } = rate_law_mut(&mut m, *reaction)? {
                    *k = value;
                } else {
                    return Err(mismatch("mass-action"));
                }
            }
            ParamTarget::Vmax { reaction } => match rate_law_mut(&mut m, *reaction)? {
                RateLaw::MichaelisMenten { vmax, .. } | RateLaw::Hill { vmax, .. } => {
                    *vmax = value;
                }
                _ => return Err(mismatch("MM/Hill")),
            },
            ParamTarget::Km { reaction } => {
                if let RateLaw::MichaelisMenten { km, .. } = rate_law_mut(&mut m, *reaction)? {
                    *km = value;
                } else {
                    return Err(mismatch("Michaelis-Menten"));
                }
            }
            ParamTarget::Kd { reaction } => {
                if let RateLaw::Hill { kd, .. } = rate_law_mut(&mut m, *reaction)? {
                    *kd = value;
                } else {
                    return Err(mismatch("Hill"));
                }
            }
            ParamTarget::HillN { reaction } => {
                if let RateLaw::Hill { n, .. } = rate_law_mut(&mut m, *reaction)? {
                    *n = value;
                } else {
                    return Err(mismatch("Hill"));
                }
            }
            ParamTarget::ConstantRate { reaction } => {
                if let RateLaw::Constant { rate } = rate_law_mut(&mut m, *reaction)? {
                    *rate = value;
                } else {
                    return Err(mismatch("constant"));
                }
            }
            ParamTarget::Global { id } => {
                let idx = m
                    .parameter_index(id)
                    .ok_or_else(|| SysbioError::invalid("parameter", "unknown global"))?;
                m.parameters[idx].value = value;
            }
            ParamTarget::InitialAmount { species } => {
                let s = m
                    .species
                    .get_mut(*species)
                    .ok_or_else(|| SysbioError::invalid("species", "species index out of range"))?;
                s.initial = value;
            }
        }
        Ok(m)
    }
}

fn mismatch(expected: &str) -> SysbioError {
    SysbioError::invalid(
        "param_target",
        format!("target does not match a {expected} rate law"),
    )
}

/// Mutable access to reaction `i`'s rate law, range-checked.
fn rate_law_mut(model: &mut Model, i: usize) -> Result<&mut RateLaw> {
    model
        .reactions
        .get_mut(i)
        .map(|r| &mut r.rate_law)
        .ok_or_else(|| SysbioError::invalid("reaction", "reaction index out of range"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Reaction, Species};

    fn mm_model() -> Model {
        let mut m = Model::new("mm");
        let a = m.add_species(Species::new("A", 10.0));
        let b = m.add_species(Species::new("B", 0.0));
        m.add_reaction(Reaction {
            id: "r".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MichaelisMenten {
                vmax: 5.0,
                km: 2.0,
                substrate: a,
            },
            reversible: false,
        });
        m.add_parameter(crate::model::Parameter::new("g", 1.5));
        m
    }

    #[test]
    fn read_and_apply_vmax() {
        let m = mm_model();
        let t = ParamTarget::Vmax { reaction: 0 };
        assert_eq!(t.read(&m).unwrap(), 5.0);
        let m2 = t.apply(&m, 9.0).unwrap();
        assert_eq!(t.read(&m2).unwrap(), 9.0);
        // Original untouched.
        assert_eq!(t.read(&m).unwrap(), 5.0);
    }

    #[test]
    fn read_and_apply_global() {
        let m = mm_model();
        let t = ParamTarget::Global { id: "g".into() };
        assert_eq!(t.read(&m).unwrap(), 1.5);
        let m2 = t.apply(&m, 3.0).unwrap();
        assert_eq!(t.read(&m2).unwrap(), 3.0);
    }

    #[test]
    fn apply_initial_amount() {
        let m = mm_model();
        let t = ParamTarget::InitialAmount { species: 0 };
        let m2 = t.apply(&m, 42.0).unwrap();
        assert_eq!(m2.species[0].initial, 42.0);
    }

    #[test]
    fn type_mismatch_is_an_error() {
        let m = mm_model();
        // Reaction 0 is MM, not mass-action.
        let t = ParamTarget::MassActionK { reaction: 0 };
        assert!(t.read(&m).is_err());
        assert!(t.apply(&m, 1.0).is_err());
    }

    #[test]
    fn out_of_range_reaction_is_an_error() {
        let m = mm_model();
        let t = ParamTarget::Vmax { reaction: 99 };
        assert!(t.read(&m).is_err());
    }
}
