//! Feature 6 — AutoDock4-class force-field scoring function.
//!
//! Where the Vina function ([`crate::score::vina`]) is an empirical
//! Gaussian model, AutoDock 4's scoring function is a physics-style
//! force field — a sum of four physically-motivated terms, each with a
//! published global weight:
//!
//! ```text
//! ΔG = w_vdw·Σ vdW(12-6) + w_hb·Σ E(t)·hbond(12-10)
//!      + w_elec·Σ Coulomb_screened + w_sol·Σ desolvation
//!      + w_tors·N_torsions
//! ```
//!
//! - **van der Waals** — a Lennard-Jones 12-6 term. Each AD4 atom
//!   pair has an equilibrium separation `r_eq` and a well depth
//!   `epsilon`; this module combines per-element parameters by the
//!   Lorentz-Berthelot rules (`r_eq` arithmetic mean, `epsilon`
//!   geometric mean).
//! - **hydrogen bonding** — a 12-10 term that replaces the 12-6 for
//!   donor/acceptor pairs, with a directional weight (here a simple
//!   distance-only model — see the v1 note).
//! - **electrostatics** — a Coulomb term with AutoDock's
//!   distance-dependent (sigmoidal) dielectric, scaled by the product
//!   of the two partial charges.
//! - **desolvation** — a volume-and-charge term: burying an atom's
//!   solvent-accessible volume against the partner costs / gains free
//!   energy depending on its solvation parameter.
//!
//! The weights are AutoDock 4.2's published values
//! (Huey, Morris, Olson & Goodsell, J. Comput. Chem. 28 (2007)
//! 1145-1152).
//!
//! ### v1 note
//!
//! The hydrogen-bond term here is distance-only — the directional
//! `E(t)` factor (which scales the bond by the cosine of the
//! donor-H···acceptor angle) needs explicit hydrogen positions and a
//! lone-pair model that the PDBQT atom list does not carry. A
//! distance-only 12-10 is the honest v1; the angular factor is a
//! documented follow-on. The per-element vdW / solvation table is a
//! representative subset of AD4.2's `AD4_parameters.dat`, not the full
//! file.

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::receptor::Receptor;

/// AutoDock 4.2 published global term weights.
pub mod weights {
    /// van der Waals (12-6) term weight.
    pub const VDW: f64 = 0.1662;
    /// Hydrogen-bond (12-10) term weight.
    pub const HBOND: f64 = 0.1209;
    /// Electrostatics term weight.
    pub const ELEC: f64 = 0.1406;
    /// Desolvation term weight.
    pub const DESOLV: f64 = 0.1322;
    /// Per-torsion entropy term weight (kcal/mol per rotatable bond).
    pub const TORS: f64 = 0.2983;
}

/// Per-atom-type van der Waals + solvation parameters. A representative
/// subset of AutoDock 4.2's `AD4_parameters.dat`.
#[derive(Clone, Copy, Debug)]
struct Ad4Params {
    /// Lennard-Jones equilibrium separation `r_eq` (Å).
    r_eq: f64,
    /// Lennard-Jones well depth `epsilon` (kcal/mol).
    epsilon: f64,
    /// Atomic solvation parameter (kcal/mol·Å⁻³).
    solpar: f64,
    /// Atomic fragmental volume (Å³).
    volume: f64,
}

/// Look up the AD4 force-field parameters for an atom type. Values
/// drawn from AutoDock 4.2's parameter file (subset).
fn params(t: Ad4AtomType) -> Ad4Params {
    use Ad4AtomType::*;
    match t {
        C | A => Ad4Params {
            r_eq: 4.00,
            epsilon: 0.150,
            solpar: -0.00143,
            volume: 33.51,
        },
        N | NA | NS => Ad4Params {
            r_eq: 3.50,
            epsilon: 0.160,
            solpar: -0.00162,
            volume: 22.45,
        },
        OA | OS => Ad4Params {
            r_eq: 3.20,
            epsilon: 0.200,
            solpar: -0.00251,
            volume: 17.16,
        },
        S | SA => Ad4Params {
            r_eq: 4.00,
            epsilon: 0.200,
            solpar: -0.00214,
            volume: 33.51,
        },
        HD | H => Ad4Params {
            r_eq: 2.00,
            epsilon: 0.020,
            solpar: 0.000510,
            volume: 0.00,
        },
        P => Ad4Params {
            r_eq: 4.20,
            epsilon: 0.200,
            solpar: -0.00110,
            volume: 38.79,
        },
        F => Ad4Params {
            r_eq: 3.09,
            epsilon: 0.080,
            solpar: -0.00110,
            volume: 15.45,
        },
        Cl => Ad4Params {
            r_eq: 4.09,
            epsilon: 0.276,
            solpar: -0.00110,
            volume: 35.82,
        },
        Br => Ad4Params {
            r_eq: 4.33,
            epsilon: 0.389,
            solpar: -0.00110,
            volume: 42.57,
        },
        I => Ad4Params {
            r_eq: 4.72,
            epsilon: 0.550,
            solpar: -0.00110,
            volume: 55.06,
        },
        Metal => Ad4Params {
            r_eq: 1.98,
            epsilon: 0.875,
            solpar: -0.00110,
            volume: 1.56,
        },
    }
}

/// `true` if the pair is a hydrogen-bonding pair (one donor, one
/// acceptor) — the 12-10 term replaces the 12-6 for these.
fn is_hbond_pair(a: Ad4AtomType, b: Ad4AtomType) -> bool {
    let pa = a.props();
    let pb = b.props();
    (pa.is_donor && pb.is_acceptor) || (pa.is_acceptor && pb.is_donor)
}

/// A term-by-term breakdown of an AutoDock4 score.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ad4Terms {
    /// Weighted Σ of the 12-6 van der Waals term.
    pub vdw: f64,
    /// Weighted Σ of the 12-10 directional hydrogen-bond term.
    pub hbond: f64,
    /// Weighted Σ of the screened-Coulomb electrostatics term.
    pub electrostatic: f64,
    /// Weighted Σ of the desolvation term.
    pub desolvation: f64,
    /// Number of torsions used in the entropy term.
    pub n_torsions: usize,
}

impl Ad4Terms {
    /// The weighted torsional-entropy term, `w_tors · n_torsions`.
    pub fn torsional(&self) -> f64 {
        weights::TORS * self.n_torsions as f64
    }

    /// The total predicted free energy of binding (kcal/mol): the four
    /// interaction terms plus the torsional-entropy term.
    pub fn total(&self) -> f64 {
        self.vdw + self.hbond + self.electrostatic + self.desolvation + self.torsional()
    }

    /// The intermolecular interaction energy only (no torsion term) —
    /// the part that lives on an AutoDock4 affinity grid.
    pub fn intermolecular(&self) -> f64 {
        self.vdw + self.hbond + self.electrostatic + self.desolvation
    }
}

/// AutoDock's distance-dependent (sigmoidal) dielectric. At short
/// range the effective dielectric is ~6; far away it climbs toward the
/// bulk-water value ~78. From Mehler & Solmajer (1991), the form
/// AutoDock 4 adopts.
fn epsilon_r(r: f64) -> f64 {
    const A: f64 = -8.5525;
    const B: f64 = 78.4 - A; // lambda·B in the Mehler-Solmajer form
    const LAMBDA: f64 = 0.003627;
    const K: f64 = 7.7839;
    A + B / (1.0 + K * (-LAMBDA * B * r).exp())
}

/// AutoDock's smoothed desolvation distance factor: a Gaussian of
/// width 3.6 Å — burial only matters at close range.
fn desolv_factor(r: f64) -> f64 {
    const SIGMA: f64 = 3.6;
    (-(r * r) / (2.0 * SIGMA * SIGMA)).exp()
}

/// Score one receptor/ligand atom pair with the four AD4 interaction
/// terms (no weights applied yet — the caller scales). `r` is the
/// centre-to-centre distance; `q_l` / `q_r` the partial charges.
fn pair_terms(
    lt: Ad4AtomType,
    rt: Ad4AtomType,
    q_l: f64,
    q_r: f64,
    r: f64,
) -> (f64, f64, f64, f64) {
    // Guard against a zero distance — overlapping atoms give a huge
    // but finite repulsion via the clamp.
    let r = r.max(0.5);
    let pl = params(lt);
    let pr = params(rt);
    // Lorentz-Berthelot combination.
    let r_eq = 0.5 * (pl.r_eq + pr.r_eq);
    let eps = (pl.epsilon * pr.epsilon).sqrt();

    let (vdw, hbond) = if is_hbond_pair(lt, rt) {
        // 12-10 hydrogen-bond term. Coefficients chosen so the well
        // depth is `eps` at `r = r_eq` (distance-only — see module v1
        // note on the missing angular factor).
        let c12 = 5.0 * eps * r_eq.powi(12);
        let c10 = 6.0 * eps * r_eq.powi(10);
        let hb = c12 / r.powi(12) - c10 / r.powi(10);
        (0.0, hb)
    } else {
        // 12-6 Lennard-Jones. c12 / c6 chosen so the minimum is -eps
        // at r = r_eq.
        let c12 = eps * r_eq.powi(12);
        let c6 = 2.0 * eps * r_eq.powi(6);
        let lj = c12 / r.powi(12) - c6 / r.powi(6);
        (lj, 0.0)
    };

    // Screened Coulomb electrostatics.
    let elec = 332.06 * q_l * q_r / (epsilon_r(r) * r);

    // Desolvation: each atom's (solpar + a charge contribution) times
    // the partner's fragmental volume, smoothed by distance.
    let charge_solpar = 0.01097; // AD4's sigma·|q| coefficient
    let sl = pl.solpar + charge_solpar * q_l.abs();
    let sr = pr.solpar + charge_solpar * q_r.abs();
    let desolv = (sl * pr.volume + sr * pl.volume) * desolv_factor(r);

    (vdw, hbond, elec, desolv)
}

/// Score a posed ligand against a receptor with the AutoDock4-class
/// force field, returning the full term breakdown.
///
/// `ligand_atoms` is `(world position, AD4 type, partial charge)` per
/// ligand atom; `n_torsions` feeds the entropy term. Pairs beyond 8 Å
/// are skipped (AutoDock's non-bonded cutoff).
pub fn score_complex(
    receptor: &Receptor,
    ligand_atoms: &[(Vector3<f64>, Ad4AtomType, f64)],
    n_torsions: usize,
) -> Ad4Terms {
    const CUTOFF_SQ: f64 = 8.0 * 8.0;
    let mut terms = Ad4Terms {
        n_torsions,
        ..Ad4Terms::default()
    };
    for &(lp, lt, q_l) in ligand_atoms {
        for ra in &receptor.atoms {
            let r2 = (lp - ra.position).norm_squared();
            if r2 > CUTOFF_SQ {
                continue;
            }
            let r = r2.sqrt();
            let (vdw, hb, elec, desolv) = pair_terms(lt, ra.ad4_type, q_l, ra.partial_charge, r);
            terms.vdw += weights::VDW * vdw;
            terms.hbond += weights::HBOND * hb;
            terms.electrostatic += weights::ELEC * elec;
            terms.desolvation += weights::DESOLV * desolv;
        }
    }
    terms
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_dock::receptor::ReceptorAtom;

    fn carbon_receptor() -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        }
    }

    #[test]
    fn vdw_term_has_a_negative_well() {
        // A carbon probe swept across distance should dip below zero
        // somewhere near the equilibrium separation (~4 Å for C-C).
        let receptor = carbon_receptor();
        let mut min_vdw = f64::INFINITY;
        let mut d = 2.0;
        while d < 8.0 {
            let ligand = vec![(Vector3::new(d, 0.0, 0.0), Ad4AtomType::C, 0.0)];
            let t = score_complex(&receptor, &ligand, 0);
            min_vdw = min_vdw.min(t.vdw);
            d += 0.1;
        }
        assert!(min_vdw < 0.0, "vdW well never went negative: {min_vdw}");
    }

    #[test]
    fn clash_gives_large_positive_vdw() {
        // Two carbons right on top of each other — strong repulsion.
        let receptor = carbon_receptor();
        let ligand = vec![(Vector3::new(0.6, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let t = score_complex(&receptor, &ligand, 0);
        assert!(
            t.vdw > 0.0,
            "overlapping atoms must repel, got vdw {}",
            t.vdw
        );
    }

    #[test]
    fn opposite_charges_give_favourable_electrostatics() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::OA,
                partial_charge: -0.5,
            }],
        };
        // A +0.5 ligand atom near the -0.5 receptor atom.
        let ligand = vec![(Vector3::new(3.0, 0.0, 0.0), Ad4AtomType::N, 0.5)];
        let t = score_complex(&receptor, &ligand, 0);
        assert!(
            t.electrostatic < 0.0,
            "opposite charges should attract, got {}",
            t.electrostatic
        );
    }

    #[test]
    fn like_charges_give_unfavourable_electrostatics() {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::N,
                partial_charge: 0.5,
            }],
        };
        let ligand = vec![(Vector3::new(3.0, 0.0, 0.0), Ad4AtomType::N, 0.5)];
        let t = score_complex(&receptor, &ligand, 0);
        assert!(t.electrostatic > 0.0, "like charges should repel");
    }

    #[test]
    fn hbond_term_fires_for_donor_acceptor() {
        // Receptor HD donor, ligand OA acceptor — the 12-10 term, not
        // the 12-6, carries the energy.
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::HD,
                partial_charge: 0.2,
            }],
        };
        let ligand = vec![(Vector3::new(1.9, 0.0, 0.0), Ad4AtomType::OA, -0.4)];
        let t = score_complex(&receptor, &ligand, 0);
        assert!(t.hbond != 0.0, "H-bond term should be non-zero");
        assert_eq!(t.vdw, 0.0, "H-bond pairs use the 12-10, not the 12-6");
    }

    #[test]
    fn torsional_term_scales_with_torsions() {
        let receptor = carbon_receptor();
        let ligand = vec![(Vector3::new(4.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let t0 = score_complex(&receptor, &ligand, 0);
        let t5 = score_complex(&receptor, &ligand, 5);
        assert_eq!(t0.torsional(), 0.0);
        assert!((t5.torsional() - weights::TORS * 5.0).abs() < 1e-12);
        // The torsional penalty makes the total less favourable.
        assert!(t5.total() > t0.total());
    }

    #[test]
    fn total_is_sum_of_all_terms() {
        let receptor = carbon_receptor();
        let ligand = vec![(Vector3::new(4.0, 0.0, 0.0), Ad4AtomType::C, 0.1)];
        let t = score_complex(&receptor, &ligand, 3);
        let manual = t.vdw + t.hbond + t.electrostatic + t.desolvation + t.torsional();
        assert!((t.total() - manual).abs() < 1e-12);
        assert!((t.intermolecular() - (t.total() - t.torsional())).abs() < 1e-12);
    }

    #[test]
    fn far_pair_contributes_nothing() {
        let receptor = carbon_receptor();
        let ligand = vec![(Vector3::new(20.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let t = score_complex(&receptor, &ligand, 0);
        assert_eq!(t.intermolecular(), 0.0);
    }

    #[test]
    fn dielectric_increases_with_distance() {
        // The sigmoidal dielectric must be monotone increasing.
        assert!(epsilon_r(2.0) < epsilon_r(10.0));
        assert!(epsilon_r(10.0) < epsilon_r(30.0));
    }
}
