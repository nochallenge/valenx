//! Nucleotide substitution models and their transition matrices.
//!
//! A continuous-time Markov model of nucleotide substitution is a 4×4
//! instantaneous-rate matrix `Q` whose stationary distribution is the
//! equilibrium base frequencies `π`. The probability of seeing base
//! `j` a time `t` after base `i` is the `(i, j)` entry of
//! `P(t) = e^{Q t}`.
//!
//! ## The models
//!
//! From least to most general — each adds free parameters:
//!
//! - **JC69** — equal base frequencies, one rate. No free parameters.
//! - **K80** — equal frequencies, separate transition / transversion
//!   rates (`κ` = ti/tv ratio).
//! - **F81** — unequal frequencies, one rate.
//! - **HKY85** — unequal frequencies *and* a `κ` ratio. (F81 + K80.)
//! - **GTR** — unequal frequencies and six free exchangeabilities
//!   (`a..f` for AC, AG, AT, CG, CT, GT). The most general
//!   time-reversible model.
//!
//! ## Matrix exponential
//!
//! Every model here is **time-reversible**, so `Q` can be written
//! `Q = Π⁻¹ S` with `S` symmetric and `Π = diag(π)`. The substitution
//! `Q' = Π^{1/2} Q Π^{-1/2}` is then a *real symmetric* matrix with the
//! same eigenvalues as `Q`; `nalgebra`'s symmetric eigendecomposition
//! diagonalises it, and `P(t) = e^{Qt}` follows from
//! `e^{Qt} = Π^{-1/2} U e^{Λt} Uᵀ Π^{1/2}`. This is numerically far
//! more robust than a general (non-symmetric) eigensolver.

use crate::error::{PhyloError, Result};
use nalgebra::{Matrix4, SymmetricEigen, Vector4};

/// A nucleotide substitution model.
///
/// Base order throughout the crate is **A, C, G, T** (indices 0..3).
#[derive(Debug, Clone, PartialEq)]
pub enum SubstModel {
    /// Jukes-Cantor 1969: equal frequencies, one rate.
    Jc69,
    /// Kimura 1980: equal frequencies, transition/transversion ratio.
    K80 {
        /// Transition / transversion rate ratio `κ` (`> 0`).
        kappa: f64,
    },
    /// Felsenstein 1981: unequal frequencies, one rate.
    F81 {
        /// Equilibrium base frequencies `(πA, πC, πG, πT)`, summing
        /// to 1.
        freqs: [f64; 4],
    },
    /// Hasegawa-Kishino-Yano 1985: unequal frequencies + `κ`.
    Hky85 {
        /// Transition / transversion ratio `κ` (`> 0`).
        kappa: f64,
        /// Equilibrium base frequencies, summing to 1.
        freqs: [f64; 4],
    },
    /// General time-reversible: unequal frequencies + six
    /// exchangeabilities.
    Gtr {
        /// Exchangeabilities `[AC, AG, AT, CG, CT, GT]`, all `> 0`.
        rates: [f64; 6],
        /// Equilibrium base frequencies, summing to 1.
        freqs: [f64; 4],
    },
}

impl SubstModel {
    /// The model's equilibrium base frequencies (A, C, G, T).
    pub fn frequencies(&self) -> [f64; 4] {
        match self {
            SubstModel::Jc69 | SubstModel::K80 { .. } => [0.25; 4],
            SubstModel::F81 { freqs }
            | SubstModel::Hky85 { freqs, .. }
            | SubstModel::Gtr { freqs, .. } => *freqs,
        }
    }

    /// Validates the model parameters.
    ///
    /// # Errors
    /// [`PhyloError::Invalid`] if a `κ` / rate is non-positive or the
    /// base frequencies do not sum to 1 (within `1e-6`) or contain a
    /// non-positive entry.
    pub fn validate(&self) -> Result<()> {
        let check_freqs = |f: &[f64; 4]| -> Result<()> {
            if f.iter().any(|&x| x <= 0.0) {
                return Err(PhyloError::invalid(
                    "freqs",
                    "every base frequency must be positive",
                ));
            }
            let sum: f64 = f.iter().sum();
            if (sum - 1.0).abs() > 1e-6 {
                return Err(PhyloError::invalid(
                    "freqs",
                    format!("base frequencies must sum to 1 (got {sum})"),
                ));
            }
            Ok(())
        };
        match self {
            SubstModel::Jc69 => Ok(()),
            SubstModel::K80 { kappa } => {
                if *kappa <= 0.0 {
                    Err(PhyloError::invalid("kappa", "must be positive"))
                } else {
                    Ok(())
                }
            }
            SubstModel::F81 { freqs } => check_freqs(freqs),
            SubstModel::Hky85 { kappa, freqs } => {
                if *kappa <= 0.0 {
                    return Err(PhyloError::invalid("kappa", "must be positive"));
                }
                check_freqs(freqs)
            }
            SubstModel::Gtr { rates, freqs } => {
                if rates.iter().any(|&r| r <= 0.0) {
                    return Err(PhyloError::invalid(
                        "rates",
                        "every exchangeability must be positive",
                    ));
                }
                check_freqs(freqs)
            }
        }
    }

    /// Builds the symmetric exchangeability matrix `S` (`Q = S Π`,
    /// `S` symmetric, zero diagonal placeholder).
    ///
    /// For every reversible model the off-diagonal rate is
    /// `Q[i][j] = S[i][j] · π[j]`; this returns `S`.
    fn exchangeabilities(&self) -> Matrix4<f64> {
        // Index order A=0, C=1, G=2, T=3.
        // Six pair slots: AC, AG, AT, CG, CT, GT.
        let s = match self {
            SubstModel::Jc69 | SubstModel::F81 { .. } => [1.0; 6],
            SubstModel::K80 { kappa } | SubstModel::Hky85 { kappa, .. } => {
                // AG (transition) and CT (transition) get κ.
                [1.0, *kappa, 1.0, 1.0, *kappa, 1.0]
            }
            SubstModel::Gtr { rates, .. } => *rates,
        };
        let mut m = Matrix4::zeros();
        let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for (k, &(i, j)) in pairs.iter().enumerate() {
            m[(i, j)] = s[k];
            m[(j, i)] = s[k];
        }
        m
    }

    /// Builds the instantaneous-rate matrix `Q`, normalised so the
    /// expected number of substitutions per unit branch length is 1
    /// (so a branch length *is* an expected substitution count).
    pub fn rate_matrix(&self) -> Matrix4<f64> {
        let pi = self.frequencies();
        let s = self.exchangeabilities();
        let mut q = Matrix4::zeros();
        // Off-diagonal: Q[i][j] = S[i][j] * pi[j].
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    q[(i, j)] = s[(i, j)] * pi[j];
                }
            }
        }
        // Diagonal: row sums to zero.
        for i in 0..4 {
            let off: f64 = (0..4).filter(|&j| j != i).map(|j| q[(i, j)]).sum();
            q[(i, i)] = -off;
        }
        // Normalise: expected rate -Σ π_i Q_ii = 1.
        let mu: f64 = (0..4).map(|i| -pi[i] * q[(i, i)]).sum();
        if mu > 0.0 {
            q /= mu;
        }
        q
    }

    /// Precomputes the eigendecomposition used to evaluate `P(t)` for
    /// many branch lengths cheaply.
    ///
    /// # Errors
    /// [`PhyloError::Invalid`] if the model parameters are invalid.
    pub fn transition_engine(&self) -> Result<TransitionMatrix> {
        self.validate()?;
        let q = self.rate_matrix();
        let pi = self.frequencies();

        // Symmetrise: Q' = Π^{1/2} Q Π^{-1/2} is symmetric for a
        // reversible model and shares Q's eigenvalues.
        let sqrt_pi = Vector4::new(
            pi[0].sqrt(),
            pi[1].sqrt(),
            pi[2].sqrt(),
            pi[3].sqrt(),
        );
        let inv_sqrt_pi = Vector4::new(
            1.0 / sqrt_pi[0],
            1.0 / sqrt_pi[1],
            1.0 / sqrt_pi[2],
            1.0 / sqrt_pi[3],
        );
        let mut qsym = Matrix4::zeros();
        for i in 0..4 {
            for j in 0..4 {
                qsym[(i, j)] = sqrt_pi[i] * q[(i, j)] * inv_sqrt_pi[j];
            }
        }
        // Force exact symmetry against tiny round-off before the
        // symmetric eigensolver.
        for i in 0..4 {
            for j in (i + 1)..4 {
                let avg = 0.5 * (qsym[(i, j)] + qsym[(j, i)]);
                qsym[(i, j)] = avg;
                qsym[(j, i)] = avg;
            }
        }
        let eig = SymmetricEigen::new(qsym);
        Ok(TransitionMatrix {
            eigenvalues: eig.eigenvalues,
            eigenvectors: eig.eigenvectors,
            sqrt_pi,
            inv_sqrt_pi,
            freqs: pi,
        })
    }
}

/// A precomputed transition-probability engine for one substitution
/// model: stores the eigendecomposition so `P(t)` is a cheap
/// matrix triple-product for any `t`.
#[derive(Debug, Clone)]
pub struct TransitionMatrix {
    /// Eigenvalues `Λ` of the (symmetrised) rate matrix.
    eigenvalues: Vector4<f64>,
    /// Orthonormal eigenvectors `U` of the symmetrised rate matrix.
    eigenvectors: Matrix4<f64>,
    /// `diag(π)^{1/2}`.
    sqrt_pi: Vector4<f64>,
    /// `diag(π)^{-1/2}`.
    inv_sqrt_pi: Vector4<f64>,
    /// Equilibrium base frequencies.
    freqs: [f64; 4],
}

impl TransitionMatrix {
    /// Equilibrium base frequencies of the underlying model.
    pub fn frequencies(&self) -> [f64; 4] {
        self.freqs
    }

    /// Transition-probability matrix `P(t) = e^{Q t}`.
    ///
    /// `P(t)[i][j]` is the probability of ending in base `j` after a
    /// branch of length `t` starting in base `i`. Each row sums to 1.
    /// A negative `t` is clamped to 0.
    pub fn p(&self, t: f64) -> Matrix4<f64> {
        let t = t.max(0.0);
        // e^{Qt} = Π^{-1/2} U e^{Λt} Uᵀ Π^{1/2}.
        let exp_lambda = Vector4::new(
            (self.eigenvalues[0] * t).exp(),
            (self.eigenvalues[1] * t).exp(),
            (self.eigenvalues[2] * t).exp(),
            (self.eigenvalues[3] * t).exp(),
        );
        // U · diag(exp_lambda) · Uᵀ — the symmetric exponential.
        let mut sym = Matrix4::zeros();
        for i in 0..4 {
            for j in 0..4 {
                let mut acc = 0.0;
                for k in 0..4 {
                    acc += self.eigenvectors[(i, k)]
                        * exp_lambda[k]
                        * self.eigenvectors[(j, k)];
                }
                sym[(i, j)] = acc;
            }
        }
        // De-symmetrise: P = Π^{-1/2} sym Π^{1/2}.
        let mut p = Matrix4::zeros();
        for i in 0..4 {
            for j in 0..4 {
                p[(i, j)] = self.inv_sqrt_pi[i] * sym[(i, j)] * self.sqrt_pi[j];
            }
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A representative non-trivial frequency vector.
    fn freqs() -> [f64; 4] {
        [0.30, 0.20, 0.25, 0.25]
    }

    #[test]
    fn validation_catches_bad_parameters() {
        assert!(SubstModel::Jc69.validate().is_ok());
        assert!(SubstModel::K80 { kappa: -1.0 }.validate().is_err());
        assert!(SubstModel::F81 {
            freqs: [0.5, 0.5, 0.5, 0.5]
        }
        .validate()
        .is_err());
        assert!(SubstModel::Hky85 {
            kappa: 2.0,
            freqs: freqs()
        }
        .validate()
        .is_ok());
    }

    #[test]
    fn rate_matrix_rows_sum_to_zero() {
        for model in [
            SubstModel::Jc69,
            SubstModel::K80 { kappa: 3.0 },
            SubstModel::Hky85 {
                kappa: 2.5,
                freqs: freqs(),
            },
            SubstModel::Gtr {
                rates: [1.0, 2.0, 0.5, 0.8, 3.0, 1.2],
                freqs: freqs(),
            },
        ] {
            let q = model.rate_matrix();
            for i in 0..4 {
                let row: f64 = (0..4).map(|j| q[(i, j)]).sum();
                assert!(row.abs() < 1e-10, "row {i} sum = {row}");
            }
        }
    }

    #[test]
    fn p_of_zero_is_identity() {
        let engine = SubstModel::Hky85 {
            kappa: 2.0,
            freqs: freqs(),
        }
        .transition_engine()
        .unwrap();
        let p0 = engine.p(0.0);
        for i in 0..4 {
            for j in 0..4 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((p0[(i, j)] - want).abs() < 1e-9, "P(0)[{i}][{j}]");
            }
        }
    }

    #[test]
    fn p_rows_are_probability_distributions() {
        let engine = SubstModel::Gtr {
            rates: [1.0, 2.5, 0.7, 1.1, 3.2, 0.9],
            freqs: freqs(),
        }
        .transition_engine()
        .unwrap();
        for &t in &[0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = engine.p(t);
            for i in 0..4 {
                let row: f64 = (0..4).map(|j| p[(i, j)]).sum();
                assert!((row - 1.0).abs() < 1e-9, "P({t}) row {i} = {row}");
                for j in 0..4 {
                    assert!(p[(i, j)] >= -1e-12, "negative probability");
                }
            }
        }
    }

    #[test]
    fn p_converges_to_equilibrium() {
        // At a very long branch every row of P(t) -> the equilibrium
        // frequency vector.
        let f = freqs();
        let engine = SubstModel::Hky85 {
            kappa: 2.0,
            freqs: f,
        }
        .transition_engine()
        .unwrap();
        let p = engine.p(100.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!((p[(i, j)] - f[j]).abs() < 1e-4, "not at equilibrium");
            }
        }
    }

    #[test]
    fn jc69_p_has_the_known_closed_form() {
        // JC69: P_ii = 1/4 + 3/4 e^{-4t/3}, P_ij = 1/4 - 1/4 e^{-4t/3}.
        let engine = SubstModel::Jc69.transition_engine().unwrap();
        let t = 0.3;
        let p = engine.p(t);
        let e = (-4.0 / 3.0 * t).exp();
        let diag = 0.25 + 0.75 * e;
        let off = 0.25 - 0.25 * e;
        for i in 0..4 {
            for j in 0..4 {
                let want = if i == j { diag } else { off };
                assert!((p[(i, j)] - want).abs() < 1e-9, "JC P[{i}][{j}]");
            }
        }
    }

    #[test]
    fn detailed_balance_holds() {
        // Reversibility: π_i Q_ij == π_j Q_ji.
        let f = freqs();
        let q = SubstModel::Gtr {
            rates: [1.3, 2.1, 0.6, 0.9, 2.8, 1.0],
            freqs: f,
        }
        .rate_matrix();
        for i in 0..4 {
            for j in 0..4 {
                let lhs = f[i] * q[(i, j)];
                let rhs = f[j] * q[(j, i)];
                assert!((lhs - rhs).abs() < 1e-10, "not reversible at {i},{j}");
            }
        }
    }
}
