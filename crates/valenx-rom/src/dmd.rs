//! Dynamic Mode Decomposition (DMD), standard and exact.
//!
//! DMD identifies the best-fit *linear* operator `A` such that `x_{k+1} ≈ A x_k`
//! from a snapshot sequence, then returns that operator's spectrum: complex
//! eigenvalues `λ_i` (the discrete-time DMD eigenvalues), the spatial **DMD
//! modes** `φ_i`, and — given the sampling interval `dt` — the continuous-time
//! **growth rates** and **frequencies**.
//!
//! ## Algorithm (Tu et al. 2014, "On Dynamic Mode Decomposition")
//!
//! Split the snapshots into `X = [x_0 … x_{m-2}]` and `X' = [x_1 … x_{m-1}]`.
//! Take the (optionally rank-`r` truncated) SVD `X = U Σ Vᵀ` and form the
//! reduced operator
//!
//! ```text
//! Ã = Uᵣᵀ X' Vᵣ Σᵣ⁻¹     (r × r)
//! ```
//!
//! Eigendecompose `Ã w_i = λ_i w_i`. The DMD eigenvalues are the `λ_i`. Two
//! mode conventions are supported:
//!
//! - **Standard / projected** ([`DmdVariant::Standard`]): `φ_i = Uᵣ w_i`. The
//!   modes live in the column space of `U`.
//! - **Exact** ([`DmdVariant::Exact`]): `φ_i = (1/λ_i) X' Vᵣ Σᵣ⁻¹ w_i`. These
//!   are the true eigenvectors of the full operator `A` (Tu et al.); they
//!   coincide with the standard modes when `X'` lies in the column space of
//!   `U`.
//!
//! ## Continuous time
//!
//! With time step `dt`, the continuous-time eigenvalue is `ω_i = ln(λ_i)/dt`.
//! Then **growth rate** `= Re(ω_i)` (per unit time; negative = decaying) and
//! **angular frequency** `= Im(ω_i)`, so the (ordinary) **frequency** is
//! `Im(ω_i) / (2π)`. A purely real positive `λ_i` has zero frequency; a
//! complex-conjugate pair encodes one oscillation.
//!
//! ## Rank truncation caveat
//!
//! `rank = None` keeps every numerically significant SVD direction of `X`
//! (those above the relative singular-value floor). An explicit `rank` smaller
//! than the data rank yields a *projected, lower-dimensional* operator — the
//! recovered eigenvalues are then approximations, not exact, and over-
//! truncation can drop genuine dynamics. The decaying-oscillator benchmark in
//! the tests is full-rank (`rank = None`) and recovers `λ` to < 1e-6.

use nalgebra::{Complex, DMatrix, DVector};

use crate::error::RomError;
use crate::snapshots::Snapshots;

/// Relative singular-value floor for the truncating SVD of `X`.
const SV_REL_FLOOR: f64 = 1e-12;

/// Which DMD mode convention to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DmdVariant {
    /// Projected modes `φ_i = Uᵣ w_i` (live in the POD subspace).
    Standard,
    /// Exact modes `φ_i = (1/λ_i) X' Vᵣ Σᵣ⁻¹ w_i` (Tu et al. 2014).
    Exact,
}

/// A fitted DMD model: discrete eigenvalues, modes, and continuous-time rates.
#[derive(Debug, Clone)]
pub struct Dmd {
    eigenvalues: Vec<Complex<f64>>,
    modes: Vec<DVector<Complex<f64>>>,
    growth_rates: Vec<f64>,
    frequencies: Vec<f64>,
    dt: f64,
    rank: usize,
    variant: DmdVariant,
}

impl Dmd {
    /// Fit standard (projected) DMD. See [`Dmd::fit_variant`].
    pub fn fit(snapshots: &Snapshots, dt: f64, rank: Option<usize>) -> Result<Self, RomError> {
        Self::fit_variant(snapshots, dt, rank, DmdVariant::Standard)
    }

    /// Fit exact DMD (Tu et al. 2014 exact modes). See [`Dmd::fit_variant`].
    pub fn fit_exact(
        snapshots: &Snapshots,
        dt: f64,
        rank: Option<usize>,
    ) -> Result<Self, RomError> {
        Self::fit_variant(snapshots, dt, rank, DmdVariant::Exact)
    }

    /// Fit DMD with an explicit mode convention.
    ///
    /// `dt` is the (uniform) sampling interval between snapshots; `rank` is the
    /// optional POD truncation rank (`None` = all significant directions).
    ///
    /// # Errors
    /// - [`RomError::BadTimeStep`] if `dt` is not finite and positive.
    /// - [`RomError::NotEnoughSamples`] if there are fewer than 2 snapshots.
    /// - [`RomError::InvalidRank`] if a supplied `rank` is `0` or exceeds the
    ///   numerical rank of `X`.
    /// - [`RomError::RankDeficient`] if `X` is effectively zero.
    /// - [`RomError::NotConverged`] if a factorisation fails.
    pub fn fit_variant(
        snapshots: &Snapshots,
        dt: f64,
        rank: Option<usize>,
        variant: DmdVariant,
    ) -> Result<Self, RomError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(RomError::BadTimeStep { value: dt });
        }
        let m = snapshots.n_time();
        if m < 2 {
            return Err(RomError::NotEnoughSamples {
                what: "DMD",
                needed: 2,
                got: m,
            });
        }
        let full = snapshots.matrix();
        let n = full.nrows();
        // X = first m-1 columns, X' = last m-1 columns.
        let x = full.columns(0, m - 1).into_owned();
        let xp = full.columns(1, m - 1).into_owned();

        // Truncated SVD of X.
        let svd = x.clone().svd(true, true);
        let u = svd.u.ok_or(RomError::NotConverged { what: "DMD SVD" })?;
        let v_t = svd.v_t.ok_or(RomError::NotConverged { what: "DMD SVD" })?;
        let sv = svd.singular_values;

        let smax = sv.iter().cloned().fold(0.0_f64, f64::max);
        let floor = smax * SV_REL_FLOOR;
        if smax <= 0.0 || !smax.is_finite() {
            return Err(RomError::RankDeficient {
                what: "DMD snapshots",
                tol: floor,
            });
        }
        let significant = sv.iter().filter(|&&s| s > floor).count();
        let r = match rank {
            None => significant,
            Some(req) => {
                if req == 0 || req > significant {
                    return Err(RomError::InvalidRank {
                        requested: req,
                        max: significant,
                    });
                }
                req
            }
        };

        // Ur (n x r), Σr⁻¹ (r x r diagonal), Vr (cols of V = rows of Vᵀ → r x (m-1)).
        let ur = u.columns(0, r).into_owned();
        let sigma_inv = DMatrix::from_fn(r, r, |i, j| if i == j { 1.0 / sv[i] } else { 0.0 });
        // Vᵀ is r' x (m-1); take first r rows then transpose → (m-1) x r.
        let vr = v_t.rows(0, r).transpose(); // (m-1) x r

        // Ã = Urᵀ X' Vr Σr⁻¹   (r x r), real.
        let a_tilde = ur.transpose() * &xp * &vr * &sigma_inv;

        // Eigen-decomposition of the small real operator → complex λ_i, w_i.
        let Eig {
            values: eigenvalues,
            vectors: eigvecs,
        } = eig_general(&a_tilde)?;

        // Build the spatial modes per the chosen convention.
        let xp_c = real_to_complex(&xp);
        let vr_c = real_to_complex(&vr);
        let sigma_inv_c = real_to_complex(&sigma_inv);
        let ur_c = real_to_complex(&ur);

        let mut modes: Vec<DVector<Complex<f64>>> = Vec::with_capacity(r);
        for (i, w) in eigvecs.iter().enumerate() {
            let phi = match variant {
                DmdVariant::Standard => &ur_c * w,
                DmdVariant::Exact => {
                    let lam = eigenvalues[i];
                    // φ = (1/λ) X' Vr Σr⁻¹ w
                    let tmp = &xp_c * &vr_c * &sigma_inv_c * w;
                    if lam.norm() <= f64::MIN_POSITIVE {
                        // λ ≈ 0: exact-mode scaling is undefined; fall back to
                        // the unscaled projection rather than dividing by ~0.
                        tmp
                    } else {
                        tmp / lam
                    }
                }
            };
            modes.push(phi);
        }

        // Continuous-time rates: ω = ln(λ)/dt.
        let mut growth_rates = Vec::with_capacity(r);
        let mut frequencies = Vec::with_capacity(r);
        for &lam in &eigenvalues {
            let omega = lam.ln() / dt; // principal branch
            growth_rates.push(omega.re);
            frequencies.push(omega.im / (2.0 * std::f64::consts::PI));
        }

        let _ = n; // state dim retained for clarity; modes already carry it
        Ok(Self {
            eigenvalues,
            modes,
            growth_rates,
            frequencies,
            dt,
            rank: r,
            variant,
        })
    }

    /// The discrete-time DMD eigenvalues `λ_i`.
    pub fn eigenvalues(&self) -> &[Complex<f64>] {
        &self.eigenvalues
    }

    /// The spatial DMD modes `φ_i` (each of length = state dimension).
    pub fn modes(&self) -> &[DVector<Complex<f64>>] {
        &self.modes
    }

    /// Continuous-time growth rates `Re(ln(λ_i)/dt)` (per unit time).
    pub fn growth_rates(&self) -> &[f64] {
        &self.growth_rates
    }

    /// Ordinary frequencies `Im(ln(λ_i)/dt) / (2π)` (cycles per unit time).
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    /// The truncation rank used.
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// The sampling interval `dt`.
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// The mode convention used.
    pub fn variant(&self) -> DmdVariant {
        self.variant
    }
}

/// Promote a real matrix to a complex one (zero imaginary part).
fn real_to_complex(m: &DMatrix<f64>) -> DMatrix<Complex<f64>> {
    DMatrix::from_fn(m.nrows(), m.ncols(), |i, j| Complex::new(m[(i, j)], 0.0))
}

/// Eigen-decomposition of a real, possibly non-symmetric square matrix:
/// complex eigenvalues paired with their unit-norm complex eigenvectors.
struct Eig {
    values: Vec<Complex<f64>>,
    vectors: Vec<DVector<Complex<f64>>>,
}

/// Eigen-decomposition of a real, possibly non-symmetric square matrix into
/// complex eigenvalues and (unit-norm) complex eigenvectors.
///
/// Eigenvalues come from nalgebra's real Schur form (`complex_eigenvalues`).
/// For each eigenvalue `λ` the eigenvector is the smallest right singular
/// vector of `(A − λI)` computed via a complex SVD — a numerically stable
/// null-space extraction that works for the small reduced operators DMD
/// produces.
fn eig_general(a: &DMatrix<f64>) -> Result<Eig, RomError> {
    let r = a.nrows();
    debug_assert_eq!(r, a.ncols());
    let lambdas: Vec<Complex<f64>> = a.complex_eigenvalues().iter().copied().collect();
    if lambdas.len() != r {
        return Err(RomError::NotConverged {
            what: "DMD eigenvalues",
        });
    }
    let a_c = real_to_complex(a);
    let id = DMatrix::<Complex<f64>>::identity(r, r);
    let mut vecs = Vec::with_capacity(r);
    for &lam in &lambdas {
        let shifted = &a_c - &id * lam; // A - λI
        let svd = shifted.svd(false, true);
        let v_t = svd.v_t.ok_or(RomError::NotConverged {
            what: "DMD eigenvector SVD",
        })?;
        // The right singular vector for the smallest singular value is the last
        // row of Vᵀ (singular values are returned in descending order). For a
        // complex SVD, the right singular vectors are the conjugate-transpose
        // rows of Vᵀ → take the row and conjugate it.
        let last = v_t.row(r - 1);
        let vec = DVector::from_fn(r, |i, _| last[i].conj());
        // Normalise (guard the degenerate zero vector).
        let norm = vec.norm();
        let vec = if norm > 0.0 {
            vec / Complex::new(norm, 0.0)
        } else {
            vec
        };
        vecs.push(vec);
    }
    Ok(Eig {
        values: lambdas,
        vectors: vecs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    /// Scalar decaying oscillator: x_{k+1} = λ x_k with a known complex λ.
    /// The 2-D real embedding rotates+scales, so DMD must recover λ's magnitude
    /// and phase, hence the frequency and growth rate, to high precision.
    #[test]
    fn recovers_known_complex_eigenvalue() {
        // Choose continuous-time ω = σ + iΩ, then λ = exp(ω dt).
        let dt = 0.05;
        let sigma = -0.3; // decay
        let omega = 2.0 * std::f64::consts::PI * 1.5; // 1.5 Hz
        let lam = (Complex::new(sigma, omega) * dt).exp();

        // Real 2x2 rotation-scaling block with eigenvalues λ, conj(λ):
        // [[Re λ, -Im λ], [Im λ, Re λ]].
        let a = DMatrix::from_row_slice(2, 2, &[lam.re, -lam.im, lam.im, lam.re]);

        // Generate a trajectory from a real initial condition.
        let mut cols: Vec<Vec<f64>> = Vec::new();
        let mut state = DVector::from_row_slice(&[1.0, 0.0]);
        let steps = 60;
        for _ in 0..steps {
            cols.push(vec![state[0], state[1]]);
            state = &a * &state;
        }
        let snaps = Snapshots::from_columns(&cols).unwrap();

        let dmd = Dmd::fit(&snaps, dt, None).unwrap();
        assert_eq!(dmd.rank(), 2);

        // One of the recovered eigenvalues must match λ in magnitude & phase.
        let best = dmd
            .eigenvalues()
            .iter()
            .map(|e| (e - lam).norm())
            .fold(f64::INFINITY, f64::min);
        assert!(best < 1e-9, "closest eigenvalue error = {best:e}");

        // Frequency and growth rate to < 1e-6.
        let target_freq = omega / (2.0 * std::f64::consts::PI);
        let got_freq = dmd
            .frequencies()
            .iter()
            .copied()
            .map(|f| (f - target_freq).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(got_freq < 1e-6, "frequency error = {got_freq:e}");

        let got_growth = dmd
            .growth_rates()
            .iter()
            .copied()
            .map(|g| (g - sigma).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(got_growth < 1e-6, "growth-rate error = {got_growth:e}");
    }

    #[test]
    fn exact_and_standard_agree_on_full_rank_data() {
        let dt = 0.1;
        let a = DMatrix::from_row_slice(2, 2, &[0.9, -0.2, 0.2, 0.9]);
        let mut cols: Vec<Vec<f64>> = Vec::new();
        let mut state = DVector::from_row_slice(&[1.0, 0.5]);
        for _ in 0..40 {
            cols.push(vec![state[0], state[1]]);
            state = &a * &state;
        }
        let snaps = Snapshots::from_columns(&cols).unwrap();
        let std = Dmd::fit(&snaps, dt, None).unwrap();
        let exact = Dmd::fit_exact(&snaps, dt, None).unwrap();
        // Eigenvalues are identical regardless of mode convention.
        let mut se: Vec<f64> = std.eigenvalues().iter().map(|e| e.norm()).collect();
        let mut ee: Vec<f64> = exact.eigenvalues().iter().map(|e| e.norm()).collect();
        se.sort_by(|a, b| a.partial_cmp(b).unwrap());
        ee.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (a, b) in se.iter().zip(ee.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn purely_real_decay_has_zero_frequency() {
        // x_{k+1} = 0.5 x_k : λ = 0.5 real, frequency 0, growth = ln(0.5)/dt.
        let dt = 0.2;
        let mut cols = Vec::new();
        let mut v = 1.0_f64;
        for _ in 0..20 {
            cols.push(vec![v]);
            v *= 0.5;
        }
        let snaps = Snapshots::from_columns(&cols).unwrap();
        let dmd = Dmd::fit(&snaps, dt, None).unwrap();
        assert_eq!(dmd.rank(), 1);
        assert!((dmd.eigenvalues()[0].re - 0.5).abs() < 1e-9);
        assert!(dmd.eigenvalues()[0].im.abs() < 1e-9);
        assert!(dmd.frequencies()[0].abs() < 1e-9);
        let expected_growth = (0.5_f64).ln() / dt;
        assert!((dmd.growth_rates()[0] - expected_growth).abs() < 1e-9);
    }

    #[test]
    fn rejects_single_snapshot() {
        let snaps = Snapshots::from_columns(&[vec![1.0, 2.0]]).unwrap();
        assert_eq!(
            Dmd::fit(&snaps, 0.1, None).unwrap_err().code(),
            "not_enough_samples"
        );
    }

    #[test]
    fn rejects_bad_dt() {
        let snaps = Snapshots::from_columns(&[vec![1.0], vec![2.0], vec![3.0]]).unwrap();
        assert_eq!(
            Dmd::fit(&snaps, 0.0, None).unwrap_err().code(),
            "bad_time_step"
        );
        assert_eq!(
            Dmd::fit(&snaps, -1.0, None).unwrap_err().code(),
            "bad_time_step"
        );
    }

    #[test]
    fn rejects_overlarge_rank() {
        let dt = 0.1;
        let a = DMatrix::from_row_slice(2, 2, &[0.9, -0.2, 0.2, 0.9]);
        let mut cols = Vec::new();
        let mut state = DVector::from_row_slice(&[1.0, 0.5]);
        for _ in 0..10 {
            cols.push(vec![state[0], state[1]]);
            state = &a * &state;
        }
        let snaps = Snapshots::from_columns(&cols).unwrap();
        assert_eq!(
            Dmd::fit(&snaps, dt, Some(99)).unwrap_err().code(),
            "invalid_rank"
        );
        assert_eq!(
            Dmd::fit(&snaps, dt, Some(0)).unwrap_err().code(),
            "invalid_rank"
        );
    }
}
