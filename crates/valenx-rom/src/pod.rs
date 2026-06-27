//! Proper Orthogonal Decomposition (POD).
//!
//! POD finds the orthonormal spatial basis that captures the most state
//! "energy" (variance, in the 2-norm) for a given number of modes. It is the
//! thin SVD of the snapshot matrix `X` (columns = state-in-time):
//!
//! ```text
//! X = U Σ Vᵀ ,   modes = columns of U,   energy(i) = σ_i²
//! ```
//!
//! The basis is truncated to the smallest `r` whose cumulative singular-value
//! energy `Σ_{i<r} σ_i² / Σ σ_i²` reaches the requested tolerance. Projecting
//! a state onto the basis (`Uᵣᵀ x`) and reconstructing (`Uᵣ a`) are the two
//! workhorse operations; the truncated reconstruction error is the tail energy
//! `sqrt(Σ_{i≥r} σ_i²)` (Eckart–Young).
//!
//! POD here is taken about the **origin** (no mean subtraction): the caller is
//! responsible for centring the data if a fluctuation basis is wanted. This
//! keeps `project`/`reconstruct` exact linear maps and makes the rank-2
//! synthetic-field benchmark in the tests hit machine precision.

use nalgebra::{DMatrix, DVector};

use crate::error::RomError;
use crate::snapshots::Snapshots;

/// Relative singular-value floor below which a mode is treated as numerical
/// noise (multiplied by the largest singular value to get the absolute floor).
const SV_REL_FLOOR: f64 = 1e-12;

/// An energy-truncated POD basis.
#[derive(Debug, Clone)]
pub struct PodBasis {
    /// Orthonormal spatial modes, one per column: `state_dim x rank`.
    modes: DMatrix<f64>,
    /// Retained singular values (length = `rank`), descending.
    singular_values: DVector<f64>,
    /// All singular values of the snapshot matrix (for energy accounting).
    all_singular_values: DVector<f64>,
}

impl PodBasis {
    /// Fit a POD basis to `snapshots`, keeping the fewest modes whose
    /// cumulative energy is at least `energy_tol` (e.g. `0.999` for 99.9 %).
    ///
    /// `energy_tol == 1.0` keeps every numerically significant mode (those
    /// above the relative singular-value floor).
    ///
    /// # Errors
    /// - [`RomError::BadEnergyTol`] if `energy_tol` is not in `(0, 1]`.
    /// - [`RomError::NotConverged`] if the SVD fails.
    /// - [`RomError::RankDeficient`] if every singular value is at/below the
    ///   numerical floor (the snapshot matrix is effectively zero).
    pub fn fit(snapshots: &Snapshots, energy_tol: f64) -> Result<Self, RomError> {
        if !(energy_tol.is_finite()) || energy_tol <= 0.0 || energy_tol > 1.0 {
            return Err(RomError::BadEnergyTol { value: energy_tol });
        }
        let x = snapshots.matrix();
        let svd = x.clone().svd(true, true);
        let u = svd.u.ok_or(RomError::NotConverged { what: "POD SVD" })?;
        let sv = svd.singular_values;

        let smax = sv.iter().cloned().fold(0.0_f64, f64::max);
        let floor = smax * SV_REL_FLOOR;
        if smax <= 0.0 || !smax.is_finite() {
            return Err(RomError::RankDeficient {
                what: "POD snapshots",
                tol: floor,
            });
        }

        // Total energy = Σ σ_i² over numerically significant modes.
        let total_energy: f64 = sv.iter().map(|s| s * s).sum();

        // Grow the basis until cumulative energy reaches the tolerance, but
        // never include a mode at/below the noise floor.
        let mut cum = 0.0_f64;
        let mut rank = 0usize;
        for &s in sv.iter() {
            if s <= floor {
                break;
            }
            cum += s * s;
            rank += 1;
            if cum / total_energy >= energy_tol {
                break;
            }
        }
        // At least one significant mode exists (smax > floor only if smax>0;
        // guard anyway).
        if rank == 0 {
            return Err(RomError::RankDeficient {
                what: "POD snapshots",
                tol: floor,
            });
        }

        let modes = u.columns(0, rank).into_owned();
        let singular_values = DVector::from_iterator(rank, sv.iter().take(rank).copied());
        Ok(Self {
            modes,
            singular_values,
            all_singular_values: sv,
        })
    }

    /// Fit a POD basis truncated to an explicit `rank` (ignores energy).
    ///
    /// # Errors
    /// - [`RomError::InvalidRank`] if `rank` is `0` or exceeds the number of
    ///   numerically significant singular values.
    /// - [`RomError::NotConverged`] if the SVD fails.
    pub fn fit_rank(snapshots: &Snapshots, rank: usize) -> Result<Self, RomError> {
        let x = snapshots.matrix();
        let svd = x.clone().svd(true, true);
        let u = svd.u.ok_or(RomError::NotConverged { what: "POD SVD" })?;
        let sv = svd.singular_values;
        let smax = sv.iter().cloned().fold(0.0_f64, f64::max);
        let floor = smax * SV_REL_FLOOR;
        let significant = sv.iter().filter(|&&s| s > floor).count();
        if rank == 0 || rank > significant {
            return Err(RomError::InvalidRank {
                requested: rank,
                max: significant,
            });
        }
        let modes = u.columns(0, rank).into_owned();
        let singular_values = DVector::from_iterator(rank, sv.iter().take(rank).copied());
        Ok(Self {
            modes,
            singular_values,
            all_singular_values: sv,
        })
    }

    /// Number of retained modes.
    pub fn rank(&self) -> usize {
        self.modes.ncols()
    }

    /// The state (spatial) dimension of each mode.
    pub fn state_dim(&self) -> usize {
        self.modes.nrows()
    }

    /// The orthonormal modes as columns (`state_dim x rank`).
    pub fn modes(&self) -> &DMatrix<f64> {
        &self.modes
    }

    /// The retained singular values (descending, length = [`PodBasis::rank`]).
    pub fn singular_values(&self) -> &DVector<f64> {
        &self.singular_values
    }

    /// Fraction of total snapshot energy captured by the retained modes,
    /// in `[0, 1]`.
    pub fn captured_energy(&self) -> f64 {
        let total: f64 = self.all_singular_values.iter().map(|s| s * s).sum();
        if total <= 0.0 {
            return 0.0;
        }
        let kept: f64 = self.singular_values.iter().map(|s| s * s).sum();
        kept / total
    }

    /// Project a full state `x` onto the reduced coordinates `a = Uᵣᵀ x`.
    ///
    /// # Errors
    /// [`RomError::DimensionMismatch`] if `x.len() != state_dim`.
    pub fn project(&self, x: &DVector<f64>) -> Result<DVector<f64>, RomError> {
        if x.len() != self.state_dim() {
            return Err(RomError::DimensionMismatch {
                what: "POD project input",
                expected: self.state_dim(),
                got: x.len(),
            });
        }
        Ok(self.modes.transpose() * x)
    }

    /// Reconstruct a full state from reduced coordinates `x ≈ Uᵣ a`.
    ///
    /// # Errors
    /// [`RomError::DimensionMismatch`] if `a.len() != rank`.
    pub fn reconstruct(&self, a: &DVector<f64>) -> Result<DVector<f64>, RomError> {
        if a.len() != self.rank() {
            return Err(RomError::DimensionMismatch {
                what: "POD reconstruct coords",
                expected: self.rank(),
                got: a.len(),
            });
        }
        Ok(&self.modes * a)
    }

    /// Project then reconstruct every column of `snapshots` and return the
    /// relative Frobenius reconstruction error `‖X − Uᵣ Uᵣᵀ X‖_F / ‖X‖_F`.
    ///
    /// For a basis fit on the same data this equals the analytic tail energy
    /// `sqrt(Σ_{i≥r} σ_i²) / sqrt(Σ σ_i²)` up to round-off (Eckart–Young).
    ///
    /// # Errors
    /// [`RomError::DimensionMismatch`] if the snapshot state dimension differs
    /// from the basis state dimension.
    pub fn reconstruction_error(&self, snapshots: &Snapshots) -> Result<f64, RomError> {
        let x = snapshots.matrix();
        if x.nrows() != self.state_dim() {
            return Err(RomError::DimensionMismatch {
                what: "reconstruction_error state dim",
                expected: self.state_dim(),
                got: x.nrows(),
            });
        }
        let proj = self.modes.transpose() * x; // rank x time
        let recon = &self.modes * proj; // state x time
        let diff = x - &recon;
        let denom = x.norm();
        if denom == 0.0 {
            return Ok(0.0);
        }
        Ok(diff.norm() / denom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    /// Build a rank-2 synthetic field: two known orthonormal spatial modes with
    /// distinct time coefficients. POD must recover exactly 2 significant modes
    /// and reconstruct to machine precision.
    fn rank2_field() -> (Snapshots, DVector<f64>, DVector<f64>) {
        let n = 8; // state dim
        let m = 20; // time
                    // Two orthonormal spatial modes (cosine modes on a grid).
        let phi1 = DVector::from_fn(n, |i, _| {
            (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos()
        });
        let phi2 = DVector::from_fn(n, |i, _| {
            (2.0 * std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos()
        });
        let phi1 = phi1.normalize();
        let phi2 = phi2.normalize();
        // Distinct, linearly independent time coefficients.
        let mut x = DMatrix::<f64>::zeros(n, m);
        for k in 0..m {
            let t = k as f64 * 0.1;
            let a1 = (t).sin() + 2.0;
            let a2 = (0.7 * t).cos() * 0.5;
            let col = a1 * &phi1 + a2 * &phi2;
            x.set_column(k, &col);
        }
        (Snapshots::from_matrix(x).unwrap(), phi1, phi2)
    }

    #[test]
    fn rank2_field_has_exactly_two_significant_modes() {
        let (snaps, _, _) = rank2_field();
        let basis = PodBasis::fit(&snaps, 1.0).unwrap();
        assert_eq!(basis.rank(), 2, "expected exactly 2 POD modes");
        // The third singular value should be at the noise floor.
        let sv = &basis.all_singular_values;
        let ratio = sv[2] / sv[0];
        assert!(ratio < 1e-12, "third sv ratio = {ratio:e} not negligible");
    }

    #[test]
    fn rank2_reconstruction_is_machine_precision() {
        let (snaps, _, _) = rank2_field();
        let basis = PodBasis::fit(&snaps, 0.9999999999).unwrap();
        let err = basis.reconstruction_error(&snaps).unwrap();
        assert!(err < 1e-10, "reconstruction error = {err:e}");
    }

    #[test]
    fn modes_span_the_known_subspace() {
        // Each true mode, projected then reconstructed, must return unchanged.
        let (snaps, phi1, phi2) = rank2_field();
        let basis = PodBasis::fit(&snaps, 1.0).unwrap();
        for phi in [&phi1, &phi2] {
            let a = basis.project(phi).unwrap();
            let back = basis.reconstruct(&a).unwrap();
            let e = (phi - &back).norm();
            assert!(e < 1e-10, "mode not in POD span: err = {e:e}");
        }
    }

    #[test]
    fn modes_are_orthonormal() {
        let (snaps, _, _) = rank2_field();
        let basis = PodBasis::fit(&snaps, 1.0).unwrap();
        let gram = basis.modes().transpose() * basis.modes();
        let id = DMatrix::<f64>::identity(basis.rank(), basis.rank());
        assert!((gram - id).norm() < 1e-12);
    }

    #[test]
    fn energy_truncation_keeps_fewer_modes_at_lower_tol() {
        let (snaps, _, _) = rank2_field();
        // The first mode dominates; a modest tolerance keeps just one.
        let basis = PodBasis::fit(&snaps, 0.5).unwrap();
        assert_eq!(basis.rank(), 1);
        assert!(basis.captured_energy() >= 0.5);
    }

    #[test]
    fn reconstruction_error_matches_tail_energy() {
        let (snaps, _, _) = rank2_field();
        let basis = PodBasis::fit_rank(&snaps, 1).unwrap();
        let measured = basis.reconstruction_error(&snaps).unwrap();
        // Analytic: sqrt(σ2²) / sqrt(σ1²+σ2²) using all singular values.
        let sv = &basis.all_singular_values;
        let tail = (sv[1] * sv[1]).sqrt();
        let total = (sv.iter().map(|s| s * s).sum::<f64>()).sqrt();
        let analytic = tail / total;
        assert!(
            (measured - analytic).abs() < 1e-10,
            "measured {measured:e} vs analytic {analytic:e}"
        );
    }

    #[test]
    fn rejects_bad_energy_tol() {
        let (snaps, _, _) = rank2_field();
        assert_eq!(
            PodBasis::fit(&snaps, 0.0).unwrap_err().code(),
            "bad_energy_tol"
        );
        assert_eq!(
            PodBasis::fit(&snaps, 1.5).unwrap_err().code(),
            "bad_energy_tol"
        );
    }

    #[test]
    fn rejects_zero_field() {
        let z = DMatrix::<f64>::zeros(4, 6);
        let snaps = Snapshots::from_matrix(z).unwrap();
        assert_eq!(
            PodBasis::fit(&snaps, 0.99).unwrap_err().code(),
            "rank_deficient"
        );
    }

    #[test]
    fn project_dimension_mismatch_errs() {
        let (snaps, _, _) = rank2_field();
        let basis = PodBasis::fit(&snaps, 1.0).unwrap();
        let bad = DVector::<f64>::zeros(basis.state_dim() + 1);
        assert_eq!(
            basis.project(&bad).unwrap_err().code(),
            "dimension_mismatch"
        );
    }
}
