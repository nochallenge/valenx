//! Docking configuration: search box, exhaustiveness, output controls.

use nalgebra::Vector3;

use crate::{DEFAULT_ENERGY_RANGE, DEFAULT_EXHAUSTIVENESS, DEFAULT_NUM_MODES};

/// User-facing knobs for [`crate::dock`].
#[derive(Clone, Debug)]
pub struct DockConfig {
    /// Centre of the search box (Å, receptor coords).
    pub center: Vector3<f64>,
    /// Edge lengths of the search box (Å).
    pub size: Vector3<f64>,
    /// Number of independent ILS chains. Each chain gets its own seed.
    pub exhaustiveness: u32,
    /// Maximum poses to return.
    pub num_modes: u32,
    /// kcal/mol cutoff above best pose for which modes are kept.
    pub energy_range: f64,
    /// Reproducibility seed.
    pub seed: u64,
    /// Grid spacing in Å (Vina default 0.375).
    pub grid_spacing: f64,
}

impl Default for DockConfig {
    fn default() -> Self {
        Self {
            center: Vector3::zeros(),
            size: Vector3::new(20.0, 20.0, 20.0),
            exhaustiveness: DEFAULT_EXHAUSTIVENESS,
            num_modes: DEFAULT_NUM_MODES,
            energy_range: DEFAULT_ENERGY_RANGE,
            seed: 0,
            grid_spacing: 0.375,
        }
    }
}

impl DockConfig {
    /// Grid dimensions derived from box edges + spacing (Vina convention:
    /// at least one cell on each side of centre).
    pub fn grid_dims(&self) -> (usize, usize, usize) {
        let nx = (self.size.x / self.grid_spacing).ceil() as usize + 1;
        let ny = (self.size.y / self.grid_spacing).ceil() as usize + 1;
        let nz = (self.size.z / self.grid_spacing).ceil() as usize + 1;
        (nx.max(2), ny.max(2), nz.max(2))
    }

    /// Grid origin (corner) = centre - size/2.
    pub fn grid_origin(&self) -> Vector3<f64> {
        self.center - self.size / 2.0
    }

    /// Validate that all knobs are physically sensible. Called by
    /// [`crate::dock`], [`crate::runner::dock_with_events`], and
    /// [`crate::dry_run::dock_dry_run`] before any work happens.
    ///
    /// Returns the first violated invariant; callers that need every
    /// problem at once should chain [`Result`]s on the field-by-field
    /// helpers they pass through (`VinaInput::from_case_dir` already
    /// surfaces a similar set of checks at the TOML boundary).
    pub fn validate(&self) -> Result<(), crate::error::DockError> {
        use crate::error::DockError;
        const AXES: [&str; 3] = ["x", "y", "z"];
        // Edge lengths: strictly positive, finite, and bounded above
        // by the largest realistic search volume Vina is asked for
        // in practice (1000 Å is already two orders of magnitude
        // beyond a normal binding-site box).
        for (axis, value) in AXES.iter().zip(self.size.iter()) {
            if !value.is_finite() || *value <= 0.0 {
                return Err(DockError::BadBox {
                    axis,
                    value: *value,
                });
            }
            if *value > 1000.0 {
                return Err(DockError::BadBox {
                    axis,
                    value: *value,
                });
            }
        }
        if !self.grid_spacing.is_finite() || self.grid_spacing <= 0.0 {
            return Err(DockError::BadSpacing(self.grid_spacing));
        }
        if self.exhaustiveness == 0 || self.exhaustiveness > 32 {
            return Err(DockError::BadExhaustiveness(self.exhaustiveness));
        }
        if self.num_modes == 0 {
            return Err(DockError::BadNumModes(self.num_modes));
        }
        if !self.energy_range.is_finite() || self.energy_range <= 0.0 {
            return Err(DockError::BadEnergyRange(self.energy_range));
        }
        // Center components only need to be finite; the box can be
        // placed anywhere in receptor space.
        for (axis, value) in AXES.iter().zip(self.center.iter()) {
            if !value.is_finite() {
                return Err(DockError::BadBox {
                    axis,
                    value: *value,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_has_vina_published_values() {
        let c = DockConfig::default();
        assert_eq!(c.exhaustiveness, 8);
        assert_eq!(c.num_modes, 9);
        assert_eq!(c.energy_range, 3.0);
        assert_eq!(c.grid_spacing, 0.375);
    }

    #[test]
    fn grid_dims_match_box_size() {
        let c = DockConfig {
            size: Vector3::new(7.5, 7.5, 7.5),
            ..DockConfig::default()
        };
        let (nx, ny, nz) = c.grid_dims();
        // 7.5 / 0.375 = 20, +1 = 21
        assert_eq!((nx, ny, nz), (21, 21, 21));
    }

    #[test]
    fn validate_accepts_default() {
        DockConfig::default()
            .validate()
            .expect("default must be valid");
    }

    #[test]
    fn validate_rejects_zero_box_edge() {
        let c = DockConfig {
            size: Vector3::new(20.0, 0.0, 20.0),
            ..DockConfig::default()
        };
        let err = c.validate().unwrap_err();
        assert!(matches!(
            err,
            crate::error::DockError::BadBox { axis: "y", .. }
        ));
    }

    #[test]
    fn validate_rejects_non_finite_box() {
        let c = DockConfig {
            size: Vector3::new(20.0, f64::NAN, 20.0),
            ..DockConfig::default()
        };
        let err = c.validate().unwrap_err();
        assert!(matches!(err, crate::error::DockError::BadBox { .. }));
    }

    #[test]
    fn validate_rejects_oversized_box() {
        let c = DockConfig {
            size: Vector3::new(20.0, 20.0, 1500.0),
            ..DockConfig::default()
        };
        let err = c.validate().unwrap_err();
        assert!(matches!(
            err,
            crate::error::DockError::BadBox { axis: "z", .. }
        ));
    }

    #[test]
    fn validate_rejects_zero_grid_spacing() {
        let c = DockConfig {
            grid_spacing: 0.0,
            ..DockConfig::default()
        };
        assert!(matches!(
            c.validate().unwrap_err(),
            crate::error::DockError::BadSpacing(_)
        ));
    }

    #[test]
    fn validate_rejects_zero_exhaustiveness() {
        let c = DockConfig {
            exhaustiveness: 0,
            ..DockConfig::default()
        };
        assert!(matches!(
            c.validate().unwrap_err(),
            crate::error::DockError::BadExhaustiveness(0)
        ));
    }

    #[test]
    fn validate_rejects_excess_exhaustiveness() {
        let c = DockConfig {
            exhaustiveness: 64,
            ..DockConfig::default()
        };
        assert!(matches!(
            c.validate().unwrap_err(),
            crate::error::DockError::BadExhaustiveness(64)
        ));
    }

    #[test]
    fn validate_rejects_zero_num_modes() {
        let c = DockConfig {
            num_modes: 0,
            ..DockConfig::default()
        };
        assert!(matches!(
            c.validate().unwrap_err(),
            crate::error::DockError::BadNumModes(0)
        ));
    }

    #[test]
    fn validate_rejects_zero_energy_range() {
        let c = DockConfig {
            energy_range: 0.0,
            ..DockConfig::default()
        };
        assert!(matches!(
            c.validate().unwrap_err(),
            crate::error::DockError::BadEnergyRange(_)
        ));
    }

    #[test]
    fn validate_rejects_non_finite_center() {
        let c = DockConfig {
            center: Vector3::new(0.0, f64::INFINITY, 0.0),
            ..DockConfig::default()
        };
        assert!(matches!(
            c.validate().unwrap_err(),
            crate::error::DockError::BadBox { .. }
        ));
    }
}
