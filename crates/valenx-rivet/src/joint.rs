//! The riveted lap / butt joint and its three failure modes.
//!
//! ## Model
//!
//! A row (or rows) of `n` rivets of diameter `d` connect two plates of
//! thickness `t`. The joint is loaded in tension across the width `w` of
//! the critical cross-section. Three independent strengths are computed,
//! each the product of an allowable stress and the area that resists it:
//!
//! - **Shear** of the rivets. Each rivet shank carries the load across
//!   its circular cross-section `A = π d² / 4`. With `s` shear planes
//!   (1 for a lap / single-cover butt joint, 2 for a double-cover butt
//!   joint) the group resists
//!
//!   `P_shear = n · s · (π d² / 4) · τ`.
//!
//! - **Bearing** (crushing) between rivet shank and plate. The contact
//!   pressure acts on the *projected* rectangle `d · t`, so
//!
//!   `P_bearing = n · d · t · σ_b`.
//!
//! - **Tension** of the plate on the **net section** — the gross width
//!   less the material punched out by the holes in the critical row:
//!
//!   `P_tension = (w − n_row · d) · t · σ_t`,
//!
//!   where `n_row` is the number of rivets in that row.
//!
//! The **joint strength** is the smallest of the three — the joint fails
//! in whichever mode reaches its allowable first. The **efficiency** is
//! that strength divided by the strength of the un-drilled solid plate,
//! `P_solid = w · t · σ_t`, and is always strictly below one because the
//! holes can only remove material, never add it.
//!
//! Conventions: all lengths are metres, stresses pascals, forces newtons.

use crate::error::{Result, RivetError};
use crate::material::Allowables;
use core::f64::consts::PI;
use serde::{Deserialize, Serialize};

/// A pattern of identical rivets in a joint.
///
/// The group is described by how many rivets sit in the critical
/// (tension-carrying) row, how many such rows there are, the rivet
/// diameter, and how many shear planes each rivet crosses.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RivetGroup {
    /// Rivet (shank / hole) diameter `d`, metres.
    pub diameter: f64,
    /// Number of rivets in the critical row, `n_row` (≥ 1). This is the
    /// count of holes that interrupt the net tension section.
    pub rivets_per_row: u32,
    /// Number of rivet rows, `n_rows` (≥ 1). The total rivet count is
    /// `rivets_per_row · rows`.
    pub rows: u32,
    /// Shear planes per rivet, `s` (≥ 1): 1 for a lap or single-cover
    /// butt joint, 2 for a double-cover (double-strap) butt joint.
    pub shear_planes: u32,
}

impl RivetGroup {
    /// Build a validated rivet group.
    ///
    /// # Errors
    ///
    /// Returns [`RivetError::NotPositive`] if `diameter` is not finite
    /// and positive, or [`RivetError::ZeroCount`] if any of
    /// `rivets_per_row`, `rows` or `shear_planes` is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_rivet::RivetGroup;
    ///
    /// // Single row of three 20 mm rivets in single shear.
    /// let g = RivetGroup::new(0.020, 3, 1, 1).unwrap();
    /// assert_eq!(g.total_rivets(), 3);
    /// ```
    pub fn new(diameter: f64, rivets_per_row: u32, rows: u32, shear_planes: u32) -> Result<Self> {
        Ok(Self {
            diameter: RivetError::require_positive("diameter", diameter)?,
            rivets_per_row: RivetError::require_count("rivets_per_row", rivets_per_row)?,
            rows: RivetError::require_count("rows", rows)?,
            shear_planes: RivetError::require_count("shear_planes", shear_planes)?,
        })
    }

    /// Total number of rivets in the group, `rivets_per_row · rows`.
    pub fn total_rivets(&self) -> u32 {
        self.rivets_per_row * self.rows
    }

    /// Cross-sectional area of a single rivet shank, `π d² / 4` (m²).
    pub fn rivet_area(&self) -> f64 {
        PI * self.diameter * self.diameter / 4.0
    }
}

/// A plate being joined, described by the cross-section that carries the
/// joint load.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Plate {
    /// Gross width of the plate at the critical row, `w` (metres).
    pub width: f64,
    /// Plate thickness, `t` (metres).
    pub thickness: f64,
}

impl Plate {
    /// Build a validated plate.
    ///
    /// # Errors
    ///
    /// Returns [`RivetError::NotPositive`] if `width` or `thickness` is
    /// not finite and positive.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_rivet::Plate;
    ///
    /// let p = Plate::new(0.150, 0.010).unwrap();
    /// assert!((p.gross_area() - 0.0015).abs() < 1e-9);
    /// ```
    pub fn new(width: f64, thickness: f64) -> Result<Self> {
        Ok(Self {
            width: RivetError::require_positive("width", width)?,
            thickness: RivetError::require_positive("thickness", thickness)?,
        })
    }

    /// Gross cross-sectional area of the un-drilled plate, `w · t` (m²).
    pub fn gross_area(&self) -> f64 {
        self.width * self.thickness
    }

    /// Net cross-sectional area remaining after `holes` rivet holes of
    /// diameter `d` are punched across the row: `(w − holes · d) · t`.
    ///
    /// # Errors
    ///
    /// Returns [`RivetError::NetSectionNonPositive`] if the holes remove
    /// at least the full width, leaving no material to carry tension.
    pub fn net_area(&self, holes: u32, diameter: f64) -> Result<f64> {
        let removed = holes as f64 * diameter;
        let net_width = self.width - removed;
        if net_width > 0.0 {
            Ok(net_width * self.thickness)
        } else {
            Err(RivetError::NetSectionNonPositive {
                width: self.width,
                holes,
                diameter,
                removed,
            })
        }
    }
}

/// Which failure mode governs (reaches its allowable load first).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FailureMode {
    /// The rivets shear off across their shanks.
    Shear,
    /// The plate or rivet crushes in bearing.
    Bearing,
    /// The plate tears across its net section.
    Tension,
}

/// One riveted joint: a [`RivetGroup`] connecting a [`Plate`] under a set
/// of [`Allowables`].
///
/// Both plates in a real lap/butt joint are assumed identical here, so a
/// single [`Plate`] describes the tension member. Build with
/// [`Joint::new`], then call [`Joint::analyze`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Joint {
    /// The rivet pattern.
    pub group: RivetGroup,
    /// The plate cross-section carrying the load.
    pub plate: Plate,
    /// The permissible working stresses.
    pub allow: Allowables,
}

impl Joint {
    /// Assemble a joint from its already-validated parts.
    ///
    /// The components carry their own invariants (positive, finite,
    /// counts ≥ 1), so this constructor cannot itself fail; it exists so
    /// the public surface is uniform and future cross-component checks
    /// have a home.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_rivet::{Allowables, Joint, Plate, RivetGroup};
    ///
    /// let group = RivetGroup::new(0.020, 3, 1, 1).unwrap();
    /// let plate = Plate::new(0.150, 0.010).unwrap();
    /// let allow = Allowables::new(80.0e6, 160.0e6, 100.0e6).unwrap();
    /// let joint = Joint::new(group, plate, allow);
    /// let r = joint.analyze().unwrap();
    /// assert!(r.efficiency < 1.0);
    /// ```
    pub fn new(group: RivetGroup, plate: Plate, allow: Allowables) -> Self {
        Self {
            group,
            plate,
            allow,
        }
    }

    /// Rivet-shear strength of the group:
    /// `P_shear = n · s · (π d² / 4) · τ` (newtons).
    pub fn shear_strength(&self) -> f64 {
        let n = self.group.total_rivets() as f64;
        let s = self.group.shear_planes as f64;
        n * s * self.group.rivet_area() * self.allow.shear
    }

    /// Plate / rivet bearing strength:
    /// `P_bearing = n · d · t · σ_b` (newtons).
    pub fn bearing_strength(&self) -> f64 {
        let n = self.group.total_rivets() as f64;
        n * self.group.diameter * self.plate.thickness * self.allow.bearing
    }

    /// Plate tension strength on the net section:
    /// `P_tension = (w − n_row · d) · t · σ_t` (newtons).
    ///
    /// # Errors
    ///
    /// Returns [`RivetError::NetSectionNonPositive`] if the holes in the
    /// critical row remove the whole width.
    pub fn tension_strength(&self) -> Result<f64> {
        let net = self
            .plate
            .net_area(self.group.rivets_per_row, self.group.diameter)?;
        Ok(net * self.allow.tension)
    }

    /// Strength of the un-drilled solid plate, `P_solid = w · t · σ_t`
    /// (newtons) — the denominator of the joint efficiency.
    pub fn solid_strength(&self) -> f64 {
        self.plate.gross_area() * self.allow.tension
    }

    /// Evaluate all three failure modes and return the governing
    /// strength, the mode, and the joint efficiency.
    ///
    /// # Errors
    ///
    /// Returns [`RivetError::NetSectionNonPositive`] (via
    /// [`Joint::tension_strength`]) if the net section is non-positive.
    pub fn analyze(&self) -> Result<JointResult> {
        let shear = self.shear_strength();
        let bearing = self.bearing_strength();
        let tension = self.tension_strength()?;

        // Governing = the minimum failure load. Compare numerically and
        // pick the mode that owns that minimum; ties resolve in the order
        // shear → bearing → tension, which is deterministic and matches
        // the usual textbook reporting order.
        let mut mode = FailureMode::Shear;
        let mut strength = shear;
        if bearing < strength {
            mode = FailureMode::Bearing;
            strength = bearing;
        }
        if tension < strength {
            mode = FailureMode::Tension;
            strength = tension;
        }

        let solid = self.solid_strength();
        let efficiency = strength / solid;

        Ok(JointResult {
            shear,
            bearing,
            tension,
            solid,
            strength,
            mode,
            efficiency,
        })
    }
}

/// The full result of analysing a [`Joint`].
///
/// Carries each individual failure load, the governing strength and
/// mode, and the efficiency, so a caller can both read the design answer
/// and see how close the other modes were.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JointResult {
    /// Rivet-shear strength, `P_shear` (N).
    pub shear: f64,
    /// Bearing strength, `P_bearing` (N).
    pub bearing: f64,
    /// Net-section tension strength, `P_tension` (N).
    pub tension: f64,
    /// Solid-plate strength, `P_solid` (N).
    pub solid: f64,
    /// The governing (minimum) joint strength, `P = min(...)` (N).
    pub strength: f64,
    /// Which mode governs.
    pub mode: FailureMode,
    /// Joint efficiency `η = P / P_solid`, in `(0, 1)`.
    pub efficiency: f64,
}
