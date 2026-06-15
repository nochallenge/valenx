//! Hydrostatic force and centre of pressure on a submerged flat plate.
//!
//! Consider a flat plate of area `A` fully submerged in a fluid of
//! density `rho`, its plane inclined at angle `theta` to the horizontal
//! free surface (`theta = 90°` is a vertical wall, `theta = 0°` a
//! horizontal plate). Let `h_c` be the vertical depth of the plate's
//! area centroid and `I_xc` the second moment of the plate area about a
//! horizontal axis through its centroid.
//!
//! **Resultant force.** Because gauge pressure varies linearly with
//! depth, the total force on the plate is the centroidal pressure times
//! the area:
//!
//! ```text
//! F = P_c * A = rho * g * h_c * A
//! ```
//!
//! **Centre of pressure.** The resultant does *not* act through the
//! centroid; it acts at the centre of pressure, which lies *below* it.
//! Measured as a vertical depth,
//!
//! ```text
//! h_cp = h_c + (I_xc * sin^2(theta)) / (h_c * A)
//! ```
//!
//! and measured as a slant distance `y = h / sin(theta)` along the
//! plane from the surface line,
//!
//! ```text
//! y_cp = y_c + I_xc / (y_c * A).
//! ```
//!
//! The extra term `I_xc * sin^2(theta) / (h_c * A)` is strictly positive
//! for any real plate, so the centre of pressure is always deeper than
//! the centroid — and it approaches the centroid as the plate is sunk
//! deeper (`h_c → ∞`), where the pressure field becomes effectively
//! uniform.
//!
//! A [`RectangularPlate`] supplies the area and centroidal second moment
//! for the common rectangle in closed form (`I_xc = b * H^3 / 12`); the
//! lower-level [`SubmergedPlate`] takes those two geometric properties
//! directly so any shape can be analysed.

use crate::error::{
    require_finite, require_non_negative, require_positive, FluidStaticsError, Result,
};
use crate::fluid::{Fluid, STANDARD_GRAVITY};
use serde::{Deserialize, Serialize};
use std::f64::consts::FRAC_PI_2;

/// The hydrostatic loading result for a submerged flat plate.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PlateLoad {
    /// Magnitude of the resultant hydrostatic force normal to the plate,
    /// in newtons: `F = rho * g * h_c * A`.
    pub force_n: f64,
    /// Vertical depth of the centre of pressure below the free surface,
    /// in metres. Always `>= centroid_depth_m` for a real plate.
    pub center_of_pressure_depth_m: f64,
    /// Slant distance of the centre of pressure from the surface line,
    /// measured along the plane of the plate, in metres.
    pub center_of_pressure_slant_m: f64,
    /// Vertical depth of the plate's area centroid below the free
    /// surface, in metres (echoed for convenience).
    pub centroid_depth_m: f64,
}

impl PlateLoad {
    /// How far the centre of pressure sits *below* the centroid, in
    /// metres of vertical depth: `h_cp - h_c`. Always `>= 0`; it shrinks
    /// toward zero as the plate is submerged ever deeper.
    pub fn cp_below_centroid_m(&self) -> f64 {
        self.center_of_pressure_depth_m - self.centroid_depth_m
    }
}

/// A flat plate of arbitrary shape submerged in a fluid, described by
/// the two geometric quantities the hydrostatics needs: its area and the
/// second moment of that area about a horizontal centroidal axis.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SubmergedPlate {
    /// Plate area, in square metres. Strictly positive.
    area_m2: f64,
    /// Second moment of area about a horizontal axis through the
    /// centroid (`I_xc`), in metres to the fourth power. Non-negative.
    second_moment_m4: f64,
    /// Inclination of the plate plane to the horizontal free surface, in
    /// radians, in the open interval `(0, pi/2]`. `pi/2` is a vertical
    /// plate.
    inclination_rad: f64,
}

impl SubmergedPlate {
    /// Construct a submerged plate from its area, centroidal second
    /// moment of area, and inclination to the horizontal.
    ///
    /// `inclination_rad` must lie in `(0, pi/2]`: a horizontal plate
    /// (`0`) has no centre-of-pressure offset and is handled by the
    /// pressure module instead, while angles beyond `pi/2` simply mirror
    /// angles below it.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `area_m2` is not strictly positive, `second_moment_m4` is
    /// negative / non-finite, or `inclination_rad` is non-finite, and
    /// [`Geometry`](crate::FluidStaticsError::Geometry)
    /// if the inclination is outside `(0, pi/2]`.
    pub fn new(area_m2: f64, second_moment_m4: f64, inclination_rad: f64) -> Result<Self> {
        let area_m2 = require_positive("area_m2", area_m2)?;
        let second_moment_m4 = require_non_negative("second_moment_m4", second_moment_m4)?;
        let inclination_rad = require_finite("inclination_rad", inclination_rad)?;
        if inclination_rad <= 0.0 || inclination_rad > FRAC_PI_2 {
            return Err(FluidStaticsError::geometry(
                "submerged plate",
                format!(
                    "inclination {inclination_rad} rad must lie in (0, pi/2]; \
                     a horizontal plate has no centre-of-pressure offset"
                ),
            ));
        }
        Ok(SubmergedPlate {
            area_m2,
            second_moment_m4,
            inclination_rad,
        })
    }

    /// Construct a *vertical* submerged plate (`inclination = pi/2`).
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `area_m2` is not strictly positive or `second_moment_m4` is
    /// negative / non-finite.
    pub fn vertical(area_m2: f64, second_moment_m4: f64) -> Result<Self> {
        SubmergedPlate::new(area_m2, second_moment_m4, FRAC_PI_2)
    }

    /// The plate area, in square metres.
    pub fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// The centroidal second moment of area, in metres to the fourth.
    pub fn second_moment_m4(&self) -> f64 {
        self.second_moment_m4
    }

    /// The inclination of the plate to the horizontal, in radians.
    pub fn inclination_rad(&self) -> f64 {
        self.inclination_rad
    }

    /// The resultant hydrostatic force on the plate, in newtons, given a
    /// centroidal depth `centroid_depth_m`: `F = rho * g * h_c * A`.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `gravity` is not strictly positive or `centroid_depth_m` is
    /// negative / non-finite.
    pub fn resultant_force(
        &self,
        fluid: &Fluid,
        gravity: f64,
        centroid_depth_m: f64,
    ) -> Result<f64> {
        let gravity = require_positive("gravity", gravity)?;
        let centroid_depth_m = require_non_negative("centroid_depth_m", centroid_depth_m)?;
        Ok(fluid.density() * gravity * centroid_depth_m * self.area_m2)
    }

    /// The full hydrostatic load — resultant force and centre of pressure
    /// — for the plate whose centroid sits at vertical depth
    /// `centroid_depth_m` below the free surface.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `gravity` is not strictly positive or `centroid_depth_m` is
    /// negative / non-finite, and
    /// [`Singular`](crate::FluidStaticsError::Singular)
    /// if `centroid_depth_m` is zero (the centroid lies on the free
    /// surface, so the resultant force is zero and the centre of pressure
    /// is indeterminate).
    pub fn load(&self, fluid: &Fluid, gravity: f64, centroid_depth_m: f64) -> Result<PlateLoad> {
        let gravity = require_positive("gravity", gravity)?;
        let centroid_depth_m = require_non_negative("centroid_depth_m", centroid_depth_m)?;
        if centroid_depth_m == 0.0 {
            return Err(FluidStaticsError::singular(
                "centre of pressure",
                "centroid lies on the free surface (depth 0); the resultant \
                 force is zero and its line of action is undefined",
            ));
        }

        let force_n = fluid.density() * gravity * centroid_depth_m * self.area_m2;

        // Slant centroidal distance from the surface line: y_c = h_c / sin(theta).
        let sin_t = self.inclination_rad.sin();
        let y_c = centroid_depth_m / sin_t;

        // Centre of pressure along the slope: y_cp = y_c + I_xc / (y_c * A).
        let y_cp = y_c + self.second_moment_m4 / (y_c * self.area_m2);

        // Back to vertical depth: h_cp = y_cp * sin(theta)
        //  = h_c + I_xc * sin^2(theta) / (h_c * A).
        let h_cp = y_cp * sin_t;

        Ok(PlateLoad {
            force_n,
            center_of_pressure_depth_m: h_cp,
            center_of_pressure_slant_m: y_cp,
            centroid_depth_m,
        })
    }
}

/// A rectangular flat plate of width `b` (horizontal, along the surface
/// line) and height `H` (down the slope), submerged in a fluid.
///
/// Supplies the rectangle's area `A = b * H` and centroidal second
/// moment `I_xc = b * H^3 / 12` in closed form, then defers to
/// [`SubmergedPlate`] for the loading.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RectangularPlate {
    /// Width along the (horizontal) surface line, in metres.
    width_m: f64,
    /// Height down the slope of the plate, in metres.
    height_m: f64,
    /// Inclination to the horizontal, in radians, in `(0, pi/2]`.
    inclination_rad: f64,
}

impl RectangularPlate {
    /// Construct an inclined rectangular plate of the given width,
    /// height and inclination to the horizontal.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `width_m` or `height_m` is not strictly positive, and
    /// [`Geometry`](crate::FluidStaticsError::Geometry)
    /// if the inclination is outside `(0, pi/2]`.
    pub fn new(width_m: f64, height_m: f64, inclination_rad: f64) -> Result<Self> {
        let width_m = require_positive("width_m", width_m)?;
        let height_m = require_positive("height_m", height_m)?;
        // Validate the inclination eagerly via SubmergedPlate so both
        // constructors share one rule.
        let _ = SubmergedPlate::new(width_m * height_m, 1.0, inclination_rad)?;
        Ok(RectangularPlate {
            width_m,
            height_m,
            inclination_rad,
        })
    }

    /// Construct a *vertical* rectangular plate (`inclination = pi/2`),
    /// the common dam-wall / tank-wall case.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `width_m` or `height_m` is not strictly positive.
    pub fn vertical(width_m: f64, height_m: f64) -> Result<Self> {
        RectangularPlate::new(width_m, height_m, FRAC_PI_2)
    }

    /// The plate width, in metres.
    pub fn width_m(&self) -> f64 {
        self.width_m
    }

    /// The plate height (down the slope), in metres.
    pub fn height_m(&self) -> f64 {
        self.height_m
    }

    /// The plate area, in square metres: `A = b * H`.
    pub fn area_m2(&self) -> f64 {
        self.width_m * self.height_m
    }

    /// The centroidal second moment of area, in metres to the fourth:
    /// `I_xc = b * H^3 / 12`.
    pub fn second_moment_m4(&self) -> f64 {
        self.width_m * self.height_m.powi(3) / 12.0
    }

    /// The equivalent [`SubmergedPlate`] with this rectangle's area,
    /// centroidal second moment and inclination.
    pub fn as_submerged_plate(&self) -> SubmergedPlate {
        // All inputs were validated in the constructor, so this cannot
        // fail; constructing here keeps the conversion total.
        SubmergedPlate {
            area_m2: self.area_m2(),
            second_moment_m4: self.second_moment_m4(),
            inclination_rad: self.inclination_rad,
        }
    }

    /// The full hydrostatic load on the rectangle whose *top edge* sits
    /// at vertical depth `top_depth_m` below the free surface.
    ///
    /// The centroid of an inclined rectangle lies half its slant height
    /// below the top edge, i.e. at vertical depth
    /// `h_c = top_depth_m + (H/2) * sin(theta)`.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `gravity` is not strictly positive or `top_depth_m` is negative
    /// / non-finite, and
    /// [`Singular`](crate::FluidStaticsError::Singular)
    /// if the resulting centroidal depth is zero.
    pub fn load_from_top_edge(
        &self,
        fluid: &Fluid,
        gravity: f64,
        top_depth_m: f64,
    ) -> Result<PlateLoad> {
        let top_depth_m = require_non_negative("top_depth_m", top_depth_m)?;
        let half_height_depth = 0.5 * self.height_m * self.inclination_rad.sin();
        let centroid_depth_m = top_depth_m + half_height_depth;
        self.as_submerged_plate()
            .load(fluid, gravity, centroid_depth_m)
    }
}

/// Hydrostatic load on a submerged plate under [`STANDARD_GRAVITY`] — a
/// convenience wrapper around [`SubmergedPlate::load`].
///
/// # Errors
///
/// Propagates the errors of [`SubmergedPlate::load`].
pub fn plate_load_standard(
    plate: &SubmergedPlate,
    fluid: &Fluid,
    centroid_depth_m: f64,
) -> Result<PlateLoad> {
    plate.load(fluid, STANDARD_GRAVITY, centroid_depth_m)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rectangle_area_and_second_moment_are_closed_form() {
        let plate = RectangularPlate::vertical(3.0, 2.0).unwrap();
        assert!((plate.area_m2() - 6.0).abs() < EPS);
        // I = b*H^3/12 = 3 * 8 / 12 = 2.0 m^4.
        assert!(
            (plate.second_moment_m4() - 2.0).abs() < EPS,
            "got {}",
            plate.second_moment_m4()
        );
    }

    #[test]
    fn vertical_wall_resultant_force_matches_hand_calc() {
        // Vertical rectangular gate, width 2 m, height 3 m, top at the
        // free surface, in water. Centroid depth = 1.5 m.
        // F = rho*g*h_c*A = 1000 * 9.80665 * 1.5 * 6 = 88259.85 N.
        let plate = RectangularPlate::vertical(2.0, 3.0).unwrap();
        let load = plate
            .load_from_top_edge(&Fluid::water(), STANDARD_GRAVITY, 0.0)
            .unwrap();
        assert!((load.centroid_depth_m - 1.5).abs() < EPS);
        assert!(
            (load.force_n - 88_259.85).abs() < 1e-2,
            "got {}",
            load.force_n
        );
    }

    #[test]
    fn vertical_wall_center_of_pressure_is_two_thirds_depth() {
        // Classic result: for a vertical rectangle with its top edge on
        // the surface, the centre of pressure is at 2/3 of the total
        // depth H. Here H = 3 m -> h_cp = 2.0 m.
        let plate = RectangularPlate::vertical(2.0, 3.0).unwrap();
        let load = plate
            .load_from_top_edge(&Fluid::water(), STANDARD_GRAVITY, 0.0)
            .unwrap();
        assert!(
            (load.center_of_pressure_depth_m - 2.0).abs() < 1e-9,
            "got {}",
            load.center_of_pressure_depth_m
        );
    }

    #[test]
    fn center_of_pressure_is_below_centroid() {
        // For any finite submergence the CP sits strictly below the
        // centroid by I_xc*sin^2/(h_c*A) > 0.
        let plate = RectangularPlate::vertical(2.0, 3.0).unwrap();
        // Top edge 5 m down: centroid at 6.5 m.
        let load = plate
            .load_from_top_edge(&Fluid::water(), STANDARD_GRAVITY, 5.0)
            .unwrap();
        assert!(load.center_of_pressure_depth_m > load.centroid_depth_m);
        // Analytic offset: I/(h_c*A) = 4.5 / (6.5 * 6) = 0.115384... m.
        let expected_offset = plate.second_moment_m4() / (load.centroid_depth_m * plate.area_m2());
        assert!(
            (load.cp_below_centroid_m() - expected_offset).abs() < 1e-9,
            "offset {} expected {expected_offset}",
            load.cp_below_centroid_m()
        );
    }

    #[test]
    fn cp_offset_shrinks_with_depth() {
        // As the plate sinks deeper, the pressure field is more uniform
        // and the CP-to-centroid gap shrinks toward zero.
        let plate = RectangularPlate::vertical(2.0, 2.0).unwrap();
        let shallow = plate
            .load_from_top_edge(&Fluid::water(), STANDARD_GRAVITY, 1.0)
            .unwrap();
        let deep = plate
            .load_from_top_edge(&Fluid::water(), STANDARD_GRAVITY, 50.0)
            .unwrap();
        assert!(
            deep.cp_below_centroid_m() < shallow.cp_below_centroid_m(),
            "deep {} shallow {}",
            deep.cp_below_centroid_m(),
            shallow.cp_below_centroid_m()
        );
        assert!(deep.cp_below_centroid_m() > 0.0);
    }

    #[test]
    fn general_plate_matches_rectangle_specialisation() {
        // Building a SubmergedPlate from the rectangle's A and I_xc must
        // reproduce the RectangularPlate result exactly.
        let rect = RectangularPlate::vertical(2.5, 4.0).unwrap();
        let general = SubmergedPlate::vertical(rect.area_m2(), rect.second_moment_m4()).unwrap();
        let h_c = 6.0;
        let from_rect = rect
            .as_submerged_plate()
            .load(&Fluid::water(), STANDARD_GRAVITY, h_c)
            .unwrap();
        let from_general = general
            .load(&Fluid::water(), STANDARD_GRAVITY, h_c)
            .unwrap();
        assert!((from_rect.force_n - from_general.force_n).abs() < 1e-6);
        assert!(
            (from_rect.center_of_pressure_depth_m - from_general.center_of_pressure_depth_m).abs()
                < 1e-12
        );
    }

    #[test]
    fn inclined_plate_slant_cp_uses_slope_geometry() {
        // 45-degree plate: the vertical-depth CP offset is reduced by
        // sin^2(45) = 0.5 relative to the same plate held vertical.
        let area = 6.0;
        let i_xc = 2.0;
        let h_c = 4.0;
        let theta = std::f64::consts::FRAC_PI_4;
        let inclined = SubmergedPlate::new(area, i_xc, theta).unwrap();
        let vertical = SubmergedPlate::vertical(area, i_xc).unwrap();
        let li = inclined
            .load(&Fluid::water(), STANDARD_GRAVITY, h_c)
            .unwrap();
        let lv = vertical
            .load(&Fluid::water(), STANDARD_GRAVITY, h_c)
            .unwrap();
        // Force depends only on h_c and A, so it is identical.
        assert!((li.force_n - lv.force_n).abs() < 1e-6);
        // Vertical-depth CP offset scales with sin^2(theta) = 0.5 here.
        let ratio = li.cp_below_centroid_m() / lv.cp_below_centroid_m();
        assert!((ratio - 0.5).abs() < 1e-9, "ratio {ratio}");
    }

    #[test]
    fn force_scales_with_fluid_density() {
        // Same gate in seawater vs fresh water -> force ratio = density ratio.
        let plate = RectangularPlate::vertical(2.0, 3.0).unwrap();
        let fresh = plate
            .load_from_top_edge(&Fluid::water(), STANDARD_GRAVITY, 1.0)
            .unwrap();
        let salt = plate
            .load_from_top_edge(&Fluid::seawater(), STANDARD_GRAVITY, 1.0)
            .unwrap();
        let ratio = salt.force_n / fresh.force_n;
        assert!((ratio - 1025.0 / 1000.0).abs() < 1e-9, "ratio {ratio}");
        // CP depth is purely geometric, independent of density.
        assert!((salt.center_of_pressure_depth_m - fresh.center_of_pressure_depth_m).abs() < 1e-12);
    }

    #[test]
    fn zero_centroid_depth_is_singular() {
        let plate = SubmergedPlate::vertical(4.0, 1.0).unwrap();
        let err = plate.load(&Fluid::water(), STANDARD_GRAVITY, 0.0);
        assert!(err.is_err());
        assert_eq!(
            err.unwrap_err().category(),
            "singular",
            "zero-depth centroid must be a singular error"
        );
    }

    #[test]
    fn constructors_reject_bad_geometry() {
        assert!(RectangularPlate::vertical(0.0, 1.0).is_err());
        assert!(RectangularPlate::vertical(1.0, -2.0).is_err());
        // Horizontal plate (inclination 0) is rejected.
        assert!(RectangularPlate::new(1.0, 1.0, 0.0).is_err());
        // Inclination beyond pi/2 is rejected.
        assert!(SubmergedPlate::new(1.0, 1.0, 2.0).is_err());
        assert!(SubmergedPlate::new(-1.0, 1.0, FRAC_PI_2).is_err());
        assert!(SubmergedPlate::new(1.0, -1.0, FRAC_PI_2).is_err());
    }

    #[test]
    fn standard_gravity_wrapper_matches_explicit() {
        let plate = SubmergedPlate::vertical(6.0, 2.0).unwrap();
        let a = plate_load_standard(&plate, &Fluid::water(), 3.0).unwrap();
        let b = plate.load(&Fluid::water(), STANDARD_GRAVITY, 3.0).unwrap();
        assert!((a.force_n - b.force_n).abs() < 1e-9);
        assert!((a.center_of_pressure_depth_m - b.center_of_pressure_depth_m).abs() < 1e-12);
    }
}
