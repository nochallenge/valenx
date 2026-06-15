//! Free-air consumption of a reciprocating pneumatic cylinder.
//!
//! ## Model
//!
//! Each powered stroke sweeps a volume equal to the effective piston area
//! times the stroke length. That volume is filled with *compressed* air at
//! the supply pressure; sizing a compressor needs the equivalent volume of
//! *free* (atmospheric) air, which is larger by the compression ratio
//! `r = p_abs / p_atm` (see [`crate::compression`]):
//!
//! ```text
//! V_swept   = A * L                       (compressed volume per stroke)
//! V_free    = V_swept * r                  (free air per stroke)
//! Q_free    = V_free * N_cycles            (free air for N cycles)
//! ```
//!
//! For a **single-acting** cylinder only the extend stroke consumes air
//! (the spring returns it), so one cycle = one swept bore volume. For a
//! **double-acting** cylinder a cycle is an extend *and* a retract stroke,
//! consuming the bore volume plus the annular volume.
//!
//! Multiply by the cycle *rate* instead of a fixed count and you get a
//! volumetric flow demand (free air per unit time) — exactly the quantity
//! a compressor is rated in (e.g. l/min or m^3/min "FAD", free air
//! delivery).
//!
//! ## Honest scope
//!
//! This counts only the air that fills the swept stroke volume. It ignores
//! the dead volume of ports, fittings and hoses (the unswept clearance
//! that must also be pressurised), valve and line leakage, and any air
//! used for pilot signals or blow-off. Real consumption is therefore
//! somewhat higher; treat this as the irreducible stroke demand.

use crate::cylinder::{Cylinder, Stroke};
use crate::error::{PneumaticsError, Result};

/// Whether a cycle counts one powered stroke (single-acting) or both an
/// extend and a retract stroke (double-acting).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Action {
    /// Single-acting: only the extend stroke consumes air.
    SingleActing,
    /// Double-acting: both extend and retract strokes consume air.
    DoubleActing,
}

/// The compressed (swept) volume of a single stroke, `A * L`, in cubic
/// metres, where `A` is the effective area for that [`Stroke`].
///
/// This is the volume *at supply pressure*; multiply by the compression
/// ratio to convert to free air (see [`free_air_per_stroke`]).
///
/// # Errors
///
/// [`PneumaticsError::NonPositive`] if `stroke_length` is not finite and
/// `> 0`.
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::cylinder::{Cylinder, Stroke};
/// use valenx_pneumatics::consumption::swept_volume;
/// let c = Cylinder::single_acting(0.05).unwrap();
/// // A = pi/4 * 0.05^2, L = 0.1 m.
/// let v = swept_volume(&c, 0.1, Stroke::Extend).unwrap();
/// assert!((v - c.bore_area() * 0.1).abs() < 1e-15);
/// ```
pub fn swept_volume(cylinder: &Cylinder, stroke_length: f64, stroke: Stroke) -> Result<f64> {
    let l = PneumaticsError::positive("stroke_length", stroke_length)?;
    Ok(cylinder.effective_area(stroke) * l)
}

/// Free-air volume consumed by a single stroke, `A * L * r`, in cubic
/// metres, where `r` is the compression ratio computed from the *gauge*
/// supply pressure and `atmospheric_pressure`.
///
/// # Errors
///
/// - [`PneumaticsError::NonPositive`] if `stroke_length` or
///   `atmospheric_pressure` is not finite and `> 0`.
/// - [`PneumaticsError::Negative`] if `gauge_pressure` is negative or
///   non-finite.
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::cylinder::{Cylinder, Stroke};
/// use valenx_pneumatics::consumption::free_air_per_stroke;
/// let c = Cylinder::single_acting(0.05).unwrap();
/// // 6 bar gauge over a 100 kPa atmosphere -> r = 7.
/// let v = free_air_per_stroke(&c, 0.1, Stroke::Extend, 600_000.0, 100_000.0).unwrap();
/// assert!((v - 7.0 * c.bore_area() * 0.1).abs() < 1e-12);
/// ```
pub fn free_air_per_stroke(
    cylinder: &Cylinder,
    stroke_length: f64,
    stroke: Stroke,
    gauge_pressure: f64,
    atmospheric_pressure: f64,
) -> Result<f64> {
    let v_swept = swept_volume(cylinder, stroke_length, stroke)?;
    let r = crate::compression::compression_ratio(gauge_pressure, atmospheric_pressure)?;
    Ok(v_swept * r)
}

/// The combined compressed swept volume of *one full cycle*, in cubic
/// metres: the extend stroke for a single-acting cylinder, or the extend
/// plus retract strokes for a double-acting one.
///
/// # Errors
///
/// [`PneumaticsError::NonPositive`] if `stroke_length` is not finite and
/// `> 0`.
pub fn swept_volume_per_cycle(
    cylinder: &Cylinder,
    stroke_length: f64,
    action: Action,
) -> Result<f64> {
    let extend = swept_volume(cylinder, stroke_length, Stroke::Extend)?;
    match action {
        Action::SingleActing => Ok(extend),
        Action::DoubleActing => {
            let retract = swept_volume(cylinder, stroke_length, Stroke::Retract)?;
            Ok(extend + retract)
        }
    }
}

/// Total free-air consumption over `cycles` full cycles, in cubic metres.
///
/// This is the headline sizing number: `Q = (per-cycle swept volume) * r *
/// cycles`. It scales linearly with both the stroke length (through the
/// swept volume) and the cycle count.
///
/// # Errors
///
/// - [`PneumaticsError::NonPositive`] if `stroke_length` or
///   `atmospheric_pressure` is not finite and `> 0`.
/// - [`PneumaticsError::Negative`] if `gauge_pressure` or `cycles` is
///   negative or non-finite (zero cycles is allowed and yields zero).
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::cylinder::Cylinder;
/// use valenx_pneumatics::consumption::{free_air_consumption, Action};
/// let c = Cylinder::single_acting(0.05).unwrap();
/// // 100 cycles, 0.1 m stroke, 6 bar gauge over 100 kPa atm (r = 7).
/// let q = free_air_consumption(&c, 0.1, Action::SingleActing, 100.0, 600_000.0, 100_000.0).unwrap();
/// let per_cycle = 7.0 * c.bore_area() * 0.1;
/// assert!((q - 100.0 * per_cycle).abs() < 1e-9);
/// ```
pub fn free_air_consumption(
    cylinder: &Cylinder,
    stroke_length: f64,
    action: Action,
    cycles: f64,
    gauge_pressure: f64,
    atmospheric_pressure: f64,
) -> Result<f64> {
    let n = PneumaticsError::non_negative("cycles", cycles)?;
    let per_cycle = swept_volume_per_cycle(cylinder, stroke_length, action)?;
    let r = crate::compression::compression_ratio(gauge_pressure, atmospheric_pressure)?;
    Ok(per_cycle * r * n)
}

/// Free-air volumetric flow *demand*, in cubic metres per unit time, for a
/// cylinder cycling at `cycles_per_unit_time`.
///
/// Identical to [`free_air_consumption`] but interpreting the cycle figure
/// as a rate; the result is therefore a flow (e.g. m^3/min if the rate is
/// cycles/min). This is the quantity to match against a compressor's
/// free-air-delivery rating.
///
/// # Errors
///
/// Same as [`free_air_consumption`].
pub fn free_air_flow_demand(
    cylinder: &Cylinder,
    stroke_length: f64,
    action: Action,
    cycles_per_unit_time: f64,
    gauge_pressure: f64,
    atmospheric_pressure: f64,
) -> Result<f64> {
    free_air_consumption(
        cylinder,
        stroke_length,
        action,
        cycles_per_unit_time,
        gauge_pressure,
        atmospheric_pressure,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for floating comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn swept_volume_is_area_times_length() {
        let c = Cylinder::single_acting(0.05).unwrap();
        let v = swept_volume(&c, 0.2, Stroke::Extend).unwrap();
        assert!((v - c.bore_area() * 0.2).abs() < EPS);
    }

    #[test]
    fn consumption_scales_linearly_with_stroke() {
        // Doubling the stroke doubles the consumption (V = A*L).
        let c = Cylinder::single_acting(0.04).unwrap();
        let q1 = free_air_consumption(&c, 0.10, Action::SingleActing, 50.0, 600_000.0, 100_000.0)
            .unwrap();
        let q2 = free_air_consumption(&c, 0.20, Action::SingleActing, 50.0, 600_000.0, 100_000.0)
            .unwrap();
        assert!((q2 - 2.0 * q1).abs() < EPS);
    }

    #[test]
    fn consumption_scales_linearly_with_cycles() {
        // Tripling the cycle count triples the consumption.
        let c = Cylinder::single_acting(0.04).unwrap();
        let q1 = free_air_consumption(&c, 0.15, Action::SingleActing, 100.0, 600_000.0, 100_000.0)
            .unwrap();
        let q3 = free_air_consumption(&c, 0.15, Action::SingleActing, 300.0, 600_000.0, 100_000.0)
            .unwrap();
        assert!((q3 - 3.0 * q1).abs() < EPS);
    }

    #[test]
    fn consumption_scales_with_compression_ratio() {
        // At p_atm = 100 kPa, 6 bar gauge gives r = 7, so the free-air
        // figure is 7x the raw swept volume over the same cycles.
        let c = Cylinder::single_acting(0.05).unwrap();
        let cycles = 100.0;
        let q = free_air_consumption(&c, 0.1, Action::SingleActing, cycles, 600_000.0, 100_000.0)
            .unwrap();
        let raw_swept = c.bore_area() * 0.1 * cycles;
        assert!((q - 7.0 * raw_swept).abs() < 1e-9);
    }

    #[test]
    fn double_acting_uses_more_air_than_single_acting() {
        // The retract stroke adds the annular volume each cycle.
        let c = Cylinder::double_acting(0.05, 0.02).unwrap();
        let single =
            free_air_consumption(&c, 0.1, Action::SingleActing, 100.0, 600_000.0, 100_000.0)
                .unwrap();
        let double =
            free_air_consumption(&c, 0.1, Action::DoubleActing, 100.0, 600_000.0, 100_000.0)
                .unwrap();
        assert!(double > single);

        // The extra is exactly the annular swept volume * r * cycles.
        let extra = 7.0 * c.annulus_area() * 0.1 * 100.0;
        assert!((double - single - extra).abs() < 1e-9);
    }

    #[test]
    fn zero_cycles_consumes_nothing() {
        let c = Cylinder::single_acting(0.05).unwrap();
        let q =
            free_air_consumption(&c, 0.1, Action::SingleActing, 0.0, 600_000.0, 100_000.0).unwrap();
        assert!((q - 0.0).abs() < EPS);
    }

    #[test]
    fn flow_demand_equals_consumption_at_unit_rate() {
        let c = Cylinder::single_acting(0.05).unwrap();
        let q = free_air_consumption(&c, 0.1, Action::DoubleActing, 30.0, 600_000.0, 100_000.0)
            .unwrap();
        let demand =
            free_air_flow_demand(&c, 0.1, Action::DoubleActing, 30.0, 600_000.0, 100_000.0)
                .unwrap();
        assert!((q - demand).abs() < EPS);
    }

    #[test]
    fn free_air_per_stroke_is_swept_times_ratio() {
        let c = Cylinder::single_acting(0.05).unwrap();
        let v = free_air_per_stroke(&c, 0.1, Stroke::Extend, 600_000.0, 100_000.0).unwrap();
        assert!((v - 7.0 * c.bore_area() * 0.1).abs() < 1e-12);
    }

    #[test]
    fn rejects_nonpositive_stroke() {
        let c = Cylinder::single_acting(0.05).unwrap();
        assert!(swept_volume(&c, 0.0, Stroke::Extend).is_err());
        assert!(swept_volume(&c, -0.1, Stroke::Extend).is_err());
    }

    #[test]
    fn rejects_negative_cycles() {
        let c = Cylinder::single_acting(0.05).unwrap();
        let err = free_air_consumption(&c, 0.1, Action::SingleActing, -1.0, 600_000.0, 100_000.0)
            .unwrap_err();
        assert_eq!(err.code(), "pneumatics.negative");
    }

    #[test]
    fn known_value_litre_scale_check() {
        // A 50 mm bore, 100 mm stroke single-acting cylinder at 6 bar
        // gauge (p_atm = 100 kPa, r = 7) per cycle:
        //   V_swept = pi/4 * 0.05^2 * 0.1 = 1.9635e-4 m^3 = 0.19635 L
        //   V_free  = 7 * 0.19635 L       = 1.37445 L of free air.
        let c = Cylinder::single_acting(0.05).unwrap();
        let v_free_m3 =
            free_air_consumption(&c, 0.1, Action::SingleActing, 1.0, 600_000.0, 100_000.0).unwrap();
        let litres = v_free_m3 * 1000.0;
        assert!(
            (litres - 1.374_446_785_310_206).abs() < 1e-9,
            "got {litres} L"
        );
    }
}
