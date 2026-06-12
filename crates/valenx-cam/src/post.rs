//! G-code postprocessors — convert a [`crate::Toolpath`] into
//! controller-specific text.
//!
//! Three flavors ship in v1:
//!
//! - [`grbl::Grbl`] — most popular hobby CNC firmware. Single-tool
//!   (no `T{n} M6` support).
//! - [`linuxcnc::LinuxCnc`] — full G-code dialect with tool change.
//! - [`fanuc::Fanuc`] — industrial controller. N-numbered lines,
//!   `T{n} M6` tool change.
//!
//! Every postprocessor implements [`Postprocessor`]. The
//! [`Postprocessor::process`] entry point composes header → moves →
//! footer into one final string.
//!
//! See [`save_nc`] for the file-output helper.

pub mod fanuc;
pub mod grbl;
pub mod linuxcnc;

// Phase 17C — templated 27 postprocessor variants.
pub mod haas;
#[macro_use]
pub mod macros;
pub mod phase17;
pub mod template;

use std::path::Path;

use nalgebra::Vector3;

use crate::{error::CamError, tool::Tool, toolpath::Toolpath};

/// Which postprocessor variant the host has selected. Used by the
/// CAM panel's dropdown and the [`process_with_kind`] dispatcher.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum PostKind {
    /// GRBL — single-tool hobby firmware.
    #[default]
    Grbl,
    /// LinuxCNC — full G-code with tool change.
    LinuxCnc,
    /// Fanuc — N-numbered industrial controller.
    Fanuc,
    // ----- Phase 17C variants -----
    /// Haas CNC controllers.
    Haas,
    /// Heidenhain TNC.
    Heidenhain,
    /// Mazak Mazatrol.
    Mazatrol,
    /// Siemens Sinumerik 840D.
    Sinumerik,
    /// Okuma OSP.
    Okuma,
    /// Mori-Seiki MAPPS.
    Mori,
    /// Makino Professional 5.
    Makino,
    /// Kitamura Mycenter.
    Kitamura,
    /// Centroid hobby controller.
    Centroid,
    /// Tormach PCNC.
    Tormach,
    /// DeepNest sheet-metal flat-pattern.
    DeepNest,
    /// Marlin 3D-printer firmware.
    Marlin,
    /// Klipper firmware.
    Klipper,
    /// Repetier firmware.
    Repetier,
    /// FluidNC firmware.
    FluidNc,
    /// SnapMaker combo firmware.
    SnapMaker,
    /// SmoothieWare firmware.
    Smoothie,
    /// Smoothieboard firmware.
    Smoothieboard,
    /// TinyG controller.
    TinyG,
    /// SourceRabbit (Mach3/4 derivative).
    SourceRabbit,
    /// HSMAdvisor cycle-time-aware post.
    HsmAdvisor,
    /// Fusion 360 branded post.
    Fusion360Brand,
    /// OpenDMG (DMG Mori clone).
    OpenDmg,
    /// VMC 4-axis (A-axis).
    VmcAAxis,
    /// VMC 4-axis (B-axis).
    VmcBAxis,
    /// Fanuc 4-axis (A-axis).
    FanucAAxis,
    /// Fanuc 5-axis (A+B).
    FanucBAxis,
}

impl PostKind {
    /// Short label for panels.
    pub fn label(self) -> &'static str {
        match self {
            PostKind::Grbl => "GRBL",
            PostKind::LinuxCnc => "LinuxCNC",
            PostKind::Fanuc => "Fanuc",
            PostKind::Haas => "Haas",
            PostKind::Heidenhain => "Heidenhain TNC",
            PostKind::Mazatrol => "Mazatrol",
            PostKind::Sinumerik => "Sinumerik 840D",
            PostKind::Okuma => "Okuma OSP",
            PostKind::Mori => "Mori MAPPS",
            PostKind::Makino => "Makino",
            PostKind::Kitamura => "Kitamura",
            PostKind::Centroid => "Centroid",
            PostKind::Tormach => "Tormach",
            PostKind::DeepNest => "DeepNest",
            PostKind::Marlin => "Marlin",
            PostKind::Klipper => "Klipper",
            PostKind::Repetier => "Repetier",
            PostKind::FluidNc => "FluidNC",
            PostKind::SnapMaker => "SnapMaker",
            PostKind::Smoothie => "SmoothieWare",
            PostKind::Smoothieboard => "Smoothieboard",
            PostKind::TinyG => "TinyG",
            PostKind::SourceRabbit => "SourceRabbit",
            PostKind::HsmAdvisor => "HSMAdvisor",
            PostKind::Fusion360Brand => "Fusion 360 (branded)",
            PostKind::OpenDmg => "OpenDMG",
            PostKind::VmcAAxis => "VMC 4-axis (A)",
            PostKind::VmcBAxis => "VMC 4-axis (B)",
            PostKind::FanucAAxis => "Fanuc 4-axis (A)",
            PostKind::FanucBAxis => "Fanuc 5-axis (A+B)",
        }
    }
}

/// Trait every CNC postprocessor implements. The host calls
/// [`Postprocessor::process`] to convert a [`Toolpath`] to G-code.
pub trait Postprocessor {
    /// Initialisation block emitted before any moves (units mode,
    /// absolute positioning, spindle on, etc.).
    fn header(&self, tool: &Tool, spindle_rpm: f64) -> String;

    /// Wind-down block emitted after the last move (spindle off,
    /// program end).
    fn footer(&self) -> String;

    /// Format a rapid traverse (G0) to `p`.
    fn move_g0(&self, p: Vector3<f64>) -> String;

    /// Format a feed move (G1) to `p` at `feed` mm/min.
    fn move_g1(&self, p: Vector3<f64>, feed: f64) -> Result<String, CamError>;

    /// Format a spindle-on (M3 S{rpm}) line.
    fn spindle_on(&self, rpm: f64) -> String;

    /// Format a spindle-off (M5) line.
    fn spindle_off(&self) -> String;

    /// Format a tool-change (`T{n} M6`) line, or empty string if the
    /// controller is single-tool.
    fn tool_change(&self, tool_id: u32) -> String;

    /// Format a coolant-on (M8) line.
    fn coolant_on(&self) -> String;

    /// Format a coolant-off (M9) line.
    fn coolant_off(&self) -> String;

    /// Convert a complete toolpath to a G-code string.
    fn process(
        &self,
        toolpath: &Toolpath,
        tool: &Tool,
        spindle_rpm: f64,
    ) -> Result<String, CamError>;

    /// Format a 5-axis feed move with tool centre at `p` and rotary
    /// axes at `a` / `b` degrees. Default impl emits a warning
    /// comment for 3-axis postprocessors.
    ///
    /// Phase 17E posts that natively support 5-axis (fanuc_a_axis /
    /// fanuc_b_axis / vmc_a_axis / vmc_b_axis) override this to emit
    /// real `A{deg} B{deg}` lines.
    fn move_g1_5ax(&self, p: Vector3<f64>, a: f64, b: f64, feed: f64) -> Result<String, CamError> {
        let _ = (p, a, b, feed);
        Ok("; 5-axis move dropped — postprocessor lacks A/B output".into())
    }

    /// Format an XY-plane circular arc (G2 = clockwise, G3 = CCW)
    /// from the previous position to `end`, with `centre_xy` the arc
    /// centre. The standard G-code form emits I/J as the *relative*
    /// offset from the arc start to the centre, so the postprocessor
    /// needs the start position too.
    ///
    /// Default impl uses `format_g23` which produces the
    /// standard-controller-agnostic `G2|G3 X Y Z I J F` line.
    fn move_g23(
        &self,
        start: Vector3<f64>,
        end: Vector3<f64>,
        centre_xy: nalgebra::Vector2<f64>,
        dir: crate::arcfit::ArcDir,
        feed: f64,
    ) -> Result<String, CamError> {
        format_g23(start, end, centre_xy, dir, feed)
    }
}

/// Convenience dispatcher — runs [`Postprocessor::process`] for the
/// chosen [`PostKind`] without the caller having to import each impl.
pub fn process_with_kind(
    kind: PostKind,
    toolpath: &Toolpath,
    tool: &Tool,
    spindle_rpm: f64,
) -> Result<String, CamError> {
    match kind {
        PostKind::Grbl => grbl::Grbl.process(toolpath, tool, spindle_rpm),
        PostKind::LinuxCnc => linuxcnc::LinuxCnc.process(toolpath, tool, spindle_rpm),
        PostKind::Fanuc => fanuc::Fanuc.process(toolpath, tool, spindle_rpm),
        PostKind::Haas => haas::Haas.process(toolpath, tool, spindle_rpm),
        PostKind::Heidenhain => {
            phase17::heidenhain::Heidenhain.process(toolpath, tool, spindle_rpm)
        }
        PostKind::Mazatrol => phase17::mazatrol::Mazatrol.process(toolpath, tool, spindle_rpm),
        PostKind::Sinumerik => phase17::sinumerik::Sinumerik.process(toolpath, tool, spindle_rpm),
        PostKind::Okuma => phase17::okuma::Okuma.process(toolpath, tool, spindle_rpm),
        PostKind::Mori => phase17::mori::Mori.process(toolpath, tool, spindle_rpm),
        PostKind::Makino => phase17::makino::Makino.process(toolpath, tool, spindle_rpm),
        PostKind::Kitamura => phase17::kitamura::Kitamura.process(toolpath, tool, spindle_rpm),
        PostKind::Centroid => phase17::centroid::Centroid.process(toolpath, tool, spindle_rpm),
        PostKind::Tormach => phase17::tormach::Tormach.process(toolpath, tool, spindle_rpm),
        PostKind::DeepNest => phase17::deepnest::DeepNest.process(toolpath, tool, spindle_rpm),
        PostKind::Marlin => phase17::marlin::Marlin.process(toolpath, tool, spindle_rpm),
        PostKind::Klipper => phase17::klipper::Klipper.process(toolpath, tool, spindle_rpm),
        PostKind::Repetier => phase17::repetier::Repetier.process(toolpath, tool, spindle_rpm),
        PostKind::FluidNc => phase17::fluid_nc::FluidNc.process(toolpath, tool, spindle_rpm),
        PostKind::SnapMaker => phase17::snap_maker::SnapMaker.process(toolpath, tool, spindle_rpm),
        PostKind::Smoothie => phase17::smoothie::Smoothie.process(toolpath, tool, spindle_rpm),
        PostKind::Smoothieboard => {
            phase17::smoothieboard::Smoothieboard.process(toolpath, tool, spindle_rpm)
        }
        PostKind::TinyG => phase17::tinyg::TinyG.process(toolpath, tool, spindle_rpm),
        PostKind::SourceRabbit => {
            phase17::sourcerabbit::SourceRabbit.process(toolpath, tool, spindle_rpm)
        }
        PostKind::HsmAdvisor => {
            phase17::hsmadvisor::HsmAdvisor.process(toolpath, tool, spindle_rpm)
        }
        PostKind::Fusion360Brand => {
            phase17::fusion360_brand::Fusion360Brand.process(toolpath, tool, spindle_rpm)
        }
        PostKind::OpenDmg => phase17::opendmg::OpenDmg.process(toolpath, tool, spindle_rpm),
        PostKind::VmcAAxis => phase17::vmc_a_axis::VmcAAxis.process(toolpath, tool, spindle_rpm),
        PostKind::VmcBAxis => phase17::vmc_b_axis::VmcBAxis.process(toolpath, tool, spindle_rpm),
        PostKind::FanucAAxis => {
            phase17::fanuc_a_axis::FanucAAxis.process(toolpath, tool, spindle_rpm)
        }
        PostKind::FanucBAxis => {
            phase17::fanuc_b_axis::FanucBAxis.process(toolpath, tool, spindle_rpm)
        }
    }
}

/// Format a 5-axis G1 move (template helper). Used by Fanuc 4/5-axis
/// + VMC 4-axis postprocessors that override [`Postprocessor::move_g1_5ax`].
pub fn format_g1_5ax(p: Vector3<f64>, a: f64, b: f64, feed: f64) -> Result<String, CamError> {
    // Round-3 fix: reject not just `feed <= 0` but also NaN / ±inf —
    // `!(feed > 0.0)` is true for negative AND for NaN but `+inf > 0.0`
    // is true and would produce `F inf` (which most controllers parse
    // as "unlimited rapid", and on some machines that translates to a
    // smashed spindle).
    if !feed.is_finite() || feed <= 0.0 {
        return Err(CamError::PostprocessorFailed {
            reason: format!("cut feed must be a finite number > 0 (got {feed})"),
        });
    }
    if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!(
                "G1_5ax target must be finite (got x={}, y={}, z={})",
                p.x, p.y, p.z
            ),
        });
    }
    if !(a.is_finite() && b.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!("G1_5ax rotary angles must be finite (got a={a}, b={b})"),
        });
    }
    Ok(format!(
        "G1 X{:.3} Y{:.3} Z{:.3} A{:.3} B{:.3} F{:.0}",
        p.x, p.y, p.z, a, b, feed
    ))
}

/// Run [`Postprocessor::process`] and write the result to `path`.
///
/// Bytes are written in UTF-8 with whatever line endings the
/// postprocessor produced (LF for all v1 impls — GRBL/LinuxCNC/Fanuc
/// all accept LF or CRLF).
pub fn save_nc(
    kind: PostKind,
    toolpath: &Toolpath,
    tool: &Tool,
    spindle_rpm: f64,
    path: &Path,
) -> Result<(), CamError> {
    let gcode = process_with_kind(kind, toolpath, tool, spindle_rpm)?;
    valenx_core::io_caps::atomic_write_str(path, &gcode)?;
    Ok(())
}

/// Shared helper: format a feed move with .3f-precision G1 syntax.
///
/// Round-3 fix: rejects not just `feed <= 0` but also NaN / ±inf, and
/// requires finite XYZ. The previous `!(feed > 0.0)` chain caught NaN
/// (any comparison with NaN is false, so the negation fires) and
/// negative values — but `+inf > 0.0` is true, so a `+inf` feed used
/// to slip through and emit `F inf` to the .nc file. Most controllers
/// parse that as "unlimited rapid", which on some machines smashes the
/// spindle.
pub(crate) fn format_g1(p: Vector3<f64>, feed: f64) -> Result<String, CamError> {
    if !feed.is_finite() || feed <= 0.0 {
        return Err(CamError::PostprocessorFailed {
            reason: format!("cut feed must be a finite number > 0 (got {feed})"),
        });
    }
    if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!(
                "G1 target must be finite (got x={}, y={}, z={})",
                p.x, p.y, p.z
            ),
        });
    }
    Ok(format!(
        "G1 X{:.3} Y{:.3} Z{:.3} F{:.0}",
        p.x, p.y, p.z, feed
    ))
}

/// Shared helper: format a rapid move with .3f precision.
///
/// Round-3 fix: when any of X/Y/Z is NaN/±inf we replace it with `0.0`
/// rather than emitting `X NaN` (which most controllers parse as zero
/// silently, sending the spindle straight through stock at rapid
/// speed). Returning a fault-tolerant string here keeps the trait
/// signature `String` and lets the Postprocessor::process loop's
/// upstream validators decide whether to bail. The fallible
/// [`format_g0_strict`] is available for callers that prefer hard
/// failures.
pub(crate) fn format_g0(p: Vector3<f64>) -> String {
    let safe = |v: f64| if v.is_finite() { v } else { 0.0 };
    format!("G0 X{:.3} Y{:.3} Z{:.3}", safe(p.x), safe(p.y), safe(p.z))
}

/// Like [`format_g0`] but returns [`CamError::PostprocessorFailed`]
/// when any coordinate is NaN / ±inf, instead of substituting zeros.
/// Use this from validators that need to surface "the toolpath is
/// corrupt" upstream rather than silently emit safe-but-wrong G-code.
///
/// Currently unused by the in-tree postprocessors (they all rely on
/// the fault-tolerant [`format_g0`]); kept available for future
/// callers and exercised by tests in this module.
#[allow(dead_code)]
pub(crate) fn format_g0_strict(p: Vector3<f64>) -> Result<String, CamError> {
    if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!(
                "G0 target must be finite (got x={}, y={}, z={})",
                p.x, p.y, p.z
            ),
        });
    }
    Ok(format!("G0 X{:.3} Y{:.3} Z{:.3}", p.x, p.y, p.z))
}

/// Shared helper: format a G2/G3 circular-arc move. `start` is the
/// previous tool position (used to compute the I/J offsets); `end`
/// is this move's target; `centre_xy` is the arc centre in XY;
/// `dir` selects G2 (CW) vs G3 (CCW).
///
/// I/J are the offset from arc *start* to centre, the universal
/// modern G-code convention (G91.1 incremental I/J interpretation
/// is the default on every contemporary controller — Fanuc,
/// Heidenhain, Siemens, LinuxCNC, GRBL).
pub fn format_g23(
    start: Vector3<f64>,
    end: Vector3<f64>,
    centre_xy: nalgebra::Vector2<f64>,
    dir: crate::arcfit::ArcDir,
    feed: f64,
) -> Result<String, CamError> {
    // Round-14 M5 (round-3 sister gap for arc moves): pre-fix the
    // `!(feed > 0.0)` check let NaN slip (since NaN compares neither
    // true nor false to 0.0, the negated comparison evaluated true
    // — but only by accident; `+inf > 0.0` was `true` so it passed),
    // and the start / end / centre_xy components were never
    // validated. A NaN or +inf in any of those flows straight into
    // the `{:.3}` format specifier as `NaN` / `inf`, which most
    // controllers parse as zero — driving the spindle to (0,0,0) at
    // a feed rate of "unlimited" if the feed itself was inf.
    if !feed.is_finite() || feed <= 0.0 {
        return Err(CamError::PostprocessorFailed {
            reason: format!("arc feed must be finite and > 0 (got {feed})"),
        });
    }
    if !start.iter().all(|v| v.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!("arc start position must be finite (got {start:?})"),
        });
    }
    if !end.iter().all(|v| v.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!("arc end position must be finite (got {end:?})"),
        });
    }
    if !centre_xy.iter().all(|v| v.is_finite()) {
        return Err(CamError::PostprocessorFailed {
            reason: format!("arc centre_xy must be finite (got {centre_xy:?})"),
        });
    }
    let code = match dir {
        crate::arcfit::ArcDir::Clockwise => "G2",
        crate::arcfit::ArcDir::Counterclockwise => "G3",
    };
    let i = centre_xy.x - start.x;
    let j = centre_xy.y - start.y;
    Ok(format!(
        "{} X{:.3} Y{:.3} Z{:.3} I{:.3} J{:.3} F{:.0}",
        code, end.x, end.y, end.z, i, j, feed
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// Round-3 fix: a positive-infinity feed used to produce `F inf`
    /// in the .nc output — most controllers parse that as "unlimited
    /// rapid" and on some machines that smashes the spindle.
    #[test]
    fn format_g1_rejects_infinity_feed() {
        let result = format_g1(Vector3::new(0.0, 0.0, 0.0), f64::INFINITY);
        let err = result.expect_err("inf feed must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    #[test]
    fn format_g1_rejects_nan_feed() {
        let result = format_g1(Vector3::new(0.0, 0.0, 0.0), f64::NAN);
        let err = result.expect_err("NaN feed must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    #[test]
    fn format_g1_rejects_negative_feed() {
        let result = format_g1(Vector3::new(0.0, 0.0, 0.0), -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn format_g1_accepts_finite_positive_feed() {
        let result = format_g1(Vector3::new(1.0, 2.0, 3.0), 500.0);
        assert_eq!(result.unwrap(), "G1 X1.000 Y2.000 Z3.000 F500");
    }

    /// Round-3 fix: NaN/inf coordinates must not slip through and emit
    /// `X NaN` / `X inf` strings (which most controllers parse as
    /// zero).
    #[test]
    fn format_g1_rejects_nan_position() {
        let result = format_g1(Vector3::new(f64::NAN, 0.0, 0.0), 500.0);
        let err = result.expect_err("NaN position must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    /// `format_g0` must NOT emit `X NaN` / `X inf` — we substitute 0.0
    /// for the offending coordinate as a fault-tolerant fallback, since
    /// the trait shape is `String` (not `Result`). The NaN substring
    /// must not survive into the .nc text.
    #[test]
    fn format_g0_substitutes_zero_for_nan_position() {
        let g0 = format_g0(Vector3::new(f64::NAN, 1.0, 2.0));
        assert!(!g0.contains("NaN"), "G0 must not emit NaN: {g0}");
        assert!(!g0.contains("inf"), "G0 must not emit inf: {g0}");
        assert!(g0.contains("X0.000"), "expected X substituted to 0: {g0}");
    }

    #[test]
    fn format_g0_substitutes_zero_for_inf_position() {
        let g0 = format_g0(Vector3::new(0.0, f64::INFINITY, f64::NEG_INFINITY));
        assert!(!g0.contains("inf"), "G0 must not emit inf: {g0}");
        assert!(g0.contains("Y0.000"));
        assert!(g0.contains("Z0.000"));
    }

    #[test]
    fn format_g0_strict_rejects_nan_position() {
        let result = format_g0_strict(Vector3::new(f64::NAN, 0.0, 0.0));
        let err = result.expect_err("NaN position must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    /// Round-14 M5 RED→GREEN: arc moves (G2 / G3) must reject an
    /// infinite feed for the same reason the round-3 G1 fix did —
    /// most controllers parse `F inf` as "unlimited rapid" and the
    /// spindle can collide with the workpiece. Pre-fix the
    /// `!(feed > 0.0)` check let +inf through because `+inf > 0.0`
    /// is `true`; the new `is_finite() && > 0.0` shape catches it.
    #[test]
    fn format_g23_rejects_infinity_feed() {
        let result = format_g23(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            nalgebra::Vector2::new(0.5, 0.5),
            crate::arcfit::ArcDir::Clockwise,
            f64::INFINITY,
        );
        let err = result.expect_err("inf feed must be rejected for arc");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    /// Round-14 M5 sister: NaN in any start / end / centre position
    /// must be rejected. Pre-fix it slipped through and the format
    /// string emitted `X NaN`, which most controllers read as 0.0
    /// — drives the cutter into the fixture.
    #[test]
    fn format_g23_rejects_nan_start_position() {
        let result = format_g23(
            Vector3::new(f64::NAN, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            nalgebra::Vector2::new(0.5, 0.5),
            crate::arcfit::ArcDir::Counterclockwise,
            500.0,
        );
        let err = result.expect_err("NaN start must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
        assert!(msg.contains("start"), "msg: {msg}");
    }
}
