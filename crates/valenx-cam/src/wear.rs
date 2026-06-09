//! Tool wear model (Phase 17F).
//!
//! v1 uses Taylor's tool-life equation: `V × T^n = C`, where:
//!
//! - `V` = cutting speed (m/min, derived from rpm + diameter).
//! - `T` = tool life (minutes).
//! - `n` = Taylor exponent (material-specific, 0.1–0.4 typical).
//! - `C` = Taylor constant (material-specific).
//!
//! Material constants are kept in [`material_constants`] keyed on a
//! lower-cased substring match (the host's free-form material
//! descriptor on the [`crate::stock::Stock`] is matched
//! case-insensitively).
//!
//! v1 limitations:
//!
//! - **No chip-load awareness** — the equation ignores feed/tooth,
//!   coolant, and depth-of-cut. The output is an order-of-magnitude
//!   estimate, not a CAM-vendor-grade tool-life prediction.
//! - **No tool-coating multipliers** — coated carbide and uncoated
//!   carbide use the same constant.

use std::time::Duration;

use serde::{Deserialize, Serialize};

use crate::tool::Tool;

/// Material-specific Taylor constants. Keyed on a lowercased
/// substring of the stock material descriptor; first match wins.
pub fn material_constants(material: &str) -> (f64, f64) {
    // Returns (C, n).
    let m = material.to_lowercase();
    if m.contains("alu") {
        (300.0, 0.30)
    } else if m.contains("brass") || m.contains("copper") {
        (250.0, 0.25)
    } else if m.contains("steel") || m.contains("iron") {
        (60.0, 0.18)
    } else if m.contains("titanium") || m.contains("ti-") {
        (35.0, 0.15)
    } else if m.contains("plastic") || m.contains("nylon") || m.contains("hdpe") {
        (500.0, 0.40)
    } else if m.contains("wood") || m.contains("mdf") || m.contains("ply") {
        (800.0, 0.45)
    } else {
        // Generic carbon-steel fallback.
        (60.0, 0.18)
    }
}

/// Cutting speed in m/min for a given diameter (mm) and rpm.
pub fn cutting_speed_m_per_min(diameter_mm: f64, rpm: f64) -> f64 {
    std::f64::consts::PI * (diameter_mm / 1000.0) * rpm
}

/// Feed per tooth (chip load) `f_z = v_f / (n · z)` (mm/tooth) — the advance of the
/// workpiece per cutting edge, the key parameter for chatter-free milling: too low rubs
/// and work-hardens, too high overloads the edge. `feed_mm_per_min` is the table feed
/// `v_f`, `rpm` the spindle speed `n`, and `teeth` the number of flutes `z`. Returns `0`
/// for non-physical input (non-finite, or `rpm` / `teeth` non-positive).
pub fn feed_per_tooth(feed_mm_per_min: f64, rpm: f64, teeth: f64) -> f64 {
    if !feed_mm_per_min.is_finite()
        || !rpm.is_finite()
        || rpm <= 0.0
        || !teeth.is_finite()
        || teeth <= 0.0
    {
        return 0.0;
    }
    feed_mm_per_min / (rpm * teeth)
}

/// Operation summary needed for the wear model.
#[derive(Clone, Copy, Debug)]
pub struct OpRunSpec {
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Estimated minutes the tool is engaged in this op.
    pub estimated_minutes: f64,
}

/// Saturating ceiling for [`tool_life`]'s Duration. f64 can only
/// represent integers exactly up to 2^53; we clamp well below that
/// so the f64-to-u64 round-trip in `Duration::from_secs_f64` lands
/// on a finite, representable value. 1e15 seconds is ~31 million
/// years — far longer than any real tool life.
const TOOL_LIFE_SECS_MAX: f64 = 1.0e15;

/// Predicted tool life (minutes) under Taylor's equation.
///
/// Round-8 hardening: the result is clamped against
/// [`Duration::from_secs_f64`]'s panic surface. `from_secs_f64`
/// panics on NaN, negatives, or seconds outside the [0, u64::MAX]
/// range — and `powf(1.0 / n)` over extreme inputs can produce any of
/// those. We clamp to `[0, TOOL_LIFE_SECS_MAX]` seconds and
/// substitute `Duration::ZERO` for NaN.
pub fn tool_life(tool: &Tool, op: OpRunSpec, material: &str) -> Duration {
    let v = cutting_speed_m_per_min(tool.diameter_mm, op.spindle_rpm).max(1e-6);
    let (c, n) = material_constants(material);
    // V * T^n = C  ⇒  T = (C / V)^(1/n).
    let t_min = (c / v).powf(1.0 / n);
    let secs = t_min * 60.0;
    if !secs.is_finite() {
        // NaN or ±inf — Duration::from_secs_f64 panics on these.
        // Substitute zero (no predicted life under degenerate inputs).
        return Duration::ZERO;
    }
    // Clamp into [0, TOOL_LIFE_SECS_MAX] so from_secs_f64 stays
    // well clear of the saturation edge on extreme inputs. We
    // guard against NaN above; clamp is safe here.
    let clamped = secs.clamp(0.0, TOOL_LIFE_SECS_MAX);
    Duration::from_secs_f64(clamped)
}

/// Warning emitted by [`check_op`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WearWarning {
    /// Code (kebab-cased) — stable across versions for LLM dispatch.
    pub code: String,
    /// Human-readable explanation.
    pub message: String,
}

/// Check whether the operation exceeds the predicted tool life.
pub fn check_op(tool: &Tool, op: OpRunSpec, material: &str) -> Vec<WearWarning> {
    let life = tool_life(tool, op, material);
    let life_min = life.as_secs_f64() / 60.0;
    let mut out = Vec::new();
    if op.estimated_minutes > life_min {
        out.push(WearWarning {
            code: "wear.life_exceeded".into(),
            message: format!(
                "operation runs {:.1} min — predicted tool life is {:.1} min",
                op.estimated_minutes, life_min
            ),
        });
    } else if op.estimated_minutes > life_min * 0.75 {
        out.push(WearWarning {
            code: "wear.life_warning".into(),
            message: format!(
                "operation runs {:.1} min — within 75% of predicted tool life ({:.1} min)",
                op.estimated_minutes, life_min
            ),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;

    #[test]
    fn cutting_speed_basic() {
        // 10mm @ 10000 rpm = pi * 0.01 * 10000 = ~314 m/min.
        let v = cutting_speed_m_per_min(10.0, 10000.0);
        assert!((v - 314.159).abs() < 0.1);
    }

    #[test]
    fn feed_per_tooth_basic() {
        // v_f = 800 mm/min, n = 10000 rpm, z = 4 flutes → f_z = 800/(10000·4) = 0.02 mm/tooth.
        let fz = feed_per_tooth(800.0, 10000.0, 4.0);
        assert!((fz - 0.02).abs() < 1e-12, "f_z = v_f/(n·z), got {fz}");
        // Non-tautological thread: the table feed reconstructs from f_z (v_f = f_z·n·z).
        assert!((fz * 10000.0 * 4.0 - 800.0).abs() < 1e-9, "v_f = f_z·n·z");
        // Inversely proportional to flute count and to rpm.
        assert!((feed_per_tooth(800.0, 10000.0, 2.0) - 2.0 * fz).abs() < 1e-12, "∝ 1/z");
        assert!((feed_per_tooth(800.0, 5000.0, 4.0) - 2.0 * fz).abs() < 1e-12, "∝ 1/n");
        // Non-physical input → 0.
        assert_eq!(feed_per_tooth(800.0, 0.0, 4.0), 0.0);
        assert_eq!(feed_per_tooth(800.0, 10000.0, 0.0), 0.0);
        assert_eq!(feed_per_tooth(f64::NAN, 10000.0, 4.0), 0.0);
    }

    #[test]
    fn material_constants_alu() {
        let (c, n) = material_constants("6061 aluminum");
        assert!(c > 200.0);
        assert!(n > 0.2 && n < 0.4);
    }

    #[test]
    fn tool_life_steel_short() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let life = tool_life(
            &tool,
            OpRunSpec {
                spindle_rpm: 10000.0,
                estimated_minutes: 0.0,
            },
            "1045 steel",
        );
        // Steel at 10000 rpm on a 6mm cutter is aggressive — life
        // should be on the order of a few minutes (not infinite).
        let life_min = life.as_secs_f64() / 60.0;
        assert!(life_min > 0.0 && life_min < 60.0, "got {life_min} min");
    }

    #[test]
    fn check_op_warns_when_exceeded() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let warns = check_op(
            &tool,
            OpRunSpec {
                spindle_rpm: 12000.0,
                estimated_minutes: 1e6,
            },
            "steel",
        );
        assert!(!warns.is_empty());
        assert_eq!(warns[0].code, "wear.life_exceeded");
    }

    #[test]
    fn check_op_clean_when_safe() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let warns = check_op(
            &tool,
            OpRunSpec {
                spindle_rpm: 5000.0,
                estimated_minutes: 0.001,
            },
            "wood",
        );
        assert!(warns.is_empty(), "expected clean check, got {warns:?}");
    }

    #[test]
    fn tool_life_extreme_inputs_do_not_panic() {
        // Round-8 RED→GREEN: pre-fix, `tool_life` could feed a NaN or
        // out-of-range f64 into `Duration::from_secs_f64`, which
        // panics — a single hostile tool spec could crash the host
        // process. With the saturating clamp in place we get back a
        // finite Duration on every f64 input.
        let tool = Tool::new(1, "EM_huge", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        // Extremely low rpm pushes t_min toward infinity — the
        // clamp must keep the result finite.
        let life = tool_life(
            &tool,
            OpRunSpec {
                spindle_rpm: 1e-300,
                estimated_minutes: 0.0,
            },
            "steel",
        );
        // Should be at least non-zero; the cap is TOOL_LIFE_SECS_MAX.
        assert!(
            life.as_secs() <= TOOL_LIFE_SECS_MAX as u64 + 1,
            "tool_life exceeded the clamp ceiling: {} s",
            life.as_secs()
        );
        // Zero RPM also pre-fix triggered panic in some material
        // configurations. Verify it now returns a finite Duration.
        let life0 = tool_life(
            &tool,
            OpRunSpec {
                spindle_rpm: 0.0,
                estimated_minutes: 0.0,
            },
            "steel",
        );
        assert!(life0.as_secs() <= TOOL_LIFE_SECS_MAX as u64 + 1);
    }
}
