//! 2D tooth-profile math. Standard involute construction:
//!
//! ```text
//! x(theta) = r_b * (cos(theta) + theta * sin(theta))
//! y(theta) = r_b * (sin(theta) - theta * cos(theta))
//! ```
//!
//! Where `r_b` is the base-circle radius. `theta = 0` lands on the
//! base circle along the +x axis, and `theta` increases tangentially.

use crate::error::GearsError;
use crate::spec::GearSpec;

/// Single point on the involute generated from `base_radius` at
/// parameter `theta` (radians).
pub fn involute_point(base_radius: f64, theta: f64) -> [f64; 2] {
    let cs = theta.cos();
    let sn = theta.sin();
    [
        base_radius * (cs + theta * sn),
        base_radius * (sn - theta * cs),
    ]
}

/// Sample one full tooth (one tooth + one space) profile in 2D.
///
/// The output is a closed polyline starting at the dedendum
/// circle, going up the leading involute flank, across the
/// addendum top-land, down the trailing flank, and across the
/// dedendum bottom-land back to the start. v1 samples 8 points
/// per flank and 3 across the addendum/dedendum lands.
pub fn tooth_profile(spec: &GearSpec) -> Result<Vec<[f64; 2]>, GearsError> {
    if spec.module_mm <= 0.0 {
        return Err(GearsError::BadParameter {
            name: "module_mm",
            reason: format!("must be > 0, got {}", spec.module_mm),
        });
    }
    if spec.teeth < 5 {
        return Err(GearsError::TeethTooFew(spec.teeth));
    }
    if !(0.0..=45.0).contains(&spec.pressure_angle_deg) {
        return Err(GearsError::BadParameter {
            name: "pressure_angle_deg",
            reason: format!(
                "must be in [0, 45], got {}",
                spec.pressure_angle_deg
            ),
        });
    }

    let r_pitch = spec.pitch_diameter_mm() * 0.5;
    let r_base = spec.base_diameter_mm() * 0.5;
    let r_add = spec.addendum_diameter_mm() * 0.5;
    let r_ded = spec.dedendum_diameter_mm() * 0.5;
    let n = spec.teeth as f64;

    // Maximum involute parameter reaches the addendum circle.
    let max_theta = ((r_add / r_base).powi(2) - 1.0).max(0.0).sqrt();

    // Half-tooth angle at the pitch circle. Classical formula:
    // pi / (2n) - inv(alpha) where alpha = pressure angle.
    let alpha = spec.pressure_angle_deg.to_radians();
    let inv_alpha = alpha.tan() - alpha;
    let half_tooth_angle = std::f64::consts::FRAC_PI_2 / n - inv_alpha;

    let steps = 8;
    let mut leading: Vec<[f64; 2]> = Vec::with_capacity(steps + 1);
    let mut trailing: Vec<[f64; 2]> = Vec::with_capacity(steps + 1);
    for i in 0..=steps {
        let t = max_theta * (i as f64 / steps as f64);
        // Leading flank starts from base circle.
        let p = involute_point(r_base, t);
        // Rotate the involute so that it crosses the pitch circle at
        // -half_tooth_angle (a standard "tooth straddles +x" layout).
        // Find pitch-crossing theta_p:
        let theta_p = ((r_pitch / r_base).powi(2) - 1.0).max(0.0).sqrt();
        let p_pitch = involute_point(r_base, theta_p);
        let phi = p_pitch[1].atan2(p_pitch[0]) + half_tooth_angle;
        let rot = -phi;
        let cs = rot.cos();
        let sn = rot.sin();
        let q = [p[0] * cs - p[1] * sn, p[0] * sn + p[1] * cs];
        leading.push(q);
        // Trailing flank = mirror across +x axis.
        trailing.push([q[0], -q[1]]);
    }

    // Build the closed tooth polyline.
    // Start at dedendum on the +x axis (between this tooth and the
    // next one).
    let mut out: Vec<[f64; 2]> = Vec::new();
    // 1) Dedendum point at the start of the leading flank.
    let theta_ded_start = -std::f64::consts::PI / n;
    if r_ded > 0.0 {
        out.push([
            r_ded * theta_ded_start.cos(),
            r_ded * theta_ded_start.sin(),
        ]);
    }
    // 2) Leading flank from base → addendum. (Trailing list is
    // mirrored about +x; we want the flank closer to the trailing
    // edge first in the layout going CCW, then back over the
    // addendum then down the leading flank.)
    for p in trailing.iter() {
        out.push(*p);
    }
    // 3) Addendum top-land — sample 3 points.
    let top_steps = 3;
    let top_start = trailing.last().unwrap();
    let top_end = leading.last().unwrap();
    let a0 = top_start[1].atan2(top_start[0]);
    let a1 = top_end[1].atan2(top_end[0]);
    for i in 1..top_steps {
        let t = i as f64 / top_steps as f64;
        let a = a0 + (a1 - a0) * t;
        out.push([r_add * a.cos(), r_add * a.sin()]);
    }
    // 4) Leading flank (already CCW order from base → addendum, so
    // reverse for the downward sweep).
    for p in leading.iter().rev() {
        out.push(*p);
    }
    // 5) Dedendum point at the start of the next tooth's leading
    // flank.
    let theta_ded_end = std::f64::consts::PI / n;
    if r_ded > 0.0 {
        out.push([r_ded * theta_ded_end.cos(), r_ded * theta_ded_end.sin()]);
    }

    Ok(out)
}

/// Replicate the tooth polyline around the gear `teeth` times to
/// build the complete 2D outline.
pub fn full_profile(spec: &GearSpec) -> Result<Vec<[f64; 2]>, GearsError> {
    let tooth = tooth_profile(spec)?;
    let n = spec.teeth as f64;
    let delta = std::f64::consts::TAU / n;
    let mut out: Vec<[f64; 2]> = Vec::with_capacity(tooth.len() * spec.teeth as usize);
    for i in 0..spec.teeth {
        let phi = i as f64 * delta;
        let cs = phi.cos();
        let sn = phi.sin();
        for p in &tooth {
            out.push([p[0] * cs - p[1] * sn, p[0] * sn + p[1] * cs]);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn involute_zero_theta_lands_on_base_circle() {
        let p = involute_point(10.0, 0.0);
        assert!((p[0] - 10.0).abs() < 1e-9);
        assert!(p[1].abs() < 1e-9);
    }

    #[test]
    fn tooth_profile_grows_with_addendum() {
        let spec = GearSpec::standard_spur(20);
        let pts = tooth_profile(&spec).unwrap();
        assert!(!pts.is_empty());
        let max_r = pts
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt())
            .fold(0.0_f64, f64::max);
        // Should reach close to the addendum radius.
        let r_add = spec.addendum_diameter_mm() * 0.5;
        assert!((max_r - r_add).abs() < 0.05 * r_add);
    }

    #[test]
    fn full_profile_has_teeth_times_tooth_length() {
        let spec = GearSpec::standard_spur(20);
        let tooth = tooth_profile(&spec).unwrap();
        let full = full_profile(&spec).unwrap();
        assert_eq!(full.len(), tooth.len() * spec.teeth as usize);
    }

    #[test]
    fn rejects_too_few_teeth() {
        let mut spec = GearSpec::standard_spur(3);
        spec.teeth = 3;
        assert!(matches!(
            tooth_profile(&spec),
            Err(GearsError::TeethTooFew(3))
        ));
    }
}
