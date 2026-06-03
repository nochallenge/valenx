//! Snell's-law refraction at an interface line.
//!
//! Phase 12B Task 26 / 12.6. Models the geometric form of Snell's
//! law: `n1 · sin(θ_i) = n2 · sin(θ_t)` where the angles are measured
//! between each ray (`ray_in`, `ray_out`) and the normal of the
//! `interface_line`.
//!
//! ## Exact residual — no small-angle approximation
//!
//! The residual is the **exact** Snell's-law mismatch
//! `n1·sin θ_i − n2·sin θ_t`. The sines are computed directly as the
//! 2D cross product of the unit ray direction with the unit interface
//! normal — `sin θ = û_ray × n̂` — which is the true sine for any
//! angle in `(−90°, 90°)`, not the `sin θ ≈ θ` linearisation. The
//! constraint therefore converges correctly for steep incidence
//! (total-internal-reflection geometry aside, which has no real
//! refracted solution and is the caller's responsibility to avoid).

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Residual r = n1 * sin(theta_i) - n2 * sin(theta_t).
///
/// `ray_in` and `ray_out` are lines whose endpoints define the
/// incoming and outgoing rays. The interface line provides the
/// surface normal (perpendicular to its direction).
pub fn residuals(
    sketch: &Sketch,
    ray_in: EntityId,
    ray_out: EntityId,
    interface_line: EntityId,
    n1: f64,
    n2: f64,
    out: &mut [f64],
) {
    let (Ok(lin), Ok(lout), Ok(iface)) = (
        sketch.line_at(ray_in),
        sketch.line_at(ray_out),
        sketch.line_at(interface_line),
    ) else {
        out[0] = 0.0;
        return;
    };
    let (din_x, din_y) = lin.direction(&sketch.vars);
    let (dout_x, dout_y) = lout.direction(&sketch.vars);
    let (dif_x, dif_y) = iface.direction(&sketch.vars);
    let nlen = (dif_x * dif_x + dif_y * dif_y).sqrt();
    if nlen < 1e-15 {
        out[0] = 0.0;
        return;
    }
    // Surface normal = perpendicular of interface direction.
    let nx = -dif_y / nlen;
    let ny = dif_x / nlen;
    // sin(angle to normal) = |cross(unit_ray, unit_normal)|.
    let l_in = (din_x * din_x + din_y * din_y).sqrt();
    let l_out = (dout_x * dout_x + dout_y * dout_y).sqrt();
    if l_in < 1e-15 || l_out < 1e-15 {
        out[0] = 0.0;
        return;
    }
    let sin_i = (din_x / l_in) * ny - (din_y / l_in) * nx;
    let sin_t = (dout_x / l_out) * ny - (dout_y / l_out) * nx;
    out[0] = n1 * sin_i - n2 * sin_t;
}

/// Jacobian: finite-difference (12 candidate variables).
pub fn jacobian(
    sketch: &Sketch,
    ray_in: EntityId,
    ray_out: EntityId,
    interface_line: EntityId,
    n1: f64,
    n2: f64,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(lin), Ok(lout), Ok(iface)) = (
        sketch.line_at(ray_in),
        sketch.line_at(ray_out),
        sketch.line_at(interface_line),
    ) else {
        return;
    };
    let vars = vec![
        lin.start.x_var,
        lin.start.y_var,
        lin.end.x_var,
        lin.end.y_var,
        lout.start.x_var,
        lout.start.y_var,
        lout.end.x_var,
        lout.end.y_var,
        iface.start.x_var,
        iface.start.y_var,
        iface.end.x_var,
        iface.end.y_var,
    ];
    let mut base = vec![0.0; 1];
    residuals(sketch, ray_in, ray_out, interface_line, n1, n2, &mut base);
    let r0 = base[0];
    let h = 1e-7;
    let mut perturbed = sketch.clone();
    for var in vars {
        let saved = perturbed.vars[var];
        perturbed.vars[var] = saved + h;
        let mut r = vec![0.0; 1];
        residuals(&perturbed, ray_in, ray_out, interface_line, n1, n2, &mut r);
        let d = (r[0] - r0) / h;
        if d.abs() > 1e-15 {
            triplets.push((0, var, d));
        }
        perturbed.vars[var] = saved;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_for_normal_incidence() {
        // Ray going straight down hits a horizontal interface; refracted
        // ray also goes straight down. theta_i = theta_t = 0 → residual
        // = 0 regardless of n1, n2.
        let mut s = Sketch::new();
        let p_in_top = s.add_point(0.0, 2.0);
        let p_hit = s.add_point(0.0, 0.0);
        let p_out_bottom = s.add_point(0.0, -2.0);
        let p_iface_l = s.add_point(-1.0, 0.0);
        let p_iface_r = s.add_point(1.0, 0.0);
        let ray_in = s.add_line(p_in_top, p_hit).unwrap();
        let ray_out = s.add_line(p_hit, p_out_bottom).unwrap();
        let iface = s.add_line(p_iface_l, p_iface_r).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, ray_in, ray_out, iface, 1.5, 1.0, &mut out);
        assert!(out[0].abs() < 1e-12, "got {}", out[0]);
    }

    #[test]
    fn residual_uses_exact_sine_at_steep_angle() {
        // 12.6: verify the residual is the *exact* sine, not sin θ ≈ θ.
        // Incident ray at 60° to the (vertical) normal of a horizontal
        // interface; refracted ray at the angle Snell's law predicts
        // for n1=1.0, n2=1.5. The residual must be ~0 — which only
        // holds if the true sine (not the linearisation) is used.
        let mut s = Sketch::new();
        let theta_i = 60f64.to_radians();
        // n1 sin θ_i = n2 sin θ_t  →  sin θ_t = (1.0/1.5) sin 60°.
        let theta_t = ((1.0 / 1.5) * theta_i.sin()).asin();
        // Interface horizontal → normal is vertical (±y). Rays point
        // downward; measure angle from the -y normal.
        let p_hit = s.add_point(0.0, 0.0);
        // Incident ray comes in from upper-left at θ_i from vertical.
        let p_in = s.add_point(-theta_i.sin(), theta_i.cos());
        // Refracted ray continues to lower-right at θ_t from vertical.
        let p_out = s.add_point(theta_t.sin(), -theta_t.cos());
        let p_iface_l = s.add_point(-1.0, 0.0);
        let p_iface_r = s.add_point(1.0, 0.0);
        let ray_in = s.add_line(p_in, p_hit).unwrap();
        let ray_out = s.add_line(p_hit, p_out).unwrap();
        let iface = s.add_line(p_iface_l, p_iface_r).unwrap();
        let mut out = vec![0.0; 1];
        residuals(&s, ray_in, ray_out, iface, 1.0, 1.5, &mut out);
        assert!(
            out[0].abs() < 1e-9,
            "exact-sine Snell residual should vanish at 60°, got {}",
            out[0]
        );
        // A small-angle (sin θ ≈ θ) model would give a residual of
        // n1·θ_i − n2·θ_t, which is markedly non-zero here:
        let small_angle_residual = 1.0 * theta_i - 1.5 * theta_t;
        assert!(
            small_angle_residual.abs() > 0.05,
            "sanity: small-angle model genuinely differs here"
        );
    }
}
