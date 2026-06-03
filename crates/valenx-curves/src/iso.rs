//! Iso-curve extraction from a NURBS surface.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_surface::{NurbsCurve, NurbsSurface};

use crate::error::CurvesError;

/// Which iso-direction to extract.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum IsoKind {
    /// Curve along `v` at fixed `u`.
    IsoU,
    /// Curve along `u` at fixed `v`.
    IsoV,
}

/// Extract an iso-u (constant `u`, sweep `v`) or iso-v (constant
/// `v`, sweep `u`) curve from a surface.
///
/// v1 implementation: sample 64 points along the iso line then fit a
/// NURBS curve with the surface's matching-direction degree. Phase
/// 27.5 will substitute a true row/column extraction that preserves
/// the surface's CP grid exactly.
pub fn extract_iso(
    surface: &NurbsSurface,
    kind: IsoKind,
    parameter: f64,
) -> Result<NurbsCurve, CurvesError> {
    let u_min = surface.u_knots[surface.u_degree];
    let u_max = surface.u_knots[surface.u_knots.len() - 1 - surface.u_degree];
    let v_min = surface.v_knots[surface.v_degree];
    let v_max = surface.v_knots[surface.v_knots.len() - 1 - surface.v_degree];
    let n_samples = 64usize;
    let mut pts: Vec<Vector3<f64>> = Vec::with_capacity(n_samples);
    let (degree, n_cps_max) = match kind {
        IsoKind::IsoU => (surface.v_degree, n_samples),
        IsoKind::IsoV => (surface.u_degree, n_samples),
    };
    if !(u_min..=u_max).contains(&parameter) && !(v_min..=v_max).contains(&parameter) {
        return Err(CurvesError::BadParameter {
            name: "parameter",
            reason: format!("not in surface u or v domain ({parameter})"),
        });
    }
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let p = match kind {
            IsoKind::IsoU => {
                let v = v_min + t * (v_max - v_min);
                surface.evaluate(parameter, v)
            }
            IsoKind::IsoV => {
                let u = u_min + t * (u_max - u_min);
                surface.evaluate(u, parameter)
            }
        };
        pts.push(p);
    }
    let n_cps = (degree + 1).max(n_cps_max.min(16));
    let fit = valenx_surface::fit::nurbs_curve_through_points(&pts, degree, n_cps)?;
    Ok(fit.curve)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iso_out_of_range_errors() {
        // Build a trivial 2x2 bilinear plane surface for the test —
        // smallest valid surface.
        let surface = NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::zeros(), Vector3::new(0.0, 1.0, 0.0)],
                vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .unwrap();
        // 2.0 is outside both [0,1] domains.
        assert!(extract_iso(&surface, IsoKind::IsoU, 2.0).is_err());
    }

    #[test]
    fn iso_u_at_half_returns_v_line() {
        let surface = NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::zeros(), Vector3::new(0.0, 1.0, 0.0)],
                vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .unwrap();
        let iso = extract_iso(&surface, IsoKind::IsoU, 0.5).unwrap();
        // Curve should be a horizontal line at x = 0.5 going y: 0..1.
        let u_min = iso.knots[iso.degree];
        let u_max = iso.knots[iso.knots.len() - 1 - iso.degree];
        let p0 = iso.evaluate(u_min);
        let p1 = iso.evaluate(u_max);
        assert!((p0.x - 0.5).abs() < 0.01);
        assert!((p1.x - 0.5).abs() < 0.01);
    }
}
