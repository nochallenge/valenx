//! Helix polyline + 3D thread profile sweep.

use std::f64::consts::PI;

use nalgebra::Vector3;

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::ThreadsProError;
use crate::spec::ThreadSpecPro;
use crate::standard::ProfileShape;

/// Sample a true helical polyline matching the thread spec. Returns a
/// list of `n_turns * samples_per_turn + 1` points spiralling along
/// +Z, starting at `(major/2, 0, 0)`.
pub fn helix_polyline(
    spec: &ThreadSpecPro,
    n_turns: u32,
    samples_per_turn: u32,
) -> Result<Vec<Vector3<f64>>, ThreadsProError> {
    if n_turns == 0 {
        return Err(ThreadsProError::BadParameter {
            name: "n_turns",
            reason: "must be >= 1".into(),
        });
    }
    if samples_per_turn < 8 {
        return Err(ThreadsProError::BadParameter {
            name: "samples_per_turn",
            reason: "must be >= 8 for a recognisable helix".into(),
        });
    }
    let r = spec.major_diameter() / 2.0;
    let total = n_turns as usize * samples_per_turn as usize + 1;
    let mut out = Vec::with_capacity(total);
    for i in 0..total {
        let t = i as f64 / samples_per_turn as f64; // turns
        let theta = t * 2.0 * PI;
        let z = t * spec.pitch;
        out.push(Vector3::new(r * theta.cos(), r * theta.sin(), z));
    }
    Ok(out)
}

/// Build a mesh-backed [`Solid`] approximating the helical thread
/// profile swept along the axis. v1 emits a coarse band of triangles
/// (one tooth per turn, profile shape encoded in radial offset) so the
/// result is visually correct but is **not** a manifold solid — the
/// purpose is panel preview + downstream mesh ops.
pub fn profile_solid(spec: &ThreadSpecPro, length: f64) -> Result<Solid, ThreadsProError> {
    if length <= 0.0 {
        return Err(ThreadsProError::BadParameter {
            name: "length",
            reason: format!("must be > 0, got {length}"),
        });
    }
    // Pitch must be finite and positive: a zero/negative/NaN pitch makes
    // `length / pitch` non-finite, and `(+inf).ceil() as u32` saturates to
    // u32::MAX -- `helix_polyline` would then request a ~u32::MAX * 32
    // element (multi-terabyte) allocation and abort. `ThreadSpecPro` is
    // `Deserialize` with no field validation, so a garbled saved spec can
    // carry `pitch = 0`.
    if !spec.pitch.is_finite() || spec.pitch <= 0.0 {
        return Err(ThreadsProError::BadParameter {
            name: "pitch",
            reason: format!("must be finite and > 0, got {}", spec.pitch),
        });
    }
    let samples_per_turn: u32 = 32;
    // Reject an absurd turn count *before* the `as u32` cast: a tiny but
    // positive pitch would otherwise still drive `n_turns` to billions and
    // OOM. A real thread is far under this cap (the longest standard
    // fastener is well under 1000 mm at >= 0.2 mm pitch -- a few thousand
    // turns).
    const MAX_TURNS: f64 = 100_000.0;
    let turns = (length / spec.pitch).ceil();
    if turns > MAX_TURNS {
        return Err(ThreadsProError::BadParameter {
            name: "pitch",
            reason: format!(
                "length/pitch = {turns} turns exceeds the {MAX_TURNS} cap (pitch too small for this length)"
            ),
        });
    }
    let n_turns = (turns as u32).max(1);
    let helix = helix_polyline(spec, n_turns, samples_per_turn)?;

    let r_major = spec.major_diameter() / 2.0;
    let r_minor = spec.minor_diameter() / 2.0;
    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(helix.len() * 2);
    for (i, p) in helix.iter().enumerate() {
        // Profile shape modulates the radial offset within one turn.
        let t = (i % samples_per_turn as usize) as f64 / samples_per_turn as f64;
        let r = profile_radial(spec.profile, t, r_minor, r_major);
        let theta = (i as f64 / samples_per_turn as f64) * 2.0 * PI;
        // Crest (outer) + root (inner) ribbon.
        nodes.push(Vector3::new(r * theta.cos(), r * theta.sin(), p.z));
        nodes.push(Vector3::new(
            r_minor * theta.cos(),
            r_minor * theta.sin(),
            p.z,
        ));
    }

    let mut mesh = Mesh::new(format!("thread_{}", spec.designation));
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    let n = helix.len();
    for i in 0..(n - 1) {
        let a = (2 * i) as u32;
        let b = (2 * i + 1) as u32;
        let c = (2 * (i + 1)) as u32;
        let d = (2 * (i + 1) + 1) as u32;
        // Two tris per quad.
        block.connectivity.extend_from_slice(&[a, b, c]);
        block.connectivity.extend_from_slice(&[b, d, c]);
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Radial offset within one turn for the given profile shape. `t` is
/// the fraction of the turn `[0, 1)`. `r_minor` / `r_major` are the
/// inside / outside radii.
fn profile_radial(profile: ProfileShape, t: f64, r_minor: f64, r_major: f64) -> f64 {
    match profile {
        // Symmetric V triangle — peak at the centre of the turn.
        ProfileShape::V => {
            let s = if t < 0.5 { 2.0 * t } else { 2.0 * (1.0 - t) };
            r_minor + (r_major - r_minor) * s
        }
        // Acme / Trapezoidal — flat-topped trapezium occupying the
        // middle 50 % of the turn.
        ProfileShape::Acme | ProfileShape::Trapezoidal => {
            if !(0.25..=0.75).contains(&t) {
                r_minor
            } else if (0.25..0.375).contains(&t) {
                let s = (t - 0.25) / 0.125;
                r_minor + (r_major - r_minor) * s
            } else if (0.625..0.75).contains(&t) {
                let s = (0.75 - t) / 0.125;
                r_minor + (r_major - r_minor) * s
            } else {
                r_major
            }
        }
        // Buttress — steep face followed by a 7° back face. Approximate
        // with an asymmetric triangle peaking at t = 0.3.
        ProfileShape::Buttress => {
            if t < 0.3 {
                let s = t / 0.3;
                r_minor + (r_major - r_minor) * s
            } else {
                let s = (1.0 - t) / 0.7;
                r_minor + (r_major - r_minor) * s.max(0.0)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::standard::ThreadStandardPro;

    #[test]
    fn helix_starts_on_positive_x_axis() {
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        let pts = helix_polyline(&spec, 2, 16).unwrap();
        assert!((pts[0].x - 4.0).abs() < 1e-9);
        assert!(pts[0].y.abs() < 1e-9);
        assert!(pts[0].z.abs() < 1e-9);
    }

    #[test]
    fn helix_advances_z_by_pitch_per_turn() {
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        let pts = helix_polyline(&spec, 3, 32).unwrap();
        // First sample past the third turn = index 3 * 32.
        let z = pts[3 * 32].z;
        assert!((z - 3.0 * 1.25).abs() < 1e-9, "got z = {z}");
    }

    #[test]
    fn helix_rejects_zero_turns() {
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        let err = helix_polyline(&spec, 0, 32).unwrap_err();
        assert!(matches!(err, ThreadsProError::BadParameter { .. }));
    }

    #[test]
    fn profile_solid_returns_non_empty_mesh_backed_solid() {
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        let solid = profile_solid(&spec, 10.0).unwrap();
        // Mesh-backed solid → faces() returns 0 but the underlying
        // mesh has triangles. Check via tessellation:
        let mesh = valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE).unwrap();
        assert!(mesh.nodes.len() > 10);
    }

    #[test]
    fn profile_solid_rejects_zero_length() {
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        let err = profile_solid(&spec, 0.0).unwrap_err();
        assert!(matches!(err, ThreadsProError::BadParameter { .. }));
    }

    #[test]
    fn profile_solid_rejects_zero_pitch() {
        // A zero pitch makes length/pitch non-finite; pre-fix the
        // `(+inf).ceil() as u32` saturated to u32::MAX and the helix
        // requested a multi-terabyte allocation (abort). It must be a
        // clean BadParameter instead. `ThreadSpecPro` is `Deserialize`
        // with no field validation, so a garbled saved spec can carry
        // pitch = 0.
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 0.0);
        let err = profile_solid(&spec, 10.0).unwrap_err();
        assert!(
            matches!(&err, ThreadsProError::BadParameter { name, .. } if *name == "pitch"),
            "expected a pitch BadParameter, got {err:?}"
        );
    }

    #[test]
    fn profile_solid_rejects_pitch_too_small_for_length() {
        // A tiny but positive pitch would drive n_turns to billions and
        // OOM the same way; the turn-count cap must reject it before the
        // `as u32` cast.
        let spec = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1e-9);
        let err = profile_solid(&spec, 10.0).unwrap_err();
        assert!(
            matches!(&err, ThreadsProError::BadParameter { name, .. } if *name == "pitch"),
            "expected a pitch BadParameter, got {err:?}"
        );
    }
}
