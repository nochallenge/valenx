//! Centreline emitters + tessellated swept solid + stiffness.

use nalgebra::Vector3;

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::SpringsError;
use crate::spec::{SpringKind, SpringSpec};

const SAMPLES_PER_COIL: usize = 24;

fn validate(spec: &SpringSpec) -> Result<(), SpringsError> {
    if spec.wire_diameter_mm <= 0.0 {
        return Err(SpringsError::BadParameter {
            name: "wire_diameter_mm",
            reason: format!("must be > 0, got {}", spec.wire_diameter_mm),
        });
    }
    if spec.mean_coil_diameter_mm <= spec.wire_diameter_mm {
        return Err(SpringsError::Degenerate(format!(
            "mean coil dia ({}) must be > wire dia ({})",
            spec.mean_coil_diameter_mm, spec.wire_diameter_mm
        )));
    }
    if spec.free_length_mm <= 0.0 {
        return Err(SpringsError::BadParameter {
            name: "free_length_mm",
            reason: format!("must be > 0, got {}", spec.free_length_mm),
        });
    }
    if spec.n_active_coils <= 0.0 {
        return Err(SpringsError::BadParameter {
            name: "n_active_coils",
            reason: format!("must be > 0, got {}", spec.n_active_coils),
        });
    }
    Ok(())
}

/// Helical centreline for a compression spring.
pub fn compression_centerline(spec: &SpringSpec) -> Result<Vec<Vector3<f64>>, SpringsError> {
    validate(spec)?;
    let r = spec.mean_coil_diameter_mm * 0.5;
    let n = spec.n_active_coils;
    let pitch = spec.pitch_mm();
    let total_samples = (SAMPLES_PER_COIL as f64 * n).ceil() as usize;
    let mut pts = Vec::with_capacity(total_samples + 1);
    for s in 0..=total_samples {
        let t = s as f64 / SAMPLES_PER_COIL as f64; // coil parameter
        let phi = t * std::f64::consts::TAU;
        let z = -spec.free_length_mm * 0.5 + t * pitch;
        pts.push(Vector3::new(r * phi.cos(), r * phi.sin(), z));
    }
    Ok(pts)
}

/// Centreline for an extension spring — same helix as compression
/// but with a hook loop at each end (semicircle in the y-z plane
/// returning the wire perpendicular to the axis).
pub fn extension_centerline(spec: &SpringSpec) -> Result<Vec<Vector3<f64>>, SpringsError> {
    let mut pts = compression_centerline(spec)?;
    let r = spec.mean_coil_diameter_mm * 0.5;
    let hook_steps = 12;
    // Start hook — semicircle in y-z plane around the first point.
    let p0 = pts[0];
    let mut start_hook = Vec::with_capacity(hook_steps + 1);
    for s in 0..=hook_steps {
        let theta = (s as f64 / hook_steps as f64) * std::f64::consts::PI;
        let z = p0.z - r + r * theta.cos();
        let y = -r * theta.sin();
        start_hook.push(Vector3::new(p0.x, y, z));
    }
    // End hook.
    let pn = *pts.last().unwrap();
    let mut end_hook = Vec::with_capacity(hook_steps + 1);
    for s in 0..=hook_steps {
        let theta = (s as f64 / hook_steps as f64) * std::f64::consts::PI;
        let z = pn.z + r - r * theta.cos();
        let y = r * theta.sin();
        end_hook.push(Vector3::new(pn.x, y, z));
    }
    let mut full = Vec::with_capacity(start_hook.len() + pts.len() + end_hook.len());
    full.extend(start_hook.into_iter().rev());
    full.append(&mut pts);
    full.extend(end_hook);
    Ok(full)
}

/// Centreline for a torsion spring — compression helix plus
/// straight radial legs at each end (length = mean coil diameter).
pub fn torsion_centerline(spec: &SpringSpec) -> Result<Vec<Vector3<f64>>, SpringsError> {
    let helix = compression_centerline(spec)?;
    let r = spec.mean_coil_diameter_mm * 0.5;
    let leg_len = r * 2.0;
    let leg_steps = 8;
    let p0 = helix[0];
    let pn = *helix.last().unwrap();
    // Start leg pointing in +x direction outward from the helix at p0.
    let dir0 = Vector3::new(p0.x, p0.y, 0.0).normalize();
    let mut start_leg = Vec::with_capacity(leg_steps + 1);
    for s in 0..=leg_steps {
        let t = s as f64 / leg_steps as f64;
        start_leg.push(p0 + dir0 * leg_len * (1.0 - t));
    }
    let dir_n = Vector3::new(pn.x, pn.y, 0.0).normalize();
    let mut end_leg = Vec::with_capacity(leg_steps + 1);
    for s in 0..=leg_steps {
        let t = s as f64 / leg_steps as f64;
        end_leg.push(pn + dir_n * leg_len * t);
    }
    let mut full = Vec::with_capacity(start_leg.len() + helix.len() + end_leg.len());
    full.extend(start_leg);
    full.extend(helix);
    full.extend(end_leg);
    Ok(full)
}

/// Sweep the wire-diameter circle along the centreline picked from
/// `spec.kind`. Returns a [`Solid::Mesh`].
pub fn to_solid(spec: &SpringSpec) -> Result<Solid, SpringsError> {
    let centreline = match spec.kind {
        SpringKind::Compression => compression_centerline(spec)?,
        SpringKind::Extension => extension_centerline(spec)?,
        SpringKind::Torsion => torsion_centerline(spec)?,
    };
    if centreline.len() < 2 {
        return Err(SpringsError::Degenerate(
            "centreline has < 2 points".into(),
        ));
    }

    // 12-sided wire cross-section.
    let n = 12usize;
    let r_wire = spec.wire_diameter_mm * 0.5;
    let mut mesh = Mesh::new(format!("springs_{}", spec.kind.label()));
    let mut block = ElementBlock::new(ElementType::Tri3);
    let mut ring_base = Vec::with_capacity(centreline.len());

    let world_up = Vector3::new(0.0, 0.0, 1.0);
    let world_alt = Vector3::new(0.0, 1.0, 0.0);
    for i in 0..centreline.len() {
        let tangent = if i == 0 {
            centreline[1] - centreline[0]
        } else if i == centreline.len() - 1 {
            centreline[i] - centreline[i - 1]
        } else {
            (centreline[i + 1] - centreline[i - 1]).normalize()
        };
        let t = tangent.normalize();
        let ref_up = if t.dot(&world_up).abs() > 0.95 {
            world_alt
        } else {
            world_up
        };
        let u = ref_up.cross(&t).normalize();
        let v = t.cross(&u).normalize();
        let base = mesh.nodes.len() as u32;
        ring_base.push(base);
        for k in 0..n {
            let phi = (k as f64 / n as f64) * std::f64::consts::TAU;
            let p = centreline[i] + u * (r_wire * phi.cos()) + v * (r_wire * phi.sin());
            mesh.nodes.push(p);
        }
    }
    for w in 0..(centreline.len() - 1) {
        let a = ring_base[w];
        let b = ring_base[w + 1];
        for k in 0..n {
            let j = (k + 1) % n;
            block.connectivity.extend_from_slice(&[
                a + k as u32,
                a + j as u32,
                b + j as u32,
            ]);
            block.connectivity.extend_from_slice(&[
                a + k as u32,
                b + j as u32,
                b + k as u32,
            ]);
        }
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Linear axial stiffness for a compression spring:
///
/// `k = G * d^4 / (8 * D^3 * n_active)`
///
/// where `G` is the shear modulus (MPa = N/mm²), `d` the wire
/// diameter, `D` the mean coil diameter, and `n_active` the active
/// coil count. Result is in N/mm.
pub fn stiffness_n_per_mm(spec: &SpringSpec) -> f64 {
    let g = spec.shear_modulus_mpa;
    let d = spec.wire_diameter_mm;
    let big_d = spec.mean_coil_diameter_mm;
    let n = spec.n_active_coils.max(1e-6);
    g * d.powi(4) / (8.0 * big_d.powi(3) * n)
}

/// Spring index `C = D / d` (ISO 26909) — the ratio of mean coil diameter to wire diameter, a
/// dimensionless measure of coil curvature (typically 4–12). It feeds the Wahl stress-correction
/// factor and buckling checks. Returns `0.0` for a non-positive or non-finite wire diameter.
pub fn spring_index(spec: &SpringSpec) -> f64 {
    if !spec.wire_diameter_mm.is_finite() || spec.wire_diameter_mm <= 0.0 {
        return 0.0;
    }
    spec.mean_coil_diameter_mm / spec.wire_diameter_mm
}

/// Wahl stress-correction factor `K_w = (4C − 1)/(4C − 4) + 0.615/C` for a helical spring of
/// index `C = D/d` (see [`spring_index`]). It corrects the shear stress for coil curvature and
/// direct shear; K_w decreases toward 1 as C grows. Returns `0.0` for `C ≤ 1` (the 4C−4
/// denominator vanishes at C = 1) or non-finite input.
pub fn wahl_factor(c: f64) -> f64 {
    if !c.is_finite() || c <= 1.0 {
        return 0.0;
    }
    (4.0 * c - 1.0) / (4.0 * c - 4.0) + 0.615 / c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compression_centreline_has_expected_length() {
        let spec = SpringSpec::default_compression();
        let pts = compression_centerline(&spec).unwrap();
        let expected = (SAMPLES_PER_COIL as f64 * spec.n_active_coils).ceil() as usize + 1;
        assert_eq!(pts.len(), expected);
    }

    #[test]
    fn helix_length_matches_developed_arc() {
        let spec = SpringSpec::default_compression();
        // n=8, D=10, p=30/8=3.75 → ℓ = 8·√((π·10)² + 3.75²) = 8·√(986.96+14.06) ≈ 253.11 mm.
        let expected =
            8.0 * ((std::f64::consts::PI * 10.0).powi(2) + (30.0_f64 / 8.0).powi(2)).sqrt();
        assert!((spec.helix_length_mm() - expected).abs() < 1e-9);
        assert!((spec.helix_length_mm() - 253.111).abs() < 0.01);
        // The coiled wire is far longer than the straight free length.
        assert!(spec.helix_length_mm() > spec.free_length_mm);
        // Non-positive active coils → 0 (guard, no panic, no division).
        let mut degen = SpringSpec::default_compression();
        degen.n_active_coils = 0.0;
        assert_eq!(degen.helix_length_mm(), 0.0);
    }

    #[test]
    fn extension_adds_two_hooks() {
        let spec = SpringSpec::default_compression();
        let pts_c = compression_centerline(&spec).unwrap();
        let pts_e = extension_centerline(&spec).unwrap();
        assert!(pts_e.len() > pts_c.len() + 20);
    }

    #[test]
    fn torsion_adds_two_legs() {
        let spec = SpringSpec::default_compression();
        let pts_c = compression_centerline(&spec).unwrap();
        let pts_t = torsion_centerline(&spec).unwrap();
        assert!(pts_t.len() > pts_c.len() + 10);
    }

    #[test]
    fn to_solid_emits_mesh_backed() {
        let spec = SpringSpec::default_compression();
        let s = to_solid(&spec).unwrap();
        match s {
            Solid::Mesh(m) => {
                assert!(!m.nodes.is_empty());
                assert!(m.total_elements() > 0);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn stiffness_matches_canonical_formula() {
        // Steel, d=1, D=10, n=8 -> k = 79300 / (8 * 1000 * 8)
        //                            = 79300 / 64000
        //                            ≈ 1.2391 N/mm
        let spec = SpringSpec::default_compression();
        let k = stiffness_n_per_mm(&spec);
        assert!((k - 1.2390625).abs() < 1e-6);
    }

    #[test]
    fn rejects_wire_bigger_than_coil() {
        let mut spec = SpringSpec::default_compression();
        spec.wire_diameter_mm = 12.0;
        assert!(matches!(
            compression_centerline(&spec),
            Err(SpringsError::Degenerate(_))
        ));
    }

    #[test]
    fn spring_index_is_ratio_of_diameters() {
        // default_compression: D = 10, d = 1 → C = 10.
        let spec = SpringSpec::default_compression();
        assert!((spring_index(&spec) - 10.0).abs() < 1e-9);
        // C scales linearly with D and inversely with d.
        let mut s2 = SpringSpec::default_compression();
        s2.mean_coil_diameter_mm = 20.0;
        assert!((spring_index(&s2) - 20.0).abs() < 1e-9);
        s2.wire_diameter_mm = 2.0;
        assert!((spring_index(&s2) - 10.0).abs() < 1e-9); // 20 / 2
        // Guard: non-positive wire diameter → 0.0.
        let mut bad = SpringSpec::default_compression();
        bad.wire_diameter_mm = 0.0;
        assert_eq!(spring_index(&bad), 0.0);
    }

    #[test]
    fn wahl_factor_corrects_for_curvature() {
        // C = 10 → K_w = 39/36 + 0.0615 ≈ 1.144833.
        assert!((wahl_factor(10.0) - 1.144_833_333).abs() < 1e-5);
        // C = 6 → 23/20 + 0.1025 = 1.2525.
        assert!((wahl_factor(6.0) - 1.2525).abs() < 1e-9);
        // K_w decreases toward 1 as C grows.
        assert!(wahl_factor(6.0) > wahl_factor(12.0));
        assert!(wahl_factor(12.0) > 1.0);
        // Guards: C ≤ 1 (denominator vanishes) or non-finite → 0.
        assert_eq!(wahl_factor(1.0), 0.0);
        assert_eq!(wahl_factor(0.5), 0.0);
        assert_eq!(wahl_factor(f64::NAN), 0.0);
    }
}
