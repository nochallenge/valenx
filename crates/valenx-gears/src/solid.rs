//! 3D gear solid emitters. Each `to_solid_*` returns a
//! [`valenx_cad::Solid::Mesh`] sized for the viewport / BOM cut-list
//! display.

use nalgebra::Vector3;

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::GearsError;
use crate::profile::full_profile;
use crate::spec::{GearKind, GearSpec};

/// Spur gear — extrude the 2D outline straight along +z.
pub fn to_solid_spur(spec: &GearSpec) -> Result<Solid, GearsError> {
    if spec.face_width_mm <= 0.0 {
        return Err(GearsError::BadParameter {
            name: "face_width_mm",
            reason: format!("must be > 0, got {}", spec.face_width_mm),
        });
    }
    let profile2d = full_profile(spec)?;
    extrude_twisted(&profile2d, spec.face_width_mm, 0.0, "gears_spur")
}

/// Helical gear — twisted extrusion of the 2D outline along +z.
/// Twist amount = helix angle × face width / pitch radius (so the
/// helix is constant-pitch along the cylinder).
pub fn to_solid_helical(spec: &GearSpec) -> Result<Solid, GearsError> {
    if !matches!(spec.kind, GearKind::Helical) {
        return Err(GearsError::BadParameter {
            name: "kind",
            reason: format!("expected Helical, got {:?}", spec.kind),
        });
    }
    if spec.helix_angle_deg.abs() < 1e-6 {
        return Err(GearsError::Degenerate(
            "helix_angle_deg = 0 — use to_solid_spur instead".into(),
        ));
    }
    let profile2d = full_profile(spec)?;
    let r_p = spec.pitch_diameter_mm() * 0.5;
    let twist_rad = spec.helix_angle_deg.to_radians() * spec.face_width_mm / r_p;
    extrude_twisted(&profile2d, spec.face_width_mm, twist_rad, "gears_helical")
}

/// Bevel gear — conical sweep. Each station along the cone slope
/// scales the 2D profile to its local pitch-cone radius.
pub fn to_solid_bevel(spec: &GearSpec) -> Result<Solid, GearsError> {
    if spec.face_width_mm <= 0.0 {
        return Err(GearsError::BadParameter {
            name: "face_width_mm",
            reason: format!("must be > 0, got {}", spec.face_width_mm),
        });
    }
    let profile2d = full_profile(spec)?;
    let r_p = spec.pitch_diameter_mm() * 0.5;
    // Standard 45° pitch cone for v1.
    let pitch_cone_angle = std::f64::consts::FRAC_PI_4;
    let mut mesh = Mesh::new("gears_bevel");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let n_layers = 6;
    let n = profile2d.len();
    let mut ring_base: Vec<u32> = Vec::with_capacity(n_layers + 1);
    for layer in 0..=n_layers {
        let t = layer as f64 / n_layers as f64; // 0..=1
        let z = -spec.face_width_mm * 0.5 + spec.face_width_mm * t;
        // Scale factor: 1 at z=−half, shrinks as +z along the cone.
        let scale = 1.0 - (spec.face_width_mm * t) / (r_p / pitch_cone_angle.tan() + 1e-9);
        let scale = scale.max(0.05);
        let base = mesh.nodes.len() as u32;
        ring_base.push(base);
        for p in &profile2d {
            mesh.nodes.push(Vector3::new(p[0] * scale, p[1] * scale, z));
        }
    }
    for layer in 0..n_layers {
        for i in 0..n {
            let j = (i + 1) % n;
            let a = ring_base[layer] + i as u32;
            let b = ring_base[layer] + j as u32;
            let c = ring_base[layer + 1] + j as u32;
            let d = ring_base[layer + 1] + i as u32;
            block.connectivity.extend_from_slice(&[a, b, c]);
            block.connectivity.extend_from_slice(&[a, c, d]);
        }
    }
    // Caps.
    cap_polygon(&mut mesh, &mut block, ring_base[0], n, true);
    cap_polygon(&mut mesh, &mut block, *ring_base.last().unwrap(), n, false);
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Worm gear — single-start helical thread sweep around the worm
/// axis. v1 emits a cylindrical shaft with a triangular thread
/// running along it (3 turns by default).
pub fn to_solid_worm(spec: &GearSpec) -> Result<Solid, GearsError> {
    if spec.face_width_mm <= 0.0 {
        return Err(GearsError::BadParameter {
            name: "face_width_mm",
            reason: format!("must be > 0, got {}", spec.face_width_mm),
        });
    }
    let r_p = spec.pitch_diameter_mm() * 0.5;
    let r_root = (r_p - spec.module_mm).max(spec.module_mm * 0.1);
    let r_tip = r_p + spec.module_mm;
    let mut mesh = Mesh::new("gears_worm");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let segments = 64;
    let turns = (spec.face_width_mm / (spec.module_mm * std::f64::consts::PI)).max(1.0);
    let helix_lead = spec.module_mm * std::f64::consts::PI; // axial advance per turn
    let total_steps = (segments as f64 * turns) as usize;
    // Two parallel polylines along the helix: root + tip.
    let mut root_base = Vec::with_capacity(total_steps);
    let mut tip_base = Vec::with_capacity(total_steps);
    for s in 0..=total_steps {
        let t = s as f64 / segments as f64;
        let phi = t * std::f64::consts::TAU;
        let z = -spec.face_width_mm * 0.5 + t * helix_lead;
        if z > spec.face_width_mm * 0.5 {
            break;
        }
        root_base.push(mesh.nodes.len() as u32);
        mesh.nodes
            .push(Vector3::new(r_root * phi.cos(), r_root * phi.sin(), z));
        tip_base.push(mesh.nodes.len() as u32);
        mesh.nodes
            .push(Vector3::new(r_tip * phi.cos(), r_tip * phi.sin(), z));
    }
    for w in 0..(root_base.len() - 1) {
        let a = root_base[w];
        let b = tip_base[w];
        let c = tip_base[w + 1];
        let d = root_base[w + 1];
        block.connectivity.extend_from_slice(&[a, b, c]);
        block.connectivity.extend_from_slice(&[a, c, d]);
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Top-level dispatcher — picks the right `to_solid_*` from the
/// spec's [`GearKind`].
pub fn to_solid(spec: &GearSpec) -> Result<Solid, GearsError> {
    match spec.kind {
        GearKind::Spur => to_solid_spur(spec),
        GearKind::Helical => to_solid_helical(spec),
        GearKind::Bevel => to_solid_bevel(spec),
        GearKind::Worm => to_solid_worm(spec),
    }
}

fn extrude_twisted(
    profile2d: &[[f64; 2]],
    height: f64,
    twist_rad: f64,
    mesh_id: &str,
) -> Result<Solid, GearsError> {
    let n = profile2d.len();
    if n < 3 {
        return Err(GearsError::Degenerate(format!("profile has {n} vertices")));
    }
    let n_layers = 12;
    let mut mesh = Mesh::new(mesh_id);
    let mut block = ElementBlock::new(ElementType::Tri3);
    let mut ring_base: Vec<u32> = Vec::with_capacity(n_layers + 1);
    for layer in 0..=n_layers {
        let t = layer as f64 / n_layers as f64;
        let z = -height * 0.5 + height * t;
        let phi = twist_rad * t;
        let cs = phi.cos();
        let sn = phi.sin();
        let base = mesh.nodes.len() as u32;
        ring_base.push(base);
        for p in profile2d {
            mesh.nodes.push(Vector3::new(
                p[0] * cs - p[1] * sn,
                p[0] * sn + p[1] * cs,
                z,
            ));
        }
    }
    for layer in 0..n_layers {
        for i in 0..n {
            let j = (i + 1) % n;
            let a = ring_base[layer] + i as u32;
            let b = ring_base[layer] + j as u32;
            let c = ring_base[layer + 1] + j as u32;
            let d = ring_base[layer + 1] + i as u32;
            block.connectivity.extend_from_slice(&[a, b, c]);
            block.connectivity.extend_from_slice(&[a, c, d]);
        }
    }
    cap_polygon(&mut mesh, &mut block, ring_base[0], n, true);
    cap_polygon(&mut mesh, &mut block, *ring_base.last().unwrap(), n, false);
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

fn cap_polygon(mesh: &mut Mesh, block: &mut ElementBlock, base: u32, n: usize, flip: bool) {
    let mut c = Vector3::zeros();
    for i in 0..n {
        c += mesh.nodes[base as usize + i];
    }
    c /= n as f64;
    mesh.nodes.push(c);
    let centroid = mesh.nodes.len() as u32 - 1;
    for i in 0..n {
        let j = (i + 1) % n;
        if flip {
            block
                .connectivity
                .extend_from_slice(&[centroid, base + j as u32, base + i as u32]);
        } else {
            block
                .connectivity
                .extend_from_slice(&[centroid, base + i as u32, base + j as u32]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spur_to_solid_emits_mesh() {
        let mut spec = GearSpec::standard_spur(20);
        spec.face_width_mm = 8.0;
        let s = to_solid_spur(&spec).unwrap();
        match s {
            Solid::Mesh(m) => {
                assert!(!m.nodes.is_empty());
                assert!(m.total_elements() > 0);
            }
            _ => panic!("expected mesh-backed"),
        }
    }

    #[test]
    fn helical_requires_nonzero_helix_angle() {
        let mut spec = GearSpec::standard_spur(20);
        spec.kind = GearKind::Helical;
        spec.helix_angle_deg = 0.0;
        assert!(matches!(
            to_solid_helical(&spec),
            Err(GearsError::Degenerate(_))
        ));
        spec.helix_angle_deg = 20.0;
        assert!(to_solid_helical(&spec).is_ok());
    }

    #[test]
    fn bevel_to_solid_meshes() {
        let mut spec = GearSpec::standard_spur(20);
        spec.kind = GearKind::Bevel;
        spec.face_width_mm = 4.0;
        let s = to_solid_bevel(&spec).unwrap();
        assert!(matches!(s, Solid::Mesh(_)));
    }

    #[test]
    fn worm_to_solid_meshes() {
        let mut spec = GearSpec::standard_spur(20);
        spec.kind = GearKind::Worm;
        spec.face_width_mm = 30.0;
        let s = to_solid_worm(&spec).unwrap();
        assert!(matches!(s, Solid::Mesh(_)));
    }

    #[test]
    fn dispatcher_picks_right_emitter() {
        let mut spec = GearSpec::standard_spur(20);
        assert!(to_solid(&spec).is_ok());
        spec.kind = GearKind::Bevel;
        assert!(to_solid(&spec).is_ok());
    }
}
