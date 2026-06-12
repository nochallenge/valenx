//! Nuts — ISO 4032 hex.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::Solid;
use valenx_feature_tree::threads::ThreadSpec;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::FastenerError;

/// One nut entry in the ISO 4032 table.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NutSpec {
    /// Nominal diameter (e.g. "M6").
    pub nominal: String,
    /// Thread spec.
    pub thread: ThreadSpec,
    /// Thickness (m, ISO 4032: ≈0.8 d for hex nuts).
    pub thickness_mm: f64,
    /// Width across flats.
    pub width_across_flats_mm: f64,
}

/// Standard ISO 4032 hex nut table.
pub fn iso4032_hex_table() -> Vec<NutSpec> {
    let entries: &[(f64, f64, f64, f64)] = &[
        // (d, pitch, thickness, WAF)
        (3.0, 0.5, 2.4, 5.5),
        (4.0, 0.7, 3.2, 7.0),
        (5.0, 0.8, 4.0, 8.0),
        (6.0, 1.0, 5.0, 10.0),
        (8.0, 1.25, 6.5, 13.0),
        (10.0, 1.5, 8.0, 16.0),
        (12.0, 1.75, 10.0, 18.0),
        (14.0, 2.0, 11.0, 21.0),
        (16.0, 2.0, 13.0, 24.0),
        (18.0, 2.5, 15.0, 27.0),
        (20.0, 2.5, 16.0, 30.0),
        (22.0, 2.5, 18.0, 34.0),
        (24.0, 3.0, 19.0, 36.0),
        (27.0, 3.0, 22.0, 41.0),
        (30.0, 3.5, 24.0, 46.0),
    ];
    entries
        .iter()
        .map(|&(d, p, t, waf)| NutSpec {
            nominal: format!("M{d:.0}"),
            thread: ThreadSpec::new(
                valenx_feature_tree::threads::ThreadStandard::IsoMetric,
                format!("M{d:.0}"),
                d,
                p,
            ),
            thickness_mm: t,
            width_across_flats_mm: waf,
        })
        .collect()
}

/// Generate a parametric solid for `spec`. Returns a
/// [`Solid::Mesh`] — hexagonal prism with a cylindrical bore at the
/// thread's tap-drill diameter.
pub fn to_solid(spec: &NutSpec) -> Result<Solid, FastenerError> {
    if spec.thickness_mm <= 0.0 || spec.width_across_flats_mm <= 0.0 {
        return Err(FastenerError::BadParameter {
            name: "thickness|waf",
            reason: format!(
                "must be > 0, got t={} waf={}",
                spec.thickness_mm, spec.width_across_flats_mm
            ),
        });
    }
    let mut mesh = Mesh::new(format!("nut_{}", spec.nominal));
    let mut block = ElementBlock::new(ElementType::Tri3);

    let r_outer = spec.width_across_flats_mm / (2.0 * (std::f64::consts::PI / 6.0).cos());
    let r_bore = spec.thread.tap_drill_diameter() * 0.5;
    let h = spec.thickness_mm;

    let n_outer = 6;
    let n_bore = 16;

    // Outer hex top + bottom rings.
    let outer_top = mesh.nodes.len() as u32;
    for i in 0..n_outer {
        let theta =
            (i as f64 / n_outer as f64) * std::f64::consts::TAU + std::f64::consts::FRAC_PI_6;
        mesh.nodes.push(Vector3::new(
            r_outer * theta.cos(),
            r_outer * theta.sin(),
            h,
        ));
    }
    let outer_bot = mesh.nodes.len() as u32;
    for i in 0..n_outer {
        let theta =
            (i as f64 / n_outer as f64) * std::f64::consts::TAU + std::f64::consts::FRAC_PI_6;
        mesh.nodes.push(Vector3::new(
            r_outer * theta.cos(),
            r_outer * theta.sin(),
            0.0,
        ));
    }
    // Inner bore top + bottom rings.
    let bore_top = mesh.nodes.len() as u32;
    for i in 0..n_bore {
        let theta = (i as f64 / n_bore as f64) * std::f64::consts::TAU;
        mesh.nodes
            .push(Vector3::new(r_bore * theta.cos(), r_bore * theta.sin(), h));
    }
    let bore_bot = mesh.nodes.len() as u32;
    for i in 0..n_bore {
        let theta = (i as f64 / n_bore as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(
            r_bore * theta.cos(),
            r_bore * theta.sin(),
            0.0,
        ));
    }
    // Outer hex side walls.
    for i in 0..n_outer {
        let j = (i + 1) % n_outer;
        let a = outer_top + i as u32;
        let b = outer_top + j as u32;
        let c = outer_bot + j as u32;
        let d = outer_bot + i as u32;
        block.connectivity.extend_from_slice(&[a, b, c, a, c, d]);
    }
    // Bore side walls (inward-facing).
    for i in 0..n_bore {
        let j = (i + 1) % n_bore;
        let a = bore_top + i as u32;
        let b = bore_top + j as u32;
        let c = bore_bot + j as u32;
        let d = bore_bot + i as u32;
        block.connectivity.extend_from_slice(&[a, c, b, a, d, c]);
    }
    // Top + bottom annular rings: fan from each hex vertex to the
    // nearest bore vertex pair. Cheap approximation — gives a
    // visually closed annulus.
    for i in 0..n_outer {
        let j = (i + 1) % n_outer;
        let bore_i = (i * n_bore / n_outer) as u32;
        let bore_j = (j * n_bore / n_outer) as u32;
        // Top face (CCW seen from +Z).
        block.connectivity.extend_from_slice(&[
            outer_top + i as u32,
            outer_top + j as u32,
            bore_top + bore_j,
        ]);
        block.connectivity.extend_from_slice(&[
            outer_top + i as u32,
            bore_top + bore_j,
            bore_top + bore_i,
        ]);
        // Bottom face (CW seen from +Z so it faces -Z).
        block.connectivity.extend_from_slice(&[
            outer_bot + j as u32,
            outer_bot + i as u32,
            bore_bot + bore_i,
        ]);
        block.connectivity.extend_from_slice(&[
            outer_bot + j as u32,
            bore_bot + bore_i,
            bore_bot + bore_j,
        ]);
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iso_nut_table_complete() {
        let t = iso4032_hex_table();
        assert_eq!(t.len(), 15);
        let m6 = t.iter().find(|n| n.nominal == "M6").unwrap();
        assert!((m6.thickness_mm - 5.0).abs() < 1e-9);
        assert!((m6.width_across_flats_mm - 10.0).abs() < 1e-9);
    }

    #[test]
    fn to_solid_m6_produces_mesh() {
        let nut = iso4032_hex_table()
            .into_iter()
            .find(|n| n.nominal == "M6")
            .unwrap();
        let s = to_solid(&nut).unwrap();
        match s {
            Solid::Mesh(m) => assert!(!m.nodes.is_empty()),
            _ => panic!("expected mesh-backed"),
        }
    }

    #[test]
    fn zero_thickness_errors() {
        let mut nut = iso4032_hex_table().into_iter().next().unwrap();
        nut.thickness_mm = 0.0;
        assert!(matches!(
            to_solid(&nut),
            Err(FastenerError::BadParameter { .. })
        ));
    }
}
