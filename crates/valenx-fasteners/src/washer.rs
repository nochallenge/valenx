//! Washers — ISO 7089.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::FastenerError;

/// One washer entry in the ISO 7089 table.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WasherSpec {
    /// Nominal bolt diameter the washer fits (e.g. "M6").
    pub nominal: String,
    /// Inner diameter (mm) — should be slightly larger than the
    /// nominal bolt diameter.
    pub inner_diameter_mm: f64,
    /// Outer diameter (mm).
    pub outer_diameter_mm: f64,
    /// Thickness (mm).
    pub thickness_mm: f64,
}

/// Standard ISO 7089 washer table (Form A — flat washer).
pub fn iso7089_table() -> Vec<WasherSpec> {
    let entries: &[(f64, f64, f64, f64)] = &[
        // (nominal_d, id, od, thickness)
        (3.0, 3.2, 7.0, 0.5),
        (4.0, 4.3, 9.0, 0.8),
        (5.0, 5.3, 10.0, 1.0),
        (6.0, 6.4, 12.0, 1.6),
        (8.0, 8.4, 16.0, 1.6),
        (10.0, 10.5, 20.0, 2.0),
        (12.0, 13.0, 24.0, 2.5),
        (14.0, 15.0, 28.0, 2.5),
        (16.0, 17.0, 30.0, 3.0),
        (18.0, 19.0, 34.0, 3.0),
        (20.0, 21.0, 37.0, 3.0),
        (22.0, 23.0, 39.0, 3.0),
        (24.0, 25.0, 44.0, 4.0),
        (27.0, 28.0, 50.0, 4.0),
        (30.0, 31.0, 56.0, 4.0),
    ];
    entries
        .iter()
        .map(|&(d, id, od, t)| WasherSpec {
            nominal: format!("M{d:.0}"),
            inner_diameter_mm: id,
            outer_diameter_mm: od,
            thickness_mm: t,
        })
        .collect()
}

/// Generate a parametric solid for `spec` — annular cylinder.
/// Returns a [`Solid::Mesh`].
pub fn to_solid(spec: &WasherSpec) -> Result<Solid, FastenerError> {
    if spec.thickness_mm <= 0.0
        || spec.outer_diameter_mm <= spec.inner_diameter_mm
        || spec.inner_diameter_mm <= 0.0
    {
        return Err(FastenerError::BadParameter {
            name: "diameters|thickness",
            reason: format!(
                "must satisfy 0 < id < od, t > 0 (got id={} od={} t={})",
                spec.inner_diameter_mm, spec.outer_diameter_mm, spec.thickness_mm
            ),
        });
    }
    let mut mesh = Mesh::new(format!("washer_{}", spec.nominal));
    let mut block = ElementBlock::new(ElementType::Tri3);
    let n = 24;
    let r_in = spec.inner_diameter_mm * 0.5;
    let r_out = spec.outer_diameter_mm * 0.5;
    let h = spec.thickness_mm;

    let out_top = mesh.nodes.len() as u32;
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(r_out * theta.cos(), r_out * theta.sin(), h));
    }
    let out_bot = mesh.nodes.len() as u32;
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(r_out * theta.cos(), r_out * theta.sin(), 0.0));
    }
    let in_top = mesh.nodes.len() as u32;
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(r_in * theta.cos(), r_in * theta.sin(), h));
    }
    let in_bot = mesh.nodes.len() as u32;
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(r_in * theta.cos(), r_in * theta.sin(), 0.0));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        // Outer wall (outward facing).
        block.connectivity.extend_from_slice(&[
            out_top + i as u32,
            out_top + j as u32,
            out_bot + j as u32,
            out_top + i as u32,
            out_bot + j as u32,
            out_bot + i as u32,
        ]);
        // Inner wall (inward facing).
        block.connectivity.extend_from_slice(&[
            in_top + i as u32,
            in_bot + j as u32,
            in_top + j as u32,
            in_top + i as u32,
            in_bot + i as u32,
            in_bot + j as u32,
        ]);
        // Top annular ring.
        block.connectivity.extend_from_slice(&[
            out_top + i as u32,
            in_top + i as u32,
            in_top + j as u32,
            out_top + i as u32,
            in_top + j as u32,
            out_top + j as u32,
        ]);
        // Bottom annular ring.
        block.connectivity.extend_from_slice(&[
            out_bot + i as u32,
            in_bot + j as u32,
            in_bot + i as u32,
            out_bot + i as u32,
            out_bot + j as u32,
            in_bot + j as u32,
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
    fn washer_table_complete() {
        let t = iso7089_table();
        assert_eq!(t.len(), 15);
    }

    #[test]
    fn to_solid_m6_produces_mesh() {
        let w = iso7089_table()
            .into_iter()
            .find(|w| w.nominal == "M6")
            .unwrap();
        let s = to_solid(&w).unwrap();
        match s {
            Solid::Mesh(m) => assert!(!m.nodes.is_empty()),
            _ => panic!("expected mesh-backed"),
        }
    }

    #[test]
    fn bad_geometry_errors() {
        let bad = WasherSpec {
            nominal: "bad".into(),
            inner_diameter_mm: 10.0,
            outer_diameter_mm: 5.0, // id > od
            thickness_mm: 1.0,
        };
        assert!(matches!(
            to_solid(&bad),
            Err(FastenerError::BadParameter { .. })
        ));
    }
}
