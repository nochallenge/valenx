//! Bolts — ISO 4017 hex + ANSI B18.2.1 hex specs + to_solid.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::Solid;
use valenx_feature_tree::threads::ThreadSpec;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::FastenerError;

/// Bolt head + thread style.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoltKind {
    /// Hex head (ISO 4017 metric, ANSI B18.2.1 imperial).
    HexHead,
    /// Internal hex (socket head cap screw — ISO 4762).
    SocketHead,
    /// Flat / countersunk head.
    FlatHead,
    /// Button head.
    ButtonHead,
    /// Hex flange.
    HexFlange,
}

impl BoltKind {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::HexHead => "Hex head",
            Self::SocketHead => "Socket head",
            Self::FlatHead => "Flat head",
            Self::ButtonHead => "Button head",
            Self::HexFlange => "Hex flange",
        }
    }
}

/// One bolt entry in a standard parts table.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BoltSpec {
    /// Head + thread style.
    pub kind: BoltKind,
    /// Nominal diameter (e.g. "M6", "1/4-20").
    pub nominal: String,
    /// Total length under the head (mm).
    pub length_mm: f64,
    /// Thread spec from Phase 13.
    pub thread: ThreadSpec,
}

impl BoltSpec {
    /// Width across flats — derived from the thread spec's nominal
    /// diameter using the ISO 4017 convention `WAF ≈ 1.5 d`.
    pub fn width_across_flats_mm(&self) -> f64 {
        1.5 * self.thread.nominal_diameter
    }

    /// Head height — `0.7 d` (ISO 4017 default).
    pub fn head_height_mm(&self) -> f64 {
        0.7 * self.thread.nominal_diameter
    }

    /// ISO 898-1 tensile stress area (mm²) — the effective thread cross-section at the stress
    /// diameter, `Aₛ = (π/4)·(d − 0.9382·P)²`, for nominal diameter `d` and pitch `P` (mm). Feeds
    /// the bolt tensile capacity. M6×1.0 → ≈ 20.1 mm².
    pub fn tensile_stress_area_mm2(&self) -> f64 {
        let d = self.thread.nominal_diameter;
        let p = self.thread.pitch;
        std::f64::consts::PI / 4.0 * (d - 0.9382 * p).powi(2)
    }

    /// ISO 261 pitch (effective) diameter (mm) — `d₂ = d − 0.6495·P` (the `0.6495 = 3√3/8`
    /// 60°-thread constant), the diameter where thread thickness equals groove width. Distinct
    /// from the major (nominal) diameter and the tensile-stress area. M6×1.0 → ≈ 5.351 mm.
    pub fn pitch_diameter_mm(&self) -> f64 {
        let d = self.thread.nominal_diameter;
        let p = self.thread.pitch;
        d - 0.6495 * p
    }
}

/// Standard ISO 4017 hex bolt table (M3 through M30 with a common
/// length per size — the call-site can override `length_mm`).
pub fn iso4017_hex_table() -> Vec<BoltSpec> {
    let entries: &[(f64, f64, f64)] = &[
        // (nominal_diameter, pitch, default_length)
        (3.0, 0.5, 16.0),
        (4.0, 0.7, 20.0),
        (5.0, 0.8, 25.0),
        (6.0, 1.0, 30.0),
        (8.0, 1.25, 40.0),
        (10.0, 1.5, 50.0),
        (12.0, 1.75, 60.0),
        (14.0, 2.0, 70.0),
        (16.0, 2.0, 80.0),
        (18.0, 2.5, 90.0),
        (20.0, 2.5, 100.0),
        (22.0, 2.5, 110.0),
        (24.0, 3.0, 120.0),
        (27.0, 3.0, 130.0),
        (30.0, 3.5, 140.0),
    ];
    entries
        .iter()
        .map(|&(d, p, l)| BoltSpec {
            kind: BoltKind::HexHead,
            nominal: format!("M{d:.0}"),
            length_mm: l,
            thread: ThreadSpec::new(
                valenx_feature_tree::threads::ThreadStandard::IsoMetric,
                format!("M{d:.0}"),
                d,
                p,
            ),
        })
        .collect()
}

/// Standard ANSI B18.2.1 hex bolt table (1/4 through 1-1/2 UNC, with
/// default thread pitches).
pub fn ansi_b18_2_1_table() -> Vec<BoltSpec> {
    let entries: &[(&str, f64, u32, f64)] = &[
        // (nominal label, diameter_mm, tpi, default_length_mm)
        ("1/4-20", 6.35, 20, 25.4),
        ("5/16-18", 7.94, 18, 31.75),
        ("3/8-16", 9.525, 16, 38.1),
        ("7/16-14", 11.11, 14, 44.45),
        ("1/2-13", 12.7, 13, 50.8),
        ("5/8-11", 15.875, 11, 63.5),
        ("3/4-10", 19.05, 10, 76.2),
        ("7/8-9", 22.225, 9, 88.9),
        ("1-8", 25.4, 8, 101.6),
        ("1-1/8-7", 28.575, 7, 114.3),
        ("1-1/4-7", 31.75, 7, 127.0),
        ("1-3/8-6", 34.925, 6, 139.7),
        ("1-1/2-6", 38.1, 6, 152.4),
    ];
    entries
        .iter()
        .map(|&(label, d, tpi, l)| {
            let pitch = 25.4 / tpi as f64;
            BoltSpec {
                kind: BoltKind::HexHead,
                nominal: label.to_string(),
                length_mm: l,
                thread: ThreadSpec::new(
                    valenx_feature_tree::threads::ThreadStandard::UnifiedNational,
                    label,
                    d,
                    pitch,
                ),
            }
        })
        .collect()
}

/// Generate a parametric solid for `spec`. Returns a
/// [`Solid::Mesh`] — a hexagonal-prism head + cylindrical shank.
/// True threaded geometry is not generated (the thread is
/// represented by `spec.thread`).
pub fn to_solid(spec: &BoltSpec) -> Result<Solid, FastenerError> {
    if spec.length_mm <= 0.0 || spec.thread.nominal_diameter <= 0.0 {
        return Err(FastenerError::BadParameter {
            name: "length|nominal_diameter",
            reason: format!(
                "must be > 0, got L={} d={}",
                spec.length_mm, spec.thread.nominal_diameter
            ),
        });
    }
    let mut mesh = Mesh::new(format!("bolt_{}", spec.nominal));
    let mut block = ElementBlock::new(ElementType::Tri3);

    let r_shank = spec.thread.nominal_diameter * 0.5;
    let waf = spec.width_across_flats_mm();
    let r_head = waf / (2.0 * (std::f64::consts::PI / 6.0).cos()); // circumscribed
    let head_h = spec.head_height_mm();

    // Hex head: 6 corner verts top + bottom.
    let n = 6;
    let head_top_base = mesh.nodes.len() as u32;
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU + std::f64::consts::FRAC_PI_6;
        mesh.nodes.push(Vector3::new(
            r_head * theta.cos(),
            r_head * theta.sin(),
            head_h,
        ));
    }
    let head_bot_base = mesh.nodes.len() as u32;
    for i in 0..n {
        let theta = (i as f64 / n as f64) * std::f64::consts::TAU + std::f64::consts::FRAC_PI_6;
        mesh.nodes.push(Vector3::new(
            r_head * theta.cos(),
            r_head * theta.sin(),
            0.0,
        ));
    }
    // Top + bottom face fans.
    let top_center = mesh.nodes.len() as u32;
    mesh.nodes.push(Vector3::new(0.0, 0.0, head_h));
    let bot_center = mesh.nodes.len() as u32;
    mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
    for i in 0..n {
        let j = (i + 1) % n;
        block.connectivity.extend_from_slice(&[
            top_center,
            head_top_base + i as u32,
            head_top_base + j as u32,
        ]);
        block.connectivity.extend_from_slice(&[
            bot_center,
            head_bot_base + j as u32,
            head_bot_base + i as u32,
        ]);
    }
    // 6 side walls.
    for i in 0..n {
        let j = (i + 1) % n;
        let a = head_top_base + i as u32;
        let b = head_top_base + j as u32;
        let c = head_bot_base + j as u32;
        let d = head_bot_base + i as u32;
        block.connectivity.extend_from_slice(&[a, b, c]);
        block.connectivity.extend_from_slice(&[a, c, d]);
    }

    // Shank cylinder — 16-sided.
    let shank_segs = 16;
    let shank_top_base = mesh.nodes.len() as u32;
    for i in 0..shank_segs {
        let theta = (i as f64 / shank_segs as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(
            r_shank * theta.cos(),
            r_shank * theta.sin(),
            0.0,
        ));
    }
    let shank_bot_base = mesh.nodes.len() as u32;
    for i in 0..shank_segs {
        let theta = (i as f64 / shank_segs as f64) * std::f64::consts::TAU;
        mesh.nodes.push(Vector3::new(
            r_shank * theta.cos(),
            r_shank * theta.sin(),
            -spec.length_mm,
        ));
    }
    let shank_bot_center = mesh.nodes.len() as u32;
    mesh.nodes.push(Vector3::new(0.0, 0.0, -spec.length_mm));
    for i in 0..shank_segs {
        let j = (i + 1) % shank_segs;
        let a = shank_top_base + i as u32;
        let b = shank_top_base + j as u32;
        let c = shank_bot_base + j as u32;
        let d = shank_bot_base + i as u32;
        block.connectivity.extend_from_slice(&[a, b, c]);
        block.connectivity.extend_from_slice(&[a, c, d]);
        // Shank bottom cap.
        block.connectivity.extend_from_slice(&[
            shank_bot_center,
            shank_bot_base + j as u32,
            shank_bot_base + i as u32,
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
    fn iso_table_covers_m3_through_m30() {
        let t = iso4017_hex_table();
        assert_eq!(t.len(), 15);
        assert_eq!(t.first().unwrap().nominal, "M3");
        assert_eq!(t.last().unwrap().nominal, "M30");
    }

    #[test]
    fn ansi_table_covers_quarter_through_1_5() {
        let t = ansi_b18_2_1_table();
        assert_eq!(t.len(), 13);
        assert_eq!(t.first().unwrap().nominal, "1/4-20");
        assert_eq!(t.last().unwrap().nominal, "1-1/2-6");
    }

    #[test]
    fn to_solid_m6_produces_mesh_backed() {
        let m6 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .unwrap();
        let s = to_solid(&m6).unwrap();
        match s {
            Solid::Mesh(m) => {
                assert!(!m.nodes.is_empty());
                assert!(m.total_elements() > 0);
            }
            _ => panic!("expected mesh-backed solid"),
        }
    }

    #[test]
    fn to_solid_zero_length_errors() {
        let mut m6 = iso4017_hex_table().into_iter().next().unwrap();
        m6.length_mm = 0.0;
        assert!(matches!(
            to_solid(&m6),
            Err(FastenerError::BadParameter { .. })
        ));
    }

    #[test]
    fn tensile_stress_area_m6() {
        // M6×1.0: (π/4)·(6 − 0.9382)² = (π/4)·5.0618² ≈ 20.12 mm² (ISO 898-1).
        let m6 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .unwrap();
        let a_s = m6.tensile_stress_area_mm2();
        assert!((a_s - 20.12).abs() < 0.05);
        // Strictly less than the major-diameter circular area (root < nominal).
        assert!(a_s < std::f64::consts::PI / 4.0 * 6.0_f64.powi(2));
    }

    #[test]
    fn root_minor_diameter_m6_iso724() {
        // ISO 724 minor (root) diameter d₃ = d − 1.2269·P. For M6×1.0:
        //   d₃ = 6 − 1.2269·1.0 = 4.7731 mm.
        // The exposed `ThreadSpec::minor_diameter()` (= d − 2·depth with
        // depth = 0.61343·P, i.e. d − 1.22686·P) realises this; pin it to the
        // published closed form. The 1.22686 vs 1.2269 constant differ by
        // 6e-5·P → use tol 1e-3 mm (covers the rounding, still far tighter
        // than the ~1.23 mm the root sits below the nominal).
        let m6 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .unwrap();
        let d3 = m6.thread.minor_diameter();
        let expected = 6.0 - 1.2269 * 1.0; // ISO 724 → 4.7731 mm
        assert!(
            (d3 - expected).abs() < 1e-3,
            "minor dia {d3} mm vs ISO 724 {expected} mm"
        );
        assert!((d3 - 4.7731).abs() < 1e-3, "minor dia {d3} mm, expected 4.7731");
        // Ordering sanity: root < pitch < major (nominal).
        assert!(d3 < m6.pitch_diameter_mm());
        assert!(m6.pitch_diameter_mm() < m6.thread.nominal_diameter);
    }

    #[test]
    fn pitch_diameter_m6_and_m8() {
        // M6×1.0: 6 − 0.6495·1 = 5.3505 mm; M8×1.25: 8 − 0.6495·1.25 = 7.188125 mm.
        let m6 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .unwrap();
        assert!((m6.pitch_diameter_mm() - 5.3505).abs() < 1e-4);
        let m8 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M8")
            .unwrap();
        assert!((m8.pitch_diameter_mm() - 7.188125).abs() < 1e-4);
        // The pitch diameter is strictly less than the nominal (major) diameter.
        assert!(m6.pitch_diameter_mm() < m6.thread.nominal_diameter);
    }
}
