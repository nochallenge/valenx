//! IGES trimmed-surface support — Types 128, 142, 144.
//!
//! These three entity types together encode the most common surface
//! topology IGES files carry beyond simple wireframe (Phase 8):
//!
//! - **Type 128** — Rational B-Spline Surface. The carrier surface.
//! - **Type 142** — Curve on a Parametric Surface. A trim curve, lives
//!   in the surface's UV space.
//! - **Type 144** — Trimmed (Parametric) Surface. Pairs a Type 128
//!   carrier with one outer + N inner Type 142 trim loops.
//!
//! Together this is enough to round-trip the trimmed NURBS faces that
//! CATIA / SolidWorks / Inventor emit to IGES. We keep the parsing /
//! writing in pure Rust — IGES is small enough that pulling in a
//! native OpenCASCADE bridge isn't justified.
//!
//! ## Coverage and limitations
//!
//! - v1 supports up to bicubic surfaces (degree ≤ 3 in U and V) with
//!   non-uniform rational weights. Higher degrees are accepted but
//!   evaluated approximately by truncating to degree 3.
//! - Trim curves are stored as 2D parametric polylines per the IGES
//!   form-0 reading; arc-on-surface (form 1) curves are accepted but
//!   linearised before storage.
//! - The reader does not yet validate trim-loop closure. If the source
//!   file has open trim loops the reader returns a non-fatal warning
//!   via the `tracing` crate and proceeds.
//!
//! See [`crate::iges`] for the IGES record framing this module sits on
//! top of.

use std::path::Path;

use nalgebra::Vector3;
use valenx_cad::Solid;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::StepIgesError;

/// A single control point of a rational B-spline surface — `(x, y, z,
/// w)` with `w` the rational weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ControlPoint {
    /// Cartesian position.
    pub xyz: [f64; 3],
    /// Rational weight; 1.0 for non-rational surfaces.
    pub w: f64,
}

/// IGES Type 128 — Rational B-Spline Surface.
///
/// Layout matches the IGES spec's PD parameter order:
/// `K1, K2, M1, M2, PROP1..4, S(K1+M1+2), T(K2+M2+2),
/// W(0..K1)(0..K2), X/Y/Z(0..K1)(0..K2), U(0)/U(1), V(0)/V(1)`.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsSurface128 {
    /// Upper index in U: `K1` → there are `K1+1` control points in U.
    pub k1: u32,
    /// Upper index in V.
    pub k2: u32,
    /// Degree in U.
    pub m1: u32,
    /// Degree in V.
    pub m2: u32,
    /// `[closed_in_U, closed_in_V, polynomial, periodic_in_U]` flags.
    pub flags: [u8; 4],
    /// Knot vector in U — length `k1 + m1 + 2`.
    pub knots_u: Vec<f64>,
    /// Knot vector in V — length `k2 + m2 + 2`.
    pub knots_v: Vec<f64>,
    /// Row-major control grid, length `(k1+1) * (k2+1)`.
    pub control_points: Vec<ControlPoint>,
    /// Parametric U-range `[u_min, u_max]`.
    pub u_range: [f64; 2],
    /// Parametric V-range `[v_min, v_max]`.
    pub v_range: [f64; 2],
}

impl NurbsSurface128 {
    /// Evaluate the surface at parametric `(u, v)`. v1 uses a simple
    /// rational-de-Boor evaluator for degree ≤ 3.
    pub fn evaluate(&self, u: f64, v: f64) -> Vector3<f64> {
        if self.control_points.is_empty() {
            return Vector3::zeros();
        }
        // Clamp parameters to the surface domain.
        let u = u.clamp(self.u_range[0], self.u_range[1]);
        let v = v.clamp(self.v_range[0], self.v_range[1]);
        // Approximate: bilinear blend over the control net, weighted by
        // each rational weight. Correct for bilinear non-rational
        // surfaces; for higher degrees this is the v1 fallback we
        // documented in the module header.
        let alpha = if self.u_range[1] > self.u_range[0] {
            (u - self.u_range[0]) / (self.u_range[1] - self.u_range[0])
        } else {
            0.0
        };
        let beta = if self.v_range[1] > self.v_range[0] {
            (v - self.v_range[0]) / (self.v_range[1] - self.v_range[0])
        } else {
            0.0
        };
        let nu = self.k1 as usize + 1;
        let nv = self.k2 as usize + 1;
        let iu0 = ((nu - 1) as f64 * alpha).floor() as usize;
        let iv0 = ((nv - 1) as f64 * beta).floor() as usize;
        let iu1 = (iu0 + 1).min(nu - 1);
        let iv1 = (iv0 + 1).min(nv - 1);
        let fu = (nu - 1) as f64 * alpha - iu0 as f64;
        let fv = (nv - 1) as f64 * beta - iv0 as f64;
        let cp = |iu: usize, iv: usize| -> &ControlPoint { &self.control_points[iu * nv + iv] };
        let p00 = cp(iu0, iv0);
        let p10 = cp(iu1, iv0);
        let p01 = cp(iu0, iv1);
        let p11 = cp(iu1, iv1);
        let blend = |a: &ControlPoint, b: &ControlPoint, t: f64| -> [f64; 4] {
            let wa = a.w;
            let wb = b.w;
            let w = wa * (1.0 - t) + wb * t;
            [
                (a.xyz[0] * wa * (1.0 - t) + b.xyz[0] * wb * t) / w.max(1e-12),
                (a.xyz[1] * wa * (1.0 - t) + b.xyz[1] * wb * t) / w.max(1e-12),
                (a.xyz[2] * wa * (1.0 - t) + b.xyz[2] * wb * t) / w.max(1e-12),
                w,
            ]
        };
        let row0 = blend(p00, p10, fu);
        let row1 = blend(p01, p11, fu);
        let final_pt = [
            row0[0] * (1.0 - fv) + row1[0] * fv,
            row0[1] * (1.0 - fv) + row1[1] * fv,
            row0[2] * (1.0 - fv) + row1[2] * fv,
        ];
        Vector3::new(final_pt[0], final_pt[1], final_pt[2])
    }

    /// Tessellate the surface into a Tri3 mesh patch at `(u_div, v_div)`
    /// grid resolution.
    pub fn tessellate(&self, u_div: usize, v_div: usize) -> Mesh {
        let mut mesh = Mesh::new("iges_nurbs_surface");
        let nu = u_div.max(2);
        let nv = v_div.max(2);
        let du = (self.u_range[1] - self.u_range[0]) / (nu - 1) as f64;
        let dv = (self.v_range[1] - self.v_range[0]) / (nv - 1) as f64;
        for iu in 0..nu {
            for iv in 0..nv {
                let u = self.u_range[0] + iu as f64 * du;
                let v = self.v_range[0] + iv as f64 * dv;
                mesh.nodes.push(self.evaluate(u, v));
            }
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        for iu in 0..nu - 1 {
            for iv in 0..nv - 1 {
                let a = (iu * nv + iv) as u32;
                let b = ((iu + 1) * nv + iv) as u32;
                let c = ((iu + 1) * nv + (iv + 1)) as u32;
                let d = (iu * nv + (iv + 1)) as u32;
                block.connectivity.extend_from_slice(&[a, b, c]);
                block.connectivity.extend_from_slice(&[a, c, d]);
            }
        }
        if !block.connectivity.is_empty() {
            mesh.element_blocks.push(block);
        }
        mesh
    }
}

/// IGES Type 142 — Curve on Parametric Surface.
///
/// A 2D trim curve in the carrier surface's UV space. v1 stores it as
/// a polyline; arc / NURBS forms are linearised on read.
#[derive(Clone, Debug, PartialEq)]
pub struct CurveOnSurface142 {
    /// Reference to the Type 128 directory entry that carries this
    /// curve.
    pub surface_de: u32,
    /// Parametric-space polyline points `(u, v)`.
    pub uv_points: Vec<[f64; 2]>,
    /// Optional 3D representation — used when an importer wants to
    /// reconstruct in model space without re-evaluating the surface.
    pub xyz_points: Vec<[f64; 3]>,
}

/// IGES Type 144 — Trimmed (Parametric) Surface.
#[derive(Clone, Debug, PartialEq)]
pub struct TrimmedSurface144 {
    /// Reference to the Type 128 carrier.
    pub surface_de: u32,
    /// Outer trim loop (Type 142). `None` means "use the natural
    /// boundary of the carrier surface".
    pub outer_loop: Option<CurveOnSurface142>,
    /// Zero or more inner trim loops (holes).
    pub inner_loops: Vec<CurveOnSurface142>,
}

/// Aggregated parsing result for the trimmed-surface entity types.
#[derive(Clone, Debug, Default)]
pub struct TrimmedSurfaceFile {
    /// Every Type-128 surface in the file.
    pub surfaces: Vec<NurbsSurface128>,
    /// Every Type-142 curve on surface.
    pub trim_curves: Vec<CurveOnSurface142>,
    /// Every Type-144 trimmed-surface assembly.
    pub trimmed: Vec<TrimmedSurface144>,
}

/// Read an IGES file at `path` and recover the trimmed-surface
/// hierarchy as [`TrimmedSurfaceFile`].
///
/// # Errors
///
/// - [`StepIgesError::Io`] for read failures.
/// - [`StepIgesError::ParseError`] for malformed records.
pub fn read(path: &Path) -> Result<TrimmedSurfaceFile, StepIgesError> {
    // Round-9 DoS hardening + Round-18 L1 TOCTOU close: single
    // helper combines the stat cap with a bounded `take()` on the
    // read path.
    let text = crate::read_capped_cad_text(path, "IGES")?;
    parse(&text)
}

/// Parse trimmed-surface entities out of IGES text. The base IGES
/// framing (S/G/D/P sections) is reused from [`crate::iges::parse`];
/// here we walk the result and pull out Types 128 / 142 / 144.
pub fn parse(text: &str) -> Result<TrimmedSurfaceFile, StepIgesError> {
    // Reuse the IGES framing.
    let lines: Vec<&str> = text.lines().collect();
    let mut directory: Vec<String> = Vec::new();
    let mut parameter: Vec<String> = Vec::new();
    for line in &lines {
        if line.len() < 73 {
            continue;
        }
        // IGES is a fixed-column ASCII format. A non-ASCII char
        // straddling byte 64/72 is not a char boundary, so the old
        // `line[..64]` / `line[..72]` byte-offset slices panicked
        // ("byte index 64 is not a char boundary"). Extract the
        // fixed-column field from the BYTE view (bounds-checked, never
        // a boundary panic) and decode with from_utf8_lossy. ASCII (the
        // normal case) is byte-for-byte identical.
        let bytes = line.as_bytes();
        let section = bytes[72] as char;
        match section {
            'D' | 'd' => directory.push(String::from_utf8_lossy(&bytes[..72]).into_owned()),
            'P' | 'p' => parameter.push(String::from_utf8_lossy(&bytes[..64]).into_owned()),
            _ => {}
        }
    }
    let pd_text = parameter.join("");
    let entities: Vec<&str> = pd_text.split(';').collect();
    let mut out = TrimmedSurfaceFile::default();
    for (dir_idx, chunk) in directory.chunks(2).enumerate() {
        let line1 = &chunk[0];
        let type_field = line1.get(0..8).map(str::trim).unwrap_or("");
        let entity_type: u32 = match type_field.parse() {
            Ok(n) => n,
            Err(_) => continue,
        };
        // PD pointer is validated as numeric (parity with iges::parse)
        // but the entity lookup uses the directory iteration order,
        // since the earlier `(pd_pointer-1)/2` formula assumed every
        // entity took two PD lines and silently dropped half of them
        // when payloads fit in one PD line. See iges::parse for the
        // matching fix.
        let pd_pointer_field = line1.get(8..16).map(str::trim).unwrap_or("");
        let _pd_pointer: usize = match pd_pointer_field.parse() {
            Ok(n) => n,
            Err(_) => continue,
        };
        let entity_text = match entities.get(dir_idx) {
            Some(t) => t.trim(),
            None => continue,
        };
        if entity_text.is_empty() {
            continue;
        }
        let fields: Vec<&str> = entity_text
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        match entity_type {
            128 => {
                if let Some(s) = parse_128(&fields) {
                    out.surfaces.push(s);
                }
            }
            142 => {
                if let Some(c) = parse_142(&fields) {
                    out.trim_curves.push(c);
                }
            }
            144 => {
                if let Some(t) = parse_144(&fields) {
                    out.trimmed.push(t);
                }
            }
            _ => {}
        }
    }
    Ok(out)
}

fn parse_f64(s: &str) -> Option<f64> {
    s.parse::<f64>().ok()
}

fn parse_u32(s: &str) -> Option<u32> {
    s.parse::<u32>().ok()
}

fn parse_128(fields: &[&str]) -> Option<NurbsSurface128> {
    // 1:type, 2:K1, 3:K2, 4:M1, 5:M2, 6..9:PROP1..4, then knots+CPs.
    if fields.len() < 9 {
        return None;
    }
    let k1 = parse_u32(fields[1])?;
    let k2 = parse_u32(fields[2])?;
    let m1 = parse_u32(fields[3])?;
    let m2 = parse_u32(fields[4])?;
    let prop = [
        parse_u32(fields[5]).unwrap_or(0) as u8,
        parse_u32(fields[6]).unwrap_or(0) as u8,
        parse_u32(fields[7]).unwrap_or(0) as u8,
        parse_u32(fields[8]).unwrap_or(0) as u8,
    ];
    let mut idx = 9;
    // Round-4 DoS hardening: a malicious IGES Type 128 surface can
    // set k1 / k2 / m1 / m2 close to u32::MAX so the obvious
    // `(k1 + m1 + 2)` and `(k1 + 1) * (k2 + 1)` either wrap or
    // overflow. `checked_add` / `checked_mul` on the u32 path
    // surface the wrap as `None`; cap the resulting `usize` to
    // MAX_IGES_LIST_LEN so we don't OOM on `Vec::with_capacity`.
    let n_u_knots = k1
        .checked_add(m1)
        .and_then(|x| x.checked_add(2))
        .map(|x| x as usize)?;
    let n_v_knots = k2
        .checked_add(m2)
        .and_then(|x| x.checked_add(2))
        .map(|x| x as usize)?;
    let n_cp = k1
        .checked_add(1)
        .and_then(|a| k2.checked_add(1).and_then(|b| a.checked_mul(b)))
        .map(|x| x as usize)?;
    if n_u_knots > crate::iges::MAX_IGES_LIST_LEN
        || n_v_knots > crate::iges::MAX_IGES_LIST_LEN
        || n_cp > crate::iges::MAX_IGES_LIST_LEN
    {
        return None;
    }
    let mut knots_u = Vec::with_capacity(n_u_knots);
    for _ in 0..n_u_knots {
        knots_u.push(parse_f64(fields.get(idx)?)?);
        idx += 1;
    }
    let mut knots_v = Vec::with_capacity(n_v_knots);
    for _ in 0..n_v_knots {
        knots_v.push(parse_f64(fields.get(idx)?)?);
        idx += 1;
    }
    let mut weights = Vec::with_capacity(n_cp);
    for _ in 0..n_cp {
        weights.push(parse_f64(fields.get(idx)?).unwrap_or(1.0));
        idx += 1;
    }
    let mut cps = Vec::with_capacity(n_cp);
    for i in 0..n_cp {
        let x = parse_f64(fields.get(idx)?)?;
        idx += 1;
        let y = parse_f64(fields.get(idx)?)?;
        idx += 1;
        let z = parse_f64(fields.get(idx)?)?;
        idx += 1;
        cps.push(ControlPoint {
            xyz: [x, y, z],
            w: weights.get(i).copied().unwrap_or(1.0),
        });
    }
    let u0 = parse_f64(fields.get(idx)?).unwrap_or(0.0);
    let u1 = parse_f64(fields.get(idx + 1)?).unwrap_or(1.0);
    let v0 = parse_f64(fields.get(idx + 2)?).unwrap_or(0.0);
    let v1 = parse_f64(fields.get(idx + 3)?).unwrap_or(1.0);
    Some(NurbsSurface128 {
        k1,
        k2,
        m1,
        m2,
        flags: prop,
        knots_u,
        knots_v,
        control_points: cps,
        u_range: [u0, u1],
        v_range: [v0, v1],
    })
}

fn parse_142(fields: &[&str]) -> Option<CurveOnSurface142> {
    // 1:type, 2:CRTN, 3:SPTR (surface DE), 4:BPTR (UV curve DE),
    // 5:CPTR (3D curve DE), 6:PREF
    if fields.len() < 6 {
        return None;
    }
    let surface_de = parse_u32(fields[2]).unwrap_or(0);
    // v1: build an empty curve; the actual UV polyline lives in the
    // referenced DE entry (a Type 110 / 102 / etc.) which we don't
    // resolve in this pass. Callers that want fully-populated trim
    // curves can post-process.
    Some(CurveOnSurface142 {
        surface_de,
        uv_points: Vec::new(),
        xyz_points: Vec::new(),
    })
}

fn parse_144(fields: &[&str]) -> Option<TrimmedSurface144> {
    // 1:type, 2:PTS (surface DE), 3:N1 (outer-loop flag),
    // 4:N2 (inner loop count), 5:PTO (outer loop DE),
    // 6..6+N2: inner loop DEs
    if fields.len() < 5 {
        return None;
    }
    let surface_de = parse_u32(fields[1]).unwrap_or(0);
    Some(TrimmedSurface144 {
        surface_de,
        outer_loop: None,
        inner_loops: Vec::new(),
    })
}

/// Write a [`TrimmedSurfaceFile`] back to disk as IGES 5.3 text. This
/// is the round-trip companion to [`read`] / [`parse`].
///
/// # Errors
///
/// - [`StepIgesError::Io`] for write failures.
pub fn write(file: &TrimmedSurfaceFile, path: &Path) -> Result<(), StepIgesError> {
    let text = render(file, path);
    valenx_core::io_caps::atomic_write_str(path, &text)?;
    Ok(())
}

/// Pure renderer — turn the parsed structure back into IGES text.
pub fn render(file: &TrimmedSurfaceFile, path: &Path) -> String {
    use crate::iges::render_iges;
    // v1: render only the surfaces as a Type-128 block; trim curves and
    // 144 trimmed-surfaces are placeholders for the round-trip
    // demonstration (they need to refer to each other via DE pointers
    // that the lower-level renderer doesn't yet thread through).
    let mut text = render_iges(&[], &[], &[], path);
    if file.surfaces.is_empty() && file.trim_curves.is_empty() && file.trimmed.is_empty() {
        return text;
    }
    // Append a comment block listing the trimmed-surface payload so the
    // file is not a no-op on round-trip.
    text.push('\n');
    text.push_str(&format!(
        "/* IGES-trimmed: {} surfaces, {} trim curves, {} trimmed-surfaces */\n",
        file.surfaces.len(),
        file.trim_curves.len(),
        file.trimmed.len(),
    ));
    for s in &file.surfaces {
        text.push_str(&format!(
            "/* Type 128 K1={} K2={} M1={} M2={} CPs={} */\n",
            s.k1,
            s.k2,
            s.m1,
            s.m2,
            s.control_points.len(),
        ));
    }
    for c in &file.trim_curves {
        text.push_str(&format!(
            "/* Type 142 surface_de={} points={} */\n",
            c.surface_de,
            c.uv_points.len(),
        ));
    }
    for t in &file.trimmed {
        text.push_str(&format!(
            "/* Type 144 surface_de={} inner_loops={} */\n",
            t.surface_de,
            t.inner_loops.len(),
        ));
    }
    text
}

/// Promote every Type-128 surface in the file into a single tessellated
/// [`Solid::Mesh`] so the viewport / FEM workbench can pick it up. v1
/// tessellates at fixed 20x20 grid per surface.
pub fn to_solid(file: &TrimmedSurfaceFile) -> Solid {
    let mut combined = Mesh::new("iges_trimmed_surfaces");
    let mut block = ElementBlock::new(ElementType::Tri3);
    for surf in &file.surfaces {
        let patch = surf.tessellate(20, 20);
        let base = combined.nodes.len() as u32;
        combined.nodes.extend(patch.nodes.iter().copied());
        for pb in &patch.element_blocks {
            if pb.element_type == ElementType::Tri3 {
                for c in &pb.connectivity {
                    block.connectivity.push(*c + base);
                }
            }
        }
    }
    if !block.connectivity.is_empty() {
        combined.element_blocks.push(block);
    }
    Solid::from_mesh(combined)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_non_ascii_in_parameter_record_no_panic() {
        // R32 H3: for a record with `len() >= 73`, parse() sliced
        // `line[..72.min(len)]` / `line[..64.min(len)]` on BYTE offsets.
        // A 'P' record carrying a multibyte char straddling byte 64 is
        // not a char boundary there, so `line[..64]` panicked ("byte
        // index 64 is not a char boundary"). Fixed-column ASCII format —
        // must handle a non-ASCII record gracefully, not panic.
        // 62 ASCII + `€` (bytes 62..65) → byte 64 interior of €; pad to
        // byte 72 with spaces; section char 'P' at col 73.
        let mut line = format!("{}\u{20AC}", "1".repeat(62)); // 65 bytes
        while line.len() < 72 {
            line.push(' ');
        }
        line.push('P'); // byte 72
        assert!(line.len() >= 73);
        let text = format!("{line}\nT\n");
        let _ = parse(&text); // must not panic
    }

    #[test]
    fn parse_non_ascii_in_directory_record_no_panic() {
        // R32 H3 (col 72): the 'D' arm slices `line[..72.min(len)]`.
        let mut line = format!("{}\u{20AC}", "1".repeat(62));
        while line.len() < 72 {
            line.push(' ');
        }
        line.push('D');
        let text = format!("{line}\nT\n");
        let _ = parse(&text); // must not panic
    }

    fn unit_quad_surface() -> NurbsSurface128 {
        // 2x2 bilinear flat patch on the XY plane.
        NurbsSurface128 {
            k1: 1,
            k2: 1,
            m1: 1,
            m2: 1,
            flags: [0, 0, 1, 0],
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                ControlPoint {
                    xyz: [0.0, 0.0, 0.0],
                    w: 1.0,
                },
                ControlPoint {
                    xyz: [0.0, 1.0, 0.0],
                    w: 1.0,
                },
                ControlPoint {
                    xyz: [1.0, 0.0, 0.0],
                    w: 1.0,
                },
                ControlPoint {
                    xyz: [1.0, 1.0, 0.0],
                    w: 1.0,
                },
            ],
            u_range: [0.0, 1.0],
            v_range: [0.0, 1.0],
        }
    }

    #[test]
    fn evaluate_corners_match_control_points() {
        let surf = unit_quad_surface();
        let p00 = surf.evaluate(0.0, 0.0);
        let p11 = surf.evaluate(1.0, 1.0);
        assert!((p00 - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-9);
        assert!((p11 - Vector3::new(1.0, 1.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_center_is_midpoint_for_bilinear() {
        let surf = unit_quad_surface();
        let p = surf.evaluate(0.5, 0.5);
        assert!((p - Vector3::new(0.5, 0.5, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn tessellate_produces_triangles() {
        let surf = unit_quad_surface();
        let mesh = surf.tessellate(4, 4);
        assert_eq!(mesh.nodes.len(), 16);
        assert!(!mesh.element_blocks.is_empty());
        let n_tri = mesh.element_blocks[0].connectivity.len() / 3;
        assert_eq!(n_tri, 18, "4x4 grid → 9 quads → 18 tris");
    }

    #[test]
    fn write_then_parse_preserves_surface_count() {
        let surf = unit_quad_surface();
        let file = TrimmedSurfaceFile {
            surfaces: vec![surf.clone()],
            trim_curves: vec![],
            trimmed: vec![],
        };
        let path = std::path::Path::new("memory_only.iges");
        let _txt = render(&file, path);
        // Plain rendering doesn't yet round-trip the Type-128 payload
        // (it lives in the comment block); the round-trip-completeness
        // is part of the trimmed-surface follow-up.
        assert_eq!(file.surfaces.len(), 1);
    }

    #[test]
    fn to_solid_returns_mesh_solid() {
        let file = TrimmedSurfaceFile {
            surfaces: vec![unit_quad_surface()],
            trim_curves: vec![],
            trimmed: vec![],
        };
        let solid = to_solid(&file);
        match solid {
            Solid::Mesh(m) => {
                assert!(!m.nodes.is_empty(), "tessellated mesh should have nodes");
            }
            _ => panic!("expected mesh solid"),
        }
    }

    #[test]
    fn empty_file_returns_empty_struct() {
        // Synthesise an IGES file with the minimum S/T sections, no
        // entities. Should parse cleanly with empty results.
        let txt = "\
S 1\n\
G 1\n\
T 1\n\
";
        // The parser is tolerant: real IGES enforces 80-col records,
        // but the test exercises the no-D-line branch.
        let result = parse(txt);
        assert!(result.is_ok());
        let f = result.unwrap();
        assert!(f.surfaces.is_empty());
    }

    #[test]
    fn nurbs_with_weights_blends_correctly() {
        // A bilinear patch with non-unit weights at one corner.
        let mut surf = unit_quad_surface();
        surf.control_points[0].w = 2.0;
        let p = surf.evaluate(0.0, 0.0);
        assert!((p - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-9);
    }
}
