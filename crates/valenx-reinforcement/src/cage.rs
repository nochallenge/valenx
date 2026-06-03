//! Rebar cages (longitudinal bars + transverse hoops) + recipe
//! generators for beams / columns / slabs.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::ReinforcementError;
use crate::rebar::{Rebar, RebarGrade, RebarShape};

/// A complete reinforcement cage — longitudinal bars + transverse
/// hoops (or spiral) + spacings + cover thickness.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RebarCage {
    /// Longitudinal bars (each carries an `origin` in cage-local
    /// coordinates — the bar itself uses [`RebarShape::Straight`] for
    /// the centreline).
    pub longitudinal: Vec<PlacedBar>,
    /// Transverse hoops (typically [`RebarShape::U`] open-top stirrups
    /// or a single [`RebarShape::Spiral`] for round columns).
    pub transverse: Vec<PlacedBar>,
    /// Centre-to-centre spacing of transverse hoops (m).
    pub hoop_spacing: f64,
    /// Concrete cover thickness (m).
    pub cover: f64,
}

/// A rebar instance placed at `origin` within a cage's local frame.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlacedBar {
    /// The rebar.
    pub bar: Rebar,
    /// Where the bar's local-origin starts.
    pub origin: Vector3<f64>,
}

/// Default cover (30 mm = 0.03 m) — common for above-ground
/// columns / beams not exposed to weather.
const DEFAULT_COVER: f64 = 0.03;

/// Generate a beam cage with `n_long_bars` evenly-spaced bars across
/// the width on both top and bottom faces, plus open-top U-stirrups
/// spaced every `hoop_spacing` along the length.
///
/// `width`, `depth`, `length` in m.
pub fn generate_beam(
    width: f64,
    depth: f64,
    length: f64,
    n_long_bars: usize,
    hoop_spacing: f64,
) -> Result<RebarCage, ReinforcementError> {
    if width <= 0.0 || depth <= 0.0 || length <= 0.0 {
        return Err(ReinforcementError::BadParameter {
            name: "width|depth|length",
            reason: format!("must be > 0, got {width}x{depth}x{length}"),
        });
    }
    if n_long_bars < 2 {
        return Err(ReinforcementError::BadParameter {
            name: "n_long_bars",
            reason: "must be >= 2".into(),
        });
    }
    if hoop_spacing <= 0.0 {
        return Err(ReinforcementError::BadParameter {
            name: "hoop_spacing",
            reason: "must be > 0".into(),
        });
    }
    let cover = DEFAULT_COVER;
    let mut cage = RebarCage {
        cover,
        hoop_spacing,
        ..Default::default()
    };
    // Longitudinal bars: top face + bottom face.
    let inner_w = (width - 2.0 * cover).max(0.001);
    let dx = inner_w / (n_long_bars as f64 - 1.0);
    for face_z in [cover, depth - cover] {
        for i in 0..n_long_bars {
            let x = cover + dx * (i as f64);
            cage.longitudinal.push(PlacedBar {
                bar: Rebar {
                    diameter_mm: 16.0,
                    length_m: length,
                    shape: RebarShape::Straight { length },
                    grade: RebarGrade::G60,
                },
                origin: Vector3::new(x, 0.0, face_z),
            });
        }
    }
    // Transverse U-stirrups along the length.
    let mut y = hoop_spacing * 0.5;
    while y <= length - hoop_spacing * 0.5 {
        let inner_h = (depth - 2.0 * cover).max(0.001);
        cage.transverse.push(PlacedBar {
            bar: Rebar {
                diameter_mm: 10.0,
                length_m: 2.0 * inner_w + 2.0 * inner_h,
                shape: RebarShape::U {
                    width: inner_w,
                    height: inner_h,
                },
                grade: RebarGrade::G60,
            },
            origin: Vector3::new(cover, y, depth - cover),
        });
        y += hoop_spacing;
    }
    Ok(cage)
}

/// Generate a column cage — `n_bars_per_face` longitudinal bars on
/// each of the 4 faces plus rectangular (U-stirrup) hoops at
/// `hoop_spacing` along the height. For round columns prefer a
/// [`RebarShape::Spiral`]-based cage in v1.5+.
///
/// `width`, `depth`, `height` in m.
pub fn generate_column(
    width: f64,
    depth: f64,
    height: f64,
    n_bars_per_face: usize,
    hoop_spacing: f64,
) -> Result<RebarCage, ReinforcementError> {
    if width <= 0.0 || depth <= 0.0 || height <= 0.0 {
        return Err(ReinforcementError::BadParameter {
            name: "width|depth|height",
            reason: format!("must be > 0, got {width}x{depth}x{height}"),
        });
    }
    if n_bars_per_face < 2 {
        return Err(ReinforcementError::BadParameter {
            name: "n_bars_per_face",
            reason: "must be >= 2".into(),
        });
    }
    if hoop_spacing <= 0.0 {
        return Err(ReinforcementError::BadParameter {
            name: "hoop_spacing",
            reason: "must be > 0".into(),
        });
    }
    let cover = DEFAULT_COVER;
    let mut cage = RebarCage {
        cover,
        hoop_spacing,
        ..Default::default()
    };
    // Helper closure to create one vertical bar.
    let push_long = |cage: &mut RebarCage, x: f64, y: f64| {
        cage.longitudinal.push(PlacedBar {
            bar: Rebar {
                diameter_mm: 20.0,
                length_m: height,
                shape: RebarShape::Straight { length: height },
                grade: RebarGrade::G60,
            },
            origin: Vector3::new(x, y, 0.0),
        });
    };
    // 4 faces — corners are shared but for v1 we just place them per
    // face (BOM-honest count = 4 corners + edge bars).
    let xs: Vec<f64> = (0..n_bars_per_face)
        .map(|i| cover + (width - 2.0 * cover) * (i as f64) / (n_bars_per_face as f64 - 1.0))
        .collect();
    let ys: Vec<f64> = (0..n_bars_per_face)
        .map(|i| cover + (depth - 2.0 * cover) * (i as f64) / (n_bars_per_face as f64 - 1.0))
        .collect();
    for &x in &xs {
        push_long(&mut cage, x, cover);
        push_long(&mut cage, x, depth - cover);
    }
    for &y in &ys {
        if y > cover + 1e-9 && y < depth - cover - 1e-9 {
            push_long(&mut cage, cover, y);
            push_long(&mut cage, width - cover, y);
        }
    }
    // Hoops along the height.
    let inner_w = (width - 2.0 * cover).max(0.001);
    let inner_d = (depth - 2.0 * cover).max(0.001);
    let mut z = hoop_spacing * 0.5;
    while z <= height - hoop_spacing * 0.5 {
        cage.transverse.push(PlacedBar {
            bar: Rebar {
                diameter_mm: 10.0,
                length_m: 2.0 * inner_w + 2.0 * inner_d,
                shape: RebarShape::U {
                    width: inner_w,
                    height: inner_d,
                },
                grade: RebarGrade::G60,
            },
            origin: Vector3::new(cover, cover, z),
        });
        z += hoop_spacing;
    }
    Ok(cage)
}

/// Generate a slab cage — straight bars in two orthogonal directions
/// at `bar_spacing_x` along Y and `bar_spacing_y` along X. Returns
/// the flat list of [`Rebar`] instances (no transverse element).
///
/// `width` along X, `length` along Y, `thickness` along Z.
pub fn generate_slab(
    width: f64,
    length: f64,
    thickness: f64,
    bar_spacing_x: f64,
    bar_spacing_y: f64,
) -> Result<Vec<PlacedBar>, ReinforcementError> {
    if width <= 0.0 || length <= 0.0 || thickness <= 0.0 {
        return Err(ReinforcementError::BadParameter {
            name: "width|length|thickness",
            reason: format!("must be > 0, got {width}x{length}x{thickness}"),
        });
    }
    if bar_spacing_x <= 0.0 || bar_spacing_y <= 0.0 {
        return Err(ReinforcementError::BadParameter {
            name: "bar_spacing_x|bar_spacing_y",
            reason: "must be > 0".into(),
        });
    }
    let cover = DEFAULT_COVER;
    let mut out = Vec::new();
    // Bars along X spaced by `bar_spacing_x` along Y.
    let mut y = cover;
    while y <= length - cover {
        out.push(PlacedBar {
            bar: Rebar {
                diameter_mm: 12.0,
                length_m: width - 2.0 * cover,
                shape: RebarShape::Straight {
                    length: width - 2.0 * cover,
                },
                grade: RebarGrade::G60,
            },
            origin: Vector3::new(cover, y, cover),
        });
        y += bar_spacing_x;
    }
    // Bars along Y spaced by `bar_spacing_y` along X.
    let mut x = cover;
    while x <= width - cover {
        out.push(PlacedBar {
            bar: Rebar {
                diameter_mm: 12.0,
                length_m: length - 2.0 * cover,
                shape: RebarShape::Straight {
                    length: length - 2.0 * cover,
                },
                grade: RebarGrade::G60,
            },
            origin: Vector3::new(x, cover, cover + 0.012), // slight z offset for 2nd layer
        });
        x += bar_spacing_y;
    }
    Ok(out)
}

/// Number of facets in the circular cross-section ring. 12 reads as a
/// smooth round bar at viewport scale without ballooning the mesh.
const SECTION_FACETS: usize = 12;

/// Tessellate a cage into one Tri3 [`Mesh`] for the viewport.
///
/// Phase 33.5: each bar is swept with a **circular cross-section**
/// (a `SECTION_FACETS`-sided ring at the bar's true nominal radius),
/// replacing the earlier 4-sided square approximation. A ring is
/// placed at every centreline vertex; the ring frame is
/// parallel-transported segment to segment so the round section does
/// not twist around bends, and consecutive rings are connected with
/// quad walls. This is a genuine circular sweep, so a Ø16 mm bar
/// renders round instead of as a square prism.
pub fn to_mesh(cage: &RebarCage) -> Mesh {
    let mut out = Mesh::new("rebar_cage");
    let mut block = ElementBlock::new(ElementType::Tri3);
    for b in &cage.longitudinal {
        push_bar_circular(b, &mut out, &mut block);
    }
    for b in &cage.transverse {
        push_bar_circular(b, &mut out, &mut block);
    }
    out.element_blocks.push(block);
    out.recompute_stats();
    out
}

/// Sweep one bar's centreline polyline as a circular-section tube.
fn push_bar_circular(bar: &PlacedBar, out: &mut Mesh, block: &mut ElementBlock) {
    // True nominal radius: diameter mm → m, halved.
    let r = bar.bar.diameter_mm * 0.5e-3;
    let pts = bar.bar.shape.to_polyline();
    if pts.len() < 2 || r <= 0.0 {
        return;
    }

    // Per-vertex tangents (averaged at interior vertices).
    let n_pts = pts.len();
    let mut tangents = Vec::with_capacity(n_pts);
    for i in 0..n_pts {
        let incoming = if i > 0 {
            (pts[i] - pts[i - 1]).try_normalize(1e-12)
        } else {
            None
        };
        let outgoing = if i + 1 < n_pts {
            (pts[i + 1] - pts[i]).try_normalize(1e-12)
        } else {
            None
        };
        let t = match (incoming, outgoing) {
            (Some(a), Some(b)) => (a + b).try_normalize(1e-12).unwrap_or(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => Vector3::x(),
        };
        tangents.push(t);
    }

    // Parallel-transport an orthonormal (u, v) frame along the curve.
    let mut frames: Vec<(Vector3<f64>, Vector3<f64>)> = Vec::with_capacity(n_pts);
    let (mut fu, mut fv) = perp_frame(tangents[0]);
    frames.push((fu, fv));
    for i in 1..n_pts {
        let rot = min_rotation(tangents[i - 1], tangents[i]);
        fu = rot * fu;
        fv = rot * fv;
        frames.push((fu, fv));
    }

    // Lay out one circular ring of `SECTION_FACETS` nodes per vertex.
    let ring_base = out.nodes.len() as u32;
    for (i, p) in pts.iter().enumerate() {
        let centre = bar.origin + p;
        let (u, v) = frames[i];
        for k in 0..SECTION_FACETS {
            let ang = k as f64 / SECTION_FACETS as f64 * std::f64::consts::TAU;
            out.nodes
                .push(centre + u * (r * ang.cos()) + v * (r * ang.sin()));
        }
    }

    // Connect consecutive rings with quad walls (2 tris each).
    let facets = SECTION_FACETS as u32;
    for seg in 0..(n_pts as u32 - 1) {
        let a = ring_base + seg * facets;
        let b = ring_base + (seg + 1) * facets;
        for k in 0..facets {
            let k1 = (k + 1) % facets;
            // quad (a+k, a+k1, b+k1, b+k) → two triangles.
            block.connectivity.push(a + k);
            block.connectivity.push(a + k1);
            block.connectivity.push(b + k1);
            block.connectivity.push(a + k);
            block.connectivity.push(b + k1);
            block.connectivity.push(b + k);
        }
    }
}

/// Two orthonormal vectors spanning the plane ⟂ unit vector `t`.
fn perp_frame(t: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let seed = if t.dot(&Vector3::z()).abs() < 0.9 {
        Vector3::z()
    } else {
        Vector3::x()
    };
    let u = t.cross(&seed).normalize();
    let v = t.cross(&u);
    (u, v)
}

/// Minimal rotation carrying unit vector `from` onto unit vector `to`.
fn min_rotation(from: Vector3<f64>, to: Vector3<f64>) -> nalgebra::UnitQuaternion<f64> {
    nalgebra::UnitQuaternion::rotation_between(&from, &to)
        .unwrap_or_else(nalgebra::UnitQuaternion::identity)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn beam_cage_has_longitudinal_and_hoops() {
        let cage = generate_beam(0.3, 0.5, 4.0, 4, 0.2).unwrap();
        // Top + bottom faces → 2 × 4 = 8 long bars.
        assert_eq!(cage.longitudinal.len(), 8);
        assert!(!cage.transverse.is_empty());
    }

    #[test]
    fn beam_bad_param_errors() {
        assert!(matches!(
            generate_beam(-1.0, 0.5, 4.0, 4, 0.2),
            Err(ReinforcementError::BadParameter { .. })
        ));
        assert!(matches!(
            generate_beam(0.3, 0.5, 4.0, 1, 0.2),
            Err(ReinforcementError::BadParameter { .. })
        ));
    }

    #[test]
    fn column_cage_has_all_4_faces() {
        let cage = generate_column(0.4, 0.4, 3.0, 3, 0.2).unwrap();
        // 4 faces with corner sharing: 4 corners + 2 mid-edge bars * 4 faces.
        // Actual count for n=3: 6 (top/bottom X faces × 3) + 2 mid-edge Y bars × 2 = 12.
        assert!(cage.longitudinal.len() >= 8);
        assert!(!cage.transverse.is_empty());
    }

    #[test]
    fn slab_two_directions() {
        let bars = generate_slab(3.0, 3.0, 0.2, 0.2, 0.2).unwrap();
        assert!(bars.len() >= 4);
    }

    #[test]
    fn to_mesh_produces_triangles() {
        let cage = generate_beam(0.3, 0.5, 4.0, 4, 0.5).unwrap();
        let m = to_mesh(&cage);
        assert!(!m.nodes.is_empty());
        assert!(m.total_elements() > 0);
    }

    // --- Phase 33.5 circular-section tests ---

    #[test]
    fn straight_bar_section_is_circular() {
        // A single straight bar: every ring vertex must sit at the
        // bar's nominal radius from the centreline — i.e. the cross
        // section is a real circle, not a square.
        let mut cage = RebarCage::default();
        let bar = Rebar {
            diameter_mm: 20.0,
            length_m: 1.0,
            shape: RebarShape::Straight { length: 1.0 },
            grade: RebarGrade::G60,
        };
        cage.longitudinal.push(PlacedBar {
            bar,
            origin: Vector3::zeros(),
        });
        let m = to_mesh(&cage);
        // 2 centreline vertices × SECTION_FACETS ring nodes.
        assert_eq!(m.nodes.len(), 2 * SECTION_FACETS);
        let r = 20.0 * 0.5e-3; // 10 mm in metres
        // The straight bar runs along +X; ring nodes lie in the YZ
        // plane offset, so their (y,z) radius equals r.
        for node in &m.nodes {
            let radial = (node.y * node.y + node.z * node.z).sqrt();
            assert!(
                (radial - r).abs() < 1e-9,
                "ring node radius {radial} should equal nominal {r}"
            );
        }
    }

    #[test]
    fn circular_section_has_more_facets_than_a_square() {
        // SECTION_FACETS sides → SECTION_FACETS quads → 2× triangles
        // per segment. A square section would give only 4 quads.
        let mut cage = RebarCage::default();
        cage.longitudinal.push(PlacedBar {
            bar: Rebar {
                diameter_mm: 16.0,
                length_m: 1.0,
                shape: RebarShape::Straight { length: 2.0 },
                grade: RebarGrade::G60,
            },
            origin: Vector3::zeros(),
        });
        let m = to_mesh(&cage);
        // One segment → SECTION_FACETS quads → 2*SECTION_FACETS tris.
        // SECTION_FACETS (12) is far more than a 4-sided square's 8
        // triangles per segment — proves the section is round.
        assert_eq!(m.total_elements(), 2 * SECTION_FACETS);
    }
}
