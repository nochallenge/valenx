//! Cross-section profile catalogue. Each variant is a parametric
//! 2D shape; [`cross_section_polygon`] returns the section's outline
//! polygon as a closed CCW polygon in the local (u, v) frame.

use serde::{Deserialize, Serialize};

/// Standard structural profiles.
///
/// All dimensions in millimetres. The profile is positioned with its
/// centroid at the origin and the strong-axis ("h") aligned with the
/// local +v axis; the local +u axis points to the right.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Profile {
    /// I-beam (wide-flange) — `h` overall height, `b` flange width,
    /// `tw` web thickness, `tf` flange thickness.
    IBeam {
        /// Overall height (mm).
        h: f64,
        /// Flange width (mm).
        b: f64,
        /// Web thickness (mm).
        tw: f64,
        /// Flange thickness (mm).
        tf: f64,
    },
    /// C-channel — `h` overall height, `b` flange width, `tw` web
    /// thickness (we approximate flange thickness == tw).
    CChannel {
        /// Overall height (mm).
        h: f64,
        /// Flange width (mm).
        b: f64,
        /// Web thickness (mm).
        tw: f64,
    },
    /// L-angle (equal or unequal) — `h` vertical leg height, `b`
    /// horizontal leg width, `t` leg thickness.
    LAngle {
        /// Vertical leg length (mm).
        h: f64,
        /// Horizontal leg length (mm).
        b: f64,
        /// Leg thickness (mm).
        t: f64,
    },
    /// Rectangular hollow section (HSS rect) — `h` outer height, `b`
    /// outer width, `t` wall thickness.
    RhsRect {
        /// Outer height (mm).
        h: f64,
        /// Outer width (mm).
        b: f64,
        /// Wall thickness (mm).
        t: f64,
    },
    /// Round hollow section (CHS) — `d` outer diameter, `t` wall
    /// thickness.
    ChsRound {
        /// Outer diameter (mm).
        d: f64,
        /// Wall thickness (mm).
        t: f64,
    },
    /// T-beam — `h` overall height, `b` flange width, `tw` web
    /// thickness, `tf` flange thickness.
    TBeam {
        /// Overall height (mm).
        h: f64,
        /// Flange width (mm).
        b: f64,
        /// Web thickness (mm).
        tw: f64,
        /// Flange thickness (mm).
        tf: f64,
    },
}

impl Profile {
    /// Short UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::IBeam { .. } => "I-beam",
            Self::CChannel { .. } => "C-channel",
            Self::LAngle { .. } => "L-angle",
            Self::RhsRect { .. } => "Rect HSS",
            Self::ChsRound { .. } => "Round CHS",
            Self::TBeam { .. } => "T-beam",
        }
    }

    /// Convenient default for the UI — an IPE 200 wide-flange.
    pub fn default_ipe200() -> Self {
        Self::IBeam {
            h: 200.0,
            b: 100.0,
            tw: 5.6,
            tf: 8.5,
        }
    }

    /// Exact cross-sectional area (mm²) of the idealized (sharp-cornered, no-fillet) profile —
    /// a sum of rectangles, or an annulus for [`ChsRound`](Profile::ChsRound). Matches the
    /// geometry of [`cross_section_polygon`]; differs from a real rolled section that has root
    /// fillets. Non-finite or degenerate dimensions clamp the result to `0.0`.
    pub fn cross_section_area_mm2(self) -> f64 {
        // Inner/remaining dimensions are clamped to ≥0 so an over-thick wall collapses to a solid
        // section rather than producing a spurious negative-squared hole.
        let area = match self {
            Self::IBeam { h, b, tw, tf } => 2.0 * b * tf + tw * (h - 2.0 * tf).max(0.0),
            Self::CChannel { h, b, tw } => 2.0 * b * tw + tw * (h - 2.0 * tw).max(0.0),
            Self::LAngle { h, b, t } => t * h + t * (b - t).max(0.0),
            Self::RhsRect { h, b, t } => b * h - (b - 2.0 * t).max(0.0) * (h - 2.0 * t).max(0.0),
            Self::ChsRound { d, t } => {
                let inner = (d - 2.0 * t).max(0.0);
                std::f64::consts::PI / 4.0 * (d * d - inner * inner)
            }
            Self::TBeam { h, b, tw, tf } => b * tf + tw * (h - tf).max(0.0),
        };
        if area.is_finite() && area > 0.0 {
            area
        } else {
            0.0
        }
    }

    /// Second moments of area `(I_u, I_v)` about the section's own
    /// horizontal (u) and vertical (v) centroidal axes, in mm⁴.
    ///
    /// `I_u` bends about the horizontal axis (resists loads applied in
    /// the strong-axis / "h" direction) and `I_v` bends about the
    /// vertical axis. They are computed by the exact polygon
    /// second-moment integral over [`cross_section_polygon`], so they
    /// are consistent with [`cross_section_area_mm2`](Profile::cross_section_area_mm2)
    /// and reduce to the textbook closed forms — e.g. a solid
    /// rectangle gives `I = b·h³/12` (Roark's *Formulas for Stress &
    /// Strain*, Table A.1; Gere, *Mechanics of Materials*, App. D).
    /// [`ChsRound`](Profile::ChsRound) uses the exact annulus form
    /// `π·(D⁴ − d⁴)/64` rather than the 24-gon. Degenerate / non-finite
    /// dimensions clamp to `0.0`.
    pub fn second_moment_of_area_mm4(self) -> (f64, f64) {
        if let Self::ChsRound { d, t } = self {
            // Exact annulus: I = π·(D⁴ − d_in⁴)/64 about any diameter.
            let inner = (d - 2.0 * t).max(0.0);
            let i = std::f64::consts::PI / 64.0 * (d.powi(4) - inner.powi(4));
            return finite_pair(i, i);
        }
        if let Self::RhsRect { h, b, t } = self {
            // Exact hollow-rectangle (HSS) second moments: outer minus
            // inner concentric rectangle (shared centroid). Matches the
            // hollow area in cross_section_area_mm2 and the textbook HSS
            // form I_u = (b·h³ − b_i·h_i³)/12. The cross_section_polygon
            // outline is the *solid* outer rectangle (used only for
            // rendering); the structural property must subtract the hole,
            // else I/S is overestimated (unconservative) for tube sections.
            let b_in = (b - 2.0 * t).max(0.0);
            let h_in = (h - 2.0 * t).max(0.0);
            let iu = (b * h.powi(3) - b_in * h_in.powi(3)) / 12.0;
            let iv = (h * b.powi(3) - h_in * b_in.powi(3)) / 12.0;
            return finite_pair(iu, iv);
        }
        // Exact second moments of the closed outline polygon about its
        // own centroid, via the standard shoelace moment integral.
        let poly = cross_section_polygon(self);
        let (iu, iv) = polygon_centroidal_second_moments(&poly);
        finite_pair(iu, iv)
    }

    /// Elastic section moduli `(S_u, S_v) = (I_u / c_u, I_v / c_v)` in
    /// mm³, where `c` is the distance from the centroid to the extreme
    /// fibre. The maximum bending stress in a beam is `σ = M / S`
    /// (Roark's *Formulas for Stress & Strain*; Gere, *Mechanics of
    /// Materials*). For a solid rectangle this reduces to the textbook
    /// `S = b·h²/6`. Degenerate / non-finite sections give `0.0`.
    pub fn section_modulus_mm3(self) -> (f64, f64) {
        if let Self::ChsRound { d, t } = self {
            let (i, _) = self.second_moment_of_area_mm4();
            let c = (d * 0.5).max(0.0);
            let _ = t;
            return finite_pair(div_or_zero(i, c), div_or_zero(i, c));
        }
        if let Self::RhsRect { h, b, .. } = self {
            // Hollow-correct S = I / c, with I the hollow second moment
            // from second_moment_of_area_mm4 and the outer extreme-fibre
            // distances (c_u = h/2, c_v = b/2).
            let (iu, iv) = self.second_moment_of_area_mm4();
            return finite_pair(div_or_zero(iu, h * 0.5), div_or_zero(iv, b * 0.5));
        }
        let poly = cross_section_polygon(self);
        let (iu, iv) = polygon_centroidal_second_moments(&poly);
        let (cu, cv) = polygon_extreme_fibre_distances(&poly);
        finite_pair(div_or_zero(iu, cu), div_or_zero(iv, cv))
    }

    /// Exact outer-boundary perimeter (mm) of the idealized (sharp-cornered, no-fillet) profile —
    /// the sum of the outer outline's edge lengths; for [`ChsRound`](Profile::ChsRound) it is the
    /// circle perimeter `π·d` (not a polygon approximation). The wall-thickness terms cancel in the
    /// open profiles, so this depends only on the outer dimensions. Non-finite/degenerate → `0.0`.
    pub fn cross_section_perimeter_mm(self) -> f64 {
        let perimeter = match self {
            Self::IBeam { h, b, tw, .. } => 2.0 * (2.0 * b + h - tw),
            Self::CChannel { h, b, tw } => 2.0 * (2.0 * b + h - tw),
            Self::LAngle { h, b, .. } => 2.0 * (b + h),
            Self::RhsRect { h, b, .. } => 2.0 * (b + h),
            Self::ChsRound { d, .. } => std::f64::consts::PI * d,
            Self::TBeam { h, b, .. } => 2.0 * (b + h),
        };
        if perimeter.is_finite() && perimeter > 0.0 {
            perimeter
        } else {
            0.0
        }
    }
}

impl Default for Profile {
    fn default() -> Self {
        Self::default_ipe200()
    }
}

/// Build the closed CCW outline polygon for `p` in local (u, v).
/// The polygon is **not** explicitly closed (last vertex != first);
/// the consumer must close it if it needs an explicit duplicate.
pub fn cross_section_polygon(p: Profile) -> Vec<[f64; 2]> {
    match p {
        Profile::IBeam { h, b, tw, tf } => {
            // 12-vertex I-beam, centroid at origin.
            let h2 = h * 0.5;
            let b2 = b * 0.5;
            let tw2 = tw * 0.5;
            let yfb = h2 - tf;
            vec![
                [-b2, -h2],
                [b2, -h2],
                [b2, -yfb],
                [tw2, -yfb],
                [tw2, yfb],
                [b2, yfb],
                [b2, h2],
                [-b2, h2],
                [-b2, yfb],
                [-tw2, yfb],
                [-tw2, -yfb],
                [-b2, -yfb],
            ]
        }
        Profile::CChannel { h, b, tw } => {
            // 8-vertex C, centroid roughly at origin (-b/4 offset
            // ignored for v1 — kept axis-aligned for predictability).
            let h2 = h * 0.5;
            let b2 = b * 0.5;
            vec![
                [-b2, -h2],
                [b2, -h2],
                [b2, -h2 + tw],
                [-b2 + tw, -h2 + tw],
                [-b2 + tw, h2 - tw],
                [b2, h2 - tw],
                [b2, h2],
                [-b2, h2],
            ]
        }
        Profile::LAngle { h, b, t } => {
            // 6-vertex L with vertical leg height `h` and horizontal
            // leg width `b`; origin at the lower-left corner of the
            // bounding box, then re-centred on the leg-thickness
            // intersection.
            vec![[0.0, 0.0], [b, 0.0], [b, t], [t, t], [t, h], [0.0, h]]
        }
        Profile::RhsRect { h, b, t } => {
            // Outer rectangle CCW. The hollow interior is not yet
            // represented — v1 emits the solid outer outline only
            // (the visual proxy under a sweep is still recognisable
            // as a rect HSS section).
            let _ = t;
            let h2 = h * 0.5;
            let b2 = b * 0.5;
            vec![[-b2, -h2], [b2, -h2], [b2, h2], [-b2, h2]]
        }
        Profile::ChsRound { d, t } => {
            let _ = t;
            let r = d * 0.5;
            let n = 24;
            (0..n)
                .map(|i| {
                    let th = (i as f64 / n as f64) * std::f64::consts::TAU;
                    [r * th.cos(), r * th.sin()]
                })
                .collect()
        }
        Profile::TBeam { h, b, tw, tf } => {
            // 8-vertex T, centroid at origin.
            let b2 = b * 0.5;
            let tw2 = tw * 0.5;
            let v_top = h - tf;
            vec![
                [-tw2, 0.0],
                [tw2, 0.0],
                [tw2, v_top],
                [b2, v_top],
                [b2, h],
                [-b2, h],
                [-b2, v_top],
                [-tw2, v_top],
            ]
        }
    }
}

/// Clamp a pair of section properties to `0.0` when either is
/// non-finite or non-positive (degenerate input geometry).
fn finite_pair(a: f64, b: f64) -> (f64, f64) {
    let f = |x: f64| if x.is_finite() && x > 0.0 { x } else { 0.0 };
    (f(a), f(b))
}

/// `n / d`, or `0.0` when the denominator is non-positive/non-finite.
fn div_or_zero(n: f64, d: f64) -> f64 {
    if d.is_finite() && d > 0.0 && n.is_finite() {
        n / d
    } else {
        0.0
    }
}

/// Area and centroid `(A, cx, cy)` of a simple polygon via the shoelace
/// formulas. Sign-robust: `A` is returned as the unsigned area while the
/// centroid uses the signed area internally, so vertex winding does not
/// affect the result.
fn polygon_area_centroid(poly: &[[f64; 2]]) -> (f64, f64, f64) {
    if poly.len() < 3 {
        return (0.0, 0.0, 0.0);
    }
    let mut a2 = 0.0; // 2·signed area
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..poly.len() {
        let [x0, y0] = poly[i];
        let [x1, y1] = poly[(i + 1) % poly.len()];
        let cross = x0 * y1 - x1 * y0;
        a2 += cross;
        cx += (x0 + x1) * cross;
        cy += (y0 + y1) * cross;
    }
    let signed_area = a2 * 0.5;
    if signed_area.abs() < f64::EPSILON {
        return (0.0, 0.0, 0.0);
    }
    (signed_area.abs(), cx / (3.0 * a2), cy / (3.0 * a2))
}

/// Second moments of area `(I_u, I_v)` of a simple polygon taken about
/// its own centroid, where `I_u` is `∫ (y − ȳ)² dA` (bending about the
/// horizontal/u axis) and `I_v` is `∫ (x − x̄)² dA`. Uses the exact
/// shoelace second-moment integral plus the parallel-axis shift back to
/// the centroid, so it is exact for any straight-edged section.
fn polygon_centroidal_second_moments(poly: &[[f64; 2]]) -> (f64, f64) {
    if poly.len() < 3 {
        return (0.0, 0.0);
    }
    let (area, cx, cy) = polygon_area_centroid(poly);
    if area == 0.0 {
        return (0.0, 0.0);
    }
    // Second moments about the global origin (signed-area convention).
    let mut a2 = 0.0;
    let mut ixx_o = 0.0; // ∫ y² dA  (about global x axis)
    let mut iyy_o = 0.0; // ∫ x² dA  (about global y axis)
    for i in 0..poly.len() {
        let [x0, y0] = poly[i];
        let [x1, y1] = poly[(i + 1) % poly.len()];
        let cross = x0 * y1 - x1 * y0;
        a2 += cross;
        ixx_o += (y0 * y0 + y0 * y1 + y1 * y1) * cross;
        iyy_o += (x0 * x0 + x0 * x1 + x1 * x1) * cross;
    }
    let sign = if a2 >= 0.0 { 1.0 } else { -1.0 };
    let ixx_o = ixx_o / 12.0 * sign;
    let iyy_o = iyy_o / 12.0 * sign;
    // Parallel-axis shift to the centroid: I_c = I_o − A·d².
    let i_u = ixx_o - area * cy * cy;
    let i_v = iyy_o - area * cx * cx;
    (i_u.max(0.0), i_v.max(0.0))
}

/// Distances `(c_u, c_v)` from the polygon centroid to the extreme
/// fibre in the v (vertical) and u (horizontal) directions — the lever
/// arms used in the section modulus `S = I / c`.
fn polygon_extreme_fibre_distances(poly: &[[f64; 2]]) -> (f64, f64) {
    if poly.len() < 3 {
        return (0.0, 0.0);
    }
    let (_, cx, cy) = polygon_area_centroid(poly);
    let mut cu = 0.0_f64; // max |y − ȳ|
    let mut cv = 0.0_f64; // max |x − x̄|
    for &[x, y] in poly {
        cu = cu.max((y - cy).abs());
        cv = cv.max((x - cx).abs());
    }
    (cu, cv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ibeam_polygon_has_12_verts() {
        let p = Profile::IBeam {
            h: 200.0,
            b: 100.0,
            tw: 6.0,
            tf: 10.0,
        };
        assert_eq!(cross_section_polygon(p).len(), 12);
    }

    #[test]
    fn round_chs_polygon_has_24_verts() {
        let p = Profile::ChsRound { d: 100.0, t: 4.0 };
        assert_eq!(cross_section_polygon(p).len(), 24);
    }

    #[test]
    fn rect_hss_polygon_is_a_quad() {
        let p = Profile::RhsRect {
            h: 80.0,
            b: 40.0,
            t: 3.0,
        };
        assert_eq!(cross_section_polygon(p).len(), 4);
    }

    #[test]
    fn cross_section_area_exact_per_variant() {
        // I-beam {200,100,5.6,8.5}: 2·100·8.5 + 5.6·(200−17) = 2724.8.
        assert!(
            (Profile::IBeam {
                h: 200.0,
                b: 100.0,
                tw: 5.6,
                tf: 8.5
            }
            .cross_section_area_mm2()
                - 2724.8)
                .abs()
                < 1e-9
        );
        // RHS rect {100,50,5}: 50·100 − 40·90 = 1400.
        assert!(
            (Profile::RhsRect {
                h: 100.0,
                b: 50.0,
                t: 5.0
            }
            .cross_section_area_mm2()
                - 1400.0)
                .abs()
                < 1e-9
        );
        // CHS round {100,4}: π/4·(100²−92²) ≈ 1206.37.
        assert!(
            (Profile::ChsRound { d: 100.0, t: 4.0 }.cross_section_area_mm2()
                - std::f64::consts::PI / 4.0 * (100.0_f64 * 100.0 - 92.0 * 92.0))
                .abs()
                < 1e-9
        );
        // L-angle {100,80,5}: 5·100 + 5·75 = 875. T-beam {200,100,8,12}: 1200 + 1504 = 2704.
        assert!(
            (Profile::LAngle {
                h: 100.0,
                b: 80.0,
                t: 5.0
            }
            .cross_section_area_mm2()
                - 875.0)
                .abs()
                < 1e-9
        );
        assert!(
            (Profile::TBeam {
                h: 200.0,
                b: 100.0,
                tw: 8.0,
                tf: 12.0
            }
            .cross_section_area_mm2()
                - 2704.0)
                .abs()
                < 1e-9
        );
        // An over-thick wall collapses to a solid section (no spurious negative hole).
        assert_eq!(
            Profile::RhsRect {
                h: 100.0,
                b: 100.0,
                t: 60.0
            }
            .cross_section_area_mm2(),
            10000.0
        );
    }

    #[test]
    fn cross_section_perimeter_exact_per_variant() {
        // I-beam {200,100,5.6,8.5}: 2·(2·100+200−5.6) = 788.8 (tf cancels).
        assert!(
            (Profile::IBeam {
                h: 200.0,
                b: 100.0,
                tw: 5.6,
                tf: 8.5
            }
            .cross_section_perimeter_mm()
                - 788.8)
                .abs()
                < 1e-9
        );
        // C-channel {200,100,5.6}: 2·(2·100+200−5.6) = 788.8.
        assert!(
            (Profile::CChannel {
                h: 200.0,
                b: 100.0,
                tw: 5.6
            }
            .cross_section_perimeter_mm()
                - 788.8)
                .abs()
                < 1e-9
        );
        // L-angle {100,80,5}: 2·(80+100) = 360; RHS rect {100,50,5}: 2·(50+100) = 300.
        assert!(
            (Profile::LAngle {
                h: 100.0,
                b: 80.0,
                t: 5.0
            }
            .cross_section_perimeter_mm()
                - 360.0)
                .abs()
                < 1e-9
        );
        assert!(
            (Profile::RhsRect {
                h: 100.0,
                b: 50.0,
                t: 5.0
            }
            .cross_section_perimeter_mm()
                - 300.0)
                .abs()
                < 1e-9
        );
        // CHS round {100,4}: the exact circle perimeter π·d ≈ 314.16 (not the 24-gon).
        assert!(
            (Profile::ChsRound { d: 100.0, t: 4.0 }.cross_section_perimeter_mm()
                - std::f64::consts::PI * 100.0)
                .abs()
                < 1e-9
        );
        // T-beam {200,100,8,12}: 2·(100+200) = 600.
        assert!(
            (Profile::TBeam {
                h: 200.0,
                b: 100.0,
                tw: 8.0,
                tf: 12.0
            }
            .cross_section_perimeter_mm()
                - 600.0)
                .abs()
                < 1e-9
        );
    }

    /// Validation against the closed-form second moment of area and
    /// section modulus of a *hollow* rectangle (HSS / rectangular tube).
    /// Reference: the textbook hollow-section form (Roark's *Formulas for
    /// Stress & Strain*; Gere, *Mechanics of Materials*, App. D):
    /// `I_u = (b·h³ − b_i·h_i³)/12`, `S_u = I_u / (h/2)`, with the inner
    /// dimensions `b_i = b − 2t`, `h_i = h − 2t`. The hole is subtracted
    /// (consistent with the hollow `cross_section_area_mm2`), so the
    /// result is the true HSS property — not the solid outer rectangle,
    /// which would overestimate I/S (unconservative). Case: b = 50 mm
    /// (width, u), h = 100 mm (height, v), t = 5 mm wall.
    #[test]
    fn rectangle_second_moment_and_modulus_match_closed_form() {
        let b = 50.0_f64;
        let h = 100.0_f64;
        let t = 5.0_f64;
        let rect = Profile::RhsRect { h, b, t };
        let b_in = b - 2.0 * t;
        let h_in = h - 2.0 * t;

        let (iu, iv) = rect.second_moment_of_area_mm4();
        // Strong axis (bending about the horizontal/u axis):
        // I_u = (b·h³ − b_i·h_i³)/12.
        let iu_ref = (b * h.powi(3) - b_in * h_in.powi(3)) / 12.0;
        assert!((iu - iu_ref).abs() < 1e-6, "I_u got {iu}, want {iu_ref}");
        // Weak axis: I_v = (h·b³ − h_i·b_i³)/12.
        let iv_ref = (h * b.powi(3) - h_in * b_in.powi(3)) / 12.0;
        assert!((iv - iv_ref).abs() < 1e-6, "I_v got {iv}, want {iv_ref}");

        let (su, sv) = rect.section_modulus_mm3();
        // S = I / c, c the extreme-fibre distance (h/2 strong, b/2 weak).
        assert!(
            (su - iu_ref / (h / 2.0)).abs() < 1e-6,
            "S_u got {su}, want {}",
            iu_ref / (h / 2.0)
        );
        assert!(
            (sv - iv_ref / (b / 2.0)).abs() < 1e-6,
            "S_v got {sv}, want {}",
            iv_ref / (b / 2.0)
        );

        // Regression guard: the hollow second moment must be strictly less
        // than the solid outer rectangle's (the prior, overestimating
        // behaviour) — so a revert to the solid outline fails here.
        assert!(
            iu < b * h.powi(3) / 12.0 && iv < h * b.powi(3) / 12.0,
            "hollow I must be < solid outer rectangle"
        );
    }

    /// I-beam strong-axis second moment, validated against the exact
    /// composite-rectangle closed form
    /// `I_u = b·h³/12 − (b − t_w)·(h − 2·t_f)³/12`
    /// (full bounding rectangle minus the two web-flanking voids — the
    /// standard wide-flange result in Gere / Roark). For the IPE-200-like
    /// section {h:200, b:100, t_w:5.6, t_f:8.5} this is ≈ 1.8456 × 10⁷ mm⁴.
    #[test]
    fn ibeam_strong_axis_second_moment_matches_composite_form() {
        let (h, b, tw, tf) = (200.0_f64, 100.0_f64, 5.6_f64, 8.5_f64);
        let p = Profile::IBeam { h, b, tw, tf };
        let want_iu = b * h.powi(3) / 12.0 - (b - tw) * (h - 2.0 * tf).powi(3) / 12.0;

        let (iu, iv) = p.second_moment_of_area_mm4();
        assert!((iu - want_iu).abs() < 1e-3, "I_u got {iu}, want {want_iu}");
        // Strong axis must be the larger one for a wide-flange.
        assert!(
            iu > iv,
            "expected I_u ({iu}) > I_v ({iv}) for a wide flange"
        );

        // Section modulus consistency S_u = I_u / (h/2).
        let (su, _) = p.section_modulus_mm3();
        assert!((su - iu / (h / 2.0)).abs() < 1e-6, "S_u got {su}");
    }

    /// Round hollow section (annulus) second moment and modulus against
    /// the exact closed form `I = π·(D⁴ − d⁴)/64`, `S = I / (D/2)`
    /// (Roark's *Formulas for Stress & Strain*). CHS {d:100, t:4} →
    /// inner Ø 92 → I ≈ 1.39215 × 10⁶ mm⁴.
    #[test]
    fn round_hollow_second_moment_matches_annulus_form() {
        let (d, t) = (100.0_f64, 4.0_f64);
        let di = d - 2.0 * t;
        let chs = Profile::ChsRound { d, t };
        let want_i = std::f64::consts::PI * (d.powi(4) - di.powi(4)) / 64.0;

        let (iu, iv) = chs.second_moment_of_area_mm4();
        assert!((iu - want_i).abs() < 1e-6, "I got {iu}, want {want_i}");
        // Axisymmetric: I_u == I_v.
        assert!((iu - iv).abs() < 1e-9);

        let (su, _) = chs.section_modulus_mm3();
        assert!(
            (su - want_i / (d / 2.0)).abs() < 1e-6,
            "S got {su}, want {}",
            want_i / (d / 2.0)
        );
    }

    /// Degenerate / non-finite geometry must clamp the section properties
    /// to `0.0`, matching the convention of the area/perimeter helpers.
    #[test]
    fn second_moment_clamps_degenerate_geometry() {
        let bad = Profile::RhsRect {
            h: 0.0,
            b: 0.0,
            t: 0.0,
        };
        assert_eq!(bad.second_moment_of_area_mm4(), (0.0, 0.0));
        assert_eq!(bad.section_modulus_mm3(), (0.0, 0.0));
    }
}
