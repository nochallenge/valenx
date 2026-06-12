//! 2D structure layout coordinates (naview / Forna-class).
//!
//! To *draw* an RNA secondary structure every base needs an `(x, y)`
//! position. This module computes such a layout: helices are drawn as
//! straight ladders and loops as circles, the standard "squiggle
//! plot" produced by naview, RNAplot and Forna.
//!
//! ## Method
//!
//! A recursive radial layout:
//!
//! 1. The exterior loop is laid out left-to-right along a baseline.
//! 2. Each base pair `(i, j)` opens a helix — the two strands are
//!    drawn as parallel ladders advancing by a fixed bond length.
//! 3. The loop closed by a pair is drawn as a circle whose radius is
//!    chosen so its `m` boundary elements (unpaired bases + the stubs
//!    of enclosed helices) are evenly spaced; the enclosed helices
//!    then recurse, each rooted on the circle and pointing radially
//!    outward.
//!
//! The result is a [`Layout`] — one point per base plus the helix
//! connectivity — that a rendering layer turns into an SVG / canvas
//! drawing. Coordinates are in arbitrary drawing units; the bond
//! length is 1.0.
//!
//! ## v1 scope
//!
//! This produces a clean, non-self-intersecting layout for the common
//! nested structures. It does not run the force-directed
//! overlap-removal pass that naview/Forna apply for very large or
//! pathological structures, and pseudoknotted pairs are drawn as
//! straight chords (their crossing is shown, not resolved). Both
//! limitations are stated plainly.

use crate::structure::Structure;
use std::f64::consts::PI;

/// The fixed distance between consecutive bases / between the two
/// strands of a helix, in drawing units.
pub const BOND_LENGTH: f64 = 1.0;

/// A 2-D point in drawing-unit space.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Point {
    /// Horizontal coordinate.
    pub x: f64,
    /// Vertical coordinate.
    pub y: f64,
}

impl Point {
    /// A point at `(x, y)`.
    pub fn new(x: f64, y: f64) -> Self {
        Point { x, y }
    }
}

/// A computed 2-D layout of a secondary structure.
#[derive(Clone, Debug, PartialEq)]
pub struct Layout {
    /// One position per base, indexed by sequence position.
    pub points: Vec<Point>,
    /// The base pairs (so a renderer can draw the pairing bonds).
    pub pairs: Vec<(usize, usize)>,
}

impl Layout {
    /// The number of bases laid out.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// `true` if the layout has no bases.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// The axis-aligned bounding box of the layout as
    /// `(min_x, min_y, max_x, max_y)`. Returns all-zero for an empty
    /// layout.
    pub fn bounding_box(&self) -> (f64, f64, f64, f64) {
        if self.points.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for p in &self.points {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
        (min_x, min_y, max_x, max_y)
    }
}

/// Computes a 2-D layout of `structure` (naview / Forna-class).
///
/// The layout never fails; a pseudoknotted structure is laid out from
/// its nested skeleton with the crossing pairs drawn as chords.
pub fn layout(structure: &Structure) -> Layout {
    let n = structure.len();
    let mut points = vec![Point::new(0.0, 0.0); n];
    if n == 0 {
        return Layout {
            points,
            pairs: Vec::new(),
        };
    }

    // Lay out the exterior loop along the x-axis, recursing into each
    // helix that springs from it.
    let mut cursor_x = 0.0;
    let mut k = 0;
    while k < n {
        match structure.partner(k) {
            Some(p) if p > k => {
                // A helix from k to p springs upward from the
                // baseline at x = cursor_x.
                let base_5 = Point::new(cursor_x, 0.0);
                let base_3 = Point::new(cursor_x + BOND_LENGTH, 0.0);
                place_helix(
                    structure,
                    k,
                    p,
                    base_5,
                    base_3,
                    0.0, // helix grows in +y
                    1.0,
                    &mut points,
                );
                cursor_x += BOND_LENGTH * 2.0;
                k = p + 1;
            }
            _ => {
                // An exterior unpaired base.
                points[k] = Point::new(cursor_x, 0.0);
                cursor_x += BOND_LENGTH;
                k += 1;
            }
        }
    }

    Layout {
        points,
        pairs: structure
            .pairs()
            .into_iter()
            .map(|bp| (bp.i, bp.j))
            .collect(),
    }
}

/// Places the helix closed by pair `(i, j)` whose two closing bases
/// sit at `base_5` (for `i`) and `base_3` (for `j`). `dir_x`/`dir_y`
/// is the unit direction the helix advances in.
#[allow(clippy::too_many_arguments)]
fn place_helix(
    s: &Structure,
    i: usize,
    j: usize,
    base_5: Point,
    base_3: Point,
    dir_x: f64,
    dir_y: f64,
    points: &mut [Point],
) {
    points[i] = base_5;
    points[j] = base_3;

    // Walk inward while the helix keeps stacking (i+1 pairs j-1, ...).
    let mut a = i;
    let mut b = j;
    let mut p5 = base_5;
    let mut p3 = base_3;
    loop {
        let next_a = a + 1;
        if next_a >= b {
            break;
        }
        let next_b = b - 1;
        if s.partner(next_a) == Some(next_b) {
            // helix continues: advance both strands by one bond length
            p5 = Point::new(p5.x + dir_x * BOND_LENGTH, p5.y + dir_y * BOND_LENGTH);
            p3 = Point::new(p3.x + dir_x * BOND_LENGTH, p3.y + dir_y * BOND_LENGTH);
            points[next_a] = p5;
            points[next_b] = p3;
            a = next_a;
            b = next_b;
        } else {
            break;
        }
    }

    // The loop closed by the innermost stacked pair (a, b): place its
    // elements on a circle.
    place_loop(s, a, b, p5, p3, dir_x, dir_y, points);
}

/// Places the loop closed by pair `(i, j)` (whose closing bases are
/// already at `p5`/`p3`) by spreading its boundary elements around a
/// circle. Enclosed helices recurse outward from the circle.
#[allow(clippy::too_many_arguments)]
fn place_loop(
    s: &Structure,
    i: usize,
    j: usize,
    p5: Point,
    p3: Point,
    dir_x: f64,
    dir_y: f64,
    points: &mut [Point],
) {
    // Collect the loop's boundary elements in 5'->3' order: either an
    // unpaired base or the 5' base of an enclosed helix.
    enum Elem {
        Unpaired(usize),
        Helix(usize, usize),
    }
    let mut elems: Vec<Elem> = Vec::new();
    let mut k = i + 1;
    while k < j {
        match s.partner(k) {
            Some(p) if p > k && p < j => {
                elems.push(Elem::Helix(k, p));
                k = p + 1;
            }
            _ => {
                elems.push(Elem::Unpaired(k));
                k += 1;
            }
        }
    }
    if elems.is_empty() {
        return; // a 0-size loop (shouldn't happen for a valid hairpin)
    }

    // The circle passes through both closing bases and carries
    // `m = elems.len() + 1` arcs (the +1 is the closing-pair arc).
    let m = elems.len() + 1;
    // Radius so consecutive boundary points are ~BOND_LENGTH apart.
    let radius = (BOND_LENGTH / (2.0 * (PI / m as f64).sin())).max(BOND_LENGTH);

    // Circle centre: offset from the midpoint of the closing pair,
    // along the helix direction, by the circle "sagitta".
    let mid = Point::new(0.5 * (p5.x + p3.x), 0.5 * (p5.y + p3.y));
    let half_chord = 0.5 * BOND_LENGTH;
    let sagitta = (radius * radius - half_chord * half_chord).max(0.0).sqrt();
    let centre = Point::new(mid.x + dir_x * sagitta, mid.y + dir_y * sagitta);

    // Angles: the closing pair occupies one arc. Place p5 and p3 at
    // their actual angles around the centre, then distribute the
    // boundary elements on the remaining arc.
    let ang_5 = (p5.y - centre.y).atan2(p5.x - centre.x);
    let ang_3 = (p3.y - centre.y).atan2(p3.x - centre.x);
    // Sweep from ang_5 to ang_3 the "long way" (through the loop).
    let mut sweep = ang_3 - ang_5;
    // Normalise so the sweep goes around the loop interior.
    while sweep <= 0.0 {
        sweep += 2.0 * PI;
    }
    let step = sweep / m as f64;

    for (idx, e) in elems.iter().enumerate() {
        let ang = ang_5 + step * (idx as f64 + 1.0);
        let pt = Point::new(centre.x + radius * ang.cos(), centre.y + radius * ang.sin());
        match e {
            Elem::Unpaired(b) => {
                points[*b] = pt;
            }
            Elem::Helix(h5, h3) => {
                // The helix springs radially outward from the circle.
                let out_x = ang.cos();
                let out_y = ang.sin();
                // Two closing bases of the sub-helix sit a bond-length
                // apart, tangent to the circle.
                let tan_x = -ang.sin();
                let tan_y = ang.cos();
                let b5 = Point::new(
                    pt.x - tan_x * 0.5 * BOND_LENGTH,
                    pt.y - tan_y * 0.5 * BOND_LENGTH,
                );
                let b3 = Point::new(
                    pt.x + tan_x * 0.5 * BOND_LENGTH,
                    pt.y + tan_y * 0.5 * BOND_LENGTH,
                );
                place_helix(s, *h5, *h3, b5, b3, out_x, out_y, points);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_structure_lays_out_empty() {
        let l = layout(&Structure::empty(0));
        assert!(l.is_empty());
        assert_eq!(l.bounding_box(), (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn unpaired_chain_is_a_straight_line() {
        let l = layout(&Structure::empty(5));
        assert_eq!(l.len(), 5);
        // all on the baseline y == 0
        for p in &l.points {
            assert!((p.y).abs() < 1e-9, "exterior bases should be on y=0");
        }
        // x increases monotonically
        for w in l.points.windows(2) {
            assert!(w[1].x > w[0].x);
        }
    }

    #[test]
    fn hairpin_has_distinct_finite_points() {
        let s = Structure::from_dot_bracket("(((....)))").unwrap();
        let l = layout(&s);
        assert_eq!(l.len(), 10);
        // every coordinate is finite
        for p in &l.points {
            assert!(p.x.is_finite() && p.y.is_finite());
        }
        // the loop bases are lifted off the baseline
        let (_, _, _, max_y) = l.bounding_box();
        assert!(max_y > 0.5, "the hairpin loop should rise above y=0");
    }

    #[test]
    fn paired_bases_are_roughly_a_bond_length_apart() {
        let s = Structure::from_dot_bracket("(((....)))").unwrap();
        let l = layout(&s);
        for &(i, j) in &l.pairs {
            let dx = l.points[i].x - l.points[j].x;
            let dy = l.points[i].y - l.points[j].y;
            let d = (dx * dx + dy * dy).sqrt();
            assert!(
                (d - BOND_LENGTH).abs() < 0.5,
                "pair ({i},{j}) spans {d}, expected ~{BOND_LENGTH}"
            );
        }
    }

    #[test]
    fn pairs_are_reported() {
        let s = Structure::from_dot_bracket("(((...)))").unwrap();
        let l = layout(&s);
        assert_eq!(l.pairs.len(), 3);
    }

    #[test]
    fn multiloop_lays_out_without_nan() {
        let s = Structure::from_dot_bracket("((((....))((....))))").unwrap();
        let l = layout(&s);
        for p in &l.points {
            assert!(
                p.x.is_finite() && p.y.is_finite(),
                "NaN in multiloop layout"
            );
        }
    }

    #[test]
    fn two_hairpins_do_not_collapse_to_a_point() {
        let s = Structure::from_dot_bracket("(((...)))(((...)))").unwrap();
        let l = layout(&s);
        let (min_x, _, max_x, _) = l.bounding_box();
        assert!(max_x - min_x > 1.0, "layout should have horizontal extent");
    }
}
