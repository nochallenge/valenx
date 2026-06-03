//! **Feature 4 — loop modelling (CCD loop closure).**
//!
//! The backbone-transfer step leaves gaps wherever the target has
//! residues the template does not — insertions, and the regions a
//! partial template does not cover. These must be built *de novo* and
//! made to connect smoothly to the anchored backbone on either side.
//!
//! This module builds loops by **cyclic coordinate descent (CCD)** —
//! the loop-closure algorithm from robotics adapted to protein
//! backbones by Canutescu & Dunbrack (2003). A loop is grown from one
//! anchor with idealised geometry, leaving a *closure gap* to the
//! other anchor. CCD then iterates: for each adjustable backbone
//! dihedral, in turn, rotate it by the angle that minimises the
//! distance between the loop's free end and the fixed target anchor.
//! The loop "reaches" toward the anchor and, on convergence, closes.
//!
//! The CCD move and the convergence test are the genuine published
//! algorithm. A production loop modeller additionally biases the
//! dihedral sampling toward the Ramachandran-allowed region and
//! screens for clashes; this v1 closes the loop and reports whether
//! the closure tolerance was met.

use nalgebra::{Point3, Rotation3, Unit, Vector3};

use crate::error::{Result, StructPredictError};
use crate::model::{ideal, place_atom, ModelResidue, ProteinModel};

/// The outcome of closing one loop.
#[derive(Clone, Debug, PartialEq)]
pub struct LoopClosure {
    /// Half-open `[start, end)` residue range that was rebuilt.
    pub range: (usize, usize),
    /// Final distance from the loop's free end to the fixed anchor
    /// (ångström). Below `tolerance` means the loop closed.
    pub gap: f64,
    /// CCD iterations performed.
    pub iterations: usize,
    /// Whether the closure tolerance was reached.
    pub closed: bool,
}

/// Builds an initial extended-geometry backbone for a loop region.
///
/// Residues `[start, end)` are grown from the residue *before*
/// `start` (the N-terminal anchor) with idealised bond lengths and
/// angles and φ = ψ = 180° (a fully extended chain). The free end
/// will not yet meet the C-terminal anchor — that is what CCD then
/// fixes.
fn build_extended_loop(model: &mut ProteinModel, start: usize, end: usize) -> Result<()> {
    if start == 0 {
        // No N-anchor: seed the very first three atoms at canonical
        // positions so the loop has a frame to grow from.
        let r = &mut model.residues[0];
        r.n = Some(Point3::new(0.0, 0.0, 0.0));
        r.ca = Some(Point3::new(ideal::N_CA, 0.0, 0.0));
        r.c = Some(place_atom(
            &Point3::new(-1.0, 1.0, 0.0),
            &r.n.unwrap(),
            &r.ca.unwrap(),
            ideal::CA_C,
            ideal::N_CA_C.to_radians(),
            std::f64::consts::PI,
        ));
    }
    for i in start..end {
        // Anchor atoms come from residue i-1 (always built by now).
        let prev = if i == 0 {
            // residue 0 was seeded above
            model.residues[0].clone()
        } else {
            model.residues[i - 1].clone()
        };
        let (pn, pca, pc) = (
            prev.n.ok_or_else(|| anchor_err(i))?,
            prev.ca.ok_or_else(|| anchor_err(i))?,
            prev.c.ok_or_else(|| anchor_err(i))?,
        );
        // N of residue i: from prev N-CA-C, ψ=180.
        let n = place_atom(
            &pn,
            &pca,
            &pc,
            ideal::C_N,
            ideal::CA_C_N.to_radians(),
            std::f64::consts::PI,
        );
        // CA of residue i: from prev CA-C-N(i), ω=180.
        let ca = place_atom(
            &pca,
            &pc,
            &n,
            ideal::N_CA,
            ideal::C_N_CA.to_radians(),
            std::f64::consts::PI,
        );
        // C of residue i: from prev C-N(i)-CA(i), φ=180.
        let c = place_atom(
            &pc,
            &n,
            &ca,
            ideal::CA_C,
            ideal::N_CA_C.to_radians(),
            std::f64::consts::PI,
        );
        // O of residue i: off the C, in the peptide plane.
        let o = place_atom(
            &n,
            &ca,
            &c,
            ideal::C_O,
            (180.0 - ideal::CA_C_N).to_radians(),
            0.0,
        );
        let r = &mut model.residues[i];
        r.n = Some(n);
        r.ca = Some(ca);
        r.c = Some(c);
        r.o = Some(o);
    }
    Ok(())
}

fn anchor_err(i: usize) -> StructPredictError {
    StructPredictError::invalid(
        "loop_anchor",
        format!("residue {i} before the loop has no backbone to grow from"),
    )
}

/// Closes a single loop region by cyclic coordinate descent.
///
/// `model` must already have residues `[start, end)` set up as a gap
/// and the residue at `end` (the C-terminal anchor) built with its
/// backbone. The loop is grown extended, then CCD rotates each loop
/// φ/ψ dihedral to drive the loop's free end (the `N` atom of residue
/// `end`, recomputed from the loop) onto the fixed anchor.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a degenerate range, a missing
/// C-anchor, or a missing N-anchor when `start > 0`.
pub fn close_loop(
    model: &mut ProteinModel,
    start: usize,
    end: usize,
    max_iterations: usize,
    tolerance: f64,
) -> Result<LoopClosure> {
    if end <= start || end > model.residues.len() {
        return Err(StructPredictError::invalid(
            "loop_range",
            format!("[{start}, {end}) is not a valid loop"),
        ));
    }
    if end >= model.residues.len() {
        return Err(StructPredictError::invalid(
            "loop_range",
            "loop has no C-terminal anchor residue",
        ));
    }
    if start > 0 && !model.residues[start - 1].has_backbone() {
        return Err(anchor_err(start));
    }
    let anchor = model.residues[end]
        .ca
        .ok_or_else(|| StructPredictError::invalid("loop_anchor", "C-anchor has no Cα"))?;

    build_extended_loop(model, start, end)?;

    // CCD adjustable dihedrals: the φ (about N-CA) and ψ (about CA-C)
    // bonds of every loop residue. Each is an axis through two atoms;
    // rotating it pivots the downstream chain.
    let mut iterations = 0;
    let mut gap = loop_gap(model, end, anchor);
    for _ in 0..max_iterations {
        iterations += 1;
        if gap < tolerance {
            break;
        }
        for i in start..end {
            // φ: rotate about the N(i)→CA(i) axis.
            ccd_rotate_about(model, i, end, BondAxis::Phi, anchor);
            // ψ: rotate about the CA(i)→C(i) axis.
            ccd_rotate_about(model, i, end, BondAxis::Psi, anchor);
        }
        gap = loop_gap(model, end, anchor);
    }
    let closed = gap < tolerance;
    Ok(LoopClosure {
        range: (start, end),
        gap,
        iterations,
        closed,
    })
}

/// Which backbone bond a CCD rotation pivots about.
#[derive(Copy, Clone)]
enum BondAxis {
    /// The N→Cα bond (the φ dihedral).
    Phi,
    /// The Cα→C bond (the ψ dihedral).
    Psi,
}

/// The current closure gap: distance from the loop's recomputed free
/// end (the N atom of the residue *just past* the loop, grown from
/// the loop's last residue) to the fixed anchor Cα.
fn loop_gap(model: &ProteinModel, end: usize, anchor: Point3<f64>) -> f64 {
    // The free end is the C atom of the loop's last residue; closing
    // means that C connects properly to the anchor. We use the last
    // loop C → anchor distance minus the ideal C–N + N–CA span as the
    // residual, but the simplest robust target is the loop's terminal
    // Cα reaching the anchor minus one virtual Cα–Cα step.
    let last = end - 1;
    if let Some(c) = model.residues[last].c {
        // Ideal: the loop's final C sits ~1.33 Å (C–N) from the
        // anchor's N, which is ~1.46 Å from the anchor Cα. Use the
        // straight C→anchorCα target distance.
        let target = ideal::C_N + ideal::N_CA;
        ((c - anchor).norm() - target).abs()
    } else {
        f64::INFINITY
    }
}

/// One CCD move: rotate the loop's downstream atoms about residue
/// `i`'s `axis` bond by the angle that best drives the moving end
/// onto the anchor.
fn ccd_rotate_about(
    model: &mut ProteinModel,
    i: usize,
    end: usize,
    axis: BondAxis,
    anchor: Point3<f64>,
) {
    let res = model.residues[i].clone();
    let (Some(n), Some(ca), Some(c)) = (res.n, res.ca, res.c) else {
        return;
    };
    let (origin, dir) = match axis {
        BondAxis::Phi => (n, (ca - n)),
        BondAxis::Psi => (ca, (c - ca)),
    };
    let dir = match Unit::try_new(dir, 1e-9) {
        Some(u) => u,
        None => return,
    };
    // The moving end: the C atom of the loop's last residue.
    let moving = match model.residues[end - 1].c {
        Some(p) => p,
        None => return,
    };
    // CCD closed-form optimal angle (Canutescu-Dunbrack): project the
    // moving point and the anchor into the plane perpendicular to the
    // rotation axis and rotate to align them.
    let r_m = perp_component(moving - origin, &dir);
    let r_t = perp_component(anchor - origin, &dir);
    if r_m.norm() < 1e-6 || r_t.norm() < 1e-6 {
        return;
    }
    let s = r_m.normalize();
    let t_hat = r_t.normalize();
    // signed angle from s to t_hat about dir
    let cos = s.dot(&t_hat).clamp(-1.0, 1.0);
    let sin = dir.dot(&s.cross(&t_hat));
    let theta = sin.atan2(cos);
    if theta.abs() < 1e-6 {
        return;
    }
    let rot = Rotation3::from_axis_angle(&dir, theta);
    // Apply to all downstream atoms: the rest of residue i past the
    // axis, plus every atom of residues i+1..end.
    rotate_downstream(model, i, end, axis, &origin, &rot);
}

/// The component of `v` perpendicular to the unit axis `dir`.
fn perp_component(v: Vector3<f64>, dir: &Unit<Vector3<f64>>) -> Vector3<f64> {
    v - dir.into_inner() * v.dot(dir)
}

/// Rotates every loop atom downstream of residue `i`'s `axis` bond.
fn rotate_downstream(
    model: &mut ProteinModel,
    i: usize,
    end: usize,
    axis: BondAxis,
    origin: &Point3<f64>,
    rot: &Rotation3<f64>,
) {
    let apply = |p: &mut Option<Point3<f64>>| {
        if let Some(point) = p {
            *point = origin + rot * (*point - origin);
        }
    };
    // Atoms of residue i past the rotation axis.
    {
        let r = &mut model.residues[i];
        match axis {
            BondAxis::Phi => {
                // φ axis is N→CA: C and O move (CA stays on the axis).
                apply(&mut r.c);
                apply(&mut r.o);
                apply(&mut r.cb);
            }
            BondAxis::Psi => {
                // ψ axis is CA→C: only O moves within residue i.
                apply(&mut r.o);
            }
        }
    }
    for j in (i + 1)..end {
        let r = &mut model.residues[j];
        apply(&mut r.n);
        apply(&mut r.ca);
        apply(&mut r.c);
        apply(&mut r.o);
        apply(&mut r.cb);
    }
}

/// Closes every gap in a model by CCD.
///
/// Runs [`close_loop`] on each gap reported by
/// [`ProteinModel::gaps`]. A gap that runs to the very end of the
/// chain (no C-terminal anchor) is grown extended but not
/// CCD-closed — there is nothing to close to — and reported with
/// `closed = false`.
///
/// # Errors
/// Propagates [`close_loop`] errors.
pub fn model_loops(
    model: &mut ProteinModel,
    max_iterations: usize,
    tolerance: f64,
) -> Result<Vec<LoopClosure>> {
    let gaps = model.gaps();
    let mut out = Vec::new();
    for (start, end) in gaps {
        if end >= model.residues.len() {
            // C-terminal tail: grow extended, no closure target.
            build_extended_loop(model, start, end)?;
            out.push(LoopClosure {
                range: (start, end),
                gap: f64::NAN,
                iterations: 0,
                closed: false,
            });
        } else {
            out.push(close_loop(model, start, end, max_iterations, tolerance)?);
        }
    }
    Ok(out)
}

/// Rebuilds an idealised carbonyl `O` for a residue whose backbone
/// `N`/`CA`/`C` are set — used after CCD perturbs a backbone.
pub fn rebuild_carbonyl_oxygen(res: &mut ModelResidue) {
    if let (Some(n), Some(ca), Some(c)) = (res.n, res.ca, res.c) {
        res.o = Some(place_atom(
            &n,
            &ca,
            &c,
            ideal::C_O,
            (180.0 - ideal::CA_C_N).to_radians(),
            0.0,
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A model with a built N-anchor, a 4-residue gap, and a built
    /// C-anchor placed within reach.
    fn loop_test_model() -> ProteinModel {
        let mut m = ProteinModel::from_sequence("AAAAAAA").expect("model");
        // N-anchor: residue 0 with a sane backbone.
        {
            let r = &mut m.residues[0];
            r.n = Some(Point3::new(0.0, 0.0, 0.0));
            r.ca = Some(Point3::new(1.458, 0.0, 0.0));
            r.c = Some(Point3::new(2.0, 1.4, 0.0));
            r.o = Some(Point3::new(1.3, 2.4, 0.0));
        }
        // C-anchor: residue 6, ~13 Å away — reachable by a 5-residue
        // loop (residues 1..6).
        {
            let r = &mut m.residues[6];
            r.n = Some(Point3::new(12.0, 3.0, 0.0));
            r.ca = Some(Point3::new(13.0, 3.5, 0.0));
            r.c = Some(Point3::new(14.0, 2.8, 0.0));
            r.o = Some(Point3::new(14.5, 1.8, 0.0));
        }
        m
    }

    #[test]
    fn ccd_reduces_the_closure_gap() {
        let mut m = loop_test_model();
        let anchor = m.residues[6].ca.unwrap();
        // Gap before closure (extended loop).
        let mut probe = m.clone();
        build_extended_loop(&mut probe, 1, 6).expect("extend");
        let before = loop_gap(&probe, 6, anchor);
        // Closure.
        let res = close_loop(&mut m, 1, 6, 500, 0.3).expect("close");
        assert!(
            res.gap <= before + 1e-6,
            "gap did not shrink: {before} -> {}",
            res.gap
        );
        assert_eq!(res.range, (1, 6));
    }

    #[test]
    fn extended_loop_has_sane_bond_lengths() {
        let mut m = loop_test_model();
        build_extended_loop(&mut m, 1, 6).expect("extend");
        // The N–CA bond of a built loop residue is ~1.458 Å.
        let r = &m.residues[3];
        let d = (r.ca.unwrap() - r.n.unwrap()).norm();
        assert!((d - ideal::N_CA).abs() < 1e-6, "N-CA = {d}");
    }

    #[test]
    fn degenerate_range_rejected() {
        let mut m = loop_test_model();
        assert!(close_loop(&mut m, 3, 3, 10, 0.1).is_err());
    }

    #[test]
    fn model_loops_handles_c_terminal_tail() {
        let mut m = ProteinModel::from_sequence("AAAAA").expect("model");
        // Only residue 0 built — residues 1..5 are a C-terminal tail.
        let r = &mut m.residues[0];
        r.n = Some(Point3::new(0.0, 0.0, 0.0));
        r.ca = Some(Point3::new(1.458, 0.0, 0.0));
        r.c = Some(Point3::new(2.0, 1.4, 0.0));
        r.o = Some(Point3::new(1.3, 2.4, 0.0));
        let closures = model_loops(&mut m, 100, 0.3).expect("loops");
        assert_eq!(closures.len(), 1);
        assert!(!closures[0].closed, "tail is not closed");
        assert!(m.is_complete(), "tail still gets a backbone");
    }
}
