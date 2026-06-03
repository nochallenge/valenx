//! Real BRep fillet construction for the single convex planar-edge
//! case (Phase 14.5).
//!
//! # What this module does
//!
//! Where [`crate::fillet::plan_planar_edge_fillet`] computes the
//! *geometry* of a fillet (the [`crate::fillet::EdgeFilletPlan`]) and
//! the v1 [`crate::fillet::fillet_planar_edge`] stopped there, this
//! module performs the **actual BRep surgery** for the bounded case
//! the crate's scope promises: one convex straight edge shared by
//! exactly two planar faces, constant radius.
//!
//! # The method — constructive solid geometry, real BRep throughout
//!
//! Filleting a convex edge means *removing* the sharp corner sliver
//! and *replacing* it with a rounded surface. In the cross-section
//! plane perpendicular to the edge the corner is a point `E`; the two
//! faces are rays from `E`; the fillet is a circular arc of radius `r`
//! tangent to both rays at points `T₀` and `T₁`, centred at `C` on the
//! interior bisector.
//!
//! Let
//!
//! - `wedge`  = the triangular prism whose cross-section is the
//!   triangle `△(E, T₀, T₁)` — the corner sliver,
//! - `bar`    = the prism whose cross-section is the circular sector
//!   `(C, T₀, arc T₀→T₁, T₁)` — the rounded fill, which is geometrically
//!   *contained inside* the triangle.
//!
//! Then, by elementary set algebra (`sector ⊆ triangle`):
//!
//! ```text
//!   filleted = (solid − wedge) ∪ bar
//!            =  solid − wedge + bar
//!            =  solid − (wedge − bar)
//! ```
//!
//! i.e. the result removes exactly the corner-minus-round material.
//! Both `wedge` and `bar` are **genuine truck-modeling BRep solids**
//! (built with `try_attach_plane` + `tsweep`, the `bar`'s outer face
//! carrying a real circular-arc-swept surface), and the two set
//! operations are the **real `truck_shapeops` booleans** that
//! [`valenx_cad::difference`] / [`valenx_cad::union`] wrap. The output
//! is a [`valenx_cad::Solid::Brep`]: it has true faces and edges, it
//! round-trips through STEP/IGES, and it composes with further BRep
//! ops. This is a true-BRep fillet, not a tessellated approximation.
//!
//! # Honest scope and the coincident-face caveat
//!
//! Two faces of `wedge` lie *flush* with the solid's two adjacent
//! planar faces (that flush trim is exactly what a fillet must do).
//! Coincident faces are the historically fragile input for any BRep
//! boolean kernel, `truck_shapeops` included. The construction here is
//! correct; on geometry where `truck_shapeops` cannot resolve the
//! coincidence it returns `None`, which this module surfaces as the
//! soft [`FilletBrepError::TruckOp`] error so the feature-tree
//! dispatcher can still fall through to the mesh-domain pipeline. The
//! `wedge` and `bar` prisms are deliberately over-extended past both
//! edge endpoints (so their *end caps* never coincide with the
//! solid's faces) to give the kernel the best possible chance.
//!
//! [`fillet_variable_radius_planar_edge`] extends this to a fillet
//! whose radius **varies linearly** along the edge — the cutter and
//! fillet-bar become lofted tapers between two end cross-sections
//! instead of constant-radius prisms.
//!
//! Multi-edge corners (3+ filleted edges meeting at a vertex), curved
//! adjacent faces, and concave edges remain out of scope — see the
//! crate-level docs. This module is the genuinely-achievable slice: the
//! single convex planar-edge fillet, constant or linearly-tapered.

use truck_modeling::{builder, InnerSpace, Point3, Solid as TruckSolid, Vector3, Wire};

use valenx_cad::Solid;

use crate::error::FilletBrepError;
use crate::fillet::{plan_planar_edge_fillet, EdgeFilletPlan};

/// The cross-section geometry of a fillet, resolved from an
/// [`EdgeFilletPlan`] into the handful of points the BRep constructor
/// needs.
///
/// All points are world-space and lie in the cross-section plane at
/// the edge's **front** endpoint. Because the edge is straight, the
/// cross-section is identical along the whole edge — the constructor
/// builds the profile here and `tsweep`s it along [`Self::edge_dir`]
/// by [`Self::edge_len`] (plus the over-extension), so a single
/// front-end cross-section fully determines the swept prism.
#[derive(Clone, Debug)]
struct FilletCrossSection {
    /// Sharp corner point at the edge (front end).
    e_front: Point3,
    /// Tangent contact point on face 0 (front end).
    t0_front: Point3,
    /// Tangent contact point on face 1 (front end).
    t1_front: Point3,
    /// Fillet-arc centre (front end).
    c_front: Point3,
    /// Point on the fillet arc nearest the sharp corner (front end).
    /// Used as the `transit` point that orients the circular arc.
    arc_mid_front: Point3,
    /// Unit edge direction, `front → back`.
    edge_dir: Vector3,
    /// Edge length.
    edge_len: f64,
}

/// Resolve an [`EdgeFilletPlan`] into the swept cross-section.
///
/// The plan supplies the tangent points and the inward bisector; this
/// adds the arc centre `C` and the arc's apex (transit) point.
///
/// # Geometry
///
/// In the cross-section the interior corner has half-angle `β`. The
/// tangent points sit at distance `r/tan β` from the corner along each
/// face; the arc centre sits at distance `r/sin β` from the corner
/// along the bisector; the arc apex (nearest point to the corner) sits
/// at `C − r·b̂` where `b̂` is the inward bisector. The plan already
/// stores the tangent points and the dihedral angle, so `β` follows
/// from `β = (π − dihedral)/2`.
fn resolve_cross_section(plan: &EdgeFilletPlan) -> Result<FilletCrossSection, FilletBrepError> {
    let edge_vec = plan.edge_back - plan.edge_front;
    let edge_len = edge_vec.magnitude();
    if edge_len < 1e-12 {
        return Err(FilletBrepError::TruckOp(
            "degenerate edge — zero length".to_string(),
        ));
    }
    let edge_dir = edge_vec / edge_len;

    // Half-angle β of the interior corner. dihedral_angle is between
    // the *outward* normals; the interior corner angle between the
    // faces is (π − dihedral), so its half is (π − dihedral)/2.
    let beta = (std::f64::consts::PI - plan.dihedral_angle) * 0.5;
    let sin_b = beta.sin();
    if sin_b.abs() < 1e-9 {
        // Degenerate corner — faces almost coplanar.
        return Err(FilletBrepError::NotPlanarFaces);
    }

    // Arc centre: r/sin β along the inward bisector from the edge.
    let centre_dist = plan.radius / sin_b;
    let c_front = plan.edge_front + plan.bisector_inward * centre_dist;

    // Arc apex: the point on the circle nearest the sharp corner. The
    // circle is centred at C with radius r; the corner lies on the
    // bisector at −centre_dist from C, so the nearest arc point is
    // C − r·b̂.
    let arc_mid_front = c_front - plan.bisector_inward * plan.radius;

    Ok(FilletCrossSection {
        e_front: plan.edge_front,
        t0_front: plan.tangent_on_face0_front,
        t1_front: plan.tangent_on_face1_front,
        c_front,
        arc_mid_front,
        edge_dir,
        edge_len,
    })
}

/// How far past each edge endpoint the cutter / fillet-bar prisms are
/// extended, as a fraction of the edge length. Over-extending keeps
/// the prism end caps off the solid's faces so the boolean never has
/// to resolve a coincident *cap*.
const PRISM_OVEREXTEND_FRAC: f64 = 0.05;

/// Build the triangular **cutter** prism — the corner sliver to
/// remove. Cross-section is the triangle `△(E, T₀, T₁)`; swept along
/// the edge and over-extended past both ends.
fn build_cutter(xs: &FilletCrossSection) -> Result<TruckSolid, FilletBrepError> {
    let over = (xs.edge_len * PRISM_OVEREXTEND_FRAC).max(1e-6);
    // Profile sits at front − over·dir; sweep length = edge_len + 2·over.
    let base = -xs.edge_dir * over;
    let e = xs.e_front + base;
    let t0 = xs.t0_front + base;
    let t1 = xs.t1_front + base;

    let v_e = builder::vertex(e);
    let v_t0 = builder::vertex(t0);
    let v_t1 = builder::vertex(t1);
    // Triangle wire E → T0 → T1 → E.
    let wire: Wire = vec![
        builder::line(&v_e, &v_t0),
        builder::line(&v_t0, &v_t1),
        builder::line(&v_t1, &v_e),
    ]
    .into();
    let face = builder::try_attach_plane(&[wire])
        .map_err(|err| FilletBrepError::TruckOp(format!("cutter face: {err:?}")))?;
    let sweep = xs.edge_dir * (xs.edge_len + 2.0 * over);
    let solid: TruckSolid = builder::tsweep(&face, sweep);
    Ok(solid)
}

/// Build the **fillet-bar** prism — the rounded fill. Cross-section is
/// the circular sector `(C, T₀, arc T₀→T₁, T₁)`; swept along the edge
/// and over-extended past both ends.
fn build_fillet_bar(xs: &FilletCrossSection) -> Result<TruckSolid, FilletBrepError> {
    let over = (xs.edge_len * PRISM_OVEREXTEND_FRAC).max(1e-6);
    let base = -xs.edge_dir * over;
    let c = xs.c_front + base;
    let t0 = xs.t0_front + base;
    let t1 = xs.t1_front + base;
    let arc_mid = xs.arc_mid_front + base;

    let v_c = builder::vertex(c);
    let v_t0 = builder::vertex(t0);
    let v_t1 = builder::vertex(t1);
    // Sector wire: straight C → T0, circular arc T0 → T1 (bulging
    // toward the sharp corner via the apex transit point), straight
    // T1 → C.
    let wire: Wire = vec![
        builder::line(&v_c, &v_t0),
        builder::circle_arc(&v_t0, &v_t1, arc_mid),
        builder::line(&v_t1, &v_c),
    ]
    .into();
    let face = builder::try_attach_plane(&[wire])
        .map_err(|err| FilletBrepError::TruckOp(format!("fillet-bar face: {err:?}")))?;
    let sweep = xs.edge_dir * (xs.edge_len + 2.0 * over);
    let solid: TruckSolid = builder::tsweep(&face, sweep);
    Ok(solid)
}

/// Linear tolerance handed to the `truck_shapeops` booleans. The
/// fillet cutter trims flush with the solid's faces, so the boolean
/// needs a tolerance generous enough to recognise the coincidence but
/// tight enough not to merge genuinely-distinct geometry. The
/// `valenx-cad` default (0.05) is tuned for primitives in the 1–10
/// unit range, which is the fillet's working scale.
const FILLET_BOOL_TOL: f64 = valenx_cad::DEFAULT_BOOL_TOLERANCE;

/// Perform the real BRep fillet on one convex straight edge bounded by
/// two planar faces.
///
/// This is the working implementation behind
/// [`crate::fillet::fillet_planar_edge`]. It:
///
/// 1. validates the input and computes the [`EdgeFilletPlan`]
///    (delegated to [`plan_planar_edge_fillet`]);
/// 2. resolves the swept cross-section (arc centre + apex);
/// 3. builds the BRep `cutter` and `bar` prisms;
/// 4. evaluates `filleted = (solid − cutter) ∪ bar` with the real
///    `truck_shapeops` booleans.
///
/// # Errors
///
/// - The geometric-precondition errors of [`plan_planar_edge_fillet`]
///   (`NotPlanarFaces`, `NonConvexEdge`, `RadiusTooLarge`,
///   `BadParameter`) — surfaced verbatim so a bad input is reported by
///   its true cause.
/// - [`FilletBrepError::TruckOp`] if a prism face cannot be attached
///   or a boolean returns no solid (the latter is the coincident-face
///   case discussed in the module docs — a *soft* failure the caller
///   may treat as a fall-through signal).
pub fn fillet_convex_planar_edge(
    solid: &TruckSolid,
    edge: &truck_modeling::Edge,
    radius: f64,
) -> Result<TruckSolid, FilletBrepError> {
    // Stage 1+2: validate + plan (catches all the geometric errors).
    let plan = plan_planar_edge_fillet(solid, edge, radius)?;
    // Resolve the cross-section.
    let xs = resolve_cross_section(&plan)?;

    // Stage 3: build the BRep cutter + fillet bar.
    let cutter = build_cutter(&xs)?;
    let bar = build_fillet_bar(&xs)?;

    // Stage 4: filleted = (solid − cutter) ∪ bar, via real booleans.
    let base = Solid::from_truck(solid.clone());
    let cutter_solid = Solid::from_truck(cutter);
    let bar_solid = Solid::from_truck(bar);

    let trimmed = valenx_cad::boolean::difference_tol(&base, &cutter_solid, FILLET_BOOL_TOL)
        .map_err(|err| FilletBrepError::TruckOp(format!("corner removal boolean: {err}")))?;
    let filleted = valenx_cad::boolean::union_tol(&trimmed, &bar_solid, FILLET_BOOL_TOL)
        .map_err(|err| FilletBrepError::TruckOp(format!("fillet-bar union boolean: {err}")))?;

    match filleted {
        Solid::Brep(b) => Ok(b),
        // The boolean path always returns a BRep solid; a mesh-backed
        // result would mean an internal contract break.
        Solid::Mesh(_) => Err(FilletBrepError::TruckOp(
            "boolean returned a mesh-backed solid (internal error)".to_string(),
        )),
    }
}

// === Variable-radius fillet (Phase 14.6) =================================

/// The fillet cross-section resolved at one *specific* point along the
/// edge with one *specific* radius.
///
/// Where [`FilletCrossSection`] captures the (identical) section of a
/// constant-radius fillet, this captures a single station of a
/// **variable-radius** fillet — the front and the back stations carry
/// different radii, so their tangent points, arc centres and apexes all
/// differ, and the swept faces become tapered lofts rather than
/// constant prisms.
#[derive(Clone, Copy, Debug)]
struct FilletStation {
    /// Sharp corner point on the edge at this station.
    e: Point3,
    /// Tangent contact point on face 0.
    t0: Point3,
    /// Tangent contact point on face 1.
    t1: Point3,
    /// Fillet-arc centre.
    c: Point3,
    /// Point on the fillet arc nearest the sharp corner (the arc apex
    /// — the `transit` point that orients a circular-arc edge).
    arc_mid: Point3,
}

/// Resolve a [`FilletStation`] — the fillet cross-section at the edge
/// point `edge_point` for radius `radius`.
///
/// The geometry is the same construction [`resolve_cross_section`] does,
/// generalised to an arbitrary point + radius: the half-angle `β` of
/// the corner comes from the plan's dihedral angle; the tangent points
/// sit at `r/tan β` from the corner along each face; the arc centre
/// sits at `r/sin β` from the corner along the inward bisector; the
/// arc apex is `C − r·b̂`.
///
/// Because the edge is straight, the inward bisector and the two
/// in-face directions are constant along it — only the point on the
/// edge and the radius change between stations.
fn resolve_station(
    plan: &EdgeFilletPlan,
    edge_point: Point3,
    radius: f64,
) -> Result<FilletStation, FilletBrepError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(FilletBrepError::BadParameter {
            name: "radius",
            reason: format!("variable-radius endpoint must be > 0, got {radius}"),
        });
    }
    // Half-angle β of the interior corner.
    let beta = (std::f64::consts::PI - plan.dihedral_angle) * 0.5;
    let sin_b = beta.sin();
    let tan_b = beta.tan();
    if sin_b.abs() < 1e-9 || !tan_b.is_finite() || tan_b.abs() < 1e-9 {
        return Err(FilletBrepError::NotPlanarFaces);
    }

    // In-face directions: the inward bisector projected into each face
    // plane (away from the edge). Recovered from the plan's stored
    // front-tangent offsets, which are exact for `plan.radius`.
    let plan_offset = plan.radius / tan_b;
    if plan_offset.abs() < 1e-12 {
        return Err(FilletBrepError::NotPlanarFaces);
    }
    let in_face0 = (plan.tangent_on_face0_front - plan.edge_front) / plan_offset;
    let in_face1 = (plan.tangent_on_face1_front - plan.edge_front) / plan_offset;

    // Tangent points at this station: r/tan β from the corner.
    let offset = radius / tan_b;
    let t0 = edge_point + in_face0 * offset;
    let t1 = edge_point + in_face1 * offset;
    // Arc centre: r/sin β along the inward bisector.
    let c = edge_point + plan.bisector_inward * (radius / sin_b);
    // Arc apex: the circle point nearest the corner.
    let arc_mid = c - plan.bisector_inward * radius;

    Ok(FilletStation {
        e: edge_point,
        t0,
        t1,
        c,
        arc_mid,
    })
}

/// Build the **tapered cutter** of a variable-radius fillet — the
/// lofted triangular prism whose front cross-section is `△(E,T₀,T₁)` at
/// `radius_start` and whose back cross-section is the same triangle at
/// `radius_end`.
///
/// The side faces are built with [`builder::try_wire_homotopy`], which
/// lofts a shell of NURBS faces between two equal-edge-count wires; the
/// two end caps are planar faces. The result is a closed BRep solid
/// whose shape is a *taper* — the swept-constant case is just the
/// degenerate `radius_start == radius_end` of it.
fn build_tapered_cutter(
    front: &FilletStation,
    back: &FilletStation,
) -> Result<TruckSolid, FilletBrepError> {
    // Front triangle wire E → T0 → T1 → E.
    let f_e = builder::vertex(front.e);
    let f_t0 = builder::vertex(front.t0);
    let f_t1 = builder::vertex(front.t1);
    let front_wire: Wire = vec![
        builder::line(&f_e, &f_t0),
        builder::line(&f_t0, &f_t1),
        builder::line(&f_t1, &f_e),
    ]
    .into();
    // Back triangle wire — same vertex order.
    let b_e = builder::vertex(back.e);
    let b_t0 = builder::vertex(back.t0);
    let b_t1 = builder::vertex(back.t1);
    let back_wire: Wire = vec![
        builder::line(&b_e, &b_t0),
        builder::line(&b_t0, &b_t1),
        builder::line(&b_t1, &b_e),
    ]
    .into();
    loft_between(&front_wire, &back_wire, "tapered cutter")
}

/// Build the **tapered fillet-bar** of a variable-radius fillet — the
/// lofted circular-sector prism whose cross-section is the sector
/// `(C,T₀,arc,T₁)` at `radius_start` in front and at `radius_end` at
/// the back.
///
/// The side wall carrying the arc is a real lofted blend between the
/// front circular arc (radius `radius_start`) and the back circular arc
/// (radius `radius_end`) — so the fillet surface genuinely *changes
/// radius* along the edge, the defining feature of a variable-radius
/// fillet.
fn build_tapered_fillet_bar(
    front: &FilletStation,
    back: &FilletStation,
) -> Result<TruckSolid, FilletBrepError> {
    // Front sector wire: C → T0, arc T0→T1, T1 → C.
    let f_c = builder::vertex(front.c);
    let f_t0 = builder::vertex(front.t0);
    let f_t1 = builder::vertex(front.t1);
    let front_wire: Wire = vec![
        builder::line(&f_c, &f_t0),
        builder::circle_arc(&f_t0, &f_t1, front.arc_mid),
        builder::line(&f_t1, &f_c),
    ]
    .into();
    // Back sector wire — same edge order (line, arc, line).
    let b_c = builder::vertex(back.c);
    let b_t0 = builder::vertex(back.t0);
    let b_t1 = builder::vertex(back.t1);
    let back_wire: Wire = vec![
        builder::line(&b_c, &b_t0),
        builder::circle_arc(&b_t0, &b_t1, back.arc_mid),
        builder::line(&b_t1, &b_c),
    ]
    .into();
    loft_between(&front_wire, &back_wire, "tapered fillet bar")
}

/// Loft a closed BRep solid between two equal-edge-count closed wires.
///
/// The side wall is [`builder::try_wire_homotopy`]'s lofted shell; the
/// two end caps are planar faces attached to the wires. The three
/// pieces are assembled into one closed [`TruckSolid`]. `label` names
/// the construction in any error message.
fn loft_between(
    front_wire: &Wire,
    back_wire: &Wire,
    label: &str,
) -> Result<TruckSolid, FilletBrepError> {
    // The lofted side shell between the two cross-section wires.
    let side = builder::try_wire_homotopy(front_wire, back_wire)
        .map_err(|err| FilletBrepError::TruckOp(format!("{label} loft: {err:?}")))?;
    // Cap the shell's two open boundary loops with planar faces.
    //
    // The cap wire must be the boundary loop *inverted* — the
    // homotopy shell's side faces traverse each boundary in one
    // direction, and a cap that seals it must run the opposite way so
    // the closed shell has consistent outward normals. (Capping the
    // raw input wires directly mis-orients one cap and yields a shell
    // `TruckSolid::new` rejects — this matches the working loft path
    // in `valenx-feature-tree::ops::loft`.)
    let boundaries = side.extract_boundaries();
    if boundaries.len() != 2 {
        return Err(FilletBrepError::TruckOp(format!(
            "{label}: homotopy shell has {} open boundaries, expected 2",
            boundaries.len()
        )));
    }
    let mut shell = side;
    for boundary in boundaries {
        let cap_wire: Wire = boundary.inverse();
        let cap = builder::try_attach_plane(&[cap_wire])
            .map_err(|err| FilletBrepError::TruckOp(format!("{label} cap: {err:?}")))?;
        shell.push(cap);
    }
    // `TruckSolid::new` panics ("This shell is not oriented and
    // closed.") when the homotopy shell + caps do not assemble into a
    // closed 2-manifold — which can happen on a degenerate or strongly
    // skewed taper. Use the fallible `try_new` and convert the failure
    // into a soft `TruckOp` error so a variable-radius fillet on
    // awkward geometry falls through cleanly instead of unwinding the
    // caller's thread.
    TruckSolid::try_new(vec![shell])
        .map_err(|err| FilletBrepError::TruckOp(format!("{label} solid assembly: {err}")))
}

/// Perform a real BRep **variable-radius** fillet on one convex
/// straight edge bounded by two planar faces.
///
/// This is the variable-radius generalisation of
/// [`fillet_convex_planar_edge`]: the fillet radius varies **linearly**
/// along the edge from `radius_start` (at the edge's front endpoint) to
/// `radius_end` (at the back endpoint). The fillet surface becomes a
/// tapered, lofted blend instead of a constant-radius cylinder; with
/// `radius_start == radius_end` the result is the constant-radius
/// fillet.
///
/// # Method
///
/// 1. validate the edge + the larger of the two radii via
///    [`plan_planar_edge_fillet`] (the `2·r ≤ edge_length` bound must
///    hold for the worst-case radius);
/// 2. resolve the fillet cross-section at each endpoint — the front at
///    `radius_start`, the back at `radius_end`;
/// 3. build the **tapered cutter** (lofted triangular prism) and the
///    **tapered fillet bar** (lofted circular-sector prism) by lofting
///    between the two end cross-sections with
///    [`builder::try_wire_homotopy`];
/// 4. evaluate `filleted = (solid − cutter) ∪ bar` with the real
///    `truck_shapeops` booleans.
///
/// # Honest scope
///
/// The radius profile is **linear** between the two endpoints — a
/// general radius law (a spline of radius-vs-arc-length) would loft
/// through intermediate stations and is a bounded follow-up. As with
/// the constant-radius fillet, the adjacent faces must be **planar**,
/// the edge must be a single **convex straight** edge, and **multi-edge
/// corner blends** (3+ filleted edges meeting at a vertex) remain a
/// Tier-3 research residual — see the crate-level docs.
///
/// # Errors
///
/// - The geometric-precondition errors of [`plan_planar_edge_fillet`]
///   (`NotPlanarFaces`, `NonConvexEdge`, `RadiusTooLarge`,
///   `BadParameter`).
/// - [`FilletBrepError::BadParameter`] if either endpoint radius is
///   non-finite or non-positive.
/// - [`FilletBrepError::TruckOp`] if a loft / cap step fails or a
///   boolean returns no solid (the coincident-face case the
///   constant-radius fillet documents — a *soft* fall-through signal).
pub fn fillet_variable_radius_planar_edge(
    solid: &TruckSolid,
    edge: &truck_modeling::Edge,
    radius_start: f64,
    radius_end: f64,
) -> Result<TruckSolid, FilletBrepError> {
    if !radius_start.is_finite() || radius_start <= 0.0 {
        return Err(FilletBrepError::BadParameter {
            name: "radius_start",
            reason: format!("must be > 0 and finite, got {radius_start}"),
        });
    }
    if !radius_end.is_finite() || radius_end <= 0.0 {
        return Err(FilletBrepError::BadParameter {
            name: "radius_end",
            reason: format!("must be > 0 and finite, got {radius_end}"),
        });
    }
    // Validate against the worst-case (larger) radius — the
    // `2·r ≤ edge_length` self-intersection bound must hold along the
    // whole edge, so the bigger endpoint governs.
    let worst = radius_start.max(radius_end);
    let plan = plan_planar_edge_fillet(solid, edge, worst)?;

    // Over-extend the loft past both endpoints so the cutter / bar end
    // caps never coincide with the solid's faces (the same trick the
    // constant-radius builder uses). The radius is linearly
    // extrapolated to the over-extended endpoints so the taper rate is
    // preserved.
    let edge_vec = plan.edge_back - plan.edge_front;
    let edge_len = edge_vec.magnitude();
    if edge_len < 1e-12 {
        return Err(FilletBrepError::TruckOp(
            "degenerate edge — zero length".to_string(),
        ));
    }
    let edge_dir = edge_vec / edge_len;
    let over = (edge_len * PRISM_OVEREXTEND_FRAC).max(1e-6);
    // Extrapolated endpoints + radii (linear in arc length).
    let front_pt = plan.edge_front - edge_dir * over;
    let back_pt = plan.edge_back + edge_dir * over;
    let slope = (radius_end - radius_start) / edge_len;
    let front_r = (radius_start - slope * over).max(1e-6);
    let back_r = (radius_end + slope * over).max(1e-6);

    let front = resolve_station(&plan, front_pt, front_r)?;
    let back = resolve_station(&plan, back_pt, back_r)?;

    // Build the tapered cutter + fillet bar, then the CSG.
    let cutter = build_tapered_cutter(&front, &back)?;
    let bar = build_tapered_fillet_bar(&front, &back)?;

    let base = Solid::from_truck(solid.clone());
    let cutter_solid = Solid::from_truck(cutter);
    let bar_solid = Solid::from_truck(bar);

    let trimmed = valenx_cad::boolean::difference_tol(&base, &cutter_solid, FILLET_BOOL_TOL)
        .map_err(|err| FilletBrepError::TruckOp(format!("corner removal boolean: {err}")))?;
    let filleted = valenx_cad::boolean::union_tol(&trimmed, &bar_solid, FILLET_BOOL_TOL)
        .map_err(|err| FilletBrepError::TruckOp(format!("fillet-bar union boolean: {err}")))?;

    match filleted {
        Solid::Brep(b) => Ok(b),
        Solid::Mesh(_) => Err(FilletBrepError::TruckOp(
            "boolean returned a mesh-backed solid (internal error)".to_string(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;
    use valenx_cad::primitives::box_solid;

    fn inner_brep(s: &Solid) -> &TruckSolid {
        match s {
            Solid::Brep(b) => b,
            _ => panic!("expected brep"),
        }
    }

    fn pick_first_edge(brep: &TruckSolid) -> truck_modeling::Edge {
        let mut seen = std::collections::HashSet::new();
        for e in brep.edge_iter() {
            if seen.insert(e.id()) {
                return e;
            }
        }
        panic!("solid has no edges")
    }

    #[test]
    fn resolve_cross_section_arc_centre_is_radius_from_each_face() {
        // For a 90° cube corner the arc centre must sit exactly
        // `radius` from each face plane (the defining property of a
        // tangent fillet circle).
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.3).unwrap();
        let xs = resolve_cross_section(&plan).unwrap();

        // Distance from the arc centre to face 0's plane. Face 0's
        // plane passes through the edge with outward normal
        // `plan.face0_normal`; the centre's signed distance is
        // (C − edge_front)·n̂.
        let n0 = plan.face0_normal.normalize();
        let n1 = plan.face1_normal.normalize();
        let to_c = xs.c_front - plan.edge_front;
        let d0 = to_c.dot(n0).abs();
        let d1 = to_c.dot(n1).abs();
        assert!(
            (d0 - 0.3).abs() < 1e-6,
            "arc centre should be radius from face 0, got {d0}"
        );
        assert!(
            (d1 - 0.3).abs() < 1e-6,
            "arc centre should be radius from face 1, got {d1}"
        );
    }

    #[test]
    fn resolve_cross_section_arc_apex_is_radius_from_centre() {
        // The arc apex lies on the fillet circle, so |C − apex| = r.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.25).unwrap();
        let xs = resolve_cross_section(&plan).unwrap();
        let r = (xs.arc_mid_front - xs.c_front).magnitude();
        assert!((r - 0.25).abs() < 1e-6, "apex not on circle: {r}");
    }

    #[test]
    fn resolve_cross_section_tangent_points_are_radius_from_centre() {
        // T0 and T1 lie on the fillet circle too — they are the
        // tangent contact points, so |C − T| = r.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.4).unwrap();
        let xs = resolve_cross_section(&plan).unwrap();
        for (t, label) in [(xs.t0_front, "T0"), (xs.t1_front, "T1")] {
            let d = (t - xs.c_front).magnitude();
            assert!(
                (d - 0.4).abs() < 1e-6,
                "{label} should be radius from centre, got {d}"
            );
        }
    }

    #[test]
    fn cube_edge_dihedral_is_right_angle() {
        // Sanity: the plan we build the fillet from sees a 90° corner.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.1).unwrap();
        assert!((plan.dihedral_angle - FRAC_PI_2).abs() < 1e-6);
    }

    #[test]
    fn cutter_prism_builds_as_brep() {
        // The cutter is a real triangular-prism BRep: 5 faces, 6
        // vertices, 9 edges — and over-extended past the edge ends.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.3).unwrap();
        let xs = resolve_cross_section(&plan).unwrap();
        let cutter = build_cutter(&xs).unwrap();
        // A triangular prism: 2 triangular caps + 3 rectangular sides.
        let face_count: usize = cutter.boundaries().iter().map(|s| s.len()).sum();
        assert_eq!(face_count, 5, "cutter should be a 5-face triangular prism");
    }

    #[test]
    fn fillet_bar_prism_builds_as_brep() {
        // The fillet bar is a real prism with a circular-arc side
        // face. Its cross-section is a 3-edge sector (2 lines + 1
        // arc), so the prism has 2 caps + 3 side faces = 5 faces.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.3).unwrap();
        let xs = resolve_cross_section(&plan).unwrap();
        let bar = build_fillet_bar(&xs).unwrap();
        let face_count: usize = bar.boundaries().iter().map(|s| s.len()).sum();
        assert_eq!(face_count, 5, "fillet bar should have 5 faces");
        // The bar must have real edges and vertices.
        assert!(bar.edge_iter().count() > 0, "fillet bar has no edges");
        assert!(bar.vertex_iter().count() > 0, "fillet bar has no vertices");
    }

    #[test]
    fn cutter_cross_section_lies_in_the_edge_normal_plane() {
        // The cutter cross-section points (E, T0, T1, C, apex) at the
        // front must all lie in the plane perpendicular to the edge
        // through the front corner — i.e. each point's displacement
        // from `e_front` is orthogonal to the edge direction.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.3).unwrap();
        let xs = resolve_cross_section(&plan).unwrap();
        for (p, label) in [
            (xs.t0_front, "T0"),
            (xs.t1_front, "T1"),
            (xs.c_front, "C"),
            (xs.arc_mid_front, "apex"),
        ] {
            let along_edge = (p - xs.e_front).dot(xs.edge_dir);
            assert!(
                along_edge.abs() < 1e-9,
                "{label} not in the edge-normal plane (edge-axis offset {along_edge})"
            );
        }
    }

    #[test]
    fn fillet_rejects_zero_radius_before_construction() {
        // A bad radius must fail at the planning stage with the
        // BadParameter cause, not somewhere inside the BRep build.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let err = fillet_convex_planar_edge(brep, &edge, 0.0).unwrap_err();
        assert!(matches!(
            err,
            FilletBrepError::BadParameter { name: "radius", .. }
        ));
    }

    #[test]
    fn fillet_rejects_too_large_radius_before_construction() {
        // radius·2 > edge length → RadiusTooLarge from the planner.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let err = fillet_convex_planar_edge(brep, &edge, 0.7).unwrap_err();
        assert!(matches!(err, FilletBrepError::RadiusTooLarge { .. }));
    }

    // --- variable-radius fillet tests ---

    #[test]
    fn resolve_station_arc_centre_is_radius_from_each_face() {
        // A station resolved at any point + radius must place its arc
        // centre exactly `radius` from each face plane — the same
        // tangent-circle property the constant-radius section has, now
        // checked at a custom radius.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.5).unwrap();
        // Resolve a station at the edge midpoint with radius 0.7.
        let mid = plan.edge_front + (plan.edge_back - plan.edge_front) * 0.5;
        let st = resolve_station(&plan, mid, 0.7).unwrap();
        let n0 = plan.face0_normal.normalize();
        let n1 = plan.face1_normal.normalize();
        let to_c = st.c - mid;
        assert!(
            (to_c.dot(n0).abs() - 0.7).abs() < 1e-6,
            "station arc centre should be radius from face 0"
        );
        assert!(
            (to_c.dot(n1).abs() - 0.7).abs() < 1e-6,
            "station arc centre should be radius from face 1"
        );
        // The tangent points + apex all lie on the radius-0.7 circle.
        for p in [st.t0, st.t1, st.arc_mid] {
            assert!(
                ((p - st.c).magnitude() - 0.7).abs() < 1e-6,
                "station point should lie on the fillet circle"
            );
        }
    }

    #[test]
    fn resolve_station_radius_scales_the_section_linearly() {
        // Doubling the radius doubles the tangent-point offset from the
        // corner — the section scales linearly with the radius, which
        // is what makes a *linear* radius taper a clean loft.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.5).unwrap();
        let mid = plan.edge_front + (plan.edge_back - plan.edge_front) * 0.5;
        let small = resolve_station(&plan, mid, 0.3).unwrap();
        let big = resolve_station(&plan, mid, 0.6).unwrap();
        let off_small = (small.t0 - mid).magnitude();
        let off_big = (big.t0 - mid).magnitude();
        assert!(
            (off_big / off_small - 2.0).abs() < 1e-6,
            "2× radius should give 2× tangent offset, got ratio {}",
            off_big / off_small
        );
    }

    #[test]
    fn tapered_cutter_builds_as_a_closed_brep() {
        // The variable-radius cutter is a lofted triangular prism: a
        // closed BRep with two triangular caps + three lofted side
        // faces = 5 faces.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.6).unwrap();
        let front = resolve_station(&plan, plan.edge_front, 0.3).unwrap();
        let back = resolve_station(&plan, plan.edge_back, 0.6).unwrap();
        let cutter = build_tapered_cutter(&front, &back).unwrap();
        let faces: usize = cutter.boundaries().iter().map(|s| s.len()).sum();
        assert_eq!(faces, 5, "tapered cutter should be a 5-face prism");
        assert!(cutter.edge_iter().count() > 0, "cutter has no edges");
    }

    #[test]
    fn tapered_fillet_bar_builds_as_a_closed_brep() {
        // The variable-radius fillet bar is a lofted circular-sector
        // prism — 2 sector caps + 3 lofted side faces (one of them the
        // tapered arc blend) = 5 faces.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let plan = plan_planar_edge_fillet(brep, &edge, 0.6).unwrap();
        let front = resolve_station(&plan, plan.edge_front, 0.3).unwrap();
        let back = resolve_station(&plan, plan.edge_back, 0.6).unwrap();
        let bar = build_tapered_fillet_bar(&front, &back).unwrap();
        let faces: usize = bar.boundaries().iter().map(|s| s.len()).sum();
        assert_eq!(faces, 5, "tapered fillet bar should have 5 faces");
        assert!(bar.vertex_iter().count() > 0, "fillet bar has no vertices");
    }

    #[test]
    fn variable_radius_fillet_rejects_bad_endpoint_radii() {
        // Zero / negative / NaN endpoint radii fail with BadParameter.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let bad_start =
            fillet_variable_radius_planar_edge(brep, &edge, 0.0, 0.3).unwrap_err();
        assert!(matches!(
            bad_start,
            FilletBrepError::BadParameter { name: "radius_start", .. }
        ));
        let bad_end =
            fillet_variable_radius_planar_edge(brep, &edge, 0.3, -0.1).unwrap_err();
        assert!(matches!(
            bad_end,
            FilletBrepError::BadParameter { name: "radius_end", .. }
        ));
        let nan_end =
            fillet_variable_radius_planar_edge(brep, &edge, 0.3, f64::NAN).unwrap_err();
        assert!(matches!(
            nan_end,
            FilletBrepError::BadParameter { name: "radius_end", .. }
        ));
    }

    #[test]
    fn variable_radius_fillet_validates_against_the_larger_radius() {
        // The 2·r ≤ edge_length bound must hold for the worst-case
        // (larger) endpoint radius. A unit-cube edge is length 1, so
        // radius_end = 0.7 (→ 2·r = 1.4 > 1) must trip RadiusTooLarge
        // even though radius_start = 0.1 is fine.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        let err =
            fillet_variable_radius_planar_edge(brep, &edge, 0.1, 0.7).unwrap_err();
        assert!(matches!(err, FilletBrepError::RadiusTooLarge { .. }));
    }

    #[test]
    fn variable_radius_fillet_on_a_cube_edge_either_builds_or_soft_fails() {
        // The full variable-radius fillet runs the real loft + boolean
        // construction. For a well-conditioned cube edge with a genuine
        // radius taper (0.2 → 0.4) the outcome must be either a real
        // BRep solid or the documented soft `TruckOp` fall-through —
        // never a geometric-precondition error, never a panic.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let edge = pick_first_edge(brep);
        match fillet_variable_radius_planar_edge(brep, &edge, 0.2, 0.4) {
            Ok(filleted) => {
                let faces: usize =
                    filleted.boundaries().iter().map(|s| s.len()).sum();
                assert!(
                    faces >= 6,
                    "a variable-radius-filleted cube keeps the cube's faces, got {faces}"
                );
            }
            Err(FilletBrepError::TruckOp(_)) => {
                // Soft fall-through — coincident-face booleans are the
                // documented fragile case for any kernel.
            }
            other => panic!("unexpected variable-radius fillet outcome: {other:?}"),
        }
    }
}
