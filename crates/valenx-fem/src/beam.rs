//! Native **3D structural beam / frame** finite-element solver
//! (Phase 24.8).
//!
//! ## What this is
//!
//! A genuine, self-contained finite-element solver for **2-node 3D beam
//! elements** — the element type frames, trusses, lattices and
//! structural skeletons are built from. Every node carries **six
//! degrees of freedom**: three translations `(u,v,w)` and three
//! rotations `(θx,θy,θz)`. One beam element therefore has a 12×12
//! stiffness matrix coupling
//!
//! - **axial** stretch (`EA/L`),
//! - **torsion** about the beam axis (`GJ/L`),
//! - **bending** in each of the two principal planes (`EI`), with
//!   **Timoshenko transverse-shear flexibility** (`GA·κ`) — so the
//!   element is correct for stocky beams as well as slender ones and
//!   does not shear-lock.
//!
//! The element is defined in a **local** frame (`x` along the beam,
//! `y`/`z` the cross-section principal axes); a 12×12 block-diagonal
//! rotation `T` built from the local triad maps it to global
//! coordinates. The solver assembles the global `K`, applies six-DOF
//! nodal constraints and loads (forces *and* moments), and solves the
//! sparse SPD system. A companion mass matrix gives a beam modal solve.
//!
//! ## Honest scope
//!
//! Real, validated v1 — it reproduces analytic cantilever / simply-
//! supported deflections and the first natural frequency of a beam (see
//! the tests and [`crate::validation`]). It is **prismatic, linear-
//! elastic, small-displacement**: one constant cross-section per
//! element, no tapering, no member loads (apply distributed load as
//! equivalent nodal forces), no geometric stiffening (`P-Δ`), no
//! offsets / rigid links / released DOFs. Those are bounded follow-ups;
//! the 90% frame-analysis case — axial + biaxial bending + torsion of
//! prismatic members — is covered.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};
use thiserror::Error;

use crate::material::FemMaterial;

/// The analytic **Euler–Bernoulli beam natural frequency**
/// `f = (β·L)²/(2π·L²)·√(E·I / (ρ·A))` (Hz) of a slender prismatic beam in
/// transverse (bending) vibration. `beta_l` is the dimensionless mode eigenvalue
/// `β·L` for the boundary condition — `π` for the simply-supported (pinned–pinned)
/// fundamental, `1.875104` for the cantilever (clamped–free) fundamental,
/// `4.730041` for clamped–clamped — `e_modulus` is `E` (Pa),
/// `area_moment_of_inertia` the section's `I` (m⁴) about the bending axis,
/// `density` `ρ` (kg/m³), `area` the cross-section `A` (m²), and `length` the
/// span `L` (m).
///
/// This is the *analytic* slender-beam reference the finite-element modal solver
/// ([`solve_beam_modal`], and the tet [`crate::modal_solver`]) converges to as
/// the mesh refines — the bending-vibration companion to the
/// [`crate::buckling::euler_critical_load`] stability reference. It gives a quick
/// hand-check without meshing: the frequency rises with the bending stiffness
/// `√(E·I)`, falls with the running mass `√(ρ·A)`, drops as `1/L²`, and scales
/// with the square of the mode eigenvalue `(β·L)²` (so the cantilever fundamental
/// sits a factor `(1.875104/π)² ≈ 0.356` below the simply-supported one). Returns
/// `0` for any non-physical input (`E`, `I`, `ρ`, `A`, or `L` non-positive, or any
/// argument non-finite).
pub fn euler_bernoulli_beam_frequency(
    beta_l: f64,
    e_modulus: f64,
    area_moment_of_inertia: f64,
    density: f64,
    area: f64,
    length: f64,
) -> f64 {
    if !beta_l.is_finite()
        || !e_modulus.is_finite()
        || e_modulus <= 0.0
        || !area_moment_of_inertia.is_finite()
        || area_moment_of_inertia <= 0.0
        || !density.is_finite()
        || density <= 0.0
        || !area.is_finite()
        || area <= 0.0
        || !length.is_finite()
        || length <= 0.0
    {
        return 0.0;
    }
    let beta = beta_l / length; // β = (β·L)/L, the bending wavenumber (1/m)
    beta * beta / (2.0 * std::f64::consts::PI)
        * (e_modulus * area_moment_of_inertia / (density * area)).sqrt()
}

/// The analytic **cantilever tip deflection** `δ = P·L³/(3·E·I)` (m) of a
/// slender Euler–Bernoulli cantilever — a prismatic beam clamped at one end and
/// loaded by a transverse point force `load` `P` (N) at the free end, with span
/// `length` `L` (m), Young's modulus `youngs_modulus` `E` (Pa), and section
/// second moment of area `second_moment_area` `I` (m⁴) about the bending axis.
///
/// This is the *analytic* slender-beam reference the finite-element beam solver
/// ([`solve_beam_static`]) converges to as the mesh refines — the static-bending
/// companion to the [`euler_bernoulli_beam_frequency`] vibration reference and
/// the [`crate::buckling::euler_critical_load`] stability reference. It gives a
/// quick hand-check without meshing: the deflection grows *linearly* with the
/// load `P` (and is sign-preserving — an upward load lifts the tip), with the
/// *cube* of the span `L` (the dominant lever: doubling the length softens the
/// tip eight-fold), and falls inversely with the flexural rigidity `E·I`. (A real
/// short/Timoshenko beam adds a small shear term `P·L/(κ·G·A)` on top; this is
/// the pure-bending part.) Returns `0` for non-physical input (`P` non-finite, or
/// `E`, `I`, or `L` non-positive or non-finite).
pub fn cantilever_tip_deflection(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * length.powi(3) / (3.0 * youngs_modulus * second_moment_area)
}

/// The analytic **free-end deflection of a cantilever under a pure end moment**
/// `δ = M₀·L² / (2·E·I)` (m) — the tip deflection of a cantilever of span `length` `L`
/// (m) loaded only by a couple `end_moment` `M₀` (N·m) applied at its free end, for
/// Young's modulus `youngs_modulus` `E` (Pa) and section second moment of area
/// `second_moment_area` `I` (m⁴). A pure end couple is carried as a *constant* bending
/// moment `M₀` along the whole span, so the curvature `κ = M₀/(EI)` ([`beam_curvature`])
/// is uniform and the beam bows into a circular arc; the small-deflection tip rise is
/// `κ·L²/2`. It is the moment-loaded companion to the force-loaded
/// [`cantilever_tip_deflection`] (`P·L³/3EI`): the deflection grows *linearly* with the
/// moment (and follows its sign), with the *square* of the span, and falls inversely with
/// the flexural rigidity `E·I`. Returns `0` for non-physical input (`M₀` non-finite, or
/// `E`, `I`, or `L` non-positive or non-finite).
pub fn cantilever_end_moment_tip_deflection(
    end_moment: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !end_moment.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    end_moment * length * length / (2.0 * youngs_modulus * second_moment_area)
}

/// The analytic **free-end slope of a cantilever under a pure end moment**
/// `θ = M₀·L / (E·I)` (rad) — the rotation of the free end of a cantilever of span
/// `length` `L` (m) loaded only by a couple `end_moment` `M₀` (N·m) at its tip, for
/// Young's modulus `youngs_modulus` `E` (Pa) and section second moment of area
/// `second_moment_area` `I` (m⁴). A pure end couple is carried as a *constant* bending
/// moment `M₀`, so the curvature `κ = M₀/(EI)` ([`beam_curvature`]) is uniform and the
/// slope accumulates linearly along the span, `θ = κ·L`. It is the rotational companion to
/// the [`cantilever_end_moment_tip_deflection`] `δ = M₀L²/2EI`, the two locked together by
/// the constant-curvature arc `δ = θ·L/2`. The slope grows *linearly* with the moment (and
/// follows its sign), *linearly* with the span, and falls inversely with the flexural
/// rigidity `E·I`. Returns `0` for non-physical input (`M₀` non-finite, or `E`, `I`, or
/// `L` non-positive or non-finite).
pub fn cantilever_end_moment_tip_slope(
    end_moment: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !end_moment.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    end_moment * length / (youngs_modulus * second_moment_area)
}

/// The analytic **strain energy of a cantilever under a pure end moment**
/// `U = M₀²·L / (2·E·I)` (J) — the elastic bending energy stored in a cantilever of span
/// `length` `L` (m) loaded only by a couple `end_moment` `M₀` (N·m) at its free end, for
/// Young's modulus `youngs_modulus` `E` (Pa) and section second moment of area
/// `second_moment_area` `I` (m⁴). A pure end couple is carried as a *constant* bending
/// moment `M₀`, so the energy integral `U = ∫M²/(2EI) dx` collapses to `M₀²L/(2EI)`. It is
/// the moment-loaded member of the cantilever energy-method family alongside the point-load
/// [`cantilever_point_load_strain_energy`] and UDL [`cantilever_udl_strain_energy`] cases,
/// and follows from **Clapeyron's theorem** `U = ½·M₀·θ` with the tip slope
/// [`cantilever_end_moment_tip_slope`] `θ = M₀L/EI`. The energy grows with the *square* of
/// the moment (so it is sign-independent), *linearly* with the span, and falls inversely
/// with the flexural rigidity `E·I`. Returns `0` for non-physical input (`M₀` non-finite,
/// or `E`, `I`, or `L` non-positive or non-finite).
pub fn cantilever_end_moment_strain_energy(
    end_moment: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !end_moment.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    end_moment * end_moment * length / (2.0 * youngs_modulus * second_moment_area)
}

/// The analytic **maximum bending moment at the fixed root of a tip-loaded
/// cantilever** `M = P·L` (N·m) — the peak moment, at the built-in (encastré) end,
/// which sets the maximum bending stress and so governs the strength design of the
/// member. `load` `P` is the transverse tip force (N) and `length` `L` the span (m).
///
/// The moment varies linearly from `P·L` at the root to zero at the free tip; this
/// root value is the design maximum. It threads the tip deflection
/// [`cantilever_tip_deflection`] (`δ = P·L³/3EI = M·L²/3EI`) and the strain energy
/// [`cantilever_point_load_strain_energy`] (`U = P²L³/6EI = M²·L/6EI`). Linear and
/// sign-preserving in `P`, linear in `L`, and — a statics result — independent of `E`
/// and `I`. Returns `0` for non-physical input (`P` non-finite, or `L` non-positive
/// or non-finite).
pub fn cantilever_point_load_root_moment(load: f64, length: f64) -> f64 {
    if !load.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load * length
}

/// The analytic **maximum bending moment at the fixed root of a cantilever under a
/// uniformly distributed load** `M = w·L²/2` (N·m) — the peak moment, at the built-in
/// (encastré) end, which sets the maximum bending stress and governs the strength
/// design of the member. `load_per_length` `w` is the load intensity (N/m) and
/// `length` `L` the span (m).
///
/// The moment grows quadratically from zero at the free tip to `w·L²/2` at the root;
/// this root value is the design maximum. It is the UDL companion to the point-load
/// [`cantilever_point_load_root_moment`] (`P·L`), and threads the tip deflection
/// [`cantilever_udl_tip_deflection`] (`δ = w·L⁴/8EI = M·L²/4EI`) and the strain energy
/// [`cantilever_udl_strain_energy`] (`U = w²L⁵/40EI = M²·L/10EI`). Quadratic in `L`,
/// linear and sign-preserving in `w`, and — a statics result — independent of `E` and
/// `I`. Returns `0` for non-physical input (`w` non-finite, or `L` non-positive or
/// non-finite).
pub fn cantilever_udl_root_moment(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load_per_length * length * length / 2.0
}

/// The analytic **maximum transverse shear force at the fixed root of a cantilever under a
/// uniformly distributed load** `V = w·L` (N) — the peak shear, carried at the built-in
/// (encastré) end, equal to the total distributed load on the span; it sets the maximum
/// transverse shear stress and governs the shear design of the member. `load_per_length`
/// `w` is the load intensity (N/m) and `length` `L` the span (m).
///
/// The shear decreases linearly from `w·L` at the root to zero at the free tip. It is the
/// shear companion to the root bending moment [`cantilever_udl_root_moment`]
/// (`M = w·L²/2 = V·L/2`). Linear and sign-preserving in `w`, linear in `L`, and — a statics
/// result — independent of `E` and `I`. Returns `0` for non-physical input (`w` non-finite,
/// or `L` non-positive or non-finite).
pub fn cantilever_udl_max_shear(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load_per_length * length
}

/// The analytic **strain energy of a tip-loaded cantilever**
/// `U = P²·L³/(6·E·I)` (J) — the elastic energy stored in bending when a slender
/// Euler–Bernoulli cantilever of span `length` `L` (m), Young's modulus
/// `youngs_modulus` `E` (Pa) and section second moment of area `second_moment_area`
/// `I` (m⁴) carries a transverse point force `load` `P` (N) at its free end.
///
/// It is the bending-energy integral `U = ∫₀^L M²/(2EI) dx` with the linearly
/// varying moment `M(x) = P·x`, and the first of the energy-method references
/// (Castigliano's theorems, unit-load deflections). By **Clapeyron's theorem** the
/// work a single static load does equals half the load times the deflection it
/// produces, so `U = ½·P·δ_tip` with the tip deflection
/// [`cantilever_tip_deflection`] `δ = P·L³/(3EI)` — a useful consistency check and
/// the basis of energy methods for deflection. The energy grows with the *square*
/// of the load (so it is sign-independent — an up- or down-load store the same
/// energy), with the *cube* of the span, and falls inversely with the flexural
/// rigidity `E·I`. Returns `0` for non-physical input (`P` non-finite, or `E`, `I`,
/// or `L` non-positive or non-finite).
pub fn cantilever_point_load_strain_energy(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * load * length.powi(3) / (6.0 * youngs_modulus * second_moment_area)
}

/// The analytic **cantilever tip slope** `θ = P·L²/(2·E·I)` (rad) — the
/// end-rotation of a slender Euler–Bernoulli cantilever clamped at one end and
/// loaded by a transverse point force `load` `P` (N) at the free end, with span
/// `length` `L` (m), Young's modulus `youngs_modulus` `E` (Pa), and section second
/// moment of area `second_moment_area` `I` (m⁴).
///
/// This is the *rotational* companion to the deflection
/// [`cantilever_tip_deflection`] — the slope of the elastic curve at the free end,
/// the quantity slope-continuity and moment-area methods track. The two are locked
/// together by `δ = (2/3)·L·θ`: integrating the curvature `M/(E·I)` once gives the
/// slope (`∝ L²`), twice the deflection (`∝ L³`), so the tip deflection is
/// two-thirds of the span times the tip slope. Like the deflection it grows
/// linearly with the load `P` (and is sign-preserving — an upward load rotates the
/// tip up), with the *square* of the span `L`, and falls inversely with the
/// flexural rigidity `E·I`. Returns `0` for non-physical input (`P` non-finite, or
/// `E`, `I`, or `L` non-positive or non-finite).
pub fn cantilever_tip_slope(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * length.powi(2) / (2.0 * youngs_modulus * second_moment_area)
}

/// The analytic **cantilever tip deflection under a uniformly distributed load**
/// `δ = w·L⁴/(8·E·I)` (m) of a slender Euler–Bernoulli cantilever — a prismatic
/// beam clamped at one end carrying a transverse load of intensity
/// `load_per_length` `w` (N/m) spread evenly over its full span `length` `L` (m),
/// with Young's modulus `youngs_modulus` `E` (Pa) and section second moment of
/// area `second_moment_area` `I` (m⁴) about the bending axis.
///
/// This is the **distributed-load companion** to the point-load
/// [`cantilever_tip_deflection`] (`P·L³/(3·E·I)`): the *same total load*
/// `W = w·L`, spread uniformly along the span instead of concentrated at the
/// free end, deflects the tip only `3/8` as far — the load near the clamp acts on
/// a short lever arm and contributes little. The deflection grows *linearly* with
/// the intensity `w` (and is sign-preserving — an upward load lifts the tip) and
/// with the *fourth* power of the span `L` (doubling the length softens the tip
/// sixteen-fold), and falls inversely with the flexural rigidity `E·I`. It is the
/// self-weight / pressure-loading member of the analytic beam-reference set the
/// finite-element solver ([`solve_beam_static`]) converges to. Returns `0` for
/// non-physical input (`w` non-finite, or `E`, `I`, or `L` non-positive or
/// non-finite).
pub fn cantilever_udl_tip_deflection(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load_per_length * length.powi(4) / (8.0 * youngs_modulus * second_moment_area)
}

/// The analytic **strain energy of a uniformly-loaded cantilever**
/// `U = w²·L⁵/(40·E·I)` (J) — the elastic energy stored in bending when a slender
/// Euler–Bernoulli cantilever of span `length` `L` (m), Young's modulus
/// `youngs_modulus` `E` (Pa) and section second moment of area `second_moment_area`
/// `I` (m⁴) carries a uniformly distributed transverse load of intensity
/// `load_per_length` `w` (N/m) over its full length.
///
/// It extends the energy-method family (Castigliano / Clapeyron) into the
/// distributed-load case — the UDL companion to the point-load
/// [`cantilever_point_load_strain_energy`] — and is the bending-energy integral
/// `U = ∫₀^L M²/(2EI) dx` with the cantilever-UDL moment `M(x) = w(L−x)²/2`. Because
/// the load is distributed, the tip deflection alone does not give the work
/// directly, but `U = (1/5)·(w·L)·δ_tip` ties it to
/// [`cantilever_udl_tip_deflection`] `δ = wL⁴/(8EI)`. The energy grows with the
/// *square* of the load intensity (so it is sign-independent), the *fifth* power of
/// the span, and falls inversely with the flexural rigidity `E·I`. Returns `0` for
/// non-physical input (`w` non-finite, or `E`, `I`, or `L` non-positive or
/// non-finite).
pub fn cantilever_udl_strain_energy(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load_per_length * load_per_length * length.powi(5)
        / (40.0 * youngs_modulus * second_moment_area)
}

/// The analytic **strain energy of a simply-supported beam under a uniformly
/// distributed load** `U = w²·L⁵/(240·E·I)` (J) — the elastic energy stored in
/// bending when a slender Euler–Bernoulli beam of span `length` `L` (m), Young's
/// modulus `youngs_modulus` `E` (Pa) and section second moment of area
/// `second_moment_area` `I` (m⁴) carries a uniformly distributed transverse load of
/// intensity `load_per_length` `w` (N/m) over its full span.
///
/// It is the **last corner of the energy-method 2×2 matrix** {cantilever,
/// simply-supported} × {point, UDL}, the simply-supported-UDL companion to the
/// cantilever [`cantilever_udl_strain_energy`]. It is the bending-energy integral
/// `U = ∫₀^L M²/(2EI) dx` with the simply-supported-UDL moment `M(x) = (w/2)·x·(L−x)`
/// (parabolic, zero at the pinned ends, peaking at mid-span). For the same `w, L, E,
/// I` a cantilever stores exactly `6×` this energy (`w²L⁵/(40EI)` vs `w²L⁵/(240EI)`)
/// — the far more compliant free end. The energy grows with the *square* of the load
/// intensity (so it is sign-independent), the *fifth* power of the span, and falls
/// inversely with the flexural rigidity `E·I`. Returns `0` for non-physical input
/// (`w` non-finite, or `E`, `I`, or `L` non-positive or non-finite).
pub fn simply_supported_udl_strain_energy(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load_per_length * load_per_length * length.powi(5)
        / (240.0 * youngs_modulus * second_moment_area)
}

/// The analytic **cantilever tip slope under a uniformly distributed load**
/// `θ = w·L³/(6·E·I)` (rad) — the end-rotation of a slender Euler–Bernoulli
/// cantilever carrying a transverse load of intensity `load_per_length` `w` (N/m)
/// over its full span `length` `L` (m), with Young's modulus `youngs_modulus` `E`
/// (Pa) and section second moment of area `second_moment_area` `I` (m⁴).
///
/// This is the *distributed-load* slope companion to the point-load
/// [`cantilever_tip_slope`] (`P·L²/(2EI)`) and the rotational companion to the UDL
/// deflection [`cantilever_udl_tip_deflection`] (`w·L⁴/(8EI)`); the two UDL
/// quantities are locked by `δ = (3/4)·L·θ`. The slope grows linearly with the
/// intensity `w` (sign-preserving), with the *cube* of the span `L`, and falls
/// inversely with the flexural rigidity `E·I`. Returns `0` for non-physical input
/// (`w` non-finite, or `E`, `I`, or `L` non-positive or non-finite).
pub fn cantilever_udl_tip_slope(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load_per_length * length.powi(3) / (6.0 * youngs_modulus * second_moment_area)
}

/// The analytic **angle of twist** `θ = T·L/(G·J)` (rad) of a prismatic shaft in
/// pure torsion — a shaft of length `length` `L` (m), shear modulus
/// `shear_modulus` `G` (Pa) and polar second moment of area (the St-Venant
/// torsion constant) `polar_moment` `J` (m⁴), twisted by an axial torque `torque`
/// `T` (N·m).
///
/// This is the *torsion* member of the analytic beam-reference set the
/// finite-element solver ([`solve_beam_static`]) converges to — the companion to
/// the bending [`cantilever_tip_deflection`], the
/// [`euler_bernoulli_beam_frequency`] vibration reference and the
/// [`crate::buckling::euler_critical_load`] stability reference, each covering a
/// different deformation mode. The twist grows *linearly* with the torque `T`
/// (and is sign-preserving — reversing the torque reverses the twist) and with
/// the length `L`, and falls inversely with the **torsional rigidity** `G·J` (a
/// stiffer material or a fatter section twists less). For a solid circular shaft
/// `J = πr⁴/2`. Returns `0` for non-physical input (`T` non-finite, or `G`, `J`,
/// or `L` non-positive or non-finite).
pub fn beam_angle_of_twist(
    torque: f64,
    length: f64,
    shear_modulus: f64,
    polar_moment: f64,
) -> f64 {
    if !torque.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !shear_modulus.is_finite()
        || shear_modulus <= 0.0
        || !polar_moment.is_finite()
        || polar_moment <= 0.0
    {
        return 0.0;
    }
    torque * length / (shear_modulus * polar_moment)
}

/// The **torsional shear stress** `τ = T·r/J` (Pa) — the shear stress a torque `torque`
/// `T` (N·m) induces at radius `radius` `r` (m) from the axis of a shaft of polar second
/// moment of area `polar_moment` `J` (m⁴). Like the bending [`bending_stress`] it varies
/// linearly across the section — zero on the axis, maximal at the outer surface
/// (`r = r_outer`), which sets the failure check — and it is the shear conjugate of the
/// [`beam_angle_of_twist`] through Hooke's law `τ = G·γ` with the shear strain
/// `γ = r·θ/L`. Returns `0` for non-physical input (`T` or `r` non-finite, or `J`
/// non-positive or non-finite).
pub fn torsional_shear_stress(torque: f64, radius: f64, polar_moment: f64) -> f64 {
    if !torque.is_finite() || !radius.is_finite() || !polar_moment.is_finite() || polar_moment <= 0.0
    {
        return 0.0;
    }
    torque * radius / polar_moment
}

/// The **torsional strain energy** `U = T²·L / (2·G·J)` (J) stored in a prismatic shaft
/// of length `length` `L` (m), shear modulus `shear_modulus` `G` (Pa) and polar second
/// moment of area `polar_moment` `J` (m⁴) carrying a torque `torque` `T` (N·m). It is the
/// elastic work done twisting the shaft — equivalently the Clapeyron form `½·T·θ` with the
/// [`beam_angle_of_twist`] `θ` — and the torsion analogue of the bending strain energy.
/// Being a square in `T` it is non-negative regardless of the torque's sign. Returns `0`
/// for non-physical input (non-finite torque, or a non-positive / non-finite length, shear
/// modulus or polar moment).
pub fn torsional_strain_energy(
    torque: f64,
    length: f64,
    shear_modulus: f64,
    polar_moment: f64,
) -> f64 {
    if !torque.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !shear_modulus.is_finite()
        || shear_modulus <= 0.0
        || !polar_moment.is_finite()
        || polar_moment <= 0.0
    {
        return 0.0;
    }
    torque * torque * length / (2.0 * shear_modulus * polar_moment)
}

/// The **polar (torsional) section modulus** `Z_p = J / r` (m³) — the torsion design
/// property that maps an applied torque to the peak surface shear stress through
/// `τ_max = T / Z_p`, where `polar_moment` `J` (m⁴) is the polar second moment of area
/// and `outer_radius` `r` (m) is the outer radius. It is the torsion analogue of the
/// bending [`elastic_section_modulus`] (`S = I/c`), and the conjugate of
/// [`torsional_shear_stress`]. Returns `0` for non-physical input (`J` or `r` non-positive
/// or non-finite).
pub fn polar_section_modulus(polar_moment: f64, outer_radius: f64) -> f64 {
    if !polar_moment.is_finite()
        || polar_moment <= 0.0
        || !outer_radius.is_finite()
        || outer_radius <= 0.0
    {
        return 0.0;
    }
    polar_moment / outer_radius
}

/// The **torsional moment (torque) capacity** `T = τ·Z_p` (N·m) — the torque a shaft
/// carries when its outer surface reaches the shear stress `shear_stress` `τ` (Pa),
/// given its [`polar_section_modulus`] `polar_modulus` `Z_p` (m³). It is the inverse of
/// the torsion formula `τ_max = T/Z_p` ([`torsional_shear_stress`] at the outer radius):
/// feeding the *allowable* shear stress gives the section's allowable torque, and feeding
/// the *yield* shear stress gives the torque at first yield. The torsion analogue of the
/// bending [`bending_moment_capacity`]. Linear in both the stress and the polar modulus.
/// Returns `0` for non-physical input (non-finite, or a non-positive polar modulus).
pub fn torsional_moment_capacity(shear_stress: f64, polar_modulus: f64) -> f64 {
    if !shear_stress.is_finite() || !polar_modulus.is_finite() || polar_modulus <= 0.0 {
        return 0.0;
    }
    shear_stress * polar_modulus
}

/// The analytic **axial extension** `δ = F·L/(E·A)` (m) of a prismatic bar in
/// pure tension or compression — a member of length `length` `L` (m), Young's
/// modulus `youngs_modulus` `E` (Pa) and cross-section area `area` `A` (m²) under
/// an axial force `force` `F` (N). This is Hooke's law for an axially-loaded bar
/// (`δ = F·L/(E·A)`, equivalently `σ = E·ε`).
///
/// It is the *axial* member of the analytic beam-reference set the finite-element
/// solver ([`solve_beam_static`]) converges to — the simplest deformation mode,
/// completing the trio with the bending [`cantilever_tip_deflection`] and the
/// torsional [`beam_angle_of_twist`] (the three ways a straight member yields,
/// under a transverse load, a torque, and an axial force). The extension grows
/// *linearly* with the force `F` (and is sign-preserving — tension lengthens,
/// compression shortens) and with the length `L`, and falls inversely with the
/// **axial rigidity** `E·A`. Returns `0` for non-physical input (`F` non-finite,
/// or `E`, `A`, or `L` non-positive or non-finite).
pub fn beam_axial_extension(force: f64, length: f64, youngs_modulus: f64, area: f64) -> f64 {
    if !force.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !area.is_finite()
        || area <= 0.0
    {
        return 0.0;
    }
    force * length / (youngs_modulus * area)
}

/// The **axial strain energy** `U = F²·L / (2·E·A)` (J) stored in a prismatic bar of
/// length `length` `L` (m), Young's modulus `youngs_modulus` `E` (Pa) and cross-section
/// area `area` `A` (m²) under an axial force `force` `F` (N). It is the elastic work done
/// stretching (or compressing) the bar — equivalently the Clapeyron form `½·F·δ` with the
/// [`beam_axial_extension`] `δ` — and the axial member of the strain-energy set
/// (axial / torsion / bending). Being a square in `F` it is non-negative: tension and
/// compression of equal magnitude store equal energy. Returns `0` for non-physical input
/// (non-finite force, or a non-positive / non-finite length, modulus or area).
pub fn axial_strain_energy(force: f64, length: f64, youngs_modulus: f64, area: f64) -> f64 {
    if !force.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !area.is_finite()
        || area <= 0.0
    {
        return 0.0;
    }
    force * force * length / (2.0 * youngs_modulus * area)
}

/// The **axial (extensional) rigidity** `EA = E·A` (N) — the product of Young's modulus
/// `youngs_modulus` `E` (Pa) and the cross-sectional area `area` `A` (m²). It is the
/// constant of proportionality between axial force and strain (`F = EA·ε`) and the
/// extensional companion to the [`flexural_rigidity`] (`EI`): the stiffness in the
/// [`beam_axial_extension`] (`δ = FL/EA`) and the [`axial_strain_energy`] (`U = F²L/2EA`).
/// A stiffer or thicker bar stretches less under the same load. Returns `0` for
/// non-physical input (`E` or `A` non-positive or non-finite).
pub fn axial_rigidity(youngs_modulus: f64, area: f64) -> f64 {
    if !youngs_modulus.is_finite() || youngs_modulus <= 0.0 || !area.is_finite() || area <= 0.0 {
        return 0.0;
    }
    youngs_modulus * area
}

/// The **torsional rigidity** `GJ = G·J` (N·m²) — the product of the shear modulus
/// `shear_modulus` `G` (Pa) and the polar moment of inertia (or torsion constant)
/// `polar_moment_of_inertia` `J` (m⁴). It is the constant of proportionality between torque
/// and twist-rate (`T = GJ·dφ/dz`) and the torsional member of the rigidity trio with the
/// [`axial_rigidity`] (`EA`) and the [`flexural_rigidity`] (`EI`): the stiffness in the
/// [`beam_angle_of_twist`] (`φ = TL/GJ`) and the [`torsional_strain_energy`] (`U = T²L/2GJ`).
/// A stiffer or fatter shaft twists less under the same torque. Returns `0` for non-physical
/// input (`G` or `J` non-positive or non-finite).
pub fn torsional_rigidity(shear_modulus: f64, polar_moment_of_inertia: f64) -> f64 {
    if !shear_modulus.is_finite()
        || shear_modulus <= 0.0
        || !polar_moment_of_inertia.is_finite()
        || polar_moment_of_inertia <= 0.0
    {
        return 0.0;
    }
    shear_modulus * polar_moment_of_inertia
}

/// The **axial normal stress** `σ = F / A` (Pa) in a prismatic bar of cross-section area
/// `area` `A` (m²) under an axial force `force` `F` (N) — positive in tension, negative in
/// compression. It is the axial member of the beam stress family alongside the bending
/// [`bending_stress`] (`σ = M·y/I`) and the [`torsional_shear_stress`] (`τ = T·r/J`), and
/// the Hooke's-law conjugate of the [`beam_axial_extension`] (`σ = E·δ/L`). Returns `0` for
/// non-physical input (`F` or `A` non-finite, or `A` non-positive).
pub fn axial_stress(force: f64, area: f64) -> f64 {
    if !force.is_finite() || !area.is_finite() || area <= 0.0 {
        return 0.0;
    }
    force / area
}

/// The **axial force capacity** `F = σ·A` (N) — the axial force a prismatic bar carries
/// when its cross-section of area `area` `A` (m²) reaches the normal stress `stress` `σ`
/// (Pa), positive in tension and negative in compression. It is the inverse of
/// [`axial_stress`] (`σ = F/A` ⟹ `F = σ·A`): feeding the *allowable* stress gives the
/// member's allowable axial load, and feeding the *yield* stress gives the squash / yield
/// load. The axial member of the capacity family alongside the bending
/// [`bending_moment_capacity`] and the [`torsional_moment_capacity`]. Linear in both the
/// stress and the area. Returns `0` for non-physical input (non-finite, or a non-positive
/// area).
pub fn axial_force_capacity(stress: f64, area: f64) -> f64 {
    if !stress.is_finite() || !area.is_finite() || area <= 0.0 {
        return 0.0;
    }
    stress * area
}

/// The **moment–curvature relation** `κ = M/(E·I)` (1/m) — the local bending curvature
/// produced in a beam section by a bending moment `moment` `M` (N·m), for Young's
/// modulus `youngs_modulus` `E` (Pa) and second moment of area `second_moment_area`
/// `I` (m⁴). This is the Euler–Bernoulli constitutive law `M = E·I·κ` that underlies
/// every deflection and slope in this module: integrating the curvature once gives the
/// slope and twice the deflection. The curvature is linear in the moment (and so
/// follows its sign — sagging vs hogging) and falls inversely with the flexural
/// rigidity `E·I`. Returns `0` for non-physical input (`M` non-finite, or `E` or `I`
/// non-positive or non-finite).
pub fn beam_curvature(moment: f64, youngs_modulus: f64, second_moment_area: f64) -> f64 {
    if !moment.is_finite()
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    moment / (youngs_modulus * second_moment_area)
}

/// The **flexural rigidity (bending stiffness)** `EI = E·I` (N·m²) — the product of
/// Young's modulus `youngs_modulus` `E` (Pa) and the section's second moment of area
/// `second_moment_area` `I` (m⁴). It is the constant of proportionality between bending
/// moment and curvature ([`beam_curvature`]: `M = EI·κ`) and the stiffness in every
/// Euler–Bernoulli deflection formula (e.g. [`cantilever_tip_deflection`] `δ = PL³/(3·EI)`):
/// a stiffer or deeper section resists bending more. Returns `0` for non-physical input
/// (`E` or `I` non-positive or non-finite).
pub fn flexural_rigidity(youngs_modulus: f64, second_moment_area: f64) -> f64 {
    if !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    youngs_modulus * second_moment_area
}

/// The **bulk modulus** `K = E / (3(1−2ν))` (Pa) — the volumetric (compressive)
/// stiffness of an isotropic elastic material: the ratio of a uniform hydrostatic
/// pressure to the fractional volume change `−p / (ΔV/V)` it produces. It is the third
/// of the standard isotropic elastic moduli alongside Young's modulus `youngs_modulus`
/// `E` (Pa) and the shear modulus `G` (the private companion the beam solver uses), and
/// the one that governs how a material resists a change of *volume* rather than of
/// shape. `poisson_ratio` `ν` is the dimensionless lateral-contraction ratio. The bulk
/// modulus rises without bound as `ν → 0.5` (an incompressible material) and reduces to
/// `E/3` at `ν = 0`; the relation inverts as `E = 3·K·(1−2ν)`. Returns `0` for
/// non-physical input (`E` non-positive or non-finite, or `ν` outside the physical
/// range `(−1, 0.5)`).
pub fn bulk_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    if !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !poisson_ratio.is_finite()
        || poisson_ratio <= -1.0
        || poisson_ratio >= 0.5
    {
        return 0.0;
    }
    youngs_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio))
}

/// The **Lamé first parameter** `λ = E·ν / ((1+ν)(1−2ν))` (Pa) — the first of the two
/// Lamé constants that parameterise the isotropic linear-elastic stiffness tensor,
/// `σ = λ·tr(ε)·I + 2μ·ε`, in which the second Lamé parameter `μ` is the shear modulus
/// `G` (the private companion the beam solver uses). `youngs_modulus` `E` (Pa) and
/// `poisson_ratio` `ν` are the engineering constants; `λ` is the coefficient that couples
/// the volumetric strain `tr(ε)` into the normal stresses. It connects to the
/// [`bulk_modulus`] `K` by `λ = K − 2G/3 = 3·K·ν/(1+ν)`. Returns `0` for non-physical
/// input (`E` non-positive or non-finite, or `ν` outside the physical range `(−1, 0.5)`);
/// note `λ` is also *validly* `0` at `ν = 0` and *negative* for an auxetic material
/// (`ν < 0`).
pub fn lames_first_parameter(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    if !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !poisson_ratio.is_finite()
        || poisson_ratio <= -1.0
        || poisson_ratio >= 0.5
    {
        return 0.0;
    }
    youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
}

/// The **shear modulus** (modulus of rigidity) `G = E / (2(1+ν))` (Pa) — the second of
/// the two Lamé constants (`μ`), relating a shear stress to the shear strain it produces.
/// `youngs_modulus` `E` (Pa) and `poisson_ratio` `ν` are the engineering constants. With
/// the [`bulk_modulus`] `K` and the [`lames_first_parameter`] `λ` it completes the
/// isotropic elastic-constant set: `λ = K − 2G/3`, and Young's modulus is recovered as
/// `E = 9KG / (3K + G)`. (A closed-form public companion to the crate's private
/// material-based shear-modulus helper.) Returns `0` for non-physical input (`E`
/// non-positive or non-finite, or `ν` outside the physical range `(−1, 0.5)`).
pub fn shear_modulus_from_youngs(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    if !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !poisson_ratio.is_finite()
        || poisson_ratio <= -1.0
        || poisson_ratio >= 0.5
    {
        return 0.0;
    }
    youngs_modulus / (2.0 * (1.0 + poisson_ratio))
}

/// The **P-wave (longitudinal / constrained) modulus** `M = E(1−ν) / ((1+ν)(1−2ν))` (Pa)
/// — the uniaxial-strain stiffness: the stress needed to compress a material along one
/// axis while it is *constrained* from straining laterally (an oedometer / 1-D
/// consolidation test). It also sets the compressional (primary) wave speed
/// `v_p = √(M/ρ)`. `youngs_modulus` `E` (Pa) and `poisson_ratio` `ν` are the engineering
/// constants. It is the capstone of the isotropic elastic-constant family, tying it all
/// together: `M = K + 4G/3 = λ + 2G`, where `K` is the [`bulk_modulus`], `G` the
/// [`shear_modulus_from_youngs`], and `λ` the [`lames_first_parameter`]. At `ν = 0` it
/// reduces to `E` (no lateral constraint). Returns `0` for non-physical input (`E`
/// non-positive or non-finite, or `ν` outside the physical range `(−1, 0.5)`).
pub fn p_wave_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    if !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !poisson_ratio.is_finite()
        || poisson_ratio <= -1.0
        || poisson_ratio >= 0.5
    {
        return 0.0;
    }
    youngs_modulus * (1.0 - poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
}

/// The **bending (flexure) stress** `σ = M·y/I` (Pa) — the longitudinal normal stress
/// a bending moment `moment` `M` (N·m) induces at distance `distance_from_neutral_axis`
/// `y` (m) from the neutral axis in a section of second moment of area
/// `second_moment_area` `I` (m⁴). It varies linearly across the section: zero at the
/// neutral axis, tensile on the convex face and compressive on the concave one, with
/// the design-critical maximum at the extreme fibre `y = c` (checked against the
/// material's yield strength). It is the stress conjugate of the curvature
/// [`beam_curvature`] through Hooke's law `σ = E·κ·y`. Returns `0` for non-physical
/// input (`M` or `y` non-finite, or `I` non-positive or non-finite).
pub fn bending_stress(moment: f64, distance_from_neutral_axis: f64, second_moment_area: f64) -> f64 {
    if !moment.is_finite()
        || !distance_from_neutral_axis.is_finite()
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    moment * distance_from_neutral_axis / second_moment_area
}

/// The **transverse (flexural) shear stress** `τ = V·Q / (I·b)` (Pa) — Jourawski's formula
/// for the shear in a beam cross-section under a transverse shear force `shear_force` `V`
/// (N), where `first_moment_area` `Q` (m³) is the first moment about the neutral axis of the
/// area beyond the cut, `second_moment_area` `I` (m⁴) the section's second moment, and
/// `width` `b` (m) the section width at the cut. It is the third classic beam stress
/// alongside the bending normal stress [`bending_stress`] (`σ = M·y/I`) and the
/// [`torsional_shear_stress`] (`τ = T·r/J`); for a rectangular section it peaks at the
/// neutral axis at `τ_max = 1.5·V/A`, half again the average shear `V/A`. Returns `0` for
/// non-physical input (`I` or `b` non-positive, or any non-finite argument).
pub fn beam_transverse_shear_stress(
    shear_force: f64,
    first_moment_area: f64,
    second_moment_area: f64,
    width: f64,
) -> f64 {
    if !shear_force.is_finite()
        || !first_moment_area.is_finite()
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
        || !width.is_finite()
        || width <= 0.0
    {
        return 0.0;
    }
    shear_force * first_moment_area / (second_moment_area * width)
}

/// The **elastic section modulus** `S = I/c` (m³) — the cross-section property (the `W`
/// or `Z` of structural beam tables) that turns a bending moment directly into the peak
/// fibre stress, `σ_max = M/S`. `second_moment_area` `I` (m⁴) is the section's second
/// moment of area and `extreme_fiber_distance` `c` (m) the distance from the neutral
/// axis to the most-stressed (outermost) fibre. A larger `S` carries a given moment at a
/// lower stress, so it is the single number beam selection optimises. Returns `0` for
/// non-physical input (`I` or `c` non-positive or non-finite).
pub fn elastic_section_modulus(second_moment_area: f64, extreme_fiber_distance: f64) -> f64 {
    if !second_moment_area.is_finite()
        || second_moment_area <= 0.0
        || !extreme_fiber_distance.is_finite()
        || extreme_fiber_distance <= 0.0
    {
        return 0.0;
    }
    second_moment_area / extreme_fiber_distance
}

/// The **second moment of area of a rectangular cross-section** about its centroidal
/// bending axis `I = b·h³/12` (m⁴), for a section of `width` `b` (m) and `height` `h`
/// (m). It is the foundational section property of the commonest beam section — the `I`
/// that feeds [`elastic_section_modulus`] (`S = I/(h/2)`), [`flexural_rigidity`] (`EI`),
/// and every Euler–Bernoulli deflection. The bending axis is normal to `height`, so the
/// depth enters cubed (a deeper beam is far stiffer). Returns `0` for non-physical input
/// (`b` or `h` non-positive or non-finite).
pub fn rectangular_second_moment_of_area(width: f64, height: f64) -> f64 {
    if !width.is_finite() || width <= 0.0 || !height.is_finite() || height <= 0.0 {
        return 0.0;
    }
    width * height.powi(3) / 12.0
}

/// The **plastic section modulus of a rectangular cross-section** `Z = b·h²/4` (m³), for a
/// section of `width` `b` (m) and `height` `h` (m) — the combined first moment of the two
/// half-areas about the equal-area (plastic neutral) axis, the section property for
/// fully-plastic limit bending (`M_plastic = σ_y·Z`). It is `1.5×` the rectangular *elastic*
/// section modulus ([`elastic_section_modulus`] of [`rectangular_second_moment_of_area`],
/// `S = b·h²/6`) — the rectangle's shape factor `Z/S = 1.5`, the reserve between first yield
/// and a full plastic hinge. Returns `0` for non-physical input (`b` or `h` non-positive or
/// non-finite).
pub fn rectangular_plastic_section_modulus(width: f64, height: f64) -> f64 {
    if !width.is_finite() || width <= 0.0 || !height.is_finite() || height <= 0.0 {
        return 0.0;
    }
    width * height * height / 4.0
}

/// The **polar second moment of area of a rectangular cross-section** about the centroidal
/// axis normal to the section `J = I_x + I_y = (b·h/12)·(b² + h²)` (m⁴), for a section of
/// `width` `b` (m) and `height` `h` (m) — by the perpendicular-axis theorem the sum of the
/// two bending [`rectangular_second_moment_of_area`]s. **NB:** this is the polar *second
/// moment of area*; for a non-circular section it is **not** the St-Venant torsion constant
/// (a rectangle's torsion constant is the separate, smaller `β·b·h³` value). `J` equals the
/// torsion constant only for a circle. Returns `0` for non-physical input (`b` or `h`
/// non-positive or non-finite).
pub fn rectangular_polar_second_moment_of_area(width: f64, height: f64) -> f64 {
    if !width.is_finite() || width <= 0.0 || !height.is_finite() || height <= 0.0 {
        return 0.0;
    }
    width * height * (width * width + height * height) / 12.0
}

/// The **second moment of area of a solid circular cross-section** about a centroidal
/// diameter `I = π·d⁴/64` (m⁴), for a round bar of diameter `diameter` `d` (m) — the
/// shaft/rod companion to [`rectangular_second_moment_of_area`], the `I` that feeds
/// [`elastic_section_modulus`] (`S = I/(d/2) = π·d³/32`) and [`flexural_rigidity`] (`EI`).
/// The diameter enters to the fourth power. Returns `0` for non-physical input (`d`
/// non-positive or non-finite).
pub fn circular_second_moment_of_area(diameter: f64) -> f64 {
    if !diameter.is_finite() || diameter <= 0.0 {
        return 0.0;
    }
    std::f64::consts::PI * diameter.powi(4) / 64.0
}

/// The **second moment of area of a hollow circular (tube/pipe) cross-section** about a
/// centroidal diameter `I = π·(D⁴ − d⁴)/64` (m⁴), for outer diameter `outer_diameter` `D`
/// (m) and inner bore `inner_diameter` `d` (m) — the annulus is the outer disc minus the
/// inner disc, the bending `I` of structural tubing (far stiffer per unit weight than a
/// solid bar, since the removed core carries little bending stress). Returns `0` for
/// non-physical input (`D` non-positive or non-finite, or the bore `d` negative, non-finite,
/// or not strictly inside the outer diameter).
pub fn hollow_circular_second_moment_of_area(outer_diameter: f64, inner_diameter: f64) -> f64 {
    if !outer_diameter.is_finite()
        || outer_diameter <= 0.0
        || !inner_diameter.is_finite()
        || inner_diameter < 0.0
        || inner_diameter >= outer_diameter
    {
        return 0.0;
    }
    std::f64::consts::PI * (outer_diameter.powi(4) - inner_diameter.powi(4)) / 64.0
}

/// The **polar second moment of area of a hollow circular shaft** about its longitudinal axis
/// `J = π·(D⁴ − d⁴)/32` (m⁴), for outer diameter `outer_diameter` `D` (m) and inner bore
/// `inner_diameter` `d` (m). For a *circular* tube this polar second moment **is** the
/// St-Venant torsion constant — the torsional stiffness of drive-shaft tubing, which (like the
/// solid shaft) is exactly twice the bending [`hollow_circular_second_moment_of_area`] by the
/// perpendicular-axis theorem, and the annulus = the outer disc's polar moment minus the
/// inner's. Returns `0` for non-physical input (`D` non-positive or non-finite, or the bore
/// `d` negative, non-finite, or not strictly inside the outer diameter).
pub fn hollow_circular_polar_second_moment_of_area(
    outer_diameter: f64,
    inner_diameter: f64,
) -> f64 {
    if !outer_diameter.is_finite()
        || outer_diameter <= 0.0
        || !inner_diameter.is_finite()
        || inner_diameter < 0.0
        || inner_diameter >= outer_diameter
    {
        return 0.0;
    }
    std::f64::consts::PI * (outer_diameter.powi(4) - inner_diameter.powi(4)) / 32.0
}

/// The **second moment of area of a hollow rectangular (box) cross-section** about its
/// centroidal bending axis `I = (b·h³ − bᵢ·hᵢ³)/12` (m⁴), for outer width `width` `b` (m),
/// outer height `height` `h` (m), inner width `inner_width` `bᵢ` (m) and inner height
/// `inner_height` `hᵢ` (m). The box tube is the canonical stiffness-to-weight section — the
/// removed core carries little bending stress — and is the rectangular companion to
/// [`hollow_circular_second_moment_of_area`]. The bending axis is normal to `height`, so the
/// depth enters cubed. Returns `0` for non-physical input (outer dimensions non-positive or
/// non-finite, inner dimensions negative or non-finite, or the inner bore not strictly inside
/// the outer envelope).
pub fn hollow_rectangular_second_moment_of_area(
    width: f64,
    height: f64,
    inner_width: f64,
    inner_height: f64,
) -> f64 {
    if !width.is_finite()
        || width <= 0.0
        || !height.is_finite()
        || height <= 0.0
        || !inner_width.is_finite()
        || inner_width < 0.0
        || !inner_height.is_finite()
        || inner_height < 0.0
        || inner_width >= width
        || inner_height >= height
    {
        return 0.0;
    }
    (width * height.powi(3) - inner_width * inner_height.powi(3)) / 12.0
}

/// The **polar second moment of area of a hollow rectangular (box) cross-section**
/// `J = (b·h·(b²+h²) − bᵢ·hᵢ·(bᵢ²+hᵢ²))/12` (m⁴) — the torsion-relevant polar moment of a
/// box tube, equal to `Ix + Iy` (the sum of the bending second moments about both centroidal
/// axes), for outer width `width` `b` (m), outer height `height` `h` (m), inner width
/// `inner_width` `bᵢ` (m) and inner height `inner_height` `hᵢ` (m). It is the polar companion
/// to [`hollow_rectangular_second_moment_of_area`] (bending) and to
/// [`hollow_circular_polar_second_moment_of_area`] (round tube), and equals `J_outer − J_inner`
/// where `J_rect` is [`rectangular_polar_second_moment_of_area`]. For a square box (`b = h`,
/// `bᵢ = hᵢ`) it is exactly twice the bending second moment (`Ix = Iy`). Returns `0` for
/// non-physical input (outer dimensions non-positive or non-finite, inner dimensions negative
/// or non-finite, or the inner bore not strictly inside the outer envelope).
pub fn hollow_rectangular_polar_second_moment_of_area(
    width: f64,
    height: f64,
    inner_width: f64,
    inner_height: f64,
) -> f64 {
    if !width.is_finite()
        || width <= 0.0
        || !height.is_finite()
        || height <= 0.0
        || !inner_width.is_finite()
        || inner_width < 0.0
        || !inner_height.is_finite()
        || inner_height < 0.0
        || inner_width >= width
        || inner_height >= height
    {
        return 0.0;
    }
    (width * height * (width * width + height * height)
        - inner_width * inner_height * (inner_width * inner_width + inner_height * inner_height))
        / 12.0
}

/// The **plastic section modulus of a solid circular section** `Z = d³/6` (m³), for a round
/// bar of diameter `diameter` `d` (m) — the section property for fully-plastic limit bending
/// of a shaft (`M_plastic = σ_y·Z`), the circular companion to
/// [`rectangular_plastic_section_modulus`]. For a solid circle the shape factor is
/// `Z/S = 16/(3π) ≈ 1.698` (`Z = d³/6` vs the elastic `S = π·d³/32`, the
/// [`elastic_section_modulus`] of [`circular_second_moment_of_area`]). Returns `0` for
/// non-physical input (`d` non-positive or non-finite).
pub fn circular_plastic_section_modulus(diameter: f64) -> f64 {
    if !diameter.is_finite() || diameter <= 0.0 {
        return 0.0;
    }
    diameter * diameter * diameter / 6.0
}

/// The **polar second moment of area of a solid circular shaft** about its longitudinal
/// axis `J = π·d⁴/32` (m⁴), for a round bar of diameter `diameter` `d` (m) — the **torsion
/// constant** that feeds [`polar_section_modulus`] (`Z_p = J/(d/2)`) and
/// [`torsional_rigidity`] (`GJ`). By the perpendicular-axis theorem it is exactly twice the
/// bending [`circular_second_moment_of_area`] (`J = I_x + I_y = 2·I`, since `I_x = I_y` for
/// a circle). Returns `0` for non-physical input (`d` non-positive or non-finite).
pub fn circular_polar_second_moment_of_area(diameter: f64) -> f64 {
    if !diameter.is_finite() || diameter <= 0.0 {
        return 0.0;
    }
    std::f64::consts::PI * diameter.powi(4) / 32.0
}

/// The **elastic bending-moment capacity** `M = σ·S` (N·m) — the bending moment a
/// section carries when its most-stressed fibre reaches the stress `stress` `σ` (Pa),
/// given the [`elastic_section_modulus`] `section_modulus` `S` (m³). It is the inverse
/// of the flexure formula `σ_max = M/S` ([`bending_stress`] at the extreme fibre):
/// feeding the *allowable* stress gives the section's allowable working moment, and
/// feeding the *yield* stress gives the first-yield moment `M_y`. Linear in both the
/// stress and the section modulus. Returns `0` for non-physical input (non-finite, or a
/// non-positive section modulus).
pub fn bending_moment_capacity(stress: f64, section_modulus: f64) -> f64 {
    if !stress.is_finite() || !section_modulus.is_finite() || section_modulus <= 0.0 {
        return 0.0;
    }
    stress * section_modulus
}

/// The analytic **simply-supported (pinned–pinned) mid-span deflection**
/// `δ = P·L³/(48·E·I)` (m) of a slender Euler–Bernoulli beam under a *central*
/// transverse point load `load` `P` (N) — span `length` `L` (m), Young's modulus
/// `youngs_modulus` `E` (Pa), section second moment of area `second_moment_area`
/// `I` (m⁴).
///
/// This is the *simply-supported* bending reference the finite-element solver
/// ([`solve_beam_static`]) converges to — the boundary-condition companion to the
/// clamped–free [`cantilever_tip_deflection`]. Same load type, different
/// supports: pinning *both* ends makes the beam **16× stiffer** at mid-span than a
/// cantilever of the same length is at its tip (the `1/48` coefficient versus
/// `1/3`), so `δ_ss = `[`cantilever_tip_deflection`]`(P, L, E, I)/16`. It grows
/// linearly with the load `P` (and is sign-preserving), with the cube of the span
/// `L`, and falls inversely with the flexural rigidity `E·I`. Returns `0` for
/// non-physical input (`P` non-finite, or `E`, `I`, or `L` non-positive or
/// non-finite).
pub fn simply_supported_center_deflection(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * length.powi(3) / (48.0 * youngs_modulus * second_moment_area)
}

/// The analytic **maximum bending moment of a simply-supported beam under a central
/// point load** `M = P·L/4` (N·m) — the peak moment, at mid-span beneath the load,
/// which sets the maximum bending stress and so governs the strength design of the
/// member. `load` `P` is the transverse load (N) and `length` `L` the span (m).
///
/// It is the simply-supported companion to the clamped–clamped
/// [`fixed_fixed_point_load_end_moment`] (`P·L/8`): building in the ends halves the
/// peak moment, redistributing it so the same magnitude `P·L/8` appears at both
/// fixed ends *and* at mid-span instead of the bare `P·L/4` of the pinned span.
/// Linear and sign-preserving in `P`, linear in `L`, and — being a statics result —
/// independent of `E` and `I`. Returns `0` for non-physical input (`P` non-finite,
/// or `L` non-positive or non-finite).
pub fn simply_supported_point_load_max_moment(load: f64, length: f64) -> f64 {
    if !load.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load * length / 4.0
}

/// The analytic **maximum bending moment of a simply-supported beam under a
/// uniformly distributed load** `M = w·L²/8` (N·m) — the peak moment, at mid-span,
/// which sets the maximum bending stress and so governs the strength design of the
/// member. `load_per_length` `w` is the load intensity (N/m) and `length` `L` the
/// span (m).
///
/// It is the UDL companion to the point-load
/// [`simply_supported_point_load_max_moment`] (`P·L/4`), completing the
/// simply-supported peak-moment pair. It is `3/2` of the clamped–clamped fixing
/// moment [`fixed_fixed_udl_end_moment`] (`w·L²/12`) — building in the ends sheds
/// a third of the peak moment. Quadratic in `L`, linear and sign-preserving in `w`,
/// and — a statics result — independent of `E` and `I`. Returns `0` for non-physical
/// input (`w` non-finite, or `L` non-positive or non-finite).
pub fn simply_supported_udl_max_moment(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load_per_length * length * length / 8.0
}

/// The analytic **fixed-fixed (clamped–clamped) beam mid-span deflection**
/// `δ = P·L³/(192·E·I)` (m) — the deflection at mid-span of a slender
/// Euler–Bernoulli beam **built in at both ends** (encastré) carrying a transverse
/// point force `load` `P` (N) at the centre, with span `length` `L` (m), Young's
/// modulus `youngs_modulus` `E` (Pa), and section second moment of area
/// `second_moment_area` `I` (m⁴).
///
/// This is the third classic single-span boundary condition, completing the set
/// with the clamped–free [`cantilever_tip_deflection`] (`P·L³/3EI`) and the
/// pinned–pinned [`simply_supported_center_deflection`] (`P·L³/48EI`). Clamping
/// *both* ends adds end-fixing moments that make the span the stiffest of the
/// three: it deflects exactly **¼** as far as the simply-supported beam, and
/// **1/64** as far as a cantilever of the same length, under the same central
/// load. The deflection is linear (and sign-preserving) in `P`, grows with the
/// cube of the span, and falls inversely with the flexural rigidity `E·I`. Returns
/// `0` for non-physical input (`P` non-finite, or `E`, `I`, or `L` non-positive or
/// non-finite).
pub fn fixed_fixed_center_deflection(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * length.powi(3) / (192.0 * youngs_modulus * second_moment_area)
}

/// The analytic **fixed-fixed (clamped–clamped) beam mid-span deflection under a
/// uniformly distributed load** `δ = w·L⁴/(384·E·I)` (m) — the mid-span deflection
/// of a slender Euler–Bernoulli beam **built in at both ends** carrying a transverse
/// load of intensity `load_per_length` `w` (N/m) over its full span `length` `L`
/// (m), with Young's modulus `youngs_modulus` `E` (Pa) and section second moment of
/// area `second_moment_area` `I` (m⁴).
///
/// It is the distributed-load companion to the point-load
/// [`fixed_fixed_center_deflection`] (`P·L³/192EI`), and the clamped–clamped
/// counterpart to the pinned–pinned [`simply_supported_udl_center_deflection`]
/// (`5·w·L⁴/384EI`). Clamping *both* ends makes the span exactly **5× stiffer** at
/// mid-span than pinning both under the same UDL. The deflection is linear (and
/// sign-preserving) in `w`, grows with the *fourth* power of the span, and falls
/// inversely with the flexural rigidity `E·I`. Returns `0` for non-physical input
/// (`w` non-finite, or `E`, `I`, or `L` non-positive or non-finite).
pub fn fixed_fixed_udl_center_deflection(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load_per_length * length.powi(4) / (384.0 * youngs_modulus * second_moment_area)
}

/// The analytic **fixed-end (fixing) moment** `M = P·L/8` (N·m) of a slender
/// Euler–Bernoulli beam **built in at both ends** (clamped–clamped) under a central
/// transverse point force `load` `P` (N) on a span `length` `L` (m) — the reaction
/// moment the clamps must supply to hold the ends at zero slope.
///
/// This is the classic **fixed-end moment** of the slope-deflection and
/// moment-distribution methods, and the structural complement to the deflection
/// [`fixed_fixed_center_deflection`]: it is exactly the end moment that, superposed
/// on a pinned span, stiffens it into a clamped one (lifting mid-span by
/// `M·L²/8EI`, the difference between the pinned and clamped centre deflections).
/// It is linear and sign-preserving in `P` and linear in `L`, and — being a pure
/// statics result — is independent of `E` and `I`. Returns `0` for non-physical
/// input (`P` non-finite, or `L` non-positive or non-finite).
pub fn fixed_fixed_point_load_end_moment(load: f64, length: f64) -> f64 {
    if !load.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load * length / 8.0
}

/// The analytic **fixed-end (fixing) moment under a uniformly distributed load**
/// `M = w·L²/12` (N·m) of a clamped–clamped beam — the reaction moment each built-in
/// end supplies to hold its slope at zero under a transverse load of intensity
/// `load_per_length` `w` (N/m) over a span `length` `L` (m).
///
/// It is the UDL companion to the point-load [`fixed_fixed_point_load_end_moment`]
/// (`P·L/8`), completing the clamped–clamped fixed-end-moment pair of the
/// slope-deflection and moment-distribution methods. Like its sibling it is the end
/// moment that, superposed on a pinned span, stiffens it into a clamped one (lifting
/// mid-span by `M·L²/8EI`, the difference between the pinned and clamped UDL centre
/// deflections). Quadratic in `L`, linear and sign-preserving in `w`, and — being a
/// pure statics result — independent of `E` and `I`. Returns `0` for non-physical
/// input (`w` non-finite, or `L` non-positive or non-finite).
pub fn fixed_fixed_udl_end_moment(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load_per_length * length * length / 12.0
}

/// The analytic **strain energy of a fixed–fixed (clamped–clamped) beam under a central
/// point load** `U = P²·L³/(384·E·I)` (J) — the elastic bending energy stored in a slender
/// Euler–Bernoulli beam built in at both ends carrying a transverse point force `load` `P`
/// (N) at mid-span, with span `length` `L` (m), Young's modulus `youngs_modulus` `E` (Pa),
/// and section second moment of area `second_moment_area` `I` (m⁴).
///
/// By Clapeyron's theorem `U = ½·P·δ` with the centre deflection
/// [`fixed_fixed_center_deflection`] `δ = P·L³/(192EI)`. It completes the point-load
/// strain-energy set: clamping both ends stores **1/64** the energy of a cantilever
/// [`cantilever_point_load_strain_energy`] (`P²L³/6EI`) and **1/4** of a simply-supported
/// beam [`simply_supported_point_load_strain_energy`] (`P²L³/96EI`) under the same central
/// load. The energy grows with the square of the load (sign-independent), the cube of the
/// span, and falls inversely with the flexural rigidity `E·I`. Returns `0` for non-physical
/// input (`P` non-finite, or `E`, `I`, or `L` non-positive or non-finite).
pub fn fixed_fixed_point_load_strain_energy(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * load * length.powi(3) / (384.0 * youngs_modulus * second_moment_area)
}

/// The analytic **prop reaction of a propped cantilever under a uniformly distributed
/// load** `R_B = 3·w·L/8` (N) — the vertical reaction at the *propped* (simple-support)
/// end of a propped cantilever (fixed at one end, simply supported at the other) carrying
/// a uniform load `load_per_length` `w` (N/m) over span `length` `L` (m). This is the
/// classic statically-indeterminate beam released by the **force (compatibility) method**:
/// removing the prop leaves a cantilever that deflects `w·L⁴/8EI` at the tip
/// ([`cantilever_udl_tip_deflection`]); the prop supplies exactly the upward point force
/// that pushes the tip back to zero ([`cantilever_tip_deflection`] `R·L³/3EI`), giving
/// `R_B = 3wL/8`. The remaining `5wL/8` goes to the fixed end. Being a pure statics result
/// it is linear and sign-preserving in `w`, linear in `L`, and independent of `E` and `I`.
/// Returns `0` for non-physical input (`w` non-finite, or `L` non-positive or non-finite).
pub fn propped_cantilever_udl_prop_reaction(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    3.0 * load_per_length * length / 8.0
}

/// The analytic **fixed-end (clamp) moment of a propped cantilever under a uniformly
/// distributed load** `M_A = w·L²/8` (N·m) — the hogging bending moment that develops at
/// the *fixed* (clamped) end of a propped cantilever (fixed at one end, simply supported
/// at the other) carrying a uniform load `load_per_length` `w` (N/m) over span `length`
/// `L` (m). It is the design-critical peak moment of the case. By statics it is the
/// complement of the [`propped_cantilever_udl_prop_reaction`] `R_B = 3wL/8`: taking
/// moments about the fixed end, `M_A = w·L²/2 − R_B·L = wL²/2 − 3wL²/8 = wL²/8`.
///
/// Note this evaluates to `w·L²/8`, numerically the same as the simply-supported midspan
/// moment [`simply_supported_udl_max_moment`] — a genuine coincidence; the two are
/// physically distinct (this is the *hogging* moment at the *clamp* of a propped
/// cantilever, not the *sagging* midspan moment of a pinned–pinned beam). Quadratic in
/// `L`, linear and sign-preserving in `w`, and — being pure statics — independent of `E`
/// and `I`. Returns `0` for non-physical input (`w` non-finite, or `L` non-positive or
/// non-finite).
pub fn propped_cantilever_udl_fixed_end_moment(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    load_per_length * length * length / 8.0
}

/// The **middle-support hogging moment of a two-span continuous beam under a uniformly
/// distributed load** `M_B = w·L²/8` (N·m) — the bending moment over the central support
/// `B` of a beam on three simple supports `A`–`B`–`C` with two equal spans `span_length`
/// `L` (m) each carrying a uniform load `load_per_length` `w` (N/m). It is the classic
/// **three-moment-theorem (Clapeyron)** result and the strength-governing moment of the
/// arrangement.
///
/// By symmetry the centre support acts as a fixed (clamped) end for each span, so every
/// span is itself a **propped cantilever** under UDL — hence `M_B` is numerically equal to
/// the [`propped_cantilever_udl_fixed_end_moment`] `w·L²/8`, a different configuration
/// reaching the same value through the symmetry of the continuous beam. Being pure statics
/// it is independent of `E` and `I`, linear and sign-preserving in `w`, and quadratic in
/// `L`. Returns `0` for non-physical input (`w` non-finite, or `L` non-positive or
/// non-finite).
pub fn two_span_continuous_beam_udl_middle_moment(load_per_length: f64, span_length: f64) -> f64 {
    if !load_per_length.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    load_per_length * span_length * span_length / 8.0
}

/// The **middle-support reaction of a two-span continuous beam under a uniformly
/// distributed load** `R_B = 5·w·L/4` (N) — the vertical reaction at the central support
/// `B` of a beam on three simple supports `A`–`B`–`C` with two equal spans `span_length`
/// `L` (m) each carrying a uniform load `load_per_length` `w` (N/m). It is the
/// heaviest-loaded support of the arrangement.
///
/// By the same symmetry that fixes the [`two_span_continuous_beam_udl_middle_moment`], each
/// span behaves as a propped cantilever clamped at `B`, so the centre support collects the
/// fixed-end reaction of *both* spans: `R_B = 2 · 5wL/8 = 5wL/4`, exactly twice the
/// [`propped_cantilever_udl_fixed_end_reaction`]. With the two outer reactions
/// `R_A = R_C = 3wL/8` (each the [`propped_cantilever_udl_prop_reaction`]) it balances the
/// total load `R_A + R_B + R_C = 2wL`. Being pure statics it is independent of `E` and `I`,
/// and linear in both `w` and `L`. Returns `0` for non-physical input (`w` non-finite, or
/// `L` non-positive or non-finite).
pub fn two_span_continuous_beam_udl_middle_reaction(load_per_length: f64, span_length: f64) -> f64 {
    if !load_per_length.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    5.0 * load_per_length * span_length / 4.0
}

/// The **outer-support reaction of a two-span continuous beam under a uniformly
/// distributed load** `R_A = R_C = 3·w·L/8` (N) — the vertical reaction at each of the two
/// end simple supports `A` and `C` of a beam on three simple supports `A`–`B`–`C` with two
/// equal spans `span_length` `L` (m) each carrying a uniform load `load_per_length` `w`
/// (N/m).
///
/// By symmetry each end support is the *propped* (simple) end of its span's equivalent
/// propped cantilever, so `R_A = R_C = 3wL/8`, exactly the
/// [`propped_cantilever_udl_prop_reaction`]. Together with the centre
/// [`two_span_continuous_beam_udl_middle_reaction`] `R_B = 5wL/4` the three reactions carry
/// the whole load: `2·R_A + R_B = 2wL`. Being pure statics it is independent of `E` and
/// `I`, and linear in both `w` and `L`. Returns `0` for non-physical input (`w` non-finite,
/// or `L` non-positive or non-finite).
pub fn two_span_continuous_beam_udl_outer_reaction(load_per_length: f64, span_length: f64) -> f64 {
    if !load_per_length.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    3.0 * load_per_length * span_length / 8.0
}

/// The **interior-support moment of a three-span continuous beam under a uniformly
/// distributed load** `M = w·L²/10` (N·m) — the bending moment over each of the two interior
/// supports of a beam on four simple supports `A`–`B`–`C`–`D` with three equal spans
/// `span_length` `L` (m), each carrying a uniform load `load_per_length` `w` (N/m).
///
/// By the three-moment (Clapeyron) theorem with `M_A = M_D = 0` and symmetry `M_B = M_C`,
/// the support equation `4·M_B·L + M_C·L = −w·L³/2` (free-BM area `wL³/12` per span) gives
/// `5·M_B·L = −wL³/2`, hence `M_B = −wL²/10`. The magnitude `wL²/10` is **smaller** than the
/// two-span [`two_span_continuous_beam_udl_middle_moment`] `wL²/8` (ratio 8:10), since the
/// extra span shares the load. Pure statics: independent of `E` and `I`, linear in `w`,
/// quadratic in `L`. Returns `0` for non-physical input (`w` non-finite, or `L` non-positive
/// or non-finite).
pub fn three_span_continuous_beam_udl_interior_moment(load_per_length: f64, span_length: f64) -> f64 {
    if !load_per_length.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    load_per_length * span_length * span_length / 10.0
}

/// The **end-support reaction of a three-span continuous beam under a uniformly distributed
/// load** `R_A = R_D = 2·w·L/5` (N) — the reaction at each of the two end simple supports
/// `A` and `D` of a beam on four simple supports `A`–`B`–`C`–`D` with three equal spans
/// `span_length` `L` (m), each carrying a uniform load `load_per_length` `w` (N/m).
///
/// From end-span statics with the interior moment `M_B = −wL²/10`
/// ([`three_span_continuous_beam_udl_interior_moment`]): `R_A = wL/2 − wL/10 = 2wL/5`. With
/// the interior [`three_span_continuous_beam_udl_interior_reaction`] `R_int = 11wL/10` the
/// four reactions balance the total load `2·R_A + 2·R_int = 3wL`. Pure statics: independent
/// of `E` and `I`, linear in `w` and `L`. Returns `0` for non-physical input (`w` non-finite,
/// or `L` non-positive or non-finite).
pub fn three_span_continuous_beam_udl_end_reaction(load_per_length: f64, span_length: f64) -> f64 {
    if !load_per_length.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    2.0 * load_per_length * span_length / 5.0
}

/// The **interior-support reaction of a three-span continuous beam under a uniformly
/// distributed load** `R_B = R_C = 11·w·L/10` (N) — the reaction at each of the two interior
/// simple supports `B` and `C` of a beam on four simple supports `A`–`B`–`C`–`D` with three
/// equal spans `span_length` `L` (m), each carrying a uniform load `load_per_length` `w`
/// (N/m).
///
/// The interior supports carry the heaviest share. By global vertical equilibrium with the
/// end reactions [`three_span_continuous_beam_udl_end_reaction`] `R_A = 2wL/5`, the four
/// reactions sum to the total load `2·R_A + 2·R_B = 3wL`, giving `R_B = 11wL/10`. Pure
/// statics: independent of `E` and `I`, linear in `w` and `L`. Returns `0` for non-physical
/// input (`w` non-finite, or `L` non-positive or non-finite).
pub fn three_span_continuous_beam_udl_interior_reaction(load_per_length: f64, span_length: f64) -> f64 {
    if !load_per_length.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    11.0 * load_per_length * span_length / 10.0
}

/// The **middle-support moment of a two-span continuous beam under a central point load in
/// one span** `M_B = 3PL/32` (N·m) — the bending moment over the central support `B` of a
/// beam on three simple supports `A`–`B`–`C` with two equal spans `span_length` `L` (m)
/// when a single transverse point force `point_load` `P` (N) acts at the **mid-span of one
/// span** (the other span unloaded).
///
/// It is the point-load companion to the distributed-load
/// [`two_span_continuous_beam_udl_middle_moment`] `wL²/8`, from the same three-moment
/// (Clapeyron) theorem: the loaded span's free bending-moment triangle (area `PL²/8`,
/// centroid `L/2` from the end support) gives `2·M_B·(2L) = −6·(PL²/8)(L/2)/L`, i.e.
/// `M_B = 3PL/32`. Equivalently it is **half** the propped-cantilever clamping moment
/// [`propped_cantilever_central_load_fixed_end_moment`] `3PL/16` — the *unloaded* equal-
/// stiffness adjacent span provides finite (not rigid) rotational restraint at `B`, halving
/// the fully-fixed value. Being pure statics it is independent of `E` and `I`, linear and
/// sign-preserving in `P`, and linear in `L`. Returns `0` for non-physical input (`P`
/// non-finite, or `L` non-positive or non-finite).
pub fn two_span_continuous_beam_central_point_load_middle_moment(
    point_load: f64,
    span_length: f64,
) -> f64 {
    if !point_load.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    3.0 * point_load * span_length / 32.0
}

/// The **loaded-span outer reaction of a two-span continuous beam under a central point
/// load in one span** `R_A = 13P/32` (N) — the reaction at the end simple support `A` of
/// the *loaded* span (`A`–`B`) of a beam on three simple supports `A`–`B`–`C` (two equal
/// spans `span_length` `L`), carrying a point force `point_load` `P` (N) at the mid-span
/// of span `A`–`B` (span `B`–`C` unloaded).
///
/// Derived from the centre moment
/// [`two_span_continuous_beam_central_point_load_middle_moment`] `M_B = 3PL/32` (hogging)
/// by span statics (`R_A·L = P·L/2 − M_B ⇒ R_A = 13P/32`). With the centre
/// [`two_span_continuous_beam_central_point_load_middle_reaction`] `R_B = 11P/16` and the
/// unloaded-span [`two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction`]
/// `R_C = −3P/32` it satisfies vertical equilibrium `R_A + R_B + R_C = P`. Pure statics:
/// independent of `E`, `I`, and `L`; linear and sign-preserving in `P`. Returns `0` for
/// non-physical input (`P` non-finite, or `L` non-positive or non-finite).
pub fn two_span_continuous_beam_central_point_load_loaded_span_outer_reaction(
    point_load: f64,
    span_length: f64,
) -> f64 {
    if !point_load.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    13.0 * point_load / 32.0
}

/// The **centre-support reaction of a two-span continuous beam under a central point load
/// in one span** `R_B = 11P/16` (N) — the reaction at the middle support `B` of a beam on
/// three simple supports `A`–`B`–`C` (two equal spans `span_length` `L`), with a point
/// force `point_load` `P` (N) at the mid-span of one span.
///
/// It is the dual of the centre moment
/// [`two_span_continuous_beam_central_point_load_middle_moment`] `M_B = 3PL/32`, and the
/// largest of the three reactions (`11/16 ≈ 69%` of `P`): the middle support collects the
/// inner shears of both spans. With the loaded-span `R_A = 13P/32` and the unloaded-span
/// `R_C = −3P/32` it satisfies `R_A + R_B + R_C = P`. Pure statics: independent of `E`,
/// `I`, and `L`; linear and sign-preserving in `P`. Returns `0` for non-physical input
/// (`P` non-finite, or `L` non-positive or non-finite).
pub fn two_span_continuous_beam_central_point_load_middle_reaction(
    point_load: f64,
    span_length: f64,
) -> f64 {
    if !point_load.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    11.0 * point_load / 16.0
}

/// The **unloaded-span outer reaction of a two-span continuous beam under a central point
/// load in one span** `R_C = −3P/32` (N) — the reaction at the end simple support `C` of
/// the *unloaded* span (`B`–`C`) of a beam on three simple supports `A`–`B`–`C` (two equal
/// spans `span_length` `L`), with a point force `point_load` `P` (N) at the mid-span of the
/// other span (`A`–`B`).
///
/// The **negative sign is uplift**: the hogging moment at `B` lifts the far support, so for
/// a downward load the unloaded end pulls *up* (a classic continuous-beam result, and why
/// the far support of a propped span can need hold-down). It is the smallest of the three
/// reactions. With the loaded-span `R_A = 13P/32` and the centre `R_B = 11P/16` it satisfies
/// `R_A + R_B + R_C = P`. Pure statics: independent of `E`, `I`, and `L`; magnitude linear
/// in `P` (so it reverses sign with the load). Returns `0` for non-physical input (`P`
/// non-finite, or `L` non-positive or non-finite).
pub fn two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction(
    point_load: f64,
    span_length: f64,
) -> f64 {
    if !point_load.is_finite() || !span_length.is_finite() || span_length <= 0.0 {
        return 0.0;
    }
    -3.0 * point_load / 32.0
}

/// The analytic **fixed-end (clamp) reaction of a propped cantilever under a uniformly
/// distributed load** `R_A = 5·w·L/8` (N) — the vertical support reaction at the *fixed*
/// (clamped) end of a propped cantilever (fixed at one end, simply supported at the other)
/// carrying a uniform load `load_per_length` `w` (N/m) over span `length` `L` (m). With the
/// [`propped_cantilever_udl_prop_reaction`] `R_B = 3wL/8` it carries the whole applied
/// load — `R_A + R_B = w·L` (vertical equilibrium) — the clamp taking the larger `5/8`
/// share because it also resists the [`propped_cantilever_udl_fixed_end_moment`]. Being a
/// pure statics result it is linear and sign-preserving in `w`, linear in `L`, and
/// independent of `E` and `I`. Returns `0` for non-physical input (`w` non-finite, or `L`
/// non-positive or non-finite).
pub fn propped_cantilever_udl_fixed_end_reaction(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    5.0 * load_per_length * length / 8.0
}

/// The analytic **maximum sagging (span) moment of a propped cantilever under a
/// uniformly distributed load** `M_sag = 9·w·L²/128` (N·m) — the largest *positive*
/// (sagging) bending moment in the span of a propped cantilever (fixed at one end,
/// simply supported at the other) under a uniform load `load_per_length` `w` (N/m) over
/// span `length` `L` (m). It occurs at the point of zero shear, `x = 5L/8` from the fixed
/// end (`= 3L/8` from the prop). Together with the fixed-end hogging moment
/// [`propped_cantilever_udl_fixed_end_moment`] `M_A = wL²/8 = 16wL²/128` it defines the
/// full bending-moment envelope: the clamp (`|M_A| = 16wL²/128`) governs strength, while
/// this span peak (`9wL²/128`) is the maximum sagging value the bottom fibre sees. By
/// statics from the prop end it is `M_sag = R_B·(3L/8) − w·(3L/8)²/2` with the prop
/// reaction [`propped_cantilever_udl_prop_reaction`] `R_B = 3wL/8`. Quadratic in `L`,
/// linear and sign-preserving in `w`, and — being pure statics — independent of `E` and
/// `I`. Returns `0` for non-physical input (`w` non-finite, or `L` non-positive or
/// non-finite).
pub fn propped_cantilever_udl_max_sagging_moment(load_per_length: f64, length: f64) -> f64 {
    if !load_per_length.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    9.0 * load_per_length * length * length / 128.0
}

/// The reaction at the **propped (simple) support of a propped cantilever under a
/// central point load** `R_B = 5P/16` (N) — the redundant support reaction for a beam
/// of span `length` `L` (m) that is built-in (clamped) at one end and simply supported
/// at the other, carrying a transverse point force `point_load` `P` (N) at mid-span
/// (`a = L/2`).
///
/// It is the point-load companion to the distributed-load
/// [`propped_cantilever_udl_prop_reaction`] `R_B = 3wL/8`. The beam is statically
/// indeterminate to the first degree; releasing the prop leaves a cantilever, and the
/// compatibility condition that the prop deflection vanish (`δ_B = 0`) gives the general
/// redundant reaction `R_B = P·a²(3L − a)/(2L³)`, which at the mid-span case `a = L/2`
/// collapses to the clean rational `R_B = 5P/16`. The fixed end then carries the
/// complement `R_A = 11P/16` and a clamping moment `M_A = 3PL/16`. Being pure statics
/// once the redundant is found, it is independent of `E` and `I`; it is linear and
/// sign-preserving in `P`, and — for the central case — independent of `L` (the `length`
/// argument is validated only by the physicality guard). Returns `0` for non-physical
/// input (`P` non-finite, or `L` non-positive or non-finite).
pub fn propped_cantilever_central_load_prop_reaction(point_load: f64, length: f64) -> f64 {
    if !point_load.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    5.0 * point_load / 16.0
}

/// The reaction at the **clamped (fixed) end of a propped cantilever under a central
/// point load** `R_A = 11P/16` (N) — the support reaction at the built-in end of a beam
/// of span `length` `L` (m) that is clamped at one end and simply supported at the
/// other, carrying a transverse point force `point_load` `P` (N) at mid-span.
///
/// It is the complement of the prop reaction
/// [`propped_cantilever_central_load_prop_reaction`] `R_B = 5P/16`: vertical equilibrium
/// of the whole beam requires `R_A + R_B = P`, so `R_A = P − 5P/16 = 11P/16`. The fixed
/// end carries the larger share of the load (and additionally resists the clamping
/// moment `M_A = 3PL/16`). Being pure statics it is independent of `E` and `I`, linear
/// and sign-preserving in `P`, and — for the central case — independent of `L` (the
/// `length` argument is validated only by the physicality guard). Returns `0` for
/// non-physical input (`P` non-finite, or `L` non-positive or non-finite).
pub fn propped_cantilever_central_load_fixed_end_reaction(point_load: f64, length: f64) -> f64 {
    if !point_load.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    11.0 * point_load / 16.0
}

/// The **clamping moment at the fixed end of a propped cantilever under a central
/// point load** `M_A = 3PL/16` (N·m) — the bending moment the built-in support must
/// resist on a beam of span `length` `L` (m) clamped at one end and simply supported at
/// the other, carrying a transverse point force `point_load` `P` (N) at mid-span. It is
/// the strength-governing quantity of the case (the largest bending moment along the
/// beam).
///
/// It completes the propped-cantilever central-load case alongside the reaction pair
/// [`propped_cantilever_central_load_prop_reaction`] `R_B = 5P/16` and
/// [`propped_cantilever_central_load_fixed_end_reaction`] `R_A = 11P/16`. By moment
/// equilibrium of the whole beam about the fixed end,
/// `M_A = P·(L/2) − R_B·L = PL/2 − 5PL/16 = 3PL/16`. Unlike the reactions it is
/// **linear in `L`** (a moment, not a force). Being pure statics it is independent of
/// `E` and `I`, and sign-preserving in `P`. Returns `0` for non-physical input (`P`
/// non-finite, or `L` non-positive or non-finite).
pub fn propped_cantilever_central_load_fixed_end_moment(point_load: f64, length: f64) -> f64 {
    if !point_load.is_finite() || !length.is_finite() || length <= 0.0 {
        return 0.0;
    }
    3.0 * point_load * length / 16.0
}

/// The analytic **strain energy of a simply-supported beam under a central point
/// load** `U = P²·L³/(96·E·I)` (J) — the elastic energy stored in bending when a
/// slender Euler–Bernoulli beam of span `length` `L` (m), Young's modulus
/// `youngs_modulus` `E` (Pa) and section second moment of area `second_moment_area`
/// `I` (m⁴) carries a transverse point force `load` `P` (N) at mid-span.
///
/// It is the simply-supported companion to the cantilever
/// [`cantilever_point_load_strain_energy`] in the energy-method family, and like it
/// follows from **Clapeyron's theorem** `U = ½·P·δ` with the mid-span deflection
/// [`simply_supported_center_deflection`] `δ = P·L³/(48EI)`. Comparing the two load
/// cases, a cantilever stores `16×` the energy of a simply-supported beam under the
/// same point load (`P²L³/(6EI)` vs `P²L³/(96EI)`) — the far stiffer propped-span
/// support condition. The energy grows with the *square* of the load (so it is
/// sign-independent), the *cube* of the span, and falls inversely with the flexural
/// rigidity `E·I`. Returns `0` for non-physical input (`P` non-finite, or `E`, `I`,
/// or `L` non-positive or non-finite).
pub fn simply_supported_point_load_strain_energy(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * load * length.powi(3) / (96.0 * youngs_modulus * second_moment_area)
}

/// The analytic **simply-supported (pinned–pinned) end slope** `θ = P·L²/(16·E·I)`
/// (rad) — the rotation of the elastic curve at *each support* of a slender
/// Euler–Bernoulli beam under a *central* transverse point load `load` `P` (N),
/// span `length` `L` (m), Young's modulus `youngs_modulus` `E` (Pa), section second
/// moment of area `second_moment_area` `I` (m⁴).
///
/// This is the *rotational* companion to the mid-span deflection
/// [`simply_supported_center_deflection`], and the simply-supported counterpart to
/// the clamped–free [`cantilever_tip_slope`]. The two simply-supported point-load
/// quantities are locked together by `δ_centre = (L/3)·θ_end`: the central
/// deflection is a third of the span times the support rotation. Like the
/// deflection it grows linearly with the load `P` (sign-preserving), with the
/// *square* of the span `L`, and falls inversely with the flexural rigidity `E·I`.
/// Returns `0` for non-physical input (`P` non-finite, or `E`, `I`, or `L`
/// non-positive or non-finite).
pub fn simply_supported_end_slope(
    load: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load * length.powi(2) / (16.0 * youngs_modulus * second_moment_area)
}

/// The analytic **simply-supported (pinned–pinned) mid-span deflection under a
/// uniformly distributed load** `δ = 5·w·L⁴/(384·E·I)` (m) of a slender
/// Euler–Bernoulli beam — span `length` `L` (m) carrying a transverse load of
/// intensity `load_per_length` `w` (N/m) spread evenly over the whole span, with
/// Young's modulus `youngs_modulus` `E` (Pa) and section second moment of area
/// `second_moment_area` `I` (m⁴).
///
/// This is the **distributed-load companion** to the point-load
/// [`simply_supported_center_deflection`] (`P·L³/(48·E·I)`), completing the
/// `{cantilever, simply-supported} × {point, UDL}` set of beam-deflection
/// references. The *same total load* `W = w·L`, spread uniformly instead of
/// concentrated at mid-span, deflects the centre to only `5/8` as far (`5/384`
/// versus `1/48`). The deflection grows linearly with the intensity `w` (and is
/// sign-preserving — an upward load lifts mid-span), with the *fourth* power of
/// the span `L`, and falls inversely with the flexural rigidity `E·I`. Returns
/// `0` for non-physical input (`w` non-finite, or `E`, `I`, or `L` non-positive
/// or non-finite).
pub fn simply_supported_udl_center_deflection(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    5.0 * load_per_length * length.powi(4) / (384.0 * youngs_modulus * second_moment_area)
}

/// The analytic **simply-supported (pinned–pinned) end slope under a uniformly
/// distributed load** `θ = w·L³/(24·E·I)` (rad) — the support rotation of a slender
/// Euler–Bernoulli beam carrying a transverse load of intensity `load_per_length`
/// `w` (N/m) over its full span `length` `L` (m), with Young's modulus
/// `youngs_modulus` `E` (Pa) and section second moment of area `second_moment_area`
/// `I` (m⁴).
///
/// This is the *distributed-load* slope companion to the point-load
/// [`simply_supported_end_slope`] (`P·L²/(16EI)`) and the rotational companion to
/// the UDL mid-span deflection [`simply_supported_udl_center_deflection`]
/// (`5wL⁴/(384EI)`); the two UDL quantities are locked by `δ_centre = (5/16)·L·θ`.
/// With it the analytic slope set spans both support conditions and both load
/// types. The slope grows linearly with `w` (sign-preserving), with the *cube* of
/// the span `L`, and falls inversely with the flexural rigidity `E·I`. Returns `0`
/// for non-physical input (`w` non-finite, or `E`, `I`, or `L` non-positive or
/// non-finite).
pub fn simply_supported_udl_end_slope(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    second_moment_area: f64,
) -> f64 {
    if !load_per_length.is_finite()
        || !length.is_finite()
        || length <= 0.0
        || !youngs_modulus.is_finite()
        || youngs_modulus <= 0.0
        || !second_moment_area.is_finite()
        || second_moment_area <= 0.0
    {
        return 0.0;
    }
    load_per_length * length.powi(3) / (24.0 * youngs_modulus * second_moment_area)
}

/// Errors from the native 3D beam solver.
#[derive(Debug, Error)]
pub enum BeamSolverError {
    /// The frame has no nodes.
    #[error("beam model has no nodes")]
    EmptyModel,
    /// The frame has no elements.
    #[error("beam model has no elements")]
    NoElements,
    /// A beam element has zero (or near-zero) length — its two end
    /// nodes coincide.
    #[error("beam element {0} has zero length (coincident end nodes)")]
    ZeroLength(usize),
    /// An element references a node index past the end of the node
    /// array.
    #[error("beam element {elem} references node {node} but the model has only {n_nodes} nodes")]
    BadConnectivity {
        /// 0-based element index.
        elem: usize,
        /// Out-of-range node index.
        node: usize,
        /// Node-array length.
        n_nodes: usize,
    },
    /// A cross-section property is non-physical (≤ 0).
    #[error("invalid cross-section for element {elem}: {what} must be positive")]
    BadSection {
        /// 0-based element index.
        elem: usize,
        /// Which property was bad.
        what: &'static str,
    },
    /// A material constant is non-physical.
    #[error("invalid material: {0}")]
    BadMaterial(String),
    /// No constraint was supplied, so the structure can float — the
    /// stiffness matrix is rigid-body singular.
    #[error("no constraint — the frame is unrestrained (rigid-body singular)")]
    Unconstrained,
    /// The linear solve failed: the assembled system was not positive-
    /// definite (an under-constrained or mechanism frame).
    #[error("linear solve failed: stiffness matrix is not positive-definite")]
    SolveFailed,
    /// A modal solve was asked for more modes than the constrained
    /// system can supply, or for zero modes.
    #[error("requested {requested} modes but the constrained system has only {available} DOFs")]
    TooManyModes {
        /// Modes requested.
        requested: usize,
        /// Free DOFs available.
        available: usize,
    },
    /// Every DOF is constrained — nothing left to vibrate / deflect.
    #[error("all degrees of freedom are constrained")]
    FullyConstrained,
    /// The reduced mass matrix was not positive-definite.
    #[error("reduced mass matrix is not positive-definite")]
    MassNotPositiveDefinite,
    /// The symmetric eigensolver did not converge.
    #[error("symmetric eigensolver failed to converge")]
    EigenFailed,
    /// The frame is too large for the **dense** beam assembly path. Each
    /// beam node carries six DOFs, so the global stiffness (and, for a
    /// modal solve, mass) is a dense `6·n_nodes × 6·n_nodes` `f64`
    /// matrix needing `8·(6·n_nodes)²` bytes. Round-2 fix: the continuum
    /// solver's [`crate::native_solver::MAX_DENSE_DOFS`] cap covered only
    /// the 3-DOF/node volumetric path (via `3·n_nodes`), so the beam
    /// path's `6·n_nodes` allocation was uncapped — a large frame would
    /// OOM the host. The DOF count is now routed through
    /// [`crate::native_solver::check_dense_dof_count`] *before* any
    /// [`nalgebra::DMatrix::zeros`] allocation.
    #[error(
        "dense beam solve needs {dofs} DOFs, the dense path supports at most {max}; \
         coarsen the frame"
    )]
    TooLarge {
        /// The frame's degree-of-freedom count `6·n_nodes` (or the
        /// overflow sentinel `usize::MAX` if `6·n_nodes` itself
        /// overflowed).
        dofs: usize,
        /// The dense-path upper bound,
        /// [`crate::native_solver::MAX_DENSE_DOFS`].
        max: usize,
    },
    /// A beam load or boundary-condition input carried a non-finite
    /// value (`NaN` or `±∞`). Round-1 fix: beam forces / moments /
    /// prescribed displacements were pushed straight into the RHS, where
    /// the Cholesky back-substitution turns a non-finite input into a
    /// silently-non-finite displacement returned as `Ok(..)`. Validating
    /// up front lets the error name the cause.
    #[error("non-finite {kind} at node {node} (NaN or infinity is not a valid load/BC)")]
    InvalidLoad {
        /// 0-based node index carrying the bad value.
        node: usize,
        /// Which input was non-finite — `"force"`, `"moment"`, or
        /// `"prescribed displacement"`.
        kind: &'static str,
    },
}

/// Cross-section properties of a prismatic beam element.
///
/// `iy` / `iz` are the second moments of area about the section's two
/// **principal** axes (local `y` and `z`); `j` is the
/// **torsion constant** (the St-Venant `J`, equal to the polar moment
/// `Iy+Iz` only for a circular section). `shear_*` are the Timoshenko
/// **shear correction factors** (`κ ≈ 5/6` for a rectangle, `≈ 0.9` for
/// a solid circle, `1.0` disables transverse-shear flexibility and
/// recovers the Euler-Bernoulli element).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BeamSection {
    /// Cross-section area `A` in m².
    pub area: f64,
    /// Second moment of area about the local `y` axis, `Iy` in m⁴
    /// (governs bending in the local `x-z` plane).
    pub iy: f64,
    /// Second moment of area about the local `z` axis, `Iz` in m⁴
    /// (governs bending in the local `x-y` plane).
    pub iz: f64,
    /// St-Venant torsion constant `J` in m⁴.
    pub j: f64,
    /// Timoshenko shear correction factor for shear along local `y`
    /// (dimensionless, `(0,1]`). `1.0` → no shear flexibility.
    pub shear_y: f64,
    /// Timoshenko shear correction factor for shear along local `z`.
    pub shear_z: f64,
}

impl BeamSection {
    /// A solid **rectangular** section, `width` along local `y`,
    /// `height` along local `z`.
    ///
    /// `Iz = w·h³/12`, `Iy = h·w³/12`, the torsion constant uses the
    /// standard thin-to-square rectangle approximation, and the shear
    /// factor is `5/6`.
    pub fn rectangle(width: f64, height: f64) -> Self {
        let a = width * height;
        // Note local convention: Iz governs x-y plane bending (about z)
        //   → depends on the y-extent `width`  → Iz = h·w³/12.
        //   Iy governs x-z plane bending (about y)
        //   → depends on the z-extent `height` → Iy = w·h³/12.
        let iz = height * width.powi(3) / 12.0;
        let iy = width * height.powi(3) / 12.0;
        // St-Venant torsion constant of a rectangle (Roark): with
        // a = ½·long side, b = ½·short side,
        //   J = a·b³·(16/3 − 3.36·(b/a)·(1 − b⁴/(12a⁴))).
        let (long, short) = if width >= height {
            (width, height)
        } else {
            (height, width)
        };
        let a_h = long / 2.0;
        let b_h = short / 2.0;
        let ratio = b_h / a_h;
        let j = a_h
            * b_h.powi(3)
            * (16.0 / 3.0 - 3.36 * ratio * (1.0 - ratio.powi(4) / 12.0));
        Self {
            area: a,
            iy,
            iz,
            j,
            shear_y: 5.0 / 6.0,
            shear_z: 5.0 / 6.0,
        }
    }

    /// A solid **circular** section of the given `radius`.
    ///
    /// `Iy = Iz = πr⁴/4`, the torsion constant is the polar moment
    /// `J = πr⁴/2`, and the shear factor is `0.9`.
    pub fn circle(radius: f64) -> Self {
        let a = std::f64::consts::PI * radius * radius;
        let i = std::f64::consts::PI * radius.powi(4) / 4.0;
        Self {
            area: a,
            iy: i,
            iz: i,
            j: 2.0 * i,
            shear_y: 0.9,
            shear_z: 0.9,
        }
    }

    /// A thin-walled **circular tube**, `outer` / `inner` radii.
    pub fn tube(outer: f64, inner: f64) -> Self {
        let a = std::f64::consts::PI * (outer * outer - inner * inner);
        let i = std::f64::consts::PI * (outer.powi(4) - inner.powi(4)) / 4.0;
        Self {
            area: a,
            iy: i,
            iz: i,
            j: 2.0 * i,
            shear_y: 0.5,
            shear_z: 0.5,
        }
    }

    /// Validate that every property is finite and positive.
    fn check(&self, elem: usize) -> Result<(), BeamSolverError> {
        let bad = |what| BeamSolverError::BadSection { elem, what };
        if !(self.area.is_finite()) || self.area <= 0.0 {
            return Err(bad("area"));
        }
        if !(self.iy.is_finite()) || self.iy <= 0.0 {
            return Err(bad("Iy"));
        }
        if !(self.iz.is_finite()) || self.iz <= 0.0 {
            return Err(bad("Iz"));
        }
        if !(self.j.is_finite()) || self.j <= 0.0 {
            return Err(bad("J"));
        }
        if !(self.shear_y.is_finite()) || self.shear_y <= 0.0 {
            return Err(bad("shear_y"));
        }
        if !(self.shear_z.is_finite()) || self.shear_z <= 0.0 {
            return Err(bad("shear_z"));
        }
        Ok(())
    }
}

/// One 2-node 3D beam element.
///
/// `nodes` are the two end-node indices. `section` carries the
/// cross-section properties. `orientation` is an optional reference
/// vector that, together with the beam axis, fixes the cross-section's
/// principal-`y` direction (so a beam can be rolled about its axis);
/// `None` uses an automatic, well-conditioned choice.
#[derive(Copy, Clone, Debug)]
pub struct BeamElement {
    /// The two end-node indices `[start, end]`.
    pub nodes: [usize; 2],
    /// Cross-section properties.
    pub section: BeamSection,
    /// Optional roll reference: a vector that is *not* parallel to the
    /// beam axis; the local `y` axis is taken in the plane it spans
    /// with the axis. `None` → an automatic choice.
    pub orientation: Option<Vector3<f64>>,
}

impl BeamElement {
    /// A beam element between two nodes with the given section and the
    /// automatic cross-section orientation.
    pub fn new(start: usize, end: usize, section: BeamSection) -> Self {
        Self {
            nodes: [start, end],
            section,
            orientation: None,
        }
    }

    /// A beam element with an explicit roll-reference vector.
    pub fn with_orientation(
        start: usize,
        end: usize,
        section: BeamSection,
        orientation: Vector3<f64>,
    ) -> Self {
        Self {
            nodes: [start, end],
            section,
            orientation: Some(orientation),
        }
    }
}

/// A single-node six-DOF constraint for the beam solver.
///
/// `fixed[0..3]` pin the three translations, `fixed[3..6]` the three
/// rotations; `Some(v)` fixes that DOF to value `v`, `None` leaves it
/// free.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BeamConstraint {
    /// 0-based node index.
    pub node: usize,
    /// Per-DOF pin: `[ux,uy,uz, θx,θy,θz]`.
    pub fixed: [Option<f64>; 6],
}

impl BeamConstraint {
    /// Fully clamp a node — all three translations and all three
    /// rotations fixed to zero (an encastré / built-in support).
    pub fn clamped(node: usize) -> Self {
        Self {
            node,
            fixed: [Some(0.0); 6],
        }
    }

    /// A **pinned** support — translations fixed, rotations free
    /// (a frictionless spherical joint / simple support).
    pub fn pinned(node: usize) -> Self {
        Self {
            node,
            fixed: [
                Some(0.0),
                Some(0.0),
                Some(0.0),
                None,
                None,
                None,
            ],
        }
    }
}

/// A single-node six-component load for the beam solver — a
/// concentrated force and/or moment.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BeamLoad {
    /// 0-based node index.
    pub node: usize,
    /// Force `[Fx,Fy,Fz]` in newtons.
    pub force: [f64; 3],
    /// Moment `[Mx,My,Mz]` in newton-metres.
    pub moment: [f64; 3],
}

impl BeamLoad {
    /// A pure force at a node (zero moment).
    pub fn force(node: usize, force: [f64; 3]) -> Self {
        Self {
            node,
            force,
            moment: [0.0; 3],
        }
    }

    /// A pure moment at a node (zero force).
    pub fn moment(node: usize, moment: [f64; 3]) -> Self {
        Self {
            node,
            force: [0.0; 3],
            moment,
        }
    }
}

/// Result of a 3D beam static solve.
#[derive(Clone, Debug)]
pub struct BeamSolution {
    /// Per-node translation `[ux,uy,uz]` in metres.
    pub translation: Vec<[f64; 3]>,
    /// Per-node rotation `[θx,θy,θz]` in radians.
    pub rotation: Vec<[f64; 3]>,
}

impl BeamSolution {
    /// Largest nodal translation magnitude.
    pub fn max_translation(&self) -> f64 {
        self.translation
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }

    /// Largest nodal rotation magnitude.
    pub fn max_rotation(&self) -> f64 {
        self.rotation
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }
}

/// One natural mode of a beam frame.
#[derive(Clone, Debug)]
pub struct BeamMode {
    /// Natural frequency in hertz.
    pub frequency_hz: f64,
    /// Angular frequency in rad/s.
    pub angular_frequency: f64,
    /// Per-node mode shape — `(translation, rotation)`, mass-normalised.
    pub translation: Vec<[f64; 3]>,
    /// Per-node rotational mode shape.
    pub rotation: Vec<[f64; 3]>,
}

/// Result of a 3D beam modal solve.
#[derive(Clone, Debug)]
pub struct BeamModalSolution {
    /// The modes, ascending in frequency.
    pub modes: Vec<BeamMode>,
}

impl BeamModalSolution {
    /// Fundamental (lowest) frequency in hertz, or `None` if empty.
    pub fn fundamental_hz(&self) -> Option<f64> {
        self.modes.first().map(|m| m.frequency_hz)
    }
}

/// Shear modulus `G = E / (2(1+ν))` of an isotropic material.
fn shear_modulus(m: &FemMaterial) -> Result<f64, BeamSolverError> {
    let e = m.youngs_modulus;
    let nu = m.poisson_ratio;
    if !e.is_finite() || e <= 0.0 {
        return Err(BeamSolverError::BadMaterial(format!(
            "Young's modulus must be finite and positive, got {e}"
        )));
    }
    if !nu.is_finite() || nu <= -1.0 || nu >= 0.5 {
        return Err(BeamSolverError::BadMaterial(format!(
            "Poisson's ratio must lie in (-1, 0.5), got {nu}"
        )));
    }
    Ok(e / (2.0 * (1.0 + nu)))
}

/// Build the 3×3 local→global rotation `R` of a beam element.
///
/// Columns of `R` are the local `x` (along the beam), local `y` and
/// local `z` axes expressed in global coordinates. `local_x` is the
/// unit beam axis; `local_y` is taken from the roll reference (or an
/// automatic choice) made orthogonal to `x`; `local_z = x × y`.
fn beam_triad(
    p_start: Vector3<f64>,
    p_end: Vector3<f64>,
    orientation: Option<Vector3<f64>>,
) -> Option<Matrix3<f64>> {
    let axis = p_end - p_start;
    let len = axis.norm();
    if len < 1.0e-12 {
        return None;
    }
    let local_x = axis / len;
    // Reference vector for the cross-section roll.
    let reference = match orientation {
        Some(v) if v.norm() > 1.0e-12 => v.normalize(),
        _ => {
            // Automatic: use global Z unless the beam is (near-)vertical,
            // in which case use global Y. This is the standard
            // well-conditioned default.
            if local_x.z.abs() < 0.9 {
                Vector3::new(0.0, 0.0, 1.0)
            } else {
                Vector3::new(0.0, 1.0, 0.0)
            }
        }
    };
    // local_y ⟂ local_x, in the plane of (x, reference).
    let mut local_y = reference - local_x * local_x.dot(&reference);
    let yn = local_y.norm();
    if yn < 1.0e-9 {
        // reference parallel to the axis — fall back.
        let alt = if local_x.x.abs() < 0.9 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };
        local_y = alt - local_x * local_x.dot(&alt);
    }
    local_y.normalize_mut();
    let local_z = local_x.cross(&local_y);
    Some(Matrix3::from_columns(&[local_x, local_y, local_z]))
}

/// The 12×12 **local** stiffness matrix of a 2-node Timoshenko beam.
///
/// DOF order is `[u1 v1 w1 θx1 θy1 θz1  u2 v2 w2 θx2 θy2 θz2]` in the
/// element's local frame (`x` along the beam). Couples axial, torsion
/// and the two bending planes; the bending blocks carry the Timoshenko
/// shear-flexibility factors `Φ` so a stocky beam is handled correctly.
fn beam_local_stiffness(
    length: f64,
    e: f64,
    g: f64,
    s: &BeamSection,
) -> DMatrix<f64> {
    let l = length;
    let mut k = DMatrix::<f64>::zeros(12, 12);

    // --- axial: DOFs 0 and 6 ---
    let ea_l = e * s.area / l;
    k[(0, 0)] += ea_l;
    k[(6, 6)] += ea_l;
    k[(0, 6)] -= ea_l;
    k[(6, 0)] -= ea_l;

    // --- torsion: DOFs 3 and 9 ---
    let gj_l = g * s.j / l;
    k[(3, 3)] += gj_l;
    k[(9, 9)] += gj_l;
    k[(3, 9)] -= gj_l;
    k[(9, 3)] -= gj_l;

    // --- bending in the local x-y plane (about local z) ---
    // Transverse v (DOFs 1,7), rotation θz (DOFs 5,11). The Timoshenko
    // shear parameter Φy = 12·E·Iz / (κy·G·A·L²).
    let phi_y = 12.0 * e * s.iz / (s.shear_y * g * s.area * l * l);
    add_bending_block(&mut k, l, e * s.iz, phi_y, 1, 5, 7, 11, 1.0);

    // --- bending in the local x-z plane (about local y) ---
    // Transverse w (DOFs 2,8), rotation θy (DOFs 4,10). The coupling of
    // a +w with a rotation θy has the opposite sign to the x-y plane
    // (right-handed frame), captured by the `sign = -1` argument.
    let phi_z = 12.0 * e * s.iy / (s.shear_z * g * s.area * l * l);
    add_bending_block(&mut k, l, e * s.iy, phi_z, 2, 4, 8, 10, -1.0);

    k
}

/// Scatter the 4×4 Timoshenko bending sub-matrix into the 12×12 local
/// stiffness for one bending plane.
///
/// `tr_a` / `tr_b` are the transverse-translation DOF indices of the
/// two nodes; `rot_a` / `rot_b` the matching bending-rotation DOF
/// indices. `ei` is the bending rigidity, `phi` the Timoshenko shear
/// parameter, `sign` the ±1 that orients the translation-rotation
/// coupling for the plane.
#[allow(clippy::too_many_arguments)]
fn add_bending_block(
    k: &mut DMatrix<f64>,
    l: f64,
    ei: f64,
    phi: f64,
    tr_a: usize,
    rot_a: usize,
    tr_b: usize,
    rot_b: usize,
    sign: f64,
) {
    // Classical Timoshenko 4×4 bending stiffness. With Φ the shear
    // parameter the common factor is EI / (L³(1+Φ)).
    let f = ei / (l * l * l * (1.0 + phi));
    let k_tt = 12.0 * f; // translation-translation
    let k_tr = 6.0 * l * f * sign; // translation-rotation
    let k_rr_near = (4.0 + phi) * l * l * f; // same-node rotation-rotation
    let k_rr_far = (2.0 - phi) * l * l * f; // cross-node rotation-rotation

    // The 4 DOFs of this plane, in order [tr_a, rot_a, tr_b, rot_b].
    let idx = [tr_a, rot_a, tr_b, rot_b];
    // The 4×4 sub-stiffness (standard Timoshenko form).
    let sub = [
        [k_tt, k_tr, -k_tt, k_tr],
        [k_tr, k_rr_near, -k_tr, k_rr_far],
        [-k_tt, -k_tr, k_tt, -k_tr],
        [k_tr, k_rr_far, -k_tr, k_rr_near],
    ];
    for (a, &ia) in idx.iter().enumerate() {
        for (b, &ib) in idx.iter().enumerate() {
            k[(ia, ib)] += sub[a][b];
        }
    }
}

/// The 12×12 **local consistent mass** matrix of a 2-node beam.
///
/// Translational inertia uses the classical cubic-Hermite consistent
/// mass; rotary (torsional + bending-rotation) inertia uses the
/// standard lumped-with-coupling form. This is the textbook beam
/// consistent mass — exact for the element kinematics and accurate for
/// modal analysis. `rho` is the density.
fn beam_local_mass(length: f64, rho: f64, s: &BeamSection) -> DMatrix<f64> {
    let l = length;
    let mut m = DMatrix::<f64>::zeros(12, 12);
    let mass = rho * s.area * l; // total element mass

    // --- axial inertia (DOFs 0, 6): consistent 2×2 [2 1;1 2]·(m/6) ---
    let a6 = mass / 6.0;
    m[(0, 0)] += 2.0 * a6;
    m[(6, 6)] += 2.0 * a6;
    m[(0, 6)] += a6;
    m[(6, 0)] += a6;

    // --- torsional inertia (DOFs 3, 9) ---
    // Polar mass moment per length = ρ·(Iy+Iz); consistent [2 1;1 2]/6.
    let polar = rho * (s.iy + s.iz) * l;
    let t6 = polar / 6.0;
    m[(3, 3)] += 2.0 * t6;
    m[(9, 9)] += 2.0 * t6;
    m[(3, 9)] += t6;
    m[(9, 3)] += t6;

    // --- bending-plane consistent mass (translational part) ---
    // The classical cubic-Hermite consistent mass for bending; the
    // dominant term for a slender beam. Rotary-inertia contributions
    // are small and added as a lumped term for robustness.
    add_bending_mass(&mut m, l, mass, 1, 5, 7, 11, 1.0);
    add_bending_mass(&mut m, l, mass, 2, 4, 8, 10, -1.0);

    m
}

/// Scatter the 4×4 cubic-Hermite consistent-mass sub-matrix for one
/// bending plane into the 12×12 local mass.
#[allow(clippy::too_many_arguments)]
fn add_bending_mass(
    m: &mut DMatrix<f64>,
    l: f64,
    mass: f64,
    tr_a: usize,
    rot_a: usize,
    tr_b: usize,
    rot_b: usize,
    sign: f64,
) {
    // Classical Euler-Bernoulli cubic-Hermite consistent mass,
    // common factor m/420.
    let f = mass / 420.0;
    let idx = [tr_a, rot_a, tr_b, rot_b];
    // 4×4 sub-mass. The off-diagonal translation-rotation terms carry
    // the plane sign so the matrix matches the stiffness convention.
    let s = sign;
    let sub = [
        [156.0 * f, 22.0 * l * f * s, 54.0 * f, -13.0 * l * f * s],
        [
            22.0 * l * f * s,
            4.0 * l * l * f,
            13.0 * l * f * s,
            -3.0 * l * l * f,
        ],
        [54.0 * f, 13.0 * l * f * s, 156.0 * f, -22.0 * l * f * s],
        [
            -13.0 * l * f * s,
            -3.0 * l * l * f,
            -22.0 * l * f * s,
            4.0 * l * l * f,
        ],
    ];
    for (a, &ia) in idx.iter().enumerate() {
        for (b, &ib) in idx.iter().enumerate() {
            m[(ia, ib)] += sub[a][b];
        }
    }
}

/// The 12×12 block-diagonal transform `T` that rotates a beam
/// element's local matrices to global coordinates.
///
/// `T` is four copies of the 3×3 triad `R` (one per node, one each for
/// the translation and rotation triplets). A local matrix `Kₗ` becomes
/// the global `Kg = Tᵀ·Kₗ·T`.
fn beam_transform(r: &Matrix3<f64>) -> DMatrix<f64> {
    let mut t = DMatrix::<f64>::zeros(12, 12);
    // R maps local→global; the transform that takes a global DOF vector
    // to local is Rᵀ applied blockwise. We build T = blkdiag(Rᵀ ×4) so
    // that  u_local = T · u_global  and  Kg = Tᵀ Kl T.
    let rt = r.transpose();
    for block in 0..4 {
        let o = 3 * block;
        for i in 0..3 {
            for j in 0..3 {
                t[(o + i, o + j)] = rt[(i, j)];
            }
        }
    }
    t
}

/// Reject a frame too large for the **dense** beam path *before* any
/// `6·n_nodes × 6·n_nodes` allocation, returning the validated DOF count.
///
/// Each beam node carries six DOFs, so the dense path's `6·n_nodes` DOF
/// count is twice the continuum solver's `3·n_nodes`. The 6× multiply is
/// guarded with [`usize::checked_mul`] (a wrap would understate the size
/// and defeat the cap), then the count is routed through the single
/// shared capping path
/// [`crate::native_solver::check_dense_dof_count`] so the beam path uses
/// the *same* [`crate::native_solver::MAX_DENSE_DOFS`] bound as the
/// volumetric path — the `NativeSolverError::TooLarge` it returns is
/// re-tagged as the beam-local [`BeamSolverError::TooLarge`]. Pure
/// `O(1)` arithmetic, no allocation.
fn check_beam_dense_dofs(n_nodes: usize) -> Result<usize, BeamSolverError> {
    let n_dof = n_nodes
        .checked_mul(6)
        .ok_or(BeamSolverError::TooLarge {
            dofs: usize::MAX,
            max: crate::native_solver::MAX_DENSE_DOFS,
        })?;
    crate::native_solver::check_dense_dof_count(n_dof).map_err(|e| match e {
        crate::native_solver::NativeSolverError::TooLarge { dofs, max } => {
            BeamSolverError::TooLarge { dofs, max }
        }
        // `check_dense_dof_count` only ever returns `TooLarge`.
        _ => BeamSolverError::TooLarge {
            dofs: n_dof,
            max: crate::native_solver::MAX_DENSE_DOFS,
        },
    })
}

/// Assemble the global beam stiffness (and optionally consistent mass)
/// as dense `6·n_nodes` matrices.
///
/// Returns `(K, optional M)`. When `with_mass` is false `M` is `None`
/// — the static solver does not need it.
#[allow(clippy::type_complexity)]
fn assemble_beam_system(
    nodes: &[Vector3<f64>],
    elements: &[BeamElement],
    material: &FemMaterial,
    with_mass: bool,
) -> Result<(DMatrix<f64>, Option<DMatrix<f64>>), BeamSolverError> {
    let n_nodes = nodes.len();
    if n_nodes == 0 {
        return Err(BeamSolverError::EmptyModel);
    }
    if elements.is_empty() {
        return Err(BeamSolverError::NoElements);
    }
    let e = material.youngs_modulus;
    if !e.is_finite() || e <= 0.0 {
        return Err(BeamSolverError::BadMaterial(format!(
            "Young's modulus must be finite and positive, got {e}"
        )));
    }
    let g = shear_modulus(material)?;
    if with_mass && (!material.density.is_finite() || material.density <= 0.0) {
        return Err(BeamSolverError::BadMaterial(format!(
            "density must be finite and positive for a modal solve, got {}",
            material.density
        )));
    }

    // Cap the dense `6·n_nodes × 6·n_nodes` allocation BEFORE it is made.
    // This is the single chokepoint both `solve_beam_static` and
    // `solve_beam_modal` flow through, so capping here covers both their
    // dense paths (the modal solve's later `n_free² ≤ n_dof²` reduction
    // is bounded by the same check).
    let n_dof = check_beam_dense_dofs(n_nodes)?;
    let mut k = DMatrix::<f64>::zeros(n_dof, n_dof);
    let mut m = if with_mass {
        Some(DMatrix::<f64>::zeros(n_dof, n_dof))
    } else {
        None
    };

    for (ei, elem) in elements.iter().enumerate() {
        for &nd in &elem.nodes {
            if nd >= n_nodes {
                return Err(BeamSolverError::BadConnectivity {
                    elem: ei,
                    node: nd,
                    n_nodes,
                });
            }
        }
        elem.section.check(ei)?;
        let p0 = nodes[elem.nodes[0]];
        let p1 = nodes[elem.nodes[1]];
        let length = (p1 - p0).norm();
        if length < 1.0e-12 {
            return Err(BeamSolverError::ZeroLength(ei));
        }
        let r = beam_triad(p0, p1, elem.orientation).ok_or(BeamSolverError::ZeroLength(ei))?;
        let t = beam_transform(&r);
        let t_t = t.transpose();

        let kl = beam_local_stiffness(length, e, g, &elem.section);
        let kg = &t_t * &kl * &t;
        scatter_beam(&mut k, &kg, &elem.nodes);

        if let Some(ref mut mm) = m {
            let ml = beam_local_mass(length, material.density, &elem.section);
            let mg = &t_t * &ml * &t;
            scatter_beam(mm, &mg, &elem.nodes);
        }
    }
    Ok((k, m))
}

/// Scatter a 12×12 element matrix into the global `6·n_nodes` system.
/// Local DOF `6a+i` of element node `a` → global DOF `6·node[a]+i`.
fn scatter_beam(global: &mut DMatrix<f64>, elem: &DMatrix<f64>, nodes: &[usize; 2]) {
    for a in 0..2 {
        for i in 0..6 {
            let gi = 6 * nodes[a] + i;
            for b in 0..2 {
                for j in 0..6 {
                    let gj = 6 * nodes[b] + j;
                    global[(gi, gj)] += elem[(6 * a + i, 6 * b + j)];
                }
            }
        }
    }
}

/// Solve a **linear-static** 3D beam-frame problem.
///
/// `nodes` are the frame's node coordinates; `elements` the beam
/// members with their cross-sections; `material` the (isotropic)
/// elastic constants; `constraints` the six-DOF supports (at least one
/// is required, or the frame floats); `loads` the concentrated nodal
/// forces and moments.
///
/// Returns the per-node translation and rotation fields.
///
/// # Method
///
/// Each element's 12×12 local Timoshenko stiffness is rotated to global
/// coordinates by the block-diagonal triad transform and scatter-added
/// into the global `6·n_nodes` system. Forces and moments go straight
/// into the load vector; constraints are imposed by the large-penalty
/// method (keeps the matrix SPD so the Cholesky path stays valid); the
/// system is factorised with [`CscCholesky`].
///
/// # Errors
///
/// See [`BeamSolverError`].
pub fn solve_beam_static(
    nodes: &[Vector3<f64>],
    elements: &[BeamElement],
    material: &FemMaterial,
    constraints: &[BeamConstraint],
    loads: &[BeamLoad],
) -> Result<BeamSolution, BeamSolverError> {
    let (k, _m) = assemble_beam_system(nodes, elements, material, false)?;
    let n_nodes = nodes.len();
    let n_dof = 6 * n_nodes;

    if constraints.is_empty() {
        return Err(BeamSolverError::Unconstrained);
    }

    // Peak diagonal → penalty scale.
    let mut max_diag = 0.0_f64;
    for i in 0..n_dof {
        max_diag = max_diag.max(k[(i, i)].abs());
    }
    if max_diag <= 0.0 {
        return Err(BeamSolverError::SolveFailed);
    }
    let penalty = max_diag * 1.0e8;

    // Load vector.
    let mut f = DVector::<f64>::zeros(n_dof);
    for load in loads {
        if load.node >= n_nodes {
            return Err(BeamSolverError::BadConnectivity {
                elem: usize::MAX,
                node: load.node,
                n_nodes,
            });
        }
        // Round-1 H4: reject a non-finite force / moment before it
        // reaches the RHS — the Cholesky solve would otherwise return a
        // silently-non-finite displacement as Ok(..).
        if load.force.iter().any(|v| !v.is_finite()) {
            return Err(BeamSolverError::InvalidLoad {
                node: load.node,
                kind: "force",
            });
        }
        if load.moment.iter().any(|v| !v.is_finite()) {
            return Err(BeamSolverError::InvalidLoad {
                node: load.node,
                kind: "moment",
            });
        }
        for i in 0..3 {
            f[6 * load.node + i] += load.force[i];
            f[6 * load.node + 3 + i] += load.moment[i];
        }
    }

    // Penalty constraints.
    let mut penalty_diag = vec![0.0_f64; n_dof];
    let mut any = false;
    for c in constraints {
        if c.node >= n_nodes {
            return Err(BeamSolverError::BadConnectivity {
                elem: usize::MAX,
                node: c.node,
                n_nodes,
            });
        }
        for (i, fixed) in c.fixed.iter().enumerate() {
            if let Some(value) = fixed {
                // Round-1 H4: a prescribed displacement / rotation is
                // folded into the RHS as `penalty·value`; reject a
                // non-finite value.
                if !value.is_finite() {
                    return Err(BeamSolverError::InvalidLoad {
                        node: c.node,
                        kind: "prescribed displacement",
                    });
                }
                let dof = 6 * c.node + i;
                penalty_diag[dof] += penalty;
                f[dof] += penalty * value;
                any = true;
            }
        }
    }
    if !any {
        return Err(BeamSolverError::Unconstrained);
    }

    // Build the sparse stiffened system and factorise.
    let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);
    for i in 0..n_dof {
        for j in 0..n_dof {
            let mut v = k[(i, j)];
            if i == j {
                v += penalty_diag[i];
            }
            if v != 0.0 {
                coo.push(i, j, v);
            }
        }
    }
    let csc = CscMatrix::from(&coo);
    let chol = CscCholesky::factor(&csc).map_err(|_| BeamSolverError::SolveFailed)?;
    let u = chol.solve(&f);
    let u = u.column(0);

    let mut translation = vec![[0.0_f64; 3]; n_nodes];
    let mut rotation = vec![[0.0_f64; 3]; n_nodes];
    for n in 0..n_nodes {
        for i in 0..3 {
            translation[n][i] = u[6 * n + i];
            rotation[n][i] = u[6 * n + 3 + i];
        }
    }
    Ok(BeamSolution {
        translation,
        rotation,
    })
}

/// Solve the **modal** (natural-frequency) eigenproblem of a 3D beam
/// frame.
///
/// Assembles the global beam stiffness `K` and consistent mass `M`,
/// eliminates the constrained DOFs, and solves the generalised
/// symmetric eigenproblem `K φ = λ M φ` for the lowest `n_modes` —
/// exactly the treatment [`crate::modal_solver`] uses for the
/// continuum solver. Returns the natural frequencies and mass-
/// normalised mode shapes.
///
/// # Errors
///
/// See [`BeamSolverError`].
pub fn solve_beam_modal(
    nodes: &[Vector3<f64>],
    elements: &[BeamElement],
    material: &FemMaterial,
    constraints: &[BeamConstraint],
    n_modes: usize,
) -> Result<BeamModalSolution, BeamSolverError> {
    if n_modes == 0 {
        return Err(BeamSolverError::TooManyModes {
            requested: 0,
            available: 0,
        });
    }
    let (k, m) = assemble_beam_system(nodes, elements, material, true)?;
    let m = m.expect("mass requested");
    let n_nodes = nodes.len();
    let n_dof = 6 * n_nodes;

    // Free-DOF set.
    let mut constrained = vec![false; n_dof];
    for c in constraints {
        if c.node >= n_nodes {
            return Err(BeamSolverError::BadConnectivity {
                elem: usize::MAX,
                node: c.node,
                n_nodes,
            });
        }
        for (i, fixed) in c.fixed.iter().enumerate() {
            if fixed.is_some() {
                constrained[6 * c.node + i] = true;
            }
        }
    }
    let free: Vec<usize> = (0..n_dof).filter(|&d| !constrained[d]).collect();
    if free.is_empty() {
        return Err(BeamSolverError::FullyConstrained);
    }
    let n_free = free.len();
    if n_modes > n_free {
        return Err(BeamSolverError::TooManyModes {
            requested: n_modes,
            available: n_free,
        });
    }

    // Reduce.
    let mut k_ff = DMatrix::<f64>::zeros(n_free, n_free);
    let mut m_ff = DMatrix::<f64>::zeros(n_free, n_free);
    for (ri, &gr) in free.iter().enumerate() {
        for (ci, &gc) in free.iter().enumerate() {
            k_ff[(ri, ci)] = k[(gr, gc)];
            m_ff[(ri, ci)] = m[(gr, gc)];
        }
    }

    // Generalised → standard via the Cholesky factor of M_ff.
    let chol = m_ff
        .clone()
        .cholesky()
        .ok_or(BeamSolverError::MassNotPositiveDefinite)?;
    let l = chol.l();
    let l_inv = invert_lower_triangular(&l).ok_or(BeamSolverError::MassNotPositiveDefinite)?;
    let l_inv_t = l_inv.transpose();
    let mut c = &l_inv * &k_ff * &l_inv_t;
    let c_t = c.transpose();
    c = (&c + &c_t) * 0.5;

    let eigen =
        nalgebra::SymmetricEigen::try_new(c, 1.0e-12, 0).ok_or(BeamSolverError::EigenFailed)?;
    let eigvals = &eigen.eigenvalues;
    let eigvecs = &eigen.eigenvectors;

    let mut order: Vec<usize> = (0..n_free).collect();
    order.sort_by(|&a, &b| {
        eigvals[a]
            .partial_cmp(&eigvals[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut modes = Vec::with_capacity(n_modes);
    for &idx in order.iter().take(n_modes) {
        let lambda = eigvals[idx].max(0.0);
        let omega = lambda.sqrt();
        let freq = omega / (2.0 * std::f64::consts::PI);
        let psi = eigvecs.column(idx).into_owned();
        let phi_free = &l_inv_t * &psi;
        let mut translation = vec![[0.0_f64; 3]; n_nodes];
        let mut rotation = vec![[0.0_f64; 3]; n_nodes];
        for (fi, &gd) in free.iter().enumerate() {
            let node = gd / 6;
            let comp = gd % 6;
            if comp < 3 {
                translation[node][comp] = phi_free[fi];
            } else {
                rotation[node][comp - 3] = phi_free[fi];
            }
        }
        modes.push(BeamMode {
            frequency_hz: freq,
            angular_frequency: omega,
            translation,
            rotation,
        });
    }
    Ok(BeamModalSolution { modes })
}

/// Invert a lower-triangular matrix by forward substitution. `None` if
/// a diagonal entry is too small to divide by.
fn invert_lower_triangular(l: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    let n = l.nrows();
    let mut inv = DMatrix::<f64>::zeros(n, n);
    for col in 0..n {
        for row in 0..n {
            let mut sum = if row == col { 1.0 } else { 0.0 };
            for k in 0..row {
                sum -= l[(row, k)] * inv[(k, col)];
            }
            let diag = l[(row, row)];
            if diag.abs() < 1.0e-300 {
                return None;
            }
            inv[(row, col)] = sum / diag;
        }
    }
    Some(inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steel() -> FemMaterial {
        FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            density: 7850.0,
            ..FemMaterial::default()
        }
    }

    #[test]
    fn rectangle_section_moments_are_correct() {
        let s = BeamSection::rectangle(0.1, 0.2);
        assert!((s.area - 0.02).abs() < 1e-12);
        // Iy = w·h³/12 = 0.1·0.008/12.
        assert!((s.iy - 0.1 * 0.2_f64.powi(3) / 12.0).abs() < 1e-15);
        assert!((s.iz - 0.2 * 0.1_f64.powi(3) / 12.0).abs() < 1e-15);
        assert!(s.j > 0.0);
    }

    #[test]
    fn circle_section_polar_moment() {
        let s = BeamSection::circle(0.05);
        let i = std::f64::consts::PI * 0.05_f64.powi(4) / 4.0;
        assert!((s.iy - i).abs() < 1e-18);
        assert!((s.j - 2.0 * i).abs() < 1e-18, "J should be the polar moment");
    }

    #[test]
    fn local_stiffness_is_symmetric() {
        let s = BeamSection::rectangle(0.1, 0.1);
        let k = beam_local_stiffness(2.0, 200e9, 80e9, &s);
        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (k[(i, j)] - k[(j, i)]).abs() < 1e-3 * k[(i, i)].abs().max(1.0),
                    "local K not symmetric at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn local_stiffness_has_six_rigid_body_modes() {
        // A free 3D beam element has six zero-energy rigid-body modes.
        let s = BeamSection::rectangle(0.1, 0.15);
        let l = 2.5;
        let k = beam_local_stiffness(l, 200e9, 80e9, &s);
        // Rigid translation along local x: both nodes' u DOF = 1.
        let mut tx = DVector::zeros(12);
        tx[0] = 1.0;
        tx[6] = 1.0;
        assert!((&k * &tx).norm() < 1e-3 * k.norm(), "axial rigid mode");
        // Rigid translation along local y.
        let mut ty = DVector::zeros(12);
        ty[1] = 1.0;
        ty[7] = 1.0;
        assert!((&k * &ty).norm() < 1e-3 * k.norm(), "transverse-y rigid mode");
        // Rigid rotation about local z: v = θz·x, so node 2 (at x=L)
        // gets v = L and both nodes get θz = 1.
        let mut rz = DVector::zeros(12);
        rz[5] = 1.0; // θz node 1
        rz[11] = 1.0; // θz node 2
        rz[7] = l; // v node 2
        assert!(
            (&k * &rz).norm() < 1e-3 * k.norm() * l,
            "rigid rotation about z gave force {}",
            (&k * &rz).norm()
        );
    }

    #[test]
    fn cantilever_tip_load_matches_analytic_euler_bernoulli() {
        // A cantilever along global X, clamped at node 0, with a
        // transverse tip load. With a slender section the Timoshenko
        // element reproduces δ = P·L³/(3·E·I) closely (shear adds a
        // small extra term P·L/(κ·G·A)).
        let mat = steel();
        let l = 4.0;
        let n_elem = 8;
        let nodes: Vec<Vector3<f64>> = (0..=n_elem)
            .map(|i| Vector3::new(l * i as f64 / n_elem as f64, 0.0, 0.0))
            .collect();
        let section = BeamSection::rectangle(0.05, 0.05);
        let elements: Vec<BeamElement> = (0..n_elem)
            .map(|i| BeamElement::new(i, i + 1, section))
            .collect();
        let p = 1000.0;
        let constraints = [BeamConstraint::clamped(0)];
        // Load in -Z at the tip.
        let loads = [BeamLoad::force(n_elem, [0.0, 0.0, -p])];
        let sol = solve_beam_static(&nodes, &elements, &mat, &constraints, &loads).unwrap();

        let tip = sol.translation[n_elem][2];
        // Bending about local y (the x-z plane). Iy = w·h³/12.
        let i = section.iy;
        let bending = cantilever_tip_deflection(p, l, mat.youngs_modulus, i);
        let g = mat.youngs_modulus / (2.0 * (1.0 + mat.poisson_ratio));
        let shear = p * l / (section.shear_z * g * section.area);
        let analytic = bending + shear;
        let rel = (tip.abs() - analytic).abs() / analytic;
        assert!(
            rel < 0.02,
            "tip deflection {} vs analytic {analytic} (rel {rel})",
            tip.abs()
        );
        assert!(tip < 0.0, "tip should deflect in -Z, got {tip}");
    }

    #[test]
    fn cantilever_udl_root_moment_matches_statics() {
        // Worked: w = 1 kN/m on a 2 m cantilever → M_root = w·L²/2 = 2000 N·m.
        let m = cantilever_udl_root_moment(1000.0, 2.0);
        assert!((m - 2000.0).abs() < 1e-9, "M_root = w·L²/2, got {m}");

        // Threads cantilever_udl_tip_deflection (δ = wL⁴/8EI = M·L²/4EI) and
        // cantilever_udl_strain_energy (U = w²L⁵/40EI = M²·L/10EI).
        for &(w, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 1.2, 120.0e9, 9.0e-8),
        ] {
            let from_moment = cantilever_udl_root_moment(w, l) * l * l / (4.0 * e * i);
            let delta = cantilever_udl_tip_deflection(w, l, e, i);
            assert!((from_moment - delta).abs() <= 1e-12 * delta.abs(), "M·L²/4EI = δ_udl");

            let from_moment_u = cantilever_udl_root_moment(w, l).powi(2) * l / (10.0 * e * i);
            let energy = cantilever_udl_strain_energy(w, l, e, i);
            assert!((from_moment_u - energy).abs() <= 1e-12 * energy.abs(), "M²·L/10EI = U");
        }

        // Quadratic in span; linear and sign-preserving in w.
        assert!(
            (cantilever_udl_root_moment(1000.0, 4.0) - 4.0 * cantilever_udl_root_moment(1000.0, 2.0))
                .abs()
                < 1e-9,
            "quadratic in L"
        );
        assert!(cantilever_udl_root_moment(-1000.0, 2.0) < 0.0, "sign follows the load");

        // Non-physical input → 0.
        assert_eq!(cantilever_udl_root_moment(f64::NAN, 2.0), 0.0);
        assert_eq!(cantilever_udl_root_moment(1000.0, 0.0), 0.0);
        assert_eq!(cantilever_udl_root_moment(1000.0, -1.0), 0.0);
    }

    #[test]
    fn cantilever_udl_max_shear_matches_statics() {
        // Worked: w = 1 kN/m on a 2 m cantilever → V_root = w·L = 2000 N.
        let v = cantilever_udl_max_shear(1000.0, 2.0);
        assert!((v - 2000.0).abs() < 1e-9, "V_root = w·L, got {v}");

        // STRONG non-tautological threads over signed (w, L): the root shear equals the total
        // distributed load, and the root moment is the shear acting at the load centroid L/2.
        for &(w, l) in &[(1000.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 1.2)] {
            let shear = cantilever_udl_max_shear(w, l);
            assert!((shear - w * l).abs() <= 1e-12 * (w * l).abs(), "V = total load w·L");
            let moment = cantilever_udl_root_moment(w, l);
            assert!((moment - shear * l / 2.0).abs() <= 1e-9 * moment.abs(), "M_root = V·L/2");
        }

        // Linear in w and L; sign follows the load.
        assert!(
            (cantilever_udl_max_shear(2000.0, 2.0) - 2.0 * cantilever_udl_max_shear(1000.0, 2.0))
                .abs()
                < 1e-9,
            "linear in w"
        );
        assert!(
            (cantilever_udl_max_shear(1000.0, 4.0) - 2.0 * cantilever_udl_max_shear(1000.0, 2.0))
                .abs()
                < 1e-9,
            "linear in L"
        );
        assert!(cantilever_udl_max_shear(-1000.0, 2.0) < 0.0, "sign follows the load");

        // Non-physical input → 0.
        assert_eq!(cantilever_udl_max_shear(f64::NAN, 2.0), 0.0);
        assert_eq!(cantilever_udl_max_shear(1000.0, 0.0), 0.0);
        assert_eq!(cantilever_udl_max_shear(1000.0, -1.0), 0.0);
    }

    #[test]
    fn cantilever_point_load_root_moment_matches_statics() {
        // Worked: P = 1 kN at the tip of a 2 m cantilever → M_root = P·L = 2000 N·m.
        let m = cantilever_point_load_root_moment(1000.0, 2.0);
        assert!((m - 2000.0).abs() < 1e-9, "M_root = P·L, got {m}");

        // Threads cantilever_tip_deflection (δ = PL³/3EI = M·L²/3EI) and
        // cantilever_point_load_strain_energy (U = P²L³/6EI = M²·L/6EI).
        for &(p, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 1.2, 120.0e9, 9.0e-8),
        ] {
            let from_moment = cantilever_point_load_root_moment(p, l) * l * l / (3.0 * e * i);
            let delta = cantilever_tip_deflection(p, l, e, i);
            assert!((from_moment - delta).abs() <= 1e-12 * delta.abs(), "M·L²/3EI = δ_tip");

            let from_moment_u =
                cantilever_point_load_root_moment(p, l).powi(2) * l / (6.0 * e * i);
            let energy = cantilever_point_load_strain_energy(p, l, e, i);
            assert!((from_moment_u - energy).abs() <= 1e-12 * energy.abs(), "M²·L/6EI = U");
        }

        // Linear and sign-preserving in P; linear in L.
        assert!(cantilever_point_load_root_moment(-1000.0, 2.0) < 0.0, "sign follows the load");
        assert!(
            (cantilever_point_load_root_moment(1000.0, 4.0)
                - 2.0 * cantilever_point_load_root_moment(1000.0, 2.0))
            .abs()
                < 1e-9,
            "linear in L"
        );

        // Non-physical input → 0.
        assert_eq!(cantilever_point_load_root_moment(f64::NAN, 2.0), 0.0);
        assert_eq!(cantilever_point_load_root_moment(1000.0, 0.0), 0.0);
        assert_eq!(cantilever_point_load_root_moment(1000.0, -1.0), 0.0);
    }

    #[test]
    fn cantilever_end_moment_tip_deflection_matches_constant_curvature() {
        // (a) WORKED: M₀ = 1 kN·m at the free end of a 2 m steel cantilever (E = 200 GPa,
        // I = 1e-6 m⁴) → δ = M₀·L²/(2EI) = 1000·4/(2·200e9·1e-6) = 0.01 m.
        assert!(
            (cantilever_end_moment_tip_deflection(1000.0, 2.0, 200.0e9, 1.0e-6) - 0.01).abs()
                <= 1e-9 * 0.01,
            "δ = M₀L²/(2EI) = 0.01 m"
        );

        // (b) THREAD beam_curvature (non-tautological): a pure end couple gives a constant
        // curvature κ = M₀/(EI), so δ = κ·L²/2.
        for &(m, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-820.0, 3.5, 70.0e9, 4.2e-7),
        ] {
            let d = cantilever_end_moment_tip_deflection(m, l, e, i);
            assert!(
                (d - beam_curvature(m, e, i) * l * l / 2.0).abs() <= 1e-9 * d.abs().max(1e-300),
                "δ = κ·L²/2"
            );
        }

        // (c) SCALING: linear and sign-preserving in M₀, quadratic in L.
        let base = cantilever_end_moment_tip_deflection(1000.0, 2.0, 200.0e9, 1.0e-6);
        assert!(
            (cantilever_end_moment_tip_deflection(2000.0, 2.0, 200.0e9, 1.0e-6) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in M₀"
        );
        assert!(
            (cantilever_end_moment_tip_deflection(1000.0, 4.0, 200.0e9, 1.0e-6) - 4.0 * base).abs()
                <= 1e-9 * 4.0 * base,
            "quadratic in L"
        );
        assert!(
            cantilever_end_moment_tip_deflection(-1000.0, 2.0, 200.0e9, 1.0e-6) < 0.0,
            "sign follows the moment"
        );

        // (d) GUARD: non-physical input → 0.
        assert_eq!(cantilever_end_moment_tip_deflection(f64::NAN, 2.0, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_tip_deflection(1000.0, 0.0, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_tip_deflection(1000.0, 2.0, 0.0, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_tip_deflection(1000.0, 2.0, 200.0e9, 0.0), 0.0);
    }

    #[test]
    fn cantilever_end_moment_tip_slope_matches_constant_curvature() {
        // (a) WORKED: M₀ = 1 kN·m at the free end of a 2 m steel cantilever (E = 200 GPa,
        // I = 1e-6 m⁴) → θ = M₀·L/(EI) = 1000·2/(200e9·1e-6) = 0.01 rad.
        assert!(
            (cantilever_end_moment_tip_slope(1000.0, 2.0, 200.0e9, 1.0e-6) - 0.01).abs()
                <= 1e-9 * 0.01,
            "θ = M₀L/(EI) = 0.01 rad"
        );

        // (b) THREAD beam_curvature (non-tautological): a pure end couple gives a constant
        // curvature κ = M₀/(EI), so θ = κ·L.
        for &(m, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-820.0, 3.5, 70.0e9, 4.2e-7),
        ] {
            let t = cantilever_end_moment_tip_slope(m, l, e, i);
            assert!(
                (t - beam_curvature(m, e, i) * l).abs() <= 1e-9 * t.abs().max(1e-300),
                "θ = κ·L"
            );
        }

        // (c) THREAD cantilever_end_moment_tip_deflection (non-tautological): the constant-
        // curvature arc locks slope and deflection by δ = θ·L/2, so θ = 2·δ/L.
        let (m, l, e, i) = (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64);
        assert!(
            (cantilever_end_moment_tip_slope(m, l, e, i)
                - 2.0 * cantilever_end_moment_tip_deflection(m, l, e, i) / l)
            .abs()
                <= 1e-9 * cantilever_end_moment_tip_slope(m, l, e, i),
            "θ = 2δ/L"
        );

        // (d) SCALING: linear and sign-preserving in M₀, linear in L.
        let base = cantilever_end_moment_tip_slope(1000.0, 2.0, 200.0e9, 1.0e-6);
        assert!(
            (cantilever_end_moment_tip_slope(2000.0, 2.0, 200.0e9, 1.0e-6) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in M₀"
        );
        assert!(
            (cantilever_end_moment_tip_slope(1000.0, 4.0, 200.0e9, 1.0e-6) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in L"
        );
        assert!(
            cantilever_end_moment_tip_slope(-1000.0, 2.0, 200.0e9, 1.0e-6) < 0.0,
            "sign follows the moment"
        );

        // (e) GUARD: non-physical input → 0.
        assert_eq!(cantilever_end_moment_tip_slope(f64::NAN, 2.0, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_tip_slope(1000.0, 0.0, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_tip_slope(1000.0, 2.0, 0.0, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_tip_slope(1000.0, 2.0, 200.0e9, 0.0), 0.0);
    }

    #[test]
    fn cantilever_end_moment_strain_energy_matches_clapeyron() {
        // (a) WORKED: M₀ = 1 kN·m at the free end of a 2 m steel cantilever (E = 200 GPa,
        // I = 1e-6 m⁴) → U = M₀²·L/(2EI) = 1000²·2/(2·200e9·1e-6) = 5.0 J.
        assert!(
            (cantilever_end_moment_strain_energy(1000.0, 2.0, 200.0e9, 1.0e-6) - 5.0).abs()
                <= 1e-9 * 5.0,
            "U = M₀²L/(2EI) = 5.0 J"
        );

        // (b) CLAPEYRON THREAD (non-tautological): the work a single static load does is
        // half the load times the conjugate displacement, U = ½·M₀·θ with the tip slope.
        for &(m, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-820.0, 3.5, 70.0e9, 4.2e-7),
        ] {
            let u = cantilever_end_moment_strain_energy(m, l, e, i);
            assert!(
                (u - 0.5 * m * cantilever_end_moment_tip_slope(m, l, e, i)).abs() <= 1e-9 * u,
                "U = ½·M₀·θ"
            );
        }

        // (c) SCALING: quadratic and sign-independent in M₀ (the M² term), linear in L.
        let base = cantilever_end_moment_strain_energy(1000.0, 2.0, 200.0e9, 1.0e-6);
        assert!(
            (cantilever_end_moment_strain_energy(2000.0, 2.0, 200.0e9, 1.0e-6) - 4.0 * base).abs()
                <= 1e-9 * 4.0 * base,
            "quadratic in M₀"
        );
        assert!(
            (cantilever_end_moment_strain_energy(-1000.0, 2.0, 200.0e9, 1.0e-6) - base).abs()
                <= 1e-9 * base,
            "sign-independent (M²)"
        );
        assert!(
            (cantilever_end_moment_strain_energy(1000.0, 4.0, 200.0e9, 1.0e-6) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in L"
        );

        // (d) GUARD: non-physical input → 0.
        assert_eq!(cantilever_end_moment_strain_energy(f64::NAN, 2.0, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_strain_energy(1000.0, 0.0, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_strain_energy(1000.0, 2.0, 0.0, 1.0e-6), 0.0);
        assert_eq!(cantilever_end_moment_strain_energy(1000.0, 2.0, 200.0e9, 0.0), 0.0);
    }

    #[test]
    fn cantilever_tip_deflection_matches_the_closed_form() {
        // Worked point: P = 1 kN at the tip of a 2 m steel cantilever,
        // E = 200 GPa, I = 1e-6 m⁴ → δ = P·L³/(3·E·I) = 8000/6e5 = 1/75 ≈ 0.01333 m.
        let (p, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let delta = cantilever_tip_deflection(p, l, e, i);
        assert!((delta - 1.0 / 75.0).abs() / delta < 1e-9, "δ = 1/75 m, got {delta}");
        // Linear in the load, and sign-preserving (an upward load lifts the tip).
        assert!((cantilever_tip_deflection(2.0 * p, l, e, i) - 2.0 * delta).abs() / delta < 1e-12);
        assert!((cantilever_tip_deflection(-p, l, e, i) + delta).abs() / delta < 1e-12, "sign-preserving");
        // Cubic in the span: double L → 8× δ.
        assert!((cantilever_tip_deflection(p, 2.0 * l, e, i) - 8.0 * delta).abs() / delta < 1e-9, "L³ scaling");
        // Inverse in the flexural rigidity E·I: double E or I → half δ.
        assert!((cantilever_tip_deflection(p, l, 2.0 * e, i) - 0.5 * delta).abs() / delta < 1e-12, "1/E");
        assert!((cantilever_tip_deflection(p, l, e, 2.0 * i) - 0.5 * delta).abs() / delta < 1e-12, "1/I");
        // Non-physical input → 0.
        assert_eq!(cantilever_tip_deflection(p, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(cantilever_tip_deflection(p, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(cantilever_tip_deflection(p, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(cantilever_tip_deflection(f64::NAN, l, e, i), 0.0); // non-finite P
        assert_eq!(cantilever_tip_deflection(p, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn cantilever_point_load_strain_energy_matches_clapeyron() {
        // Worked point: P = 1 kN at the tip of a 2 m steel cantilever, E = 200 GPa,
        // I = 1e-6 m⁴ → U = P²·L³/(6·E·I) = 8e6/1.2e6 = 20/3 ≈ 6.6667 J.
        let (p, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let u = cantilever_point_load_strain_energy(p, l, e, i);
        assert!((u - 20.0 / 3.0).abs() / u < 1e-9, "U = 20/3 J, got {u}");
        // Quadratic in the load → SIGN-INDEPENDENT: an up- or down-load store equal
        // energy (unlike deflection, which is sign-preserving).
        assert!((cantilever_point_load_strain_energy(2.0 * p, l, e, i) - 4.0 * u).abs() / u < 1e-9, "P² scaling");
        assert!((cantilever_point_load_strain_energy(-p, l, e, i) - u).abs() / u < 1e-12, "sign-independent");
        // Cubic in the span: double L → 8× U.
        assert!((cantilever_point_load_strain_energy(p, 2.0 * l, e, i) - 8.0 * u).abs() / u < 1e-9, "L³ scaling");
        // Inverse in the flexural rigidity E·I: double E or I → half U.
        assert!((cantilever_point_load_strain_energy(p, l, 2.0 * e, i) - 0.5 * u).abs() / u < 1e-12, "1/E");
        assert!((cantilever_point_load_strain_energy(p, l, e, 2.0 * i) - 0.5 * u).abs() / u < 1e-12, "1/I");
        // STRONG non-tautological cross-check via Clapeyron's theorem: the work a
        // single static load does is half the load times the deflection it makes,
        // U = ½·P·δ_tip. The energy impl is ∫M²/(2EI) → P²L³/(6EI); the check uses the
        // independent deflection fn (P·L³/(3EI)) — different path, same value.
        for &(pp, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (-800.0, 3.0), (4200.0, 0.8)] {
            let energy = cantilever_point_load_strain_energy(pp, ll, e, i);
            let half_p_delta = 0.5 * pp * cantilever_tip_deflection(pp, ll, e, i);
            assert!(
                (energy - half_p_delta).abs() / energy < 1e-12,
                "Clapeyron U = ½·P·δ at P={pp}, L={ll}: {energy} vs {half_p_delta}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(cantilever_point_load_strain_energy(p, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(cantilever_point_load_strain_energy(p, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(cantilever_point_load_strain_energy(p, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(cantilever_point_load_strain_energy(f64::NAN, l, e, i), 0.0); // non-finite P
        assert_eq!(cantilever_point_load_strain_energy(p, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn simply_supported_point_load_strain_energy_matches_clapeyron() {
        // Worked point: P = 1 kN central load on a 2 m simply-supported steel beam,
        // E = 200 GPa, I = 1e-6 m⁴ → U = P²·L³/(96·E·I) = 8e6/1.92e7 = 5/12 ≈ 0.4167 J.
        let (p, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let u = simply_supported_point_load_strain_energy(p, l, e, i);
        assert!((u - 5.0 / 12.0).abs() / u < 1e-9, "U = 5/12 J, got {u}");
        // Quadratic in the load → SIGN-INDEPENDENT; cubic in span; inverse in E·I.
        assert!((simply_supported_point_load_strain_energy(2.0 * p, l, e, i) - 4.0 * u).abs() / u < 1e-9, "P² scaling");
        assert!((simply_supported_point_load_strain_energy(-p, l, e, i) - u).abs() / u < 1e-12, "sign-independent");
        assert!((simply_supported_point_load_strain_energy(p, 2.0 * l, e, i) - 8.0 * u).abs() / u < 1e-9, "L³ scaling");
        assert!((simply_supported_point_load_strain_energy(p, l, 2.0 * e, i) - 0.5 * u).abs() / u < 1e-12, "1/E");
        assert!((simply_supported_point_load_strain_energy(p, l, e, 2.0 * i) - 0.5 * u).abs() / u < 1e-12, "1/I");
        // STRONG cross-check (1) via Clapeyron's theorem U = ½·P·δ_centre, threading
        // the independent simply_supported_center_deflection (δ = P·L³/(48EI)).
        for &(pp, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (-800.0, 3.0), (4200.0, 0.8)] {
            let energy = simply_supported_point_load_strain_energy(pp, ll, e, i);
            let half_p_delta = 0.5 * pp * simply_supported_center_deflection(pp, ll, e, i);
            assert!(
                (energy - half_p_delta).abs() / energy < 1e-12,
                "Clapeyron U = ½·P·δ at P={pp}, L={ll}: {energy} vs {half_p_delta}"
            );
        }
        // STRONG cross-check (2): for equal P,L,E,I a cantilever stores exactly 16× the
        // energy of the propped span (P²L³/(6EI) vs P²L³/(96EI)), threading #236.
        for &(pp, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (4200.0, 0.8)] {
            let ss = simply_supported_point_load_strain_energy(pp, ll, e, i);
            let cant = cantilever_point_load_strain_energy(pp, ll, e, i);
            assert!((cant - 16.0 * ss).abs() / cant < 1e-12, "cantilever = 16× SS at P={pp}, L={ll}");
        }
        // Non-physical input → 0.
        assert_eq!(simply_supported_point_load_strain_energy(p, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(simply_supported_point_load_strain_energy(p, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(simply_supported_point_load_strain_energy(p, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(simply_supported_point_load_strain_energy(f64::NAN, l, e, i), 0.0); // non-finite P
        assert_eq!(simply_supported_point_load_strain_energy(p, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn fixed_fixed_point_load_strain_energy_completes_the_set() {
        // Worked: P = 1 kN, L = 4 m, E = 200 GPa, I = 1e-6 m⁴ →
        // U = P²L³/(384EI) = 1e6·64/(384·2e5) = 0.8333… J.
        let (p, l, e, i) = (1000.0_f64, 4.0_f64, 200.0e9_f64, 1.0e-6_f64);
        let u = fixed_fixed_point_load_strain_energy(p, l, e, i);
        assert!((u - 64.0e6 / 7.68e7).abs() <= 1e-9 * u, "U = P²L³/384EI");

        // STRONG non-tautological threads: (i) Clapeyron U = ½·P·δ with the centre-deflection
        // fn; (ii) ratio identities — 1/64 of the cantilever and 1/4 of the simply-supported
        // strain energy under the same central load.
        assert!(
            (u - 0.5 * p * fixed_fixed_center_deflection(p, l, e, i)).abs() <= 1e-9 * u,
            "U = ½·P·δ (Clapeyron, threads centre deflection)",
        );
        assert!(
            (u - cantilever_point_load_strain_energy(p, l, e, i) / 64.0).abs() <= 1e-9 * u,
            "fixed-fixed = 1/64 × cantilever",
        );
        assert!(
            (u - simply_supported_point_load_strain_energy(p, l, e, i) / 4.0).abs() <= 1e-9 * u,
            "fixed-fixed = 1/4 × simply-supported",
        );

        // Scaling: ∝ P², ∝ L³; sign-independent.
        assert!(
            (fixed_fixed_point_load_strain_energy(2.0 * p, l, e, i) - 4.0 * u).abs() <= 1e-9 * u,
            "∝ P²",
        );
        assert!(
            (fixed_fixed_point_load_strain_energy(p, 2.0 * l, e, i) - 8.0 * u).abs() <= 1e-9 * u,
            "∝ L³",
        );
        assert!(
            (fixed_fixed_point_load_strain_energy(-p, l, e, i) - u).abs() <= 1e-12 * u,
            "sign-independent",
        );

        // Guards: non-physical input → 0.
        assert_eq!(fixed_fixed_point_load_strain_energy(p, l, 0.0, i), 0.0);
        assert_eq!(fixed_fixed_point_load_strain_energy(p, l, e, -1.0e-6), 0.0);
        assert_eq!(fixed_fixed_point_load_strain_energy(p, -1.0, e, i), 0.0);
        assert_eq!(fixed_fixed_point_load_strain_energy(f64::NAN, l, e, i), 0.0);
    }

    #[test]
    fn cantilever_tip_slope_matches_closed_form() {
        // Worked point: P = 1 kN at the tip of a 2 m steel cantilever, E = 200 GPa,
        // I = 1e-6 m⁴ → θ = P·L²/(2·E·I) = 4000/4e5 = 0.01 rad.
        let (p, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let theta = cantilever_tip_slope(p, l, e, i);
        assert!((theta - 0.01).abs() / theta < 1e-12, "θ = 0.01 rad, got {theta}");
        // Linear in the load, and sign-preserving (an upward load rotates the tip up).
        assert!((cantilever_tip_slope(2.0 * p, l, e, i) - 2.0 * theta).abs() / theta < 1e-12);
        assert!((cantilever_tip_slope(-p, l, e, i) + theta).abs() / theta < 1e-12, "sign-preserving");
        // Quadratic in the span: double L → 4× θ.
        assert!((cantilever_tip_slope(p, 2.0 * l, e, i) - 4.0 * theta).abs() / theta < 1e-12, "L² scaling");
        // Inverse in the flexural rigidity E·I: double E or I → half θ.
        assert!((cantilever_tip_slope(p, l, 2.0 * e, i) - 0.5 * theta).abs() / theta < 1e-12, "1/E");
        assert!((cantilever_tip_slope(p, l, e, 2.0 * i) - 0.5 * theta).abs() / theta < 1e-12, "1/I");
        // STRONG non-tautological cross-check tying it to #176: for a tip-loaded
        // cantilever the tip deflection is exactly two-thirds of the span times the
        // tip slope, δ = (2/3)·L·θ. The slope impl uses L²/(2EI); the check uses the
        // independent deflection fn L³/(3EI) — a known beam relation, different path.
        for &(pp, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (-800.0, 3.0), (4200.0, 0.8)] {
            let slope = cantilever_tip_slope(pp, ll, e, i);
            let defl = cantilever_tip_deflection(pp, ll, e, i);
            assert!(
                (defl - 2.0 / 3.0 * ll * slope).abs() / defl.abs() < 1e-12,
                "δ = (2/3)·L·θ at P={pp}, L={ll}: {defl} vs {slope}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(cantilever_tip_slope(p, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(cantilever_tip_slope(p, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(cantilever_tip_slope(p, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(cantilever_tip_slope(f64::NAN, l, e, i), 0.0); // non-finite P
        assert_eq!(cantilever_tip_slope(p, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn cantilever_udl_tip_deflection_matches_closed_form() {
        // Worked point: w = 1 kN/m uniformly along a 2 m steel cantilever,
        // E = 200 GPa, I = 1e-6 m⁴ → δ = w·L⁴/(8·E·I) = 16000/1.6e6 = 0.01 m.
        let (w, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let delta = cantilever_udl_tip_deflection(w, l, e, i);
        assert!((delta - 0.01).abs() / delta < 1e-9, "δ = 0.01 m, got {delta}");
        // Linear in the load intensity, and sign-preserving (an upward load lifts the tip).
        assert!((cantilever_udl_tip_deflection(2.0 * w, l, e, i) - 2.0 * delta).abs() / delta < 1e-12);
        assert!((cantilever_udl_tip_deflection(-w, l, e, i) + delta).abs() / delta < 1e-12, "sign-preserving");
        // Quartic in the span: double L → 16× δ.
        assert!((cantilever_udl_tip_deflection(w, 2.0 * l, e, i) - 16.0 * delta).abs() / delta < 1e-9, "L⁴ scaling");
        // Inverse in the flexural rigidity E·I: double E or I → half δ.
        assert!((cantilever_udl_tip_deflection(w, l, 2.0 * e, i) - 0.5 * delta).abs() / delta < 1e-12, "1/E");
        assert!((cantilever_udl_tip_deflection(w, l, e, 2.0 * i) - 0.5 * delta).abs() / delta < 1e-12, "1/I");
        // STRONG non-tautological cross-check: the same TOTAL load W = w·L, spread as a
        // UDL, deflects the tip to exactly 3/8 of the same total concentrated at the
        // tip. The UDL impl uses w·L⁴/(8EI); the check uses the independent point-load
        // fn cantilever_tip_deflection (P·L³/(3EI)) with P = w·L — a known
        // structural-mechanics ratio, a different code path.
        for &(ww, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (-800.0, 3.0), (4200.0, 0.8)] {
            let udl = cantilever_udl_tip_deflection(ww, ll, e, i);
            let point = cantilever_tip_deflection(ww * ll, ll, e, i);
            assert!(
                (udl - 3.0 / 8.0 * point).abs() / point.abs() < 1e-12,
                "UDL = 3/8 × point-load(W=w·L) at w={ww}, L={ll}: {udl} vs {point}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(cantilever_udl_tip_deflection(w, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(cantilever_udl_tip_deflection(w, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(cantilever_udl_tip_deflection(w, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(cantilever_udl_tip_deflection(f64::NAN, l, e, i), 0.0); // non-finite w
        assert_eq!(cantilever_udl_tip_deflection(w, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn cantilever_udl_strain_energy_matches_the_moment_integral() {
        // Worked point: w = 1 kN/m UDL on a 2 m steel cantilever, E = 200 GPa,
        // I = 1e-6 m⁴ → U = w²·L⁵/(40·E·I) = 3.2e7/8e6 = 4.0 J.
        let (w, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let u = cantilever_udl_strain_energy(w, l, e, i);
        assert!((u - 4.0).abs() / u < 1e-9, "U = 4.0 J, got {u}");
        // Quadratic in w → SIGN-INDEPENDENT; L⁵ scaling; inverse in E·I.
        assert!((cantilever_udl_strain_energy(2.0 * w, l, e, i) - 4.0 * u).abs() / u < 1e-9, "w² scaling");
        assert!((cantilever_udl_strain_energy(-w, l, e, i) - u).abs() / u < 1e-12, "sign-independent");
        assert!((cantilever_udl_strain_energy(w, 2.0 * l, e, i) - 32.0 * u).abs() / u < 1e-9, "L⁵ scaling");
        assert!((cantilever_udl_strain_energy(w, l, 2.0 * e, i) - 0.5 * u).abs() / u < 1e-12, "1/E");
        assert!((cantilever_udl_strain_energy(w, l, e, 2.0 * i) - 0.5 * u).abs() / u < 1e-12, "1/I");
        // STRONG cross-check (1): a first-principles NUMERICAL Riemann integral of the
        // bending energy ∫₀^L M(x)²/(2EI) dx with the cantilever-UDL moment
        // M(x) = w(L−x)²/2 — an independent derivation of the closed form.
        let n = 100_000;
        let dx = l / n as f64;
        let mut energy_sum = 0.0;
        for k in 0..n {
            let x = (k as f64 + 0.5) * dx; // midpoint
            let m = w * (l - x) * (l - x) / 2.0;
            energy_sum += m * m / (2.0 * e * i) * dx;
        }
        assert!((u - energy_sum).abs() / u < 1e-5, "U = ∫M²/(2EI)dx numerically: {u} vs {energy_sum}");
        // STRONG cross-check (2): U = (1/5)·(w·L)·δ_tip, threading the independent
        // cantilever_udl_tip_deflection (δ = wL⁴/(8EI)).
        for &(ww, ll) in &[(1000.0_f64, 2.0_f64), (300.0, 1.5), (-750.0, 3.0)] {
            let energy = cantilever_udl_strain_energy(ww, ll, e, i);
            let from_defl = 0.2 * ww * ll * cantilever_udl_tip_deflection(ww, ll, e, i);
            assert!(
                (energy - from_defl).abs() / energy < 1e-12,
                "U = (1/5)·wL·δ_tip at w={ww}, L={ll}: {energy} vs {from_defl}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(cantilever_udl_strain_energy(w, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(cantilever_udl_strain_energy(w, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(cantilever_udl_strain_energy(w, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(cantilever_udl_strain_energy(f64::NAN, l, e, i), 0.0); // non-finite w
        assert_eq!(cantilever_udl_strain_energy(w, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn simply_supported_udl_strain_energy_completes_the_energy_matrix() {
        // Worked point: w = 1 kN/m UDL on a 2 m simply-supported steel beam, E = 200
        // GPa, I = 1e-6 m⁴ → U = w²·L⁵/(240·E·I) = 3.2e7/4.8e7 = 2/3 ≈ 0.6667 J.
        let (w, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let u = simply_supported_udl_strain_energy(w, l, e, i);
        assert!((u - 2.0 / 3.0).abs() / u < 1e-9, "U = 2/3 J, got {u}");
        // Quadratic in w → SIGN-INDEPENDENT; L⁵ scaling; inverse in E·I.
        assert!((simply_supported_udl_strain_energy(2.0 * w, l, e, i) - 4.0 * u).abs() / u < 1e-9, "w² scaling");
        assert!((simply_supported_udl_strain_energy(-w, l, e, i) - u).abs() / u < 1e-12, "sign-independent");
        assert!((simply_supported_udl_strain_energy(w, 2.0 * l, e, i) - 32.0 * u).abs() / u < 1e-9, "L⁵ scaling");
        assert!((simply_supported_udl_strain_energy(w, l, 2.0 * e, i) - 0.5 * u).abs() / u < 1e-12, "1/E");
        assert!((simply_supported_udl_strain_energy(w, l, e, 2.0 * i) - 0.5 * u).abs() / u < 1e-12, "1/I");
        // STRONG cross-check (1): a first-principles NUMERICAL Riemann integral of the
        // bending energy ∫₀^L M(x)²/(2EI) dx with the SS-UDL moment M(x) = (w/2)·x·(L−x).
        let n = 100_000;
        let dx = l / n as f64;
        let mut energy_sum = 0.0;
        for k in 0..n {
            let x = (k as f64 + 0.5) * dx; // midpoint
            let m = 0.5 * w * x * (l - x);
            energy_sum += m * m / (2.0 * e * i) * dx;
        }
        assert!((u - energy_sum).abs() / u < 1e-5, "U = ∫M²/(2EI)dx numerically: {u} vs {energy_sum}");
        // STRONG cross-check (2): completing the energy 2×2 matrix — for equal w,L,E,I
        // a cantilever stores exactly 6× the simply-supported UDL energy (threads #248).
        for &(ww, ll) in &[(1000.0_f64, 2.0_f64), (300.0, 1.5), (-750.0, 3.0)] {
            let ss = simply_supported_udl_strain_energy(ww, ll, e, i);
            let cant = cantilever_udl_strain_energy(ww, ll, e, i);
            assert!((cant - 6.0 * ss).abs() / cant < 1e-12, "cantilever = 6× SS UDL at w={ww}, L={ll}");
        }
        // Non-physical input → 0.
        assert_eq!(simply_supported_udl_strain_energy(w, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(simply_supported_udl_strain_energy(w, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(simply_supported_udl_strain_energy(w, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(simply_supported_udl_strain_energy(f64::NAN, l, e, i), 0.0); // non-finite w
        assert_eq!(simply_supported_udl_strain_energy(w, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn cantilever_udl_tip_slope_matches_closed_form() {
        // Worked point: w = 1 kN/m uniformly along a 2 m steel cantilever, E = 200 GPa,
        // I = 1e-6 m⁴ → θ = w·L³/(6·E·I) = 8000/1.2e6 = 1/150 rad.
        let (w, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let theta = cantilever_udl_tip_slope(w, l, e, i);
        assert!((theta - 1.0 / 150.0).abs() / theta < 1e-12, "θ = 1/150 rad, got {theta}");
        // Linear in the load intensity (sign-preserving); cubic in span; inverse in E·I.
        assert!((cantilever_udl_tip_slope(2.0 * w, l, e, i) - 2.0 * theta).abs() / theta < 1e-12);
        assert!((cantilever_udl_tip_slope(-w, l, e, i) + theta).abs() / theta < 1e-12, "sign-preserving");
        assert!((cantilever_udl_tip_slope(w, 2.0 * l, e, i) - 8.0 * theta).abs() / theta < 1e-12, "L³ scaling");
        assert!((cantilever_udl_tip_slope(w, l, 2.0 * e, i) - 0.5 * theta).abs() / theta < 1e-12, "1/E");
        assert!((cantilever_udl_tip_slope(w, l, e, 2.0 * i) - 0.5 * theta).abs() / theta < 1e-12, "1/I");
        // Cross-check (a) tying it to the UDL deflection #200: δ = (3/4)·L·θ.
        for &(ww, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (-800.0, 3.0)] {
            let slope = cantilever_udl_tip_slope(ww, ll, e, i);
            let defl = cantilever_udl_tip_deflection(ww, ll, e, i);
            assert!(
                (defl - 3.0 / 4.0 * ll * slope).abs() / defl.abs() < 1e-12,
                "δ = (3/4)·L·θ at w={ww}, L={ll}"
            );
        }
        // Cross-check (b) tying it to the point-load slope #212: the same TOTAL load
        // W = w·L spread as a UDL rotates the tip to 1/3 of the concentrated case.
        for &(ww, ll) in &[(1000.0_f64, 2.0_f64), (250.0, 1.5), (-800.0, 3.0)] {
            let udl = cantilever_udl_tip_slope(ww, ll, e, i);
            let point = cantilever_tip_slope(ww * ll, ll, e, i);
            assert!(
                (udl - point / 3.0).abs() / udl.abs() < 1e-12,
                "θ_udl = θ_point(W=w·L)/3 at w={ww}, L={ll}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(cantilever_udl_tip_slope(w, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(cantilever_udl_tip_slope(w, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(cantilever_udl_tip_slope(w, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(cantilever_udl_tip_slope(f64::NAN, l, e, i), 0.0); // non-finite w
        assert_eq!(cantilever_udl_tip_slope(w, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn cantilever_axial_load_matches_analytic() {
        // Axial extension δ = F·L/(E·A).
        let mat = steel();
        let l = 3.0;
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(l, 0.0, 0.0),
        ];
        let section = BeamSection::circle(0.02);
        let elements = [BeamElement::new(0, 1, section)];
        let f = 5.0e4;
        let constraints = [BeamConstraint::clamped(0)];
        let loads = [BeamLoad::force(1, [f, 0.0, 0.0])];
        let sol = solve_beam_static(&nodes, &elements, &mat, &constraints, &loads).unwrap();
        let analytic = beam_axial_extension(f, l, mat.youngs_modulus, section.area);
        let rel = (sol.translation[1][0] - analytic).abs() / analytic;
        assert!(rel < 1e-6, "axial δ {} vs {analytic}", sol.translation[1][0]);
    }

    #[test]
    fn elastic_section_modulus_threads_the_flexure_formula() {
        // Worked: S = I/c = 1e-6 / 0.05 = 2e-5 m³.
        let s = elastic_section_modulus(1.0e-6, 0.05);
        assert!((s - 2.0e-5).abs() <= 1e-12 * 2.0e-5, "S = I/c, got {s}");

        // Threads bending_stress: the peak stress is the moment over the section modulus.
        for &(m, c, i) in &[
            (2000.0_f64, 0.05_f64, 1.0e-6_f64),
            (-450.0, 0.02, 4.2e-7),
            (8200.0, 0.12, 9.0e-8),
        ] {
            let from_modulus = m / elastic_section_modulus(i, c);
            assert!(
                (bending_stress(m, c, i) - from_modulus).abs() <= 1e-12 * from_modulus.abs(),
                "σ_max = M/S"
            );
        }

        // Linear in I, inverse in c.
        let base = elastic_section_modulus(1.0e-6, 0.05);
        assert!(
            (elastic_section_modulus(2.0e-6, 0.05) - 2.0 * base).abs() <= 1e-12 * base,
            "linear in I"
        );
        assert!(
            (elastic_section_modulus(1.0e-6, 0.10) - 0.5 * base).abs() <= 1e-12 * base,
            "inverse in c"
        );

        // Non-physical input → 0.
        assert_eq!(elastic_section_modulus(0.0, 0.05), 0.0);
        assert_eq!(elastic_section_modulus(1.0e-6, 0.0), 0.0);
        assert_eq!(elastic_section_modulus(1.0e-6, -0.05), 0.0);
        assert_eq!(elastic_section_modulus(f64::NAN, 0.05), 0.0);
    }

    #[test]
    fn rectangular_second_moment_of_area_is_bh3_over_12() {
        // (a) WORKED: b=0.1, h=0.2 → I = b·h³/12 = 6.6667e-5 m⁴.
        let i = rectangular_second_moment_of_area(0.1, 0.2);
        assert!((i - 0.1 * 0.2_f64.powi(3) / 12.0).abs() <= 1e-9 * i, "I = b·h³/12");

        // (b) THREAD elastic_section_modulus (#380) (non-tautological): the rectangular
        // section modulus is S = I/(h/2) = b·h²/6.
        let (b, h) = (0.1_f64, 0.2_f64);
        assert!(
            (elastic_section_modulus(rectangular_second_moment_of_area(b, h), h / 2.0)
                - b * h * h / 6.0)
            .abs()
                <= 1e-9 * (b * h * h / 6.0),
            "S = I/(h/2) = b·h²/6"
        );

        // (c) THREAD flexural_rigidity (#398): EI of a steel rectangle.
        let ei = 200.0e9 * rectangular_second_moment_of_area(0.1, 0.2);
        assert!(
            (flexural_rigidity(200.0e9, rectangular_second_moment_of_area(0.1, 0.2)) - ei).abs()
                <= 1e-6 * ei,
            "EI = E·I"
        );

        // (d) SCALING: linear in width, cubic in height.
        let base = rectangular_second_moment_of_area(0.1, 0.2);
        assert!(
            (rectangular_second_moment_of_area(0.2, 0.2) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in width"
        );
        assert!(
            (rectangular_second_moment_of_area(0.1, 0.4) - 8.0 * base).abs() <= 1e-9 * 8.0 * base,
            "cubic in height"
        );

        // (e) GUARD: non-positive or non-finite → 0.
        assert_eq!(rectangular_second_moment_of_area(0.0, 0.2), 0.0);
        assert_eq!(rectangular_second_moment_of_area(0.1, 0.0), 0.0);
        assert_eq!(rectangular_second_moment_of_area(f64::NAN, 0.2), 0.0);
        assert_eq!(rectangular_second_moment_of_area(0.1, -0.2), 0.0);
    }

    #[test]
    fn rectangular_plastic_section_modulus_is_bh2_over_4() {
        // (a) WORKED: b=0.1, h=0.2 → Z = b·h²/4 = 0.001 m³.
        let z = rectangular_plastic_section_modulus(0.1, 0.2);
        assert!((z - 0.1 * 0.2_f64 * 0.2 / 4.0).abs() <= 1e-9 * z, "Z = b·h²/4");

        // (b) THREAD elastic_section_modulus + rectangular_second_moment_of_area (#415)
        // (non-tautological): the rectangular elastic modulus is S = I/(h/2) = b·h²/6, and the
        // plastic modulus is exactly 1.5·S — the rectangle's shape factor Z/S = 1.5.
        let (b, h) = (0.1_f64, 0.2_f64);
        let s = elastic_section_modulus(rectangular_second_moment_of_area(b, h), h / 2.0);
        assert!(
            (rectangular_plastic_section_modulus(b, h) - 1.5 * s).abs()
                <= 1e-9 * rectangular_plastic_section_modulus(b, h),
            "Z = 1.5·S"
        );
        assert!(
            (rectangular_plastic_section_modulus(b, h) / s - 1.5).abs() <= 1e-9,
            "shape factor Z/S = 1.5"
        );

        // (c) SCALING: linear in width, quadratic in height.
        let base = rectangular_plastic_section_modulus(0.1, 0.2);
        assert!(
            (rectangular_plastic_section_modulus(0.2, 0.2) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in width"
        );
        assert!(
            (rectangular_plastic_section_modulus(0.1, 0.4) - 4.0 * base).abs() <= 1e-9 * 4.0 * base,
            "quadratic in height"
        );

        // (d) GUARD: non-positive or non-finite → 0.
        assert_eq!(rectangular_plastic_section_modulus(0.0, 0.2), 0.0);
        assert_eq!(rectangular_plastic_section_modulus(0.1, 0.0), 0.0);
        assert_eq!(rectangular_plastic_section_modulus(f64::NAN, 0.2), 0.0);
        assert_eq!(rectangular_plastic_section_modulus(0.1, -0.2), 0.0);
    }

    #[test]
    fn rectangular_polar_second_moment_of_area_is_ix_plus_iy() {
        // (a) WORKED: b=0.1, h=0.2 → J = (b·h/12)(b²+h²) = (0.02/12)(0.05) ≈ 8.3333e-5 m⁴.
        let j = rectangular_polar_second_moment_of_area(0.1, 0.2);
        assert!(
            (j - 0.1 * 0.2_f64 * (0.1 * 0.1 + 0.2 * 0.2) / 12.0).abs() <= 1e-9 * j,
            "J = (b·h/12)(b²+h²)"
        );

        // (b) THREAD rectangular_second_moment_of_area (#415) (non-tautological, perpendicular-
        // axis theorem): J = I_x + I_y = I(b,h) + I(h,b).
        let (b, h) = (0.1_f64, 0.2_f64);
        assert!(
            (rectangular_polar_second_moment_of_area(b, h)
                - (rectangular_second_moment_of_area(b, h)
                    + rectangular_second_moment_of_area(h, b)))
            .abs()
                <= 1e-9 * rectangular_polar_second_moment_of_area(b, h),
            "J = I_x + I_y"
        );

        // (c) SQUARE: for b = h = a the polar moment is a⁴/6.
        assert!(
            (rectangular_polar_second_moment_of_area(0.1, 0.1) - 0.1_f64.powi(4) / 6.0).abs()
                <= 1e-9 * rectangular_polar_second_moment_of_area(0.1, 0.1),
            "square → a⁴/6"
        );

        // (d) SYMMETRY: the polar moment is symmetric in the two dimensions, J(b,h) = J(h,b).
        assert!(
            (rectangular_polar_second_moment_of_area(0.1, 0.2)
                - rectangular_polar_second_moment_of_area(0.2, 0.1))
            .abs()
                <= 1e-9 * rectangular_polar_second_moment_of_area(0.1, 0.2),
            "J(b,h) = J(h,b)"
        );

        // (e) GUARD: non-positive or non-finite → 0.
        assert_eq!(rectangular_polar_second_moment_of_area(0.0, 0.2), 0.0);
        assert_eq!(rectangular_polar_second_moment_of_area(0.1, 0.0), 0.0);
        assert_eq!(rectangular_polar_second_moment_of_area(f64::NAN, 0.2), 0.0);
        assert_eq!(rectangular_polar_second_moment_of_area(0.1, -0.2), 0.0);
    }

    #[test]
    fn circular_second_moment_of_area_is_pi_d4_over_64() {
        // (a) WORKED: d = 0.1 → I = π·d⁴/64 ≈ 4.9087e-6 m⁴.
        let i = circular_second_moment_of_area(0.1);
        assert!(
            (i - std::f64::consts::PI * 0.1_f64.powi(4) / 64.0).abs() <= 1e-9 * i,
            "I = π·d⁴/64"
        );

        // (b) THREAD elastic_section_modulus (#380) (non-tautological): the circular section
        // modulus S = I/(d/2) = π·d³/32.
        let d = 0.1_f64;
        assert!(
            (elastic_section_modulus(circular_second_moment_of_area(d), d / 2.0)
                - std::f64::consts::PI * d.powi(3) / 32.0)
            .abs()
                <= 1e-9 * (std::f64::consts::PI * d.powi(3) / 32.0),
            "S = I/(d/2) = π·d³/32"
        );

        // (c) THREAD flexural_rigidity (#398): EI of a steel round bar.
        let ei = 200.0e9 * circular_second_moment_of_area(0.1);
        assert!(
            (flexural_rigidity(200.0e9, circular_second_moment_of_area(0.1)) - ei).abs()
                <= 1e-6 * ei,
            "EI = E·I"
        );

        // (d) SCALING: quartic in diameter.
        let base = circular_second_moment_of_area(0.1);
        assert!(
            (circular_second_moment_of_area(0.2) - 16.0 * base).abs() <= 1e-9 * 16.0 * base,
            "quartic in diameter"
        );

        // (e) GUARD: non-positive or non-finite → 0.
        assert_eq!(circular_second_moment_of_area(0.0), 0.0);
        assert_eq!(circular_second_moment_of_area(f64::NAN), 0.0);
        assert_eq!(circular_second_moment_of_area(-0.1), 0.0);
    }

    #[test]
    fn hollow_circular_second_moment_of_area_is_pi_d4_minus_d4_over_64() {
        // (a) WORKED: D=0.1, d=0.05 → I = π(D⁴−d⁴)/64 ≈ 4.6019e-6 m⁴.
        let i = hollow_circular_second_moment_of_area(0.1, 0.05);
        assert!(
            (i - std::f64::consts::PI * (0.1_f64.powi(4) - 0.05_f64.powi(4)) / 64.0).abs()
                <= 1e-9 * i,
            "I = π(D⁴−d⁴)/64"
        );

        // (b) THREAD circular_second_moment_of_area (#419) (non-tautological): the annulus is
        // the outer disc minus the inner disc.
        assert!(
            (hollow_circular_second_moment_of_area(0.1, 0.05)
                - (circular_second_moment_of_area(0.1) - circular_second_moment_of_area(0.05)))
            .abs()
                <= 1e-9 * hollow_circular_second_moment_of_area(0.1, 0.05),
            "I = I(D) − I(d)"
        );

        // (c) SOLID LIMIT: a zero bore is a solid circle.
        assert!(
            (hollow_circular_second_moment_of_area(0.1, 0.0) - circular_second_moment_of_area(0.1))
                .abs()
                <= 1e-9 * circular_second_moment_of_area(0.1),
            "zero bore → solid circle"
        );

        // (d) MONOTONICITY: a larger bore removes more material, so I drops.
        assert!(
            hollow_circular_second_moment_of_area(0.1, 0.08)
                < hollow_circular_second_moment_of_area(0.1, 0.05),
            "larger bore → smaller I"
        );

        // (e) GUARD: non-physical input → 0.
        assert_eq!(hollow_circular_second_moment_of_area(0.0, 0.0), 0.0);
        assert_eq!(hollow_circular_second_moment_of_area(0.1, 0.1), 0.0); // d ≥ D
        assert_eq!(hollow_circular_second_moment_of_area(0.1, 0.2), 0.0); // d > D
        assert_eq!(hollow_circular_second_moment_of_area(f64::NAN, 0.05), 0.0);
        assert_eq!(hollow_circular_second_moment_of_area(0.1, -0.01), 0.0);
    }

    #[test]
    fn hollow_circular_polar_second_moment_of_area_is_pi_d4_minus_d4_over_32() {
        // (a) WORKED: D=0.1, d=0.05 → J = π(D⁴−d⁴)/32 ≈ 9.2038e-6 m⁴.
        let j = hollow_circular_polar_second_moment_of_area(0.1, 0.05);
        assert!(
            (j - std::f64::consts::PI * (0.1_f64.powi(4) - 0.05_f64.powi(4)) / 32.0).abs()
                <= 1e-9 * j,
            "J = π(D⁴−d⁴)/32"
        );

        // (b) THREAD hollow_circular_second_moment_of_area (#435) (non-tautological,
        // perpendicular-axis theorem): J = 2·I_hollow.
        assert!(
            (hollow_circular_polar_second_moment_of_area(0.1, 0.05)
                - 2.0 * hollow_circular_second_moment_of_area(0.1, 0.05))
            .abs()
                <= 1e-9 * hollow_circular_polar_second_moment_of_area(0.1, 0.05),
            "J = 2·I_hollow"
        );

        // (c) THREAD circular_polar_second_moment_of_area (#423) (non-tautological, annulus):
        // J = J(D) − J(d).
        assert!(
            (hollow_circular_polar_second_moment_of_area(0.1, 0.05)
                - (circular_polar_second_moment_of_area(0.1)
                    - circular_polar_second_moment_of_area(0.05)))
            .abs()
                <= 1e-9 * hollow_circular_polar_second_moment_of_area(0.1, 0.05),
            "J = J(D) − J(d)"
        );

        // (d) SOLID LIMIT: a zero bore is a solid shaft.
        assert!(
            (hollow_circular_polar_second_moment_of_area(0.1, 0.0)
                - circular_polar_second_moment_of_area(0.1))
            .abs()
                <= 1e-9 * circular_polar_second_moment_of_area(0.1),
            "zero bore → solid shaft"
        );

        // (e) GUARD: non-physical input → 0.
        assert_eq!(hollow_circular_polar_second_moment_of_area(0.0, 0.0), 0.0);
        assert_eq!(hollow_circular_polar_second_moment_of_area(0.1, 0.1), 0.0); // d ≥ D
        assert_eq!(hollow_circular_polar_second_moment_of_area(0.1, 0.2), 0.0); // d > D
        assert_eq!(hollow_circular_polar_second_moment_of_area(f64::NAN, 0.05), 0.0);
        assert_eq!(hollow_circular_polar_second_moment_of_area(0.1, -0.01), 0.0);
    }

    #[test]
    fn hollow_rectangular_second_moment_of_area_is_outer_minus_inner() {
        // WORKED: 2×2 outer, 1×1 inner → I = (2·2³ − 1·1³)/12 = (16 − 1)/12 = 1.25 m⁴.
        assert!(
            (hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.0, 1.0) - 1.25).abs() < 1e-12,
            "I = (bh³ − bᵢhᵢ³)/12 = 1.25",
        );

        // STRONG non-tautological thread: the box is the solid outer minus the solid inner,
        // threading the existing rectangular_second_moment_of_area.
        for &(b, h, bi, hi) in &[(0.1_f64, 0.2_f64, 0.06, 0.12), (3.0, 1.0, 2.0, 0.5)] {
            let hollow = hollow_rectangular_second_moment_of_area(b, h, bi, hi);
            let outer = rectangular_second_moment_of_area(b, h);
            let inner = rectangular_second_moment_of_area(bi, hi);
            assert!(
                (hollow - (outer - inner)).abs() <= 1e-12 * outer,
                "hollow = I_outer − I_inner",
            );
        }

        // A larger bore removes more material → smaller I.
        assert!(
            hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.5, 1.5)
                < hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.0, 1.0),
            "larger bore → smaller I",
        );
        // Cubic in outer depth: doubling h (with proportional bore) scales I by 8.
        let base = hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.0, 1.0);
        assert!(
            (hollow_rectangular_second_moment_of_area(2.0, 4.0, 1.0, 2.0) - 8.0 * base).abs()
                <= 1e-9 * (8.0 * base),
            "∝ h³ (with proportional bore)",
        );

        // Guards: non-physical / bore-not-inside → 0.
        assert_eq!(hollow_rectangular_second_moment_of_area(0.0, 2.0, 1.0, 1.0), 0.0);
        assert_eq!(hollow_rectangular_second_moment_of_area(2.0, 2.0, 2.0, 1.0), 0.0); // bᵢ ≥ b
        assert_eq!(hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.0, 2.0), 0.0); // hᵢ ≥ h
        assert_eq!(hollow_rectangular_second_moment_of_area(f64::NAN, 2.0, 1.0, 1.0), 0.0);
        assert_eq!(hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.0, -0.1), 0.0);
    }

    #[test]
    fn hollow_rectangular_polar_second_moment_of_area_is_outer_minus_inner() {
        // WORKED: 2×2 outer, 1×1 bore → J = (b·h·(b²+h²) − bᵢ·hᵢ·(bᵢ²+hᵢ²))/12
        //         = (2·2·8 − 1·1·2)/12 = (32 − 2)/12 = 2.5 m⁴.
        assert!(
            (hollow_rectangular_polar_second_moment_of_area(2.0, 2.0, 1.0, 1.0) - 2.5).abs() < 1e-12,
            "J = 2.5",
        );

        // STRONG non-tautological thread: J_hollow = J_outer − J_inner via the solid polar fn.
        for &(b, h, bi, hi) in &[(0.1_f64, 0.2_f64, 0.06, 0.12), (3.0, 1.0, 2.0, 0.5)] {
            let hollow = hollow_rectangular_polar_second_moment_of_area(b, h, bi, hi);
            let outer = rectangular_polar_second_moment_of_area(b, h);
            let inner = rectangular_polar_second_moment_of_area(bi, hi);
            assert!(
                (hollow - (outer - inner)).abs() <= 1e-12 * outer,
                "J = J_outer − J_inner",
            );
        }

        // For a SQUARE box (b=h, bᵢ=hᵢ) the polar J = Ix + Iy = 2·(bending I).
        let j_square = hollow_rectangular_polar_second_moment_of_area(2.0, 2.0, 1.0, 1.0);
        let i_square = hollow_rectangular_second_moment_of_area(2.0, 2.0, 1.0, 1.0);
        assert!(
            (j_square - 2.0 * i_square).abs() <= 1e-12 * j_square,
            "square box: J = 2·I",
        );

        // Guards: non-physical / bore-not-inside → 0.
        assert_eq!(hollow_rectangular_polar_second_moment_of_area(0.0, 2.0, 1.0, 1.0), 0.0);
        assert_eq!(hollow_rectangular_polar_second_moment_of_area(2.0, 2.0, 2.0, 1.0), 0.0); // bᵢ ≥ b
        assert_eq!(hollow_rectangular_polar_second_moment_of_area(2.0, 2.0, 1.0, 2.0), 0.0); // hᵢ ≥ h
        assert_eq!(hollow_rectangular_polar_second_moment_of_area(f64::NAN, 2.0, 1.0, 1.0), 0.0);
        assert_eq!(hollow_rectangular_polar_second_moment_of_area(2.0, 2.0, 1.0, -0.1), 0.0);
    }

    #[test]
    fn circular_plastic_section_modulus_is_d3_over_6() {
        // (a) WORKED: d = 0.1 → Z = d³/6 = 0.001/6 ≈ 1.6667e-4 m³.
        let z = circular_plastic_section_modulus(0.1);
        assert!((z - 0.1_f64.powi(3) / 6.0).abs() <= 1e-9 * z, "Z = d³/6");

        // (b) THREAD elastic_section_modulus + circular_second_moment_of_area (#419)
        // (non-tautological): the circular elastic modulus is S = I/(d/2) = π·d³/32, and the
        // plastic modulus is exactly (16/(3π))·S — the solid-circle shape factor.
        let d = 0.1_f64;
        let s = elastic_section_modulus(circular_second_moment_of_area(d), d / 2.0);
        let shape = 16.0 / (3.0 * std::f64::consts::PI);
        assert!(
            (circular_plastic_section_modulus(d) - shape * s).abs()
                <= 1e-9 * circular_plastic_section_modulus(d),
            "Z = (16/3π)·S"
        );
        assert!(
            (circular_plastic_section_modulus(d) / s - shape).abs() <= 1e-9 * shape,
            "shape factor Z/S = 16/(3π)"
        );

        // (c) SCALING: cubic in diameter.
        let base = circular_plastic_section_modulus(0.1);
        assert!(
            (circular_plastic_section_modulus(0.2) - 8.0 * base).abs() <= 1e-9 * 8.0 * base,
            "cubic in diameter"
        );

        // (d) GUARD: non-positive or non-finite → 0.
        assert_eq!(circular_plastic_section_modulus(0.0), 0.0);
        assert_eq!(circular_plastic_section_modulus(f64::NAN), 0.0);
        assert_eq!(circular_plastic_section_modulus(-0.1), 0.0);
    }

    #[test]
    fn circular_polar_second_moment_of_area_is_pi_d4_over_32() {
        // (a) WORKED: d = 0.1 → J = π·d⁴/32 ≈ 9.8175e-6 m⁴.
        let j = circular_polar_second_moment_of_area(0.1);
        assert!(
            (j - std::f64::consts::PI * 0.1_f64.powi(4) / 32.0).abs() <= 1e-9 * j,
            "J = π·d⁴/32"
        );

        // (b) THREAD circular_second_moment_of_area (#419) (non-tautological, perpendicular-axis
        // theorem): J = I_x + I_y = 2·I for a circle.
        assert!(
            (circular_polar_second_moment_of_area(0.1) - 2.0 * circular_second_moment_of_area(0.1))
                .abs()
                <= 1e-9 * j,
            "J = 2·I (perpendicular-axis theorem)"
        );

        // (c) THREAD polar_section_modulus: the polar section modulus Z_p = J/(d/2) = π·d³/16.
        let d = 0.1_f64;
        assert!(
            (polar_section_modulus(circular_polar_second_moment_of_area(d), d / 2.0)
                - std::f64::consts::PI * d.powi(3) / 16.0)
            .abs()
                <= 1e-9 * (std::f64::consts::PI * d.powi(3) / 16.0),
            "Z_p = J/(d/2) = π·d³/16"
        );

        // (d) THREAD torsional_rigidity (#410): GJ of a steel shaft.
        let gj = 80.0e9 * circular_polar_second_moment_of_area(0.1);
        assert!(
            (torsional_rigidity(80.0e9, circular_polar_second_moment_of_area(0.1)) - gj).abs()
                <= 1e-6 * gj,
            "GJ = G·J"
        );

        // (e) SCALING: quartic in diameter.
        let base = circular_polar_second_moment_of_area(0.1);
        assert!(
            (circular_polar_second_moment_of_area(0.2) - 16.0 * base).abs() <= 1e-9 * 16.0 * base,
            "quartic in diameter"
        );

        // (f) GUARD: non-positive or non-finite → 0.
        assert_eq!(circular_polar_second_moment_of_area(0.0), 0.0);
        assert_eq!(circular_polar_second_moment_of_area(f64::NAN), 0.0);
        assert_eq!(circular_polar_second_moment_of_area(-0.1), 0.0);
    }

    #[test]
    fn beam_transverse_shear_stress_is_jourawski() {
        // Rectangular section b×h, transverse shear V: at the neutral axis Q = b·h²/8,
        // I = b·h³/12, so τ_max = 1.5·V/A (the classic 3/2 factor over the average V/A).
        let (b, h) = (0.05_f64, 0.10_f64);
        let q = b * h * h / 8.0;
        let i = b * h.powi(3) / 12.0;
        let v = 10000.0;
        let tau = beam_transverse_shear_stress(v, q, i, b);
        let area = b * h;
        assert!((tau - 1.5 * v / area).abs() <= 1e-9 * tau, "τ_max = 1.5·V/A");
        assert!(tau > v / area, "peak exceeds the average shear V/A");

        // Linear in V and Q; inverse in I and b.
        assert!(
            (beam_transverse_shear_stress(2.0 * v, q, i, b) - 2.0 * tau).abs() <= 1e-9 * (2.0 * tau),
            "∝ V"
        );
        assert!(
            (beam_transverse_shear_stress(v, 3.0 * q, i, b) - 3.0 * tau).abs() <= 1e-9 * (3.0 * tau),
            "∝ Q"
        );
        assert!(
            (beam_transverse_shear_stress(v, q, 2.0 * i, b) - 0.5 * tau).abs() <= 1e-9 * (0.5 * tau),
            "∝ 1/I"
        );
        assert!(
            (beam_transverse_shear_stress(v, q, i, 2.0 * b) - 0.5 * tau).abs() <= 1e-9 * (0.5 * tau),
            "∝ 1/b"
        );

        // 0-sentinel for non-physical input.
        assert_eq!(beam_transverse_shear_stress(v, q, 0.0, b), 0.0);
        assert_eq!(beam_transverse_shear_stress(v, q, i, 0.0), 0.0);
        assert_eq!(beam_transverse_shear_stress(f64::NAN, q, i, b), 0.0);
    }

    #[test]
    fn bending_stress_is_the_flexure_formula() {
        // Worked: σ = M·y/I = 2000·0.05/1e-6 = 100 MPa.
        let s = bending_stress(2000.0, 0.05, 1.0e-6);
        assert!((s - 1.0e8).abs() <= 1e-12 * 1.0e8, "σ = M·y/I, got {s}");

        // Threads beam_curvature via Hooke's law σ = E·κ·y (the E cancels).
        let e = 200.0e9;
        for &(m, y, i) in &[
            (2000.0_f64, 0.05_f64, 1.0e-6_f64),
            (-450.0, -0.02, 4.2e-7),
            (8200.0, 0.12, 9.0e-8),
        ] {
            let from_curvature = e * beam_curvature(m, e, i) * y;
            assert!(
                (bending_stress(m, y, i) - from_curvature).abs() <= 1e-12 * from_curvature.abs(),
                "σ = E·κ·y"
            );
        }

        // Threads the cantilever moment family: σ_root = M_root·c/I = P·L·c/I.
        for &(p, l, c, i) in &[(1000.0_f64, 2.0_f64, 0.05_f64, 1.0e-6_f64), (-450.0, 3.5, 0.03, 4.2e-7)]
        {
            let sigma_root = bending_stress(cantilever_point_load_root_moment(p, l), c, i);
            assert!(
                (sigma_root - p * l * c / i).abs() <= 1e-12 * (p * l * c / i).abs(),
                "σ_root = PLc/I"
            );
        }

        // Opposite faces carry opposite (tension/compression) stress; inverse in I.
        assert!(
            (bending_stress(2000.0, -0.05, 1.0e-6) + bending_stress(2000.0, 0.05, 1.0e-6)).abs() < 1e-3,
            "σ(−y) = −σ(y)"
        );
        assert!(
            (bending_stress(2000.0, 0.05, 2.0e-6) - 0.5 * bending_stress(2000.0, 0.05, 1.0e-6)).abs()
                < 1e-3,
            "inverse in I"
        );

        // Non-physical input → 0.
        assert_eq!(bending_stress(f64::NAN, 0.05, 1.0e-6), 0.0);
        assert_eq!(bending_stress(2000.0, f64::NAN, 1.0e-6), 0.0);
        assert_eq!(bending_stress(2000.0, 0.05, 0.0), 0.0);
        assert_eq!(bending_stress(2000.0, 0.05, -1.0e-6), 0.0);
    }

    #[test]
    fn bending_moment_capacity_inverts_the_flexure_formula() {
        // (a) WORKED: M = σ·S = 1.0e8·2.0e-5 = 2000 N·m.
        assert!(
            (bending_moment_capacity(1.0e8, 2.0e-5) - 2000.0).abs() <= 1e-9 * 2000.0,
            "M = σ·S = 2000 N·m"
        );

        // (b) ROUND-TRIP threading bending_stress + elastic_section_modulus
        // (non-tautological): σ·S = (M/S)·S = M at the extreme fibre y = c.
        for &(m, c, i) in &[(2000.0_f64, 0.05_f64, 1.0e-6_f64), (8200.0, 0.12, 9.0e-8)] {
            let s = elastic_section_modulus(i, c);
            let sigma = bending_stress(m, c, i);
            assert!(
                (bending_moment_capacity(sigma, s) - m).abs() <= 1e-9 * m,
                "σ·S recovers M for (M,c,I)=({m},{c},{i})"
            );
        }

        // (c) FIRST-YIELD round-trip threading all three: a section's first-yield
        // moment M_y = σ_y·S maps back to exactly σ_y at the extreme fibre.
        let (sy, c, i) = (250.0e6_f64, 0.05_f64, 1.0e-6_f64);
        let s = elastic_section_modulus(i, c);
        let my = bending_moment_capacity(sy, s);
        assert!(
            (bending_stress(my, c, i) - sy).abs() <= 1e-9 * sy,
            "M_y = σ_y·S maps back to σ_y"
        );

        // (d) LINEARITY + SIGN: linear in stress and in section modulus; a negative
        // stress yields a negative (reversed) moment.
        let base = bending_moment_capacity(1.0e8, 2.0e-5);
        assert!(
            (bending_moment_capacity(2.0e8, 2.0e-5) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in stress"
        );
        assert!(
            (bending_moment_capacity(1.0e8, 4.0e-5) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in section modulus"
        );
        assert!(
            bending_moment_capacity(-1.0e8, 2.0e-5) < 0.0,
            "negative stress → negative moment"
        );

        // (e) GUARD: non-positive section modulus or non-finite input → 0 sentinel.
        assert_eq!(bending_moment_capacity(1.0e8, 0.0), 0.0);
        assert_eq!(bending_moment_capacity(1.0e8, -1.0e-5), 0.0);
        assert_eq!(bending_moment_capacity(f64::NAN, 2.0e-5), 0.0);
        assert_eq!(bending_moment_capacity(1.0e8, f64::NAN), 0.0);
    }

    #[test]
    fn beam_curvature_is_the_moment_curvature_law() {
        // Worked: κ = M/(EI) = 2000/(200e9·1e-6) = 0.01 1/m.
        let k = beam_curvature(2000.0, 200.0e9, 1.0e-6);
        assert!((k - 0.01).abs() < 1e-12, "κ = M/EI, got {k}");

        // Moment–curvature inverse: M = E·I·κ.
        for &(m, e, i) in &[
            (2000.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 70.0e9, 4.2e-7),
            (8200.0, 120.0e9, 9.0e-8),
        ] {
            assert!((beam_curvature(m, e, i) * e * i - m).abs() <= 1e-12 * m.abs(), "M = E·I·κ");
        }

        // Threads the cantilever family: δ_tip = κ_root·L²/3 (since δ = PL³/3EI and
        // κ_root = M_root/EI = PL/EI).
        for &(p, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 1.2, 120.0e9, 9.0e-8),
        ] {
            let from_curvature =
                beam_curvature(cantilever_point_load_root_moment(p, l), e, i) * l * l / 3.0;
            let delta = cantilever_tip_deflection(p, l, e, i);
            assert!((from_curvature - delta).abs() <= 1e-12 * delta.abs(), "δ_tip = κ_root·L²/3");
        }

        // Linear & sign-preserving in M; inverse in E·I.
        assert!(beam_curvature(-2000.0, 200.0e9, 1.0e-6) < 0.0, "sign follows the moment");
        assert!(
            (beam_curvature(2000.0, 400.0e9, 2.0e-6)
                - 0.25 * beam_curvature(2000.0, 200.0e9, 1.0e-6))
            .abs()
                < 1e-12,
            "inverse in E·I"
        );

        // Non-physical input → 0.
        assert_eq!(beam_curvature(f64::NAN, 200.0e9, 1.0e-6), 0.0);
        assert_eq!(beam_curvature(2000.0, 0.0, 1.0e-6), 0.0);
        assert_eq!(beam_curvature(2000.0, 200.0e9, -1.0e-6), 0.0);
    }

    #[test]
    fn flexural_rigidity_is_the_bending_stiffness() {
        // (a) WORKED: EI = E·I = 200e9·1e-6 = 2e5 N·m².
        assert!(
            (flexural_rigidity(200.0e9, 1.0e-6) - 2.0e5).abs() <= 1e-9 * 2.0e5,
            "EI = E·I = 2e5 N·m²"
        );

        // (b) THREAD beam_curvature (non-tautological): EI = M/κ.
        for &(m, e, i) in &[(2000.0_f64, 200.0e9_f64, 1.0e-6_f64), (8200.0, 70.0e9, 9.0e-8)] {
            assert!(
                (flexural_rigidity(e, i) - m / beam_curvature(m, e, i)).abs()
                    <= 1e-9 * flexural_rigidity(e, i),
                "EI = M/κ"
            );
        }

        // (c) THREAD cantilever_tip_deflection (non-tautological): EI = P·L³/(3·δ).
        let (p, l, e, i) = (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64);
        assert!(
            (flexural_rigidity(e, i)
                - p * l * l * l / (3.0 * cantilever_tip_deflection(p, l, e, i)))
            .abs()
                <= 1e-9 * flexural_rigidity(e, i),
            "EI = PL³/(3δ)"
        );

        // (d) LINEARITY: linear in both E and I.
        let base = flexural_rigidity(200.0e9, 1.0e-6);
        assert!(
            (flexural_rigidity(2.0 * 200.0e9, 1.0e-6) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in E"
        );
        assert!(
            (flexural_rigidity(200.0e9, 2.0e-6) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in I"
        );

        // (e) GUARD: non-positive E or I, or non-finite → 0 sentinel.
        assert_eq!(flexural_rigidity(0.0, 1.0e-6), 0.0);
        assert_eq!(flexural_rigidity(200.0e9, 0.0), 0.0);
        assert_eq!(flexural_rigidity(f64::NAN, 1.0e-6), 0.0);
        assert_eq!(flexural_rigidity(200.0e9, -1.0e-6), 0.0);
    }

    #[test]
    fn bulk_modulus_is_the_volumetric_stiffness() {
        // (a) WORKED (exact): ν=0.25 → 1−2ν=0.5, K = E/(3·0.5) = E/1.5; E=3e9 → 2e9.
        assert!(
            (bulk_modulus(3.0e9, 0.25) - 2.0e9).abs() <= 1e-9 * 2.0e9,
            "K = E/(3(1−2ν)) = 2 GPa"
        );

        // (b) INVERSE round-trip (non-tautological): E = 3·K·(1−2ν) recovers Young's modulus.
        for &(e, nu) in &[(200.0e9_f64, 0.3_f64), (70.0e9, 0.33), (110.0e9, 0.21)] {
            let k = bulk_modulus(e, nu);
            assert!(
                (3.0 * k * (1.0 - 2.0 * nu) - e).abs() <= 1e-9 * e,
                "E = 3K(1−2ν)"
            );
        }

        // (c) ν=0 limit: no lateral coupling → K = E/3.
        assert!(
            (bulk_modulus(200.0e9, 0.0) - 200.0e9 / 3.0).abs() <= 1e-9 * (200.0e9 / 3.0),
            "K(ν=0) = E/3"
        );

        // (d) INCOMPRESSIBLE trend: K → ∞ as ν → 0.5.
        assert!(
            bulk_modulus(200.0e9, 0.49) > 10.0 * bulk_modulus(200.0e9, 0.0),
            "K grows without bound toward ν=0.5"
        );

        // (e) GUARD: non-physical input → 0 sentinel.
        assert_eq!(bulk_modulus(0.0, 0.3), 0.0);
        assert_eq!(bulk_modulus(f64::NAN, 0.3), 0.0);
        assert_eq!(bulk_modulus(200.0e9, 0.5), 0.0);
        assert_eq!(bulk_modulus(200.0e9, -1.0), 0.0);
        assert_eq!(bulk_modulus(200.0e9, f64::NAN), 0.0);
    }

    #[test]
    fn lames_first_parameter_couples_volumetric_strain() {
        // (a) WORKED: E=200e9, ν=0.3 → λ = E·ν/((1+ν)(1−2ν)) = 60e9/0.52 ≈ 1.153846e11.
        let lam = lames_first_parameter(200.0e9, 0.3);
        assert!(
            (lam - 200.0e9 * 0.3 / (1.3 * 0.4)).abs() <= 1e-9 * lam,
            "λ = Eν/((1+ν)(1−2ν))"
        );

        // (b) THREAD bulk_modulus (non-tautological): λ = 3·K·ν/(1+ν).
        for &(e, nu) in &[(200.0e9_f64, 0.3_f64), (70.0e9, 0.33)] {
            let lam = lames_first_parameter(e, nu);
            assert!(
                (lam - 3.0 * bulk_modulus(e, nu) * nu / (1.0 + nu)).abs() <= 1e-9 * lam,
                "λ = 3Kν/(1+ν)"
            );
        }

        // (c) ν=0 → λ=0 (a valid zero, not a sentinel): no volumetric coupling.
        assert_eq!(lames_first_parameter(200.0e9, 0.0), 0.0);

        // (d) AUXETIC ν<0 → λ<0: a negative Poisson's ratio flips the coupling sign.
        assert!(lames_first_parameter(200.0e9, -0.2) < 0.0);

        // (e) GUARD: non-physical input → 0 sentinel.
        assert_eq!(lames_first_parameter(0.0, 0.3), 0.0);
        assert_eq!(lames_first_parameter(f64::NAN, 0.3), 0.0);
        assert_eq!(lames_first_parameter(200.0e9, 0.5), 0.0);
        assert_eq!(lames_first_parameter(200.0e9, -1.0), 0.0);
        assert_eq!(lames_first_parameter(200.0e9, f64::NAN), 0.0);
    }

    #[test]
    fn shear_modulus_from_youngs_completes_the_elastic_constant_set() {
        // (a) WORKED: E=200e9, ν=0.3 → G = E/(2(1+ν)) = 200e9/2.6 ≈ 7.6923e10.
        let g = shear_modulus_from_youngs(200.0e9, 0.3);
        assert!((g - 200.0e9 / 2.6).abs() <= 1e-9 * g, "G = E/(2(1+ν))");

        // (b) LAMÉ IDENTITY (thread bulk_modulus + lames_first_parameter, non-tautological):
        // K = λ + (2/3)G, so G = 1.5·(K − λ).
        for &(e, nu) in &[(200.0e9_f64, 0.3_f64), (70.0e9, 0.33)] {
            let g = shear_modulus_from_youngs(e, nu);
            assert!(
                (g - 1.5 * (bulk_modulus(e, nu) - lames_first_parameter(e, nu))).abs() <= 1e-9 * g,
                "G = 1.5(K − λ)"
            );
        }

        // (c) INTERCONVERSION (thread bulk_modulus): E = 9KG/(3K+G) recovers Young's modulus.
        let (e, nu) = (200.0e9_f64, 0.3_f64);
        let g = shear_modulus_from_youngs(e, nu);
        let k = bulk_modulus(e, nu);
        assert!((9.0 * k * g / (3.0 * k + g) - e).abs() <= 1e-9 * e, "E = 9KG/(3K+G)");

        // (d) ν=0 → G=E/2: no lateral coupling.
        assert!(
            (shear_modulus_from_youngs(200.0e9, 0.0) - 100.0e9).abs() <= 1e-9 * 100.0e9,
            "G(ν=0) = E/2"
        );

        // (e) GUARD: non-physical input → 0 sentinel.
        assert_eq!(shear_modulus_from_youngs(0.0, 0.3), 0.0);
        assert_eq!(shear_modulus_from_youngs(f64::NAN, 0.3), 0.0);
        assert_eq!(shear_modulus_from_youngs(200.0e9, 0.5), 0.0);
        assert_eq!(shear_modulus_from_youngs(200.0e9, -1.0), 0.0);
        assert_eq!(shear_modulus_from_youngs(200.0e9, f64::NAN), 0.0);
    }

    #[test]
    fn p_wave_modulus_is_the_capstone_of_the_elastic_constants() {
        // (a) WORKED: E=200e9, ν=0.3 → M = E(1−ν)/((1+ν)(1−2ν)) = 140e9/0.52 ≈ 2.6923e11.
        let m = p_wave_modulus(200.0e9, 0.3);
        assert!((m - 200.0e9 * 0.7 / (1.3 * 0.4)).abs() <= 1e-9 * m, "M = E(1−ν)/((1+ν)(1−2ν))");

        // (b) IDENTITY M = K + 4G/3 (thread bulk_modulus + shear_modulus_from_youngs, exact).
        for &(e, nu) in &[(200.0e9_f64, 0.3_f64), (70.0e9, 0.33)] {
            let m = p_wave_modulus(e, nu);
            assert!(
                (m - (bulk_modulus(e, nu) + 4.0 / 3.0 * shear_modulus_from_youngs(e, nu))).abs()
                    <= 1e-9 * m,
                "M = K + 4G/3"
            );
        }

        // (c) IDENTITY M = λ + 2G (thread lames_first_parameter + shear, exact — ties all
        // the constants together).
        let (e, nu) = (70.0e9_f64, 0.33_f64);
        let m = p_wave_modulus(e, nu);
        assert!(
            (m - (lames_first_parameter(e, nu) + 2.0 * shear_modulus_from_youngs(e, nu))).abs()
                <= 1e-9 * m,
            "M = λ + 2G"
        );

        // (d) ν=0 → M=E: no lateral constraint.
        assert!(
            (p_wave_modulus(200.0e9, 0.0) - 200.0e9).abs() <= 1e-9 * 200.0e9,
            "M(ν=0) = E"
        );

        // (e) GUARD: non-physical input → 0 sentinel.
        assert_eq!(p_wave_modulus(0.0, 0.3), 0.0);
        assert_eq!(p_wave_modulus(f64::NAN, 0.3), 0.0);
        assert_eq!(p_wave_modulus(200.0e9, 0.5), 0.0);
        assert_eq!(p_wave_modulus(200.0e9, -1.0), 0.0);
        assert_eq!(p_wave_modulus(200.0e9, f64::NAN), 0.0);
    }

    #[test]
    fn axial_rigidity_is_the_extensional_stiffness() {
        // (a) WORKED: EA = E·A = 200e9·1e-4 = 2e7 N.
        assert!(
            (axial_rigidity(200.0e9, 1.0e-4) - 2.0e7).abs() <= 1e-9 * 2.0e7,
            "EA = E·A = 2e7 N"
        );

        // (b) THREAD beam_axial_extension (non-tautological): EA = F·L/δ.
        let (f, l, e, a) = (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-4_f64);
        assert!(
            (axial_rigidity(e, a) - f * l / beam_axial_extension(f, l, e, a)).abs()
                <= 1e-9 * axial_rigidity(e, a),
            "EA = F·L/δ"
        );

        // (c) THREAD axial_strain_energy (non-tautological): EA = F²·L/(2·U).
        assert!(
            (axial_rigidity(e, a) - f * f * l / (2.0 * axial_strain_energy(f, l, e, a))).abs()
                <= 1e-9 * axial_rigidity(e, a),
            "EA = F²·L/(2U)"
        );

        // (d) LINEARITY: linear in both E and A.
        let base = axial_rigidity(200.0e9, 1.0e-4);
        assert!(
            (axial_rigidity(2.0 * 200.0e9, 1.0e-4) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in E"
        );
        assert!(
            (axial_rigidity(200.0e9, 2.0e-4) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in A"
        );

        // (e) GUARD: non-positive E or A, or non-finite → 0 sentinel.
        assert_eq!(axial_rigidity(0.0, 1.0e-4), 0.0);
        assert_eq!(axial_rigidity(200.0e9, 0.0), 0.0);
        assert_eq!(axial_rigidity(f64::NAN, 1.0e-4), 0.0);
        assert_eq!(axial_rigidity(200.0e9, -1.0e-4), 0.0);
    }

    #[test]
    fn torsional_rigidity_is_the_torsional_stiffness() {
        // (a) WORKED: GJ = G·J = 80e9·2e-8 = 1600 N·m².
        assert!(
            (torsional_rigidity(80.0e9, 2.0e-8) - 1600.0).abs() <= 1e-9 * 1600.0,
            "GJ = G·J = 1600 N·m²"
        );

        // (b) THREAD beam_angle_of_twist (non-tautological): GJ = T·L/φ.
        let (t, l, g, j) = (100.0_f64, 2.0_f64, 80.0e9_f64, 2.0e-8_f64);
        assert!(
            (torsional_rigidity(g, j) - t * l / beam_angle_of_twist(t, l, g, j)).abs()
                <= 1e-9 * torsional_rigidity(g, j),
            "GJ = T·L/φ"
        );

        // (c) THREAD torsional_strain_energy (non-tautological): GJ = T²·L/(2·U).
        assert!(
            (torsional_rigidity(g, j) - t * t * l / (2.0 * torsional_strain_energy(t, l, g, j)))
                .abs()
                <= 1e-9 * torsional_rigidity(g, j),
            "GJ = T²·L/(2U)"
        );

        // (d) LINEARITY: linear in both G and J.
        let base = torsional_rigidity(80.0e9, 2.0e-8);
        assert!(
            (torsional_rigidity(2.0 * 80.0e9, 2.0e-8) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in G"
        );
        assert!(
            (torsional_rigidity(80.0e9, 4.0e-8) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in J"
        );

        // (e) GUARD: non-positive G or J, or non-finite → 0 sentinel.
        assert_eq!(torsional_rigidity(0.0, 2.0e-8), 0.0);
        assert_eq!(torsional_rigidity(80.0e9, 0.0), 0.0);
        assert_eq!(torsional_rigidity(f64::NAN, 2.0e-8), 0.0);
        assert_eq!(torsional_rigidity(80.0e9, -2.0e-8), 0.0);
    }

    #[test]
    fn axial_stress_is_force_over_area() {
        // Threads beam_axial_extension via Hooke's law σ = E·δ/L (E and L cancel, since
        // δ = F·L/E·A).
        let (l, e) = (2.0, 200.0e9);
        for &(f, a) in &[(10000.0_f64, 1.0e-4_f64), (-5000.0, 3.0e-4), (25000.0, 5.0e-5)] {
            let from_strain = e * beam_axial_extension(f, l, e, a) / l;
            assert!(
                (axial_stress(f, a) - from_strain).abs() <= 1e-12 * from_strain.abs(),
                "σ = E·δ/L = F/A"
            );
        }

        // Worked: σ = F/A = 1e4 / 1e-4 = 100 MPa.
        assert!((axial_stress(10000.0, 1.0e-4) - 1.0e8).abs() <= 1e-12 * 1.0e8, "F/A = 100 MPa");

        // Sign-preserving (tension +, compression −), linear in F, inverse in A.
        assert_eq!(
            axial_stress(-10000.0, 1.0e-4),
            -axial_stress(10000.0, 1.0e-4),
            "compression is negative"
        );
        assert!(
            (axial_stress(20000.0, 1.0e-4) - 2.0 * axial_stress(10000.0, 1.0e-4)).abs() < 1e-3,
            "linear in F"
        );
        assert!(
            (axial_stress(10000.0, 2.0e-4) - 0.5 * axial_stress(10000.0, 1.0e-4)).abs() < 1e-3,
            "inverse in A"
        );

        // 0 sentinel for non-physical input.
        assert_eq!(axial_stress(10000.0, 0.0), 0.0);
        assert_eq!(axial_stress(10000.0, -1.0e-4), 0.0);
        assert_eq!(axial_stress(f64::NAN, 1.0e-4), 0.0);
    }

    #[test]
    fn axial_force_capacity_inverts_the_axial_stress() {
        // (a) WORKED: F = σ·A = 2.0e8·5.0e-4 = 1.0e5 N.
        assert!(
            (axial_force_capacity(2.0e8, 5.0e-4) - 1.0e5).abs() <= 1e-9 * 1.0e5,
            "F = σ·A = 100 kN"
        );

        // (b) ROUND-TRIP threading axial_stress (non-tautological): the force that
        // produces stress σ is recovered from σ and the area.
        for &(f, a) in &[(10000.0_f64, 1.0e-4_f64), (-5000.0, 3.0e-4)] {
            assert!(
                (axial_force_capacity(axial_stress(f, a), a) - f).abs() <= 1e-9 * f.abs().max(1.0),
                "σ·A recovers F"
            );
        }

        // (c) YIELD round-trip: a member's squash load F_y = σ_y·A maps back to σ_y.
        let (sy, a) = (250.0e6_f64, 5.0e-4_f64);
        let fy = axial_force_capacity(sy, a);
        assert!((axial_stress(fy, a) - sy).abs() <= 1e-9 * sy, "F_y = σ_y·A maps back to σ_y");

        // (d) LINEARITY + SIGN: linear in stress and area; a negative (compressive)
        // stress yields a negative force.
        let base = axial_force_capacity(2.0e8, 5.0e-4);
        assert!(
            (axial_force_capacity(2.0 * 2.0e8, 5.0e-4) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in stress"
        );
        assert!(
            (axial_force_capacity(2.0e8, 1.0e-3) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in area"
        );
        assert!(
            axial_force_capacity(-2.0e8, 5.0e-4) < 0.0,
            "compressive stress → negative force"
        );

        // (e) GUARD: non-positive area or non-finite input → 0 sentinel.
        assert_eq!(axial_force_capacity(2.0e8, 0.0), 0.0);
        assert_eq!(axial_force_capacity(2.0e8, -1.0e-4), 0.0);
        assert_eq!(axial_force_capacity(f64::NAN, 5.0e-4), 0.0);
        assert_eq!(axial_force_capacity(2.0e8, f64::NAN), 0.0);
    }

    #[test]
    fn axial_strain_energy_is_the_clapeyron_work() {
        // Threads beam_axial_extension via Clapeyron U = ½·F·δ (exact, incl. compression).
        for &(f, l, e, a) in &[
            (10000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-4_f64),
            (-5000.0, 1.5, 70.0e9, 3.0e-4),
            (25000.0, 4.0, 200.0e9, 5.0e-5),
        ] {
            let from_delta = 0.5 * f * beam_axial_extension(f, l, e, a);
            assert!(
                (axial_strain_energy(f, l, e, a) - from_delta).abs() <= 1e-12 * from_delta.abs(),
                "U = ½·F·δ"
            );
        }

        // Worked: U = F²L/2EA = 10000²·2/(2·200e9·1e-4) = 5 J.
        assert!(
            (axial_strain_energy(10000.0, 2.0, 200.0e9, 1.0e-4) - 5.0).abs() <= 1e-12 * 5.0,
            "U = F²L/2EA = 5 J"
        );

        // Energy is non-negative and quadratic in F (tension and compression store equal).
        assert_eq!(
            axial_strain_energy(-10000.0, 2.0, 200.0e9, 1.0e-4),
            axial_strain_energy(10000.0, 2.0, 200.0e9, 1.0e-4),
            "even in F"
        );
        assert!(
            (axial_strain_energy(20000.0, 2.0, 200.0e9, 1.0e-4)
                - 4.0 * axial_strain_energy(10000.0, 2.0, 200.0e9, 1.0e-4))
            .abs()
                < 1e-9,
            "quadratic in F"
        );
        assert!(axial_strain_energy(10000.0, 2.0, 200.0e9, 1.0e-4) > 0.0, "non-negative");

        // 0 sentinel for non-physical input.
        assert_eq!(axial_strain_energy(10000.0, 0.0, 200.0e9, 1.0e-4), 0.0);
        assert_eq!(axial_strain_energy(10000.0, 2.0, 0.0, 1.0e-4), 0.0);
        assert_eq!(axial_strain_energy(10000.0, 2.0, 200.0e9, 0.0), 0.0);
        assert_eq!(axial_strain_energy(f64::NAN, 2.0, 200.0e9, 1.0e-4), 0.0);
    }

    #[test]
    fn beam_axial_extension_matches_the_closed_form() {
        // Worked point: F = 10 kN on a 2 m steel bar, E = 200 GPa, A = 1e-4 m²
        // → δ = F·L/(E·A) = 20000 / 2e7 = 0.001 m = 1 mm.
        let (f, l, e, a) = (10_000.0, 2.0, 200.0e9, 1.0e-4);
        let delta = beam_axial_extension(f, l, e, a);
        assert!((delta - 0.001).abs() / delta < 1e-9, "δ = 1 mm, got {delta}");
        // Linear in the force, and sign-preserving (tension lengthens, compression shortens).
        assert!((beam_axial_extension(2.0 * f, l, e, a) - 2.0 * delta).abs() / delta < 1e-12);
        assert!(
            (beam_axial_extension(-f, l, e, a) + delta).abs() / delta < 1e-12,
            "compression shortens"
        );
        // Linear in the length: double L → double δ.
        assert!(
            (beam_axial_extension(f, 2.0 * l, e, a) - 2.0 * delta).abs() / delta < 1e-12,
            "L scaling"
        );
        // Inverse in the axial rigidity E·A: double E or A → half δ.
        assert!((beam_axial_extension(f, l, 2.0 * e, a) - 0.5 * delta).abs() / delta < 1e-12, "1/E");
        assert!((beam_axial_extension(f, l, e, 2.0 * a) - 0.5 * delta).abs() / delta < 1e-12, "1/A");
        // Non-physical input → 0.
        assert_eq!(beam_axial_extension(f, l, e, -1.0e-4), 0.0); // A ≤ 0
        assert_eq!(beam_axial_extension(f, l, 0.0, a), 0.0); // E ≤ 0
        assert_eq!(beam_axial_extension(f, -1.0, e, a), 0.0); // L ≤ 0
        assert_eq!(beam_axial_extension(f64::NAN, l, e, a), 0.0); // non-finite F
        assert_eq!(beam_axial_extension(f, l, f64::INFINITY, a), 0.0); // non-finite E
    }

    #[test]
    fn cantilever_torque_matches_analytic_twist() {
        // Twist of a shaft under an end torque: φ = T·L/(G·J).
        let mat = steel();
        let l = 2.0;
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(l, 0.0, 0.0),
        ];
        let section = BeamSection::circle(0.03);
        let elements = [BeamElement::new(0, 1, section)];
        let torque = 800.0;
        let constraints = [BeamConstraint::clamped(0)];
        // Moment about the beam axis (global X).
        let loads = [BeamLoad::moment(1, [torque, 0.0, 0.0])];
        let sol = solve_beam_static(&nodes, &elements, &mat, &constraints, &loads).unwrap();
        let g = mat.youngs_modulus / (2.0 * (1.0 + mat.poisson_ratio));
        let analytic = beam_angle_of_twist(torque, l, g, section.j);
        let rel = (sol.rotation[1][0] - analytic).abs() / analytic;
        assert!(rel < 1e-6, "twist {} vs analytic {analytic}", sol.rotation[1][0]);
    }

    #[test]
    fn polar_section_modulus_is_the_torsion_design_property() {
        // Worked: Z_p = J/r = 1e-6/0.05 = 2e-5 m³.
        assert!(
            (polar_section_modulus(1.0e-6, 0.05) - 2.0e-5).abs() <= 1e-12 * 2.0e-5,
            "Z_p = J/r"
        );

        // Threads torsional_shear_stress: τ_max = T / Z_p (the conjugate design relation).
        for &(tq, r, j) in &[
            (1000.0_f64, 0.05_f64, 1.0e-6_f64),
            (-450.0, 0.02, 4.2e-7),
            (8200.0, 0.12, 9.0e-8),
        ] {
            let from_zp = tq / polar_section_modulus(j, r);
            assert!(
                (torsional_shear_stress(tq, r, j) - from_zp).abs() <= 1e-12 * from_zp.abs(),
                "τ_max = T / Z_p"
            );
        }

        // Linear in J, inverse in r.
        assert!(
            (polar_section_modulus(2.0e-6, 0.05) - 2.0 * polar_section_modulus(1.0e-6, 0.05)).abs()
                < 1e-15,
            "linear in J"
        );
        assert!(
            (polar_section_modulus(1.0e-6, 0.10) - 0.5 * polar_section_modulus(1.0e-6, 0.05)).abs()
                < 1e-15,
            "inverse in r"
        );

        // Non-physical input → 0.
        assert_eq!(polar_section_modulus(0.0, 0.05), 0.0);
        assert_eq!(polar_section_modulus(1.0e-6, 0.0), 0.0);
        assert_eq!(polar_section_modulus(-1.0e-6, 0.05), 0.0);
        assert_eq!(polar_section_modulus(f64::NAN, 0.05), 0.0);
    }

    #[test]
    fn torsional_shear_stress_is_the_torsion_formula() {
        // Worked: τ = T·r/J = 1000·0.05/1e-6 = 50 MPa.
        let t = torsional_shear_stress(1000.0, 0.05, 1.0e-6);
        assert!((t - 5.0e7).abs() <= 1e-12 * 5.0e7, "τ = T·r/J, got {t}");

        // Threads beam_angle_of_twist via Hooke's law in shear τ = G·θ·r/L (G and L
        // cancel, since θ = T·L/(G·J)).
        let g = 80.0e9;
        for &(tq, r, j, l) in &[
            (1000.0_f64, 0.05_f64, 1.0e-6_f64, 2.0_f64),
            (-450.0, 0.02, 4.2e-7, 0.8),
            (8200.0, 0.12, 9.0e-8, 3.5),
        ] {
            let from_twist = g * beam_angle_of_twist(tq, l, g, j) * r / l;
            assert!(
                (torsional_shear_stress(tq, r, j) - from_twist).abs() <= 1e-12 * from_twist.abs(),
                "τ = G·θ·r/L"
            );
        }

        // Linear & sign-preserving in T and r, inverse in J.
        assert!(torsional_shear_stress(-1000.0, 0.05, 1.0e-6) < 0.0, "sign follows the torque");
        assert!(
            (torsional_shear_stress(1000.0, 0.05, 2.0e-6)
                - 0.5 * torsional_shear_stress(1000.0, 0.05, 1.0e-6))
            .abs()
                < 1e-3,
            "inverse in J"
        );

        // Non-physical input → 0.
        assert_eq!(torsional_shear_stress(f64::NAN, 0.05, 1.0e-6), 0.0);
        assert_eq!(torsional_shear_stress(1000.0, 0.05, 0.0), 0.0);
        assert_eq!(torsional_shear_stress(1000.0, 0.05, -1.0e-6), 0.0);
    }

    #[test]
    fn torsional_moment_capacity_inverts_the_torsion_formula() {
        // (a) WORKED: T = τ·Z_p = 5.0e7·2.0e-5 = 1000 N·m.
        assert!(
            (torsional_moment_capacity(5.0e7, 2.0e-5) - 1000.0).abs() <= 1e-9 * 1000.0,
            "T = τ·Z_p = 1000 N·m"
        );

        // (b) ROUND-TRIP threading torsional_shear_stress + polar_section_modulus
        // (non-tautological): at the outer radius τ = T·r/J = T/Z_p, so τ·Z_p = T.
        for &(t, r, j) in &[(1000.0_f64, 0.05_f64, 1.0e-6_f64), (450.0, 0.03, 4.2e-7)] {
            let zp = polar_section_modulus(j, r);
            let tau = torsional_shear_stress(t, r, j);
            assert!(
                (torsional_moment_capacity(tau, zp) - t).abs() <= 1e-9 * t,
                "τ·Z_p recovers T for (T,r,J)=({t},{r},{j})"
            );
        }

        // (c) FIRST-YIELD round-trip threading all three: a shaft's yield torque
        // T_y = τ_y·Z_p maps back to exactly τ_y at the outer surface.
        let (ty, r, j) = (120.0e6_f64, 0.05_f64, 1.0e-6_f64);
        let zp = polar_section_modulus(j, r);
        let ty_torque = torsional_moment_capacity(ty, zp);
        assert!(
            (torsional_shear_stress(ty_torque, r, j) - ty).abs() <= 1e-9 * ty,
            "T_y = τ_y·Z_p maps back to τ_y"
        );

        // (d) LINEARITY + SIGN: linear in stress and polar modulus; a negative shear
        // stress yields a negative (reversed) torque.
        let base = torsional_moment_capacity(5.0e7, 2.0e-5);
        assert!(
            (torsional_moment_capacity(2.0 * 5.0e7, 2.0e-5) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in shear stress"
        );
        assert!(
            (torsional_moment_capacity(5.0e7, 4.0e-5) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in polar modulus"
        );
        assert!(
            torsional_moment_capacity(-5.0e7, 2.0e-5) < 0.0,
            "negative shear stress → negative torque"
        );

        // (e) GUARD: non-positive polar modulus or non-finite input → 0 sentinel.
        assert_eq!(torsional_moment_capacity(5.0e7, 0.0), 0.0);
        assert_eq!(torsional_moment_capacity(5.0e7, -1.0e-5), 0.0);
        assert_eq!(torsional_moment_capacity(f64::NAN, 2.0e-5), 0.0);
        assert_eq!(torsional_moment_capacity(5.0e7, f64::NAN), 0.0);
    }

    #[test]
    fn torsional_strain_energy_is_the_clapeyron_work() {
        // Threads beam_angle_of_twist via Clapeyron U = ½·T·θ (exact, incl. negative T).
        for &(tq, l, g, j) in &[
            (100.0_f64, 2.0_f64, 80.0e9_f64, 1.0e-6_f64),
            (-450.0, 0.8, 27.0e9, 4.2e-7),
            (8200.0, 3.5, 80.0e9, 9.0e-8),
        ] {
            let from_twist = 0.5 * tq * beam_angle_of_twist(tq, l, g, j);
            assert!(
                (torsional_strain_energy(tq, l, g, j) - from_twist).abs() <= 1e-12 * from_twist.abs(),
                "U = ½·T·θ"
            );
        }

        // Worked: U = T²L/2GJ = 100²·2/(2·80e9·1e-6) = 0.125 J.
        assert!(
            (torsional_strain_energy(100.0, 2.0, 80.0e9, 1.0e-6) - 0.125).abs() <= 1e-12 * 0.125,
            "U = T²L/2GJ = 0.125 J"
        );

        // Energy is non-negative and quadratic in T (sign-independent).
        assert_eq!(
            torsional_strain_energy(-100.0, 2.0, 80.0e9, 1.0e-6),
            torsional_strain_energy(100.0, 2.0, 80.0e9, 1.0e-6),
            "even in T"
        );
        assert!(
            (torsional_strain_energy(200.0, 2.0, 80.0e9, 1.0e-6)
                - 4.0 * torsional_strain_energy(100.0, 2.0, 80.0e9, 1.0e-6))
            .abs()
                < 1e-12,
            "quadratic in T"
        );
        assert!(torsional_strain_energy(100.0, 2.0, 80.0e9, 1.0e-6) > 0.0, "non-negative");

        // 0 sentinel for non-physical input.
        assert_eq!(torsional_strain_energy(100.0, 0.0, 80.0e9, 1.0e-6), 0.0);
        assert_eq!(torsional_strain_energy(100.0, 2.0, 0.0, 1.0e-6), 0.0);
        assert_eq!(torsional_strain_energy(100.0, 2.0, 80.0e9, 0.0), 0.0);
        assert_eq!(torsional_strain_energy(f64::NAN, 2.0, 80.0e9, 1.0e-6), 0.0);
    }

    #[test]
    fn beam_angle_of_twist_matches_the_closed_form() {
        // Worked point: T = 100 N·m over a 2 m shaft, G = 80 GPa steel,
        // J = 1e-8 m⁴ → θ = T·L/(G·J) = 200 / 800 = 0.25 rad.
        let (t, l, g, j) = (100.0, 2.0, 80.0e9, 1.0e-8);
        let theta = beam_angle_of_twist(t, l, g, j);
        assert!((theta - 0.25).abs() / theta < 1e-9, "θ = 0.25 rad, got {theta}");
        // Linear in the torque, and sign-preserving (reverse T → reverse twist).
        assert!((beam_angle_of_twist(2.0 * t, l, g, j) - 2.0 * theta).abs() / theta < 1e-12);
        assert!(
            (beam_angle_of_twist(-t, l, g, j) + theta).abs() / theta < 1e-12,
            "sign-preserving"
        );
        // Linear in the length: double L → double θ.
        assert!(
            (beam_angle_of_twist(t, 2.0 * l, g, j) - 2.0 * theta).abs() / theta < 1e-12,
            "L scaling"
        );
        // Inverse in the torsional rigidity G·J: double G or J → half θ.
        assert!((beam_angle_of_twist(t, l, 2.0 * g, j) - 0.5 * theta).abs() / theta < 1e-12, "1/G");
        assert!((beam_angle_of_twist(t, l, g, 2.0 * j) - 0.5 * theta).abs() / theta < 1e-12, "1/J");
        // Non-physical input → 0.
        assert_eq!(beam_angle_of_twist(t, l, g, -1.0e-8), 0.0); // J ≤ 0
        assert_eq!(beam_angle_of_twist(t, l, 0.0, j), 0.0); // G ≤ 0
        assert_eq!(beam_angle_of_twist(t, -1.0, g, j), 0.0); // L ≤ 0
        assert_eq!(beam_angle_of_twist(f64::NAN, l, g, j), 0.0); // non-finite T
        assert_eq!(beam_angle_of_twist(t, l, f64::INFINITY, j), 0.0); // non-finite G
    }

    #[test]
    fn simply_supported_beam_centre_deflection_matches_analytic() {
        // A simply-supported beam, central point load P: the analytic
        // mid-span deflection is δ = P·L³/(48·E·I).
        let mat = steel();
        let l = 6.0;
        let n_elem = 12; // even so a node sits at mid-span
        let nodes: Vec<Vector3<f64>> = (0..=n_elem)
            .map(|i| Vector3::new(l * i as f64 / n_elem as f64, 0.0, 0.0))
            .collect();
        let section = BeamSection::rectangle(0.08, 0.08);
        let elements: Vec<BeamElement> = (0..n_elem)
            .map(|i| BeamElement::new(i, i + 1, section))
            .collect();
        let p = 2000.0;
        // Pinned at both ends. To keep the frame from spinning about
        // its own axis / sliding axially, fix axial + torsion at the
        // left support too.
        let mut left = BeamConstraint::pinned(0);
        left.fixed[0] = Some(0.0); // ux
        left.fixed[3] = Some(0.0); // θx (torsion)
        let constraints = [left, BeamConstraint::pinned(n_elem)];
        let mid = n_elem / 2;
        let loads = [BeamLoad::force(mid, [0.0, 0.0, -p])];
        let sol = solve_beam_static(&nodes, &elements, &mat, &constraints, &loads).unwrap();
        let i = section.iy;
        let bending = simply_supported_center_deflection(p, l, mat.youngs_modulus, i);
        // Allow a few % for the Timoshenko shear contribution.
        let centre = sol.translation[mid][2].abs();
        let rel = (centre - bending).abs() / bending;
        assert!(
            rel < 0.05,
            "mid-span deflection {centre} vs Euler-Bernoulli {bending} (rel {rel})"
        );
    }

    #[test]
    fn fixed_fixed_center_deflection_is_the_stiff_clamped_clamped_case() {
        // Worked point: P = 1 kN central load on a 4 m clamped–clamped steel beam,
        // E = 200 GPa, I = 1e-6 m⁴ → δ = P·L³/(192·E·I) = 64000/3.84e7 = 1/600 m.
        let d = fixed_fixed_center_deflection(1000.0, 4.0, 200.0e9, 1.0e-6);
        assert!((d - 1.0 / 600.0).abs() / (1.0 / 600.0) < 1e-12, "worked δ = 1/600 m, got {d}");

        // STRONG cross-checks threading the other two boundary conditions: clamping
        // both ends is exactly 4× stiffer at mid-span than pinning both (δ = δ_ss/4),
        // and 64× stiffer than a cantilever tip (δ = δ_cant/64), under the same load.
        for &(p, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 0.8, 120.0e9, 9.0e-8),
        ] {
            let ff = fixed_fixed_center_deflection(p, l, e, i);
            let ss = simply_supported_center_deflection(p, l, e, i);
            let cant = cantilever_tip_deflection(p, l, e, i);
            assert!((ff - ss / 4.0).abs() <= 1e-12 * (ss / 4.0).abs(), "δ_ff = δ_ss/4");
            assert!((ff - cant / 64.0).abs() <= 1e-12 * (cant / 64.0).abs(), "δ_ff = δ_cant/64");
            // Sign-preserving: a downward (negative) load gives a downward deflection.
            assert!(ff * p > 0.0, "deflection follows the load sign");
        }

        // Cubic in span: doubling L multiplies the deflection by 8.
        let d1 = fixed_fixed_center_deflection(1000.0, 1.0, 200.0e9, 1.0e-6);
        let d2 = fixed_fixed_center_deflection(1000.0, 2.0, 200.0e9, 1.0e-6);
        assert!((d2 - 8.0 * d1).abs() / (8.0 * d1) < 1e-12, "δ ∝ L³");

        // Non-physical input → 0.
        assert_eq!(fixed_fixed_center_deflection(f64::NAN, 4.0, 200.0e9, 1.0e-6), 0.0); // P NaN
        assert_eq!(fixed_fixed_center_deflection(1000.0, 0.0, 200.0e9, 1.0e-6), 0.0); // L = 0
        assert_eq!(fixed_fixed_center_deflection(1000.0, 4.0, -1.0, 1.0e-6), 0.0); // E < 0
        assert_eq!(fixed_fixed_center_deflection(1000.0, 4.0, 200.0e9, 0.0), 0.0); // I = 0
    }

    #[test]
    fn simply_supported_udl_max_moment_matches_statics() {
        // Worked point: w = 1 kN/m on a 4 m span → M_max = w·L²/8 = 2000 N·m.
        let m = simply_supported_udl_max_moment(1000.0, 4.0);
        assert!((m - 2000.0).abs() < 1e-9, "M_max = w·L²/8, got {m}");

        // Threads fixed_fixed_udl_end_moment (#278): the simply-supported UDL mid-span
        // moment wL²/8 is 3/2 the clamped-clamped fixing moment wL²/12.
        assert!(
            (m - 1.5 * fixed_fixed_udl_end_moment(1000.0, 4.0)).abs() < 1e-9,
            "M_ss = 1.5 · M_ff_end"
        );

        // Threads simply_supported_udl_center_deflection: δ = 5wL⁴/384EI = 5·M·L²/(48EI).
        for &(w, l, e, i) in &[
            (1000.0_f64, 4.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 0.8, 120.0e9, 9.0e-8),
        ] {
            let from_moment = 5.0 * simply_supported_udl_max_moment(w, l) * l * l / (48.0 * e * i);
            let delta = simply_supported_udl_center_deflection(w, l, e, i);
            assert!((from_moment - delta).abs() <= 1e-12 * delta.abs(), "5·M·L²/48EI = δ_udl");
        }

        // Quadratic in span; linear and sign-preserving in w.
        assert!(
            (simply_supported_udl_max_moment(1000.0, 8.0)
                - 4.0 * simply_supported_udl_max_moment(1000.0, 4.0))
            .abs()
                < 1e-9,
            "quadratic in L"
        );
        assert!(simply_supported_udl_max_moment(-1000.0, 4.0) < 0.0, "sign follows the load");

        // Non-physical input → 0.
        assert_eq!(simply_supported_udl_max_moment(f64::NAN, 4.0), 0.0);
        assert_eq!(simply_supported_udl_max_moment(1000.0, 0.0), 0.0);
        assert_eq!(simply_supported_udl_max_moment(1000.0, -1.0), 0.0);
    }

    #[test]
    fn simply_supported_point_load_max_moment_matches_statics() {
        // Worked point: P = 1 kN central load on a 4 m span → M_max = P·L/4 = 1000 N·m.
        let m = simply_supported_point_load_max_moment(1000.0, 4.0);
        assert!((m - 1000.0).abs() < 1e-9, "M_max = P·L/4, got {m}");

        // Classic "clamping halves the peak moment": the simply-supported mid-span
        // moment P·L/4 is exactly twice the clamped-clamped fixing moment P·L/8 (the
        // fixed-fixed beam shares the same magnitude between mid-span and both ends).
        assert!(
            (m - 2.0 * fixed_fixed_point_load_end_moment(1000.0, 4.0)).abs() < 1e-9,
            "M_ss = 2 · M_ff_end"
        );

        // Deflection relation threading simply_supported_center_deflection: the
        // mid-span deflection δ = P·L³/48EI equals M·L²/(12EI).
        for &(p, l, e, i) in &[
            (1000.0_f64, 4.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-650.0, 2.5, 70.0e9, 4.2e-7),
            (8200.0, 1.2, 120.0e9, 9.0e-8),
        ] {
            let from_moment =
                simply_supported_point_load_max_moment(p, l) * l * l / (12.0 * e * i);
            let delta = simply_supported_center_deflection(p, l, e, i);
            assert!((from_moment - delta).abs() <= 1e-12 * delta.abs(), "M·L²/12EI = δ_center");
        }

        // Linear and sign-preserving in P; linear in L.
        assert!(simply_supported_point_load_max_moment(-1000.0, 4.0) < 0.0, "sign follows the load");
        assert!(
            (simply_supported_point_load_max_moment(1000.0, 8.0)
                - 2.0 * simply_supported_point_load_max_moment(1000.0, 4.0))
            .abs()
                < 1e-9,
            "linear in L"
        );

        // Non-physical input → 0.
        assert_eq!(simply_supported_point_load_max_moment(f64::NAN, 4.0), 0.0);
        assert_eq!(simply_supported_point_load_max_moment(1000.0, 0.0), 0.0);
        assert_eq!(simply_supported_point_load_max_moment(1000.0, -1.0), 0.0);
    }

    #[test]
    fn simply_supported_center_deflection_matches_the_closed_form() {
        // Worked point: P = 1 kN central load on a 4 m simply-supported steel
        // beam, E = 200 GPa, I = 1e-6 m⁴ → δ = P·L³/(48·E·I) = 64000/9.6e6 = 1/150 m.
        let (p, l, e, i) = (1000.0, 4.0, 200.0e9, 1.0e-6);
        let delta = simply_supported_center_deflection(p, l, e, i);
        assert!((delta - 1.0 / 150.0).abs() / delta < 1e-9, "δ = 1/150 m, got {delta}");
        // 16× stiffer than a cantilever of the same span/section under the same load.
        assert!(
            (delta - cantilever_tip_deflection(p, l, e, i) / 16.0).abs() / delta < 1e-12,
            "δ_ss = δ_cantilever / 16"
        );
        // Linear in the load (and sign-preserving).
        assert!((simply_supported_center_deflection(2.0 * p, l, e, i) - 2.0 * delta).abs() / delta < 1e-12);
        assert!(
            (simply_supported_center_deflection(-p, l, e, i) + delta).abs() / delta < 1e-12,
            "sign-preserving"
        );
        // Cubic in the span; inverse in the flexural rigidity E·I.
        assert!(
            (simply_supported_center_deflection(p, 2.0 * l, e, i) - 8.0 * delta).abs() / delta < 1e-9,
            "L³ scaling"
        );
        assert!((simply_supported_center_deflection(p, l, 2.0 * e, i) - 0.5 * delta).abs() / delta < 1e-12, "1/E");
        assert!((simply_supported_center_deflection(p, l, e, 2.0 * i) - 0.5 * delta).abs() / delta < 1e-12, "1/I");
        // Non-physical input → 0.
        assert_eq!(simply_supported_center_deflection(p, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(simply_supported_center_deflection(p, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(simply_supported_center_deflection(p, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(simply_supported_center_deflection(f64::NAN, l, e, i), 0.0); // non-finite P
        assert_eq!(simply_supported_center_deflection(p, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn simply_supported_end_slope_matches_closed_form() {
        // Worked point: P = 1 kN central load on a 4 m simply-supported steel beam,
        // E = 200 GPa, I = 1e-6 m⁴ → θ = P·L²/(16·E·I) = 16000/3.2e6 = 0.005 rad.
        let (p, l, e, i) = (1000.0, 4.0, 200.0e9, 1.0e-6);
        let theta = simply_supported_end_slope(p, l, e, i);
        assert!((theta - 0.005).abs() / theta < 1e-12, "θ = 0.005 rad, got {theta}");
        // Linear in the load (and sign-preserving).
        assert!((simply_supported_end_slope(2.0 * p, l, e, i) - 2.0 * theta).abs() / theta < 1e-12);
        assert!((simply_supported_end_slope(-p, l, e, i) + theta).abs() / theta < 1e-12, "sign-preserving");
        // Quadratic in the span; inverse in the flexural rigidity E·I.
        assert!((simply_supported_end_slope(p, 2.0 * l, e, i) - 4.0 * theta).abs() / theta < 1e-12, "L² scaling");
        assert!((simply_supported_end_slope(p, l, 2.0 * e, i) - 0.5 * theta).abs() / theta < 1e-12, "1/E");
        assert!((simply_supported_end_slope(p, l, e, 2.0 * i) - 0.5 * theta).abs() / theta < 1e-12, "1/I");
        // STRONG non-tautological cross-check tying it to #194: for a centrally loaded
        // simply-supported beam the mid-span deflection is a third of the span times
        // the support rotation, δ = (L/3)·θ. The slope impl uses L²/(16EI); the check
        // uses the independent deflection fn L³/(48EI) — a known beam relation.
        for &(pp, ll) in &[(1000.0_f64, 4.0_f64), (250.0, 1.5), (-800.0, 3.0), (4200.0, 0.8)] {
            let slope = simply_supported_end_slope(pp, ll, e, i);
            let defl = simply_supported_center_deflection(pp, ll, e, i);
            assert!(
                (defl - ll / 3.0 * slope).abs() / defl.abs() < 1e-12,
                "δ = (L/3)·θ at P={pp}, L={ll}: {defl} vs {slope}"
            );
        }
        // The support rotation is 1/8 of the cantilever tip slope (#212) for the same
        // P, L, E, I (P·L²/16EI vs P·L²/2EI) — ties the two slope references.
        assert!(
            (theta - cantilever_tip_slope(p, l, e, i) / 8.0).abs() / theta < 1e-12,
            "θ_ss = θ_cantilever / 8"
        );
        // Non-physical input → 0.
        assert_eq!(simply_supported_end_slope(p, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(simply_supported_end_slope(p, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(simply_supported_end_slope(p, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(simply_supported_end_slope(f64::NAN, l, e, i), 0.0); // non-finite P
        assert_eq!(simply_supported_end_slope(p, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn fixed_fixed_udl_end_moment_matches_superposition() {
        // Worked point: w = 1 kN/m UDL on a 4 m clamped–clamped beam → the fixing
        // moment at each end is M = w·L²/12 = 16000/12 = 4000/3 N·m.
        let m = fixed_fixed_udl_end_moment(1000.0, 4.0);
        assert!((m - 4000.0 / 3.0).abs() < 1e-9, "M_end = w·L²/12, got {m}");

        // STRONG cross-check threading the two UDL centre-deflection formulas via
        // superposition: two equal end moments lift a simply-supported span's
        // mid-span by M·L²/(8EI), which must equal δ_ss_udl − δ_ff_udl
        // (5wL⁴/384EI − wL⁴/384EI = wL⁴/96EI = (wL²/12)·L²/8EI).
        for &(w, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 0.8, 120.0e9, 9.0e-8),
        ] {
            let lift = fixed_fixed_udl_end_moment(w, l) * l * l / (8.0 * e * i);
            let diff = simply_supported_udl_center_deflection(w, l, e, i)
                - fixed_fixed_udl_center_deflection(w, l, e, i);
            assert!((lift - diff).abs() <= 1e-12 * diff.abs(), "M·L²/8EI = δ_ss_udl − δ_ff_udl");
        }

        // Quadratic in span; linear and sign-preserving in w.
        assert!(
            (fixed_fixed_udl_end_moment(1000.0, 2.0)
                - 4.0 * fixed_fixed_udl_end_moment(1000.0, 1.0))
            .abs()
                < 1e-9,
            "quadratic in L"
        );
        assert!(fixed_fixed_udl_end_moment(-1000.0, 4.0) < 0.0, "sign follows the load");

        // Non-physical input → 0.
        assert_eq!(fixed_fixed_udl_end_moment(f64::NAN, 4.0), 0.0); // w NaN
        assert_eq!(fixed_fixed_udl_end_moment(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(fixed_fixed_udl_end_moment(1000.0, -1.0), 0.0); // L < 0
    }

    #[test]
    fn propped_cantilever_udl_prop_reaction_matches_compatibility() {
        // (a) WORKED: w = 1 kN/m UDL on a 2 m propped cantilever → R_B = 3·w·L/8 =
        // 3·1000·2/8 = 750 N.
        assert!(
            (propped_cantilever_udl_prop_reaction(1000.0, 2.0) - 750.0).abs() <= 1e-9 * 750.0,
            "R_B = 3wL/8 = 750 N"
        );

        // (b) COMPATIBILITY THREAD (non-tautological, force method): the prop reaction is
        // exactly the tip point-force that cancels the cantilever's UDL tip deflection, so
        // cantilever_tip_deflection(R, L, E, I) == cantilever_udl_tip_deflection(w, L, E, I).
        for &(w, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 0.8, 120.0e9, 9.0e-8),
        ] {
            let r = propped_cantilever_udl_prop_reaction(w, l);
            let udl = cantilever_udl_tip_deflection(w, l, e, i);
            assert!(
                (cantilever_tip_deflection(r, l, e, i) - udl).abs() <= 1e-9 * udl.abs(),
                "prop reaction zeroes the cantilever tip deflection"
            );
        }

        // (c) LINEARITY: linear and sign-preserving in w, linear in L.
        let base = propped_cantilever_udl_prop_reaction(1000.0, 2.0);
        assert!(
            (propped_cantilever_udl_prop_reaction(2000.0, 2.0) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in w"
        );
        assert!(
            (propped_cantilever_udl_prop_reaction(1000.0, 4.0) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "linear in L"
        );
        assert!(propped_cantilever_udl_prop_reaction(-1000.0, 2.0) < 0.0, "sign follows the load");

        // (d) Non-physical input → 0.
        assert_eq!(propped_cantilever_udl_prop_reaction(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(propped_cantilever_udl_prop_reaction(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_udl_prop_reaction(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn propped_cantilever_udl_fixed_end_moment_matches_equilibrium() {
        // (a) WORKED: w = 1 kN/m UDL on a 2 m propped cantilever → M_A = w·L²/8 =
        // 1000·4/8 = 500 N·m.
        assert!(
            (propped_cantilever_udl_fixed_end_moment(1000.0, 2.0) - 500.0).abs() <= 1e-9 * 500.0,
            "M_A = wL²/8 = 500 N·m"
        );

        // (b) MOMENT-EQUILIBRIUM THREAD (non-tautological): moments about the fixed end
        // give M_A = w·L²/2 − R_B·L, threading the prop reaction R_B = 3wL/8.
        for &(w, l) in &[(1000.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let m = propped_cantilever_udl_fixed_end_moment(w, l);
            let equil = w * l * l / 2.0 - propped_cantilever_udl_prop_reaction(w, l) * l;
            assert!((m - equil).abs() <= 1e-9 * m.abs().max(1.0), "M_A = wL²/2 − R_B·L");
        }

        // (c) SCALING: quadratic in L, linear and sign-preserving in w.
        let base = propped_cantilever_udl_fixed_end_moment(1000.0, 2.0);
        assert!(
            (propped_cantilever_udl_fixed_end_moment(1000.0, 4.0) - 4.0 * base).abs()
                <= 1e-9 * 4.0 * base,
            "quadratic in L"
        );
        assert!(
            (propped_cantilever_udl_fixed_end_moment(2000.0, 2.0) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in w"
        );
        assert!(propped_cantilever_udl_fixed_end_moment(-1000.0, 2.0) < 0.0, "sign follows load");

        // (d) Non-physical input → 0.
        assert_eq!(propped_cantilever_udl_fixed_end_moment(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(propped_cantilever_udl_fixed_end_moment(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_udl_fixed_end_moment(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn propped_cantilever_udl_max_sagging_moment_matches_statics() {
        // (a) WORKED: w = 1 kN/m UDL on a 2 m propped cantilever → M_sag = 9·w·L²/128 =
        // 9·1000·4/128 = 281.25 N·m.
        assert!(
            (propped_cantilever_udl_max_sagging_moment(1000.0, 2.0) - 281.25).abs() <= 1e-9 * 281.25,
            "M_sag = 9wL²/128 = 281.25 N·m"
        );

        // (b) STATICS THREAD (non-tautological): the max sagging moment is at the zero-
        // shear section, 3L/8 from the prop; moments there from the prop end give
        // M_sag = R_B·(3L/8) − w·(3L/8)²/2, threading the prop reaction R_B = 3wL/8.
        for &(w, l) in &[(1000.0_f64, 2.0_f64), (8200.0, 0.8), (450.0, 3.5)] {
            let a = 3.0 * l / 8.0;
            let m = propped_cantilever_udl_max_sagging_moment(w, l);
            let stat = propped_cantilever_udl_prop_reaction(w, l) * a - w * a * a / 2.0;
            assert!((m - stat).abs() <= 1e-9 * m.abs().max(1.0), "M_sag = R_B·(3L/8) − w·(3L/8)²/2");
        }

        // (c) SCALING: quadratic in L, linear and sign-preserving in w.
        let base = propped_cantilever_udl_max_sagging_moment(1000.0, 2.0);
        assert!(
            (propped_cantilever_udl_max_sagging_moment(1000.0, 4.0) - 4.0 * base).abs()
                <= 1e-9 * 4.0 * base,
            "quadratic in L"
        );
        assert!(
            (propped_cantilever_udl_max_sagging_moment(2000.0, 2.0) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in w"
        );
        assert!(propped_cantilever_udl_max_sagging_moment(-1000.0, 2.0) < 0.0, "sign follows load");

        // (d) Non-physical input → 0.
        assert_eq!(propped_cantilever_udl_max_sagging_moment(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(propped_cantilever_udl_max_sagging_moment(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_udl_max_sagging_moment(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn two_span_continuous_beam_udl_middle_moment_matches_the_three_moment_theorem() {
        // WORKED: w = 1 kN/m UDL over both spans of a two-span continuous beam with
        // equal spans L = 2 m → the middle-support moment is M_B = wL²/8 = 500 N·m.
        assert!(
            (two_span_continuous_beam_udl_middle_moment(1000.0, 2.0) - 500.0).abs() <= 1e-9 * 500.0,
            "M_B = wL²/8 = 500 N·m",
        );

        // STRONG non-tautological CROSS-CONFIGURATION thread: by symmetry each span of the
        // two-span beam behaves as a propped cantilever fixed at the centre support, so M_B
        // equals the existing propped-cantilever fixed-end moment w·L²/8 — same value,
        // different structure, reached independently.
        for &(w, l) in &[(1200.0_f64, 3.5_f64), (8200.0, 0.8), (-450.0, 2.0)] {
            let m_b = two_span_continuous_beam_udl_middle_moment(w, l);
            assert!(
                (m_b - propped_cantilever_udl_fixed_end_moment(w, l)).abs() <= 1e-9 * m_b.abs(),
                "M_B must equal the propped-cantilever fixed-end moment by symmetry",
            );
        }

        // Linear in w, quadratic in L.
        let base = two_span_continuous_beam_udl_middle_moment(1000.0, 2.0);
        assert!(
            (two_span_continuous_beam_udl_middle_moment(2000.0, 2.0) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in w",
        );
        assert!(
            (two_span_continuous_beam_udl_middle_moment(1000.0, 4.0) - 4.0 * base).abs()
                <= 1e-9 * 4.0 * base,
            "quadratic in L",
        );
        assert!(
            two_span_continuous_beam_udl_middle_moment(-1000.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(two_span_continuous_beam_udl_middle_moment(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(two_span_continuous_beam_udl_middle_moment(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(two_span_continuous_beam_udl_middle_moment(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn two_span_continuous_beam_udl_middle_reaction_collects_both_spans() {
        // WORKED: w = 1 kN/m UDL over both spans, equal spans L = 2 m → the middle-support
        // reaction is R_B = 5wL/4 = 2500 N.
        assert!(
            (two_span_continuous_beam_udl_middle_reaction(1000.0, 2.0) - 2500.0).abs()
                <= 1e-9 * 2500.0,
            "R_B = 5wL/4 = 2500 N",
        );

        // STRONG non-tautological symmetry thread: the centre support collects the fixed-end
        // reaction of BOTH propped-cantilever spans, R_B = 2 · (5wL/8); and global
        // equilibrium R_A + R_B + R_C = 2wL with the outer reactions = the prop reaction.
        for &(w, l) in &[(1200.0_f64, 3.5_f64), (8200.0, 0.8), (-450.0, 2.0)] {
            let r_b = two_span_continuous_beam_udl_middle_reaction(w, l);
            assert!(
                (r_b - 2.0 * propped_cantilever_udl_fixed_end_reaction(w, l)).abs()
                    <= 1e-9 * r_b.abs(),
                "R_B must equal twice the propped-cantilever fixed-end reaction",
            );
            let r_outer = propped_cantilever_udl_prop_reaction(w, l);
            assert!(
                (2.0 * r_outer + r_b - 2.0 * w * l).abs() <= 1e-9 * (2.0 * w * l).abs(),
                "R_A + R_B + R_C = total load 2wL",
            );
        }

        // Linear in both w and L.
        let base = two_span_continuous_beam_udl_middle_reaction(1000.0, 2.0);
        assert!(
            (two_span_continuous_beam_udl_middle_reaction(2000.0, 2.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in w",
        );
        assert!(
            (two_span_continuous_beam_udl_middle_reaction(1000.0, 4.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in L",
        );
        assert!(
            two_span_continuous_beam_udl_middle_reaction(-1000.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(two_span_continuous_beam_udl_middle_reaction(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(two_span_continuous_beam_udl_middle_reaction(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(two_span_continuous_beam_udl_middle_reaction(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn two_span_continuous_beam_udl_outer_reaction_completes_the_reaction_set() {
        // WORKED: w = 1 kN/m UDL over both spans, equal spans L = 2 m → each outer
        // (end-support) reaction is R_A = 3wL/8 = 750 N.
        assert!(
            (two_span_continuous_beam_udl_outer_reaction(1000.0, 2.0) - 750.0).abs()
                <= 1e-9 * 750.0,
            "R_A = 3wL/8 = 750 N",
        );

        // STRONG non-tautological DUAL thread over three signed (w, L) cases: (i) the outer
        // support is the propped (simple) end of each span, so R_A equals the existing
        // propped-cantilever prop reaction; (ii) global equilibrium 2·R_A + R_B = 2wL with
        // the #469 centre reaction.
        for &(w, l) in &[(1200.0_f64, 3.5_f64), (8200.0, 0.8), (-450.0, 2.0)] {
            let ra = two_span_continuous_beam_udl_outer_reaction(w, l);
            assert!(
                (ra - propped_cantilever_udl_prop_reaction(w, l)).abs() <= 1e-9 * ra.abs(),
                "R_A must equal the propped-cantilever prop reaction",
            );
            assert!(
                (2.0 * ra + two_span_continuous_beam_udl_middle_reaction(w, l) - 2.0 * w * l).abs()
                    <= 1e-9 * (2.0 * w * l).abs(),
                "2·R_A + R_B = total load 2wL",
            );
        }

        // Linear in both w and L.
        let base = two_span_continuous_beam_udl_outer_reaction(1000.0, 2.0);
        assert!(
            (two_span_continuous_beam_udl_outer_reaction(2000.0, 2.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in w",
        );
        assert!(
            (two_span_continuous_beam_udl_outer_reaction(1000.0, 4.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in L",
        );
        assert!(
            two_span_continuous_beam_udl_outer_reaction(-1000.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(two_span_continuous_beam_udl_outer_reaction(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(two_span_continuous_beam_udl_outer_reaction(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(two_span_continuous_beam_udl_outer_reaction(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn three_span_continuous_beam_udl_interior_moment_matches_clapeyron() {
        // WORKED: w = 1 kN/m UDL over all three spans, equal spans L = 2 m → each interior
        // support moment is M = wL²/10 = 400 N·m.
        assert!(
            (three_span_continuous_beam_udl_interior_moment(1000.0, 2.0) - 400.0).abs()
                <= 1e-12 * 400.0,
            "M = wL²/10 = 400 N·m",
        );

        // STRONG non-tautological CROSS-CONFIGURATION thread: the three-span interior moment
        // is EXACTLY 0.8× the two-span middle moment (wL²/10 vs wL²/8 = ratio 8/10), the extra
        // span sharing the load — threads the existing two-span fn.
        for &(w, l) in &[(1200.0_f64, 3.5_f64), (8200.0, 0.8), (-450.0, 2.0)] {
            let m3 = three_span_continuous_beam_udl_interior_moment(w, l);
            let m2 = two_span_continuous_beam_udl_middle_moment(w, l);
            assert!(
                (m3 - 0.8 * m2).abs() <= 1e-9 * m3.abs(),
                "three-span M = 0.8 × two-span M",
            );
        }

        // Linear in w, quadratic in L.
        let base = three_span_continuous_beam_udl_interior_moment(1000.0, 2.0);
        assert!(
            (three_span_continuous_beam_udl_interior_moment(2000.0, 2.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in w",
        );
        assert!(
            (three_span_continuous_beam_udl_interior_moment(1000.0, 4.0) - 4.0 * base).abs()
                < 1e-9 * base,
            "quadratic in L",
        );
        assert!(
            three_span_continuous_beam_udl_interior_moment(-1000.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(three_span_continuous_beam_udl_interior_moment(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(three_span_continuous_beam_udl_interior_moment(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(three_span_continuous_beam_udl_interior_moment(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn three_span_continuous_beam_udl_reactions_balance() {
        // WORKED (w = 1 kN/m, L = 2 m): R_end = 2wL/5 = 800 N, R_int = 11wL/10 = 2200 N.
        assert!(
            (three_span_continuous_beam_udl_end_reaction(1000.0, 2.0) - 800.0).abs()
                <= 1e-12 * 800.0,
            "R_end = 2wL/5 = 800 N",
        );
        assert!(
            (three_span_continuous_beam_udl_interior_reaction(1000.0, 2.0) - 2200.0).abs()
                <= 1e-12 * 2200.0,
            "R_int = 11wL/10 = 2200 N",
        );

        // STRONG non-tautological vertical-equilibrium thread over signed (w, L): the four
        // reactions (2 end + 2 interior) carry the whole 3wL load.
        for &(w, l) in &[(1200.0_f64, 3.5_f64), (8200.0, 0.8), (-450.0, 2.0)] {
            let r_end = three_span_continuous_beam_udl_end_reaction(w, l);
            let r_int = three_span_continuous_beam_udl_interior_reaction(w, l);
            assert!(
                (2.0 * r_end + 2.0 * r_int - 3.0 * w * l).abs() <= 1e-9 * (3.0 * w * l).abs(),
                "2·R_end + 2·R_int = 3wL (vertical equilibrium)",
            );
        }

        // Linear in w and L.
        let base = three_span_continuous_beam_udl_end_reaction(1000.0, 2.0);
        assert!(
            (three_span_continuous_beam_udl_end_reaction(2000.0, 2.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in w",
        );
        assert!(
            (three_span_continuous_beam_udl_end_reaction(1000.0, 4.0) - 2.0 * base).abs()
                < 1e-9 * base,
            "linear in L",
        );
        assert!(
            three_span_continuous_beam_udl_interior_reaction(-1000.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(three_span_continuous_beam_udl_end_reaction(f64::NAN, 2.0), 0.0);
        assert_eq!(three_span_continuous_beam_udl_interior_reaction(1000.0, 0.0), 0.0);
        assert_eq!(three_span_continuous_beam_udl_end_reaction(1000.0, -2.0), 0.0);
    }

    #[test]
    fn two_span_continuous_beam_central_point_load_middle_moment_matches_clapeyron() {
        // WORKED: P = 32 kN at the mid-span of one span of a two-span continuous beam,
        // equal spans L = 2 m → the middle-support moment is M_B = 3PL/32 = 6 kN·m.
        assert!(
            (two_span_continuous_beam_central_point_load_middle_moment(32.0, 2.0) - 6.0).abs()
                <= 1e-9 * 6.0,
            "M_B = 3PL/32 = 6 kN·m",
        );

        // STRONG non-tautological DUAL thread over three signed (P, L): (i) it is HALF the
        // propped-cantilever clamping moment (#466), the unloaded equal span halving the
        // restraint; (ii) an independent three-moment-theorem recompute from the loaded
        // span's free-BM triangle (area PL²/8, centroid L/2).
        for &(p, l) in &[(7.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let m = two_span_continuous_beam_central_point_load_middle_moment(p, l);
            assert!(
                (m - 0.5 * propped_cantilever_central_load_fixed_end_moment(p, l)).abs()
                    <= 1e-9 * m.abs(),
                "M_B = ½ · propped-cantilever clamping moment",
            );
            let a1 = p * l * l / 8.0; // free-BM triangle area of the loaded span
            let xbar = l / 2.0; // its centroid from the end support
            let mb3 = -6.0 * a1 * xbar / (l * 2.0 * (l + l)); // three-moment theorem (signed)
            assert!(
                (m + mb3).abs() <= 1e-9 * m.abs(),
                "M_B magnitude must equal the three-moment-theorem result",
            );
        }

        // Linear in both P and L.
        let base = two_span_continuous_beam_central_point_load_middle_moment(10.0, 2.0);
        assert!(
            (two_span_continuous_beam_central_point_load_middle_moment(20.0, 2.0) - 2.0 * base)
                .abs()
                < 1e-9 * base,
            "linear in P",
        );
        assert!(
            (two_span_continuous_beam_central_point_load_middle_moment(10.0, 4.0) - 2.0 * base)
                .abs()
                < 1e-9 * base,
            "linear in L",
        );
        assert!(
            two_span_continuous_beam_central_point_load_middle_moment(-10.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(
            two_span_continuous_beam_central_point_load_middle_moment(f64::NAN, 2.0),
            0.0
        ); // P NaN
        assert_eq!(two_span_continuous_beam_central_point_load_middle_moment(32.0, 0.0), 0.0); // L = 0
        assert_eq!(two_span_continuous_beam_central_point_load_middle_moment(32.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn two_span_continuous_beam_central_point_load_reactions_balance() {
        // WORKED (P = 32 kN, L = 2 m): R_A = 13P/32 = 13, R_B = 11P/16 = 22,
        // R_C = −3P/32 = −3 kN — summing to 32 kN = P.
        assert!(
            (two_span_continuous_beam_central_point_load_loaded_span_outer_reaction(32.0, 2.0)
                - 13.0)
                .abs()
                < 1e-12,
            "R_A = 13P/32 = 13",
        );
        assert!(
            (two_span_continuous_beam_central_point_load_middle_reaction(32.0, 2.0) - 22.0).abs()
                < 1e-12,
            "R_B = 11P/16 = 22",
        );
        assert!(
            (two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction(32.0, 2.0)
                + 3.0)
                .abs()
                < 1e-12,
            "R_C = −3P/32 = −3",
        );

        // STRONG non-tautological DUAL thread over signed (P, L): (i) vertical equilibrium
        // R_A + R_B + R_C = P; (ii) loaded-span statics R_A·L − P·L/2 = −M_B recovers the
        // centre hogging moment (#473), threading the existing moment fn.
        for &(p, l) in &[(7.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let ra = two_span_continuous_beam_central_point_load_loaded_span_outer_reaction(p, l);
            let rb = two_span_continuous_beam_central_point_load_middle_reaction(p, l);
            let rc = two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction(p, l);
            assert!(
                (ra + rb + rc - p).abs() <= 1e-9 * p.abs(),
                "R_A + R_B + R_C = P (vertical equilibrium)",
            );
            let m_from_statics = ra * l - p * (l / 2.0);
            let m_b = -two_span_continuous_beam_central_point_load_middle_moment(p, l);
            assert!(
                (m_from_statics - m_b).abs() <= 1e-9 * m_b.abs(),
                "R_A·L − P·L/2 = −M_B (hogging) — threads the centre moment",
            );
        }

        // Linear in P (these reactions are L-independent); R_C reverses sign with the load.
        let base_rb = two_span_continuous_beam_central_point_load_middle_reaction(10.0, 2.0);
        assert!(
            (two_span_continuous_beam_central_point_load_middle_reaction(20.0, 2.0) - 2.0 * base_rb)
                .abs()
                < 1e-12,
            "R_B linear in P",
        );
        assert!(
            two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction(-10.0, 2.0)
                > 0.0,
            "R_C reverses sign with the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(
            two_span_continuous_beam_central_point_load_middle_reaction(f64::NAN, 2.0),
            0.0
        );
        assert_eq!(
            two_span_continuous_beam_central_point_load_loaded_span_outer_reaction(32.0, 0.0),
            0.0
        );
        assert_eq!(
            two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction(32.0, -2.0),
            0.0
        );
    }

    #[test]
    fn propped_cantilever_udl_fixed_end_reaction_completes_the_reaction_pair() {
        // (a) WORKED: w = 1 kN/m UDL on a 2 m propped cantilever → R_A = 5·w·L/8 =
        // 5·1000·2/8 = 1250 N.
        assert!(
            (propped_cantilever_udl_fixed_end_reaction(1000.0, 2.0) - 1250.0).abs() <= 1e-9 * 1250.0,
            "R_A = 5wL/8 = 1250 N"
        );

        // (b) VERTICAL-EQUILIBRIUM THREAD (non-tautological): the two reactions carry the
        // whole UDL, R_A + R_B = w·L, threading the prop reaction R_B = 3wL/8.
        for &(w, l) in &[(1000.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let ra = propped_cantilever_udl_fixed_end_reaction(w, l);
            let rb = propped_cantilever_udl_prop_reaction(w, l);
            assert!((ra + rb - w * l).abs() <= 1e-9 * (w * l).abs().max(1.0), "R_A + R_B = wL");
        }

        // (c) SCALING: linear and sign-preserving in w, linear in L.
        let base = propped_cantilever_udl_fixed_end_reaction(1000.0, 2.0);
        assert!(
            (propped_cantilever_udl_fixed_end_reaction(2000.0, 2.0) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in w"
        );
        assert!(
            (propped_cantilever_udl_fixed_end_reaction(1000.0, 4.0) - 2.0 * base).abs()
                <= 1e-9 * 2.0 * base,
            "linear in L"
        );
        assert!(propped_cantilever_udl_fixed_end_reaction(-1000.0, 2.0) < 0.0, "sign follows load");

        // (d) Non-physical input → 0.
        assert_eq!(propped_cantilever_udl_fixed_end_reaction(f64::NAN, 2.0), 0.0); // w NaN
        assert_eq!(propped_cantilever_udl_fixed_end_reaction(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_udl_fixed_end_reaction(1000.0, -2.0), 0.0); // L < 0
    }

    #[test]
    fn propped_cantilever_central_load_prop_reaction_matches_compatibility() {
        // Worked point: P = 16 kN central load on a propped cantilever → the prop
        // (simple-support) reaction is R_B = 5P/16 = 5 kN, independent of span.
        assert!(
            (propped_cantilever_central_load_prop_reaction(16.0, 3.0) - 5.0).abs() < 1e-12,
            "R_B = 5P/16 = 5 for P = 16",
        );

        // STRONG non-tautological thread: derive 5/16 independently from the general
        // first-degree-redundant reaction R_B = P·a²(3L − a)/(2L³) (force method: prop
        // released → cantilever, compatibility δ_B = 0) evaluated at the mid-span case
        // a = L/2. The closed form must agree with the 5P/16 the function returns.
        for &(p, l) in &[(7.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let a = l / 2.0;
            let rb = p * a * a * (3.0 * l - a) / (2.0 * l * l * l);
            let got = propped_cantilever_central_load_prop_reaction(p, l);
            assert!(
                (got - rb).abs() <= 1e-9 * rb.abs(),
                "R_B must equal P·a²(3L−a)/(2L³)|a=L/2; got {got}, want {rb}",
            );
        }

        // Linear and sign-preserving in P; independent of L for the central case.
        let base = propped_cantilever_central_load_prop_reaction(10.0, 2.0);
        assert!(
            (propped_cantilever_central_load_prop_reaction(20.0, 2.0) - 2.0 * base).abs() < 1e-12,
            "doubling P doubles R_B",
        );
        assert!(
            propped_cantilever_central_load_prop_reaction(-10.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(propped_cantilever_central_load_prop_reaction(8.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_central_load_prop_reaction(8.0, -2.0), 0.0); // L < 0
        assert_eq!(propped_cantilever_central_load_prop_reaction(f64::NAN, 2.0), 0.0); // P NaN
    }

    #[test]
    fn propped_cantilever_central_load_fixed_end_reaction_completes_the_reaction_pair() {
        // Worked point: P = 16 kN central load → the clamped-end reaction is
        // R_A = 11P/16 = 11 kN, independent of span.
        assert!(
            (propped_cantilever_central_load_fixed_end_reaction(16.0, 3.0) - 11.0).abs() < 1e-12,
            "R_A = 11P/16 = 11 for P = 16",
        );

        // STRONG non-tautological thread through the EXISTING prop reaction (#462):
        // vertical equilibrium of the whole beam requires R_A + R_B = P.
        for &(p, l) in &[(7.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let ra = propped_cantilever_central_load_fixed_end_reaction(p, l);
            let rb = propped_cantilever_central_load_prop_reaction(p, l);
            assert!(
                (ra + rb - p).abs() <= 1e-9 * p.abs(),
                "R_A + R_B must equal P; got {ra} + {rb} vs {p}",
            );
        }

        // Linear and sign-preserving in P; independent of L for the central case.
        let base = propped_cantilever_central_load_fixed_end_reaction(10.0, 2.0);
        assert!(
            (propped_cantilever_central_load_fixed_end_reaction(20.0, 2.0) - 2.0 * base).abs()
                < 1e-12,
            "doubling P doubles R_A",
        );
        assert!(
            propped_cantilever_central_load_fixed_end_reaction(-10.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(propped_cantilever_central_load_fixed_end_reaction(8.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_central_load_fixed_end_reaction(8.0, -2.0), 0.0); // L < 0
        assert_eq!(propped_cantilever_central_load_fixed_end_reaction(f64::NAN, 2.0), 0.0); // P NaN
    }

    #[test]
    fn propped_cantilever_central_load_fixed_end_moment_completes_the_case() {
        // Worked point: P = 16 kN central load on a 2 m beam → the clamping moment is
        // M_A = 3PL/16 = 3·16·2/16 = 6 kN·m.
        assert!(
            (propped_cantilever_central_load_fixed_end_moment(16.0, 2.0) - 6.0).abs() < 1e-12,
            "M_A = 3PL/16 = 6 for P = 16, L = 2",
        );

        // STRONG non-tautological statics thread through the EXISTING prop reaction
        // (#462): moment equilibrium of the whole beam about the fixed end requires
        // M_A = P·(L/2) − R_B·L.
        for &(p, l) in &[(7.0_f64, 2.0_f64), (-450.0, 3.5), (8200.0, 0.8)] {
            let ma = propped_cantilever_central_load_fixed_end_moment(p, l);
            let statics = p * (l / 2.0) - propped_cantilever_central_load_prop_reaction(p, l) * l;
            assert!(
                (ma - statics).abs() <= 1e-9 * ma.abs(),
                "M_A must equal P·(L/2) − R_B·L; got {ma} vs {statics}",
            );
        }

        // Linear in P and (unlike the reactions) linear in L.
        let base = propped_cantilever_central_load_fixed_end_moment(10.0, 2.0);
        assert!(
            (propped_cantilever_central_load_fixed_end_moment(20.0, 2.0) - 2.0 * base).abs() < 1e-12,
            "doubling P doubles M_A",
        );
        assert!(
            (propped_cantilever_central_load_fixed_end_moment(10.0, 4.0) - 2.0 * base).abs() < 1e-12,
            "doubling L doubles M_A",
        );
        assert!(
            propped_cantilever_central_load_fixed_end_moment(-10.0, 2.0) < 0.0,
            "sign follows the load",
        );

        // Guards: non-physical input → 0.
        assert_eq!(propped_cantilever_central_load_fixed_end_moment(8.0, 0.0), 0.0); // L = 0
        assert_eq!(propped_cantilever_central_load_fixed_end_moment(8.0, -2.0), 0.0); // L < 0
        assert_eq!(propped_cantilever_central_load_fixed_end_moment(f64::NAN, 2.0), 0.0); // P NaN
    }

    #[test]
    fn fixed_fixed_point_load_end_moment_matches_superposition() {
        // Worked point: P = 1 kN central load on a 4 m clamped–clamped beam → the
        // fixing moment at each end is M = P·L/8 = 500 N·m.
        let m = fixed_fixed_point_load_end_moment(1000.0, 4.0);
        assert!((m - 500.0).abs() < 1e-9, "M_end = P·L/8 = 500, got {m}");

        // STRONG cross-check threading the two centre-deflection formulas via
        // superposition: the end fixing moment is exactly what stiffens a pinned span
        // into a clamped one. Two equal end moments lift a simply-supported beam's
        // mid-span by M·L²/(8EI), which must equal δ_ss − δ_ff
        // (PL³/48EI − PL³/192EI = PL³/64EI = (PL/8)·L²/8EI).
        for &(p, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 0.8, 120.0e9, 9.0e-8),
        ] {
            let lift = fixed_fixed_point_load_end_moment(p, l) * l * l / (8.0 * e * i);
            let diff = simply_supported_center_deflection(p, l, e, i)
                - fixed_fixed_center_deflection(p, l, e, i);
            assert!((lift - diff).abs() <= 1e-12 * diff.abs(), "M·L²/8EI = δ_ss − δ_ff");
        }

        // Linear in P (sign-preserving) and linear in L.
        assert!(
            (fixed_fixed_point_load_end_moment(2000.0, 4.0)
                - 2.0 * fixed_fixed_point_load_end_moment(1000.0, 4.0))
            .abs()
                < 1e-9,
            "linear in P"
        );
        assert!(fixed_fixed_point_load_end_moment(-1000.0, 4.0) < 0.0, "sign follows the load");

        // Non-physical input → 0.
        assert_eq!(fixed_fixed_point_load_end_moment(f64::NAN, 4.0), 0.0); // P NaN
        assert_eq!(fixed_fixed_point_load_end_moment(1000.0, 0.0), 0.0); // L = 0
        assert_eq!(fixed_fixed_point_load_end_moment(1000.0, -1.0), 0.0); // L < 0
    }

    #[test]
    fn fixed_fixed_udl_center_deflection_is_the_stiff_clamped_clamped_udl() {
        // Worked point: w = 1 kN/m UDL on a 4 m clamped–clamped steel beam,
        // E = 200 GPa, I = 1e-6 m⁴ → δ = w·L⁴/(384·E·I) = 256000/7.68e7 = 1/300 m.
        let d = fixed_fixed_udl_center_deflection(1000.0, 4.0, 200.0e9, 1.0e-6);
        assert!((d - 1.0 / 300.0).abs() / (1.0 / 300.0) < 1e-12, "worked δ = 1/300 m, got {d}");

        // STRONG cross-check threading simply_supported_udl_center_deflection: clamping
        // both ends is exactly 5× stiffer at mid-span than pinning both under a UDL
        // (w·L⁴/384EI vs 5·w·L⁴/384EI), so δ_ff = δ_ss/5.
        for &(w, l, e, i) in &[
            (1000.0_f64, 2.0_f64, 200.0e9_f64, 1.0e-6_f64),
            (-450.0, 3.5, 70.0e9, 4.2e-7),
            (8200.0, 0.8, 120.0e9, 9.0e-8),
        ] {
            let ff = fixed_fixed_udl_center_deflection(w, l, e, i);
            let ss = simply_supported_udl_center_deflection(w, l, e, i);
            assert!((ff - ss / 5.0).abs() <= 1e-12 * (ss / 5.0).abs(), "δ_ff_udl = δ_ss_udl/5");
            // Sign-preserving: a downward (negative) load gives a downward deflection.
            assert!(ff * w > 0.0, "deflection follows the load sign");
        }

        // Quartic in span: doubling L multiplies the deflection by 16.
        let d1 = fixed_fixed_udl_center_deflection(1000.0, 1.0, 200.0e9, 1.0e-6);
        let d2 = fixed_fixed_udl_center_deflection(1000.0, 2.0, 200.0e9, 1.0e-6);
        assert!((d2 - 16.0 * d1).abs() / (16.0 * d1) < 1e-12, "δ ∝ L⁴");

        // Non-physical input → 0.
        assert_eq!(fixed_fixed_udl_center_deflection(f64::NAN, 4.0, 200.0e9, 1.0e-6), 0.0); // w NaN
        assert_eq!(fixed_fixed_udl_center_deflection(1000.0, 0.0, 200.0e9, 1.0e-6), 0.0); // L = 0
        assert_eq!(fixed_fixed_udl_center_deflection(1000.0, 4.0, -1.0, 1.0e-6), 0.0); // E < 0
        assert_eq!(fixed_fixed_udl_center_deflection(1000.0, 4.0, 200.0e9, 0.0), 0.0); // I = 0
    }

    #[test]
    fn simply_supported_udl_center_deflection_matches_closed_form() {
        // Worked point: w = 1 kN/m UDL on a 2 m simply-supported steel beam,
        // E = 200 GPa, I = 1e-6 m⁴ → δ = 5·w·L⁴/(384·E·I) = 80000/7.68e7 = 1/960 m.
        let (w, l, e, i) = (1000.0, 2.0, 200.0e9, 1.0e-6);
        let delta = simply_supported_udl_center_deflection(w, l, e, i);
        assert!((delta - 1.0 / 960.0).abs() / delta < 1e-9, "δ = 1/960 m, got {delta}");
        // Linear in the load intensity, and sign-preserving (an upward load lifts mid-span).
        assert!((simply_supported_udl_center_deflection(2.0 * w, l, e, i) - 2.0 * delta).abs() / delta < 1e-12);
        assert!(
            (simply_supported_udl_center_deflection(-w, l, e, i) + delta).abs() / delta < 1e-12,
            "sign-preserving"
        );
        // Quartic in the span: double L → 16× δ.
        assert!(
            (simply_supported_udl_center_deflection(w, 2.0 * l, e, i) - 16.0 * delta).abs() / delta < 1e-9,
            "L⁴ scaling"
        );
        // Inverse in the flexural rigidity E·I: double E or I → half δ.
        assert!((simply_supported_udl_center_deflection(w, l, 2.0 * e, i) - 0.5 * delta).abs() / delta < 1e-12, "1/E");
        assert!((simply_supported_udl_center_deflection(w, l, e, 2.0 * i) - 0.5 * delta).abs() / delta < 1e-12, "1/I");
        // STRONG non-tautological cross-check: the same TOTAL load W = w·L, spread as a
        // UDL, deflects mid-span to exactly 5/8 of the same total as a central point
        // load. The UDL impl uses 5wL⁴/(384EI); the check uses the independent point-
        // load fn simply_supported_center_deflection (P·L³/(48EI)) with P = w·L — a
        // known structural-mechanics ratio, a different code path.
        for &(ww, ll) in &[(1000.0_f64, 2.0_f64), (300.0, 1.5), (-750.0, 3.0), (4200.0, 0.8)] {
            let udl = simply_supported_udl_center_deflection(ww, ll, e, i);
            let point = simply_supported_center_deflection(ww * ll, ll, e, i);
            assert!(
                (udl - 5.0 / 8.0 * point).abs() / point.abs() < 1e-12,
                "UDL = 5/8 × point-load(W=w·L) at w={ww}, L={ll}: {udl} vs {point}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(simply_supported_udl_center_deflection(w, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(simply_supported_udl_center_deflection(w, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(simply_supported_udl_center_deflection(w, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(simply_supported_udl_center_deflection(f64::NAN, l, e, i), 0.0); // non-finite w
        assert_eq!(simply_supported_udl_center_deflection(w, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn simply_supported_udl_end_slope_matches_closed_form() {
        // Worked point: w = 1 kN/m UDL on a 4 m simply-supported steel beam, E = 200
        // GPa, I = 1e-6 m⁴ → θ = w·L³/(24·E·I) = 64000/4.8e6 = 1/75 rad.
        let (w, l, e, i) = (1000.0, 4.0, 200.0e9, 1.0e-6);
        let theta = simply_supported_udl_end_slope(w, l, e, i);
        assert!((theta - 1.0 / 75.0).abs() / theta < 1e-12, "θ = 1/75 rad, got {theta}");
        // Linear in w (sign-preserving); cubic in span; inverse in E·I.
        assert!((simply_supported_udl_end_slope(2.0 * w, l, e, i) - 2.0 * theta).abs() / theta < 1e-12);
        assert!((simply_supported_udl_end_slope(-w, l, e, i) + theta).abs() / theta < 1e-12, "sign-preserving");
        assert!((simply_supported_udl_end_slope(w, 2.0 * l, e, i) - 8.0 * theta).abs() / theta < 1e-12, "L³ scaling");
        assert!((simply_supported_udl_end_slope(w, l, 2.0 * e, i) - 0.5 * theta).abs() / theta < 1e-12, "1/E");
        assert!((simply_supported_udl_end_slope(w, l, e, 2.0 * i) - 0.5 * theta).abs() / theta < 1e-12, "1/I");
        // Cross-check (a) tying it to the UDL centre deflection #206: δ = (5/16)·L·θ.
        for &(ww, ll) in &[(1000.0_f64, 4.0_f64), (300.0, 1.5), (-750.0, 3.0)] {
            let slope = simply_supported_udl_end_slope(ww, ll, e, i);
            let defl = simply_supported_udl_center_deflection(ww, ll, e, i);
            assert!(
                (defl - 5.0 / 16.0 * ll * slope).abs() / defl.abs() < 1e-12,
                "δ = (5/16)·L·θ at w={ww}, L={ll}"
            );
        }
        // Cross-check (b) tying it to the point-load end slope #218: the same TOTAL
        // load W = w·L spread as a UDL rotates the supports to 2/3 of the point case.
        for &(ww, ll) in &[(1000.0_f64, 4.0_f64), (300.0, 1.5), (-750.0, 3.0)] {
            let udl = simply_supported_udl_end_slope(ww, ll, e, i);
            let point = simply_supported_end_slope(ww * ll, ll, e, i);
            assert!(
                (udl - 2.0 / 3.0 * point).abs() / udl.abs() < 1e-12,
                "θ_udl = (2/3)·θ_point(W=w·L) at w={ww}, L={ll}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(simply_supported_udl_end_slope(w, l, e, -1.0e-6), 0.0); // I ≤ 0
        assert_eq!(simply_supported_udl_end_slope(w, l, 0.0, i), 0.0); // E ≤ 0
        assert_eq!(simply_supported_udl_end_slope(w, -1.0, e, i), 0.0); // L ≤ 0
        assert_eq!(simply_supported_udl_end_slope(f64::NAN, l, e, i), 0.0); // non-finite w
        assert_eq!(simply_supported_udl_end_slope(w, l, f64::INFINITY, i), 0.0); // non-finite E
    }

    #[test]
    fn solve_rejects_unconstrained_and_empty() {
        let mat = steel();
        let nodes = vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)];
        let elements = [BeamElement::new(0, 1, BeamSection::circle(0.01))];
        assert!(matches!(
            solve_beam_static(&nodes, &elements, &mat, &[], &[]),
            Err(BeamSolverError::Unconstrained)
        ));
        assert!(matches!(
            solve_beam_static(&[], &elements, &mat, &[BeamConstraint::clamped(0)], &[]),
            Err(BeamSolverError::EmptyModel)
        ));
    }

    #[test]
    fn solve_rejects_zero_length_element() {
        let mat = steel();
        let nodes = vec![Vector3::zeros(), Vector3::zeros()];
        let elements = [BeamElement::new(0, 1, BeamSection::circle(0.01))];
        assert!(matches!(
            solve_beam_static(&nodes, &elements, &mat, &[BeamConstraint::clamped(0)], &[]),
            Err(BeamSolverError::ZeroLength(0))
        ));
    }

    #[test]
    fn cantilever_first_natural_frequency_matches_analytic() {
        // First bending natural frequency of a clamped-free beam:
        //   f₁ = (β₁L)²/(2π)·√(E·I/(ρ·A·L⁴)),  (β₁L)² = 3.51602.
        // The beam consistent-mass element reproduces this closely.
        let mat = steel();
        let l = 5.0;
        let n_elem = 16;
        let nodes: Vec<Vector3<f64>> = (0..=n_elem)
            .map(|i| Vector3::new(l * i as f64 / n_elem as f64, 0.0, 0.0))
            .collect();
        let section = BeamSection::rectangle(0.05, 0.05);
        let elements: Vec<BeamElement> = (0..n_elem)
            .map(|i| BeamElement::new(i, i + 1, section))
            .collect();
        let constraints = [BeamConstraint::clamped(0)];
        let sol = solve_beam_modal(&nodes, &elements, &mat, &constraints, 4).unwrap();
        let f_fe = sol.fundamental_hz().unwrap();

        // The analytic cantilever-fundamental reference via this module's helper
        // (β₁L = 1.875104 for clamped–free), which the FE modal solver converges to.
        let f_analytic = euler_bernoulli_beam_frequency(
            1.875_104,
            mat.youngs_modulus,
            section.iy,
            mat.density,
            section.area,
            l,
        );
        let rel = (f_fe - f_analytic).abs() / f_analytic;
        assert!(
            rel < 0.03,
            "FE fundamental {f_fe} Hz vs analytic {f_analytic} Hz (rel {rel})"
        );
    }

    #[test]
    fn euler_bernoulli_beam_frequency_matches_the_closed_form_and_modes() {
        use std::f64::consts::PI;
        // A simply-supported steel beam (β₁L = π): f = π²/(2π·L²)·√(EI/(ρA)).
        let (e, i, rho, area, l) = (200.0e9, 1.0e-8, 7850.0, 1.0e-4, 1.0);
        let f_ss = euler_bernoulli_beam_frequency(PI, e, i, rho, area, l);
        let omega = PI.powi(2) / (l * l) * (e * i / (rho * area)).sqrt();
        assert!((f_ss - omega / (2.0 * PI)).abs() < 1e-9, "f = ω/2π");
        assert!((f_ss - 79.28).abs() < 0.5, "SS fundamental ≈ 79.3 Hz, got {f_ss}");
        // The boundary/mode enters only through (β·L)²: the cantilever fundamental
        // (β₁L = 1.875104) sits a factor (1.875104/π)² below the simply-supported one.
        let f_cant = euler_bernoulli_beam_frequency(1.875_104, e, i, rho, area, l);
        assert!(
            (f_cant / f_ss - (1.875_104_f64 / PI).powi(2)).abs() < 1e-9,
            "f ∝ (β·L)²"
        );
        // Scaling: ∝ √E, ∝ √I, ∝ 1/√ρ, ∝ 1/L².
        assert!((euler_bernoulli_beam_frequency(PI, 4.0 * e, i, rho, area, l) - 2.0 * f_ss).abs() < 1e-6, "∝ √E");
        assert!((euler_bernoulli_beam_frequency(PI, e, 4.0 * i, rho, area, l) - 2.0 * f_ss).abs() < 1e-6, "∝ √I");
        assert!((euler_bernoulli_beam_frequency(PI, e, i, 4.0 * rho, area, l) - f_ss / 2.0).abs() < 1e-6, "∝ 1/√ρ");
        assert!((euler_bernoulli_beam_frequency(PI, e, i, rho, area, 2.0 * l) - f_ss / 4.0).abs() < 1e-6, "∝ 1/L²");
        // Non-physical input → 0.
        assert_eq!(euler_bernoulli_beam_frequency(PI, 0.0, i, rho, area, l), 0.0);
        assert_eq!(euler_bernoulli_beam_frequency(PI, e, i, rho, area, 0.0), 0.0);
        assert_eq!(euler_bernoulli_beam_frequency(PI, e, i, 0.0, area, l), 0.0);
        assert_eq!(euler_bernoulli_beam_frequency(f64::NAN, e, i, rho, area, l), 0.0);
    }

    #[test]
    fn beam_modal_modes_are_ascending_and_positive() {
        let mat = steel();
        let l = 4.0;
        let n_elem = 10;
        let nodes: Vec<Vector3<f64>> = (0..=n_elem)
            .map(|i| Vector3::new(l * i as f64 / n_elem as f64, 0.0, 0.0))
            .collect();
        let section = BeamSection::circle(0.03);
        let elements: Vec<BeamElement> = (0..n_elem)
            .map(|i| BeamElement::new(i, i + 1, section))
            .collect();
        let sol = solve_beam_modal(&nodes, &elements, &mat, &[BeamConstraint::clamped(0)], 5)
            .unwrap();
        assert_eq!(sol.modes.len(), 5);
        for w in sol.modes.windows(2) {
            assert!(
                w[1].frequency_hz >= w[0].frequency_hz - 1e-6,
                "modes not ascending"
            );
        }
        assert!(sol.modes[0].frequency_hz > 0.0);
    }

    #[test]
    fn portal_frame_assembles_and_solves() {
        // A simple 2D portal frame in 3D space: two columns + a beam.
        //   node0 (base L) — node1 (top L) — node2 (top R) — node3 (base R)
        let mat = steel();
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0), // base left
            Vector3::new(0.0, 0.0, 3.0), // top left
            Vector3::new(4.0, 0.0, 3.0), // top right
            Vector3::new(4.0, 0.0, 0.0), // base right
        ];
        let section = BeamSection::rectangle(0.2, 0.3);
        let elements = [
            BeamElement::new(0, 1, section), // left column
            BeamElement::new(1, 2, section), // top beam
            BeamElement::new(2, 3, section), // right column
        ];
        let constraints = [BeamConstraint::clamped(0), BeamConstraint::clamped(3)];
        // A lateral load pushing the top of the frame sideways (+X).
        let loads = [BeamLoad::force(1, [1.0e4, 0.0, 0.0])];
        let sol = solve_beam_static(&nodes, &elements, &mat, &constraints, &loads).unwrap();
        // The loaded top corner must sway in +X; the clamped bases stay.
        assert!(sol.translation[1][0] > 0.0, "frame should sway +X");
        assert!(sol.translation[0][0].abs() < 1e-9, "clamped base moved");
        assert!(sol.translation[3][0].abs() < 1e-9, "clamped base moved");
        // The top beam is stiff axially, so node 2 sways a similar
        // amount to node 1.
        assert!(
            (sol.translation[2][0] - sol.translation[1][0]).abs()
                < 0.5 * sol.translation[1][0].abs(),
            "portal top beam should carry the sway across"
        );
    }

    // ----- Round-1 H4: non-finite beam load / BC rejection -------------

    /// A minimal valid cantilever: two nodes along X, one element,
    /// node 0 clamped. Shared by the H4 rejection tests below.
    fn simple_cantilever() -> (Vec<Vector3<f64>>, Vec<BeamElement>, Vec<BeamConstraint>) {
        let nodes = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        let elements = vec![BeamElement::new(0, 1, BeamSection::rectangle(0.1, 0.1))];
        let constraints = vec![BeamConstraint::clamped(0)];
        (nodes, elements, constraints)
    }

    #[test]
    fn beam_rejects_nan_force() {
        // A NaN beam force pushed straight into the RHS would yield a
        // silently-NaN displacement returned as Ok(..). Reject it.
        let (nodes, elements, constraints) = simple_cantilever();
        let loads = vec![BeamLoad::force(1, [f64::NAN, 0.0, 0.0])];
        let err = solve_beam_static(&nodes, &elements, &steel(), &constraints, &loads)
            .unwrap_err();
        assert!(
            matches!(err, BeamSolverError::InvalidLoad { .. }),
            "expected InvalidLoad, got {err:?}"
        );
    }

    #[test]
    fn beam_rejects_infinite_moment() {
        let (nodes, elements, constraints) = simple_cantilever();
        let loads = vec![BeamLoad::moment(1, [0.0, f64::INFINITY, 0.0])];
        let err = solve_beam_static(&nodes, &elements, &steel(), &constraints, &loads)
            .unwrap_err();
        assert!(matches!(err, BeamSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn beam_rejects_non_finite_prescribed_displacement() {
        // A prescribed constraint value is folded into the RHS as
        // `penalty·value`; a non-finite value corrupts the solve.
        let (nodes, elements, _c) = simple_cantilever();
        let constraints = vec![
            BeamConstraint::clamped(0),
            BeamConstraint {
                node: 1,
                fixed: [Some(f64::NAN), None, None, None, None, None],
            },
        ];
        let err = solve_beam_static(&nodes, &elements, &steel(), &constraints, &[])
            .unwrap_err();
        assert!(matches!(err, BeamSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn beam_still_accepts_finite_loads() {
        // The validation must not reject a normal finite load.
        let (nodes, elements, constraints) = simple_cantilever();
        let loads = vec![BeamLoad::force(1, [0.0, 0.0, -1.0e3])];
        let sol = solve_beam_static(&nodes, &elements, &steel(), &constraints, &loads)
            .unwrap();
        assert!(sol.translation.iter().all(|t| t.iter().all(|c| c.is_finite())));
        assert!(sol.translation[1][2] < 0.0, "tip should deflect downward");
    }

    // ----- Round-2 F1: dense beam allocation cap (6 DOF/node) ----------

    #[test]
    fn check_beam_dense_dofs_uses_six_dof_per_node() {
        use crate::native_solver::MAX_DENSE_DOFS;
        // 6·n_nodes, returned untouched up to the cap.
        assert_eq!(check_beam_dense_dofs(2).unwrap(), 12);
        let at = MAX_DENSE_DOFS / 6;
        assert!(check_beam_dense_dofs(at).is_ok());
        // The motivating regression: a node count whose 3·n is under the
        // continuum cap but whose 6·n is over it must be rejected by the
        // beam path — pure arithmetic, instant, no allocation.
        let n = MAX_DENSE_DOFS / 6 + 1;
        let err = check_beam_dense_dofs(n).unwrap_err();
        match err {
            BeamSolverError::TooLarge { dofs, max } => {
                assert_eq!(dofs, n * 6);
                assert_eq!(max, MAX_DENSE_DOFS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
        // A 6·n that overflows usize surfaces as TooLarge, never a wrap.
        assert!(matches!(
            check_beam_dense_dofs(usize::MAX).unwrap_err(),
            BeamSolverError::TooLarge { .. }
        ));
    }

    /// A long straight chain of `n_nodes` collinear nodes with a beam
    /// element between each adjacent pair, node 0 clamped. Cheap to
    /// build (the `O(n_dof²)` cost is the dense matrix the cap prevents).
    fn straight_chain(n_nodes: usize) -> (Vec<Vector3<f64>>, Vec<BeamElement>, Vec<BeamConstraint>) {
        let nodes: Vec<Vector3<f64>> = (0..n_nodes)
            .map(|i| Vector3::new(i as f64, 0.0, 0.0))
            .collect();
        let elements: Vec<BeamElement> = (0..n_nodes.saturating_sub(1))
            .map(|i| BeamElement::new(i, i + 1, BeamSection::rectangle(0.1, 0.1)))
            .collect();
        let constraints = vec![BeamConstraint::clamped(0)];
        (nodes, elements, constraints)
    }

    #[test]
    fn solve_beam_static_rejects_oversized_frame_without_allocating() {
        use crate::native_solver::MAX_DENSE_DOFS;
        // 6·n_nodes just over the cap. The cap fires before the
        // `6·n_nodes × 6·n_nodes` (~16002²·8 B ≈ 2 GB) matrix is touched;
        // if it did not, this test would OOM rather than return Err.
        let n = MAX_DENSE_DOFS / 6 + 1;
        let (nodes, elements, constraints) = straight_chain(n);
        let err = solve_beam_static(&nodes, &elements, &steel(), &constraints, &[])
            .unwrap_err();
        assert!(
            matches!(err, BeamSolverError::TooLarge { .. }),
            "expected TooLarge, got {err:?}"
        );
    }

    #[test]
    fn solve_beam_modal_rejects_oversized_frame_without_allocating() {
        use crate::native_solver::MAX_DENSE_DOFS;
        let n = MAX_DENSE_DOFS / 6 + 1;
        let (nodes, elements, constraints) = straight_chain(n);
        let err = solve_beam_modal(&nodes, &elements, &steel(), &constraints, 2)
            .unwrap_err();
        assert!(
            matches!(err, BeamSolverError::TooLarge { .. }),
            "expected TooLarge, got {err:?}"
        );
    }

    #[test]
    fn small_beam_still_solves_after_cap() {
        // A normal small frame is comfortably under the cap and solves
        // exactly as before — the cap must not perturb valid inputs.
        let (nodes, elements, constraints) = simple_cantilever();
        let loads = vec![BeamLoad::force(1, [0.0, 0.0, -1.0e3])];
        let sol = solve_beam_static(&nodes, &elements, &steel(), &constraints, &loads)
            .expect("small frame must still solve");
        assert!(sol.translation[1][2] < 0.0);
        // And the modal path likewise.
        let modal = solve_beam_modal(&nodes, &elements, &steel(), &constraints, 1)
            .expect("small modal must still solve");
        assert_eq!(modal.modes.len(), 1);
    }
}
