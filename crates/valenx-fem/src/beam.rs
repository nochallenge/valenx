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
