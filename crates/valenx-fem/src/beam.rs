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
        let analytic = f * l / (mat.youngs_modulus * section.area);
        let rel = (sol.translation[1][0] - analytic).abs() / analytic;
        assert!(rel < 1e-6, "axial δ {} vs {analytic}", sol.translation[1][0]);
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
        let analytic = torque * l / (g * section.j);
        let rel = (sol.rotation[1][0] - analytic).abs() / analytic;
        assert!(rel < 1e-6, "twist {} vs analytic {analytic}", sol.rotation[1][0]);
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
        let bending = p * l.powi(3) / (48.0 * mat.youngs_modulus * i);
        // Allow a few % for the Timoshenko shear contribution.
        let centre = sol.translation[mid][2].abs();
        let rel = (centre - bending).abs() / bending;
        assert!(
            rel < 0.05,
            "mid-span deflection {centre} vs Euler-Bernoulli {bending} (rel {rel})"
        );
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
