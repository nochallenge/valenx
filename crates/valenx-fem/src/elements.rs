//! Generic **continuum solid element library** for the native FEA
//! solvers (Phase 24.8).
//!
//! ## What this is
//!
//! The original native solvers ([`crate::native_solver`],
//! [`crate::modal_solver`], …) shipped a single element — the 4-node
//! constant-strain tetrahedron ([`ElementType::Tet4`]). Commercial FEA
//! (ANSYS, Abaqus) carries dozens of element types because no single
//! element is good everywhere: the linear tet is cheap but notoriously
//! **over-stiff in bending** ("shear / volumetric locking"), so a
//! coarse Tet4 mesh badly under-predicts the deflection of a beam.
//!
//! This module is the element-library layer. It defines a common
//! [`SolidElement`] trait — every 3-translational-DOF-per-node
//! continuum element implements it — and ships three implementations:
//!
//! - [`Tet4`] — the 4-node linear tetrahedron (the original element,
//!   re-expressed against the trait so the generic assembly can consume
//!   it).
//! - [`Hex8`] — the 8-node trilinear hexahedron ("brick"), the
//!   workhorse of commercial solid FEA. 2×2×2 Gauss quadrature.
//! - [`Tet10`] — the 10-node **quadratic** tetrahedron. A quadratic
//!   displacement field captures bending far better per DOF than the
//!   linear tet, so it converges to beam theory dramatically faster.
//!
//! All three are **isoparametric**: the same shape functions
//! interpolate the geometry and the displacement field, and the element
//! integrals are evaluated by Gauss-Legendre quadrature on a reference
//! element. The [`SolidElement`] trait exposes
//! [`SolidElement::stiffness`] (`Kₑ = ∫ Bᵀ D B dV`),
//! [`SolidElement::consistent_mass`] (`Mₑ = ∫ ρ Nᵀ N dV`) and
//! [`SolidElement::recover_stress`] (`σ = D B uₑ` at the element
//! centroid), all assembled by the *same* Gauss loop.
//!
//! ## Honest scope
//!
//! These are **isotropic linear-elastic small-strain** continuum
//! elements — exactly the regime the existing solvers target. They are
//! real, validated elements (each passes the constant-strain **patch
//! test** to solver precision — see the tests and
//! [`crate::validation`]); they are not the full commercial element
//! zoo (no reduced-integration / hourglass-stabilised bricks, no
//! incompatible-modes elements, no anisotropy). The 3D structural
//! **beam** element lives in [`crate::beam`] because it carries
//! rotational DOFs and so does not fit the 3-DOF-per-node continuum
//! trait.

use nalgebra::{DMatrix, DVector, Matrix3, Matrix6, Vector3};

use valenx_mesh::element::ElementType;

/// 6×6 isotropic elasticity matrix `D` in Voigt notation
/// `[σxx σyy σzz σxy σyz σzx] = D·[εxx εyy εzz γxy γyz γzx]`.
///
/// Uses engineering shear strain (`γ = 2ε`); the lower-right block is
/// therefore the shear modulus `G` (not `2G`). This is the same
/// formulation the `native_solver`'s `elasticity_matrix` helper uses —
/// the element library takes the pre-built matrix so it stays free of
/// the material type.
pub type ElasticityMatrix = Matrix6<f64>;

/// A single Gauss-Legendre quadrature point: a position in the
/// element's reference (natural) coordinates plus its integration
/// weight.
#[derive(Copy, Clone, Debug)]
pub struct GaussPoint {
    /// Natural coordinates of the point. The meaning of the three
    /// components depends on the element family (`ξ,η,ζ ∈ [-1,1]` for
    /// the hex; volume / area coordinates for the simplex elements).
    pub coords: [f64; 3],
    /// Integration weight.
    pub weight: f64,
}

/// Evaluated shape-function data at one quadrature point.
///
/// Produced by [`SolidElement::shape`]: the shape-function values `N`
/// and their gradients with respect to the *physical* coordinates,
/// together with the Jacobian determinant that turns a reference-space
/// integral into a physical-space one.
#[derive(Clone, Debug)]
pub struct ShapeEval {
    /// Shape-function values `Nₐ` at the point, one per element node.
    pub n: Vec<f64>,
    /// Physical gradients `∂Nₐ/∂x`, `∂Nₐ/∂y`, `∂Nₐ/∂z`, one
    /// [`Vector3`] per node.
    pub grad: Vec<Vector3<f64>>,
    /// Determinant of the isoparametric Jacobian `∂x/∂ξ` at the point.
    /// The physical volume element is `det(J)·dξdηdζ`.
    pub detj: f64,
}

/// A 3-translational-DOF-per-node **continuum solid element**.
///
/// Every implementor (Tet4 / Hex8 / Tet10) is isoparametric: it knows
/// its reference shape functions and a Gauss quadrature rule, and the
/// trait's default methods assemble the element stiffness, consistent
/// mass and recovered stress from those two pieces. An element with
/// `n` nodes has `3n` DOFs ordered node-major: local DOF `3a+i` is
/// component `i ∈ {x,y,z}` of node `a`.
pub trait SolidElement {
    /// The mesh [`ElementType`] this element implements.
    fn element_type(&self) -> ElementType;

    /// Number of nodes (and hence `3·n` DOFs).
    fn n_nodes(&self) -> usize;

    /// Node coordinates in physical space, in connectivity order.
    fn node_coords(&self) -> &[Vector3<f64>];

    /// The Gauss-Legendre quadrature rule used to integrate this
    /// element's **stiffness** `Kₑ = ∫ BᵀDB dV`.
    fn quadrature(&self) -> Vec<GaussPoint>;

    /// The quadrature rule used to integrate this element's
    /// **consistent mass** `Mₑ = ∫ ρNᵀN dV`.
    ///
    /// Defaults to [`SolidElement::quadrature`] — correct for the
    /// hexahedral and quadratic-tet elements, whose stiffness rule
    /// already integrates the mass integrand exactly. The
    /// constant-strain [`Tet4`] overrides it: its stiffness integrand
    /// is *constant* so a single point suffices there, but the mass
    /// integrand `NᵢNⱼ` is *quadratic*, so a one-point rule would give
    /// a rank-deficient (non-positive-definite) mass matrix. Tet4
    /// therefore uses a 4-point rule for the mass.
    fn mass_quadrature(&self) -> Vec<GaussPoint> {
        self.quadrature()
    }

    /// Evaluate the shape functions and their *physical* gradients at a
    /// natural-coordinate point. Returns `None` if the isoparametric
    /// Jacobian is singular (a degenerate / inside-out element).
    fn shape(&self, natural: [f64; 3]) -> Option<ShapeEval>;

    /// The 6×(3n) strain-displacement matrix `B` at a pre-evaluated
    /// quadrature point — maps the `3n` nodal displacement components
    /// to the 6 Voigt strain components. The default builds the
    /// standard layout from the shape gradients; no element needs to
    /// override it.
    fn strain_displacement(&self, shape: &ShapeEval) -> DMatrix<f64> {
        let n = self.n_nodes();
        let mut b = DMatrix::<f64>::zeros(6, 3 * n);
        for (a, g) in shape.grad.iter().enumerate() {
            let (bx, by, bz) = (g.x, g.y, g.z);
            let col = 3 * a;
            // εxx = ∂u/∂x
            b[(0, col)] = bx;
            // εyy = ∂v/∂y
            b[(1, col + 1)] = by;
            // εzz = ∂w/∂z
            b[(2, col + 2)] = bz;
            // γxy = ∂u/∂y + ∂v/∂x
            b[(3, col)] = by;
            b[(3, col + 1)] = bx;
            // γyz = ∂v/∂z + ∂w/∂y
            b[(4, col + 1)] = bz;
            b[(4, col + 2)] = by;
            // γzx = ∂w/∂x + ∂u/∂z
            b[(5, col)] = bz;
            b[(5, col + 2)] = bx;
        }
        b
    }

    /// The `3n × 3n` element stiffness `Kₑ = ∫ Bᵀ·D·B dV`, integrated
    /// by the element's Gauss rule. `None` for a degenerate element.
    fn stiffness(&self, d: &ElasticityMatrix) -> Option<DMatrix<f64>> {
        let ndof = 3 * self.n_nodes();
        let mut ke = DMatrix::<f64>::zeros(ndof, ndof);
        let mut any = false;
        for gp in self.quadrature() {
            let sh = self.shape(gp.coords)?;
            if !(sh.detj.is_finite()) || sh.detj <= 0.0 {
                return None;
            }
            let b = self.strain_displacement(&sh);
            let bt_d = b.transpose() * d;
            ke += (bt_d * &b) * (sh.detj * gp.weight);
            any = true;
        }
        if !any {
            return None;
        }
        Some(ke)
    }

    /// The `3n × 3n` **consistent** element mass `Mₑ = ∫ ρ·Nᵀ·N dV`,
    /// integrated by the element's Gauss rule. The three translational
    /// DOFs of a node are uncoupled, so each `(a,b)` node pair
    /// contributes `(∫ρ Nₐ Nᵦ dV)·I₃`. `None` for a degenerate element.
    fn consistent_mass(&self, density: f64) -> Option<DMatrix<f64>> {
        let n = self.n_nodes();
        let ndof = 3 * n;
        let mut me = DMatrix::<f64>::zeros(ndof, ndof);
        let mut any = false;
        for gp in self.mass_quadrature() {
            let sh = self.shape(gp.coords)?;
            if !(sh.detj.is_finite()) || sh.detj <= 0.0 {
                return None;
            }
            let scale = density * sh.detj * gp.weight;
            for a in 0..n {
                for b in 0..n {
                    let coeff = scale * sh.n[a] * sh.n[b];
                    for i in 0..3 {
                        me[(3 * a + i, 3 * b + i)] += coeff;
                    }
                }
            }
            any = true;
        }
        if !any {
            return None;
        }
        Some(me)
    }

    /// The element volume `∫ dV`, integrated by the Gauss rule. Useful
    /// for checking a mesh and for lumped-mass fallbacks. `None` for a
    /// degenerate element.
    fn volume(&self) -> Option<f64> {
        let mut v = 0.0;
        let mut any = false;
        for gp in self.quadrature() {
            let sh = self.shape(gp.coords)?;
            v += sh.detj * gp.weight;
            any = true;
        }
        if !any {
            return None;
        }
        Some(v)
    }

    /// Recover the Voigt stress `σ = D·B·uₑ` evaluated at the element
    /// **centroid**, given the `3n`-long element displacement vector.
    /// `None` for a degenerate element.
    fn recover_stress(&self, d: &ElasticityMatrix, ue: &DVector<f64>) -> Option<[f64; 6]> {
        let sh = self.shape(self.centroid_natural())?;
        let b = self.strain_displacement(&sh);
        let strain = &b * ue;
        let sigma = d * strain;
        Some([sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5]])
    }

    /// The element centroid in natural coordinates — where
    /// [`SolidElement::recover_stress`] samples the stress.
    fn centroid_natural(&self) -> [f64; 3];
}

// ===========================================================================
// Shared helpers
// ===========================================================================

/// One point of a tetrahedron quadrature orbit of the type
/// `(a, b, b, b)` — one barycentric coordinate is `a`, the other three
/// are `b`. `which` selects which of the four barycentric slots holds
/// `a`. Returns the natural coordinates `(ξ,η,ζ) = (L1,L2,L3)`.
fn tet_orbit_abbb(which: usize, a: f64, b: f64) -> [f64; 3] {
    let mut l = [b; 4];
    l[which] = a;
    [l[1], l[2], l[3]]
}

/// One point of a tetrahedron quadrature orbit of the type
/// `(a, a, b, b)` — two barycentric coordinates are `a`, the other two
/// are `b`. `edge` (0..6) selects which pair of slots holds `a`,
/// following the [`TET10_EDGES`] pairing.
fn tet_orbit_aabb(edge: usize, a: f64, b: f64) -> [f64; 3] {
    // The six unordered pairs of {0,1,2,3}.
    const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let (p, q) = PAIRS[edge];
    let mut l = [b; 4];
    l[p] = a;
    l[q] = a;
    [l[1], l[2], l[3]]
}

/// Invert a 3×3 isoparametric Jacobian, returning `(J⁻¹, det J)`.
/// `None` if the determinant is below a tiny tolerance.
fn invert_jacobian(j: &Matrix3<f64>) -> Option<(Matrix3<f64>, f64)> {
    let det = j.determinant();
    if !det.is_finite() || det.abs() < 1.0e-20 {
        return None;
    }
    let inv = j.try_inverse()?;
    Some((inv, det))
}

/// Assemble the 3×3 isoparametric Jacobian `J = Σ xₐ ⊗ ∂Nₐ/∂ξ`.
///
/// The columns of `J` are `∂x/∂ξ`, `∂x/∂η`, `∂x/∂ζ`. `coords` and
/// `reference_grad` must be the same length (one entry per element
/// node). Shared by every isoparametric element so the Jacobian loop is
/// written once.
fn assemble_jacobian(coords: &[Vector3<f64>], reference_grad: &[Vector3<f64>]) -> Matrix3<f64> {
    let mut j = Matrix3::zeros();
    for (x, g) in coords.iter().zip(reference_grad.iter()) {
        for c in 0..3 {
            j[(0, c)] += x.x * g[c];
            j[(1, c)] += x.y * g[c];
            j[(2, c)] += x.z * g[c];
        }
    }
    j
}

// ===========================================================================
// Tet4 — 4-node linear tetrahedron
// ===========================================================================

/// The 4-node **linear tetrahedron** (`C3D4`), expressed against the
/// [`SolidElement`] trait.
///
/// The four shape functions are the barycentric (volume) coordinates
/// `N₀ = 1−ξ−η−ζ`, `N₁ = ξ`, `N₂ = η`, `N₃ = ζ`; their gradients are
/// constant over the element, so the strain field is constant
/// ("constant-strain tetrahedron"). A single-point quadrature
/// integrates the constant integrand exactly. This is the original
/// native-solver element — re-expressed here so the generic
/// mixed-element assembly can consume a Tet4 the same way it consumes a
/// Hex8 or a Tet10.
#[derive(Clone, Debug)]
pub struct Tet4 {
    coords: [Vector3<f64>; 4],
}

impl Tet4 {
    /// Build a Tet4 from its four physical node coordinates.
    pub fn new(coords: [Vector3<f64>; 4]) -> Self {
        Self { coords }
    }
}

impl SolidElement for Tet4 {
    fn element_type(&self) -> ElementType {
        ElementType::Tet4
    }
    fn n_nodes(&self) -> usize {
        4
    }
    fn node_coords(&self) -> &[Vector3<f64>] {
        &self.coords
    }
    fn centroid_natural(&self) -> [f64; 3] {
        [0.25, 0.25, 0.25]
    }
    fn quadrature(&self) -> Vec<GaussPoint> {
        // The integrand (Bᵀ D B) is constant over a linear tet, so the
        // exact integral is value × volume. The reference tet has
        // volume 1/6; one centroid point with weight 1/6 integrates it.
        vec![GaussPoint {
            coords: [0.25, 0.25, 0.25],
            weight: 1.0 / 6.0,
        }]
    }
    fn mass_quadrature(&self) -> Vec<GaussPoint> {
        // The consistent-mass integrand NᵢNⱼ is quadratic, so a single
        // centroid point gives a rank-1 (non-SPD) mass. The degree-2
        // 4-point tet rule integrates it exactly → the classical
        // ρV/20·(2 on diag, 1 off) consistent mass.
        let a = 0.585_410_196_624_968_5;
        let b = 0.138_196_601_125_010_5;
        let w = (1.0 / 6.0) / 4.0;
        vec![
            GaussPoint {
                coords: [a, b, b],
                weight: w,
            },
            GaussPoint {
                coords: [b, a, b],
                weight: w,
            },
            GaussPoint {
                coords: [b, b, a],
                weight: w,
            },
            GaussPoint {
                coords: [b, b, b],
                weight: w,
            },
        ]
    }
    fn shape(&self, natural: [f64; 3]) -> Option<ShapeEval> {
        let (xi, eta, zeta) = (natural[0], natural[1], natural[2]);
        let n = vec![1.0 - xi - eta - zeta, xi, eta, zeta];
        // Reference-coordinate gradients ∂N/∂(ξ,η,ζ) — constant.
        let dn = [
            Vector3::new(-1.0, -1.0, -1.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        // Jacobian J = Σ xₐ ⊗ ∂Nₐ/∂ξ  (columns are ∂x/∂ξ, ∂x/∂η, ∂x/∂ζ).
        let j = assemble_jacobian(&self.coords, &dn);
        let (jinv, detj) = invert_jacobian(&j)?;
        // Physical gradient ∂Nₐ/∂x = J⁻ᵀ · ∂Nₐ/∂ξ.
        let jinv_t = jinv.transpose();
        let grad = dn.iter().map(|g| jinv_t * g).collect();
        Some(ShapeEval { n, grad, detj })
    }
}

// ===========================================================================
// Hex8 — 8-node trilinear hexahedron
// ===========================================================================

/// The 8-node **trilinear hexahedron** (`C3D8`), the workhorse solid
/// element of commercial FEA.
///
/// The reference element is the cube `[-1,1]³`; the eight shape
/// functions are the trilinear products
/// `Nₐ = ⅛(1+ξₐξ)(1+ηₐη)(1+ζₐζ)` where `(ξₐ,ηₐ,ζₐ)` is corner `a`'s
/// reference position. The element is integrated by **2×2×2 Gauss-
/// Legendre quadrature** (8 points at `±1/√3`) — full integration,
/// exact for the trilinear stiffness integrand on a parallelepiped.
///
/// Node ordering follows the standard hex convention used by
/// [`crate::native_solver::structured_box_mesh`]'s cells:
/// `0..4` are the `ζ=-1` face counter-clockwise, `4..8` the `ζ=+1`
/// face. A Hex8 is far less stiff in bending than the constant-strain
/// Tet4: a single layer of bricks already captures a linear bending
/// stress, so a Hex8 mesh converges to beam theory with far fewer DOFs.
#[derive(Clone, Debug)]
pub struct Hex8 {
    coords: [Vector3<f64>; 8],
}

/// The eight reference-cube corner positions, in standard hex node
/// order — `0..4` is the bottom (`ζ=-1`) face CCW, `4..8` the top.
const HEX8_CORNERS: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];

impl Hex8 {
    /// Build a Hex8 from its eight physical node coordinates, in the
    /// standard hex node order — the `ζ=-1` face counter-clockwise
    /// (`0..4`) then the `ζ=+1` face (`4..8`).
    pub fn new(coords: [Vector3<f64>; 8]) -> Self {
        Self { coords }
    }
}

impl SolidElement for Hex8 {
    fn element_type(&self) -> ElementType {
        ElementType::Hex8
    }
    fn n_nodes(&self) -> usize {
        8
    }
    fn node_coords(&self) -> &[Vector3<f64>] {
        &self.coords
    }
    fn centroid_natural(&self) -> [f64; 3] {
        [0.0, 0.0, 0.0]
    }
    fn quadrature(&self) -> Vec<GaussPoint> {
        // 2×2×2 Gauss-Legendre: points at ±1/√3, every weight 1.
        let g = 1.0 / 3.0_f64.sqrt();
        let mut pts = Vec::with_capacity(8);
        for &z in &[-g, g] {
            for &y in &[-g, g] {
                for &x in &[-g, g] {
                    pts.push(GaussPoint {
                        coords: [x, y, z],
                        weight: 1.0,
                    });
                }
            }
        }
        pts
    }
    fn shape(&self, natural: [f64; 3]) -> Option<ShapeEval> {
        let (xi, eta, zeta) = (natural[0], natural[1], natural[2]);
        let mut n = vec![0.0; 8];
        // Reference gradients ∂N/∂(ξ,η,ζ).
        let mut dn = [Vector3::zeros(); 8];
        for (a, corner) in HEX8_CORNERS.iter().enumerate() {
            let (xa, ya, za) = (corner[0], corner[1], corner[2]);
            n[a] = 0.125 * (1.0 + xa * xi) * (1.0 + ya * eta) * (1.0 + za * zeta);
            dn[a] = Vector3::new(
                0.125 * xa * (1.0 + ya * eta) * (1.0 + za * zeta),
                0.125 * ya * (1.0 + xa * xi) * (1.0 + za * zeta),
                0.125 * za * (1.0 + xa * xi) * (1.0 + ya * eta),
            );
        }
        let j = assemble_jacobian(&self.coords, &dn);
        let (jinv, detj) = invert_jacobian(&j)?;
        let jinv_t = jinv.transpose();
        let grad = dn.iter().map(|g| jinv_t * g).collect();
        Some(ShapeEval { n, grad, detj })
    }
}

// ===========================================================================
// Tet10 — 10-node quadratic tetrahedron
// ===========================================================================

/// The 10-node **quadratic tetrahedron** (`C3D10`).
///
/// Four corner nodes (`0..4`) plus six mid-edge nodes (`4..10`) carry a
/// **quadratic** displacement field, so the strain varies *linearly*
/// over the element — a single Tet10 already represents the linear
/// bending stress that defeats the constant-strain Tet4. Edge ordering
/// for the mid-side nodes is the standard `C3D10` convention:
///
/// ```text
///   node 4 : edge 0-1     node 7 : edge 0-3
///   node 5 : edge 1-2     node 8 : edge 1-3
///   node 6 : edge 2-0     node 9 : edge 2-3
/// ```
///
/// Integrated by the **4-point** Gauss rule for a tetrahedron (exact
/// for the quadratic-times-quadratic stiffness integrand on a
/// straight-edged tet).
#[derive(Clone, Debug)]
pub struct Tet10 {
    coords: [Vector3<f64>; 10],
}

/// Mid-side node → its two end corner indices, in standard `C3D10`
/// edge order.
pub const TET10_EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];

impl Tet10 {
    /// Build a Tet10 from its ten physical node coordinates — four
    /// corners then six mid-edge nodes (see [`TET10_EDGES`]).
    pub fn new(coords: [Vector3<f64>; 10]) -> Self {
        Self { coords }
    }
}

impl SolidElement for Tet10 {
    fn element_type(&self) -> ElementType {
        ElementType::Tet10
    }
    fn n_nodes(&self) -> usize {
        10
    }
    fn node_coords(&self) -> &[Vector3<f64>] {
        &self.coords
    }
    fn centroid_natural(&self) -> [f64; 3] {
        [0.25, 0.25, 0.25]
    }
    fn quadrature(&self) -> Vec<GaussPoint> {
        // 4-point Gauss rule for a tetrahedron — degree-2 exact.
        // The Tet10 strain field is linear, so BᵀDB is quadratic and
        // this rule integrates the stiffness exactly. Points at the
        // (a,b,b,b) permutations with a = 0.5854102, b = 0.1381966;
        // each weight = (1/6)/4 (they sum to the reference-tet
        // volume 1/6).
        let a = 0.585_410_196_624_968_5;
        let b = 0.138_196_601_125_010_5;
        let w = (1.0 / 6.0) / 4.0;
        vec![
            GaussPoint {
                coords: [a, b, b],
                weight: w,
            },
            GaussPoint {
                coords: [b, a, b],
                weight: w,
            },
            GaussPoint {
                coords: [b, b, a],
                weight: w,
            },
            GaussPoint {
                coords: [b, b, b],
                weight: w,
            },
        ]
    }
    fn mass_quadrature(&self) -> Vec<GaussPoint> {
        // The Tet10 consistent-mass integrand NᵢNⱼ is *quartic*, so the
        // degree-2 stiffness rule under-integrates it (rank-deficient
        // mass). Keast's degree-4, 15-point tetrahedron rule integrates
        // it exactly. The weights below are quoted for the reference
        // tetrahedron of volume 1/6 — they sum to 1/6 directly.
        let mut pts: Vec<GaussPoint> = Vec::with_capacity(15);
        // Centroid.
        pts.push(GaussPoint {
            coords: [0.25, 0.25, 0.25],
            weight: 0.030_283_678_097_089,
        });
        // Orbit at (a, b, b, b) — 4 permutations.
        {
            let a = 0.0;
            let b = 1.0 / 3.0;
            let w = 0.006_026_785_714_286;
            for corner in 0..4 {
                pts.push(GaussPoint {
                    coords: tet_orbit_abbb(corner, a, b),
                    weight: w,
                });
            }
        }
        // Orbit at (a, b, b, b) with a = 8/11 — 4 permutations.
        {
            let a = 8.0 / 11.0;
            let b = 1.0 / 11.0;
            let w = 0.011_645_249_086_029;
            for corner in 0..4 {
                pts.push(GaussPoint {
                    coords: tet_orbit_abbb(corner, a, b),
                    weight: w,
                });
            }
        }
        // Orbit at (a, a, b, b) — 6 permutations (the edge orbit).
        {
            let a = 0.066_550_153_573_664;
            let b = 0.433_449_846_426_336;
            let w = 0.010_949_141_561_386;
            for edge in 0..6 {
                pts.push(GaussPoint {
                    coords: tet_orbit_aabb(edge, a, b),
                    weight: w,
                });
            }
        }
        pts
    }
    fn shape(&self, natural: [f64; 3]) -> Option<ShapeEval> {
        // Volume coordinates: L0 = 1−ξ−η−ζ, L1 = ξ, L2 = η, L3 = ζ.
        let (xi, eta, zeta) = (natural[0], natural[1], natural[2]);
        let l0 = 1.0 - xi - eta - zeta;
        let l = [l0, xi, eta, zeta];
        // Quadratic shape functions.
        //   corner a:  Nₐ = Lₐ(2Lₐ−1)
        //   mid-edge between corners (p,q): N = 4·Lₚ·L_q
        let mut n = vec![0.0; 10];
        for a in 0..4 {
            n[a] = l[a] * (2.0 * l[a] - 1.0);
        }
        for (m, &(p, q)) in TET10_EDGES.iter().enumerate() {
            n[4 + m] = 4.0 * l[p] * l[q];
        }
        // Reference gradients. dL0/dξ = (-1,-1,-1); dL1 = e_x; dL2 =
        // e_y; dL3 = e_z. Chain-rule each shape function.
        let dl = [
            Vector3::new(-1.0, -1.0, -1.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut dn = [Vector3::zeros(); 10];
        for a in 0..4 {
            // d/dξ [ Lₐ(2Lₐ−1) ] = (4Lₐ−1)·dLₐ
            dn[a] = dl[a] * (4.0 * l[a] - 1.0);
        }
        for (m, &(p, q)) in TET10_EDGES.iter().enumerate() {
            // d/dξ [ 4 Lₚ L_q ] = 4(L_q·dLₚ + Lₚ·dL_q)
            dn[4 + m] = (dl[p] * l[q] + dl[q] * l[p]) * 4.0;
        }
        let j = assemble_jacobian(&self.coords, &dn);
        let (jinv, detj) = invert_jacobian(&j)?;
        let jinv_t = jinv.transpose();
        let grad = dn.iter().map(|g| jinv_t * g).collect();
        Some(ShapeEval { n, grad, detj })
    }
}

/// von Mises equivalent stress from a Voigt stress vector
/// `[σxx σyy σzz σxy σyz σzx]`. A convenience identical to the
/// `native_solver`'s `von_mises_from_voigt` helper.
pub fn von_mises(s: &[f64; 6]) -> f64 {
    let (sx, sy, sz, sxy, syz, szx) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2))
        + 3.0 * (sxy * sxy + syz * syz + szx * szx))
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A reference unit tet (corners at the origin + the three axes).
    fn ref_tet4() -> Tet4 {
        Tet4::new([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ])
    }

    /// A 2×3×4 axis-aligned brick.
    fn brick_hex8(lx: f64, ly: f64, lz: f64) -> Hex8 {
        Hex8::new([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ])
    }

    /// A straight-edged unit Tet10: corners + the six edge midpoints.
    fn ref_tet10() -> Tet10 {
        let c = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut all = [Vector3::zeros(); 10];
        all[..4].copy_from_slice(&c);
        for (m, &(p, q)) in TET10_EDGES.iter().enumerate() {
            all[4 + m] = (c[p] + c[q]) * 0.5;
        }
        Tet10::new(all)
    }

    /// A simple isotropic D matrix for E = 1, ν = 0 (so it is trivial
    /// to reason about): D = diag(1,1,1, 0.5,0.5,0.5).
    fn d_simple() -> ElasticityMatrix {
        let mut d = Matrix6::zeros();
        for i in 0..3 {
            d[(i, i)] = 1.0;
        }
        for i in 3..6 {
            d[(i, i)] = 0.5;
        }
        d
    }

    #[test]
    fn tet4_shape_functions_sum_to_one() {
        let t = ref_tet4();
        for p in [[0.25, 0.25, 0.25], [0.1, 0.2, 0.3], [0.0, 0.0, 0.0]] {
            let sh = t.shape(p).unwrap();
            let sum: f64 = sh.n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "ΣN = {sum} ≠ 1");
        }
    }

    #[test]
    fn hex8_shape_functions_sum_to_one_and_partition() {
        let h = brick_hex8(2.0, 3.0, 4.0);
        for p in [[0.0, 0.0, 0.0], [0.3, -0.7, 0.5], [-1.0, -1.0, -1.0]] {
            let sh = h.shape(p).unwrap();
            let sum: f64 = sh.n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "ΣN = {sum} ≠ 1");
            // Shape gradients must sum to zero (a constant field has
            // zero gradient — the partition-of-unity derivative).
            let gsum: Vector3<f64> = sh.grad.iter().sum();
            assert!(gsum.norm() < 1e-10, "Σ∇N = {gsum} ≠ 0");
        }
    }

    #[test]
    fn tet10_shape_functions_sum_to_one_and_partition() {
        let t = ref_tet10();
        for p in [[0.25, 0.25, 0.25], [0.1, 0.2, 0.3], [0.5, 0.0, 0.0]] {
            let sh = t.shape(p).unwrap();
            let sum: f64 = sh.n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "ΣN = {sum} ≠ 1");
            let gsum: Vector3<f64> = sh.grad.iter().sum();
            assert!(gsum.norm() < 1e-10, "Σ∇N = {gsum} ≠ 0");
        }
    }

    #[test]
    fn tet10_corner_nodes_interpolate_kronecker_delta() {
        // Nₐ(node_b) = δₐᵦ — the defining property. Corner natural
        // coords + the six edge-midpoint natural coords.
        let t = ref_tet10();
        let mut nat = [[0.0; 3]; 10];
        let corners = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        nat[..4].copy_from_slice(&corners);
        for (m, &(p, q)) in TET10_EDGES.iter().enumerate() {
            for c in 0..3 {
                nat[4 + m][c] = 0.5 * (corners[p][c] + corners[q][c]);
            }
        }
        for (b, &node) in nat.iter().enumerate() {
            let sh = t.shape(node).unwrap();
            for (a, &na) in sh.n.iter().enumerate() {
                let want = if a == b { 1.0 } else { 0.0 };
                assert!((na - want).abs() < 1e-10, "N{a}(node{b}) = {na} ≠ {want}");
            }
        }
    }

    #[test]
    fn hex8_volume_is_exact() {
        let h = brick_hex8(2.0, 3.0, 4.0);
        let v = h.volume().unwrap();
        assert!((v - 24.0).abs() < 1e-9, "hex volume {v} ≠ 24");
    }

    #[test]
    fn tet4_and_tet10_volume_match_for_straight_edges() {
        // A straight-edged Tet10 spans the same region as its Tet4.
        let v4 = ref_tet4().volume().unwrap();
        let v10 = ref_tet10().volume().unwrap();
        assert!((v4 - 1.0 / 6.0).abs() < 1e-12, "tet4 vol {v4}");
        assert!((v10 - 1.0 / 6.0).abs() < 1e-12, "tet10 vol {v10}");
    }

    #[test]
    fn stiffness_matrices_are_symmetric() {
        let d = d_simple();
        for ke in [
            ref_tet4().stiffness(&d).unwrap(),
            brick_hex8(1.0, 1.0, 1.0).stiffness(&d).unwrap(),
            ref_tet10().stiffness(&d).unwrap(),
        ] {
            let n = ke.nrows();
            for i in 0..n {
                for j in 0..n {
                    let asym = (ke[(i, j)] - ke[(j, i)]).abs();
                    assert!(
                        asym < 1e-9 * ke[(i, i)].abs().max(1.0),
                        "Kₑ not symmetric at ({i},{j})"
                    );
                }
            }
        }
    }

    #[test]
    fn stiffness_has_six_zero_energy_rigid_modes() {
        // Each element's Kₑ must annihilate every rigid-body motion:
        // 3 translations + 3 (small) rotations → 6 zero-energy modes.
        let d = d_simple();
        let cases: Vec<(Box<dyn SolidElement>, usize)> = vec![
            (Box::new(ref_tet4()), 4),
            (Box::new(brick_hex8(1.5, 2.0, 1.0)), 8),
            (Box::new(ref_tet10()), 10),
        ];
        for (elem, n) in cases {
            let ke = elem.stiffness(&d).unwrap();
            let coords = elem.node_coords();
            // Translations.
            for axis in 0..3 {
                let mut t = DVector::zeros(3 * n);
                for a in 0..n {
                    t[3 * a + axis] = 1.0;
                }
                let f = &ke * &t;
                assert!(
                    f.norm() < 1e-7 * ke.norm(),
                    "{:?}: rigid translation {axis} produced force {}",
                    elem.element_type(),
                    f.norm()
                );
            }
            // Infinitesimal rotation about Z: u = (-y, x, 0).
            let mut r = DVector::zeros(3 * n);
            for a in 0..n {
                r[3 * a] = -coords[a].y;
                r[3 * a + 1] = coords[a].x;
            }
            let f = &ke * &r;
            assert!(
                f.norm() < 1e-7 * ke.norm() * r.amax().max(1.0),
                "{:?}: rigid rotation produced force {}",
                elem.element_type(),
                f.norm()
            );
        }
    }

    #[test]
    fn consistent_mass_distributes_total_mass() {
        // Sum of all entries of Mₑ = 3·ρ·V (the total mass counted once
        // per spatial direction).
        let density = 7850.0;
        let cases: Vec<Box<dyn SolidElement>> = vec![
            Box::new(ref_tet4()),
            Box::new(brick_hex8(2.0, 1.0, 3.0)),
            Box::new(ref_tet10()),
        ];
        for elem in cases {
            let me = elem.consistent_mass(density).unwrap();
            let vol = elem.volume().unwrap();
            let total: f64 = me.iter().sum();
            let want = 3.0 * density * vol;
            assert!(
                (total - want).abs() < 1e-6 * want,
                "{:?}: mass sum {total} ≠ 3ρV {want}",
                elem.element_type()
            );
        }
    }

    #[test]
    fn consistent_mass_is_positive_definite() {
        let cases: Vec<Box<dyn SolidElement>> = vec![
            Box::new(ref_tet4()),
            Box::new(brick_hex8(1.0, 2.0, 1.5)),
            Box::new(ref_tet10()),
        ];
        for elem in cases {
            let me = elem.consistent_mass(2700.0).unwrap();
            assert!(
                me.clone().cholesky().is_some(),
                "{:?}: Mₑ not positive-definite",
                elem.element_type()
            );
        }
    }

    #[test]
    fn degenerate_elements_are_rejected() {
        // A flat (zero-volume) tet.
        let flat = Tet4::new([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ]);
        let d = d_simple();
        assert!(flat.stiffness(&d).is_none());
        // A collapsed hex (all top nodes on the bottom face).
        let collapsed = Hex8::new([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ]);
        assert!(collapsed.stiffness(&d).is_none());
    }

    #[test]
    fn von_mises_basic_values() {
        // Pure hydrostatic → 0.
        assert!(von_mises(&[5.0, 5.0, 5.0, 0.0, 0.0, 0.0]) < 1e-9);
        // Uniaxial σ → σ.
        assert!((von_mises(&[3.0, 0.0, 0.0, 0.0, 0.0, 0.0]) - 3.0).abs() < 1e-12);
    }
}
