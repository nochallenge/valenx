//! Element-validation suite for the native FEA element library
//! (Phase 24.8).
//!
//! ## What this is
//!
//! The element library ([`crate::elements`], [`crate::beam`]) ships
//! three new element types (Hex8, Tet10, the 3D beam). This module is
//! their **validation harness** — the standard finite-element
//! correctness tests, each asserting a genuine analytic value, never a
//! weakened tolerance:
//!
//! - **Patch test.** The fundamental FE correctness test. A small mesh
//!   has a *constant-strain* displacement field prescribed on its
//!   boundary; a correct element must reproduce that field exactly at
//!   every interior node and recover the *constant* analytic stress
//!   everywhere. An element that fails the patch test does not
//!   converge. [`patch_test_hex8`], [`patch_test_tet10`],
//!   [`patch_test_tet4`] run it for each continuum element.
//! - **Beam-bending convergence.** A cantilever is solved with the
//!   Tet4, Hex8 and Tet10 elements at increasing mesh density; the tip
//!   deflection is compared to the Euler-Bernoulli analytic value. The
//!   test demonstrates that the constant-strain Tet4 is badly
//!   over-stiff while the Hex8 and Tet10 converge fast.
//!   [`beam_convergence_study`] returns the numbers.
//! - **Beam-element benchmarks.** The 3D beam element is checked
//!   against analytic cantilever / simply-supported deflections and the
//!   first natural frequency — those live as `#[test]`s in
//!   [`crate::beam`].
//!
//! The functions are public so a caller (a `tests/` harness, the
//! desktop QA panel) can run the suite and inspect the numbers; the
//! `#[test]`s at the bottom assert the pass criteria.

use nalgebra::Vector3;

use valenx_mesh::Mesh;

use crate::elements::von_mises;
use crate::material::FemMaterial;
use crate::native_solver::{solve_linear_static_mixed, NodalConstraint, NodalForce};

/// A constant (homogeneous) strain state, as the 3×3 displacement-
/// gradient matrix `A` of the linear field `u(x) = A·x`.
///
/// The strain is the symmetric part of `A`; an antisymmetric part is a
/// rigid rotation. The patch test uses a non-trivial `A` with both
/// normal and shear strain.
#[derive(Copy, Clone, Debug)]
pub struct ConstantStrainField {
    /// The 3×3 displacement-gradient matrix `A` (row-major).
    pub gradient: [[f64; 3]; 3],
}

impl ConstantStrainField {
    /// The displacement `u = A·x` at a point.
    pub fn displacement(&self, x: Vector3<f64>) -> [f64; 3] {
        let a = &self.gradient;
        [
            a[0][0] * x.x + a[0][1] * x.y + a[0][2] * x.z,
            a[1][0] * x.x + a[1][1] * x.y + a[1][2] * x.z,
            a[2][0] * x.x + a[2][1] * x.y + a[2][2] * x.z,
        ]
    }

    /// The Voigt strain `[εxx εyy εzz γxy γyz γzx]` of the field
    /// (engineering shear strain `γ = 2ε`).
    pub fn voigt_strain(&self) -> [f64; 6] {
        let a = &self.gradient;
        [
            a[0][0],
            a[1][1],
            a[2][2],
            a[0][1] + a[1][0],
            a[1][2] + a[2][1],
            a[2][0] + a[0][2],
        ]
    }
}

/// The result of one element patch test.
#[derive(Clone, Debug)]
pub struct PatchTestResult {
    /// Element family under test (a human label).
    pub element: &'static str,
    /// Largest displacement error at any node, relative to the field
    /// amplitude. A passing element drives this to solver precision.
    pub max_displacement_error: f64,
    /// Largest stress error at any node, relative to the analytic
    /// constant stress magnitude.
    pub max_stress_error: f64,
    /// Number of nodes that were *not* prescribed (the genuine
    /// interior the solver had to recover). A patch test with zero
    /// interior nodes only checks representability; with interior
    /// nodes it also checks the assembled equilibrium.
    pub interior_nodes: usize,
}

impl PatchTestResult {
    /// Did the element pass — both errors below the given tolerance?
    pub fn passed(&self, tol: f64) -> bool {
        self.max_displacement_error < tol && self.max_stress_error < tol
    }
}

/// Run the constant-strain **patch test** on an arbitrary solid mesh.
///
/// `mesh` must carry a supported continuum-element block. The boundary
/// nodes (those on the outer faces of the mesh's axis-aligned bounding
/// box) have the field `u = A·x` prescribed; every *interior* node is
/// left free. A correct element reproduces the linear field exactly at
/// the interior nodes and recovers the constant analytic stress
/// `σ = D·ε` everywhere.
///
/// Returns the worst-case displacement and stress errors — both must
/// be at solver precision for the element to pass.
pub fn run_patch_test(
    element_label: &'static str,
    mesh: &Mesh,
    material: &FemMaterial,
    field: &ConstantStrainField,
) -> PatchTestResult {
    // Bounding box.
    let mut lo = Vector3::new(f64::MAX, f64::MAX, f64::MAX);
    let mut hi = Vector3::new(f64::MIN, f64::MIN, f64::MIN);
    for p in &mesh.nodes {
        lo = lo.inf(p);
        hi = hi.sup(p);
    }
    let tol_geom = 1.0e-9 * (hi - lo).norm().max(1.0);

    // Prescribe every boundary node; leave interior nodes free.
    let mut constraints = Vec::new();
    let mut interior = 0usize;
    for (n, p) in mesh.nodes.iter().enumerate() {
        let on_boundary = (p.x - lo.x).abs() < tol_geom
            || (p.x - hi.x).abs() < tol_geom
            || (p.y - lo.y).abs() < tol_geom
            || (p.y - hi.y).abs() < tol_geom
            || (p.z - lo.z).abs() < tol_geom
            || (p.z - hi.z).abs() < tol_geom;
        if on_boundary {
            let u = field.displacement(*p);
            constraints.push(NodalConstraint::displaced(n, u));
        } else {
            interior += 1;
        }
    }

    let no_force: [NodalForce; 0] = [];
    let sol = solve_linear_static_mixed(mesh, material, &constraints, &no_force)
        .expect("patch-test solve failed");

    // Field amplitude for relative errors.
    let amp = mesh
        .nodes
        .iter()
        .map(|p| {
            let u = field.displacement(*p);
            (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt()
        })
        .fold(0.0_f64, f64::max)
        .max(1.0e-30);

    // Displacement error at every node.
    let mut max_disp_err = 0.0_f64;
    for (n, p) in mesh.nodes.iter().enumerate() {
        let want = field.displacement(*p);
        let got = sol.displacement[n];
        for k in 0..3 {
            max_disp_err = max_disp_err.max((got[k] - want[k]).abs() / amp);
        }
    }

    // Stress error: the analytic constant stress is σ = D·ε.
    let d = crate::native_solver::elasticity_matrix(material).expect("bad material");
    let eps = field.voigt_strain();
    let eps_v = nalgebra::Vector6::from_row_slice(&eps);
    let sigma_analytic = d * eps_v;
    let sigma_mag = sigma_analytic.norm().max(1.0e-30);
    let mut max_stress_err = 0.0_f64;
    for s in &sol.stress {
        for k in 0..6 {
            max_stress_err = max_stress_err.max((s[k] - sigma_analytic[k]).abs() / sigma_mag);
        }
    }

    PatchTestResult {
        element: element_label,
        max_displacement_error: max_disp_err,
        max_stress_error: max_stress_err,
        interior_nodes: interior,
    }
}

/// A non-trivial constant-strain field for the patch test — normal
/// strain on all three axes plus all three shear components, so the
/// test exercises every entry of the strain-displacement matrix.
pub fn patch_test_field() -> ConstantStrainField {
    // u = A·x. Small amplitudes keep it a small-strain state.
    ConstantStrainField {
        gradient: [
            [1.0e-3, 0.4e-3, 0.2e-3],
            [0.3e-3, 0.8e-3, 0.5e-3],
            [0.1e-3, 0.6e-3, 1.2e-3],
        ],
    }
}

/// Run the patch test for the **Hex8** element on a 2×2×2 brick mesh
/// (which has one genuine interior node).
pub fn patch_test_hex8() -> PatchTestResult {
    // Round-10: structured_hex_mesh is now fallible. The hard-coded
    // 2×2×2/1.0/1.3/0.7 args are within bounds — expect() is the
    // honest "this can't fail" assertion.
    let mesh = crate::meshgen::structured_hex_mesh(1.0, 1.3, 0.7, 2, 2, 2)
        .expect("patch_test_hex8: 2x2x2 box dimensions are within allocation bounds");
    run_patch_test("Hex8", &mesh, &patch_material(), &patch_test_field())
}

/// Run the patch test for the **Tet10** element on a 2×2×2 box mesh.
pub fn patch_test_tet10() -> PatchTestResult {
    // Round-10: structured_tet10_mesh is now fallible.
    let mesh = crate::meshgen::structured_tet10_mesh(1.0, 1.3, 0.7, 2, 2, 2)
        .expect("patch_test_tet10: 2x2x2 box dimensions are within allocation bounds");
    run_patch_test("Tet10", &mesh, &patch_material(), &patch_test_field())
}

/// Run the patch test for the **Tet4** element on a 2×2×2 box mesh —
/// the original element, re-checked through the generic path.
pub fn patch_test_tet4() -> PatchTestResult {
    // Round-5: structured_box_mesh is now fallible. Hard-coded
    // 2×2×2/1.0/1.3/0.7 params are within bounds — expect() is the
    // honest "this can't fail with these specific args" assertion.
    let mesh = crate::native_solver::structured_box_mesh(1.0, 1.3, 0.7, 2, 2, 2)
        .expect("patch_test_tet4: 2x2x2 box dimensions are within allocation bounds");
    run_patch_test("Tet4", &mesh, &patch_material(), &patch_test_field())
}

/// The material the patch tests use — a generic steel with a non-zero
/// Poisson ratio (so the lateral coupling is exercised too).
fn patch_material() -> FemMaterial {
    FemMaterial {
        youngs_modulus: 200.0e9,
        poisson_ratio: 0.3,
        ..FemMaterial::default()
    }
}

/// One row of a beam-bending convergence study: an element family at a
/// given mesh density and its predicted tip deflection.
#[derive(Clone, Debug)]
pub struct ConvergencePoint {
    /// Element family label.
    pub element: &'static str,
    /// Number of elements along the beam's length.
    pub elements_along_length: usize,
    /// Total degrees of freedom of the model.
    pub dof: usize,
    /// Predicted tip deflection magnitude (metres).
    pub tip_deflection: f64,
    /// Tip deflection as a fraction of the Euler-Bernoulli analytic
    /// value — `1.0` is exact, `< 1` is over-stiff.
    pub fraction_of_analytic: f64,
}

/// Result of the cantilever beam-bending convergence study.
#[derive(Clone, Debug)]
pub struct BeamConvergenceStudy {
    /// The Euler-Bernoulli analytic tip deflection `δ = P·L³/(3·E·I)`.
    pub analytic_deflection: f64,
    /// One [`ConvergencePoint`] per (element family, mesh density),
    /// appended coarse-mesh-first within each family.
    pub points: Vec<ConvergencePoint>,
}

impl BeamConvergenceStudy {
    /// Every convergence point of one element family, in the order they
    /// were computed (coarsest mesh first).
    pub fn family(&self, element: &str) -> Vec<&ConvergencePoint> {
        self.points
            .iter()
            .filter(|p| p.element == element)
            .collect()
    }

    /// The finest-mesh result for an element family, as a fraction of
    /// the analytic deflection.
    pub fn finest_fraction(&self, element: &str) -> Option<f64> {
        self.family(element)
            .last()
            .map(|p| p.fraction_of_analytic)
    }
}

/// Run the **cantilever beam-bending convergence study**.
///
/// A slender cantilever (length `L`, square `b×b` section) is clamped
/// at one end and given a transverse tip load. It is solved with the
/// Tet4, Hex8 and Tet10 elements at increasing mesh densities; each
/// tip deflection is compared to the Euler-Bernoulli analytic value
/// `δ = P·L³/(3·E·I)`.
///
/// The study demonstrates the headline result: the constant-strain
/// Tet4 is badly over-stiff in bending (a coarse Tet4 mesh predicts a
/// small fraction of the true deflection), while the Hex8 and the
/// quadratic Tet10 converge to the analytic value quickly.
pub fn beam_convergence_study() -> BeamConvergenceStudy {
    let (lx, ly, lz) = (10.0_f64, 1.0_f64, 1.0_f64);
    let material = FemMaterial {
        youngs_modulus: 200.0e9,
        poisson_ratio: 0.3,
        ..FemMaterial::default()
    };
    let load = 1.0e4;

    // Euler-Bernoulli analytic tip deflection.
    let i_section = ly * lz.powi(3) / 12.0;
    let analytic = load * lx.powi(3) / (3.0 * material.youngs_modulus * i_section);

    let mut points = Vec::new();
    // (length subdivisions, transverse subdivisions) — coarse → fine.
    let densities = [(4usize, 1usize), (8, 1), (16, 2), (24, 2)];

    for family in ["Tet4", "Hex8", "Tet10"] {
        for &(nx, nt) in &densities {
            let mesh = match family {
                // Round-5: structured_box_mesh is now fallible.
                // Iteration densities top out at 24×2×2 = 96 cells —
                // well within the cap.
                "Tet4" => crate::native_solver::structured_box_mesh(lx, ly, lz, nx, nt, nt)
                    .expect("convergence study densities are within allocation bounds"),
                // Round-10: both fallible. Iteration densities top
                // out at the same scale as Tet4 — within the cap.
                "Hex8" => crate::meshgen::structured_hex_mesh(lx, ly, lz, nx, nt, nt)
                    .expect("convergence study densities are within allocation bounds"),
                "Tet10" => crate::meshgen::structured_tet10_mesh(lx, ly, lz, nx, nt, nt)
                    .expect("convergence study densities are within allocation bounds"),
                _ => unreachable!(),
            };
            let tip = cantilever_tip(&mesh, &material, nx, nt, nt, load);
            points.push(ConvergencePoint {
                element: match family {
                    "Tet4" => "Tet4",
                    "Hex8" => "Hex8",
                    _ => "Tet10",
                },
                elements_along_length: nx,
                dof: 3 * mesh.nodes.len(),
                tip_deflection: tip,
                fraction_of_analytic: tip / analytic,
            });
        }
    }

    BeamConvergenceStudy {
        analytic_deflection: analytic,
        points,
    }
}

/// Solve one cantilever case and return the mean tip deflection
/// magnitude.
///
/// The mesh's corner-node grid is laid out identically by all three
/// structured generators, so the same `(i,j,k) → node` arithmetic
/// picks the clamped face (`x=0`) and the loaded face (`x=L`).
fn cantilever_tip(
    mesh: &Mesh,
    material: &FemMaterial,
    nx: usize,
    ny: usize,
    nz: usize,
    load: f64,
) -> f64 {
    // Corner-node id of grid point (i,j,k).
    let nid = |i: usize, j: usize, k: usize| i + (nx + 1) * j + (nx + 1) * (ny + 1) * k;

    // Clamp every corner node on the x=0 face.
    let mut constraints = Vec::new();
    for k in 0..=nz {
        for j in 0..=ny {
            constraints.push(NodalConstraint::fixed(nid(0, j, k)));
        }
    }
    // Distribute a downward (-Z) load over the corner nodes of the x=L
    // face.
    let tip_nodes: Vec<usize> = (0..=nz)
        .flat_map(|k| (0..=ny).map(move |j| nid(nx, j, k)))
        .collect();
    let per = -load / tip_nodes.len() as f64;
    let forces: Vec<NodalForce> = tip_nodes
        .iter()
        .map(|&n| NodalForce {
            node: n,
            force: [0.0, 0.0, per],
        })
        .collect();

    let sol = solve_linear_static_mixed(mesh, material, &constraints, &forces)
        .expect("cantilever solve failed");
    tip_nodes
        .iter()
        .map(|&n| sol.displacement[n][2].abs())
        .sum::<f64>()
        / tip_nodes.len() as f64
}

/// von Mises stress of an analytic uniaxial state — a tiny convenience
/// the validation `#[test]`s use.
pub fn uniaxial_von_mises(sigma: f64) -> f64 {
    von_mises(&[sigma, 0.0, 0.0, 0.0, 0.0, 0.0])
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- PATCH TESTS -------------------------------------------------

    #[test]
    fn hex8_passes_the_constant_strain_patch_test() {
        // The fundamental FE correctness test for the Hex8 element: a
        // constant-strain field prescribed on the boundary must be
        // reproduced EXACTLY at the interior node and the recovered
        // stress must be the constant analytic σ = D·ε everywhere.
        let r = patch_test_hex8();
        assert!(r.interior_nodes >= 1, "patch test had no interior node");
        assert!(
            r.passed(1e-8),
            "Hex8 patch test FAILED: disp err {:.3e}, stress err {:.3e}",
            r.max_displacement_error,
            r.max_stress_error
        );
    }

    #[test]
    fn tet10_passes_the_constant_strain_patch_test() {
        let r = patch_test_tet10();
        assert!(r.interior_nodes >= 1, "patch test had no interior node");
        assert!(
            r.passed(1e-8),
            "Tet10 patch test FAILED: disp err {:.3e}, stress err {:.3e}",
            r.max_displacement_error,
            r.max_stress_error
        );
    }

    #[test]
    fn tet4_passes_the_constant_strain_patch_test() {
        // The original element, re-verified through the generic
        // mixed-element path.
        let r = patch_test_tet4();
        assert!(
            r.passed(1e-8),
            "Tet4 patch test FAILED: disp err {:.3e}, stress err {:.3e}",
            r.max_displacement_error,
            r.max_stress_error
        );
    }

    // --- BEAM-BENDING CONVERGENCE ------------------------------------

    /// The beam-bending convergence study is moderately expensive and
    /// three tests inspect it — compute it once and share.
    fn shared_study() -> &'static BeamConvergenceStudy {
        use std::sync::OnceLock;
        static STUDY: OnceLock<BeamConvergenceStudy> = OnceLock::new();
        STUDY.get_or_init(beam_convergence_study)
    }

    #[test]
    fn hex8_and_tet10_converge_much_faster_than_tet4_in_bending() {
        // The headline element-library result. A cantilever solved
        // with each element family, refined; Tet4 is badly over-stiff
        // (constant-strain locking), Hex8 and Tet10 converge to the
        // Euler-Bernoulli analytic deflection.
        let study = shared_study();
        assert!(study.analytic_deflection > 0.0);

        let tet4_best = study.finest_fraction("Tet4").unwrap();
        let hex8_best = study.finest_fraction("Hex8").unwrap();
        let tet10_best = study.finest_fraction("Tet10").unwrap();

        // Tet4 is measurably over-stiff: even the finest Tet4 mesh
        // here recovers well under the true deflection.
        assert!(
            tet4_best < 0.85,
            "Tet4 should be over-stiff in bending, got fraction {tet4_best:.3}"
        );
        // Hex8 converges close to the analytic value.
        assert!(
            hex8_best > 0.9 && hex8_best < 1.15,
            "Hex8 should converge to beam theory, got fraction {hex8_best:.3}"
        );
        // Tet10 converges close to the analytic value.
        assert!(
            tet10_best > 0.9 && tet10_best < 1.15,
            "Tet10 should converge to beam theory, got fraction {tet10_best:.3}"
        );
        // The new elements are dramatically better than Tet4.
        assert!(
            hex8_best > tet4_best + 0.2,
            "Hex8 ({hex8_best:.3}) should beat Tet4 ({tet4_best:.3}) by a wide margin"
        );
        assert!(
            tet10_best > tet4_best + 0.2,
            "Tet10 ({tet10_best:.3}) should beat Tet4 ({tet4_best:.3}) by a wide margin"
        );
    }

    #[test]
    fn each_element_monotonically_converges_under_refinement() {
        // Refining a mesh must move the tip deflection *toward* the
        // analytic value — for every element family. A displacement-
        // based element under a load converges monotonically from
        // below, so the fraction sequence is non-decreasing.
        let study = shared_study();
        for family in ["Tet4", "Hex8", "Tet10"] {
            let fr: Vec<f64> = study
                .family(family)
                .iter()
                .map(|p| p.fraction_of_analytic)
                .collect();
            assert!(fr.len() >= 2, "{family}: not enough convergence points");
            for w in fr.windows(2) {
                assert!(
                    w[1] >= w[0] - 1e-6,
                    "{family}: refinement decreased the deflection {} → {}",
                    w[0],
                    w[1]
                );
            }
            // The finest mesh is at least as close to analytic as the
            // coarsest (monotone approach from below).
            let coarse = *fr.first().unwrap();
            let fine = *fr.last().unwrap();
            assert!(
                fine >= coarse - 1e-6,
                "{family}: finest mesh {fine} not closer than coarsest {coarse}"
            );
        }
    }

    // --- BENCHMARKS --------------------------------------------------

    #[test]
    fn hex8_uniaxial_tension_recovers_hookes_law() {
        // Standard benchmark: a Hex8 bar in uniaxial tension must
        // reproduce σ = E·ε to solver precision (a constant-strain
        // state — the patch test in disguise, displacement driven).
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        // Round-10: structured_hex_mesh is now fallible.
        let mesh = crate::meshgen::structured_hex_mesh(lx, ly, lz, nx, ny, nz)
            .expect("uniaxial-tension hex mesh dims are within bounds");
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.0, // pure uniaxial, no lateral coupling
            ..FemMaterial::default()
        };
        let nid = |i: usize, j: usize, k: usize| i + (nx + 1) * j + (nx + 1) * (ny + 1) * k;
        let strain = 1.0e-3;
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                // x=0 face fully clamped — kills every rigid-body mode.
                // With ν = 0 there is no Poisson contraction to be
                // restrained, so the σxx field stays exactly uniform.
                constraints.push(NodalConstraint::fixed(nid(0, j, k)));
                // x=L face displaced by ε·L in X (transverse DOFs free
                // so the bar can extend without spurious stress).
                constraints.push(NodalConstraint {
                    node: nid(nx, j, k),
                    fixed: [Some(strain * lx), None, None],
                });
            }
        }
        let no_force: [NodalForce; 0] = [];
        let sol = solve_linear_static_mixed(&mesh, &mat, &constraints, &no_force).unwrap();
        // σxx = E·ε everywhere.
        let want = mat.youngs_modulus * strain;
        for s in &sol.stress {
            let rel = (s[0] - want).abs() / want;
            assert!(rel < 1e-6, "Hex8 σxx {} vs E·ε {want} (rel {rel})", s[0]);
        }
    }

    #[test]
    fn tet10_uniaxial_tension_recovers_hookes_law() {
        // The same uniaxial benchmark for the quadratic Tet10. Every
        // node (corners AND mid-edge) is prescribed from the analytic
        // linear field u = ε·x — for a quadratic element the mid-edge
        // nodes must be on the field too, and the structured generator
        // put them at the true edge midpoints so u = ε·x is exact.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (2, 1, 1);
        // Round-10: structured_tet10_mesh is now fallible.
        let mesh = crate::meshgen::structured_tet10_mesh(lx, ly, lz, nx, ny, nz)
            .expect("uniaxial-tension tet10 mesh dims are within bounds");
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.0,
            ..FemMaterial::default()
        };
        let strain = 1.0e-3;
        let mut constraints = Vec::new();
        for (n, p) in mesh.nodes.iter().enumerate() {
            constraints.push(NodalConstraint::displaced(n, [strain * p.x, 0.0, 0.0]));
        }
        let no_force: [NodalForce; 0] = [];
        let sol = solve_linear_static_mixed(&mesh, &mat, &constraints, &no_force).unwrap();
        let want = mat.youngs_modulus * strain;
        for s in &sol.stress {
            let rel = (s[0] - want).abs() / want;
            assert!(rel < 1e-6, "Tet10 σxx {} vs E·ε {want} (rel {rel})", s[0]);
        }
    }

    #[test]
    fn convergence_study_dof_counts_are_sane() {
        // A coarse Tet4 cantilever has fewer DOF than the finest
        // Tet10 — the study really does refine. And the Tet10 mesh
        // carries more DOF than the Tet4 mesh at the same cell count
        // (the mid-edge nodes).
        let study = shared_study();
        let tet4 = study.family("Tet4");
        let tet10 = study.family("Tet10");
        assert!(tet4.first().unwrap().dof < tet4.last().unwrap().dof);
        // Same cell density, Tet10 has the mid-edge nodes → more DOF.
        assert!(
            tet10[0].dof > tet4[0].dof,
            "Tet10 should carry more DOF than Tet4 at equal cell count"
        );
    }

    #[test]
    fn uniaxial_von_mises_equals_the_stress() {
        assert!((uniaxial_von_mises(5.0e6) - 5.0e6).abs() < 1.0);
    }
}


