//! # valenx-fem
//!
//! Native FEM workbench for Valenx (Phase 24).
//!
//! Provides the in-process pre/post pipeline that sits *on top of*
//! the existing subprocess adapters in
//! `crates/valenx-adapters/fea/valenx-adapter-{calculix,elmer,code-aster}`.
//! The adapters do the actual solving; this crate owns the workbench
//! model:
//!
//! - **Pre.** [`FemAnalysis`] captures the source solid, the
//!   material, the constraints, the loads, the mesh parameters, and
//!   the desired solver. [`FemAnalysis::write_input_file`] emits a
//!   solver-specific input deck.
//! - **Solver dispatch.** [`FemAnalysis::launch_adapter_id`] is a thin
//!   description of which adapter ID to invoke; the actual subprocess
//!   call is left to the desktop shell so this crate stays free of
//!   tokio / process / network deps. v1 documents the dispatch
//!   contract via [`FemSolverChoice::adapter_id`].
//! - **Post.** [`FemResult`] parses solver output (displacement +
//!   stress field), [`FemResult::tessellate_visualization`] returns a
//!   coloured mesh for the viewport.
//!
//! ## Native solvers
//!
//! Eight genuine in-process FEA solvers ship in this crate — no
//! subprocess, no input deck. The original eight all assemble on the
//! constant-strain linear-tetrahedron element and
//! [`native_solver::structured_box_mesh`], a tetrahedral box mesh so
//! every solver is usable end-to-end without an external mesher.
//!
//! ## Element library (Phase 24.8)
//!
//! Commercial FEA carries dozens of element types because the linear
//! tetrahedron — the original sole element — is notoriously **over-
//! stiff in bending**. The [`elements`] module is the generic element
//! layer: a [`elements::SolidElement`] trait every 3-DOF-per-node
//! continuum element implements, with three implementations —
//! [`elements::Tet4`] (the linear tet), [`elements::Hex8`] (the
//! 8-node trilinear brick, the workhorse of commercial solid FEA) and
//! [`elements::Tet10`] (the 10-node *quadratic* tet, far more accurate
//! per DOF). [`assembly`] assembles a **mixed-element** mesh — Tet4,
//! Hex8 and Tet10 in any combination — into one coupled global system;
//! [`native_solver::solve_linear_static_mixed`] and
//! [`modal_solver::solve_modal_mixed`] are the mixed-element static and
//! modal solvers. The [`beam`] module is a separate native solver for
//! **2-node 3D Timoshenko beam** elements (6 DOF per node — axial,
//! biaxial bending, torsion) — the element frames and lattices are
//! built from. [`meshgen`] supplies structured Hex8 / Tet10 box meshes;
//! [`validation`] is the element-validation suite (constant-strain
//! patch tests + beam-bending convergence).
//!
//! - **Linear static** ([`native_solver`]). Assembles the global
//!   stiffness matrix, applies displacement / force boundary
//!   conditions, solves the sparse SPD system with a Cholesky
//!   factorisation, and recovers the nodal displacement field plus von
//!   Mises stress. See [`native_solver::solve_linear_static`].
//! - **Modal** ([`modal_solver`]). Pairs the elastic stiffness with a
//!   consistent mass matrix and solves the generalised symmetric
//!   eigenproblem `K φ = λ M φ` for the lowest `n` natural frequencies
//!   + mode shapes. See [`modal_solver::solve_modal`].
//! - **Steady-state thermal** ([`thermal_solver`]). Solves the steady
//!   heat-conduction equation `−∇·(k∇T) = q` with Dirichlet
//!   (fixed-temperature) and Neumann (heat-flux) boundary conditions,
//!   recovering the nodal temperature field and heat flux. See
//!   [`thermal_solver::solve_steady_thermal`].
//! - **Geometrically-nonlinear static** ([`nonlinear_solver`]). Wraps
//!   the linear Tet4 element in a **corotational** kinematic
//!   description + a **Newton-Raphson** iteration, so the solver stays
//!   correct under **large displacements and rotations** — a
//!   cantilever under a large tip load comes out stiffer than the
//!   linear prediction (geometric stiffening). See
//!   [`nonlinear_solver::solve_nonlinear_static`].
//! - **Material plasticity** ([`plasticity`]). **von Mises (J2)
//!   plasticity** with linear isotropic hardening — a
//!   **radial-return** (closest-point) stress update, the **consistent
//!   elastoplastic tangent**, and an incremental Newton-Raphson load
//!   loop. A metal modelled with it is elastic up to yield then flows
//!   plastically, and carries a permanent set when unloaded. See
//!   [`plasticity::solve_plastic`].
//! - **Contact** ([`contact`]). **Penalty-method node-to-surface
//!   contact** — gap detection against a rigid plane, a penalty
//!   stiffness opposing penetration, and the contact force assembled
//!   into the Newton residual. Two bodies pushed together do not
//!   interpenetrate. See [`contact::solve_contact`].
//! - **Transient structural dynamics** ([`dynamics`]). **Newmark-β**
//!   implicit time integration on the consistent mass + stiffness, so
//!   the solver can march a structure's transient vibration. A
//!   single-DOF spring-mass oscillator reproduces its analytic natural
//!   period. See [`dynamics::solve_transient_dynamics`].
//! - **Linear buckling** ([`buckling`]). **Eigenvalue buckling** — a
//!   **geometric (stress) stiffness** matrix assembled from a reference
//!   linear-static stress state, and the generalised eigenproblem
//!   `(K + λ·K_g)·φ = 0` solved for the critical load factors. A
//!   slender column's lowest buckling load approaches the Euler
//!   critical load. See [`buckling::solve_buckling`].
//!
//! Each module documents its honest scope (Tet4, isotropic, linear
//! isotropic hardening, penalty rigid-plane contact; none is a
//! CalculiX replacement — that gap stays documented in Tier 3).
//!
//! ## v1 limitations
//!
//! - The solver-specific input-deck emitters are minimal:
//!   `write_input_file` produces a CalculiX `.inp` plus a parallel
//!   Elmer / Code_Aster placeholder that the adapter's existing
//!   writer can refine. The 90% case (linear static, isotropic
//!   elastic, one material, surface BCs) is covered.
//! - Result parsing reads displacement / stress out of the canonical
//!   `valenx_fields::Results` returned by the adapter's
//!   `collect()` step — this lets us reuse the existing parser path.
//! - Animation playback (deformation interpolation) is a pure-math
//!   helper here; the egui slider lives in `valenx-app`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod analysis;
pub mod assembly;
pub mod beam;
pub mod buckling;
pub mod constraints;
pub mod contact;
pub mod dynamics;
pub mod elements;
pub mod loads;
pub mod material;
pub mod mesh_params;
pub mod meshgen;
pub mod modal_solver;
pub mod native_solver;
pub mod nonlinear_solver;
pub mod ordering;
pub mod plasticity;
pub mod result;
pub mod solver;
pub mod thermal_solver;
pub mod validation;

pub use analysis::{FemAnalysis, FemAnalysisError};
pub use assembly::{
    assemble_global_stiffness, assemble_global_stiffness_mass_mixed, has_solid_elements,
    recover_nodal_stress_mixed,
};
pub use beam::{
    axial_force_capacity, axial_rigidity, axial_strain_energy, axial_stress, beam_angle_of_twist,
    beam_axial_extension,
    beam_curvature, beam_transverse_shear_stress,
    bending_moment_capacity,
    bending_stress, bulk_modulus,
    cantilever_point_load_root_moment,
    cantilever_point_load_strain_energy, cantilever_tip_deflection, cantilever_tip_slope,
    cantilever_udl_root_moment, cantilever_udl_strain_energy, cantilever_udl_tip_deflection,
    cantilever_udl_tip_slope, circular_plastic_section_modulus,
    circular_polar_second_moment_of_area,
    circular_second_moment_of_area, elastic_section_modulus,
    euler_bernoulli_beam_frequency, fixed_fixed_center_deflection,
    fixed_fixed_point_load_end_moment, fixed_fixed_udl_center_deflection,
    fixed_fixed_udl_end_moment, flexural_rigidity, hollow_circular_polar_second_moment_of_area,
    hollow_circular_second_moment_of_area,
    lames_first_parameter,
    p_wave_modulus,
    polar_section_modulus,
    propped_cantilever_udl_fixed_end_moment, propped_cantilever_udl_prop_reaction,
    rectangular_plastic_section_modulus, rectangular_polar_second_moment_of_area,
    rectangular_second_moment_of_area, shear_modulus_from_youngs,
    simply_supported_center_deflection, simply_supported_end_slope,
    simply_supported_point_load_max_moment, simply_supported_point_load_strain_energy,
    simply_supported_udl_center_deflection,
    simply_supported_udl_end_slope, simply_supported_udl_max_moment,
    simply_supported_udl_strain_energy, solve_beam_modal,
    solve_beam_static, torsional_moment_capacity, torsional_rigidity, torsional_shear_stress,
    torsional_strain_energy, BeamConstraint,
    BeamElement, BeamLoad,
    BeamModalSolution, BeamMode,
    BeamSection, BeamSolution, BeamSolverError,
};
pub use buckling::{
    critical_buckling_stress, euler_critical_load, section_radius_of_gyration,
    slenderness_ratio, solve_buckling, BucklingMode, BucklingSolution, BucklingSolverError,
};
pub use constraints::FemConstraint;
pub use contact::{
    solve_contact, ContactControls, ContactPlane, ContactSolution,
};
pub use dynamics::{
    solve_transient_dynamics, DynamicsControls, DynamicsSolution, NewmarkParameters,
    NodalInitialState,
};
pub use elements::{Hex8, SolidElement, Tet10, Tet4};
pub use loads::FemLoad;
pub use meshgen::{structured_hex_mesh, structured_tet10_mesh};
pub use material::{material_library, FemMaterial, PlasticProperties};
pub use mesh_params::{ElementOrder, FemMeshParams};
pub use modal_solver::{
    solve_modal, solve_modal_mixed, ModalSolution, ModalSolverError, VibrationMode,
};
pub use native_solver::{
    check_dense_dofs, solve_linear_static, solve_linear_static_mixed, structured_box_mesh,
    NativeSolution, NativeSolverError, NodalConstraint, NodalForce, MAX_DENSE_DOFS,
};
pub use nonlinear_solver::{
    corotational_element, solve_nonlinear_static, NonlinearControls, NonlinearSolution,
};
pub use plasticity::{
    consistent_tangent, radial_return, solve_plastic, PlasticControls, PlasticSolution,
    PlasticState, ReturnResult,
};
pub use result::{FemResult, FemResultError};
pub use solver::FemSolverChoice;
pub use thermal_solver::{
    solve_steady_thermal, FixedTemperature, HeatLoad, ThermalSolution, ThermalSolverError,
};
