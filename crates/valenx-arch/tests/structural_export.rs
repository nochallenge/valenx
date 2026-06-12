//! Integration test for the structural-model export.
//!
//! Builds a 2-column-1-beam portal frame, exports it to a
//! [`StructuralModel`], translates that model into the FEM beam
//! solver's input vectors, and verifies the solve runs end-to-end.

use nalgebra::Vector3;
use valenx_arch::{
    export_structural_model, ArchDocument, ArchEntity, BeamParams, BeamSection, ColumnParams,
    ColumnSection, StructuralMaterial, StructuralMember, StructuralModelOptions, SupportKind,
};
use valenx_fem::{
    material::FemMaterial, solve_beam_static, BeamConstraint, BeamElement, BeamLoad,
    BeamSection as FemBeamSection,
};

fn portal_frame_doc() -> ArchDocument {
    let mut d = ArchDocument::new("Portal");
    d.add_entity(ArchEntity::Column(ColumnParams {
        base: Vector3::new(0.0, 0.0, 0.0),
        height: 3.0,
        cross_section: ColumnSection::Rectangle {
            width: 0.3,
            depth: 0.3,
        },
        material: "Steel".into(),
        structural: Some(StructuralMember {
            material: StructuralMaterial::SteelS355,
            support: SupportKind::Clamped,
            applied_force: [0.0; 3],
            applied_moment: [0.0; 3],
            self_weight_load: false,
        }),
    }));
    d.add_entity(ArchEntity::Column(ColumnParams {
        base: Vector3::new(5.0, 0.0, 0.0),
        height: 3.0,
        cross_section: ColumnSection::Rectangle {
            width: 0.3,
            depth: 0.3,
        },
        material: "Steel".into(),
        structural: Some(StructuralMember {
            material: StructuralMaterial::SteelS355,
            support: SupportKind::Clamped,
            applied_force: [0.0; 3],
            applied_moment: [0.0; 3],
            self_weight_load: false,
        }),
    }));
    d.add_entity(ArchEntity::Beam(BeamParams {
        start: Vector3::new(0.0, 0.0, 3.0),
        end: Vector3::new(5.0, 0.0, 3.0),
        cross_section: BeamSection::IBeam {
            width: 0.2,
            depth: 0.4,
            flange_thickness: 0.02,
            web_thickness: 0.01,
        },
        orientation_angle: 0.0,
        material: "Steel".into(),
        structural: Some(StructuralMember {
            material: StructuralMaterial::SteelS355,
            support: SupportKind::Free,
            applied_force: [0.0, 0.0, -10_000.0],
            applied_moment: [0.0; 3],
            self_weight_load: false,
        }),
    }));
    d
}

#[test]
fn portal_frame_solves_through_fem_beam_solver() {
    let doc = portal_frame_doc();
    let model = export_structural_model(&doc, &StructuralModelOptions::default()).unwrap();

    assert_eq!(model.elements.len(), 3, "expected 3 beam elements");
    assert_eq!(model.nodes.len(), 4, "expected 4 unique nodes");
    assert_eq!(model.supports.len(), 2, "expected 2 clamped supports");
    assert_eq!(model.loads.len(), 1, "expected 1 nodal load");

    // Translate the StructuralModel to the FEM solver's vectors.
    let nodes: Vec<Vector3<f64>> = model.nodes.iter().map(|n| n.position).collect();
    let mat = &model.materials[0];
    let fem_mat = FemMaterial {
        name: mat.label().into(),
        youngs_modulus: mat.youngs_modulus(),
        poisson_ratio: mat.poisson_ratio(),
        density: mat.density(),
        thermal_conductivity: 50.0,
        plasticity: None,
    };
    let elements: Vec<BeamElement> = model
        .elements
        .iter()
        .map(|e| {
            BeamElement::new(
                e.start_node,
                e.end_node,
                FemBeamSection {
                    area: e.section.area,
                    iy: e.section.iy,
                    iz: e.section.iz,
                    j: e.section.j,
                    shear_y: 5.0 / 6.0,
                    shear_z: 5.0 / 6.0,
                },
            )
        })
        .collect();
    let constraints: Vec<BeamConstraint> = model
        .supports
        .iter()
        .map(|s| BeamConstraint {
            node: s.node,
            fixed: s.fixed,
        })
        .collect();
    let loads: Vec<BeamLoad> = model
        .loads
        .iter()
        .map(|l| BeamLoad {
            node: l.node,
            force: l.force,
            moment: l.moment,
        })
        .collect();

    let sol = solve_beam_static(&nodes, &elements, &fem_mat, &constraints, &loads).unwrap();

    assert_eq!(sol.translation.len(), nodes.len());
    assert_eq!(sol.rotation.len(), nodes.len());

    // Find the loaded crown node and check it deflects downward.
    let load_node = model.loads[0].node;
    let max_z = sol.translation[load_node][2];
    assert!(
        max_z < 0.0,
        "loaded crown should deflect downward, got {max_z}"
    );

    // Clamped bases should not move.
    for sup in &model.supports {
        let t = sol.translation[sup.node];
        let mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
        assert!(mag < 1.0e-6, "clamped node {} moved by {mag}", sup.node);
    }

    // The solve must produce a finite, non-trivial max translation.
    let max = sol.max_translation();
    assert!(max.is_finite() && max > 1.0e-12, "got max disp {max}");
}

#[test]
fn portal_frame_dof_counts_match_expected() {
    let doc = portal_frame_doc();
    let model = export_structural_model(&doc, &StructuralModelOptions::default()).unwrap();
    // 4 nodes × 6 DOF = 24 total DOF.
    assert_eq!(model.dof_count(), 24);
    // 2 clamped supports × 6 DOF = 12 constrained DOF.
    assert_eq!(model.constrained_dof_count(), 12);
}
