//! Integration tests for the IFC4 coverage expansion.
//!
//! Verifies the IFC4 writer emits well-formed STEP-21 lines for every
//! new entity kind (covering, curtain wall, footing, pile, railing,
//! ramp, chimney, furnishing, opening + rel-voids, space boundary,
//! duct, pipe, cable, conduit, equipment) and that Psets attach via
//! IfcRelDefinesByProperties.

use nalgebra::Vector3;
use valenx_arch::ifc::{
    emit_pset, emit_rel_space_boundary, emit_rel_voids_element, ifc_guid_v4, write_cable,
    write_chimney, write_conduit, write_covering, write_curtain_wall, write_document, write_duct,
    write_footing, write_furnishing, write_mep_equipment, write_opening_for_door,
    write_opening_for_window, write_pile, write_pipe, write_railing, write_ramp, IfcWriter,
    PropValue,
};
use valenx_arch::{
    ArchDocument, ArchEntity, BeamParams, BeamSection, CableSegmentParams, ColumnParams,
    ColumnSection, ConduitSegmentParams, DoorParams, DoorStyle, DuctSegmentParams, DuctShape,
    EquipmentKind, FlowDirection, MepEquipmentParams, PipeSegmentParams, Side, SlabParams,
    SpaceParams, WallParams, WindowParams, WindowStyle,
};

fn fixture_wall() -> WallParams {
    WallParams {
        start: Vector3::zeros(),
        end: Vector3::new(5.0, 0.0, 0.0),
        height: 2.7,
        thickness: 0.2,
        material: "Brick".into(),
    }
}

fn fixture_slab() -> SlabParams {
    SlabParams {
        boundary: vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(4.0, 3.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
        ],
        thickness: 0.2,
        material: "Concrete".into(),
        structural: None,
    }
}

fn fixture_column() -> ColumnParams {
    ColumnParams {
        base: Vector3::zeros(),
        height: 3.0,
        cross_section: ColumnSection::Rectangle {
            width: 0.3,
            depth: 0.3,
        },
        material: "Steel".into(),
        structural: None,
    }
}

#[test]
fn ifc_covering_emits_predefined_type() {
    let mut w = IfcWriter::new();
    let id = write_covering(&mut w, &fixture_slab(), ".FLOORING.");
    let s = w.finish();
    assert!(s.contains("IFCCOVERING"));
    assert!(s.contains(".FLOORING."));
    assert!(id > 0);
}

#[test]
fn ifc_curtain_wall_emits_correct_class() {
    let mut w = IfcWriter::new();
    write_curtain_wall(&mut w, &fixture_wall());
    let s = w.finish();
    assert!(s.contains("IFCCURTAINWALL"));
    assert!(s.contains("Curtain Wall"));
}

#[test]
fn ifc_footing_emits_strip_footing_predef() {
    let mut w = IfcWriter::new();
    write_footing(&mut w, &fixture_slab());
    let s = w.finish();
    assert!(s.contains("IFCFOOTING"));
    assert!(s.contains(".STRIP_FOOTING."));
}

#[test]
fn ifc_pile_emits_bored_predef() {
    let mut w = IfcWriter::new();
    write_pile(&mut w, &fixture_column());
    let s = w.finish();
    assert!(s.contains("IFCPILE"));
    assert!(s.contains(".BORED."));
}

#[test]
fn ifc_railing_emits_handrail_predef() {
    let mut w = IfcWriter::new();
    write_railing(
        &mut w,
        Vector3::zeros(),
        Vector3::new(3.0, 0.0, 0.0),
        1.1,
        "Steel",
    );
    let s = w.finish();
    assert!(s.contains("IFCRAILING"));
    assert!(s.contains(".HANDRAIL."));
}

#[test]
fn ifc_ramp_emits_straight_run_predef() {
    let mut w = IfcWriter::new();
    write_ramp(
        &mut w,
        &valenx_arch::StairParams {
            base: Vector3::zeros(),
            direction: Vector3::new(1.0, 0.0, 0.0),
            total_rise: 1.0,
            total_run: 8.0,
            num_steps: 1,
            width: 1.5,
        },
    );
    let s = w.finish();
    assert!(s.contains("IFCRAMP"));
    assert!(s.contains(".STRAIGHT_RUN_RAMP."));
}

#[test]
fn ifc_chimney_emits_with_material_in_name() {
    let mut w = IfcWriter::new();
    write_chimney(&mut w, &fixture_column());
    let s = w.finish();
    assert!(s.contains("IFCCHIMNEY"));
    assert!(s.contains("Steel"));
}

#[test]
fn ifc_furnishing_emits_with_tag() {
    let mut w = IfcWriter::new();
    let e = MepEquipmentParams {
        position: Vector3::zeros(),
        size: [0.8, 0.4, 1.6],
        kind: EquipmentKind::LightFitting,
        tag: "F-101".into(),
        description: "Cabinet".into(),
    };
    write_furnishing(&mut w, &e);
    let s = w.finish();
    assert!(s.contains("IFCFURNISHINGELEMENT"));
    assert!(s.contains("F-101"));
}

#[test]
fn ifc_duct_with_oval_shape_picks_oval_label() {
    let mut w = IfcWriter::new();
    let d = DuctSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(2.0, 0.0, 0.0),
        shape: DuctShape::Oval {
            width: 0.4,
            height: 0.2,
        },
        material: "Galv steel".into(),
        flow_direction: FlowDirection::Bidirectional,
    };
    write_duct(&mut w, &d);
    let s = w.finish();
    assert!(s.contains("IFCDUCTSEGMENT"));
    assert!(s.contains("Oval"));
    assert!(s.contains("Bidir"));
}

#[test]
fn ifc_pipe_includes_fluid_and_material_in_name() {
    let mut w = IfcWriter::new();
    let p = PipeSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(1.5, 0.0, 0.0),
        diameter: 0.05,
        material: "Copper".into(),
        fluid: "Cold water".into(),
        operating_pressure: 4.0e5,
    };
    write_pipe(&mut w, &p);
    let s = w.finish();
    assert!(s.contains("IFCPIPESEGMENT"));
    assert!(s.contains("Copper"));
    assert!(s.contains("Cold water"));
}

#[test]
fn ifc_cable_includes_voltage_in_name() {
    let mut w = IfcWriter::new();
    let c = CableSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(8.0, 0.0, 0.0),
        diameter: 0.012,
        conductor_csa_mm2: 4.0,
        voltage: 230.0,
        material: "PVC".into(),
    };
    write_cable(&mut w, &c);
    let s = w.finish();
    assert!(s.contains("IFCCABLESEGMENT"));
    assert!(s.contains("230V"));
    assert!(s.contains("4.0mm2"));
}

#[test]
fn ifc_conduit_emits_cable_carrier_segment_with_predefined() {
    let mut w = IfcWriter::new();
    let c = ConduitSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(3.0, 0.0, 0.0),
        outer_diameter: 0.025,
        inner_diameter: 0.020,
        material: "EMT".into(),
    };
    write_conduit(&mut w, &c);
    let s = w.finish();
    assert!(s.contains("IFCCABLECARRIERSEGMENT"));
    assert!(s.contains(".CONDUITSEGMENT."));
}

#[test]
fn ifc_mep_equipment_picks_kind_specific_entity() {
    let kinds = [
        (EquipmentKind::AirHandlingUnit, "IFCAIRTERMINALBOX"),
        (EquipmentKind::Pump, "IFCPUMP"),
        (EquipmentKind::Valve, "IFCVALVE"),
        (EquipmentKind::SprinklerHead, "IFCFIRESUPPRESSIONTERMINAL"),
        (
            EquipmentKind::ElectricalPanel,
            "IFCELECTRICDISTRIBUTIONBOARD",
        ),
        (EquipmentKind::LightFitting, "IFCLIGHTFIXTURE"),
    ];
    for (kind, expected) in kinds {
        let mut w = IfcWriter::new();
        let e = MepEquipmentParams {
            position: Vector3::zeros(),
            size: [0.5, 0.5, 0.5],
            kind,
            tag: "X-1".into(),
            description: "".into(),
        };
        write_mep_equipment(&mut w, &e);
        let s = w.finish();
        assert!(
            s.contains(expected),
            "kind {kind:?} missing {expected} in:\n{s}"
        );
    }
}

#[test]
fn pset_attaches_via_rel_defines_by_properties() {
    let mut w = IfcWriter::new();
    // Stub an entity so the Pset has something to attach to.
    let stub = w.add("IFCWALL('guid',$,'W',$,$,$,$,$,$)");
    emit_pset(
        &mut w,
        stub,
        "Pset_WallCommon",
        &[
            ("LoadBearing", PropValue::Bool(true)),
            ("ThermalTransmittance", PropValue::Real(0.21)),
            ("FireRating", PropValue::Label("F60".into())),
        ],
    );
    let s = w.finish();
    assert!(s.contains("IFCPROPERTYSET"));
    assert!(s.contains("Pset_WallCommon"));
    assert!(s.contains("IFCRELDEFINESBYPROPERTIES"));
    assert!(s.contains("IFCPROPERTYSINGLEVALUE"));
    assert!(s.contains("LoadBearing"));
    assert!(s.contains("ThermalTransmittance"));
    assert!(s.contains("FireRating"));
    // IfcReal + IfcBoolean + IfcLabel measure types.
    assert!(s.contains("IFCBOOLEAN"));
    assert!(s.contains("IFCREAL"));
    assert!(s.contains("IFCLABEL"));
}

#[test]
fn rel_voids_element_links_host_and_opening() {
    let mut w = IfcWriter::new();
    let wall = fixture_wall();
    let win = WindowParams {
        host: 1,
        position_along_wall: 2.5,
        position_height: 1.0,
        width: 1.0,
        height: 1.0,
        frame_thickness: 0.05,
        style: WindowStyle::Casement,
    };
    let opening = write_opening_for_window(&mut w, &win, &wall);
    // Fabricate a host id.
    let host = w.add("IFCWALL('h',$,'HW',$,$,$,$,$,$)");
    let rel = emit_rel_voids_element(&mut w, host, opening);
    let s = w.finish();
    assert!(s.contains("IFCOPENINGELEMENT"));
    assert!(s.contains("IFCRELVOIDSELEMENT"));
    assert!(rel > opening);
}

#[test]
fn rel_space_boundary_emits_with_internal_kind() {
    let mut w = IfcWriter::new();
    let sp = w.add("IFCSPACE('g',$,'S',$,$,$,$,$,.ELEMENT.,.INTERNAL.,$)");
    let el = w.add("IFCWALL('g',$,'W',$,$,$,$,$,$)");
    emit_rel_space_boundary(&mut w, sp, el, ".INTERNAL.");
    let s = w.finish();
    assert!(s.contains("IFCRELSPACEBOUNDARY"));
    assert!(s.contains(".INTERNAL."));
    assert!(s.contains(".PHYSICAL."));
}

#[test]
fn write_door_for_opening_emits_door_opening_label() {
    let mut w = IfcWriter::new();
    let wall = fixture_wall();
    let door = DoorParams {
        host: 1,
        position_along_wall: 2.5,
        width: 0.9,
        height: 2.1,
        style: DoorStyle::Single,
        hinge_side: Side::Left,
    };
    write_opening_for_door(&mut w, &door, &wall);
    let s = w.finish();
    assert!(s.contains("IFCOPENINGELEMENT"));
    assert!(s.contains("DoorOpening"));
}

#[test]
fn ifc_guid_is_22_chars() {
    let g = ifc_guid_v4();
    assert_eq!(g.len(), 22);
}

#[test]
fn full_document_writes_psets_for_each_kind() {
    let mut d = ArchDocument::new("integration-test");
    let wall_id = d.add_entity(ArchEntity::Wall(fixture_wall()));
    d.add_entity(ArchEntity::Slab(fixture_slab()));
    d.add_entity(ArchEntity::Column(fixture_column()));
    d.add_entity(ArchEntity::Beam(BeamParams {
        start: Vector3::zeros(),
        end: Vector3::new(4.0, 0.0, 0.0),
        cross_section: BeamSection::Rectangle {
            width: 0.2,
            depth: 0.3,
        },
        orientation_angle: 0.0,
        material: "Steel".into(),
        structural: None,
    }));
    d.add_entity(ArchEntity::Window(WindowParams {
        host: wall_id,
        position_along_wall: 2.0,
        position_height: 1.0,
        width: 1.0,
        height: 1.0,
        frame_thickness: 0.05,
        style: WindowStyle::Casement,
    }));
    d.add_entity(ArchEntity::Door(DoorParams {
        host: wall_id,
        position_along_wall: 4.0,
        width: 0.9,
        height: 2.1,
        style: DoorStyle::Single,
        hinge_side: Side::Left,
    }));
    d.add_entity(ArchEntity::Space(SpaceParams {
        boundary: vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(4.0, 3.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
        ],
        ceiling_height: 2.7,
        space_name: "Living".into(),
    }));
    d.add_entity(ArchEntity::DuctSegment(DuctSegmentParams {
        start: Vector3::new(0.0, 0.0, 2.4),
        end: Vector3::new(4.0, 0.0, 2.4),
        shape: DuctShape::Round { diameter: 0.25 },
        material: "Galv steel".into(),
        flow_direction: FlowDirection::SourceToSink,
    }));
    d.add_entity(ArchEntity::PipeSegment(PipeSegmentParams {
        start: Vector3::new(0.0, 1.0, 0.5),
        end: Vector3::new(4.0, 1.0, 0.5),
        diameter: 0.05,
        material: "Copper".into(),
        fluid: "Cold water".into(),
        operating_pressure: 4.0e5,
    }));
    d.add_entity(ArchEntity::CableSegment(CableSegmentParams {
        start: Vector3::new(0.0, 0.0, 2.5),
        end: Vector3::new(4.0, 0.0, 2.5),
        diameter: 0.012,
        conductor_csa_mm2: 4.0,
        voltage: 230.0,
        material: "PVC".into(),
    }));
    d.add_entity(ArchEntity::ConduitSegment(ConduitSegmentParams {
        start: Vector3::new(0.0, 2.0, 2.5),
        end: Vector3::new(4.0, 2.0, 2.5),
        outer_diameter: 0.025,
        inner_diameter: 0.020,
        material: "EMT".into(),
    }));
    d.add_entity(ArchEntity::MepEquipment(MepEquipmentParams {
        position: Vector3::new(0.5, 1.0, 0.0),
        size: [1.5, 0.8, 1.2],
        kind: EquipmentKind::AirHandlingUnit,
        tag: "AHU-101".into(),
        description: "Roof-top air handler".into(),
    }));
    let tmp = std::env::temp_dir().join("valenx_arch_full_ifc4.ifc");
    write_document(&d, &tmp).unwrap();
    let s = std::fs::read_to_string(&tmp).unwrap();
    let _ = std::fs::remove_file(&tmp);

    // Per-entity Psets must show up.
    assert!(s.contains("Pset_WallCommon"));
    assert!(s.contains("Pset_SlabCommon"));
    assert!(s.contains("Pset_ColumnCommon"));
    assert!(s.contains("Pset_BeamCommon"));
    assert!(s.contains("Pset_WindowCommon"));
    assert!(s.contains("Pset_DoorCommon"));
    assert!(s.contains("Pset_SpaceCommon"));
    assert!(s.contains("Pset_DuctSegmentTypeCommon"));
    assert!(s.contains("Pset_PipeSegmentTypeCommon"));
    assert!(s.contains("Pset_CableSegmentTypeCommon"));
    assert!(s.contains("Pset_CableCarrierSegmentTypeCommon"));
    assert!(s.contains("Pset_DistributionElementCommon"));

    // Window + door must produce OpeningElement + RelVoidsElement +
    // RelFillsElement.
    assert!(s.contains("IFCOPENINGELEMENT"));
    assert!(s.contains("IFCRELVOIDSELEMENT"));
    assert!(s.contains("IFCRELFILLSELEMENT"));

    // Space boundary linking the wall and the space (the wall
    // midpoint at (2.5, 0, 1.35) is on the rectangle edge — the
    // boundary AABB is inclusive, so the boundary fires).
    assert!(s.contains("IFCRELSPACEBOUNDARY"));

    // The MEP class names + furnishing IFC4 entity prefix.
    assert!(s.contains("IFCDUCTSEGMENT"));
    assert!(s.contains("IFCPIPESEGMENT"));
    assert!(s.contains("IFCCABLESEGMENT"));
    assert!(s.contains("IFCCABLECARRIERSEGMENT"));
    assert!(s.contains("IFCAIRTERMINALBOX"));
}
