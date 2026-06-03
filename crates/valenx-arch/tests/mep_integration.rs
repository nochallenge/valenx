//! Integration tests for the MEP entities — verifies they flow
//! through document tessellation, schedule grouping, and the
//! summary/kind hooks.

use nalgebra::Vector3;
use valenx_arch::{
    ArchDocument, ArchEntity, ArchEntityKind, CableSegmentParams, ConduitSegmentParams,
    DuctSegmentParams, DuctShape, EquipmentKind, FlowDirection, MepEquipmentParams,
    PipeSegmentParams, Schedule,
};

#[test]
fn document_tessellates_through_mep_entities() {
    let mut d = ArchDocument::new("MEP");
    d.add_entity(ArchEntity::DuctSegment(DuctSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(3.0, 0.0, 0.0),
        shape: DuctShape::Rectangular {
            width: 0.4,
            height: 0.3,
        },
        material: "Galv steel".into(),
        flow_direction: FlowDirection::SourceToSink,
    }));
    d.add_entity(ArchEntity::PipeSegment(PipeSegmentParams {
        start: Vector3::new(0.0, 1.0, 0.5),
        end: Vector3::new(3.0, 1.0, 0.5),
        diameter: 0.05,
        material: "Copper".into(),
        fluid: "Cold water".into(),
        operating_pressure: 4.0e5,
    }));
    d.add_entity(ArchEntity::MepEquipment(MepEquipmentParams {
        position: Vector3::new(1.0, 2.0, 0.0),
        size: [0.8, 0.8, 2.0],
        kind: EquipmentKind::ElectricalPanel,
        tag: "P-1".into(),
        description: "Main panel".into(),
    }));

    let mesh = d.tessellate_all(0.1).unwrap();
    // Every MEP segment / equipment tessellates to 12 triangles each.
    assert_eq!(mesh.total_elements(), 36);
    assert!(d.bbox().is_some(), "document should have a bounding box");
}

#[test]
fn schedule_groups_mep_kinds_and_aggregates_lengths() {
    let mut d = ArchDocument::new("MEP");
    d.add_entity(ArchEntity::DuctSegment(DuctSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(3.0, 0.0, 0.0),
        shape: DuctShape::Round { diameter: 0.25 },
        material: "Galv".into(),
        flow_direction: FlowDirection::SourceToSink,
    }));
    d.add_entity(ArchEntity::DuctSegment(DuctSegmentParams {
        start: Vector3::new(3.0, 0.0, 0.0),
        end: Vector3::new(3.0, 4.0, 0.0),
        shape: DuctShape::Round { diameter: 0.25 },
        material: "Galv".into(),
        flow_direction: FlowDirection::SourceToSink,
    }));
    d.add_entity(ArchEntity::PipeSegment(PipeSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(5.0, 0.0, 0.0),
        diameter: 0.05,
        material: "Copper".into(),
        fluid: "Water".into(),
        operating_pressure: 4e5,
    }));
    d.add_entity(ArchEntity::CableSegment(CableSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(10.0, 0.0, 0.0),
        diameter: 0.012,
        conductor_csa_mm2: 2.5,
        voltage: 230.0,
        material: "PVC".into(),
    }));
    d.add_entity(ArchEntity::ConduitSegment(ConduitSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(7.0, 0.0, 0.0),
        outer_diameter: 0.025,
        inner_diameter: 0.020,
        material: "EMT".into(),
    }));
    d.add_entity(ArchEntity::MepEquipment(MepEquipmentParams {
        position: Vector3::zeros(),
        size: [1.0, 0.5, 0.5],
        kind: EquipmentKind::AirHandlingUnit,
        tag: "AHU".into(),
        description: "".into(),
    }));

    let s = Schedule::from_document(&d);

    let duct = s.entries.get(&ArchEntityKind::DuctSegment).unwrap();
    assert_eq!(duct.count, 2);
    assert!((duct.linear_m - 7.0).abs() < 1e-9);

    let pipe = s.entries.get(&ArchEntityKind::PipeSegment).unwrap();
    assert_eq!(pipe.count, 1);
    assert!((pipe.linear_m - 5.0).abs() < 1e-9);

    let cable = s.entries.get(&ArchEntityKind::CableSegment).unwrap();
    assert!((cable.linear_m - 10.0).abs() < 1e-9);

    let conduit = s.entries.get(&ArchEntityKind::ConduitSegment).unwrap();
    assert!((conduit.linear_m - 7.0).abs() < 1e-9);

    let eq = s.entries.get(&ArchEntityKind::MepEquipment).unwrap();
    assert_eq!(eq.count, 1);
    // Equipment volume = 1.0 × 0.5 × 0.5 = 0.25 m³.
    assert!((eq.volume_m3 - 0.25).abs() < 1e-9);

    let csv = s.to_csv();
    assert!(csv.contains("Duct,2,"));
    assert!(csv.contains("Pipe,1,"));
    assert!(csv.contains("Cable,1,"));
    assert!(csv.contains("Conduit,1,"));
    assert!(csv.contains("Equipment,1,"));
}

#[test]
fn entity_summary_includes_mep_descriptors() {
    let ent = ArchEntity::DuctSegment(DuctSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(2.0, 0.0, 0.0),
        shape: DuctShape::Round { diameter: 0.3 },
        material: "Galv steel".into(),
        flow_direction: FlowDirection::SourceToSink,
    });
    let s = ent.summary();
    assert!(s.starts_with("Duct"));
    assert!(s.contains("Round"));
    assert!(s.contains("Supply"));

    let ent = ArchEntity::PipeSegment(PipeSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(2.0, 0.0, 0.0),
        diameter: 0.05,
        material: "Copper".into(),
        fluid: "Hot water".into(),
        operating_pressure: 4e5,
    });
    let s = ent.summary();
    assert!(s.starts_with("Pipe"));
    assert!(s.contains("Hot water"));
    assert!(s.contains("Copper"));
}

#[test]
fn persist_round_trips_mep_entities() {
    use valenx_arch::ArchFile;
    let mut d = ArchDocument::new("Persist");
    d.add_entity(ArchEntity::DuctSegment(DuctSegmentParams {
        start: Vector3::zeros(),
        end: Vector3::new(2.0, 0.0, 0.0),
        shape: DuctShape::Oval {
            width: 0.4,
            height: 0.2,
        },
        material: "PVC".into(),
        flow_direction: FlowDirection::Return,
    }));
    d.add_entity(ArchEntity::MepEquipment(MepEquipmentParams {
        position: Vector3::new(1.0, 0.0, 0.0),
        size: [1.0, 1.0, 1.0],
        kind: EquipmentKind::Pump,
        tag: "P-101".into(),
        description: "Booster".into(),
    }));
    let ron = ArchFile::from_document(&d).to_ron().unwrap();
    let back = ArchFile::from_ron(&ron).unwrap();
    assert_eq!(back.document, d);
}
