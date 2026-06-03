//! Integration tests for multi-op chaining + stock containment.
//!
//! Tasks 30-32 in the Phase 10 plan: ensure ops stay inside the
//! stock AABB and that face + pocket + drill chain together via
//! [`valenx_cam::Toolpath::concatenate`].

use nalgebra::Vector3;
use valenx_cam::op;
use valenx_cam::operation::{DrillParams, FaceParams, PocketParams};
use valenx_cam::stock::Stock;
use valenx_cam::tool::{Tool, ToolKind};
use valenx_cam::toolpath::Toolpath;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

fn cube_mesh(size: f64) -> Mesh {
    let s = size * 0.5;
    let nodes = vec![
        Vector3::new(-s, -s, -s),
        Vector3::new(s, -s, -s),
        Vector3::new(s, s, -s),
        Vector3::new(-s, s, -s),
        Vector3::new(-s, -s, s),
        Vector3::new(s, -s, s),
        Vector3::new(s, s, s),
        Vector3::new(-s, s, s),
    ];
    let conn: Vec<u32> = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5,
    ];
    let mut mesh = Mesh::new("cube");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh
}

#[test]
fn profile_xy_positions_stay_within_a_reasonable_envelope() {
    // We don't enforce strict stock AABB containment (the toolpath
    // legitimately rises above the stock for safe-Z rapids) but the
    // XY positions of the cut moves should stay close to the cube
    // boundary they were generated from.
    let mesh = cube_mesh(10.0);
    let stock = Stock::new(
        Vector3::new(-20.0, -20.0, -5.0),
        Vector3::new(40.0, 40.0, 10.0),
        "wood",
    )
    .unwrap();
    let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
    let params = valenx_cam::operation::ProfileParams {
        tool_id: 1,
        step_down: 2.0,
        depth: 2.0,
        ..Default::default()
    };
    let tp = op::profile::generate(&stock, &mesh, &params, &tool).unwrap();
    let (min, max) = tp.bounding_box().unwrap();
    // Profile offsets outward by tool radius (3) so XY range is
    // ~[-8, +8]. Z range is [stock_top - depth, stock_top + safe_z].
    assert!(
        min.x >= -10.0 && max.x <= 10.0,
        "XY range out of bounds: x in [{}, {}]",
        min.x,
        max.x
    );
    assert!(
        min.y >= -10.0 && max.y <= 10.0,
        "XY range out of bounds: y in [{}, {}]",
        min.y,
        max.y
    );
    // Z stays above stock_bottom and below safe_z.
    assert!(min.z >= -5.0, "Z too low: {}", min.z);
    assert!(max.z <= stock.top_z() + 10.0, "Z too high: {}", max.z);
}

#[test]
fn face_plus_pocket_plus_drill_chain_concatenates() {
    let mesh = cube_mesh(10.0);
    let stock = Stock::new(
        Vector3::new(-6.0, -6.0, -5.0),
        Vector3::new(12.0, 12.0, 10.0),
        "wood",
    )
    .unwrap();
    let face_tool = Tool::new(1, "FM10", ToolKind::FaceMill, 10.0, 30.0, 4, "carbide").unwrap();
    let em_tool = Tool::new(2, "EM2", ToolKind::EndMill, 2.0, 25.0, 2, "carbide").unwrap();
    let drill_tool = Tool::new(3, "Drill3", ToolKind::Drill, 3.0, 30.0, 2, "HSS").unwrap();

    let face = op::face::generate(
        &stock,
        &FaceParams {
            tool_id: 1,
            step_over: 4.0,
            step_down: 0.5,
            depth: 0.5,
            ..Default::default()
        },
        &face_tool,
    )
    .unwrap();
    let pocket = op::pocket::generate(
        &stock,
        &mesh,
        &PocketParams {
            tool_id: 2,
            step_over: 0.8,
            step_down: 1.0,
            depth: 2.0,
            ..Default::default()
        },
        &em_tool,
    )
    .unwrap();
    let drill = op::drill::generate(
        &stock,
        &DrillParams {
            tool_id: 3,
            peck_depth: 1.0,
            total_depth: 3.0,
            hole_positions: vec![Vector3::new(-3.0, -3.0, 0.0), Vector3::new(3.0, 3.0, 0.0)],
            ..Default::default()
        },
        &drill_tool,
    )
    .unwrap();

    let mut chain = Toolpath::new();
    chain.concatenate(&face);
    chain.concatenate(&pocket);
    chain.concatenate(&drill);

    assert_eq!(chain.len(), face.len() + pocket.len() + drill.len());
    assert!(!chain.is_empty());
}

#[test]
fn safe_z_helper_matches_stock_top_plus_clearance() {
    let s = Stock::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(10.0, 10.0, 5.0),
        "x",
    )
    .unwrap();
    let z = op::safe_z_for(&s, 2.0);
    assert!((z - 7.0).abs() < 1e-9);
}
