//! Postprocessor tests — known-pattern toolpath → expected G-code
//! string match.

use nalgebra::Vector3;
use valenx_cam::post::{
    fanuc::Fanuc, grbl::Grbl, linuxcnc::LinuxCnc, save_nc, PostKind, Postprocessor,
};
use valenx_cam::tool::{Tool, ToolKind};
use valenx_cam::toolpath::{Move, MoveKind, Toolpath};

fn make_toolpath() -> Toolpath {
    let mut tp = Toolpath::new();
    tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(10.0, 0.0, 5.0),
        0.0,
    ));
    tp.push(Move::new(
        MoveKind::Plunge,
        Vector3::new(10.0, 0.0, 0.0),
        200.0,
    ));
    tp.push(Move::new(
        MoveKind::Cut,
        Vector3::new(10.0, 10.0, 0.0),
        500.0,
    ));
    tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
    tp
}

fn make_tool() -> Tool {
    Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap()
}

#[test]
fn grbl_emits_expected_g_code() {
    let g = Grbl
        .process(&make_toolpath(), &make_tool(), 12000.0)
        .unwrap();
    assert!(g.contains("G21"));
    assert!(g.contains("G90"));
    assert!(g.contains("M3 S12000"));
    assert!(g.contains("G0 X10.000 Y0.000 Z5.000"));
    assert!(g.contains("G1 X10.000 Y0.000 Z0.000 F200"));
    assert!(g.contains("G1 X10.000 Y10.000 Z0.000 F500"));
    assert!(g.contains("M5"));
    assert!(g.contains("M30"));
    assert!(!g.contains("T1 M6"), "GRBL should not emit tool change");
}

#[test]
fn linuxcnc_emits_tool_change() {
    let g = LinuxCnc
        .process(&make_toolpath(), &make_tool(), 12000.0)
        .unwrap();
    assert!(g.contains("T1 M6"));
    assert!(g.contains("G0 X10.000 Y0.000 Z5.000"));
    assert!(g.contains("G1 X10.000 Y10.000 Z0.000 F500"));
}

#[test]
fn fanuc_emits_numbered_lines_and_program_brackets() {
    let g = Fanuc
        .process(&make_toolpath(), &make_tool(), 12000.0)
        .unwrap();
    assert!(g.starts_with("%"));
    assert!(g.contains("O1000"));
    assert!(g.contains("N10 G0 X0.000 Y0.000 Z5.000"));
    assert!(g.contains("N20 G0 X10.000"));
    assert!(g.contains("N30 G1 X10.000 Y0.000 Z0.000 F200"));
    assert!(g.trim_end().ends_with("%"));
}

#[test]
fn empty_toolpath_errors() {
    let err = Grbl
        .process(&Toolpath::new(), &make_tool(), 12000.0)
        .unwrap_err();
    assert_eq!(err.code(), "cam.empty_toolpath");
}

#[test]
fn cut_with_zero_feed_errors() {
    let mut tp = Toolpath::new();
    tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
    tp.push(Move::new(MoveKind::Cut, Vector3::new(1.0, 1.0, 0.0), 0.0));
    let err = Grbl.process(&tp, &make_tool(), 12000.0).unwrap_err();
    assert_eq!(err.code(), "cam.postprocessor_failed");
}

#[test]
fn save_nc_writes_file_round_trip() {
    let tp = make_toolpath();
    let tool = make_tool();
    let path = std::env::temp_dir().join("valenx_cam_save_nc_test.nc");
    save_nc(PostKind::Grbl, &tp, &tool, 12000.0, &path).unwrap();
    let read = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = read.lines().take(10).collect();
    assert!(lines[0].contains("valenx-cam GRBL"));
    assert!(read.contains("G1 X10.000 Y10.000 Z0.000 F500"));
    let _ = std::fs::remove_file(&path);
}

#[test]
fn post_kind_label() {
    assert_eq!(PostKind::Grbl.label(), "GRBL");
    assert_eq!(PostKind::LinuxCnc.label(), "LinuxCNC");
    assert_eq!(PostKind::Fanuc.label(), "Fanuc");
}
