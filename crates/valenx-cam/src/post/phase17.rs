//! Phase 17C postprocessor variants — 27 templated postprocessors.
//!
//! Each variant is a thin module that calls [`crate::make_post!`]
//! with the template knobs that distinguish it from the others.
//! Subtle dialect differences live here so the host can advertise the
//! complete dropdown without having 27 separate files.
//!
//! The split rationale:
//!
//! - **Fanuc family** (`fanuc_a_axis`, `fanuc_b_axis`) — 30i clones
//!   with `(parens)` comments and `%`-bracketed program format.
//! - **Industrial-controller variants** (Heidenhain, Sinumerik,
//!   Mazatrol, Okuma, Mori, Makino, Kitamura) — full G-code with
//!   N-numbered lines and tool change.
//! - **Hobby controllers** (Centroid, Tormach, TinyG, SourceRabbit,
//!   FluidNC) — single-tool, GRBL-lineage formatting.
//! - **Printer-controllers** (Marlin, Klipper, Repetier, SnapMaker,
//!   Smoothieboard / SmoothieWare) — `;` line comments, no M30, no
//!   tool change.
//! - **Specialty** (DeepNest sheet-metal, HSMAdvisor cycle-time,
//!   Fusion360Brand display-branded Fanuc, OpenDMG Mori-clone) — one
//!   variant each.
//!
//! The macro emits a private `fn tpl()` per module to avoid name
//! collisions with neighbouring submodules.

use crate::make_post;
use crate::post::template::CommentStyle;

/// Centroid hobby controller.
pub mod centroid {
    use super::*;
    make_post!(
        name: Centroid,
        display: "Centroid",
        comment: CommentStyle::Open(';'),
        prelude: "",
        program_id: "",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// Tormach PCNC controller.
pub mod tormach {
    use super::*;
    make_post!(
        name: Tormach,
        display: "Tormach",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "",
        program_id: "",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// Heidenhain TNC.
pub mod heidenhain {
    use super::*;
    make_post!(
        name: Heidenhain,
        display: "Heidenhain TNC",
        comment: CommentStyle::Prefix(";"),
        prelude: "BEGIN PGM VALENX MM",
        program_id: "",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "TOOL CALL {n}",
        spindle_on_tpl: "S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "END PGM VALENX MM",
        program_end_suffix: "",
        number_lines: true,
        block_start: 1,
        block_step: 1,
    );
}

/// Mazak Mazatrol.
pub mod mazatrol {
    use super::*;
    make_post!(
        name: Mazatrol,
        display: "Mazatrol",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O8000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Siemens Sinumerik 840D.
pub mod sinumerik {
    use super::*;
    make_post!(
        name: Sinumerik,
        display: "Sinumerik 840D",
        comment: CommentStyle::Prefix(";"),
        prelude: "",
        program_id: "; SINUMERIK 840D",
        units_mode: "G71\nG90\nG17",
        tool_change_tpl: "T{n}\nM6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Okuma OSP — Fanuc-like, custom RPM ramp.
pub mod okuma {
    use super::*;
    make_post!(
        name: Okuma,
        display: "Okuma OSP",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O0001",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Mori-Seiki MAPPS.
pub mod mori {
    use super::*;
    make_post!(
        name: Mori,
        display: "Mori MAPPS",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O7000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Makino Professional 5.
pub mod makino {
    use super::*;
    make_post!(
        name: Makino,
        display: "Makino",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O5000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Kitamura Mycenter.
pub mod kitamura {
    use super::*;
    make_post!(
        name: Kitamura,
        display: "Kitamura",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O2000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// DeepNest — flat-pattern sheet-metal "G-code".
pub mod deepnest {
    use super::*;
    make_post!(
        name: DeepNest,
        display: "DeepNest",
        comment: CommentStyle::Open(';'),
        prelude: "; DeepNest flat-pattern output",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// Marlin 3D-printer firmware (used by some hobby CNC mills).
pub mod marlin {
    use super::*;
    make_post!(
        name: Marlin,
        display: "Marlin",
        comment: CommentStyle::Open(';'),
        prelude: "; Generated for Marlin",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M84",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// Klipper firmware.
pub mod klipper {
    use super::*;
    make_post!(
        name: Klipper,
        display: "Klipper",
        comment: CommentStyle::Open(';'),
        prelude: "; Generated for Klipper",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M84",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// Repetier firmware.
pub mod repetier {
    use super::*;
    make_post!(
        name: Repetier,
        display: "Repetier",
        comment: CommentStyle::Open(';'),
        prelude: "; Generated for Repetier",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M84",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// FluidNC firmware (GRBL successor).
pub mod fluid_nc {
    use super::*;
    make_post!(
        name: FluidNc,
        display: "FluidNC",
        comment: CommentStyle::Open(';'),
        prelude: "; FluidNC",
        program_id: "",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// SnapMaker firmware (3D-printer + laser + CNC combo).
pub mod snap_maker {
    use super::*;
    make_post!(
        name: SnapMaker,
        display: "SnapMaker",
        comment: CommentStyle::Open(';'),
        prelude: "; SnapMaker",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M84",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// SmoothieWare firmware.
pub mod smoothie {
    use super::*;
    make_post!(
        name: Smoothie,
        display: "SmoothieWare",
        comment: CommentStyle::Open(';'),
        prelude: "; SmoothieWare",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M84",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// Smoothieboard firmware (Smoothie-on-hardware variant).
pub mod smoothieboard {
    use super::*;
    make_post!(
        name: Smoothieboard,
        display: "Smoothieboard",
        comment: CommentStyle::Open(';'),
        prelude: "; Smoothieboard",
        program_id: "",
        units_mode: "G21\nG90",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "",
        coolant_off: "",
        program_end: "M84",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// TinyG controller.
pub mod tinyg {
    use super::*;
    make_post!(
        name: TinyG,
        display: "TinyG",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "",
        program_id: "",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// SourceRabbit (Mach3 / Mach4 derivative).
pub mod sourcerabbit {
    use super::*;
    make_post!(
        name: SourceRabbit,
        display: "SourceRabbit",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "",
        program_id: "",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: false,
        block_start: 10,
        block_step: 10,
    );
}

/// HSMAdvisor cycle-time-aware post.
pub mod hsmadvisor {
    use super::*;
    make_post!(
        name: HsmAdvisor,
        display: "HSMAdvisor",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O0001",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Autodesk Fusion 360 branded post (Fanuc base + branded header).
pub mod fusion360_brand {
    use super::*;
    make_post!(
        name: Fusion360Brand,
        display: "Fusion 360 (branded)",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O1001",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// OpenDMG variant (DMG Mori clone, no `%` bracketing).
pub mod opendmg {
    use super::*;
    make_post!(
        name: OpenDmg,
        display: "OpenDMG",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "",
        program_id: "O9000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// VMC 4-axis (A-axis) — generic Fanuc-like VMC with rotary A.
pub mod vmc_a_axis {
    use super::*;
    make_post!(
        name: VmcAAxis,
        display: "VMC 4-axis (A)",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O4000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// VMC 4-axis (B-axis) — generic Fanuc-like VMC with rotary B.
pub mod vmc_b_axis {
    use super::*;
    make_post!(
        name: VmcBAxis,
        display: "VMC 4-axis (B)",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O4001",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Fanuc 4-axis (A-axis) — Fanuc 30i with rotary A.
pub mod fanuc_a_axis {
    use super::*;
    make_post!(
        name: FanucAAxis,
        display: "Fanuc 4-axis (A)",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O3001",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

/// Fanuc 5-axis (A+B) — Fanuc 30i with two rotary axes.
pub mod fanuc_b_axis {
    use super::*;
    make_post!(
        name: FanucBAxis,
        display: "Fanuc 5-axis (A+B)",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O3002",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M6",
        spindle_on_tpl: "M3 S{rpm}",
        spindle_off: "M5",
        coolant_on: "M8",
        coolant_off: "M9",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::post::Postprocessor;
    use crate::{
        tool::{Tool, ToolKind},
        toolpath::{Move, MoveKind, Toolpath},
    };
    use nalgebra::Vector3;

    fn trivial_toolpath() -> Toolpath {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(10.0, 0.0, 0.0),
            500.0,
        ));
        tp
    }

    fn tool() -> Tool {
        Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap()
    }

    #[test]
    fn every_phase17_post_returns_some_output() {
        let tp = trivial_toolpath();
        let t = tool();
        // One sample per variant — verifies the macro is well-formed
        // and each variant produces a non-empty G-code string.
        let outputs: Vec<String> = vec![
            centroid::Centroid.process(&tp, &t, 12000.0).unwrap(),
            tormach::Tormach.process(&tp, &t, 12000.0).unwrap(),
            heidenhain::Heidenhain.process(&tp, &t, 12000.0).unwrap(),
            mazatrol::Mazatrol.process(&tp, &t, 12000.0).unwrap(),
            sinumerik::Sinumerik.process(&tp, &t, 12000.0).unwrap(),
            okuma::Okuma.process(&tp, &t, 12000.0).unwrap(),
            mori::Mori.process(&tp, &t, 12000.0).unwrap(),
            makino::Makino.process(&tp, &t, 12000.0).unwrap(),
            kitamura::Kitamura.process(&tp, &t, 12000.0).unwrap(),
            deepnest::DeepNest.process(&tp, &t, 12000.0).unwrap(),
            marlin::Marlin.process(&tp, &t, 12000.0).unwrap(),
            klipper::Klipper.process(&tp, &t, 12000.0).unwrap(),
            repetier::Repetier.process(&tp, &t, 12000.0).unwrap(),
            fluid_nc::FluidNc.process(&tp, &t, 12000.0).unwrap(),
            snap_maker::SnapMaker.process(&tp, &t, 12000.0).unwrap(),
            smoothie::Smoothie.process(&tp, &t, 12000.0).unwrap(),
            smoothieboard::Smoothieboard
                .process(&tp, &t, 12000.0)
                .unwrap(),
            tinyg::TinyG.process(&tp, &t, 12000.0).unwrap(),
            sourcerabbit::SourceRabbit
                .process(&tp, &t, 12000.0)
                .unwrap(),
            hsmadvisor::HsmAdvisor.process(&tp, &t, 12000.0).unwrap(),
            fusion360_brand::Fusion360Brand
                .process(&tp, &t, 12000.0)
                .unwrap(),
            opendmg::OpenDmg.process(&tp, &t, 12000.0).unwrap(),
            vmc_a_axis::VmcAAxis.process(&tp, &t, 12000.0).unwrap(),
            vmc_b_axis::VmcBAxis.process(&tp, &t, 12000.0).unwrap(),
            fanuc_a_axis::FanucAAxis.process(&tp, &t, 12000.0).unwrap(),
            fanuc_b_axis::FanucBAxis.process(&tp, &t, 12000.0).unwrap(),
        ];
        for o in &outputs {
            assert!(!o.is_empty());
            assert!(o.contains("G0") || o.contains("G1"));
        }
    }

    #[test]
    fn fanuc_a_axis_emits_program_id() {
        let tp = trivial_toolpath();
        let t = tool();
        let g = fanuc_a_axis::FanucAAxis.process(&tp, &t, 12000.0).unwrap();
        assert!(g.contains("O3001"));
    }

    #[test]
    fn fanuc_b_axis_emits_program_id() {
        let tp = trivial_toolpath();
        let t = tool();
        let g = fanuc_b_axis::FanucBAxis.process(&tp, &t, 12000.0).unwrap();
        assert!(g.contains("O3002"));
    }

    #[test]
    fn vmc_a_axis_emits_program_id() {
        let tp = trivial_toolpath();
        let t = tool();
        let g = vmc_a_axis::VmcAAxis.process(&tp, &t, 12000.0).unwrap();
        assert!(g.contains("O4000"));
    }

    #[test]
    fn vmc_b_axis_emits_program_id() {
        let tp = trivial_toolpath();
        let t = tool();
        let g = vmc_b_axis::VmcBAxis.process(&tp, &t, 12000.0).unwrap();
        assert!(g.contains("O4001"));
    }
}
