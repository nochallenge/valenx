//! Haas postprocessor — Fanuc-like with `M06 T{n}` order swap.

use nalgebra::Vector3;

use crate::{
    error::CamError,
    post::{
        template::{process_template, CommentStyle, PostTemplate},
        Postprocessor,
    },
    tool::Tool,
    toolpath::Toolpath,
};

/// Haas CNC postprocessor.
#[derive(Clone, Copy, Debug, Default)]
pub struct Haas;

fn tpl() -> PostTemplate {
    PostTemplate {
        name: "Haas",
        comment: CommentStyle::Pair('(', ')'),
        prelude: "%",
        program_id: "O1000",
        units_mode: "G21\nG90\nG17",
        tool_change_tpl: "T{n} M06",
        spindle_on_tpl: "M03 S{rpm}",
        spindle_off: "M05",
        coolant_on: "M08",
        coolant_off: "M09",
        program_end: "M30",
        program_end_suffix: "%\n",
        number_lines: true,
        block_start: 10,
        block_step: 10,
    }
}

impl Postprocessor for Haas {
    fn header(&self, tool: &Tool, spindle_rpm: f64) -> String {
        tpl().header(tool, spindle_rpm)
    }
    fn footer(&self) -> String {
        tpl().footer()
    }
    fn move_g0(&self, p: Vector3<f64>) -> String {
        crate::post::format_g0(p)
    }
    fn move_g1(&self, p: Vector3<f64>, feed: f64) -> Result<String, CamError> {
        crate::post::format_g1(p, feed)
    }
    fn spindle_on(&self, rpm: f64) -> String {
        tpl().spindle_on(rpm)
    }
    fn spindle_off(&self) -> String {
        tpl().spindle_off.into()
    }
    fn tool_change(&self, tool_id: u32) -> String {
        tpl().tool_change(tool_id)
    }
    fn coolant_on(&self) -> String {
        tpl().coolant_on.into()
    }
    fn coolant_off(&self) -> String {
        tpl().coolant_off.into()
    }
    fn process(&self, tp: &Toolpath, tool: &Tool, spindle_rpm: f64) -> Result<String, CamError> {
        process_template(&tpl(), tp, tool, spindle_rpm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        tool::ToolKind,
        toolpath::{Move, MoveKind},
    };

    #[test]
    fn haas_emits_m06_tool_change() {
        let t = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        let g = Haas.process(&tp, &t, 12000.0).unwrap();
        assert!(g.contains("T1 M06"));
        assert!(g.contains("M03 S12000"));
    }
}
