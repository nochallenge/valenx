//! GRBL postprocessor.
//!
//! GRBL is single-tool (no `T{n} M6`); the header sets mm-units,
//! absolute mode, XY-plane, and starts the spindle.
//!
//! Move format: `G0 X{x:.3} Y{y:.3} Z{z:.3}` / `G1 X{...} F{feed}`.
//! Three-decimal precision is GRBL-friendly.

use nalgebra::Vector3;

use crate::{
    error::CamError,
    post::{format_g0, format_g1, Postprocessor},
    tool::Tool,
    toolpath::{MoveKind, Toolpath},
};

/// GRBL postprocessor — emits the most-portable subset of G-code
/// that GRBL accepts.
#[derive(Clone, Copy, Debug, Default)]
pub struct Grbl;

impl Postprocessor for Grbl {
    fn header(&self, tool: &Tool, spindle_rpm: f64) -> String {
        // GRBL is single-tool, so tool number goes in a comment only.
        format!(
            "; valenx-cam GRBL\n; Tool: T{} {} D{:.2}mm\nG21\nG90\nG17\n{}\n",
            tool.id,
            tool.name,
            tool.diameter_mm,
            self.spindle_on(spindle_rpm).trim()
        )
    }

    fn footer(&self) -> String {
        format!("{}\nM30\n", self.spindle_off().trim())
    }

    fn move_g0(&self, p: Vector3<f64>) -> String {
        format_g0(p)
    }

    fn move_g1(&self, p: Vector3<f64>, feed: f64) -> Result<String, CamError> {
        format_g1(p, feed)
    }

    fn spindle_on(&self, rpm: f64) -> String {
        format!("M3 S{rpm:.0}")
    }

    fn spindle_off(&self) -> String {
        "M5".into()
    }

    fn tool_change(&self, _tool_id: u32) -> String {
        // GRBL has no tool change.
        String::new()
    }

    fn coolant_on(&self) -> String {
        "M8".into()
    }

    fn coolant_off(&self) -> String {
        "M9".into()
    }

    fn process(
        &self,
        toolpath: &Toolpath,
        tool: &Tool,
        spindle_rpm: f64,
    ) -> Result<String, CamError> {
        if toolpath.is_empty() {
            return Err(CamError::EmptyToolpath);
        }
        let mut out = String::new();
        out.push_str(&self.header(tool, spindle_rpm));
        for (idx, m) in toolpath.moves.iter().enumerate() {
            let line = match m.kind {
                MoveKind::Rapid => self.move_g0(m.position),
                MoveKind::Cut | MoveKind::Plunge => self.move_g1(m.position, m.feed)?,
                MoveKind::Arc { centre_xy, dir } => {
                    let start = if idx == 0 {
                        m.position
                    } else {
                        toolpath.moves[idx - 1].position
                    };
                    self.move_g23(start, m.position, centre_xy, dir, m.feed)?
                }
            };
            out.push_str(&line);
            out.push('\n');
        }
        out.push_str(&self.footer());
        Ok(out)
    }
}
