//! HeeksCNC-compatible `.nc` writer. Writes a G-code variant that
//! mirrors HeeksCAD's post-processor: metric, absolute coordinates,
//! G0/G1 with per-move feed-rate annotation, and a `; valenx-heekscad`
//! header comment.

use std::path::Path;

use crate::cam::{Move, Toolpath};
use crate::error::HeeksCadError;

/// Write the G-code for `toolpath` to a string.
pub fn write_heeks_string(toolpath: &Toolpath) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "; valenx-heekscad op `{}` tool D={:.3}mm feed={:.1}\n",
        toolpath.op_name, toolpath.tool.diameter, toolpath.tool.feed_rate
    ));
    s.push_str("G21\nG90\nG17\n");
    s.push_str(&format!("G0 Z{:.4}\n", toolpath.clearance_z));
    for m in &toolpath.moves {
        match m {
            Move::Rapid { x, y } => s.push_str(&format!("G0 X{x:.4} Y{y:.4}\n")),
            Move::Plunge { z } => {
                s.push_str(&format!("G1 Z{z:.4} F{:.1}\n", toolpath.tool.plunge_rate));
            }
            Move::Retract { z } => s.push_str(&format!("G0 Z{z:.4}\n")),
            Move::Feed { x, y } => s.push_str(&format!(
                "G1 X{x:.4} Y{y:.4} F{:.1}\n",
                toolpath.tool.feed_rate
            )),
        }
    }
    s.push_str("M30\n");
    s
}

/// Write the G-code for `toolpath` to `path` on disk.
pub fn write_heeks(toolpath: &Toolpath, path: &Path) -> Result<(), HeeksCadError> {
    valenx_core::io_caps::atomic_write_str(path, &write_heeks_string(toolpath))
        .map_err(|e| HeeksCadError::Io(e.to_string()))
}
