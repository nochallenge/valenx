//! Shared postprocessor scaffolding.
//!
//! Many controllers differ only in: header comment marker, units /
//! coordinate-mode preamble, tool-change syntax, spindle/coolant
//! M-codes, line numbering, and program-end syntax. This module
//! captures all of that in a single [`PostTemplate`] struct and
//! provides a [`process_template`] driver that the 27 variants
//! introduced in Phase 17C wire up via thin newtype wrappers.

use nalgebra::Vector3;

use crate::{
    error::CamError,
    post::{format_g0, format_g1, format_g23},
    tool::Tool,
    toolpath::{MoveKind, Toolpath},
};

/// Comment / line-number / dialect knobs shared by every templated
/// postprocessor.
#[derive(Clone, Debug)]
pub struct PostTemplate {
    /// Display name embedded into the header comment.
    pub name: &'static str,
    /// Character used for a single-line comment (`;`, `(`, `--`, etc.).
    /// `Open(c)` opens with `c`; `Pair(o, c)` wraps in `o ... c` (e.g.
    /// `(comment)` on Fanuc-family controllers).
    pub comment: CommentStyle,
    /// Program-start prologue ("%" on Fanuc, empty on GRBL).
    pub prelude: &'static str,
    /// Program ID line emitted right after the prelude (e.g. `O1000`).
    /// Empty string skips.
    pub program_id: &'static str,
    /// Units & coordinate mode header lines (e.g. `"G21\nG90\nG17"`).
    pub units_mode: &'static str,
    /// Tool-change template — `{n}` substitutes the tool id.
    /// Empty string ⇒ single-tool controller.
    pub tool_change_tpl: &'static str,
    /// Spindle-on template — `{rpm}` substitutes the integer RPM.
    pub spindle_on_tpl: &'static str,
    /// Spindle-off line.
    pub spindle_off: &'static str,
    /// Coolant-on line.
    pub coolant_on: &'static str,
    /// Coolant-off line.
    pub coolant_off: &'static str,
    /// Program-end line.
    pub program_end: &'static str,
    /// Program-end suffix (`"%\n"` on Fanuc, empty elsewhere).
    pub program_end_suffix: &'static str,
    /// `true` to prefix each move line with `N{n}` block numbers.
    pub number_lines: bool,
    /// Starting block number when `number_lines` is `true`.
    pub block_start: u32,
    /// Step between block numbers.
    pub block_step: u32,
}

/// Comment-syntax flavour.
#[derive(Clone, Copy, Debug)]
pub enum CommentStyle {
    /// Line comment opened by a single character (`;`, `#`).
    Open(char),
    /// Paired comment delimiter (e.g. `(comment)` on Fanuc).
    Pair(char, char),
    /// Custom prefix string (e.g. `";"`, `"//"`).
    Prefix(&'static str),
}

impl CommentStyle {
    /// Wrap `text` in the appropriate comment delimiters.
    ///
    /// Round-12 M4 (machine-safety): user-controlled strings reach
    /// this method via `Tool.name` (G-code header) and op labels. A
    /// hostile name like `"bit) G0 Z-99 ("` would break out of a
    /// `(...)` comment and inject a real motion command into the
    /// post-processed program — bad enough on a hobby CNC, lethal on
    /// a 5-axis machining centre. We scrub newlines (which would
    /// always escape a comment) and, for paired styles, the close
    /// delimiter (which would escape that specific style).
    pub fn wrap(self, text: &str) -> String {
        let cleaned = self.sanitize(text);
        match self {
            CommentStyle::Open(c) => format!("{c} {cleaned}"),
            CommentStyle::Pair(o, c) => format!("{o}{cleaned}{c}"),
            CommentStyle::Prefix(p) => format!("{p} {cleaned}"),
        }
    }

    /// Strip every character that could let `text` break out of the
    /// comment in this flavour:
    /// - `\n` / `\r` always (line-break escapes any comment style),
    /// - the close-delim character for [`CommentStyle::Pair`].
    fn sanitize(self, text: &str) -> String {
        let mut out = String::with_capacity(text.len());
        for ch in text.chars() {
            match ch {
                '\n' | '\r' => out.push(' '),
                c if matches!(self, CommentStyle::Pair(_, close) if c == close) => {
                    out.push('_');
                }
                c => out.push(c),
            }
        }
        out
    }
}

impl PostTemplate {
    /// Format the tool-change line per the template, returning an
    /// empty string for single-tool controllers.
    pub fn tool_change(&self, tool_id: u32) -> String {
        if self.tool_change_tpl.is_empty() {
            String::new()
        } else {
            self.tool_change_tpl.replace("{n}", &tool_id.to_string())
        }
    }

    /// Format the spindle-on line.
    pub fn spindle_on(&self, rpm: f64) -> String {
        self.spindle_on_tpl.replace("{rpm}", &format!("{rpm:.0}"))
    }

    /// Build the header block: prelude → program id → metadata
    /// comments → units → tool change → spindle on.
    pub fn header(&self, tool: &Tool, spindle_rpm: f64) -> String {
        let mut out = String::new();
        if !self.prelude.is_empty() {
            out.push_str(self.prelude);
            out.push('\n');
        }
        if !self.program_id.is_empty() {
            out.push_str(self.program_id);
            out.push('\n');
        }
        out.push_str(&self.comment.wrap(&format!("valenx-cam {}", self.name)));
        out.push('\n');
        out.push_str(&self.comment.wrap(&format!(
            "Tool: T{} {} D{:.2}mm",
            tool.id, tool.name, tool.diameter_mm
        )));
        out.push('\n');
        out.push_str(self.units_mode);
        out.push('\n');
        let tc = self.tool_change(tool.id);
        if !tc.is_empty() {
            out.push_str(&tc);
            out.push('\n');
        }
        out.push_str(&self.spindle_on(spindle_rpm));
        out.push('\n');
        out
    }

    /// Build the footer block: spindle off → program end → suffix.
    pub fn footer(&self) -> String {
        let mut out = String::new();
        out.push_str(self.spindle_off);
        out.push('\n');
        out.push_str(self.program_end);
        out.push('\n');
        if !self.program_end_suffix.is_empty() {
            out.push_str(self.program_end_suffix);
        }
        out
    }
}

/// Convert `toolpath` to G-code text using the supplied template.
///
/// 5-axis moves (Move with non-zero A/B in the future) are not
/// handled here — the templated postprocessors are 3-axis. Phase 17E
/// introduces a `move_g1_5ax` extension that defaults to a warning
/// comment line for templated posts.
pub fn process_template(
    tpl: &PostTemplate,
    toolpath: &Toolpath,
    tool: &Tool,
    spindle_rpm: f64,
) -> Result<String, CamError> {
    if toolpath.is_empty() {
        return Err(CamError::EmptyToolpath);
    }
    let mut out = String::new();
    out.push_str(&tpl.header(tool, spindle_rpm));
    let mut n = tpl.block_start;
    for (idx, m) in toolpath.moves.iter().enumerate() {
        // Negative feed encodes a dwell marker (peck_drill_full). We
        // emit a G4 line for those, where supported.
        if m.kind == MoveKind::Cut && m.feed < 0.0 {
            let seconds = -m.feed;
            let line = format!("G4 P{seconds:.3}");
            push_line(&mut out, tpl, &mut n, &line);
            continue;
        }
        let line = match m.kind {
            MoveKind::Rapid => format_g0(m.position),
            MoveKind::Cut | MoveKind::Plunge => format_g1(m.position, m.feed)?,
            MoveKind::Arc { centre_xy, dir } => {
                let start = if idx == 0 {
                    m.position
                } else {
                    toolpath.moves[idx - 1].position
                };
                format_g23(start, m.position, centre_xy, dir, m.feed)?
            }
        };
        push_line(&mut out, tpl, &mut n, &line);
    }
    out.push_str(&tpl.footer());
    Ok(out)
}

fn push_line(out: &mut String, tpl: &PostTemplate, n: &mut u32, line: &str) {
    if tpl.number_lines {
        out.push_str(&format!("N{n} {line}\n"));
        *n = n.saturating_add(tpl.block_step);
    } else {
        out.push_str(line);
        out.push('\n');
    }
}

/// Templated 5-axis fallback emitter — emits a warning comment for
/// posts that have no native A/B support. Returns the formatted line.
pub fn warn_five_axis_unsupported(
    tpl: &PostTemplate,
    _p: Vector3<f64>,
    _a: f64,
    _b: f64,
    _feed: f64,
) -> String {
    tpl.comment
        .wrap("5-axis move dropped — postprocessor lacks A/B output")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;
    use crate::toolpath::{Move, MoveKind};

    fn fanuc_like_tpl() -> PostTemplate {
        PostTemplate {
            name: "TestFanuc",
            comment: CommentStyle::Pair('(', ')'),
            prelude: "%",
            program_id: "O1000",
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
        }
    }

    #[test]
    fn header_and_footer_include_prelude_and_suffix() {
        let tpl = fanuc_like_tpl();
        let t = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let h = tpl.header(&t, 12000.0);
        assert!(h.starts_with("%\nO1000\n"), "header missing prelude: {h}");
        assert!(h.contains("(valenx-cam TestFanuc)"));
        assert!(h.contains("T1 M6"));
        assert!(h.contains("M3 S12000"));
        let f = tpl.footer();
        assert!(f.contains("M5"));
        assert!(f.ends_with("%\n"));
    }

    #[test]
    fn process_template_numbers_lines() {
        let tpl = fanuc_like_tpl();
        let t = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(10.0, 0.0, 0.0),
            500.0,
        ));
        let g = process_template(&tpl, &tp, &t, 12000.0).unwrap();
        assert!(g.contains("N10 G0"));
        assert!(g.contains("N20 G1"));
    }

    #[test]
    fn empty_toolpath_errors() {
        let tpl = fanuc_like_tpl();
        let t = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let tp = Toolpath::new();
        assert!(process_template(&tpl, &tp, &t, 12000.0).is_err());
    }

    #[test]
    fn dwell_encoded_as_g4() {
        let tpl = fanuc_like_tpl();
        let t = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Cut, Vector3::new(0.0, 0.0, 0.0), -0.5));
        let g = process_template(&tpl, &tp, &t, 12000.0).unwrap();
        assert!(g.contains("G4 P0.500"), "expected dwell marker, got:\n{g}");
    }

    #[test]
    fn comment_styles() {
        assert_eq!(CommentStyle::Open(';').wrap("hi"), "; hi");
        assert_eq!(CommentStyle::Pair('(', ')').wrap("hi"), "(hi)");
        assert_eq!(CommentStyle::Prefix("//").wrap("hi"), "// hi");
    }

    /// Round-12 M4 RED→GREEN: `CommentStyle::Pair` is used by every
    /// Fanuc-flavour post (Haas, Mach3, LinuxCNC, Tormach) and a
    /// hostile `Tool.name` like `")\nG0 Z-99\n("` would otherwise
    /// break out of the `(...)` comment and inject a real motion
    /// command into the G-code. The sanitiser replaces the close
    /// delimiter with `_` and newlines with spaces.
    #[test]
    fn paired_comment_neutralises_close_delim_and_newlines() {
        let attack = ")\nG0 Z-99\n(";
        let out = CommentStyle::Pair('(', ')').wrap(attack);
        // No newlines anywhere — they would always escape a comment
        // regardless of style.
        assert!(!out.contains('\n'), "newline leaked through: {out}");
        assert!(!out.contains('\r'), "carriage return leaked: {out}");
        // The outer wrap is `(...)` — outer `(` at index 0, outer `)`
        // at the last index. Both are legitimate parts of the wrap
        // itself. Carve out the inner content and assert it has no
        // close-delim left to escape the comment with.
        assert!(out.starts_with('('), "must still wrap with paren: {out}");
        assert!(out.ends_with(')'), "must still close with paren: {out}");
        let inner = &out[1..out.len() - 1];
        assert!(
            !inner.contains(')'),
            "inner text must not contain unescaped close-delim: {inner}"
        );
        // Pin the exact sanitised payload so a regression is obvious.
        // Newlines become spaces, `)` becomes `_`, and `(` passes
        // through unmolested (the open-delim is harmless inside a
        // paired comment — only the close escapes).
        assert_eq!(inner, "_ G0 Z-99 (");
    }

    /// Round-12 M4: open-style and prefix-style comments only need
    /// newline scrubbing — there's no close-delim to escape.
    #[test]
    fn open_and_prefix_comments_neutralise_newlines() {
        let attack = "first\nG0 Z-99";
        let out_open = CommentStyle::Open(';').wrap(attack);
        let out_prefix = CommentStyle::Prefix("//").wrap(attack);
        assert!(
            !out_open.contains('\n'),
            "; comment leaked newline: {out_open}"
        );
        assert!(
            !out_prefix.contains('\n'),
            "// comment leaked newline: {out_prefix}"
        );
    }
}
