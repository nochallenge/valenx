//! Ergonomic RS-274 G-code emission from a [`Toolpath`].
//!
//! The [`crate::post`] module hosts the full controller-specific
//! postprocessor family (GRBL / LinuxCNC / Fanuc / Haas / …), each of
//! which needs a [`Tool`](crate::tool::Tool) and a spindle RPM. This
//! module is the *minimal* counterpart: a single in-house function,
//! [`to_gcode`], that turns a bare [`Toolpath`] into a standard
//! ([RS-274 / ISO 6983]) G-code string with a sane header and footer
//! and **no** mandatory tool/spindle metadata. It is the path most
//! callers want when they simply have a toolpath and need plain
//! controller-agnostic text:
//!
//! ```text
//! let nc = valenx_cam::gcode::to_gcode(&toolpath);
//! ```
//!
//! ## What it emits
//!
//! - A header: `%` start, a `; valenx-cam` banner comment, then the
//!   modal setup block `G21 G90 G94 G17` (millimetres, absolute,
//!   feed-per-minute, XY plane).
//! - Optional spindle start (`M3 S{rpm}`) when an RPM is supplied.
//! - One line per [`Move`]:
//!   - [`MoveKind::Rapid`] → `G0 X Y Z`
//!   - [`MoveKind::Cut`] / [`MoveKind::Plunge`] → `G1 X Y Z F{feed}`
//!   - [`MoveKind::Arc`] → `G2`/`G3 X Y Z I J F{feed}`
//! - A footer: spindle off (`M5`, if it was started), program end
//!   (`M30`), and the `%` terminator.
//!
//! All coordinate / feed formatting (and the NaN / ±inf guards that
//! stop a corrupt toolpath from emitting a `F inf` "unlimited rapid")
//! is shared with the postprocessors via [`crate::post`]'s
//! `format_g0` / `format_g1` / `format_g23` helpers — so `to_gcode`
//! and a postprocessor always agree on the move syntax.
//!
//! This is a fully **in-house** generator. The `gcode` crate on
//! crates.io *parses* G-code; it does not generate it, so it cannot be
//! used for the emission half. The round-trip tests re-parse the
//! emitted text with a small in-tree tokenizer to prove every cut
//! move survives a generate → parse round trip with its coordinates
//! intact.
//!
//! [RS-274 / ISO 6983]: https://en.wikipedia.org/wiki/G-code

use crate::{
    error::CamError,
    post::{format_g0, format_g1, format_g23},
    toolpath::{Move, MoveKind, Toolpath},
};

/// Options controlling [`to_gcode_with`] output. Defaults to a
/// metric, absolute, spindle-less program — call [`GcodeOptions::rpm`]
/// to add a spindle start/stop pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GcodeOptions {
    /// Spindle speed in RPM. `Some(rpm)` emits `M3 S{rpm}` in the
    /// header and `M5` in the footer; `None` omits both (the toolpath
    /// is assumed to run under an externally-commanded spindle, e.g.
    /// a laser or a pre-started job).
    pub spindle_rpm: Option<f64>,
    /// When `true`, emit the `%` program start/stop sentinels that
    /// tape-style controllers (Fanuc, Haas) require. GRBL and most
    /// hobby firmware ignore them harmlessly. Default `true`.
    pub program_markers: bool,
}

impl Default for GcodeOptions {
    fn default() -> Self {
        Self {
            spindle_rpm: None,
            program_markers: true,
        }
    }
}

impl GcodeOptions {
    /// Builder: set the spindle RPM (emits `M3 S{rpm}` / `M5`).
    #[must_use]
    pub fn rpm(mut self, rpm: f64) -> Self {
        self.spindle_rpm = Some(rpm);
        self
    }

    /// Builder: toggle the `%` program-start/stop markers.
    #[must_use]
    pub fn program_markers(mut self, on: bool) -> Self {
        self.program_markers = on;
        self
    }
}

/// Emit standard RS-274 G-code for `toolpath` with default options
/// (metric, absolute, no spindle command).
///
/// Returns `; (empty toolpath)\n` for an empty toolpath rather than
/// erroring — the result is always a syntactically valid (if inert)
/// program. Callers that want a hard error on empty input should use
/// [`to_gcode_checked`].
///
/// # Example
/// ```
/// use nalgebra::Vector3;
/// use valenx_cam::gcode::to_gcode;
/// use valenx_cam::toolpath::{Move, MoveKind, Toolpath};
///
/// let mut tp = Toolpath::new();
/// tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
/// tp.push(Move::new(MoveKind::Cut, Vector3::new(10.0, 0.0, 5.0), 600.0));
/// let nc = to_gcode(&tp);
/// assert!(nc.contains("G1 X10.000 Y0.000 Z5.000 F600"));
/// ```
#[must_use]
pub fn to_gcode(toolpath: &Toolpath) -> String {
    // Default options never trigger the feed/coordinate validation
    // failure path for well-formed toolpaths; on a corrupt move the
    // checked emitter would error, but here we degrade to a comment.
    match to_gcode_with(toolpath, GcodeOptions::default()) {
        Ok(s) => s,
        Err(e) => format!("; (G-code emission failed: {e})\n"),
    }
}

/// Like [`to_gcode`] but returns [`CamError::EmptyToolpath`] for an
/// empty path and propagates any [`CamError::PostprocessorFailed`]
/// raised by a non-finite coordinate / feed in a move.
pub fn to_gcode_checked(toolpath: &Toolpath) -> Result<String, CamError> {
    if toolpath.is_empty() {
        return Err(CamError::EmptyToolpath);
    }
    to_gcode_with(toolpath, GcodeOptions::default())
}

/// Emit RS-274 G-code for `toolpath` under explicit [`GcodeOptions`].
///
/// An empty toolpath yields a valid no-op program (header + footer
/// only). A move with a non-finite cut feed or coordinate returns
/// [`CamError::PostprocessorFailed`] (rapids degrade their bad
/// coordinates to `0.000` via the internal `format_g0` helper, matching
/// postprocessor behaviour).
pub fn to_gcode_with(toolpath: &Toolpath, opts: GcodeOptions) -> Result<String, CamError> {
    let mut out = String::new();

    // ---- header ----
    if opts.program_markers {
        out.push_str("%\n");
    }
    out.push_str("; valenx-cam G-code (RS-274)\n");
    // G21 mm · G90 absolute · G94 feed-per-minute · G17 XY plane.
    out.push_str("G21 G90 G94 G17\n");
    if let Some(rpm) = opts.spindle_rpm {
        if !rpm.is_finite() || rpm < 0.0 {
            return Err(CamError::PostprocessorFailed {
                reason: format!("spindle rpm must be a finite, non-negative number (got {rpm})"),
            });
        }
        out.push_str(&format!("M3 S{rpm:.0}\n"));
    }

    // ---- moves ----
    for (idx, m) in toolpath.moves.iter().enumerate() {
        let line = format_move(toolpath, idx, m)?;
        out.push_str(&line);
        out.push('\n');
    }

    // ---- footer ----
    if opts.spindle_rpm.is_some() {
        out.push_str("M5\n");
    }
    out.push_str("M30\n");
    if opts.program_markers {
        out.push_str("%\n");
    }
    Ok(out)
}

/// Format a single move to its G-code line, resolving the arc start
/// from the previous move when needed.
fn format_move(toolpath: &Toolpath, idx: usize, m: &Move) -> Result<String, CamError> {
    Ok(match m.kind {
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
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn p(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    /// A closed rectangular contour at Z = depth: rapid in above the
    /// start corner, plunge to depth, cut the four sides back to the
    /// start, then rapid clear. Feed in mm/min.
    fn rect_contour(w: f64, h: f64, depth: f64, feed: f64) -> Toolpath {
        let mut t = Toolpath::new();
        t.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        t.push(Move::new(MoveKind::Plunge, p(0.0, 0.0, depth), feed));
        t.push(Move::new(MoveKind::Cut, p(w, 0.0, depth), feed));
        t.push(Move::new(MoveKind::Cut, p(w, h, depth), feed));
        t.push(Move::new(MoveKind::Cut, p(0.0, h, depth), feed));
        t.push(Move::new(MoveKind::Cut, p(0.0, 0.0, depth), feed)); // close
        t.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        t
    }

    /// A minimal RS-274 line parser sufficient to re-read what
    /// `to_gcode` emits: pulls the leading G/M word plus any
    /// X/Y/Z/I/J/F address values. Proves the emitter's output is
    /// re-parseable (a generate → parse round trip) without depending
    /// on an external crate (the crates.io `gcode` crate parses but
    /// does not generate, so it can validate but cannot produce our
    /// output).
    #[derive(Debug, Default, PartialEq)]
    struct ParsedLine {
        word: String,
        x: Option<f64>,
        y: Option<f64>,
        z: Option<f64>,
        f: Option<f64>,
    }

    fn parse_line(line: &str) -> Option<ParsedLine> {
        // Strip a trailing `; comment`.
        let code = line.split(';').next().unwrap_or("").trim();
        if code.is_empty() || code == "%" {
            return None;
        }
        let mut pl = ParsedLine::default();
        for (i, tok) in code.split_whitespace().enumerate() {
            let letter = tok.chars().next().unwrap();
            let rest = &tok[letter.len_utf8()..];
            match letter {
                'G' | 'M' if i == 0 => pl.word = tok.to_string(),
                'X' => pl.x = rest.parse().ok(),
                'Y' => pl.y = rest.parse().ok(),
                'Z' => pl.z = rest.parse().ok(),
                'F' => pl.f = rest.parse().ok(),
                _ => {}
            }
        }
        if pl.word.is_empty() {
            None
        } else {
            Some(pl)
        }
    }

    #[test]
    fn emits_header_and_footer() {
        let nc = to_gcode(&rect_contour(10.0, 20.0, -2.0, 600.0));
        assert!(nc.starts_with("%\n"), "missing program start marker:\n{nc}");
        assert!(nc.contains("; valenx-cam"), "missing banner:\n{nc}");
        assert!(nc.contains("G21 G90 G94 G17"), "missing modal setup:\n{nc}");
        assert!(nc.contains("M30"), "missing program end:\n{nc}");
        assert!(nc.trim_end().ends_with('%'), "missing end marker:\n{nc}");
    }

    #[test]
    fn rectangular_contour_has_correct_g1_moves_and_coords() {
        let nc = to_gcode(&rect_contour(10.0, 20.0, -2.0, 600.0));
        // The four cut sides + the closing plunge are all G1 lines; the
        // four side end-points must appear verbatim with .3f coords and
        // the F600 feed.
        for expected in [
            "G1 X10.000 Y0.000 Z-2.000 F600",  // bottom edge
            "G1 X10.000 Y20.000 Z-2.000 F600", // right edge
            "G1 X0.000 Y20.000 Z-2.000 F600",  // top edge
            "G1 X0.000 Y0.000 Z-2.000 F600",   // closing edge
        ] {
            assert!(nc.contains(expected), "missing `{expected}` in:\n{nc}");
        }
        // The retract / approach are rapids, not feeds.
        assert!(
            nc.contains("G0 X0.000 Y0.000 Z5.000"),
            "missing rapid:\n{nc}"
        );
        // Exactly five G1 lines (one plunge + four cuts).
        let g1_count = nc
            .lines()
            .filter(|l| l.trim_start().starts_with("G1 "))
            .count();
        assert_eq!(g1_count, 5, "expected 5 G1 lines, got {g1_count}:\n{nc}");
        // Exactly two rapids (approach + retract).
        let g0_count = nc
            .lines()
            .filter(|l| l.trim_start().starts_with("G0 "))
            .count();
        assert_eq!(g0_count, 2, "expected 2 G0 lines, got {g0_count}:\n{nc}");
    }

    #[test]
    fn round_trips_through_parser() {
        let depth = -2.0;
        let feed = 600.0;
        let tp = rect_contour(10.0, 20.0, depth, feed);
        let nc = to_gcode(&tp);

        // Re-parse every motion line and pair it back with the source
        // move; coordinates and feed must survive the round trip.
        let parsed: Vec<ParsedLine> = nc.lines().filter_map(parse_line).collect();
        let motion: Vec<&ParsedLine> = parsed
            .iter()
            .filter(|p| p.word == "G0" || p.word == "G1" || p.word == "G2" || p.word == "G3")
            .collect();
        assert_eq!(
            motion.len(),
            tp.moves.len(),
            "motion-line count must match move count"
        );

        for (mv, pl) in tp.moves.iter().zip(motion.iter()) {
            // Word matches the move kind.
            let expect_word = match mv.kind {
                MoveKind::Rapid => "G0",
                MoveKind::Cut | MoveKind::Plunge => "G1",
                MoveKind::Arc { .. } => unreachable!("contour has no arcs"),
            };
            assert_eq!(pl.word, expect_word, "wrong word for {mv:?}");
            // Coordinates round-trip to 1e-3 (the emission precision).
            assert!(
                (pl.x.unwrap() - mv.position.x).abs() < 1e-3,
                "x mismatch: {pl:?} vs {mv:?}"
            );
            assert!(
                (pl.y.unwrap() - mv.position.y).abs() < 1e-3,
                "y mismatch: {pl:?} vs {mv:?}"
            );
            assert!(
                (pl.z.unwrap() - mv.position.z).abs() < 1e-3,
                "z mismatch: {pl:?} vs {mv:?}"
            );
            // Feed only present on G1.
            if expect_word == "G1" {
                assert!((pl.f.unwrap() - feed).abs() < 1e-6, "feed mismatch: {pl:?}");
            } else {
                assert!(pl.f.is_none(), "rapid should not carry feed: {pl:?}");
            }
        }
    }

    #[test]
    fn spindle_option_emits_m3_m5() {
        let tp = rect_contour(10.0, 10.0, -1.0, 300.0);
        let nc = to_gcode_with(&tp, GcodeOptions::default().rpm(12000.0)).unwrap();
        assert!(nc.contains("M3 S12000"), "missing spindle on:\n{nc}");
        assert!(nc.contains("M5"), "missing spindle off:\n{nc}");
        // Default (no rpm) omits both.
        let plain = to_gcode(&tp);
        assert!(
            !plain.contains("M3 S"),
            "plain output must not start spindle"
        );
    }

    #[test]
    fn arc_move_emits_g2_g3_with_ij() {
        use crate::arcfit::ArcDir;
        let mut t = Toolpath::new();
        t.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 0.0), 0.0));
        // Quarter-circle CW from (0,0) to (1,1) about centre (1,0):
        // I = 1-0 = 1, J = 0-0 = 0.
        t.push(Move::new(
            MoveKind::Arc {
                centre_xy: nalgebra::Vector2::new(1.0, 0.0),
                dir: ArcDir::Clockwise,
            },
            p(1.0, 1.0, 0.0),
            400.0,
        ));
        let nc = to_gcode(&t);
        assert!(
            nc.contains("G2 X1.000 Y1.000 Z0.000 I1.000 J0.000 F400"),
            "arc line wrong:\n{nc}"
        );
    }

    #[test]
    fn program_markers_toggle() {
        let tp = rect_contour(5.0, 5.0, -1.0, 200.0);
        let nc = to_gcode_with(&tp, GcodeOptions::default().program_markers(false)).unwrap();
        assert!(!nc.contains('%'), "markers should be suppressed:\n{nc}");
        assert!(nc.contains("G1 X5.000 Y0.000 Z-1.000 F200"));
    }

    #[test]
    fn empty_toolpath_is_valid_noop_but_checked_errors() {
        let empty = Toolpath::new();
        // Lenient emitter: still a syntactically valid program.
        let nc = to_gcode(&empty);
        assert!(nc.contains("G21 G90 G94 G17"));
        assert!(nc.contains("M30"));
        // No motion lines.
        assert_eq!(
            nc.lines()
                .filter(|l| l.starts_with("G0") || l.starts_with("G1"))
                .count(),
            0
        );
        // Strict emitter: errors.
        assert!(matches!(
            to_gcode_checked(&empty),
            Err(CamError::EmptyToolpath)
        ));
    }

    #[test]
    fn nonfinite_cut_feed_is_rejected_by_checked() {
        let mut t = Toolpath::new();
        t.push(Move::new(MoveKind::Cut, p(1.0, 0.0, 0.0), f64::INFINITY));
        let err = to_gcode_checked(&t).expect_err("inf feed must error");
        assert!(matches!(err, CamError::PostprocessorFailed { .. }));
        assert_eq!(err.code(), "cam.postprocessor_failed");
    }
}
