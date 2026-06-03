//! Shared egui helpers for the Genetics workbench panels.
//!
//! The 13 bio-crate panels all follow the same shape — a tool palette
//! header, input forms (multiline sequence / numeric drag-values /
//! file pickers), a Run button and a results area. These helpers keep
//! that idiom consistent and the per-panel modules short.
//!
//! The polish pass extended this module with:
//!
//! - **Icon-prefixed Run / Undo / Redo buttons** — the `▶ Run` /
//!   `↶ Undo` / `↷ Redo` cadence comes from `valenx_icons` so every
//!   panel reads the same shape.
//! - **`undo_redo_inline`** — a one-call helper to render the small
//!   `↶ ↷` button pair next to an input.
//! - **`friendly_error`** — a name-based recovery-hint generator for
//!   the broader set of bio panels. The Sequence / Aero panels keep
//!   their bespoke mappers (which carry domain knowledge), but every
//!   other panel can call this one for a "what went wrong + what to
//!   try" line in one call.

use eframe::egui;
use valenx_icons::{edit, run as icon_run, status};

/// Render a panel section header (matches `mesh_toolbox`'s strong-text
/// section labels).
pub fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(2.0);
    ui.label(egui::RichText::new(title).strong());
}

/// A monospace, read-only, multi-line results box. Used for alignment
/// renderings, dot-bracket strings, Newick text, ASCII trees, etc.
pub fn mono_output(ui: &mut egui::Ui, id: &str, text: &str, rows: usize) {
    let mut buf = text.to_string();
    ui.add(
        egui::TextEdit::multiline(&mut buf)
            .id_source(id)
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(rows)
            .interactive(false),
    );
}

/// A labelled error line — red, only drawn when `err` is `Some`.
///
/// Uses the canonical `accent::ERROR` colour from `valenx-design-tokens`
/// so the rendered red matches whatever theme variant is active
/// (Dark / Light / High-Contrast).
///
/// Adds an auto-generated recovery hint when the error message matches
/// a known pattern (see [`friendly_error`]). Panels with bespoke
/// mappers (Sequence's `friendly_error`, Aero's `friendly_aero_error`)
/// still surface their domain-specific text — the bare `e` is passed
/// through unchanged when no pattern matches.
pub fn error_line(ui: &mut egui::Ui, err: &Option<String>) {
    if let Some(e) = err {
        // Show the message itself.
        ui.colored_label(error_color(), format!("{} {e}", status::ERROR));
        // Plus a heuristic recovery hint (silent if no pattern matched).
        let (_, hint) = friendly_error(e, "");
        if !hint.is_empty() {
            ui.weak(format!("→ {hint}"));
        }
    }
}

/// Like [`error_line`] but with a follow-up "recovery hint" weak
/// line. The intent of the polish pass: every error tells the user
/// what to try next so they're not staring at a raw `Err(...)`
/// debug-print.
pub fn error_line_with_hint(ui: &mut egui::Ui, err: &Option<String>, hint: &str) {
    if let Some(e) = err {
        ui.colored_label(error_color(), format!("{} {e}", status::ERROR));
        if !hint.is_empty() {
            ui.weak(format!("→ {hint}"));
        }
    }
}

/// Map a raw `Err(...)` from a Run action into a friendly message +
/// recovery hint. The pattern is the same one `sequence::friendly_error`
/// + `aero::friendly_aero_error` follow, but reusable from any
///   panel that doesn't have its own bespoke mapper.
///
/// Recognised patterns (case-insensitive substring match on `raw`):
///
/// - `"empty"` / `"no input"` — "needs a non-empty input".
/// - `"invalid"` / `"unsupported"` — "input has unsupported characters".
/// - `"not found"` / `"unknown"` — "couldn't find a referenced item".
/// - `"timeout"` / `"timed out"` — "took too long, try smaller input".
/// - `"out of memory"` / `"oom"` — "input too large, narrow the inputs".
/// - `"converge"` — "didn't converge — try a different parameter range".
///
/// Anything we don't recognise passes through unchanged so we never
/// lose information.
pub fn friendly_error(raw: &str, action: &str) -> (String, String) {
    let low = raw.to_ascii_lowercase();
    if low.contains("empty") || low.contains("no input") {
        (
            format!("{action} needs a non-empty input."),
            "Paste a sequence / select a file / fill the form, then try again.".to_string(),
        )
    } else if low.contains("invalid") || low.contains("unsupported") {
        (
            format!("{action}: {raw}"),
            "Check the input characters — only the documented alphabet is accepted.".to_string(),
        )
    } else if low.contains("not found") || low.contains("unknown") {
        (
            format!("{action}: {raw}"),
            "Check spelling — names are case-sensitive.".to_string(),
        )
    } else if low.contains("timeout") || low.contains("timed out") {
        (
            format!("{action} timed out."),
            "Reduce the input size or tighten the iteration count, then retry.".to_string(),
        )
    } else if low.contains("out of memory") || low.contains("oom") {
        (
            format!("{action} ran out of memory."),
            "The input is too large — narrow the window / split into chunks.".to_string(),
        )
    } else if low.contains("converge") {
        (
            format!("{action} did not converge: {raw}"),
            "Try a different parameter range, more iterations, or a coarser tolerance.".to_string(),
        )
    } else if low.contains("parse") || low.contains("malformed") {
        (
            format!("{action}: {raw}"),
            "Check the input format — see the panel help (F1) for examples.".to_string(),
        )
    } else {
        // Unrecognised — pass the message through with a generic hint.
        (format!("{action}: {raw}"), String::new())
    }
}

/// Token-driven error colour — `accent::ERROR` from
/// `valenx-design-tokens` converted to egui's `Color32`. The previous
/// pre-polish hex literal (`#E5787C`) drifted from the token-defined
/// palette over time; this routes back through the canonical source.
fn error_color() -> egui::Color32 {
    let hex = valenx_design_tokens::color::accent::ERROR;
    egui::Color32::from_rgb(
        ((hex >> 16) & 0xFF) as u8,
        ((hex >> 8) & 0xFF) as u8,
        (hex & 0xFF) as u8,
    )
}

/// Token-driven success / "OK" colour. Mirrors [`error_color`] for
/// the green half of the success / error pair.
fn success_color() -> egui::Color32 {
    let hex = valenx_design_tokens::color::accent::SUCCESS;
    egui::Color32::from_rgb(
        ((hex >> 16) & 0xFF) as u8,
        ((hex >> 8) & 0xFF) as u8,
        (hex & 0xFF) as u8,
    )
}

/// Run-button helper that decorates the standard button with the
/// canonical Ctrl+R / F5 tooltip + the play-glyph icon. Replaces
/// bare `common::run_button` at the call sites that have been
/// polish-passed.
pub fn run_button_with_tt(ui: &mut egui::Ui, label: &str, tooltip: &str) -> bool {
    let combined = format!("{tooltip}\nShortcut: Ctrl+R / F5");
    ui.add_sized(
        [ui.available_width(), 24.0],
        egui::Button::new(format!("{} {label}", icon_run::PLAY)),
    )
    .on_hover_text(combined)
    .clicked()
}

/// A green "ok / done" status line.
pub fn ok_line(ui: &mut egui::Ui, msg: &str) {
    ui.colored_label(success_color(), format!("{} {msg}", status::SUCCESS));
}

/// A primary "Run" button — slightly wider so it reads as the panel's
/// main action. Returns `true` on click.
///
/// All Run buttons in the Genetics workbench share the Ctrl+R / F5
/// keyboard hint via this helper's hover-tooltip — added in the
/// frontend-polish pass so every panel surfaces the shortcut without
/// each call-site having to remember.
///
/// The label is prefixed with the [`valenx_icons::run::PLAY`] glyph
/// so the Run button is visually consistent across every panel.
pub fn run_button(ui: &mut egui::Ui, label: &str) -> bool {
    ui.add_sized(
        [ui.available_width(), 24.0],
        egui::Button::new(format!("{} {label}", icon_run::PLAY)),
    )
    .on_hover_text("Execute the active sub-tool's primary action.\nShortcut: Ctrl+R / F5")
    .clicked()
}

/// Render the inline undo / redo button pair next to an input.
///
/// Returns `(undo_clicked, redo_clicked)` — the caller calls
/// `panel.undo_edit()` / `panel.redo_edit()` on a `true` to drive
/// the panel's `History<T>` snapshot stack. The buttons share the
/// same `↶` / `↷` glyphs from `valenx_icons::edit` so every panel
/// presents the same affordance.
///
/// `can_undo` / `can_redo` are the panel's `History::can_undo()` /
/// `History::can_redo()` — they disable the button when nothing's
/// left to undo / redo so the user gets immediate feedback that
/// they've reached the end of the timeline.
pub fn undo_redo_inline(
    ui: &mut egui::Ui,
    can_undo: bool,
    can_redo: bool,
) -> (bool, bool) {
    let undo_clicked = ui
        .add_enabled(can_undo, egui::Button::new(edit::UNDO))
        .on_hover_text("Undo last edit (Ctrl+Z)")
        .clicked();
    let redo_clicked = ui
        .add_enabled(can_redo, egui::Button::new(edit::REDO))
        .on_hover_text("Redo last edit (Ctrl+Y / Ctrl+Shift+Z)")
        .clicked();
    (undo_clicked, redo_clicked)
}

/// Render a small "Undo / Redo" row right under the section header
/// of a panel, calling back into the panel's snapshot stack on click.
///
/// Returns `(undo_clicked, redo_clicked)`. Callers route those through
/// `panel.undo_edit()` / `panel.redo_edit()`. Pairs with
/// [`undo_redo_inline`] (which is the wrapped-horizontal variant
/// callers paste into an existing `ui.horizontal(|ui| { ... })`).
pub fn undo_redo_row(
    ui: &mut egui::Ui,
    can_undo: bool,
    can_redo: bool,
) -> (bool, bool) {
    let mut undo_clicked = false;
    let mut redo_clicked = false;
    ui.horizontal(|ui| {
        let (u, r) = undo_redo_inline(ui, can_undo, can_redo);
        undo_clicked = u;
        redo_clicked = r;
        ui.label(
            egui::RichText::new("Ctrl+Z / Ctrl+Y")
                .weak()
                .small(),
        );
    });
    (undo_clicked, redo_clicked)
}

/// Trait every editable-state panel implements so the host's
/// `Ctrl+Z` / `Ctrl+Y` dispatcher can drive undo / redo uniformly.
///
/// The 3 + 11 + … panels each own a `History<Snapshot>` over their
/// own snapshot type. This trait abstracts that so the host doesn't
/// have to `match` on every panel variant.
///
/// Implementations are mechanical — see
/// [`crate::genetics::sequence::SequencePanel`] /
/// [`crate::genetics::alignment::AlignmentPanel`] for the canonical
/// shape (snapshot + restore + record-on-mutate).
pub trait PanelHistory {
    /// Undo the most recent edit. Returns `true` if a snapshot was
    /// restored.
    fn undo_edit(&mut self) -> bool;
    /// Redo the most recently undone edit. Returns `true` if a
    /// snapshot was restored.
    fn redo_edit(&mut self) -> bool;
    /// Whether an undo would do something.
    fn can_undo(&self) -> bool;
    /// Whether a redo would do something.
    fn can_redo(&self) -> bool;
}

/// A two-column key/value row inside a grid.
pub fn kv(ui: &mut egui::Ui, key: &str, value: impl Into<String>) {
    ui.label(key);
    ui.label(value.into());
    ui.end_row();
}

/// A multiline sequence input with a label above it.
///
/// The text-edit carries a hover tooltip explaining the input
/// conventions (FASTA-friendly, IUPAC codes, whitespace stripped)
/// so users don't have to guess what the panel will accept.
pub fn seq_input(ui: &mut egui::Ui, id: &str, label: &str, buf: &mut String, rows: usize) {
    ui.label(label);
    ui.add(
        egui::TextEdit::multiline(buf)
            .id_source(id)
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(rows)
            .hint_text("paste a sequence here"),
    )
    .on_hover_text(
        "Paste a sequence as bare residue letters or as a FASTA record. \
         Whitespace, digits, and FASTA header lines (starting with `>` or `;`) \
         are stripped before validation.",
    );
}

/// Strip whitespace / FASTA headers / digits from a pasted sequence,
/// returning the bare residue letters uppercased. Lines beginning with
/// `>` or `;` (FASTA / comment) are dropped — so a user can paste a
/// whole FASTA record and the panel uses the first sequence body.
pub fn clean_sequence(raw: &str) -> String {
    raw.lines()
        .filter(|l| !l.trim_start().starts_with('>') && !l.trim_start().starts_with(';'))
        .flat_map(|l| l.chars())
        .filter(|c| c.is_ascii_alphabetic() || *c == '*' || *c == '-')
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

/// Parse a multi-FASTA blob into `(label, sequence)` pairs. A run of
/// lines with no leading `>` before the first header is taken as one
/// unnamed record. Sequence bodies are whitespace-stripped but NOT
/// case-folded (callers fold as needed).
pub fn parse_fasta(raw: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    let mut label = String::new();
    let mut body = String::new();
    let mut seen_header = false;
    for line in raw.lines() {
        let t = line.trim_end();
        if let Some(rest) = t.strip_prefix('>') {
            if seen_header || !body.is_empty() {
                out.push((std::mem::take(&mut label), std::mem::take(&mut body)));
            }
            label = rest.trim().to_string();
            seen_header = true;
        } else {
            body.extend(t.chars().filter(|c| !c.is_whitespace()));
        }
    }
    if seen_header || !body.is_empty() {
        out.push((label, body));
    }
    out.into_iter().filter(|(_, s)| !s.is_empty()).collect()
}
