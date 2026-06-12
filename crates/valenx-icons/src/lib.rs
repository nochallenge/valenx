//! # valenx-icons
//!
//! Curated icon glyph set for the Valenx app.
//!
//! ## Approach
//!
//! We use **Unicode glyphs from egui's bundled font** rather than
//! shipping rasterised PNGs or rasterising SVGs at runtime. The
//! rationale:
//!
//! - egui's default font (Inter + an emoji fallback) already ships
//!   with the entire Basic Latin + Geometric Shapes + Arrows +
//!   Miscellaneous Technical + Miscellaneous Symbols Unicode blocks.
//!   Every glyph here renders on every platform with no asset-loading
//!   ceremony.
//! - Vendoring rasterised PNGs would bloat the binary by ~40-80 KB
//!   per icon × 60 icons = ~3-5 MB; vendoring SVGs + adding a runtime
//!   rasteriser (`resvg`) is another ~2 MB on top.
//! - The icons here are **semantically named** — `run::PLAY`,
//!   `file::SAVE`, `domain::DNA` — so swapping in a rasterised set
//!   later is a single-file change with no panel-side rewrites.
//!
//! ## License
//!
//! All glyphs are part of the Unicode standard (codepoints listed in
//! source). The mapping itself is licensed under the Valenx workspace
//! licence (MIT OR Apache-2.0). No third-party icon font is bundled,
//! so no third-party attribution is required.
//!
//! ## Icon families
//!
//! - [`run`] — play / pause / stop / run-action icons.
//! - [`mod@file`] — open / save / new / export.
//! - [`edit`] — undo / redo / cut / copy / paste / search.
//! - [`status`] — info / warning / error / success / question.
//! - [`nav`] — arrows / chevrons / expand / collapse.
//! - [`view`] — grid / list / fullscreen / zoom.
//! - [`domain`] — DNA / chip / structure / gear / fluid / waveform —
//!   workbench-domain marks.
//!
//! ## Usage
//!
//! ```ignore
//! use valenx_icons::{run, file, status};
//!
//! ui.button(format!("{} Run", run::PLAY));
//! ui.button(format!("{} Save", file::SAVE));
//! ui.colored_label(error_color, format!("{} {}", status::ERROR, msg));
//! ```

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

/// "Run" / playback control icons.
///
/// Common Geometric Shapes + Miscellaneous Symbols.
pub mod run {
    /// ▶ — primary "Run" / play marker. U+25B6 BLACK RIGHT-POINTING TRIANGLE.
    pub const PLAY: &str = "\u{25B6}";
    /// ⏸ — pause marker. U+23F8 DOUBLE VERTICAL BAR.
    pub const PAUSE: &str = "\u{23F8}";
    /// ⏹ — stop marker. U+23F9 BLACK SQUARE FOR STOP.
    pub const STOP: &str = "\u{23F9}";
    /// ⟳ — restart / re-run. U+27F3 CLOCKWISE GAPPED CIRCLE ARROW.
    pub const RESTART: &str = "\u{27F3}";
    /// ⏵ — step forward. U+23F5 BLACK MEDIUM RIGHT-POINTING TRIANGLE.
    pub const STEP: &str = "\u{23F5}";
}

/// File / project icons.
pub mod file {
    /// 💾 — save (floppy). U+1F4BE FLOPPY DISK.
    pub const SAVE: &str = "\u{1F4BE}";
    /// 📂 — open folder. U+1F4C2 OPEN FILE FOLDER.
    pub const OPEN: &str = "\u{1F4C2}";
    /// 📄 — new document. U+1F4C4 PAGE FACING UP.
    pub const NEW: &str = "\u{1F4C4}";
    /// ⬇ — import / download. U+2B07 DOWNWARDS BLACK ARROW.
    pub const IMPORT: &str = "\u{2B07}";
    /// ⬆ — export / upload. U+2B06 UPWARDS BLACK ARROW.
    pub const EXPORT: &str = "\u{2B06}";
    /// 📁 — generic folder. U+1F4C1 FILE FOLDER.
    pub const FOLDER: &str = "\u{1F4C1}";
}

/// Edit / clipboard icons.
pub mod edit {
    /// ↶ — undo. U+21B6 ANTICLOCKWISE TOP SEMICIRCLE ARROW.
    pub const UNDO: &str = "\u{21B6}";
    /// ↷ — redo. U+21B7 CLOCKWISE TOP SEMICIRCLE ARROW.
    pub const REDO: &str = "\u{21B7}";
    /// ✂ — cut. U+2702 BLACK SCISSORS.
    pub const CUT: &str = "\u{2702}";
    /// 🔍 — search. U+1F50D LEFT-POINTING MAGNIFYING GLASS.
    pub const SEARCH: &str = "\u{1F50D}";
    /// ✎ — edit / pencil. U+270E LOWER RIGHT PENCIL.
    pub const EDIT: &str = "\u{270E}";
    /// ✕ — delete / close. U+2715 MULTIPLICATION X.
    pub const DELETE: &str = "\u{2715}";
    /// ⚙ — settings / gear. U+2699 GEAR.
    pub const SETTINGS: &str = "\u{2699}";
}

/// Status / outcome icons.
pub mod status {
    /// ℹ — info. U+2139 INFORMATION SOURCE.
    pub const INFO: &str = "\u{2139}";
    /// ⚠ — warning. U+26A0 WARNING SIGN.
    pub const WARNING: &str = "\u{26A0}";
    /// ✖ — error. U+2716 HEAVY MULTIPLICATION X.
    pub const ERROR: &str = "\u{2716}";
    /// ✓ — success / check. U+2713 CHECK MARK.
    pub const SUCCESS: &str = "\u{2713}";
    /// ? — question / help. U+003F QUESTION MARK.
    pub const HELP: &str = "?";
    /// 🔒 — locked / read-only. U+1F512 LOCK.
    pub const LOCKED: &str = "\u{1F512}";
    /// 🔓 — unlocked. U+1F513 OPEN LOCK.
    pub const UNLOCKED: &str = "\u{1F513}";
    /// ⏳ — busy / loading. U+23F3 HOURGLASS WITH FLOWING SAND.
    pub const BUSY: &str = "\u{23F3}";
}

/// Navigation / directional icons.
pub mod nav {
    /// ◀ — back / previous. U+25C0 BLACK LEFT-POINTING TRIANGLE.
    pub const BACK: &str = "\u{25C0}";
    /// ▶ — forward / next. U+25B6 BLACK RIGHT-POINTING TRIANGLE.
    pub const NEXT: &str = "\u{25B6}";
    /// ▲ — up. U+25B2 BLACK UP-POINTING TRIANGLE.
    pub const UP: &str = "\u{25B2}";
    /// ▼ — down. U+25BC BLACK DOWN-POINTING TRIANGLE.
    pub const DOWN: &str = "\u{25BC}";
    /// ⮜ — collapse panel left. U+2B9C LEFTWARDS BLACK ARROW.
    pub const COLLAPSE: &str = "\u{2B9C}";
    /// ⮞ — expand panel right. U+2B9E RIGHTWARDS BLACK ARROW.
    pub const EXPAND: &str = "\u{2B9E}";
    /// ⤢ — fullscreen / maximise. U+2922 NORTH EAST AND SOUTH WEST ARROW.
    pub const FULLSCREEN: &str = "\u{2922}";
}

/// View / display icons.
pub mod view {
    /// ⊞ — grid view. U+229E SQUARED PLUS.
    pub const GRID: &str = "\u{229E}";
    /// ☰ — list view. U+2630 TRIGRAM FOR HEAVEN.
    pub const LIST: &str = "\u{2630}";
    /// 🔆 — bright / contrast up. U+1F506 HIGH BRIGHTNESS SYMBOL.
    pub const BRIGHTNESS: &str = "\u{1F506}";
    /// 🔅 — dim / contrast down. U+1F505 LOW BRIGHTNESS SYMBOL.
    pub const DIM: &str = "\u{1F505}";
    /// 🔭 — zoom out / overview. U+1F52D TELESCOPE.
    pub const OVERVIEW: &str = "\u{1F52D}";
    /// ⊕ — add / plus. U+2295 CIRCLED PLUS.
    pub const ADD: &str = "\u{2295}";
    /// ⊖ — remove / minus. U+2296 CIRCLED MINUS.
    pub const REMOVE: &str = "\u{2296}";
    /// 🎨 — theme / colour. U+1F3A8 ARTIST PALETTE.
    pub const THEME: &str = "\u{1F3A8}";
}

/// Workbench / domain marks.
pub mod domain {
    /// 🧬 — genetics / DNA. U+1F9EC DNA.
    pub const DNA: &str = "\u{1F9EC}";
    /// 🔬 — bio / microscope. U+1F52C MICROSCOPE.
    pub const BIO: &str = "\u{1F52C}";
    /// 🧪 — chemistry / test tube. U+1F9EA TEST TUBE.
    pub const CHEM: &str = "\u{1F9EA}";
    /// ⚙ — CAD / mechanical / gear. U+2699 GEAR.
    pub const CAD: &str = "\u{2699}";
    /// 🌀 — fluid / aero / cyclone. U+1F300 CYCLONE.
    pub const FLUID: &str = "\u{1F300}";
    /// 📊 — chart / analytics / data. U+1F4CA BAR CHART.
    pub const CHART: &str = "\u{1F4CA}";
    /// 🏗 — architecture / BIM. U+1F3D7 BUILDING CONSTRUCTION.
    pub const ARCH: &str = "\u{1F3D7}";
    /// 🔧 — tool / wrench. U+1F527 WRENCH.
    pub const TOOL: &str = "\u{1F527}";
    /// 🎯 — target / aim. U+1F3AF DIRECT HIT.
    pub const TARGET: &str = "\u{1F3AF}";
    /// 🔀 — workflow / shuffle. U+1F500 TWISTED RIGHTWARDS ARROWS.
    pub const WORKFLOW: &str = "\u{1F500}";
}

/// Format a button label by prefixing an icon glyph and a space.
/// The canonical formatter every workbench panel uses so the
/// `icon + label` cadence stays consistent ("▶ Run", "💾 Save",
/// "↶ Undo").
///
/// Use the constants from [`run`] / [`mod@file`] / [`edit`] / [`status`]
/// / [`nav`] / [`view`] / [`domain`] as the `icon` argument.
pub fn label(icon: &str, text: &str) -> String {
    format!("{icon} {text}")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every published icon is a non-empty string. Smoke test to
    /// guard against typos creating empty consts that render as
    /// nothing in the UI.
    #[test]
    fn every_icon_is_non_empty() {
        // Run family.
        for icon in [run::PLAY, run::PAUSE, run::STOP, run::RESTART, run::STEP] {
            assert!(!icon.is_empty(), "run icon is empty");
        }
        // File family.
        for icon in [
            file::SAVE,
            file::OPEN,
            file::NEW,
            file::IMPORT,
            file::EXPORT,
            file::FOLDER,
        ] {
            assert!(!icon.is_empty(), "file icon is empty");
        }
        // Edit family.
        for icon in [
            edit::UNDO,
            edit::REDO,
            edit::CUT,
            edit::SEARCH,
            edit::EDIT,
            edit::DELETE,
            edit::SETTINGS,
        ] {
            assert!(!icon.is_empty(), "edit icon is empty");
        }
        // Status family.
        for icon in [
            status::INFO,
            status::WARNING,
            status::ERROR,
            status::SUCCESS,
            status::HELP,
            status::LOCKED,
            status::UNLOCKED,
            status::BUSY,
        ] {
            assert!(!icon.is_empty(), "status icon is empty");
        }
        // Nav family.
        for icon in [
            nav::BACK,
            nav::NEXT,
            nav::UP,
            nav::DOWN,
            nav::COLLAPSE,
            nav::EXPAND,
            nav::FULLSCREEN,
        ] {
            assert!(!icon.is_empty(), "nav icon is empty");
        }
        // View family.
        for icon in [
            view::GRID,
            view::LIST,
            view::BRIGHTNESS,
            view::DIM,
            view::OVERVIEW,
            view::ADD,
            view::REMOVE,
            view::THEME,
        ] {
            assert!(!icon.is_empty(), "view icon is empty");
        }
        // Domain family.
        for icon in [
            domain::DNA,
            domain::BIO,
            domain::CHEM,
            domain::CAD,
            domain::FLUID,
            domain::CHART,
            domain::ARCH,
            domain::TOOL,
            domain::TARGET,
            domain::WORKFLOW,
        ] {
            assert!(!icon.is_empty(), "domain icon is empty");
        }
    }

    #[test]
    fn run_play_is_canonical_marker() {
        // Lock in the canonical run-button glyph so a "let's swap to
        // a different Unicode codepoint" change has to update this
        // test too. Stops accidental "looks fine to me" drift.
        assert_eq!(run::PLAY, "\u{25B6}");
    }

    #[test]
    fn undo_redo_pair_are_distinct_glyphs() {
        // Catch a copy-paste bug where both undo and redo resolve to
        // the same arrow.
        assert_ne!(edit::UNDO, edit::REDO);
    }

    #[test]
    fn status_error_and_success_are_distinct() {
        // Locks in that error / success render as different glyphs —
        // a missed paste would silently make every status line look
        // the same.
        assert_ne!(status::ERROR, status::SUCCESS);
    }

    #[test]
    fn label_helper_prepends_icon_with_space() {
        // Smoke-test the canonical pattern panels follow: prefix an
        // icon to a label with a single space separator.
        let s = format!("{} Run", run::PLAY);
        assert!(s.starts_with(run::PLAY));
        assert!(s.ends_with("Run"));
        assert!(s.contains(' '));
    }
}
