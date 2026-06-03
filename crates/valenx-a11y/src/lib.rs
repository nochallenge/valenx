//! # valenx-a11y
//!
//! Accessibility helpers — WCAG 2.1 colour-contrast computation
//! and the AA / AAA gate predicates the accessibility audit lane
//! needs.
//!
//! ## Scope
//!
//! - Relative-luminance computation per WCAG 2.x (sRGB → linearised
//!   → weighted sum).
//! - Contrast ratio between two colours (1:1 to 21:1).
//! - AA / AAA gate predicates for normal text (4.5 / 7) and large
//!   text (3 / 4.5).
//! - A bulk audit helper that takes a list of fg/bg pairs and
//!   returns one `ContrastReport` per pair with pass/fail flags.
//!
//! Out of scope for v0.1.0: APCA contrast (the WCAG 3 successor),
//! reduced-motion preferences, screen-reader narration. Those land
//! once the core gates ship.
//!
//! ## Why a separate crate
//!
//! The contrast logic is generic over any RGB triple — it doesn't
//! depend on egui, eframe, or `valenx-design-tokens`. Keeping it
//! standalone means the WCAG arithmetic gets unit-tested without
//! an event loop, and any downstream tool (the design-tokens
//! build script, an external linter) can pull just this crate.

#![forbid(unsafe_code)]

use serde::{Deserialize, Serialize};

/// Sized 8-bit RGB triple — the shape every UI library converges
/// on. Alpha intentionally omitted: contrast is computed against
/// an opaque backdrop, so transparency is the caller's problem to
/// flatten.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Rgb {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Rgb {
    /// Construct an `Rgb` from raw 0–255 channel values.
    pub const fn new(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b }
    }

    /// Black — handy reference colour for contrast tests.
    pub const BLACK: Self = Self::new(0, 0, 0);
    /// White — same.
    pub const WHITE: Self = Self::new(255, 255, 255);
}

/// Compute relative luminance per
/// <https://www.w3.org/TR/WCAG21/#dfn-relative-luminance>.
///
/// Returns a value in `[0.0, 1.0]` — 0 for pure black, 1 for pure
/// white. The formula:
///
/// ```text
/// channel_lin = if c8/255 ≤ 0.04045 { (c8/255) / 12.92 }
///               else                 { ((c8/255 + 0.055) / 1.055).powf(2.4) }
/// L = 0.2126*R_lin + 0.7152*G_lin + 0.0722*B_lin
/// ```
pub fn relative_luminance(c: Rgb) -> f64 {
    fn linearise(c8: u8) -> f64 {
        let n = (c8 as f64) / 255.0;
        if n <= 0.04045 {
            n / 12.92
        } else {
            ((n + 0.055) / 1.055).powf(2.4)
        }
    }
    0.2126 * linearise(c.r) + 0.7152 * linearise(c.g) + 0.0722 * linearise(c.b)
}

/// WCAG contrast ratio between two colours. Always `≥ 1.0` — the
/// brighter colour is the numerator regardless of argument order.
///
/// Anchors:
/// - Black on white  → 21:1
/// - Black on black  →  1:1
/// - #767676 on white → ≈ 4.54:1 (the canonical AA threshold for
///   normal text)
pub fn contrast_ratio(a: Rgb, b: Rgb) -> f64 {
    let la = relative_luminance(a);
    let lb = relative_luminance(b);
    let (lighter, darker) = if la >= lb { (la, lb) } else { (lb, la) };
    (lighter + 0.05) / (darker + 0.05)
}

/// WCAG 2.1 conformance level for a fg/bg pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum WcagLevel {
    /// Insufficient contrast — fails AA for normal text.
    Fail,
    /// Passes AA for normal text (≥ 4.5:1) but not AAA.
    Aa,
    /// Passes AAA for normal text (≥ 7:1).
    Aaa,
}

/// Conformance for "large text" — bigger gates because larger
/// glyphs tolerate weaker contrast.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum WcagLevelLarge {
    /// Insufficient contrast — fails AA for large text.
    Fail,
    /// Passes AA for large text (≥ 3:1).
    Aa,
    /// Passes AAA for large text (≥ 4.5:1).
    Aaa,
}

/// Threshold for normal-text AA contrast.
pub const AA_NORMAL: f64 = 4.5;
/// Threshold for normal-text AAA contrast.
pub const AAA_NORMAL: f64 = 7.0;
/// Threshold for large-text AA contrast.
pub const AA_LARGE: f64 = 3.0;
/// Threshold for large-text AAA contrast.
pub const AAA_LARGE: f64 = 4.5;

/// Classify a contrast ratio against the normal-text gates.
pub fn classify_normal(ratio: f64) -> WcagLevel {
    if ratio >= AAA_NORMAL {
        WcagLevel::Aaa
    } else if ratio >= AA_NORMAL {
        WcagLevel::Aa
    } else {
        WcagLevel::Fail
    }
}

/// Classify a contrast ratio against the large-text gates.
pub fn classify_large(ratio: f64) -> WcagLevelLarge {
    if ratio >= AAA_LARGE {
        WcagLevelLarge::Aaa
    } else if ratio >= AA_LARGE {
        WcagLevelLarge::Aa
    } else {
        WcagLevelLarge::Fail
    }
}

/// Audit row for one fg/bg pair.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ContrastReport {
    /// Caller-provided label so the report is self-describing.
    pub label: String,
    pub fg: Rgb,
    pub bg: Rgb,
    pub ratio: f64,
    pub normal: WcagLevel,
    pub large: WcagLevelLarge,
}

impl ContrastReport {
    /// Build a report for the `(fg, bg)` pair labelled `label`. Computes
    /// the WCAG 2.1 contrast ratio and classifies it for both normal-text
    /// and large-text thresholds.
    pub fn for_pair(label: impl Into<String>, fg: Rgb, bg: Rgb) -> Self {
        let ratio = contrast_ratio(fg, bg);
        Self {
            label: label.into(),
            fg,
            bg,
            ratio,
            normal: classify_normal(ratio),
            large: classify_large(ratio),
        }
    }

    /// `true` if the pair passes AA for normal text — the gate the
    /// audit fails on.
    pub fn passes_aa_normal(&self) -> bool {
        !matches!(self.normal, WcagLevel::Fail)
    }
}

/// Audit a slice of `(label, fg, bg)` triples, producing one
/// [`ContrastReport`] per row. Caller decides how to surface the
/// results — log them, render a UI table, or fail a CI gate when
/// any row's `passes_aa_normal()` is `false`.
pub fn audit(pairs: &[(String, Rgb, Rgb)]) -> Vec<ContrastReport> {
    pairs
        .iter()
        .map(|(label, fg, bg)| ContrastReport::for_pair(label, *fg, *bg))
        .collect()
}

/// Subset of [`audit`] that returns only the failing rows. Useful
/// for CI gates that want a non-empty failure list to fail on.
pub fn audit_failures(pairs: &[(String, Rgb, Rgb)]) -> Vec<ContrastReport> {
    audit(pairs)
        .into_iter()
        .filter(|r| !r.passes_aa_normal())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn black_on_white_is_21_to_1() {
        let r = contrast_ratio(Rgb::BLACK, Rgb::WHITE);
        assert!(approx_eq(r, 21.0, 0.001), "got: {r}");
    }

    #[test]
    fn contrast_ratio_is_symmetric() {
        let a = contrast_ratio(Rgb::BLACK, Rgb::WHITE);
        let b = contrast_ratio(Rgb::WHITE, Rgb::BLACK);
        assert!(approx_eq(a, b, 1e-9));
    }

    #[test]
    fn identical_colours_have_one_to_one() {
        let r = contrast_ratio(Rgb::new(50, 50, 50), Rgb::new(50, 50, 50));
        assert!(approx_eq(r, 1.0, 1e-9));
    }

    #[test]
    fn middle_grey_anchor_passes_aa_normal() {
        // #767676 on white is the canonical near-AA-threshold case.
        // The actual ratio is ≈ 4.54 — just over the 4.5 gate.
        let mid = Rgb::new(0x76, 0x76, 0x76);
        let r = contrast_ratio(mid, Rgb::WHITE);
        assert!(r >= AA_NORMAL, "got: {r}");
        assert!(r < AAA_NORMAL, "got: {r}");
        assert_eq!(classify_normal(r), WcagLevel::Aa);
    }

    #[test]
    fn classify_normal_anchors() {
        // Just below AA → Fail.
        assert_eq!(classify_normal(4.49), WcagLevel::Fail);
        // Exact AA → Aa.
        assert_eq!(classify_normal(4.50), WcagLevel::Aa);
        // Just below AAA → Aa.
        assert_eq!(classify_normal(6.99), WcagLevel::Aa);
        // Exact AAA → Aaa.
        assert_eq!(classify_normal(7.00), WcagLevel::Aaa);
        // Far above AAA → Aaa.
        assert_eq!(classify_normal(15.0), WcagLevel::Aaa);
    }

    #[test]
    fn classify_large_anchors() {
        assert_eq!(classify_large(2.99), WcagLevelLarge::Fail);
        assert_eq!(classify_large(3.00), WcagLevelLarge::Aa);
        assert_eq!(classify_large(4.49), WcagLevelLarge::Aa);
        assert_eq!(classify_large(4.50), WcagLevelLarge::Aaa);
    }

    #[test]
    fn relative_luminance_anchors() {
        assert!(approx_eq(relative_luminance(Rgb::BLACK), 0.0, 1e-12));
        assert!(approx_eq(relative_luminance(Rgb::WHITE), 1.0, 1e-12));
        // Mid-grey #777777: linear value ≈ 0.182.
        let lum = relative_luminance(Rgb::new(0x77, 0x77, 0x77));
        assert!(lum > 0.17 && lum < 0.20, "got: {lum}");
    }

    #[test]
    fn contrast_report_round_trips_through_for_pair() {
        let r = ContrastReport::for_pair(
            "primary text on s1",
            Rgb::new(0xe0, 0xe0, 0xe0),
            Rgb::new(0x16, 0x16, 0x16),
        );
        assert!(r.passes_aa_normal());
        assert!(matches!(r.normal, WcagLevel::Aa | WcagLevel::Aaa));
        assert!(r.ratio > AA_NORMAL);
    }

    #[test]
    fn audit_failures_isolates_failing_rows() {
        let pairs = vec![
            ("good".to_string(), Rgb::BLACK, Rgb::WHITE),
            (
                "bad".to_string(),
                Rgb::new(200, 200, 200),
                Rgb::new(220, 220, 220),
            ),
        ];
        let all = audit(&pairs);
        assert_eq!(all.len(), 2);
        let bad = audit_failures(&pairs);
        assert_eq!(bad.len(), 1);
        assert_eq!(bad[0].label, "bad");
    }

    #[test]
    fn poor_contrast_pair_classifies_as_fail() {
        let r = contrast_ratio(Rgb::new(220, 220, 220), Rgb::WHITE);
        assert!(r < AA_NORMAL);
        assert_eq!(classify_normal(r), WcagLevel::Fail);
    }

    #[test]
    fn rgb_round_trips_through_serde() {
        let c = Rgb::new(0x12, 0x34, 0x56);
        let s = serde_json::to_string(&c).unwrap();
        let back: Rgb = serde_json::from_str(&s).unwrap();
        assert_eq!(back, c);
    }
}
