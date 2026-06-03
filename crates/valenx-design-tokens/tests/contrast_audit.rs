//! WCAG 2.1 AA contrast audit on the canonical foreground /
//! background pairs the GUI uses. Drift in `tokens.json` that
//! drops any pair below the AA threshold fails this test, so a
//! visual-regression slip is caught before it ships.
//!
//! ## Pairs covered
//!
//! - **Primary text on every surface tier (S0..S5)** —
//!   text::T1 must read everywhere the user sees the body
//!   typography, including the deepest surface (S0) and the
//!   shallowest hover-state (S5).
//! - **Secondary text on the deepest surface** — text::T2 on
//!   surface::S0 is the canonical "muted caption" pairing; AA-
//!   normal at minimum.
//! - **Accent + physics colours on surface::S1** — the surface
//!   tier most accent badges sit on. We require AA-large
//!   (≥ 3:1) since these typically appear on chips / pills /
//!   icons rather than body text.
//!
//! ## Out of scope
//!
//! - **text::T3** ("dimmed" / disabled) is intentionally not
//!   gated against AA-normal — it's a dimmed-on-purpose tier
//!   and only needs to clear AA-large.
//! - **Cross-tier pairings** (e.g. T1 on accent::PRIMARY) ride
//!   along with the future button-spec audit.
//!
//! Adding a new pair: append to the `PAIRS` slice. Changing a
//! token: tweak `tokens.json` and observe whether CI flags the
//! diff.

use valenx_a11y::{audit_failures, classify_large, contrast_ratio, Rgb, AA_LARGE};
use valenx_design_tokens::color::{accent, physics, surface, text};

fn rgb_from_u32(n: u32) -> Rgb {
    Rgb::new(
        ((n >> 16) & 0xFF) as u8,
        ((n >> 8) & 0xFF) as u8,
        (n & 0xFF) as u8,
    )
}

#[test]
fn primary_text_passes_aa_normal_on_every_surface_tier() {
    let t1 = rgb_from_u32(text::T1);
    let surfaces = [
        ("S0", surface::S0),
        ("S1", surface::S1),
        ("S2", surface::S2),
        ("S3", surface::S3),
        ("S4", surface::S4),
        ("S5", surface::S5),
    ];
    let pairs: Vec<(String, Rgb, Rgb)> = surfaces
        .iter()
        .map(|(name, hex)| {
            (
                format!("text::T1 on surface::{name}"),
                t1,
                rgb_from_u32(*hex),
            )
        })
        .collect();
    let failures = audit_failures(&pairs);
    assert!(
        failures.is_empty(),
        "WCAG AA-normal failures (need ≥ 4.5:1):\n{}",
        failures
            .iter()
            .map(|r| format!("  · {} = {:.2}:1", r.label, r.ratio))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn secondary_text_passes_aa_normal_on_deepest_surface() {
    let t2 = rgb_from_u32(text::T2);
    let s0 = rgb_from_u32(surface::S0);
    let pairs = vec![("text::T2 on surface::S0".to_string(), t2, s0)];
    let failures = audit_failures(&pairs);
    assert!(
        failures.is_empty(),
        "secondary text fails AA-normal on the deepest surface — \
         either lighten T2 or darken S0 in tokens.json"
    );
}

#[test]
fn accent_colours_pass_aa_large_on_surface_s1() {
    // Accent badges typically render at ≥ 18 pt or bold weight —
    // the AA-large gate (3:1) is the operative threshold.
    let s1 = rgb_from_u32(surface::S1);
    let accents = [
        ("accent::PRIMARY", accent::PRIMARY),
        ("accent::SUCCESS", accent::SUCCESS),
        ("accent::WARNING", accent::WARNING),
        ("accent::ERROR", accent::ERROR),
        ("accent::INFO", accent::INFO),
    ];
    let mut failures: Vec<String> = Vec::new();
    for (name, hex) in accents {
        let ratio = contrast_ratio(rgb_from_u32(hex), s1);
        if ratio < AA_LARGE {
            failures.push(format!("{name} on S1 = {ratio:.2}:1 (need ≥ {AA_LARGE})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn physics_badge_colours_pass_aa_large_on_surface_s1() {
    // Physics-domain badges follow the same gate as accents —
    // they appear on the same chip surfaces.
    let s1 = rgb_from_u32(surface::S1);
    let physics_colours = [
        ("physics::CFD", physics::CFD),
        ("physics::FEA", physics::FEA),
        ("physics::EM", physics::EM),
        ("physics::CHEM", physics::CHEM),
        ("physics::MD", physics::MD),
        ("physics::BATTERY", physics::BATTERY),
    ];
    let mut failures: Vec<String> = Vec::new();
    for (name, hex) in physics_colours {
        let ratio = contrast_ratio(rgb_from_u32(hex), s1);
        let level = classify_large(ratio);
        if matches!(level, valenx_a11y::WcagLevelLarge::Fail) {
            failures.push(format!("{name} on S1 = {ratio:.2}:1 (Fail)"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn rgb_from_u32_round_trip() {
    // Sanity-check the helper: the build script's hex-string →
    // u32 path produces the same byte layout this test consumes.
    assert_eq!(rgb_from_u32(0xFF0000), Rgb::new(0xFF, 0x00, 0x00));
    assert_eq!(rgb_from_u32(0x00FF00), Rgb::new(0x00, 0xFF, 0x00));
    assert_eq!(rgb_from_u32(0x0000FF), Rgb::new(0x00, 0x00, 0xFF));
    assert_eq!(rgb_from_u32(0x123456), Rgb::new(0x12, 0x34, 0x56));
}

/// High-contrast palette — every text tier must clear AAA (≥ 7:1) on
/// every HC surface tier. This is the WCAG-AAA gold standard for
/// normal text and is the entire point of the high-contrast variant.
#[test]
fn high_contrast_text_passes_aaa_on_every_hc_surface() {
    use valenx_a11y::AAA_NORMAL;
    use valenx_design_tokens::color::{hc_surface, hc_text};

    let hc_t1 = rgb_from_u32(hc_text::T1);
    let hc_t2 = rgb_from_u32(hc_text::T2);
    let hc_surfaces = [
        ("S0", hc_surface::S0),
        ("S1", hc_surface::S1),
        ("S2", hc_surface::S2),
        ("S3", hc_surface::S3),
        ("S4", hc_surface::S4),
        ("S5", hc_surface::S5),
    ];

    let mut failures: Vec<String> = Vec::new();
    for tier in &[("T1", hc_t1), ("T2", hc_t2)] {
        for (s_name, s_hex) in &hc_surfaces {
            let ratio = contrast_ratio(tier.1, rgb_from_u32(*s_hex));
            if ratio < AAA_NORMAL {
                failures.push(format!(
                    "hc_text::{} on hc_surface::{} = {:.2}:1 (need ≥ {AAA_NORMAL} for AAA)",
                    tier.0, s_name, ratio
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// High-contrast accent palette — every accent must clear AA (≥ 4.5)
/// against the canonical HC surface S1 since they're often used as
/// inline text-colour (status messages), not just chips.
#[test]
fn high_contrast_accents_pass_aa_normal_on_hc_surface_s1() {
    use valenx_a11y::AA_NORMAL;
    use valenx_design_tokens::color::{hc_accent, hc_surface};

    let s1 = rgb_from_u32(hc_surface::S1);
    let accents = [
        ("hc_accent::PRIMARY", hc_accent::PRIMARY),
        ("hc_accent::SUCCESS", hc_accent::SUCCESS),
        ("hc_accent::WARNING", hc_accent::WARNING),
        ("hc_accent::ERROR", hc_accent::ERROR),
        ("hc_accent::INFO", hc_accent::INFO),
    ];
    let mut failures: Vec<String> = Vec::new();
    for (name, hex) in accents {
        let ratio = contrast_ratio(rgb_from_u32(hex), s1);
        if ratio < AA_NORMAL {
            failures.push(format!("{name} on hc_surface::S1 = {ratio:.2}:1"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Light-theme palette — primary text must pass AA on every light
/// surface tier (so the rendered Light variant doesn't drop below
/// WCAG-AA contrast).
#[test]
fn light_text_passes_aa_normal_on_every_light_surface() {
    use valenx_design_tokens::color::{light_surface, light_text};

    let t1 = rgb_from_u32(light_text::T1);
    let surfaces = [
        ("S0", light_surface::S0),
        ("S1", light_surface::S1),
        ("S2", light_surface::S2),
        ("S3", light_surface::S3),
        ("S4", light_surface::S4),
        ("S5", light_surface::S5),
    ];
    let pairs: Vec<(String, Rgb, Rgb)> = surfaces
        .iter()
        .map(|(name, hex)| {
            (
                format!("light_text::T1 on light_surface::{name}"),
                t1,
                rgb_from_u32(*hex),
            )
        })
        .collect();
    let failures = audit_failures(&pairs);
    assert!(
        failures.is_empty(),
        "{:#?}",
        failures
            .iter()
            .map(|r| format!("{} = {:.2}:1", r.label, r.ratio))
            .collect::<Vec<_>>()
    );
}
