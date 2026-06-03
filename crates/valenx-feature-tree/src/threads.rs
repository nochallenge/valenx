//! Standard thread tables for Hole features.
//!
//! Phase 13A — Task 1-5. Provides the canonical thread-spec tables used
//! by the Hole feature to attach thread metadata to drilled pockets.
//! The metadata is consumed downstream by TechDraw (callouts) and
//! exporters; v1 does NOT model the helical thread geometry itself —
//! see [`crate::feature::HoleParams`] for the limitation.
//!
//! Four standards are supported:
//!
//! - **ISO metric** (M1.6 — M100, per ISO 261 coarse-pitch series).
//! - **Unified National** (UN, ANSI/ASME B1.1 — #0 through 4 inch).
//! - **British Standard Pipe Parallel** (BSPP, BS EN ISO 228).
//! - **NPT** (American National Standard Pipe Taper Thread, ANSI/ASME B1.20.1).
//!
//! All measurements are in **millimetres** (even for the imperial
//! standards — the tables convert at construction). Diameters and
//! pitches are the canonical published values; if you need a
//! non-standard pitch use [`ThreadSpec`] directly with custom values
//! rather than extending the tables.

use serde::{Deserialize, Serialize};

/// One of the four supported standard families.
///
/// The enum is a *family tag*; the actual pitch + diameter come from
/// the per-spec table entries returned by [`iso_metric_table`] /
/// [`un_table`] / [`bspp_table`] / [`npt_table`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ThreadStandard {
    /// ISO 261 metric coarse-pitch series. Designations like `M8` or
    /// `M16x1.5` (the fine-pitch variants are not yet in the table).
    IsoMetric,
    /// ANSI/ASME B1.1 Unified National (UN, UNC, UNF — table holds
    /// the UNC coarse series).
    UnifiedNational,
    /// BS EN ISO 228 British Standard Pipe Parallel (BSPP).
    BSPP,
    /// ANSI/ASME B1.20.1 American National Standard Taper Pipe
    /// Thread (NPT).
    NPT,
}

impl ThreadStandard {
    /// Short label for UI display (TechDraw callouts, panel headers).
    pub fn short_label(self) -> &'static str {
        match self {
            ThreadStandard::IsoMetric => "ISO Metric",
            ThreadStandard::UnifiedNational => "UN",
            ThreadStandard::BSPP => "BSPP",
            ThreadStandard::NPT => "NPT",
        }
    }
}

/// Optional class designation (fit / tolerance).
///
/// Per-standard interpretation:
/// - ISO: `6g` (external) / `6H` (internal), etc.
/// - UN: `1A`/`2A`/`3A` (external), `1B`/`2B`/`3B` (internal).
/// - BSPP: `A`/`B` (looser/tighter parallel fit).
/// - NPT: not generally used — kept as a free-form label for callouts.
pub type ThreadClass = String;

/// A single thread specification entry in a table.
///
/// Each entry is the canonical published "what to call out" for one
/// nominal-diameter row. Tap-drill, major, and minor diameters are
/// derived per-standard from the pitch with the conventions baked in
/// to the table functions; see [`Self::tap_drill_diameter`] for the
/// per-standard formulas.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ThreadSpec {
    /// Which family this spec is for.
    pub standard: ThreadStandard,
    /// Display designation (e.g. `"M8"`, `"1/4-20 UNC"`, `"G 1/2"`,
    /// `"NPT 1/4"`). Used verbatim in UI labels and callouts.
    pub designation: String,
    /// Nominal diameter in mm — the "outside" diameter the standard
    /// names the spec after. For metric M8 this is 8.0; for imperial
    /// 1/4 it's 6.35.
    pub nominal_diameter: f64,
    /// Thread pitch in mm (distance between adjacent thread crests).
    /// For imperial standards, computed from TPI as 25.4 / tpi.
    pub pitch: f64,
    /// Thread depth in mm — height of the truncated triangle from
    /// minor to major. Conservative value used for tap-drill
    /// computation; per the ISO/UN standards this is `0.61343 *
    /// pitch` (0.6P for a sharp-V) — pitch * 0.61343 for a
    /// truncated-V external thread.
    pub depth: f64,
    /// Optional fit/class designator (see [`ThreadClass`]).
    pub class: Option<ThreadClass>,
}

impl ThreadSpec {
    /// Construct a spec for a custom (non-tabular) thread.
    pub fn new(
        standard: ThreadStandard,
        designation: impl Into<String>,
        nominal_diameter: f64,
        pitch: f64,
    ) -> Self {
        let depth = pitch * 0.61343;
        Self {
            standard,
            designation: designation.into(),
            nominal_diameter,
            pitch,
            depth,
            class: None,
        }
    }

    /// The diameter of the *tap-drill* — i.e. the hole you drill
    /// before tapping a female thread. Per ASME / ISO convention this
    /// is `major - 1.0825 * pitch` for a 75 %-engagement thread. For
    /// NPT (tapered) we return the nominal minor at the small end.
    pub fn tap_drill_diameter(&self) -> f64 {
        match self.standard {
            ThreadStandard::IsoMetric | ThreadStandard::UnifiedNational | ThreadStandard::BSPP => {
                // 75 % engagement formula common to ISO + UN tap drills.
                self.nominal_diameter - 1.0825 * self.pitch
            }
            ThreadStandard::NPT => {
                // NPT is tapered; the "tap drill" is normally the
                // small-end ID at full engagement. Approximate as
                // nominal - 1.5 * depth which lands close to the
                // published table values for the 1/8 — 4 inch range.
                self.nominal_diameter - 1.5 * self.depth
            }
        }
    }

    /// Major (outside) diameter of the thread. Equals
    /// `nominal_diameter` for metric / UN; equals the published OD for
    /// pipe threads (which is larger than the nominal designation —
    /// e.g. G 1/2 has OD ≈ 20.955 mm).
    pub fn major_diameter(&self) -> f64 {
        self.nominal_diameter
    }

    /// Minor (root) diameter of the thread. Equals
    /// `major - 2 * depth` for a parallel thread.
    pub fn minor_diameter(&self) -> f64 {
        self.nominal_diameter - 2.0 * self.depth
    }
}

/// ISO 261 metric coarse-pitch series, M1.6 through M100.
///
/// Pitch values are from the standard's "first choice" column for
/// general engineering. Designations are the bare `M<diameter>` form
/// — fine-pitch variants (`M8x1`, etc.) are not in this table.
pub fn iso_metric_table() -> Vec<ThreadSpec> {
    // (nominal mm, coarse pitch mm). Numbers from ISO 261 Table 1.
    let rows: &[(f64, f64)] = &[
        (1.6, 0.35),
        (2.0, 0.4),
        (2.5, 0.45),
        (3.0, 0.5),
        (4.0, 0.7),
        (5.0, 0.8),
        (6.0, 1.0),
        (8.0, 1.25),
        (10.0, 1.5),
        (12.0, 1.75),
        (14.0, 2.0),
        (16.0, 2.0),
        (18.0, 2.5),
        (20.0, 2.5),
        (22.0, 2.5),
        (24.0, 3.0),
        (27.0, 3.0),
        (30.0, 3.5),
        (33.0, 3.5),
        (36.0, 4.0),
        (39.0, 4.0),
        (42.0, 4.5),
        (45.0, 4.5),
        (48.0, 5.0),
        (52.0, 5.0),
        (56.0, 5.5),
        (60.0, 5.5),
        (64.0, 6.0),
        (68.0, 6.0),
        (72.0, 6.0),
        (76.0, 6.0),
        (80.0, 6.0),
        (85.0, 6.0),
        (90.0, 6.0),
        (95.0, 6.0),
        (100.0, 6.0),
    ];
    rows.iter()
        .map(|(d, p)| {
            // Format trailing zero for "nice" designations: 1.6 → "M1.6",
            // 6.0 → "M6", not "M6.0".
            let label = if (d - d.trunc()).abs() < 1e-9 {
                format!("M{}", *d as i64)
            } else {
                format!("M{d}")
            };
            ThreadSpec::new(ThreadStandard::IsoMetric, label, *d, *p)
        })
        .collect()
}

/// ANSI/ASME B1.1 Unified National coarse-pitch (UNC) series.
///
/// Covers numbered sizes #0 — #10 (where #N nominal = (0.060 + 0.013*N)
/// inch) and fractional sizes 1/4 — 4 inch. The "TPI" (threads-per-inch)
/// → pitch conversion is `25.4 / TPI` mm.
pub fn un_table() -> Vec<ThreadSpec> {
    // (designation, nominal inch, TPI)
    let rows: &[(&str, f64, f64)] = &[
        ("#0-80", 0.060, 80.0),
        ("#2-56", 0.086, 56.0),
        ("#4-40", 0.112, 40.0),
        ("#6-32", 0.138, 32.0),
        ("#8-32", 0.164, 32.0),
        ("#10-24", 0.190, 24.0),
        ("1/4-20", 0.250, 20.0),
        ("5/16-18", 0.3125, 18.0),
        ("3/8-16", 0.375, 16.0),
        ("7/16-14", 0.4375, 14.0),
        ("1/2-13", 0.500, 13.0),
        ("5/8-11", 0.625, 11.0),
        ("3/4-10", 0.750, 10.0),
        ("7/8-9", 0.875, 9.0),
        ("1-8", 1.000, 8.0),
        ("1-1/8-7", 1.125, 7.0),
        ("1-1/4-7", 1.250, 7.0),
        ("1-3/8-6", 1.375, 6.0),
        ("1-1/2-6", 1.500, 6.0),
        ("1-3/4-5", 1.750, 5.0),
        ("2-4-1/2", 2.000, 4.5),
        ("2-1/4-4-1/2", 2.250, 4.5),
        ("2-1/2-4", 2.500, 4.0),
        ("2-3/4-4", 2.750, 4.0),
        ("3-4", 3.000, 4.0),
        ("3-1/4-4", 3.250, 4.0),
        ("3-1/2-4", 3.500, 4.0),
        ("3-3/4-4", 3.750, 4.0),
        ("4-4", 4.000, 4.0),
    ];
    rows.iter()
        .map(|(label, in_dia, tpi)| {
            ThreadSpec::new(
                ThreadStandard::UnifiedNational,
                *label,
                *in_dia * 25.4,
                25.4 / *tpi,
            )
        })
        .collect()
}

/// BSPP (BS EN ISO 228) parallel-pipe thread table.
///
/// BSPP designations are written `G <size>`. The "size" is the
/// historic *bore* of the pipe in inches, NOT the thread OD; the table
/// captures the actual major diameter the standard publishes.
pub fn bspp_table() -> Vec<ThreadSpec> {
    // (designation, major-OD mm, TPI). TPI = 28 for 1/16-1/4, 19 for
    // 3/8-1/2, 14 for 5/8-1, 11 for >=1-1/4.
    let rows: &[(&str, f64, f64)] = &[
        ("G 1/16", 7.723, 28.0),
        ("G 1/8", 9.728, 28.0),
        ("G 1/4", 13.157, 19.0),
        ("G 3/8", 16.662, 19.0),
        ("G 1/2", 20.955, 14.0),
        ("G 5/8", 22.911, 14.0),
        ("G 3/4", 26.441, 14.0),
        ("G 7/8", 30.201, 14.0),
        ("G 1", 33.249, 11.0),
        ("G 1-1/4", 41.910, 11.0),
        ("G 1-1/2", 47.803, 11.0),
        ("G 2", 59.614, 11.0),
        ("G 2-1/2", 75.184, 11.0),
        ("G 3", 87.884, 11.0),
        ("G 4", 113.030, 11.0),
        ("G 5", 138.430, 11.0),
        ("G 6", 163.830, 11.0),
    ];
    rows.iter()
        .map(|(label, dia, tpi)| ThreadSpec::new(ThreadStandard::BSPP, *label, *dia, 25.4 / *tpi))
        .collect()
}

/// NPT (ANSI/ASME B1.20.1) tapered-pipe thread table.
///
/// Same nominal-size scheme as BSPP but with a 1° 47′ 24″ taper. We
/// surface only the nominal OD + pitch — the taper itself is metadata
/// the downstream consumer can apply if it wants to model the cone.
pub fn npt_table() -> Vec<ThreadSpec> {
    // (designation, major-OD mm at large end, TPI).
    let rows: &[(&str, f64, f64)] = &[
        ("NPT 1/16", 7.895, 27.0),
        ("NPT 1/8", 10.272, 27.0),
        ("NPT 1/4", 13.716, 18.0),
        ("NPT 3/8", 17.145, 18.0),
        ("NPT 1/2", 21.336, 14.0),
        ("NPT 3/4", 26.670, 14.0),
        ("NPT 1", 33.401, 11.5),
        ("NPT 1-1/4", 42.164, 11.5),
        ("NPT 1-1/2", 48.260, 11.5),
        ("NPT 2", 60.325, 11.5),
        ("NPT 2-1/2", 73.025, 8.0),
        ("NPT 3", 88.900, 8.0),
        ("NPT 4", 114.300, 8.0),
    ];
    rows.iter()
        .map(|(label, dia, tpi)| ThreadSpec::new(ThreadStandard::NPT, *label, *dia, 25.4 / *tpi))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iso_metric_table_includes_m8_with_pitch_1_25() {
        let table = iso_metric_table();
        let m8 = table
            .iter()
            .find(|s| s.designation == "M8")
            .expect("M8 in ISO table");
        assert!((m8.pitch - 1.25).abs() < 1e-9);
        assert!((m8.nominal_diameter - 8.0).abs() < 1e-9);
        // Tap drill for M8x1.25 ≈ 8 - 1.0825*1.25 ≈ 6.647 mm (close to
        // the standard 6.8 mm drill bit).
        assert!((m8.tap_drill_diameter() - (8.0 - 1.0825 * 1.25)).abs() < 1e-9);
    }

    #[test]
    fn un_table_quarter_twenty_has_canonical_pitch() {
        let table = un_table();
        let qt = table
            .iter()
            .find(|s| s.designation == "1/4-20")
            .expect("1/4-20");
        // Pitch = 25.4 / 20 = 1.27 mm.
        assert!((qt.pitch - 1.27).abs() < 1e-9);
        // Nominal = 0.25 inch = 6.35 mm.
        assert!((qt.nominal_diameter - 6.35).abs() < 1e-9);
    }

    #[test]
    fn bspp_g_half_uses_published_od() {
        let table = bspp_table();
        let g_half = table
            .iter()
            .find(|s| s.designation == "G 1/2")
            .expect("G 1/2");
        // Standard published OD for G 1/2 = 20.955 mm.
        assert!((g_half.nominal_diameter - 20.955).abs() < 1e-9);
        // TPI 14 → pitch 25.4/14 ≈ 1.8143 mm.
        assert!((g_half.pitch - 25.4 / 14.0).abs() < 1e-9);
    }

    #[test]
    fn npt_table_present_with_quarter_inch_entry() {
        let table = npt_table();
        assert!(table.iter().any(|s| s.designation == "NPT 1/4"));
        // NPT 1/4 has TPI = 18.
        let npt_q = table.iter().find(|s| s.designation == "NPT 1/4").unwrap();
        assert!((npt_q.pitch - 25.4 / 18.0).abs() < 1e-9);
    }

    #[test]
    fn major_and_minor_diameters_are_related_by_2_x_depth() {
        let spec = ThreadSpec::new(ThreadStandard::IsoMetric, "M10", 10.0, 1.5);
        // depth = 1.5 * 0.61343 ≈ 0.92014
        // minor = 10 - 2*0.92014 ≈ 8.16
        assert!((spec.major_diameter() - 10.0).abs() < 1e-9);
        assert!((spec.minor_diameter() - (10.0 - 2.0 * 1.5 * 0.61343)).abs() < 1e-9);
    }

    #[test]
    fn thread_standard_short_labels() {
        assert_eq!(ThreadStandard::IsoMetric.short_label(), "ISO Metric");
        assert_eq!(ThreadStandard::UnifiedNational.short_label(), "UN");
        assert_eq!(ThreadStandard::BSPP.short_label(), "BSPP");
        assert_eq!(ThreadStandard::NPT.short_label(), "NPT");
    }

    #[test]
    fn iso_table_has_expected_extent_endpoints() {
        let table = iso_metric_table();
        // First entry should be M1.6, last should be M100.
        assert_eq!(table.first().unwrap().designation, "M1.6");
        assert_eq!(table.last().unwrap().designation, "M100");
    }
}
