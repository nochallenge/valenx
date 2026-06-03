//! Thread tables — full ISO metric, UN/UNC/UNF/UNEF, BSPP/BSPT,
//! NPT/NPTF/NPS, Acme, Trapezoidal, Whitworth.

use crate::spec::ThreadSpecPro;
use crate::standard::ThreadStandardPro;

/// Full ISO 261 metric coarse-pitch series, M0.5 through M300.
///
/// Pitches follow the published ISO 261 first-choice series.
pub fn metric_full_table() -> Vec<ThreadSpecPro> {
    // (designation suffix in mm, pitch in mm). One row per nominal.
    // Source: ISO 261:1998 Table 1 first-choice column.
    let rows: &[(f64, f64)] = &[
        (0.5, 0.125),
        (0.6, 0.15),
        (0.8, 0.2),
        (1.0, 0.25),
        (1.2, 0.25),
        (1.4, 0.3),
        (1.6, 0.35),
        (1.8, 0.35),
        (2.0, 0.4),
        (2.5, 0.45),
        (3.0, 0.5),
        (3.5, 0.6),
        (4.0, 0.7),
        (5.0, 0.8),
        (6.0, 1.0),
        (7.0, 1.0),
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
        (105.0, 6.0),
        (110.0, 6.0),
        (115.0, 6.0),
        (120.0, 6.0),
        (130.0, 6.0),
        (140.0, 6.0),
        (150.0, 6.0),
        (160.0, 6.0),
        (170.0, 6.0),
        (180.0, 6.0),
        (190.0, 6.0),
        (200.0, 6.0),
        (220.0, 6.0),
        (240.0, 6.0),
        (260.0, 6.0),
        (280.0, 6.0),
        (300.0, 6.0),
    ];
    rows.iter()
        .map(|&(d, p)| {
            ThreadSpecPro::new(
                ThreadStandardPro::IsoMetric,
                format!("M{}", short_num(d)),
                d,
                p,
            )
        })
        .collect()
}

/// Metric fine pitches (subset — common engineering sizes).
pub fn metric_fine_table() -> Vec<ThreadSpecPro> {
    let rows: &[(f64, f64)] = &[
        (4.0, 0.5),
        (5.0, 0.5),
        (6.0, 0.75),
        (8.0, 1.0),
        (10.0, 1.0),
        (12.0, 1.25),
        (14.0, 1.5),
        (16.0, 1.5),
        (18.0, 1.5),
        (20.0, 1.5),
        (22.0, 1.5),
        (24.0, 2.0),
        (27.0, 2.0),
        (30.0, 2.0),
        (33.0, 2.0),
        (36.0, 3.0),
        (39.0, 3.0),
        (42.0, 3.0),
        (45.0, 3.0),
        (48.0, 3.0),
        (52.0, 4.0),
        (56.0, 4.0),
        (60.0, 4.0),
        (64.0, 4.0),
        (72.0, 4.0),
        (80.0, 4.0),
        (90.0, 4.0),
        (100.0, 4.0),
    ];
    rows.iter()
        .map(|&(d, p)| {
            ThreadSpecPro::new(
                ThreadStandardPro::MetricFine,
                format!("M{}x{}", short_num(d), short_num(p)),
                d,
                p,
            )
        })
        .collect()
}

/// Metric extra-fine.
pub fn metric_extra_fine_table() -> Vec<ThreadSpecPro> {
    let rows: &[(f64, f64)] = &[
        (8.0, 0.5),
        (10.0, 0.5),
        (10.0, 0.75),
        (12.0, 0.5),
        (12.0, 0.75),
        (12.0, 1.0),
        (14.0, 1.0),
        (16.0, 1.0),
        (18.0, 1.0),
        (20.0, 1.0),
        (22.0, 1.0),
        (24.0, 1.0),
        (24.0, 1.5),
    ];
    rows.iter()
        .map(|&(d, p)| {
            ThreadSpecPro::new(
                ThreadStandardPro::MetricExtraFine,
                format!("M{}x{}", short_num(d), short_num(p)),
                d,
                p,
            )
        })
        .collect()
}

/// Unified National Coarse (UNC) — #0 through 4 in.
pub fn unc_table() -> Vec<ThreadSpecPro> {
    // (display, nominal in mm, tpi).
    let rows: &[(&str, f64, f64)] = &[
        ("#0-80", 1.524, 80.0),
        ("#1-64", 1.854, 64.0),
        ("#2-56", 2.184, 56.0),
        ("#3-48", 2.515, 48.0),
        ("#4-40", 2.845, 40.0),
        ("#5-40", 3.175, 40.0),
        ("#6-32", 3.505, 32.0),
        ("#8-32", 4.166, 32.0),
        ("#10-24", 4.826, 24.0),
        ("#12-24", 5.486, 24.0),
        ("1/4-20", 6.35, 20.0),
        ("5/16-18", 7.938, 18.0),
        ("3/8-16", 9.525, 16.0),
        ("7/16-14", 11.113, 14.0),
        ("1/2-13", 12.7, 13.0),
        ("9/16-12", 14.288, 12.0),
        ("5/8-11", 15.875, 11.0),
        ("3/4-10", 19.05, 10.0),
        ("7/8-9", 22.225, 9.0),
        ("1-8", 25.4, 8.0),
        ("1-1/8-7", 28.575, 7.0),
        ("1-1/4-7", 31.75, 7.0),
        ("1-3/8-6", 34.925, 6.0),
        ("1-1/2-6", 38.1, 6.0),
        ("1-3/4-5", 44.45, 5.0),
        ("2-4.5", 50.8, 4.5),
        ("2-1/4-4.5", 57.15, 4.5),
        ("2-1/2-4", 63.5, 4.0),
        ("2-3/4-4", 69.85, 4.0),
        ("3-4", 76.2, 4.0),
        ("3-1/4-4", 82.55, 4.0),
        ("3-1/2-4", 88.9, 4.0),
        ("3-3/4-4", 95.25, 4.0),
        ("4-4", 101.6, 4.0),
    ];
    rows.iter()
        .map(|&(name, d, tpi)| {
            ThreadSpecPro::new(
                ThreadStandardPro::UnifiedNationalCoarse,
                format!("{name} UNC"),
                d,
                25.4 / tpi,
            )
        })
        .collect()
}

/// Unified National Fine (UNF).
pub fn unf_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("#0-80", 1.524, 80.0),
        ("#1-72", 1.854, 72.0),
        ("#2-64", 2.184, 64.0),
        ("#3-56", 2.515, 56.0),
        ("#4-48", 2.845, 48.0),
        ("#5-44", 3.175, 44.0),
        ("#6-40", 3.505, 40.0),
        ("#8-36", 4.166, 36.0),
        ("#10-32", 4.826, 32.0),
        ("#12-28", 5.486, 28.0),
        ("1/4-28", 6.35, 28.0),
        ("5/16-24", 7.938, 24.0),
        ("3/8-24", 9.525, 24.0),
        ("7/16-20", 11.113, 20.0),
        ("1/2-20", 12.7, 20.0),
        ("9/16-18", 14.288, 18.0),
        ("5/8-18", 15.875, 18.0),
        ("3/4-16", 19.05, 16.0),
        ("7/8-14", 22.225, 14.0),
        ("1-12", 25.4, 12.0),
        ("1-1/8-12", 28.575, 12.0),
        ("1-1/4-12", 31.75, 12.0),
        ("1-3/8-12", 34.925, 12.0),
        ("1-1/2-12", 38.1, 12.0),
    ];
    rows.iter()
        .map(|&(name, d, tpi)| {
            ThreadSpecPro::new(
                ThreadStandardPro::UnifiedNationalFine,
                format!("{name} UNF"),
                d,
                25.4 / tpi,
            )
        })
        .collect()
}

/// Unified National Extra-Fine (UNEF).
pub fn unef_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("#12-32", 5.486, 32.0),
        ("1/4-32", 6.35, 32.0),
        ("5/16-32", 7.938, 32.0),
        ("3/8-32", 9.525, 32.0),
        ("7/16-28", 11.113, 28.0),
        ("1/2-28", 12.7, 28.0),
        ("9/16-24", 14.288, 24.0),
        ("5/8-24", 15.875, 24.0),
        ("11/16-24", 17.463, 24.0),
        ("3/4-20", 19.05, 20.0),
        ("13/16-20", 20.638, 20.0),
        ("7/8-20", 22.225, 20.0),
        ("15/16-20", 23.813, 20.0),
        ("1-20", 25.4, 20.0),
        ("1-1/16-18", 26.988, 18.0),
        ("1-1/8-18", 28.575, 18.0),
        ("1-3/16-18", 30.163, 18.0),
        ("1-1/4-18", 31.75, 18.0),
        ("1-5/16-18", 33.338, 18.0),
        ("1-3/8-18", 34.925, 18.0),
        ("1-7/16-18", 36.513, 18.0),
        ("1-1/2-18", 38.1, 18.0),
    ];
    rows.iter()
        .map(|&(name, d, tpi)| {
            ThreadSpecPro::new(
                ThreadStandardPro::UnifiedNationalExtraFine,
                format!("{name} UNEF"),
                d,
                25.4 / tpi,
            )
        })
        .collect()
}

/// BSPP G1/16 through G6 (parallel pipe thread).
pub fn bspp_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("G 1/16", 7.723, 0.907),
        ("G 1/8", 9.728, 0.907),
        ("G 1/4", 13.157, 1.337),
        ("G 3/8", 16.662, 1.337),
        ("G 1/2", 20.955, 1.814),
        ("G 5/8", 22.911, 1.814),
        ("G 3/4", 26.441, 1.814),
        ("G 7/8", 30.201, 1.814),
        ("G 1", 33.249, 2.309),
        ("G 1-1/8", 37.897, 2.309),
        ("G 1-1/4", 41.910, 2.309),
        ("G 1-3/8", 44.323, 2.309),
        ("G 1-1/2", 47.803, 2.309),
        ("G 1-3/4", 53.746, 2.309),
        ("G 2", 59.614, 2.309),
        ("G 2-1/4", 65.710, 2.309),
        ("G 2-1/2", 75.184, 2.309),
        ("G 2-3/4", 81.534, 2.309),
        ("G 3", 87.884, 2.309),
        ("G 3-1/2", 100.330, 2.309),
        ("G 4", 113.030, 2.309),
        ("G 4-1/2", 125.730, 2.309),
        ("G 5", 138.430, 2.309),
        ("G 5-1/2", 151.130, 2.309),
        ("G 6", 163.830, 2.309),
    ];
    rows.iter()
        .map(|&(name, d, p)| ThreadSpecPro::new(ThreadStandardPro::BSPP, name, d, p))
        .collect()
}

/// BSPT R1/16 through R6 (taper pipe thread). Uses major diameter at
/// the gauging plane.
pub fn bspt_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("R 1/16", 7.723, 0.907),
        ("R 1/8", 9.728, 0.907),
        ("R 1/4", 13.157, 1.337),
        ("R 3/8", 16.662, 1.337),
        ("R 1/2", 20.955, 1.814),
        ("R 3/4", 26.441, 1.814),
        ("R 1", 33.249, 2.309),
        ("R 1-1/4", 41.910, 2.309),
        ("R 1-1/2", 47.803, 2.309),
        ("R 2", 59.614, 2.309),
        ("R 2-1/2", 75.184, 2.309),
        ("R 3", 87.884, 2.309),
        ("R 4", 113.030, 2.309),
        ("R 5", 138.430, 2.309),
        ("R 6", 163.830, 2.309),
    ];
    rows.iter()
        .map(|&(name, d, p)| ThreadSpecPro::new(ThreadStandardPro::BSPT, name, d, p))
        .collect()
}

/// NPT (American National Standard Pipe Taper).
pub fn npt_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("NPT 1/16", 7.895, 25.4 / 27.0),
        ("NPT 1/8", 10.287, 25.4 / 27.0),
        ("NPT 1/4", 13.716, 25.4 / 18.0),
        ("NPT 3/8", 17.145, 25.4 / 18.0),
        ("NPT 1/2", 21.336, 25.4 / 14.0),
        ("NPT 3/4", 26.670, 25.4 / 14.0),
        ("NPT 1", 33.401, 25.4 / 11.5),
        ("NPT 1-1/4", 42.164, 25.4 / 11.5),
        ("NPT 1-1/2", 48.260, 25.4 / 11.5),
        ("NPT 2", 60.325, 25.4 / 11.5),
        ("NPT 2-1/2", 73.025, 25.4 / 8.0),
        ("NPT 3", 88.900, 25.4 / 8.0),
        ("NPT 3-1/2", 101.600, 25.4 / 8.0),
        ("NPT 4", 114.300, 25.4 / 8.0),
    ];
    rows.iter()
        .map(|&(name, d, p)| ThreadSpecPro::new(ThreadStandardPro::NPT, name, d, p))
        .collect()
}

/// NPTF (dryseal) — identical dimensions to NPT; the difference is
/// crest/root truncation, captured in the family tag.
pub fn nptf_table() -> Vec<ThreadSpecPro> {
    npt_table()
        .into_iter()
        .map(|mut s| {
            s.standard = ThreadStandardPro::NPTF;
            s.designation = s.designation.replace("NPT", "NPTF");
            s
        })
        .collect()
}

/// NPS straight pipe thread — same diameters / pitches as NPT but
/// parallel.
pub fn nps_table() -> Vec<ThreadSpecPro> {
    npt_table()
        .into_iter()
        .map(|mut s| {
            s.standard = ThreadStandardPro::NPS;
            s.designation = s.designation.replace("NPT", "NPS");
            s
        })
        .collect()
}

/// Acme general-purpose thread (ASME B1.5 — selected sizes).
pub fn acme_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("1/4-16 Acme", 6.35, 25.4 / 16.0),
        ("5/16-14 Acme", 7.938, 25.4 / 14.0),
        ("3/8-12 Acme", 9.525, 25.4 / 12.0),
        ("7/16-12 Acme", 11.113, 25.4 / 12.0),
        ("1/2-10 Acme", 12.7, 25.4 / 10.0),
        ("5/8-8 Acme", 15.875, 25.4 / 8.0),
        ("3/4-6 Acme", 19.05, 25.4 / 6.0),
        ("7/8-6 Acme", 22.225, 25.4 / 6.0),
        ("1-5 Acme", 25.4, 25.4 / 5.0),
        ("1-1/4-5 Acme", 31.75, 25.4 / 5.0),
        ("1-1/2-4 Acme", 38.1, 25.4 / 4.0),
        ("2-4 Acme", 50.8, 25.4 / 4.0),
        ("2-1/2-3 Acme", 63.5, 25.4 / 3.0),
        ("3-2 Acme", 76.2, 25.4 / 2.0),
    ];
    rows.iter()
        .map(|&(name, d, p)| ThreadSpecPro::new(ThreadStandardPro::Acme, name, d, p))
        .collect()
}

/// ISO 2904 / DIN 103 metric trapezoidal (selected sizes).
pub fn trapezoidal_iso_table() -> Vec<ThreadSpecPro> {
    let rows: &[(f64, f64)] = &[
        (8.0, 1.5),
        (10.0, 2.0),
        (12.0, 2.0),
        (14.0, 3.0),
        (16.0, 3.0),
        (18.0, 4.0),
        (20.0, 4.0),
        (22.0, 5.0),
        (24.0, 5.0),
        (26.0, 5.0),
        (28.0, 5.0),
        (30.0, 6.0),
        (32.0, 6.0),
        (36.0, 6.0),
        (40.0, 7.0),
        (44.0, 7.0),
        (48.0, 8.0),
        (50.0, 8.0),
        (60.0, 9.0),
        (70.0, 10.0),
        (80.0, 10.0),
        (100.0, 12.0),
    ];
    rows.iter()
        .map(|&(d, p)| {
            ThreadSpecPro::new(
                ThreadStandardPro::TrapezoidalIso,
                format!("Tr {}x{}", short_num(d), short_num(p)),
                d,
                p,
            )
        })
        .collect()
}

/// BSW Whitworth (BS 84).
pub fn whitworth_bsw_table() -> Vec<ThreadSpecPro> {
    let rows: &[(&str, f64, f64)] = &[
        ("1/8 W", 3.175, 25.4 / 40.0),
        ("3/16 W", 4.763, 25.4 / 24.0),
        ("1/4 W", 6.35, 25.4 / 20.0),
        ("5/16 W", 7.938, 25.4 / 18.0),
        ("3/8 W", 9.525, 25.4 / 16.0),
        ("7/16 W", 11.113, 25.4 / 14.0),
        ("1/2 W", 12.7, 25.4 / 12.0),
        ("5/8 W", 15.875, 25.4 / 11.0),
        ("3/4 W", 19.05, 25.4 / 10.0),
        ("7/8 W", 22.225, 25.4 / 9.0),
        ("1 W", 25.4, 25.4 / 8.0),
    ];
    rows.iter()
        .map(|&(name, d, p)| ThreadSpecPro::new(ThreadStandardPro::WhitworthBSW, name, d, p))
        .collect()
}

fn short_num(v: f64) -> String {
    let i = v as i64;
    if (v - i as f64).abs() < 1e-9 {
        i.to_string()
    } else {
        // Always emit at most 3 decimals, then trim trailing zeros.
        let mut s = format!("{v:.3}");
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metric_full_table_has_many_entries() {
        let t = metric_full_table();
        assert!(t.len() >= 60);
        assert!(t.iter().any(|s| s.designation == "M8"));
        assert!(t.iter().any(|s| s.designation == "M300"));
    }

    #[test]
    fn un_tables_combined_size_meets_target() {
        let n = unc_table().len() + unf_table().len() + unef_table().len();
        assert!(n >= 75, "got {n}");
    }

    #[test]
    fn bspp_starts_at_1_16_and_ends_at_6() {
        let t = bspp_table();
        assert_eq!(t.first().unwrap().designation, "G 1/16");
        assert_eq!(t.last().unwrap().designation, "G 6");
    }

    #[test]
    fn npt_pitch_matches_tpi() {
        let t = npt_table();
        let half = t.iter().find(|s| s.designation == "NPT 1/2").unwrap();
        assert!((half.pitch - (25.4 / 14.0)).abs() < 1e-6);
    }

    #[test]
    fn acme_table_returns_acme_profile() {
        let t = acme_table();
        for s in t {
            assert_eq!(s.standard, ThreadStandardPro::Acme);
            assert_eq!(s.profile, crate::standard::ProfileShape::Acme);
        }
    }

    #[test]
    fn trapezoidal_designation_has_tr_prefix() {
        let t = trapezoidal_iso_table();
        for s in t {
            assert!(s.designation.starts_with("Tr "));
        }
    }
}
