//! NPS → outside-diameter table + schedule wall-thickness factors.

/// Look up the outside diameter (in **inches**) for a given Nominal
/// Pipe Size designation. Returns `None` for unknown sizes.
///
/// NPS is *not* equal to the bore — for 1/8" through 12" they're
/// distinct, beyond that NPS = OD in inches (ASME B36.10).
pub fn nominal_to_od_in(nps: &str) -> Option<f64> {
    match nps {
        "1/8" => Some(0.405),
        "1/4" => Some(0.540),
        "3/8" => Some(0.675),
        "1/2" => Some(0.840),
        "3/4" => Some(1.050),
        "1" => Some(1.315),
        "1-1/4" => Some(1.660),
        "1-1/2" => Some(1.900),
        "2" => Some(2.375),
        "2-1/2" => Some(2.875),
        "3" => Some(3.500),
        "3-1/2" => Some(4.000),
        "4" => Some(4.500),
        "5" => Some(5.563),
        "6" => Some(6.625),
        "8" => Some(8.625),
        "10" => Some(10.750),
        "12" => Some(12.750),
        "14" => Some(14.000),
        "16" => Some(16.000),
        "18" => Some(18.000),
        "20" => Some(20.000),
        "24" => Some(24.000),
        "30" => Some(30.000),
        "36" => Some(36.000),
        _ => None,
    }
}

/// Look up the outside diameter in **millimetres**.
pub fn nominal_to_od_mm(nps: &str) -> Option<f64> {
    nominal_to_od_in(nps).map(|inches| inches * 25.4)
}

/// Wall-thickness in **inches** for an NPS + Schedule. Returns `None`
/// when the pair isn't in the abbreviated table. v1 ships Sch 40/80/
/// 160 for the common sizes; downstream callers can extend via the
/// `Schedule::Custom(t_in)` variant.
pub fn wall_thickness_in(nps: &str, sched: Schedule) -> Option<f64> {
    if let Schedule::Custom(t) = sched {
        return Some(t);
    }
    let rows: &[(&str, f64, f64, f64)] = &[
        // (nps, sch40, sch80, sch160) inches
        ("1/8", 0.068, 0.095, 0.0),
        ("1/4", 0.088, 0.119, 0.0),
        ("3/8", 0.091, 0.126, 0.0),
        ("1/2", 0.109, 0.147, 0.187),
        ("3/4", 0.113, 0.154, 0.218),
        ("1", 0.133, 0.179, 0.250),
        ("1-1/4", 0.140, 0.191, 0.250),
        ("1-1/2", 0.145, 0.200, 0.281),
        ("2", 0.154, 0.218, 0.343),
        ("2-1/2", 0.203, 0.276, 0.375),
        ("3", 0.216, 0.300, 0.437),
        ("4", 0.237, 0.337, 0.531),
        ("6", 0.280, 0.432, 0.718),
        ("8", 0.322, 0.500, 0.906),
        ("10", 0.365, 0.593, 1.125),
        ("12", 0.375, 0.687, 1.312),
    ];
    let row = rows.iter().find(|(n, _, _, _)| *n == nps)?;
    let v = match sched {
        Schedule::Sch40 => row.1,
        Schedule::Sch80 => row.2,
        Schedule::Sch160 => row.3,
        Schedule::Custom(_) => unreachable!(),
    };
    if v > 0.0 {
        Some(v)
    } else {
        None
    }
}

/// Pipe schedule — wall thickness selector.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Schedule {
    /// Standard wall (Sch 40).
    Sch40,
    /// Extra-strong (Sch 80).
    Sch80,
    /// Double extra-strong (Sch 160).
    Sch160,
    /// Custom wall thickness in inches.
    Custom(f64),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nominal_to_od_matches_canonical_sizes() {
        assert_eq!(nominal_to_od_in("1/2"), Some(0.840));
        assert_eq!(nominal_to_od_in("2"), Some(2.375));
        assert_eq!(nominal_to_od_in("24"), Some(24.000));
        assert_eq!(nominal_to_od_in("99"), None);
    }

    #[test]
    fn od_mm_converts_correctly() {
        let v = nominal_to_od_mm("1").unwrap();
        assert!((v - 33.401).abs() < 0.001);
    }

    #[test]
    fn schedule_wall_thickness_picks_correct_column() {
        assert!((wall_thickness_in("2", Schedule::Sch40).unwrap() - 0.154).abs() < 1e-6);
        assert!((wall_thickness_in("2", Schedule::Sch80).unwrap() - 0.218).abs() < 1e-6);
        assert!((wall_thickness_in("2", Schedule::Sch160).unwrap() - 0.343).abs() < 1e-6);
    }

    #[test]
    fn schedule_custom_passes_through() {
        assert_eq!(
            wall_thickness_in("anything", Schedule::Custom(0.5)),
            Some(0.5)
        );
    }
}
