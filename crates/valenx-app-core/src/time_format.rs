//! UI rendering of [`valenx_fields::TimeKey`] values. Used by the
//! time-series slider in the Results pane to label snapshots in a
//! way the user can read at a glance.

/// Render a [`valenx_fields::TimeKey`] as a short, readable label
/// for the time-series slider. Steady runs read as `"steady"`,
/// iteration-indexed as `"iter 500"`, time-stamped as `"t=0.005 s"`.
pub fn format_time_key(time: valenx_fields::TimeKey) -> String {
    use valenx_fields::TimeKey;
    match time {
        TimeKey::Steady => "steady".to_string(),
        TimeKey::Iteration(n) => format!("iter {n}"),
        TimeKey::Time { value, units } => {
            let suffix = units.display.unwrap_or("");
            if suffix.is_empty() {
                format!("t={value:.4}")
            } else {
                format!("t={value:.4} {suffix}")
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_time_key_renders_each_variant() {
        use valenx_fields::TimeKey;
        assert_eq!(format_time_key(TimeKey::Steady), "steady");
        assert_eq!(format_time_key(TimeKey::Iteration(500)), "iter 500");
        // Time-stamped: with units suffix.
        let with_units = TimeKey::Time {
            value: 0.005,
            units: valenx_fields::units::SECOND,
        };
        let s = format_time_key(with_units);
        assert!(s.starts_with("t=0.0050"));
        assert!(s.contains("s"), "expected SI second suffix in {s:?}");
    }
}
