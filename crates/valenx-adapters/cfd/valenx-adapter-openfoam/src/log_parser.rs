//! Parser for the subset of OpenFOAM solver stdout that we care
//! about: iteration time markers and per-field initial residuals.
//!
//! OpenFOAM solver output is deliberately unstructured (it's a
//! fifty-year-old Fortran-flavoured log format), but the shapes we
//! want are consistent across every solver in the `v2206+` range:
//!
//! ```text
//! Time = 5
//!
//! smoothSolver:  Solving for Ux, Initial residual = 0.0123, Final residual = 1.23e-4, No Iterations 5
//! smoothSolver:  Solving for Uy, Initial residual = 0.0078, Final residual = 9.80e-5, No Iterations 4
//! GAMG:  Solving for p, Initial residual = 0.234, Final residual = 1.50e-3, No Iterations 12
//! smoothSolver:  Solving for k, Initial residual = 0.005, Final residual = 5.00e-5, No Iterations 4
//! smoothSolver:  Solving for omega, Initial residual = 0.004, Final residual = 4.00e-5, No Iterations 4
//! ExecutionTime = 0.12 s  ClockTime = 0 s
//! ```
//!
//! We parse the `Time = N` marker so the run loop knows which
//! iteration each residual belongs to, and extract the initial
//! residual per field because the UI's residual chart plots the
//! history of initial residuals over time.
//!
//! No regex — the format is fixed enough to parse with
//! `str::find` / `strip_prefix`, which keeps the adapter's
//! dependency list tiny.

use valenx_core::ResidualSample;

/// What kind of signal a single log line carries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LogSignal<'a> {
    /// The solver advanced to time `t` (for steady-state simpleFoam
    /// this is `Time = N` with `t = N` as a pseudo-time iteration
    /// counter; for transient pimpleFoam / icoFoam it's the real time
    /// in seconds, e.g. `Time = 0.0005`).
    ///
    /// We carry it as `f64` rather than `u64` because transient
    /// `Time = 0.005` would truncate to zero — the run loop never
    /// makes progress and the residual chart parks all samples on
    /// iteration 0. The run loop turns this back into an integer
    /// step counter for `ResidualSample::iteration`.
    Time(f64),
    /// Initial residual reported for a named field at the current
    /// iteration. The caller combines this with the most recent step
    /// counter to produce a [`ResidualSample`].
    Residual { field: &'a str, initial: f64 },
    /// A wall-time report — `ExecutionTime = X s`. Useful for ETA
    /// estimates.
    Execution { seconds: f64 },
    /// None of the above — keep on reading.
    Other,
}

/// Parse a single log line. Whitespace-insensitive on the left.
pub fn parse_line(line: &str) -> LogSignal<'_> {
    let trimmed = line.trim();

    // Time marker — "Time = 5" or "Time = 5.0000" or "Time = 0.0005".
    if let Some(rest) = trimmed.strip_prefix("Time = ") {
        if let Ok(v) = rest.trim().parse::<f64>() {
            return LogSignal::Time(v.max(0.0));
        }
    }

    // Residual line — "<solver>:  Solving for <field>, Initial residual = X, Final residual = Y, No Iterations Z".
    if let Some((field, initial)) = parse_residual(trimmed) {
        return LogSignal::Residual { field, initial };
    }

    // Execution time — "ExecutionTime = 0.12 s  ClockTime = 0 s".
    if let Some(rest) = trimmed.strip_prefix("ExecutionTime = ") {
        if let Some(space) = rest.find(' ') {
            if let Ok(secs) = rest[..space].parse::<f64>() {
                return LogSignal::Execution { seconds: secs };
            }
        }
    }

    LogSignal::Other
}

fn parse_residual(line: &str) -> Option<(&str, f64)> {
    let solving_idx = line.find("Solving for ")?;
    let after_solving = &line[solving_idx + "Solving for ".len()..];
    let comma_idx = after_solving.find(',')?;
    let field = after_solving[..comma_idx].trim();

    let initial_idx = line.find("Initial residual = ")?;
    let after_initial = &line[initial_idx + "Initial residual = ".len()..];
    let end_idx = after_initial.find(',').unwrap_or(after_initial.len());
    let value: f64 = after_initial[..end_idx].trim().parse().ok()?;

    Some((field, value))
}

/// Turn a known field name into the static-str form that
/// [`ResidualSample::field`] expects. OpenFOAM's log uses a bounded
/// set of field names; anything outside the set is dropped rather
/// than smuggled through.
pub fn intern_field(name: &str) -> Option<&'static str> {
    match name {
        "Ux" => Some("Ux"),
        "Uy" => Some("Uy"),
        "Uz" => Some("Uz"),
        "U" => Some("U"),
        "p" => Some("p"),
        "p_rgh" => Some("p_rgh"),
        "k" => Some("k"),
        "omega" => Some("omega"),
        "epsilon" => Some("epsilon"),
        "nuTilda" => Some("nuTilda"),
        "T" => Some("T"),
        "h" => Some("h"),
        _ => None,
    }
}

/// Convenience: if `signal` is a residual for a recognised field,
/// produce a ready-to-send `ResidualSample` under the given iteration.
pub fn signal_to_sample(signal: &LogSignal<'_>, iteration: u64) -> Option<ResidualSample> {
    if let LogSignal::Residual { field, initial } = signal {
        let interned = intern_field(field)?;
        Some(ResidualSample {
            iteration,
            field: interned,
            value: *initial,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_marker_integer() {
        match parse_line("Time = 5") {
            LogSignal::Time(t) => assert!((t - 5.0).abs() < 1e-12),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn time_marker_with_decimals() {
        match parse_line("   Time = 5.0000e+00   ") {
            LogSignal::Time(t) => assert!((t - 5.0).abs() < 1e-12),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn time_marker_subsecond_is_preserved_for_transient() {
        // Regression: previously `Time = 0.0005` truncated to 0,
        // pinning the progress bar at zero forever. The parser must
        // now return the real f64.
        match parse_line("Time = 0.0005") {
            LogSignal::Time(t) => assert!((t - 5e-4).abs() < 1e-12),
            other => panic!("{other:?}"),
        }
        match parse_line("Time = 0.5") {
            LogSignal::Time(t) => assert!((t - 0.5).abs() < 1e-12),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn residual_line_ux() {
        let line = "smoothSolver:  Solving for Ux, Initial residual = 0.0123, Final residual = 1.23e-4, No Iterations 5";
        match parse_line(line) {
            LogSignal::Residual { field, initial } => {
                assert_eq!(field, "Ux");
                assert!((initial - 0.0123).abs() < 1e-9);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn residual_line_gamg_pressure() {
        let line = "GAMG:  Solving for p, Initial residual = 0.234, Final residual = 1.50e-3, No Iterations 12";
        match parse_line(line) {
            LogSignal::Residual { field, initial } => {
                assert_eq!(field, "p");
                assert!((initial - 0.234).abs() < 1e-9);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn execution_time_line() {
        match parse_line("ExecutionTime = 0.12 s  ClockTime = 0 s") {
            LogSignal::Execution { seconds } => {
                assert!((seconds - 0.12).abs() < 1e-9);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn unknown_line_is_other() {
        assert!(matches!(parse_line(""), LogSignal::Other));
        assert!(matches!(parse_line("some header line"), LogSignal::Other));
    }

    #[test]
    fn intern_field_rejects_unknown() {
        assert_eq!(intern_field("Ux"), Some("Ux"));
        assert_eq!(intern_field("p"), Some("p"));
        assert_eq!(intern_field("unknown-field"), None);
    }

    #[test]
    fn signal_to_sample_roundtrip() {
        let signal = parse_line(
            "smoothSolver:  Solving for omega, Initial residual = 0.002, Final residual = 1e-6, No Iterations 3",
        );
        let sample = signal_to_sample(&signal, 42).expect("sample");
        assert_eq!(sample.iteration, 42);
        assert_eq!(sample.field, "omega");
        assert!((sample.value - 0.002).abs() < 1e-9);
    }

    #[test]
    fn residual_for_unknown_field_drops() {
        let signal = parse_line(
            "smoothSolver:  Solving for mystery, Initial residual = 0.01, Final residual = 1e-5, No Iterations 1",
        );
        // parse_line accepts "mystery" at the parse layer...
        assert!(matches!(
            signal,
            LogSignal::Residual {
                field: "mystery",
                ..
            }
        ));
        // ...but the interner strips it out before it reaches Results.
        assert!(signal_to_sample(&signal, 1).is_none());
    }
}
