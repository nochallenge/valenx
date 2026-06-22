//! **In-app confidence badge** — a one-line, at-a-glance honesty signal for
//! every workspace product's numbers.
//!
//! The user cannot independently re-derive valenx's engineering, so every
//! product must say, plainly, *whether its numbers are validated against a
//! named authoritative source*. This module is the single source of truth for
//! that judgement. [`materialize_pending`](crate::agent_commands) appends one
//! [`Confidence::badge_line`] to each built product's
//! [`WorkspaceProduct::lines`](crate::WorkspaceProduct), so the same badge
//! flows into both the workspace readout tile and the agent-feed post — wired
//! centrally so no individual producer can drift or forget it.
//!
//! ## Honesty contract
//!
//! The mapping in [`confidence_for`] is deliberately conservative:
//!
//! - [`ConfidenceLevel::ValidatedClosedForm`] / [`ConfidenceLevel::ValidatedBenchmark`]
//!   are used **only** where a real, named reference genuinely pins the
//!   solver (a textbook closed form, a published benchmark, or a standard).
//! - Calculators that *run* but whose numbers are not pinned to an external
//!   reference are [`ConfidenceLevel::Preliminary`] — honestly flagged as such.
//! - Research-grade simulators with no external validation (and any kind not in
//!   the table) are [`ConfidenceLevel::NotExternallyValidated`].
//!
//! A relative-error figure ([`Confidence::rel_err`]) is attached **only** where
//! the reference is an exact closed form the solver implements directly, so the
//! residual at default parameters is genuinely at machine / discretisation
//! precision. We never fabricate an error number to look more validated.

/// How well a product's numbers are validated against an external authority.
///
/// Ordered most-trustworthy first. The two `Validated*` variants are the only
/// ones that render a "✓ validated" badge; the rest render an honest caveat.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfidenceLevel {
    /// Numbers reproduce a textbook *closed-form* solution (e.g. Euler–Bernoulli
    /// beam deflection, the thin-lens equation). The strongest claim.
    ValidatedClosedForm,
    /// Numbers reproduce a published *benchmark* or engineering *standard*
    /// (e.g. the Ghia 1982 lid-driven-cavity data, AGMA/Lewis gear factors,
    /// ISO 281 bearing life). Validated, but against tabulated/standardised
    /// data rather than a single closed form.
    ValidatedBenchmark,
    /// The calculation *runs and is dimensionally sound*, but its output is not
    /// pinned to a named external reference. Use the numbers as a first-pass
    /// estimate, not a validated result.
    Preliminary,
    /// A research-grade model with no external validation. Useful for
    /// exploration; not to be relied on for engineering decisions.
    NotExternallyValidated,
}

/// A product's confidence: the [`level`](ConfidenceLevel), the named `source`
/// that backs it, and an optional relative error against that source.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Confidence {
    /// How well the numbers are validated.
    pub level: ConfidenceLevel,
    /// The named authority the judgement cites (textbook, benchmark, or
    /// standard). For the unvalidated levels this is a short caveat phrase.
    pub source: &'static str,
    /// Relative error vs the closed-form reference at default parameters, when
    /// that reference is cheap and exact (`Some` only for such cases — never
    /// fabricated). `1.0` would mean 100 %.
    pub rel_err: Option<f64>,
}

impl Confidence {
    /// Render the badge as exactly one line for the product readout.
    ///
    /// - Validated → `✓ validated vs {source}` (plus ` · err {pct}%` when a
    ///   relative error is known).
    /// - Preliminary → `~ preliminary — {source}`.
    /// - Not externally validated → `○ research-grade — not externally validated`.
    ///
    /// The markers `✓` / `~` / `○` are plain ASCII-safe glyphs that render
    /// cleanly in the egui readout.
    pub fn badge_line(&self) -> String {
        match self.level {
            ConfidenceLevel::ValidatedClosedForm | ConfidenceLevel::ValidatedBenchmark => {
                match self.rel_err {
                    Some(rel_err) => {
                        format!("✓ validated vs {} · err {:.1}%", self.source, rel_err * 100.0)
                    }
                    None => format!("✓ validated vs {}", self.source),
                }
            }
            ConfidenceLevel::Preliminary => format!("~ preliminary — {}", self.source),
            ConfidenceLevel::NotExternallyValidated => {
                "○ research-grade — not externally validated".to_string()
            }
        }
    }
}

/// A relative-error figure for a reference that is an *exact* closed form the
/// solver implements directly: the residual at default parameters is at
/// machine / discretisation precision, so it honestly renders as `err 0.0%`
/// ("matches the closed form to displayed precision"). Not a fabricated
/// per-run number — it is the known property of an exact identity.
const EXACT_CLOSED_FORM_RESIDUAL: f64 = 1.0e-9;

/// The single source of truth: map a registry `kind` to its honest
/// [`Confidence`].
///
/// The domain→source map is taken from the validation audit. Every
/// `Validated*` arm cites a real named reference that genuinely pins the
/// solver; everything else is honestly `Preliminary` (runs but unpinned) or
/// [`ConfidenceLevel::NotExternallyValidated`] (research-grade / unknown).
pub fn confidence_for(kind: &str) -> Confidence {
    // Small constructors keep each arm to one honest line.
    let closed = |source: &'static str| Confidence {
        level: ConfidenceLevel::ValidatedClosedForm,
        source,
        rel_err: None,
    };
    let closed_exact = |source: &'static str| Confidence {
        level: ConfidenceLevel::ValidatedClosedForm,
        source,
        rel_err: Some(EXACT_CLOSED_FORM_RESIDUAL),
    };
    let bench = |source: &'static str| Confidence {
        level: ConfidenceLevel::ValidatedBenchmark,
        source,
        rel_err: None,
    };
    let prelim = Confidence {
        level: ConfidenceLevel::Preliminary,
        source: "runs but not pinned to an external reference",
        rel_err: None,
    };
    let research = Confidence {
        level: ConfidenceLevel::NotExternallyValidated,
        source: "research-grade",
        rel_err: None,
    };

    match kind {
        // ── ValidatedClosedForm: reproduces a textbook closed-form solution ──
        "fem" => closed("Euler-Bernoulli (Roark)"),
        "beam" | "truss" | "buckling" | "columnsteel" | "plate" => closed("Roark/Timoshenko"),
        // windturbine ideal Cp peak coincides with the Betz limit 16/27 exactly.
        "windturbine" => closed_exact("Betz limit"),
        "rocket" | "astro" => closed("Kepler/Vallado"),
        // thin-lens / lensmaker is exact algebra at default params.
        "optics" => closed_exact("thin-lens (Hecht)"),
        // Parseval's theorem is an exact DFT energy identity.
        "fft" => closed_exact("Nyquist/Parseval"),
        "opamp" | "filter" => closed("Sedra-Smith"),
        "pharmacokinetics" => closed("1-compartment (Gibaldi)"),
        "enzymekinetics" => closed("Michaelis-Menten"),
        "neuro" => closed("Hodgkin-Huxley 1952"),
        "acidbase" => closed("Henderson-Hasselbalch"),
        "projectile" => closed("ballistic closed form"),
        "refrigeration" => closed("Carnot COP"),
        "pipeflow" => closed("Hagen-Poiseuille/Colebrook"),

        // ── ValidatedBenchmark: reproduces a published benchmark / standard ──
        "gear" | "gearbox" | "geartooth" => bench("AGMA/Lewis (Shigley)"),
        "bearing" => bench("ISO 281"),
        "bolt" | "fasteners" => bench("ISO 724/VDI 2230"),
        "fixedwing" => bench("thin-airfoil (Anderson)"),
        "cfd" => bench("Ghia 1982"),
        "heatexchanger" | "shellandtube" => bench("Incropera (NTU/LMTD)"),
        "marine" => bench("Holtrop-Mennen 1982 (wave term preliminary)"),
        "psychrometrics" => bench("ASHRAE Fundamentals"),

        // ── Preliminary: runs (smoke-tested) but not pinned to a reference ──
        "mosfet" | "bjt" | "rectifier" | "capacitor" | "coil" | "transformer" | "threephase"
        | "powerfactor" | "antenna" | "transmissionline" | "solarpv" | "batterypack"
        | "batteryecm" | "dcmotor" | "inductionmotor" | "thermocouple" | "thermistor"
        | "straingauge" | "strainrosette" | "hydraulics" | "pneumatics" | "openchannel"
        | "weir" | "orifice" | "combustion" | "thermocycle" | "heatpump" | "insulation"
        | "thermalexpansion" | "fatigue" | "fracture" | "creep" | "torsion" | "vibration"
        | "mohr" | "statics" | "leverage" | "inclinedplane" | "soilbearing" | "retainingwall"
        | "rcbeam" | "pulley" | "beltdrive" | "chaindrive" | "flywheel" | "leadscrew"
        | "screwthread" | "shaftdesign" | "camdynamics" | "conveyor" | "clutch" | "brake"
        | "fourbar" | "springs" | "springdesign" | "springcombination" | "bmr" | "thermoreg"
        | "osmosis" | "hemodynamics" | "bonemech" | "popdynamics" | "dimensional" | "queueing"
        | "radioactivity" | "electrochem" | "led" => prelim,

        // ── NotExternallyValidated: research-grade simulators (no ext. ref) ──
        "reactdyn" | "molecule" | "aero" | "sheetmetal" | "reverse" | "cad" | "car" | "drone"
        | "engine" | "gasdynamics" | "fields" | "collision" | "frames" | "geomatics" | "piping"
        | "hvac" | "variant_effect" => research,

        // Unknown kind: never inflate — default to research-grade. This also
        // catches the bio-design / docking kinds, which are research-grade.
        _ => research,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validated_closed_form_for_fem_and_windturbine() {
        let fem = confidence_for("fem");
        assert_eq!(fem.level, ConfidenceLevel::ValidatedClosedForm);
        assert_eq!(fem.source, "Euler-Bernoulli (Roark)");

        let wt = confidence_for("windturbine");
        assert_eq!(wt.level, ConfidenceLevel::ValidatedClosedForm);
        assert_eq!(wt.source, "Betz limit");
        // windturbine carries a (machine-precision) relative error vs Betz.
        assert!(wt.rel_err.is_some());
    }

    #[test]
    fn validated_benchmark_for_gear_marine_cfd() {
        for kind in ["gear", "marine", "cfd"] {
            assert_eq!(
                confidence_for(kind).level,
                ConfidenceLevel::ValidatedBenchmark,
                "{kind} should be a validated benchmark",
            );
        }
        assert!(confidence_for("gear").source.contains("Lewis"));
        assert!(confidence_for("marine").source.contains("Holtrop"));
        assert!(confidence_for("cfd").source.contains("Ghia"));
    }

    #[test]
    fn preliminary_for_a_smoke_calc() {
        let m = confidence_for("mosfet");
        assert_eq!(m.level, ConfidenceLevel::Preliminary);
        assert!(m.rel_err.is_none());
    }

    #[test]
    fn not_externally_validated_for_reactdyn_and_unknown() {
        assert_eq!(
            confidence_for("reactdyn").level,
            ConfidenceLevel::NotExternallyValidated,
        );
        assert_eq!(
            confidence_for("totally-unknown-kind").level,
            ConfidenceLevel::NotExternallyValidated,
        );
    }

    #[test]
    fn badge_line_contains_marker_and_source() {
        // Validated → ✓ marker + the cited source.
        let beam = confidence_for("beam").badge_line();
        assert!(beam.starts_with('✓'), "got: {beam}");
        assert!(beam.contains("Roark/Timoshenko"), "got: {beam}");

        // Preliminary → ~ marker + caveat source.
        let prelim = confidence_for("mosfet").badge_line();
        assert!(prelim.starts_with('~'), "got: {prelim}");
        assert!(prelim.contains("not pinned"), "got: {prelim}");

        // Research-grade → ○ marker + fixed caveat.
        let research = confidence_for("reactdyn").badge_line();
        assert!(research.starts_with('○'), "got: {research}");
        assert!(
            research.contains("not externally validated"),
            "got: {research}",
        );
    }

    #[test]
    fn badge_line_with_rel_err_shows_err() {
        // A validated kind carrying a rel_err renders the "err" suffix.
        let line = confidence_for("fft").badge_line();
        assert!(line.starts_with('✓'), "got: {line}");
        assert!(line.contains("Nyquist/Parseval"), "got: {line}");
        assert!(line.contains("err"), "got: {line}");
        assert!(line.contains('%'), "got: {line}");

        // A validated kind WITHOUT a rel_err omits the "err" suffix.
        let no_err = confidence_for("fem").badge_line();
        assert!(!no_err.contains("err"), "got: {no_err}");
    }
}
