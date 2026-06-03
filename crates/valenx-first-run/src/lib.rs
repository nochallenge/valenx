//! # valenx-first-run
//!
//! First-launch wizard state machine. Pure logic — the egui shell
//! that renders the dialog lives in `valenx-app`. Splitting the
//! decision layer from the rendering layer keeps this crate
//! unit-testable without spinning up an event loop.
//!
//! ## What the wizard does
//!
//! On a fresh install the user has no idea which solvers Valenx
//! supports or whether the ones they care about are already
//! installed on `PATH`. The wizard:
//!
//! 1. Walks the registered adapter list and probes each.
//! 2. Renders a per-adapter status row (Installed / Missing /
//!    Outdated / Broken / Disabled) with a per-OS install hint
//!    next to the missing ones.
//! 3. Lets the user check / uncheck which adapters they expect
//!    to use. The selection persists to `settings.json` so the
//!    GUI's case-browser can grey out un-selected adapters
//!    instead of nagging on every run.
//! 4. Records that the wizard ran, so it never re-opens unless
//!    the user explicitly invokes it from the command palette.
//!
//! ## Why this is its own crate
//!
//! - `valenx-app`'s test pipeline currently writes a `settings.json`
//!   on every test run (see the JSON-file constraint in CONTRIBUTING).
//!   Keeping the wizard logic out of that crate means we can test
//!   it freely.
//! - The wizard touches the project loader (case names) and the
//!   adapter registry (probe results). Both already live downstream
//!   of `valenx-app` in the dep graph; pulling the logic up keeps
//!   the dep direction acyclic.

#![forbid(unsafe_code)]

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// Per-adapter probe outcome the wizard renders.
///
/// Mirrors the canonical `valenx-core::AdapterStatus` enum but lives
/// in this crate so we don't drag the full adapter dep tree in. The
/// caller (typically `valenx-app::init_registry`) maps the canonical
/// status to one of these variants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum AdapterAvailability {
    /// Probe succeeded — tool is on PATH and the version matches.
    Installed,
    /// Tool not on PATH.
    Missing,
    /// Tool installed but at a version outside our supported range.
    Outdated,
    /// Probe call returned an error other than "not on PATH" (e.g.
    /// permission denied, broken install).
    Broken,
    /// Adapter is registered but disabled in settings.
    Disabled,
}

impl AdapterAvailability {
    /// `true` for an adapter that's ready to run a case today.
    pub fn is_ready(self) -> bool {
        matches!(self, Self::Installed)
    }

    /// Short user-facing label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Installed => "Installed",
            Self::Missing => "Missing",
            Self::Outdated => "Outdated",
            Self::Broken => "Broken",
            Self::Disabled => "Disabled",
        }
    }
}

/// One adapter's row in the wizard report.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct AdapterStatus {
    /// Stable id (`"openfoam"`, `"gmsh"`, …).
    pub id: String,
    /// Human display name (`"OpenFOAM"`).
    pub display_name: String,
    /// Probe result.
    pub availability: AdapterAvailability,
    /// Detected version when [`AdapterAvailability::Installed`] /
    /// `Outdated` — empty otherwise.
    #[serde(default)]
    pub detected_version: Option<String>,
    /// Per-OS install hint shown when the adapter is missing.
    /// `None` when the adapter doesn't need a hint (built-in,
    /// already covered, etc.).
    #[serde(default)]
    pub install_hint: Option<String>,
}

/// Aggregated environment report.
///
/// The wizard renders this directly: rows in declaration order,
/// summary count at the top.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct EnvironmentReport {
    pub adapters: Vec<AdapterStatus>,
}

impl EnvironmentReport {
    /// `true` when at least one adapter is ready — the user has
    /// something runnable today.
    pub fn has_any_ready(&self) -> bool {
        self.adapters.iter().any(|a| a.availability.is_ready())
    }

    /// Count of ready adapters. Useful for the wizard's headline
    /// ("12 of 17 adapters ready to run").
    pub fn ready_count(&self) -> usize {
        self.adapters
            .iter()
            .filter(|a| a.availability.is_ready())
            .count()
    }

    /// Adapter count broken down by availability — feeds the
    /// summary block at the top of the wizard.
    pub fn count_by_availability(&self) -> BTreeMap<AdapterAvailability, usize> {
        let mut out: BTreeMap<AdapterAvailability, usize> = BTreeMap::new();
        for a in &self.adapters {
            *out.entry(a.availability).or_insert(0) += 1;
        }
        out
    }
}

// AdapterAvailability needs Ord for BTreeMap key use above.
impl PartialOrd for AdapterAvailability {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for AdapterAvailability {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Stable display order: Installed first, then Outdated /
        // Broken / Missing / Disabled.
        let rank = |a: &AdapterAvailability| match a {
            AdapterAvailability::Installed => 0,
            AdapterAvailability::Outdated => 1,
            AdapterAvailability::Broken => 2,
            AdapterAvailability::Missing => 3,
            AdapterAvailability::Disabled => 4,
        };
        rank(self).cmp(&rank(other))
    }
}

/// User's selection at the end of the wizard. Persisted to
/// `settings.json` under a top-level `first_run` key.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct FirstRunDecision {
    /// Whether the wizard has run (and therefore shouldn't auto-
    /// open again).
    #[serde(default)]
    pub completed: bool,
    /// Per-adapter "the user expects to use this" flag. Adapters
    /// not in the map fall back to "respect installed status" —
    /// a missing adapter is dimmed, an installed one is enabled.
    #[serde(default)]
    pub adapter_enabled: BTreeMap<String, bool>,
}

impl FirstRunDecision {
    /// Build a fresh decision from an [`EnvironmentReport`] —
    /// every Installed adapter starts checked, everything else
    /// starts unchecked. The user can flip individual rows in the
    /// wizard before clicking Done.
    pub fn from_report(report: &EnvironmentReport) -> Self {
        let mut adapter_enabled: BTreeMap<String, bool> = BTreeMap::new();
        for a in &report.adapters {
            adapter_enabled.insert(a.id.clone(), a.availability.is_ready());
        }
        Self {
            completed: false,
            adapter_enabled,
        }
    }

    /// `true` if this adapter was ticked in the wizard. Adapters
    /// the user never saw (added in a later release) default to
    /// `true` — the GUI surfaces them and the user can disable in
    /// Settings if they want.
    pub fn is_adapter_enabled(&self, id: &str) -> bool {
        self.adapter_enabled.get(id).copied().unwrap_or(true)
    }

    /// Mark the wizard complete. Called when the user clicks
    /// "Done" or "Skip".
    pub fn mark_completed(&mut self) {
        self.completed = true;
    }
}

/// Should the GUI auto-open the wizard at launch?
///
/// `decision` is what we read from `settings.json`. Returns `true`
/// only when the wizard has never run — every subsequent launch
/// skips. The user can re-open it manually from the command
/// palette ("Settings → Run first-launch wizard").
pub fn should_auto_show(decision: &FirstRunDecision) -> bool {
    !decision.completed
}

/// Per-OS install hint catalogue.
///
/// Keeping this in its own function (rather than hard-coded into
/// every adapter) means the wizard can render the right hint for
/// the current OS without each adapter knowing about Linux / macOS
/// / Windows package managers. Hints prefer the package manager
/// most users on that OS will recognise — apt on Linux,
/// homebrew on macOS, winget on Windows.
pub fn install_hint_for(adapter_id: &str) -> Option<String> {
    let os = std::env::consts::OS;
    install_hint_for_with_os(adapter_id, os)
}

/// Pure helper that takes the OS as an argument so unit tests can
/// pin it without monkey-patching `std::env::consts::OS`.
pub fn install_hint_for_with_os(adapter_id: &str, os: &str) -> Option<String> {
    let hint = match (adapter_id, os) {
        ("openfoam", "linux") => "sudo apt install openfoam",
        ("openfoam", "macos") => "brew install openfoam",
        ("openfoam", _) => "https://www.openfoam.com/download/",
        ("gmsh", "linux") => "sudo apt install gmsh",
        ("gmsh", "macos") => "brew install gmsh",
        ("gmsh", "windows") => "winget install gmsh",
        ("gmsh", _) => "https://gmsh.info/",
        ("netgen", _) => "https://ngsolve.org/downloads",
        ("freecad", "linux") => "sudo apt install freecad",
        ("freecad", "macos") => "brew install --cask freecad",
        ("freecad", "windows") => "winget install FreeCAD.FreeCAD",
        ("freecad", _) => "https://www.freecad.org/downloads.php",
        ("calculix", "linux") => "sudo apt install calculix-ccx",
        ("calculix", _) => "http://www.calculix.de/",
        ("elmer", "linux") => "sudo apt install elmerfem-csc",
        ("elmer", _) => "https://www.elmerfem.org/blog/binaries/",
        ("cantera", _) => "pip install cantera",
        ("lammps", "linux") => "sudo apt install lammps",
        ("lammps", "macos") => "brew install lammps",
        ("lammps", _) => "https://lammps.org/download.html",
        ("gromacs", "linux") => "sudo apt install gromacs",
        ("gromacs", "macos") => "brew install gromacs",
        ("gromacs", _) => "https://www.gromacs.org/Downloads",
        ("openems", _) => "https://openems.de/start/",
        ("meep", _) => "pip install meep",
        ("pybamm", _) => "pip install pybamm",
        ("mujoco", _) => "https://github.com/google-deepmind/mujoco/releases",
        ("precice", _) => "https://precice.org/installation-overview.html",
        ("su2", _) => "https://su2code.github.io/download.html",
        ("openradioss", _) => "https://github.com/OpenRadioss/OpenRadioss/releases",
        ("code-aster", _) => "https://www.code-aster.org/V2/spip.php?article272",
        ("occt", "linux") => "sudo apt install libocct-foundation-dev",
        ("occt", _) => "https://dev.opencascade.org/release",
        _ => return None,
    };
    Some(hint.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report_for(rows: &[(&str, AdapterAvailability)]) -> EnvironmentReport {
        EnvironmentReport {
            adapters: rows
                .iter()
                .map(|(id, av)| AdapterStatus {
                    id: (*id).to_string(),
                    display_name: id.to_uppercase(),
                    availability: *av,
                    detected_version: None,
                    install_hint: None,
                })
                .collect(),
        }
    }

    #[test]
    fn availability_is_ready_only_for_installed() {
        assert!(AdapterAvailability::Installed.is_ready());
        assert!(!AdapterAvailability::Missing.is_ready());
        assert!(!AdapterAvailability::Outdated.is_ready());
        assert!(!AdapterAvailability::Broken.is_ready());
        assert!(!AdapterAvailability::Disabled.is_ready());
    }

    #[test]
    fn availability_label_is_human_readable() {
        assert_eq!(AdapterAvailability::Installed.label(), "Installed");
        assert_eq!(AdapterAvailability::Outdated.label(), "Outdated");
    }

    #[test]
    fn report_has_any_ready_only_when_one_is_installed() {
        let none = report_for(&[
            ("a", AdapterAvailability::Missing),
            ("b", AdapterAvailability::Broken),
        ]);
        assert!(!none.has_any_ready());

        let some = report_for(&[
            ("a", AdapterAvailability::Missing),
            ("b", AdapterAvailability::Installed),
        ]);
        assert!(some.has_any_ready());
        assert_eq!(some.ready_count(), 1);
    }

    #[test]
    fn report_count_by_availability_groups_correctly() {
        let r = report_for(&[
            ("a", AdapterAvailability::Installed),
            ("b", AdapterAvailability::Installed),
            ("c", AdapterAvailability::Missing),
            ("d", AdapterAvailability::Broken),
        ]);
        let counts = r.count_by_availability();
        assert_eq!(counts[&AdapterAvailability::Installed], 2);
        assert_eq!(counts[&AdapterAvailability::Missing], 1);
        assert_eq!(counts[&AdapterAvailability::Broken], 1);
    }

    #[test]
    fn decision_from_report_pre_ticks_installed_only() {
        let r = report_for(&[
            ("openfoam", AdapterAvailability::Installed),
            ("gmsh", AdapterAvailability::Missing),
            ("calculix", AdapterAvailability::Installed),
        ]);
        let d = FirstRunDecision::from_report(&r);
        assert!(d.is_adapter_enabled("openfoam"));
        assert!(!d.is_adapter_enabled("gmsh"));
        assert!(d.is_adapter_enabled("calculix"));
        // Adapter the user never saw — defaults to enabled so a
        // future release doesn't silently hide its new adapters.
        assert!(d.is_adapter_enabled("brand-new-adapter"));
    }

    #[test]
    fn decision_default_is_not_completed() {
        let d = FirstRunDecision::default();
        assert!(!d.completed);
        assert!(should_auto_show(&d));
    }

    #[test]
    fn decision_mark_completed_suppresses_auto_show() {
        let mut d = FirstRunDecision::default();
        d.mark_completed();
        assert!(d.completed);
        assert!(!should_auto_show(&d));
    }

    #[test]
    fn install_hint_for_with_os_returns_per_os_string() {
        assert_eq!(
            install_hint_for_with_os("openfoam", "linux"),
            Some("sudo apt install openfoam".to_string())
        );
        assert_eq!(
            install_hint_for_with_os("gmsh", "macos"),
            Some("brew install gmsh".to_string())
        );
        assert_eq!(
            install_hint_for_with_os("gmsh", "windows"),
            Some("winget install gmsh".to_string())
        );
    }

    #[test]
    fn install_hint_for_with_os_falls_back_to_url() {
        // An OS we don't have a package-manager hint for falls
        // through to the upstream URL.
        let hint = install_hint_for_with_os("openfoam", "freebsd").unwrap();
        assert!(hint.starts_with("http"), "got: {hint}");
    }

    #[test]
    fn install_hint_for_with_os_returns_none_for_unknown_adapter() {
        assert!(install_hint_for_with_os("not-a-real-adapter", "linux").is_none());
    }

    #[test]
    fn decision_round_trips_through_serde() {
        let r = report_for(&[
            ("openfoam", AdapterAvailability::Installed),
            ("gmsh", AdapterAvailability::Missing),
        ]);
        let mut d = FirstRunDecision::from_report(&r);
        d.mark_completed();
        let s = serde_json::to_string(&d).unwrap();
        let back: FirstRunDecision = serde_json::from_str(&s).unwrap();
        assert_eq!(back, d);
    }

    #[test]
    fn report_round_trips_through_serde() {
        let r = report_for(&[
            ("openfoam", AdapterAvailability::Installed),
            ("calculix", AdapterAvailability::Outdated),
        ]);
        let s = serde_json::to_string(&r).unwrap();
        let back: EnvironmentReport = serde_json::from_str(&s).unwrap();
        assert_eq!(back, r);
    }

    #[test]
    fn install_hints_cover_every_workspace_adapter() {
        // Drift guard — adding a new adapter to the workspace
        // without giving it an install hint here should surface
        // here so the wizard isn't left with a blank cell.
        let workspace_adapters = [
            "openfoam",
            "su2",
            "gmsh",
            "netgen",
            "freecad",
            "calculix",
            "elmer",
            "code-aster",
            "openradioss",
            "cantera",
            "lammps",
            "gromacs",
            "openems",
            "meep",
            "pybamm",
            "mujoco",
            "precice",
            "occt",
        ];
        for id in workspace_adapters {
            let hint = install_hint_for_with_os(id, "linux");
            assert!(hint.is_some(), "no Linux install hint for `{id}`");
        }
    }
}
