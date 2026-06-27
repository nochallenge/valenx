//! Baked-in **product self-test** — one fast, token-cheap command that
//! exercises every valenx product's real compute path with representative
//! inputs, reads its output, and checks it is sane (finite, non-error), with
//! exact known-value assertions where ground truth exists.
//!
//! # Why this exists
//!
//! Verifying valenx's 56 products used to be hand-scripted each run (open a
//! workbench, drive it, eyeball the readout). That is slow and burns tokens.
//! This module bakes the driving + checking into the binary so a single
//! trigger emits a **compact** per-product report (one line each:
//! `id  PASS|FAIL|SKIP  <key value or reason>`) plus a summary tally — cheap to
//! read, no walls of text.
//!
//! # How to invoke
//!
//! Headless CLI, no GUI window and no `rfd` file dialog (both would hang a
//! head-less / CI run):
//!
//! ```text
//! valenx --self-test                  # every product
//! valenx --self-test --group Simulation
//! valenx --self-test --id thermo
//! ```
//!
//! Wired in [`crate::headless::run_headless`] under the `self-test` task.
//!
//! # Mechanism
//!
//! Each check runs against a fresh [`ValenxApp::default()`] — the SAME pure-CPU
//! construction the workbench unit tests use (no GPU, no window). Two tiers:
//!
//! * **Deep** (exemplars with known ground truth): drive the workbench's real
//!   compute (`thermo`/`quantum`/`optics`/`acoustics` run paths,
//!   `waveform` parse) with representative inputs, read the structured
//!   `agent_readout()`, and assert the key number against the documented value
//!   with a tolerance.
//! * **Generic** (every other drivable product): open the product's tab and
//!   render its panel in a throwaway headless `accesskit` frame (the SAME probe
//!   `read_text` / the AI-drivability tests use, via
//!   [`crate::agent_commands::probe_active_workbench_text`]). This executes the
//!   workbench's real draw/compute path; PASS iff it renders substantive text
//!   with no `NaN`/`inf`/error token.
//! * **Skip**: a product that genuinely cannot self-verify head-less (GPU
//!   render, external tool, interactive scan/file input) is reported `SKIP`
//!   with a one-word reason — never a faked pass.
//!
//! # Extending — adding a product check is **one registry entry**
//!
//! Add a row to [`product_checks`]:
//!
//! ```ignore
//! ProductCheck { id: "myproduct", kind: TabKind::MyProduct, mode: CheckMode::Generic }
//! ```
//!
//! …or `CheckMode::Deep(my_deep_fn)` for a known-value assertion, or
//! `CheckMode::Skip("reason")`. The [`registry_covers_every_template`] test
//! guarantees the registry stays in lock-step with `TabKind::TEMPLATES` (the
//! authoritative 56), so a newly-added product fails CI until it has a row.

use crate::project_tabs::TabKind;
use crate::ValenxApp;

/// PASS / FAIL / SKIP for one product check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// The compute path ran and its output passed the check.
    Pass,
    /// The compute path ran but its output failed (empty, non-finite, wrong
    /// value, or an error line).
    Fail,
    /// The product cannot self-verify head-less (GPU / external tool / file
    /// input). Not a pass and not a failure.
    Skip,
}

impl Status {
    /// Three-letter tag for the compact report.
    fn tag(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Fail => "FAIL",
            Status::Skip => "SKIP",
        }
    }
}

/// The result of one product check: a [`Status`] plus a short detail (the key
/// value on PASS, the reason on FAIL/SKIP).
#[derive(Debug, Clone)]
pub struct CheckOutcome {
    pub status: Status,
    pub detail: String,
}

impl CheckOutcome {
    fn pass(detail: impl Into<String>) -> Self {
        Self {
            status: Status::Pass,
            detail: detail.into(),
        }
    }
    fn fail(detail: impl Into<String>) -> Self {
        Self {
            status: Status::Fail,
            detail: detail.into(),
        }
    }
    fn skip(reason: impl Into<String>) -> Self {
        Self {
            status: Status::Skip,
            detail: reason.into(),
        }
    }
}

/// How a product is checked.
pub enum CheckMode {
    /// A known-value assertion against documented ground truth. The closure
    /// gets a fresh app and returns the outcome.
    Deep(fn(&mut ValenxApp) -> CheckOutcome),
    /// Open the tab, probe the rendered panel, PASS iff it renders substantive
    /// finite text with no error token.
    Generic,
    /// Cannot self-verify head-less; the `&str` is the one-word reason.
    Skip(&'static str),
}

/// One product's self-test registration: its canonical id (the
/// `TabKind::from_id` string, authoritative in `docs/PRODUCTS.md`), the tab
/// kind it opens, and how it is checked.
pub struct ProductCheck {
    pub id: &'static str,
    pub kind: TabKind,
    pub mode: CheckMode,
}

/// One line of the final report.
#[derive(Debug, Clone)]
pub struct ReportLine {
    pub id: String,
    pub status: Status,
    pub detail: String,
}

/// The full self-test report: one [`ReportLine`] per checked product, in
/// registry order, plus the pass/fail/skip tally.
#[derive(Debug, Clone, Default)]
pub struct SelfTestReport {
    pub lines: Vec<ReportLine>,
    pub passed: usize,
    pub failed: usize,
    pub skipped: usize,
}

impl SelfTestReport {
    /// Render the compact report as text: one `id  STATUS  detail` line per
    /// product (id left-padded to a stable column) then a one-line tally.
    /// Deliberately terse — cheap to read, no walls of text.
    pub fn render(&self) -> String {
        let mut out = String::new();
        // Stable id column width so the STATUS column lines up.
        let w = self
            .lines
            .iter()
            .map(|l| l.id.len())
            .max()
            .unwrap_or(0)
            .max(2);
        for l in &self.lines {
            out.push_str(&format!(
                "{:<width$}  {}  {}\n",
                l.id,
                l.status.tag(),
                l.detail,
                width = w
            ));
        }
        out.push_str(&format!(
            "— {} product(s): {} PASS · {} FAIL · {} SKIP\n",
            self.lines.len(),
            self.passed,
            self.failed,
            self.skipped
        ));
        out
    }

    /// `true` iff no product check FAILED (skips are not failures). The exit
    /// signal a CI gate keys on.
    pub fn ok(&self) -> bool {
        self.failed == 0
    }
}

/// Which products to run, parsed from the CLI flags after `self-test`.
#[derive(Debug, Clone, Default)]
pub enum Filter {
    /// Every product in the registry.
    #[default]
    All,
    /// Only products whose [`TabKind::group`] matches (case-insensitive).
    Group(String),
    /// Only the single product with this canonical id.
    Id(String),
}

impl Filter {
    /// Parse `--group <G>` / `--id <id>` out of the CLI tokens that follow
    /// `self-test`. An explicit `--id` wins over `--group`; neither ⇒ [`All`].
    pub fn from_args(args: &[String]) -> Filter {
        let mut group: Option<String> = None;
        let mut id: Option<String> = None;
        let mut i = 0;
        while i < args.len() {
            match args[i].as_str() {
                "--group" | "-g" => {
                    group = args.get(i + 1).cloned();
                    i += 2;
                }
                "--id" | "-i" => {
                    id = args.get(i + 1).cloned();
                    i += 2;
                }
                _ => i += 1,
            }
        }
        if let Some(id) = id {
            Filter::Id(id)
        } else if let Some(g) = group {
            Filter::Group(g)
        } else {
            Filter::All
        }
    }

    /// Does product `pc` pass this filter?
    fn matches(&self, pc: &ProductCheck) -> bool {
        match self {
            Filter::All => true,
            Filter::Group(g) => pc.kind.group().eq_ignore_ascii_case(g.trim()),
            Filter::Id(id) => pc.id.eq_ignore_ascii_case(id.trim()),
        }
    }
}

// ===========================================================================
// The registry — one entry per product. THE place to add a product check.
// ===========================================================================

/// Every product's self-test registration, in `TabKind::TEMPLATES` order.
///
/// The single source of truth for what the self-test drives. Keep one row per
/// `TabKind::TEMPLATES` entry; [`registry_covers_every_template`] enforces it.
pub fn product_checks() -> Vec<ProductCheck> {
    use CheckMode::{Deep, Generic, Skip};
    use TabKind as K;
    vec![
        // -- Aerospace --
        // rocket: DEEP — Tsiolkovsky Δv of the LV-1 vs the analytic 2-stage sum.
        ProductCheck {
            id: "rocket",
            kind: K::Rocket,
            mode: Deep(check_rocket),
        },
        // engine: DEEP — characteristic velocity c* vs 1-D Vandenkerckhove theory.
        ProductCheck {
            id: "engine",
            kind: K::Engine,
            mode: Deep(check_engine),
        },
        // astro: DEEP — Hohmann LEO→GEO total Δv vs closed-form Kepler/vis-viva.
        ProductCheck {
            id: "astro",
            kind: K::Astro,
            mode: Deep(check_astro),
        },
        // note: aero stays Generic — the aero WORKBENCH drives an async (background
        // thread) RANS CFD solve on a demo box, whose Cl is a numerical box result,
        // NOT the thin-airfoil 2π·α analytic. The `2π` thin-airfoil slope is a
        // valenx-aero *crate* benchmark, not the product's drive path, so a clean
        // known-value assertion would test the library, not the product. Honest
        // call: keep the generic panel-renders-finite check here.
        ProductCheck {
            id: "aero",
            kind: K::Aero,
            mode: Generic,
        },
        // gasdynamics: DEEP — isentropic A/A* + normal-shock M2, p2/p1 (exact).
        ProductCheck {
            id: "gasdynamics",
            kind: K::Gasdynamics,
            mode: Deep(check_gasdynamics),
        },
        // rotor: DEEP — BEMT hover figure of merit vs momentum theory P_ideal/P.
        ProductCheck {
            id: "rotor",
            kind: K::Rotor,
            mode: Deep(check_rotor),
        },
        // uas: DEEP — multirotor hover power vs disk-loading momentum theory.
        ProductCheck {
            id: "uas",
            kind: K::Uas,
            mode: Deep(check_uas),
        },
        // -- Astrophysics --
        // blackhole: DEEP — Schwarzschild r₊=2M, photon=3M, ISCO=6M, shadow=3√3 M.
        ProductCheck {
            id: "blackhole",
            kind: K::BlackHole,
            mode: Deep(check_blackhole),
        },
        // -- Simulation --
        // note: cfd stays Generic — the default case is the lid-driven cavity at
        // Re=100, whose validation is the *tabular* Ghia-1982 centreline-u
        // benchmark (no clean closed form), and the case selector (cavity vs the
        // plane-Poiseuille channel, which DOES have an analytic parabola) is not
        // exposed via the agent control API to switch to it. Honest call: keep
        // the generic panel-renders-finite check rather than assert a weak bound.
        ProductCheck {
            id: "cfd",
            kind: K::Cfd,
            mode: Generic,
        },
        // fem: DEEP — 3-D linear-static cantilever tip δ vs Euler–Bernoulli PL³/3EI.
        ProductCheck {
            id: "fem",
            kind: K::Fem,
            mode: Deep(check_fem),
        },
        ProductCheck {
            id: "topopt",
            kind: K::TopOpt,
            mode: Generic,
        },
        // nodegraph: DEEP — default Constant(2)+Constant(2)→Add→Output = 4.0.
        ProductCheck {
            id: "nodegraph",
            kind: K::NodeGraph,
            mode: Deep(check_nodegraph),
        },
        ProductCheck {
            id: "bondgraph",
            kind: K::BondGraph,
            mode: Generic,
        },
        ProductCheck {
            id: "surrogate",
            kind: K::Surrogate,
            mode: Generic,
        },
        // reactdyn: DEEP — H2 AIMD NVE energy conservation |ΔE/E₀| < 1e-2.
        ProductCheck {
            id: "reactdyn",
            kind: K::Reactdyn,
            mode: Deep(check_reactdyn),
        },
        // thermo: DEEP — CO₂ Peng-Robinson vapor Z + Psat against NIST range.
        ProductCheck {
            id: "thermo",
            kind: K::Thermo,
            mode: Deep(check_thermo),
        },
        // quantum: DEEP — Bell pair p(|00>) = p(|11>) ≈ 0.5.
        ProductCheck {
            id: "quantum",
            kind: K::Quantum,
            mode: Deep(check_quantum),
        },
        // fields: DEEP — descriptive stats of 1..5 (mean=3, σ_pop=√2) exact.
        ProductCheck {
            id: "fields",
            kind: K::Fields,
            mode: Deep(check_fields),
        },
        ProductCheck {
            id: "fluids",
            kind: K::Fluids,
            mode: Generic,
        },
        ProductCheck {
            id: "ocean",
            kind: K::Ocean,
            mode: Generic,
        },
        ProductCheck {
            id: "rom",
            kind: K::Rom,
            mode: Generic,
        },
        ProductCheck {
            id: "uq",
            kind: K::Uq,
            mode: Generic,
        },
        ProductCheck {
            id: "missionsim",
            kind: K::MissionSim,
            mode: Generic,
        },
        ProductCheck {
            id: "missionplanner",
            kind: K::MissionPlanner,
            mode: Generic,
        },
        ProductCheck {
            id: "survivability",
            kind: K::Survivability,
            mode: Generic,
        },
        // cosim: external FMI/HELICS co-simulation tool.
        ProductCheck {
            id: "cosim",
            kind: K::Cosim,
            mode: Skip("external-tool"),
        },
        // mbd: DEEP — pendulum period vs 2π√(L/g) + NVE energy conservation.
        ProductCheck {
            id: "mbd",
            kind: K::Mbd,
            mode: Deep(check_mbd),
        },
        // optics: DEEP — thin-lens magnification sign (object at 2f ⇒ m ≈ -1).
        ProductCheck {
            id: "optics",
            kind: K::Optics,
            mode: Deep(check_optics),
        },
        // acoustics: DEEP — monopole 1/r law (p at 2r is half p at r).
        ProductCheck {
            id: "acoustics",
            kind: K::Acoustics,
            mode: Deep(check_acoustics),
        },
        // waveform: DEEP — sample VCD parses to clk + cnt (2 signals).
        ProductCheck {
            id: "waveform",
            kind: K::Waveform,
            mode: Deep(check_waveform),
        },
        // -- CAD & mesh --
        // cad: DEEP — parametric CSG volume of unit box − Ø0.5 hole = 1−π/16 (exact).
        ProductCheck {
            id: "cad",
            kind: K::Cad,
            mode: Deep(check_cad),
        },
        // brep: DEEP — boolean difference volume ≈ 1−π/16 (tessellated, 3% band).
        ProductCheck {
            id: "brep",
            kind: K::BrepCad,
            mode: Deep(check_brep),
        },
        // mesh: DEEP — canonical LV-1 mesh AABB extents (9.6×9.6×34.7).
        ProductCheck {
            id: "mesh",
            kind: K::MeshToolbox,
            mode: Deep(check_mesh),
        },
        // sheetmetal: DEEP — bend allowance (π·θ/180)·(R+K·t) = 2.2619 mm (exact).
        ProductCheck {
            id: "sheetmetal",
            kind: K::Sheetmetal,
            mode: Deep(check_sheetmetal),
        },
        // reverse: DEEP — unit-sphere recon (600 pts, bbox ≈ 2.0³).
        ProductCheck {
            id: "reverse",
            kind: K::Reverse,
            mode: Deep(check_reverse),
        },
        // photogrammetry: STAYS SKIP — needs an imported image set / scan.
        ProductCheck {
            id: "photogrammetry",
            kind: K::Photogrammetry,
            mode: Skip("scan-input"),
        },
        // draft2d: DEEP — exact 2-D geometry (4 entities, extent 60×40 units).
        ProductCheck {
            id: "draft2d",
            kind: K::Draft2d,
            mode: Deep(check_draft2d),
        },
        // render: STAYS SKIP — path-traced raster output (GPU pixels, not a readout).
        ProductCheck {
            id: "render",
            kind: K::Render,
            mode: Skip("gpu-render"),
        },
        // animate: DEEP — keyframe interpolation (Linear t=1s = π/2; endpoints 0,π).
        ProductCheck {
            id: "animate",
            kind: K::Animate,
            mode: Deep(check_animate),
        },
        // -- Machine design --
        // springs: DEEP — helical spring rate k = G·d⁴/(8·D³·n) (closed form).
        ProductCheck {
            id: "springs",
            kind: K::Springs,
            mode: Deep(check_springs),
        },
        // gears: DEEP — spur pitch diameter m·z/cosβ + meshing ratio z2/z1.
        ProductCheck {
            id: "gears",
            kind: K::Gears,
            mode: Deep(check_gears),
        },
        // fasteners: DEEP — M6 pitch dia d−0.6495P + ISO-898 tensile stress area.
        ProductCheck {
            id: "fasteners",
            kind: K::Fasteners,
            mode: Deep(check_fasteners),
        },
        // frames: DEEP — IPE 200 cross-section area 2·b·tf+tw·(h−2tf) (the panel
        // computes area + perimeter, not I/S — see check_frames note).
        ProductCheck {
            id: "frames",
            kind: K::Frames,
            mode: Deep(check_frames),
        },
        // collision: DEEP — AABB overlap both ways (disjoint gap=10 / overlap sep=0).
        ProductCheck {
            id: "collision",
            kind: K::Collision,
            mode: Deep(check_collision),
        },
        // -- Civil & AEC --
        // piping: DEEP — ASME B36.10 NPS 2 Sch 40 OD=60.325mm + bore (exact).
        ProductCheck {
            id: "piping",
            kind: K::Piping,
            mode: Deep(check_piping),
        },
        // hvac: DEEP — Darcy–Weisbach duct ΔP = f·(L/D)·½ρV² (15.05 Pa).
        ProductCheck {
            id: "hvac",
            kind: K::Hvac,
            mode: Deep(check_hvac),
        },
        // note: reinforcement stays Generic — the workbench is a parametric
        // rebar-CAGE GEOMETRY generator (width/depth/length/n_bars/hoop-spacing →
        // a 3-D cage mesh), NOT a structural calculator: it exposes no As, fy,
        // f'c, bar diameter, nominal moment Mn, or reinforcement ratio ρ, so there
        // is no closed-form structural quantity to assert against. Honest call:
        // keep the generic panel-renders-finite check.
        ProductCheck {
            id: "reinforcement",
            kind: K::Reinforcement,
            mode: Generic,
        },
        // interior: DEEP — floor-plan bookkeeping (6×4 m room, 2 pieces) exact.
        ProductCheck {
            id: "interior",
            kind: K::Interior,
            mode: Deep(check_interior),
        },
        // geomatics: DEEP — Cambridge→Paris great-circle (haversine) ≈ 404.3 km.
        ProductCheck {
            id: "geomatics",
            kind: K::Geomatics,
            mode: Deep(check_geomatics),
        },
        // -- Life sciences --
        // genetics: DEEP — exact sequence GC content + length (ATGC → 50%, 4).
        ProductCheck {
            id: "genetics",
            kind: K::Genetics,
            mode: Deep(check_genetics),
        },
        // neuro: DEEP — unmyelinated conduction velocity v = k·√d (cable theory).
        ProductCheck {
            id: "neuro",
            kind: K::Neuro,
            mode: Deep(check_neuro),
        },
        // variant: DEEP — HGVS parse p.R273H → ProteinSub R(273)H (exact).
        ProductCheck {
            id: "variant",
            kind: K::VariantEffect,
            mode: Deep(check_variant),
        },
        // ppi: DEEP — shortest-path BFS GUARD→EFF-A = 1 hop (exact graph result).
        ProductCheck {
            id: "ppi",
            kind: K::Ppi,
            mode: Deep(check_ppi),
        },
        // morphogenesis: DEEP — Gray–Scott bounded-stability invariant U,V∈[0,1].
        ProductCheck {
            id: "morphogenesis",
            kind: K::Morphogenesis,
            mode: Deep(check_morphogenesis),
        },
        // -- Sensors --
        // sensors: DEEP — LiDAR angular resolution = FOV/(N−1) = 90°/31 (exact).
        ProductCheck {
            id: "sensors",
            kind: K::Sensors,
            mode: Deep(check_sensors),
        },
        // autonomy: DEEP — V&V min clearance = point-to-path distance (3.5 m) + PASS.
        ProductCheck {
            id: "autonomy",
            kind: K::Autonomy,
            mode: Deep(check_autonomy),
        },
    ]
}

// ===========================================================================
// Runner + compact report
// ===========================================================================

/// Run the product self-tests selected by `filter` and return the compact
/// [`SelfTestReport`]. Each product gets a fresh [`ValenxApp::default()`] so one
/// product's state can never leak into another's check.
pub fn run_self_tests(filter: &Filter) -> SelfTestReport {
    let mut report = SelfTestReport::default();
    for pc in product_checks().into_iter().filter(|pc| filter.matches(pc)) {
        let outcome = run_one(&pc);
        match outcome.status {
            Status::Pass => report.passed += 1,
            Status::Fail => report.failed += 1,
            Status::Skip => report.skipped += 1,
        }
        report.lines.push(ReportLine {
            id: pc.id.to_string(),
            status: outcome.status,
            detail: outcome.detail,
        });
    }
    report
}

/// Run a single product's check on a fresh app. `Deep` runs its closure;
/// `Generic` opens the tab and probes the panel; `Skip` reports the reason
/// without constructing anything.
fn run_one(pc: &ProductCheck) -> CheckOutcome {
    match &pc.mode {
        CheckMode::Skip(reason) => CheckOutcome::skip(*reason),
        CheckMode::Deep(f) => {
            let mut app = ValenxApp::default();
            f(&mut app)
        }
        CheckMode::Generic => {
            let mut app = ValenxApp::default();
            generic_check(&mut app, pc.kind)
        }
    }
}

/// The universal product check: open `kind`'s tab (so it is the active
/// workbench), render its panel head-less, and inspect the readable text.
///
/// This executes the workbench's real draw/compute path. PASS iff the panel
/// renders **substantive** text (so the product actually loaded + drew) with no
/// `NaN` / `inf` / error token. An empty / panicked panel ⇒ FAIL.
fn generic_check(app: &mut ValenxApp, kind: TabKind) -> CheckOutcome {
    app.tab_bar.open(kind);
    let Some(text) = crate::agent_commands::probe_active_workbench_text(app) else {
        return CheckOutcome::fail("panel emitted no accessible text");
    };
    // Require some substance: a stub that drew nothing useful is not a pass.
    let joined = text.join(" • ");
    if joined.trim().is_empty() || text.len() < 2 {
        return CheckOutcome::fail(format!("panel text too sparse ({} node(s))", text.len()));
    }
    if let Some(bad) = non_finite_or_error_token(&text) {
        return CheckOutcome::fail(format!("non-finite/error in panel: {bad:?}"));
    }
    // PASS — quote a short, representative slice of the rendered text as the key
    // value so the report line carries evidence the panel really computed.
    CheckOutcome::pass(short_evidence(&joined))
}

/// Scan probed panel text for a clearly-bad signal and return the offending
/// fragment (or `None` if clean). Two classes:
///
/// * a **non-finite numeric token** — a standalone `NaN` / `inf` / `infinity`
///   value (word-boundaried so benign substrings like "**inf**ormation" don't
///   trip it);
/// * a **compute-failure phrase** — `failed`, `error:` (with the colon, as the
///   workbenches' `"… failed: {e}"` / `"error: {e}"` lines emit) or `panicked`.
///
/// Critically, the **bare word "error"** is NOT a failure marker: many panels
/// legitimately label a metric "Relative error", "reproj error", "coupling
/// error" etc., and flagging those would false-FAIL a working product. Only the
/// `error:`-with-colon failure format counts.
fn non_finite_or_error_token(parts: &[String]) -> Option<String> {
    for p in parts {
        let lower = p.to_ascii_lowercase();
        // Compute-failure phrases (substring is fine — these are unambiguous).
        if lower.contains("failed") || lower.contains("error:") || lower.contains("panicked") {
            return Some(p.clone());
        }
        // Non-finite numeric tokens, word-boundaried.
        for tok in lower.split(|c: char| !(c.is_ascii_alphanumeric() || c == '.' || c == '_')) {
            let t = tok.trim_matches('.');
            if matches!(t, "nan" | "inf" | "-inf" | "infinity") {
                return Some(p.clone());
            }
        }
    }
    None
}

/// A short evidence slice for a PASS report line: the first chunk of the panel
/// text, length-capped so the report stays cheap to read.
fn short_evidence(joined: &str) -> String {
    const MAX: usize = 90;
    let one_line = joined.replace('\n', " ");
    if one_line.chars().count() > MAX {
        let truncated: String = one_line.chars().take(MAX).collect();
        format!("panel ok: {truncated}…")
    } else {
        format!("panel ok: {one_line}")
    }
}

// ===========================================================================
// Deep known-value checks (exemplars with documented ground truth)
// ===========================================================================

/// Pull the float that follows `key` in a structured readout line, e.g.
/// `parse_keyed("…Z_vap=0.94321 …", "Z_vap=")` ⇒ `Some(0.94321)`. Stops at the
/// first non-numeric char (so a trailing unit / separator is ignored).
fn parse_keyed(readout: &str, key: &str) -> Option<f64> {
    let rest = &readout[readout.find(key)? + key.len()..];
    let num: String = rest
        .chars()
        .take_while(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | '+' | 'e' | 'E'))
        .collect();
    num.parse::<f64>().ok()
}

/// Pull the float from a `label … : value` line in a multi-section readout —
/// e.g. `parse_after_colon("… A/A*             : 1.687500\n…", "A/A*")` ⇒
/// `Some(1.6875)`. Finds `label`, then the next `:` after it, then the first
/// numeric run. Used for the monospace gas-dynamics report whose columns are
/// padded with spaces between the label and the `:`.
fn parse_after_colon(readout: &str, label: &str) -> Option<f64> {
    let from_label = &readout[readout.find(label)? + label.len()..];
    let after_colon = &from_label[from_label.find(':')? + 1..];
    let trimmed = after_colon.trim_start();
    let num: String = trimmed
        .chars()
        .take_while(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | '+' | 'e' | 'E'))
        .collect();
    num.parse::<f64>().ok()
}

// ---------------------------------------------------------------------------
// Aerospace deep checks — each asserts a real physics/analytic result against
// the EXACT quantity the workbench computes (not just "runs + sane").
// ---------------------------------------------------------------------------

/// **gasdynamics** — 1-D compressible-flow relations are textbook-exact. At the
/// default `M = 2.0`, `γ = 1.4` the workbench's report carries (NACA-1135):
/// isentropic `A/A* = 1.6875`, and the normal-shock `M2 = 0.57735`,
/// `p2/p1 = 4.5`. Asserts all three against the closed-form values (tol 1e-3),
/// driving the SAME `run_gasdynamics` the panel's Compute button fires.
fn check_gasdynamics(app: &mut ValenxApp) -> CheckOutcome {
    crate::gasdynamics_workbench::run(app);
    let Some(r) = app.gasdynamics.agent_readout() else {
        return CheckOutcome::fail("no readout after compute");
    };
    let area = parse_after_colon(&r, "A/A*");
    let m2 = parse_after_colon(&r, "M2");
    let p2 = parse_after_colon(&r, "p2/p1");
    let (Some(area), Some(m2), Some(p2)) = (area, m2, p2) else {
        return CheckOutcome::fail(format!("could not parse A/A*, M2, p2/p1 from: {r}"));
    };
    let bad =
        (area - 1.6875).abs() >= 1e-3 || (m2 - 0.577_350).abs() >= 1e-3 || (p2 - 4.5).abs() >= 1e-3;
    if bad || !(area.is_finite() && m2.is_finite() && p2.is_finite()) {
        return CheckOutcome::fail(format!(
            "M=2,γ=1.4 expected A/A*=1.6875 M2=0.5774 p2/p1=4.5, got {area:.4}/{m2:.4}/{p2:.4}"
        ));
    }
    CheckOutcome::pass(format!(
        "M2 A/A*={area:.4} shock M2={m2:.4} p2/p1={p2:.3} (NACA-1135 exact)"
    ))
}

/// **astro** — the default Hohmann transfer (LEO 300 km → GEO 35 786 km) is
/// closed-form orbital mechanics. Independently recomputes the total Δv from the
/// crate's own constants (`μ⊕`, `R⊕`) — `Δv = |v_peri − v_circ1| +
/// |v_circ2 − v_apo|` with vis-viva on the transfer ellipse — and asserts the
/// workbench's reported `total Δv` matches (tol 2 m/s, the readout rounds to
/// whole m/s). Drives the real `astro_product()` builder.
fn check_astro(_app: &mut ValenxApp) -> CheckOutcome {
    use valenx_astro::constants::{MU_EARTH, R_EARTH};
    // The default planner altitudes the product builder uses.
    let r1 = R_EARTH + 300.0 * 1_000.0;
    let r2 = R_EARTH + 35_786.0 * 1_000.0;
    // Closed-form Hohmann total Δv (same formulae as valenx_astro::hohmann_transfer).
    let v1 = (MU_EARTH / r1).sqrt();
    let v2 = (MU_EARTH / r2).sqrt();
    let a_t = 0.5 * (r1 + r2);
    let v_peri = (MU_EARTH * (2.0 / r1 - 1.0 / a_t)).sqrt();
    let v_apo = (MU_EARTH * (2.0 / r2 - 1.0 / a_t)).sqrt();
    let expected = (v_peri - v1).abs() + (v2 - v_apo).abs();

    // The workbench's genuine product readout (built from valenx_astro::hohmann_transfer).
    let product = crate::astro_workbench::astro_product();
    let line = product
        .lines
        .iter()
        .find(|l| l.contains("total Δv"))
        .cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!(
            "no 'total Δv' line in astro product: {:?}",
            product.lines
        ));
    };
    let Some(reported) = parse_after_colon(&line, "total Δv") else {
        return CheckOutcome::fail(format!("could not parse total Δv from: {line}"));
    };
    if !reported.is_finite() || (reported - expected).abs() > 2.0 {
        return CheckOutcome::fail(format!(
            "Hohmann total Δv: workbench {reported:.0} vs analytic {expected:.0} m/s (>2 m/s)"
        ));
    }
    CheckOutcome::pass(format!(
        "Hohmann LEO→GEO Δv={reported:.0} m/s (Kepler/vis-viva {expected:.0})"
    ))
}

/// **rocket** — Tsiolkovsky `Δv = Σ Isp_vac·g₀·ln(m₀/m_f)`. Independently
/// recomputes the LV-1's ideal Δv budget from its (stable, source-documented)
/// two-stage mass + Isp data, using the crate's `g₀`, and asserts the rocket
/// product's reported `Δv` (the live `valenx_astro::simulate_ascent` ideal Δv of
/// the LV-1) matches (tol 3 m/s, the readout rounds to whole m/s). Drives the
/// real `rocket_product()` builder.
fn check_rocket(_app: &mut ValenxApp) -> CheckOutcome {
    use valenx_astro::constants::G0;
    // LV-1 stages (see rocket_workbench::lv1_vehicle): payload 2000 kg. Tsiolkovsky
    // uses the VACUUM Isp (the crate's Stage::ideal_delta_v uses isp_vac).
    // Stage 0: dry 6000, prop 90000, Isp_vac 311. Stage 1: dry 1500, prop 12000, Isp_vac 345.
    let payload: f64 = 2_000.0;
    // Stage 1 (upper) burns with only payload above it.
    let (s1_dry, s1_prop, s1_isp): (f64, f64, f64) = (1_500.0, 12_000.0, 345.0);
    let upper1 = payload;
    let dv1 = s1_isp * G0 * ((upper1 + s1_dry + s1_prop) / (upper1 + s1_dry)).ln();
    // Stage 0 (booster) carries the full wet upper stage + payload above it.
    let (s0_dry, s0_prop, s0_isp): (f64, f64, f64) = (6_000.0, 90_000.0, 311.0);
    let upper0 = payload + s1_dry + s1_prop;
    let dv0 = s0_isp * G0 * ((upper0 + s0_dry + s0_prop) / (upper0 + s0_dry)).ln();
    let expected = dv0 + dv1;

    let product = crate::rocket_workbench::rocket_product();
    let line = product.lines.iter().find(|l| l.contains("Δv")).cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!("no Δv line in rocket product: {:?}", product.lines));
    };
    // Line looks like "Δv 9603 m/s · max-Q … · peak … g".
    let Some(reported) = parse_keyed(&line, "Δv ") else {
        return CheckOutcome::fail(format!("could not parse Δv from: {line}"));
    };
    if !reported.is_finite() || (reported - expected).abs() > 3.0 {
        return CheckOutcome::fail(format!(
            "LV-1 ideal Δv: workbench {reported:.0} vs Tsiolkovsky {expected:.0} m/s (>3 m/s)"
        ));
    }
    CheckOutcome::pass(format!(
        "LV-1 ideal Δv={reported:.0} m/s (Tsiolkovsky 2-stage {expected:.0})"
    ))
}

/// **engine** — characteristic velocity `c* = √(R̄/M · T_c) / Γ(γ)` is exact
/// 1-D chamber/nozzle theory (Γ the Vandenkerckhove choked-throat function).
/// Independently recomputes c* from the default kerolox design (`T_c = 3500 K`,
/// `γ = 1.2`, `M = 22 g/mol`) and asserts the engine product's reported `c*`
/// matches (tol 2 m/s, readout rounds to whole m/s). Drives the real
/// `engine_product()` builder (`analyze`).
fn check_engine(_app: &mut ValenxApp) -> CheckOutcome {
    // Default EngineDesign (see engine_workbench::Default).
    let (tc, gamma, molar_mass_g) = (3_500.0_f64, 1.2_f64, 22.0_f64);
    let r_universal = 8_314.462_618_f64; // J/(kmol·K) form (g/mol → kg/kmol) ⇒ R̄/M with M in g/mol.
    let r_specific = r_universal / molar_mass_g;
    let gam = gamma.sqrt() * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let expected = (r_specific * tc).sqrt() / gam;

    let product = crate::engine_workbench::engine_product();
    let line = product.lines.iter().find(|l| l.contains("c*")).cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!("no c* line in engine product: {:?}", product.lines));
    };
    let Some(reported) = parse_after_colon(&line, "c*") else {
        return CheckOutcome::fail(format!("could not parse c* from: {line}"));
    };
    if !reported.is_finite() || (reported - expected).abs() > 2.0 {
        return CheckOutcome::fail(format!(
            "engine c*: workbench {reported:.0} vs Vandenkerckhove {expected:.0} m/s (>2 m/s)"
        ));
    }
    CheckOutcome::pass(format!(
        "kerolox c*={reported:.0} m/s (1-D nozzle theory {expected:.0})"
    ))
}

/// **rotor** — momentum-theory hover: ideal induced power `P_ideal = T·v_h` with
/// `v_h = √(T/(2ρA))`, and the figure of merit `FM = P_ideal/P_actual` (≤ 1 — a
/// real rotor can't beat the ideal). Solves the default BEMT rotor **in hover**
/// (freestream V = 0) and asserts the solver's reported `figure_of_merit` equals
/// the momentum-theory `T^1.5/√(2ρA)/P` recomputed from its own
/// thrust/power/disk-area (tol 1e-6), and that FM ∈ (0, 1].
fn check_rotor(app: &mut ValenxApp) -> CheckOutcome {
    // Hover: the FM is only defined at V = 0 (see valenx_rotor::bemt).
    app.rotor.freestream_v = 0.0;
    let perf = match app.rotor.solve() {
        Ok(p) => p,
        Err(e) => return CheckOutcome::fail(format!("rotor hover solve failed: {e}")),
    };
    let (t, p, a, rho, fm) = (
        perf.thrust_n,
        perf.power_w,
        perf.disk_area_m2,
        perf.air_density,
        perf.figure_of_merit,
    );
    if !(t > 0.0 && p > 0.0 && a > 0.0 && t.is_finite() && p.is_finite()) {
        return CheckOutcome::fail(format!("non-physical hover solve T={t} P={p} A={a}"));
    }
    let ideal = t.powf(1.5) / (2.0 * rho * a).sqrt();
    let fm_expected = ideal / p;
    if (fm - fm_expected).abs() >= 1e-6 {
        return CheckOutcome::fail(format!(
            "FM mismatch: solver {fm:.6} vs momentum-theory {fm_expected:.6}"
        ));
    }
    if !(fm > 0.0 && fm <= 1.0) {
        return CheckOutcome::fail(format!("hover FM={fm:.4} not in (0,1] (beats ideal?)"));
    }
    CheckOutcome::pass(format!(
        "hover FM={fm:.4}=P_ideal/P (momentum theory, T={t:.3} N)"
    ))
}

/// **uas** — multirotor hover disk-loading momentum theory: ideal hover power
/// `P_ideal = W^1.5/√(2ρA)` and actual shaft power `P = P_ideal/FM`. Runs the
/// real `app.uas.run()` performance pipeline, reads the multirotor performance
/// rows, and asserts (a) the disk-loading identity `disk_loading = m·g₀/A` and
/// (b) the hover power equals the momentum-theory `W^1.5/√(2ρA)/FM` recomputed
/// from the all-up mass, disk area, the crate's sea-level ρ and the default
/// FM = 0.70 (tol 0.5 %).
fn check_uas(app: &mut ValenxApp) -> CheckOutcome {
    use valenx_uas::SEA_LEVEL_AIR_DENSITY;
    const G0: f64 = 9.806_65;
    let fm = app.uas.params.figure_of_merit; // default 0.70
    let result = match app.uas.run() {
        Ok(r) => r,
        Err(e) => return CheckOutcome::fail(format!("uas pipeline failed: {e}")),
    };
    // Pull the typed numbers back out of the multirotor performance rows.
    let row = |label: &str| -> Option<f64> {
        result
            .perf_rows
            .iter()
            .find(|r| r.label == label)
            .and_then(|r| parse_keyed(&r.value, ""))
    };
    let (Some(mass), Some(area), Some(dl), Some(hover_p)) = (
        row("all-up mass"),
        row("disk area"),
        row("disk loading"),
        row("hover power"),
    ) else {
        return CheckOutcome::fail(format!(
            "missing multirotor perf rows: {:?}",
            result
                .perf_rows
                .iter()
                .map(|r| &r.label)
                .collect::<Vec<_>>()
        ));
    };
    if !(mass > 0.0 && area > 0.0 && hover_p > 0.0) {
        return CheckOutcome::fail(format!("non-physical perf m={mass} A={area} P={hover_p}"));
    }
    let weight = mass * G0;
    // Disk-loading identity W/A.
    let dl_expected = weight / area;
    if (dl - dl_expected).abs() > 0.005 * dl_expected {
        return CheckOutcome::fail(format!("disk loading {dl:.2} vs W/A {dl_expected:.2} N/m²"));
    }
    // Momentum-theory hover power: P_actual = W^1.5/√(2ρA) / FM.
    let ideal = weight.powf(1.5) / (2.0 * SEA_LEVEL_AIR_DENSITY * area).sqrt();
    let p_expected = ideal / fm;
    if (hover_p - p_expected).abs() > 0.005 * p_expected {
        return CheckOutcome::fail(format!(
            "hover power {hover_p:.1} vs momentum theory {p_expected:.1} W (FM={fm:.2})"
        ));
    }
    CheckOutcome::pass(format!(
        "hover P={hover_p:.0} W = W^1.5/√(2ρA)/FM (FM={fm:.2}); DL={dl:.1} N/m²"
    ))
}

// ---------------------------------------------------------------------------
// Machine-design deep checks — each asserts a closed-form result against the
// EXACT quantity the workbench computes (not just "runs + sane").
// ---------------------------------------------------------------------------

/// **springs** — helical compression spring rate `k = G·d⁴/(8·D³·n_a)`.
/// Recomputes k from the panel's own default inputs (d = 1 mm, D = 10 mm,
/// n_a = 8, G = 79 300 MPa ⇒ k = 1.2391 N/mm) and asserts the workbench's
/// reported `stiffness k` matches (tol 1e-3 N/mm). Drives the SAME `run_springs`
/// the panel's Analyze button fires.
fn check_springs(app: &mut ValenxApp) -> CheckOutcome {
    // Default geometry/material (see springs_workbench::Default). Units: mm, MPa
    // ⇒ N/mm directly (MPa·mm⁴/mm³ = N/mm).
    let (d, dia, n, g): (f64, f64, f64, f64) = (1.0, 10.0, 8.0, 79_300.0);
    let expected = g * d.powi(4) / (8.0 * dia.powi(3) * n);
    crate::springs_workbench::run(app);
    let Some(r) = app.springs.agent_readout() else {
        return CheckOutcome::fail("no readout after analyze");
    };
    let Some(k) = parse_after_colon(&r, "stiffness k") else {
        return CheckOutcome::fail(format!("could not parse stiffness k from: {r}"));
    };
    if !k.is_finite() || (k - expected).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "spring rate k={k:.4} vs G·d⁴/(8·D³·n)={expected:.4} N/mm"
        ));
    }
    CheckOutcome::pass(format!("k={k:.4} N/mm = G·d⁴/(8·D³·n) ({expected:.4})"))
}

/// **gears** — spur geometry: pitch diameter `d = m·z/cosβ` and meshing ratio
/// `= z_mate/z`. Recomputes from the panel's default inputs (spur m = 1 mm,
/// z = 20, β = 0, z_mate = 40 ⇒ d = 20 mm, ratio = 2.0) and asserts the
/// workbench's reported `pitch diameter` and `gear ratio` match (tol 1e-3).
/// Drives the SAME `run_gears` the panel's Analyze button fires.
fn check_gears(app: &mut ValenxApp) -> CheckOutcome {
    // Default: spur (β = 0), m = 1, z = 20, mate = 40 (see gears_workbench::Default).
    let (module, z, beta_deg, z_mate): (f64, f64, f64, f64) = (1.0, 20.0, 0.0, 40.0);
    let d_expected = module * z / beta_deg.to_radians().cos();
    let ratio_expected = z_mate / z;
    crate::gears_workbench::run(app);
    let Some(r) = app.gears.agent_readout() else {
        return CheckOutcome::fail("no readout after analyze");
    };
    let pitch = parse_after_colon(&r, "pitch diameter");
    let ratio = parse_after_colon(&r, "gear ratio");
    let (Some(pitch), Some(ratio)) = (pitch, ratio) else {
        return CheckOutcome::fail(format!(
            "could not parse pitch diameter / gear ratio from: {r}"
        ));
    };
    if (pitch - d_expected).abs() >= 1e-3 || (ratio - ratio_expected).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "gears expected d={d_expected:.3} ratio={ratio_expected:.3}, got {pitch:.3}/{ratio:.3}"
        ));
    }
    CheckOutcome::pass(format!(
        "pitch d={pitch:.3} mm (m·z/cosβ) · ratio={ratio:.3} (z2/z1)"
    ))
}

/// **fasteners** — ISO metric bolt: pitch (flank) diameter `d2 = d − 0.6495·P`
/// and the ISO-898 tensile stress area `A_t`. For the default **M6** coarse
/// (P = 1.0): `d2 = 6 − 0.6495 = 5.3505 mm` and the standard `A_t = 20.12 mm²`.
/// Asserts the workbench's reported `pitch diameter` (tol 1e-3 mm) and
/// `tensile stress area` (tol 0.1 mm² of the published value). Drives the SAME
/// `run_fasteners` the panel's Compute button fires.
fn check_fasteners(app: &mut ValenxApp) -> CheckOutcome {
    // M6 coarse, pitch P = 1.0 mm: d2 = d − 0.6495·P (ISO 724).
    let d2_expected = 6.0 - 0.6495 * 1.0;
    let at_published = 20.12_f64; // ISO 898 tensile-stress-area for M6.
    crate::fasteners_workbench::run(app);
    let Some(r) = app.fasteners.agent_readout() else {
        return CheckOutcome::fail("no readout after compute");
    };
    let pitch = parse_after_colon(&r, "pitch diameter");
    let at = parse_after_colon(&r, "tensile stress area");
    let (Some(pitch), Some(at)) = (pitch, at) else {
        return CheckOutcome::fail(format!("could not parse pitch dia / stress area from: {r}"));
    };
    if (pitch - d2_expected).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "M6 pitch dia {pitch:.4} vs d−0.6495·P {d2_expected:.4} mm"
        ));
    }
    if !at.is_finite() || (at - at_published).abs() >= 0.1 {
        return CheckOutcome::fail(format!("M6 A_t {at:.3} vs ISO-898 {at_published:.2} mm²"));
    }
    CheckOutcome::pass(format!(
        "M6 d2={pitch:.4} mm (d−0.6495P) · A_t={at:.2} mm² (ISO 898)"
    ))
}

/// **frames** — cross-section area of the default **IPE 200 I-beam** is the
/// exact built-up sum `A = 2·b·t_f + t_w·(h − 2·t_f)`. (The panel reports area +
/// perimeter — not I/S — so the area is the closed-form quantity to assert.)
/// For h = 200, b = 100, t_w = 5.6, t_f = 8.5 ⇒ `A = 2724.8 mm²`. Asserts the
/// reported `cross-section area` (tol 1e-2 mm²). Drives the SAME `run_frames`
/// the panel's Compute button fires.
///
/// note: moment of inertia / section modulus are NOT exposed by the frames
/// panel (it computes area + perimeter only), so the deep check asserts the
/// area's closed form rather than `I = b·h³/12`.
fn check_frames(app: &mut ValenxApp) -> CheckOutcome {
    // Default IPE 200 I-beam (see frames_workbench::Default).
    let (h, b, tw, tf): (f64, f64, f64, f64) = (200.0, 100.0, 5.6, 8.5);
    let area_expected = 2.0 * b * tf + tw * (h - 2.0 * tf);
    crate::frames_workbench::run(app);
    let Some(r) = app.frames.agent_readout() else {
        return CheckOutcome::fail("no readout after compute");
    };
    let Some(area) = parse_after_colon(&r, "cross-section area") else {
        return CheckOutcome::fail(format!("could not parse cross-section area from: {r}"));
    };
    if !area.is_finite() || (area - area_expected).abs() >= 1e-2 {
        return CheckOutcome::fail(format!(
            "I-beam area {area:.2} vs 2·b·tf+tw·(h−2tf) {area_expected:.2} mm²"
        ));
    }
    CheckOutcome::pass(format!(
        "IPE200 area={area:.2} mm² = 2·b·tf+tw·(h−2tf) ({area_expected:.2})"
    ))
}

/// **collision** — AABB overlap predicate, asserted both ways. The default pair
/// (A = [0,0,0]→[10,20,30], B = [20,0,0]→[30,10,10]) is disjoint with an exact
/// L2 gap of 10 along x ⇒ `intersect: no`, `separation = 10`. Then an
/// overlapping pair (B moved to [5,5,5]→[15,15,15] via the agent control) ⇒
/// `intersect: yes`, `separation = 0`. Drives the SAME `run_collision` the
/// panel's Compute button fires.
fn check_collision(app: &mut ValenxApp) -> CheckOutcome {
    use crate::agent_commands::AgentValue;
    // (1) Default disjoint pair.
    crate::collision_workbench::run(app);
    let Some(r1) = app.collision.agent_readout() else {
        return CheckOutcome::fail("no readout (disjoint)");
    };
    // The verdict line is `intersect  : no — disjoint`; note "overlapping" also
    // appears in the separation line's annotation, so match the verdict, not a
    // bare "overlap" substring.
    let disjoint_ok = r1.contains("disjoint") && r1.contains("intersect");
    let Some(sep) = parse_after_colon(&r1, "separation") else {
        return CheckOutcome::fail(format!("could not parse separation from: {r1}"));
    };
    if !disjoint_ok || (sep - 10.0).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "default pair: expected disjoint, gap 10; got sep={sep:.4} / {r1}"
        ));
    }

    // (2) Overlapping pair: move Box B to straddle A.
    for (cap, v) in [
        ("Box B min x", 5.0_f64),
        ("Box B min y", 5.0),
        ("Box B min z", 5.0),
        ("Box B max x", 15.0),
        ("Box B max y", 15.0),
        ("Box B max z", 15.0),
    ] {
        if let Err(e) = app.collision.agent_set(cap, &AgentValue::Float(v)) {
            return CheckOutcome::fail(format!("set {cap} failed: {e}"));
        }
    }
    crate::collision_workbench::run(app);
    let Some(r2) = app.collision.agent_readout() else {
        return CheckOutcome::fail("no readout (overlap)");
    };
    let Some(sep2) = parse_after_colon(&r2, "separation") else {
        return CheckOutcome::fail(format!("could not parse separation from: {r2}"));
    };
    // Verdict is `yes — boxes overlap`; match that, not the annotation's
    // "overlapping" which is present in both cases.
    if !r2.contains("boxes overlap") || sep2.abs() >= 1e-9 {
        return CheckOutcome::fail(format!(
            "overlap pair: expected overlap, sep 0; got sep={sep2:.4} / {r2}"
        ));
    }
    CheckOutcome::pass("AABB: disjoint gap=10 (no), overlap sep=0 (yes)".to_string())
}

// ---------------------------------------------------------------------------
// Simulation deep checks (batch 1) — each asserts a closed-form / conservation
// result against the EXACT quantity the workbench computes.
// ---------------------------------------------------------------------------

/// **fem** — Euler–Bernoulli cantilever tip deflection `δ = P·L³/(3·E·I)`. The
/// canonical FEM product is a steel cantilever (L = 1 m, section 50×100 mm,
/// E = 200 GPa, P = 5 kN; `I = b·h³/12 = 4.1667e-6 m⁴` ⇒ `δ_EB = 2.000 mm`).
/// Two assertions:
///   1. the **analytic** PL³/3EI the panel prints equals the independently
///      recomputed `δ_EB = 2.000 mm` exactly (tol 1e-3 mm) — the exact known
///      value, validating the closed-form reference; and
///   2. the **FE** tip deflection from the real `valenx_fem` 3-D linear-static
///      solver sits in the physically-correct band relative to it.
///
/// note: the FE mesh is coarse linear tetrahedra (Tet4, 16×4×2), which **shear-
/// lock** in bending and so report a *stiffer* (smaller) deflection than the
/// shear-free Euler–Bernoulli value — here ≈ 0.6·δ_EB. That is a documented
/// discretization signature, not an error, so the FE result is asserted to be in
/// `(0.4·δ_EB, δ_EB)` (stiffer than, but the same order as, beam theory) rather
/// than matched exactly. The panel prints both numbers side by side.
fn check_fem(_app: &mut ValenxApp) -> CheckOutcome {
    // Closed-form Euler–Bernoulli δ for the fem_product cantilever (constants
    // mirror fem_workbench::fem_product: LX=1, LY=0.10, LZ=0.05, E=200 GPa, P=5 kN).
    let (lx, ly, lz, e, p): (f64, f64, f64, f64, f64) = (1.0, 0.10, 0.05, 200.0e9, 5000.0);
    let i_area = lz * ly.powi(3) / 12.0;
    let delta_eb_mm = p * lx.powi(3) / (3.0 * e * i_area) * 1.0e3;

    let product = crate::fem_workbench::fem_product();
    // The card line is "max deflection (FE): X mm   vs analytical PL³/3EI = Y mm".
    let line = product
        .lines
        .iter()
        .find(|l| l.contains("max deflection (FE)"))
        .cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!(
            "no FE deflection line in fem product: {:?}",
            product.lines
        ));
    };
    let fe_mm = parse_after_colon(&line, "max deflection (FE)");
    let analytic_mm = parse_keyed(&line, "PL³/3EI = ");
    let (Some(fe_mm), Some(analytic_mm)) = (fe_mm, analytic_mm) else {
        return CheckOutcome::fail(format!(
            "could not parse FE / analytic deflection from: {line}"
        ));
    };
    // (1) The panel's analytic PL³/3EI must equal the independent closed form.
    if (analytic_mm - delta_eb_mm).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "panel analytic δ={analytic_mm:.3} vs PL³/3EI {delta_eb_mm:.3} mm"
        ));
    }
    // (2) The FE result must be stiffer than beam theory (Tet4 locking) but the
    // same order — physically ordered, not a wild number.
    if !fe_mm.is_finite() || fe_mm <= 0.4 * delta_eb_mm || fe_mm >= delta_eb_mm {
        return CheckOutcome::fail(format!(
            "FE tip δ={fe_mm:.3} mm not in expected locked band (0.4·δ_EB, δ_EB)=({:.3},{:.3})",
            0.4 * delta_eb_mm,
            delta_eb_mm
        ));
    }
    CheckOutcome::pass(format!(
        "PL³/3EI={analytic_mm:.3} mm (exact); FE δ={fe_mm:.3} mm (Tet4 locked, 0.4–1×)"
    ))
}

/// **fields** — descriptive statistics are exact. Sets a known dataset
/// (`1 2 3 4 5`) via the agent control, runs the compute, and asserts the
/// reported mean (= 3), population std dev (= √2 ≈ 1.41421), and min/max (1/5)
/// match exactly (tol 1e-4). Drives the SAME `run_fields` the Compute button
/// fires.
fn check_fields(app: &mut ValenxApp) -> CheckOutcome {
    use crate::agent_commands::AgentValue;
    if let Err(e) = app
        .fields
        .agent_set("numbers", &AgentValue::Str("1 2 3 4 5".into()))
    {
        return CheckOutcome::fail(format!("set numbers failed: {e}"));
    }
    crate::fields_workbench::run(app);
    let Some(r) = app.fields.agent_readout() else {
        return CheckOutcome::fail("no readout after compute");
    };
    let mean = parse_after_colon(&r, "mean");
    let std = parse_after_colon(&r, "std dev");
    let (Some(mean), Some(std)) = (mean, std) else {
        return CheckOutcome::fail(format!("could not parse mean / std dev from: {r}"));
    };
    let std_expected = 2.0_f64.sqrt(); // population σ of 1..5 = √2.
    if (mean - 3.0).abs() >= 1e-4 || (std - std_expected).abs() >= 1e-4 {
        return CheckOutcome::fail(format!(
            "stats of 1..5: mean {mean:.5} (want 3), std {std:.5} (want {std_expected:.5})"
        ));
    }
    CheckOutcome::pass(format!("mean={mean:.4} σ_pop={std:.5} (=√2) for 1,2,3,4,5"))
}

/// **nodegraph** — the default graph `Constant(2) + Constant(2) → Add → Output`
/// evaluates to exactly `4`. Drives the real topological evaluator and asserts
/// the Output node's value is 4.0 (tol 1e-6). Drives the SAME `run_eval` the
/// Evaluate button fires.
fn check_nodegraph(app: &mut ValenxApp) -> CheckOutcome {
    crate::nodegraph_workbench::run(app);
    let Some(r) = app.nodegraph.agent_readout() else {
        return CheckOutcome::fail("no readout after evaluate");
    };
    // The readout ends "outputs: #3=4.0000"; pull the value after the first "=".
    let Some(after) = r.split("outputs:").nth(1) else {
        return CheckOutcome::fail(format!("no outputs in readout: {r}"));
    };
    let Some(eq) = after.find('=') else {
        return CheckOutcome::fail(format!("no output value in: {r}"));
    };
    let Some(out) = parse_keyed(&after[eq..], "=") else {
        return CheckOutcome::fail(format!("could not parse output value from: {r}"));
    };
    if !out.is_finite() || (out - 4.0).abs() >= 1e-6 {
        return CheckOutcome::fail(format!("Output={out} expected 2+2=4"));
    }
    CheckOutcome::pass(format!("Output node = {out:.4} (Constant 2 + Constant 2)"))
}

/// **reactdyn** — NVE (microcanonical) energy conservation is the integrator's
/// correctness signal. Runs a short H2 AIMD trajectory (the default RHF/STO-3G
/// velocity-Verlet, steps trimmed for speed) **synchronously** and asserts the
/// relative total-energy drift `max|E(t) − E₀|/|E₀|` stays below 1e-2 — i.e. the
/// velocity-Verlet + analytic forces conserve energy.
fn check_reactdyn(app: &mut ValenxApp) -> CheckOutcome {
    // Default is AIMD / H2 / NVE; trim the step count so the QM MD is fast.
    app.reactdyn = Default::default();
    if let Err(e) = app
        .reactdyn
        .agent_set("steps", &crate::agent_commands::AgentValue::Int(20))
    {
        return CheckOutcome::fail(format!("set steps failed: {e}"));
    }
    crate::reactdyn_workbench::run(app);
    // A failed run leaves no trajectory ⇒ no drift; the status carries the reason.
    let Some(drift) = app.reactdyn.last_energy_rel_drift() else {
        return CheckOutcome::fail(format!(
            "AIMD produced no trajectory ({})",
            app.reactdyn.status_line()
        ));
    };
    if !drift.is_finite() || drift >= 1.0e-2 {
        return CheckOutcome::fail(format!("NVE energy drift |ΔE/E₀|={drift:.2e} (≥1e-2)"));
    }
    CheckOutcome::pass(format!(
        "H2 AIMD NVE: |ΔE/E₀|={drift:.2e} (<1e-2, energy conserved)"
    ))
}

/// **mbd** — the default articulated-rod pendulum (L = 1 m, small release 0.2
/// rad) swings at the small-angle period `T = 2π√(L/g)`. Runs the real
/// `valenx_mbd` integrator, measures the period from the angle-trace zero
/// crossings, and asserts it matches `2π√(L/g) = 2.006 s` within 3 % (covers the
/// finite-amplitude correction + sample resolution); also asserts the NVE energy
/// drift is small. Drives the SAME `MbdWorkbenchState::run` the Run button fires.
fn check_mbd(app: &mut ValenxApp) -> CheckOutcome {
    let result = match app.mbd.run() {
        Ok(r) => r,
        Err(e) => return CheckOutcome::fail(format!("mbd run failed: {e}")),
    };
    // Small-angle period of the default 1 m pendulum at g = 9.81.
    let (l, g) = (app.mbd.params.rod_length, app.mbd.params.gravity);
    let t_analytic = std::f64::consts::TAU * (l / g).sqrt();

    // Measure the period from the angle trace: the bob starts at +θ₀, so the
    // gaps between consecutive sign-changes of `angle` are each a half-period.
    let s = &result.samples;
    let mut crossings: Vec<f64> = Vec::new();
    for w in s.windows(2) {
        if w[0].angle == 0.0 {
            continue;
        }
        if w[0].angle.signum() != w[1].angle.signum() {
            // Linear interpolation to the zero crossing time.
            let (t0, a0, t1, a1) = (w[0].t, w[0].angle, w[1].t, w[1].angle);
            let frac = a0 / (a0 - a1);
            crossings.push(t0 + frac * (t1 - t0));
        }
    }
    if crossings.len() < 2 {
        return CheckOutcome::fail(format!(
            "too few zero-crossings ({}) to measure period",
            crossings.len()
        ));
    }
    let t_measured = 2.0 * (crossings[1] - crossings[0]);
    if !t_measured.is_finite() || (t_measured - t_analytic).abs() / t_analytic > 0.03 {
        return CheckOutcome::fail(format!(
            "pendulum period {t_measured:.4} s vs 2π√(L/g) {t_analytic:.4} s (>3%)"
        ));
    }
    // NVE conservation: the undamped pendulum should hold energy.
    if !result.energy_rel_drift.is_finite() || result.energy_rel_drift >= 1.0e-2 {
        return CheckOutcome::fail(format!(
            "pendulum energy drift |ΔE/E₀|={:.2e} (≥1e-2)",
            result.energy_rel_drift
        ));
    }
    CheckOutcome::pass(format!(
        "pendulum T={t_measured:.3} s ≈ 2π√(L/g)={t_analytic:.3} s; drift {:.1e}",
        result.energy_rel_drift
    ))
}

// ---------------------------------------------------------------------------
// Life-sciences deep checks — each asserts an exact / invariant result against
// the EXACT quantity the workbench computes.
// ---------------------------------------------------------------------------

/// **genetics** — sequence statistics are exact. Feeds a known DNA sequence
/// (`ATGC`) to the Sequence panel's real `valenx_bioseq` Analyze and asserts the
/// reported length (= 4) and GC content (= 50.00 %, since G+C = 2 of 4). Drives
/// the SAME `run_analyze` the panel's Analyze button fires.
fn check_genetics(app: &mut ValenxApp) -> CheckOutcome {
    use valenx_bioseq::SeqKind;
    // ATGC: one each of A,T,G,C ⇒ length 4, GC = (G+C)/4 = 50 %.
    let readout = app.genetics.sequence.agent_analyze("ATGC", SeqKind::Dna);
    let len = parse_after_colon(&readout, "Length");
    let gc = parse_after_colon(&readout, "GC content");
    let (Some(len), Some(gc)) = (len, gc) else {
        return CheckOutcome::fail(format!(
            "could not parse Length / GC content from: {readout}"
        ));
    };
    if (len - 4.0).abs() >= 1e-6 || (gc - 50.0).abs() >= 1e-4 {
        return CheckOutcome::fail(format!(
            "ATGC: length {len} (want 4), GC {gc:.4}% (want 50)"
        ));
    }
    CheckOutcome::pass(format!("ATGC: length={len:.0}, GC={gc:.2}% (G+C=2/4)"))
}

/// **neuro** — unmyelinated axon conduction velocity is cable theory's
/// `v = k·√d`. For the default fiber (d = 1 µm, k = 2.0 m/s per √µm) the velocity
/// is `2.0·√1 = 2.0 m/s`. Recomputes from the panel's own inputs via the SAME
/// `valenx_neuro::unmyelinated_conduction_velocity` the panel calls and asserts
/// the match (tol 1e-9); also checks the myelinated Hursh `v = k·d`.
fn check_neuro(app: &mut ValenxApp) -> CheckOutcome {
    let (v_myel, v_unmyel) = app.neuro.conduction_velocities();
    // Default state: d = 1.0 µm, k_unmyel = 2.0, k_myel = HURSH (6.0).
    let d: f64 = 1.0;
    let (k_unmyel, k_myel): (f64, f64) = (2.0, valenx_neuro::HURSH_FACTOR_M_PER_S_PER_UM);
    let v_unmyel_expected = k_unmyel * d.sqrt();
    let v_myel_expected = k_myel * d;
    if !v_unmyel.is_finite() || (v_unmyel - v_unmyel_expected).abs() >= 1e-9 {
        return CheckOutcome::fail(format!(
            "unmyelinated v={v_unmyel} vs k·√d={v_unmyel_expected} m/s"
        ));
    }
    if (v_myel - v_myel_expected).abs() >= 1e-9 {
        return CheckOutcome::fail(format!(
            "myelinated v={v_myel} vs k·d={v_myel_expected} m/s"
        ));
    }
    CheckOutcome::pass(format!(
        "unmyel v={v_unmyel:.3} m/s = k·√d; myel v={v_myel:.1} = k·d (Hursh)"
    ))
}

/// **variant** — the HGVS parser is exact. Sets `p.R273H` via the agent control,
/// runs the real `valenx_variant_effect` parse, and asserts the readout reports
/// the protein substitution `R273H` at residue 273 (ref R, pos 273, alt H).
fn check_variant(app: &mut ValenxApp) -> CheckOutcome {
    use crate::agent_commands::AgentValue;
    if let Err(e) = app.variant_effect.agent_set(
        "Variants (one per line)",
        &AgentValue::Str("p.R273H".into()),
    ) {
        return CheckOutcome::fail(format!("set variant failed: {e}"));
    }
    crate::variant_effect_workbench::run(app);
    let Some(r) = app.variant_effect.agent_readout() else {
        return CheckOutcome::fail("no readout after parse");
    };
    // The describe() line reads "protein substitution R273H — residue 273".
    if r.contains("error:") {
        return CheckOutcome::fail(format!("parse error: {r}"));
    }
    if !(r.contains("R273H") && r.contains("273")) {
        return CheckOutcome::fail(format!("expected R273H / residue 273, got: {r}"));
    }
    CheckOutcome::pass("HGVS p.R273H → ProteinSub R(273)H".to_string())
}

/// **ppi** — protein-interaction shortest path is an exact BFS graph result. The
/// default demo network plants a **direct** GUARD→EFF-A edge, so the shortest
/// path between node 0 (GUARD) and node 3 (EFF-A, the first pathogen) is exactly
/// one hop. Sets the analysis to ShortestPath, runs the real `valenx_ppi`
/// screen + BFS, and asserts `path = [0, 3]` (1 hop, GUARD→EFF-A).
fn check_ppi(app: &mut ValenxApp) -> CheckOutcome {
    use crate::ppi_workbench::PpiAnalysis;
    // Default endpoints are path_from=0 (GUARD), path_to=3 (EFF-A); select the
    // shortest-path analysis (default is degree centrality).
    app.ppi.params.analysis = PpiAnalysis::ShortestPath;
    let result = match app.ppi.run() {
        Ok(r) => r,
        Err(e) => return CheckOutcome::fail(format!("ppi run failed: {e}")),
    };
    let Some(hops) = result.path_hops() else {
        return CheckOutcome::fail("no shortest path (GUARD↔EFF-A should be connected)");
    };
    let path = result.path.clone().unwrap_or_default();
    if hops != 1 || path.first() != Some(&0) || path.last() != Some(&3) {
        return CheckOutcome::fail(format!(
            "expected 1-hop path [0,3], got {hops} hops {path:?}"
        ));
    }
    // Confirm the endpoint labels are GUARD and EFF-A.
    let from = result.nodes.first().map(|n| n.name.as_str()).unwrap_or("?");
    let to = result.nodes.get(3).map(|n| n.name.as_str()).unwrap_or("?");
    if from != "GUARD" || to != "EFF-A" {
        return CheckOutcome::fail(format!("endpoints {from}→{to} (want GUARD→EFF-A)"));
    }
    CheckOutcome::pass(format!(
        "shortest path {from}→{to} = {hops} hop (direct edge)"
    ))
}

/// **morphogenesis** — Gray–Scott reaction–diffusion has no clean scalar ground
/// truth (it's an emergent pattern), so this asserts the solver's documented
/// **bounded-stability invariant**: the `U`/`V` morphogen concentrations stay
/// finite and in `[0, 1]` every step, and the field actually evolves from its
/// seeded germ (mean V changes). Steps the real field and checks the invariant.
fn check_morphogenesis(app: &mut ValenxApp) -> CheckOutcome {
    let mean_v0 = app.morphogenesis.field.mean_v();
    if mean_v0 <= 0.0 || !mean_v0.is_finite() {
        return CheckOutcome::fail(format!("seeded field has no germ (mean V = {mean_v0})"));
    }
    // Advance several frames of the real explicit-Euler Gray–Scott solver.
    for _ in 0..10 {
        crate::morphogenesis_workbench::run(app);
    }
    let f = &app.morphogenesis.field;
    // Every U and V must be finite and clamped to [0, 1] (the documented guard).
    let bad_u =
        f.u.iter()
            .find(|&&x| !x.is_finite() || !(0.0..=1.0).contains(&x));
    let bad_v =
        f.v.iter()
            .find(|&&x| !x.is_finite() || !(0.0..=1.0).contains(&x));
    if let Some(b) = bad_u.or(bad_v) {
        return CheckOutcome::fail(format!("U/V out of [0,1] or non-finite: {b}"));
    }
    let (vmin, vmax) = f.field_minmax();
    let mean_v1 = f.mean_v();
    // The reaction must have run (the field changed from the static seed).
    if (mean_v1 - mean_v0).abs() < 1e-9 {
        return CheckOutcome::fail(format!("field did not evolve (mean V stayed {mean_v0})"));
    }
    CheckOutcome::pass(format!(
        "Gray–Scott U,V∈[0,1] (V {vmin:.3}..{vmax:.3}); evolved {} steps",
        f.steps
    ))
}

// ---------------------------------------------------------------------------
// Civil & AEC deep checks — each asserts an exact / standard value against the
// quantity the workbench computes (via its deterministic *_product builder).
// ---------------------------------------------------------------------------

/// **geomatics** — the great-circle (haversine) distance is exact spherical
/// trigonometry. The default product is the canonical Cambridge (52.205, 0.119)
/// → Paris (48.857, 2.351) worked example ≈ **404.3 km**. Drives the real
/// `geomatics_product()` (which calls `valenx_geomatics::haversine_distance`,
/// `R = 6 371 008.8 m`) and asserts the reported great-circle matches the
/// documented value (tol 0.5 km).
fn check_geomatics(_app: &mut ValenxApp) -> CheckOutcome {
    let product = crate::geomatics_workbench::geomatics_product();
    let line = product
        .lines
        .iter()
        .find(|l| l.contains("great-circle"))
        .cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!(
            "no great-circle line in geomatics: {:?}",
            product.lines
        ));
    };
    let Some(km) = parse_after_colon(&line, "great-circle") else {
        return CheckOutcome::fail(format!("could not parse great-circle km from: {line}"));
    };
    if !km.is_finite() || (km - 404.3).abs() > 0.5 {
        return CheckOutcome::fail(format!(
            "Cambridge→Paris great-circle {km:.3} km vs haversine 404.3 km (>0.5)"
        ));
    }
    CheckOutcome::pass(format!(
        "Cambridge→Paris great-circle = {km:.2} km (haversine)"
    ))
}

/// **piping** — ASME B36.10 pipe dimensions are exact. The default product is
/// NPS 2 Sch 40, whose outer diameter is `2.375 in × 25.4 = 60.325 mm` and whose
/// bore (`OD − 2·wall`, wall 3.91 mm) is ≈ 52.5 mm. Drives the real
/// `piping_product()` and asserts the reported OD = 60.325 mm (tol 1e-3) and the
/// ID is in the ASME bore range (tol 0.1 mm of 52.5).
fn check_piping(_app: &mut ValenxApp) -> CheckOutcome {
    let product = crate::piping_workbench::piping_product();
    let od_line = product.lines.iter().find(|l| l.contains("outer diameter"));
    let id_line = product.lines.iter().find(|l| l.contains("inner diameter"));
    let (Some(od_line), Some(id_line)) = (od_line, id_line) else {
        return CheckOutcome::fail(format!("no OD/ID lines in piping: {:?}", product.lines));
    };
    let od = parse_after_colon(od_line, "outer diameter");
    let id = parse_after_colon(id_line, "inner diameter");
    let (Some(od), Some(id)) = (od, id) else {
        return CheckOutcome::fail("could not parse OD/ID".to_string());
    };
    // NPS 2 OD = 2.375 in × 25.4 (exact ASME B36.10).
    if (od - 60.325).abs() >= 1e-3 {
        return CheckOutcome::fail(format!("NPS 2 OD {od:.3} mm vs ASME 60.325 mm"));
    }
    // NPS 2 Sch 40 wall 3.91 mm ⇒ ID ≈ 52.5 mm.
    if (id - 52.5).abs() > 0.1 || id >= od {
        return CheckOutcome::fail(format!("NPS 2 Sch40 ID {id:.3} mm vs ASME ≈52.5 mm"));
    }
    CheckOutcome::pass(format!(
        "NPS 2 Sch 40 OD={od:.3} mm, ID={id:.3} mm (ASME B36.10)"
    ))
}

/// **hvac** — Darcy–Weisbach duct pressure drop `ΔP = f·(L/D)·½·ρ·V²`. The
/// default product is a 200 mm round duct, 10 m long at 5 m/s, `f = 0.02`,
/// `ρ = 1.204 kg/m³` (air 20 °C) ⇒ ΔP = `0.02·(10/0.2)·0.5·1.204·25 = 15.05 Pa`.
/// Drives the real `hvac_product()` and asserts the reported pressure drop
/// matches the recomputed Darcy–Weisbach value (tol 0.05 Pa).
fn check_hvac(_app: &mut ValenxApp) -> CheckOutcome {
    // Recompute ΔP with the crate's own air density + Darcy–Weisbach.
    let (d, l, v, f, rho): (f64, f64, f64, f64, f64) = (0.2, 10.0, 5.0, 0.02, 1.204);
    let dp_expected = f * (l / d) * 0.5 * rho * v * v;

    let product = crate::hvac_workbench::hvac_product();
    let line = product
        .lines
        .iter()
        .find(|l| l.contains("pressure drop"))
        .cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!(
            "no pressure-drop line in hvac: {:?}",
            product.lines
        ));
    };
    let Some(dp) = parse_keyed(&line, "pressure drop ") else {
        return CheckOutcome::fail(format!("could not parse ΔP from: {line}"));
    };
    if !dp.is_finite() || (dp - dp_expected).abs() > 0.05 {
        return CheckOutcome::fail(format!(
            "duct ΔP {dp:.2} Pa vs Darcy–Weisbach {dp_expected:.2} Pa (>0.05)"
        ));
    }
    CheckOutcome::pass(format!(
        "duct ΔP={dp:.2} Pa = f·(L/D)·½ρV² ({dp_expected:.2})"
    ))
}

/// **interior** — floor-plan bookkeeping is exact. The default product is a
/// 6 m × 4 m living room with 2 furniture pieces (a sofa + a table). Drives the
/// real `interior_product()` and asserts the reported room extent (6.0 × 4.0 m)
/// and piece count (2).
fn check_interior(_app: &mut ValenxApp) -> CheckOutcome {
    let product = crate::interior_workbench::interior_product();
    let joined = product.lines.join(" • ");
    // Lines: "1 room · 2 pieces" and "6.0 × 4.0 m".
    let dims_ok = joined.contains("6.0 × 4.0 m");
    let pieces_ok = joined.contains("2 pieces");
    if !(dims_ok && pieces_ok) {
        return CheckOutcome::fail(format!("expected 6.0×4.0 m room / 2 pieces, got: {joined}"));
    }
    CheckOutcome::pass("floor plan 6.0 × 4.0 m (24 m²), 2 pieces".to_string())
}

// ---------------------------------------------------------------------------
// CAD & mesh deep checks — each asserts an exact geometric / bookkeeping value
// against the quantity the workbench computes.
// ---------------------------------------------------------------------------

/// **cad** — parametric CSG solid volume is exact. The default feature tree is a
/// unit box (size 1) cut by a Ø0.5 through-cylinder (r = 0.25), so the exact CSG
/// volume is `1³ − π·0.25²·1 = 1 − π/16 ≈ 0.80365 u³`. Drives the real
/// `cad_product()` (which uses `valenx_cad::solid_volume`, analytic) and asserts
/// the reported volume matches (tol 1e-3).
fn check_cad(_app: &mut ValenxApp) -> CheckOutcome {
    let v_expected = 1.0 - std::f64::consts::PI * 0.25_f64.powi(2) * 1.0;
    let product = crate::cad_workbench::cad_product();
    let line = product.lines.iter().find(|l| l.contains("volume")).cloned();
    let Some(line) = line else {
        return CheckOutcome::fail(format!("no volume line in cad: {:?}", product.lines));
    };
    // Line: "volume 0.8036 u³ · surface area …".
    let Some(vol) = parse_keyed(&line, "volume ") else {
        return CheckOutcome::fail(format!("could not parse volume from: {line}"));
    };
    if !vol.is_finite() || (vol - v_expected).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "CSG volume {vol:.4} vs 1−π·0.25² {v_expected:.4} u³"
        ));
    }
    CheckOutcome::pass(format!(
        "punched-box volume {vol:.4} u³ = 1−π/16 ({v_expected:.4})"
    ))
}

/// **brep** — boolean solid volume. The default is the SAME punched cube as cad
/// (unit box − Ø0.5 cylinder), so the analytic difference volume is
/// `1 − π·0.25²·1 ≈ 0.804`. Drives the real `valenx-truck-cad` boolean + tessellation
/// (via the `brep.build` path) and reads `agent_readout`'s reported volume.
///
/// note: unlike cad's analytic `solid_volume`, brep measures the **tessellated
/// mesh** volume, so the faceted cylinder makes it deviate slightly from the
/// analytic value — hence a 3 % band rather than an exact match.
fn check_brep(app: &mut ValenxApp) -> CheckOutcome {
    let v_analytic = 1.0 - std::f64::consts::PI * 0.25_f64.powi(2) * 1.0;
    crate::brep_workbench::run(app);
    let Some(r) = app.brep.agent_readout() else {
        return CheckOutcome::fail("no readout after build");
    };
    if r.contains("build failed") {
        return CheckOutcome::fail(format!("brep build failed: {r}"));
    }
    // Readout: "… · volume 0.8036 · watertight · bbox …".
    let Some(vol) = parse_keyed(&r, "volume ") else {
        return CheckOutcome::fail(format!("could not parse volume from: {r}"));
    };
    if !vol.is_finite() || (vol - v_analytic).abs() / v_analytic > 0.03 {
        return CheckOutcome::fail(format!(
            "brep volume {vol:.4} vs analytic difference {v_analytic:.4} (>3%)"
        ));
    }
    CheckOutcome::pass(format!(
        "boolean difference volume {vol:.4} ≈ 1−π/16 ({v_analytic:.4}, tessellated)"
    ))
}

/// **mesh** — AABB of the canonical valenx LV-1 mesh, the SAME extents the Mesh
/// Toolbox Inspector reports. Asserts the LV-1 mesh's measured bounding-box
/// extents (the documented ≈ 9.600 × 9.600 × 34.700) — finite, positive, and a
/// tall body (Z ≫ X ≈ Y). Drives the real `valenx_mesh` AABB computation.
fn check_mesh(_app: &mut ValenxApp) -> CheckOutcome {
    let Some((dx, dy, dz)) = crate::mesh_toolbox::lv1_canonical_aabb_extents() else {
        return CheckOutcome::fail("LV-1 mesh has no AABB (empty?)");
    };
    if !(dx.is_finite() && dy.is_finite() && dz.is_finite() && dx > 0.0 && dy > 0.0 && dz > 0.0) {
        return CheckOutcome::fail(format!("non-physical AABB {dx} × {dy} × {dz}"));
    }
    // Documented LV-1 extents ≈ 9.6 × 9.6 × 34.7; assert within 0.1.
    if (dx - 9.6).abs() > 0.1 || (dy - 9.6).abs() > 0.1 || (dz - 34.7).abs() > 0.1 {
        return CheckOutcome::fail(format!(
            "LV-1 AABB {dx:.3}×{dy:.3}×{dz:.3} vs documented 9.6×9.6×34.7"
        ));
    }
    CheckOutcome::pass(format!("LV-1 mesh AABB {dx:.3} × {dy:.3} × {dz:.3}"))
}

/// **sheetmetal** — bend allowance is exact: `BA = (π·θ/180)·(R + K·t)`. For the
/// default 90° bend (t = 1 mm, R = 1 mm, K = 0.44) the neutral-axis arc is
/// `(π/2)·(1 + 0.44) = 1.44·π/2 ≈ 2.2619 mm`. Drives the real `run_sheetmetal`
/// and asserts the reported bend allowance matches (tol 1e-3 mm).
fn check_sheetmetal(app: &mut ValenxApp) -> CheckOutcome {
    // Default: θ=90°, R=1, K=0.44, t=1 ⇒ BA = (π·90/180)·(1 + 0.44·1).
    let (theta_deg, r, k, t): (f64, f64, f64, f64) = (90.0, 1.0, 0.44, 1.0);
    let ba_expected = (std::f64::consts::PI * theta_deg / 180.0) * (r + k * t);
    crate::sheetmetal_workbench::run(app);
    let Some(readout) = app.sheetmetal.agent_readout() else {
        return CheckOutcome::fail("no readout after compute");
    };
    let Some(ba) = parse_after_colon(&readout, "bend allowance") else {
        return CheckOutcome::fail(format!("could not parse bend allowance from: {readout}"));
    };
    if !ba.is_finite() || (ba - ba_expected).abs() >= 1e-3 {
        return CheckOutcome::fail(format!(
            "bend allowance {ba:.4} vs (π·θ/180)·(R+K·t) {ba_expected:.4} mm"
        ));
    }
    CheckOutcome::pass(format!(
        "BA={ba:.4} mm = (π·θ/180)·(R+K·t) ({ba_expected:.4})"
    ))
}

/// **draft2d** — exact 2-D geometry/bookkeeping. The default demo drawing is a
/// closed 60×40 polyline + a circle + a diagonal + a text label (4 entities), so
/// its bounding extent is exactly 60 × 40 units. Drives the real
/// `draft2d_product()` and asserts the reported extent (60 × 40) and entity
/// count (4).
fn check_draft2d(_app: &mut ValenxApp) -> CheckOutcome {
    let product = crate::draft2d_workbench::draft2d_product();
    let joined = product.lines.join(" • ");
    // Lines: "4 entities", "extent 60 × 40 units".
    if !(joined.contains("4 entities") && joined.contains("extent 60 × 40 units")) {
        return CheckOutcome::fail(format!("expected 4 entities / extent 60×40, got: {joined}"));
    }
    CheckOutcome::pass("2-D drawing: 4 entities, extent 60 × 40 units".to_string())
}

/// **animate** — keyframe interpolation. The demo sweep is joint 0: 0 → π over
/// 2 s. The keyframe endpoints are exact and easing-independent (start = 0,
/// end = π), and with **Linear** easing the midpoint (t = 1 s) is exactly the
/// linear interpolation `π/2`. Asserts the product's start/end and the linear
/// midpoint (tol 1e-6).
fn check_animate(_app: &mut ValenxApp) -> CheckOutcome {
    use valenx_animate::TweenMode;
    // Endpoints from the product card (EaseInOut default): start 0.000, end 3.142.
    let product = crate::animate_workbench::animate_product();
    let joined = product.lines.join(" • ");
    if !(joined.contains("t=0.00 s  →  0.000") && joined.contains("t=2.00 s  →  3.142")) {
        return CheckOutcome::fail(format!("expected endpoints 0.000 / 3.142, got: {joined}"));
    }
    // Linear interpolation midpoint: at t=1 s of a 0→π/2s sweep, value = π/2.
    let mid_lin = crate::animate_workbench::sample_demo_joint0(TweenMode::Linear, 1.0);
    let mid_expected = std::f64::consts::FRAC_PI_2;
    if (mid_lin - mid_expected).abs() >= 1e-6 {
        return CheckOutcome::fail(format!(
            "Linear midpoint {mid_lin:.6} vs π/2 {mid_expected:.6} rad"
        ));
    }
    CheckOutcome::pass(format!(
        "keyframe sweep 0→π; Linear t=1s = {mid_lin:.5} rad (=π/2)"
    ))
}

/// **reverse** — point-cloud → mesh reconstructs a known primitive. The default
/// is a unit-sphere cloud (density 24 ⇒ 25·24 = 600 sample points), so the
/// reconstructed surface's bounding box ≈ the unit sphere (each extent ≈ 2.0).
/// Drives the real `reverse_product()` and asserts the reported sample count
/// (600) and the reconstructed-mesh AABB extents (≈ 2.0, tol 0.05).
fn check_reverse(_app: &mut ValenxApp) -> CheckOutcome {
    let product = crate::reverse_workbench::reverse_product();
    let joined = product.lines.join(" • ");
    // "600 sampled points · k = 8 neighbours".
    if !joined.contains("600 sampled points") {
        return CheckOutcome::fail(format!("expected 600 sample points, got: {joined}"));
    }
    // The reconstructed mesh's bbox should approximate the unit sphere (±1).
    let Some(loaded) = product.mesh.as_ref() else {
        return CheckOutcome::fail("reverse product has no mesh".to_string());
    };
    let nodes = &loaded.mesh.nodes;
    if nodes.is_empty() {
        return CheckOutcome::fail("reconstructed mesh is empty".to_string());
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for n in nodes {
        for i in 0..3 {
            lo[i] = lo[i].min(n[i]);
            hi[i] = hi[i].max(n[i]);
        }
    }
    for i in 0..3 {
        let ext = hi[i] - lo[i];
        if !ext.is_finite() || (ext - 2.0).abs() > 0.05 {
            return CheckOutcome::fail(format!(
                "unit-sphere reconstruction axis {i} extent {ext:.3} vs 2.0"
            ));
        }
    }
    CheckOutcome::pass("unit-sphere recon: 600 pts, bbox ≈ 2.0³".to_string())
}

// ---------------------------------------------------------------------------
// Astrophysics + Sensors deep checks.
// ---------------------------------------------------------------------------

/// **blackhole** — Schwarzschild GR is exact closed form (geometrized units, the
/// panel's `M = 1` default with spin a = 0): event horizon `r₊ = 2M`, photon
/// sphere `r_ph = 3M`, ISCO `r_isco = 6M`, shadow radius `r_sh = 3√3 M ≈ 5.196 M`.
/// Drives the real `valenx_relativity` observables (`compute_observables`) and
/// asserts all four against their closed forms (tol 1e-3 M).
fn check_blackhole(app: &mut ValenxApp) -> CheckOutcome {
    let m = app.blackhole.mass; // default 1.0
    let readout = match crate::blackhole_workbench::compute_observables(&app.blackhole) {
        Ok(r) => r,
        Err(e) => return CheckOutcome::fail(format!("observables failed: {e}")),
    };
    let rp = parse_after_colon(&readout, "event horizon r+");
    let ph = parse_after_colon(&readout, "photon sphere");
    let isco = parse_after_colon(&readout, "ISCO (prograde)");
    let sh = parse_after_colon(&readout, "shadow radius");
    let (Some(rp), Some(ph), Some(isco), Some(sh)) = (rp, ph, isco, sh) else {
        return CheckOutcome::fail(format!("could not parse GR observables from: {readout}"));
    };
    let shadow_expected = 3.0 * 3.0_f64.sqrt() * m; // 3√3 M.
    if (rp - 2.0 * m).abs() >= 1e-3
        || (ph - 3.0 * m).abs() >= 1e-3
        || (isco - 6.0 * m).abs() >= 1e-3
        || (sh - shadow_expected).abs() >= 1e-3
    {
        return CheckOutcome::fail(format!(
            "Schwarzschild M={m}: r+={rp:.3} (2M), r_ph={ph:.3} (3M), ISCO={isco:.3} (6M), \
             shadow={sh:.3} (3√3 M={shadow_expected:.3})"
        ));
    }
    CheckOutcome::pass(format!(
        "Schwarzschild r+={rp:.2}M r_ph={ph:.2}M ISCO={isco:.2}M shadow={sh:.3}M (=3√3)"
    ))
}

/// **sensors** — LiDAR fan geometry is deterministic. The default scan lays
/// `N = 32` azimuth beams uniformly across the `90°` H-FOV in `[−45°, +45°]`
/// inclusive, so the angular resolution between adjacent beams is exactly
/// `FOV/(N−1) = 90/31 ≈ 2.9032°`. Drives the real `scan_lidar()`, reads the beam
/// azimuths, and asserts the measured step matches (tol 1e-4°) and the full span
/// is the 90° FOV.
fn check_sensors(app: &mut ValenxApp) -> CheckOutcome {
    let n = app.sensors.lidar.azimuth_steps; // 32
    let fov = app.sensors.lidar.h_fov_deg; // 90
    if n < 2 {
        return CheckOutcome::fail(format!("need ≥2 azimuth beams, got {n}"));
    }
    let result = match app.sensors.scan_lidar() {
        Ok(r) => r,
        Err(e) => return CheckOutcome::fail(format!("lidar scan failed: {e}")),
    };
    // Distinct azimuths (rad → deg), sorted; the fan is elevation-major so the
    // azimuth pattern repeats — collect the unique set.
    let mut az: Vec<f64> = result.scan.iter().map(|b| b.azimuth.to_degrees()).collect();
    az.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    az.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    if az.len() != n {
        return CheckOutcome::fail(format!("expected {n} distinct azimuths, got {}", az.len()));
    }
    let step = az[1] - az[0];
    let step_expected = fov / (n as f64 - 1.0);
    let span = az[az.len() - 1] - az[0];
    if (step - step_expected).abs() >= 1e-4 || (span - fov).abs() >= 1e-4 {
        return CheckOutcome::fail(format!(
            "azimuth step {step:.4}° vs FOV/(N−1) {step_expected:.4}°, span {span:.3}° vs {fov}°"
        ));
    }
    CheckOutcome::pass(format!(
        "LiDAR Δaz={step:.4}° = {fov}°/{}(N−1) over [{:.1}°,{:.1}°]",
        n - 1,
        az[0],
        az[az.len() - 1]
    ))
}

/// **autonomy** — scenario V&V closest-approach is exact geometry. The default
/// scenario coasts a vehicle from the origin east at 3 m/s for 25 ticks of
/// 0.1 s (⇒ 7.5 m, ending at x = 7.5) past a sphere obstacle at (12, 0) r = 1.
/// The closest path point is the endpoint, so the min clearance to the obstacle
/// surface is `|12 − 7.5| − 1 = 3.5 m`. Drives the real `valenx_autonomy_vnv`
/// run and asserts the reported `min_clearance_achieved` matches the independent
/// point-to-path distance (tol 1e-2 m) and that the scenario PASSES.
fn check_autonomy(app: &mut ValenxApp) -> CheckOutcome {
    let p = &app.autonomy.params;
    // Straight coast: end x = start_x + speed·dt·steps (heading 0° = +x).
    let end_x = p.start_x + p.speed * p.dt * p.steps as f64;
    // Closest approach of the [start_x, end_x] segment (on y=0) to the obstacle
    // centre (obstacle on y=0 too). Distance = |obstacle_x − clamp(obstacle_x to segment)|.
    let nearest_x = p
        .obstacle_x
        .clamp(p.start_x.min(end_x), p.start_x.max(end_x));
    let dist_centre =
        ((p.obstacle_x - nearest_x).powi(2) + (p.obstacle_y - 0.0_f64).powi(2)).sqrt();
    let clearance_expected = (dist_centre - p.obstacle_radius).max(0.0);

    let result = match app.autonomy.run() {
        Ok(r) => r,
        Err(e) => return CheckOutcome::fail(format!("autonomy run failed: {e}")),
    };
    if !result.overall_pass {
        return CheckOutcome::fail(
            "default scenario should PASS (clears the obstacle)".to_string(),
        );
    }
    let clr = result.min_clearance_achieved;
    if !clr.is_finite() || (clr - clearance_expected).abs() > 1e-2 {
        return CheckOutcome::fail(format!(
            "min clearance {clr:.3} m vs point-to-path {clearance_expected:.3} m"
        ));
    }
    CheckOutcome::pass(format!(
        "min clearance {clr:.2} m = |obst−path|−r ({clearance_expected:.2}); V&V PASS"
    ))
}

/// **thermo** — CO₂ with Peng–Robinson at 350 K / 5 MPa (the default seed).
/// Asserts the vapor compressibility `Z` is in the documented (0.7, 1.0) range
/// and that the saturation pressure at a subcritical T matches the NIST value
/// (CO₂ Psat at 273.15 K ≈ 3.49 MPa, within 5 %).
fn check_thermo(app: &mut ValenxApp) -> CheckOutcome {
    // Run the SAME compute the in-panel Compute button / `thermo.compute`
    // bridge fire, against the default CO₂/PR/350 K/5 MPa seed.
    crate::thermo_workbench::run(app);
    let Some(readout) = app.thermo.agent_readout() else {
        return CheckOutcome::fail("no readout after compute");
    };
    let Some(z) = parse_keyed(&readout, "Z_vap=") else {
        return CheckOutcome::fail(format!("could not parse Z_vap from: {readout}"));
    };
    if !(z.is_finite() && z > 0.7 && z < 1.0) {
        return CheckOutcome::fail(format!("Z_vap={z} outside NIST range (0.7,1.0)"));
    }

    // Second app: Psat at the subcritical 273.15 K via the agent control.
    let mut app2 = ValenxApp::default();
    if let Err(e) = app2.thermo.agent_set(
        "Temperature [K]",
        &crate::agent_commands::AgentValue::Float(273.15),
    ) {
        return CheckOutcome::fail(format!("set T failed: {e}"));
    }
    crate::thermo_workbench::run(&mut app2);
    let Some(r2) = app2.thermo.agent_readout() else {
        return CheckOutcome::fail("no Psat readout");
    };
    let Some(psat) = parse_keyed(&r2, "Psat(T)=") else {
        return CheckOutcome::fail(format!("could not parse Psat from: {r2}"));
    };
    if !psat.is_finite() || (psat - 3.49e6).abs() / 3.49e6 >= 0.05 {
        return CheckOutcome::fail(format!("Psat={psat:.4e} off NIST 3.49 MPa (>5%)"));
    }
    CheckOutcome::pass(format!(
        "Z_vap={z:.4} (NIST 0.7–1.0); Psat={psat:.3e}Pa≈3.49MPa"
    ))
}

/// **quantum** — the 2-qubit Bell pair `H(0)·CNOT(0→1)` (the default seed)
/// yields `p(|00>) = p(|11>) = 0.5`. Asserts both basis-state probabilities are
/// ≈ 0.5 (tol 1e-3) and the other two basis states carry ≈ 0.
fn check_quantum(app: &mut ValenxApp) -> CheckOutcome {
    crate::quantum_workbench::run(app);
    let Some(readout) = app.quantum.agent_readout() else {
        return CheckOutcome::fail("no readout after run");
    };
    let p00 = parse_keyed(&readout, "|00>=");
    let p11 = parse_keyed(&readout, "|11>=");
    match (p00, p11) {
        (Some(a), Some(b))
            if a.is_finite()
                && b.is_finite()
                && (a - 0.5).abs() < 1e-3
                && (b - 0.5).abs() < 1e-3 =>
        {
            CheckOutcome::pass(format!("Bell p(|00>)={a:.4} p(|11>)={b:.4} (≈0.5/0.5)"))
        }
        _ => CheckOutcome::fail(format!("Bell probs not ~0.5/0.5 in: {readout}")),
    }
}

/// **optics** — a thin lens with the object at 2 f (default `do=0.20`,
/// `f=0.10`) forms a real, inverted, unit-magnification image: the transverse
/// magnification `m ≈ -1` (the sign is the load-bearing physics — a real image
/// is inverted). Asserts `m` is finite, negative, and ≈ -1 (tol 1e-3).
fn check_optics(app: &mut ValenxApp) -> CheckOutcome {
    // Optics' readout is computed live from (do, f); run() also fills the panel
    // string. Either way the agent_readout carries `m=`.
    crate::optics_workbench::run(app);
    let Some(readout) = app.optics.agent_readout() else {
        return CheckOutcome::fail("no optics readout");
    };
    let Some(m) = parse_keyed(&readout, "m=") else {
        return CheckOutcome::fail(format!("could not parse m from: {readout}"));
    };
    if !m.is_finite() || m >= 0.0 {
        return CheckOutcome::fail(format!(
            "magnification m={m} not negative (image not inverted)"
        ));
    }
    if (m - (-1.0)).abs() >= 1e-3 {
        return CheckOutcome::fail(format!("m={m:.5} not ≈ -1 for object at 2f"));
    }
    CheckOutcome::pass(format!("thin-lens m={m:.5} (real inverted 1:1, sign<0)"))
}

/// **acoustics** — a free-field monopole (pulsating sphere) obeys the `1/r`
/// pressure law: doubling the observer distance halves the radiated pressure.
/// Runs the default seed at `r = 1 m`, then again at `r = 2 m` (via the agent
/// control), and asserts `p(2 m) / p(1 m) ≈ 0.5` (tol 2 %).
fn check_acoustics(app: &mut ValenxApp) -> CheckOutcome {
    crate::acoustics_workbench::run(app);
    let Some(r1) = app.acoustics.agent_readout() else {
        return CheckOutcome::fail("no readout at r=1m");
    };
    let Some(p1) = parse_keyed(&r1, "p=") else {
        return CheckOutcome::fail(format!("could not parse p at 1m: {r1}"));
    };

    let mut app2 = ValenxApp::default();
    if let Err(e) = app2.acoustics.agent_set(
        "observer distance (m)",
        &crate::agent_commands::AgentValue::Float(2.0),
    ) {
        return CheckOutcome::fail(format!("set distance failed: {e}"));
    }
    crate::acoustics_workbench::run(&mut app2);
    let Some(r2) = app2.acoustics.agent_readout() else {
        return CheckOutcome::fail("no readout at r=2m");
    };
    let Some(p2) = parse_keyed(&r2, "p=") else {
        return CheckOutcome::fail(format!("could not parse p at 2m: {r2}"));
    };
    if !(p1.is_finite() && p2.is_finite() && p1 > 0.0) {
        return CheckOutcome::fail(format!("non-finite pressures p1={p1} p2={p2}"));
    }
    let ratio = p2 / p1;
    if (ratio - 0.5).abs() >= 0.02 {
        return CheckOutcome::fail(format!("p(2m)/p(1m)={ratio:.4} not ≈0.5 (1/r law)"));
    }
    CheckOutcome::pass(format!("monopole p(2m)/p(1m)={ratio:.4} (1/r law)"))
}

/// **waveform** — the seeded clock+counter VCD parses to exactly the two
/// documented signals (`clk` + `cnt`). Asserts the readout reports `2 signal`
/// and names the clock — the parser's ground-truth trace.
fn check_waveform(app: &mut ValenxApp) -> CheckOutcome {
    crate::waveform_workbench::run(app);
    let Some(readout) = app.waveform.agent_readout() else {
        return CheckOutcome::fail("no readout after parse");
    };
    if readout.contains("parse failed") {
        return CheckOutcome::fail(format!("VCD parse failed: {readout}"));
    }
    if readout.contains("2 signal") && readout.contains("clk") {
        CheckOutcome::pass("VCD parsed 2 signals (clk + cnt)".to_string())
    } else {
        CheckOutcome::fail(format!("expected 2 signals incl clk, got: {readout}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The registry must cover **exactly** the authoritative `TabKind::TEMPLATES`
    /// (the 56 products) — no missing rows, no duplicates, and every id must
    /// round-trip through `TabKind::from_id` to its declared kind. This is the
    /// guard that keeps the self-test in lock-step with the product list: add a
    /// product to `TEMPLATES` and this fails until you add its check row.
    #[test]
    fn registry_covers_every_template() {
        let checks = product_checks();
        // 1:1 with the templates count (excludes the non-product Blank).
        assert_eq!(
            checks.len(),
            TabKind::TEMPLATES.len(),
            "self-test registry must have one row per product template"
        );

        // Every registered kind is a real template, and the id resolves to it.
        for pc in &checks {
            assert!(
                TabKind::TEMPLATES.contains(&pc.kind),
                "{} maps to a non-template kind {:?}",
                pc.id,
                pc.kind
            );
            assert_eq!(
                TabKind::from_id(pc.id),
                Some(pc.kind),
                "id {:?} must from_id to {:?}",
                pc.id,
                pc.kind
            );
        }

        // No duplicate ids and no duplicate kinds.
        let mut ids: Vec<&str> = checks.iter().map(|c| c.id).collect();
        ids.sort_unstable();
        let n = ids.len();
        ids.dedup();
        assert_eq!(ids.len(), n, "duplicate id in the self-test registry");

        let mut kinds: Vec<TabKind> = checks.iter().map(|c| c.kind).collect();
        kinds.sort_unstable_by_key(|k| *k as usize);
        let kn = kinds.len();
        kinds.dedup();
        assert_eq!(kinds.len(), kn, "duplicate kind in the self-test registry");

        // Every template has a row (set equality, given the count + uniqueness).
        for k in TabKind::TEMPLATES {
            assert!(
                checks.iter().any(|c| c.kind == k),
                "template {k:?} has no self-test row"
            );
        }
    }

    /// `parse_keyed` pulls the float that follows a key and stops at the unit.
    #[test]
    fn parse_keyed_extracts_the_value() {
        let r = "Thermo · CO2 · pengrobinson · T=350.00 K P=5.0000e6 Pa · \
                 Z_vap=0.93210 Z_liq=0.01234 · single real root · Psat(T)=6.4000e6 Pa";
        assert!((parse_keyed(r, "Z_vap=").unwrap() - 0.93210).abs() < 1e-6);
        assert!((parse_keyed(r, "T=").unwrap() - 350.0).abs() < 1e-9);
        assert!((parse_keyed(r, "Psat(T)=").unwrap() - 6.4e6).abs() < 1.0);
        assert!(parse_keyed(r, "nope=").is_none());
        // Negative + exponent forms parse.
        assert!((parse_keyed("m=-1.00000 · real", "m=").unwrap() - (-1.0)).abs() < 1e-9);
    }

    /// The bad-token scanner trips on non-finite numbers and the colon-form
    /// failure phrases, but NOT on a benign "… error" metric label or on
    /// substrings like "information" — so a working panel that displays a
    /// "Relative error" metric still passes.
    #[test]
    fn non_finite_scanner_is_word_boundaried() {
        // Clean.
        assert!(non_finite_or_error_token(&["Z = 0.93".into(), "ok".into()]).is_none());
        assert!(non_finite_or_error_token(&["more information here".into()]).is_none());
        // Benign metric labels must NOT trip (no colon ⇒ not a failure).
        assert!(non_finite_or_error_token(&["Relative error 0.4%".into()]).is_none());
        assert!(non_finite_or_error_token(&["reproj error: 0.3 px".into()]).is_some()); // this DOES (colon form)
                                                                                        // Non-finite numerics.
        assert!(non_finite_or_error_token(&["result = NaN".into()]).is_some());
        assert!(non_finite_or_error_token(&["p = inf Pa".into()]).is_some());
        // Failure phrases.
        assert!(non_finite_or_error_token(&["compute error: bad input".into()]).is_some());
        assert!(non_finite_or_error_token(&["Quantum run failed: nan".into()]).is_some());
    }

    /// Deep exemplar: thermo must PASS (CO₂ PR Z in NIST range + Psat match).
    #[test]
    fn deep_thermo_passes_known_value() {
        let mut app = ValenxApp::default();
        let out = check_thermo(&mut app);
        assert_eq!(
            out.status,
            Status::Pass,
            "thermo deep check: {}",
            out.detail
        );
        assert!(
            out.detail.contains("Z_vap"),
            "detail names Z: {}",
            out.detail
        );
    }

    /// Deep exemplar: quantum Bell pair must PASS at ≈0.5/0.5.
    #[test]
    fn deep_quantum_passes_bell_state() {
        let mut app = ValenxApp::default();
        let out = check_quantum(&mut app);
        assert_eq!(
            out.status,
            Status::Pass,
            "quantum deep check: {}",
            out.detail
        );
    }

    /// Deep exemplar: optics thin-lens magnification sign (m ≈ -1, negative).
    #[test]
    fn deep_optics_passes_magnification_sign() {
        let mut app = ValenxApp::default();
        let out = check_optics(&mut app);
        assert_eq!(
            out.status,
            Status::Pass,
            "optics deep check: {}",
            out.detail
        );
    }

    /// Deep exemplar: acoustics monopole 1/r law (p halves at 2×distance).
    #[test]
    fn deep_acoustics_passes_inverse_r_law() {
        let mut app = ValenxApp::default();
        let out = check_acoustics(&mut app);
        assert_eq!(
            out.status,
            Status::Pass,
            "acoustics deep check: {}",
            out.detail
        );
    }

    /// Deep exemplar: waveform VCD parse (clk + cnt).
    #[test]
    fn deep_waveform_passes_vcd_parse() {
        let mut app = ValenxApp::default();
        let out = check_waveform(&mut app);
        assert_eq!(
            out.status,
            Status::Pass,
            "waveform deep check: {}",
            out.detail
        );
    }

    /// `parse_after_colon` pulls the value from a padded `label … : value` line.
    #[test]
    fn parse_after_colon_handles_padded_columns() {
        let r = "ISENTROPIC\n  A/A*             : 1.687500\nNORMAL SHOCK (M1 = 2.0000)\n  M2               : 0.577350\n  p2/p1            : 4.500000\n";
        assert!((parse_after_colon(r, "A/A*").unwrap() - 1.6875).abs() < 1e-6);
        assert!((parse_after_colon(r, "M2").unwrap() - 0.57735).abs() < 1e-6);
        assert!((parse_after_colon(r, "p2/p1").unwrap() - 4.5).abs() < 1e-6);
        assert!(parse_after_colon(r, "nope").is_none());
    }

    /// Aerospace DEEP: gas dynamics M=2,γ=1.4 → A/A*=1.6875, shock M2=0.5774, p2/p1=4.5.
    #[test]
    fn deep_gasdynamics_passes_naca1135() {
        let mut app = ValenxApp::default();
        let out = check_gasdynamics(&mut app);
        assert_eq!(out.status, Status::Pass, "gasdynamics: {}", out.detail);
        assert!(out.detail.contains("A/A*"), "detail: {}", out.detail);
    }

    /// Aerospace DEEP: astro Hohmann LEO→GEO total Δv vs closed-form Kepler/vis-viva.
    #[test]
    fn deep_astro_passes_hohmann_delta_v() {
        let mut app = ValenxApp::default();
        let out = check_astro(&mut app);
        assert_eq!(out.status, Status::Pass, "astro: {}", out.detail);
        assert!(out.detail.contains("Δv"), "detail: {}", out.detail);
    }

    /// Aerospace DEEP: rocket LV-1 ideal Δv vs Tsiolkovsky 2-stage sum.
    #[test]
    fn deep_rocket_passes_tsiolkovsky() {
        let mut app = ValenxApp::default();
        let out = check_rocket(&mut app);
        assert_eq!(out.status, Status::Pass, "rocket: {}", out.detail);
    }

    /// Aerospace DEEP: engine c* vs 1-D Vandenkerckhove nozzle theory.
    #[test]
    fn deep_engine_passes_cstar() {
        let mut app = ValenxApp::default();
        let out = check_engine(&mut app);
        assert_eq!(out.status, Status::Pass, "engine: {}", out.detail);
    }

    /// Aerospace DEEP: rotor BEMT hover figure of merit vs momentum theory.
    #[test]
    fn deep_rotor_passes_momentum_theory() {
        let mut app = ValenxApp::default();
        let out = check_rotor(&mut app);
        assert_eq!(out.status, Status::Pass, "rotor: {}", out.detail);
    }

    /// Aerospace DEEP: uas multirotor hover power vs disk-loading momentum theory.
    #[test]
    fn deep_uas_passes_hover_power() {
        let mut app = ValenxApp::default();
        let out = check_uas(&mut app);
        assert_eq!(out.status, Status::Pass, "uas: {}", out.detail);
    }

    /// Machine-design DEEP: springs helical rate k = G·d⁴/(8·D³·n).
    #[test]
    fn deep_springs_passes_rate() {
        let mut app = ValenxApp::default();
        let out = check_springs(&mut app);
        assert_eq!(out.status, Status::Pass, "springs: {}", out.detail);
        assert!(out.detail.contains("N/mm"), "detail: {}", out.detail);
    }

    /// Machine-design DEEP: gears spur pitch dia m·z + ratio z2/z1.
    #[test]
    fn deep_gears_passes_geometry() {
        let mut app = ValenxApp::default();
        let out = check_gears(&mut app);
        assert_eq!(out.status, Status::Pass, "gears: {}", out.detail);
    }

    /// Machine-design DEEP: fasteners M6 pitch dia + ISO-898 tensile stress area.
    #[test]
    fn deep_fasteners_passes_stress_area() {
        let mut app = ValenxApp::default();
        let out = check_fasteners(&mut app);
        assert_eq!(out.status, Status::Pass, "fasteners: {}", out.detail);
    }

    /// Machine-design DEEP: frames IPE 200 cross-section area closed form.
    #[test]
    fn deep_frames_passes_area() {
        let mut app = ValenxApp::default();
        let out = check_frames(&mut app);
        assert_eq!(out.status, Status::Pass, "frames: {}", out.detail);
    }

    /// Machine-design DEEP: collision AABB overlap both ways (disjoint / overlap).
    #[test]
    fn deep_collision_passes_aabb_overlap() {
        let mut app = ValenxApp::default();
        let out = check_collision(&mut app);
        assert_eq!(out.status, Status::Pass, "collision: {}", out.detail);
    }

    /// Simulation DEEP: fem 3-D cantilever tip δ vs Euler–Bernoulli PL³/3EI.
    #[test]
    fn deep_fem_passes_beam_deflection() {
        let mut app = ValenxApp::default();
        let out = check_fem(&mut app);
        assert_eq!(out.status, Status::Pass, "fem: {}", out.detail);
        assert!(out.detail.contains("PL³/3EI"), "detail: {}", out.detail);
    }

    /// Simulation DEEP: fields descriptive stats of 1..5 (mean=3, σ=√2).
    #[test]
    fn deep_fields_passes_descriptive_stats() {
        let mut app = ValenxApp::default();
        let out = check_fields(&mut app);
        assert_eq!(out.status, Status::Pass, "fields: {}", out.detail);
    }

    /// Simulation DEEP: nodegraph default graph evaluates to 4.0.
    #[test]
    fn deep_nodegraph_passes_arithmetic() {
        let mut app = ValenxApp::default();
        let out = check_nodegraph(&mut app);
        assert_eq!(out.status, Status::Pass, "nodegraph: {}", out.detail);
    }

    /// Simulation DEEP: reactdyn H2 AIMD NVE energy conservation.
    #[test]
    fn deep_reactdyn_passes_energy_conservation() {
        let mut app = ValenxApp::default();
        let out = check_reactdyn(&mut app);
        assert_eq!(out.status, Status::Pass, "reactdyn: {}", out.detail);
    }

    /// Simulation DEEP: mbd pendulum period vs 2π√(L/g) + energy conservation.
    #[test]
    fn deep_mbd_passes_pendulum_period() {
        let mut app = ValenxApp::default();
        let out = check_mbd(&mut app);
        assert_eq!(out.status, Status::Pass, "mbd: {}", out.detail);
    }

    /// Life-sciences DEEP: genetics sequence GC content + length (ATGC → 50%, 4).
    #[test]
    fn deep_genetics_passes_gc_content() {
        let mut app = ValenxApp::default();
        let out = check_genetics(&mut app);
        assert_eq!(out.status, Status::Pass, "genetics: {}", out.detail);
        assert!(out.detail.contains("GC="), "detail: {}", out.detail);
    }

    /// Life-sciences DEEP: neuro unmyelinated conduction velocity v = k·√d.
    #[test]
    fn deep_neuro_passes_conduction_velocity() {
        let mut app = ValenxApp::default();
        let out = check_neuro(&mut app);
        assert_eq!(out.status, Status::Pass, "neuro: {}", out.detail);
    }

    /// Life-sciences DEEP: variant HGVS parse p.R273H → ProteinSub R(273)H.
    #[test]
    fn deep_variant_passes_hgvs_parse() {
        let mut app = ValenxApp::default();
        let out = check_variant(&mut app);
        assert_eq!(out.status, Status::Pass, "variant: {}", out.detail);
    }

    /// Life-sciences DEEP: ppi shortest-path GUARD→EFF-A = 1 hop.
    #[test]
    fn deep_ppi_passes_shortest_path() {
        let mut app = ValenxApp::default();
        let out = check_ppi(&mut app);
        assert_eq!(out.status, Status::Pass, "ppi: {}", out.detail);
    }

    /// Life-sciences DEEP: morphogenesis Gray–Scott U,V∈[0,1] invariant + evolves.
    #[test]
    fn deep_morphogenesis_passes_bounded_invariant() {
        let mut app = ValenxApp::default();
        let out = check_morphogenesis(&mut app);
        assert_eq!(out.status, Status::Pass, "morphogenesis: {}", out.detail);
    }

    /// Civil & AEC DEEP: geomatics Cambridge→Paris great-circle ≈ 404.3 km.
    #[test]
    fn deep_geomatics_passes_haversine() {
        let mut app = ValenxApp::default();
        let out = check_geomatics(&mut app);
        assert_eq!(out.status, Status::Pass, "geomatics: {}", out.detail);
    }

    /// Civil & AEC DEEP: piping ASME B36.10 NPS 2 Sch 40 OD = 60.325 mm.
    #[test]
    fn deep_piping_passes_asme_dimensions() {
        let mut app = ValenxApp::default();
        let out = check_piping(&mut app);
        assert_eq!(out.status, Status::Pass, "piping: {}", out.detail);
        assert!(out.detail.contains("60.325"), "detail: {}", out.detail);
    }

    /// Civil & AEC DEEP: hvac Darcy–Weisbach duct ΔP = 15.05 Pa.
    #[test]
    fn deep_hvac_passes_pressure_drop() {
        let mut app = ValenxApp::default();
        let out = check_hvac(&mut app);
        assert_eq!(out.status, Status::Pass, "hvac: {}", out.detail);
    }

    /// Civil & AEC DEEP: interior floor-plan bookkeeping (6×4 m, 2 pieces).
    #[test]
    fn deep_interior_passes_floor_plan() {
        let mut app = ValenxApp::default();
        let out = check_interior(&mut app);
        assert_eq!(out.status, Status::Pass, "interior: {}", out.detail);
    }

    /// CAD & mesh DEEP: cad parametric CSG volume = 1−π/16.
    #[test]
    fn deep_cad_passes_csg_volume() {
        let mut app = ValenxApp::default();
        let out = check_cad(&mut app);
        assert_eq!(out.status, Status::Pass, "cad: {}", out.detail);
    }

    /// CAD & mesh DEEP: brep boolean difference volume ≈ 1−π/16 (tessellated).
    #[test]
    fn deep_brep_passes_boolean_volume() {
        let mut app = ValenxApp::default();
        let out = check_brep(&mut app);
        assert_eq!(out.status, Status::Pass, "brep: {}", out.detail);
    }

    /// CAD & mesh DEEP: mesh canonical LV-1 AABB extents.
    #[test]
    fn deep_mesh_passes_lv1_aabb() {
        let mut app = ValenxApp::default();
        let out = check_mesh(&mut app);
        assert_eq!(out.status, Status::Pass, "mesh: {}", out.detail);
    }

    /// CAD & mesh DEEP: sheetmetal bend allowance (π·θ/180)·(R+K·t).
    #[test]
    fn deep_sheetmetal_passes_bend_allowance() {
        let mut app = ValenxApp::default();
        let out = check_sheetmetal(&mut app);
        assert_eq!(out.status, Status::Pass, "sheetmetal: {}", out.detail);
    }

    /// CAD & mesh DEEP: draft2d exact 2-D geometry (4 entities, 60×40 extent).
    #[test]
    fn deep_draft2d_passes_geometry() {
        let mut app = ValenxApp::default();
        let out = check_draft2d(&mut app);
        assert_eq!(out.status, Status::Pass, "draft2d: {}", out.detail);
    }

    /// CAD & mesh DEEP: animate keyframe interpolation (Linear t=1s = π/2).
    #[test]
    fn deep_animate_passes_interpolation() {
        let mut app = ValenxApp::default();
        let out = check_animate(&mut app);
        assert_eq!(out.status, Status::Pass, "animate: {}", out.detail);
    }

    /// CAD & mesh DEEP: reverse unit-sphere reconstruction (600 pts, bbox ≈ 2³).
    #[test]
    fn deep_reverse_passes_sphere_recon() {
        let mut app = ValenxApp::default();
        let out = check_reverse(&mut app);
        assert_eq!(out.status, Status::Pass, "reverse: {}", out.detail);
    }

    /// Astrophysics DEEP: blackhole Schwarzschild r₊/photon/ISCO/shadow.
    #[test]
    fn deep_blackhole_passes_schwarzschild() {
        let mut app = ValenxApp::default();
        let out = check_blackhole(&mut app);
        assert_eq!(out.status, Status::Pass, "blackhole: {}", out.detail);
        assert!(out.detail.contains("ISCO"), "detail: {}", out.detail);
    }

    /// Sensors DEEP: LiDAR angular resolution = FOV/(N−1) = 90°/31.
    #[test]
    fn deep_sensors_passes_angular_resolution() {
        let mut app = ValenxApp::default();
        let out = check_sensors(&mut app);
        assert_eq!(out.status, Status::Pass, "sensors: {}", out.detail);
    }

    /// Sensors DEEP: autonomy V&V min clearance = point-to-path distance + PASS.
    #[test]
    fn deep_autonomy_passes_min_clearance() {
        let mut app = ValenxApp::default();
        let out = check_autonomy(&mut app);
        assert_eq!(out.status, Status::Pass, "autonomy: {}", out.detail);
    }

    /// The `--id` filter selects exactly one product and runs it end-to-end.
    #[test]
    fn id_filter_runs_a_single_product() {
        let rep = run_self_tests(&Filter::Id("thermo".into()));
        assert_eq!(rep.lines.len(), 1);
        assert_eq!(rep.lines[0].id, "thermo");
        assert_eq!(rep.lines[0].status, Status::Pass, "{}", rep.lines[0].detail);
        assert!(rep.ok());
    }

    /// The `--group` filter selects every product in a group (and only those).
    #[test]
    fn group_filter_selects_the_group() {
        let rep = run_self_tests(&Filter::Group("Aerospace".into()));
        assert!(!rep.lines.is_empty());
        for l in &rep.lines {
            let kind = TabKind::from_id(&l.id).unwrap();
            assert_eq!(kind.group(), "Aerospace", "{} not Aerospace", l.id);
        }
    }

    /// Skips are reported as SKIP (not pass, not fail) with their reason.
    #[test]
    fn skips_are_honest() {
        let rep = run_self_tests(&Filter::Id("render".into()));
        assert_eq!(rep.lines.len(), 1);
        assert_eq!(rep.lines[0].status, Status::Skip);
        assert_eq!(rep.lines[0].detail, "gpu-render");
        assert_eq!(rep.skipped, 1);
    }

    /// Filter parsing: `--id` wins over `--group`; bare ⇒ All.
    #[test]
    fn filter_parsing() {
        assert!(matches!(Filter::from_args(&[]), Filter::All));
        assert!(matches!(
            Filter::from_args(&["--group".into(), "Simulation".into()]),
            Filter::Group(g) if g == "Simulation"
        ));
        assert!(matches!(
            Filter::from_args(&["--id".into(), "thermo".into()]),
            Filter::Id(i) if i == "thermo"
        ));
        // --id wins.
        assert!(matches!(
            Filter::from_args(&["--group".into(), "X".into(), "--id".into(), "fem".into()]),
            Filter::Id(i) if i == "fem"
        ));
    }

    /// The compact report renders one line per product plus a tally, and the
    /// tally arithmetic is internally consistent.
    #[test]
    fn report_renders_compact_and_consistent() {
        let rep = run_self_tests(&Filter::Group("Aerospace".into()));
        let text = rep.render();
        // One body line per product + one tally line.
        assert_eq!(text.lines().count(), rep.lines.len() + 1);
        assert!(text.contains("PASS") || text.contains("SKIP") || text.contains("FAIL"));
        assert_eq!(rep.passed + rep.failed + rep.skipped, rep.lines.len());
        assert!(text.contains("product(s):"));
    }
}
