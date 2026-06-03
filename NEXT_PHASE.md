# Next phase — scope

**Goal:** ship `0.1.0` (the first tagged public alpha) by closing
the remaining Phase 10 acceptance items, then start Phase 11.

This document is the actionable bridge between
[STATUS.md](./STATUS.md) (current snapshot) and
[ROADMAP.md](./ROADMAP.md) (16-phase long arc). It scopes the
work into commit-shaped milestones with rough effort estimates,
identifies the critical path, and calls out the side quests that
deserve their own focused sessions.

> **Calibration.** Effort estimates are "engineer-days" — calendar
> time depends on parallelism, review cadence, and how many of the
> known unknowns surface. Treat the totals as planning anchors, not
> guarantees.

---

## TL;DR

| Lane | Effort | Status | Critical for 0.1.0? |
|---|---|---|---|
| **A. i18n pipeline** | ~3 days | 🟡 library landed (`valenx-i18n` crate, en-US baseline, pseudo-locale, `format_with`); UI wiring TBD | optional |
| **B. First-run wizard** | ~4 days | 🟡 logic landed (`valenx-first-run` crate, `EnvironmentReport` / `FirstRunDecision`, per-OS install hints); egui shell TBD | **yes** — onboarding |
| **C. In-app crash reporter** | ~3 days | ✅ shipped (`valenx-crash-reporter` crate, panic hook, sanitiser, subprocess integration test) | **yes** — unattended data |
| **D. Installer signing pipeline** | ~5 days | ✅ shipped (`.github/workflows/release.yml`, deb / rpm / .app / .msi configs, RELEASING.md). Certs deferred to v0.2.0 — pipeline degrades gracefully to unsigned artefacts, which is the convention for pre-alpha OSS releases | optional for 0.1.0 (was: yes — re-classified) |
| **E. Theme snapshot tests** | ~2 days | ⚪ planned | optional |
| **F. Accessibility audit** | ~3 days | 🟡 contrast subset landed (`valenx-a11y` crate, WCAG 2.1 ratios + AA / AAA gates, audit helpers); accesskit shell TBD | partial |
| **G. Visual polish pass** | ~2 days | ⚪ planned | optional |
| **H. Per-pick viewport readouts** | ~7 days | ⚪ planned | defer to v0.2 |
| **I. ssh-tail log streaming** | ~5 days | ⚪ planned | defer to Phase 11 |
| **J. OCCT FFI shim** | ~15-20 days | ⚪ planned | defer — multi-week |
| **K. Real-cluster CI** | ~10 days | ⚪ planned | defer to Phase 11 |

> 🟢 = done · 🟡 = library / CI landed, GUI shell still to wire ·
> ⚪ = not started

**Library / CI work for the four critical-path lanes (B, C, D, F)
is now landed.** The remaining work to reach 0.1.0 is:
- Wire `valenx-first-run` into the GUI as an egui dialog (lane B GUI half)
- Wire `valenx-crash-reporter::install_panic_hook` into `main.rs` + add a Settings → Privacy toggle (lane C GUI half)
- Procure the Apple Developer ID + Authenticode certs (lane D — kicked off; ~1 week procurement lead time)
- Wire the design-token colour pairs through `valenx-a11y::audit` in a CI test (lane F gate)
- A → E → G optional polish lanes

**Single-engineer critical path to 0.1.0:** `C → B → F (partial) → tag` —
about 1-2 weeks of focused work. The cert procurement that was
flagged as a 1-week parallel item turns out to be optional for an
alpha (`release.yml` ships unsigned artefacts cleanly when secrets
are unset). Sign-up happens in v0.2.0 when the user base justifies
the $300+/yr ongoing cost.

**Three-engineer parallel path:** ship 0.1.0 in ~2 calendar weeks,
plus headroom to start Phase 11 work in parallel.

---

## Phase 10 — close the acceptance checklist

The remaining items from
[`docs/src/phases/phase-10-ux-polish.md`](./docs/src/phases/phase-10-ux-polish.md):

### A. i18n pipeline (~3 days)

- Pick a runtime crate (`fluent-rs` is the obvious choice — used
  by Firefox, mature, supports plurals + ICU MessageFormat).
- Externalise every user-visible string into `locales/en-US/*.ftl`.
- Wire a `t!("key")` macro that resolves at compile time when the
  locale is fixed and at runtime otherwise.
- Add a CI lint that scans Rust UI code for string literals
  matching `\.label\("[A-Z]"`-style patterns and rejects them.
- Pseudo-locale (`Pseudo`) for dev builds — wraps every string in
  brackets to surface hard-coded leaks visually.

**Acceptance:** zero hard-coded UI strings, every label routed
through the catalogue, pseudo-locale visibly catches new ones.

**Risk:** if no good lint exists, may need a small proc-macro.
Mitigation: ship the catalogue + manual review for v0.1.0, lint
in v0.2.0.

### B. First-run wizard (~4 days)

- New `first_run::detect_environment()` helper that calls each
  registered adapter's `probe()` and returns a per-adapter status
  vector.
- New `FirstRunWizard` egui modal that renders the status, surfaces
  install hints (per OS), and lets the user check / uncheck which
  ones they expect to use.
- Persist the answer to `settings.json` so the wizard runs at most
  once per user.
- Skippable via a "Don't show this again" checkbox or `--skip-wizard`
  CLI arg.

**Acceptance:** fresh install on a clean VM with no solvers
shows the wizard; ticking "I have OpenFOAM" but failing the probe
surfaces an actionable error; the wizard never re-opens.

**Touches:** `valenx-app/src/first_run.rs` (new), `ValenxApp::new`
init flow, `Settings::first_run_complete: bool`.

### C. In-app crash reporter (~3 days)

- Set a `std::panic::set_hook` that captures the panic payload +
  stack + sanitised metadata (no file paths, no project IDs, just
  adapter / version / platform).
- Write the report to `~/.local/state/valenx/crashes/<ts>.json`
  on every crash.
- Opt-in upload to a configurable endpoint via Settings → Privacy.
- Display "Send crash report" prompt on next launch if any
  unsent reports queued.

**Acceptance:** trigger a panic in a test build, confirm the
report lands on disk; opt-in upload sends correctly to a stub
endpoint; opt-out path leaves the file local-only.

**Risk:** stack symbolication on Windows requires `.pdb`s shipped
with the installer. Mitigate by including them in the package step
(adds a few MB to the download).

### D. Installer signing pipeline (~5 days) — pipeline shipped, certs deferred

- **macOS:** `cargo-bundle` → `.app` → `codesign` → `notarytool`
  submit → `stapler staple`. Needs Apple Developer ID
  certificate + an Apple ID with notarisation entitlement.
- **Windows:** `cargo-wix` → `.msi` → `signtool` with an
  Authenticode certificate (DigiCert / Sectigo / GlobalSign).
- **Linux:** `cargo-deb` for `.deb`, `cargo-generate-rpm` for
  `.rpm`. Reproducible-build flags so checksums match across
  builders. AppImage as a fallback for distros we don't pin.
- CI matrix in `.github/workflows/release.yml`: build all three
  platforms via `gh workflow run release.yml`, upload signed
  artefacts to the GitHub Release. (Round-8 removed the tag-push
  trigger; release builds are manual `workflow_dispatch` only so an
  accidental retag never re-runs the signing pipeline.)

**Acceptance:** tag `0.1.0`, run `gh workflow run release.yml`,
watch CI produce signed artefacts on all three platforms, install
on a fresh machine and confirm no "unidentified developer" /
"untrusted publisher" warnings.

**Procurement is optional for v0.1.0-alpha.1.** The pipeline
already ships unsigned artefacts when the cert secrets are unset
— that's the conventional shape for OSS pre-alphas (rustup, uv,
ripgrep all started this way). Users see one Gatekeeper /
SmartScreen warning the first time they run, then it's gone.
Trade-off: $0/yr unsigned vs ~$99/yr Apple + ~$100-300/yr
Authenticode for the smoother out-of-the-box UX. Defer to v0.2.0
when the user base justifies the cost or a sponsor covers it.

### E. Theme snapshot tests (~2 days)

- Capture rendered images of every panel in light + dark themes
  using `egui_kittest` or a custom harness over the existing
  app frames.
- Pixel-diff CI gate; intentional changes update the baseline
  with a `cargo test -- --update-snapshots` flow.
- Three reference resolutions (1280×720 / 1920×1080 / 2560×1440)
  to catch DPI-dependent layout drift.

**Acceptance:** changing any pixel in any panel without updating
the baseline fails the test; updating the baseline produces a
diff that the reviewer can sanity-check.

**Touches:** `crates/valenx-app/tests/theme_snapshots/` (new),
plus a new integration crate or feature flag that doesn't pull
the JSON-file persistence path.

### F. Accessibility audit (~3 days)

- Wire `accesskit` (egui has integration) for screen-reader
  narration on every interactive control.
- Add a colour-contrast test: every fg/bg pair in the design-tokens
  module passes WCAG 2.1 AA (4.5:1 for normal text, 3:1 for large).
- Keyboard nav coverage: every button reachable via Tab, Space /
  Enter activates, Escape closes modals.

**Acceptance:** automated contrast test passes, keyboard nav test
visits every focusable control on the case-browser + viewport +
results pane.

### G. Visual polish pass (~2 days)

- Tooltip timing tokens (300 ms hover delay, 100 ms hide grace).
- Animation budget: nothing runs at < 60 fps, motion respects
  `prefers-reduced-motion` (via the system setting on macOS / Windows).
- Hairline / spacing audit: every panel uses tokens from
  `valenx-design-tokens`, no magic numbers.

**Acceptance:** review every visible surface in the app once with
the tokens cheat sheet; commit any drift that surfaces.

---

## Post-0.1.0 side quests

These items are documented as known gaps in STATUS.md. They don't
block the alpha tag but each cleanly slots into an upcoming phase:

### H. Per-pick numerical readouts (~7 days)

Closes the most-asked-for viewport gap. Click on a mesh node →
tooltip shows the field value at that node.

- Screen-space → world-space ray construction (camera + projection
  matrix already exists in `viewport::Camera`).
- Acceleration structure on the canonical Mesh: BVH or octree over
  element bounding boxes.
- Per-pick lookup: ray-vs-element intersection → barycentric
  interpolation of the active scalar.
- Tooltip rendering: world-space anchor in egui's painter.

Belongs to Phase 10 visual but lands as its own milestone because
the BVH work is non-trivial. Defer to v0.2.0.

### I. ssh-tail log streaming (~5 days) — Phase 11 foundation

The `valenx-executor-slurm` crate already has the after-the-fact
log reader (`read_slurm_log_tail`) and the ssh argv constructors
(`build_ssh_wrapped_command`). The remaining piece is a persistent
ssh session that streams `tail -f slurm-<id>.out` line-by-line into
the run's progress channel.

- New `slurm::ssh_tail::start(host, native_id)` that spawns
  `ssh <host> tail -f --pid=<sentinel> ...` and yields a
  `Receiver<String>`.
- Wire the receiver into `RunContext::log` so the GUI's log panel
  shows live remote output.
- Cleanup: kill the ssh process when the run ends or the user
  cancels.

**Risk:** persistent ssh sessions need careful kill-on-drop
discipline so a cancelled run doesn't leak processes on the
client.

### J. OCCT FFI shim (~15-20 days) — Phase 2 closure

Genuinely multi-week. The shape:

- New `crates/occt-sys/` crate with a `build.rs` that probes for
  OpenCASCADE via `pkg-config` (Linux) or environment vars
  (Windows / macOS).
- `cxx`-based bindings for the BRep kernel: `BRepPrimAPI_MakeBox`,
  `BRepBuilderAPI_MakeShape`, the STEP / IGES readers (`STEPControl_Reader`,
  `IGESControl_Reader`), and the mesh discretiser
  (`BRepMesh_IncrementalMesh`).
- Wire them into `valenx-adapter-occt::prepare()` / `run()` /
  `collect()`, replacing the current `not_implemented` stubs.
- Test fixtures: a STEP file, an IGES file, a BREP file. Each
  imports cleanly, exports as STL, round-trips through canonical
  Mesh.

**Hard prerequisites:** OpenCASCADE 7.6+ installed on the build
machine. Multi-platform CI needs OCCT in every runner image.

**Owner profile:** dedicated C++/Rust FFI contributor, ideally
familiar with OpenCASCADE's quirky lifecycle model.

### K. Real-cluster CI integration (~10 days) — Phase 11 acceptance

End-to-end runs of SU2 / OpenFOAM / GROMACS against checked-in
fixtures, on real cluster hardware, on every CI run.

- Self-hosted runner with SLURM + the three solvers installed.
- Cluster fixture: `tests/fixtures/cluster.valenx/` with a
  cavity-flow case sized for ~1 minute on 4 cores.
- New CI workflow `cluster-e2e.yml` that runs on a cron schedule
  (not every PR — too expensive).

**Risk:** self-hosted runners are an SRE / DevOps commitment. May
warrant a sponsorship arrangement with an HPC partner.

---

## Phase 11+ — what's after 0.1.0

Once the alpha ships, the natural next-phase split (per ROADMAP.md):

**Phase 11 — HPC / cluster execution** (~4-6 weeks)
- ssh-tail (I above) lands first — foundational.
- Remote viewport: headless render on cluster, pixel stream back.
- Credential management with OS keychain.
- Result-cache pinning rules.

**Phase 12 — Optimisation + adjoint workflows** (~6-8 weeks)
- DAKOTA integration for parameter sweeps.
- Adjoint solver wiring (SU2 has the strongest adjoint story).
- Pareto-front visualisation.

**Phase 13 — ML-assisted surrogates** (~8-12 weeks)
- The `valenx-export` crate already ships dataset-export helpers.
- Phase 13 wraps an ONNX-runtime inference layer + a "run
  surrogate then refine with high-fidelity" workflow.

**Phases 14-16** (plugin marketplace / enterprise / stewardship)
remain documentation-only. Plenty of time to get there.

---

## Recommended ordering

### Solo engineer (sequential)

```
Week 1:  D (signing) + crash-reporter foundation
Week 2:  C (crash reporter complete) + B (first-run wizard)
Week 3:  F (accessibility) + E (theme snapshots)
Week 4:  A (i18n) + G (polish) + tag 0.1.0
Week 5:  Begin I (ssh-tail) — Phase 11 foundation
Week 6+: H (per-pick) or J (OCCT FFI) depending on pull
```

### Two-engineer parallel

```
Engineer A: D → C → B → tag 0.1.0   (release-pipeline + onboarding)
Engineer B: A + E + F + G            (UX polish)
Tag converges Week 2-3
Then both: I + H or J in parallel
```

### Three-engineer parallel

```
Engineer A: D → C → B
Engineer B: A → E → F → G
Engineer C: I (ssh-tail starting now — independent of 0.1.0)
Tag 0.1.0 Week 2
Engineers A+B: H (per-pick)
Engineer C: continue Phase 11
```

---

## Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| Apple Developer ID enrolment delay | Low | Deferred to v0.2.0 — pipeline ships unsigned for 0.1.0-alpha |
| Authenticode certificate lead time | Low | Deferred to v0.2.0 — same rationale |
| OpenCASCADE FFI complexity (lane J) | High | Treat as its own multi-week project; don't block 0.1.0 |
| SLURM CI infrastructure cost (lane K) | High | Defer to Phase 11; consider HPC partnership |
| Pseudo-locale lint may need proc-macro | Low | Ship catalogue + manual review for 0.1.0; lint in 0.2.0 |
| Theme snapshot test flakiness on mixed DPI | Low | Pin reference resolutions; allow ±1 px tolerance |

---

## Open questions

1. **Telemetry policy.** Crash reporter (lane C) needs an upload
   endpoint. Decision: self-hosted vs Sentry vs custom? Affects
   privacy doc + Settings UI.

2. **i18n scope at 0.1.0.** Ship with en-US only and a "translations
   welcome" CONTRIBUTING note, or pre-seed one or two locales?
   Recommendation: en-US only, accept community PRs.

3. **Installer cadence.** GitHub Releases on every tag, or also a
   nightly channel? Recommendation: tags only for 0.1.0; nightly
   for 0.2.0+.

4. **OCCT FFI ownership.** Find a dedicated contributor with
   OpenCASCADE experience, or scope the scaffold-only crate (no
   actual binding) for 0.1.0 and treat the binding as a Phase 2.5
   project? Recommendation: scaffold-only, write a clear
   "wanted: OCCT contributor" section in CONTRIBUTING.md.

---

## Tracking

This doc lives next to STATUS.md and gets updated as lanes land.
Each lane should reference a tracking issue (or RFC if the design
warrants one). When a lane completes, mark it done here and move
the entry into the matching CHANGELOG entry.
