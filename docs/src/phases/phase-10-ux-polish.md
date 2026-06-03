# Phase 10 — UX polish + first public release

**Status:** 🟡 Partially landed during Phase 1 — the workflow loop
(prepare → edit → run → inspect) is live, including click-to-run /
prepare-only / open-workdir-in-host-browser / per-case run-history
badges / adapter-status badges / field-coloured wireframe overlay
with colour-bar legend / time-series scrubbing for transient runs /
clickable field picker in the Results pane. The bigger Phase 10
items below (i18n, signed installers, crash reporter, first-run
wizard, theme snapshot tests) are still planned.

## Goal

Move from "everything works" to "everything feels right" — and ship
the first stable public release.

## Capability inventory

- Visual polish pass: hairline consistency, tooltip timing,
  animation budget respected, motion tokens enforced.
- i18n: all user-visible strings externalised, pseudo-locale enabled
  in dev builds to catch hard-coded strings.
- Theme parity: light/dark render pixel-identical snapshot tests.
- First-run wizard: detects installed OSS tools, offers one-click
  install, writes per-user preferences.
- In-app crash reporter (opt-in) with sanitised stack + project
  header.
- Installer signing + notarisation for macOS, code-signing on
  Windows, reproducible deb/rpm builds for Linux.

## What landed early (during Phase 1)

These shipped during Phase 1's foundation work and the results-
rendering arc — not the full Phase 10 polish, but enough that the
workflow loop is genuinely usable today:

- [x] Click-to-run case workflow (browser → Run, F5 keystroke,
      command-palette entry).
- [x] Click-to-prepare-without-running (write the dict tree,
      inspect / edit by hand, then "Run from prepared workdir").
- [x] "Open prepared / run workdir in host file browser"
      cross-platform (Explorer / Finder / xdg-open).
- [x] Adapter status badges in the case browser
      (Ready / Missing / Outdated / Broken / Disabled / Unregistered).
- [x] Per-case run-history badges (✓ / ✗ / · with hover tooltip).
- [x] Field-coloured wireframe overlay in the viewport with a five-
      stop cool-to-warm divergent colour ramp.
- [x] Colour-bar legend (field name + min/max + ramp + timestep).
- [x] Time-series slider for transient runs (scrub through every
      VTU snapshot in the Field catalog).
- [x] Clickable field picker in the Results pane (switch which
      scalar drives the overlay).

## Acceptance checklist

- [ ] Keyboard navigation reaches every interactive control.
- [ ] Colour contrast AA on every surface/text pair.
- [ ] Screen-reader narration matches visible labels.
- [ ] No string hard-coded in Rust UI code (lint enforces).
- [ ] All three platforms produce signed, notarised installers in CI.

## Success metrics

| Metric                                      | Target          |
|---------------------------------------------|-----------------|
| `valenx --version` on a fresh install       | < 1 s           |
| Novice "hello airfoil" time-to-first-solve  | < 15 min        |
| Crash rate in telemetry                     | < 0.1 % / hour  |

## Leads into

[Phase 11 — HPC / cluster execution](./phase-11-hpc.md).
