# RFC 0008 — Reports and export

- **Status:** Draft
- **Authors:** @valenx-maintainers
- **Created:** 2026-04
- **Target phase:** 10

## Summary

Add a first-class **Report** concept: a reproducible, typeset
document assembled from a case's `Results` (plus supporting artifacts
and screenshots) that re-renders on demand from canonical data.
Reports are not markdown exports with embedded PNGs — they're live
templates whose contents update when a run is re-executed.

## Motivation

Engineers who use Valenx still need to hand something to a
colleague, regulator, or customer. The current options — screenshot
the viewport, paste into Word, hope the numbers don't drift — are
fragile. A native reporting surface lets the app own the entire
lifecycle: author once, re-run, re-render, archive.

## Design

### The Report document model

A report is a YAML / TOML document that references:

- One or more cases (by stable ID).
- A template (built-in or user-supplied) describing the structure:
  title, sections, captions, what fields / scalars to pull.
- Layout hints: page size, header / footer, table-of-contents
  behaviour, figure numbering scheme.
- Export targets: PDF, HTML, docx-compatible (ooxml), slide deck
  (ODP/PPTX-via-ooxml).

```toml
[report]
format = "1.0"
name   = "Airfoil drag report"
cases  = ["cfd-steady"]
theme  = "default"

[[report.section]]
title = "Summary"
body  = """
The {{case.cfd-steady.scalars.C_D | fmt("%.4f")}} drag coefficient was
computed using {{case.cfd-steady.header.solver}} with a
{{case.cfd-steady.flow.turbulence}} model.
"""

[[report.section]]
title = "Residual history"
figure = { kind = "residual-chart", case = "cfd-steady" }

[[report.section]]
title = "Pressure contour on surface"
figure = { kind = "viewport-snapshot", case = "cfd-steady",
            view = "iso", field = "p", range = "auto" }
```

### Typesetting

- Built-in rendering pipeline emits PDF directly (via a pure-Rust
  backend such as `typst-ref-parser` + `typst` bindings, or
  `printpdf` as a fallback).
- HTML output for archival in run bundles.
- docx / pptx paths go through a minimal OOXML writer specifically
  to avoid any LGPL OOXML libraries.

### Live figures

Figures are not cached PNGs — they're specs the renderer resolves
against the current `Results`:

- `residual-chart`: parameterised by case ID; pulls the `ScalarRecord`
  time series named `residual.*`.
- `viewport-snapshot`: the 3D view renders headlessly with a
  declared camera (named ViewCube direction or custom) against the
  current mesh + field.
- `contour-plot` / `vector-plot` / `streamline-plot`: 2D slices with
  declared orientation + offset.
- `table`: columns pulled from `ScalarCatalog`.

Every figure stamps its provenance (which run produced it) on the
output so users can verify the chain.

### Interaction with the run pipeline

- A report is just another artifact of a run: it lives under
  `results/<run-id>/report/`, with its source YAML and final PDF.
- Re-running the case rebuilds the report automatically.
- An export command (`File → Export Report`) produces a standalone
  PDF with the full data set embedded so recipients can open it
  offline.

### Default templates

Shipped with Valenx:

- **Quick summary** — title, one figure per physics, headline
  scalars.
- **Comparison** — side-by-side two cases with diff tables.
- **Study** — parameter sweep table with scatter plot + summary
  statistics.
- **Regulatory** — boilerplate structure for disclosures (loads,
  assumptions, standards cited, provenance chain).

### UI surface

- New left-panel node: *Reports*, sibling to *Cases* and *Results*.
- Template gallery on "new report".
- WYSIWYG-ish preview: the report renders live as the user edits
  captions; figures render on save because they're expensive.
- Version control: saving a report bumps its version and records the
  trigger (manual save vs post-run auto-refresh).

## Drawbacks

- Building a reporting engine is non-trivial — PDF generation, table
  layout, figure rendering, template evaluation.
- Users may still want to copy-paste into their own tooling; we must
  keep that path easy.

## Alternatives

- Just emit a markdown dump with PNGs. Works today, but doesn't
  regenerate when the run changes.
- Punt to Pandoc. Adds a build-time dep, licensing isn't our
  headache, but installer size and the "pandoc crashed on obscure
  input" failure mode land on us.

## Prior art

- Typst (pure-Rust typesetting, excellent output quality).
- Simulia Abaqus report generator (proprietary, but the model of
  declarative templates + live data is proven).
- R Markdown / Quarto for scientific reports — Valenx can learn from
  the data-binding pattern without adopting the Pandoc stack.

## Unresolved

- Template language: Liquid-style mustache, or a Rust-specific
  minimal engine? Leaning towards mustache-compatible so people can
  hand-edit outside the app.
- Offline HTML output bundles: do we embed everything or leave
  artifact paths relative? Leaning embed for archival use.
- Accessibility: tagged PDF output for screen readers is
  non-trivial — plan to ship in Phase 10.5 or later.
