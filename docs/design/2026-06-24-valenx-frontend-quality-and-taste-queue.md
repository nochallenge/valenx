# valenx Front-End Quality & Taste — build queue

**Date:** 2026-06-24
**Sources (user):** `thedaviddias/Front-End-Checklist` (run/apply it) and `Leonxlnx/taste-skill`
(AI-UI design-taste skills, MIT).

## Honest framing (read first)

Both inputs target **web** front-ends (HTML/CSS/JS). valenx's main UI is **native egui/wgpu — no HTML,
no CSS, no SEO/meta-tags**, so a literal application is mostly N/A and would produce fake findings.
They land in two honest ways:

1. **The `valenx-remote` served web page** (the phone-as-2nd-screen — a *real* HTML/JS surface): the
   checklist + taste guidance apply **directly** here.
2. **Translatable principles** for the native egui UI: accessibility, performance, responsive layout,
   security/input-validation, and *design taste* (hierarchy, spacing, restraint, consistent theming,
   intentional polish) — applied as the native analog, not as literal HTML items.

Front-End-Checklist categories (11): HTML(25)/CSS(32)/JS(26)/SEO(94)/Images(25) are web-page-specific
(~200 rules — N/A to egui, apply to the remote page); **Accessibility(95) / Performance(43) /
Security(22) / Testing(13) / Privacy(5) / i18n(5)** are the cross-cutting principles that translate.

## The track

### F1 — `valenx-remote` web-page audit + polish  (the checklist + taste apply directly)
Audit/harden the served page against the applicable items: `<meta charset=utf-8>` + viewport, semantic
HTML, **accessibility** (keyboard, ARIA, contrast, screen-reader), **performance** (the frame-poll
JPEG path, lazy/caching), **security** (the PIN gate, LAN-binding, HTTPS — required anyway for
`getUserMedia` in the scanning track — input validation, no secrets in the page), **privacy**. Then a
**taste pass** (taste-skill principles) on the page's look: type scale, spacing, restraint, a coherent
theme — it's the one genuinely-web surface valenx has.

### F2 — native egui UI-quality + taste pass  (translated principles)
A UI-quality sweep of the desktop app, framed natively:
- **Accessibility:** verify accesskit coverage — every interactive control has a clear accessible name
  (this *is* the AI-drivable-first gate), keyboard navigability, sufficient contrast in the theme.
- **Performance:** frame timing on heavy workbenches; the `rust-lld` build-speed win is already in.
- **Responsive/adaptive layout:** panels/workbenches resize cleanly; the tab/dock chrome holds up.
- **Taste:** consistent spacing/hierarchy, the painter-icon workbench chrome, a coherent theme +
  View→Text-size; apply taste-skill's "intentional, not generic" principles to egui idioms.
- **Testing/Privacy/i18n:** headless egui-logic tests for new panels; no PII in logs/feed; string
  externalization readiness (don't hard-block, just note).

Each finding becomes a small, scoped, gated commit (per AGENTS.md) — not one mega-PR.

### Note (separate from the build queue)
`taste-skill` can also be **installed as a Claude skill** (`SKILL.md`) so future UI-building turns
load its design guidance automatically — a tooling choice, independent of these build tasks.

## Loop integration
Lower-urgency than the connective-eval + MD/viewer tracks (it's polish, not new capability), but F1
(remote page) is self-contained and F2 dovetails with the standing AI-drivable-first gate. Scheduled
as slots free; F2 touches `valenx-app` (no new dep). Each chunk: scoped-gated green → email-safe
commit (leak 0) → review → next. GitHub HELD (local).
