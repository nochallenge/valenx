# Valenx Design

The design plan for Valenx — why we do things the way we do, what needs
to get designed, in what order, and how we know it's working.

See also:
- [ARCHITECTURE.md](./ARCHITECTURE.md) — how the code fits together
- [ROADMAP.md](./ROADMAP.md) — the 20-year direction this design serves
- [rfcs/](./rfcs/) — specs; this doc is the plan, RFCs are commitments

---

## 1. How to read this doc

This is the **plan** — living, opinionated, meant to guide decisions
and be argued with. It isn't a spec; specs live in RFCs. When a part
of this plan matures enough to commit, it turns into an RFC.

If the ROADMAP and this doc disagree, the ROADMAP is direction; this
doc is method. If this doc and an RFC disagree, the RFC wins — it's
been reviewed and merged; this file can drift.

Two realities this plan honors up front:

- **We are solo.** For at least Year 1, one person does design and
  code and everything else. Every validation / process line has a
  "cheap mode" for what we can actually afford right now.
- **We are long-term.** 20 years means this file should age well.
  Over-committing to specific tools or vendors is worse than
  under-committing and iterating.

---

## 2. Design principles — the non-negotiables

Five commitments. When trade-offs come up, these break ties.

**1. Polish before scope.** One physics done beautifully beats five
done roughly. We ship CFD to Fusion quality before we ship FEA
half-done. This is a deliberate order, not a rejection of breadth.

**2. Dense, not busy.** Simulation engineers live in numbers. Thirty
data points on screen is fine if they're organized; three data points
with wasted space is not. **Fusion 360** is the aesthetic reference —
dense, but never crowded.

**3. Two speeds, same tool.** A newcomer runs their first airfoil
case by following prompts. An expert does it in 12 keystrokes via
the command palette. Neither experience is bolted on; both are
primary. No "pro mode" toggle.

**4. The app is the product.** No modal upsells, no signup walls, no
"complete your profile." No telemetry by default. No marketing copy
inside the UI. The app earns attention by being good.

**5. Quiet.** Animations are subtle. The app does not make sounds.
Colors are restrained. The user's model and results are loud; the
app is quiet around them.

These get extracted into a shorter `DESIGN_PRINCIPLES.md` (< 1 page)
referenced from every UI PR.

---

## 3. What design owns (scope)

| Domain | Deliverables |
|---|---|
| Visual language | Color, typography, icons, illustration |
| Design system | Tokens → components → patterns, in code and mockups |
| Interaction | Keyboard model, mouse/trackpad conventions, drag/drop, gestures |
| Information architecture | Screen hierarchy, navigation, search, what lives where |
| Motion | Transitions, feedback timing, progress behavior |
| Data visualization | Viewport rendering style, plot conventions, result overlays |
| Accessibility | Keyboard-complete, screen-reader-friendly, high-contrast, colorblind-safe |
| Themes + localization | Dark / light / high-contrast; i18n infrastructure |
| Voice + copy | Microcopy, error messages, empty states, documentation tone |
| Onboarding + learning | Tours, tutorials, contextual help, progressive disclosure |
| Design process | How decisions are made, reviewed, and shipped |

Code owns the implementation; design owns the spec.

---

## 4. Reality check: two modes

Everything that costs money has a **cheap** default and a **funded**
upgrade. Lines marked **(funded)** are deferred until a grant,
sponsorship, or TSC-directed spend materializes.

| Activity | Cheap mode (default) | Funded mode |
|---|---|---|
| Accessibility audit | Self-check against WCAG list; NVDA smoke test | External audit, annually |
| Usability testing | Hallway tests with 3 friends per release | Recruited users, 5 per persona, 60-min sessions |
| Visual audit | Peer review inside team | Outside designer review, annually |
| Design tooling | SVG in the repo, free Figma / Penpot / Inkscape | Paid design-tool library seats |
| Fonts and icons | Open-licensed only | Licensed if warranted |

The rest of the plan assumes cheap mode until stated otherwise.

---

## 5. The design system — three layers

### Layer A — Tokens

Machine-readable primitives. Single source of truth: a JSON file at
`crates/valenx-design-tokens/tokens.json`, with a small build script
that emits:

- Rust `const`s for `egui` code
- An SVG color palette for mockups
- Optional JSON mirrors for any design tool that imports tokens

Token categories:

| Category | Examples |
|---|---|
| Color | surface-0..5, text-1..3, accent-primary / success / warning / error / info, physics-cfd / fea / em / chem / md / battery |
| Type | size-xs..3xl, weight-regular / medium / semibold, line heights, mono variants |
| Space | space-0..12 (4-pixel base) |
| Radius | radius-none / sm / md / lg / full |
| Shadow | elev-0..3 |
| Motion | duration-fast / base / slow, easing-standard / emphasized / decelerate |
| Border | border-hairline / default / strong |
| Z-index | z-base / docked / overlay / modal / tooltip |

Principle: tokens are named by role (`surface-1`), never by value
(`gray-900`). A theme swap changes the mapping; components don't know.

**Typography specifics.**

- **UI body** — **Inter**, SIL Open Font License. Weights: regular
  (400), medium (500), semibold (600).
- **Numerics, code, logs** — **JetBrains Mono**, Apache 2.0.
  Tabular figures forced on for aligned columns of numbers.
- Both fonts live in `crates/valenx-fonts/` and are **embedded in
  the binary**. No system-font fallback for core UI — avoids
  cross-OS layout drift and keeps rendering identical for pixel
  snapshot tests.
- Font sizes resolve from the type-size tokens; the whole scale is
  multiplied by a user-settable font-size factor (Settings →
  Appearance) for low-vision users.

### Layer B — Components

UI building blocks. Year 1 ships ~15; the full library (~40–50)
matures across Year 2–3.

**Year 1 set:**
- Buttons — primary, secondary, ghost, icon
- Inputs — text, number, dropdown, checkbox, slider
- Cards and panels
- Dialog and sheet
- Toast and snackbar
- Tree
- Tab-bar and rail
- Ribbon — tab, group, button
- Command palette
- Status dot / badge
- Progress bar and spinner

Each with states: default / hover / active / focus / disabled /
error / loading / selected. Documented via a visual index generated
from snapshot tests.

### Layer C — Patterns

Opinionated compositions for recurring problems.

**Year 1 patterns** (the shell needs these):
- Ribbon tab composition (how adapters contribute)
- Browser tree conventions
- Log viewer
- Dialog composition rules (modal vs. sheet vs. popover)
- Long-running operation feedback (progress + ETA + cancel)
- Validation + error surfacing
- Command palette invocation
- Units display + input

**Year 2+:** property panel, timeline, residual pane, expression
editing (once its RFC lands), results comparison, parametric sweep
view.

Each pattern gets `docs/design/patterns/<name>.md` with: problem,
when to use, when not to, worked example, accessibility notes.

---

## 6. Screen inventory

**Year 1 — minimum viable product.** Seven screens. Cut hard from
the original eighteen because a solo developer can't ship that many
at Fusion quality in a year.

1. Home — recent / templates / status / learn
2. Workspace shell — ribbon + browser tree + viewport + log
3. Command palette (overlay)
4. Project dialogs — new / open / save-as
5. Run overlay — launch + progress + cancel
6. Settings — appearance, shortcuts, tools
7. Welcome tour — first launch only, dismissible forever

**Year 2 — polish and breadth**
- Tool Manager
- Plugin Manager
- Residual / convergence pane
- Plot editor
- Tabular results view
- Export dialogs
- Keyboard shortcut cheat sheet
- In-app user manual

**Year 3–5 — depth**
- Sketcher + constraints panel
- Meshing workflow UI
- Material library
- Post-processor editor (streamlines, iso-surfaces, slices)
- Case diff viewer
- Validation gallery browser
- Parametric sweep view
- Results comparison (A/B)

**Year 5+ — scale**
- Plugin registry browser
- Collaboration / review mode
- Cloud compute setup
- Enterprise admin
- Classroom mode

Each screen gets three artifacts: mockup (SVG in repo), pattern
doc, and reference implementation in `valenx-app`.

---

## 7. Core interactions — commitments

**Keyboard is first-class.** Every action has a keyboard path.
Everything reachable in ≤3 keystrokes through the command palette.

**Shortcut assignment convention.** Not every action gets a direct
shortcut — too many to remember. Direct shortcuts are reserved for:

1. **OS-standard actions** — follow the host platform verbatim
   (`Cmd/Ctrl+S`, `Cmd/Ctrl+Z`, `Cmd/Ctrl+C/V`, `Cmd/Ctrl+W`).
2. **Top-level navigation** — open project (`Cmd/Ctrl+O`), new
   project (`Cmd/Ctrl+N`), command palette (`Cmd/Ctrl+P`),
   settings (`Cmd/Ctrl+,`), toggle log (`Cmd/Ctrl+\``), help
   (`F1`).
3. **Frequent workflow actions** — run (`Cmd/Ctrl+R`), save as
   (`Cmd/Ctrl+Shift+S`), toggle ribbon (`Cmd/Ctrl+F1`).
4. **Single-key verbs in the viewport** — Fusion-style when
   focus is on the 3-D canvas: `E` extrude, `F` fit, `S` sketch,
   `M` measure. Ambiguous ones (`S` for sketch vs. save) are
   context-scoped, not global.

Everything else flows through the command palette. Users rebind
anything via Settings → Keyboard. Shortcut hints appear in tooltips
after >500 ms hover (not for everyone, not always).

**Keyboard-layout handling.** Shortcuts bind to **physical keys**,
not characters — so `Cmd+Z` on AZERTY stays on the bottom-left
corner of the keyboard, where English-layout users' muscle memory
lives. Users on alternative layouts (Dvorak, Colemak) can opt for
character-bound shortcuts in Settings if they prefer.

**File dialogs.** OS-native via the `rfd` crate. Familiar look,
right drag-and-drop integration, remembers the last directory per
OS convention. No custom in-app file chooser until a specific use
case justifies one (and that goes through an RFC).

**Undo has a defined boundary.**
- *In the undo stack:* project-model edits, case edits, timeline
  actions, text edits in persistent fields
- *Not in the undo stack:* viewport camera (has a separate camera
  history), transient UI state (panel layout), running a solver,
  script execution, file I/O
- One global stack for what's in it; `Ctrl/Cmd+Z` rolls back the
  last *project* action regardless of which pane it came from.
- **UI feedback on undo:** brief toast bottom-right — *"Undid: Set
  inlet velocity"* with a Redo link; auto-dismiss at 2 s. No
  animation, no sound. Stack depth shown discreetly in the status
  bar (`↶ 24 ↷ 0`) for users who want it.

**Inline editing over dialogs.** Click a number, edit, tab away.
Dialogs are for operations that can't be single-value edits.

**Direct manipulation in the viewport.** Drag a BC handle; selection
is first-class.

**Save is implicit + explicit.** Auto-save project state every 30 s
to a crash-safe side file; `Ctrl/Cmd+S` writes to the canonical file.
Viewport / panel state: 5 s cadence.

**Log is always one keystroke away.** `Ctrl/Cmd+\`` toggles.

**Background runs are first-class.** Start a CFD run, immediately
open another project; the run keeps going; a status chip in the rail
reports progress. (Engineering cost acknowledged in Section 26.)

**Right-click does something useful.** Context menus are where power
lives.

**Gestures need keyboard equivalents,** *except* for industry
standards: pinch-to-zoom, two-finger pan, trackpad rotate — these
are discoverable from muscle memory and don't need a key.

**Expressions in numeric fields** — originally listed here as a
commitment, now moved to its own RFC. The parser, safety model, units
system, and scoping work is real engineering that wasn't budgeted in
Year 1. Target: Phase 2+.

---

## 8. Viewport design

The viewport is where most user time goes. Pinning enough now that
it stops being "TBD." Full spec lands as its own RFC in Phase 0.

**Camera model.**
- Orbit around a *pivot point*. Pivot is user-settable (middle-click
  anywhere to re-pivot); defaults to scene center on reset.
- Roll disabled by default (prevents disorientation); togglable in
  Settings for power users.
- Zoom anchors at the cursor, not the pivot — matches Fusion and
  Rhino convention.
- Camera state persists per-project via the `.valenx` bundle.

**Selection.**
- Single-click selects the topmost visible face / edge / body per
  the selection-filter setting in the status bar (All / Face / Edge
  / Vertex / Body).
- `Shift+click` adds; `Ctrl/Cmd+click` toggles.
- **Box-select** — drag a rectangle on empty canvas; selects all
  entities whose bounding box intersects. Left-to-right = fully
  contained; right-to-left = intersect (matching AutoCAD / Fusion
  convention).
- **Lasso-select** — `Shift+drag` on empty canvas draws a lasso;
  selects entities inside the closed shape.
- **Tree ↔ viewport sync** — selecting a node in the browser tree
  highlights in the viewport; the tree scrolls and expands to
  show the selected node when the viewport selection changes.
  Two-way mirror, never drifts.
- Highlight: 60 %-opacity accent overlay plus a darker edge stroke.
- Persists through rotation; clears on `Esc` or click-on-empty.

**Manipulators** (Year 2):
- Translation / rotation gizmos when an editable transform is
  selected.
- Snap to grid / axis / existing geometry, toggled in the status bar.

**Overlays.**
- *Year 1:* coordinate axes (bottom-left), ground grid, ViewCube
  (top-right).
- *Year 2:* measure tool, section planes.
- *Year 3+:* dimension annotations, callout labels.

**Rendering styles.**
- Flat-shaded with edges (default)
- Shaded only (clean renders)
- Hidden-line
- X-ray (50 % transparency for picking through geometry)
- Section view (when section planes ship)

**Performance budget.**
- 30 fps minimum at 1 M triangles on mid-range integrated GPU
- 60 fps target at 1 M on discrete GPU
- LOD + GPU-resident meshes for 10 M+

---

## 9. Data visualization

The other surface engineers live in: plots and tables. Full RFC in
Phase 1.

**Plots.**
- Library: build on `egui_plot` for 2-D; extend where needed;
  own code only where the library blocks us.
- Aesthetic: close to matplotlib's `seaborn-whitegrid` — minimal
  grid, thin axes, no chart junk.
- Defaults: drag to pan, wheel to zoom on the axis under cursor,
  hover crosshair + values, click-drag region zooms to it, `R`
  resets view.
- Color cycling: colorblind-safe **Okabe-Ito** palette (8 hues,
  deuteranopia / protanopia / tritanopia friendly).
- Typography uses the same token family as the rest of the app.
- Export: SVG and PNG at arbitrary DPI; paper-publishable by default.

**Residual / convergence plot.**
- Log-y by default; switchable to linear.
- Multiple series, colored by field (U, p, k, ω, …).
- Vertical gridline at current iteration.
- Wall-time on a secondary axis.

**Tables.**
- Virtualized scrolling for large results.
- Sortable, filterable, CSV-copyable.
- Numeric columns right-aligned, monospace.
- Sparkline columns for per-row trends.

**Viewport overlays for results.**
- Contour maps; adjustable colormap (viridis, coolwarm, grayscale,
  custom).
- Legend floating, repositionable, unit-aware.
- Vector fields (glyphs / LIC).
- Iso-surfaces.
- Streamlines (Year 2+).

**Colour and monitor calibration.** Plots and contours target
**sRGB** rendering for predictability — `wgpu`'s output surface is
configured for sRGB and all tokens / palettes are authored in
sRGB. We don't attempt ICC-profile-aware colour management in v1
(good enough for 99% of engineering use, too much complexity for
1% who need certified colour). Users who need publication-grade
colour fidelity should export to SVG / PDF via `typst` and let
their publishing pipeline handle colour management.

---

## 10. Error and recovery UX

Errors happen constantly in simulation. Design them first-class, not
as afterthoughts. Taxonomy aligned with `AdapterError` in RFC 0002.

| Error | Where surfaced | Remedy hint |
|---|---|---|
| ToolNotInstalled | Modal with "Install now" | Opens Tool Manager |
| ToolVersionMismatch | Inline banner at the case | "Update tools.lock?" / "Use override" |
| InvalidCase | Inline red ring at the offending field | Tooltip explains why |
| Translate | Dialog on run-start | Shows offending parameter + adapter-docs link |
| Run (exit ≠ 0) | Log viewer, auto-scrolled to last error | "Copy error", "Open workdir" buttons |
| ParseOutput | Warning banner in Results | Partial results still loadable |
| Cancelled | Status chip in the rail | Silent — user knew |
| IO | Toast | Generic retry |

**Principles:**
- Never a bare "something went wrong." Every error says what failed
  and where.
- Solver crashes don't crash the app — subprocess isolation means we
  show the log and let the user retry.
- Recovery is always one click from the error surface.
- Messages are in the user's voice, not the developer's. "The mesh
  has 0 cells" — not "mesh.cell_count == 0 assertion failed."

**Corrupt project file on load.** If a `.valenx` directory is
present but `project.toml` fails to parse, the app shows a
recovery dialog:

- Diagnostic — exact line / column of the parse error
- Options: *Open read-only* (load what we can; refuse to save
  over the original), *Restore from last auto-save*, *Show in
  Explorer / Finder*, *Cancel*
- "Open read-only" is the default and never destroys data

Never silently skip corrupt fields — partial-load paths are a
source of downstream bugs. Always explicit.

**Crash recovery.** If the app detects an unexpected exit on
next launch (a lock file left behind by the previous session),
the Home screen shows a banner:

> Valenx closed unexpectedly. Your last auto-save is 2 minutes
> ago. *Restore session* · *Discard*

Restore opens the last auto-saved project(s) at the viewport /
tab state they had. Discard dismisses the banner and carries
on.

**"Nothing works" home state.** If every adapter's `probe()`
returns `ToolNotInstalled` or `Broken` — e.g., user hasn't
completed the first-run wizard, or a system update broke the
tool directory — the status strip on Home turns red:

> ✗ No solvers available. *Run first-run setup* · *Open Tool
> Manager*

The New Project and Open Project actions remain enabled — users
can still look at existing projects even if they can't run
anything. The UI degrades gracefully, not catastrophically.

---

## 11. Documentation, microcopy, tutorials

Design owns the words.

**Microcopy.**
- Direct, calm, short. "Run" — not "Execute solver."
- Second person. "Your mesh has 2.3 M cells."
- Exclamation points sparingly and only with intent — never for
  emphasis or default cheer.
- No emoji inside the app — they're for casual chat, not
  simulation tools. Exception: unicode physics glyphs where they
  carry meaning.
- Numbers first. "50 m/s inlet velocity" — not "inlet velocity:
  50 m/s."

**Documentation structure.**
- *In-app:* `F1` opens a context-sensitive manual pane.
- *Shipped:* mdBook at `docs/src/`, embedded in the binary.
- *External:* the public website mirrors the same mdBook.

**Help discovery beyond F1.**
- **Help menu** in the title bar / application menu: *Getting
  started*, *User manual*, *Keyboard shortcuts*, *Validation
  gallery*, *Community*, *Check for updates*, *About*.
- **Command palette** indexes every help page; typing "turbulence"
  finds both the feature and the manual chapter about it.
- **Tooltips with "Learn more" links** for moderately complex
  concepts (mesh refinement, turbulence models, BC types) —
  tooltip shows a one-sentence definition plus a link that opens
  the manual pane to the right chapter.
- **Contextual "?" affordances** in dense panels (property
  editors, solver settings) sit next to parameters and open the
  relevant manual section.

**Tutorials.**
- Every tutorial has three phases: setup (pre-made project), task
  (things the user does), discussion (why it worked, what to tweak).
- Tutorials live in the Samples tab on Home and open as real projects
  with a side panel of instructions that advance on actions.

**Empty states.**
- Every empty list has a short sentence plus one affordance telling
  the user how to fill it.
- No decorative illustrations — one icon is enough.

---

## 12. Novice-to-expert journey

Sim users come from two directions: physics-strong / tool-weak
(researchers) and tool-strong / physics-weak (engineering
generalists). Design supports both.

**First session.**
- Welcome tour optional, dismissible forever.
- Open a template; all parameters named in the user's vocabulary,
  not the solver's internal terms.
- One-click "Run on sample data" before they author anything.

**First ten sessions.**
- Breadcrumb suggestions surface as they explore ("Try a
  mesh-independence check?").
- Command palette learns what they've searched; promotes frequent
  commands. History stored locally at
  `~/.valenx/command-history.toml`, never transmitted.
- Tutorials adapt: skip sections they've clearly mastered.

**After ten sessions.**
- UI stops suggesting; user owns the learning.
- Power features (palette macros, expressions, scripting) become
  visible in menus without ever being pushed.
- Keyboard-shortcut hints appear in tooltips only after a >500 ms
  hover. We don't teach experts what they already know.

**Past the hundredth session.** Experts benefit from:
- **Workspace layouts** — save named panel arrangements per
  workflow (CFD run, FEA modal, debugging a failing solve),
  switch between them instantly
- **Command macros** — sequences recorded or scripted, bound to
  shortcuts; runs N actions on a keystroke
- **Parameter sets** — saved configurations of solver / case
  settings, reusable across projects
- **Custom templates** — promote a project to a template the
  user reuses; optionally shareable with teammates

None of these are discoverable on first launch — they appear in
menus once the user's palette-history or usage-count reaches a
threshold. Experts find them; novices never feel them.

No difficulty setting; the journey is implicit in interaction
history, which is local and never leaves the machine.

---

## 13. Iconography

Custom-drawn, one coherent set, SVG source at 24 × 24.

- Three weights: regular / bold / filled.
- Categories: physics, tool, action, status.
- Shipped via `crates/valenx-icons/` as embedded assets with
  type-safe names generated from the SVG inventory.
- Contribution: SVG files in the repo, PR-reviewed like code. A
  design-steward role (Section 24) reviews icon PRs for visual
  consistency.

**Icon inventory.** Full inventory is a Phase 0 deliverable,
scheduled after Year-1 mockups stabilize — it's the list of every
distinct role on every committed screen, no fabrication.

Pre-inventory rough order of magnitude: tens at Year 1, low
hundreds by Year 5. Firm numbers land with the inventory PR, not
in this doc.

---

## 14. Motion and sound

### Motion
Subtle, consistent, purposeful. Motion explains, never decorates.
Respect `prefers-reduced-motion` — cut durations 75 % when set.

| Situation | Duration | Easing |
|---|---|---|
| Hover / focus | 120 ms | standard |
| Menu open / close | 150 ms | decelerate |
| Dialog / sheet enter | 200 ms | emphasized |
| View transition | 180 ms | standard |
| Panel resize | 250 ms | standard |
| Viewport camera tween | 400 ms | emphasized |
| Error shake | 300 ms | custom spring |

What does **not** animate: text in lists, live plots, log lines, any
value being actively edited.

### Sound
**None, ever.** Not on startup. Not on errors. Not on completion.

- Saves disk
- No localization complication
- Welcome in open-plan offices and libraries by default
- Deaf and hard-of-hearing users aren't disadvantaged

If users later want completion chimes, that's a sound-system RFC
starting from zero. The default is silence.

---

## 15. Themes, accessibility, localization

### Themes
- **Dark** (default) — near-black background, cool accents
- **Light** — warm off-white, same accent roles
- **High-contrast** — WCAG AAA
- **Colorblind-safe variant** — accent remapping for deuteranopia /
  protanopia

### Accessibility
- **Target:** WCAG 2.1 AA. Cheap-mode verification covers what
  automated tools and checklists can see — colour contrast,
  keyboard order, focus rings, label presence, heading structure.
  Full AA certification requires a funded-mode external audit.
- Every action keyboard-reachable
- Focus ring visible everywhere
- Screen-reader support via `egui`'s AccessKit integration
- *Cheap mode:* self-audit against WCAG checklist every release;
  NVDA / VoiceOver smoke test
- *Funded:* external audit yearly; usability testing with disabled
  users

### Localization
- English only at 1.0
- Infrastructure from day one — all strings externalized via a
  single Fluent file (`strings.ftl`); no baked-in English
- Post-1.0 languages in priority order of contributor availability
- RTL deferred until a contributor commits to Arabic or Hebrew

---

## 16. Branding identity

A 20-year project has a brand, whether you design it or not.
Design it.

**Elements.**
- **Logomark** — wordmark plus symbol. Must read at 16 px (taskbar
  icon) and at 512 px (about dialog). First draft in Phase 0,
  formalized in Phase 1.
- **App icon** — platform-native bundles: `.ico` on Windows, `.icns`
  on macOS, `.png` set for AppImage / `.deb` / `.rpm`. Respects OS
  convention (squircle on macOS, flatter on Windows 11).
- **Accent colour family** — one primary accent, mapped from
  `accent-primary` token; status colours from success / warning /
  error / info. Changes only with a design RFC.
- **Splash screen** — none. Cold-start straight to the main window
  is the launch experience.
- **Favicon / docs branding** — mirrors the app icon.

**Voice and personality.** Calm, technical, respectful of the user's
time. Not playful. Not academic-stiff either. Valenx is a tool —
confident, direct, assumes the user is an adult engineer.

**Marketing / website.** Explicitly out of scope of this doc. When
a website stands up, it's a sibling effort that uses the same
tokens, icons, and fonts, and is reviewed by the same stewards.

---

## 17. Settings

The surface users tune to make the app theirs. Designed, not just
dumped into a panel.

**Panes at 1.0** (Year 1 minimum):
- **Appearance** — theme, font-size factor, motion preference
  (including `prefers-reduced-motion` respect toggle)
- **Keyboard & shortcuts** — full shortcut map, searchable,
  per-action rebindable, conflict detection
- **Tools** — installed solvers, versions, update channel per
  tool, user overrides, add / remove (the Tool Manager)
- **Units** — default unit system (SI / Imperial / custom);
  per-physics overrides
- **General** — auto-save interval (default 30 s), logging
  verbosity, working directory default, update channel

**Panes at Year 2+:**
- **Plugins** — installed plugins, permissions, registry access
- **Scripting** — Python / Lua runtime selection, module path
- **Network** — proxy, offline mode, registry endpoints, per-
  category opt-in for user-initiated remote calls
- **Advanced** — experimental flags, developer tools toggle,
  tracing level

**Persistence location** (OS conventions):

| Platform | Path |
|---|---|
| Linux | `$XDG_CONFIG_HOME/valenx/settings.toml` (default `~/.config/valenx/`) |
| macOS | `~/Library/Application Support/Valenx/settings.toml` |
| Windows | `%APPDATA%\Valenx\settings.toml` |

**Secrets** (API keys, cloud credentials, registry tokens) go in
the OS keychain — Secret Service / Keychain / Credential Vault —
never plaintext on disk.

**Sync across machines** is explicitly **not a feature at 1.0**.
Settings stay local. Users who want sync can put the config
directory under `git` or `syncthing` themselves.

---

## 18. Scripting UX

Valenx embeds Python (pyo3) and Lua (mlua) per LANGUAGES.md and
ARCHITECTURE.md. Design positions them as follows.

**Surfaces.**
- **Script Console** — dockable pane. Default shortcut
  `Cmd/Ctrl+Shift+\``; rebindable. REPL on the embedded
  interpreter, full access to the `valenx.scripting` module, the
  current project, and the result cache.
- **Script editor** — edit a `.py` or `.lua` file inside the app;
  syntax-highlight via the same tokens; run against the current
  project from a play button or shortcut.
- **Project scripts** — scripts inside a `.valenx` project under
  `scripts/`; surfaced in the browser tree and runnable from the
  Scripts ribbon tab.

**Trust model.** Scripts run with the capability model of the
plugin API (RFC 0003). By default: read-only in the project;
explicit grants for file writes outside the project, network,
subprocess spawn. Prompted on first attempt; saved per project
or globally.

**Language picker.** Per-project setting. If you don't want Python,
it's off; never initialized, no interpreter start-up cost.

**Collaboration.** Scripts travel with the project — a shared
`.valenx` includes them. Running a received project's script
requires explicit user approval on first run.

---

## 19. Reports and export

Simulation produces artifacts. Design the export story.

**Export targets (Year 1):**
- Plot → SVG, PNG (arbitrary DPI)
- Table → CSV, TSV, clipboard
- Mesh + fields → VTK (legacy + XML formats)
- Full report → PDF, rendered via **typst** (fast, reproducible,
  LaTeX-quality output, Apache 2.0)

**Export targets (Year 2+):**
- Word `.docx` for corporate workflows
- LaTeX source for thesis / paper writers
- Native solver input deck (for verifying the translation layer
  against the upstream tool — debug feature)
- Animation (MP4 / WebM) for transient cases

**Template-based reports.** A built-in set for each physics type:
CFD convergence study, FEA linear-static summary, modal analysis,
validation-case write-up. Templates are Markdown-with-placeholders
files rendered against project data. Users can author their own
templates in `~/.valenx/templates/` or inside a project's
`reports/` folder.

**Headless rendering.** `valenx render project.valenx --report
convergence --out pdf` produces a report without opening the UI.
Lets CI build nightly validation reports; lets organizations run
Valenx on a headless server for batch reporting.

**No "Share" button.** Export writes a file; sharing is the user's
workflow (email, Slack, Drive, whatever). We don't bake in an
integration that puts their data on any particular server.

---

## 20. App update UX

**Check is manual, not automatic.** `Help → Check for updates`.
The app does not background-poll.

**Channels.** Stable, LTS, Nightly. User picks in Settings →
General. Default: Stable.

**Update mechanism.** Signed packages per platform:
- **Windows** — MSIX with update feed, or the installer re-run
- **macOS** — Sparkle-compatible feed, with notarization
- **Linux** — distro packages (`.deb`, `.rpm`) and `AppImageUpdate`
  for AppImage

**Release notes inline.** Before the user accepts an update, they
see what changed — bullet list from the CHANGELOG, with a link to
the full entry. No update proceeds on a click-through; the user
reads what they're getting.

**Rollback.** The previous installed version is kept for one major
release so `Help → Roll back` works if a new version misbehaves.
Two majors back, the rollback bundle is garbage-collected.

**Offline installs always work.** Users can always download the
installer for any supported platform from the releases page and
run it manually. No online check required to install.

---

## 21. Multi-document and multi-window

**One instance, many projects.** Default model: tabs at the top of
the window, one per open project. Familiar from VS Code and Fusion.

- `Cmd/Ctrl+N` new project in a new tab
- `Cmd/Ctrl+O` open existing project in a new tab
- `Cmd/Ctrl+W` close current tab
- `Cmd/Ctrl+Tab` / `Cmd/Ctrl+Shift+Tab` cycle tabs

**Tab detach.** Drag a tab off the tab bar and it pops into its own
OS window. Same app process, same registry, same adapter state —
just a second viewport and a second set of panes for multi-monitor
users.

**Multiple app instances.** Launching Valenx twice is allowed but
discouraged — each instance has its own settings cache, no shared
state, and concurrent access to the same `.valenx` project is
undefined. Power users only; not an officially supported workflow.

**Background solver runs** persist across tab switches and window
changes. A tab whose solver is running shows a status dot in its
tab header. Closing that tab prompts:
- *Keep running in background* (default)
- *Pause and keep state*
- *Cancel the run*

**State isolation.** Each project has its own undo stack, its own
camera history, its own timeline. Switching tabs is a clean
context switch; nothing leaks across.

---

## 22. Mockup and source-of-truth tooling

Decision: **no proprietary tool is a dependency**.

- **Canonical source:** SVG mockups in `docs/design/mockups/`,
  versioned with the code. Edited in any SVG editor.
- **Working tools:** Figma, Penpot, Inkscape — use what you like.
  Export to SVG and commit. Figma is not required and is not the
  source of truth.
- **Tokens:** JSON at `crates/valenx-design-tokens/tokens.json`,
  mirrored to Rust `const`s and optional design-tool palettes.
- **Sharing:** public. Anyone can browse `docs/design/` on GitHub.

A BDFL-or-maintainer can use whichever tool they want for a personal
draft. When a draft becomes a design decision, it lands as SVG +
pattern doc in the repo, reviewed in a PR.

---

## 23. How design ships — process

Workflow for a UI change:

```
1. Problem sketch         (issue or discussion)
2. Mockup                 (SVG committed, or any tool's export)
3. Design RFC             (if pattern-level; otherwise skip)
4. Implementation PR      (Rust + snapshot tests + docs)
5. Merge
```

Design RFCs follow the normal RFC process
([rfcs/README.md](./rfcs/README.md)) — minimum 10 days open,
consensus-seeking, BDFL tie-break.

**PR checklist for UI changes** (lives in `PULL_REQUEST_TEMPLATE.md`):
- Mockup link or screenshot
- Snapshot test for the new / changed state
- Tokens used — no hard-coded colors, spacings, motions
- Keyboard-reachable
- String key used — no hard-coded English

**Visual regression.**
- *Structural / tree-based* assertions via `egui_kittest` for the
  bulk of components — not flaky across platforms.
- *Pixel diffs* only for a small "hero" set (Home thumbnail, ribbon
  header, command palette), and only on Linux CI. Reason: our
  Linux runners use a pinned FreeType plus our bundled Inter /
  JetBrains Mono, so the same draw call produces the same pixels
  run after run. Windows (ClearType, DPI-variable) and macOS
  (Apple's rasterizer) render the same call differently on the same
  OS across point releases; pixel diffs there are too noisy to be
  useful as gates.

---

## 24. Contribution flow and bus factor

### Contributing to design
- *Mockups:* PR with SVG plus pattern doc
- *Icons:* PR with SVG plus entry in the icons index
- *Tokens:* PR with JSON change plus rationale (usually an RFC)
- *Copy / translations:* PR to `strings.ftl` or a locale file
- Everything reviewed by the design-steward role

### Attribution
- Community contributions are credited via Git history — no
  separate CONTRIBUTORS file needed for design work; `git blame`
  and `git log` are the canonical record.
- Icons, mockups, and templates contributed by community members
  may carry a short `author` field in their source file (`<!--
  author: handle -->` in SVG, `author = "handle"` in template
  TOML) when the contributor wants visible credit.
- **We do not ask contributors to assign copyright** — Apache 2.0
  contributor-license terms apply (see [POLICIES.md](./POLICIES.md)).
- The MAINTAINERS.md file lists design-stewards and emeritus
  stewards by year.

### Bus factor
- **Year 1–2:** BDFL is design-steward. Acknowledged single point
  of failure. Mitigations: this document, mockups versioned in the
  repo, no proprietary tool lock-in.
- **Year 2+:** As soon as a contributor has landed 5+ quality
  design PRs, promote them to co-steward. Two-person design
  stewardship becomes the minimum.
- **Year 5+:** Formal design-steward role with TSC membership;
  mirrors the code-maintainer rule of ≥2 per subsystem.

---

## 25. Telemetry and network silence

**No telemetry. No automatic server reach.**

The app does not phone home for:
- Usage metrics
- Crash reports (by default — see opt-in note below)
- Feature-flag experimentation
- A/B tests
- Marketing
- "Check-in" pings of any kind

**User-initiated network is a different category and is fine.**
When the user clicks *Check for Updates*, opens the Plugin Registry,
downloads a tool from the first-run wizard, exports to a cloud
endpoint they chose, or triggers any explicit remote action — we
make the request. Every such request:

- Is initiated by a visible user action
- Surfaces a status indicator while it runs
- Is logged in the activity log
- Can be disabled per-category in Settings → Network

If the TSC ever wants *opt-in* crash reporting — user explicitly
enables it in Settings and can inspect the payload before send —
that's an RFC. Until then, the default is silence.

---

## 26. Engineering prerequisites

Design depends on real engineering. Budgeting here so it doesn't
disappear from the roadmap.

| Infrastructure | Rough scope | Phase |
|---|---|---|
| Token pipeline (JSON → Rust consts + SVG palette) | ~1 week | 0 |
| Icon pipeline (SVG → embedded, typed names) | ~1 week | 0 |
| Snapshot-test harness (`egui_kittest` wiring) | ~1 week | 1 |
| String externalization (Fluent + loader) | ~2 weeks | 1 |
| Theme system (token binding + hot-swap) | ~1 week | 1 |
| Command-palette infrastructure | ~2 weeks | 1 |
| Expression evaluator | uncertain — 6–10+ weeks realistic | 2+, via its own RFC |
| Viewport LOD + culling | 3+ weeks | 2+ |

Each gets a tracking issue in the relevant phase. Missing any blocks
the design work that builds on it.

---

## 27. Design RFC queue

Numbers get assigned at merge time. The sequencing below is *intent*,
not commitment, and will shuffle.

- **Design principles** — short RFC (< 1 page) locking the five
  non-negotiables
- **Token system and pipeline** — separate, deeper RFC covering
  JSON schema, Rust-const generation, SVG-palette export, theme
  binding
- **Results and fields data model** — affects plotting, tables,
  export
- **Ribbon composition protocol** — how adapters contribute tabs
- **Command palette protocol** — search infrastructure
- **Icon system** — SVG pipeline, naming, contribution
- **Viewport interaction model** — selection, manipulators, gestures
- **Plot / chart system** — charting library + aesthetic + API
- **Expression evaluation** — deferred from Year 1; parser, safety,
  units, scoping

---

## 28. Timeline — design deliverables per phase

Aligned with ROADMAP phases.

**Phase 0 (months 0–6) — Foundation.**
- `DESIGN_PRINCIPLES.md` extracted (< 1 page)
- Icon inventory PR — one row per icon on the committed Year-1
  screens; counts fall out of it
- Tokens JSON + Rust pipeline working
- Icon pipeline + first icon set (count = inventory)
- Mockups: Home, Workspace shell, command palette (SVGs in repo)
- Two design RFCs opened, discussed, accepted: one for principles
  (short), one for the token schema + pipeline (deeper)

**Phase 1 (months 3–18) — Shell + adapters.**
- Workspace shell implemented to match mockups
- Home view implemented
- Ribbon composition RFC accepted and implemented
- Year-1 component set shipped with snapshot tests
- Dark theme finished; light theme draft
- Keyboard-shortcut system shipped
- String externalization infrastructure
- *Cheap-mode* accessibility audit (self-check + hallway testing)
- *Funded:* first external accessibility audit, if budget allows

**Phase 2–3 (months 9–36) — Geometry + meshing.**
- Sketcher + constraints panel patterns
- Mesh preview pattern
- Material library UI
- Light theme finished; high-contrast theme shipped
- Second self-audit pass

**Phase 4–5 (months 18–72) — CFD basics + advanced.**
- Residual pane finalized
- Plot editor
- Parametric sweep UI
- Results comparison

**Phase 6–10 (years 2–12) — all physics + multi-physics.**
- Physics-specific ribbon tabs per vertical
- Localization to top 3 non-English languages
- Design system 1.0 — documented, hardened for plugin authors
- Plugin UI conventions RFC

**Phase 11–15 (years 7–20) — maturity.**
- Enterprise / classroom modes
- Cloud-compute setup UX
- Design-system retrospective at Year 10

---

## 29. Success criteria

How we know the design is working per phase. Floors we have to
clear, not KPIs to optimize.

**Year 1**
- First-run to first solved case: < 15 minutes for a novice,
  measured with 3 hallway-test subjects
- Home → running a simpleFoam airfoil: ≤ 8 clicks or ≤ 15 keystrokes
- All Year-1 actions reachable through the command palette in ≤ 3
  keystrokes
- Zero hard-coded colors, spacings, strings in the shell code
- WCAG AA self-audit passes on all Year-1 screens

**Year 3**
- Five complete physics workflows (CFD steady, CFD transient,
  linear-static FEA, modal FEA, antenna FDTD), each 8-click
  reachable
- *Funded:* external accessibility audit passes at AA on all screens
- Three non-English translations at ≥ 80 % string coverage
- 100+ community-contributed icons / mockups / copy edits merged

**Year 5**
- Adoption in ≥ 2 academic programs for teaching
- Design system 1.0 stabilized — breaking changes < 1/year
- Plugin authors shipping UI through the documented protocol
- Two-person design stewardship in place

**Year 10**
- Recognizable visual brand ("Valenx looks like Valenx")
- Certification-oriented workflows (DO-178C, ISO) designed and shipped
- Localization to ≥ 8 languages

---

## 30. Validation

**Cheap mode (always):**
- `egui_kittest` snapshot tests every PR
- Token lint / no-hard-coded-values lint every PR
- Accessibility smoke (focus order, labels) every PR
- Hallway test with 3 friends per release
- Keyboard-only walkthrough of core workflows per release

**Funded mode (when budget allows):**
- External accessibility audit, annually
- Recruited user testing — 5 per persona, 60-min sessions, recorded
  and analyzed
- Outside designer for visual / brand audit, annually

Neither mode uses telemetry. Validation is deliberately effortful,
not automated surveillance.

---

## 31. First actions — what the plan kicks off

### Already landed

- ✅ `DESIGN_PRINCIPLES.md` extracted (< 1 page)
- ✅ RFC 0005 (principles) and RFC 0006 (token schema + pipeline)
  merged
- ✅ `crates/valenx-design-tokens/` scaffolded with `tokens.json`,
  JSON schema, and build-time Rust codegen
- ✅ UI-PR checklist added to `.github/PULL_REQUEST_TEMPLATE.md`
- ✅ Icon inventory stubbed at `docs/design/icon-inventory.md`

### Still open — Phase 0 / early Phase 1

1. Commit the first mockup SVGs (Home + Workspace shell + command
   palette) to `docs/design/mockups/`
2. Draw the first ~50 icons listed in
   [`docs/design/icon-inventory.md`](./docs/design/icon-inventory.md)
   and drop them into `crates/valenx-icons/assets/`
3. Build a tokens-lint that fails CI on hard-coded colours / spacings
   / durations in UI code (Section 17 promises this; today it's
   human-enforced only)
4. Open tracking issues for the engineering prerequisites in
   Section 26 (token pipeline is done; icon pipeline, snapshot-test
   harness, string externalization, theme hot-swap, command-palette
   infrastructure remain)
5. Author the first pattern doc under `docs/design/patterns/`
   (ribbon-tab composition, since it's the first thing an adapter
   author will touch)

All Phase 0 work; doable in weeks.

---

Updates to this doc go through the RFC process
([rfcs/README.md](./rfcs/README.md)). Drift happens; keep the file
honest.
