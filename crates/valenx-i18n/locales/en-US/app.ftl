# Valenx en-US baseline string catalogue.
#
# Format: `key = value`, one per line. `#` starts a comment.
# Multi-line values are not supported in v0.1.0 — split into
# multiple keys if you need a long block.
#
# Keys use dot-separated namespaces:
#   ribbon.<button>          — ribbon labels
#   browser.<element>        — case-browser tree elements
#   dialog.<dialog>.<field>  — dialog labels
#   menu.<menu>.<item>       — menu items
#   palette.<command>        — command-palette entries
#   status.<state>           — status-bar messages
#   tooltip.<element>        — tooltip text
#   error.<code>             — user-facing error messages
#
# When adding a new key, add the same key to every other locale
# file in this directory or the loader will fall back to the
# en-US value at lookup time. The pseudo-locale build (when we
# add it) is generated programmatically — no .ftl file required.

# ---------------------------------------------------------------
# Ribbon
# ---------------------------------------------------------------
ribbon.run = Run
ribbon.prepare = Prepare
ribbon.cancel = Cancel
ribbon.open-project = Open project
ribbon.import-stl = Import STL
ribbon.import-mesh = Load canonical mesh
ribbon.export-html = Export HTML report
ribbon.export-csv = Export scalars CSV

# ---------------------------------------------------------------
# Browser
# ---------------------------------------------------------------
browser.section.project = Project
browser.section.cases = Cases
browser.section.geometry = Geometry
browser.empty-project = (no project loaded)
browser.empty-cases = (no cases yet)
browser.run-selected = Run selected
browser.prepare-selected = Prepare selected case (no execute)
browser.adapter-status.ready = Ready
browser.adapter-status.missing = Missing
browser.adapter-status.outdated = Outdated
browser.adapter-status.broken = Broken
browser.adapter-status.disabled = Disabled
browser.adapter-status.unregistered = Unregistered
browser.run-history.success = Last run succeeded
browser.run-history.failure = Last run failed
browser.run-history.never = Not yet run

# ---------------------------------------------------------------
# Ribbon menus
# ---------------------------------------------------------------
menu.file = File
menu.file.open-project = Open project…
menu.file.import-stl = Import STL…
menu.file.load-mesh = Load canonical mesh…
menu.file.exit = Exit
menu.view = View
menu.view.shaded = Shaded
menu.view.wireframe = Wireframe
menu.view.frame = Frame geometry (F)
menu.run = Run
menu.run.selected = Run selected case (F5)
menu.run.prepare = Prepare selected case (no execute)
menu.run.prepare-tooltip = Write the solver input deck to a temp workdir without spawning the solver. Useful for inspecting generated files or working without the underlying tool installed.
menu.run.from-prepared = Run from prepared workdir
menu.run.from-prepared-tooltip = Run the solver against the workdir from the last Prepare, including any edits you've made to the generated dicts. Skips the prepare step so your edits survive.
menu.run.sweep = Sweep selected case (materialise only)
menu.run.sweep-tooltip = Read the case's [sweep] block and materialise N derived case.toml files in a temp workdir. Doesn't execute the runs — opens the workdir for inspection.
menu.run.first = Run first case
menu.run.cancel = Cancel
menu.settings = Settings
menu.settings.preferences = Preferences…
menu.settings.reprobe = Re-probe adapters
menu.help = Help
menu.help.palette = Command palette (Ctrl+P)
menu.help.about = About Valenx

# ---------------------------------------------------------------
# Status bar
# ---------------------------------------------------------------
status.idle = Idle
status.preparing = Preparing case…
status.running = Running solver…
status.collecting = Collecting results…
status.cancelled = Cancelled
status.failed = Run failed
status.done = Run complete

# ---------------------------------------------------------------
# Dialogs
# ---------------------------------------------------------------
dialog.about.title = About Valenx
dialog.about.heading = Valenx
dialog.about.version = Version { $version }
dialog.about.tagline = A native open-source desktop simulation suite.
dialog.about.licence = Dual-licensed under MIT OR Apache-2.0. See LICENSE-MIT and LICENSE-APACHE.
dialog.about.unifies-1 = Unifies CFD, FEA, EM, chemistry, MD, battery, and
dialog.about.unifies-2 = multibody physics behind one shell — no browser,
dialog.about.unifies-3 = no subscription, no vendor lock-in.
dialog.about.tip = Tip: press Ctrl+P to open the command palette.

dialog.first-run.title-window = Welcome to Valenx
dialog.first-run.title = First-launch setup
dialog.first-run.subtitle = Detected adapters on this machine
dialog.first-run.skip = Skip
dialog.first-run.done = Done
dialog.first-run.summary-ready = Detected { $count } of { $total } adapters ready to run on this machine.
dialog.first-run.col-adapter = Adapter
dialog.first-run.col-status = Status
dialog.first-run.col-hint = Hint

dialog.settings.title = Settings
dialog.settings.section.appearance = Appearance
dialog.settings.section.viewport = Viewport
dialog.settings.section.residuals = Residuals
dialog.settings.section.adapters = Adapters
dialog.settings.section.privacy = Privacy
dialog.settings.theme-label = Theme:
dialog.settings.theme.auto = Auto
dialog.settings.theme.dark = Dark
dialog.settings.theme.light = Light
dialog.settings.shading-label = Default shading:
dialog.settings.shading.shaded = Shaded
dialog.settings.shading.wireframe = Wireframe
dialog.settings.residual-scale-label = Y-axis scale:
dialog.settings.residual-scale.log10 = log₁₀
dialog.settings.residual-scale.linear = linear
dialog.settings.reprobe-on-close = Re-probe adapters when this window closes
dialog.settings.crash-report-opt-in = Upload crash reports to the maintainers
dialog.settings.crash-report-explainer = Reports stay on disk regardless. With this off, the next launch asks before sending. With this on, the next launch sends automatically. Reports are sanitised — usernames, paths, UUIDs, and SHA-256 hashes are redacted before write.
dialog.settings.crash-report-open-folder = Open crashes folder…
dialog.settings.crash-report-folder-tooltip = Reveals { $path } in your file browser. Inspect or delete reports manually before sending.
dialog.settings.reset-defaults = Reset to defaults
dialog.settings.persistence-note = Settings persist to <state_dir>/settings.json on close.

dialog.crash-report.title = A previous run crashed
dialog.crash-report.body = Valenx detected one or more unsent crash reports. Send them to the maintainers?
dialog.crash-report.send = Send
dialog.crash-report.discard = Discard
dialog.crash-report.privacy-link = What's in a crash report?

# ---------------------------------------------------------------
# Command palette
# ---------------------------------------------------------------
palette.placeholder = Type a command…
palette.command.run-selected = Run selected case
palette.command.prepare-selected = Prepare selected case
palette.command.cancel-run = Cancel current run
palette.command.open-project = Open project
palette.command.reprobe-adapters = Re-probe adapters
palette.command.toggle-shading = Toggle shading mode
palette.command.frame-mesh = Frame mesh in viewport
palette.command.show-first-run = Show first-launch wizard
palette.command.show-settings = Settings…

# ---------------------------------------------------------------
# Tooltips
# ---------------------------------------------------------------
tooltip.adapter-status = Adapter probe status — hover for details
tooltip.history-glyph = Run history for this case
tooltip.field-picker = Switch which scalar drives the field overlay
tooltip.time-slider = Scrub through transient snapshots

# ---------------------------------------------------------------
# Error messages — user-facing surface; technical details go to logs
# ---------------------------------------------------------------
error.tool-not-installed = Tool `{ $tool }` is not installed or not on PATH.
error.case-malformed = Case `{ $case }` is malformed: { $reason }.
error.run-failed = Run failed with exit code { $code }.
error.io = Could not access `{ $path }`: { $reason }.
