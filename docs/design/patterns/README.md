# Patterns

Each pattern is one Markdown file documenting an opinionated
composition the whole app uses. Structure, per
[DESIGN.md § 5 Layer C](../../../DESIGN.md):

- **Problem** the pattern solves
- **When to use**
- **When not to use**
- **Worked example** (with a mockup link if visual)
- **Accessibility notes**

## Year-1 patterns (per DESIGN.md § 5)

- `ribbon-tab.md` — ribbon tab composition; how adapters contribute
- `browser-tree.md` — browser-tree conventions
- `log-viewer.md` — log viewer
- `dialog.md` — dialog composition rules (modal / sheet / popover)
- `long-running-op.md` — long-running operation feedback
- `validation-error.md` — validation + error surfacing
- `command-palette.md` — command palette invocation
- `units.md` — units display + input

Each file lands with the pattern's implementation PR, not before.
This README is a placeholder until the first pattern ships.
