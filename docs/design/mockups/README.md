# Mockups

SVG source files for every Valenx screen and component state.

## Conventions

- **One SVG file per screen or component state.**
- **Name by role:** `home-view.svg`, `workspace-shell.svg`,
  `command-palette.svg`, `dialog-new-project.svg`.
- **Include a `viewBox`** matching the logical pixel dimensions the
  app renders at (commonly 1440×900 for default; 1920×1080 for
  large-screen checks; 1024×700 for minimum supported).
- **Reference tokens by role in comments** (`<!-- bg: surface-0 -->`)
  rather than embedding raw hex when possible. The palette at
  `target/tokens-export/palette.svg` has every role.
- **Dark theme is the default.** A matching `-light` variant ships
  for each screen that significantly differs in light mode.

## Status

This directory is empty at repository initialization. Mockups land
in Phase 0 per [DESIGN.md § 28 Timeline](../../../DESIGN.md).
Initial set:

- `home-view.svg`
- `workspace-shell.svg`
- `command-palette.svg`

Subsequent phases add mockups for screens in their phase inventory.
