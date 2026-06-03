# Design sources

Everything under `docs/design/` is **canonical design source material**
— versioned with the code, editable in any SVG/Markdown tool.

Layout:

```
docs/design/
├── README.md              ← this file
├── icon-inventory.md      ← every icon, role by role
├── mockups/               ← SVG mockups of every screen
│   └── README.md
└── patterns/              ← pattern documentation (one per pattern)
    └── README.md
```

## How to contribute

- **Mockups** — SVG at the size the app renders (include a viewBox).
  One file per screen or component state. Commit the working file.
- **Icons** — SVG at 24×24. Named after the role, not the shape
  (`icon-run.svg` not `icon-play-triangle.svg`).
- **Pattern docs** — Markdown in `patterns/`, one per pattern, with
  problem / when to use / when not to / example / accessibility notes
  (see the design-ships process in [DESIGN.md](../../DESIGN.md)).

Reviews follow the normal PR flow; UI PRs require a mockup link or
reference in the PR body (see `.github/PULL_REQUEST_TEMPLATE.md`).

## Tools

- **Figma / Penpot / Inkscape** — whatever you're comfortable with.
  Export to SVG and commit. Mockups are not the source of truth in any
  proprietary tool — SVG-in-repo is.
- **Tokens** — color / spacing / type values come from
  `crates/valenx-design-tokens/tokens.json`. The tokens crate emits
  an SVG palette under `target/tokens-export/palette.svg` on build
  that you can drop into your design tool.
