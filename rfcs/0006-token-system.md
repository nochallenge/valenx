# RFC 0006: Design Token System and Pipeline

- **Status:** Accepted (initial schema; additive changes expected)
- **Author(s):** BDFL
- **Created:** 2026-04-23
- **Discussion PR:** (this commit)
- **Tracking issue:** TBD

---

## Summary

Define the design-token schema, the JSON source of truth, and the
build pipeline that produces Rust `const` values for the app and
export artifacts for mockup tools. Tokens are the single source of
truth for every design primitive — color, typography, spacing,
radius, motion, shadow, borders, z-index.

---

## Motivation

Without tokens, every UI component hard-codes values, themes
become impossible to swap, and rebranding later is a rewrite.
With tokens, a theme switch is a mapping swap and the UI never
knows the difference.

Token-based design is industry-standard (Material, Polaris,
Fluent, Apple HIG all do it). What this RFC commits to is
Valenx-specific:

- **One JSON source of truth** in the repo — not Figma-locked, not
  multi-source
- **Generated Rust consts** for compile-time use in `egui` code —
  no runtime JSON parsing, no string keys
- **Generated SVG / palette exports** for mockup tools, so
  mockups can match the actual app
- **Strict naming by role**, never by value — so a theme swap is
  a mapping swap
- **Build is deterministic** — same JSON → same Rust bytes, every
  time, for reproducible binaries and snapshot tests

---

## Guide-level explanation

### The JSON source

Lives at `crates/valenx-design-tokens/tokens.json`. Edited like
any other source file; reviewed in PRs; changes that alter
token semantics (role renames, theme splits) require an RFC.

Structure:

```json
{
  "$schema": "./tokens.schema.json",
  "version": "0.1.0",

  "color": {
    "surface": { "0": "#0E1116", "1": "#161A22", "2": "#1D232C",
                  "3": "#252C37", "4": "#2E3642", "5": "#38414F" },
    "text":    { "1": "#E6EBF1", "2": "#A8B2BF", "3": "#6B7482" },
    "accent":  { "primary": "#4B9EFF",
                  "success": "#3FC27A",
                  "warning": "#E0A43D",
                  "error":   "#E55B5B",
                  "info":    "#4B9EFF" },
    "physics": { "cfd": "#4B9EFF", "fea": "#E0A43D",
                  "em":  "#9D7DF0", "chem": "#3FC27A",
                  "md":  "#E55B5B", "battery": "#F07DA6" }
  },

  "type": {
    "family": {
      "ui":   "Inter",
      "mono": "JetBrains Mono"
    },
    "size": {
      "xs": 11, "sm": 12, "base": 14,
      "lg": 16, "xl": 18, "2xl": 24, "3xl": 32
    },
    "weight": {
      "regular":  400,
      "medium":   500,
      "semibold": 600
    },
    "lineHeight": {
      "tight":  1.2,
      "normal": 1.45,
      "loose":  1.6
    }
  },

  "space":  { "0": 0, "1": 4, "2": 8, "3": 12, "4": 16, "5": 20,
               "6": 24, "7": 32, "8": 40, "9": 48, "10": 64,
               "11": 80, "12": 96 },

  "radius": { "none": 0, "sm": 2, "md": 6, "lg": 10, "full": 9999 },

  "shadow": { "0": "none",
               "1": "0 1px 2px rgba(0,0,0,0.15)",
               "2": "0 4px 8px rgba(0,0,0,0.20)",
               "3": "0 8px 24px rgba(0,0,0,0.30)" },

  "motion": {
    "duration": { "fast": 120, "base": 180, "slow": 300,
                   "camera": 400 },
    "easing":   { "standard":    "cubic-bezier(0.4, 0, 0.2, 1)",
                   "decelerate":  "cubic-bezier(0.0, 0, 0.2, 1)",
                   "emphasized":  "cubic-bezier(0.2, 0, 0, 1)" }
  },

  "border":  { "hairline": 0.5, "default": 1, "strong": 2 },

  "z": { "base": 0, "docked": 10, "overlay": 100,
          "modal": 1000, "tooltip": 2000 },

  "themes": {
    "dark":  { "extends": null, "overrides": {} },
    "light": {
      "extends": "dark",
      "overrides": {
        "color.surface.0": "#F5F6F8",
        "color.surface.1": "#FFFFFF",
        "color.surface.2": "#ECEEF2",
        "color.text.1":    "#0E1116",
        "color.text.2":    "#4A525D"
      }
    },
    "high-contrast": {
      "extends": "dark",
      "overrides": {
        "color.surface.0": "#000000",
        "color.surface.1": "#000000",
        "color.text.1":    "#FFFFFF",
        "color.accent.primary": "#FFFF00"
      }
    }
  }
}
```

### The pipeline

A `build.rs` in `valenx-design-tokens` runs at compile time:

1. Reads `tokens.json`
2. Validates against `tokens.schema.json` (a JSON Schema file also in the crate)
3. Generates `src/generated.rs` with Rust `const` values
4. Generates `target/tokens-export/` with palette SVG, Figma-importable JSON, and documentation reference

Consumers in other crates `use valenx_design_tokens::color::surface::S1;` — no runtime lookup, no string keys.

---

## Reference-level explanation

### Crate layout

```
crates/valenx-design-tokens/
├── Cargo.toml
├── build.rs                        # generates Rust from tokens.json
├── tokens.json                     # THE SOURCE OF TRUTH
├── tokens.schema.json              # JSON Schema for validation
├── src/
│   ├── lib.rs                      # re-exports, theme-switching API
│   ├── generated.rs                # generated by build.rs (git-ignored)
│   └── theme.rs                    # runtime theme application
└── tests/
    ├── schema.rs                   # tokens.json validates against schema
    ├── stability.rs                # snapshot of the generated API; prevents accidental renames
    └── theme_fallback.rs           # ensures every theme provides every role
```

### Generated Rust API

For color:

```rust
pub mod color {
    pub mod surface {
        pub const S0: Color = Color::from_hex(0x0E1116);
        pub const S1: Color = Color::from_hex(0x161A22);
        // ...
    }
    pub mod text {
        pub const T1: Color = Color::from_hex(0xE6EBF1);
        // ...
    }
    // ... accent, physics
}
```

For spacing, type, etc., similar module-scoped consts.

### Runtime theming

Even though the default theme compiles as `const`, theme swap is
runtime:

```rust
pub struct Theme { /* ... */ }

impl Theme {
    pub const DARK:   Theme = /* generated from tokens.json */;
    pub const LIGHT:  Theme = /* ... */;
    pub const HIGH_CONTRAST: Theme = /* ... */;
}

pub fn set_active_theme(theme: &'static Theme);
pub fn active() -> &'static Theme;
```

Components call `active().color.surface.s1` at paint time; theme
swap is atomic (a single pointer swap in a `RwLock`).

### Build-script behavior

- Deterministic output (sorted map keys, stable formatting)
- Errors on unknown token categories, malformed colors, missing
  theme overrides
- Emits `cargo:rerun-if-changed=tokens.json` and
  `cargo:rerun-if-changed=tokens.schema.json`
- If `tokens.json` is syntactically invalid, the build fails with
  a precise line + column

### Exports to mockup tools

Also written during `build.rs`:

```
target/tokens-export/
├── palette.svg               # visual palette, one swatch per role
├── tokens.figma.json         # importable by the Figma Tokens plugin (and Penpot)
├── tokens.css                # CSS custom properties, for the website
└── TOKENS.md                 # generated Markdown reference for docs
```

These make the mockup side of work (which often happens in Figma
or Penpot) track the canonical source automatically. The mockup
tool imports the Figma JSON; the mockups use live tokens; no
drift.

### Versioning the tokens themselves

The `"version"` field in `tokens.json` follows SemVer independently
from the crate version:

- **MAJOR** bump: a role was removed or renamed; consumers must
  update
- **MINOR** bump: a role was added; old consumers keep working
- **PATCH** bump: values changed but roles didn't (e.g., accent
  color shifted a few steps)

An additional stability-snapshot test (`tests/stability.rs`)
captures the generated API surface and fails if it changes
without a corresponding `tokens.json` version bump. Prevents
stealth renames.

### Themes

Themes are declared under `"themes"` in `tokens.json`. Each has
an `extends` (usually `"dark"`) and an `overrides` map. Every
generated theme is verified to cover every role — if a theme
omits a role and doesn't inherit one, `build.rs` fails.

Adding a theme is adding an entry; no code changes elsewhere.

### Naming discipline

**Tokens are named by role, never by value.** `surface.s1` —
not `gray-800`. `accent.primary` — not `blue`. This is enforced
by review, and the role names are aligned with the component
and pattern docs in `docs/design/patterns/`.

### Validation

- **Schema validation** via `jsonschema` in `build.rs`
- **Hex color parsing** at build time (not runtime) — invalid
  colors fail the build
- **Motion durations** constrained to 0–2000 ms
- **Space values** constrained to multiples of 4 (soft
  convention — warns if violated)

---

## Drawbacks

- **One more build step.** `build.rs` adds a few hundred
  milliseconds to cold compiles. Cached between changes; fine.
- **A JSON file is slightly less ergonomic to edit than a TOML
  file.** We chose JSON because it's what every design-tool
  import/export speaks natively. TOML would require a
  bidirectional converter.
- **Theme `overrides` using dotted paths** (`"color.surface.0"`)
  is string-typed and can silently typo. Mitigation: the build
  script validates every override path against the base theme.
- **Generated Rust code under `src/generated.rs`** is normal but
  can confuse readers who expect everything in `src/` to be
  authored. A comment at the top of the generated file clarifies.

---

## Rationale and alternatives

**TOML as the source.**
Considered. Rejected because every mockup tool imports / exports
JSON natively, and we want the mockup pipeline to be
zero-friction.

**YAML.**
Rejected. YAML's whitespace-sensitivity causes subtle bugs; no
native import in any design tool.

**Rust as the source** (put the tokens directly in `lib.rs` as
`const`).
Rejected. Mockup tools can't read Rust. Tokens-in-code also
makes theme definitions awkward (nested consts don't express
inheritance cleanly).

**Style Dictionary** (Amazon's toolchain).
Considered. Powerful but brings a Node.js dependency and more
ceremony than we need. Our `build.rs` does ~10% of Style
Dictionary's job and is ~100 lines.

**Runtime JSON loading on app start** instead of build-time
codegen.
Rejected. Adds a startup cost; loses type-safety; makes snapshot
tests flakier.

---

## Prior art

- **Material Design tokens** — the canonical reference for
  role-based naming. We borrow the pattern heavily.
- **GitHub Primer** — industrial-strength token system, good
  reference for the dark / light theme split.
- **Radix Colors** — scale-based palette (0–12 per hue) that's
  worth aligning with for a future expansion; our current 0–5
  surface scale is a subset.
- **Tailwind's design tokens** — influential format; our spacing
  scale is Tailwind-flavoured (4-pixel base).
- **Figma Variables** — the modern way to consume tokens in
  Figma. Our `tokens.figma.json` export targets this format.

---

## Unresolved questions

- **Accent-color user preference.** Do we let users pick their
  own accent like Slack / macOS? Low-cost to add; not in v1.
- **Density multiplier** — one global space-scale multiplier for
  "compact / comfortable / spacious" like Gmail. Straightforward
  extension; parked for post-1.0.
- **Dark / light auto-switching.** Follow OS preference? Time of
  day? Manual? Current plan: manual in Settings, with an
  explicit "follow OS" option; no automatic time switching.
- **Per-physics accent usage rules.** When does a component use
  `accent.primary` versus `color.physics.cfd`? Needs a style
  guide section in the component docs, not an RFC.

---

## Future possibilities

- **Token lint** — a custom clippy lint that rejects hard-coded
  colors / spacing / font sizes in UI code. Turns the "use
  tokens" principle into a compile-time check.
- **Preview server** — `cargo xtask tokens-preview` opens a
  native window showing every token in every theme.
- **Export to other tools as they appear** — Penpot has its own
  JSON format variant; we can add an adapter in `build.rs`
  without touching the source.
- **Semantic layer on top of roles** — e.g., `danger-action-
  background`, `neutral-surface-elevated`. Adds ergonomics but
  also complexity; evaluate after Year 1.
- **Brand-partner tokens** — a subfolder of tokens specifically
  for co-branded or institutional deployments (university or
  company color overrides). Supports theming of the app for
  classroom or enterprise contexts.
