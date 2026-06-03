# Installers

Per-platform installer packaging scaffolds. Empty at pre-alpha; filled
in during Phase 1 as the native app becomes runnable.

## Planned layout

- `windows/` — MSI (via `cargo-wix`) and MSIX bundle
- `macos/` — `.app` bundle, DMG, Apple notarization scripts
- `linux/`
  - `deb/` — Debian/Ubuntu package
  - `rpm/` — RHEL/Fedora package
  - `appimage/` — portable AppImage
  - `flatpak/` — Flatpak manifest
- `tools/` — first-run tool-picker wizard bundle (downloaded separately
  to keep the base installer small per ROADMAP Section 10)

## Build flow

```powershell
# Per-platform builds driven by cargo-xtask (to be added in Phase 1):
cargo xtask installer windows
cargo xtask installer macos
cargo xtask installer linux --format deb
```

Binaries are signed per PLATFORM conventions (see
[SECURITY.md](../SECURITY.md) for signing + provenance).
