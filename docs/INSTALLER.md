# Valenx installer guide

End-to-end runbook for the three OS-native installers produced by the
release pipeline (`.github/workflows/release.yml`):

* **Windows** — `.msi` via [cargo-wix](https://volks73.github.io/cargo-wix/)
* **macOS** — `.dmg` containing a signed `.app` via [cargo-bundle](https://github.com/burtonageo/cargo-bundle) + [create-dmg](https://github.com/create-dmg/create-dmg)
* **Linux** — `.deb` via [cargo-deb](https://github.com/kornelski/cargo-deb) and `.rpm` via [cargo-generate-rpm](https://github.com/cat-in-136/cargo-generate-rpm)

Every installer ships the same branded icon (the white "V" on a
`#4B9EFF` rounded square — generated from a single source design that
feeds Windows `.ico`, macOS `.icns`, and Linux `.png`+`.svg`), and
every installer puts Valenx into the OS-native launcher so users can
pin it for one-click access.

---

## Quick start (end users)

| OS | Download | Install | Launch |
| --- | --- | --- | --- |
| Windows 10 / 11 | `Valenx-<ver>-x86_64.msi` from [Releases](https://github.com/nochallenge/valenx/releases) | Double-click → next → next → finish | Start Menu → type "Valenx" |
| macOS 11+ | `Valenx-<ver>.dmg` from Releases | Open → drag **Valenx.app** onto the **Applications** alias | Launchpad → Valenx (or Spotlight) |
| Ubuntu / Debian | `valenx_<ver>_amd64.deb` | `sudo apt install ./valenx_<ver>_amd64.deb` | Activities → "Valenx" |
| Fedora / RHEL | `valenx-<ver>.x86_64.rpm` | `sudo dnf install ./valenx-<ver>.x86_64.rpm` | App grid → "Valenx" |

---

## Pinning Valenx for quick launch

OS-specific recipes — same outcome: a one-click launcher entry that
lives outside the start screen / launcher grid.

### Windows

Windows enforces **you can only pin shortcuts, not raw `.exe` files**.
The installer creates two shortcuts so all three pinning targets work:

1. **Pin to taskbar** — Open Start Menu → type `Valenx` → click to
   launch → right-click the running taskbar icon → **Pin to taskbar**.
   (Windows 11 puts this on an inner submenu.) Alternative: right-click
   the Start Menu **Valenx** entry → **More → Pin to taskbar**.
2. **Pin to Start** — Right-click the Start Menu **Valenx** entry →
   **Pin to Start**.
3. **Desktop shortcut** — The installer creates `Desktop\Valenx.lnk`
   by default. To opt out, uncheck *"Create a desktop shortcut"* on the
   install dialog, or install silently with
   `msiexec /i Valenx-<ver>-x86_64.msi INSTALLDESKTOPSHORTCUT=0 /quiet`.

### macOS

1. Open the `.dmg` and drag **Valenx.app** onto the **Applications**
   alias (standard macOS install flow).
2. Launch Valenx once (Spotlight → "Valenx", or double-click in
   `/Applications`). The dock icon appears.
3. **Right-click the dock icon → Options → Keep in Dock.** That's
   macOS's permanent-pin equivalent — the icon stays in the dock even
   when Valenx isn't running.
4. To put it on the desktop, ⌥-drag from `/Applications/Valenx.app`
   to the desktop (creates a `.app` alias).

### Linux (GNOME, KDE Plasma, XFCE, Cinnamon, Unity, …)

Because the `.deb` / `.rpm` ships a freedesktop `.desktop` file and
hicolor icons, Valenx appears in every standards-compliant launcher
automatically (no manual `desktop-file-install` needed).

* **GNOME** — Press *Activities* → search "Valenx" → right-click →
  **Add to Favorites** (pins to the left dock).
* **KDE Plasma** — Open the application menu → search "Valenx" →
  right-click → **Add to Favorites** (kickoff menu) or **Add to Task
  Manager** (panel).
* **XFCE / MATE / Cinnamon** — Right-click the panel → **Add to panel**
  → **Application Launcher** → pick Valenx.
* **Unity / Pantheon (elementary)** — Drag the Valenx icon from the
  applications launcher onto the dock.
* **Desktop file copy** (works everywhere) —
  `cp /usr/share/applications/valenx.desktop ~/Desktop/ && chmod +x ~/Desktop/valenx.desktop`
  drops a clickable icon directly on the desktop.

---

## Building installers locally

The repo ships `scripts/build-installer.{ps1,sh}` as one-command
wrappers per host OS. Pick the section that matches the machine you're
building on.

### Windows (.msi)

Prereqs:

```powershell
rustup toolchain install stable
cargo install cargo-wix --locked
# Then install WiX Toolset 3.x from https://wixtoolset.org/docs/wix3/
# (provides candle.exe + light.exe; cargo-wix auto-discovers it).
```

Build:

```powershell
pwsh scripts/build-installer.ps1
# Or manually:
cargo build --release -p valenx-app
cargo wix -p valenx-app --no-build --nocapture
```

Output: `target/wix/Valenx-<version>-x86_64.msi`.

### macOS (.dmg)

Prereqs:

```sh
rustup toolchain install stable
cargo install cargo-bundle --locked
brew install create-dmg jq
```

Build:

```sh
bash scripts/build-installer.sh
# Or manually:
cargo bundle --release -p valenx-app
VERSION=$(cargo metadata --no-deps --format-version 1 \
  | jq -r '.packages[] | select(.name=="valenx-app") | .version')
mkdir -p dist
create-dmg \
  --volname "Valenx" \
  --window-size 600 400 \
  --icon "Valenx.app" 150 200 \
  --hide-extension "Valenx.app" \
  --app-drop-link 450 200 \
  --hdiutil-quiet \
  "dist/Valenx-${VERSION}.dmg" \
  "target/release/bundle/osx/"
```

Output: `dist/Valenx-<version>.dmg` containing `Valenx.app` and an
alias to `/Applications`.

### Linux (.deb + .rpm)

Prereqs (Debian/Ubuntu):

```sh
sudo apt-get install -y libgl1-mesa-dev libwayland-dev \
    libxkbcommon-dev libasound2-dev pkg-config
cargo install --locked cargo-deb cargo-generate-rpm
```

Prereqs (Fedora/RHEL/openSUSE):

```sh
sudo dnf install -y mesa-libGL-devel libxkbcommon-devel \
    wayland-devel alsa-lib-devel
cargo install --locked cargo-deb cargo-generate-rpm
```

Build (one command, both formats):

```sh
bash scripts/build-installer.sh
# Or manually:
cargo build --release -p valenx-app -p valenx-core -p valenx-fields \
    -p valenx-mesh -p valenx-audit -p valenx-export
cargo deb -p valenx-app --no-build
cargo generate-rpm -p crates/valenx-app
```

Output:

* `target/debian/valenx_<version>_amd64.deb`
* `target/generate-rpm/valenx-<version>-1.x86_64.rpm`

---

## Testing the installers

### Windows

```cmd
msiexec /i target\wix\Valenx-<ver>-x86_64.msi /l*v install.log
```

Then verify:

1. Open Start, type `Valenx` — Valenx entry with the blue "V" icon.
2. `Valenx.lnk` is on the desktop (unless you unchecked the box).
3. `appwiz.cpl` shows Valenx with the blue icon, the repository URL as
   the "About" link, and only an **Uninstall** button.
4. Launch Valenx → right-click taskbar icon → **Pin to taskbar** is
   enabled.

Uninstall: `appwiz.cpl → Valenx → Uninstall`, or
`msiexec /x {7530132B-431C-4F23-B9B3-2744F5BE49E1} /qb`.

### macOS

```sh
hdiutil attach dist/Valenx-<ver>.dmg
ls /Volumes/Valenx/   # should list Valenx.app + Applications alias
cp -R /Volumes/Valenx/Valenx.app /Applications/
hdiutil detach /Volumes/Valenx/
open -a Valenx
```

Verify:

1. Launchpad shows the Valenx icon (the brand "V").
2. Spotlight ⌘-Space → "Valenx" finds it.
3. Right-click the dock icon → **Options → Keep in Dock** is offered.

Uninstall: drag `/Applications/Valenx.app` to the Trash, then empty.

### Linux

```sh
# Debian/Ubuntu
sudo apt install ./target/debian/valenx_*_amd64.deb
which valenx                                     # /usr/bin/valenx
test -f /usr/share/applications/valenx.desktop && echo "desktop OK"
test -f /usr/share/icons/hicolor/256x256/apps/valenx.png && echo "icon OK"
update-desktop-database -q ~/.local/share/applications 2>/dev/null
desktop-file-validate /usr/share/applications/valenx.desktop  # silent = pass

# Fedora/RHEL
sudo dnf install ./target/generate-rpm/valenx-*.x86_64.rpm
# Same verification steps as above
```

Launch from the application menu, confirm the icon renders, then
pin per the recipes in the *Pinning* section above.

Uninstall:

```sh
sudo apt remove valenx    # Debian/Ubuntu
sudo dnf remove valenx    # Fedora/RHEL
```

---

## Internals & cross-references

| File | Role |
| --- | --- |
| `crates/valenx-app/wix/main.wxs` | Hand-authored WiX 3.x template (Windows). Defines install layout, shortcuts, ARP entries, MajorUpgrade. |
| `crates/valenx-app/wix/valenx.ico` | Multi-resolution Windows app icon (16/32/48/256). |
| `crates/valenx-app/build.rs` | On Windows, embeds `wix/valenx.ico` into `valenx.exe` via `embed-resource`. No-op on Linux + macOS. |
| `crates/valenx-app/Cargo.toml` `[package.metadata.wix]` | Tells cargo-wix to skip the EULA dialog and pins the `upgrade-guid`. |
| `crates/valenx-app/Cargo.toml` `[package.metadata.bundle]` | cargo-bundle config. `icon = ["../../packaging/macos/valenx.icns"]` wires the macOS icon. |
| `crates/valenx-app/Cargo.toml` `[package.metadata.deb]` | cargo-deb config. `assets` list installs `valenx.desktop` + hicolor icons. |
| `crates/valenx-app/Cargo.toml` `[package.metadata.generate-rpm.assets]` | cargo-generate-rpm config. Same layout as `.deb`. |
| `packaging/linux/valenx.desktop` | freedesktop XDG launcher entry. `Categories=Science;Education;Engineering;…`. |
| `packaging/linux/valenx-256.png` | 256×256 PNG for older launchers. |
| `packaging/linux/valenx.svg` | Scalable vector icon for HiDPI launchers. |
| `packaging/macos/valenx.icns` | Apple Icon Image format. Multi-resolution (16 → 1024). |
| `scripts/build_icon.ps1` | Regenerates `valenx.ico` from the brand colour (Windows-only, uses .NET `System.Drawing`). |
| `scripts/build_icns.py` | Regenerates `valenx.icns` from the source PNG (Python + Pillow). |
| `scripts/build-installer.ps1` | One-command local MSI build (Windows). |
| `scripts/build-installer.sh` | One-command local `.deb` / `.rpm` / `.dmg` build — dispatches by `uname -s`. |
| `.github/workflows/release.yml` | CI: manual `workflow_dispatch` only — builds + signs + uploads all four artifacts when triggered via `gh workflow run release.yml`. Tag-push triggers were removed in round-8 to prevent accidental re-builds of historical tags. |

### GUIDs (must stay stable across Windows releases)

| Element | GUID |
| --- | --- |
| UpgradeCode (product identity) | `7530132B-431C-4F23-B9B3-2744F5BE49E1` |
| Component `valenx.exe` | `63D5237D-A0CA-4114-9B6B-B0AC2A6D6CC2` |
| Component Start Menu shortcut | `AC45D57D-2565-41E7-A3C6-6C0DE877D3D8` |
| Component Desktop shortcut | `55163810-C628-4BF4-8ACC-ADB6B02FA8AA` |

Re-generating any of these would orphan existing installs (the old
version wouldn't upgrade cleanly), so they're permanent.

### Regenerating icons

When the brand mark changes, regenerate all three formats from the
same design:

```powershell
# Windows .ico
pwsh scripts/build_icon.ps1
```

```sh
# macOS .icns (requires Python + Pillow)
pip install Pillow
python scripts/build_icns.py
```

The Linux `.png` + `.svg` are committed directly to
`packaging/linux/` — re-render with any vector tool and commit.

### Why a custom WiX template instead of `cargo wix init`?

`cargo wix init`'s default template installs `valenx.exe` into
`Program Files\Valenx\` and stops there — no Start Menu entry, no
shortcut, no app icon in Add/Remove Programs. Without a Start Menu
shortcut, **Windows refuses to let users pin the app**. Our custom
template adds:

* `ProgramMenuFolder → Valenx → Valenx.lnk` (mandatory — the pinning source)
* `DesktopFolder → Valenx.lnk` gated by `INSTALLDESKTOPSHORTCUT`
* `ARPPRODUCTICON` so the Add/Remove Programs entry isn't a generic
  installer icon
* `MajorUpgrade` so re-runs upgrade in place instead of installing
  side-by-side
* `WixUI_InstallDir` flow so users can pick the install location
  (still defaults to `Program Files`)

---

## Code signing

Signing is **gated on per-OS secrets** in CI. Without secrets,
installers ship unsigned and end-users see OS-specific "untrusted
publisher" warnings on first launch which they can bypass.

| OS | Secrets | What gets signed | First-launch UX without signing |
| --- | --- | --- | --- |
| Windows | `WINDOWS_CERT` (base64 PFX) + `WINDOWS_CERT_PASSWORD` | The `.msi` via Authenticode (`signtool sign /fd sha256 /tr http://timestamp.digicert.com /td sha256`) | SmartScreen "Windows protected your PC" → **More info → Run anyway** |
| macOS | `APPLE_ID` + `APPLE_TEAM_ID` + `APPLE_APP_PASSWORD` + `APPLE_DEVELOPER_ID_CERT` (base64 P12) + `APPLE_DEVELOPER_ID_CERT_PASSWORD` | The `.app` (codesign --options runtime --timestamp), then the `.dmg`, then notarytool-notarise + staple | Gatekeeper "Valenx can't be opened because Apple cannot check it" → System Preferences → Security → **Open Anyway** |
| Linux | None | n/a — distros sign their package indices, not individual packages | apt / dnf "untrusted source" warning if not installed from a configured repo |

To Authenticode-sign locally on Windows, set `WINDOWS_CERT` +
`WINDOWS_CERT_PASSWORD` and re-run the signtool block from
`release.yml` manually after `cargo wix`. Apple signing locally
requires a Developer ID Application certificate in Keychain.
