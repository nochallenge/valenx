#!/usr/bin/env bash
#
# Build the Valenx installer for the host OS — dispatches by `uname -s`:
#
#   Linux   → .deb (via cargo-deb)  +  .rpm (via cargo-generate-rpm)
#   Darwin  → .dmg (via cargo-bundle → create-dmg)
#   *       → fall back; run scripts/build-installer.ps1 on Windows
#
# The output(s) land under either target/debian/, target/generate-rpm/,
# or dist/ (.dmg). The CI flows in .github/workflows/release.yml run
# the same commands.
#
# Prereqs:
#
#   Linux:
#     cargo install --locked cargo-deb cargo-generate-rpm
#     sudo apt-get install -y libgl1-mesa-dev libwayland-dev \
#         libxkbcommon-dev libasound2-dev pkg-config       # Debian/Ubuntu
#     sudo dnf install -y mesa-libGL-devel libxkbcommon-devel \
#         wayland-devel alsa-lib-devel                     # Fedora/RHEL
#
#   macOS:
#     cargo install --locked cargo-bundle
#     brew install create-dmg jq
#
#   Windows: use scripts/build-installer.ps1 (this script will tell you).

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/.." && pwd)"
cd "${repo_root}"

# ---------------------------------------------------------------------
# Linux: .deb + .rpm out of one release build.
# ---------------------------------------------------------------------
build_linux() {
    echo "[1/4] Building valenx release binaries..."
    cargo build --release -p valenx-app -p valenx-core -p valenx-fields \
        -p valenx-mesh -p valenx-audit -p valenx-export

    echo "[2/4] Packaging .deb via cargo-deb..."
    cargo deb -p valenx-app --no-build

    echo "[3/4] Packaging .rpm via cargo-generate-rpm..."
    cargo generate-rpm -p crates/valenx-app

    echo "[4/4] Installer artefacts:"
    for deb in "${repo_root}/target/debian/"*.deb; do
        [[ -f "${deb}" ]] || continue
        size_kb=$(( $(stat -c%s "${deb}" 2>/dev/null || stat -f%z "${deb}") / 1024 ))
        echo "  ${deb}  (${size_kb} KB)"
    done
    for rpm in "${repo_root}/target/generate-rpm/"*.rpm; do
        [[ -f "${rpm}" ]] || continue
        size_kb=$(( $(stat -c%s "${rpm}" 2>/dev/null || stat -f%z "${rpm}") / 1024 ))
        echo "  ${rpm}  (${size_kb} KB)"
    done
}

# ---------------------------------------------------------------------
# macOS: .app (cargo-bundle) wrapped in a .dmg with the drag-to-
# Applications UX (via Homebrew's `create-dmg`).
# ---------------------------------------------------------------------
build_macos() {
    echo "[1/3] Building .app bundle via cargo-bundle..."
    cargo bundle --release -p valenx-app

    version=$(cargo metadata --no-deps --format-version 1 \
        | jq -r '.packages[] | select(.name=="valenx-app") | .version')
    bundle_dir="target/release/bundle/osx"
    if [[ ! -d "${bundle_dir}" ]]; then
        echo "ERROR: cargo-bundle output missing at ${bundle_dir}" >&2
        exit 1
    fi
    mkdir -p dist
    dmg_path="dist/Valenx-${version}.dmg"
    rm -f "${dmg_path}"  # create-dmg refuses to overwrite

    echo "[2/3] Wrapping .app in .dmg via create-dmg (version ${version})..."
    create-dmg \
        --volname "Valenx" \
        --window-pos 200 120 \
        --window-size 600 400 \
        --icon-size 100 \
        --icon "Valenx.app" 150 200 \
        --hide-extension "Valenx.app" \
        --app-drop-link 450 200 \
        --hdiutil-quiet \
        "${dmg_path}" \
        "${bundle_dir}/"

    echo "[3/3] Installer artefact:"
    size_mb=$(( $(stat -f%z "${dmg_path}") / 1024 / 1024 ))
    echo "  ${dmg_path}  (${size_mb} MB)"
}

case "$(uname -s)" in
    Linux)
        build_linux
        ;;
    Darwin)
        build_macos
        ;;
    MINGW*|MSYS*|CYGWIN*)
        echo "Windows host detected — run scripts/build-installer.ps1 instead." >&2
        exit 1
        ;;
    *)
        echo "Unsupported host OS: $(uname -s)" >&2
        echo "Supported: Linux (.deb + .rpm), Darwin (.dmg), Windows (use .ps1)." >&2
        exit 1
        ;;
esac
