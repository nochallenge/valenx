//! Build script for `valenx-app`.
//!
//! On Windows targets, embed `wix/valenx.ico` into the `valenx.exe`
//! Windows resource section so the running process picks up the
//! Valenx app icon in:
//!
//! * the taskbar / Alt+Tab switcher,
//! * File Explorer when browsing to `valenx.exe`, and
//! * Start Menu shortcuts that the WiX installer creates (the .lnk
//!   files also carry an explicit `Icon` attribute pointing at the
//!   same .ico, but a stripped `valenx.exe` shipped standalone still
//!   gets the right Explorer icon thanks to this step).
//!
//! On every other target this is a no-op so the workspace continues
//! to build cleanly on Linux + macOS without pulling Windows-only
//! tooling into the build graph.
//!
//! Implementation uses the `embed-resource` crate which invokes
//! `windres` / `llvm-rc` / MSVC `rc.exe` (whichever is on the host)
//! against a tiny generated `.rc` script. The dependency itself is
//! gated to `cfg(windows)` in `Cargo.toml` under
//! `[target.'cfg(windows)'.build-dependencies]`, so non-Windows
//! builds never even pull it.

#[cfg(target_os = "windows")]
fn main() {
    use std::path::PathBuf;

    let manifest_dir =
        std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR is always set by cargo");
    let out_dir = std::env::var("OUT_DIR").expect("OUT_DIR is always set by cargo");

    let ico_path: PathBuf = [manifest_dir.as_str(), "wix", "valenx.ico"]
        .iter()
        .collect();
    let rc_path: PathBuf = [out_dir.as_str(), "valenx.rc"].iter().collect();

    // Resource ID 1 + ICON convention: explorer.exe picks the
    // numerically lowest icon resource as the file icon. The path
    // gets double-backslash-escaped for the .rc string literal.
    let rc_body = format!(
        "1 ICON \"{}\"\r\n",
        ico_path.display().to_string().replace('\\', "\\\\"),
    );
    std::fs::write(&rc_path, rc_body).expect("failed to write valenx.rc");

    // Re-run only when these sources change.
    println!("cargo:rerun-if-changed=wix/valenx.ico");
    println!("cargo:rerun-if-changed=build.rs");

    embed_resource::compile(&rc_path, embed_resource::NONE)
        .manifest_required()
        .expect("failed to embed Windows resources");
}

#[cfg(not(target_os = "windows"))]
fn main() {
    // Nothing to do on non-Windows targets.
}
