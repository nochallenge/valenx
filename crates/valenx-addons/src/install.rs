//! Install / update / uninstall helpers.
//!
//! v1 ships **manual install by source directory** — the user picks a
//! directory that contains a `valenx-addon.toml` and the install
//! function recursively copies it to `target_dir/{manifest.name}/`.
//! Git-clone-based install is the natural follow-up once the GitHub-
//! search story lands (currently parked, see lib.rs).

use std::fs;
use std::io;
use std::path::Path;

use crate::error::AddonError;
use crate::manifest::read_manifest_at;
use crate::registry::LocalAddon;

/// Per-file size cap applied while copying an add-on source tree.
/// Round-6 hardening: pre-fix `fs::copy(&from, &to)` would happily
/// stream a 100 GiB asset out of the source dir, filling the user's
/// disk before any error surfaced. 64 MiB is far more than any
/// honest Python / WASM asset (an add-on with a 64 MiB single file
/// is already suspect); the cap surfaces a typed error so the user
/// can spot the misbehaving source tree.
///
/// Round-14 H1: kept for back-compat exports / docs only — the
/// actual cap enforced at copy time is
/// [`valenx_core::adapter_helpers::MAX_COPY_DIR_FILE_BYTES`] (8 GiB)
/// since install.rs delegates the recursive copy to the shared
/// helper for the round-6 hardening (symlink rejection, depth cap,
/// per-file cap) and to keep the security policy single-sourced.
pub const MAX_ADDON_FILE_BYTES: u64 = 64 * 1024 * 1024;

/// Maximum recursion depth while copying an add-on source tree.
/// 32 levels is far past any sensible source layout; the cap stops
/// a symlink loop or pathological deep nesting from blowing the
/// stack via the recursive copy function.
///
/// Round-14 H1: mirrors
/// [`valenx_core::adapter_helpers::MAX_COPY_DIR_DEPTH`] — the actual
/// depth check fires inside the shared helper.
pub const MAX_ADDON_COPY_DEPTH: usize = 32;

/// Install an add-on by copying `source_dir` to
/// `target_dir/{manifest.name}/`.
///
/// # Errors
///
/// - [`AddonError::SourceNotFound`] when `source_dir` does not exist.
/// - [`AddonError::Manifest`] / [`AddonError::InvalidManifest`] when
///   the manifest is malformed.
/// - [`AddonError::AlreadyInstalled`] if a same-named add-on already
///   exists under `target_dir`.
/// - [`AddonError::Io`] for any filesystem error.
pub fn install_from_dir(source_dir: &Path, target_dir: &Path) -> Result<LocalAddon, AddonError> {
    if !source_dir.is_dir() {
        return Err(AddonError::SourceNotFound(source_dir.display().to_string()));
    }
    let manifest_path = source_dir.join("valenx-addon.toml");
    let manifest = read_manifest_at(&manifest_path)?;
    let install_path = target_dir.join(&manifest.name);
    if install_path.exists() {
        return Err(AddonError::AlreadyInstalled(manifest.name.clone()));
    }
    fs::create_dir_all(target_dir)?;
    copy_dir_recursive(source_dir, &install_path)?;
    Ok(LocalAddon {
        manifest,
        path: install_path,
    })
}

/// Re-install on top of an existing add-on. Internally: uninstall +
/// install. Useful for "Update" actions in the UI.
///
/// # Errors
///
/// Same as [`install_from_dir`] plus [`AddonError::Io`] from the
/// implicit uninstall step.
pub fn update_from_dir(source_dir: &Path, target_dir: &Path) -> Result<LocalAddon, AddonError> {
    let manifest = read_manifest_at(&source_dir.join("valenx-addon.toml"))?;
    let install_path = target_dir.join(&manifest.name);
    if install_path.exists() {
        fs::remove_dir_all(&install_path)?;
    }
    install_from_dir(source_dir, target_dir)
}

/// Uninstall an add-on by removing its install directory.
///
/// Returns `Ok(false)` if the directory didn't exist.
///
/// # Errors
///
/// - [`AddonError::Io`] for any deletion error other than `NotFound`.
pub fn uninstall(install_path: &Path) -> Result<bool, AddonError> {
    match fs::remove_dir_all(install_path) {
        Ok(()) => Ok(true),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(false),
        Err(e) => Err(AddonError::Io(e)),
    }
}

/// Recursive directory copy for add-on installs.
///
/// Round-14 H1: delegates to
/// [`valenx_core::adapter_helpers::copy_dir_recursive`], the round-6
/// hardened helper the rest of the adapter ecosystem already uses
/// (precice adapter migration sister-fix from round-7). Pre-fix the
/// addons crate kept its own copy of the walker; that worked but
/// meant future hardening to the shared helper (e.g. tightening the
/// per-file cap or adding a per-tree size cap) would silently NOT
/// propagate to the addon install path. Single-sourcing the policy
/// avoids that drift class permanently.
///
/// The shared helper raises `AdapterError::Other`; we map every
/// helper error into [`AddonError::InvalidManifest`] (the closest
/// addon-side variant for "the source tree violates a copy policy")
/// so callers see typed errors with the same shape as the pre-fix
/// implementation.
fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<(), AddonError> {
    valenx_core::adapter_helpers::copy_dir_recursive(src, dst)
        .map_err(|e| AddonError::InvalidManifest(format!("{e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn tmpdir(name: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("valenx_addons_test_{name}"));
        let _ = fs::remove_dir_all(&p);
        p
    }

    fn make_addon_dir(dir: &Path, name: &str) {
        fs::create_dir_all(dir).unwrap();
        let manifest = format!(
            r#"
                name = "{name}"
                description = "Test addon"
                version = "0.1.0"

                [entry_point]
                kind = "python"
                module = "main"
            "#
        );
        fs::write(dir.join("valenx-addon.toml"), manifest).unwrap();
        fs::write(dir.join("main.py"), "def hello(): return 1").unwrap();
    }

    #[test]
    fn install_copies_source_to_target() {
        let src = tmpdir("install_src");
        let dst = tmpdir("install_dst");
        make_addon_dir(&src, "ok-addon");
        let local = install_from_dir(&src, &dst).unwrap();
        assert_eq!(local.manifest.name, "ok-addon");
        assert!(local.path.join("valenx-addon.toml").exists());
        assert!(local.path.join("main.py").exists());
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
    }

    #[test]
    fn install_into_already_existing_returns_already_installed() {
        let src = tmpdir("aei_src");
        let dst = tmpdir("aei_dst");
        make_addon_dir(&src, "dup");
        install_from_dir(&src, &dst).unwrap();
        let err = install_from_dir(&src, &dst).unwrap_err();
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
        assert!(matches!(err, AddonError::AlreadyInstalled(_)));
    }

    #[test]
    fn install_missing_source_returns_source_not_found() {
        let src = tmpdir("missing_src");
        let _ = fs::remove_dir_all(&src);
        let dst = tmpdir("missing_dst");
        let err = install_from_dir(&src, &dst).unwrap_err();
        assert!(matches!(err, AddonError::SourceNotFound(_)));
    }

    #[test]
    fn update_replaces_existing_install() {
        let src = tmpdir("upd_src");
        let dst = tmpdir("upd_dst");
        make_addon_dir(&src, "to-update");
        install_from_dir(&src, &dst).unwrap();
        // Modify source.
        fs::write(src.join("new_file.py"), "x = 2").unwrap();
        let local = update_from_dir(&src, &dst).unwrap();
        assert!(local.path.join("new_file.py").exists());
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
    }

    #[test]
    fn uninstall_removes_directory() {
        let src = tmpdir("uni_src");
        let dst = tmpdir("uni_dst");
        make_addon_dir(&src, "to-rm");
        let local = install_from_dir(&src, &dst).unwrap();
        assert!(local.path.exists());
        assert!(uninstall(&local.path).unwrap());
        assert!(!local.path.exists());
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
    }

    #[test]
    fn uninstall_missing_returns_false() {
        let p = std::env::temp_dir().join("valenx_addons_no_such_addon");
        let _ = fs::remove_dir_all(&p);
        assert!(!uninstall(&p).unwrap());
    }

    #[cfg(unix)]
    #[test]
    fn install_rejects_symlinked_source_directory() {
        // Round-6 RED→GREEN: a poisoned source dir that contains
        // a symlink to /etc could otherwise let the recursive copy
        // walk outside the user's chosen source tree. The fix
        // refuses any symlink the walker meets.
        let src = tmpdir("install_symlink_src");
        let dst = tmpdir("install_symlink_dst");
        make_addon_dir(&src, "symlink-addon");
        // Plant a symlink inside the source pointing at a system
        // directory. `std::os::unix::fs::symlink` is the standard
        // helper.
        let symlink_target = std::path::PathBuf::from("/etc");
        let symlink_at = src.join("evil_link");
        std::os::unix::fs::symlink(&symlink_target, &symlink_at).unwrap();

        let err = install_from_dir(&src, &dst).unwrap_err();
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
        match err {
            AddonError::InvalidManifest(msg) => {
                assert!(msg.contains("symlink"), "msg: {msg}");
            }
            other => panic!("expected InvalidManifest for symlink, got {other:?}"),
        }
    }

    #[test]
    fn install_rejects_oversized_source_file() {
        // Round-6 RED→GREEN: an addon source dir with a single
        // 100 MiB file would otherwise let the recursive copy
        // happily fill the user's disk. The size cap stops the
        // copy with a typed error.
        //
        // Round-14 H1: cap now lives on the shared helper
        // (`valenx_core::adapter_helpers::MAX_COPY_DIR_FILE_BYTES`)
        // at 8 GiB. We sparse-allocate one byte past that to
        // confirm the helper still rejects oversized files. The
        // sparse-file `set_len` call doesn't actually reserve disk;
        // ZFS/NTFS/ext4 all honour sparse semantics.
        let src = tmpdir("install_big_src");
        let dst = tmpdir("install_big_dst");
        make_addon_dir(&src, "big-addon");
        let big = src.join("payload.bin");
        let oversize_bytes =
            valenx_core::adapter_helpers::MAX_COPY_DIR_FILE_BYTES + 1;
        let f = std::fs::File::create(&big).unwrap();
        // set_len on a freshly-opened file allocates a sparse hole on
        // every modern filesystem; the 8 GiB doesn't hit the disk.
        f.set_len(oversize_bytes).unwrap();

        let err = install_from_dir(&src, &dst).unwrap_err();
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
        match err {
            AddonError::InvalidManifest(msg) => {
                assert!(msg.contains("exceeds"), "msg: {msg}");
                assert!(
                    msg.contains(
                        &valenx_core::adapter_helpers::MAX_COPY_DIR_FILE_BYTES.to_string()
                    ),
                    "msg: {msg}"
                );
            }
            other => panic!("expected InvalidManifest for oversized file, got {other:?}"),
        }
    }

    /// Round-14 H1 RED→GREEN: confirm that delegating to the shared
    /// helper preserves the symlink rejection on Unix. Pre-fix the
    /// addons crate had its own copy of the walker, so any future
    /// hardening to `valenx_core::adapter_helpers::copy_dir_recursive`
    /// would silently NOT propagate. This test re-runs the round-6
    /// symlink scenario through the delegated path to lock in that
    /// the policy is enforced end-to-end via the shared helper.
    #[cfg(unix)]
    #[test]
    fn install_rejects_symlinked_source_directory_via_shared_helper() {
        let src = tmpdir("install_symlink_shared_src");
        let dst = tmpdir("install_symlink_shared_dst");
        make_addon_dir(&src, "symlink-shared-addon");
        let symlink_target = std::path::PathBuf::from("/etc");
        let symlink_at = src.join("evil_link");
        std::os::unix::fs::symlink(&symlink_target, &symlink_at).unwrap();

        let err = install_from_dir(&src, &dst)
            .expect_err("shared helper must still reject symlinks");
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
        match err {
            AddonError::InvalidManifest(msg) => {
                assert!(msg.contains("symlink"), "msg: {msg}");
                // The shared helper namespaces its error message with
                // its own function name. Cross-checking the substring
                // confirms the delegated path is what's running.
                assert!(
                    msg.contains("copy_dir_recursive"),
                    "shared helper marker missing from msg: {msg}"
                );
            }
            other => panic!("expected InvalidManifest for shared-helper symlink, got {other:?}"),
        }
    }
}
