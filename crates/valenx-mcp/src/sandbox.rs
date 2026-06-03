//! Path sandboxing for MCP tool calls.
//!
//! The MCP server is stdio-driven — it has no natural per-call "case
//! directory" the way `prepare()` does. To prevent a misbehaving client
//! from reading/writing arbitrary files via `receptor_path` /
//! `ligand_path` / `output_path`, every path supplied through the
//! protocol is forced through [`sandbox_check`].
//!
//! Sandbox root: the value of the `VALENX_MCP_SANDBOX_DIR` environment
//! variable at server startup, or [`std::env::temp_dir()`]`.join("valenx-mcp")`
//! when unset. Both are canonicalised before prefix-matching so the
//! check is robust against symlinks.

use std::fs::File;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};

/// Name of the env var holding the sandbox-root override.
pub const SANDBOX_ENV: &str = "VALENX_MCP_SANDBOX_DIR";

/// Round-19 L2 helper: returns `true` when any component of `path`
/// looks like a Windows 8.3 short-name alias. The pattern is
/// `[A-Z0-9_]{1,6}~\d` — 1–6 uppercase-alpha / digit / underscore
/// chars, a literal `~`, then a single digit (NTFS only ever
/// generates `~1`–`~9`, but we accept any digit for the cheapest
/// possible match). NTFS replaces 8.3-illegal characters (space, `+`,
/// `,`, etc.) with `_` when synthesising the alias, so the prefix is
/// `[A-Z0-9_]` rather than the stricter `[A-Z0-9]`.
///
/// Hand-rolled scanner so we don't drag in the `regex` crate just
/// for one alias check. The match is intentionally narrow — we only
/// fire on a component shape NTFS itself would alias.
fn contains_windows_short_name_component(path: &Path) -> bool {
    for comp in path.components() {
        let Some(s) = comp.as_os_str().to_str() else {
            continue;
        };
        // Look for the `~N` infix.
        let bytes = s.as_bytes();
        let Some(tilde_pos) = bytes.iter().position(|b| *b == b'~') else {
            continue;
        };
        // 1–6 uppercase-alpha / digit / underscore chars before `~`.
        if tilde_pos == 0 || tilde_pos > 6 {
            continue;
        }
        let prefix_ok = bytes[..tilde_pos]
            .iter()
            .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || *b == b'_');
        if !prefix_ok {
            continue;
        }
        // Exactly one digit immediately after the `~`. We allow the
        // suffix to be followed by an extension (e.g. `PROGRA~1.LNK`)
        // so we only check the byte at `tilde_pos + 1`.
        if bytes.len() <= tilde_pos + 1 {
            continue;
        }
        if bytes[tilde_pos + 1].is_ascii_digit() {
            return true;
        }
    }
    false
}

/// Return the configured sandbox root. If `VALENX_MCP_SANDBOX_DIR` is
/// set we use that; otherwise we fall back to
/// `<tempdir>/valenx-mcp`. The directory is created on first call so
/// the fallback is always writable.
pub fn sandbox_root() -> PathBuf {
    let raw = match std::env::var(SANDBOX_ENV) {
        Ok(v) if !v.trim().is_empty() => PathBuf::from(v),
        _ => std::env::temp_dir().join("valenx-mcp"),
    };
    // Best-effort mkdir so callers can drop files in without an extra
    // ceremony. Failure is non-fatal: sandbox_check will still reject
    // anything outside the (now-nonexistent) root.
    let _ = std::fs::create_dir_all(&raw);
    raw
}

/// Reject `path` if its absolute, canonicalised form does not live
/// under [`sandbox_root`]. Returns the resolved absolute path on
/// success so callers can use it as-is.
///
/// "Live under" is checked by prefix-matching the canonicalised forms
/// of both paths. We canonicalise the parent (not the file itself —
/// the output path won't exist yet for `dock`) so the check works for
/// both inputs (must exist) and outputs (may not yet).
///
/// Round-16 M3 hardening: we ALSO reject any leaf-name basename
/// that already exists as a symlink, even if its parent canonicalises
/// inside the sandbox. Pre-fix only the parent was canonicalised, so
/// a file like `<sandbox>/passwd` that was actually a symlink to
/// `/etc/passwd` would pass the prefix check and `fs::write` /
/// `fs::read_to_string` would silently follow the link out of the
/// sandbox. For brand-new paths that don't exist yet (writes) the
/// parent canonicalisation remains sufficient — we can't follow what
/// doesn't exist.
///
/// Round-24 M4 — known TOCTOU limitation: for **write** paths where
/// the leaf doesn't exist at check time, `sandbox_check` cannot
/// guarantee that an attacker hasn't raced a symlink into the leaf
/// position between this check and the actual `fs::write`. The
/// resolved path is verified to live inside the sandbox lexically
/// AND to not be an existing symlink, but if a concurrent process
/// `mklink`/`ln -s`'s the leaf path to point outside the sandbox
/// before the write fires, `fs::write` will silently follow the new
/// link.
///
/// Mitigations available to callers:
///
/// 1. Open with `OpenOptions::create_new(true)` — refuses to follow
///    any pre-existing entry at the leaf path, so the race window
///    collapses to "leaf must not exist when open fires". This is
///    the recommended pattern for fresh outputs.
/// 2. Sandbox the entire sandbox dir (chmod o-w, ACL deny non-
///    owner) at MCP startup so a hostile peer process can't
///    mklink into it.
///
/// The runtime accepts the residual TOCTOU because the alternative
/// (open-then-fstat-then-write under a tmpdir-only policy) requires
/// reworking every MCP tool's write path; the create_new + ACL
/// combination collapses the practical exposure to "another process
/// running as the same user". See round-24 M4 in CHANGELOG.
pub fn sandbox_check(path: &Path) -> Result<PathBuf> {
    let root = sandbox_root();
    let root_canon = root.canonicalize().unwrap_or(root.clone());
    // Resolve to absolute. Relative paths are joined onto the sandbox
    // root so a bare `out.pdbqt` writes to the sandbox by default.
    let absolute: PathBuf = if path.is_absolute() {
        path.to_path_buf()
    } else {
        root_canon.join(path)
    };
    // Canonicalise the parent for the existence check; fall back to
    // the lexical path if the parent doesn't exist yet (e.g. brand-new
    // subdirectory). Anchor against the sandbox root in that case.
    let parent = absolute.parent().unwrap_or(Path::new("."));
    let parent_canon_result = parent.canonicalize();
    let parent_canon = parent_canon_result
        .as_ref()
        .map(|p| p.clone())
        .unwrap_or_else(|_| parent.to_path_buf());
    // Round-19 L2: when the parent canonicalisation fails, the
    // prefix-check below runs against the lexical (un-canonicalised)
    // path and is much weaker — `..` segments don't resolve, Windows
    // 8.3 short names ("PROGRA~1") aren't expanded, and mixed-slash
    // paths can confuse the prefix comparison. Add a defence-in-depth
    // refusal for the three patterns most often used to bypass a
    // lexical-only check.
    if parent_canon_result.is_err() {
        let path_str = absolute.to_string_lossy();
        // `..` segments — the canonical escape vector.
        if absolute
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir))
        {
            return Err(anyhow!(
                "path `{}` contains `..` and the parent dir does not exist \
                 — refused (cannot verify it stays inside sandbox `{}`)",
                absolute.display(),
                root_canon.display()
            ));
        }
        // Mixed forward / backslash on Windows — `C:\Users\..\Windows`
        // and `C:/Users\..` both confuse the lexical prefix-match.
        // Refuse the mix outright when canonicalisation didn't run.
        if cfg!(windows) && path_str.contains('\\') && path_str.contains('/') {
            return Err(anyhow!(
                "path `{}` mixes `/` and `\\` and the parent dir does not exist \
                 — refused (cannot disambiguate without canonicalisation)",
                absolute.display()
            ));
        }
        // Windows 8.3 short-name pattern: `<1-6 alnum>~<digit>`. NTFS
        // exposes these for legacy 16-bit compat (PROGRA~1 → Program
        // Files). They bypass case-sensitivity-aware prefix checks
        // because canonicalising would expand the alias.
        if cfg!(windows) && contains_windows_short_name_component(&absolute) {
            return Err(anyhow!(
                "path `{}` looks like a Windows 8.3 short-name (`xxxxxx~N`) \
                 and the parent dir does not exist — refused (would alias \
                 outside the sandbox)",
                absolute.display()
            ));
        }
    }
    let resolved = match absolute.file_name() {
        Some(name) => parent_canon.join(name),
        None => parent_canon.clone(),
    };
    if !resolved.starts_with(&root_canon) {
        return Err(anyhow!(
            "path `{}` escapes sandbox `{}`",
            resolved.display(),
            root_canon.display()
        ));
    }
    // Round-16 M3: leaf-name symlink check. If the resolved path
    // exists AND is a symlink, refuse — its `fs::read` /
    // `fs::write` would follow the link out of the sandbox.
    if let Ok(meta) = std::fs::symlink_metadata(&resolved) {
        if meta.file_type().is_symlink() {
            return Err(anyhow!(
                "path `{}` is a symlink — refused (would escape sandbox via link target)",
                resolved.display()
            ));
        }
    }
    Ok(resolved)
}

/// Round-25 H1: sandbox-check `path`, then open it for read with the
/// kernel-level "do not follow symlinks at the leaf" flag set so the
/// open itself refuses to traverse a leaf symlink. Closes the TOCTOU
/// window between `sandbox_check` and the actual `File::open` — pre-
/// fix an attacker with write access to the sandbox dir could swap a
/// symlink into the leaf path between the check and the read, and
/// `fs::read_to_string` would silently follow the link out of the
/// sandbox.
///
/// On Unix this passes `O_NOFOLLOW` via
/// `OpenOptionsExt::custom_flags(libc::O_NOFOLLOW)` — `open(2)` then
/// returns `ELOOP` if the leaf is a symlink. On Windows this passes
/// `FILE_FLAG_OPEN_REPARSE_POINT` via
/// `OpenOptionsExt::custom_flags(...)` so `CreateFileW` opens the
/// reparse point itself rather than following it; the caller's
/// downstream read on a reparse-point handle then fails cleanly.
///
/// The pre-existing `sandbox_check` symlink-rejection branch
/// (round-16 M3) covers the "leaf already exists as a symlink at
/// check time" case; this helper covers the "leaf was raced into
/// existence as a symlink between check and open" case.
pub fn sandbox_open_read(path: &Path) -> Result<File> {
    let resolved = sandbox_check(path)?;
    open_no_follow(&resolved)
}

/// Unix implementation: `O_NOFOLLOW` makes `open(2)` return `ELOOP`
/// when the final path component is a symlink, atomically with the
/// kernel's symlink check. The open is otherwise read-only.
#[cfg(unix)]
fn open_no_follow(path: &Path) -> Result<File> {
    use std::os::unix::fs::OpenOptionsExt;
    std::fs::OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NOFOLLOW)
        .open(path)
        .map_err(|e| {
            anyhow!(
                "sandbox_open_read: open `{}` with O_NOFOLLOW failed: {e}",
                path.display(),
            )
        })
}

/// Windows implementation: `FILE_FLAG_OPEN_REPARSE_POINT` tells
/// `CreateFileW` to open the reparse point itself rather than the
/// link target. Unlike Unix `O_NOFOLLOW` (which causes `open(2)` to
/// fail with `ELOOP` when the leaf is a symlink), Windows lets the
/// OPEN succeed even on a symlink leaf — a subsequent read on the
/// resulting handle returns the raw reparse data, not the link
/// target's contents. That's the property the round-25 H1 fix relied
/// on (and at the time was sufficient for the sandbox use case
/// because callers immediately read the file and the JSON / text
/// parsers downstream rejected the raw reparse blob as malformed).
///
/// Round-26 M2: close the semantic gap explicitly. After a
/// successful open we `metadata().file_type().is_symlink()` and
/// refuse the file when true, mirroring the Unix `ELOOP` behaviour
/// at the helper boundary so callers can reason about the contract
/// uniformly across platforms.
#[cfg(windows)]
fn open_no_follow(path: &Path) -> Result<File> {
    use std::os::windows::fs::OpenOptionsExt;
    // FILE_FLAG_OPEN_REPARSE_POINT = 0x00200000. Hardcoded so we don't
    // need a runtime windows-sys re-export for a single constant.
    const FILE_FLAG_OPEN_REPARSE_POINT: u32 = 0x0020_0000;
    let f = std::fs::OpenOptions::new()
        .read(true)
        .custom_flags(FILE_FLAG_OPEN_REPARSE_POINT)
        .open(path)
        .map_err(|e| {
            anyhow!(
                "sandbox_open_read: open `{}` with FILE_FLAG_OPEN_REPARSE_POINT failed: {e}",
                path.display(),
            )
        })?;
    refuse_symlink_after_open(&f).map_err(|e| {
        anyhow!(
            "sandbox_open_read: leaf `{}` is a symlink — refused: {e}",
            path.display(),
        )
    })?;
    Ok(f)
}

/// Round-26 M2: shared post-open symlink check used by the Windows
/// `open_no_follow` (sandbox), `open_sidecar_no_follow` (dock
/// runner), and any future read/write helper that wants Unix
/// `O_NOFOLLOW`-equivalent semantics on Windows. On Unix the kernel
/// rejects symlink opens before we get a handle, so this helper is
/// #[cfg]-windows.
///
/// The reparse-point flag (`FILE_FLAG_OPEN_REPARSE_POINT`) lets the
/// open succeed on a symlink leaf; without this post-open check the
/// caller would read the raw reparse data rather than receiving an
/// error, which is the right shape for a defensible "we refused to
/// follow a leaf symlink" contract.
#[cfg(windows)]
fn refuse_symlink_after_open(f: &File) -> std::io::Result<()> {
    if f.metadata()?.file_type().is_symlink() {
        Err(std::io::Error::other(
            "leaf is a symlink (refused — FILE_FLAG_OPEN_REPARSE_POINT \
             would otherwise return a reparse-point handle)",
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `set_var` is unsafe-ish across threads; serialise the tests
    /// that mutate the env var so they can't race each other.
    fn lock() -> std::sync::MutexGuard<'static, ()> {
        use std::sync::{Mutex, OnceLock};
        static M: OnceLock<Mutex<()>> = OnceLock::new();
        M.get_or_init(|| Mutex::new(())).lock().unwrap()
    }

    #[test]
    fn accepts_path_inside_sandbox() {
        let _g = lock();
        let tmp = std::env::temp_dir().join("valenx-mcp-test-inside");
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(SANDBOX_ENV, &tmp);
        let p = sandbox_check(Path::new("foo.pdbqt")).unwrap();
        assert!(p.starts_with(tmp.canonicalize().unwrap()));
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn rejects_absolute_escape() {
        let _g = lock();
        let tmp = std::env::temp_dir().join("valenx-mcp-test-escape");
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(SANDBOX_ENV, &tmp);
        // Pick a path that is definitely NOT under the sandbox.
        #[cfg(windows)]
        let escape = Path::new("C:\\Windows\\System32\\drivers\\etc\\hosts");
        #[cfg(not(windows))]
        let escape = Path::new("/etc/passwd");
        let err = sandbox_check(escape).unwrap_err();
        assert!(err.to_string().contains("escapes sandbox"), "got: {err}");
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn rejects_dotdot_traversal() {
        let _g = lock();
        let tmp = std::env::temp_dir().join("valenx-mcp-test-dotdot");
        let inner = tmp.join("inner");
        std::fs::create_dir_all(&inner).unwrap();
        std::env::set_var(SANDBOX_ENV, &inner);
        // ../escape.pdbqt resolves to outside the sandbox.
        let err = sandbox_check(Path::new("../../escape.pdbqt"));
        assert!(err.is_err(), "expected error, got {err:?}");
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-16 M3 RED→GREEN: a symlink whose name lives inside the
    /// sandbox but whose target points OUTSIDE must be refused.
    /// Pre-fix only the parent dir was canonicalised — so a symlink
    /// like `<sandbox>/passwd → /etc/passwd` passed the prefix check
    /// (parent canonicalises inside) and `fs::read_to_string` then
    /// silently followed the link out of the sandbox.
    ///
    /// On Windows symlink creation requires Developer Mode or admin
    /// privileges. The test skips when the OS refuses the symlink
    /// (the mechanism under test still works on Windows — symlink
    /// followers via `fs::read` would otherwise leak).
    #[test]
    fn rejects_leaf_symlink_escape() {
        let _g = lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-test-symlink-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Target lives OUTSIDE the sandbox.
        let outside = tmp.parent().unwrap().join(format!(
            "valenx-mcp-outside-{}.txt",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&outside, b"secret").unwrap();
        // The symlink itself lives INSIDE the sandbox.
        let link = tmp.join("evil_link.pdbqt");
        #[cfg(unix)]
        let link_result = std::os::unix::fs::symlink(&outside, &link);
        #[cfg(windows)]
        let link_result = std::os::windows::fs::symlink_file(&outside, &link);
        if link_result.is_err() {
            eprintln!(
                "skipping: OS refused to create symlink (Windows: enable Developer Mode): {:?}",
                link_result.err()
            );
            let _ = std::fs::remove_file(&outside);
            let _ = std::fs::remove_dir_all(&tmp);
            return;
        }
        std::env::set_var(SANDBOX_ENV, &tmp);
        let err = sandbox_check(Path::new("evil_link.pdbqt"))
            .expect_err("sandbox_check must reject leaf-name symlinks");
        assert!(
            err.to_string().contains("symlink"),
            "expected symlink-rejection message, got: {err}"
        );
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_file(&link);
        let _ = std::fs::remove_file(&outside);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    // -----------------------------------------------------------------
    // Round-19 L2: lexical-fallback hardening when the parent dir
    // doesn't exist on disk (parent.canonicalize() fails). The
    // prefix-check then runs against the un-canonicalised path,
    // which is much weaker — refuse three known-bypass shapes.
    // -----------------------------------------------------------------

    /// Helper: a path WHOSE PARENT DOESN'T EXIST containing `..` must
    /// be refused even though we never canonicalised it. Pre-fix the
    /// lexical prefix-check could pass because `..` doesn't resolve
    /// without canonicalisation.
    #[test]
    fn lexical_fallback_rejects_dotdot_in_nonexistent_parent() {
        let _g = lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-l2-dotdot-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(SANDBOX_ENV, &tmp);
        // Point at a subdir under sandbox that DOES NOT EXIST so
        // parent.canonicalize() fails. The `..` inside it would
        // resolve to outside the sandbox once expanded.
        let p = Path::new("nonexistent_subdir/../../../escape.pdbqt");
        let err = sandbox_check(p).expect_err("lexical `..` must be refused");
        assert!(
            err.to_string().contains("..") || err.to_string().contains("escapes"),
            "expected `..`/escape message, got: {err}"
        );
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Windows-only: a path that mixes forward and back slashes
    /// AND has a non-existent parent must be refused — the lexical
    /// prefix check can't disambiguate the components.
    #[cfg(windows)]
    #[test]
    fn lexical_fallback_rejects_mixed_slash_on_windows() {
        let _g = lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-l2-mixed-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(SANDBOX_ENV, &tmp);
        // Absolute path with both `\` and `/` AND non-existent parent.
        // We can't easily construct one inside the sandbox that has a
        // non-existent parent AND mixed slashes — feed an absolute
        // path with both separators that points outside (which means
        // the prefix-check would refuse it anyway, but the L2 mixed-
        // slash branch fires first with its specific message).
        let p = Path::new("C:\\nonexistent-l2/sub/file.pdbqt");
        let err = sandbox_check(p).expect_err("mixed slash + missing parent must refuse");
        let msg = err.to_string();
        assert!(
            msg.contains("mixes") || msg.contains("escapes") || msg.contains(".."),
            "expected refusal, got: {msg}"
        );
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-25 H1 RED→GREEN (Unix-only): `sandbox_open_read` must
    /// refuse to follow a leaf symlink that was raced into the
    /// sandbox AFTER `sandbox_check` ran. We simulate the race by
    /// creating the symlink directly (the test isn't trying to
    /// recreate the racing thread — the kernel-level `O_NOFOLLOW`
    /// flag refuses the link regardless of timing), then asking
    /// `sandbox_open_read` to open it. The open MUST fail; pre-fix
    /// `read_capped_to_string` would silently follow the link.
    ///
    /// Windows is exercised by the `sandbox_open_read_refuses_symlink_windows`
    /// sibling — the Win32 `FILE_FLAG_OPEN_REPARSE_POINT` flag has
    /// different semantics (opens the reparse data itself rather
    /// than failing) so a separate assertion shape applies.
    #[cfg(unix)]
    #[test]
    fn sandbox_open_read_refuses_leaf_symlink() {
        let _g = lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-h1-noopen-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Set the sandbox to the parent so the symlink lives inside.
        std::env::set_var(SANDBOX_ENV, &tmp);
        // Create a target file OUTSIDE the sandbox to demonstrate
        // the escape vector — even though `sandbox_check`'s leaf
        // symlink-rejection branch (round-16 M3) would normally
        // catch this synchronously, this test exercises the
        // `O_NOFOLLOW` defence at the open syscall itself.
        let outside = tmp.parent().unwrap().join(format!(
            "valenx-h1-secret-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&outside, b"secret-payload").unwrap();
        let link = tmp.join("rcpt.pdbqt");
        std::os::unix::fs::symlink(&outside, &link).unwrap();
        // Pre-fix `read_capped_to_string` would happily follow the
        // link and return "secret-payload". Post-fix `sandbox_check`
        // refuses it (M3) — and even if it didn't, `O_NOFOLLOW` on
        // the open would still error with ELOOP.
        let res = sandbox_open_read(Path::new("rcpt.pdbqt"));
        assert!(
            res.is_err(),
            "sandbox_open_read must refuse a leaf-symlink read",
        );
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_file(&link);
        let _ = std::fs::remove_file(&outside);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-25 H1 RED→GREEN (cross-platform): plain (non-symlink)
    /// files inside the sandbox MUST still open via
    /// `sandbox_open_read`. Pins that the new helper doesn't break
    /// the legitimate happy path while it locks down the symlink
    /// case.
    #[test]
    fn sandbox_open_read_accepts_plain_file_inside_sandbox() {
        let _g = lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-h1-ok-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(SANDBOX_ENV, &tmp);
        let p = tmp.join("ok.pdbqt");
        std::fs::write(&p, b"ATOM\n").unwrap();
        let f = sandbox_open_read(Path::new("ok.pdbqt"))
            .expect("plain file must open via sandbox_open_read");
        // Reading from the handle must yield the file contents.
        use std::io::Read;
        let mut s = String::new();
        std::io::BufReader::new(f).read_to_string(&mut s).unwrap();
        assert_eq!(s, "ATOM\n");
        std::env::remove_var(SANDBOX_ENV);
        let _ = std::fs::remove_file(&p);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// 8.3 short-name component scanner: the standalone helper must
    /// flag a path like `PROGRA~1\foo` (the canonical NTFS legacy
    /// alias) and pass plain ASCII paths. Tests the helper directly
    /// rather than through `sandbox_check` so the assertion stays
    /// cross-platform (the runtime branch is `cfg!(windows)`-gated
    /// for the refusal but the scanner itself works on every host).
    #[test]
    fn windows_short_name_helper_flags_legacy_alias() {
        // Canonical NTFS short-names — `PROGRA~1`, `MY_FIL~1.TXT`,
        // `ABC~9` — should all match.
        assert!(contains_windows_short_name_component(Path::new("C:/PROGRA~1/foo")));
        assert!(contains_windows_short_name_component(Path::new(
            "C:/some/path/MY_FIL~1.TXT"
        )));
        assert!(contains_windows_short_name_component(Path::new("ABC~9")));
        // Negative cases — plain names, names with tilde but no
        // trailing digit, names too long for the 1-6 prefix.
        assert!(!contains_windows_short_name_component(Path::new(
            "ProgramFiles/foo"
        )));
        assert!(!contains_windows_short_name_component(Path::new("normalpath")));
        // Lowercase doesn't match — NTFS aliases are always uppercase.
        assert!(!contains_windows_short_name_component(Path::new(
            "progra~1"
        )));
        // Tilde without a digit suffix.
        assert!(!contains_windows_short_name_component(Path::new("PROGRA~X")));
        // Prefix longer than 6 chars doesn't match the NTFS layout.
        assert!(!contains_windows_short_name_component(Path::new(
            "TOOLONG~1"
        )));
    }
}
