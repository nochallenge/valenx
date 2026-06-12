//! Cross-platform "show this directory in the host's file browser".
//! The user clicks "Open in file browser" anywhere the app shows a
//! workdir / artifact path; this is the single entry point that
//! routes through `explorer.exe` / `open` / `xdg-open`.

/// Cross-platform "show this directory in the host's file browser".
///
/// Picks the per-OS launcher: `explorer.exe` on Windows, `open` on
/// macOS, `xdg-open` on Linux/BSD. We spawn the launcher and detach;
/// the child is allowed to outlive us (no [`std::process::Child`]
/// kept around to wait on or kill).
///
/// Returns `Err(reason)` only when the launcher fails to spawn —
/// e.g. `xdg-open` isn't on a headless Linux box. Successfully
/// launching but the launcher then failing to actually open anything
/// (e.g. a misconfigured default app) is invisible from here; that's
/// the host's problem, not Valenx's.
/// Build (but do NOT spawn) the per-OS file-browser launcher command
/// for `path`. Split out from [`open_path_in_file_browser`] so the
/// command-construction logic is unit-testable WITHOUT ever spawning a
/// real launcher — spawning `explorer.exe` during `cargo test` pops a
/// user-visible window on Windows, which is exactly the popup-regression
/// class the project's kill-switch guards against.
///
/// Windows: `explorer.exe <file>` would open the file with its default
/// association (e.g. VS Code for .jsonl), which is wrong here — we want
/// to REVEAL the file in Explorer, not launch its editor.
/// `explorer.exe /select,<path>` highlights the file in its containing
/// folder. For directories we pass the path as-is so Explorer opens the
/// folder directly.
///
/// macOS: `open -R <file>` reveals + selects in Finder; `open <dir>`
/// opens the dir.
///
/// Linux/BSD: `xdg-open` doesn't have a reveal flag, so for a file path
/// we route through the parent directory instead — the user can't be
/// "popped into the file's editor" because xdg-open would do exactly
/// that.
fn build_file_browser_command(path: &std::path::Path) -> std::process::Command {
    if cfg!(target_os = "windows") {
        let mut c = std::process::Command::new("explorer.exe");
        if path.is_file() {
            let mut arg = std::ffi::OsString::from("/select,");
            arg.push(path.as_os_str());
            c.arg(arg);
        } else {
            c.arg(path.as_os_str());
        }
        c
    } else if cfg!(target_os = "macos") {
        let mut c = std::process::Command::new("open");
        if path.is_file() {
            c.arg("-R");
        }
        c.arg(path.as_os_str());
        c
    } else {
        let mut c = std::process::Command::new("xdg-open");
        let target: &std::path::Path = if path.is_file() {
            path.parent().unwrap_or(path)
        } else {
            path
        };
        c.arg(target.as_os_str());
        c
    }
}

pub fn open_path_in_file_browser(path: &std::path::Path) -> Result<(), String> {
    let mut cmd = build_file_browser_command(path);
    // Round-24 L5: `cmd.spawn().map(|_child| ())` drops the Child
    // immediately. On POSIX this leaves the child as a zombie until
    // the parent process exits — `Drop` for `std::process::Child`
    // does NOT reap (the docs say "the resources of the child
    // process are not released until the parent process explicitly
    // waits for the child"). For a short-lived launcher like `open`
    // / `xdg-open` the leak is small per invocation, but a session
    // that opens dozens of workdirs accumulates dozens of zombies
    // pinned to the GUI's PID table.
    //
    // Fix: hand the Child off to a detached background thread that
    // calls `wait()` so the OS reaps when the launcher exits. The
    // thread is fire-and-forget — we don't care about the launcher's
    // exit code (a failed `open` would have errored at spawn time;
    // a successful spawn + later launcher misconfiguration is the
    // host's problem). On Windows `Child::wait()` works the same
    // way; the thread costs ~4 KiB of stack and dies as soon as the
    // launcher exits (usually within a second).
    match cmd.spawn() {
        Ok(child) => {
            std::thread::spawn(move || {
                let mut child = child;
                let _ = child.wait();
            });
            Ok(())
        }
        Err(e) => Err(format!("spawn file browser: {e}")),
    }
}

/// Sentinel `Err(_)` returned by [`open_path_or_copy`] when the user
/// has opted in to the global "no file-browser popups" kill-switch.
/// Callers detect this with [`str::contains`] and surface it as a
/// neutral status message rather than a red `last_error` — the popup
/// suppression is the user's affirmative choice, not a failure.
pub const POPUP_DISABLED_PREFIX: &str = "File browser disabled — path: ";

/// Kill-switched wrapper around [`open_path_in_file_browser`]: when
/// `disable_popups` is true, this becomes a no-op that returns
/// `Err(POPUP_DISABLED_PREFIX + path)`. The string-error contract lets
/// the existing call sites — which already format `last_error` /
/// `status` from the launcher's `Err(reason)` — surface the path to
/// the user without growing a new return type.
///
/// **Why a status-message fallback rather than a clipboard copy?**
/// `arboard` is in the transitive dep graph via `egui-winit`, but
/// pulling it in as a *direct* dep just to copy a path here would
/// turn off a Cargo feature flag's deduplication safety (the
/// egui-winit selection could be re-resolved) and add a sync /
/// init step (`Clipboard::new()` is fallible on headless Linux,
/// CI, locked-screen sessions). The existing pattern in
/// `run_actions.rs` already includes the path verbatim in the
/// error message ("Couldn't open file browser: ... The workdir is
/// at `<path>`"), so a status-line path is consistent with what
/// users already see when the launcher fails. If a future revision
/// wants real clipboard integration, swap the `Err(...)` branch
/// here for an `arboard::Clipboard::set_text` call — no call-site
/// changes needed.
pub fn open_path_or_copy(path: &std::path::Path, disable_popups: bool) -> Result<(), String> {
    if disable_popups {
        return Err(format!("{POPUP_DISABLED_PREFIX}{}", path.display()));
    }
    open_path_in_file_browser(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// RED→GREEN: when `disable_popups = true`, the helper MUST return
    /// `Err(_)` with the kill-switch prefix — no `explorer.exe` /
    /// `open` / `xdg-open` spawn under any circumstance. Pre-fix,
    /// every call site spawned unconditionally; this test pins the
    /// new kill-switch behaviour so a future refactor that drops the
    /// `if disable_popups` branch fails here instead of silently
    /// re-popping File Explorer on users who opted out.
    #[test]
    fn open_path_or_copy_with_disable_does_not_spawn_and_signals_disabled() {
        let path = std::path::Path::new("/tmp/valenx-test-no-such-path");
        let result = open_path_or_copy(path, true);
        let err = result.expect_err("disable_popups=true must return Err");
        assert!(
            err.starts_with(POPUP_DISABLED_PREFIX),
            "expected the disabled-sentinel prefix, got: {err}"
        );
        // The path must round-trip through the error message so the
        // call site can surface it to the user as a status line.
        assert!(
            err.contains("valenx-test-no-such-path"),
            "expected the path to round-trip into the error, got: {err}"
        );
    }

    /// R29 C1: assert the per-OS launcher command is constructed
    /// correctly WITHOUT ever spawning it. The previous test here
    /// (`open_path_or_copy_with_disable_off_does_not_use_disabled_prefix`)
    /// called `open_path_or_copy(path, false)` which delegates to
    /// `open_path_in_file_browser` → `cmd.spawn()` — on Windows that
    /// REALLY launched `explorer.exe` during `cargo test`, popping a
    /// visible Explorer window. That is the exact popup-regression
    /// the kill-switch exists to prevent, so the spawning test is
    /// gone; we now inspect the built [`std::process::Command`] via
    /// `get_program()` / `get_args()` (stable since Rust 1.57) and
    /// never spawn anything.
    #[cfg(target_os = "windows")]
    #[test]
    fn build_command_selects_reveal_flag_per_os() {
        let dir = std::env::temp_dir().join("valenx-fb-test-dir-c1");
        std::fs::create_dir_all(&dir).expect("create tempdir");

        // Directory branch: program == explorer.exe, args == [<dir>].
        let cmd = build_file_browser_command(&dir);
        assert_eq!(cmd.get_program(), std::ffi::OsStr::new("explorer.exe"));
        let args: Vec<_> = cmd.get_args().collect();
        assert_eq!(args.len(), 1, "dir branch passes exactly one arg");
        assert_eq!(args[0], dir.as_os_str());

        // File branch: a real temp file → first arg starts `/select,`.
        let file = dir.join("artifact.jsonl");
        std::fs::write(&file, b"{}").expect("write temp file");
        let cmd = build_file_browser_command(&file);
        assert_eq!(cmd.get_program(), std::ffi::OsStr::new("explorer.exe"));
        let args: Vec<_> = cmd.get_args().collect();
        assert_eq!(args.len(), 1, "file branch passes exactly one arg");
        let arg0 = args[0].to_string_lossy();
        assert!(
            arg0.starts_with("/select,"),
            "file branch must reveal via /select, got: {arg0}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[cfg(target_os = "macos")]
    #[test]
    fn build_command_selects_reveal_flag_per_os() {
        let dir = std::env::temp_dir().join("valenx-fb-test-dir-c1");
        std::fs::create_dir_all(&dir).expect("create tempdir");

        // Directory branch: program == open, single arg <dir>.
        let cmd = build_file_browser_command(&dir);
        assert_eq!(cmd.get_program(), std::ffi::OsStr::new("open"));
        let args: Vec<_> = cmd.get_args().collect();
        assert_eq!(args.len(), 1, "dir branch passes exactly one arg");
        assert_eq!(args[0], dir.as_os_str());

        // File branch: `open -R <file>` reveals in Finder.
        let file = dir.join("artifact.jsonl");
        std::fs::write(&file, b"{}").expect("write temp file");
        let cmd = build_file_browser_command(&file);
        assert_eq!(cmd.get_program(), std::ffi::OsStr::new("open"));
        let args: Vec<_> = cmd.get_args().collect();
        assert_eq!(
            args.first().map(|a| a.to_string_lossy()).as_deref(),
            Some("-R")
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[cfg(not(any(target_os = "windows", target_os = "macos")))]
    #[test]
    fn build_command_selects_reveal_flag_per_os() {
        let dir = std::env::temp_dir().join("valenx-fb-test-dir-c1");
        std::fs::create_dir_all(&dir).expect("create tempdir");

        // Directory branch: program == xdg-open, single arg <dir>.
        let cmd = build_file_browser_command(&dir);
        assert_eq!(cmd.get_program(), std::ffi::OsStr::new("xdg-open"));
        let args: Vec<_> = cmd.get_args().collect();
        assert_eq!(args.len(), 1, "dir branch passes exactly one arg");
        assert_eq!(args[0], dir.as_os_str());

        // File branch: xdg-open has no reveal flag, so we route to the
        // parent directory instead of the file itself.
        let file = dir.join("artifact.jsonl");
        std::fs::write(&file, b"{}").expect("write temp file");
        let cmd = build_file_browser_command(&file);
        assert_eq!(cmd.get_program(), std::ffi::OsStr::new("xdg-open"));
        let args: Vec<_> = cmd.get_args().collect();
        assert_eq!(args.len(), 1, "file branch routes to parent, one arg");
        assert_eq!(args[0], dir.as_os_str(), "file branch uses parent dir");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
