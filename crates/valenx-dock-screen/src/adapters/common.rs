//! Shared scaffolding for the neural-network-tool subprocess adapters.
//!
//! Every adapter in [`crate::adapters`] needs the same two things:
//!
//! - a way to find an external binary on `PATH`
//!   ([`find_executable`]) — with a Windows `.exe` / `PATHEXT`
//!   fallback so conda / scoop / chocolatey shims are visible;
//! - a typed representation of the subprocess command line that
//!   *would* run the tool ([`AdapterCommand`]).
//!
//! No subprocess is ever launched from this module — see the
//! [`crate::adapters`] module docs for why. [`AdapterCommand`] is a
//! plain data structure; the Valenx job runner (or a test) decides
//! whether and when to execute it.

use std::path::PathBuf;

/// The availability of an external tool on the host.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ToolStatus {
    /// The tool's binary was found on `PATH`.
    Available {
        /// Absolute path of the located binary.
        path: PathBuf,
    },
    /// The tool's binary was not found.
    Missing,
}

impl ToolStatus {
    /// `true` if the tool is available.
    pub fn is_available(&self) -> bool {
        matches!(self, ToolStatus::Available { .. })
    }

    /// The located binary path, if the tool is available.
    pub fn path(&self) -> Option<&PathBuf> {
        match self {
            ToolStatus::Available { path } => Some(path),
            ToolStatus::Missing => None,
        }
    }
}

/// A fully-specified subprocess command line for an external tool.
///
/// This is what an adapter's `run_*` returns when the tool is present:
/// the program path plus its arguments. Executing it (a long-running,
/// often GPU-bound job) is the caller's responsibility — a library
/// function must not block on it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterCommand {
    /// Absolute path of the program to run.
    pub program: PathBuf,
    /// The program's command-line arguments, in order.
    pub args: Vec<String>,
    /// A human-readable description of what the command does.
    pub description: String,
}

impl AdapterCommand {
    /// Build an adapter command.
    pub fn new(
        program: impl Into<PathBuf>,
        args: Vec<String>,
        description: impl Into<String>,
    ) -> Self {
        AdapterCommand {
            program: program.into(),
            args,
            description: description.into(),
        }
    }

    /// Render the command as a single shell-style string (for logging
    /// / display). Arguments containing whitespace are quoted.
    pub fn to_display_string(&self) -> String {
        let mut s = self.program.display().to_string();
        for a in &self.args {
            s.push(' ');
            if a.contains(char::is_whitespace) {
                s.push('"');
                s.push_str(a);
                s.push('"');
            } else {
                s.push_str(a);
            }
        }
        s
    }

    /// Construct a [`std::process::Command`] from this adapter command.
    ///
    /// **This builds the command — it does not run it.** The caller
    /// (the Valenx job runner) decides when to spawn it. Provided so a
    /// caller does not have to re-thread `program` / `args` by hand.
    pub fn to_process_command(&self) -> std::process::Command {
        let mut cmd = std::process::Command::new(&self.program);
        cmd.args(&self.args);
        cmd
    }
}

/// Candidate filenames for `name` on the current platform.
///
/// On non-Windows: just `name`. On Windows: the bare name first (it
/// may already carry an extension), then `name + ext` for each
/// extension in `PATHEXT` — this catches the `.bat` / `.cmd` shims
/// conda / scoop / chocolatey produce.
fn platform_candidates(name: &str) -> Vec<String> {
    if !cfg!(windows) {
        return vec![name.to_string()];
    }
    let pathext = std::env::var("PATHEXT").unwrap_or_else(|_| ".COM;.EXE;.BAT;.CMD".to_string());
    let mut out = vec![name.to_string()];
    let name_lower = name.to_ascii_lowercase();
    let mut seen: Vec<String> = Vec::new();
    for ext in pathext.split(';').map(str::trim).filter(|s| !s.is_empty()) {
        let ext_lower = ext.to_ascii_lowercase();
        if seen.contains(&ext_lower) {
            continue;
        }
        seen.push(ext_lower.clone());
        if name_lower.ends_with(&ext_lower) {
            continue;
        }
        out.push(format!("{name}{ext_lower}"));
    }
    out
}

/// Locate an executable on `PATH`.
///
/// Tries each name in `names` (in order) against every `PATH`
/// directory, applying the platform candidate rules. Returns the
/// absolute path of the first hit, or [`ToolStatus::Missing`] if no
/// candidate exists.
///
/// This is the single PATH-probing primitive every adapter uses; it
/// mirrors `valenx_core::adapter_helpers::find_on_path` so the
/// behaviour is consistent with the rest of the Valenx adapter stack,
/// kept local here to avoid pulling `valenx-core` into this crate.
pub fn find_executable(names: &[&str]) -> ToolStatus {
    let Some(path) = std::env::var_os("PATH") else {
        return ToolStatus::Missing;
    };
    for dir in std::env::split_paths(&path) {
        for name in names {
            for candidate in platform_candidates(name) {
                let full = dir.join(&candidate);
                if full.is_file() {
                    return ToolStatus::Available { path: full };
                }
            }
        }
    }
    ToolStatus::Missing
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_executable_returns_missing_for_a_nonexistent_tool() {
        // A name that cannot possibly be on any PATH.
        let status = find_executable(&["valenx-no-such-tool-9q8w7e6r"]);
        assert_eq!(status, ToolStatus::Missing);
        assert!(!status.is_available());
        assert!(status.path().is_none());
    }

    #[test]
    fn tool_status_available_exposes_its_path() {
        let p = PathBuf::from("/usr/bin/example");
        let status = ToolStatus::Available { path: p.clone() };
        assert!(status.is_available());
        assert_eq!(status.path(), Some(&p));
    }

    #[test]
    fn adapter_command_display_quotes_whitespace_args() {
        let cmd = AdapterCommand::new(
            PathBuf::from("/usr/bin/tool"),
            vec!["--in".into(), "a path/with space".into()],
            "run the tool",
        );
        let s = cmd.to_display_string();
        assert!(s.contains("/usr/bin/tool"));
        assert!(s.contains("\"a path/with space\""), "got: {s}");
        assert!(s.contains("--in"));
    }

    #[test]
    fn adapter_command_to_process_command_carries_program_and_args() {
        let cmd = AdapterCommand::new(PathBuf::from("/usr/bin/tool"), vec!["--flag".into()], "x");
        let proc = cmd.to_process_command();
        // get_program is the program path.
        assert_eq!(proc.get_program(), std::ffi::OsStr::new("/usr/bin/tool"));
        let args: Vec<_> = proc.get_args().collect();
        assert_eq!(args, vec![std::ffi::OsStr::new("--flag")]);
    }

    #[test]
    fn platform_candidates_includes_bare_name() {
        // On every platform the bare name is a candidate.
        let cands = platform_candidates("alphafold");
        assert!(cands.contains(&"alphafold".to_string()));
        if cfg!(windows) {
            // On Windows there should be extra extension candidates.
            assert!(cands.len() > 1);
        } else {
            assert_eq!(cands.len(), 1);
        }
    }
}
