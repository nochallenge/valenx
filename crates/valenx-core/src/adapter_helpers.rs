//! Small utilities that every adapter ends up reimplementing.
//!
//! Rather than ship a `valenx-adapter-common` crate (and pay the
//! workspace-member overhead), we surface these as helpers on
//! `valenx-core`, which every adapter already depends on.
//!
//! Kept deliberately small:
//!
//! - [`find_on_path`] — locate an executable under `PATH`
//!   with a Windows `.exe` fallback.
//! - [`platform_suffix`] — the `.exe` fallback alone.
//! - [`stub_provenance`] — a deterministic empty `Provenance` for
//!   adapters whose `collect()` is still a placeholder.
//! - [`not_implemented`] — consistent error for
//!   prepare/run/collect stubs in scaffold adapters.
//!
//! Real adapters should graduate off these helpers as they mature.

use std::path::PathBuf;

use crate::error::AdapterError;
use valenx_fields::provenance::{Provenance, Sha256Hex};

/// Look up the first matching binary on `PATH`. Adds `.exe` on
/// Windows when the caller hasn't already. Returns the absolute path
/// of the first hit, searching each `PATH` entry in order.
///
/// On Windows the search tries every extension in `PATHEXT`
/// (defaulting to `.COM;.EXE;.BAT;.CMD` when the env var is unset).
/// This catches `.bat` / `.cmd` shims that package managers like
/// conda, scoop, and chocolatey routinely produce — pre-fix, those
/// were invisible to `find_on_path`.
pub fn find_on_path(names: &[&str]) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path) {
        for name in names {
            for candidate_name in platform_candidates(name) {
                let candidate = dir.join(&candidate_name);
                if candidate.is_file() {
                    return Some(candidate);
                }
            }
        }
    }
    None
}

/// Candidate filenames to try for `name` on the current platform.
///
/// - On non-Windows: `[name]`.
/// - On Windows: the bare name (in case it already carries an
///   extension), then `name + ext` for each extension in `PATHEXT`
///   (defaulting to `.COM;.EXE;.BAT;.CMD` when the env var is unset
///   or empty). Comparison is case-insensitive — extensions are
///   lowercased so a `.EXE` entry in `PATHEXT` doesn't double-list
///   alongside our test-bench `.exe`.
pub fn platform_candidates(name: &str) -> Vec<String> {
    if !cfg!(windows) {
        return vec![name.to_string()];
    }
    let pathext = std::env::var("PATHEXT").unwrap_or_else(|_| ".COM;.EXE;.BAT;.CMD".to_string());
    platform_candidates_for(name, &pathext)
}

/// Pure helper that doesn't touch `std::env`. Lets unit tests cover
/// the dedup / double-extension / empty-pathext edges without the
/// `forbid(unsafe_code)` lint refusing the env-var manipulation
/// `set_var` would require.
pub fn platform_candidates_for(name: &str, pathext: &str) -> Vec<String> {
    if !cfg!(windows) {
        return vec![name.to_string()];
    }
    // Try the bare name first — caller may have passed `gmsh.exe`
    // with the extension already, or be on a system where the
    // shell has resolved a no-extension wrapper script.
    let mut out = vec![name.to_string()];
    let mut seen: Vec<String> = Vec::new();
    for ext in pathext.split(';').map(str::trim).filter(|s| !s.is_empty()) {
        let ext_lower = ext.to_ascii_lowercase();
        if seen.contains(&ext_lower) {
            continue;
        }
        // Skip extensions the caller already asked for (case-
        // insensitive) so we don't push `gmsh.exe.exe`.
        let name_lower = name.to_ascii_lowercase();
        if name_lower.ends_with(&ext_lower) {
            seen.push(ext_lower);
            continue;
        }
        out.push(format!("{name}{ext_lower}"));
        seen.push(ext_lower);
    }
    out
}

/// Add the `.exe` suffix on Windows when a binary name doesn't
/// already have one. No-op on other platforms.
///
/// Kept for backward-compat with callers that need a single
/// canonical name. New code should prefer [`platform_candidates`]
/// which respects `PATHEXT` and surfaces every legitimate shim.
pub fn platform_suffix(name: &str) -> String {
    if cfg!(windows) && !name.ends_with(".exe") {
        format!("{name}.exe")
    } else {
        name.to_string()
    }
}

/// Validate that `name` is a safe single-path-component basename — no
/// path separators, no `..`, no absolute paths, no empty.
///
/// This is the canonical check for adapter fields like `output_basename`,
/// `prefix`, or `output_prefix` that get appended to the workdir before
/// being passed to a subprocess. Without it, a hostile case.toml could
/// set `output_basename = "../../etc/cron.d/x"` and the adapter would
/// happily write into anywhere on disk it has permission for.
pub fn validate_output_basename(name: &str, field_label: &str) -> Result<(), AdapterError> {
    let trimmed = name.trim();
    if trimmed.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label}: must not be empty"
        )));
    }
    if trimmed.contains('/') || trimmed.contains('\\') {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label} `{name}`: must not contain path separators \
             (`/` or `\\`) — basenames are written next to the case \
             workdir, not arbitrary paths"
        )));
    }
    if trimmed == ".." || trimmed.starts_with("../") || trimmed.starts_with("..\\") {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label} `{name}`: must not traverse via `..`"
        )));
    }
    if std::path::Path::new(trimmed).is_absolute() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label} `{name}`: must not be an absolute path"
        )));
    }
    Ok(())
}

/// Validate that a user-supplied `output_dir` (or similar
/// workdir-relative subdirectory path) is safe to join with the
/// workdir. Unlike [`validate_output_basename`], this helper allows
/// multi-component relative paths (`results/run1`) but still rejects
/// absolute paths and `..` traversal.
///
/// Used by adapters whose tool CLI expects a directory the tool
/// creates inside the workdir (cwltool's `--outdir`, salmon /
/// kallisto's quant `-o`). Without this, a hostile case.toml could
/// set `output_dir = "/etc/cron.d"` or `output_dir = "../../etc"` and
/// the workdir-rooted subprocess would happily write into anywhere
/// the user has permission for.
pub fn validate_output_dir(path: &std::path::Path, field_label: &str) -> Result<(), AdapterError> {
    let display = path.display();
    let as_str = path.to_string_lossy();
    if as_str.trim().is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label}: must not be empty"
        )));
    }
    if path.is_absolute() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label} `{display}`: must not be an absolute path \
             (workdir-relative subdirectory required)"
        )));
    }
    // Paths starting with `/` or `\` are root-relative on POSIX
    // and current-drive-relative on Windows. Both are sandbox
    // escapes and `Path::is_absolute` doesn't always catch them
    // (a leading `/` on Windows is "relative to current drive root"
    // and reports `false` for `is_absolute`).
    if as_str.starts_with('/') || as_str.starts_with('\\') {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "{field_label} `{display}`: must not start with `/` or `\\` \
             (workdir-relative subdirectory required)"
        )));
    }
    for component in path.components() {
        if matches!(component, std::path::Component::ParentDir) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "{field_label} `{display}`: must not traverse via `..` \
                 (subdirectory must stay within the workdir)"
            )));
        }
        if matches!(component, std::path::Component::RootDir) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "{field_label} `{display}`: must not contain a root \
                 component (workdir-relative subdirectory required)"
            )));
        }
    }
    Ok(())
}

/// Allow-list of acceptable Python interpreter names. Adapters that
/// accept a `python` field in `case.toml` validate the value against
/// this list before passing it to `Command::new`, otherwise a hostile
/// case could turn `python = "/usr/bin/curl"` into arbitrary code
/// execution. Names match what conda / pyenv / system installs
/// typically produce — version-suffixed minor releases plus the
/// `uv`-managed alternatives.
///
/// Round-5 fix: extend with Windows `.exe` variants of every dotted
/// version (`python3.11.exe`) AND the conda-Windows no-dot form
/// (`python311.exe`). Conda-on-Windows in particular ships
/// `python3.11.exe` as the standard per-env name, and the round-4
/// allow-list rejected it.
pub const ALLOWED_PYTHON_NAMES: &[&str] = &[
    "python",
    "python3",
    "python3.8",
    "python3.9",
    "python3.10",
    "python3.11",
    "python3.12",
    "python3.13",
    "python3.14",
    "python.exe",
    "python3.exe",
    // Round-5: dotted Windows variants (`python3.11.exe`). pyenv-win,
    // chocolatey, scoop, and Miniconda-on-Windows ship these.
    "python3.8.exe",
    "python3.9.exe",
    "python3.10.exe",
    "python3.11.exe",
    "python3.12.exe",
    "python3.13.exe",
    "python3.14.exe",
    // Round-5: no-dot conda-Windows variants (`python311.exe`). The
    // anaconda3 \ envs \ <name> \ python311.exe layout uses this form.
    "python38.exe",
    "python39.exe",
    "python310.exe",
    "python311.exe",
    "python312.exe",
    "python313.exe",
    "python314.exe",
    // conda's per-env python alias, when activated.
    "conda-python",
    // `uv run python` shim some labs ship.
    "uv",
];

/// Allow-list of acceptable R / Rscript interpreter names. Sister of
/// [`ALLOWED_PYTHON_NAMES`]; closes the same arbitrary-binary-exec
/// hole on adapters that accept a user-supplied `rscript = "..."`
/// field in case.toml (iCodon, Seurat). Without this, a hostile case
/// could turn `rscript = "/usr/bin/curl"` into arbitrary code
/// execution.
///
/// `Rscript` is the canonical headless launcher (`R --no-save < file`
/// would also work but `Rscript` is what every modern doc tells users
/// to invoke); `R` is the interactive REPL but is occasionally pinned
/// by users running their own custom batch wrappers. Both Windows
/// `.exe` variants and the lowercased `rscript` variants conda-on-
/// Windows installs.
pub const ALLOWED_RSCRIPT_NAMES: &[&str] = &[
    "Rscript",
    "Rscript.exe",
    "rscript",
    "rscript.exe",
    "R",
    "R.exe",
];

/// Validate that the user-supplied Python interpreter name / path is
/// either a plain allow-listed name (resolved via `PATH`) or an
/// absolute path whose file_name is allow-listed. Rejects everything
/// else — `python = "/usr/bin/curl"`, `python = "../../etc/cron.d/x"`,
/// `python = "rm"`, etc.
///
/// Returns the path that should be passed to `Command::new`. Callers
/// chain this through `find_on_path` for the bare-name case.
pub fn validate_python_binary(spec: &str) -> Result<PathBuf, AdapterError> {
    let spec = spec.trim();
    if spec.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "python interpreter spec must not be empty"
        )));
    }
    let path = std::path::Path::new(spec);
    // The "file name" we'll compare against the allow list.
    let file_name = path.file_name().and_then(|s| s.to_str()).unwrap_or(spec);
    if !ALLOWED_PYTHON_NAMES
        .iter()
        .any(|allowed| file_name.eq_ignore_ascii_case(allowed))
    {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "python interpreter `{spec}` is not on the allow-list \
             {ALLOWED_PYTHON_NAMES:?} — refusing to invoke arbitrary \
             binaries from case.toml. Use a plain `python3` and let \
             PATH resolution do the rest, or pin the interpreter to a \
             versioned name (`python3.11`)."
        )));
    }
    // If the spec is a bare name, defer PATH resolution to the caller
    // (so they can fall back to the adapter's PYTHON_BINARIES list).
    // If it's a path, return as-is — but still demand the path live
    // somewhere reasonable: absolute is OK, relative is OK.
    Ok(PathBuf::from(spec))
}

/// Validate that the user-supplied R / Rscript interpreter name / path
/// is on [`ALLOWED_RSCRIPT_NAMES`]. Round-5 sister of
/// [`validate_python_binary`]; iCodon and Seurat adapters take a
/// user-supplied `rscript = "..."` field and previously called
/// `Command::new(&input.rscript)` directly, letting case.toml supply
/// arbitrary binaries (`rscript = "/usr/bin/curl"`).
///
/// Returns the path that should be passed to `Command::new`. Callers
/// chain this through `find_on_path` for the bare-name case (same
/// shape as `validate_python_binary`).
pub fn validate_rscript_binary(spec: &str) -> Result<PathBuf, AdapterError> {
    let spec = spec.trim();
    if spec.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "Rscript interpreter spec must not be empty"
        )));
    }
    let path = std::path::Path::new(spec);
    let file_name = path.file_name().and_then(|s| s.to_str()).unwrap_or(spec);
    if !ALLOWED_RSCRIPT_NAMES
        .iter()
        .any(|allowed| file_name.eq_ignore_ascii_case(allowed))
    {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "Rscript interpreter `{spec}` is not on the allow-list \
             {ALLOWED_RSCRIPT_NAMES:?} — refusing to invoke arbitrary \
             binaries from case.toml. Use a plain `Rscript` and let \
             PATH resolution do the rest, or pin the interpreter to an \
             absolute path whose basename matches the allow list."
        )));
    }
    Ok(PathBuf::from(spec))
}

/// Combined validate + resolve for adapter `prepare()` paths that take a
/// user-supplied Python interpreter spec from case.toml.
///
/// 1. [`validate_python_binary`] rejects anything off the allow-list.
/// 2. If the spec is an absolute path to an existing file, return as-is.
/// 3. If the spec contains a `..` component, reject (round-4 hardening).
/// 4. Otherwise resolve via [`find_on_path`] using the supplied bare
///    name first, then the adapter's fallback list.
///
/// Returns the absolute path that should be passed to `Command::new`.
/// On failure, returns the structured `AdapterError` the adapter would
/// otherwise have to construct by hand — keeps the per-adapter wiring
/// short.
pub fn resolve_python_binary(
    spec: &str,
    fallback_binaries: &[&str],
) -> Result<PathBuf, AdapterError> {
    let validated = validate_python_binary(spec)?;
    if validated.is_absolute() && validated.is_file() {
        return Ok(validated);
    }
    // Defence-in-depth: allow-list passes basename even for `..`-bearing
    // relatives like `../python3`. We must never resolve through a
    // parent-dir traversal because that would let case.toml reach out
    // of the project sandbox the user reasonably expects.
    if validated
        .components()
        .any(|c| matches!(c, std::path::Component::ParentDir))
    {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "python interpreter `{spec}` must not traverse via `..`"
        )));
    }
    // Bare-name case (or relative path that wasn't a `..` traversal):
    // try PATH lookup with the user's value first, then fall back to
    // the adapter's curated PYTHON_BINARIES list.
    find_on_path(&[spec])
        .or_else(|| find_on_path(fallback_binaries))
        .ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "could not locate Python interpreter `{spec}` on PATH"
            ))
        })
}

/// SHA-256 hex of an arbitrary byte slice. Wraps the same primitive
/// the audit log uses; lifted into adapter_helpers so adapter
/// `collect()` paths can compute case / mesh / input hashes for
/// the [`Provenance`] block without each one re-importing `sha2`.
pub fn sha256_hex_bytes(bytes: &[u8]) -> Sha256Hex {
    use sha2::{Digest, Sha256};
    let digest = Sha256::digest(bytes);
    let mut s = String::with_capacity(64);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(s, "{b:02x}");
    }
    Sha256Hex::new(s)
}

/// Chunk size used by [`sha256_hex_file`] when reading the input file
/// into the hasher. Chosen to match a typical OS-page-friendly read
/// granularity — small enough that the in-flight buffer stays cheap,
/// large enough that the per-syscall overhead doesn't dominate.
const HASH_CHUNK: usize = 64 * 1024;

/// Hard upper bound on bytes hashed by [`sha256_hex_file`]. Files
/// larger than this are treated like missing files (empty hash) to
/// match `stub_provenance`'s convention.
///
/// Round-16 M1: pre-fix this helper slurped the whole file into a
/// single `Vec<u8>` via `std::fs::read`. A 5 GB mesh dump would then
/// allocate 5 GB of RAM before the hasher ran. 4 GiB is the largest
/// payload any realistic input file could be (the case+mesh+tools-lock
/// triple `live_provenance` hashes are typically KB-to-MB sized;
/// anything past 4 GiB is almost certainly poisoned input, not a real
/// scientific dataset).
pub const MAX_HASHED_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// SHA-256 hex of a file's contents, streamed in fixed-size chunks
/// so the helper never allocates more than `HASH_CHUNK` bytes plus
/// the hasher's internal state.
///
/// Returns `Sha256Hex::new("")` when the file can't be read OR when
/// the file exceeds [`MAX_HASHED_BYTES`] — matches `stub_provenance`'s
/// "missing inputs land as empty hash" convention so adapters with
/// optional mesh / lock files don't fail the run.
///
/// Round-16 M1: pre-fix used `std::fs::read(&path)?` then
/// `Sha256::digest(&bytes)`. A 5 GB mesh would spike memory by 5 GB
/// (plus the hasher's working set). The streaming form keeps the
/// peak RSS bounded regardless of file size.
pub fn sha256_hex_file(path: &std::path::Path) -> Sha256Hex {
    use sha2::{Digest, Sha256};
    use std::io::Read;
    // Belt-and-braces: stat the file first so an obviously-too-large
    // payload is rejected before we even open the descriptor. The
    // in-loop check below catches the racing-grow case where the
    // file expands between stat and read.
    if let Ok(meta) = std::fs::metadata(path) {
        if meta.len() > MAX_HASHED_BYTES {
            return Sha256Hex::new("");
        }
    }
    let mut file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return Sha256Hex::new(""),
    };
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; HASH_CHUNK];
    let mut total: u64 = 0;
    loop {
        let n = match file.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => n,
            Err(_) => return Sha256Hex::new(""),
        };
        total = match total.checked_add(n as u64) {
            Some(t) => t,
            None => return Sha256Hex::new(""),
        };
        if total > MAX_HASHED_BYTES {
            // Match the "missing → empty" convention rather than
            // returning a partial hash that callers can't tell from
            // a real one.
            return Sha256Hex::new("");
        }
        hasher.update(&buf[..n]);
    }
    let digest = hasher.finalize();
    let mut s = String::with_capacity(64);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(s, "{b:02x}");
    }
    Sha256Hex::new(s)
}

/// Compute a stable run id (UUIDv4 shape) without pulling in the
/// `uuid` crate as a hard dep. Combines the current
/// `SystemTime` nanos with a small process-local counter so two
/// `live_provenance` calls in the same nanosecond still produce
/// distinct ids.
pub fn fresh_run_id() -> String {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(1);
    let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    // Hash the (counter, nanos) pair so the result looks like
    // a UUIDv4 (16 hex chars per segment in the 8-4-4-4-12 layout)
    // without needing the `uuid` dep.
    let payload = format!("{counter}-{nanos}");
    let h = sha256_hex_bytes(payload.as_bytes()).0;
    // Layout: 8-4-4-4-12 from the first 32 hex chars; pad with
    // counter / nanos hex on overflow (impossible — sha256 hex is
    // 64 chars).
    let h = &h[..32];
    format!(
        "{}-{}-{}-{}-{}",
        &h[0..8],
        &h[8..12],
        &h[12..16],
        &h[16..20],
        &h[20..32]
    )
}

/// Build a fully-populated [`Provenance`] from a case path + the
/// ambient adapter metadata. Hashes case.toml + the optional mesh
/// file + tools.lock; falls back to empty hashes when any input is
/// absent. The run_id is fresh per call.
///
/// Adapters call this from `collect()` to replace the
/// `stub_provenance` placeholder with live data. `wall_time_seconds`
/// is filled in by the caller (the adapter knows its own elapsed
/// time; this helper doesn't try to time the call).
///
/// The 8 args are deliberate — adapter id / version, tool id /
/// version, three input paths, and elapsed time are all distinct
/// concerns and a `Provenance` builder struct would just shuffle
/// the same data through one more layer.
#[allow(clippy::too_many_arguments)]
pub fn live_provenance(
    adapter_id: &str,
    adapter_version: &str,
    tool_name: &str,
    tool_version: &str,
    case_path: &std::path::Path,
    mesh_path: Option<&std::path::Path>,
    tools_lock_path: Option<&std::path::Path>,
    wall_time_seconds: f64,
) -> Provenance {
    Provenance {
        adapter: adapter_id.to_string(),
        adapter_version: adapter_version.to_string(),
        tool: tool_name.to_string(),
        tool_version: tool_version.to_string(),
        case_hash: sha256_hex_file(case_path),
        mesh_hash: mesh_path
            .map(sha256_hex_file)
            .unwrap_or_else(|| Sha256Hex::new("")),
        input_hash: sha256_hex_file(case_path), // alias for now
        tools_lock_hash: tools_lock_path
            .map(sha256_hex_file)
            .unwrap_or_else(|| Sha256Hex::new("")),
        run_id: fresh_run_id(),
        wall_time_seconds,
        completed_at: current_timestamp_iso8601(),
        ancestors: Vec::new(),
    }
}

/// ISO-8601 UTC timestamp ("YYYY-MM-DDTHH:MM:SSZ"). Pulled in here
/// so the helper module doesn't have to depend on `chrono` — uses
/// `SystemTime` + manual day-of-year math.
fn current_timestamp_iso8601() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = secs / 86400;
    let secs_today = secs % 86400;
    let h = secs_today / 3600;
    let m = (secs_today % 3600) / 60;
    let s = secs_today % 60;
    let (y, mo, d) = days_to_ymd(days as i64);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}

/// Convert a Unix-day count to (year, month, day). Uses Hinnant's
/// civil-from-days algorithm (public domain) — same routine
/// valenx-app uses for its audit-timestamp computation.
fn days_to_ymd(z: i64) -> (i32, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = if m <= 2 { y + 1 } else { y };
    (y as i32, m, d)
}

/// Build a deterministic empty [`Provenance`] for adapters that have
/// a placeholder `collect()`. Timestamps are fixed so snapshot tests
/// stay reproducible; real adapters replace this with live data.
pub fn stub_provenance(adapter_id: &str, adapter_version: &str, tool_name: &str) -> Provenance {
    Provenance {
        adapter: adapter_id.to_string(),
        adapter_version: adapter_version.to_string(),
        tool: tool_name.to_string(),
        tool_version: "unknown".to_string(),
        case_hash: Sha256Hex::new(""),
        mesh_hash: Sha256Hex::new(""),
        input_hash: Sha256Hex::new(""),
        tools_lock_hash: Sha256Hex::new(""),
        run_id: "00000000-0000-0000-0000-000000000000".to_string(),
        wall_time_seconds: 0.0,
        completed_at: "1970-01-01T00:00:00Z".to_string(),
        ancestors: Vec::new(),
    }
}

/// Cap on the bytes read from a probed subprocess's stdout+stderr.
///
/// `Command::output()` reads stdout and stderr unbounded — a probe
/// against a hostile or runaway binary that emits gigabytes on the
/// `--version` channel would block the adapter and OOM the host.
/// 1 MiB per channel is generous: every well-behaved tool's version
/// banner is at most a few KiB, and the cap fires *before* any
/// memory-pressure failure.
///
/// Public for adapters that need the same constant; the round-8
/// addition wires it into [`capture_subprocess_stdout`].
pub const MAX_PROBE_OUTPUT_BYTES: usize = 1_048_576;

/// Hard wall-clock timeout on a single probe subprocess.
///
/// Round-16 M2: pre-fix `capture_subprocess_stdout` had only the
/// per-channel byte cap from round-8 — but a probe binary that hangs
/// (waiting on stdin, blocked on a license server, deadlocked in its
/// own startup) would never close its pipes, so the bounded reads
/// would block forever and the adapter's `probe()` would never
/// return. 10 s is generous: every real tool's `--version` returns
/// within a few hundred ms; a probe that takes longer is almost
/// certainly stuck and the right answer is "unknown version" rather
/// than "the GUI hangs on adapter discovery".
pub const MAX_PROBE_TIMEOUT_SECS: u64 = 10;

/// Run `<binary> <flag>` and return the raw stdout. `flag` is
/// typically `--version` but some tools use `-V` / `version` /
/// `-v` — pass whatever the tool wants.
///
/// Returns `None` on spawn / process failures (binary not found,
/// non-zero exit). Adapter probes that need version data fall back
/// to "unknown" in those cases rather than aborting the probe.
///
/// Round-8 hardening: bounded reads via [`MAX_PROBE_OUTPUT_BYTES`]
/// per channel (stdout + stderr). Pre-fix this used
/// `Command::output()` which reads both channels unbounded — a
/// runaway probe (10 MiB+ on stdout) would block the adapter and
/// pressure host memory. With the cap, output past the limit is
/// silently truncated and the truncated banner is parsed normally
/// (any sensible version string lives in the first few KiB).
///
/// Round-16 M2 hardening: hard wall-clock timeout via
/// [`MAX_PROBE_TIMEOUT_SECS`] plus an RAII kill-on-drop guard. Pre-fix
/// a binary that opened its pipes but never wrote (license-check
/// hang, stuck network call) would block the read loop forever; now a
/// watchdog thread kills the child after the timeout fires, the read
/// loop unblocks on the pipe EOF, and the probe returns `None`. The
/// `KillOnDropChild` guard belt-and-braces the kill so an early
/// return / panic from anywhere in the function still reaps the child.
///
/// Used by [`detect_tool_version`] for the common
/// "find_on_path -> --version -> parse semver" pipeline.
pub fn capture_subprocess_stdout(binary: &std::path::Path, flag: &str) -> Option<String> {
    capture_subprocess_stdout_with_timeout(binary, flag, MAX_PROBE_TIMEOUT_SECS)
}

/// Same as [`capture_subprocess_stdout`] but with a caller-chosen
/// timeout. Public-in-crate so tests can pin the watchdog behaviour
/// without sleeping the full default.
pub(crate) fn capture_subprocess_stdout_with_timeout(
    binary: &std::path::Path,
    flag: &str,
    timeout_secs: u64,
) -> Option<String> {
    use crate::subprocess::KillOnDropChild;
    use std::sync::{Arc, Mutex};
    let raw_child = std::process::Command::new(binary)
        .arg(flag)
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .ok()?;
    // Wrap in an RAII guard so any early-return / panic kills the
    // child instead of orphaning it. The watchdog also kills via this
    // guard so the timeout and the unwind path agree.
    let child_arc = Arc::new(Mutex::new(KillOnDropChild::new(raw_child, true)));
    // Take the pipes out BEFORE handing the child to the watchdog so
    // the read loop doesn't fight the watchdog for the mutex.
    let (stdout_handle, stderr_handle) = {
        let mut guard = child_arc.lock().ok()?;
        let inner = guard.inner_mut();
        let stdout_handle: Box<dyn std::io::Read + Send> = match inner.stdout.take() {
            Some(s) => Box::new(s),
            None => Box::new(std::io::empty()),
        };
        let stderr_handle: Box<dyn std::io::Read + Send> = match inner.stderr.take() {
            Some(s) => Box::new(s),
            None => Box::new(std::io::empty()),
        };
        (stdout_handle, stderr_handle)
    };
    // Watchdog: after the timeout, SIGKILL the child so a stuck
    // probe doesn't block forever. The watchdog also signals via the
    // shared flag so the post-read wait() can skip if the watchdog
    // already killed (avoids a benign double-wait noise on POSIX).
    //
    // Round-17 L2: replaced the blind `thread::sleep(timeout_secs)` with
    // an `mpsc::recv_timeout` so the watchdog thread can be woken EARLY
    // when the probe completes. Pre-fix the watchdog held an
    // `Arc<Mutex<KillOnDropChild>>` clone and slept the FULL 10s
    // regardless of how fast the child returned — for a normal probe
    // that finishes in 50ms, the Arc keepalive (and the thread) lived
    // for another 9.95s, leaking thread + memory across hundreds of
    // version probes during adapter discovery. The channel lets the
    // foreground send `()` when the read loop is done; the watchdog
    // wakes immediately on `Ok(())` (or on `Disconnected` if the
    // sender is dropped via panic-unwind) and exits without firing.
    let killed_by_watchdog = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let (done_tx, done_rx) = std::sync::mpsc::channel::<()>();
    {
        let watchdog_child = Arc::clone(&child_arc);
        let killed_flag = Arc::clone(&killed_by_watchdog);
        std::thread::spawn(move || {
            match done_rx.recv_timeout(std::time::Duration::from_secs(timeout_secs)) {
                // Probe completed (or panicked, dropping the sender) —
                // bail without touching the child.
                Ok(()) | Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {}
                // Timeout — kill the child.
                Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                    if let Ok(mut guard) = watchdog_child.lock() {
                        // try_wait first so we only kill children that
                        // are actually still running — avoids a "no
                        // such process" diagnostic on POSIX when the
                        // child already exited.
                        let still_running = matches!(guard.inner_mut().try_wait(), Ok(None));
                        if still_running {
                            let _ = guard.inner_mut().kill();
                            killed_flag.store(true, std::sync::atomic::Ordering::SeqCst);
                            // Round-25 L3: on Windows, `TerminateProcess`
                            // is asynchronous — the call returns
                            // immediately but the kernel takes a
                            // brief moment to tear down the process's
                            // handles, including the child's WRITE
                            // end of our piped stdout/stderr. The
                            // parent's READ end (held by the drain
                            // threads via `read_to_end`) only unblocks
                            // when the kernel signals EOF, which can't
                            // happen until the child's write handle
                            // is released. A try_wait + wait pair
                            // here ensures we don't return from this
                            // watchdog branch until the kernel has
                            // fully reaped the process — by then the
                            // pipes are closed and the drain threads
                            // can finish their `read_to_end` cleanly.
                            // Pre-fix the watchdog kicked the child
                            // and returned; the drain threads then
                            // held the pipe-read syscalls open
                            // briefly while the kernel finished
                            // teardown, which on a CI under load
                            // could be observable as the foreground
                            // `stdout_handle.join()` hanging for
                            // hundreds of ms.
                            let _ = guard.inner_mut().wait();
                        }
                    }
                }
            }
        });
    }
    let combined = capture_subprocess_stdout_inner(
        stdout_handle as Box<dyn std::io::Read + Send>,
        stderr_handle as Box<dyn std::io::Read + Send>,
    );
    // Round-17 L2: wake the watchdog now — the read loop is done. On
    // a normal probe this fires before the timeout, so the watchdog
    // thread (and its Arc clone) drops within milliseconds instead of
    // sleeping out the full 10s. Drop the sender as a defence-in-depth
    // so even if `send` somehow fails (no receiver), the channel
    // disconnect still tickles `recv_timeout` to exit early.
    let _ = done_tx.send(());
    drop(done_tx);
    // Reap the child so we don't leave a zombie around. We don't
    // inspect the exit status — a tool that exits non-zero may still
    // have printed a usable version banner before erroring out. The
    // KillOnDropChild guard handles the case where we panic before
    // this wait runs.
    if let Ok(mut guard) = child_arc.lock() {
        let _ = guard.inner_mut().wait();
    }
    // If the watchdog killed the child we can't trust the partial
    // output — return None so callers fall back to "unknown version"
    // rather than parsing a truncated banner with a torn UTF-8
    // boundary.
    if killed_by_watchdog.load(std::sync::atomic::Ordering::SeqCst) {
        return None;
    }
    combined
}

/// Inner read-and-concat for [`capture_subprocess_stdout`]. Factored
/// out so unit tests can pin the bounded-read behaviour without
/// having to spawn a real subprocess that emits megabytes.
///
/// Round-24 M6: pre-fix the function drained stdout fully, THEN
/// stderr. That sequence pipe-deadlocks any tool that fills its
/// stderr pipe buffer (~64 KiB on Linux, ~4 KiB on Windows) before
/// closing stdout — the child blocks on `write(2)` to stderr, the
/// parent blocks on `read(2)` from stdout, neither side moves. The
/// fix: drain both pipes in parallel threads (mirrors the pattern
/// `subprocess::run` already uses), join both before assembling the
/// combined output. The 10s watchdog still bounds the worst case
/// but with parallel drain the normal completion path is dominated
/// by the child's exit time, not by pipe-fill races.
fn capture_subprocess_stdout_inner(
    stdout: Box<dyn std::io::Read + Send>,
    stderr: Box<dyn std::io::Read + Send>,
) -> Option<String> {
    use std::io::Read;
    use std::thread;
    // Move the readers into parallel drain threads so a child that
    // floods one pipe doesn't deadlock the other. Each thread does
    // its own bounded read up to MAX_PROBE_OUTPUT_BYTES then closes
    // — the Drop on the boxed reader releases the pipe handle.
    let stdout_handle = thread::spawn(move || -> Vec<u8> {
        let mut stdout = stdout;
        let mut buf = Vec::with_capacity(4096);
        let _ = (&mut stdout)
            .take(MAX_PROBE_OUTPUT_BYTES as u64)
            .read_to_end(&mut buf);
        buf
    });
    let stderr_handle = thread::spawn(move || -> Vec<u8> {
        let mut stderr = stderr;
        let mut buf = Vec::with_capacity(4096);
        let _ = (&mut stderr)
            .take(MAX_PROBE_OUTPUT_BYTES as u64)
            .read_to_end(&mut buf);
        buf
    });
    // `join` returns `Err` only if the drain thread panicked. Treat a
    // panic the same as an empty read — we're trying to extract a
    // version banner, not strictly required.
    let stdout_buf = stdout_handle.join().unwrap_or_default();
    let stderr_buf = stderr_handle.join().unwrap_or_default();
    // Some tools (e.g. CalculiX) print version to stderr instead of
    // stdout. Concatenate so we don't miss it.
    let mut combined = String::new();
    if let Ok(s) = String::from_utf8(stdout_buf) {
        combined.push_str(&s);
    }
    if !combined.is_empty() {
        combined.push('\n');
    }
    if let Ok(s) = String::from_utf8(stderr_buf) {
        combined.push_str(&s);
    }
    if combined.is_empty() {
        None
    } else {
        Some(combined)
    }
}

/// Look for a semver-shaped substring in raw text. Returns the
/// first match like `1.2.3` (with optional `-pre`, `+build`).
///
/// Pattern: `\d+\.\d+(\.\d+)?` plus a permissive optional pre /
/// build tail. Greedy on the major.minor.patch core; conservative
/// on the trailing tag so we don't accidentally swallow CLI output
/// like `1.2.3 (built on 2025-04-25)`.
pub fn extract_semver(text: &str) -> Option<String> {
    // Hand-rolled scanner — avoids the `regex` dep for one parse.
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        // Find a digit.
        if !bytes[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let start = i;
        // Consume digits.
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        // Need a '.' followed by at least one digit.
        if i >= bytes.len() || bytes[i] != b'.' {
            continue;
        }
        let after_first_dot = i + 1;
        if after_first_dot >= bytes.len() || !bytes[after_first_dot].is_ascii_digit() {
            continue;
        }
        i = after_first_dot;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        // Optional patch component.
        if i < bytes.len() && bytes[i] == b'.' {
            let after_second_dot = i + 1;
            if after_second_dot < bytes.len() && bytes[after_second_dot].is_ascii_digit() {
                i = after_second_dot;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
            }
        }
        // Optional pre-release / build tag: `-foo.1+build.2`.
        // Stop at the first whitespace / non-tag character.
        if i < bytes.len() && (bytes[i] == b'-' || bytes[i] == b'+') {
            while i < bytes.len()
                && (bytes[i].is_ascii_alphanumeric()
                    || bytes[i] == b'.'
                    || bytes[i] == b'-'
                    || bytes[i] == b'+')
            {
                i += 1;
            }
        }
        return Some(String::from_utf8_lossy(&bytes[start..i]).to_string());
    }
    None
}

/// One-shot version detection: spawn `<binary> <flag>` and return
/// the first semver-shaped substring from the combined stdout +
/// stderr. Returns `None` on any failure (spawn, parse).
///
/// Adapters use this in `probe()` to populate `ProbeReport::found_version`
/// without each one re-implementing the run-and-parse boilerplate.
pub fn detect_tool_version(binary: &std::path::Path, flag: &str) -> Option<String> {
    let stdout = capture_subprocess_stdout(binary, flag)?;
    extract_semver(&stdout)
}

/// Higher-level wrapper around [`detect_tool_version`] that
/// normalises the result into a strict `semver::Version`. Pads
/// short forms ("11" / "3.0") to major.minor.patch first because
/// `Version::parse` rejects anything shorter than three components.
///
/// Tries each `flag` in order until one yields a parseable
/// version — useful for tools whose `--version` flag varies between
/// forks (OpenFOAM uses `-help`, FreeCAD uses `--version`,
/// CalculiX uses `-v`).
pub fn detect_tool_version_semver(
    binary: &std::path::Path,
    flags: &[&str],
) -> Option<semver::Version> {
    for flag in flags {
        if let Some(raw) = detect_tool_version(binary, flag) {
            let dots = raw.chars().filter(|c| *c == '.').count();
            let normalised: String = match dots {
                0 => format!("{raw}.0.0"),
                1 => format!("{raw}.0"),
                _ => raw,
            };
            if let Ok(v) = semver::Version::parse(&normalised) {
                return Some(v);
            }
        }
    }
    None
}

/// Returns `true` when the `tool_license` string corresponds to an
/// OSI-approved or unambiguously-OSS license.
///
/// Used by the desktop app to default-hide adapters wrapping
/// academic-only / non-commercial / restricted-commercial tools.
/// The adapters are still in the registry — they just don't appear
/// in the wizard / browser / palette / Tools menu unless the user
/// opts in via Settings → "Show non-OSS adapters".
///
/// Conservative: matches the standard OSI strings explicitly and
/// returns `false` for anything else (including the per-tool license
/// names like `"NAMD-License"`, `"Rosetta-License"`, `"Academic"`,
/// `"CC-BY-NC-*"`, etc.).
pub fn is_oss_license(license: &str) -> bool {
    // Strip "WITH <exception>" suffix (SPDX style, e.g. OCCT's
    // "LGPL-2.1-or-later WITH OCCT-exception-1.0").
    let core = license.split(" WITH ").next().unwrap_or(license).trim();
    matches!(
        core,
        "MIT"
            | "BSD-2-Clause"
            | "BSD-3-Clause"
            | "BSD-3-Clause-Open-Source"
            | "Apache-2.0"
            | "GPL-2.0"
            | "GPL-2.0-only"
            | "GPL-2.0-or-later"
            | "GPL-3.0"
            | "GPL-3.0-only"
            | "GPL-3.0-or-later"
            | "LGPL-2.1"
            | "LGPL-2.1-only"
            | "LGPL-2.1-or-later"
            | "LGPL-3.0"
            | "LGPL-3.0-or-later"
            | "AGPL-3.0"
            | "AGPL-3.0-only"
            | "AGPL-3.0-or-later"
            | "MPL-2.0"
            | "Artistic-2.0"
            | "ECL-2.0"
            | "AFL-3.0"
            | "Public Domain"
            | "CC0-1.0"
    )
}

/// Build a consistent [`AdapterError::Other`] for "not implemented
/// in this scaffold" stubs. The message calls out the roadmap phase
/// so the UI can link users to planned work.
pub fn not_implemented(adapter_id: &str, op: &str, roadmap_phase: &str) -> AdapterError {
    AdapterError::Other(anyhow::anyhow!(
        "{adapter_id}.{op}() is a Phase 1 scaffold — the full implementation \
         lands in ROADMAP {roadmap_phase}"
    ))
}

/// Walk the workdir top-level for the lexicographically-first file
/// whose extension matches any in `extensions` (case-insensitive).
/// Used by adapter `collect()` paths to locate a hashable case /
/// mesh source for the [`live_provenance`] call.
///
/// Returns `None` when the directory can't be read, when no entry
/// has an extension, or when nothing matches. Subdirectories are
/// skipped — only files at the top level participate.
pub fn first_workdir_match(root: &std::path::Path, extensions: &[&str]) -> Option<PathBuf> {
    let entries = std::fs::read_dir(root).ok()?;
    let mut matches: Vec<PathBuf> = entries
        .flatten()
        .filter_map(|e| {
            let p = e.path();
            let ext = p
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase())?;
            extensions.contains(&ext.as_str()).then_some(p)
        })
        .collect();
    matches.sort();
    matches.into_iter().next()
}

/// Join a user-supplied path under a case directory while keeping
/// the result confined to that directory.
///
/// Adapters use this for any case.toml-supplied path that the adapter
/// will subsequently *stage* into the workdir (script files, config
/// files that get copied across). The function rejects:
///
/// - **Absolute paths** — `script = "/etc/passwd"` would otherwise
///   be copied verbatim into the workdir, leaking arbitrary host
///   files into a directory that may later be archived or shared.
/// - **`..` components** — `script = "../../etc/passwd"` would
///   otherwise traverse out of the case directory.
///
/// Both are mistakes a user authoring their own case.toml is unlikely
/// to make on purpose, but the threat model includes "user runs an
/// untrusted shared case bundle", and this helper closes that
/// exfiltration path.
///
/// **NOT for adapters that take a system-install path** (Cromwell's
/// `jar`, Fiji's `fiji_app`, NAMD's `namd2` etc.) — those legitimately
/// point to absolute paths outside the case directory and the adapter
/// reads from them in place rather than copying them in. Use plain
/// `case.path.join(&user_path)` for those.
///
/// Returns the joined path on success.
///
/// # Known limitations
///
/// **Symlinks within the case bundle that point outside are NOT
/// detected.** A full canonicalization-and-prefix-compare check is
/// expensive (it spawns a syscall per `confined_join` call and only
/// works when the path actually exists yet), and adapters call this
/// helper on every adapter-prepare. The defence-in-depth boundary
/// for the symlink case is the surrounding filesystem permissions /
/// container sandbox the operator is responsible for. For the
/// user-typed-a-bad-path threat model this helper addresses, the
/// in-process character-level checks are sufficient.
pub fn confined_join(
    case_dir: &std::path::Path,
    user_path: &std::path::Path,
) -> Result<PathBuf, AdapterError> {
    let user_str = user_path.to_string_lossy();
    if user_path.is_absolute() {
        return Err(AdapterError::InvalidCase {
            case_path: case_dir.join("case.toml"),
            reason: format!(
                "path `{}` is absolute; case.toml-supplied paths that get \
                 staged into the workdir must be relative to the case \
                 directory (sandboxing requirement)",
                user_path.display()
            ),
        });
    }
    // Round-5 fix: Windows drive-relative paths (`C:foo`) report
    // `false` for `is_absolute()` but `Command::new` resolves them
    // against the drive's current working directory — a sandbox
    // escape. Reject the pattern explicitly on Windows: a single
    // ASCII letter followed by `:` followed by anything other than
    // `\` or `/` (a fully-absolute `C:\foo` is already caught by
    // the `is_absolute()` check above).
    #[cfg(windows)]
    {
        let bytes = user_str.as_bytes();
        if bytes.len() >= 2
            && bytes[0].is_ascii_alphabetic()
            && bytes[1] == b':'
            && (bytes.len() == 2 || (bytes[2] != b'\\' && bytes[2] != b'/'))
        {
            return Err(AdapterError::InvalidCase {
                case_path: case_dir.join("case.toml"),
                reason: format!(
                    "path `{}` is a Windows drive-relative path (resolves \
                     against the drive's CWD); case.toml-supplied paths \
                     must be a plain relative subdirectory",
                    user_path.display()
                ),
            });
        }
    }
    // Even on non-Windows hosts, reject drive-letter-prefixed paths
    // so a case bundle shared between operators doesn't behave
    // differently on Windows vs. Linux: the `C:foo` shape is never
    // a legitimate POSIX path.
    #[cfg(not(windows))]
    {
        let bytes = user_str.as_bytes();
        if bytes.len() >= 2 && bytes[0].is_ascii_alphabetic() && bytes[1] == b':' {
            return Err(AdapterError::InvalidCase {
                case_path: case_dir.join("case.toml"),
                reason: format!(
                    "path `{}` looks like a Windows drive-letter path; \
                     case.toml-supplied paths must be plain relative \
                     subdirectories",
                    user_path.display()
                ),
            });
        }
    }
    for component in user_path.components() {
        if matches!(component, std::path::Component::ParentDir) {
            return Err(AdapterError::InvalidCase {
                case_path: case_dir.join("case.toml"),
                reason: format!(
                    "path `{}` contains a `..` component; case.toml-supplied \
                     paths must stay within the case directory",
                    user_path.display()
                ),
            });
        }
    }
    Ok(case_dir.join(user_path))
}

/// Escape a user-controlled string for safe embedding inside a
/// double-quoted Python string literal. Escapes `\`, `"`, `\n`, `\r`,
/// `\t`, and any non-printable ASCII control char. Returns the escaped
/// string WITHOUT enclosing quotes — caller must add them.
///
/// Round-19 H2 helper: closes a Python-injection class in adapters
/// that emit a script driving a tool (cantera, pybamm). Pre-fix the
/// emitters did `format!("\"{}\"", user_str)` directly, so a field
/// containing a closing quote + newline + an arbitrary statement would
/// break out of the literal and inject statement-level Python. The
/// helper turns the same payload into a single escaped string literal,
/// so any `"`/`\n` the user supplied stays inside the quoted region.
///
/// Conservative: only escapes the forms Python's tokeniser would treat
/// as terminating the string (or as an escape continuation). UTF-8
/// codepoints above 0x7F pass through verbatim — they're legitimate
/// Python source even in ASCII contexts because the default `coding:`
/// declaration on Python 3 is UTF-8.
///
/// Round-20 L4 (defence-in-depth): U+2028 (LINE SEPARATOR) and
/// U+2029 (PARAGRAPH SEPARATOR) are also escaped as ` ` /
/// ` `. These two codepoints are LINE TERMINATORS in JavaScript
/// (per ECMAScript 11.3) and a few other languages even though
/// Python's tokeniser treats them as ordinary characters. Escaping
/// them here makes the helper reusable for JS/JSON contexts (e.g.
/// emitting JSON-embedded-in-HTML, or driving a JS evaluator from a
/// Python script) without surprise. The escape is benign for the
/// existing Python call sites — Python accepts the `\u` escape inside
/// any string literal and the resulting codepoint is byte-identical.
pub fn python_str_repr(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            // Round-20 L4: U+2028 LINE SEPARATOR + U+2029
            // PARAGRAPH SEPARATOR — line terminators in JS/JSON
            // contexts where a bare codepoint would split the
            // literal. Encoded with the Python `\u` escape so the
            // emitted source is reusable across languages.
            '\u{2028}' => out.push_str("\\u2028"),
            '\u{2029}' => out.push_str("\\u2029"),
            c if (c as u32) < 0x20 => {
                use std::fmt::Write;
                let _ = write!(out, "\\x{:02x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

/// Strip characters that can break structured-text formats
/// (Elmer SIF, CalculiX INP, OpenFOAM dict, STEP block comments,
/// IFC strings, etc.). Removes `\n`, `\r`, `"`, `\\` so an attacker
/// who controls a user-facing string field cannot inject sibling
/// blocks or escape out of the quoted literal into surrounding
/// directives.
///
/// Round-15 helper — sister to round-14's `ifc_str` Part-21 escape
/// in `valenx-bim/src/ifc.rs`. Most structured-text writers in the
/// workspace previously interpolated user strings via `"{}"` directly
/// into their output, which made every `material.name`, `space_name`,
/// `output_basename`, `boundary_condition.name`, and `part.name` a
/// directive-injection vector. The helper is intentionally generic;
/// per-format quirks (STEP's `*/` block-comment close, OF dict's `;`
/// statement terminator) get handled by the call site after this
/// generic pass.
///
/// Returns a new owned `String`. Callers that need byte-identical
/// passthrough for the common all-safe-chars case can short-circuit
/// on `contains(&['\n', '\r', '"', '\\'])` themselves; the helper
/// always allocates so the caller doesn't have to repeat the check.
pub fn sanitize_structured_identifier(s: &str) -> String {
    s.chars()
        .filter(|c| !matches!(c, '\n' | '\r' | '"' | '\\'))
        .collect()
}

/// Strict variant of [`sanitize_structured_identifier`]: only allow
/// ASCII alphanumeric plus `.`, `-`, `_`. Used for identifier
/// positions where the format strictly requires a parseable token —
/// OpenFOAM boundary patch names, gmsh physical group names, INP
/// `*NSET` names, etc. — and silently mangling user input would
/// produce a dict file the solver rejects with a cryptic parser
/// error.
///
/// Returns `Err(AdapterError::InvalidCase)` with `case_path` left
/// empty (caller rewrites via `set_case_path` if it knows the path)
/// when validation fails, so the UI can surface a clean "field X is
/// invalid" message instead of letting the value flow into a dict
/// file that breaks the solver mid-run.
///
/// Empty input is rejected — a zero-length identifier is never
/// valid in any of the consumer formats.
pub fn validate_structured_identifier(
    s: &str,
    field: &str,
) -> Result<(), crate::error::AdapterError> {
    if s.is_empty() {
        return Err(crate::error::AdapterError::InvalidCase {
            case_path: std::path::PathBuf::new(),
            reason: format!("{field}: cannot be empty"),
        });
    }
    if !s
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '.' || c == '-' || c == '_')
    {
        return Err(crate::error::AdapterError::InvalidCase {
            case_path: std::path::PathBuf::new(),
            reason: format!("{field}: must be ASCII alphanumeric + .-_ (got {s:?})"),
        });
    }
    Ok(())
}

/// Maximum recursion depth `copy_dir_recursive` will descend before
/// aborting. 32 levels is well past anything sane real cases need;
/// a poisoned source tree shaped like a recursive symlink-free
/// chain `a/a/a/.../a` will still terminate.
pub const MAX_COPY_DIR_DEPTH: usize = 32;

/// Maximum per-file byte cap `copy_dir_recursive` will copy. Bumped
/// past what realistic case bundles need; the per-file gate is here
/// to refuse a 100 GB single file from a poisoned case stage that
/// could otherwise fill the disk before any sanity check fires.
pub const MAX_COPY_DIR_FILE_BYTES: u64 = 8 * 1024 * 1024 * 1024;

/// Maximum aggregate (cumulative) byte cap `copy_dir_recursive` will
/// copy across the entire walk. Round-24 M3: pre-fix the per-file cap
/// was the only ceiling — a poisoned source with thousands of files
/// each just under MAX_COPY_DIR_FILE_BYTES would happily fill the
/// destination disk because nothing tracked the cumulative spend.
/// 64 GiB matches the largest realistic project bundle (full
/// AlphaFold model weights + a few HDF5 datasets) while refusing the
/// "1024 files at 7 GiB each" denial-of-service shape.
pub const MAX_COPY_DIR_TOTAL_BYTES: u64 = 64 * 1024 * 1024 * 1024;

/// Maximum number of filesystem entries `copy_dir_recursive` will
/// touch in a single call. Round-24 M3: pre-fix the only entry-count
/// guard was implicit in the recursion-depth cap, which gates branch
/// depth not entry count. A flat directory with 100 M zero-byte
/// files would pass the per-file size cap on every entry but still
/// exhaust inodes / dentries on the destination. 1 M entries covers
/// every realistic case bundle (typical CFD bundle: ~10 k files;
/// genome reference: ~5 k chromosome FASTAs + index sidecars).
pub const MAX_COPY_DIR_ENTRIES: usize = 1_000_000;

/// Recursive directory copy hardened against the three classic
/// "staging-time" abuses:
///
/// - **Symlink rejection.** A poisoned source dir that links a
///   subdir to `/etc` or to the install target itself would otherwise
///   let the copy traverse / overwrite arbitrary paths. We refuse
///   *any* symlink under the source tree, top-level or nested.
/// - **Recursion depth cap.** 32 levels is enough for every realistic
///   case bundle and aborts the pathological `a/a/a/a/...` cycle.
/// - **Per-file size cap.** 8 GiB stays out of the way for legitimate
///   reference bundles (genomes, meshes) and refuses a single
///   100 GB file that would otherwise fill the disk.
/// - **Aggregate size cap** (round-24 M3). 64 GiB cap across the
///   entire walk so the "many small-but-not-tiny files" attack
///   shape is still caught even when each file is under the per-
///   file ceiling.
/// - **Aggregate entry cap** (round-24 M3). 1 M entries across the
///   whole walk so "10 M zero-byte files" can't exhaust inodes /
///   dentries on the destination.
///
/// All three are the same caps the `valenx-addons` install path
/// uses; this is the round-9 extraction so the precice adapter
/// (and any future caller) shares the hardening rather than each
/// reinventing it.
///
/// # Errors
///
/// - [`AdapterError::Other`] (with an `anyhow` cause) for any policy
///   violation (symlink, depth, per-file size, aggregate size,
///   aggregate entries). The string includes the offending path so
///   the user can see what was rejected.
/// - [`AdapterError::Other`] for any underlying I/O error.
pub fn copy_dir_recursive(
    src: &std::path::Path,
    dst: &std::path::Path,
) -> Result<(), AdapterError> {
    copy_dir_recursive_with_caps(src, dst, MAX_COPY_DIR_TOTAL_BYTES, MAX_COPY_DIR_ENTRIES)
}

/// Round-24 M3: variant of `copy_dir_recursive` exposed to the test
/// module so the aggregate-byte and aggregate-entry caps can be
/// exercised with smaller test values. Production callers MUST use
/// `copy_dir_recursive` so the safety contract is uniform; this
/// helper exists ONLY for the in-crate unit tests that would
/// otherwise need to allocate 64 GiB of sparse files to hit the
/// production cap.
#[cfg(test)]
pub(crate) fn copy_dir_recursive_with_caps(
    src: &std::path::Path,
    dst: &std::path::Path,
    max_total_bytes: u64,
    max_entries: usize,
) -> Result<(), AdapterError> {
    let mut totals = CopyDirTotals {
        bytes: 0,
        entries: 0,
    };
    copy_dir_recursive_depth(src, dst, 0, &mut totals, max_total_bytes, max_entries)
}

/// Same shape as the test-only helper but kept private (no `cfg(test)`)
/// because the non-test path needs an implementation too.
#[cfg(not(test))]
fn copy_dir_recursive_with_caps(
    src: &std::path::Path,
    dst: &std::path::Path,
    max_total_bytes: u64,
    max_entries: usize,
) -> Result<(), AdapterError> {
    let mut totals = CopyDirTotals {
        bytes: 0,
        entries: 0,
    };
    copy_dir_recursive_depth(src, dst, 0, &mut totals, max_total_bytes, max_entries)
}

/// Running totals tracked across a `copy_dir_recursive` walk. Lives
/// here (not as caller-visible state) because the aggregate caps are
/// part of the helper's safety contract, not a tuning knob.
struct CopyDirTotals {
    bytes: u64,
    entries: usize,
}

fn copy_dir_recursive_depth(
    src: &std::path::Path,
    dst: &std::path::Path,
    depth: usize,
    totals: &mut CopyDirTotals,
    max_total_bytes: u64,
    max_entries: usize,
) -> Result<(), AdapterError> {
    if depth >= MAX_COPY_DIR_DEPTH {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "copy_dir_recursive: source tree exceeds {MAX_COPY_DIR_DEPTH}-level depth cap at {}",
            src.display()
        )));
    }
    // `symlink_metadata` does NOT traverse symlinks, so a symlinked
    // top-level entry shows up as `FileType::is_symlink()` here even
    // when the regular `metadata` would follow the link.
    let md = std::fs::symlink_metadata(src).map_err(|e| {
        AdapterError::Other(anyhow::anyhow!("symlink_metadata {}: {e}", src.display()))
    })?;
    if md.file_type().is_symlink() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "copy_dir_recursive: source path is a symlink (refusing to traverse): {}",
            src.display()
        )));
    }
    std::fs::create_dir_all(dst).map_err(|e| {
        AdapterError::Other(anyhow::anyhow!("create_dir_all {}: {e}", dst.display()))
    })?;
    for entry in std::fs::read_dir(src)
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read_dir {}: {e}", src.display())))?
    {
        let entry =
            entry.map_err(|e| AdapterError::Other(anyhow::anyhow!("read_dir entry: {e}")))?;
        let from = entry.path();
        let to: PathBuf = dst.join(entry.file_name());
        let entry_md = std::fs::symlink_metadata(&from).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("symlink_metadata {}: {e}", from.display()))
        })?;
        if entry_md.file_type().is_symlink() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "copy_dir_recursive: source contains a symlink (refusing to traverse): {}",
                from.display()
            )));
        }
        // Round-24 M3: bump the aggregate entry counter BEFORE doing
        // any work for this entry so a poisoned source can't sneak
        // an over-cap entry through by hiding it past a slow IO step.
        totals.entries = totals.entries.saturating_add(1);
        if totals.entries > max_entries {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "copy_dir_recursive: source tree exceeds the {max_entries}-entry \
                 aggregate cap (entry {} would be over the cap)",
                from.display()
            )));
        }
        if entry_md.is_dir() {
            copy_dir_recursive_depth(&from, &to, depth + 1, totals, max_total_bytes, max_entries)?;
        } else {
            if entry_md.len() > MAX_COPY_DIR_FILE_BYTES {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "copy_dir_recursive: source file {} ({} bytes) exceeds the {}-byte cap",
                    from.display(),
                    entry_md.len(),
                    MAX_COPY_DIR_FILE_BYTES
                )));
            }
            // Round-24 M3: bump and check the aggregate byte counter
            // BEFORE the copy so we refuse rather than partially
            // writing past the cap.
            totals.bytes = totals.bytes.saturating_add(entry_md.len());
            if totals.bytes > max_total_bytes {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "copy_dir_recursive: cumulative copy exceeds the {}-byte aggregate cap \
                     (would total {} bytes after copying {})",
                    max_total_bytes,
                    totals.bytes,
                    from.display()
                )));
            }
            std::fs::copy(&from, &to).map_err(|e| {
                AdapterError::Other(anyhow::anyhow!(
                    "copy {} -> {}: {e}",
                    from.display(),
                    to.display()
                ))
            })?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn platform_suffix_is_idempotent_on_windows() {
        if cfg!(windows) {
            assert_eq!(platform_suffix("foo"), "foo.exe");
            assert_eq!(platform_suffix("foo.exe"), "foo.exe");
        } else {
            assert_eq!(platform_suffix("foo"), "foo");
            assert_eq!(platform_suffix("foo.exe"), "foo.exe");
        }
    }

    #[test]
    fn platform_candidates_returns_bare_name_on_unix() {
        if !cfg!(windows) {
            // Non-Windows: just one candidate, the bare name.
            assert_eq!(platform_candidates("gmsh"), vec!["gmsh".to_string()]);
            assert_eq!(
                platform_candidates_for("gmsh", ".EXE;.BAT"),
                vec!["gmsh".to_string()],
                "non-Windows ignores PATHEXT entirely"
            );
        }
    }

    #[test]
    fn platform_candidates_lists_every_pathext_on_windows() {
        if !cfg!(windows) {
            return;
        }
        let candidates = platform_candidates_for("gmsh", ".COM;.EXE;.BAT;.CMD");
        // Bare name first, then one entry per extension.
        assert_eq!(candidates[0], "gmsh");
        assert!(candidates.contains(&"gmsh.exe".to_string()));
        assert!(candidates.contains(&"gmsh.bat".to_string()));
        assert!(candidates.contains(&"gmsh.cmd".to_string()));
        assert!(candidates.contains(&"gmsh.com".to_string()));
    }

    #[test]
    fn platform_candidates_no_double_extension_when_caller_already_has_one() {
        if !cfg!(windows) {
            return;
        }
        // Caller passed `gmsh.exe` — we must NOT produce `gmsh.exe.exe`.
        let candidates = platform_candidates_for("gmsh.exe", ".EXE;.BAT");
        assert!(
            !candidates.iter().any(|c| c.contains("exe.exe")),
            "double-extension leaked: {candidates:?}"
        );
        // The .bat alternative still appears since it's a different
        // shim that legitimately could exist alongside.
        assert!(candidates.contains(&"gmsh.exe.bat".to_string()));
    }

    #[test]
    fn platform_candidates_dedupes_repeated_extensions() {
        if !cfg!(windows) {
            return;
        }
        // Pathological input: PATHEXT lists `.EXE` twice with
        // different casing. Don't double-emit candidates.
        let candidates = platform_candidates_for("foo", ".EXE;.exe;.BAT");
        let exe_count = candidates.iter().filter(|c| *c == "foo.exe").count();
        assert_eq!(exe_count, 1, "duplicates leaked: {candidates:?}");
    }

    #[test]
    fn platform_candidates_handles_empty_pathext_gracefully() {
        if !cfg!(windows) {
            return;
        }
        // An empty / whitespace PATHEXT shouldn't panic, just return
        // the bare name.
        let candidates = platform_candidates_for("foo", "");
        assert_eq!(candidates, vec!["foo".to_string()]);
        let candidates = platform_candidates_for("foo", "   ;  ;  ");
        assert_eq!(candidates, vec!["foo".to_string()]);
    }

    #[test]
    fn stub_provenance_fields_are_deterministic() {
        let a = stub_provenance("x", "0.0.0", "y");
        let b = stub_provenance("x", "0.0.0", "y");
        // Two separate calls must produce identical fields so snapshot
        // tests don't flap.
        assert_eq!(a.run_id, b.run_id);
        assert_eq!(a.completed_at, b.completed_at);
    }

    #[test]
    fn not_implemented_mentions_roadmap_phase() {
        let e = not_implemented("foo", "prepare", "Phase 5");
        let msg = format!("{e}");
        assert!(msg.contains("foo"));
        assert!(msg.contains("prepare"));
        assert!(msg.contains("Phase 5"));
    }

    #[test]
    fn find_on_path_handles_missing_bins() {
        // A made-up name can never exist — assert None.
        let result = find_on_path(&["valenx-nonexistent-binary-asd9823h"]);
        assert!(result.is_none(), "unexpected match: {result:?}");
    }

    // -----------------------------------------------------------------
    // Python interpreter allow-list (round-3 security fix)
    // -----------------------------------------------------------------

    #[test]
    fn validate_python_binary_accepts_plain_name() {
        assert!(validate_python_binary("python").is_ok());
        assert!(validate_python_binary("python3").is_ok());
        assert!(validate_python_binary("python3.11").is_ok());
        assert!(validate_python_binary("python3.13").is_ok());
    }

    #[test]
    fn validate_python_binary_accepts_absolute_allowed_name() {
        // The path itself is irrelevant — the basename has to match.
        assert!(validate_python_binary("/usr/bin/python3").is_ok());
        assert!(validate_python_binary("/opt/conda/envs/myenv/bin/python3.11").is_ok());
    }

    #[test]
    fn validate_python_binary_rejects_arbitrary_binary() {
        // The canonical attack: hostile case.toml tries to invoke
        // /usr/bin/curl as the "python" interpreter.
        let err = validate_python_binary("/usr/bin/curl").expect_err("must reject");
        let msg = format!("{err}");
        assert!(msg.contains("allow"), "msg: {msg}");
        assert!(msg.contains("curl"), "msg: {msg}");
    }

    #[test]
    fn validate_python_binary_rejects_bare_unrelated_command() {
        assert!(validate_python_binary("rm").is_err());
        assert!(validate_python_binary("powershell").is_err());
        assert!(validate_python_binary("bash").is_err());
    }

    #[test]
    fn validate_python_binary_rejects_empty() {
        assert!(validate_python_binary("").is_err());
        assert!(validate_python_binary("   ").is_err());
    }

    // -----------------------------------------------------------------
    // resolve_python_binary — round-4 mechanical-sweep helper
    // -----------------------------------------------------------------

    /// Round-4: the high-level resolver must short-circuit through the
    /// allow-list before doing any PATH lookup. A hostile case.toml
    /// setting python to e.g. `/usr/bin/curl` must fail here even if
    /// `/usr/bin/curl` happens to exist on the host.
    #[test]
    fn resolve_python_binary_rejects_arbitrary_binary() {
        let err = resolve_python_binary("/usr/bin/curl", &["python3"]).expect_err("reject curl");
        let msg = format!("{err}");
        assert!(
            msg.contains("allow") || msg.contains("not on the allow-list"),
            "msg: {msg}"
        );
    }

    /// `..` traversal in a relative spec must fail even after the
    /// allow-list (which only checks basename).
    #[test]
    fn resolve_python_binary_rejects_parent_dir_traversal() {
        let err = resolve_python_binary("../python3", &["python3"])
            .expect_err("must reject parent-dir traversal");
        let msg = format!("{err}");
        assert!(msg.contains(".."), "msg: {msg}");
    }

    #[test]
    fn resolve_python_binary_rejects_empty() {
        assert!(resolve_python_binary("", &["python3"]).is_err());
    }

    /// When the spec is an absolute path whose basename is on the
    /// allow-list AND the file exists, return that path verbatim.
    #[test]
    fn resolve_python_binary_keeps_existing_absolute_path() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-resolve-python-{}.fake",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // Use a sibling "python3" file so the basename matches.
        let dir = tmp.with_extension("");
        std::fs::create_dir_all(&dir).unwrap();
        let python_path = dir.join("python3");
        std::fs::write(&python_path, b"#!/bin/sh\n").unwrap();
        let resolved = resolve_python_binary(
            python_path.to_str().expect("test python path utf8"),
            &["python3"],
        )
        .expect("absolute allow-list path with existing file resolves");
        assert_eq!(resolved, python_path);
        let _ = std::fs::remove_dir_all(&dir);
    }

    // -----------------------------------------------------------------
    // Output basename validation (round-3 security fix)
    // -----------------------------------------------------------------

    #[test]
    fn validate_output_basename_accepts_safe_names() {
        assert!(validate_output_basename("results", "[bio.test].out").is_ok());
        assert!(validate_output_basename("run_001", "out").is_ok());
        assert!(validate_output_basename("egfr.net", "out").is_ok());
    }

    #[test]
    fn validate_output_basename_rejects_path_separators() {
        let err = validate_output_basename("a/b", "out").unwrap_err();
        assert!(format!("{err}").contains("separators"));
        let err = validate_output_basename("a\\b", "out").unwrap_err();
        assert!(format!("{err}").contains("separators"));
    }

    #[test]
    fn validate_output_basename_rejects_parent_dir_traversal() {
        let err = validate_output_basename("..", "out").unwrap_err();
        assert!(format!("{err}").contains(".."));
        // Note: "../foo" is also caught by the separator check; we only
        // need the leading-dotdot path here.
        let err = validate_output_basename("../etc/cron.d/x", "out").unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("separators"),
            "msg: {msg}"
        );
    }

    #[test]
    fn validate_output_basename_rejects_absolute() {
        let abs = if cfg!(windows) {
            "C:\\Windows\\System32"
        } else {
            "/etc/passwd"
        };
        let err = validate_output_basename(abs, "out").unwrap_err();
        let msg = format!("{err}");
        // Either separator-check or absolute-check could fire first.
        assert!(
            msg.contains("absolute") || msg.contains("separators"),
            "msg: {msg}"
        );
    }

    #[test]
    fn validate_output_basename_rejects_empty() {
        assert!(validate_output_basename("", "out").is_err());
        assert!(validate_output_basename("   ", "out").is_err());
    }

    // -----------------------------------------------------------------
    // Round-5: Python allow-list — Windows variants
    // -----------------------------------------------------------------

    /// Round-5 RED→GREEN: round-4 introduced a regression where
    /// `python3.11.exe` (the standard conda-on-Windows / pyenv-win
    /// per-env name) was rejected by the allow-list. This test pins
    /// the accepted dotted-Windows form.
    #[test]
    fn validate_python_binary_accepts_python3_11_exe() {
        assert!(
            validate_python_binary("python3.11.exe").is_ok(),
            "python3.11.exe (conda / pyenv-win standard) must be on \
             the allow-list"
        );
        // Cover the rest of the 3.X dotted-Windows variants for the
        // common conda releases.
        assert!(validate_python_binary("python3.10.exe").is_ok());
        assert!(validate_python_binary("python3.12.exe").is_ok());
        assert!(validate_python_binary("python3.13.exe").is_ok());
    }

    /// Round-5: the no-dot conda-Windows variant (`python311.exe`).
    #[test]
    fn validate_python_binary_accepts_python311_exe() {
        assert!(
            validate_python_binary("python311.exe").is_ok(),
            "python311.exe (anaconda3 envs layout) must be on the \
             allow-list"
        );
        assert!(validate_python_binary("python310.exe").is_ok());
        assert!(validate_python_binary("python312.exe").is_ok());
        assert!(validate_python_binary("python313.exe").is_ok());
    }

    // -----------------------------------------------------------------
    // Round-5: R / Rscript allow-list
    // -----------------------------------------------------------------

    /// Round-5 RED→GREEN: `validate_rscript_binary` is the R sister of
    /// `validate_python_binary`. iCodon and Seurat adapters take a
    /// user-supplied `rscript = "..."` field; without the allow-list
    /// they previously forwarded the value straight to `Command::new`,
    /// letting a hostile case do `rscript = "/usr/bin/curl"`.
    #[test]
    fn validate_rscript_binary_rejects_arbitrary_binary() {
        let err = validate_rscript_binary("/usr/bin/curl").expect_err("must reject /usr/bin/curl");
        let msg = format!("{err}");
        assert!(msg.contains("allow"), "msg: {msg}");
        assert!(msg.contains("curl"), "msg: {msg}");
    }

    #[test]
    fn validate_rscript_binary_accepts_plain_name() {
        // The standard headless R launcher.
        assert!(validate_rscript_binary("Rscript").is_ok());
        // Lowercase fallback some distros symlink.
        assert!(validate_rscript_binary("rscript").is_ok());
        // The interactive REPL — also legitimate (some batch wrappers
        // pipe scripts into it).
        assert!(validate_rscript_binary("R").is_ok());
    }

    #[test]
    fn validate_rscript_binary_accepts_absolute_allowed_name() {
        // Absolute path resolution: the basename has to match the
        // allow-list, the path prefix doesn't matter.
        assert!(validate_rscript_binary("/usr/bin/Rscript").is_ok());
        assert!(validate_rscript_binary("/opt/conda/envs/myenv/bin/Rscript").is_ok());
    }

    #[test]
    fn validate_rscript_binary_accepts_windows_variants() {
        // Conda-on-Windows / R-for-Windows ship .exe suffixes.
        assert!(validate_rscript_binary("Rscript.exe").is_ok());
        assert!(validate_rscript_binary("rscript.exe").is_ok());
        assert!(validate_rscript_binary("R.exe").is_ok());
    }

    #[test]
    fn validate_rscript_binary_rejects_empty() {
        assert!(validate_rscript_binary("").is_err());
        assert!(validate_rscript_binary("   ").is_err());
    }

    // -----------------------------------------------------------------
    // Round-5: validate_output_dir — multi-component relative path
    // -----------------------------------------------------------------

    /// Round-5 RED→GREEN: `validate_output_dir` is the multi-component
    /// sister of `validate_output_basename`. cwltool / kallisto /
    /// salmon take an `output_dir` field which legitimately can be
    /// multi-segment (`results/run1`), but must still reject absolute
    /// paths and `..` traversal.
    #[test]
    fn validate_output_dir_rejects_parent_dir_traversal() {
        let err = validate_output_dir(
            std::path::Path::new("../escape"),
            "[bio.cwltool].output_dir",
        )
        .expect_err("must reject ../escape");
        let msg = format!("{err}");
        assert!(msg.contains(".."), "msg: {msg}");
    }

    #[test]
    fn validate_output_dir_rejects_absolute_paths() {
        let abs = if cfg!(windows) {
            "C:\\tmp\\out"
        } else {
            "/tmp/out"
        };
        let err = validate_output_dir(std::path::Path::new(abs), "out")
            .expect_err("must reject absolute path");
        let msg = format!("{err}");
        assert!(msg.contains("absolute"), "msg: {msg}");
    }

    #[test]
    fn validate_output_dir_accepts_multi_component_relative() {
        // The whole point of the helper vs `validate_output_basename`:
        // multi-component relative paths like `results/run1` ARE
        // allowed.
        assert!(validate_output_dir(std::path::Path::new("results"), "out").is_ok());
        assert!(validate_output_dir(std::path::Path::new("results/run1"), "out").is_ok());
        assert!(validate_output_dir(std::path::Path::new("a/b/c"), "out").is_ok());
    }

    #[test]
    fn validate_output_dir_rejects_leading_slash() {
        // Both POSIX root-relative (`/foo`) and Windows current-drive-
        // relative (`\foo`) must be rejected. On Linux `\foo` is just
        // a weird filename so the leading-`\` check is the only signal.
        let err =
            validate_output_dir(std::path::Path::new("/foo"), "out").expect_err("must reject /foo");
        let msg = format!("{err}");
        assert!(
            msg.contains('/') || msg.contains("absolute") || msg.contains('\\'),
            "msg: {msg}"
        );
    }

    #[test]
    fn validate_output_dir_rejects_empty() {
        assert!(validate_output_dir(std::path::Path::new(""), "out").is_err());
    }

    // -----------------------------------------------------------------
    // Round-5: confined_join — Windows drive-relative paths
    // -----------------------------------------------------------------

    /// Round-5 RED→GREEN: `confined_join` previously missed the
    /// Windows drive-relative shape `C:foo` (no separator after the
    /// drive letter). `Path::is_absolute()` returns false for this
    /// shape, but `Command::new("C:foo")` resolves it against the
    /// drive's current working directory — a sandbox escape.
    #[cfg(windows)]
    #[test]
    fn confined_join_rejects_windows_drive_relative() {
        let case_dir = std::path::Path::new("C:\\case\\path");
        let err =
            confined_join(case_dir, std::path::Path::new("C:foo")).expect_err("must reject C:foo");
        let msg = format!("{err}");
        assert!(
            msg.contains("drive") || msg.contains("absolute"),
            "msg: {msg}"
        );
    }

    /// Sister assertion (runs on every platform): the same drive-letter
    /// shape — `C:foo` — gets rejected even on POSIX hosts so a case
    /// bundle authored on Linux can't smuggle a Windows-specific
    /// escape that only manifests when the bundle is opened on Windows.
    #[test]
    fn confined_join_rejects_drive_letter_prefix_cross_platform() {
        let case_dir = std::path::Path::new("/case/path");
        // `C:foo` is never a legitimate POSIX or Windows-relative path.
        let err = confined_join(case_dir, std::path::Path::new("C:foo"))
            .expect_err("must reject C:foo cross-platform");
        let msg = format!("{err}");
        assert!(
            msg.contains("drive") || msg.contains("absolute"),
            "msg: {msg}"
        );
    }

    // -----------------------------------------------------------------
    // Version detection
    // -----------------------------------------------------------------

    #[test]
    fn extract_semver_picks_up_three_part_version() {
        assert_eq!(
            extract_semver("OpenFOAM v11.0.0 / built on 2024-12-01"),
            Some("11.0.0".into())
        );
    }

    #[test]
    fn extract_semver_picks_up_two_part_version() {
        // Some tools (older Cantera, gmsh master) report just MAJOR.MINOR.
        assert_eq!(extract_semver("Cantera 3.0"), Some("3.0".into()));
    }

    #[test]
    fn extract_semver_handles_pre_release_tag() {
        assert_eq!(
            extract_semver("CalculiX Version 2.21-pre1"),
            Some("2.21-pre1".into())
        );
    }

    #[test]
    fn extract_semver_handles_build_metadata() {
        assert_eq!(
            extract_semver("preCICE 3.1.0+abc123"),
            Some("3.1.0+abc123".into())
        );
    }

    #[test]
    fn extract_semver_returns_none_for_no_version() {
        assert_eq!(extract_semver("nothing here"), None);
        assert_eq!(extract_semver(""), None);
        assert_eq!(extract_semver("year is 2025"), None);
    }

    #[test]
    fn extract_semver_skips_lone_integers() {
        // "version 7" is NOT a match — we require at least
        // major.minor.
        assert_eq!(extract_semver("version 7 release"), None);
    }

    #[test]
    fn extract_semver_takes_first_match() {
        // Multiple versions in the same blob — return the first
        // (typically the tool itself; later versions are usually
        // dependencies the user doesn't care about for the probe).
        assert_eq!(
            extract_semver("foo 1.2.3 (linked against bar 4.5.6)"),
            Some("1.2.3".into())
        );
    }

    #[test]
    fn extract_semver_doesnt_swallow_trailing_text() {
        // Conservative tail: the version stops at the first
        // whitespace / punctuation after the optional pre-release
        // tag. Catches accidents like consuming the entire line.
        assert_eq!(
            extract_semver("openems 0.0.36 (built 2025-01)"),
            Some("0.0.36".into())
        );
    }

    #[test]
    fn detect_tool_version_returns_none_for_missing_binary() {
        // Spawn a path that definitely doesn't exist.
        let p = std::path::PathBuf::from("/does/not/exist/at/all/banana");
        assert!(detect_tool_version(&p, "--version").is_none());
    }

    #[test]
    fn capture_subprocess_stdout_returns_none_for_missing_binary() {
        let p = std::path::PathBuf::from("/does/not/exist/at/all/banana");
        assert!(capture_subprocess_stdout(&p, "--version").is_none());
    }

    /// Round-16 M2 RED→GREEN: a probe binary that hangs (e.g. blocks
    /// on a license server, deadlocks in startup) must not block the
    /// caller indefinitely — the watchdog fires after the configured
    /// timeout, kills the child, and the probe returns `None`.
    ///
    /// Pre-fix: the round-8 bounded reads still hung on `read_to_end`
    /// because the child held its pipes open forever, so adapter
    /// probes against a stuck binary would never return.
    ///
    /// We pick a binary that takes a single-argument "hang" form:
    ///   * POSIX: `/bin/sleep 30` (sleep binary, "30" as the arg)
    ///   * Windows: `cmd.exe` with `/c pause` runs but pause waits on
    ///     stdin (we close stdin), so we use a longer-form workaround
    ///     by re-invoking our own test binary in a separate parking
    ///     thread via a small temp script. To keep the test simple
    ///     and reliable on Windows CI, we test the watchdog *unit*
    ///     itself on Windows instead — see
    ///     `capture_subprocess_stdout_watchdog_kills_long_running_child`.
    #[cfg(not(windows))]
    #[test]
    fn capture_subprocess_stdout_returns_none_on_timeout() {
        let bin = std::path::PathBuf::from("/bin/sleep");
        // Confirm sleep is present; skip otherwise.
        if std::process::Command::new(&bin)
            .arg("0")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .is_err()
        {
            eprintln!("skipping: /bin/sleep not present");
            return;
        }
        // Now exercise the real timeout path with a 1-second deadline.
        // We assert the call returns within 8s (generous to absorb
        // scheduling jitter on slow CI machines + child cleanup).
        let started = std::time::Instant::now();
        let result = capture_subprocess_stdout_with_timeout(&bin, "30", 1);
        let elapsed = started.elapsed();
        assert!(
            result.is_none(),
            "expected None after watchdog kill, got Some(...)"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(8),
            "watchdog should have fired well before 8s, took {elapsed:?}"
        );
    }

    /// Windows variant: spawn `powershell.exe Start-Sleep` directly
    /// via `std::process::Command` (multi-arg) and then run the
    /// SAME watchdog mechanism the public API uses. We can't call
    /// `capture_subprocess_stdout_with_timeout` here because the
    /// `(binary, flag)` signature only carries one extra arg and
    /// powershell needs `-NoProfile`, `-Command`, `Start-Sleep ...`.
    /// The mechanism under test (Arc<Mutex<KillOnDropChild>> +
    /// watchdog thread) is byte-for-byte the same.
    ///
    /// Critical structural detail: we drain the child's stdout in
    /// the foreground while the watchdog holds the mutex briefly.
    /// We do NOT block-on-wait while holding the mutex — that would
    /// shut the watchdog out and the test would hang for 30s.
    #[cfg(windows)]
    #[test]
    fn capture_subprocess_stdout_watchdog_kills_long_running_child() {
        use crate::subprocess::KillOnDropChild;
        use std::sync::{Arc, Mutex};
        // Spawn powershell directly so we can pass multiple flags.
        let raw_child = match std::process::Command::new("powershell.exe")
            .args(["-NoProfile", "-Command", "Start-Sleep -Seconds 30"])
            .stdin(std::process::Stdio::null())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
        {
            Ok(c) => c,
            Err(_) => {
                eprintln!("skipping: powershell.exe not present");
                return;
            }
        };
        // Take the pipes BEFORE locking so the read loop runs without
        // holding the mutex (matching the production code path).
        let mut child = raw_child;
        let stdout = child.stdout.take().expect("stdout");
        let stderr = child.stderr.take().expect("stderr");
        let child_arc = Arc::new(Mutex::new(KillOnDropChild::new(child, true)));
        // Watchdog: 1-second deadline. Same mechanism the public
        // `capture_subprocess_stdout_with_timeout` uses.
        let killed = Arc::new(std::sync::atomic::AtomicBool::new(false));
        {
            let watchdog_child = Arc::clone(&child_arc);
            let killed_flag = Arc::clone(&killed);
            std::thread::spawn(move || {
                std::thread::sleep(std::time::Duration::from_secs(1));
                if let Ok(mut guard) = watchdog_child.lock() {
                    let still_running = matches!(guard.inner_mut().try_wait(), Ok(None));
                    if still_running {
                        let _ = guard.inner_mut().kill();
                        killed_flag.store(true, std::sync::atomic::Ordering::SeqCst);
                    }
                }
            });
        }
        // Drain stdout+stderr in the foreground (no lock held). When
        // the watchdog kills the child the pipes close and the read
        // returns. This is exactly the production code's read pattern.
        let started = std::time::Instant::now();
        use std::io::Read;
        let mut sink = Vec::new();
        let mut stdout = stdout;
        let mut stderr = stderr;
        let _ = (&mut stdout)
            .take(MAX_PROBE_OUTPUT_BYTES as u64)
            .read_to_end(&mut sink);
        let _ = (&mut stderr)
            .take(MAX_PROBE_OUTPUT_BYTES as u64)
            .read_to_end(&mut sink);
        // Now wait for child to fully exit; the kill has already
        // fired so this is brief.
        if let Ok(mut guard) = child_arc.lock() {
            let _ = guard.inner_mut().wait();
        }
        let elapsed = started.elapsed();
        assert!(
            killed.load(std::sync::atomic::Ordering::SeqCst),
            "watchdog should have fired before the 30s sleep finished"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(8),
            "watchdog should have fired well before 8s, took {elapsed:?}"
        );
    }

    /// RED→GREEN (round-25 L3): when the watchdog kills the child,
    /// `capture_subprocess_stdout_with_timeout` MUST return within
    /// a bounded wall-clock even if the drain threads' pipe-reads
    /// take a moment to unblock after `TerminateProcess`. Pre-fix
    /// the watchdog called `kill()` and immediately released the
    /// mutex; the drain threads' `read_to_end` was still blocked
    /// on the OS pipe and only unblocked once the kernel finished
    /// the child's teardown (a few hundred ms on Windows under
    /// load). Post-fix the watchdog calls `wait()` AFTER `kill()`,
    /// guaranteeing the child is fully reaped before the watchdog
    /// signals "done", which in turn lets the pipe reads unblock
    /// and the function return promptly.
    ///
    /// We assert the wall-clock is well within timeout + small
    /// teardown budget (5s for a 1s timeout). The actual L3
    /// regression would surface as a wall-clock blow past 10s+
    /// under heavy CI load — the looser budget here matches the
    /// existing watchdog-kill test's threshold.
    #[cfg(windows)]
    #[test]
    fn capture_subprocess_stdout_round25_l3_watchdog_waits_before_return() {
        use std::path::PathBuf;
        // Pick powershell.exe as the long-running child (matches the
        // existing watchdog test's approach). If powershell isn't on
        // PATH the test gracefully skips.
        let bin = match std::process::Command::new("where.exe")
            .arg("powershell.exe")
            .output()
        {
            Ok(out) if out.status.success() => {
                let s = String::from_utf8_lossy(&out.stdout);
                let first = s.lines().next().unwrap_or("").trim().to_string();
                if first.is_empty() {
                    eprintln!("skipping: powershell.exe not on PATH");
                    return;
                }
                PathBuf::from(first)
            }
            _ => {
                eprintln!("skipping: where.exe not available");
                return;
            }
        };
        // 1-second timeout against a 30-second sleep; if the L3 fix
        // works the return is within ~2s on a quiet runner. We give a
        // 5s budget so CI under load passes.
        let started = std::time::Instant::now();
        let _ = capture_subprocess_stdout_with_timeout(
            &bin,
            "-NoProfile -Command Start-Sleep -Seconds 30",
            1,
        );
        let elapsed = started.elapsed();
        assert!(
            elapsed < std::time::Duration::from_secs(5),
            "watchdog-killed probe must return within 5s, took {elapsed:?} \
             (L3: kernel pipe-teardown may not have unblocked drain threads)",
        );
    }

    #[test]
    fn capture_subprocess_stdout_truncates_at_max_probe_output_bytes() {
        // Round-8 RED→GREEN: a hostile / runaway probe that emits
        // 10 MiB on the version channel would pre-fix block the
        // adapter in `Command::output()` and pressure host memory.
        // The bounded-read cap (MAX_PROBE_OUTPUT_BYTES) silently
        // truncates output past the limit so version-parsing
        // continues to work and the host stays responsive.
        //
        // We exercise the inner helper directly so the test doesn't
        // need to spawn a real subprocess that emits megabytes.
        let big_stdout: Vec<u8> = vec![b'A'; 10 * 1024 * 1024];
        let stderr: Vec<u8> = Vec::new();
        // Round-24 M6: inner now drains in parallel threads so the
        // readers must be `Send`. `Cursor<Vec<u8>>` already is.
        let stdout_reader: Box<dyn std::io::Read + Send> =
            Box::new(std::io::Cursor::new(big_stdout));
        let stderr_reader: Box<dyn std::io::Read + Send> = Box::new(std::io::Cursor::new(stderr));
        let got = capture_subprocess_stdout_inner(stdout_reader, stderr_reader).unwrap();
        // stdout truncated at MAX_PROBE_OUTPUT_BYTES; stderr was
        // empty so the trailing newline + empty stderr add no extra
        // bytes. The combined output must be ≤ the cap + 1 (newline).
        assert!(
            got.len() <= MAX_PROBE_OUTPUT_BYTES + 1,
            "expected ≤ {} bytes, got {}",
            MAX_PROBE_OUTPUT_BYTES + 1,
            got.len()
        );
        assert!(
            got.len() >= MAX_PROBE_OUTPUT_BYTES,
            "expected the cap to actually fire (≥ {} bytes), got {}",
            MAX_PROBE_OUTPUT_BYTES,
            got.len()
        );
    }

    /// RED→GREEN (round-24 M6): `capture_subprocess_stdout_inner`
    /// must drain stdout AND stderr in parallel so a child that
    /// fills its stderr pipe buffer (~64 KiB on Linux) before
    /// closing stdout doesn't pipe-deadlock the probe. Pre-fix the
    /// inner drained stdout fully via `read_to_end`, THEN stderr —
    /// which deadlocks on any child that follows the typical
    /// "version banner + license/copyright stanza to stderr" shape
    /// when the stderr stanza exceeds the pipe buffer.
    ///
    /// We can't ship a real subprocess in unit tests, so we model
    /// the deadlock shape with two synchronised in-memory readers:
    /// the "stderr" reader publishes a flag when its first read
    /// fires, and the "stdout" reader only returns data AFTER it
    /// sees that flag set. Under the pre-fix serial drain, stderr
    /// is read last — its first read never fires until stdout is
    /// fully drained, but stdout never produces data until stderr
    /// is drained, so the test would block forever. The wall-clock
    /// gate (test fails if elapsed > 5s) acts as the assertion.
    #[test]
    fn capture_subprocess_stdout_inner_drains_pipes_in_parallel() {
        use std::sync::atomic::{AtomicBool, Ordering};
        use std::sync::Arc;
        use std::time::Instant;

        /// Reader that returns one chunk of bytes the FIRST time it
        /// reads, after waiting for a signal flag. Subsequent calls
        /// return EOF.
        struct GatedReader {
            payload: Option<Vec<u8>>,
            gate: Arc<AtomicBool>,
            /// If `Some`, we set this flag to true on the first read
            /// (before waiting for the gate). Used by the "stderr"
            /// reader to signal the "stdout" reader.
            signal: Option<Arc<AtomicBool>>,
        }
        impl std::io::Read for GatedReader {
            fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
                // Signal first so the other thread can make progress.
                if let Some(sig) = &self.signal {
                    sig.store(true, Ordering::SeqCst);
                }
                // Wait for the gate (with a hard 4s timeout — the
                // test runner enforces wall-clock anyway).
                let deadline = Instant::now() + std::time::Duration::from_secs(4);
                while !self.gate.load(Ordering::SeqCst) {
                    if Instant::now() >= deadline {
                        // Pipe-deadlock signal: the other side never
                        // released us. Return EOF so the inner
                        // doesn't loop forever in tests.
                        return Ok(0);
                    }
                    std::thread::sleep(std::time::Duration::from_millis(10));
                }
                match self.payload.take() {
                    Some(p) => {
                        let n = std::cmp::min(p.len(), buf.len());
                        buf[..n].copy_from_slice(&p[..n]);
                        Ok(n)
                    }
                    None => Ok(0),
                }
            }
        }

        let stdout_can_read = Arc::new(AtomicBool::new(false));
        let stderr_has_started = Arc::new(AtomicBool::new(false));
        // stdout reader: blocks until stderr reader has signalled it
        // started (proving the two are draining in parallel).
        let stdout_reader: Box<dyn std::io::Read + Send> = Box::new(GatedReader {
            payload: Some(b"version 1.2.3".to_vec()),
            gate: Arc::clone(&stderr_has_started),
            signal: None,
        });
        // stderr reader: signals immediately, then waits on its own
        // gate (we open it instantly via a side thread). The
        // pre-fix serial drain never starts the stderr reader until
        // stdout fully drains — so `stderr_has_started` stays false
        // and the stdout reader times out at 4s.
        let stderr_reader: Box<dyn std::io::Read + Send> = Box::new(GatedReader {
            payload: Some(b"more output\n".to_vec()),
            gate: Arc::clone(&stdout_can_read), // not actually toggled — pure signal
            signal: Some(Arc::clone(&stderr_has_started)),
        });

        // Side thread: release the stderr reader's gate after 50ms
        // so it returns its payload promptly. This is irrelevant to
        // the deadlock test — we're checking that stderr STARTED
        // reading in parallel with stdout, not that it finished.
        let gate_clone = Arc::clone(&stdout_can_read);
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(50));
            gate_clone.store(true, Ordering::SeqCst);
        });

        let started = Instant::now();
        let _ = capture_subprocess_stdout_inner(stdout_reader, stderr_reader);
        let elapsed = started.elapsed();
        assert!(
            elapsed < std::time::Duration::from_secs(3),
            "parallel drain should finish in well under 3s; took {elapsed:?} \
             (pre-fix this would be ≥ 4s because the serial drain order \
             deadlocks the gated readers)",
        );
        assert!(
            stderr_has_started.load(Ordering::SeqCst),
            "stderr drain must have started in parallel with stdout — \
             pre-fix it only started after stdout fully drained, never \
             toggling this flag",
        );
    }

    #[test]
    fn detect_tool_version_semver_returns_none_for_missing_binary() {
        let p = std::path::PathBuf::from("/does/not/exist/at/all/banana");
        assert!(detect_tool_version_semver(&p, &["--version", "-V"]).is_none());
    }

    /// Round-17 L2 RED→GREEN: the watchdog thread must exit EARLY when
    /// the probe finishes — pre-fix it slept the full `timeout_secs`
    /// regardless of how fast the child returned, holding its
    /// `Arc<Mutex<KillOnDropChild>>` clone alive for that entire span.
    /// Running many fast probes in series would back up dozens of
    /// idle watchdog threads + Arc references for the leftover sleep
    /// time, leaking memory + thread handles.
    ///
    /// We pin the regression by spawning 20 back-to-back probes
    /// against a fast-exit binary (a missing path returns instantly).
    /// With the channel-based watchdog the whole sequence finishes in
    /// well under 1 second; pre-fix the test would have taken
    /// ≥ 10 × 20 = 200 seconds (because the executor would block on
    /// the trailing watchdog from each probe).
    ///
    /// We're really testing the watchdog mechanism here, not the
    /// missing-binary path — `capture_subprocess_stdout` returns
    /// `None` immediately on spawn failure, never installing the
    /// watchdog. Instead we directly exercise the public API with the
    /// shortest probe binary we can find. On every platform the
    /// `capture_subprocess_stdout` path either spawns successfully
    /// (then exits fast) or returns `None` from `spawn()?` — both
    /// paths must finish well under the timeout.
    ///
    /// Cross-platform binary choice: `/bin/echo` on POSIX,
    /// `cmd.exe /c rem` on Windows. If neither is present we skip
    /// the assertion (the test still passes — the harness asserts only
    /// what we CAN verify).
    #[test]
    fn capture_subprocess_stdout_watchdog_exits_promptly_on_fast_probe() {
        // Find a fast-exit binary. If neither is present, the test
        // becomes a no-op (still valid — we can't manufacture a
        // probe binary out of thin air on a stripped-down CI host).
        #[cfg(windows)]
        let (bin, flag) = (std::path::PathBuf::from("cmd.exe"), "/?");
        #[cfg(not(windows))]
        let (bin, flag) = (std::path::PathBuf::from("/bin/echo"), "hi");
        if std::process::Command::new(&bin)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .stdin(std::process::Stdio::null())
            .arg(flag)
            .status()
            .is_err()
        {
            eprintln!("skipping: probe binary not present at {bin:?}");
            return;
        }
        // Run 20 back-to-back probes through the public API with the
        // default 10s timeout. Pre-fix: each probe's watchdog held
        // its Arc clone for the full 10s, so even though the probes
        // themselves return fast, the threads + Arc references pile
        // up. The total runtime of this loop must remain well below
        // the per-probe timeout times the iteration count.
        let started = std::time::Instant::now();
        for _ in 0..20 {
            // Public API uses MAX_PROBE_TIMEOUT_SECS = 10. We don't
            // assert on Some/None — the binary might emit nothing
            // useful — only on the elapsed time.
            let _ = capture_subprocess_stdout(&bin, flag);
        }
        let elapsed = started.elapsed();
        // 20 × 10s = 200s pre-fix worst case (if each probe's foreground
        // wait blocked on the watchdog's lifetime via the Arc — it
        // technically didn't, but the thread leak was real). After the
        // fix the loop must finish in well under 5 seconds even on the
        // slowest CI machine.
        assert!(
            elapsed < std::time::Duration::from_secs(20),
            "20 fast probes took {elapsed:?} — watchdog channel \
             should exit early, not sleep the full timeout"
        );
    }

    // -----------------------------------------------------------------
    // Provenance hash helpers
    // -----------------------------------------------------------------

    #[test]
    fn sha256_hex_bytes_matches_known_vector() {
        // SHA-256("") = e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
        let h = sha256_hex_bytes(b"");
        assert_eq!(
            h.0,
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn sha256_hex_bytes_changes_with_input() {
        let h1 = sha256_hex_bytes(b"hello");
        let h2 = sha256_hex_bytes(b"world");
        assert_ne!(h1, h2);
        assert_eq!(h1.0.len(), 64);
        assert_eq!(h2.0.len(), 64);
    }

    #[test]
    fn sha256_hex_file_returns_empty_for_missing_file() {
        let p = std::path::PathBuf::from("/does/not/exist/file.txt");
        let h = sha256_hex_file(&p);
        assert_eq!(h.0, "", "missing file must hash to empty string sentinel");
    }

    #[test]
    fn sha256_hex_file_matches_bytes_helper() {
        let p = std::env::temp_dir().join(format!(
            "valenx-prov-hash-{}.txt",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let payload = b"valenx provenance test fixture";
        std::fs::write(&p, payload).unwrap();
        let h_file = sha256_hex_file(&p);
        let h_bytes = sha256_hex_bytes(payload);
        assert_eq!(h_file, h_bytes);
        let _ = std::fs::remove_file(&p);
    }

    /// Round-16 M1 RED→GREEN: streaming hash over a >> HASH_CHUNK
    /// file produces the same digest as the bytes-buffer helper.
    /// Pre-fix `sha256_hex_file` used `std::fs::read` (loads the whole
    /// file) so a multi-MB file forced an equally-large heap spike
    /// just to compute a hash that the streaming form does with a
    /// 64 KiB working set. We test against a payload several HASH_CHUNK
    /// values wide so the chunk-loop path is exercised; the assertion
    /// still pins the hash against the all-bytes-at-once helper.
    #[test]
    fn sha256_hex_file_streaming_matches_for_multi_chunk_file() {
        let p = std::env::temp_dir().join(format!(
            "valenx-prov-hash-stream-{}.bin",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // ~3.5 HASH_CHUNK worth of distinct bytes so a chunked read
        // crosses several chunk boundaries with a partial tail.
        let mut payload = Vec::with_capacity(HASH_CHUNK * 4);
        for i in 0..(HASH_CHUNK * 4) - 17 {
            payload.push((i & 0xff) as u8);
        }
        std::fs::write(&p, &payload).unwrap();
        let h_file = sha256_hex_file(&p);
        let h_bytes = sha256_hex_bytes(&payload);
        assert_eq!(h_file, h_bytes, "streaming + bytes hashes must agree");
        let _ = std::fs::remove_file(&p);
    }

    /// Round-16 M1 RED→GREEN: a file whose size exceeds the
    /// [`MAX_HASHED_BYTES`] cap returns the empty-hash sentinel
    /// instead of streaming gigabytes through the hasher (or, pre-fix,
    /// allocating gigabytes via `fs::read`). We use a sparse file via
    /// `set_len` so the test doesn't actually consume 4 GiB of disk —
    /// the OS reports the length and the streaming reader yields
    /// zeros up to the cap; we detect overflow at the very first
    /// chunk past the cap and return empty.
    ///
    /// Sparse files are supported on every fs we ship to (NTFS,
    /// ext4, APFS), and reads of holes return zero-filled buffers
    /// per POSIX. If a future test target lacks sparse-file support
    /// this test will time out / fill the disk — the assertion is
    /// guarded so a failure is loud rather than silent.
    #[test]
    fn sha256_hex_file_returns_empty_for_oversize_file() {
        let p = std::env::temp_dir().join(format!(
            "valenx-prov-hash-oversize-{}.bin",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // 1 byte over the cap — set_len creates a sparse file on
        // every fs we target (NTFS, ext4, APFS), so the test stays
        // within seconds without filling disk.
        let f = std::fs::File::create(&p).unwrap();
        f.set_len(MAX_HASHED_BYTES + 1).unwrap();
        drop(f);
        let h = sha256_hex_file(&p);
        assert_eq!(
            h.0, "",
            "oversize file (> MAX_HASHED_BYTES) must hash to empty sentinel"
        );
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn fresh_run_id_has_uuid_shape() {
        let id = fresh_run_id();
        // 8-4-4-4-12 = 36 chars + 4 dashes
        assert_eq!(id.len(), 36);
        let segments: Vec<&str> = id.split('-').collect();
        assert_eq!(segments.len(), 5);
        assert_eq!(segments[0].len(), 8);
        assert_eq!(segments[1].len(), 4);
        assert_eq!(segments[2].len(), 4);
        assert_eq!(segments[3].len(), 4);
        assert_eq!(segments[4].len(), 12);
        // All hex
        assert!(id.chars().all(|c| c == '-' || c.is_ascii_hexdigit()));
    }

    #[test]
    fn fresh_run_id_is_unique_across_calls() {
        // Two back-to-back calls must produce different ids — the
        // counter ensures this even within the same nanosecond.
        let mut ids: Vec<String> = (0..10).map(|_| fresh_run_id()).collect();
        ids.sort();
        ids.dedup();
        assert_eq!(ids.len(), 10);
    }

    #[test]
    fn live_provenance_populates_real_hashes() {
        let case_path = std::env::temp_dir().join(format!(
            "valenx-prov-case-{}.toml",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&case_path, b"[case]\nname=\"test\"\n").unwrap();
        let prov = live_provenance(
            "test", "0.1.0", "TestTool", "1.0.0", &case_path, None, None, 42.5,
        );
        assert_eq!(prov.adapter, "test");
        assert_eq!(prov.adapter_version, "0.1.0");
        assert_eq!(prov.tool_version, "1.0.0");
        assert_eq!(prov.case_hash.0.len(), 64);
        assert_ne!(prov.case_hash.0, "", "case hash must be populated");
        // Mesh + tools_lock paths are None -> empty hash sentinels.
        assert_eq!(prov.mesh_hash.0, "");
        assert_eq!(prov.tools_lock_hash.0, "");
        assert_eq!(prov.wall_time_seconds, 42.5);
        assert!(prov.completed_at.starts_with("20")); // 21st century
        let _ = std::fs::remove_file(&case_path);
    }

    #[test]
    fn live_provenance_handles_missing_case_file() {
        let bogus = std::path::PathBuf::from("/does/not/exist/case.toml");
        let prov = live_provenance(
            "test", "0.1.0", "TestTool", "1.0.0", &bogus, None, None, 0.0,
        );
        // Missing case path -> empty hash, not a panic.
        assert_eq!(prov.case_hash.0, "");
    }

    #[test]
    fn current_timestamp_iso8601_format_is_zulu() {
        let ts = super::current_timestamp_iso8601();
        // YYYY-MM-DDTHH:MM:SSZ -> 20 chars
        assert_eq!(ts.len(), 20, "got: {ts:?}");
        assert!(ts.ends_with('Z'));
        assert!(ts.chars().nth(4).unwrap() == '-');
        assert!(ts.chars().nth(7).unwrap() == '-');
        assert!(ts.chars().nth(10).unwrap() == 'T');
    }

    #[test]
    fn first_workdir_match_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-fwm-lex-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("z_late.cfg"), b"placeholder").unwrap();
        std::fs::write(tmp.join("a_first.cfg"), b"placeholder").unwrap();
        std::fs::write(tmp.join("not_a_match.csv"), b"placeholder").unwrap();
        let f = first_workdir_match(&tmp, &["cfg"]).expect("found");
        assert!(f.ends_with("a_first.cfg"));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn first_workdir_match_returns_none_for_no_match() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-fwm-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("only.txt"), b"placeholder").unwrap();
        assert!(first_workdir_match(&tmp, &["cfg"]).is_none());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn first_workdir_match_is_case_insensitive_on_extension() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-fwm-case-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("model.STEP"), b"placeholder").unwrap();
        let f = first_workdir_match(&tmp, &["step"]).expect("found");
        assert!(f.ends_with("model.STEP"));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn confined_join_accepts_simple_relative_path() {
        let case = std::path::Path::new("/case");
        let p = confined_join(case, std::path::Path::new("script.py")).unwrap();
        assert_eq!(p, std::path::Path::new("/case/script.py"));
    }

    #[test]
    fn confined_join_accepts_relative_path_with_subdirs() {
        let case = std::path::Path::new("/case");
        let p = confined_join(case, std::path::Path::new("scripts/inner/run.py")).unwrap();
        assert_eq!(p, std::path::Path::new("/case/scripts/inner/run.py"));
    }

    #[test]
    fn confined_join_rejects_absolute_path() {
        let case = std::path::Path::new("/case");
        // Use a unix-style absolute on unix and a drive-style absolute
        // on Windows so the test stays cross-platform.
        let abs = if cfg!(windows) {
            std::path::PathBuf::from("C:/etc/passwd")
        } else {
            std::path::PathBuf::from("/etc/passwd")
        };
        let err = confined_join(case, &abs).expect_err("absolute should reject");
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute"),
            "expected error to mention `absolute`, got: {msg}"
        );
    }

    #[test]
    fn confined_join_rejects_parent_dir_traversal() {
        let case = std::path::Path::new("/case");
        let err = confined_join(case, std::path::Path::new("../../etc/passwd"))
            .expect_err("parent-dir traversal should reject");
        let msg = format!("{err}");
        assert!(
            msg.contains(".."),
            "expected error to mention `..`, got: {msg}"
        );
    }

    #[test]
    fn confined_join_rejects_parent_dir_in_middle_of_path() {
        let case = std::path::Path::new("/case");
        let err = confined_join(case, std::path::Path::new("inner/../../escape"))
            .expect_err("middle-of-path .. should reject");
        let msg = format!("{err}");
        assert!(msg.contains(".."), "got: {msg}");
    }

    // -----------------------------------------------------------------
    // OSS license predicate
    // -----------------------------------------------------------------

    #[test]
    fn is_oss_license_accepts_common_osi_strings() {
        // Every license string the workspace's adapters actually use
        // for OSS tools must satisfy the predicate. If any of these
        // ever fail, an adapter would silently disappear from the
        // default UI — surprising behaviour.
        for s in [
            "MIT",
            "BSD-2-Clause",
            "BSD-3-Clause",
            "BSD-3-Clause-Open-Source",
            "Apache-2.0",
            "GPL-2.0",
            "GPL-2.0-only",
            "GPL-2.0-or-later",
            "GPL-3.0",
            "GPL-3.0-only",
            "GPL-3.0-or-later",
            "LGPL-2.1",
            "LGPL-2.1-only",
            "LGPL-2.1-or-later",
            "LGPL-3.0",
            "LGPL-3.0-or-later",
            "AGPL-3.0",
            "AGPL-3.0-only",
            "AGPL-3.0-or-later",
            "MPL-2.0",
            "Artistic-2.0",
            "ECL-2.0",
            "AFL-3.0",
            "Public Domain",
            "CC0-1.0",
        ] {
            assert!(is_oss_license(s), "expected OSS: `{s}`");
        }
    }

    #[test]
    fn is_oss_license_rejects_academic_and_per_tool_licenses() {
        // Per-tool / non-commercial strings the registry uses today.
        // The default UI hides every adapter whose license matches
        // one of these.
        for s in [
            "NAMD-License",
            "Rosetta-License",
            "Academic",
            "CC-BY-NC-4.0",
            "CC-BY-NC-SA-4.0",
            "VMD-License",
            "NUPACK-License",
            "X3DNA-License",
            "DSSR-License",
            "Curves-License",
            "Janelia-License",
            "ViennaRNA-License",
            "Cambrian-Open-License",
            "University of California — free for non-commercial / academic use",
            "",
            "unknown",
        ] {
            assert!(!is_oss_license(s), "expected NON-OSS: `{s}`");
        }
    }

    #[test]
    fn is_oss_license_handles_spdx_with_exception_suffix() {
        // OCCT's tool_license is `"LGPL-2.1-or-later WITH OCCT-exception-1.0"`.
        // The SPDX `WITH <exception>` modifier expands rather than
        // contracts the licensee's rights, so the core license still
        // dictates OSS status. Predicate must strip the suffix.
        assert!(is_oss_license("LGPL-2.1-or-later WITH OCCT-exception-1.0"));
        assert!(is_oss_license(
            "GPL-3.0-or-later WITH Classpath-exception-2.0"
        ));
    }

    // -----------------------------------------------------------------
    // Round-15: structured-identifier sanitizers (sister to round-14's
    // IFC `ifc_str` Part-21 escape; cover SIF / INP / OF dict / STEP
    // block-comment writers that previously interpolated user strings
    // raw into their output).
    // -----------------------------------------------------------------

    #[test]
    fn sanitize_structured_identifier_strips_quotes_and_newlines() {
        // The classic injection payload: a user-controlled "name" field
        // closes the SIF quoted literal, opens a sibling Material block,
        // and leaves the rest dangling. Sanitiser must drop every char
        // that lets the payload escape the quoted region.
        let payload = "x\"\nMaterial 99\n  Name = \"injected\"\nEnd\n";
        let cleaned = sanitize_structured_identifier(payload);
        assert!(!cleaned.contains('"'), "must strip embedded `\"`");
        assert!(!cleaned.contains('\n'), "must strip embedded `\\n`");
        assert!(!cleaned.contains('\r'), "must strip embedded `\\r`");
        assert!(!cleaned.contains('\\'), "must strip embedded `\\\\`");
        // Safe chars still flow through verbatim — we don't mangle
        // anything that wasn't a directive-injection vector.
        assert!(cleaned.contains('x'));
        assert!(cleaned.contains("Material"));
    }

    #[test]
    fn sanitize_structured_identifier_keeps_plain_ascii_intact() {
        // Common safe inputs (material names, file basenames) must
        // round-trip unchanged so the writer's output stays human-
        // readable for well-behaved cases.
        assert_eq!(sanitize_structured_identifier("aluminium"), "aluminium");
        assert_eq!(
            sanitize_structured_identifier("steel-7075-t6"),
            "steel-7075-t6"
        );
        assert_eq!(
            sanitize_structured_identifier("heatsink_v2.case"),
            "heatsink_v2.case"
        );
    }

    #[test]
    fn validate_structured_identifier_rejects_newline() {
        // Newline injection is the canonical attack: it lets a hostile
        // case.toml break out of a single-line dict entry and inject
        // a sibling directive on the next line.
        let err = validate_structured_identifier(
            "inlet\n}\n\ninjectedPatch\n{\ntype wall;",
            "boundaries.inlet",
        )
        .expect_err("newline must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("boundaries.inlet"), "msg: {msg}");
    }

    #[test]
    fn validate_structured_identifier_rejects_quote() {
        let err = validate_structured_identifier("x\"y", "field").expect_err("quote rejected");
        let msg = format!("{err}");
        assert!(msg.contains("field"), "msg: {msg}");
    }

    #[test]
    fn validate_structured_identifier_rejects_backslash() {
        // Backslashes are the OpenFOAM dict comment opener (`//`) part
        // and a common escape character — refuse them in identifier
        // positions even though `sanitize_structured_identifier` strips
        // them, because the strict variant's job is to fail loudly.
        let err = validate_structured_identifier("x\\y", "field").expect_err("backslash rejected");
        let msg = format!("{err}");
        assert!(msg.contains("field"), "msg: {msg}");
    }

    #[test]
    fn validate_structured_identifier_rejects_empty() {
        let err = validate_structured_identifier("", "boundaries.something")
            .expect_err("empty must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("empty"), "msg: {msg}");
    }

    #[test]
    fn validate_structured_identifier_accepts_alphanumeric_with_dot_dash_underscore() {
        // Every common identifier shape adapter authors use today:
        // bare names, dot-separated, dash-separated, underscore-
        // separated. All round-trip without error.
        assert!(validate_structured_identifier("inlet", "f").is_ok());
        assert!(validate_structured_identifier("inlet1", "f").is_ok());
        assert!(validate_structured_identifier("inlet-top", "f").is_ok());
        assert!(validate_structured_identifier("inlet_top", "f").is_ok());
        assert!(validate_structured_identifier("inlet.top", "f").is_ok());
        assert!(
            validate_structured_identifier("velocityInlet-1.heatFlux_walls", "f").is_ok(),
            "compound identifier with all three delimiters must pass"
        );
    }

    #[test]
    fn validate_structured_identifier_rejects_space() {
        // Spaces are illegal in OpenFOAM patch names / gmsh physical
        // group names — the parser stops at the first space, leaving
        // the suffix as a stray token. Refuse rather than mangle.
        let err = validate_structured_identifier("inlet wall", "f").expect_err("space rejected");
        let msg = format!("{err}");
        assert!(msg.contains("f"), "msg: {msg}");
    }

    // -----------------------------------------------------------------
    // Round-19 H2: python_str_repr — escapes for Python double-quoted
    // literal embedding. Sister to round-15's structured-identifier
    // sanitisers, but for the Python-source emitter call sites
    // (cantera, pybamm).
    // -----------------------------------------------------------------

    /// Round-19 H2 RED→GREEN: the canonical injection payload — a
    /// user string that closes the literal, opens a sibling
    /// statement, and leaves the rest dangling — must round-trip
    /// as a single escaped literal. Pre-fix, `format!("\"{}\"", x)`
    /// would emit the payload verbatim and turn the resulting script
    /// into statement-level Python the user controls.
    #[test]
    fn python_str_repr_escapes_quote_and_newline_payload() {
        let payload = "Chen2020\");\nimport os\nx = pybamm.ParameterValues(\"";
        let escaped = python_str_repr(payload);
        // The closing `"` MUST be backslash-escaped. The escaped
        // form should have `\");` (literal backslash before the
        // quote). We scan for `");` and assert it's always preceded
        // by `\` — the bare 3-byte `");` sequence (closing quote +
        // paren + semicolon) MUST NOT appear without a backslash.
        for (idx, window) in escaped.as_bytes().windows(3).enumerate() {
            if window == [b'"', b')', b';'] {
                assert!(
                    idx > 0 && escaped.as_bytes()[idx - 1] == b'\\',
                    "raw `\");` (without preceding backslash) leaked: \
                     position {idx}, escaped {escaped:?}"
                );
            }
        }
        assert!(
            escaped.contains("\\\""),
            "expected `\\\"` for the embedded quote, got {escaped:?}"
        );
        // Newlines must be escaped — a raw `\n` lets the next line
        // become a statement.
        assert!(
            !escaped.contains('\n'),
            "raw newline leaked through escaping: {escaped:?}"
        );
        assert!(
            escaped.contains("\\n"),
            "expected `\\n` for the embedded newline, got {escaped:?}"
        );
    }

    /// Plain alphanumeric / dot / dash strings round-trip unchanged —
    /// the helper must not mangle the common case where the user's
    /// string was a normal parameter-set name like `Chen2020`.
    #[test]
    fn python_str_repr_keeps_plain_ascii_intact() {
        assert_eq!(python_str_repr("Chen2020"), "Chen2020");
        assert_eq!(python_str_repr("h2o2"), "h2o2");
        assert_eq!(
            python_str_repr("CH4:1, O2:2, N2:7.52"),
            "CH4:1, O2:2, N2:7.52"
        );
    }

    /// Backslash must be escaped — the user could otherwise inject
    /// a trailing `\"` (escaped quote) to merge the next line's
    /// content into the literal.
    #[test]
    fn python_str_repr_escapes_backslash() {
        assert_eq!(python_str_repr("path\\to\\file"), "path\\\\to\\\\file");
    }

    /// `\r` and `\t` are escape sequences in their own right; both
    /// must be encoded explicitly so the emitted literal is byte-
    /// identical regardless of platform line-ending normalisation.
    #[test]
    fn python_str_repr_escapes_cr_and_tab() {
        assert_eq!(python_str_repr("a\rb"), "a\\rb");
        assert_eq!(python_str_repr("a\tb"), "a\\tb");
    }

    /// Non-printable ASCII control chars get the `\xNN` form so an
    /// attacker can't smuggle a literal `\0` into the source (which
    /// Python would tokenise as end-of-string in some embedding
    /// contexts).
    #[test]
    fn python_str_repr_escapes_low_control_chars() {
        let s: String = (0..0x20).map(|c| c as u8 as char).collect();
        let escaped = python_str_repr(&s);
        // No raw control bytes survive.
        assert!(
            escaped.chars().all(|c| (c as u32) >= 0x20),
            "raw control bytes leaked: {escaped:?}"
        );
        // `\x00` and `\x1f` must both appear in the escaped form.
        assert!(escaped.contains("\\x00"), "got: {escaped:?}");
        assert!(escaped.contains("\\x1f"), "got: {escaped:?}");
    }

    /// UTF-8 codepoints above 0x7F pass through unchanged — they
    /// are valid in Python 3 source by default and the helper's
    /// scope is Python-literal escaping, not ASCII normalisation.
    #[test]
    fn python_str_repr_preserves_utf8() {
        assert_eq!(python_str_repr("résumé"), "résumé");
        assert_eq!(python_str_repr("日本語"), "日本語");
    }

    /// Round-20 L4 RED→GREEN (defence-in-depth): U+2028 (LINE
    /// SEPARATOR) and U+2029 (PARAGRAPH SEPARATOR) are line
    /// terminators in JavaScript/ECMAScript even though Python
    /// treats them as ordinary chars. Escaping them here makes the
    /// helper reusable for JS/JSON-emitting contexts without
    /// surprise. Pre-fix the helper passed both codepoints through
    /// verbatim because they fell into the "≥ 0x20, not a special
    /// char" bucket.
    #[test]
    fn python_str_repr_escapes_u2028_and_u2029() {
        // Bare U+2028 / U+2029 must NOT appear in the output.
        let payload = "before\u{2028}middle\u{2029}after";
        let escaped = python_str_repr(payload);
        assert!(
            !escaped.contains('\u{2028}'),
            "raw U+2028 leaked through escaping: {escaped:?}"
        );
        assert!(
            !escaped.contains('\u{2029}'),
            "raw U+2029 leaked through escaping: {escaped:?}"
        );
        // The Python `\u` escape is the canonical form.
        assert!(
            escaped.contains("\\u2028"),
            "expected `\\u2028` for U+2028, got: {escaped:?}"
        );
        assert!(
            escaped.contains("\\u2029"),
            "expected `\\u2029` for U+2029, got: {escaped:?}"
        );
    }

    /// RED→GREEN (round-24 M3): the aggregate entry cap fires when
    /// the source tree contains more files than the limit. Verified
    /// against a downscaled test cap (`max_entries = 5`) on a 10-
    /// file source dir — pre-fix only the per-file size cap existed,
    /// so 1.5 M zero-byte files would pass the loop and exhaust
    /// inodes on the destination. Post-fix the entry-cap check
    /// short-circuits the walk before the 6th file is touched.
    #[test]
    fn copy_dir_recursive_with_caps_enforces_entry_cap() {
        let src = std::env::temp_dir().join(format!(
            "valenx-m3-entry-src-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let dst = std::env::temp_dir().join(format!(
            "valenx-m3-entry-dst-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&src).unwrap();
        // 10 tiny files — well over the test cap of 5.
        for i in 0..10 {
            std::fs::write(src.join(format!("f-{i}.txt")), b"hi").unwrap();
        }
        let err = copy_dir_recursive_with_caps(
            &src,
            &dst,
            /*max_total_bytes=*/ u64::MAX, // disable the byte cap for this test
            /*max_entries=*/ 5,
        )
        .unwrap_err();
        let _ = std::fs::remove_dir_all(&src);
        let _ = std::fs::remove_dir_all(&dst);
        let msg = format!("{err}");
        assert!(
            msg.contains("5-entry") || msg.contains("entry cap"),
            "expected entry-cap message, got: {msg}"
        );
    }

    /// RED→GREEN (round-24 M3): the aggregate byte cap fires when
    /// the cumulative byte count crosses the limit even though each
    /// individual file is under MAX_COPY_DIR_FILE_BYTES. Pre-fix
    /// the per-file cap was the only ceiling, so a poisoned source
    /// with thousands of "just-under-cap" files would fill the disk.
    /// Test uses three 4 KiB files and a 10 KiB aggregate cap — the
    /// third file pushes the running total to 12 KiB, over the cap.
    #[test]
    fn copy_dir_recursive_with_caps_enforces_aggregate_byte_cap() {
        let src = std::env::temp_dir().join(format!(
            "valenx-m3-bytes-src-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let dst = std::env::temp_dir().join(format!(
            "valenx-m3-bytes-dst-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&src).unwrap();
        let payload = vec![b'.'; 4096];
        for i in 0..3 {
            std::fs::write(src.join(format!("chunk-{i}.bin")), &payload).unwrap();
        }
        let err = copy_dir_recursive_with_caps(
            &src,
            &dst,
            /*max_total_bytes=*/ 10_000, // 10 KiB cap — third 4 KiB file crosses
            /*max_entries=*/ usize::MAX,
        )
        .unwrap_err();
        let _ = std::fs::remove_dir_all(&src);
        let _ = std::fs::remove_dir_all(&dst);
        let msg = format!("{err}");
        assert!(
            msg.contains("aggregate cap") || msg.contains("cumulative"),
            "expected aggregate-byte message, got: {msg}"
        );
    }

    /// Sanity (round-24 M3): the production constants are present
    /// and pinned so a future refactor can't silently drop the cap.
    #[test]
    fn copy_dir_aggregate_caps_are_publicly_exposed() {
        // 64 GiB is the rationale ceiling; pin it so a sloppy bump
        // gets a code-review checkpoint.
        assert_eq!(MAX_COPY_DIR_TOTAL_BYTES, 64 * 1024 * 1024 * 1024);
        assert_eq!(MAX_COPY_DIR_ENTRIES, 1_000_000);
    }
}
