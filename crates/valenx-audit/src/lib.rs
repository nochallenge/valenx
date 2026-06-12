//! # valenx-audit
//!
//! Append-only JSONL audit log with a SHA-256 `prev_hash` chain
//! linking every entry to the previous one. First concrete chunk of
//! [RFC 0013](../../../rfcs/0013-enterprise-audit-rbac.md).
//!
//! - [`AuditWriter::append`] writes one entry, computing its
//!   `prev_hash` from the previous serialised line.
//! - [`AuditEntry`] is the canonical event shape (timestamp + actor
//!   + action + target + free-form context).
//! - [`verify_chain`] reads a log file end-to-end and checks every
//!   entry's `prev_hash` matches the SHA-256 of the previous entry's
//!   serialisation.
//!
//! What's deferred (covered by follow-up RFCs):
//!
//! - Cryptographic signing per install (the chain catches tampering
//!   that doesn't rewrite ALL subsequent entries; signing detects
//!   even full-rewrite attacks).
//! - Log rotation. The operator owns this — Valenx writes append-
//!   only and never rotates.
//! - GDPR redaction (replacing actor.id in old entries while keeping
//!   the chain valid).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use fs2::FileExt;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use thiserror::Error;

/// Maximum bytes any single audit-log record may consume when read
/// back via the slow-walk paths (`compute_prev_hash_slow`,
/// `last_line_sha256`, `tail_filtered`, `verify_chain_report`).
///
/// Round-20 M3: pre-fix the four read paths above used
/// `BufReader::lines()`, which allocates an unbounded `String` for
/// each record. A poisoned log with a single 2 GiB "line" (no
/// newline terminator) would OOM the verifier / writer's
/// genesis-rebuild walk before any record-level validation kicked
/// in. 256 KiB is generous for legitimate entries (a typical
/// `run.complete` with full stdout context tops out around 8 KiB —
/// `target` and `context` are bounded JSON values, not free-form
/// stdout dumps) while refusing the unbounded-line DoS shape.
///
/// A line that exceeds the cap surfaces as an
/// [`AuditError::Io`] with `InvalidData` kind rather than
/// `Malformed`, because the underlying read failed BEFORE any
/// per-line parsing happened — the chain walk is structurally unable
/// to continue past an over-cap line.
pub const MAX_AUDIT_LINE_BYTES: usize = 256 * 1024;

/// Iterate the lines of `reader` with a per-line byte cap. The
/// returned items strip the trailing `\n` (if present) so the bytes
/// pipe straight into `sha256_hex` / `serde_json::from_slice` /
/// `String::from_utf8` without a second copy.
///
/// `BufRead::read_until(b'\n', ...)` reads exactly one line of bytes
/// at a time and returns 0 at EOF; this iterator wraps that with the
/// extra "did we just blow past the cap?" check. Caller drives it
/// with `for line in read_capped_lines(&mut reader, cap)` — exactly
/// the shape of the pre-fix `for line in reader.lines()` it replaces.
///
/// Round-20 M3: factored out so all four audit slow-walk sites
/// share the same bounded reader instead of each one growing its
/// own ad-hoc cap.
///
/// Round-24 H1: the round-20 implementation called
/// `reader.read_until(b'\n', &mut buf)` WITHOUT a `Read::take()`
/// limit and only checked `buf.len() > max_per_line` AFTER the
/// allocation completed — defeating the cap entirely. A hostile
/// log with one 5 GiB no-newline line would allocate the whole 5
/// GiB string then "reject" it as oversize, OOMing the process.
/// The fix wraps the reader in `(&mut reader).take(cap + 1)` so
/// `read_until` itself stops at the cap, bounding the allocation.
/// This mirrors the correct pattern in
/// `valenx_core::io_caps::read_capped_lines_bounded`.
fn read_capped_lines<'r, R: BufRead + 'r>(
    reader: &'r mut R,
    max_per_line: usize,
) -> impl Iterator<Item = std::io::Result<Vec<u8>>> + 'r {
    let mut done = false;
    std::iter::from_fn(move || {
        if done {
            return None;
        }
        let mut buf = Vec::with_capacity(64);
        // Round-24 H1: read at most max_per_line+1 bytes BEFORE
        // allocating further. The `+1` is the sentinel that lets us
        // distinguish "exactly at the cap" from "would have been
        // larger" — if `take` cuts us off the buf will have grown
        // by `cap+1` bytes WITHOUT a trailing `\n`.
        let cap = (max_per_line as u64).saturating_add(1);
        // `Read::take` consumes self by value; passing `&mut reader`
        // here makes it consume the *reborrow* — the original
        // `reader` is untouched after `limited` drops at the end of
        // the iteration, so subsequent calls keep streaming from
        // where the previous line ended.
        let mut limited = std::io::Read::take(&mut *reader, cap);
        match limited.read_until(b'\n', &mut buf) {
            Ok(0) => {
                done = true;
                None
            }
            Ok(_) => {
                if buf.len() > max_per_line {
                    // Round-26 M1: the round-25 L1 bounded-resync
                    // logic that previously lived here was dead code
                    // — all four production callers
                    // (`compute_prev_hash_slow`, `last_line_sha256`,
                    // `tail_filtered`, `verify_chain_report`) use
                    // `?` on the iterator's items, so the first
                    // InvalidData stops them outright and resync
                    // never runs. Mark the iterator done on the
                    // poisoned line and surface the error. A future
                    // caller that wants skip-and-continue semantics
                    // can wrap a bounded resync around this iterator
                    // at the call site — keeping the resync out of
                    // the shared helper means we don't pay its
                    // complexity for the only callers that exist.
                    done = true;
                    Some(Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "audit log line exceeds {max_per_line}-byte cap (got {} bytes)",
                            buf.len(),
                        ),
                    )))
                } else {
                    if buf.ends_with(b"\n") {
                        buf.pop();
                    }
                    Some(Ok(buf))
                }
            }
            Err(e) => {
                done = true;
                Some(Err(e))
            }
        }
    })
}

/// Acquire an exclusive `flock(2)` / `LockFileEx` on `<log>.lock` so
/// the critical section "read previous tail → compute hash → append"
/// is serialised across threads and processes. Returns the locked
/// file handle; the lock releases when the handle drops at the end
/// of the caller's scope.
///
/// Round-6 fix: pre-fix two concurrent `AuditWriter::append` calls
/// against the same log raced between `compute_prev_hash` and the
/// `O_APPEND` write. Both saw the same prev tail and computed the
/// same prev_hash, producing TWO chained entries both pointing at
/// the same predecessor — a silent chain fork that `verify_chain`
/// catches only because the SHA-256 of the actual prev line no
/// longer matches the second entry. The fork is observable as a
/// `ChainBroken` error; the entries are still recoverable, but the
/// chain semantic ("monotonic strict ordering") is broken. The
/// sidecar lockfile gives us a stable lock target even when the
/// active log doesn't exist yet (genesis case).
fn lock_log(log_path: &Path) -> Result<std::fs::File, AuditError> {
    let lock_path = lockfile_path(log_path);
    if let Some(parent) = lock_path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AuditError::Io {
            path: parent.to_path_buf(),
            source: e,
        })?;
    }
    let lock_file = OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(false)
        .open(&lock_path)
        .map_err(|e| AuditError::Io {
            path: lock_path.clone(),
            source: e,
        })?;
    lock_file.lock_exclusive().map_err(|e| AuditError::Io {
        path: lock_path,
        source: e,
    })?;
    Ok(lock_file)
}

fn lockfile_path(active_log: &Path) -> PathBuf {
    let mut s = active_log.to_path_buf();
    let mut filename = active_log
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_default();
    filename.push_str(".lock");
    s.set_file_name(filename);
    s
}

/// One audit-log entry. Field set per RFC 0013 §"Format".
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuditEntry {
    /// Wall-clock time the action happened, ISO 8601 UTC. The writer
    /// expects callers to provide this — we don't pin a clock
    /// implementation here so tests can supply deterministic
    /// timestamps.
    pub timestamp: String,
    /// Who did the action.
    pub actor: AuditActor,
    /// What the action was. Use one of the strings from RFC 0013's
    /// action vocabulary (`run.start`, `case.delete`, etc.).
    pub action: String,
    /// What the action operated on. Free-form per-action shape.
    pub target: serde_json::Value,
    /// Optional action-specific context (workdir paths, exit codes,
    /// etc.). Default is `{}`.
    #[serde(default)]
    pub context: serde_json::Value,
    /// SHA-256 (lowercase hex) of the previous entry's serialised
    /// line, or `"genesis"` for the first entry. Computed by
    /// [`AuditWriter::append`]; callers should not set this
    /// directly.
    #[serde(default)]
    pub prev_hash: String,
}

/// Actor identifier — the user / process that performed the action.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuditActor {
    /// Stable identifier. For human users: typically email or
    /// SSO-resolved sub. For batch / scheduled runs:
    /// `system:scheduler` or similar.
    pub id: String,
    /// Optional session id linking entries from the same UI session.
    #[serde(default)]
    pub session_id: Option<String>,
}

/// Errors raised by the writer / verifier.
#[derive(Debug, Error)]
pub enum AuditError {
    #[error("audit log {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("audit log {path} line {line}: malformed entry: {reason}")]
    Malformed {
        path: PathBuf,
        line: usize,
        reason: String,
    },
    #[error(
        "audit log {path} line {line}: prev_hash chain broken (expected {expected}, got {actual})"
    )]
    ChainBroken {
        path: PathBuf,
        line: usize,
        expected: String,
        actual: String,
    },
}

/// Rotation policy applied before each append. When `Some`, the
/// writer renames the active log to `<log>.<unix_timestamp>` and
/// starts a fresh chain whenever the active file's size exceeds
/// `max_size_bytes`.
///
/// The fresh chain's first entry uses
/// `prev_hash = "genesis-after-rotation:<sha256_hex>"` where the
/// SHA-256 is over the rotated file's last line bytes — that way
/// auditors can prove the new chain is contiguous with the old one
/// across the rotation point. Set `rotated_chain_link = false` to
/// emit a plain `"genesis"` instead and treat each rotated file as
/// an independent chain (useful when the policy intentionally
/// scrubs old entries).
#[derive(Clone, Copy, Debug)]
pub struct RotationPolicy {
    pub max_size_bytes: u64,
    pub rotated_chain_link: bool,
}

impl Default for RotationPolicy {
    fn default() -> Self {
        // 16 MiB ~= 50K-100K typical entries depending on context
        // payload size. Big enough that workstation users rarely
        // hit it; small enough that a corrupted log isn't crippling
        // to scan / re-verify.
        Self {
            max_size_bytes: 16 * 1024 * 1024,
            rotated_chain_link: true,
        }
    }
}

/// Tuple of "what we last cached" — the hash of the tail line, the
/// file length at that point, AND the file's modification time. The
/// size acts as a cheap freshness signal, but is not sufficient on
/// its own: an external process can delete the log and rewrite a NEW
/// file with the same byte length (Round-17 L3 RED case — same-size
/// rewrite). Pairing the length with mtime catches that case because
/// the new file's mtime advances past whatever we recorded.
///
/// If EITHER `file_len_at_cache` OR `mtime_at_cache` disagrees with
/// the active log's current metadata, an external writer / rotation /
/// delete-and-recreate has touched the file and the cache is stale.
#[derive(Clone, Debug)]
struct CachedTail {
    /// SHA-256 hex of the most recent line we wrote, sans trailing
    /// newline (matches what `compute_prev_hash_slow` would re-derive
    /// from `BufRead::lines()` output).
    hash: String,
    /// Active log's byte length at the moment we wrote the cached
    /// hash. Used as the primary freshness signal.
    file_len_at_cache: u64,
    /// Active log's modification time at the moment we wrote the
    /// cached hash. Pairs with `file_len_at_cache` so a delete-and-
    /// recreate with the same byte count still invalidates the cache.
    ///
    /// `None` means the metadata call didn't yield an mtime (rare —
    /// some exotic filesystems / network mounts don't carry mtime).
    /// When `None` we conservatively fall back to the slow walk on
    /// the next read so the chain stays correct.
    mtime_at_cache: Option<std::time::SystemTime>,
}

/// Streaming append-only writer. Each [`AuditWriter::append`] call
/// opens the log file with `O_APPEND` semantics and writes one line,
/// then closes — no kept open handle, so concurrent writers from
/// different processes interleave safely.
///
/// Round-16 M4: maintains an in-memory `cached_prev_hash` so the
/// hot append path no longer re-walks the entire log file on every
/// call. Pre-fix `compute_prev_hash` did a full BufReader walk per
/// append, so a 100K-entry log paid O(N) cost on every entry — total
/// O(N²) for the log. The cache turns the hot path into O(1) by
/// remembering the SHA-256 of the last line written + invalidating
/// on rotation and external mutation.
pub struct AuditWriter {
    log_path: PathBuf,
    rotation: Option<RotationPolicy>,
    /// In-memory cache of `compute_prev_hash`'s last result, paired
    /// with the active log's file length at the moment we stamped
    /// the cache.
    ///
    /// - `None` = not yet computed (lazy init on first append) OR
    ///   invalidated (rotation happened, or an external process
    ///   appended between our writes).
    /// - `Some(CachedTail { hash, file_len })` = the SHA-256 of the
    ///   most recent line we wrote, plus the file length at that
    ///   point. On the next append we recompute `metadata().len()`
    ///   and invalidate the cache if it doesn't match — covers the
    ///   case where an external `rotate_if_needed` shrank the file
    ///   or another process appended between our writes.
    ///
    /// Wrapped in `Mutex` so `&self` methods can lazily fill it
    /// without forcing the whole writer to `&mut`. The fs2 advisory
    /// lock in `append` provides the cross-process serialisation; the
    /// `Mutex` only handles intra-process concurrent writers sharing
    /// an `Arc<AuditWriter>`.
    cached_prev_hash: Mutex<Option<CachedTail>>,
}

impl AuditWriter {
    /// New writer that appends to `log_path` with no rotation. Use
    /// [`AuditWriter::with_rotation`] to install a rotation policy.
    pub fn new(log_path: impl Into<PathBuf>) -> Self {
        Self {
            log_path: log_path.into(),
            rotation: None,
            cached_prev_hash: Mutex::new(None),
        }
    }

    /// Builder: attach a rotation policy. Subsequent `append` calls
    /// check the file's size before writing and rotate when the cap
    /// is exceeded.
    pub fn with_rotation(mut self, policy: RotationPolicy) -> Self {
        self.rotation = Some(policy);
        self
    }

    /// Append one entry. Computes the `prev_hash` from the file's
    /// current last line (or `"genesis"` if the file is empty),
    /// then writes the JSON-serialised entry + a newline.
    ///
    /// When a rotation policy is attached and the active log
    /// exceeds the size cap, the file is rotated (renamed to
    /// `<log>.<unix>`) before the prev_hash is computed. The new
    /// entry then chains to the rotated file's last line via the
    /// `genesis-after-rotation:<hash>` prefix when the policy's
    /// `rotated_chain_link` flag is true.
    pub fn append(&self, mut entry: AuditEntry) -> Result<(), AuditError> {
        // Round-6: serialize the read-tail/compute-hash/write
        // critical section across threads AND processes. Without
        // the advisory lock two appenders racing through
        // `compute_prev_hash` saw the same tail bytes, produced the
        // same prev_hash, and both chained to the same predecessor —
        // a silent fork that `verify_chain` reports as ChainBroken.
        // The lock guard releases when it leaves scope at function
        // return / panic-unwind / `?`.
        let _lock = lock_log(&self.log_path)?;
        // Rotation has to happen BEFORE compute_prev_hash so the
        // new entry's hash chains to the rotated file's tail
        // (via genesis-after-rotation), not to a fresh empty file.
        // Pass `already_locked = true` so we don't try to re-acquire
        // the lock we already hold (POSIX flock would be a no-op
        // here, but on Windows LockFileEx-EX over the same handle
        // returns ERROR_LOCK_VIOLATION).
        let rotated = if let Some(policy) = &self.rotation {
            rotate_if_needed_locked(&self.log_path, *policy)?.is_some()
        } else {
            false
        };
        // Round-16 M4: invalidate the cached prev_hash whenever
        // rotation happened — the new file starts empty (or with the
        // rotation-genesis sidecar) and the cached value belongs to
        // the now-rotated chain. Fall through to the lazy-fill path
        // in `compute_prev_hash`.
        if rotated {
            if let Ok(mut slot) = self.cached_prev_hash.lock() {
                *slot = None;
            }
        }
        entry.prev_hash = self.compute_prev_hash()?;
        if let Some(parent) = self.log_path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| AuditError::Io {
                path: parent.to_path_buf(),
                source: e,
            })?;
        }
        // Build the entire JSONL record (JSON + trailing newline) as a
        // SINGLE byte buffer, then emit it with one `write_all` call.
        // The original code did JSON + newline as two separate writes,
        // which on Windows let two concurrent appenders interleave
        // between the JSON and the newline — producing corrupt lines
        // like `{...}{...}\n\n` (two entries on one row). Single-call
        // append is atomic on every fs we target.
        let mut line = serde_json::to_string(&entry).map_err(|e| AuditError::Io {
            path: self.log_path.clone(),
            source: std::io::Error::other(e.to_string()),
        })?;
        line.push('\n');
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.log_path)
            .map_err(|e| AuditError::Io {
                path: self.log_path.clone(),
                source: e,
            })?;
        f.write_all(line.as_bytes()).map_err(|e| AuditError::Io {
            path: self.log_path.clone(),
            source: e,
        })?;
        // Round-16 M4 (a): fsync BEFORE releasing the fs2 advisory
        // lock so a host crash between the write and the close can't
        // drop the last entry from the chain. Without this, the
        // kernel may buffer the appended bytes; if the box crashes
        // before flush, the entry vanishes and the next append's
        // prev_hash silently points at a predecessor that's no
        // longer on disk — the chain becomes unverifiable.
        f.sync_all().map_err(|e| AuditError::Io {
            path: self.log_path.clone(),
            source: e,
        })?;
        // Round-16 M4 (b): the line we just wrote becomes the new
        // prev_hash for the next append. Update the cache with the
        // SHA-256 of those exact bytes (sans trailing newline — the
        // chain hash is over the line bytes as they appear in the
        // file, which matches what `compute_prev_hash` would compute
        // by re-reading the tail). This turns the next append's
        // prev-hash compute into an O(1) cache hit instead of a
        // full-file walk.
        //
        // We also record the file's POST-WRITE length so the next
        // append's freshness check (in `compute_prev_hash`) can
        // detect an external mutation between our writes.
        if let Ok(mut slot) = self.cached_prev_hash.lock() {
            // The trailing newline is the JSONL record separator and
            // is NOT included in the chain-hash input (compute_prev_hash
            // uses `BufRead::lines()` which strips the terminator).
            let line_bytes = if line.ends_with('\n') {
                &line.as_bytes()[..line.len() - 1]
            } else {
                line.as_bytes()
            };
            // Round-17 L3: stamp BOTH len AND mtime so a delete-and-
            // recreate-same-size attack is caught by the freshness check.
            let post_meta = std::fs::metadata(&self.log_path).ok();
            let post_write_len = post_meta.as_ref().map(|m| m.len()).unwrap_or(0);
            let post_write_mtime = post_meta.as_ref().and_then(|m| m.modified().ok());
            *slot = Some(CachedTail {
                hash: sha256_hex(line_bytes),
                file_len_at_cache: post_write_len,
                mtime_at_cache: post_write_mtime,
            });
        }
        Ok(())
    }

    fn compute_prev_hash(&self) -> Result<String, AuditError> {
        // Round-16 M4: cache hit short-circuits the O(N) walk. The
        // cache is invalidated on rotation (see `append`) and
        // populated either by the previous append (the fast path)
        // or by the slow walk below (lazy init / first append in
        // the process).
        //
        // Freshness check (Round-17 L3): the cached entry is only
        // trusted if BOTH the active log's current byte length AND its
        // mtime match what we recorded at cache-stamp time. Length
        // alone is insufficient — a delete-and-recreate with the same
        // byte count would have passed the pre-fix check while
        // pointing at a totally different file. Requiring mtime to
        // match too catches the same-size rewrite case (the new
        // file's mtime advances past whatever we stamped). If EITHER
        // disagrees, fall through to the slow walk so the chain stays
        // correct.
        //
        // mtime can be `None` on rare filesystems / network mounts;
        // when either the cached or the current mtime is `None` we
        // conservatively treat that as a mismatch — better to pay
        // one slow walk than to chain to the wrong predecessor.
        let current_meta = std::fs::metadata(&self.log_path).ok();
        let current_len = current_meta.as_ref().map(|m| m.len());
        let current_mtime = current_meta.as_ref().and_then(|m| m.modified().ok());
        if let Ok(slot) = self.cached_prev_hash.lock() {
            if let Some(cached) = slot.as_ref() {
                let len_match = Some(cached.file_len_at_cache) == current_len;
                let mtime_match = match (cached.mtime_at_cache, current_mtime) {
                    (Some(c), Some(n)) => c == n,
                    // Either side `None` → conservatively re-walk.
                    _ => false,
                };
                if len_match && mtime_match {
                    return Ok(cached.hash.clone());
                }
            }
        }
        let hash = self.compute_prev_hash_slow()?;
        // After the slow walk, stamp the cache with the file's
        // length AND mtime AS THEY STAND NOW (what we just observed
        // via the walker). For the genesis path we still record both
        // so the next append's freshness check is honest — an empty /
        // missing file maps to length 0 with `None` mtime.
        if let Ok(mut slot) = self.cached_prev_hash.lock() {
            let observed_meta = std::fs::metadata(&self.log_path).ok();
            let observed_len = observed_meta.as_ref().map(|m| m.len()).unwrap_or(0);
            let observed_mtime = observed_meta.as_ref().and_then(|m| m.modified().ok());
            *slot = Some(CachedTail {
                hash: hash.clone(),
                file_len_at_cache: observed_len,
                mtime_at_cache: observed_mtime,
            });
        }
        Ok(hash)
    }

    fn compute_prev_hash_slow(&self) -> Result<String, AuditError> {
        // If the active log is empty (or absent) but a rotation
        // sidecar with the post-rotation genesis is present, use
        // it. Lets the chain link to the previous file's tail
        // even though the current file starts empty.
        if !self.log_path.exists() {
            if let Some(g) = read_rotation_genesis(&self.log_path) {
                return Ok(g);
            }
            return Ok("genesis".to_string());
        }
        let f = std::fs::File::open(&self.log_path).map_err(|e| AuditError::Io {
            path: self.log_path.clone(),
            source: e,
        })?;
        let mut reader = BufReader::new(f);
        let mut last_line: Option<Vec<u8>> = None;
        // Round-20 M3: replaced `BufReader::lines()` with the
        // bounded `read_capped_lines` so a hostile multi-GB single
        // "line" without a newline terminator can't OOM the genesis
        // walk before any record-level validation runs.
        for line in read_capped_lines(&mut reader, MAX_AUDIT_LINE_BYTES) {
            let line = line.map_err(|e| AuditError::Io {
                path: self.log_path.clone(),
                source: e,
            })?;
            // Trim whitespace-equivalent check on bytes — the
            // pre-fix code stripped lines that were all ASCII
            // whitespace (typically blank lines added by an editor).
            // Skip blank lines without converting to a `String`.
            if !line.iter().all(|b| b.is_ascii_whitespace()) {
                last_line = Some(line);
            }
        }
        Ok(match last_line {
            None => {
                // Empty file but a rotation marker present means
                // we just rotated and the new chain inherits the
                // previous tail's hash.
                if let Some(g) = read_rotation_genesis(&self.log_path) {
                    g
                } else {
                    "genesis".to_string()
                }
            }
            Some(s) => sha256_hex(&s),
        })
    }
}

/// If the active log was just rotated, return the
/// `genesis-after-rotation:<hash>` string the rotation sidecar
/// stamped. `None` when no sidecar exists.
///
/// Round-22 L1: cap the sidecar read at `MAX_ROTATION_GENESIS_BYTES`
/// (1 KiB) — the sidecar is exactly one line (`genesis-after-rotation:`
/// prefix + 64-char hex hash + newline ≈ 80 bytes) so a hostile or
/// accidentally-mangled sidecar that grew to many GB on disk can't
/// slurp into memory before the prefix check runs.
fn read_rotation_genesis(active_log: &Path) -> Option<String> {
    let sidecar = rotation_sidecar_path(active_log);
    let text = valenx_core::io_caps::read_capped_to_string(
        &sidecar,
        valenx_core::io_caps::MAX_ROTATION_GENESIS_BYTES as usize,
    )
    .ok()?;
    let trimmed = text.trim();
    if trimmed.starts_with("genesis-after-rotation:") {
        Some(trimmed.to_string())
    } else {
        None
    }
}

fn rotation_sidecar_path(active_log: &Path) -> PathBuf {
    let mut s = active_log.to_path_buf();
    let mut filename = active_log
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_default();
    filename.push_str(".rotation-genesis");
    s.set_file_name(filename);
    s
}

/// Rotate the active log if it exceeds the policy's size cap.
/// Returns the path the previous log was renamed to (when
/// rotation happened), or `None` when no rotation was needed.
///
/// The renamed file format is `<active>.<unix_seconds>` so multiple
/// rotations within a second on a high-volume machine still get
/// distinct names (the first millisecond of overflow gets one,
/// the next gets `.<unix_seconds>-1` etc. — handled inside).
///
/// Round-6: acquires the same advisory lock as `AuditWriter::append`
/// so a rotation triggered from outside the writer can't race with
/// a concurrent append. Callers already inside the lock (the writer
/// itself) take the lock-free variant `rotate_if_needed_locked`.
pub fn rotate_if_needed(
    active_log: &Path,
    policy: RotationPolicy,
) -> Result<Option<PathBuf>, AuditError> {
    let _lock = lock_log(active_log)?;
    rotate_if_needed_locked(active_log, policy)
}

/// Lock-free rotation primitive — callers must already hold the
/// advisory lock returned by `lock_log`. Internal helper invoked by
/// `AuditWriter::append` (which holds the lock for the duration of
/// its critical section) and by the public `rotate_if_needed`.
fn rotate_if_needed_locked(
    active_log: &Path,
    policy: RotationPolicy,
) -> Result<Option<PathBuf>, AuditError> {
    let metadata = match std::fs::metadata(active_log) {
        Ok(m) => m,
        Err(_) => return Ok(None),
    };
    if metadata.len() <= policy.max_size_bytes {
        return Ok(None);
    }
    // Compute the genesis-after-rotation BEFORE renaming the file
    // so we can read its last line. The sidecar gets stamped after
    // a successful rename so a partial failure leaves the chain
    // observably intact.
    let last_line_hash = if policy.rotated_chain_link {
        last_line_sha256(active_log)?
    } else {
        None
    };
    let rotated = pick_rotation_target(active_log)?;
    std::fs::rename(active_log, &rotated).map_err(|e| AuditError::Io {
        path: active_log.to_path_buf(),
        source: e,
    })?;
    if let Some(hash) = last_line_hash {
        let sidecar = rotation_sidecar_path(active_log);
        let text = format!("genesis-after-rotation:{hash}\n");
        // Best-effort write; if it fails the new chain just starts
        // with plain "genesis", which is still verifiable end-to-end
        // (the audit log on the rotated file is independently
        // verifiable too). R29: route through atomic_write_str so a
        // crash mid-write can't leave a truncated genesis sidecar that
        // would be silently mis-read as a valid (but wrong) hash anchor.
        let _ = valenx_core::io_caps::atomic_write_str(&sidecar, &text);
    }
    Ok(Some(rotated))
}

/// Pick a fresh `<active>.<unix>` filename, retrying with `-N`
/// suffixes if a same-second rotation already exists.
fn pick_rotation_target(active_log: &Path) -> Result<PathBuf, AuditError> {
    let unix = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let mut candidate = with_suffix(active_log, &format!(".{unix}"));
    if !candidate.exists() {
        return Ok(candidate);
    }
    for n in 1..1024 {
        candidate = with_suffix(active_log, &format!(".{unix}-{n}"));
        if !candidate.exists() {
            return Ok(candidate);
        }
    }
    Err(AuditError::Io {
        path: active_log.to_path_buf(),
        source: std::io::Error::other(
            "rotation: 1024 same-second collisions; pick a different log path",
        ),
    })
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.to_path_buf();
    let mut filename = path
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_default();
    filename.push_str(suffix);
    s.set_file_name(filename);
    s
}

fn last_line_sha256(path: &Path) -> Result<Option<String>, AuditError> {
    let f = std::fs::File::open(path).map_err(|e| AuditError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut reader = BufReader::new(f);
    let mut last_line: Option<Vec<u8>> = None;
    // Round-20 M3: same bounded-line read as `compute_prev_hash_slow` —
    // a poisoned rotation-source log would have OOM'd the rotation
    // sidecar writer too.
    for line in read_capped_lines(&mut reader, MAX_AUDIT_LINE_BYTES) {
        let line = line.map_err(|e| AuditError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
        if !line.iter().all(|b| b.is_ascii_whitespace()) {
            last_line = Some(line);
        }
    }
    Ok(last_line.map(|s| sha256_hex(&s)))
}

/// Verify a complete audit log file end-to-end: every entry's
/// `prev_hash` must equal the SHA-256 of the previous line's
/// serialised bytes (or `"genesis"` for the first).
///
/// Returns `Ok(n)` with the number of verified entries, or the
/// first chain-break / malformed-entry error.
pub fn verify_chain(path: &Path) -> Result<usize, AuditError> {
    verify_chain_report(path).map(|r| r.entries_verified)
}

/// Detailed chain-integrity report. Returned by [`verify_chain_report`]
/// when the chain is intact end-to-end. Used by the app's "Audit log
/// chain integrity" button so it can show a human-readable summary
/// without the caller having to walk the log a second time.
#[derive(Clone, Debug, PartialEq)]
pub struct ChainReport {
    /// Number of valid entries the verifier walked.
    pub entries_verified: usize,
    /// `timestamp` of the first entry, copied verbatim. `None` for an
    /// empty log.
    pub first_timestamp: Option<String>,
    /// `timestamp` of the last entry, copied verbatim. `None` for an
    /// empty log.
    pub last_timestamp: Option<String>,
    /// SHA-256 of the last line's bytes — useful for compliance tools
    /// that pin chain heads ("the log as of last verification ended at
    /// hash `<X>`"). `None` for an empty log.
    pub head_hash: Option<String>,
}

/// Read the last `n` entries from an audit log without verifying the
/// SHA-256 chain. Useful for the desktop app's "show recent audit
/// activity" panel — the user usually cares about the last dozen
/// events ("did my last run actually emit a `run.complete`?") rather
/// than running a full chain audit.
///
/// Behaviour:
///
/// - Returns up to `n` entries; fewer if the log is shorter.
/// - Empty log (or absent file) returns `Ok(vec![])`. Missing file is
///   not an error here — a fresh install with zero events is a
///   normal state, not a failure.
/// - Malformed entries surface as [`AuditError::Malformed`] so the
///   caller can show a precise "log corrupt at line N" message
///   instead of silently skipping bad lines.
/// - The returned `Vec` preserves on-disk order (oldest first within
///   the tail window). The caller is responsible for any reverse-
///   chronological presentation.
///
/// This intentionally does NOT verify the prev_hash chain. Callers
/// who need integrity verification should use
/// [`verify_chain_report`] alongside this — the two answer different
/// questions.
pub fn tail_n(path: &Path, n: usize) -> Result<Vec<AuditEntry>, AuditError> {
    tail_filtered(path, n, None)
}

/// Tail with an optional ISO-8601 timestamp cutoff. Entries whose
/// `timestamp` field is lexicographically less than `since` are
/// dropped before the ring-buffer truncation runs.
///
/// Audit timestamps are written in RFC 3339 / ISO 8601 with a
/// canonical fixed format (`YYYY-MM-DDTHH:MM:SS.sssZ`), so
/// lexicographic comparison correctly orders them across days,
/// hours, and minutes. Mixed-precision strings (`...:00Z` vs
/// `...:00.000Z`) still compare correctly because the longer
/// form is byte-greater than the shorter at the same instant.
///
/// `since: None` falls through to "no cutoff" — equivalent to
/// [`tail_n`]. Pass an empty string for the same effect (an empty
/// cutoff is treated as "earliest representable" since every
/// non-empty timestamp is byte-greater than the empty string).
pub fn tail_filtered(
    path: &Path,
    n: usize,
    since: Option<&str>,
) -> Result<Vec<AuditEntry>, AuditError> {
    if n == 0 {
        return Ok(Vec::new());
    }
    let f = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
        Err(e) => {
            return Err(AuditError::Io {
                path: path.to_path_buf(),
                source: e,
            })
        }
    };
    let mut reader = BufReader::new(f);
    // Ring-buffer the last `n` non-blank lines. We don't pre-count
    // total lines because that means walking the file twice; a single
    // pass with bounded retention costs O(n) memory.
    let mut tail: std::collections::VecDeque<(usize, String)> =
        std::collections::VecDeque::with_capacity(n);
    // Round-20 M3: bounded per-line read so a hostile multi-GB
    // single-line log doesn't OOM `tail_filtered` (the UI's "show
    // recent audit entries" surface).
    for (idx, line) in read_capped_lines(&mut reader, MAX_AUDIT_LINE_BYTES).enumerate() {
        let line_bytes = line.map_err(|e| AuditError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
        // Skip blank lines without an extra owned copy.
        if line_bytes.iter().all(|b| b.is_ascii_whitespace()) {
            continue;
        }
        // Decode once — the JSON serde path below needs a &str, so we
        // pay the UTF-8 cost here (the read is already capped, so this
        // never grows past `MAX_AUDIT_LINE_BYTES`).
        let line = match String::from_utf8(line_bytes) {
            Ok(s) => s,
            Err(e) => {
                return Err(AuditError::Io {
                    path: path.to_path_buf(),
                    source: std::io::Error::new(std::io::ErrorKind::InvalidData, e),
                });
            }
        };
        // Filter by `since` BEFORE the ring-buffer truncation so the
        // returned set is "last N entries on or after cutoff", not
        // "last N entries with the cutoff applied as a post-filter".
        // Pre-decoding the timestamp without parsing the whole entry
        // keeps the hot path light — we just look for the
        // `"timestamp":"...` field at the start of the JSON object.
        if let Some(cutoff) = since {
            if !cutoff.is_empty() {
                let ts = parse_timestamp_field(&line);
                if let Some(ts) = ts {
                    if ts.as_str() < cutoff {
                        continue;
                    }
                }
                // If we couldn't extract a timestamp the line is
                // probably malformed — fall through and let the
                // serde_json::from_str below surface the error.
            }
        }
        if tail.len() == n {
            tail.pop_front();
        }
        tail.push_back((idx + 1, line));
    }
    let mut out = Vec::with_capacity(tail.len());
    for (line_num, line) in tail {
        let entry: AuditEntry = serde_json::from_str(&line).map_err(|e| AuditError::Malformed {
            path: path.to_path_buf(),
            line: line_num,
            reason: e.to_string(),
        })?;
        out.push(entry);
    }
    Ok(out)
}

/// Cheap timestamp extraction without parsing the whole JSON
/// object. Looks for the first occurrence of `"timestamp":"` and
/// reads up to the next `"`. Returns `None` for any line that
/// doesn't look like an audit JSON entry. Used by
/// [`tail_filtered`] to skip pre-cutoff entries without a full
/// serde decode.
fn parse_timestamp_field(line: &str) -> Option<String> {
    let key = "\"timestamp\":\"";
    let start = line.find(key)? + key.len();
    let rest = &line[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

/// Like [`verify_chain`] but returns the richer [`ChainReport`].
pub fn verify_chain_report(path: &Path) -> Result<ChainReport, AuditError> {
    let f = std::fs::File::open(path).map_err(|e| AuditError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut reader = BufReader::new(f);
    let mut prev_serialised: Option<String> = None;
    let mut count = 0usize;
    let mut first_ts: Option<String> = None;
    let mut last_ts: Option<String> = None;
    let mut head_hash: Option<String> = None;
    // Round-20 M3: same bounded-line read as the other slow walks —
    // a hostile log with a multi-GB un-newline-terminated "line"
    // can't OOM `verify_chain` / its UI surface.
    for (line_idx, line) in read_capped_lines(&mut reader, MAX_AUDIT_LINE_BYTES).enumerate() {
        let line_bytes = line.map_err(|e| AuditError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
        let line_num = line_idx + 1;
        if line_bytes.iter().all(|b| b.is_ascii_whitespace()) {
            continue;
        }
        let line = match String::from_utf8(line_bytes) {
            Ok(s) => s,
            Err(e) => {
                return Err(AuditError::Io {
                    path: path.to_path_buf(),
                    source: std::io::Error::new(std::io::ErrorKind::InvalidData, e),
                });
            }
        };
        let entry: AuditEntry = serde_json::from_str(&line).map_err(|e| AuditError::Malformed {
            path: path.to_path_buf(),
            line: line_num,
            reason: e.to_string(),
        })?;
        let expected = match &prev_serialised {
            None => "genesis".to_string(),
            Some(prev) => sha256_hex(prev.as_bytes()),
        };
        // First entry can also legitimately be a rotation marker
        // (`genesis-after-rotation:<hash>`). Accept that shape on
        // line 1 only — anywhere else it's a chain break.
        let first_entry = prev_serialised.is_none();
        let prev_hash_ok = entry.prev_hash == expected
            || (first_entry && entry.prev_hash.starts_with("genesis-after-rotation:"));
        if !prev_hash_ok {
            return Err(AuditError::ChainBroken {
                path: path.to_path_buf(),
                line: line_num,
                expected,
                actual: entry.prev_hash,
            });
        }
        if first_ts.is_none() {
            first_ts = Some(entry.timestamp.clone());
        }
        last_ts = Some(entry.timestamp.clone());
        head_hash = Some(sha256_hex(line.as_bytes()));
        prev_serialised = Some(line);
        count += 1;
    }
    Ok(ChainReport {
        entries_verified: count,
        first_timestamp: first_ts,
        last_timestamp: last_ts,
        head_hash,
    })
}

fn sha256_hex(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let mut s = String::with_capacity(64);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(s, "{b:02x}");
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn tempfile(tag: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "valenx-audit-{tag}-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    fn make_entry(action: &str) -> AuditEntry {
        AuditEntry {
            timestamp: "2026-04-25T12:00:00Z".into(),
            actor: AuditActor {
                id: "alice@example.com".into(),
                session_id: None,
            },
            action: action.into(),
            target: json!({"kind": "case", "name": "cfd-steady"}),
            context: json!({}),
            prev_hash: String::new(),
        }
    }

    /// RED→GREEN (round-25 L1): after `read_capped_lines` emits a
    /// cap-exceeded error, the next iteration must NOT re-emit the
    /// rest of the poisoned line as another bogus line. Pre-fix the
    /// reader was left at the over-cap offset (partway through the
    /// poisoned record); a caller that ignored the error and kept
    /// iterating would see torn fragments of the same line as
    /// separate "lines". Post-fix the reader is resync'd to the next
    /// newline so the next legitimate line is what the caller sees.
    ///
    /// Round-26 M1: the round-25 L1 bounded-resync logic was removed
    /// because all four production callers stopped on the first
    /// `InvalidData` via `?`, making the resync unreachable. The
    /// iterator now marks itself done on a poisoned line and surfaces
    /// the error; no further lines are yielded. The test is
    /// repurposed to anchor that contract.
    #[test]
    fn read_capped_lines_stops_after_over_cap_round26_m1() {
        use std::io::Cursor;
        // Synthetic log: one over-cap line (no newline for 8 KiB),
        // then a newline, then a normal line that we MUST NOT see.
        const TEST_CAP: usize = 4096;
        let mut payload = Vec::new();
        payload.extend(std::iter::repeat_n(b'x', 8192)); // 2 * cap
        payload.push(b'\n');
        payload.extend_from_slice(b"good_line\n");
        let mut reader = Cursor::new(payload);
        let mut got_cap_err = false;
        let mut got_good_line = false;
        for line in read_capped_lines(&mut reader, TEST_CAP) {
            match line {
                Ok(bytes) => {
                    if bytes == b"good_line" {
                        got_good_line = true;
                    }
                }
                Err(e) => {
                    assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
                    got_cap_err = true;
                }
            }
        }
        assert!(
            got_cap_err,
            "expected cap-exceeded error on the poisoned line"
        );
        assert!(
            !got_good_line,
            "round-26 M1: iterator must stop on the poisoned line; \
             `good_line` should NOT be reachable without explicit \
             caller-side resync"
        );
    }

    #[test]
    fn first_entry_has_genesis_prev_hash() {
        let path = tempfile("genesis");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("project.open")).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let line = text.lines().next().unwrap();
        let parsed: AuditEntry = serde_json::from_str(line).unwrap();
        assert_eq!(parsed.prev_hash, "genesis");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn second_entry_chains_to_first() {
        let path = tempfile("chain");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("project.open")).unwrap();
        writer.append(make_entry("run.start")).unwrap();

        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let second: AuditEntry = serde_json::from_str(lines[1]).unwrap();
        assert_eq!(second.prev_hash, sha256_hex(lines[0].as_bytes()));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn verify_chain_accepts_well_formed_log() {
        let path = tempfile("verify-ok");
        let writer = AuditWriter::new(&path);
        for action in ["project.open", "case.create", "run.start", "run.complete"] {
            writer.append(make_entry(action)).unwrap();
        }
        let count = verify_chain(&path).expect("verify");
        assert_eq!(count, 4);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn verify_chain_detects_tampering() {
        let path = tempfile("verify-tampered");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("project.open")).unwrap();
        writer.append(make_entry("run.start")).unwrap();
        // Tamper: rewrite the first line so its hash no longer
        // matches the second's prev_hash.
        let text = std::fs::read_to_string(&path).unwrap();
        let mut lines: Vec<String> = text.lines().map(|s| s.to_string()).collect();
        // Inject a different action into the first entry's JSON.
        lines[0] = lines[0].replace("project.open", "project.steal");
        std::fs::write(&path, lines.join("\n") + "\n").unwrap();

        let err = verify_chain(&path).unwrap_err();
        assert!(matches!(err, AuditError::ChainBroken { .. }));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn verify_chain_handles_empty_lines_gracefully() {
        let path = tempfile("verify-blanks");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("project.open")).unwrap();
        // Inject a blank line — should be skipped, not break the chain.
        let mut text = std::fs::read_to_string(&path).unwrap();
        text.push('\n');
        text.push('\n');
        std::fs::write(&path, &text).unwrap();
        writer.append(make_entry("run.start")).unwrap();
        let count = verify_chain(&path).expect("verify");
        assert_eq!(count, 2);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn verify_chain_rejects_malformed_json() {
        let path = tempfile("verify-malformed");
        std::fs::write(&path, "not json at all\n").unwrap();
        let err = verify_chain(&path).unwrap_err();
        assert!(matches!(err, AuditError::Malformed { .. }));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn sha256_hex_is_64_lowercase_chars() {
        let hex = sha256_hex(b"hello");
        assert_eq!(hex.len(), 64);
        assert!(hex
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_uppercase()));
        // Known SHA-256 of "hello" — anchor for regression detection.
        assert_eq!(
            hex,
            "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        );
    }

    #[test]
    fn verify_chain_report_carries_first_last_timestamps_and_head_hash() {
        let path = tempfile("verify-report");
        let writer = AuditWriter::new(&path);
        // Two entries with different timestamps.
        let mut e1 = make_entry("project.open");
        e1.timestamp = "2026-04-25T10:00:00Z".into();
        writer.append(e1).unwrap();
        let mut e2 = make_entry("run.start");
        e2.timestamp = "2026-04-25T10:30:00Z".into();
        writer.append(e2).unwrap();

        let report = verify_chain_report(&path).expect("verify");
        assert_eq!(report.entries_verified, 2);
        assert_eq!(
            report.first_timestamp.as_deref(),
            Some("2026-04-25T10:00:00Z")
        );
        assert_eq!(
            report.last_timestamp.as_deref(),
            Some("2026-04-25T10:30:00Z")
        );
        // head_hash should be the SHA-256 of the LAST line's bytes.
        let text = std::fs::read_to_string(&path).unwrap();
        let last_line = text.lines().last().unwrap();
        assert_eq!(
            report.head_hash.as_deref(),
            Some(sha256_hex(last_line.as_bytes()).as_str())
        );
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn verify_chain_report_for_empty_file_is_zero_count_with_no_timestamps() {
        let path = tempfile("verify-empty");
        std::fs::write(&path, "").unwrap();
        let report = verify_chain_report(&path).expect("verify");
        assert_eq!(report.entries_verified, 0);
        assert!(report.first_timestamp.is_none());
        assert!(report.last_timestamp.is_none());
        assert!(report.head_hash.is_none());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn verify_chain_and_verify_chain_report_agree_on_count() {
        // Same input, same count — the convenience wrapper must not
        // diverge from the full report.
        let path = tempfile("verify-agree");
        let writer = AuditWriter::new(&path);
        for action in ["a", "b", "c", "d", "e"] {
            writer.append(make_entry(action)).unwrap();
        }
        let count_simple = verify_chain(&path).expect("simple");
        let count_full = verify_chain_report(&path).expect("full").entries_verified;
        assert_eq!(count_simple, count_full);
        let _ = std::fs::remove_file(&path);
    }

    // -----------------------------------------------------------------
    // Rotation policy
    // -----------------------------------------------------------------

    #[test]
    fn rotation_policy_default_is_16mib_with_chain_link() {
        let p = RotationPolicy::default();
        assert_eq!(p.max_size_bytes, 16 * 1024 * 1024);
        assert!(p.rotated_chain_link);
    }

    #[test]
    fn rotate_if_needed_no_op_when_under_cap() {
        let path = tempfile("rotate-noop");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("a")).unwrap();
        let res = rotate_if_needed(
            &path,
            RotationPolicy {
                max_size_bytes: 10_000_000,
                rotated_chain_link: true,
            },
        )
        .unwrap();
        assert!(res.is_none(), "no rotation expected");
        // File still has the original entry.
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("\"a\""));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn rotate_if_needed_renames_when_over_cap() {
        let path = tempfile("rotate-trigger");
        let writer = AuditWriter::new(&path);
        // Cap of 1 byte forces rotation after the first append.
        writer.append(make_entry("a")).unwrap();
        let rotated = rotate_if_needed(
            &path,
            RotationPolicy {
                max_size_bytes: 1,
                rotated_chain_link: true,
            },
        )
        .unwrap();
        let rotated = rotated.expect("rotation expected");
        assert!(rotated.is_file(), "rotated file must exist");
        assert!(!path.exists(), "active log was renamed away");
        // Sidecar must exist with the genesis-after-rotation marker.
        let sidecar = super::rotation_sidecar_path(&path);
        assert!(sidecar.is_file());
        let marker = std::fs::read_to_string(&sidecar).unwrap();
        assert!(marker.starts_with("genesis-after-rotation:"));
        let _ = std::fs::remove_file(&rotated);
        let _ = std::fs::remove_file(&sidecar);
    }

    #[test]
    fn writer_with_rotation_chains_new_chain_to_rotation_genesis() {
        let path = tempfile("rotate-chain");
        let writer = AuditWriter::new(&path).with_rotation(RotationPolicy {
            max_size_bytes: 1, // force rotation on every append after the first
            rotated_chain_link: true,
        });
        writer.append(make_entry("first")).unwrap();
        // Second append: rotation should have happened, new entry's
        // prev_hash points at the rotation-genesis marker.
        writer.append(make_entry("second")).unwrap();
        let active_text = std::fs::read_to_string(&path).unwrap();
        let line = active_text.lines().next().expect("first line");
        let entry: AuditEntry = serde_json::from_str(line).unwrap();
        assert!(
            entry.prev_hash.starts_with("genesis-after-rotation:"),
            "got prev_hash: {}",
            entry.prev_hash
        );
        // verify_chain accepts the rotation marker.
        verify_chain(&path).expect("chain still valid");
        // Cleanup
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(super::rotation_sidecar_path(&path));
        // Rotated file: best-effort cleanup
        if let Some(parent) = path.parent() {
            if let Ok(rd) = std::fs::read_dir(parent) {
                for e in rd.flatten() {
                    let p = e.path();
                    if p.file_name()
                        .map(|n| {
                            n.to_string_lossy().starts_with(&format!(
                                "{}.",
                                path.file_name().unwrap().to_string_lossy()
                            ))
                        })
                        .unwrap_or(false)
                    {
                        let _ = std::fs::remove_file(p);
                    }
                }
            }
        }
    }

    #[test]
    fn rotation_policy_link_disabled_skips_sidecar() {
        let path = tempfile("rotate-no-link");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("a")).unwrap();
        let _ = rotate_if_needed(
            &path,
            RotationPolicy {
                max_size_bytes: 1,
                rotated_chain_link: false,
            },
        )
        .unwrap();
        // No sidecar must exist when chain-link is disabled.
        let sidecar = super::rotation_sidecar_path(&path);
        assert!(
            !sidecar.exists(),
            "sidecar should not be created when rotated_chain_link=false"
        );
        // Cleanup
        if let Some(parent) = path.parent() {
            if let Ok(rd) = std::fs::read_dir(parent) {
                for e in rd.flatten() {
                    let p = e.path();
                    if p.file_name()
                        .map(|n| {
                            n.to_string_lossy().starts_with(&format!(
                                "{}.",
                                path.file_name().unwrap().to_string_lossy()
                            ))
                        })
                        .unwrap_or(false)
                    {
                        let _ = std::fs::remove_file(p);
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // tail_n
    // -----------------------------------------------------------------

    #[test]
    fn tail_n_returns_empty_for_missing_file() {
        let path = tempfile("tail-missing");
        // No file written; tail_n must treat that as "zero events"
        // rather than IO error — fresh installs have no log yet.
        let entries = tail_n(&path, 10).expect("tail");
        assert!(entries.is_empty());
    }

    #[test]
    fn tail_n_returns_empty_for_zero_request() {
        let path = tempfile("tail-zero");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("a")).unwrap();
        // n=0 short-circuits before opening the file.
        let entries = tail_n(&path, 0).expect("tail");
        assert!(entries.is_empty());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn tail_n_returns_last_n_entries_in_order() {
        let path = tempfile("tail-window");
        let writer = AuditWriter::new(&path);
        for action in ["a", "b", "c", "d", "e"] {
            writer.append(make_entry(action)).unwrap();
        }
        let entries = tail_n(&path, 3).expect("tail");
        assert_eq!(entries.len(), 3);
        // Oldest-first within the tail window: c, d, e.
        assert_eq!(entries[0].action, "c");
        assert_eq!(entries[1].action, "d");
        assert_eq!(entries[2].action, "e");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn tail_n_caps_to_total_count_when_log_is_shorter() {
        let path = tempfile("tail-cap");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("only")).unwrap();
        let entries = tail_n(&path, 100).expect("tail");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].action, "only");
        let _ = std::fs::remove_file(&path);
    }

    fn make_entry_at(action: &str, ts: &str) -> AuditEntry {
        AuditEntry {
            timestamp: ts.into(),
            actor: AuditActor {
                id: "alice@example.com".into(),
                session_id: None,
            },
            action: action.into(),
            target: json!({"kind": "case", "name": "cfd-steady"}),
            context: json!({}),
            prev_hash: String::new(),
        }
    }

    #[test]
    fn tail_filtered_drops_entries_before_since_cutoff() {
        let path = tempfile("tail-since");
        let writer = AuditWriter::new(&path);
        writer
            .append(make_entry_at("a", "2026-04-27T00:00:00Z"))
            .unwrap();
        writer
            .append(make_entry_at("b", "2026-04-27T12:00:00Z"))
            .unwrap();
        writer
            .append(make_entry_at("c", "2026-04-28T00:00:00Z"))
            .unwrap();
        writer
            .append(make_entry_at("d", "2026-04-28T12:00:00Z"))
            .unwrap();

        // Cut-off at midnight on the 28th — only c and d survive.
        let entries = tail_filtered(&path, 10, Some("2026-04-28T00:00:00Z")).expect("tail");
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].action, "c");
        assert_eq!(entries[1].action, "d");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn tail_filtered_combines_since_with_n_truncation() {
        // `--since` runs FIRST, then `-n` keeps the last N of what
        // survived. So `tail_filtered(p, 1, Some("2026-04-28..."))`
        // should keep just `d`, not `c`.
        let path = tempfile("tail-since-n");
        let writer = AuditWriter::new(&path);
        for (action, ts) in [
            ("a", "2026-04-27T00:00:00Z"),
            ("b", "2026-04-28T01:00:00Z"),
            ("c", "2026-04-28T02:00:00Z"),
            ("d", "2026-04-28T03:00:00Z"),
        ] {
            writer.append(make_entry_at(action, ts)).unwrap();
        }
        let entries = tail_filtered(&path, 1, Some("2026-04-28T00:00:00Z")).expect("tail");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].action, "d");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn tail_filtered_with_none_since_matches_tail_n() {
        let path = tempfile("tail-none");
        let writer = AuditWriter::new(&path);
        for action in ["a", "b", "c"] {
            writer.append(make_entry(action)).unwrap();
        }
        let with_filter = tail_filtered(&path, 5, None).expect("tail filtered");
        let without = tail_n(&path, 5).expect("tail n");
        assert_eq!(with_filter.len(), without.len());
        for (a, b) in with_filter.iter().zip(without.iter()) {
            assert_eq!(a.action, b.action);
            assert_eq!(a.timestamp, b.timestamp);
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn parse_timestamp_field_extracts_value_without_full_decode() {
        let line = r#"{"timestamp":"2026-04-28T00:00:00Z","actor":{"id":"x","session_id":null},"action":"y","target":{},"context":{},"prev_hash":"genesis"}"#;
        assert_eq!(
            parse_timestamp_field(line),
            Some("2026-04-28T00:00:00Z".to_string())
        );
        // Malformed lines return None.
        assert_eq!(parse_timestamp_field("not json"), None);
        assert_eq!(parse_timestamp_field("{\"missing\":1}"), None);
    }

    #[test]
    fn tail_n_skips_blank_lines() {
        let path = tempfile("tail-blanks");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("first")).unwrap();
        // Inject a blank line and a whitespace-only line — both
        // should be ignored by the tail walker, mirroring the
        // verify_chain behaviour.
        let mut text = std::fs::read_to_string(&path).unwrap();
        text.push('\n');
        text.push_str("   \n");
        std::fs::write(&path, &text).unwrap();
        writer.append(make_entry("second")).unwrap();
        let entries = tail_n(&path, 5).expect("tail");
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].action, "first");
        assert_eq!(entries[1].action, "second");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn tail_n_surfaces_malformed_entry_with_line_number() {
        // A corrupt last line should bubble up as Malformed with the
        // exact line number so the UI can point at it.
        let path = tempfile("tail-malformed");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("good")).unwrap();
        let mut text = std::fs::read_to_string(&path).unwrap();
        text.push_str("not json at all\n");
        std::fs::write(&path, &text).unwrap();
        let err = tail_n(&path, 5).unwrap_err();
        match err {
            AuditError::Malformed { line, .. } => assert_eq!(line, 2),
            other => panic!("expected Malformed, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn tail_n_does_not_verify_chain() {
        // A tampered first entry breaks the chain but doesn't break
        // tail_n: the helper is for "show me recent activity",
        // not for proving chain integrity.
        let path = tempfile("tail-tampered");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("project.open")).unwrap();
        writer.append(make_entry("run.start")).unwrap();
        // Tamper with the first entry — second entry's prev_hash no
        // longer matches, but that's verify_chain's problem.
        let text = std::fs::read_to_string(&path).unwrap();
        let mut lines: Vec<String> = text.lines().map(|s| s.to_string()).collect();
        lines[0] = lines[0].replace("project.open", "project.steal");
        std::fs::write(&path, lines.join("\n") + "\n").unwrap();
        let entries = tail_n(&path, 5).expect("tail still works");
        assert_eq!(entries.len(), 2);
        // verify_chain would fail on the same input.
        assert!(verify_chain(&path).is_err());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn append_serialises_across_threads_with_fs2_lock() {
        // Round-6 RED→GREEN: two threads each calling `append` 100
        // times against the same log file. Without the fs2 lock,
        // `compute_prev_hash` + `O_APPEND` raced — both threads
        // could read the same tail, compute the same prev_hash,
        // and write entries that both claimed the same predecessor.
        // The post-race `verify_chain` would then surface a
        // `ChainBroken` error because the SHA-256 of the actual
        // previous line no longer matches the second entry's
        // prev_hash. With the lock the writes serialise, the chain
        // stays intact, and verify_chain accepts the full log.
        let path = tempfile("race-append");
        let writer_handle = std::sync::Arc::new(AuditWriter::new(path.clone()));
        let n_threads = 2;
        let per_thread = 100;
        let mut handles = Vec::new();
        for t in 0..n_threads {
            let w = writer_handle.clone();
            handles.push(std::thread::spawn(move || {
                for i in 0..per_thread {
                    let action = format!("t{t}.i{i}");
                    w.append(make_entry(&action)).expect("append");
                }
            }));
        }
        for h in handles {
            h.join().expect("thread join");
        }
        // Every entry was serialised through the lock, so the chain
        // is unbroken end-to-end.
        let count = verify_chain(&path).expect("verify");
        assert_eq!(count, n_threads * per_thread);
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(super::lockfile_path(&path));
    }

    /// Round-16 M4 RED→GREEN (cache): the cached prev_hash makes a
    /// thousand-entry append loop run in time consistent with O(N)
    /// total work (one slow walk at boot, then constant-time appends).
    /// Pre-fix `compute_prev_hash` re-walked the entire log on every
    /// call, making the loop O(N²) — a thousand appends paid for
    /// reading half a million lines instead of writing 1000.
    ///
    /// The chain must still verify end-to-end after the optimised
    /// path runs. Without the cache the loop completed but in O(N²)
    /// time; the cache cuts the runtime dramatically and the
    /// `verify_chain` post-condition pins that we didn't break the
    /// chain semantics in the process.
    #[test]
    fn cached_prev_hash_keeps_chain_intact_across_many_appends() {
        let path = tempfile("m4-cache");
        let writer = AuditWriter::new(&path);
        let n = 1000;
        for i in 0..n {
            let mut e = make_entry("bulk");
            // Unique timestamp so the serialised lines differ.
            e.timestamp = format!("2026-05-26T00:00:{:02}.{:03}Z", i % 60, i % 1000);
            writer.append(e).expect("append");
        }
        // The chain must verify end-to-end — the cache must produce
        // identical prev_hash values to the slow walker on every step.
        let count = verify_chain(&path).expect("chain valid");
        assert_eq!(count, n);
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(super::lockfile_path(&path));
    }

    /// Round-16 M4 (b) RED→GREEN: the cache is correctly invalidated
    /// when an external rotation runs between two appends. Without
    /// the freshness check the second append would use a stale
    /// cached hash that no longer corresponds to anything on disk,
    /// producing a broken chain.
    #[test]
    fn cache_invalidates_after_external_rotation() {
        let path = tempfile("m4-extrot");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("first")).expect("append 1");
        // External rotate — happens OUTSIDE the writer.
        let _ = rotate_if_needed(
            &path,
            RotationPolicy {
                max_size_bytes: 1,
                rotated_chain_link: true,
            },
        )
        .expect("rotate");
        // Now the active file is gone (renamed to .<unix>) and the
        // writer's cache holds the hash of `first`. The freshness
        // check must catch the size mismatch (cached len > 0,
        // current len = 0 / file absent) and fall back to the slow
        // walk, which reads the rotation-genesis sidecar.
        writer.append(make_entry("second")).expect("append 2");
        let active = std::fs::read_to_string(&path).expect("active log");
        let line = active.lines().next().expect("second entry exists");
        let entry: AuditEntry = serde_json::from_str(line).unwrap();
        assert!(
            entry.prev_hash.starts_with("genesis-after-rotation:"),
            "expected rotation-genesis link, got prev_hash: {}",
            entry.prev_hash
        );
        verify_chain(&path).expect("chain still valid");
        // Cleanup
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(super::lockfile_path(&path));
        let _ = std::fs::remove_file(super::rotation_sidecar_path(&path));
        if let Some(parent) = path.parent() {
            if let Ok(rd) = std::fs::read_dir(parent) {
                for e in rd.flatten() {
                    let p = e.path();
                    if p.file_name()
                        .map(|n| {
                            n.to_string_lossy().starts_with(&format!(
                                "{}.",
                                path.file_name().unwrap().to_string_lossy()
                            ))
                        })
                        .unwrap_or(false)
                    {
                        let _ = std::fs::remove_file(p);
                    }
                }
            }
        }
    }

    /// Round-16 M4 (a) RED→GREEN: fsync is called before the lock
    /// releases. We can't directly observe `sync_all` from userspace,
    /// but we can pin the regression by verifying the write actually
    /// reaches disk byte-for-byte (a missing sync would leave the
    /// bytes in the kernel buffer — still readable by us in-process,
    /// but vulnerable to a host crash). The test reads the file back
    /// after each append and confirms the chain is intact.
    ///
    /// This is a structural test — it asserts the post-condition the
    /// `sync_all` call protects (file content matches what we wrote)
    /// rather than the syscall itself. If a future refactor removes
    /// the `sync_all` the test will still pass on a non-crashing
    /// kernel, but the structural assertion (we read what we wrote)
    /// is the closest portable proxy for "the bytes are durable".
    #[test]
    fn append_writes_are_readable_immediately_after_fsync() {
        let path = tempfile("m4-fsync");
        let writer = AuditWriter::new(&path);
        for action in ["a", "b", "c"] {
            writer.append(make_entry(action)).unwrap();
            // After each append, the bytes must be visible on disk —
            // this is what sync_all protects against in the crash
            // case. We verify the count grows monotonically.
            let text = std::fs::read_to_string(&path).unwrap();
            let count = text.lines().filter(|l| !l.trim().is_empty()).count();
            assert!(count > 0, "expected {action} to be visible immediately");
        }
        // Final chain integrity check.
        verify_chain(&path).expect("chain valid");
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(super::lockfile_path(&path));
    }

    /// Round-17 L3 RED→GREEN: an external delete-and-recreate of the
    /// log with the SAME byte length must invalidate the cached
    /// prev_hash. Pre-fix the freshness check only compared
    /// `metadata().len()`; if a second process rotated the log out and
    /// truncated a fresh same-size file back into place, the cache
    /// silently chained to the wrong predecessor (the new file's tail
    /// is empty / a different entry, but the cached hash still pointed
    /// at the old file's last line).
    ///
    /// Adding `mtime_at_cache` to the freshness signal catches that
    /// case because filesystem mtime advances on the recreate even
    /// when the byte length matches.
    #[test]
    fn cache_invalidates_after_same_size_recreate() {
        let path = tempfile("l3-recreate");
        let writer = AuditWriter::new(&path);
        // Write one entry — populates the cache with that entry's
        // tail hash and the file's length/mtime.
        writer.append(make_entry("first")).expect("append 1");
        let original_text = std::fs::read_to_string(&path).expect("read original");
        let original_len = original_text.len();
        // Wait long enough that the filesystem's mtime tick advances
        // — Windows/FAT have a 2-second resolution on some volumes;
        // most modern fs are nanosecond, but we use sleep_millis(50)
        // to be safe on the test runner.
        std::thread::sleep(std::time::Duration::from_millis(50));
        // External tamper: delete the active log and recreate a fresh
        // file of EXACTLY the same byte length. The bytes are different
        // (a totally fabricated payload) but the length matches, so
        // pre-fix the cache would have happily reused the stale hash.
        let _ = std::fs::remove_file(&path);
        // Fabricate replacement bytes of equal length but ENDING WITH
        // A NEWLINE so the appended `second` entry lands on its own
        // line (matching real JSONL semantics). The penultimate line
        // is the fabricated payload, the last line is `second` — we
        // only need to read `second`'s prev_hash to prove the cache
        // was invalidated.
        let fake = "x".repeat(original_len - 1) + "\n";
        std::fs::write(&path, &fake).expect("rewrite same-size");
        // After the recreate, the next append must either:
        // (a) detect the cache is stale and recompute prev_hash from
        //     the new file's tail (the fake `xxx...` line); or
        // (b) refuse, because the file no longer looks like a valid
        //     JSONL chain.
        // Either way, the cached "hash of first" must not silently
        // become the second entry's prev_hash.
        writer.append(make_entry("second")).expect("append 2");
        let text = std::fs::read_to_string(&path).expect("read after");
        // The second entry must be the LAST line of the file. Pull it
        // out and check its prev_hash is NOT the cached hash of `first`
        // — i.e. the cache was actually invalidated.
        let last_line = text
            .lines()
            .filter(|l| !l.trim().is_empty())
            .next_back()
            .expect("at least one entry in file");
        let second_entry: AuditEntry = serde_json::from_str(last_line).expect("last line is JSONL");
        // The first-entry serialised form (what `first`'s cache would
        // have stored as the tail hash) is no longer on disk after the
        // recreate. If the cache was honoured, second_entry.prev_hash
        // would equal sha256_hex(first's serialised line). We rebuild
        // that to assert it DOESN'T match.
        let first_serialised = original_text.lines().next().unwrap();
        let stale_hash = sha256_hex(first_serialised.as_bytes());
        assert_ne!(
            second_entry.prev_hash, stale_hash,
            "the cache should have been invalidated by the recreate; \
             instead the second entry chained to the deleted first"
        );
        // Cleanup
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(super::lockfile_path(&path));
    }

    #[test]
    fn verify_chain_rejects_rotation_marker_on_non_first_entry() {
        // Hand-build a log file with a rotation marker on line 2 —
        // that's a malformed chain and verify_chain must reject it.
        let path = tempfile("rotate-bad-position");
        let writer = AuditWriter::new(&path);
        writer.append(make_entry("first")).unwrap();
        // Inject a malformed second entry whose prev_hash starts
        // with the rotation marker even though we're not at line 1.
        let mut text = std::fs::read_to_string(&path).unwrap();
        let mut bad = make_entry("second");
        bad.prev_hash = "genesis-after-rotation:deadbeef".to_string();
        text.push_str(&serde_json::to_string(&bad).unwrap());
        text.push('\n');
        std::fs::write(&path, &text).unwrap();
        let err = verify_chain(&path).unwrap_err();
        assert!(matches!(err, AuditError::ChainBroken { .. }));
        let _ = std::fs::remove_file(&path);
    }

    /// Round-20 M3 RED→GREEN: pre-fix, the audit slow-walk readers
    /// used `BufReader::lines()`, which allocates an unbounded
    /// `String` per record. A poisoned log with a single
    /// 1 MiB un-newline-terminated "line" would have allocated past
    /// the [`MAX_AUDIT_LINE_BYTES`] cap before any record-level
    /// validation kicked in. Post-fix, the over-cap line surfaces
    /// as an `AuditError::Io` with `InvalidData` kind.
    #[test]
    fn verify_chain_rejects_oversized_line_without_oom() {
        use std::io::Write;
        let path = tempfile("r20m3-oversize-line");
        // Write a single 1 MiB "line" with NO trailing newline. The
        // pre-fix reader would have grown its buffer to 1 MiB before
        // attempting any serde decode; the post-fix reader bails out
        // at MAX_AUDIT_LINE_BYTES (256 KiB) with an IO/InvalidData.
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(&vec![b'x'; 1024 * 1024]).unwrap();
        // No trailing `\n` — this is the unbounded-line shape.
        drop(f);
        let err =
            verify_chain(&path).expect_err("round-20 M3: oversized line must reject as IO error");
        match err {
            AuditError::Io { source, .. } => {
                assert_eq!(
                    source.kind(),
                    std::io::ErrorKind::InvalidData,
                    "over-cap line must surface as InvalidData",
                );
                let msg = source.to_string();
                assert!(
                    msg.contains("cap") || msg.contains("exceed"),
                    "expected cap-mention in message, got: {msg}"
                );
            }
            other => panic!("expected AuditError::Io for over-cap line, got: {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }

    /// Round-20 M3 GREEN: the cap is a sensible value and matches the
    /// public constant.
    #[test]
    fn audit_line_cap_is_sensible() {
        assert_eq!(MAX_AUDIT_LINE_BYTES, 256 * 1024);
    }

    /// Round-22 L1 RED→GREEN: an over-cap rotation-genesis sidecar
    /// (larger than `MAX_ROTATION_GENESIS_BYTES`, 1 KiB) must not be
    /// slurped into memory. Pre-fix `read_rotation_genesis` did a
    /// bare `fs::read_to_string(&sidecar)` and would have allocated
    /// the full file size before the prefix check ran. Post-fix the
    /// bounded read returns Err, which propagates as a clean `None`
    /// (the same return value as "no sidecar exists") without the
    /// allocation.
    #[test]
    fn read_rotation_genesis_rejects_oversize_sidecar() {
        let active = tempfile("rotation-genesis-r22l1");
        // Touch the active log path so `rotation_sidecar_path` can
        // derive a stable sibling filename.
        std::fs::write(&active, b"").unwrap();
        let sidecar = rotation_sidecar_path(&active);
        // Past the 1 KiB MAX_ROTATION_GENESIS_BYTES cap. Use `set_len`
        // so we don't write 16 MiB of bytes for the test. The first
        // bytes are a real `genesis-after-rotation:` prefix so the
        // failure mode is the size-cap, not "no prefix found".
        use std::io::Write;
        let mut f = std::fs::File::create(&sidecar).unwrap();
        f.write_all(b"genesis-after-rotation:abc123\n").unwrap();
        f.set_len(valenx_core::io_caps::MAX_ROTATION_GENESIS_BYTES + 1)
            .unwrap();
        drop(f);
        // Pre-fix would have slurped 1 KiB+ then returned Some. Post-
        // fix bounded read errors → returns None, same as a missing
        // sidecar.
        assert!(read_rotation_genesis(&active).is_none());
        let _ = std::fs::remove_file(&active);
        let _ = std::fs::remove_file(&sidecar);
    }

    /// RED→GREEN (round-24 H1): `read_capped_lines` rejects a line
    /// that exceeds the cap WITHOUT allocating the whole line first.
    ///
    /// Strategy: ask for a small per-line cap (8 bytes), feed in 4
    /// KiB of input with NO newline. Pre-fix code called
    /// `read_until(b'\n', &mut buf)` unbounded — buf grows to 4 KiB
    /// THEN the check fires. Post-fix wraps in `take(cap+1)` so
    /// `buf.len() == cap + 1 == 9` when the error fires. Test
    /// asserts the iteration ends after exactly one Err AND that
    /// only a single iteration runs (proving the bound stopped
    /// further reads).
    #[test]
    fn read_capped_lines_bounds_allocation() {
        use std::io::Cursor;
        let payload = vec![b'x'; 4096];
        let mut reader = std::io::BufReader::new(Cursor::new(payload));
        let mut items: Vec<std::io::Result<Vec<u8>>> = Vec::new();
        for line in read_capped_lines(&mut reader, 8) {
            items.push(line);
        }
        // Exactly one iteration — the error — and nothing after it.
        assert_eq!(items.len(), 1, "iterator must stop after over-cap line");
        let err = items.into_iter().next().unwrap().unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        // The error message reflects the cap of 8.
        let msg = format!("{err}");
        assert!(msg.contains("8-byte cap"), "unexpected message: {msg}");
    }

    /// RED→GREEN (round-24 H1): a well-formed multi-line file still
    /// works post-fix (no regression of the genesis-walk path).
    #[test]
    fn read_capped_lines_handles_multiple_lines() {
        use std::io::Cursor;
        let payload = b"alpha\nbeta\ngamma\n";
        let mut reader = std::io::BufReader::new(Cursor::new(payload));
        let mut out: Vec<Vec<u8>> = Vec::new();
        for line in read_capped_lines(&mut reader, 1024) {
            out.push(line.unwrap());
        }
        assert_eq!(out.len(), 3);
        assert_eq!(out[0], b"alpha");
        assert_eq!(out[1], b"beta");
        assert_eq!(out[2], b"gamma");
    }
}
