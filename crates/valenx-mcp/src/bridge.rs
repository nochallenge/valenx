//! **File-bridge to a *running* valenx GUI.**
//!
//! valenx is AI-drivable through a newline-delimited-JSON *command file* that
//! the live app polls each frame (see `crates/valenx-app/src/agent_commands.rs`):
//! an external agent appends one JSON command per line and valenx executes it
//! through the same vetted tab/dock/workbench methods a user click would. This
//! module is the **writer/reader half** of that bridge, used by the
//! [`crate::bridge_tools`] MCP tools so *any* MCP-capable client can drive a
//! local valenx without re-implementing the wire format.
//!
//! # Two channels
//!
//! - **Global** command file `<base>/valenx_chat_cmd.jsonl` (no suffix). valenx
//!   polls this on every frame regardless of how many units exist, and it is the
//!   **only** channel that honours the `new_unit` command (the bootstrap that
//!   opens a Workbench+Agent unit before any unit exists). See [`global_cmd_path`].
//! - **Per-unit** command file `<base>/valenx_chat_cmd_u{n}.jsonl` for unit `n`.
//!   Every *other* command (open_workbench, set_control, run_command, …) is a
//!   per-unit command and must be written here. See [`unit_cmd_path`].
//!
//! `<base>` is the directory of `$VALENX_ASSISTANT_INBOX` (or the OS temp dir
//! when that env var is unset) — exactly the derivation valenx's
//! `agent_commands::cmd_path` uses, so the file the MCP server writes is the
//! file valenx reads.
//!
//! # Acks / results
//!
//! valenx posts every ack/result as a *feed* line `{"title","detail","kind"}`
//! into the matching feed file — the **base** feed `$VALENX_ASSISTANT_FEED`
//! (or `<state_dir>/assistant_feed.jsonl`, with an OS-temp fallback) for the
//! global channel, and `<base-feed>` with `_u{n}` inserted before `.jsonl` for
//! unit `n`. [`read_feed_since`] reads the lines appended *after* a recorded
//! byte offset so a tool can return just the new ack(s) it triggered. The feed
//! path is resolved independently of the command path (different env var), so
//! the launch convention that points both at the same temp dir keeps the bridge
//! whole.
//!
//! # Security posture — LOCAL ONLY
//!
//! This bridge speaks to a valenx running **on the same machine** purely through
//! **local files** in a user-writable directory. It opens **no socket**, binds
//! **no port**, and makes **no network call**; nothing here is reachable from
//! another host. The MCP transport itself is stdio (a child process of the MCP
//! client), so the whole path — client → MCP server → command file → valenx — is
//! confined to the local user session. An agent can only drive a valenx the
//! local user has already launched and pointed at the same `$VALENX_ASSISTANT_*`
//! directory.

use std::path::{Path, PathBuf};

/// Fixed stem for the command files (mirrors `agent_commands::CMD_STEM`). The
/// global file is `<stem>.jsonl`; a per-unit file is `<stem>_u{n}.jsonl`.
const CMD_STEM: &str = "valenx_chat_cmd";

/// Resolve the command-file **base directory** the live valenx polls: the
/// directory of `$VALENX_ASSISTANT_INBOX` if that env var is set+non-empty,
/// else the OS temp dir. This is the exact derivation
/// `valenx_app::agent_commands::cmd_path` uses, so the file written here is the
/// file valenx reads.
pub fn base_dir() -> PathBuf {
    if let Ok(p) = std::env::var("VALENX_ASSISTANT_INBOX") {
        if !p.is_empty() {
            if let Some(parent) = Path::new(&p).parent() {
                // `parent()` of a bare filename is `""`; treat that as "current
                // dir" by falling through to temp only when there is truly no
                // parent component.
                if !parent.as_os_str().is_empty() {
                    return parent.to_path_buf();
                }
            }
        }
    }
    std::env::temp_dir()
}

/// The **global** command file `<base>/valenx_chat_cmd.jsonl` (no `_u` suffix).
/// The only channel that honours the `new_unit` command.
pub fn global_cmd_path() -> PathBuf {
    base_dir().join(format!("{CMD_STEM}.jsonl"))
}

/// The **per-unit** command file `<base>/valenx_chat_cmd_u{n}.jsonl` for unit
/// `n`. Every non-`new_unit` command is written here.
pub fn unit_cmd_path(n: usize) -> PathBuf {
    base_dir().join(format!("{CMD_STEM}_u{n}.jsonl"))
}

/// Resolve the **base feed file** valenx writes acks/results into:
/// `$VALENX_ASSISTANT_FEED` if set+non-empty, else `<state_dir>/assistant_feed.jsonl`,
/// with an OS-temp fallback. Mirrors `assistant_workbench::assistant_feed_path`.
///
/// `valenx-mcp` deliberately avoids depending on `valenx-app` (which would pull
/// in the GUI), so the `<state_dir>` default is reproduced here with the same
/// `directories`-free, env-driven logic: when `$VALENX_ASSISTANT_FEED` is unset
/// we cannot see valenx's true state dir without the GUI crate, so we fall back
/// to the OS temp dir — which is also where the documented launch convention
/// points `$VALENX_ASSISTANT_FEED`. In practice an agent driving valenx sets
/// both `$VALENX_ASSISTANT_INBOX` and `$VALENX_ASSISTANT_FEED` to the same temp
/// directory, so command and feed files sit side by side.
pub fn base_feed_path() -> PathBuf {
    if let Ok(p) = std::env::var("VALENX_ASSISTANT_FEED") {
        if !p.is_empty() {
            return PathBuf::from(p);
        }
    }
    std::env::temp_dir().join("valenx_assistant_feed.jsonl")
}

/// Insert `_u{n}` before a feed path's final `.jsonl` (mirrors
/// `assistant_workbench::per_unit_path`): `…/valenx_chat_feed.jsonl` →
/// `…/valenx_chat_feed_u3.jsonl`. A path that does not end in `.jsonl` gets the
/// suffix appended whole, keeping the parent directory unchanged.
fn per_unit_feed(base: &Path, n: usize) -> PathBuf {
    let file = base
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let renamed = match file.strip_suffix(".jsonl") {
        Some(stem) => format!("{stem}_u{n}.jsonl"),
        None => format!("{file}_u{n}"),
    };
    match base.parent() {
        Some(dir) if !dir.as_os_str().is_empty() => dir.join(renamed),
        _ => PathBuf::from(renamed),
    }
}

/// The **per-unit feed file** for unit `n` (`_u{n}` before `.jsonl` on
/// [`base_feed_path`]) — where valenx posts that unit's acks/results.
pub fn unit_feed_path(n: usize) -> PathBuf {
    per_unit_feed(&base_feed_path(), n)
}

/// The feed file a command on `channel` produces acks into: the base feed for
/// the global channel (`None`), else unit `n`'s feed.
pub fn feed_path_for(channel: Option<usize>) -> PathBuf {
    match channel {
        None => base_feed_path(),
        Some(n) => unit_feed_path(n),
    }
}

/// **Append one JSON command line** to `path`, creating it if absent. This is
/// the write half of the bridge: valenx polls `path` and runs every appended
/// line from the first (stale files are wiped by valenx at launch, so appending
/// is always "new commands now"). `value` is serialised compactly and a single
/// trailing `\n` is added so each command is its own JSONL line.
///
/// Returns the serialised line (without the newline) so a caller can echo the
/// exact wire text it wrote, and surfaces any I/O error to the MCP layer.
pub fn write_command(path: &Path, value: &serde_json::Value) -> std::io::Result<String> {
    use std::io::Write;
    let line = serde_json::to_string(value)?;
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    writeln!(f, "{line}")?;
    Ok(line)
}

/// The current byte length of `path` (0 if it does not exist yet). Recorded
/// *before* a command is written so [`read_feed_since`] can return only the feed
/// lines valenx appends in response.
pub fn feed_len(path: &Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

/// One parsed feed entry — valenx's ack/result line shape
/// (`assistant_workbench::FeedEntry`): a `title` (the poster, e.g. `"Claude"`),
/// a `detail` (the message body), and an accent `kind`
/// (`build`/`result`/`ship`/`warn`/`user`/…).
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
pub struct FeedEntry {
    /// Who posted the line (valenx's bridge uses `"Claude"`).
    #[serde(default)]
    pub title: String,
    /// The message body — the human-readable ack/result text.
    #[serde(default)]
    pub detail: String,
    /// Accent tag: `build` / `result` / `ship` / `warn` / `user` / other.
    #[serde(default)]
    pub kind: String,
}

/// **Read the feed lines appended after `since` bytes** of `path`, parsed into
/// [`FeedEntry`]s. Blank and malformed lines are skipped (a half-written final
/// line while valenx is appending must not break parsing), exactly like
/// valenx's own `parse_feed`. A missing/unreadable file yields an empty list.
///
/// The read is bounded: at most [`MAX_FEED_READ_BYTES`] are pulled from the tail
/// so a runaway feed can't blow up the MCP server. `since` past the current
/// length (file truncated/rotated) is clamped to read from the start of what
/// remains.
pub fn read_feed_since(path: &Path, since: u64) -> Vec<FeedEntry> {
    use std::io::{Read, Seek, SeekFrom};
    let Ok(mut f) = std::fs::File::open(path) else {
        return Vec::new();
    };
    let len = f.metadata().map(|m| m.len()).unwrap_or(0);
    // Clamp the start: if the file shrank since we recorded `since`, read from
    // the start of the remaining bytes rather than seeking past EOF.
    let start = since.min(len);
    if f.seek(SeekFrom::Start(start)).is_err() {
        return Vec::new();
    }
    let mut buf = Vec::new();
    if f.take(MAX_FEED_READ_BYTES as u64)
        .read_to_end(&mut buf)
        .is_err()
    {
        return Vec::new();
    }
    let body = String::from_utf8_lossy(&buf);
    body.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .filter_map(|l| serde_json::from_str::<FeedEntry>(l).ok())
        .collect()
}

/// Upper bound on how many tail bytes [`read_feed_since`] pulls from a feed file
/// in one call — generous for the handful of ack lines a single tool triggers,
/// while bounding a hostile/runaway feed.
pub const MAX_FEED_READ_BYTES: usize = 1024 * 1024;

/// Render a slice of feed entries as a compact human-readable block for an MCP
/// `content` text reply: one `kind: detail` line each. Empty input → a short
/// "no ack yet" sentinel so a caller always gets *some* feedback even if valenx
/// is not running (or had not polled within the wait window).
pub fn render_feed(entries: &[FeedEntry]) -> String {
    if entries.is_empty() {
        return "(no response from valenx yet — is it running and pointed at the same \
                $VALENX_ASSISTANT_INBOX/$VALENX_ASSISTANT_FEED directory?)"
            .to_string();
    }
    entries
        .iter()
        .map(|e| format!("[{}] {}", e.kind, e.detail))
        .collect::<Vec<_>>()
        .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::env_lock;

    #[test]
    fn cmd_paths_follow_inbox_env_dir() {
        let _g = env_lock();
        let dir = std::env::temp_dir().join("valenx-mcp-bridge-pathtest");
        std::env::set_var(
            "VALENX_ASSISTANT_INBOX",
            dir.join("valenx_chat_inbox.jsonl"),
        );
        assert_eq!(global_cmd_path(), dir.join("valenx_chat_cmd.jsonl"));
        assert_eq!(unit_cmd_path(3), dir.join("valenx_chat_cmd_u3.jsonl"));
        std::env::remove_var("VALENX_ASSISTANT_INBOX");
    }

    #[test]
    fn feed_paths_follow_feed_env_and_insert_unit_suffix() {
        let _g = env_lock();
        let base = std::env::temp_dir().join("valenx-mcp-bridge-feedtest");
        let feed = base.join("valenx_chat_feed.jsonl");
        std::env::set_var("VALENX_ASSISTANT_FEED", &feed);
        assert_eq!(base_feed_path(), feed);
        assert_eq!(
            unit_feed_path(7),
            base.join("valenx_chat_feed_u7.jsonl"),
            "unit feed must insert _u{{n}} before .jsonl"
        );
        assert_eq!(feed_path_for(None), feed);
        assert_eq!(
            feed_path_for(Some(7)),
            base.join("valenx_chat_feed_u7.jsonl")
        );
        std::env::remove_var("VALENX_ASSISTANT_FEED");
    }

    #[test]
    fn write_command_appends_one_jsonl_line() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-mcp-bridge-write-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cmd.jsonl");
        let a = write_command(&path, &serde_json::json!({"cmd":"new_unit"})).unwrap();
        let b = write_command(&path, &serde_json::json!({"cmd":"note","text":"hi"})).unwrap();
        assert_eq!(a, r#"{"cmd":"new_unit"}"#);
        assert_eq!(b, r#"{"cmd":"note","text":"hi"}"#);
        // Two distinct newline-delimited lines on disk.
        let body = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], r#"{"cmd":"new_unit"}"#);
        assert_eq!(lines[1], r#"{"cmd":"note","text":"hi"}"#);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_feed_since_returns_only_new_lines() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-mcp-bridge-feed-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let feed = dir.join("feed.jsonl");
        // Pre-existing line the tool should NOT see.
        std::fs::write(
            &feed,
            "{\"title\":\"You\",\"detail\":\"old\",\"kind\":\"user\"}\n",
        )
        .unwrap();
        let mark = feed_len(&feed);
        // valenx appends two ack lines after the mark, plus a half-written line.
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(&feed)
            .unwrap();
        writeln!(
            f,
            "{{\"title\":\"Claude\",\"detail\":\"ran view.front\",\"kind\":\"result\"}}"
        )
        .unwrap();
        writeln!(
            f,
            "{{\"title\":\"Claude\",\"detail\":\"Unit 1 ready\",\"kind\":\"ship\"}}"
        )
        .unwrap();
        write!(f, "{{\"title\":\"Claude\",\"detail\":\"half").unwrap(); // malformed tail
        let got = read_feed_since(&feed, mark);
        assert_eq!(got.len(), 2, "only the two complete new lines, got {got:?}");
        assert_eq!(got[0].detail, "ran view.front");
        assert_eq!(got[0].kind, "result");
        assert_eq!(got[1].detail, "Unit 1 ready");
        // The rendered block shows both acks, not the old line.
        let rendered = render_feed(&got);
        assert!(rendered.contains("ran view.front"));
        assert!(rendered.contains("Unit 1 ready"));
        assert!(!rendered.contains("old"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_feed_since_missing_file_is_empty() {
        let p = std::env::temp_dir().join("valenx-mcp-bridge-does-not-exist.jsonl");
        let _ = std::fs::remove_file(&p);
        assert!(read_feed_since(&p, 0).is_empty());
        assert_eq!(feed_len(&p), 0);
    }

    #[test]
    fn render_feed_empty_is_sentinel() {
        let s = render_feed(&[]);
        assert!(s.contains("no response from valenx"));
    }
}
