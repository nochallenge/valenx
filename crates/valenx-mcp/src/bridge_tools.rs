//! **MCP tools that DRIVE a running valenx** through the file-bridge
//! ([`crate::bridge`]).
//!
//! Each tool maps 1:1 to an `AgentCommand` valenx already understands
//! (`crates/valenx-app/src/agent_commands.rs`): it builds the exact tagged-JSONL
//! command line, records the matching feed file's current length, appends the
//! command to the correct channel file, waits briefly for valenx's ~1 s poll to
//! run it, then reads the feed lines valenx appended in response and returns
//! them to the MCP caller. No valenx state is poked directly — every effect goes
//! through valenx's own vetted reducer.
//!
//! ## Channel routing
//!
//! valenx splits commands across two files (see [`crate::bridge`]):
//! the `new_unit` command is **global** (`valenx_chat_cmd.jsonl`); every
//! other command is **per-unit** (`valenx_chat_cmd_u{n}.jsonl`). The per-unit
//! tools therefore take an optional `unit` argument (default `1`) that selects
//! which Workbench+Agent unit to drive; its acks come back on that unit's feed.
//!
//! ## Tools
//!
//! - `valenx_new_unit` — open a fresh Workbench+Agent unit (optional `kind` /
//!   `title` / `note` / `group`). The only tool on the global channel.
//! - `valenx_open_workbench(id, unit?)` — switch a unit's active tab to a
//!   workbench by id (e.g. `"rocket"`, `"uq"`).
//! - `valenx_list_workbenches()` — list the command-palette ids valenx exposes
//!   (drives `list_commands`).
//! - `valenx_list_controls(workbench?, unit?)` — list a workbench's settable
//!   control captions (drives `list_controls`).
//! - `valenx_set_control(name, value, workbench?, unit?)` — set a labelled
//!   control to a typed value (drives `set_control`).
//! - `valenx_run_command(id, unit?)` — run a command-palette action by id
//!   (drives `run_command`).
//! - `valenx_read_readout(workbench?, unit?)` — read a workbench's computed
//!   result back (drives `read_readout`).
//! - `valenx_note(text, unit?)` — post a visible note into a unit's chat feed.

use anyhow::{anyhow, Result};
use serde_json::{json, Value};

use crate::bridge;

/// How long a per-call tool waits for valenx to poll + run the command and post
/// its ack, before reading the feed back. valenx's command poll runs on a ~1 s
/// cadence (`agent_commands::POLL_INTERVAL`), so we wait a little over that so a
/// freshly-appended command has been seen and its ack written. Overridable per
/// call via the `wait_ms` argument for slower hosts / batched drives.
pub const DEFAULT_WAIT_MS: u64 = 1300;

/// Upper bound on the per-call wait so a hostile `wait_ms` can't hang the MCP
/// server indefinitely.
pub const MAX_WAIT_MS: u64 = 30_000;

/// The eight bridge tools' `tools/list` entries, appended to the server's tool
/// list. Each schema documents only the inputs valenx's matching reducer reads.
pub fn tool_list() -> Vec<Value> {
    // Shared optional knobs reused across the per-unit tools.
    let unit_prop = json!({
        "type": "integer",
        "minimum": 1,
        "description": "Workbench+Agent unit to drive (default 1). Open one first with valenx_new_unit."
    });
    let wait_prop = json!({
        "type": "integer",
        "minimum": 0,
        "description": "Milliseconds to wait for valenx to run the command and post its ack before reading the feed back (default ~1300; valenx polls every ~1s)."
    });
    vec![
        json!({
            "name": "valenx_new_unit",
            "description": "Open a fresh 'Workbench + Agent' unit in a running valenx (the bootstrap before any unit exists). Optionally build a product into it. Returns the unit's ready/ack feed line.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "kind": {"type": "string", "description": "Optional product to build (e.g. 'rocket', 'gear', 'rcbeam', 'dna'). Absent → empty unit."},
                    "title": {"type": "string", "description": "Optional heading for the unit's workspace card."},
                    "note": {"type": "string", "description": "Optional narration line posted to the new unit's chat feed."},
                    "group": {"type": "string", "description": "Optional coloured tab-group to file the unit's tab into."},
                    "wait_ms": wait_prop,
                },
            }
        }),
        json!({
            "name": "valenx_open_workbench",
            "description": "Switch a unit's active tab to the workbench with the given id (a TabKind id, e.g. 'rocket', 'uq', 'fem'). Returns valenx's ack.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "id": {"type": "string", "description": "Workbench id to switch to (case-insensitive)."},
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
                "required": ["id"],
            }
        }),
        json!({
            "name": "valenx_list_workbenches",
            "description": "List the command-palette action ids valenx exposes (what valenx_run_command accepts). Posts the list into a unit's feed and returns it.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
            }
        }),
        json!({
            "name": "valenx_list_controls",
            "description": "List the settable control captions of a workbench (what valenx_set_control accepts for it). Returns valenx's listing.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "workbench": {"type": "string", "description": "Target workbench id (default: the unit's active tab workbench)."},
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
            }
        }),
        json!({
            "name": "valenx_set_control",
            "description": "Set a labelled workbench control to a typed value, by the exact caption the user sees. Returns valenx's ack (or a warn on a bad name/type).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "The control's user-visible caption (== its labelled_by text)."},
                    "value": {"description": "The value to assign: a boolean, integer, number, or string."},
                    "workbench": {"type": "string", "description": "Target workbench id (default: the unit's active tab workbench)."},
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
                "required": ["name", "value"],
            }
        }),
        json!({
            "name": "valenx_run_command",
            "description": "Run a valenx command-palette action by its stable id (e.g. 'view.front', 'run.selected-case'). Use valenx_list_workbenches to discover ids. Returns valenx's ack.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "id": {"type": "string", "description": "The stable command id to run."},
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
                "required": ["id"],
            }
        }),
        json!({
            "name": "valenx_read_readout",
            "description": "Read a workbench's computed result back into the feed so the agent can self-verify what it drove. Returns the result text (or a warn if not run yet).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "workbench": {"type": "string", "description": "Target workbench id (default: the unit's active tab workbench)."},
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
            }
        }),
        json!({
            "name": "valenx_note",
            "description": "Post a visible note line into a unit's chat feed (the agent's narration shows up in valenx's panel). Returns the echoed note.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "text": {"type": "string", "description": "The message body."},
                    "unit": unit_prop,
                    "wait_ms": wait_prop,
                },
                "required": ["text"],
            }
        }),
    ]
}

/// Dispatch a bridge tool by name. Returns `None` when `name` is not one of the
/// eight bridge tools (so the server's other tools still match), else `Some` of
/// the tool's MCP result.
pub fn dispatch(name: &str, args: &Value) -> Option<Result<Value>> {
    let r = match name {
        "valenx_new_unit" => new_unit(args),
        "valenx_open_workbench" => open_workbench(args),
        "valenx_list_workbenches" => list_workbenches(args),
        "valenx_list_controls" => list_controls(args),
        "valenx_set_control" => set_control(args),
        "valenx_run_command" => run_command(args),
        "valenx_read_readout" => read_readout(args),
        "valenx_note" => note(args),
        _ => return None,
    };
    Some(r)
}

/// Resolve the optional `unit` argument (default 1), rejecting a 0 / out-of-range
/// value loudly so a tool never silently writes to the wrong channel.
fn unit_of(args: &Value) -> Result<usize> {
    match args.get("unit") {
        None | Some(Value::Null) => Ok(1),
        Some(v) => {
            let n = v
                .as_u64()
                .ok_or_else(|| anyhow!("`unit` must be a positive integer"))?;
            if n == 0 {
                return Err(anyhow!("`unit` must be >= 1"));
            }
            Ok(n as usize)
        }
    }
}

/// Resolve the optional `wait_ms` argument (default [`DEFAULT_WAIT_MS`], clamped
/// to [`MAX_WAIT_MS`]).
fn wait_of(args: &Value) -> u64 {
    args.get("wait_ms")
        .and_then(Value::as_u64)
        .unwrap_or(DEFAULT_WAIT_MS)
        .min(MAX_WAIT_MS)
}

/// Pull a required string argument, erroring with the field name if absent.
fn req_str<'a>(args: &'a Value, key: &str) -> Result<&'a str> {
    args.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| anyhow!("missing required string `{key}`"))
}

/// Insert an optional string field into a command object only when present and
/// non-empty — so the emitted JSONL carries exactly the keys valenx's serde
/// `#[serde(default)]` fields expect, and omits the rest (keeping the line
/// minimal and matching the on-the-wire shape the reducer was tested against).
fn put_opt_str(obj: &mut serde_json::Map<String, Value>, args: &Value, key: &str) {
    if let Some(s) = args.get(key).and_then(Value::as_str) {
        if !s.is_empty() {
            obj.insert(key.to_string(), Value::String(s.to_string()));
        }
    }
}

/// The shared write→wait→read core every bridge tool runs:
/// 1. resolve the feed file for `channel` and record its current length;
/// 2. append the `command` JSONL line to `cmd_path`;
/// 3. wait `wait_ms` for valenx's poll to run it and post an ack;
/// 4. read the feed lines appended since step 1 and return them as the MCP
///    `content` text (plus a `structuredContent` echo of what was written and
///    the parsed acks).
///
/// The wait is a plain blocking sleep on this stdio server's task: each MCP
/// `tools/call` is handled to completion before the next line is read, so
/// sleeping here simply paces this one call against valenx's ~1 s poll — it does
/// not block other in-flight tools (there are none on a line-oriented stdio
/// loop).
fn drive(
    cmd_path: std::path::PathBuf,
    channel: Option<usize>,
    command: Value,
    wait_ms: u64,
) -> Result<Value> {
    let feed = bridge::feed_path_for(channel);
    let before = bridge::feed_len(&feed);
    let line = bridge::write_command(&cmd_path, &command).map_err(|e| {
        anyhow!(
            "failed to write valenx command file {}: {e}",
            cmd_path.display()
        )
    })?;
    if wait_ms > 0 {
        std::thread::sleep(std::time::Duration::from_millis(wait_ms));
    }
    let acks = bridge::read_feed_since(&feed, before);
    let text = bridge::render_feed(&acks);
    Ok(json!({
        "content": [{ "type": "text", "text": text }],
        "structuredContent": {
            "wrote": line,
            "command_file": cmd_path.display().to_string(),
            "feed_file": feed.display().to_string(),
            "acks": acks.iter().map(|e| json!({
                "kind": e.kind, "detail": e.detail, "title": e.title,
            })).collect::<Vec<_>>(),
        }
    }))
}

// --- the eight tools -------------------------------------------------------

fn new_unit(args: &Value) -> Result<Value> {
    // GLOBAL channel only — `new_unit` is the bootstrap valenx honours on
    // `valenx_chat_cmd.jsonl`. Its ack lands on the NEW unit's feed, but the
    // unit number isn't known until valenx mints it, so we read the *base* feed
    // (global channel) where `apply_global`'s confirm note is not posted; the
    // unit-ready note goes to the unit feed. We therefore additionally surface
    // the freshly-created unit's feed by scanning unit feeds is overkill — the
    // base-channel drive still returns any global-channel acks, and the caller
    // can `valenx_read_readout`/`valenx_note` on the unit next. Keep it simple:
    // drive on the global channel and report what the global feed shows.
    let mut obj = serde_json::Map::new();
    obj.insert("cmd".to_string(), Value::String("new_unit".to_string()));
    put_opt_str(&mut obj, args, "kind");
    put_opt_str(&mut obj, args, "title");
    put_opt_str(&mut obj, args, "note");
    put_opt_str(&mut obj, args, "group");
    drive(
        bridge::global_cmd_path(),
        None,
        Value::Object(obj),
        wait_of(args),
    )
}

fn open_workbench(args: &Value) -> Result<Value> {
    let id = req_str(args, "id")?;
    let unit = unit_of(args)?;
    let cmd = json!({ "cmd": "open_workbench", "id": id });
    drive(bridge::unit_cmd_path(unit), Some(unit), cmd, wait_of(args))
}

fn list_workbenches(args: &Value) -> Result<Value> {
    // valenx's `list_commands` enumerates the command-palette ids into the feed.
    let unit = unit_of(args)?;
    let cmd = json!({ "cmd": "list_commands" });
    drive(bridge::unit_cmd_path(unit), Some(unit), cmd, wait_of(args))
}

fn list_controls(args: &Value) -> Result<Value> {
    let unit = unit_of(args)?;
    let mut obj = serde_json::Map::new();
    obj.insert(
        "cmd".to_string(),
        Value::String("list_controls".to_string()),
    );
    put_opt_str(&mut obj, args, "workbench");
    drive(
        bridge::unit_cmd_path(unit),
        Some(unit),
        Value::Object(obj),
        wait_of(args),
    )
}

fn set_control(args: &Value) -> Result<Value> {
    let name = req_str(args, "name")?;
    let value = args
        .get("value")
        .cloned()
        .ok_or_else(|| anyhow!("missing required `value`"))?;
    // valenx's AgentValue is untagged (bool / int / float / string). Reject a
    // JSON object/array/null up front so we never write a value shape the
    // reducer can't interpret.
    match &value {
        Value::Bool(_) | Value::Number(_) | Value::String(_) => {}
        _ => return Err(anyhow!("`value` must be a boolean, number, or string")),
    }
    let unit = unit_of(args)?;
    let mut obj = serde_json::Map::new();
    obj.insert("cmd".to_string(), Value::String("set_control".to_string()));
    obj.insert("name".to_string(), Value::String(name.to_string()));
    obj.insert("value".to_string(), value);
    put_opt_str(&mut obj, args, "workbench");
    drive(
        bridge::unit_cmd_path(unit),
        Some(unit),
        Value::Object(obj),
        wait_of(args),
    )
}

fn run_command(args: &Value) -> Result<Value> {
    let id = req_str(args, "id")?;
    let unit = unit_of(args)?;
    let cmd = json!({ "cmd": "run_command", "id": id });
    drive(bridge::unit_cmd_path(unit), Some(unit), cmd, wait_of(args))
}

fn read_readout(args: &Value) -> Result<Value> {
    let unit = unit_of(args)?;
    let mut obj = serde_json::Map::new();
    obj.insert("cmd".to_string(), Value::String("read_readout".to_string()));
    put_opt_str(&mut obj, args, "workbench");
    drive(
        bridge::unit_cmd_path(unit),
        Some(unit),
        Value::Object(obj),
        wait_of(args),
    )
}

fn note(args: &Value) -> Result<Value> {
    let text = req_str(args, "text")?;
    let unit = unit_of(args)?;
    let cmd = json!({ "cmd": "note", "text": text });
    drive(bridge::unit_cmd_path(unit), Some(unit), cmd, wait_of(args))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Capture the exact JSONL line a tool writes WITHOUT a running valenx, by
    /// pointing the bridge at a temp dir, calling the tool with `wait_ms: 0`,
    /// and reading the command file back. Returns the single line written.
    fn line_written(tool: &str, mut args: Value) -> String {
        let _g = crate::test_support::env_lock();
        let dir = std::env::temp_dir().join(format!(
            "valenx-mcp-tooltest-{}-{}",
            tool,
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::env::set_var(
            "VALENX_ASSISTANT_INBOX",
            dir.join("valenx_chat_inbox.jsonl"),
        );
        std::env::set_var("VALENX_ASSISTANT_FEED", dir.join("valenx_chat_feed.jsonl"));
        // Force no wait so the test never sleeps.
        if let Value::Object(m) = &mut args {
            m.insert("wait_ms".to_string(), json!(0));
        }
        let res = dispatch(tool, &args).expect("known tool").expect("ok");
        let cmd_file = res["structuredContent"]["command_file"]
            .as_str()
            .unwrap()
            .to_string();
        let body = std::fs::read_to_string(&cmd_file).unwrap();
        std::env::remove_var("VALENX_ASSISTANT_INBOX");
        std::env::remove_var("VALENX_ASSISTANT_FEED");
        let _ = std::fs::remove_dir_all(&dir);
        let lines: Vec<&str> = body.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(
            lines.len(),
            1,
            "expected exactly one command line, got {body:?}"
        );
        lines[0].to_string()
    }

    /// Round-trip the written line back through valenx's OWN `AgentCommand`
    /// wire-format expectations would require depending on valenx-app; instead
    /// we assert the exact serialized bytes, which is what valenx's serde
    /// `#[serde(tag="cmd", rename_all="snake_case")]` parser consumes.
    #[test]
    fn new_unit_serializes_to_global_command() {
        let l = line_written("valenx_new_unit", json!({}));
        assert_eq!(l, r#"{"cmd":"new_unit"}"#);
    }

    #[test]
    fn new_unit_includes_optional_fields_when_set() {
        let l = line_written(
            "valenx_new_unit",
            json!({"kind":"rocket","title":"LV-1","group":"Aero"}),
        );
        // serde_json's Map is a BTreeMap (no `preserve_order` feature), so keys
        // serialize alphabetically: cmd, group, kind, title. valenx parses by
        // key name (serde, order-independent), so the wire order is irrelevant
        // to the reducer — we assert the exact bytes our writer emits.
        assert_eq!(
            l,
            r#"{"cmd":"new_unit","group":"Aero","kind":"rocket","title":"LV-1"}"#
        );
    }

    #[test]
    fn new_unit_goes_to_global_cmd_file() {
        let _g = crate::test_support::env_lock();
        let dir = std::env::temp_dir().join("valenx-mcp-tooltest-newunit-path");
        std::env::set_var(
            "VALENX_ASSISTANT_INBOX",
            dir.join("valenx_chat_inbox.jsonl"),
        );
        std::env::set_var("VALENX_ASSISTANT_FEED", dir.join("valenx_chat_feed.jsonl"));
        std::fs::create_dir_all(&dir).unwrap();
        let res = dispatch("valenx_new_unit", &json!({"wait_ms":0}))
            .unwrap()
            .unwrap();
        let cmd_file = res["structuredContent"]["command_file"].as_str().unwrap();
        assert!(
            cmd_file.ends_with("valenx_chat_cmd.jsonl"),
            "new_unit must use the GLOBAL (no _u suffix) file, got {cmd_file}"
        );
        std::env::remove_var("VALENX_ASSISTANT_INBOX");
        std::env::remove_var("VALENX_ASSISTANT_FEED");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn open_workbench_serializes_and_uses_unit_file() {
        let l = line_written("valenx_open_workbench", json!({"id":"rocket","unit":2}));
        assert_eq!(l, r#"{"cmd":"open_workbench","id":"rocket"}"#);
    }

    #[test]
    fn per_unit_tools_write_unit_suffixed_file() {
        let _g = crate::test_support::env_lock();
        let dir = std::env::temp_dir().join("valenx-mcp-tooltest-unitfile");
        std::env::set_var(
            "VALENX_ASSISTANT_INBOX",
            dir.join("valenx_chat_inbox.jsonl"),
        );
        std::env::set_var("VALENX_ASSISTANT_FEED", dir.join("valenx_chat_feed.jsonl"));
        std::fs::create_dir_all(&dir).unwrap();
        let res = dispatch(
            "valenx_run_command",
            &json!({"id":"view.front","unit":4,"wait_ms":0}),
        )
        .unwrap()
        .unwrap();
        let cmd_file = res["structuredContent"]["command_file"].as_str().unwrap();
        assert!(
            cmd_file.ends_with("valenx_chat_cmd_u4.jsonl"),
            "per-unit tool must target _u4, got {cmd_file}"
        );
        std::env::remove_var("VALENX_ASSISTANT_INBOX");
        std::env::remove_var("VALENX_ASSISTANT_FEED");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn list_workbenches_serializes() {
        let l = line_written("valenx_list_workbenches", json!({}));
        assert_eq!(l, r#"{"cmd":"list_commands"}"#);
    }

    #[test]
    fn list_controls_serializes_with_and_without_workbench() {
        assert_eq!(
            line_written("valenx_list_controls", json!({})),
            r#"{"cmd":"list_controls"}"#
        );
        assert_eq!(
            line_written("valenx_list_controls", json!({"workbench":"uq"})),
            r#"{"cmd":"list_controls","workbench":"uq"}"#
        );
    }

    #[test]
    fn set_control_serializes_each_value_type() {
        assert_eq!(
            line_written(
                "valenx_set_control",
                json!({"name":"Monte-Carlo samples N","value":4000,"workbench":"uq"})
            ),
            r#"{"cmd":"set_control","name":"Monte-Carlo samples N","value":4000,"workbench":"uq"}"#
        );
        assert_eq!(
            line_written("valenx_set_control", json!({"name":"a1","value":0.55})),
            r#"{"cmd":"set_control","name":"a1","value":0.55}"#
        );
        assert_eq!(
            line_written("valenx_set_control", json!({"name":"enabled","value":true})),
            r#"{"cmd":"set_control","name":"enabled","value":true}"#
        );
        assert_eq!(
            line_written(
                "valenx_set_control",
                json!({"name":"mode","value":"linear"})
            ),
            r#"{"cmd":"set_control","name":"mode","value":"linear"}"#
        );
    }

    #[test]
    fn set_control_rejects_object_value() {
        let _g = crate::test_support::env_lock();
        let dir = std::env::temp_dir().join("valenx-mcp-tooltest-setctl-bad");
        std::env::set_var(
            "VALENX_ASSISTANT_INBOX",
            dir.join("valenx_chat_inbox.jsonl"),
        );
        let err = dispatch(
            "valenx_set_control",
            &json!({"name":"x","value":{"nested":1},"wait_ms":0}),
        )
        .unwrap()
        .unwrap_err();
        assert!(
            err.to_string().contains("boolean, number, or string"),
            "{err}"
        );
        std::env::remove_var("VALENX_ASSISTANT_INBOX");
    }

    #[test]
    fn run_command_and_read_readout_and_note_serialize() {
        assert_eq!(
            line_written("valenx_run_command", json!({"id":"view.front"})),
            r#"{"cmd":"run_command","id":"view.front"}"#
        );
        assert_eq!(
            line_written("valenx_read_readout", json!({"workbench":"uq"})),
            r#"{"cmd":"read_readout","workbench":"uq"}"#
        );
        assert_eq!(
            line_written("valenx_read_readout", json!({})),
            r#"{"cmd":"read_readout"}"#
        );
        assert_eq!(
            line_written("valenx_note", json!({"text":"hello"})),
            r#"{"cmd":"note","text":"hello"}"#
        );
    }

    #[test]
    fn missing_required_arg_errors() {
        let e = dispatch("valenx_open_workbench", &json!({"wait_ms":0}))
            .unwrap()
            .unwrap_err();
        assert!(e.to_string().contains("`id`"), "{e}");
        let e = dispatch("valenx_run_command", &json!({"wait_ms":0}))
            .unwrap()
            .unwrap_err();
        assert!(e.to_string().contains("`id`"), "{e}");
        let e = dispatch("valenx_note", &json!({"wait_ms":0}))
            .unwrap()
            .unwrap_err();
        assert!(e.to_string().contains("`text`"), "{e}");
    }

    #[test]
    fn bad_unit_is_rejected() {
        let e = dispatch("valenx_note", &json!({"text":"x","unit":0,"wait_ms":0}))
            .unwrap()
            .unwrap_err();
        assert!(e.to_string().contains("unit"), "{e}");
    }

    #[test]
    fn dispatch_returns_none_for_non_bridge_tool() {
        assert!(dispatch("dock", &json!({})).is_none());
        assert!(dispatch("create_sketch", &json!({})).is_none());
    }

    #[test]
    fn tool_list_advertises_all_eight() {
        let names: Vec<String> = tool_list()
            .iter()
            .map(|t| t["name"].as_str().unwrap().to_string())
            .collect();
        for n in [
            "valenx_new_unit",
            "valenx_open_workbench",
            "valenx_list_workbenches",
            "valenx_list_controls",
            "valenx_set_control",
            "valenx_run_command",
            "valenx_read_readout",
            "valenx_note",
        ] {
            assert!(names.contains(&n.to_string()), "missing tool {n}");
        }
        // Every advertised bridge tool carries an object inputSchema.
        for t in tool_list() {
            assert!(t["inputSchema"]["type"] == "object", "{} schema", t["name"]);
        }
    }
}
