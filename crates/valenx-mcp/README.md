# valenx-mcp

An [MCP](https://modelcontextprotocol.io) (Model Context Protocol) server that
exposes **valenx** to any MCP-capable AI agent over **stdio**. It advertises two
families of tools:

1. **Capability tools** (in-process, GUI-free): native docking (`dock`,
   `dry_run`, `describe`, `list_adapters`) and a generative-design / parametric
   CAD set (`create_sketch`, `pad`, `pocket`, `revolve`, `fillet`, …).
2. **Drive-a-running-valenx tools** (this document's focus): `valenx_*` tools
   that steer a **live valenx GUI** through its existing **file-bridge** — the
   same newline-delimited-JSON command channel valenx polls every frame. An
   agent can open workbench tabs, set tool inputs, run actions, and read results
   back, entirely by name.

The MCP layer is **hand-rolled** (no third-party MCP SDK): a small,
well-specified JSON-RPC 2.0 line loop over stdin/stdout implementing
`initialize`, `tools/list`, and `tools/call`. This keeps the dependency surface
to `serde`/`serde_json`/`anyhow`/`tokio` (all already in the workspace) and
avoids pulling a heavier SDK tree into a tool that mostly does small file I/O.

## The drive-a-running-valenx tools

Each tool maps 1:1 onto an `AgentCommand` valenx already understands
(`crates/valenx-app/src/agent_commands.rs`). A tool call:

1. records the current length of the relevant **feed** file,
2. appends the exact tagged-JSONL command line to the correct **command** file,
3. waits briefly (~1.3 s by default; valenx polls on a ~1 s cadence) for valenx
   to run it and post an ack,
4. reads the feed lines valenx appended in response and returns them.

| MCP tool | valenx command | Channel | Purpose |
|---|---|---|---|
| `valenx_new_unit` | `new_unit` | global | Open a fresh "Workbench + Agent" unit (optional `kind`/`title`/`note`/`group`). The bootstrap before any unit exists. |
| `valenx_open_workbench(id, unit?)` | `open_workbench` | per-unit | Switch a unit's active tab to a workbench by id (`rocket`, `uq`, `fem`, …). |
| `valenx_list_workbenches(unit?)` | `list_commands` | per-unit | List the command-palette action ids `valenx_run_command` accepts. |
| `valenx_list_controls(workbench?, unit?)` | `list_controls` | per-unit | List a workbench's settable control captions. |
| `valenx_set_control(name, value, workbench?, unit?)` | `set_control` | per-unit | Set a labelled control to a typed value (bool/int/float/string). |
| `valenx_run_command(id, unit?)` | `run_command` | per-unit | Run a command-palette action by id (e.g. `view.front`). |
| `valenx_read_readout(workbench?, unit?)` | `read_readout` | per-unit | Read a workbench's computed result back to self-verify. |
| `valenx_note(text, unit?)` | `note` | per-unit | Post a visible note into a unit's chat feed. |

`new_unit` is the **only** command valenx honours on the *global* channel
(`valenx_chat_cmd.jsonl`); every other tool writes a **per-unit** channel
(`valenx_chat_cmd_u{n}.jsonl`) selected by the optional `unit` argument
(default `1`). Open a unit with `valenx_new_unit` first, then drive `unit: 1`
(or whatever number valenx reports as ready).

Every tool also accepts an optional `wait_ms` to tune how long it waits for the
ack (clamped to 30 s).

### Where the files live

The command + feed files are resolved from the **same** environment variables
valenx uses, so the server writes the file valenx reads:

- command files: in the **directory of `$VALENX_ASSISTANT_INBOX`** (or the OS
  temp dir if unset);
- feed files: at **`$VALENX_ASSISTANT_FEED`** (or an OS-temp default), with
  `_u{n}` inserted before `.jsonl` for unit `n`.

Launch valenx and the MCP server with both variables pointing at the **same
directory**, e.g.:

```
VALENX_ASSISTANT_INBOX=<dir>/valenx_chat_inbox.jsonl
VALENX_ASSISTANT_FEED=<dir>/valenx_chat_feed.jsonl
```

## Building

```sh
cargo build -p valenx-mcp --release
# binary: target/release/valenx-mcp(.exe)
```

## Registering with an MCP client

The server speaks MCP over stdio, so any stdio MCP client can launch it. Example
**Claude Desktop** `claude_desktop_config.json` stanza (point `command` at the
built binary; set the env so it shares valenx's bridge directory):

```json
{
  "mcpServers": {
    "valenx": {
      "command": "C:/Users/you/valenx/target/release/valenx-mcp.exe",
      "env": {
        "VALENX_ASSISTANT_INBOX": "C:/Users/you/AppData/Local/Temp/valenx_chat_inbox.jsonl",
        "VALENX_ASSISTANT_FEED": "C:/Users/you/AppData/Local/Temp/valenx_chat_feed.jsonl",
        "VALENX_MCP_SANDBOX_DIR": "C:/Users/you/AppData/Local/Temp/valenx-mcp"
      }
    }
  }
}
```

On macOS/Linux use the matching paths (`target/release/valenx-mcp`,
`$TMPDIR`/`/tmp`). Start valenx with the **same** `VALENX_ASSISTANT_INBOX` /
`VALENX_ASSISTANT_FEED` so the `valenx_*` tools reach the running window; then
in the client call `valenx_new_unit`, `valenx_open_workbench`, etc.

`VALENX_MCP_SANDBOX_DIR` only scopes the *capability* tools' file paths (docking
receptors/ligands/outputs); the drive-a-running-valenx tools do not read
arbitrary files.

## Security posture — LOCAL ONLY

This server and its bridge are confined to the **local user session**:

- **No network.** The bridge opens no socket, binds no port, and makes no
  network call. It drives valenx purely by reading and writing **local files**
  in a user-writable directory. Nothing here is reachable from another host.
- **Stdio transport.** The MCP server is a child process of the MCP client
  (stdin/stdout), not a listening service.
- **Indirection through valenx's own reducers.** A `valenx_*` tool can only
  append a command line; valenx executes it through the **same vetted
  tab/dock/workbench methods a user click would**, never a raw field poke. An
  agent can only drive a valenx the local user has already launched and pointed
  at the same `$VALENX_ASSISTANT_*` directory.
- **Bounded reads.** Command-file and feed reads are size-capped so a corrupt or
  runaway file can't exhaust memory.
- **Capability tools are path-sandboxed.** Every `*_path` argument to the
  docking tools must resolve under `$VALENX_MCP_SANDBOX_DIR` (TOCTOU-resistant
  open); see the crate docs.

## Tests

The test suite runs **without a running valenx**: it points the bridge at a
temp directory and asserts (a) each tool serializes to the exact JSONL command
line valenx's parser consumes, (b) per-unit vs. global channel file selection,
(c) feed parsing returns only the lines appended after the recorded offset
(skipping malformed/half-written lines), and (d) the `tools/list` schema
advertises every tool with an object input schema.

```sh
cargo test -p valenx-mcp
```
