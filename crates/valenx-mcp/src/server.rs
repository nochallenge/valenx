//! MCP stdio server loop.
//!
//! We implement the minimal subset of MCP needed for tool calls:
//! `initialize`, `tools/list`, `tools/call`. The transport is line-
//! oriented JSON-RPC 2.0 over stdin/stdout.

use anyhow::Result;
use serde_json::{json, Value};
use tokio::io::{AsyncBufReadExt, AsyncReadExt, AsyncWrite, AsyncWriteExt, BufReader};

use crate::tools;

/// Hard cap on the byte length of a single JSON-RPC line we'll accept
/// from stdin.
///
/// Round-16 H2: pre-fix the server used `BufReader::lines()` which has
/// no per-line cap — a misbehaving (or malicious) MCP client could
/// stream gigabytes onto a single line and OOM the server before the
/// JSON parser saw anything. 10 MiB is generous for JSON-RPC (the
/// largest message MCP itself emits is `tools/list` at ~tens of KiB;
/// even a `tools/call` with embedded binary data tops out in the
/// low-MB range) while still being small enough to absorb without
/// host pressure.
pub const MAX_MCP_REQUEST_BYTES: usize = 10 * 1024 * 1024;

/// Run the MCP server on stdin / stdout. Blocks until stdin closes.
pub async fn serve_stdio() -> Result<()> {
    let stdin = tokio::io::stdin();
    let stdout = tokio::io::stdout();
    serve_stdio_with(BufReader::new(stdin), stdout).await
}

/// Inner serve loop, generic over the input + output transports. The
/// public [`serve_stdio`] wraps real stdin/stdout; tests inject
/// in-memory pipes so they can feed oversize lines without spawning a
/// subprocess.
pub async fn serve_stdio_with<R, W>(mut reader: BufReader<R>, mut writer: W) -> Result<()>
where
    R: tokio::io::AsyncRead + Unpin,
    W: AsyncWrite + Unpin,
{
    let mut buf: Vec<u8> = Vec::with_capacity(8 * 1024);
    loop {
        buf.clear();
        // Bounded read: a malicious client streaming gigabytes onto
        // one line can no longer OOM the server. We `take()` the
        // reader to one byte over the cap so we can DETECT overflow
        // (the reader yields exactly cap+1 bytes when the line is at
        // least that long) and reply with a JSON-RPC parse error.
        // Using `read_until` directly would still allocate up to the
        // cap; that's intended — the cap is the budget we accept.
        let limit_with_marker = (MAX_MCP_REQUEST_BYTES as u64) + 1;
        let mut limited = (&mut reader).take(limit_with_marker);
        let n = limited.read_until(b'\n', &mut buf).await?;
        if n == 0 {
            break; // EOF
        }
        // Detect cap overflow: we read more bytes than the cap permits,
        // OR we hit the marker byte without finding a newline (the
        // line is longer than cap + the marker). Either way we drain
        // the rest of the line so the next iteration starts on a fresh
        // JSON-RPC frame and emit a parse-error response.
        let overflowed_cap = buf.len() > MAX_MCP_REQUEST_BYTES;
        let no_newline_at_cap = buf.len() == (MAX_MCP_REQUEST_BYTES + 1)
            && !buf.ends_with(b"\n");
        if overflowed_cap || no_newline_at_cap {
            // Drain to end-of-line so the client's next frame lands
            // cleanly. We don't accumulate the drained bytes — just
            // skip past the newline (if any) so we resync.
            drain_to_newline(&mut reader).await?;
            let resp = json!({
                "jsonrpc": "2.0",
                "id": null,
                "error": {
                    "code": -32600,
                    "message": format!(
                        "request line exceeds {MAX_MCP_REQUEST_BYTES}-byte cap"
                    )
                }
            });
            let out = serde_json::to_string(&resp)?;
            writer.write_all(out.as_bytes()).await?;
            writer.write_all(b"\n").await?;
            writer.flush().await?;
            continue;
        }
        // Strip the trailing newline if present (read_until keeps it).
        if buf.ends_with(b"\n") {
            buf.pop();
            if buf.ends_with(b"\r") {
                buf.pop();
            }
        }
        // Skip blank lines (matches the pre-fix behaviour of
        // `next_line()` + `line.trim().is_empty()`).
        let line = match std::str::from_utf8(&buf) {
            Ok(s) => s,
            Err(_) => {
                let resp = json!({
                    "jsonrpc": "2.0",
                    "id": null,
                    "error": { "code": -32700, "message": "request line is not valid UTF-8" }
                });
                let out = serde_json::to_string(&resp)?;
                writer.write_all(out.as_bytes()).await?;
                writer.write_all(b"\n").await?;
                writer.flush().await?;
                continue;
            }
        };
        if line.trim().is_empty() {
            continue;
        }
        let response = match serde_json::from_str::<Value>(line) {
            Ok(req) => handle_request(req).await,
            Err(e) => json!({
                "jsonrpc": "2.0",
                "id": null,
                "error": { "code": -32700, "message": format!("parse: {e}") }
            }),
        };
        let out = serde_json::to_string(&response)?;
        writer.write_all(out.as_bytes()).await?;
        writer.write_all(b"\n").await?;
        writer.flush().await?;
    }
    Ok(())
}

/// After we reject an oversize line, drain the rest of the client's
/// line bytes (up to and including the next newline) so the next
/// iteration starts on a fresh JSON-RPC frame. Reads at most
/// `MAX_MCP_REQUEST_BYTES` extra bytes of garbage — anything more
/// and we treat the client as dead and stop draining.
///
/// Uses the BufReader's `fill_buf` + `consume` pair so we can scan
/// for the newline chunk-by-chunk without ever pulling bytes past
/// the newline out of the buffer. When we spot '\n' in the current
/// chunk, we consume exactly up to (and including) that byte —
/// every byte AFTER the newline stays in the BufReader's internal
/// buffer, where the next iteration of the main `serve_stdio_with`
/// loop picks them up as the start of the next JSON-RPC frame.
async fn drain_to_newline<R>(reader: &mut BufReader<R>) -> Result<()>
where
    R: tokio::io::AsyncRead + Unpin,
{
    let max_drain = MAX_MCP_REQUEST_BYTES;
    let mut drained = 0usize;
    loop {
        let chunk = reader.fill_buf().await?;
        if chunk.is_empty() {
            // EOF mid-line. Next iteration of the main loop sees
            // EOF and shuts down cleanly.
            return Ok(());
        }
        if let Some(pos) = chunk.iter().position(|&b| b == b'\n') {
            // Consume up to AND INCLUDING the newline. Bytes after
            // `pos` stay in the buffer for the main loop.
            let consume_n = pos + 1;
            reader.consume(consume_n);
            return Ok(());
        }
        // No newline in this chunk — consume it all and keep
        // looking, but only up to the drain budget.
        let chunk_len = chunk.len();
        reader.consume(chunk_len);
        drained = drained.saturating_add(chunk_len);
        if drained >= max_drain {
            return Ok(());
        }
    }
}

async fn handle_request(req: Value) -> Value {
    let id = req.get("id").cloned().unwrap_or(Value::Null);
    let method = req.get("method").and_then(|v| v.as_str()).unwrap_or("");
    let params = req.get("params").cloned().unwrap_or(json!({}));
    let result = match method {
        "initialize" => json!({
            "protocolVersion": "2024-11-05",
            "capabilities": { "tools": {} },
            "serverInfo": { "name": "valenx-mcp", "version": env!("CARGO_PKG_VERSION") }
        }),
        "tools/list" => json!({ "tools": tools::list() }),
        "tools/call" => match tools::call(&params).await {
            Ok(v) => v,
            Err(e) => {
                return json!({
                    "jsonrpc": "2.0",
                    "id": id,
                    "error": { "code": -32000, "message": format!("{e}") }
                });
            }
        },
        other => {
            return json!({
                "jsonrpc": "2.0",
                "id": id,
                "error": { "code": -32601, "message": format!("method not found: {other}") }
            });
        }
    };
    json!({ "jsonrpc": "2.0", "id": id, "result": result })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-16 H2 RED→GREEN: a single JSON-RPC line longer than the
    /// cap is rejected with a parse-error response, the server does
    /// NOT OOM, and a subsequent well-formed request on a fresh line
    /// is still processed correctly.
    ///
    /// Pre-fix: the server used `BufReader::lines()` (unbounded) so a
    /// hostile client emitting 20 MiB on one line would pre-allocate
    /// the whole line into memory before any parser saw it.
    #[tokio::test]
    async fn oversize_line_returns_parse_error_and_continues() {
        // 20 MiB of garbage followed by a newline, then a well-formed
        // initialize request on its own line. Post-fix the server
        // must reject frame #1 with a parse-error response and still
        // produce the normal initialize response for frame #2.
        let mut input: Vec<u8> = Vec::with_capacity(20 * 1024 * 1024 + 256);
        input.extend(std::iter::repeat_n(b'X', 20 * 1024 * 1024));
        input.push(b'\n');
        input.extend_from_slice(
            br#"{"jsonrpc":"2.0","id":1,"method":"initialize"}"#,
        );
        input.push(b'\n');

        let reader = BufReader::new(std::io::Cursor::new(input));
        let mut output: Vec<u8> = Vec::new();
        serve_stdio_with(reader, &mut output).await.unwrap();

        let text = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines.len(),
            2,
            "expected one error response + one initialize response, got: {text:?}"
        );
        // First response: parse error from the oversize frame.
        let first: Value = serde_json::from_str(lines[0]).unwrap();
        assert_eq!(first["error"]["code"], -32600);
        let msg = first["error"]["message"].as_str().unwrap();
        assert!(
            msg.contains("exceeds") && msg.contains("cap"),
            "got error message: {msg}"
        );
        // Second response: the normal initialize result. Confirms the
        // server resynced past the rejected frame and is still alive.
        let second: Value = serde_json::from_str(lines[1]).unwrap();
        assert_eq!(second["id"], 1);
        assert_eq!(second["result"]["protocolVersion"], "2024-11-05");
    }

    /// Sanity check: under-cap lines still flow through normally.
    #[tokio::test]
    async fn under_cap_line_processed_normally() {
        let input: &[u8] = br#"{"jsonrpc":"2.0","id":7,"method":"initialize"}
"#;
        let reader = BufReader::new(std::io::Cursor::new(input));
        let mut output: Vec<u8> = Vec::new();
        serve_stdio_with(reader, &mut output).await.unwrap();
        let text = String::from_utf8(output).unwrap();
        let v: Value = serde_json::from_str(text.trim()).unwrap();
        assert_eq!(v["id"], 7);
        assert!(v["result"].is_object());
    }
}
