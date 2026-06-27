//! valenx-remote — phone-as-2nd-screen companion for Valenx.
//!
//! Serves the live Valenx desktop window to a phone browser over the LAN and
//! relays the phone's touches back as operating-system mouse input, turning a
//! handset into a small touch-driven second screen for the running app.
//!
//! # How it works
//!
//! A tiny hand-rolled HTTP server (`std::net::TcpListener`, no web framework)
//! exposes three routes, every one of which is gated on a shared `--pin`:
//!
//! | Route             | Method | Purpose                                            |
//! |-------------------|--------|----------------------------------------------------|
//! | `/`               | GET    | The single-page client ([`PAGE`]).                 |
//! | `/frame.jpg`      | GET    | The current screen as a JPEG (needs `live-capture`).|
//! | `/input`          | POST   | A normalized `{nx, ny, kind}` pointer event.        |
//!
//! The phone page polls `/frame.jpg` continuously (each loaded frame schedules
//! the next request) and `POST`s pointer events as the user touches the image.
//!
//! # Default build is std-only (CI-green on headless Linux)
//!
//! The platform-specific parts — capturing the screen and injecting OS mouse
//! events — live behind the **off-by-default `live-capture`** Cargo feature so
//! the crate compiles on a headless Linux CI box with no GUI libraries. The
//! default build still contains the entire HTTP server, the served page, the
//! `--pin` authentication, and the full server-side input *validation*; only
//! the capture/injection back-ends are absent (and `/frame.jpg` then answers
//! `501 Not Implemented`). The user turns `live-capture` on when building on a
//! desktop that actually has a display. This mirrors the `binary-fmu` pattern
//! in `valenx-adapter-fmi`.
//!
//! # Security
//!
//! LAN-only and **no TLS**. Anyone on the same network who has the PIN can
//! control this PC, so do not expose the port to the internet and stop the
//! server (Ctrl-C) when finished. HTTPS (needed for some phone browser APIs)
//! is a documented follow-up — it requires pulling in a TLS dependency and is
//! intentionally out of scope here.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};

#[cfg(feature = "live-capture")]
mod capture;
#[cfg(feature = "live-capture")]
mod inject;

/// Default TCP port the server listens on.
pub const DEFAULT_PORT: u16 = 7333;
/// Default bind address (all interfaces, so the phone on the LAN can reach it).
pub const DEFAULT_BIND: &str = "0.0.0.0";
/// Default case-insensitive window-title substring used to pick the window to
/// capture (only relevant under `live-capture`).
pub const DEFAULT_TITLE: &str = "valenx";
/// Default frame rate cap, in frames per second.
pub const DEFAULT_FPS: u32 = 8;
/// Default JPEG quality (1..=100).
pub const DEFAULT_QUALITY: u8 = 60;

/// Runtime configuration, assembled from CLI flags with environment fallbacks.
#[derive(Debug, Clone)]
pub struct Config {
    /// Shared secret required on every request (the only authentication).
    pub pin: String,
    /// TCP port to listen on.
    pub port: u16,
    /// Bind address.
    pub bind: String,
    /// Case-insensitive window-title substring to capture.
    pub title: String,
    /// Frame rate cap (1..=60).
    pub fps: u32,
    /// JPEG quality (1..=100).
    pub quality: u8,
}

/// The single-page phone client, served verbatim at `GET /`.
///
/// It carries the `?pin=` from its own URL onto every sub-request, polls
/// `/frame.jpg` continuously (each `load`/`error` schedules the next fetch,
/// with a `?t=` cache-bust and a 300 ms backoff on error), and relays touch /
/// mouse events to `POST /input` as a normalized `{nx, ny, kind}` JSON body.
///
/// Accessibility / web-audit (F1) fixes baked in: `<html lang="en">`,
/// `<meta charset="utf-8">` first, a `role="status"` + `aria-live="polite"`
/// status hint, and an `alt` on the live-view image. No secrets are embedded
/// in the markup — the PIN arrives only via the URL the user opens.
pub const PAGE: &str = r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no, viewport-fit=cover">
<title>Valenx Remote</title>
<style>
  html, body { margin: 0; padding: 0; height: 100%; background: #16181f; overflow: hidden; }
  body { display: flex; align-items: center; justify-content: center; }
  #screen {
    max-width: 100vw; max-height: 100vh;
    width: auto; height: auto;
    touch-action: none; -webkit-user-select: none; user-select: none;
    -webkit-touch-callout: none;
  }
  #hint {
    position: fixed; left: 0; right: 0; bottom: 0;
    font: 12px -apple-system, system-ui, sans-serif;
    color: #8a90a0; text-align: center; padding: 4px;
    background: rgba(0,0,0,0.35); pointer-events: none;
  }
</style>
</head>
<body>
<img id="screen" alt="Valenx live view" draggable="false">
<div id="hint" role="status" aria-live="polite">Valenx Remote — tap &amp; drag to control the PC</div>
<script>
(function () {
  // Carry the PIN from this page's URL to every sub-request.
  var params = new URLSearchParams(window.location.search);
  var pin = params.get("pin") || "";
  var img = document.getElementById("screen");

  // Continuous polling: each loaded (or errored) frame schedules the next fetch.
  function nextFrame() {
    img.src = "/frame.jpg?pin=" + encodeURIComponent(pin) + "&t=" + Date.now() + "_" + Math.random();
  }
  img.addEventListener("load", function () { requestAnimationFrame(nextFrame); });
  img.addEventListener("error", function () { setTimeout(nextFrame, 300); });
  nextFrame();

  // Map a touch point to normalized [0,1] coords within the displayed image.
  function norm(touch) {
    var r = img.getBoundingClientRect();
    var nx = (touch.clientX - r.left) / Math.max(r.width, 1);
    var ny = (touch.clientY - r.top) / Math.max(r.height, 1);
    return { nx: Math.min(1, Math.max(0, nx)), ny: Math.min(1, Math.max(0, ny)) };
  }

  function send(kind, touch) {
    var p = norm(touch);
    fetch("/input?pin=" + encodeURIComponent(pin), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ nx: p.nx, ny: p.ny, kind: kind }),
      keepalive: true
    }).catch(function () {});
  }

  img.addEventListener("touchstart", function (e) {
    e.preventDefault();
    if (e.touches.length > 0) send("down", e.touches[0]);
  }, { passive: false });

  img.addEventListener("touchmove", function (e) {
    e.preventDefault();
    if (e.touches.length > 0) send("move", e.touches[0]);
  }, { passive: false });

  img.addEventListener("touchend", function (e) {
    e.preventDefault();
    // touchend has no active touches; use the last known changedTouches point.
    if (e.changedTouches.length > 0) send("up", e.changedTouches[0]);
  }, { passive: false });

  // Mouse fallback so the page is also usable from a desktop browser for testing.
  var mouseDown = false;
  img.addEventListener("mousedown", function (e) { mouseDown = true; send("down", e); });
  img.addEventListener("mousemove", function (e) { if (mouseDown) send("move", e); });
  window.addEventListener("mouseup", function (e) { if (mouseDown) { mouseDown = false; send("up", e); } });
})();
</script>
</body>
</html>
"#;

const HELP: &str = "valenx-remote — serve the live Valenx window to a phone browser; relay touches as mouse input.

USAGE:
    valenx-remote --pin <STR> [--port <U16>] [--bind <ADDR>] [--title <STR>] [--fps <N>] [--quality <1-100>]

FLAGS (env fallback in parentheses):
    --pin <STR>        REQUIRED. Shared secret required on every request.   (VALENX_REMOTE_PIN)
    --port <U16>       Port to listen on.            [default: 7333]        (VALENX_REMOTE_PORT)
    --bind <ADDR>      Bind address.                 [default: 0.0.0.0]     (VALENX_REMOTE_BIND)
    --title <STR>      Case-insensitive window-title substring to capture.
                                                     [default: valenx]     (VALENX_REMOTE_TITLE)
    --fps <N>          Frame rate (1-60).            [default: 8]           (VALENX_REMOTE_FPS)
    --quality <1-100>  JPEG quality.                 [default: 60]          (VALENX_REMOTE_QUALITY)
    -h, --help         Print this help.

EXAMPLE:
    valenx-remote --pin 1234

SECURITY: LAN-only, no TLS. Anyone on this network with the PIN can control this
PC. Do not expose the port to the internet. Stop with Ctrl-C when done.
";

/// An outcome of CLI parsing: either a ready [`Config`], or a request to print
/// help and exit cleanly.
#[derive(Debug)]
enum Cli {
    Run(Config),
    Help,
}

/// Parse argv (with environment fallbacks) into a [`Cli`].
///
/// `--pin` is mandatory: if it is absent both as a flag and as
/// `VALENX_REMOTE_PIN`, this returns an error (the server refuses to start
/// without authentication). Numeric flags are range-checked and clamped to
/// valid bounds so a typo can never produce a degenerate server.
fn parse_args(
    args: impl IntoIterator<Item = String>,
    env: impl Fn(&str) -> Option<String>,
) -> Result<Cli, String> {
    let mut pin = env("VALENX_REMOTE_PIN");
    let mut port = env("VALENX_REMOTE_PORT");
    let mut bind = env("VALENX_REMOTE_BIND");
    let mut title = env("VALENX_REMOTE_TITLE");
    let mut fps = env("VALENX_REMOTE_FPS");
    let mut quality = env("VALENX_REMOTE_QUALITY");

    let mut it = args.into_iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => return Ok(Cli::Help),
            "--pin" => pin = Some(next_value(&mut it, "--pin")?),
            "--port" => port = Some(next_value(&mut it, "--port")?),
            "--bind" => bind = Some(next_value(&mut it, "--bind")?),
            "--title" => title = Some(next_value(&mut it, "--title")?),
            "--fps" => fps = Some(next_value(&mut it, "--fps")?),
            "--quality" => quality = Some(next_value(&mut it, "--quality")?),
            other => return Err(format!("unknown argument: {other}")),
        }
    }

    let pin = pin
        .filter(|p| !p.is_empty())
        .ok_or_else(|| "--pin <STR> is required (the only auth for this server)".to_string())?;

    let port = match port {
        Some(s) => s
            .parse::<u16>()
            .map_err(|_| format!("invalid --port value: {s}"))?,
        None => DEFAULT_PORT,
    };
    let bind = bind.unwrap_or_else(|| DEFAULT_BIND.to_string());
    let title = title.unwrap_or_else(|| DEFAULT_TITLE.to_string());
    let fps = match fps {
        Some(s) => s
            .parse::<u32>()
            .map_err(|_| format!("invalid --fps value: {s}"))?
            .clamp(1, 60),
        None => DEFAULT_FPS,
    };
    let quality = match quality {
        Some(s) => s
            .parse::<u8>()
            .map_err(|_| format!("invalid --quality value: {s}"))?
            .clamp(1, 100),
        None => DEFAULT_QUALITY,
    };

    Ok(Cli::Run(Config {
        pin,
        port,
        bind,
        title,
        fps,
        quality,
    }))
}

/// Pull the value that must follow a value-taking flag, or error if missing.
fn next_value(it: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    it.next().ok_or_else(|| format!("missing value for {flag}"))
}

// ---------------------------------------------------------------------------
// Input model + validation (ALWAYS ON — not behind any feature)
// ---------------------------------------------------------------------------

/// The kind of pointer event the phone reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputKind {
    /// Pointer pressed down.
    Down,
    /// Pointer moved while pressed.
    Move,
    /// Pointer released.
    Up,
    /// A discrete tap (press + release at one point).
    Tap,
}

impl InputKind {
    /// Parse the `kind` field. Returns `None` for any unrecognized value.
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "down" => Some(Self::Down),
            "move" => Some(Self::Move),
            "up" => Some(Self::Up),
            "tap" => Some(Self::Tap),
            _ => None,
        }
    }
}

/// A validated pointer event: normalized coordinates in `[0, 1]` plus a kind.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InputEvent {
    /// Normalized horizontal position, guaranteed finite and in `[0, 1]`.
    pub nx: f64,
    /// Normalized vertical position, guaranteed finite and in `[0, 1]`.
    pub ny: f64,
    /// The pointer event kind.
    pub kind: InputKind,
}

/// Parse and **validate** an `/input` JSON body into an [`InputEvent`].
///
/// This is the security-critical, always-on guard required by the F1 web
/// audit: it is fail-loud and never panics on hostile input. It rejects, with
/// a descriptive error:
///
/// * a body that is not a JSON object with the three expected fields,
/// * `nx`/`ny` that are not finite (NaN / ±∞) or fall outside `[0, 1]`,
/// * a `kind` that is not one of `down` / `move` / `up` / `tap`.
///
/// A deliberately small hand-rolled parser is used (no `serde` dependency) so
/// the always-on core stays dependency-free; it tolerates whitespace and key
/// ordering but is otherwise strict.
pub fn parse_input(body: &str) -> Result<InputEvent, String> {
    let nx = extract_number(body, "nx").ok_or("missing or invalid 'nx'")?;
    let ny = extract_number(body, "ny").ok_or("missing or invalid 'ny'")?;
    let kind_str = extract_string(body, "kind").ok_or("missing or invalid 'kind'")?;

    if !nx.is_finite() || !ny.is_finite() {
        return Err("'nx'/'ny' must be finite numbers".to_string());
    }
    if !(0.0..=1.0).contains(&nx) || !(0.0..=1.0).contains(&ny) {
        return Err("'nx'/'ny' must be within [0, 1]".to_string());
    }
    let kind =
        InputKind::from_str(&kind_str).ok_or_else(|| format!("unknown 'kind': {kind_str:?}"))?;

    Ok(InputEvent { nx, ny, kind })
}

/// Extract a JSON number value for `"key"` from a flat object body.
///
/// Returns `None` if the key is absent or the value does not parse as `f64`.
/// (`f64::from_str` accepts `"NaN"`/`"inf"`; those non-finite cases are caught
/// and rejected by [`parse_input`].)
fn extract_number(body: &str, key: &str) -> Option<f64> {
    let raw = extract_raw_value(body, key)?;
    // A number value ends at the next comma or closing brace.
    let end = raw.find([',', '}']).unwrap_or(raw.len());
    raw[..end].trim().parse::<f64>().ok()
}

/// Extract a JSON string value for `"key"` from a flat object body.
fn extract_string(body: &str, key: &str) -> Option<String> {
    let raw = extract_raw_value(body, key)?;
    let raw = raw.trim_start();
    let rest = raw.strip_prefix('"')?;
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

/// Return the slice of `body` immediately after `"key" :`, or `None`.
fn extract_raw_value<'a>(body: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\"");
    let pos = body.find(&needle)?;
    let after = &body[pos + needle.len()..];
    let after = after.trim_start();
    let after = after.strip_prefix(':')?;
    Some(after.trim_start())
}

// ---------------------------------------------------------------------------
// HTTP routing (pure, testable — no socket involved)
// ---------------------------------------------------------------------------

/// A minimal HTTP response produced by [`route`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Response {
    /// Numeric status code (e.g. 200, 400, 403, 404, 501).
    pub status: u16,
    /// Human-readable reason phrase.
    pub reason: &'static str,
    /// `Content-Type` header value.
    pub content_type: &'static str,
    /// Response body bytes.
    pub body: Vec<u8>,
}

impl Response {
    fn text(status: u16, reason: &'static str, body: &str) -> Self {
        Response {
            status,
            reason,
            content_type: "text/plain; charset=utf-8",
            body: body.as_bytes().to_vec(),
        }
    }

    /// `403 Forbidden` for a missing or wrong PIN — the single auth failure.
    fn forbidden() -> Self {
        Self::text(403, "Forbidden", "forbidden: bad or missing pin")
    }

    /// `404 Not Found` for an unknown route.
    fn not_found() -> Self {
        Self::text(404, "Not Found", "not found")
    }
}

/// Does the request carry the correct PIN in its `?pin=` query?
///
/// Constant-ish string comparison against the configured secret. The PIN is
/// required on *every* route — there is no unauthenticated surface.
fn pin_ok(query: &str, expected: &str) -> bool {
    query_param(query, "pin").as_deref() == Some(expected)
}

/// Pull a single query-string parameter value (URL-decoded) by name.
fn query_param(query: &str, key: &str) -> Option<String> {
    for pair in query.split('&') {
        let mut kv = pair.splitn(2, '=');
        let k = kv.next()?;
        if k == key {
            let v = kv.next().unwrap_or("");
            return Some(url_decode(v));
        }
    }
    None
}

/// Minimal `application/x-www-form-urlencoded` decode (`%XX` + `+` → space).
fn url_decode(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'+' => {
                out.push(b' ');
                i += 1;
            }
            b'%' if i + 2 < bytes.len() => {
                let hi = (bytes[i + 1] as char).to_digit(16);
                let lo = (bytes[i + 2] as char).to_digit(16);
                if let (Some(hi), Some(lo)) = (hi, lo) {
                    out.push((hi * 16 + lo) as u8);
                    i += 3;
                } else {
                    out.push(bytes[i]);
                    i += 1;
                }
            }
            b => {
                out.push(b);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// Route a parsed request to a [`Response`].
///
/// This is the whole HTTP surface, expressed as a pure function so it can be
/// unit-tested without binding a socket. Every route validates the PIN first;
/// an absent/wrong PIN yields `403` regardless of path. `cfg` is the live
/// configuration; `capture` produces the current screen JPEG when the
/// `live-capture` feature is enabled.
pub fn route(
    method: &str,
    path: &str,
    query: &str,
    body: &str,
    cfg: &Config,
    capture: &mut dyn FnMut(&Config) -> Result<Vec<u8>, String>,
) -> Response {
    // Authentication gate: applies to every route.
    if !pin_ok(query, &cfg.pin) {
        return Response::forbidden();
    }

    match (method, path) {
        ("GET", "/") => Response {
            status: 200,
            reason: "OK",
            content_type: "text/html; charset=utf-8",
            body: PAGE.as_bytes().to_vec(),
        },
        ("GET", "/favicon.ico") => Response {
            status: 200,
            reason: "OK",
            content_type: "image/x-icon",
            body: Vec::new(),
        },
        ("GET", "/frame.jpg") => match capture(cfg) {
            Ok(jpeg) => Response {
                status: 200,
                reason: "OK",
                content_type: "image/jpeg",
                body: jpeg,
            },
            // Without `live-capture` the closure reports it is unavailable; we
            // surface a clear 501 rather than pretending to serve a frame.
            Err(msg) => Response::text(501, "Not Implemented", &msg),
        },
        ("POST", "/input") => match parse_input(body) {
            Ok(ev) => {
                #[cfg(feature = "live-capture")]
                {
                    // Best-effort: an injection failure must not crash the
                    // server or leak details to the client.
                    let _ = inject::dispatch(ev, cfg);
                }
                #[cfg(not(feature = "live-capture"))]
                {
                    let _ = ev; // validated, but nothing to inject without the feature
                }
                Response::text(200, "OK", "ok")
            }
            // F1 fix: malformed / out-of-range input is rejected fail-loud
            // with a 400 and a reason, never a panic.
            Err(msg) => Response::text(400, "Bad Request", &msg),
        },
        _ => Response::not_found(),
    }
}

// ---------------------------------------------------------------------------
// Socket plumbing (thin shell around `route`)
// ---------------------------------------------------------------------------

/// A parsed request line + body, extracted from a raw HTTP/1.x stream.
struct RawRequest {
    method: String,
    path: String,
    query: String,
    body: String,
}

/// Read and minimally parse one HTTP/1.x request from `stream`.
///
/// Reads headers up to the blank line, then reads exactly `Content-Length`
/// further bytes as the body (capped to a sane maximum so a hostile
/// `Content-Length` cannot exhaust memory). Returns `None` on a malformed or
/// empty request rather than panicking.
fn read_request(stream: &mut TcpStream) -> Option<RawRequest> {
    const MAX_HEADER: usize = 16 * 1024;
    const MAX_BODY: usize = 64 * 1024;

    let mut buf = Vec::new();
    let mut tmp = [0u8; 4096];
    // Read until we have the end-of-headers marker (or hit the header cap).
    let header_end = loop {
        if let Some(pos) = find_subslice(&buf, b"\r\n\r\n") {
            break pos + 4;
        }
        if buf.len() > MAX_HEADER {
            return None;
        }
        let n = stream.read(&mut tmp).ok()?;
        if n == 0 {
            // Connection closed before headers completed.
            if let Some(pos) = find_subslice(&buf, b"\r\n\r\n") {
                break pos + 4;
            }
            return None;
        }
        buf.extend_from_slice(&tmp[..n]);
    };

    let head = String::from_utf8_lossy(&buf[..header_end]);
    let mut lines = head.split("\r\n");
    let request_line = lines.next()?;
    let mut parts = request_line.split_whitespace();
    let method = parts.next()?.to_string();
    let target = parts.next()?.to_string();

    let (path, query) = match target.split_once('?') {
        Some((p, q)) => (p.to_string(), q.to_string()),
        None => (target, String::new()),
    };

    // Find Content-Length (case-insensitive).
    let mut content_length = 0usize;
    for line in lines {
        if let Some((k, v)) = line.split_once(':') {
            if k.trim().eq_ignore_ascii_case("content-length") {
                content_length = v.trim().parse::<usize>().unwrap_or(0).min(MAX_BODY);
            }
        }
    }

    // Body bytes already buffered, plus whatever else we still need to read.
    let mut body = buf[header_end..].to_vec();
    while body.len() < content_length {
        let n = stream.read(&mut tmp).ok()?;
        if n == 0 {
            break;
        }
        body.extend_from_slice(&tmp[..n]);
        if body.len() > MAX_BODY {
            break;
        }
    }
    body.truncate(content_length.min(MAX_BODY));

    Some(RawRequest {
        method,
        path,
        query,
        body: String::from_utf8_lossy(&body).into_owned(),
    })
}

/// Find the first index of `needle` within `hay`.
fn find_subslice(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// Serialize and write a [`Response`] back to the client.
///
/// `Cache-Control: no-store` is sent on every response so the phone never
/// caches a stale frame or the page.
fn write_response(stream: &mut TcpStream, resp: &Response) {
    let header = format!(
        "HTTP/1.1 {} {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nCache-Control: no-store, no-cache, must-revalidate\r\nConnection: close\r\n\r\n",
        resp.status,
        resp.reason,
        resp.content_type,
        resp.body.len(),
    );
    // Best-effort write; a dropped phone connection is not fatal.
    let _ = stream.write_all(header.as_bytes());
    let _ = stream.write_all(&resp.body);
    let _ = stream.flush();
}

/// Produce the current screen as a JPEG, or an error string.
///
/// With `live-capture` this delegates to the `capture` module; without it,
/// it always reports that capture is unavailable (driving a `501`).
fn capture_frame(cfg: &Config) -> Result<Vec<u8>, String> {
    #[cfg(feature = "live-capture")]
    {
        capture::capture_jpeg(cfg)
    }
    #[cfg(not(feature = "live-capture"))]
    {
        let _ = cfg;
        Err("frame capture unavailable: build with --features live-capture".to_string())
    }
}

/// Handle a single accepted connection end-to-end.
fn handle_connection(mut stream: TcpStream, cfg: &Config) {
    if let Some(req) = read_request(&mut stream) {
        let resp = route(
            &req.method,
            &req.path,
            &req.query,
            &req.body,
            cfg,
            &mut capture_frame,
        );
        write_response(&mut stream, &resp);
    }
}

/// Bind the listener and serve forever (one thread per connection).
fn serve(cfg: Config) -> std::io::Result<()> {
    let addr = format!("{}:{}", cfg.bind, cfg.port);
    let listener = TcpListener::bind(&addr)?;

    eprintln!("valenx-remote starting");
    eprintln!(
        "Open on your phone (same Wi-Fi):  http://<this-pc-ip>:{}/?pin=<your-pin>",
        cfg.port
    );
    eprintln!(
        "SECURITY: Anyone on this network with the PIN can control this PC; stop the server (Ctrl-C) when done. No TLS."
    );
    #[cfg(not(feature = "live-capture"))]
    eprintln!(
        "NOTE: built without `live-capture` — /frame.jpg returns 501 and touches are validated but not injected."
    );

    let cfg = std::sync::Arc::new(cfg);
    for stream in listener.incoming() {
        match stream {
            Ok(stream) => {
                let cfg = std::sync::Arc::clone(&cfg);
                std::thread::spawn(move || handle_connection(stream, &cfg));
            }
            Err(e) => eprintln!("accept error: {e}"),
        }
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(args, |k| std::env::var(k).ok()) {
        Ok(Cli::Help) => {
            print!("{HELP}");
        }
        Ok(Cli::Run(cfg)) => {
            if let Err(e) = serve(cfg) {
                eprintln!("error: {e}");
                std::process::exit(1);
            }
        }
        Err(msg) => {
            eprintln!("error: {msg}\n");
            eprint!("{HELP}");
            std::process::exit(2);
        }
    }
}

// ---------------------------------------------------------------------------
// Tests (always-on core — no feature required)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn test_cfg() -> Config {
        Config {
            pin: "secret".to_string(),
            port: DEFAULT_PORT,
            bind: DEFAULT_BIND.to_string(),
            title: DEFAULT_TITLE.to_string(),
            fps: DEFAULT_FPS,
            quality: DEFAULT_QUALITY,
        }
    }

    /// A capture stub used by routing tests so they never touch a display.
    fn ok_capture(_: &Config) -> Result<Vec<u8>, String> {
        Ok(vec![0xFF, 0xD8, 0xFF, 0xD9]) // tiny fake JPEG SOI..EOI
    }

    // --- PIN enforcement -----------------------------------------------------

    #[test]
    fn missing_pin_is_rejected_on_every_route() {
        let cfg = test_cfg();
        for (m, p) in [("GET", "/"), ("GET", "/frame.jpg"), ("POST", "/input")] {
            let r = route(m, p, "", "", &cfg, &mut ok_capture);
            assert_eq!(r.status, 403, "{m} {p} with no pin must be 403");
            assert_eq!(r.body, b"forbidden: bad or missing pin");
        }
    }

    #[test]
    fn wrong_pin_is_rejected() {
        let cfg = test_cfg();
        let r = route("GET", "/", "pin=nope", "", &cfg, &mut ok_capture);
        assert_eq!(r.status, 403);
    }

    #[test]
    fn correct_pin_serves_page() {
        let cfg = test_cfg();
        let r = route("GET", "/", "pin=secret", "", &cfg, &mut ok_capture);
        assert_eq!(r.status, 200);
        assert_eq!(r.content_type, "text/html; charset=utf-8");
        assert!(r.body.starts_with(b"<!doctype html>"));
    }

    #[test]
    fn pin_is_url_decoded() {
        let mut cfg = test_cfg();
        cfg.pin = "a b".to_string();
        // "a b" percent-encoded as "a%20b" must authenticate.
        let r = route("GET", "/", "pin=a%20b", "", &cfg, &mut ok_capture);
        assert_eq!(r.status, 200);
    }

    // --- /input validation (F1: fail-loud, no panic) ------------------------

    #[test]
    fn input_accepts_well_formed_event() {
        let ev = parse_input(r#"{"nx":0.25,"ny":0.75,"kind":"down"}"#).unwrap();
        assert_eq!(ev.kind, InputKind::Down);
        assert!((ev.nx - 0.25).abs() < 1e-12);
        assert!((ev.ny - 0.75).abs() < 1e-12);
    }

    #[test]
    fn input_accepts_all_kinds() {
        for (s, k) in [
            ("down", InputKind::Down),
            ("move", InputKind::Move),
            ("up", InputKind::Up),
            ("tap", InputKind::Tap),
        ] {
            let body = format!(r#"{{"nx":0.5,"ny":0.5,"kind":"{s}"}}"#);
            assert_eq!(parse_input(&body).unwrap().kind, k);
        }
    }

    #[test]
    fn input_rejects_nan_without_panicking() {
        assert!(parse_input(r#"{"nx":NaN,"ny":0.5,"kind":"down"}"#).is_err());
    }

    #[test]
    fn input_rejects_infinity() {
        assert!(parse_input(r#"{"nx":inf,"ny":0.5,"kind":"down"}"#).is_err());
        assert!(parse_input(r#"{"nx":0.5,"ny":-inf,"kind":"down"}"#).is_err());
    }

    #[test]
    fn input_rejects_out_of_range_coords() {
        assert!(parse_input(r#"{"nx":1.5,"ny":0.5,"kind":"move"}"#).is_err());
        assert!(parse_input(r#"{"nx":0.5,"ny":-0.01,"kind":"move"}"#).is_err());
        assert!(parse_input(r#"{"nx":-0.0001,"ny":0.5,"kind":"move"}"#).is_err());
    }

    #[test]
    fn input_rejects_bad_kind() {
        assert!(parse_input(r#"{"nx":0.5,"ny":0.5,"kind":"explode"}"#).is_err());
        assert!(parse_input(r#"{"nx":0.5,"ny":0.5,"kind":""}"#).is_err());
    }

    #[test]
    fn input_rejects_missing_fields() {
        assert!(parse_input(r#"{"nx":0.5,"ny":0.5}"#).is_err());
        assert!(parse_input(r#"{"kind":"down"}"#).is_err());
        assert!(parse_input("").is_err());
        assert!(parse_input("not json at all").is_err());
        assert!(parse_input("{}").is_err());
    }

    #[test]
    fn input_route_rejects_malformed_with_400_not_panic() {
        let cfg = test_cfg();
        let r = route(
            "POST",
            "/input",
            "pin=secret",
            r#"{"nx":2,"ny":0.5,"kind":"down"}"#,
            &cfg,
            &mut ok_capture,
        );
        assert_eq!(r.status, 400);
    }

    #[test]
    fn input_route_accepts_valid_with_200() {
        let cfg = test_cfg();
        let r = route(
            "POST",
            "/input",
            "pin=secret",
            r#"{"nx":0.5,"ny":0.5,"kind":"tap"}"#,
            &cfg,
            &mut ok_capture,
        );
        assert_eq!(r.status, 200);
        assert_eq!(r.body, b"ok");
    }

    // --- Routing edges -------------------------------------------------------

    #[test]
    fn unknown_route_is_404() {
        let cfg = test_cfg();
        let r = route("GET", "/nope", "pin=secret", "", &cfg, &mut ok_capture);
        assert_eq!(r.status, 404);
    }

    #[test]
    fn frame_route_surfaces_capture_error_as_501() {
        let cfg = test_cfg();
        let mut failing = |_: &Config| Err("capture unavailable".to_string());
        let r = route("GET", "/frame.jpg", "pin=secret", "", &cfg, &mut failing);
        assert_eq!(r.status, 501);
    }

    #[test]
    fn frame_route_serves_jpeg_on_success() {
        let cfg = test_cfg();
        let r = route("GET", "/frame.jpg", "pin=secret", "", &cfg, &mut ok_capture);
        assert_eq!(r.status, 200);
        assert_eq!(r.content_type, "image/jpeg");
    }

    // --- Served HTML contains the required (F1) tags ------------------------

    #[test]
    fn page_has_required_meta_and_accessibility_tags() {
        // charset first, viewport, lang, title, alt text.
        assert!(PAGE.contains(r#"<meta charset="utf-8">"#));
        assert!(PAGE.contains(r#"<meta name="viewport""#));
        assert!(PAGE.contains(r#"<html lang="en">"#));
        assert!(PAGE.contains("<title>Valenx Remote</title>"));
        assert!(PAGE.contains(r#"alt="Valenx live view""#));
        // F1: status hint is announced to assistive tech.
        assert!(PAGE.contains(r#"role="status""#));
        assert!(PAGE.contains(r#"aria-live="polite""#));
        // No secret baked into the page.
        assert!(!PAGE.contains("secret"));
    }

    #[test]
    fn charset_meta_appears_before_title() {
        let charset = PAGE.find("charset").expect("charset present");
        let title = PAGE.find("<title>").expect("title present");
        assert!(charset < title, "charset must come before <title>");
    }

    // --- CLI parsing ---------------------------------------------------------

    fn no_env(_: &str) -> Option<String> {
        None
    }

    #[test]
    fn cli_requires_pin() {
        let r = parse_args(Vec::<String>::new(), no_env);
        assert!(r.is_err(), "no --pin and no env must error");
    }

    #[test]
    fn cli_empty_pin_is_rejected() {
        let args = vec!["--pin".to_string(), String::new()];
        assert!(parse_args(args, no_env).is_err());
    }

    #[test]
    fn cli_parses_pin_and_defaults() {
        let args = vec!["--pin".to_string(), "1234".to_string()];
        match parse_args(args, no_env).unwrap() {
            Cli::Run(cfg) => {
                assert_eq!(cfg.pin, "1234");
                assert_eq!(cfg.port, DEFAULT_PORT);
                assert_eq!(cfg.bind, DEFAULT_BIND);
                assert_eq!(cfg.fps, DEFAULT_FPS);
                assert_eq!(cfg.quality, DEFAULT_QUALITY);
            }
            Cli::Help => panic!("expected Run"),
        }
    }

    #[test]
    fn cli_env_fallback_supplies_pin() {
        let env = |k: &str| (k == "VALENX_REMOTE_PIN").then(|| "envpin".to_string());
        match parse_args(Vec::<String>::new(), env).unwrap() {
            Cli::Run(cfg) => assert_eq!(cfg.pin, "envpin"),
            Cli::Help => panic!("expected Run"),
        }
    }

    #[test]
    fn cli_flag_overrides_and_clamps() {
        let args = vec![
            "--pin".to_string(),
            "p".to_string(),
            "--port".to_string(),
            "9000".to_string(),
            "--fps".to_string(),
            "999".to_string(),
            "--quality".to_string(),
            "0".to_string(),
        ];
        match parse_args(args, no_env).unwrap() {
            Cli::Run(cfg) => {
                assert_eq!(cfg.port, 9000);
                assert_eq!(cfg.fps, 60, "fps clamped to 60");
                assert_eq!(cfg.quality, 1, "quality clamped to >=1");
            }
            Cli::Help => panic!("expected Run"),
        }
    }

    #[test]
    fn cli_bad_port_errors() {
        let args = vec![
            "--pin".to_string(),
            "p".to_string(),
            "--port".to_string(),
            "not-a-number".to_string(),
        ];
        assert!(parse_args(args, no_env).is_err());
    }

    #[test]
    fn cli_help_flag() {
        assert!(matches!(
            parse_args(vec!["-h".to_string()], no_env).unwrap(),
            Cli::Help
        ));
        assert!(matches!(
            parse_args(vec!["--help".to_string()], no_env).unwrap(),
            Cli::Help
        ));
    }

    #[test]
    fn cli_unknown_arg_errors() {
        let args = vec!["--pin".to_string(), "p".to_string(), "--bogus".to_string()];
        assert!(parse_args(args, no_env).is_err());
    }

    // --- url_decode ----------------------------------------------------------

    #[test]
    fn url_decode_handles_percent_and_plus() {
        assert_eq!(url_decode("a%20b"), "a b");
        assert_eq!(url_decode("a+b"), "a b");
        assert_eq!(url_decode("plain"), "plain");
        // Malformed trailing percent is passed through, not panicked on.
        assert_eq!(url_decode("bad%"), "bad%");
        assert_eq!(url_decode("bad%2"), "bad%2");
    }
}
