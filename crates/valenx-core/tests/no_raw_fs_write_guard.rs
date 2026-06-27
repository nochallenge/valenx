//! R30 structural guard: no raw durable write in production source.
//!
//! Every round since R20 hand-grepped for un-migrated `std::fs::write`
//! calls and routed them through the crash-safe
//! `valenx_core::io_caps::atomic_write_{str,bytes}` helpers (unique
//! sidecar → fsync → rename). R29 finally automated the grep with a
//! source-walking test — but that test matched only the literal text
//! `fs::write(`. R30 review caught the predictable sister gap: two
//! adapters (ctffind, curves) wrote subprocess input decks via
//! `File::create(path)? + write_all(..)?`, a shape the textual scan
//! never saw. A torn/concurrent write fed a truncated deck to the child.
//!
//! This rewrite parses every `crates/**/src/**/*.rs` file with `syn`
//! (a real Rust AST, so comments / string literals / `cfg(test)` are
//! handled correctly — no naive text stripping) and flags the WHOLE
//! raw-durable-write surface:
//!
//!   - `fs::write(..)` / `std::fs::write(..)`
//!   - `File::create(..)` / `File::create_new(..)` (and the
//!     `std::fs::File::create*` spellings)
//!   - any `OpenOptions::new()....open(..)` builder chain
//!
//! Anything the visitor flags must be either MIGRATED to
//! `atomic_write_*` (durable small-state persistence) or registered in
//! [`ALLOWLIST`] with a one-line reason (legit non-atomic site:
//! streaming export, subprocess stdout/stderr redirect sink,
//! append-only log, or a canonical/local atomic-write implementation).
//!
//! Design notes:
//!
//! - `#[cfg(test)]` modules/functions and `#[test]` functions are
//!   skipped (the visitor does not recurse into them), so test fixtures
//!   that legitimately `File::create(..).unwrap()` never trip the guard.
//!   There is deliberately **no** `.unwrap()` / `.expect(..)` escape
//!   hatch: a *production* `fs::write(..).expect(..)` is exactly the
//!   kind of regression R30 review wanted closed, and the previous
//!   textual guard would have let it slip.
//! - A `syn::parse_file` failure PANICS with the path rather than being
//!   silently skipped — a file the guard can't parse is a file the guard
//!   is blind to, and that blindness is itself a defect.
//! - The allowlist matches by `(path_suffix, fn_name)` so it survives
//!   line drift; the enclosing function name comes from the AST visitor.
//! - `build.rs` scripts live at the crate root, not under `src/`, so the
//!   `src/`-scoped walk excludes them: a build script writing generated
//!   code into `OUT_DIR` is build-private, not a runtime durability
//!   concern.

use std::path::{Path, PathBuf};

use syn::visit::Visit;

/// Genuine, reviewed exceptions: `(path_suffix, fn_name, reason)`.
///
/// A flag is suppressed iff its file path ends with `path_suffix` AND
/// the enclosing function name equals `fn_name`. `fn_name` is `"<file>"`
/// for a flag at item scope (outside any fn) — none exist today.
///
/// Every entry was read at its site and confirmed to be a legitimate
/// non-atomic write. Categories:
///
///   - **Subprocess stdout/stderr redirect sinks** — the `File` is moved
///     into `Stdio::from` / `cmd.stdout(..)`; the OS writes the child's
///     output, not us. Atomic-rename is meaningless for a live append
///     target the kernel owns.
///   - **Append-only audit hash-chain** — atomic rename would clobber
///     the prior chain; append is the correct, intended semantics.
///   - **Streaming binary exports** — intentionally streamed, not
///     buffered into memory. A torn export is a re-runnable artifact, not
///     corrupt durable state; buffering a multi-GB mesh just to rename it
///     is strictly worse.
///   - **Canonical / local atomic-write implementations** — these ARE
///     the atomic write (or the documented local copy for crates that
///     can't depend on valenx-core because of the core→fields→mesh
///     cycle).
const ALLOWLIST: &[(&str, &str, &str)] = &[
    // ── Local atomic-write primitive + intentional append sinks ────────
    (
        "crates/valenx-recipes/src/lib.rs",
        "atomic_write",
        "local atomic-write primitive: sidecar temp + sync_all + rename, mirrors io_caps::atomic_write_str (keeps this crate pure-std/serde, no valenx-core dep)",
    ),
    (
        "crates/valenx-app/src/assistant_workbench.rs",
        "append_line",
        "append-only .jsonl agent-chat bridge: an atomic full-file replace would truncate the prior appended lines; incremental append is the intended semantics",
    ),
    (
        "crates/valenx-mcp/src/bridge.rs",
        "write_command",
        "append-only .jsonl command queue to a LIVE valenx GUI that polls the file by byte offset each frame: an atomic temp+rename would swap the inode out from under the running reader and truncate prior queued commands; incremental append is the intended (and required) wire semantics, mirroring assistant_workbench::append_line",
    ),
    // ── Subprocess stdout/stderr redirect sinks ────────────────────────
    // The created File handle is moved into `Stdio::from(..)` /
    // `cmd.stdout(..)`; the kernel writes the child's output to it. There
    // is no durable state WE author, and atomic-rename of a live append
    // target the OS owns is meaningless.
    (
        "crates/valenx-core/src/executor.rs",
        "submit",
        "subprocess redirect sinks: stdout.log/stderr.log FDs handed to cmd.stdout/stderr; OS writes, not us",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-badread/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-fasttree/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-xtb/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-mafft/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout (the alignment) streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-nwchem/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-viennarna/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-linearfold/src/lib.rs",
        "run",
        "subprocess redirect sink: child stdout streamed to a file via Stdio::from(out_file)",
    ),
    (
        "crates/valenx-adapters/bio/valenx-adapter-samtools/src/lib.rs",
        "run_capture_stdout",
        "subprocess redirect sink: child stdout (flagstat) streamed to a file via Stdio::from(out_file)",
    ),
    // ── Re-creatable marker file (not durable content) ──────────────────
    (
        "crates/valenx-adapters/bio/valenx-adapter-ctffind/src/lib.rs",
        "prepare",
        "marker-file touch: `let _ = File::create(&output_txt)` so collect() finds a path even on partial runs; CTFFIND overwrites it on success — no content authored here",
    ),
    // ── Audit log: append-only hash-chain + advisory lock file ──────────
    (
        "crates/valenx-audit/src/lib.rs",
        "lock_log",
        "advisory lock-file open: OpenOptions create+read+write+truncate(false) for fs2 lock_exclusive; not a content write",
    ),
    (
        "crates/valenx-audit/src/lib.rs",
        "append",
        "append-only audit hash-chain: OpenOptions create+append + single-call write; atomic rename would break append semantics and the prev_hash chain",
    ),
    // ── Hardened read-only open (flagged because it is an OpenOptions
    //    chain, but `.read(true)` only — no write/create/append) ─────────
    (
        "crates/valenx-mcp/src/sandbox.rs",
        "open_no_follow",
        "read-only sandbox open: OpenOptions.read(true) + O_NOFOLLOW / FILE_FLAG_OPEN_REPARSE_POINT; not a write at all",
    ),
    // ── Streaming geometry exports (re-runnable, not durable state) ─────
    // These stream element-by-element to a BufWriter; atomic_write_bytes
    // would buffer the whole (potentially multi-GB) mesh in memory. A
    // torn export is a re-runnable artifact, not corrupt durable state.
    (
        "crates/valenx-mesh/src/stl_write.rs",
        "write_stl_binary",
        "streaming export — atomic_write would buffer whole mesh in memory; torn file is re-runnable, not durable state",
    ),
    (
        "crates/valenx-mesh/src/format/obj.rs",
        "write_path",
        "streaming OBJ export — atomic_write would buffer whole mesh in memory; torn file is re-runnable, not durable state",
    ),
    (
        "crates/valenx-mesh/src/format/ply.rs",
        "write_path",
        "streaming PLY export — atomic_write would buffer whole mesh in memory; torn file is re-runnable, not durable state",
    ),
    (
        "crates/valenx-occt-exchange/src/ply_writer_extended.rs",
        "ply_writer_extended",
        "streaming PLY-with-annotations export via write_extended(&mut BufWriter); atomic_write would buffer whole mesh; torn file is re-runnable",
    ),
    (
        "crates/valenx-occt-exchange/src/stl_writer_extended.rs",
        "write_ascii",
        "streaming ASCII-STL export to a BufWriter; atomic_write would buffer whole mesh; torn file is re-runnable",
    ),
    (
        "crates/valenx-app/src/mesh_toolbox.rs",
        "write_triangle_mesh_stl",
        "streaming export — atomic_write would buffer whole mesh in memory; torn file is re-runnable, not durable state",
    ),
    (
        "crates/valenx-app/src/mesh_toolbox.rs",
        "write_binary_stl",
        "streaming export — atomic_write would buffer whole mesh in memory; torn file is re-runnable, not durable state",
    ),
    // ── Canonical / local atomic-write implementations ─────────────────
    (
        "crates/valenx-core/src/io_caps.rs",
        "atomic_write_bytes",
        "canonical atomic-write implementation: the sidecar File + fsync + rename IS the atomic write",
    ),
    (
        "crates/valenx-core/src/io_caps.rs",
        "atomic_write_streaming",
        "canonical atomic-write implementation: streaming sidecar File + fsync + rename",
    ),
    (
        "crates/valenx-fields/src/pvd.rs",
        "atomic_write_local",
        "local atomic-write implementation: valenx-fields can't dep valenx-core (core→fields→mesh cycle)",
    ),
];

/// Resolve the workspace root from this crate's manifest dir.
/// `CARGO_MANIFEST_DIR` is `<root>/crates/valenx-core`, so `../../`
/// climbs back to `<root>`.
fn workspace_root() -> PathBuf {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest
        .parent() // <root>/crates
        .and_then(Path::parent) // <root>
        .expect("CARGO_MANIFEST_DIR has two ancestors")
        .to_path_buf()
}

/// Collect every `*.rs` file under `crates/**/src/**`, skipping any
/// directory named `tests`, `target`, `vendor`, `examples`, or
/// `benches`, and any path that runs through a `tests/` segment.
///
/// `examples` and `benches` hold non-production code (demo binaries and
/// criterion harnesses) that legitimately `File::create(..)` scratch
/// output; they are excluded for the same reason `tests` is. No such
/// directory lives under any crate `src/` today, so this is a latent
/// safeguard against a future example/bench tripping the guard.
fn collect_src_rs(root: &Path) -> Vec<PathBuf> {
    let crates_dir = root.join("crates");
    let mut out = Vec::new();
    let mut stack = vec![crates_dir];
    while let Some(dir) = stack.pop() {
        let entries = match std::fs::read_dir(&dir) {
            Ok(e) => e,
            Err(_) => continue,
        };
        for entry in entries.flatten() {
            let path = entry.path();
            let name = entry.file_name();
            let name = name.to_string_lossy();
            if path.is_dir() {
                if matches!(
                    name.as_ref(),
                    "tests" | "target" | "vendor" | "examples" | "benches"
                ) {
                    continue;
                }
                stack.push(path);
            } else if name.ends_with(".rs") {
                let rel = path.strip_prefix(root).unwrap_or(&path);
                let rel_str = rel.to_string_lossy().replace('\\', "/");
                if rel_str.contains("/src/") && !rel_str.contains("/tests/") {
                    out.push(path);
                }
            }
        }
    }
    out
}

/// A flagged raw-write site.
#[derive(Debug)]
struct Flag {
    rel_path: String,
    line: usize,
    enclosing_fn: String,
    snippet: String,
}

/// True if any attribute is `#[test]`, `#[cfg(test)]`, or a cfg whose
/// predicate mentions `test` (e.g. `#[cfg(all(test, ...))]`).
fn has_test_attr(attrs: &[syn::Attribute]) -> bool {
    attrs.iter().any(|attr| {
        let path = attr.path();
        if path.is_ident("test") {
            return true;
        }
        if path.is_ident("cfg") {
            // Render the cfg(...) tokens and look for a `test` ident.
            // Robust to `cfg(test)`, `cfg(all(test, feature = ..))`, etc.
            if let syn::Meta::List(list) = &attr.meta {
                return list
                    .tokens
                    .clone()
                    .into_iter()
                    .any(|tt| matches!(tt, proc_macro2::TokenTree::Ident(ref id) if id == "test"))
                    || list
                        .tokens
                        .to_string()
                        .split(|c: char| !c.is_alphanumeric() && c != '_')
                        .any(|w| w == "test");
            }
        }
        false
    })
}

/// Last segment ident of a path, e.g. `std::fs::write` → `"write"`.
fn last_segment(path: &syn::Path) -> Option<String> {
    path.segments.last().map(|s| s.ident.to_string())
}

/// The segment ident immediately before the last, e.g.
/// `std::fs::File::create` → `"File"`, `std::fs::write` → `"fs"`.
fn penultimate_segment(path: &syn::Path) -> Option<String> {
    let n = path.segments.len();
    if n >= 2 {
        Some(path.segments[n - 2].ident.to_string())
    } else {
        None
    }
}

/// Does this call path name a raw durable-write free function/associated
/// fn we care about? Returns a short label for the snippet if so.
fn raw_write_call_label(path: &syn::Path) -> Option<&'static str> {
    let last = last_segment(path)?;
    let prev = penultimate_segment(path);
    match (prev.as_deref(), last.as_str()) {
        // fs::write / std::fs::write
        (Some("fs"), "write") => Some("fs::write"),
        // File::create / File::create_new / std::fs::File::create*
        (Some("File"), "create") => Some("File::create"),
        (Some("File"), "create_new") => Some("File::create_new"),
        _ => None,
    }
}

/// Walk a method-call receiver chain to its root expression, returning
/// true if the chain is rooted at an `OpenOptions` builder
/// (`OpenOptions::new()` or a bare `OpenOptions` path).
fn receiver_roots_at_open_options(mut expr: &syn::Expr) -> bool {
    loop {
        match expr {
            syn::Expr::MethodCall(mc) => {
                expr = &mc.receiver;
            }
            // `OpenOptions::new()` — an ExprCall whose func path ends in
            // `OpenOptions :: new` (or `new` with `OpenOptions` before).
            syn::Expr::Call(call) => {
                if let syn::Expr::Path(p) = &*call.func {
                    let last = last_segment(&p.path);
                    let prev = penultimate_segment(&p.path);
                    if prev.as_deref() == Some("OpenOptions") && last.as_deref() == Some("new") {
                        return true;
                    }
                    // Also handle the rarer `fs::OpenOptions::new()` where
                    // the builder ident is the penultimate segment.
                    if last.as_deref() == Some("new")
                        && p.path.segments.iter().any(|s| s.ident == "OpenOptions")
                    {
                        return true;
                    }
                }
                return false;
            }
            // A bare `OpenOptions` path (unusual, but be safe).
            syn::Expr::Path(p) => {
                return p.path.segments.iter().any(|s| s.ident == "OpenOptions");
            }
            // Parenthesised / referenced / grouped wrappers: unwrap.
            syn::Expr::Paren(p) => expr = &p.expr,
            syn::Expr::Group(g) => expr = &g.expr,
            syn::Expr::Reference(r) => expr = &r.expr,
            _ => return false,
        }
    }
}

/// AST visitor that records raw-write flags, skipping `#[cfg(test)]` /
/// `#[test]` items and tracking the enclosing function name.
struct WriteVisitor<'a> {
    rel_path: &'a str,
    fn_stack: Vec<String>,
    flags: Vec<Flag>,
}

impl<'a> WriteVisitor<'a> {
    fn enclosing_fn(&self) -> String {
        self.fn_stack
            .last()
            .cloned()
            .unwrap_or_else(|| "<file>".to_string())
    }

    fn push_flag(&mut self, line: usize, label: &str) {
        self.flags.push(Flag {
            rel_path: self.rel_path.to_string(),
            line,
            enclosing_fn: self.enclosing_fn(),
            snippet: label.to_string(),
        });
    }
}

impl<'ast, 'a> Visit<'ast> for WriteVisitor<'a> {
    fn visit_item_fn(&mut self, node: &'ast syn::ItemFn) {
        if has_test_attr(&node.attrs) {
            return; // don't recurse into test fns
        }
        self.fn_stack.push(node.sig.ident.to_string());
        syn::visit::visit_item_fn(self, node);
        self.fn_stack.pop();
    }

    fn visit_impl_item_fn(&mut self, node: &'ast syn::ImplItemFn) {
        if has_test_attr(&node.attrs) {
            return;
        }
        self.fn_stack.push(node.sig.ident.to_string());
        syn::visit::visit_impl_item_fn(self, node);
        self.fn_stack.pop();
    }

    fn visit_item_mod(&mut self, node: &'ast syn::ItemMod) {
        if has_test_attr(&node.attrs) {
            return; // don't recurse into #[cfg(test)] modules
        }
        syn::visit::visit_item_mod(self, node);
    }

    fn visit_expr_call(&mut self, node: &'ast syn::ExprCall) {
        if let syn::Expr::Path(p) = &*node.func {
            if let Some(label) = raw_write_call_label(&p.path) {
                let line = p
                    .path
                    .segments
                    .last()
                    .map_or(0, |s| s.ident.span().start().line);
                self.push_flag(line, label);
            }
        }
        syn::visit::visit_expr_call(self, node);
    }

    fn visit_expr_method_call(&mut self, node: &'ast syn::ExprMethodCall) {
        if node.method == "open" && receiver_roots_at_open_options(&node.receiver) {
            let line = node.method.span().start().line;
            self.push_flag(line, "OpenOptions::open");
        }
        syn::visit::visit_expr_method_call(self, node);
    }
}

/// True if `flag` matches an [`ALLOWLIST`] entry.
///
/// Setting `GUARD_NO_ALLOWLIST=1` disables the allowlist so the test
/// prints EVERY flagged raw-write site (migrate-vs-allowlist triage for
/// a future round). This is fail-safe: it can only make the guard fail
/// *harder* (every legit site is reported too) — it can never let a real
/// regression pass — so it is safe to leave wired in but must never be
/// set in CI.
fn is_allowlisted(flag: &Flag) -> bool {
    if std::env::var("GUARD_NO_ALLOWLIST").is_ok() {
        return false; // triage mode: surface EVERY flag
    }
    ALLOWLIST.iter().any(|(suffix, fn_name, _)| {
        flag.rel_path.ends_with(suffix) && flag.enclosing_fn == *fn_name
    })
}

#[test]
fn no_raw_durable_write_in_production_source() {
    let root = workspace_root();
    let files = collect_src_rs(&root);
    assert!(
        files.len() > 100,
        "guard walked only {} files — workspace root resolution is probably wrong (root = {})",
        files.len(),
        root.display()
    );

    let mut offenders: Vec<String> = Vec::new();

    for file in &files {
        let rel = file
            .strip_prefix(&root)
            .unwrap_or(file)
            .to_string_lossy()
            .replace('\\', "/");

        let src = match std::fs::read_to_string(file) {
            Ok(s) => s,
            // A non-UTF-8 *.rs file can't contain Rust we'd flag; skip.
            Err(_) => continue,
        };

        // A parse failure means the guard is BLIND to this file — that is
        // itself a defect, so panic with the path rather than skip.
        let ast = syn::parse_file(&src).unwrap_or_else(|e| {
            panic!(
                "no_raw_fs_write_guard: failed to parse {rel} with syn ({e}). \
                 The guard cannot see writes in a file it can't parse; fix the \
                 parse or adjust the guard."
            )
        });

        let mut visitor = WriteVisitor {
            rel_path: &rel,
            fn_stack: Vec::new(),
            flags: Vec::new(),
        };
        visitor.visit_file(&ast);

        for flag in &visitor.flags {
            if is_allowlisted(flag) {
                continue;
            }
            offenders.push(format!(
                "{}:{} — {} (in fn `{}`)",
                flag.rel_path, flag.line, flag.snippet, flag.enclosing_fn
            ));
        }
    }

    assert!(
        offenders.is_empty(),
        "Found {} raw production durable-write call(s) — route each through \
         `valenx_core::io_caps::atomic_write_str` / `atomic_write_bytes` (or a \
         documented local helper for valenx-fields/valenx-mesh), or register a \
         genuine exception in ALLOWLIST with a one-line reason:\n{}",
        offenders.len(),
        offenders.join("\n")
    );
}
