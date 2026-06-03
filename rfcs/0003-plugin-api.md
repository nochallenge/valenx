# RFC 0003: Plugin API (WIT-based WASM Sandbox)

- **Status:** Accepted (initial design; subject to refinement once prototyped)
- **Author(s):** BDFL
- **Created:** 2026-04-21
- **Discussion PR:** (this commit)
- **Tracking issue:** TBD

---

## Summary

Define how third-party code extends Valenx. Plugins compile to
**WebAssembly** against **WIT (WebAssembly Interface Types)**
interfaces published by Valenx, load into the app at runtime via
`wasmtime`, and run in a sandbox with only the capabilities the user
grants. Plugins are distinct from in-tree adapters (RFC 0002); they're
for untrusted / third-party code.

---

## Motivation

A serious simulation suite needs extensibility the maintainers don't
directly control:

- A university lab wants a custom boundary condition for their research
- An industrial user wants a custom post-processor pulling from their
  internal systems
- A community contributor wants to wrap an OSS tool we haven't
  integrated yet
- A company wants proprietary extensions on top of open-source Valenx

Requirements:

1. **Safe** — a malicious plugin must not be able to read arbitrary
   files, exfiltrate data, or crash the host
2. **Language-agnostic** — Rust, C, C++, AssemblyScript, Go, Zig, all
   should work
3. **Versioned** — plugins compiled against API v1 work on host v1.x
4. **Hot-loadable** — load a plugin without restarting the app
5. **Discoverable** — a plugin manager UI; no manual config-file hacks
6. **Doesn't mix license classes** — a GPL plugin can exist without
   making Valenx GPL

Anti-requirements:

- Not a full OS. Plugins don't get arbitrary filesystem, arbitrary
  network, arbitrary subprocess spawning
- Not a scripting language. We have Python/Lua for that (embedded
  interpreters, different trust model)
- Not native dynamic libraries. `.so`/`.dll`/`.dylib` plugins are a
  security nightmare and would re-introduce ABI problems

---

## Guide-level explanation

### What is a plugin

A plugin is:

- A single `.wasm` file (or `.wasm` plus metadata, as a `.valenx-plugin` bundle)
- Compiled against a published WIT world
- Loaded by Valenx at startup (or runtime, via the plugin manager)
- Runs inside a `wasmtime` instance, sandboxed from the host filesystem
  except for capabilities the user granted

### Plugin kinds (the "worlds")

We publish several WIT worlds — interfaces a plugin can implement:

| World | What it does | Examples |
|-------|--------------|----------|
| `valenx:adapter@1.0` | A user-space adapter for a tool | Wrap an OSS tool not yet in-tree |
| `valenx:postprocess@1.0` | Read results, produce derived data | Compute custom aerodynamic coefficients |
| `valenx:case-template@1.0` | Parametric case generator | "Generate N airfoil meshes at sweep of angles" |
| `valenx:exporter@1.0` | Write results to a custom format | Export to proprietary LMS format |
| `valenx:importer@1.0` | Read an external format into Valenx | Import CAD from an obscure vendor format |
| `valenx:ui-panel@1.0` | Add a custom panel to the UI | Custom parameter editor |

A single plugin file can implement multiple worlds.

### Plugin package format

```
my-plugin.valenx-plugin/
├── plugin.toml         # manifest — name, author, version, capabilities
├── plugin.wasm         # the compiled WASM
├── README.md           # user-facing description
├── icon.png            # 256x256 plugin icon
└── assets/             # static assets the plugin reads (read-only)
```

The `.valenx-plugin` extension is a directory (like `.valenx`
projects). Users install a plugin by dropping the directory into
`~/.valenx/plugins/` or using the plugin manager.

### Manifest

`plugin.toml`:

```toml
[plugin]
id = "org.example.my-cfd-tool"     # reverse-DNS
name = "My CFD Tool"
version = "0.3.2"
api_version = "1.0"                 # WIT world version compiled against
authors = ["Jane Doe <jane@example.com>"]
license = "Apache-2.0"
description = "Wraps the My-CFD-Tool solver for Valenx"
homepage = "https://example.com/my-cfd-tool"
repository = "https://github.com/example/my-cfd-plugin"

[plugin.worlds]
# Which WIT worlds this plugin implements
adapter = true
exporter = false

[plugin.capabilities]
# Requested at install time; user approves or denies
read_files = ["assets/"]            # relative to plugin dir
write_files = []
network = []                         # list of hostnames, [] means none
subprocess = []                      # list of binaries, [] means none
env_vars = []                        # list of env-var names, [] means none

[plugin.physics]
# What this plugin handles (for UI categorization)
kind = ["cfd"]
```

### Capability model

Plugins cannot do anything except:

- Call back into the host via **exported host functions** defined in
  the WIT world
- Read files it was granted access to at install time
- Use the CPU and whatever memory was allocated to its instance

Everything else is denied by default: no filesystem roaming, no raw
syscalls, no network except through host-mediated APIs, no spawning
processes.

The user sees the requested capabilities at install and can:

- Approve all
- Approve selectively (e.g., allow file reads but deny network)
- Decline — plugin doesn't install

### Example plugin (Rust)

```rust
// Cargo.toml: crate-type = ["cdylib"], deps: wit-bindgen

wit_bindgen::generate!({
    world: "adapter",
    path: "../wit/valenx-adapter@1.0",
});

struct MyCfdAdapter;

impl Guest for MyCfdAdapter {
    fn info() -> AdapterInfo {
        AdapterInfo {
            id: "my-cfd-tool",
            display_name: "My CFD Tool",
            physics: vec![Physics::Cfd],
            version_range: "1.0..2.0".to_string(),
            license_mode: LicenseMode::Subprocess,
        }
    }

    fn probe() -> Result<ProbeReport, AdapterError> {
        let binary = host::find_binary("my-cfd-tool")?;
        Ok(ProbeReport {
            ok: true,
            binary_path: Some(binary),
            ..Default::default()
        })
    }

    fn prepare(case: Case, workdir: String) -> Result<PreparedJob, AdapterError> {
        // Write input file
        let input = translate(&case)?;
        host::write_file(&format!("{workdir}/input.txt"), &input)?;
        Ok(PreparedJob {
            workdir,
            native_command: vec!["my-cfd-tool".into(), "input.txt".into()],
            ..Default::default()
        })
    }

    // ... run(), collect()
}

export!(MyCfdAdapter);
```

Compiled with `cargo build --target wasm32-wasip2`, producing a
`.wasm` file that is put into the plugin package.

---

## Reference-level explanation

### Choice of toolkit

- **Runtime:** `wasmtime` 20+
  - Mature, production-quality, Bytecode Alliance governed
  - Supports WASM Component Model + WASI Preview 2
  - Fast compile-cache, fast instantiation
- **Interface description:** WIT (WebAssembly Interface Types)
  - Language-agnostic interfaces
  - Tooling: `wit-bindgen` for Rust/C/C++/Go/Python guests
- **Guest languages:**
  - First-class: Rust (examples, tutorials, reference plugins)
  - Supported: C, C++, AssemblyScript, Go (TinyGo)
  - Best-effort: Python (componentize-py), JS

### WIT file organization

Published under `wit/` in the main repo and versioned:

```
wit/
├── valenx-core@1.0/
│   ├── types.wit        # shared types (Case, Results, etc.)
│   ├── host.wit         # what the host exposes to plugins
│   └── world.wit
├── valenx-adapter@1.0/
│   └── world.wit        # adapter interface
├── valenx-postprocess@1.0/
├── valenx-case-template@1.0/
├── valenx-exporter@1.0/
├── valenx-importer@1.0/
└── valenx-ui-panel@1.0/
```

Every world is tagged with a SemVer. Plugins declare which version
they were compiled against in their manifest; Valenx enforces
compatibility at load time.

### Host capabilities (what plugins can call)

```wit
// wit/valenx-core@1.0/host.wit
package valenx:core@1.0;

interface host {
    // Logging
    log: func(level: log-level, message: string);

    // Progress reporting
    progress: func(pct: f32, message: string);

    // Filesystem (subject to capability grants)
    read-file: func(path: string) -> result<list<u8>, fs-error>;
    write-file: func(path: string, data: list<u8>) -> result<_, fs-error>;
    exists: func(path: string) -> bool;

    // Subprocess (subject to capability grants + subprocess sandbox)
    spawn: func(command: list<string>, env: list<tuple<string, string>>) -> result<process-handle, subprocess-error>;
    wait: func(handle: process-handle) -> result<exit-code, subprocess-error>;
    cancel: func(handle: process-handle);

    // Network (subject to capability grants; HTTP only, not raw sockets)
    http-get: func(url: string) -> result<list<u8>, net-error>;
    http-post: func(url: string, body: list<u8>, headers: list<tuple<string, string>>) -> result<list<u8>, net-error>;

    // Find a known binary on the system
    find-binary: func(name: string) -> result<string, fs-error>;
}
```

Every host function checks capabilities before executing. A plugin
that didn't request `subprocess` will get a capability-denied error if
it tries to `spawn`.

### Sandbox model

Per-plugin instance:

- **Memory:** 256 MB default, configurable up to 4 GB per plugin
- **CPU:** cooperative; plugins may be killed if they don't yield or
  they exceed a time budget
- **Filesystem:** only paths the user granted, plus a `/tmp`-like
  writable area under `~/.valenx/plugins/<id>/cache/`
- **Network:** HTTP/HTTPS only, only to hostnames in the approved list
- **Subprocess:** only binaries in the approved list, no shell
  interpretation
- **Environment variables:** only those in the approved list (each
  read is explicit)

### Determinism and caching

Plugins compile once on first load (`wasmtime`'s compile cache).
Subsequent loads are fast — a typical plugin starts in ~5 ms.

If a plugin's `.wasm` changes, its cache entry is invalidated.

### Version compatibility

| Host version | Plugin API version | Plugin loads? |
|--------------|--------------------|---------------|
| 1.5 | 1.0 | Yes — SemVer minor compat |
| 1.5 | 1.5 | Yes |
| 1.5 | 1.7 | No — host too old |
| 1.5 | 0.x | No — pre-stable |
| 2.0 | 1.0 | Only via compat shim (TBD) |

We publish a "compat layer" when making a breaking change, so 1.x
plugins can load in 2.x for at least one major version.

### Plugin lifecycle

1. **Discover** — scan `~/.valenx/plugins/` at startup, check manifests
2. **Load** — instantiate `wasmtime` module, call `init()` export
3. **Register** — host queries `info()` exports, adds to registries
4. **Execute** — calls happen during user actions
5. **Unload** — on app close or user-triggered disable; plugin's
   `cleanup()` called

Hot-reload: the user can disable a plugin from the plugin manager,
edit its `.wasm`, and re-enable. Re-instantiation is fast enough
(~50 ms including re-register).

### Plugin signing (future)

For plugins distributed through an official registry, signing is
required. For sideloaded plugins, a big orange warning appears at
install:

> **Unsigned plugin.** This plugin has not been verified by the Valenx
> team. Only install if you trust the source. Permissions requested:
> [list].

Detailed signing scheme deferred to a future RFC once we stand up a
registry.

### CLI support

```bash
valenx plugin install ./my-plugin.valenx-plugin
valenx plugin list
valenx plugin enable  org.example.my-cfd-tool
valenx plugin disable org.example.my-cfd-tool
valenx plugin remove  org.example.my-cfd-tool
valenx plugin info    org.example.my-cfd-tool
```

### Testing plugins

A test harness crate `valenx-plugin-test` lets plugin authors run
their plugin against a fake host, exercising every host function and
verifying capability checks.

---

## Drawbacks

- **WIT + component model is still young.** The ecosystem changes
  faster than our RFC cycle. We pin a `wasmtime` and `wit-bindgen`
  version and upgrade deliberately.
- **WASM has real overhead** for math-heavy work — ~1.5-2x slower than
  native. For plugins doing heavy numerics, we recommend they
  subprocess to a native tool rather than compute in WASM.
- **Capability model has a learning curve.** Plugin authors have to
  think about what they ask for; users have to understand what they
  approve.
- **Language support is uneven.** Rust and C work great; Go is
  usable but large; Python plugins are 20+ MB even for simple code.
  We document best practice (Rust for serious plugins).
- **We now own a plugin API** — breaking changes hurt downstream. We
  commit to the deprecation cycle in POLICIES.md.

---

## Rationale and alternatives

**Native dynamic libraries (`.so` / `.dll` plugins):**
Rejected.
- Security: arbitrary code in the process, no sandbox
- ABI: Rust ABI not stable; each host version requires rebuild
- Portability: per-OS binaries per plugin
- License risk: GPL plugins would taint the host

**Lua / Python scripting only:**
Rejected as the only option.
- Too slow for anything beyond scripting
- We already have embedded scripting; that's separate (smaller scope:
  user-written scripts in their own project, not redistributable
  plugins)
- No language choice for the plugin author

**Node.js-style NPM packages (JS):**
Rejected. Node is the wrong runtime for a native desktop app; WASM
covers the same "run third-party code safely" use case with better
performance and language choice.

**gRPC-based plugins:**
Considered. Out-of-process plugin via gRPC gives you sandbox-for-free
(OS-level process isolation), but startup cost is higher and the IPC
protocol becomes another versioned surface. Reasonable fallback if
WASM turns out to be too constraining; we'd scope a separate RFC.

---

## Prior art

- **Blender Python addons** — huge ecosystem, proves plugin systems
  work in 3D apps. Pythonic, but unsandboxed; addons can do anything.
- **VS Code extensions** — WASM-not-required, JS in Electron; we
  learn from the extension marketplace UX but don't copy the runtime
- **Steam Workshop / plugins in games** — signed content, marketplace,
  user ratings. Good model for our eventual plugin registry.
- **Wasmtime's own use in Shopify Functions, Fastly Compute@Edge,
  Fermyon Spin** — production validation that the stack works
- **Rust's proc-macros** — WASM for proc-macros is being discussed
  upstream for determinism; similar motivation
- **FreeCAD's workbench system** — Python-based, great extensibility,
  but no sandbox; users run arbitrary code with full permissions

---

## Unresolved questions

- **Cross-plugin communication.** Plugin A wants to call Plugin B. Is
  this allowed? Safe? Current plan: no for v1 — plugins can only talk
  to the host. Revisit later.
- **GPU access from plugins.** Viz plugins might want shaders. WebGPU
  through WASI is emerging; defer.
- **Long-running plugin state.** Does a plugin keep state across
  invocations, or is every call fresh? Plan: plugin chooses (init
  allocates state; each invocation gets it), but state is in-memory
  only, not persisted.
- **Signing, distribution, registry.** Scoped separately. For v1, load
  from local disk only.
- **Python plugin size.** Even a hello-world Python-to-WASM component
  is 20+ MB. Do we offer a shared Python runtime across plugins?
  Possible but adds coupling risk.
- **UI plugins' rendering.** `valenx:ui-panel` is listed but the
  rendering protocol is nontrivial — egui immediate-mode doesn't
  translate cleanly across the WASM boundary. May end up restricting
  UI plugins to a declarative widget set (like a retained-mode
  description the host renders).

---

## Future possibilities

- **Plugin registry** at `plugins.valenx.org` — curated, signed,
  reviewed; one-click install from inside the app
- **Plugin reviews / ratings** in the registry
- **Plugin composition** — compose multiple plugins into a workflow
  template
- **In-app plugin creation wizard** — "I want to add support for X;
  generate a plugin skeleton for me"
- **Commercial plugin support** — third parties selling Valenx
  extensions; the WASM sandbox + capability model makes this
  practical
- **Enterprise policy** — IT admins restricting which plugins can be
  loaded
