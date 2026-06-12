//! Tool registration + dispatch.

use std::io::Read;

use anyhow::{anyhow, Result};
use serde_json::{json, Value};

use crate::sandbox::{sandbox_check, sandbox_open_read};

/// Hard cap on the bytes either of `dock`/`dry_run`'s `receptor_path` or
/// `ligand_path` may consume.
///
/// Round-16 H1: pre-fix `dock_tool` and `dry_run_tool` read the receptor
/// and ligand files via bare `std::fs::read_to_string` — an MCP client
/// (which is LLM-untrusted input) could point either path at a
/// gigabytes-large file inside the sandbox root and OOM the MCP server.
/// 64 MiB is generous for chemistry (even a PDBQT for a multi-chain
/// receptor with explicit waters is in the low-MB range) while still
/// being well under any host's memory budget.
pub const MAX_PDBQT_FILE_BYTES: usize = 64 * 1024 * 1024;

/// Static `tools/list` response.
///
/// The docking tools are listed first, then the generative-design tool
/// set ([`crate::design::tool_list`]) is appended — one flat array of
/// every tool the server advertises.
pub fn list() -> Value {
    let mut tools = vec![
        json!({
            "name": "dock",
            "description": "Run a native AutoDock Vina docking job. Returns ranked poses.",
            "inputSchema": valenx_dock::schema::dock_tool_schema(),
        }),
        json!({
            "name": "dry_run",
            "description": "Pre-flight a docking job: validate inputs and return an execution plan.",
            "inputSchema": valenx_dock::schema::dry_run_tool_schema(),
        }),
        json!({
            "name": "describe",
            "description": "Return JSON Schema for DockConfig and Pose.",
            "inputSchema": { "type": "object" }
        }),
        json!({
            "name": "list_adapters",
            "description": "Return the full Valenx adapter catalog.",
            "inputSchema": { "type": "object" }
        }),
    ];
    // Append the generative-design (parametric CAD) tools.
    tools.extend(crate::design::tool_list());
    Value::Array(tools)
}

/// Dispatch a `tools/call`.
pub async fn call(params: &Value) -> Result<Value> {
    let name = params
        .get("name")
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow!("tools/call: missing `name`"))?;
    let args = params.get("arguments").cloned().unwrap_or(json!({}));
    // The generative-design tools are synchronous, in-process CAD
    // operations; route them first. `dispatch` returns `None` when the
    // name is not a design tool, so the docking tools below still match.
    if let Some(result) = crate::design::dispatch(name, &args) {
        return result;
    }
    match name {
        "dock" => dock_tool(&args).await,
        "dry_run" => dry_run_tool(&args).await,
        "describe" => Ok(describe_tool()),
        "list_adapters" => Ok(list_adapters_tool()),
        other => Err(anyhow!("unknown tool: {other}")),
    }
}

async fn dock_tool(args: &Value) -> Result<Value> {
    let receptor_raw = args["receptor_path"]
        .as_str()
        .ok_or_else(|| anyhow!("receptor_path"))?;
    let ligand_raw = args["ligand_path"]
        .as_str()
        .ok_or_else(|| anyhow!("ligand_path"))?;
    let output_raw = args["output_path"]
        .as_str()
        .ok_or_else(|| anyhow!("output_path"))?;
    // Output stays with the lexical-only `sandbox_check` — the write
    // path is the round-24 M4 dock runner's `OpenOptions + O_NOFOLLOW`
    // (round-25 H2) flow, so the TOCTOU surface is the same as the
    // read-side fix below but lives in `valenx-dock::write_output`.
    let output_path = sandbox_check(std::path::Path::new(output_raw))?;
    let cfg = parse_config(args)?;
    // Round-25 H1: open receptor + ligand via the TOCTOU-resistant
    // `sandbox_open_read` helper, which does the lexical check AND
    // passes `O_NOFOLLOW` / `FILE_FLAG_OPEN_REPARSE_POINT` to the
    // kernel so a leaf-symlink raced in between check and open fails
    // synchronously. Pre-fix `sandbox_check` + bare
    // `read_capped_to_string` (round-16 H1) was a TOCTOU pair: an
    // attacker who could write to the sandbox dir could swap a symlink
    // into the leaf path after the check returned and `read_to_string`
    // would silently follow it. The size cap from H1 still applies —
    // we `take(MAX_PDBQT_FILE_BYTES + 1)` against the handle so an
    // oversize read still surfaces as `InvalidData`.
    let receptor = read_handle_capped_to_string(
        sandbox_open_read(std::path::Path::new(receptor_raw))?,
        MAX_PDBQT_FILE_BYTES,
    )?;
    let ligand = read_handle_capped_to_string(
        sandbox_open_read(std::path::Path::new(ligand_raw))?,
        MAX_PDBQT_FILE_BYTES,
    )?;
    let poses = valenx_dock::dock(&receptor, &ligand, &cfg, &output_path, None)?;
    Ok(json!({
        "content": [{
            "type": "text",
            "text": format!("{} poses written to {}", poses.len(), output_path.display())
        }],
        "structuredContent": {
            "poses": poses.iter().enumerate().map(|(i, (_, score))| json!({
                "rank": i + 1,
                "score": score,
            })).collect::<Vec<_>>()
        }
    }))
}

async fn dry_run_tool(args: &Value) -> Result<Value> {
    let receptor_raw = args["receptor_path"]
        .as_str()
        .ok_or_else(|| anyhow!("receptor_path"))?;
    let ligand_raw = args["ligand_path"]
        .as_str()
        .ok_or_else(|| anyhow!("ligand_path"))?;
    let cfg = parse_config(args)?;
    // Round-25 H1: same TOCTOU-resistant read-open as `dock_tool`.
    // The dry-run tool reads the same receptor/ligand pair to build
    // an execution plan, so it inherits the exact same race window.
    let receptor = read_handle_capped_to_string(
        sandbox_open_read(std::path::Path::new(receptor_raw))?,
        MAX_PDBQT_FILE_BYTES,
    )?;
    let ligand = read_handle_capped_to_string(
        sandbox_open_read(std::path::Path::new(ligand_raw))?,
        MAX_PDBQT_FILE_BYTES,
    )?;
    let plan = valenx_dock::dock_dry_run(&receptor, &ligand, &cfg)?;
    let body = serde_json::to_value(&plan)?;
    Ok(json!({
        "content": [{ "type": "text", "text": format!("{body:#}") }],
        "structuredContent": body
    }))
}

/// Round-25 H1: read at most `cap` bytes from an already-open file
/// handle and decode as UTF-8. Mirrors
/// `valenx_core::io_caps::read_capped_to_string` but takes a `File`
/// instead of a path — necessary because the H1 fix opens the file
/// via `sandbox_open_read` (with `O_NOFOLLOW`) and then passes the
/// handle straight to the bounded read; re-opening by path here would
/// reintroduce the TOCTOU window. The `take(cap + 1)` + post-check
/// pattern matches the upstream helper: a file that exactly matched
/// the cap reads `cap` bytes; a file that outgrew the cap reads
/// `cap + 1` and we reject it.
fn read_handle_capped_to_string(file: std::fs::File, cap: usize) -> Result<String> {
    let mut buf = Vec::new();
    file.take(cap as u64 + 1).read_to_end(&mut buf)?;
    if buf.len() > cap {
        return Err(anyhow!(
            "MCP receptor/ligand exceeds {}-byte cap (read {} bytes)",
            cap,
            buf.len(),
        ));
    }
    String::from_utf8(buf).map_err(|e| anyhow!("MCP file not valid UTF-8: {e}"))
}

fn describe_tool() -> Value {
    json!({
        "content": [{ "type": "text", "text": "JSON Schema for DockConfig and Pose" }],
        "structuredContent": {
            "dock_config": valenx_dock::schema::dock_config_schema(),
            "pose": valenx_dock::schema::pose_schema(),
        }
    })
}

fn list_adapters_tool() -> Value {
    // The MCP crate intentionally avoids depending on valenx-app (which
    // would pull the GUI). Instead we instantiate each native-capable
    // adapter directly via its `info()` method and build the catalog
    // entry from [`valenx_core::AdapterDescriptor`]. Today the only
    // native adapter is Vina; the loop is structured so additional
    // native adapters can be appended one line at a time.
    use valenx_core::{Adapter, AdapterDescriptor};
    let adapters: Vec<Box<dyn Adapter>> = vec![Box::new(valenx_adapter_vina::VinaAdapter::new())];
    let descriptors: Vec<AdapterDescriptor> = adapters
        .iter()
        .map(|a| AdapterDescriptor::from_info(&a.info()))
        .collect();
    json!({
        "content": [{
            "type": "text",
            "text": format!("{} native adapters available", descriptors.len())
        }],
        "structuredContent": {
            "adapters": descriptors,
        }
    })
}

fn parse_config(args: &Value) -> Result<valenx_dock::DockConfig> {
    let mut cfg = valenx_dock::DockConfig::default();
    if let Some(c) = args.get("center").and_then(|v| v.as_array()) {
        if c.len() == 3 {
            cfg.center = nalgebra::Vector3::new(
                c[0].as_f64().unwrap_or(0.0),
                c[1].as_f64().unwrap_or(0.0),
                c[2].as_f64().unwrap_or(0.0),
            );
        }
    }
    if let Some(s) = args.get("size").and_then(|v| v.as_array()) {
        if s.len() == 3 {
            cfg.size = nalgebra::Vector3::new(
                s[0].as_f64().unwrap_or(20.0),
                s[1].as_f64().unwrap_or(20.0),
                s[2].as_f64().unwrap_or(20.0),
            );
        }
    }
    // u64→u32 narrowing checks: fail at the parse boundary so the
    // error mentions the field name. `DockConfig::validate()` also
    // catches these, but it can only complain about the final value.
    if let Some(v) = args.get("exhaustiveness").and_then(|v| v.as_u64()) {
        if v == 0 || v > 32 {
            return Err(anyhow!("exhaustiveness must be in 1..=32, got {v}"));
        }
        cfg.exhaustiveness = v as u32;
    }
    if let Some(v) = args.get("num_modes").and_then(|v| v.as_u64()) {
        if v == 0 || v > u32::MAX as u64 {
            return Err(anyhow!("num_modes must be in 1..=u32::MAX, got {v}"));
        }
        cfg.num_modes = v as u32;
    }
    if let Some(v) = args.get("energy_range").and_then(|v| v.as_f64()) {
        cfg.energy_range = v;
    }
    if let Some(v) = args.get("seed").and_then(|v| v.as_u64()) {
        cfg.seed = v;
    }
    Ok(cfg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn list_advertises_the_docking_and_design_tools() {
        let v = list();
        let arr = v.as_array().unwrap();
        let names: Vec<&str> = arr.iter().map(|t| t["name"].as_str().unwrap()).collect();
        // The four docking tools.
        assert!(names.contains(&"dock"));
        assert!(names.contains(&"dry_run"));
        assert!(names.contains(&"describe"));
        assert!(names.contains(&"list_adapters"));
        // ...and the generative-design tool set is appended after them.
        for design_tool in [
            "create_sketch",
            "add_sketch_line",
            "add_sketch_circle",
            "add_constraint",
            "pad",
            "pocket",
            "revolve",
            "fillet",
            "boolean",
            "evaluate_design",
            "export_design",
            "reset_design",
        ] {
            assert!(
                names.contains(&design_tool),
                "tool `{design_tool}` should be advertised"
            );
        }
        // Every advertised tool carries a non-empty inputSchema.
        for t in arr {
            assert!(
                t["inputSchema"].is_object(),
                "{} lacks an inputSchema",
                t["name"]
            );
        }
    }

    #[test]
    fn list_dock_schema_requires_paths() {
        let v = list();
        let dock = v
            .as_array()
            .unwrap()
            .iter()
            .find(|t| t["name"] == "dock")
            .unwrap();
        let req = dock["inputSchema"]["required"].as_array().unwrap();
        assert!(req.iter().any(|x| x == "receptor_path"));
        assert!(req.iter().any(|x| x == "ligand_path"));
        assert!(req.iter().any(|x| x == "output_path"));
    }

    #[test]
    fn list_dry_run_schema_requires_input_paths() {
        let v = list();
        let dry = v
            .as_array()
            .unwrap()
            .iter()
            .find(|t| t["name"] == "dry_run")
            .unwrap();
        let req = dry["inputSchema"]["required"].as_array().unwrap();
        assert!(req.iter().any(|x| x == "receptor_path"));
        assert!(req.iter().any(|x| x == "ligand_path"));
        assert!(
            !req.iter().any(|x| x == "output_path"),
            "dry_run must not require output_path"
        );
    }

    #[tokio::test]
    async fn describe_tool_returns_both_schemas() {
        let v = describe_tool();
        let sc = &v["structuredContent"];
        assert!(sc["dock_config"].is_object());
        assert!(sc["pose"].is_object());
    }

    #[test]
    fn parse_config_rejects_excess_exhaustiveness() {
        let args = json!({ "exhaustiveness": 64 });
        let err = parse_config(&args).unwrap_err();
        assert!(err.to_string().contains("1..=32"), "{err}");
    }

    #[test]
    fn parse_config_rejects_zero_exhaustiveness() {
        let args = json!({ "exhaustiveness": 0 });
        let err = parse_config(&args).unwrap_err();
        assert!(err.to_string().contains("1..=32"), "{err}");
    }

    #[test]
    fn parse_config_rejects_zero_num_modes() {
        let args = json!({ "num_modes": 0 });
        let err = parse_config(&args).unwrap_err();
        assert!(err.to_string().contains("num_modes"), "{err}");
    }

    #[test]
    fn parse_config_rejects_oversize_num_modes() {
        let args = json!({ "num_modes": (u32::MAX as u64) + 1 });
        let err = parse_config(&args).unwrap_err();
        assert!(err.to_string().contains("num_modes"), "{err}");
    }

    /// Serialise sandbox env-var-mutating tests so they can't race.
    fn sandbox_lock() -> std::sync::MutexGuard<'static, ()> {
        use std::sync::{Mutex, OnceLock};
        static M: OnceLock<Mutex<()>> = OnceLock::new();
        M.get_or_init(|| Mutex::new(())).lock().unwrap()
    }

    /// Round-16 H1 RED→GREEN: a receptor file larger than the cap is
    /// rejected before the docker reads it into memory. Pre-fix the
    /// bare `std::fs::read_to_string` would slurp the entire file
    /// regardless of size — a hostile MCP client (LLM-untrusted input)
    /// could OOM the server by writing a 5 GB receptor.pdbqt into the
    /// sandbox and then asking `dock` to use it.
    ///
    /// The `MutexGuard` is held across the `.await` deliberately —
    /// without it the env-var mutation (sandbox dir) races other
    /// tests in the file. Single-threaded tokio test runtime makes
    /// this safe (no executor stall).
    #[allow(clippy::await_holding_lock)]
    #[tokio::test]
    async fn dock_tool_rejects_oversize_receptor() {
        let _g = sandbox_lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-h1-dock-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(crate::sandbox::SANDBOX_ENV, &tmp);
        let receptor = tmp.join("big.pdbqt");
        // Write 1 byte over the cap so the stat-and-take path rejects.
        let oversize = vec![b'X'; MAX_PDBQT_FILE_BYTES + 1];
        std::fs::write(&receptor, &oversize).unwrap();
        let ligand = tmp.join("ligand.pdbqt");
        std::fs::write(&ligand, b"minimal\n").unwrap();
        let output = tmp.join("out.pdbqt");
        let args = json!({
            "receptor_path": receptor.to_string_lossy(),
            "ligand_path": ligand.to_string_lossy(),
            "output_path": output.to_string_lossy(),
        });
        let res = dock_tool(&args).await;
        assert!(
            res.is_err(),
            "expected oversize receptor read to be rejected, got Ok"
        );
        std::env::remove_var(crate::sandbox::SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-16 H1 RED→GREEN: same cap on the dry-run tool. Pre-fix
    /// `dry_run_tool` had the same bare `std::fs::read_to_string`
    /// pattern as `dock_tool` — closing the cap there too means an
    /// attacker can't OOM the server via the cheaper pre-flight tool.
    ///
    /// `#[allow(clippy::await_holding_lock)]`: see sibling test's
    /// rationale — env-var serialisation requires the lock; single-
    /// threaded tokio test runtime makes it safe.
    #[allow(clippy::await_holding_lock)]
    #[tokio::test]
    async fn dry_run_tool_rejects_oversize_ligand() {
        let _g = sandbox_lock();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mcp-h1-dry-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::env::set_var(crate::sandbox::SANDBOX_ENV, &tmp);
        let receptor = tmp.join("receptor.pdbqt");
        std::fs::write(&receptor, b"minimal\n").unwrap();
        let ligand = tmp.join("big.pdbqt");
        let oversize = vec![b'X'; MAX_PDBQT_FILE_BYTES + 1];
        std::fs::write(&ligand, &oversize).unwrap();
        let args = json!({
            "receptor_path": receptor.to_string_lossy(),
            "ligand_path": ligand.to_string_lossy(),
        });
        let res = dry_run_tool(&args).await;
        assert!(
            res.is_err(),
            "expected oversize ligand read to be rejected, got Ok"
        );
        std::env::remove_var(crate::sandbox::SANDBOX_ENV);
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
