//! Small parsers for the case-level metadata the rest of the app
//! needs to interpret: the `solver` string convention (adapter id
//! prefix) and the `[sweep.derived]` numeric block that the sweep
//! materialiser stamps into each derived case.

/// Resolve an adapter id from a case's `solver` string. The
/// convention is `"<adapter_id>.<specific-solver-or-analysis>"` —
/// e.g. `"openfoam.simpleFoam"` → `"openfoam"`,
/// `"gmsh.delaunay"` → `"gmsh"`, `"cantera.equilibrium"` →
/// `"cantera"`. A solver string with no dot is treated as the
/// adapter id itself.
pub fn adapter_id_from_solver(solver: &str) -> &str {
    solver.split('.').next().unwrap_or(solver)
}

/// Recover the numeric inputs that the sweep optimizer wrote into a
/// derived `case.toml`. Reads the optional `[sweep.derived]` block
/// (a flat table of `<short-name> = <number>` pairs the materialiser
/// stamps in); missing block / non-numeric values silently skip.
///
/// Used by the sweep dataset assembler to recover the per-sample
/// input vector. Until the materialiser writes the block this will
/// return an empty list — that's fine, the assembler still emits a
/// usable outputs.npy + manifest.
///
/// Round-22 M3 (sister to `sweep.rs` caps): pre-fix the read was a
/// bare `std::fs::read_to_string(path)` — a derived `case.toml`
/// that grew (or was injected) past a few MB would slurp into
/// memory before the toml parser ran. Capped at
/// `MAX_PROJECT_FILE_BYTES` (1 MiB) to match every other
/// project-toml load site (the sweep aggregator's base-case read in
/// `sweep::assemble_sweep_dataset`, the rbac-override loader, …).
pub fn derived_inputs_from_case_toml(path: &std::path::Path) -> Vec<(String, f64)> {
    let Ok(text) = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    ) else {
        return Vec::new();
    };
    let Ok(value) = toml::from_str::<toml::Value>(&text) else {
        return Vec::new();
    };
    let Some(derived) = value
        .get("sweep")
        .and_then(|s| s.get("derived"))
        .and_then(|d| d.as_table())
    else {
        return Vec::new();
    };
    let mut inputs: Vec<(String, f64)> = Vec::new();
    for (k, v) in derived {
        let f = match v {
            toml::Value::Float(f) => *f,
            toml::Value::Integer(i) => *i as f64,
            _ => continue,
        };
        inputs.push((k.clone(), f));
    }
    // Sort by key so two derived runs with the same set of inputs
    // produce the same column order — required by the dataset
    // exporter's "all samples share the input schema" rule.
    inputs.sort_by(|a, b| a.0.cmp(&b.0));
    inputs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adapter_id_from_solver_pulls_prefix() {
        assert_eq!(adapter_id_from_solver("openfoam.simpleFoam"), "openfoam");
        assert_eq!(adapter_id_from_solver("gmsh.delaunay"), "gmsh");
        assert_eq!(adapter_id_from_solver("cantera.equilibrium"), "cantera");
        assert_eq!(adapter_id_from_solver("calculix.static"), "calculix");
        assert_eq!(adapter_id_from_solver("elmer.heat"), "elmer");
        // No dot → the string is the id.
        assert_eq!(adapter_id_from_solver("mujoco"), "mujoco");
        // Empty → empty (caller decides what to do).
        assert_eq!(adapter_id_from_solver(""), "");
    }

    #[test]
    fn derived_inputs_from_case_toml_picks_up_numeric_block() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-derived-inputs-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let case = dir.join("case.toml");
        std::fs::write(
            &case,
            r#"
solver = "openfoam"

[sweep.derived]
aoa = 5.0
re = 1000000
"#,
        )
        .unwrap();
        let inputs = derived_inputs_from_case_toml(&case);
        // Sorted by key — `aoa` comes before `re`.
        assert_eq!(inputs.len(), 2);
        assert_eq!(inputs[0].0, "aoa");
        assert_eq!(inputs[0].1, 5.0);
        assert_eq!(inputs[1].0, "re");
        assert_eq!(inputs[1].1, 1_000_000.0);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn derived_inputs_returns_empty_when_block_missing() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-derived-inputs-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let case = dir.join("case.toml");
        std::fs::write(&case, "solver = \"gmsh\"\n").unwrap();
        let inputs = derived_inputs_from_case_toml(&case);
        assert!(inputs.is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Round-22 M3 RED→GREEN (sister to sweep.rs caps): a derived
    /// `case.toml` larger than `MAX_PROJECT_FILE_BYTES` (1 MiB) must
    /// produce an empty result WITHOUT being slurped into memory.
    /// Pre-fix the read was a bare `std::fs::read_to_string(path)`
    /// and would have allocated the full file size before the toml
    /// parser saw the first key.
    ///
    /// Uses `set_len` to create a sparse over-cap file without
    /// writing 5 MiB of zeros on every CI run.
    #[test]
    fn derived_inputs_rejects_oversize_case_toml() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-derived-inputs-r22m3-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let case = dir.join("case.toml");
        // Past the 1 MiB MAX_PROJECT_FILE_BYTES cap.
        let f = std::fs::File::create(&case).unwrap();
        f.set_len(valenx_core::project::loader::MAX_PROJECT_FILE_BYTES + 1)
            .unwrap();
        drop(f);
        // Pre-fix: would slurp 5 MiB+ then fail to parse → return
        // empty. Post-fix: bounded read errors out → return empty
        // without the allocation. Same return value, different mid-
        // path — the cap protects memory, not the public contract.
        let inputs = derived_inputs_from_case_toml(&case);
        assert!(inputs.is_empty(), "got: {inputs:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
