//! JSON Schema documents for [`crate::DockConfig`] inputs and
//! [`crate::Pose`] outputs. Designed for LLM self-documentation:
//! a client hits [`dock_config_schema`] before calling [`crate::dock`]
//! to learn what fields exist and what their types/ranges are.
//!
//! We hand-write the schemas (vs. deriving) because schemars-derived
//! schemas for `nalgebra::Vector3` etc. don't quite match the
//! semantic types we want to expose ("3-element array of f64 in Å").

use serde_json::{json, Value};

/// JSON Schema for [`crate::DockConfig`].
pub fn dock_config_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "DockConfig",
        "type": "object",
        "required": ["center", "size"],
        "properties": {
            "center": {
                "type": "array",
                "description": "Centre of the search box (Å, receptor coords)",
                "items": { "type": "number" },
                "minItems": 3,
                "maxItems": 3
            },
            "size": {
                "type": "array",
                "description": "Edge lengths of the search box (Å)",
                "items": { "type": "number", "exclusiveMinimum": 0 },
                "minItems": 3,
                "maxItems": 3
            },
            "exhaustiveness": {
                "type": "integer",
                "description": "Number of independent ILS chains",
                "minimum": 1,
                "maximum": 32,
                "default": 8
            },
            "num_modes": {
                "type": "integer",
                "description": "Max poses to return",
                "minimum": 1,
                "default": 9
            },
            "energy_range": {
                "type": "number",
                "description": "kcal/mol cutoff above best pose",
                "exclusiveMinimum": 0,
                "default": 3.0
            },
            "seed": {
                "type": "integer",
                "description": "Reproducibility seed (u64)",
                "default": 0
            },
            "grid_spacing": {
                "type": "number",
                "description": "Grid spacing in Å",
                "exclusiveMinimum": 0,
                "default": 0.375
            }
        }
    })
}

/// JSON Schema for the MCP `dock` tool call — wraps [`dock_config_schema`]
/// with the receptor/ligand/output path fields the MCP handler requires.
/// All paths are sandboxed to `$VALENX_MCP_SANDBOX_DIR` (or a temp-dir
/// fallback if unset); see `valenx-mcp` crate docs.
pub fn dock_tool_schema() -> Value {
    let mut s = dock_config_schema();
    let props = s["properties"].as_object_mut().unwrap();
    props.insert(
        "receptor_path".into(),
        json!({
            "type": "string",
            "description": "Path to receptor PDBQT file (sandboxed under VALENX_MCP_SANDBOX_DIR)"
        }),
    );
    props.insert(
        "ligand_path".into(),
        json!({
            "type": "string",
            "description": "Path to ligand PDBQT file (sandboxed under VALENX_MCP_SANDBOX_DIR)"
        }),
    );
    props.insert(
        "output_path".into(),
        json!({
            "type": "string",
            "description": "Output PDBQT path (sandboxed under VALENX_MCP_SANDBOX_DIR)"
        }),
    );
    let req = s["required"].as_array_mut().unwrap();
    req.push(json!("receptor_path"));
    req.push(json!("ligand_path"));
    req.push(json!("output_path"));
    s
}

/// JSON Schema for the MCP `dry_run` tool — wraps [`dock_config_schema`]
/// with the receptor/ligand path fields (no `output_path`; dry-run
/// writes nothing). Paths are sandboxed; see [`dock_tool_schema`].
pub fn dry_run_tool_schema() -> Value {
    let mut s = dock_config_schema();
    let props = s["properties"].as_object_mut().unwrap();
    props.insert(
        "receptor_path".into(),
        json!({
            "type": "string",
            "description": "Path to receptor PDBQT file (sandboxed under VALENX_MCP_SANDBOX_DIR)"
        }),
    );
    props.insert(
        "ligand_path".into(),
        json!({
            "type": "string",
            "description": "Path to ligand PDBQT file (sandboxed under VALENX_MCP_SANDBOX_DIR)"
        }),
    );
    let req = s["required"].as_array_mut().unwrap();
    req.push(json!("receptor_path"));
    req.push(json!("ligand_path"));
    s
}

/// JSON Schema for one [`crate::Pose`] (as returned to the caller).
pub fn pose_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "Pose",
        "type": "object",
        "required": ["translation", "orientation_axis_angle", "torsions", "score"],
        "properties": {
            "translation": {
                "type": "array",
                "items": { "type": "number" },
                "minItems": 3,
                "maxItems": 3
            },
            "orientation_axis_angle": {
                "type": "array",
                "description": "Rodrigues vector (axis × angle); identity = [0,0,0]",
                "items": { "type": "number" },
                "minItems": 3,
                "maxItems": 3
            },
            "torsions": {
                "type": "array",
                "items": { "type": "number" }
            },
            "score": {
                "type": "number",
                "description": "Inter-molecular energy (kcal/mol)"
            },
            "rmsd_to_top": {
                "type": "number",
                "description": "Heavy-atom RMSD vs the top pose (0 for rank 1)"
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dock_config_schema_has_required_fields() {
        let s = dock_config_schema();
        assert_eq!(s["title"], "DockConfig");
        let req = s["required"].as_array().unwrap();
        assert!(req.iter().any(|v| v == "center"));
        assert!(req.iter().any(|v| v == "size"));
    }

    #[test]
    fn dock_config_schema_exhaustiveness_is_bounded() {
        let s = dock_config_schema();
        assert_eq!(s["properties"]["exhaustiveness"]["maximum"], 32);
    }

    #[test]
    fn pose_schema_has_score() {
        let s = pose_schema();
        assert!(s["properties"]["score"].is_object());
    }

    #[test]
    fn dock_tool_schema_requires_paths() {
        let s = dock_tool_schema();
        let req = s["required"].as_array().unwrap();
        assert!(req.iter().any(|v| v == "receptor_path"));
        assert!(req.iter().any(|v| v == "ligand_path"));
        assert!(req.iter().any(|v| v == "output_path"));
        let props = s["properties"].as_object().unwrap();
        assert!(props.contains_key("receptor_path"));
        assert!(props.contains_key("ligand_path"));
        assert!(props.contains_key("output_path"));
    }

    #[test]
    fn dry_run_tool_schema_requires_input_paths_only() {
        let s = dry_run_tool_schema();
        let req = s["required"].as_array().unwrap();
        assert!(req.iter().any(|v| v == "receptor_path"));
        assert!(req.iter().any(|v| v == "ligand_path"));
        assert!(
            !req.iter().any(|v| v == "output_path"),
            "dry_run must not require output_path"
        );
        let props = s["properties"].as_object().unwrap();
        assert!(!props.contains_key("output_path"));
    }
}
