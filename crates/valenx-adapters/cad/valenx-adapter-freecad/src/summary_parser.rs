//! Parse the `summary.json` that the generated FreeCAD Python
//! script writes back out. Fields map 1:1 to what the script
//! emits; anything extra is tolerated so future script revisions
//! stay backwards-compatible.

use std::path::Path;

use serde::{Deserialize, Serialize};

use valenx_core::AdapterError;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct GeometrySummary {
    #[serde(default)]
    pub valenx_adapter: String,
    #[serde(default)]
    pub source: String,
    #[serde(default)]
    pub parts: Vec<PartSummary>,
    #[serde(default)]
    pub feature_tree: Vec<FeatureNode>,
    #[serde(default)]
    pub bounding_box: Option<BoundingBox>,
    #[serde(default)]
    pub volume: Option<f64>,
    #[serde(default)]
    pub area: Option<f64>,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct PartSummary {
    #[serde(default)]
    pub name: String,
    #[serde(default, rename = "type")]
    pub type_id: String,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct FeatureNode {
    #[serde(default)]
    pub name: String,
    #[serde(default, rename = "type")]
    pub type_id: String,
    #[serde(default)]
    pub visible: bool,
    #[serde(default)]
    pub has_shape: bool,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize)]
pub struct BoundingBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

impl BoundingBox {
    pub fn size(&self) -> [f64; 3] {
        [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ]
    }

    pub fn center(&self) -> [f64; 3] {
        [
            (self.max[0] + self.min[0]) * 0.5,
            (self.max[1] + self.min[1]) * 0.5,
            (self.max[2] + self.min[2]) * 0.5,
        ]
    }
}

pub fn parse_file(path: &Path) -> Result<GeometrySummary, AdapterError> {
    // Round-23 sweep: bound the summary read at
    // MAX_FREECAD_SUMMARY_BYTES (1 MiB) — sister to Cantera / MuJoCo
    // summary caps. FreeCAD summaries are tiny (bbox + area / volume /
    // centre-of-mass scalars); 1 MiB is generous.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_FREECAD_SUMMARY_BYTES as usize,
    )?;
    parse_str(&text, path)
}

pub fn parse_str(text: &str, for_path: &Path) -> Result<GeometrySummary, AdapterError> {
    serde_json::from_str::<GeometrySummary>(text).map_err(|e| AdapterError::ParseOutput {
        file: for_path.to_path_buf(),
        reason: format!("summary.json: {e}"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_happy_path() {
        let json = r#"{
          "valenx_adapter": "freecad",
          "source": "bracket.step",
          "parts": [{"name": "Body001", "type": "Part::Feature"}],
          "feature_tree": [
             {"name": "Body001", "type": "Part::Feature", "visible": true, "has_shape": true}
          ],
          "bounding_box": {"min": [0, 0, 0], "max": [1, 2, 3]},
          "volume": 6.0,
          "area": 22.0
        }"#;
        let s = parse_str(json, std::path::Path::new("summary.json")).unwrap();
        assert_eq!(s.parts.len(), 1);
        assert_eq!(s.parts[0].name, "Body001");
        let bb = s.bounding_box.unwrap();
        assert_eq!(bb.size(), [1.0, 2.0, 3.0]);
        assert_eq!(bb.center(), [0.5, 1.0, 1.5]);
        assert_eq!(s.volume.unwrap(), 6.0);
    }

    #[test]
    fn tolerates_missing_fields() {
        // A barely-filled summary — only `valenx_adapter` — should
        // still parse so partial FreeCAD runs don't break the flow.
        let json = r#"{ "valenx_adapter": "freecad" }"#;
        let s = parse_str(json, std::path::Path::new("summary.json")).unwrap();
        assert_eq!(s.valenx_adapter, "freecad");
        assert!(s.parts.is_empty());
        assert!(s.bounding_box.is_none());
    }

    #[test]
    fn malformed_json_is_parse_output_error() {
        let json = "not json at all";
        let err = parse_str(json, std::path::Path::new("summary.json")).unwrap_err();
        assert!(matches!(err, AdapterError::ParseOutput { .. }));
    }
}
