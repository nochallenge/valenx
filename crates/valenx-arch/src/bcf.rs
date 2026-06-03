//! BCF (BIM Collaboration Format) stub.
//!
//! BCF is an open OASIS/buildingSMART standard for exchanging
//! issue-tracker entries that reference parts of an IFC model. A real
//! BCF file is a ZIP archive containing a tree of XML files (one
//! topic per directory, plus image attachments + viewpoints).
//!
//! v1 ships an in-memory model and a directory-form writer that
//! emits one `markup.bcf` XML file per topic plus a top-level
//! `bcf.version` file. A future Phase 15.5 will pull a `zip` crate
//! into the workspace and bundle the directory into a true `.bcf` /
//! `.bcfzip` file.
//!
//! ## Schema
//!
//! v1 emits BCF 2.1 markup files:
//! - `bcf.version` — schema version metadata.
//! - `<topic-guid>/markup.bcf` — issue metadata.
//! - `<topic-guid>/viewpoint.bcfv` — viewpoint payload (camera +
//!   selection).
//!
//! This is enough for a downstream BCF viewer to import the
//! directory (most viewers accept either a `.bcf` file or a
//! directory).

use serde::{Deserialize, Serialize};

use crate::error::ArchError;

/// Open / Closed status of a BCF issue.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BcfStatus {
    /// Issue is open and awaiting action.
    Open,
    /// Issue is being worked on.
    InProgress,
    /// Issue is resolved.
    Resolved,
    /// Issue is closed (no longer requires action).
    Closed,
}

impl BcfStatus {
    /// Stable BCF 2.1 string for the `<TopicStatus>` element.
    pub fn label(self) -> &'static str {
        match self {
            BcfStatus::Open => "Open",
            BcfStatus::InProgress => "InProgress",
            BcfStatus::Resolved => "Resolved",
            BcfStatus::Closed => "Closed",
        }
    }
}

/// A 3D camera + selection — referenced from the issue's markup.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BcfViewpoint {
    /// Camera position.
    pub camera_position: [f64; 3],
    /// Camera direction (unit vector).
    pub camera_direction: [f64; 3],
    /// Camera up vector.
    pub camera_up: [f64; 3],
    /// Field of view in degrees.
    pub field_of_view_deg: f64,
    /// Optional list of IFC GUIDs that should be highlighted.
    pub selected_components: Vec<String>,
}

impl Default for BcfViewpoint {
    fn default() -> Self {
        Self {
            camera_position: [0.0, 0.0, 5.0],
            camera_direction: [0.0, 0.0, -1.0],
            camera_up: [0.0, 1.0, 0.0],
            field_of_view_deg: 60.0,
            selected_components: Vec::new(),
        }
    }
}

/// A single BCF issue.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BcfIssue {
    /// Stable issue identifier — used as the directory name. Use a
    /// UUID or short slug.
    pub id: String,
    /// Short title shown in the issue list.
    pub title: String,
    /// Long description.
    pub description: String,
    /// Status.
    pub status: BcfStatus,
    /// Optional list of viewpoints.
    pub viewpoints: Vec<BcfViewpoint>,
}

/// A BCF document — list of issues plus a schema version.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bcf {
    /// Detail-version of the BCF schema we emit (always `"2.1"` in
    /// v1).
    pub version: String,
    /// All issues in the project.
    pub issues: Vec<BcfIssue>,
}

impl Default for Bcf {
    /// Empty BCF at schema version `"2.1"`.
    fn default() -> Self {
        Self {
            version: "2.1".to_string(),
            issues: Vec::new(),
        }
    }
}

impl Bcf {
    /// Write the BCF as a directory tree at `dir`. Creates the
    /// directory if missing. Each issue becomes a sub-directory.
    ///
    /// v1 emits XML files only — no ZIP envelope. Downstream tools
    /// that need a `.bcf` file can `zip -r project.bcf .` on the
    /// emitted directory.
    pub fn write_to_directory(&self, dir: &std::path::Path) -> Result<(), ArchError> {
        std::fs::create_dir_all(dir)
            .map_err(|e| ArchError::BcfWriteFailed(format!("mkdir {}: {e}", dir.display())))?;
        // bcf.version.
        let version_xml = format!(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
             <Version VersionId=\"{v}\" />\n",
            v = escape_xml(&self.version)
        );
        valenx_core::io_caps::atomic_write_str(&dir.join("bcf.version"), &version_xml)
            .map_err(|e| ArchError::BcfWriteFailed(format!("write bcf.version: {e}")))?;

        for issue in &self.issues {
            if issue.id.is_empty() {
                return Err(ArchError::BcfWriteFailed("issue has empty id".into()));
            }
            // SECURITY: issue.id is deserialized from user-supplied BCF data.
            // Sanitize to a single safe path component so that crafted IDs
            // like "../../etc/passwd" cannot escape the target directory.
            let safe_id = sanitize_id(&issue.id);
            let topic_dir = dir.join(&safe_id);
            std::fs::create_dir_all(&topic_dir).map_err(|e| {
                ArchError::BcfWriteFailed(format!("mkdir {}: {e}", topic_dir.display()))
            })?;
            valenx_core::io_caps::atomic_write_str(&topic_dir.join("markup.bcf"), &markup_for(issue))
                .map_err(|e| ArchError::BcfWriteFailed(format!("write markup.bcf: {e}")))?;
            if let Some(vp) = issue.viewpoints.first() {
                valenx_core::io_caps::atomic_write_str(
                    &topic_dir.join("viewpoint.bcfv"),
                    &viewpoint_xml(vp),
                )
                .map_err(|e| ArchError::BcfWriteFailed(format!("write viewpoint.bcfv: {e}")))?;
            }
        }
        Ok(())
    }
}

fn markup_for(issue: &BcfIssue) -> String {
    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    s.push_str("<Markup>\n");
    s.push_str("  <Topic Guid=\"");
    s.push_str(&escape_xml(&issue.id));
    s.push_str("\" TopicStatus=\"");
    s.push_str(issue.status.label());
    s.push_str("\">\n");
    s.push_str("    <Title>");
    s.push_str(&escape_xml(&issue.title));
    s.push_str("</Title>\n");
    s.push_str("    <Description>");
    s.push_str(&escape_xml(&issue.description));
    s.push_str("</Description>\n");
    s.push_str("  </Topic>\n");
    if !issue.viewpoints.is_empty() {
        s.push_str("  <Viewpoints Guid=\"vp-0\">\n");
        s.push_str("    <Viewpoint>viewpoint.bcfv</Viewpoint>\n");
        s.push_str("  </Viewpoints>\n");
    }
    s.push_str("</Markup>\n");
    s
}

fn viewpoint_xml(vp: &BcfViewpoint) -> String {
    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    s.push_str("<VisualizationInfo>\n");
    s.push_str("  <PerspectiveCamera>\n");
    s.push_str(&format!(
        "    <CameraViewPoint><X>{}</X><Y>{}</Y><Z>{}</Z></CameraViewPoint>\n",
        vp.camera_position[0], vp.camera_position[1], vp.camera_position[2]
    ));
    s.push_str(&format!(
        "    <CameraDirection><X>{}</X><Y>{}</Y><Z>{}</Z></CameraDirection>\n",
        vp.camera_direction[0], vp.camera_direction[1], vp.camera_direction[2]
    ));
    s.push_str(&format!(
        "    <CameraUpVector><X>{}</X><Y>{}</Y><Z>{}</Z></CameraUpVector>\n",
        vp.camera_up[0], vp.camera_up[1], vp.camera_up[2]
    ));
    s.push_str(&format!(
        "    <FieldOfView>{}</FieldOfView>\n",
        vp.field_of_view_deg
    ));
    s.push_str("  </PerspectiveCamera>\n");
    if !vp.selected_components.is_empty() {
        s.push_str("  <Components>\n");
        s.push_str("    <Selection>\n");
        for guid in &vp.selected_components {
            s.push_str(&format!(
                "      <Component IfcGuid=\"{}\" />\n",
                escape_xml(guid)
            ));
        }
        s.push_str("    </Selection>\n");
        s.push_str("  </Components>\n");
    }
    s.push_str("</VisualizationInfo>\n");
    s
}

/// Sanitize a BCF issue identifier into a safe single-component directory
/// name. Any character outside the ASCII alphanumeric set + `-` + `_` is
/// replaced with `_`. This blocks path-traversal payloads such as
/// `../../etc/passwd`, embedded NUL bytes, and platform path separators.
fn sanitize_id(id: &str) -> String {
    let out: String = id
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if out.is_empty() {
        "_".to_string()
    } else {
        out
    }
}

fn escape_xml(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '&' => out.push_str("&amp;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            other => out.push(other),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_bcf_has_version() {
        let b = Bcf::default();
        assert_eq!(b.version, "2.1");
        assert!(b.issues.is_empty());
    }

    #[test]
    fn write_two_issues_to_directory() {
        let mut b = Bcf::default();
        b.issues.push(BcfIssue {
            id: "00000000-0000-4000-8000-000000000001".into(),
            title: "Door clash".into(),
            description: "Door overlaps with column.".into(),
            status: BcfStatus::Open,
            viewpoints: vec![BcfViewpoint::default()],
        });
        b.issues.push(BcfIssue {
            id: "00000000-0000-4000-8000-000000000002".into(),
            title: "Missing slab".into(),
            description: "Floor 2 has no slab modelled.".into(),
            status: BcfStatus::Resolved,
            viewpoints: vec![],
        });
        let dir = std::env::temp_dir().join("valenx_arch_bcf_test");
        let _ = std::fs::remove_dir_all(&dir);
        b.write_to_directory(&dir).unwrap();
        // Expect bcf.version and two topic dirs.
        assert!(dir.join("bcf.version").exists());
        assert!(dir
            .join("00000000-0000-4000-8000-000000000001")
            .join("markup.bcf")
            .exists());
        assert!(dir
            .join("00000000-0000-4000-8000-000000000001")
            .join("viewpoint.bcfv")
            .exists());
        assert!(dir
            .join("00000000-0000-4000-8000-000000000002")
            .join("markup.bcf")
            .exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn empty_id_rejected() {
        let mut b = Bcf::default();
        b.issues.push(BcfIssue {
            id: "".into(),
            title: "x".into(),
            description: "x".into(),
            status: BcfStatus::Open,
            viewpoints: vec![],
        });
        let dir = std::env::temp_dir().join("valenx_arch_bcf_bad_id");
        let _ = std::fs::remove_dir_all(&dir);
        assert!(b.write_to_directory(&dir).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn path_traversal_id_is_sanitized() {
        let mut b = Bcf::default();
        b.issues.push(BcfIssue {
            id: "../../etc/passwd".into(),
            title: "x".into(),
            description: "x".into(),
            status: BcfStatus::Open,
            viewpoints: vec![],
        });
        let dir = std::env::temp_dir().join("valenx_arch_bcf_traversal");
        let _ = std::fs::remove_dir_all(&dir);
        b.write_to_directory(&dir).unwrap();
        // Sanitized name: "..", "/", "\" all become "_".
        let expected = sanitize_id("../../etc/passwd");
        let topic_dir = dir.join(&expected);
        assert!(topic_dir.exists(), "expected {}", topic_dir.display());
        // Confirm the topic is a *direct* child of dir — i.e. it did not
        // escape.
        let canon_dir = std::fs::canonicalize(&dir).unwrap();
        let canon_topic = std::fs::canonicalize(&topic_dir).unwrap();
        assert!(
            canon_topic.starts_with(&canon_dir),
            "topic dir escaped: {} not under {}",
            canon_topic.display(),
            canon_dir.display()
        );
        // Sanity-check the sanitized name has no traversal chars.
        assert!(!expected.contains('/'));
        assert!(!expected.contains('\\'));
        assert!(!expected.contains(".."));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn sanitize_id_replaces_unsafe_chars() {
        assert_eq!(sanitize_id("foo/bar"), "foo_bar");
        assert_eq!(sanitize_id("../etc"), "___etc");
        assert_eq!(sanitize_id("foo\\bar"), "foo_bar");
        assert_eq!(sanitize_id("ok-id_123"), "ok-id_123");
        assert_eq!(sanitize_id("uuid:abc"), "uuid_abc");
        assert_eq!(sanitize_id(""), "_");
    }

    #[test]
    fn markup_xml_includes_title_and_status() {
        let issue = BcfIssue {
            id: "x".into(),
            title: "T1".into(),
            description: "D1".into(),
            status: BcfStatus::InProgress,
            viewpoints: vec![],
        };
        let x = markup_for(&issue);
        assert!(x.contains("T1"));
        assert!(x.contains("InProgress"));
        assert!(x.contains("D1"));
    }
}
