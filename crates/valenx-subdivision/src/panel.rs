//! UI panel-state envelope for the Subdivision workbench.

use crate::mesh::SubdivMesh;

/// Selected subdivision scheme.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Scheme {
    /// Catmull-Clark (general polygon → quads).
    CatmullClark,
    /// Loop (triangles only).
    Loop,
}

impl Default for Scheme {
    fn default() -> Self {
        Self::CatmullClark
    }
}

/// Workbench panel state.
pub struct SubdivPanelState {
    /// Source mesh.
    pub source: SubdivMesh,
    /// Subdivided mesh (last result).
    pub result: SubdivMesh,
    /// Current scheme.
    pub scheme: Scheme,
    /// Iteration count slider.
    pub iterations: u32,
    /// Top-menu entry label.
    pub menu_label: String,
    /// Status message.
    pub last_status: Option<String>,
    /// Error message.
    pub last_error: Option<String>,
}

impl Default for SubdivPanelState {
    fn default() -> Self {
        Self {
            source: SubdivMesh::new(),
            result: SubdivMesh::new(),
            scheme: Scheme::default(),
            iterations: 2,
            menu_label: "Modeling > Subdivision".into(),
            last_status: None,
            last_error: None,
        }
    }
}

impl SubdivPanelState {
    /// New empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the source mesh and reset the result.
    pub fn set_source(&mut self, m: SubdivMesh) {
        self.last_status = Some(format!(
            "loaded source: {} verts, {} faces",
            m.n_verts(),
            m.n_faces()
        ));
        self.last_error = None;
        self.result = m.clone();
        self.source = m;
    }

    /// Status setter.
    pub fn set_status(&mut self, s: impl Into<String>) {
        self.last_status = Some(s.into());
        self.last_error = None;
    }

    /// Error setter.
    pub fn set_error(&mut self, s: impl Into<String>) {
        self.last_error = Some(s.into());
        self.last_status = None;
    }
}
