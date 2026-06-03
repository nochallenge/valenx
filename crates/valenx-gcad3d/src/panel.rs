//! UI panel envelope.

use crate::arc::Arc3d;
use crate::line::Line3d;
use crate::plane::Plane3d;
use crate::surface::RuledSurface;
use crate::text::TextSolid;

/// Workbench panel state.
#[derive(Default)]
pub struct Gcad3dPanelState {
    /// Planes added by the user.
    pub planes: Vec<Plane3d>,
    /// Lines added.
    pub lines: Vec<Line3d>,
    /// Arcs added.
    pub arcs: Vec<Arc3d>,
    /// Ruled surfaces.
    pub surfaces: Vec<RuledSurface>,
    /// Extruded text solids.
    pub text_solids: Vec<TextSolid>,
    /// Status message.
    pub last_status: Option<String>,
    /// Error message.
    pub last_error: Option<String>,
}

impl Gcad3dPanelState {
    /// New empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// Push a plane.
    pub fn add_plane(&mut self, p: Plane3d) {
        self.planes.push(p);
        self.last_status = Some("added plane".into());
        self.last_error = None;
    }

    /// Push a line.
    pub fn add_line(&mut self, l: Line3d) {
        self.lines.push(l);
        self.last_status = Some("added line".into());
        self.last_error = None;
    }

    /// Push an arc.
    pub fn add_arc(&mut self, a: Arc3d) {
        self.arcs.push(a);
        self.last_status = Some("added arc".into());
        self.last_error = None;
    }

    /// Push a ruled surface.
    pub fn add_surface(&mut self, s: RuledSurface) {
        self.surfaces.push(s);
        self.last_status = Some("added ruled surface".into());
        self.last_error = None;
    }

    /// Push an extruded text solid.
    pub fn add_text(&mut self, t: TextSolid) {
        let n = t.glyphs.len();
        self.text_solids.push(t);
        self.last_status = Some(format!("added text ({n} glyphs)"));
        self.last_error = None;
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
