//! UI envelope — operator palette.

/// Available operator kinds (matches the in-app op-palette).
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum BlenderOp {
    /// Extrude region.
    Extrude,
    /// Bevel edges.
    Bevel,
    /// Inset faces.
    Inset,
    /// Loop cut.
    LoopCut,
    /// Bridge edge loops.
    Bridge,
    /// Boolean modifier.
    Boolean,
    /// Solidify.
    Solidify,
}

impl BlenderOp {
    /// All seven op kinds — drives the palette layout.
    pub fn all() -> &'static [BlenderOp] {
        &[
            Self::Extrude,
            Self::Bevel,
            Self::Inset,
            Self::LoopCut,
            Self::Bridge,
            Self::Boolean,
            Self::Solidify,
        ]
    }

    /// Hot-key matching Blender (free-form; informational).
    pub fn hotkey(self) -> &'static str {
        match self {
            Self::Extrude => "E",
            Self::Bevel => "Ctrl+B",
            Self::Inset => "I",
            Self::LoopCut => "Ctrl+R",
            Self::Bridge => "Ctrl+E > Bridge",
            Self::Boolean => "Modifier > Boolean",
            Self::Solidify => "Modifier > Solidify",
        }
    }
}

/// Workbench panel state.
#[derive(Default)]
pub struct BlenderOpPanelState {
    /// Selected operator.
    pub selected: Option<BlenderOp>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl BlenderOpPanelState {
    /// New empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Select an op.
    pub fn select(&mut self, op: BlenderOp) {
        self.selected = Some(op);
        self.last_status = Some(format!("selected `{op:?}` ({})", op.hotkey()));
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
