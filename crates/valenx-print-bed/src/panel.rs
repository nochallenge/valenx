//! UI state envelope for the Print Bed Layout panel.

use crate::printer::Printer;

/// Workbench-panel state for print-bed layout.
#[derive(Clone, Debug)]
pub struct PrintBedPanelState {
    /// Active printer.
    pub printer: Printer,
    /// Number of parts currently on the bed.
    pub last_parts: usize,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl Default for PrintBedPanelState {
    fn default() -> Self {
        Self {
            printer: Printer::new(
                (220.0, 220.0, 250.0),
                crate::printer::BedType::Heated,
                crate::printer::BedMaterial::Pei,
            ),
            last_parts: 0,
            last_status: None,
            last_error: None,
        }
    }
}

impl PrintBedPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record success.
    pub fn set_status(&mut self, msg: impl Into<String>, parts: usize) {
        self.last_status = Some(msg.into());
        self.last_error = None;
        self.last_parts = parts;
    }

    /// Record failure.
    pub fn set_error(&mut self, msg: impl Into<String>) {
        self.last_error = Some(msg.into());
        self.last_status = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_printer_is_220_bed() {
        let s = PrintBedPanelState::new();
        assert!((s.printer.bed_size.0 - 220.0).abs() < 1e-9);
    }
}
