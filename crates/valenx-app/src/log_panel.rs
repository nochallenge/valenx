//! Log panel — the second tab in the bottom dock, alongside
//! residuals. Shows every `RunEvent::LogLine` the adapter streamed,
//! colour-coded by level, with an autoscroll toggle and per-level
//! filters so noisy logs stay legible.
//!
//! The log is a ring buffer capped at [`MAX_LOG_LINES`] so a
//! long-running solver can't blow up the UI's memory. When the cap is
//! reached the oldest lines are dropped from the front.

use std::collections::VecDeque;

use eframe::egui;
use valenx_core::LogLevel;

/// How many log lines to retain. A 500-iteration simpleFoam run with
/// ~8 solver lines per iteration + banner text is well under this.
pub const MAX_LOG_LINES: usize = 20_000;

/// One captured log line.
#[derive(Clone, Debug)]
pub struct LogLine {
    pub level: LogLevel,
    pub text: String,
}

/// In-memory log state. `Default` gives an empty log with autoscroll
/// on and every level visible.
#[derive(Debug)]
pub struct LogPanel {
    pub lines: VecDeque<LogLine>,
    pub autoscroll: bool,
    pub show_trace: bool,
    pub show_debug: bool,
    pub show_info: bool,
    pub show_warn: bool,
    pub show_error: bool,
}

impl Default for LogPanel {
    fn default() -> Self {
        Self {
            lines: VecDeque::with_capacity(1024),
            autoscroll: true,
            // Trace/Debug off by default — these are overwhelming on
            // a real solver; Info + above is the common case.
            show_trace: false,
            show_debug: false,
            show_info: true,
            show_warn: true,
            show_error: true,
        }
    }
}

impl LogPanel {
    /// Append a line, dropping the oldest if we've hit the cap.
    pub fn push(&mut self, level: LogLevel, text: String) {
        if self.lines.len() >= MAX_LOG_LINES {
            self.lines.pop_front();
        }
        self.lines.push_back(LogLine { level, text });
    }

    /// Drop every buffered log line. Filter switches are preserved.
    pub fn clear(&mut self) {
        self.lines.clear();
    }

    /// Is this level visible given the current filter switches?
    pub fn level_visible(&self, level: LogLevel) -> bool {
        match level {
            LogLevel::Trace => self.show_trace,
            LogLevel::Debug => self.show_debug,
            LogLevel::Info => self.show_info,
            LogLevel::Warn => self.show_warn,
            LogLevel::Error => self.show_error,
        }
    }

    /// Number of lines currently passing the filter.
    pub fn filtered_count(&self) -> usize {
        self.lines
            .iter()
            .filter(|l| self.level_visible(l.level))
            .count()
    }
}

/// Header row with level-filter checkboxes + autoscroll + clear.
pub fn header(ui: &mut egui::Ui, panel: &mut LogPanel) {
    ui.checkbox(&mut panel.show_info, "info");
    ui.checkbox(&mut panel.show_warn, "warn");
    ui.checkbox(&mut panel.show_error, "error");
    ui.separator();
    ui.checkbox(&mut panel.show_debug, "debug");
    ui.checkbox(&mut panel.show_trace, "trace");
    ui.separator();
    ui.checkbox(&mut panel.autoscroll, "auto-scroll");
    if ui.button("Clear").clicked() {
        panel.clear();
    }
}

/// Render the log body (scrollable list).
pub fn show(ui: &mut egui::Ui, panel: &LogPanel) {
    if panel.lines.is_empty() {
        ui.vertical_centered(|ui| {
            ui.add_space(8.0);
            ui.label("Log lines from running solvers will appear here.");
            ui.add_space(8.0);
        });
        return;
    }

    let mut scroll = egui::ScrollArea::vertical().auto_shrink([false, false]);
    if panel.autoscroll {
        scroll = scroll.stick_to_bottom(true);
    }
    scroll.show(ui, |ui| {
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        let font = egui::FontId::monospace(12.0);
        for line in &panel.lines {
            if !panel.level_visible(line.level) {
                continue;
            }
            let color = color_for(line.level);
            ui.colored_label(color, egui::RichText::new(&line.text).font(font.clone()));
        }
    });
}

fn color_for(level: LogLevel) -> egui::Color32 {
    match level {
        LogLevel::Trace => egui::Color32::from_gray(120),
        LogLevel::Debug => egui::Color32::from_gray(170),
        LogLevel::Info => egui::Color32::from_rgb(200, 210, 220),
        LogLevel::Warn => egui::Color32::from_rgb(240, 200, 110),
        LogLevel::Error => egui::Color32::from_rgb(240, 130, 130),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_respects_cap() {
        let mut panel = LogPanel::default();
        for i in 0..(MAX_LOG_LINES + 100) {
            panel.push(LogLevel::Info, format!("line {i}"));
        }
        assert_eq!(panel.lines.len(), MAX_LOG_LINES);
        // Oldest 100 lines should have been dropped.
        assert_eq!(panel.lines.front().unwrap().text, "line 100");
    }

    #[test]
    fn filters_honoured() {
        let mut panel = LogPanel::default();
        panel.push(LogLevel::Info, "info line".into());
        panel.push(LogLevel::Warn, "warn line".into());
        panel.push(LogLevel::Error, "error line".into());
        assert_eq!(panel.filtered_count(), 3);
        panel.show_info = false;
        assert_eq!(panel.filtered_count(), 2);
        panel.show_warn = false;
        panel.show_error = false;
        assert_eq!(panel.filtered_count(), 0);
    }

    #[test]
    fn clear_empties_buffer() {
        let mut panel = LogPanel::default();
        panel.push(LogLevel::Info, "a".into());
        panel.push(LogLevel::Info, "b".into());
        panel.clear();
        assert_eq!(panel.lines.len(), 0);
    }
}
