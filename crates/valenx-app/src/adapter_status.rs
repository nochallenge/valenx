//! Status-pill colours + text for adapter rows in the browser. Each
//! [`AdapterStatus`] gets a (fill, text) colour pair and a short
//! human label.

use eframe::egui;

use valenx_core::AdapterStatus;

pub(crate) fn status_color(status: &AdapterStatus) -> (egui::Color32, egui::Color32) {
    match status {
        AdapterStatus::Ready { .. } => (
            egui::Color32::from_rgb(120, 220, 140),
            egui::Color32::from_rgb(220, 230, 230),
        ),
        AdapterStatus::Missing { .. } => (
            egui::Color32::from_rgb(200, 200, 200),
            egui::Color32::from_rgb(160, 160, 160),
        ),
        AdapterStatus::Outdated { .. } => (
            egui::Color32::from_rgb(240, 190, 110),
            egui::Color32::from_rgb(220, 220, 220),
        ),
        AdapterStatus::Broken { .. } => (
            egui::Color32::from_rgb(230, 120, 120),
            egui::Color32::from_rgb(220, 220, 220),
        ),
        AdapterStatus::Disabled => (
            egui::Color32::from_rgb(150, 150, 180),
            egui::Color32::from_rgb(160, 160, 160),
        ),
    }
}

pub(crate) fn status_label(status: &AdapterStatus) -> &'static str {
    match status {
        AdapterStatus::Ready { .. } => "Ready",
        AdapterStatus::Missing { .. } => "Missing",
        AdapterStatus::Outdated { .. } => "Outdated",
        AdapterStatus::Broken { .. } => "Broken",
        AdapterStatus::Disabled => "Disabled",
    }
}
