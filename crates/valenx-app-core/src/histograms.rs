//! Mesh-quality histogram rendering — aspect-ratio + skewness rows
//! with Unicode block bars proportional to per-bucket counts.

use eframe::egui;

/// Render an aspect-ratio histogram as text rows with Unicode block
/// bars sized proportionally to the per-bucket count. Rows look like:
/// `≤ 1.50: ████████ 1234`. Uncategorised + overflow get their own
/// summary line at the bottom when non-zero.
pub fn render_aspect_histogram(ui: &mut egui::Ui, hist: &valenx_mesh::AspectRatioHistogram) {
    ui.label("Aspect ratio");
    let max_count = hist
        .counts
        .iter()
        .copied()
        .chain(std::iter::once(hist.overflow))
        .max()
        .unwrap_or(0);
    for (i, &edge) in hist.buckets.iter().enumerate() {
        let count = hist.counts[i];
        let bar = histogram_bar(count, max_count, 16);
        ui.label(egui::RichText::new(format!("≤ {edge:>6.2}: {bar} {count}")).monospace());
    }
    if hist.overflow > 0 {
        let bar = histogram_bar(hist.overflow, max_count, 16);
        ui.label(egui::RichText::new(format!("> max  : {bar} {n}", n = hist.overflow)).monospace());
    }
    if hist.uncategorised > 0 {
        ui.label(format!("uncategorised: {}", hist.uncategorised));
    }
}

/// Render a skewness histogram as text rows with quality-band labels
/// (excellent / good / acceptable / poor / very poor).
pub fn render_skewness_histogram(ui: &mut egui::Ui, hist: &valenx_mesh::SkewnessHistogram) {
    ui.label("Skewness");
    let labels = ["excellent", "good", "acceptable", "poor", "very poor"];
    let max_count = hist.counts.iter().copied().max().unwrap_or(0);
    for (i, &edge) in hist.buckets.iter().enumerate() {
        let count = hist.counts[i];
        let bar = histogram_bar(count, max_count, 16);
        let label = labels.get(i).copied().unwrap_or("");
        ui.label(
            egui::RichText::new(format!("≤ {edge:.2} ({label:<11}): {bar} {count}")).monospace(),
        );
    }
    if hist.uncategorised > 0 {
        ui.label(format!("uncategorised: {}", hist.uncategorised));
    }
}

/// Render `count` as a Unicode block-bar of width up to `max_width`
/// proportional to `max`. Always returns at least one cell when
/// count > 0 so non-zero buckets remain visible.
pub fn histogram_bar(count: u64, max: u64, max_width: usize) -> String {
    if max == 0 || count == 0 {
        return " ".repeat(max_width);
    }
    let n = ((count as f64 / max as f64) * max_width as f64).ceil() as usize;
    let n = n.max(1).min(max_width);
    let bar: String = "█".repeat(n);
    let pad: String = " ".repeat(max_width - n);
    format!("{bar}{pad}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_bar_zero_count_is_all_spaces() {
        let s = histogram_bar(0, 100, 8);
        assert_eq!(s, "        ");
    }

    #[test]
    fn histogram_bar_zero_max_is_all_spaces() {
        let s = histogram_bar(5, 0, 8);
        assert_eq!(s, "        ");
    }

    #[test]
    fn histogram_bar_count_equals_max_fills_width() {
        let s = histogram_bar(10, 10, 4);
        assert_eq!(s, "████");
    }

    #[test]
    fn histogram_bar_small_nonzero_count_shows_at_least_one_block() {
        // count/max = 1/100; ceil(1/100 * 8) = 1 — but the .max(1)
        // floor guarantees we never drop a non-zero bucket.
        let s = histogram_bar(1, 100, 8);
        assert!(s.starts_with("█"), "expected leading block, got {s:?}");
        // Exactly 1 block + 7 spaces.
        let blocks = s.chars().filter(|c| *c == '█').count();
        assert_eq!(blocks, 1);
        assert_eq!(s.chars().filter(|c| *c == ' ').count(), 7);
    }

    #[test]
    fn histogram_bar_proportional_scaling_at_half() {
        // count=50/100 -> ceil(0.5 * 8) = 4 blocks.
        let s = histogram_bar(50, 100, 8);
        let blocks = s.chars().filter(|c| *c == '█').count();
        assert_eq!(blocks, 4);
    }
}
