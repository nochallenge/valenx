//! SVG rendering. The PNG entry-point is a deferred stub.

use std::fmt::Write;

use crate::error::PlotError;
use crate::plot::Plot;
use crate::series::SeriesStyle;

/// Default colour palette used when a [`crate::Series`] doesn't carry
/// its own colour.
const PALETTE: &[&str] = &[
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
];

const MARGIN: f64 = 50.0;

/// Render the plot as a standalone SVG document.
///
/// `width` / `height` are the SVG viewport in user units (≈ pixels).
/// Errors:
/// - [`PlotError::Empty`] when the plot has no series.
/// - [`PlotError::BadParameter`] when `width` or `height` is non-
///   positive.
pub fn to_svg(plot: &Plot, width: u32, height: u32) -> Result<String, PlotError> {
    if plot.series.is_empty() {
        return Err(PlotError::Empty);
    }
    if width == 0 || height == 0 {
        return Err(PlotError::BadParameter {
            name: if width == 0 { "width" } else { "height" },
            reason: "must be > 0".into(),
        });
    }
    let w = width as f64;
    let h = height as f64;

    let (x0, x1) = plot.resolved_x_range();
    let (y0, y1) = plot.resolved_y_range();

    let plot_l = MARGIN;
    let plot_r = w - MARGIN * 0.6;
    let plot_t = MARGIN * 0.8;
    let plot_b = h - MARGIN;

    let mapx =
        |x: f64| -> f64 { plot_l + (x - x0) / (x1 - x0).max(f64::EPSILON) * (plot_r - plot_l) };
    let mapy =
        |y: f64| -> f64 { plot_b - (y - y0) / (y1 - y0).max(f64::EPSILON) * (plot_b - plot_t) };

    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(
        out,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">",
    );
    let _ = writeln!(
        out,
        "  <rect width=\"{w}\" height=\"{h}\" fill=\"#ffffff\"/>",
    );

    // Title.
    if !plot.title.is_empty() {
        let cx = w * 0.5;
        let title = escape_xml(&plot.title);
        let _ = writeln!(
            out,
            "  <text x=\"{cx}\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\" fill=\"#222\">{title}</text>",
        );
    }

    // Plot area.
    let rw = plot_r - plot_l;
    let rh = plot_b - plot_t;
    let _ = writeln!(
        out,
        "  <rect x=\"{plot_l}\" y=\"{plot_t}\" width=\"{rw}\" height=\"{rh}\" fill=\"#fafafa\" stroke=\"#888\" stroke-width=\"1\"/>",
    );

    // Ticks (5 each).
    for i in 0..=5 {
        let frac = i as f64 / 5.0;
        // X.
        let x_val = x0 + (x1 - x0) * frac;
        let xp = mapx(x_val);
        let xy2 = plot_b + 5.0;
        let _ = writeln!(
            out,
            "  <line x1=\"{xp}\" y1=\"{plot_b}\" x2=\"{xp}\" y2=\"{xy2}\" stroke=\"#bbb\"/>",
        );
        let xty = plot_b + 18.0;
        let _ = writeln!(
            out,
            "  <text x=\"{xp}\" y=\"{xty}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#444\">{x_val:.2}</text>",
        );
        // Y.
        let y_val = y0 + (y1 - y0) * frac;
        let yp = mapy(y_val);
        let ytx = plot_l - 5.0;
        let _ = writeln!(
            out,
            "  <line x1=\"{ytx}\" y1=\"{yp}\" x2=\"{plot_l}\" y2=\"{yp}\" stroke=\"#bbb\"/>",
        );
        let ytx2 = plot_l - 8.0;
        let _ = writeln!(
            out,
            "  <text x=\"{ytx2}\" y=\"{yp}\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#444\" dominant-baseline=\"middle\">{y_val:.2}</text>",
        );
    }

    // Axis labels.
    if !plot.x_label.is_empty() {
        let cx = (plot_l + plot_r) * 0.5;
        let ly = h - 8.0;
        let label = escape_xml(&plot.x_label);
        let _ = writeln!(
            out,
            "  <text x=\"{cx}\" y=\"{ly}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#222\">{label}</text>",
        );
    }
    if !plot.y_label.is_empty() {
        let cy = (plot_t + plot_b) * 0.5;
        let label = escape_xml(&plot.y_label);
        let _ = writeln!(
            out,
            "  <text x=\"14\" y=\"{cy}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#222\" transform=\"rotate(-90 14 {cy})\">{label}</text>",
        );
    }

    // Series.
    for (i, s) in plot.series.iter().enumerate() {
        let colour = s
            .colour
            .clone()
            .unwrap_or_else(|| PALETTE[i % PALETTE.len()].into());
        match s.style {
            SeriesStyle::Line | SeriesStyle::Area => {
                if s.data.is_empty() {
                    continue;
                }
                // Polyline points.
                let mut pts = String::new();
                for &(x, y) in &s.data {
                    let _ = write!(pts, "{},{} ", mapx(x), mapy(y));
                }
                if matches!(s.style, SeriesStyle::Area) {
                    // Close down to y = 0 (or bottom of range if 0 is
                    // outside it).
                    let baseline_y = if (y0..=y1).contains(&0.0) { 0.0 } else { y0 };
                    let bp = mapy(baseline_y);
                    let first_x = mapx(s.data[0].0);
                    let last_x = mapx(s.data.last().unwrap().0);
                    let mut poly = String::new();
                    let _ = write!(poly, "{first_x},{bp} ");
                    poly.push_str(&pts);
                    let _ = write!(poly, "{last_x},{bp} ");
                    let _ = writeln!(
                        out,
                        "  <polygon points=\"{poly}\" fill=\"{colour}\" fill-opacity=\"0.18\" stroke=\"none\"/>",
                    );
                }
                let _ = writeln!(
                    out,
                    "  <polyline points=\"{pts}\" fill=\"none\" stroke=\"{colour}\" stroke-width=\"1.8\"/>",
                );
            }
            SeriesStyle::Scatter => {
                for &(x, y) in &s.data {
                    let cx = mapx(x);
                    let cy = mapy(y);
                    let _ = writeln!(
                        out,
                        "  <circle cx=\"{cx}\" cy=\"{cy}\" r=\"3\" fill=\"{colour}\"/>",
                    );
                }
            }
            SeriesStyle::Bar => {
                let bar_w = if s.data.len() > 1 {
                    ((plot_r - plot_l) / s.data.len() as f64) * 0.7
                } else {
                    14.0
                };
                let baseline_y = if (y0..=y1).contains(&0.0) { 0.0 } else { y0 };
                for &(x, y) in &s.data {
                    let xp = mapx(x);
                    let y_top = mapy(y.max(baseline_y));
                    let y_bot = mapy(y.min(baseline_y));
                    let rx = xp - bar_w * 0.5;
                    let rh_local = (y_bot - y_top).max(0.5);
                    let _ = writeln!(
                        out,
                        "  <rect x=\"{rx}\" y=\"{y_top}\" width=\"{bar_w}\" height=\"{rh_local}\" fill=\"{colour}\" fill-opacity=\"0.7\" stroke=\"{colour}\"/>",
                    );
                }
            }
        }
    }

    // Legend top-right.
    let legend_x = plot_r - 110.0;
    let mut legend_y = plot_t + 8.0;
    for (i, s) in plot.series.iter().enumerate() {
        let colour = s
            .colour
            .clone()
            .unwrap_or_else(|| PALETTE[i % PALETTE.len()].into());
        let _ = writeln!(
            out,
            "  <rect x=\"{legend_x}\" y=\"{legend_y}\" width=\"10\" height=\"10\" fill=\"{colour}\"/>",
        );
        let tx = legend_x + 14.0;
        let ty = legend_y + 9.0;
        let n = escape_xml(&s.name);
        let _ = writeln!(
            out,
            "  <text x=\"{tx}\" y=\"{ty}\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#222\">{n}</text>",
        );
        legend_y += 14.0;
    }

    let _ = writeln!(out, "</svg>");
    Ok(out)
}

/// PNG rendering — deferred to a follow-up phase.
///
/// In v1 this always returns [`PlotError::PngDeferred`]. The desktop
/// shell can render the SVG via [`to_svg`] and let `egui` /
/// `resvg` rasterise it instead.
pub fn to_png(_plot: &Plot, _width: u32, _height: u32) -> Result<Vec<u8>, PlotError> {
    Err(PlotError::PngDeferred)
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
            _ => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::series::{Series, SeriesStyle};

    fn sample_plot() -> Plot {
        let mut p = Plot::new("Demo");
        p.x_label = "x".into();
        p.y_label = "y".into();
        let mut s = Series::new("sine", SeriesStyle::Line);
        for i in 0..32 {
            let x = (i as f64) / 31.0 * std::f64::consts::TAU;
            s.push(x, x.sin());
        }
        p.add_series(s);
        p
    }

    #[test]
    fn to_svg_basic() {
        let p = sample_plot();
        let svg = to_svg(&p, 480, 360).unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Demo"));
        assert!(svg.contains("polyline"));
    }

    #[test]
    fn to_svg_empty_errors() {
        let p = Plot::new("empty");
        assert!(matches!(to_svg(&p, 100, 100), Err(PlotError::Empty)));
    }

    #[test]
    fn to_svg_bad_width() {
        let p = sample_plot();
        assert!(matches!(
            to_svg(&p, 0, 100),
            Err(PlotError::BadParameter { name: "width", .. })
        ));
    }

    #[test]
    fn to_png_deferred() {
        let p = sample_plot();
        assert!(matches!(to_png(&p, 100, 100), Err(PlotError::PngDeferred)));
    }
}
