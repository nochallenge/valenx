//! Plot — title + axis labels + series + optional explicit ranges.

use serde::{Deserialize, Serialize};

use crate::series::Series;

/// Inclusive axis range. `Auto` derives `(min, max)` from the series
/// extents; `Manual(lo, hi)` overrides.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum AxisRange {
    /// Fit to series data.
    Auto,
    /// Explicit `(lo, hi)`.
    Manual(f64, f64),
}

impl Default for AxisRange {
    fn default() -> Self {
        Self::Auto
    }
}

/// A 2D plot — one or more [`Series`] sharing the same x / y axes.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Plot {
    /// Display title.
    pub title: String,
    /// X axis label.
    pub x_label: String,
    /// Y axis label.
    pub y_label: String,
    /// Series.
    pub series: Vec<Series>,
    /// X axis range.
    pub x_range: AxisRange,
    /// Y axis range.
    pub y_range: AxisRange,
}

impl Plot {
    /// Empty plot with the given title.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            ..Default::default()
        }
    }

    /// Append a series.
    pub fn add_series(&mut self, s: Series) {
        self.series.push(s);
    }

    /// Resolve `x_range` to concrete `(lo, hi)` against the series
    /// data. Empty plot returns `(0.0, 1.0)`.
    pub fn resolved_x_range(&self) -> (f64, f64) {
        match self.x_range {
            AxisRange::Manual(a, b) => (a.min(b), a.max(b)),
            AxisRange::Auto => {
                let (mut x0, mut x1) = (f64::INFINITY, f64::NEG_INFINITY);
                let mut any = false;
                for s in &self.series {
                    if let Some((sx0, sx1, _, _)) = s.extent() {
                        if sx0 < x0 {
                            x0 = sx0;
                        }
                        if sx1 > x1 {
                            x1 = sx1;
                        }
                        any = true;
                    }
                }
                if !any || !x0.is_finite() || !x1.is_finite() {
                    (0.0, 1.0)
                } else if (x1 - x0).abs() < f64::EPSILON {
                    (x0 - 0.5, x0 + 0.5)
                } else {
                    (x0, x1)
                }
            }
        }
    }

    /// Resolve `y_range` to concrete `(lo, hi)`.
    pub fn resolved_y_range(&self) -> (f64, f64) {
        match self.y_range {
            AxisRange::Manual(a, b) => (a.min(b), a.max(b)),
            AxisRange::Auto => {
                let (mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY);
                let mut any = false;
                for s in &self.series {
                    if let Some((_, _, sy0, sy1)) = s.extent() {
                        if sy0 < y0 {
                            y0 = sy0;
                        }
                        if sy1 > y1 {
                            y1 = sy1;
                        }
                        any = true;
                    }
                }
                if !any || !y0.is_finite() || !y1.is_finite() {
                    (0.0, 1.0)
                } else if (y1 - y0).abs() < f64::EPSILON {
                    (y0 - 0.5, y0 + 0.5)
                } else {
                    // Pad 5% so the data isn't flush with the axes.
                    let pad = (y1 - y0) * 0.05;
                    (y0 - pad, y1 + pad)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::series::SeriesStyle;

    #[test]
    fn auto_range_fits_series() {
        let mut p = Plot::new("test");
        let mut s = Series::new("a", SeriesStyle::Line);
        s.push(0.0, 0.0);
        s.push(2.0, 4.0);
        p.add_series(s);
        let (x0, x1) = p.resolved_x_range();
        assert!((x0 - 0.0).abs() < 1e-9 && (x1 - 2.0).abs() < 1e-9);
        let (y0, y1) = p.resolved_y_range();
        // Padded.
        assert!(y0 < 0.0 && y1 > 4.0);
    }

    #[test]
    fn manual_range_honoured() {
        let mut p = Plot::new("test");
        p.x_range = AxisRange::Manual(-10.0, 10.0);
        p.y_range = AxisRange::Manual(-1.0, 1.0);
        assert_eq!(p.resolved_x_range(), (-10.0, 10.0));
        assert_eq!(p.resolved_y_range(), (-1.0, 1.0));
    }
}
