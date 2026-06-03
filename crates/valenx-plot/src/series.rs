//! Series — one labelled line/scatter/bar/area of `(x, y)` data points.

use serde::{Deserialize, Serialize};

/// How a series is drawn on the plot canvas.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum SeriesStyle {
    /// Connected polyline.
    #[default]
    Line,
    /// Discrete circular markers.
    Scatter,
    /// Vertical bars from y=0 (or the bottom of the y axis range) up
    /// to each data point.
    Bar,
    /// Polyline plus a translucent fill down to y=0.
    Area,
}

impl SeriesStyle {
    /// UI dropdown label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Line => "Line",
            Self::Scatter => "Scatter",
            Self::Bar => "Bar",
            Self::Area => "Area",
        }
    }
}

/// One named series of `(x, y)` data points.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Series {
    /// Display name (legend + tooltip).
    pub name: String,
    /// Data points sorted by x for Line / Area / Bar to make sense.
    pub data: Vec<(f64, f64)>,
    /// Draw style.
    pub style: SeriesStyle,
    /// Stroke / marker / fill colour as `#RRGGBB`. `None` → engine
    /// picks from a default palette by series index.
    pub colour: Option<String>,
}

impl Series {
    /// Construct a series.
    pub fn new(name: impl Into<String>, style: SeriesStyle) -> Self {
        Self {
            name: name.into(),
            data: Vec::new(),
            style,
            colour: None,
        }
    }

    /// Push a point.
    pub fn push(&mut self, x: f64, y: f64) {
        self.data.push((x, y));
    }

    /// Inclusive `(min_x, max_x, min_y, max_y)` over `data`. Returns
    /// `None` for an empty series.
    pub fn extent(&self) -> Option<(f64, f64, f64, f64)> {
        let (mut x0, mut x1) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY);
        if self.data.is_empty() {
            return None;
        }
        for &(x, y) in &self.data {
            if x < x0 {
                x0 = x;
            }
            if x > x1 {
                x1 = x;
            }
            if y < y0 {
                y0 = y;
            }
            if y > y1 {
                y1 = y;
            }
        }
        Some((x0, x1, y0, y1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extent_basic() {
        let mut s = Series::new("a", SeriesStyle::Line);
        s.push(0.0, 1.0);
        s.push(1.0, 5.0);
        s.push(2.0, -2.0);
        let e = s.extent().unwrap();
        assert_eq!(e, (0.0, 2.0, -2.0, 5.0));
    }

    #[test]
    fn extent_empty() {
        let s = Series::new("a", SeriesStyle::Line);
        assert!(s.extent().is_none());
    }

    #[test]
    fn style_labels() {
        assert_eq!(SeriesStyle::Bar.label(), "Bar");
    }
}
