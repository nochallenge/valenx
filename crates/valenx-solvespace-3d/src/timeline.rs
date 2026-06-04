//! Parametric feature timeline — Fusion-style feature history.
//!
//! A [`FeatureTimeline`] is an ordered list of modelling [`Feature`]s whose
//! dimensions are **parameter expressions**, resolved against a
//! [`ParameterTable`]. [`FeatureTimeline::rebuild`] resolves every feature's
//! expressions and produces the solids; editing a parameter and rebuilding
//! re-drives the whole model — the parametric-history loop at the heart of
//! Fusion / SolidWorks. The geometry kernel is `valenx-cad` (box / cylinder
//! primitives + profile extrusion).

use valenx_cad::{box_solid, cylinder, prism, CadError, Solid};

use crate::parameters::{ParamError, ParameterTable};

/// One modelling operation. Dimensions are parameter-expression strings
/// (e.g. `"width"`, `"base * 2"`, `"40"`) resolved against a [`ParameterTable`].
pub enum Feature {
    /// An axis-aligned box `dx × dy × dz`.
    Box {
        /// Width expression (x).
        dx: String,
        /// Depth expression (y).
        dy: String,
        /// Height expression (z).
        dz: String,
    },
    /// A cylinder of the given radius and height.
    Cylinder {
        /// Radius expression.
        radius: String,
        /// Height expression.
        height: String,
    },
    /// A prism: a fixed 2-D profile extruded to the given height.
    Extrude {
        /// Closed 2-D profile, as `(x, y)` points.
        profile: Vec<(f64, f64)>,
        /// Extrusion-height expression.
        height: String,
    },
}

/// An ordered list of parametric features.
#[derive(Default)]
pub struct FeatureTimeline {
    /// The features, in build order.
    pub features: Vec<Feature>,
}

/// A failure while rebuilding the timeline.
#[derive(Debug)]
pub enum TimelineError {
    /// A dimension expression failed to resolve.
    Param(ParamError),
    /// A geometry-kernel operation failed (e.g. a degenerate dimension).
    Cad(CadError),
}

impl std::fmt::Display for TimelineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TimelineError::Param(e) => write!(f, "parameter: {e}"),
            TimelineError::Cad(e) => write!(f, "geometry: {e}"),
        }
    }
}

impl std::error::Error for TimelineError {}

impl FeatureTimeline {
    /// An empty timeline.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a feature.
    pub fn push(&mut self, feature: Feature) {
        self.features.push(feature);
    }

    /// Number of features.
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// Whether the timeline has no features.
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// Rebuild every feature against `params`, returning one [`Solid`] per
    /// feature. Editing a parameter and rebuilding re-drives the model.
    pub fn rebuild(&self, params: &ParameterTable) -> Result<Vec<Solid>, TimelineError> {
        let resolve = |expr: &str| params.compute(expr).map_err(TimelineError::Param);
        let mut solids = Vec::with_capacity(self.features.len());
        for feature in &self.features {
            let solid = match feature {
                Feature::Box { dx, dy, dz } => {
                    box_solid(resolve(dx)?, resolve(dy)?, resolve(dz)?).map_err(TimelineError::Cad)?
                }
                Feature::Cylinder { radius, height } => {
                    cylinder(resolve(radius)?, resolve(height)?).map_err(TimelineError::Cad)?
                }
                Feature::Extrude { profile, height } => {
                    prism(profile, resolve(height)?).map_err(TimelineError::Cad)?
                }
            };
            solids.push(solid);
        }
        Ok(solids)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> ParameterTable {
        let mut p = ParameterTable::new();
        p.set("base", "10");
        p.set("width", "base * 3"); // 30
        p.set("height", "base * 2"); // 20
        p
    }

    #[test]
    fn rebuilds_parametric_features_into_solids() {
        let mut tl = FeatureTimeline::new();
        tl.push(Feature::Box { dx: "width".into(), dy: "height".into(), dz: "base".into() });
        tl.push(Feature::Cylinder { radius: "base / 2".into(), height: "height".into() });
        let solids = tl.rebuild(&params()).expect("rebuild");
        assert_eq!(solids.len(), 2);
        assert_eq!(solids[0].faces(), 6, "a box has six faces");
    }

    #[test]
    fn extrudes_a_profile() {
        let mut tl = FeatureTimeline::new();
        // A unit-square profile extruded to a parameter-driven height.
        tl.push(Feature::Extrude {
            profile: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
            height: "base".into(),
        });
        let solids = tl.rebuild(&params()).expect("extrude");
        assert_eq!(solids.len(), 1);
        assert!(solids[0].faces() >= 5, "an extruded quad is a box-like solid");
    }

    #[test]
    fn editing_a_parameter_redrives_the_build() {
        let mut tl = FeatureTimeline::new();
        tl.push(Feature::Box { dx: "size".into(), dy: "size".into(), dz: "size".into() });
        let mut p = ParameterTable::new();
        p.set("size", "10");
        assert!(tl.rebuild(&p).is_ok(), "a valid size builds");
        // Edit the parameter to a broken expression — re-resolved on rebuild,
        // so the edit surfaces as an error (proving the value is re-read).
        p.set("size", "10 +");
        assert!(
            matches!(tl.rebuild(&p), Err(TimelineError::Param(_))),
            "the edited parameter is re-resolved each rebuild"
        );
    }

    #[test]
    fn undefined_parameter_errors() {
        let mut tl = FeatureTimeline::new();
        tl.push(Feature::Cylinder { radius: "missing".into(), height: "1".into() });
        assert!(matches!(
            tl.rebuild(&ParameterTable::new()),
            Err(TimelineError::Param(_))
        ));
    }
}
