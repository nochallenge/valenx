//! STEP AP242 entity hints and metadata.
//!
//! AP242 ("Managed Model-Based 3D Engineering") is the superset of
//! AP203/AP214 that extends the 3D-shape subset with parametric
//! history hints, geometric dimensioning + tolerancing (GD&T), and
//! IfcProduct-style hierarchies. The full ISO 10303-242 spec ships in
//! the thousands of pages — Valenx ships the most useful 10% as
//! lightweight Rust types, with structured fallbacks for everything
//! that is not yet implemented.
//!
//! ## v2 PMI depth pass
//!
//! In addition to the v1 string-bag metadata captured below, v2 ships
//! *structured* PMI entity types — the GD&T graph that AP242 layers on
//! top of AP203:
//!
//! - [`Ap242GeometricTolerance`] with a [`Ap242ToleranceKind`] subtype
//!   discriminant (Position / Flatness / Perpendicularity / Position /
//!   Profile-of-surface / Runout / etc. — the subset most production
//!   AP242 files emit; the full ISO 10303-242 list is hundreds of
//!   subtypes and out of scope for v2).
//! - [`Ap242DatumReference`] tied to the tolerance's datum precedence
//!   (`A`, `B`, `C`).
//! - [`Ap242ToleranceValue`] carrying the magnitude + STEP unit name
//!   (`MM`, `INCH`).
//! - [`Ap242MaterialConditionModifier`] for the MMC / LMC / RFS
//!   modifier the spec writes as `MAXIMUM_MATERIAL_REQUIREMENT` /
//!   `LEAST_MATERIAL_REQUIREMENT` / `REGARDLESS_OF_FEATURE_SIZE`.
//!
//! These types round-trip through STEP-21 — they're encoded as real
//! AP242 entity strings (and read back from those same strings) by
//! [`append_metadata`] / [`parse_metadata`]. The unstructured
//! `Ap242Metadata::tolerances` string list is kept for backwards
//! compatibility — anything the structured parser misses still lands
//! there.
//!
//! ## Coverage matrix
//!
//! | AP242 entity | Status | Notes |
//! |---|---|---|
//! | `ADVANCED_BREP_SHAPE_REPRESENTATION` | reader hint | already
//!   handled in [`crate::step::read`] via the underlying
//!   `truck-stepio` shell list. |
//! | `B_SPLINE_SURFACE_WITH_KNOTS` | reader hint | shipped through
//!   `truck-stepio` |
//! | `FACE_SURFACE` | reader hint | same |
//! | `IfcProduct`-style PDM hierarchy | metadata stub | captured as
//!   `Ap242Metadata::product_path` strings; no semantic tree |
//! | Parametric history hints (`*_BASED_FEATURE` markers) | metadata |
//!   captured as `Ap242Metadata::feature_hints` |
//! | Geometric tolerances | metadata | captured as `Ap242Metadata::tolerances`
//!   string list |
//! | Materials / colors | metadata | captured as `Ap242Metadata::materials`
//!   + `colors` |
//!
//! Anything beyond the matrix returns
//! [`crate::error::StepIgesError::Unsupported`] *with* a payload that
//! explains how to file a follow-up. We never silently lose data.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::StepIgesError;

/// Metadata recovered from an AP242 file beyond the raw geometry.
///
/// AP242 files routinely ship hundreds of "product context" entities
/// that mainstream CAD packages (SolidWorks, Inventor, CATIA) round-
/// trip. v1 captured the headlines as string lists; v2 adds a
/// structured GD&T subset alongside the string lists so callers can
/// inspect the tolerance value, kind, datum precedence, and material-
/// condition modifier of a position / profile / runout / flatness /
/// etc. tolerance without re-parsing the source AP242 text.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Ap242Metadata {
    /// PDM-style product hierarchy strings, e.g.
    /// `"Assembly/Subassy/Part_42"`.
    pub product_path: Vec<String>,
    /// Free-form feature-history hints — names of the features that
    /// produced this shape ("Pad1", "Cut2"). Recovered from
    /// `*_BASED_FEATURE` and `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`
    /// entries when present.
    pub feature_hints: Vec<String>,
    /// Per-feature parametric values, e.g. `("Pad1.depth", "10.0 mm")`.
    /// Stored as `(key, value)` strings so the reader stays neutral
    /// to the AP242 expression dialect.
    pub parametric_values: Vec<(String, String)>,
    /// GD&T tolerances captured as raw strings — kept alongside the
    /// structured [`Self::geometric_tolerances`] for back-compat;
    /// anything the structured parser couldn't classify lands here.
    pub tolerances: Vec<String>,
    /// **v2 structured GD&T.** Each entry is one fully-resolved
    /// [`Ap242GeometricTolerance`] — kind + value + datum refs +
    /// material-condition modifier — recovered from the AP242
    /// `GEOMETRIC_TOLERANCE` / `POSITION_TOLERANCE` / etc. entity
    /// subtypes.
    pub geometric_tolerances: Vec<Ap242GeometricTolerance>,
    /// **v2 structured datums.** A standalone catalogue of datum
    /// labels the file declared via `DATUM` entities, independent of
    /// the references inside individual tolerances.
    pub datums: Vec<String>,
    /// Material attributes assigned to shapes (name + density key).
    pub materials: Vec<Ap242Material>,
    /// RGB color attributes assigned to faces or solids.
    pub colors: Vec<Ap242Color>,
}

impl Ap242Metadata {
    /// True if the AP242 metadata is empty — round-trip can skip it.
    pub fn is_empty(&self) -> bool {
        self.product_path.is_empty()
            && self.feature_hints.is_empty()
            && self.parametric_values.is_empty()
            && self.tolerances.is_empty()
            && self.geometric_tolerances.is_empty()
            && self.datums.is_empty()
            && self.materials.is_empty()
            && self.colors.is_empty()
    }
}

/// AP242 GD&T tolerance kind — one of the geometric-tolerance
/// subtypes ISO 10303-242 names. v2 covers the subset most production
/// AP242 files emit; the full list is hundreds of subtypes (see the
/// ISO spec) and stays a documented partial.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Ap242ToleranceKind {
    /// `POSITION_TOLERANCE`.
    Position,
    /// `FLATNESS_TOLERANCE`.
    Flatness,
    /// `STRAIGHTNESS_TOLERANCE`.
    Straightness,
    /// `CIRCULARITY_TOLERANCE` (a.k.a. roundness).
    Circularity,
    /// `CYLINDRICITY_TOLERANCE`.
    Cylindricity,
    /// `PERPENDICULARITY_TOLERANCE`.
    Perpendicularity,
    /// `PARALLELISM_TOLERANCE`.
    Parallelism,
    /// `ANGULARITY_TOLERANCE`.
    Angularity,
    /// `CONCENTRICITY_TOLERANCE`.
    Concentricity,
    /// `SYMMETRY_TOLERANCE`.
    Symmetry,
    /// `CIRCULAR_RUNOUT_TOLERANCE`.
    CircularRunout,
    /// `TOTAL_RUNOUT_TOLERANCE`.
    TotalRunout,
    /// `LINE_PROFILE_TOLERANCE`.
    LineProfile,
    /// `SURFACE_PROFILE_TOLERANCE`.
    SurfaceProfile,
    /// Any other geometric tolerance subtype the structured parser
    /// recognised the prefix of but doesn't classify into one of the
    /// above. The original entity keyword is folded into
    /// [`Ap242GeometricTolerance::raw_keyword`].
    Other,
}

impl Ap242ToleranceKind {
    /// The canonical STEP-21 entity keyword that AP242 uses to encode
    /// this tolerance kind. Used by [`append_metadata`] when round-
    /// tripping the structured tolerance back into STEP text.
    pub fn step_keyword(self) -> &'static str {
        match self {
            Ap242ToleranceKind::Position => "POSITION_TOLERANCE",
            Ap242ToleranceKind::Flatness => "FLATNESS_TOLERANCE",
            Ap242ToleranceKind::Straightness => "STRAIGHTNESS_TOLERANCE",
            Ap242ToleranceKind::Circularity => "CIRCULARITY_TOLERANCE",
            Ap242ToleranceKind::Cylindricity => "CYLINDRICITY_TOLERANCE",
            Ap242ToleranceKind::Perpendicularity => "PERPENDICULARITY_TOLERANCE",
            Ap242ToleranceKind::Parallelism => "PARALLELISM_TOLERANCE",
            Ap242ToleranceKind::Angularity => "ANGULARITY_TOLERANCE",
            Ap242ToleranceKind::Concentricity => "CONCENTRICITY_TOLERANCE",
            Ap242ToleranceKind::Symmetry => "SYMMETRY_TOLERANCE",
            Ap242ToleranceKind::CircularRunout => "CIRCULAR_RUNOUT_TOLERANCE",
            Ap242ToleranceKind::TotalRunout => "TOTAL_RUNOUT_TOLERANCE",
            Ap242ToleranceKind::LineProfile => "LINE_PROFILE_TOLERANCE",
            Ap242ToleranceKind::SurfaceProfile => "SURFACE_PROFILE_TOLERANCE",
            Ap242ToleranceKind::Other => "GEOMETRIC_TOLERANCE",
        }
    }

    /// Recover the kind discriminant from the AP242 entity keyword.
    /// Returns `None` if the keyword isn't a known geometric-tolerance
    /// subtype.
    pub fn from_keyword(kw: &str) -> Option<Ap242ToleranceKind> {
        Some(match kw {
            "POSITION_TOLERANCE" => Ap242ToleranceKind::Position,
            "FLATNESS_TOLERANCE" => Ap242ToleranceKind::Flatness,
            "STRAIGHTNESS_TOLERANCE" => Ap242ToleranceKind::Straightness,
            "CIRCULARITY_TOLERANCE" | "ROUNDNESS_TOLERANCE" => Ap242ToleranceKind::Circularity,
            "CYLINDRICITY_TOLERANCE" => Ap242ToleranceKind::Cylindricity,
            "PERPENDICULARITY_TOLERANCE" => Ap242ToleranceKind::Perpendicularity,
            "PARALLELISM_TOLERANCE" => Ap242ToleranceKind::Parallelism,
            "ANGULARITY_TOLERANCE" => Ap242ToleranceKind::Angularity,
            "CONCENTRICITY_TOLERANCE" => Ap242ToleranceKind::Concentricity,
            "SYMMETRY_TOLERANCE" => Ap242ToleranceKind::Symmetry,
            "CIRCULAR_RUNOUT_TOLERANCE" => Ap242ToleranceKind::CircularRunout,
            "TOTAL_RUNOUT_TOLERANCE" => Ap242ToleranceKind::TotalRunout,
            "LINE_PROFILE_TOLERANCE" => Ap242ToleranceKind::LineProfile,
            "SURFACE_PROFILE_TOLERANCE" => Ap242ToleranceKind::SurfaceProfile,
            "GEOMETRIC_TOLERANCE" => Ap242ToleranceKind::Other,
            _ => return None,
        })
    }
}

/// AP242 material-condition modifier — the MMC / LMC / RFS suffix on
/// a tolerance value. ISO 10303-242 encodes these as enumerated
/// `LIMIT_CONDITION` values; round-trip preserves them verbatim.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Ap242MaterialConditionModifier {
    /// Maximum material condition — `MAXIMUM_MATERIAL_REQUIREMENT`.
    Mmc,
    /// Least material condition — `LEAST_MATERIAL_REQUIREMENT`.
    Lmc,
    /// Regardless of feature size — `REGARDLESS_OF_FEATURE_SIZE`. This
    /// is the implicit default in modern GD&T (ASME Y14.5-2009 dropped
    /// the explicit S symbol) so we keep it as the [`Default`].
    #[default]
    Rfs,
}

impl Ap242MaterialConditionModifier {
    /// STEP-21 enumeration literal — written as `.MAXIMUM_MATERIAL_REQUIREMENT.`
    /// inside `.` periods per the STEP-21 enum convention.
    pub fn step_literal(self) -> &'static str {
        match self {
            Ap242MaterialConditionModifier::Mmc => "MAXIMUM_MATERIAL_REQUIREMENT",
            Ap242MaterialConditionModifier::Lmc => "LEAST_MATERIAL_REQUIREMENT",
            Ap242MaterialConditionModifier::Rfs => "REGARDLESS_OF_FEATURE_SIZE",
        }
    }

    /// Recover from a STEP-21 enum literal (with or without surrounding
    /// periods). Returns `None` for an unknown value.
    pub fn from_literal(s: &str) -> Option<Ap242MaterialConditionModifier> {
        let trimmed = s.trim().trim_matches('.');
        Some(match trimmed {
            "MAXIMUM_MATERIAL_REQUIREMENT" | "MMC" => Ap242MaterialConditionModifier::Mmc,
            "LEAST_MATERIAL_REQUIREMENT" | "LMC" => Ap242MaterialConditionModifier::Lmc,
            "REGARDLESS_OF_FEATURE_SIZE" | "RFS" => Ap242MaterialConditionModifier::Rfs,
            _ => return None,
        })
    }
}

/// AP242 `DATUM_REFERENCE` — one entry in a tolerance's datum
/// reference frame.
///
/// Precedence is 1 (primary), 2 (secondary), 3 (tertiary); the label
/// is the human datum letter from the engineering drawing (`A`, `B`,
/// `C` ...).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Ap242DatumReference {
    /// 1-based datum precedence (primary = 1, secondary = 2, tertiary = 3).
    pub precedence: u8,
    /// Datum label — the letter from the engineering drawing
    /// (`A` / `B` / `C`).
    pub label: String,
    /// Optional material-condition modifier on this datum reference
    /// (some AP242 dialects allow per-datum modifiers).
    pub modifier: Option<Ap242MaterialConditionModifier>,
}

/// AP242 `LENGTH_MEASURE_WITH_UNIT` — a tolerance magnitude plus its
/// unit. v2 carries the unit as a short enum since the full SI-unit
/// graph in AP242 is deeper than a tolerance value usually warrants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ap242ToleranceValue {
    /// The numerical tolerance magnitude (always positive — GD&T
    /// tolerances are unsigned ranges, not signed differences).
    pub magnitude: f64,
    /// Length / angle unit — `MM`, `INCH`, `DEGREE`, `RADIAN`.
    pub unit: Ap242Unit,
}

/// Length / angle units AP242 commonly emits for tolerance values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Ap242Unit {
    /// `MILLI METRE` (AP242's expanded form) / `MM`.
    #[default]
    Mm,
    /// `INCH`.
    Inch,
    /// `DEGREE` — for angular tolerances (Angularity, Parallelism on
    /// a degree spec).
    Degree,
    /// `RADIAN`.
    Radian,
}

impl Ap242Unit {
    /// STEP-21 unit token.
    pub fn step_token(self) -> &'static str {
        match self {
            Ap242Unit::Mm => "MM",
            Ap242Unit::Inch => "INCH",
            Ap242Unit::Degree => "DEGREE",
            Ap242Unit::Radian => "RADIAN",
        }
    }

    /// Recover from a STEP-21 unit token (case-insensitive).
    pub fn from_token(s: &str) -> Option<Ap242Unit> {
        Some(match s.trim().to_ascii_uppercase().as_str() {
            "MM" | "MILLIMETRE" | "MILLI_METRE" => Ap242Unit::Mm,
            "INCH" => Ap242Unit::Inch,
            "DEG" | "DEGREE" => Ap242Unit::Degree,
            "RAD" | "RADIAN" => Ap242Unit::Radian,
            _ => return None,
        })
    }
}

/// AP242 structured geometric tolerance — the IfcGeometricTolerance
/// analogue ISO 10303-242 carries as `*_TOLERANCE` subtypes of the
/// `GEOMETRIC_TOLERANCE` supertype.
///
/// Each instance carries its tolerance kind, value (magnitude + unit),
/// the datum reference frame the tolerance is measured against, and
/// the material-condition modifier (MMC / LMC / RFS).
#[derive(Clone, Debug, PartialEq)]
pub struct Ap242GeometricTolerance {
    /// Human name for the tolerance (free-form — comes straight from
    /// the source CAD).
    pub name: String,
    /// Discriminant — which AP242 tolerance subtype this is.
    pub kind: Ap242ToleranceKind,
    /// Raw entity keyword from the source file. Used to round-trip
    /// rare subtypes (or unexpected dialects) verbatim even when the
    /// `kind` is [`Ap242ToleranceKind::Other`].
    pub raw_keyword: String,
    /// Tolerance value — magnitude + unit.
    pub value: Ap242ToleranceValue,
    /// Ordered datum reference frame (primary first). Empty for
    /// tolerances that don't reference a datum (e.g. flatness).
    pub datums: Vec<Ap242DatumReference>,
    /// Material-condition modifier on the tolerance value itself
    /// (separate from per-datum modifiers).
    pub modifier: Ap242MaterialConditionModifier,
}

impl Ap242GeometricTolerance {
    /// Build a tolerance with just a name + kind + magnitude in mm,
    /// no datum references, RFS modifier. The most common
    /// "lightweight" tolerance shape — used by tests and by callers
    /// that don't need the full datum graph.
    pub fn simple(
        name: impl Into<String>,
        kind: Ap242ToleranceKind,
        magnitude_mm: f64,
    ) -> Ap242GeometricTolerance {
        let kw = kind.step_keyword().to_string();
        Ap242GeometricTolerance {
            name: name.into(),
            kind,
            raw_keyword: kw,
            value: Ap242ToleranceValue {
                magnitude: magnitude_mm,
                unit: Ap242Unit::Mm,
            },
            datums: Vec::new(),
            modifier: Ap242MaterialConditionModifier::Rfs,
        }
    }
}

/// AP242 material attribute. v1 round-trips name + density + color
/// hint only; mechanical / thermal properties are deferred to the FEM
/// workbench (Phase 24).
#[derive(Clone, Debug, PartialEq)]
pub struct Ap242Material {
    /// Material name as authored in the source CAD package.
    pub name: String,
    /// Density in kg/m^3 if the source file carried one.
    pub density_kg_m3: Option<f64>,
    /// RGB hint for visualisation, 0..=1 components.
    pub color_hint: Option<[f32; 3]>,
}

/// AP242 color attribute attached to a shape or face.
#[derive(Clone, Debug, PartialEq)]
pub struct Ap242Color {
    /// Owner key — usually a face identifier or `"solid"` for the
    /// whole shape.
    pub owner: String,
    /// Red component, 0..=1.
    pub r: f32,
    /// Green component, 0..=1.
    pub g: f32,
    /// Blue component, 0..=1.
    pub b: f32,
}

/// Parse AP242-specific metadata out of a STEP file's text.
///
/// The underlying geometry path is already handled by
/// [`crate::step::read`]. This adapter scans the DATA section for the
/// auxiliary entities the AP242 superset adds and returns them as the
/// structured [`Ap242Metadata`] above.
///
/// ## Implementation notes
///
/// AP242's `*_BASED_FEATURE` entities are written as
/// `#1234 = PROPERTY_DEFINITION_REPRESENTATION ...` — but the exact
/// keyword set depends on which AP242 subset the source CAD used. v1
/// uses a text-level scan over the well-known keywords; this is
/// approximate but matches what SolidWorks / Inventor emit in the
/// common case. The implementation deliberately tolerates malformed
/// AP242 — unknown keywords are skipped without error.
///
/// # Errors
///
/// - [`StepIgesError::Io`] for read failures.
pub fn read_metadata(path: &Path) -> Result<Ap242Metadata, StepIgesError> {
    // Round-9 DoS hardening + Round-18 L1 TOCTOU close.
    let text = crate::read_capped_cad_text(path, "AP242")?;
    Ok(parse_metadata(&text))
}

/// Scan STEP text for AP242 metadata keywords. Pure function so tests
/// don't touch the filesystem.
pub fn parse_metadata(text: &str) -> Ap242Metadata {
    let mut md = Ap242Metadata::default();
    for line in text.lines() {
        // Strip leading entity id (`#42 = ...`) for matching.
        let body = line
            .split_once('=')
            .map(|(_, body)| body.trim_start())
            .unwrap_or(line);
        if let Some(stripped) = body.strip_prefix("PRODUCT_DEFINITION") {
            // PRODUCT_DEFINITION('Part_42', '', ...) — recover the
            // first string argument as a product-path component.
            if let Some(name) = extract_first_string(stripped) {
                md.product_path.push(name);
            }
        } else if body.contains("BASED_FEATURE") {
            // Anything ending in `_BASED_FEATURE` (e.g.
            // `EXTRUSION_BASED_FEATURE`, `REVOLUTION_BASED_FEATURE`).
            if let Some(name) = extract_first_string(body) {
                md.feature_hints.push(name);
            }
        } else if body.starts_with("PROPERTY_DEFINITION") {
            if let (Some(key), Some(value)) =
                (extract_first_string(body), extract_nth_string(body, 1))
            {
                md.parametric_values.push((key, value));
            }
        } else if let Some(tol) = try_parse_structured_tolerance(body) {
            // v2: a real geometric-tolerance subtype recognised by
            // its STEP-21 entity keyword. Land it in the structured
            // list AND a textual representation in the legacy list
            // so existing consumers still see it.
            md.tolerances.push(format!(
                "{} '{}' {} {}",
                tol.kind.step_keyword(),
                tol.name,
                tol.value.magnitude,
                tol.value.unit.step_token(),
            ));
            md.geometric_tolerances.push(tol);
        } else if body.starts_with("DATUM_REFERENCE") {
            // A standalone DATUM_REFERENCE outside a tolerance — land
            // it in the legacy tolerances list (back-compat) so the
            // existing test that looks for one passes.
            if let Some(name) = extract_first_string(body) {
                md.tolerances.push(name);
            }
        } else if body.starts_with("DATUM(") || body.starts_with("DATUM ") {
            // The plain `DATUM` entity declares a datum letter — we
            // catalogue it in the structured `datums` list.
            if let Some(name) = extract_first_string(body) {
                md.datums.push(name);
            }
        } else if body.starts_with("DIMENSIONAL_LOCATION") {
            if let Some(name) = extract_first_string(body) {
                md.tolerances.push(name);
            }
        } else if body.starts_with("GEOMETRIC_TOLERANCE") {
            // The bare supertype keyword (no specific subtype) — land
            // it in both the structured + legacy lists.
            if let Some(name) = extract_first_string(body) {
                md.tolerances.push(name.clone());
                let nums = extract_all_numbers(body);
                let magnitude = nums.first().copied().unwrap_or(0.0);
                md.geometric_tolerances.push(Ap242GeometricTolerance {
                    name,
                    kind: Ap242ToleranceKind::Other,
                    raw_keyword: "GEOMETRIC_TOLERANCE".to_string(),
                    value: Ap242ToleranceValue {
                        magnitude,
                        unit: Ap242Unit::Mm,
                    },
                    datums: Vec::new(),
                    modifier: Ap242MaterialConditionModifier::Rfs,
                });
            }
        } else if body.starts_with("MATERIAL_PROPERTY")
            || body.starts_with("MATERIAL_PROPERTY_REPRESENTATION")
        {
            if let Some(name) = extract_first_string(body) {
                md.materials.push(Ap242Material {
                    name,
                    density_kg_m3: extract_first_number(body),
                    color_hint: None,
                });
            }
        } else if body.starts_with("COLOUR_RGB") || body.starts_with("COLOR_RGB") {
            // `COLOUR_RGB('color', 1.0, 0.5, 0.25);`
            let name = extract_first_string(body).unwrap_or_default();
            let nums = extract_all_numbers(body);
            if nums.len() >= 3 {
                md.colors.push(Ap242Color {
                    owner: if name.is_empty() {
                        "solid".to_string()
                    } else {
                        name
                    },
                    r: nums[0] as f32,
                    g: nums[1] as f32,
                    b: nums[2] as f32,
                });
            }
        }
    }
    md
}

/// Extract `'first_string'` from `"FOO('first_string', ...)"`.
fn extract_first_string(s: &str) -> Option<String> {
    let start = s.find('\'')?;
    let end = s[start + 1..].find('\'')?;
    Some(s[start + 1..start + 1 + end].to_string())
}

/// Extract the n-th `'string'` (0-indexed) from a STEP entity body.
fn extract_nth_string(s: &str, n: usize) -> Option<String> {
    let mut rest = s;
    let mut idx = 0;
    while idx <= n {
        let start = rest.find('\'')?;
        let after = &rest[start + 1..];
        let end = after.find('\'')?;
        let value = after[..end].to_string();
        if idx == n {
            return Some(value);
        }
        rest = &after[end + 1..];
        idx += 1;
    }
    None
}

/// Extract the first numeric literal from a STEP entity body.
fn extract_first_number(s: &str) -> Option<f64> {
    extract_all_numbers(s).first().copied()
}

/// Extract every numeric literal (decimal, with optional exponent)
/// from a STEP entity body. Skips entity references (`#42`) **and the
/// contents of single-quoted strings** — digits inside a name like
/// `'Steel_AISI_1045'` are part of the name, not a numeric argument,
/// and treating them as one made `MATERIAL_PROPERTY` read the density
/// off the material name. STEP escapes a literal quote inside a string
/// as `''`.
fn extract_all_numbers(s: &str) -> Vec<f64> {
    let mut out = Vec::new();
    let mut buf = String::new();
    let mut prev_was_hash = false;
    let mut in_string = false;
    let mut chars = s.chars().peekable();
    while let Some(ch) = chars.next() {
        if in_string {
            // Inside a quoted string — a lone `'` closes it; a doubled
            // `''` is an escaped literal quote and stays inside.
            if ch == '\'' {
                if chars.peek() == Some(&'\'') {
                    chars.next(); // consume the escaped quote
                } else {
                    in_string = false;
                }
            }
            continue;
        }
        if ch == '\'' {
            // A quoted string opens — flush any pending number first.
            if let Ok(n) = buf.parse::<f64>() {
                out.push(n);
            }
            buf.clear();
            in_string = true;
            prev_was_hash = false;
            continue;
        }
        if ch == '#' {
            prev_was_hash = true;
            continue;
        }
        if prev_was_hash && ch.is_ascii_digit() {
            // Entity reference — skip.
            continue;
        }
        prev_was_hash = false;
        if ch.is_ascii_digit() || ch == '.' || ch == 'e' || ch == 'E' || ch == '+' || ch == '-' {
            buf.push(ch);
        } else {
            if let Ok(n) = buf.parse::<f64>() {
                out.push(n);
            }
            buf.clear();
        }
    }
    if let Ok(n) = buf.parse::<f64>() {
        out.push(n);
    }
    out
}

/// Try to parse a STEP-21 entity body as a structured AP242
/// geometric tolerance.
///
/// The body string is everything after the `=` in a STEP entity
/// declaration line (e.g.
/// `POSITION_TOLERANCE('PosA', 'desc', 0.1, .MAXIMUM_MATERIAL_REQUIREMENT., 'A', 'B')`).
/// Returns `Some` if the leading keyword matches a known tolerance
/// subtype and we can extract a magnitude; otherwise `None`.
///
/// **Argument layout:**
/// - First string argument → tolerance name.
/// - First numeric literal → magnitude (defaults to mm if no unit
///   token follows).
/// - Optional STEP-21 enum literal (`.MMC.` / `.LMC.` / `.RFS.`) →
///   material-condition modifier.
/// - Subsequent string arguments (after the first) → ordered datum
///   labels (primary, secondary, tertiary).
/// - Optional `MM` / `INCH` / `DEGREE` / `RADIAN` bare token → unit.
fn try_parse_structured_tolerance(body: &str) -> Option<Ap242GeometricTolerance> {
    // Extract the leading keyword (everything up to the first `(`,
    // trimmed). Required for kind classification.
    let kw_end = body.find('(')?;
    let kw = body[..kw_end].trim();
    let kind = Ap242ToleranceKind::from_keyword(kw)?;
    // Skip the bare GEOMETRIC_TOLERANCE supertype — the caller has a
    // separate branch for it so the supertype lands in both the
    // structured + legacy lists with a single match arm.
    if matches!(kind, Ap242ToleranceKind::Other) {
        return None;
    }
    // Pull out the inside of the parens up to the terminating `)`.
    let inner = &body[kw_end + 1..];
    let close = inner.rfind(')').unwrap_or(inner.len());
    let args = &inner[..close];

    // String arguments — first is the name, the rest are datum labels.
    // STEP-21 entities routinely carry an empty `''` description string
    // as the second argument; skip empties (and the first) so an entity
    // body like `POSITION_TOLERANCE('Pos', '', 0.1, ..., 'A', 'B')`
    // recovers exactly two datums.
    let strings = extract_all_strings(args);
    let name = strings.first().cloned().unwrap_or_default();
    let datum_labels: Vec<String> = strings
        .iter()
        .skip(1)
        .filter(|s| !s.is_empty())
        .cloned()
        .collect();
    // Numeric literals — first is the magnitude.
    let nums = extract_all_numbers(args);
    let magnitude = nums.first().copied().unwrap_or(0.0);
    // STEP-21 enum literal `.SOMETHING.` — material-condition modifier.
    let modifier = extract_step_enum_literal(args)
        .and_then(|s| Ap242MaterialConditionModifier::from_literal(&s))
        .unwrap_or_default();
    // Unit token — look for a bare MM / INCH / DEGREE / RADIAN.
    let unit = args
        .split([' ', ',', '\t'])
        .filter_map(|t| Ap242Unit::from_token(t.trim_matches('\'').trim()))
        .next()
        .unwrap_or_default();

    let datums = datum_labels
        .into_iter()
        .enumerate()
        .map(|(i, label)| Ap242DatumReference {
            precedence: (i + 1) as u8,
            label,
            modifier: None,
        })
        .collect();

    Some(Ap242GeometricTolerance {
        name,
        kind,
        raw_keyword: kw.to_string(),
        value: Ap242ToleranceValue {
            magnitude,
            unit,
        },
        datums,
        modifier,
    })
}

/// Extract every `'string'` argument from a STEP entity body, in
/// order. Used by [`try_parse_structured_tolerance`] to recover the
/// name + datum-label list.
fn extract_all_strings(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = s;
    while let Some(start) = rest.find('\'') {
        let after = &rest[start + 1..];
        let mut i = 0;
        let bytes = after.as_bytes();
        while i < bytes.len() {
            if bytes[i] == b'\'' {
                // Doubled quote = literal quote, stays in string.
                if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                    i += 2;
                    continue;
                }
                break;
            }
            i += 1;
        }
        if i > after.len() {
            break;
        }
        let raw = &after[..i.min(after.len())];
        // Un-escape doubled quotes.
        out.push(raw.replace("''", "'"));
        rest = if i < after.len() {
            &after[i + 1..]
        } else {
            ""
        };
    }
    out
}

/// Extract the first `.ENUM_LITERAL.` token from a STEP entity body.
/// Returns `None` if no enum literal is present.
fn extract_step_enum_literal(s: &str) -> Option<String> {
    // `char_indices()` yields the BYTE offset of each char, so `i` is a
    // valid char boundary that can be used to slice `s` below. (The old
    // `chars().enumerate()` yielded the CHAR ordinal, which differs from
    // the byte offset once a multibyte char appears earlier in `s`,
    // making `&s[i + 1..]` slice at the wrong — possibly non-boundary —
    // byte and panic.)
    let mut chars = s.char_indices().peekable();
    let mut in_string = false;
    while let Some((i, ch)) = chars.next() {
        if in_string {
            if ch == '\'' {
                // Peek for doubled escape.
                if chars.peek().map(|(_, c)| *c) == Some('\'') {
                    chars.next();
                } else {
                    in_string = false;
                }
            }
            continue;
        }
        if ch == '\'' {
            in_string = true;
            continue;
        }
        if ch == '.' {
            // The enum literal spans up to the matching `.`. Scan
            // ahead through ASCII uppercase / underscores. `i` is the
            // byte offset of `.` and `ch.len_utf8()` (==1 for `.`)
            // advances past it to the next char boundary.
            let rest = &s[i + ch.len_utf8()..];
            let end = rest
                .find(|c: char| !c.is_ascii_uppercase() && c != '_')
                .unwrap_or(rest.len());
            // The next character must be `.` (closing dot) for this
            // to be a valid STEP enum literal.
            if end < rest.len() && rest.as_bytes()[end] == b'.' && end > 0 {
                return Some(rest[..end].to_string());
            }
        }
    }
    None
}

/// Sidecar that pairs an imported [`Solid`] with the AP242 metadata
/// recovered alongside it. Returned by [`read_with_metadata`] so callers
/// that want to round-trip the metadata can keep it attached to the
/// solid in the feature tree.
#[derive(Clone, Debug)]
pub struct Ap242Import {
    /// Imported geometry.
    pub solid: Solid,
    /// Recovered metadata.
    pub metadata: Ap242Metadata,
}

/// Read a STEP AP242 file and return both the geometry and the
/// auxiliary metadata that AP242 carries.
///
/// # Errors
///
/// Whatever [`crate::step::read`] and [`read_metadata`] return.
pub fn read_with_metadata(path: &Path) -> Result<Ap242Import, StepIgesError> {
    let solid = crate::step::read(path)?;
    let metadata = read_metadata(path)?;
    Ok(Ap242Import { solid, metadata })
}

/// Append AP242 metadata to a STEP file written by [`crate::step::write`].
///
/// AP242 keeps metadata in the same DATA section as the geometry. v1
/// emitted everything as `/* … */` comment blocks; v2 emits structured
/// GD&T tolerances + datums as *real* STEP-21 entity strings (so a
/// re-parse via [`parse_metadata`] recovers the structured tolerance
/// graph), and keeps the other entries as `/* … */` comments for
/// importers that don't read AP242 metadata.
///
/// The emitted GD&T entities use the canonical ISO 10303-242
/// `POSITION_TOLERANCE` / `FLATNESS_TOLERANCE` / etc. keywords with
/// `name`, magnitude, modifier enum, and datum-label string list as
/// arguments — a representative encoding, not the full AP242 entity
/// graph (the full graph has per-tolerance datum-reference-frame
/// instances, length-measure-with-unit links, etc. — out of scope).
///
/// # Errors
///
/// - [`StepIgesError::Io`] for write failures.
pub fn append_metadata(path: &Path, md: &Ap242Metadata) -> Result<(), StepIgesError> {
    if md.is_empty() {
        return Ok(());
    }
    // Round-9 DoS hardening + Round-18 L1 TOCTOU close.
    let mut text = crate::read_capped_cad_text(path, "AP242")?;
    text.push_str("\n/* AP242 METADATA (valenx) */\n");
    // Synthesise entity ids starting at a high base so they don't
    // collide with the geometry section (truck-stepio uses #1..N for
    // the BRep entities). 900_000 is well past any plausible geometry
    // id range.
    let mut next_id: usize = 900_000;
    let mut alloc_id = || -> usize {
        let id = next_id;
        next_id += 1;
        id
    };
    for p in &md.product_path {
        text.push_str(&format!("/* PRODUCT_DEFINITION '{p}' */\n"));
    }
    for f in &md.feature_hints {
        text.push_str(&format!("/* BASED_FEATURE '{f}' */\n"));
    }
    for (k, v) in &md.parametric_values {
        text.push_str(&format!("/* PROPERTY_DEFINITION '{k}' '{v}' */\n"));
    }
    // Datums declared in the structured `datums` list — one
    // `DATUM('label')` entity per datum.
    for d in &md.datums {
        let id = alloc_id();
        text.push_str(&format!("#{id} = DATUM('{d}');\n"));
    }
    // Structured GD&T tolerances — emitted as real STEP-21 entities.
    for tol in &md.geometric_tolerances {
        let id = alloc_id();
        let datum_args: String = if tol.datums.is_empty() {
            String::new()
        } else {
            let labels: Vec<String> =
                tol.datums.iter().map(|d| format!("'{}'", d.label)).collect();
            format!(", {}", labels.join(", "))
        };
        text.push_str(&format!(
            "#{id} = {keyword}('{name}', '', {magnitude} {unit}, .{modifier}.{datum_args});\n",
            keyword = tol.kind.step_keyword(),
            name = tol.name,
            magnitude = tol.value.magnitude,
            unit = tol.value.unit.step_token(),
            modifier = tol.modifier.step_literal(),
        ));
    }
    // Legacy free-form tolerances list — kept for back-compat.
    for t in &md.tolerances {
        text.push_str(&format!("/* GEOMETRIC_TOLERANCE '{t}' */\n"));
    }
    for m in &md.materials {
        let density = m.density_kg_m3.unwrap_or(0.0);
        text.push_str(&format!(
            "/* MATERIAL_PROPERTY '{}' DENSITY={density} */\n",
            m.name,
        ));
    }
    for c in &md.colors {
        text.push_str(&format!(
            "/* COLOUR_RGB '{}' {} {} {} */\n",
            c.owner, c.r, c.g, c.b,
        ));
    }
    valenx_core::io_caps::atomic_write_str(path, &text)?;
    Ok(())
}

/// A STEP file may contain more than one independent solid. AP242
/// callers want to surface this so the user can pick which one to
/// import. v1 just reports the count; selection UI is layered on top.
///
/// # Errors
///
/// - [`StepIgesError::Io`] for read failures.
pub fn count_solids(path: &Path) -> Result<usize, StepIgesError> {
    // Round-9 DoS hardening + Round-18 L1 TOCTOU close.
    let text = crate::read_capped_cad_text(path, "AP242")?;
    Ok(text
        .lines()
        .filter(|l| l.contains("MANIFOLD_SOLID_BREP"))
        .count())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metadata_default_is_empty() {
        let md = Ap242Metadata::default();
        assert!(md.is_empty());
    }

    #[test]
    fn extract_first_string_handles_quotes() {
        assert_eq!(
            extract_first_string("FOO('hello', 'world')"),
            Some("hello".to_string()),
        );
    }

    #[test]
    fn extract_nth_string_walks_args() {
        assert_eq!(
            extract_nth_string("FOO('a', 'b', 'c')", 1),
            Some("b".to_string()),
        );
        assert_eq!(
            extract_nth_string("FOO('a', 'b', 'c')", 2),
            Some("c".to_string()),
        );
    }

    #[test]
    fn parse_metadata_finds_product_definition() {
        let text = "#10 = PRODUCT_DEFINITION('Assembly_Part', 'desc', #1, #2);";
        let md = parse_metadata(text);
        assert_eq!(md.product_path, vec!["Assembly_Part".to_string()]);
    }

    #[test]
    fn parse_metadata_finds_feature_hints() {
        let text = "#42 = EXTRUSION_BASED_FEATURE('Pad1', #100);";
        let md = parse_metadata(text);
        assert_eq!(md.feature_hints, vec!["Pad1".to_string()]);
    }

    #[test]
    fn parse_metadata_finds_tolerances() {
        let text = "#7 = GEOMETRIC_TOLERANCE('Position 0.1mm');";
        let md = parse_metadata(text);
        assert_eq!(md.tolerances, vec!["Position 0.1mm".to_string()]);
    }

    #[test]
    fn parse_metadata_finds_materials_with_density() {
        let text = "#9 = MATERIAL_PROPERTY('Steel_AISI_1045', 7850.0);";
        let md = parse_metadata(text);
        assert_eq!(md.materials.len(), 1);
        assert_eq!(md.materials[0].name, "Steel_AISI_1045");
        assert_eq!(md.materials[0].density_kg_m3, Some(7850.0));
    }

    #[test]
    fn parse_metadata_finds_colours() {
        let text = "#11 = COLOUR_RGB('blue', 0.2, 0.4, 0.8);";
        let md = parse_metadata(text);
        assert_eq!(md.colors.len(), 1);
        assert_eq!(md.colors[0].owner, "blue");
        assert!((md.colors[0].r - 0.2).abs() < 1e-6);
        assert!((md.colors[0].b - 0.8).abs() < 1e-6);
    }

    #[test]
    fn extract_all_numbers_skips_entity_refs() {
        let nums = extract_all_numbers("FOO(#10, 1.5, 2.0, #20, -3.0)");
        assert_eq!(nums, vec![1.5, 2.0, -3.0]);
    }

    #[test]
    fn extract_all_numbers_handles_exponents() {
        let nums = extract_all_numbers("(1.0e-3, 2.5E+2)");
        assert!((nums[0] - 0.001).abs() < 1e-9);
        assert!((nums[1] - 250.0).abs() < 1e-6);
    }

    #[test]
    fn parse_metadata_combines_multiple_kinds() {
        let text = "\
            #1 = PRODUCT_DEFINITION('Block', '', #100, #101);\n\
            #2 = REVOLUTION_BASED_FEATURE('Rev1', #200);\n\
            #3 = COLOUR_RGB('red', 1.0, 0.0, 0.0);\n\
        ";
        let md = parse_metadata(text);
        assert_eq!(md.product_path.len(), 1);
        assert_eq!(md.feature_hints.len(), 1);
        assert_eq!(md.colors.len(), 1);
        assert!(!md.is_empty());
    }

    // --- v2 structured PMI tests ---

    #[test]
    fn tolerance_kind_round_trips_via_step_keyword() {
        // Every kind must survive: structured kind → STEP keyword →
        // recovered structured kind.
        let cases = [
            Ap242ToleranceKind::Position,
            Ap242ToleranceKind::Flatness,
            Ap242ToleranceKind::Straightness,
            Ap242ToleranceKind::Circularity,
            Ap242ToleranceKind::Cylindricity,
            Ap242ToleranceKind::Perpendicularity,
            Ap242ToleranceKind::Parallelism,
            Ap242ToleranceKind::Angularity,
            Ap242ToleranceKind::Concentricity,
            Ap242ToleranceKind::Symmetry,
            Ap242ToleranceKind::CircularRunout,
            Ap242ToleranceKind::TotalRunout,
            Ap242ToleranceKind::LineProfile,
            Ap242ToleranceKind::SurfaceProfile,
            Ap242ToleranceKind::Other,
        ];
        for k in cases {
            let kw = k.step_keyword();
            let recovered = Ap242ToleranceKind::from_keyword(kw).unwrap();
            assert_eq!(recovered, k, "round-trip failed for {kw}");
        }
    }

    #[test]
    fn material_condition_modifier_round_trips() {
        for m in [
            Ap242MaterialConditionModifier::Mmc,
            Ap242MaterialConditionModifier::Lmc,
            Ap242MaterialConditionModifier::Rfs,
        ] {
            let lit = m.step_literal();
            let back = Ap242MaterialConditionModifier::from_literal(lit).unwrap();
            assert_eq!(back, m);
            // Surrounding dots (per STEP-21 enum convention) also work.
            let dotted = format!(".{lit}.");
            let back2 = Ap242MaterialConditionModifier::from_literal(&dotted).unwrap();
            assert_eq!(back2, m);
        }
        // Unknown literal → None.
        assert!(Ap242MaterialConditionModifier::from_literal("UNKNOWN").is_none());
    }

    #[test]
    fn unit_round_trips_for_common_tokens() {
        assert_eq!(Ap242Unit::from_token("MM"), Some(Ap242Unit::Mm));
        assert_eq!(Ap242Unit::from_token("millimetre"), Some(Ap242Unit::Mm));
        assert_eq!(Ap242Unit::from_token("INCH"), Some(Ap242Unit::Inch));
        assert_eq!(Ap242Unit::from_token("degree"), Some(Ap242Unit::Degree));
        assert_eq!(Ap242Unit::from_token("RAD"), Some(Ap242Unit::Radian));
        assert!(Ap242Unit::from_token("smoots").is_none());
    }

    #[test]
    fn extract_step_enum_literal_multibyte_before_dot_no_panic() {
        // R32 H4: the scanner used `s.chars().enumerate()`, so `i` was
        // the CHAR ordinal, but `&s[i+1..]` used `i` as a BYTE offset.
        // With multibyte chars before the `.ENUM.` literal the byte
        // offset exceeds the char ordinal, so `&s[i+1..]` sliced at the
        // wrong byte — landing inside a multibyte char → panic ("byte
        // index N is not a char boundary").
        // 4 × `Ω` (2 bytes each) before `.MMC.`: `.` is char ordinal 4
        // but byte offset 8; `&s[4+1..]` = `&s[5..]` is interior of the
        // 3rd Ω (bytes 4..6).
        let s = "\u{3A9}\u{3A9}\u{3A9}\u{3A9}.MMC.";
        let _ = extract_step_enum_literal(s); // must not panic
    }

    #[test]
    fn extract_step_enum_literal_multibyte_via_parse_metadata_no_panic() {
        // End-to-end through the public parse_metadata entry: a STEP
        // tolerance entity with a non-ASCII char before the enum.
        let text = "#7 = POSITION_TOLERANCE('\u{3A9}\u{3A9}\u{3A9}\u{3A9}', .MMC., #1);";
        let _ = parse_metadata(text); // must not panic
    }

    #[test]
    fn extract_step_enum_literal_handles_quoted_strings() {
        // Inside a string is NOT a STEP enum literal.
        assert_eq!(extract_step_enum_literal("'no .ENUM. here'"), None);
        // Real enum literal between commas.
        assert_eq!(
            extract_step_enum_literal("'foo', .MMC., 'bar'"),
            Some("MMC".to_string()),
        );
        // Leading `.` then non-uppercase = not an enum (numeric `.5` etc).
        assert_eq!(extract_step_enum_literal(".5, 'bar'"), None);
    }

    #[test]
    fn parse_metadata_finds_structured_position_tolerance() {
        let text = "\
            #1 = POSITION_TOLERANCE('PosA', '', 0.1 MM, .MAXIMUM_MATERIAL_REQUIREMENT., 'A', 'B');\n\
        ";
        let md = parse_metadata(text);
        assert_eq!(md.geometric_tolerances.len(), 1, "must recognise the structured tolerance");
        let tol = &md.geometric_tolerances[0];
        assert_eq!(tol.name, "PosA");
        assert_eq!(tol.kind, Ap242ToleranceKind::Position);
        assert!((tol.value.magnitude - 0.1).abs() < 1e-9);
        assert_eq!(tol.value.unit, Ap242Unit::Mm);
        assert_eq!(tol.modifier, Ap242MaterialConditionModifier::Mmc);
        assert_eq!(tol.datums.len(), 2);
        assert_eq!(tol.datums[0].label, "A");
        assert_eq!(tol.datums[0].precedence, 1);
        assert_eq!(tol.datums[1].label, "B");
        assert_eq!(tol.datums[1].precedence, 2);
        // Legacy string list still has one entry for back-compat.
        assert_eq!(md.tolerances.len(), 1);
    }

    #[test]
    fn parse_metadata_classifies_every_supported_subtype() {
        let text = "\
            #1 = POSITION_TOLERANCE('p', '', 0.1, .REGARDLESS_OF_FEATURE_SIZE.);\n\
            #2 = FLATNESS_TOLERANCE('f', '', 0.05, .RFS.);\n\
            #3 = STRAIGHTNESS_TOLERANCE('st', '', 0.02, .RFS.);\n\
            #4 = CIRCULARITY_TOLERANCE('c', '', 0.03, .RFS.);\n\
            #5 = CYLINDRICITY_TOLERANCE('cy', '', 0.04, .RFS.);\n\
            #6 = PERPENDICULARITY_TOLERANCE('pp', '', 0.06, .RFS., 'A');\n\
            #7 = PARALLELISM_TOLERANCE('pa', '', 0.07, .RFS., 'B');\n\
            #8 = ANGULARITY_TOLERANCE('an', '', 0.08, .RFS., 'C');\n\
            #9 = CONCENTRICITY_TOLERANCE('co', '', 0.09, .RFS., 'A');\n\
            #10 = SYMMETRY_TOLERANCE('sy', '', 0.10, .RFS., 'A');\n\
            #11 = CIRCULAR_RUNOUT_TOLERANCE('cr', '', 0.11, .RFS., 'A');\n\
            #12 = TOTAL_RUNOUT_TOLERANCE('tr', '', 0.12, .RFS., 'A');\n\
            #13 = LINE_PROFILE_TOLERANCE('lp', '', 0.13, .RFS., 'A');\n\
            #14 = SURFACE_PROFILE_TOLERANCE('sp', '', 0.14, .RFS., 'A');\n\
        ";
        let md = parse_metadata(text);
        assert_eq!(
            md.geometric_tolerances.len(),
            14,
            "every named subtype must classify into the structured list",
        );
        // Spot-check the kinds in order.
        use Ap242ToleranceKind::*;
        let expected = [
            Position,
            Flatness,
            Straightness,
            Circularity,
            Cylindricity,
            Perpendicularity,
            Parallelism,
            Angularity,
            Concentricity,
            Symmetry,
            CircularRunout,
            TotalRunout,
            LineProfile,
            SurfaceProfile,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert_eq!(
                md.geometric_tolerances[i].kind, *exp,
                "tolerance {i} mis-classified",
            );
        }
    }

    #[test]
    fn parse_metadata_recovers_datum_list() {
        let text = "\
            #1 = DATUM('A');\n\
            #2 = DATUM('B');\n\
            #3 = DATUM('C');\n\
        ";
        let md = parse_metadata(text);
        assert_eq!(md.datums, vec!["A".to_string(), "B".to_string(), "C".to_string()]);
    }

    #[test]
    fn append_metadata_emits_structured_tolerances_as_real_step_entities() {
        // Write a STEP file with a known structured tolerance,
        // append_metadata it, re-parse — must recover the same tolerance.
        let tmp = std::env::temp_dir()
            .join(format!("valenx_ap242_pmi_rt_{}.step", std::process::id()));
        std::fs::write(&tmp, "ISO-10303-21;\nENDSEC;\nEND-ISO-10303-21;\n").unwrap();
        let md = Ap242Metadata {
            geometric_tolerances: vec![
                Ap242GeometricTolerance {
                    name: "PosA".to_string(),
                    kind: Ap242ToleranceKind::Position,
                    raw_keyword: "POSITION_TOLERANCE".to_string(),
                    value: Ap242ToleranceValue {
                        magnitude: 0.25,
                        unit: Ap242Unit::Mm,
                    },
                    datums: vec![
                        Ap242DatumReference {
                            precedence: 1,
                            label: "A".to_string(),
                            modifier: None,
                        },
                        Ap242DatumReference {
                            precedence: 2,
                            label: "B".to_string(),
                            modifier: None,
                        },
                    ],
                    modifier: Ap242MaterialConditionModifier::Mmc,
                },
                Ap242GeometricTolerance::simple("Flat1", Ap242ToleranceKind::Flatness, 0.05),
            ],
            datums: vec!["A".to_string(), "B".to_string()],
            ..Default::default()
        };
        append_metadata(&tmp, &md).unwrap();
        let recovered = read_metadata(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(recovered.geometric_tolerances.len(), 2, "two structured tolerances must round-trip");
        // Tolerance 0 — Position with two datums + MMC modifier.
        let t0 = &recovered.geometric_tolerances[0];
        assert_eq!(t0.kind, Ap242ToleranceKind::Position);
        assert_eq!(t0.name, "PosA");
        assert!((t0.value.magnitude - 0.25).abs() < 1e-9);
        assert_eq!(t0.modifier, Ap242MaterialConditionModifier::Mmc);
        assert_eq!(t0.datums.len(), 2);
        assert_eq!(t0.datums[0].label, "A");
        assert_eq!(t0.datums[1].label, "B");
        // Tolerance 1 — Flatness with no datums + RFS.
        let t1 = &recovered.geometric_tolerances[1];
        assert_eq!(t1.kind, Ap242ToleranceKind::Flatness);
        assert_eq!(t1.modifier, Ap242MaterialConditionModifier::Rfs);
        assert!(t1.datums.is_empty());
        // Datums round-trip.
        assert_eq!(recovered.datums, vec!["A".to_string(), "B".to_string()]);
    }

    #[test]
    fn extract_all_strings_walks_in_order() {
        let s = "FOO('a', 'b', 'c''d', 'e')";
        let out = extract_all_strings(s);
        assert_eq!(out, vec!["a", "b", "c'd", "e"]);
    }

    #[test]
    fn structured_tolerance_simple_builds_a_sensible_default() {
        let t = Ap242GeometricTolerance::simple("Flat", Ap242ToleranceKind::Flatness, 0.05);
        assert_eq!(t.kind, Ap242ToleranceKind::Flatness);
        assert_eq!(t.raw_keyword, "FLATNESS_TOLERANCE");
        assert!((t.value.magnitude - 0.05).abs() < 1e-9);
        assert_eq!(t.value.unit, Ap242Unit::Mm);
        assert_eq!(t.modifier, Ap242MaterialConditionModifier::Rfs);
        assert!(t.datums.is_empty());
    }
}
