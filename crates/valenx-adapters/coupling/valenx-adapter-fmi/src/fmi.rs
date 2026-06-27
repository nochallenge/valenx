//! FMI **co-simulation** import: parse an FMU's `modelDescription.xml`.
//!
//! Scope (honest): this is a *model-description importer for co-simulation
//! FMUs* — it extracts the model name and the scalar interface variables
//! (name, value reference, causality) so the [`crate::cosim`] master knows
//! the FMU's input/output ports. It is **not** a model-exchange importer,
//! and the DEFAULT crate path does not require a real binary FMU at all:
//! the native in-process [`crate::cosim::Subsystem`] is what every test
//! drives. Loading and stepping an actual compiled FMU binary
//! (`.so` / `.dll` / `.dylib`) lives behind the off-by-default `binary-fmu`
//! cargo feature (the `binary` module, present only when that feature is
//! enabled).
//!
//! ## Why hand-parse?
//!
//! Per `AGENTS.md` (native-first: prefer extending in-house Rust over
//! adding a dependency), the parser is a small, self-contained scanner
//! rather than a new XML crate in the workspace tree. It understands the
//! two FMI generations:
//!
//! * **FMI 2.0** — `<ScalarVariable name="…" valueReference="…"
//!   causality="…">` children of `<ModelVariables>`.
//! * **FMI 3.0** — typed variable elements (`<Float64>`, `<Int32>`,
//!   `<Boolean>`, …) with the same `name` / `valueReference` / `causality`
//!   attributes, also under `<ModelVariables>`.
//!
//! Anything it cannot make sense of is a fail-loud
//! [`FmiError::ModelDescriptionParse`], never a silent empty list.

use crate::error::{FmiError, Result};

/// Causality of an FMI scalar variable (the subset relevant to wiring a
/// co-simulation coupling). Unknown / unsupported causality strings are a
/// parse error rather than a silent default.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Causality {
    /// Free input the importer drives.
    Input,
    /// Computed output the importer reads.
    Output,
    /// Tunable / fixed parameter.
    Parameter,
    /// Calculated parameter (FMI: `calculatedParameter`).
    CalculatedParameter,
    /// Internal variable exposed for inspection.
    Local,
    /// Independent variable (typically simulation time).
    Independent,
}

impl Causality {
    fn parse(s: &str) -> Result<Self> {
        match s {
            "input" => Ok(Causality::Input),
            "output" => Ok(Causality::Output),
            "parameter" => Ok(Causality::Parameter),
            "calculatedParameter" => Ok(Causality::CalculatedParameter),
            "local" => Ok(Causality::Local),
            "independent" => Ok(Causality::Independent),
            other => Err(FmiError::ModelDescriptionParse(format!(
                "unknown causality {other:?}"
            ))),
        }
    }
}

/// One scalar interface variable of an FMU.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScalarVariable {
    /// Variable name as declared in the FMU.
    pub name: String,
    /// FMI value reference (unique per type within the FMU).
    pub value_reference: u32,
    /// Variable causality.
    pub causality: Causality,
}

/// The parsed interface of an FMU's `modelDescription.xml`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModelDescription {
    /// `modelName` attribute of `<fmiModelDescription>`.
    pub model_name: String,
    /// `fmiVersion` attribute (e.g. `"2.0"` or `"3.0"`), if present.
    pub fmi_version: Option<String>,
    /// Scalar variables, in document order.
    pub variables: Vec<ScalarVariable>,
}

impl ModelDescription {
    /// Parse a `modelDescription.xml` document from a string.
    ///
    /// Fail-loud: malformed XML, a missing `modelName`, an unparsable
    /// `valueReference`, or an unknown `causality` all return
    /// [`FmiError::ModelDescriptionParse`] — never a partial-but-plausible
    /// result.
    pub fn parse(xml: &str) -> Result<Self> {
        let tags = tokenize(xml)?;

        // Find the root <fmiModelDescription ...> and pull modelName +
        // fmiVersion off its attributes.
        let root = tags
            .iter()
            .find(|t| t.name == "fmiModelDescription" && t.kind != TagKind::Close)
            .ok_or_else(|| {
                FmiError::ModelDescriptionParse(
                    "missing <fmiModelDescription> root element".to_string(),
                )
            })?;
        let model_name = root.attr("modelName").ok_or_else(|| {
            FmiError::ModelDescriptionParse(
                "<fmiModelDescription> has no modelName attribute".to_string(),
            )
        })?;
        let fmi_version = root.attr("fmiVersion");

        // Scope variable parsing to the <ModelVariables> … </ModelVariables>
        // span so we don't mistake, e.g., a <ScalarVariable> inside some
        // other (illegal) context. Absence of <ModelVariables> means an
        // FMU with no interface variables — that is unusual enough that we
        // treat it as a parse error to stay fail-loud.
        let mv_open = tags
            .iter()
            .position(|t| t.name == "ModelVariables" && t.kind != TagKind::Close);
        let Some(mv_open) = mv_open else {
            return Err(FmiError::ModelDescriptionParse(
                "missing <ModelVariables> section".to_string(),
            ));
        };
        // If the <ModelVariables> tag is self-closing, the section is empty.
        let mv_close = if tags[mv_open].kind == TagKind::SelfClosing {
            mv_open
        } else {
            tags[mv_open + 1..]
                .iter()
                .position(|t| t.name == "ModelVariables" && t.kind == TagKind::Close)
                .map(|p| p + mv_open + 1)
                .ok_or_else(|| {
                    FmiError::ModelDescriptionParse(
                        "unterminated <ModelVariables> section".to_string(),
                    )
                })?
        };

        let mut variables = Vec::new();
        for tag in &tags[mv_open + 1..mv_close] {
            if !is_variable_element(&tag.name, tag.kind) {
                continue;
            }
            // A variable element carries the interface attributes whether
            // it is self-closing (FMI 3.0 typed scalar, FMI 2.0 with no
            // nested <Real>) or an open tag (FMI 2.0 with a nested type
            // element). Either way the attributes live on the OPEN tag.
            if tag.kind == TagKind::Close {
                continue;
            }
            let name = tag.attr("name").ok_or_else(|| {
                FmiError::ModelDescriptionParse(format!(
                    "<{}> variable has no name attribute",
                    tag.name
                ))
            })?;
            let vr_str = tag.attr("valueReference").ok_or_else(|| {
                FmiError::ModelDescriptionParse(format!(
                    "variable {name:?} has no valueReference attribute"
                ))
            })?;
            let value_reference: u32 = vr_str.trim().parse().map_err(|_| {
                FmiError::ModelDescriptionParse(format!(
                    "variable {name:?} has non-integer valueReference {vr_str:?}"
                ))
            })?;
            // FMI defaults causality to "local" when the attribute is
            // omitted; honour that default explicitly rather than dropping
            // the variable.
            let causality = match tag.attr("causality") {
                Some(c) => Causality::parse(&c)?,
                None => Causality::Local,
            };
            variables.push(ScalarVariable {
                name,
                value_reference,
                causality,
            });
        }

        Ok(ModelDescription {
            model_name,
            fmi_version,
            variables,
        })
    }

    /// Names of the input variables, in document order.
    pub fn inputs(&self) -> Vec<&str> {
        self.variables
            .iter()
            .filter(|v| v.causality == Causality::Input)
            .map(|v| v.name.as_str())
            .collect()
    }

    /// Names of the output variables, in document order.
    pub fn outputs(&self) -> Vec<&str> {
        self.variables
            .iter()
            .filter(|v| v.causality == Causality::Output)
            .map(|v| v.name.as_str())
            .collect()
    }
}

/// Is `name` a `<ModelVariables>` child that declares a scalar variable?
///
/// FMI 2.0 uses `ScalarVariable`; FMI 3.0 uses one element per scalar type.
fn is_variable_element(name: &str, kind: TagKind) -> bool {
    if kind == TagKind::Close {
        return false;
    }
    matches!(
        name,
        // FMI 2.0
        "ScalarVariable"
        // FMI 3.0 typed scalars
        | "Float32" | "Float64"
        | "Int8" | "UInt8" | "Int16" | "UInt16"
        | "Int32" | "UInt32" | "Int64" | "UInt64"
        | "Boolean" | "String" | "Enumeration"
    )
}

// ---------------------------------------------------------------------------
// Minimal XML tag scanner.
//
// This is intentionally small: it splits the document into element tags and
// reads their attributes, which is all the FMI model description needs. It
// is NOT a general XML parser (it does not build a DOM, resolve entities
// beyond the common five, or validate nesting beyond what we use), but it
// is fail-loud on the malformations that matter: an unterminated tag, an
// unterminated quoted attribute, or junk where a tag should be.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TagKind {
    /// `<name …>`
    Open,
    /// `</name>`
    Close,
    /// `<name … />`
    SelfClosing,
}

#[derive(Debug)]
struct Tag {
    name: String,
    kind: TagKind,
    attrs: Vec<(String, String)>,
}

impl Tag {
    fn attr(&self, key: &str) -> Option<String> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.clone())
    }
}

/// Split `xml` into a flat list of element tags, skipping the prolog
/// (`<?xml …?>`), comments (`<!-- … -->`), DOCTYPE/CDATA-style `<!…>`
/// markup, and text content. Fail-loud on a tag that never closes or a
/// quoted attribute value that is never terminated.
fn tokenize(xml: &str) -> Result<Vec<Tag>> {
    let bytes = xml.as_bytes();
    let mut i = 0usize;
    let n = bytes.len();
    let mut tags = Vec::new();

    while i < n {
        // Advance to the next '<'.
        match bytes[i] {
            b'<' => {}
            _ => {
                i += 1;
                continue;
            }
        }

        // Classify the markup that opens at `i`.
        if xml[i..].starts_with("<?") {
            // Processing instruction / XML declaration: skip to "?>".
            let end = xml[i..].find("?>").ok_or_else(|| {
                FmiError::ModelDescriptionParse("unterminated <? … ?> declaration".to_string())
            })?;
            i += end + 2;
            continue;
        }
        if xml[i..].starts_with("<!--") {
            let end = xml[i..].find("-->").ok_or_else(|| {
                FmiError::ModelDescriptionParse("unterminated <!-- … --> comment".to_string())
            })?;
            i += end + 3;
            continue;
        }
        if xml[i..].starts_with("<!") {
            // DOCTYPE / CDATA / other declaration: skip to the next '>'.
            let end = xml[i..].find('>').ok_or_else(|| {
                FmiError::ModelDescriptionParse("unterminated <! … > declaration".to_string())
            })?;
            i += end + 1;
            continue;
        }

        // A real element tag. Find its terminating '>', respecting quotes
        // so a '>' inside an attribute value doesn't end the tag early.
        let tag_start = i + 1;
        let mut j = tag_start;
        let mut quote: Option<u8> = None;
        let mut closed = false;
        while j < n {
            let c = bytes[j];
            match quote {
                Some(q) => {
                    if c == q {
                        quote = None;
                    }
                }
                None => match c {
                    b'"' | b'\'' => quote = Some(c),
                    b'>' => {
                        closed = true;
                        break;
                    }
                    _ => {}
                },
            }
            j += 1;
        }
        if quote.is_some() {
            return Err(FmiError::ModelDescriptionParse(
                "unterminated quoted attribute value".to_string(),
            ));
        }
        if !closed {
            return Err(FmiError::ModelDescriptionParse(
                "unterminated element tag (missing '>')".to_string(),
            ));
        }

        let inner = &xml[tag_start..j]; // between '<' and '>'
        i = j + 1;
        let tag = parse_tag(inner)?;
        tags.push(tag);
    }

    Ok(tags)
}

/// Parse the inside of a tag (everything between `<` and `>`).
fn parse_tag(inner: &str) -> Result<Tag> {
    let inner = inner.trim();
    if inner.is_empty() {
        return Err(FmiError::ModelDescriptionParse("empty <> tag".to_string()));
    }

    let (kind, body) = if let Some(rest) = inner.strip_prefix('/') {
        (TagKind::Close, rest.trim())
    } else if let Some(rest) = inner.strip_suffix('/') {
        (TagKind::SelfClosing, rest.trim())
    } else {
        (TagKind::Open, inner)
    };

    // Element name = everything up to the first whitespace.
    let name_split = body
        .char_indices()
        .find(|(_, c)| c.is_whitespace())
        .map(|(k, _)| k)
        .unwrap_or(body.len());
    let name = body[..name_split].to_string();
    if name.is_empty() {
        return Err(FmiError::ModelDescriptionParse(
            "tag with empty element name".to_string(),
        ));
    }
    let attr_str = body[name_split..].trim();
    let attrs = parse_attrs(attr_str)?;

    Ok(Tag { name, kind, attrs })
}

/// Parse `key="value"` / `key='value'` attribute pairs from the remainder
/// of a tag. Fail-loud on a key without a quoted value.
fn parse_attrs(s: &str) -> Result<Vec<(String, String)>> {
    let mut attrs = Vec::new();
    let bytes = s.as_bytes();
    let mut i = 0usize;
    let n = bytes.len();

    while i < n {
        // Skip whitespace.
        while i < n && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= n {
            break;
        }
        // Read key up to '=' or whitespace.
        let key_start = i;
        while i < n && bytes[i] != b'=' && !bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        let key = s[key_start..i].to_string();
        if key.is_empty() {
            return Err(FmiError::ModelDescriptionParse(format!(
                "malformed attributes near {:?}",
                &s[key_start..]
            )));
        }
        // Skip whitespace then require '='.
        while i < n && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= n || bytes[i] != b'=' {
            return Err(FmiError::ModelDescriptionParse(format!(
                "attribute {key:?} has no '=' value"
            )));
        }
        i += 1; // consume '='
        while i < n && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        // Require a quote.
        if i >= n || (bytes[i] != b'"' && bytes[i] != b'\'') {
            return Err(FmiError::ModelDescriptionParse(format!(
                "attribute {key:?} value is not quoted"
            )));
        }
        let quote = bytes[i];
        i += 1;
        let val_start = i;
        while i < n && bytes[i] != quote {
            i += 1;
        }
        if i >= n {
            return Err(FmiError::ModelDescriptionParse(format!(
                "attribute {key:?} has an unterminated quoted value"
            )));
        }
        let raw = &s[val_start..i];
        i += 1; // consume closing quote
        attrs.push((key, unescape(raw)));
    }

    Ok(attrs)
}

/// Resolve the five predefined XML entities in an attribute value.
fn unescape(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    s.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        // Ampersand last so we don't double-decode the ones above.
        .replace("&amp;", "&")
}

#[cfg(test)]
mod tests {
    use super::*;

    const FMI2: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<fmiModelDescription fmiVersion="2.0" modelName="SpringDamper" guid="abc">
  <CoSimulation modelIdentifier="SpringDamper"/>
  <ModelVariables>
    <ScalarVariable name="force_in" valueReference="0" causality="input">
      <Real start="0.0"/>
    </ScalarVariable>
    <ScalarVariable name="position" valueReference="1" causality="output">
      <Real/>
    </ScalarVariable>
    <ScalarVariable name="stiffness" valueReference="2" causality="parameter">
      <Real start="100.0"/>
    </ScalarVariable>
  </ModelVariables>
</fmiModelDescription>"#;

    const FMI3: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<fmiModelDescription fmiVersion="3.0" modelName="Mass3">
  <CoSimulation modelIdentifier="Mass3"/>
  <ModelVariables>
    <Float64 name="u" valueReference="10" causality="input"/>
    <Float64 name="y" valueReference="11" causality="output"/>
    <Int32 name="mode" valueReference="12" causality="parameter"/>
  </ModelVariables>
</fmiModelDescription>"#;

    #[test]
    fn parses_fmi2_model_description() {
        let md = ModelDescription::parse(FMI2).expect("parse fmi2");
        assert_eq!(md.model_name, "SpringDamper");
        assert_eq!(md.fmi_version.as_deref(), Some("2.0"));
        assert_eq!(md.variables.len(), 3);

        assert_eq!(md.variables[0].name, "force_in");
        assert_eq!(md.variables[0].value_reference, 0);
        assert_eq!(md.variables[0].causality, Causality::Input);

        assert_eq!(md.variables[1].name, "position");
        assert_eq!(md.variables[1].value_reference, 1);
        assert_eq!(md.variables[1].causality, Causality::Output);

        assert_eq!(md.variables[2].causality, Causality::Parameter);

        assert_eq!(md.inputs(), vec!["force_in"]);
        assert_eq!(md.outputs(), vec!["position"]);
    }

    #[test]
    fn parses_fmi3_typed_variables() {
        let md = ModelDescription::parse(FMI3).expect("parse fmi3");
        assert_eq!(md.model_name, "Mass3");
        assert_eq!(md.fmi_version.as_deref(), Some("3.0"));
        assert_eq!(md.variables.len(), 3);
        assert_eq!(md.variables[0].name, "u");
        assert_eq!(md.variables[0].value_reference, 10);
        assert_eq!(md.variables[0].causality, Causality::Input);
        assert_eq!(md.outputs(), vec!["y"]);
    }

    #[test]
    fn malformed_xml_is_fail_loud_err() {
        // Unterminated tag.
        let bad = r#"<fmiModelDescription modelName="x" <ModelVariables>"#;
        assert!(matches!(
            ModelDescription::parse(bad),
            Err(FmiError::ModelDescriptionParse(_))
        ));
    }

    #[test]
    fn unterminated_quote_is_err() {
        let bad = r#"<?xml version="1.0"?><fmiModelDescription modelName="oops>
        <ModelVariables></ModelVariables></fmiModelDescription>"#;
        assert!(matches!(
            ModelDescription::parse(bad),
            Err(FmiError::ModelDescriptionParse(_))
        ));
    }

    #[test]
    fn missing_model_name_is_err() {
        let bad = r#"<fmiModelDescription fmiVersion="2.0">
          <ModelVariables></ModelVariables>
        </fmiModelDescription>"#;
        let err = ModelDescription::parse(bad).unwrap_err();
        assert!(matches!(err, FmiError::ModelDescriptionParse(ref m) if m.contains("modelName")));
    }

    #[test]
    fn missing_model_variables_is_err() {
        let bad = r#"<fmiModelDescription modelName="x"></fmiModelDescription>"#;
        assert!(matches!(
            ModelDescription::parse(bad),
            Err(FmiError::ModelDescriptionParse(_))
        ));
    }

    #[test]
    fn non_integer_value_reference_is_err() {
        let bad = r#"<fmiModelDescription modelName="x">
          <ModelVariables>
            <Float64 name="u" valueReference="not-a-number" causality="input"/>
          </ModelVariables>
        </fmiModelDescription>"#;
        assert!(matches!(
            ModelDescription::parse(bad),
            Err(FmiError::ModelDescriptionParse(_))
        ));
    }

    #[test]
    fn unknown_causality_is_err() {
        let bad = r#"<fmiModelDescription modelName="x">
          <ModelVariables>
            <Float64 name="u" valueReference="0" causality="sideways"/>
          </ModelVariables>
        </fmiModelDescription>"#;
        assert!(matches!(
            ModelDescription::parse(bad),
            Err(FmiError::ModelDescriptionParse(_))
        ));
    }

    #[test]
    fn entity_unescape_in_attribute() {
        let xml = r#"<fmiModelDescription modelName="a &amp; b">
          <ModelVariables>
            <Float64 name="x" valueReference="0" causality="local"/>
          </ModelVariables>
        </fmiModelDescription>"#;
        let md = ModelDescription::parse(xml).unwrap();
        assert_eq!(md.model_name, "a & b");
    }
}
