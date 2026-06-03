//! SBML-subset reader and writer — features 2 and 3.
//!
//! SBML (the Systems Biology Markup Language) is the lingua franca of
//! the COPASI / Tellurium / iBioSim ecosystem. The full Level-3
//! specification with its package extensions is enormous; this module
//! handles the **core subset** that round-trips a kinetic
//! reaction-network model:
//!
//! - `<listOfCompartments>` → [`Compartment`] (`id`, `size`).
//! - `<listOfSpecies>` → [`Species`] (`id`, `compartment`,
//!   `initialAmount` / `initialConcentration`, `constant` /
//!   `boundaryCondition`).
//! - `<listOfParameters>` → [`Parameter`] (`id`, `value`).
//! - `<listOfReactions>` → [`Reaction`] with `<listOfReactants>`,
//!   `<listOfProducts>` and a `<kineticLaw>`.
//!
//! Because Round 6 keeps a heavy XML crate off the dependency tree,
//! the parser is a small hand-rolled tag scanner. It is deliberately
//! *forgiving* about attribute order and whitespace but *strict*
//! about structure — a reaction referencing an undeclared species is
//! a [`SysbioError::Parse`].
//!
//! ## v1 caveats
//!
//! The kinetic law is **not** a full MathML evaluator. Instead the
//! writer emits, and the reader recognises, a compact
//! `sbml:rateLawKind` annotation naming one of the four
//! [`RateLaw`] variants plus its constants. A foreign SBML file
//! whose `<kineticLaw>` carries only `<math>` will still parse — its
//! reactions get a `Constant { rate: 0.0 }` placeholder law and the
//! reader records that in the returned [`SbmlReadReport`] so the
//! caller is never silently misled.

use std::fmt::Write as _;

use crate::error::{Result, SysbioError};
use crate::model::events::{
    AssignmentRule, EventAssignment, RateRule, SbmlEvent, SbmlRules, VarRef,
};
use crate::model::expr::Expr;
use crate::model::{Compartment, Model, Parameter, RateLaw, Reaction, Species};

/// Diagnostics from [`read_sbml`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SbmlReadReport {
    /// Number of reactions whose kinetic law could not be reconstructed
    /// from a `sbml:rateLawKind` annotation and were given a
    /// zero-flux placeholder.
    pub placeholder_laws: usize,
    /// Ids of those reactions.
    pub placeholder_reaction_ids: Vec<String>,
}

/// Parse an SBML-subset document into a [`Model`].
///
/// Returns the model plus an [`SbmlReadReport`] flagging any reaction
/// whose kinetics had to be approximated.
pub fn read_sbml(xml: &str) -> Result<(Model, SbmlReadReport)> {
    let toks = tokenize(xml)?;
    let mut model = Model {
        id: String::new(),
        compartments: Vec::new(),
        species: Vec::new(),
        reactions: Vec::new(),
        parameters: Vec::new(),
        events: Vec::new(),
        rules: crate::model::SbmlRules::default(),
    };
    let mut report = SbmlReadReport::default();

    // Pass 1: compartments, species, parameters, model id.
    for t in &toks {
        match t.name.as_str() {
            "model" => {
                if let Some(id) = t.attr("id") {
                    model.id = id.to_string();
                }
            }
            "compartment" => {
                let id = t.attr("id").ok_or_else(|| {
                    SysbioError::parse("sbml", "<compartment> without id")
                })?;
                let size = t
                    .attr("size")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(1.0);
                model.compartments.push(Compartment {
                    id: id.to_string(),
                    size,
                });
            }
            _ => {}
        }
    }
    if model.compartments.is_empty() {
        model.compartments.push(Compartment::new("default"));
    }
    for t in &toks {
        if t.name == "species" {
            let id = t
                .attr("id")
                .ok_or_else(|| SysbioError::parse("sbml", "<species> without id"))?;
            let comp = t
                .attr("compartment")
                .and_then(|c| model.compartments.iter().position(|x| x.id == c))
                .unwrap_or(0);
            let initial = t
                .attr("initialAmount")
                .or_else(|| t.attr("initialConcentration"))
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0);
            let constant = t.attr("constant") == Some("true")
                || t.attr("boundaryCondition") == Some("true");
            model.species.push(Species {
                id: id.to_string(),
                compartment: comp,
                initial,
                constant,
            });
        } else if t.name == "parameter" {
            // Skip parameters that live inside a <kineticLaw> — those
            // are local and handled with the reaction. A top-level
            // listOfParameters parameter has no enclosing reaction.
            if t.in_kinetic_law {
                continue;
            }
            let id = t
                .attr("id")
                .ok_or_else(|| SysbioError::parse("sbml", "<parameter> without id"))?;
            let value = t
                .attr("value")
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0);
            model.parameters.push(Parameter {
                id: id.to_string(),
                value,
            });
        }
    }

    // Pass 2: reactions. Reactant / product / kineticLaw tokens carry
    // the enclosing reaction id (stamped by the tokenizer).
    let resolve = |sid: &str| -> Result<usize> {
        model
            .species
            .iter()
            .position(|s| s.id == sid)
            .ok_or_else(|| {
                SysbioError::parse("sbml", format!("reaction references unknown species `{sid}`"))
            })
    };
    let mut current: Option<Reaction> = None;
    let mut reactions: Vec<Reaction> = Vec::new();
    for t in &toks {
        match t.name.as_str() {
            "reaction" => {
                if let Some(r) = current.take() {
                    reactions.push(r);
                }
                let id = t
                    .attr("id")
                    .ok_or_else(|| SysbioError::parse("sbml", "<reaction> without id"))?;
                current = Some(Reaction {
                    id: id.to_string(),
                    reactants: Vec::new(),
                    products: Vec::new(),
                    rate_law: RateLaw::Constant { rate: 0.0 },
                    reversible: t.attr("reversible") != Some("false"),
                });
            }
            "speciesReference" => {
                let Some(r) = current.as_mut() else { continue };
                let sid = t.attr("species").ok_or_else(|| {
                    SysbioError::parse("sbml", "<speciesReference> without species")
                })?;
                let idx = resolve(sid)?;
                let stoich = t
                    .attr("stoichiometry")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(1.0);
                if t.side_is_product {
                    r.products.push((idx, stoich));
                } else {
                    r.reactants.push((idx, stoich));
                }
            }
            "kineticLaw" => {
                let Some(r) = current.as_mut() else { continue };
                match parse_rate_law_annotation(t, &model) {
                    Some(law) => r.rate_law = law,
                    None => {
                        report.placeholder_laws += 1;
                        report.placeholder_reaction_ids.push(r.id.clone());
                    }
                }
            }
            _ => {}
        }
    }
    if let Some(r) = current.take() {
        reactions.push(r);
    }
    model.reactions = reactions;

    // Pass 3: events and rules. The compact-annotation encoding
    // (`sbml:expr`, `sbml:target`, …) is symmetrical with the
    // rate-law approach: serialised faithfully by [`write_sbml`] and
    // round-tripped by this reader without going through a MathML
    // evaluator.
    let resolve_var = |target: &str| -> Result<VarRef> {
        if let Some(rest) = target.strip_prefix("species:") {
            return model
                .species_index(rest)
                .map(VarRef::Species)
                .ok_or_else(|| {
                    SysbioError::parse("sbml", format!("unknown species `{rest}` in target"))
                });
        }
        if let Some(rest) = target.strip_prefix("parameter:") {
            return model
                .parameter_index(rest)
                .map(VarRef::Parameter)
                .ok_or_else(|| {
                    SysbioError::parse("sbml", format!("unknown parameter `{rest}` in target"))
                });
        }
        Err(SysbioError::parse(
            "sbml",
            format!("malformed target `{target}` (expected `species:<id>` or `parameter:<id>`)"),
        ))
    };

    let mut events: Vec<SbmlEvent> = Vec::new();
    let mut current_ev: Option<SbmlEvent> = None;
    let mut current_ev_assignments: Vec<EventAssignment> = Vec::new();
    let mut rules = SbmlRules::default();
    for t in &toks {
        match t.name.as_str() {
            "event" => {
                if let Some(mut ev) = current_ev.take() {
                    ev.assignments = std::mem::take(&mut current_ev_assignments);
                    events.push(ev);
                }
                let id = t
                    .attr("id")
                    .ok_or_else(|| SysbioError::parse("sbml", "<event> without id"))?
                    .to_string();
                let trig_str = t.attr("sbml:trigger").ok_or_else(|| {
                    SysbioError::parse("sbml", format!("event `{id}` missing sbml:trigger"))
                })?;
                let trigger = Expr::parse(trig_str).ok_or_else(|| {
                    SysbioError::parse(
                        "sbml",
                        format!("event `{id}` has unparseable trigger `{trig_str}`"),
                    )
                })?;
                let delay = t.attr("sbml:delay").and_then(|v| v.parse().ok());
                let priority = t
                    .attr("sbml:priority")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0.0);
                let use_trigger_values = t.attr("sbml:useValuesFromTriggerTime") != Some("false");
                current_ev = Some(SbmlEvent {
                    id,
                    trigger,
                    assignments: Vec::new(),
                    delay,
                    priority,
                    use_values_from_trigger_time: use_trigger_values,
                    initial_value: false,
                });
            }
            "eventAssignment" => {
                let target = t.attr("sbml:target").ok_or_else(|| {
                    SysbioError::parse("sbml", "<eventAssignment> missing sbml:target")
                })?;
                let formula_str = t.attr("sbml:expr").ok_or_else(|| {
                    SysbioError::parse("sbml", "<eventAssignment> missing sbml:expr")
                })?;
                let target_ref = resolve_var(target)?;
                let formula = Expr::parse(formula_str).ok_or_else(|| {
                    SysbioError::parse(
                        "sbml",
                        format!("unparseable assignment formula `{formula_str}`"),
                    )
                })?;
                current_ev_assignments.push(EventAssignment {
                    target: target_ref,
                    formula,
                });
            }
            "assignmentRule" => {
                let target = t.attr("sbml:target").ok_or_else(|| {
                    SysbioError::parse("sbml", "<assignmentRule> missing sbml:target")
                })?;
                let formula_str = t.attr("sbml:expr").ok_or_else(|| {
                    SysbioError::parse("sbml", "<assignmentRule> missing sbml:expr")
                })?;
                let target_ref = resolve_var(target)?;
                let formula = Expr::parse(formula_str).ok_or_else(|| {
                    SysbioError::parse(
                        "sbml",
                        format!("unparseable assignment-rule formula `{formula_str}`"),
                    )
                })?;
                rules.assignments.push(AssignmentRule {
                    target: target_ref,
                    formula,
                });
            }
            "rateRule" => {
                let target = t.attr("sbml:target").ok_or_else(|| {
                    SysbioError::parse("sbml", "<rateRule> missing sbml:target")
                })?;
                let formula_str = t.attr("sbml:expr").ok_or_else(|| {
                    SysbioError::parse("sbml", "<rateRule> missing sbml:expr")
                })?;
                let target_ref = resolve_var(target)?;
                let formula = Expr::parse(formula_str).ok_or_else(|| {
                    SysbioError::parse(
                        "sbml",
                        format!("unparseable rate-rule formula `{formula_str}`"),
                    )
                })?;
                rules.rates.push(RateRule {
                    target: target_ref,
                    formula,
                });
            }
            _ => {}
        }
    }
    if let Some(mut ev) = current_ev.take() {
        ev.assignments = std::mem::take(&mut current_ev_assignments);
        events.push(ev);
    }
    model.events = events;
    model.rules = rules;

    model.validate().map_err(|e| {
        SysbioError::parse("sbml", format!("document parsed but model invalid: {e}"))
    })?;
    Ok((model, report))
}

/// Reconstruct a [`RateLaw`] from the `sbml:rateLawKind` attribute set
/// that [`write_sbml`] emits. Returns `None` for a foreign file with
/// only `<math>` content.
fn parse_rate_law_annotation(t: &Token, model: &Model) -> Option<RateLaw> {
    let kind = t.attr("sbml:rateLawKind")?;
    let f = |k: &str| t.attr(k).and_then(|v| v.parse::<f64>().ok());
    let sp = |k: &str| {
        t.attr(k)
            .and_then(|id| model.species.iter().position(|s| s.id == id))
    };
    match kind {
        "constant" => Some(RateLaw::Constant { rate: f("rate")? }),
        "mass_action" => {
            let k = f("k")?;
            // Reactant indices/orders encoded as "id:order;id:order".
            let mut reactants = Vec::new();
            if let Some(spec) = t.attr("reactants") {
                for part in spec.split(';').filter(|p| !p.is_empty()) {
                    let (id, ord) = part.split_once(':')?;
                    let idx = model.species.iter().position(|s| s.id == id)?;
                    reactants.push((idx, ord.parse().ok()?));
                }
            }
            Some(RateLaw::MassAction { k, reactants })
        }
        "michaelis_menten" => Some(RateLaw::MichaelisMenten {
            vmax: f("vmax")?,
            km: f("km")?,
            substrate: sp("substrate")?,
        }),
        "hill" => Some(RateLaw::Hill {
            vmax: f("vmax")?,
            kd: f("kd")?,
            n: f("n")?,
            regulator: sp("regulator")?,
            repress: t.attr("repress") == Some("true"),
        }),
        _ => None,
    }
}

/// Serialise a [`Model`] back to an SBML-subset document.
///
/// The output is valid SBML Level-3 core structure; the kinetic laws
/// are carried in `sbml:rateLawKind` attributes (plus, for human
/// readers and foreign tools, a best-effort `<math>` placeholder is
/// *not* emitted — the annotation is the source of truth and
/// [`read_sbml`] reads it back exactly).
pub fn write_sbml(model: &Model) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(
        "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version2/core\" \
         level=\"3\" version=\"2\">\n",
    );
    let _ = writeln!(out, "  <model id=\"{}\">", esc(&model.id));

    out.push_str("    <listOfCompartments>\n");
    for c in &model.compartments {
        let _ = writeln!(
            out,
            "      <compartment id=\"{}\" size=\"{}\" constant=\"true\"/>",
            esc(&c.id),
            c.size
        );
    }
    out.push_str("    </listOfCompartments>\n");

    out.push_str("    <listOfSpecies>\n");
    for s in &model.species {
        let comp = model
            .compartments
            .get(s.compartment)
            .map(|c| c.id.as_str())
            .unwrap_or("default");
        let _ = writeln!(
            out,
            "      <species id=\"{}\" compartment=\"{}\" initialAmount=\"{}\" \
             hasOnlySubstanceUnits=\"true\" boundaryCondition=\"{}\" constant=\"{}\"/>",
            esc(&s.id),
            esc(comp),
            s.initial,
            s.constant,
            s.constant
        );
    }
    out.push_str("    </listOfSpecies>\n");

    if !model.parameters.is_empty() {
        out.push_str("    <listOfParameters>\n");
        for p in &model.parameters {
            let _ = writeln!(
                out,
                "      <parameter id=\"{}\" value=\"{}\" constant=\"true\"/>",
                esc(&p.id),
                p.value
            );
        }
        out.push_str("    </listOfParameters>\n");
    }

    out.push_str("    <listOfReactions>\n");
    for r in &model.reactions {
        let _ = writeln!(
            out,
            "      <reaction id=\"{}\" reversible=\"{}\" fast=\"false\">",
            esc(&r.id),
            r.reversible
        );
        if !r.reactants.is_empty() {
            out.push_str("        <listOfReactants>\n");
            for &(i, c) in &r.reactants {
                emit_species_ref(&mut out, model, i, c);
            }
            out.push_str("        </listOfReactants>\n");
        }
        if !r.products.is_empty() {
            out.push_str("        <listOfProducts>\n");
            for &(i, c) in &r.products {
                emit_species_ref(&mut out, model, i, c);
            }
            out.push_str("        </listOfProducts>\n");
        }
        out.push_str("        <kineticLaw");
        emit_rate_law_attrs(&mut out, model, &r.rate_law);
        out.push_str("/>\n");
        out.push_str("      </reaction>\n");
    }
    out.push_str("    </listOfReactions>\n");

    // SBML L3 rules (assignment + rate). The compact `sbml:expr`
    // annotation carries the formula in the same parenthesised infix
    // form the `Expr::to_string_compact` writer produces; targets are
    // tagged `species:<id>` or `parameter:<id>` so the reader can
    // disambiguate without consulting the model.
    if !model.rules.assignments.is_empty() || !model.rules.rates.is_empty() {
        out.push_str("    <listOfRules>\n");
        for rule in &model.rules.assignments {
            let _ = writeln!(
                out,
                "      <assignmentRule sbml:target=\"{}\" sbml:expr=\"{}\"/>",
                esc(&target_label(model, &rule.target)),
                esc(&rule.formula.to_string_compact()),
            );
        }
        for rule in &model.rules.rates {
            let _ = writeln!(
                out,
                "      <rateRule sbml:target=\"{}\" sbml:expr=\"{}\"/>",
                esc(&target_label(model, &rule.target)),
                esc(&rule.formula.to_string_compact()),
            );
        }
        out.push_str("    </listOfRules>\n");
    }

    // SBML L3 events.
    if !model.events.is_empty() {
        out.push_str("    <listOfEvents>\n");
        for ev in &model.events {
            let _ = write!(
                out,
                "      <event id=\"{}\" sbml:trigger=\"{}\" sbml:priority=\"{}\" \
                 sbml:useValuesFromTriggerTime=\"{}\"",
                esc(&ev.id),
                esc(&ev.trigger.to_string_compact()),
                ev.priority,
                ev.use_values_from_trigger_time,
            );
            if let Some(d) = ev.delay {
                let _ = write!(out, " sbml:delay=\"{d}\"");
            }
            out.push_str(">\n");
            out.push_str("        <listOfEventAssignments>\n");
            for a in &ev.assignments {
                let _ = writeln!(
                    out,
                    "          <eventAssignment sbml:target=\"{}\" sbml:expr=\"{}\"/>",
                    esc(&target_label(model, &a.target)),
                    esc(&a.formula.to_string_compact()),
                );
            }
            out.push_str("        </listOfEventAssignments>\n");
            out.push_str("      </event>\n");
        }
        out.push_str("    </listOfEvents>\n");
    }

    out.push_str("  </model>\n");
    out.push_str("</sbml>\n");
    out
}

/// Format a [`VarRef`] for the `sbml:target` attribute - `species:<id>`
/// or `parameter:<id>`.
fn target_label(model: &Model, vref: &VarRef) -> String {
    match vref {
        VarRef::Species(i) => format!(
            "species:{}",
            model.species.get(*i).map(|s| s.id.as_str()).unwrap_or("?")
        ),
        VarRef::Parameter(i) => format!(
            "parameter:{}",
            model
                .parameters
                .get(*i)
                .map(|p| p.id.as_str())
                .unwrap_or("?")
        ),
    }
}

fn emit_species_ref(out: &mut String, model: &Model, idx: usize, coeff: f64) {
    let id = model
        .species
        .get(idx)
        .map(|s| s.id.as_str())
        .unwrap_or("?");
    let _ = writeln!(
        out,
        "          <speciesReference species=\"{}\" stoichiometry=\"{}\" constant=\"true\"/>",
        esc(id),
        coeff
    );
}

fn emit_rate_law_attrs(out: &mut String, model: &Model, law: &RateLaw) {
    let sid = |i: usize| {
        model
            .species
            .get(i)
            .map(|s| s.id.clone())
            .unwrap_or_else(|| "?".into())
    };
    match law {
        RateLaw::Constant { rate } => {
            let _ = write!(out, " sbml:rateLawKind=\"constant\" rate=\"{rate}\"");
        }
        RateLaw::MassAction { k, reactants } => {
            let spec: String = reactants
                .iter()
                .map(|&(i, o)| format!("{}:{};", sid(i), o))
                .collect();
            let _ = write!(
                out,
                " sbml:rateLawKind=\"mass_action\" k=\"{k}\" reactants=\"{}\"",
                esc(&spec)
            );
        }
        RateLaw::MichaelisMenten {
            vmax,
            km,
            substrate,
        } => {
            let _ = write!(
                out,
                " sbml:rateLawKind=\"michaelis_menten\" vmax=\"{vmax}\" km=\"{km}\" \
                 substrate=\"{}\"",
                esc(&sid(*substrate))
            );
        }
        RateLaw::Hill {
            vmax,
            kd,
            n,
            regulator,
            repress,
        } => {
            let _ = write!(
                out,
                " sbml:rateLawKind=\"hill\" vmax=\"{vmax}\" kd=\"{kd}\" n=\"{n}\" \
                 regulator=\"{}\" repress=\"{repress}\"",
                esc(&sid(*regulator))
            );
        }
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn unesc(s: &str) -> String {
    s.replace("&quot;", "\"")
        .replace("&gt;", ">")
        .replace("&lt;", "<")
        .replace("&amp;", "&")
}

// --- tiny XML tag scanner --------------------------------------------

/// One parsed start / empty XML tag plus position context.
struct Token {
    name: String,
    attrs: Vec<(String, String)>,
    /// `true` if this tag lay inside a `<listOfProducts>`.
    side_is_product: bool,
    /// `true` if this tag lay inside a `<kineticLaw>` element.
    in_kinetic_law: bool,
}

impl Token {
    fn attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

/// Scan an XML string into a flat list of start / empty tags. Closing
/// tags are consumed only to maintain the `products` / `kineticLaw`
/// context flags. Comments, processing instructions and `<math>`
/// subtrees are skipped wholesale.
fn tokenize(xml: &str) -> Result<Vec<Token>> {
    let bytes = xml.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    let mut in_products = false;
    let mut kinetic_depth = 0usize;

    while i < bytes.len() {
        if bytes[i] != b'<' {
            i += 1;
            continue;
        }
        // Comment / CDATA / processing instruction / declaration.
        if xml[i..].starts_with("<!--") {
            let end = xml[i..]
                .find("-->")
                .ok_or_else(|| SysbioError::parse("sbml", "unterminated comment"))?;
            i += end + 3;
            continue;
        }
        if xml[i..].starts_with("<?") {
            let end = xml[i..]
                .find("?>")
                .ok_or_else(|| SysbioError::parse("sbml", "unterminated processing instruction"))?;
            i += end + 2;
            continue;
        }
        // Find the matching '>'.
        let close = xml[i..]
            .find('>')
            .ok_or_else(|| SysbioError::parse("sbml", "unterminated tag"))?;
        let raw = &xml[i + 1..i + close];
        i += close + 1;

        if let Some(stripped) = raw.strip_prefix('/') {
            // Closing tag — update context.
            let name = stripped.trim();
            if name == "listOfProducts" || name == "listOfReactants" {
                in_products = false;
            } else if name == "kineticLaw" {
                kinetic_depth = kinetic_depth.saturating_sub(1);
            }
            continue;
        }
        if raw.starts_with('!') {
            continue; // <!DOCTYPE …> etc.
        }
        let empty = raw.ends_with('/');
        let body = raw.trim_end_matches('/').trim();
        let (name, attr_str) = match body.find(char::is_whitespace) {
            Some(p) => (&body[..p], &body[p..]),
            None => (body, ""),
        };
        // Strip any namespace prefix for the *element* name only.
        let local = name.rsplit(':').next().unwrap_or(name).to_string();
        let attrs = parse_attrs(attr_str)?;
        if local == "listOfProducts" {
            in_products = true;
        } else if local == "listOfReactants" {
            in_products = false;
        }
        let token = Token {
            name: local.clone(),
            attrs,
            side_is_product: in_products,
            in_kinetic_law: kinetic_depth > 0,
        };
        if local == "kineticLaw" && !empty {
            kinetic_depth += 1;
        }
        out.push(token);
    }
    if out.is_empty() {
        return Err(SysbioError::parse("sbml", "no XML tags found"));
    }
    Ok(out)
}

/// Parse `key="value"` pairs out of an attribute string.
fn parse_attrs(s: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        if chars[i].is_whitespace() {
            i += 1;
            continue;
        }
        let key_start = i;
        while i < chars.len() && chars[i] != '=' && !chars[i].is_whitespace() {
            i += 1;
        }
        let key: String = chars[key_start..i].iter().collect();
        while i < chars.len() && chars[i].is_whitespace() {
            i += 1;
        }
        if i >= chars.len() || chars[i] != '=' {
            // Valueless attribute — skip it.
            continue;
        }
        i += 1; // consume '='
        while i < chars.len() && chars[i].is_whitespace() {
            i += 1;
        }
        if i >= chars.len() || (chars[i] != '"' && chars[i] != '\'') {
            return Err(SysbioError::parse(
                "sbml",
                format!("attribute `{key}` has no quoted value"),
            ));
        }
        let quote = chars[i];
        i += 1;
        let val_start = i;
        while i < chars.len() && chars[i] != quote {
            i += 1;
        }
        if i >= chars.len() {
            return Err(SysbioError::parse(
                "sbml",
                format!("attribute `{key}` value not closed"),
            ));
        }
        let val: String = chars[val_start..i].iter().collect();
        i += 1; // consume closing quote
        out.push((key, unesc(&val)));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_model() -> Model {
        let mut m = Model::new("glycolysis_toy");
        m.compartments[0].size = 2.0;
        let a = m.add_species(Species::new("A", 10.0));
        let b = m.add_species(Species::new("B", 0.0));
        let e = m.add_species(Species::new("E", 1.0).constant());
        m.add_parameter(Parameter::new("k1", 0.3));
        m.add_reaction(Reaction {
            id: "conv".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MichaelisMenten {
                vmax: 5.0,
                km: 2.0,
                substrate: a,
            },
            reversible: false,
        });
        let _ = e;
        m
    }

    #[test]
    fn write_then_read_roundtrips() {
        let m = sample_model();
        let xml = write_sbml(&m);
        let (back, report) = read_sbml(&xml).expect("re-parse");
        assert_eq!(report.placeholder_laws, 0);
        assert_eq!(back.id, m.id);
        assert_eq!(back.species.len(), m.species.len());
        assert_eq!(back.reactions.len(), 1);
        assert_eq!(back.compartments[0].size, 2.0);
        // The MM law survived the round-trip.
        match &back.reactions[0].rate_law {
            RateLaw::MichaelisMenten { vmax, km, .. } => {
                assert_eq!(*vmax, 5.0);
                assert_eq!(*km, 2.0);
            }
            other => panic!("rate law not preserved: {other:?}"),
        }
        // Boundary species kept its constant flag.
        assert!(back.species.iter().any(|s| s.id == "E" && s.constant));
    }

    #[test]
    fn mass_action_roundtrips_with_reactant_orders() {
        let mut m = Model::new("ma");
        let a = m.add_species(Species::new("A", 4.0));
        let b = m.add_species(Species::new("B", 0.0));
        m.add_reaction(Reaction {
            id: "r".into(),
            reactants: vec![(a, 2.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 0.7,
                reactants: vec![(a, 2.0)],
            },
            reversible: false,
        });
        let (back, _) = read_sbml(&write_sbml(&m)).unwrap();
        match &back.reactions[0].rate_law {
            RateLaw::MassAction { k, reactants } => {
                assert!((k - 0.7).abs() < 1e-12);
                assert_eq!(reactants, &vec![(0usize, 2.0)]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn foreign_kinetic_law_gets_placeholder() {
        // A <kineticLaw> with only <math> — no rateLawKind annotation.
        let xml = r#"<?xml version="1.0"?>
        <sbml level="3" version="2">
          <model id="foreign">
            <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
            <listOfSpecies>
              <species id="X" compartment="c" initialAmount="1"/>
              <species id="Y" compartment="c" initialAmount="0"/>
            </listOfSpecies>
            <listOfReactions>
              <reaction id="r1">
                <listOfReactants>
                  <speciesReference species="X" stoichiometry="1"/>
                </listOfReactants>
                <listOfProducts>
                  <speciesReference species="Y" stoichiometry="1"/>
                </listOfProducts>
                <kineticLaw>
                  <math xmlns="http://www.w3.org/1998/Math/MathML"><ci>X</ci></math>
                </kineticLaw>
              </reaction>
            </listOfReactions>
          </model>
        </sbml>"#;
        let (model, report) = read_sbml(xml).expect("parse foreign");
        assert_eq!(model.reactions.len(), 1);
        assert_eq!(report.placeholder_laws, 1);
        assert_eq!(report.placeholder_reaction_ids, vec!["r1"]);
    }

    #[test]
    fn reaction_with_unknown_species_is_rejected() {
        let xml = r#"<sbml><model id="m">
          <listOfSpecies><species id="A" compartment="c" initialAmount="1"/></listOfSpecies>
          <listOfReactions><reaction id="r">
            <listOfReactants><speciesReference species="GHOST"/></listOfReactants>
          </reaction></listOfReactions>
        </model></sbml>"#;
        assert!(read_sbml(xml).is_err());
    }

    #[test]
    fn empty_document_errors() {
        assert!(read_sbml("").is_err());
        assert!(read_sbml("   \n  ").is_err());
    }

    #[test]
    fn events_round_trip_through_sbml() {
        // Build a model with one event that has a delay, priority,
        // multiple assignments and a state-based trigger.
        let mut m = Model::new("evt");
        let a = m.add_species(Species::new("A", 10.0));
        let _b = m.add_species(Species::new("B", 0.0));
        m.add_parameter(Parameter::new("k", 0.7));
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        m.add_event(
            SbmlEvent::new(
                "pulse",
                Expr::and(
                    Expr::ge(Expr::Time, Expr::k(2.0)),
                    Expr::lt(Expr::var(a), Expr::k(100.0)),
                ),
                vec![
                    EventAssignment {
                        target: VarRef::Species(a),
                        formula: Expr::k(50.0),
                    },
                    EventAssignment {
                        target: VarRef::Parameter(0),
                        formula: Expr::mul(Expr::param(0), Expr::k(2.0)),
                    },
                ],
            )
            .with_delay(0.5)
            .with_priority(3.5)
            .evaluate_at_execution_time(),
        );
        let xml = write_sbml(&m);
        let (back, _) = read_sbml(&xml).expect("re-parse");
        assert_eq!(back.events.len(), 1);
        let ev = &back.events[0];
        assert_eq!(ev.id, "pulse");
        assert_eq!(ev.delay, Some(0.5));
        assert_eq!(ev.priority, 3.5);
        assert!(!ev.use_values_from_trigger_time);
        assert_eq!(ev.assignments.len(), 2);
        // Trigger and assignment formulas evaluate identically before
        // and after the round trip.
        let y = vec![5.0, 0.0];
        let p = vec![0.7];
        assert_eq!(
            ev.trigger.value(&y, &p, 3.0),
            m.events[0].trigger.value(&y, &p, 3.0)
        );
        assert_eq!(
            ev.assignments[1].formula.value(&y, &p, 3.0),
            m.events[0].assignments[1].formula.value(&y, &p, 3.0)
        );
    }

    #[test]
    fn rules_round_trip_through_sbml() {
        let mut m = Model::new("rules");
        let _a = m.add_species(Species::new("A", 0.0));
        let _b = m.add_species(Species::new("B", 0.0));
        m.add_parameter(Parameter::new("k", 1.5));
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        m.add_rate_rule(RateRule {
            target: VarRef::Species(0),
            formula: Expr::param(0),
        });
        m.add_assignment_rule(AssignmentRule {
            target: VarRef::Species(1),
            formula: Expr::mul(Expr::var(0), Expr::k(2.0)),
        });
        let xml = write_sbml(&m);
        let (back, _) = read_sbml(&xml).unwrap();
        assert_eq!(back.rules.rates.len(), 1);
        assert_eq!(back.rules.assignments.len(), 1);
        assert!(matches!(
            back.rules.rates[0].target,
            VarRef::Species(0)
        ));
        assert!(matches!(
            back.rules.assignments[0].target,
            VarRef::Species(1)
        ));
    }
}
