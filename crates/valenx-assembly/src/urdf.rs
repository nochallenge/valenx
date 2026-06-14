//! Native **URDF importer** — turn a Unified Robot Description Format file
//! into a valenx [`Assembly`] of parts and joints, entirely in-house.
//!
//! This is the no-adapter path: a robot described in URDF (the standard ROS
//! robot-description format) is parsed and mapped directly onto valenx's own
//! assembly model — each `<link>` becomes a [`Part`] with primitive geometry,
//! and each `<joint>` becomes a valenx [`Joint`]:
//!
//! - `revolute` / `continuous` → [`JointKind::Revolute`]
//! - `prismatic` → [`JointKind::Prismatic`]
//! - `fixed` → [`JointKind::Fixed`]
//! - `planar` → [`JointKind::Planar`]
//! - anything else → `Fixed` with a recorded warning.
//!
//! The link tree is then forward-kinematically placed (each child posed by its
//! joint origin relative to its parent) so the imported robot renders already
//! assembled, and the root link is marked `fixed`. The result is a live valenx
//! assembly: its joint parameters can be driven by the assembly kinematics to
//! pose the robot — no MuJoCo, no external tool, no extra dependency (a small
//! in-crate XML reader does the parsing).
//!
//! ## Honest scope
//!
//! A **first-order kinematic import**: link geometry is mapped to the nearest
//! valenx-cad primitive (`box` / `cylinder` / `sphere`; `mesh` and unknown
//! geometry fall back to a small placeholder box with a warning), per-`<visual>`
//! origin offsets *within* a link are not baked into the solid, and URDF
//! features beyond the kinematic tree (joint limits, `mimic`, dynamics,
//! materials, `safety_controller`) are not modelled. It captures the robot's
//! structure and kinematics — the parts, the joints, and how they move.

use std::collections::{HashMap, HashSet, VecDeque};

use nalgebra::{UnitQuaternion, Vector3};

use crate::joint::{Joint, JointKind};
use crate::{Assembly, Part, PartTransform};

/// Failure modes of [`import_urdf`].
#[derive(Debug, thiserror::Error)]
pub enum UrdfError {
    /// The XML was malformed (unbalanced tags, unexpected end of input, …).
    #[error("malformed URDF XML: {0}")]
    Xml(String),
    /// The document's root element is not `<robot>`.
    #[error("URDF root element is not <robot>")]
    NoRobot,
    /// A `<link>` element has no `name` attribute.
    #[error("a <link> element is missing its name")]
    UnnamedLink,
    /// A joint references a link that no `<link>` element defines.
    #[error("joint `{joint}` references unknown link `{link}`")]
    UnknownLink {
        /// The offending joint's name.
        joint: String,
        /// The undefined link name it referenced.
        link: String,
    },
    /// A link's geometry could not be realized as a valenx-cad solid.
    #[error("link `{link}` geometry is invalid: {reason}")]
    Geometry {
        /// The link whose geometry failed.
        link: String,
        /// Why the solid could not be built.
        reason: String,
    },
}

/// A robot imported from URDF: the native assembly plus the lookup tables that
/// connect it back to the URDF names.
pub struct UrdfRobot {
    /// The robot's `name` attribute.
    pub name: String,
    /// The native assembly — parts (links) + joints, forward-kinematically
    /// posed, with the root link marked `fixed`.
    pub assembly: Assembly,
    /// Map from URDF link name to the assigned [`Part`] id.
    pub link_ids: HashMap<String, usize>,
    /// Joint names, in document order (parallel to `assembly.joints`).
    pub joint_names: Vec<String>,
    /// Non-fatal notes (unsupported geometry/joint types mapped to fallbacks).
    pub warnings: Vec<String>,
}

/// Parse a URDF document and build a native valenx [`Assembly`].
///
/// # Errors
///
/// Returns [`UrdfError`] if the XML is malformed, the root is not `<robot>`,
/// a link is unnamed, a joint references an undefined link, or a link's
/// geometry cannot be built into a valenx-cad solid.
pub fn import_urdf(xml: &str) -> Result<UrdfRobot, UrdfError> {
    let root = parse_xml(xml)?;
    if root.tag != "robot" {
        return Err(UrdfError::NoRobot);
    }
    let name = root.attr("name").unwrap_or("urdf_robot").to_string();

    let mut assembly = Assembly::new();
    let mut link_ids: HashMap<String, usize> = HashMap::new();
    let mut warnings: Vec<String> = Vec::new();

    // 1. Links → Parts.
    for lk in root.children.iter().filter(|c| c.tag == "link") {
        let lname = lk.attr("name").ok_or(UrdfError::UnnamedLink)?;
        let solid = link_solid(lname, lk, &mut warnings)?;
        let pid = assembly.add_part(Part::new(0, lname, solid));
        link_ids.insert(lname.to_string(), pid);
    }

    // 2. Joints → Joints. Collect the parent→child edges for forward kinematics.
    let mut joint_names: Vec<String> = Vec::new();
    let mut children: HashSet<String> = HashSet::new();
    let mut edges: HashMap<String, Vec<(String, PartTransform)>> = HashMap::new();

    for j in root.children.iter().filter(|c| c.tag == "joint") {
        let jname = j.attr("name").unwrap_or("joint").to_string();
        let jtype = j.attr("type").unwrap_or("fixed");
        let parent = j
            .child("parent")
            .and_then(|n| n.attr("link"))
            .ok_or_else(|| UrdfError::Xml(format!("joint `{jname}` has no <parent link=…>")))?
            .to_string();
        let child = j
            .child("child")
            .and_then(|n| n.attr("link"))
            .ok_or_else(|| UrdfError::Xml(format!("joint `{jname}` has no <child link=…>")))?
            .to_string();

        let pa = *link_ids
            .get(&parent)
            .ok_or_else(|| UrdfError::UnknownLink {
                joint: jname.clone(),
                link: parent.clone(),
            })?;
        let pb = *link_ids.get(&child).ok_or_else(|| UrdfError::UnknownLink {
            joint: jname.clone(),
            link: child.clone(),
        })?;

        let origin = j.child("origin");
        let xyz = origin.and_then(|o| o.attr("xyz")).map_or([0.0; 3], parse3);
        let rpy = origin.and_then(|o| o.attr("rpy")).map_or([0.0; 3], parse3);
        // URDF default joint axis is (1, 0, 0); express it in the parent frame.
        let axis = j
            .child("axis")
            .and_then(|a| a.attr("xyz"))
            .map_or([1.0, 0.0, 0.0], parse3);
        let q = UnitQuaternion::from_euler_angles(rpy[0], rpy[1], rpy[2]);
        let axis_dir = q * Vector3::new(axis[0], axis[1], axis[2]);
        let axis_origin = Vector3::new(xyz[0], xyz[1], xyz[2]);

        let kind = match jtype {
            "revolute" | "continuous" => JointKind::Revolute {
                part_a: pa,
                part_b: pb,
                axis_origin,
                axis_dir,
            },
            "prismatic" => JointKind::Prismatic {
                part_a: pa,
                part_b: pb,
                axis_dir,
            },
            "fixed" => JointKind::Fixed {
                part_a: pa,
                part_b: pb,
            },
            "planar" => JointKind::Planar {
                part_a: pa,
                part_b: pb,
                plane_origin: axis_origin,
                plane_normal: axis_dir,
            },
            other => {
                warnings.push(format!("joint `{jname}`: type `{other}` mapped to Fixed"));
                JointKind::Fixed {
                    part_a: pa,
                    part_b: pb,
                }
            }
        };
        assembly.add_joint(Joint::new(0, kind));
        joint_names.push(jname);
        children.insert(child.clone());
        edges
            .entry(parent)
            .or_default()
            .push((child, rigid(xyz, rpy)));
    }

    // 3. Forward kinematics: pose every link from each root (a link that is
    //    never a child). Roots are the grounded "world" parts.
    let mut pose: HashMap<String, PartTransform> = HashMap::new();
    let mut queue: VecDeque<String> = VecDeque::new();
    for (lname, &pid) in &link_ids {
        if !children.contains(lname) {
            pose.insert(lname.clone(), PartTransform::identity());
            if let Ok(p) = assembly.get_part_mut(pid) {
                p.fixed = true;
                p.transform = PartTransform::identity();
            }
            queue.push_back(lname.clone());
        }
    }
    while let Some(pname) = queue.pop_front() {
        let pt = pose
            .get(&pname)
            .cloned()
            .unwrap_or_else(PartTransform::identity);
        if let Some(kids) = edges.get(&pname) {
            for (child, local) in kids {
                let ct = compose(&pt, local);
                if let Some(&cid) = link_ids.get(child) {
                    if let Ok(p) = assembly.get_part_mut(cid) {
                        p.transform = ct.clone();
                    }
                }
                pose.insert(child.clone(), ct);
                queue.push_back(child.clone());
            }
        }
    }

    Ok(UrdfRobot {
        name,
        assembly,
        link_ids,
        joint_names,
        warnings,
    })
}

/// Tessellate every part of an assembly at its current pose and merge the
/// triangles into one [`valenx_mesh::Mesh`] — a renderable, exportable mesh of
/// the whole robot, built entirely in-house.
pub fn assembly_to_mesh(assembly: &Assembly, tessellation_tolerance: f64) -> valenx_mesh::Mesh {
    use valenx_mesh::{ElementBlock, ElementType, Mesh};

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<u32> = Vec::new();
    for part in &assembly.parts {
        let Ok(m) = valenx_cad::solid_to_mesh(&part.solid, tessellation_tolerance) else {
            continue;
        };
        let base = nodes.len() as u32;
        for n in &m.nodes {
            nodes.push(part.transform.apply_point(*n));
        }
        for b in &m.element_blocks {
            for &idx in &b.connectivity {
                tris.push(base + idx);
            }
        }
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris;
    let mut mesh = Mesh::new("urdf-robot");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

/// Build a valenx-cad solid for a link's first `<visual>`/`<collision>`
/// geometry, falling back to a small box (with a warning) for mesh/unknown
/// geometry or a geometry-less link.
fn link_solid(
    lname: &str,
    link: &XmlNode,
    warnings: &mut Vec<String>,
) -> Result<valenx_cad::Solid, UrdfError> {
    let geometry = link
        .children
        .iter()
        .filter(|c| c.tag == "visual" || c.tag == "collision")
        .find_map(|vc| vc.child("geometry"));

    let mkerr = |reason: String| UrdfError::Geometry {
        link: lname.to_string(),
        reason,
    };
    let clamp = |v: f64| v.max(1.0e-4);

    let Some(g) = geometry else {
        warnings.push(format!(
            "link `{lname}` has no geometry; using a placeholder box"
        ));
        return valenx_cad::box_solid(0.01, 0.01, 0.01).map_err(|e| mkerr(e.to_string()));
    };
    if let Some(b) = g.child("box") {
        let s = b.attr("size").map_or([0.01; 3], parse3);
        valenx_cad::box_solid(clamp(s[0]), clamp(s[1]), clamp(s[2]))
            .map_err(|e| mkerr(e.to_string()))
    } else if let Some(c) = g.child("cylinder") {
        let r = c
            .attr("radius")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        let h = c
            .attr("length")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        valenx_cad::cylinder(clamp(r), clamp(h)).map_err(|e| mkerr(e.to_string()))
    } else if let Some(s) = g.child("sphere") {
        let r = s
            .attr("radius")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.01);
        valenx_cad::sphere(clamp(r)).map_err(|e| mkerr(e.to_string()))
    } else {
        warnings.push(format!(
            "link `{lname}` uses mesh/unknown geometry; using a placeholder box"
        ));
        valenx_cad::box_solid(0.01, 0.01, 0.01).map_err(|e| mkerr(e.to_string()))
    }
}

/// A rigid transform from a URDF `xyz` translation + `rpy` (fixed-axis XYZ)
/// rotation.
fn rigid(xyz: [f64; 3], rpy: [f64; 3]) -> PartTransform {
    PartTransform {
        translation: Vector3::new(xyz[0], xyz[1], xyz[2]),
        orientation: UnitQuaternion::from_euler_angles(rpy[0], rpy[1], rpy[2]),
    }
}

/// Compose two rigid transforms: `parent ∘ local` (apply `local` in `parent`'s
/// frame).
fn compose(parent: &PartTransform, local: &PartTransform) -> PartTransform {
    PartTransform {
        orientation: parent.orientation * local.orientation,
        translation: parent.orientation * local.translation + parent.translation,
    }
}

/// Parse a whitespace-separated triple like `"0.04 0 0.004"` into `[f64; 3]`,
/// filling missing components with 0.
fn parse3(s: &str) -> [f64; 3] {
    let mut out = [0.0; 3];
    for (slot, tok) in out.iter_mut().zip(s.split_whitespace()) {
        if let Ok(v) = tok.parse::<f64>() {
            *slot = v;
        }
    }
    out
}

// ── A minimal in-crate XML reader (no external dependency) ─────────────────

/// A parsed XML element: tag, attributes, and child elements. Text content is
/// not retained (URDF carries its data in attributes).
struct XmlNode {
    tag: String,
    attrs: Vec<(String, String)>,
    children: Vec<XmlNode>,
}

impl XmlNode {
    /// First attribute value for `key`.
    fn attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// First direct child element with tag `tag`.
    fn child(&self, tag: &str) -> Option<&XmlNode> {
        self.children.iter().find(|c| c.tag == tag)
    }
}

/// Parse an XML document into its root element.
fn parse_xml(input: &str) -> Result<XmlNode, UrdfError> {
    let mut p = XmlParser {
        c: input.chars().collect(),
        i: 0,
    };
    p.skip_misc();
    let node = p.parse_element()?;
    Ok(node)
}

struct XmlParser {
    c: Vec<char>,
    i: usize,
}

impl XmlParser {
    fn peek(&self) -> Option<char> {
        self.c.get(self.i).copied()
    }

    fn starts(&self, pat: &str) -> bool {
        let pc: Vec<char> = pat.chars().collect();
        self.i + pc.len() <= self.c.len() && self.c[self.i..self.i + pc.len()] == pc[..]
    }

    fn skip_ws(&mut self) {
        while self.peek().is_some_and(|ch| ch.is_whitespace()) {
            self.i += 1;
        }
    }

    /// Advance past `pat`; if absent, advance to end of input.
    fn skip_until(&mut self, pat: &str) {
        let pn = pat.chars().count();
        while self.i < self.c.len() && !self.starts(pat) {
            self.i += 1;
        }
        if self.starts(pat) {
            self.i += pn;
        }
    }

    /// Skip whitespace, XML declarations, comments and doctype.
    fn skip_misc(&mut self) {
        loop {
            self.skip_ws();
            if self.starts("<?") {
                self.skip_until("?>");
            } else if self.starts("<!--") {
                self.skip_until("-->");
            } else if self.starts("<!") {
                self.skip_until(">");
            } else {
                break;
            }
        }
    }

    fn take_while(&mut self, pred: impl Fn(char) -> bool) -> String {
        let mut s = String::new();
        while let Some(ch) = self.peek() {
            if pred(ch) {
                s.push(ch);
                self.i += 1;
            } else {
                break;
            }
        }
        s
    }

    fn parse_attrs(&mut self) -> Vec<(String, String)> {
        let mut attrs = Vec::new();
        loop {
            self.skip_ws();
            match self.peek() {
                Some('>') | Some('/') | None => break,
                Some(_) => {
                    let name = self.take_while(|ch| {
                        ch != '=' && !ch.is_whitespace() && ch != '>' && ch != '/'
                    });
                    self.skip_ws();
                    if self.peek() == Some('=') {
                        self.i += 1;
                        self.skip_ws();
                        let q = self.peek();
                        if q == Some('"') || q == Some('\'') {
                            let quote = q.unwrap();
                            self.i += 1;
                            let val = self.take_while(|ch| ch != quote);
                            if self.peek() == Some(quote) {
                                self.i += 1;
                            }
                            attrs.push((name, val));
                        } else {
                            let val =
                                self.take_while(|ch| !ch.is_whitespace() && ch != '>' && ch != '/');
                            attrs.push((name, val));
                        }
                    } else {
                        attrs.push((name, String::new()));
                    }
                }
            }
        }
        attrs
    }

    fn parse_element(&mut self) -> Result<XmlNode, UrdfError> {
        if self.peek() != Some('<') {
            return Err(UrdfError::Xml("expected `<`".into()));
        }
        self.i += 1;
        let tag = self.take_while(|ch| !ch.is_whitespace() && ch != '>' && ch != '/');
        if tag.is_empty() {
            return Err(UrdfError::Xml("empty tag name".into()));
        }
        let attrs = self.parse_attrs();
        self.skip_ws();
        let mut node = XmlNode {
            tag,
            attrs,
            children: Vec::new(),
        };
        if self.starts("/>") {
            self.i += 2;
            return Ok(node);
        }
        if self.peek() != Some('>') {
            return Err(UrdfError::Xml(format!("malformed open tag `{}`", node.tag)));
        }
        self.i += 1;

        loop {
            self.skip_ws();
            if self.starts("<!--") {
                self.skip_until("-->");
            } else if self.starts("</") {
                self.i += 2;
                let _close = self.take_while(|ch| ch != '>');
                if self.peek() == Some('>') {
                    self.i += 1;
                }
                break;
            } else if self.starts("<") {
                node.children.push(self.parse_element()?);
            } else if self.peek().is_none() {
                return Err(UrdfError::Xml(format!("unclosed element `{}`", node.tag)));
            } else {
                // Text content — skip to the next tag.
                self.take_while(|ch| ch != '<');
            }
        }
        Ok(node)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a 16-DOF hand URDF (4 fingers × 4 revolute joints) in code, so the
    /// test owns its fixture without an external file.
    fn hand_urdf_16dof() -> String {
        let mut s = String::from("<?xml version=\"1.0\"?>\n<robot name=\"test_hand\">\n");
        s.push_str(
            "<link name=\"palm\"><visual><geometry><box size=\"0.08 0.09 0.02\"/>\
             </geometry></visual></link>\n",
        );
        let segs = [
            ("knuckle", 0.012),
            ("proximal", 0.045),
            ("medial", 0.030),
            ("distal", 0.026),
        ];
        for f in ["index", "middle", "ring", "thumb"] {
            for (seg, len) in segs {
                s.push_str(&format!(
                    "<link name=\"{f}_{seg}\"><visual><geometry>\
                     <box size=\"{len} 0.015 0.015\"/></geometry></visual></link>\n"
                ));
            }
            s.push_str(&format!(
                "<joint name=\"{f}_j0\" type=\"revolute\"><parent link=\"palm\"/>\
                 <child link=\"{f}_knuckle\"/><origin xyz=\"0.04 0 0\"/><axis xyz=\"0 0 1\"/></joint>\n"
            ));
            s.push_str(&format!(
                "<joint name=\"{f}_j1\" type=\"revolute\"><parent link=\"{f}_knuckle\"/>\
                 <child link=\"{f}_proximal\"/><origin xyz=\"0.012 0 0\"/><axis xyz=\"0 1 0\"/></joint>\n"
            ));
            s.push_str(&format!(
                "<joint name=\"{f}_j2\" type=\"revolute\"><parent link=\"{f}_proximal\"/>\
                 <child link=\"{f}_medial\"/><origin xyz=\"0.045 0 0\"/><axis xyz=\"0 1 0\"/></joint>\n"
            ));
            s.push_str(&format!(
                "<joint name=\"{f}_j3\" type=\"revolute\"><parent link=\"{f}_medial\"/>\
                 <child link=\"{f}_distal\"/><origin xyz=\"0.030 0 0\"/><axis xyz=\"0 1 0\"/></joint>\n"
            ));
        }
        s.push_str("</robot>\n");
        s
    }

    #[test]
    fn imports_16dof_hand_as_a_native_assembly() {
        let robot = import_urdf(&hand_urdf_16dof()).expect("import");
        assert_eq!(robot.name, "test_hand");
        assert_eq!(robot.assembly.parts.len(), 17, "palm + 16 finger links");
        assert_eq!(robot.assembly.joints.len(), 16);
        let revolute = robot
            .assembly
            .joints
            .iter()
            .filter(|j| matches!(j.kind, JointKind::Revolute { .. }))
            .count();
        assert_eq!(revolute, 16, "16 revolute DOF");
        // The palm is the tree root → grounded/fixed.
        let palm = robot.link_ids["palm"];
        assert!(robot.assembly.get_part(palm).unwrap().fixed);
        // Forward kinematics advanced the fingertip well down +x from the palm.
        let tip = robot.link_ids["index_distal"];
        let x = robot
            .assembly
            .get_part(tip)
            .unwrap()
            .transform
            .translation
            .x;
        assert!(x > 0.05, "fingertip should be advanced along +x, got {x}");
    }

    #[test]
    fn builds_a_renderable_mesh_in_house() {
        let robot = import_urdf(&hand_urdf_16dof()).expect("import");
        let mesh = assembly_to_mesh(&robot.assembly, 0.01);
        assert!(!mesh.nodes.is_empty(), "the merged robot mesh has vertices");
        let conn = &mesh.element_blocks[0].connectivity;
        assert!(!conn.is_empty() && conn.len() % 3 == 0, "triangulated");
    }

    #[test]
    fn maps_fixed_prismatic_and_geometry_kinds() {
        let xml = "<robot name=\"r\">\
            <link name=\"a\"><visual><geometry><box size=\"1 1 1\"/></geometry></visual></link>\
            <link name=\"b\"><collision><geometry><cylinder radius=\"0.1\" length=\"0.5\"/>\
            </geometry></collision></link>\
            <joint name=\"slide\" type=\"prismatic\"><parent link=\"a\"/><child link=\"b\"/>\
            <axis xyz=\"1 0 0\"/></joint></robot>";
        let r = import_urdf(xml).expect("import");
        assert_eq!(r.assembly.parts.len(), 2);
        assert!(matches!(
            r.assembly.joints[0].kind,
            JointKind::Prismatic { .. }
        ));
    }

    #[test]
    fn rejects_unknown_link_and_non_robot_root() {
        let dangling = "<robot name=\"r\">\
            <link name=\"a\"><visual><geometry><box size=\"1 1 1\"/></geometry></visual></link>\
            <joint name=\"j\" type=\"fixed\"><parent link=\"a\"/><child link=\"ghost\"/></joint></robot>";
        assert!(matches!(
            import_urdf(dangling),
            Err(UrdfError::UnknownLink { .. })
        ));
        assert!(matches!(
            import_urdf("<notrobot/>"),
            Err(UrdfError::NoRobot)
        ));
    }

    #[test]
    fn xml_reader_handles_comments_and_self_closing_tags() {
        let xml = "<!-- a hand --><robot name=\"r\">\
            <link name=\"a\"><visual><geometry><sphere radius=\"0.02\"/></geometry></visual></link>\
            </robot>";
        let r = import_urdf(xml).expect("import");
        assert_eq!(r.assembly.parts.len(), 1);
        assert_eq!(r.name, "r");
    }
}
