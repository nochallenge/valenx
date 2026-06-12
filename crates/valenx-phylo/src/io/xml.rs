//! PhyloXML and NeXML writers.
//!
//! Both are XML serialisations of a phylogenetic tree. They are
//! *writers only* — Valenx reads trees from Newick / NEXUS, and the XML
//! formats are emitted for interchange with web tools (phyloXML
//! viewers, the NeXML ecosystem).
//!
//! - **PhyloXML** (<http://www.phyloxml.org>) nests `<clade>` elements,
//!   each carrying an optional `<name>` and `<branch_length>`.
//! - **NeXML** (<http://www.nexml.org>) is flatter: a `<otus>` block
//!   declares the taxa, then a `<tree>` lists every `<node>` and
//!   `<edge>` by id.
//!
//! The output is hand-built with explicit XML escaping; no XML
//! dependency is pulled in for two small writers.

use crate::tree::{NodeId, Tree};

/// Escapes the five XML special characters in element text / attribute
/// values.
fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

/// Serialises a [`Tree`] to a PhyloXML document.
///
/// The document declares one `<phylogeny>` whose `rooted` attribute
/// mirrors [`Tree::rooted`]; `<clade>` elements nest to encode the
/// topology.
pub fn write_phyloxml(tree: &Tree) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(
        "<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
xmlns=\"http://www.phyloxml.org\">\n",
    );
    out.push_str(&format!(
        "  <phylogeny rooted=\"{}\">\n",
        if tree.rooted { "true" } else { "false" }
    ));
    write_clade(tree, tree.root(), 2, &mut out);
    out.push_str("  </phylogeny>\n");
    out.push_str("</phyloxml>\n");
    out
}

fn write_clade(tree: &Tree, id: NodeId, depth: usize, out: &mut String) {
    let pad = "  ".repeat(depth);
    let node = tree.node(id);
    out.push_str(&format!("{pad}<clade>\n"));
    if let Some(name) = &node.label {
        out.push_str(&format!("{pad}  <name>{}</name>\n", xml_escape(name)));
    }
    if let Some(bl) = node.branch_length {
        out.push_str(&format!("{pad}  <branch_length>{bl}</branch_length>\n"));
    }
    for &c in &node.children {
        write_clade(tree, c, depth + 1, out);
    }
    out.push_str(&format!("{pad}</clade>\n"));
}

/// Serialises a [`Tree`] to a NeXML document.
///
/// All leaves become `<otu>` elements inside one `<otus>` block; the
/// `<tree>` then lists every node and edge by stable `n<id>` / `e<id>`
/// identifiers.
pub fn write_nexml(tree: &Tree) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(
        "<nex:nexml version=\"0.9\" \
xmlns:nex=\"http://www.nexml.org/2009\" \
xmlns=\"http://www.nexml.org/2009\" \
xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n",
    );

    // OTUs — one per leaf.
    out.push_str("  <otus id=\"taxa1\">\n");
    for &leaf in &tree.leaves() {
        let label = tree
            .node(leaf)
            .label
            .clone()
            .unwrap_or_else(|| format!("leaf{leaf}"));
        out.push_str(&format!(
            "    <otu id=\"otu{leaf}\" label=\"{}\"/>\n",
            xml_escape(&label)
        ));
    }
    out.push_str("  </otus>\n");

    // The tree itself.
    out.push_str("  <trees id=\"trees1\" otus=\"taxa1\">\n");
    out.push_str("    <tree id=\"tree1\" xsi:type=\"nex:FloatTree\">\n");
    for id in tree.preorder() {
        let node = tree.node(id);
        let is_root = id == tree.root();
        let root_attr = if is_root { " root=\"true\"" } else { "" };
        let otu_attr = if node.is_leaf() {
            format!(" otu=\"otu{id}\"")
        } else {
            String::new()
        };
        let label_attr = node
            .label
            .as_ref()
            .map(|l| format!(" label=\"{}\"", xml_escape(l)))
            .unwrap_or_default();
        out.push_str(&format!(
            "      <node id=\"n{id}\"{otu_attr}{label_attr}{root_attr}/>\n"
        ));
    }
    for id in tree.preorder() {
        let node = tree.node(id);
        if let Some(parent) = node.parent {
            let len = node.branch_length.unwrap_or(0.0);
            out.push_str(&format!(
                "      <edge id=\"e{id}\" source=\"n{parent}\" target=\"n{id}\" length=\"{len}\"/>\n"
            ));
        }
    }
    out.push_str("    </tree>\n");
    out.push_str("  </trees>\n");
    out.push_str("</nex:nexml>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn phyloxml_has_expected_structure() {
        let t = read_newick("((A:0.1,B:0.2):0.3,C:0.5);").unwrap();
        let xml = write_phyloxml(&t);
        assert!(xml.contains("<phyloxml"));
        assert!(xml.contains("<phylogeny rooted=\"true\">"));
        assert!(xml.contains("<name>A</name>"));
        assert!(xml.contains("<branch_length>0.1</branch_length>"));
        // Three leaves => at least three <name> elements.
        assert_eq!(xml.matches("<name>").count(), 3);
    }

    #[test]
    fn nexml_declares_one_otu_per_leaf() {
        let t = read_newick("((A,B),(C,D));").unwrap();
        let xml = write_nexml(&t);
        assert!(xml.contains("<nex:nexml"));
        assert_eq!(xml.matches("<otu ").count(), 4);
        // One node per arena entry, one edge per non-root node.
        assert_eq!(xml.matches("<node ").count(), t.node_count());
        assert_eq!(xml.matches("<edge ").count(), t.node_count() - 1);
        assert!(xml.contains("root=\"true\""));
    }

    #[test]
    fn xml_escaping_handles_specials() {
        let t = read_newick("('A & B':1,C:1);").unwrap();
        let xml = write_phyloxml(&t);
        assert!(xml.contains("A &amp; B"), "got: {xml}");
        assert!(!xml.contains("A & B"));
    }
}
