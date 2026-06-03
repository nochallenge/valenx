//! Tree rendering — ASCII cladograms and SVG phylograms.
//!
//! Two layouts:
//!
//! - [`render_ascii`] draws a text cladogram (FigTree's text view,
//!   `Bio.Phylo.draw_ascii`-style). Internal structure is drawn with
//!   the box-drawing characters `|`, `_` and `-`; leaves are written at
//!   the right margin.
//! - [`render_svg`] draws a rectangular phylogram: horizontal extent is
//!   proportional to branch length, vertical position spreads the
//!   leaves evenly.
//!
//! Both share one internal layout pass that assigns every node an
//! `(x, y)` coordinate; the ASCII renderer then snaps to a character
//! grid while the SVG renderer uses the floating-point coordinates
//! directly.

use crate::tree::Tree;

/// Per-node layout coordinates produced by [`layout`].
#[derive(Debug, Clone, Copy, PartialEq)]
struct Coord {
    /// Horizontal position — cumulative branch length from the root.
    x: f64,
    /// Vertical position — leaf index for tips, child mean for internals.
    y: f64,
}

/// Computes `(x, y)` coordinates for every node.
///
/// `x` is the root-to-node distance (cumulative branch length, missing
/// lengths counting as `unit`). `y` numbers leaves `0, 1, 2, …` top to
/// bottom; an internal node sits at the mean `y` of its children.
fn layout(tree: &Tree, unit: f64) -> Vec<Coord> {
    let n = tree.node_count();
    let mut coords = vec![Coord { x: 0.0, y: 0.0 }; n];

    // x: a preorder pass accumulates branch lengths from the parent.
    for &id in &tree.preorder() {
        if let Some(p) = tree.node(id).parent {
            let bl = tree.node(id).branch_length.unwrap_or(unit);
            coords[id].x = coords[p].x + bl.max(0.0);
        }
    }

    // y: leaves get successive integers; a postorder pass averages
    // each internal node over its children.
    let mut next_leaf = 0.0;
    for &id in &tree.postorder() {
        let node = tree.node(id);
        if node.is_leaf() {
            coords[id].y = next_leaf;
            next_leaf += 1.0;
        } else {
            let sum: f64 = node.children.iter().map(|&c| coords[c].y).sum();
            coords[id].y = sum / node.children.len() as f64;
        }
    }
    coords
}

/// Renders a tree as a text cladogram.
///
/// `width` is the target inner drawing width in characters (the actual
/// string is wider — leaf labels extend past it). A reasonable default
/// is 60.
pub fn render_ascii(tree: &Tree, width: usize) -> String {
    if tree.is_trivial() {
        return tree
            .node(tree.root())
            .label
            .clone()
            .unwrap_or_else(|| "(empty)".into());
    }
    let width = width.max(10);
    let coords = layout(tree, 1.0);
    let max_x = coords.iter().map(|c| c.x).fold(0.0_f64, f64::max).max(1e-9);
    let n_leaves = tree.leaf_count();
    // One text row per leaf (leaves occupy y = 0..n_leaves-1).
    let rows = n_leaves;
    let mut grid: Vec<Vec<char>> = vec![vec![' '; width + 1]; rows];

    // Snap a layout x to a column.
    let col = |x: f64| -> usize {
        ((x / max_x) * width as f64).round() as usize
    };
    // Snap a layout y to a row.
    let row = |y: f64| -> usize { y.round() as usize };

    let mut labels: Vec<(usize, String)> = Vec::new();

    for &id in &tree.preorder() {
        let node = tree.node(id);
        let r = row(coords[id].y).min(rows - 1);
        let c = col(coords[id].x).min(width);
        if node.is_leaf() {
            // Horizontal branch from the parent to this leaf.
            if let Some(p) = node.parent {
                let pc = col(coords[p].x).min(width);
                for x in pc..=c {
                    if grid[r][x] == ' ' {
                        grid[r][x] = '_';
                    }
                }
            }
            let name = node.label.clone().unwrap_or_else(|| format!("n{id}"));
            labels.push((r, name));
        } else {
            // Vertical connector spanning the children's rows.
            let child_rows: Vec<usize> = node
                .children
                .iter()
                .map(|&ch| row(coords[ch].y).min(rows - 1))
                .collect();
            let (top, bot) = (
                *child_rows.iter().min().unwrap(),
                *child_rows.iter().max().unwrap(),
            );
            for line in grid.iter_mut().take(bot + 1).skip(top) {
                if line[c] == ' ' || line[c] == '_' {
                    line[c] = '|';
                }
            }
            // Horizontal stub from the parent into this node.
            if let Some(p) = node.parent {
                let pc = col(coords[p].x).min(width);
                for x in pc..c {
                    if grid[r][x] == ' ' {
                        grid[r][x] = '_';
                    }
                }
            }
        }
    }

    // Assemble: each grid row, then its leaf label if any.
    let mut out = String::new();
    for (r, line) in grid.iter().enumerate() {
        let text: String = line.iter().collect();
        out.push_str(text.trim_end());
        if let Some((_, name)) = labels.iter().find(|(lr, _)| *lr == r) {
            out.push(' ');
            out.push_str(name);
        }
        out.push('\n');
    }
    out
}

/// Renders a tree as a rectangular-phylogram SVG document.
///
/// `width` and `height` are the SVG canvas dimensions in pixels. Branch
/// lengths drive the horizontal extent; the result is a self-contained
/// `<svg>` string suitable for embedding or writing to a `.svg` file.
pub fn render_svg(tree: &Tree, width: u32, height: u32) -> String {
    let margin = 10.0;
    let label_space = 110.0;
    let w = width as f64;
    let h = height as f64;
    let plot_w = (w - 2.0 * margin - label_space).max(50.0);
    let plot_h = (h - 2.0 * margin).max(20.0);

    if tree.is_trivial() {
        return format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\">\
<text x=\"{margin}\" y=\"20\">{}</text></svg>",
            tree.node(tree.root()).label.as_deref().unwrap_or("(empty)")
        );
    }

    let coords = layout(tree, 1.0);
    let max_x = coords.iter().map(|c| c.x).fold(0.0_f64, f64::max).max(1e-9);
    let n_leaves = (tree.leaf_count().max(2) - 1) as f64;

    let px = |x: f64| margin + (x / max_x) * plot_w;
    let py = |y: f64| margin + (y / n_leaves) * plot_h;

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" \
viewBox=\"0 0 {width} {height}\">\n"
    ));
    svg.push_str("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");
    svg.push_str("  <g stroke=\"black\" stroke-width=\"1.5\" fill=\"none\">\n");

    // One horizontal segment per edge plus a vertical segment joining
    // each internal node's children.
    for &id in &tree.preorder() {
        let node = tree.node(id);
        if let Some(p) = node.parent {
            // Horizontal edge: from the parent's x to this node's x,
            // at this node's y.
            svg.push_str(&format!(
                "    <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\"/>\n",
                px(coords[p].x),
                py(coords[id].y),
                px(coords[id].x),
                py(coords[id].y),
            ));
        }
        if node.is_internal() {
            // Vertical connector at this node's x, spanning children.
            let ys: Vec<f64> = node.children.iter().map(|&c| coords[c].y).collect();
            let (top, bot) = (
                ys.iter().cloned().fold(f64::INFINITY, f64::min),
                ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            );
            svg.push_str(&format!(
                "    <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\"/>\n",
                px(coords[id].x),
                py(top),
                px(coords[id].x),
                py(bot),
            ));
        }
    }
    svg.push_str("  </g>\n");

    // Leaf labels at the right margin.
    svg.push_str("  <g font-family=\"sans-serif\" font-size=\"12\" fill=\"black\">\n");
    for &leaf in &tree.leaves() {
        let name = tree
            .node(leaf)
            .label
            .clone()
            .unwrap_or_else(|| format!("n{leaf}"));
        let escaped = name
            .replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;");
        svg.push_str(&format!(
            "    <text x=\"{:.2}\" y=\"{:.2}\" dominant-baseline=\"middle\">{}</text>\n",
            px(coords[leaf].x) + 4.0,
            py(coords[leaf].y),
            escaped,
        ));
    }
    svg.push_str("  </g>\n");
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn ascii_lists_every_leaf() {
        let t = read_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let art = render_ascii(&t, 40);
        for taxon in ["A", "B", "C", "D"] {
            assert!(art.contains(taxon), "missing {taxon} in:\n{art}");
        }
        // One text row per leaf.
        assert_eq!(art.lines().count(), 4);
    }

    #[test]
    fn ascii_handles_trivial_tree() {
        let t = read_newick("A;").unwrap();
        assert_eq!(render_ascii(&t, 40).trim(), "A");
    }

    #[test]
    fn svg_is_well_formed_and_labelled() {
        let t = read_newick("((A:0.2,B:0.3):0.1,C:0.5);").unwrap();
        let svg = render_svg(&t, 400, 300);
        assert!(svg.starts_with("<svg"));
        assert!(svg.trim_end().ends_with("</svg>"));
        assert!(svg.contains(">A<") && svg.contains(">B<") && svg.contains(">C<"));
        // Each non-root edge draws a line; internals add verticals.
        assert!(svg.matches("<line").count() >= t.node_count() - 1);
    }

    #[test]
    fn svg_escapes_label_specials() {
        let t = read_newick("('x<y':1,z:1);").unwrap();
        let svg = render_svg(&t, 300, 200);
        assert!(svg.contains("x&lt;y"));
    }
}
