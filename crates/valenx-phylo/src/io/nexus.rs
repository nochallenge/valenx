//! NEXUS file format — reader and writer.
//!
//! NEXUS is the block-structured format used by PAUP*, MrBayes and
//! FigTree. A file is `#NEXUS` followed by `BEGIN <name>; … END;`
//! blocks. This module reads and writes the three blocks a
//! phylogenetics pipeline needs:
//!
//! - **`TAXA`** — a `DIMENSIONS NTAX=n;` and a `TAXLABELS a b c;` list.
//! - **`TREES`** — optional `TRANSLATE` table plus `TREE name = …;`
//!   lines whose value is a [Newick](super::newick) string.
//! - **`DATA`** (a.k.a. `CHARACTERS`) — a `MATRIX` of taxon → sequence
//!   rows; read and written minimally (interleaved matrices are *not*
//!   supported — one row per taxon).
//!
//! The reader is forgiving: it is case-insensitive on keywords, treats
//! `[...]` as comments, and tolerates blocks it does not recognise by
//! skipping to the matching `END;`.

use crate::error::{PhyloError, Result};
use crate::io::newick::{read_newick, write_newick};
use crate::tree::Tree;

/// The parsed contents of a NEXUS file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NexusFile {
    /// Taxon labels from the `TAXA` block, in declared order.
    pub taxa: Vec<String>,
    /// Named trees from the `TREES` block (name, tree) pairs.
    pub trees: Vec<(String, Tree)>,
    /// Aligned sequence rows from the `DATA` / `CHARACTERS` block:
    /// `(taxon, sequence)` pairs.
    pub data: Vec<(String, String)>,
}

impl NexusFile {
    /// Convenience: the first tree in the file, if any.
    pub fn first_tree(&self) -> Option<&Tree> {
        self.trees.first().map(|(_, t)| t)
    }
}

/// Removes `[...]` NEXUS comment blocks (which may nest).
fn strip_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth = 0u32;
    for c in s.chars() {
        match c {
            '[' => depth += 1,
            ']' if depth > 0 => depth -= 1,
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

/// Parses a NEXUS document.
///
/// # Errors
/// [`PhyloError::Parse`] if the `#NEXUS` header is missing, a block is
/// not terminated, or an embedded Newick tree fails to parse.
pub fn read_nexus(input: &str) -> Result<NexusFile> {
    let cleaned = strip_comments(input);
    let upper = cleaned.to_uppercase();
    if !upper.trim_start().starts_with("#NEXUS") {
        return Err(PhyloError::parse("nexus", "missing #NEXUS header"));
    }
    let mut file = NexusFile::default();
    // Split into `BEGIN … END;` blocks. All keyword scanning is done
    // case-insensitively *on the original `cleaned` string itself*
    // (never on a separate `to_lowercase()` copy), so every byte offset
    // is valid in `cleaned` — a lowercased copy can change UTF-8 byte
    // lengths (e.g. `İ` U+0130 → `i̇`), and offsets taken from such a copy
    // then applied back to the original overshoot or split a char and
    // panic. `cleaned` stays the single source of truth, so block bodies
    // also retain their original case for labels/sequences.
    let mut idx = 0usize;
    while let Some(rel) = find_ci(&cleaned[idx..], "begin ") {
        let block_start = idx + rel;
        let after_begin = block_start + "begin ".len();
        let semi = cleaned[after_begin..]
            .find(';')
            .ok_or_else(|| PhyloError::parse("nexus", "BEGIN with no ';'"))?;
        let block_name = cleaned[after_begin..after_begin + semi].trim().to_lowercase();
        let body_start = after_begin + semi + 1;
        let end_rel = find_ci(&cleaned[body_start..], "end;")
            .or_else(|| find_ci(&cleaned[body_start..], "endblock;"))
            .ok_or_else(|| {
                PhyloError::parse("nexus", format!("block `{block_name}` not terminated"))
            })?;
        let body = &cleaned[body_start..body_start + end_rel];
        match block_name.as_str() {
            "taxa" => parse_taxa_block(body, &mut file),
            "trees" => parse_trees_block(body, &mut file)?,
            "data" | "characters" => parse_data_block(body, &mut file),
            _ => { /* unknown block — skipped intentionally */ }
        }
        idx = body_start + end_rel;
    }
    Ok(file)
}

/// Finds the first case-insensitive occurrence of the ASCII `needle` in
/// `haystack`, returning a byte offset that is valid in `haystack`.
///
/// Unlike searching a `to_lowercase()` copy (whose byte offsets can drift
/// from the original when a char's lowercase has a different UTF-8
/// length), this scans `haystack` directly, so the returned index always
/// lands on a char boundary of `haystack`. `needle` must be ASCII and
/// non-empty.
fn find_ci(haystack: &str, needle: &str) -> Option<usize> {
    debug_assert!(needle.is_ascii() && !needle.is_empty());
    let hay = haystack.as_bytes();
    let need = needle.as_bytes();
    if need.len() > hay.len() {
        return None;
    }
    // ASCII needle ⇒ every match starts and ends on a char boundary
    // (an ASCII byte can never be a UTF-8 continuation byte).
    (0..=hay.len() - need.len()).find(|&i| hay[i..i + need.len()].eq_ignore_ascii_case(need))
}

/// Splits a block body into `;`-terminated commands.
fn commands(body: &str) -> Vec<String> {
    body.split(';')
        .map(|c| c.trim().to_string())
        .filter(|c| !c.is_empty())
        .collect()
}

/// If `s` begins with the ASCII `keyword` (case-insensitively), returns
/// the remainder of `s` *in its original case*.
///
/// Matching is done directly on `s`'s own bytes — never on a separate
/// `to_lowercase()` copy — so the returned slice boundary is always a
/// valid char boundary in `s`. (A `to_lowercase()` copy can change UTF-8
/// byte lengths, e.g. `İ` U+0130 → `i̇`; offsets taken from such a copy
/// and applied back to the original overshoot or split a char and
/// panic.) `keyword` must be ASCII; if `s` matches it, the leading
/// `keyword.len()` bytes of `s` are therefore all ASCII and form a valid
/// boundary.
fn strip_keyword_ci<'a>(s: &'a str, keyword: &str) -> Option<&'a str> {
    debug_assert!(keyword.is_ascii());
    let n = keyword.len();
    let head = s.as_bytes().get(..n)?;
    if head.eq_ignore_ascii_case(keyword.as_bytes()) {
        // `head` is ASCII (it case-insensitively equals an ASCII
        // keyword), so byte index `n` is a char boundary in `s`.
        Some(&s[n..])
    } else {
        None
    }
}

fn parse_taxa_block(body: &str, file: &mut NexusFile) {
    for cmd in commands(body) {
        if let Some(labels) = strip_keyword_ci(&cmd, "taxlabels") {
            file.taxa = split_tokens(labels);
        }
    }
}

fn parse_trees_block(body: &str, file: &mut NexusFile) -> Result<()> {
    // An optional TRANSLATE table maps short ids to taxon labels.
    let mut translate: Vec<(String, String)> = Vec::new();
    for cmd in commands(body) {
        let lower = cmd.to_lowercase();
        if let Some(table) = strip_keyword_ci(&cmd, "translate") {
            for pair in table.split(',') {
                let toks = split_tokens(pair);
                if toks.len() == 2 {
                    translate.push((toks[0].clone(), toks[1].clone()));
                }
            }
        } else if lower.starts_with("tree ") || lower.starts_with("utree ") {
            // `tree NAME = (newick);`  — split on the first '='.
            let eq = cmd
                .find('=')
                .ok_or_else(|| PhyloError::parse("nexus", "TREE command without '='"))?;
            let name_part = cmd[..eq].trim();
            let name = name_part
                .split_whitespace()
                .nth(1)
                .unwrap_or("tree")
                .trim_start_matches('*')
                .to_string();
            let newick = cmd[eq + 1..].trim();
            let tree = read_newick(newick)?;
            let tree = apply_translation(tree, &translate);
            file.trees.push((name, tree));
        }
    }
    Ok(())
}

/// Rewrites leaf labels through a TRANSLATE table.
fn apply_translation(mut tree: Tree, translate: &[(String, String)]) -> Tree {
    if translate.is_empty() {
        return tree;
    }
    for i in 0..tree.node_count() {
        let relabel = {
            let node = tree.node(i);
            node.label.as_ref().and_then(|lbl| {
                translate
                    .iter()
                    .find(|(id, _)| id == lbl)
                    .map(|(_, name)| name.clone())
            })
        };
        if let Some(name) = relabel {
            tree.node_mut(i).label = Some(name);
        }
    }
    tree
}

fn parse_data_block(body: &str, file: &mut NexusFile) {
    // Find the MATRIX command and read whitespace-separated
    // taxon/sequence pairs from it.
    for cmd in commands(body) {
        if let Some(matrix) = strip_keyword_ci(&cmd, "matrix") {
            for line in matrix.lines() {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let mut it = line.split_whitespace();
                if let (Some(taxon), Some(seq)) = (it.next(), it.next()) {
                    file.data.push((
                        taxon.trim_matches('\'').to_string(),
                        seq.to_string(),
                    ));
                }
            }
        }
    }
}

/// Splits a token list honouring single-quoted multi-word labels.
fn split_tokens(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    for c in s.chars() {
        match c {
            '\'' => quoted = !quoted,
            c if c.is_whitespace() && !quoted => {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
            }
            c => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// Serialises a [`NexusFile`] to a NEXUS document string.
///
/// Emits a `TAXA` block when `taxa` is non-empty, a `DATA` block when
/// `data` is non-empty, and a `TREES` block when `trees` is non-empty.
pub fn write_nexus(file: &NexusFile) -> String {
    let mut out = String::from("#NEXUS\n\n");

    if !file.taxa.is_empty() {
        out.push_str("BEGIN TAXA;\n");
        out.push_str(&format!("    DIMENSIONS NTAX={};\n", file.taxa.len()));
        out.push_str("    TAXLABELS");
        for t in &file.taxa {
            out.push(' ');
            out.push_str(&quote_if_needed(t));
        }
        out.push_str(";\nEND;\n\n");
    }

    if !file.data.is_empty() {
        let nchar = file.data.first().map(|(_, s)| s.len()).unwrap_or(0);
        out.push_str("BEGIN DATA;\n");
        out.push_str(&format!(
            "    DIMENSIONS NTAX={} NCHAR={};\n",
            file.data.len(),
            nchar
        ));
        out.push_str("    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n");
        out.push_str("    MATRIX\n");
        for (taxon, seq) in &file.data {
            out.push_str(&format!("        {}  {}\n", quote_if_needed(taxon), seq));
        }
        out.push_str("    ;\nEND;\n\n");
    }

    if !file.trees.is_empty() {
        out.push_str("BEGIN TREES;\n");
        for (name, tree) in &file.trees {
            out.push_str(&format!("    TREE {} = {}\n", name, write_newick(tree)));
        }
        out.push_str("END;\n");
    }
    out
}

/// Single-quotes a token if it contains whitespace or NEXUS specials.
fn quote_if_needed(s: &str) -> String {
    if s.is_empty()
        || s.chars()
            .any(|c| c.is_whitespace() || matches!(c, '(' | ')' | ',' | ';' | ':' | '\'' | '['))
    {
        format!("'{}'", s.replace('\'', "''"))
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=3;
    TAXLABELS A B C;
END;
BEGIN TREES;
    TREE t1 = ((A:0.1,B:0.2):0.3,C:0.5);
END;
"#;

    #[test]
    fn reads_taxa_and_trees() {
        let f = read_nexus(SAMPLE).unwrap();
        assert_eq!(f.taxa, vec!["A", "B", "C"]);
        assert_eq!(f.trees.len(), 1);
        assert_eq!(f.trees[0].0, "t1");
        assert_eq!(f.first_tree().unwrap().leaf_count(), 3);
    }

    #[test]
    fn round_trip() {
        let f = read_nexus(SAMPLE).unwrap();
        let text = write_nexus(&f);
        let f2 = read_nexus(&text).unwrap();
        assert_eq!(f.taxa, f2.taxa);
        assert_eq!(f2.trees.len(), 1);
        assert_eq!(f2.first_tree().unwrap().leaf_labels(), vec!["A", "B", "C"]);
    }

    #[test]
    fn translate_table_is_applied() {
        let src = r#"#NEXUS
BEGIN TREES;
    TRANSLATE 1 Alpha, 2 Beta, 3 Gamma;
    TREE x = ((1:1,2:1):1,3:1);
END;
"#;
        let f = read_nexus(src).unwrap();
        let t = f.first_tree().unwrap();
        assert!(t.find("Alpha").is_some());
        assert!(t.find("Gamma").is_some());
        assert!(t.find("1").is_none());
    }

    #[test]
    fn data_block_round_trips() {
        let src = r#"#NEXUS
BEGIN DATA;
    DIMENSIONS NTAX=2 NCHAR=4;
    FORMAT DATATYPE=DNA;
    MATRIX
        A  ACGT
        B  ACGA
    ;
END;
"#;
        let f = read_nexus(src).unwrap();
        assert_eq!(f.data.len(), 2);
        assert_eq!(f.data[0], ("A".to_string(), "ACGT".to_string()));
        let text = write_nexus(&f);
        let f2 = read_nexus(&text).unwrap();
        assert_eq!(f.data, f2.data);
    }

    #[test]
    fn rejects_missing_header() {
        assert!(read_nexus("BEGIN TAXA; END;").is_err());
    }

    #[test]
    fn comments_are_stripped() {
        let src = "#NEXUS\n[ a comment ]\nBEGIN TAXA;\n TAXLABELS A B;\nEND;\n";
        let f = read_nexus(src).unwrap();
        assert_eq!(f.taxa, vec!["A", "B"]);
    }

    #[test]
    fn unicode_special_case_char_does_not_panic() {
        // `İ` (U+0130, LATIN CAPITAL LETTER I WITH DOT ABOVE) lowercases
        // to `i̇` (U+0069 U+0307) — TWO bytes become THREE. The parser
        // computed byte offsets on a `to_lowercase()` copy and then
        // indexed the ORIGINAL string with them; once a case char of a
        // different UTF-8 length appears in a block name or a command
        // body, those offsets overshoot or land mid-char and the string
        // slice panics. Every block must parse (or error) gracefully —
        // never abort the process.
        //
        // A `İ` in: (1) a TAXLABELS body, (2) a TRANSLATE body,
        // (3) a MATRIX body, (4) a block name.
        let srcs = [
            "#NEXUS\nBEGIN TAXA;\n TAXLABELS A\u{130} B;\nEND;\n",
            "#NEXUS\nBEGIN TREES;\n TRANSLATE 1 A\u{130}, 2 B;\n TREE t = (1:1,2:1);\nEND;\n",
            "#NEXUS\nBEGIN DATA;\n MATRIX\n  A\u{130} ACGT\n  B ACGA\n ;\nEND;\n",
            "#NEXUS\nBEGIN F\u{130}OO;\n SOMETHING x;\nEND;\n",
        ];
        for src in srcs {
            // The contract is "does not panic"; a clean Ok or Err is fine.
            let _ = read_nexus(src);
        }

        // And a concrete positive: a `İ` taxon label is read with its
        // case preserved (the original-case text is not corrupted).
        let f =
            read_nexus("#NEXUS\nBEGIN TAXA;\n TAXLABELS Foo\u{130}bar Baz;\nEND;\n").unwrap();
        assert_eq!(f.taxa, vec!["Foo\u{130}bar".to_string(), "Baz".to_string()]);
    }
}
