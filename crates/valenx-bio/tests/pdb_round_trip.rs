use std::path::PathBuf;
use valenx_bio::format::pdb;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("biology")
        .join("1ubq-tiny.pdb")
}

#[test]
fn read_tiny_ubq_extract() {
    let text = std::fs::read_to_string(fixture()).unwrap();
    let s = pdb::read("1ubq-tiny", &text).unwrap();
    assert_eq!(s.id, "1ubq-tiny");
    assert_eq!(s.chains.len(), 1);
    assert_eq!(s.chains[0].id, 'A');
    // 2 residues from the snippet (MET, GLN).
    assert_eq!(s.chains[0].residues.len(), 2);
    assert_eq!(s.atom_count(), 9);
}

#[test]
fn accepts_66_char_atom_line_without_element_symbol() {
    // Phenix-style 66-char ATOM record: no element symbol at cols
    // 77-78. Pre-fix this would reject as "too short (66 < 78)".
    // Columns:                  111111111122222222223333333333444444444455555555556666666666
    //                  123456789012345678901234567890123456789012345678901234567890123456
    let line = "ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00 49.05";
    assert_eq!(line.len(), 66);
    let s = pdb::read("p", line).expect("66-char line should parse");
    assert_eq!(s.atom_count(), 1);
    let atom = &s.chains[0].residues[0].atoms[0];
    assert_eq!(atom.name, "N");
    assert!(atom.element.is_empty(), "no element when line < 78 wide");
    assert!((atom.position.x - 27.340).abs() < 1e-6);
}

#[test]
fn rejects_too_short_atom_line() {
    // 50-char line — below the 66-char minimum (no B-factor / temp).
    let line = "ATOM      1  N   MET A   1      27.340  24.430";
    assert!(pdb::read("p", line).is_err());
}

#[test]
fn dedups_alt_loc_keeping_first_conformer() {
    // Two ATOM records, same atom name (" CA "), different altLoc
    // (A then B). Pre-fix both got appended, inflating the count.
    // Post-fix we keep only the first.
    let pdb = "\
ATOM      1  CA AMET A   1      27.340  24.430   2.614  1.00 49.05           C
ATOM      2  CA BMET A   1      27.500  24.600   2.700  1.00 50.00           C
ATOM      3  N   MET A   1      28.000  25.000   3.000  1.00 50.00           N
";
    let s = pdb::read("alt", pdb).unwrap();
    // CA should appear once (altLoc 'A' kept), plus N — total 2.
    assert_eq!(s.atom_count(), 2);
    let ca = s.chains[0].residues[0]
        .atoms
        .iter()
        .find(|a| a.name == "CA")
        .expect("CA present");
    // First conformer's coords (27.340).
    assert!((ca.position.x - 27.340).abs() < 1e-6, "got: {ca:?}");
}
