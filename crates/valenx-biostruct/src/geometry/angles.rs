//! Bond angles, dihedral / torsion angles, and protein backbone
//! (φ/ψ/ω) and sidechain (χ) torsions.
//!
//! All angles are returned in **degrees**. Dihedrals follow the
//! IUPAC convention: the angle between the plane of atoms `(1,2,3)`
//! and the plane of atoms `(2,3,4)`, in `(-180, 180]`.

use crate::structure::{Chain, Residue};
use nalgebra::{Point3, Vector3};

/// Bond angle at the central atom `b` of the triple `a–b–c`, degrees.
/// Returns `None` if any two atoms coincide.
pub fn bond_angle(a: &Point3<f64>, b: &Point3<f64>, c: &Point3<f64>) -> Option<f64> {
    let ba = a - b;
    let bc = c - b;
    let nba = ba.norm();
    let nbc = bc.norm();
    if nba < 1e-9 || nbc < 1e-9 {
        return None;
    }
    let cos = (ba.dot(&bc) / (nba * nbc)).clamp(-1.0, 1.0);
    Some(cos.acos().to_degrees())
}

/// Dihedral (torsion) angle of the four atoms `a–b–c–d`, in degrees,
/// in `(-180, 180]`. Returns `None` for a degenerate geometry.
pub fn dihedral(a: &Point3<f64>, b: &Point3<f64>, c: &Point3<f64>, d: &Point3<f64>) -> Option<f64> {
    let b1 = b - a;
    let b2 = c - b;
    let b3 = d - c;
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let n1n = n1.norm();
    let n2n = n2.norm();
    let b2n = b2.norm();
    if n1n < 1e-9 || n2n < 1e-9 || b2n < 1e-9 {
        return None;
    }
    // Atan2 form is numerically stable across the full range.
    let m = n1.cross(&(b2 / b2n));
    let x = n1.dot(&n2);
    let y = m.dot(&n2);
    Some((-y).atan2(x).to_degrees())
}

/// Improper-dihedral helper — identical maths to [`dihedral`] but
/// named for the chirality / planarity use-case where the four atoms
/// are not a consecutive bonded chain.
pub fn improper_dihedral(
    a: &Point3<f64>,
    b: &Point3<f64>,
    c: &Point3<f64>,
    d: &Point3<f64>,
) -> Option<f64> {
    dihedral(a, b, c, d)
}

/// The φ/ψ/ω backbone torsions of one residue.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct BackboneTorsions {
    /// φ — `C(i-1)–N(i)–CA(i)–C(i)`. `None` for the chain N-terminus.
    pub phi: Option<f64>,
    /// ψ — `N(i)–CA(i)–C(i)–N(i+1)`. `None` for the C-terminus.
    pub psi: Option<f64>,
    /// ω — `CA(i)–C(i)–N(i+1)–CA(i+1)`. `None` for the C-terminus.
    pub omega: Option<f64>,
}

/// Compute φ/ψ/ω for residue at index `i` within `chain`'s residue
/// list. Uses the highest-occupancy `N`, `CA`, `C` atoms of the
/// relevant residues; any missing atom yields `None` for the torsions
/// that need it.
pub fn backbone_torsions(chain: &Chain, i: usize) -> BackboneTorsions {
    let residues = &chain.residues;
    let mut out = BackboneTorsions::default();
    let cur = match residues.get(i) {
        Some(r) if r.is_amino_acid() => r,
        _ => return out,
    };
    let n = bb_atom(cur, "N");
    let ca = bb_atom(cur, "CA");
    let c = bb_atom(cur, "C");

    // φ needs C of the previous residue.
    if let (Some(prev), Some(n), Some(ca), Some(c)) =
        (i.checked_sub(1).and_then(|p| residues.get(p)), n, ca, c)
    {
        if prev.is_amino_acid() {
            if let Some(c_prev) = bb_atom(prev, "C") {
                out.phi = dihedral(&c_prev, &n, &ca, &c);
            }
        }
    }

    // ψ and ω need N (and CA) of the next residue.
    if let (Some(next), Some(n), Some(ca), Some(c)) = (residues.get(i + 1), n, ca, c) {
        if next.is_amino_acid() {
            if let Some(n_next) = bb_atom(next, "N") {
                out.psi = dihedral(&n, &ca, &c, &n_next);
                if let Some(ca_next) = bb_atom(next, "CA") {
                    out.omega = dihedral(&ca, &c, &n_next, &ca_next);
                }
            }
        }
    }
    out
}

/// φ/ψ for every residue of a chain, paired with the residue index.
/// Residues with both angles defined are the ones a Ramachandran plot
/// uses.
pub fn chain_phi_psi(chain: &Chain) -> Vec<(usize, BackboneTorsions)> {
    (0..chain.residues.len())
        .map(|i| (i, backbone_torsions(chain, i)))
        .collect()
}

/// Highest-occupancy backbone atom coordinate by name.
fn bb_atom(r: &Residue, name: &str) -> Option<Point3<f64>> {
    r.primary_atom(name).map(|a| a.coord)
}

/// The χ sidechain-torsion atom quadruples for each amino acid, by
/// three-letter code. Each inner slice is one χ angle's four atom
/// names in `χ1, χ2, …` order. Standard IUPAC definitions.
fn chi_definitions(resname: &str) -> &'static [&'static [&'static str]] {
    match resname {
        "ARG" => &[
            &["N", "CA", "CB", "CG"],
            &["CA", "CB", "CG", "CD"],
            &["CB", "CG", "CD", "NE"],
            &["CG", "CD", "NE", "CZ"],
        ],
        "ASN" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "OD1"]],
        "ASP" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "OD1"]],
        "CYS" => &[&["N", "CA", "CB", "SG"]],
        "GLN" => &[
            &["N", "CA", "CB", "CG"],
            &["CA", "CB", "CG", "CD"],
            &["CB", "CG", "CD", "OE1"],
        ],
        "GLU" => &[
            &["N", "CA", "CB", "CG"],
            &["CA", "CB", "CG", "CD"],
            &["CB", "CG", "CD", "OE1"],
        ],
        "HIS" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "ND1"]],
        "ILE" => &[&["N", "CA", "CB", "CG1"], &["CA", "CB", "CG1", "CD1"]],
        "LEU" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "CD1"]],
        "LYS" => &[
            &["N", "CA", "CB", "CG"],
            &["CA", "CB", "CG", "CD"],
            &["CB", "CG", "CD", "CE"],
            &["CG", "CD", "CE", "NZ"],
        ],
        "MET" => &[
            &["N", "CA", "CB", "CG"],
            &["CA", "CB", "CG", "SD"],
            &["CB", "CG", "SD", "CE"],
        ],
        "PHE" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "CD1"]],
        "PRO" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "CD"]],
        "SER" => &[&["N", "CA", "CB", "OG"]],
        "THR" => &[&["N", "CA", "CB", "OG1"]],
        "TRP" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "CD1"]],
        "TYR" => &[&["N", "CA", "CB", "CG"], &["CA", "CB", "CG", "CD1"]],
        "VAL" => &[&["N", "CA", "CB", "CG1"]],
        // ALA and GLY have no rotatable sidechain torsions.
        _ => &[],
    }
}

/// Compute the χ sidechain torsions of `residue`, in degrees.
/// `chi[k]` is χ(k+1). A χ whose atoms are not all present is
/// omitted, so a residue with a missing tip atom yields a shorter
/// vector.
pub fn sidechain_chi(residue: &Residue) -> Vec<f64> {
    let defs = chi_definitions(&residue.name);
    let mut out = Vec::new();
    for quad in defs {
        let pts: Option<Vec<Point3<f64>>> = quad
            .iter()
            .map(|name| residue.primary_atom(name).map(|a| a.coord))
            .collect();
        match pts {
            Some(p) if p.len() == 4 => {
                if let Some(d) = dihedral(&p[0], &p[1], &p[2], &p[3]) {
                    out.push(d);
                } else {
                    break; // degenerate — stop the χ chain
                }
            }
            _ => break, // a missing atom truncates the χ chain
        }
    }
    out
}

/// Number of χ torsions a residue *type* defines (independent of
/// which atoms are present).
pub fn chi_count(resname: &str) -> usize {
    chi_definitions(resname).len()
}

/// Bond length of a residue's named atom pair, ångström.
pub fn bond_length(residue: &Residue, a: &str, b: &str) -> Option<f64> {
    let pa = residue.primary_atom(a)?;
    let pb = residue.primary_atom(b)?;
    Some(pa.distance(pb))
}

/// Convenience: a `Vector3` cross-product helper used by callers that
/// want a face normal from three atom coordinates.
pub fn plane_normal(a: &Point3<f64>, b: &Point3<f64>, c: &Point3<f64>) -> Option<Vector3<f64>> {
    let n = (b - a).cross(&(c - a));
    let len = n.norm();
    if len < 1e-9 {
        None
    } else {
        Some(n / len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Atom;

    fn at(name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom::new(name, "C", Point3::new(x, y, z))
    }

    #[test]
    fn right_angle() {
        let a = Point3::new(1.0, 0.0, 0.0);
        let b = Point3::new(0.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        assert!((bond_angle(&a, &b, &c).unwrap() - 90.0).abs() < 1e-9);
    }

    #[test]
    fn dihedral_cis_trans() {
        // planar trans arrangement -> ~180 degrees
        let a = Point3::new(0.0, 1.0, 0.0);
        let b = Point3::new(0.0, 0.0, 0.0);
        let c = Point3::new(1.0, 0.0, 0.0);
        let d = Point3::new(1.0, -1.0, 0.0);
        let v = dihedral(&a, &b, &c, &d).unwrap();
        assert!(v.abs() > 179.0, "trans dihedral was {v}");

        // 90-degree twist
        let d2 = Point3::new(1.0, 0.0, 1.0);
        let v2 = dihedral(&a, &b, &c, &d2).unwrap();
        assert!((v2.abs() - 90.0).abs() < 1e-6, "twisted dihedral was {v2}");
    }

    #[test]
    fn degenerate_returns_none() {
        let o = Point3::origin();
        assert!(bond_angle(&o, &o, &Point3::new(1.0, 0.0, 0.0)).is_none());
        assert!(dihedral(&o, &o, &o, &o).is_none());
    }

    #[test]
    fn backbone_torsions_middle_residue() {
        // 3-residue chain; build an arbitrary but non-degenerate
        // geometry so the middle residue has all of phi/psi/omega.
        let mut chain = Chain::new("A");
        let coords = [
            // residue 1: N, CA, C
            ("N", 0.0, 0.0, 0.0),
            ("CA", 1.5, 0.0, 0.0),
            ("C", 2.0, 1.4, 0.0),
            // residue 2
            ("N", 3.3, 1.6, 0.3),
            ("CA", 4.0, 2.8, 0.7),
            ("C", 5.4, 2.6, 1.2),
            // residue 3
            ("N", 6.0, 3.7, 1.6),
            ("CA", 7.4, 3.8, 2.0),
            ("C", 8.0, 5.0, 2.3),
        ];
        for (ri, group) in coords.chunks(3).enumerate() {
            let mut r = Residue::new("ALA", ri as i32 + 1);
            for (name, x, y, z) in group {
                r.atoms.push(at(name, *x, *y, *z));
            }
            chain.residues.push(r);
        }
        let mid = backbone_torsions(&chain, 1);
        assert!(mid.phi.is_some(), "middle phi missing");
        assert!(mid.psi.is_some(), "middle psi missing");
        assert!(mid.omega.is_some(), "middle omega missing");

        // The N-terminus has no phi; the C-terminus has no psi.
        assert!(backbone_torsions(&chain, 0).phi.is_none());
        assert!(backbone_torsions(&chain, 2).psi.is_none());
    }

    #[test]
    fn chi_definitions_count() {
        assert_eq!(chi_count("ARG"), 4);
        assert_eq!(chi_count("SER"), 1);
        assert_eq!(chi_count("ALA"), 0);
        assert_eq!(chi_count("GLY"), 0);
    }

    #[test]
    fn sidechain_chi_serine() {
        let mut r = Residue::new("SER", 1);
        r.atoms.push(at("N", 0.0, 0.0, 0.0));
        r.atoms.push(at("CA", 1.5, 0.0, 0.0));
        r.atoms.push(at("CB", 2.0, 1.4, 0.0));
        r.atoms.push(at("OG", 3.0, 1.4, 1.0));
        let chi = sidechain_chi(&r);
        assert_eq!(chi.len(), 1);
        assert!(chi[0].is_finite());
    }

    #[test]
    fn sidechain_chi_truncates_on_missing_atom() {
        // ARG with only the first chi's atoms -> exactly one chi.
        let mut r = Residue::new("ARG", 1);
        r.atoms.push(at("N", 0.0, 0.0, 0.0));
        r.atoms.push(at("CA", 1.5, 0.0, 0.0));
        r.atoms.push(at("CB", 2.0, 1.4, 0.0));
        r.atoms.push(at("CG", 3.0, 1.4, 1.0));
        assert_eq!(sidechain_chi(&r).len(), 1);
    }

    #[test]
    fn bond_length_works() {
        let mut r = Residue::new("ALA", 1);
        r.atoms.push(at("N", 0.0, 0.0, 0.0));
        r.atoms.push(at("CA", 1.46, 0.0, 0.0));
        assert!((bond_length(&r, "N", "CA").unwrap() - 1.46).abs() < 1e-9);
    }
}
