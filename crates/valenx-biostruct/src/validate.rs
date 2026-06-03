//! Structure-validation report.
//!
//! [`validate_structure`] bundles the individual geometry checks of
//! this crate into one [`ValidationReport`]: steric clashes,
//! Ramachandran outliers, bond-length / bond-angle outliers, and
//! missing-atom detection. It is the MolProbity-style "is this model
//! geometrically sane" summary.
//!
//! For the higher-level per-structure *analysis* bundle (chain
//! summary, composition, secondary-structure content, radius of
//! gyration and a validation verdict) plus multi-structure batch
//! utilities, see the [`analyze`](crate::analyze) module — that
//! report embeds a [`ValidationReport`].
//!
//! ## Scope of this v1
//!
//! Bond-length / bond-angle outlier thresholds use simple ideal
//! values with fixed tolerances, not the full Engh-Huber
//! residue-specific parameter set. Missing-atom detection compares
//! against a per-residue expected heavy-atom list. The report is a
//! real, useful screen; it is not a substitute for a full
//! MolProbity / wwPDB validation pipeline.

use crate::error::Result;
use crate::geometry::angles::bond_angle;
use crate::geometry::contacts::{detect_clashes, Clash};
use crate::geometry::ramachandran::summarize;
use crate::structure::{Model, Residue, Structure};

/// One flagged geometry outlier.
#[derive(Clone, Debug, PartialEq)]
pub struct GeometryOutlier {
    /// `(chain_id, residue_seq, residue_name)` the outlier belongs to.
    pub residue: (String, i32, String),
    /// What is wrong — a short human-readable description.
    pub description: String,
    /// The measured value (a length, angle, or count).
    pub measured: f64,
    /// The ideal / expected value.
    pub expected: f64,
}

/// A complete structure-validation report.
#[derive(Clone, Debug, PartialEq)]
pub struct ValidationReport {
    /// The structure id the report was run on.
    pub structure_id: String,
    /// Steric clashes found in the first model.
    pub clashes: Vec<Clash>,
    /// Bond-length and bond-angle outliers.
    pub geometry_outliers: Vec<GeometryOutlier>,
    /// Residues whose `(φ, ψ)` lands in a Ramachandran outlier
    /// region, as `(chain_id, residue_seq)`.
    pub ramachandran_outliers: Vec<(String, i32)>,
    /// Residues missing one or more expected heavy atoms, with the
    /// missing atom names.
    pub missing_atoms: Vec<(String, i32, Vec<String>)>,
    /// Total residue count across the model.
    pub residue_count: usize,
    /// Total atom count across the model.
    pub atom_count: usize,
}

impl ValidationReport {
    /// A single 0–100 quality score: starts at 100 and is docked for
    /// each clash, geometry outlier, Ramachandran outlier and
    /// incomplete residue, scaled by the structure size. Higher is
    /// better.
    pub fn quality_score(&self) -> f64 {
        let n = self.residue_count.max(1) as f64;
        let penalties = self.clashes.len() as f64 * 2.0
            + self.geometry_outliers.len() as f64
            + self.ramachandran_outliers.len() as f64 * 1.5
            + self.missing_atoms.len() as f64 * 0.5;
        (100.0 - 100.0 * penalties / n / 4.0).clamp(0.0, 100.0)
    }

    /// Whether the structure passed every check with no flags.
    pub fn is_clean(&self) -> bool {
        self.clashes.is_empty()
            && self.geometry_outliers.is_empty()
            && self.ramachandran_outliers.is_empty()
            && self.missing_atoms.is_empty()
    }

    /// A short multi-line human-readable summary.
    pub fn summary(&self) -> String {
        format!(
            "structure {}: {} residues, {} atoms\n  clashes: {}\n  \
             geometry outliers: {}\n  Ramachandran outliers: {}\n  \
             incomplete residues: {}\n  quality score: {:.1}/100",
            self.structure_id,
            self.residue_count,
            self.atom_count,
            self.clashes.len(),
            self.geometry_outliers.len(),
            self.ramachandran_outliers.len(),
            self.missing_atoms.len(),
            self.quality_score(),
        )
    }
}

/// Run every validation check on a structure's first model.
///
/// `clash_tolerance` is the van der Waals overlap allowance (≈ 0.4 Å
/// is the MolProbity default).
pub fn validate_structure(
    structure: &Structure,
    clash_tolerance: f64,
) -> Result<ValidationReport> {
    let model = structure.first_model();

    let clashes = detect_clashes(model, clash_tolerance)?;
    let geometry_outliers = check_geometry(model);
    let ramachandran_outliers = check_ramachandran(model);
    let missing_atoms = check_missing_atoms(model);

    Ok(ValidationReport {
        structure_id: structure.id.clone(),
        clashes,
        geometry_outliers,
        ramachandran_outliers,
        missing_atoms,
        residue_count: model.residues().count(),
        atom_count: model.atom_count(),
    })
}

/// Bond-length and bond-angle outlier check over backbone geometry.
fn check_geometry(model: &Model) -> Vec<GeometryOutlier> {
    // Ideal protein backbone geometry, with generous tolerances.
    const N_CA: (f64, f64) = (1.46, 0.10);
    const CA_C: (f64, f64) = (1.52, 0.10);
    const C_O: (f64, f64) = (1.23, 0.10);
    const N_CA_C: (f64, f64) = (111.0, 8.0); // angle, degrees

    let mut outliers = Vec::new();
    for chain in &model.chains {
        for residue in &chain.residues {
            if !residue.is_amino_acid() {
                continue;
            }
            let id = (chain.id.clone(), residue.seq_num, residue.name.clone());
            // Bond lengths.
            for (a, b, (ideal, tol), label) in [
                ("N", "CA", N_CA, "N-CA bond"),
                ("CA", "C", CA_C, "CA-C bond"),
                ("C", "O", C_O, "C=O bond"),
            ] {
                if let (Some(pa), Some(pb)) =
                    (residue.primary_atom(a), residue.primary_atom(b))
                {
                    let d = pa.distance(pb);
                    if (d - ideal).abs() > tol {
                        outliers.push(GeometryOutlier {
                            residue: id.clone(),
                            description: format!("{label} length out of range"),
                            measured: d,
                            expected: ideal,
                        });
                    }
                }
            }
            // N-CA-C backbone angle.
            if let (Some(n), Some(ca), Some(c)) = (
                residue.primary_atom("N"),
                residue.primary_atom("CA"),
                residue.primary_atom("C"),
            ) {
                if let Some(angle) = bond_angle(&n.coord, &ca.coord, &c.coord) {
                    if (angle - N_CA_C.0).abs() > N_CA_C.1 {
                        outliers.push(GeometryOutlier {
                            residue: id.clone(),
                            description: "N-CA-C backbone angle out of range".to_string(),
                            measured: angle,
                            expected: N_CA_C.0,
                        });
                    }
                }
            }
        }
    }
    outliers
}

/// Ramachandran-outlier check over every protein chain.
fn check_ramachandran(model: &Model) -> Vec<(String, i32)> {
    let mut out = Vec::new();
    for chain in &model.chains {
        // summarize gives counts; for the residue list re-run the
        // per-residue classification.
        let _ = summarize(chain);
        for point in crate::geometry::ramachandran::chain_ramachandran(chain) {
            if !point.region.is_allowed() {
                if let Some(r) = chain.residues.get(point.residue_index) {
                    out.push((chain.id.clone(), r.seq_num));
                }
            }
        }
    }
    out
}

/// The expected heavy-atom names of a standard amino acid.
fn expected_heavy_atoms(resname: &str) -> &'static [&'static str] {
    match resname {
        "GLY" => &["N", "CA", "C", "O"],
        "ALA" => &["N", "CA", "C", "O", "CB"],
        "SER" => &["N", "CA", "C", "O", "CB", "OG"],
        "CYS" => &["N", "CA", "C", "O", "CB", "SG"],
        "THR" => &["N", "CA", "C", "O", "CB", "OG1", "CG2"],
        "VAL" => &["N", "CA", "C", "O", "CB", "CG1", "CG2"],
        "LEU" => &["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
        "ILE" => &["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
        "PRO" => &["N", "CA", "C", "O", "CB", "CG", "CD"],
        "MET" => &["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
        "ASP" => &["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
        "ASN" => &["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
        "GLU" => &["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
        "GLN" => &["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
        "LYS" => &["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
        "ARG" => &[
            "N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2",
        ],
        "HIS" => &["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        "PHE" => &[
            "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
        ],
        "TYR" => &[
            "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH",
        ],
        "TRP" => &[
            "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3",
            "CZ2", "CZ3", "CH2",
        ],
        _ => &[],
    }
}

/// Missing-atom check: compare each standard amino acid's atoms
/// against its expected heavy-atom list.
fn check_missing_atoms(model: &Model) -> Vec<(String, i32, Vec<String>)> {
    let mut out = Vec::new();
    for chain in &model.chains {
        for residue in &chain.residues {
            let expected = expected_heavy_atoms(&residue.name);
            if expected.is_empty() {
                continue; // not a standard residue we track
            }
            let missing: Vec<String> = expected
                .iter()
                .filter(|name| !has_atom(residue, name))
                .map(|s| s.to_string())
                .collect();
            if !missing.is_empty() {
                out.push((chain.id.clone(), residue.seq_num, missing));
            }
        }
    }
    out
}

/// Whether a residue contains an atom of the given name.
fn has_atom(residue: &Residue, name: &str) -> bool {
    residue.atoms.iter().any(|a| a.name == name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Chain};
    use nalgebra::Point3;

    /// A geometrically clean alanine with correct backbone bonds.
    fn good_alanine(seq: i32, origin: Point3<f64>) -> Residue {
        let mut r = Residue::new("ALA", seq);
        let n = origin;
        let ca = n + nalgebra::Vector3::new(1.46, 0.0, 0.0);
        // CA-C 1.52 A; place C so N-CA-C ~111 deg.
        let c = ca
            + nalgebra::Vector3::new(
                -(1.52 * (111.0_f64.to_radians()).cos()),
                1.52 * (111.0_f64.to_radians()).sin(),
                0.0,
            );
        let o = c + nalgebra::Vector3::new(0.0, 1.23, 0.0);
        let cb = ca + nalgebra::Vector3::new(0.0, -1.53, 0.0);
        r.atoms.push(Atom::new("N", "N", n));
        r.atoms.push(Atom::new("CA", "C", ca));
        r.atoms.push(Atom::new("C", "C", c));
        r.atoms.push(Atom::new("O", "O", o));
        r.atoms.push(Atom::new("CB", "C", cb));
        r
    }

    #[test]
    fn clean_structure_reports_clean() {
        let mut s = Structure::new("CLEAN");
        let mut chain = Chain::new("A");
        chain
            .residues
            .push(good_alanine(1, Point3::new(0.0, 0.0, 0.0)));
        chain
            .residues
            .push(good_alanine(2, Point3::new(20.0, 0.0, 0.0)));
        s.first_model_mut().chains.push(chain);

        let report = validate_structure(&s, 0.4).unwrap();
        assert!(report.clashes.is_empty());
        assert!(report.geometry_outliers.is_empty());
        assert!(report.missing_atoms.is_empty());
        assert_eq!(report.residue_count, 2);
        assert!(report.quality_score() > 95.0);
    }

    #[test]
    fn detects_a_bad_bond_length() {
        let mut s = Structure::new("BADBOND");
        let mut chain = Chain::new("A");
        let mut bad = Residue::new("ALA", 1);
        // N-CA stretched to 2.5 A — far outside tolerance.
        bad.atoms.push(Atom::new("N", "N", Point3::new(0.0, 0.0, 0.0)));
        bad.atoms
            .push(Atom::new("CA", "C", Point3::new(2.5, 0.0, 0.0)));
        bad.atoms
            .push(Atom::new("C", "C", Point3::new(3.5, 1.4, 0.0)));
        bad.atoms
            .push(Atom::new("O", "O", Point3::new(3.5, 2.6, 0.0)));
        bad.atoms
            .push(Atom::new("CB", "C", Point3::new(2.5, -1.5, 0.0)));
        chain.residues.push(bad);
        s.first_model_mut().chains.push(chain);

        let report = validate_structure(&s, 0.4).unwrap();
        assert!(
            report
                .geometry_outliers
                .iter()
                .any(|o| o.description.contains("N-CA")),
            "expected an N-CA bond outlier"
        );
        assert!(!report.is_clean());
    }

    #[test]
    fn detects_missing_atoms() {
        let mut s = Structure::new("INCOMPLETE");
        let mut chain = Chain::new("A");
        // A TRP with only the backbone — many sidechain atoms missing.
        let mut trp = Residue::new("TRP", 1);
        trp.atoms.push(Atom::new("N", "N", Point3::new(0.0, 0.0, 0.0)));
        trp.atoms
            .push(Atom::new("CA", "C", Point3::new(1.46, 0.0, 0.0)));
        trp.atoms
            .push(Atom::new("C", "C", Point3::new(2.0, 1.4, 0.0)));
        trp.atoms
            .push(Atom::new("O", "O", Point3::new(2.0, 2.6, 0.0)));
        chain.residues.push(trp);
        s.first_model_mut().chains.push(chain);

        let report = validate_structure(&s, 0.4).unwrap();
        assert_eq!(report.missing_atoms.len(), 1);
        let (_, seq, missing) = &report.missing_atoms[0];
        assert_eq!(*seq, 1);
        assert!(missing.contains(&"CB".to_string()));
        assert!(missing.contains(&"CG".to_string()));
    }

    #[test]
    fn quality_score_is_bounded() {
        let mut s = Structure::new("X");
        let mut chain = Chain::new("A");
        chain.residues.push(good_alanine(1, Point3::origin()));
        s.first_model_mut().chains.push(chain);
        let report = validate_structure(&s, 0.4).unwrap();
        let q = report.quality_score();
        assert!((0.0..=100.0).contains(&q));
    }

    #[test]
    fn summary_mentions_the_structure() {
        let mut s = Structure::new("MYSTRUCT");
        let mut chain = Chain::new("A");
        chain.residues.push(good_alanine(1, Point3::origin()));
        s.first_model_mut().chains.push(chain);
        let report = validate_structure(&s, 0.4).unwrap();
        assert!(report.summary().contains("MYSTRUCT"));
        assert!(report.summary().contains("quality score"));
    }
}
