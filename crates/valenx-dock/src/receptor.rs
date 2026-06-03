//! Rigid receptor: flat list of atoms with AD4 types and partial charges.

use nalgebra::Vector3;
use valenx_bio::format::pdbqt::{parse, PdbqtRecord};

use crate::atom_type::Ad4AtomType;
use crate::error::DockError;

/// One receptor atom in dock-engine form.
#[derive(Clone, Debug)]
pub struct ReceptorAtom {
    /// Position in world coordinates (Å).
    pub position: Vector3<f64>,
    /// Vina/AD4 atom type.
    pub ad4_type: Ad4AtomType,
    /// Gasteiger partial charge.
    pub partial_charge: f64,
}

/// Rigid receptor — no flexibility in v1.
#[derive(Clone, Debug, Default)]
pub struct Receptor {
    /// All receptor atoms.
    pub atoms: Vec<ReceptorAtom>,
}

impl Receptor {
    /// Parse a PDBQT receptor file.
    pub fn from_pdbqt(text: &str) -> Result<Self, DockError> {
        let records = parse(text).map_err(DockError::from)?;
        let mut atoms = Vec::new();
        for r in records {
            if let PdbqtRecord::Atom(a) = r {
                let ad4: Ad4AtomType = a.ad4_type.parse().map_err(DockError::AtomType)?;
                atoms.push(ReceptorAtom {
                    position: a.position,
                    ad4_type: ad4,
                    partial_charge: a.partial_charge,
                });
            }
        }
        if atoms.is_empty() {
            return Err(DockError::EmptyReceptor);
        }
        Ok(Self { atoms })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_two_atom_receptor() {
        let pdbqt = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ATOM      2  OG1 SER A   2       2.500   0.000   0.000  1.00  0.00    -0.300 OA
";
        let r = Receptor::from_pdbqt(pdbqt).unwrap();
        assert_eq!(r.atoms.len(), 2);
        assert_eq!(r.atoms[0].ad4_type, Ad4AtomType::C);
        assert_eq!(r.atoms[1].ad4_type, Ad4AtomType::OA);
    }

    #[test]
    fn rejects_empty() {
        let err = Receptor::from_pdbqt("REMARK only a comment\n").unwrap_err();
        assert!(matches!(err, DockError::EmptyReceptor));
    }
}
