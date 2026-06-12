//! The protein-model data type used across the crate.
//!
//! Prediction, refinement and design all operate on a
//! [`ProteinModel`]: an ordered list of [`ModelResidue`]s, each
//! carrying its amino-acid identity, its four backbone atom
//! coordinates (`N`, `CA`, `C`, `O`) and an optional `CB`
//! pseudo-/real sidechain centroid. This is a deliberately small,
//! self-contained model — the homology and ab-initio protocols build
//! and mutate it directly, and it converts to and from the richer
//! [`valenx_biostruct::Structure`] hierarchy at the I/O boundary.
//!
//! Coordinates are in ångström. A residue whose backbone has not yet
//! been built carries `None` for its atom slots; [`ProteinModel`]
//! reports such gaps so the loop-modelling code knows what to fill.

use nalgebra::Point3;
use serde::{Deserialize, Serialize};

use crate::aa::{one_to_three, three_to_one};
use crate::error::{Result, StructPredictError};

/// Idealised peptide-geometry constants (textbook engineering values,
/// ångström / degrees).
pub mod ideal {
    /// N–Cα bond length.
    pub const N_CA: f64 = 1.458;
    /// Cα–C bond length.
    pub const CA_C: f64 = 1.525;
    /// C–N peptide bond length.
    pub const C_N: f64 = 1.329;
    /// C=O carbonyl bond length.
    pub const C_O: f64 = 1.231;
    /// Cα–Cβ bond length.
    pub const CA_CB: f64 = 1.530;
    /// N–Cα–C backbone bond angle, degrees.
    pub const N_CA_C: f64 = 111.0;
    /// Cα–C–N bond angle, degrees.
    pub const CA_C_N: f64 = 116.2;
    /// C–N–Cα bond angle, degrees.
    pub const C_N_CA: f64 = 121.7;
    /// The peptide-bond ω dihedral (trans), degrees.
    pub const OMEGA_TRANS: f64 = 180.0;
    /// Virtual Cα–Cα distance for a trans peptide unit.
    pub const CA_CA: f64 = 3.80;
}

/// One residue of a [`ProteinModel`]: an amino-acid identity plus its
/// backbone (and optional Cβ) coordinates.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModelResidue {
    /// One-letter amino-acid code (`A`..`V`; `X` for an unknown
    /// position whose identity is not yet assigned).
    pub aa: char,
    /// Backbone amide-nitrogen coordinate, if built.
    pub n: Option<Point3<f64>>,
    /// Backbone α-carbon coordinate, if built.
    pub ca: Option<Point3<f64>>,
    /// Backbone carbonyl-carbon coordinate, if built.
    pub c: Option<Point3<f64>>,
    /// Backbone carbonyl-oxygen coordinate, if built.
    pub o: Option<Point3<f64>>,
    /// β-carbon coordinate (a real Cβ when copied from a template, a
    /// rebuilt pseudo-Cβ otherwise). `None` for glycine or an
    /// unbuilt residue.
    pub cb: Option<Point3<f64>>,
}

impl ModelResidue {
    /// A residue with the given amino-acid identity and no
    /// coordinates (a gap to be built).
    pub fn empty(aa: char) -> Self {
        ModelResidue {
            aa,
            n: None,
            ca: None,
            c: None,
            o: None,
            cb: None,
        }
    }

    /// `true` when the residue has its full `N`, `CA`, `C`, `O`
    /// backbone built.
    pub fn has_backbone(&self) -> bool {
        self.n.is_some() && self.ca.is_some() && self.c.is_some() && self.o.is_some()
    }

    /// The residue's three-letter PDB name (`"ALA"`, …; `"UNK"` for
    /// an unknown identity).
    pub fn resname(&self) -> &'static str {
        one_to_three(self.aa).unwrap_or("UNK")
    }
}

/// A protein model: an ordered chain of [`ModelResidue`]s.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProteinModel {
    /// Residues in N-to-C chain order.
    pub residues: Vec<ModelResidue>,
}

impl ProteinModel {
    /// An empty model with no residues.
    pub fn new() -> Self {
        ProteinModel {
            residues: Vec::new(),
        }
    }

    /// Builds a coordinate-free model from a one-letter sequence —
    /// every residue is a gap awaiting a backbone.
    ///
    /// # Errors
    /// [`StructPredictError::Invalid`] if the sequence is empty.
    pub fn from_sequence(seq: &str) -> Result<Self> {
        let seq = seq.trim();
        if seq.is_empty() {
            return Err(StructPredictError::invalid("sequence", "empty"));
        }
        let residues = seq
            .chars()
            .map(|c| {
                let up = c.to_ascii_uppercase();
                let aa = if crate::aa::aa_index(up).is_some() {
                    up
                } else {
                    'X'
                };
                ModelResidue::empty(aa)
            })
            .collect();
        Ok(ProteinModel { residues })
    }

    /// Residue count.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// `true` when the model has no residues.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// The model's one-letter sequence.
    pub fn sequence(&self) -> String {
        self.residues.iter().map(|r| r.aa).collect()
    }

    /// Every built Cα coordinate, in chain order. Residues without a
    /// Cα are skipped.
    pub fn ca_trace(&self) -> Vec<Point3<f64>> {
        self.residues.iter().filter_map(|r| r.ca).collect()
    }

    /// `true` when every residue has a complete backbone.
    pub fn is_complete(&self) -> bool {
        !self.residues.is_empty() && self.residues.iter().all(|r| r.has_backbone())
    }

    /// Half-open `[start, end)` index ranges of residues lacking a
    /// backbone — the gaps the loop modeller must close.
    pub fn gaps(&self) -> Vec<(usize, usize)> {
        let mut gaps = Vec::new();
        let mut start: Option<usize> = None;
        for (i, r) in self.residues.iter().enumerate() {
            match (r.has_backbone(), start) {
                (false, None) => start = Some(i),
                (true, Some(s)) => {
                    gaps.push((s, i));
                    start = None;
                }
                _ => {}
            }
        }
        if let Some(s) = start {
            gaps.push((s, self.residues.len()));
        }
        gaps
    }

    /// Translates every built atom so the model's Cα centroid sits at
    /// the origin. A no-op on a backbone-free model.
    pub fn center(&mut self) {
        let trace = self.ca_trace();
        if trace.is_empty() {
            return;
        }
        let mut c = nalgebra::Vector3::zeros();
        for p in &trace {
            c += p.coords;
        }
        c /= trace.len() as f64;
        for r in &mut self.residues {
            for p in [&mut r.n, &mut r.ca, &mut r.c, &mut r.o, &mut r.cb]
                .into_iter()
                .flatten()
            {
                *p -= c;
            }
        }
    }

    /// Builds a model from a chain of a [`valenx_biostruct::Structure`].
    ///
    /// Only the polymer amino-acid residues of the named chain are
    /// taken; their `N`/`CA`/`C`/`O`/`CB` atoms populate the model.
    ///
    /// # Errors
    /// [`StructPredictError::NotFound`] if the chain is absent or has
    /// no amino-acid residues.
    pub fn from_structure_chain(
        structure: &valenx_biostruct::Structure,
        chain_id: &str,
    ) -> Result<Self> {
        let model = structure.first_model();
        let chain = model
            .chain(chain_id)
            .ok_or_else(|| StructPredictError::not_found("chain", format!("chain `{chain_id}`")))?;
        let mut residues = Vec::new();
        for res in &chain.residues {
            if !res.is_amino_acid() {
                continue;
            }
            let aa = three_to_one(&res.name).unwrap_or('X');
            let pos = |name: &str| res.atom(name).map(|a| a.coord);
            residues.push(ModelResidue {
                aa,
                n: pos("N"),
                ca: pos("CA"),
                c: pos("C"),
                o: pos("O"),
                cb: pos("CB"),
            });
        }
        if residues.is_empty() {
            return Err(StructPredictError::not_found(
                "chain",
                format!("chain `{chain_id}` has no amino-acid residues"),
            ));
        }
        Ok(ProteinModel { residues })
    }

    /// Renders the model as PDB `ATOM` records (single chain `A`).
    /// Only built atoms are written; unbuilt residues contribute no
    /// records. The output round-trips through
    /// [`valenx_biostruct::read_structure`].
    pub fn to_pdb(&self) -> String {
        let mut out = String::new();
        let mut serial = 1;
        for (i, r) in self.residues.iter().enumerate() {
            let resseq = (i + 1) as i32;
            let resname = r.resname();
            let mut emit = |name: &str, p: &Point3<f64>, serial: &mut i32, elem: &str| {
                out.push_str(&format!(
                    "ATOM  {:>5} {:<4} {:<3} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00          {:>2}\n",
                    serial, name, resname, resseq, p.x, p.y, p.z, elem
                ));
                *serial += 1;
            };
            if let Some(p) = &r.n {
                emit(" N", p, &mut serial, "N");
            }
            if let Some(p) = &r.ca {
                emit(" CA", p, &mut serial, "C");
            }
            if let Some(p) = &r.c {
                emit(" C", p, &mut serial, "C");
            }
            if let Some(p) = &r.o {
                emit(" O", p, &mut serial, "O");
            }
            if let Some(p) = &r.cb {
                emit(" CB", p, &mut serial, "C");
            }
        }
        out.push_str("END\n");
        out
    }
}

impl Default for ProteinModel {
    fn default() -> Self {
        ProteinModel::new()
    }
}

/// Builds a full backbone for a model from a per-residue list of
/// `(phi, psi)` dihedral angles (degrees), using idealised bond
/// lengths and angles and a trans ω peptide bond.
///
/// `torsions[i]` is the `(φ, ψ)` of residue `i`. The chain is grown
/// N-to-C from a canonical seed; every residue receives `N`, `CA`,
/// `C` and `O`. Cβ atoms are *not* placed here — call the rotamer /
/// sidechain code for those.
///
/// # Errors
/// [`StructPredictError::Invalid`] if `torsions.len()` disagrees with
/// the model's residue count, or the model is empty.
pub fn build_backbone_from_torsions(
    model: &mut ProteinModel,
    torsions: &[(f64, f64)],
) -> Result<()> {
    let n = model.residues.len();
    if n == 0 {
        return Err(StructPredictError::invalid("model", "empty"));
    }
    if torsions.len() != n {
        return Err(StructPredictError::invalid(
            "torsions",
            format!("length {} disagrees with {n} residues", torsions.len()),
        ));
    }
    // Seed residue 0: N at origin, CA along +x, C placed by the
    // backbone bond angle.
    let n0 = Point3::new(0.0, 0.0, 0.0);
    let ca0 = Point3::new(ideal::N_CA, 0.0, 0.0);
    // A virtual atom before N0 so the seed C has a dihedral frame.
    let virt = Point3::new(-1.0, 1.0, 0.0);
    let c0 = place_atom(
        &virt,
        &n0,
        &ca0,
        ideal::CA_C,
        ideal::N_CA_C.to_radians(),
        torsions[0].1.to_radians(), // ψ of residue 0
    );
    {
        let r = &mut model.residues[0];
        r.n = Some(n0);
        r.ca = Some(ca0);
        r.c = Some(c0);
    }
    for (i, &(phi, _psi)) in torsions.iter().enumerate().skip(1) {
        let prev = model.residues[i - 1].clone();
        let (pn, pca, pc) = (
            prev.n.expect("prev N"),
            prev.ca.expect("prev CA"),
            prev.c.expect("prev C"),
        );
        // Standard NeRF backbone placement — each atom's dihedral is the
        // four-atom torsion ending at that atom:
        //   N(i)  : dihedral N(i-1)-CA(i-1)-C(i-1)-N(i)  = ψ(i-1)
        //   CA(i) : dihedral CA(i-1)-C(i-1)-N(i)-CA(i)   = ω  (trans, 180)
        //   C(i)  : dihedral C(i-1)-N(i)-CA(i)-C(i)      = φ(i)
        // (The earlier code drove N(i) by ω, CA(i) by φ and C(i) by ψ —
        // shifted by one atom — so the built backbone did not have the
        // requested φ/ψ at all.)
        //
        // `place_atom` uses the opposite dihedral-sign convention to
        // `refine::ramachandran::dihedral_deg` (the one the rest of the
        // crate measures φ/ψ with), so the requested torsions are negated
        // here to make build → measure a faithful round-trip. ω = 180 is
        // its own negative, so it is unaffected.
        let psi_prev = torsions[i - 1].1;
        let n_i = place_atom(
            &pn,
            &pca,
            &pc,
            ideal::C_N,
            ideal::CA_C_N.to_radians(),
            (-psi_prev).to_radians(),
        );
        let ca_i = place_atom(
            &pca,
            &pc,
            &n_i,
            ideal::N_CA,
            ideal::C_N_CA.to_radians(),
            ideal::OMEGA_TRANS.to_radians(),
        );
        let c_i = place_atom(
            &pc,
            &n_i,
            &ca_i,
            ideal::CA_C,
            ideal::N_CA_C.to_radians(),
            (-phi).to_radians(),
        );
        let r = &mut model.residues[i];
        r.n = Some(n_i);
        r.ca = Some(ca_i);
        r.c = Some(c_i);
    }
    // Carbonyl oxygens, in the peptide plane.
    for r in &mut model.residues {
        if let (Some(ni), Some(cai), Some(ci)) = (r.n, r.ca, r.c) {
            r.o = Some(place_atom(
                &ni,
                &cai,
                &ci,
                ideal::C_O,
                (180.0 - ideal::CA_C_N).to_radians(),
                0.0,
            ));
        }
    }
    Ok(())
}

/// Places a fourth point `d` given three anchor points and the
/// internal coordinates (bond length `d`–`c`, bond angle `b`–`c`–`d`,
/// dihedral `a`–`b`–`c`–`d` in radians). The standard NeRF
/// (Natural-Extension Reference Frame) atom-placement formula used to
/// grow a backbone from φ/ψ angles.
pub fn place_atom(
    a: &Point3<f64>,
    b: &Point3<f64>,
    c: &Point3<f64>,
    bond: f64,
    angle: f64,
    dihedral: f64,
) -> Point3<f64> {
    let bc = (c - b).normalize();
    let ab = b - a;
    // n is perpendicular to the a-b-c plane.
    let n = ab.cross(&bc).normalize();
    let m = n.cross(&bc);
    // Local coordinate: bond projected by the angle and dihedral.
    let d2 = nalgebra::Vector3::new(
        -bond * angle.cos(),
        bond * angle.sin() * dihedral.cos(),
        bond * angle.sin() * dihedral.sin(),
    );
    // Basis [bc, m, n] maps local → world.
    let world = bc * d2.x + m * d2.y + n * d2.z;
    c + world
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_sequence_marks_unknowns() {
        let m = ProteinModel::from_sequence("ACDXZ").expect("model");
        assert_eq!(m.len(), 5);
        // Z is not standard → X; X stays X.
        assert_eq!(m.sequence(), "ACDXX");
        assert!(!m.is_complete());
    }

    #[test]
    fn empty_sequence_rejected() {
        assert!(ProteinModel::from_sequence("   ").is_err());
    }

    #[test]
    fn gaps_detected() {
        let mut m = ProteinModel::from_sequence("AAAAA").expect("model");
        // build residues 0 and 4, leave 1..4 as a gap
        for &i in &[0usize, 4] {
            let r = &mut m.residues[i];
            r.n = Some(Point3::origin());
            r.ca = Some(Point3::origin());
            r.c = Some(Point3::origin());
            r.o = Some(Point3::origin());
        }
        assert_eq!(m.gaps(), vec![(1, 4)]);
    }

    #[test]
    fn build_backbone_from_torsions_produces_sane_geometry() {
        let mut m = ProteinModel::from_sequence("AAAAAAAA").expect("model");
        // All-helical torsions.
        let torsions = vec![(-63.0, -42.0); 8];
        build_backbone_from_torsions(&mut m, &torsions).expect("build");
        assert!(m.is_complete());
        // The N–CA bond length is idealised.
        let r = &m.residues[4];
        let d = (r.ca.unwrap() - r.n.unwrap()).norm();
        assert!((d - ideal::N_CA).abs() < 1e-6, "N-CA = {d}");
        // An α-helix has ~1.5 Å rise per residue → the chain is
        // markedly shorter end-to-end than an extended chain.
        let end_to_end = (m.residues[7].ca.unwrap() - m.residues[0].ca.unwrap()).norm();
        assert!(end_to_end < 7.0 * CA_GAP, "helix is compact: {end_to_end}");
    }

    const CA_GAP: f64 = ideal::CA_CA;

    #[test]
    fn build_backbone_rejects_length_mismatch() {
        let mut m = ProteinModel::from_sequence("AAA").expect("model");
        assert!(build_backbone_from_torsions(&mut m, &[(-63.0, -42.0)]).is_err());
    }

    #[test]
    fn place_atom_respects_bond_length() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.5, 0.0, 0.0);
        let c = Point3::new(3.0, 1.0, 0.0);
        let d = place_atom(&a, &b, &c, 1.33, 2.0, 1.0);
        let len = (d - c).norm();
        assert!((len - 1.33).abs() < 1e-9, "bond length {len}");
    }
}
