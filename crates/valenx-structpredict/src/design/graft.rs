//! **Feature 21 — loop / functional-motif grafting.**
//!
//! A common design task is to take a known functional motif — an
//! active-site loop, a binding epitope, a metal-coordinating helix —
//! and **graft** it onto a stable scaffold protein, so the scaffold
//! displays the motif's function. This is the design strategy behind
//! many engineered enzymes and binders.
//!
//! The geometric heart of grafting is **anchor superposition**: the
//! motif's flanking residues (the few residues either side of the
//! functional part) are superposed — by Kabsch — onto the scaffold
//! residues they will replace. If those anchors superpose with low
//! RMSD the motif "fits" the scaffold there; the motif is then placed
//! by that transform and its functional residues spliced into the
//! scaffold sequence and coordinates.
//!
//! This module performs that superposition-driven graft and reports
//! the anchor RMSD — the standard graftability metric. A real
//! grafting protocol would additionally scan many scaffold positions
//! and relax the junction; this v1 grafts at a caller-chosen position
//! and reports how well it fits.

use serde::{Deserialize, Serialize};
use valenx_biostruct::kabsch;

use crate::error::{Result, StructPredictError};
use crate::model::{ModelResidue, ProteinModel};

/// The outcome of grafting a motif onto a scaffold.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GraftResult {
    /// The scaffold with the motif grafted in.
    pub model: ProteinModel,
    /// Cα-RMSD of the motif anchors onto the scaffold anchors after
    /// superposition — the graftability metric. Below ~1 Å is a
    /// clean graft.
    pub anchor_rmsd: f64,
    /// Half-open `[start, end)` scaffold range the motif replaced.
    pub grafted_range: (usize, usize),
}

impl GraftResult {
    /// Whether the graft is geometrically clean (anchor RMSD below a
    /// 1.5 Å tolerance).
    pub fn is_clean(&self) -> bool {
        self.anchor_rmsd < 1.5
    }
}

/// Grafts a motif onto a scaffold at a chosen position.
///
/// `scaffold` is the host structure; `motif` is the structural motif
/// to graft; `scaffold_range` is the half-open `[start, end)` range of
/// scaffold residues the motif replaces; `anchor` is how many
/// residues at each end of the motif are flanking *anchors* (used for
/// superposition, not themselves functional).
///
/// The motif's `anchor` flanking residues at each end are Kabsch-
/// superposed onto the scaffold residues just outside
/// `scaffold_range`; the whole motif is moved by that transform and
/// its residues spliced in.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a degenerate range, too few
/// anchor residues, or motifs/scaffolds lacking the needed Cα atoms.
pub fn graft_motif(
    scaffold: &ProteinModel,
    motif: &ProteinModel,
    scaffold_range: (usize, usize),
    anchor: usize,
) -> Result<GraftResult> {
    let (start, end) = scaffold_range;
    if end <= start || end > scaffold.residues.len() {
        return Err(StructPredictError::invalid(
            "scaffold_range",
            "not a valid range",
        ));
    }
    if anchor == 0 {
        return Err(StructPredictError::invalid(
            "anchor",
            "need at least 1 anchor residue per side",
        ));
    }
    if motif.residues.len() < 2 * anchor + 1 {
        return Err(StructPredictError::invalid(
            "motif",
            "motif is shorter than its two anchors plus a body",
        ));
    }
    // Scaffold anchors: the `anchor` residues just before `start` and
    // just after `end-1`.
    if start < anchor || end + anchor > scaffold.residues.len() {
        return Err(StructPredictError::invalid(
            "scaffold_range",
            "not enough scaffold residues flanking the graft site for anchors",
        ));
    }

    // Collect the anchor Cα pairs: motif N-anchor ↔ scaffold before,
    // motif C-anchor ↔ scaffold after.
    let mut motif_ca = Vec::new();
    let mut scaffold_ca = Vec::new();
    for k in 0..anchor {
        // N-side.
        let m_idx = k;
        let s_idx = start - anchor + k;
        push_pair(motif, scaffold, m_idx, s_idx, &mut motif_ca, &mut scaffold_ca)?;
        // C-side.
        let m_idx = motif.residues.len() - anchor + k;
        let s_idx = end + k;
        push_pair(motif, scaffold, m_idx, s_idx, &mut motif_ca, &mut scaffold_ca)?;
    }

    let sup = kabsch(&motif_ca, &scaffold_ca)
        .map_err(|e| StructPredictError::invalid("graft_superpose", e.to_string()))?;

    // Build the grafted model: scaffold[..start] + moved motif body +
    // scaffold[end..].
    let mut residues: Vec<ModelResidue> = scaffold.residues[..start].to_vec();
    // The motif body is everything between the anchors.
    for body_res in &motif.residues[anchor..motif.residues.len() - anchor] {
        let mut r = body_res.clone();
        // Move every atom by the anchor superposition transform.
        for p in [&mut r.n, &mut r.ca, &mut r.c, &mut r.o, &mut r.cb]
            .into_iter()
            .flatten()
        {
            *p = sup.transform.apply(p);
        }
        residues.push(r);
    }
    residues.extend_from_slice(&scaffold.residues[end..]);

    let grafted_end = start + (motif.residues.len() - 2 * anchor);
    Ok(GraftResult {
        model: ProteinModel { residues },
        anchor_rmsd: sup.rmsd,
        grafted_range: (start, grafted_end),
    })
}

/// Pushes a motif/scaffold Cα pair onto the superposition lists.
fn push_pair(
    motif: &ProteinModel,
    scaffold: &ProteinModel,
    m_idx: usize,
    s_idx: usize,
    motif_ca: &mut Vec<nalgebra::Point3<f64>>,
    scaffold_ca: &mut Vec<nalgebra::Point3<f64>>,
) -> Result<()> {
    let m_ca = motif
        .residues
        .get(m_idx)
        .and_then(|r| r.ca)
        .ok_or_else(|| StructPredictError::invalid("motif", "anchor residue lacks a Cα"))?;
    let s_ca = scaffold
        .residues
        .get(s_idx)
        .and_then(|r| r.ca)
        .ok_or_else(|| {
            StructPredictError::invalid("scaffold", "anchor residue lacks a Cα")
        })?;
    motif_ca.push(m_ca);
    scaffold_ca.push(s_ca);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn linear_model(seq: &str, dx: f64, dy: f64) -> ProteinModel {
        let mut m = ProteinModel::from_sequence(seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new(i as f64 * 3.8 + dx, dy, 0.0));
            r.n = Some(Point3::new(i as f64 * 3.8 + dx - 0.7, dy + 0.5, 0.0));
            r.c = Some(Point3::new(i as f64 * 3.8 + dx + 0.7, dy + 0.5, 0.0));
            r.o = Some(Point3::new(i as f64 * 3.8 + dx + 0.7, dy + 1.5, 0.0));
        }
        m
    }

    #[test]
    fn grafting_a_matching_motif_is_clean() {
        // Scaffold and motif are both linear; the motif's anchors
        // line up with the scaffold's graft-site flanks.
        let scaffold = linear_model("AAAAAAAAAAAA", 0.0, 0.0);
        // Motif: 2 anchors + 3 body + 2 anchors = 7 residues, geometry
        // matching the scaffold span it replaces.
        let motif = linear_model("GGWWWGG", 3.0 * 3.8, 0.0);
        let res = graft_motif(&scaffold, &motif, (5, 8), 2).expect("graft");
        assert!(res.is_clean(), "anchor RMSD {}", res.anchor_rmsd);
        // The grafted region carries the motif's body residues.
        let grafted: String = res.model.sequence();
        assert!(grafted.contains("WWW"), "motif body grafted: {grafted}");
    }

    #[test]
    fn graft_replaces_the_target_range() {
        let scaffold = linear_model("AAAAAAAAAAAA", 0.0, 0.0);
        let motif = linear_model("GGWWWGG", 3.0 * 3.8, 0.0);
        let res = graft_motif(&scaffold, &motif, (5, 8), 2).expect("graft");
        // Body is 3 residues; replaced span was 3 → length unchanged.
        assert_eq!(res.model.len(), scaffold.len());
        assert_eq!(res.grafted_range, (5, 8));
    }

    #[test]
    fn degenerate_inputs_rejected() {
        let scaffold = linear_model("AAAAAAAA", 0.0, 0.0);
        let motif = linear_model("GGWGG", 0.0, 0.0);
        // bad range
        assert!(graft_motif(&scaffold, &motif, (5, 5), 2).is_err());
        // zero anchors
        assert!(graft_motif(&scaffold, &motif, (3, 5), 0).is_err());
        // motif too short for 3 anchors per side
        assert!(graft_motif(&scaffold, &motif, (3, 5), 3).is_err());
    }
}
