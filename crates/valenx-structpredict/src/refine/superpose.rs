//! **Feature 17 — model-vs-reference superposition + RMSD / GDT.**
//!
//! To judge a predicted model you superpose it on the experimental
//! reference structure and measure how far the atoms are apart. This
//! module reuses [`valenx_biostruct`]'s Kabsch optimal superposition
//! and reports the two standard accuracy metrics:
//!
//! - **RMSD** — root-mean-square deviation of the Cα atoms after
//!   optimal rigid superposition. Simple, but dominated by the
//!   worst-modelled regions.
//! - **GDT** (Global Distance Test) — the CASP metric: the fraction
//!   of Cα atoms within a distance cutoff of the reference after
//!   superposition. **GDT-TS** averages the fractions under 1, 2, 4
//!   and 8 Å cutoffs. GDT is robust to a few badly-placed residues
//!   and is the standard structure-prediction accuracy score.
//!
//! Both operate on the equal-length Cα traces of the model and the
//! reference (residue `i` of one corresponds to residue `i` of the
//! other).

use nalgebra::Point3;
use valenx_biostruct::kabsch;

use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;

/// A model-vs-reference superposition report.
#[derive(Clone, Debug, PartialEq)]
pub struct SuperposeReport {
    /// Cα-RMSD after optimal superposition, ångström.
    pub rmsd: f64,
    /// GDT-TS score in `[0, 100]` — higher is more accurate.
    pub gdt_ts: f64,
    /// GDT-HA (high-accuracy) score in `[0, 100]` — the stricter
    /// 0.5/1/2/4 Å variant.
    pub gdt_ha: f64,
    /// Number of corresponded Cα pairs.
    pub n: usize,
}

/// A pair of equal-length coordinate lists.
type TracePair = (Vec<Point3<f64>>, Vec<Point3<f64>>);

/// Extracts the two equal-length Cα traces of a model and a
/// reference, paired by residue index.
///
/// A residue contributes a pair only if *both* the model and the
/// reference have a Cα at that index. Returns the paired coordinate
/// lists.
fn paired_traces(model: &ProteinModel, reference: &ProteinModel) -> Result<TracePair> {
    let n = model.residues.len().min(reference.residues.len());
    let mut a = Vec::new();
    let mut b = Vec::new();
    for i in 0..n {
        if let (Some(ca_m), Some(ca_r)) = (model.residues[i].ca, reference.residues[i].ca) {
            a.push(ca_m);
            b.push(ca_r);
        }
    }
    if a.len() < 3 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 3 corresponded Cα pairs to superpose",
        ));
    }
    Ok((a, b))
}

/// Cα-RMSD of a model against a reference after optimal Kabsch
/// superposition.
///
/// # Errors
/// [`StructPredictError::Invalid`] if fewer than 3 Cα atoms
/// correspond.
pub fn ca_rmsd_superposed(model: &ProteinModel, reference: &ProteinModel) -> Result<f64> {
    let (a, b) = paired_traces(model, reference)?;
    let sup =
        kabsch(&a, &b).map_err(|e| StructPredictError::invalid("superpose", e.to_string()))?;
    Ok(sup.rmsd)
}

/// GDT-TS (Global Distance Test, Total Score) of a model against a
/// reference.
///
/// After Kabsch superposition, GDT-TS is the mean of the fractions of
/// Cα atoms within 1, 2, 4 and 8 Å of the reference, scaled to
/// `[0, 100]`.
///
/// # Errors
/// [`StructPredictError::Invalid`] if fewer than 3 Cα atoms
/// correspond.
pub fn gdt_ts(model: &ProteinModel, reference: &ProteinModel) -> Result<f64> {
    let report = model_vs_reference(model, reference)?;
    Ok(report.gdt_ts)
}

/// Full model-vs-reference superposition: RMSD, GDT-TS and GDT-HA.
///
/// # Errors
/// [`StructPredictError::Invalid`] if fewer than 3 Cα atoms
/// correspond.
pub fn model_vs_reference(
    model: &ProteinModel,
    reference: &ProteinModel,
) -> Result<SuperposeReport> {
    let (a, b) = paired_traces(model, reference)?;
    let sup =
        kabsch(&a, &b).map_err(|e| StructPredictError::invalid("superpose", e.to_string()))?;
    // Apply the transform to the model trace, then measure distances.
    let moved = sup.transform.apply_all(&a);
    let n = moved.len();
    let deviations: Vec<f64> = moved.iter().zip(&b).map(|(p, q)| (p - q).norm()).collect();

    let fraction_within = |cutoff: f64| -> f64 {
        deviations.iter().filter(|&&d| d <= cutoff).count() as f64 / n as f64
    };

    // GDT-TS: cutoffs 1, 2, 4, 8 Å.
    let gdt_ts = 100.0
        * (fraction_within(1.0)
            + fraction_within(2.0)
            + fraction_within(4.0)
            + fraction_within(8.0))
        / 4.0;
    // GDT-HA: stricter cutoffs 0.5, 1, 2, 4 Å.
    let gdt_ha = 100.0
        * (fraction_within(0.5)
            + fraction_within(1.0)
            + fraction_within(2.0)
            + fraction_within(4.0))
        / 4.0;

    Ok(SuperposeReport {
        rmsd: sup.rmsd,
        gdt_ts,
        gdt_ha,
        n,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Rotation3, Vector3};

    fn helix_model(n: usize) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            let t = i as f64 * 1.745; // ~100° per residue
            r.ca = Some(Point3::new(2.3 * t.cos(), 2.3 * t.sin(), 1.5 * i as f64));
        }
        m
    }

    #[test]
    fn identical_model_has_zero_rmsd_full_gdt() {
        let m = helix_model(20);
        let report = model_vs_reference(&m, &m).expect("superpose");
        assert!(report.rmsd < 1e-6, "rmsd {}", report.rmsd);
        assert!((report.gdt_ts - 100.0).abs() < 1e-6);
        assert!((report.gdt_ha - 100.0).abs() < 1e-6);
    }

    #[test]
    fn rigid_motion_does_not_change_rmsd() {
        let m = helix_model(20);
        let mut moved = m.clone();
        let rot = Rotation3::from_euler_angles(0.5, 1.1, -0.3);
        let trans = Vector3::new(7.0, -4.0, 12.0);
        for r in moved.residues.iter_mut() {
            if let Some(ca) = r.ca.as_mut() {
                *ca = Point3::from(rot * ca.coords + trans);
            }
        }
        // Kabsch undoes the rigid motion → RMSD ≈ 0.
        let rmsd = ca_rmsd_superposed(&moved, &m).expect("rmsd");
        assert!(rmsd < 1e-6, "rmsd after rigid motion {rmsd}");
    }

    #[test]
    fn perturbed_model_loses_gdt() {
        let m = helix_model(20);
        let mut noisy = m.clone();
        for (i, r) in noisy.residues.iter_mut().enumerate() {
            if let Some(ca) = r.ca.as_mut() {
                ca.x += if i % 2 == 0 { 3.0 } else { -3.0 };
            }
        }
        let report = model_vs_reference(&noisy, &m).expect("superpose");
        assert!(report.rmsd > 0.5, "rmsd {}", report.rmsd);
        assert!(report.gdt_ts < 100.0, "gdt {}", report.gdt_ts);
        // GDT-HA is the stricter metric → never above GDT-TS here.
        assert!(report.gdt_ha <= report.gdt_ts + 1e-9);
    }

    #[test]
    fn too_short_rejected() {
        let m = helix_model(2);
        assert!(model_vs_reference(&m, &m).is_err());
    }
}
