//! Net charge versus pH and the isoelectric point (Henderson-Hasselbalch).

use crate::aa::{negative_sidechain_pka, positive_sidechain_pka, PKA_C_TERM, PKA_N_TERM};
use crate::error::DevelopabilityError;

fn validate(seq: &str) -> Result<(), DevelopabilityError> {
    if seq.is_empty() {
        return Err(DevelopabilityError::EmptySequence);
    }
    if let Some((pos, residue)) = crate::aa::first_invalid(seq) {
        return Err(DevelopabilityError::InvalidResidue { residue, pos });
    }
    Ok(())
}

/// Fraction of a basic (positively-ionizable) group that is protonated (`+1`).
fn frac_protonated(pka: f64, ph: f64) -> f64 {
    1.0 / (1.0 + 10f64.powf(ph - pka))
}

/// Fraction of an acidic (negatively-ionizable) group that is deprotonated
/// (`-1`).
fn frac_deprotonated(pka: f64, ph: f64) -> f64 {
    1.0 / (1.0 + 10f64.powf(pka - ph))
}

/// Net charge of the sequence at a given `ph`, summing the protonated basic
/// groups (N-terminus, His/Lys/Arg) as `+` and the deprotonated acidic groups
/// (C-terminus, Cys/Asp/Glu/Tyr) as `-`.
pub fn net_charge_at_ph(seq: &str, ph: f64) -> Result<f64, DevelopabilityError> {
    validate(seq)?;
    if !ph.is_finite() {
        return Err(DevelopabilityError::NonFinite { what: "ph" });
    }
    let mut positive = frac_protonated(PKA_N_TERM, ph);
    let mut negative = frac_deprotonated(PKA_C_TERM, ph);
    for b in seq.bytes() {
        if let Some(pk) = positive_sidechain_pka(b) {
            positive += frac_protonated(pk, ph);
        }
        if let Some(pk) = negative_sidechain_pka(b) {
            negative += frac_deprotonated(pk, ph);
        }
    }
    Ok(positive - negative)
}

/// The isoelectric point: the pH at which [`net_charge_at_ph`] is zero, found by
/// bisection on `[0, 14]` (net charge is monotone decreasing in pH).
pub fn isoelectric_point(seq: &str) -> Result<f64, DevelopabilityError> {
    validate(seq)?;
    let (mut lo, mut hi) = (0.0_f64, 14.0_f64);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if net_charge_at_ph(seq, mid)? > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lysine_is_positive_aspartate_is_negative_at_neutral_ph() {
        assert!(net_charge_at_ph("K", 7.0).unwrap() > 0.0);
        assert!(net_charge_at_ph("D", 7.0).unwrap() < 0.0);
    }

    #[test]
    fn net_charge_decreases_with_ph() {
        let seq = "KDEYHR";
        assert!(net_charge_at_ph(seq, 3.0).unwrap() > net_charge_at_ph(seq, 11.0).unwrap());
    }

    #[test]
    fn pi_bands_make_sense() {
        // termini-only -> pI near (8.6 + 3.6)/2 = 6.1
        let neutral = isoelectric_point("AAAA").unwrap();
        assert!((5.0..7.0).contains(&neutral), "pI = {neutral}");
        // poly-Lys is strongly basic, poly-Glu strongly acidic
        assert!(isoelectric_point("KKKK").unwrap() > 9.5);
        assert!(isoelectric_point("EEEE").unwrap() < 4.5);
    }

    #[test]
    fn net_charge_is_zero_at_pi() {
        let seq = "KDEYHRAA";
        let pi = isoelectric_point(seq).unwrap();
        assert!(net_charge_at_ph(seq, pi).unwrap().abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            net_charge_at_ph("", 7.0).unwrap_err().code(),
            "empty_sequence"
        );
        assert_eq!(
            net_charge_at_ph("AZ", 7.0).unwrap_err().code(),
            "invalid_residue"
        );
        assert_eq!(
            net_charge_at_ph("A", f64::NAN).unwrap_err().code(),
            "non_finite"
        );
    }
}
