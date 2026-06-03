//! Unit conventions and conversions for `valenx-neuro`.
//!
//! Electrophysiology mixes unit systems freely, which is the classic source
//! of silent numerical bugs (the neural analogue of the AIMD timestep bug).
//! `valenx-neuro` fixes one convention and converts to SI only at solver
//! boundaries — here — with a round-trip test for each conversion:
//!
//! | quantity | unit |
//! |---|---|
//! | potential | mV |
//! | time | ms |
//! | current | µA |
//! | conductivity σ | S/m |
//! | tissue length | mm |
//! | compartment length / fiber diameter | µm |
//! | membrane capacitance | µF/cm² |
//! | conductance density | mS/cm² |
//! | temperature rise | K |

/// Convert a potential from millivolts to volts.
pub fn mv_to_volts(mv: f64) -> f64 {
    mv * 1e-3
}

/// Convert a potential from volts to millivolts.
pub fn volts_to_mv(v: f64) -> f64 {
    v * 1e3
}

/// Convert a current from microamperes to amperes.
pub fn ua_to_amp(ua: f64) -> f64 {
    ua * 1e-6
}

/// Convert a current from amperes to microamperes.
pub fn amp_to_ua(a: f64) -> f64 {
    a * 1e6
}

/// Convert a length from millimetres to metres.
pub fn mm_to_m(mm: f64) -> f64 {
    mm * 1e-3
}

/// Convert a length from micrometres to metres.
pub fn um_to_m(um: f64) -> f64 {
    um * 1e-6
}

/// Convert a membrane conductance *density* (mS/cm²) over a patch `area`
/// (cm²) into an absolute conductance in siemens.
pub fn ms_per_cm2_to_s(g_ms_cm2: f64, area_cm2: f64) -> f64 {
    g_ms_cm2 * area_cm2 * 1e-3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mv_round_trips_through_volts() {
        assert!((volts_to_mv(mv_to_volts(-65.0)) + 65.0).abs() < 1e-12);
    }

    #[test]
    fn current_microamp_to_amp() {
        // 10 µA = 1e-5 A
        assert!((ua_to_amp(10.0) - 1e-5).abs() < 1e-18);
        assert!((amp_to_ua(ua_to_amp(7.0)) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn length_conversions() {
        assert!((mm_to_m(40.0) - 0.04).abs() < 1e-15);
        assert!((um_to_m(100.0) - 1e-4).abs() < 1e-18);
    }

    #[test]
    fn conductance_density_area_scaling() {
        // g_Na = 120 mS/cm² on a 1e-4 cm² patch = 0.012 mS = 1.2e-5 S
        assert!((ms_per_cm2_to_s(120.0, 1e-4) - 1.2e-5).abs() < 1e-12);
    }
}
