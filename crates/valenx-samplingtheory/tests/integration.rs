//! Cross-module integration tests for `valenx-samplingtheory`.
//!
//! These exercise the public surface end to end: a sampling-chain
//! sanity check that ties the Nyquist, aliasing and quantization
//! modules together, and a JSON round-trip of the serializable analysis
//! bundle (the only place the crate touches serde at the boundary).

use valenx_samplingtheory::aliasing::{alias_frequency, is_aliased};
use valenx_samplingtheory::nyquist::{nyquist_frequency, satisfies_nyquist};
use valenx_samplingtheory::quantization::{ideal_snr_db, quant_step, QuantizationAnalysis};

/// Tolerance for floating-point comparisons.
const EPS: f64 = 1e-9;

#[test]
fn full_sampling_chain_is_consistent() {
    // A 12 kHz tone sampled at 8 kHz. The Nyquist frequency is 4 kHz,
    // the tone is above it (so it aliases), and it folds to
    // |12000 - round(1.5)*8000| = |12000 - 16000| = 4000 Hz.
    let fs = 8_000.0;
    let tone = 12_000.0;

    let fnyq = nyquist_frequency(fs).unwrap();
    assert!((fnyq - 4_000.0).abs() < EPS, "got {fnyq}");

    // The tone exceeds the Nyquist frequency, so it cannot be captured.
    assert!(!satisfies_nyquist(fs, tone).unwrap());
    assert!(is_aliased(tone, fs).unwrap());

    let a = alias_frequency(tone, fs).unwrap();
    assert!((a - 4_000.0).abs() < EPS, "got {a}");
    // The alias always lands within the baseband [0, fs/2].
    assert!(a <= fnyq + EPS, "alias {a} exceeds Nyquist {fnyq}");
}

#[test]
fn quantization_figures_agree_across_helpers() {
    // The analysis bundle must agree with the standalone helpers.
    let range = 5.0;
    let bits = 14;

    let bundle = QuantizationAnalysis::new(range, bits).unwrap();
    assert!((bundle.step - quant_step(range, bits).unwrap()).abs() < EPS);
    assert!((bundle.ideal_snr_db - ideal_snr_db(bits).unwrap()).abs() < EPS);
    assert_eq!(bundle.levels, 1u64 << bits);
}

#[test]
fn analysis_bundle_survives_json_round_trip() {
    let original = QuantizationAnalysis::new(3.3, 12).unwrap();

    let json = serde_json::to_string(&original).expect("serialize");
    let restored: QuantizationAnalysis = serde_json::from_str(&json).expect("deserialize");

    assert_eq!(original, restored);
    // Spot-check a field survived the trip with full precision.
    assert!(
        (restored.step - 3.3 / 4096.0).abs() < EPS,
        "got {}",
        restored.step
    );
}
