//! Ground-truth analytic tests for the op-amp models.
//!
//! Each test checks a closed-form result against a hand-computed value
//! using an absolute-difference tolerance — never an exact float
//! equality — so the arithmetic is verified without being brittle.

use valenx_opamp::{Gbw, Inverting, NonInverting, OpAmpError, SummingAmplifier, VoltageFollower};

/// Absolute tolerance for the dimensionless / volt-scale assertions.
const EPS: f64 = 1e-9;
/// Looser tolerance for hertz-scale assertions (values up to ~1e6).
const EPS_HZ: f64 = 1e-3;

fn close(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

// ---------------------------------------------------------------------
// Inverting amplifier: G = -Rf / Rin
// ---------------------------------------------------------------------

#[test]
fn inverting_gain_is_negative_ratio() {
    // Rf = 100k, Rin = 10k -> G = -10.
    let stage = Inverting::new(10_000.0, 100_000.0).expect("valid resistors");
    let g = stage.gain();
    assert!(
        close(g, -10.0, EPS),
        "inverting gain should be -Rf/Rin = -10, got {g}"
    );
    assert!(g < 0.0, "inverting gain must be negative, got {g}");
}

#[test]
fn inverting_unity_when_resistors_equal() {
    // Equal resistors give a precise -1 inverter.
    let stage = Inverting::new(4_700.0, 4_700.0).expect("valid resistors");
    let g = stage.gain();
    assert!(
        close(g, -1.0, EPS),
        "equal-R inverter should be -1, got {g}"
    );
    let mag = stage.gain_magnitude();
    assert!(close(mag, 1.0, EPS), "magnitude should be 1, got {mag}");
}

#[test]
fn inverting_output_scales_input() {
    let stage = Inverting::new(1_000.0, 22_000.0).expect("valid resistors");
    // G = -22; Vin = 0.1 V -> Vout = -2.2 V.
    let v_out = stage.output(0.1);
    assert!(
        close(v_out, -2.2, EPS),
        "Vout should be G*Vin = -2.2, got {v_out}"
    );
}

// ---------------------------------------------------------------------
// Non-inverting amplifier: G = 1 + Rf / Rin, always >= 1
// ---------------------------------------------------------------------

#[test]
fn non_inverting_gain_is_one_plus_ratio() {
    // Rf = 90k, Rin = 10k -> G = 1 + 9 = 10.
    let stage = NonInverting::new(10_000.0, 90_000.0).expect("valid resistors");
    let g = stage.gain();
    assert!(
        close(g, 10.0, EPS),
        "non-inverting gain should be 1 + Rf/Rin = 10, got {g}"
    );
}

#[test]
fn non_inverting_gain_never_below_one() {
    // Sweep a range of resistor ratios; every gain must be >= 1.
    for r_in in [100.0, 1_000.0, 47_000.0, 1_000_000.0] {
        for r_f in [1.0, 100.0, 10_000.0, 5_000_000.0] {
            let stage = NonInverting::new(r_in, r_f).expect("valid resistors");
            let g = stage.gain();
            assert!(
                g >= 1.0,
                "non-inverting gain must be >= 1 for r_in={r_in}, r_f={r_f}, got {g}"
            );
        }
    }
}

#[test]
fn non_inverting_approaches_unity_for_small_feedback() {
    // Tiny Rf relative to Rin -> gain just above 1.
    let stage = NonInverting::new(1_000_000.0, 1.0).expect("valid resistors");
    let g = stage.gain();
    // 1 + 1/1e6 = 1.000001.
    assert!(
        close(g, 1.000_001, EPS),
        "gain should approach 1 from above, got {g}"
    );
    assert!(g > 1.0, "gain must stay strictly above 1, got {g}");
}

// ---------------------------------------------------------------------
// Voltage follower: G = 1 exactly
// ---------------------------------------------------------------------

#[test]
fn voltage_follower_is_unity() {
    let buf = VoltageFollower::new();
    let g = buf.gain();
    assert!(close(g, 1.0, EPS), "follower gain should be 1, got {g}");
    let v_out = buf.output(3.3);
    assert!(
        close(v_out, 3.3, EPS),
        "follower output should equal input, got {v_out}"
    );
}

#[test]
fn voltage_follower_matches_zero_feedback_non_inverting() {
    // A non-inverting stage with vanishing Rf tends to the follower.
    let buf = VoltageFollower::new();
    let near = NonInverting::new(1.0, 1e-12).expect("valid resistors");
    let diff = (near.gain() - buf.gain()).abs();
    assert!(
        diff < 1e-6,
        "near-zero feedback should approach unity follower, diff {diff}"
    );
}

// ---------------------------------------------------------------------
// Summing amplifier: Vout = -Rf * Σ(Vi / Ri)
// ---------------------------------------------------------------------

#[test]
fn summing_two_equal_branches() {
    // Rf = Ri = 10k for both branches -> Vout = -(V1 + V2).
    let amp =
        SummingAmplifier::new(10_000.0, [(1.0, 10_000.0), (2.0, 10_000.0)]).expect("valid amp");
    let v_out = amp.output();
    assert!(
        close(v_out, -3.0, EPS),
        "equal-weight sum of 1 V + 2 V should be -3 V, got {v_out}"
    );
}

#[test]
fn summing_weighted_branches() {
    // Rf = 20k. Branch A: 1 V / 10k -> weight 2. Branch B: 1 V / 5k -> weight 4.
    // Vout = -20k * (1/10k + 1/5k) = -(2 + 4) = -6 V.
    let amp =
        SummingAmplifier::new(20_000.0, [(1.0, 10_000.0), (1.0, 5_000.0)]).expect("valid amp");
    let v_out = amp.output();
    assert!(
        close(v_out, -6.0, EPS),
        "weighted sum should be -6 V, got {v_out}"
    );
}

#[test]
fn summing_single_branch_reduces_to_inverting() {
    // One branch is just an inverting amplifier: -Rf/Ri * Vi.
    let r_in = 2_000.0;
    let r_f = 8_000.0;
    let v_in = 0.5;
    let amp = SummingAmplifier::new(r_f, [(v_in, r_in)]).expect("valid amp");
    let inv = Inverting::new(r_in, r_f).expect("valid resistors");
    let diff = (amp.output() - inv.output(v_in)).abs();
    assert!(
        diff < EPS,
        "single-branch summer must equal inverting stage, diff {diff}"
    );
}

#[test]
fn summing_branch_gain_matches_definition() {
    let amp =
        SummingAmplifier::new(30_000.0, [(0.0, 10_000.0), (0.0, 15_000.0)]).expect("valid amp");
    // Branch 0: -30k/10k = -3. Branch 1: -30k/15k = -2.
    let g0 = amp.branch_gain(0).expect("branch 0 exists");
    let g1 = amp.branch_gain(1).expect("branch 1 exists");
    assert!(close(g0, -3.0, EPS), "branch 0 gain should be -3, got {g0}");
    assert!(close(g1, -2.0, EPS), "branch 1 gain should be -2, got {g1}");
    assert!(
        amp.branch_gain(2).is_none(),
        "out-of-range branch must be None"
    );
}

// ---------------------------------------------------------------------
// Gain-bandwidth relations
// ---------------------------------------------------------------------

#[test]
fn gbw_from_gain_bandwidth_product() {
    // A part with closed-loop gain 10 at 100 kHz has GBW = 1 MHz.
    let gbw = Gbw::from_gain_bandwidth(10.0, 100_000.0).expect("valid pair");
    assert!(
        close(gbw.hz(), 1_000_000.0, EPS_HZ),
        "GBW = gain*bandwidth should be 1 MHz, got {hz}",
        hz = gbw.hz()
    );
}

#[test]
fn closed_loop_bandwidth_is_gbw_over_gain() {
    // 1 MHz GBW at gain 100 -> 10 kHz bandwidth.
    let gbw = Gbw::new(1_000_000.0).expect("valid gbw");
    let bw = gbw.closed_loop_bandwidth(100.0).expect("valid gain");
    assert!(
        close(bw, 10_000.0, EPS_HZ),
        "closed-loop BW = GBW/gain should be 10 kHz, got {bw}"
    );
}

#[test]
fn unity_gain_bandwidth_equals_gbw() {
    let gbw = Gbw::new(3_000_000.0).expect("valid gbw");
    let ugbw = gbw.unity_gain_bandwidth();
    assert!(
        close(ugbw, gbw.hz(), EPS_HZ),
        "unity-gain BW must equal GBW, got {ugbw} vs {hz}",
        hz = gbw.hz()
    );
    // And it must equal the bandwidth at gain = 1.
    let at_unity = gbw.closed_loop_bandwidth(1.0).expect("valid gain");
    assert!(
        close(at_unity, gbw.hz(), EPS_HZ),
        "closed-loop BW at gain 1 must equal GBW, got {at_unity}"
    );
}

#[test]
fn gbw_round_trips_through_gain() {
    // GBW derived from (gain, bw) must reproduce that bw at that gain.
    let gain = 47.0;
    let bw_in = 21_276.6;
    let gbw = Gbw::from_gain_bandwidth(gain, bw_in).expect("valid pair");
    let bw_out = gbw.closed_loop_bandwidth(gain).expect("valid gain");
    let diff = (bw_in - bw_out).abs();
    assert!(diff < EPS_HZ, "round-trip bandwidth mismatch, diff {diff}");
}

// ---------------------------------------------------------------------
// Validation / error paths
// ---------------------------------------------------------------------

#[test]
fn rejects_non_positive_resistance() {
    let err = Inverting::new(0.0, 100.0).expect_err("zero r_in must fail");
    assert!(
        matches!(err, OpAmpError::NonPositive { name: "r_in", .. }),
        "expected NonPositive for r_in, got {err:?}"
    );
    assert_eq!(err.code(), "opamp.non-positive", "stable code mismatch");

    let err2 = NonInverting::new(100.0, -5.0).expect_err("negative r_f must fail");
    assert!(
        matches!(err2, OpAmpError::NonPositive { name: "r_f", .. }),
        "expected NonPositive for r_f, got {err2:?}"
    );
}

#[test]
fn rejects_non_finite_input() {
    let err = Inverting::new(f64::NAN, 100.0).expect_err("NaN r_in must fail");
    assert!(
        matches!(err, OpAmpError::NotFinite { name: "r_in", .. }),
        "expected NotFinite for r_in, got {err:?}"
    );
    assert_eq!(err.code(), "opamp.not-finite", "stable code mismatch");

    let err2 = Gbw::new(f64::INFINITY).expect_err("infinite GBW must fail");
    assert!(
        matches!(err2, OpAmpError::NotFinite { name: "gbw_hz", .. }),
        "expected NotFinite for gbw_hz, got {err2:?}"
    );
}

#[test]
fn rejects_empty_summing_inputs() {
    let empty: [(f64, f64); 0] = [];
    let err = SummingAmplifier::new(10_000.0, empty).expect_err("no inputs must fail");
    assert!(
        matches!(err, OpAmpError::NoInputs),
        "expected NoInputs, got {err:?}"
    );
    assert_eq!(err.code(), "opamp.no-inputs", "stable code mismatch");
}

#[test]
fn summing_rejects_bad_branch_resistance() {
    let err =
        SummingAmplifier::new(10_000.0, [(1.0, 5_000.0), (2.0, 0.0)]).expect_err("zero branch R");
    assert!(
        matches!(err, OpAmpError::NonPositive { name: "r", .. }),
        "expected NonPositive for branch r, got {err:?}"
    );
}
