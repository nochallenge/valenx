//! `valenx-quantum` — quantum circuit simulation: state-vector qubits,
//! gates, and projective measurement.
//!
//! This is an in-house, dependency-light state-vector simulator. An
//! `n`-qubit register is stored as `2^n` complex amplitudes
//! ([`StateVector`]); gates are applied directly to that vector and
//! [`StateVector::probabilities`] returns the Born-rule distribution.
//!
//! # Conventions
//!
//! * Qubit `0` is the **least-significant** bit of a basis-state index, so
//!   the index `b` of computational basis state `|q_{n-1} … q_1 q_0>` has
//!   bit `q_k = (b >> k) & 1`.
//! * Amplitudes are normalised so that `sum |amplitude|^2 == 1`; every gate
//!   preserves this norm exactly (up to floating-point round-off).
//! * Out-of-range qubit indices fail loud with a [`QuantumError`] rather
//!   than panicking or silently producing a wrong answer.
//!
//! # Example — a Bell state
//!
//! ```
//! use valenx_quantum::{Circuit, Gate};
//!
//! let probs = Circuit::new(2)
//!     .h(0)
//!     .cnot(0, 1)
//!     .run()
//!     .unwrap()
//!     .probabilities();
//!
//! // P(00) == P(11) == 0.5, the cross terms vanish.
//! assert!((probs[0b00] - 0.5).abs() < 1e-12);
//! assert!((probs[0b11] - 0.5).abs() < 1e-12);
//! assert!(probs[0b01].abs() < 1e-12);
//! assert!(probs[0b10].abs() < 1e-12);
//! let _ = Gate::H; // gates are also exposed directly.
//! ```

use num_complex::Complex64;
use thiserror::Error;

/// A complex probability amplitude (64-bit real and imaginary parts).
pub type Amp = Complex64;

/// Errors returned by the simulator instead of panicking.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum QuantumError {
    /// A register with zero qubits was requested. At least one qubit is
    /// required because the state vector would otherwise be empty.
    #[error("a quantum register needs at least one qubit")]
    ZeroQubits,

    /// `n` qubits would need `2^n` amplitudes, which overflows the address
    /// space (or the supported maximum). Reported with the offending count.
    #[error("{0} qubits is too large to simulate as a dense state vector")]
    TooManyQubits(usize),

    /// A gate referenced a qubit index outside `0..num_qubits`.
    #[error("qubit index {index} is out of range for a {num_qubits}-qubit register")]
    QubitOutOfRange {
        /// The offending qubit index.
        index: usize,
        /// The number of qubits in the register.
        num_qubits: usize,
    },

    /// A two-qubit gate was given the same index for both operands.
    #[error("a two-qubit gate needs two distinct qubits, got {0} twice")]
    DuplicateQubit(usize),
}

/// The largest register this dense simulator will allocate.
///
/// `2^26` amplitudes is already a gigabyte of `Complex64`; refusing larger
/// registers keeps a mistaken qubit count from exhausting memory.
pub const MAX_QUBITS: usize = 26;

/// A pure quantum state of `n` qubits, stored as a dense vector of `2^n`
/// complex amplitudes.
#[derive(Debug, Clone, PartialEq)]
pub struct StateVector {
    num_qubits: usize,
    amps: Vec<Amp>,
}

impl StateVector {
    /// Create the all-zero computational basis state `|0…0>` of `n` qubits.
    ///
    /// # Errors
    ///
    /// Returns [`QuantumError::ZeroQubits`] if `n == 0`, or
    /// [`QuantumError::TooManyQubits`] if `n > `[`MAX_QUBITS`].
    pub fn new(n: usize) -> Result<Self, QuantumError> {
        if n == 0 {
            return Err(QuantumError::ZeroQubits);
        }
        if n > MAX_QUBITS {
            return Err(QuantumError::TooManyQubits(n));
        }
        let dim = 1usize << n;
        let mut amps = vec![Amp::new(0.0, 0.0); dim];
        amps[0] = Amp::new(1.0, 0.0);
        Ok(Self {
            num_qubits: n,
            amps,
        })
    }

    /// The number of qubits in this register.
    #[must_use]
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    /// The dimension of the state vector, i.e. `2^n`.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.amps.len()
    }

    /// The raw amplitude vector, indexed by computational basis state.
    #[must_use]
    pub fn amplitudes(&self) -> &[Amp] {
        &self.amps
    }

    /// The Born-rule probability distribution over computational basis
    /// states: entry `b` is `|amplitude_b|^2`.
    ///
    /// The returned vector sums to `1` up to floating-point round-off.
    #[must_use]
    pub fn probabilities(&self) -> Vec<f64> {
        self.amps
            .iter()
            .map(num_complex::Complex::norm_sqr)
            .collect()
    }

    /// The total probability `sum |amplitude|^2`. Equals `1` for a valid
    /// normalised state; handy as a test/diagnostic invariant.
    #[must_use]
    pub fn norm_sqr(&self) -> f64 {
        self.amps.iter().map(num_complex::Complex::norm_sqr).sum()
    }

    /// Validate that `q` is a legal qubit index for this register.
    fn check(&self, q: usize) -> Result<(), QuantumError> {
        if q >= self.num_qubits {
            return Err(QuantumError::QubitOutOfRange {
                index: q,
                num_qubits: self.num_qubits,
            });
        }
        Ok(())
    }

    /// Apply a single-qubit gate given by its `2x2` unitary `m` (row-major:
    /// `[[m00, m01], [m10, m11]]`) to qubit `q`.
    ///
    /// # Errors
    ///
    /// Returns [`QuantumError::QubitOutOfRange`] if `q` is not a valid qubit.
    pub fn apply_1q(&mut self, q: usize, m: [[Amp; 2]; 2]) -> Result<(), QuantumError> {
        self.check(q)?;
        let bit = 1usize << q;
        // Iterate over the half of basis states with bit `q` cleared; each
        // pairs with its bit-`q`-set partner. This visits every amplitude
        // exactly once.
        for base in 0..self.amps.len() {
            if base & bit != 0 {
                continue;
            }
            let i0 = base;
            let i1 = base | bit;
            let a0 = self.amps[i0];
            let a1 = self.amps[i1];
            self.amps[i0] = m[0][0] * a0 + m[0][1] * a1;
            self.amps[i1] = m[1][0] * a0 + m[1][1] * a1;
        }
        Ok(())
    }

    /// Apply a [`Gate`] to qubit `q`.
    ///
    /// # Errors
    ///
    /// Propagates [`QuantumError::QubitOutOfRange`] from [`Self::apply_1q`].
    pub fn apply(&mut self, gate: Gate, q: usize) -> Result<(), QuantumError> {
        self.apply_1q(q, gate.matrix())
    }

    /// Apply a controlled-NOT: flip `target` iff `control` is `|1>`.
    ///
    /// # Errors
    ///
    /// Returns [`QuantumError::DuplicateQubit`] if `control == target`, or
    /// [`QuantumError::QubitOutOfRange`] for an invalid index.
    pub fn cnot(&mut self, control: usize, target: usize) -> Result<(), QuantumError> {
        self.check(control)?;
        self.check(target)?;
        if control == target {
            return Err(QuantumError::DuplicateQubit(control));
        }
        let cbit = 1usize << control;
        let tbit = 1usize << target;
        for base in 0..self.amps.len() {
            // Visit each (control=1, target=0) state once and swap with its
            // target=1 partner.
            if base & cbit != 0 && base & tbit == 0 {
                let partner = base | tbit;
                self.amps.swap(base, partner);
            }
        }
        Ok(())
    }

    /// Apply a controlled-Z: multiply the amplitude by `-1` iff both
    /// `control` and `target` are `|1>`. (CZ is symmetric in its operands.)
    ///
    /// # Errors
    ///
    /// Returns [`QuantumError::DuplicateQubit`] if the two indices coincide,
    /// or [`QuantumError::QubitOutOfRange`] for an invalid index.
    pub fn cz(&mut self, control: usize, target: usize) -> Result<(), QuantumError> {
        self.check(control)?;
        self.check(target)?;
        if control == target {
            return Err(QuantumError::DuplicateQubit(control));
        }
        let mask = (1usize << control) | (1usize << target);
        for (b, amp) in self.amps.iter_mut().enumerate() {
            if b & mask == mask {
                *amp = -*amp;
            }
        }
        Ok(())
    }

    /// Sample one computational basis outcome (full-register measurement)
    /// using the supplied seeded PRNG, **without** collapsing the state.
    ///
    /// Sampling is deterministic for a given seed, which keeps tests
    /// reproducible. Returns the measured basis-state index in `0..2^n`.
    pub fn sample(&self, rng: &mut SeededRng) -> usize {
        let r = rng.next_f64();
        let mut acc = 0.0;
        let last = self.amps.len() - 1;
        for (b, amp) in self.amps.iter().enumerate() {
            acc += amp.norm_sqr();
            if r < acc {
                return b;
            }
        }
        // Floating-point round-off can leave `r` just past the final
        // cumulative bound; attribute it to the last basis state.
        last
    }

    /// Draw `shots` independent measurement samples with the supplied seeded
    /// PRNG and return a histogram of length `2^n` counting each outcome.
    ///
    /// The state is not collapsed; each shot samples the same distribution.
    #[must_use]
    pub fn measure_counts(&self, shots: usize, rng: &mut SeededRng) -> Vec<usize> {
        let mut counts = vec![0usize; self.amps.len()];
        for _ in 0..shots {
            counts[self.sample(rng)] += 1;
        }
        counts
    }
}

/// The built-in single-qubit gates.
///
/// Each variant maps to a fixed `2x2` unitary via [`Gate::matrix`]. The
/// rotation gates carry their angle in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Gate {
    /// Hadamard — maps `|0>` to the equal superposition `(|0>+|1>)/√2`.
    H,
    /// Pauli-X (bit flip / NOT).
    X,
    /// Pauli-Y.
    Y,
    /// Pauli-Z (phase flip).
    Z,
    /// Phase gate `S = diag(1, i)` (a quarter turn about Z).
    S,
    /// `T = diag(1, e^{iπ/4})` (an eighth turn about Z).
    T,
    /// Rotation about the X axis by `θ` radians.
    Rx(f64),
    /// Rotation about the Y axis by `θ` radians.
    Ry(f64),
    /// Rotation about the Z axis by `θ` radians.
    Rz(f64),
}

impl Gate {
    /// The `2x2` unitary matrix for this gate, row-major
    /// `[[m00, m01], [m10, m11]]`.
    #[must_use]
    pub fn matrix(self) -> [[Amp; 2]; 2] {
        let re = |x: f64| Amp::new(x, 0.0);
        let i = Amp::new(0.0, 1.0);
        match self {
            Gate::H => {
                let s = std::f64::consts::FRAC_1_SQRT_2;
                [[re(s), re(s)], [re(s), re(-s)]]
            }
            Gate::X => [[re(0.0), re(1.0)], [re(1.0), re(0.0)]],
            Gate::Y => [[re(0.0), -i], [i, re(0.0)]],
            Gate::Z => [[re(1.0), re(0.0)], [re(0.0), re(-1.0)]],
            Gate::S => [[re(1.0), re(0.0)], [re(0.0), i]],
            Gate::T => {
                let phase = (i * std::f64::consts::FRAC_PI_4).exp();
                [[re(1.0), re(0.0)], [re(0.0), phase]]
            }
            Gate::Rx(theta) => {
                let c = re((theta / 2.0).cos());
                let s = (theta / 2.0).sin();
                let ms = Amp::new(0.0, -s);
                [[c, ms], [ms, c]]
            }
            Gate::Ry(theta) => {
                let c = re((theta / 2.0).cos());
                let s = re((theta / 2.0).sin());
                [[c, -s], [s, c]]
            }
            Gate::Rz(theta) => {
                let neg = (-i * (theta / 2.0)).exp();
                let pos = (i * (theta / 2.0)).exp();
                [[neg, re(0.0)], [re(0.0), pos]]
            }
        }
    }
}

/// A single instruction in a [`Circuit`].
#[derive(Debug, Clone, Copy, PartialEq)]
enum Op {
    Single(Gate, usize),
    Cnot(usize, usize),
    Cz(usize, usize),
}

/// A small builder that records a sequence of gates and applies them to a
/// fresh `|0…0>` register.
///
/// Builder methods are chainable; [`Circuit::run`] executes the recorded
/// program and returns the resulting [`StateVector`].
#[derive(Debug, Clone)]
pub struct Circuit {
    num_qubits: usize,
    ops: Vec<Op>,
}

impl Circuit {
    /// Start an empty circuit on `n` qubits.
    ///
    /// The qubit count is validated lazily when [`Circuit::run`] is called.
    #[must_use]
    pub fn new(n: usize) -> Self {
        Self {
            num_qubits: n,
            ops: Vec::new(),
        }
    }

    /// The number of qubits this circuit targets.
    #[must_use]
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    /// Append an arbitrary single-qubit [`Gate`] on qubit `q`.
    #[must_use]
    pub fn gate(mut self, gate: Gate, q: usize) -> Self {
        self.ops.push(Op::Single(gate, q));
        self
    }

    /// Append a Hadamard on qubit `q`.
    #[must_use]
    pub fn h(self, q: usize) -> Self {
        self.gate(Gate::H, q)
    }

    /// Append a Pauli-X on qubit `q`.
    #[must_use]
    pub fn x(self, q: usize) -> Self {
        self.gate(Gate::X, q)
    }

    /// Append a Pauli-Y on qubit `q`.
    #[must_use]
    pub fn y(self, q: usize) -> Self {
        self.gate(Gate::Y, q)
    }

    /// Append a Pauli-Z on qubit `q`.
    #[must_use]
    pub fn z(self, q: usize) -> Self {
        self.gate(Gate::Z, q)
    }

    /// Append a phase gate `S` on qubit `q`.
    #[must_use]
    pub fn s(self, q: usize) -> Self {
        self.gate(Gate::S, q)
    }

    /// Append a `T` gate on qubit `q`.
    #[must_use]
    pub fn t(self, q: usize) -> Self {
        self.gate(Gate::T, q)
    }

    /// Append an X-rotation by `theta` radians on qubit `q`.
    #[must_use]
    pub fn rx(self, q: usize, theta: f64) -> Self {
        self.gate(Gate::Rx(theta), q)
    }

    /// Append a Y-rotation by `theta` radians on qubit `q`.
    #[must_use]
    pub fn ry(self, q: usize, theta: f64) -> Self {
        self.gate(Gate::Ry(theta), q)
    }

    /// Append a Z-rotation by `theta` radians on qubit `q`.
    #[must_use]
    pub fn rz(self, q: usize, theta: f64) -> Self {
        self.gate(Gate::Rz(theta), q)
    }

    /// Append a CNOT with the given `control` and `target` qubits.
    #[must_use]
    pub fn cnot(mut self, control: usize, target: usize) -> Self {
        self.ops.push(Op::Cnot(control, target));
        self
    }

    /// Append a CZ with the given `control` and `target` qubits.
    #[must_use]
    pub fn cz(mut self, control: usize, target: usize) -> Self {
        self.ops.push(Op::Cz(control, target));
        self
    }

    /// Execute the recorded program on a fresh `|0…0>` register.
    ///
    /// # Errors
    ///
    /// Propagates any [`QuantumError`] raised while building the state or
    /// applying a gate (bad qubit count or out-of-range / duplicate index).
    pub fn run(&self) -> Result<StateVector, QuantumError> {
        let mut sv = StateVector::new(self.num_qubits)?;
        for op in &self.ops {
            match *op {
                Op::Single(g, q) => sv.apply(g, q)?,
                Op::Cnot(c, t) => sv.cnot(c, t)?,
                Op::Cz(c, t) => sv.cz(c, t)?,
            }
        }
        Ok(sv)
    }
}

/// A tiny deterministic PRNG (`SplitMix64`) used for measurement sampling.
///
/// Seeding it with a fixed value makes [`StateVector::sample`] and
/// [`StateVector::measure_counts`] reproducible, so tests stay deterministic
/// without taking a dependency on an external RNG crate.
#[derive(Debug, Clone)]
pub struct SeededRng {
    state: u64,
}

impl SeededRng {
    /// Create a PRNG seeded with `seed`.
    #[must_use]
    pub fn seed_from_u64(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next raw 64-bit value (the `SplitMix64` step).
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Next value in the half-open interval `[0, 1)`.
    fn next_f64(&mut self) -> f64 {
        // Use the top 53 bits for a uniform double in [0, 1).
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    #[test]
    fn zero_qubits_fails_loud() {
        assert_eq!(StateVector::new(0), Err(QuantumError::ZeroQubits));
    }

    #[test]
    fn fresh_state_is_ket_zero() {
        let sv = StateVector::new(3).unwrap();
        assert_eq!(sv.dim(), 8);
        assert!(close(sv.probabilities()[0], 1.0));
        assert!(close(sv.norm_sqr(), 1.0));
    }

    #[test]
    fn h_on_zero_is_fifty_fifty() {
        let p = Circuit::new(1).h(0).run().unwrap().probabilities();
        assert!(close(p[0], 0.5));
        assert!(close(p[1], 0.5));
    }

    #[test]
    fn x_flips_zero_to_one() {
        let p = Circuit::new(1).x(0).run().unwrap().probabilities();
        assert!(close(p[0], 0.0));
        assert!(close(p[1], 1.0));
    }

    #[test]
    fn h_is_involution() {
        let sv = Circuit::new(1).h(0).h(0).run().unwrap();
        let p = sv.probabilities();
        assert!(close(p[0], 1.0));
        assert!(close(p[1], 0.0));
        // Amplitude, not just probability, returns to |0>.
        assert!(close(sv.amplitudes()[0].re, 1.0));
        assert!(close(sv.amplitudes()[0].im, 0.0));
    }

    #[test]
    fn bell_state_correlations() {
        let p = Circuit::new(2)
            .h(0)
            .cnot(0, 1)
            .run()
            .unwrap()
            .probabilities();
        assert!(close(p[0b00], 0.5));
        assert!(close(p[0b11], 0.5));
        assert!(close(p[0b01], 0.0));
        assert!(close(p[0b10], 0.0));
    }

    #[test]
    fn ghz_three_qubit() {
        let p = Circuit::new(3)
            .h(0)
            .cnot(0, 1)
            .cnot(1, 2)
            .run()
            .unwrap()
            .probabilities();
        assert!(close(p[0b000], 0.5));
        assert!(close(p[0b111], 0.5));
        // Every other basis state has zero probability.
        for (b, pb) in p.iter().enumerate() {
            if b != 0b000 && b != 0b111 {
                assert!(close(*pb, 0.0), "basis {b:03b} should be zero");
            }
        }
    }

    #[test]
    fn cz_phase_flip_on_eleven() {
        // |11> picks up a -1 phase under CZ; probabilities are unchanged but
        // the amplitude sign flips. Prepare |11> with two X gates.
        let sv = Circuit::new(2).x(0).x(1).cz(0, 1).run().unwrap();
        assert!(close(sv.amplitudes()[0b11].re, -1.0));
        assert!(close(sv.norm_sqr(), 1.0));
    }

    #[test]
    fn pauli_z_phase_and_s_t_diagonal() {
        // Z|1> = -|1>.
        let z = Circuit::new(1).x(0).z(0).run().unwrap();
        assert!(close(z.amplitudes()[1].re, -1.0));
        // S|1> = i|1>; T|1> = e^{iπ/4}|1>.
        let s = Circuit::new(1).x(0).s(0).run().unwrap();
        assert!(close(s.amplitudes()[1].re, 0.0) && close(s.amplitudes()[1].im, 1.0));
        let t = Circuit::new(1).x(0).t(0).run().unwrap();
        assert!(close(t.amplitudes()[1].re, std::f64::consts::FRAC_1_SQRT_2));
        assert!(close(t.amplitudes()[1].im, std::f64::consts::FRAC_1_SQRT_2));
    }

    #[test]
    fn rotations_match_analytic() {
        // Rx(π) maps |0> to -i|1>.
        let rx = Circuit::new(1).rx(0, std::f64::consts::PI).run().unwrap();
        assert!(close(rx.amplitudes()[1].im, -1.0));
        assert!(close(rx.probabilities()[1], 1.0));
        // Ry(π/2) on |0> gives equal real superposition.
        let ry = Circuit::new(1)
            .ry(0, std::f64::consts::FRAC_PI_2)
            .run()
            .unwrap();
        assert!(close(ry.probabilities()[0], 0.5));
        assert!(close(ry.probabilities()[1], 0.5));
        // Rz leaves computational-basis probabilities unchanged.
        let rz = Circuit::new(1).rz(0, 1.234).run().unwrap();
        assert!(close(rz.probabilities()[0], 1.0));
    }

    #[test]
    fn normalization_preserved_through_a_program() {
        let sv = Circuit::new(3)
            .h(0)
            .rx(1, 0.7)
            .ry(2, 1.9)
            .cnot(0, 2)
            .cz(1, 2)
            .t(0)
            .s(1)
            .run()
            .unwrap();
        assert!(close(sv.norm_sqr(), 1.0));
    }

    #[test]
    fn out_of_range_qubit_fails_loud() {
        let mut sv = StateVector::new(2).unwrap();
        assert_eq!(
            sv.apply(Gate::H, 2),
            Err(QuantumError::QubitOutOfRange {
                index: 2,
                num_qubits: 2,
            })
        );
        assert_eq!(
            sv.cnot(0, 5).unwrap_err(),
            QuantumError::QubitOutOfRange {
                index: 5,
                num_qubits: 2,
            }
        );
        // The circuit builder surfaces the same error from run().
        assert!(Circuit::new(2).h(3).run().is_err());
    }

    #[test]
    fn duplicate_qubit_two_qubit_gate_fails() {
        let mut sv = StateVector::new(2).unwrap();
        assert_eq!(sv.cnot(1, 1), Err(QuantumError::DuplicateQubit(1)));
        assert_eq!(sv.cz(0, 0), Err(QuantumError::DuplicateQubit(0)));
    }

    #[test]
    fn seeded_sampling_is_deterministic_and_tracks_born_rule() {
        // Bell state: only 00 and 11 ever appear, ~50/50, and the exact
        // counts are reproducible for a fixed seed.
        let sv = Circuit::new(2).h(0).cnot(0, 1).run().unwrap();
        let mut a = SeededRng::seed_from_u64(42);
        let mut b = SeededRng::seed_from_u64(42);
        let counts_a = sv.measure_counts(10_000, &mut a);
        let counts_b = sv.measure_counts(10_000, &mut b);
        assert_eq!(counts_a, counts_b, "same seed -> identical samples");
        assert_eq!(counts_a[0b01], 0);
        assert_eq!(counts_a[0b10], 0);
        // Within a few percent of 50/50.
        let frac00 = counts_a[0b00] as f64 / 10_000.0;
        assert!((frac00 - 0.5).abs() < 0.05, "got {frac00}");
        assert_eq!(counts_a.iter().sum::<usize>(), 10_000);
    }
}
