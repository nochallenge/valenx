# Phase 13 — ML-assisted surrogates

**Status:** ⚪ Planned.

## Goal

Speed up design loops by training surrogates on validated runs —
transparently, with honest uncertainty reporting so nobody mistakes
an interpolation for a solve.

## Capability inventory

- Surrogate types: Gaussian process, random forest, multi-layer
  perceptron, physics-informed neural net.
- Training-data manager tied to provenance hashes — a surrogate
  declares exactly which runs trained it.
- Active learning: the surrogate asks for the next most-informative
  full-fidelity run.
- Uncertainty surface: every prediction ships an error bar; the UI
  uses shade / width to convey it.
- Model cards: a README-like page per surrogate with its inputs,
  outputs, training set, accuracy metrics, and known failure modes.

## Acceptance checklist

- [ ] Train a GP surrogate on a parameter sweep; predict new samples
      with uncertainty.
- [ ] Active-learning loop that converges to a target tolerance.
- [ ] Surrogate badge visible wherever predicted (not simulated)
      values are shown.
- [ ] Export surrogate as an ONNX model for reuse outside Valenx.

## Leads into

[Phase 14 — Plugin marketplace](./phase-14-plugins.md).
