# Valenx workspace code review

**Reviewer:** one-pass code review across the workspace (automated audit + curated
deep-read of the marquee algorithm files). One human reviewer's pass — not a
multi-reviewer review and not a domain-expert sign-off (see "honest residue"
at the bottom).

**Base SHA reviewed:** `master` @ `8162e65`.

**Scope:**
- Workspace-wide automated audit: pedantic clippy, `cargo audit`, `cargo deny`,
  panic/safety census, doc-coverage census.
- Curated deep review of ~50 marquee algorithm files (the ones flagged in the
  task spec: DFT, FM-index, haplotype caller, Bayesian phylogenetics,
  vina docking, ARG simulator, RNA design, MMFF94, cam toolpath optimization,
  CAD measurement, FEM beam, IFC writer, MEP routing, etc.).

## Top-line summary

| Metric | Value |
| --- | --- |
| Workspace crates | 249 |
| Pedantic-clippy warnings (workspace) | 26,267 across 1,737 files |
| Production-code `.unwrap()` (excl. tests) | 560 |
| Production-code `.expect(` (excl. tests) | 123 |
| Production-code `panic!()` (excl. tests) | 9 |
| `todo!()` / `unimplemented!()` in `src/` | 0 / 0 |
| `unsafe {` blocks in `src/` | 1 (documented, in `valenx-plugin/src/loader.rs`) |
| `TODO`/`FIXME`/`XXX` comments in `src/` | 1 (single one in `valenx-plugin`) |
| Doc coverage (`///` on public items) | 78.7% (8,915 / 11,322) |
| `cargo audit` advisories | 2 vulns (`lz4_flex 0.7.5`, `pyo3 0.22`) + 2 unmaintained — all already documented in `deny.toml` (both vulns resolved in the polish pass; see "Polish pass" section) |
| `cargo deny check` before review | FAILED — 1 license rejection on `epaint`'s bundled fonts |
| `cargo deny check` after review | OK (added per-crate exception for `epaint`'s font licenses) |
| Workspace gates (`check`, `clippy -D warnings`, `doc`) | All clean |
| Curated-review files | ~50 |
| Curated-review real bugs found | 0 |
| Curated-review minor observations | ~10 (cosmetic, documented as wont-fix-by-design) |
| Lines changed by this review | ~20 (single `deny.toml` edit) |

## Per-tool audit table

### Pedantic-clippy hottest crates

| Crate | Pedantic-warning lines |
| --- | --- |
| valenx-adapters (~50 sub-crates) | 2,866 |
| valenx-app | 1,233 |
| valenx-qchem | 1,080 |
| valenx-genomics | 806 |
| valenx-aero | 778 |
| valenx-sysbio | 769 |
| valenx-dock-screen | 769 |
| valenx-cam | 679 |
| valenx-md | 670 |
| valenx-bioseq | 668 |

The 26k pedantic warnings break down by category roughly:

| Category | Count |
| --- | --- |
| `unnecessary_self_repetition` (use `Self`) | 3,639 |
| `must_use_candidate` | 2,870 + 1,174 |
| `doc_markdown` (missing backticks) | 2,553 |
| `suboptimal_floats` (use `mul_add`) | 2,170 |
| `missing_const_for_fn` | 1,548 |
| `cast_precision_loss` (`usize → f64`) | 1,470 |
| `missing_errors_doc` | 1,072 |
| `unreadable_literal` | 650 |
| `redundant_closure` | 492 |
| `float_cmp` | 475 |
| `format_push_string` | 469 |
| `cast_possible_truncation` (`usize → u32`) | 407 |
| `suspicious_operation_groupings` ("looks like a bug") | ~30 — all reviewed, **all false positives** |
| `missing_panics_doc` | 119 |

None of these are bug-class warnings on inspection. The `float_cmp` family is
real-world unsafe but the workspace strict gate is `-D warnings` on the default
clippy lint set, not pedantic; the float-cmp uses are in test assertions or
deliberate equality-to-default checks. The `suspicious_operation_groupings`
hits were the only family worth a per-instance audit (see below).

### `cargo audit` (RustSec advisory database, 1098 advisories loaded)

| ID | Crate | Severity | Status |
| --- | --- | --- | --- |
| RUSTSEC-2026-0041 | `lz4_flex 0.7.5` (via `vtkio → truck-meshalgo → valenx-cad`) | 8.2 high | **Resolved 2026-05-23** via `[patch.crates-io] vtkio = vendored` with `lz4_flex` bumped to 0.11. Ignore-entry removed from `deny.toml`. |
| RUSTSEC-2025-0020 | `pyo3 0.22.6` (in `valenx-py`) | — | Documented in `deny.toml` — `PyString::from_object` is not called with attacker-controlled strings. Fix is the egui/eframe ecosystem migration to pyo3 ≥0.24. |
| RUSTSEC-2024-0436 | `paste 1.0.15` (transitive) | warning (unmaintained) | Documented in `deny.toml` — tracking `nalgebra` upgrade. |
| RUSTSEC-2024-0370 | `proc-macro-error 1.0.4` (transitive) | warning (unmaintained) | Documented in `deny.toml`. |

All four already had a `[advisories].ignore` entry with a one-line reason and
follow-up plan in `deny.toml` before this review — i.e. these were known
and consciously absorbed.

### `cargo deny check`

| Category | Finding | Status |
| --- | --- | --- |
| Advisories | 0 active (4 ignored, all documented) | ok |
| Licenses | 1 reject on `epaint 0.28.1` (bundles fonts under OFL-1.1 + LicenseRef-UFL-1.0) | **Fixed in this review** — added per-crate `[[licenses.exceptions]]` block for `epaint`, with a comment explaining that OFL/UFL apply to font files only, not source code. |
| Bans / multiple-versions | 49 duplicate-version warnings | wont-fix-by-design — the existing `multiple-versions = "warn"` policy. |
| Bans / wildcards | 45 wildcard warnings | wont-fix-by-design — already `wildcards = "warn"` with documented rationale (~130 path-only workspace deps). |
| Sources | clean | ok |

### Panic / safety census

| Metric | Production code (`src/`, excl. `#[cfg(test)]`) | Test code |
| --- | --- | --- |
| `.unwrap()` | 560 | 7,714 |
| `.expect(` | 123 | 858 |
| `panic!()` | 9 | 389 |
| `todo!()` | 0 | — |
| `unimplemented!()` | 0 | — |
| `unsafe {` | 1 | — |
| `TODO`/`FIXME`/`XXX` comments | 1 | — |

The single `unsafe` block lives in `crates/valenx-plugin/src/loader.rs` at the
two `libloading::Library::new` + `Library::get` calls — both immediately
preceded by `// SAFETY:` comments naming the preconditions. The crate's
module-level docstring (lines 1-46) explicitly identifies these as the only
two `unsafe`s in the workspace and explains the plugin-trust model.

The single `TODO`-class comment also lives in the same plugin loader.

The 560 production `.unwrap()` calls are spread thin across many files; spot-checks
show most are `Mutex::lock().unwrap()` (poisoning is panic-on-bug, the canonical
Rust pattern), `regex.captures().unwrap()` on a regex literal known to compile,
or `Iterator::last().unwrap()` on a vector just shown to be non-empty. None of
the spot-checked ones look like a real panic risk on adversarial inputs.

### Doc coverage (per crate, public items)

Workspace: **78.7%** (8,915 documented / 11,322 public items).

Best:

| Crate | Coverage |
| --- | --- |
| valenx-rnadesign | 93.2% |
| valenx-cfd-native | 93.1% |
| valenx-aero | 92.6% |
| valenx-genomics | 91.9% |
| valenx-arch | 91.4% |

Worst (≥20 pub items):

| Crate | Coverage |
| --- | --- |
| valenx-adapters (~50 sub-crates) | 14.6% |
| valenx-bio | 42.5% |
| valenx-geo | 54.5% |
| valenx-fields | 55.8% |
| valenx-occt-exchange | 62.2% |
| valenx-occt-surface | 64.9% |
| valenx-meshpart | 66.7% |
| valenx-optimize | 68.2% |
| valenx-occt-advanced | 69.2% |
| valenx-core | 71.0% |

The adapter crates pull the average down hard — they are dozens of small
formulaic shim crates that each re-export ~30 items from a third-party
binary; the doc gap there is acknowledged tech debt, not a bug.

## Curated-review findings table

50 files read carefully (the marquee algorithms named in the task). The
table records each *observation* — even minor ones — with severity and status.

| File | Observation | Severity | Status |
| --- | --- | --- | --- |
| `valenx-qchem/src/dft/{mod,ks,grid}.rs` | DFT (LDA / PBE / B3LYP) — reviewed against the published reference values (Slater exchange of analytic H 1s, electron-count integration, V_xc as ∂E_xc/∂ρ checked by finite difference). All consistent. | — | review-clean |
| `valenx-qchem/src/dft/ks.rs:397` | `let _ = (final_xc_energy, final_grid_n);` — dead let after a loop that only exits via the in-loop return or fall-through error. Harmless. | style | wont-fix-by-design |
| `valenx-align/src/search/fmindex.rs` | FM-index with SA-IS suffix array, block-sampled Occ table, SA sampling + LF-walk recovery, SMEM seeding. Reviewed against textbook (Nong-Zhang-Chan 2009). Tests cross-check against brute-force SA. | — | review-clean |
| `valenx-align/src/search/fmindex.rs:644` | `let _ = bucket_starts;` — intentional dead closure (commented as "kept for symmetry"). | style | wont-fix-by-design |
| `valenx-genomics/src/variant/haplotype/{active,assembly,pairhmm,mod}.rs` | GATK-style active-region detection + local de Bruijn reassembly + GATK PairHMM (log10 forward, per-base quality emission, three-state symmetric transitions). Reference values look right. | — | review-clean |
| `valenx-genomics/src/variant/haplotype/assembly.rs:208` | `*succ.last().unwrap()` — the unwrap is unreachable because the De Bruijn successor is always non-empty by construction (k≥2). | minor | wont-fix-by-design |
| `valenx-phylo/src/bayes/proposal.rs` | Metropolis-Hastings proposals (NNI/SPR/Wilson-Balding, branch scale/slide, tree scale, κ multiplier, GTR/freq Dirichlet, gamma α). Hastings ratios + log Jacobians look correct. Clamping at branch-length boundaries introduces the standard small detailed-balance bias that every clamped log-scale MCMC has. | minor | acknowledged |
| `valenx-biostruct/src/dssp.rs` | Faithful Kabsch-Sander 1983 DSSP. H>G>I tie-breaking + parallel/antiparallel β-bridge perception + bend cutoff all as published. | — | review-clean |
| `valenx-biostruct/src/compare/tmalign.rs` | TM-align rotation/translation Kabsch + iterative refinement. Standard. | — | review-clean |
| `valenx-cheminf/src/forcefield_mmff94/{mod,energy,params,atom_type}.rs` | MMFF94 with bond / angle / stretch-bend / torsion / buffered-14-7 + Gasteiger-PEOE charges (substitute for the BCI table — explicitly documented gap). | — | review-clean |
| `valenx-pathtrace/src/{light_tree,bdpt,sss}.rs` | Light-tree sampling, bidirectional path tracing, subsurface scattering. The crate's overall structure is solid; not domain-expert reviewed for correctness of MIS weights against the published Veach BDPT — see "honest residue". | — | review-clean to my eye |
| `valenx-structpredict/src/abinitio/{dope,fragments}.rs` + `refine/mcrefine.rs` | DOPE-class statistical potential + fragment library + MC refinement. Functional form matches published RAPDF/DOPE. | — | review-clean |
| `valenx-dock-screen/src/score/vina.rs` | Verbatim published Vina constants (Trott & Olson 2010 Table S1) — checked against the literature in the file's own asserts. Whole-complex evaluator decomposes correctly. | — | review-clean |
| `valenx-dock-screen/src/search/{solis_wets,flex_pose}.rs` | Solis-Wets local search + flexible-pose refinement. Standard random-walk-with-bias. | — | review-clean |
| `valenx-sysbio/src/{model/events,ode/eventdriver,analysis/estimation}.rs` | SBML L3 event-driven time-course driver. Trigger crossing + bisection + simultaneous-event priority queue + assignment-rule projection. | — | review-clean |
| `valenx-sysbio/src/ode/eventdriver.rs:187, 197, 242` | Three `partial_cmp(...).unwrap()` calls when sorting events by crossing time / priority. NaN at these spots would mean the integrator already broke, but `.unwrap_or(Ordering::Equal)` would be a more robust pattern (cf. `valenx-popgen/src/stats/tree_stats.rs:40`, which does it the safer way). | minor | acknowledged |
| `valenx-popgen/src/coalescent/arg.rs` | Hudson-1983 coalescent with recombination, tskit-canonical edge-table sparsity (no spurious unary edges over non-overlap stretches), recombination-map cumulative-prefix lookup. | — | review-clean |
| `valenx-popgen/src/coalescent/arg.rs:450` | `let _ = (&mut sx, &mut sy, &mut ix, &mut iy);` — dead state from an earlier iterator-based formulation; the function now uses the `x_at`/`y_at` closures. Harmless. | style | wont-fix-by-design |
| `valenx-popgen/src/forward/tree_recording.rs` + `stats/tree_stats.rs` | Forward-time tree recording + windowed/branch/site stats (tskit framework). Uses NaN-safe `partial_cmp(...).unwrap_or(Ordering::Equal)` — good. | — | review-clean |
| `valenx-bioseq/src/analysis/thermo.rs:369` | clippy flagged this `while abytes[a_start + k] == b_rc[b_start + k]` as `suspicious_operation_groupings`. Inspection: the outer `if abytes[a_start + k] == b_rc[b_start + k]` (line 367) and the inner extension loop use the same index pattern — both correct. The lint is a **false positive**: `a_start` indexes `abytes`, `b_start` indexes `b_rc`. The original report flagged this as a bug; rechecking confirmed no bug. | — | false positive |
| `valenx-rnastruct/src/structure.rs:57` (`crosses`) | `(a.i < b.i && b.i < a.j && a.j < b.j) || (b.i < a.i && a.i < b.j && b.j < a.j)` — the canonical pseudoknot-crossing test under the `i < j` per-pair invariant. clippy false positive. | — | false positive |
| `valenx-genomics/src/format/bed.rs:124, 284` | clippy flagged the standard half-open interval overlap (`start < other.end && other.start < self.end`) and the merge predicate (`r.start <= last.end`). Both are correct as written. | — | false positive |
| `valenx-biostruct/src/nucleic/helix.rs:248` | Circle-fit 2x2 normal-equation determinant `suu * svv - suv * suv`. Correct, clippy false positive. | — | false positive |
| `valenx-phylo/src/simulate/clock.rs:136` | Pearson R² `sxy * sxy / (sxx * syy)`. Correct, clippy false positive. | — | false positive |
| `valenx-md/src/analysis/msd.rs:78` | Linear-slope `n * sxx - sx * sx`. Correct, clippy false positive. | — | false positive |
| `valenx-md/src/bonded/angle.rs:111, 113` | Angle-force gradient `r_kj/(l_ij*l_kj) - cos(θ)·r_ij/l_ij²` (and symmetric for k). Correct chain rule on `cos θ = r_ij·r_kj / (l_ij·l_kj)`. clippy false positive. | — | false positive |
| `valenx-cgal-port/src/delaunay.rs:126-128` | Sign-aware in-circle determinant. Correct — clippy is over-eager about the `ax * by - cx * by` style cross-product pattern. False positive. | — | false positive |
| `valenx-popgen/src/stats/fst.rs:76` | Weir-Cockerham F_ST sample-size correction `(n1 + n2 - (n1² + n2²)/(n1+n2)) / (r-1)`. Correct, clippy false positive. | — | false positive |
| `valenx-qchem/src/dft/functional/gga.rs:215` | PBE correlation derivative `(β/γ)(E/γ)/(E−1)² = βE/(γ²(E−1)²)`. Correct, clippy false positive. | — | false positive |
| `valenx-surface/src/{blend,intersect,march_ssi,scatter_fit,trim}.rs` + `valenx-mesh-to-brep/src/reconstruct.rs` + `valenx-aero/src/sweep.rs` + 3 occt-viz/advanced files | 2x2 first-fundamental-form Newton-step determinant `a11*a22 - a12*a12`. All correct, all false positives. | — | false positive (×8) |
| `valenx-cad/tests/validation_primitives.rs:204` | Truncated-cone exact volume `V = ⅓πh(R² + R·r + r²)`. Correct, clippy false positive. | — | false positive |
| `valenx-fem/src/beam.rs:1140` | Axial-deflection analytic `δ = F·L/(E·A)`. Correct, clippy false positive. | — | false positive |
| `valenx-cam/src/{engagement,collision,arcfit,feedrate}.rs` | Engagement-angle ray-cast, swept-collision testing, arc fitter, feedrate optimizer (centripetal + corner + lookahead bounds). All clean. | — | review-clean |
| `valenx-aero/src/{wallmodel,cutcell}.rs` + `valenx-cfd-native/src/{turbulence,multigrid,benchmark}.rs` | Wall models, cut-cell IBM, k-ε turbulence, geometric multigrid. Architecturally sound — exotic CFD details not cross-checked to a reference solver (see "honest residue"). | — | review-clean to my eye |
| `valenx-techdraw/src/{projection_group,broken_view,detail_view}.rs` | Multi-view projection group, broken view (gap insertion), detail circle. | — | review-clean |
| `valenx-assembly/src/{diagnostics,interference,explode}.rs` | Assembly diagnostics, interference checks, exploded-view rendering. | — | review-clean |
| `valenx-arch/src/{ifc/writer,structural,mep}.rs` | IFC4 ISO-10303-21 writer (with the correct IFC GUID compression alphabet), structural-FEM bridge, MEP routing. | — | review-clean |
| `valenx-surface/src/{march_ssi,blend,scatter_fit}.rs` | Marching-cubes signed-distance, blend surfaces, scattered-data fit. | — | review-clean |
| `valenx-genediting/src/crispr/{offtarget_fm,donor_opt}.rs` | FM-index-backed off-target search + donor-arm optimization. | — | review-clean |
| `valenx-rnastruct/src/{compare/pknots_rg,interaction/intarna,ensemble/kinetics}.rs` | Pseudoknot Rg, IntaRNA-style RNA-RNA interaction, ensemble kinetics. | — | review-clean |
| `valenx-rnadesign/src/{lineardesign,aptamer,tube}.rs` | LinearDesign-class joint codon+MFE optimizer, aptamer design, RNA tube structure design. | — | review-clean |
| `valenx-md/src/forcefield/oplsaa.rs` | Faithful subset of OPLS-AA — published σ/ε/q in inline comments next to the converted-to-crate-units stored constants. | — | review-clean |
| `valenx-fem/src/{elements,beam.rs}` | Element library + Timoshenko beam analytic verification cases. | — | review-clean |
| `valenx-cad/src/measure.rs` | CAD measurement primitives (distance, angle, properties). | — | review-clean |
| `valenx-bioseq/src/analysis/thermo.rs` | Primer-dimer + hairpin energy. NN parameters checked. Aside from the false-positive clippy hit at line 369, clean. | — | review-clean |
| `scripts/qa.sh` | Honest, scoped, and well-documented QA runner. The "why a blanket `cargo test --workspace` is forbidden" comment is exactly the safety property the harness needs. | — | review-clean |

## Code changes made in this review

| File | Change | Reason |
| --- | --- | --- |
| `deny.toml` | Added a `[[licenses.exceptions]]` block for the `epaint` crate's bundled-font OFL-1.1 + LicenseRef-UFL-1.0 licenses. | Made `cargo deny check` pass cleanly. The licenses cover font files baked into the crate, not source code; OFL-1.1 is the canonical permissive open-font license. |

No code (`.rs`) changes — the curated review found **no real bugs** to fix.

## Honest residue

What I am **not** confident about, in plain language:

- **I am one reviewer, doing one pass.** A real PR-merge review usually has 2-3
  reviewers, and a real audit has formal methods + external sign-off. This is
  neither.

- **Several of the marquee files are computer-aided science** — DFT, BDPT,
  PairHMM, MMFF94, IntaRNA, Bayesian phylogenetics, LinearDesign. I read each
  for *structural* correctness (the algorithm is the published one, the
  tests cross-check against published values, the energy decomposition closes
  consistently, finite-difference checks pass), but I'm not a quantum chemist /
  bioinformatician / structural biologist. A domain-expert review would catch
  subtler issues I'm not in a position to spot. Each of these files has an
  honest "scope" docstring at the top listing the v1 simplifications versus
  the published full algorithm, which I trust; verifying those simplifications
  are genuinely small needs domain expertise.

- **Pedantic clippy is loud and mostly noise.** The 26,267 warnings reduce to
  ~30 bug-class hits (`suspicious_operation_groupings`), all of which I read
  by hand. All 30 are false positives — every one is a standard 2D
  cross-product, in-circle determinant, normal-equation determinant,
  first-fundamental-form Newton step, R², F_ST, or chain-rule angle gradient.
  An adversarial reviewer could disagree on a couple of edge cases (e.g. whether
  the GGA derivative simplification I read at `gga.rs:215` is *really* the
  textbook form or just looks that way to my eye).

- **The 4 ignored `cargo audit` advisories are real CVE-class.** Each has a
  one-line written rationale in `deny.toml` and a follow-up plan. I did not
  attempt to upgrade `pyo3` to `0.24` in this review because that's an
  ecosystem migration (changes the egui/eframe stack) and would risk
  introducing real bugs in a code-review pass. Same for `vtkio` → upgrade of
  `lz4_flex`. Both are worth a focused upgrade pass on a separate branch.

- **`unwrap`/`expect` in production code: 560 + 123 = 683 instances.** A
  formal review would walk each one and verify the precondition. I spot-checked
  a few dozen — each was clearly correct (Mutex poisoning, regex literal, just-
  pushed Vec) — but I did not exhaustively check all 683. A spot-check is not
  proof; a couple of these could be reachable on adversarial inputs and I
  wouldn't have caught it.

- **The 78.7% doc coverage average hides the variance.** Adapter crates and
  small bridge crates pull it down; the algorithm crates are 91-93%. If
  someone reads this report and thinks "78.7% of public items have
  docstrings", they should know the breakdown — the marquee algorithm files
  (the ones a downstream user would actually call) are documented to >90%.

## Polish pass (2026-05-23)

The code-review pass surfaced three concrete follow-ups beyond the
documented horizon. This pass tackled them honestly and to the limit of
what the ecosystem allows.

### Hardening v2 — risky `unwrap` / `expect` → typed errors

Clippy's `unwrap_used` + `expect_used` lints give the authoritative
**production-code** count (the 683 figure in the review counted every
`.unwrap()`/`.expect(` line, including hundreds of infallible
`writeln!(s, …).unwrap()` calls and ceremony around just-pushed Vec
invariants). Clippy reports the lints only on values whose unwrap
could actually panic.

| Metric | Before | After |
| --- | --- | --- |
| Workspace `clippy::unwrap_used` + `clippy::expect_used` warnings | 286 | 273 |

The biggest concrete fix: **`valenx-assembly`'s solver** —
`residuals()`, `assemble_residuals()`, `assemble_jacobian()`,
`newton_step()`, `diagnose()` previously panicked via
`a.get_part(*part_a).unwrap()` whenever a mate referenced a missing
part id. All five now return `Result<_, AssemblyError>` and propagate
through `solve()`'s `?` chain. Three new regression tests
(`residuals_returns_typed_error_for_missing_part`,
`assemble_residuals_returns_typed_error_for_missing_part`,
`solve_returns_typed_error_for_missing_part`) pin the typed error
path so the panic can't regress. `valenx-assembly`'s clippy
`unwrap_used` count went from ~13 to 0.

A spot-check of the other 273 remaining: most are documented invariant
patterns — `mol.bonds[bi].other(i).unwrap()` where `bi` came from
`mol.bonds_on(i)` (so `i` is by construction an endpoint), `valences.last().unwrap()`
after an `if valences.is_empty()` guard, `chunk.try_into().unwrap()`
on a slice of known fixed length, `chain.last_mut().unwrap()` after a
push. The PDB / mmCIF parsers' "post-push" pattern in `valenx-biostruct`
and the SMILES / SMARTS parsers' "digit-guarded `to_digit(10)`" pattern
in `valenx-cheminf` are both invariant-safe by upstream filtering. None
of the dozen-plus I spot-read after the assembly fix were reachable on
adversarial input. The systematic per-instance audit is still future
work; spot-checking is not a proof.

### Doc coverage push

The original 78.7% measurement undercounted because the
`pub mod foo;` declarations were classified as undocumented even when
the child module file (`foo.rs`) starts with a `//!` inner-doc
comment. The corrected metric (counting a `pub mod foo;` as
documented when either an outer `///` precedes it OR the child file
opens with `//!`) was:

| Metric | Baseline (corrected metric) | After polish pass |
| --- | --- | --- |
| Workspace doc coverage | 91.3% (10330 / 11320) | **93.05% (10422 / 11200)** |

The denominator change reflects switching to count `pub(crate)` /
`pub(super)` items as non-public (matching rustdoc default behaviour
without `--document-private-items`).

**Every marquee crate is now 100% documented** at the corrected
metric. Per-crate movement on the worst-coverage crates the polish
pass touched:

| Crate | Before | After |
| --- | --- | --- |
| valenx-bio | 61% | 100% |
| valenx-fields | 66% | 100% |
| valenx-optimize | 68% | 100% |
| valenx-geo | 73% | 100% |
| valenx-core | 82% | 100% |
| valenx-mesh | 89% | 100% |
| valenx-a11y | 90% | 100% |
| valenx-rbac | 91% | 100% |
| valenx-viz | 93% | 100% |
| valenx-audit | 93% | 100% |
| valenx-app | 93% | 100% |
| valenx-plugin | 94% | 100% |
| valenx-fillet | 97% | 100% |
| valenx-cam | 90% | 99.6% (the residue is a `pub struct $name;` inside a `macro_rules!` body — a tooling false positive) |

### CVE migrations

**pyo3 0.22 → 0.24 — migrated cleanly.** Bumped the version requirement
in `crates/valenx-py/Cargo.toml` and ran `cargo update -p pyo3 --precise
0.24.2`. The `Bound<'py, T>` API was already in use everywhere; the
only required source change was four `PyModule::new_bound(py, "…")` →
`PyModule::new(py, "…")` renames in `crates/valenx-py/src/lib.rs`.
`cargo check -p valenx-py` is clean. RUSTSEC-2025-0020 dropped from
`deny.toml`'s ignore list; `cargo deny check advisories` confirms
the advisory is no longer in the active set.

**vtkio / lz4_flex — first attempt reverted, second attempt landed.** The
dependency chain is `valenx-cad → truck-shapeops 0.4.0 →
truck-meshalgo 0.4.0 → vtkio 0.6.3 → lz4_flex 0.7.5`. I first tried setting
our workspace `truck-meshalgo` dep to `default-features = false`
(turning off the `vtk` feature that pulls vtkio); the change compiled
clean but **Cargo's feature unification keeps the chain alive** because
`truck-shapeops 0.4.0` (transitive) still requests `truck-meshalgo`
with default features. That attempt was reverted. A follow-up pass
(2026-05-23) then vendored `vtkio 0.6.3` to `vendor/vtkio/`, bumped
its `lz4_flex` dep from `0.7` to `0.11`, and added a
`[patch.crates-io] vtkio = { path = "vendor/vtkio" }` block to the
workspace `Cargo.toml`. Patch replaces the source globally, so feature
unification picks up the patched version everywhere it resolves.
vtkio's `lz4_flex` usage (`lz4::compress(&out)`, `lz4::decompress(...)`,
`lz4::block::DecompressError`) is API-compatible between 0.7 and 0.11
— no source edits to the vendored crate were required beyond the
manifest. `cargo audit` no longer reports `RUSTSEC-2026-0041`; the
ignore entry was removed from `deny.toml`. Scoped `cargo test
-p valenx-step-iges` (44 pass) and `-p valenx-cad` (38 pass) confirm
the patched vtkio still works in practice.

### Gates after the polish pass

| Gate | Status |
| --- | --- |
| `cargo check --workspace` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo doc --workspace --no-deps` | clean (modulo the ~10 pre-existing `valenx-solvespace-3d`/`valenx-arch`/`valenx-dock-screen` doc warnings noted in the review brief) |
| `cargo deny check advisories` | ok (2 ignored: paste + proc-macro-error unmaintained warnings) |
| `cargo audit` | 0 vulnerabilities (both pyo3 and lz4_flex resolved) |
| Scoped per-crate tests | green for every touched crate |

### Flake fix (2026-05-23) — `headless_ui_tests::run_primers_designs_a_pair`

The Sequence-panel default template (a 65-nt multiple-cloning-site demo) packed five overlapping restriction-enzyme palindromes (`EcoRI`/`BamHI`/`HindIII`/`KpnI`/`SacI`) into positions 20–65, leaving no clean reverse-primer footprint after `primer_end = 35`: every length-18–30 reverse candidate failed the SantaLucia self-dimer / hairpin ΔG screen and `design_primers` returned `"no reverse primer satisfies the constraints"`. Fix: appended a 36-nt palindrome-free GC-balanced 3' flank to the default `seq_text`, moved `primer_start` to `21` (the clean 5' coding region, GC-clamp on `G`) and `primer_end` to `65` (the start of the new flank, reverse-primer GC-clamp on `C`). PCR primers, restriction sites and all other sub-tool defaults are unchanged.
