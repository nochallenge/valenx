# valenx-neuro Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `valenx-neuro` — a native-Rust neural-interface / BCI simulation suite (extracellular FEM field, Hodgkin–Huxley cable, Rattay activating function, Pennes bioheat, electrode–tissue impedance) plus a `valenx-app` workbench, each module validated against a closed-form/textbook result.

**Architecture:** New crate `crates/valenx-neuro` (mirrors `valenx-reactdyn`). The extracellular field and bioheat reuse `valenx-fem`'s `solve_steady_thermal` (the operator −∇·(σ∇φ)=I is identical to −∇·(k∇T)=q). The HH cable is a hand-rolled RK4 integrator. The activating function couples the 3-D field to the 1-D cable by sampling φ along each fiber. Compute-then-visualize, on branch `feat/valenx-neuro` (draft PR #9).

**Tech Stack:** Rust, `valenx-fem` (FEM solver reuse), `valenx-core` (no hard dep expected), `egui`/`egui_plot` (workbench), `valenx-app` (host). Tests scoped `cargo test -p valenx-neuro` / `-p valenx-app` — **never** `cargo test --workspace` (hangs; see docs/QA.md).

**Reference conventions (RFC 0011 §Units):** potential **mV**, time **ms**, current **µA**, conductivity σ **S/m**, tissue length **mm**, compartment length / fiber diameter **µm**, C_m **µF/cm²**, conductances **mS/cm²**, ΔT **K**.

**Commit identity:** local git config is locked to `nochallenge <…noreply…>` — a plain `git commit` authors correctly. Commit each task individually.

---

## File Structure

**Create — `crates/valenx-neuro/`:**
- `Cargo.toml` — mirrors reactdyn; deps `valenx-fem = { path = "../valenx-fem" }`; `[lints.rust] missing_docs = "warn"`.
- `src/lib.rs` — crate docs + `pub mod` declarations + re-exports.
- `src/units.rs` — unit constants + converters + round-trip tests (Task 0).
- `src/error.rs` — `NeuroError` enum.
- `src/cable.rs` — HH gating + single & multi-compartment cable + RK4 (Tasks 1–2).
- `src/field.rs` — extracellular field wrapping `valenx_fem::solve_steady_thermal` (Task 3).
- `src/scene.rs` — procedural tissue grid + electrode + axon paths (Task 3 helper).
- `src/activating.rs` — sample φ along a fiber + activating function (Task 4).
- `src/engine.rs` — orchestrate coupled run → `Trajectory`/`Recruitment` (Task 5).
- `src/bioheat.rs` — Pennes steady solve (Task 6).
- `src/impedance.rs` — disk-electrode R_a + CPE (Task 7).

**Modify — workspace + app:**
- `Cargo.toml` (workspace root) — add `crates/valenx-neuro` to `members`.
- `crates/valenx-app/Cargo.toml` — add `valenx-neuro = { path = "../valenx-neuro" }`.
- `crates/valenx-app/src/lib.rs` — `mod neuro_workbench;` + `show_neuro_workbench: bool` field on the App struct.
- `crates/valenx-app/src/update.rs` — View-menu checkbox (near line 374) + draw dispatch (near line 862).
- `crates/valenx-app/src/neuro_workbench.rs` (create) — panel + 3-D playback + plots (Task 8).

---

## Task 0: Crate scaffold + units

**Files:** Create `crates/valenx-neuro/Cargo.toml`, `src/lib.rs`, `src/units.rs`, `src/error.rs`; modify workspace `Cargo.toml`.

- [ ] **Step 1 — scaffold.** Create `Cargo.toml` mirroring `crates/valenx-reactdyn/Cargo.toml` (workspace-inherited fields, `description` for the neuro suite, `[dependencies] valenx-fem = { path = "../valenx-fem" }`, `[lints.rust] missing_docs = "warn"`). Add `"crates/valenx-neuro"` to the root `Cargo.toml` `members`. Create `src/lib.rs` with `#![doc = "..."]` and `pub mod units; pub mod error;`. Create `src/error.rs` with `#[derive(Debug, thiserror-free)] pub enum NeuroError { BadScene(&'static str), Solver(String) }` + `Display`/`std::error::Error` (mirror reactdyn's `error.rs` — no thiserror).

- [ ] **Step 2 — failing units test.** In `src/units.rs`:
```rust
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
    }
    #[test]
    fn conductance_density_area_scaling() {
        // g_Na = 120 mS/cm² on 1e-4 cm² patch = 0.012 mS = 1.2e-5 S
        assert!((ms_per_cm2_to_s(120.0, 1e-4) - 1.2e-5).abs() < 1e-12);
    }
}
```

- [ ] **Step 3 — run, expect FAIL.** `cargo test -p valenx-neuro units` → fails (functions undefined).

- [ ] **Step 4 — implement** `mv_to_volts`, `volts_to_mv`, `ua_to_amp`, `amp_to_ua`, `ms_per_cm2_to_s(g_ms_cm2, area_cm2)`, plus `pub const` factors, all `///`-documented.

- [ ] **Step 5 — run, expect PASS.** `cargo test -p valenx-neuro` (all green) + `cargo clippy -p valenx-neuro --all-targets -- -D warnings`.

- [ ] **Step 6 — commit.** `git add crates/valenx-neuro Cargo.toml && git commit -m "feat(neuro): crate scaffold + unit conversions"`

---

## Task 1: Hodgkin–Huxley single compartment → textbook action potential

**Files:** `crates/valenx-neuro/src/cable.rs` (create); `src/lib.rs` (+`pub mod cable;`).

**Pinned constants (HH 1952, modern mV convention, V_rest = −65 mV):**
`C_m = 1.0` µF/cm²; `g_Na = 120.0`, `g_K = 36.0`, `g_L = 0.3` mS/cm²; `E_Na = 50.0`, `E_K = -77.0`, `E_L = -54.4` mV.
Rate functions (V in mV, rates in 1/ms):
```
α_m = 0.1·(V+40)/(1 − e^{−(V+40)/10})      β_m = 4·e^{−(V+65)/18}
α_h = 0.07·e^{−(V+65)/20}                   β_h = 1/(1 + e^{−(V+35)/10})
α_n = 0.01·(V+55)/(1 − e^{−(V+55)/10})      β_n = 0.125·e^{−(V+65)/80}
```
(Guard the two `x/(1−e^{−x})` forms against 0/0 with a series limit near the singularity.)

- [ ] **Step 1 — failing test:**
```rust
#[test]
fn suprathreshold_stimulus_fires_action_potential() {
    let mut c = HhCompartment::at_rest();          // V≈-65, m,h,n at steady state
    // 10 µA/cm² for 0.5 ms then off, integrate 20 ms at dt=0.005 ms (RK4)
    let trace = c.run(StimPulse { amp_ua_cm2: 10.0, start_ms: 1.0, width_ms: 0.5 },
                      20.0, 0.005);
    let vmax = trace.iter().cloned().fold(f64::MIN, f64::max);
    assert!(vmax > 20.0, "AP should overshoot 0 mV; vmax={vmax}");   // ~+40 mV
    assert!(vmax < 60.0, "but not blow up; vmax={vmax}");
}
#[test]
fn subthreshold_stimulus_does_not_fire() {
    let mut c = HhCompartment::at_rest();
    let trace = c.run(StimPulse { amp_ua_cm2: 2.0, start_ms: 1.0, width_ms: 0.5 },
                      20.0, 0.005);
    let vmax = trace.iter().cloned().fold(f64::MIN, f64::max);
    assert!(vmax < -40.0, "weak stim must not trigger an AP; vmax={vmax}");
}
#[test]
fn resting_state_is_stable() {
    let mut c = HhCompartment::at_rest();
    let trace = c.run(StimPulse { amp_ua_cm2: 0.0, start_ms: 0.0, width_ms: 0.0 },
                      10.0, 0.005);
    let drift = trace.last().unwrap() - (-65.0);
    assert!(drift.abs() < 1.0, "rest must not drift; drift={drift}");
}
```

- [ ] **Step 2 — run, expect FAIL.** `cargo test -p valenx-neuro cable` → fails (types undefined).

- [ ] **Step 3 — implement** `HhCompartment { v, m, h, n }`, `at_rest()` (set gates to α/(α+β) at −65), the ionic currents, `StimPulse`, and `run(...)` doing RK4 on the 4-state vector `(V,m,h,n)` with the pinned constants. Choose `dt = 0.005 ms` (conservative for HH stiffness; the test uses it). Return `Vec<f64>` of V samples.

- [ ] **Step 4 — run, expect PASS.** `cargo test -p valenx-neuro cable`.

- [ ] **Step 5 — refractory test** (add, then confirm green):
```rust
#[test]
fn absolute_refractory_blocks_immediate_second_spike() {
    let mut c = HhCompartment::at_rest();
    // two identical pulses 2 ms apart — the second lands in the refractory period
    let trace = c.run_two(StimPulse{amp_ua_cm2:10.0,start_ms:1.0,width_ms:0.5},
                          StimPulse{amp_ua_cm2:10.0,start_ms:3.0,width_ms:0.5},
                          20.0, 0.005);
    // count peaks above 0 mV → exactly 1
    assert_eq!(count_spikes(&trace, 0.0), 1);
}
```
Implement `run_two` + `count_spikes` (a peak-above-threshold counter), confirm `cargo test -p valenx-neuro` green + clippy clean.

- [ ] **Step 6 — commit.** `git commit -am "feat(neuro): Hodgkin-Huxley compartment — validated action potential"`

---

## Task 2: Multi-compartment cable → propagating AP + conduction velocity

**Files:** `crates/valenx-neuro/src/cable.rs` (extend).

**Axial coupling:** between compartments of length `Δx` (µm), radius `a` (µm), axial resistivity `R_i` (Ω·cm): `I_axial,k = (a/(2 R_i Δx²))·(V_{k-1} − 2V_k + V_{k+1})` (consistent units via `units.rs`). Sealed ends (no-flux) at both tips.

- [ ] **Step 1 — failing test:**
```rust
#[test]
fn action_potential_propagates_with_finite_velocity() {
    // 200 compartments, Δx=100 µm, a=238 µm (squid), R_i=35.4 Ω·cm
    let mut cable = HhCable::uniform(200, 100.0, 238.0, 35.4);
    // stimulate compartment 0; record AP peak time at comp 50 and comp 150
    let r = cable.stimulate_end(StimPulse{amp_ua_cm2:50.0,start_ms:1.0,width_ms:0.5},
                                30.0, 0.005);
    let t50  = r.peak_time_ms(50).expect("comp 50 fires");
    let t150 = r.peak_time_ms(150).expect("comp 150 fires");
    assert!(t150 > t50, "AP must travel 50→150 in +x");
    let dist_mm = (150 - 50) as f64 * 0.1;          // 10 mm
    let vel_m_s = (dist_mm / (t150 - t50)) * 1.0;    // mm/ms == m/s
    assert!((1.0..100.0).contains(&vel_m_s), "squid-ish CV; got {vel_m_s} m/s");
}
```

- [ ] **Step 2 — run, expect FAIL.** `cargo test -p valenx-neuro propagat`.

- [ ] **Step 3 — implement** `HhCable { comps: Vec<HhState>, dx_um, a_um, ri_ohm_cm }`, `uniform(...)`, `stimulate_end(...)` (RK4 on the whole state vector each step, adding `I_axial` per compartment), and a `CableRun` with `peak_time_ms(idx)` (argmax of V at that compartment, returns `Some(t)` only if it crosses 0 mV).

- [ ] **Step 4 — run, expect PASS** + clippy clean.

- [ ] **Step 5 — commit.** `git commit -am "feat(neuro): multi-compartment cable — propagating AP + conduction velocity"`

---

## Task 3: Extracellular field (reuse valenx-fem) → point-source φ = I/(4πσr)

**Files:** `crates/valenx-neuro/src/scene.rs` (create), `src/field.rs` (create); `src/lib.rs` (+mods).

- [ ] **Step 1 — confirm the fem API.** Read `crates/valenx-fem/src/thermal_solver.rs:150-240` to record the exact signature of `solve_steady_thermal(...)` and the fields of `FixedTemperature`, `HeatLoad`, `ThermalSolution`, and the node/Tet4 mesh types it expects. (Recorded surface: `pub fn solve_steady_thermal(...)`, structs `FixedTemperature`, `HeatLoad`, `ThermalSolution`, enum `ThermalSolverError`.) Write `field.rs` to call it with σ in place of k.

- [ ] **Step 2 — failing test:**
```rust
#[test]
fn point_source_matches_inverse_r_law() {
    // 21×21×21 structured tet grid, 40 mm cube, σ=0.2 S/m (gray matter)
    let mesh = TissueGrid::cube(40.0, 21, 0.2);
    // inject I=100 µA at the center node, ground the outer boundary
    let field = solve_extracellular(&mesh, ElectrodeSource::point_center(100.0))
        .expect("solve");
    // φ should fall as 1/r: φ(r)/φ(2r) ≈ 2 along +x away from the source
    let p1 = field.potential_mv_at_radius_x(5.0);   // 5 mm
    let p2 = field.potential_mv_at_radius_x(10.0);   // 10 mm
    assert!((p1 / p2 - 2.0).abs() < 0.4, "1/r law: {p1}/{p2} = {}", p1/p2);
    // absolute magnitude vs φ=I/(4πσr): at 5 mm, analytic ≈ ... (assert ±25%)
    let analytic_mv = analytic_point_source_mv(100.0, 0.2, 5.0);
    assert!((p1/analytic_mv - 1.0).abs() < 0.25, "got {p1}, analytic {analytic_mv}");
}
```

- [ ] **Step 3 — run, expect FAIL.** `cargo test -p valenx-neuro field`.

- [ ] **Step 4 — implement** `TissueGrid::cube(side_mm, n, sigma)` (structured node grid → 6 tets/cell), `ElectrodeSource::point_center` (a `HeatLoad`-equivalent nodal current injection), `solve_extracellular(...)` (assemble σ-conductivity via `solve_steady_thermal`, boundary nodes as `FixedTemperature{value:0}`), `ExtracellularField` with `potential_mv_at_radius_x` (nearest-node or trilinear sample) and the free fn `analytic_point_source_mv(i_ua, sigma_s_m, r_mm)` = `I/(4πσr)` in mV.

- [ ] **Step 5 — run, expect PASS** + clippy clean.

- [ ] **Step 6 — commit.** `git commit -am "feat(neuro): extracellular FEM field via valenx-fem — point-source validated"`

---

## Task 4: Activating function → cathodic<anodic, strength–distance ∝ r²

**Files:** `crates/valenx-neuro/src/activating.rs` (create); `src/lib.rs` (+mod).

**Activating function** along a fiber sampled at points `x_k`: `f_k ∝ (Vₑ(x_{k-1}) − 2Vₑ(x_k) + Vₑ(x_{k+1}))/Δx²`. Cathodic source = negative current.

- [ ] **Step 1 — failing test:**
```rust
#[test]
fn cathodic_depolarizes_under_electrode() {
    // straight fiber 1 mm below a point electrode; sample Ve along it
    let f_cath = activating_profile_for(/*current_ua*/ -100.0, /*depth_mm*/ 1.0);
    let f_anod = activating_profile_for(/* +100 */         100.0,             1.0);
    // node nearest the electrode (center of the profile)
    let mid = f_cath.len()/2;
    assert!(f_cath[mid] > 0.0, "cathodic → depolarizing at the near node");
    assert!(f_anod[mid] < 0.0, "anodic → hyperpolarizing at the near node");
}
#[test]
fn recruitment_threshold_scales_with_distance_squared() {
    // find threshold current to fire a fiber at distances 0.5,1,2 mm
    let mut log_r = vec![]; let mut log_i = vec![];
    for r in [0.5_f64, 1.0, 2.0] {
        let i_th = threshold_current_ua(r);
        log_r.push(r.ln()); log_i.push(i_th.ln());
    }
    // slope of ln(I_th) vs ln(r) ≈ 2
    let slope = lin_slope(&log_r, &log_i);
    assert!((slope - 2.0).abs() < 0.5, "strength-distance exponent={slope}");
}
```

- [ ] **Step 2 — run, expect FAIL.** `cargo test -p valenx-neuro activating`.

- [ ] **Step 3 — implement** `activating_profile_for(current_ua, depth_mm)` (solve field, sample Ve along the fiber, second difference), `threshold_current_ua(r_mm)` (bisection: smallest |I| that makes the coupled cable fire — depends on Task 5's `engine`, so stub a minimal couple-and-check here or land Task 5 first; see note), and `lin_slope`.

  **Note:** Task 4's `threshold_current_ua` needs the coupled cable solve. Build Task 5's `engine::stimulate` first if cleaner, then return here — keep the *sign* test (no cable) in Task 4 and the *r²* test after Task 5. Adjust order at execution.

- [ ] **Step 4 — run, expect PASS** + clippy clean.

- [ ] **Step 5 — commit.** `git commit -am "feat(neuro): activating function — cathodic/anodic sign + strength-distance r^2"`

---

## Task 5: Coupled stimulation → recruitment

**Files:** `crates/valenx-neuro/src/engine.rs` (create); `src/lib.rs` (+mod).

- [ ] **Step 1 — failing test:**
```rust
#[test]
fn electrode_recruits_near_axon_above_threshold() {
    let scene = Scene::single_axon(/*depth_mm*/ 1.0);
    let below = stimulate(&scene, StimPulse{amp_ua_cm2:0.0,..}, /*electrode_ua*/ 20.0);
    let above = stimulate(&scene, StimPulse{amp_ua_cm2:0.0,..}, /*electrode_ua*/ 400.0);
    assert!(!below.any_fired(), "20 µA must not recruit");
    assert!(above.any_fired(), "400 µA must recruit the near axon");
}
#[test]
fn recruitment_curve_is_monotonic() {
    let scene = Scene::bundle(/*n*/ 20, /*spread_mm*/ 2.0);
    let frac: Vec<f64> = [50.0,100.0,200.0,400.0,800.0].iter()
        .map(|&i| stimulate(&scene, StimPulse::default(), i).recruited_fraction())
        .collect();
    for w in frac.windows(2) { assert!(w[1] >= w[0], "recruitment must not decrease"); }
}
```

- [ ] **Step 2 — run, expect FAIL.** `cargo test -p valenx-neuro engine`.

- [ ] **Step 3 — implement** `Scene` (tissue grid + electrode + `Vec<axon path>`), `stimulate(scene, intracellular_stim, electrode_ua) -> Recruitment` (solve field once for unit current, scale by `electrode_ua`, for each axon sample Ve → drive the `HhCable` with the extracellular activating term → record fired/not), `Recruitment { fired: Vec<bool> }` with `any_fired()` and `recruited_fraction()`.

- [ ] **Step 4 — run, expect PASS** + clippy clean. Then return to Task 4's `r²` test (now buildable) and confirm green.

- [ ] **Step 5 — commit.** `git commit -am "feat(neuro): coupled field↔cable stimulation + recruitment curve"`

---

## Task 6: Bioheat (Pennes) → ΔT vs analytic point source

**Files:** `crates/valenx-neuro/src/bioheat.rs` (create); `src/lib.rs` (+mod).

**Pennes steady:** `∇·(k∇T) − ω_b ρ_b c_b (T−T_a) + Q = 0`. Reuse `solve_steady_thermal` (conduction) + perfusion as a diagonal reaction term + source `Q = σ|∇φ|²` (Joule). Brain: `k≈0.5` W/m·K, `ω_b≈0.008` 1/s, `ρ_b c_b≈3.6e6` J/m³K, `T_a=37 °C`.

- [ ] **Step 1 — failing test:**
```rust
#[test]
fn point_heat_source_matches_analytic_no_perfusion() {
    let mesh = TissueGrid::cube(40.0, 21, /*sigma unused here*/ 0.2);
    let dt = solve_bioheat(&mesh, PointHeat{ power_w: 0.01 }, Perfusion::off(), 0.5);
    let num = dt.delta_t_k_at_radius_x(5.0);          // 5 mm
    let ana = 0.01 / (4.0*std::f64::consts::PI*0.5*0.005);  // Q/(4πk r), r in m
    assert!((num/ana - 1.0).abs() < 0.25, "ΔT num={num} analytic={ana}");
}
#[test]
fn perfusion_reduces_far_field_temperature() {
    let mesh = TissueGrid::cube(40.0, 21, 0.2);
    let off = solve_bioheat(&mesh, PointHeat{power_w:0.01}, Perfusion::off(), 0.5)
        .delta_t_k_at_radius_x(8.0);
    let on  = solve_bioheat(&mesh, PointHeat{power_w:0.01}, Perfusion::brain(), 0.5)
        .delta_t_k_at_radius_x(8.0);
    assert!(on < off, "perfusion must cool the far field: on={on} off={off}");
}
```

- [ ] **Step 2 — run, expect FAIL.** `cargo test -p valenx-neuro bioheat`.

- [ ] **Step 3 — implement** `solve_bioheat(mesh, source, perfusion, k)` (assemble conduction via the fem path, add `ω_b ρ_b c_b · V_node` to the diagonal and to the RHS via `T_a`, inject `Q`), `PointHeat`, `Perfusion::{off, brain}`, `BioheatSolution::delta_t_k_at_radius_x`.

- [ ] **Step 4 — run, expect PASS** + clippy clean.

- [ ] **Step 5 — commit.** `git commit -am "feat(neuro): Pennes bioheat — analytic point-source validated"`

---

## Task 7: Electrode–tissue impedance → R_a = 1/(4σa) + CPE Bode

**Files:** `crates/valenx-neuro/src/impedance.rs` (create); `src/lib.rs` (+mod).

- [ ] **Step 1 — failing test:**
```rust
#[test]
fn disk_access_resistance_matches_formula() {
    // disk radius a=50 µm, σ=0.2 S/m → R_a = 1/(4σa)
    let z = ElectrodeImpedance::disk(/*a_um*/ 50.0, /*sigma*/ 0.2, Cpe::default());
    let r_a = z.access_resistance_ohm();
    let expect = 1.0 / (4.0 * 0.2 * 50e-6);
    assert!((r_a/expect - 1.0).abs() < 1e-9, "R_a={r_a} expect={expect}");
}
#[test]
fn impedance_is_capacitive_low_resistive_high() {
    let z = ElectrodeImpedance::disk(50.0, 0.2, Cpe::default());
    let lo = z.magnitude_ohm(1.0);        // 1 Hz
    let hi = z.magnitude_ohm(1.0e5);      // 100 kHz
    assert!(lo > hi, "low-f capacitive > high-f resistive: lo={lo} hi={hi}");
    assert!((hi - z.access_resistance_ohm()).abs()/hi < 0.1, "hi-f → R_a plateau");
}
```

- [ ] **Step 2 — run, expect FAIL.** `cargo test -p valenx-neuro impedance`.

- [ ] **Step 3 — implement** `ElectrodeImpedance::disk(a_um, sigma, Cpe)`, `access_resistance_ohm()` = `1/(4σa)`, `Cpe { q, n }` with `Z_CPE = 1/(Q(jω)^n)`, `magnitude_ohm(freq_hz)` = `|R_a + Z_CPE|` (complex). Use a tiny inline complex (re, im) — no new dep.

- [ ] **Step 4 — run, expect PASS** + clippy clean.

- [ ] **Step 5 — commit.** `git commit -am "feat(neuro): electrode impedance — access resistance + CPE"`

---

## Task 8: Workbench wiring (`valenx-app`)

**Files:** modify `crates/valenx-app/Cargo.toml`, `src/lib.rs`, `src/update.rs`; create `src/neuro_workbench.rs`.

- [ ] **Step 1 — dep + module.** Add `valenx-neuro = { path = "../valenx-neuro" }` to `crates/valenx-app/Cargo.toml`. In `src/lib.rs` add `mod neuro_workbench;` and a `show_neuro_workbench: bool` field on the App struct (default `false`), mirroring `show_reactdyn_workbench`.

- [ ] **Step 2 — View toggle + dispatch.** In `src/update.rs` near line 374 add `ui.checkbox(&mut self.show_neuro_workbench, "Neural Interface");` and near line 862 add `crate::neuro_workbench::draw_neuro_workbench(self, ctx);` (guarded internally by the bool), mirroring the reactdyn lines.

- [ ] **Step 3 — failing headless test.** In `neuro_workbench.rs`:
```rust
#[cfg(test)]
mod tests {
    #[test]
    fn default_scene_runs_and_produces_a_trajectory() {
        let scene = valenx_neuro::Scene::single_axon(1.0);
        let r = valenx_neuro::stimulate(&scene, Default::default(), 400.0);
        assert!(r.any_fired(), "default demo scene should recruit");
    }
    #[test]
    fn run_with_zero_axons_is_handled() {
        let scene = valenx_neuro::Scene::bundle(0, 1.0);
        let r = valenx_neuro::stimulate(&scene, Default::default(), 400.0);
        assert_eq!(r.recruited_fraction(), 0.0);   // no panic, empty result
    }
}
```

- [ ] **Step 4 — run, expect FAIL** then implement `draw_neuro_workbench(app, ctx)` (setup panel: tissue σ preset, electrode amplitude/width/polarity, axon count/depth; **Run** spawns the `stimulate` sweep on a background thread like reactdyn; 3-D viewport reusing the reactdyn playback scaffolding — axons colored by V, field heatmap slice, time scrubber; `egui_plot` panels: V(t), recruitment curve, ΔT, |Z|(ω)). Extract draw/run logic into `run_*` free fns for testability (mirror reactdyn).

- [ ] **Step 5 — run, expect PASS.** `cargo test -p valenx-app neuro` + `cargo test -p valenx-neuro` + `cargo clippy -p valenx-app -p valenx-neuro --all-targets -- -D warnings`.

- [ ] **Step 6 — commit.** `git commit -am "feat(neuro): neural-interface workbench — setup, run, 3-D playback, plots"`

---

## Task 9: Docs + RFC status + PR ready

- [ ] **Step 1.** Update RFC 0011 status Draft→Implemented note (Amendments section); add `valenx-neuro` to README Native engines (Reaction dynamics & graphics group → "Neuroengineering") + a Validation-table row per validated module (AP, point-source φ, r² recruitment, bioheat ΔT, R_a).
- [ ] **Step 2.** `cargo test -p valenx-neuro && cargo test -p valenx-app neuro` green; clippy clean; `cargo doc -p valenx-neuro --no-deps` warning-free.
- [ ] **Step 3 — commit + ready the PR.** `git commit -am "docs(neuro): README + RFC 0011 status + validation rows"`; `git push`; `gh pr ready 9`.
- [ ] **Step 4 — résumé bullet upgrade** (out of repo): the "extending toward neural-interface/BCI" line becomes a concrete validated "Neuroengineering" group.

---

## Self-review notes

- **Spec coverage:** all five RFC modules map to Tasks 3 (field), 1–2 (cable), 4 (activating), 6 (bioheat), 7 (impedance); coupling = Task 5; workbench = Task 8; units = Task 0. ✓
- **Ordering caveat (flagged in Task 4):** `threshold_current_ua` (r² test) depends on Task 5's `stimulate`; build the sign test in Task 4, land Task 5, then the r² test. Keep both green by the end of Task 5.
- **Numbers are pinned** (HH constants, σ, k, R_a) so results are reproducible — wrong numbers are worse than a crash.
- **fem API:** Task 3 Step 1 reads the exact `solve_steady_thermal` signature before wrapping — no invented signatures.
- **Tolerances are loose where mesh error is real** (point-source φ and ΔT: ±25%; r² slope: ±0.5) — validating the *law*, not solver-precision matching.
