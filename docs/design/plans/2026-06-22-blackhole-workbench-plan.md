# Black-Hole / Relativity Workbench — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `valenx-app` workbench that surfaces the (already complete + validated) `valenx-relativity` engine — black-hole observables, a ray-traced shadow/lensing image, and geodesic orbits in the 3-D viewport.

**Architecture:** Thin UI over a done engine. The workbench holds `(spacetime, M, a, Q, observer)`, builds a `KerrNewman`, and calls the engine's `Result`-returning observables/thermodynamics/shadow functions; the shadow `ShadowImage` (a per-pixel bool mask) maps to an RGB buffer shown via egui's `ColorImage::from_rgb → load_texture` (the same pattern `render_workbench.rs` uses for the path tracer). No new physics — wiring + presentation only.

**Tech Stack:** Rust · `valenx-relativity` (KerrNewman, observables, thermo, shadow, geodesics) · eframe/egui (panel + texture) · the existing viewport line-draw for geodesics.

**Design/spec:** `docs/design/2026-06-22-valenx-capability-roadmap.md` → Track A.

**Commit identity (every commit):**
```
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "<msg>"
```
After each: `git log -1 --format='%ae' | grep -c forceblue` must print `0`.

**Pattern to mirror for all UI wiring:** `crates/valenx-app/src/render_workbench.rs` (panel + `egui::TextureHandle` + `ColorImage::from_rgb` + `ctx.load_texture`) and any existing `*_workbench.rs` for the show-flag + View-menu + `update.rs` dispatch registration. Follow those verbatim for the boilerplate; this plan gives exact code only for the engine-calling + image-mapping parts.

---

## File Structure

- **Create** `crates/valenx-app/src/blackhole_workbench.rs` — state + panel + the pure compute/render helpers.
- **Modify** `crates/valenx-app/src/lib.rs` — add `mod blackhole_workbench;`, a `show_blackhole_workbench: bool` field on the app, and (mirroring a sibling workbench) a `BlackHoleWorkbenchState` field.
- **Modify** `crates/valenx-app/src/update.rs` — a View-menu toggle + a draw dispatch (mirror an existing workbench's two lines).
- **Modify** `crates/valenx-app/Cargo.toml` — add `valenx-relativity = { path = "../valenx-relativity" }` if not already a dep.

---

## Task 0: Scaffold the workbench + register it

**Files:** create `blackhole_workbench.rs`; modify `lib.rs`, `update.rs`, `Cargo.toml`.

- [ ] **Step 1: Confirm/add the dep** — in `crates/valenx-app/Cargo.toml`, ensure `valenx-relativity = { path = "../valenx-relativity" }` is under `[dependencies]`.

- [ ] **Step 2: Create `blackhole_workbench.rs` with state + an empty panel:**

```rust
//! Black-hole / relativity workbench — surfaces the valenx-relativity engine:
//! observables, a ray-traced shadow image, and geodesic orbits.

use valenx_relativity::spacetimes::KerrNewman;

/// Which named member of the Kerr–Newman family the user is exploring.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SpacetimeKind {
    #[default]
    Schwarzschild,
    Kerr,
    ReissnerNordstrom,
    KerrNewman,
}

/// Workbench UI state.
pub struct BlackHoleWorkbenchState {
    pub kind: SpacetimeKind,
    pub mass: f64,
    pub spin: f64,
    pub charge: f64,
    /// Observer radius (units of M) for the shadow render.
    pub r_obs: f64,
    /// Observer polar angle (degrees) for the shadow render.
    pub theta_obs_deg: f64,
    /// Shadow image resolution (pixels, square).
    pub img_size: usize,
    /// Formatted observables readout (or an error string).
    pub readout: String,
    /// The rendered shadow texture, if any.
    pub texture: Option<egui::TextureHandle>,
}

impl Default for BlackHoleWorkbenchState {
    fn default() -> Self {
        Self {
            kind: SpacetimeKind::Schwarzschild,
            mass: 1.0,
            spin: 0.0,
            charge: 0.0,
            r_obs: 50.0,
            theta_obs_deg: 80.0,
            img_size: 192,
            readout: String::new(),
            texture: None,
        }
    }
}

impl BlackHoleWorkbenchState {
    /// Build the parameter struct, zeroing spin/charge for the simpler members
    /// so the UI can't produce a contradictory hole.
    pub fn hole(&self) -> KerrNewman {
        match self.kind {
            SpacetimeKind::Schwarzschild => KerrNewman { mass: self.mass, spin: 0.0, charge: 0.0 },
            SpacetimeKind::Kerr => KerrNewman { mass: self.mass, spin: self.spin, charge: 0.0 },
            SpacetimeKind::ReissnerNordstrom => {
                KerrNewman { mass: self.mass, spin: 0.0, charge: self.charge }
            }
            SpacetimeKind::KerrNewman => {
                KerrNewman { mass: self.mass, spin: self.spin, charge: self.charge }
            }
        }
    }
}
```

- [ ] **Step 3: Add a minimal panel fn** in the same file (mirror the signature of an existing `show_*_workbench`; the exact `&mut ValenxApp`/`egui::Context` shape comes from a sibling):

```rust
/// Draw the black-hole workbench panel. Mirrors render_workbench.rs's panel shape.
pub fn show_blackhole_workbench(state: &mut BlackHoleWorkbenchState, ui: &mut egui::Ui) {
    ui.heading("Black hole");
    ui.label("Relativity engine front-end (valenx-relativity).");
    // Observables + shadow controls are added in Tasks 1–2.
}
```

- [ ] **Step 4: Register** — in `lib.rs` add `mod blackhole_workbench;`, a `blackhole: blackhole_workbench::BlackHoleWorkbenchState` field (default), and a `show_blackhole_workbench: bool`; in `update.rs` add a View-menu checkbox and a dispatch call to `show_blackhole_workbench(&mut self.blackhole, ui)` — copying the two registration lines from a sibling workbench.

- [ ] **Step 5: Build**

Run: `cargo build -p valenx-app --lib`
Expected: compiles.

- [ ] **Step 6: Commit**

```bash
git add crates/valenx-app/src/blackhole_workbench.rs crates/valenx-app/src/lib.rs crates/valenx-app/src/update.rs crates/valenx-app/Cargo.toml
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(app): scaffold black-hole workbench (task 0)"
```

---

## Task 1: Observables panel

**Files:** modify `blackhole_workbench.rs`.

- [ ] **Step 1: Write the failing test** — append to `blackhole_workbench.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schwarzschild_observables_match_closed_form() {
        let s = BlackHoleWorkbenchState { mass: 1.0, ..Default::default() };
        let text = compute_observables(&s).expect("valid hole");
        // Event horizon r+ = 2M, photon sphere 3M, ISCO 6M for M=1.
        assert!(text.contains("2.000"), "horizon r+=2M:\n{text}");
        assert!(text.contains("3.000"), "photon sphere 3M:\n{text}");
        assert!(text.contains("6.000"), "ISCO 6M:\n{text}");
    }

    #[test]
    fn super_extremal_is_rejected_not_faked() {
        // a² + Q² > M² is a naked singularity — must error, not invent numbers.
        let s = BlackHoleWorkbenchState {
            kind: SpacetimeKind::Kerr,
            mass: 1.0,
            spin: 2.0,
            ..Default::default()
        };
        assert!(compute_observables(&s).is_err());
    }
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p valenx-app --lib blackhole`
Expected: FAIL — `compute_observables` not defined.

- [ ] **Step 3: Implement `compute_observables`** (add to `blackhole_workbench.rs`; confirm each function name against `valenx_relativity::observables` — they are re-exported from the crate root):

```rust
use std::fmt::Write as _;
use valenx_relativity::{
    gravitational_redshift, horizons, isco, photon_sphere, shadow_radius, thermodynamics,
    OrbitSense, RelativityError,
};

/// Format the full observables readout, or return the engine's error verbatim
/// (a super-extremal hole / non-positive mass must surface honestly, not fake
/// a number).
pub fn compute_observables(s: &BlackHoleWorkbenchState) -> Result<String, RelativityError> {
    let bh = s.hole();
    let h = horizons(&bh)?;
    let ph = photon_sphere(&bh, OrbitSense::Prograde)?;
    let risco = isco(&bh, OrbitSense::Prograde)?;
    let rsh = shadow_radius(&bh)?;
    let thermo = thermodynamics(&bh)?;
    // Redshift from the ISCO out to infinity, as a representative figure.
    let z = gravitational_redshift(&bh, risco)?;

    let mut out = String::new();
    let _ = writeln!(out, "mass M            : {:.3}", bh.mass);
    let _ = writeln!(out, "spin a            : {:.3}", bh.spin);
    let _ = writeln!(out, "charge Q          : {:.3}", bh.charge);
    let _ = writeln!(out, "event horizon r+  : {:.3} M", h.outer);
    let _ = writeln!(out, "inner horizon r-  : {:.3} M", h.inner);
    let _ = writeln!(out, "photon sphere     : {:.3} M", ph);
    let _ = writeln!(out, "ISCO (prograde)   : {:.3} M", risco);
    let _ = writeln!(out, "shadow radius     : {:.3} M", rsh);
    let _ = writeln!(out, "redshift z @ISCO  : {:.4}", z);
    let _ = writeln!(out, "Hawking T_H       : {:.4e}", thermo.hawking_temperature);
    let _ = writeln!(out, "entropy S         : {:.4e}", thermo.entropy);
    let _ = writeln!(out, "horizon Ω_H       : {:.4}", thermo.horizon_angular_velocity);
    Ok(out)
}
```

> NOTE: confirm the exact argument lists of `photon_sphere` / `isco` / `shadow_radius` / `gravitational_redshift` in `crates/valenx-relativity/src/observables.rs` (some take an `OrbitSense`, some take a radius). The shapes above match the documented closed forms; adjust the calls to the real signatures — the engine functions themselves are correct + tested.

- [ ] **Step 4: Wire into the panel** — in `show_blackhole_workbench`, add the spacetime picker + `DragValue`s for M/a/Q, and a "Compute" button that sets `state.readout = match compute_observables(state) { Ok(t) => t, Err(e) => format!("⚠ {e}") };`, then show `state.readout` monospaced. (Disable the spin field for Schwarzschild/RN and charge for Schwarzschild/Kerr — cosmetic.)

- [ ] **Step 5: Run the tests**

Run: `cargo test -p valenx-app --lib blackhole`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/valenx-app/src/blackhole_workbench.rs
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(app): black-hole observables panel — horizons/ISCO/photon-sphere/shadow/Hawking (task 1)"
```

---

## Task 2: Shadow / lensing image

**Files:** modify `blackhole_workbench.rs`.

- [ ] **Step 1: Write the failing test** (the pure bool-mask → RGB mapping, no GPU/UI):

```rust
    #[test]
    fn shadow_mask_maps_to_black_inside_and_nonblack_outside() {
        use valenx_relativity::shadow::ShadowImage;
        let img = ShadowImage {
            width: 2,
            height: 1,
            half_extent: 8.0,
            shadow: vec![true, false], // pixel 0 captured, pixel 1 escapes
        };
        let rgb = shadow_to_rgb(&img);
        assert_eq!(rgb.len(), 2 * 1 * 3);
        assert_eq!(&rgb[0..3], &[0, 0, 0], "shadow pixel is black");
        assert!(rgb[3..6].iter().any(|&c| c > 0), "sky pixel is not black");
    }
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p valenx-app --lib blackhole`
Expected: FAIL — `shadow_to_rgb` not defined.

- [ ] **Step 3: Implement the mapping + render trigger:**

```rust
use valenx_relativity::shadow::{render_shadow, ShadowImage};

/// Map a shadow mask to an RGB8 buffer (row-major, `w*h*3` bytes): the captured
/// region is black; escaping pixels get a simple radial "sky" tint so the
/// photon ring reads as the boundary. Good enough for the panel; richer shading
/// (accretion disk) is a later phase.
pub fn shadow_to_rgb(img: &ShadowImage) -> Vec<u8> {
    let (w, h) = (img.width, img.height);
    let mut rgb = vec![0u8; w * h * 3];
    for row in 0..h {
        for col in 0..w {
            let i = (row * w + col) * 3;
            if img.is_shadow(col, row) {
                continue; // black
            }
            // Sky tint by distance from image centre (purely cosmetic).
            let (cx, cy) = (w as f64 / 2.0, h as f64 / 2.0);
            let d = (((col as f64 - cx).powi(2) + (row as f64 - cy).powi(2)).sqrt()
                / (cx.max(cy)))
            .clamp(0.0, 1.0);
            rgb[i] = (40.0 + 120.0 * d) as u8;
            rgb[i + 1] = (60.0 + 120.0 * d) as u8;
            rgb[i + 2] = (120.0 + 135.0 * d) as u8;
        }
    }
    rgb
}

/// Render the shadow for the current state into an egui texture.
pub fn render_shadow_texture(
    s: &BlackHoleWorkbenchState,
    ctx: &egui::Context,
) -> Result<egui::TextureHandle, RelativityError> {
    let bh = s.hole();
    let n = s.img_size.clamp(16, 512);
    // Confirm render_shadow's signature in shadow.rs — it takes the hole, the
    // observer (r_obs, theta_obs in radians), the image-plane half-extent, and
    // the pixel dimensions, and returns Result<ShadowImage>.
    let img: ShadowImage = render_shadow(
        &bh,
        s.r_obs,
        s.theta_obs_deg.to_radians(),
        12.0, // image-plane half-extent in M (≈ a few × shadow radius)
        n,
        n,
    )?;
    let rgb = shadow_to_rgb(&img);
    let color = egui::ColorImage::from_rgb([img.width, img.height], &rgb);
    Ok(ctx.load_texture("blackhole_shadow", color, egui::TextureOptions::LINEAR))
}
```

- [ ] **Step 4: Wire into the panel** — add a "Render shadow" button: `state.texture = match render_shadow_texture(state, ui.ctx()) { Ok(t) => Some(t), Err(e) => { state.readout = format!("⚠ {e}"); None } };` then, if `Some(tex)`, show it scaled-to-fit (copy the texture-display block from `render_workbench.rs`). Add preset buttons that set `(mass, spin)` for **Sgr A\*** and **M87\*** (use `M = 1`, `a ≈ 0.5–0.9` as illustrative spins — they're shown in geometrized units). Cap `img_size` ≤ 256 in the UI so a click can't hang the frame; for larger sizes, move the render onto the background-job thread (the run pipeline) in a later pass.

- [ ] **Step 5: Run the tests**

Run: `cargo test -p valenx-app --lib blackhole`
Expected: PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/valenx-app/src/blackhole_workbench.rs
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(app): black-hole shadow/lensing render into the image panel (task 2)"
```

---

## Task 3: Geodesic orbits in the 3-D viewport

**Files:** modify `blackhole_workbench.rs` (+ whatever viewport line-draw entry the app exposes).

- [ ] **Step 1:** Add a pure helper that integrates an equatorial orbit and returns a polyline of `[x,y,z]` points (convert Boyer–Lindquist `(r, φ)` at `θ = π/2` to Cartesian):

```rust
use valenx_relativity::{equatorial_state, integrate_geodesic, GeodesicOptions, Kind};

/// Integrate one equatorial geodesic and return its Cartesian polyline (M units).
/// `kind` selects a photon (null) or massive (timelike) orbit.
pub fn orbit_polyline(
    s: &BlackHoleWorkbenchState,
    r0: f64,
    kind: Kind,
) -> Result<Vec<[f64; 3]>, RelativityError> {
    let bh = s.hole();
    let state0 = equatorial_state(&bh, r0, kind)?;        // confirm signature
    let traj = integrate_geodesic(&bh, state0, &GeodesicOptions::default())?;
    Ok(traj
        .states
        .iter()
        .map(|st| {
            let r = st.position[1];
            let phi = st.position[3];
            [r * phi.cos(), r * phi.sin(), 0.0]
        })
        .collect())
}
```

> Confirm `equatorial_state` / `Trajectory` field names against `geodesics.rs` (the module is re-exported from the crate root). The conversion above assumes equatorial `θ = π/2`.

- [ ] **Step 2:** A unit test that a bound massive orbit at `r0 = 10M` (Schwarzschild) returns a non-empty polyline whose radii stay finite and roughly bounded:

```rust
    #[test]
    fn massive_orbit_returns_a_bounded_polyline() {
        let s = BlackHoleWorkbenchState { mass: 1.0, ..Default::default() };
        let line = orbit_polyline(&s, 10.0, valenx_relativity::Kind::Timelike).unwrap();
        assert!(line.len() > 10);
        assert!(line.iter().all(|p| p.iter().all(|c| c.is_finite())));
    }
```

- [ ] **Step 3:** Run → fail → implement → pass (`cargo test -p valenx-app --lib blackhole`).

- [ ] **Step 4:** In the panel, an "Add orbit" control that pushes the polyline into the 3-D viewport's debug-line list (reuse the existing viewport polyline-draw path the app already uses for other overlays). Draw the event-horizon sphere (radius `horizons().outer`) as a shaded sphere for context.

- [ ] **Step 5: Commit**

```bash
git add crates/valenx-app/src/blackhole_workbench.rs
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(app): trace geodesic orbits around the hole in the 3-D viewport (task 3)"
```

---

## Task 4: AI-drivable command + time-dilation extra

- [ ] **Step 1:** Add an agent-bridge command `BlackHole { spacetime, mass, spin, charge }` that sets the state + runs `compute_observables` (mirror an existing agent command in the bridge dispatch). Test by driving it through the bridge as the other workbenches are.
- [ ] **Step 2:** Add a small "clock at r vs ∞" readout: for the current hole, `dτ/dt = √(−g_tt)` at a user radius — use `gravitational_redshift` (already called) to show time-dilation factor. One label, no new physics.
- [ ] **Step 3:** `cargo clippy -p valenx-app --lib -- -D warnings` clean; `cargo fmt`.
- [ ] **Step 4: Commit**

```bash
git add -A crates/valenx-app/src
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(app): black-hole agent command + time-dilation readout; clippy/fmt clean (task 4)"
```

---

## Done — what exists

A Black-Hole workbench: pick a spacetime + (M, a, Q); read horizons / ergosphere / photon-sphere / ISCO / shadow / redshift / Hawking thermodynamics (super-extremal holes rejected honestly); render the ray-traced shadow/lensing image; trace geodesic orbits in the 3-D viewport; AI-drivable. All on the existing validated engine — no new physics.

## Execution
Deferred. When building: subagent-driven (recommended) or inline. The engine functions are tested already; the risk surface here is only the egui wiring + confirming the exact `observables`/`shadow`/`geodesics` signatures (noted inline) — so a quick checkpoint after Task 1 (it proves the engine calls compile + return the right numbers) de-risks the rest.
