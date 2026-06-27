//! In-house **Turing morphogenesis** reaction–diffusion model (Gray–Scott).
//!
//! This is the pure, headless compute core behind the Morphogenesis workbench.
//! It evolves two virtual "morphogens" `U` and `V` on a rectangular grid under
//! the Gray–Scott reaction–diffusion system — the same class of model Alan
//! Turing proposed in 1952 ("The Chemical Basis of Morphogenesis") as the
//! mechanism by which developing organisms break symmetry into spots, stripes,
//! and other organic patterns. Different feed/kill rates settle into visibly
//! different biological-looking régimes (spots, coral, mitosis, mazes).
//!
//! ## Model
//!
//! Per explicit Euler step (`dt` baked at `1.0`), with `lap` the 5-point
//! Laplacian on a **toroidal** (wrap-around) grid:
//!
//! ```text
//! U' = U + (Du·lap(U) − U·V·V + F·(1 − U)) · dt
//! V' = V + (Dv·lap(V) + U·V·V − (F + k)·V) · dt
//! ```
//!
//! `Du`/`Dv` are diffusion rates (`U` diffuses faster than `V`), `F` is the
//! feed rate (replenishes `U`), and `k` the kill rate (removes `V`). Values are
//! clamped to `[0, 1]` and guarded against non-finite drift every step so the
//! field can be rendered directly as a height/colour map.
//!
//! The struct is deliberately framework-free (no `egui`, no app state) so it is
//! cheap to unit-test and reuse. The workbench
//! ([`crate::morphogenesis_workbench`]) owns the GUI, the real-time stepping
//! loop, and the 3-D rendering of the `V` field.

/// A live Gray–Scott reaction–diffusion field on a `w × h` toroidal grid.
///
/// `u` and `v` are row-major (`idx = y·w + x`) concentration buffers in
/// `[0, 1]`. `step` advances the system with explicit Euler integration; the
/// reaction parameters (`du`, `dv`, `feed`, `kill`) select the pattern régime
/// (see [`Morphogenesis::preset`]).
#[derive(Clone, Debug)]
pub struct Morphogenesis {
    /// Grid width in cells.
    pub w: usize,
    /// Grid height in cells.
    pub h: usize,
    /// Morphogen `U` concentration, row-major, length `w·h`, in `[0, 1]`.
    pub u: Vec<f32>,
    /// Morphogen `V` concentration, row-major, length `w·h`, in `[0, 1]`.
    pub v: Vec<f32>,
    /// Diffusion rate of `U` (`Du`). `U` diffuses faster than `V`.
    pub du: f32,
    /// Diffusion rate of `V` (`Dv`).
    pub dv: f32,
    /// Feed rate `F` (replenishes `U`).
    pub feed: f32,
    /// Kill rate `k` (removes `V`).
    pub kill: f32,
    /// Total number of explicit-Euler steps taken since the last reseed.
    pub steps: u64,
}

impl Morphogenesis {
    /// Build a fresh field on a `w × h` grid (each dimension clamped to a sane
    /// `[16, 512]` range), seeded with `U = 1.0` everywhere and a small square
    /// blob of `V` in the centre — the symmetry-breaking germ the pattern grows
    /// from. Default parameters are the classic "coral"-ish balance
    /// (`Du = 0.16`, `Dv = 0.08`, the `mitosis` feed/kill pair) which evolves a
    /// clearly organic pattern.
    pub fn new(w: usize, h: usize) -> Self {
        let w = w.clamp(16, 512);
        let h = h.clamp(16, 512);
        let (feed, kill) = Self::preset("mitosis");
        let mut m = Self {
            w,
            h,
            u: vec![1.0; w * h],
            v: vec![0.0; w * h],
            du: 0.16,
            dv: 0.08,
            feed,
            kill,
            steps: 0,
        };
        m.seed_center();
        m
    }

    /// Mean `V` concentration over the whole grid (the pattern "mass"); `0.0`
    /// for an empty grid. Handy as a one-number progress readout.
    pub fn mean_v(&self) -> f32 {
        if self.v.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.v.iter().map(|&x| x as f64).sum();
        (sum / self.v.len() as f64) as f32
    }

    /// `(min, max)` of the `V` field, used by the renderer to normalise the
    /// height/colour map. Returns `(0.0, 1.0)` for an empty grid.
    pub fn field_minmax(&self) -> (f32, f32) {
        if self.v.is_empty() {
            return (0.0, 1.0);
        }
        let mut lo = f32::INFINITY;
        let mut hi = f32::NEG_INFINITY;
        for &x in &self.v {
            if x.is_finite() {
                lo = lo.min(x);
                hi = hi.max(x);
            }
        }
        if !(lo.is_finite() && hi.is_finite()) {
            return (0.0, 1.0);
        }
        (lo, hi)
    }

    /// Read the `V` concentration at `(x, y)` (row-major), clamped to the grid
    /// bounds. The field the workbench renders as a 3-D surface.
    #[inline]
    pub fn v_at(&self, x: usize, y: usize) -> f32 {
        if self.w == 0 || self.h == 0 {
            return 0.0;
        }
        let x = x.min(self.w - 1);
        let y = y.min(self.h - 1);
        self.v[y * self.w + x]
    }

    /// Seed `U = 1`, `V = 0` everywhere, then drop a small centred square of
    /// `V = 0.5` / `U = 0.25` — the germ the Turing pattern grows from.
    fn seed_center(&mut self) {
        for x in self.u.iter_mut() {
            *x = 1.0;
        }
        for x in self.v.iter_mut() {
            *x = 0.0;
        }
        // A blob roughly 1/8 of the smaller dimension, at least a few cells.
        let r = (self.w.min(self.h) / 8).max(3);
        let cx = self.w / 2;
        let cy = self.h / 2;
        for dy in 0..(2 * r) {
            for dx in 0..(2 * r) {
                let x = (cx + dx).wrapping_sub(r);
                let y = (cy + dy).wrapping_sub(r);
                if x < self.w && y < self.h {
                    let i = y * self.w + x;
                    self.v[i] = 0.5;
                    self.u[i] = 0.25;
                }
            }
        }
    }

    /// Reset the field to the centred seed and zero the step counter (keeps the
    /// current grid size and reaction parameters). Use this to restart growth
    /// after changing `feed`/`kill` to watch a different régime emerge.
    pub fn reseed(&mut self) {
        self.seed_center();
        self.steps = 0;
    }

    /// Resize the grid to `w × h` (clamped to `[16, 512]`) and reseed. Safe for
    /// any input; reallocates the buffers and restarts the pattern.
    pub fn resize(&mut self, w: usize, h: usize) {
        let w = w.clamp(16, 512);
        let h = h.clamp(16, 512);
        if w == self.w && h == self.h {
            self.reseed();
            return;
        }
        self.w = w;
        self.h = h;
        self.u = vec![1.0; w * h];
        self.v = vec![0.0; w * h];
        self.reseed();
    }

    /// The named Gray–Scott `(feed, kill)` parameter pairs that settle into
    /// distinct biological-looking patterns. Unknown names fall back to the
    /// `"mitosis"` pair (never panics):
    ///
    /// | name      | F       | k       | look                         |
    /// |-----------|---------|---------|------------------------------|
    /// | `spots`   | 0.035   | 0.065   | isolated dots                |
    /// | `coral`   | 0.055   | 0.062   | branching coral / fingerprint|
    /// | `mitosis` | 0.0367  | 0.0649  | self-replicating cells       |
    /// | `maze`    | 0.029   | 0.057   | labyrinthine stripes         |
    pub fn preset(name: &str) -> (f32, f32) {
        match name {
            "spots" => (0.035, 0.065),
            "coral" => (0.055, 0.062),
            "mitosis" => (0.0367, 0.0649),
            "maze" => (0.029, 0.057),
            _ => (0.0367, 0.0649),
        }
    }

    /// Apply one of the named [`preset`](Self::preset) régimes by setting
    /// `feed`/`kill` and reseeding so the new pattern grows from scratch.
    pub fn apply_preset(&mut self, name: &str) {
        let (f, k) = Self::preset(name);
        self.feed = f;
        self.kill = k;
        self.reseed();
    }

    /// 5-point Laplacian of buffer `b` at `(x, y)` with **toroidal** wrap
    /// (so the pattern tiles seamlessly and has no hard boundary):
    /// `b[x-1] + b[x+1] + b[y-1] + b[y+1] − 4·b[x,y]`.
    #[inline]
    fn laplacian(b: &[f32], w: usize, h: usize, x: usize, y: usize) -> f32 {
        let xm = if x == 0 { w - 1 } else { x - 1 };
        let xp = if x + 1 == w { 0 } else { x + 1 };
        let ym = if y == 0 { h - 1 } else { y - 1 };
        let yp = if y + 1 == h { 0 } else { y + 1 };
        let c = b[y * w + x];
        b[y * w + xm] + b[y * w + xp] + b[ym * w + x] + b[yp * w + x] - 4.0 * c
    }

    /// Advance the system by `substeps` explicit-Euler steps (`dt = 1`). A
    /// `substeps` of `0` is a **no-op**. Each step computes the Gray–Scott
    /// reaction–diffusion update into scratch buffers, then swaps them in;
    /// every written value is clamped to `[0, 1]` and any non-finite result is
    /// coerced to `0.0`, so the field stays renderable no matter the parameters.
    pub fn step(&mut self, substeps: u32) {
        if substeps == 0 || self.w == 0 || self.h == 0 {
            return;
        }
        let (w, h) = (self.w, self.h);
        let mut nu = vec![0.0f32; w * h];
        let mut nv = vec![0.0f32; w * h];
        for _ in 0..substeps {
            for y in 0..h {
                for x in 0..w {
                    let i = y * w + x;
                    let u = self.u[i];
                    let v = self.v[i];
                    let lu = Self::laplacian(&self.u, w, h, x, y);
                    let lv = Self::laplacian(&self.v, w, h, x, y);
                    let uvv = u * v * v;
                    // dt baked at 1.0.
                    let mut un = u + (self.du * lu - uvv + self.feed * (1.0 - u));
                    let mut vn = v + (self.dv * lv + uvv - (self.feed + self.kill) * v);
                    if !un.is_finite() {
                        un = 0.0;
                    }
                    if !vn.is_finite() {
                        vn = 0.0;
                    }
                    nu[i] = un.clamp(0.0, 1.0);
                    nv[i] = vn.clamp(0.0, 1.0);
                }
            }
            std::mem::swap(&mut self.u, &mut nu);
            std::mem::swap(&mut self.v, &mut nv);
            self.steps += 1;
        }
    }
}

impl Default for Morphogenesis {
    fn default() -> Self {
        Self::new(96, 96)
    }
}

// ---------------------------------------------------------------------------
// Tests (pure)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_seeds_u_full_and_a_v_blob() {
        let m = Morphogenesis::new(64, 64);
        assert_eq!(m.u.len(), 64 * 64);
        assert_eq!(m.v.len(), 64 * 64);
        // U starts at 1 in the corners (away from the seed blob).
        assert!((m.v_at(0, 0)).abs() < 1e-6, "corner V should start at 0");
        assert!((m.u[0] - 1.0).abs() < 1e-6, "corner U should start at 1");
        // The centre has a non-zero V seed.
        assert!(m.v_at(32, 32) > 0.0, "centre V seed should be > 0");
        assert_eq!(m.steps, 0);
    }

    #[test]
    fn new_clamps_extreme_grid_sizes() {
        let tiny = Morphogenesis::new(1, 1);
        assert!(tiny.w >= 16 && tiny.h >= 16);
        let huge = Morphogenesis::new(100_000, 100_000);
        assert!(huge.w <= 512 && huge.h <= 512);
        assert_eq!(huge.u.len(), huge.w * huge.h);
    }

    #[test]
    fn seeded_field_evolves_mass_changes() {
        // After stepping, the total V "mass" must differ from the seed — the
        // reaction is doing something (the core liveness check).
        let mut m = Morphogenesis::new(80, 80);
        let before: f64 = m.v.iter().map(|&x| x as f64).sum();
        m.step(200);
        let after: f64 = m.v.iter().map(|&x| x as f64).sum();
        assert!(m.steps == 200);
        assert!(
            (after - before).abs() > 1e-3,
            "V mass should change after stepping (before={before}, after={after})"
        );
    }

    #[test]
    fn all_values_stay_finite_and_in_unit_range() {
        let mut m = Morphogenesis::new(64, 64);
        m.step(300);
        for (&u, &v) in m.u.iter().zip(m.v.iter()) {
            assert!(u.is_finite() && v.is_finite(), "values must stay finite");
            assert!((0.0..=1.0).contains(&u), "U out of [0,1]: {u}");
            assert!((0.0..=1.0).contains(&v), "V out of [0,1]: {v}");
        }
    }

    #[test]
    fn fixed_seed_is_reproducible() {
        // The model is fully deterministic: two identically-built+stepped
        // fields must match bit-for-bit.
        let mut a = Morphogenesis::new(48, 48);
        let mut b = Morphogenesis::new(48, 48);
        a.step(150);
        b.step(150);
        assert_eq!(a.v, b.v, "deterministic V fields must be identical");
        assert_eq!(a.u, b.u, "deterministic U fields must be identical");
    }

    #[test]
    fn presets_return_expected_pairs() {
        assert_eq!(Morphogenesis::preset("spots"), (0.035, 0.065));
        assert_eq!(Morphogenesis::preset("coral"), (0.055, 0.062));
        assert_eq!(Morphogenesis::preset("mitosis"), (0.0367, 0.0649));
        assert_eq!(Morphogenesis::preset("maze"), (0.029, 0.057));
        // Unknown falls back to the mitosis pair (no panic).
        assert_eq!(Morphogenesis::preset("nope"), (0.0367, 0.0649));
    }

    #[test]
    fn step_zero_substeps_is_a_noop() {
        let mut m = Morphogenesis::new(32, 32);
        let v0 = m.v.clone();
        m.step(0);
        assert_eq!(m.v, v0, "step(0) must not change the field");
        assert_eq!(m.steps, 0, "step(0) must not advance the counter");
    }

    #[test]
    fn resize_is_safe_and_reseeds() {
        let mut m = Morphogenesis::new(64, 64);
        m.step(50);
        m.resize(40, 100);
        assert_eq!(m.w, 40);
        assert_eq!(m.h, 100);
        assert_eq!(m.u.len(), 40 * 100);
        assert_eq!(m.v.len(), 40 * 100);
        assert_eq!(m.steps, 0, "resize reseeds and resets the step counter");
        // Extreme resize is clamped, not panicking.
        m.resize(1, 1);
        assert!(m.w >= 16 && m.h >= 16);
    }

    #[test]
    fn apply_preset_sets_params_and_reseeds() {
        let mut m = Morphogenesis::new(48, 48);
        m.step(20);
        m.apply_preset("maze");
        assert_eq!((m.feed, m.kill), (0.029, 0.057));
        assert_eq!(m.steps, 0, "applying a preset reseeds");
    }

    #[test]
    fn minmax_and_mean_are_well_formed() {
        let mut m = Morphogenesis::new(48, 48);
        m.step(100);
        let (lo, hi) = m.field_minmax();
        assert!(lo.is_finite() && hi.is_finite() && lo <= hi);
        let mean = m.mean_v();
        assert!((0.0..=1.0).contains(&mean), "mean V in [0,1]: {mean}");
    }
}
