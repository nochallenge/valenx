//! Panel 8 — **Macromolecular Structure** (`valenx-biostruct`).
//!
//! Load a macromolecular structure (PDB or mmCIF), compute a per-chain
//! secondary-structure / composition summary, classify the
//! Ramachandran φ/ψ distribution, and Kabsch-superpose two structures'
//! Cα atoms for an RMSD — all native `valenx-biostruct` calls.

use eframe::egui;
use nalgebra::Point3;

use valenx_biostruct::analyze::StructureReport;
use valenx_biostruct::geometry::ramachandran::summarize as rama_summarize;
use valenx_biostruct::io::read_structure;
use valenx_biostruct::structure::Structure;
use valenx_biostruct::superpose::{kabsch, rmsd};

use super::common;
use super::molecule_view::{self, ViewAtom, ViewMolecule};
use crate::molviz::{self, AtomAttr, BackbonePoint, ColorScheme, MolvizParams, Representation};
use crate::ValenxApp;

/// MD coordinates are nanometres; the molecular viewer (biostruct / cheminf)
/// works in ångström. A loaded [`valenx_md::io::trajectory::Trajectory`] is
/// scaled by this factor on attach so every frame shares the viewer length
/// scale, exactly as [`ViewMolecule::from_md_system`] does.
const ANGSTROM_PER_NM: f32 = 10.0;

/// MD-trajectory playback state for the molecular viewer — **roadmap
/// feature**: animate the currently-shown structure across a series of
/// coordinate frames (the MolecularNodes-style playback).
///
/// The structure is parsed **once** into a base [`ViewMolecule`] (atoms in the
/// exact [`ViewMolecule::from_biostruct`] order) plus its per-atom [`AtomAttr`]
/// and Cα backbone, so each frame only overwrites the atom positions and
/// re-meshes — the connectivity/colour metadata is reused. Frames are stored in
/// **ångström** (`[f32; 3]` per atom), so applying a frame is a unit-correct
/// position copy regardless of whether it came from the synthetic generator or
/// a loaded `valenx-md` trajectory (which is converted from nm on attach).
///
/// Fail-loud, never panic: a frame whose atom count differs from the base
/// structure sets [`note`](Self::note) and is skipped; an empty trajectory is a
/// no-op (nothing to attach / play).
#[derive(Clone, Debug, Default)]
pub(crate) struct TrajectoryPlayback {
    /// The base molecule (atoms in `from_biostruct` order) the frames animate.
    /// `None` until a trajectory is attached.
    base: Option<ViewMolecule>,
    /// Per-atom colour attributes for the base molecule, in lockstep with
    /// `base.atoms` — reused for every frame's coloured rebuild.
    attrs: Vec<AtomAttr>,
    /// Cα backbone control points (DSSP-tagged) of the base structure — reused
    /// for the cartoon / ribbon representations. (Static across frames: only
    /// atom positions move per frame, the secondary-structure track is from the
    /// reference structure.)
    backbone: Vec<BackbonePoint>,
    /// The coordinate frames, each `base.atoms.len()` positions in **ångström**.
    frames: Vec<Vec<[f32; 3]>>,
    /// Currently-displayed frame index (clamped to `0..frames.len()`).
    frame: usize,
    /// Whether playback is advancing each repaint.
    playing: bool,
    /// Frames advanced per second while playing.
    speed: f32,
    /// Fractional-frame accumulator so a non-integer `speed × dt` advances
    /// smoothly (carries the leftover between repaints).
    accum: f32,
    /// A human-readable label for where the trajectory came from (synthetic /
    /// loaded), shown next to the controls.
    source: String,
    /// In-panel note (atom-count mismatch, etc.) — surfaced instead of a panic.
    note: Option<String>,
    /// Diagnostics of the last in-house ENM-MD run, when the attached
    /// trajectory came from [`Self::attach_enm_md`] — surfaced in the panel and
    /// the agent readout (frames / springs / temperature / RMSD / energy). The
    /// synthetic-wobble and loaded-file paths leave this `None`.
    md: Option<EnmMd>,
}

impl TrajectoryPlayback {
    /// Number of attached frames (0 when nothing is attached).
    pub(crate) fn n_frames(&self) -> usize {
        self.frames.len()
    }

    /// Whether a trajectory is attached and has at least one frame.
    pub(crate) fn is_attached(&self) -> bool {
        self.base.is_some() && !self.frames.is_empty()
    }

    /// Drop any attached trajectory and reset playback.
    fn clear(&mut self) {
        *self = TrajectoryPlayback::default();
    }

    /// Attach an explicit set of ångström frames to a base molecule.
    ///
    /// Fail-loud: an empty `frames` list is a no-op (nothing is attached); any
    /// frame whose atom count differs from `base.atoms` sets [`note`] and the
    /// whole attach is rejected (so playback never indexes a short frame).
    fn attach(
        &mut self,
        base: ViewMolecule,
        attrs: Vec<AtomAttr>,
        backbone: Vec<BackbonePoint>,
        frames: Vec<Vec<[f32; 3]>>,
        source: impl Into<String>,
    ) {
        if frames.is_empty() {
            // Empty trajectory -> no-op, leave any prior attachment untouched.
            return;
        }
        let n = base.atoms.len();
        if let Some(bad) = frames.iter().position(|f| f.len() != n) {
            self.note = Some(format!(
                "trajectory not attached: frame {} has {} atoms but the structure \
                 has {} — atom counts must match",
                bad,
                frames[bad].len(),
                n,
            ));
            return;
        }
        self.base = Some(base);
        self.attrs = attrs;
        self.backbone = backbone;
        self.frames = frames;
        self.frame = 0;
        self.playing = false;
        self.accum = 0.0;
        self.source = source.into();
        self.note = None;
        self.md = None; // cleared here; an MD attach sets it afterwards
    }

    /// Build the base molecule + frames from a parsed [`Structure`] by
    /// generating a small **synthetic wobble** — each atom oscillates about its
    /// reference position with a per-atom phase, so there is a real, deterministic
    /// trajectory to play with no external file. `n_frames` frames over one full
    /// period; `amplitude` is the ångström displacement.
    fn attach_synthetic(&mut self, s: &Structure, n_frames: usize, amplitude: f32) {
        let base = ViewMolecule::from_biostruct(s);
        if base.atoms.is_empty() {
            self.note = Some("no atoms to animate".to_string());
            return;
        }
        let attrs = structure_atom_attrs(s);
        let backbone = ca_backbone(s);
        let frames = synthetic_wobble_frames(&base, n_frames.max(2), amplitude);
        self.attach(base, attrs, backbone, frames, "synthetic wobble");
    }

    /// Attach a loaded `valenx-md` [`Trajectory`](valenx_md::io::trajectory::Trajectory)
    /// to a parsed [`Structure`]. The MD positions (nm) are converted to ångström.
    ///
    /// Fail-loud: an empty trajectory is a no-op; an atom-count mismatch between
    /// the structure and the trajectory's per-frame atom count sets [`note`] and
    /// attaches nothing.
    fn attach_md(&mut self, s: &Structure, traj: &valenx_md::io::trajectory::Trajectory) {
        let base = ViewMolecule::from_biostruct(s);
        if base.atoms.is_empty() {
            self.note = Some("no atoms to animate".to_string());
            return;
        }
        if traj.is_empty() {
            return; // empty trajectory -> no-op
        }
        if traj.n_atoms() != base.atoms.len() {
            self.note = Some(format!(
                "trajectory not attached: it has {} atoms per frame but the \
                 structure has {} — atom counts must match",
                traj.n_atoms(),
                base.atoms.len(),
            ));
            return;
        }
        let attrs = structure_atom_attrs(s);
        let backbone = ca_backbone(s);
        let frames: Vec<Vec<[f32; 3]>> = (0..traj.len())
            .filter_map(|i| traj.frame(i))
            .map(|frame| {
                frame
                    .iter()
                    .map(|p| {
                        [
                            p.x as f32 * ANGSTROM_PER_NM,
                            p.y as f32 * ANGSTROM_PER_NM,
                            p.z as f32 * ANGSTROM_PER_NM,
                        ]
                    })
                    .collect()
            })
            .collect();
        self.attach(base, attrs, backbone, frames, "loaded valenx-md trajectory");
    }

    /// Build the base molecule + frames from a parsed [`Structure`] by running
    /// the in-house **Elastic Network Model MD** ([`enm_md`]) — a real,
    /// physically-motivated thermal-vibration trajectory of the protein (the
    /// large-scale collective "breathing" of its low-frequency normal modes),
    /// not a synthetic per-atom wobble.
    ///
    /// Deterministic for a fixed seed; stores the run diagnostics
    /// ([`EnmMd`]) so the panel / agent can report frames / springs /
    /// temperature / RMSD / energy. Fail-loud: a structure with no atoms (or
    /// fewer than two) sets [`note`] and attaches nothing.
    fn attach_enm_md(&mut self, s: &Structure, params: EnmParams) {
        let base = ViewMolecule::from_biostruct(s);
        if base.atoms.len() < 2 {
            self.note = Some("need ≥ 2 atoms to run MD".to_string());
            return;
        }
        let attrs = structure_atom_attrs(s);
        let backbone = ca_backbone(s);
        let md = enm_md(&base, params);
        if md.frames.is_empty() {
            self.note = Some("ENM MD produced no frames".to_string());
            return;
        }
        let summary = md.summary();
        self.attach(base, attrs, backbone, md.frames.clone(), summary);
        self.md = Some(md);
    }

    /// Diagnostics of the last in-house ENM-MD run, if the attached trajectory
    /// came from [`Self::attach_enm_md`].
    pub(crate) fn md(&self) -> Option<&EnmMd> {
        self.md.as_ref()
    }

    /// The base molecule with the positions of frame `index` applied, or `None`
    /// if nothing is attached / `index` is out of range. The returned molecule
    /// has **freshly detected bonds** for that frame's geometry, so a trajectory
    /// in which atoms move apart/together re-bonds correctly each frame.
    fn molecule_at(&self, index: usize) -> Option<ViewMolecule> {
        let base = self.base.as_ref()?;
        let frame = self.frames.get(index)?;
        if frame.len() != base.atoms.len() {
            return None; // guarded at attach, but never index a short frame
        }
        let atoms: Vec<ViewAtom> = base
            .atoms
            .iter()
            .zip(frame)
            .map(|(a, &pos)| ViewAtom::new(pos, a.element.clone()))
            .collect();
        let bonds = molecule_view::detect_bonds(&atoms);
        Some(ViewMolecule { atoms, bonds })
    }
}

/// Generate a deterministic synthetic wobble: `n_frames` frames in which every
/// atom oscillates about its reference position. The displacement is
/// `amplitude · sin(t) · dir(i)` where `t` sweeps one full period over the
/// frames and `dir(i)` is a fixed per-atom direction (decorrelated by the atom
/// index) — so the cloud *breathes* (atoms move along different axes) rather
/// than translating rigidly, giving a visibly animated but reproducible
/// trajectory with no external data.
///
/// Because the envelope is `sin(t)` (not `sin(t + phase)`), **frame 0 is the
/// reference structure exactly** (`sin 0 == 0`), so attaching a synthetic
/// trajectory shows the input structure before playback starts.
fn synthetic_wobble_frames(
    base: &ViewMolecule,
    n_frames: usize,
    amplitude: f32,
) -> Vec<Vec<[f32; 3]>> {
    let n_frames = n_frames.max(1);
    (0..n_frames)
        .map(|f| {
            let t = 2.0 * std::f32::consts::PI * f as f32 / n_frames as f32;
            let env = amplitude * t.sin(); // 0 at frame 0; same scalar for all atoms
            base.atoms
                .iter()
                .enumerate()
                .map(|(i, atom)| {
                    // Per-atom unit-ish direction, fixed across frames so each
                    // atom oscillates along its own axis (decorrelated breathing).
                    let a = i as f32 * 0.7;
                    let dir = [a.sin(), (a + 2.094).sin(), (a + 4.189).sin()];
                    [
                        atom.pos[0] + env * dir[0],
                        atom.pos[1] + env * dir[1],
                        atom.pos[2] + env * dir[2],
                    ]
                })
                .collect()
        })
        .collect()
}

// ===========================================================================
// In-house Elastic Network Model (ENM) molecular dynamics
// ===========================================================================
//
// A real, physically-motivated MD trajectory of the loaded protein — the
// thermal "breathing" you see in a CASP/PyMOL normal-mode movie — built
// in-house with no external engine and no force field. It is the classic
// **Tirion / anisotropic elastic-network model**: every pair of atoms closer
// than a cutoff in the *reference* structure is joined by a harmonic spring
// at its equilibrium length, the atoms are given Maxwell-Boltzmann thermal
// velocities at a chosen temperature, and the system is propagated with the
// symplectic **velocity-Verlet** integrator under light Langevin damping.
//
// Why ENM (and not the full valenx-md / OPLS force field): the ENM needs only
// the heavy-atom coordinates already in the [`ViewMolecule`] — no bond
// perception, no atom-typing, no partial charges — yet it reproduces the
// large-scale collective motion (the low-frequency normal modes) that
// dominate a protein's thermal fluctuation. That makes it the simplest path
// that is *robust*: it cannot blow up (the potential is purely harmonic about
// the input minimum, so the energy is bounded), it is fully deterministic for
// a fixed seed, and every atom stays within a fraction of an ångström of its
// start — exactly the gentle vibration the viewer should animate.
//
// All maths is done in f64 in a self-consistent reduced unit system
// (length = Å, the spring constant `k` carries the energy/time scale); the
// resulting per-frame Å coordinates feed straight into [`TrajectoryPlayback`].

/// One full ENM-MD result: the per-frame coordinates plus the run diagnostics
/// the panel / agent readout report. Pure data, no rendering.
#[derive(Clone, Debug, Default)]
pub(crate) struct EnmMd {
    /// Per-frame atom coordinates in **ångström**, in `ViewMolecule::atoms`
    /// order — fed directly to [`TrajectoryPlayback::attach`].
    pub(crate) frames: Vec<Vec<[f32; 3]>>,
    /// The requested temperature (reduced units; higher ⇒ larger wobble).
    pub(crate) temperature: f64,
    /// Number of harmonic springs in the network (pairs within the cutoff).
    pub(crate) n_springs: usize,
    /// Total ENM energy (kinetic + spring potential) at the first and last
    /// frame — equal-ish (bounded) is the stability signal.
    pub(crate) energy_first: f64,
    pub(crate) energy_last: f64,
    /// Peak all-atom RMSD (Å) of any frame from the reference (frame 0). A
    /// physical thermal run keeps this well under ~1 Å.
    pub(crate) rmsd_max: f64,
}

impl EnmMd {
    /// A compact one-line readout for the panel / agent bridge.
    pub(crate) fn summary(&self) -> String {
        format!(
            "ENM MD: {} frames · {} springs · T={:.2} · RMSD≤{:.3} Å · E {:.3}→{:.3} (Δ {:+.1}%)",
            self.frames.len(),
            self.n_springs,
            self.temperature,
            self.rmsd_max,
            self.energy_first,
            self.energy_last,
            if self.energy_first.abs() > 1e-12 {
                100.0 * (self.energy_last - self.energy_first) / self.energy_first
            } else {
                0.0
            },
        )
    }
}

/// A tiny deterministic PRNG (SplitMix64) so the thermal velocities are
/// reproducible for a fixed seed without pulling the `rand` crate — the same
/// dependency-light stance as `valenx-md`'s in-crate PRNG.
struct SplitMix64(u64);
impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform f64 in `[0, 1)`.
    fn next_f64(&mut self) -> f64 {
        // Top 53 bits → exact double in [0,1).
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// One standard-normal sample (Box-Muller; the paired value is discarded
    /// for simplicity — we never need bit-exact rand parity, only determinism).
    fn next_gaussian(&mut self) -> f64 {
        // Guard u1 away from 0 so ln() is finite.
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// Atomic mass (u) for the elements crambin contains, with a carbon-weight
/// fallback so an exotic element never yields a zero/NaN acceleration.
fn element_mass(symbol: &str) -> f64 {
    match symbol.trim().to_ascii_uppercase().as_str() {
        "H" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "S" => 32.06,
        "P" => 30.974,
        "FE" => 55.845,
        _ => 12.011,
    }
}

/// Parameters for an ENM-MD run. Defaults are tuned so crambin (327 atoms)
/// produces a visibly vibrating but rock-stable trajectory.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EnmParams {
    /// Spring cutoff in ångström — pairs closer than this in the reference
    /// structure are connected (8 Å is the standard ANM value).
    pub(crate) cutoff: f64,
    /// Harmonic spring constant (reduced energy / Å²). Sets the stiffness /
    /// vibration frequency.
    pub(crate) k: f64,
    /// Target temperature (reduced units). Scales the thermal kinetic energy
    /// and hence the vibration amplitude.
    pub(crate) temperature: f64,
    /// Integration time step (reduced). Kept well inside the stability limit
    /// `dt < 2/ω_max` for the stiffest mode.
    pub(crate) dt: f64,
    /// Langevin velocity-damping coefficient per step — a light drag that
    /// removes the slow numerical energy creep so a long run stays bounded
    /// without killing the visible motion.
    pub(crate) gamma: f64,
    /// Number of trajectory frames to emit.
    pub(crate) n_frames: usize,
    /// Verlet sub-steps integrated between emitted frames (so the motion
    /// advances meaningfully per displayed frame without storing every step).
    pub(crate) substeps: usize,
    /// PRNG seed — fixes the thermal velocities, making the whole run
    /// deterministic.
    pub(crate) seed: u64,
}

impl Default for EnmParams {
    fn default() -> Self {
        EnmParams {
            cutoff: 8.0,
            k: 1.0,
            temperature: 0.30,
            dt: 0.05,
            // Very light per-substep drag: velocity-Verlet on a harmonic
            // potential is already energy-conserving (symplectic, no drift), so
            // this is only insurance against slow numerical creep — small enough
            // that the protein keeps visibly vibrating across the whole reel
            // rather than freezing after a few frames.
            gamma: 0.0005,
            n_frames: 60,
            substeps: 6,
            seed: 0x5EED,
        }
    }
}

/// Run an in-house **ENM velocity-Verlet** MD on a molecule's atoms and return
/// the trajectory + diagnostics. Pure and deterministic for a fixed seed.
///
/// Steps:
/// 1. Build the harmonic spring list: every atom pair within `cutoff` of each
///    other in the reference coordinates, storing the equilibrium length `r0`.
/// 2. Assign Maxwell-Boltzmann thermal velocities `v ~ N(0, kT/m)` from the
///    seeded PRNG, then **remove net linear momentum** (no rigid drift) and
///    rescale so the instantaneous temperature equals `temperature` exactly.
/// 3. Propagate with velocity-Verlet (the symplectic integrator) under light
///    Langevin damping, emitting `n_frames` snapshots `substeps` apart.
///
/// Returns an empty result (no frames) for a molecule with < 2 atoms — nothing
/// to vibrate. Never panics: masses fall back to carbon, an atom with no
/// spring still integrates (it just drifts thermally and is damped back).
//
// The fixed-size-3 coordinate loops (`for c in 0..3`) below index several
// parallel per-atom buffers (`pos`/`vel`/`acc`/`f`) in lockstep; an iterator
// rewrite across 2-3 arrays is less clear than the explicit component index,
// so the range-loop lint is allowed for this routine.
#[allow(clippy::needless_range_loop)]
fn enm_md(mol: &ViewMolecule, p: EnmParams) -> EnmMd {
    let n = mol.atoms.len();
    if n < 2 {
        return EnmMd {
            temperature: p.temperature,
            ..Default::default()
        };
    }

    // --- positions (Å, f64), reference, masses --------------------------
    let r0_ref: Vec<[f64; 3]> = mol
        .atoms
        .iter()
        .map(|a| [a.pos[0] as f64, a.pos[1] as f64, a.pos[2] as f64])
        .collect();
    let mass: Vec<f64> = mol.atoms.iter().map(|a| element_mass(&a.element)).collect();

    // --- spring network: pairs within the cutoff (store r0) -------------
    let cutoff2 = p.cutoff * p.cutoff;
    let mut springs: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = [
                r0_ref[i][0] - r0_ref[j][0],
                r0_ref[i][1] - r0_ref[j][1],
                r0_ref[i][2] - r0_ref[j][2],
            ];
            let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if d2 <= cutoff2 && d2 > 1e-12 {
                springs.push((i, j, d2.sqrt()));
            }
        }
    }

    // --- harmonic spring potential energy + forces ----------------------
    // U = Σ ½k(|rij| - r0)² ;  F on i = -k(|rij| - r0) · (rij / |rij|).
    let forces = |pos: &[[f64; 3]]| -> (Vec<[f64; 3]>, f64) {
        let mut f = vec![[0.0; 3]; n];
        let mut energy = 0.0;
        for &(i, j, r0) in &springs {
            let d = [
                pos[i][0] - pos[j][0],
                pos[i][1] - pos[j][1],
                pos[i][2] - pos[j][2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if r < 1e-9 {
                continue;
            }
            let dr = r - r0;
            energy += 0.5 * p.k * dr * dr;
            // Force magnitude / r, projected onto the bond unit vector.
            let g = -p.k * dr / r;
            for c in 0..3 {
                let fc = g * d[c];
                f[i][c] += fc;
                f[j][c] -= fc;
            }
        }
        (f, energy)
    };

    // --- Maxwell-Boltzmann velocities, momentum-removed, T-rescaled -----
    let mut rng = SplitMix64(p.seed);
    let mut vel = vec![[0.0; 3]; n];
    for (i, v) in vel.iter_mut().enumerate() {
        let sigma = (p.temperature / mass[i]).sqrt(); // √(kT/m), k_B ≡ 1
        for c in 0..3 {
            v[c] = sigma * rng.next_gaussian();
        }
    }
    // Remove net linear momentum so the whole molecule does not translate.
    let mtot: f64 = mass.iter().sum();
    let mut p_net = [0.0; 3];
    for i in 0..n {
        for c in 0..3 {
            p_net[c] += mass[i] * vel[i][c];
        }
    }
    for i in 0..n {
        for c in 0..3 {
            vel[i][c] -= p_net[c] / mtot;
        }
    }
    // Rescale to the requested temperature: T_inst = 2·KE / (3N·k_B).
    let kinetic = |v: &[[f64; 3]]| -> f64 {
        0.5 * (0..n)
            .map(|i| mass[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]))
            .sum::<f64>()
    };
    let dof = 3.0 * n as f64 - 3.0; // momentum removed ⇒ 3 fewer DOF
    let t_inst = 2.0 * kinetic(&vel) / dof; // k_B ≡ 1
    if t_inst > 1e-12 {
        let scale = (p.temperature / t_inst).sqrt();
        for v in &mut vel {
            for c in 0..3 {
                v[c] *= scale;
            }
        }
    }

    // --- velocity-Verlet propagation with light Langevin damping --------
    let mut pos = r0_ref.clone();
    let (mut acc, e0_pot) = forces(&pos);
    for (i, a) in acc.iter_mut().enumerate() {
        for c in 0..3 {
            a[c] /= mass[i];
        }
    }
    let energy_first = kinetic(&vel) + e0_pot;

    let mut frames: Vec<Vec<[f32; 3]>> = Vec::with_capacity(p.n_frames.max(1));
    let snapshot = |pos: &[[f64; 3]]| -> Vec<[f32; 3]> {
        pos.iter()
            .map(|r| [r[0] as f32, r[1] as f32, r[2] as f32])
            .collect()
    };
    frames.push(snapshot(&pos)); // frame 0 = the reference structure exactly

    // Per-step velocity-damping factor (explicit, unconditionally stable for
    // 0 ≤ gamma·dt < 1): v ← v·(1 - gamma) after the kick.
    let damp = (1.0 - p.gamma).clamp(0.0, 1.0);
    let mut energy_last = energy_first;
    for _ in 1..p.n_frames.max(1) {
        for _ in 0..p.substeps.max(1) {
            // v(t+½dt) = v + ½dt·a ; r += dt·v ; a' = F(r)/m ; v += ½dt·a'.
            for i in 0..n {
                for c in 0..3 {
                    vel[i][c] += 0.5 * p.dt * acc[i][c];
                    pos[i][c] += p.dt * vel[i][c];
                }
            }
            let (f, _e) = forces(&pos);
            for i in 0..n {
                for c in 0..3 {
                    acc[i][c] = f[i][c] / mass[i];
                    vel[i][c] += 0.5 * p.dt * acc[i][c];
                    vel[i][c] *= damp; // light Langevin drag
                }
            }
        }
        let (_f, pe) = forces(&pos);
        energy_last = kinetic(&vel) + pe;
        frames.push(snapshot(&pos));
    }

    // --- peak RMSD from the reference (frame 0) -------------------------
    let mut rmsd_max = 0.0_f64;
    for frame in &frames {
        let mut s = 0.0;
        for (a, r) in frame.iter().zip(&r0_ref) {
            let dx = a[0] as f64 - r[0];
            let dy = a[1] as f64 - r[1];
            let dz = a[2] as f64 - r[2];
            s += dx * dx + dy * dy + dz * dz;
        }
        rmsd_max = rmsd_max.max((s / n as f64).sqrt());
    }

    EnmMd {
        frames,
        temperature: p.temperature,
        n_springs: springs.len(),
        energy_first,
        energy_last,
        rmsd_max,
    }
}

/// Which structure sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Analyze,
    Ramachandran,
    Superpose,
}

/// Snapshot of every editable input the Biostruct panel owns.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) struct BiostructSnapshot {
    pub(crate) structure_a: String,
    pub(crate) structure_b: String,
    pub(crate) clash_tolerance: f64,
}

/// Form + result state for the Macromolecular Structure panel.
pub struct BiostructPanel {
    tool: Tool,
    /// Structure text for the analyze / Ramachandran tools, and the
    /// mobile structure of the superpose tool.
    structure_a: String,
    /// Reference structure text for the superpose tool.
    structure_b: String,
    /// Steric-clash tolerance (Å) for the analysis report.
    clash_tolerance: f64,
    error: Option<String>,
    result: String,
    /// Undo / redo over both structure-text inputs + the clash tol.
    history: crate::undo::History<BiostructSnapshot>,
    /// Which molecular-viewer representation the "Show in 3D viewport" button
    /// builds. Reactive picker; default [`Representation::BallAndStick`].
    pub(crate) representation: Representation,
    /// How the viewer colours atoms. Reactive picker; default
    /// [`ColorScheme::Element`] (the CPK palette). When non-default
    /// (chain / residue / B-factor), "Show in 3D viewport" builds a
    /// per-vertex-coloured mesh and uploads the colours so the viewport
    /// renders the scheme instead of monochrome.
    pub(crate) color_scheme: ColorScheme,
    /// Per-representation mesh-generation tunables (surface grid resolution,
    /// cartoon tube, ball/stick radii). Mutated by the picker's sliders.
    pub(crate) molviz_params: MolvizParams,
    /// MD-trajectory playback over the currently-shown structure — animate the
    /// atoms across coordinate frames (synthetic wobble or a loaded `valenx-md`
    /// trajectory). Empty until a trajectory is attached.
    pub(crate) trajectory: TrajectoryPlayback,
    /// Target temperature (reduced units) for the in-house **Run MD (ENM)**
    /// control — scales the thermal vibration amplitude. Reactive slider;
    /// also settable through the agent `agent_set` bridge as `"MD temperature"`.
    pub(crate) md_temperature: f64,
    /// Which built-in demo structure the "Demo structure" picker last loaded
    /// into `structure_a`. Default [`DemoStructure::Crambin`] (the real
    /// 46-residue protein, also `structure_a`'s default text).
    pub(crate) demo: DemoStructure,
}

impl BiostructPanel {
    fn snapshot(&self) -> BiostructSnapshot {
        BiostructSnapshot {
            structure_a: self.structure_a.clone(),
            structure_b: self.structure_b.clone(),
            clash_tolerance: self.clash_tolerance,
        }
    }
    fn restore(&mut self, s: BiostructSnapshot) {
        self.structure_a = s.structure_a;
        self.structure_b = s.structure_b;
        self.clash_tolerance = s.clash_tolerance;
    }
    pub fn undo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(prev) = self.history.undo(current) {
            self.restore(prev);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn redo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(next) = self.history.redo(current) {
            self.restore(next);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }

    /// Load a built-in [`DemoStructure`] into `structure_a` (the viewer / analysis
    /// input), recording the prior state for undo. Used by both the "Demo
    /// structure" picker and the agent `SetControl` bridge, so selecting a demo
    /// through the accessibility tree and selecting it by hand take the same path.
    /// Clears any prior error and detaches a stale trajectory (the new structure
    /// has a different atom count).
    pub(crate) fn load_demo(&mut self, demo: DemoStructure) {
        let snap = self.snapshot();
        self.history.record(snap);
        self.demo = demo;
        self.structure_a = demo.pdb().to_string();
        self.error = None;
        self.trajectory.clear();
    }

    /// Run the in-house **ENM velocity-Verlet MD** on the structure currently in
    /// `structure_a` and attach the resulting real thermal-vibration trajectory
    /// to the viewer's playback (so pressing Play animates the protein
    /// breathing). Uses the default ENM parameters with the panel's
    /// `md_temperature`. Shared by the "Run MD (ENM)" button and the agent
    /// `RunMd` bridge so both take the same path.
    ///
    /// Fail-loud: an unparseable structure sets the trajectory note and attaches
    /// nothing (never panics). Returns the run summary on success for the agent.
    pub(crate) fn run_md(&mut self) -> Result<String, String> {
        match read_structure(&self.structure_a, "viewer") {
            Ok(s) => {
                let params = EnmParams {
                    temperature: self.md_temperature,
                    ..EnmParams::default()
                };
                self.trajectory.attach_enm_md(&s, params);
                // Sensible default playback speed for the ~60-frame run.
                self.trajectory.speed = 12.0;
                match self.trajectory.md() {
                    Some(md) => Ok(md.summary()),
                    None => Err(self
                        .trajectory
                        .note
                        .clone()
                        .unwrap_or_else(|| "MD did not attach".to_string())),
                }
            }
            Err(e) => {
                let msg = format!("could not parse structure: {e}");
                self.trajectory.note = Some(msg.clone());
                Err(msg)
            }
        }
    }

    /// A one-line MD readout for the agent bridge: the last ENM-MD run summary
    /// (frames / springs / temperature / RMSD / energy), or a hint that no run
    /// has happened yet.
    pub(crate) fn md_readout(&self) -> String {
        match self.trajectory.md() {
            Some(md) => md.summary(),
            None => format!(
                "no MD run yet (Run MD (ENM) to animate; target T={:.2})",
                self.md_temperature
            ),
        }
    }

    /// A one-line readout of the structure currently loaded in `structure_a`,
    /// for the agent `agent_readout` bridge: the demo name (if it matches a
    /// built-in) plus the parsed residue / atom counts. Fail-soft — an
    /// unparseable buffer reports the parse error rather than panicking.
    pub(crate) fn structure_readout(&self) -> String {
        let src = if self.structure_a == CRAMBIN_PDB {
            "Crambin (1CRN)"
        } else if self.structure_a == DEMO_PDB {
            "Triglycine peptide"
        } else {
            "custom / loaded structure"
        };
        match read_structure(&self.structure_a, "readout") {
            Ok(s) => {
                let (mut res, mut atoms) = (0usize, 0usize);
                for chain in &s.first_model().chains {
                    res += chain.residues.len();
                    for r in &chain.residues {
                        atoms += r.atoms.len();
                    }
                }
                format!("{src}: {res} residues, {atoms} atoms")
            }
            Err(e) => format!("{src}: unparseable ({e})"),
        }
    }
}

/// A minimal 3-residue glycine peptide PDB — enough for the analysis +
/// Ramachandran tools to produce real output without file I/O. Tiny by design
/// (12 atoms): a smoke structure, not something that looks like a protein. The
/// default loaded structure is [`CRAMBIN_PDB`]; this stays as a second demo
/// option so the trivial case is still one click away.
const DEMO_PDB: &str = "\
ATOM      1  N   GLY A   1      -1.204   1.045   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       1.250   0.881   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       1.150   2.100   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       2.430   0.300   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       3.720   0.960   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       4.880   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       4.770  -1.220   0.000  1.00  0.00           O
ATOM      9  N   GLY A   3       6.050   0.620   0.000  1.00  0.00           N
ATOM     10  CA  GLY A   3       7.310  -0.080   0.000  1.00  0.00           C
ATOM     11  C   GLY A   3       8.500   0.870   0.000  1.00  0.00           C
ATOM     12  O   GLY A   3       8.380   2.090   0.000  1.00  0.00           O
END
";

/// **Crambin (PDB 1CRN)** — a real 46-residue, 327-atom plant-seed protein, the
/// default loaded structure for the molecular viewer. It is a classic small
/// fold that carries *both* an α-helix (two helices, res 7–19 and 23–30) and a
/// β-sheet (two strands, res 1–4 and 32–35), so the ribbon / cartoon renderers
/// and the secondary-structure (DSSP) colour scheme all have something real to
/// show — unlike the 12-atom [`DEMO_PDB`], which looks like nothing.
///
/// The atom/`HELIX`/`SHEET`/`TER` records are the original experimental
/// coordinates from Teeter & Hendrickson (1981); the existing
/// [`read_structure`] PDB parser consumes them as-is (it already handles
/// `HELIX`, `SHEET`, `TER`, and the `OXT` terminal-oxygen atom). PDB coordinate
/// data is public-domain scientific data.
const CRAMBIN_PDB: &str = "\
HEADER    PLANT PROTEIN                           30-APR-81   1CRN
HELIX    1  H1 ILE A    7  PRO A   19  13/10 CONFORMATION RES 17,19       13
HELIX    2  H2 GLU A   23  THR A   30  1DISTORTED 3/10 AT RES 30           8
SHEET    1  S1 2 THR A   1  CYS A   4  0
SHEET    2  S1 2 CYS A  32  ILE A  35 -1
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O
ATOM      7  CG2 THR A   1      18.150  11.546   6.304  1.00 14.23           C
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  8.31           C
ATOM     10  C   THR A   2      14.164  10.785   7.379  1.00  5.80           C
ATOM     11  O   THR A   2      14.993   9.862   7.443  1.00  6.94           O
ATOM     12  CB  THR A   2      12.732  10.711   5.261  1.00 10.32           C
ATOM     13  OG1 THR A   2      13.308   9.439   4.926  1.00 12.81           O
ATOM     14  CG2 THR A   2      12.484  11.442   3.895  1.00 11.90           C
ATOM     15  N   CYS A   3      13.488  11.241   8.417  1.00  5.24           N
ATOM     16  CA  CYS A   3      13.660  10.707   9.787  1.00  5.39           C
ATOM     17  C   CYS A   3      12.269  10.431  10.323  1.00  4.45           C
ATOM     18  O   CYS A   3      11.393  11.308  10.185  1.00  6.54           O
ATOM     19  CB  CYS A   3      14.368  11.748  10.691  1.00  5.99           C
ATOM     20  SG  CYS A   3      15.885  12.426  10.016  1.00  7.01           S
ATOM     21  N   CYS A   4      12.019   9.272  10.928  1.00  3.90           N
ATOM     22  CA  CYS A   4      10.646   8.991  11.408  1.00  4.24           C
ATOM     23  C   CYS A   4      10.654   8.793  12.919  1.00  3.72           C
ATOM     24  O   CYS A   4      11.659   8.296  13.491  1.00  5.30           O
ATOM     25  CB  CYS A   4      10.057   7.752  10.682  1.00  4.41           C
ATOM     26  SG  CYS A   4       9.837   8.018   8.904  1.00  4.72           S
ATOM     27  N   PRO A   5       9.561   9.108  13.563  1.00  3.96           N
ATOM     28  CA  PRO A   5       9.448   9.034  15.012  1.00  4.25           C
ATOM     29  C   PRO A   5       9.288   7.670  15.606  1.00  4.96           C
ATOM     30  O   PRO A   5       9.490   7.519  16.819  1.00  7.44           O
ATOM     31  CB  PRO A   5       8.230   9.957  15.345  1.00  5.11           C
ATOM     32  CG  PRO A   5       7.338   9.786  14.114  1.00  5.24           C
ATOM     33  CD  PRO A   5       8.366   9.804  12.958  1.00  5.20           C
ATOM     34  N   SER A   6       8.875   6.686  14.796  1.00  4.83           N
ATOM     35  CA  SER A   6       8.673   5.314  15.279  1.00  4.45           C
ATOM     36  C   SER A   6       8.753   4.376  14.083  1.00  4.99           C
ATOM     37  O   SER A   6       8.726   4.858  12.923  1.00  4.61           O
ATOM     38  CB  SER A   6       7.340   5.121  15.996  1.00  5.05           C
ATOM     39  OG  SER A   6       6.274   5.220  15.031  1.00  6.39           O
ATOM     40  N   ILE A   7       8.881   3.075  14.358  1.00  4.94           N
ATOM     41  CA  ILE A   7       8.912   2.083  13.258  1.00  6.33           C
ATOM     42  C   ILE A   7       7.581   2.090  12.506  1.00  5.32           C
ATOM     43  O   ILE A   7       7.670   2.031  11.245  1.00  6.85           O
ATOM     44  CB  ILE A   7       9.207   0.677  13.924  1.00  8.43           C
ATOM     45  CG1 ILE A   7      10.714   0.702  14.312  1.00  9.78           C
ATOM     46  CG2 ILE A   7       8.811  -0.477  12.969  1.00 11.70           C
ATOM     47  CD1 ILE A   7      11.185  -0.516  15.142  1.00  9.92           C
ATOM     48  N   VAL A   8       6.458   2.162  13.159  1.00  5.02           N
ATOM     49  CA  VAL A   8       5.145   2.209  12.453  1.00  6.93           C
ATOM     50  C   VAL A   8       5.115   3.379  11.461  1.00  5.39           C
ATOM     51  O   VAL A   8       4.664   3.268  10.343  1.00  6.30           O
ATOM     52  CB  VAL A   8       3.995   2.354  13.478  1.00  9.64           C
ATOM     53  CG1 VAL A   8       2.716   2.891  12.869  1.00 13.85           C
ATOM     54  CG2 VAL A   8       3.758   1.032  14.208  1.00 11.97           C
ATOM     55  N   ALA A   9       5.606   4.546  11.941  1.00  3.73           N
ATOM     56  CA  ALA A   9       5.598   5.767  11.082  1.00  3.56           C
ATOM     57  C   ALA A   9       6.441   5.527   9.850  1.00  4.13           C
ATOM     58  O   ALA A   9       6.052   5.933   8.744  1.00  4.36           O
ATOM     59  CB  ALA A   9       6.022   6.977  11.891  1.00  4.80           C
ATOM     60  N   ARG A  10       7.647   4.909  10.005  1.00  3.73           N
ATOM     61  CA  ARG A  10       8.496   4.609   8.837  1.00  3.38           C
ATOM     62  C   ARG A  10       7.798   3.609   7.876  1.00  3.47           C
ATOM     63  O   ARG A  10       7.878   3.778   6.651  1.00  4.67           O
ATOM     64  CB  ARG A  10       9.847   4.020   9.305  1.00  3.95           C
ATOM     65  CG  ARG A  10      10.752   3.607   8.149  1.00  4.55           C
ATOM     66  CD  ARG A  10      11.226   4.699   7.244  1.00  5.89           C
ATOM     67  NE  ARG A  10      12.143   5.571   8.035  1.00  6.20           N
ATOM     68  CZ  ARG A  10      12.758   6.609   7.443  1.00  7.52           C
ATOM     69  NH1 ARG A  10      12.539   6.932   6.158  1.00 10.68           N
ATOM     70  NH2 ARG A  10      13.601   7.322   8.202  1.00  9.48           N
ATOM     71  N   SER A  11       7.186   2.582   8.445  1.00  5.19           N
ATOM     72  CA  SER A  11       6.500   1.584   7.565  1.00  4.60           C
ATOM     73  C   SER A  11       5.382   2.313   6.773  1.00  4.84           C
ATOM     74  O   SER A  11       5.213   2.016   5.557  1.00  5.84           O
ATOM     75  CB  SER A  11       5.908   0.462   8.400  1.00  5.91           C
ATOM     76  OG  SER A  11       6.990  -0.272   9.012  1.00  8.38           O
ATOM     77  N   ASN A  12       4.648   3.182   7.446  1.00  3.54           N
ATOM     78  CA  ASN A  12       3.545   3.935   6.751  1.00  4.57           C
ATOM     79  C   ASN A  12       4.107   4.851   5.691  1.00  4.14           C
ATOM     80  O   ASN A  12       3.536   5.001   4.617  1.00  5.52           O
ATOM     81  CB  ASN A  12       2.663   4.677   7.748  1.00  6.42           C
ATOM     82  CG  ASN A  12       1.802   3.735   8.610  1.00  8.25           C
ATOM     83  OD1 ASN A  12       1.567   2.613   8.165  1.00 12.72           O
ATOM     84  ND2 ASN A  12       1.394   4.252   9.767  1.00  9.92           N
ATOM     85  N   PHE A  13       5.259   5.498   6.005  1.00  3.43           N
ATOM     86  CA  PHE A  13       5.929   6.358   5.055  1.00  3.49           C
ATOM     87  C   PHE A  13       6.304   5.578   3.799  1.00  3.40           C
ATOM     88  O   PHE A  13       6.136   6.072   2.653  1.00  4.07           O
ATOM     89  CB  PHE A  13       7.183   6.994   5.754  1.00  5.48           C
ATOM     90  CG  PHE A  13       7.884   8.006   4.883  1.00  5.57           C
ATOM     91  CD1 PHE A  13       8.906   7.586   4.027  1.00  6.99           C
ATOM     92  CD2 PHE A  13       7.532   9.373   4.983  1.00  6.52           C
ATOM     93  CE1 PHE A  13       9.560   8.539   3.194  1.00  8.20           C
ATOM     94  CE2 PHE A  13       8.176  10.281   4.145  1.00  6.34           C
ATOM     95  CZ  PHE A  13       9.141   9.845   3.292  1.00  6.84           C
ATOM     96  N   ASN A  14       6.900   4.390   3.989  1.00  3.64           N
ATOM     97  CA  ASN A  14       7.331   3.607   2.791  1.00  4.31           C
ATOM     98  C   ASN A  14       6.116   3.210   1.915  1.00  3.98           C
ATOM     99  O   ASN A  14       6.240   3.144   0.684  1.00  6.22           O
ATOM    100  CB  ASN A  14       8.145   2.404   3.240  1.00  5.81           C
ATOM    101  CG  ASN A  14       9.555   2.856   3.730  1.00  6.82           C
ATOM    102  OD1 ASN A  14      10.013   3.895   3.323  1.00  9.43           O
ATOM    103  ND2 ASN A  14      10.120   1.956   4.539  1.00  8.21           N
ATOM    104  N   VAL A  15       4.993   2.927   2.571  1.00  3.76           N
ATOM    105  CA  VAL A  15       3.782   2.599   1.742  1.00  3.98           C
ATOM    106  C   VAL A  15       3.296   3.871   1.004  1.00  3.80           C
ATOM    107  O   VAL A  15       2.947   3.817  -0.189  1.00  4.85           O
ATOM    108  CB  VAL A  15       2.698   1.953   2.608  1.00  4.71           C
ATOM    109  CG1 VAL A  15       1.384   1.826   1.806  1.00  6.67           C
ATOM    110  CG2 VAL A  15       3.174   0.533   3.005  1.00  6.26           C
ATOM    111  N   CYS A  16       3.321   4.987   1.720  1.00  3.79           N
ATOM    112  CA  CYS A  16       2.890   6.285   1.126  1.00  3.54           C
ATOM    113  C   CYS A  16       3.687   6.597  -0.111  1.00  3.48           C
ATOM    114  O   CYS A  16       3.200   7.147  -1.103  1.00  4.63           O
ATOM    115  CB  CYS A  16       3.039   7.369   2.240  1.00  4.58           C
ATOM    116  SG  CYS A  16       2.559   9.014   1.649  1.00  5.66           S
ATOM    117  N   ARG A  17       4.997   6.227  -0.100  1.00  3.99           N
ATOM    118  CA  ARG A  17       5.895   6.489  -1.213  1.00  3.83           C
ATOM    119  C   ARG A  17       5.738   5.560  -2.409  1.00  3.79           C
ATOM    120  O   ARG A  17       6.228   5.901  -3.507  1.00  5.39           O
ATOM    121  CB  ARG A  17       7.370   6.507  -0.731  1.00  4.11           C
ATOM    122  CG  ARG A  17       7.717   7.687   0.206  1.00  4.69           C
ATOM    123  CD  ARG A  17       7.949   8.947  -0.615  1.00  5.10           C
ATOM    124  NE  ARG A  17       9.212   8.856  -1.337  1.00  4.71           N
ATOM    125  CZ  ARG A  17       9.537   9.533  -2.431  1.00  5.28           C
ATOM    126  NH1 ARG A  17       8.659  10.350  -3.032  1.00  6.67           N
ATOM    127  NH2 ARG A  17      10.793   9.491  -2.899  1.00  6.41           N
ATOM    128  N   LEU A  18       5.051   4.411  -2.204  1.00  4.70           N
ATOM    129  CA  LEU A  18       4.933   3.431  -3.326  1.00  5.46           C
ATOM    130  C   LEU A  18       4.397   4.014  -4.620  1.00  5.13           C
ATOM    131  O   LEU A  18       4.988   3.755  -5.687  1.00  5.55           O
ATOM    132  CB  LEU A  18       4.196   2.184  -2.863  1.00  6.47           C
ATOM    133  CG  LEU A  18       4.960   1.178  -1.991  1.00  7.43           C
ATOM    134  CD1 LEU A  18       3.907   0.097  -1.634  1.00  8.70           C
ATOM    135  CD2 LEU A  18       6.129   0.606  -2.768  1.00  9.39           C
ATOM    136  N   PRO A  19       3.329   4.795  -4.543  1.00  4.28           N
ATOM    137  CA  PRO A  19       2.792   5.376  -5.797  1.00  5.38           C
ATOM    138  C   PRO A  19       3.573   6.540  -6.322  1.00  6.30           C
ATOM    139  O   PRO A  19       3.260   7.045  -7.422  1.00  9.62           O
ATOM    140  CB  PRO A  19       1.358   5.766  -5.472  1.00  5.87           C
ATOM    141  CG  PRO A  19       1.223   5.694  -3.993  1.00  6.47           C
ATOM    142  CD  PRO A  19       2.421   4.941  -3.408  1.00  6.45           C
ATOM    143  N   GLY A  20       4.565   7.047  -5.559  1.00  4.94           N
ATOM    144  CA  GLY A  20       5.366   8.191  -6.018  1.00  5.39           C
ATOM    145  C   GLY A  20       5.007   9.481  -5.280  1.00  5.03           C
ATOM    146  O   GLY A  20       5.535  10.510  -5.730  1.00  7.34           O
ATOM    147  N   THR A  21       4.181   9.438  -4.262  1.00  4.10           N
ATOM    148  CA  THR A  21       3.767  10.609  -3.513  1.00  3.94           C
ATOM    149  C   THR A  21       5.017  11.397  -3.042  1.00  3.96           C
ATOM    150  O   THR A  21       5.947  10.757  -2.523  1.00  5.82           O
ATOM    151  CB  THR A  21       2.992  10.188  -2.225  1.00  4.13           C
ATOM    152  OG1 THR A  21       2.051   9.144  -2.623  1.00  5.45           O
ATOM    153  CG2 THR A  21       2.260  11.349  -1.551  1.00  5.41           C
ATOM    154  N   PRO A  22       4.971  12.703  -3.176  1.00  5.04           N
ATOM    155  CA  PRO A  22       6.143  13.513  -2.696  1.00  4.69           C
ATOM    156  C   PRO A  22       6.400  13.233  -1.225  1.00  4.19           C
ATOM    157  O   PRO A  22       5.485  13.061  -0.382  1.00  4.47           O
ATOM    158  CB  PRO A  22       5.703  14.969  -2.920  1.00  7.12           C
ATOM    159  CG  PRO A  22       4.676  14.893  -3.996  1.00  7.03           C
ATOM    160  CD  PRO A  22       3.964  13.567  -3.811  1.00  4.90           C
ATOM    161  N   GLU A  23       7.728  13.297  -0.921  1.00  5.16           N
ATOM    162  CA  GLU A  23       8.114  13.103   0.500  1.00  5.31           C
ATOM    163  C   GLU A  23       7.427  14.073   1.410  1.00  4.11           C
ATOM    164  O   GLU A  23       7.036  13.682   2.540  1.00  5.11           O
ATOM    165  CB  GLU A  23       9.648  13.285   0.660  1.00  6.16           C
ATOM    166  CG  GLU A  23      10.440  12.093   0.063  1.00  7.48           C
ATOM    167  CD  GLU A  23      11.941  12.170   0.391  1.00  9.40           C
ATOM    168  OE1 GLU A  23      12.416  13.225   0.681  1.00 10.40           O
ATOM    169  OE2 GLU A  23      12.539  11.070   0.292  1.00 13.32           O
ATOM    170  N   ALA A  24       7.212  15.334   0.966  1.00  4.56           N
ATOM    171  CA  ALA A  24       6.614  16.317   1.913  1.00  4.49           C
ATOM    172  C   ALA A  24       5.212  15.936   2.350  1.00  4.10           C
ATOM    173  O   ALA A  24       4.782  16.166   3.495  1.00  5.64           O
ATOM    174  CB  ALA A  24       6.605  17.695   1.246  1.00  5.80           C
ATOM    175  N   ILE A  25       4.445  15.318   1.405  1.00  4.37           N
ATOM    176  CA  ILE A  25       3.074  14.894   1.756  1.00  5.44           C
ATOM    177  C   ILE A  25       3.085  13.643   2.645  1.00  4.32           C
ATOM    178  O   ILE A  25       2.315  13.523   3.578  1.00  4.72           O
ATOM    179  CB  ILE A  25       2.204  14.637   0.462  1.00  6.42           C
ATOM    180  CG1 ILE A  25       1.815  16.048  -0.129  1.00  7.50           C
ATOM    181  CG2 ILE A  25       0.903  13.864   0.811  1.00  7.65           C
ATOM    182  CD1 ILE A  25       0.756  16.761   0.757  1.00  7.80           C
ATOM    183  N   CYS A  26       4.032  12.764   2.313  1.00  3.92           N
ATOM    184  CA  CYS A  26       4.180  11.549   3.187  1.00  4.37           C
ATOM    185  C   CYS A  26       4.632  11.944   4.596  1.00  3.95           C
ATOM    186  O   CYS A  26       4.227  11.252   5.547  1.00  4.74           O
ATOM    187  CB  CYS A  26       5.038  10.518   2.539  1.00  4.63           C
ATOM    188  SG  CYS A  26       4.349   9.794   1.022  1.00  5.61           S
ATOM    189  N   ALA A  27       5.408  13.012   4.694  1.00  3.89           N
ATOM    190  CA  ALA A  27       5.879  13.502   6.026  1.00  4.43           C
ATOM    191  C   ALA A  27       4.696  13.908   6.882  1.00  4.26           C
ATOM    192  O   ALA A  27       4.528  13.422   8.025  1.00  5.44           O
ATOM    193  CB  ALA A  27       6.880  14.615   5.830  1.00  5.36           C
ATOM    194  N   THR A  28       3.827  14.802   6.358  1.00  4.53           N
ATOM    195  CA  THR A  28       2.691  15.221   7.194  1.00  5.08           C
ATOM    196  C   THR A  28       1.672  14.132   7.434  1.00  4.62           C
ATOM    197  O   THR A  28       0.947  14.112   8.468  1.00  7.80           O
ATOM    198  CB  THR A  28       1.986  16.520   6.614  1.00  6.03           C
ATOM    199  OG1 THR A  28       1.664  16.221   5.230  1.00  7.19           O
ATOM    200  CG2 THR A  28       2.914  17.739   6.700  1.00  7.34           C
ATOM    201  N   TYR A  29       1.621  13.190   6.511  1.00  5.01           N
ATOM    202  CA  TYR A  29       0.715  12.045   6.657  1.00  6.60           C
ATOM    203  C   TYR A  29       1.125  11.125   7.815  1.00  4.92           C
ATOM    204  O   TYR A  29       0.286  10.632   8.545  1.00  7.13           O
ATOM    205  CB  TYR A  29       0.755  11.229   5.322  1.00  9.66           C
ATOM    206  CG  TYR A  29      -0.203  10.044   5.354  1.00 11.56           C
ATOM    207  CD1 TYR A  29      -1.547  10.337   5.645  1.00 12.85           C
ATOM    208  CD2 TYR A  29       0.193   8.750   5.100  1.00 14.44           C
ATOM    209  CE1 TYR A  29      -2.496   9.329   5.673  1.00 16.61           C
ATOM    210  CE2 TYR A  29      -0.801   7.705   5.156  1.00 17.11           C
ATOM    211  CZ  TYR A  29      -2.079   8.031   5.430  1.00 19.99           C
ATOM    212  OH  TYR A  29      -3.097   7.057   5.458  1.00 28.98           O
ATOM    213  N   THR A  30       2.470  10.984   7.995  1.00  5.31           N
ATOM    214  CA  THR A  30       2.986   9.994   8.950  1.00  5.70           C
ATOM    215  C   THR A  30       3.609  10.505  10.230  1.00  6.28           C
ATOM    216  O   THR A  30       3.766   9.715  11.186  1.00  8.77           O
ATOM    217  CB  THR A  30       4.076   9.103   8.225  1.00  6.55           C
ATOM    218  OG1 THR A  30       5.125  10.027   7.824  1.00  6.57           O
ATOM    219  CG2 THR A  30       3.493   8.324   7.035  1.00  7.29           C
ATOM    220  N   GLY A  31       3.984  11.764  10.241  1.00  4.99           N
ATOM    221  CA  GLY A  31       4.769  12.336  11.360  1.00  5.50           C
ATOM    222  C   GLY A  31       6.255  12.243  11.106  1.00  4.19           C
ATOM    223  O   GLY A  31       7.037  12.750  11.954  1.00  6.12           O
ATOM    224  N   CYS A  32       6.710  11.631   9.992  1.00  4.30           N
ATOM    225  CA  CYS A  32       8.140  11.694   9.635  1.00  4.89           C
ATOM    226  C   CYS A  32       8.500  13.141   9.206  1.00  5.50           C
ATOM    227  O   CYS A  32       7.581  13.949   8.944  1.00  5.82           O
ATOM    228  CB  CYS A  32       8.504  10.686   8.530  1.00  4.66           C
ATOM    229  SG  CYS A  32       8.048   8.987   8.881  1.00  5.33           S
ATOM    230  N   ILE A  33       9.793  13.410   9.173  1.00  6.02           N
ATOM    231  CA  ILE A  33      10.280  14.760   8.823  1.00  5.24           C
ATOM    232  C   ILE A  33      11.346  14.658   7.743  1.00  5.16           C
ATOM    233  O   ILE A  33      11.971  13.583   7.552  1.00  7.19           O
ATOM    234  CB  ILE A  33      10.790  15.535  10.085  1.00  5.49           C
ATOM    235  CG1 ILE A  33      12.059  14.803  10.671  1.00  6.85           C
ATOM    236  CG2 ILE A  33       9.684  15.686  11.138  1.00  6.45           C
ATOM    237  CD1 ILE A  33      12.733  15.676  11.781  1.00  8.94           C
ATOM    238  N   ILE A  34      11.490  15.773   7.038  1.00  5.52           N
ATOM    239  CA  ILE A  34      12.552  15.877   6.036  1.00  6.82           C
ATOM    240  C   ILE A  34      13.590  16.917   6.560  1.00  6.92           C
ATOM    241  O   ILE A  34      13.168  18.006   6.945  1.00  9.22           O
ATOM    242  CB  ILE A  34      11.987  16.360   4.681  1.00  8.11           C
ATOM    243  CG1 ILE A  34      10.914  15.338   4.163  1.00  9.59           C
ATOM    244  CG2 ILE A  34      13.131  16.517   3.629  1.00  9.73           C
ATOM    245  CD1 ILE A  34      10.151  16.024   2.938  1.00 13.41           C
ATOM    246  N   ILE A  35      14.856  16.493   6.536  1.00  7.06           N
ATOM    247  CA  ILE A  35      15.930  17.454   6.941  1.00  7.52           C
ATOM    248  C   ILE A  35      16.913  17.550   5.819  1.00  6.63           C
ATOM    249  O   ILE A  35      17.097  16.660   4.970  1.00  7.90           O
ATOM    250  CB  ILE A  35      16.622  16.995   8.285  1.00  8.07           C
ATOM    251  CG1 ILE A  35      17.360  15.651   8.067  1.00  9.41           C
ATOM    252  CG2 ILE A  35      15.592  16.974   9.434  1.00  9.46           C
ATOM    253  CD1 ILE A  35      18.298  15.206   9.219  1.00  9.85           C
ATOM    254  N   PRO A  36      17.664  18.669   5.806  1.00  8.07           N
ATOM    255  CA  PRO A  36      18.635  18.861   4.738  1.00  8.78           C
ATOM    256  C   PRO A  36      19.925  18.042   4.949  1.00  8.31           C
ATOM    257  O   PRO A  36      20.593  17.742   3.945  1.00  9.09           O
ATOM    258  CB  PRO A  36      18.945  20.364   4.783  1.00  9.67           C
ATOM    259  CG  PRO A  36      18.238  20.937   5.908  1.00 10.15           C
ATOM    260  CD  PRO A  36      17.371  19.900   6.596  1.00  9.53           C
ATOM    261  N   GLY A  37      20.172  17.730   6.217  1.00  8.48           N
ATOM    262  CA  GLY A  37      21.452  16.969   6.513  1.00  9.20           C
ATOM    263  C   GLY A  37      21.143  15.478   6.427  1.00 10.41           C
ATOM    264  O   GLY A  37      20.138  15.023   5.878  1.00 12.06           O
ATOM    265  N   ALA A  38      22.055  14.701   7.032  1.00  9.24           N
ATOM    266  CA  ALA A  38      22.019  13.242   7.020  1.00  9.24           C
ATOM    267  C   ALA A  38      21.944  12.628   8.396  1.00  9.60           C
ATOM    268  O   ALA A  38      21.869  11.387   8.435  1.00 13.65           O
ATOM    269  CB  ALA A  38      23.246  12.697   6.275  1.00 10.43           C
ATOM    270  N   THR A  39      21.894  13.435   9.436  1.00  8.70           N
ATOM    271  CA  THR A  39      21.936  12.911  10.809  1.00  9.46           C
ATOM    272  C   THR A  39      20.615  13.191  11.521  1.00  8.32           C
ATOM    273  O   THR A  39      20.357  14.317  11.948  1.00  9.89           O
ATOM    274  CB  THR A  39      23.131  13.601  11.593  1.00 10.72           C
ATOM    275  OG1 THR A  39      24.284  13.401  10.709  1.00 11.66           O
ATOM    276  CG2 THR A  39      23.340  12.935  12.962  1.00 11.81           C
ATOM    277  N   CYS A  40      19.827  12.110  11.642  1.00  7.64           N
ATOM    278  CA  CYS A  40      18.504  12.312  12.298  1.00  8.05           C
ATOM    279  C   CYS A  40      18.684  12.451  13.784  1.00  7.63           C
ATOM    280  O   CYS A  40      19.533  11.718  14.362  1.00  9.64           O
ATOM    281  CB  CYS A  40      17.582  11.117  11.996  1.00  7.80           C
ATOM    282  SG  CYS A  40      17.199  10.929  10.237  1.00  7.30           S
ATOM    283  N   PRO A  41      17.880  13.266  14.426  1.00  8.00           N
ATOM    284  CA  PRO A  41      17.924  13.421  15.877  1.00  8.96           C
ATOM    285  C   PRO A  41      17.392  12.206  16.594  1.00  9.06           C
ATOM    286  O   PRO A  41      16.652  11.368  16.033  1.00  8.82           O
ATOM    287  CB  PRO A  41      17.076  14.658  16.145  1.00 10.39           C
ATOM    288  CG  PRO A  41      16.098  14.689  14.997  1.00 10.99           C
ATOM    289  CD  PRO A  41      16.859  14.150  13.779  1.00 10.49           C
ATOM    290  N   GLY A  42      17.728  12.124  17.884  1.00  7.55           N
ATOM    291  CA  GLY A  42      17.334  10.956  18.691  1.00  8.00           C
ATOM    292  C   GLY A  42      15.875  10.688  18.871  1.00  7.22           C
ATOM    293  O   GLY A  42      15.434   9.550  19.166  1.00  8.41           O
ATOM    294  N   ASP A  43      15.036  11.747  18.715  1.00  5.54           N
ATOM    295  CA  ASP A  43      13.564  11.573  18.836  1.00  5.85           C
ATOM    296  C   ASP A  43      12.936  11.227  17.470  1.00  5.87           C
ATOM    297  O   ASP A  43      11.720  11.040  17.428  1.00  7.29           O
ATOM    298  CB  ASP A  43      12.933  12.737  19.580  1.00  6.72           C
ATOM    299  CG  ASP A  43      13.140  14.094  18.958  1.00  8.59           C
ATOM    300  OD1 ASP A  43      14.109  14.303  18.212  1.00  9.59           O
ATOM    301  OD2 ASP A  43      12.267  14.963  19.265  1.00 11.45           O
ATOM    302  N   TYR A  44      13.725  11.174  16.425  1.00  5.22           N
ATOM    303  CA  TYR A  44      13.257  10.745  15.081  1.00  5.56           C
ATOM    304  C   TYR A  44      14.275   9.687  14.612  1.00  4.61           C
ATOM    305  O   TYR A  44      14.930   9.862  13.568  1.00  6.04           O
ATOM    306  CB  TYR A  44      13.200  11.914  14.071  1.00  5.41           C
ATOM    307  CG  TYR A  44      12.000  12.819  14.399  1.00  5.34           C
ATOM    308  CD1 TYR A  44      12.119  13.853  15.332  1.00  6.59           C
ATOM    309  CD2 TYR A  44      10.775  12.617  13.762  1.00  5.94           C
ATOM    310  CE1 TYR A  44      11.045  14.675  15.610  1.00  5.97           C
ATOM    311  CE2 TYR A  44       9.676  13.433  14.048  1.00  5.17           C
ATOM    312  CZ  TYR A  44       9.802  14.456  14.996  1.00  5.96           C
ATOM    313  OH  TYR A  44       8.740  15.265  15.269  1.00  8.60           O
ATOM    314  N   ALA A  45      14.342   8.640  15.422  1.00  4.76           N
ATOM    315  CA  ALA A  45      15.445   7.667  15.246  1.00  5.89           C
ATOM    316  C   ALA A  45      15.171   6.533  14.280  1.00  6.67           C
ATOM    317  O   ALA A  45      16.093   5.705  14.039  1.00  7.56           O
ATOM    318  CB  ALA A  45      15.680   7.099  16.682  1.00  6.82           C
ATOM    319  N   ASN A  46      13.966   6.502  13.739  1.00  5.80           N
ATOM    320  CA  ASN A  46      13.512   5.395  12.878  1.00  6.15           C
ATOM    321  C   ASN A  46      13.311   5.853  11.455  1.00  6.61           C
ATOM    322  O   ASN A  46      13.733   6.929  11.026  1.00  7.18           O
ATOM    323  CB  ASN A  46      12.266   4.769  13.501  1.00  7.27           C
ATOM    324  CG  ASN A  46      12.538   4.304  14.922  1.00  7.98           C
ATOM    325  OD1 ASN A  46      11.982   4.849  15.886  1.00 11.00           O
ATOM    326  ND2 ASN A  46      13.407   3.298  15.015  1.00 10.32           N
ATOM    327  OXT ASN A  46      12.703   4.973  10.746  1.00  7.86           O
TER     328      ASN A  46
END
";

/// Which built-in demo structure the molecular viewer loads — picked by the
/// "Demo structure" combo and settable by an agent. [`Crambin`](Self::Crambin)
/// (PDB 1CRN, a real 46-residue mixed α/β protein) is the default; the tiny
/// glycine peptide stays as the trivial smoke case.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub(crate) enum DemoStructure {
    /// Real 46-residue, 327-atom crambin (PDB 1CRN) — the default. Has both an
    /// α-helix and a β-sheet, so the ribbon / cartoon + DSSP colour scheme show
    /// real secondary structure.
    #[default]
    Crambin,
    /// The 12-atom triglycine peptide ([`DEMO_PDB`]) — a minimal smoke
    /// structure that intentionally looks like almost nothing.
    GlyPeptide,
}

impl DemoStructure {
    /// Every demo structure, in picker order.
    pub(crate) const ALL: [DemoStructure; 2] = [DemoStructure::Crambin, DemoStructure::GlyPeptide];

    /// The user-visible / agent-facing caption for this demo structure.
    pub(crate) fn label(self) -> &'static str {
        match self {
            DemoStructure::Crambin => "Crambin (1CRN, 46 res)",
            DemoStructure::GlyPeptide => "Triglycine peptide (3 res)",
        }
    }

    /// The embedded PDB text for this demo structure — parses via
    /// [`read_structure`].
    pub(crate) fn pdb(self) -> &'static str {
        match self {
            DemoStructure::Crambin => CRAMBIN_PDB,
            DemoStructure::GlyPeptide => DEMO_PDB,
        }
    }

    /// Resolve a demo structure from its display [`label`](Self::label),
    /// case-insensitively. Accepts a few short aliases (`"crambin"`, `"1crn"`,
    /// `"glycine"`, `"peptide"`) so an agent need not quote the full caption.
    /// Fail-loud: an unrecognised name returns `Err`.
    pub(crate) fn from_label(s: &str) -> Result<DemoStructure, String> {
        let n = s.trim().to_ascii_lowercase();
        DemoStructure::ALL
            .into_iter()
            .find(|d| d.label().to_ascii_lowercase() == n)
            .or(match n.as_str() {
                "crambin" | "1crn" | "crn" => Some(DemoStructure::Crambin),
                "glycine" | "triglycine" | "peptide" | "gly" => Some(DemoStructure::GlyPeptide),
                _ => None,
            })
            .ok_or_else(|| format!("unknown demo structure: {s:?}"))
    }
}

impl Default for BiostructPanel {
    fn default() -> Self {
        BiostructPanel {
            tool: Tool::Analyze,
            // Default to the real 46-residue crambin so the viewer shows an
            // actual protein (ribbon + secondary structure), not the 12-atom
            // smoke peptide. `structure_b` (the superpose reference) starts as
            // the same structure so a default RMSD run is well-defined.
            structure_a: CRAMBIN_PDB.to_string(),
            structure_b: CRAMBIN_PDB.to_string(),
            clash_tolerance: 0.4,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
            representation: Representation::default(),
            color_scheme: ColorScheme::default(),
            molviz_params: MolvizParams::default(),
            trajectory: TrajectoryPlayback::default(),
            md_temperature: EnmParams::default().temperature,
            demo: DemoStructure::default(),
        }
    }
}

/// Collect every Cα atom coordinate from the first model of a structure.
fn ca_coords(s: &Structure) -> Vec<Point3<f64>> {
    let mut out = Vec::new();
    for chain in &s.first_model().chains {
        for res in &chain.residues {
            if let Some(ca) = res.ca() {
                out.push(ca.coord);
            }
        }
    }
    out
}

/// Render the Macromolecular Structure panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.biostruct;

    common::section(ui, "Tool");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Analyze, "Structure analysis")
            .on_hover_text("Detect secondary structure, contacts, clashes, and chains.");
        ui.selectable_value(&mut p.tool, Tool::Ramachandran, "Ramachandran")
            .on_hover_text(
                "Compute φ/ψ backbone dihedrals and classify into Ramachandran regions.",
            );
        ui.selectable_value(&mut p.tool, Tool::Superpose, "RMSD / superpose")
            .on_hover_text("Superpose two structures via Kabsch rotation + report RMSD.");
        ui.separator();
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u {
            p.undo_edit();
        }
        if r {
            p.redo_edit();
        }
    });
    ui.separator();

    match p.tool {
        Tool::Analyze | Tool::Ramachandran => draw_single(p, ui),
        Tool::Superpose => draw_superpose(p, ui),
    }

    common::error_line(ui, &p.error);

    // --- 3-D viewport integration ---------------------------------
    // The structure-A text feeds the viewer for every tool (it is the
    // analyse / Ramachandran input and the mobile superpose
    // structure). A reactive representation picker selects how it is
    // meshed (ball-and-stick / sticks / spacefill / cartoon / surface),
    // and "Show in 3D viewport" pushes that mesh into the app's wgpu 3-D
    // viewport. The picker is a row of named `selectable_value` widgets,
    // so an agent can switch representation through the accessibility
    // tree by widget name.
    if !app.genetics.biostruct.structure_a.trim().is_empty() {
        draw_representation_picker(app, ui);
        if ui.button("Show in 3D viewport").clicked() {
            show_in_viewport(app);
        }
        draw_trajectory_controls(app, ui);
    }

    let p = &app.genetics.biostruct;
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "biostruct_result", &p.result, 16);
    }
}

/// Draw the reactive **representation picker** for the 3-D viewer: a row of
/// named `selectable_value` buttons (one per [`Representation`]) plus the
/// per-representation tuning controls. The named widgets are what makes the
/// representation AI-drivable — an agent flips representation by selecting the
/// button by its `label()` through the accessibility tree.
fn draw_representation_picker(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.biostruct;
    common::section(ui, "3-D representation");
    ui.horizontal_wrapped(|ui| {
        for rep in Representation::ALL {
            ui.selectable_value(&mut p.representation, rep, rep.label())
                .on_hover_text(match rep {
                    Representation::BallAndStick => "Element spheres + bond cylinders.",
                    Representation::Sticks => "Bonds only (licorice).",
                    Representation::Spacefill => "Full van-der-Waals spheres (CPK).",
                    Representation::Cartoon => {
                        "Smooth Catmull-Rom round tube through the Cα backbone \
                         (helices/strands fatten via DSSP)."
                    }
                    Representation::Ribbon => {
                        "Flat elliptical ribbon swept along the Cα spline — a wide, \
                         thin band (helices/strands widen via DSSP)."
                    }
                    Representation::Surface => "Marching-cubes molecular surface (union-of-balls).",
                    Representation::Density => {
                        "Marching-cubes isosurface of a Gaussian electron-density-like \
                         field (sum of per-atom Gaussians; a smooth-blob model, not real QM)."
                    }
                });
        }
    });

    // Colour-scheme picker: a row of named `selectable_value` buttons (one per
    // [`ColorScheme`]). A non-default scheme (chain / residue / B-factor) makes
    // "Show in 3D viewport" build a per-vertex-coloured mesh and upload the
    // colours so the viewport renders the scheme. Each button is `labelled_by`
    // the "Colour scheme" caption so it carries an unambiguous accessibility /
    // UI-Automation Name (AI-drivable: an agent flips the scheme by selecting
    // the button by its `label()`), and the caption itself names the group.
    ui.horizontal_wrapped(|ui| {
        let caption = ui.label("Colour scheme");
        for scheme in ColorScheme::ALL {
            ui.selectable_value(&mut p.color_scheme, scheme, scheme.label())
                .labelled_by(caption.id)
                .on_hover_text(match scheme {
                    ColorScheme::Element => "CPK by element (C grey, N blue, O red, S yellow, …).",
                    ColorScheme::Chain => "A distinct hue per chain.",
                    ColorScheme::Residue => "Rainbow ramp by residue index (N→C terminus).",
                    ColorScheme::BFactor => "Blue→white→red ramp by B-factor (low→high).",
                    ColorScheme::SecondaryStructure => {
                        "By DSSP secondary structure: helix magenta-red, sheet \
                         yellow, coil/loop grey."
                    }
                });
        }
    });
    if p.color_scheme != ColorScheme::Element {
        // Sticks / cartoon / ribbon / surface / density have no per-atom colour
        // builder, so they take a single scheme-derived tint rather than
        // per-atom colour. Tell the user which they'll get.
        let per_atom = matches!(
            p.representation,
            Representation::BallAndStick | Representation::Spacefill
        );
        if !per_atom {
            ui.label(
                egui::RichText::new(format!(
                    "note: {} renders a single {}-derived tint (per-atom colour needs \
                     ball-and-stick or spacefill)",
                    p.representation.label().to_ascii_lowercase(),
                    p.color_scheme.label().to_ascii_lowercase(),
                ))
                .weak()
                .italics(),
            );
        }
    }

    // Per-representation tuning, only shown when relevant.
    match p.representation {
        Representation::Surface => {
            // Probe-based surface mode: a row of named `selectable_value` buttons
            // (vdW / SAS / SES), each `labelled_by` the "Surface" caption so it
            // carries an unambiguous accessibility / UI-Automation Name — an agent
            // flips the surface type by selecting the button by its `label()`.
            ui.horizontal_wrapped(|ui| {
                let caption = ui.label("Surface type");
                for mode in molviz::SurfaceMode::ALL {
                    ui.selectable_value(&mut p.molviz_params.surface_mode, mode, mode.label())
                        .labelled_by(caption.id)
                        .on_hover_text(match mode {
                            molviz::SurfaceMode::Vdw => {
                                "Van-der-Waals surface: union of the bare atom spheres \
                                 (probe ignored)."
                            }
                            molviz::SurfaceMode::Sas => {
                                "Solvent-accessible surface (Lee–Richards): union of \
                                 probe-inflated spheres — the path of the probe centre."
                            }
                            molviz::SurfaceMode::Ses => {
                                "Solvent-excluded / Connolly surface (Richards): the smooth \
                                 re-entrant surface, built by eroding the SAS solid by the \
                                 probe radius on the grid."
                            }
                        });
                }
            });
            ui.horizontal(|ui| {
                ui.label("Surface grid:");
                ui.add(egui::Slider::new(&mut p.molviz_params.grid_max, 16..=128).text("cells"))
                    .on_hover_text(
                        "Marching-cubes resolution along the longest axis — higher is \
                     smoother but costs O(n³). SES re-entrant detail sharpens with resolution.",
                    );
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.probe_radius)
                        .speed(0.05)
                        .range(0.0..=3.0)
                        .prefix("probe Å "),
                )
                .on_hover_text(
                    "Rolling-probe radius (1.4 Å ≈ water). Inflates each atom for SAS \
                     and is the erosion radius for SES; ignored for vdW.",
                );
            });
        }
        Representation::Density => {
            ui.horizontal(|ui| {
                ui.label("Density grid:");
                ui.add(
                    egui::Slider::new(&mut p.molviz_params.density_grid_max, 16..=128)
                        .text("cells"),
                )
                .on_hover_text(
                    "Marching-cubes resolution along the longest axis — higher is \
                     smoother but costs O(n³).",
                );
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.density_sigma)
                        .speed(0.05)
                        .range(0.2..=4.0)
                        .prefix("σ Å "),
                )
                .on_hover_text(
                    "Gaussian width per atom — larger σ gives fatter, more-merged blobs.",
                );
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.density_iso)
                        .speed(0.02)
                        .range(0.05..=0.95)
                        .prefix("iso "),
                )
                .on_hover_text(
                    "Iso-level as a fraction of one atom's peak — lower wraps more of \
                     each Gaussian tail. (Phenomenological blob, not real electron density.)",
                );
            });
            ui.checkbox(
                &mut p.molviz_params.density_weight_by_element,
                "Weight density by element",
            )
            .on_hover_text(
                "Scale each atom's Gaussian by a crude electron count so heavy atoms read denser.",
            );
        }
        Representation::Cartoon => {
            ui.horizontal(|ui| {
                ui.label("Tube radius (Å):");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.tube_radius)
                        .speed(0.02)
                        .range(0.1..=2.0),
                );
            });
        }
        Representation::Ribbon => {
            ui.horizontal(|ui| {
                ui.label("Ribbon width (Å):");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.ribbon_width)
                        .speed(0.05)
                        .range(0.2..=4.0),
                )
                .on_hover_text("Half-width of the flat band (the wide axis).");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.ribbon_thickness)
                        .speed(0.02)
                        .range(0.05..=1.0)
                        .prefix("thick Å "),
                )
                .on_hover_text("Half-thickness of the band (the thin axis).");
            });
        }
        Representation::BallAndStick | Representation::Sticks => {
            ui.horizontal(|ui| {
                ui.label("Stick radius (Å):");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.stick_radius)
                        .speed(0.01)
                        .range(0.02..=0.6),
                );
            });
        }
        Representation::Spacefill => {}
    }
}

/// Build the structure-A text into the selected [`Representation`]'s mesh and
/// push it into the app's 3-D viewport. For the cartoon, the Cα backbone and a
/// per-residue DSSP secondary-structure track are extracted first; for the
/// other modes the atom/bond [`ViewMolecule`] is meshed.
fn show_in_viewport(app: &mut ValenxApp) {
    let rep = app.genetics.biostruct.representation;
    let scheme = app.genetics.biostruct.color_scheme;
    let params = app.genetics.biostruct.molviz_params;
    match read_structure(&app.genetics.biostruct.structure_a, "viewer") {
        Ok(s) => {
            let backbone = if rep.needs_backbone() {
                ca_backbone(&s)
            } else {
                Vec::new()
            };
            // A cartoon / ribbon needs a real Cα trace — fail loudly rather
            // than silently showing nothing.
            if rep.needs_backbone() && backbone.len() < 2 {
                app.genetics.biostruct.error = Some(format!(
                    "{} needs a protein backbone (≥ 2 Cα atoms) — this \
                     structure has none; try ball-and-stick or surface",
                    rep.label().to_ascii_lowercase()
                ));
                return;
            }
            let view = ViewMolecule::from_biostruct(&s);
            let label = format!("structure.{}.{}", rep.token(), scheme.token());
            // Default (Element / CPK) keeps the original monochrome-metal path
            // (`show_molecule`). A non-default scheme builds the colour-aware
            // mesh and uploads per-vertex colours so the viewport renders the
            // scheme instead of a single material.
            let result = if scheme == ColorScheme::Element {
                let mesh = molviz::build_mesh(&view, rep, &backbone, &params);
                molecule_view::show_molecule(app, mesh, &label)
            } else {
                // Per-atom annotations (chain / residue / B-factor) in lockstep
                // with `view.atoms`, read from the structure the viewer parsed.
                let attrs = structure_atom_attrs(&s);
                let (mesh, per_tri_colors) =
                    molviz::build_mesh_colored(&view, rep, &backbone, &params, scheme, &attrs);
                molecule_view::show_molecule_colored(app, mesh, &per_tri_colors, &label)
            };
            match result {
                Ok(_) => app.genetics.biostruct.error = None,
                Err(e) => app.genetics.biostruct.error = Some(e),
            }
        }
        Err(e) => app.genetics.biostruct.error = Some(e.to_string()),
    }
}

/// Draw the **MD-trajectory playback** controls (the MolecularNodes-style
/// animation) and run the per-frame update.
///
/// All controls carry an accessibility Name and are `labelled_by` their caption
/// so the playback is AI-drivable (an agent finds the frame slider / Play /
/// speed by name in the accessibility tree). When playing, the panel advances
/// the frame index by `speed × dt` each repaint and requests another repaint so
/// the animation is self-driving from inside the panel's `update`.
fn draw_trajectory_controls(app: &mut ValenxApp, ui: &mut egui::Ui) {
    common::section(ui, "MD trajectory playback");

    // --- attach buttons -------------------------------------------------
    ui.horizontal_wrapped(|ui| {
        if ui
            .button("Run MD (ENM)")
            .on_hover_text(
                "In-house Elastic Network Model molecular dynamics: harmonic \
                 springs between nearby atoms, Maxwell-Boltzmann thermal \
                 velocities, velocity-Verlet integration. Produces a real, \
                 deterministic trajectory of the protein vibrating with thermal \
                 motion — press Play to animate the breathing.",
            )
            .clicked()
        {
            // Result is surfaced via the trajectory note / readout below.
            let _ = app.genetics.biostruct.run_md();
            if app.genetics.biostruct.trajectory.is_attached() {
                render_trajectory_frame(app);
            }
        }
        // AI-drivable MD temperature: a labelled slider an agent can set by
        // name (`"MD temperature"`) to scale the vibration amplitude.
        let temp_cap = ui.label("MD temperature");
        ui.add(egui::Slider::new(&mut app.genetics.biostruct.md_temperature, 0.0..=1.5).text("T"))
            .labelled_by(temp_cap.id)
            .on_hover_text(
                "Reduced-unit temperature for Run MD (ENM): higher ⇒ larger thermal \
             vibration. Re-run MD after changing it.",
            );
        if ui
            .button("Generate demo trajectory")
            .on_hover_text(
                "Build a small synthetic wobble of the current atoms over 60 \
                 frames — a real, deterministic trajectory to play with no \
                 external file.",
            )
            .clicked()
        {
            attach_synthetic_trajectory(app);
        }
        if ui
            .button("Load valenx-md trajectory…")
            .on_hover_text(
                "Load a binary (.dcd-class) or framed-text valenx-md trajectory \
                 and animate the current structure across its frames (atom counts \
                 must match).",
            )
            .clicked()
        {
            load_md_trajectory(app);
        }
        if app.genetics.biostruct.trajectory.is_attached()
            && ui.button("Clear trajectory").clicked()
        {
            app.genetics.biostruct.trajectory.clear();
        }
    });

    // --- playback transport (only once a trajectory is attached) --------
    if app.genetics.biostruct.trajectory.is_attached() {
        let n_frames = app.genetics.biostruct.trajectory.n_frames();
        let last = n_frames.saturating_sub(1);

        // Advance the clock while playing BEFORE drawing the slider so the
        // slider reflects the frame we are about to render this repaint.
        if app.genetics.biostruct.trajectory.playing && n_frames > 1 {
            // `stable_dt` is the real frame time; fall back to ~60 FPS on the
            // first frame (when egui reports 0) so playback still advances.
            let dt = {
                let d = ui.ctx().input(|i| i.stable_dt);
                if d.is_finite() && d > 0.0 {
                    d.min(0.1) // clamp a long stall so we don't jump the whole reel
                } else {
                    1.0 / 60.0
                }
            };
            let tj = &mut app.genetics.biostruct.trajectory;
            tj.accum += tj.speed * dt;
            while tj.accum >= 1.0 {
                tj.accum -= 1.0;
                tj.frame = (tj.frame + 1) % n_frames; // loop
            }
        }

        ui.horizontal(|ui| {
            let caption = ui.label("Frame");
            let mut frame = app.genetics.biostruct.trajectory.frame;
            let changed = ui
                .add(egui::Slider::new(&mut frame, 0..=last).clamp_to_range(true))
                .labelled_by(caption.id)
                .changed();
            let dv_changed = ui
                .add(egui::DragValue::new(&mut frame).range(0..=last))
                .labelled_by(caption.id)
                .changed();
            ui.label(format!("/ {n_frames}"));
            if changed || dv_changed {
                // Scrubbing pauses auto-play so the user keeps control.
                app.genetics.biostruct.trajectory.playing = false;
            }
            app.genetics.biostruct.trajectory.frame = frame.min(last);
        });

        ui.horizontal(|ui| {
            let caption = ui.label("Transport");
            let playing = app.genetics.biostruct.trajectory.playing;
            let label = if playing { "Pause" } else { "Play" };
            if ui
                .add(egui::Button::new(label))
                .labelled_by(caption.id)
                .on_hover_text("Start / stop animating the structure across frames.")
                .clicked()
            {
                app.genetics.biostruct.trajectory.playing = !playing;
            }
            if ui
                .add(egui::Button::new("⏮ Reset"))
                .labelled_by(caption.id)
                .on_hover_text("Jump back to frame 0.")
                .clicked()
            {
                app.genetics.biostruct.trajectory.frame = 0;
                app.genetics.biostruct.trajectory.accum = 0.0;
            }
            ui.separator();
            let speed_caption = ui.label("Speed (fps)");
            ui.add(
                egui::Slider::new(&mut app.genetics.biostruct.trajectory.speed, 0.5..=60.0)
                    .text("fps"),
            )
            .labelled_by(speed_caption.id)
            .on_hover_text("Frames advanced per second during playback.");
        });

        let tj = &app.genetics.biostruct.trajectory;
        ui.label(
            egui::RichText::new(format!(
                "{} · frame {}/{}",
                tj.source,
                tj.frame + 1,
                n_frames
            ))
            .weak(),
        );

        // Render the current frame into the viewport every repaint (so the
        // slider and the auto-advancing clock both show live geometry).
        render_trajectory_frame(app);

        // Keep the animation alive: a playing trajectory must keep repainting
        // even with no input events.
        if app.genetics.biostruct.trajectory.playing {
            ui.ctx().request_repaint();
        }
    }

    // Surface any attach-time / playback note (atom-count mismatch, no atoms…).
    if let Some(note) = &app.genetics.biostruct.trajectory.note {
        ui.label(
            egui::RichText::new(format!("trajectory: {note}"))
                .italics()
                .color(egui::Color32::from_rgb(0xB0, 0x60, 0x20)),
        );
    }
}

/// Parse the current structure text and attach a synthetic-wobble trajectory.
fn attach_synthetic_trajectory(app: &mut ValenxApp) {
    match read_structure(&app.genetics.biostruct.structure_a, "viewer") {
        Ok(s) => {
            app.genetics.biostruct.trajectory.speed = 12.0;
            // 60 frames, 0.4 Å wobble — enough motion to read on screen.
            app.genetics
                .biostruct
                .trajectory
                .attach_synthetic(&s, 60, 0.4);
            // Show frame 0 immediately so the viewport isn't stale.
            if app.genetics.biostruct.trajectory.is_attached() {
                render_trajectory_frame(app);
            }
        }
        Err(e) => {
            app.genetics.biostruct.trajectory.note =
                Some(format!("could not parse structure: {e}"));
        }
    }
}

/// Pick a `valenx-md` trajectory file (binary or framed-text), parse it, and
/// attach it to the current structure. Fail-loud: a bad file / atom-count
/// mismatch sets an in-panel note, never panics.
fn load_md_trajectory(app: &mut ValenxApp) {
    use valenx_md::io::trajectory::{read_binary, read_text};

    let Some(path) = rfd::FileDialog::new()
        .add_filter("valenx-md trajectory", &["dcd", "trj", "traj", "txt"])
        .pick_file()
    else {
        return;
    };
    // Cap the read so a huge picked file can't OOM the renderer before parsing.
    let bytes = match valenx_core::io_caps::read_capped_to_bytes(
        &path,
        valenx_core::io_caps::MAX_GENETICS_FILE_BYTES,
    ) {
        Ok(b) => b,
        Err(e) => {
            app.genetics.biostruct.trajectory.note = Some(format!("read: {e}"));
            return;
        }
    };
    // Try the binary format first, then fall back to the framed-text format.
    let traj = read_binary(&bytes).or_else(|_| {
        std::str::from_utf8(&bytes)
            .map_err(|_| valenx_md::error::MdError::parse("trajectory", "not UTF-8 text"))
            .and_then(read_text)
    });
    match traj {
        Ok(traj) => match read_structure(&app.genetics.biostruct.structure_a, "viewer") {
            Ok(s) => {
                app.genetics.biostruct.trajectory.attach_md(&s, &traj);
                if app.genetics.biostruct.trajectory.is_attached() {
                    render_trajectory_frame(app);
                }
            }
            Err(e) => {
                app.genetics.biostruct.trajectory.note =
                    Some(format!("could not parse structure: {e}"));
            }
        },
        Err(e) => {
            app.genetics.biostruct.trajectory.note =
                Some(format!("could not read trajectory: {e}"));
        }
    }
}

/// Mesh the current trajectory frame with the panel's selected representation +
/// colour scheme and push it into the 3-D viewport.
///
/// This is the per-frame update the slider / playback drive: it overwrites the
/// base molecule's atom positions with the current frame (via
/// [`TrajectoryPlayback::molecule_at`], which also re-detects bonds), rebuilds
/// the representation mesh, and shows it through the same
/// [`molecule_view::show_molecule`] / [`molecule_view::show_molecule_colored`]
/// path as the static "Show in 3D viewport" button — so a frame renders
/// identically to a freshly-shown structure at those coordinates.
///
/// Fail-loud, never panic: nothing attached / a short frame is a silent no-op;
/// a backbone representation with too few Cα points sets the in-panel note (the
/// secondary-structure backbone is the reference structure's, static across
/// frames).
fn render_trajectory_frame(app: &mut ValenxApp) {
    let rep = app.genetics.biostruct.representation;
    let scheme = app.genetics.biostruct.color_scheme;
    let params = app.genetics.biostruct.molviz_params;
    let frame = app.genetics.biostruct.trajectory.frame;

    let Some(view) = app.genetics.biostruct.trajectory.molecule_at(frame) else {
        return; // nothing attached, or a guarded short frame
    };
    let backbone = if rep.needs_backbone() {
        app.genetics.biostruct.trajectory.backbone.clone()
    } else {
        Vec::new()
    };
    if rep.needs_backbone() && backbone.len() < 2 {
        app.genetics.biostruct.trajectory.note = Some(format!(
            "{} needs a protein backbone (≥ 2 Cα atoms) — this structure has \
             none; pick ball-and-stick or surface to animate it",
            rep.label().to_ascii_lowercase()
        ));
        return;
    }
    let label = format!("trajectory.{}.{}.frame{frame}", rep.token(), scheme.token());
    let result = if scheme == ColorScheme::Element {
        let mesh = molviz::build_mesh(&view, rep, &backbone, &params);
        molecule_view::show_molecule(app, mesh, &label)
    } else {
        let attrs = &app.genetics.biostruct.trajectory.attrs;
        let (mesh, per_tri_colors) =
            molviz::build_mesh_colored(&view, rep, &backbone, &params, scheme, attrs);
        molecule_view::show_molecule_colored(app, mesh, &per_tri_colors, &label)
    };
    if let Err(e) = result {
        app.genetics.biostruct.trajectory.note = Some(e);
    }
}

/// Build the per-atom [`AtomAttr`] (chain id, residue index, B-factor,
/// secondary structure) for a structure's first model, **in the exact atom
/// order** [`ViewMolecule::from_biostruct`] walks (chain-major, then residue,
/// then atom) so the returned vec is in lockstep with `ViewMolecule::atoms` and
/// the colour schemes line up atom-for-atom.
///
/// Chain id, residue index and B-factor come straight from the
/// `valenx-biostruct` model, which carries them on every atom (`Atom::b_factor`,
/// `Chain::id`, `Residue::seq_num`): no field needs a fallback. The residue
/// index is a **monotone counter across the whole model** (incremented per
/// residue in iteration order) rather than the raw `seq_num`, so the residue
/// rainbow runs cleanly N→C even when `seq_num` has gaps, insertion codes, or
/// resets between chains.
///
/// The **secondary structure** is the same per-chain DSSP assignment the
/// cartoon backbone uses ([`valenx_biostruct::dssp::assign_chain`]), keyed by
/// residue position: every atom of residue `i` in a chain inherits that chain's
/// `states[i]`, collapsed to the three-state [`SsKind`]. A residue past the end
/// of the DSSP track (defensive — `assign_chain` covers every residue) gets
/// `ss = None`, which the SS colour scheme renders as coil grey (fail-loud, no
/// panic).
fn structure_atom_attrs(s: &Structure) -> Vec<AtomAttr> {
    use crate::molviz::SsKind;
    use valenx_biostruct::dssp;

    let mut out = Vec::new();
    let mut residue_index: i32 = 0;
    for chain in &s.first_model().chains {
        // Per-chain DSSP state, one entry per residue in chain order — the same
        // assignment `ca_backbone` tags the cartoon Cα trace with, so the
        // per-atom SS colour and the cartoon's tube modulation agree.
        let states = dssp::assign_chain(chain).states;
        for (i, res) in chain.residues.iter().enumerate() {
            let ss = states.get(i).map(|st| SsKind::from_dssp_code(st.code()));
            for atom in &res.atoms {
                out.push(
                    AtomAttr::new(chain.id.clone(), residue_index, atom.b_factor as f32)
                        .with_ss(ss),
                );
            }
            residue_index += 1;
        }
    }
    out
}

/// Extract the per-residue Cα backbone control points of a structure's first
/// model, tagged with their DSSP secondary-structure code, for the cartoon
/// representation. Residues without a Cα (ligands, waters, missing atoms) are
/// skipped; the DSSP track is computed per chain via [`valenx_biostruct::dssp`]
/// and indexed by residue position so each kept Cα carries its own state.
fn ca_backbone(s: &valenx_biostruct::structure::Structure) -> Vec<BackbonePoint> {
    use valenx_biostruct::dssp;
    let mut out = Vec::new();
    for chain in &s.first_model().chains {
        // Per-chain DSSP state, one entry per residue in chain order.
        let states = dssp::assign_chain(chain).states;
        for (i, res) in chain.residues.iter().enumerate() {
            if let Some(ca) = res.ca() {
                let ss = states.get(i).map(|st| st.code());
                out.push(BackbonePoint::new(
                    [ca.coord.x as f32, ca.coord.y as f32, ca.coord.z as f32],
                    ss,
                ));
            }
        }
    }
    out
}

fn structure_text_input(
    ui: &mut egui::Ui,
    id: &str,
    label: &str,
    buf: &mut String,
) -> Option<String> {
    let mut err = None;
    common::section(ui, label);
    ui.horizontal(|ui| {
        if ui.small_button("Load PDB / mmCIF…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("Structure", &["pdb", "ent", "cif", "mmcif"])
                .pick_file()
            {
                // Round-21 H1: file-dialog paths flow straight to a
                // bare `fs::read_to_string` pre-fix. A user (or a
                // stale dialog state) pointing at a multi-GB file
                // would OOM the renderer before the parser saw a
                // single byte. `read_capped_to_string` rejects
                // anything past `MAX_GENETICS_FILE_BYTES` (64 MiB).
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => *buf = t,
                    Err(e) => err = Some(format!("read: {e}")),
                }
            }
        }
        ui.label(format!("{} lines", buf.lines().count()));
    });
    ui.add(
        egui::TextEdit::multiline(buf)
            .id_source(id)
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(6),
    );
    err
}

/// Draw the **"Demo structure"** picker — a labelled row of `selectable_value`
/// chips, one per built-in [`DemoStructure`]. Selecting one loads its embedded
/// PDB into `structure_a` (via [`BiostructPanel::load_demo`]). Each chip is
/// `labelled_by` its caption and carries the demo's display name as its
/// accessibility Name (the same always-visible pattern as the representation /
/// colour-scheme pickers), so an agent can switch the loaded structure through
/// the accessibility tree by name without first opening a popup; the same path
/// backs the host `agent_set`.
fn demo_structure_picker(p: &mut BiostructPanel, ui: &mut egui::Ui) {
    ui.horizontal_wrapped(|ui| {
        let cap = ui.label("Demo structure");
        for d in DemoStructure::ALL {
            let selected = p.demo == d;
            let resp = ui
                .selectable_label(selected, d.label())
                .on_hover_text(match d {
                    DemoStructure::Crambin => {
                        "Real 46-residue protein (PDB 1CRN) with an α-helix and a \
                         β-sheet — shows secondary structure in the ribbon view."
                    }
                    DemoStructure::GlyPeptide => "Tiny 12-atom triglycine smoke peptide.",
                })
                .labelled_by(cap.id);
            if resp.clicked() && !selected {
                p.load_demo(d);
            }
        }
    });
}

fn draw_single(p: &mut BiostructPanel, ui: &mut egui::Ui) {
    demo_structure_picker(p, ui);
    if let Some(e) = structure_text_input(
        ui,
        "biostruct_input_a",
        "Structure (PDB / mmCIF)",
        &mut p.structure_a,
    ) {
        p.error = Some(e);
    }
    if p.tool == Tool::Analyze {
        ui.horizontal(|ui| {
            ui.label("Clash tolerance (Å):");
            ui.add(
                egui::DragValue::new(&mut p.clash_tolerance)
                    .speed(0.05)
                    .range(0.0..=2.0),
            );
        });
    }
    let label = if p.tool == Tool::Analyze {
        "Analyze structure"
    } else {
        "Classify Ramachandran"
    };
    if common::run_button(ui, label) {
        let snap = p.snapshot();
        p.history.record(snap);
        run_single(p);
    }
}

/// Run the structure-analysis / Ramachandran tool — extracted from the
/// button closure so it is callable from the headless UI tests.
fn run_single(p: &mut BiostructPanel) {
    p.error = None;
    match read_structure(&p.structure_a, "input") {
        Ok(s) => match p.tool {
            Tool::Analyze => match StructureReport::analyze(&s, p.clash_tolerance) {
                Ok(r) => {
                    let mut out = format!(
                        "title          : {}\nmodels         : {}\n\
                             chains         : {}  ({} protein, {} nucleic)\n\
                             residues       : {}\natoms          : {}\n\
                             water / ligand : {} / {}\nradius of gyr. : {:.2} Å\n\
                             mean helix     : {:.1} %\nmean sheet     : {:.1} %\n\n\
                             -- per chain --\n",
                        r.title,
                        r.model_count,
                        r.chains.len(),
                        r.protein_chain_count(),
                        r.nucleic_chain_count(),
                        r.residue_count,
                        r.atom_count,
                        r.water_count,
                        r.ligand_count,
                        r.radius_of_gyration,
                        r.mean_helix_fraction() * 100.0,
                        r.mean_sheet_fraction() * 100.0,
                    );
                    for ch in &r.chains {
                        out.push_str(&format!(
                            "  {} {:<10} {:>4} res  H {:>4.0}% E {:>4.0}% C {:>4.0}%\n",
                            ch.id,
                            format!("{:?}", ch.kind),
                            ch.residue_count,
                            ch.secondary.helix * 100.0,
                            ch.secondary.sheet * 100.0,
                            ch.secondary.coil * 100.0,
                        ));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            },
            Tool::Ramachandran => {
                let mut out = String::new();
                for chain in &s.first_model().chains {
                    let summary = rama_summarize(chain);
                    if summary.total == 0 {
                        continue;
                    }
                    out.push_str(&format!(
                        "chain {} — {} phi/psi points\n  alpha-helix : {}\n  \
                             beta-sheet  : {}\n  left-alpha  : {}\n  bridge      : {}\n  \
                             outliers    : {}\n  allowed     : {:.1} %\n\n",
                        chain.id,
                        summary.total,
                        summary.alpha,
                        summary.beta,
                        summary.left_alpha,
                        summary.bridge,
                        summary.outliers,
                        summary.allowed_fraction() * 100.0,
                    ));
                }
                if out.is_empty() {
                    out = "no residues with a defined φ/ψ pair (need ≥ 3 \
                               consecutive amino acids)"
                        .to_string();
                }
                p.result = out;
            }
            Tool::Superpose => unreachable!(),
        },
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_superpose(p: &mut BiostructPanel, ui: &mut egui::Ui) {
    demo_structure_picker(p, ui);
    if let Some(e) = structure_text_input(
        ui,
        "biostruct_input_mob",
        "Mobile structure",
        &mut p.structure_a,
    ) {
        p.error = Some(e);
    }
    if let Some(e) = structure_text_input(
        ui,
        "biostruct_input_ref",
        "Reference structure",
        &mut p.structure_b,
    ) {
        p.error = Some(e);
    }
    if common::run_button(ui, "Kabsch superpose (Cα)") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_superpose(p);
    }
}

/// Run the Kabsch Cα superposition — extracted for the headless UI
/// tests.
fn run_superpose(p: &mut BiostructPanel) {
    p.error = None;
    match (
        read_structure(&p.structure_a, "mobile"),
        read_structure(&p.structure_b, "reference"),
    ) {
        (Ok(mob), Ok(reference)) => {
            let ca_m = ca_coords(&mob);
            let ca_r = ca_coords(&reference);
            let n = ca_m.len().min(ca_r.len());
            if n < 3 {
                p.error = Some(format!(
                    "need ≥ 3 paired Cα atoms (mobile {}, reference {})",
                    ca_m.len(),
                    ca_r.len(),
                ));
                return;
            }
            let m = &ca_m[..n];
            let r = &ca_r[..n];
            let pre = rmsd(m, r);
            match kabsch(m, r) {
                Ok(sup) => {
                    p.result = format!(
                        "paired Cα atoms : {}\nRMSD before     : {}\n\
                         RMSD after fit  : {:.4} Å\n\n\
                         optimal rotation + translation found by the \
                         Kabsch algorithm.",
                        n,
                        pre.map(|v| format!("{v:.4} Å"))
                            .unwrap_or_else(|_| "(n/a)".into()),
                        sup.rmsd,
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        (Err(e), _) => p.error = Some(format!("mobile: {e}")),
        (_, Err(e)) => p.error = Some(format!("reference: {e}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_pdb_parses() {
        let s = read_structure(DEMO_PDB, "demo").expect("demo PDB must parse");
        assert_eq!(ca_coords(&s).len(), 3);
    }

    #[test]
    fn demo_superposes_to_zero_rmsd() {
        // The same structure superposed on itself has ~zero RMSD.
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let ca = ca_coords(&s);
        let sup = kabsch(&ca, &ca).unwrap();
        assert!(sup.rmsd < 1.0e-6);
    }

    #[test]
    fn crambin_pdb_parses_as_a_real_protein() {
        // The embedded 1CRN is a real 46-residue / 327-atom protein. The
        // existing `read_structure` PDB parser must consume it as-is (HELIX /
        // SHEET / TER / the OXT terminal oxygen included).
        let s = read_structure(CRAMBIN_PDB, "crambin").expect("crambin (1CRN) must parse");
        let (mut res, mut atoms) = (0usize, 0usize);
        for chain in &s.first_model().chains {
            res += chain.residues.len();
            for r in &chain.residues {
                atoms += r.atoms.len();
            }
        }
        assert_eq!(res, 46, "crambin has 46 residues");
        assert_eq!(atoms, 327, "crambin has 327 atoms");
        // 46 Cα atoms (one per residue) — enough for a ribbon and a meaningful
        // RMSD, unlike the 3-Cα demo peptide.
        assert_eq!(ca_coords(&s).len(), 46);
    }

    // ----- in-house ENM molecular dynamics --------------------------------

    /// Helper: build crambin's ViewMolecule for the MD tests.
    fn crambin_view() -> ViewMolecule {
        let s = read_structure(CRAMBIN_PDB, "crambin").unwrap();
        ViewMolecule::from_biostruct(&s)
    }

    #[test]
    fn enm_md_runs_on_crambin_and_is_well_formed() {
        let mol = crambin_view();
        let md = enm_md(&mol, EnmParams::default());
        // One frame per requested frame, every frame the right atom count.
        assert_eq!(md.frames.len(), EnmParams::default().n_frames);
        for f in &md.frames {
            assert_eq!(f.len(), mol.atoms.len());
            for c in f {
                assert!(c.iter().all(|x| x.is_finite()), "no NaN/Inf coordinates");
            }
        }
        // Crambin is dense → the 8 Å network has many springs.
        assert!(md.n_springs > 1000, "ENM built {} springs", md.n_springs);
    }

    #[test]
    fn enm_md_frame0_is_the_reference_structure() {
        // Velocity-Verlet starts from the input coordinates, so frame 0 must be
        // the reference exactly (the viewer shows the input before playback).
        let mol = crambin_view();
        let md = enm_md(&mol, EnmParams::default());
        for (a, atom) in md.frames[0].iter().zip(&mol.atoms) {
            assert_eq!(*a, atom.pos);
        }
    }

    #[test]
    fn enm_md_is_deterministic_for_a_fixed_seed() {
        let mol = crambin_view();
        let a = enm_md(&mol, EnmParams::default());
        let b = enm_md(&mol, EnmParams::default());
        assert_eq!(a.frames, b.frames, "same seed ⇒ identical trajectory");
        // A different seed gives a *different* (still valid) trajectory.
        let c = enm_md(
            &mol,
            EnmParams {
                seed: 0xC0FFEE,
                ..EnmParams::default()
            },
        );
        assert_ne!(a.frames, c.frames, "different seed ⇒ different trajectory");
    }

    #[test]
    fn enm_md_energy_is_bounded_no_blowup() {
        // The symplectic integrator + light damping keep the total energy in a
        // bounded band — it must not drift up by orders of magnitude (the
        // signature of an unstable run). We integrate a long run and assert the
        // final energy stays close to the initial.
        let mol = crambin_view();
        let md = enm_md(
            &mol,
            EnmParams {
                n_frames: 200,
                substeps: 8,
                ..EnmParams::default()
            },
        );
        assert!(md.energy_first.is_finite() && md.energy_last.is_finite());
        assert!(md.energy_first > 0.0);
        // Final within a factor of ~2 of initial (damping pulls it *down*, never
        // an exponential blow-up).
        let ratio = md.energy_last / md.energy_first;
        assert!(
            (0.25..=2.0).contains(&ratio),
            "energy ratio {ratio} out of bounds (first {}, last {})",
            md.energy_first,
            md.energy_last
        );
    }

    #[test]
    fn enm_md_atoms_stay_near_equilibrium() {
        // A physical thermal vibration keeps every atom within a fraction of an
        // ångström of its start — the RMSD must stay small (no explosion).
        let mol = crambin_view();
        let md = enm_md(&mol, EnmParams::default());
        assert!(md.rmsd_max > 0.0, "the structure must actually move");
        assert!(
            md.rmsd_max < 1.0,
            "peak RMSD {:.3} Å too large — not a gentle thermal wobble",
            md.rmsd_max
        );
    }

    #[test]
    fn enm_md_higher_temperature_gives_larger_motion() {
        // Monotone physics check: a hotter run vibrates more (larger RMSD).
        let mol = crambin_view();
        let cold = enm_md(
            &mol,
            EnmParams {
                temperature: 0.1,
                ..EnmParams::default()
            },
        );
        let hot = enm_md(
            &mol,
            EnmParams {
                temperature: 0.8,
                ..EnmParams::default()
            },
        );
        assert!(
            hot.rmsd_max > cold.rmsd_max,
            "hot RMSD {:.3} should exceed cold {:.3}",
            hot.rmsd_max,
            cold.rmsd_max
        );
    }

    #[test]
    fn enm_md_short_molecule_is_a_no_op() {
        // Fewer than two atoms → nothing to vibrate, empty result, no panic.
        let one = ViewMolecule {
            atoms: vec![ViewAtom::new([0.0, 0.0, 0.0], "C")],
            bonds: vec![],
        };
        let md = enm_md(&one, EnmParams::default());
        assert!(md.frames.is_empty());
    }

    #[test]
    fn run_md_attaches_a_playable_trajectory_and_readout() {
        // The panel-level path used by the button + agent bridge: run MD on the
        // default (crambin) structure and confirm it attaches and reports.
        let mut p = BiostructPanel::default();
        let summary = p.run_md().expect("MD must attach on the default crambin");
        assert!(p.trajectory.is_attached(), "trajectory must be attached");
        assert!(p.trajectory.n_frames() >= 2);
        assert!(summary.contains("ENM MD"));
        // The readout reflects the run, not the "no run yet" hint.
        let readout = p.md_readout();
        assert!(readout.contains("ENM MD") && readout.contains("springs"));
    }

    #[test]
    fn enm_md_temperature_zero_is_static_but_stable() {
        // T = 0 ⇒ no thermal velocities ⇒ the structure sits at its minimum and
        // does not move; still finite, still bounded (a degenerate-but-valid run).
        let mol = crambin_view();
        let md = enm_md(
            &mol,
            EnmParams {
                temperature: 0.0,
                ..EnmParams::default()
            },
        );
        assert!(md.rmsd_max < 1e-6, "T=0 keeps the structure put");
        assert!(md.energy_last.is_finite());
    }

    #[test]
    fn crambin_backbone_has_both_helix_and_sheet() {
        // The whole point of crambin as the default: its DSSP-tagged Cα trace
        // shows *both* an α-helix and a β-sheet, so the ribbon / cartoon and
        // the secondary-structure colour scheme have real structure to render.
        use crate::molviz::SsKind;
        let s = read_structure(CRAMBIN_PDB, "crambin").unwrap();
        let bb = super::ca_backbone(&s);
        assert_eq!(bb.len(), 46, "one Cα control point per residue");
        let kinds: Vec<Option<SsKind>> = bb
            .iter()
            .map(|p| p.ss.map(SsKind::from_dssp_code))
            .collect();
        assert!(
            kinds.contains(&Some(SsKind::Helix)),
            "crambin must contain α-helix backbone points"
        );
        assert!(
            kinds.contains(&Some(SsKind::Sheet)),
            "crambin must contain β-sheet backbone points"
        );
    }

    #[test]
    fn default_panel_loads_crambin_not_the_tiny_peptide() {
        // The viewer's default structure is the real protein, so opening the
        // panel shows something — not the 12-atom smoke peptide.
        let p = BiostructPanel::default();
        assert_eq!(p.demo, DemoStructure::Crambin);
        assert_eq!(p.structure_a, CRAMBIN_PDB);
        let s = read_structure(&p.structure_a, "default").unwrap();
        assert_eq!(ca_coords(&s).len(), 46);
    }

    #[test]
    fn demo_structure_from_label_round_trips_and_rejects_garbage() {
        for d in DemoStructure::ALL {
            assert_eq!(DemoStructure::from_label(d.label()), Ok(d));
            assert_eq!(DemoStructure::from_label(&d.label().to_uppercase()), Ok(d));
        }
        // Short aliases resolve too.
        assert_eq!(
            DemoStructure::from_label("crambin"),
            Ok(DemoStructure::Crambin)
        );
        assert_eq!(
            DemoStructure::from_label("1CRN"),
            Ok(DemoStructure::Crambin)
        );
        assert_eq!(
            DemoStructure::from_label("triglycine"),
            Ok(DemoStructure::GlyPeptide)
        );
        assert!(DemoStructure::from_label("hemoglobin").is_err());
    }

    #[test]
    fn load_demo_swaps_the_structure_and_is_undoable() {
        let mut p = BiostructPanel::default(); // crambin
        p.load_demo(DemoStructure::GlyPeptide);
        assert_eq!(p.demo, DemoStructure::GlyPeptide);
        assert_eq!(p.structure_a, DEMO_PDB);
        assert!(p.can_undo(), "loading a demo records an undo step");
        // Undo restores the prior (crambin) structure text.
        assert!(p.undo_edit());
        assert_eq!(p.structure_a, CRAMBIN_PDB);
    }

    #[test]
    fn structure_readout_names_the_loaded_structure_with_counts() {
        let mut p = BiostructPanel::default();
        let r = p.structure_readout();
        assert!(r.contains("Crambin"), "readout names crambin: {r}");
        assert!(r.contains("46 residues"), "readout has residue count: {r}");
        assert!(r.contains("327 atoms"), "readout has atom count: {r}");
        p.load_demo(DemoStructure::GlyPeptide);
        let r = p.structure_readout();
        assert!(r.contains("Triglycine"), "readout names the peptide: {r}");
        assert!(r.contains("3 residues"), "readout has residue count: {r}");
    }

    /// Round-21 H1 RED→GREEN: the genetics-workbench file loaders
    /// route through [`valenx_core::io_caps::read_capped_to_string`]
    /// with the [`valenx_core::io_caps::MAX_GENETICS_FILE_BYTES`]
    /// cap. Pre-fix the loader did a bare `fs::read_to_string`, so a
    /// user-picked multi-GB file would OOM. We exercise the helper
    /// here against an oversized scratch file (allocating 100 MiB on
    /// disk would slow CI, so the test uses a small cap so a
    /// modestly-sized file trips it the same way 64 MiB+ would trip
    /// the production cap).
    #[test]
    fn oversize_genetics_load_returns_invalid_data() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join("valenx_r21_biostruct_oversize.pdb");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&vec![b'A'; 4096]).unwrap();
        drop(f);
        // Cap of 1 KiB simulates the 64 MiB production cap shape —
        // a file larger than the cap is rejected with InvalidData.
        let err = valenx_core::io_caps::read_capped_to_string(&tmp, 1024).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let _ = std::fs::remove_file(&tmp);
        // Sanity: the constant the loader actually uses is the
        // production cap (proves we didn't downgrade by accident).
        assert_eq!(
            valenx_core::io_caps::MAX_GENETICS_FILE_BYTES,
            64u64 * 1024 * 1024
        );
    }
}

/// Headless egui UI-logic tests for the Macromolecular Structure panel.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::genetics_workbench::GeneticsPanel;
    use crate::ValenxApp;

    fn draw_headless(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
    }

    /// Draw the panel once with the accesskit tree enabled and return its
    /// nodes — the harness the AI-drivability (`labelled_by`) assertions use.
    fn draw_and_collect_nodes(
        app: &mut ValenxApp,
    ) -> Vec<(egui::accesskit::NodeId, egui::accesskit::Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::MacromolecularStructure;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [Tool::Analyze, Tool::Ramachandran, Tool::Superpose] {
            let mut app = app_with_panel();
            app.genetics.biostruct.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.biostruct.result = "chains : 1\nresidues : 3\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.biostruct.error = Some("could not parse PDB".to_string());
        draw_headless(&mut app);
        // Empty structure text — the "Show in 3D" affordance is gated.
        let mut app = app_with_panel();
        app.genetics.biostruct.structure_a.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn run_single_analyzes_the_demo_structure() {
        // The demo glycine peptide → the real valenx-biostruct
        // analyzer produces a correctly-formatted report.
        let mut p = BiostructPanel {
            tool: Tool::Analyze,
            ..BiostructPanel::default()
        };
        run_single(&mut p);
        assert!(p.error.is_none(), "analyze errored: {:?}", p.error);
        assert!(p.result.contains("chains"));
        assert!(p.result.contains("residues"));
    }

    #[test]
    fn run_single_classifies_ramachandran() {
        let mut p = BiostructPanel {
            tool: Tool::Ramachandran,
            ..BiostructPanel::default()
        };
        run_single(&mut p);
        assert!(p.error.is_none(), "Ramachandran errored: {:?}", p.error);
        assert!(!p.result.is_empty());
    }

    #[test]
    fn run_superpose_aligns_identical_structures() {
        // The demo structure superposed on itself → near-zero RMSD.
        let mut p = BiostructPanel {
            tool: Tool::Superpose,
            ..BiostructPanel::default()
        };
        run_superpose(&mut p);
        assert!(p.error.is_none(), "superpose errored: {:?}", p.error);
        assert!(p.result.contains("paired Cα atoms"));
        assert!(p.result.contains("RMSD after fit"));
    }

    #[test]
    fn run_actions_surface_errors_on_bad_input() {
        // A non-PDB string is malformed structure input.
        let mut p = BiostructPanel {
            tool: Tool::Analyze,
            structure_a: "not a structure".to_string(),
            ..BiostructPanel::default()
        };
        run_single(&mut p);
        assert!(p.error.is_some(), "analyze should error on malformed input");
        // Superpose with too few Cα atoms (an empty mobile structure).
        let mut p = BiostructPanel {
            tool: Tool::Superpose,
            structure_a: "END\n".to_string(),
            ..BiostructPanel::default()
        };
        run_superpose(&mut p);
        assert!(p.error.is_some(), "superpose should error with no Cα atoms");
    }

    #[test]
    fn draws_with_every_representation_selected_without_panic() {
        // The representation picker + its per-mode tuning row must draw for
        // every representation (the demo PDB is loaded by default, so the
        // picker is shown).
        for rep in Representation::ALL {
            let mut app = app_with_panel();
            app.genetics.biostruct.representation = rep;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn show_in_viewport_meshes_the_demo_for_each_representation() {
        // The demo glycine peptide is a real protein backbone, so every
        // representation — including cartoon (needs Cα) and surface (marching
        // cubes) — produces a non-empty mesh in the viewport without error.
        for rep in Representation::ALL {
            let mut app = app_with_panel();
            app.genetics.biostruct.representation = rep;
            // Keep the surface grid small so the test stays fast.
            app.genetics.biostruct.molviz_params.grid_max = 24;
            super::show_in_viewport(&mut app);
            assert!(
                app.genetics.biostruct.error.is_none(),
                "{rep:?} errored: {:?}",
                app.genetics.biostruct.error
            );
            assert!(
                app.stl.is_some(),
                "{rep:?} should have pushed a mesh to the viewport"
            );
        }
    }

    #[test]
    fn cartoon_errors_on_a_structure_with_no_backbone() {
        // A hetero-only structure (a lone zinc ion) has no Cα trace, so the
        // cartoon representation must fail loudly rather than show nothing.
        let mut app = app_with_panel();
        app.genetics.biostruct.representation = Representation::Cartoon;
        app.genetics.biostruct.structure_a = "\
HETATM    1 ZN    ZN A   1       0.000   0.000   0.000  1.00  0.00          ZN
END
"
        .to_string();
        super::show_in_viewport(&mut app);
        assert!(
            app.genetics.biostruct.error.is_some(),
            "cartoon with no backbone should surface an error"
        );
    }

    #[test]
    fn ca_backbone_extracts_dssp_tagged_trace() {
        // The demo 3-glycine peptide yields 3 Cα control points, each carrying
        // a DSSP secondary-structure code (so the cartoon tube can modulate).
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let bb = super::ca_backbone(&s);
        assert_eq!(bb.len(), 3, "three Cα control points");
        assert!(
            bb.iter().all(|p| p.ss.is_some()),
            "every backbone point carries a DSSP code"
        );
    }

    #[test]
    fn draws_with_every_color_scheme_selected_without_panic() {
        // The colour-scheme picker (+ its per-scheme note row) must draw for
        // every scheme, in combination with each representation (the demo PDB
        // is loaded by default, so the picker is shown).
        for scheme in ColorScheme::ALL {
            for rep in Representation::ALL {
                let mut app = app_with_panel();
                app.genetics.biostruct.color_scheme = scheme;
                app.genetics.biostruct.representation = rep;
                draw_headless(&mut app);
            }
        }
    }

    #[test]
    fn structure_atom_attrs_lockstep_with_view_atoms() {
        // The per-atom attribute slice must be in lockstep with the
        // `ViewMolecule::from_biostruct` atom order so the colour schemes line
        // up atom-for-atom. The demo 3-glycine peptide has 12 atoms across one
        // chain "A" and three residues.
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let view = ViewMolecule::from_biostruct(&s);
        let attrs = super::structure_atom_attrs(&s);
        assert_eq!(
            attrs.len(),
            view.atoms.len(),
            "one AtomAttr per ViewMolecule atom, in order"
        );
        assert!(attrs.iter().all(|a| a.chain == "A"), "demo is all chain A");
        // Residue index is a monotone 0-based counter across the model: three
        // residues → indices 0, 1, 2 present.
        let max_res = attrs.iter().map(|a| a.residue_index).max().unwrap();
        assert_eq!(max_res, 2, "three residues → max residue index 2");
    }

    #[test]
    fn structure_atom_attrs_carry_secondary_structure_in_lockstep() {
        // Each atom's per-atom SS must equal its residue's DSSP state — the
        // same per-chain assignment `ca_backbone` tags the Cα trace with —
        // keyed by residue position, in the exact `ViewMolecule` atom order.
        use crate::molviz::SsKind;
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let view = ViewMolecule::from_biostruct(&s);
        let attrs = super::structure_atom_attrs(&s);
        assert_eq!(attrs.len(), view.atoms.len(), "one attr per atom, in order");
        // Every atom of the single-chain demo carries an SS state (the demo is
        // all amino acids, so DSSP assigns every residue — none is `None`).
        assert!(
            attrs.iter().all(|a| a.ss.is_some()),
            "every demo atom inherits its residue's DSSP state"
        );

        // Cross-check against the independent `ca_backbone` assignment: the SS
        // of each residue's Cα (in backbone order) must match the SS the atom
        // attrs carry for that residue's atoms. Build the per-residue expected
        // SS from `ca_backbone` and confirm the atom attrs agree residue-for-
        // residue.
        let backbone = super::ca_backbone(&s);
        let expected_per_res: Vec<Option<SsKind>> = backbone
            .iter()
            .map(|bp| bp.ss.map(SsKind::from_dssp_code))
            .collect();
        // Walk the structure in the same chain→residue→atom order and compare.
        // `residue_index` is the monotone all-residue counter (the attrs' key);
        // `ca_pos` is the backbone (Cα-only) position used to index the expected
        // SS — kept separate so the comparison is correct even if some residue
        // lacks a Cα. For the demo every residue has a Cα, so they coincide.
        let mut residue_index = 0usize;
        let mut ca_pos = 0usize;
        for chain in &s.first_model().chains {
            for res in &chain.residues {
                if res.ca().is_some() {
                    let want = expected_per_res[ca_pos];
                    for atom_attr in attrs
                        .iter()
                        .filter(|a| a.residue_index == residue_index as i32)
                    {
                        assert_eq!(
                            atom_attr.ss, want,
                            "atom SS must match its residue's backbone DSSP state"
                        );
                    }
                    ca_pos += 1;
                }
                residue_index += 1;
            }
        }
    }

    #[test]
    fn ss_color_scheme_reads_the_plumbed_per_atom_ss() {
        // End-to-end: the SS colour scheme paints each atom by the per-atom SS
        // that `structure_atom_attrs` plumbed in — helix atoms helix-red, sheet
        // atoms sheet-yellow, coil atoms grey, a missing SS grey (no panic).
        // We assert the scheme's colour for every atom equals the canonical SS
        // colour for that atom's plumbed state (so the chain of evidence —
        // residue DSSP → per-atom ss → rendered colour — is unbroken).
        use crate::molviz::{atom_color, ColorContext, ColorScheme, SsKind};
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let attrs = super::structure_atom_attrs(&s);
        let ctx = ColorContext::build(&attrs);
        // The standard convention colours (mirrors molviz::ss_color).
        let conv = |ss: Option<SsKind>| -> [f32; 3] {
            match ss {
                Some(SsKind::Helix) => [0.90, 0.18, 0.55],
                Some(SsKind::Sheet) => [0.95, 0.85, 0.18],
                Some(SsKind::Coil) | None => [0.80, 0.80, 0.80],
            }
        };
        for a in &attrs {
            let got = atom_color(ColorScheme::SecondaryStructure, "C", a, &ctx);
            assert_eq!(
                got,
                conv(a.ss),
                "SS scheme must colour the atom by its plumbed SS state"
            );
            assert!(got
                .iter()
                .all(|&x| x.is_finite() && (0.0..=1.0).contains(&x)));
        }
        // An atom whose SS is None colours as coil grey, never a panic.
        let none = crate::molviz::AtomAttr::new("A", 0, 0.0);
        assert_eq!(
            atom_color(ColorScheme::SecondaryStructure, "C", &none, &ctx),
            [0.80, 0.80, 0.80],
            "a missing SS must render as coil grey"
        );
    }

    #[test]
    fn ss_color_scheme_is_named_in_the_picker_for_ai_driving() {
        // The "Secondary structure" colour-scheme option must be an AI-drivable
        // node: it appears by its label as an accesskit Name and is labelled_by
        // its "Colour scheme" caption (so an agent / screen reader can select
        // it by name through the accessibility tree).
        use egui::accesskit::Node;
        let mut app = app_with_panel();
        let nodes = draw_and_collect_nodes(&mut app);
        let label = ColorScheme::SecondaryStructure.label(); // "Secondary structure"
        let matching: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.name() == Some(label))
            .collect();
        assert!(
            !matching.is_empty(),
            "the '{label}' option must appear as a named accesskit node"
        );
        assert!(
            matching.iter().any(|n| !n.labelled_by().is_empty()),
            "the '{label}' option must be labelled_by its caption (AI-drivable)"
        );
    }

    #[test]
    fn show_in_viewport_uploads_per_vertex_colors_for_each_non_default_scheme() {
        // A non-default colour scheme builds a colour-aware mesh and uploads
        // per-vertex colours: the viewport STL carries a `colors` buffer whose
        // length is exactly 3 × the triangle count (one colour per surface
        // vertex), for every representation. The default Element scheme keeps
        // the monochrome path (no colour buffer).
        for scheme in ColorScheme::ALL {
            for rep in Representation::ALL {
                let mut app = app_with_panel();
                app.genetics.biostruct.color_scheme = scheme;
                app.genetics.biostruct.representation = rep;
                // Keep the surface/density grids small so the test stays fast.
                app.genetics.biostruct.molviz_params.grid_max = 24;
                app.genetics.biostruct.molviz_params.density_grid_max = 24;
                super::show_in_viewport(&mut app);
                assert!(
                    app.genetics.biostruct.error.is_none(),
                    "{scheme:?}/{rep:?} errored: {:?}",
                    app.genetics.biostruct.error
                );
                let stl = app
                    .stl
                    .as_ref()
                    .unwrap_or_else(|| panic!("{scheme:?}/{rep:?} pushed no mesh"));
                let tri_count = stl.mesh.triangle_count();
                match scheme {
                    ColorScheme::Element => {
                        // Default scheme → no per-vertex colour override.
                        assert!(
                            stl.colors.is_none(),
                            "{rep:?}: Element scheme must keep the monochrome path"
                        );
                    }
                    _ => {
                        let colors = stl.colors.as_ref().unwrap_or_else(|| {
                            panic!("{scheme:?}/{rep:?} must carry per-vertex colours")
                        });
                        assert_eq!(
                            colors.len(),
                            tri_count * 3,
                            "{scheme:?}/{rep:?}: colours must equal 3 × triangle count"
                        );
                        // Every colour component is a finite 0..=1 value.
                        assert!(
                            colors.iter().all(|c| c
                                .iter()
                                .all(|&x| x.is_finite() && (0.0..=1.0).contains(&x))),
                            "{scheme:?}/{rep:?}: colour components must be finite in [0,1]"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn demo_structure_picker_is_labelled_for_ai_driving() {
        // The "Demo structure" picker must be AI-drivable: each demo option
        // carries its label as an accessibility Name, and the combo is
        // `labelled_by` its "Demo structure" caption — so an agent can switch
        // the loaded structure (e.g. to crambin) through the accessibility tree.
        use egui::accesskit::Node;
        let mut app = app_with_panel();
        let nodes = draw_and_collect_nodes(&mut app);

        // The caption itself appears as a named node.
        assert!(
            nodes
                .iter()
                .any(|(_, n)| n.name() == Some("Demo structure")),
            "the 'Demo structure' caption must be a named accesskit node"
        );
        // Every demo option appears as a named node, and at least one carries a
        // `labelled_by` association (the combo body links to its caption).
        for d in DemoStructure::ALL {
            let label = d.label();
            let matching: Vec<&Node> = nodes
                .iter()
                .map(|(_, n)| n)
                .filter(|n| n.name() == Some(label))
                .collect();
            assert!(
                !matching.is_empty(),
                "demo structure '{label}' must appear as a named node in the accesskit tree"
            );
        }
        assert!(
            nodes.iter().any(|(_, n)| !n.labelled_by().is_empty()),
            "the demo-structure combo must expose a labelled_by association (AI-drivable)"
        );
    }

    #[test]
    fn color_scheme_picker_is_labelled_for_ai_driving() {
        // The colour-scheme picker must be AI-drivable: each scheme button
        // carries the scheme label as its accessibility Name AND is
        // `labelled_by` the "Colour scheme" caption (so an agent / screen
        // reader can find and select it by name). Assert the four scheme labels
        // appear as named, `labelled_by`-associated nodes in the accesskit tree.
        use egui::accesskit::Node;
        let mut app = app_with_panel();
        let nodes = draw_and_collect_nodes(&mut app);

        for scheme in ColorScheme::ALL {
            let label = scheme.label();
            let matching: Vec<&Node> = nodes
                .iter()
                .map(|(_, n)| n)
                .filter(|n| n.name() == Some(label))
                .collect();
            assert!(
                !matching.is_empty(),
                "colour scheme '{label}' must appear as a named node in the accesskit tree"
            );
            assert!(
                matching.iter().any(|n| !n.labelled_by().is_empty()),
                "colour scheme '{label}' button must be labelled_by its caption (AI-drivable)"
            );
        }
    }

    // ---- MD-trajectory playback ----------------------------------------

    /// Build a 3-frame `valenx-md` trajectory whose atom count matches the demo
    /// glycine peptide (12 atoms), so it attaches cleanly. Frame `f` shifts
    /// every atom by `f` Å on x, giving frames that are trivially distinguishable.
    fn matching_traj(natoms: usize) -> valenx_md::io::trajectory::Trajectory {
        use nalgebra::Vector3;
        let mut traj = valenx_md::io::trajectory::Trajectory::new(natoms, 0.002).unwrap();
        for f in 0..3 {
            let frame: Vec<Vector3<f64>> = (0..natoms)
                // nm here; ×10 → ångström on attach. `f` nm → 10·f Å shift.
                .map(|i| Vector3::new(f as f64 + 0.1 * i as f64, 0.0, 0.0))
                .collect();
            traj.push_frame(frame).unwrap();
        }
        traj
    }

    #[test]
    fn synthetic_trajectory_attaches_and_frame0_is_the_reference() {
        // Generating the demo trajectory attaches a multi-frame wobble whose
        // frame 0 is the reference structure exactly (sin 0 == 0).
        let mut app = app_with_panel();
        super::attach_synthetic_trajectory(&mut app);
        let tj = &app.genetics.biostruct.trajectory;
        assert!(tj.is_attached(), "synthetic trajectory must attach");
        assert!(tj.n_frames() >= 2, "needs ≥ 2 frames to animate");
        assert!(tj.note.is_none(), "clean attach has no note: {:?}", tj.note);
        // Frame 0 reproduces the base atom positions exactly.
        let base = tj.base.as_ref().unwrap();
        let f0 = tj.molecule_at(0).unwrap();
        assert_eq!(f0.atoms.len(), base.atoms.len());
        for (a, b) in f0.atoms.iter().zip(&base.atoms) {
            assert_eq!(a.pos, b.pos, "frame 0 must equal the reference structure");
        }
    }

    #[test]
    fn stepping_the_frame_index_changes_atom_positions_exactly() {
        // The core playback contract: stepping the frame index sets the
        // displayed atom positions to that frame's coordinates, exactly.
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let natoms = ViewMolecule::from_biostruct(&s).atoms.len();
        let traj = matching_traj(natoms);
        let mut tj = TrajectoryPlayback::default();
        tj.attach_md(&s, &traj);
        assert!(tj.is_attached());
        assert_eq!(tj.n_frames(), 3);

        // Frame f shifts atom i to x = 10·(f + 0.1·i) Å (nm→Å on attach).
        for f in 0..3 {
            let mol = tj.molecule_at(f).expect("frame in range");
            assert_eq!(mol.atoms.len(), natoms);
            for (i, atom) in mol.atoms.iter().enumerate() {
                let expect_x = 10.0 * (f as f32 + 0.1 * i as f32);
                assert!(
                    (atom.pos[0] - expect_x).abs() < 1e-3,
                    "frame {f} atom {i}: x {} != {expect_x}",
                    atom.pos[0]
                );
                assert_eq!(atom.pos[1], 0.0);
                assert_eq!(atom.pos[2], 0.0);
            }
        }
        // Consecutive frames differ (the motion is real, not static).
        let a = tj.molecule_at(0).unwrap();
        let b = tj.molecule_at(1).unwrap();
        assert_ne!(a.atoms[0].pos, b.atoms[0].pos);
        // Out-of-range frame is None, never a panic.
        assert!(tj.molecule_at(99).is_none());
    }

    #[test]
    fn atom_count_mismatch_is_handled_without_panic() {
        // A trajectory whose per-frame atom count differs from the structure
        // must NOT attach and must leave an in-panel note (no panic).
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let natoms = ViewMolecule::from_biostruct(&s).atoms.len();
        let wrong = matching_traj(natoms + 5); // deliberately too many atoms
        let mut tj = TrajectoryPlayback::default();
        tj.attach_md(&s, &wrong);
        assert!(!tj.is_attached(), "a mismatched trajectory must not attach");
        assert!(
            tj.note
                .as_deref()
                .unwrap_or("")
                .contains("atom counts must match"),
            "mismatch must leave an explanatory note: {:?}",
            tj.note
        );

        // The same guard via the explicit-frames path: a short frame is rejected.
        let base = ViewMolecule::from_biostruct(&s);
        let mut tj2 = TrajectoryPlayback::default();
        tj2.attach(
            base,
            Vec::new(),
            Vec::new(),
            vec![vec![[0.0; 3]; natoms], vec![[0.0; 3]; natoms - 1]],
            "test",
        );
        assert!(!tj2.is_attached());
        assert!(tj2.note.is_some());
    }

    #[test]
    fn empty_trajectory_is_a_noop() {
        // An empty `valenx-md` trajectory attaches nothing (no-op), no panic.
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let natoms = ViewMolecule::from_biostruct(&s).atoms.len();
        let empty = valenx_md::io::trajectory::Trajectory::new(natoms, 0.002).unwrap();
        let mut tj = TrajectoryPlayback::default();
        tj.attach_md(&s, &empty);
        assert!(!tj.is_attached(), "empty trajectory must not attach");
    }

    #[test]
    fn render_trajectory_frame_pushes_the_frames_geometry_to_the_viewport() {
        // The per-frame update meshes the current frame and pushes it into the
        // viewport; switching frames re-renders without panic.
        let mut app = app_with_panel();
        super::attach_synthetic_trajectory(&mut app);
        assert!(app.genetics.biostruct.trajectory.is_attached());
        // Render a couple of distinct frames — each must populate the viewport.
        for f in [0usize, 1, 30] {
            app.genetics.biostruct.trajectory.frame =
                f.min(app.genetics.biostruct.trajectory.n_frames() - 1);
            super::render_trajectory_frame(&mut app);
            assert!(
                app.stl.is_some(),
                "frame {f} should have pushed a mesh to the viewport"
            );
            assert!(
                app.genetics.biostruct.trajectory.note.is_none(),
                "frame {f} render errored: {:?}",
                app.genetics.biostruct.trajectory.note
            );
        }
    }

    #[test]
    fn draws_at_a_mid_frame_without_panic() {
        // The full panel (with the playback transport) must draw at a mid-frame,
        // both paused and playing, without panic.
        let mut app = app_with_panel();
        super::attach_synthetic_trajectory(&mut app);
        let mid = app.genetics.biostruct.trajectory.n_frames() / 2;
        app.genetics.biostruct.trajectory.frame = mid;
        draw_headless(&mut app); // paused
        app.genetics.biostruct.trajectory.playing = true;
        draw_headless(&mut app); // playing (advances the clock + repaints)
    }

    #[test]
    fn playback_controls_are_labelled_for_ai_driving() {
        // The frame slider, Play/Pause and speed controls must be AI-drivable:
        // each carries its caption as an accessibility Name and is `labelled_by`
        // that caption. Attach a trajectory first so the transport is shown.
        use egui::accesskit::Node;
        let mut app = app_with_panel();
        super::attach_synthetic_trajectory(&mut app);
        let nodes = draw_and_collect_nodes(&mut app);

        // The captions that name the playback groups.
        for caption in ["Frame", "Transport", "Speed (fps)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "playback caption '{caption}' must appear as a named accesskit node"
            );
        }
        // At least some control is `labelled_by` a caption (the slider / drag /
        // buttons are associated to their captions) — the AI-drivable handle.
        let labelled: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| !n.labelled_by().is_empty())
            .collect();
        assert!(
            !labelled.is_empty(),
            "playback transport must expose labelled_by-associated controls"
        );
        // The Play button is reachable by name.
        assert!(
            nodes
                .iter()
                .any(|(_, n)| n.name() == Some("Play") || n.name() == Some("Pause")),
            "the Play/Pause button must be a named accesskit node"
        );
    }
}
