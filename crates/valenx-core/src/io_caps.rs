//! Size caps for file I/O.
//!
//! Round-12 hardening: every `fs::read_to_string` site in the
//! workbench `persist.rs` family (sketch, cam, arch, feature-tree,
//! techdraw, surface, spreadsheet, macro, draft, lattice, assembly)
//! and the app-level mesh-load + sweep-reload sites went through a
//! bare unbounded read that would slurp an arbitrarily large file
//! into memory before any parser saw it. The helpers here enforce a
//! stat-and-bounded-take pattern (TOCTOU defence-in-depth — the stat
//! might lie if the file grew between metadata and open) so a
//! multi-GB payload is rejected before allocation.
//!
//! Mirrors the per-file cap in
//! `valenx_core::project::loader::read_capped` (which uses its own
//! 1 MiB project-scoped cap) and the
//! `valenx_app::rbac_io::rbac_override_from_project_toml` cap.

use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
#[cfg(unix)]
use std::os::unix::fs::OpenOptionsExt as _;
#[cfg(windows)]
use std::os::windows::fs::OpenOptionsExt as _;

/// Round-27 STRUCTURAL: process-monotonic counter appended to the
/// sidecar `<basename>.tmp.<pid>.<counter>` name so two concurrent
/// writers in the SAME PID can't collide on the per-process clock.
/// `(pid, counter)` is strictly unique across the process lifetime,
/// guaranteeing each writer owns a distinct sidecar path even when
/// many threads race the same nanosecond.
///
/// Replaces the (pid, nanos, counter) pattern from R25 M3 (which
/// used three tuple components for distinctness) — the counter
/// alone, combined with the PID, is already unique-across-process,
/// so the nanosecond component is redundant for collision-safety
/// purposes and only added cosmetic distinctness.
static ATOMIC_WRITE_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Windows `FILE_FLAG_OPEN_REPARSE_POINT` — opens the reparse point
/// itself rather than following it. Hardcoded so this crate doesn't
/// need a `windows-sys` dependency for a single constant. Mirrors
/// the dock-runner pattern (R25 H2 / R26 M2).
#[cfg(windows)]
const FILE_FLAG_OPEN_REPARSE_POINT: u32 = 0x0020_0000;

/// Default cap for in-workbench document files (`.ron` for
/// sketch / arch / cam / techdraw / lattice / etc.). 16 MiB is well
/// above any realistic file in those domains — even a heavily
/// populated assembly is at most a few hundred KB of RON — while
/// still being small enough to allocate without OOMing a 32-bit
/// process.
pub const MAX_DOC_FILE_BYTES: usize = 16 * 1024 * 1024;

/// Cap for the `.json` mesh-load path in `valenx-app`. Mesh files
/// (especially debug exports of subdivision results) can be larger
/// than ordinary documents, but 64 MiB of JSON is already absurd — a
/// 64 MiB JSON mesh is millions of vertices, well past anything an
/// interactive viewer can handle.
pub const MAX_MESH_JSON_BYTES: usize = 64 * 1024 * 1024;

/// Cap for the per-derived-case `results.json` files the sweep
/// aggregator reads back in `valenx_app::sweep::assemble_sweep_dataset`.
///
/// Round-14 H3 (R13 carry-over): pre-fix the sweep aggregator read
/// each subdir's `results.json` with a bare `fs::read_to_string`, so
/// a poisoned multi-GB file written under a derived workdir could be
/// slurped into memory before serde_json saw it. 64 MiB is the same
/// shape as [`MAX_MESH_JSON_BYTES`] — a sweep result with that much
/// JSON is unrealistic (results.json is typically a few KB; even
/// adapters that dump every solver residual top out in the low
/// hundreds of KB).
pub const MAX_RESULTS_JSON_BYTES: usize = 64 * 1024 * 1024;

/// Cap for `.pvd` time-series manifest files. Real-world `.pvd`
/// manifests reference a few hundred to a few thousand `.vtu`
/// snapshots — even a transient run with millions of time steps
/// rarely tops a few hundred KiB. 16 MiB is generous while
/// rejecting a hostile multi-GB manifest that would slurp into
/// memory before the XML parser saw any of it.
///
/// Round-18 M2: the auto-loader's `latest_entry_from_pvd` was
/// reading the manifest with a bare `fs::read_to_string`.
pub const MAX_PVD_FILE_BYTES: usize = 16 * 1024 * 1024;

/// Cap for CalculiX `.frd` ASCII result files. Real-world `.frd`
/// dumps for a 1M-DOF model with several time steps + stress
/// components top out around 100 MiB; 256 MiB is generous while
/// rejecting a hostile multi-GB file that would slurp into memory
/// before the parser saw the first line.
///
/// Round-18 M1: both `valenx_fem::parse_frd` and the CalculiX
/// adapter's `load_frd_fields_into_results` previously did a bare
/// `fs::read_to_string` on artifact paths — a corrupted (or
/// adversarial) workdir with a multi-GB `.frd` would OOM the
/// renderer before any post-processing happened.
pub const MAX_FRD_FILE_BYTES: usize = 256 * 1024 * 1024;

/// Cap on the number of sibling subdirs the sweep aggregator will
/// walk under a single sweep parent directory.
///
/// Round-14 H3 (R13 carry-over): pre-fix
/// `assemble_sweep_dataset` accepted whatever `read_dir(parent)`
/// yielded, then sorted + walked the full list. A poisoned (or just
/// pathological) sweep parent dir with millions of subdirs would
/// allocate the full Vec before any reasonability check fired. 100k
/// is well past any realistic sweep (a 10-parameter grid sweep at the
/// grid-cell cap is 10M cells, but those are split across many
/// sweep-parent dirs; a single parent rarely tops a few thousand
/// derived cases).
pub const MAX_SWEEP_SIBLINGS: usize = 100_000;

/// Cap on the bytes a `.pdbqt` receptor / ligand / dock-output file
/// may consume when read by the Desktop Dock panel.
///
/// Round-20 H1: the Dock panel's `run_dock_now` + `push_pose_as_mesh`
/// previously did three bare `std::fs::read_to_string` against
/// user-chosen paths (receptor, ligand, dock output). A user (or a
/// stale path in the panel's saved state) pointing at a multi-GB file
/// would OOM the renderer process before the docker ever ran. 64 MiB
/// matches the MCP-side cap (see `valenx_mcp::tools::MAX_PDBQT_FILE_BYTES`)
/// and is generous for chemistry: even a PDBQT for a multi-chain
/// receptor with explicit waters is in the low-MB range.
pub const MAX_PDBQT_FILE_BYTES: usize = 64 * 1024 * 1024;

/// Cap on the bytes a Radiance `.hdr` environment-map file may
/// consume when loaded by the render-bridge environment loader.
///
/// Round-20 M2: the render-bridge `EnvironmentRef::load` did a bare
/// `std::fs::read(&self.hdr_path)` against a serialised path that
/// could be anywhere on disk. HDR maps can legitimately be larger
/// than ordinary documents (an 8K equirectangular HDR is ~256 MiB of
/// float32 pixels, but the Radiance `.hdr` RGBE encoding is ~32 MiB
/// — well under the cap), but a hostile or stale path could still
/// point at a multi-GB file. 256 MiB is generous for legitimate
/// 8K-resolution HDR environment maps while refusing the
/// `cat /dev/zero > big.hdr` denial of service.
///
/// Typed `u64` to match the [`read_capped_to_bytes`] cap parameter
/// (file metadata `len()` is `u64`).
pub const MAX_HDR_FILE_BYTES: u64 = 256 * 1024 * 1024;

/// Cap on the bytes a `.vtu` (XML) or `.vtk` (legacy binary) result
/// file may consume when read by the CFD adapters (SU2, OpenFOAM)
/// from a solver workdir.
///
/// Round-20 L1: pre-fix the SU2 + OpenFOAM `collect` paths walked
/// the workdir for VTK output and did a bare `std::fs::read(&path)`
/// on each match. A corrupted (or adversarial) workdir with a
/// multi-GB VTK file would slurp into memory before `vtk_dispatch`
/// saw the magic bytes. 4 GiB is generous — production CFD VTU
/// snapshots top out around 1 GiB for an HPC mesh with 100M cells
/// and a few fields per snapshot — while still refusing a hostile
/// or runaway-write artefact that would OOM the renderer.
pub const MAX_VTK_FILE_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// Cap on the bytes the `valenx_params.json` file an adapter stages
/// into a workdir may consume when read back during `collect()`.
///
/// Round-21 L3: ~35 adapters stage the case parameters into the
/// workdir as `valenx_params.json` and read that file back during
/// `collect()` to recover `output_basename` / similar look-ups. The
/// pre-fix reads were bare `fs::read_to_string(workdir.join(...)).ok()?`
/// — a corrupted or poisoned workdir with a multi-GB file would slurp
/// into memory before any parser saw the first byte. 1 MiB is the
/// same shape as [`crate::project::loader::MAX_PROJECT_FILE_BYTES`]
/// and well past anything realistic (a typical staged params file is
/// in the low-KB range).
pub const MAX_ADAPTER_PARAMS_BYTES: u64 = 1024 * 1024;

/// Cap on the bytes a Genetics-Workbench file-loader may consume
/// when the user picks a structure / sequence / tree file via the
/// native file dialog.
///
/// Round-21 H1: eight `std::fs::read_to_string(&path)` sites under
/// `crates/valenx-app/src/genetics/` (biostruct, docking ×2,
/// genomics, md, phylogenetics, qchem, sequence) accepted whatever
/// the rfd-picked path pointed at. A user who picked a multi-GB
/// `.pdb` (or whose saved dialog state pointed at one) would OOM
/// the renderer process before any parser ran. 64 MiB is generous
/// — a 1 M-atom multi-chain PDB ASCII dump is in the low hundreds
/// of MB *only* with multiple solvent shells; honest in-workbench
/// loads are under 10 MB — while still rejecting a hostile or
/// stale path that would OOM the egui frame.
pub const MAX_GENETICS_FILE_BYTES: u64 = 64 * 1024 * 1024;

/// Cap on the bytes a bio-CLI inspector (valenx-fasta /
/// valenx-fastq / valenx-pdb-info / valenx-sam-info /
/// valenx-vcf-info / valenx-blast) may consume from either stdin or
/// a positional path argument.
///
/// Round-21 M4: every bio-inspector CLI pre-fix did a bare
/// `std::io::stdin().read_to_string(&mut buf)` / `fs::read_to_string`
/// against caller-supplied input. A `cat /dev/zero | valenx-fasta
/// inspect -` would OOM the process before any sequence got
/// parsed; the file branch had the same shape against a stale
/// CI-runner path. 256 MiB is generous — even a chromosome-scale
/// FASTA / FASTQ is well under that — while refusing hostile or
/// runaway input.
pub const MAX_BIO_CLI_BYTES: u64 = 256 * 1024 * 1024;

/// Cap on the bytes a `valenx-plugin.toml` manifest may consume
/// when discovered by [`crate`]-side plugin loader.
///
/// Round-21 M5 sister to `valenx_addons::manifest::MAX_ADDON_MANIFEST_BYTES`:
/// `valenx-plugin::load_manifest` walked the per-user / per-system
/// plugin search paths and did a bare `fs::read_to_string(manifest_path)`
/// on every candidate. A poisoned plugins directory containing a
/// multi-MB `valenx-plugin.toml` would OOM during discovery before
/// any plugin actually loaded. 256 KiB matches the addons cap and is
/// well past any honest manifest (which is a few hundred bytes of
/// id / version / capability metadata).
pub const MAX_PLUGIN_MANIFEST_BYTES: u64 = 256 * 1024;

/// Cap on the bytes a `.gltf` JSON manifest may consume.
///
/// Round-21 M1 / L4: `valenx_occt_exchange::gltf2_reader` did a bare
/// `std::fs::read_to_string(path)` on the manifest before parsing.
/// glTF JSON is structured metadata that points at base64-encoded
/// buffers (which can be large); the manifest itself is typically
/// under 1 MiB. 64 MiB is generous while refusing a hostile
/// multi-GB file that would OOM `serde_json::from_str` before the
/// allocator's safety limits triggered. Smaller than
/// `valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES` (256 MiB)
/// because JSON is denser than STEP / IGES text.
pub const MAX_GLTF_JSON_BYTES: u64 = 64 * 1024 * 1024;

/// Cap on the bytes the openEMS adapter may read from a per-probe
/// `.csv` time-series file during `collect()`.
///
/// Round-21 L1: openEMS writes one CSV per probe into the workdir;
/// the adapter pre-fix did `fs::read_to_string(&path)` on each
/// matching `.csv`. A poisoned workdir with a multi-GB CSV would
/// slurp into memory before the line iterator parsed the header. 64 MiB
/// is generous (a typical probe CSV is in the low MB range — one
/// row per FDTD time step) while refusing runaway artefacts.
pub const MAX_OPENEMS_CSV_BYTES: u64 = 64 * 1024 * 1024;

/// Cap on the bytes the MuJoCo adapter may read from the per-run
/// trajectory JSONL written by the staged Python driver.
///
/// Round-21 L2: `valenx-adapter-mujoco::collect` pre-fix did
/// `fs::read_to_string(&ts_path)` on the trajectory file. A long
/// rollout with many DOFs can in principle produce multi-GB JSONL;
/// a corrupted run or hostile workdir could exceed memory before
/// the line parser ran. 64 MiB is generous for honest rollouts
/// (one JSONL row per timestep, ≤ tens of DOFs per row, ≤ hundreds
/// of bytes per row) while refusing runaway captures.
pub const MAX_MUJOCO_TIMESTEP_BYTES: u64 = 64 * 1024 * 1024;

/// Cap on the bytes a `.pdb` (or `.cif`) protein-structure file may
/// consume when read for the collect-label path of a bio adapter
/// (chimerax / alphafold2 / colabfold / biopython / rfdiffusion /
/// rfantibody / proteinmpnn / esmfold / chroma / openmm).
///
/// Round-22 M2: each adapter's `collect()` walks the workdir for
/// `.pdb` (or `.fasta`) artefacts and does a bare
/// `fs::read_to_string(&path)` so it can populate a rich label
/// ("N atoms, M residues") via `valenx_bio::format::pdb::read`. A
/// poisoned (or just runaway) workdir with a multi-GB `.pdb` would
/// slurp into memory before the parser saw the first `ATOM` line.
/// 256 MiB is generous — even a 1-million-atom ribosome with
/// alt-confs is in the low hundreds of MB of PDB text — while
/// rejecting a hostile or runaway artefact.
pub const MAX_PDB_FILE_BYTES: u64 = 256 * 1024 * 1024;

/// Cap on the bytes a `.dcd` (binary trajectory) frame file may
/// consume when read by the collect-label path of a molecular-
/// dynamics adapter (currently `valenx-adapter-mdanalysis`).
///
/// Round-22 M2: pre-fix the MDAnalysis adapter's `.dcd` branch did
/// a bare `fs::read(&path)` on every trajectory in the workdir. A
/// long production run legitimately produces a multi-GB DCD (the
/// CHARMM binary trajectory format packs `nframes × natoms × 12`
/// bytes plus a tiny header), so the cap must accommodate honest
/// MD output. 4 GiB matches the workbench's largest VTK cap and
/// rejects only truly pathological / corrupted files.
pub const MAX_DCD_FRAME_FILE_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// Cap on the bytes the audit-log rotation-genesis sidecar may
/// consume when read back during `read_rotation_genesis`.
///
/// Round-22 L1: pre-fix `valenx_audit::read_rotation_genesis` did
/// a bare `fs::read_to_string(&sidecar)` to recover the
/// `genesis-after-rotation:<sha256>` line the rotator stamped. The
/// sidecar is exactly one line — a 7-byte `genesis-after-rotation:`
/// prefix plus a 64-char hex hash plus a trailing newline ≈ 80
/// bytes — so 1 KiB is more than 10x the honest size while still
/// refusing a hostile or accidentally-mangled sidecar that grew to
/// many GB on disk.
pub const MAX_ROTATION_GENESIS_BYTES: u64 = 1024;

/// Cap on the bytes the LAMMPS adapter may read from `log.lammps`
/// during `collect()` for thermo-column extraction.
///
/// Round-23 named finding: pre-fix the LAMMPS adapter did
/// `fs::read_to_string(&log_path)` on the per-run thermo log. A long
/// run with many thermo lines (or a hostile / runaway artefact) can
/// legitimately push the log into the hundreds of MiB; 256 MiB is
/// generous (LAMMPS thermo blocks are typically a few MiB even for
/// million-atom production runs) while refusing pathological
/// inputs that would OOM the renderer when the report layer tries
/// to chart energy / temperature / pressure vs step.
pub const MAX_LAMMPS_LOG_BYTES: u64 = 256 * 1024 * 1024;

/// Cap on the bytes the PyBaMM adapter may read from the discharge
/// time-series CSV during `collect()`.
///
/// Round-23 named finding: pre-fix the PyBaMM adapter did
/// `std::fs::read_to_string(&ts_path)` on `discharge.csv`. Battery
/// discharge time series at fine sampling can grow large, but
/// honest PyBaMM rollouts (one row per second × several columns of
/// f64 ASCII) top out in the low tens of MiB even for multi-hour
/// runs; 256 MiB is generous while refusing runaway / hostile
/// artefacts that would OOM the charting layer.
pub const MAX_PYBAMM_TIMESERIES_BYTES: u64 = 256 * 1024 * 1024;

/// Cap on the bytes a Gmsh `.msh` mesh file may consume when read
/// by `crate::adapters::mesh::gmsh::msh_parser::parse_file`.
///
/// Round-23 named finding: pre-fix the gmsh msh parser did a bare
/// `fs::read_to_string(path)` on the entire mesh. Production
/// engineering meshes are legitimately huge — a 100M-element
/// tetrahedral mesh in version-4.1 ASCII format easily passes 1 GiB
/// — so the cap must accommodate honest HPC output. 4 GiB matches
/// the workbench's largest VTK cap and rejects only pathological /
/// corrupted files that would OOM before parsing began.
pub const MAX_MSH_FILE_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// Cap on the bytes a Netgen `.vol` mesh file may consume when read
/// by `crate::adapters::mesh::netgen::vol_parser::parse_file`.
///
/// Round-23 named finding: sister to [`MAX_MSH_FILE_BYTES`] —
/// pre-fix the netgen vol parser did a bare
/// `fs::read_to_string(path)` on the mesh. Netgen `.vol` ASCII is
/// denser per element than gmsh, but a high-resolution adaptive
/// mesh can still cross several GiB. 4 GiB matches the gmsh cap.
pub const MAX_VOL_FILE_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// Cap on the bytes the Cantera adapter may read from
/// `summary.json` during `collect()`.
///
/// Round-23 named finding: pre-fix the Cantera summary parser did
/// `fs::read_to_string(path)` on the staged summary. Cantera
/// summaries are tiny — adapter ID, analysis tag, a handful of
/// thermo frames (T/P/H/S/rho per row) — so honest output is well
/// under 100 KiB. 1 MiB is the same shape as
/// [`MAX_ADAPTER_PARAMS_BYTES`] and rejects any poisoned or runaway
/// summary that would OOM `serde_json::from_str`.
pub const MAX_CANTERA_SUMMARY_BYTES: u64 = 1024 * 1024;

/// Cap on the bytes the MuJoCo adapter may read from `summary.json`
/// during `collect()`.
///
/// Round-23 named finding: sister to [`MAX_CANTERA_SUMMARY_BYTES`]
/// — pre-fix the MuJoCo summary parser did `fs::read_to_string(path)`
/// on the staged summary. MuJoCo summaries record only model
/// metadata (name, duration, timestep, step count, nq / nv / nu) —
/// always under 10 KiB. 1 MiB matches the Cantera cap and refuses
/// any hostile or runaway summary.
pub const MAX_MUJOCO_SUMMARY_BYTES: u64 = 1024 * 1024;

/// Cap on the bytes the FreeCAD adapter may read from
/// `summary.json` during `collect()`.
///
/// Round-23 workspace sweep: sister to [`MAX_CANTERA_SUMMARY_BYTES`]
/// — pre-fix the FreeCAD geometry-summary parser did
/// `fs::read_to_string(path)` on the staged summary. FreeCAD
/// summaries record bbox extents + a few scalar geometry metrics
/// (volume / area / centre-of-mass); always under 10 KiB. 1 MiB
/// matches the Cantera cap.
pub const MAX_FREECAD_SUMMARY_BYTES: u64 = 1024 * 1024;

/// Cap on the bytes an `.obj` / `.mtl` mesh-exchange file may
/// consume when read by the round-trip writer / reader in
/// `valenx_occt_exchange::obj_*_extended`.
///
/// Round-23 workspace sweep: pre-fix the obj writer's
/// "splice mtllib into freshly written OBJ" round-trip and the obj
/// reader's path-form variant did bare `fs::read_to_string(path)` on
/// the OBJ + MTL. OBJ for a million-vertex mesh in ASCII is in the
/// low hundreds of MiB; 1 GiB is generous while refusing the
/// `cat /dev/zero` denial of service that would OOM the parser
/// before the line iterator saw any content.
pub const MAX_OBJ_FILE_BYTES: u64 = 1024 * 1024 * 1024;

/// Cap on the bytes a 2D DXF drawing file may consume when read by
/// `valenx_librecad_2d::dxf::read_full`.
///
/// Round-23 workspace sweep: pre-fix the LibreCAD-2D DXF reader did
/// `std::fs::read_to_string(path)` on the entire drawing. AutoCAD
/// R12 ASCII DXF for a complex multi-layer drawing can legitimately
/// be hundreds of MiB; 1 GiB matches the OBJ cap and refuses
/// pathological inputs.
pub const MAX_DXF_FILE_BYTES: u64 = 1024 * 1024 * 1024;

/// Cap on the bytes a `.kicad_pcb` S-expression file may consume
/// when read by `valenx_kicad::parse::import_kicad_pcb`.
///
/// Round-23 workspace sweep: pre-fix the kicad PCB parser did
/// `fs::read_to_string(path)` on the board file. Production PCBs
/// with thousands of footprints can cross 100 MiB of S-expression
/// text; 1 GiB matches the DXF cap.
pub const MAX_KICAD_FILE_BYTES: u64 = 1024 * 1024 * 1024;

/// Cap on the bytes a `.ron` part-library index file may consume
/// when read by `valenx_partlib::library::load_index`.
///
/// Round-23 workspace sweep: pre-fix `load_index` did
/// `fs::read_to_string(&path)` on the index. The index file
/// catalogues every part in a library — for a large org's parts
/// library this can grow but stays in the low MiB. 16 MiB matches
/// [`MAX_DOC_FILE_BYTES`] and refuses pathological inputs.
pub const MAX_PARTLIB_INDEX_BYTES: u64 = 16 * 1024 * 1024;

/// Cap on the bytes an ASCII `.ply` point-cloud file may consume
/// when read by `valenx_reverse::pointcloud::from_ply`.
///
/// Round-23 workspace sweep: pre-fix `from_ply` did
/// `std::fs::read_to_string(path)` on the point cloud. ASCII PLY
/// for a million-vertex scan is in the low hundreds of MiB; 1 GiB
/// matches the OBJ / DXF caps.
pub const MAX_PLY_ASCII_BYTES: u64 = 1024 * 1024 * 1024;

/// Cap on the bytes a JT (Siemens Jupiter Tessellation) file may
/// consume when read by `valenx_occt_exchange::jt_reader::read_jt_model`.
///
/// Round-23 workspace sweep: pre-fix `read_jt_model` did a bare
/// `std::fs::read(path)` on the entire binary container. JT is a
/// dense binary CAD interchange format with the LSG assembly tree,
/// tessellated mesh, and optional XT B-rep payload all in one file
/// — production CAD assemblies (e.g. car body-in-white) can push
/// past 1 GiB. 2 GiB is generous while refusing the cat /dev/zero
/// denial of service.
pub const MAX_JT_FILE_BYTES: u64 = 2 * 1024 * 1024 * 1024;

/// Cap on the bytes any single line of an ASCII OBJ file may consume
/// when parsed by `valenx_mesh::format::obj::read_path`.
///
/// Round-24 H3: pre-fix the OBJ reader did
/// `BufReader::lines().map_while(Result::ok)`, an iterator that
/// allocates an unbounded `String` per `\n`-delimited record and
/// silently drops IO errors. A poisoned file with a single 4 GiB
/// "line" (no newline terminator) would OOM the import before any
/// per-line parsing happened. OBJ lines are usually under 1 KiB
/// (single `v x y z` or `f i j k`); 4 MiB is generous for hand-
/// authored CAD output with very wide `f` polygon lines while
/// refusing the unbounded-line DoS shape.
pub const MAX_OBJ_LINE_BYTES: usize = 4 * 1024 * 1024;

/// Cap on the bytes any single line of a SLURM `.out` / `.err` log
/// may consume when read by `valenx_executor_slurm::read_slurm_log_tail`.
///
/// Round-24 H4: pre-fix the tail walker used
/// `BufReader::lines().map_while(Result::ok)` directly, with no per-
/// line cap. SLURM logs include solver stdout — a runaway MPI rank
/// printing without newlines for the entire wall-clock is a real
/// failure mode. 4 MiB is generous for any single legitimate log
/// entry while refusing multi-GiB pathological records.
pub const MAX_SLURM_LOG_LINE_BYTES: usize = 4 * 1024 * 1024;

/// Read an entire file to a `String`, refusing to allocate more than
/// `cap` bytes regardless of what the filesystem reports.
///
/// Belt-and-braces: both the stat and the reader are bounded, because
/// the file size at stat time can disagree with the file size at read
/// time on a hostile / racing filesystem.
///
/// # Errors
///
/// - [`std::io::ErrorKind::InvalidData`] when the file exceeds `cap`
///   (size known up-front) or when the bytes aren't valid UTF-8.
/// - The underlying IO error for stat / open / read failures.
pub fn read_capped_to_string(path: &Path, cap: usize) -> Result<String, std::io::Error> {
    let meta = std::fs::metadata(path)?;
    if meta.len() > cap as u64 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "file {} exceeds {}-byte cap (actual: {})",
                path.display(),
                cap,
                meta.len()
            ),
        ));
    }
    let mut buf = Vec::new();
    std::fs::File::open(path)?
        .take(cap as u64)
        .read_to_end(&mut buf)?;
    String::from_utf8(buf)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

/// Read an entire file to a `Vec<u8>`, refusing to allocate more than
/// `cap` bytes regardless of what the filesystem reports.
///
/// The binary counterpart to [`read_capped_to_string`] — same
/// belt-and-braces stat + bounded-take pattern, but skips the UTF-8
/// validation step so binary formats (HDR, VTK legacy, .vtu XML with
/// embedded base64 / binary appendices) can be read safely without
/// the helper rejecting them as "not valid UTF-8". Use the string
/// variant for ASCII-only formats; use this for everything else.
///
/// Round-20 M2 + L1: the render-bridge `EnvironmentRef::load` and
/// the SU2 / OpenFOAM CFD `collect` paths were the original callers
/// that motivated this helper — both did bare `std::fs::read(&path)`
/// on user-supplied paths and would have OOM'd on a multi-GB input.
///
/// # Errors
///
/// - [`std::io::ErrorKind::InvalidData`] when the file exceeds `cap`
///   (either at stat time or mid-read on a racing / hostile FS).
/// - The underlying IO error for stat / open / read failures.
pub fn read_capped_to_bytes(path: &Path, cap: u64) -> Result<Vec<u8>, std::io::Error> {
    let meta = std::fs::metadata(path)?;
    if meta.len() > cap {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "file {} exceeds {}-byte cap (actual: {})",
                path.display(),
                cap,
                meta.len()
            ),
        ));
    }
    // Bounded `take(cap + 1)` so a file that grew between stat and
    // open by more than the cap still produces a clean error rather
    // than slurping unbounded. The `+ 1` is the tell — if we got
    // exactly `cap + 1` bytes back, the file outgrew its cap.
    let mut buf = Vec::new();
    std::fs::File::open(path)?
        .take(cap.saturating_add(1))
        .read_to_end(&mut buf)?;
    if buf.len() as u64 > cap {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "file {} grew past {}-byte cap mid-read",
                path.display(),
                cap,
            ),
        ));
    }
    Ok(buf)
}

/// Read at most `cap` bytes from stdin to a `String`. The CLI
/// counterpart to [`read_capped_to_string`] — same bounded-take
/// pattern, but reads from `std::io::stdin().lock()` so a piped
/// `cat /dev/zero | valenx-fasta inspect -` produces a clean error
/// rather than slurping unbounded into memory.
///
/// # Errors
///
/// - [`std::io::ErrorKind::InvalidData`] when stdin produces more
///   than `cap` bytes (we read `cap + 1` and reject any read whose
///   final length exceeds the cap) or when the bytes aren't valid
///   UTF-8.
/// - The underlying IO error for read failures.
pub fn read_capped_stdin_to_string(cap: u64) -> Result<String, std::io::Error> {
    let mut buf = Vec::new();
    std::io::stdin()
        .lock()
        .take(cap.saturating_add(1))
        .read_to_end(&mut buf)?;
    if buf.len() as u64 > cap {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("stdin exceeded {cap}-byte cap"),
        ));
    }
    String::from_utf8(buf).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

/// Read at most `cap` bytes from stdin to a `Vec<u8>`. Binary
/// counterpart to [`read_capped_stdin_to_string`] — same
/// bounded-take pattern, no UTF-8 validation.
///
/// # Errors
///
/// - [`std::io::ErrorKind::InvalidData`] when stdin exceeds `cap`.
/// - The underlying IO error for read failures.
pub fn read_capped_stdin_to_bytes(cap: u64) -> Result<Vec<u8>, std::io::Error> {
    let mut buf = Vec::new();
    std::io::stdin()
        .lock()
        .take(cap.saturating_add(1))
        .read_to_end(&mut buf)?;
    if buf.len() as u64 > cap {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("stdin exceeded {cap}-byte cap"),
        ));
    }
    Ok(buf)
}

/// Iterator over lines of a `BufRead`, with each line bounded at
/// `max_per_line` bytes — the audit-log helper that lets a line
/// reader survive a pathological input file with no `\n` for many
/// GB.
///
/// Mirrors the pattern audit-log scanning uses (sister to the
/// round-20 M3 `read_capped_lines` migration) so callers that don't
/// want the audit-specific framing can reuse the same shape. Each
/// `next()` yields up to `max_per_line` bytes from the next line
/// (including any embedded `\n` until the boundary). Lines longer
/// than `max_per_line` produce a single `Err(InvalidData)` and the
/// iterator stops — callers can still surface a partial result for
/// everything read before the offending line.
///
/// Use [`read_capped_to_string`] for files that fit in memory; this
/// helper is for cases where the file *might* be larger than RAM
/// (long-running audit logs, streaming pipeline output) but each
/// individual line is bounded.
pub fn read_capped_lines_bounded<R: BufRead>(
    mut reader: R,
    max_per_line: usize,
) -> impl Iterator<Item = std::io::Result<Vec<u8>>> {
    let mut done = false;
    std::iter::from_fn(move || {
        if done {
            return None;
        }
        let mut buf = Vec::with_capacity(64);
        // Read up to max_per_line+1 bytes, stop at newline.
        let cap = (max_per_line as u64).saturating_add(1);
        let mut limited = (&mut reader).take(cap);
        match limited.read_until(b'\n', &mut buf) {
            Ok(0) => {
                done = true;
                None
            }
            Ok(_n) => {
                if buf.len() > max_per_line {
                    done = true;
                    return Some(Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "line exceeded {max_per_line}-byte cap (possible \
                             missing newline or hostile input)"
                        ),
                    )));
                }
                Some(Ok(buf))
            }
            Err(e) => {
                done = true;
                Some(Err(e))
            }
        }
    })
}

/// Round-27 STRUCTURAL: canonical atomic-write helper used by every
/// crate in the workspace that publishes a file via the
/// write-sidecar → fsync → rename-over-target pattern.
///
/// Replaces 4 near-identical inlined copies (`valenx-app::state_paths::atomic_write`,
/// `valenx-dock::runner::atomic_write_pdbqt`, `valenx-crash-reporter::atomic_write_bytes`,
/// `valenx-render-bridge::persist::write_to`). All 4 sites now
/// delegate here, so a single fix to the crash-safety invariants
/// (parent-fsync, O_NOFOLLOW, counter-based sidecar names,
/// pre-rename fsync) lands everywhere at once instead of having to
/// be back-ported across 4 copies.
///
/// ## Crash-safety invariants
///
/// 1. **Unique sidecar name** — `<basename>.tmp.<pid>.<counter>`
///    where `counter` is a process-monotonic [`AtomicU64`]. Two
///    concurrent writers (even at the same nanosecond) own distinct
///    sidecar paths. Solves R24 M1 + R25 M3 + R26 H1.
/// 2. **`O_NOFOLLOW` (Unix) / `FILE_FLAG_OPEN_REPARSE_POINT`
///    (Windows)** on the sidecar open — defence in depth against an
///    attacker who pre-seeds the sidecar path as a symlink.
/// 3. **`create_new(true)`** — kernel refuses to open a pre-existing
///    sidecar. Belt-and-braces against the unique-counter naming.
/// 4. **`sync_all()` before the rename** — ensures the file data is
///    on durable storage before the rename publishes the new
///    contents. Without this a power loss between `write_all` and
///    `rename` could land a zero-length file as the durable state.
/// 5. **`fs::rename` for atomic publication** — POSIX `rename(2)` is
///    atomic; Windows `MoveFileEx`/`MOVEFILE_REPLACE_EXISTING` (set
///    implicitly by Rust's `fs::rename` when the destination
///    exists) is atomic too. A concurrent reader sees either the
///    old file or the new file, never a partial write.
/// 6. **Parent directory fsync (Unix)** — after rename, fsync the
///    parent directory so the dentry update (the new file name
///    binding) is durable, not just the file data. NTFS journals
///    metadata via the USN journal so the Windows path is a no-op.
/// 7. **Sidecar cleanup on rename failure** — best-effort
///    `remove_file` on the tmp path so a failed rename doesn't
///    accumulate orphans under the parent dir.
///
/// ## Errors
///
/// - [`std::io::ErrorKind::InvalidInput`] when the target path has
///   no filename component (e.g. a bare `/` path).
/// - The underlying IO error for parent-dir creation, sidecar open,
///   write, sync, or rename failures.
pub fn atomic_write_bytes(path: &Path, contents: &[u8]) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp_path = canonical_sidecar_path(path)?;
    let mut f = open_sidecar_no_follow(&tmp_path)?;
    // Round-28 L1 / R29 M2: RAII guard so a panic inside `write_all`
    // (rare but possible — a `Vec<u8>::extend` allocator failure, or
    // any future closure path) unwinds through `Drop` and removes the
    // sidecar. Pre-fix the manual `let _ = remove_file(...)`-on-every-
    // error-path approach left the sidecar behind on panic. R29 M2:
    // construct the guard on the line IMMEDIATELY after the open so
    // the sidecar is never on disk without an armed cleanup guard —
    // closing even the tiny no-op window between open and guard.
    let mut guard = SidecarGuard {
        tmp_path: &tmp_path,
        disarmed: false,
    };
    f.write_all(contents)?;
    f.sync_all()?;
    drop(f);
    std::fs::rename(&tmp_path, path)?;
    guard.disarmed = true;
    parent_dir_fsync(path);
    Ok(())
}

/// Round-28 L1: RAII sidecar cleanup. `Drop` removes the sidecar
/// unless [`Self::disarmed`] is set. Set `disarmed = true` only
/// after the rename succeeds — until then any path out of the
/// caller (early return, panic, `?` unwinding) must remove the
/// orphan sidecar.
///
/// Lifetime'd against the caller's `tmp_path` so the guard can hold
/// a borrow without an owned `PathBuf` allocation.
struct SidecarGuard<'a> {
    tmp_path: &'a Path,
    disarmed: bool,
}

impl Drop for SidecarGuard<'_> {
    fn drop(&mut self) {
        if !self.disarmed {
            let _ = std::fs::remove_file(self.tmp_path);
        }
    }
}

/// String-counterpart of [`atomic_write_bytes`] — same crash-safety
/// invariants, takes a `&str` instead of `&[u8]`.
pub fn atomic_write_str(path: &Path, contents: &str) -> std::io::Result<()> {
    atomic_write_bytes(path, contents.as_bytes())
}

/// Streaming counterpart of [`atomic_write_bytes`] — same crash-
/// safety invariants (unique sidecar, O_NOFOLLOW, fsync-before-
/// rename, parent fsync), but the contents are produced lazily by
/// the supplied closure instead of being materialised in memory
/// up-front. Use this when the file is large enough that a
/// monolithic `Vec<u8>` would be wasteful (e.g. multi-MiB IFC
/// exports).
///
/// The closure receives a `&mut BufWriter<File>` and is expected to
/// emit the file contents via `write_all` etc. The 8 KiB BufWriter
/// window is the only memory overhead; no monolithic buffer is
/// allocated.
///
/// The BufWriter is flushed before the file is `sync_all`'d and
/// renamed, so any IO error surfaced inside the closure short-
/// circuits the publication.
pub fn atomic_write_streaming<F>(path: &Path, write_fn: F) -> std::io::Result<()>
where
    F: FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
{
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp_path = canonical_sidecar_path(path)?;
    let f = open_sidecar_no_follow(&tmp_path)?;
    // Round-28 L1 / R29 M2: RAII guard so a panic inside `write_fn` (a
    // careless caller's IFC writer is the headline concern) unwinds
    // through `Drop` and removes the sidecar. Pre-fix the manual
    // `let _ = remove_file(...)`-on-every-error-path approach left
    // the sidecar behind on panic. The guard's `Drop` also covers
    // the `BufWriter::into_inner` failure path (returns
    // `IntoInnerError` if a queued flush errored). R29 M2: arm the
    // guard on the line IMMEDIATELY after the open — before wrapping
    // `f` in a BufWriter — so the sidecar is never on disk without an
    // armed cleanup guard.
    //
    // Round-28 L2: the previous `let _ = bw.flush();` inside the
    // closure error path was cosmetic — `drop(bw)` on the next line
    // already flushes via `BufWriter`'s `Drop` impl (ignoring
    // errors, same as `let _ = bw.flush()` did explicitly). It has
    // been removed.
    let mut guard = SidecarGuard {
        tmp_path: &tmp_path,
        disarmed: false,
    };
    let mut bw = BufWriter::new(f);
    write_fn(&mut bw)?;
    bw.flush()?;
    let f = bw.into_inner().map_err(|e| e.into_error())?;
    f.sync_all()?;
    drop(f);
    std::fs::rename(&tmp_path, path)?;
    guard.disarmed = true;
    parent_dir_fsync(path);
    Ok(())
}

/// Compute the canonical sidecar path
/// `<basename>.tmp.<pid>.<counter>` for `path`. Errors with
/// [`std::io::ErrorKind::InvalidInput`] when `path` has no filename
/// component.
fn canonical_sidecar_path(path: &Path) -> std::io::Result<std::path::PathBuf> {
    let name = match path.file_name().and_then(|s| s.to_str()) {
        Some(s) => s.to_string(),
        None => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "atomic_write: target path has no filename component",
            ));
        }
    };
    let pid = std::process::id();
    let counter = ATOMIC_WRITE_COUNTER.fetch_add(1, Ordering::Relaxed);
    let tmp_name = format!("{name}.tmp.{pid}.{counter}");
    Ok(path.with_file_name(tmp_name))
}

/// Open the canonical sidecar for write with `O_NOFOLLOW | O_EXCL`
/// (Unix) or `FILE_FLAG_OPEN_REPARSE_POINT | create_new` (Windows).
///
/// On Unix, `open(2)` returns `ELOOP` if the leaf is a symlink — the
/// kernel refuses the open outright.
///
/// On Windows, `FILE_FLAG_OPEN_REPARSE_POINT` makes CreateFileW
/// open the reparse point itself rather than following it; we then
/// run [`refuse_symlink_after_open`] to mirror Unix's `ELOOP`
/// behaviour for symlinked sidecars (the Windows flag's semantic
/// gap closed by R26 M2).
fn open_sidecar_no_follow(path: &Path) -> std::io::Result<File> {
    let mut opts = OpenOptions::new();
    opts.write(true).create_new(true);
    #[cfg(unix)]
    opts.custom_flags(libc::O_NOFOLLOW);
    #[cfg(windows)]
    opts.custom_flags(FILE_FLAG_OPEN_REPARSE_POINT);
    let f = opts.open(path)?;
    #[cfg(windows)]
    refuse_symlink_after_open(&f)?;
    Ok(f)
}

/// Round-26 M2 (lifted from dock-runner): on Windows
/// `FILE_FLAG_OPEN_REPARSE_POINT` opens the reparse point itself but
/// doesn't fail with an `ELOOP`-equivalent when the leaf IS a
/// symlink — it just returns a handle whose reads/writes target the
/// raw reparse data. Close the semantic gap with a post-open
/// metadata check: if the leaf was a symlink we refuse it with an
/// `Other`-kind error mirroring Unix's `ELOOP` behaviour. On Unix
/// the kernel already rejected the open with `ELOOP` before we got
/// here, so the helper is `#[cfg]`'d out entirely.
#[cfg(windows)]
fn refuse_symlink_after_open(f: &File) -> std::io::Result<()> {
    if f.metadata()?.file_type().is_symlink() {
        Err(std::io::Error::other(
            "leaf is a symlink (refused — FILE_FLAG_OPEN_REPARSE_POINT \
             would otherwise return a reparse-point handle)",
        ))
    } else {
        Ok(())
    }
}

/// Round-26 M3 (lifted from state_paths): fsync the parent directory
/// after the rename so the directory entry update (the new file
/// name binding) is durable, not just the file data.
///
/// Without this a power loss immediately after rename could leave
/// the file data on disk but the dentry update lost in the journal
/// — readers would see the OLD file at the target path, not the new
/// one.
///
/// On Windows the equivalent semantic (NTFS metadata journalled
/// alongside file data via the USN journal) is satisfied by the
/// rename itself, so we skip the directory fsync. The `File::open`
/// / `sync_all` pattern below works on Unix only — Windows doesn't
/// let you open a directory as a regular `File` handle for
/// synchronisation.
///
/// Best-effort: a failed sync of the directory entry is not
/// actionable here (rename already succeeded, the file content is
/// durable). Operators that care about crash-resilience-at-rename-
/// boundary should use a journalled filesystem.
fn parent_dir_fsync(path: &Path) {
    #[cfg(unix)]
    if let Some(parent) = path.parent() {
        if let Ok(dir) = File::open(parent) {
            let _ = dir.sync_all();
        }
    }
    #[cfg(not(unix))]
    {
        let _ = path;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// RED→GREEN: a file larger than the cap is rejected with
    /// `InvalidData` instead of being slurped into memory.
    #[test]
    fn rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_io_caps_oversize.dat");
        let mut f = std::fs::File::create(&tmp).unwrap();
        // 1 KiB cap, 4 KiB file → must reject.
        f.write_all(&vec![b'x'; 4096]).unwrap();
        drop(f);
        let err = read_capped_to_string(&tmp, 1024).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn accepts_undersize_file() {
        let tmp = std::env::temp_dir().join("valenx_io_caps_ok.dat");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(b"hello world").unwrap();
        drop(f);
        let s = read_capped_to_string(&tmp, 1024).unwrap();
        assert_eq!(s, "hello world");
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn rejects_non_utf8() {
        let tmp = std::env::temp_dir().join("valenx_io_caps_bad_utf8.dat");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&[0xFF, 0xFE, 0xFD]).unwrap();
        drop(f);
        let err = read_capped_to_string(&tmp, 1024).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn caps_are_sensible() {
        assert_eq!(MAX_DOC_FILE_BYTES, 16 * 1024 * 1024);
        assert_eq!(MAX_MESH_JSON_BYTES, 64 * 1024 * 1024);
        assert_eq!(MAX_RESULTS_JSON_BYTES, 64 * 1024 * 1024);
        assert_eq!(MAX_FRD_FILE_BYTES, 256 * 1024 * 1024);
        assert_eq!(MAX_PVD_FILE_BYTES, 16 * 1024 * 1024);
        assert_eq!(MAX_SWEEP_SIBLINGS, 100_000);
        // Round-20 caps.
        assert_eq!(MAX_PDBQT_FILE_BYTES, 64 * 1024 * 1024);
        assert_eq!(MAX_HDR_FILE_BYTES, 256u64 * 1024 * 1024);
        assert_eq!(MAX_VTK_FILE_BYTES, 4u64 * 1024 * 1024 * 1024);
        // Round-21 caps.
        assert_eq!(MAX_ADAPTER_PARAMS_BYTES, 1024u64 * 1024);
        assert_eq!(MAX_GENETICS_FILE_BYTES, 64u64 * 1024 * 1024);
        assert_eq!(MAX_BIO_CLI_BYTES, 256u64 * 1024 * 1024);
        assert_eq!(MAX_PLUGIN_MANIFEST_BYTES, 256u64 * 1024);
        assert_eq!(MAX_GLTF_JSON_BYTES, 64u64 * 1024 * 1024);
        assert_eq!(MAX_OPENEMS_CSV_BYTES, 64u64 * 1024 * 1024);
        assert_eq!(MAX_MUJOCO_TIMESTEP_BYTES, 64u64 * 1024 * 1024);
        // Round-22 caps.
        assert_eq!(MAX_PDB_FILE_BYTES, 256u64 * 1024 * 1024);
        assert_eq!(MAX_DCD_FRAME_FILE_BYTES, 4u64 * 1024 * 1024 * 1024);
        assert_eq!(MAX_ROTATION_GENESIS_BYTES, 1024u64);
        // Round-23 caps.
        assert_eq!(MAX_LAMMPS_LOG_BYTES, 256u64 * 1024 * 1024);
        assert_eq!(MAX_PYBAMM_TIMESERIES_BYTES, 256u64 * 1024 * 1024);
        assert_eq!(MAX_MSH_FILE_BYTES, 4u64 * 1024 * 1024 * 1024);
        assert_eq!(MAX_VOL_FILE_BYTES, 4u64 * 1024 * 1024 * 1024);
        assert_eq!(MAX_CANTERA_SUMMARY_BYTES, 1024u64 * 1024);
        assert_eq!(MAX_MUJOCO_SUMMARY_BYTES, 1024u64 * 1024);
        assert_eq!(MAX_FREECAD_SUMMARY_BYTES, 1024u64 * 1024);
        assert_eq!(MAX_OBJ_FILE_BYTES, 1024u64 * 1024 * 1024);
        assert_eq!(MAX_DXF_FILE_BYTES, 1024u64 * 1024 * 1024);
        assert_eq!(MAX_KICAD_FILE_BYTES, 1024u64 * 1024 * 1024);
        assert_eq!(MAX_PARTLIB_INDEX_BYTES, 16u64 * 1024 * 1024);
        assert_eq!(MAX_PLY_ASCII_BYTES, 1024u64 * 1024 * 1024);
        assert_eq!(MAX_JT_FILE_BYTES, 2u64 * 1024 * 1024 * 1024);
    }

    /// Round-21: `read_capped_lines_bounded` yields each line up to
    /// the per-line cap. Pathological "no newline for 10 GB" input
    /// is refused mid-iteration without slurping unbounded.
    #[test]
    fn read_capped_lines_bounded_caps_long_line() {
        use std::io::Cursor;
        // 4-byte cap, one well-formed line then a 16-byte runaway.
        let input = b"abc\nxxxxxxxxxxxxxxxxx";
        let mut it = read_capped_lines_bounded(Cursor::new(&input[..]), 4);
        let first = it.next().unwrap().unwrap();
        assert_eq!(first, b"abc\n");
        let second = it.next().unwrap();
        assert!(second.is_err(), "second line must trip the cap");
        assert_eq!(
            second.err().unwrap().kind(),
            std::io::ErrorKind::InvalidData
        );
        // Iterator stops after the error.
        assert!(it.next().is_none());
    }

    /// Round-21: bounded-line iterator passes through honest input
    /// (every line under the cap, mixed final-newline shape).
    #[test]
    fn read_capped_lines_bounded_passes_clean_input() {
        use std::io::Cursor;
        let input = b"alpha\nbeta\ngamma";
        let lines: Vec<Vec<u8>> = read_capped_lines_bounded(Cursor::new(&input[..]), 1024)
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], b"alpha\n");
        assert_eq!(lines[1], b"beta\n");
        assert_eq!(lines[2], b"gamma");
    }

    /// Round-20: binary variant rejects oversized files exactly like
    /// the string variant — same stat-then-bounded-take pattern.
    #[test]
    fn read_capped_to_bytes_rejects_oversize() {
        let tmp = std::env::temp_dir().join("valenx_io_caps_oversize_bin.dat");
        let mut f = std::fs::File::create(&tmp).unwrap();
        // 1 KiB cap, 4 KiB file → must reject.
        f.write_all(&vec![0xFF; 4096]).unwrap();
        drop(f);
        let err = read_capped_to_bytes(&tmp, 1024).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-20: binary variant accepts non-UTF-8 bytes (unlike the
    /// string variant which rejects them) — this is the whole point
    /// of having a separate `read_capped_to_bytes`.
    #[test]
    fn read_capped_to_bytes_accepts_non_utf8() {
        let tmp = std::env::temp_dir().join("valenx_io_caps_binary.dat");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&[0xFF, 0xFE, 0xFD, 0xFC]).unwrap();
        drop(f);
        let bytes = read_capped_to_bytes(&tmp, 1024).unwrap();
        assert_eq!(bytes, vec![0xFF, 0xFE, 0xFD, 0xFC]);
        let _ = std::fs::remove_file(&tmp);
    }

    /// RED→GREEN (round-27 STRUCTURAL): the canonical
    /// `atomic_write_bytes` helper round-trips a small payload.
    /// Smoke-tests that the parent-dir auto-create, sidecar open,
    /// fsync, rename, and parent-fsync legs all land cleanly.
    #[test]
    fn atomic_write_bytes_round_trips_round27() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-r27-aw-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let path = dir.join("nested").join("data.bin");
        // Parent dir doesn't exist yet — the helper must create it.
        atomic_write_bytes(&path, b"hello world").expect("atomic_write_bytes");
        assert_eq!(std::fs::read(&path).unwrap(), b"hello world");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// RED→GREEN (round-27 STRUCTURAL): `atomic_write_str` is the
    /// string-flavoured wrapper — must produce the same bytes the
    /// `&[u8]` form would.
    #[test]
    fn atomic_write_str_round_trips_round27() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-r27-aws-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("data.txt");
        atomic_write_str(&path, "ünïcødé").expect("atomic_write_str");
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "ünïcødé");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// RED→GREEN (round-27 STRUCTURAL): 10 threads concurrently
    /// call `atomic_write_bytes` against the same target with
    /// distinct payloads. Pre-canonical inlining used a per-callsite
    /// counter; the helper here is workspace-wide and must remain
    /// safe under concurrent use. Each writer must own a unique
    /// sidecar (no `AlreadyExists` errors from `create_new`) and
    /// the final file must contain exactly ONE writer's content
    /// (no interleaved bytes).
    #[test]
    fn atomic_write_bytes_concurrent_writes_pick_winner_round27() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        let dir = std::env::temp_dir().join(format!(
            "valenx-r27-concurrent-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let target = dir.join("contended.bin");
        const N: usize = 10;
        let payloads: Vec<Vec<u8>> =
            (0..N).map(|i| format!("payload-{i}").into_bytes()).collect();
        let barrier = Arc::new(Barrier::new(N));
        let mut handles = Vec::with_capacity(N);
        for payload in payloads.clone() {
            let target = target.clone();
            let barrier = Arc::clone(&barrier);
            handles.push(thread::spawn(move || {
                barrier.wait();
                atomic_write_bytes(&target, &payload)
            }));
        }
        let mut ok = 0;
        for h in handles {
            match h.join().unwrap() {
                Ok(()) => ok += 1,
                Err(e) => panic!("concurrent atomic_write_bytes failed: {e}"),
            }
        }
        assert_eq!(ok, N, "all {N} writers must succeed");
        let final_contents = std::fs::read(&target).unwrap();
        assert!(
            payloads.iter().any(|p| p == &final_contents),
            "final contents must equal exactly ONE input (no interleaving); got {} bytes",
            final_contents.len(),
        );
        // No leaked sidecars.
        let mut orphans = Vec::new();
        for entry in std::fs::read_dir(&dir).unwrap().flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            if name.contains(".tmp.") {
                orphans.push(name);
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
        assert!(orphans.is_empty(), "orphaned tmp files: {orphans:?}");
    }

    /// RED→GREEN (round-27 STRUCTURAL): 1000 atomic writes in a
    /// tight loop, each to a distinct path in the same dir. The
    /// process-monotonic counter guarantees a strictly unique
    /// sidecar per call, even at the same nanosecond.
    #[test]
    fn atomic_write_bytes_counter_prevents_collision_round27() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-r27-counter-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        const N: usize = 1000;
        for i in 0..N {
            let p = dir.join(format!("file-{i}.bin"));
            atomic_write_bytes(&p, format!("payload-{i}").as_bytes())
                .expect("atomic_write_bytes in loop");
        }
        // Spot-check a handful of payloads.
        for i in [0, 1, 500, 999] {
            let body = std::fs::read(dir.join(format!("file-{i}.bin"))).unwrap();
            assert_eq!(body, format!("payload-{i}").as_bytes());
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// RED→GREEN (round-27 STRUCTURAL): the streaming variant
    /// produces the same bytes the in-memory variant would, but
    /// without materialising a monolithic buffer.
    #[test]
    fn atomic_write_streaming_round_trips_round27() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-r27-stream-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("streamed.txt");
        atomic_write_streaming(&path, |bw| {
            bw.write_all(b"line 1\n")?;
            bw.write_all(b"line 2\n")?;
            bw.write_all(b"line 3\n")?;
            Ok(())
        })
        .expect("atomic_write_streaming");
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "line 1\nline 2\nline 3\n",
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// RED→GREEN (round-27 STRUCTURAL, Unix only): pre-creating the
    /// sidecar path as a SYMLINK must cause `atomic_write_bytes` to
    /// return an error — the `O_NOFOLLOW` flag refuses to traverse a
    /// symlinked sidecar (the sidecar shouldn't exist at all, since
    /// the counter-based name is unique, but defence-in-depth
    /// catches an attacker who races a symlink in).
    #[cfg(unix)]
    #[test]
    fn atomic_write_bytes_refuses_symlink_sidecar_round27() {
        use std::os::unix::fs::symlink;
        let dir = std::env::temp_dir().join(format!(
            "valenx-r27-symlink-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let target = dir.join("data.bin");
        // Pre-seed the FIRST sidecar path the helper will compute
        // for this call — we can't know the counter ahead of time
        // because it's process-global, but we can run a probe call
        // first to advance it, then read the next-counter sidecar
        // shape and plant a symlink there.
        //
        // Simpler test surrogate: call `atomic_write_bytes` once to
        // observe the failure surface — the cleanest demonstration
        // is to plant the symlink at the exact name the helper will
        // try next. We compute it the same way the helper does.
        let counter_before = ATOMIC_WRITE_COUNTER.load(Ordering::Relaxed);
        let pid = std::process::id();
        let sidecar_name = format!("data.bin.tmp.{pid}.{counter_before}");
        let sidecar = dir.join(&sidecar_name);
        let bait = dir.join("bait.bin");
        std::fs::write(&bait, b"bait").unwrap();
        // Plant a symlink at the exact path `atomic_write_bytes`
        // will try to create_new — open(2) with O_NOFOLLOW must
        // refuse it.
        symlink(&bait, &sidecar).unwrap();
        let err = atomic_write_bytes(&target, b"payload").expect_err(
            "atomic_write_bytes must refuse a symlinked sidecar via O_NOFOLLOW",
        );
        // O_NOFOLLOW returns ELOOP on Linux; on macOS it's EMLINK.
        // Either way the kernel surfaces an `Err`; we don't pin the
        // specific kind because POSIX permits a few variations.
        assert!(
            !matches!(err.kind(), std::io::ErrorKind::NotFound),
            "got {err:?} ({:?}); expected a follow-related error",
            err.kind(),
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// RED→GREEN (round-27 STRUCTURAL): `atomic_write_bytes` returns
    /// `InvalidInput` when the target path has no filename
    /// component. Pins the helper's defensive case so a future
    /// refactor that drops the check fails this test instead of
    /// silently writing to a degenerate path.
    #[test]
    fn atomic_write_bytes_rejects_dirless_path_round27() {
        #[cfg(unix)]
        let p = std::path::Path::new("/");
        #[cfg(windows)]
        let p = std::path::Path::new("C:\\");
        let err = atomic_write_bytes(p, b"x").expect_err("dirless path must fail");
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    }

    /// Round-28 L1 RED→GREEN — pre-fix the manual `let _ =
    /// remove_file(...)` cleanup at each error path inside
    /// `atomic_write_streaming` was only reached on `Err`-via-`?`.
    /// A panic inside `write_fn` (e.g. an over-eager `unwrap()` in a
    /// caller's IFC writer) skipped the cleanup and left the
    /// sidecar on disk forever. Post-fix the `SidecarGuard`'s
    /// `Drop` runs during unwinding too, so we expect the sidecar
    /// to be gone after `catch_unwind` returns.
    ///
    /// Observability: we count the dir entries containing `.tmp.`
    /// in their name before and after. Pre-fix there's exactly one
    /// orphan; post-fix there are zero.
    #[test]
    fn atomic_write_streaming_cleans_sidecar_on_panic_round28_l1() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-r28-l1-stream-panic-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("payload.txt");

        let count_sidecars = || {
            std::fs::read_dir(&dir)
                .unwrap()
                .filter_map(|e| e.ok())
                .filter(|e| {
                    e.file_name()
                        .to_string_lossy()
                        .contains(".tmp.")
                })
                .count()
        };

        let before = count_sidecars();
        assert_eq!(before, 0, "no sidecars before the panic'ing write");

        let path_clone = path.clone();
        let result = std::panic::catch_unwind(move || {
            let _ = atomic_write_streaming(&path_clone, |bw| {
                bw.write_all(b"some bytes before the panic\n")?;
                // Simulate a careless caller — panic midway
                // through the write. Pre-fix the sidecar leaks;
                // post-fix the RAII guard cleans it up during
                // unwinding.
                panic!("panic from inside write_fn — sidecar must be cleaned up");
            });
        });
        assert!(result.is_err(), "catch_unwind must surface the panic");

        let after = count_sidecars();
        assert_eq!(
            after, 0,
            "RAII guard must remove the sidecar during unwinding (found {after} orphans)",
        );
        // The target file must not exist (we never reached rename).
        assert!(
            !path.exists(),
            "target file must not exist after a panic during the write",
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Round-28 L1 sister — same RAII guarantee for
    /// `atomic_write_bytes`. A panic between sidecar open and
    /// rename (rare — `Vec::write_all` doesn't naturally panic, but
    /// future closure variants might) must clean up the sidecar.
    ///
    /// We can't naturally panic inside `atomic_write_bytes` since
    /// it doesn't take a closure. Instead this test exercises the
    /// happy path AND uses the absence of sidecar orphans on the
    /// streaming variant's panic test as the structural guarantee
    /// for both. We pin a structural invariant: after a successful
    /// write, no `.tmp.` sidecars survive in the parent dir.
    #[test]
    fn atomic_write_bytes_leaves_no_sidecar_on_success_round28_l1() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-r28-l1-bytes-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("payload.bin");
        atomic_write_bytes(&path, b"hello").expect("happy-path write");
        // No `.tmp.` survivors.
        let orphans: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().contains(".tmp."))
            .collect();
        assert!(
            orphans.is_empty(),
            "no `.tmp.` sidecar survivors after success — found {} orphans",
            orphans.len(),
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
