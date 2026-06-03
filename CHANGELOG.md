# Changelog

All notable changes to Valenx are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
(see [POLICIES.md](./POLICIES.md) for what SemVer covers here).

---

## [Unreleased]

### Security

**Rounds 15–33 — sustained robustness & untrusted-input hardening
(summary).** A long series of code-review rounds drove the
review loop toward zero open findings. Rather than enumerate every
round, the major classes addressed were:

- **Atomic-write migration.** Every persistence write surface was moved
  onto the canonical `valenx_core::io_caps::atomic_write_*` helpers
  (sidecar-temp with O_NOFOLLOW-style protection, fsync-before-rename,
  parent-dir fsync on Unix), replacing bare `std::fs::write` calls that
  silently followed leaf symlinks and were non-atomic. A workspace-wide
  machine-guard test pins the contract so new bare writes are caught.
- **Path / symlink / TOCTOU hardening.** Confined-join helpers, symlink
  rejection, and shared-cap constants were applied consistently across
  adapters and loaders so a poisoned project or case directory cannot
  escape its sandbox or race a double-read.
- **Parser robustness.** Recursion-depth caps, char-boundary-safe string
  slicing (multibyte-safe `split_at` / column extraction), and
  declared-count allocation caps were swept across the text and binary
  parsers (STEP/IGES, NEXUS/PHYLIP, OBJ/PLY/MSH, VTK, and the bio
  sequence/structure formats) so a hostile header or non-ASCII payload
  is rejected instead of panicking or OOMing.
- **Validated-deserialize.** Loaded documents whose internal index /
  handle references point out of range are now rejected at load with a
  typed error rather than panicking "index out of bounds" the first
  time they are consumed (e.g. sketch entity variable handles validated
  on `.valenx` project / sketch load; gmsh `.msh` element node tags
  bounds-checked at parse).

**Round 14 — addons copy_dir helper bypass + vtk_legacy DoS cap +
sweep.rs 3 caps (R13 carry-over) + IFC string injection + format_g23
NaN/inf check + LHS uncapped n_samples + addons manifest cap tighten
+ rbac_io double-read shared cap + LocalExecutor kill_on_drop +
CITATION.cff regression fix + AF2 absl form + 2 doc tag-push fixes.**

**High-severity:**

- **H1 — Addons `copy_dir_recursive` helper bypass (R8 sister gap).**
  `valenx-addons/src/install.rs` kept its own copy of the recursive
  copy walker. Round-6 had hardened the shared
  `valenx_core::adapter_helpers::copy_dir_recursive` (symlink
  rejection, depth cap, per-file size cap) but the addons crate never
  migrated. Now delegates to the shared helper so any future
  hardening propagates automatically. Same migration the precice
  adapter did in round-7.
- **H2 — `vtk_legacy.rs` DoS cap (R11 sister gap).** New
  `MAX_VTK_LEGACY_POINTS = 256_000_000` cap (matches the round-11
  `MAX_VTU_POINTS` on the XML VTK parser). `read_points_block`, the
  CELLS path, the CELL_TYPES path, and `read_typed_block` now
  `checked_mul` + cap-check before any allocation. A hostile header
  like `POINTS 10000000000 float` no longer overflows usize on 64-bit
  hosts. New `ParseError::TooLarge { what, count, max }` variant.
- **H3 — `sweep.rs` `assemble_sweep_dataset` 3 caps (R13 carry-over).**
  Three sites in the post-sweep aggregator: (a) the case.toml re-read
  now uses `read_capped_to_string(MAX_PROJECT_FILE_BYTES)` (1 MiB),
  matching the round-12 M6 fix on `sweep_selected_case`; (b) the
  per-derived-case `results.json` read now uses
  `read_capped_to_string(MAX_RESULTS_JSON_BYTES)` (64 MiB); (c) the
  sibling-dir enumeration caps at `MAX_SWEEP_SIBLINGS = 100_000` so a
  poisoned sweep parent dir with millions of subdirs can't fill the
  Vec before any reasonability check fires.

**Medium-severity:**

- **M4 — IFC `ifc_str` newline / backslash injection.** Pre-fix only
  `'` was escaped; a hostile name like `"X');\n#9999=IFCPROJECT(..."`
  survived single-quote escape and let the embedded `);\n#9999=`
  break out of the entity wrapper and inject sibling IFC entities.
  Now strips `\n` / `\r` and escapes `\\` per Part 21. Mirrors the
  round-12 CAM `CommentStyle::wrap` sanitiser pattern.
- **M5 — `format_g23` NaN / inf check (R3 sister gap, machine safety).**
  Round-3 fixed `format_g1`, `format_g0`, `format_g0_strict`,
  `format_g1_5ax`. `format_g23` (arc moves G2/G3) used `!(feed > 0.0)`
  which let `f64::INFINITY` through (`+inf > 0.0` is `true`) and
  never validated start / end / centre_xy components. Now uses
  `!feed.is_finite() || feed <= 0.0` for the feed plus per-vector
  `is_finite()` checks. Same controller-safety class as the round-3
  set.
- **M6 — LHS uncapped `n_samples` (R4 sister gap).** New
  `MAX_LHS_SAMPLES = 1_000_000` cap on the LHS optimizer's
  `n_samples`. Pre-fix `n_samples = 10_000_000_000` flowed straight
  into the strata-Vec allocation and OOMed. New
  `OptimizerError::TooManySamples { optimizer, requested, cap }`
  variant.
- **M7 — Addons `read_manifest_at` cap tightened.** Lowered from
  `MAX_MANIFEST_BYTES = 1 MiB` to `MAX_ADDON_MANIFEST_BYTES = 256 KiB`.
  Manifest schema is a dozen fields — even a maxed-out addon fits
  under 32 KiB. Old name is `#[deprecated]` aliased so existing
  callers keep working through one release cycle.
- **M8 — `rbac_io::rbac_override_from_project_toml` double-read TOCTOU.**
  Pre-fix the overlay reader used a local `MAX_PROJECT_TOML_BYTES`
  constant that happened to match the project loader's cap; staying
  in sync was manual. Now imports
  `valenx_core::project::loader::MAX_PROJECT_FILE_BYTES` directly so
  loader + overlay share a single source of truth. Cleaner long-term
  fix (refactor signature to take a `&toml::Value` from the parsed
  loader output) deferred — documented in code comments as option (b).
- **M9 — `LocalExecutor::Drop` orphans children.** Pre-fix the
  executor's `children: HashMap<String, Child>` table stored bare
  `Child` handles; dropping the executor (e.g. UI exit while a sweep
  was still running) tore down the map without signalling anything
  and every outstanding subprocess outlived the parent. Now wraps
  each entry in `KillOnDropChild` (the round-6 RAII guard, factored
  out to `pub` for executor consumption). Same shape the
  `subprocess::run` path uses.

**Low-severity / doc cluster:**

- **L10 — CITATION.cff license + abstract regression fix.** Round-12
  L12 rewrote CITATION.cff and accidentally dropped MIT from the
  `license:` list + collapsed "dual-licensed MIT OR Apache-2.0" to
  "Apache-2.0 only" in the abstract. Restored both — Valenx remains
  dual-licensed.
- **L11 — AF2 absl `--flag=value` joined form.** Mirrors AF3's
  round-12 fix: `--fasta_paths`, `--output_dir`, `--data_dir` now
  use the separated `--flag value` form so a path containing `=`
  (legal POSIX) doesn't get mis-parsed by absl::flags. Non-path
  flags (`--max_template_date`, `--model_preset`) stay compact.
- **L12 — Stale tag-push doc claims (R8/R9 sister sweep).**
  `docs/INSTALLER.md` table row + `NEXT_PHASE.md` D-section both
  said the release CI runs "on `v*` tag push". Round-8 had removed
  the tag-push trigger; release builds are now `workflow_dispatch`
  only via `gh workflow run release.yml`. Both files updated to
  match reality.

12 new RED→GREEN tests pin every fix.

---

**Round 12 — 11 persist.rs read-cap sweep + spreadsheet parser
depth-cap completeness + voxel resolution overflow + CAM toolname
G-code injection + load_mesh JSON cap + sweep case.toml cap + 5
Python adapter shared-helper migration + TOML injection via folder
name + main→master sweep + ARCHITECTURE.md 75→141 update + CITATION
+ macOS bundle biology mention.**

**Medium-severity:**

- **M1 — 11 *persist.rs `fs::read_to_string` unbounded read.** Sister
  to round-11 R11-2 on the project-loader axis: sketch / cam / arch /
  feature-tree / techdraw / surface / spreadsheet / macro / draft /
  lattice / assembly all migrated to a `MAX_DOC_FILE_BYTES = 16 MiB`
  cap (stat + bounded `take()` for TOCTOU defence-in-depth). A
  multi-GB hostile `.ron` no longer slurps into RAM before the parser
  sees it. New `valenx_core::io_caps` module exports the canonical
  helper for any future callers.
- **M2 — Spreadsheet parser depth-cap bypass via `-` and `^` chains.**
  Round-8 wired `MAX_PARSE_DEPTH = 100` into `parse_primary`'s LParen
  arm but `parse_factor` (right-associative `^`) and `parse_power`
  (unary `-`) recursed without ever consuming a paren — `=---...x`
  and `=2^2^...^2` could blow the stack despite the cap. All three
  recursive entry points now bump + check uniformly.
- **M3 + M8 — Voxel `from_aabb` + `to_mesh_surface_nets` resolution
  overflow.** New `MAX_VOXEL_CELLS = 100_000_000` cap plus
  `checked_mul` chain on `nx*ny*nz`. New `VoxelError` enum
  (`Overflow`, `TooManyCells`) — both methods now return `Result`. All
  internal callers (`simulate::voxel_from_stock` / `animate` /
  `final_state`, `camotics_sim::Animation::{frame, frame_smooth,
  frames, frame_metadata, side_by_side, fresh_grid}`, plus the
  mesh_toolbox call site in valenx-app) migrated.
- **M4 — CAM `CommentStyle::wrap` G-code injection (machine safety).**
  `Tool.name` (user-controlled) flows into the G-code header via
  `(Tool: T1 …)` on every Fanuc-flavour post (Haas, Mach3, LinuxCNC,
  Tormach). A hostile name like `")\nG0 Z-99\n("` would otherwise
  escape the comment and inject a real motion command. `wrap` now
  scrubs `\n` / `\r` for every flavour and the close-delim character
  for the `Pair(open, close)` flavour.
- **M5 — `load_mesh` `.json` branch unbounded `read_to_string`.**
  Now uses `valenx_core::io_caps::read_capped_to_string` with
  `MAX_MESH_JSON_BYTES = 64 MiB` (meshes can be larger than ordinary
  workbench docs but 64 MiB is still well past anything an
  interactive viewer can render).
- **M6 — `sweep_selected_case` case.toml re-read uncapped.** A
  case.toml swapped between project-load and sweep-button click
  would re-read uncapped. Now uses the same cap helper at
  `MAX_PROJECT_FILE_BYTES = 1 MiB` the round-11 project loader uses.
- **M7 — 5 round-3 Python adapters + alphafold3 migrated to
  `resolve_python_binary`.** alphamissense, anndata, be-designer,
  esmfold, rfdiffusion all had the same hand-rolled
  `validate_python_binary` + `find_on_path` pattern. The `..`-traversal
  guard was inconsistent across them (some had it, some didn't).
  Migrated to the round-4 shared helper, which bundles allow-list +
  absolute-path acceptance + `..`-traversal rejection + PATH
  resolution in a single call. alphafold3 — which had the inline
  `ParentDir` check — also migrated for consistency.

**Low-severity / doc cluster:**

- **L9 — TOML injection via folder name in `render_project_toml`.**
  A folder name like `evil"\n[rbac]\ndefault_role = "viewer"\n#`
  would let an attacker who controls the scaffold destination inject
  an arbitrary `[rbac]` block into the rendered `project.toml`. New
  `sanitize_project_name` allow-list (ASCII alphanumeric + `_`,
  `.`, `-`) gates both the CLI `scaffold_project` and the GUI
  `new_case_for_adapter` paths. `render_project_toml` now returns
  `Result<String, String>`.
- **L10 — `main` → `master` branch-name sweep.** 4 missed files
  (SECURITY.md, POLICIES.md × 2 sites, STATUS.md, CHANGELOG.md
  comment block) updated to match the canonical branch name.
- **L11 — ARCHITECTURE.md `~75` / `75` → `141`.** 4 sites updated;
  added a sentence acknowledging biology as the largest surface
  (123 of 141 adapters) and expanded the domain list to include
  alignment / CRISPR / structure prediction.
- **L12 — CITATION.cff + macOS bundle `long_description` biology
  mention.** CITATION abstract rewritten to surface the bio
  ecosystem; added `bioinformatics` / `crispr` / `protein-design`
  / `genomics` / `structure-prediction` keywords. The
  `[package.metadata.bundle].long_description` on macOS now mirrors
  the `.deb extended-description` text — acknowledges 141 adapters
  + biology.

---

**Round 10 — Lattice DoS cap + print-bed traversal + 11-adapter
output-path sweep + FEM meshgen Result + pathtrace render trio Result +
load_project RBAC/audit + file-browser popup kill-switch + lazy
crashes dir + doc cluster.**

**Privacy / UX (commit `823279d`):**

- New Settings → Privacy → "Never open the system file browser"
  kill-switch suppresses every "Open in file browser" action across
  the UI (Settings → "Open crashes folder", Run → Open prepared / run
  workdir, Audit → Open audit log location, plus the four command-
  palette entries that dispatch into the same methods). On the kill-
  switch path the helper returns the path as a neutral status string
  instead of spawning `explorer.exe` / `open` / `xdg-open`. Defaults
  OFF so behaviour is preserved; Serde-default'd so existing
  `settings.json` files load unchanged.
- `install_crash_reporter` no longer pre-creates `<state_dir>/crashes/`
  on startup — `CrashReport::write_to_disk` already runs
  `create_dir_all` on a real crash.

**High-severity (commit `b71e269`):**

- **H1 — valenx-lattice `MAX_LATTICE_PLACEMENTS` cap + `checked_mul` on
  grid / on_surface / bezier placement axes (DoS class).** A hostile
  case requesting `1e9 × 1e9 × 1e9` lattice cells would overflow into a
  small positive `Vec::with_capacity` and either OOM or panic; the cap
  plus `checked_mul` rejects at validation time.
- **H2 — valenx-print-bed `export_layout` `Part.name` validation
  (path-traversal).** Pre-fix the output filename was taken verbatim
  from the user-controlled `Part.name`; a name like `../../etc/passwd`
  would land outside the workdir. Routed through
  `validate_output_basename`.
- **H3 — 11 adapters' output-path traversal closed.** Sister to the
  round-3/4/5 sweep on the `workdir.join` axis: badread, cas-offinder,
  ctffind, dssr, kallisto, nwchem, openbabel, psi4, rnastructure,
  salmon, wgsim all migrated to `validate_output_basename` /
  `validate_output_dir` for every adapter-emitted output filename.

**Medium-severity (commit `b71e269`):**

- **M4 — valenx-fem `structured_hex_mesh` + `structured_tet10_mesh` now
  fallible `Result` variant with `MAX_NODES` cap.** Closes a round-5
  sister gap — the surface-mesh generators were already capped, but
  the volumetric structured-grid generators were not. 5 internal
  callers migrated.
- **M5 — valenx-pathtrace `render_mis` + `render_bdpt` + `render_volume`
  → `Result<HdrFramebuffer, FramebufferError>`.** Closes the round-9
  sister gap — `render` was already fallible after the round-9 tracer
  pass, but the three sibling renderers still went through the
  panicking `HdrFramebuffer::new` path.
- **M6 — `load_project` now calls `rbac_check(Action::ProjectOpen)`
  and emits an `emit_audit("project.open")` line.** Closes a round-3
  RBAC sister gap; pre-fix every project load was un-audited so a
  shared-host environment had no record of which user opened which
  case bundle.

**Documentation (commit `b71e269`):**

- **L7 — `scripts/qa.sh` + `scripts/qa.ps1` `--help` text now lists
  `deny` in the `--gates` description** so contributors know the
  available gate names.
- **L8 — CHANGELOG documents Round 9 + Round 10 entries.** Round-9
  fix agent forgot to add its own entry; round-10 catches up + adds
  round-10 itself.

Round-10 ships 26 new RED→GREEN tests pinning every load-bearing fix.
(Commits `823279d` + `b71e269`.)

**Round 9 — 18-adapter sibling-field sweep + tracer Result + RBAC cap +
sync_channel deadlock + ResidualHistory VecDeque + STEP/IGES/STL caps +
precice helper migration + doc cluster.**

**H1-18 — sibling-field comprehensive sweep (18 user-data sites):**

- bio: alphafold2 (`run_script`), alphafold3 (`run_script`), blast
  (`database`), cellprofiler (`input_dir`), cromwell (`jar` relative
  form), diamond (`database`), fiji (`input_image`), foldseek
  (`database`), ilastik (`project` + `input_images[]`), kallisto
  (`transcriptome`), mmseqs2 (`target`), salmon (`transcriptome`),
  star (`genome_dir` + `reference` + `sjdb_gtf`).
- non-bio: cfd/su2 (`mesh`), chemistry/cantera (`Mechanism::External`
  path), coupling/precice (participant `case_dir`), dynamics/mujoco
  (`model.path`), fea/elmer (`mesh_dir`), mesh/gmsh
  (`Domain::MergeFile` path), md/lammps (`ReadData{path}` +
  `Eam{path}` in `stage_external_files`).
- Intentional KEEP sites (admin-managed paths) documented in-place.
  Verification grep `case.path.join(&input.*)` now shows only 9
  remaining sites, each with a round-9 classification comment.

**H19 — tracer::render → Result.** `tracer::render` migrated to
`Result<HdrFramebuffer, FramebufferError>`; the fallible
`HdrFramebuffer::try_new` propagates so a hostile scene with a
100k×100k camera returns `FramebufferError::TooLarge` instead of
panicking inside `HdrFramebuffer::new`. 6 internal callers (tests
in `tracer.rs` / `mis.rs` / `bdpt.rs` + the lib.rs doctest) migrated
with `.expect("render small framebuffer")`.

**H20 — RBAC + project.toml file-size cap.**
`MAX_RBAC_FILE_BYTES = 1 MiB` applied to `valenx_rbac::load` (new
`FileTooLarge` variant) and to the project.toml rbac-override loader
in `valenx-app::rbac_io` (1 MiB local const with `tracing::warn` on
overflow). Pre-fix a multi-GiB `rbac.json` would force the loader to
allocate the whole file before the JSON parser could reject it.

**M21 — sync_channel deadlock fix.** `ChannelLogSink` +
`ChannelProgressSink` now use `try_send` + `SinkDropCounter`
(`Arc<AtomicUsize>`) instead of blocking `send()`. One-shot lifecycle
events (Starting / Finished / Failed / Collected) still use `send()`
since they're O(1) per run and dropping would break the UI state
machine.

**M22 — ResidualHistory Vec → VecDeque.** O(1) `pop_front` eviction
replaces O(n) `Vec::remove(0)` at cap. Same migration for the
`by_field_log10` mirror. `show()` pulls `PlotPoints` via
`from_iter(samples.iter().copied())`.

**M23 — STEP / IGES / AP242 file-size cap.**
`MAX_CAD_INTERCHANGE_FILE_BYTES = 256 MiB` gates `fs::read_to_string`
in `step::read`, `iges::read`, `iges_trimmed::read`, and
`ap242::{read_metadata, append_metadata, count_solids}`. New
`FileTooLarge` variant.

**M24 — STL file-size cap.** `MAX_STL_FILE_BYTES = 512 MiB` (legit
STL exports can hit hundreds of MiB) gates `valenx_viz::stl::load`.

**M25 — precice copy_dir_recursive helper migration.** Extracted the
round-6 hardened `copy_dir_recursive` (symlink rejection + depth cap
+ per-file size cap) from `valenx-addons` into
`valenx-core::adapter_helpers` (public). The precice adapter dropped
its local copy in favour of the shared helper.

**Doc cluster (D26-28):** CHANGELOG gains Round 7 (review-only) +
Round 8 entries. `RELEASING.md` Hotfix + "Tag was wrong" sections
rewritten around `gh workflow run release.yml -f tag=...` instead of
`git tag + push`. `docs/QA.md` adds step 7 (`cargo deny check`).
`release.yml` header comment trimmed of the obsolete "Builds on
every tag push (v*)" claim.

Round-9 ships 22 new RED→GREEN tests pinning every load-bearing
fix. (Commit `da319d3`.)

**Round 8 — 8-adapter sibling-field confined_join sweep + CIGAR op-count
cap + 7 mediums + 4 lows.** Round 8 finds eight more adapters where a
sibling field of an already-confined-join'd field was still using bare
`case.path.join` (bcftools / cromwell / planemo / pksim / rosetta /
snakemake / nextflow — nextflow's `inputs` got the `confined_join`
treatment for the first time). `CIGAR` parser gets a second cap
(`MAX_CIGAR_OPS = 1_000_000`) sibling to the round-6
`MAX_CIGAR_OP_LEN` so a 4 M-op CIGAR string with each op = 1 still
exits cleanly. `capture_subprocess_stdout` (the helper every probe
goes through) gets `MAX_PROBE_OUTPUT_BYTES = 1 MiB` bounded read. CAM
adapter's `Duration::from_secs_f64(tool_life_hours * 3600.0)` now
clamps to `Duration::MAX` (no more panic when a hostile case sets
tool life to f64::INFINITY). Spreadsheet parser gets
`MAX_PARSE_DEPTH = 100` (sister to the round-3 evaluator cap). FEM
dynamics gets `MAX_TIME_STEPS = 1_000_000` cap. Palette cache clone
gated on `palette.open` so the per-frame `clone()` doesn't churn at
~3 MB/s when the palette is closed. ResidualHistory gets
`MAX_RESIDUAL_SAMPLES_PER_FIELD = 50_000` with FIFO drop. RunEvent +
SweepEvent channels switched to `sync_channel(4096)` back-pressure
(sister to the round-6 subprocess fix). `Solid::rotated` zero-axis
`is_finite` check (sister to round-3 mirrored). Three sister
state-loaders (settings / history / first-run) get
`MAX_STATE_FILE_BYTES = 10 MiB` cap. `scripts/qa.sh` + `qa.ps1` now
run `cargo deny check` (matches SECURITY.md claim). `RELEASING.md`
tag-push procedure replaced with `workflow_dispatch` flow.
`STATUS.md` OCCT contradiction resolved. `CHANGELOG` round-6 test
count fixed. `.deb` `extended-description` ambiguous `141+123`
phrasing fixed. Round-8 ships 16 new RED→GREEN tests pinning every
load-bearing fix. (Commit `0a45c05`.)

**Round 7 — review pass.** Round 7 was a review-only pass — no fix
commit landed. Findings flowed into the round-8 fix queue.

**Round 6 — CIGAR DoS cap + 6 aligners confined_join + audit log fs2 lock +
glTF accessor cap + 8 mediums + 7-file doc sweep.** Round-6 sweeps the
remaining DoS / traversal classes the round-1..5 work either missed or
introduced.

**Critical (4):**

- **CIGAR amplification DoS** (`valenx-genomics::format::sam::Cigar::parse`).
  Pre-fix a 100-byte SAM with `CIGAR=4294967295M` triggered ~4.3 B
  `BTreeMap` insertions in `build_pileup`. New
  `MAX_CIGAR_OP_LEN: u32 = 1_000_000` rejects oversized ops at parse
  time; defensive `MAX_PILEUP_SPAN: usize = 100_000_000` guard in
  `build_pileup` catches hand-rolled `SamRecord` values that bypass
  the parser.
- **Reads-loop traversal class** (6 aligner adapters: bowtie2, hisat2,
  kallisto, minimap2, salmon, STAR). The pre-round-6 `case.path.join(read)`
  fall-through accepted absolute paths and `..` traversal; replaced with
  the BWA-style `confined_join(&case.path, read)?` pattern. Class-sweep
  grep verifies `read.is_absolute()` is gone from
  `crates/valenx-adapters/bio/`.
- **Audit log race condition** (`valenx-audit::AuditWriter::append`).
  Two concurrent appends could read the same tail, compute the same
  prev_hash, and silently fork the chain. New `fs2`-backed advisory
  lock on `<log>.lock` serialises the read-tail/compute-hash/write
  critical section across threads AND processes. `rotate_if_needed`
  gets the same lock; internal `_locked` variant for callers that
  already hold it.
- **glTF accessor count unbounded** (`valenx-occt-exchange::gltf2_reader`).
  A manifest with `"count": 18446744073709551615` would let
  `Vec::with_capacity(count as usize)` ask the allocator for tens
  of exabytes. New `MAX_GLTF_ACCESSOR_COUNT: usize = 1_000_000`
  rejects in `accessor_bytes` before any allocation.

**Medium (8):**

- **`kill_on_drop` wiring** — pre-fix the field on `PreparedJob` was
  set across 140+ adapter sites but no path actually honoured it.
  New `KillOnDropChild` RAII wrapper in `valenx-core::subprocess`
  SIGKILLs the child on Drop when the flag is set. The executor
  path stays advisory (it already kills explicitly via the
  submit/poll/cancel lifecycle).
- **BackgroundRun cancel** (`valenx-app::genetics::rna_designer`).
  RNA Designer's inverse-design / LinearDesign / λ-sweep workers
  now carry an `Arc<AtomicBool>` cancel flag wired through the
  App's `on_exit`; future workers that poll the flag exit early
  on window close (the v1 workers don't yet poll inside their
  `valenx-rnadesign` call tree — documented limitation).
- **Subprocess unbounded mpsc** (`valenx-core::subprocess`).
  Channel between the stdout/stderr pump threads and the main loop
  is now `mpsc::sync_channel(SUBPROCESS_CHANNEL_CAPACITY = 4096)`
  so a chatty solver back-pressures the child rather than OOMing
  the parent. Also adds `MAX_LINE_BYTES = 1 MiB` cap via a
  `read_capped_line` helper that truncates pathological lines.
- **CrashReport size cap** (`valenx-crash-reporter::load_all`).
  A 20 MiB `*.json` planted in the crash directory used to slurp
  in full via `std::fs::read`. New `MAX_REPORT_BYTES = 10 MiB` cap
  via stat-then-take guards the read path; `read_capped` helper
  bounds memory growth even when metadata is unreliable.
- **copy_dir symlink-naive** (`valenx-addons::install`).
  Recursive copy now uses `fs::symlink_metadata`, refuses every
  symlink encountered (top-level + per-entry), caps per-file at
  64 MiB, and limits recursion to 32 levels. Manifest reader
  (`read_manifest_at`) gets a matching `MAX_MANIFEST_BYTES = 1 MiB`
  cap.
- **Solid::translated / Solid::rotated finite validation**
  (`valenx-cad::solid`). Both methods now return `Result<Self, CadError>`
  and reject NaN/inf components (otherwise the BRep matrix silently
  inherited the bad values and downstream tessellation panicked).
  All 46 callers of `translated` and 16 of `rotated` migrated;
  `valenx-py::cad::Solid::translated` now surfaces `PyValueError`.
- **HdrFramebuffer checked_mul** (`valenx-pathtrace::framebuffer`).
  `HdrFramebuffer::new(65536, 65536)` used to silently wrap the
  `width * height` multiplication. New `try_new` returns
  `FramebufferError::TooLarge` for overflow or > 8K² pixels;
  `new` is now a panicking thin wrapper around `try_new`.
- **EntryPoint::Python.module validation**
  (`valenx-addons::manifest`). A manifest `module = "evil; __import__('os')"`
  could pass straight to `importlib` and execute. New
  `validate_python_module` allow-list restricts the field to
  `[a-zA-Z0-9_.]`.

**Doc cluster (7-file sweep):**

- README + STATUS bio-adapter count corrected from 121 → 123
  (matches `ls crates/valenx-adapters/bio/ | wc -l`).
- POLICIES.md "Tier 1 tested in CI every commit" claim rewritten
  to acknowledge the manual-trigger-only CI state since 2026-05-24,
  with the rationale link to `docs/CI.md`.
- SECURITY.md `cargo audit` / `cargo deny` "block merges in CI"
  claim rewritten — cargo-audit was removed, CI is manual-trigger.
  Now: "Locally, run `bash scripts/qa.sh` which includes
  `cargo deny check`."
- RELEASING.md `tag push fires release.yml` and
  `cargo wix init --force clobbers template` sections deleted —
  both untrue post-round-1/2 (release.yml is manual-trigger and
  wix template lives at the standard cargo-wix path).
- CONTRIBUTING.md "pre-alpha note: workspace hasn't been scaffolded
  yet" deleted (years stale); `main` → `master` fixed at the two
  remaining mentions.
- TESTING.md / CONTRIBUTING.md / LANGUAGES.md non-existent crate
  references: `valenx-cfd` → `valenx-cfd-native`, `valenx-bench`
  → `valenx-cfd-native` (the actual benchmark home), `valenx-fea`
  → `valenx-fem`.
- `valenx-app/Cargo.toml`'s `.deb` `extended-description` rewritten
  to acknowledge the 141 live adapters + the biology surface
  (was claiming only 17 OSS solvers).

Round-6 ships 20 new RED→GREEN tests pinning every load-bearing fix.

**Round 5 — Tier-1 DoS class sweep + R-script allow-list + Windows /
macOS LICENSE shipping.** Round 5 sweeps the DoS class to
`valenx-fields` (FRD / VTK-legacy / VTU / KdTree readers got the
same `MAX_*_LIST_LEN` cap + bounds-checked iters + monotonic-offset
validation + iterative kd-tree descend treatment as the round-1
PLY / round-3 DCD / round-4 IGES readers). R-script
arbitrary-binary-exec sister-class fix lands `validate_rscript_binary`
helper in `valenx-core::adapter_helpers` and wires it into the icodon
and seurat adapters (the two R-driven cases in v1). Windows `.msi`
+ macOS `.app` bundle now ship `LICENSE-MIT` + `LICENSE-APACHE`
per MIT §1 + Apache §4(a) compliance. Python allow-list extended
for `python3.X.exe` + `python3X.exe` (Windows conda / pyenv-win
regression). Doc sweep completed across TESTING.md /
CONTRIBUTING.md sub-blocks / LANGUAGES.md. `output_dir` validation
landed in cwltool, kallisto, salmon. cwltool `inputs` field wrapped
with `confined_join`. `structured_box_mesh` panicking variant
deleted, all 45 internal callers migrated to fallible
`try_structured_box_mesh` Result variant. OBJ fan-triangulation
`MAX_OBJ_FACE_VERTICES = 4096` cap. `confined_join` Windows
drive-relative pattern detection + symlink-escape limitation
documented. 38-adapter broken Python error template swept. mRNA
designer codon-table sum claim + 5′UTR provenance comment made
honest. Round-5 ships 22 new RED→GREEN tests pinning every
load-bearing fix.

**Round 4 — mechanical class-pattern sweeps + IGES DoS + spawn_sweep parity
+ doc cluster.** Previous rounds fixed sample adapters but the same vulnerable
patterns repeated across 100+ adapter crates. Round 4 grep-finds every
occurrence of each pattern class and applies the same fix systematically.

**Sweeps (mechanical, across 100+ files):**

- **Sweep 1 — Python interpreter allow-list (40 adapters).** Round-3
  added `validate_python_binary()` but only landed in 8 adapters; the
  vulnerable `if input.python == "python3" || input.python == "python"`
  block survived in 40 more. Round 4 introduces
  `valenx_core::adapter_helpers::resolve_python_binary()` (validate +
  resolve in one call) and the sweep wires it into every affected
  adapter: alphafold2, be-hive, biopython, cellprofiler, chopchop,
  chroma, crispor, crispritz, deepchem, dnachisel, esm3, esmc, esm-if,
  eternafold, forecast, hoomd, indelphi, mdanalysis, mdtraj, molstar,
  msprime, ngl, nupack, omegafold, openfold, openmm, pegfinder,
  primedesign, prody, proteinmpnn, pydna, pyrosetta, pysbol, rdkit,
  rfantibody, rosettafold, scanpy, scvi, tskit, occt.
  Verification: `grep -rln 'input.python ==' crates/valenx-adapters/`
  is empty.

- **Sweep 2 — `confined_join` for staged-file inputs (74 sites across
  64 files).** Round-2/round-3 fixed 7 adapters; round 4 sweeps every
  remaining call site. The script classifies each `case.path.join(&input.X)`
  by field name: STAGED fields (the adapter reads + copies the file
  into the workdir — `input`, `query`, `script`, `xml`, `reference`,
  `topology`, …) get `confined_join` which rejects absolute paths and
  `..` traversal; SYSTEM-PATH fields (the adapter reads the file in
  place from a user-managed install — `jar`, `model_dir`, `db_dir`,
  `fiji_app`, …) deliberately keep `case.path.join` since they
  legitimately need to be absolute. Affected adapters span bio
  (amber-sander, art, autodock4, badread, beast2, blast, bowtie2,
  cas-offinder, cello, cellprofiler, clustalo, copasi, cpptraj,
  cromwell, ctffind, curves, cwltool, deepvariant, diamond, dssr,
  eman2, fasttree, fiji, foldseek, hisat2, igv, iqtree, j5, jalview,
  lineardesign, linearfold, mafft, mcell, mfold, minimap2, mmseqs2,
  mrbayes, muscle, namd, nwchem, omegafold, openbabel, physicell,
  planemo, plumed, psi4, raxml-ng, rnastructure, rosetta, samtools,
  simrna, slim, smoldyn, snakemake, tcoffee, wgsim, x3dna, xtb),
  cad (freecad), cfd (su2), coupling (precice), fea (code-aster,
  openradioss), and mesh (netgen). Verification: post-sweep grep
  shows only the 14 SYSTEM-path field names remain.

- **Sweep 3 — `extra_args` ordering (8 adapters).** Pre-fix, BWA,
  beast2, cwltool, hmmer, mafft, minimap2, nwchem, and slim all
  pushed `extra_args` BEFORE the positional source paths, which
  meant a hostile case.toml could supply `extra_args = ["--help"]`
  and shadow the source path positional. The sweep moves the
  `for arg in &input.extra_args` loop to AFTER every positional in
  each affected `prepare()`. nwchem also gets a refactor to apply
  the rule to both the MPI and non-MPI branches.

- **Sweep 4 — `validate_output_basename` (62 adapters).** Round-3
  fixed 4 (bionetgen, iqtree, art, fasttree); round 4 sweeps every
  remaining adapter that takes an `output_basename` / `prefix` field
  in case.toml. The injection point lands right after the
  `let input = ...::from_case_dir(&case.path)?;` line in `prepare()`.
  Affected: alphamissense, amber-sander, anndata, be-designer, be-hive,
  cello, cellprofiler, chopchop, chroma, clustalo, crispor, crispritz,
  cromwell, curves, dnachisel, eman2, esm-if, esm3, esmc, eternafold,
  fiji, foldseek, forecast, hoomd, icodon, ilastik, indelphi, j5,
  jalview, lineardesign, linearfold, mdtraj, mfold, molstar, msprime,
  ngl, nupack, omegafold, pegfinder, pksim, planemo, plumed,
  primedesign, prody, proteinmpnn, pydna, pyrosetta, pysbol, raxml-ng
  (`prefix`), relion, rfantibody, rfdiffusion, rosetta, rosettafold,
  seurat, simrna, slim, tcoffee, tskit, x3dna, occt, elmer.

**Critical individual fixes:**

- **IGES DoS cluster (4 sites in iges.rs + NURBS surface math in
  iges_trimmed.rs).** Add `MAX_IGES_LIST_LEN = 1_000_000` constant
  + `StepIgesError::ListTooLarge { count, max }` + `ArithmeticOverflow`
  variant. Each `Vec::with_capacity(count)` call site now bounds-checks
  count first; NURBS Type 128 surfaces use `checked_add` / `checked_mul`
  on the (k1 + m1 + 2) and (k1+1) * (k2+1) arithmetic. RED→GREEN
  test crafts an IGES file with `count = usize::MAX`, confirms
  reader returns `ListTooLarge` instead of OOM.

- **spawn_sweep catch_unwind parity** — round 3 wrapped `spawn_inner`'s
  thread body in `std::panic::catch_unwind`; `spawn_sweep` was
  missing the equivalent guard so a panicking sweep would unwind
  the worker silently and leave the progress bar stuck forever.
  Round 4 mirrors the spawn_inner pattern + adds a RED→GREEN test
  that panics from a sweep adapter's `prepare()` and verifies the
  outer `catch_unwind` surfaces evidence in the `SweepEvent` stream.

- **DCD natom signed check** — pre-fix, `i32_le(&natom_rec[0..4])? as usize`
  silently wrapped negative values to ~18 quintillion, then
  `Vec::with_capacity(natom)` OOMed. Round 4 adds `MAX_NATOM = 25_000_000`
  (about 100 MiB of f64 coordinates) + an explicit signed check
  before the cast. Two RED→GREEN tests pin both the negative and
  positive-overflow paths.

- **Crash report atomic write** — `CrashReport::write_to_disk` now
  goes through an inline `atomic_write_bytes` helper (sidecar `.tmp`
  + rename) rather than `std::fs::write`. Crash reports are written
  from the panic hook, so a double-panic mid-write would otherwise
  leave a truncated JSON behind for `load_all` to silently drop.

**Medium fixes:**

- **valenx-optimize cartesian-product overflow** — `GridOptimizer::plan`
  now uses `try_fold` + `checked_mul` for the cell-count product and
  caps at `MAX_GRID_CELLS = 10_000_000` (with a structured
  `InvalidConfig` error rather than a silent wrap or OOM).
- **valenx-arch Column segments cap** — `MAX_COLUMN_SEGMENTS = 4096`
  cap on `ColumnSection::Circular { segments }` to stop a malicious
  BCF from driving `Vec::with_capacity(u32::MAX as usize)` to OOM.
- **valenx-geo `try_regular_polygon_xy`** — new fallible variant
  with `MAX_REGULAR_POLYGON_SIDES = 100_000` cap; the existing
  panicking `regular_polygon_xy` is preserved for call sites that
  know their input is bounded.
- **AeroRunHandle::cancel()** — cancellation token wired through the
  aero sweep loop with per-angle checkpoints; App `on_exit` now also
  cancels `aero.run` (not just `run_handle` / `sweep_handle`).
- **`pump_sweep_events` per-frame cap** — adds the same
  `MAX_EVENTS_PER_FRAME = 256` budget as `pump_run_events` so a
  chatty 10,000-case sweep doesn't starve the UI thread.
- **opener crate promoted to `[workspace.dependencies]`** — was
  pinned in `valenx-app/Cargo.toml` directly.
- **docs/CI.md crate count** — 235 → 249 (workspace has grown).
- **vendor/vtkio/CHANGES_VS_UPSTREAM.md** — adds an explicit
  verification command that proves the zero-source-diff invariant
  vs upstream 0.6.3.

**Documentation sweep:**

- `cargo test --workspace` references in CONTRIBUTING.md (×2),
  QUICKSTART.md, STATUS.md, PR template, tests/README.md, and
  `.github/workflows/ci.yml` all replaced with the canonical
  `bash scripts/qa.sh` / `pwsh scripts/qa.ps1` invocation
  (the scoped harness). The forbidden-pattern reference in
  `docs/QA.md` and the historical narrative in CHANGELOG / STATUS
  / plans is intentionally untouched.
- `Valenx.app.zip` references in RELEASING.md (×3) replaced with
  the `.dmg` artefact name that release.yml actually produces.
- `docs/QA.md` self-contradiction fixed — was claiming "no active
  GitHub Actions workflow" while ci.yml / ci-nightly.yml /
  release.yml exist as `workflow_dispatch`-gated workflows.
  Now acknowledges the manual-only workflows + the fact that
  the active CI invokes the scoped harness.
- Dead docstring in `valenx-app/src/lib.rs` for `first_run_open`
  updated to match the runtime behaviour (wizard never auto-opens
  on first launch — that branch was deliberately disabled).

Round-3 code-review hardening pass with explicit RED→GREEN
test-first verification (no fix lands without a failing test the
fix turns green). Round-1 and round-2 fixes addressed `cargo
deny check licenses` configuration, BCF path traversal, PLY DOS,
RBAC fail-closed, crash-reporter UTF-8 mojibake, `confined_join`
across seven bio adapters, SAM seq/qual + VCF column-count
validators (round-2 reader wiring confirmed by round-3 tests),
PDB 66-char minimum, atomic `state_paths::atomic_write`,
adaptive_clearing NaN-finite gates, AlphaFold3 path arg
splitting, RELION MPI procs cap, and FEM modal solver size
cap. Round-3 adds:

- **valenx-macro Python-export injection** — sanitise `\n`/`\r` in
  `Macro::name`, `Macro::description`, and `Step.label` before they
  hit `format!(# Macro: {} ...)`. Without the sanitisation a newline
  in any of these fields would break out of the comment context and
  let an attacker land an executable Python statement at column 0
  of the generated script. Also escape `"""` in multi-line
  `run_python` script bodies so a user-supplied script can't break
  out of the triple-quoted raw-string literal.
- **valenx-app worker panic** — wrap the run-pipeline thread body
  in `std::panic::catch_unwind`. Previously a panicking adapter
  unwound the worker silently, leaving the UI stuck on "Starting…"
  forever with no `Failed` event ever delivered.
- **valenx-app `on_exit`** — cooperatively cancel `run_handle` /
  `sweep_handle` when the window closes, so an in-flight
  subprocess (e.g. `simpleFoam.exe`) doesn't get orphaned. Adds a
  200 ms grace pause for the worker to honour the cancellation
  token before eframe tears down the GL context.
- **valenx-i18n UTF-8 mojibake** — rewrite `format_placeholders` to
  iterate codepoints rather than bytes. Pre-fix, every byte of a
  multibyte UTF-8 codepoint (`Größe`, `日本語`, 🎉) became its own
  Latin-1 `char`. Same fix applied to `valenx-openscad-import`,
  `valenx-genomics::gff3_decode`, and
  `valenx-occt-exchange::renumber_entities`.
- **valenx-rbac silent demote** — change `default_role` from `Role`
  (serde-default `Runner`) to `Option<Role>`. A project override
  that only adds a per-user mapping no longer silently demotes a
  hardened global `Viewer` to `Runner` just because the on-wire
  shape made "absent" and "Runner" indistinguishable.
- **valenx-bio DCD DOS cap** — refuse to allocate when the DCD
  header declares `nframes = i32::MAX` or a single Fortran record
  claims a 4 GB length. New `DcdError::TooLarge` short-circuits
  before any `Vec::with_capacity` / `vec![0u8; N]` call.
- **PhysiCell adapter** — restrict `[bio.physicell].binary` to a
  relative path inside the case directory. Pre-fix a hostile
  case.toml could set `binary = "/usr/bin/curl"` and turn "Run
  case" into arbitrary code execution.
- **AlphaFold3 + AlphaMissense + AnnData + Be-Designer + ESMFold +
  RFdiffusion adapters** — allow-list `input.python` against
  `valenx_core::adapter_helpers::ALLOWED_PYTHON_NAMES`. Same threat
  shape as PhysiCell — `python = "/usr/bin/curl"` no longer turns
  the case run into arbitrary exec.
- **bionetgen / iqtree / art / fasttree adapters** — validate
  `output_basename` / `prefix` / `output_prefix` / `output` as
  single-component path basenames via the new
  `valenx_core::adapter_helpers::validate_output_basename` helper.
  Pre-fix a value like `"../../etc/cron.d/x"` could let the
  subprocess write outside the workdir.
- **wgsim / fasttree / mrbayes positional shift** — move
  `extra_args` to be appended AFTER the positional arguments, so a
  hostile case.toml can't slip a phantom positional in via
  `extra_args = ["phantom"]` and shift the real positionals onto
  different argument slots.
- **valenx-cam G-code NaN/inf** — promote the `feed > 0.0` check
  in `format_g1` / `format_g1_5ax` to `feed.is_finite() && > 0.0`,
  add finite checks on the XYZ target. Pre-fix `+inf` would
  produce `F inf` which most controllers parse as "unlimited
  rapid" (and smash the spindle). `format_g0` substitutes 0.0 for
  any non-finite axis as a fault-tolerant fallback.
- **valenx-cam n_passes infinite loop** — cap
  `(depth / step_down).ceil() as usize` at 10 000 via the shared
  `op::compute_n_passes` helper. Pre-fix
  `step_down = f64::MIN_POSITIVE` would let `n_passes` saturate to
  `usize::MAX` and the per-pass loop run ~2^64 times. Applied
  uniformly across all 10 Z-stepping cam ops.
- **valenx-mesh PLY signed→u32** — replace `value as u32` with the
  new `checked_index` helper that rejects negatives / NaN / ±inf /
  overflow. Pre-fix a negative index in a signed-int face list
  silently saturated to 0, producing wrong-but-not-out-of-bounds
  connectivity that pointed at random vertices. New
  `PlyError::InvalidIndex` variant.
- **valenx-fem panic** — replace the `assert!()` calls in
  `structured_box_mesh` with a fallible `try_structured_box_mesh`
  returning `NativeSolverError::InvalidParams`. Convenience
  non-fallible wrapper retains the original panicking signature
  for tests with hard-coded inputs. Also adds `checked_mul` on
  `(nx+1)*(ny+1)*(nz+1)` plus a 200 M-node cap.
- **valenx-cad mirror panic** — `Solid::mirrored` was a hard
  `assert!(n_len > 0.0, ...)` on the normal vector. Round-3 changes
  the signature to `Result<Self, CadError>` and returns
  `InvalidParam` for zero, NaN, or ±inf components. Feature-tree
  evaluators map the error onto `FeatureError::CadError`.
- **valenx-spreadsheet recursion cap** — add
  `MAX_EVAL_DEPTH = 100` to `Evaluator::evaluate_expr` so a
  pathologically nested expression (`1*(1*(1*(...)))`) can no
  longer blow the OS stack before the per-cell cycle detector
  fires. Adds a `depth: usize` field on the `Evaluator` struct.

### Changed
- **deny.toml license allow-list** — dropped four unused licenses
  (MIT-0, MPL-2.0, Unicode-DFS-2016, OFL-1.1 — the last is still
  accepted for epaint via the per-crate exception block). Brings
  `cargo deny check licenses` to a warning-free pass.
- **License:** moved from the custom "Valenx Software License
  Agreement" (proprietary) to **dual MIT OR Apache-2.0** — the Rust
  ecosystem standard, matching rustc, Tokio, Serde, Bevy, and most
  of the Rust crate ecosystem. `LICENSE.txt` removed; `LICENSE-MIT`
  and `LICENSE-APACHE` added at the repo root. The SPDX expression
  `MIT OR Apache-2.0` applies workspace-wide via
  `[workspace.package].license`.
- **README:** rewrote the opening section to be public-facing —
  replaced the phase-status wall with a project overview, install
  instructions, a slimmed "Supported solvers" grouping, and a clean
  dual-license note. Existing capability tables remain accessible
  via [STATUS.md](./STATUS.md).
- **OS metadata:** updated cargo-deb assets, cargo-generate-rpm
  assets + `license` field, cargo-bundle copyright string, WiX
  installer comments, and the `valenx-i18n` About-dialog locale
  string to reference the new dual-license layout.
- **Internal docs:** `POLICIES.md`, `CONTRIBUTING.md`,
  `CITATION.cff`, and the `valenx-icons` module docstring updated
  to reflect MIT OR Apache-2.0.

---

**🎯 Phase 44.5 RNA folding expansion (3 bio adapters across 1
phase, sister-expansion of Phase 28).** Phase 44.5 (RNA folding
expansion — mfold + EternaFold + LinearFold) ships 3 more bio
adapters as a sister-expansion of the Phase 28 ViennaRNA /
RNAstructure / NUPACK trio, taking the bio-adapter count to
**111** and the headline live-adapter total to **131**. mfold/
UNAFold is Michael Zuker's classic Zuker / Stiegler dynamic-
programming RNA folder (academic-license, single-binary CLI
sister to Phase 18 BWA with `KEY=VALUE`-style `mfold SEQ=
<sequence> NA=RNA T=<temperature>` invocation, surfaces an
`"academic"` / `"non-commercial"`-keyworded license-awareness
warning sister to ViennaRNA / NUPACK); EternaFold is the Eterna
project's MIT-licensed ML-aware folder via the Das lab's `arnie`
Python wrapper (Python-script subprocess sister to Phase 17
Biopython / Phase 28 NUPACK); LinearFold is Baidu / Oregon
State's Apache-2.0 beam-search linear-time folder (the folding-
only sister to Phase 43 LinearDesign — same beam-search core,
same group, applied to the inverse problem; non-standard stdin
contract — sequence on stdin, structure on stdout).

**🎯 Phase 35.5 base + prime editing design (4 bio adapters
across 1 phase, sister-expansion of Phase 35).** Phase 35.5
(base + prime editing design — BE-Designer + BE-Hive +
PrimeDesign + pegFinder) ships 4 more bio adapters as a sister-
expansion of the Phase 35 CHOPCHOP / CRISPOR / Cas-OFFinder
trio with the modern non-cleavage editing tools the Phase 35
Cas9-cut-focused adapters don't cover, taking the bio-adapter
count to **115** and the headline live-adapter total to **135**.
BE-Designer is the Komor lab's MIT-licensed base-editor guide
design tool; BE-Hive is the Liu lab's MIT-licensed base-editing
outcome predictor (the canonical Python module is `be_predict`);
PrimeDesign is the Liu lab's MIT-licensed pegRNA designer for
the Anzalone / Liu prime-editing system; pegFinder is the
Komor lab's MIT-licensed alternative pegRNA finder with a
different scoring model emphasising pegRNA secondary-structure
stability + RT-template-length tradeoffs. All four ride the
established Python-script subprocess pattern (sister to Phase
17 Biopython, Phase 35 CHOPCHOP / CRISPOR, Phase 41 pydna,
Phase 43 DNA Chisel / iCodon).

**🎯 Phase 35.6 edit-outcome prediction (4 bio adapters across
1 phase, sister-expansion of Phase 35 + 35.5).** Phase 35.6
(edit-outcome prediction — inDelphi + FORECasT + AlphaMissense
+ CRISPRitz) ships 4 more bio adapters that close the design →
predict-outcome → off-target loop, taking the bio-adapter count
to **119** and the headline live-adapter total to **139**.
inDelphi is the Liu lab's MIT-licensed Cas9-cut indel pattern
predictor; FORECasT is the Sanger Institute's Apache-2.0
alternative indel predictor (the Python module is `selftarget`,
named after Allen's data-collection assay rather than the
predictor's published name); AlphaMissense is DeepMind's
missense-effect predictor extending the AlphaFold lineage —
released under **CC-BY-NC-SA-4.0 / academic non-commercial
weights** with a mandatory probe warning surfaced whenever
Python is on PATH (sister to AlphaFold 3 / mfold / ViennaRNA /
NUPACK / VMD / NAMD); CRISPRitz is the Pinello lab's MIT-
licensed variant-aware off-target genome-wide search (sister to
Phase 35 Cas-OFFinder with the distinguishing property of
walking population VCFs for off-target sites in specific
haplotypes).

**🎯 Phase 45 pharmacokinetics + RNA tertiary structure (2 bio
adapters across 1 phase, opens 2 new domains).** Phase 45 (PK-
Sim + SimRNA) ships 2 more bio adapters opening **two new
domains** in Valenx — the **first PK/PD pharmacokinetics modeling
category** (PK-Sim) and the **first RNA tertiary 3D structure
prediction category** (SimRNA) — taking the bio-adapter count to
**121** and the headline live-adapter total to **141**. PK-Sim
is the Open Systems Pharmacology suite's GPL-2.0 physiologically-
based PK (PBPK) simulator, the de-facto open-source PBPK modeling
tool descended from the Bayer internal pharmacokinetic simulator
opened to the community via the OSP Initiative; models whole-body
drug ADME using a physiologically-grounded compartmental
representation. SimRNA is the Bujnicki group's GPL-3.0 coarse-
grained Monte Carlo RNA tertiary-structure predictor — the only
adapter in Valenx that predicts the full 3D Cartesian backbone of
an RNA (the Phase 28 + 44.5 folders predict 2D secondary
structure only); five-bead per-nucleotide coarse-graining +
replica-exchange Monte Carlo. Both adapters ride the single-
binary subprocess shape sister to Phase 18 BWA / Phase 32.5
Smoldyn / Phase 5 GROMACS / Phase 43 LinearDesign.

**🎯 Bio ecosystem complete + Phase 43 mRNA design (108 bio
adapters across 39 phases).** Phase 43 (mRNA design — DNA Chisel
+ LinearDesign + iCodon) ships 3 more bio adapters on top of the
bio-ecosystem-complete milestone reached at Phase 22.5 + 42,
**opening the first mRNA / vaccine therapeutic design domain in
Valenx** — the codon-optimization + joint-design half of the
mRNA workflow that the existing Phase 28 RNA structure-prediction
stack and Phase 33 / 41 synthetic-biology composition stack leave
incomplete. Phase 43 takes the bio-adapter count to **108** and
the headline live-adapter total to **127**. DNA Chisel is the
Edinburgh Genome Foundry's MIT-licensed constraint-driven Python
codon optimizer (Python-script subprocess sister to Phase 17
Biopython / Phase 41 pydna / Phase 42 Mol* / NGL); LinearDesign
is Baidu Research's Apache-2.0 single-binary CLI that jointly
optimizes codon usage and mRNA secondary-structure stability —
the modern mRNA-vaccine design workhorse since the 2021 _Nature_
paper; iCodon is the Vejnar lab's GPL-3.0 R-based codon-level
mRNA stability predictor (Rscript subprocess sister to Phase
19.6 Seurat).

**🎯 Bio ecosystem complete — every category from the original
/review list now covered (105 bio adapters across 38 phases).**
Phase 22.5 (workflow expansion — planemo + Cromwell + cwltool) and
Phase 42 (web visualization — Mol* + NGL Viewer) ship 5 more bio
adapters on top of the prior set, taking the bio-adapter count to
**105** and the headline live-adapter total to **124**. The bio
surface spans alignment / sequence editors / cheminformatics /
cryo-EM / CRISPR / DNA geometry / docking / MD analysis / MD
engines / microscopy / phylogenetics / population genetics /
protein design / quantum chemistry / RNA structure / sequence read
simulators / single-cell / spatial stochastic / structure
prediction / structure search / synthetic biology / systems
biology / variant calling / viewers (desktop + web) / web
visualization / workflow managers. Phase 22.5 sister-expands Phase
22 workflow managers with three more workflow languages (Galaxy /
WDL / CWL); Phase 42 opens the first modern web 3D molecular
visualization category in Valenx.

**🎉 100 bio adapters across 36 biology / biotech / chemistry
phases.** Phase 32.5 (spatial stochastic — Smoldyn + MCell), Phase
40 (microscopy — Fiji + CellProfiler + Ilastik), and Phase 41
(sequence editors — pydna + Jalview) ship 7 more bio adapters on
top of the prior set, taking the bio-adapter count to **100** and
the headline live-adapter total to **119**. Phase 32.5 sister-
expands Phase 32 systems biology with cell-scale spatial
stochastic simulators; Phase 40 opens the first microscopy /
bioimage analysis category in Valenx; Phase 41 opens the first
plasmid-design / alignment-viewer category.

### Added

- **DNA Chisel + LinearDesign + iCodon adapters (Phase 43).**
  Opens the **first mRNA / vaccine therapeutic design domain**
  in Valenx — the codon-optimization + joint-design half of the
  mRNA workflow that the existing Phase 28 RNA structure-
  prediction stack and Phase 33 / 41 synthetic-biology
  composition stack leave incomplete. Three new adapters land
  under `crates/valenx-adapters/bio/` spanning the codon-
  optimization + joint-design tradeoff space:
  `valenx-adapter-dnachisel` (Edinburgh Genome Foundry's
  constraint-driven codon-optimization library — MIT, version
  range `3.0.0..4.0.0`, ribbon `bio.dnachisel.optimize`; the de-
  facto Python choice for end-to-end synthetic-gene design
  pipelines feeding into Phase 33 j5 assembly + Phase 41 pydna
  cloning workflows; Python-script subprocess shape sister to
  Phase 17 Biopython / Phase 19.5 Scanpy / Phase 33 pySBOL /
  Phase 41 pydna / Phase 42 Mol* / NGL Viewer; knobs `script`
  (`.py` enforced) / `python` (default `"python3"`) /
  `input_fasta` (`Option<PathBuf>` — optional starting `.fa` /
  `.fasta` FASTA) / `output_basename`; `prepare()` enforces
  `.py`, routes script + optional input_fasta through
  `confined_join` to stage them safely in the workdir, writes
  `valenx_params.json` with `output_basename` always plus
  `input_fasta` (staged filename) only when set — key omitted
  entirely when `None` rather than emitted as `null`, matching
  the hand-rolled JSON convention the rest of the bio adapters
  use; collects `<output_basename>*.fasta` (`Native`, "DNA
  Chisel optimized FASTA"), `<output_basename>*.gb` /
  `.genbank` (`Native`, "DNA Chisel GenBank"),
  `<output_basename>*.json` (`Tabular`, "DNA Chisel constraint
  report"), `<output_basename>*.png` (`Native`, "DNA Chisel
  plot"), `*.log`; probe via Python on PATH then `<python> -c
  "import dnachisel"` — on import failure surface as a
  `ProbeReport.warnings` entry, not error),
  `valenx-adapter-lineardesign` (Baidu Research's joint codon
  + secondary-structure mRNA design tool — Apache-2.0, version
  range `1.0.0..2.0.0`, ribbon `bio.lineardesign.design`; the
  modern mRNA-vaccine design workhorse since the 2021
  _Nature_ paper that demonstrated dramatic stability /
  expression gains for mRNA vaccines designed under the joint
  CAI + MFE objective; single-binary CLI subprocess shape
  sister to Phase 18 BWA / Phase 32.5 Smoldyn / Phase 5
  GROMACS with `lineardesign --aa <protein> --lambda
  <lambda_param> --codon_usage <codon_usage> --output_basename
  <basename> [extras...]`; knobs `protein` (path to protein
  FASTA — read in place, no staging) / `output_basename` /
  `lambda_param` (`f64`, finite and ≥ 0.0; default 1.0; the
  Rust field is `lambda_param` because `lambda` is a Rust
  reserved keyword — the CLI emits `--lambda <value>` regardless;
  tunable Lagrangian tradeoff between codon-adaptation-index
  and predicted mRNA secondary-structure stability) /
  `codon_usage` (default `"human"` — selectable from the
  LinearDesign-shipped set: `"human"` / `"mouse"` / `"yeast"` /
  `"ecoli"` / etc.) / `extra_args`; `prepare()` validates
  `lambda_param` is finite and ≥ 0.0 (returns `InvalidCase`
  when negative or NaN), resolves `protein` against the case
  directory when relative, validates the file exists on disk;
  collects `<output_basename>*.fasta` (`Native`, "LinearDesign
  optimized mRNA"), `<output_basename>*.txt` (`Tabular`,
  "LinearDesign report"), `*.log`; probe via
  `find_on_path(&["lineardesign"])` — when the `lineardesign`
  binary isn't found but Python is on PATH the probe surfaces
  a targeted `"clone https://github.com/LinearDesignSoftware/
  LinearDesign and add the bin directory to PATH"` warning),
  and `valenx-adapter-icodon` (Vejnar lab's codon-level mRNA
  stability prediction tool — GPL-3.0, version range
  `1.0.0..2.0.0`, ribbon `bio.icodon.predict`; the canonical
  R-based mRNA stability predictor; ships as a
  `devtools::install_github('santiago1234/iCodon')` R package;
  **Rscript subprocess pattern** sister to Phase 19.6 Seurat —
  the user supplies an `.R` script that loads
  `library(iCodon)` and reads `valenx_params.json` for the
  parsed knobs via `jsonlite::fromJSON`; knobs `script` (`.R`
  enforced) / `rscript` (default `"Rscript"`) / `input_fasta`
  (`Option<PathBuf>`) / `output_basename`; `prepare()` enforces
  the `.R` extension, routes script + optional input_fasta
  through `confined_join` to stage them safely in the workdir,
  writes `valenx_params.json` with the same hand-rolled JSON
  shape as DNA Chisel (key omitted when `None`); collects
  `<output_basename>*.csv` / `*.tsv` (`Tabular`, "iCodon
  stability table"), `<output_basename>*.rds` (`Native`,
  "iCodon R object (RDS)"), `<output_basename>*.png` (`Native`,
  "iCodon plot"), `*.log`; probe via
  `find_on_path(&["Rscript"])` — does not attempt to confirm
  iCodon itself is installed because that would require running
  R, an expensive multi-second startup at probe time; the
  `ToolNotInstalled` install hint mentions the canonical
  `devtools::install_github('santiago1234/iCodon')` install
  path). Three new `valenx-init` templates (`dnachisel`,
  `lineardesign`, `icodon`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **127 of 128**
  fully live, taking the headline number from 124 to 127 total.

- **planemo + Cromwell + cwltool adapters (Phase 22.5).** Sister-
  adapter expansion of the Phase 22 Nextflow + Snakemake workflow-
  manager pair. Three new adapters land under
  `crates/valenx-adapters/bio/` rounding out the bio workflow-
  orchestration surface with three more canonical workflow
  languages: `valenx-adapter-planemo` (the Galaxy project's
  official command-line companion for tool development + workflow
  execution outside a full Galaxy server — AFL-3.0, version range
  `0.75.0..1.0.0`, ribbon `bio.planemo.run`; the same `planemo`
  binary lints tool wrappers, runs Galaxy workflow tests, and
  executes `.ga` workflow files / `.gxwf.yml` Galaxy-flavoured
  workflows; single-binary subprocess shape sister to Phase 22
  Nextflow / Snakemake with `planemo <action> <workflow> [inputs]
  [extras...]`; knobs `workflow` (`.ga` / `.gxwf.yml`; required) /
  `inputs` (`Option<PathBuf>`) / `output_basename` / `action`
  (default `"run"`; rejected at parse time if not in `{run, test,
  lint}`) / `extra_args`; `prepare()` resolves both files against
  case dir, validates each exists; `collect()` walks for
  `<output_basename>*.html` (`Native`, "Planemo report"), `*.json`
  (`Tabular`, "Planemo run JSON"), `*.log` (`Log`); probe via
  `find_on_path(&["planemo"])`), `valenx-adapter-cromwell` (the
  Broad Institute's canonical Workflow Description Language (WDL)
  engine — BSD-3-Clause, version range `80.0.0..100.0.0`, ribbon
  `bio.cromwell.run`; powers most production GATK + Terra
  pipelines; **JAR-distributed** — no `cromwell` launcher binary
  on PATH; the user supplies `[bio.cromwell].jar` absolute path;
  single-binary subprocess shape sister to Phase 33 j5 / Cello /
  Phase 41 Jalview with `java -jar <jar> <action> <workflow> [-i
  <inputs>] [extras...]`; knobs `jar` / `workflow` (`.wdl`) /
  `inputs` (`Option<PathBuf>` — emitted as TWO separate args
  `-i` + `<inputs>` only when `Some`) / `output_basename` /
  `action` (default `"run"`; rejected at parse time if not in
  `{run, submit, validate}`) / `extra_args`; `prepare()` resolves
  all three paths against case dir, validates each exists;
  `collect()` walks **the top level only** of the workdir for
  `<output_basename>*.json` (`Tabular`, "Cromwell metadata"),
  `*.log`; probe via `find_on_path(&["java"])` — Cromwell's
  version comes from the jar itself, not from `java`, so we
  surface no version here, same shape as Phase 33 j5 / Cello /
  Phase 41 Jalview), and `valenx-adapter-cwltool` (the reference
  implementation of the Common Workflow Language — Apache-2.0,
  version range `3.1.0..4.0.0`, ribbon `bio.cwltool.run`; CWL is
  the cross-tool standard for describing analytical workflows in
  YAML / JSON; single-binary subprocess shape sister to Phase 22
  Snakemake with `cwltool --outdir <output_dir> [extras...]
  <workflow> [inputs]`; knobs `workflow` (`.cwl`) / `inputs`
  (`Option<PathBuf>` — JSON or YAML CWL input object) /
  `output_dir` / `extra_args`; `prepare()` resolves both files
  against case dir; `collect()` walks **one level deep** into
  `<output_dir>/` for any file (`Native`, "cwltool output"),
  top-level `*.log`; probe prefers `cwltool` console-script with
  Python-on-PATH fallback + `"cwltool not found on PATH; install
  via pip install cwltool"` warning when only Python is
  present). Three new `valenx-init` templates (`planemo`,
  `cromwell`, `cwltool`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **122 of 123**
  fully live alongside the Phase 42 web-visualization pair,
  taking the headline number to 124 total.

- **Mol* + NGL Viewer adapters (Phase 42).** Opens the **first
  modern web 3D molecular visualization domain** in Valenx. Two
  new adapters land under `crates/valenx-adapters/bio/` spanning
  the web-visualization tradeoff space:
  `valenx-adapter-molstar` (the EMBL-EBI / RCSB-led modern WebGL
  molecular viewer — MIT, version range `3.0.0..5.0.0`, ribbon
  `bio.molstar.view`; de-facto modern molecular viewer embedded
  in the PDB / PDBe / AlphaFold DB / ESM Atlas web properties
  since the late 2010s; wrapped via the `molstar` Python binding
  so it slots into the existing Python-script subprocess pattern
  sister to Phase 17 Biopython / Phase 19.5 Scanpy / Phase 33
  pySBOL / Phase 41 pydna; knobs `script` (`.py` enforced) /
  `python` (default `"python3"`) / `input_structure`
  (`Option<PathBuf>` — optional `.pdb` / `.cif` / `.mmcif`
  structure file) / `output_basename`; `prepare()` enforces
  `.py`, stages script + optional input_structure, writes
  `valenx_params.json` with `output_basename` always plus
  `input_structure` (staged filename) only when set — key
  omitted entirely when `None` rather than emitted as `null`,
  matching the hand-rolled JSON convention the rest of the bio
  adapters use; `collect()` walks for `<output_basename>*.html`
  (`Native`, "Mol* viewer HTML"), `<output_basename>*.molj`
  (`Native`, "Mol* state file" — the JSON state format that
  captures the entire viewer state for reproducible replay),
  `<output_basename>*.png` (`Native`, "Mol* rendered image"),
  `*.log`; probe via Python on PATH then `<python> -c "import
  molstar"` — on import failure surface as a
  `ProbeReport.warnings` entry, not error — sister to the Phase
  19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6 HOOMD-blue
  / Phase 5.7 MDTraj / Phase 41 pydna probe convention) and
  `valenx-adapter-ngl` (the Rose lab's high-performance WebGL
  framework for molecular visualization — MIT, version range
  `3.0.0..5.0.0`, ribbon `bio.ngl.view`; predated Mol* and still
  powers a large fraction of the Jupyter-friendly notebook
  visualization ecosystem via its `nglview` Python binding;
  wrapped via the `nglview` Python binding so it slots into the
  existing Python-script subprocess pattern sister to Mol*;
  knobs `script` (`.py` enforced) / `python` (default
  `"python3"`) / `input_structure` (`Option<PathBuf>`) /
  `output_basename`; `prepare()` mirrors Mol* shape exactly;
  `collect()` walks for `<output_basename>*.html` (`Native`,
  "NGL viewer HTML"), `<output_basename>*.png` (`Native`, "NGL
  rendered image"), `<output_basename>*.json` (`Tabular`, "NGL
  state JSON"), `*.log`; probe via Python on PATH then `<python>
  -c "import nglview"` — sister probe to Mol* with `nglview`-
  specific warning). Two new `valenx-init` templates (`molstar`,
  `ngl`). Each adapter wired into `valenx-app::init_registry`.
  Adapter inventory: **124 of 125** fully live alongside the
  Phase 22.5 workflow-expansion trio. **This completes the bio
  ecosystem from the user's original /review list** — every
  major category called out is now covered, with the bio-adapter
  count totaling **105 bio adapters across 38 biology / biotech
  / chemistry phases**.

- **Smoldyn + MCell adapters (Phase 32.5).** Sister-adapter
  expansion of the Phase 32 systems-biology trio (COPASI /
  BioNetGen / PhysiCell). Two new single-binary CLI subprocess
  adapters land under `crates/valenx-adapters/bio/` rounding out
  the systems-biology / multiscale modeling surface with the
  canonical **spatial stochastic / cell-scale reaction-diffusion**
  simulators that Phase 32 explicitly deferred:
  `valenx-adapter-smoldyn` (Steve Andrews's spatial stochastic
  reaction-diffusion simulator — LGPL-2.1, version range
  `2.70.0..3.0.0`, ribbon `bio.smoldyn.simulate`; resolves
  individual molecules as particles diffusing and reacting in
  continuous 3D space (no lattice discretisation), the canonical
  choice when the question is "where does each molecule actually
  end up over time" rather than "what is the well-mixed
  concentration vs. t" Phase 32 COPASI's ODE / SSA covers; single-
  binary subprocess shape sister to Phase 18 BWA with `smoldyn
  <config> [extras...]`; knobs `config` (Smoldyn `.txt`
  configuration describing simulation geometry — boundaries,
  surfaces, compartments — plus molecule species + diffusion
  coefficients and per-pair / per-surface reactions; required) /
  `extra_args`; `prepare()` resolves `config` against the case
  directory when relative, validates it exists on disk;
  `collect()` walks for `*.txt` (`Tabular`, "Smoldyn output
  table" — Smoldyn's per-step particle / reaction tables), `*.dat`
  (`Tabular`, "Smoldyn data" — reaction-event / molecule-position
  dumps the config may direct here), `*.log` (`Log`, "Smoldyn
  log"); probe via `find_on_path(&["smoldyn"])`) and
  `valenx-adapter-mcell` (Salk Institute / Stiles, Bartol's cell-
  scale Monte Carlo spatial stochastic simulator — GPL-2.0,
  version range `4.0.0..5.0.0`, ribbon `bio.mcell.simulate`;
  walks the user's `.mdl` (Model Description Language) model —
  geometry built from triangle meshes, molecule species with
  diffusion coefficients, surface / volume reactions, release
  patterns, observation counts — and runs Brownian-dynamics
  particle trajectories with Monte Carlo reaction sampling;
  canonical use case is sub-cellular signaling — synaptic
  transmission, calcium dynamics, receptor binding — where
  geometry is intricate enough that Smoldyn's continuous-space
  mode would be overkill but a well-mixed COPASI / BioNetGen
  treatment misses the spatial structure; single-binary
  subprocess shape sister to Smoldyn with `mcell [-seed <N>]
  <mdl> [extras...]` (the `-seed` flag and its integer argument
  emitted as TWO separate OsString tokens, only when `seed` is
  `Some(_)`); knobs `mdl` (`.mdl` MCell model description file;
  required) / `seed` (`Option<u32>` — when `Some(n)` MCell uses
  that seed; when `None` MCell picks its own seed and prints it
  on the run banner — same shape as the Phase 29 SLiM `-s` and
  Phase 30.5 BEAST 2 `-seed` knobs) / `extra_args`; `prepare()`
  resolves `mdl` against the case directory when relative,
  validates it exists on disk; `collect()` walks for `*.dat`
  (`Tabular`, "MCell reaction data" — per-observation count
  tables MCell writes from the model's REACTION_DATA_OUTPUT
  block), `*.dx` (`Native`, "MCell visualization data" — DReAMM /
  OpenDX visualization frames), `*.log`; probe via
  `find_on_path(&["mcell"])`). Two new `valenx-init` templates
  (`smoldyn`, `mcell`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **114 of 115**
  fully live alongside the Phase 40 microscopy trio and Phase 41
  sequence-editors pair, taking the headline number to 119 total.

- **Fiji + CellProfiler + Ilastik adapters (Phase 40).** Opens
  the **first microscopy / bioimage analysis domain** in Valenx.
  Three new adapters land under `crates/valenx-adapters/bio/`
  spanning the bioimage analysis tradeoff space:
  `valenx-adapter-fiji` (the [Fiji Is Just ImageJ](https://fiji.sc/)
  distribution of NIH ImageJ — Schindelin et al, GPL-3.0, version
  range `2.0.0..3.0.0`, ribbon `bio.fiji.process`; bundles
  ImageJ2 + a curated set of plugins for biological image
  processing — channel splitting, thresholding, particle analysis,
  deconvolution, registration, segmentation, the TrackMate single-
  particle tracker, the BoneJ trabecular-bone toolkit, the entire
  ImageJ macro / Jython / scripting surface; per-platform launcher
  binaries `ImageJ-linux64`, `ImageJ.exe`,
  `Contents/MacOS/ImageJ-macosx`; app-launcher subprocess shape
  sister to Phase 36 RELION / EMAN2 with `<fiji_app> --headless
  --console -macro <macro_file> [extras...]`; knobs `fiji_app`
  (absolute path to per-platform Fiji launcher; required) /
  `macro_file` (`.ijm` Fiji macro; required) / `input_image`
  (`Option<PathBuf>` — optional input image the macro will operate
  on, typically picked up via `getArgument()`; the macro is
  responsible for opening it) / `output_basename` / `extra_args`;
  `prepare()` resolves `fiji_app`, `macro_file`, and the optional
  `input_image` against the case directory when relative,
  validates each existing file exists on disk; `collect()` walks
  for `<output_basename>*.tif` / `.tiff` (`Native`, "Fiji image
  (TIFF)"), `<output_basename>*.png` (`Native`, "Fiji image
  (PNG)"), `<output_basename>*.csv` (`Tabular`, "Fiji
  measurements"), `*.log`; probe via
  `find_on_path(&["ImageJ-linux64", "ImageJ-macosx",
  "ImageJ.exe", "fiji"])` — surfaces a Java-fallback warning
  whenever Fiji isn't on PATH but `java` is),
  `valenx-adapter-cellprofiler` (Broad Institute's pipeline-
  driven cell segmentation + measurement suite — BSD-3-Clause,
  version range `4.0.0..5.0.0`, ribbon `bio.cellprofiler.segment`;
  canonical tool for high-content screening — the user authors a
  `.cppipe` pipeline in the GUI and the CLI runs that pipeline
  over a directory of input images, emitting per-object
  measurement CSVs + segmented label-image overlays; Python-CLI
  subprocess shape with bundled `cellprofiler` launcher as the
  primary entry point and `<python> -m cellprofiler ...` as the
  fallback when the launcher isn't on PATH but Python is;
  `cellprofiler -c -r -p <pipeline> -i <input_dir> -o <basename>
  [extras...]`; knobs `pipeline` (`.cppipe` / `.cpproj`;
  required) / `input_dir` (input image directory; required — the
  adapter validates it is a directory at prepare time) /
  `output_basename` / `python` (default `"python3"`) /
  `extra_args`; `prepare()` resolves `pipeline` and `input_dir`
  against the case directory when relative, validates `pipeline`
  exists on disk and `input_dir` is a directory, looks up the
  `cellprofiler` binary on PATH first then falls back to `<python>
  -m cellprofiler ...`; `collect()` walks **one level deep** into
  `<output_basename>/` for `*.csv` (`Tabular`, "CellProfiler
  measurements"), `*.tif` / `.tiff` (`Native`, "CellProfiler
  segmented image"), `*.png` (`Native`, "CellProfiler plot");
  top-level `*.log` (`Log`, "CellProfiler log"); probe via
  `find_on_path(&["cellprofiler", "python3", "python"])` with
  warning when `cellprofiler` itself isn't on PATH but Python is),
  and `valenx-adapter-ilastik` (Hamprecht lab's interactive-ML
  pixel / object classification suite — GPL-3.0, version range
  `1.4.0..2.0.0`, ribbon `bio.ilastik.classify`; leans on user-
  trained random-forest classifiers — the user paints a few
  foreground / background strokes per image in the GUI to teach
  the classifier, saves the resulting `.ilp` project file, and
  then runs the headless CLI to apply that trained classifier to
  a batch of new images; canonical use case is hard segmentation
  tasks where rule-based pipelines (CellProfiler) or threshold-
  driven macros (Fiji) struggle — light-sheet imagery, tissue
  cross-sections, anything with low contrast or irregular
  textures; app-launcher subprocess shape sister to Fiji with
  `<ilastik_app> --headless --project=<project>
  --output_filename_format=<basename>_{nickname}.h5
  <input_images...> [extras...]`; the `--project=` and
  `--output_filename_format=` flags are emitted as single OsString
  args each so `=` and the value travel together, and the literal
  `{nickname}` substring in the format string is Ilastik's per-
  image nickname placeholder — must reach Ilastik unmodified for
  per-input-image output disambiguation; knobs `ilastik_app`
  (absolute path to per-platform Ilastik launcher; required) /
  `project` (`.ilp` Ilastik project file containing the trained
  classifier; required) / `input_images` (`Vec<PathBuf>` — must
  contain ≥ 1 entry; the adapter rejects an empty vector at
  prepare time) / `output_basename` / `workflow` (default `"Pixel
  Classification"` — selectable from Ilastik's set: `"Pixel
  Classification"`, `"Object Classification"`, etc.) /
  `extra_args`; `collect()` walks for `<output_basename>*.h5`
  (`Native`, "Ilastik probability map (HDF5)"),
  `<output_basename>*.tif` (`Native`, "Ilastik segmentation"),
  `*.log`; probe via `find_on_path(&["ilastik", "run_ilastik.sh",
  "ilastik.exe"])` with warning when nothing matches but still
  returns `ok = true` since the user can supply the launcher via
  `case.toml` — sister to Phase 32 PhysiCell's per-project-binary
  probe convention). Three new `valenx-init` templates (`fiji`,
  `cellprofiler`, `ilastik`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **117 of 118**
  fully live alongside the Phase 32.5 spatial-stochastic pair and
  Phase 41 sequence-editors pair, taking the headline number to
  119 total.

- **pydna + Jalview adapters (Phase 41).** Opens the **first
  plasmid-design / alignment-viewer domain** in Valenx. Two new
  adapters land under `crates/valenx-adapters/bio/` spanning the
  sequence-editor tradeoff space: `valenx-adapter-pydna` (Bjorn
  Johansson's Python plasmid / clone-design library — BSD-3-Clause,
  version range `5.0.0..7.0.0`, ribbon `bio.pydna.design`;
  handles the long tail of cloning operations programmatically —
  PCR primer design, restriction-enzyme digests, Gibson assembly,
  Golden-Gate assembly, sequence-overlap detection, ligation
  simulation, primer Tm calculation, cloning-strategy validation;
  canonical use case is "I need to assemble these N parts into
  this target construct — what primers should I order, what
  enzymes should I cut with, and does the end-product match my
  target sequence?" — pydna replaces hours of manual work in ApE /
  SnapGene / Vector NTI with a few dozen lines of Python; Python-
  script subprocess shape sister to Phase 17 Biopython / Phase
  19.5 Scanpy / Phase 33 pySBOL; knobs `script` (`.py` enforced)
  / `python` (default `"python3"`) / `input_genbank`
  (`Option<PathBuf>` — optional starting GenBank file the script
  can use as the parent / template construct; `None` when the
  script generates the design from scratch) / `output_basename`;
  `prepare()` enforces `.py`, stages script + optional
  input_genbank, writes `valenx_params.json` with `output_basename`
  always plus `input_genbank` (staged filename) only when set —
  key omitted entirely when `None` rather than emitted as `null`,
  matching the hand-rolled JSON convention the rest of the bio
  adapters use (Phase 19.6 Seurat / AnnData, Phase 27.5 ESM-IF);
  `collect()` walks for `<output_basename>*.gb` / `.genbank`
  (`Native`, "pydna GenBank file"), `<output_basename>*.fasta`
  (`Native`, "pydna FASTA"), `<output_basename>*.csv` (`Tabular`,
  "pydna table"), `*.log`; probe via Python on PATH with `import
  pydna` check — on import failure surface as a
  `ProbeReport.warnings` entry, not error — sister to the Phase
  19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6 HOOMD-blue
  / Phase 5.7 MDTraj probe convention) and
  `valenx-adapter-jalview` (the Barton group's Java alignment
  viewer — GPL-3.0, version range `2.11.0..3.0.0`, ribbon
  `bio.jalview.view`; the reference alignment viewer in molecular
  biology labs since the 2000s — multiple-sequence alignment
  viewing + editing, conservation / consensus / occupancy plots,
  structural overlays via Jmol / Chimera links, per-column
  annotations; ships a **headless mode** (`-nodisplay`) for batch
  image / format conversion — the user feeds an alignment in,
  picks an output format (PNG image, HTML report, SVG vector
  graphic, FASTA / Clustal alignment re-export), and Jalview
  writes the requested artifact without opening its GUI; **JAR-
  distributed** — single-binary subprocess shape sister to Phase
  33 j5 / Cello with `java -jar <jar> -nodisplay -open <input>
  -<output_format> <basename>.<ext> [extras...]`; knobs `jar`
  (absolute path to the Jalview jar; required) / `input`
  (alignment input — `.fa` / `.aln` / `.clustal` / `.stockholm`
  and friends Jalview reads natively; required) /
  `output_basename` / `output_format` (default `"png"` —
  selectable from `"png"` / `"html"` / `"svg"` / `"fasta"` /
  `"clustal"`) / `extra_args`; `prepare()` derives the output
  extension from `output_format` (png → .png, html → .html, svg
  → .svg, fasta → .fasta, clustal → .aln, default → use the
  format string itself as extension); `collect()` walks for
  `<output_basename>*.png` (`Native`, "Jalview alignment image"),
  `<output_basename>*.svg` (`Native`, "Jalview SVG"),
  `<output_basename>*.html` (`Native`, "Jalview HTML"),
  `<output_basename>*.fasta` (`Native`, "Jalview FASTA"),
  `<output_basename>*.aln` (`Tabular`, "Jalview alignment"),
  `*.log`; probe via `find_on_path(&["java"])` — Jalview's
  version comes from the jar itself, not from `java`, so we
  surface no version here — the user pins the Jalview release
  implicitly by the jar they point at; same shape as Phase 33
  j5 / Cello). Two new `valenx-init` templates (`pydna`,
  `jalview`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **119 of 120**
  fully live alongside the Phase 32.5 spatial-stochastic pair
  and Phase 40 microscopy trio. **This crosses the 100-bio-
  adapter milestone** — bio-adapter count now totals 100 across
  36 biology / biotech / chemistry phases.

- **NAMD + AmberTools sander + HOOMD-blue adapters (Phase 5.6).**
  Sister-domain expansion of the Phase 5 GROMACS / LAMMPS MD
  engine beachhead. Three new adapters land under
  `crates/valenx-adapters/bio/` rounding out the all-atom + GPU-
  native MD-engine surface alongside the Phase 17 OpenMM Python-
  native engine: `valenx-adapter-namd` (UIUC's flagship all-atom
  MD engine — custom NAMD-License — academic / non-commercial use
  only, version range `2.14.0..4.0.0`, ribbon
  `bio.namd.simulate`; the de-facto choice in biomolecular MD
  pedagogy and a workhorse on every academic HPC cluster; NAMD
  2.x ships an SMP-threaded CHARMM-style integrator, NAMD 3.x
  adds GPU-resident kernels; single-binary subprocess shape with
  `<binary> +p<processors> <config> [extras...]` where
  `<binary>` is `namd2` or `namd3` (probe accepts either) and
  `+p<N>` (no space — NAMD's own flag syntax) is NAMD's threading
  flag emitted as a single OsString so the flag and value travel
  together exactly as NAMD parses them; knobs `config` (NAMD
  `.namd` / `.conf` configuration; required) / `processors`
  (`u32`, default 1) / `extra_args`; `prepare()` resolves
  `config` against the case directory when relative; `collect()`
  walks for `*.dcd` (`Native`, "NAMD trajectory (DCD)"), `*.coor`
  (`Native`, "NAMD coordinates"), `*.vel` (`Native`, "NAMD
  velocities"), `*.xsc` (`Tabular`, "NAMD extended system"),
  `*.log`; probe via `find_on_path(&["namd2", "namd3"])` pushes
  an `"academic"`-keyworded warning containing both `"academic"`
  and `"non-commercial"` substrings into `ProbeReport.warnings`
  whenever the binary is detected, and `tool_license` surfaces as
  `"NAMD-License"` rather than mislabeling the custom UIUC NAMD
  terms as a recognised SPDX identifier),
  `valenx-adapter-amber-sander` (AmberTools' OSS MD engine
  portion — sander itself is GPL-3.0 OSS; the proprietary
  `pmemd.cuda` GPU engine is NOT wrapped, version range
  `22.0.0..26.0.0`, ribbon `bio.sander.simulate`; sister to Phase
  5.5 cpptraj — installing AmberTools installs both; single-
  binary subprocess shape with `sander -O -i <config> -p
  <topology> -c <coordinates> -o <basename>.out -r <basename>.rst
  -x <basename>.nc [extras...]` (the `-O` flag overwrites
  existing outputs — standard re-run convention); knobs
  `topology` (`.prmtop` / `.parm7`; required) / `coordinates`
  (`.inpcrd` / `.rst7`; required) / `config` (`.in` / `.mdin`;
  required) / `output_basename` / `extra_args`; `prepare()`
  resolves all three input paths against the case directory when
  relative, validates each file exists on disk; `collect()`
  walks for `<output_basename>*.out` (`Log`, "sander mdout"),
  `<output_basename>*.nc` (`Native`, "sander NetCDF trajectory"),
  `<output_basename>*.rst` (`Native`, "sander restart
  coordinates"), `<output_basename>*.mdinfo` (`Log`, "sander
  mdinfo"); probe via `find_on_path(&["sander"])` — no academic-
  license caveat), and `valenx-adapter-hoomd` (Glotzer lab's GPU-
  native particle simulator — BSD-3-Clause, version range
  `3.0.0..6.0.0`, ribbon `bio.hoomd.simulate`; HOOMD-blue v3+ is
  fully Python-scripted (no native CLI) — the user supplies a
  `.py` script that does `import hoomd` and runs the simulation;
  canonical engine for soft-matter / coarse-grained particle
  systems, polymers, colloids, rigid-body assemblies — sister to
  LAMMPS in the particle-MD surface but GPU-first by design;
  Python-script subprocess shape sister to Phase 17 OpenMM;
  knobs `script` (`.py` enforced) / `python` (default
  `"python3"`) / `output_basename`; `prepare()` enforces the
  `.py` extension, stages the script into the workdir under its
  original filename, writes a flat `valenx_params.json`
  containing `output_basename`, builds `<python>
  <staged_script>`; `collect()` walks for
  `<output_basename>*.gsd` (`Native`, "HOOMD trajectory (GSD)"),
  `<output_basename>*.h5` (`Native`, "HOOMD HDF5 output"),
  `*.log`; probe via `find_on_path(&["python3", "python"])` then
  `<python> -c "import hoomd"` — on import failure surface as a
  `ProbeReport.warnings` entry not error). Three new
  `valenx-init` templates (`namd`, `sander`, `hoomd`). Each
  adapter wired into `valenx-app::init_registry`. Adapter
  inventory: **108 of 109** fully live.

- **MDTraj adapter (Phase 5.7).** Single-adapter sister to the
  Phase 17 MDAnalysis adapter and the Phase 5.5 PLUMED / ProDy /
  cpptraj analysis trio. Single-adapter phases are a precedent
  in Valenx — when an established tool fills a clearly-defined
  corner of an existing surface without requiring new
  infrastructure, the phase ships as a single adapter. New
  adapter lands under `crates/valenx-adapters/bio/`:
  `valenx-adapter-mdtraj` (Pande / VanderSpoel / Beauchamp lab's
  Python MD trajectory analysis library — LGPL-2.1, version range
  `1.9.0..2.0.0`, ribbon `bio.mdtraj.analyze`; the second-most-
  used Python MD trajectory analyzer alongside MDAnalysis with
  wider format support (`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr`
  / `.binpos` / `.lh5` / `.amber` / `.gromacs`) and deeper
  integration with the OpenMM ecosystem (the Pande / Beauchamp
  lab is co-located with the OpenMM developers — MDTraj's HDF5
  trajectory format is OpenMM's native streaming output) plus a
  pandas-friendly per-frame property API; Python-script
  subprocess shape sister to Phase 17 Biopython, Phase 5.5
  ProDy, Phase 17 OpenMM; knobs `script` (`.py` enforced) /
  `python` (default `"python3"`) / `trajectory` (`.xtc` /
  `.dcd` / `.h5` / `.nc` / `.trr` / `.binpos` / `.lh5` MDTraj-
  supported trajectory; required) / `topology` (`.pdb` /
  `.prmtop` / `.gro` / `.psf` topology MDTraj uses for atom +
  residue + chain metadata; required) / `output_basename`;
  `prepare()` enforces the `.py` extension on the script,
  resolves all three input paths against the case directory when
  relative, stages script + trajectory + topology into the
  workdir under their original filenames so the script can
  resolve them via relative paths, then writes a flat
  `valenx_params.json` containing `output_basename`, the bare
  `trajectory` filename, and the bare `topology` filename, builds
  `<python> <staged_script>`; `collect()` walks for
  `<output_basename>*.csv` (`Tabular`, "MDTraj analysis table"),
  `<output_basename>*.npz` (`Native`, "MDTraj numpy archive"),
  `<output_basename>*.h5` (`Native`, "MDTraj HDF5 output"),
  `<output_basename>*.png` (`Native`, "MDTraj plot"), `*.log`;
  probe via `find_on_path(&["python3", "python"])` then
  `<python> -c "import mdtraj"` — on import failure surface as a
  `ProbeReport.warnings` entry not error). One new `valenx-init`
  template (`mdtraj`). Adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **109 of 110**
  fully live alongside the Phase 5.6 bio MD-engine trio.

- **RoseTTAFold + OmegaFold + FoldSeek adapters (Phase 17.7).**
  Sister-adapter expansion of the Phase 17.5 structure-prediction
  beachhead and the Phase 17 ColabFold adapter. Three new
  adapters land under `crates/valenx-adapters/bio/` rounding out
  the protein structure prediction + structure search surface
  that Phase 17.5 ESMFold / OpenFold / AlphaFold 2 / AlphaFold 3
  + Phase 17 ColabFold opened: `valenx-adapter-rosettafold`
  (Baker lab's original 3-track structure-prediction network —
  MIT, version range `1.0.0..3.0.0`, ribbon
  `bio.rosettafold.predict`; three concurrent attention tracks
  over the 1D sequence, the 2D pair-distance map, and the 3D
  Cartesian backbone with cross-track message passing refining
  all three jointly; canonical pre-AlphaFold-3 sibling that
  established the 3-track SE(3)-equivariant attention pattern;
  Python-script subprocess shape sister to Phase 17.5 ESMFold;
  knobs `script` (`.py` enforced) / `python` (default
  `"python3"`) / `fasta` (input FASTA query sequence; required) /
  `output_basename`; `prepare()` enforces `.py`, stages script +
  fasta, writes `valenx_params.json` with `output_basename` +
  bare `fasta` filename, builds `<python> <staged_script>`;
  `collect()` walks for `<output_basename>*.pdb` (`Native`,
  "RoseTTAFold predicted structure" — pLDDT-style per-residue
  confidence in the B-factor column lifted by the existing Phase
  17 PDB reader without any structure-prediction-specific code
  path), `<output_basename>*.npz` (`Native`, "RoseTTAFold
  confidence arrays"), `*.log`; probe via
  `find_on_path(&["python3", "python"])` — deliberately doesn't
  try `import rosettafold` (RoseTTAFold is not a pip package);
  pushes a probe warning whenever Python is detected: "RoseTTAFold
  model weights + dependencies not bundled — clone
  https://github.com/RosettaCommons/RoseTTAFold and follow the
  install README"), `valenx-adapter-omegafold` (HelixonAI's
  single-sequence protein-structure predictor — Apache-2.0,
  version range `1.0.0..2.0.0`, ribbon `bio.omegafold.predict`;
  MSA-free like ESMFold but uses a larger pre-trained transformer
  backbone trained on a much wider sequence corpus; works on
  single sequences but routinely matches AlphaFold-2-with-MSA
  quality on orphan / synthetic / fast-evolving sequences where
  MSA-based methods struggle; ships its own CLI binary
  (`omegafold <fasta> <output_dir>`) and falls back to `<python>
  -m omegafold ...` when the CLI launcher isn't on PATH but
  Python is; knobs `fasta` / `output_basename` (workdir-relative
  output directory name) / `python` (default `"python3"`; used
  only as fallback) / `model_dir` (`Option<PathBuf>` — optional
  pre-downloaded model checkpoint directory); `prepare()` builds
  `omegafold <fasta> <output_basename> [--model <model_dir>]`
  with the FASTA passed by absolute path (NOT staged into the
  workdir); `collect()` walks one level deep into the
  `<output_basename>/` subdirectory for `*.pdb` (`Native`,
  "OmegaFold predicted structure") and `*.json` (`Log`,
  "OmegaFold metadata"), plus the workdir-top-level `*.log`;
  probe via `find_on_path(&["omegafold", "python3", "python"])`
  — surfaces a warning if `omegafold` itself isn't on PATH but
  Python is ("OmegaFold CLI not found on PATH; install via pip
  install git+https://github.com/HeliXonProtein/OmegaFold.git")),
  and `valenx-adapter-foldseek` (Steinegger lab's protein-
  structure search via the 3Di alphabet — GPL-3.0, version range
  `8.0.0..10.0.0`, ribbon `bio.foldseek.search`; the **3D
  analogue of Phase 18.5 MMseqs2** with both built on the same
  fast many-vs-many search core but FoldSeek encoding the per-
  residue 3D geometry as a custom "3Di alphabet" (a 20-letter
  alphabet over local backbone geometry patterns, designed so
  structural matches have high 3Di alphabet identity) and running
  3Di-vs-3Di comparisons at sequence-search speed — finds
  structural homologs at PDB-scale in seconds rather than the
  hours / days HMM-based or geometry-based search tools take;
  single-binary subprocess shape sister to Phase 18.5 MMseqs2 /
  Phase 18 BWA with `foldseek easy-search <query> <database>
  <basename>.m8 tmp_<basename> --threads <N> [extras...]` (the
  `tmp_<basename>` is a per-run temp directory FoldSeek
  requires); knobs `query` (`.pdb` / `.cif` query structure;
  required) / `database` (FoldSeek database path prefix — the
  user supplies the path stem and FoldSeek resolves the
  `<prefix>_*` sidecar files itself; required) / `output_basename`
  / `threads` (`u32`, default 1) / `extra_args`; `prepare()`
  resolves `query` against the case directory when relative,
  validates the `database` parent directory exists on disk (the
  database files themselves use the prefix convention so we
  cannot validate them by name — same shape as Phase 18.7
  BLAST+'s `database` validation); `collect()` walks for
  `<output_basename>.m8` (`Tabular`, "FoldSeek search results"
  — the canonical BLAST-style M8 hit table format) and `*.log`;
  the temp directory is not surfaced in artifacts — it's
  intermediate; probe via `find_on_path(&["foldseek"])`). Three
  new `valenx-init` templates (`rosettafold`, `omegafold`,
  `foldseek`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **112 of 113**
  fully live alongside the Phase 5.6 bio MD-engine trio and
  Phase 5.7 MDTraj single-adapter.

- **BLAST+ + Clustal Omega + T-Coffee adapters (Phase 18.7).**
  Sister-adapter expansion of Phase 18 / 18.5 / 18.6 alignment
  beachhead. Three new single-binary CLI subprocess adapters land
  under `crates/valenx-adapters/bio/` rounding out the foundational
  sequence-alignment surface that the existing BWA / minimap2 /
  MAFFT / MUSCLE adapters explicitly left out: `valenx-adapter-blast`
  (NCBI BLAST+ — Public Domain US-government work, version range
  `2.10.0..3.0.0`, ribbon `bio.blast.search`; ships five user-facing
  search programs `blastn` / `blastp` / `blastx` / `tblastn` /
  `tblastx` covering every nucleotide / protein search direction;
  knobs `program` / `query` / `database` (path prefix) / `evalue`
  (default 10.0) / `outfmt` (default 0) / `threads` (default 1);
  `prepare()` validates the database parent directory exists on
  disk and looks up the per-program binary at prepare time;
  `collect()` walks for `blast_results.txt` (`Tabular`, "BLAST
  search results") + `*.log`), `valenx-adapter-clustalo` (Clustal
  Omega — GPL-2.0, version range `1.2.0..2.0.0`, ribbon
  `bio.clustalo.align`; modern HMM-driven progressive multiple-
  sequence aligner; knobs `input` / `output_basename` / `outfmt`
  default `"clustal"` / `threads` default 1; `prepare()` derives
  `<ext>` from `outfmt` (clustal → .aln, fasta → .fasta, phylip →
  .phy, vienna → .vie, nexus → .nex); `collect()` walks for
  `<output_basename>*` + `*.log`), and `valenx-adapter-tcoffee`
  (T-Coffee — GPL-2.0, version range `13.0.0..14.0.0`, ribbon
  `bio.tcoffee.align`; library-based consistency-weighted multiple-
  sequence aligner; knobs `input` / `output_basename` / `outfmt`
  default `"clustalw"` / `mode` `Option<String>` — omit for default
  progressive mode; output always pinned to `.aln` via T-Coffee's
  `=`-style flag form; `collect()` walks for `<output_basename>*` +
  `*.dnd` (Newick guide tree) + `*.log`; probe via
  `find_on_path(&["t_coffee"])` — note the underscore: T-Coffee
  installs as `t_coffee`, not `t-coffee`). Three new `valenx-init`
  templates (`blast`, `clustalo`, `tcoffee`). Each adapter wired
  into `valenx-app::init_registry`. Adapter inventory: **103 of
  104** fully live. Cross-binary roundtrip test sweeps all 101
  templates clean alongside the Phase 19.6 single-cell-expansion
  pair.

- **Seurat + AnnData adapters (Phase 19.6).** Sister-adapter
  expansion of the Phase 19.5 single-cell Scanpy + scVI beachhead.
  Two new adapters land under `crates/valenx-adapters/bio/`
  rounding out the single-cell genomics surface that Phase 19.5
  explicitly deferred: `valenx-adapter-seurat` (Satija lab's
  R-based single-cell analysis toolkit — MIT, version range
  `4.0.0..6.0.0`, ribbon `bio.seurat.analyze`; introduces the
  **Rscript subprocess pattern** to Valenx — the R analogue of the
  Python-script pattern that Phase 17 Biopython / Phase 19.5
  Scanpy / Phase 33 pySBOL established; user supplies an `.R`
  script that loads `library(Seurat)` and reads
  `valenx_params.json` for the parsed knobs via `jsonlite::
  fromJSON`; knobs `script` (`.R` enforced) / `rscript` default
  `"Rscript"` / `input_data` `Option<PathBuf>` — supports `.h5` /
  `.mtx` / `.rds` so users can drop in 10x HDF5, sparse Matrix
  Market, or pre-saved Seurat object formats / `output_basename`;
  `prepare()` enforces the `.R` extension, stages script + optional
  input_data, writes `valenx_params.json` with `output_basename`
  and `input_data` (staged filename when set; key omitted entirely
  when `None` rather than emitted as `null`, matching the hand-
  rolled JSON convention the rest of the bio adapters use);
  `collect()` walks for `<output_basename>*.rds` (`Native`, "Seurat
  object (RDS)" — canonical R-serialised Seurat object format
  consumed by every downstream Seurat / signac / Azimuth pipeline)
  / `.csv` / `.png` / `*.log`; probe via `find_on_path(&
  ["Rscript"])` — does not attempt to confirm Seurat itself is
  installed because that would require running R, an expensive
  multi-second startup at probe time) and `valenx-adapter-anndata`
  (scverse's Python single-cell HDF5 data container library —
  BSD-3-Clause, version range `0.9.0..1.0.0`, ribbon
  `bio.anndata.process`; the canonical container that ties the
  entire scverse Python ecosystem together — scanpy / scvi /
  scirpy / squidpy / muon all read and write `.h5ad`; Python-script
  subprocess shape sister to Phase 19.5 Scanpy / scVI; knobs
  `script` (`.py` enforced) / `python` default `"python3"` /
  `input_h5ad` `Option<PathBuf>` — supports `.h5ad` / `.h5` /
  `output_basename`; `prepare()` enforces the `.py` extension,
  stages, writes `valenx_params.json` with the same hand-rolled
  shape as Seurat (key omitted when `None`); `collect()` walks for
  `<output_basename>*.h5ad` (`Native`, "AnnData h5ad file") /
  `.csv` / `.png` / `*.log`; probe via Python on PATH then
  `<python> -c "import anndata"` — on import failure surface as a
  `ProbeReport.warnings` entry, not error). Two new `valenx-init`
  templates (`seurat`, `anndata`). Each adapter wired into
  `valenx-app::init_registry`. Adapter inventory: **105 of 106**
  fully live. Cross-binary roundtrip test sweeps all 101 templates
  clean alongside the Phase 18.7 alignment-toolkit-expansion trio.

**🎉 100-adapter milestone.** With this Phase 5.5 + Phase 33 docs
pass, **Valenx now ships 100 live adapters** spanning the physics-
domain phases 1-9 plus the entire Phase 5.5 / 17 → 39 biology /
biotech / chemistry expansion. Adapter inventory at 100 of 101
fully live (only `occt` remains stub-only — needs the `occt-sys`
C++ FFI shim). Phase 5.5 MD analysis expansion (PLUMED + ProDy +
cpptraj) and Phase 33 synthetic biology (pySBOL + j5 + Cello)
land together in this combined phase pair, taking the headline
adapter count from 94 to 100 and the cross-binary roundtrip from
90 to 96 templates.

### Phase 33 — Synthetic biology

Open the **first synthetic biology / genetic-circuit design
domain** in Valenx with three established open-source tools that
span the synthetic-biology tradeoff space — canonical SBOL-
standard Python composition (pySBOL), DNA assembly automation that
plans the optimal Gibson / Golden-Gate / SLIC / SLIM strategy from
a target circuit + parts library (j5), and genetic-circuit DNA
compilation from a Verilog netlist describing the desired logic
function (Cello v2): pySBOL (the Python implementation (pySBOL3)
of the Synthetic Biology Open Language standard, Apache-2.0 — SBOL
captures components, sequences, interactions, constraints, and the
full provenance of a synthetic design as RDF/XML or JSON-LD that
round-trips with every SBOL-conformant tool — j5, Cello, SynBioHub,
iBioSim; Python-script subprocess shape sister to Phase 17
Biopython), j5 (JBEI's canonical DNA-assembly automation tool,
BSD-3-Clause — consumes a target circuit design (CSV row per
cassette) plus a parts library (CSV row per part / oligo), then
plans the optimal Gibson / Golden-Gate / SLIC / SLIM assembly
strategy and writes the per-step protocol + GenBank construct
files; **JAR-distributed** — the user supplies the absolute path
to `j5.jar` via case input, and we probe `java` itself rather than
the JAR), and Cello (CIDAR's canonical genetic-circuit DNA
compiler — Cello v2, BSD-3-Clause — consumes a Verilog netlist
describing the desired logic function plus a triplet of JSON
constraint files (a user constraint file pinning the chassis /
library, an input sensor file pinning the input promoters, an
output device file pinning the reporter) and emits a fully
assembled DNA construct that implements the logic in a living cell,
running a simulated-annealing optimization over the gate-assignment
problem and outputting a Graphviz `.dot` netlist + circuit diagram
PNG + human-readable report; **JAR-distributed** — same shape as
j5).

**This is the first synthetic biology / genetic-circuit design
domain to land in Valenx.** The biology adapter family started
with Phase 17 (foundation — sequence / structure / trajectory
canonical types + classical MD + cheminformatics scripts) and
expanded through Phase 5.5 / 17.5 / 18 / 18.5 / 18.6 / 19 / 19.5
/ 20 / 22 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5
/ 31 / 32 / 34 / 35 / 36 / 38 / 39 to cover MD-trajectory analysis
expansion, sequence prediction, alignment, RNA-seq, variant
calling, single-cell, transcript quantification, workflow
orchestration, molecular viewers, cheminformatics, quantum
chemistry, protein design, EvolutionaryScale models, RNA structure,
population genetics, phylogenetics, Bayesian phylogenetics,
sequencing read simulation, systems biology, small-molecule
docking, CRISPR design, cryo-EM reconstruction, Rosetta protein
modeling, and DNA structural geometry — but until Phase 33 the
synthetic-biology / genetic-circuit-design surface (SBOL-standard
composition, DNA assembly automation, Verilog → DNA circuit
compilation) was absent. Phase 33 closes that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-pysbol` — the Python implementation (pySBOL3) of
  the Synthetic Biology Open Language standard (Apache-2.0). SBOL
  captures components, sequences, interactions, constraints, and
  the full provenance of a synthetic design as RDF/XML or JSON-LD
  that round-trips with every SBOL-conformant tool — j5, Cello,
  SynBioHub, iBioSim. Python-script subprocess shape (sister to
  Phase 17 Biopython): the user supplies a Python script
  referenced from `[bio.pysbol].script` in `case.toml` that imports
  `sbol3` and reads `valenx_params.json` for the parsed knobs.
  `[bio.pysbol]` knobs: `script` (path to user-supplied Python
  script; required), `python` (interpreter name; default
  `"python3"`), `input_sbol` (optional starting SBOL XML document;
  `None` when the script generates the design from scratch),
  `output_basename` (filename stem the user's script uses for
  outputs — surfaced here so collect() can label artifacts
  uniformly; required, non-empty). `prepare()` stages the script +
  optional input SBOL into the workdir under their original
  filenames so the script can resolve them via relative paths,
  then writes a flat `valenx_params.json` containing `input_sbol`
  (staged filename or literal `null`) and `output_basename`.
  `collect()` walks the workdir for `<output_basename>*.xml`
  (`Tabular`, "pySBOL document") and `<output_basename>*.json`
  (`Log`, "pySBOL composition log"). Probe via Python on PATH with
  an `import sbol3` check (returns `ok = true` with a warning when
  import fails so non-standard installs aren't blocked). Version
  range `3.0.0..4.0.0` (pySBOL3 is the modern Python rewrite — the
  older 2.x line is deprecated; 3.0 is the floor; upper bound 4.0
  reserves room for an eventual major bump). The init alias `sbol`
  resolves to the same template as the canonical `pysbol` name.
  `bio.pysbol.compose` ribbon capability.
- `valenx-adapter-j5` — JBEI's canonical DNA-assembly automation
  tool (BSD-3-Clause). j5 consumes a target circuit design (CSV
  row per cassette) plus a parts library (CSV row per part /
  oligo), then plans the optimal Gibson / Golden-Gate / SLIC /
  SLIM assembly strategy and writes the per-step protocol +
  GenBank construct files. **JAR-distributed** — no `j5` launcher
  binary on PATH; the user supplies the absolute path to `j5.jar`
  via `[bio.j5].jar` in `case.toml`. Single-binary subprocess
  shape (sister to Phase 18 BWA): the CLI is `java -jar <jar> -d
  <design_csv> -p <parts_csv> -o <output_basename> [extras...]`.
  `[bio.j5]` knobs: `jar` (absolute path to `j5.jar`; required),
  `design_csv` (j5 design CSV; required), `parts_csv` (parts list
  CSV; required), `output_basename` (filename stem the user
  expects j5 to produce; required, non-empty), `extra_args`.
  `prepare()` resolves all three input paths against the case
  directory when relative, validates each file exists on disk
  (returns `InvalidCase` with a helpful message when missing), and
  composes the `java -jar` invocation. `collect()` walks the
  workdir for `<output_basename>*.csv` (`Tabular`, "j5 assembly
  plan") and `<output_basename>*.gb` (`Native`, "j5 GenBank
  output"). Probe via `find_on_path(&["java"])` — j5's version
  comes from the jar itself, not from `java`, so we surface no
  version here; the user pins the j5 release implicitly by the jar
  they point at. Version range `1.0.0..2.0.0` (j5 has been on a
  1.x line for over a decade; upper bound 2.0 reserves room for an
  eventual major bump). `bio.j5.assemble` ribbon capability.
- `valenx-adapter-cello` — CIDAR's canonical genetic-circuit DNA
  compiler (Cello v2, BSD-3-Clause). Cello consumes a Verilog
  netlist describing the desired logic function plus a triplet of
  JSON constraint files (a user constraint file pinning the
  chassis / library, an input sensor file pinning the input
  promoters, an output device file pinning the reporter), and
  emits a fully assembled DNA construct that implements the logic
  in a living cell. The compiler runs a simulated-annealing
  optimization over the gate-assignment problem and outputs a
  Graphviz `.dot` netlist, a circuit diagram PNG, and a human-
  readable report. **JAR-distributed** — no `cello` launcher
  binary on PATH; the user supplies the absolute path to the jar
  via `[bio.cello].jar` in `case.toml`. Single-binary subprocess
  shape (sister to j5): the CLI is `java -jar <jar> -inputNetlist
  <verilog> -targetDataFile <user_constraints> -inputSensorFile
  <input_sensors> -outputDeviceFile <output_devices> -outputDir
  <output_basename> [extras...]`. `[bio.cello]` knobs: `jar`
  (absolute path to the Cello jar; required), `verilog` (`.v`
  Verilog circuit description; required), `user_constraints`
  (`.UCF` user constraints file pinning the chassis / library;
  required), `input_sensors` (`.input.json` pinning the input
  promoters; required), `output_devices` (`.output.json` pinning
  the reporter; required), `output_basename` (filename stem Cello
  uses for the output directory; required, non-empty),
  `extra_args`. `prepare()` resolves all five input paths against
  the case directory when relative, validates each file exists on
  disk, and composes the `java -jar` invocation. `collect()` walks
  the workdir for `<output_basename>*.txt` (`Log`, "Cello report"),
  `<output_basename>*.png` (`Native`, "Cello circuit diagram"),
  and `<output_basename>*.dot` (`Native`, "Cello Graphviz
  netlist"). Probe via `find_on_path(&["java"])` (same JAR-
  versioning shape as j5). Version range `2.0.0..3.0.0` (Cello v2
  is the modern Java rewrite (2020+); the v1 line was Python and
  is deprecated; upper bound 3.0 reserves room for an eventual
  major bump). `bio.cello.compile` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (pySBOL Python composition
scripts + optional starting SBOL XML, j5 design + parts CSVs + the
j5 jar, Cello Verilog netlist + UCF + input-sensor + output-device
JSONs + the Cello jar) and emit user-readable artifacts (pySBOL
XML / JSON SBOL documents, j5 CSV assembly plans + GenBank `.gb`
constructs, Cello Graphviz `.dot` netlists + circuit-diagram PNGs
+ human-readable text reports) that the unchanged
`Results.artifacts` collection model surfaces directly. A first-
class synthetic-biology canonical type — a typed SBOL-document
representation spanning pySBOL output + j5 GenBank + Cello
netlists, with parsed component / sequence / interaction graphs —
defers to a future phase along with circuit-diagram visualizers
and per-construct interactive overlays.

Three new `valenx-init` templates ship: `pysbol` with alias `sbol`
(`pysbol-compose`), `j5` (`j5-assemble`), and `cello`
(`cello-compile`). Cross-binary roundtrip test sweeps all 96
templates clean.

Adapter inventory: **100 of 101** fully live (only `occt` remains
stub-only) — **crosses the 100-adapter milestone** alongside the
Phase 5.5 MD-analysis-expansion trio.

What's not in this phase: libSBOL (the C++ implementation of SBOL —
sister to pySBOL with a different language surface; defer to Phase
33.5), iBioSim (Bridges lab's SBOL-conformant genetic-circuit
modeling + simulation environment; defer), SBOLDesigner (Anderson
lab's drag-and-drop GUI for SBOL composition; the GUI shape doesn't
fit the headless adapter pattern; defer), SynBioHub (Watson lab's
online repository for sharing SBOL designs; the web-service shape
doesn't fit the local-binary adapter pattern; defer), Tellurium
(Sauro lab's Python environment for systems / synthetic biology;
adjacent to Phase 32 systems biology; defer),
GeneticCircuitGenerator (CIDAR's combinatorial circuit-library
enumeration; defer to 33.5).

The full plan lives at
`docs/superpowers/plans/2026-04-30-synthetic-biology.md`.

### Phase 5.5 — MD analysis expansion

Sister-adapter expansion of the existing Phase 17 MDAnalysis
adapter. Round out the post-MD analysis surface that Phase 17
MDAnalysis opened with three more established open-source tools
that span the post-MD analysis tradeoff space — enhanced-sampling
collective-variable evaluation + free-energy reweighting (PLUMED,
the de-facto plug-in that wraps every major MD engine for biased-
simulation / reweighting work; LGPL-3.0; the `plumed driver`
sub-command runs PLUMED standalone over a pre-computed trajectory:
read frames, evaluate the collective variables defined in
`plumed.dat`, write COLVAR / bias / HILLS files; single-binary
subprocess shape sister to Phase 18 BWA), protein-dynamics
elastic-network / normal-mode analysis (ProDy, the canonical
Python toolkit for ENM / GNM / ANM and ensemble PCA; MIT; ships
elastic-network models, normal-mode analysis, ensemble PCA, the
NMD trajectory format consumed by VMD's NMWiz plug-in, and
integrations with the BLAST / DALI / PDB databases; Python-script
subprocess shape sister to Phase 17 Biopython), and canonical
AmberTools trajectory analysis via cpptraj's domain language
(cpptraj, the reference workhorse for `rms` / `radgyr` / `hbond`
/ `clustering` over Amber-format trajectories; GPL-3.0; reads
Amber `.prmtop` / `.parm7` topologies plus `.nc` / `.dcd` /
`.mdcrd` trajectories, runs an analysis script authored in
cpptraj's domain language, and writes results into the workdir as
`.dat` per-frame tables, `.agr` XmGrace plot data, or `.gnu`
gnuplot scripts; single-binary subprocess shape sister to PLUMED).

**This rounds out the post-MD analysis surface that Phase 17
MDAnalysis opened.** Phase 17 MDAnalysis (the de-facto Python
library for trajectory I/O + per-frame property calculation)
covered the standard trajectory-analysis surface; Phase 5.5
(PLUMED + ProDy + cpptraj) covers the corners MDAnalysis doesn't
reach: PLUMED for biased / enhanced-sampling work and free-energy
reweighting, ProDy for elastic-network / normal-mode protein
dynamics, cpptraj for the canonical AmberTools trajectory analysis
CLI workflow with the long tail of `rms` / `radgyr` / `hbond` /
`clustering` analyses. With Phase 5.5 the post-MD analysis surface
in Valenx covers all four canonical shapes — Python library API
(MDAnalysis + ProDy), enhanced-sampling plug-in CLI (PLUMED), and
AmberTools domain-language CLI (cpptraj).

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-plumed` — the de-facto enhanced-sampling and
  free-energy plug-in that wraps every major MD engine (GROMACS,
  LAMMPS, AMBER, NAMD, OpenMM) (LGPL-3.0). PLUMED defines
  collective variables (RMSD, dihedrals, distances, contact maps),
  biases (metadynamics, well-tempered metad, umbrella sampling,
  ABF), and a reweighting framework that turns biased trajectories
  back into unbiased free-energy surfaces. The `plumed driver`
  sub-command runs PLUMED standalone over a pre-computed
  trajectory: read frames from `--mf_xtc <traj>`, evaluate the
  collective variables defined in `--plumed <plumed.dat>`, write
  COLVAR / bias / HILLS files into the workdir. Single-binary
  subprocess shape (sister to Phase 18 BWA): the CLI is `plumed
  driver --plumed <plumed_dat> --mf_xtc <trajectory> --kt <kt>
  [extras...]`. `[bio.plumed]` knobs: `plumed_dat` (PLUMED input
  file describing the collective variables and bias to compute;
  required), `trajectory` (XTC trajectory; required — users running
  DCD / TRR can swap to `--mf_dcd` / `--mf_trr` via `extra_args`),
  `output_basename` (filename stem PLUMED uses for COLVAR / bias
  outputs; required, non-empty), `kt` (`f64`, > 0.0 and finite;
  PLUMED's `k_B T` in its energy units — kJ/mol by default;
  default 2.494 = room temperature 300 K; a zero or NaN `kt` would
  crash PLUMED's reweighting on the first frame), `extra_args`.
  `prepare()` resolves both paths against the case directory when
  relative, validates each file exists on disk, and composes the
  invocation. `run()` streams PLUMED's `PLUMED: PLUMED is
  starting` startup banner / periodic per-frame status / `PLUMED:
  Finishing` end-of-run sentinels into progress hints. `collect()`
  walks the workdir for `<output_basename>*.dat` (`Tabular`,
  "PLUMED COLVAR output") and `<output_basename>*.bias` (`Tabular`,
  "PLUMED bias"). Probe via `find_on_path(&["plumed"])`. Version
  range `2.9.0..3.0.0` (PLUMED 2.9 (2023) is the modern stable
  line — the `driver` sub-command, the metadynamics / OPES bias
  family, and the Python interface are all mature; upper bound 3.0
  reserves room for the long-promised next major).
  `bio.plumed.analyze` ribbon capability.
- `valenx-adapter-prody` — Bahar lab's canonical Python library
  for protein dynamics (MIT). ProDy ships elastic-network models
  (ENM / GNM / ANM), normal-mode analysis, ensemble PCA, the NMD
  trajectory format consumed by VMD's NMWiz plug-in, and
  integrations with the BLAST / DALI / PDB databases. Python-
  script subprocess shape (sister to Phase 17 Biopython): the user
  supplies a Python script referenced from `[bio.prody].script` in
  `case.toml` that imports `prody` and reads `valenx_params.json`
  for the parsed knobs. `[bio.prody]` knobs: `script` (path to
  user-supplied Python script; required), `python` (interpreter
  name; default `"python3"`), `input_pdb` (input PDB; required),
  `output_basename` (filename stem ProDy uses for ENM / mode / NMD
  outputs; required, non-empty), `num_modes` (`u32`, ≥ 1; number
  of normal modes to compute; default 20), `cutoff` (`f64`, > 0.0
  and finite; ENM contact cutoff in Å; default 15.0). `prepare()`
  stages the script + input PDB into the workdir under their
  original filenames so the script can resolve them via relative
  paths, then writes a flat `valenx_params.json` containing
  `input_pdb` (staged filename), `output_basename`, `num_modes`,
  and `cutoff`. `collect()` walks the workdir for
  `<output_basename>*.npz` (`Native`, "ProDy ENM modes"),
  `<output_basename>*.nmd` (`Native`, "ProDy NMD trajectory" — the
  NMD format consumed by VMD's NMWiz plug-in for normal-mode
  visualisation), and `<output_basename>*.csv` (`Tabular`, "ProDy
  table"). Probe via Python on PATH with an `import prody` check
  (returns `ok = true` with a warning when import fails so non-
  standard installs aren't blocked). Version range `2.4.0..3.0.0`
  (ProDy 2.x is the modern stable line; 2.4 is the floor; upper
  bound 3.0 reserves room for an eventual major bump).
  `bio.prody.analyze` ribbon capability.
- `valenx-adapter-cpptraj` — AmberTools' canonical trajectory
  analysis tool (GPL-3.0). cpptraj reads Amber `.prmtop` /
  `.parm7` topologies plus `.nc` / `.dcd` / `.mdcrd` trajectories,
  runs an analysis script authored in cpptraj's domain language
  (`trajin`, `rms`, `radgyr`, `hbond`, `volume`, `clustering`,
  ...), and writes results into the workdir as `.dat` (per-frame
  tables), `.agr` (XmGrace plot data), or `.gnu` (gnuplot
  scripts). Single-binary subprocess shape (sister to PLUMED): the
  CLI is `cpptraj -p <topology> -i <script> [extras...]`.
  `[bio.cpptraj]` knobs: `script` (`.ptraj` / `.cpptraj` analysis
  script; required), `topology` (Amber `.prmtop` / `.parm7`;
  required), `extra_args`. `prepare()` resolves both paths against
  the case directory when relative, validates each file exists on
  disk, and composes the invocation. `collect()` walks the
  workdir for `*.dat` (`Tabular`, "cpptraj analysis output"),
  `*.agr` (`Tabular`, "cpptraj XmGrace plot"), and `*.gnu` (`Log`,
  "cpptraj gnuplot script"). Probe via
  `find_on_path(&["cpptraj"])`. Version range `6.0.0..7.0.0`
  (cpptraj 6.x is the modern stable line shipped with AmberTools
  23+ (2023); upper bound 7.0 reserves room for the next major
  bump). `bio.cpptraj.analyze` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (PLUMED collective-variable
scripts + MD trajectories, ProDy Python analysis scripts + input
PDBs, cpptraj domain-language scripts + Amber topologies +
trajectories) and emit user-readable artifacts (PLUMED COLVAR
`.dat` / `.bias` tables, ProDy `.npz` ENM-mode arrays + `.nmd`
NMD-format trajectories + `.csv` analysis tables, cpptraj `.dat`
per-frame tables + `.agr` XmGrace plots + `.gnu` gnuplot scripts)
that the unchanged `Results.artifacts` collection model surfaces
directly. A first-class MD-analysis canonical type — a typed
collective-variable / normal-mode / per-frame-statistics
representation spanning all three back-ends and the existing
Phase 17 MDAnalysis adapter — defers to a future phase along with
COLVAR plotters, normal-mode visualizers, and per-statistic time-
series viewers.

Three new `valenx-init` templates ship: `plumed` (`plumed-analyze`),
`prody` (`prody-analyze`), and `cpptraj` (`cpptraj-analyze`).
Cross-binary roundtrip test sweeps all 96 templates clean.

Adapter inventory: **100 of 101** fully live (only `occt` remains
stub-only) — **crosses the 100-adapter milestone** alongside the
Phase 33 synthetic-biology trio.

What's not in this phase: AMBER `pmemd.cuda` (already part of
AmberTools but a different shape — a real MD engine rather than a
trajectory analyzer; defer to a future MD-engine phase), MDTraj
(sister Python trajectory analyzer to MDAnalysis; defer to Phase
5.6), nMOLDYN (neutron-scattering observables from MD
trajectories; defer), Lomap2 (alchemical perturbations / free-
energy GPU engines; defer to a docking / free-energy phase),
CHARMM-GUI (web-fronted CHARMM input generator; defer), HOOMD-blue
(GPU-native particle simulator; defer to a future MD-engine
phase).

The full plan lives at
`docs/superpowers/plans/2026-04-30-md-analysis.md`.

**6 newly-live adapters** — `prepare()` and `run()` filled in
end-to-end on the previously-stub-only adapters that had a
"single binary on a user-provided input" pattern:

- **SU2** (`SU2_CFD wing.cfg`) with optional mesh staging + OMP
  threading via `OMP_NUM_THREADS`.
- **Netgen** (`netgen -batchmode -geofile=… -meshfile=… [-meshsize=…]`)
  for batch CSG / BREP meshing.
- **OpenRadioss** (`engine_<arch> -i <_0001.rad> -nspmd N -nthread M`)
  for the engine phase. Starter→engine conversion stays user-side.
- **Code_Aster** (`as_run case.export`) with companion `.comm` /
  `.med` / `.py` staged automatically.
- **Meep** (Python or legacy Scheme `.ctl`) with optional `mpirun
  -np N` wrapping.
- **GROMACS** (`gmx mdrun -s <tpr> -deffnm <name> [-nt N]`) for the
  mdrun phase. Grompp stays user-side.

Each comes with progress hints on the tool's banner conventions,
output classification in `collect()`, and 8-12 tests covering
case-input parsing + prepare staging + collect classification.
Adapter inventory: 17 of 18 fully live (only `occt` remains
stub-only pending an `occt-sys` C++ FFI shim).

**Mesh quality module** — the long-deferred
"orthogonality / non-orthogonality / skewness" trio now ships:

- Equiangle skewness on every element type (Tri3 / Quad4 / Tet4 /
  Hex8 / Pyr5 / Prism6 + quadratic reduce-to-linear), exposed as
  `equiangle_skewness(t, pts) -> Option<f64>` and rolled into
  `QualityReport.max_skewness`.
- Cell-face orthogonality (FV-style: cosine angle between
  cell-centre vector and face normal) via a new
  `valenx_mesh::adjacency` module. Both face adjacency (3D) and
  edge adjacency (2D) ship.
- Aspect-ratio + skewness histograms on default CFD-prep buckets
  with cumulative-fraction helpers.
- `Mesh::recompute_quality_stats(&mut self) -> QualityReport`
  rolls up the four scalars (size, AR, skewness, orthogonality)
  into `mesh.stats` in one walk; gmsh / VTU / VTK Legacy mesh
  importers now call it instead of the bare `recompute_stats`.
- `Display` impl on `QualityReport` for CLI / log output.

**`valenx-mesh-info` CLI** for headless mesh-quality inspection.
Text + `--format json` output modes. `--check METRIC=THRESHOLD`
gates (max-skew / max-aspect / inverted / min-orthogonality) for
CI pipelines: exit 0 on pass, 4 on threshold breach. Subprocess
integration tests exercise the binary end-to-end.

**`valenx-audit tail` subcommand** for headless audit-log
inspection (mirrors the GUI command-palette caller). Default
50-entry tail, `-n` / `--lines` for explicit count, `--json`
for parseable output. Subprocess integration tests included.

**SLURM end-to-end remote-cluster submission** — `valenx-executor-slurm`
graduates from "local sbatch only" to true remote submit:

- `StagingMode::Rsync { host, remote_workdir_root }` opts into
  remote workflow.
- `submit()` rsyncs the workdir before sbatch; `poll()` /
  `cancel()` route `squeue` / `sacct` / `scancel` via ssh.
- `SlurmExecutor::fetch_results(handle)` pulls the cluster-side
  workdir back after `RunStatus::Completed`.
- `sacct` fallback so terminal state isn't optimistic-Completed
  when squeue empty.
- GPU resource declarations (`--gres`), multi-node `srun`
  wrapping when `nodes > 1` or `ntasks_per_node > 1`.
- `read_slurm_log_tail(workdir, native_id, n)` for after-the-fact
  stdout inspection.

**`valenx-init` per-adapter templates** — surfaces the newly-live
adapters as one-command starter projects. `--template su2` /
`openradioss` / `code-aster` / `netgen` / `meep` / `gromacs`
each scaffold a runnable case.toml; the text-driven adapters
(netgen / meep / su2) also drop a sample input file alongside.

**Vector + tensor field rendering** — `Field::magnitude_field()`
projects Vector{dim:1..=3} to ||v||₂ and Tensor{rows×cols} to the
Frobenius norm. The viewport's renderable filter accepts all
three field kinds; the field-list panel labels them as
"scalar" / "vector |v|" / "tensor |T|_F".

**Histograms in the GUI Quality panel** — Unicode-bar AR and
skewness histograms inside a `Quality histograms` collapsible
section, using the same data the `valenx-mesh-info --format json`
output emits.

**Adapter live_provenance rollout (18/18)** — every adapter's
`collect()` now emits a real `Provenance` block (case_hash,
mesh_hash, run_id, tool/adapter version, timestamp) instead of
the stub. The duplicated `first_workdir_match` helper across
12 adapters got DRY'd into `valenx_core::adapter_helpers`.

**OpenFOAM rhoSimpleFoam (compressible CFD)** — fourth solver in the
OpenFOAM family, completing the incompressible/compressible ×
steady/transient grid. New `Thermo` struct + optional `[flow.thermo]`
block lets users override the air-at-room-conditions defaults.
Compressible writes `thermophysicalProperties` (`hePsiThermo` +
`perfectGas` + Sutherland) + `0/T` instead of `transportProperties`,
plus the SIMPLE block grows `rhoMin`/`rhoMax` + the energy equation
in residualControl. `is_compressible()` + `is_steady()` helpers on
`SolverKind` keep the dispatching honest.

**Per-case run-history badges in the case browser** — each case row
gains a small ✓/✗/· glyph between the adapter-status dot and the
row label. Tracks the last run's success + wall time + convergence
in a `BTreeMap<String, RunHistoryEntry>` on `ValenxApp`. Hover gives
the wall time on success or "last run failed" / "not yet run".
Persistence across app restarts is a Phase 10 polish item; for now
the history lives only in memory.

**Netgen `.vol` → canonical Mesh parser** — closes the meshing UX
loop for the second mesh adapter. Netgen's `collect()` now
parses any produced `.vol` file (Tet4 / Pyr5 / Prism6 / Hex8 in 3D;
Tri3 / Quad4 in 2D) and writes `mesh.canonical.json` next to it,
mirroring the gmsh convention. The app's `on_run_finished` matches
both `gmsh` and `netgen`, so users see their generated mesh in the
viewport the moment netgen finishes.

**`valenx-validate` CLI** — structural pre-flight on a `.valenx`
project bundle (manifest, tools.lock, every case in `[cases].order`).
Text or JSON output, exit-code-driven for CI gates: 0 clean, 1
structural issue, 2 usage error, 3 IO. Doesn't depend on the full
adapter zoo, so the binary stays small and fast to build.

**`valenx-results` CLI** — headless inspector for the `results.json`
sidecar that ValenxApp writes next to every finished run. Lists
fields (with timestep counts), scalars (with sample counts), and
artifacts; provenance header includes adapter / tool versions and
the run UUID. `--format json` re-emits the file pretty-printed for
downstream tooling.

**`valenx-report` CLI** — headless HTML / Markdown / CSV exporter
for adapter `results.json` files. Wraps the existing
`write_html_report` / `write_scalars_csv` helpers and adds a
GitHub-flavoured Markdown renderer (`render_markdown_report` +
`write_markdown_report`) so the CLI can emit a PR-comment-friendly
summary alongside the self-contained HTML and the flat scalar
history CSV. Refuses to no-op silently — at least one of `--html`
/ `--markdown` / `--csv` must be supplied.

**3 new `valenx-init` templates** — `gmsh` / `lammps` /
`elmer-heat`. Closes the most-noticed gaps where a popular
adapter (gmsh, the most-used mesher; LAMMPS, classical MD;
Elmer's heat-equation mode) had no one-command starter. Aliases
include `delaunay` / `lj` / `classical-md` / `elmer` / `heat`.
Sample case directories: `box-mesh` / `lj-fluid` / `heat-cube`.

**Cross-binary roundtrip integration test** — every
`valenx-init` template (now 13) produces a project that
`valenx-validate` accepts. Locks in the contract that the two
CLIs stay in sync.

**`-V` / `--version` on every CLI** — uniform version-print
convention across init / validate / mesh-info / audit / results /
report. Each binary prints `<name> v<CARGO_PKG_VERSION>` and
exits 0. Useful when a CI recipe needs to pin tool versions.

**`valenx-init` next-step hints** — successful scaffold prints
two follow-up commands a user can run immediately
(`valenx-validate <dir>` and `$EDITOR <cases>`).

**`valenx-init --list-templates`** — quick discovery flag that
prints the full template catalogue (canonical name + brief
description + default case directory) without scaffolding
anything. Aliases: `--list-templates`, `-l`, `list-templates`.
Drift guard: a unit test asserts every catalogued name
round-trips through `Template::from_str`.

**`valenx-audit tail --since <TS>`** — ISO-8601 timestamp filter.
Drops entries older than the cutoff before the ring-buffer
truncation runs, so `tail -n 50 --since 2026-04-28T00:00:00Z`
returns "the last 50 entries on or after that instant" rather
than "the last 50 entries period, post-filtered". Library
generalisation lands as `tail_filtered(path, n, since)`.

**Windows `find_on_path` honours `PATHEXT`** — pre-fix, the
adapter probe layer only looked for `<name>.exe` on Windows,
missing `.bat` / `.cmd` / `.com` shims that conda, scoop, and
chocolatey routinely produce. New `platform_candidates(name)`
helper iterates the full `PATHEXT` and dedupes against the
caller-provided extension. Fixes a class of "tool installed but
adapter status badge is gray" bugs.

**Two new project fixtures** — `tests/fixtures/minimal.valenx/`
gains `heat-cube` (Elmer steady heat conduction with two pinned
faces) and `netgen-cylinder` (Netgen CSG primitive with sibling
`cylinder.geo`). The pre-existing `box-mesh` (gmsh.delaunay) is
also surfaced in `cases.order`. The `project_roundtrip` integration
test asserts all six cases load cleanly.

### Fixed

**`valenx-init` `cases.order` mismatch** — pre-fix, every scaffold
emitted `cases.order = ["case-1"]` regardless of template, even
though the case actually landed in `cases/cavity/`,
`cases/cantilever/`, etc. The mismatch left the loader to fall back
on filesystem scanning. Now `render_project_toml` takes the case
directory name as a parameter; both `render_project_toml_case_order_matches_dir_name`
and `scaffold_project_creates_full_skeleton` lock in the fix.

### Phase 10 release prep — release pipeline + library substrate

**`valenx-crash-reporter` crate** — in-app crash reporter. Panic
hook captures every panic, sanitises the payload (home-dir
paths, UUIDs, SHA-256 hashes redacted, message capped at 4 KiB),
and writes a JSON report to `<state_dir>/crashes/<ts>.json` before
chaining to the previously-installed hook. Reports always land
on disk; network egress is gated on the new Settings → Privacy
opt-in toggle. 19 unit tests + 1 subprocess integration test
that spawns a panic-bin and verifies the report shape.

**`valenx-first-run` crate** — first-launch wizard logic. Pure
state machine — egui rendering rides in `valenx-app::first_run`.
Builds an `EnvironmentReport` from adapter probe results, surfaces
per-OS install hints (apt / brew / winget per adapter, falling
through to upstream URLs), persists a `FirstRunDecision` to
`<state_dir>/first-run.json`. Wizard auto-opens on first install,
re-openable from the command palette afterwards. Drift guard test
asserts every workspace adapter has a Linux install hint.

**`valenx-i18n` crate** — lightweight string catalogue. Loads
`.ftl` files (key=value pairs, comment + blank-line tolerance),
supports `{ $name }` placeholder substitution, swappable for
fluent-rs in v0.2.0. Pseudo-locale (`⟦…⟧` wrap) for dev builds
to make hard-coded strings visually obvious. en-US baseline
covers ribbon / browser / status / dialog / palette / tooltip /
error namespaces. About / first-run / Settings dialogs already
wired through the catalogue.

**`valenx-a11y` crate** — WCAG 2.1 contrast computation. Relative
luminance per the W3C formula, contrast ratios always ≥ 1.0,
classification into Fail / Aa / Aaa for normal and large text.
Audit helpers + a CI gate test in `valenx-design-tokens` that
runs every documented foreground/background pair through the
audit on every build. Drift in tokens.json that drops any pair
below AA fails CI before it ships.

**Signed-installer release pipeline** — new
`.github/workflows/release.yml` builds `.deb` / `.rpm` /
`.app.zip` / `.msi` on every `v*` tag push. Signs when
APPLE_ID / WINDOWS_CERT secrets are set; falls through to
unsigned artefacts with a CI warning otherwise. Maintainer guide
at `RELEASING.md` — covers the unsigned shortcut for pre-alpha
(rustup / uv / ripgrep convention; saves $300+/yr until user
base justifies certs) plus the full signed path for v0.2.0+.

**Settings → Privacy → "Upload crash reports"** opt-in toggle.
Default OFF (privacy-preserving). Reports always write to disk;
the flag controls whether the next launch's "submit?" prompt
auto-sends or asks. Three new tests cover the default, the
serde round-trip, and forward-compat with old settings.json
files lacking the new key.

**Per-panel i18n wiring** — About dialog, first-launch wizard,
and Settings panel route every user-visible string through
`catalogue.lookup` / `format_with`. Catalogue baked into the
binary via `include_str!` — no external file resolution at
startup. Embedded en-US baseline ships ~75 keys covering every
wired panel.

### Quality

Workspace-wide hygiene: zero clippy warnings, zero rustdoc
warnings. 5 new workspace crates land alongside ~85 new tests
(19+1 valenx-crash-reporter, 13 valenx-first-run, 14 valenx-i18n,
11 valenx-a11y, 5 design-tokens contrast audit, 3 valenx-app
first_run + Settings privacy round-trip, plus 4 lib_panel-shim
tests).

The four-of-four critical-path Phase 10 lanes (B / C / D / F) are
now end-to-end complete; lane A (i18n) is library-complete with
3 panels wired and the rest pending mechanical translation.
`v0.1.0-alpha.1` is ready to tag from the current state — see
NEXT_PHASE.md for the lane status matrix.

### Phase 17 — Biology + biotech foundation

The first non-physics domain ships in Valenx. New `valenx-bio`
crate hosts canonical Sequence / Structure / Trajectory types
plus FASTA / PDB / DCD format readers (mmCIF + structured
trajectory parsing deferred to Phase 17.5). Seven first-class
adapters cover the user-visible workflows: Biopython, RDKit,
OpenMM (Python-native MD), ChimeraX (3D viz), oxDNA (CG DNA),
MDAnalysis (trajectory analysis), ColabFold (protein structure
prediction). Three new headless CLIs — `valenx-fasta`,
`valenx-pdb-info`, `valenx-blast` — round out the workflow
loop.

Adapter inventory: 24 of 25 fully live (only `occt` remains
stub-only pending an `occt-sys` C++ FFI shim). `valenx-init`
gains 7 new templates; the cross-binary roundtrip test sweeps
all 20 templates clean.

The full plan + future-phase scope-out lives at
`docs/superpowers/plans/2026-04-30-biology-foundation.md`.

### Phase 18 — Sequence alignment toolkit

The second non-physics phase ships in Valenx, extending the
biology beachhead with the most-used alignment + read-mapping
toolset. `valenx-bio` gains two new canonical types
(`FastqRecord` — sequence + per-base quality, length-validated;
`Alignment` — multiple-sequence alignment as a list of named
gapped sequences) plus new format readers for FASTQ (4-line)
and a minimal SAM-text reader (header + records sufficient for
summary inspection). BAM (binary BGZF) deferred to Phase 18.5.

Six new adapter crates land under `crates/valenx-adapters/bio/`:
`valenx-adapter-bwa` (short-read alignment via `bwa mem` /
`bwa aln`), `valenx-adapter-minimap2` (long-read + spliced +
asm-vs-asm with selectable preset), `valenx-adapter-mafft`
(`mafft --auto`-driven multiple-sequence alignment),
`valenx-adapter-muscle` (alternate MSA back-end with v3 / v5+
flag dispatch), `valenx-adapter-hmmer` (profile-HMM search
covering `hmmbuild` / `hmmsearch` / `phmmer` / `jackhmmer`),
and `valenx-adapter-samtools` (SAM/BAM utilities — `flagstat` /
`view` / `sort` / `index` / `stats`). Each wired into
`valenx-app::init_registry`.

Two new headless CLIs round out the workflow loop:
`valenx-fastq` (FASTQ inspect / validate, text + JSON output,
stdin via `-`) and `valenx-sam-info` (alignment summary —
record count, mapped / unmapped tally, reference list, average
MAPQ — from a SAM file).

Six new `valenx-init` templates ship: `bwa` (`bwa-align`),
`minimap2` (`minimap2-align`), `mafft` (`mafft-msa`), `muscle`
(`muscle-msa`), `hmmer` (`hmmer-search`), and `samtools`
(`samtools-flagstat`). Cross-binary roundtrip test sweeps all
26 templates clean.

Adapter inventory: 30 of 31 fully live (only `occt` remains
stub-only).

What's not in this phase: BAM (binary BGZF) reading and the
remaining aligners that share BWA / minimap2's adapter shape
(Bowtie2, HISAT2, STAR, MMseqs2, DIAMOND). These land in Phase
18.5; variant calling is Phase 19.

The full plan + future-phase scope-out lives at
`docs/superpowers/plans/2026-04-30-sequence-alignment-toolkit.md`.

### Phase 17.5 — Structure prediction expansion

The Phase 17 ColabFold adapter expands into the full set of
open-source sibling tools. Four new adapter crates land under
`crates/valenx-adapters/bio/`:

- `valenx-adapter-esmfold` — Meta's protein language model for
  single-sequence structure prediction. No MSA, no separate
  database step. Probes via `python -c "import esm"`.
- `valenx-adapter-openfold` — PyTorch reimplementation of
  AlphaFold 2 with the full preset family
  (`model_1` through `model_5_multimer_v3`) validated at the
  case-input layer. Optional `use_templates` knob.
- `valenx-adapter-alphafold2` — DeepMind reference AF2
  (`run_alphafold.py`). Validates `model_preset ∈ {monomer,
  monomer_ptm, multimer}` and that `max_template_date` matches
  `\d{4}-\d{2}-\d{2}`. MSA / template database stays user-provided.
- `valenx-adapter-alphafold3` — DeepMind's all-atom complex
  predictor (protein + nucleic acid + ligand). Consumes AF3's
  JSON job-spec format rather than a FASTA. The probe surfaces
  a non-commercial-weights warning into `ProbeReport.warnings`
  because AF3's model weights are released under CC-BY-NC-4.0.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type changes.** pLDDT rides the existing
`Atom.b_factor` field that the Phase 17 PDB reader already lifts
cleanly; the predicted PDBs round-trip through the unchanged
`valenx_bio::format::pdb` reader. No new format readers, no new
CLIs — `valenx-pdb-info` (Phase 17) inspects the predictions.

Four new `valenx-init` templates ship: `esmfold`
(`esmfold-predict`), `openfold` (`openfold-predict`), `alphafold2`
with alias `af2` (`alphafold2-predict`), and `alphafold3` with
alias `af3` (`alphafold3-predict`). Cross-binary roundtrip test
sweeps all 30 templates clean.

Adapter inventory: 34 of 35 fully live (only `occt` remains
stub-only).

The full plan + future-phase scope-out lives at
`docs/superpowers/plans/2026-04-30-structure-prediction-expansion.md`.

### Phase 19 — Variant calling toolkit

The next link of the genomics workflow after the Phase 18 alignment
beachhead: variant calling from aligned reads. `valenx-bio` gains
two new canonical types (`Vcf` — `##` header lines verbatim, sample
IDs from the `#CHROM` column header, list of records; `VcfRecord` —
single variant row with chrom / pos / optional id / ref / comma-split
alt / optional Phred qual / `;`-split filter / raw info / optional
format + per-sample columns) plus a minimal VCF text reader
(`valenx_bio::format::vcf`). BCF (binary) and bgzf-compressed VCF
deferred to Phase 19.5 — convert with `bcftools view` first.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-bcftools` — VCF/BCF multitool with per-action
  dispatch (`view` / `call` / `filter` / `concat`). `call` defaults
  to the multiallelic-caller variants-only flow (`-m -v -O v`),
  with the reference FASTA staged + `--threads` plumbed through.
- `valenx-adapter-gatk` — Broad Institute reference variant caller.
  Wraps `gatk --java-options "-Xmx<heap>" HaplotypeCaller` with
  reference + sorted-indexed BAM staging and an optional intervals
  (BED) restriction. Java heap validated against the conventional
  `8g` / `16g` style suffix.
- `valenx-adapter-deepvariant` — Google ML-based variant caller.
  Wraps `run_deepvariant` with a typed `model_type` ∈ `{WGS, WES,
  PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA}` and `num_shards` knob.
  Probe hint mentions both the direct binary and the Docker /
  Singularity wrapper paths — the adapter does not manage container
  runtimes, the user brings their own.

Each adapter wired into `valenx-app::init_registry`.

One new headless CLI rounds out the workflow loop: `valenx-vcf-info`
(header-line count, sample count, total records, PASS / FAIL split,
no-ALT count from a VCF file or stdin via `-`; text + JSON output;
mirrors `valenx-sam-info`).

Three new `valenx-init` templates ship: `bcftools` (`bcftools-call`),
`gatk` with alias `hc` (`gatk-haplotype`), and `deepvariant` with
alias `dv` (`deepvariant-call`). Cross-binary roundtrip test sweeps
all 33 templates clean.

Adapter inventory: 37 of 38 fully live (only `occt` remains
stub-only).

The full plan + future-phase scope-out lives at
`docs/superpowers/plans/2026-04-30-variant-calling-toolkit.md`.

### Phase 23 — Molecular viewers

Round out the visualization surface for everything Valenx's biology
stack produces. Phase 17 shipped ChimeraX as the first script-driven
molecular renderer; Phase 23 ships its three most-used siblings,
following the same shape (script-driven subprocess, headless mode,
output-in-workdir).

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-pymol` — open-source PyMOL build (the Schrödinger
  fork is proprietary; we wrap the BSD-licensed open-source line).
  Drives off `.pml` (Python-style) command files; defaults to
  `pymol -c -q <script.pml>` for headless quiet rendering.
  Collects the `.png` / `.pse` / `.cif` / `.pdb` outputs the
  script generates. `bio.pymol.render` ribbon capability.
- `valenx-adapter-vmd` — Tcl-scripted MD trajectory viewer
  (`vmd -dispdev text -e <script>`). Optional `structure` field in
  the case-input schema lets a `.pdb` / `.gro` / `.psf` get
  pre-loaded as a positional arg without the script having to
  know the path. Collects `.png` / `.tga` / `.bmp` renders, `.pdb`
  / `.gro` exported structures, and `.dat` analysis output. The
  probe pushes a license-awareness warning (containing the keyword
  `"academic"`) into `ProbeReport.warnings` because VMD ships
  under a custom non-OSS-but-free-for-academic-use license.
  `bio.vmd.render` ribbon capability.
- `valenx-adapter-igv` — `igvtools` wrapper for headless BAM /
  VCF / WIG indexing + tile generation. Per-action dispatch on
  `action ∈ {index, count, sort, tile}` — `index` writes the
  `.bai` / `.idx` sidecar next to the input (no `output` field),
  the other three actions consume an explicit `output` path.
  `count` exposes the conventional 25-bp default `window_size`.
  The companion GUI viewer is out of scope. `bio.igv.index`
  ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
viewers consume existing Phase 17 / 18 / 19 inputs and emit
user-readable artifacts (PNG / PSE / TGA renders, BAI / IDX index
sidecars) that the unchanged `Results.artifacts` collection model
surfaces directly.

Three new `valenx-init` templates ship: `pymol` (`pymol-render`),
`vmd` (`vmd-render`), and `igv` with alias `igvtools` (`igv-index`).
Cross-binary roundtrip test sweeps all 36 templates clean.

Adapter inventory: 40 of 41 fully live (only `occt` remains
stub-only).

The full plan lives at
`docs/superpowers/plans/2026-04-30-molecular-viewers.md`.

### Phase 27 — Protein design

Pair the structure-prediction adapters Valenx already ships
(Phases 17 + 17.5: ColabFold, ESMFold, OpenFold, AlphaFold 2/3)
with their de novo design counterparts. Phase 27 ships RFdiffusion
(GPU-driven protein backbone generation) and ProteinMPNN (sequence
design from backbone). Together with the prediction stack, this
gives Valenx the complete **design → predict → validate** loop in
one shell.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-rfdiffusion` — GPU-driven protein backbone
  generation. Drives off a user-supplied Python entry script that
  imports `rfdiffusion` and reads `valenx_params.json` (written
  by the adapter into the workdir) for config knobs:
  `mode ∈ {motif, binder, unconditional, partial-diffusion}`,
  `num_designs` (default 8), `diffusion_steps` (default 50,
  RFdiffusion's recommended value), `output_basename`, and the
  staged `input_pdb` filename. Sampled designs land at
  `<output_basename>_0.pdb`, `<output_basename>_1.pdb`, … and are
  surfaced as typed `Native` artifacts via
  `valenx_bio::format::pdb::read` (RFdiffusion writes pLDDT into
  the B-factor column too — same as the prediction tools).
  BSD-3-Clause licensed. `bio.rfdiffusion.design` ribbon
  capability.
- `valenx-adapter-proteinmpnn` — sequence design from a backbone
  PDB. Same Python-script-subprocess pattern as RFdiffusion;
  takes a backbone PDB and emits FASTA sequences (one per design)
  with per-residue probabilities. `valenx_params.json` carries
  `model_variant ∈ {vanilla, soluble, ca-only}`, `temperature`
  (default 0.1), `num_seq_per_target` (default 8),
  `output_basename`, and the staged `input_pdb` filename. Output
  FASTA lands at `<output_basename>.fa` and is parsed via
  `valenx_bio::format::fasta::read_str` for a richer
  `"ProteinMPNN · N sequences"` artifact label. MIT licensed.
  `bio.proteinmpnn.design` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both
adapters consume the existing Phase 17 PDB inputs and emit
user-readable artifacts (PDB backbones, FASTA sequences) that
the unchanged `Results.artifacts` collection model surfaces
directly. RFdiffusion's PDB outputs are inspectable through the
existing `valenx-pdb-info` CLI; ProteinMPNN's FASTA outputs
through the existing `valenx-fasta` CLI.

Two new `valenx-init` templates ship: `rfdiffusion` with alias
`rfd` (`rfdiffusion-design`), and `proteinmpnn` with alias
`mpnn` (`proteinmpnn-design`). Cross-binary roundtrip test
sweeps all 38 templates clean.

Adapter inventory: 42 of 43 fully live (only `occt` remains
stub-only).

The full plan lives at
`docs/superpowers/plans/2026-04-30-protein-design.md`.

### Phase 34 — Molecular docking

Add the de-facto open-source small-molecule docking pair to
Valenx's biology / chemistry stack: AutoDock Vina (the modern
single-binary docker) and AutoDock 4 (the older two-stage
`autogrid4 → autodock4` workflow still common in pharma
teaching + tutorials). Both adapters follow the established
Phase 18 BWA shape — single-action subprocess, file in / file
out, no GPU required. AutoDock 4's two-stage `prepare()`
mirrors BWA's `bwa index` → `bwa mem` pattern.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-vina` — AutoDock Vina, the modern single-
  binary small-molecule docker. Takes a receptor PDBQT +
  ligand PDBQT (prepared upstream via `prepare_receptor4.py` /
  Open Babel — adapter does NOT manage prep) and writes ranked
  binding poses to a user-named output PDBQT. Required search
  space knobs `center [x, y, z]` and `size [x, y, z]` in Å;
  `exhaustiveness` (default 8, range 1..=32) tunes search
  depth; `num_modes` (default 9) controls the number of poses
  surfaced; `energy_range` (default 3.0 kcal/mol) bounds the
  energy window between best and worst returned poses; `cpu`
  (default 0 = auto-detect) selects thread count. Output PDBQT
  collected as a `Native` artifact with label
  `"AutoDock Vina docked poses"`. Apache-2.0 licensed.
  `bio.vina.dock` ribbon capability.
- `valenx-adapter-autodock4` — AutoDock 4, the older two-stage
  docker. Stage 1: `autogrid4 -p <receptor>.gpf -l <grid_log>`
  writes the grid maps. Stage 2: `autodock4 -p <ligand>.dpf -l
  <dock_log>` reads the maps + the docking parameter file and
  runs the docking. Adapter mirrors BWA's two-stage shape:
  stage 1 runs synchronously inside `prepare()`, stage 2 lands
  as the `PreparedJob.native_command` for the shared subprocess
  runner. `skip_grid` (default `false`) lets users reuse pre-
  generated grid maps; `grid_log` (default `"autogrid4.glg"`)
  and `dock_log` (default `"autodock4.dlg"`) name the per-stage
  log files inside the workdir. Probe surfaces a warning if
  `autogrid4` is missing from PATH while `autodock4` is present
  (since the full workflow needs both binaries unless
  `skip_grid` is on). GPL-2.0-or-later licensed.
  `bio.autodock4.dock` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both
adapters consume PDBQT inputs and write PDBQT / `.dlg` /
`.glg` outputs that the unchanged `Results.artifacts`
collection model surfaces directly. PDBQT is a PDB-extension
format the existing `valenx_bio::format::pdb` reader already
inspects; Vina's docked-pose PDBQT outputs and AutoDock 4's
`.dlg` / `.pdbqt` outputs are accessible through the existing
`valenx-pdb-info` CLI.

Two new `valenx-init` templates ship: `vina` with alias
`autodock-vina` (`vina-dock`), and `autodock4` with alias
`ad4` (`autodock4-dock`). Cross-binary roundtrip test sweeps
all 40 templates clean.

Adapter inventory: 44 of 45 fully live (only `occt` remains
stub-only).

The full plan lives at
`docs/superpowers/plans/2026-04-30-docking.md`.

### Phase 24 — Cheminformatics expansion

Round out the cheminformatics surface that Phase 17's RDKit
adapter started. Phase 24 ships three sister adapters
(DeepChem + Open Babel + Avogadro 2) that together with RDKit
(already shipped) and the Phase 34 docking pair give Valenx
the complete small-molecule + cheminformatics stack. All three
follow established patterns: DeepChem mirrors RDKit's Python-
script subprocess shape, Open Babel uses BWA's single-binary
CLI shape (`obabel <in> -O <out>`), and Avogadro 2 mirrors
ChimeraX's script-driven-headless pattern. Phase 24 sits
numerically between Phase 23 and Phase 27 but ships
chronologically after Phase 34 — same convention as Phase 17.5.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-deepchem` — PyTorch-backed deep-learning
  cheminformatics; sister to RDKit's classical chemistry. Drives
  off a user-provided Python script that imports `deepchem` and
  reads `valenx_params.json` (written by the adapter into the
  workdir) for config knobs: optional inline `smiles` list
  (passed through the params file for the script to consume),
  optional `dataset_csv` (staged into the workdir), optional
  `checkpoint` model path. Output classification walks the
  workdir for `.csv` (kind `Tabular`, label
  `"DeepChem analysis output"`), `.png` (kind `Native`, label
  `"DeepChem plot"`), and `.pkl` / `.pt` (kind `Native`, label
  `"DeepChem model checkpoint"`). MIT licensed.
  `bio.deepchem.script` ribbon capability.
- `valenx-adapter-openbabel` — the de-facto open-source
  chemistry-format converter; `obabel` translates between ~120
  file formats (SMILES, MOL, MOL2, PDB, SDF, XYZ, …). Single-
  binary CLI shape: `obabel <input> -O <output> [-i
  <input_format>] [-o <output_format>] [--gen3D] [-h]
  [extras…]`. `gen_3d` (default `false`) toggles `--gen3D` for
  2D → 3D coordinate generation; `add_hydrogens` (default
  `false`) toggles the `-h` hydrogen-adding flag; explicit
  `input_format` / `output_format` overrides let users pin a
  format that the extension would mis-detect. Output collected
  as a `Native` artifact with label
  `"Open Babel converted file"`. GPL-2.0 licensed.
  `bio.openbabel.convert` ribbon capability.
- `valenx-adapter-avogadro` — Python-scriptable chemistry
  editor with a small-molecule rendering pipeline. Drives off a
  user-supplied Python script via `avogadro2 --script
  <script.py>`; an optional `structure` field (`.cml` / `.mol`
  / `.xyz` / `.pdb`) gets staged + passed as a positional arg
  so the script doesn't need to know the path. `headless`
  (default `true`) toggles `--no-gui` for batch / CI use.
  Output classification walks the workdir for `.png` (label
  `"Avogadro 2 render"`), `.cml` / `.mol` / `.xyz` (label
  `"Avogadro 2 exported structure"`). GPL-2.0-or-later
  licensed. `bio.avogadro.render` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume existing Phase 17 / 18 / 19 inputs (PDB / SDF
/ SMILES / CSV) and emit user-readable artifacts (CSV tables,
PNG renders, exported chemistry structures, PyTorch / Pickle
model checkpoints) that the unchanged `Results.artifacts`
collection model surfaces directly.

Three new `valenx-init` templates ship: `deepchem` with alias
`dc` (`deepchem-screen`), `openbabel` with alias `obabel`
(`openbabel-convert`), and `avogadro` with alias `avogadro2`
(`avogadro-render`). Cross-binary roundtrip test sweeps all 43
templates clean.

Adapter inventory: 47 of 48 fully live (only `occt` remains
stub-only).

The full plan lives at
`docs/superpowers/plans/2026-04-30-cheminformatics-expansion.md`.

### Phase 22 — Workflow managers

Add the two de-facto bioinformatics workflow orchestrators to
Valenx: Nextflow (the DSL-driven pipeline language behind
nf-core) and Snakemake (the Python-flavoured rule-based
orchestrator). Unlike the per-tool adapters Phase 17 / 18 / 19 /
23 / 24 ship, these are meta-tools — they invoke pipelines that
themselves call other bio adapters' underlying binaries. Both
follow the established Phase 18 BWA single-binary CLI shape:
probe / prepare / run / collect, output-in-workdir. Phase 22
sits numerically before Phase 23 but ships chronologically after
Phase 24 — same convention as Phase 17.5.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-nextflow` — DSL-driven pipeline orchestrator
  behind nf-core. Single-binary CLI shape: `nextflow run
  <pipeline> [-c <config>] [-profile <profile>] [-resume]
  [--<key> <value>...] [extras...]`. The `pipeline` field
  accepts a local `.nf` filename, a relative/absolute path, or a
  registry identifier like `nf-core/rnaseq`. `params` maps as
  `--<key> <value>` on the command line (values stringified;
  numeric / bool conversions happen in the Nextflow DSL).
  `profile` selects a config profile (`-profile <name>`);
  `resume` toggles `-resume` for incremental re-runs; optional
  `config` path passes through `-c <file>`. The native_command
  lives at the workdir's parent so Nextflow writes its `work/`
  and `.nextflow/` directories there. `collect()` reports the
  workdir as a `Native` artifact with label
  `"Nextflow run workdir"` and walks for `report.html` /
  `timeline.html` / `dag.svg` (Nextflow's standard observability
  outputs) surfacing them as `Log` artifacts. Apache-2.0
  licensed. `bio.nextflow.run` ribbon capability.
- `valenx-adapter-snakemake` — Python-flavoured rule-based
  pipeline orchestrator. Single-binary CLI shape: `snakemake -s
  <snakefile> --cores N [--use-conda] [-n] [--configfile <path>]
  [<targets>...] [extras...]`. `snakefile` points at the
  canonical `Snakefile` (default name; relative to the case dir
  or absolute); `targets` lists specific rules to build (empty =
  all default targets); `cores` (default 1, must be ≥ 1) sets
  `--cores N` for parallel-rule execution; `use_conda` toggles
  `--use-conda` for managed environments; `dry_run` toggles `-n`
  for plan-only inspection; optional `config_file` passes
  through `--configfile`. `collect()` reports the workdir as a
  `Native` artifact with label `"Snakemake run workdir"` and
  walks `.snakemake/log/*.log` if present, surfacing the most-
  recent log file as a `Log` artifact. MIT licensed.
  `bio.snakemake.run` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Workflow
managers are meta-orchestrators — they don't produce a single
canonical artifact of their own. The pipelines they invoke
produce whatever the underlying tools do (BAM via BWA, VCF via
bcftools, FASTA via ColabFold, …), and the unchanged
`Results.artifacts` collection model surfaces them through their
respective adapters' canonical types.

Two new `valenx-init` templates ship: `nextflow` with alias `nf`
(`nextflow-pipeline`), and `snakemake` with alias `smk`
(`snakemake-pipeline`). Cross-binary roundtrip test sweeps all 45
templates clean.

Adapter inventory: 49 of 50 fully live (only `occt` remains
stub-only).

The full plan lives at
`docs/superpowers/plans/2026-04-30-workflow-managers.md`.

### Phase 19.5 — Single-cell genomics

Open the single-cell genomics domain in Valenx with the two
most-used Python tools: Scanpy (the de-facto single-cell analysis
library — clustering, dimensionality reduction, marker discovery)
and scVI (probabilistic deep-learning models for single-cell data
via the `scvi-tools` package). Both adapters follow the established
Phase 17 Biopython / RDKit pattern — Python-script subprocess where
the user's script imports `scanpy` or `scvi` and reads
`valenx_params.json` (auto-written by the adapter) for config knobs.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-scanpy` — de-facto Python single-cell analysis
  library (BSD-3-Clause). `valenx_params.json` knobs `input_h5ad`,
  `output_h5ad`, `n_top_genes` (default 2000), `n_pcs` (default 50),
  `n_neighbors` (default 15), `resolution` (default 1.0). `collect()`
  walks `.h5ad` ("Scanpy AnnData output"), `.png` / `.pdf`
  ("Scanpy plot"), `.csv` / `.tsv` (`Tabular`, "Scanpy table").
  Probe via `find_on_path(&["python3", "python"])` then `python -c
  "import scanpy"` — surfaces an install hint when Python is on
  PATH but `scanpy` isn't importable.
- `valenx-adapter-scvi` — probabilistic deep-learning models for
  single-cell data via `scvi-tools` (BSD-3-Clause). Same Python-
  script subprocess shape as Scanpy with `valenx_params.json` knobs
  `input_h5ad`, `output_h5ad`, typed `model` ∈ `{scvi, scanvi,
  totalvi, linear-scvi}` (default `scvi`), `n_latent` (default 10),
  `n_layers` (default 2), `max_epochs` (default 400), optional
  `batch_key`. Collects under "scVI" labels.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** AnnData
reader as a canonical type (`.h5ad` is HDF5-backed and needs the
`hdf5` crate, a non-trivial C-library dep) is deferred to Phase
19.6 along with the Seurat R-runtime work.

Two new `valenx-init` templates ship: `scanpy` (`scanpy-analyse`)
and `scvi` with alias `scvi-tools` (`scvi-train`). Cross-binary
roundtrip test sweeps all 47 templates clean.

Adapter inventory: 51 of 52 fully live (only `occt` remains
stub-only).

What's not in this phase: Seurat (the dominant R single-cell
library — needs an R-runtime adapter pattern (`Rscript`-based
subprocess) that's new to Valenx; the runtime infrastructure
ships alongside the adapter), AnnData reader as a canonical type
(`.h5ad` HDF5-backed format reader; the `hdf5` crate brings a
non-trivial C-library dep), `scanpy-spatial` (niche enough to
defer), and CellxGene visualization (viewer concern, slot into
Phase 23.5). These land in Phase 19.6.

The full plan lives at
`docs/superpowers/plans/2026-04-30-single-cell-genomics.md`.

### Phase 27.5 — Protein design expansion

Sister-adapter expansion of Phase 27. Add three more open-source
protein design tools to round out the de novo design surface:
Chroma (Generate Biomedicines' joint backbone + sequence
diffusion model), ESM-IF (Meta's GVP-based inverse-folding
sequence designer — alternative to ProteinMPNN), and RFantibody
(RosettaCommons antibody-specific RFdiffusion fork). All three
follow the established Phase 27 RFdiffusion / ProteinMPNN pattern
— Python-script subprocess where the user's script imports the
relevant package and reads `valenx_params.json` (auto-written by
the adapter) for config knobs. No new infrastructure.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-chroma` — Generate Biomedicines' joint
  backbone-and-sequence diffusion model (Apache-2.0). Drives off
  a user-supplied Python entry script that imports `chroma` and
  reads `valenx_params.json` (written by the adapter into the
  workdir) for config knobs: `num_samples` (default 4),
  `length`, `temperature` (default 1.0), and `output_basename`.
  Sampled designs land at `<output_basename>_N.pdb` /
  `<output_basename>_N.fa`. `collect()` walks
  `<output_basename>*.pdb` (`Native`, "Chroma design") and
  `<output_basename>*.fa` (`Tabular`, "Chroma sequence").
  `bio.chroma.design` ribbon capability.
- `valenx-adapter-esm-if` — Meta's GVP-based inverse-folding
  sequence designer via the `fair-esm` package (MIT). Same
  Python-script-subprocess shape; the user's script imports
  `esm` (the same package as ESMFold). `valenx_params.json`
  carries `input_pdb`, `model` (default
  `esm_if1_gvp4_t16_142M_UR50` — not whitelisted because ESM-IF
  model identifiers evolve fast and the upstream package
  validates), `temperature` (default 1.0), `num_samples`
  (default 8), and `output_basename`. Output FASTA lands at
  `<output_basename>.fa` and is parsed via
  `valenx_bio::format::fasta::read` for a richer
  `"ESM-IF · N sequences"` artifact label, falling back to
  `"ESM-IF designed sequences"` on parse failure (ProteinMPNN
  pattern). `bio.esm-if.design` ribbon capability.
- `valenx-adapter-rfantibody` — RosettaCommons antibody-specific
  fork of RFdiffusion (BSD-3-Clause). Adds antibody-aware modes
  and a CDR-loop-focused sampling protocol. Same Python-script-
  subprocess shape; the user's script imports `rfantibody`.
  `valenx_params.json` carries `framework_pdb` (antibody
  framework), `target_pdb` (target antigen), `design_loops`
  (non-empty subset of `["H1", "H2", "H3", "L1", "L2", "L3"]`),
  `num_designs` (default 8), `diffusion_steps` (default 50), and
  `output_basename`. Sampled designs land at
  `<output_basename>_N.pdb` and are surfaced as `Native`
  artifacts labelled "RFantibody design".
  `bio.rfantibody.design` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume the existing Phase 17 PDB inputs and emit
user-readable artifacts (PDB backbones, FASTA sequences) that
the unchanged `Results.artifacts` collection model surfaces
directly. Chroma + RFantibody PDB outputs are inspectable
through the existing `valenx-pdb-info` CLI; ESM-IF FASTA
outputs through the existing `valenx-fasta` CLI.

Three new `valenx-init` templates ship: `chroma` (`chroma-design`),
`esm-if` with aliases `esmif` / `inverse-folding`
(`esm-if-design`), and `rfantibody` with alias `rfab`
(`rfantibody-design`). Cross-binary roundtrip test sweeps all
50 templates clean.

Adapter inventory: 54 of 55 fully live (only `occt` remains
stub-only).

What's not in this phase: framediff and Genie (alternative
diffusion-based design models — niche enough to defer further;
slot into Phase 27.6 if user demand surfaces),
AlphaFold-Multimer-Design (different shape, would need direct
AlphaFold integration), Hallucination / TrDesign-style design
(different shape, separate phase). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-protein-design-expansion.md`.

### Phase 20 — Transcript quantification

Sister-domain expansion of Phase 18 / 18.5 / 18.6 closing the
transcript-level quantification gap that Phase 18.6 explicitly
deferred. Add the two de-facto transcript-level quantification tools
to Valenx: Salmon (Rob Patro's quasi-mapping plus two-phase EM
quantifier — the GTEx / TCGA / nf-core / GENCODE reference
quantifier) and Kallisto (Lior Pachter's pseudoalignment-based
quantifier — the original "skip the alignment" approach). Both
pseudo-align reads to a transcriptome and report TPM / count per
transcript without producing intermediate SAM / BAM files; they are
shape-distinct from the Phase 18.6 RNA-seq aligners (HISAT2, STAR)
because they emit per-transcript abundance tables rather than aligned
reads. Both adapters mirror the established Phase 18 BWA two-stage
shape — single-binary CLI subprocess, file in / file out, `index →
quant` pipeline. No new infrastructure.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-salmon` — Rob Patro's quasi-mapping plus two-phase
  EM transcript-level quantification tool (GPL-3.0). Two-stage
  `salmon index → salmon quant` pipeline mirrors BWA's `bwa index →
  bwa mem` pattern, Bowtie2's `bowtie2-build → bowtie2` pattern, and
  HISAT2's `hisat2-build → hisat2` pattern. `[bio.salmon]` knobs:
  `transcriptome` (FASTA, required — used to build the index when
  `skip_index = false`), `index_dir` (the index directory salmon
  writes to / reads from; required), `reads` (1 entry single-end, 2
  entries paired-end), `output_dir` (`salmon quant -o`; required),
  `threads` (≥ 1, default 1), `skip_index` (default false; set true
  to reuse a pre-built salmon index), `libtype` (default `"A"`;
  Salmon's library-type DSL — `"A"` auto-detects orientation, `"U"`
  unstranded, `"ISF"` / `"ISR"` paired-end stranded forward / reverse,
  `"IU"` paired-end unstranded — left non-whitelisted because the
  libtype DSL has many valid combos), `extra_args`. `prepare()`
  synchronously runs `salmon index -t <transcriptome> -i <index_dir>
  -p <threads>` unless `skip_index` is set, then composes the quant
  command — single-end: `salmon quant -i <index_dir> -l <libtype>
  -p <threads> -o <output_dir> -r <reads[0]> [extras...]`; paired-
  end: `salmon quant -i <index_dir> -l <libtype> -p <threads> -o
  <output_dir> -1 <reads[0]> -2 <reads[1]> [extras...]`. `collect()`
  walks `<output_dir>` for `quant.sf` (`Tabular`, "Salmon transcript
  quantification") and `cmd_info.json` (`Log`, "Salmon command
  info"). Probe via `find_on_path(&["salmon"])`. `bio.salmon.quant`
  ribbon capability.
- `valenx-adapter-kallisto` — Lior Pachter's pseudoalignment-based
  transcript quantifier (BSD-2-Clause). Two-stage `kallisto index →
  kallisto quant` pipeline. Kallisto's index is a single `.idx` file
  (not a directory) — the only shape difference from Salmon.
  `[bio.kallisto]` knobs: `transcriptome` (FASTA, required), `index`
  (single `.idx` file path — kallisto convention), `reads` (1 or 2
  entries), `output_dir` (required), `threads` (≥ 1, default 1),
  `skip_index` (default false), `fragment_length` (optional `f64` —
  required for single-end reads only; `kallisto quant -l`),
  `fragment_sd` (optional `f64` — required for single-end reads only;
  `kallisto quant -s`), `extra_args`. Validation: when `reads.len() ==
  1`, both `fragment_length` and `fragment_sd` must be present,
  finite, and `> 0.0` (kallisto auto-detects fragment statistics from
  paired-end reads but cannot for single-end). `prepare()`
  synchronously runs `kallisto index -i <index> <transcriptome>`
  unless `skip_index` is set, then composes the quant command —
  paired-end: `kallisto quant -i <index> -o <output_dir> -t <threads>
  <reads[0]> <reads[1]> [extras...]`; single-end: `kallisto quant -i
  <index> -o <output_dir> -t <threads> --single -l <fragment_length>
  -s <fragment_sd> <reads[0]> [extras...]`. `collect()` walks
  `<output_dir>` for `abundance.tsv` (`Tabular`, "Kallisto transcript
  abundance"), `abundance.h5` (`Native`, "Kallisto HDF5 abundance"),
  and `run_info.json` (`Log`, "Kallisto run info"). Probe via
  `find_on_path(&["kallisto"])`. `bio.kallisto.quant` ribbon
  capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both adapters
consume the existing Phase 17 FASTA + Phase 18 FASTQ inputs and emit
per-transcript abundance tables (Salmon `quant.sf`, Kallisto
`abundance.tsv`) that the unchanged `Results.artifacts` collection
model surfaces directly through the `Tabular` artifact kind.
Kallisto's `abundance.h5` HDF5 sidecar needs an external tool
(`h5dump` / Python `h5py`) to inspect — the canonical H5 reader as a
Valenx CLI defers to Phase 19.6 along with the Seurat / AnnData
R-runtime work.

Two new `valenx-init` templates ship: `salmon` (`salmon-quant`) and
`kallisto` (`kallisto-quant`). Cross-binary roundtrip test sweeps all
57 templates clean.

Adapter inventory: 61 of 62 fully live (only `occt` remains
stub-only).

What's not in this phase: StringTie / Cufflinks (transcript assembly
— different workflow shape; defer to Phase 20.5 if user demand
surfaces), Tximport / DESeq2 / edgeR (downstream differential
expression — R-runtime territory; slot into Phase 19.6 along with the
Seurat R-runtime work). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-transcript-quantification.md`.

### Phase 18.6 — RNA-seq alignment

Sister-adapter expansion of Phase 18 / 18.5 closing the splice-aware
RNA-seq alignment gap that Phase 18 and Phase 18.5 explicitly
deferred. Add the two de-facto RNA-seq aligners to Valenx: HISAT2
(Daehwan Kim's graph-based splice-aware aligner — successor to
TopHat) and STAR (Alex Dobin's most-used spliced aligner, the
reference RNA-seq mapper backing GTEx / TCGA / ENCODE pipelines and
the only Phase 18.x aligner that doubles as a chromatin-
conformation tool). Both are spliced extensions of Phase 18's BWA /
Phase 18.5's Bowtie2 — they handle reads that span exon-exon
junctions, where the linear short-read aligners would soft-clip or
misalign. Both adapters mirror the established Phase 18 BWA
two-stage shape — single-binary CLI subprocess, file in / file out,
`index → align` pipeline. STAR has a heavier index step (genomic
+ splice-junction database) but the same overall shape. No new
infrastructure.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-hisat2` — Daehwan Kim's graph-based splice-aware
  RNA-seq aligner (GPL-3.0). Two-stage `hisat2-build → hisat2`
  pipeline mirrors BWA's `bwa index → bwa mem` and Bowtie2's
  `bowtie2-build → bowtie2` patterns. `[bio.hisat2]` knobs:
  `reference` (FASTA, required), `reads` (1 entry single-end, 2
  entries paired-end), `threads` (≥ 1, default 1), `skip_index`
  (default false; set true to reuse a pre-built HFM index),
  `strandness` (default `"unstranded"`; whitelist
  `["unstranded", "F", "R", "FR", "RF"]` — F/R variants match
  Illumina TruSeq stranded library prep conventions), `extra_args`.
  `prepare()` synchronously runs `hisat2-build` unless `skip_index`
  is set, then composes `hisat2 -x <ref_basename> -p <threads>
  -S out.sam [--rna-strandness <strandness>] [-U <single-read>
  | -1 <r1> -2 <r2>] [extras...]`. The `--rna-strandness` flag is
  omitted when `strandness = "unstranded"` because HISAT2 treats
  unstranded data as the default. `collect()` walks `out.sam`
  (`Tabular`, "HISAT2 aligned reads"). Probe via
  `find_on_path(&["hisat2"])`. `bio.hisat2.align` ribbon
  capability.
- `valenx-adapter-star` — Alex Dobin's spliced RNA-seq aligner
  (MIT). Note the capitalized binary name —
  `find_on_path(&["STAR"])`, not `star`. Two-stage
  `--runMode genomeGenerate → --runMode alignReads` pipeline;
  STAR's index step is heavier than BWA / Bowtie2 / HISAT2 — it
  builds a suffix-array-indexed genome under `genome_dir/` and
  optionally a splice-junction database from a GTF — but the
  adapter shape is the same. `[bio.star]` knobs: `genome_dir` (the
  pre-built STAR index directory, or where the adapter writes one
  if `skip_index = false`; required), `reference` (FASTA; required
  only when generating the index), `reads` (1 or 2 entries),
  `threads` (≥ 1, default 1), `skip_index` (default false),
  `output_type` (default `"BAM_SortedByCoordinate"`; whitelist
  `["BAM_Unsorted", "BAM_SortedByCoordinate", "SAM"]`), `sjdb_gtf`
  (optional GTF for splice-junction-database-aware indexing),
  `extra_args`. The `output_type` underscore-delimited canonical
  names map to STAR's two-arg `--outSAMtype` form:
  `"BAM_Unsorted"` → `--outSAMtype BAM Unsorted`,
  `"BAM_SortedByCoordinate"` → `--outSAMtype BAM SortedByCoordinate`,
  `"SAM"` → `--outSAMtype SAM`. `prepare()` synchronously runs
  `STAR --runMode genomeGenerate --genomeDir <genome_dir>
  --genomeFastaFiles <reference> --runThreadN N
  [--sjdbGTFfile <sjdb_gtf>] [--sjdbOverhang 100]` unless
  `skip_index` is set, then composes `STAR --runMode alignReads
  --genomeDir <genome_dir> --readFilesIn <reads...> --runThreadN N
  --outSAMtype <output_type spec> --outFileNamePrefix star_
  [extras...]`. `collect()` walks for `star_Aligned.out.{bam,sam}`
  (`Tabular` for SAM, `Native` for BAM, "STAR aligned reads") and
  `star_Log.final.out` (`Log`, "STAR alignment summary"). Validation:
  when `skip_index == false`, `reference` is required.
  `bio.star.align` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both adapters
consume the existing Phase 17 FASTA + Phase 18 FASTQ inputs and emit
SAM (HISAT2 + STAR `SAM` mode) or BAM (STAR `BAM_*` modes) outputs
that the unchanged `Results.artifacts` collection model surfaces
directly. HISAT2 SAM outputs are inspectable through the existing
`valenx-sam-info` CLI; STAR BAM outputs need the existing samtools
adapter (`samtools view`) to convert to SAM before `valenx-sam-info`
can read them; STAR's `star_Log.final.out` is plain text and surfaces
directly through the `Log` artifact kind.

Two new `valenx-init` templates ship: `hisat2` with alias `hisat`
(`hisat2-align`) and `star` (`star-align`). Cross-binary roundtrip
test sweeps all 55 templates clean.

Adapter inventory: 59 of 60 fully live (only `occt` remains
stub-only).

What's not in this phase: Salmon / Kallisto (transcript
quantification, not alignment — different shape, k-mer-based
pseudoalignment, no genome index; defer to Phase 20), TopHat
(deprecated; HISAT2 is the successor — skip), Cufflinks (assembly,
not alignment — defer to Phase 20.5 if user demand surfaces). Out
of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-rna-seq-alignment.md`.

### Phase 18.5 — Aligners expansion

Sister-adapter expansion of Phase 18. Add three more aligners
covering distinct user-facing use cases: Bowtie2 (Langmead &
Salzberg's gapped FM-index short-read aligner — alternative to
BWA), MMseqs2 (Söding lab's many-vs-many protein search +
clustering toolkit — fast alternative to BLAST), and DIAMOND
(Buchfink, Reuter & Drost's ultra-fast BLAST-protocol-compatible
protein aligner). All three follow the established Phase 18 BWA
pattern — single-binary CLI subprocess, file in / file out.
Bowtie2 mirrors BWA's two-stage `index → align` shape; MMseqs2
and DIAMOND dispatch per-action via the bcftools-style
`build_command(...) -> Result<Vec<OsString>, AdapterError>`
helper. No new infrastructure.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-bowtie2` — Langmead & Salzberg's gapped FM-
  index short-read aligner (GPL-3.0). Two-stage `bowtie2-build →
  bowtie2` pipeline mirrors BWA's `bwa index → bwa mem` pattern.
  `[bio.bowtie2]` knobs: `reference` (FASTA, required), `reads`
  (1 entry single-end, 2 entries paired-end), `threads` (≥ 1,
  default 1), `skip_index` (default false; set true to reuse a
  pre-built FM-index), `preset` (default `"sensitive"`; whitelist
  `["very-fast", "fast", "sensitive", "very-sensitive"]`),
  `extra_args`. `prepare()` synchronously runs `bowtie2-build`
  unless `skip_index` is set, then composes `bowtie2 -x
  <ref_basename> --<preset> -p <threads> -S out.sam [-U <single>
  | -1 <r1> -2 <r2>] [extras...]`. `collect()` walks `out.sam`
  (`Tabular`, "Bowtie2 aligned reads") and `.log` files (`Log`,
  "Bowtie2 log"). `bio.bowtie2.align` ribbon capability.
- `valenx-adapter-mmseqs2` — Söding lab's many-vs-many protein
  search + clustering toolkit (MIT). Per-action dispatch on
  `action ∈ ["easy-search", "easy-cluster", "easy-linsearch"]`
  via `build_command(...) -> Result<...>` (post-fix bcftools
  shape — `InvalidCase` on schema drift, never panics).
  `[bio.mmseqs2]` knobs: `query` (required), `target` (required
  for search modes, ignored for cluster), `output`, `sensitivity`
  (range `1.0..=7.5`, finite-checked, default 7.5 = max
  sensitivity), `threads` (≥ 1), `extra_args`. `easy-search` →
  `mmseqs easy-search <query> <target> <output> tmp -s
  <sensitivity> --threads N [extras...]`; `easy-linsearch` →
  same shape minus `-s`; `easy-cluster` → `mmseqs easy-cluster
  <query> <output_prefix> tmp -s <sensitivity> --threads N
  [extras...]`. `collect()` reports `output` as `Tabular` with
  per-action label ("MMseqs2 easy-search hits" /
  "MMseqs2 easy-linsearch hits" / "MMseqs2 easy-cluster output").
  Probe via `find_on_path(&["mmseqs"])` — the on-disk binary is
  just `mmseqs` (no `2` suffix). MMseqs2 versions are git-hash-
  tagged (e.g. `14-7e284`); `version_range` spans `14.0.0..17.0.0`
  to cover the current major lines. `bio.mmseqs2.search` ribbon
  capability.
- `valenx-adapter-diamond` — Buchfink, Reuter & Drost's ultra-
  fast BLAST-protocol-compatible protein aligner (GPL-3.0). Per-
  action dispatch on `action ∈ ["blastp", "blastx", "makedb"]`
  via `build_command(...)`. `[bio.diamond]` knobs: `query`,
  `database`, `output` (all required), `sensitivity` (whitelist
  `["default", "fast", "sensitive", "more-sensitive",
  "very-sensitive", "ultra-sensitive"]`), `threads` (≥ 1),
  `extra_args`. The `--default` sensitivity flag is omitted when
  `sensitivity = "default"` because DIAMOND's out-of-the-box
  default has no flag; in `makedb` mode the schema field roles
  flip — `query` is the input FASTA and `database` is the output
  DB basename (DIAMOND appends `.dmnd`). `collect()` reports the
  `output` for `blastp` / `blastx` (`Tabular`, "DIAMOND <action>
  hits"); for `makedb` reports `<database>.dmnd` (`Native`,
  "DIAMOND .dmnd database"). `bio.diamond.search` ribbon
  capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume the existing Phase 17 FASTA + Phase 18 FASTQ
inputs and emit SAM (Bowtie2) or tabular hit-table (MMseqs2 /
DIAMOND BLAST format-8) outputs that the unchanged
`Results.artifacts` collection model surfaces directly. Bowtie2
SAM outputs are inspectable through the existing `valenx-sam-info`
CLI; MMseqs2 + DIAMOND tabular hits surface via the existing
`Tabular` artifact kind.

Three new `valenx-init` templates ship: `bowtie2` with alias
`bt2` (`bowtie2-align`), `mmseqs2` with alias `mmseqs`
(`mmseqs2-search`), and `diamond` with alias `dmnd`
(`diamond-search`). Cross-binary roundtrip test sweeps all 53
templates clean.

Adapter inventory: 57 of 58 fully live (only `occt` remains
stub-only).

What's not in this phase: HISAT2 and STAR (RNA-seq-specific
splice-aware aligners — different shape; defer to Phase 18.6
RNA-seq toolkit), LAST (niche pairwise aligner), BAM (binary)
reader (needs BGZF + same scope reason as Phase 18). Out of
scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-aligners-expansion.md`.

### Phase 30 — Phylogenetics

Open the molecular phylogenetics domain in Valenx with the three
most-used maximum-likelihood tree inference tools: IQ-TREE (Bui
Quang Minh & Robert Lanfear's de-facto modern ML tree builder —
ModelFinder + UFBoot ultrafast bootstrap), RAxML-NG (Alexey Kozlov's
next-generation RAxML rewrite — successor to classical `raxmlHPC`),
and FastTree (Morgan Price's approximate-ML inference, optimized for
very large trees — sub-quadratic in alignment size). All three
follow the established Phase 18 BWA single-binary CLI pattern: input
alignment in, tree out. No two-stage index step (the alignment is
the input; the tree is the output). No new infrastructure.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-iqtree` — Bui Quang Minh & Robert Lanfear's de-
  facto modern maximum-likelihood phylogenetic tree builder
  (GPL-2.0). Single-binary subprocess shape; alignment in, tree
  out. `[bio.iqtree]` knobs: `alignment` (FASTA / PHYLIP / NEXUS /
  CLUSTAL; required), `model` (default `"MFP"` — `"TEST"` / `"MFP"`
  trigger ModelFinder's automatic model selection; otherwise pass
  e.g. `"GTR+G"` / `"WAG+I+G"` verbatim; required non-empty),
  `bootstrap` (UFBoot ultrafast bootstrap replicates; `0` disables;
  default 1000), `threads` (default `"AUTO"` — IQ-TREE's auto-
  detect, otherwise an integer count; validated against
  `^(AUTO|\d+)$`), `prefix` (output file prefix; required non-
  empty), `extra_args`. `prepare()` builds `iqtree2 -s <alignment>
  -m <model> -B <bootstrap> -T <threads> --prefix <prefix>
  [extras...]` (omitting `-B` when `bootstrap == 0`). `collect()`
  walks for `<prefix>.treefile` (`Native`, "IQ-TREE ML tree"),
  `<prefix>.iqtree` (`Log`), `<prefix>.log` (`Log`). Probe via
  `find_on_path(&["iqtree2", "iqtree"])` — newer 2.x ships as
  `iqtree2`; older 1.x as `iqtree`. `bio.iqtree.tree` ribbon
  capability.
- `valenx-adapter-raxml-ng` — Alexey Kozlov's next-generation RAxML
  rewrite (AGPL-3.0). Successor to classical `raxmlHPC`. Single-
  binary subprocess with mode dispatch. `[bio.raxml-ng]` knobs:
  `alignment` (required), `model` (substitution model — `"GTR+G"`
  / `"WAG+I+G"` / etc.; required non-empty), `mode` (`"search"`
  single-tree ML, `"all"` search + bootstrap, or `"bootstrap"`
  bootstrap-only on existing tree), `bootstrap` (replicates —
  required ≥ 1 when `mode ∈ {all, bootstrap}`, ignored otherwise),
  `threads` (≥ 1, default 1), `prefix` (required non-empty),
  `extra_args`. `prepare()` builds `raxml-ng --<mode> --msa
  <alignment> --model <model> --threads <N> --prefix <prefix>
  [--bs-trees <bootstrap> if mode in {all, bootstrap}]
  [extras...]`. `collect()` walks for `<prefix>.raxml.bestTree`
  (`Native`, "RAxML-NG ML tree"), `<prefix>.raxml.support`
  (`Native`), `<prefix>.raxml.log` (`Log`). Probe via
  `find_on_path(&["raxml-ng"])`. `bio.raxml-ng.tree` ribbon
  capability.
- `valenx-adapter-fasttree` — Morgan Price's approximate-ML
  phylogenetic inference tool (GPL-2.0). Optimized for very large
  trees: sub-quadratic in alignment size. Single-binary subprocess
  shape; writes Newick to stdout (the MAFFT-style stdout-redirect
  pattern captures stdout to the `output` path). `[bio.fasttree]`
  knobs: `alignment` (required), `output` (Newick tree path;
  required), `seq_type` (`"nt"` nucleotide or `"aa"` amino acid),
  `use_gtr` (default `true` — uses GTR for nucleotides, ignored for
  amino acid; FastTree's default is JC without this flag), `gamma`
  (gamma rate-variation model toggle; default `false`),
  `extra_args`. `prepare()` builds — nucleotide: `FastTree [-nt]
  [-gtr if use_gtr] [-gamma if gamma] <alignment>` → stdout; amino-
  acid: `FastTree [-gamma if gamma] <alignment>` → stdout.
  `collect()` reports `output` as a `Native` artifact "FastTree
  Newick tree". Probe via `find_on_path(&["FastTree", "fasttree"])`
  — binary name varies by distro. `bio.fasttree.tree` ribbon
  capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume the existing Phase 18 / 20 multiple-sequence
alignment inputs (FASTA / PHYLIP / NEXUS / CLUSTAL) and emit Newick-
format tree files that the unchanged `Results.artifacts` collection
model surfaces directly through the `Native` artifact kind. A
first-class `Tree` canonical type with a Newick reader as a Valenx
CLI defers to a future phase along with visualization integrations.

Three new `valenx-init` templates ship: `iqtree` with alias
`iqtree2` (`iqtree-build`), `raxml-ng` with alias `raxml`
(`raxml-ng-build`), and `fasttree` (`fasttree-build`). Cross-binary
roundtrip test sweeps all 60 templates clean.

Adapter inventory: 64 of 65 fully live (only `occt` remains
stub-only).

What's not in this phase: BEAST 2 / MrBayes / RevBayes (Bayesian
phylogenetics — different shape, MCMC convergence-monitoring story;
defer to Phase 30.5), PhyML (niche; defer to user demand),
ModelTest / jModelTest (model selection — workflow-orchestration
concern, slot into the workflow-manager surface), tree visualization
(FigTree, Dendroscope, TreeViewer — viewer concern, slot into
Phase 23.5 if user demand surfaces). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-phylogenetics.md`.

### Phase 28 — RNA structure

Open the RNA secondary-structure prediction domain in Valenx with
three established tools: ViennaRNA (the most-cited RNA secondary-
structure suite — `RNAfold` minimum-free-energy folding), RNAstructure
(Mathews lab's classic RNA folding toolkit — `Fold` is the flagship,
BSD-3-Clause), and NUPACK (Caltech's nucleic-acid package — academic-
license-only, surfaces a license-awareness warning à la VMD /
AlphaFold 3 / ChimeraX). ViennaRNA follows the MAFFT-style stdout-
redirect pattern (RNAfold writes to stdout); RNAstructure follows the
BWA single-binary CLI pattern with explicit `-o`-style output; NUPACK
follows the OpenMM / Scanpy Python-script-subprocess pattern (NUPACK 4
is Python-driven — the 3.x CLI is deprecated). Phase 28 sits
numerically between Phase 27.5 and Phase 30 and ships chronologically
right after the Phase 30 phylogenetics beachhead — the same
chronological-vs-numerical convention as Phase 17.5 sits numerically
between Phase 17 and Phase 18.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-viennarna` — the most-cited RNA secondary-structure
  suite (Ivo Hofacker et al., custom non-commercial / academic-use
  license). Single-binary subprocess with stdout-redirect: RNAfold
  writes the dot-bracket structure to stdout, captured to `output`
  via the MAFFT-style stdout-capture pattern. `[bio.viennarna]` knobs:
  `input` (FASTA file containing the sequence(s) to fold; required),
  `output` (dot-bracket output filename relative to workdir;
  required), `temperature` (Celsius; default 37.0; finite),
  `partition_function` (default false — toggles `-p` for partition
  function + base-pair probabilities), `allow_gu` (default true —
  `--noGU` disables GU pairs), `extra_args`. `prepare()` builds
  `RNAfold -i <input> -T <temperature> [-p] [--noGU if !allow_gu]
  [extras...]` → stdout captured to `output`. `collect()` reports
  `output` as a `Native` artifact "ViennaRNA secondary structure".
  Probe via `find_on_path(&["RNAfold"])` (capital R-N-A — that's
  ViennaRNA's binary name). **Academic-license-only** — when probe
  finds RNAfold on PATH, an `"academic"`-keyworded license-awareness
  warning gets pushed into `ProbeReport.warnings` reminding users
  that ViennaRNA is licensed for non-commercial / academic use only;
  confirm your use case complies with the upstream license before
  redistributing folds or derived data. The init aliases `vienna`
  and `rnafold` resolve to the same template. `bio.viennarna.fold`
  ribbon capability.
- `valenx-adapter-rnastructure` — Mathews lab's classic RNA folding
  toolkit (BSD-3-Clause). Single-binary subprocess (binary literally
  named `Fold`, mirroring the capital-F naming of `FastTree` and
  `STAR`). `[bio.rnastructure]` knobs: `input` (FASTA or `.seq`
  RNAstructure-native format; required), `output` (`.ct` connection-
  table file; required), `max_structures` (number of structures to
  predict; default 20; ≥ 1), `max_percent` (energy difference cutoff
  as % of MFE; default 10; in `0..=100`), `temperature` (Kelvin —
  RNAstructure's convention; default 310.15; finite, > 0.0),
  `extra_args`. `prepare()` builds `Fold <input> <output> -m
  <max_structures> -p <max_percent> -t <temperature> [extras...]`.
  `collect()` reports `output` as a `Native` artifact "RNAstructure
  connectivity table". Probe via `find_on_path(&["Fold"])`.
  `bio.rnastructure.fold` ribbon capability.
- `valenx-adapter-nupack` — Caltech's nucleic-acid package (Niles
  Pierce lab, custom academic-only license). Python-script subprocess
  shape: NUPACK 4 is Python-driven (the 3.x CLI is deprecated), so
  the user supplies a Python script that imports `nupack` and reads
  `valenx_params.json` for the config knobs. `[bio.nupack]` knobs:
  `script` (required Python file), `python` (default `python3`),
  `input` (optional FASTA / `.npc` NUPACK config), `output_basename`
  (script reads from `valenx_params.json` and writes outputs prefixed
  with this; required non-empty), `temperature` (Celsius; default
  37.0; finite), `sodium` (salt concentration in molar — NUPACK's
  `sodium` parameter; default 1.0; > 0.0 and finite). `prepare()`
  stages the script + optional input, writes `valenx_params.json`
  with the staged filename / `output_basename` / `temperature` /
  `sodium`, builds `native_command = [python, script]`. `collect()`
  walks `<output_basename>*` (`Native`, "NUPACK output") and `.npc`
  / `.json` files (`Tabular` / `Log`). Probe via
  `find_on_path(&["python3", "python"])` then `python -c "import
  nupack; print(nupack.__version__)"` — surfaces an install hint
  when Python is on PATH but `nupack` isn't importable. **Academic-
  license-only** — when probe succeeds (Python on PATH and `nupack`
  importable), an `"academic"`-keyworded license-awareness warning
  gets pushed into `ProbeReport.warnings` reminding users that
  Caltech's NUPACK license restricts redistribution + commercial
  use; confirm your use case complies before publishing analyses.
  `bio.nupack.analyze` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume the existing Phase 17 FASTA inputs and emit user-
readable artifacts (dot-bracket structures, `.ct` connection tables,
NUPACK output files) that the unchanged `Results.artifacts`
collection model surfaces directly. A first-class RNA-secondary-
structure canonical type with a dot-bracket or `.ct` reader as a
Valenx CLI defers to a future phase along with visualization
integrations.

Three new `valenx-init` templates ship: `viennarna` with aliases
`vienna` / `rnafold` (`viennarna-fold`), `rnastructure`
(`rnastructure-fold`), and `nupack` (`nupack-analyze`). Cross-
binary roundtrip test sweeps all 63 templates clean.

**Academic-license callouts.** Two of the three adapters wrap tools
that ship under non-commercial / academic-use-only licenses —
ViennaRNA (custom non-commercial license from the University of
Vienna) and NUPACK (custom Caltech academic-only license). Both
surface a license-awareness warning through their `probe()` call so
the user sees it in the registry status before they ship folds or
analyses downstream. This mirrors the existing VMD (Phase 23) /
AlphaFold 3 (Phase 17.5) / ChimeraX (Phase 27.5 expansion) probe-
warning pattern. RNAstructure ships under BSD-3-Clause and needs no
analogous callout.

Adapter inventory: 67 of 68 fully live (only `occt` remains stub-
only).

What's not in this phase: ContraFold / IPknot / ProbKnot (sub-tools
of RNA suites — niche enough to defer until user demand surfaces;
slot into Phase 28.5), LocARNA (alignment-based RNA structure
prediction — different shape, separate phase), SimRNA (3D RNA
structure — different shape, would slot alongside the Phase 17.5
protein-prediction stack), mfold / UNAFold (predecessor to
RNAstructure / superseded — defer). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-rna-structure.md`.

### Phase 25 — Quantum chemistry

Open the **first quantum-chemistry domain** in Valenx with three
established open-source tools that span the quantum-chemistry tradeoff
space — semiempirical at the fast-and-approximate end, general-purpose
HF/DFT/post-HF in the middle, massively-parallel ab initio at the high
end: Psi4 (HF/DFT/post-HF general-purpose ab initio quantum chemistry,
LGPL-3.0; Psithon-scriptable input), NWChem (Pacific Northwest National
Lab's massively-parallel ab initio + plane-wave DFT package, ECL-2.0;
its own `.nw` input format with optional `mpirun` launcher when
`mpi_procs > 1`), and xTB (Stefan Grimme's extended tight-binding
semiempirical method, LGPL-3.0; reads `.xyz` coordinates directly with
all options on the CLI). Psi4 follows the Phase 18 BWA single-binary
CLI shape with explicit `-i`/`-o` arguments; NWChem follows the same
BWA shape with optional MPI wrapping plus the MAFFT-style stdout-
redirect pattern (NWChem writes its run report to stdout, captured to
`output`); xTB follows the BWA shape with stdout captured to `xtb.log`
via the same MAFFT-style stdout-redirect pattern. Phase 25 sits
numerically between Phase 24 and Phase 27 but ships chronologically
right after Phase 28 — the same chronological-vs-numerical convention
used for Phase 17.5 / 24 / 28.

**This is the first quantum-chemistry domain to land in Valenx.** The
biology adapter family started with Phase 17 (foundation — sequence /
structure / trajectory canonical types + classical MD + cheminformatics
scripts) and expanded through Phase 17.5 / 18 / 18.5 / 18.6 / 19 /
19.5 / 20 / 22 / 23 / 24 / 27 / 27.5 / 28 / 30 / 34 to cover sequence
prediction, alignment, RNA-seq, variant calling, single-cell,
transcript quantification, workflow orchestration, molecular viewers,
cheminformatics, protein design, RNA structure, phylogenetics, and
small-molecule docking — but until Phase 25 the quantum-mechanics
surface (HF / DFT / post-HF / semiempirical methods) was absent.
Phase 25 closes that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-psi4` — open-source HF/DFT/post-HF quantum chemistry
  (Justin Turney et al., LGPL-3.0). Single-binary subprocess shape:
  Psi4 reads a Psithon (Python-scriptable) input file and writes its
  run report to a user-named output file. `[bio.psi4]` knobs: `input`
  (`.in` / `.dat` Psithon script; required), `output` (output filename
  relative to workdir; required), `threads` (default 1; ≥ 1), `memory`
  (default `"1 gb"`; matches `^\d+\s*(mb|gb|MB|GB)$` via the
  `is_valid_memory` helper), `extra_args`. `prepare()` builds
  `psi4 -i <input> -o <output> -n <threads> [-m <memory>] [extras...]`.
  The `-m` flag is only emitted when the user asked for something
  other than the documented `"1 gb"` default — passing `-m` every
  time would override Psi4's internal `"500 mb"` default with our
  fixed value even when the user didn't ask for one. `collect()`
  reports `output` as a `Log` artifact "Psi4 output" and walks the
  workdir for `.fchk` (`Native`, "Psi4 formatted checkpoint") and
  `.molden` (`Native`, "Psi4 Molden orbital data") files. Probe via
  `find_on_path(&["psi4"])`. `bio.psi4.compute` ribbon capability.
- `valenx-adapter-nwchem` — Pacific Northwest National Laboratory's
  massively-parallel ab initio + plane-wave DFT package (ECL-2.0).
  Single-binary subprocess shape with optional MPI wrapping: NWChem
  reads its own `.nw` input format and writes its run report to
  stdout — captured to a user-named output file via the MAFFT-style
  stdout-redirect pattern. `[bio.nwchem]` knobs: `input` (`.nw`
  NWChem-format script; required), `output` (output filename relative
  to workdir; required), `mpi_procs` (default 1; ≥ 1), `extra_args`.
  `prepare()` builds — serial: `nwchem [extras...] <input>`;
  parallel: `mpirun -n <mpi_procs> [extras...] nwchem <input>`. When
  `mpi_procs > 1`, prepare resolves `mpirun` via `find_on_path` and
  fails with a helpful install-hint `InvalidCase` if it's missing
  (`apt install openmpi-bin` / `apt install mpich`) rather than
  letting the child fail later with a less obvious "command not
  found". The output path is stashed in
  `PreparedJob.environment[VALENX_NWCHEM_OUTPUT]` so `run()` can
  redirect stdout to it via the MAFFT-style stdout-redirect pattern
  without re-parsing the case TOML. `collect()` reports `output` as
  a `Log` artifact "NWChem output". Probe via
  `find_on_path(&["nwchem"])`. `bio.nwchem.compute` ribbon capability.
- `valenx-adapter-xtb` — Stefan Grimme's extended tight-binding
  semiempirical quantum chemistry package (LGPL-3.0). Single-binary
  subprocess shape with stdout-redirect: xTB reads `.xyz` coordinates
  directly and writes its run report to stdout — captured to
  `xtb.log` via the MAFFT-style stdout-redirect pattern. `[bio.xtb]`
  knobs: `input` (`.xyz` geometry; required), `mode` ∈
  `{single-point, opt, ohess, hess, md}` (default `"single-point"`),
  `charge` (electron-balance `i32`; default 0), `uhf` (xTB's
  multiplicity convention — number of unpaired electrons, `u32`;
  default 0), `gfn` ∈ `{0, 1, 2}` (GFN method; default 2 — GFN2-xTB
  is the modern default), `solvent` (optional ALPB solvent name e.g.
  `"water"` / `"thf"`; `None` = gas phase), `extra_args`. `prepare()`
  builds `xtb <input> --gfn <gfn> --chrg <charge> --uhf <uhf>
  [--<mode> if mode != "single-point"] [--alpb <solvent> if Some]
  [extras...]`. `single-point` is xTB's default run type so it gets
  no flag; every other mode maps to `--<mode>`. Charge, multiplicity,
  and the GFN parameter set are always emitted so the invocation is
  unambiguous regardless of whether xTB's own defaults match ours.
  `collect()` reports `xtb.log` as a `Log` artifact "xTB stdout log"
  and walks the workdir for `xtbopt.xyz` (`Native`, "xTB optimized
  geometry"), `xtbopt.log` (`Log`), `gradient` / `hessian` files
  (`Native`). Probe via `find_on_path(&["xtb"])`. `bio.xtb.compute`
  ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied input files (Psithon `.in` / NWChem
`.nw` / xyz `.xyz` coordinates) and emit user-readable artifacts
(text output reports, formatted checkpoints, Molden orbital data,
optimized xyz geometries, gradient / hessian files) that the unchanged
`Results.artifacts` collection model surfaces directly. A first-class
quantum-chemistry canonical type — a generic energy / geometry /
orbital data type spanning all three back-ends — defers to a future
phase along with `.fchk` / Molden / `.cube` reader CLIs and
visualization integrations.

Three new `valenx-init` templates ship: `psi4` (`psi4-compute`),
`nwchem` (`nwchem-compute`), and `xtb` (`xtb-compute`). Canonical
names only — no aliases beyond the canonical names themselves.
Cross-binary roundtrip test sweeps all 66 templates clean.

Adapter inventory: 70 of 71 fully live (only `occt` remains stub-
only).

What's not in this phase: CP2K / Quantum ESPRESSO / GAMESS-US
(different shape — plane-wave / massively parallel; defer to a
sister-adapter expansion phase Phase 25.5), DFTB+ / ABINIT / Octopus
(niche; defer to user demand), ORCA (proprietary binary, free-tier-
only license; would need a separate license-mode flag like AlphaFold
3), PySCF (Python library, fits Phase 24 cheminformatics pattern;
can be added there). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-quantum-chemistry.md`.

### Phase 27.6 — EvolutionaryScale models

Complete EvolutionaryScale's open-source ESM lineup in Valenx.
Phase 17.5 + 27.5 already shipped ESMFold (single-sequence
structure prediction) and ESM-IF (GVP-based inverse-folding
sequence design). Phase 27.6 adds the remaining two: ESM3
(EvolutionaryScale's flagship generative multi-modal protein model
— joint reasoning over sequence + structure + function tracks,
modes `design` / `inverse-fold` / `scaffold` / `predict`) and ESM
Cambrian / ESMC (the smaller-faster protein representation model
for embedding-driven downstream ML, two open checkpoints
`esmc-300m` / `esmc-600m`). Both follow the established Phase 17.5
ESMFold / Phase 27.5 ESM-IF pattern — Python-script subprocess
where the user's script imports the relevant package and reads
`valenx_params.json` (auto-written by the adapter) for config
knobs. No new infrastructure.

**This phase closes out the open-source EvolutionaryScale ESM
lineup at 4 of 4 tools** — ESMFold (Phase 17.5, structure
prediction) + ESM-IF (Phase 27.5, inverse-folding sequence design)
+ ESM3 (Phase 27.6, generative multi-modal joint reasoning) + ESMC
(Phase 27.6, protein representation embeddings). All four ride the
same EvolutionaryScale `esm` Python package under the hood —
installing one installs them all, and the probe surfaces a single
unified "esm is importable" gate via the shared
`detect_esm_version` helper.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-esm3` — EvolutionaryScale's flagship generative
  multi-modal protein model (Cambrian-Open-License — open weights
  for the smaller checkpoints, non-commercial for the largest
  Forge-only variants). Drives off a user-supplied Python entry
  script that imports `esm` and reads `valenx_params.json` (written
  by the adapter into the workdir) for config knobs:
  `model_variant` ∈ `{open, open-multimer, small}` (the open-
  weight variants — larger Forge-only variants are not supported
  by this adapter), `mode` ∈
  `{design, inverse-fold, scaffold, predict}`, `num_samples`
  (default 4, ≥ 1), `input_pdb` (optional PDB — required for
  `inverse-fold` and `scaffold`), `input_fasta` (optional FASTA —
  required for `predict`), `temperature` (default 1.0, > 0 and
  finite), `output_basename`. Sampled outputs land at
  `<output_basename>*.pdb` (generated structures) and
  `<output_basename>*.fa` (generated sequences). `collect()` walks
  `<output_basename>*.pdb` (`Native`, "ESM3 generated structure")
  and `<output_basename>*.fa` (`Tabular`, "ESM3 generated
  sequence"). Probe via `find_on_path(&["python3", "python"])`
  then `python -c "import esm; print(esm.__version__)"` —
  surfaces an install hint when Python is on PATH but `esm` isn't
  importable. `bio.esm3.generate` ribbon capability.
- `valenx-adapter-esmc` — EvolutionaryScale's open-weight protein
  representation model (Cambrian-Open-License). Same Python-script-
  subprocess shape; the user's script imports `esm` (the same
  package as ESMFold / ESM-IF / ESM3). `valenx_params.json`
  carries `input_fasta` (required), `model_variant` ∈
  `{esmc-300m, esmc-600m}` (the two open release sizes — 300M fits
  on a consumer GPU, 600M for larger / better representations),
  `pooling` ∈ `{per-residue, mean}`, and `output_basename`.
  Embedding tables land at
  `<output_basename>.{npy,npz,parquet}`. `collect()` walks for
  those (`Tabular`, "ESMC embeddings"). Probe identical to ESM3
  (`import esm`). `bio.esmc.embed` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both adapters
consume the existing Phase 17 PDB / FASTA inputs and emit user-
readable artifacts (PDB backbones, FASTA sequences, NumPy `.npy` /
`.npz` and Parquet embedding tables) that the unchanged
`Results.artifacts` collection model surfaces directly. ESM3 PDB
outputs are inspectable through the existing `valenx-pdb-info`
CLI; ESM3 FASTA outputs through the existing `valenx-fasta` CLI;
ESMC's embedding sidecars feed straight into the user's downstream
Python pipeline. A canonical embeddings CLI defers to a future
phase along with HDF5 / Arrow reader work.

Two new `valenx-init` templates ship: `esm3` (`esm3-generate`) and
`esmc` with alias `esm-cambrian` (`esmc-embed`). Cross-binary
roundtrip test sweeps all 68 templates clean.

Adapter inventory: 72 of 73 fully live (only `occt` remains stub-
only).

What's not in this phase: ESM3 commercial Forge API (would need an
HTTP-API client, not subprocess; out of scope for the open-weights
adapter), multi-chain ESM3 design (recently released — works
through the same Python script entry; the schema is forward-
compatible), other EvolutionaryScale Biohub items (CELL×GENE,
CryoET Data Portal — different shape; tracked as future work).
Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-evolutionaryscale.md`.

### Phase 32 — Systems biology

Open the **first systems-biology / multiscale modeling domain** in
Valenx with three established open-source tools that span the systems-
biology tradeoff space — biochemical pathway / ODE simulation at the
deterministic end (COPASI), rule-based combinatorial signaling
networks in the middle (BioNetGen), and agent-based multicellular
tissue simulation at the spatial / multiscale end (PhysiCell): COPASI
(the COmplex PAthway SImulator — de-facto biochemical pathway / ODE-
based systems-biology suite descended from the Gepasi lineage,
Artistic-2.0; reads `.cps` native or SBML `.xml`), BioNetGen (rule-
based modeling language + tool suite for combinatorially-complex
signaling networks, MIT; Perl driver `BNG2.pl` reads `.bngl` rule-
based models and emits `<basename>.net` / `<basename>.gdat` /
`<basename>.cdat` outputs), and PhysiCell (Paul Macklin's agent-
based, off-lattice multicellular simulator — tens to hundreds of
thousands of individual cells coupled to a reaction-diffusion
microenvironment for substrates like oxygen and drugs; canonical use
case is tumour growth + immunology; BSD-3-Clause; models compile per-
project to a project-specific C++ binary). All three follow the
established Phase 18 BWA single-binary CLI pattern: model file in,
results out.

**This is the first systems-biology / multiscale modeling domain to
land in Valenx.** The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 30 / 34 to cover sequence prediction, alignment,
RNA-seq, variant calling, single-cell, transcript quantification,
workflow orchestration, molecular viewers, cheminformatics, quantum
chemistry, protein design, EvolutionaryScale models, RNA structure,
phylogenetics, and small-molecule docking — but until Phase 32 the
systems-biology / multiscale-modeling surface (biochemical pathway
ODE simulation, rule-based signaling networks, agent-based
multicellular tissue simulation) was absent. Phase 32 closes that
gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-copasi` — the COmplex PAthway SImulator
  (Artistic-2.0). The de-facto desktop suite for biochemical pathway
  and ODE-based systems-biology models, descended from the Gepasi
  lineage. Single-binary subprocess shape: COPASI's headless CLI is
  `CopasiSE` (capital `C-S-E`, "Self-Executing"), a task runner that
  reads a COPASI native `.cps` archive (or an SBML `.xml`) and
  executes the simulation / scan / fitting tasks defined inside.
  `[bio.copasi]` knobs: `model` (`.cps` / `.sbml` / `.xml`; required),
  `report` (optional `--save <report>` target so the run output lands
  at a known path collect() can find without walking), `run_all`
  (default `false`; when `true` adds `--scheduled`, executing every
  task in the file rather than just the primary one), `extra_args`.
  `prepare()` composes `CopasiSE [--save <report>] <model>
  [--scheduled] [extras...]`. `collect()` reports the explicit
  `report` path as `Tabular` ("COPASI report") when supplied, else
  walks the workdir top-level for `.csv` / `.txt` files (COPASI's
  tabular outputs). Probe via `find_on_path(&["CopasiSE"])`. Version
  range `4.40.0..5.0.0` (4.x is the long-running stable line; 4.40
  is a recent floor that ships SBML L3v2 + the task scheduler every
  Phase 32 model relies on). `bio.copasi.simulate` ribbon
  capability.
- `valenx-adapter-bionetgen` — the rule-based modeling language and
  tool suite for combinatorially-complex signaling networks (MIT).
  The user writes BNGL (BioNetGen Language) files describing
  molecular species, sites, and reaction *rules*, and BioNetGen
  expands the rules into the underlying reaction network and
  (optionally) integrates it deterministically (ODE) or
  stochastically (SSA). Single-binary subprocess shape: `BNG2.pl` is
  the canonical Perl driver. `[bio.bionetgen]` knobs: `model`
  (`.bngl`; required), `output_basename` (required — becomes the
  `-o` prefix every output file inherits so collect() walks
  deterministically), `generate_only` (default `false`; when `true`
  adds `--no-execute`, skipping simulate / scan / fitting actions and
  emitting just the expanded reaction network), `extra_args`.
  `prepare()` builds `BNG2.pl [--no-execute if generate_only] -o
  <output_basename> <model> [extras...]`. `collect()` walks the
  workdir top-level for `<output_basename>*.net` (`Native`,
  "BioNetGen reaction network"), `<output_basename>*.gdat`
  (`Tabular`, "BioNetGen species trajectories"), and
  `<output_basename>*.cdat` (`Tabular`, "BioNetGen concentrations")
  — `parameter_scan` per-trial variants share the basename prefix
  (e.g. `<basename>_001.gdat`) so the prefix-restricted walk picks
  them up too. Probe via `find_on_path(&["BNG2.pl"])`. Version range
  `2.8.0..3.0.0`. The `valenx-init` template ships with the alias
  `bng` alongside the canonical `bionetgen`.
  `bio.bionetgen.simulate` ribbon capability.
- `valenx-adapter-physicell` — Paul Macklin's agent-based, off-
  lattice multicellular simulator (BSD-3-Clause). PhysiCell models
  tens to hundreds of thousands of individual cells (each an agent
  with state, mechanics, secretion, and phenotype) coupled to a
  reaction-diffusion microenvironment for substrates like oxygen and
  drugs. The canonical use case is tumour growth and immunology.
  Unlike a typical CLI tool, PhysiCell models compile to a project-
  specific C++ executable: the user clones the framework, edits the
  project's `custom_modules/` source, runs `make`, and ends up with
  e.g. `./project` next to the project directory. The adapter
  therefore takes both a `binary` path and the run-time XML
  configuration. `[bio.physicell]` knobs: `binary` (the per-project
  compiled executable; required), `config` (the `.xml` settings file
  PhysiCell binaries accept as a positional argument; required),
  `extra_args`. `prepare()` validates `binary` and `config` exist on
  disk (returns `InvalidCase` with a helpful "PhysiCell models
  compile per-project — clone the framework, edit the project's
  `custom_modules/` source, run `make`, and point this field at the
  resulting executable." message if missing), then builds `<binary>
  <config> [extras...]`. `collect()` walks `output/`. PhysiCell drops
  a stack of per-snapshot files there: `output<N>.xml` (manifest),
  `output<N>_*.mat` (cell + microenvironment state in MATLAB v4
  binary), and optional `*.csv` scalar summaries — typed `Native`
  ("PhysiCell tissue snapshot") for `.xml` / `.mat` and `Tabular`
  ("PhysiCell scalar table") for `.csv`. Probe via
  `find_on_path(&["physicell"])` returns `ok = true` either way (most
  installs won't have a generic `physicell` binary on PATH — the
  per-project build pattern means there isn't a canonical one) and
  attaches a warning that the real validation happens in `prepare()`
  against the user's `binary` field. Version range `1.13.0..2.0.0`.
  `bio.physicell.simulate` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (COPASI `.cps` / SBML `.xml`
archives, BioNetGen `.bngl` rule-based models, PhysiCell per-project
compiled binaries + XML config) and emit user-readable artifacts
(CSV / TXT tabular reports, reaction-network `.net` files, species-
trajectory `.gdat` / concentration `.cdat` tabular files, `.xml` /
`.mat` per-snapshot tissue state, per-cell scalar `.csv` summaries)
that the unchanged `Results.artifacts` collection model surfaces
directly. A first-class systems-biology canonical type — a generic
SBML / BNGL / per-cell state type spanning all three back-ends —
defers to a future phase along with SBML readers and tissue-snapshot
visualizers.

Three new `valenx-init` templates ship: `copasi` (`copasi-simulate`),
`bionetgen` with alias `bng` (`bionetgen-simulate`), and `physicell`
(`physicell-simulate`). Cross-binary roundtrip test sweeps all 71
templates clean.

Adapter inventory: 75 of 76 fully live (only `occt` remains stub-
only).

What's not in this phase: Tellurium / libRoadRunner (Python library,
fits the Biopython subprocess pattern; defer to Phase 32.5), VCell
(Java GUI app; `vcell-cli` exists but workflow is heavy; defer),
E-Cell / Morpheus / CompuCell3D (niche; defer to Phase 32.5),
Smoldyn / MCell (particle-based simulators; defer), StochPy /
libSBML / PySB (Python libraries; future Phase 32.5). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-systems-biology.md`.

### Phase 39 — DNA structural geometry

Open the **first DNA structural-geometry domain** in Valenx with
three established open-source tools that span the structural-
geometry tradeoff space — base-pair / base-step parameter
calculation (X3DNA), helical-axis curvature analysis with groove
geometry (Curves+), and structural-feature annotation as a single
machine-readable JSON summary (DSSR, the modern Python-fronted
X3DNA-family tool): X3DNA (Wilma Olson and Xiang-Jun Lu's
reference toolkit for DNA / RNA structural-geometry analysis —
de-facto reference for canonical helical-step parameters: twist,
roll, tilt, slide, shift, rise plus per-base intra-pair
parameters: buckle, propeller, opening, shear, stretch, stagger;
custom X3DNA-License — academic / non-commercial use only;
single-binary subprocess shape sister to Phase 18 BWA), Curves+
(Richard Lavery's reference toolkit for DNA helical-axis analysis
— fits a curvilinear helical axis through a nucleic-acid
structure and reports per-base axis-curvature, base-pair
parameters relative to that axis, and a `.cda` file describing
the axis itself for downstream visualisation; the canonical tool
for "is this DNA bent, and if so, how" questions in protein-DNA /
drug-DNA structural studies; custom Curves-License — academic /
non-commercial use only; single-binary subprocess shape with
stdin-piped namelist parameters sister to Phase 36 CTFFIND), and
DSSR (Dissecting the Spatial Structure of RNA / DNA — the modern
Python-fronted X3DNA-family tool that reads a nucleic-acid PDB and
emits a single JSON file enumerating every detected feature: base
pairs (Watson-Crick, Hoogsteen, sugar-edge, ...), multiplets,
double helices, stems, hairpin / internal / junction loops,
kissing loops, A-minor motifs, ribose zippers, pseudoknots, and
more; the standard machine-readable feature-extraction step in
modern RNA-structure pipelines; custom DSSR-License — academic /
non-commercial use only; single-binary subprocess shape sister to
X3DNA). All three adapters surface their respective non-OSS
academic / non-commercial-use terms accurately via `tool_license =
"X3DNA-License"` / `"Curves-License"` / `"DSSR-License"` and emit
a probe warning whenever each binary is detected, with the literal
`"academic"` string baked into the warning as a stable anchor for
license-aware filters and tests.

**This is the first DNA structural-geometry domain to land in
Valenx.** The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27
/ 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 34 / 35 / 36 / 38
to cover sequence prediction, alignment, RNA-seq, variant calling,
single-cell, transcript quantification, workflow orchestration,
molecular viewers, cheminformatics, quantum chemistry, protein
design, EvolutionaryScale models, RNA structure, population
genetics, phylogenetics, Bayesian phylogenetics, sequencing read
simulation, systems biology, small-molecule docking, CRISPR
design, cryo-EM reconstruction, and Rosetta protein modeling — but
until Phase 39 the DNA structural-geometry surface (canonical
helical parameters, helical-axis curvature analysis, machine-
readable structural-feature annotation) was absent. Phase 39
closes that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-x3dna` — Wilma Olson and Xiang-Jun Lu's
  reference toolkit for DNA / RNA structural-geometry analysis
  (custom X3DNA-License — academic / non-commercial use only).
  Reads a nucleic-acid PDB, identifies base pairs, and computes
  the canonical helical-step parameters (twist, roll, tilt,
  slide, shift, rise) plus per-base intra-pair parameters
  (buckle, propeller, opening, shear, stretch, stagger). It is
  the workhorse behind structural-bioinformatics pipelines that
  need quantitative DNA geometry — bending studies, drug-DNA /
  protein-DNA complex analysis, RNA tertiary-structure
  annotation. Single-binary subprocess shape (sister to Phase 18
  BWA): the CLI is `analyze <input_pdb> [extras...]`. `analyze`
  is positional-only — it derives every output filename from the
  input basename, so the adapter just hands it the PDB and any
  user-supplied extras. `[bio.x3dna]` knobs: `input_pdb` (input
  PDB; required), `output_basename` (filename stem the user
  expects X3DNA to produce — surfaced here so collect() can label
  artefacts uniformly without scraping `analyze`'s filename
  heuristics; required, non-empty), `extra_args`. `prepare()`
  resolves the input PDB against the case directory when relative,
  validates it exists on disk (returns `InvalidCase` with a
  helpful message when missing), and composes the positional
  invocation. `collect()` walks the workdir for
  `<output_basename>*.par` (`Tabular`, "X3DNA base-step
  parameters") and `*.out` (`Log`, the per-run log `analyze`
  writes alongside). Probe via `find_on_path(&["analyze"])`
  (X3DNA's main analysis binary is literally named `analyze`).
  **Academic-license-only** — probe always pushes an
  `"academic"`-keyworded warning into `ProbeReport.warnings`
  whenever the binary is detected, and `tool_license` surfaces as
  `"X3DNA-License"` rather than mislabeling the custom X3DNA
  terms as a recognised SPDX identifier. Version range
  `2.4.0..3.0.0` (X3DNA 2.4 (2020) is the modern stable release
  and the floor we test against; upper bound 3.0 reserves room
  for an eventual major bump). `bio.x3dna.analyze` ribbon
  capability.
- `valenx-adapter-curves` — Richard Lavery's reference toolkit
  for DNA helical-axis analysis (custom Curves-License —
  academic / non-commercial use only). Fits a curvilinear helical
  axis through a nucleic-acid structure and reports per-base
  axis-curvature, base-pair parameters relative to that axis,
  and a `.cda` file describing the axis itself for downstream
  visualisation. It is the canonical tool for "is this DNA bent,
  and if so, how" questions in protein-DNA / drug-DNA structural
  studies. Single-binary subprocess shape with stdin-piped
  parameters (sister to Phase 36 CTFFIND): Curves+ takes its
  parameters as a Fortran-style `&inp ... &end` namelist block on
  stdin followed by strand / axis residue cards. The adapter
  authors that block at `prepare()` time and pipes it into `Cur+`'s
  stdin at `run()` time via `Stdio::from(file)` — the shared
  `subprocess::run` helper closes stdin which makes Curves+ read
  EOF before parsing its first parameter and exit, so the custom
  `run()` opens the parameters file with `File::open()` and hands
  its FD to the child via `Stdio::from(file)` (the custom run
  path mirrors the MAFFT stdout-redirect pattern but for stdin,
  same shape Phase 36 CTFFIND uses). `[bio.curves]` knobs:
  `input_pdb` (input PDB; required), `output_basename` (filename
  stem Curves+ uses for outputs — `<basename>.lis`,
  `<basename>.cda`, etc.; required, non-empty), `first_residue`
  (`u32` — first inclusive residue index in the strand to
  analyse; required), `last_residue` (`u32`, ≥ `first_residue` —
  a reverse range is rejected up front with a helpful message;
  required), `extra_args`. `prepare()` resolves the input PDB
  against the case directory when relative, validates it exists
  on disk, writes `curves_params.txt` containing the namelist
  body + residue-range cards, stashes the filename under the
  sentinel env var `VALENX_CURVES_PARAMS_FILE`, and the custom
  `run()` recovers the filename, strips the sentinel from the
  env table so Curves+ doesn't see it, opens the params file, and
  pipes its contents into the child. `collect()` walks the
  workdir for `<output_basename>*.lis` (`Log`, "Curves+ helical
  analysis") and `<output_basename>*.cda` (`Tabular`, "Curves+
  axis curve data"). Probe via `find_on_path(&["Cur+"])` (the
  binary name uses a literal `+`). **Academic-license-only** —
  probe always pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever the binary is detected, and
  `tool_license` surfaces as `"Curves-License"` rather than
  mislabeling the custom Curves+ terms as a recognised SPDX
  identifier. Version range `2.0.0..3.0.0` (Curves+ 2.x is the
  modern stable line; 2.0 is the floor; upper bound 3.0 reserves
  room for an eventual major bump). `bio.curves.analyze` ribbon
  capability.
- `valenx-adapter-dssr` — Dissecting the Spatial Structure of
  RNA / DNA, the modern Python-fronted X3DNA-family tool (custom
  DSSR-License — academic / non-commercial use only). Reads a
  nucleic-acid PDB and emits a single JSON file enumerating every
  detected feature: base pairs (Watson-Crick, Hoogsteen, sugar-
  edge, ...), multiplets, double helices, stems, hairpin /
  internal / junction loops, kissing loops, A-minor motifs,
  ribose zippers, pseudoknots, splayed-apart conformations, and
  more. It is the standard machine-readable feature-extraction
  step in modern RNA-structure pipelines. Single-binary
  subprocess shape (sister to X3DNA): the CLI is `x3dna-dssr
  -i=<input_pdb> -o=<output_json> --json [extras...]` (DSSR uses
  `key=value` flag form on its short-form options — no space
  between flag and value). `[bio.dssr]` knobs: `input_pdb` (input
  PDB; required), `output_json` (output JSON path; required),
  `extra_args`. `prepare()` resolves the input PDB against the
  case directory when relative, scopes the output JSON path to
  the workdir when relative, validates the input exists on disk,
  and composes the flag-form invocation. `collect()` reports the
  configured `output_json` file as a single `Tabular` artifact
  ("DSSR analysis (JSON)") — DSSR's JSON is the canonical
  machine-readable summary; tagged `Tabular` rather than `Native`
  so downstream serdes can key off a consistent kind. Probe via
  `find_on_path(&["x3dna-dssr"])`. **Academic-license-only** —
  probe always pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever the binary is detected, and
  `tool_license` surfaces as `"DSSR-License"` rather than
  mislabeling the inherited X3DNA-family terms as a recognised
  SPDX identifier. Version range `2.0.0..3.0.0` (DSSR 2.x is the
  modern stable line that ships with X3DNA 2.4+; upper bound 3.0
  reserves room for an eventual major bump). `bio.dssr.analyze`
  ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (X3DNA / Curves+ / DSSR all
take nucleic-acid PDBs, plus the Curves+ residue-range knobs) and
emit user-readable artifacts (X3DNA `.par` base-step parameter
tables and `.out` per-run logs, Curves+ `.lis` helical-analysis
logs and `.cda` axis-curve data files, DSSR JSON structural-
feature summaries) that the unchanged `Results.artifacts`
collection model surfaces directly. The existing
`valenx_bio::format::pdb` reader inspects collected PDB inputs
for chain / residue / atom counts. A first-class DNA-geometry
canonical type — a typed helical-parameter representation
spanning all three back-ends, with parsed per-step parameter
tables and a typed structural-feature summary — defers to a
future phase along with helical-axis visualizers and per-feature
interactive overlays.

Three new `valenx-init` templates ship: `x3dna` with alias `3dna`
(`x3dna-analyze`), `curves` with alias `curves+` (`curves-
analyze`; carries an inline academic-license note in the
scaffolded `case.toml`), and `dssr` (`dssr-analyze`; carries an
inline academic-license note in the scaffolded `case.toml`).
Cross-binary roundtrip test sweeps all 90 templates clean.

Adapter inventory: 94 of 95 fully live (only `occt` remains stub-
only).

What's not in this phase: web3DNA (web-fronted X3DNA companion;
defer to Phase 39.5), 3D-DART (DNA-axis transformation tool;
defer), MC-Sym / MC-Annotate (Major / Cedergren group's RNA
structural annotation; defer to 39.5), Madbend (DNA-bending
statistical analysis; defer), DNAtools (Berman lab's DNA-
conformation toolkit; defer).

The full plan lives at
`docs/superpowers/plans/2026-04-30-dna-geometry.md`.

### Phase 30.5 — Bayesian phylogenetics

Sister-domain expansion of Phase 30. Round out the molecular
phylogenetics surface with the two de-facto Bayesian phylogenetic
inference engines — BEAST 2 (Bayesian Evolutionary Analysis by
Sampling Trees v2 — the cross-platform XML-driven MCMC framework
with a sprawling package ecosystem covering tip-dated trees,
relaxed molecular clocks, coalescent demographic models, birth-
death speciation models, and the BDSKY / MASCOT / BEASTling /
StarBEAST3 universe of extensions; LGPL-2.1; single-binary
subprocess shape sister to Phase 18 BWA), and MrBayes (the long-
standing Bayesian MCMC tree inference engine that remains the
de-facto choice alongside BEAST 2 for posterior tree sampling
across nucleotide / amino-acid / morphological datasets, with its
own NEXUS-embedded model-and-mcmc command language and built-in
Metropolis-coupled MCMC ("MC^3") chain swapping; GPL-3.0; single-
binary subprocess shape sister to BEAST 2). Both adapters share
the established Phase 18 BWA single-binary CLI pattern: a
user-authored model description (BEAST 2 XML or MrBayes NEXUS
file) in, posterior tree + parameter samples out — the same
shape the Phase 30 ML tools use, with the inputs swapped from
multiple-sequence alignments to MCMC model files.

**This rounds out the molecular phylogenetics surface that Phase
30 opened from the maximum-likelihood side.** Phase 30 (IQ-TREE +
RAxML-NG + FastTree) covered the ML tradeoff space — modern ML
with ModelFinder + UFBoot bootstrap, the next-generation RAxML
rewrite, and approximate-ML for very large trees. Phase 30.5
(BEAST 2 + MrBayes) covers the Bayesian MCMC side: BEAST 2 for
time-calibrated trees with relaxed clocks, demographic priors,
and the sprawling package ecosystem; MrBayes for the historic
NEXUS-embedded MCMC workhorse with built-in MC^3 chain swapping.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-beast2` — the cross-platform Bayesian
  Evolutionary Analysis by Sampling Trees v2 engine (LGPL-2.1).
  BEAST 2 is the canonical Bayesian MCMC framework for time-
  calibrated phylogenetics: tip-dated trees, relaxed molecular
  clocks, coalescent demographic models, birth-death speciation
  models, and the ever-growing universe of BEAST 2 packages
  (BDSKY, MASCOT, BEASTling, StarBEAST3, ...). It complements the
  maximum-likelihood Phase 30 tools (IQ-TREE, RAxML-NG, FastTree)
  with a full posterior over tree topologies and parameters.
  Single-binary subprocess shape (sister to Phase 18 BWA): the CLI
  is `beast [-seed <N>] -threads <N> [-overwrite] <xml>
  [extras...]`. The user authors / generates the model XML
  (typically through BEAUti) and references it from
  `[bio.beast2].xml` — the adapter doesn't generate XML.
  `[bio.beast2]` knobs: `xml` (BEAUti-generated XML model file;
  required), `seed` (optional `u64`; passed via `beast -seed <N>`
  when present, otherwise BEAST picks its own seed and prints it
  on the run banner), `threads` (`u32`, ≥ 1, default 1; maps to
  `beast -threads N` for tree-likelihood evaluation parallelism),
  `overwrite` (default `false`; toggles `-overwrite` so an
  existing output set from a previous run is replaced rather than
  triggering a fail), `extra_args`. `prepare()` resolves the XML
  against the case directory when relative, validates it exists on
  disk (returns `InvalidCase` with a helpful message when
  missing), and composes the invocation with `seed` injected
  before `-threads` and the XML positional last so BEAST treats it
  as the model file rather than the value of an earlier flag.
  `run()` streams BEAST's `Random number seed` / `BEAST v2`
  startup banner / periodic `Sample` / `posterior` chain-status
  lines / `End likelihood` / `Total calculation time` end-of-run
  sentinels into progress hints. `collect()` walks the workdir
  for `*.log` (`Log`, "BEAST 2 trace log" — the parameter trace
  Tracer reads) and `*.trees` (`Native`, "BEAST 2 sampled trees"
  — the sampled tree posterior TreeAnnotator / DensiTree
  consumes); the adapter doesn't try to predict the exact
  filenames since BEAST writes whatever the XML's
  `<log fileName="...">` sites configure. Probe via
  `find_on_path(&["beast"])`; the generic version detector tries
  the conventional `--version` and BEAST's own `-version` form.
  Version range `2.7.0..3.0.0` (the modern stable line is the
  2.7.x series from 2022+ that introduced modern threading + the
  package manager; upper bound 3.0 reserves room for an eventual
  major bump). `bio.beast2.mcmc` ribbon capability.
- `valenx-adapter-mrbayes` — the long-standing Bayesian MCMC
  phylogenetic inference engine (GPL-3.0). MrBayes is the historic
  workhorse for Bayesian phylogenetics: alongside BEAST 2 it
  remains the de-facto choice for posterior tree sampling across
  nucleotide / amino-acid / morphological datasets, with its own
  NEXUS-embedded model-and-mcmc command language and built-in
  Metropolis-coupled MCMC ("MC^3") chain swapping. Single-binary
  subprocess shape (sister to BEAST 2): the CLI is `mb [-i]
  <nexus> [extras...]`. The binary is literally named `mb` (the
  project's own convention). The user authors a NEXUS file with a
  DATA block plus a MRBAYES block embedding the model / MCMC
  parameters and `mcmc` command and references it from
  `[bio.mrbayes].nexus`. `[bio.mrbayes]` knobs: `nexus` (NEXUS
  data file with embedded MRBAYES block; required), `batch`
  (default `false`; toggles `-i` so MrBayes runs the embedded
  commands non-interactively and exits cleanly rather than
  waiting on stdin at the prompt — the right default for non-
  interactive automation), `extra_args`. `prepare()` resolves the
  NEXUS path against the case directory when relative, validates
  it exists on disk, and composes the invocation with the NEXUS
  positional last so MrBayes treats it as the model file rather
  than the value of an earlier flag. `run()` streams MrBayes's
  `MrBayes v` / `Initializing` startup banner / periodic
  `Generation NNNN` / `Avg standard deviation of split
  frequencies` chain-status lines / `Analysis completed` /
  `Continue with analysis` end-of-run sentinels into progress
  hints. `collect()` walks the workdir for `*.t` (`Native`,
  "MrBayes tree samples"), `*.p` (`Tabular`, "MrBayes parameter
  samples"), and `*.con.tre` (`Native`, "MrBayes consensus
  tree"). Probe via `find_on_path(&["mb"])`. Version range
  `3.2.0..4.0.0` (the long-running stable 3.2.x line that every
  distro ships covers every release through 3.2.7; upper bound
  4.0 reserves room for an eventual major bump).
  `bio.mrbayes.mcmc` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both
adapters consume user-supplied inputs (BEAST 2 XML model files,
MrBayes NEXUS data files with embedded MRBAYES blocks) and emit
user-readable artifacts (BEAST 2 `.log` trace files + `.trees`
sampled-tree posteriors, MrBayes `.t` tree samples + `.p`
parameter samples + `.con.tre` consensus trees) that the
unchanged `Results.artifacts` collection model surfaces directly.
A first-class Bayesian-phylogenetics canonical type — a typed
posterior representation spanning both back-ends, with parsed
per-generation parameter traces and per-sample tree topologies
plus convergence-diagnostic helpers (effective sample size,
Gelman-Rubin) — defers to a future phase along with trace
visualizers, tree-density plots, and consensus-tree viewers.

Two new `valenx-init` templates ship: `beast2` with alias `beast`
(`beast2-mcmc`) and `mrbayes` with alias `mb` (`mrbayes-mcmc`).
Cross-binary roundtrip test sweeps all 90 templates clean.

Adapter inventory: 94 of 95 fully live (only `occt` remains stub-
only).

What's not in this phase: RevBayes (Sebastian Höhna's flexible
probabilistic graphical-model phylogenetics framework with a Rev
language; sister to BEAST 2 / MrBayes with a different scripting
surface; defer to Phase 30.6), BEAST 1.x (the long-lived
predecessor / sibling line to BEAST 2; defer), PhyloBayes (Nicolas
Lartillot's site-heterogeneous CAT model Bayesian phylogenetics;
defer), MIGRATE-N (Peter Beerli's Bayesian population-genetics
inference of migration rates and effective population sizes; sits
adjacent to Phase 29 population genetics; defer).

The full plan lives at
`docs/superpowers/plans/2026-04-30-bayesian-phylogenetics.md`.

### Phase 38 — Rosetta family

Open the **first Rosetta protein-modeling family** in Valenx with
the two most-used entry points into the RosettaCommons code base —
`rosetta_scripts` (the XML-driven protocol runner that's the
de-facto Rosetta entry point in production: every `relax` / `dock`
/ `abinitio` / FastDesign / enzyme-design pipeline lives as an XML
protocol fed to this binary) and PyRosetta (Python bindings
exposing the same C++ core through a Pythonic API for users who
prefer scripting Rosetta from `.py` rather than authoring XML
protocols): Rosetta (RosettaCommons' flagship modeling suite —
drives protein design, structure prediction, docking, ligand
binding, and a long tail of related modeling tasks through
`rosetta_scripts`, which reads an XML protocol describing the
modeling pipeline and applies it to an input `.pdb`; custom
Rosetta-License — academic / non-commercial use only; single-binary
subprocess shape sister to Phase 18 BWA), and PyRosetta (Python
bindings to the Rosetta C++ core — exposes the entire Rosetta
modeling pipeline (movers, filters, scorefunctions, task-
operations) through a Pythonic API; inherits Rosetta's custom
non-OSS license; Python-script subprocess shape sister to Phase 17
Biopython). Both adapters surface the RosettaCommons license
accurately via `tool_license = "Rosetta-License"` (a custom non-OSS
license — not a recognised SPDX identifier) and emit a probe
warning whenever the binary / bindings are detected, with the
literal `"academic"` string baked into the warning as a stable
anchor for license-aware filters and tests.

**This is the first Rosetta protein-modeling family to land in
Valenx.** The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 29 / 30 / 31 / 32 / 34 / 35 / 36 to cover
sequence prediction, alignment, RNA-seq, variant calling, single-
cell, transcript quantification, workflow orchestration, molecular
viewers, cheminformatics, quantum chemistry, protein design,
EvolutionaryScale models, RNA structure, population genetics,
phylogenetics, sequencing read simulation, systems biology, small-
molecule docking, CRISPR design, and cryo-EM reconstruction — but
until Phase 38 the canonical Rosetta surface (XML-protocol-driven
modeling, Python-bindings access to the core) was absent. Phase 38
closes that gap.

Two new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-rosetta` — RosettaCommons' flagship modeling
  suite (custom Rosetta-License — academic / non-commercial use
  only without a separate commercial agreement). Drives protein
  design, structure prediction, docking, ligand binding, and a
  long tail of related modeling tasks through `rosetta_scripts`,
  which reads an XML protocol describing the modeling pipeline
  (filters, movers, scorefunctions) and applies it to an input
  `.pdb`. Single-binary subprocess shape (sister to Phase 18 BWA).
  `[bio.rosetta]` knobs: `protocol` (XML protocol script;
  required), `input_pdb` (input PDB Rosetta will operate on;
  required), `output_basename` (filename stem the binary uses to
  label output decoys — `<basename>_0001.pdb` etc.; required,
  non-empty), `nstruct` (number of independent decoys to generate;
  `u32`, ≥ 1), `database` (path to the Rosetta `database/`
  directory — required because every `rosetta_scripts` invocation
  needs `-database <path>` pointing at the energy tables /
  fragment libraries / etc. bundled with the source distribution),
  `extra_args`. `prepare()` resolves the protocol, input PDB, and
  database paths against the case directory when relative,
  validates the protocol + PDB exist on disk (returns
  `InvalidCase` with a helpful message when missing), and composes
  `rosetta_scripts -database <path> -parser:protocol <protocol>
  -in:file:s <input_pdb> -out:prefix <output_basename> -nstruct <N>
  [extras...]`. `run()` streams Rosetta's `protocols.jd2` startup
  banner / `apply` per-mover lines / `Finished` / `successfully
  completed` end-of-run sentinels into progress hints. `collect()`
  walks the workdir for `<output_basename>*.pdb` (`Native`,
  "Rosetta designed structure") plus the canonical `score.sc`
  scorefile (`Tabular`, "Rosetta scores"). Probe via
  `find_on_path(&["rosetta_scripts",
  "rosetta_scripts.linuxgccrelease",
  "rosetta_scripts.macosclangrelease"])` — Rosetta source builds
  emit platform-suffixed names by default, conda / packaged
  distributions install a bare `rosetta_scripts` shim, and the
  probe covers all three. **Academic-license-only** — probe always
  pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever the binary is detected, and
  `tool_license` surfaces as `"Rosetta-License"` rather than
  mislabeling the custom RosettaCommons terms as a recognised SPDX
  identifier. Version range `3.13.0..4.0.0` (the stable 3.x line
  landed at 3.13 in 2021; upper bound 4.0 reserves room for an
  eventual major bump). `bio.rosetta.protocol` ribbon capability.
- `valenx-adapter-pyrosetta` — Python bindings to the Rosetta C++
  core (Rosetta-License — inherits the same academic / non-
  commercial use terms as the upstream Rosetta distribution).
  Exposes the entire Rosetta modeling pipeline (movers, filters,
  scorefunctions, task-operations) through a Pythonic API, letting
  users drive Rosetta from regular `.py` scripts rather than
  authoring XML protocols. Python-script subprocess shape (sister
  to Phase 17 Biopython). `[bio.pyrosetta]` knobs: `script` (path
  to user-authored Python script; required), `python` (interpreter
  name; default `"python3"`), `input_pdb` (optional input PDB the
  script will operate on — None when the script generates
  structures de novo; surfaced in `valenx_params.json` so the
  script can read it without re-parsing case.toml),
  `output_basename` (filename stem; required, non-empty).
  `prepare()` stages the script (and PDB, when present) into the
  workdir under their original filenames so the script can resolve
  them via relative paths, then writes a flat
  `valenx_params.json` with `input_pdb` (staged filename or
  literal `null` so user scripts can always do
  `params["input_pdb"]` without an `in` check) and
  `output_basename`. `run()` invokes `python <script>` via the
  shared subprocess runner; the script can emit a sentinel
  `[valenx] pyrosetta done` line on stdout to signal completion
  before exit (lifted to a 95% progress tick). `collect()` walks
  the workdir for `<output_basename>*.pdb` (`Native`, "PyRosetta
  designed structure") and `*.sc` files (`Tabular`, "PyRosetta
  scores"). Probe via Python on PATH with an `import pyrosetta`
  check — when the import fails the probe still returns
  `ok = true` with a warning so users with PyRosetta installed
  under a different interpreter (referenced via the case-level
  `python` override) aren't blocked. **Academic-license-only** —
  probe always pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever Python is detected (regardless
  of whether `pyrosetta` itself is importable, since the user is
  either about to install it or has it installed and needs
  reminding), and `tool_license` surfaces as `"Rosetta-License"`
  rather than mislabeling the inherited RosettaCommons terms as
  MIT / BSD. Version range `4.0.0..5.0.0` (the modern release line
  is the 4.x series with weekly nightly drops post-2017; upper
  bound 5.0 reserves room for an eventual major bump).
  `bio.pyrosetta.script` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** Both
adapters consume user-supplied inputs (Rosetta XML protocols +
input PDBs + `database/` data directories, PyRosetta Python
scripts + optional input PDBs) and emit user-readable artifacts
(`<basename>*.pdb` design decoys, `score.sc` / `*.sc` scorefiles)
that the unchanged `Results.artifacts` collection model surfaces
directly. The existing `valenx_bio::format::pdb` reader inspects
collected PDB artifacts for chain / residue / atom counts. A
first-class Rosetta canonical type — a generic protocol +
scorefile pair spanning both back-ends, parsed into a typed
scorefile model with per-decoy energy terms — defers to a future
phase along with score-distribution visualizers and per-mutation
Δ-energy heatmap viewers.

Two new `valenx-init` templates ship: `rosetta`
(`rosetta-protocol`) and `pyrosetta` (`pyrosetta-script`).
Cross-binary roundtrip test sweeps all 85 templates clean.

Adapter inventory: 89 of 90 fully live (only `occt` remains stub-
only).

What's not in this phase: specific Rosetta app adapters that wrap
individual Rosetta apps (`relax`, `dock_protocol`, `AbinitioRelax`,
`enzyme_design`, `loopmodel`) directly rather than going through
`rosetta_scripts` XML protocols (defer to Phase 38.5 — the
`rosetta_scripts` protocol runner covers all of them via XML
protocols). Rosetta@home / RosettaCommons cluster execution
(different shape — distributed work-unit dispatch rather than
single-host modeling; defer to a future phase).

The full plan lives at
`docs/superpowers/plans/2026-04-30-rosetta.md`.

### Phase 29 — Population genetics

Open the **first population-genetics / evolutionary-simulation
domain** in Valenx with three established open-source tools that
span the population-genetics tradeoff space — forward-time
individual-based simulation under arbitrary selection / demography
/ mating-system specifications (SLiM), coalescent backward-time
simulation of sample ancestries under configurable demographies
(msprime), and tree-sequence analysis / statistics on the succinct
tree-sequence outputs both simulators emit (tskit): SLiM (Philipp
Messer's forward-time population-genetics simulator — evolves a
finite-population model generation by generation under a user-
defined Eidos script (mutation rates, selection coefficients,
recombination maps, demographic events, migrations, mating
systems); tree-sequence recording (`treeSeqOutput()` family) feeds
straight into tskit / msprime downstream; GPL-3.0; single-binary
subprocess shape sister to Phase 18 BWA), msprime (Jerome
Kelleher's coalescent backwards-in-time population-genetics
simulator — speed-of-light coalescent simulator (millions of
samples per minute on a workstation); the canonical companion to
SLiM and tskit; GPL-3.0; Python-script subprocess shape sister to
Phase 17 Biopython), and tskit (the canonical tree-sequence
analysis library, MIT — built around the succinct tree-sequence
data structure pioneered by msprime; computes population-genetics
statistics (π, Tajima's D, Fst, site-frequency spectra, IBD
shares); the workhorse downstream of every Phase 29 simulator —
msprime emits `.trees`, SLiM emits `.trees`, tskit consumes them).
SLiM follows the established Phase 18 BWA single-binary CLI
pattern (script positional last). msprime + tskit follow the
established Phase 17 Biopython Python-script subprocess pattern
(user authors a `.py` driver; the adapter stages script + writes
`valenx_params.json`; run() invokes `python <script>`).

**This is the first population-genetics / evolutionary-simulation
domain to land in Valenx.** The biology adapter family started
with Phase 17 (foundation — sequence / structure / trajectory
canonical types + classical MD + cheminformatics scripts) and
expanded through Phase 17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 /
22 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 28 / 30 / 31 / 32 / 34 / 35
/ 36 / 38 to cover sequence prediction, alignment, RNA-seq,
variant calling, single-cell, transcript quantification, workflow
orchestration, molecular viewers, cheminformatics, quantum
chemistry, protein design, EvolutionaryScale models, RNA
structure, phylogenetics, sequencing read simulation, systems
biology, small-molecule docking, CRISPR design, cryo-EM
reconstruction, and Rosetta protein modeling — but until Phase 29
the population-genetics surface (forward-time individual-based
evolutionary simulation, coalescent backward-time ancestry
simulation, tree-sequence analysis) was absent. Phase 29 closes
that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-slim` — Philipp Messer's forward-time population-
  genetics simulator (GPL-3.0). Evolves a finite-population model
  generation by generation under a user-defined Eidos script:
  mutation rates, selection coefficients, recombination maps,
  demographic events, migrations, mating systems. The state is
  sampled at any generation the script asks for, and tree-sequence
  recording (the `treeSeqOutput()` family) feeds straight into
  tskit / msprime downstream. Single-binary subprocess shape
  (sister to Phase 18 BWA): the CLI is `slim [-s <seed>]
  [extras...] <script>`. `[bio.slim]` knobs: `script` (`.slim`
  Eidos model file; required), `seed` (optional `u64`; passed via
  `slim -s <N>` when present, otherwise SLiM picks its own seed
  and prints it on the run banner), `output_basename` (filename
  stem the user's script uses for outputs — surfaced here so
  collect() can label artefacts uniformly even though SLiM scripts
  choose their own output paths; required, non-empty),
  `extra_args` (additional CLI arguments appended after the
  script path; `-d KEY=VALUE` is the canonical way to inject Eidos
  constants from outside the script). `prepare()` resolves the
  script against the case directory when relative, validates it
  exists on disk, and composes the invocation with `seed` injected
  before any extras and the script positional last so SLiM treats
  it as the model file rather than the value of an earlier flag.
  `run()` streams SLiM's `// Initial random seed` banner /
  periodic `// generation N` lines / `// Run finished` end-of-run
  sentinel into progress hints. `collect()` walks the workdir for
  any `<output_basename>*.trees` (`Native`, "SLiM tree sequence")
  and `<output_basename>*.log` (`Log`) the script emitted — the
  adapter doesn't try to predict the exact filenames the script
  will write, since SLiM scripts choose their own output paths via
  `writeFile()` / `treeSeqOutput()` calls. Probe via
  `find_on_path(&["slim"])` (conda-forge / source / Homebrew all
  install under the canonical lowercase `slim` name); the generic
  version detector tries both the conventional `--version` and the
  SLiM-native `-version` form. Version range `4.0.0..5.0.0` (the
  modern release line is the 4.x series from 2022+; 4.0
  introduced the streamlined Eidos type system and the
  `treeSeqOutput()` helpers we rely on for tskit interop; upper
  bound 5.0 reserves room for an eventual major bump).
  `bio.slim.simulate` ribbon capability.
- `valenx-adapter-msprime` — Jerome Kelleher's coalescent
  backwards-in-time population-genetics simulator (GPL-3.0).
  Simulates the ancestry of a sample under a configurable
  demography and recombination map, then layers mutations onto the
  resulting tree sequence. It is the speed-of-light coalescent
  simulator (millions of samples per minute on a workstation) and
  the canonical companion to SLiM (forward-time) and tskit (tree-
  sequence analysis). Python-script subprocess shape (sister to
  Phase 17 Biopython): the user authors a `simulate.py` referenced
  from `[bio.msprime].script` in `case.toml`. `[bio.msprime]`
  knobs: `script` (path to user-authored Python script; required),
  `python` (interpreter name; default `"python3"`),
  `population_size` (`u32`, ≥ 1), `num_samples` (`u32`, ≥ 1),
  `recombination_rate` (`f64`, ≥ 0.0 and finite — per-site per-
  generation rate), `mutation_rate` (`f64`, ≥ 0.0 and finite —
  per-site per-generation rate), `output_basename` (filename stem;
  required, non-empty). `prepare()` stages the script into the
  workdir under its original filename and writes a flat
  `valenx_params.json` containing `population_size`,
  `num_samples`, `recombination_rate` (emitted via `{:e}` so
  Python's `json.load` parses it back as a float), `mutation_rate`
  (same), and `output_basename`. `run()` invokes `python <script>`
  via the shared subprocess runner. `collect()` walks the workdir
  for `<output_basename>.trees` (`Native`, "msprime tree
  sequence"), `<output_basename>.vcf` (`Tabular`, "msprime VCF"),
  and `<output_basename>.csv` (`Tabular`, "msprime per-sample
  summary") — user scripts emit any combination of these via
  msprime / tskit's tabular APIs. Probe via Python on PATH with an
  `import msprime` check — when the import fails the probe still
  returns `ok = true` with a warning so users with msprime
  installed under a different interpreter (referenced via the
  case-level `python` override) aren't blocked. Version range
  `1.3.0..2.0.0` (the modern `sim_ancestry()` /
  `sim_mutations()` split landed in 1.3 in 2024, paired with the
  tskit 0.5+ tree-sequence format we surface in collect(); upper
  bound 2.0 reserves room for an eventual major bump).
  `bio.msprime.simulate` ribbon capability.
- `valenx-adapter-tskit` — the canonical tree-sequence analysis
  library (MIT), built around the succinct tree-sequence data
  structure pioneered by msprime. Computes population-genetics
  statistics (π, Tajima's D, Fst, site-frequency spectra, IBD
  shares), exposes per-tree iteration across the genome, converts
  between tree-sequence and VCF / Newick / table formats, and
  renders phylogenetic plots. It's the workhorse downstream of
  every Phase 29 simulator — msprime emits `.trees`, SLiM emits
  `.trees`, tskit consumes them. Python-script subprocess shape
  (sister to msprime). `[bio.tskit]` knobs: `script` (path to
  user-authored Python script; required), `python` (interpreter
  name; default `"python3"`), `input_trees` (`.trees` file from
  SLiM or msprime; required), `output_basename` (filename stem;
  required, non-empty). `prepare()` stages script + tree-sequence
  file into the workdir under their original filenames so the
  script can resolve them via relative paths, then writes a flat
  `valenx_params.json` containing `input_trees` (staged filename)
  and `output_basename`. `run()` invokes `python <script>` via the
  shared subprocess runner. `collect()` walks the workdir for
  `<output_basename>*.csv` / `<output_basename>*.tsv` (`Tabular`,
  "tskit statistics") and `*.png` (`Native`, "tskit plot") — user
  scripts emit any combination of statistics tables and rendered
  plots. Probe via Python on PATH with an `import tskit` check —
  same `ok = true` + warning fallback as msprime. Version range
  `0.5.0..1.0.0` (tskit 0.5+ ships the modern `Statistics` API
  surface and the v3 tree-sequence file format msprime 1.3+
  writes; upper bound 1.0 reserves room for the long-promised
  1.0 release). `bio.tskit.analyze` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (SLiM `.slim` Eidos scripts;
msprime + tskit Python scripts; tskit input `.trees` tree-sequence
files emitted by SLiM or msprime) and emit user-readable artifacts
(`.trees` tree sequences, `.vcf` genotype calls, `.csv` / `.tsv`
statistics tables, `.png` rendered plots, `.log` run logs) that
the unchanged `Results.artifacts` collection model surfaces
directly. A first-class population-genetics canonical type — a
typed tree-sequence representation spanning all three back-ends,
with parsed per-tree edge / node / mutation tables and a typed
statistics-table representation — defers to a future phase along
with tree-sequence visualizers and per-population allele-frequency-
spectrum viewers.

Three new `valenx-init` templates ship: `slim` (`slim-simulate`),
`msprime` (`msprime-simulate`), and `tskit` (`tskit-analyze`).
Cross-binary roundtrip test sweeps all 85 templates clean.

Adapter inventory: 89 of 90 fully live (only `occt` remains stub-
only).

What's not in this phase: fwdpy11 (Kevin Thornton's Python-driven
forward-time population-genetics simulator; sister to SLiM with a
different scripting surface; defer to Phase 29.5), simuPOP (Bo
Peng's Python forward simulator with a long history; defer),
pyslim (the SLiM tree-sequence Python interop layer that bridges
SLiM `.trees` outputs into msprime / tskit; defer to 29.5),
stdpopsim (the standardised population-genetics simulation library
that wraps msprime / SLiM under a catalog of curated demographic
models; defer), demes (the human-readable demographic-model
specification format; defer to 29.5).

The full plan lives at
`docs/superpowers/plans/2026-04-30-population-genetics.md`.

### Phase 36 — Cryo-EM

Open the **first cryo-electron microscopy reconstruction domain** in
Valenx with three established open-source tools that span the cryo-EM
pipeline — Bayesian 3D reconstruction at the core (RELION), broad-
spectrum image processing across the full single-particle workflow
(EMAN2), and per-micrograph contrast transfer function (CTF)
estimation as the canonical preprocessing step (CTFFIND): RELION
(Sjors Scheres' REgularised LIkelihood OptimisatioN suite — de-facto
Bayesian 3D reconstruction workhorse in cryo-EM facilities worldwide,
GPL-2.0; single-binary `relion_refine` for the single-process path
or `mpirun -n <N> relion_refine_mpi` for multi-rank, since RELION
ships separate `_mpi`-suffixed binaries), EMAN2 (Steve Ludtke's
broad-spectrum cryo-EM image-processing package — "Swiss army knife"
of single-particle cryo-EM, BSD-3-Clause; high-level driver
`e2refine_easy.py` orchestrates particle picking, 2D classification,
initial-model building, and 3D refinement across the sprawling
`e2*.py` toolkit), and CTFFIND (Niko Grigorieff's CTF estimation
tool — gold standard for fitting per-micrograph CTF parameters that
RELION, cryoSPARC, EMAN2, and most automated pipelines all wrap as
a preprocessing step; **academic-license-flagged** under the Janelia
Research Campus non-commercial license, surfaced as `Janelia-License`
and flagged via mandatory `"academic"`-keyworded probe warning;
single-binary `ctffind` with stdin-piped parameters since the CLI is
interactive). All three follow the established Phase 18 BWA single-
binary CLI pattern: particles / micrographs / reference maps in,
typed artifacts out.

**This is the first cryo-electron microscopy reconstruction domain
to land in Valenx.** The biology adapter family started with Phase
17 (foundation — sequence / structure / trajectory canonical types
+ classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 30 / 32 / 34 to cover sequence prediction,
alignment, RNA-seq, variant calling, single-cell, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design, EvolutionaryScale
models, RNA structure, phylogenetics, systems biology, and small-
molecule docking — but until Phase 36 the cryo-electron microscopy
reconstruction surface (Bayesian 3D refinement, broad-spectrum image
processing, CTF estimation) was absent. Phase 36 closes that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-relion` — Sjors Scheres' REgularised LIkelihood
  OptimisatioN suite (GPL-2.0). The de-facto Bayesian 3D
  reconstruction workhorse in cryo-EM facilities worldwide —
  particle classification, 3D refinement, CTF correction, post-
  processing. Single-binary subprocess shape with optional MPI
  wrapping: `relion_refine` for the single-process path,
  `mpirun -n <N> relion_refine_mpi` for multi-rank runs (RELION
  ships these as separate `_mpi`-suffixed binaries so the launcher
  knows which transport to use). `[bio.relion]` knobs: `particles`
  (`*_data.star`; required), `reference` (`.mrc`; required),
  `output_basename` (becomes the `--o` prefix every output inherits
  so collect() walks deterministically; required), `angpix` (pixel
  size in Angstroms; required, > 0.0 and finite), `mpi_procs`
  (default 1, ≥ 1; > 1 switches to the MPI binary), `threads`
  (OpenMP threads per MPI rank, default 1, ≥ 1), `extra_args`.
  `prepare()` dispatches on `mpi_procs`: single-rank composes
  `relion_refine --i <particles> --ref <reference> --o
  <output_basename> --angpix <angpix> --j <threads> [extras...]`;
  multi-rank prepends `mpirun -n <mpi_procs> relion_refine_mpi ...`
  and surfaces a helpful install-hint `InvalidCase` ("install
  OpenMPI (`apt install openmpi-bin`, `brew install open-mpi`) or
  MPICH (`apt install mpich`) to enable multi-rank RELION runs") if
  `mpirun` isn't on PATH. `collect()` walks the workdir recursively
  for `<output_basename>*_class*.mrc` (`Native`, "RELION class
  average"), `<output_basename>*_data.star` (`Tabular`, "RELION
  particle assignments"), and `<output_basename>*_model.star`
  (`Log`, "RELION model summary"). Probe via
  `find_on_path(&["relion_refine"])`. Version range `4.0.0..6.0.0`
  (4.0 is the current stable line, predecessor 3.1; upper bound
  6.0 reserves room for the next major). `bio.relion.refine`
  ribbon capability.
- `valenx-adapter-eman2` — Steve Ludtke's broad-spectrum cryo-EM
  image-processing package (BSD-3-Clause). The "Swiss army knife"
  of single-particle cryo-EM: particle picking, 2D classification,
  initial-model building, 3D refinement (CTF corrected, with
  simultaneous tilt-pair handling), and a sprawling Python toolkit
  (`e2*.py`) for everything in between. Single-binary subprocess
  shape: the adapter wraps `e2refine_easy.py`, EMAN2's high-level
  orchestrator that drives the rest of the toolkit. `[bio.eman2]`
  knobs: `particles` (particle stack `.bdb` / `.hdf` / `.mrcs`;
  required), `model` (initial 3D model `.hdf` / `.mrc`; required),
  `output_basename` (becomes the `--path` argument; EMAN2 turns
  this into a `<basename>_NN/` results directory under the workdir;
  required), `target_resolution` (target resolution in Angstroms;
  required, > 0.0 and finite), `symmetry` (point group — `"c1"` /
  `"d2"` / `"icos"` / etc.; required, default `"c1"`), `threads`
  (default 1, ≥ 1), `extra_args`. `prepare()` builds
  `e2refine_easy.py --input <particles> --model <model> --path
  <output_basename> --targetres <target_resolution> --sym <symmetry>
  --threads <threads> [extras...]`. `collect()` walks recursively
  for `<output_basename>_*/threed_*.hdf` (`Native`, "EMAN2
  reconstruction") and `<output_basename>_*/log.txt` (`Log`, "EMAN2
  log"). Probe via `find_on_path(&["e2refine_easy.py"])`. Version
  range `2.99.0..3.0.0` (the 2.99 line is the current pre-3.0
  stable release; upper bound 3.0 reserves room for the long-
  rumoured 3.x line). The `valenx-init` template ships with the
  alias `eman` alongside the canonical `eman2`.
  `bio.eman2.refine` ribbon capability.
- `valenx-adapter-ctffind` — Niko Grigorieff's contrast transfer
  function (CTF) estimation tool (Janelia Research Campus non-
  commercial / academic-only license, surfaced as `Janelia-License`
  rather than mislabeling as MIT / BSD). The gold standard for
  fitting per-micrograph CTF parameters (defocus, astigmatism, phase
  shift) in single-particle cryo-EM workflows; RELION, cryoSPARC,
  EMAN2, and most automated pipelines all wrap CTFFIND under the
  hood as a preprocessing step. Single-binary subprocess shape with
  stdin-piped parameters: CTFFIND's CLI is interactive and prompts
  the user line-by-line for each microscope parameter on startup,
  so the adapter writes a parameters text file in the workdir
  during `prepare()` and uses a custom `run()` that pipes the file
  into the child's stdin via `Stdio::from(file)`. The shared
  `subprocess::run` helper closes stdin (`Stdio::null()`); CTFFIND
  on a closed stdin reads EOF before its first prompt and exits
  with an error. The custom run path mirrors the MAFFT stdout-
  redirect pattern but for stdin. `[bio.ctffind]` knobs:
  `micrograph` (input micrograph `.mrc`; required),
  `output_diagnostic` (output diagnostic image `.mrc`; required),
  `output_txt` (output text file with CTF parameters; required),
  `pixel_size` (Angstroms; required, > 0.0 and finite), `voltage`
  (acceleration voltage in kV; default 300.0, > 0.0), `cs`
  (spherical aberration in mm; default 2.7, > 0.0),
  `amplitude_contrast` (fraction; required, in `0.0..=1.0` — 0.07
  typical for cryo, 0.1 for negative stain), `extra_args`.
  `prepare()` writes `ctffind_params.txt` containing one parameter
  per CTFFIND-v4.1 prompt in order (input image, output diagnostic,
  pixel size, voltage, Cs, amplitude contrast, plus standard
  defaults for box size / min res / max res / defocus search /
  expert sub-prompts) and stashes the filename under a sentinel env
  var (`VALENX_CTFFIND_PARAMS_FILE`). The custom `run()` recovers
  the filename, strips the sentinel from the env table so CTFFIND
  doesn't see it, opens the params file with `File::open()`, and
  hands the FD to the child — CTFFIND sees a pipe pre-loaded with
  one parameter per prompt and responds as if a human had typed
  each line. `collect()` reports `output_diagnostic` (`Native`,
  "CTFFIND diagnostic image") and `output_txt` (`Tabular`, "CTFFIND
  parameters"). Probe via `find_on_path(&["ctffind"])`. Version
  range `4.1.0..5.0.0` (CTFFIND4 is the long-running stable line;
  upper bound 5.0 reserves room for the announced CTFFIND5 line).
  `bio.ctffind.estimate` ribbon capability.
  **Academic-license-only — non-OSS Janelia Research Campus terms.**
  Probe pushes the literal string `"academic"` (the asserted
  anchor) into `ProbeReport.warnings` with the full reminder:
  "CTFFIND is licensed for non-commercial / academic use only.
  Confirm your use case complies with the Janelia license before
  redistributing CTF estimates or derived data." Tool license
  surfaces as `Janelia-License` rather than mislabeling as MIT /
  BSD. Downstream license-aware tooling and tests key off the
  literal `"academic"` anchor.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (RELION particle STAR files +
reference MRC volumes, EMAN2 particle stacks `.bdb` / `.hdf` /
`.mrcs` + initial 3D models, CTFFIND micrograph `.mrc` files) and
emit user-readable artifacts (RELION class-average MRC volumes,
particle-assignment STAR files, model-summary STAR files, EMAN2
`threed_*.hdf` reconstructions plus per-run log files, CTFFIND
diagnostic-image MRC plus per-micrograph parameter text files) that
the unchanged `Results.artifacts` collection model surfaces directly.
A first-class cryo-EM canonical type — a generic `.mrc` volume /
particle-stack / micrograph type spanning all three back-ends —
defers to a future phase along with MRC readers and reconstruction
visualizers.

Three new `valenx-init` templates ship: `relion` (`relion-refine`),
`eman2` with alias `eman` (`eman2-refine`), and `ctffind`
(`ctffind-estimate`; the scaffolded `case.toml` carries an inline
academic-license note). Cross-binary roundtrip test sweeps all 74
templates clean.

Adapter inventory: 78 of 79 fully live (only `occt` remains stub-
only).

What's not in this phase: cisTEM (full single-particle cryo-EM
pipeline UI; defer to Phase 36.5), SPHIRE (TransPhire / SPHIRE
pipeline framework; defer), IMOD (cryo-ET reconstruction; different
shape; defer), Bsoft (broad cryo-EM + electron crystallography
toolkit; defer), Scipion (full cryo-EM pipeline orchestrator akin
to Nextflow / Snakemake but cryo-EM-specific; defer), Frealign (the
cisTEM predecessor; defer), motion correction (MotionCor2 /
RELION's own `relion_motioncorr`; defer to 36.5), particle picking
(Topaz / crYOLO; defer), tomography (TomoBEAR / nextPYP; different
shape — slot under cryo-ET separately). Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-cryo-em.md`.

### Phase 35 — CRISPR design

Open the **first CRISPR guide-RNA design domain** in Valenx with
three established open-source tools that span the CRISPR-design
tradeoff space — popular ranked guide design with off-target
scoring (CHOPCHOP), comprehensive guide design plus rigorous off-
target prediction across many enzymes (CRISPOR), and pure off-
target searching used as a primitive by most other CRISPR-design
web services and pipelines (Cas-OFFinder): CHOPCHOP (University of
Bergen's web-and-script CRISPR guide-RNA design tool — de-facto
first stop in academic CRISPR workflows; scores candidate gRNAs
against a target sequence under a configurable nuclease (Cas9,
Cas12a, Cas13) or TALEN design pass and ranks by efficiency /
specificity / off-target risk; MIT; Python-script subprocess shape
sister to Phase 17 Biopython), CRISPOR (Maximilian Haeussler's
CRISPR guide-RNA design + off-target prediction tool behind the
public crispor.org service — distinguishing feature is the rigorous
off-target pass via the CFD scoring model and MIT-style specificity
scores per guide; supports many more enzymes / PAMs than CHOPCHOP;
GPL-3.0; Python-script subprocess shape sister to CHOPCHOP), and
Cas-OFFinder (Bae / Park / Kim group's CRISPR off-target searching
tool from Hanyang / Seoul National University — fast, OpenCL-
accelerated scanner that walks a reference genome and reports every
position whose sequence matches one of the input guides within the
configured Hamming distance; the workhorse off-target scanner
sitting under most CRISPR design web services and pipelines;
BSD-3-Clause; single-binary subprocess shape sister to Phase 18 BWA
with fixed-shape positional CLI `cas-offinder <input> {C|G|A}
<output> [extras...]` — no `-i` / `-o` flags, the order is fixed).
CHOPCHOP + CRISPOR follow the established Phase 17 Biopython Python-
script subprocess pattern: the user supplies a Python script
referenced from `[bio.<adapter>].script` that imports the upstream
package and reads a flat `valenx_params.json` for the parsed knobs.
Cas-OFFinder follows the established Phase 18 BWA single-binary CLI
pattern.

**This is the first CRISPR guide-RNA design domain to land in
Valenx.** The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 30 / 31 / 32 / 34 / 36 to cover sequence
prediction, alignment, RNA-seq, variant calling, single-cell,
transcript quantification, workflow orchestration, molecular
viewers, cheminformatics, quantum chemistry, protein design,
EvolutionaryScale models, RNA structure, phylogenetics, sequencing
read simulation, systems biology, small-molecule docking, and cryo-
EM reconstruction — but until Phase 35 the CRISPR-design surface
(guide-RNA scoring, off-target prediction, off-target searching) was
absent. Phase 35 closes that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-chopchop` — University of Bergen's web-and-script
  CRISPR guide-RNA design tool (MIT). De-facto first stop for "I
  have a gene, what should I cut" in academic CRISPR workflows:
  scores candidate gRNAs against a target sequence under a
  configurable nuclease (Cas9, Cas12a, Cas13) or TALEN design pass,
  ranks by efficiency / specificity / off-target risk, and emits
  both a guide-ranking TSV and a guide-location BED. Python-script
  subprocess shape (sister to Phase 17 Biopython): the user supplies
  a Python script referenced from `[bio.chopchop].script` in
  `case.toml` that imports `chopchop` (or invokes `chopchop.py`)
  and reads `valenx_params.json` for the parsed knobs. `[bio.chopchop]`
  knobs: `script` (path to user-supplied Python script; required),
  `python` (interpreter name; default `"python3"`), `target` (target
  sequence FASTA; required), `genome` (CHOPCHOP-installed genome
  name — `"hg38"` / `"mm10"` / etc.; required), `cas_variant` (one
  of `"Cas9"` / `"Cas12a"` / `"Cas13"` / `"TALEN"`; required), `pam`
  (PAM sequence — `"NGG"` for Cas9, `"TTTV"` for Cas12a, etc.;
  required), `output_basename` (filename stem; required, non-empty).
  `prepare()` stages the script + target FASTA into the workdir,
  writes a flat `valenx_params.json` containing `target` (staged
  filename), `genome`, `cas_variant`, `pam`, `output_basename`, and
  composes `python <script_filename>` as the native command.
  `collect()` walks the workdir for `<output_basename>*.tsv`
  (`Tabular`, "CHOPCHOP guide rankings") and `<output_basename>*.bed`
  (`Tabular`, "CHOPCHOP guide locations"). Probe via Python on PATH
  with an `import chopchop` check — when the import fails the probe
  still returns `ok = true` with a warning so users with CHOPCHOP
  installed under a non-standard module name (`crispr_chopchop`,
  `chopchop_v3`) or invoked via the user-supplied script aren't
  blocked. Version range `3.0.0..4.0.0` (the modern web / script
  split landed in 3.0; upper bound 4.0 reserves room for an
  eventual major bump). `bio.chopchop.design` ribbon capability.
- `valenx-adapter-crispor` — Maximilian Haeussler's CRISPR guide-RNA
  design + off-target prediction tool (GPL-3.0). CRISPOR's
  distinguishing feature is the rigorous off-target pass: scores
  candidate guides against a reference genome assembly with the CFD
  scoring model and reports an MIT-style specificity score per
  guide. It powers the public crispor.org service and is also
  distributed as a standalone Python script for batch / pipeline
  use, supporting many more enzymes / PAMs than CHOPCHOP. Python-
  script subprocess shape (sister to CHOPCHOP). `[bio.crispor]`
  knobs: `script` (path to user-supplied Python script; required),
  `python` (interpreter name; default `"python3"`), `target` (target
  sequence FASTA; required), `genome` (CRISPOR-supported genome
  name; required), `pam` (PAM motif — `"NGG"` / `"NG"` / `"TTTV"` /
  etc.; required), `batch_id` (optional string; CRISPOR caches
  partial results by batch so passing the same `batch_id` resumes a
  previously-interrupted run), `output_basename` (filename stem;
  required, non-empty). `prepare()` stages the script + target
  FASTA, writes a flat `valenx_params.json` containing `target`
  (staged filename), `genome`, `pam`, `batch_id` (JSON string or
  literal `null` so user scripts can always do `params["batch_id"]`
  without an `in` check), `output_basename`, and composes `python
  <script_filename>` as the native command. `collect()` walks the
  workdir for `<output_basename>*.tsv` (`Tabular`, "CRISPOR guide
  rankings") and `<output_basename>*.txt` (`Log`). Probe via Python
  on PATH with an `import crispor` check (same `ok = true` +
  warning fallback as CHOPCHOP). Version range `5.0.0..6.0.0` (the
  modern Python 3 / batch-mode rewrite landed in 5.0; upper bound
  6.0 reserves room for an eventual major bump). `bio.crispor.design`
  ribbon capability.
- `valenx-adapter-cas-offinder` — Bae / Park / Kim group's CRISPR
  off-target searching tool from Hanyang / Seoul National University
  (BSD-3-Clause). Cas-OFFinder is a fast, OpenCL-accelerated
  scanner: given a list of guide sequences + PAM patterns + mismatch
  budget in a plain-text input file, it walks a reference genome
  and reports every position whose sequence matches one of the
  guides within the configured Hamming distance. It's the workhorse
  off-target scanner sitting under most CRISPR design web services
  (CRISPOR, CRISPRdirect, …) and pipelines. Single-binary subprocess
  shape (sister to Phase 18 BWA): the CLI is fixed-shape `cas-offinder
  <input> {C|G|A} <output> [extras...]`. `<input>` is a 3+-line text
  file with the reference genome path, the PAM pattern, and one
  guide-sequence row per query. The middle positional argument
  selects the OpenCL device class — `C` (CPU), `G` (GPU), or `A`
  (auto-pick fastest at runtime). `[bio.cas-offinder]` knobs: `input`
  (Cas-OFFinder input file; required), `output` (output text file;
  required), `backend` (one of `"C"` / `"G"` / `"A"`; required),
  `extra_args`. `prepare()` resolves both paths against the case
  directory (when relative) and composes the invocation positionally
  — no `-i` / `-o` flags, the order is fixed. `collect()` reports
  the configured `output` file as a single `Tabular` artifact
  ("Cas-OFFinder off-target hits"). Probe via
  `find_on_path(&["cas-offinder"])`. Version range `2.4.0..3.0.0`
  (the modern OpenCL device-selection CLI stabilised at 2.4; upper
  bound 3.0 reserves room for an eventual major bump). The init
  alias `cas-off` resolves to the same template as the canonical
  `cas-offinder`. `bio.cas-offinder.search` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied inputs (CHOPCHOP + CRISPOR target
FASTAs and design Python scripts, Cas-OFFinder fixed-shape input
files specifying the genome path + PAM + per-guide query lines)
and emit user-readable artifacts (CHOPCHOP guide-ranking TSV +
guide-location BED, CRISPOR guide-ranking TSV + log TXT,
Cas-OFFinder ranked off-target hit TSV) that the unchanged
`Results.artifacts` collection model surfaces directly. The existing
`valenx_bio::format::fasta` reader already inspects target FASTAs
for sequence count + identifiers + alphabets. A first-class CRISPR-
design canonical type — a generic guide / off-target / scoring type
spanning all three back-ends — defers to a future phase along with
guide-ranking visualizers and off-target heatmap viewers.

Three new `valenx-init` templates ship: `chopchop` (`chopchop-design`),
`crispor` (`crispor-design`), and `cas-offinder` with alias `cas-off`
(`cas-offinder-search`). Cross-binary roundtrip test sweeps all
80 templates clean.

Adapter inventory: 84 of 85 fully live (only `occt` remains stub-
only).

What's not in this phase: CRISPRitz (Pinello lab in-silico off-target
search with variant-aware scoring; sister to Cas-OFFinder; defer to
Phase 35.5), FlashFry (Aaron McKenna's high-throughput guide design
+ scoring; sister to CHOPCHOP; defer), E-CRISP (Boutros lab guide
design with conservation scoring; defer), CRISPRdirect (Naito lab
web-service guide selector; defer), Guidescan (off-target
enumeration via specificity scoring; defer to 35.5), CRISPResso2
(post-editing analysis from sequencing data — different shape:
mapped read alignment + indel calling rather than guide-RNA design;
defer to a future phase).

The full plan lives at
`docs/superpowers/plans/2026-04-30-crispr-design.md`.

### Phase 31 — Sequencing read simulators

Open the **first sequencing read-simulation domain** in Valenx with
three established open-source tools that span all three major
sequencing-technology classes — per-platform empirical-error-profile
Illumina short reads (ART), the simple-uniform-error classic short-
read baseline that ships with samtools (wgsim), and realistic-error-
profile Nanopore long reads (Badread): ART (Weichun Huang's NIEHS
Illumina-platform read simulator — de-facto choice for synthesising
FASTQs that match per-platform empirical error profiles for HiSeq
2500 / HiSeq X / MiSeq v3 / NextSeq 500 / MiniSeq, GPL-3.0; single-
binary `art_illumina` with paired-end dispatch via `-p -m <mean> -s
<sd>`), wgsim (Heng Li's classic Whole-Genome SIMulator that ships
alongside samtools, MIT — always paired-end, always position-uniform,
deliberately simple under a uniform sequencing-error model; the
canonical "small + classic" simulator for fast smoke-testing of
mappers and variant callers when realistic error spectra are not
required; single-binary `wgsim` with positional output arguments
after the reference, no stdout-redirect needed), and Badread (Ryan
Wick's long-read simulator with realistic Nanopore + PacBio CLR
error profiles, GPL-3.0 — random / chimeric / adapter / glitch read
types, junk-read injection, identity drift, length distributions
calibrated against actual sequencer output; single-binary `badread
simulate` that writes its FASTQ to stdout via the MAFFT-style
stdout-redirect-to-file pattern). All three follow the established
Phase 18 BWA single-binary CLI pattern: reference FASTA in,
simulated FASTQ(s) out.

**This is the first sequencing read-simulation domain to land in
Valenx.** The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 30 / 32 / 34 / 36 to cover sequence prediction,
alignment, RNA-seq, variant calling, single-cell, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design,
EvolutionaryScale models, RNA structure, phylogenetics, systems
biology, small-molecule docking, and cryo-EM reconstruction — but
until Phase 31 the read-simulation surface (synthetic FASTQ
generation across the three major sequencing-technology classes)
was absent. Phase 31 closes that gap.

Three new adapter crates land under `crates/valenx-adapters/bio/`:

- `valenx-adapter-art` — Weichun Huang's NIEHS Illumina-platform
  read simulator (GPL-3.0). The de-facto choice for synthesising
  FASTQs that match the empirical error profile of a given Illumina
  sequencing system (HiSeq 2500, HiSeq X, MiSeq v3, NextSeq 500,
  MiniSeq) so downstream pipelines can be validated against a
  known-truth reference at controlled coverage and read length.
  Single-binary subprocess shape: the adapter wraps `art_illumina`,
  the workhorse of the ART family (companion `art_454` / `art_SOLiD`
  binaries cover platforms this adapter does not surface).
  `[bio.art]` knobs: `reference` (FASTA; required), `output_prefix`
  (filename stem; ART writes `<prefix>.fq` for single-end or
  `<prefix>1.fq` + `<prefix>2.fq` for paired-end; required, non-
  empty), `sequencing_system` (one of `"HS25"` / `"HSXt"` / `"MSv3"`
  / `"NS50"` / `"MinS"`; required), `read_length` (≥ 1; required),
  `fold_coverage` (> 0.0; required), `paired_end` (default `false`),
  `fragment_mean` (mean insert size for paired-end; default 200.0,
  > 0.0 when `paired_end`), `fragment_sd` (insert-size stddev for
  paired-end; default 10.0, > 0.0 when `paired_end`), `extra_args`.
  `prepare()` builds `art_illumina -ss <sequencing_system> -i
  <reference> -l <read_length> -f <fold_coverage> -o <output_prefix>
  [-p -m <fragment_mean> -s <fragment_sd> if paired_end]
  [extras...]`. `collect()` walks the workdir top-level for
  `<output_prefix>*.fq` (`Tabular`, "ART simulated reads") and
  `<output_prefix>*.aln` (`Log`, "ART alignment record" — the per-
  read alignment record ART writes alongside the FASTQ, useful for
  validating aligner accuracy against the simulated truth). Probe
  via `find_on_path(&["art_illumina"])`. Version range
  `2.5.0..3.0.0` (the long-running ChocolateCherryCake `2.5.x`
  series since 2016; Bioconda + Homebrew ship 2.5.8). The
  `valenx-init` template ships with the alias `art-illumina`
  alongside the canonical `art`. `bio.art.simulate` ribbon
  capability.
- `valenx-adapter-wgsim` — Heng Li's classic Whole-Genome SIMulator
  that ships alongside samtools (MIT). Always paired-end, always
  position-uniform, deliberately simple under a uniform sequencing-
  error model with configurable insert size, read length, and per-
  base error rate. Unlike ART (which models per-platform empirical
  error profiles), wgsim is the canonical "small + classic"
  simulator for fast smoke-testing of mappers and variant callers
  when realistic error spectra are not required. Single-binary
  subprocess shape: `wgsim` takes the reference and both output
  FASTQs as positional arguments (no stdout-redirect needed).
  `[bio.wgsim]` knobs: `reference` (FASTA; required), `output1`
  (FASTQ for read 1; required, non-empty), `output2` (FASTQ for
  read 2; required, non-empty — wgsim is paired-end only),
  `num_pairs` (≥ 1; required), `length1` (read 1 length, default
  70, ≥ 1), `length2` (read 2 length, default 70, ≥ 1),
  `fragment_size` (outer fragment length, default 500, > 0),
  `error_rate` (per-base error rate in `0.0..=1.0`, default 0.02 —
  typical Illumina baseline), `extra_args`. `prepare()` builds
  `wgsim -N <num_pairs> -1 <length1> -2 <length2> -d <fragment_size>
  -e <error_rate> <reference> <output1> <output2> [extras...]`.
  `collect()` reports `output1` and `output2` as `Tabular` artifacts
  ("wgsim simulated reads"). Probe via `find_on_path(&["wgsim"])`.
  Version range `1.0.0..2.0.0` (wgsim is versioned alongside the
  parent samtools 1.x line; the binary historically prints the
  matching samtools tag on startup). `bio.wgsim.simulate` ribbon
  capability.
- `valenx-adapter-badread` — Ryan Wick's long-read simulator with
  realistic Nanopore (and PacBio CLR) error profiles (GPL-3.0).
  Badread's per-platform error models are calibrated against actual
  sequencer output: random / chimeric / adapter / glitch read
  types, junk-read injection, identity drift, and length
  distributions that match what users see from a live flowcell.
  The de-facto choice for stress-testing long-read pipelines under
  realistic conditions. Single-binary subprocess shape with stdout-
  redirect: Badread writes its simulated FASTQ to stdout (no `-o`
  flag), so `run()` borrows MAFFT's stdout-redirect-to-file pattern
  — spawn the child directly, attach stdout to a `File` via
  `Stdio::from(file)`, stream stderr through the line handler.
  `[bio.badread]` knobs: `reference` (FASTA; required), `output`
  (FASTQ output path; required, non-empty), `quantity` (Badread's
  `--quantity` literal — one or more decimal digits followed by an
  optional `K` / `M` / `G` / `T` SI suffix, e.g. `"100M"` for 100
  megabases or `"5G"` for 5 gigabases; validated via the
  `is_valid_quantity` helper), `error_model` (one of
  `"nanopore2018"` / `"nanopore2020"` / `"nanopore2023"` /
  `"pacbio2016"`; required — selects the per-platform error profile
  baked into the Badread distribution), `identity_mean` (read
  identity mean as a percentage in `0.0..=100.0`; default 87.5),
  `length_mean` (read length mean in bases; default 15000.0,
  > 0.0), `length_sd` (read length stddev in bases; default
  13000.0, > 0.0), `extra_args`. `prepare()` builds `badread
  simulate --reference <reference> --quantity <quantity>
  --error_model <error_model> --identity <identity_mean> --length
  <length_mean>,<length_sd> [extras...]` → stdout, captured to
  `output` via the MAFFT-style stdout-redirect pattern. `collect()`
  reports `output` as a single `Tabular` artifact ("Badread
  simulated reads"). Probe via `find_on_path(&["badread"])`. Version
  range `0.4.0..1.0.0` (the long-running 0.4.x stable series; a
  1.0 cut hasn't happened yet but the upper bound reserves room for
  it). `bio.badread.simulate` ribbon capability.

Each adapter wired into `valenx-app::init_registry`.

**No canonical-type, format-reader, or CLI changes.** All three
adapters consume user-supplied reference FASTAs (the existing
`valenx_bio::format::fasta` reader already inspects sequence count +
identifiers + alphabets) and emit FASTQ files that the existing
Phase 18 `valenx-fastq` CLI inspects for record count, base quality
distributions, and read-length statistics. The unchanged
`Results.artifacts` collection model surfaces every emitted FASTQ +
ART alignment record directly. A first-class read-simulation
provenance type — recording which simulator produced which FASTQ
under which error model — defers to a future phase along with
simulator-aware pipeline stitching.

Three new `valenx-init` templates ship: `art` with alias
`art-illumina` (`art-simulate`), `wgsim` (`wgsim-simulate`), and
`badread` (`badread-simulate`). Cross-binary roundtrip test sweeps
all 77 templates clean.

Adapter inventory: 81 of 82 fully live (only `occt` remains stub-
only).

What's not in this phase: DWGSIM (Nils Homer's wgsim fork with
structural-variant injection; sister to wgsim; defer to Phase 31.5),
pIRS (BGI's profile-based Illumina simulator; sister to ART; defer),
InSilicoSeq (HMM-based Illumina + ONT simulator; defer), Mason
(UCSC's single-binary read simulator covering Illumina + 454;
defer), CuReSim (PCR / amplicon-aware simulator; defer), pbsim2 /
pbsim3 (PacBio HiFi simulator, sister to Badread; defer to 31.5),
NanoSim (Nanopore simulator with model training; different shape —
the per-flowcell training step requires a separate phase).
Out of scope.

The full plan lives at
`docs/superpowers/plans/2026-04-30-read-simulators.md`.

---

## [0.1.0-alpha.1] — 2026-04-25

First tagged version. Pre-alpha — no shipping installer, but the
end-to-end CFD workflow loop is genuinely usable on the command line
today: load a `.valenx` project, click a case, hit Run, watch the
mesh paint by field value with a colour-bar legend and a time-step
slider for transient snapshots.

This release covers everything in the rebuild from clean-slate Rust
through the results-rendering arc. See [STATUS.md](./STATUS.md) for
a snapshot of what works end-to-end and [QUICKSTART.md](./QUICKSTART.md)
for the five-minute walkthrough.

### Highlights

- **11 live adapters** covering every physics-domain phase 1-9:
  OpenFOAM (4 solvers), gmsh, FreeCAD, CalculiX (5 analyses), Elmer
  (steady + transient heat), Cantera, LAMMPS, openEMS, PyBaMM,
  MuJoCo, preCICE.
- **End-to-end visual results loop** for OpenFOAM: VTU parser →
  canonical `Mesh` + `Field` → auto-`collect()` after run →
  auto-load mesh into viewport → field-coloured wireframe overlay
  → colour-bar legend → clickable field picker → time-series
  scrubbing slider.
- **Workflow loop** for any live adapter: click-to-run / click-to-
  prepare-only / open-workdir-in-host-browser / run-from-prepared-
  workdir / per-case run-history badges / adapter-status badges.
- **Transient time-stepping** in three adapters (OpenFOAM /
  CalculiX / Elmer) using a consistent `[<equation>.transient]`
  block in `case.toml`.
- **318 passing tests**, zero clippy warnings, `cargo check
  --workspace --all-targets` clean.

### Known gaps

- Non-OpenFOAM adapters don't surface `Results.fields` yet —
  CalculiX `.frd` and Elmer `.vtu` need their own collect() wiring.
- Wireframe overlay only; filled-triangle field rendering is a
  follow-up (needs surface extraction + wgpu shader update).
- Probe-at-point / streamlines / vector glyphs are documented as
  Phase 10+ polish.
- 7 alt-capability scaffolds (OCCT / Netgen / Code_Aster /
  OpenRadioss / SU2 / GROMACS / Meep) stay as honest scaffolds —
  their primary partners cover the same physics, and rushing them
  would erode the "every live adapter actually works" promise.
- Phases 11-16 (HPC / optimization / ML / plugins / enterprise /
  stewardship) are planning docs only.

See [POLICIES.md](./POLICIES.md) for the SemVer / pre-alpha contract.

---

### Added (rolled into 0.1.0-alpha.1)

**App: time-series stepping for transient runs:**

Step 9 of the results-rendering arc — closes the gap that's been
sitting since the OpenFOAM transient solvers landed. Multi-timestep
runs (pimpleFoam writing `cavity_100.vtu`, `cavity_200.vtu`,
`cavity_300.vtu`, …) used to dump every snapshot into the Field
catalog but only ever render the first one. Now users can scrub
through the time series.

The pieces:

- **`selected_time_index: usize`** on `ValenxApp`. Driven by a new
  slider in the Results pane that shows step N-of-M plus the
  decoded `TimeKey` (`steady` / `iter 500` / `t=0.0050 s`).
- **Field-overlay resolver** picks the field at `selected_time_index`
  within the chosen field's `time_series`. Index gets clamped each
  frame so the slider can't outrun reality after a project switch.
- **Auto-reset on field switch** — clicking a different field name
  in the Results pane resets the slider to step 0. Different fields
  may have different snapshot counts, and a stale index would have
  dropped the user mid-series silently.
- **Colour-bar legend** now shows the timestep label below the field
  name when the field isn't steady. Users see exactly which snapshot
  they're looking at without leaving the viewport.
- **Slider hidden** when the field has only one snapshot (steady
  runs) so the pane doesn't grow noise.

`format_time_key(TimeKey) -> String` is the new helper rendering
each variant: `Steady` → `"steady"`, `Iteration(n)` → `"iter N"`,
`Time { value, units }` → `"t=0.0050 s"` with the SI suffix from
`Units::display`. Used by both the slider label and the legend.

End-to-end: run a pimpleFoam case → 600 timesteps land in the
catalog → slider shows "Time step 1 of 600 — iter 1" → drag the
slider, the wireframe re-paints with the new snapshot's field
values + the legend's min/max + timestep update accordingly.

Tests: 315 → 316 (+1: format_time_key across Steady / Iteration /
Time-with-units variants; selected_time_index gets a default-state
assertion in the existing default-state test).

**App: clickable field picker in the Results pane:**

Step 8 of the results-rendering arc. The previous commit auto-picked
the first scalar OnNode field for the wireframe overlay, but users
had no way to switch — `p` was the default and they were stuck
with it even when the catalog also held `T`, `Ux`, `Uy`, `Uz`, `k`,
`omega`, etc.

The field-name list in the Results pane is no longer read-only:

- Renderable fields (scalar, OnNode, length matches mesh node count)
  are now `selectable_label` rows. Clicking one sets it as the
  viewport's overlay; the wireframe re-paints next frame.
- Non-renderable fields (vector U, tensor stress, per-cell data)
  still appear in the list but stay `weak`-styled with a hover
  tooltip explaining "vector / tensor / cell-data fields aren't
  renderable on the wireframe overlay yet" — so users see the data
  exists without thinking the click would work.
- Click-the-already-selected-field clears the explicit pick. The
  overlay falls back to the auto-default.

The selection lives in a new `selected_field_name: Option<String>`
field on `ValenxApp`, queried first by the field-overlay resolver.
A stale selection (e.g. user picked `p` then loaded a different
project that doesn't have `p`) silently falls back to auto rather
than blanking the overlay — keeps the post-run "I see colours"
moment intact across project switches.

Tests: still 315 — the new `selected_field_name` field gets a
default-state assertion in the existing
default_state_has_no_prepared_job test.

**App: colour-bar legend in the viewport:**

Step 7 of the results-rendering arc — completes the read-the-colour
loop. The previous commit painted wireframe edges by field value,
but users had no way to know which field was being shown or what
the value range was. The new legend in the bottom-right corner
fills both gaps.

Layout:

- Field name above the strip (top label).
- 14×140 px vertical gradient strip with 32 stacked slices of the
  cool-to-warm ramp — top = max (warm), bottom = min (cool),
  matching scientific-plot convention.
- Min / max value labels alongside the strip's bottom and top.
- Semi-transparent dark backing card (`rgba 20,22,26,200`) so the
  legend stays readable over busy meshes.

A new `format_field_value(f64)` helper picks the readable shape
per value: integer-looking for whole numbers, up to 4 decimal
places (trimmed) for small fractionals, scientific notation
outside `[1e-3, 1e6)` so labels never overflow the 90 px label_room
budget. Covers most CFD pressures (10⁵ Pa) and thermal
temperatures (10² K) without scientific noise, falls back to `e`
notation for residuals and species mass fractions.

The legend renders only when `state.field_overlay.is_some()` —
no overlay, no legend, no visual debt for users not running cases.

Tests: 314 → 315 (+1: format_field_value covers integer / small
fractional / sub-1e-3 / over-1e6 shapes).

**App: field-coloured wireframe overlay in the viewport:**

Step 6 of the results-rendering arc — the first one that visibly
changes what the user sees on screen. After a run completes and
the field catalog is populated, the mesh wireframe in the viewport
now colour-codes each edge by the average of its two endpoints'
field values via a five-stop **cool-to-warm** divergent ramp.

The pieces:

- **`valenx_fields::colormap`** — new module with
  `cool_to_warm(t: f32) -> [u8; 3]` and
  `cool_to_warm_in_range(v, min, max) -> [u8; 3]`. Five stops
  (dark blue → light blue → near-white → orange → dark red),
  piecewise-linear interpolation, clamped at the endpoints.
  Degenerate ranges (`min == max`) collapse to the midpoint
  colour rather than dividing by zero.
- **`viewport::ViewportState.field_overlay: Option<&Field>`** —
  new optional input. When `Some` and the field is scalar +
  OnNode + matches the mesh node count, every edge in
  `draw_mesh_wireframe` gets a per-line stroke colour from the
  ramp. Otherwise the existing flat-blue wireframe stays.
- **App-side selection** — `update()` picks the first scalar
  OnNode field from `last_run_results` whose data length matches
  the loaded mesh's node count. Most CFD runs include `p`
  (pressure), most thermal runs `T` (temperature) — both satisfy
  the filter and give the user a colour-coded wireframe
  immediately after a successful run.

This is intentionally the simplest viable path. Filled-triangle
field rendering needs surface extraction (boundary faces of 3D
elements → triangle list) plus a wgpu shader update with a
per-vertex field attribute and colour-bar legend. Both arrive in
follow-ups; this commit proves the data path works end-to-end and
gives users their first "I can see my CFD result" moment.

End-to-end: load a project → click cfd-transient → Run → solver
writes `cavity_500.vtu` → collect parses fields → mesh auto-loads
→ wireframe edges paint cool-to-warm by `p` (or whatever the
first scalar field is). All without manual UI clicks.

Tests: 308 → 314 (+6 in `colormap`: endpoints match stops,
midpoint near white, out-of-range clamps, cool-to-warm-half
brightness monotonic, degenerate range returns midpoint, normalised
range maps endpoints exactly).

**App: auto-load OpenFOAM mesh into viewport when a run finishes:**

Step 5 of the results-rendering arc. The previous commit got the
parsed Field catalog into app state; this one gets the matching
geometry visible in the viewport so the user sees what they
simulated, without clicking anything.

When an `openfoam` run completes, `on_run_finished` walks the
workdir for `.vtu` files, picks the lexicographically last one
(which is the highest-time-step snapshot for OpenFOAM's
`<case>_<N>.vtu` naming), parses its mesh half via the new
`load_mesh_from_vtu` helper, and calls `apply_mesh()` to put the
geometry on screen.

Two helpers do the lifting:

- `latest_vtu_in_workdir(workdir) -> Option<PathBuf>` walks the
  whole workdir tree (not just the top level — OpenFOAM nests
  VTKs under `VTK/`) and returns the latest `.vtu` by sort order.
- `load_mesh_from_vtu(path) -> Result<Mesh, String>` reads, parses,
  converts, and returns just the canonical `Mesh` half. The Field
  half is dropped — `Results.fields` already owns those, populated
  by `collect()`.

`load_mesh` got refactored to extract an `apply_mesh(mesh,
source_path)` helper that does the post-parse work (recompute_stats,
quality_report, status update, frame_current_mesh). `load_mesh`
is now "parse + delegate"; `apply_mesh` is the entry point for
callers that already have a `Mesh` in hand. Same UX, cleaner
internal split.

End-to-end visible result: run a `pimpleFoam` case → solver writes
`cavity_500.vtu` → run finishes → viewport snaps to a fresh tet
mesh framed at the right zoom, with the field catalog populated
in the Results pane. The next commit paints triangles by field
value.

Tests: 305 → 308 (+3: latest_vtu_in_workdir picks the highest
case_N, returns None on empty dirs, load_mesh_from_vtu round-trips
a one-tet through to a real `Mesh` with the right ElementType).

**App: auto-collect after run + show field catalog in Results pane:**

Step 4 of the results-rendering arc. Before this commit, every
adapter's `collect()` was unreachable from the run pipeline — the
worker thread sent `Finished(report)` and quit. Users would have
seen a "fields supported" message but no actual fields, because
`collect()` was only ever called from tests.

The worker now calls `adapter.collect(&prepared)` automatically
right after a successful `run()` and ships the result via a new
`RunEvent::Collected(Box<Results>)` variant. Failures of `collect()`
land as a `LogLine` warning rather than a separate event so the
"did this run finish" check stays simple — partial results are
still better than no results.

App side:

- New `last_run_results: Option<Box<Results>>` field on `ValenxApp`,
  populated when the worker sends Collected.
- `pump_run_events`' Finished handler no longer sets `finished =
  true` — that was forcing the run handle to drop before Collected
  could arrive. The natural `Disconnected` detection (when the
  worker thread exits) cleans up after both events have been
  processed.
- Results pane gains a field catalog summary when the run produces
  one: `"6 field samples in 3 unique fields"` followed by a per-
  field line like `"· U  (3 timesteps)"`. Hidden when the catalog
  is empty (preCICE meta-runs, etc.) so the pane doesn't grow noise.

End-to-end path is now live for the OpenFOAM stack:

> Run pimpleFoam → foamToVTK writes cavity_*.vtu → worker calls
> collect() → parser populates Results.fields → app receives
> Collected event → field catalog renders in Results pane.

Next commit in the arc is the wgpu colour-bar that paints triangles
by field value — the data is now there waiting for it.

Tests: still 305 passing — the new last_run_results field gets a
default-state assertion in the existing default_state_has_no_prepared_job
test rather than a new test. The NoopAdapter test fixture in run.rs
also gained a real (empty) collect() impl since the worker now calls
it on every successful run.

**OpenFOAM `collect()` parses `.vtu` artifacts into the Field catalog:**

Step 3 of the results-rendering arc. The previous two commits gave
us a parser + canonical-type converter; this commit hooks them into
OpenFOAM's `collect()` so finished runs actually populate
`Results.fields` instead of just listing artifact paths.

After every run, `collect()` now:

1. Discovers all artifacts (existing behaviour).
2. Filters for `.vtu` files (case-insensitive extension match).
3. Sorts by path so OpenFOAM's `<case>_<N>.vtu` naming gives
   chronological order.
4. For each file: `parse_ascii` → `to_canonical(stem)` → insert
   every Field into `results.fields` with a derived `TimeKey`.
5. Per-file failures (broken VTK, partial write) are
   logged-and-skipped — a single bad file shouldn't wipe out
   everything else.

`time_key_from_filename(stem)` derives the time index by splitting
on the last `_` and parsing the trailing token as `u64`. So
`cavity_500.vtu` → `TimeKey::Iteration(500)`,
`flow_around_cylinder_120.vtu` → `TimeKey::Iteration(120)`,
unparseable names fall back to `TimeKey::Steady`. This makes
multi-timestep transient runs land in the catalog as a real time
series the report layer can iterate.

The end-to-end path is now:

> Run pimpleFoam → `foamToVTK` writes `cavity_*.vtu` → `collect()`
> parses each → `Results.fields` carries `(p, U, T, …)` per
> timestep → ready for the wgpu colour-bar pipeline (next commit).

Tests: 303 → 305 (+2: collect_loads_vtk_fields_into_catalog spins
up a fake workdir with a `cavity_500.vtu`, runs collect(), asserts
the `p` field lands at `TimeKey::Iteration(500)` and as `Scalar` /
`OnNode`; time_key_from_filename_handles_canonical_shapes covers
`<case>_<N>` / multi-underscore / no-trailing-int / empty edges).

**Foundation: `.vtu` → canonical `Mesh` + `Field` converters:**

Step 2 of the results-rendering arc. The previous commit added the
ASCII parser; this commit converts the resulting `VtuData` into the
canonical types the rest of Valenx uses:

- `VtuData::to_canonical(mesh_id) -> (valenx_mesh::Mesh, Vec<Field>)`
  is the one-shot entry point. Coordinates land as `Vector3<f64>`
  nodes; cells are grouped by canonical `ElementType` into one
  `ElementBlock` each (so a mixed tet+hex VTU produces two blocks);
  point and cell fields land as a flat `Vec<Field>` with `Location::
  OnNode` / `OnCell` set appropriately.
- `vtu_to_element_type(VtuCellType)` is the public mapping function
  for callers that want to walk cells themselves. Vertex + Unknown
  return `None` and are silently skipped by `to_canonical` —
  dropping them keeps the mesh consistent (sum-of-block-counts
  matches surviving cell count).
- Field component count picks the canonical kind: 1 → `Scalar`,
  3 → `Vector { dim: 3 }`, 9 → `Tensor { rows: 3, cols: 3 }`,
  anything else → `Vector { dim: n }` so weird-but-real data isn't
  silently lost.
- `Field.range` is cached at conversion time (min/max across every
  f64 in the buffer) so the colour-bar widget can default-stretch
  without re-scanning the data.

VTU doesn't carry units or time-step indices; we default to
`DIMENSIONLESS` + `TimeKey::Steady`. Adapters that know the physics
(CFD = velocity/pressure/temperature, FEA = displacement/stress)
can re-stamp the units in a follow-up pass before the field reaches
the report layer.

`valenx-fields` now depends on `valenx-mesh` (acyclic — mesh stays
standalone). Dep direction matches the conceptual one: fields are
defined on meshes, so fields knows about mesh.

Tests: 298 → 303 (+5: one-tet round-trip with both mesh and fields,
mixed tet+hex grouping, vertex+unknown skipping, vtu_to_element_type
across the supported set, field_range edges including empty / single
/ constant / mixed-sign).

**Foundation: `.vtu` ASCII parser in `valenx-fields`:**

The first piece of the "viewport renders simulation results" arc.
The `vtu` module parses the ASCII `UnstructuredGrid` shape that
OpenFOAM's `foamToVTK`, Elmer's `Post File`, and most other VTK
writers emit by default. Output is a plain `VtuData` struct with:

- `VtuMesh` — points + flat connectivity + per-cell offsets +
  per-cell type codes (mapped to a `VtuCellType` enum: Tet, Hex,
  Tri6, Tet10, etc., with `Unknown(code)` for anything outside
  the supported set so callers can decide whether to skip).
- `point_fields: Vec<VtuField>` — named fields stored at mesh
  nodes (e.g. `U`, `p`, `T` from CFD), each with a component
  count (1 = scalar, 3 = vector) and a flat data buffer.
- `cell_fields: Vec<VtuField>` — same shape, but per-cell.

What's intentionally NOT supported (returned as
`ParseError::Unsupported` rather than wrong data):

- `<AppendedData …>` (binary appended sections)
- `compressor=` (zlib/lz4 compressed payloads)
- `PUnstructuredGrid` (parallel decomposition)
- `format != "ascii"` on any DataArray (binary / base64)

These all need a real VTK library; the contract here is "ASCII in,
real data out, nothing else." The error message points users at
the right gap so they know to convert with `vtkXMLUnstructuredGridReader`
or to ask the writer for ASCII.

No new dependencies — the format is rigid enough that hand-rolled
`find` / `strip_prefix` parsing handles it cleanly. Keeps
`valenx-fields`' license perimeter tight.

Tests: 9 new (one-tet round-trip with scalar `p` and vector `U`,
self-closed `<PointData/>`, `<CellData>` parsing, count-mismatch
diagnostic, all four Unsupported variants, cell-type code round-
trip across the supported set).

**App: adapter status badges in the case browser:**

The case browser used to render every case as `<name> · <solver>`
with no indication of whether the solver's adapter was actually
ready to run. Users had to visit the Adapters panel and cross-
reference the solver string to figure out if "click Run" would
work — annoying for projects with many cases.

Each row now starts with a coloured ● badge matching the convention
already used in the Adapters panel:

| Colour       | Status        |
|--------------|---------------|
| green        | Ready         |
| gray         | Missing       |
| yellow       | Outdated      |
| red          | Broken        |
| blue-purple  | Disabled      |
| dark gray    | Unregistered  |

Hover gives the full reason: ``` `openfoam` adapter: Ready ``` /
`` `calculix` adapter: Missing `` / etc. The new "Unregistered"
status is for solver strings whose adapter id isn't in the registry
at all (e.g. a typo or a not-yet-shipped adapter).

The selectable label and click/double-click behaviour are unchanged
— the badge is purely informational. Implementation is a small
horizontal layout per row: `[badge] [selectable name + solver]`.

No new tests — the change is purely visual. The compile-clean +
clippy-clean verification on the existing 289-test suite is the
contract.

**App: run from prepared workdir (skip prepare):**

Closes the prepare → edit → run workflow loop. Before this commit,
users could click "Prepare", "Open in file browser", and edit dicts
by hand — but then the only way to run was via "Run selected case",
which calls `prepare()` again and overwrites their edits.

`run_from_prepared_workdir()` re-uses the `PreparedJob` captured from
the most recent successful `prepare_selected_case()` and spawns the
solver against it directly, skipping the prepare step entirely. The
user's hand-edits to the dicts in the workdir survive.

`run.rs` got a new `spawn_prepared(adapter, prepared)` entry point.
Internally it shares the run-loop body with `spawn(adapter, case,
workdir)` via a private `RunSpec::{Fresh, Prepared}` enum, so the
progress / log / cancel plumbing stays single-source.

A new `last_prepared_job: Option<(String, PreparedJob)>` field on
`ValenxApp` stashes the adapter id alongside the job (the spawn
needs to look the adapter back up in the registry; the id is the
key). Cleared when prepare fails so users don't accidentally run
a stale job.

Wired through:
- Run menu: "Run from prepared workdir" between Prepare and the
  legacy "Run first case", with a tooltip explaining the workflow.
- Command palette: `run.from-prepared-workdir` (palette: 22 → 23).
- Disabled when no prepare has run yet, or when a run is in progress.

Tests: 287 → 289 (+2: run-from-prepared-without-prepare rejects
cleanly, default state has no prepared job).

**App: open prepared / run workdir in host file browser:**

Natural follow-up to the prepare-selected-case workflow. After the
user clicks "Prepare" (or finishes a run), the right-pane Results
section now shows an "Open in file browser" button next to each
workdir path. The button shells out to the platform-native launcher:

| Host    | Launcher       |
|---------|----------------|
| Windows | `explorer.exe` |
| macOS   | `open`         |
| Linux   | `xdg-open`     |

The launcher is spawned and detached — Valenx doesn't wait on the
child or kill it. If the launcher itself can't spawn (e.g.
headless Linux without `xdg-open`), the error includes the
workdir path so users have a fallback they can copy-paste.

The run pipeline got a new `last_run_workdir: Option<PathBuf>`
companion to `last_prepare_workdir`. It's set when the run handle
drops at the end of `pump_run_events`, regardless of whether the
run succeeded or failed — failed runs typically still leave a
partial dict tree + log on disk that users want to dig into.

Wired through:
- Right-pane Results: "Open in file browser" button under both the
  Prepare and Run workdir text edits, with hover tooltips.
- Command palette: `view.open-prepare-workdir` +
  `view.open-run-workdir` (palette array grew 20 → 22).

Tests: 285 → 287 (+2: open-prepare-without-prepare rejects cleanly,
open-run-without-run rejects cleanly).

**App: prepare-selected-case workflow (no execute):**

A new `prepare_selected_case()` method on `ValenxApp` mirrors
`run_selected_case()` up to adapter resolution but stops after
calling `adapter.prepare()` — no solver process is spawned.

Useful for two cases:

1. The user has the case but not the underlying tool installed and
   wants to see what the adapter would write. OpenFOAM, CalculiX,
   Elmer all emit their dict / .inp / .sif trees during prepare;
   the user can open the workdir and inspect (or hand off to a
   third-party tool, or copy onto a remote HPC node).

2. The user wants to inspect / edit the generated files before
   running. Click "Prepare", inspect the workdir, then run
   manually with the modified files.

The workdir lands in `std::env::temp_dir()` with a
`valenx-prepare-<case>-<unix>` prefix to keep it distinct from
run workdirs (`valenx-run-…`). The path is stashed in a new
`last_prepare_workdir: Option<PathBuf>` field and surfaced in the
right-pane Results section as a copyable read-only text edit.

When prepare fails — e.g. `ToolNotInstalled` because `gmshToFoam`
isn't on PATH — we still set `last_prepare_workdir` because most
adapters write the dict tree BEFORE the tool lookup. The error
message includes the workdir path so users can find what was
generated despite the failure.

Wired through:
- Run menu: "Prepare selected case (no execute)" between Run and
  Cancel, with a tooltip explaining the use cases.
- Command palette: `run.prepare-selected-case` →
  "Run: Prepare selected case (no execute)" (palette array grew
  19 → 20).

Tests: prepare-without-project rejects cleanly with a structured
error pointing at the missing project, plus a default-state
regression on the new field.

**Elmer transient heat equation — BDF-2 marching with initial conditions:**

The Elmer adapter previously hardcoded `Simulation Type = Steady State`
in every emitted SIF. A new `TimeMode` enum on `ElmerInput` mirrors
the OpenFOAM / CalculiX shape: `Steady` (default) keeps the existing
behaviour; `Transient { end_time, delta_t }` selects BDF-2 marching.

Triggered by an optional `[heat.transient]` block in `case.toml`:

```toml
[heat.transient]
end_time = 60.0
delta_t  = 0.1
```

The writer dispatches on `time.simulation_type()`:

- Steady SIF: `Simulation Type = Steady State`, no `Timestep …` lines
  (older Elmer versions warn when those land in steady configs;
  unconditional emission was a small landmine).
- Transient SIF: `Simulation Type = Transient` + `Timestepping
  Method = BDF` + `BDF Order = 2` + `Timestep Intervals = ceil(end_time
  / delta_t)` + `Timestep Sizes = delta_t`. Round-up on the interval
  count so we always reach `end_time` rather than under-shoot it.

Optional `initial_temperature` (or `T0`) on the `[heat]` block emits
an `Initial Condition 1` SIF block with `Temperature = T0`, and the
`Body 1` block grows an `Initial Condition = 1` reference. Without
that reference Elmer ignores the IC entirely and starts from zero,
which is rarely physical for transient cool-down / heat-up studies.

Tests: 8 new (transient parsing, partial-config rejection,
`delta_t > end_time` rejection, T0 alias, BDF intervals exact +
round-up, IC block shape + body reference, steady-with-T0 still
emits the IC block).

**CalculiX transient FEA — `*DYNAMIC` + transient `*HEAT TRANSFER`:**

The CalculiX adapter previously only emitted steady analyses
(`*STATIC`, `*FREQUENCY`, `*HEAT TRANSFER, STEADY STATE`). Two
transient variants now join them:

- **`LinearDynamic`** (`*DYNAMIC`) — implicit linear dynamic with
  the Hilber-Hughes-Taylor integrator. Good for impact, drop-tests,
  harmonic excitation. Selected by `analysis = "linear-dynamic"`
  (or `"dynamic"`, `"transient-dynamic"`) in `[structural]`.
- **`ThermalTransient`** (`*HEAT TRANSFER` without the
  `, STEADY STATE` qualifier) — time-marches the temperature field.
  Good for cool-down / heat-up studies. Selected by
  `analysis = "thermal-transient"` (or `"heat-transient"`).

The existing `Step` struct already carried `time_total` +
`time_increment`; the only writer change is that `*DYNAMIC` and
transient `*HEAT TRANSFER` now also emit the `<delta_t>, <total_t>`
data row that `*STATIC` was already writing. Two new helper methods
on `AnalysisKind` keep the dispatching honest:

- `needs_increment_line()` — true for `*STATIC` + `*DYNAMIC` +
  transient `*HEAT TRANSFER`; false for `*FREQUENCY` and steady
  `*HEAT TRANSFER` (which would misparse the data row).
- `is_time_marching()` — true for `*DYNAMIC` and transient
  `*HEAT TRANSFER` (real-time march); false for `*STATIC` (uses
  pseudo-time for load ramping) and `*FREQUENCY` (one-shot eigen).

Tests cover both new analysis paths against fresh `.inp` decks,
plus the steady-thermal regression: the data row after `*HEAT
TRANSFER, STEADY STATE` must NOT be a numeric line, otherwise CCX
treats it as a property and fails.

**OpenFOAM transient convergence semantic — `None` instead of `Some(false)`:**

A leftover from the transient solver work: `RunReport.converged` was
always set to `Some(last_residual_below(...))`, which makes sense for
steady simpleFoam but is meaningfully wrong for pimpleFoam / icoFoam.
Transient runs march to a fixed `end_time` regardless of where the
residuals land at the final step — saying `Some(false)` after a
successful run mislabelled it as "did not converge" in the UI.

The run loop now detects transient mode by watching for fractional
`Time = …` values during streaming (steady simpleFoam never emits
those), and reports `converged: None` for transient runs. The UI
already renders `None` as "convergence unknown", which is the honest
answer — there is no convergence criterion for a transient run that
ran to completion.

A new `is_transient_time(f64)` heuristic with unit tests handles the
edge cases (whole-1.0 vs whole-1.0+ε vs sub-1.0).

**OpenFOAM transient log parser — honest progress for sub-second times:**

The OpenFOAM log parser used to expose the `Time = …` marker as
`LogSignal::Iteration(u64)`, casting the f64 value to u64 in the
process. That was correct for steady simpleFoam (`Time = 5` →
`Iteration(5)`) but truncated every transient marker to zero
(`Time = 0.0005` → `Iteration(0)`), which pinned the progress bar
at zero and snapped every residual sample on the chart to step 0.

The variant is now `LogSignal::Time(f64)` and the run loop tracks
two values: a monotonic step counter for `ResidualSample::iteration`
(canonical type still `u64`), plus the real time for the progress
label. Steady runs show `"Time = 250"`, transient runs show
`"Time = 0.0005"` — the same shape the solver itself printed. The
residual chart in `valenx-app` now plots at the real time too, so
transient runs land on a meaningful x-axis rather than collapsing
every sample onto step 0.

A new helper `format_time_label(f64)` makes the formatting decision:
integer-valued ≥ 1 → whole number, sub-second ≥ 1e-4 → trimmed
decimal, microsecond and below → scientific notation. The chart
reuses the parser's `intern_field` so the field allow-list stays
single-source.

**OpenFOAM transient solvers — pimpleFoam + icoFoam:**

The OpenFOAM adapter previously only emitted steady-state `simpleFoam`
cases. It now also handles two transient incompressible solvers:

- **pimpleFoam** — transient PIMPLE (merged PISO + SIMPLE). Accepts
  laminar or RANS turbulence, arbitrary time-step size, with `Final`
  pressure/velocity correctors set to `relTol = 0` for accuracy on
  the last sweep of every time step.
- **icoFoam** — strictly laminar transient PISO. Lighter dict tree
  (no `0/k`, `0/omega`, `0/nut`); the case parser refuses any RANS
  turbulence model with a clear error rather than silently dropping
  it.

The dict writer dispatches on a new `SolverKind` enum (`SimpleFoam` |
`PimpleFoam` | `IcoFoam`) and a `TimeMode` enum (`Steady` |
`Transient { end_time, delta_t, write_interval }`), parsed from the
`case.solver` string and an optional `[solve.transient]` block. Steady
cases reject a stray `[solve.transient]` block; transient cases accept
defaults (1 s / 1 ms / 100 ms snapshots) when the block is omitted.

Generated dicts switch shape per solver: `controlDict` writes real
seconds + `adjustableRunTime` for transient; `fvSchemes` swaps
`steadyState` for `Euler`; `fvSolution` emits a `SIMPLE` /
`PIMPLE` / `PISO` block as appropriate, drops `relaxationFactors`
for transient runs (the time derivative does the stabilisation).

`prepare()` resolves the requested solver binary first, falls back to
any OpenFOAM binary on PATH for a structured "case requested X but I
found Y" error, then keeps the dict files on disk so users can fix
PATH and retry without regenerating the case. `SOLVER_BINARIES`
expanded to include `pimpleFoam` + `icoFoam`.

A second fixture `tests/fixtures/minimal.valenx/cases/cfd-transient`
ships alongside `cfd-steady`. The minimal-project round-trip and the
adapter's prepare-against-fixture test cover both, so the new code
paths exercise on every CI run regardless of whether OpenFOAM itself
is installed.

**Phases 7 + 8 + 9 kickoff — PyBaMM + MuJoCo + preCICE:**

Every physics-domain phase (1 through 9) now has at least one live
adapter. The remaining scaffolds (OCCT / Netgen / Code_Aster /
OpenRadioss / SU2 / GROMACS / Meep) are alt-capability adapters
within already-covered phases — kept honest as scaffolds with
working probes rather than fake-live.

- **PyBaMM** (`valenx-adapter-pybamm`, Phase 7 battery): Python-
  driven, similar to Cantera. `case_input.rs` parses `[battery]`
  into a typed `BatteryInput` with `ModelKind` ∈ {Spm, Spme, Dfn},
  `Protocol` ∈ {CcDischarge, Cccv}, initial SOC, time horizon.
  `python_script.rs` emits `valenx_pybamm.py` that loads the
  chosen lithium-ion model + `ParameterValues("Chen2020")`,
  builds a `pybamm.Experiment` from the protocol, solves it,
  writes a decimated CSV time series + JSON summary.
- **MuJoCo** (`valenx-adapter-mujoco`, Phase 8 multibody): typed
  `DynamicsInput` with `ModelSource` ∈ {Mjcf, Urdf} (extension-
  sniffed), duration, optional timestep override, per-actuator
  constant control map, initial qpos / qvel. Script emits a
  full mj_step loop recording qpos / qvel / ctrl per sample into
  a JSONL trajectory. Accepts either `physics = "robotics"` or
  `"dynamics"`.
- **preCICE** (`valenx-adapter-precice`, Phase 9 coupling):
  meta-adapter. Parses `[coupling]` + `[[coupling.participant]]`
  arrays, stages `precice-config.xml` + each participant's case
  dir into the workdir, emits a `valenx_coupling.json` manifest
  the Phase-9-tail run orchestrator reads, and runs
  `precice-tools check` for config validation. Full
  concurrent-orchestration of participating solvers is defined
  in RFC 0007 and is the Phase 9 tail.

The app registry now lists **11 live adapters** (OpenFOAM, gmsh,
FreeCAD, CalculiX, Elmer, Cantera, LAMMPS, openEMS, PyBaMM,
MuJoCo, preCICE). Every physics-domain phase 1-9 has at least
one live adapter.

**Phases 4 + 5 + 6 kickoff — Cantera + LAMMPS + openEMS:**

Three new physics domains come online in one push, all following
the pattern the earlier adapters established (typed case input +
deterministic file emitter + subprocess run + structured artifact
collection).

- **Cantera** (`valenx-adapter-cantera`): Python-driven equilibrium
  calculator. `case_input.rs` parses `[chemistry]` with `Mechanism`
  ∈ {Bundled, External}, `Analysis` ∈ {EquilibriumTP, EquilibriumHP,
  EquilibriumUV}, `ThermoState` (T in K, P in Pa, Cantera-flavoured
  composition string). `python_script.rs` emits a deterministic
  `valenx_cantera.py` that runs `gas.equilibrate()` and writes a
  filtered mole-fraction JSON summary. `summary_parser.rs` is a
  tolerant serde deserialiser (including `final` → `final_`
  rename). `lib.rs::prepare()` stages external mechanism files,
  `run()` spawns `python3 valenx_cantera.py` through the shared
  subprocess runner, `collect()` attaches summary + script. 15 new
  tests.
- **LAMMPS** (`valenx-adapter-lammps`): subprocess adapter for
  classical molecular dynamics. `case_input.rs` covers `Units` ∈
  {Lj, Metal, Real, Si}, `BoundaryCondition` per-axis, typed
  `Initialization` ∈ {LjFccBox, ReadData}, `Potential` ∈ {LjCut,
  Eam}, `Ensemble` ∈ {Nve, NvtNose, NptParrinelloRahman}.
  `input_writer.rs` emits a self-contained LAMMPS deck with
  `units`, `lattice`, `pair_style`, `fix`, `dump`, `thermo`, `run`
  — everything needed for an NVE Lennard-Jones smoke case.
  `lib.rs::prepare()` stages external `read_data` / EAM files,
  `run()` spawns `lmp -in in.lammps`, `collect()` harvests
  `traj.lammpstrj`, `log.lammps`, `thermo.dat`. 8 new tests.
- **openEMS** (`valenx-adapter-openems`): subprocess adapter
  driving openEMS via a generated Octave script. `case_input.rs`
  parses `[em]` with `Domain::Box`, `Excitation` ∈ {Gauss, Sine},
  `BoundaryCondition` ∈ {Mur, Pml, Pec}, `Probe` list.
  `octave_script.rs` emits a self-contained `valenx_openems.m`
  calling `InitFDTD` + `SetGaussExcite` + `SetBoundaryCond` +
  `DefineRectGrid` + `AddProbe` + `WriteOpenEMS` + `RunOpenEMS`.
  `lib.rs::run()` dispatches between Octave and MATLAB based on
  the resolved binary, `collect()` walks the `sim/` subdirectory
  recursively for `.h5`, `.xml`, `.vtr/.vts/.vtu`, `.dat/.csv`
  artifacts. 9 new tests.

All three registered in `ValenxApp::init_registry` — the browser's
Adapters collapsible now lists **8 live probes** (OpenFOAM, gmsh,
FreeCAD, CalculiX, Elmer, Cantera, LAMMPS, openEMS).

**Phase 2 — prism layers in gmsh:**
- `valenx-adapter-gmsh::mesh_input::BoundaryLayer` + corresponding
  `[mesh.boundary_layer]` case.toml section (first-cell thickness,
  growth rate, layer count, optional surface names).
- The `.geo` writer now emits a matching `Field[1] = BoundaryLayer;`
  block with `hwall_n`, `ratio`, and a computed total `thickness`
  (geometric-series sum), then `BoundaryLayer Field = 1;` to
  activate it. Closes the last Phase 2 meshing checklist item.
- 3 new geo_writer tests covering presence, absence, and ratio=1.

**Phase 3 — Elmer heat-equation adapter lives:**
- `valenx-adapter-elmer` graduates from scaffold to a real
  steady-state heat equation driver.
  - `case_input.rs` parses `[heat]` section into a typed
    `ElmerInput` with `Equation`, `Material` (ρ / Cₚ / k),
    `BoundaryCondition` ∈ {Temperature, HeatFlux}, `Simulation`
    control, output basename. Either `physics = "fea"` or
    `"multi-physics"` is accepted.
  - `sif_writer.rs` emits a self-contained `case.sif` covering
    `Header` + `Simulation` + `Constants` + `Body` + `Solver`
    (with the BiCGStab+ILU0 linear-solver stack Elmer's HeatSolver
    uses by default) + `Equation` + `Material` + per-BC
    `Boundary Condition N`.
  - `lib.rs::prepare()` writes the SIF, shallow-copies the user's
    `mesh_dir` into the workdir so `Mesh DB "." "mesh"` resolves.
    `run()` spawns `ElmerSolver case.sif` through the shared
    subprocess helper with progress hints keyed to Elmer's stdout
    banners. `collect()` walks for `.vtu`, `.pvtu`, `.ep`,
    `.result`, `.sif`, and `.log` artifacts.
- Registered in `ValenxApp::init_registry`, so the browser's
  Adapters collapsible now shows 5 live probes.
- 12 new tests across the adapter (4 case parser, 4 SIF writer,
  3 lib sanity, 1 copy_dir_shallow).

**Phase 2 — FreeCAD adapter lives + E2E integration test:**
- `valenx-adapter-freecad` graduates from scaffold to real import.
  New `case_input.rs` (typed `GeometryImportInput` with source path,
  format detection, exports, STL deviation, feature-tree toggle),
  `python_script.rs` (deterministic `FreeCADCmd` Python emitter
  that opens the source file, walks `ActiveDocument.Objects`,
  exports to STL / BREP / STEP as requested, and writes a JSON
  summary), and `summary_parser.rs` (serde parser for the emitted
  JSON with tolerant defaults for partial FreeCAD runs). `prepare()`
  stages the source alongside the script, `run()` spawns
  FreeCADCmd through the shared subprocess helper with coarse
  progress hints, `collect()` parses `summary.json` and attaches
  every produced artifact to `Results`.
- **End-to-end pipeline test** at
  `crates/valenx-app/tests/pipeline_e2e.rs` drives gmsh → OpenFOAM
  prepare in one test, tolerant of three machine states (neither
  tool installed, gmsh only, both installed) so it runs meaning­
  fully in CI and on dev boxes.

**Phase 3 kickoff — CalculiX adapter lives:**
- `valenx-adapter-calculix` graduates from scaffold to linear-
  static FEA. New `case_input.rs` parses the `[structural]`
  section (typed `LinearStaticInput` with `AnalysisKind` ∈
  {LinearStatic, Modal, Thermal}, `Material`, `Boundary`, `Load`,
  `Step`, `OutputField`). `inp_writer.rs` emits a self-contained
  Abaqus-flavoured `job.inp` from a canonical `valenx_mesh::Mesh`
  plus the typed input — `*NODE` / `*ELEMENT` (C3D4 / C3D10 /
  C3D8 / C3D20 / C3D6 / C3D5 / CPS3 / CPS6 / CPS4), `*MATERIAL`
  + `*ELASTIC` + optional `*DENSITY`, `*SOLID SECTION`,
  `*BOUNDARY`, `*STEP` + card for the analysis kind, `*CLOAD`,
  and `*NODE PRINT` / `*EL PRINT` / `*NODE FILE` / `*EL FILE`
  based on requested output fields. NSET names are sanitised;
  unknown element types become `**` skip comments instead of
  broken blocks.
- `run()` spawns `ccx -i job` through the shared subprocess
  runner with progress hints mapped to ccx's stdout banners.
- `collect()` walks the workdir for `.frd`, `.dat`, `.sta`,
  `.cvg`, `.12d`, and `.inp` artifacts.

**Phase 2 — closing the gmsh → OpenFOAM pipeline:**
- `valenx-adapter-openfoam::prepare()` now looks for a `mesh.msh`
  in either the case directory or the workdir and, when found,
  shells out to `gmshToFoam` to materialise `constant/polyMesh/`
  before emitting its solver dicts. Skipped silently when a
  polyMesh already exists or no .msh is present — that's the
  "users have their own meshing pipeline" escape hatch. Missing
  `gmshToFoam` surfaces as `AdapterError::ToolNotInstalled`;
  non-zero `gmshToFoam` exits surface as `AdapterError::Run` with
  the stderr tail. 3 unit tests covering all three paths.

**Phase 2 — shared subprocess runner + auto-load UX:**
- **`valenx_core::subprocess`** — one canonical
  `prepare → spawn → stream → cancel → finalize` loop shared by
  every Subprocess-mode adapter. Adapters pass in a closure
  returning `Hint { progress, warning }` to plug in their own
  per-line parsing. OpenFOAM's driver shrank from 170 lines to
  ~40; gmsh's from 150 to ~20.
- **Auto-load mesh after a gmsh run.** `RunHandle` now carries
  `adapter_id` + `workdir`; `ValenxApp::on_run_finished` loads the
  produced `mesh.canonical.json` (or fallback `.msh`) straight
  into the viewport so the user doesn't have to drag it in.

**Phase 2 kickoff — gmsh adapter goes live:**
- `valenx-adapter-gmsh::prepare()` emits a real `.geo` from the
  `[mesh]` section of `case.toml`: procedural `Box` and `Sphere`
  primitives via the OpenCASCADE factory, or `Merge` for external
  STL / BRep / STEP / IGES files (copied into the workdir so the
  emitted path is relative). New modules:
  - `mesh_input.rs` — typed `MeshSpec` with `Domain`, `MeshSizes`,
    `Algorithm2D`, `Algorithm3D`, `MeshDim`, lenient string parsing
    for algorithm names.
  - `geo_writer.rs` — `.geo` emitter with all relevant `Mesh.*`
    options, OpenCASCADE primitives, Merge + `CreateTopology`
    volume wrapping, and the trailing `Mesh N; Save "mesh.msh";`.
  - `msh_parser.rs` — first-party parser for gmsh `.msh` v4.1
    (`$MeshFormat` + `$Nodes` + `$Elements`) that produces a
    canonical `valenx_mesh::Mesh` for `Line2`, `Tri3`, `Quad4`,
    `Tet4`, `Pyr5`, `Prism6`, `Hex8` plus their quadratic
    variants. Unknown blocks are tolerated; unsupported versions
    and reserved tag=0 usage error loudly.
- `run()` spawns gmsh with the same threaded stdout/stderr
  streaming and cancellation the OpenFOAM adapter uses, plus
  coarse progress (1D → 2D → 3D → write) derived from stdout
  tokens.
- `collect()` parses the emitted `.msh` into a canonical `Mesh`,
  writes `mesh.canonical.json` alongside it (serde-pretty),
  attaches both files + the generated `.geo` as `Artifact`s on
  the `Results` bundle.
- `tests/fixtures/minimal.valenx/cases/box-mesh/case.toml` — a
  unit-cube meshing fixture for reproducible round-trips.
- **`valenx-mesh::quality`** — new module with `signed_size` (length
  / area / signed volume, covering `Line2`, `Tri3`, `Quad4`,
  `Tet4`, `Hex8`, plus the quadratic → linear reduction path) and
  `aspect_ratio` (max edge / min edge for simplicial elements,
  body-diagonal ratio for `Hex8`). Rolls both into a
  `QualityReport` via `quality::report(mesh)` with total count,
  min/max/mean size, max aspect ratio, and an `inverted_count`
  flag for negatively-oriented elements. **154 tests** pass
  workspace-wide (up from 131 in the previous push).

**Phase 1 closed — Year-1 shell complete end to end:**
- **wgpu render pass.** eframe switched from glow to the wgpu
  backend; new `valenx-app::wgpu_renderer` builds an offscreen
  colour + `Depth32Float` depth target, flat-shaded triangle
  pipeline with back-face culling, `bytemuck`-safe `Vertex` + `Uniforms`
  types, and a growable vertex buffer. `WgpuRenderer::render()`
  returns an `egui::TextureId` that `viewport::show` paints into
  the viewport rect through `painter.image(id, rect, uv, …)`. The
  previous painter-only Lambert stays as a fallback for builds
  without wgpu. 4 unit tests (Vertex/Uniform sizes, MVP finite-ness,
  `triangles_to_vertices` 3-per-triangle invariant).
- **Tabbed bottom panel.** Bottom dock now flips between
  **Residuals** and **Log** (`valenx-app::log_panel`): 20k-line
  ring buffer, per-level filter switches (info / warn / error on by
  default, debug / trace off), autoscroll toggle, Clear button,
  colour per level. The log is fed from the same `RunEvent::LogLine`
  stream the residual parser already drains. 3 unit tests (ring-
  buffer cap, filter count, clear).
- **Settings window.** `valenx-app::settings` — Theme (Auto / Dark /
  Light), default shading mode, residual Y-axis scale (log10 /
  linear), "re-probe adapters on close" toggle, Reset-to-defaults.
  Settings apply the theme to the egui context on change and
  optionally kick off a registry re-probe when closed.
- **Adapter re-probe.** A `Re-probe` button in the browser's
  **Adapters** collapsible and a menu item in `Settings →
  Re-probe adapters` call `AdapterRegistry::probe_all()`, so a user
  who just installed OpenFOAM can pick it up without restarting.
  Also surfaced as two new command-palette entries.
- **OpenFOAM `collect()` is real.** `discover_artifacts()` walks
  the workdir, classifying `.vtk` / `.vtu` / `.vtp` / `.pvd`
  (VizData), `log.*` files (Log), `.csv` / `.tsv` / `.dat`
  (Tabular), `.png` / `.svg` (Image). Populates `Results.artifacts`
  deterministically. 2 unit tests (classification + recursive walk).
- **Residual chart honours settings.** `residuals::show` takes a
  `ResidualScale` and either log10-transforms or plots linear.
- **131 tests passing, zero clippy warnings, MSRV 1.88.**

**Phase 1 (earlier) — full end-to-end Year-1 shell:**
- **Shaded viewport.** `valenx-app::viewport` now has two render
  modes (`ShadingMode::Shaded` / `Wireframe`, toggleable from the
  View menu or command palette). Shaded draws flat-shaded filled
  triangles via `egui::Mesh` with a Lambert-on-a-dot-product light,
  back-face culling against the orbit camera's view vector, and a
  painter's-algorithm depth sort. Wireframe remains as a fallback.
  The shared `valenx-viz::projection` math means the future `wgpu`
  swap is a rasteriser change, not an API rewrite.
- **Live residual chart.** New `valenx-app::residuals` module +
  `egui_plot = 0.28` dep. Residuals are extracted from the solver's
  stdout via `valenx-adapter-openfoam::log_parser`, stored
  per-field on a log-scale axis, and plotted in a resizable bottom
  panel. 3 unit tests covering ingest + clear.
- **Run orchestration with a background thread.** New
  `valenx-app::run` module — `spawn()` moves `prepare → run` onto a
  `std::thread`, uses `mpsc::channel` to stream `RunEvent`s
  (Starting / Progress / LogLine / Finished / Failed) back to the
  UI. `ChannelProgressSink` + `ChannelLogSink` adapt the sink
  traits into channel sends. `CancellationToken` is cloned so the
  UI's Cancel button actually terminates the solver. 1 integration
  test against a dummy adapter.
- **Adapter registry wired into the browser tree.** `ValenxApp`
  owns an `AdapterRegistry`, probes on startup, and the left
  sidebar's **Adapters** collapsible shows every registered adapter
  with its display name + status dot + `Ready`/`Missing`/`Outdated`/
  `Broken`/`Disabled` label. OpenFOAM registers automatically; future
  adapters join by adding a line in `ValenxApp::init_registry`.
- **Command palette (Ctrl+P).** New `valenx-app::commands` module
  with 14 commands (File: Open project / Import STL; View: Front /
  Back / Right / Left / Top / Bottom / Iso / Frame / Toggle
  shading; Run: First case / Cancel; Help: About). Fuzzy
  subsequence matching, arrow-key navigation, Enter to invoke, Esc
  to close. Every command has a stable `CommandId` so future
  keybindings + analytics can reference them. 4 unit tests
  (matching, case-insensitivity, unique IDs, non-empty labels).
- **End-to-end CFD thread works.** Open `.valenx` project → Run →
  solver spawns on background thread → residual chart fills live →
  status bar shows final convergence. The only piece not working
  today is actually having OpenFOAM installed; the adapter's
  `Missing` status is honest when `simpleFoam` isn't on PATH.
- **MSRV + clippy cleanup.** Pre-existing `needless_range_loop`,
  `result_large_err` (×3), `field_reassign_with_default`, and
  `manual_contains` (×2) warnings across `valenx-core` and
  `valenx-fields` all fixed. Workspace is now `cargo clippy
  --workspace --all-targets`-clean with zero warnings.
- **121 tests passing** (up from 90 at the start of this push).

**Phase 1 (earlier) — live OpenFOAM solve + interactive viewport:**
- **`valenx-adapter-openfoam::run()` is live.** Spawns the prepared
  solver binary via `std::process::Command`, owns two reader threads
  for stdout/stderr, streams every line back through `RunContext`'s
  `LogSink`, and publishes progress (`0..100` based on the declared
  iteration budget) + per-field residuals (`ResidualSample`) as the
  solve progresses. Cooperative cancellation via `CancellationToken`
  kills the child before propagating `AdapterError::Cancelled`.
  Non-zero exits return `AdapterError::Run { exit_code, stderr,
  phase: RunPhase::Solve }` with the tail of stderr for diagnostics.
- **`valenx-adapter-openfoam::log_parser`** — self-contained parser
  for the subset of OpenFOAM solver stdout that matters:
  - `Time = N` iteration markers
  - `<solver>: Solving for <field>, Initial residual = X, Final
    residual = Y, No Iterations Z`
  - `ExecutionTime = T s ClockTime = …`
  - `LogSignal::Other` + warning capture for `--> FOAM Warning`
    and `--> FOAM FATAL` blocks
  - `intern_field` interns the known OpenFOAM field names into
    `&'static str` so `ResidualSample.field` stays zero-allocation
    on the hot path
  - 8 unit tests covering each line shape + unknown-field drop.
- **`last_residual_below`** convergence check: per-field, uses the
  most recent residual for every tracked field.

**Phase 1 (in progress) — rfd file dialogs + interactive viewport:**
- **`rfd = 0.14`** added to workspace deps; `valenx-app`'s File menu
  now opens real native file dialogs via `rfd::FileDialog::new()`
  for both "Open project…" (folder picker) and "Import STL…" (file
  picker filtered to `.stl`/`.STL`).
- **Interactive wireframe viewport** in `valenx-app::viewport`:
  - Drag with primary or middle mouse to **orbit** (0.5°/px).
  - Scroll to **zoom** (~1 %/notch).
  - Double-click to **frame** the mesh bounding box.
  - Number keys **1..7** snap to ViewCube directions (Front / Back /
    Right / Left / Top / Bottom / Iso); F reframes.
  - Triangles are projected through `valenx-viz::projection` and
    painted with a painter's-algorithm back-to-front depth sort.
  - Bounding-box wireframe overlay in amber.
  - Camera readout (azimuth / elevation / distance) in the
    bottom-left corner.
- **`valenx-viz::projection`** — new module with `project_point`,
  `project_triangle`, `ScreenPoint`. Shares the exact math the
  forthcoming `wgpu` pass will use, so the render-swap is a
  rasteriser change, not a math rewrite. 3 unit tests (origin maps
  near centre, behind-camera culls, zero-size viewport returns
  `None`).
- **`ValenxApp` state tests** — 4 tests covering default construction,
  STL load against the fixture (camera gets framed, error stays
  `None`), STL load against a missing path (sets `last_error`), and
  project load against a missing path.

**Phase 1 (in progress) — native window lands:**
- **`valenx-app` now opens a real native window** via `eframe` +
  `egui` + a `glow` rendering backend. Shell layout per
  [DESIGN.md § 6]:
  - Ribbon / menu bar (File / View / Help) with Exit wired to
    `egui::ViewportCommand::Close`.
  - Left browser panel with `Project` / `Cases` / `Geometry` /
    `Results` collapsibles that render live data from the loaded
    `.valenx` project and STL mesh.
  - Central viewport placeholder showing the loaded STL's format,
    triangle count, and axis-aligned bounding box; the wgpu render
    pass lands next.
  - Bottom status bar that surfaces the last error or a ready hint.
  - About dialog under Help.
- **Drag-and-drop loading.** Drop a `.stl` file onto the window to
  load it through `valenx-viz::stl::load`; drop a `.valenx`
  directory to load it through `valenx-core::LoadedProject::load`.
- **CLI arg loading.** `valenx some-model.stl` loads the file on
  first frame so screenshot tests can exercise the viewport
  reproducibly.
- **`eframe = 0.28` + `egui = 0.28`** added as workspace
  dependencies. The `glow` backend ships by default; the `wgpu`
  backend is a one-line feature swap when the viewport graduates.

**Phase 1 (in progress) — first physics thread lands:**
- **`valenx-adapter-openfoam::prepare()`** now emits a complete,
  runnable `simpleFoam` case tree from a canonical `CaseDef`:
  `system/controlDict`, `system/fvSchemes`, `system/fvSolution`,
  `constant/transportProperties`, `constant/turbulenceProperties`,
  `0/U`, `0/p`, plus turbulence fields (`k`, `omega` or `epsilon` or
  `nuTilda`, `nut`) for RAS models. Split across new modules:
  - `src/dict.rs` — OpenFOAM dict writer with `FoamFile` header,
    named-block nesting, indented-entry emission.
  - `src/case_input.rs` — canonical `case.toml` → typed
    `SimpleFoamInput` with `TurbulenceModel`, `SchemePreset`, `Fluid`,
    `Boundary` enums.
  - `src/simple_foam.rs` — dict generators for every file, with
    default turbulence estimates (`k = 1.5(I|U|)²`, `omega = k^0.5 /
    (C_μ^0.25 · L)`, `epsilon = C_μ^0.75 · k^1.5 / L`).
- `valenx-core::adapter_helpers` — shared helpers (`find_on_path`,
  `platform_suffix`, `stub_provenance`, `not_implemented`) so
  adapter crates don't each reimplement PATH lookup + provenance
  boilerplate.
- **`valenx-viz::stl`** — self-contained ASCII + binary STL loader
  (`TriangleMesh`, `StlTriangle`, `StlFormat`, `StlError`) with
  auto-detection heuristic, bounding-box computation, and 5 unit
  tests. No wgpu dep; callers can run pure-geometry tests.
- **`valenx-viz::camera`** — turntable `OrbitCamera` with `target`,
  `distance`, `azimuth_deg`, `elevation_deg`, `fov_y_deg` plus
  `ViewDirection` canonical angles for the ViewCube (Front / Back /
  Top / Bottom / Left / Right / Iso). `frame_bounds`, orbit
  clamping, RHS view/projection matrices via `nalgebra`.

**Phase 2-9 structural scaffolding — 17 new adapter crates:**
- `valenx-adapter-freecad` — CAD (LGPL-2.1 via subprocess).
- `valenx-adapter-occt` — OpenCASCADE BRep kernel (LGPL-2.1
  dynamic-linked).
- `valenx-adapter-gmsh` — meshing (GPL-2.0 via subprocess).
- `valenx-adapter-netgen` — curved-geometry meshing (LGPL-2.1).
- `valenx-adapter-calculix` — FEA (GPL-2.0 via subprocess).
- `valenx-adapter-code-aster` — industrial FEA (GPL-3.0).
- `valenx-adapter-elmer` — multi-physics FE (GPL-2.0).
- `valenx-adapter-openradioss` — explicit crash (AGPL-3.0).
- `valenx-adapter-cantera` — kinetics + thermochem (BSD-3-Clause,
  dynamic-linked).
- `valenx-adapter-lammps` — materials MD (GPL-2.0 via subprocess).
- `valenx-adapter-gromacs` — biomolecular MD (LGPL-2.1).
- `valenx-adapter-openems` — FDTD EM for antennas (GPL-3.0).
- `valenx-adapter-meep` — FDTD EM for photonics (GPL-2.0).
- `valenx-adapter-pybamm` — battery DFN/SPM/SPMe (BSD-3-Clause).
- `valenx-adapter-mujoco` — multibody + contact (Apache-2.0,
  dynamic-linked).
- `valenx-adapter-precice` — multi-physics coupling meta-adapter
  (LGPL-3.0).
- `valenx-adapter-su2` — compressible CFD + adjoint (LGPL-2.1).

  Every skeleton adapter: live `probe()` (PATH lookup for the tool's
  expected binaries), honest `NotImplemented` errors for
  `prepare/run/collect` scoped to the adapter's ROADMAP phase, full
  capability list + ribbon contribution names, full `AdapterInfo`
  (display name, version range, SPDX license, docs / homepage URLs,
  license mode).

**Per-phase acceptance docs under `docs/src/phases/`:**
- `README.md` — index of all 16 phases with status.
- `phase-01-foundation.md` through `phase-09-coupling.md` —
  goal / capability inventory / integrated tools / acceptance
  checklist / success metrics / next-phase pointer for the nine
  physics-domain phases.
- `phase-10-ux-polish.md` through `phase-16-stewardship.md` —
  briefer coverage of UX polish, HPC, optimisation, ML surrogates,
  plugin marketplace, enterprise deployment, and long-term
  governance.
- Linked from `docs/src/SUMMARY.md` so mdBook renders them.

**New RFCs:**
- **RFC 0007** ([rfcs/0007-coupling-adapters.md](./rfcs/0007-coupling-adapters.md))
  — the `CouplingAdapter` trait and participant model that Phase 9
  builds on.
- **RFC 0008** ([rfcs/0008-reports-and-export.md](./rfcs/0008-reports-and-export.md))
  — a declarative reports surface (live figures, PDF / HTML / OOXML
  export, default template gallery).

**Testing / verification:**
- Workspace now has **22 member crates**; `cargo check --workspace
  --all-targets` and `cargo test --workspace` both green.
- **Test count: 90 passing** (up from 47 before this session).
- End-to-end test: `valenx_adapter_openfoam::prepare()` against the
  `tests/fixtures/minimal.valenx/cases/cfd-steady` fixture emits the
  full dict tree (opportunistically skips the `native_command`
  assertion when OpenFOAM isn't on PATH, so CI still exercises the
  writer).

**Earlier this release (pre-session) — see git history:**
- **[ROADMAP.md](./ROADMAP.md)** — 20-year plan across 16 phases, with
  capability inventory, integrated tool registry (~75 OSS tools),
  governance evolution, and success metrics per year.
- **[ARCHITECTURE.md](./ARCHITECTURE.md)** — standalone overview of how
  the pieces fit together (layer cake, workspace, canonical types,
  registry, workflow DAG, coupling, license firewall, reproducibility).
- **[DESIGN.md](./DESIGN.md)** — design plan covering principles, the
  three-layer design system (tokens → components → patterns), screen
  inventory, viewport and data-viz mini-specs, error/recovery,
  documentation tone, novice-to-expert journey, iconography, motion,
  themes/accessibility/i18n, branding identity, settings, scripting
  UX, reports and export, app updates, multi-document model, tooling,
  process, contribution flow, bus-factor, telemetry stance,
  engineering prerequisites, RFC queue, per-phase timeline, success
  criteria, and validation.
- **[DESIGN_PRINCIPLES.md](./DESIGN_PRINCIPLES.md)** — the five
  non-negotiables in short form, referenced from every UI PR.
- **RFC 0004** ([rfcs/0004-results-and-fields.md](./rfcs/0004-results-and-fields.md))
  — canonical `Results`, `Field`, `ScalarRecord`, units, provenance,
  serialization to VTK / CGNS / HDF5 / JSON.
- **RFC 0005** ([rfcs/0005-design-principles.md](./rfcs/0005-design-principles.md))
  — formal adoption of the five design principles as binding.
- **RFC 0006** ([rfcs/0006-token-system.md](./rfcs/0006-token-system.md))
  — design-token schema, JSON source of truth, Rust-const generation
  pipeline, theme inheritance, export to mockup tools.
- **Rust workspace scaffold** — root `Cargo.toml` with 10 member
  crates (`valenx-app`, `valenx-core`, `valenx-geo`, `valenx-mesh`,
  `valenx-fields`, `valenx-viz`, `valenx-icons`, `valenx-fonts`,
  `valenx-design-tokens`, `valenx-adapter-openfoam`), workspace-wide
  shared metadata and dependency versions, release / dev / bench
  profiles, and `#![forbid(unsafe_code)]` on every crate by default.
- **`valenx-design-tokens` crate** — canonical `tokens.json`,
  JSON Schema, and a `build.rs` that generates typed Rust consts
  for color / spacing / typography / radius / border / z-index /
  motion tokens at compile time.
- **mdBook skeleton** — `docs/book.toml`, `docs/src/SUMMARY.md`,
  intro / contributing / changelog stubs.
- **`docs/design/`** — README, icon-inventory stub listing ~53
  Year-1 icon roles, `mockups/` placeholder, `patterns/`
  placeholder.
- **UI PR checklist** added to `.github/PULL_REQUEST_TEMPLATE.md`
  for PRs that touch user-facing UI (mockup link, snapshot test,
  tokens-only, keyboard reachable, localization key).
- **[TESTING.md](./TESTING.md)** — how to develop and test a native
  Rust desktop app (no browser, no localhost; `cargo run` opens a
  window).
- **[LANGUAGES.md](./LANGUAGES.md)** — language and crate choices
  (Rust primary, C for FFI only, no C++/JS/Go in tree).
- **[CONTRIBUTING.md](./CONTRIBUTING.md)** — contributor guide with
  dev setup, workflow, commit style, testing expectations.
- **[MAINTAINERS.md](./MAINTAINERS.md)** — maintainer list and
  promotion criteria.
- **[CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md)** — adopting the
  Contributor Covenant v2.1.
- **[SECURITY.md](./SECURITY.md)** — vulnerability disclosure policy,
  response SLAs, scope, coordinated disclosure process.
- **[POLICIES.md](./POLICIES.md)** — SemVer commitments, deprecation
  policy, LTS cadence, supported versions, MSRV policy, supported
  platforms, dependency licensing rules.
- **[rfcs/](./rfcs/)** — RFC process and template; initial RFCs:
  - [RFC 0001](./rfcs/0001-project-file-format.md) — `.valenx` project
    file format
  - [RFC 0002](./rfcs/0002-adapter-contract.md) — Adapter contract
  - [RFC 0003](./rfcs/0003-plugin-api.md) — WIT-based plugin API
  - `rejected/` — archive slot for rejected/withdrawn RFCs
- **`.github/workflows/ci.yml`** — CI skeleton: rustfmt, clippy,
  cross-platform test matrix, rustdoc, mdBook, `cargo-deny`,
  `cargo-audit`.
- **`deny.toml`** — license allow-list and supply-chain policy.
- **`rust-toolchain.toml`** — pinned Rust toolchain (stable, with
  rustfmt/clippy/rust-src).
- **`rustfmt.toml`** — formatter config, near-defaults.
- **`.editorconfig`** — editor consistency across contributors.
- **Issue and PR templates** — bug report, feature request, PR
  checklist.
- **`.github/FUNDING.yml`** — placeholder sponsor-link config
  (populated when accounts stand up).
- **`.github/dependabot.yml`** — weekly cargo + github-actions
  dependency update cadence.
- **`CITATION.cff`** — machine-readable academic citation metadata;
  added so Valenx is citable once a release ships.
- **First substantive canonical-types implementation** under Phase 1:
  - `valenx-fields`: Units (SI 7-tuple dimension algebra), TimeKey,
    Field, FieldKind, Location, RegionRef, ScalarRecord,
    Provenance, Artifact, FieldCatalog, ScalarCatalog, Results.
    Arithmetic + conversion unit tests ship with the crate.
  - `valenx-core`: AdapterError taxonomy, Physics / Capability
    enums, LicenseMode guard (including `assert_spawn_allowed`),
    the Adapter trait with AdapterInfo / VersionRange / ProbeReport
    / PreparedJob / RunReport / ResidualSample / RunContext /
    Capabilities / Case, and an AdapterRegistry with status
    classification (Ready / Missing / Outdated / Broken /
    Disabled) and status-count summary.
  - `valenx-geo`: BoundingBox, SourceFormat, Geometry, BRepHandle.
  - `valenx-mesh`: ElementType (line/tri/quad/tet/pyr/prism/hex
    including quadratic variants), ElementBlock, Region,
    BoundaryGroup, MeshStats, Mesh.
  - `valenx-adapter-openfoam`: Adapter impl with a real probe
    (PATH lookup for `simpleFoam` / `foamRun`, Windows `.exe`
    suffix), version range (v2306..v2506), capability list,
    ribbon contribution names, and scaffold `prepare` / `run` /
    `collect` returning honest "not implemented" errors plus
    empty-Results provenance for integration tests.
- **`tests/fixtures/minimal.valenx/`** — RFC-0001-conformant
  fixture for project-load tests: `project.toml`, `tools.lock`,
  `cases/cfd-steady/case.toml`, `geometry/box.stl`.
- **`valenx-core::project`** — full `.valenx` loader and writer
  implementing RFC 0001:
  - `Project`, `ProjectHeader`, `UnitsConfig`, `GeometrySection`,
    `GeometryEntry`, `MeshEntry`, `CasesSection`, `UiSection` —
    canonical manifest shape
  - `ToolsLock`, `ToolEntry`, `LockedIntegrationMode` — per-project
    tool pinning
  - `CaseDef`, `CaseHeader` — per-case `case.toml` with
    physics-specific sections stored as `toml::Value` for
    adapter-side validation
  - `LoadedProject::load()` — loader with path-safety checks
    (absolute paths + `..`-escapes rejected), format-major version
    gate, and honest structured errors via `ProjectLoadError`
  - `LoadedProject::save()` — atomic per-file writes (temp-file +
    rename) so a crash mid-save can't corrupt the manifest
  - Integration test under
    `crates/valenx-core/tests/project_roundtrip.rs` loading
    `tests/fixtures/minimal.valenx`, verifying manifest / tools
    lock / case fields, and round-tripping save→reload
- **`valenx-core::workflow`** — DAG orchestration types per
  ARCHITECTURE.md § 6:
  - `Workflow`, `WorkflowNode`, `WorkflowEdge`, `PortType` (the
    canonical edge types: Geometry, Mesh, Case, Results, Raw)
  - `Workflow::validate()` — checks unknown-node/port references
    and rejects type mismatches between ports
  - `Workflow::topo_order()` — Kahn's algorithm, cycle-detecting
  - 5 unit tests covering valid DAGs, unknown-node errors, port
    type mismatches, and cycle detection
- **`tests/README.md`** and **`installers/README.md`** — directory
  explainers so the empty dirs don't disappear from git.
- **[legacy-reference/](./legacy-reference/)** — preserved knowledge
  from the previous web-app iteration (case generators, subprocess
  bridges, status history) to inform Rust ports without being in the
  build.

### Changed
- **MSRV bumped 1.85 → 1.88** (`rust-toolchain.toml` +
  `Cargo.toml`'s `rust-version` + `POLICIES.md` pointer). Driven by
  `eframe 0.28`'s transitive `image@0.25.10` requirement. The
  "stable - 3" policy is unchanged; this bump tracks the floor as
  stable moves.
- **Clippy cleanup pass.** Workspace now runs `cargo clippy
  --workspace --all-targets` with zero warnings:
  - `valenx-adapter-openfoam::prepare()`: replaced `Vec::new()`-plus-
    `push()` with `vec![…]` literal; dropped an `u64 as u64` cast.
  - `valenx-fields::units`: `Mul` / `Div` impls use
    `dim.iter_mut().zip(…)` instead of `for i in 0..7`.
  - `valenx-core::project::loader::ProjectLoadError::Parse::source`
    now boxes `toml::de::Error` so `Result<_, ProjectLoadError>`
    stays small (fixes `clippy::result_large_err` in 3 places).
  - `valenx-core::registry`: `ready_for_physics` / `ready_for_capability`
    use slice `contains()` instead of `iter().any(|&x| x == …)`.
  - `valenx-core::registry` test: struct-update syntax instead of
    default + field reassign.
- Repository reset from the previous React / FastAPI / Tauri web-app
  codebase to a clean slate prepared for the native Rust rewrite.
- New [README.md](./README.md) oriented around the native-desktop
  direction.
- `.gitignore` rewritten for the Rust workspace (removing Node / venv /
  Tauri entries that no longer apply).
- `POLICIES.md` release cadence set to quarterly minor (previously
  claimed both "monthly" and "quarterly" in different sections).
- `SECURITY.md` support matrix consolidated into `POLICIES.md`;
  `SECURITY.md` now focuses on reporting and response mechanics.

### Removed
- `node_modules/`, `venv_311/`, `target/`, `dist/`, `build/`,
  `frontend/`, `backend/`, `src-tauri/`, and related artifacts from the
  web-app era.
- Orphan directories (`datasets/`, `installer/`, `release/`,
  `scripts/`, `training_output/`, `valenx-sdk/`, `valenx-website/`) —
  not deleted lightly; anything worth preserving was moved to
  `legacy-reference/`.

### Deprecated
- (nothing)

### Fixed
- (nothing — pre-alpha)

### Security
- (nothing — pre-alpha; security policy now in place for future reports)

---

## Versioning scheme

This project will follow SemVer 2.0.0 from 1.0.0 onward. Until then
the rules in [POLICIES.md](./POLICIES.md) apply:

- `0.MINOR.PATCH` during pre-alpha
- Breaking changes bump `MINOR` pre-1.0
- Each release notes go under a heading `## [X.Y.Z] - YYYY-MM-DD`
- Each change goes under `Added / Changed / Deprecated / Removed / Fixed / Security`

### Channels

Once releases begin:

- **Stable** — `X.Y.Z`
- **LTS** — `X.Y.Z-lts` (every ~18 months, supported 24 months)
- **Nightly** — `X.Y.Z-nightly.YYYYMMDD` (from `main`)

---

<!--
  When the first tag lands, this becomes:
    [Unreleased]: https://github.com/<org>/valenx/compare/v0.1.0...HEAD
    [0.1.0]:      https://github.com/<org>/valenx/releases/tag/v0.1.0
  Until then, the Unreleased link points at `master`.
-->
[Unreleased]: https://github.com/nochallenge/valenx/tree/master
