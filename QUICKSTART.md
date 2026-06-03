# Valenx — Quickstart

Five-minute walkthrough that runs a transient CFD case end-to-end
and shows the result in the viewport. Aimed at someone who's
cloned the repo for the first time and has Rust installed.

> **Pre-alpha.** No installers yet — you build from source.

## Prerequisites

| Tool | Why | How |
|---|---|---|
| **Rust** ≥ 1.88 | builds Valenx itself | https://rustup.rs |
| **OpenFOAM** ≥ v2306 (optional) | actually runs the CFD case; without it, you can still **prepare** the case and inspect the dict tree | https://www.openfoam.com/download/ |
| **gmsh** ≥ 4.12 (optional) | mesh generation; the bundled fixture skips this | https://gmsh.info |

You don't need any of the optional tools to verify the build —
the scoped `scripts/qa.sh` (or `scripts/qa.ps1` on Windows)
runs **10,000+ tests** that pass with no external dependencies,
and the workspace has zero clippy + zero rustdoc warnings.
A blanket `cargo test --workspace` is intentionally forbidden in this
repo — see `docs/QA.md` for the rationale.

For headless mesh-quality inspection, the project also ships:

- `cargo run --bin valenx-mesh-info -- mesh.json` — text + JSON
  output, with `--check max-skew=0.9 --check inverted=0` for
  CI gates.
- `cargo run --bin valenx-audit -- verify <log.jsonl>` and
  `cargo run --bin valenx-audit -- tail -n 50 <log.jsonl>` —
  audit-log integrity check + recent-activity tail. The tail
  command also supports `--since 2026-04-28T00:00:00Z` to drop
  entries before a cutoff timestamp.
- `cargo run --bin valenx-results -- <workdir>/results.json` —
  inspect the post-run sidecar (fields, scalars, artifacts,
  provenance) without firing up the GUI; `--format json` for the
  full pretty-printed envelope.
- `cargo run --bin valenx-report -- <workdir>/results.json --html report.html --markdown report.md --csv scalars.csv` —
  write a self-contained HTML report, a GitHub-flavoured Markdown
  summary, and/or a flat scalar history CSV from a finished run.
  CI-friendly: at least one of `--html` / `--markdown` / `--csv` is
  required, exit-code-driven.

For biology workflows (Phase 17):

- `cargo run --bin valenx-fasta -- inspect <file.fa>` — inspect /
  validate / extract sequences from FASTA files. Text + JSON output;
  `-` reads from stdin.
- `cargo run --bin valenx-pdb-info -- <file.pdb>` — structural
  summary (chain / residue / atom counts + element tally) from a
  PDB file.
- `cargo run --bin valenx-blast -- query <query.fa> <db>` — thin
  BLAST+ wrapper that auto-detects alphabet from the query and
  routes to blastp / blastn.

For sequence-alignment workflows (Phase 18):

- `cargo run --bin valenx-fastq -- inspect <file.fq>` — inspect /
  validate FASTQ files (4-line format). Text + JSON output;
  `-` reads from stdin.
- `cargo run --bin valenx-sam-info -- <file.sam>` — alignment
  summary (record count, mapped / unmapped tally, reference list,
  average MAPQ) from a SAM file.
- `cargo run --bin valenx-vcf-info -- <file.vcf>` — VCF summary
  (header-line count, sample count, total records, PASS / FAIL
  split, no-ALT count) from a VCF file (Phase 19); `-` reads
  from stdin.

For project scaffolding:

- `cargo run --bin valenx-init -- my-case --template <T>` where
  `<T>` is one of 136 names: `empty / cfd / fea / chemistry / su2 /
  openradioss / code-aster / netgen / meep / gromacs / gmsh /
  lammps / elmer-heat / biopython / rdkit / openmm / chimerax /
  oxdna / mdanalysis / colabfold / bwa / minimap2 / mafft /
  muscle / hmmer / samtools / esmfold / openfold / alphafold2
  (alias `af2`) / alphafold3 (alias `af3`) / bcftools / gatk
  (alias `hc`) / deepvariant (alias `dv`) / scanpy / scvi (alias
  `scvi-tools`) / nextflow (alias `nf`) / snakemake (alias `smk`)
  / pymol / vmd / igv (alias `igvtools`) / deepchem (alias `dc`)
  / openbabel (alias `obabel`) / avogadro (alias `avogadro2`) /
  psi4 / nwchem / xtb / rfdiffusion (alias `rfd`) / proteinmpnn
  (alias `mpnn`) / chroma / esm-if (aliases `esmif` /
  `inverse-folding`) / rfantibody (alias `rfab`) / esm3 / esmc
  (alias `esm-cambrian`) / bowtie2 (alias `bt2`) / mmseqs2 (alias
  `mmseqs`) / diamond (alias `dmnd`) / hisat2 (alias `hisat`) /
  star / salmon / kallisto / iqtree (alias `iqtree2`) / raxml-ng
  (alias `raxml`) / fasttree / viennarna (aliases `vienna` /
  `rnafold`) / rnastructure / nupack / copasi / bionetgen (alias
  `bng`) / physicell / vina (alias `autodock-vina`) / autodock4
  (alias `ad4`) / relion / eman2 (alias `eman`) / ctffind / art
  (alias `art-illumina`) / wgsim / badread / chopchop / crispor /
  cas-offinder (alias `cas-off`) / slim / msprime / tskit /
  rosetta / pyrosetta / beast2 (alias `beast`) / mrbayes (alias
  `mb`) / x3dna (alias `3dna`) / curves (alias `curves+`) / dssr
  / plumed / prody / cpptraj / pysbol (alias `sbol`) / j5 / cello /
  blast / clustalo / tcoffee / seurat / anndata / namd (aliases
  `namd2` / `namd3`) / sander (aliases `amber-sander` /
  `ambertools-sander`) / hoomd (aliases `hoomd-blue` /
  `hoomdblue`) / mdtraj / rosettafold (alias `rf`) / omegafold
  (alias `of`) / foldseek / smoldyn / mcell / pydna / jalview /
  fiji / cellprofiler / ilastik / planemo / cromwell / cwltool /
  molstar / ngl / dnachisel / lineardesign / icodon / mfold /
  eternafold / linearfold / be-designer / be-hive / primedesign /
  pegfinder / indelphi / forecast / alphamissense / crispritz /
  pksim / simrna`
  (run `--help` for full alias list).
- `cargo run --bin valenx-validate -- path/to/project.valenx` —
  pre-flight structural check on a project bundle (manifest,
  tools.lock, every case in `[cases].order`). Exits 0 on clean,
  1 on a structural issue. `--format json` for CI consumption.

## Build + run

```bash
git clone <this-repo> valenx
cd valenx
cargo build --workspace
cargo run -p valenx-app
```

First build takes a few minutes (the wgpu / eframe / nalgebra
stacks are not small). Subsequent builds are fast.

## Drive the workflow loop

Once the window opens:

1. **Open the bundled fixture project.**
   File → Open project → pick `tests/fixtures/minimal.valenx/`.
   Six cases appear in the browser, one per live adapter family:
   - `box-mesh` (gmsh — Delaunay tet mesh of a unit cube)
   - `cfd-steady` (simpleFoam — steady RANS over a box)
   - `cfd-transient` (pimpleFoam — transient PIMPLE)
   - `fea-cantilever` (CalculiX linear-static — tip-loaded beam)
   - `heat-cube` (Elmer steady heat conduction with two pinned faces)
   - `netgen-cylinder` (Netgen CSG — unit-radius cylinder)

2. **Check adapter status.**
   Each row has a coloured ● dot:
   - **green** = the adapter's tool is installed and ready to run
   - **gray** = tool not on PATH (you can still Prepare)
   - hover for the full status reason

3. **Prepare without running** (if you don't have OpenFOAM
   installed). Click `cfd-transient` → Run → "Prepare selected
   case (no execute)". The status bar shows the temp workdir.
   Click "Open in file browser" in the Results pane to see the
   generated `system/`, `constant/`, `0/` dict tree — exactly what
   `simpleFoam` would consume.

4. **Run the case** (if you have OpenFOAM). With `cfd-transient`
   selected, hit **F5**. The residual chart updates live; the
   status bar tracks `"Time = 0.0005"` etc. as the solver writes
   stdout.

5. **See the results.** When the run finishes:
   - The mesh auto-loads into the viewport.
   - The wireframe edges paint by the first scalar field
     (typically `p`) using a five-stop blue→red ramp.
   - The bottom-right colour-bar shows the field name + min/max.
   - The Results pane lists all fields the run produced.

6. **Switch fields.** Click `T` (or `Ux`, `k`, etc.) in the
   Results pane → wireframe re-paints with that field's range.

7. **Scrub through time** (transient runs only). The time-step
   slider in the Results pane lets you drag through every snapshot;
   the wireframe + legend update live.

8. **Edit and re-run.** Click "Open prepared workdir in file
   browser" → edit `system/controlDict` to change `endTime` or
   `deltaT` → close → "Run from prepared workdir". The solver
   runs against your edits, no re-prepare.

## Try a different case

The bundled fixture covers the six most-used adapters end-to-end.
To exercise an adapter that doesn't ship a fixture (e.g. SU2,
LAMMPS, GROMACS, Code_Aster), drop a `case.toml` like the one below
into `cases/<your-case>/` and the workflow loop above works
identically — adapter status badge, click-to-run, results pane.

```toml
# Minimal CalculiX dynamic-analysis case.toml. See the
# `fea-cantilever` fixture for the linear-static shape.

[case]
format  = "1.0"
name    = "cantilever"
physics = "fea"
solver  = "calculix.dynamic"
mesh    = "primary"

[structural]
analysis    = "linear-dynamic"
mesh_source = "mesh.canonical.json"

[structural.material]
name    = "steel"
E       = 210e9
nu      = 0.3
density = 7850.0

[[structural.boundaries]]
nset      = "fixed"
dof_start = 1
dof_end   = 3
value     = 0.0

[[structural.loads]]
nset  = "tip"
dof   = 2
force = -1000.0

[structural.step]
time_total     = 0.01
time_increment = 1e-4
output_fields  = ["U", "S"]
```

Drop that into a `cases/cantilever-dynamic/case.toml` inside any
`.valenx` project, generate the mesh with gmsh adapter (or hand-
write `mesh.canonical.json`), and the workflow loop above works
identically — adapter status badge, click-to-run, results pane.

The full case-schema reference for each adapter lives next to the
adapter source in `crates/valenx-adapters/<domain>/<adapter>/src/case_input.rs`.

## When something goes wrong

| Symptom | Likely cause | Fix |
|---|---|---|
| Adapter dot is gray | Tool not on PATH | Install the tool, then Settings → "Re-probe adapters" |
| "case requested X but found Y" on Run | Solver named in `case.toml` is different from the binary on PATH | Either install the right binary or edit `case.solver` |
| Wireframe paints blue but no colour variation | The selected field has constant values across the mesh (or only one node has data) | Pick a different field, or check the solver actually wrote varying output |
| Workdir filled but no fields in Results | Solver didn't produce `.vtu` (e.g. didn't run `foamToVTK`) | Run `foamToVTK -case <workdir>` manually, then re-`collect()` |
| App freezes during a run | Shouldn't happen — the run thread is fully isolated | File a bug with the log panel contents |

## What's next

- See `STATUS.md` for what works end-to-end today and what's
  scaffolded vs. planning-doc only.
- See `docs/src/phases/` for the per-phase backlog.
- See `CHANGELOG.md` for the running list of what shipped when.
