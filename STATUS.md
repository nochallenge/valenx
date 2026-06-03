# Valenx — implementation status

A snapshot of what works end-to-end today, what's scaffolded, and
what's planning-doc only. Updated alongside `CHANGELOG.md` as
features land.

> **Pre-alpha.** No shipping release yet. The project is on
> `master` only; the first tagged version (`0.1.0-alpha.1`)
> covers what this document describes.

> **🎯 Bio ecosystem complete + Phases 44.5 + 35.5 + 35.6 + 45
> (123 adapters across 43 phases).** Phase 44.5 (RNA folding
> expansion — mfold + EternaFold + LinearFold), Phase 35.5 (base
> + prime editing design — BE-Designer + BE-Hive + PrimeDesign +
> pegFinder), Phase 35.6 (edit-outcome prediction — inDelphi +
> FORECasT + AlphaMissense + CRISPRitz), and Phase 45
> (pharmacokinetics + RNA tertiary structure — PK-Sim + SimRNA)
> ship 13 more bio adapters on top of the bio-ecosystem-complete
> milestone reached at Phase 22.5 + 42 and the Phase 43 mRNA-design
> add-on, **opening the first PK/PD pharmacokinetics modeling
> domain in Valenx (PK-Sim) plus the first RNA tertiary 3D
> structure prediction domain in Valenx (SimRNA)** while
> sister-expanding the Phase 28 RNA secondary-structure trio
> (44.5) and the Phase 35 CRISPR design trio (35.5 + 35.6). The
> bio surface now spans alignment / base + prime editing /
> cheminformatics / CRISPR / cryo-EM / DNA geometry / docking /
> edit-outcome prediction / MD analysis / MD engines / microscopy
> / mRNA design / pharmacokinetics / phylogenetics / population
> genetics / protein design / quantum chemistry / RNA structure
> (2D + 3D) / sequence editors / sequence read simulators /
> single-cell / spatial stochastic / structure prediction /
> structure search / synthetic biology / systems biology / variant
> calling / viewers (desktop + web) / web visualization / workflow
> managers. With Phases 44.5 + 35.5 + 35.6 + 45 live, **Valenx now
> ships 141 live adapters** total spanning the physics-domain
> phases 1-9 plus the entire Phase 5.5 / 17 → 45 biology / biotech
> / chemistry expansion; the **bio-adapter count now totals 123
> across 43 biology / biotech / chemistry phases**. The headline
> number is *141 fully live* (the OCCT stub was rewritten as a
> pythonocc-core subprocess wrapper, dropping the last remaining
> stub).

## End-to-end workflows that work today

These are the user-visible flows you can drive in the app right now,
without writing Rust:

1. **Open a `.valenx` project** → cases appear in the browser, each
   tagged with the **adapter status** (green = Ready, gray = Missing,
   yellow = Outdated, red = Broken, dark gray = Unregistered) and a
   **run-history badge** (✓ / ✗ / ·).
2. **Click a case** → it becomes selected. Double-click runs.
3. **Run** the selected case (`F5`, Run menu, command palette,
   or "Run selected" button in the browser). Background thread
   spawns the solver; UI stays responsive.
4. **Live progress** — residual chart + log panel + status-bar
   percentage update as the solver writes stdout.
5. **Cancel** at any time via the Run menu or command palette.
6. **Run completes** → workdir captured, fields parsed, **mesh
   auto-loads** into the viewport.
7. **Field-coloured wireframe** paints the mesh edges by the first
   scalar OnNode field via a five-stop cool-to-warm divergent ramp.
8. **Colour-bar legend** (bottom-right) shows field name, min/max,
   and (for transient runs) the timestep.
9. **Click any other scalar** in the Results pane → wireframe re-
   paints with the new field's range.
10. **Time-series slider** (transient runs only) → scrub through
    every snapshot in the field's series; wireframe + legend update.
11. **Open the run / prepare workdir in the host file browser**
    (Explorer / Finder / xdg-open) from the Results pane buttons.

Plus an **expert workflow** for hand-tuned cases:

1. **Prepare** the selected case (no execute) → dict tree lands in
   `<temp>/valenx-prepare-<case>-<unix>/` without spawning the solver.
2. **Open the prepared workdir** → inspect / edit the generated
   files manually.
3. **Run from prepared workdir** → solver runs against your edits;
   the prepare step is skipped so your changes survive.

## Live adapters

| Adapter | Phase | Capability |
|---|---|---|
| **OpenFOAM** | 1 (CFD) | simpleFoam (steady incompressible RANS) + pimpleFoam (transient PIMPLE, laminar/RANS) + icoFoam (transient laminar PISO) + rhoSimpleFoam (steady compressible RANS) |
| **SU2** | 1 (CFD) | SU2_CFD batch on a user-provided .cfg + .su2 mesh, OMP threading |
| **gmsh** | 2 (CAD/mesh) | tet/hex/prism mesh generation from `.geo` scripts; canonical Mesh export |
| **Netgen** | 2 (mesh) | batch CSG / BREP meshing on `.geo` / `.geo2d` / `.step` / `.iges` / `.brep` |
| **FreeCAD** | 2 (CAD) | STEP/IGES → STL via Python script; primitives via the FreeCAD scripting API |
| **CalculiX** | 3 (FEA) | linear-static + linear-dynamic + modal + steady-thermal + transient-thermal `.inp` emission |
| **Elmer** | 3 (FEA/heat) | steady-state + transient (BDF-2) heat equation with optional initial conditions |
| **OpenRadioss** | 3 (FEA) | engine phase (`engine_<arch> -i <_0001.rad>`) on a user-prepared starter conversion |
| **Code_Aster** | 3 (FEA) | `as_run case.export` with companion `.comm` / `.med` / `.py` staged automatically |
| **Cantera** | 4 (chemistry) | equilibrium TP / HP / UV + 0-D batch reactor via Python script |
| **LAMMPS** | 5 (MD) | LJ + EAM potentials with NVE / NVT / NPT ensembles |
| **GROMACS** | 5 (MD) | `gmx mdrun` on a user-built `.tpr` (grompp left to the user) |
| **openEMS** | 6 (EM) | FDTD via Octave / MATLAB script |
| **Meep** | 6 (EM) | Python (default) or legacy Scheme `.ctl`, MPI-aware |
| **PyBaMM** | 7 (battery) | DFN / SPM / SPMe with CC / CCCV protocols |
| **MuJoCo** | 8 (multibody) | MJCF / URDF playback with constant-control policy and trajectory capture |
| **preCICE** | 9 (coupling) | meta-adapter staging participant configs + running `precice-tools check` |

### Biology (Phase 17)

| Adapter      | Capability |
|---|---|
| **Biopython**  | user-provided Python script subprocess |
| **RDKit**      | user-provided Python script subprocess |
| **OpenMM**     | Python-native MD via user script + DCD output |
| **ChimeraX**   | `.cxc` command scripts; renders to .png / .cxs |
| **oxDNA**      | DNA/RNA coarse-grained MD on `input.dat` |
| **MDAnalysis** | trajectory analysis Python script + DCD parse-on-collect |
| **ColabFold**  | protein structure prediction from FASTA |

### Sequence alignment (Phase 18)

| Adapter      | Capability |
|---|---|
| **BWA**        | short-read alignment (`bwa mem` / `bwa aln`) against a reference index |
| **minimap2**   | long-read + spliced + asm-vs-asm alignment with selectable preset |
| **MAFFT**      | multiple-sequence alignment (`--auto` by default; algorithm overrides exposed) |
| **MUSCLE**     | alternate MSA back-end (`muscle -align` / `-super5`; v3 fallback flags) |
| **HMMER**      | profile-HMM search (`hmmbuild` / `hmmsearch` / `phmmer` / `jackhmmer`) |
| **samtools**   | SAM/BAM utilities (`flagstat` / `view` / `sort` / `index` / `stats`) |

### Aligners expansion (Phase 18.5)

| Adapter      | Capability |
|---|---|
| **Bowtie2**    | Langmead & Salzberg's gapped FM-index short-read aligner (GPL-3.0); two-stage `bowtie2-build` → `bowtie2` pipeline (mirrors BWA `index` → `mem`); knobs `reference` / `reads` (1 or 2 entries — single / paired-end) / `threads` / `skip_index` / `preset` ∈ `{very-fast, fast, sensitive, very-sensitive}`; collects `out.sam` (`Tabular`, "Bowtie2 aligned reads") + `.log` files |
| **MMseqs2**    | Söding lab's many-vs-many protein search + clustering toolkit (MIT); per-action dispatch on `action ∈ {easy-search, easy-cluster, easy-linsearch}` via `build_command(...) -> Result<...>`; knobs `query` / `target` (required for search, ignored for cluster) / `output` / `sensitivity` ∈ `1.0..=7.5` (default 7.5) / `threads`; collects `output` as `Tabular` with per-action label; binary is `mmseqs` (no `2` suffix) |
| **DIAMOND**    | Buchfink, Reuter & Drost's ultra-fast BLAST-protocol-compatible protein aligner (GPL-3.0); per-action dispatch on `action ∈ {blastp, blastx, makedb}`; knobs `query` / `database` / `output` / `sensitivity` ∈ `{default, fast, sensitive, more-sensitive, very-sensitive, ultra-sensitive}` (`--default` flag omitted when `default` since DIAMOND has no such flag) / `threads`; in `makedb` the field roles flip (input FASTA → `query`, output DB basename → `database`); collects `output` (`Tabular`, "DIAMOND <action> hits") for search modes, `<database>.dmnd` (`Native`, "DIAMOND .dmnd database") for `makedb` |

### RNA-seq alignment (Phase 18.6)

| Adapter      | Capability |
|---|---|
| **HISAT2**     | Daehwan Kim's graph-based splice-aware RNA-seq aligner (GPL-3.0); two-stage `hisat2-build` → `hisat2` pipeline (mirrors BWA `index` → `mem`); knobs `reference` (FASTA, required) / `reads` (1 or 2 entries — single / paired-end) / `threads` / `skip_index` / `strandness` ∈ `{unstranded, F, R, FR, RF}` (`--rna-strandness` flag omitted when `unstranded`); collects `out.sam` (`Tabular`, "HISAT2 aligned reads") |
| **STAR**       | Alex Dobin's spliced RNA-seq aligner (MIT); capitalized binary name (`find_on_path(&["STAR"])`); two-stage `--runMode genomeGenerate` → `--runMode alignReads` pipeline; knobs `genome_dir` (required) / `reference` (required when generating index) / `reads` (1 or 2 entries) / `threads` / `skip_index` / `output_type` ∈ `{BAM_Unsorted, BAM_SortedByCoordinate, SAM}` (underscore canonical names map to STAR's two-arg `--outSAMtype` form, e.g. `BAM_SortedByCoordinate` → `--outSAMtype BAM SortedByCoordinate`) / `sjdb_gtf` (optional GTF for splice-junction-aware indexing); collects `star_Aligned.out.{bam,sam}` (`Tabular` / `Native`, "STAR aligned reads") + `star_Log.final.out` (`Log`, "STAR alignment summary") |

### Alignment toolkit expansion (Phase 18.7)

Sister-adapter expansion of the Phase 18 / 18.5 / 18.6 alignment
beachhead — BLAST+ / Clustal Omega / T-Coffee round out the
foundational sequence-alignment surface that the existing BWA /
minimap2 / MAFFT / MUSCLE / HMMER / samtools / Bowtie2 / MMseqs2
/ DIAMOND / HISAT2 / STAR set deferred, covering canonical
sequence-database search (BLAST+, the seminal NCBI tool every
wet-lab pipeline reaches for first), modern HMM-driven
progressive multiple-sequence alignment (Clustal Omega), and
consistency-weighted library-based multiple-sequence alignment
(T-Coffee, the canonical choice for distantly-related sequences).

| Adapter      | Capability |
|---|---|
| **BLAST+**     | NCBI's seminal sequence-database search tool (Public Domain — US government work); ships five user-facing search programs (`blastn` / `blastp` / `blastx` / `tblastn` / `tblastx`) covering every nucleotide / protein search direction; single-binary subprocess shape (sister to Phase 18 BWA) with per-program CLI `<program> -query <query> -db <database> -out blast_results.txt -evalue <evalue> -outfmt <outfmt> -num_threads <threads> [extras...]`; knobs `program` (one of the five — required) / `query` (FASTA query file; required) / `database` (BLAST database path prefix — the user supplies the path stem and BLAST resolves the `<prefix>.nhr` / `<prefix>.phr` / etc. sidecars itself; required) / `evalue` (`f64`, default 10.0 — BLAST's own default) / `outfmt` (`u8`, default 0 — pairwise text; 6 = tabular, 7 = tabular-with-comments, 11 = BLAST archive) / `threads` (`usize`, default 1) / `extra_args`; `prepare()` resolves `query` against the case directory when relative, validates the database parent directory exists on disk (the database files themselves use the prefix convention so we can't validate them by name), looks up the per-program binary via `find_on_path(&[&input.program])`, and pins the output filename to `blast_results.txt`; collects `blast_results.txt` (`Tabular`, "BLAST search results") + `*.log` (`Log`); probe via `find_on_path(&["blastn", "blastp"])` — at least one BLAST+ binary on PATH counts as installed; version range `2.10.0..3.0.0` (BLAST+ 2.10 (2019) is the modern stable line that ships every NCBI release; upper bound 3.0 reserves room for an eventual major bump); `bio.blast.search` ribbon capability |
| **Clustal Omega** | Sievers / Higgins' modern HMM-driven progressive multiple-sequence aligner (GPL-2.0); modern successor to ClustalW that scales to thousands of sequences; single-binary subprocess shape (sister to Phase 18 MAFFT) with `clustalo -i <input> -o <basename>.<ext> --outfmt=<outfmt> --threads=<N> [extras...]`; knobs `input` (FASTA multi-sequence input; required) / `output_basename` (filename stem; required) / `outfmt` (default `"clustal"`; whitelist `clustal` / `fasta` / `phylip` / `vienna` / `nexus`) / `threads` (`usize`, default 1) / `extra_args`; `prepare()` derives `<ext>` from `outfmt` (`clustal` → `.aln`, `fasta` → `.fasta`, `phylip` → `.phy`, `vienna` → `.vie`, `nexus` → `.nex`, default `.aln`); collects `<output_basename>*` (`Tabular`, "Clustal Omega alignment") + `*.log` (`Log`); probe via `find_on_path(&["clustalo"])`; version range `1.2.0..2.0.0` (Clustal Omega 1.2 is the modern stable release line); `bio.clustalo.align` ribbon capability |
| **T-Coffee**   | Notredame / Higgins' library-based consistency-weighted multiple-sequence aligner (GPL-2.0); combines pairwise alignments from many sources into a single MSA — the canonical choice for difficult distantly-related sequences where progressive aligners lose accuracy; supports specialised modes (`expresso` / `psicoffee` / `mcoffee`) selectable through the `mode` knob; single-binary subprocess shape (sister to Clustal Omega) with `t_coffee <input> -output=<outfmt> -outfile=<basename>.aln [-mode=<mode>] [extras...]` (note T-Coffee's `=`-style flag form); knobs `input` / `output_basename` / `outfmt` (default `"clustalw"`; T-Coffee's own naming follows ClustalW conventions: `clustalw` / `fasta_aln` / `phylip` / `msf`) / `mode` (`Option<String>` — omit for default progressive mode) / `extra_args`; output always pinned to `.aln`; collects `<output_basename>*` (`Tabular`, "T-Coffee alignment") + `*.dnd` (`Native`, "T-Coffee guide tree" — Newick guide tree consumed by downstream phylogenetics) + `*.log` (`Log`); probe via `find_on_path(&["t_coffee"])` (note the underscore — T-Coffee installs as `t_coffee`, not `t-coffee`); version range `13.0.0..14.0.0` (T-Coffee 13.x is the modern stable line); `bio.tcoffee.align` ribbon capability |

### Transcript quantification (Phase 20)

| Adapter      | Capability |
|---|---|
| **Salmon**     | Rob Patro's quasi-mapping plus two-phase EM transcript-level quantification tool (GPL-3.0); two-stage `salmon index` → `salmon quant` pipeline (mirrors BWA `index` → `mem`); knobs `transcriptome` (FASTA, required) / `index_dir` (required) / `reads` (1 or 2 entries — single / paired-end) / `output_dir` / `threads` / `skip_index` / `libtype` (default `"A"` — Salmon's libtype DSL: `"A"` auto, `"U"` unstranded, `"ISF"` / `"ISR"` paired-end stranded forward / reverse, `"IU"` paired-end unstranded — non-whitelisted because the libtype DSL has many valid combos); collects `quant.sf` (`Tabular`, "Salmon transcript quantification") + `cmd_info.json` (`Log`, "Salmon command info") |
| **Kallisto**   | Lior Pachter's pseudoalignment-based transcript quantifier (BSD-2-Clause); two-stage `kallisto index` → `kallisto quant` pipeline; index is a single `.idx` file (kallisto convention — not a directory like Salmon / STAR); knobs `transcriptome` (FASTA, required) / `index` (single `.idx` file) / `reads` (1 or 2 entries) / `output_dir` / `threads` / `skip_index` / `fragment_length` + `fragment_sd` (`f64`, both required when `reads.len() == 1` — kallisto auto-detects fragment statistics from paired-end reads but cannot for single-end); collects `abundance.tsv` (`Tabular`, "Kallisto transcript abundance") + `abundance.h5` (`Native`, "Kallisto HDF5 abundance") + `run_info.json` (`Log`, "Kallisto run info") |

### Structure prediction (Phase 17.5)

| Adapter         | Capability |
|---|---|
| **ESMFold**     | Meta protein language model — single-sequence structure prediction (no MSA / database) |
| **OpenFold**    | PyTorch AF2 reimplementation; full `model_*` preset family validated at case-input layer |
| **AlphaFold 2** | DeepMind reference AF2 (`run_alphafold.py`) — `monomer` / `monomer_ptm` / `multimer` presets |
| **AlphaFold 3** | DeepMind all-atom complex predictor (JSON job spec); non-commercial weights flagged via probe warning |

### Structure tools expansion (Phase 17.7)

Sister-adapter expansion of the Phase 17.5 structure-prediction
beachhead and the Phase 17 ColabFold adapter — RoseTTAFold +
OmegaFold + FoldSeek round out the protein structure prediction +
structure search surface that Phase 17.5 ESMFold / OpenFold /
AlphaFold 2 / AlphaFold 3 + Phase 17 ColabFold opened, covering
the Baker-lab original 3-track network (RoseTTAFold), the MSA-
free larger-backbone HelixonAI predictor (OmegaFold), and the
Steinegger-lab 3D analogue of MMseqs2 (FoldSeek).

| Adapter      | Capability |
|---|---|
| **RoseTTAFold** | Baker lab's original 3-track structure-prediction network (MIT) — three concurrent attention tracks over the 1D sequence, the 2D pair-distance map, and the 3D Cartesian backbone with cross-track message passing refining all three jointly; canonical pre-AlphaFold-3 sibling that established the 3-track SE(3)-equivariant attention pattern; Python-script subprocess shape (sister to Phase 17.5 ESMFold / Phase 17 Biopython); knobs `script` (path to user-supplied Python script; required, `.py` enforced) / `python` (interpreter name; default `"python3"`) / `fasta` (input FASTA query sequence; required) / `output_basename`; `prepare()` enforces the `.py` extension, resolves `script` and `fasta` against the case directory when relative, stages both into the workdir under their original filenames, then writes a flat `valenx_params.json` containing `output_basename` and the bare `fasta` filename, and builds `<python> <staged_script>`; collects `<output_basename>*.pdb` (`Native`, "RoseTTAFold predicted structure" — pLDDT-style per-residue confidence in the B-factor column), `<output_basename>*.npz` (`Native`, "RoseTTAFold confidence arrays"), `*.log` (`Log`); probe via `find_on_path(&["python3", "python"])` — deliberately doesn't try `import rosettafold` (RoseTTAFold is not a pip package); pushes a probe warning whenever Python is detected: "RoseTTAFold model weights + dependencies not bundled — clone https://github.com/RosettaCommons/RoseTTAFold and follow the install README"; version range `1.0.0..3.0.0` (RoseTTAFold 1.x is the original 2021 release; RoseTTAFold 2 / RoseTTAFold All-Atom is the late-2023 follow-up that adds nucleic-acid + small-molecule support) |
| **OmegaFold** | HelixonAI's single-sequence protein-structure predictor (Apache-2.0); MSA-free like ESMFold but uses a larger pre-trained transformer backbone trained on a much wider sequence corpus — works on single sequences (sister to ESMFold from Phase 17.5) but routinely matches AlphaFold-2-with-MSA quality on orphan / synthetic / fast-evolving sequences where MSA-based methods struggle; ships its own CLI binary (`omegafold <fasta> <output_dir>`) and falls back to `<python> -m omegafold ...` when the CLI launcher isn't on PATH but Python is; knobs `fasta` (input FASTA query sequence; required) / `output_basename` (workdir-relative output directory name; required, non-empty) / `python` (interpreter name; default `"python3"`; used only as fallback) / `model_dir` (`Option<PathBuf>` — optional pre-downloaded model checkpoint directory; OmegaFold defaults to `~/.cache/omegafold_ckpt` when omitted); `prepare()` builds `omegafold <fasta> <output_basename> [--model <model_dir>]` with the FASTA passed by absolute path (NOT staged into the workdir); collects walks one level deep into the `<output_basename>/` subdirectory for `*.pdb` (`Native`, "OmegaFold predicted structure" — per-residue confidence in the B-factor column) and `*.json` (`Log`, "OmegaFold metadata"), plus the workdir-top-level `*.log` (`Log`); probe via `find_on_path(&["omegafold", "python3", "python"])` — surfaces a warning if `omegafold` itself isn't on PATH but Python is ("OmegaFold CLI not found on PATH; install via pip install git+https://github.com/HeliXonProtein/OmegaFold.git"); no academic-license caveat; version range `1.0.0..2.0.0` |
| **FoldSeek** | Steinegger lab's protein-structure search via the 3Di alphabet (GPL-3.0); the **3D analogue of Phase 18.5 MMseqs2** — both from the Steinegger lab, both built on the same fast many-vs-many search core, but FoldSeek encodes the per-residue 3D geometry as a custom "3Di alphabet" (a 20-letter alphabet over local backbone geometry patterns, designed so structural matches have high 3Di alphabet identity) and runs 3Di-vs-3Di comparisons at sequence-search speed; finds structural homologs at PDB-scale in seconds rather than the hours / days HMM-based or geometry-based search tools take; single-binary subprocess shape (sister to Phase 18.5 MMseqs2 / Phase 18 BWA) with `foldseek easy-search <query> <database> <basename>.m8 tmp_<basename> --threads <N> [extras...]` (`tmp_<basename>` is a per-run temp directory FoldSeek requires); knobs `query` (`.pdb` / `.cif` query structure; required) / `database` (FoldSeek database path prefix — the user supplies the path stem and FoldSeek resolves the `<prefix>_*` sidecar files itself; required) / `output_basename` / `threads` (`u32`, default 1) / `extra_args`; `prepare()` resolves `query` against the case directory when relative, validates the `database` parent directory exists on disk (the database files themselves use the prefix convention so we cannot validate them by name — same shape as Phase 18.7 BLAST+'s `database` validation), composes the invocation with the per-run temp directory pinned to `tmp_<basename>`; collects `<output_basename>.m8` (`Tabular`, "FoldSeek search results" — the canonical BLAST-style M8 hit table format), `*.log` (`Log`); the temp directory is not surfaced in artifacts; probe via `find_on_path(&["foldseek"])`; version range `8.0.0..10.0.0` |

### Variant calling (Phase 19)

| Adapter         | Capability |
|---|---|
| **bcftools**    | VCF/BCF multitool (`view` / `call` / `filter` / `concat`) with per-action dispatch |
| **GATK**        | Broad Institute reference variant caller — `HaplotypeCaller` on a sorted BAM with optional intervals (BED) |
| **DeepVariant** | Google ML-based variant caller — `run_deepvariant` with typed `model_type` ∈ `{WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA}` |

### Single-cell genomics (Phase 19.5)

| Adapter      | Capability |
|---|---|
| **Scanpy**   | de-facto Python single-cell analysis library (BSD-3-Clause); user-provided Python script imports `scanpy as sc` and reads `valenx_params.json` knobs (`input_h5ad`, `output_h5ad`, `n_top_genes`, `n_pcs`, `n_neighbors`, `resolution`); collects `.h5ad` ("Scanpy AnnData output"), `.png` / `.pdf` ("Scanpy plot"), `.csv` / `.tsv` (`Tabular`, "Scanpy table") |
| **scVI**     | probabilistic deep-learning models for single-cell data via `scvi-tools` (BSD-3-Clause); same Python-script subprocess shape as Scanpy with `valenx_params.json` knobs (`input_h5ad`, `output_h5ad`, typed `model` ∈ `{scvi, scanvi, totalvi, linear-scvi}`, `n_latent`, `n_layers`, `max_epochs`, optional `batch_key`); collects `.h5ad` ("scVI AnnData output"), `.png` / `.pdf` ("scVI plot"), `.csv` / `.tsv` (`Tabular`, "scVI table") |

### Single-cell expansion (Phase 19.6)

Sister-adapter expansion of the Phase 19.5 single-cell Scanpy +
scVI beachhead — Seurat / AnnData round out the single-cell
genomics surface that Phase 19.5 explicitly deferred, covering
the dominant R-based single-cell library (Seurat, the Satija
lab's reference toolkit that drives roughly half of single-cell
papers worldwide) and the canonical Python data-container
library (AnnData, the scverse foundation library that scanpy /
scvi / scirpy / squidpy / muon all read and write). Phase 19.6
introduces the **Rscript subprocess pattern** to Valenx — the R
analogue of the Python-script pattern that Phase 17 Biopython /
RDKit / OpenMM, Phase 19.5 Scanpy / scVI, and Phase 33 pySBOL
established — so future R-based bioinformatics tools
(Bioconductor, edgeR, DESeq2, monocle3, scran) can ride the
same infrastructure without further plumbing.

| Adapter      | Capability |
|---|---|
| **Seurat**     | Satija lab's R-based single-cell analysis toolkit (MIT) — dominant single-cell analysis library on the R side of bioinformatics (clustering, dimensionality reduction, integration, marker discovery, spatial transcriptomics, plus the broader Satija lab tooling Azimuth / signac / BPCells); introduces the **Rscript subprocess pattern** to Valenx; user supplies an `.R` script referenced from `[bio.seurat].script` in `case.toml` that loads `library(Seurat)` and reads `valenx_params.json` for the parsed knobs via `jsonlite::fromJSON`; knobs `script` (path to user-supplied R script; required, `.R` enforced) / `rscript` (binary name; default `"Rscript"`) / `input_data` (`Option<PathBuf>` — optional input matrix; supports `.h5` / `.mtx` / `.rds` so users can drop in 10x HDF5, sparse Matrix Market, or pre-saved Seurat object formats) / `output_basename`; `prepare()` enforces the `.R` extension, stages script + optional input_data into the workdir, writes `valenx_params.json` with `output_basename` and `input_data` (staged filename when set; key omitted entirely when `None` rather than emitted as `null`, matching the hand-rolled JSON convention the rest of the bio adapters use); collects `<output_basename>*.rds` (`Native`, "Seurat object (RDS)" — canonical R-serialised Seurat object format consumed by every downstream Seurat / signac / Azimuth pipeline) / `<output_basename>*.csv` (`Tabular`, "Seurat output table") / `<output_basename>*.png` (`Native`, "Seurat plot") / `*.log` (`Log`); probe via `find_on_path(&["Rscript"])` — surfaces a warning when Rscript is missing rather than failing (same shape as Phase 17 Biopython probe); the probe deliberately does **not** attempt to confirm Seurat itself is installed because that would require running R (an expensive multi-second startup at probe time that conflicts with the rest of the registry's snappy PATH-lookup probes); version range `4.0.0..6.0.0` (Seurat 4.x is the modern stable line; Seurat 5 (2024) is the current major; upper bound 6.0 reserves room for an eventual major bump); `bio.seurat.analyze` ribbon capability |
| **AnnData**    | scverse's Python single-cell HDF5 data container library (BSD-3-Clause) — canonical container that ties the entire scverse Python ecosystem together (scanpy / scvi / scirpy / squidpy / muon all read and write `.h5ad`); the standalone AnnData adapter lets users preprocess, convert, or inspect `.h5ad` files independent of the analysis pipelines that consume them; Python-script subprocess shape (sister to Phase 19.5 Scanpy / scVI); user supplies a Python script referenced from `[bio.anndata].script` in `case.toml` that imports `anndata` and reads `valenx_params.json` for the parsed knobs; knobs `script` (path to user-supplied Python script; required, `.py` enforced) / `python` (interpreter name; default `"python3"`) / `input_h5ad` (`Option<PathBuf>` — optional input single-cell file; supports `.h5ad` (canonical AnnData) and `.h5` (10x HDF5)) / `output_basename`; `prepare()` enforces the `.py` extension, stages script + optional input_h5ad into the workdir, writes `valenx_params.json` with the same hand-rolled shape as Seurat (`output_basename` plus `input_h5ad` (staged filename when set, key omitted when `None`)); collects `<output_basename>*.h5ad` (`Native`, "AnnData h5ad file") / `<output_basename>*.csv` (`Tabular`, "AnnData output table") / `<output_basename>*.png` (`Native`, "AnnData plot") / `*.log` (`Log`); probe via `find_on_path(&["python3", "python"])` then `<python> -c "import anndata"` — on import failure surface as a `ProbeReport.warnings` entry (not error) so non-standard installs aren't blocked; version range `0.9.0..1.0.0` (AnnData 0.9 is the modern stable line that pairs with scanpy 1.9 / scvi-tools 1.x; upper bound 1.0 reserves room for the eventual 1.0 stabilisation); `bio.anndata.process` ribbon capability |

### Molecular viewers (Phase 23)

| Adapter      | Capability |
|---|---|
| **PyMOL**    | open-source PyMOL `.pml` script-driven renderer (`pymol -c -q <script>`); collects `.png` / `.pse` / `.cif` / `.pdb` outputs |
| **VMD**      | Tcl-scripted MD trajectory viewer (`vmd -dispdev text -e <script>`); academic-license probe warning surfaces in `ProbeReport.warnings` |
| **IGV**      | `igvtools` wrapper for headless BAM/VCF indexing — per-action dispatch on `action ∈ {index, count, sort, tile}` (companion GUI viewer is out of scope) |

### Protein design (Phase 27)

| Adapter          | Capability |
|---|---|
| **RFdiffusion**  | GPU-driven protein backbone generation; user-provided Python script + `valenx_params.json` knobs (`mode` / `num_designs` / `diffusion_steps`); collects `<basename>_*.pdb` typed via `valenx_bio::format::pdb` |
| **ProteinMPNN**  | sequence design from backbone PDB; user-provided Python script + `valenx_params.json` knobs (`model_variant` ∈ `{vanilla, soluble, ca-only}` / `temperature` / `num_seq_per_target`); collects `<basename>.fa` parsed via `valenx_bio::format::fasta` |

### Protein design expansion (Phase 27.5)

| Adapter         | Capability |
|---|---|
| **Chroma**      | Generate Biomedicines' joint backbone + sequence diffusion model (Apache-2.0); user-provided Python script + `valenx_params.json` knobs (`num_samples` / `length` / `temperature`); collects `<basename>*.pdb` ("Chroma design") and `<basename>*.fa` (`Tabular`, "Chroma sequence") |
| **ESM-IF**      | Meta's GVP-based inverse-folding sequence designer via the `fair-esm` package (MIT) — sister to ProteinMPNN with a different model; user-provided Python script + `valenx_params.json` knobs (`input_pdb` / `model` / `temperature` / `num_samples`); collects `<basename>.fa` parsed via `valenx_bio::format::fasta` for a richer `"ESM-IF · N sequences"` label |
| **RFantibody**  | RosettaCommons antibody-specific fork of RFdiffusion (BSD-3-Clause); user-provided Python script + `valenx_params.json` knobs (`framework_pdb` / `target_pdb` / `design_loops` ⊆ `{H1, H2, H3, L1, L2, L3}` / `num_designs` / `diffusion_steps`); collects `<basename>*.pdb` ("RFantibody design") |

### EvolutionaryScale models (Phase 27.6)

Closes out the open-source EvolutionaryScale ESM lineup at 4 of 4
tools — ESMFold (Phase 17.5, structure prediction) + ESM-IF
(Phase 27.5, inverse-folding sequence design) + ESM3 (Phase 27.6,
generative multi-modal joint reasoning) + ESMC (Phase 27.6, protein
representation embeddings). All four ride the same EvolutionaryScale
`esm` Python package under the hood — installing one installs them
all.

| Adapter      | Capability |
|---|---|
| **ESM3**     | EvolutionaryScale's flagship generative multi-modal protein model (Cambrian-Open-License); Python-script subprocess shape — user's script imports `esm` and reads `valenx_params.json` knobs (`model_variant` ∈ `{open, open-multimer, small}` / `mode` ∈ `{design, inverse-fold, scaffold, predict}` / `num_samples` / `input_pdb` (required for `inverse-fold` / `scaffold`) / `input_fasta` (required for `predict`) / `temperature` / `output_basename`); collects `<basename>*.pdb` (`Native`, "ESM3 generated structure") and `<basename>*.fa` (`Tabular`, "ESM3 generated sequence"); probe via `find_on_path(&["python3", "python"])` then `python -c "import esm"` — surfaces an install hint when Python is on PATH but `esm` isn't importable |
| **ESM Cambrian** (ESMC) | EvolutionaryScale's open-weight protein representation model (Cambrian-Open-License); same Python-script subprocess shape as ESM3; `valenx_params.json` knobs (`input_fasta` / `model_variant` ∈ `{esmc-300m, esmc-600m}` / `pooling` ∈ `{per-residue, mean}` / `output_basename`); collects `<basename>.{npy,npz,parquet}` (`Tabular`, "ESMC embeddings"); probe identical to ESM3 (`import esm`); the init alias `esm-cambrian` resolves to the same template as `esmc` |

### RNA structure (Phase 28)

| Adapter           | Capability |
|---|---|
| **ViennaRNA**      | most-cited RNA secondary-structure suite (custom non-commercial / academic-use license); single-binary `RNAfold -i <input> -T <temperature> [-p] [--noGU if !allow_gu] [extras...]` shape (writes dot-bracket structure to stdout, captured to `output` via the MAFFT-style stdout-redirect pattern); knobs `input` (FASTA) / `output` / `temperature` (Celsius, default 37.0) / `partition_function` (toggles `-p` for partition function + base-pair probabilities) / `allow_gu` (default true; `--noGU` disables GU pairs); collects `output` as `Native` artifact "ViennaRNA secondary structure"; probe via `find_on_path(&["RNAfold"])` (capital R-N-A); **academic-license-only** — probe pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` |
| **RNAstructure**   | Mathews lab's classic RNA folding toolkit (BSD-3-Clause); single-binary `Fold <input> <output> -m <max_structures> -p <max_percent> -t <temperature> [extras...]` shape (binary literally named `Fold`); knobs `input` (FASTA or `.seq`) / `output` (`.ct` connection-table) / `max_structures` (default 20, ≥ 1) / `max_percent` (% of MFE, default 10, in `0..=100`) / `temperature` (Kelvin — RNAstructure's convention, default 310.15, > 0.0); collects `output` as `Native` artifact "RNAstructure connectivity table"; probe via `find_on_path(&["Fold"])` |
| **NUPACK**         | Caltech's nucleic-acid package (custom academic-only license); Python-script subprocess shape (NUPACK 4 is Python-driven; the 3.x CLI is deprecated); user-provided Python script imports `nupack` and reads `valenx_params.json` knobs (`input` / `output_basename` / `temperature` (Celsius, default 37.0) / `sodium` (molar, default 1.0)); collects `<output_basename>*` (`Native`, "NUPACK output") and `.npc` / `.json` files (`Tabular` / `Log`); probe via `find_on_path(&["python3", "python"])` then `python -c "import nupack"`; **academic-license-only** — probe pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` |

### RNA folding expansion (Phase 44.5)

Sister-adapter expansion of the Phase 28 RNA structure trio
(ViennaRNA / RNAstructure / NUPACK) — mfold + EternaFold +
LinearFold round out the RNA secondary-structure folding surface
with three more canonical RNA folders that span the modern
tradeoff space: the original Zuker / Stiegler dynamic-programming
folder that defined the field (mfold/UNAFold, academic-license),
the Eterna game's ML-aware folder reachable via the `arnie`
Python wrapper (EternaFold, MIT), and Baidu / Oregon State's
beam-search linear-time folder (LinearFold, Apache-2.0 — the
folding-only sister to Phase 43 LinearDesign that scales to
viral-genome-length sequences). mfold + LinearFold ride the
single-binary subprocess shape (sister to Phase 18 BWA);
EternaFold rides the Python-script subprocess shape (sister to
Phase 17 Biopython, Phase 28 NUPACK).

| Adapter      | Capability |
|---|---|
| **mfold**       | Michael Zuker's classic RNA / DNA secondary-structure folder (academic-use license) — original dynamic-programming Zuker / Stiegler folder that defined the field; minimum-free-energy folding plus a configurable suboptimal-structure ensemble within an energy budget of the MFE optimum, using the canonical Turner / Mathews thermodynamic parameters; modern `mfold` / `UNAFold.pl` driver consumes a single-sequence FASTA / `.seq` input and writes a classic connect-table `.ct`, plus a PostScript / PDF structure plot and a per-run `.out` log; single-binary subprocess shape sister to Phase 18 BWA with mfold's `KEY=VALUE`-style invocation `mfold SEQ=<sequence> NA=RNA T=<temperature> [extras...]`; knobs `sequence` / `output_basename` / `temperature` (`f64`, finite; default 37.0 — physiological) / `extra_args`; collects `*.ct` (`Tabular`, "mfold connect-table"), `*.ps` / `*.pdf` (`Native`, "mfold structure plot"), `*.out` (`Log`, "mfold log"); probe via `find_on_path(&["mfold", "UNAFold.pl"])` (modern UNAFold distribution renames the launcher to `UNAFold.pl` while keeping the original `mfold` symlink); **academic-license-only** — probe pushes an `"academic"` / `"non-commercial"`-keyworded license-awareness warning in `ProbeReport.warnings` (sister to Phase 28 ViennaRNA / NUPACK / Phase 23 VMD / Phase 5.6 NAMD pattern); version range `3.8.0..4.0.0`; `bio.mfold.fold` ribbon capability |
| **EternaFold**  | the Eterna project's ML-aware RNA folder (MIT) — rides on a half-decade of crowd-sourced Eterna gameplay-puzzle data plus modern thermodynamic + machine-learning corpora to train a maximum-expected-accuracy (MEA) predictor competitive with ViennaRNA + RNAstructure; canonical interface is the Das lab's `arnie` Python wrapper (a unified RNA-folder front-end that proxies to ViennaRNA / RNAstructure / NUPACK / EternaFold / LinearFold under a single Python API); Python-script subprocess shape sister to Phase 17 Biopython / Phase 28 NUPACK / Phase 41 pydna / Phase 43 DNA Chisel; knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_fasta` (`Option<PathBuf>`) / `output_basename`; `prepare()` enforces `.py`, routes script + optional input_fasta through `confined_join`, writes `valenx_params.json` (key omitted when `None`); collects `<output_basename>*.ct` (`Tabular`, "EternaFold connect-table"), `<output_basename>*.dot` (`Native`, "EternaFold dot-bracket"), `<output_basename>*.csv` (`Tabular`, "EternaFold MEA / probabilities"), `*.log`; probe via Python on PATH then `<python> -c "import arnie"` — on import failure surface as a `ProbeReport.warnings` entry, not error; version range `1.3.0..2.0.0`; `bio.eternafold.fold` ribbon capability |
| **LinearFold**  | Baidu Research / Oregon State's beam-search linear-time RNA folder (Apache-2.0) — folding-only sister to Phase 43 LinearDesign, same beam-search core from the same group applied to the inverse problem of "given a sequence, find the secondary structure" rather than "given a target protein, find the optimized mRNA"; decisive property is linear (`O(N)` or `O(N · beam_size)`) per-nucleotide complexity, contrasted with the cubic `O(N^3)` cost of the classical Zuker / RNAstructure / ViennaRNA folders — scales to viral-genome-length sequences (~30 kb SARS-CoV-2, ~9 kb HIV, full-length pre-mRNAs) without polynomial blowup; ships two model back-ends (`C` CONTRAfold-style ML, `V` Vienna-style thermodynamic) selectable through `model`; single-binary subprocess shape with non-standard stdin contract — LinearFold reads the sequence from stdin and writes the predicted dot-bracket / energy lines to stdout; knobs `sequence` (read in place, no staging) / `output_basename` / `model` (default `"C"`) / `beam_size` (`u32`, ≥ 1, default 100) / `extra_args`; `prepare()` resolves `sequence` against the case directory when relative, validates the file exists on disk, composes `linearfold -V <beam_size> [extras...]` (or `-C`) with stdin redirected from `sequence` and stdout redirected to `<output_basename>.txt`; collects `<output_basename>*.txt` (`Tabular`, "LinearFold structure output"), `*.log`; probe via `find_on_path(&["linearfold"])` — surfaces a `"clone https://github.com/LinearFold/LinearFold and add the bin directory to PATH"` warning when Python is on PATH but `linearfold` is missing; version range `1.0.0..2.0.0`; `bio.linearfold.fold` ribbon capability |

### Molecular docking (Phase 34)

| Adapter          | Capability |
|---|---|
| **AutoDock Vina** | modern single-binary small-molecule docker; receptor PDBQT + ligand PDBQT in, ranked-pose PDBQT out; `center` / `size` / `exhaustiveness` / `num_modes` / `energy_range` / `cpu` knobs; collects `Native` artifact `"AutoDock Vina docked poses"` |
| **AutoDock 4**    | two-stage docker (`autogrid4` writes grid maps from `.gpf`, then `autodock4` reads grid + `.dpf` to dock); `skip_grid` toggle to reuse pre-generated maps; per-stage log filenames (`grid_log` / `dock_log`); probe warns if `autogrid4` missing from PATH |

### Phylogenetics (Phase 30)

| Adapter      | Capability |
|---|---|
| **IQ-TREE**   | de-facto modern maximum-likelihood phylogenetic tree builder (GPL-2.0); single-binary `iqtree2 -s <alignment> -m <model> -B <bootstrap> -T <threads> --prefix <prefix>` shape (`-B` omitted when `bootstrap == 0`); `model` defaults `"MFP"` triggering ModelFinder's automatic model selection; `threads` default `"AUTO"` validated against `^(AUTO|\d+)$`; `bootstrap` default 1000 UFBoot ultrafast bootstrap replicates; collects `<prefix>.treefile` (`Native`, "IQ-TREE ML tree") + `.iqtree` / `.log` (`Log`); probe via `find_on_path(&["iqtree2", "iqtree"])` (newer 2.x ships as `iqtree2`; older 1.x as `iqtree`) |
| **RAxML-NG**  | next-generation RAxML rewrite (AGPL-3.0); single-binary `raxml-ng --<mode> --msa <alignment> --model <model> --threads <N> --prefix <prefix>` with mode dispatch on `mode ∈ {search, all, bootstrap}` (`--bs-trees <bootstrap>` appended when mode in `{all, bootstrap}`); `bootstrap ≥ 1` required when `mode ∈ {all, bootstrap}`; collects `<prefix>.raxml.bestTree` (`Native`, "RAxML-NG ML tree") + `<prefix>.raxml.support` (`Native`) + `<prefix>.raxml.log` (`Log`); probe via `find_on_path(&["raxml-ng"])` |
| **FastTree**  | approximate-ML phylogenetic inference (GPL-2.0), optimized for very large trees — sub-quadratic in alignment size; single-binary `FastTree [-nt] [-gtr if use_gtr] [-gamma if gamma] <alignment>` (writes Newick to stdout, captured to `output` via the MAFFT-style stdout-redirect pattern); `seq_type ∈ {nt, aa}` toggles `-nt`; `use_gtr` defaults `true` for nt (FastTree's default is JC without `-gtr`), ignored for aa; collects `output` as `Native` artifact "FastTree Newick tree"; probe via `find_on_path(&["FastTree", "fasttree"])` (binary name varies by distro) |

### Bayesian phylogenetics (Phase 30.5)

Rounds out the molecular phylogenetics surface that Phase 30
opened from the maximum-likelihood side — BEAST 2 + MrBayes
cover the Bayesian MCMC tradeoff space alongside the Phase 30
IQ-TREE / RAxML-NG / FastTree maximum-likelihood beachhead.
Both adapters share the established Phase 18 BWA single-binary
CLI pattern: a user-authored model description (BEAST 2 XML or
MrBayes NEXUS file) in, posterior tree + parameter samples out.

| Adapter      | Capability |
|---|---|
| **BEAST 2**    | the cross-platform Bayesian Evolutionary Analysis by Sampling Trees v2 engine (LGPL-2.1) — canonical Bayesian MCMC framework for time-calibrated phylogenetics: tip-dated trees, relaxed molecular clocks, coalescent demographic models, birth-death speciation models, and the ever-growing universe of BEAST 2 packages (BDSKY, MASCOT, BEASTling, StarBEAST3, ...); single-binary subprocess shape (sister to Phase 18 BWA) with `beast [-seed <N>] -threads <N> [-overwrite] <xml> [extras...]` (XML positional last so BEAST treats it as the model file rather than the value of an earlier flag); knobs `xml` (BEAUti-generated XML model file; required) / `seed` (optional `u64`; passed via `beast -seed <N>` when present, otherwise BEAST picks its own seed and prints it on the run banner) / `threads` (`u32`, ≥ 1, default 1; maps to `-threads N` for tree-likelihood evaluation parallelism) / `overwrite` (default `false`; toggles `-overwrite` so an existing output set from a previous run is replaced rather than triggering a fail) / `extra_args`; collects `*.log` (`Log`, "BEAST 2 trace log" — the parameter trace Tracer reads) and `*.trees` (`Native`, "BEAST 2 sampled trees" — the sampled tree posterior TreeAnnotator / DensiTree consumes); the adapter doesn't try to predict the exact filenames since BEAST writes whatever the XML's `<log fileName="...">` sites configure; probe via `find_on_path(&["beast"])` (the generic version detector tries the conventional `--version` and BEAST's own `-version` form) |
| **MrBayes**    | the long-standing Bayesian MCMC phylogenetic inference engine (GPL-3.0) — historic workhorse for Bayesian phylogenetics; alongside BEAST 2 the de-facto choice for posterior tree sampling across nucleotide / amino-acid / morphological datasets, with its own NEXUS-embedded model-and-mcmc command language and built-in Metropolis-coupled MCMC ("MC^3") chain swapping; single-binary subprocess shape (sister to BEAST 2) with `mb [-i if batch] <nexus> [extras...]` (the binary is literally named `mb` — the project's own convention); knobs `nexus` (NEXUS data file with embedded MRBAYES block driving the run; required) / `batch` (default `false`; toggles `-i` so MrBayes runs the embedded commands non-interactively and exits cleanly rather than waiting on stdin at the prompt — the right default for non-interactive automation) / `extra_args`; collects `*.t` (`Native`, "MrBayes tree samples"), `*.p` (`Tabular`, "MrBayes parameter samples"), and `*.con.tre` (`Native`, "MrBayes consensus tree"); probe via `find_on_path(&["mb"])` |

### Cheminformatics expansion (Phase 24)

| Adapter        | Capability |
|---|---|
| **DeepChem**   | PyTorch-backed deep-learning cheminformatics; user-provided Python script + `valenx_params.json` knobs (inline `smiles` list, optional `dataset_csv`, optional `checkpoint`); collects `.csv` (`Tabular`, "DeepChem analysis output"), `.png` ("DeepChem plot"), `.pkl` / `.pt` ("DeepChem model checkpoint") |
| **Open Babel** | de-facto open-source chemistry-format converter (~120 formats); `obabel <in> -O <out>` single-binary CLI; explicit `input_format` / `output_format` overrides, `gen_3d` toggles `--gen3D` (2D → 3D coords), `add_hydrogens` toggles `-h`; collects converted output as `Native` artifact "Open Babel converted file" |
| **Avogadro 2** | Python-scriptable chemistry editor + small-molecule rendering pipeline; `avogadro2 --script <script.py>` with optional `structure` (`.cml` / `.mol` / `.xyz` / `.pdb`) staged as positional arg; `headless` (default true) toggles `--no-gui`; collects `.png` ("Avogadro 2 render") and `.cml` / `.mol` / `.xyz` ("Avogadro 2 exported structure") |

### Quantum chemistry (Phase 25)

First quantum-chemistry domain to ship in Valenx — Psi4 / NWChem /
xTB span the tradeoff space from semiempirical (xTB, fast and
approximate) through general-purpose HF/DFT/post-HF (Psi4) to
massively-parallel ab initio + plane-wave DFT (NWChem).

| Adapter      | Capability |
|---|---|
| **Psi4**     | open-source HF/DFT/post-HF quantum chemistry (LGPL-3.0); single-binary `psi4 -i <input> -o <output> -n <threads> [-m <memory>] [extras...]` shape (Psithon-scriptable input; `-m` only emitted when `memory` is non-default to preserve Psi4's own `"500 mb"` internal default); knobs `input` (`.in` / `.dat` Psithon script) / `output` / `threads` (default 1, ≥ 1) / `memory` (default `"1 gb"`, matches `^\d+\s*(mb|gb|MB|GB)$` via `is_valid_memory` helper); collects `output` as `Log` artifact "Psi4 output" plus `.fchk` (`Native`, "Psi4 formatted checkpoint") and `.molden` (`Native`, "Psi4 Molden orbital data"); probe via `find_on_path(&["psi4"])` |
| **NWChem**   | PNNL's massively-parallel ab initio + plane-wave DFT package (ECL-2.0); single-binary subprocess with optional MPI wrapping — serial `nwchem [extras...] <input>` or parallel `mpirun -n <mpi_procs> [extras...] nwchem <input>`, both stdout-redirected to `output` (MAFFT-style); knobs `input` (`.nw`) / `output` / `mpi_procs` (default 1, ≥ 1; `prepare()` resolves `mpirun` and fails with a helpful install-hint `InvalidCase` when `mpi_procs > 1` and `mpirun` is missing rather than letting the child fail later); collects `output` as `Log` artifact "NWChem output"; probe via `find_on_path(&["nwchem"])` |
| **xTB**      | Stefan Grimme's extended tight-binding semiempirical method (LGPL-3.0); single-binary `xtb <input> --gfn <gfn> --chrg <charge> --uhf <uhf> [--<mode> if mode != "single-point"] [--alpb <solvent> if Some] [extras...]` shape (writes its run report to stdout, captured to `xtb.log` via the MAFFT-style stdout-redirect pattern); knobs `input` (`.xyz`) / `mode` ∈ `{single-point, opt, ohess, hess, md}` (default `"single-point"` — xTB's default run type, no flag emitted) / `charge` (`i32`, default 0) / `uhf` (xTB's multiplicity convention — number of unpaired electrons, default 0) / `gfn` ∈ `{0, 1, 2}` (default 2 — GFN2-xTB) / `solvent` (optional ALPB solvent name; `None` = gas phase); collects `xtb.log` as `Log` plus `xtbopt.xyz` (`Native`, "xTB optimized geometry"), `xtbopt.log` (`Log`), `gradient` / `hessian` files (`Native`); probe via `find_on_path(&["xtb"])` |

### Workflow managers (Phase 22)

| Adapter        | Capability |
|---|---|
| **Nextflow**   | DSL-driven pipeline orchestrator behind nf-core; `nextflow run <pipeline> [-c <config>] [-profile <profile>] [-resume] [--<key> <value>...]`; `pipeline` accepts a local `.nf`, an absolute path, or a registry identifier like `nf-core/rnaseq`; `params` map maps to `--<key> <value>`; collects the workdir as `Native` artifact `"Nextflow run workdir"` and walks `report.html` / `timeline.html` / `dag.svg` as `Log` artifacts |
| **Snakemake**  | Python-flavoured rule-based pipeline orchestrator; `snakemake -s <snakefile> --cores N [--use-conda] [-n] [--configfile <path>] [<targets>...]`; `targets` lists specific rules to build (empty = default); `cores` (≥ 1) sets parallelism; `use_conda` toggles managed envs; `dry_run` toggles plan-only `-n`; collects the workdir as `Native` artifact `"Snakemake run workdir"` and walks `.snakemake/log/*.log` |

### Workflow expansion (Phase 22.5)

Sister-adapter expansion of the Phase 22 Nextflow + Snakemake
workflow-manager pair — planemo / Cromwell / cwltool round out
the bio workflow-orchestration surface that Phase 22 opened,
covering three more canonical workflow languages: Galaxy
ecosystem CLI for tool development + workflow execution outside
a full Galaxy server (planemo), Broad Institute Workflow
Description Language engine (Cromwell, JAR-distributed), and
Common Workflow Language reference runner (cwltool, Python).

| Adapter      | Capability |
|---|---|
| **planemo**    | Galaxy project's official CLI companion for tool development + workflow execution outside a full Galaxy server (AFL-3.0); single-binary subprocess shape sister to Phase 22 Nextflow / Snakemake with `planemo <action> <workflow> [inputs] [extras...]`; knobs `workflow` (`.ga` / `.gxwf.yml`; required) / `inputs` (`Option<PathBuf>`) / `output_basename` / `action` (default `"run"`; rejected at parse time if not in `{run, test, lint}`) / `extra_args`; `prepare()` resolves both files against case dir, validates each exists; collects `<output_basename>*.html` (`Native`, "Planemo report"), `*.json` (`Tabular`, "Planemo run JSON"), `*.log` (`Log`); probe via `find_on_path(&["planemo"])`; version range `0.75.0..1.0.0`; `bio.planemo.run` ribbon capability |
| **Cromwell**   | Broad Institute's canonical Workflow Description Language (WDL) engine (BSD-3-Clause); **JAR-distributed** — no `cromwell` launcher binary on PATH; the user supplies `[bio.cromwell].jar` absolute path; single-binary subprocess shape sister to Phase 33 j5 / Cello / Phase 41 Jalview with `java -jar <jar> <action> <workflow> [-i <inputs>] [extras...]`; knobs `jar` / `workflow` (`.wdl`) / `inputs` (`Option<PathBuf>` — emitted as TWO separate args `-i` + `<inputs>` only when `Some`) / `output_basename` / `action` (default `"run"`; rejected at parse time if not in `{run, submit, validate}`) / `extra_args`; collects top-level `<output_basename>*.json` (`Tabular`, "Cromwell metadata"), `*.log` (`Log`); probe via `find_on_path(&["java"])`; version range `80.0.0..100.0.0`; `bio.cromwell.run` ribbon capability |
| **cwltool**    | reference implementation of the Common Workflow Language (Apache-2.0); single-binary subprocess shape sister to Phase 22 Snakemake with `cwltool --outdir <output_dir> [extras...] <workflow> [inputs]`; knobs `workflow` (`.cwl`; required) / `inputs` (`Option<PathBuf>` — JSON or YAML CWL input object) / `output_dir` / `extra_args`; `prepare()` resolves both files against case dir; `collect()` walks ONE LEVEL deep into `<output_dir>/` for any file (`Native`, "cwltool output"), top-level `*.log` (`Log`); probe prefers `cwltool` console-script with Python-on-PATH fallback + "cwltool not found on PATH; install via `pip install cwltool`" warning when only Python is present; version range `3.1.0..4.0.0`; `bio.cwltool.run` ribbon capability |

### Systems biology (Phase 32)

First systems-biology / multiscale modeling domain to ship in Valenx —
COPASI / BioNetGen / PhysiCell span the tradeoff space from biochemical
pathway / ODE simulation (COPASI, deterministic) through rule-based
combinatorially-complex signaling networks (BioNetGen, rule expansion +
ODE / SSA) to spatial agent-based multicellular tissue simulation
(PhysiCell, multiscale).

| Adapter      | Capability |
|---|---|
| **COPASI**     | the COmplex PAthway SImulator — de-facto biochemical pathway / ODE-based systems-biology suite, descended from the Gepasi lineage (Artistic-2.0); single-binary `CopasiSE [--save <report>] <model> [--scheduled if run_all] [extras...]` shape (capital `C-S-E` "Self-Executing" CLI; reads `.cps` native or SBML `.xml`); knobs `model` (`.cps` / SBML) / `report` (optional `--save` target so collect() finds the run output without walking) / `run_all` (default false; toggles `--scheduled` to execute every defined task) / `extra_args`; collects the explicit `report` path as `Tabular` ("COPASI report") when supplied, else walks the workdir top-level for `.csv` / `.txt` files; probe via `find_on_path(&["CopasiSE"])` |
| **BioNetGen**  | rule-based modeling language + tool suite for combinatorially-complex signaling networks (MIT); single-binary `BNG2.pl [--no-execute if generate_only] -o <output_basename> <model> [extras...]` shape (Perl driver); knobs `model` (`.bngl`) / `output_basename` (becomes the `-o` prefix every output inherits so collect() walks deterministically) / `generate_only` (default false; adds `--no-execute` to skip simulate / scan / fitting actions and emit just the expanded reaction network) / `extra_args`; collects `<output_basename>*.net` (`Native`, "BioNetGen reaction network"), `<output_basename>*.gdat` (`Tabular`, "BioNetGen species trajectories"), `<output_basename>*.cdat` (`Tabular`, "BioNetGen concentrations") — `parameter_scan` per-trial variants share the basename prefix (e.g. `<basename>_001.gdat`); the init alias `bng` resolves to the same template as the canonical `bionetgen` name; probe via `find_on_path(&["BNG2.pl"])` |
| **PhysiCell**  | Paul Macklin's agent-based, off-lattice multicellular simulator (BSD-3-Clause) — tens to hundreds of thousands of individual cells (each an agent with state, mechanics, secretion, phenotype) coupled to a reaction-diffusion microenvironment for substrates like oxygen and drugs; canonical use case is tumour growth + immunology; PhysiCell models compile per-project to a project-specific C++ binary, so the adapter takes both the user's compiled `binary` path and the run-time XML configuration; knobs `binary` (the per-project compiled executable) / `config` (the `.xml` settings file PhysiCell binaries accept as a positional argument) / `extra_args`; `prepare()` validates `binary` and `config` exist on disk (returns `InvalidCase` with a helpful "PhysiCell models compile per-project — clone the framework, edit the project's `custom_modules/` source, run `make`, and point this field at the resulting executable." message if missing); collects `output/*.xml` and `output/*.mat` (`Native`, "PhysiCell tissue snapshot") plus `output/*.csv` (`Tabular`, "PhysiCell scalar table"); probe via `find_on_path(&["physicell"])` returns `ok = true` either way (most installs won't have a generic `physicell` binary on PATH — the per-project build pattern means there isn't a canonical one) and attaches a warning that the real validation happens at prepare time |

### Synthetic biology (Phase 33)

First synthetic biology / genetic-circuit design domain to ship in
Valenx — pySBOL / j5 / Cello span the synthetic-biology tradeoff
space from canonical SBOL-standard Python composition (pySBOL, the
reference implementation of the Synthetic Biology Open Language for
capturing genetic designs as round-trippable RDF/XML) through DNA
assembly automation that plans the optimal Gibson / Golden-Gate /
SLIC / SLIM strategy from a target circuit + parts library (j5,
JBEI's canonical assembly automator) to genetic-circuit DNA
compilation from a Verilog netlist (Cello v2, the canonical CIDAR
genetic-circuit DNA compiler). j5 + Cello are JAR-distributed —
the user supplies the absolute path to the JAR via case input and
we probe `java` itself rather than the JAR.

| Adapter      | Capability |
|---|---|
| **pySBOL**     | the Python implementation (pySBOL3) of the Synthetic Biology Open Language standard (Apache-2.0); SBOL captures components, sequences, interactions, constraints, and the full provenance of a synthetic design as RDF/XML or JSON-LD that round-trips with every SBOL-conformant tool (j5, Cello, SynBioHub, iBioSim, ...); Python-script subprocess shape (sister to Phase 17 Biopython); knobs `script` / `python` (default `"python3"`) / `input_sbol` (optional starting SBOL XML; `None` when the script generates the design from scratch) / `output_basename`; `prepare()` stages script + optional input SBOL into the workdir under their original filenames so the script can resolve them via relative paths, then writes a flat `valenx_params.json` with `input_sbol` (staged filename or literal `null`) and `output_basename`; collects `<output_basename>*.xml` (`Tabular`, "pySBOL document") and `<output_basename>*.json` (`Log`, "pySBOL composition log"); probe via Python on PATH with an `import sbol3` check (returns `ok = true` with a warning when import fails so non-standard installs aren't blocked); the init alias `sbol` resolves to the same template as the canonical `pysbol` name |
| **j5**         | JBEI's canonical DNA-assembly automation tool (BSD-3-Clause); j5 consumes a target circuit design (CSV row per cassette) plus a parts library (CSV row per part / oligo), then plans the optimal Gibson / Golden-Gate / SLIC / SLIM assembly strategy and writes the per-step protocol + GenBank construct files; **JAR-distributed** — no `j5` launcher binary on PATH; the user supplies the absolute path to `j5.jar` via `[bio.j5].jar` in `case.toml`; single-binary subprocess shape (sister to Phase 18 BWA) with `java -jar <jar> -d <design_csv> -p <parts_csv> -o <output_basename> [extras...]`; knobs `jar` (absolute path to `j5.jar`; required) / `design_csv` (j5 design CSV; required) / `parts_csv` (parts list CSV; required) / `output_basename` / `extra_args`; `prepare()` resolves all three input paths against the case directory when relative, validates each file exists on disk (returns `InvalidCase` with a helpful message when missing), and composes the `java -jar` invocation; collects `<output_basename>*.csv` (`Tabular`, "j5 assembly plan") and `<output_basename>*.gb` (`Native`, "j5 GenBank output"); probe via `find_on_path(&["java"])` (j5's version comes from the jar itself, not from `java`, so we surface no version here — the user pins the j5 release implicitly by the jar they point at) |
| **Cello**      | CIDAR's canonical genetic-circuit DNA compiler (Cello v2, BSD-3-Clause); Cello consumes a Verilog netlist describing the desired logic function plus a triplet of JSON constraint files (a user constraint file pinning the chassis / library, an input sensor file pinning the input promoters, an output device file pinning the reporter) and emits a fully assembled DNA construct that implements the logic in a living cell, running a simulated-annealing optimization over the gate-assignment problem and outputting a Graphviz `.dot` netlist + circuit diagram PNG + human-readable report; **JAR-distributed** — no `cello` launcher binary on PATH; the user supplies the absolute path to the jar via `[bio.cello].jar` in `case.toml`; single-binary subprocess shape (sister to j5) with `java -jar <jar> -inputNetlist <verilog> -targetDataFile <user_constraints> -inputSensorFile <input_sensors> -outputDeviceFile <output_devices> -outputDir <output_basename> [extras...]`; knobs `jar` / `verilog` / `user_constraints` (`.UCF`) / `input_sensors` (`.input.json`) / `output_devices` (`.output.json`) / `output_basename` / `extra_args`; `prepare()` resolves all five input paths against the case directory when relative, validates each file exists on disk, and composes the `java -jar` invocation; collects `<output_basename>*.txt` (`Log`, "Cello report"), `<output_basename>*.png` (`Native`, "Cello circuit diagram"), and `<output_basename>*.dot` (`Native`, "Cello Graphviz netlist"); probe via `find_on_path(&["java"])` (same JAR-versioning shape as j5) |

### Spatial stochastic (Phase 32.5)

Sister-adapter expansion of the Phase 32 systems-biology trio
(COPASI / BioNetGen / PhysiCell) — Smoldyn / MCell round out the
systems-biology / multiscale modeling surface with the canonical
**spatial stochastic / cell-scale reaction-diffusion** simulators
that Phase 32 explicitly deferred. Smoldyn resolves individual
molecules as particles diffusing and reacting in continuous 3D
space; MCell does the same for intricate triangle-mesh geometry
via its own MDL model description language. Both adapters follow
the established Phase 18 BWA single-binary CLI pattern: model
file in, reaction-data + trajectory artifacts out.

| Adapter      | Capability |
|---|---|
| **Smoldyn**    | Steve Andrews's spatial stochastic reaction-diffusion simulator (LGPL-2.1); resolves individual molecules as particles diffusing and reacting in continuous 3D space (no lattice discretisation), the canonical choice when the question is "where does each molecule actually end up over time" rather than "what is the well-mixed concentration vs. t" Phase 32 COPASI's ODE / SSA covers; single-binary subprocess shape (sister to Phase 18 BWA) with `smoldyn <config> [extras...]`; knobs `config` (Smoldyn `.txt` configuration file describing simulation geometry — boundaries, surfaces, compartments — plus molecule species + diffusion coefficients and per-pair / per-surface reactions; required) / `extra_args`; `prepare()` resolves `config` against the case directory when relative, validates it exists on disk; collects `*.txt` (`Tabular`, "Smoldyn output table" — Smoldyn's per-step particle / reaction tables), `*.dat` (`Tabular`, "Smoldyn data" — reaction-event / molecule-position dumps the config may direct here), `*.log` (`Log`, "Smoldyn log"); probe via `find_on_path(&["smoldyn"])`; version range `2.70.0..3.0.0` (Smoldyn 2.70 (2023) is the modern stable line shipping the contemporary surface-reaction model + lattice-Monte-Carlo mode); `bio.smoldyn.simulate` ribbon capability |
| **MCell**      | Salk Institute / Stiles, Bartol's cell-scale Monte Carlo spatial stochastic simulator (GPL-2.0); walks the user's `.mdl` (Model Description Language) model — geometry built from triangle meshes, molecule species with diffusion coefficients, surface / volume reactions, release patterns, observation counts — and runs Brownian-dynamics particle trajectories with Monte Carlo reaction sampling; canonical use case is sub-cellular signaling (synaptic transmission, calcium dynamics, receptor binding) where geometry is intricate enough that Smoldyn's continuous-space mode would be overkill but a well-mixed COPASI / BioNetGen treatment misses the spatial structure; single-binary subprocess shape (sister to Smoldyn) with `mcell [-seed <N>] <mdl> [extras...]` (the `-seed` flag and its integer argument emitted as TWO separate OsString tokens, only when `seed` is `Some(_)`); knobs `mdl` (`.mdl` MCell model description file; required) / `seed` (`Option<u32>` — when `Some(n)` MCell uses that seed; when `None` MCell picks its own seed and prints it on the run banner — same shape as the Phase 29 SLiM `-s` and Phase 30.5 BEAST 2 `-seed` knobs) / `extra_args`; `prepare()` resolves `mdl` against the case directory when relative, validates it exists on disk; collects `*.dat` (`Tabular`, "MCell reaction data" — per-observation count tables MCell writes from the model's REACTION_DATA_OUTPUT block), `*.dx` (`Native`, "MCell visualization data" — DReAMM / OpenDX visualization frames), `*.log` (`Log`, "MCell log"); probe via `find_on_path(&["mcell"])`; version range `4.0.0..5.0.0` (MCell 4.0 (2022) is the modern Python-friendly C++ rewrite); `bio.mcell.simulate` ribbon capability |

### MD analysis expansion (Phase 5.5)

Sister-adapter expansion of the Phase 17 MDAnalysis adapter —
PLUMED / ProDy / cpptraj round out the post-MD analysis surface
that MDAnalysis opened, covering the corners MDAnalysis doesn't
reach: enhanced-sampling collective-variable evaluation + free-
energy reweighting (PLUMED, the de-facto plug-in that wraps every
major MD engine for biased-simulation work), protein-dynamics
elastic-network / normal-mode analysis (ProDy, the canonical Python
toolkit for ENM / GNM / ANM and ensemble PCA), and canonical
AmberTools trajectory analysis via cpptraj's domain language
(cpptraj, the reference workhorse for `rms` / `radgyr` / `hbond`
/ `clustering` over Amber-format trajectories).

| Adapter      | Capability |
|---|---|
| **PLUMED**     | the de-facto enhanced-sampling and free-energy plug-in that wraps every major MD engine — GROMACS, LAMMPS, AMBER, NAMD, OpenMM (LGPL-3.0); defines collective variables (RMSD, dihedrals, distances, contact maps), biases (metadynamics, well-tempered metad, umbrella sampling, ABF), and a reweighting framework that turns biased trajectories back into unbiased free-energy surfaces; the `plumed driver` sub-command runs PLUMED standalone over a pre-computed trajectory: read frames, evaluate the collective variables defined in `plumed.dat`, write COLVAR / bias / HILLS files; single-binary subprocess shape (sister to Phase 18 BWA) with `plumed driver --plumed <plumed_dat> --mf_xtc <trajectory> --kt <kt> [extras...]`; knobs `plumed_dat` (PLUMED input file; required) / `trajectory` (XTC trajectory; required — users running DCD / TRR can swap to `--mf_dcd` / `--mf_trr` via `extra_args`) / `output_basename` / `kt` (`f64`, > 0.0 and finite; PLUMED's `k_B T` in its energy units — kJ/mol by default; default 2.494 = room temperature 300 K; a zero or NaN `kt` would crash PLUMED's reweighting on the first frame) / `extra_args`; collects `<output_basename>*.dat` (`Tabular`, "PLUMED COLVAR output") and `<output_basename>*.bias` (`Tabular`, "PLUMED bias"); probe via `find_on_path(&["plumed"])` |
| **ProDy**      | Bahar lab's canonical Python library for protein dynamics (MIT); ships elastic-network models (ENM / GNM / ANM), normal-mode analysis, ensemble PCA, the NMD trajectory format consumed by VMD's NMWiz plug-in, and integrations with the BLAST / DALI / PDB databases; Python-script subprocess shape (sister to Phase 17 Biopython); knobs `script` / `python` (default `"python3"`) / `input_pdb` / `output_basename` / `num_modes` (`u32`, ≥ 1; number of normal modes to compute; default 20) / `cutoff` (`f64`, > 0.0 and finite; ENM contact cutoff in Å; default 15.0); `prepare()` stages script + input PDB into the workdir under their original filenames so the script can resolve them via relative paths, then writes a flat `valenx_params.json` containing `input_pdb` (staged filename), `output_basename`, `num_modes`, and `cutoff`; collects `<output_basename>*.npz` (`Native`, "ProDy ENM modes"), `<output_basename>*.nmd` (`Native`, "ProDy NMD trajectory" — the NMD format consumed by VMD's NMWiz plug-in for normal-mode visualisation), and `<output_basename>*.csv` (`Tabular`, "ProDy table"); probe via Python on PATH with an `import prody` check (returns `ok = true` with a warning when import fails so non-standard installs aren't blocked) |
| **cpptraj**    | AmberTools' canonical trajectory analysis tool (GPL-3.0); reads Amber `.prmtop` / `.parm7` topologies plus `.nc` / `.dcd` / `.mdcrd` trajectories, runs an analysis script authored in cpptraj's domain language (`trajin`, `rms`, `radgyr`, `hbond`, `volume`, `clustering`, ...), and writes results into the workdir as `.dat` (per-frame tables), `.agr` (XmGrace plot data), or `.gnu` (gnuplot scripts); single-binary subprocess shape (sister to PLUMED) with `cpptraj -p <topology> -i <script> [extras...]`; knobs `script` (`.ptraj` / `.cpptraj` analysis script; required) / `topology` (Amber `.prmtop` / `.parm7`; required) / `extra_args`; `prepare()` resolves both paths against the case directory when relative, validates each file exists on disk, and composes the invocation; collects `*.dat` (`Tabular`, "cpptraj analysis output"), `*.agr` (`Tabular`, "cpptraj XmGrace plot"), and `*.gnu` (`Log`, "cpptraj gnuplot script"); probe via `find_on_path(&["cpptraj"])` |

### Bio MD engines (Phase 5.6)

Sister-domain expansion of the Phase 5 GROMACS / LAMMPS MD engine
beachhead — NAMD / AmberTools sander / HOOMD-blue round out the
all-atom + GPU-native MD-engine surface that Phase 5 GROMACS /
LAMMPS opened, alongside the Phase 17 OpenMM Python-native engine.
NAMD is the de-facto UIUC academic all-atom engine (NAMD-License —
academic / non-commercial use only, flagged via probe warning);
sander is the OSS portion of AMBER's MD engine + canonical
companion to the Phase 5.5 cpptraj analyzer; HOOMD-blue is the
Glotzer-lab GPU-native particle simulator for soft-matter / coarse-
grained work.

| Adapter      | Capability |
|---|---|
| **NAMD**       | UIUC's flagship all-atom MD engine (custom NAMD-License — academic / non-commercial use only, surfaced as `NAMD-License` and flagged via mandatory `"academic"`-keyworded probe warning containing both `"academic"` and `"non-commercial"` substrings) — the de-facto choice in biomolecular MD pedagogy and a workhorse on every academic HPC cluster; NAMD 2.x ships an SMP-threaded CHARMM-style integrator, NAMD 3.x adds GPU-resident kernels; single-binary subprocess shape (sister to Phase 5 LAMMPS / GROMACS) with `<binary> +p<processors> <config> [extras...]` where `<binary>` is `namd2` or `namd3` (probe accepts either) and `+p<N>` (no space — NAMD's own flag syntax) is NAMD's threading flag (multi-threaded SMP build uses it for thread count, MPI build uses MPI-rank count from the launcher); knobs `config` (NAMD `.namd` / `.conf` configuration file; required) / `processors` (`u32`, default 1; emitted as the single OsString `+p<N>` so the flag and value travel together exactly as NAMD parses them) / `extra_args`; `prepare()` resolves `config` against the case directory when relative; collects `*.dcd` (`Native`, "NAMD trajectory (DCD)"), `*.coor` (`Native`, "NAMD coordinates"), `*.vel` (`Native`, "NAMD velocities"), `*.xsc` (`Tabular`, "NAMD extended system"), `*.log` (`Log`); probe via `find_on_path(&["namd2", "namd3"])`; version range `2.14.0..4.0.0` (NAMD 2.14 (2020) is the long-stable line; NAMD 3.x (2022+) is the GPU-resident rewrite) |
| **AmberTools sander** | AMBER's OSS MD engine portion of AmberTools (GPL-3.0 — sander itself is OSS; the proprietary `pmemd.cuda` GPU engine is NOT wrapped here); sander reads an Amber `.prmtop` / `.parm7` topology, an `.inpcrd` / `.rst7` coordinate file, and a `.in` / `.mdin` simulation control file, runs the integrator, and emits an `.out` mdout log + `.rst` restart + `.nc` NetCDF trajectory; sister to Phase 5.5 cpptraj (also AmberTools — installing AmberTools installs both); single-binary subprocess shape (sister to Phase 18 BWA) with `sander -O -i <config> -p <topology> -c <coordinates> -o <basename>.out -r <basename>.rst -x <basename>.nc [extras...]` (the `-O` flag overwrites existing outputs — standard re-run convention); knobs `topology` (`.prmtop` / `.parm7`; required) / `coordinates` (`.inpcrd` / `.rst7`; required) / `config` (`.in` / `.mdin`; required) / `output_basename` (filename stem the adapter pins for the three sander output flags so collect() walks deterministically; required, non-empty) / `extra_args`; `prepare()` resolves all three input paths against the case directory when relative, validates each file exists on disk; collects `<output_basename>*.out` (`Log`, "sander mdout"), `<output_basename>*.nc` (`Native`, "sander NetCDF trajectory"), `<output_basename>*.rst` (`Native`, "sander restart coordinates"), `<output_basename>*.mdinfo` (`Log`, "sander mdinfo"); probe via `find_on_path(&["sander"])` — no academic-license caveat (sander itself is GPL-3.0 OSS); version range `22.0.0..26.0.0` (AmberTools 22 (2022) is the floor we test against) |
| **HOOMD-blue** | Glotzer lab's GPU-native particle simulator (BSD-3-Clause); HOOMD-blue v3+ is fully Python-scripted (no native CLI) — the user supplies a `.py` script that does `import hoomd` and runs the simulation; HOOMD's GPU-resident kernels handle the per-step force evaluation transparently; canonical engine for soft-matter / coarse-grained particle systems, polymers, colloids, rigid-body assemblies — sister to LAMMPS in the particle-MD surface but GPU-first by design; Python-script subprocess shape (sister to Phase 17 OpenMM); knobs `script` (path to user-supplied Python script; required, `.py` enforced) / `python` (interpreter name; default `"python3"`) / `output_basename` (filename stem the user's script uses for outputs — surfaced here so collect() can label artefacts uniformly; required, non-empty); `prepare()` enforces the `.py` extension, stages the script into the workdir under its original filename, then writes a flat `valenx_params.json` containing `output_basename`, and builds `<python> <staged_script>`; collects `<output_basename>*.gsd` (`Native`, "HOOMD trajectory (GSD)"), `<output_basename>*.h5` (`Native`, "HOOMD HDF5 output"), `*.log` (`Log`); probe via `find_on_path(&["python3", "python"])` then `<python> -c "import hoomd"` — on import failure surface as a `ProbeReport.warnings` entry (not error); version range `3.0.0..6.0.0` (HOOMD-blue 3.x (2022) is the modern Python-first rewrite) |

### MDTraj (Phase 5.7)

Single-adapter sister to the Phase 17 MDAnalysis adapter and the
Phase 5.5 PLUMED / ProDy / cpptraj analysis trio — MDTraj rounds
out the post-MD analysis surface with the second-most-used Python
MD trajectory analyzer, with wider format support than MDAnalysis
and deeper integration with the OpenMM ecosystem. Single-adapter
phases are a precedent in Valenx — when an established tool fills
a clearly-defined corner of an existing surface without requiring
new infrastructure, the phase ships as a single adapter.

| Adapter      | Capability |
|---|---|
| **MDTraj**     | Pande / VanderSpoel / Beauchamp lab's Python MD trajectory analysis library (LGPL-2.1) — the second-most-used Python MD trajectory analyzer alongside MDAnalysis; wider format support (`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr` / `.binpos` / `.lh5` / `.amber` / `.gromacs`), deeper integration with the OpenMM ecosystem (the Pande / Beauchamp lab is co-located with the OpenMM developers — MDTraj's HDF5 trajectory format is OpenMM's native streaming output), pandas-friendly per-frame property API; Python-script subprocess shape (sister to Phase 17 Biopython, Phase 5.5 ProDy, Phase 17 OpenMM); knobs `script` (path to user-supplied Python script; required, `.py` enforced) / `python` (interpreter name; default `"python3"`) / `trajectory` (`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr` / `.binpos` / `.lh5` MDTraj-supported trajectory; required) / `topology` (`.pdb` / `.prmtop` / `.gro` / `.psf` topology MDTraj uses for atom + residue + chain metadata; required) / `output_basename`; `prepare()` enforces the `.py` extension on the script, resolves all three input paths against the case directory when relative, stages script + trajectory + topology into the workdir under their original filenames so the script can resolve them via relative paths, then writes a flat `valenx_params.json` containing `output_basename`, the bare `trajectory` filename, and the bare `topology` filename, and builds `<python> <staged_script>`; collects `<output_basename>*.csv` (`Tabular`, "MDTraj analysis table"), `<output_basename>*.npz` (`Native`, "MDTraj numpy archive"), `<output_basename>*.h5` (`Native`, "MDTraj HDF5 output"), `<output_basename>*.png` (`Native`, "MDTraj plot"), `*.log` (`Log`); probe via `find_on_path(&["python3", "python"])` then `<python> -c "import mdtraj"` — on import failure surface as a `ProbeReport.warnings` entry (not error); version range `1.9.0..2.0.0` (MDTraj 1.9 (2022) is the modern stable line that pairs with OpenMM 8.x) |

### Cryo-EM (Phase 36)

First cryo-electron microscopy reconstruction domain to ship in Valenx —
RELION / EMAN2 / CTFFIND span the cryo-EM pipeline from per-micrograph
CTF estimation (CTFFIND, the canonical preprocessing step) through
broad-spectrum image processing (EMAN2, the "Swiss army knife" across
particle picking / 2D classification / initial-model building / 3D
refinement) to Bayesian 3D reconstruction (RELION, the de-facto
single-particle workhorse).

| Adapter      | Capability |
|---|---|
| **RELION**     | Sjors Scheres' REgularised LIkelihood OptimisatioN suite (GPL-2.0) — de-facto Bayesian 3D reconstruction workhorse in cryo-EM facilities worldwide; single-binary subprocess shape with optional MPI wrapping — single-rank `relion_refine --i <particles> --ref <reference> --o <output_basename> --angpix <angpix> --j <threads> [extras...]`, multi-rank `mpirun -n <mpi_procs> relion_refine_mpi ...` (RELION ships separate `_mpi`-suffixed binaries so the launcher knows which transport to use); knobs `particles` (`*_data.star`) / `reference` (`.mrc`) / `output_basename` (becomes the `--o` prefix every output inherits so collect() walks deterministically) / `angpix` (Å, > 0.0 and finite) / `mpi_procs` (default 1, ≥ 1; > 1 switches to MPI binary and prepends `mpirun -n <N>`; surfaces helpful install-hint `InvalidCase` if `mpirun` missing) / `threads` (OpenMP threads per MPI rank, default 1) / `extra_args`; collects `<output_basename>*_class*.mrc` (`Native`, "RELION class average"), `<output_basename>*_data.star` (`Tabular`, "RELION particle assignments"), `<output_basename>*_model.star` (`Log`, "RELION model summary"); probe via `find_on_path(&["relion_refine"])` |
| **EMAN2**      | Steve Ludtke's broad-spectrum cryo-EM image-processing package (BSD-3-Clause) — "Swiss army knife" of single-particle cryo-EM (particle picking, 2D classification, initial-model building, 3D refinement, sprawling Python toolkit `e2*.py`); single-binary subprocess shape wrapping `e2refine_easy.py` (EMAN2's high-level orchestrator that drives the rest of the toolkit); knobs `particles` (`.bdb` / `.hdf` / `.mrcs`) / `model` (`.hdf` / `.mrc` initial 3D model) / `output_basename` (becomes the `--path` argument; EMAN2 turns this into a `<basename>_NN/` results directory under the workdir) / `target_resolution` (Å, > 0.0 and finite) / `symmetry` (point group — `"c1"` / `"d2"` / `"icos"` / etc.; default `"c1"`) / `threads` (default 1) / `extra_args`; collects `<output_basename>_*/threed_*.hdf` (`Native`, "EMAN2 reconstruction") and `<output_basename>_*/log.txt` (`Log`, "EMAN2 log"); the init alias `eman` resolves to the same template as the canonical `eman2` name; probe via `find_on_path(&["e2refine_easy.py"])` |
| **CTFFIND**    | Niko Grigorieff's contrast transfer function (CTF) estimation tool (Janelia non-commercial / academic-only license) — gold standard for fitting per-micrograph CTF parameters (defocus, astigmatism, phase shift); RELION, cryoSPARC, EMAN2, and most automated pipelines all wrap CTFFIND as a preprocessing step; single-binary subprocess shape with stdin-piped parameters (CTFFIND's CLI is interactive and prompts line-by-line for each microscope parameter; the adapter writes a `ctffind_params.txt` file in the workdir during `prepare()` and uses a custom `run()` that pipes the file into the child via `Stdio::from(file)` — the shared `subprocess::run` helper closes stdin which makes CTFFIND read EOF before its first prompt and exit; the custom run path mirrors the MAFFT stdout-redirect pattern but for stdin); knobs `micrograph` (input `.mrc`) / `output_diagnostic` (output diagnostic `.mrc`) / `output_txt` (output text file with CTF parameters) / `pixel_size` (Å, > 0.0 and finite) / `voltage` (kV, default 300.0, > 0.0) / `cs` (spherical aberration mm, default 2.7, > 0.0) / `amplitude_contrast` (fraction in `0.0..=1.0`; 0.07 typical for cryo, 0.1 for negative stain) / `extra_args`; collects `output_diagnostic` (`Native`, "CTFFIND diagnostic image") and `output_txt` (`Tabular`, "CTFFIND parameters"); probe via `find_on_path(&["ctffind"])`; **academic-license-only** — probe pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` and `tool_license` surfaces as `Janelia-License` rather than mislabeling as MIT / BSD |

### Sequencing read simulators (Phase 31)

First sequencing read-simulation domain to ship in Valenx —
ART / wgsim / Badread span all three major sequencing-technology
classes from per-platform empirical-error-profile Illumina short
reads (ART, the de-facto choice) through the simple-uniform-error
classic short-read baseline that ships with samtools (wgsim) to
realistic-error-profile Nanopore long reads (Badread, with calibrated
chimeric / adapter / glitch / identity-drift error models).

| Adapter      | Capability |
|---|---|
| **ART**        | Weichun Huang's NIEHS Illumina-platform read simulator (GPL-3.0) — de-facto choice for synthesising FASTQs that match per-platform empirical error profiles; single-binary subprocess shape wrapping `art_illumina` (companion `art_454` / `art_SOLiD` cover platforms this adapter does not surface); knobs `reference` (FASTA) / `output_prefix` (filename stem; ART writes `<prefix>.fq` single-end or `<prefix>1.fq` + `<prefix>2.fq` paired-end) / `sequencing_system` ∈ `{HS25, HSXt, MSv3, NS50, MinS}` (HiSeq 2500 / HiSeq X TruSeq / MiSeq v3 / NextSeq 500 / MiniSeq) / `read_length` (≥ 1) / `fold_coverage` (> 0.0) / `paired_end` (default `false`) / `fragment_mean` (mean insert size, default 200.0, used iff `paired_end`) / `fragment_sd` (insert-size stddev, default 10.0, used iff `paired_end`) / `extra_args`; `prepare()` builds `art_illumina -ss <sequencing_system> -i <reference> -l <read_length> -f <fold_coverage> -o <output_prefix> [-p -m <fragment_mean> -s <fragment_sd> if paired_end] [extras...]`; collects `<output_prefix>*.fq` (`Tabular`, "ART simulated reads") and `<output_prefix>*.aln` (`Log`, "ART alignment record" — the per-read alignment record ART writes alongside, useful for validating aligner accuracy against the simulated truth); the init alias `art-illumina` resolves to the same template as the canonical `art` name; probe via `find_on_path(&["art_illumina"])` |
| **wgsim**      | Heng Li's classic Whole-Genome SIMulator that ships alongside samtools (MIT) — always paired-end, always position-uniform, deliberately simple under a uniform sequencing-error model; the canonical "small + classic" simulator for fast smoke-testing of mappers and variant callers when realistic error spectra are not required; single-binary subprocess shape (`wgsim` takes the reference and both output FASTQs as positional arguments — no stdout-redirect needed); knobs `reference` (FASTA) / `output1` (FASTQ for read 1) / `output2` (FASTQ for read 2; required — wgsim is paired-end only) / `num_pairs` (≥ 1) / `length1` (read 1 length, default 70) / `length2` (read 2 length, default 70) / `fragment_size` (outer fragment length, default 500) / `error_rate` (per-base error rate in `0.0..=1.0`, default 0.02 — typical Illumina baseline) / `extra_args`; `prepare()` builds `wgsim -N <num_pairs> -1 <length1> -2 <length2> -d <fragment_size> -e <error_rate> <reference> <output1> <output2> [extras...]`; collects `output1` and `output2` as `Tabular` artifacts ("wgsim simulated reads"); probe via `find_on_path(&["wgsim"])` |
| **Badread**    | Ryan Wick's long-read simulator with realistic Nanopore (and PacBio CLR) error profiles (GPL-3.0) — per-platform error models calibrated against actual sequencer output (random / chimeric / adapter / glitch read types, junk-read injection, identity drift, length distributions matching live-flowcell behaviour); single-binary subprocess shape with stdout-redirect (Badread writes its simulated FASTQ to stdout — captured to `output` via the MAFFT-style stdout-redirect-to-file pattern; spawn the child directly and attach stdout to a `File` via `Stdio::from(file)`); knobs `reference` (FASTA) / `output` (FASTQ output path) / `quantity` (Badread `--quantity` literal — one or more decimal digits with optional `K` / `M` / `G` / `T` SI suffix, e.g. `"100M"` for 100 megabases or `"5G"` for 5 gigabases; validated via `is_valid_quantity` helper) / `error_model` ∈ `{nanopore2018, nanopore2020, nanopore2023, pacbio2016}` (per-platform error profile baked into the Badread distribution) / `identity_mean` (read identity mean as a percentage in `0.0..=100.0`, default 87.5) / `length_mean` (read length mean in bases, default 15000.0) / `length_sd` (read length stddev in bases, default 13000.0) / `extra_args`; `prepare()` builds `badread simulate --reference <reference> --quantity <quantity> --error_model <error_model> --identity <identity_mean> --length <length_mean>,<length_sd> [extras...]` → stdout, captured to `output` via the MAFFT-style stdout-redirect pattern; collects `output` as a single `Tabular` artifact ("Badread simulated reads"); probe via `find_on_path(&["badread"])` |

### CRISPR design (Phase 35)

First CRISPR guide-RNA design domain to ship in Valenx —
CHOPCHOP / CRISPOR / Cas-OFFinder span the CRISPR-design tradeoff
space from popular ranked guide design with off-target scoring
(CHOPCHOP, the de-facto first stop in academic CRISPR workflows)
through comprehensive guide design plus rigorous off-target
prediction across many enzymes (CRISPOR, behind the public
crispor.org service) to pure off-target searching used as a
primitive by most CRISPR-design web services and pipelines
(Cas-OFFinder, the OpenCL-accelerated workhorse scanner).

| Adapter      | Capability |
|---|---|
| **CHOPCHOP**     | University of Bergen's web-and-script CRISPR guide-RNA design tool (MIT) — de-facto first stop for "I have a gene, what should I cut" in academic CRISPR workflows; scores candidate gRNAs against a target sequence under a configurable nuclease (Cas9, Cas12a, Cas13) or TALEN design pass, ranks by efficiency / specificity / off-target risk; Python-script subprocess shape (sister to Phase 17 Biopython); knobs `script` (path to user-supplied Python script that imports `chopchop` and reads `valenx_params.json`) / `python` (interpreter name; default `"python3"`) / `target` (target sequence FASTA) / `genome` (CHOPCHOP-installed genome name — `"hg38"` / `"mm10"` / etc.) / `cas_variant` ∈ `{Cas9, Cas12a, Cas13, TALEN}` / `pam` (PAM sequence — `"NGG"` for Cas9, `"TTTV"` for Cas12a, etc.) / `output_basename`; `prepare()` stages script + target FASTA, writes `valenx_params.json` with `target` (staged filename) / `genome` / `cas_variant` / `pam` / `output_basename`, builds `python <script_filename>`; collects `<output_basename>*.tsv` (`Tabular`, "CHOPCHOP guide rankings") and `<output_basename>*.bed` (`Tabular`, "CHOPCHOP guide locations"); probe via Python on PATH with an `import chopchop` check (returns `ok = true` with a warning when import fails so non-standard installs aren't blocked) |
| **CRISPOR**      | Maximilian Haeussler's CRISPR guide-RNA design + off-target prediction tool (GPL-3.0) — distinguishing feature is the rigorous off-target pass via the CFD scoring model and MIT-style specificity scores per guide; powers the public crispor.org service and is also distributed as a standalone Python script for batch / pipeline use; supports many more enzymes / PAMs than CHOPCHOP; Python-script subprocess shape (sister to CHOPCHOP); knobs `script` / `python` / `target` / `genome` / `pam` / `batch_id` (optional — CRISPOR caches partial results by batch so passing the same `batch_id` resumes a previously-interrupted run) / `output_basename`; `prepare()` stages script + target FASTA, writes `valenx_params.json` with `target` (staged filename) / `genome` / `pam` / `batch_id` (JSON string or literal `null`) / `output_basename`, builds `python <script_filename>`; collects `<output_basename>*.tsv` (`Tabular`, "CRISPOR guide rankings") and `<output_basename>*.txt` (`Log`); probe via Python on PATH with an `import crispor` check (same `ok = true` + warning fallback as CHOPCHOP) |
| **Cas-OFFinder** | Bae / Park / Kim group's CRISPR off-target searching tool from Hanyang / Seoul National University (BSD-3-Clause) — fast, OpenCL-accelerated scanner that walks a reference genome and reports every position whose sequence matches one of the input guides within the configured Hamming distance; the workhorse off-target scanner sitting under most CRISPR design web services (CRISPOR, CRISPRdirect, …) and pipelines; single-binary subprocess shape (sister to Phase 18 BWA) with fixed-shape positional CLI `cas-offinder <input> {C\|G\|A} <output> [extras...]` (no `-i` / `-o` flags — the order is fixed); knobs `input` (Cas-OFFinder input file — 3+-line text file with reference genome path, PAM pattern, and one guide-sequence row per query) / `output` (output text file) / `backend` ∈ `{C, G, A}` (OpenCL device class — CPU / GPU / auto-pick fastest at runtime) / `extra_args`; `prepare()` resolves both paths against the case directory when relative and composes the invocation positionally; collects the configured `output` file as a single `Tabular` artifact ("Cas-OFFinder off-target hits"); the init alias `cas-off` resolves to the same template as the canonical `cas-offinder` name; probe via `find_on_path(&["cas-offinder"])` |

### Base + prime editing design (Phase 35.5)

Sister-adapter expansion of the Phase 35 CRISPR design trio
(CHOPCHOP / CRISPOR / Cas-OFFinder) — BE-Designer / BE-Hive /
PrimeDesign / pegFinder round out the CRISPR-design surface with
the modern non-cleavage editing tools that Phase 35 explicitly
deferred: base-editor guide design (BE-Designer), base-editing
outcome prediction (BE-Hive), prime-editing pegRNA design via the
Liu lab (PrimeDesign) and the Komor lab (pegFinder). All four
ride the established Python-script subprocess pattern (sister to
Phase 17 Biopython, Phase 35 CHOPCHOP / CRISPOR, Phase 41 pydna,
Phase 43 DNA Chisel / iCodon).

| Adapter      | Capability |
|---|---|
| **BE-Designer**  | Komor lab's base-editor guide design tool (MIT) — de-facto first stop in modern base-editing workflows; given a target genome region with a desired C→T or A→G base change, enumerates candidate gRNAs, scores each against the editing window of the requested base-editor variant (BE3, BE4max, ABE7.10, ABEmax), and emits a guide table with per-guide editing-window predictions, off-target predictions, and PAM-compatibility filtering; Python-script subprocess shape sister to Phase 35 CHOPCHOP; knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_fasta` (`Option<PathBuf>`) / `output_basename`; `prepare()` enforces `.py`, routes script + optional input_fasta through `confined_join`, writes `valenx_params.json` with `output_basename` always plus `input_fasta` (staged filename) only when set — key omitted when `None`; collects `<output_basename>*.csv` (`Tabular`, "BE-Designer guide table"), `<output_basename>*.fasta` (`Native`, "BE-Designer designed sequences"), `*.log`; probe via Python on PATH then `<python> -c "import bedesigner"`; version range `1.0.0..2.0.0`; `bio.be_designer.design` ribbon capability |
| **BE-Hive**      | Liu lab's base-editing outcome predictor (MIT) — answers the canonical sister question to "what guides should I order": given a designed guide / target / base-editor chassis, what fraction of edits will land on the right base, what fraction will have unintended bystander edits, where in the editing window will edits concentrate; the Liu lab's large-scale base-editing outcome dataset trains a CNN-style model that predicts per-position editing efficiency for every major base-editor chassis (BE3 / BE4 / ABE / CBE families); Python-script subprocess shape sister to BE-Designer; knobs and `prepare()` shape identical to BE-Designer; collects `<output_basename>*.csv` (`Tabular`, "BE-Hive efficiency predictions"), `<output_basename>*.png` (`Native`, "BE-Hive plot"), `*.log`; probe via Python on PATH then `<python> -c "import be_predict"` (the canonical pip-installable package on PyPI is `be_predict`, not `behive`); version range `1.0.0..2.0.0`; `bio.be_hive.predict` ribbon capability |
| **PrimeDesign**  | Liu lab's prime-editing design tool (MIT) — canonical pegRNA designer for the Anzalone / Liu prime-editing system; given a desired edit at a target locus (point mutation, small insertion, small deletion, or combinations up to a few dozen base pairs) it designs the pegRNA — the chimeric guide RNA with a 3′ extension encoding the reverse-transcriptase template + primer-binding site — plus the optional secondary nicking guide that boosts editing efficiency under PE3 / PE3b; Python-script subprocess shape sister to BE-Designer / BE-Hive; knobs and `prepare()` shape identical to BE-Designer; collects `<output_basename>*.csv` (`Tabular`, "PrimeDesign pegRNA table"), `<output_basename>*.txt` (`Tabular`, "PrimeDesign report"), `*.log`; probe via Python on PATH then `<python> -c "import primedesign"`; version range `1.0.0..2.0.0`; `bio.primedesign.design` ribbon capability |
| **pegFinder**    | Komor lab's alternative pegRNA finder (MIT) — sister to PrimeDesign with a different scoring model that emphasizes pegRNA secondary-structure stability + RT-template-length tradeoffs; the Komor lab's analysis showed pegRNA secondary structure substantially impacts prime-editing efficiency, and pegFinder explicitly scores candidate pegRNAs on predicted RT-template + PBS thermodynamic stability; having both PrimeDesign and pegFinder lets users cross-check pegRNA recommendations between the Liu and Komor labs' design philosophies; Python-script subprocess shape sister to PrimeDesign; knobs and `prepare()` shape identical to BE-Designer; collects `<output_basename>*.csv` (`Tabular`, "pegFinder pegRNA candidates"), `<output_basename>*.txt` (`Tabular`, "pegFinder summary"), `*.log`; probe via Python on PATH then `<python> -c "import pegfinder"`; version range `1.0.0..2.0.0`; `bio.pegfinder.design` ribbon capability |

### Edit-outcome prediction (Phase 35.6)

Sister-adapter expansion of the Phase 35 / 35.5 CRISPR-design +
editing surface — inDelphi / FORECasT / AlphaMissense / CRISPRitz
close the design → predict-outcome → off-target loop with four
canonical outcome predictors: Cas9-cut indel pattern prediction
(inDelphi + FORECasT for cross-checking between two independent
labs), missense pathogenicity prediction (AlphaMissense for "will
this missense edit cause disease"), and variant-aware off-target
genome-wide search (CRISPRitz). AlphaMissense ships under the
**CC-BY-NC-SA-4.0 academic / non-commercial weights license**
(sister to AlphaFold 3) — the probe surfaces a mandatory
`"academic"` / `"non-commercial"`-keyworded license-awareness
warning whenever Python is detected, regardless of whether
`import alphamissense` succeeds.

| Adapter         | Capability |
|---|---|
| **inDelphi**     | Liu lab's Cas9-cut indel pattern predictor (MIT) — de-facto first stop in modern Cas9-cut-outcome workflows; given a designed gRNA + target site + cell-line chassis (mESC / U2OS / HCT116 / HEK293), predicts the per-indel-pattern frequency distribution that will result from the Cas9 double-strand break and the cell's own end-joining repair machinery (NHEJ / MMEJ); Python-script subprocess shape sister to Phase 35 CHOPCHOP / Phase 35.5 BE-Designer; knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_fasta` (`Option<PathBuf>`) / `output_basename`; `prepare()` enforces `.py`, routes script + optional input_fasta through `confined_join`, writes `valenx_params.json` (key omitted when `None`); collects `<output_basename>*.csv` (`Tabular`, "inDelphi indel predictions"), `<output_basename>*.png` (`Native`, "inDelphi plot"), `*.log`; probe via Python on PATH then `<python> -c "import inDelphi"`; version range `1.0.0..2.0.0`; `bio.indelphi.predict` ribbon capability |
| **FORECasT**     | Sanger Institute's alternative Cas9-cut indel predictor (Apache-2.0) — different model architecture, trained on a different corpus (Felicity Allen's Sanger SelfTarget library), validated against a different assay design; lets users cross-check Cas9-cut outcome predictions across two independent groups' models — the canonical use case is "inDelphi predicted X% frame-shift; what does FORECasT predict for the same guide / target / chassis combination, and how do the two predictors agree on the dominant indel pattern"; the Python module is published under the SelfTarget name (`import selftarget`) — FORECasT is the predictor's published name in the original 2018 _Nature Biotechnology_ paper but the GitHub repository the canonical pip-install lives in is named `SelfTarget` after Allen's data-collection assay; Python-script subprocess shape sister to inDelphi; knobs and `prepare()` shape identical to inDelphi; collects `<output_basename>*.csv` (`Tabular`, "FORECasT indel predictions"), `<output_basename>*.txt` (`Tabular`, "FORECasT summary"), `*.log`; probe via Python on PATH then `<python> -c "import selftarget"`; version range `1.0.0..2.0.0`; `bio.forecast.predict` ribbon capability |
| **AlphaMissense** | DeepMind's missense-effect predictor (CC-BY-NC-SA-4.0 / academic non-commercial) — extends the AlphaFold structural-prediction lineage to score per-position missense-mutation pathogenicity: given a protein sequence + a missense change, predicts the pathogenicity score on a 0-to-1 scale that calibrates against ClinVar pathogenic / benign labels; the canonical readout in modern CRISPR-editing workflows is "given my designed CRISPR edit causing a missense change at this protein position, what's the predicted pathogenicity"; Python-script subprocess shape sister to inDelphi; knobs and `prepare()` shape identical to inDelphi; collects `<output_basename>*.csv` / `*.tsv` (`Tabular`, "AlphaMissense pathogenicity scores"), `<output_basename>*.png` (`Native`, "AlphaMissense plot"), `*.log`; **academic-license-only** — probe pushes a mandatory `"academic"` / `"non-commercial"`-keyworded warning into `ProbeReport.warnings` whenever Python is on PATH, regardless of whether `import alphamissense` succeeds (sister to the Phase 17.5 AlphaFold 3 mandatory-probe-warning pattern); version range `1.0.0..2.0.0`; `bio.alphamissense.predict` ribbon capability |
| **CRISPRitz**    | Pinello lab's variant-aware off-target genome-wide search (MIT) — sister to Phase 35 Cas-OFFinder with a different scoring model and the distinguishing property of **variant-aware off-target searching**: given a reference genome plus a population VCF (1000 Genomes, gnomAD), CRISPRitz searches for off-target sites that exist only in specific haplotypes / specific population sub-groups rather than just the reference assembly; the canonical readout for "is my CRISPR guide safe across the human population, or are there off-targets that would only manifest in certain genetic backgrounds"; Python-script subprocess shape sister to inDelphi / FORECasT; knobs and `prepare()` shape identical to inDelphi; collects `<output_basename>*.txt` (`Tabular`, "CRISPRitz off-target table"), `<output_basename>*.bed` (`Tabular`, "CRISPRitz off-target BED" — UCSC BED-format off-target hits ready for genome-browser visualisation), `*.log`; probe via Python on PATH then `<python> -c "import crispritz"`; version range `2.6.0..3.0.0`; `bio.crispritz.search` ribbon capability |

### Population genetics (Phase 29)

First population-genetics / evolutionary-simulation domain to ship
in Valenx — SLiM / msprime / tskit span the population-genetics
tradeoff space from forward-time individual-based simulation under
arbitrary selection / demography / mating-system specifications
(SLiM, the de-facto forward simulator) through coalescent
backward-time simulation of sample ancestries (msprime, the de-
facto coalescent simulator and companion to SLiM) to tree-sequence
analysis / statistics on the succinct tree-sequence outputs both
simulators emit (tskit, the canonical analysis library).

| Adapter      | Capability |
|---|---|
| **SLiM**       | Philipp Messer's forward-time population-genetics simulator (GPL-3.0) — evolves a finite-population model generation by generation under a user-defined Eidos script (mutation rates, selection coefficients, recombination maps, demographic events, migrations, mating systems); tree-sequence recording (`treeSeqOutput()` family) feeds straight into tskit / msprime downstream; single-binary subprocess shape (sister to Phase 18 BWA) with `slim [-s <seed>] [extras...] <script>` (script positional last so SLiM treats it as the model file rather than the value of an earlier flag); knobs `script` (`.slim` Eidos model file; required) / `seed` (optional `u64`; passed via `slim -s <N>` when present, otherwise SLiM picks its own seed and prints it on the run banner) / `output_basename` (filename stem the user's script uses for outputs — surfaced here so collect() can label artefacts uniformly even though SLiM scripts choose their own output paths) / `extra_args` (additional CLI arguments appended after the script path; `-d KEY=VALUE` is the canonical way to inject Eidos constants from outside the script); collects `<output_basename>*.trees` (`Native`, "SLiM tree sequence") and `<output_basename>*.log` (`Log`); probe via `find_on_path(&["slim"])` |
| **msprime**    | Jerome Kelleher's coalescent backwards-in-time population-genetics simulator (GPL-3.0) — speed-of-light coalescent simulator (millions of samples per minute on a workstation); the canonical companion to SLiM (forward-time) and tskit (tree-sequence analysis); Python-script subprocess shape (sister to Phase 17 Biopython); knobs `script` / `python` (default `"python3"`) / `population_size` (`u32`, ≥ 1) / `num_samples` (`u32`, ≥ 1) / `recombination_rate` (`f64`, ≥ 0.0 and finite — per-site per-generation rate) / `mutation_rate` (`f64`, ≥ 0.0 and finite) / `output_basename`; `prepare()` stages the script + writes `valenx_params.json` containing `population_size` / `num_samples` / `recombination_rate` (emitted as `{:e}` so Python `json.load` parses as float) / `mutation_rate` (same) / `output_basename`; collects `<output_basename>.trees` (`Native`, "msprime tree sequence"), `<output_basename>.vcf` (`Tabular`, "msprime VCF"), `<output_basename>.csv` (`Tabular`, "msprime per-sample summary"); probe via Python on PATH with an `import msprime` check (returns `ok = true` with a warning when import fails so non-standard installs aren't blocked) |
| **tskit**      | the canonical tree-sequence analysis library (MIT) — built around the succinct tree-sequence data structure pioneered by msprime; computes population-genetics statistics (π, Tajima's D, Fst, site-frequency spectra, IBD shares), exposes per-tree iteration across the genome, converts between tree-sequence and VCF / Newick / table formats, renders phylogenetic plots; the workhorse downstream of every Phase 29 simulator — msprime emits `.trees`, SLiM emits `.trees`, tskit consumes them; Python-script subprocess shape (sister to msprime); knobs `script` / `python` / `input_trees` (`.trees` file from SLiM or msprime; required) / `output_basename`; `prepare()` stages script + tree-sequence file under their original filenames so the script can resolve them via relative paths, then writes `valenx_params.json` with `input_trees` (staged filename) / `output_basename`; collects `<output_basename>*.csv` / `<output_basename>*.tsv` (`Tabular`, "tskit statistics") and `*.png` (`Native`, "tskit plot"); probe via Python on PATH with an `import tskit` check (same `ok = true` + warning fallback as msprime) |

### Rosetta family (Phase 38)

First Rosetta protein-modeling family to ship in Valenx —
Rosetta + PyRosetta open the canonical RosettaCommons code base
through the two most-used entry points: `rosetta_scripts` (the
XML-driven protocol runner that's the de-facto Rosetta entry
point in production — every `relax` / `dock` / `abinitio` /
FastDesign / enzyme-design pipeline lives as an XML protocol fed
to this binary) and PyRosetta (Python bindings exposing the same
C++ core through a Pythonic API for users who prefer scripting
Rosetta from `.py` rather than authoring XML protocols). Both
adapters are **academic-license-flagged** — the RosettaCommons
license is a custom non-OSS academic / non-commercial-use
agreement that Valenx surfaces accurately via
`tool_license = "Rosetta-License"` and a mandatory
`"academic"`-keyworded probe warning.

| Adapter      | Capability |
|---|---|
| **Rosetta**    | RosettaCommons' flagship modeling suite (custom Rosetta-License — academic / non-commercial use only) — drives protein design, structure prediction, docking, ligand binding, and a long tail of related modeling tasks through `rosetta_scripts`, which reads an XML protocol describing the modeling pipeline (filters, movers, scorefunctions) and applies it to an input `.pdb`; single-binary subprocess shape (sister to Phase 18 BWA) with `rosetta_scripts -database <path> -parser:protocol <protocol> -in:file:s <input_pdb> -out:prefix <output_basename> -nstruct <N> [extras...]`; knobs `protocol` (XML protocol script; required) / `input_pdb` / `output_basename` / `nstruct` (`u32`, ≥ 1 — number of independent decoys to generate) / `database` (path to the Rosetta `database/` directory, required because every `rosetta_scripts` invocation needs `-database <path>` pointing at the energy tables / fragment libraries / etc. bundled with the source distribution) / `extra_args`; collects `<output_basename>*.pdb` (`Native`, "Rosetta designed structure") and the canonical `score.sc` scorefile (`Tabular`, "Rosetta scores"); probe via `find_on_path(&["rosetta_scripts", "rosetta_scripts.linuxgccrelease", "rosetta_scripts.macosclangrelease"])` (Rosetta source builds emit platform-suffixed names by default, conda / packaged distributions install a bare `rosetta_scripts` shim — the probe covers all three); **academic-license-only** — probe pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` whenever the binary is detected, and `tool_license` surfaces as `"Rosetta-License"` rather than mislabeling the custom RosettaCommons terms as a recognised SPDX identifier |
| **PyRosetta**  | Python bindings to the Rosetta C++ core (Rosetta-License — inherits the same academic / non-commercial use terms) — exposes the entire Rosetta modeling pipeline (movers, filters, scorefunctions, task-operations) through a Pythonic API, letting users drive Rosetta from regular `.py` scripts rather than authoring XML protocols; Python-script subprocess shape (sister to Phase 17 Biopython); knobs `script` (path to user-authored Python script; required) / `python` (interpreter name; default `"python3"`) / `input_pdb` (optional input PDB the script will operate on — None when the script generates structures de novo; surfaced in `valenx_params.json` so the script can read it without re-parsing case.toml) / `output_basename`; `prepare()` stages script + optional input PDB into the workdir under their original filenames, then writes a flat `valenx_params.json` with `input_pdb` (staged filename or literal `null`) and `output_basename`; collects `<output_basename>*.pdb` (`Native`, "PyRosetta designed structure") and `*.sc` files (`Tabular`, "PyRosetta scores"); probe via Python on PATH with an `import pyrosetta` check (returns `ok = true` with a warning when import fails so non-standard installs aren't blocked); **academic-license-only** — probe always pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` whenever Python is detected (regardless of whether `pyrosetta` itself is importable, since the user is either about to install it or has it installed and needs reminding), and `tool_license` surfaces as `"Rosetta-License"` rather than mislabeling the inherited terms as MIT / BSD |

### DNA structural geometry (Phase 39)

First DNA structural-geometry domain to ship in Valenx — X3DNA /
Curves+ / DSSR span the structural-geometry tradeoff space from
canonical base-pair / base-step parameter calculation (X3DNA, the
de-facto reference for twist / roll / tilt / slide / shift / rise
plus the per-base intra-pair parameters) through helical-axis
curvature analysis (Curves+, the canonical "is this DNA bent"
tool) to structural-feature annotation as a single machine-
readable JSON summary (DSSR, the modern X3DNA-family Python-
fronted tool). All three are **academic-license-flagged** — the
X3DNA family (X3DNA + DSSR) ships under custom non-OSS academic
terms, Curves+ ships under a custom non-OSS academic license —
Valenx surfaces these accurately via `tool_license =
"X3DNA-License"` / `"Curves-License"` / `"DSSR-License"` and a
mandatory `"academic"`-keyworded probe warning whenever each
binary is detected.

| Adapter      | Capability |
|---|---|
| **X3DNA**      | Wilma Olson and Xiang-Jun Lu's reference toolkit for DNA / RNA structural-geometry analysis (custom X3DNA-License — academic / non-commercial use only) — reads a nucleic-acid PDB, identifies base pairs, and computes the canonical helical-step parameters (twist, roll, tilt, slide, shift, rise) plus per-base intra-pair parameters (buckle, propeller, opening, shear, stretch, stagger); single-binary subprocess shape (sister to Phase 18 BWA) with `analyze <input_pdb> [extras...]` (positional-only — `analyze` derives every output filename from the input basename, so the adapter just hands it the PDB and any user-supplied extras); knobs `input_pdb` (required) / `output_basename` (filename stem the user expects X3DNA to produce — surfaced here so collect() can label artefacts uniformly without scraping `analyze`'s filename heuristics) / `extra_args`; collects `<output_basename>*.par` (`Tabular`, "X3DNA base-step parameters") and `*.out` (`Log`, the per-run log `analyze` writes alongside); probe via `find_on_path(&["analyze"])` (X3DNA's main analysis binary is literally named `analyze`); **academic-license-only** — probe always pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` whenever the binary is detected, and `tool_license` surfaces as `"X3DNA-License"` rather than mislabeling the custom X3DNA terms as a recognised SPDX identifier |
| **Curves+**    | Richard Lavery's reference toolkit for DNA helical-axis analysis (custom Curves-License — academic / non-commercial use only) — fits a curvilinear helical axis through a nucleic-acid structure and reports per-base axis-curvature, base-pair parameters relative to that axis, and a `.cda` file describing the axis itself for downstream visualisation; the canonical tool for "is this DNA bent, and if so, how" questions in protein-DNA / drug-DNA structural studies; single-binary subprocess shape with stdin-piped parameters (sister to Phase 36 CTFFIND) — Curves+ takes its parameters as a Fortran-style `&inp ... &end` namelist block on stdin followed by strand / axis residue cards; the adapter writes a `curves_params.txt` file in the workdir during `prepare()` and uses a custom `run()` that opens the file and pipes its contents into the child via `Stdio::from(file)` — the shared `subprocess::run` helper closes stdin which makes Curves+ read EOF before parsing its first parameter and exit; knobs `input_pdb` (required) / `output_basename` (filename stem Curves+ uses for outputs — `<basename>.lis`, `<basename>.cda`, etc.) / `first_residue` (`u32` — first inclusive residue index in the strand to analyse) / `last_residue` (`u32`, ≥ `first_residue` — a reverse range is rejected up front with a helpful message) / `extra_args`; collects `<output_basename>*.lis` (`Log`, "Curves+ helical analysis") and `<output_basename>*.cda` (`Tabular`, "Curves+ axis curve data"); probe via `find_on_path(&["Cur+"])` (the binary name uses a literal `+`); **academic-license-only** — probe always pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` whenever the binary is detected, and `tool_license` surfaces as `"Curves-License"` rather than mislabeling the custom Curves+ terms as a recognised SPDX identifier |
| **DSSR**       | Dissecting the Spatial Structure of RNA / DNA — the modern Python-fronted X3DNA-family tool (custom DSSR-License — academic / non-commercial use only) — reads a nucleic-acid PDB and emits a single JSON file enumerating every detected feature: base pairs (Watson-Crick, Hoogsteen, sugar-edge, ...), multiplets, double helices, stems, hairpin / internal / junction loops, kissing loops, A-minor motifs, ribose zippers, pseudoknots, splayed-apart conformations, and more; the standard machine-readable feature-extraction step in modern RNA-structure pipelines; single-binary subprocess shape (sister to X3DNA) with `x3dna-dssr -i=<input_pdb> -o=<output_json> --json [extras...]` (DSSR uses `key=value` flag form on its short-form options — no space between flag and value); knobs `input_pdb` (required) / `output_json` (output JSON path; required) / `extra_args`; collects the configured `output_json` file as a single `Tabular` artifact ("DSSR analysis (JSON)") — DSSR's JSON is the canonical machine-readable summary; tagged `Tabular` rather than `Native` so downstream serdes can key off a consistent kind; probe via `find_on_path(&["x3dna-dssr"])`; **academic-license-only** — probe always pushes an `"academic"`-keyworded warning into `ProbeReport.warnings` whenever the binary is detected, and `tool_license` surfaces as `"DSSR-License"` rather than mislabeling the inherited X3DNA-family terms as a recognised SPDX identifier |

### Microscopy (Phase 40)

First microscopy / bioimage analysis domain to ship in Valenx —
Fiji / CellProfiler / Ilastik span the bioimage analysis tradeoff
space from script-driven general-purpose image processing in
headless mode (Fiji, the canonical ImageJ distribution) through
pipeline-driven cell segmentation + measurement (CellProfiler,
the Broad Institute pipeline-driven workhorse that powers most
high-content-screening assays) to interactive-ML pixel / object
classification (Ilastik, the Hamprecht lab tool that leans on
user-trained random-forest classifiers for hard segmentation
tasks where rule-based pipelines struggle).

| Adapter      | Capability |
|---|---|
| **Fiji**       | the [Fiji Is Just ImageJ](https://fiji.sc/) distribution of NIH ImageJ (Schindelin et al, GPL-3.0); bundles ImageJ2 + a curated set of plugins for biological image processing — channel splitting, thresholding, particle analysis, deconvolution, registration, segmentation, the TrackMate single-particle tracker, the BoneJ trabecular-bone toolkit, the entire ImageJ macro / Jython / scripting surface; per-platform launcher binaries (`ImageJ-linux64`, `ImageJ.exe`, `Contents/MacOS/ImageJ-macosx`); app-launcher subprocess shape (sister to Phase 36 RELION / EMAN2) with `<fiji_app> --headless --console -macro <macro_file> [extras...]`; knobs `fiji_app` (absolute path to per-platform Fiji launcher; required) / `macro_file` (`.ijm` Fiji macro; required) / `input_image` (`Option<PathBuf>` — optional input image the macro will operate on, typically picked up via `getArgument()`; the macro is responsible for opening it) / `output_basename` / `extra_args`; collects `<output_basename>*.tif` / `.tiff` (`Native`, "Fiji image (TIFF)"), `<output_basename>*.png` (`Native`, "Fiji image (PNG)"), `<output_basename>*.csv` (`Tabular`, "Fiji measurements"), `*.log`; probe via `find_on_path(&["ImageJ-linux64", "ImageJ-macosx", "ImageJ.exe", "fiji"])` — surfaces a Java-fallback warning whenever Fiji isn't on PATH but `java` is; version range `2.0.0..3.0.0` (Fiji's ImageJ2-based 2.x line is the modern stable shipping the curated bundle); `bio.fiji.process` ribbon capability |
| **CellProfiler** | Broad Institute's pipeline-driven cell segmentation + measurement suite (BSD-3-Clause); canonical tool for high-content screening — the user authors a `.cppipe` pipeline in the GUI (a chain of modules like LoadImages → IdentifyPrimaryObjects → MeasureObjectShape → ExportToSpreadsheet) and the CLI runs that pipeline over a directory of input images, emitting per-object measurement CSVs + segmented label-image overlays; Python-CLI subprocess shape with bundled `cellprofiler` launcher as the primary entry point and `<python> -m cellprofiler ...` as the fallback when the launcher isn't on PATH but Python is; `cellprofiler -c -r -p <pipeline> -i <input_dir> -o <basename> [extras...]` (`-c` = run without GUI, `-r` = run pipeline immediately, `-p` = pipeline file, `-i` = input directory, `-o` = output directory); knobs `pipeline` (`.cppipe` / `.cpproj`; required) / `input_dir` (input image directory; required — the adapter validates it is a directory at prepare time) / `output_basename` / `python` (default `"python3"`) / `extra_args`; `prepare()` resolves `pipeline` and `input_dir` against the case directory when relative, validates `pipeline` exists on disk and `input_dir` is a directory; collects walks **one level deep** into `<output_basename>/` for `*.csv` (`Tabular`, "CellProfiler measurements"), `*.tif` / `.tiff` (`Native`, "CellProfiler segmented image"), `*.png` (`Native`, "CellProfiler plot"); top-level `*.log` (`Log`, "CellProfiler log"); probe via `find_on_path(&["cellprofiler", "python3", "python"])` with warning when `cellprofiler` isn't on PATH but Python is; version range `4.0.0..5.0.0` (CellProfiler 4.x is the modern Python 3 line); `bio.cellprofiler.segment` ribbon capability |
| **Ilastik**    | Hamprecht lab's interactive-ML pixel / object classification suite (GPL-3.0); leans on user-trained random-forest classifiers — the user paints a few foreground / background strokes per image in the GUI to teach the classifier, saves the resulting `.ilp` project file, and then runs the headless CLI to apply that trained classifier to a batch of new images; canonical use case is hard segmentation tasks where rule-based pipelines (CellProfiler) or threshold-driven macros (Fiji) struggle — light-sheet imagery, tissue cross-sections, anything with low contrast or irregular textures; app-launcher subprocess shape (sister to Fiji) with `<ilastik_app> --headless --project=<project> --output_filename_format=<basename>_{nickname}.h5 <input_images...> [extras...]`; the `--project=` and `--output_filename_format=` flags are emitted as single OsString args each (so `=` and the value travel together, matching Ilastik's own argv parser), and the literal `{nickname}` substring in the format string is Ilastik's per-image nickname placeholder — must reach Ilastik unmodified for per-input-image output disambiguation; knobs `ilastik_app` (absolute path to per-platform Ilastik launcher; required) / `project` (`.ilp` Ilastik project file containing the trained classifier; required) / `input_images` (`Vec<PathBuf>` — must contain ≥ 1 entry; the adapter rejects an empty vector at prepare time) / `output_basename` / `workflow` (default `"Pixel Classification"` — selectable from Ilastik's set: `"Pixel Classification"`, `"Object Classification"`, etc.) / `extra_args`; collects `<output_basename>*.h5` (`Native`, "Ilastik probability map (HDF5)"), `<output_basename>*.tif` (`Native`, "Ilastik segmentation"), `*.log`; probe via `find_on_path(&["ilastik", "run_ilastik.sh", "ilastik.exe"])` with warning when nothing matches but still returns `ok = true` since the user can supply the launcher via `case.toml` (sister to Phase 32 PhysiCell's per-project-binary probe convention); version range `1.4.0..2.0.0` (Ilastik 1.4 is the modern stable line shipping the contemporary headless mode + workflow set); `bio.ilastik.classify` ribbon capability |

### Sequence editors (Phase 41)

First plasmid-design / alignment-viewer domain to ship in Valenx
— pydna / Jalview span the sequence-editor tradeoff space from a
Python plasmid-design library that handles PCR primer design,
restriction-enzyme digests, and Gibson / Golden-Gate assembly
programmatically (pydna, Bjorn Johansson's BSD-3-Clause library
that's the de-facto Python choice for cloning automation) to the
canonical Java alignment viewer with a headless mode for batch
image / format conversion (Jalview, the Barton group's GPL-3.0
viewer that's been the reference alignment viewer in molecular
biology labs since the 2000s and supports headless operation for
unattended pipeline integration).

| Adapter      | Capability |
|---|---|
| **pydna**      | Bjorn Johansson's Python plasmid / clone-design library (BSD-3-Clause); handles the long tail of cloning operations programmatically — PCR primer design, restriction-enzyme digests, Gibson assembly, Golden-Gate assembly, sequence-overlap detection, ligation simulation, primer Tm calculation, cloning-strategy validation; canonical use case is "I need to assemble these N parts into this target construct — what primers should I order, what enzymes should I cut with, and does the end-product match my target sequence?" — pydna replaces hours of manual work in ApE / SnapGene / Vector NTI with a few dozen lines of Python; Python-script subprocess shape (sister to Phase 17 Biopython, Phase 19.5 Scanpy, Phase 33 pySBOL); knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_genbank` (`Option<PathBuf>` — optional starting GenBank file the script can use as the parent / template construct; `None` when the script generates the design from scratch) / `output_basename`; `prepare()` enforces `.py`, stages script + optional input_genbank, writes `valenx_params.json` with `output_basename` always plus `input_genbank` (staged filename) only when set — key omitted entirely when `None` rather than emitted as `null`, matching the hand-rolled JSON convention the rest of the bio adapters use (Phase 19.6 Seurat / AnnData, Phase 27.5 ESM-IF); collects `<output_basename>*.gb` (`Native`, "pydna GenBank file"), `<output_basename>*.genbank` (`Native`, "pydna GenBank file" — alternate extension), `<output_basename>*.fasta` (`Native`, "pydna FASTA"), `<output_basename>*.csv` (`Tabular`, "pydna table"), `*.log`; probe via Python on PATH with `import pydna` check (returns `ok = true` with warning when import fails so non-standard installs aren't blocked); version range `5.0.0..7.0.0` (pydna 5.x pairs with Biopython 1.8x; pydna 6.x ships in 2024 with the contemporary Gibson / Golden-Gate assembly improvements); `bio.pydna.design` ribbon capability |
| **Jalview**    | the Barton group's Java alignment viewer (GPL-3.0); the reference alignment viewer in molecular biology labs since the 2000s — multiple-sequence alignment viewing + editing, conservation / consensus / occupancy plots, structural overlays via Jmol / Chimera links, per-column annotations, the canonical interactive front-end for the MSA outputs every Phase 18 / 18.5 / 18.6 / 18.7 aligner (BWA / minimap2 / MAFFT / MUSCLE / Clustal Omega / T-Coffee) emits; ships a **headless mode** (`-nodisplay`) for batch image / format conversion — the user feeds an alignment in, picks an output format (PNG image, HTML report, SVG vector graphic, FASTA / Clustal alignment re-export), and Jalview writes the requested artifact without opening its GUI; **JAR-distributed** — no `jalview` launcher binary on PATH; the user supplies the absolute path to the JAR via `[bio.jalview].jar` in `case.toml`; single-binary subprocess shape (sister to Phase 33 j5 / Cello) with `java -jar <jar> -nodisplay -open <input> -<output_format> <basename>.<ext> [extras...]`; knobs `jar` (absolute path to the Jalview jar; required) / `input` (alignment input — `.fa` / `.aln` / `.clustal` / `.stockholm` and friends Jalview reads natively; required) / `output_basename` / `output_format` (default `"png"` — selectable from `"png"` / `"html"` / `"svg"` / `"fasta"` / `"clustal"`) / `extra_args`; `prepare()` derives the output extension from `output_format` (png → .png, html → .html, svg → .svg, fasta → .fasta, clustal → .aln, default → use the format string itself as extension); collects `<output_basename>*.png` (`Native`, "Jalview alignment image"), `<output_basename>*.svg` (`Native`, "Jalview SVG"), `<output_basename>*.html` (`Native`, "Jalview HTML"), `<output_basename>*.fasta` (`Native`, "Jalview FASTA"), `<output_basename>*.aln` (`Tabular`, "Jalview alignment"), `*.log`; probe via `find_on_path(&["java"])` (Jalview's version comes from the jar itself, not from `java`, so we surface no version here — the user pins the Jalview release implicitly by the jar they point at; same shape as Phase 33 j5 / Cello); version range `2.11.0..3.0.0` (Jalview 2.11 (2022) is the modern stable line shipping the contemporary headless improvements); `bio.jalview.view` ribbon capability |

### Web visualization (Phase 42)

First modern web 3D molecular visualization domain to ship in
Valenx — Mol* / NGL Viewer span the web-visualization tradeoff
space from the canonical PDBe / RCSB modern viewer that powers
the structural-biology web (Mol*, the EMBL-EBI / RCSB-led MIT-
licensed WebGL toolkit) to the Rose lab's WebGL framework that
predated Mol* and still powers a large fraction of the Jupyter-
friendly notebook visualization ecosystem (NGL Viewer, MIT). Both
are JavaScript browser libraries wrapped via their official
Python bindings (`molstar` / `nglview`) so they slot into the
existing Python-script subprocess pattern (sister to Phase 17
Biopython, Phase 19.5 Scanpy, Phase 33 pySBOL, Phase 41 pydna).

| Adapter      | Capability |
|---|---|
| **Mol***       | EMBL-EBI / RCSB-led modern WebGL molecular viewer (MIT); de-facto modern molecular viewer embedded in the PDB / PDBe / AlphaFold DB / ESM Atlas web properties since the late 2010s; wrapped via the `molstar` Python binding so it slots into the existing Python-script subprocess pattern (sister to Phase 17 Biopython / Phase 19.5 Scanpy / Phase 33 pySBOL / Phase 41 pydna); knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_structure` (`Option<PathBuf>` — optional `.pdb` / `.cif` / `.mmcif` structure file the script can render; `None` when the script fetches from the PDB / generates the structure inline) / `output_basename`; `prepare()` enforces `.py`, stages script + optional input_structure, writes `valenx_params.json` with `output_basename` always plus `input_structure` (staged filename) only when set — key omitted entirely when `None` rather than emitted as `null`, matching the hand-rolled JSON convention the rest of the bio adapters use; collects `<output_basename>*.html` (`Native`, "Mol* viewer HTML"), `<output_basename>*.molj` (`Native`, "Mol* state file" — the JSON state format that captures the entire viewer state for reproducible replay), `<output_basename>*.png` (`Native`, "Mol* rendered image"), `*.log`; probe via Python on PATH then `<python> -c "import molstar"` — on import failure surface as a `ProbeReport.warnings` entry, not error — sister to the Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6 HOOMD-blue / Phase 5.7 MDTraj / Phase 41 pydna probe convention; version range `3.0.0..5.0.0`; `bio.molstar.view` ribbon capability |
| **NGL Viewer** | Rose lab's high-performance WebGL framework for molecular visualization (MIT); predated Mol* and still powers a large fraction of the Jupyter-friendly notebook visualization ecosystem via its `nglview` Python binding; wrapped via the `nglview` Python binding so it slots into the existing Python-script subprocess pattern (sister to Mol*); knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_structure` (`Option<PathBuf>` — optional `.pdb` / `.cif` / `.mmcif` structure file the script can render; `None` when the script fetches from the PDB / generates the structure inline) / `output_basename`; `prepare()` mirrors Mol* shape exactly — enforces `.py`, stages script + optional input_structure, writes `valenx_params.json` with the same hand-rolled JSON convention (key omitted when `None`); collects `<output_basename>*.html` (`Native`, "NGL viewer HTML"), `<output_basename>*.png` (`Native`, "NGL rendered image"), `<output_basename>*.json` (`Tabular`, "NGL state JSON"), `*.log`; probe via Python on PATH then `<python> -c "import nglview"` — on import failure surface as a `ProbeReport.warnings` entry, not error — sister to the Mol* probe convention; version range `3.0.0..5.0.0`; `bio.ngl.view` ribbon capability |

### mRNA design (Phase 43)

First mRNA / vaccine therapeutic design domain to ship in Valenx
— DNA Chisel / LinearDesign / iCodon span the codon-optimization
+ joint-design tradeoff space from the Edinburgh Genome Foundry's
constraint-driven Python codon optimizer (DNA Chisel, MIT —
codon optimization, restriction-site avoidance, repeat scanning,
GC-content tuning, arbitrary user constraints) to Baidu Research's
joint codon + secondary-structure mRNA design tool (LinearDesign,
Apache-2.0 — the modern mRNA-vaccine design workhorse since the
2021 _Nature_ paper, joint CAI + MFE optimization tunable via
`lambda_param`) and the Vejnar lab's R-based codon-level stability
predictor (iCodon, GPL-3.0 — per-position codon contributions to
mRNA half-life). DNA Chisel + iCodon ride the established
Python-script / Rscript subprocess patterns (sister to Phase 17
Biopython / Phase 19.6 Seurat); LinearDesign rides the single-
binary CLI pattern sister to Phase 18 BWA / Phase 32.5 Smoldyn.

| Adapter      | Capability |
|---|---|
| **DNA Chisel**   | Edinburgh Genome Foundry constraint-driven codon optimization (MIT); de-facto Python choice for end-to-end synthetic-gene design pipelines feeding into Phase 33 j5 assembly + Phase 41 pydna cloning workflows; Python-script subprocess shape sister to Phase 17 Biopython / Phase 19.5 Scanpy / Phase 33 pySBOL / Phase 41 pydna / Phase 42 Mol* / NGL Viewer; knobs `script` (`.py` enforced) / `python` (default `"python3"`) / `input_fasta` (`Option<PathBuf>` — optional starting `.fa` / `.fasta` FASTA) / `output_basename`; `prepare()` enforces `.py`, routes script + optional input_fasta through `confined_join` to stage them safely in the workdir, writes `valenx_params.json` with `output_basename` always plus `input_fasta` (staged filename) only when set — key omitted entirely when `None` rather than emitted as `null`, matching the hand-rolled JSON convention the rest of the bio adapters use; collects `<output_basename>*.fasta` (`Native`, "DNA Chisel optimized FASTA"), `<output_basename>*.gb` / `.genbank` (`Native`, "DNA Chisel GenBank"), `<output_basename>*.json` (`Tabular`, "DNA Chisel constraint report" — the canonical machine-readable per-constraint-pass / per-constraint-fail report DNA Chisel emits for downstream automation), `<output_basename>*.png` (`Native`, "DNA Chisel plot"), `*.log`; probe via Python on PATH then `<python> -c "import dnachisel"` — on import failure surface as a `ProbeReport.warnings` entry, not error — sister to the Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6 HOOMD-blue / Phase 5.7 MDTraj / Phase 41 pydna / Phase 42 Mol* / NGL probe convention; version range `3.0.0..4.0.0`; `bio.dnachisel.optimize` ribbon capability |
| **LinearDesign** | Baidu Research joint codon + secondary-structure mRNA design (Apache-2.0); landed as the modern mRNA-vaccine design workhorse following the 2021 _Nature_ paper that demonstrated dramatic stability / expression gains for mRNA vaccines designed under the joint CAI + MFE objective; single-binary CLI subprocess shape sister to Phase 18 BWA / Phase 32.5 Smoldyn / Phase 5 GROMACS with `lineardesign --aa <protein> --lambda <lambda_param> --codon_usage <codon_usage> --output_basename <basename> [extras...]`; knobs `protein` (path to protein FASTA; required — read in place, no staging) / `output_basename` / `lambda_param` (`f64`, finite and ≥ 0.0; default 1.0; the Rust field is `lambda_param` because `lambda` is a Rust reserved keyword — the CLI emits `--lambda <value>` regardless; tunable Lagrangian tradeoff between codon-adaptation-index and predicted mRNA secondary-structure stability — `0.0` = pure MFE-optimal, large = pure CAI-optimal, intermediate values like default `1.0` hit the joint sweet spot demonstrated in the paper) / `codon_usage` (default `"human"` — selectable from the LinearDesign-shipped set: `"human"` / `"mouse"` / `"yeast"` / `"ecoli"` / etc.) / `extra_args`; `prepare()` validates `lambda_param` is finite and ≥ 0.0 (returns `InvalidCase` when negative or NaN — LinearDesign's optimizer would either crash or silently collapse to MFE-only design on invalid input), resolves `protein` against the case directory when relative, validates the file exists on disk; collects `<output_basename>*.fasta` (`Native`, "LinearDesign optimized mRNA"), `<output_basename>*.txt` (`Tabular`, "LinearDesign report" — the canonical per-design summary report LinearDesign writes alongside the FASTA output), `*.log`; probe via `find_on_path(&["lineardesign"])` — when the `lineardesign` binary isn't found but Python is on PATH the probe surfaces a targeted `"clone https://github.com/LinearDesignSoftware/LinearDesign and add the bin directory to PATH"` warning so users see the install hint immediately; version range `1.0.0..2.0.0`; `bio.lineardesign.design` ribbon capability |
| **iCodon**       | Vejnar lab's codon-level mRNA stability prediction (GPL-3.0); predicts per-position codon contributions to mRNA half-life given a target organism — the canonical readout for "given this mRNA sequence, which codons are dragging down stability and where would a codon swap help"; the canonical R-based mRNA stability predictor; ships as a `devtools::install_github('santiago1234/iCodon')` R package; **Rscript subprocess pattern** sister to Phase 19.6 Seurat — the user supplies an `.R` script that loads `library(iCodon)` and reads `valenx_params.json` for the parsed knobs via `jsonlite::fromJSON`; knobs `script` (`.R` enforced) / `rscript` (default `"Rscript"`) / `input_fasta` (`Option<PathBuf>` — optional input mRNA FASTA the script can score; `None` when the script generates the sequence inline) / `output_basename`; `prepare()` enforces the `.R` extension, routes script + optional input_fasta through `confined_join` to stage them safely in the workdir, writes `valenx_params.json` with `output_basename` always plus `input_fasta` (staged filename) only when set — key omitted entirely when `None` rather than emitted as `null` (same hand-rolled JSON shape as DNA Chisel + every other Phase 19.6+ adapter that takes an optional input file); collects `<output_basename>*.csv` / `*.tsv` (`Tabular`, "iCodon stability table"), `<output_basename>*.rds` (`Native`, "iCodon R object (RDS)"), `<output_basename>*.png` (`Native`, "iCodon plot"), `*.log`; probe via `find_on_path(&["Rscript"])` — does not attempt to confirm iCodon itself is installed because that would require running R, an expensive multi-second startup at probe time (same shape as Phase 19.6 Seurat); the `ToolNotInstalled` install hint mentions the canonical `devtools::install_github('santiago1234/iCodon')` install path; version range `1.0.0..2.0.0`; `bio.icodon.predict` ribbon capability |

### Pharmacokinetics + RNA tertiary (Phase 45)

Opens **two new domains** in Valenx with two single-canonical-
adapter beachheads — PK/PD pharmacokinetics (PK-Sim) and RNA
tertiary 3D structure prediction (SimRNA). The two domains are
unrelated biologically but Phase 45 ships them together because
each is a single canonical adapter and each opens a never-before-
seen domain in Valenx. PK-Sim is the **first PK/PD modeling
category in Valenx** — distinct from the Phase 32 / 32.5 systems-
biology cellular-scale ODE / spatial-stochastic modeling, PK-Sim
covers whole-body absorption / distribution / metabolism /
excretion (ADME) drug pharmacokinetics. SimRNA is the **first RNA
tertiary 3D structure prediction category in Valenx** — distinct
from the Phase 28 + 44.5 RNA secondary (2D base-pairing) folders,
SimRNA predicts the full 3D Cartesian backbone via coarse-grained
Monte Carlo replica-exchange sampling. Both adapters ride the
single-binary subprocess shape sister to Phase 18 BWA / Phase
32.5 Smoldyn / Phase 5 GROMACS.

| Adapter      | Capability |
|---|---|
| **PK-Sim**     | Open Systems Pharmacology suite's physiologically-based PK (PBPK) simulator (GPL-2.0) — de-facto open-source PBPK modeling tool, descended from the Bayer internal pharmacokinetic simulator opened to the community via the Open Systems Pharmacology Initiative; models whole-body drug absorption / distribution / metabolism / excretion (ADME) using a physiologically-grounded compartmental representation (every major organ — liver, kidney, lung, gut, adipose, muscle, brain — is a compartment with its own blood flow, volume, partition coefficient, metabolic capacity); user supplies a `.pksim5` project file (XML-based, authored in PK-Sim GUI or programmatically through the OSP Python API) describing drug + dosing protocol + simulated population; single-binary subprocess shape sister to Phase 18 BWA with `pksim --project <project> --output <output_basename> [extras...]`; knobs `project` (`.pksim5` project file; required — read in place from the case directory, no staging) / `output_basename` / `extra_args`; `prepare()` resolves `project` against the case directory when relative, validates it exists on disk; collects `<output_basename>*.csv` (`Tabular`, "PK-Sim simulation results"), `<output_basename>*.json` (`Tabular`, "PK-Sim metadata"), `*.log`; probe via `find_on_path(&["pksim", "PKSim.CLI"])` (modern OSP distribution ships both the generic `pksim` launcher and the .NET-style `PKSim.CLI` invoker for Windows installs); version range `11.0.0..13.0.0`; `bio.pksim.simulate` ribbon capability |
| **SimRNA**     | Bujnicki group's coarse-grained Monte Carlo RNA tertiary-structure predictor (GPL-3.0) — predicts the full 3D Cartesian backbone of an RNA from its sequence (NOT just the 2D base-pairing pattern that the Phase 28 ViennaRNA / RNAstructure / NUPACK and Phase 44.5 mfold / EternaFold / LinearFold folders predict, but the actual 3D structure with per-residue Cartesian coordinates suitable for ChimeraX / PyMOL visualisation, structural alignment against PDB, MD simulation via Phase 5 GROMACS); model represents each nucleotide as five coarse-grained beads (one per phosphate / sugar / base ring) and runs replica-exchange Monte Carlo over the resulting reduced-coordinate energy landscape, sampling the conformational ensemble around the predicted minimum; single-binary subprocess shape sister to PK-Sim / Phase 18 BWA / Phase 32.5 Smoldyn with `SimRNA -c <config> -s <sequence> -o <output_basename> -R <n_replicas> [extras...]`; knobs `config` (SimRNA configuration file) / `sequence` (`.seq` SimRNA-format sequence file) / `output_basename` / `n_replicas` (`u32`, ≥ 1, default 1 — replica exchange helps escape local minima in the rugged RNA tertiary landscape) / `extra_args`; `prepare()` resolves both `config` and `sequence` against the case directory (read in place, no staging), validates each file exists on disk; collects `<output_basename>*.pdb` (`Native`, "SimRNA tertiary structure" — predicted 3D Cartesian backbone in PDB format, lifted by the existing Phase 17 PDB reader without RNA-3D-specific code), `<output_basename>*.trafl` (`Native`, "SimRNA trajectory"), `<output_basename>*.txt` (`Tabular`, "SimRNA energy log"), `*.log`; probe via `find_on_path(&["SimRNA", "simrna"])`; version range `3.20.0..4.0.0`; `bio.simrna.fold` ribbon capability |

That's **🎯 141 live adapters — bio ecosystem complete + Phase 43 mRNA design + Phase 44.5 RNA folding expansion + Phase 35.5 base + prime editing design + Phase 35.6 edit-outcome prediction + Phase 45 pharmacokinetics + RNA tertiary structure**: every physics-domain phase 1-9 plus the 43-phase biology / biotech / chemistry expansion (Phase 5.5 MD analysis expansion, Phase 5.6 bio MD engines, Phase 5.7 MDTraj, Phase 17 biology, Phase 17.5 structure prediction, Phase 17.7 structure tools, Phase 18 sequence alignment, Phase 18.5 aligners expansion, Phase 18.6 RNA-seq alignment, Phase 18.7 alignment-toolkit expansion, Phase 19 variant calling, Phase 19.5 single-cell, Phase 19.6 single-cell expansion, Phase 20 transcript quantification, Phase 22 workflow managers, Phase 22.5 workflow expansion — sister-expansion of Phase 22, Phase 23 molecular viewers, Phase 24 cheminformatics, Phase 25 quantum chemistry, Phase 27 protein design, Phase 27.5 protein-design expansion, Phase 27.6 EvolutionaryScale models, Phase 28 RNA structure, Phase 29 population genetics, Phase 30 phylogenetics, Phase 30.5 Bayesian phylogenetics, Phase 31 read simulators, Phase 32 systems biology, Phase 32.5 spatial stochastic — sister-expansion of Phase 32, Phase 33 synthetic biology, Phase 34 docking, Phase 35 CRISPR design, Phase 35.5 base + prime editing design — sister-expansion of Phase 35, Phase 35.6 edit-outcome prediction — sister-expansion of Phase 35 / 35.5, Phase 36 cryo-EM, Phase 38 Rosetta family, Phase 39 DNA structural geometry, Phase 40 microscopy — first bioimage analysis category in Valenx, Phase 41 sequence editors — first plasmid-design / alignment-viewer category in Valenx, Phase 42 web visualization — first modern web 3D molecular visualization category in Valenx, Phase 43 mRNA design — first mRNA / vaccine therapeutic design category in Valenx, Phase 44.5 RNA folding expansion — sister-expansion of Phase 28, Phase 45 pharmacokinetics + RNA tertiary — first PK/PD modeling category + first RNA tertiary 3D structure prediction category in Valenx). The bio surface now spans alignment / base + prime editing / cheminformatics / CRISPR / cryo-EM / DNA geometry / docking / edit-outcome prediction / MD analysis / MD engines / microscopy / mRNA design / pharmacokinetics / phylogenetics / population genetics / protein design / quantum chemistry / RNA structure (2D + 3D) / sequence editors / sequence read simulators / single-cell / spatial stochastic / structure prediction / structure search / synthetic biology / systems biology / variant calling / viewers (desktop + web) / web visualization / workflow managers. The bio-adapter count totals **123 bio adapters across 43 biology / biotech / chemistry phases**; the headline live-adapter number is **141 fully live**.

## Scaffolds (probe-only, no real prepare/run yet)

These return `not_implemented` from `prepare()` but have working
`probe()` against the underlying tool — they exist so the registry
shows the right status when the user has the tool installed:

> *(none at present — the prior `valenx-adapter-occt` stub was
> rewritten as a pythonocc-core subprocess wrapper, dropping the
> last remaining scaffold.)*

## Documentation-only (not yet started)

- **Phases 11-16** — HPC / cluster execution, optimization,
  ML surrogates, plugin marketplace, enterprise features,
  stewardship. All have planning docs; none have code.
- **Most of Phase 10 polish** — i18n pipeline, signed installers,
  theme snapshot tests, first-run wizard. The workflow-loop
  pieces of Phase 10 landed early during Phase 1; the bigger
  ship-readiness pieces are still planned.

## What you can build with this today

A complete (pre-alpha) CFD pipeline:

```
gmsh adapter → mesh.canonical.json
  ↓
OpenFOAM adapter (simpleFoam / pimpleFoam / icoFoam / rhoSimpleFoam)
  ↓
foamToVTK → cavity_*.vtu
  ↓
collect() → Results.fields
  ↓
viewport: mesh + colour-bar + time slider
```

…with the same workflow loop available for the other live adapters
(FEA, heat, chemistry, MD, EM, battery, multibody, coupling, plus the
3 MD-analysis-expansion adapters from Phase 5.5, 3 bio MD-engine
adapters from Phase 5.6, 1 MDTraj adapter from Phase 5.7, 7 biology /
biotech adapters from Phase 17, 4 structure-prediction adapters from
Phase 17.5, 3 structure-tools-expansion adapters from Phase 17.7,
6 sequence-alignment adapters from Phase 18, 3
aligners-expansion adapters from Phase 18.5, 2 RNA-seq alignment
adapters from Phase 18.6, 3 alignment-toolkit-expansion adapters
from Phase 18.7, 3 variant-calling adapters from Phase 19,
2 single-cell genomics adapters from Phase 19.5, 2
single-cell-expansion adapters from Phase 19.6, 2
transcript-quantification adapters from Phase 20, 2
workflow-manager adapters from Phase 22, 3 workflow-expansion
adapters from Phase 22.5, 3 molecular-viewer
adapters from Phase 23, 3 cheminformatics-expansion adapters from
Phase 24, 3 quantum-chemistry adapters from Phase 25, 2
protein-design adapters from Phase 27, 3 protein-design-expansion
adapters from Phase 27.5, 2 EvolutionaryScale-model adapters from
Phase 27.6, 3 RNA-structure adapters from Phase 28, 3
population-genetics adapters from Phase 29, 3 phylogenetics
adapters from Phase 30, 2 Bayesian-phylogenetics adapters from
Phase 30.5, 3 read-simulator adapters from Phase 31, 3
systems-biology adapters from Phase 32, 2 spatial-stochastic
adapters from Phase 32.5, 3 synthetic-biology adapters from
Phase 33, 2 molecular-docking adapters from Phase 34,
3 CRISPR-design adapters from Phase 35, 4 base + prime editing
design adapters from Phase 35.5, 4 edit-outcome-prediction
adapters from Phase 35.6, 3 cryo-EM adapters from Phase 36, 2
Rosetta-family adapters from Phase 38, 3 DNA-structural-geometry
adapters from Phase 39, 3 microscopy adapters from Phase 40, 2
sequence-editor adapters from Phase 41, 2 web-visualization
adapters from Phase 42, 3 mRNA-design adapters from Phase 43,
3 RNA-folding-expansion adapters from Phase 44.5, and 2
pharmacokinetics + RNA-tertiary adapters from Phase 45) —
the results-rendering
pipeline is OpenFOAM-specific today;
other adapters populate `Results.artifacts` only and need their own
VTK / .frd / .vtu parsers wired into `collect()` to get the full
visual loop.

## Phase 10 release substrate

The pieces that gate cutting a public alpha:

- **Crash reporter** (`valenx-crash-reporter`) — panic hook
  installed at app startup; reports land sanitised in
  `<state_dir>/crashes/`. Network egress gated on Settings →
  Privacy opt-in (default OFF).
- **First-launch wizard** (`valenx-first-run` + `valenx-app::first_run`) —
  auto-opens on a fresh install with adapter probe results +
  per-OS install hints; persists `<state_dir>/first-run.json`.
- **String catalogue** (`valenx-i18n`) — `.ftl`-shaped key=value
  catalogue with `{ $name }` placeholder substitution. en-US
  baseline ships ~75 keys; About + first-run + Settings panels
  wired through. Pseudo-locale `to_pseudo()` for dev builds.
- **Accessibility gate** (`valenx-a11y` + `valenx-design-tokens`
  CI test) — WCAG 2.1 contrast ratios + AA / AAA classifiers;
  every documented foreground/background pair runs through the
  audit on every CI build.
- **Signed-installer pipeline** (`.github/workflows/release.yml` +
  `RELEASING.md`) — deb / rpm / .app / .msi build matrix that
  signs when `APPLE_ID` / `WINDOWS_CERT` secrets are set, falls
  through to unsigned artefacts otherwise. Cert procurement is
  optional for v0.1.0-alpha (rustup / uv / ripgrep convention).

## Headless CLI tooling

Twelve CLIs ship for users who want to script Valenx without the GUI:

- `valenx-init` — scaffold a new project from a template.
- `valenx-validate` — structural pre-flight on a project bundle
  (manifest, tools.lock, every case in `[cases].order`). Text or
  JSON output, exit-code-driven for CI.
- `valenx-mesh-info` — inspect a `mesh.canonical.json` (text or
  JSON output). `--check max-skew=0.9 --check inverted=0` style
  threshold gates for CI.
- `valenx-audit verify` and `valenx-audit tail` — offline
  audit-log integrity check + recent-activity tail. `tail`
  supports `--since <ISO-8601>` for time-window filtering.
- `valenx-results` — inspect the `results.json` sidecar a finished
  run leaves on disk (fields, scalars, artifacts, provenance).
- `valenx-report` — write a self-contained HTML report, a
  GitHub-flavoured Markdown summary, and/or a flat scalar history
  CSV from a `results.json`. CI-friendly: refuses to no-op without
  `--html` / `--markdown` / `--csv`.

### Biology CLIs (Phase 17)

- `valenx-fasta` — inspect / validate / extract sequences from
  FASTA files. Text + JSON output, stdin via `-`.
- `valenx-pdb-info` — structural summary (chain / residue / atom
  counts + element tally) from a PDB file.
- `valenx-blast` — thin BLAST+ wrapper that auto-detects alphabet
  from the query and routes to blastp / blastn.

### Sequence alignment CLIs (Phase 18)

- `valenx-fastq` — inspect / validate FASTQ files. Text + JSON
  output, stdin via `-`.
- `valenx-sam-info` — alignment summary (record count, mapped /
  unmapped tally, reference list, average MAPQ) from a SAM file.

### Variant calling CLIs (Phase 19)

- `valenx-vcf-info` — VCF summary (header-line count, sample count,
  total records, PASS / FAIL split, no-ALT count) from a VCF file or
  stdin via `-`. Text + JSON output.

All twelve offer JSON output (or self-contained HTML for `valenx-report`)
for downstream tooling and live in `crates/<crate>/src/bin/` so
`cargo run --bin <name>` works without extra setup.

## Test posture

The scoped `scripts/qa.sh` (or `scripts/qa.ps1`) ships **10,000+ passing
tests** with **zero clippy warnings and zero rustdoc warnings** across
the workspace. A blanket `cargo test --workspace` is intentionally
forbidden in this repo — see `docs/QA.md` for the rationale.
The harness covers:

- Every adapter's case-input parser
- Every adapter's deck / dict / SIF / .inp / .py emitter
- The shared subprocess runner + cancellation token
- The canonical Mesh / Field / Results types
- The `.vtu` parser + canonical converter + colour ramp
- The Netgen `.vol` parser (Tet / Pyr / Prism / Hex; 2D Tri / Quad)
- The app's run pipeline (spawn → progress → finished → collected)
- All workflow actions (run / prepare / run-from-prepared /
  open-workdir / field selection / time stepping)
- Adapter resolution from solver string conventions
- Default-state regressions on every new ValenxApp field
- Subprocess integration tests for every CLI (mesh-info / audit /
  validate / results) — assert exit codes + stdout shapes against
  the real compiled binary

Plus a per-adapter probe smoke test on every CI run, and the
`tests/fixtures/minimal.valenx` fixture round-trip across six
demo cases (CFD steady + transient, FEA cantilever, heat cube,
gmsh box mesh, Netgen CSG cylinder).

## What's missing for "first usable release"

- **Probe at point** — click on the mesh, see the field value at
  that node / cell. Wireframe colouring + filled triangles + the
  Quality panel ship today; per-pick numerical readouts don't.
- **Live `tail -f` of cluster log files during the run** —
  `read_slurm_log_tail` reads after-the-fact; ssh-tail-style
  streaming during the run is a follow-up.
- **End-to-end showcase against real solvers** — the per-adapter
  tests exercise prepare + collect classification graceful paths;
  a CI cluster fixture that actually runs SU2 / OpenFOAM /
  GROMACS to completion is its own infra build-out.

Realistic effort to close those: **~1-2 months of focused work**
for a small team. See individual phase docs in `docs/src/phases/`
for the full per-domain backlog.
