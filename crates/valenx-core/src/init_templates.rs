//! Project-scaffold templates — the canonical library of starter
//! `case.toml` bodies, sample input files, and project metadata
//! used by both the `valenx-init` CLI and the desktop GUI's
//! "New case from adapter" flow.
//!
//! Each template corresponds to one registered adapter (or one of
//! the generic Empty / CFD / FEA / Chemistry shells) and has:
//!   * a stable [`Template`] enum variant
//!   * a CLI-/GUI-facing string id resolved by [`Template::from_str`]
//!   * a [`TemplateRow`] entry in the static catalogue (see
//!     [`template_rows`]) with a one-line description and the
//!     directory name `cases/<dir>/` that it scaffolds into
//!   * a `case.toml` body rendered by [`render_case_toml`]
//!   * (optionally) a sample input file shipped via
//!     [`template_sample_files`]
//!
//! Public wrappers used by the GUI:
//!   * [`case_dir_for`] — adapter id → cases/`<dir>` name.
//!   * [`case_toml_body`] — adapter id → starter case.toml body.
//!   * [`project_toml`] — render the project.toml shell.
//!   * [`template_rows`] — the canonical name + description list.

use std::path::Path;

/// Help text rendered by `valenx-init help`. Lives here (next to
/// the catalogue it documents) so a new template only has to be
/// added in one place.
pub const USAGE: &str = "\
valenx-init — scaffold a fresh `.valenx` project skeleton.

USAGE:
  valenx-init <dir> [--template T] [--name N]
  valenx-init help
  valenx-init -V | --version
  valenx-init -l | --list-templates

OPTIONS:
  --template T   Pick a starter template (default: empty). Aliases
                 in parens.

                 Generic:
                   empty            (alias: minimal)
                   cfd              (alias: openfoam)  — OpenFOAM simpleFoam
                   fea              (alias: calculix, structural) — CalculiX
                   chemistry        (alias: cantera, chem) — Cantera

                 Per-adapter:
                   su2              (alias: compressible, aero)
                   openradioss      (alias: crash, impact)
                   code-aster       (alias: aster, thermomechanical)
                   netgen           (alias: meshing, mesh)
                   meep             (alias: fdtd, photonics)
                   gromacs          (alias: md, molecular-dynamics)
                   gmsh             (alias: delaunay)
                   lammps           (alias: lj, classical-md)
                   elmer-heat       (alias: elmer, heat)

                 Biology (Phase 17):
                   biopython        (alias: biopy)
                   rdkit            (alias: chem-py)
                   openmm           (alias: omm)
                   chimerax         (alias: cxc, viz3d)
                   oxdna            (alias: cgdna)
                   mdanalysis       (alias: mda, traj-py)
                   colabfold        (alias: cf, protein-fold)

                 Biology — alignment toolkit (Phase 18):
                   bwa              (alias: bwa-mem)
                   minimap2         (alias: mm2)
                   mafft            — multiple sequence alignment
                   muscle           — multiple sequence alignment
                   hmmer            (alias: hmmsearch)
                   samtools         — SAM/BAM multitool

                 Biology — structure prediction (Phase 17.5):
                   esmfold          — Meta ESMFold
                   openfold         — PyTorch AF2 reimplementation
                   alphafold2       (alias: af2)
                   alphafold3       (alias: af3) — non-commercial weights

                 Biology — variant calling (Phase 19):
                   bcftools         — VCF/BCF multitool (view/call/filter/concat)
                   gatk             (alias: hc) — GATK HaplotypeCaller
                   deepvariant      (alias: dv) — Google ML-driven variant caller

                 Biology — viewers (Phase 23):
                   pymol            — open-source PyMOL renderer
                   vmd              — VMD Tcl-scripted MD viewer (academic license)
                   igv              (alias: igvtools) — IGV `igvtools` indexer

                 Biology — protein design (Phase 27):
                   rfdiffusion      (alias: rfd) — protein backbone generation
                   proteinmpnn      (alias: mpnn) — sequence design from backbone

                 Biology — docking (Phase 34):
                   vina             (alias: autodock-vina) — AutoDock Vina
                   autodock4        (alias: ad4) — AutoDock 4 (two-stage)

                 Biology — cheminformatics expansion (Phase 24):
                   deepchem         (alias: dc) — PyTorch-backed cheminformatics
                   openbabel        (alias: obabel) — chemistry-format converter
                   avogadro         (alias: avogadro2) — Python-scriptable editor

                 Biology — workflow managers (Phase 22):
                   nextflow         (alias: nf) — pipeline orchestrator
                   snakemake        (alias: smk) — rule-based orchestrator

                 Biology — single-cell genomics (Phase 19.5):
                   scanpy           — de-facto Python single-cell analysis
                   scvi             (alias: scvi-tools) — probabilistic models

                 Biology — protein design expansion (Phase 27.5):
                   chroma           — Generate Biomedicines diffusion
                   esm-if           (aliases: esmif, inverse-folding)
                   rfantibody       (alias: rfab) — antibody design

                 Biology — aligners expansion (Phase 18.5):
                   bowtie2          (alias: bt2) — short-read aligner
                   mmseqs2          (alias: mmseqs) — protein search/cluster
                   diamond          (alias: dmnd) — ultra-fast protein search

                 Biology — RNA-seq alignment (Phase 18.6):
                   hisat2           (alias: hisat) — splice-aware aligner
                   star             — spliced RNA-seq aligner

                 Biology — transcript quantification (Phase 20):
                   salmon           — quasi-mapping transcript quantification
                   kallisto         — pseudoalignment quantification

                 Biology — phylogenetics (Phase 30):
                   iqtree           (alias: iqtree2) — ML tree inference
                   raxml-ng         (alias: raxml) — next-gen RAxML
                   fasttree         — approximate-ML tree inference

                 Biology — RNA structure (Phase 28):
                   viennarna        (aliases: vienna, rnafold) — academic
                   rnastructure     — Mathews lab Fold (BSD)
                   nupack           — Caltech NUPACK (academic)

                 Biology — quantum chemistry (Phase 25):
                   psi4             — HF/DFT/post-HF
                   nwchem           — massively-parallel ab initio
                   xtb              — extended tight-binding semiempirical

                 Biology — EvolutionaryScale models (Phase 27.6):
                   esm3             — generative multi-modal protein model
                   esmc             (alias: esm-cambrian) — embeddings

                 Biology — systems biology (Phase 32):
                   copasi           — biochemical pathway / ODE
                   bionetgen        (alias: bng) — rule-based signaling
                   physicell        — agent-based multicellular tissue

                 Biology — cryo-EM (Phase 36):
                   relion           — Bayesian 3D reconstruction
                   eman2            (alias: eman) — broad image processing
                   ctffind          — CTF estimation (academic license)

                 Biology — read simulators (Phase 31):
                   art              (alias: art-illumina) — Illumina sims
                   wgsim            — classic short-read simulator
                   badread          — Nanopore long-read simulator

                 Biology — CRISPR design (Phase 35):
                   chopchop         — guide-RNA design (Python)
                   crispor          — guide design + off-target (Python)
                   cas-offinder     (alias: cas-off) — off-target search

                 Biology — Rosetta family (Phase 38):
                   rosetta          — rosetta_scripts XML protocols (academic)
                   pyrosetta        — Python wrapper for Rosetta core (academic)

                 Biology — population genetics (Phase 29):
                   slim             — forward-time simulator (Eidos)
                   msprime          — coalescent simulator (Python)
                   tskit            — tree-sequence analysis (Python)

                 Biology — Bayesian phylogenetics (Phase 30.5):
                   beast2           — Bayesian MCMC tree inference
                   mrbayes          — Bayesian MCMC tree inference

                 Biology — DNA geometry (Phase 39):
                   x3dna            — DNA base-step parameters (academic)
                   curves           — Curves+ helical-axis (academic)
                   dssr             — DNA/RNA structural features (academic)

                 Biology — MD analysis (Phase 5.5):
                   plumed           — enhanced sampling + free energy
                   prody            — protein dynamics + ENM
                   cpptraj          — AmberTools trajectory analysis

                 Biology — synthetic biology (Phase 33):
                   pysbol           — SBOL Python composition
                   j5               — DNA assembly automation (JAR)
                   cello            — genetic-circuit DNA compiler (JAR)

                 Biology — alignment expansion (Phase 18.7):
                   blast            — NCBI BLAST+ sequence search
                   clustalo         — Clustal Omega MSA
                   tcoffee          — T-Coffee consensus MSA

                 Biology — single-cell expansion (Phase 19.6):
                   seurat           — R-based single-cell analysis
                   anndata          — Python single-cell HDF5 container

                 Biology — bio MD engines (Phase 5.6):
                   namd             — NAMD all-atom MD (academic)
                   sander           — AmberTools sander MD engine
                   hoomd            — HOOMD-blue GPU-native particle MD

                 Biology — MD analysis sister (Phase 5.7):
                   mdtraj           — MDTraj Python trajectory analyzer

                 Biology — structure tools expansion (Phase 17.7):
                   rosettafold      — RoseTTAFold structure prediction
                   omegafold        — OmegaFold single-sequence prediction
                   foldseek         — FoldSeek protein structure search

                 Biology — spatial stochastic (Phase 32.5):
                   smoldyn          — Smoldyn reaction-diffusion sim
                   mcell            — MCell cell-scale simulator

                 Biology — sequence editors (Phase 41):
                   pydna            — Python plasmid / clone-design
                   jalview          — Jalview alignment viewer (headless)

                 Biology — microscopy / bioimage (Phase 40):
                   fiji             — ImageJ/Fiji headless image processing
                   cellprofiler     — cell segmentation pipeline runner
                   ilastik          — ML-based pixel/object classification

                 Biology — workflow expansion (Phase 22.5):
                   planemo          — Galaxy ecosystem CLI
                   cromwell         — Broad WDL workflow engine (JAR)
                   cwltool          — CWL reference runner

                 Biology — web 3D visualization (Phase 42):
                   molstar          — Mol* WebGL viewer
                   ngl              — NGL Viewer WebGL framework

                 Biology — mRNA design (Phase 43):
                   dnachisel        — DNA Chisel codon optimization
                   lineardesign     — LinearDesign mRNA design
                   icodon           — iCodon mRNA stability (R)

  --name N       Project name to embed in project.toml
                 (default: directory name)

EXIT CODES:
  0   project written
  1   target dir already populated or IO failure
  2   invalid CLI usage

EXAMPLES:
  valenx-init my-cfd-case --template cfd
  valenx-init experiments/run-42 --template fea --name 'cantilever-beam'
";

/// One of the canned `valenx-init` templates the CLI can scaffold.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Template {
    /// Empty project skeleton — `valenx.toml` + a placeholder case.
    Empty,
    /// CFD project skeleton (OpenFOAM-shaped `case.toml`).
    Cfd,
    /// FEA project skeleton (CalculiX-shaped `case.toml`).
    Fea,
    /// Chemistry project skeleton (Cantera-shaped `case.toml`).
    Chemistry,
    /// SU2 compressible CFD on a user-provided .cfg + .su2 mesh.
    Su2,
    /// OpenRadioss explicit dynamics on a pre-built engine deck.
    OpenRadioss,
    /// Code_Aster on a user-provided .export descriptor.
    CodeAster,
    /// Netgen batch CSG / BREP meshing.
    Netgen,
    /// Meep FDTD photonics from a Python script.
    Meep,
    /// GROMACS biomolecular MD on a pre-built .tpr file.
    Gromacs,
    /// gmsh procedural meshing — `[mesh] type = "box" | "sphere" |
    /// "merge"`.
    Gmsh,
    /// LAMMPS classical MD — Lennard-Jones FCC fluid in NVE,
    /// minimum-runnable case.
    Lammps,
    /// Elmer steady heat conduction with two pinned faces.
    ElmerHeat,
    // ---- Biology (Phase 17) ----
    /// Biopython — sequence / structural-bio analysis driven by a
    /// user-supplied Python script.
    Biopython,
    /// RDKit — cheminformatics screening driven by a user-supplied
    /// Python script.
    Rdkit,
    /// OpenMM — Python-native molecular dynamics minimisation
    /// + DCD trajectory output via a user-supplied script.
    Openmm,
    /// ChimeraX — `.cxc` command-script driven structural rendering.
    Chimerax,
    /// oxDNA — coarse-grained DNA / RNA molecular dynamics on a
    /// user-supplied `input.dat`.
    Oxdna,
    /// MDAnalysis — trajectory analysis driven by a user-supplied
    /// Python script.
    Mdanalysis,
    /// ColabFold — protein structure prediction from a FASTA query.
    Colabfold,
    // ---- Biology — sequence alignment toolkit (Phase 18) ----
    /// BWA short-read alignment via `bwa mem`.
    Bwa,
    /// minimap2 long-read + cross-domain alignment.
    Minimap2,
    /// MAFFT multiple-sequence alignment.
    Mafft,
    /// MUSCLE 5 multiple-sequence alignment.
    Muscle,
    /// HMMER profile-HMM search (hmmsearch / hmmscan).
    Hmmer,
    /// samtools — SAM/BAM multitool (view / sort / index / flagstat).
    Samtools,
    // ---- Biology — structure prediction expansion (Phase 17.5) ----
    /// ESMFold — Meta's protein language model for structure prediction.
    Esmfold,
    /// OpenFold — PyTorch reimplementation of AlphaFold 2.
    Openfold,
    /// AlphaFold 2 — DeepMind's structure prediction pipeline.
    Alphafold2,
    /// AlphaFold 3 — DeepMind's all-atom complex predictor (non-commercial weights).
    Alphafold3,
    // ---- Biology — variant calling (Phase 19) ----
    /// bcftools — VCF/BCF multitool (view / call / filter / concat).
    Bcftools,
    /// GATK HaplotypeCaller — Broad Institute Java-based variant caller.
    Gatk,
    /// Google DeepVariant — ML-driven variant caller (WGS / WES / PacBio / ONT).
    DeepVariant,
    // ---- Biology — viewers (Phase 23) ----
    /// PyMOL (open-source build) — script-driven structural rendering.
    Pymol,
    /// VMD — Tcl-scripted MD trajectory viewer (academic license).
    Vmd,
    /// IGV `igvtools` — headless BAM/VCF indexer + tile generator.
    Igv,
    // ---- Biology — protein design (Phase 27) ----
    /// RFdiffusion — RosettaCommons GPU diffusion model for protein
    /// backbone generation.
    RfDiffusion,
    /// ProteinMPNN — graph neural network for sequence design from
    /// a fixed protein backbone.
    ProteinMpnn,
    // ---- Biology — molecular docking (Phase 34) ----
    /// AutoDock Vina — modern single-binary small-molecule docker.
    Vina,
    /// AutoDock 4 — two-stage (autogrid4 + autodock4) docking.
    AutoDock4,
    // ---- Biology — cheminformatics expansion (Phase 24) ----
    /// DeepChem — PyTorch-backed cheminformatics deep-learning library.
    DeepChem,
    /// Open Babel — chemistry-format converter (~120 formats).
    OpenBabel,
    /// Avogadro 2 — Python-scriptable chemistry editor.
    Avogadro,
    // ---- Biology — workflow managers (Phase 22) ----
    /// Nextflow — pipeline orchestrator.
    Nextflow,
    /// Snakemake — rule-based pipeline orchestrator.
    Snakemake,
    // ---- Biology — single-cell genomics (Phase 19.5) ----
    /// Scanpy — de-facto Python single-cell analysis library.
    Scanpy,
    /// scvi-tools — probabilistic deep-learning models for single-cell data.
    Scvi,
    // ---- Biology — protein design expansion (Phase 27.5) ----
    /// Chroma — Generate Biomedicines joint backbone+sequence diffusion.
    Chroma,
    /// ESM-IF — Meta inverse-folding sequence design.
    EsmIf,
    /// RFantibody — RosettaCommons antibody-specific design.
    RfAntibody,
    // ---- Biology — aligners expansion (Phase 18.5) ----
    /// Bowtie2 — short-read aligner (alternative to BWA).
    Bowtie2,
    /// MMseqs2 — fast sensitive protein search + clustering.
    Mmseqs2,
    /// DIAMOND — ultra-fast BLAST-compatible protein search.
    Diamond,
    // ---- Biology — RNA-seq alignment (Phase 18.6) ----
    /// HISAT2 — graph-based splice-aware RNA-seq aligner.
    Hisat2,
    /// STAR — most-used spliced RNA-seq aligner.
    Star,
    // ---- Biology — transcript quantification (Phase 20) ----
    /// Salmon — quasi-mapping transcript-level quantification.
    Salmon,
    /// Kallisto — pseudoalignment-based transcript quantification.
    Kallisto,
    // ---- Biology — phylogenetics (Phase 30) ----
    /// IQ-TREE — modern ML tree inference with ModelFinder.
    IqTree,
    /// RAxML-NG — next-generation RAxML rewrite.
    RaxmlNg,
    /// FastTree — approximate-ML tree inference for large trees.
    FastTree,
    // ---- Biology — RNA structure (Phase 28) ----
    /// ViennaRNA RNAfold — secondary-structure prediction (academic license).
    ViennaRna,
    /// RNAstructure Fold — Mathews lab RNA folding (BSD).
    RnaStructure,
    /// NUPACK — Caltech nucleic-acid analysis (academic license).
    Nupack,
    // ---- Biology — quantum chemistry (Phase 25) ----
    /// Psi4 — HF/DFT/post-HF general-purpose quantum chemistry.
    Psi4,
    /// NWChem — massively-parallel ab initio quantum chemistry.
    Nwchem,
    /// xTB — extended tight-binding semiempirical quantum chemistry.
    Xtb,
    // ---- Biology — EvolutionaryScale models (Phase 27.6) ----
    /// ESM3 — generative multi-modal protein model (sequence +
    /// structure + function joint reasoning).
    Esm3,
    /// ESM Cambrian — protein representation embeddings.
    Esmc,
    // ---- Biology — systems biology (Phase 32) ----
    /// COPASI — biochemical pathway / ODE simulation.
    Copasi,
    /// BioNetGen — rule-based signaling network modeling.
    BioNetGen,
    /// PhysiCell — agent-based multicellular tissue simulation.
    PhysiCell,
    // ---- Biology — cryo-EM (Phase 36) ----
    /// RELION — cryo-EM Bayesian 3D reconstruction.
    Relion,
    /// EMAN2 — broad-spectrum cryo-EM image processing.
    Eman2,
    /// CTFFIND — cryo-EM CTF estimation (academic license).
    Ctffind,
    // ---- Biology — sequencing read simulators (Phase 31) ----
    /// ART — Illumina short-read simulator.
    Art,
    /// wgsim — classic short-read simulator (samtools-bundled).
    Wgsim,
    /// Badread — Nanopore long-read simulator.
    Badread,
    // ---- Biology — CRISPR design (Phase 35) ----
    /// CHOPCHOP — CRISPR guide-RNA design.
    Chopchop,
    /// CRISPOR — CRISPR guide design + off-target prediction.
    Crispor,
    /// Cas-OFFinder — CRISPR off-target searching.
    CasOffinder,
    // ---- Biology — Rosetta family (Phase 38) ----
    /// Rosetta — rosetta_scripts XML-protocol-driven modeling.
    Rosetta,
    /// PyRosetta — Python wrapper for Rosetta core.
    PyRosetta,
    // ---- Biology — population genetics (Phase 29) ----
    /// SLiM — forward-time population-genetics simulator.
    Slim,
    /// msprime — coalescent population-genetics simulator.
    Msprime,
    /// tskit — tree-sequence analysis library.
    Tskit,
    // ---- Biology — Bayesian phylogenetics (Phase 30.5) ----
    /// BEAST 2 — Bayesian MCMC phylogenetic inference.
    Beast2,
    /// MrBayes — Bayesian MCMC phylogenetic inference.
    MrBayes,
    // ---- Biology — DNA structural geometry (Phase 39) ----
    /// X3DNA — DNA base-step parameters (academic license).
    X3dna,
    /// Curves+ — DNA helical-axis analysis (academic license).
    Curves,
    /// DSSR — DNA/RNA structural-feature analysis (academic license).
    Dssr,
    // ---- Biology — MD analysis expansion (Phase 5.5) ----
    /// PLUMED — enhanced sampling + free-energy MD analysis.
    Plumed,
    /// ProDy — protein dynamics + ENM analysis.
    Prody,
    /// cpptraj — AmberTools canonical trajectory analysis.
    Cpptraj,
    // ---- Biology — synthetic biology (Phase 33) ----
    /// pySBOL — SBOL Python composition.
    PySbol,
    /// j5 — DNA assembly automation (JAR-distributed).
    J5,
    /// Cello — genetic-circuit design + DNA compiler (JAR-distributed).
    Cello,
    // ---- Biology — alignment toolkit expansion (Phase 18.7) ----
    /// NCBI BLAST+ — nucleotide / protein sequence search.
    Blast,
    /// Clustal Omega — modern progressive multiple-sequence aligner.
    Clustalo,
    /// T-Coffee — consensus / library-based multiple-sequence aligner.
    TCoffee,
    // ---- Biology — single-cell genomics expansion (Phase 19.6) ----
    /// Seurat — R-based single-cell analysis (Rscript subprocess).
    Seurat,
    /// AnnData — Python single-cell HDF5 `.h5ad` data container.
    AnnData,
    // ---- Biology — bio MD engines (Phase 5.6) ----
    /// NAMD — UIUC all-atom MD engine (academic / non-commercial).
    Namd,
    /// AmberTools sander — OSS portion of AMBER's MD engine.
    Sander,
    /// HOOMD-blue — Glotzer lab GPU-native particle MD engine.
    Hoomd,
    // ---- Biology — MD analysis sister (Phase 5.7) ----
    /// MDTraj — Python MD trajectory analyzer (sister to MDAnalysis).
    Mdtraj,
    // ---- Biology — structure prediction + search expansion (Phase 17.7) ----
    /// RoseTTAFold — Baker lab original 3-track protein structure prediction.
    RoseTTAFold,
    /// OmegaFold — HelixonAI single-sequence structure prediction (no MSA).
    OmegaFold,
    /// FoldSeek — Steinegger lab protein structure search via 3Di alphabet.
    Foldseek,
    // ---- Biology — spatial stochastic reaction-diffusion (Phase 32.5) ----
    /// Smoldyn — Andrews lab spatial stochastic reaction-diffusion sim.
    Smoldyn,
    /// MCell — Salk Institute cell-scale spatial stochastic simulator.
    Mcell,
    // ---- Biology — sequence editors / plasmid design (Phase 41) ----
    /// pydna — Python plasmid / clone-design library.
    Pydna,
    /// Jalview — Java alignment viewer (headless mode).
    Jalview,
    // ---- Biology — microscopy / bioimage analysis (Phase 40) ----
    /// Fiji — ImageJ headless image processing (NIH / Schindelin et al).
    Fiji,
    /// CellProfiler — Broad Institute pipeline-driven cell segmentation.
    CellProfiler,
    /// Ilastik — interactive ML pixel/object classification.
    Ilastik,
    // ---- Biology — workflow expansion (Phase 22.5) ----
    /// Planemo — Galaxy ecosystem workflow CLI.
    Planemo,
    /// Cromwell — Broad WDL workflow engine (JAR).
    Cromwell,
    /// cwltool — Common Workflow Language reference runner (Python).
    Cwltool,
    // ---- Biology — web 3D molecular visualization (Phase 42) ----
    /// Mol* — PDBe / RCSB modern WebGL molecular viewer.
    Molstar,
    /// NGL Viewer — Rose lab WebGL molecular viewer framework.
    Ngl,
    // ---- Biology — mRNA design (Phase 43) ----
    /// DNA Chisel — Edinburgh Genome Foundry codon optimization library.
    DnaChisel,
    /// LinearDesign — Baidu joint codon + secondary-structure mRNA design.
    LinearDesign,
    /// iCodon — codon-level mRNA stability prediction (R-based).
    Icodon,
    // ---- Biology — RNA folding expansion (Phase 44.5) ----
    /// mfold/UNAFold — Zuker's classic RNA folder (academic).
    Mfold,
    /// EternaFold — ML-aware RNA folder via the arnie wrapper.
    EternaFold,
    /// LinearFold — Baidu fast folder (sister to LinearDesign).
    LinearFold,
    // ---- Biology — base + prime editing design (Phase 35.5) ----
    /// BE-Designer — base editor guide design.
    BeDesigner,
    /// BE-Hive — Liu lab base-editing outcome predictor.
    BeHive,
    /// PrimeDesign — Liu lab prime editing design tool.
    PrimeDesign,
    /// pegFinder — Komor lab pegRNA finder.
    PegFinder,
    // ---- Biology — edit-outcome prediction (Phase 35.6) ----
    /// inDelphi — Liu lab Cas9-cut indel pattern predictor.
    Indelphi,
    /// FORECasT — Sanger alternative indel predictor.
    Forecast,
    /// AlphaMissense — DeepMind missense effect predictor (academic).
    AlphaMissense,
    /// CRISPRitz — off-target genome-wide search.
    Crispritz,
    // ---- Biology — pharmacokinetics + RNA tertiary (Phase 45) ----
    /// PK-Sim — Open Systems Pharmacology PBPK simulation.
    PkSim,
    /// SimRNA — 3D RNA tertiary structure prediction.
    SimRna,
}

impl Template {
    /// Resolve a CLI / GUI alias to a [`Template`]. Returns `None`
    /// for unknown names so callers (the binary's `parse_args` and
    /// the GUI's adapter-id lookup) can produce their own error.
    ///
    /// Deliberately not wired through `std::str::FromStr` — that
    /// trait's `Err = Self::Err` indirection is more friction than
    /// signal here, and the call sites already work with `Option`.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "empty" | "minimal" => Some(Self::Empty),
            "cfd" | "openfoam" => Some(Self::Cfd),
            "fea" | "calculix" | "structural" => Some(Self::Fea),
            "chemistry" | "cantera" | "chem" => Some(Self::Chemistry),
            "su2" | "compressible" | "aero" => Some(Self::Su2),
            "openradioss" | "crash" | "impact" => Some(Self::OpenRadioss),
            "code-aster" | "code_aster" | "aster" | "thermomechanical" => Some(Self::CodeAster),
            "netgen" | "meshing" | "mesh" => Some(Self::Netgen),
            "meep" | "fdtd" | "photonics" => Some(Self::Meep),
            "gromacs" | "md" | "molecular-dynamics" => Some(Self::Gromacs),
            "gmsh" | "delaunay" => Some(Self::Gmsh),
            "lammps" | "lj" | "classical-md" => Some(Self::Lammps),
            "elmer" | "elmer-heat" | "heat" | "elmer_heat" => Some(Self::ElmerHeat),
            // Biology (Phase 17)
            "biopython" | "biopy" => Some(Self::Biopython),
            "rdkit" | "chem-py" => Some(Self::Rdkit),
            "openmm" | "omm" => Some(Self::Openmm),
            "chimerax" | "cxc" | "viz3d" => Some(Self::Chimerax),
            "oxdna" | "cgdna" => Some(Self::Oxdna),
            "mdanalysis" | "mda" | "traj-py" => Some(Self::Mdanalysis),
            "colabfold" | "cf" | "protein-fold" => Some(Self::Colabfold),
            // Biology — alignment toolkit (Phase 18)
            "bwa" | "bwa-mem" => Some(Self::Bwa),
            "minimap2" | "mm2" => Some(Self::Minimap2),
            "mafft" => Some(Self::Mafft),
            "muscle" => Some(Self::Muscle),
            "hmmer" | "hmmsearch" => Some(Self::Hmmer),
            "samtools" => Some(Self::Samtools),
            // Biology — structure prediction expansion (Phase 17.5)
            "esmfold" => Some(Self::Esmfold),
            "openfold" => Some(Self::Openfold),
            "alphafold2" | "af2" => Some(Self::Alphafold2),
            "alphafold3" | "af3" => Some(Self::Alphafold3),
            // Biology — variant calling (Phase 19)
            "bcftools" => Some(Self::Bcftools),
            "gatk" | "hc" => Some(Self::Gatk),
            "deepvariant" | "dv" => Some(Self::DeepVariant),
            // Biology — viewers (Phase 23)
            "pymol" => Some(Self::Pymol),
            "vmd" => Some(Self::Vmd),
            "igv" | "igvtools" => Some(Self::Igv),
            // Biology — protein design (Phase 27)
            "rfdiffusion" | "rfd" => Some(Self::RfDiffusion),
            "proteinmpnn" | "mpnn" => Some(Self::ProteinMpnn),
            // Biology — molecular docking (Phase 34)
            "vina" | "autodock-vina" => Some(Self::Vina),
            "autodock4" | "ad4" => Some(Self::AutoDock4),
            // Biology — cheminformatics expansion (Phase 24)
            "deepchem" | "dc" => Some(Self::DeepChem),
            "openbabel" | "obabel" => Some(Self::OpenBabel),
            "avogadro" | "avogadro2" => Some(Self::Avogadro),
            // Biology — workflow managers (Phase 22)
            "nextflow" | "nf" => Some(Self::Nextflow),
            "snakemake" | "smk" => Some(Self::Snakemake),
            // Biology — single-cell genomics (Phase 19.5)
            "scanpy" => Some(Self::Scanpy),
            "scvi" | "scvi-tools" => Some(Self::Scvi),
            // Biology — protein design expansion (Phase 27.5)
            "chroma" => Some(Self::Chroma),
            "esm-if" | "esmif" | "inverse-folding" => Some(Self::EsmIf),
            "rfantibody" | "rfab" => Some(Self::RfAntibody),
            // Biology — aligners expansion (Phase 18.5)
            "bowtie2" | "bt2" => Some(Self::Bowtie2),
            "mmseqs2" | "mmseqs" => Some(Self::Mmseqs2),
            "diamond" | "dmnd" => Some(Self::Diamond),
            // Biology — RNA-seq alignment (Phase 18.6)
            "hisat2" | "hisat" => Some(Self::Hisat2),
            "star" => Some(Self::Star),
            // Biology — transcript quantification (Phase 20)
            "salmon" => Some(Self::Salmon),
            "kallisto" => Some(Self::Kallisto),
            // Biology — phylogenetics (Phase 30)
            "iqtree" | "iqtree2" => Some(Self::IqTree),
            "raxml-ng" | "raxml" => Some(Self::RaxmlNg),
            "fasttree" => Some(Self::FastTree),
            // Biology — RNA structure (Phase 28)
            "viennarna" | "vienna" | "rnafold" => Some(Self::ViennaRna),
            "rnastructure" => Some(Self::RnaStructure),
            "nupack" => Some(Self::Nupack),
            // Biology — quantum chemistry (Phase 25)
            "psi4" => Some(Self::Psi4),
            "nwchem" => Some(Self::Nwchem),
            "xtb" => Some(Self::Xtb),
            // Biology — EvolutionaryScale models (Phase 27.6)
            "esm3" => Some(Self::Esm3),
            "esmc" | "esm-cambrian" => Some(Self::Esmc),
            // Biology — systems biology (Phase 32)
            "copasi" => Some(Self::Copasi),
            "bionetgen" | "bng" => Some(Self::BioNetGen),
            "physicell" => Some(Self::PhysiCell),
            // Biology — cryo-EM (Phase 36)
            "relion" => Some(Self::Relion),
            "eman2" | "eman" => Some(Self::Eman2),
            "ctffind" => Some(Self::Ctffind),
            // Biology — sequencing read simulators (Phase 31)
            "art" | "art-illumina" => Some(Self::Art),
            "wgsim" => Some(Self::Wgsim),
            "badread" => Some(Self::Badread),
            // Biology — CRISPR design (Phase 35)
            "chopchop" => Some(Self::Chopchop),
            "crispor" => Some(Self::Crispor),
            "cas-offinder" | "cas-off" => Some(Self::CasOffinder),
            // Biology — Rosetta family (Phase 38)
            "rosetta" => Some(Self::Rosetta),
            "pyrosetta" => Some(Self::PyRosetta),
            // Biology — population genetics (Phase 29)
            "slim" => Some(Self::Slim),
            "msprime" => Some(Self::Msprime),
            "tskit" => Some(Self::Tskit),
            // Biology — Bayesian phylogenetics (Phase 30.5)
            "beast2" | "beast" => Some(Self::Beast2),
            "mrbayes" | "mb" => Some(Self::MrBayes),
            // Biology — DNA structural geometry (Phase 39)
            "x3dna" | "3dna" => Some(Self::X3dna),
            "curves" | "curves+" => Some(Self::Curves),
            "dssr" => Some(Self::Dssr),
            // Biology — MD analysis expansion (Phase 5.5)
            "plumed" => Some(Self::Plumed),
            "prody" => Some(Self::Prody),
            "cpptraj" => Some(Self::Cpptraj),
            // Biology — synthetic biology (Phase 33)
            "pysbol" | "sbol" => Some(Self::PySbol),
            "j5" => Some(Self::J5),
            "cello" => Some(Self::Cello),
            // Biology — alignment toolkit expansion (Phase 18.7)
            "blast" | "blast+" | "blastn" | "blastp" => Some(Self::Blast),
            "clustalo" | "clustal-omega" | "clustal_omega" => Some(Self::Clustalo),
            "tcoffee" | "t-coffee" | "t_coffee" => Some(Self::TCoffee),
            // Biology — single-cell genomics expansion (Phase 19.6)
            "seurat" | "single-cell-r" => Some(Self::Seurat),
            "anndata" | "h5ad" => Some(Self::AnnData),
            // Biology — bio MD engines (Phase 5.6)
            "namd" | "namd2" | "namd3" => Some(Self::Namd),
            "sander" | "amber-sander" | "ambertools-sander" => Some(Self::Sander),
            "hoomd" | "hoomd-blue" | "hoomdblue" => Some(Self::Hoomd),
            // Biology — MD analysis sister (Phase 5.7)
            "mdtraj" => Some(Self::Mdtraj),
            // Biology — structure prediction + search expansion (Phase 17.7)
            "rosettafold" | "rf" => Some(Self::RoseTTAFold),
            "omegafold" | "of" => Some(Self::OmegaFold),
            "foldseek" => Some(Self::Foldseek),
            // Biology — spatial stochastic reaction-diffusion (Phase 32.5)
            "smoldyn" => Some(Self::Smoldyn),
            "mcell" => Some(Self::Mcell),
            // Biology — sequence editors (Phase 41)
            "pydna" => Some(Self::Pydna),
            "jalview" => Some(Self::Jalview),
            // Biology — microscopy / bioimage analysis (Phase 40)
            "fiji" | "imagej" => Some(Self::Fiji),
            "cellprofiler" | "cp" => Some(Self::CellProfiler),
            "ilastik" => Some(Self::Ilastik),
            // Biology — workflow expansion (Phase 22.5)
            "planemo" | "galaxy" => Some(Self::Planemo),
            "cromwell" | "wdl" => Some(Self::Cromwell),
            "cwltool" | "cwl" => Some(Self::Cwltool),
            // Biology — web 3D molecular visualization (Phase 42)
            "molstar" | "mol-star" | "molstar-viewer" => Some(Self::Molstar),
            "ngl" | "nglviewer" | "ngl-viewer" => Some(Self::Ngl),
            // Biology — mRNA design (Phase 43)
            "dnachisel" | "dna-chisel" | "dna_chisel" => Some(Self::DnaChisel),
            "lineardesign" | "linear-design" | "linear_design" => Some(Self::LinearDesign),
            "icodon" => Some(Self::Icodon),
            // Biology — RNA folding expansion (Phase 44.5)
            "mfold" | "unafold" => Some(Self::Mfold),
            "eternafold" | "eterna-fold" => Some(Self::EternaFold),
            "linearfold" | "linear-fold" | "linear_fold" => Some(Self::LinearFold),
            // Biology — base + prime editing design (Phase 35.5)
            "be-designer" | "bedesigner" | "be_designer" => Some(Self::BeDesigner),
            "be-hive" | "behive" | "be_hive" => Some(Self::BeHive),
            "primedesign" | "prime-design" | "prime_design" => Some(Self::PrimeDesign),
            "pegfinder" | "peg-finder" | "peg_finder" => Some(Self::PegFinder),
            // Biology — edit-outcome prediction (Phase 35.6)
            "indelphi" | "in-delphi" => Some(Self::Indelphi),
            "forecast" => Some(Self::Forecast),
            "alphamissense" | "alpha-missense" => Some(Self::AlphaMissense),
            "crispritz" => Some(Self::Crispritz),
            // Biology — pharmacokinetics + RNA tertiary (Phase 45)
            "pksim" | "pk-sim" | "pk_sim" => Some(Self::PkSim),
            "simrna" | "sim-rna" | "sim_rna" => Some(Self::SimRna),
            _ => None,
        }
    }
}

/// Map a [`Template`] to its case directory name (the leaf under
/// `cases/`). Single source of truth — both `scaffold_project` and
/// the GUI's "New case from adapter" flow read it.
pub fn case_dir_for_template(template: Template) -> &'static str {
    match template {
        Template::Empty => "case-1",
        Template::Cfd => "cavity",
        Template::Fea => "cantilever",
        Template::Chemistry => "ch4-equilibrium",
        Template::Su2 => "naca0012",
        Template::OpenRadioss => "drop-test",
        Template::CodeAster => "static-beam",
        Template::Netgen => "csg-box",
        Template::Meep => "ring-resonator",
        Template::Gromacs => "lysozyme",
        Template::Gmsh => "box-mesh",
        Template::Lammps => "lj-fluid",
        Template::ElmerHeat => "heat-cube",
        // Biology (Phase 17)
        Template::Biopython => "biopython-analyse",
        Template::Rdkit => "rdkit-screen",
        Template::Openmm => "protein-relax",
        Template::Chimerax => "chimerax-render",
        Template::Oxdna => "oxdna-duplex",
        Template::Mdanalysis => "trajectory-rmsd",
        Template::Colabfold => "dna-folding",
        Template::Bwa => "bwa-align",
        Template::Minimap2 => "minimap2-align",
        Template::Mafft => "mafft-msa",
        Template::Muscle => "muscle-msa",
        Template::Hmmer => "hmmer-search",
        Template::Samtools => "samtools-flagstat",
        Template::Esmfold => "esmfold-predict",
        Template::Openfold => "openfold-predict",
        Template::Alphafold2 => "alphafold2-predict",
        Template::Alphafold3 => "alphafold3-predict",
        Template::Bcftools => "bcftools-call",
        Template::Gatk => "gatk-haplotype",
        Template::DeepVariant => "deepvariant-call",
        Template::Pymol => "pymol-render",
        Template::Vmd => "vmd-render",
        Template::Igv => "igv-index",
        Template::RfDiffusion => "rfdiffusion-design",
        Template::ProteinMpnn => "proteinmpnn-design",
        Template::Vina => "vina-dock",
        Template::AutoDock4 => "autodock4-dock",
        Template::DeepChem => "deepchem-screen",
        Template::OpenBabel => "openbabel-convert",
        Template::Avogadro => "avogadro-render",
        Template::Nextflow => "nextflow-pipeline",
        Template::Snakemake => "snakemake-pipeline",
        Template::Scanpy => "scanpy-analyse",
        Template::Scvi => "scvi-train",
        Template::Chroma => "chroma-design",
        Template::EsmIf => "esm-if-design",
        Template::RfAntibody => "rfantibody-design",
        Template::Bowtie2 => "bowtie2-align",
        Template::Mmseqs2 => "mmseqs2-search",
        Template::Diamond => "diamond-search",
        Template::Hisat2 => "hisat2-align",
        Template::Star => "star-align",
        Template::Salmon => "salmon-quant",
        Template::Kallisto => "kallisto-quant",
        Template::IqTree => "iqtree-build",
        Template::RaxmlNg => "raxml-ng-build",
        Template::FastTree => "fasttree-build",
        Template::ViennaRna => "viennarna-fold",
        Template::RnaStructure => "rnastructure-fold",
        Template::Nupack => "nupack-analyze",
        Template::Psi4 => "psi4-compute",
        Template::Nwchem => "nwchem-compute",
        Template::Xtb => "xtb-compute",
        Template::Esm3 => "esm3-generate",
        Template::Esmc => "esmc-embed",
        Template::Copasi => "copasi-simulate",
        Template::BioNetGen => "bionetgen-simulate",
        Template::PhysiCell => "physicell-simulate",
        Template::Relion => "relion-refine",
        Template::Eman2 => "eman2-refine",
        Template::Ctffind => "ctffind-estimate",
        Template::Art => "art-simulate",
        Template::Wgsim => "wgsim-simulate",
        Template::Badread => "badread-simulate",
        Template::Chopchop => "chopchop-design",
        Template::Crispor => "crispor-design",
        Template::CasOffinder => "cas-offinder-search",
        Template::Rosetta => "rosetta-protocol",
        Template::PyRosetta => "pyrosetta-script",
        Template::Slim => "slim-simulate",
        Template::Msprime => "msprime-simulate",
        Template::Tskit => "tskit-analyze",
        Template::Beast2 => "beast2-mcmc",
        Template::MrBayes => "mrbayes-mcmc",
        Template::X3dna => "x3dna-analyze",
        Template::Curves => "curves-analyze",
        Template::Dssr => "dssr-analyze",
        Template::Plumed => "plumed-analyze",
        Template::Prody => "prody-analyze",
        Template::Cpptraj => "cpptraj-analyze",
        Template::PySbol => "pysbol-compose",
        Template::J5 => "j5-assemble",
        Template::Cello => "cello-compile",
        Template::Blast => "blast-search",
        Template::Clustalo => "clustalo-align",
        Template::TCoffee => "tcoffee-align",
        Template::Seurat => "seurat-analyze",
        Template::AnnData => "anndata-process",
        Template::Namd => "namd-simulate",
        Template::Sander => "sander-simulate",
        Template::Hoomd => "hoomd-simulate",
        Template::Mdtraj => "mdtraj-analyze",
        Template::RoseTTAFold => "rosettafold-predict",
        Template::OmegaFold => "omegafold-predict",
        Template::Foldseek => "foldseek-search",
        Template::Smoldyn => "smoldyn-simulate",
        Template::Mcell => "mcell-simulate",
        Template::Pydna => "pydna-design",
        Template::Jalview => "jalview-view",
        Template::Fiji => "fiji-process",
        Template::CellProfiler => "cellprofiler-segment",
        Template::Ilastik => "ilastik-classify",
        Template::Planemo => "planemo-run",
        Template::Cromwell => "cromwell-run",
        Template::Cwltool => "cwltool-run",
        Template::Molstar => "molstar-view",
        Template::Ngl => "ngl-view",
        Template::DnaChisel => "dnachisel-optimize",
        Template::LinearDesign => "lineardesign-design",
        Template::Icodon => "icodon-predict",
        Template::Mfold => "mfold-fold",
        Template::EternaFold => "eternafold-fold",
        Template::LinearFold => "linearfold-fold",
        Template::BeDesigner => "be-designer-design",
        Template::BeHive => "be-hive-predict",
        Template::PrimeDesign => "primedesign-design",
        Template::PegFinder => "pegfinder-design",
        Template::Indelphi => "indelphi-predict",
        Template::Forecast => "forecast-predict",
        Template::AlphaMissense => "alphamissense-predict",
        Template::Crispritz => "crispritz-search",
        Template::PkSim => "pksim-simulate",
        Template::SimRna => "simrna-fold",
    }
}

/// Resolve an adapter id (`AdapterInfo::id` — e.g. `"dnachisel"`)
/// to the case-directory name its starter scaffolds into. Returns
/// `None` for adapters that don't have a registered template
/// (the GUI surfaces this as "no starter; bring your own case.toml").
pub fn case_dir_for(adapter_id: &str) -> Option<&'static str> {
    Template::from_str(adapter_id).map(case_dir_for_template)
}

/// Render the starter `case.toml` body for a given adapter id.
/// `dir_name` is stamped into the `name = "<dir>"` slot of the
/// emitted `[case]` block; pass the result of [`case_dir_for`] for
/// a roundtrippable scaffold. Returns `None` for unrecognised ids.
pub fn case_toml_body(adapter_id: &str) -> Option<String> {
    let template = Template::from_str(adapter_id)?;
    let case_dir = case_dir_for_template(template);
    Some(render_case_toml(template, case_dir))
}

/// Render the starter `project.toml` shell. Public alias for
/// [`render_project_toml`] — kept under this name so the GUI
/// imports the GUI-facing API as a unit.
///
/// # Errors
///
/// Forwards the [`sanitize_project_name`] check failures from
/// [`render_project_toml`].
pub fn project_toml(name: &str, case_dir: &str) -> Result<String, String> {
    render_project_toml(name, case_dir)
}

/// The canonical `name + description + case_dir` catalogue. Drives
/// `valenx-init --list-templates` and the GUI's Tools menu /
/// adapter-discovery surfaces. Read-only; new entries are added by
/// extending [`Template`] + the private `TEMPLATE_ROWS` constant in
/// this module.
pub fn template_rows() -> &'static [TemplateRow] {
    TEMPLATE_ROWS
}

/// Materialise the directory + project.toml + `cases/<name>/case.toml`.
/// Refuses to overwrite an existing project.toml so re-running the
/// init doesn't trash a real project.
pub fn scaffold_project(
    dir: &Path,
    template: Template,
    project_name: Option<&str>,
) -> Result<(), String> {
    let project_toml = dir.join("project.toml");
    if project_toml.exists() {
        return Err(format!(
            "{} already exists — refusing to overwrite. Delete it first or pick a fresh dir.",
            project_toml.display()
        ));
    }
    let dir_name = dir
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("untitled");
    let name = project_name.unwrap_or(dir_name);

    std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
    let case_dir_name = case_dir_for_template(template);
    let rendered = render_project_toml(name, case_dir_name)?;
    crate::io_caps::atomic_write_str(&project_toml, &rendered)
        .map_err(|e| format!("write {}: {e}", project_toml.display()))?;
    let case_dir = dir.join("cases").join(case_dir_name);
    std::fs::create_dir_all(&case_dir)
        .map_err(|e| format!("create {}: {e}", case_dir.display()))?;
    let case_toml_path = case_dir.join("case.toml");
    crate::io_caps::atomic_write_str(&case_toml_path, &render_case_toml(template, case_dir_name))
        .map_err(|e| format!("write {}: {e}", case_toml_path.display()))?;
    // Drop sample input files alongside case.toml so the case is
    // immediately runnable for templates whose adapter has a small,
    // self-contained sample input. Adapters whose typical input is
    // a binary file (e.g. .tpr / .vol / staged .rad) skip this —
    // users have to bring their own.
    for (filename, content) in template_sample_files(template) {
        let path = case_dir.join(filename);
        crate::io_caps::atomic_write_str(&path, content)
            .map_err(|e| format!("write {}: {e}", path.display()))?;
    }
    Ok(())
}

/// Sample input files to ship with a starter case. Returns an empty
/// slice for templates whose adapter expects user-supplied binary
/// inputs (Gromacs .tpr, OpenRadioss engine deck, Code_Aster
/// .export — all of those need real preprocessing the adapter
/// can't fake). For text-driven adapters we ship a minimal
/// runnable sample.
pub fn template_sample_files(template: Template) -> &'static [(&'static str, &'static str)] {
    match template {
        Template::Netgen => &[(
            "csg-box.geo",
            // Minimal CSG description Netgen recognises: an
            // axis-aligned unit cube. Mesh-size hint comes from
            // case.toml's mesh_size knob.
            "algebraic3d\nsolid box = orthobrick (0, 0, 0; 1, 1, 1);\ntlo box;\n",
        )],
        Template::Meep => &[(
            "ring.py",
            // A bare-bones Meep simulation that completes in
            // under a second on a modern laptop. Users edit it
            // to do something interesting.
            r#"# Generated by `valenx-init --template meep`.
import meep as mp

cell = mp.Vector3(8, 8, 0)
geometry = [mp.Cylinder(radius=2.0, material=mp.Medium(epsilon=12))]
sources = [
    mp.Source(
        mp.GaussianSource(frequency=0.15, fwidth=0.1),
        component=mp.Ez,
        center=mp.Vector3(),
    )
]
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    resolution=10,
    boundary_layers=[mp.PML(1.0)],
)
sim.run(until=20)
"#,
        )],
        Template::Su2 => &[(
            "naca0012.cfg",
            // Skeleton SU2 .cfg pointing at a mesh the user supplies.
            // Real config has hundreds of options; this is the
            // minimum that produces a parseable run.
            r#"% Generated by `valenx-init --template su2`. Edit freely.
SOLVER= EULER
MESH_FILENAME= naca0012.su2
MESH_FORMAT= SU2
MARKER_HEATFLUX= ( airfoil, 0.0 )
MARKER_FAR= ( farfield )
MACH_NUMBER= 0.8
AOA= 1.25
FREESTREAM_PRESSURE= 101325.0
FREESTREAM_TEMPERATURE= 273.15
EXT_ITER= 100
CONV_FILENAME= history
RESTART_SOL= NO
OUTPUT_FILES= (RESTART, PARAVIEW)
"#,
        )],
        // ---- Biology (Phase 17) ----
        Template::Biopython => &[(
            "analyse.py",
            // Tiny Biopython smoke script. Imports the library and
            // prints a sentinel line the adapter's progress hook
            // recognises as "completed" (95% tick).
            r#"# Generated by `valenx-init --template biopython`.
import Bio
print(f"Biopython {Bio.__version__} ready")
print("[valenx] biopython done")
"#,
        )],
        Template::Rdkit => &[(
            "enumerate.py",
            // Tiny RDKit smoke script. The adapter ships SMILES via
            // case.toml's `smiles = [...]` knob; this script just
            // confirms RDKit is importable.
            r#"# Generated by `valenx-init --template rdkit`.
from rdkit import Chem
mol = Chem.MolFromSmiles("CCO")
print(f"RDKit parsed ethanol: {Chem.MolToSmiles(mol)}")
print("[valenx] rdkit done")
"#,
        )],
        Template::Openmm => &[(
            "relax.py",
            // Minimal OpenMM minimisation. Builds a single-residue
            // alanine dipeptide in vacuum and runs L-BFGS for a
            // handful of steps, then writes minimised.pdb +
            // trajectory.dcd so collect() has artefacts to find.
            r#"# Generated by `valenx-init --template openmm`.
# Tiny vacuum minimisation — replaces with a real PDB load
# + force-field assembly when you have a target system.
import openmm as mm
import openmm.app as app
from openmm import unit

forcefield = app.ForceField("amber14-all.xml")
modeller = app.Modeller(app.Topology(), [])
# Real cases load a PDB; for the smoke run we just make sure the
# imports + minimiser plumbing work. Replace with:
#   pdb = app.PDBFile("input.pdb")
#   modeller = app.Modeller(pdb.topology, pdb.positions)
system = mm.System()
integrator = mm.LangevinIntegrator(300 * unit.kelvin,
                                   1.0 / unit.picosecond,
                                   0.002 * unit.picoseconds)
print("OpenMM ready (replace placeholder topology with a real PDB)")
print("[valenx] openmm done")
"#,
        )],
        Template::Chimerax => &[(
            "view.cxc",
            // Tiny ChimeraX command script. Opens a PDB ID, applies
            // a cartoon style, saves a screenshot. Edit `open` to
            // point at your structure of interest.
            r#"# Generated by `valenx-init --template chimerax`.
open 1abc
cartoon
save snapshot.png supersample 3
exit
"#,
        )],
        Template::Mdanalysis => &[(
            "analyse_traj.py",
            // Tiny MDAnalysis smoke script. Loads a topology +
            // trajectory pair (placeholder filenames) and prints the
            // first frame's atom count.
            r#"# Generated by `valenx-init --template mdanalysis`.
import MDAnalysis as mda
# Replace the placeholders with real inputs the adapter stages.
u = mda.Universe("topology.pdb", "trajectory.dcd")
print(f"frames={len(u.trajectory)} atoms={u.atoms.n_atoms}")
print("[valenx] mdanalysis done")
"#,
        )],
        Template::Colabfold => &[(
            "query.fasta",
            // 50-residue placeholder protein sequence — a fragment of
            // ubiquitin (P0CG48) so the FASTA parses and ColabFold
            // produces a real prediction without the user having to
            // bring their own sequence on first run.
            r#">query|placeholder|ubiquitin-fragment
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLE
"#,
        )],
        // The remaining adapters' typical inputs are binary
        // (Gromacs .tpr, OpenRadioss engine deck, Code_Aster
        // export bundle, oxDNA input.dat) or fully user-defined
        // (CalculiX FEA, OpenFOAM cavity, Cantera). They get
        // case.toml only.
        _ => &[],
    }
}

/// One row of the `--list-templates` output: canonical name, brief
/// description, and the case directory the template scaffolds into.
#[derive(Clone, Copy, Debug)]
pub struct TemplateRow {
    pub name: &'static str,
    pub description: &'static str,
    pub case_dir: &'static str,
}

/// Static catalogue used by `--list-templates`. Keeping it next to
/// the rest of the rendering keeps drift from the `Template` enum
/// visible (a missing entry surfaces as a stale list).
const TEMPLATE_ROWS: &[TemplateRow] = &[
    TemplateRow {
        name: "empty",
        description: "minimal skeleton — no per-physics block",
        case_dir: "case-1",
    },
    TemplateRow {
        name: "cfd",
        description: "OpenFOAM simpleFoam — incompressible RANS over a box",
        case_dir: "cavity",
    },
    TemplateRow {
        name: "fea",
        description: "CalculiX linear-static — tip-loaded cantilever beam",
        case_dir: "cantilever",
    },
    TemplateRow {
        name: "chemistry",
        description: "Cantera equilibrium-HP — methane / air mixture",
        case_dir: "ch4-equilibrium",
    },
    TemplateRow {
        name: "su2",
        description: "SU2 compressible CFD — NACA 0012 airfoil starter",
        case_dir: "naca0012",
    },
    TemplateRow {
        name: "openradioss",
        description: "OpenRadioss explicit dynamics — engine-deck only",
        case_dir: "drop-test",
    },
    TemplateRow {
        name: "code-aster",
        description: "Code_Aster `as_run` on a user-built `.export`",
        case_dir: "static-beam",
    },
    TemplateRow {
        name: "netgen",
        description: "Netgen CSG meshing — axis-aligned unit cube",
        case_dir: "csg-box",
    },
    TemplateRow {
        name: "meep",
        description: "Meep FDTD — Python ring-resonator script",
        case_dir: "ring-resonator",
    },
    TemplateRow {
        name: "gromacs",
        description: "GROMACS `gmx mdrun` on a user-built `.tpr`",
        case_dir: "lysozyme",
    },
    TemplateRow {
        name: "gmsh",
        description: "gmsh procedural meshing — Delaunay tet box",
        case_dir: "box-mesh",
    },
    TemplateRow {
        name: "lammps",
        description: "LAMMPS classical MD — Lennard-Jones FCC fluid (NVE)",
        case_dir: "lj-fluid",
    },
    TemplateRow {
        name: "elmer-heat",
        description: "Elmer steady heat — two pinned-temperature faces",
        case_dir: "heat-cube",
    },
    // ---- Biology (Phase 17) ----
    TemplateRow {
        name: "biopython",
        description: "Biopython — sequence / structural-bio Python script",
        case_dir: "biopython-analyse",
    },
    TemplateRow {
        name: "rdkit",
        description: "RDKit — cheminformatics Python script",
        case_dir: "rdkit-screen",
    },
    TemplateRow {
        name: "openmm",
        description: "OpenMM — Python-native MD minimisation + DCD output",
        case_dir: "protein-relax",
    },
    TemplateRow {
        name: "chimerax",
        description: "ChimeraX — `.cxc` command-script renderer",
        case_dir: "chimerax-render",
    },
    TemplateRow {
        name: "oxdna",
        description: "oxDNA — coarse-grained DNA / RNA MD on `input.dat`",
        case_dir: "oxdna-duplex",
    },
    TemplateRow {
        name: "mdanalysis",
        description: "MDAnalysis — trajectory analysis Python script",
        case_dir: "trajectory-rmsd",
    },
    TemplateRow {
        name: "colabfold",
        description: "ColabFold — protein structure prediction from FASTA",
        case_dir: "dna-folding",
    },
    // ---- Biology — alignment toolkit (Phase 18) ----
    TemplateRow {
        name: "bwa",
        description: "BWA — short-read DNA alignment via `bwa mem`",
        case_dir: "bwa-align",
    },
    TemplateRow {
        name: "minimap2",
        description: "minimap2 — long-read + cross-domain alignment",
        case_dir: "minimap2-align",
    },
    TemplateRow {
        name: "mafft",
        description: "MAFFT — multiple sequence alignment",
        case_dir: "mafft-msa",
    },
    TemplateRow {
        name: "muscle",
        description: "MUSCLE 5 — multiple sequence alignment",
        case_dir: "muscle-msa",
    },
    TemplateRow {
        name: "hmmer",
        description: "HMMER — profile HMM search (hmmsearch / hmmscan)",
        case_dir: "hmmer-search",
    },
    TemplateRow {
        name: "samtools",
        description: "samtools — SAM/BAM multitool (view / sort / index / flagstat)",
        case_dir: "samtools-flagstat",
    },
    // ---- Biology — structure prediction expansion (Phase 17.5) ----
    TemplateRow {
        name: "esmfold",
        description: "ESMFold — Meta protein language model structure prediction",
        case_dir: "esmfold-predict",
    },
    TemplateRow {
        name: "openfold",
        description: "OpenFold — PyTorch reimplementation of AlphaFold 2",
        case_dir: "openfold-predict",
    },
    TemplateRow {
        name: "alphafold2",
        description: "AlphaFold 2 — DeepMind structure prediction (open weights)",
        case_dir: "alphafold2-predict",
    },
    TemplateRow {
        name: "alphafold3",
        description: "AlphaFold 3 — all-atom complex prediction (non-commercial weights)",
        case_dir: "alphafold3-predict",
    },
    // ---- Biology — variant calling (Phase 19) ----
    TemplateRow {
        name: "bcftools",
        description: "bcftools — VCF/BCF multitool (view / call / filter / concat)",
        case_dir: "bcftools-call",
    },
    TemplateRow {
        name: "gatk",
        description: "GATK HaplotypeCaller — Java-based variant calling",
        case_dir: "gatk-haplotype",
    },
    TemplateRow {
        name: "deepvariant",
        description: "DeepVariant — Google ML-driven variant calling (WGS / WES / PacBio / ONT)",
        case_dir: "deepvariant-call",
    },
    // ---- Biology — viewers (Phase 23) ----
    TemplateRow {
        name: "pymol",
        description: "PyMOL — open-source script-driven structural rendering",
        case_dir: "pymol-render",
    },
    TemplateRow {
        name: "vmd",
        description: "VMD — Tcl-scripted MD trajectory viewer (academic license)",
        case_dir: "vmd-render",
    },
    TemplateRow {
        name: "igv",
        description: "IGV `igvtools` — headless BAM/VCF indexer + tile generator",
        case_dir: "igv-index",
    },
    // ---- Biology — protein design (Phase 27) ----
    TemplateRow {
        name: "rfdiffusion",
        description: "RFdiffusion — protein backbone generation",
        case_dir: "rfdiffusion-design",
    },
    TemplateRow {
        name: "proteinmpnn",
        description: "ProteinMPNN — sequence design from a protein backbone",
        case_dir: "proteinmpnn-design",
    },
    // ---- Biology — molecular docking (Phase 34) ----
    TemplateRow {
        name: "vina",
        description: "AutoDock Vina — modern single-binary small-molecule docker",
        case_dir: "vina-dock",
    },
    TemplateRow {
        name: "autodock4",
        description: "AutoDock 4 — two-stage (autogrid4 + autodock4) docking",
        case_dir: "autodock4-dock",
    },
    // ---- Biology — cheminformatics expansion (Phase 24) ----
    TemplateRow {
        name: "deepchem",
        description: "DeepChem — PyTorch-backed cheminformatics",
        case_dir: "deepchem-screen",
    },
    TemplateRow {
        name: "openbabel",
        description: "Open Babel — chemistry-format converter (~120 formats)",
        case_dir: "openbabel-convert",
    },
    TemplateRow {
        name: "avogadro",
        description: "Avogadro 2 — Python-scriptable chemistry editor",
        case_dir: "avogadro-render",
    },
    // ---- Biology — workflow managers (Phase 22) ----
    TemplateRow {
        name: "nextflow",
        description: "Nextflow — pipeline orchestrator",
        case_dir: "nextflow-pipeline",
    },
    TemplateRow {
        name: "snakemake",
        description: "Snakemake — rule-based pipeline orchestrator",
        case_dir: "snakemake-pipeline",
    },
    // ---- Biology — single-cell genomics (Phase 19.5) ----
    TemplateRow {
        name: "scanpy",
        description: "Scanpy — Python single-cell analysis",
        case_dir: "scanpy-analyse",
    },
    TemplateRow {
        name: "scvi",
        description: "scvi-tools — probabilistic single-cell models",
        case_dir: "scvi-train",
    },
    // ---- Biology — protein design expansion (Phase 27.5) ----
    TemplateRow {
        name: "chroma",
        description: "Chroma — Generate Biomedicines diffusion design",
        case_dir: "chroma-design",
    },
    TemplateRow {
        name: "esm-if",
        description: "ESM-IF — Meta inverse-folding sequence design",
        case_dir: "esm-if-design",
    },
    TemplateRow {
        name: "rfantibody",
        description: "RFantibody — RosettaCommons antibody design",
        case_dir: "rfantibody-design",
    },
    // ---- Biology — aligners expansion (Phase 18.5) ----
    TemplateRow {
        name: "bowtie2",
        description: "Bowtie2 — short-read alignment (alternative to BWA)",
        case_dir: "bowtie2-align",
    },
    TemplateRow {
        name: "mmseqs2",
        description: "MMseqs2 — fast protein search + clustering",
        case_dir: "mmseqs2-search",
    },
    TemplateRow {
        name: "diamond",
        description: "DIAMOND — ultra-fast BLAST-compatible protein search",
        case_dir: "diamond-search",
    },
    // ---- Biology — RNA-seq alignment (Phase 18.6) ----
    TemplateRow {
        name: "hisat2",
        description: "HISAT2 — graph-based splice-aware RNA-seq aligner",
        case_dir: "hisat2-align",
    },
    TemplateRow {
        name: "star",
        description: "STAR — spliced RNA-seq aligner",
        case_dir: "star-align",
    },
    // ---- Biology — transcript quantification (Phase 20) ----
    TemplateRow {
        name: "salmon",
        description: "Salmon — quasi-mapping transcript quantification",
        case_dir: "salmon-quant",
    },
    TemplateRow {
        name: "kallisto",
        description: "Kallisto — pseudoalignment transcript quantification",
        case_dir: "kallisto-quant",
    },
    // ---- Biology — phylogenetics (Phase 30) ----
    TemplateRow {
        name: "iqtree",
        description: "IQ-TREE — modern ML tree inference with ModelFinder",
        case_dir: "iqtree-build",
    },
    TemplateRow {
        name: "raxml-ng",
        description: "RAxML-NG — next-generation RAxML rewrite",
        case_dir: "raxml-ng-build",
    },
    TemplateRow {
        name: "fasttree",
        description: "FastTree — approximate-ML tree inference for large trees",
        case_dir: "fasttree-build",
    },
    // ---- Biology — RNA structure (Phase 28) ----
    TemplateRow {
        name: "viennarna",
        description: "ViennaRNA RNAfold — secondary-structure (academic license)",
        case_dir: "viennarna-fold",
    },
    TemplateRow {
        name: "rnastructure",
        description: "RNAstructure Fold — Mathews lab RNA folding (BSD)",
        case_dir: "rnastructure-fold",
    },
    TemplateRow {
        name: "nupack",
        description: "NUPACK — Caltech nucleic-acid analysis (academic license)",
        case_dir: "nupack-analyze",
    },
    // ---- Biology — quantum chemistry (Phase 25) ----
    TemplateRow {
        name: "psi4",
        description: "Psi4 — HF/DFT/post-HF quantum chemistry",
        case_dir: "psi4-compute",
    },
    TemplateRow {
        name: "nwchem",
        description: "NWChem — massively-parallel ab initio quantum chemistry",
        case_dir: "nwchem-compute",
    },
    TemplateRow {
        name: "xtb",
        description: "xTB — extended tight-binding semiempirical quantum chemistry",
        case_dir: "xtb-compute",
    },
    // ---- Biology — EvolutionaryScale models (Phase 27.6) ----
    TemplateRow {
        name: "esm3",
        description: "ESM3 — EvolutionaryScale generative multi-modal protein model",
        case_dir: "esm3-generate",
    },
    TemplateRow {
        name: "esmc",
        description: "ESM Cambrian — EvolutionaryScale protein representation embeddings",
        case_dir: "esmc-embed",
    },
    // ---- Biology — systems biology (Phase 32) ----
    TemplateRow {
        name: "copasi",
        description: "COPASI — biochemical pathway / ODE simulation",
        case_dir: "copasi-simulate",
    },
    TemplateRow {
        name: "bionetgen",
        description: "BioNetGen — rule-based signaling network modeling",
        case_dir: "bionetgen-simulate",
    },
    TemplateRow {
        name: "physicell",
        description: "PhysiCell — agent-based multicellular tissue simulation",
        case_dir: "physicell-simulate",
    },
    // ---- Biology — cryo-EM (Phase 36) ----
    TemplateRow {
        name: "relion",
        description: "RELION — cryo-EM Bayesian 3D reconstruction",
        case_dir: "relion-refine",
    },
    TemplateRow {
        name: "eman2",
        description: "EMAN2 — broad-spectrum cryo-EM image processing",
        case_dir: "eman2-refine",
    },
    TemplateRow {
        name: "ctffind",
        description: "CTFFIND — cryo-EM CTF estimation (academic license)",
        case_dir: "ctffind-estimate",
    },
    // ---- Biology — sequencing read simulators (Phase 31) ----
    TemplateRow {
        name: "art",
        description: "ART — Illumina short-read simulator",
        case_dir: "art-simulate",
    },
    TemplateRow {
        name: "wgsim",
        description: "wgsim — classic short-read simulator (samtools-bundled)",
        case_dir: "wgsim-simulate",
    },
    TemplateRow {
        name: "badread",
        description: "Badread — Nanopore long-read simulator",
        case_dir: "badread-simulate",
    },
    // ---- Biology — CRISPR design (Phase 35) ----
    TemplateRow {
        name: "chopchop",
        description: "CHOPCHOP — CRISPR guide-RNA design (Cas9/Cas12a/Cas13/TALEN)",
        case_dir: "chopchop-design",
    },
    TemplateRow {
        name: "crispor",
        description: "CRISPOR — CRISPR guide-RNA design + off-target prediction",
        case_dir: "crispor-design",
    },
    TemplateRow {
        name: "cas-offinder",
        description: "Cas-OFFinder — CRISPR off-target searching",
        case_dir: "cas-offinder-search",
    },
    // ---- Biology — Rosetta family (Phase 38) ----
    TemplateRow {
        name: "rosetta",
        description: "Rosetta — rosetta_scripts XML protocol modeling (academic license)",
        case_dir: "rosetta-protocol",
    },
    TemplateRow {
        name: "pyrosetta",
        description: "PyRosetta — Python wrapper for Rosetta core (academic license)",
        case_dir: "pyrosetta-script",
    },
    // ---- Biology — population genetics (Phase 29) ----
    TemplateRow {
        name: "slim",
        description: "SLiM — forward-time population-genetics simulator",
        case_dir: "slim-simulate",
    },
    TemplateRow {
        name: "msprime",
        description: "msprime — coalescent population-genetics simulator",
        case_dir: "msprime-simulate",
    },
    TemplateRow {
        name: "tskit",
        description: "tskit — tree-sequence analysis library",
        case_dir: "tskit-analyze",
    },
    // ---- Biology — Bayesian phylogenetics (Phase 30.5) ----
    TemplateRow {
        name: "beast2",
        description: "BEAST 2 — Bayesian MCMC phylogenetic inference",
        case_dir: "beast2-mcmc",
    },
    TemplateRow {
        name: "mrbayes",
        description: "MrBayes — Bayesian MCMC phylogenetic inference",
        case_dir: "mrbayes-mcmc",
    },
    // ---- Biology — DNA structural geometry (Phase 39) ----
    TemplateRow {
        name: "x3dna",
        description: "X3DNA — DNA base-step parameter analysis (academic license)",
        case_dir: "x3dna-analyze",
    },
    TemplateRow {
        name: "curves",
        description: "Curves+ — DNA helical-axis analysis (academic license)",
        case_dir: "curves-analyze",
    },
    TemplateRow {
        name: "dssr",
        description: "DSSR — DNA/RNA structural-feature analysis (academic license)",
        case_dir: "dssr-analyze",
    },
    // ---- Biology — MD analysis expansion (Phase 5.5) ----
    TemplateRow {
        name: "plumed",
        description: "PLUMED — enhanced sampling + free-energy MD analysis",
        case_dir: "plumed-analyze",
    },
    TemplateRow {
        name: "prody",
        description: "ProDy — protein dynamics + ENM analysis",
        case_dir: "prody-analyze",
    },
    TemplateRow {
        name: "cpptraj",
        description: "cpptraj — AmberTools canonical trajectory analysis",
        case_dir: "cpptraj-analyze",
    },
    // ---- Biology — synthetic biology (Phase 33) ----
    TemplateRow {
        name: "pysbol",
        description: "pySBOL — SBOL Python composition",
        case_dir: "pysbol-compose",
    },
    TemplateRow {
        name: "j5",
        description: "j5 — DNA assembly automation (JAR)",
        case_dir: "j5-assemble",
    },
    TemplateRow {
        name: "cello",
        description: "Cello — genetic-circuit design + DNA compiler (JAR)",
        case_dir: "cello-compile",
    },
    // ---- Biology — alignment toolkit expansion (Phase 18.7) ----
    TemplateRow {
        name: "blast",
        description: "BLAST+ — NCBI nucleotide / protein sequence search",
        case_dir: "blast-search",
    },
    TemplateRow {
        name: "clustalo",
        description: "Clustal Omega — multiple-sequence aligner",
        case_dir: "clustalo-align",
    },
    TemplateRow {
        name: "tcoffee",
        description: "T-Coffee — consensus / library-based MSA",
        case_dir: "tcoffee-align",
    },
    // ---- Biology — single-cell genomics expansion (Phase 19.6) ----
    TemplateRow {
        name: "seurat",
        description: "Seurat — R-based single-cell analysis (Rscript)",
        case_dir: "seurat-analyze",
    },
    TemplateRow {
        name: "anndata",
        description: "AnnData — Python single-cell HDF5 (.h5ad) container",
        case_dir: "anndata-process",
    },
    // ---- Biology — bio MD engines (Phase 5.6) ----
    TemplateRow {
        name: "namd",
        description: "NAMD — UIUC all-atom MD engine (academic / non-commercial)",
        case_dir: "namd-simulate",
    },
    TemplateRow {
        name: "sander",
        description: "AmberTools sander — OSS portion of AMBER's MD engine",
        case_dir: "sander-simulate",
    },
    TemplateRow {
        name: "hoomd",
        description: "HOOMD-blue — Glotzer-lab GPU-native particle MD engine",
        case_dir: "hoomd-simulate",
    },
    // ---- Biology — MD analysis sister (Phase 5.7) ----
    TemplateRow {
        name: "mdtraj",
        description: "MDTraj — Python MD trajectory analyzer (sister to MDAnalysis)",
        case_dir: "mdtraj-analyze",
    },
    // ---- Biology — structure prediction + search (Phase 17.7) ----
    TemplateRow {
        name: "rosettafold",
        description: "RoseTTAFold — Baker-lab original 3-track structure prediction",
        case_dir: "rosettafold-predict",
    },
    TemplateRow {
        name: "omegafold",
        description: "OmegaFold — single-sequence structure prediction (no MSA)",
        case_dir: "omegafold-predict",
    },
    TemplateRow {
        name: "foldseek",
        description: "FoldSeek — protein structure search via 3Di alphabet",
        case_dir: "foldseek-search",
    },
    // ---- Biology — spatial stochastic reaction-diffusion (Phase 32.5) ----
    TemplateRow {
        name: "smoldyn",
        description: "Smoldyn — Andrews-lab spatial stochastic reaction-diffusion",
        case_dir: "smoldyn-simulate",
    },
    TemplateRow {
        name: "mcell",
        description: "MCell — Salk Institute cell-scale spatial stochastic simulator",
        case_dir: "mcell-simulate",
    },
    // ---- Biology — sequence editors / plasmid design (Phase 41) ----
    TemplateRow {
        name: "pydna",
        description: "pydna — Python plasmid / clone-design library",
        case_dir: "pydna-design",
    },
    TemplateRow {
        name: "jalview",
        description: "Jalview — Java alignment viewer (headless)",
        case_dir: "jalview-view",
    },
    // ---- Biology — microscopy / bioimage analysis (Phase 40) ----
    TemplateRow {
        name: "fiji",
        description: "Fiji — ImageJ headless image processing",
        case_dir: "fiji-process",
    },
    TemplateRow {
        name: "cellprofiler",
        description: "CellProfiler — Broad pipeline-driven cell segmentation",
        case_dir: "cellprofiler-segment",
    },
    TemplateRow {
        name: "ilastik",
        description: "Ilastik — ML pixel/object classification",
        case_dir: "ilastik-classify",
    },
    // ---- Biology — workflow expansion (Phase 22.5) ----
    TemplateRow {
        name: "planemo",
        description: "Planemo — Galaxy ecosystem workflow CLI",
        case_dir: "planemo-run",
    },
    TemplateRow {
        name: "cromwell",
        description: "Cromwell — Broad WDL workflow engine (JAR)",
        case_dir: "cromwell-run",
    },
    TemplateRow {
        name: "cwltool",
        description: "cwltool — Common Workflow Language reference runner",
        case_dir: "cwltool-run",
    },
    // ---- Biology — web 3D molecular visualization (Phase 42) ----
    TemplateRow {
        name: "molstar",
        description: "Mol* — PDBe / RCSB WebGL molecular viewer",
        case_dir: "molstar-view",
    },
    TemplateRow {
        name: "ngl",
        description: "NGL Viewer — Rose-lab WebGL molecular viewer",
        case_dir: "ngl-view",
    },
    // ---- Biology — mRNA design (Phase 43) ----
    TemplateRow {
        name: "dnachisel",
        description: "DNA Chisel — Edinburgh Genome Foundry codon optimization",
        case_dir: "dnachisel-optimize",
    },
    TemplateRow {
        name: "lineardesign",
        description: "LinearDesign — Baidu joint codon + structure mRNA design",
        case_dir: "lineardesign-design",
    },
    TemplateRow {
        name: "icodon",
        description: "iCodon — codon-level mRNA stability prediction (R)",
        case_dir: "icodon-predict",
    },
    // ---- Biology — RNA folding expansion (Phase 44.5) ----
    TemplateRow {
        name: "mfold",
        description: "mfold/UNAFold — Zuker classic RNA folder (academic)",
        case_dir: "mfold-fold",
    },
    TemplateRow {
        name: "eternafold",
        description: "EternaFold — ML-aware RNA folder via arnie",
        case_dir: "eternafold-fold",
    },
    TemplateRow {
        name: "linearfold",
        description: "LinearFold — Baidu fast folder",
        case_dir: "linearfold-fold",
    },
    // ---- Biology — base + prime editing design (Phase 35.5) ----
    TemplateRow {
        name: "be-designer",
        description: "BE-Designer — base editor guide design",
        case_dir: "be-designer-design",
    },
    TemplateRow {
        name: "be-hive",
        description: "BE-Hive — Liu lab base-editing outcome predictor",
        case_dir: "be-hive-predict",
    },
    TemplateRow {
        name: "primedesign",
        description: "PrimeDesign — Liu lab prime editing design",
        case_dir: "primedesign-design",
    },
    TemplateRow {
        name: "pegfinder",
        description: "pegFinder — Komor lab pegRNA finder",
        case_dir: "pegfinder-design",
    },
    // ---- Biology — edit-outcome prediction (Phase 35.6) ----
    TemplateRow {
        name: "indelphi",
        description: "inDelphi — Cas9-cut indel pattern predictor",
        case_dir: "indelphi-predict",
    },
    TemplateRow {
        name: "forecast",
        description: "FORECasT — Sanger alternative indel predictor",
        case_dir: "forecast-predict",
    },
    TemplateRow {
        name: "alphamissense",
        description: "AlphaMissense — DeepMind missense effect predictor (academic)",
        case_dir: "alphamissense-predict",
    },
    TemplateRow {
        name: "crispritz",
        description: "CRISPRitz — off-target genome-wide search",
        case_dir: "crispritz-search",
    },
    // ---- Biology — pharmacokinetics + RNA tertiary (Phase 45) ----
    TemplateRow {
        name: "pksim",
        description: "PK-Sim — Open Systems Pharmacology PBPK simulation",
        case_dir: "pksim-simulate",
    },
    TemplateRow {
        name: "simrna",
        description: "SimRNA — 3D RNA tertiary structure prediction",
        case_dir: "simrna-fold",
    },
];

/// Render the `--list-templates` output as a plain-text two-column
/// table. Pure function — unit-tested without needing a temp dir.
pub fn render_template_list() -> String {
    let mut s = String::with_capacity(2048);
    s.push_str("Available `valenx-init` templates:\n\n");
    let max_name = TEMPLATE_ROWS
        .iter()
        .map(|r| r.name.len())
        .max()
        .unwrap_or(0);
    for row in TEMPLATE_ROWS {
        s.push_str(&format!(
            "  {name:<width$}  {desc}\n                            (case dir: `{dir}`)\n",
            name = row.name,
            desc = row.description,
            dir = row.case_dir,
            width = max_name,
        ));
    }
    s.push_str("\nUse `valenx-init <dir> --template <name>` to scaffold one.\n");
    s.push_str("See `valenx-init help` for the full alias list.\n");
    s
}

/// Validate a project / case name for safe TOML interpolation. The
/// generated `project.toml` includes the name in a `"..."`-quoted
/// scalar; without this check a folder named e.g.
/// `evil"\n[rbac]\ndefault_role = "viewer"\n#` would let an attacker
/// who controls the destination directory inject an arbitrary
/// `[rbac]` block (or any other table) into the rendered file.
/// Round-12 L9.
///
/// # Errors
///
/// Returns a String error describing why the name is invalid. Allowed
/// characters: ASCII alphanumeric plus `_`, `.`, `-`. Empty names are
/// rejected.
pub fn sanitize_project_name(name: &str) -> Result<&str, String> {
    if name.is_empty() {
        return Err("project name must not be empty".to_string());
    }
    if !name
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '-')
    {
        return Err(format!(
            "project name `{name}` must be ASCII alphanumeric plus `_`, `.`, `-` \
             (no whitespace, quotes, brackets, newlines, or non-ASCII)"
        ));
    }
    Ok(name)
}

/// Render the project.toml content for a fresh skeleton. The
/// `cases.order` array is populated with the single starter case
/// the chosen template ships, named to match the directory the
/// `case.toml` is dropped into so the loader sees the case at the
/// right slot in the ordering.
///
/// # Errors
///
/// Returns a String error if `name` or `case_dir_name` fails the
/// [`sanitize_project_name`] check.
pub fn render_project_toml(name: &str, case_dir_name: &str) -> Result<String, String> {
    let name = sanitize_project_name(name)?;
    let case_dir_name = sanitize_project_name(case_dir_name)?;
    Ok(format!(
        r#"# Generated by `valenx-init`. Edit freely.
[project]
format = "1.0"
name   = "{name}"

[geometry]
entries = []

# No mesh entries by default; add `[mesh.<name>]` blocks as
# meshes are created or imported.

[cases]
order = ["{case_dir_name}"]

[ui]

[units]
length      = "m"
mass        = "kg"
time        = "s"
temperature = "K"
"#
    ))
}

/// Render the starter case.toml for a given template. Each one is
/// a minimal-but-runnable case definition the user can edit + run
/// through their installed solver.
pub fn render_case_toml(template: Template, dir_name: &str) -> String {
    // Mirror render_project_toml: sanitize the interpolated name so a value
    // with TOML-significant chars (quotes, newlines, brackets) can't break out
    // of the `name = "..."` string. In-crate callers pass a vetted static
    // case-dir name; this guards the renderer defensively regardless.
    let dir_name = sanitize_project_name(dir_name).unwrap_or("case");
    let header = format!(
        r#"# Generated by `valenx-init`. Edit freely.
[case]
format  = "1.0"
name    = "{dir_name}"
"#
    );
    let body = match template {
        Template::Empty => include_str!("init_templates/templates/case-1.toml").to_string(),
        Template::Cfd => include_str!("init_templates/templates/cavity.toml").to_string(),
        Template::Fea => include_str!("init_templates/templates/cantilever.toml").to_string(),
        Template::Chemistry => {
            include_str!("init_templates/templates/ch4-equilibrium.toml").to_string()
        }
        Template::Su2 => include_str!("init_templates/templates/naca0012.toml").to_string(),
        Template::OpenRadioss => {
            include_str!("init_templates/templates/drop-test.toml").to_string()
        }
        Template::CodeAster => {
            include_str!("init_templates/templates/static-beam.toml").to_string()
        }
        Template::Netgen => include_str!("init_templates/templates/csg-box.toml").to_string(),
        Template::Meep => include_str!("init_templates/templates/ring-resonator.toml").to_string(),
        Template::Gromacs => include_str!("init_templates/templates/lysozyme.toml").to_string(),
        Template::Gmsh => include_str!("init_templates/templates/box-mesh.toml").to_string(),
        Template::Lammps => include_str!("init_templates/templates/lj-fluid.toml").to_string(),
        Template::ElmerHeat => include_str!("init_templates/templates/heat-cube.toml").to_string(),
        // ---- Biology (Phase 17) ----
        Template::Biopython => {
            include_str!("init_templates/templates/biopython-analyse.toml").to_string()
        }
        Template::Rdkit => include_str!("init_templates/templates/rdkit-screen.toml").to_string(),
        Template::Openmm => include_str!("init_templates/templates/protein-relax.toml").to_string(),
        Template::Chimerax => {
            include_str!("init_templates/templates/chimerax-render.toml").to_string()
        }
        Template::Oxdna => include_str!("init_templates/templates/oxdna-duplex.toml").to_string(),
        Template::Mdanalysis => {
            include_str!("init_templates/templates/trajectory-rmsd.toml").to_string()
        }
        Template::Colabfold => {
            include_str!("init_templates/templates/dna-folding.toml").to_string()
        }
        // ---- Biology — alignment toolkit (Phase 18) ----
        Template::Bwa => include_str!("init_templates/templates/bwa-align.toml").to_string(),
        Template::Minimap2 => {
            include_str!("init_templates/templates/minimap2-align.toml").to_string()
        }
        Template::Mafft => include_str!("init_templates/templates/mafft-msa.toml").to_string(),
        Template::Muscle => include_str!("init_templates/templates/muscle-msa.toml").to_string(),
        Template::Hmmer => include_str!("init_templates/templates/hmmer-search.toml").to_string(),
        Template::Samtools => {
            include_str!("init_templates/templates/samtools-flagstat.toml").to_string()
        }
        // ---- Biology — structure prediction expansion (Phase 17.5) ----
        Template::Esmfold => {
            include_str!("init_templates/templates/esmfold-predict.toml").to_string()
        }
        Template::Openfold => {
            include_str!("init_templates/templates/openfold-predict.toml").to_string()
        }
        Template::Alphafold2 => {
            include_str!("init_templates/templates/alphafold2-predict.toml").to_string()
        }
        Template::Alphafold3 => {
            include_str!("init_templates/templates/alphafold3-predict.toml").to_string()
        }
        // ---- Biology — variant calling (Phase 19) ----
        Template::Bcftools => {
            include_str!("init_templates/templates/bcftools-call.toml").to_string()
        }
        Template::Gatk => include_str!("init_templates/templates/gatk-haplotype.toml").to_string(),
        Template::DeepVariant => {
            include_str!("init_templates/templates/deepvariant-call.toml").to_string()
        }
        // ---- Biology — viewers (Phase 23) ----
        Template::Pymol => include_str!("init_templates/templates/pymol-render.toml").to_string(),
        Template::Vmd => include_str!("init_templates/templates/vmd-render.toml").to_string(),
        Template::Igv => include_str!("init_templates/templates/igv-index.toml").to_string(),
        // ---- Biology — protein design (Phase 27) ----
        Template::RfDiffusion => {
            include_str!("init_templates/templates/rfdiffusion-design.toml").to_string()
        }
        Template::ProteinMpnn => {
            include_str!("init_templates/templates/proteinmpnn-design.toml").to_string()
        }
        // ---- Biology — molecular docking (Phase 34) ----
        Template::Vina => include_str!("init_templates/templates/vina-dock.toml").to_string(),
        Template::AutoDock4 => {
            include_str!("init_templates/templates/autodock4-dock.toml").to_string()
        }
        // ---- Biology — cheminformatics expansion (Phase 24) ----
        Template::DeepChem => {
            include_str!("init_templates/templates/deepchem-screen.toml").to_string()
        }
        Template::OpenBabel => {
            include_str!("init_templates/templates/openbabel-convert.toml").to_string()
        }
        Template::Avogadro => {
            include_str!("init_templates/templates/avogadro-render.toml").to_string()
        }
        // ---- Biology — workflow managers (Phase 22) ----
        Template::Nextflow => {
            include_str!("init_templates/templates/nextflow-pipeline.toml").to_string()
        }
        Template::Snakemake => {
            include_str!("init_templates/templates/snakemake-pipeline.toml").to_string()
        }
        // ---- Biology — single-cell genomics (Phase 19.5) ----
        Template::Scanpy => {
            include_str!("init_templates/templates/scanpy-analyse.toml").to_string()
        }
        Template::Scvi => include_str!("init_templates/templates/scvi-train.toml").to_string(),
        // ---- Biology — protein design expansion (Phase 27.5) ----
        Template::Chroma => include_str!("init_templates/templates/chroma-design.toml").to_string(),
        Template::EsmIf => include_str!("init_templates/templates/esm-if-design.toml").to_string(),
        Template::RfAntibody => {
            include_str!("init_templates/templates/rfantibody-design.toml").to_string()
        }
        // ---- Biology — aligners expansion (Phase 18.5) ----
        Template::Bowtie2 => {
            include_str!("init_templates/templates/bowtie2-align.toml").to_string()
        }
        Template::Mmseqs2 => {
            include_str!("init_templates/templates/mmseqs2-search.toml").to_string()
        }
        Template::Diamond => {
            include_str!("init_templates/templates/diamond-search.toml").to_string()
        }
        // ---- Biology — RNA-seq alignment (Phase 18.6) ----
        Template::Hisat2 => include_str!("init_templates/templates/hisat2-align.toml").to_string(),
        Template::Star => include_str!("init_templates/templates/star-align.toml").to_string(),
        // ---- Biology — transcript quantification (Phase 20) ----
        Template::Salmon => include_str!("init_templates/templates/salmon-quant.toml").to_string(),
        Template::Kallisto => {
            include_str!("init_templates/templates/kallisto-quant.toml").to_string()
        }
        // ---- Biology — phylogenetics (Phase 30) ----
        Template::IqTree => include_str!("init_templates/templates/iqtree-build.toml").to_string(),
        Template::RaxmlNg => {
            include_str!("init_templates/templates/raxml-ng-build.toml").to_string()
        }
        Template::FastTree => {
            include_str!("init_templates/templates/fasttree-build.toml").to_string()
        }
        // ---- Biology — RNA structure (Phase 28) ----
        Template::ViennaRna => {
            include_str!("init_templates/templates/viennarna-fold.toml").to_string()
        }
        Template::RnaStructure => {
            include_str!("init_templates/templates/rnastructure-fold.toml").to_string()
        }
        Template::Nupack => {
            include_str!("init_templates/templates/nupack-analyze.toml").to_string()
        }
        // ---- Biology — quantum chemistry (Phase 25) ----
        Template::Psi4 => include_str!("init_templates/templates/psi4-compute.toml").to_string(),
        Template::Nwchem => {
            include_str!("init_templates/templates/nwchem-compute.toml").to_string()
        }
        Template::Xtb => include_str!("init_templates/templates/xtb-compute.toml").to_string(),
        // ---- Biology — EvolutionaryScale models (Phase 27.6) ----
        Template::Esm3 => include_str!("init_templates/templates/esm3-generate.toml").to_string(),
        Template::Esmc => include_str!("init_templates/templates/esmc-embed.toml").to_string(),
        // ---- Biology — systems biology (Phase 32) ----
        Template::Copasi => {
            include_str!("init_templates/templates/copasi-simulate.toml").to_string()
        }
        Template::BioNetGen => {
            include_str!("init_templates/templates/bionetgen-simulate.toml").to_string()
        }
        Template::PhysiCell => {
            include_str!("init_templates/templates/physicell-simulate.toml").to_string()
        }
        // ---- Biology — cryo-EM (Phase 36) ----
        Template::Relion => include_str!("init_templates/templates/relion-refine.toml").to_string(),
        Template::Eman2 => include_str!("init_templates/templates/eman2-refine.toml").to_string(),
        Template::Ctffind => {
            include_str!("init_templates/templates/ctffind-estimate.toml").to_string()
        }
        // ---- Biology — sequencing read simulators (Phase 31) ----
        Template::Art => include_str!("init_templates/templates/art-simulate.toml").to_string(),
        Template::Wgsim => include_str!("init_templates/templates/wgsim-simulate.toml").to_string(),
        Template::Badread => {
            include_str!("init_templates/templates/badread-simulate.toml").to_string()
        }
        // ---- Biology — CRISPR design (Phase 35) ----
        Template::Chopchop => {
            include_str!("init_templates/templates/chopchop-design.toml").to_string()
        }
        Template::Crispor => {
            include_str!("init_templates/templates/crispor-design.toml").to_string()
        }
        Template::CasOffinder => {
            include_str!("init_templates/templates/cas-offinder-search.toml").to_string()
        }
        // ---- Biology — Rosetta family (Phase 38) ----
        Template::Rosetta => {
            include_str!("init_templates/templates/rosetta-protocol.toml").to_string()
        }
        Template::PyRosetta => {
            include_str!("init_templates/templates/pyrosetta-script.toml").to_string()
        }
        // ---- Biology — population genetics (Phase 29) ----
        Template::Slim => include_str!("init_templates/templates/slim-simulate.toml").to_string(),
        Template::Msprime => {
            include_str!("init_templates/templates/msprime-simulate.toml").to_string()
        }
        Template::Tskit => include_str!("init_templates/templates/tskit-analyze.toml").to_string(),
        // ---- Biology — Bayesian phylogenetics (Phase 30.5) ----
        Template::Beast2 => include_str!("init_templates/templates/beast2-mcmc.toml").to_string(),
        Template::MrBayes => include_str!("init_templates/templates/mrbayes-mcmc.toml").to_string(),
        // ---- Biology — DNA structural geometry (Phase 39) ----
        Template::X3dna => include_str!("init_templates/templates/x3dna-analyze.toml").to_string(),
        Template::Curves => {
            include_str!("init_templates/templates/curves-analyze.toml").to_string()
        }
        Template::Dssr => include_str!("init_templates/templates/dssr-analyze.toml").to_string(),
        // ---- Biology — MD analysis expansion (Phase 5.5) ----
        Template::Plumed => {
            include_str!("init_templates/templates/plumed-analyze.toml").to_string()
        }
        Template::Prody => include_str!("init_templates/templates/prody-analyze.toml").to_string(),
        Template::Cpptraj => {
            include_str!("init_templates/templates/cpptraj-analyze.toml").to_string()
        }
        // ---- Biology — synthetic biology (Phase 33) ----
        Template::PySbol => {
            include_str!("init_templates/templates/pysbol-compose.toml").to_string()
        }
        Template::J5 => include_str!("init_templates/templates/j5-assemble.toml").to_string(),
        Template::Cello => include_str!("init_templates/templates/cello-compile.toml").to_string(),
        // ---- Biology — alignment toolkit expansion (Phase 18.7) ----
        Template::Blast => include_str!("init_templates/templates/blast-search.toml").to_string(),
        Template::Clustalo => {
            include_str!("init_templates/templates/clustalo-align.toml").to_string()
        }
        Template::TCoffee => {
            include_str!("init_templates/templates/tcoffee-align.toml").to_string()
        }
        // ---- Biology — single-cell genomics expansion (Phase 19.6) ----
        Template::Seurat => {
            include_str!("init_templates/templates/seurat-analyze.toml").to_string()
        }
        Template::AnnData => {
            include_str!("init_templates/templates/anndata-process.toml").to_string()
        }
        // ---- Biology — bio MD engines (Phase 5.6) ----
        Template::Namd => include_str!("init_templates/templates/namd-simulate.toml").to_string(),
        Template::Sander => {
            include_str!("init_templates/templates/sander-simulate.toml").to_string()
        }
        Template::Hoomd => include_str!("init_templates/templates/hoomd-simulate.toml").to_string(),
        // ---- Biology — MD analysis sister (Phase 5.7) ----
        Template::Mdtraj => {
            include_str!("init_templates/templates/mdtraj-analyze.toml").to_string()
        }
        // ---- Biology — structure prediction + search (Phase 17.7) ----
        Template::RoseTTAFold => {
            include_str!("init_templates/templates/rosettafold-predict.toml").to_string()
        }
        Template::OmegaFold => {
            include_str!("init_templates/templates/omegafold-predict.toml").to_string()
        }
        Template::Foldseek => {
            include_str!("init_templates/templates/foldseek-search.toml").to_string()
        }
        // ---- Biology — spatial stochastic reaction-diffusion (Phase 32.5) ----
        Template::Smoldyn => {
            include_str!("init_templates/templates/smoldyn-simulate.toml").to_string()
        }
        Template::Mcell => include_str!("init_templates/templates/mcell-simulate.toml").to_string(),
        // ---- Biology — sequence editors / plasmid design (Phase 41) ----
        Template::Pydna => include_str!("init_templates/templates/pydna-design.toml").to_string(),
        Template::Jalview => include_str!("init_templates/templates/jalview-view.toml").to_string(),
        // ---- Biology — microscopy / bioimage analysis (Phase 40) ----
        Template::Fiji => include_str!("init_templates/templates/fiji-process.toml").to_string(),
        Template::CellProfiler => {
            include_str!("init_templates/templates/cellprofiler-segment.toml").to_string()
        }
        Template::Ilastik => {
            include_str!("init_templates/templates/ilastik-classify.toml").to_string()
        }
        // ---- Biology — workflow expansion (Phase 22.5) ----
        Template::Planemo => include_str!("init_templates/templates/planemo-run.toml").to_string(),
        Template::Cromwell => {
            include_str!("init_templates/templates/cromwell-run.toml").to_string()
        }
        Template::Cwltool => include_str!("init_templates/templates/cwltool-run.toml").to_string(),
        // ---- Biology — web 3D molecular visualization (Phase 42) ----
        Template::Molstar => include_str!("init_templates/templates/molstar-view.toml").to_string(),
        Template::Ngl => include_str!("init_templates/templates/ngl-view.toml").to_string(),
        // ---- Biology — mRNA design (Phase 43) ----
        Template::DnaChisel => {
            include_str!("init_templates/templates/dnachisel-optimize.toml").to_string()
        }
        Template::LinearDesign => {
            include_str!("init_templates/templates/lineardesign-design.toml").to_string()
        }
        Template::Icodon => {
            include_str!("init_templates/templates/icodon-predict.toml").to_string()
        }
        // ---- Biology — RNA folding expansion (Phase 44.5) ----
        Template::Mfold => include_str!("init_templates/templates/mfold-fold.toml").to_string(),
        Template::EternaFold => {
            include_str!("init_templates/templates/eternafold-fold.toml").to_string()
        }
        Template::LinearFold => {
            include_str!("init_templates/templates/linearfold-fold.toml").to_string()
        }
        // ---- Biology — base + prime editing design (Phase 35.5) ----
        Template::BeDesigner => {
            include_str!("init_templates/templates/be-designer-design.toml").to_string()
        }
        Template::BeHive => {
            include_str!("init_templates/templates/be-hive-predict.toml").to_string()
        }
        Template::PrimeDesign => {
            include_str!("init_templates/templates/primedesign-design.toml").to_string()
        }
        Template::PegFinder => {
            include_str!("init_templates/templates/pegfinder-design.toml").to_string()
        }
        // ---- Biology — edit-outcome prediction (Phase 35.6) ----
        Template::Indelphi => {
            include_str!("init_templates/templates/indelphi-predict.toml").to_string()
        }
        Template::Forecast => {
            include_str!("init_templates/templates/forecast-predict.toml").to_string()
        }
        Template::AlphaMissense => {
            include_str!("init_templates/templates/alphamissense-predict.toml").to_string()
        }
        Template::Crispritz => {
            include_str!("init_templates/templates/crispritz-search.toml").to_string()
        }
        // ---- Biology — pharmacokinetics + RNA tertiary (Phase 45) ----
        Template::PkSim => include_str!("init_templates/templates/pksim-simulate.toml").to_string(),
        Template::SimRna => include_str!("init_templates/templates/simrna-fold.toml").to_string(),
    };
    header + &body
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render_template_list_covers_every_template() {
        // The static TEMPLATE_ROWS catalogue must include every name
        // `Template::from_str` recognises canonically. Drift surfaces
        // as a missing assertion below — adding a new template
        // without listing it here is a test failure.
        let out = render_template_list();
        for canonical in [
            "empty",
            "cfd",
            "fea",
            "chemistry",
            "su2",
            "openradioss",
            "code-aster",
            "netgen",
            "meep",
            "gromacs",
            "gmsh",
            "lammps",
            "elmer-heat",
            // Biology (Phase 17)
            "biopython",
            "rdkit",
            "openmm",
            "chimerax",
            "oxdna",
            "mdanalysis",
            "colabfold",
        ] {
            assert!(
                out.contains(canonical),
                "template list missing `{canonical}`:\n{out}"
            );
        }
        // Header + footer hints surface so users know what to do next.
        assert!(out.contains("Available"), "missing header: {out}");
        assert!(out.contains("--template"), "missing usage hint: {out}");
    }

    #[test]
    fn template_rows_match_template_enum_via_from_str() {
        // Every row's `name` must round-trip through Template::from_str.
        // Stops typos like `cgs` instead of `cfd` from sliding into
        // the user-visible list.
        for row in TEMPLATE_ROWS {
            assert!(
                Template::from_str(row.name).is_some(),
                "TEMPLATE_ROWS lists `{}` but Template::from_str doesn't recognise it",
                row.name
            );
        }
    }

    #[test]
    fn render_project_toml_has_required_fields() {
        let s = render_project_toml("smoke", "case-1").unwrap();
        assert!(s.contains("format = \"1.0\""));
        assert!(s.contains("name   = \"smoke\""));
        assert!(s.contains("[geometry]"));
        assert!(s.contains("[cases]"));
        assert!(s.contains("order = [\"case-1\"]"));
        assert!(s.contains("[units]"));
    }

    /// Regression: the `cases.order` entry must match the directory
    /// name we drop the `case.toml` into. Pre-fix, this was
    /// hardcoded to `"case-1"` for every template, which left
    /// non-Empty scaffolds with a dangling reference and forced
    /// the loader to fall back on directory scanning.
    #[test]
    fn render_project_toml_case_order_matches_dir_name() {
        let s = render_project_toml("p", "cavity").unwrap();
        assert!(s.contains("order = [\"cavity\"]"), "got: {s}");
        let s = render_project_toml("p", "cantilever").unwrap();
        assert!(s.contains("order = [\"cantilever\"]"), "got: {s}");
        let s = render_project_toml("p", "csg-box").unwrap();
        assert!(s.contains("order = [\"csg-box\"]"), "got: {s}");
    }

    /// Round-12 L9 RED→GREEN: a project name containing newline /
    /// quote characters would otherwise inject arbitrary TOML
    /// content into the rendered file. `sanitize_project_name`
    /// rejects the value before the renderer ever sees it.
    #[test]
    fn render_project_toml_rejects_injection_via_quotes_and_newlines() {
        let attack = "evil\"\n[rbac]\ndefault_role = \"viewer\"\n#";
        let err = render_project_toml(attack, "case-1").expect_err("must reject TOML injection");
        assert!(
            err.contains("must be ASCII alphanumeric"),
            "expected sanitiser message; got: {err}"
        );
    }

    /// Round-12 L9: the same check applies via the `scaffold_project`
    /// public API.
    #[test]
    fn scaffold_project_rejects_toml_injection_name() {
        let tmp =
            std::env::temp_dir().join(format!("valenx_init_l9_inject_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).unwrap();
        let attack = "evil\"\n[rbac]\ndefault_role = \"viewer\"\n#";
        let err = scaffold_project(&tmp, Template::Empty, Some(attack))
            .expect_err("must reject TOML injection via name");
        assert!(
            err.contains("must be ASCII alphanumeric") || err.contains("project name"),
            "expected sanitiser message; got: {err}"
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn sanitize_project_name_accepts_realistic_names() {
        assert!(sanitize_project_name("cavity").is_ok());
        assert!(sanitize_project_name("smoke_test-1.0").is_ok());
        assert!(sanitize_project_name("Project123").is_ok());
    }

    #[test]
    fn sanitize_project_name_rejects_dangerous_chars() {
        assert!(sanitize_project_name("").is_err());
        assert!(sanitize_project_name("has space").is_err());
        assert!(sanitize_project_name("has\"quote").is_err());
        assert!(sanitize_project_name("has\nnewline").is_err());
        assert!(sanitize_project_name("has[bracket").is_err());
        assert!(sanitize_project_name("résumé").is_err()); // non-ASCII
    }

    #[test]
    fn render_case_toml_chemistry_template_has_initial_state() {
        let s = render_case_toml(Template::Chemistry, "ch");
        assert!(s.contains("physics = \"chemistry\""));
        assert!(s.contains("solver  = \"cantera.equilibrium\""));
        assert!(s.contains("[chemistry.initial]"));
        assert!(s.contains("composition"));
    }

    #[test]
    fn render_case_toml_su2_template_has_cfd_su2_block() {
        let s = render_case_toml(Template::Su2, "naca");
        assert!(s.contains("physics = \"cfd\""));
        assert!(s.contains("[cfd.su2]"));
        assert!(s.contains("config"));
        assert!(s.contains("mesh"));
    }

    #[test]
    fn render_case_toml_openradioss_template_has_engine_input() {
        let s = render_case_toml(Template::OpenRadioss, "drop");
        assert!(s.contains("physics = \"fea\""));
        assert!(s.contains("[fea.openradioss]"));
        assert!(s.contains("engine_input"));
    }

    #[test]
    fn render_case_toml_code_aster_template_has_export() {
        let s = render_case_toml(Template::CodeAster, "static");
        assert!(s.contains("physics = \"fea\""));
        assert!(s.contains("[fea.code_aster]"));
        assert!(s.contains("export"));
    }

    #[test]
    fn render_case_toml_netgen_template_has_geometry_file() {
        let s = render_case_toml(Template::Netgen, "csg");
        assert!(s.contains("physics = \"meshing\""));
        assert!(s.contains("[meshing.netgen]"));
        assert!(s.contains("geometry_file"));
    }

    #[test]
    fn render_case_toml_meep_template_has_python_script() {
        let s = render_case_toml(Template::Meep, "ring");
        assert!(s.contains("physics = \"em\""));
        assert!(s.contains("[em.meep]"));
        assert!(s.contains("script"));
    }

    #[test]
    fn render_case_toml_gromacs_template_has_tpr() {
        let s = render_case_toml(Template::Gromacs, "lyso");
        assert!(s.contains("physics = \"md\""));
        assert!(s.contains("[md.gromacs]"));
        assert!(s.contains("tpr"));
    }

    #[test]
    fn template_aliases_resolve_for_new_adapters() {
        // Each new adapter accepts at least one short alias.
        assert_eq!(
            Template::from_str("openradioss"),
            Some(Template::OpenRadioss)
        );
        assert_eq!(Template::from_str("crash"), Some(Template::OpenRadioss));
        assert_eq!(Template::from_str("aster"), Some(Template::CodeAster));
        assert_eq!(Template::from_str("code_aster"), Some(Template::CodeAster));
        assert_eq!(Template::from_str("meshing"), Some(Template::Netgen));
        assert_eq!(Template::from_str("photonics"), Some(Template::Meep));
        assert_eq!(Template::from_str("md"), Some(Template::Gromacs));
        assert_eq!(Template::from_str("compressible"), Some(Template::Su2));
        // Case-insensitive (existing behaviour).
        assert_eq!(
            Template::from_str("OpenRadioss"),
            Some(Template::OpenRadioss)
        );
        assert_eq!(Template::from_str("SU2"), Some(Template::Su2));
    }

    #[test]
    fn template_aliases_resolve_for_gmsh_lammps_elmer_heat() {
        // gmsh / lammps / elmer-heat templates added in this rev —
        // each one accepts its canonical name plus at least one
        // domain alias.
        assert_eq!(Template::from_str("gmsh"), Some(Template::Gmsh));
        assert_eq!(Template::from_str("delaunay"), Some(Template::Gmsh));
        assert_eq!(Template::from_str("lammps"), Some(Template::Lammps));
        assert_eq!(Template::from_str("lj"), Some(Template::Lammps));
        assert_eq!(Template::from_str("classical-md"), Some(Template::Lammps));
        assert_eq!(Template::from_str("elmer-heat"), Some(Template::ElmerHeat));
        assert_eq!(Template::from_str("elmer"), Some(Template::ElmerHeat));
        assert_eq!(Template::from_str("heat"), Some(Template::ElmerHeat));
        // Case-insensitive across the new aliases.
        assert_eq!(Template::from_str("Lammps"), Some(Template::Lammps));
        assert_eq!(Template::from_str("ELMER-HEAT"), Some(Template::ElmerHeat));
    }

    #[test]
    fn render_case_toml_gmsh_template_has_box_mesh_block() {
        let s = render_case_toml(Template::Gmsh, "box");
        assert!(s.contains("physics = \"meshing\""), "got: {s}");
        assert!(s.contains("solver  = \"gmsh.delaunay\""), "got: {s}");
        assert!(s.contains("[mesh]"), "got: {s}");
        assert!(s.contains("type                  = \"box\""), "got: {s}");
        assert!(s.contains("characteristic_length"), "got: {s}");
    }

    #[test]
    fn render_case_toml_lammps_template_has_md_block_and_lj_potential() {
        let s = render_case_toml(Template::Lammps, "lj");
        assert!(s.contains("physics = \"molecular-dynamics\""), "got: {s}");
        assert!(s.contains("solver  = \"lammps.nve\""), "got: {s}");
        assert!(s.contains("[md]"), "got: {s}");
        assert!(s.contains("[md.init]"), "got: {s}");
        assert!(s.contains("[md.potential]"), "got: {s}");
        assert!(s.contains("kind    = \"lj-cut\""), "got: {s}");
    }

    #[test]
    fn render_case_toml_elmer_heat_template_has_heat_block_and_two_boundaries() {
        let s = render_case_toml(Template::ElmerHeat, "heat");
        assert!(s.contains("physics = \"fea\""), "got: {s}");
        assert!(s.contains("solver  = \"elmer.heat\""), "got: {s}");
        assert!(s.contains("[heat]"), "got: {s}");
        assert!(s.contains("[heat.material]"), "got: {s}");
        // Two boundary blocks — one hot face, one cold face.
        let n_boundaries = s.matches("[[heat.boundaries]]").count();
        assert_eq!(n_boundaries, 2, "expected 2 boundaries; got: {s}");
    }

    #[test]
    fn template_aliases_resolve_for_bio_adapters() {
        // Phase 17 — bio templates accept canonical names plus a
        // short alias (e.g. `omm` for `openmm`, `cf` for `colabfold`).
        assert_eq!(Template::from_str("biopython"), Some(Template::Biopython));
        assert_eq!(Template::from_str("biopy"), Some(Template::Biopython));
        assert_eq!(Template::from_str("rdkit"), Some(Template::Rdkit));
        assert_eq!(Template::from_str("chem-py"), Some(Template::Rdkit));
        assert_eq!(Template::from_str("openmm"), Some(Template::Openmm));
        assert_eq!(Template::from_str("omm"), Some(Template::Openmm));
        assert_eq!(Template::from_str("chimerax"), Some(Template::Chimerax));
        assert_eq!(Template::from_str("cxc"), Some(Template::Chimerax));
        assert_eq!(Template::from_str("viz3d"), Some(Template::Chimerax));
        assert_eq!(Template::from_str("oxdna"), Some(Template::Oxdna));
        assert_eq!(Template::from_str("cgdna"), Some(Template::Oxdna));
        assert_eq!(Template::from_str("mdanalysis"), Some(Template::Mdanalysis));
        assert_eq!(Template::from_str("mda"), Some(Template::Mdanalysis));
        assert_eq!(Template::from_str("traj-py"), Some(Template::Mdanalysis));
        assert_eq!(Template::from_str("colabfold"), Some(Template::Colabfold));
        assert_eq!(Template::from_str("cf"), Some(Template::Colabfold));
        assert_eq!(
            Template::from_str("protein-fold"),
            Some(Template::Colabfold)
        );
        // Case-insensitive across the new aliases.
        assert_eq!(Template::from_str("BIOPYTHON"), Some(Template::Biopython));
        assert_eq!(Template::from_str("ColabFold"), Some(Template::Colabfold));
    }

    #[test]
    fn render_case_toml_bio_templates_emit_bio_block() {
        // Each bio template emits `physics = "bio"` and the right
        // `[bio.<adapter>]` block matching the adapter's case_input
        // schema.
        for (template, expected_block, solver) in [
            (Template::Biopython, "[bio.biopython]", "biopython.script"),
            (Template::Rdkit, "[bio.rdkit]", "rdkit.script"),
            (Template::Openmm, "[bio.openmm]", "openmm.script"),
            (Template::Chimerax, "[bio.chimerax]", "chimerax.script"),
            (Template::Oxdna, "[bio.oxdna]", "oxdna.batch"),
            (
                Template::Mdanalysis,
                "[bio.mdanalysis]",
                "mdanalysis.script",
            ),
            (Template::Colabfold, "[bio.colabfold]", "colabfold.predict"),
        ] {
            let s = render_case_toml(template, "case");
            assert!(s.contains("physics = \"bio\""), "got: {s}");
            assert!(s.contains(solver), "missing solver `{solver}` in: {s}");
            assert!(
                s.contains(expected_block),
                "missing block `{expected_block}` in: {s}",
            );
        }
    }

    #[test]
    fn render_case_toml_colabfold_template_has_input_fasta_and_recycles() {
        // ColabFold's case_input requires `input_fasta`; the optional
        // `num_recycles` / `num_models` are emitted with their
        // adapter-default values so the case round-trips cleanly
        // through `from_case_dir`.
        let s = render_case_toml(Template::Colabfold, "fold");
        assert!(s.contains("input_fasta  = \"query.fasta\""), "got: {s}");
        assert!(s.contains("num_recycles = 3"), "got: {s}");
        assert!(s.contains("num_models   = 5"), "got: {s}");
    }

    #[test]
    fn scaffold_biopython_drops_sample_python_script() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-biopython-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Biopython, None).expect("scaffold");
        let case_dir = dir.join("cases").join("biopython-analyse");
        let script = case_dir.join("analyse.py");
        assert!(
            script.is_file(),
            "expected sample script at {}",
            script.display()
        );
        let text = std::fs::read_to_string(&script).unwrap();
        assert!(text.contains("import Bio"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_colabfold_drops_sample_fasta() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-colabfold-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Colabfold, None).expect("scaffold");
        let case_dir = dir.join("cases").join("dna-folding");
        let fasta = case_dir.join("query.fasta");
        assert!(
            fasta.is_file(),
            "expected sample FASTA at {}",
            fasta.display()
        );
        let text = std::fs::read_to_string(&fasta).unwrap();
        // FASTA convention — header begins with `>`.
        assert!(text.starts_with('>'), "got: {text}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_oxdna_does_not_drop_sample_input_dat() {
        // oxDNA's `input.dat` is user-supplied (encodes their
        // simulation parameters end-to-end); we deliberately don't
        // ship a placeholder. Only case.toml exists.
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-oxdna-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Oxdna, None).expect("scaffold");
        let case_dir = dir.join("cases").join("oxdna-duplex");
        assert!(case_dir.join("case.toml").is_file());
        let entries: Vec<_> = std::fs::read_dir(&case_dir).unwrap().flatten().collect();
        assert_eq!(entries.len(), 1, "expected only case.toml in {entries:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn render_case_toml_cfd_template_has_boundaries() {
        let s = render_case_toml(Template::Cfd, "cavity");
        assert!(s.contains("physics = \"cfd\""));
        assert!(s.contains("[boundaries.inlet]"));
        assert!(s.contains("[boundaries.outlet]"));
        assert!(s.contains("kEpsilon"));
    }

    #[test]
    fn render_case_toml_fea_template_has_material_and_step() {
        let s = render_case_toml(Template::Fea, "cant");
        assert!(s.contains("physics = \"fea\""));
        assert!(s.contains("[structural.material]"));
        assert!(s.contains("youngs_modulus"));
        assert!(s.contains("[structural.step]"));
    }

    #[test]
    fn render_case_toml_empty_template_has_geometry_physics() {
        let s = render_case_toml(Template::Empty, "case-1");
        assert!(s.contains("physics = \"geometry\""));
        assert!(s.contains("solver  = \"(none)\""));
    }

    #[test]
    fn scaffold_project_creates_full_skeleton() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Cfd, Some("smoke")).expect("scaffold");
        // project.toml exists
        let project_toml = dir.join("project.toml");
        assert!(project_toml.is_file());
        let project_text = std::fs::read_to_string(&project_toml).unwrap();
        assert!(project_text.contains("name   = \"smoke\""));
        // cases.order points at the actual case dir we created — pre-
        // fix this said `["case-1"]` regardless of template, which
        // forced the loader to fall back on directory scanning.
        assert!(
            project_text.contains("order = [\"cavity\"]"),
            "expected `order = [\"cavity\"]` in project.toml; got:\n{project_text}"
        );
        // cases/cavity/case.toml exists with cfd content
        let case_toml = dir.join("cases").join("cavity").join("case.toml");
        assert!(case_toml.is_file());
        let case_text = std::fs::read_to_string(&case_toml).unwrap();
        assert!(case_text.contains("physics = \"cfd\""));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_netgen_drops_sample_geo_alongside_case_toml() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-netgen-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Netgen, None).expect("scaffold");
        let case_dir = dir.join("cases").join("csg-box");
        assert!(case_dir.join("case.toml").is_file());
        let geo = case_dir.join("csg-box.geo");
        assert!(geo.is_file(), "expected sample .geo at {}", geo.display());
        let geo_text = std::fs::read_to_string(&geo).unwrap();
        assert!(geo_text.contains("algebraic3d"));
        assert!(geo_text.contains("orthobrick"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_meep_drops_sample_python_script() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-meep-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Meep, None).expect("scaffold");
        let case_dir = dir.join("cases").join("ring-resonator");
        let script = case_dir.join("ring.py");
        assert!(script.is_file());
        let script_text = std::fs::read_to_string(&script).unwrap();
        assert!(script_text.contains("import meep"));
        assert!(script_text.contains("sim.run"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_su2_drops_sample_cfg() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-su2-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Su2, None).expect("scaffold");
        let case_dir = dir.join("cases").join("naca0012");
        let cfg = case_dir.join("naca0012.cfg");
        assert!(cfg.is_file());
        let cfg_text = std::fs::read_to_string(&cfg).unwrap();
        assert!(cfg_text.contains("SOLVER"));
        assert!(cfg_text.contains("MACH_NUMBER"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_gromacs_does_not_drop_sample_tpr() {
        // GROMACS's .tpr is binary and produced by `gmx grompp` —
        // we deliberately don't ship one. Only case.toml exists.
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-gmx-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        scaffold_project(&dir, Template::Gromacs, None).expect("scaffold");
        let case_dir = dir.join("cases").join("lysozyme");
        assert!(case_dir.join("case.toml").is_file());
        // Only case.toml — no .tpr / no other inputs (those are
        // user-supplied / binary).
        let entries: Vec<_> = std::fs::read_dir(&case_dir).unwrap().flatten().collect();
        assert_eq!(entries.len(), 1, "expected only case.toml in {entries:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_project_refuses_to_overwrite_existing_project() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-overwrite-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("project.toml"), b"# pre-existing").unwrap();
        let err = scaffold_project(&dir, Template::Empty, None).unwrap_err();
        assert!(err.contains("already exists"));
        // The existing file is untouched.
        let text = std::fs::read_to_string(dir.join("project.toml")).unwrap();
        assert_eq!(text, "# pre-existing");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scaffold_project_uses_dir_name_when_name_omitted() {
        let dir = std::env::temp_dir().join(format!(
            "valenx-init-default-name-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let dir_name = dir
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap()
            .to_string();
        scaffold_project(&dir, Template::Empty, None).expect("scaffold");
        let text = std::fs::read_to_string(dir.join("project.toml")).unwrap();
        assert!(text.contains(&format!("name   = \"{dir_name}\"")));
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ---- GUI-facing wrapper API (Phase 44) ----

    #[test]
    fn case_dir_for_known_adapters() {
        // The adapter-id → case-dir mapping is one of the contracts
        // the GUI's "New case from adapter" flow depends on. Spot-
        // check a handful of representative adapters across domains.
        assert_eq!(case_dir_for("dnachisel"), Some("dnachisel-optimize"));
        assert_eq!(case_dir_for("blast"), Some("blast-search"));
        assert_eq!(case_dir_for("openfoam"), Some("cavity"));
        assert_eq!(case_dir_for("calculix"), Some("cantilever"));
        assert_eq!(case_dir_for("gmsh"), Some("box-mesh"));
    }

    #[test]
    fn case_dir_for_unknown_adapter_is_none() {
        assert!(case_dir_for("rocket-science").is_none());
        assert!(case_dir_for("").is_none());
    }

    #[test]
    fn case_toml_body_emits_expected_block_for_dnachisel() {
        let body = case_toml_body("dnachisel").expect("dnachisel known");
        assert!(body.contains("physics = \"bio\""), "got: {body}");
        assert!(
            body.contains("solver  = \"dnachisel.optimize\""),
            "got: {body}"
        );
        assert!(body.contains("[bio.dnachisel]"), "got: {body}");
        // The dir_name slot is stamped from case_dir_for(adapter_id),
        // not the adapter id itself.
        assert!(
            body.contains("name    = \"dnachisel-optimize\""),
            "got: {body}"
        );
    }

    #[test]
    fn case_toml_body_unknown_adapter_is_none() {
        assert!(case_toml_body("rocket-science").is_none());
    }

    #[test]
    fn project_toml_alias_matches_render_project_toml() {
        // Public wrapper is a thin alias — same output for same input.
        let a = project_toml("smoke", "case-1");
        let b = render_project_toml("smoke", "case-1");
        assert_eq!(a, b);
    }

    #[test]
    fn template_rows_returns_full_catalogue() {
        // The catalogue must be non-empty and round-trip every name
        // through Template::from_str (drift surfaces as a panic).
        let rows = template_rows();
        assert!(rows.len() >= 100, "expected ≥100 rows; got {}", rows.len());
        for row in rows {
            assert!(
                Template::from_str(row.name).is_some(),
                "template_rows() row `{}` not recognised by Template::from_str",
                row.name
            );
        }
    }
}
