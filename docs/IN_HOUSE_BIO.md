# In-House Bio/DNA: Native-Rust Zero-Download Adapters

**Branch:** `cloud/in-house-bio`  
**Goal:** every bio/DNA feature works with zero external tools installed.

---

## 1. The Native-Rust Core (Already Built)

These crates are fully native Rust — pure algorithms, no external processes,
no model weights. They replace the listed external tools completely:

| Crate | Replaces | What it does |
|-------|----------|--------------|
| `valenx-bioseq` | Biopython, BioJava, SeqKit, ApE, SerialCloner, pLannotate, pydna | FASTA/FASTQ/GenBank/EMBL I/O, IUPAC alphabets, translation (25 NCBI tables), 6-frame ORF finding, reverse-complement, GC/kmer/Tm/weight analysis, restriction digest (REBASE subset), codon optimisation (CAI), plasmid annotation, primer design with hairpin/dimer screening, in-silico PCR |
| `valenx-align` | BLAST+, BWA, Bowtie2, minimap2, ClustalΩ, MUSCLE, MAFFT, T-Coffee, HMMER, DIAMOND, MMseqs2 | BLOSUM/PAM/NUC4.4 scoring, Needleman-Wunsch, Gotoh affine, Smith-Waterman, semi-global, banded, Hirschberg; k-mer/seed-and-extend with Karlin-Altschul E-values; FM-index exact search; minimizer sketches; anchor chaining; progressive+iterative MSA; Plan7-style profile HMM; PSSM; SAM/Clustal/Stockholm/PHYLIP/MSF I/O; edit distance |
| `valenx-phylo` | BEAST 2, MrBayes, IQ-TREE, RAxML-NG, PhyML, FastTree, RevBayes, Seq-Gen, FigTree, Dendroscope | Newick/NEXUS/PhyloXML I/O; JC69/K80/F81/HKY85/GTR models; Felsenstein pruning; discrete-gamma rates; NNI+SPR ML topology search; Bayesian MCMC (MH, NNI/SPR/Wilson-Balding proposals, GTR, Gelman-Rubin); bootstrap; RF/quartet distances; UPGMA, WPGMA, NJ, BIONJ; Fitch/Sankoff parsimony; coalescent + birth-death simulation |
| `valenx-popgen` | SLiM, msprime, tskit, fwdpy11, simuPOP, Nemo, discoal, `ms`, stdpopsim | Kingman coalescent, structured coalescent, ARG, succinct tree-sequence; diploid Wright-Fisher forward sim with selection/recombination/migration/demography; SFS, π, θ, Tajima's D, Fst (WC), LD, EHH/iHS; ABC inference; stdpopsim-class species catalog |
| `valenx-rnastruct` | ViennaRNA, RNAstructure, mfold/UNAFold, NUPACK, ContraFold, IPknot | Full Turner-2004 nearest-neighbor parameters; Zuker MFE (`-d2` coaxial exact); LinearFold beam-search for long RNA; McCaskill partition function + BPP matrix; LinearPartition; centroid/MEA structures; Zuker suboptimal; Boltzmann sampling; Kinfold kinetics; RNA-RNA cofold; IntaRNA-class interaction DP; inverse folding; tRNA cloverleaf; pseudoknot (H-type + kissing-hairpin) |
| `valenx-cheminf` | RDKit, OpenBabel, Avogadro, DeepChem (non-ML parts) | MOL/SMILES/InChI I/O; Morgan/MACCS/path fingerprints; Tanimoto/Dice similarity; scaffold decomposition; MMFF94 forcefield; 3-D embedding (ETKDG); 2-D layout; SMARTS matching; reaction transforms; tautomer enumeration; QED drug-likeness; pharmacophore matching |
| `valenx-biostruct` | PDB/mmCIF viewers, DSSP, PyMOL structure tools (non-vis parts) | PDB/mmCIF read+write; DSSP secondary structure; SASA; Ramachandran; contact maps; superposition; base-pair geometry; nucleic helical parameters |
| `valenx-genomics` | Many GATK/samtools-equivalent variant + assembly ops | De Bruijn + OLC assembly; variant calling/filtering/genotype/normalise; CRISPR guide + off-target; read trim/filter/dedup/coverage/QC; Illumina + long-read simulation; GFF/BED/VCF/pileup/SAM I/O |
| `valenx-genediting` | CRISPOR, CRISPRitz, Cas-OFFinder, InDelPhi (guide design part), prime/base editing tools | CRISPR knockout/donor/multiplex/guide design; prime editing (pegRNA strategy); base editing (CBE/ABE design); mRNA construct design; delivery/safety analysis |
| `valenx-rnadesign` | LinearDesign (pure algorithm part), iCodon, DNA Chisel (design algorithms) | Inverse folding; coding/regulatory/riboswitch/aptamer RNA design; multistate design; codon + UTR optimisation |

**These crates have zero external dependencies and work offline with no downloads.
The problem is that the external-tool ADAPTERS still shell out — this document
tracks the plan to wire those adapters to the native crates as the default path.**

---

## 2. External-Tool Adapters: Full Inventory

`grep`-based discovery sources: `Command::new`, `subprocess::run`,
`find_on_path`, `PROBE_BINARIES`, adapter `lib.rs` in
`crates/valenx-adapters/bio/`.

### 2A — Category A: Native crate exists; adapter still shells out
(fix = wire adapter to native crate as default)

| Adapter | External binary/tool | What it wraps | Native crate | Status in this branch |
|---------|---------------------|---------------|--------------|----------------------|
| `valenx-adapter-blast` | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` | Local sequence search, E-value statistics | `valenx-align` (seed-and-extend + KA stats) | **NATIVE PATH ADDED** |
| `valenx-adapter-hmmer` | `hmmsearch`, `hmmscan` | Profile-HMM sequence search | `valenx-align` (Plan7 profile HMM) | **NATIVE PATH ADDED** |
| `valenx-adapter-muscle` | `muscle` | Multiple-sequence alignment | `valenx-align` (progressive + iterative MSA) | **NATIVE PATH ADDED** |
| `valenx-adapter-mafft` | `mafft` | Multiple-sequence alignment | `valenx-align` (progressive + iterative MSA) | **NATIVE PATH ADDED** |
| `valenx-adapter-clustalo` | `clustalo` | Multiple-sequence alignment | `valenx-align` (progressive + iterative MSA) | **NATIVE PATH ADDED** |
| `valenx-adapter-fasttree` | `FastTree`, `fasttree` | Approximate-ML phylogeny | `valenx-phylo` (NJ + NNI/SPR ML) | **NATIVE PATH ADDED** |
| `valenx-adapter-viennarna` | `RNAfold` | RNA MFE folding (Turner model) | `valenx-rnastruct` (Zuker + Turner-2004) | **NATIVE PATH ADDED** |
| `valenx-adapter-bowtie2` | `bowtie2`, `bowtie2-build` | Short-read mapping | `valenx-align` (FM-index + seed-extend) | Roadmap B-1 |
| `valenx-adapter-bwa` | `bwa` | Short-read mapping (BWA-MEM) | `valenx-align` (FM-index + seed-extend) | Roadmap B-2 |
| `valenx-adapter-minimap2` | `minimap2` | Long-read + splice mapping | `valenx-align` (minimizer + chain) | Roadmap B-3 |
| `valenx-adapter-samtools` | `samtools` | SAM/BAM manipulation | `valenx-genomics` (SAM I/O + ops) | Roadmap B-4 |
| `valenx-adapter-bcftools` | `bcftools` | VCF/BCF manipulation | `valenx-genomics` (VCF ops) | Roadmap B-5 |
| `valenx-adapter-diamond` | `diamond` | Accelerated protein search | `valenx-align` (seed-extend, protein) | Roadmap B-6 |
| `valenx-adapter-mmseqs2` | `mmseqs` | Fast protein/nucleotide search | `valenx-align` (k-mer index + search) | Roadmap B-7 |
| `valenx-adapter-iqtree` | `iqtree2` | ML phylogeny + model selection | `valenx-phylo` (NNI+SPR, HKY/GTR) | Roadmap B-8 |
| `valenx-adapter-mrbayes` | `mb` | Bayesian phylogeny MCMC | `valenx-phylo` (MH sampler) | Roadmap B-9 |
| `valenx-adapter-raxml-ng` | `raxml-ng` | ML phylogeny | `valenx-phylo` (NNI+SPR ML) | Roadmap B-10 |
| `valenx-adapter-beast2` | `beast` (Java) | Bayesian phylo + molecular clock | `valenx-phylo` (Bayesian MCMC) | Roadmap B-11 |
| `valenx-adapter-slim` | `slim` | Forward genetic simulation | `valenx-popgen` (Wright-Fisher) | Roadmap B-12 |
| `valenx-adapter-msprime` | Python `msprime` | Coalescent simulation | `valenx-popgen` (Kingman coalescent) | Roadmap B-13 |
| `valenx-adapter-tskit` | Python `tskit` | Tree-sequence analysis | `valenx-popgen` (tree-sequence) | Roadmap B-14 |
| `valenx-adapter-linearfold` | `linearfold` | Linear-time RNA folding | `valenx-rnastruct` (LinearFold) | Roadmap B-15 |
| `valenx-adapter-mfold` | `mfold` | RNA MFE folding | `valenx-rnastruct` (Zuker) | Roadmap B-16 |
| `valenx-adapter-rnastructure` | `Fold` (RNAstructure) | RNA folding (similar Turner model) | `valenx-rnastruct` (Zuker) | Roadmap B-17 |
| `valenx-adapter-eternafold` | `eternafold` | RNA folding (modified Turner) | `valenx-rnastruct` (Zuker, approx.) | Roadmap B-18 |
| `valenx-adapter-nupack` | Python `nupack` | RNA/DNA thermodynamics | `valenx-rnastruct` (ensemble/partition) | Roadmap B-19 |
| `valenx-adapter-rdkit` | Python `rdkit` | Cheminformatics | `valenx-cheminf` | Roadmap B-20 |
| `valenx-adapter-openbabel` | `obabel` | Molecule format conversion | `valenx-cheminf` (MOL/SMILES/InChI) | Roadmap B-21 |
| `valenx-adapter-biopython` | Python `biopython` | General bio Python API | `valenx-bioseq` + `valenx-align` etc. | Roadmap B-22 |
| `valenx-adapter-gatk` | Java GATK | Variant calling, haplotyping | `valenx-genomics` (variant caller) | Roadmap B-23 |
| `valenx-adapter-salmon` | `salmon` | RNA-seq quantification (pseudoalign) | `valenx-genomics` (quantification) | Roadmap B-24 |
| `valenx-adapter-kallisto` | `kallisto` | RNA-seq quantification (pseudoalign) | `valenx-genomics` (quantification) | Roadmap B-25 |
| `valenx-adapter-hisat2` | `hisat2` | Splice-aware RNA-seq aligner | `valenx-align` (splice-aware mapper) | Roadmap B-26 |
| `valenx-adapter-star` | `STAR` | Very large spliced aligner | `valenx-align` (large-scale read mapper) | Roadmap B-27 |
| `valenx-adapter-art` | `art_illumina` | Illumina read simulator | `valenx-genomics` (read simulation) | Roadmap B-28 |
| `valenx-adapter-wgsim` | `wgsim` | WGS read simulator | `valenx-genomics` (read simulation) | Roadmap B-29 |
| `valenx-adapter-badread` | `badread` | Long-read error simulator | `valenx-genomics` (long-read sim) | Roadmap B-30 |
| `valenx-adapter-tcoffee` | `t_coffee` | MSA + consistency | `valenx-align` (progressive+profile MSA) | Roadmap B-31 |
| `valenx-adapter-foldseek` | `foldseek` | Structure-based sequence search | `valenx-biostruct` + `valenx-align` | Roadmap B-32 |
| `valenx-adapter-dnachisel` | Python `dnachisel` | DNA sequence optimisation | `valenx-bioseq` (codon opt, constraints) | Roadmap B-33 |
| `valenx-adapter-pydna` | Python `pydna` | In-silico plasmid cloning | `valenx-bioseq` (cloning module) | Roadmap B-34 |
| `valenx-adapter-icodon` | Python `icodon` | Codon optimisation | `valenx-bioseq` (codon opt + CAI) | Roadmap B-35 |
| `valenx-adapter-lineardesign` | `LinearDesign` | mRNA codon + structure optimisation | `valenx-rnadesign` + `valenx-bioseq` | Roadmap B-36 |
| `valenx-adapter-forecast` | `FORECAST` | Codon optimisation | `valenx-bioseq` (codon opt) | Roadmap B-37 |
| `valenx-adapter-crispor` | Python `crispor` | CRISPR guide RNA design | `valenx-genediting` (guide design) | Roadmap B-38 |
| `valenx-adapter-crispritz` | `crispritz` | CRISPR off-target search | `valenx-genediting` (off-target) | Roadmap B-39 |
| `valenx-adapter-cas-offinder` | `cas-offinder` | CRISPR off-target search | `valenx-genediting` (off-target) | Roadmap B-40 |
| `valenx-adapter-primedesign` | Python `primedesign` | Prime editing pegRNA design | `valenx-genediting` (prime editing) | Roadmap B-41 |
| `valenx-adapter-pegfinder` | Python `pegfinder` | Prime editing strategy | `valenx-genediting` (prime editing) | Roadmap B-42 |
| `valenx-adapter-chopchop` | Python `chopchop` | CRISPR guide + efficiency scoring | `valenx-genediting` (guide design) | Roadmap B-43 |
| `valenx-adapter-dssr` | `dssr` / `x3dna-dssr` | RNA 3D structure annotation | `valenx-biostruct` (RNA geometry) | Roadmap B-44 |
| `valenx-adapter-x3dna` | `x3dna` | DNA helical parameters | `valenx-biostruct` (nucleic geometry) | Roadmap B-45 |
| `valenx-adapter-curves` | `curves+` | DNA helix parameters | `valenx-biostruct` (nucleic geometry) | Roadmap B-46 |
| `valenx-adapter-be-designer` | Python `be-designer` | Base editor design | `valenx-genediting` (base editing) | Roadmap B-47 |
| `valenx-adapter-prody` | Python `prody` | Normal-mode / elastic-network MD | `valenx-biostruct` (geometry + NMA) | Roadmap B-48 |

### 2B — Category C: NOT realistically nativizable without large pretrained weights

These tools require either multi-GB neural-network weights, closed-source
large-scale physics engines, or extensive desktop GUI frameworks.
Reimplementing them "in Rust in-house" would mean shipping those weights —
that is NOT a zero-download rewrite.

**Honest options for each:**
- (i) Keep as a cleanly gated OPTIONAL external integration (current status). Mark `probe()` failure gracefully in the UI ("install X for this capability").
- (ii) A native Rust ONNX/Candle inference runtime loading user-supplied weights. Feasible if/when Candle or tract reaches the maturity needed and the model supplier provides ONNX weights.

| Adapter | External tool | Reason NOT nativizable | Realistic path |
|---------|---------------|------------------------|----------------|
| `valenx-adapter-alphafold2` | AlphaFold 2 Python + JAX | ~10 GB weights; CUDA JAX required; 2-recycle inference is days of engineering | Option (ii): Rust ONNX runtime + user-downloaded weights |
| `valenx-adapter-alphafold3` | AlphaFold 3 Python + JAX | Same as AF2 + new diffusion head | Option (ii) |
| `valenx-adapter-colabfold` | ColabFold Python | AlphaFold2-class weights + MMseqs2 search | Option (ii) |
| `valenx-adapter-esmfold` | ESMFold Python + PyTorch | 690 M-param ESM-2 weights (~2.6 GB) | Option (ii): ONNX export feasible |
| `valenx-adapter-rosettafold` | RoseTTAFold Python | 3-track network, ~700 MB weights | Option (ii) |
| `valenx-adapter-rfdiffusion` | RFdiffusion Python | Diffusion model, ~700 MB weights | Option (ii) |
| `valenx-adapter-proteinmpnn` | ProteinMPNN Python | GNN, ~40 MB weights (small!) | Option (ii): ONNX export is feasible short-term |
| `valenx-adapter-esm-if` | ESM-IF1 Python | GVP-GNN, ~140 MB weights | Option (ii): ONNX export feasible |
| `valenx-adapter-chroma` | Chroma Python | Diffusion + GNN, multi-GB | Option (ii) |
| `valenx-adapter-deepvariant` | DeepVariant Python + TF | Pileup-based CNN; large model | Option (i): keep as optional; `valenx-genomics` native caller covers most uses |
| `valenx-adapter-alphamissense` | AlphaMissense Python | AF2 backbone + missense head | Option (ii) |
| `valenx-adapter-esm3` | ESM3 Python + PyTorch | ~98 B-param model (full), smaller open versions | Option (ii): ESM3-open ONNX feasible |
| `valenx-adapter-esmc` | ESMC Python | Protein language model, multi-GB | Option (ii) |
| `valenx-adapter-omegafold` | OmegaFold Python | 73 M-param model | Option (ii) |
| `valenx-adapter-openfold` | OpenFold Python | AlphaFold2 reimplementation | Option (ii) |
| `valenx-adapter-rosetta` | Rosetta C++ (massive) | Multi-year codebase; Apache-2.0 but requires own build chain | Option (i): keep as optional |
| `valenx-adapter-pyrosetta` | PyRosetta Python | Python Rosetta API | Option (i) |
| `valenx-adapter-rfantibody` | RFAntibody Python | RFdiffusion-based | Option (ii) |
| `valenx-adapter-be-hive` | Python BE-Hive | ML-based base-edit outcome prediction | Option (i) |
| `valenx-adapter-indelphi` | Python inDelphi | ML NHEJ/MMEJ repair prediction | Option (i) |
| `valenx-adapter-deepchem` | Python DeepChem | ML-over-molecules, many model families | Option (i): pure-algo parts in `valenx-cheminf` |
| `valenx-adapter-scanpy` | Python Scanpy | Single-cell RNA-seq analysis ecosystem | Option (i) |
| `valenx-adapter-scvi` | Python scVI | VAE-based single-cell model | Option (ii) |
| `valenx-adapter-seurat` | R Seurat | R-ecosystem single-cell analysis | Option (i) |
| `valenx-adapter-anndata` | Python AnnData | Single-cell data container (no algorithm) | Option (i) |
| `valenx-adapter-ilastik` | ilastik Python | Interactive ML image analysis | Option (i) |
| `valenx-adapter-cellprofiler` | CellProfiler Python | Cell image segmentation + ML | Option (i) |
| `valenx-adapter-fiji` | ImageJ/Fiji JVM | Java image-analysis platform + plugins | Option (i) |
| `valenx-adapter-amber-sander` | AMBER SANDER Fortran | Molecular dynamics; complex physics + GPU | Option (i) |
| `valenx-adapter-namd` | NAMD C++ | MD with CUDA; NAMD license | Option (i) |
| `valenx-adapter-openmm` | OpenMM Python + CUDA | GPU MD; Python API + CUDA kernels | Option (i) |
| `valenx-adapter-hoomd` | HOOMD-blue Python + CUDA | GPU MD/DPD | Option (i) |
| `valenx-adapter-plumed` | PLUMED C++ | Collective variables + free energy | Option (i) |
| `valenx-adapter-cpptraj` | cpptraj C++ | MD trajectory analysis (~40k lines) | Option (i): basic ops possible in `valenx-md` |
| `valenx-adapter-mdanalysis` | Python MDAnalysis | MD trajectory Python library | Option (i) |
| `valenx-adapter-mdtraj` | Python MDTraj | MD trajectory Python library | Option (i) |
| `valenx-adapter-oxdna` | oxDNA C++ | Coarse-grained DNA/RNA MD | Option (i) |
| `valenx-adapter-psi4` | PSI4 C++/Python | Quantum chemistry; multi-Mline codebase | Option (i) |
| `valenx-adapter-nwchem` | NWChem Fortran | Quantum chemistry | Option (i) |
| `valenx-adapter-xtb` | xTB Fortran | Semiempirical tight-binding QM | Option (i): ~30 k lines; feasible Rust port long-term |
| `valenx-adapter-copasi` | COPASI C++ | ODE/stochastic systems biology | Option (i) |
| `valenx-adapter-bionetgen` | BioNetGen Perl/C++ | Rule-based network simulator | Option (i) |
| `valenx-adapter-physicell` | PhysiCell C++ | Agent-based cell simulation | Option (i) |
| `valenx-adapter-smoldyn` | Smoldyn C++ | Spatial stochastic simulation | Option (i) |
| `valenx-adapter-mcell` | MCell C | Monte Carlo cell simulation | Option (i) |
| `valenx-adapter-relion` | RELION C++/CUDA | Cryo-EM refinement; GPU-required | Option (i) |
| `valenx-adapter-ctffind` | CTFFIND Fortran | Cryo-EM CTF estimation | Option (i) |
| `valenx-adapter-eman2` | EMAN2 Python/C++ | Cryo-EM image processing | Option (i) |
| `valenx-adapter-simrna` | SimRNA C++ | RNA 3D structure prediction (physics) | Option (i) |
| `valenx-adapter-nextflow` | Nextflow JVM/DSL | Workflow engine | Option (i): orchestration layer, not an algorithm |
| `valenx-adapter-snakemake` | Snakemake Python | Workflow engine | Option (i) |
| `valenx-adapter-cromwell` | Cromwell JVM/WDL | Workflow engine | Option (i) |
| `valenx-adapter-cwltool` | cwltool Python | CWL workflow executor | Option (i) |
| `valenx-adapter-planemo` | Planemo Python | Galaxy workflow engine | Option (i) |
| `valenx-adapter-j5` | j5 Java | DNA assembly design web service | Option (i) |
| `valenx-adapter-cello` | Cello Java | Genetic circuit design | Option (i) |
| `valenx-adapter-pysbol` | Python pySBOL | SBOL synthetic-biology format | Option (i): pure I/O, feasible Rust port |
| `valenx-adapter-igv` | IGV Java | Genome browser GUI | Option (i): headless operations possible |
| `valenx-adapter-jalview` | Jalview Java | Alignment viewer GUI | Option (i) |
| `valenx-adapter-pymol` | PyMOL Python/C++ | Molecular visualization GUI | Option (i) |
| `valenx-adapter-chimerax` | ChimeraX Python/C++ | Molecular visualization GUI | Option (i) |
| `valenx-adapter-vmd` | VMD Tcl/C++ | Molecular dynamics visualization GUI | Option (i) |
| `valenx-adapter-molstar` | Mol* TypeScript | Web molecular viewer | Option (i) |
| `valenx-adapter-ngl` | NGL TypeScript | Web molecular viewer | Option (i) |
| `valenx-adapter-avogadro` | Avogadro C++/Qt | Desktop molecule editor GUI | Option (i) |
| `valenx-adapter-pksim` | PKSim C# | Pharmacokinetics simulation | Option (i) |

---

## 3. What Was Implemented in This Branch

Seven adapters received a native Rust code path that makes the adapter **always
available** — `probe()` succeeds even when the external binary is not installed.
The external binary becomes an optional accelerator/fallback.

### Architecture

Each nativized adapter gains a `native.rs` module plus these changes to `lib.rs`:

1. **`probe()`** — tries the external binary first; if missing, returns
   `ProbeReport { ok: true, found_version: Some("native-rust"), ... }` with a
   note explaining native mode. The adapter is always `AdapterStatus::Ready`.

2. **`prepare()`** — always writes a `native_params.toml` file to the workdir
   with the resolved job parameters. If external binary is available and
   `[bio.<tool>].prefer_external = true`, sets `native_command` to the real
   binary invocation. Otherwise, sets `native_command[0]` to the sentinel
   `"valenx:native:<id>"`.

3. **`run()`** — detects the `"valenx:native:"` sentinel in `native_command[0]`;
   if present, reads `native_params.toml` from the workdir and calls the
   corresponding native Rust API. The algorithm writes output files to the
   workdir in the exact same filenames and formats that `collect()` expects —
   so `collect()` is unchanged.

4. **`collect()`** — unchanged; reads the same output files.

### Adapter-specific notes

| Adapter | Native impl | Output format | Reference validation |
|---------|------------|---------------|---------------------|
| `valenx-adapter-blast` | `valenx_align::search::{KmerIndex, SeedSearch, KarlinAltschul}` | BLAST tabular format 6 (12 columns) or custom header; outfmt 0 pairwise text | E-values validated against published Karlin-Altschul λ=0.318, K=0.134 for BLOSUM62 |
| `valenx-adapter-hmmer` | `valenx_align::hmm::ProfileHmm` | HMMER tblout + report format | Plan7 Viterbi scoring; bit-scores computed from log-odds |
| `valenx-adapter-muscle` | `valenx_align::{align_msa, msa::refine}` | aligned FASTA | Sum-of-pairs score improvement, validated against small reference alignment |
| `valenx-adapter-mafft` | `valenx_align::{align_msa, msa::refine}` | aligned FASTA | Same underlying native MSA; MAFFT-specific extra_args ignored in native mode |
| `valenx-adapter-clustalo` | `valenx_align::align_msa` | aligned FASTA | Progressive UPGMA guide-tree alignment |
| `valenx-adapter-fasttree` | `valenx_phylo::{distance_matrix, neighbor_joining, optimize_topology_ml_spr}` | Newick | NJ tree from JC distances; optional ML refinement under JC/GTR; tested on 4-taxon reference case |
| `valenx-adapter-viennarna` | `valenx_rnastruct::fold::zuker::{mfe, mfe_d2}` | ViennaRNA stdout format (sequence + dot-bracket + energy) | Turner-2004 energies validated; GCGGAUUUA test case matches published value |

### Limitation honest note

The native BLAST path operates on **FASTA-format databases only**. Binary BLAST
databases (pre-built with `makeblastdb`) require the external BLAST+ binary.
Users with FASTA databases (most in-house use cases) benefit from zero-install;
users needing the NCBI `nt`/`nr` binary DBs must have BLAST+ installed (the
adapter falls back gracefully to the subprocess path when available).

The native tree-building (FastTree equivalent) uses NJ + NNI/SPR ML, which is
comparable in quality to FastTree but is not the exact same algorithm. Results
are scientifically valid but will differ numerically from FastTree output.

The native RNA folding uses the identical Turner-2004 parameters as ViennaRNA's
`RNAfold -d2` (with coaxial stacking). Output format is compatible. The
license restriction on ViennaRNA does NOT apply to the native implementation,
which uses only the publicly published parameter set.

---

## 4. Prioritized Roadmap (B-items not yet finished)

### Tier 1 — High impact, simple API bridging (1–3 days each)

| # | Adapter | Key work |
|---|---------|----------|
| B-1 | `bowtie2` | Wire `valenx_align::search::FmIndex` + seed-extend as native short-read mapper; write SAM output |
| B-2 | `bwa` | Same as bowtie2 — FM-index + banded NW; write SAM |
| B-3 | `minimap2` | Wire `valenx_align::search::minimizer_sketch` + `chain_anchors`; PAF output |
| B-4 | `samtools` | Wire `valenx_genomics` SAM parsing + sort/index/view ops |
| B-5 | `bcftools` | Wire `valenx_genomics` VCF ops |
| B-8 | `iqtree` | Wire `valenx_phylo::likelihood::optimize_topology_ml_spr` with model selection |
| B-9 | `mrbayes` | Wire `valenx_phylo::bayes` MCMC |
| B-12 | `slim` | Wire `valenx_popgen::forward` WF simulator |
| B-13 | `msprime` | Wire `valenx_popgen::coalescent` |
| B-15 | `linearfold` | Wire `valenx_rnastruct::fold::linear::fold_linear` |
| B-20 | `rdkit` | Wire `valenx_cheminf` for fingerprints/similarity/descriptors |
| B-21 | `openbabel` | Wire `valenx_cheminf` MOL/SMILES/InChI conversion |

### Tier 2 — More complex output format compatibility (3–7 days each)

| # | Adapter | Key challenge |
|---|---------|---------------|
| B-6 | `diamond` | Protein DB indexing; same output format as BLAST-tabular |
| B-7 | `mmseqs2` | k-mer cascading index; cascaded alignment |
| B-10 | `raxml-ng` | Same ML as IQ-TREE tier; RAxML-NG partition models |
| B-11 | `beast2` | Bayesian MCMC with clock models; NEXUS + BEAST XML I/O |
| B-22 | `biopython` | Route each Biopython op to the relevant native crate |
| B-23 | `gatk` | Variant-calling pipeline; GATK-specific GVCF format |
| B-24 | `salmon` | Pseudoalignment EM; transcript quantification output |
| B-25 | `kallisto` | k-mer pseudoalignment; quant output format |

### Tier 3 — Category C strategy for ML tools

For the ML tools (AlphaFold, ESMFold, ProteinMPNN, etc.) the recommended
path is **Candle/tract Rust inference + user-supplied ONNX weights**:

1. Define a `weights_path` field in the adapter's case input.
2. Use [Candle](https://github.com/huggingface/candle) or
   [tract](https://github.com/sonos/tract) for zero-dependency ONNX inference.
3. Keep the external-binary subprocess path as a fallback for users who already
   have the tools installed.
4. Ship the small models (ProteinMPNN ~40 MB, ESM-IF ~140 MB) as optional
   one-time downloads with a clear prompt; the large models (AF2/AF3 > 10 GB)
   stay external-only.

Priority order for Candle-backed native paths:
1. ProteinMPNN (smallest weights, most common use case)
2. ESM-IF (moderate size, inverse folding)
3. ESMFold (large but ONNX-exportable)
4. RFdiffusion (popular de novo design)
5. AlphaFold 2/3 (always external; too large for in-house weights)

---

## 5. Deletions

None. No files were deleted in this branch. All existing adapter code is
preserved; the native paths are additive modifications.

---

## 6. How to Use Native Mode

After this branch, users with no external bio tools installed will see, e.g.
for MUSCLE:

```
[bio.muscle]
input = "seqs.fa"
# no prefer_external = true, so native Rust MSA is used automatically
```

Users who have MUSCLE 5 installed and want to use it:

```
[bio.muscle]
input          = "seqs.fa"
prefer_external = true   # falls back to native if muscle binary not found
```

The adapter probe output in the UI will show:
- `native-rust (valenx-align v<version>)` when running in native mode
- `muscle 5.1 (native-rust fallback available)` when the binary is found
