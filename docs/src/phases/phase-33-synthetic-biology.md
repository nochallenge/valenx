# Phase 33 — Synthetic biology

**Status:** 🟢 Live — pySBOL + j5 + Cello open the
**first synthetic biology / genetic-circuit design domain** in
Valenx alongside the Phase 5.5 / 17 / 17.5 / 18 / 25 / 27 / 27.5 /
27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 34 / 35 / 36 / 38 / 39
MD-analysis-expansion + biology + structure-prediction + protein-
design + RNA-structure + population-genetics + phylogenetics +
Bayesian-phylogenetics + read-simulator + systems-biology +
docking + CRISPR-design + cryo-EM + Rosetta-family + DNA-geometry
beachheads and the Phase 24 cheminformatics expansion.

## Goal

Open the synthetic biology / genetic-circuit design domain in
Valenx with three established open-source tools that span the
synthetic-biology tradeoff space — a canonical SBOL-standard
Python composition library (pySBOL, the reference implementation
of the Synthetic Biology Open Language for capturing genetic
designs as round-trippable RDF/XML or JSON-LD), DNA assembly
automation that plans the optimal Gibson / Golden-Gate / SLIC /
SLIM strategy from a target circuit + parts library (j5, JBEI's
canonical assembly automator), and genetic-circuit DNA compilation
from a Verilog netlist describing the desired logic function (Cello
v2, the canonical CIDAR genetic-circuit DNA compiler that runs
simulated-annealing optimization over the gate-assignment problem).
pySBOL follows the established Phase 17 Biopython Python-script
subprocess shape: the user supplies a Python script that imports
the upstream package and reads `valenx_params.json` for the parsed
knobs. j5 + Cello are JAR-distributed (no `j5` / `cello` launcher
binary on PATH); the user supplies the absolute path to the JAR
via case input, and we probe `java` itself rather than the JAR —
different sites pin different j5 / Cello releases under different
paths. Phase 33 sits numerically between Phase 32 systems biology
and Phase 34 molecular docking and ships chronologically right
after Phase 39 DNA structural geometry — same chronological-vs-
numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39.

## Capability inventory

### Live adapters (3)

- **pySBOL** — the Python implementation (pySBOL3) of the
  [Synthetic Biology Open Language (SBOL)](https://sbolstandard.org/)
  (Apache-2.0). SBOL captures components, sequences, interactions,
  constraints, and the full provenance of a synthetic design as
  RDF/XML or JSON-LD that round-trips with every SBOL-conformant
  tool: j5, Cello, SynBioHub, iBioSim, ... Python-script
  subprocess shape (sister to Phase 17 Biopython): the user supplies
  a Python script referenced from `[bio.pysbol].script` in
  `case.toml` that imports `sbol3` and reads `valenx_params.json`
  for the parsed knobs. Schema knobs: `script` (path to user-
  supplied Python script; required), `python` (interpreter name;
  default `"python3"`), `input_sbol` (optional starting SBOL XML
  document; `None` when the script generates the design from
  scratch), `output_basename` (filename stem the user's script uses
  for outputs — surfaced here so collect() can label artifacts
  uniformly; required, non-empty). `prepare()` stages the script +
  optional input SBOL into the workdir under their original
  filenames so the script can resolve them via relative paths,
  then writes a flat `valenx_params.json` containing `input_sbol`
  (staged filename or literal `null`) and `output_basename`.
  `collect()` walks the workdir for `<output_basename>*.xml`
  (`Tabular`, "pySBOL document") and `<output_basename>*.json`
  (`Log`, "pySBOL composition log"). Probe via Python on PATH
  with an `import sbol3` check — when the import fails the probe
  still returns `ok = true` with a warning so users with pySBOL
  installed under a non-standard module name aren't blocked. Version
  range `3.0.0..4.0.0` (pySBOL3 is the modern Python rewrite — the
  older 2.x line is deprecated; 3.0 is the floor; upper bound 4.0
  reserves room for an eventual major bump). The init alias `sbol`
  resolves to the same template as the canonical `pysbol` name.
  `bio.pysbol.compose` ribbon capability.
- **j5** — JBEI's canonical DNA-assembly automation tool
  (BSD-3-Clause). j5 consumes a target circuit design (CSV row per
  cassette) plus a parts library (CSV row per part / oligo), then
  plans the optimal Gibson / Golden-Gate / SLIC / SLIM assembly
  strategy and writes the per-step protocol + GenBank construct
  files. Single-binary subprocess shape (sister to Phase 18 BWA)
  but **JAR-distributed**: no `j5` launcher binary on PATH; the
  user supplies the absolute path to `j5.jar` via `[bio.j5].jar`
  in `case.toml`. The CLI is `java -jar <jar> -d <design_csv> -p
  <parts_csv> -o <output_basename> [extras...]`. Schema knobs:
  `jar` (absolute path to `j5.jar`; required), `design_csv` (j5
  design CSV with parts, oligos, target; required), `parts_csv`
  (parts list CSV; required), `output_basename` (filename stem the
  user expects j5 to produce; required, non-empty), `extra_args`.
  `prepare()` resolves all three input paths (jar + design CSV +
  parts CSV) against the case directory when relative, validates
  each file exists on disk (returns `InvalidCase` with a helpful
  message when missing), and composes the `java -jar` invocation.
  `collect()` walks the workdir for `<output_basename>*.csv`
  (`Tabular`, "j5 assembly plan") and `<output_basename>*.gb`
  (`Native`, "j5 GenBank output"). Probe via
  `find_on_path(&["java"])` — j5's version comes from the jar
  itself, not from `java`, so we surface no version here; the user
  pins the j5 release implicitly by the jar they point at. Version
  range `1.0.0..2.0.0` (j5 has been on a 1.x line for over a
  decade; upper bound 2.0 reserves room for an eventual major bump).
  `bio.j5.assemble` ribbon capability.
- **Cello** — CIDAR's canonical genetic-circuit DNA compiler
  ([Cello v2](https://github.com/CIDARLAB/Cello-v2), BSD-3-Clause).
  Cello consumes a Verilog netlist describing the desired logic
  function plus a triplet of JSON constraint files (a user
  constraint file pinning the chassis / library, an input sensor
  file pinning the input promoters, an output device file pinning
  the reporter), and emits a fully assembled DNA construct that
  implements the logic in a living cell. The compiler runs a
  simulated-annealing optimization over the gate-assignment problem
  and outputs a Graphviz `.dot` netlist, a circuit diagram PNG, and
  a human-readable report. Single-binary subprocess shape (sister
  to j5) but **JAR-distributed**: no `cello` launcher binary on
  PATH; the user supplies the absolute path to the jar via
  `[bio.cello].jar` in `case.toml`. The CLI is `java -jar <jar>
  -inputNetlist <verilog> -targetDataFile <user_constraints>
  -inputSensorFile <input_sensors> -outputDeviceFile
  <output_devices> -outputDir <output_basename> [extras...]`.
  Schema knobs: `jar` (absolute path to the Cello jar; required),
  `verilog` (`.v` Verilog circuit description; required),
  `user_constraints` (`.UCF` user constraints file pinning the
  chassis / library; required), `input_sensors` (`.input.json`
  pinning the input promoters; required), `output_devices`
  (`.output.json` pinning the reporter; required),
  `output_basename` (filename stem Cello uses for the output
  directory; required, non-empty), `extra_args`. `prepare()`
  resolves all five input paths against the case directory when
  relative, validates each file exists on disk, and composes the
  `java -jar` invocation. `collect()` walks the workdir for
  `<output_basename>*.txt` (`Log`, "Cello report"),
  `<output_basename>*.png` (`Native`, "Cello circuit diagram"), and
  `<output_basename>*.dot` (`Native`, "Cello Graphviz netlist").
  Probe via `find_on_path(&["java"])` — same JAR-versioning shape
  as j5: Cello's version comes from the jar itself. Version range
  `2.0.0..3.0.0` (Cello v2 is the modern Java rewrite (2020+); the
  v1 line was Python and is deprecated; upper bound 3.0 reserves
  room for an eventual major bump). `bio.cello.compile` ribbon
  capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (pySBOL Python composition scripts + optional
starting SBOL XML, j5 design + parts CSVs + the j5 jar, Cello
Verilog netlist + UCF + input-sensor + output-device JSONs + the
Cello jar) and emit user-readable artifacts (pySBOL XML / JSON
documents, j5 assembly-plan CSVs + GenBank `.gb` constructs, Cello
Graphviz `.dot` netlists + circuit-diagram PNGs + human-readable
text reports) that the unchanged `Results.artifacts` collection
model surfaces directly. A first-class synthetic-biology canonical
type — a typed SBOL-document representation spanning pySBOL output
+ j5 GenBank + Cello netlists, with parsed component / sequence /
interaction graphs — defers to a future phase along with circuit-
diagram visualizers and per-construct interactive overlays.

### Headless CLIs

**No new CLIs.** pySBOL's XML / JSON SBOL documents, j5's CSV
assembly plans + GenBank construct files, and Cello's Graphviz
`.dot` netlists + PNG circuit diagrams + text reports are all
standard formats inspectable in any editor or through the user's
downstream Python pipeline (`sbol3.Document`, `Bio.SeqIO`,
`graphviz`). A canonical synthetic-biology CLI — SBOL-document
inspection, j5 plan diffing, Cello netlist comparison — defers to
a future phase along with the canonical type.

## Domain milestone

Phase 33 is the **first synthetic biology / genetic-circuit design
domain** to land in Valenx. The biology adapter family started
with Phase 17 (foundation — sequence / structure / trajectory
canonical types + classical MD + cheminformatics scripts) and
expanded through Phase 5.5 / 17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 /
20 / 22 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5 /
31 / 32 / 34 / 35 / 36 / 38 / 39 to cover MD-trajectory analysis
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
compilation) was absent. Phase 33 closes that gap with three
established open-source tools spanning the synthetic-biology
tradeoff space — pySBOL at the SBOL-standard Python composition
end, j5 for the canonical DNA-assembly automation surface, and
Cello as the genetic-circuit DNA compiler that turns Verilog into
working living-cell logic.

## What landed early

The implementation landed across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-pysbol` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-input-sbol / rejects-empty-output-basename, plus
      the Python-script subprocess shape that stages script +
      optional input SBOL and writes `valenx_params.json` with
      `input_sbol` (staged filename or literal `null`) and
      `output_basename`
- [x] `valenx-adapter-j5` adapter ships with case-input parser +
      4 lib tests + 3 case-input tests covering parses-minimal /
      rejects-empty-jar / rejects-empty-output-basename, plus the
      JAR-distributed single-binary subprocess shape that probes
      `java` on PATH and composes `java -jar <jar> -d <design_csv>
      -p <parts_csv> -o <output_basename> [extras...]`
- [x] `valenx-adapter-cello` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      rejects-empty-jar / rejects-empty-verilog / rejects-empty-
      output-basename, plus the JAR-distributed single-binary
      subprocess shape that probes `java` on PATH and composes
      `java -jar <jar> -inputNetlist <verilog> -targetDataFile
      <user_constraints> -inputSensorFile <input_sensors>
      -outputDeviceFile <output_devices> -outputDir
      <output_basename> [extras...]`
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 94 to **100** (alongside the
      Phase 5.5 MD-analysis-expansion trio), opening the first
      synthetic biology / genetic-circuit design domain to ship in
      Valenx and crossing the **100-adapter milestone**
- [x] 3 synthetic-biology templates in `valenx-init` (`pysbol`
      with alias `sbol`, `j5`, `cello`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps **96
      templates** clean alongside the Phase 5.5 MD-analysis-
      expansion trio)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 33.5** — sister-adapter expansion of Phase 33:
      libSBOL (the C++ implementation of SBOL — sister to pySBOL
      with a different language surface; defer to sister-adapter
      expansion phase), iBioSim (Bridges lab's SBOL-conformant
      genetic-circuit modeling + simulation environment; defer),
      SBOLDesigner (Anderson lab's drag-and-drop GUI for SBOL
      composition; the GUI shape doesn't fit the headless adapter
      pattern; defer), SynBioHub (Watson lab's online repository
      for sharing SBOL designs; the web-service shape doesn't fit
      the local-binary adapter pattern; defer), Tellurium (Sauro
      lab's Python environment for systems / synthetic biology;
      sister to PySCeS but adjacent to Phase 32 systems biology;
      defer), GeneticCircuitGenerator (CIDAR's combinatorial
      circuit-library enumeration; defer to 33.5). Out of scope
      for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New synthetic-biology adapter (template + tests)      | 1 day per       |
| Compose-SBOL → assemble-DNA → compile-circuit loop across 3 tools | < tool baseline |

## Leads into

Phase 33 opens the synthetic biology / genetic-circuit design
domain that the user's bio / chemistry spec called out alongside
the Phase 17 / 17.5 biology + structure-prediction stack, the Phase
27 / 27.5 / 27.6 protein-design beachhead, and the Phase 32
systems-biology surface. Combined with the existing simulate-MD →
analyze-trajectory → reweight-free-energy → predict-structure →
fold-RNA → analyze-DNA-geometry → infer-tree-ML → infer-tree-
Bayesian → simulate-popgen → analyze-trees → simulate-pathway →
reconstruct-3D → design-protein → validate loop, the **compose-
SBOL → assemble-DNA → compile-circuit → simulate-MD → analyze-
trajectory → predict-structure → fold-RNA → analyze-DNA-geometry
→ infer-tree-ML → infer-tree-Bayesian → simulate-popgen →
analyze-trees → simulate-pathway → reconstruct-3D → design-protein
→ validate** loop now spans three synthetic-biology tools (pySBOL,
j5, Cello) feeding into the existing Phase 5.5 MD-analysis trio
(PLUMED, ProDy, cpptraj), the Phase 17 / 17.5 prediction stack
(ESMFold, OpenFold, AlphaFold 2/3, ColabFold), the Phase 28
RNA-structure tools (ViennaRNA, RNAstructure, NUPACK), the Phase
29 population-genetics trio (SLiM, msprime, tskit), the Phase 30
phylogenetic-tree builders (IQ-TREE, RAxML-NG, FastTree), the
Phase 30.5 Bayesian-phylogenetics pair (BEAST 2, MrBayes), the
Phase 32 systems-biology surface (COPASI, BioNetGen, PhysiCell),
the Phase 34 docking pair (AutoDock Vina, AutoDock 4), the Phase
35 CRISPR-design tools (CHOPCHOP, CRISPOR, Cas-OFFinder), the
Phase 36 cryo-EM reconstruction tools (RELION, EMAN2, CTFFIND),
the Phase 38 Rosetta-family adapters (Rosetta, PyRosetta), and the
Phase 39 DNA-structural-geometry tools (X3DNA, Curves+, DSSR) —
all in one Valenx shell with no glue code beyond the existing
case-toml / prepare / run / collect path.

The natural follow-up is **Phase 33.5** — the deferred synthetic-
biology work called out above (libSBOL as the C++ SBOL sister to
pySBOL, iBioSim for SBOL-conformant modeling + simulation,
SBOLDesigner / SynBioHub for the design-sharing surface,
Tellurium for the Python systems / synthetic biology environment,
GeneticCircuitGenerator for combinatorial circuit-library
enumeration), slotting in alongside the existing pySBOL + j5 +
Cello adapters with the same Python-script subprocess shape (pySBOL
sister tools) or the JAR-distributed single-binary subprocess shape
(j5 / Cello sister tools).
