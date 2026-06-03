# Summary

[Introduction](./introduction.md)

-----------

- [Reference](./reference/README.md)
  - [CLI reference](./reference/cli.md)

- [Per-phase acceptance docs](./phases/README.md)
  - [Phase 1 — Foundation + first CFD thread](./phases/phase-01-foundation.md)
  - [Phase 2 — CAD + meshing integration](./phases/phase-02-cad-meshing.md)
  - [Phase 3 — Finite-element analysis](./phases/phase-03-fea.md)
  - [Phase 4 — Reaction kinetics + combustion](./phases/phase-04-chemistry.md)
  - [Phase 5 — Molecular dynamics](./phases/phase-05-md.md)
  - [Phase 5.5 — MD analysis expansion](./phases/phase-5-5-md-analysis.md)
  - [Phase 5.6 — Bio MD engines](./phases/phase-5-6-md-engines.md)
  - [Phase 5.7 — MDTraj](./phases/phase-5-7-mdtraj.md)
  - [Phase 6 — Electromagnetics](./phases/phase-06-em.md)
  - [Phase 7 — Battery modelling](./phases/phase-07-battery.md)
  - [Phase 8 — Multibody + robotics](./phases/phase-08-dynamics.md)
  - [Phase 9 — Multi-physics coupling](./phases/phase-09-coupling.md)
  - [Phase 10 — UX polish + first public release](./phases/phase-10-ux-polish.md)
  - [Phase 11 — HPC / cluster execution](./phases/phase-11-hpc.md)
  - [Phase 12 — Optimisation + adjoint workflows](./phases/phase-12-optimisation.md)
  - [Phase 13 — ML-assisted surrogates](./phases/phase-13-ml.md)
  - [Phase 14 — Plugin marketplace](./phases/phase-14-plugins.md)
  - [Phase 15 — Enterprise deployment](./phases/phase-15-enterprise.md)
  - [Phase 16 — Stewardship + long-term governance](./phases/phase-16-stewardship.md)
  - [Phase 17 — Biology + biotech foundation](./phases/phase-17-biology.md)
  - [Phase 17.5 — Structure prediction expansion](./phases/phase-17-5-structure-prediction.md)
  - [Phase 17.7 — Structure tools expansion](./phases/phase-17-7-structure-tools.md)
  - [Phase 18 — Sequence alignment toolkit](./phases/phase-18-alignment.md)
  - [Phase 18.5 — Aligners expansion](./phases/phase-18-5-aligners.md)
  - [Phase 18.6 — RNA-seq alignment](./phases/phase-18-6-rna-seq.md)
  - [Phase 18.7 — Alignment toolkit expansion](./phases/phase-18-7-blast-alignment.md)
  - [Phase 19 — Variant calling toolkit](./phases/phase-19-variant-calling.md)
  - [Phase 19.5 — Single-cell genomics](./phases/phase-19-5-single-cell.md)
  - [Phase 19.6 — Single-cell expansion](./phases/phase-19-6-single-cell-expansion.md)
  - [Phase 20 — Transcript quantification](./phases/phase-20-transcript-quantification.md)
  - [Phase 22 — Workflow managers](./phases/phase-22-workflow-managers.md)
  - [Phase 22.5 — Workflow expansion](./phases/phase-22-5-workflow-expansion.md)
  - [Phase 23 — Molecular viewers](./phases/phase-23-viewers.md)
  - [Phase 24 — Cheminformatics expansion](./phases/phase-24-cheminformatics.md)
  - [Phase 25 — Quantum chemistry](./phases/phase-25-quantum-chemistry.md)
  - [Phase 27 — Protein design](./phases/phase-27-protein-design.md)
  - [Phase 27.5 — Protein design expansion](./phases/phase-27-5-protein-design-expansion.md)
  - [Phase 27.6 — EvolutionaryScale models](./phases/phase-27-6-evolutionaryscale.md)
  - [Phase 28 — RNA structure](./phases/phase-28-rna-structure.md)
  - [Phase 29 — Population genetics](./phases/phase-29-population-genetics.md)
  - [Phase 30 — Phylogenetics](./phases/phase-30-phylogenetics.md)
  - [Phase 30.5 — Bayesian phylogenetics](./phases/phase-30-5-bayesian-phylogenetics.md)
  - [Phase 31 — Sequencing read simulators](./phases/phase-31-read-simulators.md)
  - [Phase 32 — Systems biology](./phases/phase-32-systems-biology.md)
  - [Phase 32.5 — Spatial stochastic](./phases/phase-32-5-spatial-stochastic.md)
  - [Phase 33 — Synthetic biology](./phases/phase-33-synthetic-biology.md)
  - [Phase 34 — Molecular docking](./phases/phase-34-docking.md)
  - [Phase 35 — CRISPR design](./phases/phase-35-crispr-design.md)
  - [Phase 35.5 — Base + prime editing design](./phases/phase-35-5-editing-design.md)
  - [Phase 35.6 — Edit-outcome prediction](./phases/phase-35-6-edit-outcomes.md)
  - [Phase 36 — Cryo-EM](./phases/phase-36-cryo-em.md)
  - [Phase 38 — Rosetta family](./phases/phase-38-rosetta.md)
  - [Phase 39 — DNA structural geometry](./phases/phase-39-dna-geometry.md)
  - [Phase 40 — Microscopy](./phases/phase-40-microscopy.md)
  - [Phase 41 — Sequence editors](./phases/phase-41-sequence-editors.md)
  - [Phase 42 — Web visualization](./phases/phase-42-web-visualization.md)
  - [Phase 43 — mRNA design](./phases/phase-43-mrna-design.md)
  - [Phase 44.5 — RNA folding expansion](./phases/phase-44-5-rna-folding-expansion.md)
  - [Phase 45 — Pharmacokinetics + RNA tertiary](./phases/phase-45-pk-rna-tertiary.md)

-----------

[Contributing](./contributing.md)
[Changelog](./changelog.md)

<!--
  Future chapters (not yet written — add them here as they're authored
  to keep them listed in the sidebar). Keep the `create-missing = false`
  setting in book.toml so we don't end up with dozens of empty files in
  the repo.

  # Getting started
  - [Installation](./getting-started/installation.md)
  - [First-run setup](./getting-started/first-run.md)
  - [Your first simulation](./getting-started/first-simulation.md)

  # Concepts
  - [Projects](./concepts/projects.md)
  - [Physics and solvers](./concepts/physics.md)
  - [Adapters and the tool registry](./concepts/adapters.md)
  - [Units and numerical conventions](./concepts/units.md)

  # Physics
  - [CFD](./physics/cfd.md)
  - [FEA](./physics/fea.md)
  - [Electromagnetics](./physics/em.md)
  - [Chemistry](./physics/chem.md)
  - [Molecular dynamics](./physics/md.md)
  - [Battery modelling](./physics/battery.md)
  - [Multi-physics coupling](./physics/multi-physics.md)

  # Workflow
  - [Geometry and CAD](./workflow/geometry.md)
  - [Meshing](./workflow/meshing.md)
  - [Boundary conditions](./workflow/boundary-conditions.md)
  - [Running a solver](./workflow/running.md)
  - [Post-processing](./workflow/post-processing.md)
  - [Parametric sweeps](./workflow/sweeps.md)
  - [Reports and export](./workflow/reports.md)

  # Reference
  - [Keyboard shortcuts](./reference/shortcuts.md)
  - [Settings](./reference/settings.md)
  - [File formats](./reference/file-formats.md)
  - [CLI](./reference/cli.md)

  # Scripting
  - [Python scripting](./scripting/python.md)
  - [Lua scripting](./scripting/lua.md)
  - [The `valenx.scripting` module](./scripting/api.md)

  # Plugins
  - [Installing plugins](./plugins/installing.md)
  - [Authoring a plugin](./plugins/authoring.md)

  # Tutorials
  - [External aerodynamics — NACA 0012 airfoil](./tutorials/airfoil.md)
  - [Internal flow — pipe with heat transfer](./tutorials/pipe-heat.md)
  - [Linear static — cantilever beam](./tutorials/cantilever.md)
  - [Modal analysis — a bracket](./tutorials/modal.md)
  - [Antenna FDTD — patch antenna](./tutorials/patch-antenna.md)

  # Validation gallery
  - [About the gallery](./validation/about.md)
  - [Ghia lid-driven cavity](./validation/ghia-cavity.md)
  - [Driver-Seegmiller backward-facing step](./validation/driver-seegmiller.md)
  - [NAFEMS linear-elasticity suite](./validation/nafems.md)
  - [GRI-Mech 3.0 adiabatic flame](./validation/gri-mech.md)
-->
