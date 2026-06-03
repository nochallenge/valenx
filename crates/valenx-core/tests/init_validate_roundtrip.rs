//! Cross-binary integration test: every template `valenx-init`
//! understands must produce a project structure `valenx-validate`
//! accepts. Prevents regressions where the two CLIs drift out of
//! sync — a fix to the project loader that tightens validation
//! would fail this test if it broke a templated scaffold.
//!
//! Strategy: spawn `valenx-init` with each template against a fresh
//! temp dir, then spawn `valenx-validate` against the resulting
//! `.valenx` directory and assert exit code 0.

use std::path::PathBuf;
use std::process::Command;

fn init_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-init"))
}

fn validate_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-validate"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "valenx-init-validate-{tag}-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    // Don't pre-create — `valenx-init` creates the dir itself, and
    // a sibling test (`refuses_to_overwrite_existing_project`)
    // depends on init's first-run-only behaviour.
    d
}

/// Templates exhaustively. Mirrors the `Template` enum in
/// `bin/valenx_init.rs`. Adding a new template here without also
/// adding it to the binary's enum (or vice-versa) would surface as
/// a test failure on the next CI run.
const TEMPLATES: &[&str] = &[
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
    // Biology — alignment toolkit (Phase 18)
    "bwa",
    "minimap2",
    "mafft",
    "muscle",
    "hmmer",
    "samtools",
    // Biology — structure prediction expansion (Phase 17.5)
    "esmfold",
    "openfold",
    "alphafold2",
    "alphafold3",
    // Biology — variant calling (Phase 19)
    "bcftools",
    "gatk",
    "deepvariant",
    // Biology — viewers (Phase 23)
    "pymol",
    "vmd",
    "igv",
    // Biology — protein design (Phase 27)
    "rfdiffusion",
    "proteinmpnn",
    // Biology — molecular docking (Phase 34)
    "vina",
    "autodock4",
    // Biology — cheminformatics expansion (Phase 24)
    "deepchem",
    "openbabel",
    "avogadro",
    // Biology — workflow managers (Phase 22)
    "nextflow",
    "snakemake",
    // Biology — single-cell genomics (Phase 19.5)
    "scanpy",
    "scvi",
    // Biology — protein design expansion (Phase 27.5)
    "chroma",
    "esm-if",
    "rfantibody",
    // Biology — aligners expansion (Phase 18.5)
    "bowtie2",
    "mmseqs2",
    "diamond",
    // Biology — RNA-seq alignment (Phase 18.6)
    "hisat2",
    "star",
    // Biology — transcript quantification (Phase 20)
    "salmon",
    "kallisto",
    // Biology — phylogenetics (Phase 30)
    "iqtree",
    "raxml-ng",
    "fasttree",
    // Biology — RNA structure (Phase 28)
    "viennarna",
    "rnastructure",
    "nupack",
    // Biology — quantum chemistry (Phase 25)
    "psi4",
    "nwchem",
    "xtb",
    // Biology — EvolutionaryScale models (Phase 27.6)
    "esm3",
    "esmc",
    // Biology — systems biology (Phase 32)
    "copasi",
    "bionetgen",
    "physicell",
    // Biology — cryo-EM (Phase 36)
    "relion",
    "eman2",
    "ctffind",
    // Biology — sequencing read simulators (Phase 31)
    "art",
    "wgsim",
    "badread",
    // Biology — CRISPR design (Phase 35)
    "chopchop",
    "crispor",
    "cas-offinder",
    // Biology — Rosetta family (Phase 38)
    "rosetta",
    "pyrosetta",
    // Biology — population genetics (Phase 29)
    "slim",
    "msprime",
    "tskit",
    // Biology — Bayesian phylogenetics (Phase 30.5)
    "beast2",
    "mrbayes",
    // Biology — DNA structural geometry (Phase 39)
    "x3dna",
    "curves",
    "dssr",
    // Biology — MD analysis expansion (Phase 5.5)
    "plumed",
    "prody",
    "cpptraj",
    // Biology — synthetic biology (Phase 33)
    "pysbol",
    "j5",
    "cello",
    // Biology — alignment toolkit expansion (Phase 18.7)
    "blast",
    "clustalo",
    "tcoffee",
    // Biology — single-cell genomics expansion (Phase 19.6)
    "seurat",
    "anndata",
    // Biology — bio MD engines (Phase 5.6)
    "namd",
    "sander",
    "hoomd",
    // Biology — MD analysis sister (Phase 5.7)
    "mdtraj",
    // Biology — structure prediction + search (Phase 17.7)
    "rosettafold",
    "omegafold",
    "foldseek",
    // Biology — spatial stochastic reaction-diffusion (Phase 32.5)
    "smoldyn",
    "mcell",
    // Biology — sequence editors / plasmid design (Phase 41)
    "pydna",
    "jalview",
    // Biology — microscopy / bioimage analysis (Phase 40)
    "fiji",
    "cellprofiler",
    "ilastik",
    // Biology — workflow expansion (Phase 22.5)
    "planemo",
    "cromwell",
    "cwltool",
    // Biology — web 3D molecular visualization (Phase 42)
    "molstar",
    "ngl",
    // Biology — mRNA design (Phase 43)
    "dnachisel",
    "lineardesign",
    "icodon",
    // Biology — RNA folding expansion (Phase 44.5)
    "mfold",
    "eternafold",
    "linearfold",
    // Biology — base + prime editing design (Phase 35.5)
    "be-designer",
    "be-hive",
    "primedesign",
    "pegfinder",
    // Biology — edit-outcome prediction (Phase 35.6)
    "indelphi",
    "forecast",
    "alphamissense",
    "crispritz",
    // Biology — pharmacokinetics + RNA tertiary (Phase 45)
    "pksim",
    "simrna",
];

#[test]
fn every_template_scaffolds_a_validate_clean_project() {
    for template in TEMPLATES {
        let dir = tempdir(template);
        // 1. init.
        let init_out = Command::new(init_binary())
            .arg(&dir)
            .arg("--template")
            .arg(template)
            .output()
            .expect("spawn init");
        assert!(
            init_out.status.success(),
            "init failed for template `{template}`; stderr: {}",
            String::from_utf8_lossy(&init_out.stderr)
        );
        // 2. validate.
        let val_out = Command::new(validate_binary())
            .arg(&dir)
            .output()
            .expect("spawn validate");
        let stderr = String::from_utf8_lossy(&val_out.stderr).into_owned();
        let stdout = String::from_utf8_lossy(&val_out.stdout).into_owned();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(
            val_out.status.success(),
            "validate failed for template `{template}`\nstdout: {stdout}\nstderr: {stderr}"
        );
    }
}

#[test]
fn cfd_template_validate_output_lists_cavity_case() {
    // Specific shape check on the CFD template: the case dir is
    // `cavity` and validate's text output should call it out by
    // name. Catches the regression where init wrote
    // `cases.order = ["case-1"]` but the case actually landed in
    // `cases/cavity/` — pre-9dfa5f2 behaviour.
    let dir = tempdir("cfd-shape");
    let init_out = Command::new(init_binary())
        .arg(&dir)
        .arg("--template")
        .arg("cfd")
        .output()
        .expect("spawn init");
    assert!(init_out.status.success());

    let val_out = Command::new(validate_binary())
        .arg(&dir)
        .output()
        .expect("spawn validate");
    let stdout = String::from_utf8_lossy(&val_out.stdout).into_owned();
    let _ = std::fs::remove_dir_all(&dir);
    assert!(val_out.status.success(), "validate failed: {stdout}");
    assert!(
        stdout.contains("cavity"),
        "validate text mode should list `cavity` case; got: {stdout}"
    );
}

#[test]
fn list_templates_flag_lists_every_canonical_name() {
    // `--list-templates` is a quick discovery flag. Spawn the binary
    // and assert every template name we expect surfaces in stdout.
    // Mirrors the in-binary unit test but exercises the compiled
    // binary end-to-end so the full `match` dispatch + print path
    // is covered.
    let out = Command::new(init_binary())
        .arg("--list-templates")
        .output()
        .expect("spawn init");
    assert!(
        out.status.success(),
        "list-templates failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    for canonical in TEMPLATES {
        assert!(
            stdout.contains(canonical),
            "stdout missing `{canonical}`:\n{stdout}"
        );
    }
}

#[test]
fn fea_template_json_envelope_lists_cantilever_case() {
    let dir = tempdir("fea-json");
    let init_out = Command::new(init_binary())
        .arg(&dir)
        .arg("--template")
        .arg("fea")
        .output()
        .expect("spawn init");
    assert!(init_out.status.success());

    let val_out = Command::new(validate_binary())
        .arg(&dir)
        .arg("--format")
        .arg("json")
        .output()
        .expect("spawn validate");
    let stdout = String::from_utf8_lossy(&val_out.stdout).into_owned();
    let _ = std::fs::remove_dir_all(&dir);
    assert!(val_out.status.success(), "validate failed: {stdout}");
    let v: serde_json::Value =
        serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not JSON: {e}\n{stdout}"));
    assert_eq!(v["ok"], serde_json::Value::Bool(true));
    let cases = v["cases"].as_array().expect("cases array");
    assert!(
        cases.iter().any(|c| c["name"] == "cantilever"),
        "expected `cantilever` in cases array; got: {cases:?}"
    );
}
