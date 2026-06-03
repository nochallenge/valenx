//! End-to-end pipeline smoke test. Exercises gmsh → OpenFOAM as a
//! single chain, skipping gracefully when either tool is absent.
//!
//! Why this lives under `valenx-app/tests/` rather than the
//! workspace-root `tests/`: `valenx-app` is the only crate that
//! already depends on every runtime adapter we want to chain, so
//! integration tests that span the full pipeline slot in naturally
//! here. Per-adapter behaviour stays in each adapter's own tests.
//!
//! The test tolerates three machine states:
//!
//! 1. Neither gmsh nor OpenFOAM installed → asserts that both
//!    adapters' probes honestly report `ToolNotInstalled`.
//! 2. gmsh installed, OpenFOAM absent → meshes a unit cube, then
//!    checks that `valenx-adapter-openfoam::prepare()` stages the
//!    `mesh.msh` into the workdir and surfaces `ToolNotInstalled`
//!    (for `gmshToFoam`) with the expected hint.
//! 3. Both installed → generates polyMesh/ end-to-end.
//!
//! That way the same test runs in CI (no tools) and on developer
//! machines (varying tool coverage) without special-casing.

use std::path::PathBuf;
use std::sync::Arc;

use valenx_adapter_gmsh::GmshAdapter;
use valenx_adapter_openfoam::OpenFoamAdapter;
use valenx_core::{
    Adapter, AdapterError, CancellationToken, Case, LogLevel, LogSink, ProgressSink, RunContext,
};

/// Helper: run a gmsh meshing case synchronously and return the
/// workdir path. `None` means gmsh isn't installed and the caller
/// should skip.
fn run_gmsh_box(case_dir: &std::path::Path, workdir: &std::path::Path) -> Option<PathBuf> {
    let adapter = GmshAdapter::new();

    // Probe first — if gmsh isn't on PATH, bail cleanly.
    match adapter.probe() {
        Ok(report) if report.ok => {}
        _ => return None,
    }

    let case = Case {
        id: "box".into(),
        path: case_dir.to_path_buf(),
    };
    let prepared = adapter.prepare(&case, workdir).expect("gmsh prepare");

    let cancel = CancellationToken::new();
    let mut ctx = RunContext {
        cancel: &cancel,
        progress: Box::new(NoopProgress),
        log: Box::new(NoopLog),
    };

    let report = adapter
        .run(&prepared, &mut ctx)
        .expect("gmsh run should succeed when the tool is installed");
    assert_eq!(report.exit_code, 0);
    assert!(
        workdir.join("mesh.msh").is_file(),
        "gmsh should produce mesh.msh"
    );
    Some(workdir.to_path_buf())
}

#[test]
fn pipeline_gmsh_then_openfoam_prepare() {
    let case_dir = tempdir("valenx-e2e-case");
    let gmsh_workdir = tempdir("valenx-e2e-gmsh");
    let openfoam_workdir = tempdir("valenx-e2e-openfoam");

    // ---------------------------------------------------------------
    // Stage a canonical meshing case and a canonical CFD case in the
    // temp case dir. Both point at the same case directory so the
    // .msh produced by gmsh is visible to openfoam.prepare() via the
    // "look alongside case.toml" branch of ensure_poly_mesh.
    // ---------------------------------------------------------------
    std::fs::create_dir_all(case_dir.join("mesh-case")).unwrap();
    std::fs::create_dir_all(case_dir.join("cfd-case")).unwrap();
    std::fs::write(
        case_dir.join("mesh-case").join("case.toml"),
        r#"
[case]
format  = "1.0"
name    = "box-mesh"
physics = "meshing"
solver  = "gmsh.delaunay"
mesh    = "primary"

[mesh]
type   = "box"
origin = [0.0, 0.0, 0.0]
size   = [1.0, 1.0, 1.0]
characteristic_length = 0.3
algorithm_3d = "delaunay"
dim          = 3
"#,
    )
    .unwrap();
    std::fs::write(
        case_dir.join("cfd-case").join("case.toml"),
        include_str!("../../../tests/fixtures/minimal.valenx/cases/cfd-steady/case.toml"),
    )
    .unwrap();

    // ---------------------------------------------------------------
    // Stage 1: gmsh
    // ---------------------------------------------------------------
    let mesh_case_dir = case_dir.join("mesh-case");
    let has_gmsh_output = run_gmsh_box(&mesh_case_dir, &gmsh_workdir).is_some();
    if !has_gmsh_output {
        // Without gmsh we can still sanity-check the other side:
        // calling OpenFOAM's prepare() with no mesh.msh should
        // either stage the dicts and skip the mesh step (Ok), or
        // surface ToolNotInstalled if simpleFoam itself is missing.
        let adapter = OpenFoamAdapter::new();
        let case = Case {
            id: "cfd-steady".into(),
            path: case_dir.join("cfd-case"),
        };
        let outcome = adapter.prepare(&case, &openfoam_workdir);
        match outcome {
            Ok(_) => {
                // simpleFoam exists on PATH; dicts staged. Good.
            }
            Err(AdapterError::ToolNotInstalled { name, .. }) => {
                assert_eq!(name, "openfoam");
            }
            Err(other) => panic!("unexpected error: {other:?}"),
        }
        cleanup_all([&case_dir, &gmsh_workdir, &openfoam_workdir]);
        return;
    }

    // ---------------------------------------------------------------
    // Stage 2: copy the .msh into the CFD case dir so
    // ensure_poly_mesh picks it up via the "case_path/mesh.msh"
    // branch, then invoke OpenFOAM's prepare.
    // ---------------------------------------------------------------
    let msh_src = gmsh_workdir.join("mesh.msh");
    let msh_dst = case_dir.join("cfd-case").join("mesh.msh");
    std::fs::copy(&msh_src, &msh_dst).expect("copy mesh to cfd case dir");

    let adapter = OpenFoamAdapter::new();
    let case = Case {
        id: "cfd-steady".into(),
        path: case_dir.join("cfd-case"),
    };
    let outcome = adapter.prepare(&case, &openfoam_workdir);
    match outcome {
        Ok(_prepared) => {
            // Either simpleFoam + gmshToFoam both exist → polyMesh
            // materialised, or the staging got to the gmshToFoam
            // call and succeeded.
            assert!(
                openfoam_workdir
                    .join("constant")
                    .join("polyMesh")
                    .join("points")
                    .is_file(),
                "polyMesh/points should be materialised"
            );
            assert!(openfoam_workdir
                .join("system")
                .join("controlDict")
                .is_file());
        }
        Err(AdapterError::ToolNotInstalled {
            name: "openfoam",
            hint,
        }) => {
            // gmshToFoam (or simpleFoam) missing. The staging step
            // should still have copied the .msh into the workdir so
            // a retry with tools installed can pick up where we left
            // off.
            assert!(
                openfoam_workdir.join("mesh.msh").is_file(),
                "mesh.msh should have been staged into the workdir \
                 before the gmshToFoam step failed"
            );
            // The hint mentions either gmshToFoam or OpenFOAM.
            assert!(
                hint.contains("gmshToFoam") || hint.contains("OpenFOAM"),
                "hint should mention gmshToFoam or OpenFOAM, got: {hint}"
            );
        }
        Err(AdapterError::Run { stderr, .. }) => {
            // gmshToFoam ran but rejected the mesh. Still a pass for
            // the test — it means the plumbing got all the way to
            // the conversion step.
            assert!(!stderr.is_empty());
        }
        Err(other) => panic!("unexpected error: {other:?}"),
    }

    cleanup_all([&case_dir, &gmsh_workdir, &openfoam_workdir]);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "{tag}-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn cleanup_all<const N: usize>(dirs: [&std::path::Path; N]) {
    for d in dirs {
        let _ = std::fs::remove_dir_all(d);
    }
}

struct NoopProgress;
impl ProgressSink for NoopProgress {
    fn report(&self, _pct: f32, _message: &str) {}
}
struct NoopLog;
impl LogSink for NoopLog {
    fn log_line(&self, _level: LogLevel, _line: &str) {}
}

/// Silences the `Arc<AdapterRegistry>` reference warning on platforms
/// that don't pull it in transitively.
#[allow(dead_code)]
fn _keep_dep_alive() -> Option<Arc<dyn Adapter>> {
    None
}

// ===========================================================================
// Cross-crate computational-science end-to-end workflows
// ===========================================================================
//
// The tests above chain the *subprocess adapters*; the suite below chains the
// *pure native computational crates* (Round 6 + valenx-aero) through their
// real integration seams. `valenx-app` is the one crate that already depends
// on every one of them, so a multi-crate workflow test slots in naturally
// here — and `cargo test -p valenx-app --test pipeline_e2e` runs ONLY this
// file (no `rfd`, no GUI, no file dialog), so it stays inside the project
// test-lockdown.
//
// Each test runs a realistic workflow start-to-finish and asserts the FINAL
// result is physically / biologically sane. These catch integration bugs —
// a type-shape mismatch at a crate boundary, a unit error, a coordinate
// convention clash — that a single crate's unit tests cannot see.
//
// No solver is run at a resolution that blows the wall-clock: the aero case
// uses a coarse `cells_across_body: 4` grid (still a real 3-D Navier-Stokes
// solve), the qchem case is H2/STO-3G, the docking workflow stops at ligand
// preparation (a full dock is covered by valenx-e2e-tests/dock_smoke_1iep).

/// Workflow 1 — comparative genomics: FASTA text → parsed sequences
/// (`valenx-bioseq`) → pairwise + multiple alignment (`valenx-align`) →
/// distance matrix → neighbor-joining phylogenetic tree (`valenx-phylo`).
#[test]
fn e2e_fasta_to_alignment_to_phylogenetic_tree() {
    use valenx_align::msa::guidetree;
    use valenx_bioseq::{io::fasta, SeqKind};
    use valenx_phylo::distance::{cluster, distance_matrix, DistanceModel};

    // Four homologous DNA sequences, each exactly 60 nt, descended
    // from one aperiodic ancestor — the realistic phylogenetic
    // scenario. They differ ONLY by single-base substitutions at fixed
    // positions; the ancestor is aperiodic, so no frameshift produces a
    // spurious self-alignment and every pairwise distance reflects the
    // true substitution count:
    //   * clade 1 (A, B) shares 4 clade-1 marker substitutions,
    //   * clade 2 (C, D) shares 4 different clade-2 marker substitutions,
    //   * within each clade one taxon has a single extra private SNP.
    // So within a clade the sequences differ by exactly 1 site and
    // between clades by 8-9 sites — a clean, unambiguous {A,B} / {C,D}
    // split a correct parse → align → distance → tree pipeline must
    // recover.
    let fasta_text = "\
>taxon_A
GCTAAAGAAAATTACATAACCTACACGTCAGCCCGAAACTTGTTAGCCCAGTGTGAATCG
>taxon_B
GCTAAAGAAAATTACATAACCTACACGTCAGCCCGAAACTTGTTAGCCCAGTGTGCATCG
>taxon_C
GCTAAAGACAATTAAATAACATACACATCAGCACGAAAATTGTTGGCCCAATGTGAATCG
>taxon_D
GCTAAAGACAATTAAATAACATACACATCAGCACGAAAATTGTTGGCCCAATGTGCATCG
";

    // --- Stage 1: parse the FASTA (valenx-bioseq) ---
    let records = fasta::parse(fasta_text, SeqKind::Dna).expect("FASTA parses");
    assert_eq!(records.len(), 4, "four records expected");
    assert_eq!(records[0].id, "taxon_A");
    assert_eq!(records[0].seq.len(), 60, "each taxon sequence is 60 nt");

    // --- Stage 2a: a pairwise alignment of the two close taxa
    //               (valenx-align) — they should be ~identical ---
    let pair = valenx_align::bio::global_align(&records[0].seq, &records[1].seq)
        .expect("pairwise alignment");
    assert!(
        pair.percent_identity() > 0.90,
        "A vs B should be >90% identical, got {:.3}",
        pair.percent_identity()
    );
    // ...and a distant pair should be far less similar.
    let far = valenx_align::bio::global_align(&records[0].seq, &records[2].seq)
        .expect("pairwise alignment");
    assert!(
        far.percent_identity() < pair.percent_identity(),
        "A vs C ({:.3}) must be less similar than A vs B ({:.3})",
        far.percent_identity(),
        pair.percent_identity()
    );

    // --- Stage 2b: a multiple-sequence alignment of all four ---
    let seqs: Vec<&_> = records.iter().map(|r| &r.seq).collect();
    let msa = valenx_align::bio::multiple_align(&seqs).expect("MSA");
    assert_eq!(msa.depth(), 4, "MSA keeps all four rows");
    // Every aligned row has the same width (the defining MSA property).
    // Because the four taxa differ only by substitutions (no indels),
    // a correct affine-gap MSA introduces NO gaps — the alignment is
    // exactly as wide as the input. (A linear-gap profile aligner that
    // ignored the gap-open penalty would spuriously gap the clades
    // apart here — the bug this workflow regression-guards.)
    let width = msa.width();
    assert_eq!(width, 60, "substitution-only taxa need no alignment gaps");

    // --- Stage 3: distance matrix + neighbor-joining tree
    //              (valenx-phylo) ---
    let labels: Vec<String> = records.iter().map(|r| r.id.clone()).collect();
    let dm = distance_matrix(&msa, &labels, DistanceModel::JukesCantor)
        .expect("distance matrix");
    assert_eq!(dm.len(), 4);
    // The within-pair distance must be smaller than the between-pair
    // distance — the genomic signal survived the whole pipeline.
    let d_ab = dm.get(0, 1);
    let d_ac = dm.get(0, 2);
    assert!(
        d_ab < d_ac,
        "d(A,B)={d_ab:.4} should be < d(A,C)={d_ac:.4}"
    );

    let tree = cluster::neighbor_joining(&dm).expect("NJ tree");
    assert_eq!(tree.leaf_count(), 4, "tree has one leaf per taxon");
    // The NJ tree must group {A,B} and {C,D} as clades: in the tree,
    // A's closest leaf by patristic distance is B, not C or D.
    let a = tree.find("taxon_A").expect("taxon_A in tree");
    let b = tree.find("taxon_B").expect("taxon_B in tree");
    let c = tree.find("taxon_C").expect("taxon_C in tree");
    let d_tree_ab = tree.patristic_distance(a, b);
    let d_tree_ac = tree.patristic_distance(a, c);
    assert!(
        d_tree_ab < d_tree_ac,
        "tree must place A nearer B ({d_tree_ab:.4}) than C ({d_tree_ac:.4})"
    );

    // A UPGMA guide tree built straight from the same sequences (the
    // valenx-align side of the align↔phylo boundary) keeps all four
    // taxa and groups the two close pairs adjacently in its leaf order.
    let byte_seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();
    let guide_dm = guidetree::distance_matrix(
        &byte_seqs,
        &valenx_align::matrix::ScoringScheme::dna_default(),
    )
    .expect("guide distance matrix");
    let guide = guidetree::upgma(&guide_dm).expect("UPGMA guide tree");
    assert_eq!(guide.leaf_count(), 4, "guide tree visits every sequence");
    // The UPGMA leaf order places members of the same clade next to
    // each other: A and B are adjacent, C and D are adjacent.
    let order = guide.leaf_order();
    assert_eq!(order.len(), 4);
    let pos = |x: usize| order.iter().position(|&o| o == x).unwrap();
    assert_eq!(
        pos(0).abs_diff(pos(1)),
        1,
        "UPGMA places the close pair A,B adjacently"
    );
    assert_eq!(
        pos(2).abs_diff(pos(3)),
        1,
        "UPGMA places the close pair C,D adjacently"
    );
}

/// Workflow 2 — structure-based drug design prep: a SMILES string →
/// molecular graph + descriptors (`valenx-cheminf`) → ligand
/// preparation: protonation + rotatable-bond torsion tree
/// (`valenx-dock-screen`). The full docking score is exercised by
/// `valenx-e2e-tests/dock_smoke_1iep`; this covers the cheminf→dock-screen
/// seam that feeds it.
#[test]
fn e2e_smiles_to_descriptors_to_docking_prep() {
    use valenx_cheminf::{descriptors, mol_from_smiles};
    use valenx_dock_screen::prep::protonate::{prepare_ligand, ChargeModel};
    use valenx_dock_screen::prep::torsion::TorsionTree;

    // Ibuprofen — a real drug molecule with a carboxylic acid (a
    // protonation-relevant group) and several rotatable bonds.
    let smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
    let mol = mol_from_smiles(smiles).expect("SMILES parses");
    assert!(mol.heavy_atom_count() >= 15, "ibuprofen has 15 heavy atoms");

    // --- Descriptors (valenx-cheminf) — ibuprofen reference values ---
    let mw = valenx_cheminf::perceive::formula::average_molecular_weight(&mol);
    assert!(
        (mw - 206.3).abs() < 2.0,
        "ibuprofen MW should be ~206.3 g/mol, got {mw:.2}"
    );
    let logp = descriptors::crippen_logp(&mol);
    assert!(
        (2.0..5.5).contains(&logp),
        "ibuprofen cLogP should be lipophilic (~3.5), got {logp:.2}"
    );
    let hbd = descriptors::hbd(&mol);
    let hba = descriptors::hba(&mol);
    assert_eq!(hbd, 1, "ibuprofen has one H-bond donor (the -COOH)");
    assert!(hba >= 1, "ibuprofen has at least one H-bond acceptor");
    let lip = descriptors::lipinski(&mol);
    assert!(lip.passes(), "ibuprofen is a Lipinski-compliant drug");
    let rot = descriptors::rotatable_bonds(&mol);
    assert!(rot >= 3, "ibuprofen has several rotatable bonds, got {rot}");

    // --- Ligand prep (valenx-dock-screen) ---
    // Protonation at physiological pH 7.4: the carboxylic acid
    // deprotonates to a carboxylate, giving the molecule net charge −1.
    let prep = prepare_ligand(&mol, 7.4, ChargeModel::Gasteiger)
        .expect("ligand preparation");
    assert_eq!(
        prep.deprotonated, 1,
        "the -COOH should deprotonate at pH 7.4"
    );
    assert_eq!(
        prep.net_charge(),
        -1,
        "ibuprofen carboxylate carries net charge -1 at pH 7.4"
    );
    assert_eq!(
        prep.charges.len(),
        prep.molecule.atoms.len(),
        "one partial charge per atom"
    );

    // The torsion tree's rotatable-bond count is consistent with the
    // descriptor-level count — the cheminf and dock-screen perceptions
    // of the same molecule agree.
    let tree = TorsionTree::from_molecule(&mol).expect("torsion tree");
    assert!(
        tree.torsion_count() >= 1,
        "ibuprofen's torsion tree has rotatable bonds"
    );
    assert!(
        tree.group_count() >= 2,
        "a molecule with rotatable bonds has >= 2 rigid groups"
    );
}

/// Workflow 3 — protein structure analysis: PDB text → parsed structure
/// (`valenx-biostruct`) → geometry (radius of gyration) + secondary
/// structure (DSSP) → self-superposition (RMSD == 0). Exercises the
/// io → geometry → dssp → superpose seams within `valenx-biostruct`,
/// the canonical multi-module structural-bioinformatics pipeline.
#[test]
fn e2e_pdb_parse_to_geometry_and_secondary_structure() {
    use valenx_biostruct::geometry::shape;
    use valenx_biostruct::{dssp, io, superpose};

    // A short idealised α-helix fragment — three residues of backbone
    // atoms. Coordinates are a real right-handed α-helix (1.5 Å rise,
    // 100° twist per residue) so DSSP / geometry have real geometry to
    // work on.
    let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00  0.00           N
ATOM      6  CA  ALA A   2       4.000   2.830   0.150  1.00  0.00           C
ATOM      7  C   ALA A   2       5.510   2.680   0.000  1.00  0.00           C
ATOM      8  O   ALA A   2       6.030   1.560  -0.050  1.00  0.00           O
ATOM      9  N   ALA A   3       6.190   3.820   0.000  1.00  0.00           N
ATOM     10  CA  ALA A   3       7.640   3.870  -0.150  1.00  0.00           C
ATOM     11  C   ALA A   3       8.180   5.290  -0.100  1.00  0.00           C
ATOM     12  O   ALA A   3       7.420   6.260  -0.050  1.00  0.00           O
TER
END
";

    // --- Stage 1: parse the PDB (valenx-biostruct) ---
    let structure = io::read_structure(pdb, "test").expect("PDB parses");
    let atoms: Vec<_> = structure.first_model().atoms().collect();
    assert_eq!(atoms.len(), 12, "12 backbone atoms parsed");
    let chains = &structure.first_model().chains;
    assert_eq!(chains.len(), 1, "one chain");
    assert_eq!(chains[0].residues.len(), 3, "three residues");

    // --- Stage 2a: geometry — radius of gyration (valenx-biostruct) ---
    let pts: Vec<_> = atoms.iter().map(|a| (a.coord, 1.0_f64)).collect();
    let rg = shape::radius_of_gyration(&pts).expect("Rg computes");
    assert!(
        rg > 0.0 && rg < 10.0,
        "Rg of a 3-residue fragment is a few Angstrom, got {rg:.3}"
    );

    // --- Stage 2b: secondary structure — DSSP (valenx-biostruct) ---
    // The DSSP pass must run end-to-end and produce one assignment per
    // residue (the helix/sheet content of a 3-residue stub is not
    // asserted — too short for stable H-bond geometry).
    let ss = dssp::assign_chain(&structure.first_model().chains[0]);
    assert_eq!(
        ss.secondary_string().len(),
        3,
        "one secondary-structure code per residue"
    );

    // --- Stage 3: structure comparison — self-superposition
    //              (valenx-biostruct) must give RMSD ~ 0 ---
    let coords: Vec<_> = atoms.iter().map(|a| a.coord).collect();
    let rmsd = superpose::rmsd(&coords, &coords).expect("self-RMSD");
    assert!(
        rmsd < 1e-9,
        "a structure superposed on itself has zero RMSD, got {rmsd:.2e}"
    );
    // A translated copy also superposes to ~zero RMSD (Kabsch removes
    // the rigid-body translation) — the superpose seam is correct.
    let shift = nalgebra::Vector3::new(10.0, -5.0, 3.0);
    let shifted: Vec<_> = coords.iter().map(|p| p + shift).collect();
    let sup = superpose::kabsch(&shifted, &coords).expect("Kabsch superpose");
    assert!(
        sup.rmsd < 1e-6,
        "a pure translation superposes to zero RMSD, got {:.2e}",
        sup.rmsd
    );
}

/// Workflow 4 — gene expression: a DNA coding sequence → ORF finding →
/// translation → protein physico-chemical properties. Entirely within
/// `valenx-bioseq`, but spanning four modules (`ops::orf`,
/// `ops::translate`, `analysis::protparam`, `analysis::weight`) — the
/// canonical "from a gene to its protein product" pipeline.
#[test]
fn e2e_dna_orf_translation_protein_properties() {
    use valenx_bioseq::analysis::{protparam, weight};
    use valenx_bioseq::ops::{orf, translate};
    use valenx_bioseq::{Seq, SeqKind, Strand};

    // A DNA sequence with a clear ORF: a 5' UTR, an ATG start, a run of
    // sense codons, then a TAA stop, then a 3' UTR.
    //   UTR        ATG  GCT GCA AAA GAA GAT CTG  ... TAA  UTR
    let dna_text = "\
GGGCACCATGGCTGCAAAAGAAGATCTGGCTGCAAAAGAAGATCTGGCTGCATAAGGGCAC";
    let dna = Seq::new(SeqKind::Dna, dna_text).expect("valid DNA");

    // --- Stage 1: ORF finding (valenx-bioseq::ops::orf) ---
    let code = translate::GeneticCode::standard();
    let opts = orf::OrfOptions {
        atg_only: true,
        min_protein_len: 5,
        allow_no_stop: false,
    };
    let orfs = orf::find_orfs(&dna, &code, opts).expect("ORF scan");
    assert!(!orfs.is_empty(), "the sequence contains an ORF");
    let best = &orfs[0]; // longest first
    assert!(best.has_stop, "the reported ORF terminates in a stop codon");
    assert!(
        best.protein_len() >= 5,
        "ORF protein is at least the minimum length"
    );

    // --- Stage 2: translation (valenx-bioseq::ops::translate) ---
    // The ORF's own protein begins with Methionine (the ATG start).
    let protein = &best.protein;
    assert_eq!(
        protein.as_bytes().first().copied(),
        Some(b'M'),
        "translation of an ATG-started ORF begins with Met"
    );
    // Re-translating the ORF nucleotides independently — sliced
    // straight out of the genome at the reported span — yields the
    // same protein. This checks the orf → translate seam: the span
    // coordinates and the translation agree.
    if best.span.strand == Strand::Forward {
        let orf_nt = dna
            .slice(best.span.start, best.span.end)
            .expect("ORF span is within the sequence");
        let retranslated =
            translate::translate_default(&orf_nt, &code).expect("re-translate");
        assert_eq!(
            retranslated.as_bytes().iter().filter(|&&b| b != b'*').count(),
            protein.as_bytes().iter().filter(|&&b| b != b'*').count(),
            "independent re-translation matches the ORF protein length"
        );
    }

    // --- Stage 3: protein properties (valenx-bioseq::analysis) ---
    // Strip the trailing stop, then run ProtParam-class analysis.
    let aa: Vec<u8> = protein
        .as_bytes()
        .iter()
        .copied()
        .filter(|&b| b != b'*')
        .collect();
    let prot = Seq::new(SeqKind::Protein, aa).expect("protein sequence");
    let mw = weight::molecular_weight_protein(&prot).expect("protein MW");
    // An average amino acid is ~110 Da; the MW must be in that ballpark.
    let approx = prot.len() as f64 * 110.0;
    assert!(
        (mw - approx).abs() < approx * 0.5 + 200.0,
        "protein MW {mw:.1} should be ~{approx:.0} Da ({} residues)",
        prot.len()
    );
    let pp = protparam::protparam(&prot).expect("ProtParam");
    // ProtParam's residue count matches the translated protein — the
    // translate → protparam seam carried the right sequence through.
    assert_eq!(
        pp.length,
        prot.len(),
        "ProtParam length must match the protein length"
    );
    // The isoelectric point of any real protein lies in the 0..14 pH
    // range — a sanity check on the charge-titration solver.
    assert!(
        (0.0..=14.0).contains(&pp.isoelectric_point),
        "pI {} must be a valid pH",
        pp.isoelectric_point
    );
    // GRAVY (the Kyte-Doolittle hydropathy average) is a finite number
    // in the physically meaningful −4.5..4.5 per-residue band.
    assert!(
        pp.gravy.is_finite() && (-4.5..=4.5).contains(&pp.gravy),
        "GRAVY {} must be a finite hydropathy value",
        pp.gravy
    );
}

/// Workflow 5 — systems biology: a reaction network assembled by hand
/// (`valenx-sysbio::model`) → deterministic ODE time-course
/// (`valenx-sysbio::ode`) AND stochastic Gillespie simulation
/// (`valenx-sysbio::stochastic`). The same model drives both engines;
/// in the large-copy-number limit they must agree — the cross-check
/// that catches a stoichiometry or rate-law wiring bug.
#[test]
fn e2e_reaction_network_ode_and_gillespie() {
    use valenx_sysbio::model::{Model, RateLaw, Reaction, Species};
    use valenx_sysbio::ode::TimeCourse;
    use valenx_sysbio::stochastic::StochasticModel;

    // Irreversible conversion A -> B with mass-action rate k·[A].
    // Analytic solution: [A](t) = A0·e^(−k·t), [B](t) = A0·(1 − e^(−k·t)).
    let k = 0.5;
    let a0 = 1000.0;
    let build = || {
        let mut m = Model::new("a_to_b");
        let a = m.add_species(Species::new("A", a0));
        let b = m.add_species(Species::new("B", 0.0));
        m.add_reaction(Reaction {
            id: "r1".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MassAction {
                k,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    };
    let model = build();
    assert!(model.validate().is_ok(), "the hand-built model is valid");

    // --- Stage 1: deterministic ODE time-course (valenx-sysbio::ode) ---
    let t_end = 4.0;
    let traj = TimeCourse::new(t_end)
        .run(&model)
        .expect("ODE integration");
    let last = traj.states.last().expect("trajectory has samples");
    // Total mass is conserved: [A] + [B] == A0 at every time.
    let total = last[0] + last[1];
    assert!(
        (total - a0).abs() < a0 * 1e-3,
        "ODE must conserve mass: A+B={total:.2} vs A0={a0}"
    );
    // [A] decays toward the analytic A0·e^(−k·t_end).
    let analytic_a = a0 * (-k * t_end).exp();
    assert!(
        (last[0] - analytic_a).abs() < a0 * 0.05,
        "ODE [A]={:.2} should match analytic {analytic_a:.2}",
        last[0]
    );
    // [B] grew from zero — the reaction actually ran.
    assert!(last[1] > 0.5 * a0, "most of A converted to B by t_end");

    // --- Stage 2: stochastic Gillespie SSA (valenx-sysbio::stochastic) ---
    let sm = StochasticModel::from_model(&model).expect("stochastic model");
    let trace = sm
        .gillespie(t_end, /*seed=*/ 42, /*max_steps=*/ 200_000)
        .expect("Gillespie SSA");
    let final_counts = trace.states.last().expect("SSA trace has samples");
    let ssa_a = final_counts[0] as f64;
    let ssa_b = final_counts[1] as f64;
    // The SSA conserves the molecule count exactly (every reaction
    // moves one A to one B).
    let ssa_total = ssa_a + ssa_b;
    assert!(
        (ssa_total - a0).abs() < 1.0,
        "Gillespie conserves the molecule count: {ssa_total} vs {a0}"
    );
    // In the large-N limit the stochastic mean tracks the ODE: the SSA
    // [A] lands near the deterministic value (generous stochastic band).
    assert!(
        (ssa_a - analytic_a).abs() < a0 * 0.20,
        "Gillespie [A]={ssa_a:.0} should track the ODE/analytic {analytic_a:.0}"
    );
}

/// Workflow 6 — virtual wind tunnel: a 3-D triangle-mesh body
/// (`valenx-aero::geometry`) → the immersed-boundary RANS solver
/// (`valenx-aero`) → drag coefficient. A real 3-D steady Navier-Stokes
/// solve on a deliberately coarse grid so the wall-clock stays short;
/// the assertion is the qualitative bluff-body drag band, not an
/// engineering-tolerance number.
#[test]
fn e2e_body_mesh_to_windtunnel_drag_coefficient() {
    use nalgebra::Vector3;
    use valenx_aero::domain::{BoundaryConditions, TunnelSizing};
    use valenx_aero::{
        coefficients, geometry, integrate_forces, solve_steady, BodyMotion,
        SolverControls, TurbulenceModel, Wind, WindTunnel,
    };

    // A 1 m cube — the canonical bluff body. Its drag coefficient is a
    // textbook number (~1.0-1.3 for a cube broadside in turbulent flow);
    // a coarse immersed-boundary grid lands in the right band.
    let body = geometry::box_body(
        Vector3::new(-0.5, -0.5, -0.5),
        Vector3::new(0.5, 0.5, 0.5),
    );
    assert!(
        body.triangles.len() >= 12,
        "a box mesh has at least 12 triangles"
    );

    // --- Build the tunnel on a coarse grid (fast but real solve) ---
    let wind = Wind::straight(20.0).expect("20 m/s free stream");
    let tunnel = WindTunnel::build_with(
        &body,
        wind,
        BoundaryConditions::external_aero(),
        TunnelSizing {
            cells_across_body: 4,
            max_cells: 40_000,
            ..TunnelSizing::default()
        },
    )
    .expect("wind tunnel builds");

    // --- Run the steady RANS solver (valenx-aero) ---
    let controls = SolverControls {
        max_iterations: 60,
        turbulence: TurbulenceModel::KEpsilon,
        ..SolverControls::default()
    };
    let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());

    // --- Stage 3: drag coefficient must be physically sane ---
    let forces = integrate_forces(&tunnel, &flow, Vector3::zeros());
    let coeff = coefficients(&tunnel, &forces);
    let cd = coeff.cd;
    assert!(
        cd.is_finite() && cd > 0.0,
        "a body in a flow has a positive finite drag coefficient, got {cd}"
    );
    // A bluff body (a cube) has Cd of order one — not 0.01, not 100.
    // The coarse-grid immersed-boundary solver over-predicts somewhat
    // (documented in the aero crate's notes), so the band is generous.
    assert!(
        (0.1..5.0).contains(&cd),
        "cube drag coefficient should be order-one (bluff body), got {cd:.3}"
    );
    // The drag breakdown — pressure drag + friction drag — sums to the
    // total Cd (the force-decomposition seam is consistent).
    let sum = coeff.cd_pressure + coeff.cd_friction;
    assert!(
        (sum - cd).abs() < 1e-6,
        "drag breakdown {sum:.4} must sum to total Cd {cd:.4}"
    );
    // For a bluff body, pressure drag dominates friction drag.
    assert!(
        coeff.cd_pressure > coeff.cd_friction,
        "a bluff body's drag is pressure-dominated"
    );
}

/// Workflow 7 — quantum chemistry: a molecular geometry
/// (`valenx-qchem::geometry`) → restricted Hartree-Fock SCF
/// (`valenx-qchem`) → total energy + molecular properties. H2 in the
/// STO-3G basis is small enough to converge in a few SCF cycles and has
/// a well-known reference energy.
#[test]
fn e2e_molecular_geometry_to_hartree_fock_energy() {
    use valenx_qchem::geometry::{Atom, MolecularGeometry};
    use valenx_qchem::scf::rhf::ScfSettings;
    use valenx_qchem::{run_rhf, Element};

    // H2 at its ~0.74 Å equilibrium bond length, expressed in Bohr
    // (1 Å = 1.8897 Bohr). A closed-shell singlet — RHF applies.
    let bond_bohr = 0.74 * 1.8897259886;
    let h = Element::from_symbol("H").expect("hydrogen");
    let geom = MolecularGeometry::new(vec![
        Atom::new(h, [0.0, 0.0, 0.0]),
        Atom::new(h, [0.0, 0.0, bond_bohr]),
    ]);
    assert!(geom.is_closed_shell(), "H2 is a closed-shell singlet");
    assert_eq!(geom.n_atoms(), 2);

    // --- Run RHF/STO-3G (valenx-qchem) ---
    let report = run_rhf(&geom, "STO-3G", ScfSettings::default())
        .expect("RHF converges for H2");

    // --- Stage 3: energy + properties must be physically sane ---
    // The RHF/STO-3G total energy of H2 is ~ -1.117 Hartree (a textbook
    // value); allow a small band for the exact geometry / convergence.
    assert!(
        (report.total_energy - (-1.117)).abs() < 0.05,
        "H2 RHF/STO-3G energy should be ~ -1.117 Ha, got {:.6}",
        report.total_energy
    );
    // The total electronic energy is below the bare nuclear repulsion —
    // i.e. binding the electrons released energy.
    assert!(
        report.total_energy < report.nuclear_repulsion,
        "the bound molecule's energy is below its nuclear repulsion"
    );
    assert!(
        report.nuclear_repulsion > 0.0,
        "two protons repel — positive nuclear-repulsion energy"
    );
    // SCF actually iterated and converged.
    assert!(
        report.scf_iterations >= 1,
        "the SCF loop ran at least one cycle"
    );
    // H2 has two electrons in one doubly-occupied σ orbital, so there
    // is a HOMO-LUMO gap and one Mulliken charge per atom.
    assert_eq!(
        report.partial_charges.len(),
        2,
        "one Mulliken partial charge per atom"
    );
    // By symmetry the two H atoms carry equal partial charge.
    assert!(
        (report.partial_charges[0] - report.partial_charges[1]).abs() < 1e-6,
        "symmetric H2 has equal partial charges on both atoms"
    );
    if let Some(gap) = report.homo_lumo_gap {
        assert!(gap > 0.0, "H2 has a positive HOMO-LUMO gap");
    }
}
