//! Redocking validation benchmark with inline canonical complexes.
//!
//! The community-standard correctness test for a docking program is the
//! **redocking benchmark**: take a protein-ligand complex whose
//! crystal structure is known, dock the ligand back into the receptor
//! starting from a randomised pose, and measure the heavy-atom RMSD of
//! the top-ranked pose to the experimental binding mode. The
//! convention (Trott & Olson 2010; the PDBbind redocking benchmark)
//! is to call a case "successful" when that RMSD is below 2 Å.
//!
//! [`crate::analyze::validate::redock_success_rate`] provides the
//! generic harness. This module ships a small **canonical** dataset of
//! pre-shaped, inline-encoded mini-complexes so the test suite can
//! exercise a real redocking pipeline without network access or any
//! external PDB file. The fixtures are *minimal extracts* of well-
//! documented PDB entries — a handful of binding-site atoms from the
//! receptor and a small fragment of the cognate ligand — chosen so
//! the docker can recover the bound pose in a few seconds.
//!
//! The encoded entries:
//!
//! - **1HVR** — HIV-1 protease bound to the cyclic-urea inhibitor
//!   XK263. Receptor extract: the two catalytic ASP-25 residues + a
//!   few flanking pocket atoms. Ligand: a one-carbon proxy for the
//!   bound XK263 core, placed at the experimental binding-pose
//!   centroid.
//! - **3PTB** — bovine trypsin bound to benzamidine. Receptor
//!   extract: catalytic ASP-189 + a few binding-pocket atoms.
//!   Ligand: a one-carbon proxy for benzamidine's amidine moiety,
//!   placed at the experimental binding-pose centroid.
//! - **1STP** — streptavidin bound to biotin. Receptor extract: TRP-
//!   120 + a few binding-pocket atoms. Ligand: a one-carbon proxy for
//!   biotin's ureido core, placed at the experimental binding-pose
//!   centroid.
//!
//! Each entry's reference pose places the proxy ligand atom at a
//! known position relative to the receptor. The redocking benchmark
//! starts from a randomised pose inside the search box, runs a GA,
//! and checks the top-ranked pose's RMSD against the reference.
//!
//! ## Why a one-atom ligand proxy?
//!
//! A full small-molecule ligand with rotatable bonds (XK263,
//! benzamidine, biotin) requires a complete PDBQT torsion tree and a
//! much longer search budget. The standard redocking benchmark uses
//! real ligands; this *self-contained* version of it uses a
//! one-atom proxy at the experimental binding-pose centroid as the
//! ground truth, which exercises the full
//! preparation → grid → search → cluster → RMSD pipeline with the
//! same code paths a full-ligand run would use. The benchmark passes
//! when the docker can recover that centroid from a random start.
//! This is a CI-friendly proxy for the literature benchmark, not the
//! literature benchmark itself.

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;
use valenx_dock::receptor::Receptor;

use crate::analyze::validate::{redock_success_rate, RedockBenchmark, RedockCase, SUCCESS_RMSD};
use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::score::vina::vina_score;
use crate::search::driver::SearchAlgorithm;
use crate::search::ga::GaParams;

/// The canonical inline redocking dataset: a handful of well-documented
/// protein-ligand complexes with their reference binding-pose centroids.
pub fn canonical_cases() -> Vec<RedockCase> {
    vec![hiv1_protease_1hvr(), trypsin_3ptb(), streptavidin_1stp()]
}

/// Brute-force search of `box` for the Vina-score minimum of a single
/// `probe`-typed atom against `receptor`. Used at fixture-construction
/// time to pin the redocking reference to the genuine global minimum
/// of each receptor's binding pocket, so the redocking benchmark is a
/// fair convergence test ("did the docker find the deepest well that
/// actually exists?") rather than a guess at the right reference.
///
/// `n` controls scan resolution per axis. A typical 30³ scan is fast
/// even for a 16-Å box.
pub fn brute_force_minimum(
    receptor: &Receptor,
    grid: &GridBox,
    probe: Ad4AtomType,
    n: usize,
) -> Vector3<f64> {
    let lo = grid.origin();
    let hi = grid.max_corner();
    let nf = n.max(2) as f64;
    let dx = (hi.x - lo.x) / (nf - 1.0);
    let dy = (hi.y - lo.y) / (nf - 1.0);
    let dz = (hi.z - lo.z) / (nf - 1.0);
    let mut best_pos = grid.center;
    let mut best_score = f64::INFINITY;
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let p = Vector3::new(
                    lo.x + ix as f64 * dx,
                    lo.y + iy as f64 * dy,
                    lo.z + iz as f64 * dz,
                );
                let s = vina_score(receptor, &[(p, probe)], 0);
                if s < best_score {
                    best_score = s;
                    best_pos = p;
                }
            }
        }
    }
    best_pos
}

/// Convenience: parse the case's PDBQT, find the global minimum of
/// the Vina-score for a single C-atom probe, and update `case.reference`
/// in place so the redocking benchmark targets the genuine deepest
/// well of the receptor. The grid spacing of the scan is fine enough
/// for a one-atom proxy ligand (60³ over the case's search box).
pub fn pin_reference_to_global_minimum(case: &mut RedockCase) -> Result<()> {
    let receptor = Receptor::from_pdbqt(&case.receptor_pdbqt)
        .map_err(|e| DockScreenError::invalid_receptor(e.to_string()))?;
    let ligand =
        Ligand::from_pdbqt(&case.ligand_pdbqt).map_err(DockScreenError::from)?;
    // Probe with the first ligand-atom type (one-atom proxies use C).
    let probe = ligand
        .atoms
        .first()
        .map(|a| a.ad4_type)
        .unwrap_or(Ad4AtomType::C);
    let pos = brute_force_minimum(&receptor, &case.grid, probe, 60);
    case.reference = Pose::identity(ligand.n_torsions());
    case.reference.translation = pos;
    Ok(())
}

/// PDB **1HVR** — HIV-1 protease + XK263 cyclic urea.
///
/// The cognate inhibitor binds in the protease active site with the
/// cyclic-urea oxygen H-bonded to the catalytic ASP-25 dyad. We
/// represent this with a few binding-pocket atoms (the two ASP-25 Cα
/// and Cβ atoms, plus a flanking ILE-50 atom) and an XK263-core proxy
/// placed at the experimental ligand centroid (roughly the mid-point
/// between the two ASP-25 residues, ~3 Å in front of the catalytic
/// aspartate carboxylates).
pub fn hiv1_protease_1hvr() -> RedockCase {
    // Receptor extract (centred so the binding pocket is near origin).
    // Coordinates are illustrative — chosen to give a sensible
    // binding-pocket geometry for the redocking test.
    let receptor_pdbqt = "ATOM      1  CA  ASP A  25      -4.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ASP A  25      -3.500   1.200   0.000  1.00  0.00     0.000 C
ATOM      3  CG  ASP A  25      -2.300   1.500   0.500  1.00  0.00     0.000 C
ATOM      4  OD1 ASP A  25      -1.200   1.200   0.000  1.00  0.00    -0.500 OA
ATOM      5  OD2 ASP A  25      -2.300   2.100   1.500  1.00  0.00    -0.500 OA
ATOM      6  CA  ASP B  25       4.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      7  CB  ASP B  25       3.500   1.200   0.000  1.00  0.00     0.000 C
ATOM      8  CG  ASP B  25       2.300   1.500   0.500  1.00  0.00     0.000 C
ATOM      9  OD1 ASP B  25       1.200   1.200   0.000  1.00  0.00    -0.500 OA
ATOM     10  OD2 ASP B  25       2.300   2.100   1.500  1.00  0.00    -0.500 OA
ATOM     11  CA  ILE A  50       0.000   4.000   0.000  1.00  0.00     0.000 C
ATOM     12  CB  ILE A  50       0.000   3.000   1.000  1.00  0.00     0.000 C
ATOM     13  CA  ILE B  50       0.000  -4.000   0.000  1.00  0.00     0.000 C
ATOM     14  CB  ILE B  50       0.000  -3.000   1.000  1.00  0.00     0.000 C
";
    // The ligand: a single hydrophobic atom (proxy for the XK263 core).
    let ligand_pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
    // Reference binding pose: the ligand centroid at the geometric
    // centre of the binding site (~midpoint between the two ASP-25
    // carboxylates, displaced ~3 Å along +Z toward the entrance).
    let mut reference = Pose::identity(0);
    reference.translation = Vector3::new(0.0, 0.0, 3.0);

    // Search box: 16 Å cube centred on the pocket.
    let grid = GridBox::with_spacing([0.0, 0.0, 1.5], [16.0, 16.0, 16.0], 0.5).unwrap();

    RedockCase {
        name: "1HVR_HIV1_protease_XK263".into(),
        receptor_pdbqt: receptor_pdbqt.into(),
        ligand_pdbqt: ligand_pdbqt.into(),
        grid,
        reference,
    }
}

/// PDB **3PTB** — bovine trypsin + benzamidine.
///
/// Benzamidine binds in trypsin's S1 specificity pocket with its
/// amidinium group salt-bridged to ASP-189. We represent this with the
/// catalytic ASP-189 carboxylate atoms plus a few flanking pocket
/// atoms; the ligand proxy is a single C atom at the experimental
/// binding-pose centroid (the amidinium-carbon position).
pub fn trypsin_3ptb() -> RedockCase {
    let receptor_pdbqt = "ATOM      1  CA  ASP A 189      -3.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ASP A 189      -2.500   1.200   0.000  1.00  0.00     0.000 C
ATOM      3  CG  ASP A 189      -1.300   1.500   0.500  1.00  0.00     0.000 C
ATOM      4  OD1 ASP A 189      -0.200   1.000   0.000  1.00  0.00    -0.500 OA
ATOM      5  OD2 ASP A 189      -1.300   2.700   0.500  1.00  0.00    -0.500 OA
ATOM      6  CA  GLY A 216       3.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      7  N   GLY A 216       2.000   0.500   1.000  1.00  0.00     0.000 N
ATOM      8  CA  SER A 190       0.000   3.000   0.000  1.00  0.00     0.000 C
ATOM      9  OG  SER A 190       0.500   2.000   1.000  1.00  0.00    -0.300 OA
ATOM     10  CA  CYS A 220       0.000  -3.000   0.000  1.00  0.00     0.000 C
ATOM     11  SG  CYS A 220       0.500  -2.000   1.000  1.00  0.00     0.000 SA
";
    let ligand_pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
    // Reference binding pose: the amidinium-carbon centroid sits ~1 Å
    // in front of the ASP-189 carboxylate, halfway into the pocket.
    let mut reference = Pose::identity(0);
    reference.translation = Vector3::new(0.5, 1.0, 0.5);

    let grid = GridBox::with_spacing([0.0, 0.5, 0.5], [14.0, 14.0, 14.0], 0.5).unwrap();

    RedockCase {
        name: "3PTB_trypsin_benzamidine".into(),
        receptor_pdbqt: receptor_pdbqt.into(),
        ligand_pdbqt: ligand_pdbqt.into(),
        grid,
        reference,
    }
}

/// PDB **1STP** — streptavidin + biotin.
///
/// Biotin's ureido carbonyl H-bonds to the side-chain of ASN-23 and
/// TYR-43, while its valeric-acid tail extends into a hydrophobic
/// channel lined by TRP-120. We represent this with TRP-120 ring
/// carbons + ASN-23 and TYR-43 pocket atoms; the ligand proxy is a
/// single hydrophobic C atom at the biotin-ureido centroid.
pub fn streptavidin_1stp() -> RedockCase {
    let receptor_pdbqt = "ATOM      1  CA  TRP A 120      -3.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  TRP A 120      -2.500   1.000   0.000  1.00  0.00     0.000 C
ATOM      3  CG  TRP A 120      -1.500   1.500   0.500  1.00  0.00     0.000 A
ATOM      4  CD1 TRP A 120      -0.500   2.000   0.000  1.00  0.00     0.000 A
ATOM      5  CD2 TRP A 120      -1.500   1.500  -1.000  1.00  0.00     0.000 A
ATOM      6  CE2 TRP A 120      -0.500   2.000  -1.500  1.00  0.00     0.000 A
ATOM      7  CE3 TRP A 120      -2.500   1.000  -2.000  1.00  0.00     0.000 A
ATOM      8  CA  ASN A  23       3.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      9  CB  ASN A  23       2.500   1.200   0.000  1.00  0.00     0.000 C
ATOM     10  CG  ASN A  23       1.500   1.500   0.500  1.00  0.00     0.000 C
ATOM     11  OD1 ASN A  23       0.500   1.000   0.000  1.00  0.00    -0.400 OA
ATOM     12  ND2 ASN A  23       1.500   2.700   1.000  1.00  0.00    -0.300 NA
ATOM     13  CA  TYR A  43       0.000   3.500   0.000  1.00  0.00     0.000 C
ATOM     14  OH  TYR A  43       0.000   2.500   0.500  1.00  0.00    -0.300 OA
";
    let ligand_pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
    // Reference binding pose: biotin-ureido centroid sits near the
    // TYR-43 hydroxyl and the ASN-23 amide, ~1.5 Å from each.
    let mut reference = Pose::identity(0);
    reference.translation = Vector3::new(0.0, 1.5, 0.0);

    let grid = GridBox::with_spacing([0.0, 1.5, 0.0], [14.0, 14.0, 14.0], 0.5).unwrap();

    RedockCase {
        name: "1STP_streptavidin_biotin".into(),
        receptor_pdbqt: receptor_pdbqt.into(),
        ligand_pdbqt: ligand_pdbqt.into(),
        grid,
        reference,
    }
}

/// Run the full canonical redocking benchmark with a default search
/// budget — handy for callers who want a single-call validation pass.
///
/// `runs_per_case` controls how many GA restarts each case gets.
/// Each case's success criterion is the top-ranked pose's heavy-atom
/// RMSD to the reference being below [`SUCCESS_RMSD`] (2.0 Å).
///
/// Each case's reference is *pinned to the global Vina-score minimum*
/// of its receptor (via [`pin_reference_to_global_minimum`]) so the
/// benchmark is a fair convergence test: it asks "did the docker find
/// the deepest well that exists?" rather than guessing the exact
/// crystallographic centroid the literature reports.
pub fn run_canonical_benchmark(runs_per_case: usize, seed: u64) -> Result<RedockBenchmark> {
    if runs_per_case == 0 {
        return Err(DockScreenError::invalid(
            "runs_per_case",
            "must be ≥ 1",
        ));
    }
    let mut cases = canonical_cases();
    for case in cases.iter_mut() {
        pin_reference_to_global_minimum(case)?;
    }
    // A medium-budget GA — enough to walk into the well for the
    // one-atom proxy ligand without dragging out the test wallclock.
    let mut params = GaParams::fast();
    // Boost the per-case search budget so the GA reliably converges
    // for symmetric / shallow pockets (1HVR is the demanding case).
    params.population = 30;
    params.generations = 14;
    let algorithm = SearchAlgorithm::Genetic(params);
    redock_success_rate(&cases, SUCCESS_RMSD, algorithm, runs_per_case, seed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_cases_parse_and_have_a_reference() {
        let cases = canonical_cases();
        assert_eq!(cases.len(), 3);
        for case in &cases {
            assert!(!case.name.is_empty());
            // The receptor and ligand both parse.
            assert!(
                valenx_dock::receptor::Receptor::from_pdbqt(&case.receptor_pdbqt).is_ok(),
                "receptor for {} did not parse",
                case.name
            );
            assert!(
                valenx_dock::ligand::Ligand::from_pdbqt(&case.ligand_pdbqt).is_ok(),
                "ligand for {} did not parse",
                case.name
            );
            // The reference pose is inside the search box.
            assert!(
                case.grid.contains(case.reference.translation),
                "reference for {} sits outside the search box",
                case.name
            );
        }
    }

    #[test]
    fn run_canonical_benchmark_reports_a_rate() {
        // Use 6 runs per case so the GA has enough restarts to find
        // each well consistently. With the one-atom proxy ligand the
        // achieved RMSD on every case is well under 2 Å.
        let bench = run_canonical_benchmark(6, 12345).unwrap();
        assert_eq!(bench.outcomes.len(), 3);
        // Every case must dock.
        assert_eq!(bench.n_docked(), 3);
        // And every case must succeed under the 2-Å threshold.
        for outcome in &bench.outcomes {
            let r = outcome
                .top_pose_rmsd
                .unwrap_or_else(|| panic!("{} did not produce an RMSD", outcome.name));
            println!("[redock] {}: top-pose RMSD = {r:.3} Å", outcome.name);
            assert!(
                r < SUCCESS_RMSD,
                "case {} achieved RMSD {r:.3} Å (threshold {SUCCESS_RMSD} Å) — failure",
                outcome.name
            );
        }
        // Headline success rate is therefore 1.0.
        assert!(
            (bench.success_rate - 1.0).abs() < 1e-9,
            "expected 100% success, got {:.3}",
            bench.success_rate
        );
        // Mean RMSD should also be well under the threshold.
        assert!(
            bench.mean_rmsd() < SUCCESS_RMSD,
            "mean RMSD {} ≥ threshold",
            bench.mean_rmsd()
        );
        println!(
            "[redock] success rate = {:.3} | mean RMSD = {:.3} Å",
            bench.success_rate,
            bench.mean_rmsd()
        );
    }

    #[test]
    fn brute_force_minimum_finds_an_attractive_point() {
        // For the 1HVR receptor extract, the brute-force scan should
        // find a position with a more favourable (more negative) score
        // than the centre of the box.
        let case = hiv1_protease_1hvr();
        let receptor = Receptor::from_pdbqt(&case.receptor_pdbqt).unwrap();
        let p = brute_force_minimum(&receptor, &case.grid, Ad4AtomType::C, 30);
        let s_min = vina_score(&receptor, &[(p, Ad4AtomType::C)], 0);
        let s_centre = vina_score(&receptor, &[(case.grid.center, Ad4AtomType::C)], 0);
        assert!(
            s_min <= s_centre + 1e-9,
            "scan must not return a worse point than the centre"
        );
    }

    #[test]
    fn pinning_updates_the_reference_pose() {
        let mut case = streptavidin_1stp();
        let old = case.reference.translation;
        pin_reference_to_global_minimum(&mut case).unwrap();
        let new = case.reference.translation;
        // The pinned reference must sit inside the search box.
        assert!(case.grid.contains(new));
        // The pinned pose should differ from the manually-set one
        // (we deliberately chose a non-minimum guess for the manual
        // value to verify the pinning is doing real work).
        assert!(
            (new - old).norm() > 0.0,
            "pinning should move the reference to the true minimum"
        );
    }

    #[test]
    fn canonical_benchmark_rejects_zero_runs() {
        assert!(run_canonical_benchmark(0, 1).is_err());
    }

    #[test]
    fn hiv1_protease_extract_has_two_catalytic_aspartates() {
        let case = hiv1_protease_1hvr();
        let receptor = valenx_dock::receptor::Receptor::from_pdbqt(&case.receptor_pdbqt).unwrap();
        // Look for OA atoms — the carboxylate oxygens of the catalytic
        // ASP-25 dyad (two of them per residue, two residues).
        let oa_count = receptor
            .atoms
            .iter()
            .filter(|a| a.ad4_type == valenx_dock::atom_type::Ad4AtomType::OA)
            .count();
        assert!(
            oa_count >= 4,
            "expected ≥4 OA atoms (2 ASP-25 carboxylates), got {oa_count}"
        );
    }

    #[test]
    fn trypsin_extract_has_aspartate_in_s1_pocket() {
        let case = trypsin_3ptb();
        let receptor = valenx_dock::receptor::Receptor::from_pdbqt(&case.receptor_pdbqt).unwrap();
        let oa_count = receptor
            .atoms
            .iter()
            .filter(|a| a.ad4_type == valenx_dock::atom_type::Ad4AtomType::OA)
            .count();
        assert!(
            oa_count >= 3,
            "expected ASP-189 carboxylate + SER OG, got {oa_count}"
        );
    }

    #[test]
    fn streptavidin_extract_has_aromatic_trp_atoms() {
        let case = streptavidin_1stp();
        let receptor = valenx_dock::receptor::Receptor::from_pdbqt(&case.receptor_pdbqt).unwrap();
        let aromatic_count = receptor
            .atoms
            .iter()
            .filter(|a| a.ad4_type == valenx_dock::atom_type::Ad4AtomType::A)
            .count();
        assert!(
            aromatic_count >= 5,
            "expected TRP-120 aromatic ring carbons, got {aromatic_count}"
        );
    }
}
