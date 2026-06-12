//! Feature 24 — docking result I/O.
//!
//! AutoDock and Vina write their docked poses as a multi-`MODEL`
//! PDBQT file: each pose is a `MODEL … ENDMDL` block, and a
//! `REMARK VINA RESULT:` line carries the pose's score. This module
//! reads and writes that format and provides a tabular summary.
//!
//! - [`read_docking_result`] parses a multi-model PDBQT result file
//!   into a [`DockingResultTable`] — one [`ResultPose`] per `MODEL`,
//!   with the parsed score and atom coordinates.
//! - [`write_docking_result`] writes a [`DockingResultTable`] back to
//!   the same format, reusing [`valenx_bio`]'s
//!   [`valenx_bio::format::pdbqt::write_pose_ensemble`].
//! - [`DockingResultTable`] is the in-memory results table — ranked
//!   poses with scores and a CSV / text-table renderer.
//!
//! The pose *geometry* round-trips: every atom of the template (its
//! name, charge, AD4 type) is preserved; only the coordinates and
//! score vary between models.

use nalgebra::Vector3;

use valenx_bio::format::pdbqt::{parse, write_pose_ensemble, PdbqtAtom, PdbqtRecord};

use crate::error::{DockScreenError, Result};

/// One pose in a docking result: a score plus per-atom coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct ResultPose {
    /// 1-based model number from the `MODEL` record.
    pub model: usize,
    /// The pose's docking score (kcal/mol) from the
    /// `REMARK VINA RESULT:` line; `0.0` if the remark was absent.
    pub score: f64,
    /// World-space coordinates of every atom, in file order.
    pub coordinates: Vec<Vector3<f64>>,
}

/// An in-memory docking-result table — the template atoms plus a
/// ranked list of poses.
#[derive(Clone, Debug, PartialEq)]
pub struct DockingResultTable {
    /// The ligand's atom records (names, charges, AD4 types). Pose
    /// coordinates are applied onto these for output.
    pub template_atoms: Vec<PdbqtAtom>,
    /// The docked poses.
    pub poses: Vec<ResultPose>,
}

impl DockingResultTable {
    /// Number of poses in the table.
    pub fn len(&self) -> usize {
        self.poses.len()
    }

    /// `true` if the table has no poses.
    pub fn is_empty(&self) -> bool {
        self.poses.is_empty()
    }

    /// The best (lowest-score) pose, if any.
    pub fn best(&self) -> Option<&ResultPose> {
        self.poses.iter().min_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Sort the poses ascending by score (best first).
    pub fn sort_by_score(&mut self) {
        self.poses.sort_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Render the table as CSV: a `rank,model,score` header and one
    /// row per pose, in current table order.
    pub fn to_csv(&self) -> String {
        let mut s = String::from("rank,model,score\n");
        for (i, p) in self.poses.iter().enumerate() {
            s.push_str(&format!("{},{},{:.3}\n", i + 1, p.model, p.score));
        }
        s
    }

    /// Render the table as a fixed-width text table for the console.
    pub fn to_text_table(&self) -> String {
        let mut s = String::from(" rank  model      score\n");
        s.push_str("-----  -----  ---------\n");
        for (i, p) in self.poses.iter().enumerate() {
            s.push_str(&format!("{:>5}  {:>5}  {:>9.3}\n", i + 1, p.model, p.score));
        }
        s
    }
}

/// Feature 24 (read) — parse a multi-model PDBQT docking-result file.
///
/// Each `MODEL … ENDMDL` block becomes a [`ResultPose`]; the score is
/// lifted from a `REMARK VINA RESULT:` line if present. The first
/// model's atom records become the table's template atoms.
///
/// Returns [`DockScreenError::Parse`] if the file has no `MODEL`
/// blocks or no atoms.
pub fn read_docking_result(text: &str) -> Result<DockingResultTable> {
    // The pdbqt parser surfaces ATOM records and ROOT / BRANCH
    // markers but treats MODEL / REMARK as `Other`. So we scan
    // line-by-line for the model structure and hand each model's atom
    // lines to the pdbqt parser.
    let mut poses: Vec<ResultPose> = Vec::new();
    let mut template: Vec<PdbqtAtom> = Vec::new();

    let mut current_model: Option<usize> = None;
    let mut current_score = 0.0;
    let mut current_atom_lines: Vec<&str> = Vec::new();

    let flush = |model: usize,
                 score: f64,
                 atom_lines: &[&str],
                 poses: &mut Vec<ResultPose>,
                 template: &mut Vec<PdbqtAtom>|
     -> Result<()> {
        let block = atom_lines.join("\n");
        let records = parse(&block)
            .map_err(|e| DockScreenError::parse("pdbqt", format!("model {model}: {e}")))?;
        let atoms: Vec<PdbqtAtom> = records
            .into_iter()
            .filter_map(|r| match r {
                PdbqtRecord::Atom(a) => Some(a),
                _ => None,
            })
            .collect();
        if atoms.is_empty() {
            return Err(DockScreenError::parse(
                "pdbqt",
                format!("model {model} has no ATOM records"),
            ));
        }
        if template.is_empty() {
            *template = atoms.clone();
        }
        poses.push(ResultPose {
            model,
            score,
            coordinates: atoms.iter().map(|a| a.position).collect(),
        });
        Ok(())
    };

    for line in text.lines() {
        let trimmed = line.trim_start();
        if let Some(rest) = trimmed.strip_prefix("MODEL") {
            current_model = Some(rest.trim().parse::<usize>().unwrap_or(poses.len() + 1));
            current_score = 0.0;
            current_atom_lines.clear();
        } else if trimmed.starts_with("ENDMDL") {
            if let Some(m) = current_model.take() {
                flush(
                    m,
                    current_score,
                    &current_atom_lines,
                    &mut poses,
                    &mut template,
                )?;
            }
            current_atom_lines.clear();
        } else if let Some(score) = parse_vina_remark(trimmed) {
            current_score = score;
        } else if trimmed.starts_with("ATOM") || trimmed.starts_with("HETATM") {
            current_atom_lines.push(line);
        }
        // ROOT / BRANCH / TORSDOF lines are ignored — coordinates are
        // all this reader needs.
    }
    // A file with a single un-MODEL'd pose (a bare PDBQT) is still
    // valid — treat the trailing atom block as model 1.
    if poses.is_empty() && !current_atom_lines.is_empty() {
        flush(
            1,
            current_score,
            &current_atom_lines,
            &mut poses,
            &mut template,
        )?;
    }
    if poses.is_empty() {
        return Err(DockScreenError::parse(
            "pdbqt",
            "no MODEL blocks or ATOM records found in docking result",
        ));
    }
    Ok(DockingResultTable {
        template_atoms: template,
        poses,
    })
}

/// Parse a `REMARK VINA RESULT:` line and return the score, or `None`
/// if the line is not a Vina-result remark.
fn parse_vina_remark(line: &str) -> Option<f64> {
    let idx = line.find("VINA RESULT:")?;
    let rest = &line[idx + "VINA RESULT:".len()..];
    rest.split_whitespace().next()?.parse::<f64>().ok()
}

/// Feature 24 (write) — serialise a docking-result table to a
/// multi-model PDBQT string.
///
/// Each pose becomes a `MODEL … ENDMDL` block with a
/// `REMARK VINA RESULT:` score line, written via [`valenx_bio`]'s
/// ensemble writer.
///
/// Returns [`DockScreenError::Invalid`] if the table has no template
/// atoms, or if any pose's coordinate count does not match the
/// template atom count.
pub fn write_docking_result(table: &DockingResultTable) -> Result<String> {
    if table.template_atoms.is_empty() {
        return Err(DockScreenError::invalid(
            "template_atoms",
            "cannot write a docking result with no template atoms",
        ));
    }
    let n = table.template_atoms.len();
    for p in &table.poses {
        if p.coordinates.len() != n {
            return Err(DockScreenError::invalid(
                "coordinates",
                format!(
                    "pose (model {}) has {} coordinates, template has {n} atoms",
                    p.model,
                    p.coordinates.len()
                ),
            ));
        }
    }
    let ensemble: Vec<(Vec<Vector3<f64>>, f64)> = table
        .poses
        .iter()
        .map(|p| (p.coordinates.clone(), p.score))
        .collect();
    Ok(write_pose_ensemble(&table.template_atoms, &ensemble))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A two-model PDBQT docking result.
    const RESULT: &str = "MODEL 1
REMARK VINA RESULT:      -8.500      0.000      0.000
ROOT
ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00     0.100 C
ATOM      2  O1  LIG A   1       2.500   2.000   3.000  1.00  0.00    -0.300 OA
ENDROOT
TORSDOF 0
ENDMDL
MODEL 2
REMARK VINA RESULT:      -7.200      1.500      2.100
ROOT
ATOM      1  C1  LIG A   1       5.000   6.000   7.000  1.00  0.00     0.100 C
ATOM      2  O1  LIG A   1       6.500   6.000   7.000  1.00  0.00    -0.300 OA
ENDROOT
TORSDOF 0
ENDMDL
";

    #[test]
    fn reads_a_two_model_result() {
        let table = read_docking_result(RESULT).unwrap();
        assert_eq!(table.len(), 2);
        assert_eq!(table.template_atoms.len(), 2);
        // Model 1's score and coordinates.
        assert!((table.poses[0].score - -8.5).abs() < 1e-6);
        assert!((table.poses[0].coordinates[0].x - 1.0).abs() < 1e-6);
        // Model 2's score.
        assert!((table.poses[1].score - -7.2).abs() < 1e-6);
    }

    #[test]
    fn best_pose_is_the_lowest_score() {
        let table = read_docking_result(RESULT).unwrap();
        let best = table.best().unwrap();
        assert!((best.score - -8.5).abs() < 1e-6);
        assert_eq!(best.model, 1);
    }

    #[test]
    fn rejects_a_file_with_no_models_or_atoms() {
        assert!(read_docking_result("REMARK just a comment\n").is_err());
        assert!(read_docking_result("").is_err());
    }

    #[test]
    fn write_then_read_roundtrips() {
        let table = read_docking_result(RESULT).unwrap();
        let text = write_docking_result(&table).unwrap();
        let reparsed = read_docking_result(&text).unwrap();
        assert_eq!(reparsed.len(), table.len());
        // Scores survive the round-trip.
        assert!((reparsed.poses[0].score - table.poses[0].score).abs() < 1e-3);
        // Coordinates survive the round-trip.
        assert!((reparsed.poses[1].coordinates[0] - table.poses[1].coordinates[0]).norm() < 1e-2);
    }

    #[test]
    fn write_rejects_a_pose_with_wrong_atom_count() {
        let mut table = read_docking_result(RESULT).unwrap();
        // Corrupt model 2: drop a coordinate.
        table.poses[1].coordinates.pop();
        assert!(write_docking_result(&table).is_err());
    }

    #[test]
    fn write_rejects_an_empty_template() {
        let table = DockingResultTable {
            template_atoms: vec![],
            poses: vec![],
        };
        assert!(write_docking_result(&table).is_err());
    }

    #[test]
    fn csv_and_text_table_have_a_row_per_pose() {
        let table = read_docking_result(RESULT).unwrap();
        let csv = table.to_csv();
        // Header + 2 rows.
        assert_eq!(csv.lines().count(), 3);
        assert!(csv.contains("rank,model,score"));
        assert!(csv.contains("-8.500"));
        let txt = table.to_text_table();
        // Two header lines + 2 rows.
        assert_eq!(txt.lines().count(), 4);
    }

    #[test]
    fn sort_by_score_orders_best_first() {
        let mut table = read_docking_result(RESULT).unwrap();
        // Reverse first so the sort has work to do.
        table.poses.reverse();
        table.sort_by_score();
        assert!(table.poses[0].score <= table.poses[1].score);
    }

    #[test]
    fn reads_a_bare_single_pose_pdbqt() {
        // No MODEL wrapper — a plain PDBQT ligand counts as model 1.
        let bare = "ROOT
ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00     0.100 C
ENDROOT
TORSDOF 0
";
        let table = read_docking_result(bare).unwrap();
        assert_eq!(table.len(), 1);
        assert_eq!(table.poses[0].model, 1);
    }
}
