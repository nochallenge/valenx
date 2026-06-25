//! Aggregate outcome metrics computed from a finished scenario.
//!
//! Pure analysis over the final entity states and the per-pair first-detection
//! times — survivors per side, the total detection count, and the time to the
//! *first* detection anywhere in the run. This is reporting only; nothing here
//! influences the simulation.

use crate::entity::{Entity, Side};

/// Aggregate metrics summarising a scenario run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutcomeMetrics {
    /// Number of blue entities still alive at the stop time.
    pub survivors_blue: usize,
    /// Number of red entities still alive at the stop time.
    pub survivors_red: usize,
    /// Total number of distinct ordered (observer → target) first-detections.
    pub detection_count: usize,
    /// Time of the earliest detection anywhere in the scenario (s), or `None` if
    /// nothing was ever detected.
    pub time_to_first_detection_s: Option<f64>,
}

impl OutcomeMetrics {
    /// Compute the metrics from the final `entities`, the flattened
    /// `first_detect_time` table (`n × n`, indexed `observer * n + target`), and
    /// the entity count `n`.
    #[must_use]
    pub fn compute(entities: &[Entity], first_detect_time: &[Option<f64>], n: usize) -> Self {
        let survivors_blue = entities
            .iter()
            .filter(|e| e.alive && e.side == Side::Blue)
            .count();
        let survivors_red = entities
            .iter()
            .filter(|e| e.alive && e.side == Side::Red)
            .count();

        debug_assert_eq!(first_detect_time.len(), n * n);
        let mut detection_count = 0;
        let mut earliest: Option<f64> = None;
        for &t in first_detect_time.iter().take(n * n).flatten() {
            detection_count += 1;
            earliest = Some(match earliest {
                Some(e) => e.min(t),
                None => t,
            });
        }

        Self {
            survivors_blue,
            survivors_red,
            detection_count,
            time_to_first_detection_s: earliest,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::Mover;
    use nalgebra::Vector3;

    fn ent(side: Side, alive: bool) -> Entity {
        let mut e = Entity::new(Vector3::zeros(), side, Mover::Static, 0.0, 0.0, 0.0).unwrap();
        e.alive = alive;
        e
    }

    #[test]
    fn counts_survivors_and_detections() {
        let entities = vec![
            ent(Side::Blue, true),
            ent(Side::Blue, false),
            ent(Side::Red, true),
        ];
        let n = entities.len();
        let mut fdt = vec![None; n * n];
        let idx = |obs: usize, tgt: usize| obs * n + tgt;
        fdt[idx(0, 2)] = Some(5.0); // blue0 detected red2 at t=5
        fdt[idx(2, 0)] = Some(3.0); // red2 detected blue0 at t=3
        let m = OutcomeMetrics::compute(&entities, &fdt, n);
        assert_eq!(m.survivors_blue, 1);
        assert_eq!(m.survivors_red, 1);
        assert_eq!(m.detection_count, 2);
        assert_eq!(m.time_to_first_detection_s, Some(3.0));
    }

    #[test]
    fn no_detections_gives_none_ttfd() {
        let entities = vec![ent(Side::Blue, true)];
        let n = entities.len();
        let fdt = vec![None; n * n];
        let m = OutcomeMetrics::compute(&entities, &fdt, n);
        assert_eq!(m.detection_count, 0);
        assert!(m.time_to_first_detection_s.is_none());
    }
}
