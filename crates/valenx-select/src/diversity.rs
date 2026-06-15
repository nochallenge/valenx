//! Diversity-aware selection so a shortlist "fails differently".

use crate::error::SelectError;

/// Euclidean distance between two equal-length feature vectors.
pub fn euclidean(a: &[f64], b: &[f64]) -> Result<f64, SelectError> {
    if a.len() != b.len() {
        return Err(SelectError::Inconsistent {
            what: "feature dimension",
        });
    }
    let mut sum = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        if !x.is_finite() || !y.is_finite() {
            return Err(SelectError::NonFinite { what: "feature" });
        }
        sum += (x - y) * (x - y);
    }
    Ok(sum.sqrt())
}

fn validate_features(features: &[Vec<f64>]) -> Result<usize, SelectError> {
    if features.is_empty() {
        return Err(SelectError::Empty { what: "features" });
    }
    let dim = features[0].len();
    if dim == 0 {
        return Err(SelectError::Empty {
            what: "feature vector",
        });
    }
    for f in features {
        if f.len() != dim {
            return Err(SelectError::Inconsistent {
                what: "feature dimension",
            });
        }
    }
    Ok(dim)
}

/// MaxMin (farthest-point) sampling: start from `start`, then repeatedly add
/// the point with the greatest distance to the nearest already-selected point.
/// Returns up to `n` indices that are maximally spread out.
pub fn farthest_point_select(
    features: &[Vec<f64>],
    n: usize,
    start: usize,
) -> Result<Vec<usize>, SelectError> {
    validate_features(features)?;
    if n == 0 {
        return Err(SelectError::ZeroN);
    }
    if start >= features.len() {
        return Err(SelectError::StartOutOfRange {
            index: start,
            len: features.len(),
        });
    }
    let total = features.len();
    let mut selected = vec![start];
    while selected.len() < n && selected.len() < total {
        let mut best_idx = None;
        let mut best_dist = -1.0_f64;
        for i in 0..total {
            if selected.contains(&i) {
                continue;
            }
            let mut nearest = f64::INFINITY;
            for &s in &selected {
                let d = euclidean(&features[i], &features[s])?;
                if d < nearest {
                    nearest = d;
                }
            }
            if nearest > best_dist {
                best_dist = nearest;
                best_idx = Some(i);
            }
        }
        match best_idx {
            Some(i) => selected.push(i),
            None => break,
        }
    }
    Ok(selected)
}

/// Score-ordered sphere exclusion (Butina-style cluster-and-top): walk
/// candidates best-`scores`-first and accept one only if it lies strictly
/// beyond `radius` of every already-accepted point. Returns up to `n` indices —
/// the best representative of each well-separated cluster.
///
/// Selecting more (`n+1`) extends the previous selection (prefix-stable).
pub fn sphere_exclusion_select(
    features: &[Vec<f64>],
    scores: &[f64],
    n: usize,
    radius: f64,
) -> Result<Vec<usize>, SelectError> {
    validate_features(features)?;
    if scores.len() != features.len() {
        return Err(SelectError::Inconsistent {
            what: "scores length",
        });
    }
    if n == 0 {
        return Err(SelectError::ZeroN);
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err(SelectError::NonPositiveRadius { value: radius });
    }
    for &s in scores {
        if !s.is_finite() {
            return Err(SelectError::NonFinite { what: "score" });
        }
    }

    let mut order: Vec<usize> = (0..features.len()).collect();
    order.sort_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });

    let mut accepted: Vec<usize> = Vec::new();
    for &i in &order {
        if accepted.len() == n {
            break;
        }
        let mut ok = true;
        for &a in &accepted {
            if euclidean(&features[i], &features[a])? <= radius {
                ok = false;
                break;
            }
        }
        if ok {
            accepted.push(i);
        }
    }
    Ok(accepted)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line() -> Vec<Vec<f64>> {
        (0..=10).map(|i| vec![i as f64]).collect()
    }

    #[test]
    fn euclidean_known() {
        assert!((euclidean(&[0.0, 0.0], &[3.0, 4.0]).unwrap() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn farthest_point_spreads_out() {
        // 11 points on a line; from 0, MaxMin should grab the far end then the middle
        let sel = farthest_point_select(&line(), 3, 0).unwrap();
        assert_eq!(sel.len(), 3);
        assert_eq!(sel[0], 0);
        assert_eq!(sel[1], 10); // farthest from 0
        assert_eq!(sel[2], 5); // maximizes min-distance to {0,10}
    }

    #[test]
    fn sphere_exclusion_picks_one_per_cluster() {
        // two tight clusters (near 0 and near 10); radius 1 excludes near-dups
        let feats = vec![vec![0.0], vec![0.1], vec![10.0], vec![10.1]];
        let scores = vec![0.9, 0.8, 0.7, 0.6];
        let sel = sphere_exclusion_select(&feats, &scores, 2, 1.0).unwrap();
        assert_eq!(sel, vec![0, 2]); // top of each cluster, not the near-duplicates
    }

    #[test]
    fn sphere_exclusion_is_prefix_stable_in_n() {
        let feats = vec![vec![0.0], vec![5.0], vec![10.0]];
        let scores = vec![0.9, 0.8, 0.7];
        let one = sphere_exclusion_select(&feats, &scores, 1, 1.0).unwrap();
        let two = sphere_exclusion_select(&feats, &scores, 2, 1.0).unwrap();
        assert_eq!(one, vec![0]);
        assert_eq!(&two[..1], &one[..]); // n=1 result is a prefix of n=2
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            farthest_point_select(&line(), 0, 0).unwrap_err().code(),
            "zero_n"
        );
        assert_eq!(
            farthest_point_select(&line(), 2, 99).unwrap_err().code(),
            "start_out_of_range"
        );
        let feats = vec![vec![0.0], vec![1.0]];
        assert_eq!(
            sphere_exclusion_select(&feats, &[0.5, 0.5], 2, 0.0)
                .unwrap_err()
                .code(),
            "non_positive_radius"
        );
    }
}
