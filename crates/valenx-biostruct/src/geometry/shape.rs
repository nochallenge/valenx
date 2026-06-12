//! Bulk shape descriptors: centre of mass, radius of gyration and
//! the principal axes of inertia.
//!
//! All functions operate on a slice of `(coordinate, mass)` pairs so
//! callers can choose all-atom, Cα-only or unit-mass weighting.

use crate::error::{BiostructError, Result};
use crate::structure::Model;
use nalgebra::{Matrix3, Point3, Vector3};

/// Approximate atomic mass of an element symbol, in daltons. Unknown
/// elements return the carbon mass.
pub fn atomic_mass(element: &str) -> f64 {
    match element.to_ascii_uppercase().as_str() {
        "H" | "D" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "F" => 18.998,
        "NA" => 22.990,
        "MG" => 24.305,
        "P" => 30.974,
        "S" => 32.06,
        "CL" => 35.45,
        "K" => 39.098,
        "CA" => 40.078,
        "MN" => 54.938,
        "FE" => 55.845,
        "CO" => 58.933,
        "NI" => 58.693,
        "CU" => 63.546,
        "ZN" => 65.38,
        "SE" => 78.971,
        "BR" => 79.904,
        "I" => 126.904,
        _ => 12.011,
    }
}

/// Mass-weighted centre of mass of a point cloud.
pub fn center_of_mass(points: &[(Point3<f64>, f64)]) -> Result<Point3<f64>> {
    if points.is_empty() {
        return Err(BiostructError::invalid("points", "empty point set"));
    }
    let total: f64 = points.iter().map(|(_, m)| m).sum();
    if total <= 0.0 {
        return Err(BiostructError::invalid("masses", "total mass is zero"));
    }
    let mut acc = Vector3::zeros();
    for (p, m) in points {
        acc += p.coords * *m;
    }
    Ok(Point3::from(acc / total))
}

/// Unweighted geometric centroid (centre of *coordinates*).
pub fn centroid(points: &[Point3<f64>]) -> Result<Point3<f64>> {
    if points.is_empty() {
        return Err(BiostructError::invalid("points", "empty point set"));
    }
    let mut acc = Vector3::zeros();
    for p in points {
        acc += p.coords;
    }
    Ok(Point3::from(acc / points.len() as f64))
}

/// Mass-weighted radius of gyration, ångström.
///
/// `Rg = sqrt( Σ mᵢ |rᵢ − r_com|² / Σ mᵢ )`.
pub fn radius_of_gyration(points: &[(Point3<f64>, f64)]) -> Result<f64> {
    let com = center_of_mass(points)?;
    let total: f64 = points.iter().map(|(_, m)| m).sum();
    let mut sum = 0.0;
    for (p, m) in points {
        sum += m * (p - com).norm_squared();
    }
    Ok((sum / total).sqrt())
}

/// The principal axes of inertia of a mass-weighted point cloud.
#[derive(Clone, Debug, PartialEq)]
pub struct PrincipalAxes {
    /// Centre of mass the inertia tensor was computed about.
    pub center: Point3<f64>,
    /// Three orthonormal principal-axis directions, ordered by
    /// ascending principal moment.
    pub axes: [Vector3<f64>; 3],
    /// The three principal moments of inertia, ascending.
    pub moments: [f64; 3],
}

impl PrincipalAxes {
    /// Asphericity-style anisotropy: ratio of the largest to the
    /// smallest principal moment. `1.0` is perfectly spherical; large
    /// values mean an elongated / flat shape.
    pub fn anisotropy(&self) -> f64 {
        if self.moments[0].abs() < 1e-12 {
            return f64::INFINITY;
        }
        self.moments[2] / self.moments[0]
    }
}

/// Diagonalise the inertia tensor of a point cloud to obtain its
/// principal axes.
///
/// The 3×3 inertia tensor is symmetric, so `nalgebra`'s symmetric
/// eigendecomposition gives real eigenvalues (the principal moments)
/// and orthonormal eigenvectors (the principal axes).
pub fn principal_axes(points: &[(Point3<f64>, f64)]) -> Result<PrincipalAxes> {
    let com = center_of_mass(points)?;
    let mut tensor = Matrix3::zeros();
    for (p, m) in points {
        let r = p - com;
        let (x, y, z) = (r.x, r.y, r.z);
        // Standard inertia tensor: diagonal Σm(r²−c²), off-diag −Σm·c·c′.
        tensor[(0, 0)] += m * (y * y + z * z);
        tensor[(1, 1)] += m * (x * x + z * z);
        tensor[(2, 2)] += m * (x * x + y * y);
        tensor[(0, 1)] -= m * x * y;
        tensor[(0, 2)] -= m * x * z;
        tensor[(1, 2)] -= m * y * z;
    }
    tensor[(1, 0)] = tensor[(0, 1)];
    tensor[(2, 0)] = tensor[(0, 2)];
    tensor[(2, 1)] = tensor[(1, 2)];

    let eig = nalgebra::SymmetricEigen::new(tensor);
    // Sort eigenpairs by ascending eigenvalue.
    let mut pairs: Vec<(f64, Vector3<f64>)> = (0..3)
        .map(|i| (eig.eigenvalues[i], eig.eigenvectors.column(i).into_owned()))
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    Ok(PrincipalAxes {
        center: com,
        axes: [pairs[0].1, pairs[1].1, pairs[2].1],
        moments: [pairs[0].0, pairs[1].0, pairs[2].0],
    })
}

/// Convenience: build the `(coord, mass)` list of every atom of a
/// model, with real element masses.
pub fn model_mass_points(model: &Model) -> Vec<(Point3<f64>, f64)> {
    model
        .atoms()
        .map(|a| (a.coord, atomic_mass(&a.element)))
        .collect()
}

/// Convenience: build the `(coord, mass)` list using unit masses
/// (every atom weighted `1.0`) — the geometric, not physical, variant.
pub fn model_unit_points(model: &Model) -> Vec<(Point3<f64>, f64)> {
    model.atoms().map(|a| (a.coord, 1.0)).collect()
}

/// Radius of gyration of a whole model using real atomic masses.
pub fn model_radius_of_gyration(model: &Model) -> Result<f64> {
    radius_of_gyration(&model_mass_points(model))
}

/// Axis-aligned bounding box `(min, max)` of a model's atoms.
pub fn bounding_box(model: &Model) -> Result<(Point3<f64>, Point3<f64>)> {
    let mut it = model.atoms();
    let first = it
        .next()
        .ok_or_else(|| BiostructError::invalid("model", "model has no atoms"))?;
    let mut lo = first.coord;
    let mut hi = first.coord;
    for a in model.atoms() {
        lo.x = lo.x.min(a.coord.x);
        lo.y = lo.y.min(a.coord.y);
        lo.z = lo.z.min(a.coord.z);
        hi.x = hi.x.max(a.coord.x);
        hi.y = hi.y.max(a.coord.y);
        hi.z = hi.z.max(a.coord.z);
    }
    Ok((lo, hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn com_of_symmetric_pair() {
        let pts = vec![
            (Point3::new(-1.0, 0.0, 0.0), 1.0),
            (Point3::new(1.0, 0.0, 0.0), 1.0),
        ];
        let com = center_of_mass(&pts).unwrap();
        assert!((com - Point3::origin()).norm() < 1e-12);
    }

    #[test]
    fn com_is_mass_weighted() {
        let pts = vec![
            (Point3::new(0.0, 0.0, 0.0), 1.0),
            (Point3::new(10.0, 0.0, 0.0), 3.0),
        ];
        let com = center_of_mass(&pts).unwrap();
        // weighted mean = (0*1 + 10*3)/4 = 7.5
        assert!((com.x - 7.5).abs() < 1e-12);
    }

    #[test]
    fn rg_of_a_known_shell() {
        // Six unit masses at +/-1 on each axis: every point is 1 A
        // from the centroid, so Rg = 1.
        let pts = vec![
            (Point3::new(1.0, 0.0, 0.0), 1.0),
            (Point3::new(-1.0, 0.0, 0.0), 1.0),
            (Point3::new(0.0, 1.0, 0.0), 1.0),
            (Point3::new(0.0, -1.0, 0.0), 1.0),
            (Point3::new(0.0, 0.0, 1.0), 1.0),
            (Point3::new(0.0, 0.0, -1.0), 1.0),
        ];
        assert!((radius_of_gyration(&pts).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rg_of_a_cube_matches_closed_form() {
        // GROUND TRUTH: the eight corners of a cube of side `a` centred
        // at the origin sit at (±a/2, ±a/2, ±a/2). By symmetry the centre
        // of mass is the origin, and every corner is at distance
        //   |r| = √((a/2)² · 3) = (a/2)·√3
        // from it. With equal masses the mass weighting drops out and
        //   Rg = sqrt( (1/8) · Σ |rᵢ|² ) = (a/2)·√3.
        // For a = 2 this is exactly √3 = 1.7320508075688772.
        //
        // Distinct from `rg_of_a_known_shell` (octahedral, Rg=1): here the
        // corners are NOT all the same trivial unit distance, so this
        // exercises the |r|² accumulation with a non-unit per-point term.
        // Tolerance 1e-12 — the identity is algebraically exact.
        let a = 2.0_f64;
        let h = a / 2.0;
        let corners = [-h, h];
        let mut pts = Vec::with_capacity(8);
        for &x in &corners {
            for &y in &corners {
                for &z in &corners {
                    pts.push((Point3::new(x, y, z), 1.0));
                }
            }
        }
        let rg = radius_of_gyration(&pts).unwrap();
        let expected = h * 3.0_f64.sqrt(); // (a/2)·√3
        assert!(
            (rg - expected).abs() < 1e-12,
            "cube Rg = {rg} ≠ (a/2)·√3 = {expected}"
        );
    }

    #[test]
    fn principal_axes_of_a_rod() {
        // A rod along x: the smallest moment is about x, the two
        // larger moments are degenerate about y and z.
        let pts: Vec<_> = (-5..=5)
            .map(|i| (Point3::new(i as f64, 0.0, 0.0), 1.0))
            .collect();
        let pa = principal_axes(&pts).unwrap();
        // moment about the rod axis is ~0, the others are large.
        assert!(pa.moments[0] < 1e-6, "rod axis moment {}", pa.moments[0]);
        assert!(pa.moments[2] > 100.0);
        // the smallest-moment axis is (anti)parallel to x.
        assert!(pa.axes[0].x.abs() > 0.99);
    }

    #[test]
    fn empty_set_errors() {
        assert!(center_of_mass(&[]).is_err());
        assert!(centroid(&[]).is_err());
    }

    #[test]
    fn atomic_masses() {
        assert!((atomic_mass("C") - 12.011).abs() < 1e-9);
        assert!((atomic_mass("zz") - 12.011).abs() < 1e-9); // fallback
    }
}
