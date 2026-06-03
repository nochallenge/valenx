//! Gaussian basis sets — the data model and the built-in libraries.
//!
//! A *basis set* expands each molecular orbital as a linear combination
//! of atom-centred basis functions. Each basis function is a
//! **contracted Cartesian Gaussian**: a fixed linear combination of
//! primitive Gaussians `x^l y^m z^n exp(-α r²)` sharing a centre and a
//! set of Cartesian angular exponents `(l, m, n)`.
//!
//! ## The type hierarchy
//!
//! - [`Primitive`] — one Gaussian: an exponent `α` and a contraction
//!   coefficient `c`.
//! - [`Shell`] — a group of primitives sharing a centre and an angular
//!   momentum [`AngularMomentum`]. An *s* shell holds one basis
//!   function, a *p* shell three (`x, y, z`), a *d* shell six (the
//!   Cartesian set `xx, yy, zz, xy, xz, yz`).
//! - [`BasisFunction`] — a single contracted Gaussian with definite
//!   Cartesian exponents, the unit the integrals actually loop over.
//! - [`BasisSet`] — the [`Shell`]s for a whole molecule plus the flat
//!   [`BasisFunction`] list expanded from them.
//!
//! ## Normalisation
//!
//! [`Shell::normalised`] folds the primitive normalisation constants
//! *and* the contraction normalisation into the stored coefficients, so
//! the integral code can treat every primitive coefficient as final.
//! Cartesian *d* functions are normalised on the `xx`-type component;
//! the off-diagonal `xy` components then carry an exact `√3` relative
//! factor, which the integral recursion reproduces automatically.
//!
//! ## v1 coverage
//!
//! Built-in sets: [`StoNg`](library::sto3g) (STO-3G),
//! [`Pople321`](library::pople_321g) (3-21G),
//! [`Pople631`](library::pople_631g) (6-31G) and
//! [`Pople631s`](library::pople_631gs) (6-31G*), each defined for
//! hydrogen through neon. Angular momentum runs s, p, d — enough for
//! the polarisation shells of 6-31G* on first-row atoms. f and higher
//! shells, and elements past neon, are out of v1 scope.

pub mod library;

use crate::error::{QchemError, Result};
use crate::geometry::MolecularGeometry;
use serde::{Deserialize, Serialize};

/// Angular momentum of a shell.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AngularMomentum {
    /// `l = 0` — one component, total degree 0.
    S,
    /// `l = 1` — three components `x, y, z`.
    P,
    /// `l = 2` — six Cartesian components `xx, yy, zz, xy, xz, yz`.
    D,
}

impl AngularMomentum {
    /// The total angular-momentum quantum number `l`.
    #[inline]
    pub fn l(self) -> u32 {
        match self {
            AngularMomentum::S => 0,
            AngularMomentum::P => 1,
            AngularMomentum::D => 2,
        }
    }

    /// The number of Cartesian components: `(l+1)(l+2)/2`.
    #[inline]
    pub fn n_cartesian(self) -> usize {
        let l = self.l();
        ((l + 1) * (l + 2) / 2) as usize
    }

    /// The Cartesian exponent triples `(i, j, k)` for this shell, in a
    /// fixed canonical order. For *p* this is `x, y, z`; for *d* it is
    /// `xx, yy, zz, xy, xz, yz`.
    pub fn cartesian_components(self) -> Vec<(u32, u32, u32)> {
        match self {
            AngularMomentum::S => vec![(0, 0, 0)],
            AngularMomentum::P => vec![(1, 0, 0), (0, 1, 0), (0, 0, 1)],
            AngularMomentum::D => vec![
                (2, 0, 0),
                (0, 2, 0),
                (0, 0, 2),
                (1, 1, 0),
                (1, 0, 1),
                (0, 1, 1),
            ],
        }
    }
}

/// A single primitive Gaussian inside a contraction.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Primitive {
    /// Gaussian exponent `α` (always `> 0`).
    pub exponent: f64,
    /// Contraction coefficient `c`. As stored in a [`Shell`] obtained
    /// from [`Shell::normalised`] this already folds in the primitive
    /// and contraction normalisation factors.
    pub coefficient: f64,
}

/// A contracted Gaussian shell centred on one atom.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Shell {
    /// Index of the atom this shell sits on (into
    /// [`MolecularGeometry::atoms`]).
    pub atom_index: usize,
    /// Shell centre in bohr — a copy of the atom's position so the
    /// integral code does not need the geometry alongside.
    pub centre: [f64; 3],
    /// The shell's angular momentum.
    pub angular: AngularMomentum,
    /// The contracted primitives.
    pub primitives: Vec<Primitive>,
}

impl Shell {
    /// Number of basis functions this shell contributes — the number of
    /// Cartesian components of its angular momentum.
    #[inline]
    pub fn n_functions(&self) -> usize {
        self.angular.n_cartesian()
    }

    /// Return a copy of this shell with the primitive coefficients
    /// rescaled so that every basis function in it is normalised
    /// (unit self-overlap on the `xx`-type Cartesian component).
    ///
    /// This folds two factors into the coefficients — the per-primitive
    /// Cartesian-Gaussian normalisation, and the contraction
    /// self-overlap normalisation computed from the closed-form overlap
    /// of two primitives — so the integral routines can treat the
    /// result as final.
    pub fn normalised(&self) -> Shell {
        let l = self.angular.l();
        // Per-primitive normalisation for the leading Cartesian
        // component (l,0,0): N = (2α/π)^{3/4} (4α)^{l/2} / sqrt((2l-1)!!).
        let dfact = double_factorial_2lm1(l);
        let mut prims: Vec<Primitive> = self
            .primitives
            .iter()
            .map(|p| {
                let a = p.exponent;
                let n = (2.0 * a / std::f64::consts::PI).powf(0.75)
                    * (4.0 * a).powf(l as f64 / 2.0)
                    / dfact.sqrt();
                Primitive {
                    exponent: a,
                    coefficient: p.coefficient * n,
                }
            })
            .collect();

        // Contraction self-overlap on the (l,0,0) component:
        // S = Σ_ij c_i c_j (π/(α_i+α_j))^{3/2} (2l-1)!! / (2(α_i+α_j))^l
        let mut s = 0.0;
        for pi in &prims {
            for pj in &prims {
                let g = pi.exponent + pj.exponent;
                s += pi.coefficient * pj.coefficient
                    * (std::f64::consts::PI / g).powf(1.5)
                    * dfact
                    / (2.0 * g).powi(l as i32);
            }
        }
        let scale = 1.0 / s.sqrt();
        for p in &mut prims {
            p.coefficient *= scale;
        }
        Shell {
            atom_index: self.atom_index,
            centre: self.centre,
            angular: self.angular,
            primitives: prims,
        }
    }
}

/// A single contracted Cartesian Gaussian basis function — one shell
/// component with definite Cartesian exponents.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BasisFunction {
    /// Index of the owning atom.
    pub atom_index: usize,
    /// Basis-function centre in bohr.
    pub centre: [f64; 3],
    /// Cartesian angular exponents `(i, j, k)` — `x^i y^j z^k`.
    pub cart: (u32, u32, u32),
    /// The (normalised) contracted primitives.
    pub primitives: Vec<Primitive>,
}

impl BasisFunction {
    /// Total angular momentum `l = i + j + k`.
    #[inline]
    pub fn l(&self) -> u32 {
        self.cart.0 + self.cart.1 + self.cart.2
    }
}

/// The basis set of a complete molecule.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BasisSet {
    /// Human-readable basis-set name (`"sto-3g"`, …).
    pub name: &'static str,
    /// All shells, grouped by atom in input order.
    pub shells: Vec<Shell>,
    /// The flat basis-function list expanded from the shells. The SCF
    /// matrices are indexed by position in this vector.
    pub functions: Vec<BasisFunction>,
}

impl BasisSet {
    /// The number of basis functions — the dimension of every SCF
    /// matrix.
    #[inline]
    pub fn n_functions(&self) -> usize {
        self.functions.len()
    }

    /// Build the molecular basis set by looking each atom's element up
    /// in a built-in library and expanding the shells into basis
    /// functions. The shells are normalised on the way in.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::BasisNotFound`] when the named set has no
    /// definition for an element in the molecule, and
    /// [`QchemError::Parse`] when `name` is not a known basis-set name.
    pub fn build(name: &str, geometry: &MolecularGeometry) -> Result<BasisSet> {
        let lib = library::resolve(name)?;
        let static_name = lib.name();
        let mut shells = Vec::new();
        for (atom_index, atom) in geometry.atoms.iter().enumerate() {
            let z = atom.element.atomic_number();
            let element_shells = lib.shells_for(z).ok_or_else(|| {
                QchemError::basis_not_found(static_name, atom.element.symbol())
            })?;
            for raw in element_shells {
                let shell = Shell {
                    atom_index,
                    centre: atom.position,
                    angular: raw.angular,
                    primitives: raw.primitives.clone(),
                }
                .normalised();
                shells.push(shell);
            }
        }
        let functions = expand_functions(&shells);
        Ok(BasisSet {
            name: static_name,
            shells,
            functions,
        })
    }
}

/// Expand every shell into its Cartesian basis functions.
fn expand_functions(shells: &[Shell]) -> Vec<BasisFunction> {
    let mut out = Vec::new();
    for shell in shells {
        for cart in shell.angular.cartesian_components() {
            out.push(BasisFunction {
                atom_index: shell.atom_index,
                centre: shell.centre,
                cart,
                primitives: shell.primitives.clone(),
            });
        }
    }
    out
}

/// `(2l - 1)!!` — the double factorial that appears in the
/// Cartesian-Gaussian normalisation of the leading component. `(2·0-1)!!`
/// and `(2·1-1)!!` are both 1; `(2·2-1)!! = 3`.
pub(crate) fn double_factorial_2lm1(l: u32) -> f64 {
    match l {
        0 | 1 => 1.0,
        2 => 3.0,
        3 => 15.0,
        _ => {
            // General (2l-1)!! for completeness.
            let mut v = 1.0;
            let mut k = 2 * l as i64 - 1;
            while k > 1 {
                v *= k as f64;
                k -= 2;
            }
            v
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;

    #[test]
    fn angular_component_counts() {
        assert_eq!(AngularMomentum::S.n_cartesian(), 1);
        assert_eq!(AngularMomentum::P.n_cartesian(), 3);
        assert_eq!(AngularMomentum::D.n_cartesian(), 6);
        assert_eq!(AngularMomentum::D.cartesian_components().len(), 6);
    }

    #[test]
    fn normalised_s_shell_has_unit_self_overlap() {
        // A single uncontracted s primitive, normalised, must overlap
        // itself to 1: S = c² (π/2α)^{3/2}.
        let shell = Shell {
            atom_index: 0,
            centre: [0.0; 3],
            angular: AngularMomentum::S,
            primitives: vec![Primitive {
                exponent: 1.3,
                coefficient: 1.0,
            }],
        }
        .normalised();
        let c = shell.primitives[0].coefficient;
        let a = shell.primitives[0].exponent;
        let s = c * c * (std::f64::consts::PI / (2.0 * a)).powf(1.5);
        assert!((s - 1.0).abs() < 1.0e-12, "self-overlap {s}");
    }

    #[test]
    fn sto3g_water_has_seven_functions() {
        let atoms = vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.757, 0.587]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.757, 0.587]).unwrap(),
        ];
        let geom = MolecularGeometry::new(atoms);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        // O: 1s, 2s, 2p (5) ; each H: 1s (1) → 5 + 1 + 1 = 7.
        assert_eq!(basis.n_functions(), 7);
    }

    #[test]
    fn unknown_basis_name_errors() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
        ]);
        assert!(BasisSet::build("cc-pvqz", &geom).is_err());
    }

    #[test]
    fn double_factorial_values() {
        assert_eq!(double_factorial_2lm1(0), 1.0);
        assert_eq!(double_factorial_2lm1(1), 1.0);
        assert_eq!(double_factorial_2lm1(2), 3.0);
    }
}
