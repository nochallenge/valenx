//! The [`Model`] abstraction every UQ routine works against.
//!
//! A UQ study never needs to know *what* a model is — only that it maps a
//! vector of real inputs to a vector of real outputs and reports its
//! dimensions. That is exactly the [`Model`] trait. Wrapping a valenx solver
//! (a CFD run, an FEM stress field, an orbit propagation, …) behind this trait
//! makes it analysable by every routine in this crate; the [`FnModel`] adapter
//! does the same for a plain closure, which is convenient in tests and for
//! quick analytic models.

/// A deterministic model `f: ℝⁿ → ℝᵐ`.
///
/// Implementations are expected to be **pure and deterministic**: evaluating
/// the same input twice yields the same output. (Sampling-based UQ on a noisy
/// model still works, but the convergence guarantees assume a deterministic
/// `f`.)
///
/// `evaluate` should return a vector of length [`Model::n_outputs`] when given
/// an input of length [`Model::n_inputs`].
pub trait Model {
    /// Evaluate the model at `inputs`, returning the output vector.
    fn evaluate(&self, inputs: &[f64]) -> Vec<f64>;

    /// The number of input dimensions the model expects.
    fn n_inputs(&self) -> usize;

    /// The number of output dimensions the model produces.
    fn n_outputs(&self) -> usize;
}

/// Adapter that turns a closure into a [`Model`].
///
/// This is the workhorse for tests and analytic models. The input/output
/// dimensions are supplied explicitly so the trait can report them without
/// evaluating the closure.
///
/// ```
/// use valenx_uq::{FnModel, Model};
///
/// // y = [x0 * x1, x0 + x1]  (2 inputs, 2 outputs)
/// let m = FnModel::new(2, 2, |x| vec![x[0] * x[1], x[0] + x[1]]);
/// assert_eq!(m.n_inputs(), 2);
/// assert_eq!(m.n_outputs(), 2);
/// assert_eq!(m.evaluate(&[3.0, 4.0]), vec![12.0, 7.0]);
/// ```
#[derive(Clone)]
pub struct FnModel<F>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    n_inputs: usize,
    n_outputs: usize,
    f: F,
}

impl<F> FnModel<F>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    /// Wrap a closure `f` declaring `n_inputs` inputs and `n_outputs` outputs.
    pub fn new(n_inputs: usize, n_outputs: usize, f: F) -> Self {
        Self {
            n_inputs,
            n_outputs,
            f,
        }
    }
}

impl<F> Model for FnModel<F>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    fn evaluate(&self, inputs: &[f64]) -> Vec<f64> {
        (self.f)(inputs)
    }

    fn n_inputs(&self) -> usize {
        self.n_inputs
    }

    fn n_outputs(&self) -> usize {
        self.n_outputs
    }
}
