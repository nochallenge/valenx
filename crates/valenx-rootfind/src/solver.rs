//! The three classical scalar root-finding algorithms.
//!
//! Each solver searches for an `x` with `f(x) = 0` for a user-supplied
//! `f: Fn(f64) -> f64`, and shares a common
//! [`Settings`] (tolerance + iteration cap) and
//! a common [`Root`] result (root, residual,
//! iterations, method). The methods differ in what extra information
//! they require and how fast they converge:
//!
//! | Method | Extra input | Guarantee | Order |
//! | --- | --- | --- | --- |
//! | [`bisection`] | a sign-change bracket `[a, b]` | always converges | linear (halves each step) |
//! | [`newton`] | the derivative `f'` | needs a good start | quadratic |
//! | [`secant`] | two seed points | needs good starts | super-linear (≈1.618) |
//!
//! ## Convergence criterion
//!
//! All three stop as soon as the residual `|f(x)|` drops below the
//! tolerance. The bracketing method additionally stops once the bracket
//! width falls below the tolerance (the root is then pinned to within
//! `tol`). The open methods additionally stop when the step
//! `|x_{n+1} - x_n|` falls below the tolerance. If none of these is met
//! within the iteration budget the solver returns
//! [`RootError::MaxIterations`].
//!
//! ## Honest scope
//!
//! These are the textbook one-dimensional methods exactly as they
//! appear in an introductory numerical-analysis course. They make no
//! attempt at the production-grade robustness of a hybrid solver such as
//! Brent's method (which would safeguard the open steps with a bracket):
//! Newton and the secant can diverge or oscillate from a poor start, and
//! a vanishing slope is reported as an error rather than recovered from.
//! Research / educational grade only.

use crate::error::{Result, RootError};
use serde::{Deserialize, Serialize};

/// Default convergence tolerance used by [`Settings::default`].
///
/// Comfortably tight for `f64` work while leaving headroom above the
/// machine epsilon (`≈2.22e-16`).
pub const DEFAULT_TOL: f64 = 1e-12;

/// Default iteration cap used by [`Settings::default`].
///
/// Bisection needs about `log2((b - a) / tol)` steps (well under 100 for
/// any reasonable bracket); the open methods converge in a handful when
/// they converge at all, so 100 is a generous ceiling.
pub const DEFAULT_MAX_ITER: usize = 100;

/// Magnitude below which a slope (Newton's `f'`, or the secant's finite
/// difference) is treated as zero and the step is refused.
///
/// Dividing by a slope smaller than this would produce a wild or
/// non-finite update, so the solver raises
/// [`RootError::ZeroDerivative`] instead.
pub const SLOPE_GUARD: f64 = 1e-300;

/// Which method produced a [`Root`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Method {
    /// The bisection (interval-halving) bracketing method.
    Bisection,
    /// The Newton-Raphson tangent-line method.
    Newton,
    /// The secant (chord) method.
    Secant,
}

impl Method {
    /// Stable lowercase identifier, e.g. `"bisection"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            Method::Bisection => "bisection",
            Method::Newton => "newton",
            Method::Secant => "secant",
        }
    }
}

/// Shared stopping criteria for every solver: a convergence tolerance
/// and an iteration cap.
///
/// Construct with [`Settings::new`] (validating) or via
/// [`Settings::default`] for the [`DEFAULT_TOL`] / [`DEFAULT_MAX_ITER`]
/// pair.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Settings {
    /// Convergence tolerance on `|f(x)|` (and on the bracket width /
    /// step). Must be strictly positive and finite.
    pub tol: f64,
    /// Maximum number of iterations before the solver gives up with
    /// [`RootError::MaxIterations`]. Must be at least one.
    pub max_iter: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            tol: DEFAULT_TOL,
            max_iter: DEFAULT_MAX_ITER,
        }
    }
}

impl Settings {
    /// Build validated [`Settings`].
    ///
    /// # Errors
    ///
    /// Returns [`RootError::Invalid`] if `tol` is not strictly positive
    /// and finite, or if `max_iter` is zero.
    pub fn new(tol: f64, max_iter: usize) -> Result<Self> {
        if !(tol.is_finite() && tol > 0.0) {
            return Err(RootError::invalid(
                "tol",
                format!("must be a finite, strictly positive number, got {tol}"),
            ));
        }
        if max_iter == 0 {
            return Err(RootError::invalid("max_iter", "must be at least 1"));
        }
        Ok(Settings { tol, max_iter })
    }
}

/// A successfully located root together with diagnostics.
///
/// Returned by every solver on success. `residual` is `|f(root)|` at the
/// reported root and is guaranteed `<= settings.tol` *unless* the method
/// converged purely on a small bracket width or step (in which case the
/// residual may still be slightly above the tolerance while the root
/// itself is pinned to within `tol`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Root {
    /// The located root `x` with `f(x) ≈ 0`.
    pub root: f64,
    /// The residual `|f(root)|` at the reported root.
    pub residual: f64,
    /// Number of iterations the solver performed to reach `root`.
    pub iterations: usize,
    /// Which method produced this result.
    pub method: Method,
}

/// Evaluate `f` at `x`, returning [`RootError::NotFinite`] if the result
/// is not a finite number.
#[inline]
fn eval_finite<F>(f: &F, x: f64, what: &'static str) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    let y = f(x);
    if y.is_finite() {
        Ok(y)
    } else {
        Err(RootError::NotFinite { what, x })
    }
}

/// Find a root of `f` on the bracket `[a, b]` by bisection.
///
/// Requires a genuine sign-change bracket: `f(a)` and `f(b)` must have
/// opposite signs (one may be exactly zero, in which case that endpoint
/// is returned immediately). The interval is repeatedly halved, always
/// keeping the half that still brackets a sign change, until either the
/// midpoint residual or the bracket width falls below `settings.tol`.
/// Convergence is guaranteed and linear — the bracket width is halved
/// every iteration — so the iterate count is bounded by
/// `log2((b - a) / tol)`.
///
/// The endpoints may be given in either order (`a > b` is accepted and
/// swapped internally).
///
/// # Errors
///
/// - [`RootError::NoSignChange`] if `f(a)` and `f(b)` share a sign and
///   neither is zero.
/// - [`RootError::NotFinite`] if `f` returns a non-finite value at an
///   endpoint or midpoint.
/// - [`RootError::MaxIterations`] if the tolerance is not met within
///   `settings.max_iter` halvings.
///
/// # Examples
///
/// ```
/// use valenx_rootfind::{bisection, Settings};
///
/// // sqrt(2) is the positive root of x^2 - 2.
/// let r = bisection(|x| x * x - 2.0, 0.0, 2.0, Settings::default()).unwrap();
/// assert!((r.root - 2.0_f64.sqrt()).abs() < 1e-10);
/// ```
pub fn bisection<F>(f: F, a: f64, b: f64, settings: Settings) -> Result<Root>
where
    F: Fn(f64) -> f64,
{
    if !a.is_finite() {
        return Err(RootError::invalid(
            "a",
            format!("bracket endpoint {a} is not finite"),
        ));
    }
    if !b.is_finite() {
        return Err(RootError::invalid(
            "b",
            format!("bracket endpoint {b} is not finite"),
        ));
    }

    // Accept either order.
    let (mut lo, mut hi) = if a <= b { (a, b) } else { (b, a) };

    let mut flo = eval_finite(&f, lo, "f(x)")?;
    // `fhi` is read only for the endpoint and sign-change checks below;
    // once iterating, the live bracket is tracked by `lo` / `flo` and the
    // fresh `fmid`, so the upper value never needs updating.
    let fhi = eval_finite(&f, hi, "f(x)")?;

    // An endpoint already on the root is the trivial answer.
    if flo == 0.0 {
        return Ok(Root {
            root: lo,
            residual: 0.0,
            iterations: 0,
            method: Method::Bisection,
        });
    }
    if fhi == 0.0 {
        return Ok(Root {
            root: hi,
            residual: 0.0,
            iterations: 0,
            method: Method::Bisection,
        });
    }

    // Bracketing requires a genuine sign change.
    if flo.signum() == fhi.signum() {
        return Err(RootError::NoSignChange {
            a: lo,
            b: hi,
            fa: flo,
            fb: fhi,
        });
    }

    let mut mid = 0.5 * (lo + hi);
    let mut fmid = 0.0_f64;

    for iter in 1..=settings.max_iter {
        mid = 0.5 * (lo + hi);
        fmid = eval_finite(&f, mid, "f(x)")?;

        // Converged: residual small, or the bracket has shrunk below
        // the tolerance so the root is pinned to within `tol`.
        if fmid.abs() <= settings.tol || 0.5 * (hi - lo) <= settings.tol {
            return Ok(Root {
                root: mid,
                residual: fmid.abs(),
                iterations: iter,
                method: Method::Bisection,
            });
        }

        // Keep the half that still straddles the sign change.
        if flo.signum() == fmid.signum() {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
        }
    }

    Err(RootError::MaxIterations {
        max_iter: settings.max_iter,
        tol: settings.tol,
        last: mid,
        residual: fmid.abs(),
    })
}

/// Find a root of `f` by the Newton-Raphson method, starting from `x0`
/// and using the supplied derivative `df = f'`.
///
/// Each step takes the tangent line at the current iterate and jumps to
/// where it crosses zero: `x_{n+1} = x_n - f(x_n) / f'(x_n)`.
/// Convergence is quadratic near a simple root, so a handful of
/// iterations typically suffice — but, unlike bisection, there is no
/// safeguard: a poor start or a flat region can send the iterate far
/// away.
///
/// # Errors
///
/// - [`RootError::ZeroDerivative`] if `|f'(x_n)|` falls below
///   [`SLOPE_GUARD`] (the tangent is too flat to step along).
/// - [`RootError::NotFinite`] if `f` or `df` returns a non-finite value.
/// - [`RootError::MaxIterations`] if the tolerance is not met within
///   `settings.max_iter` steps.
///
/// # Examples
///
/// ```
/// use valenx_rootfind::{newton, Settings};
///
/// // sqrt(2): f = x^2 - 2, f' = 2x.
/// let r = newton(|x| x * x - 2.0, |x| 2.0 * x, 1.0, Settings::default()).unwrap();
/// assert!((r.root - 2.0_f64.sqrt()).abs() < 1e-10);
/// ```
pub fn newton<F, DF>(f: F, df: DF, x0: f64, settings: Settings) -> Result<Root>
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    if !x0.is_finite() {
        return Err(RootError::invalid(
            "x0",
            format!("start point {x0} is not finite"),
        ));
    }

    let mut x = x0;
    let mut fx = eval_finite(&f, x, "f(x)")?;

    for iter in 1..=settings.max_iter {
        // Already a root.
        if fx.abs() <= settings.tol {
            return Ok(Root {
                root: x,
                residual: fx.abs(),
                iterations: iter - 1,
                method: Method::Newton,
            });
        }

        let dfx = eval_finite(&df, x, "f'(x)")?;
        if dfx.abs() < SLOPE_GUARD {
            return Err(RootError::ZeroDerivative {
                x,
                slope: dfx,
                threshold: SLOPE_GUARD,
            });
        }

        let x_next = x - fx / dfx;
        if !x_next.is_finite() {
            return Err(RootError::NotFinite {
                what: "iterate",
                x: x_next,
            });
        }

        let step = (x_next - x).abs();
        let f_next = eval_finite(&f, x_next, "f(x)")?;

        // Converged on a small residual or a small step.
        if f_next.abs() <= settings.tol || step <= settings.tol {
            return Ok(Root {
                root: x_next,
                residual: f_next.abs(),
                iterations: iter,
                method: Method::Newton,
            });
        }

        x = x_next;
        fx = f_next;
    }

    Err(RootError::MaxIterations {
        max_iter: settings.max_iter,
        tol: settings.tol,
        last: x,
        residual: fx.abs(),
    })
}

/// Find a root of `f` by the secant method, starting from two seed
/// points `x0` and `x1`.
///
/// The secant method replaces Newton's analytic derivative with the
/// finite-difference slope through the two most recent iterates, so it
/// needs **no derivative** — only the function itself:
/// `x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))`.
/// Convergence is super-linear (order ≈ 1.618, the golden ratio) — a
/// touch slower than Newton per step but cheaper, since it makes a
/// single function call per iteration and asks nothing of the caller
/// beyond `f`.
///
/// # Errors
///
/// - [`RootError::Invalid`] if `x0 == x1` (the first chord would have
///   zero width).
/// - [`RootError::ZeroDerivative`] if `|f(x_n) - f(x_{n-1})|` falls below
///   [`SLOPE_GUARD`] (the chord is horizontal).
/// - [`RootError::NotFinite`] if `f` returns a non-finite value.
/// - [`RootError::MaxIterations`] if the tolerance is not met within
///   `settings.max_iter` steps.
///
/// # Examples
///
/// ```
/// use valenx_rootfind::{secant, Settings};
///
/// // sqrt(2) with no derivative supplied.
/// let r = secant(|x| x * x - 2.0, 1.0, 2.0, Settings::default()).unwrap();
/// assert!((r.root - 2.0_f64.sqrt()).abs() < 1e-10);
/// ```
pub fn secant<F>(f: F, x0: f64, x1: f64, settings: Settings) -> Result<Root>
where
    F: Fn(f64) -> f64,
{
    if !x0.is_finite() {
        return Err(RootError::invalid("x0", format!("seed {x0} is not finite")));
    }
    if !x1.is_finite() {
        return Err(RootError::invalid("x1", format!("seed {x1} is not finite")));
    }
    if x0 == x1 {
        return Err(RootError::invalid(
            "x0/x1",
            "the two secant seeds must be distinct",
        ));
    }

    let mut x_prev = x0;
    let mut x_curr = x1;
    let mut f_prev = eval_finite(&f, x_prev, "f(x)")?;
    let mut f_curr = eval_finite(&f, x_curr, "f(x)")?;

    // If a seed is already a root, report it.
    if f_prev.abs() <= settings.tol {
        return Ok(Root {
            root: x_prev,
            residual: f_prev.abs(),
            iterations: 0,
            method: Method::Secant,
        });
    }

    for iter in 1..=settings.max_iter {
        if f_curr.abs() <= settings.tol {
            return Ok(Root {
                root: x_curr,
                residual: f_curr.abs(),
                iterations: iter - 1,
                method: Method::Secant,
            });
        }

        let denom = f_curr - f_prev;
        if denom.abs() < SLOPE_GUARD {
            return Err(RootError::ZeroDerivative {
                x: x_curr,
                slope: denom,
                threshold: SLOPE_GUARD,
            });
        }

        let x_next = x_curr - f_curr * (x_curr - x_prev) / denom;
        if !x_next.is_finite() {
            return Err(RootError::NotFinite {
                what: "iterate",
                x: x_next,
            });
        }

        let step = (x_next - x_curr).abs();
        let f_next = eval_finite(&f, x_next, "f(x)")?;

        if f_next.abs() <= settings.tol || step <= settings.tol {
            return Ok(Root {
                root: x_next,
                residual: f_next.abs(),
                iterations: iter,
                method: Method::Secant,
            });
        }

        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_next;
        f_curr = f_next;
    }

    Err(RootError::MaxIterations {
        max_iter: settings.max_iter,
        tol: settings.tol,
        last: x_curr,
        residual: f_curr.abs(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `f(x) = x^2 - 2`, whose positive root is sqrt(2).
    fn quad(x: f64) -> f64 {
        x * x - 2.0
    }
    /// Its derivative `f'(x) = 2x`.
    fn quad_prime(x: f64) -> f64 {
        2.0 * x
    }

    const SQRT2: f64 = std::f64::consts::SQRT_2;
    const EPS: f64 = 1e-10;

    // --- settings validation ------------------------------------------

    #[test]
    fn settings_reject_bad_tolerance() {
        assert!(Settings::new(0.0, 50).is_err());
        assert!(Settings::new(-1e-6, 50).is_err());
        assert!(Settings::new(f64::NAN, 50).is_err());
        assert!(Settings::new(f64::INFINITY, 50).is_err());
    }

    #[test]
    fn settings_reject_zero_budget() {
        assert!(Settings::new(1e-9, 0).is_err());
    }

    #[test]
    fn settings_accept_good_values() {
        let s = Settings::new(1e-9, 25).unwrap();
        assert_eq!(s.max_iter, 25);
        assert!((s.tol - 1e-9).abs() < 1e-18);
    }

    #[test]
    fn default_settings_are_sane() {
        let s = Settings::default();
        assert_eq!(s.tol, DEFAULT_TOL);
        assert_eq!(s.max_iter, DEFAULT_MAX_ITER);
    }

    // --- all three find sqrt(2) ---------------------------------------

    #[test]
    fn bisection_finds_sqrt2() {
        let r = bisection(quad, 0.0, 2.0, Settings::default()).unwrap();
        assert!((r.root - SQRT2).abs() < EPS, "root = {}", r.root);
        assert_eq!(r.method, Method::Bisection);
    }

    #[test]
    fn newton_finds_sqrt2() {
        let r = newton(quad, quad_prime, 1.0, Settings::default()).unwrap();
        assert!((r.root - SQRT2).abs() < EPS, "root = {}", r.root);
        assert_eq!(r.method, Method::Newton);
    }

    #[test]
    fn secant_finds_sqrt2() {
        let r = secant(quad, 1.0, 2.0, Settings::default()).unwrap();
        assert!((r.root - SQRT2).abs() < EPS, "root = {}", r.root);
        assert_eq!(r.method, Method::Secant);
    }

    // --- bracketing requirement ---------------------------------------

    #[test]
    fn bisection_requires_sign_change() {
        // f > 0 on the whole interval [2, 3]: no bracket.
        let err = bisection(quad, 2.0, 3.0, Settings::default()).unwrap_err();
        assert!(matches!(err, RootError::NoSignChange { .. }));
        assert_eq!(err.code(), "rootfind.no_sign_change");
    }

    #[test]
    fn bisection_accepts_reversed_bracket() {
        // Endpoints in descending order still work.
        let r = bisection(quad, 2.0, 0.0, Settings::default()).unwrap();
        assert!((r.root - SQRT2).abs() < EPS, "root = {}", r.root);
    }

    #[test]
    fn bisection_returns_endpoint_on_exact_root() {
        // For a linear f, x = 2 is an *exact* (bit-for-bit zero) root, so
        // when it is handed in as an endpoint the solver short-circuits
        // before iterating rather than reporting a missing sign change.
        let line = |x: f64| x - 2.0;
        assert_eq!(line(2.0), 0.0, "premise: endpoint is an exact root");
        let r = bisection(line, 2.0, 5.0, Settings::default()).unwrap();
        assert!((r.root - 2.0).abs() < EPS, "root = {}", r.root);
        assert!(r.residual.abs() < EPS);
        assert_eq!(r.iterations, 0);
    }

    // --- convergence speed: Newton is fast ----------------------------

    #[test]
    fn newton_converges_in_few_iterations() {
        // From x0 = 1, quadratic convergence reaches 1e-12 well inside
        // ten iterations — far fewer than bisection needs.
        let r = newton(quad, quad_prime, 1.0, Settings::default()).unwrap();
        assert!(r.iterations <= 8, "took {} iterations", r.iterations);
    }

    #[test]
    fn newton_beats_bisection_on_iteration_count() {
        let n = newton(quad, quad_prime, 1.0, Settings::default()).unwrap();
        let b = bisection(quad, 0.0, 2.0, Settings::default()).unwrap();
        assert!(
            n.iterations < b.iterations,
            "newton {} vs bisection {}",
            n.iterations,
            b.iterations
        );
    }

    // --- secant needs no derivative -----------------------------------

    #[test]
    fn secant_uses_only_f() {
        // The closure captures nothing derivative-like; if it compiles
        // and converges, the method genuinely needs no f'.
        let r = secant(|x| x.cos() - x, 0.0, 1.0, Settings::default()).unwrap();
        // Dottie number: the unique real root of cos(x) = x.
        assert!(
            (r.root - 0.739_085_133_215_160_6).abs() < EPS,
            "root = {}",
            r.root
        );
    }

    #[test]
    fn secant_distinct_seeds_required() {
        let err = secant(quad, 1.0, 1.0, Settings::default()).unwrap_err();
        assert!(matches!(err, RootError::Invalid { .. }));
    }

    // --- max-iteration error ------------------------------------------

    #[test]
    fn bisection_exhausts_budget() {
        // One halving cannot reach 1e-12 across [0, 2].
        let s = Settings::new(1e-12, 1).unwrap();
        let err = bisection(quad, 0.0, 2.0, s).unwrap_err();
        assert!(matches!(err, RootError::MaxIterations { .. }));
        assert_eq!(err.code(), "rootfind.max_iterations");
    }

    #[test]
    fn newton_exhausts_budget_on_tight_tol() {
        // A single Newton step from x0 = 1 lands at 1.5 (f = 0.25),
        // nowhere near 1e-12, so a one-step budget overruns.
        let s = Settings::new(1e-12, 1).unwrap();
        let err = newton(quad, quad_prime, 1.0, s).unwrap_err();
        match err {
            RootError::MaxIterations { max_iter, .. } => assert_eq!(max_iter, 1),
            other => panic!("expected MaxIterations, got {other:?}"),
        }
    }

    #[test]
    fn secant_exhausts_budget_on_tight_tol() {
        let s = Settings::new(1e-12, 1).unwrap();
        let err = secant(quad, 1.0, 2.0, s).unwrap_err();
        assert!(matches!(err, RootError::MaxIterations { .. }));
    }

    // --- zero-derivative guard ----------------------------------------

    #[test]
    fn newton_flags_zero_derivative() {
        // f' = 2x is zero at x0 = 0; the tangent is horizontal.
        let err = newton(quad, quad_prime, 0.0, Settings::default()).unwrap_err();
        assert!(matches!(err, RootError::ZeroDerivative { .. }));
        assert_eq!(err.code(), "rootfind.zero_derivative");
    }

    #[test]
    fn secant_flags_flat_chord() {
        // f(1) = f(-1) for the even function x^2 - 2, so the very first
        // chord through x0 = -1, x1 = 1 is horizontal.
        let err = secant(quad, -1.0, 1.0, Settings::default()).unwrap_err();
        assert!(matches!(err, RootError::ZeroDerivative { .. }));
    }

    // --- non-finite evaluation ----------------------------------------

    #[test]
    fn non_finite_evaluation_is_reported() {
        let err = bisection(|_| f64::NAN, 0.0, 1.0, Settings::default()).unwrap_err();
        assert!(matches!(err, RootError::NotFinite { .. }));
        assert_eq!(err.code(), "rootfind.not_finite");
    }

    // --- residual quality ---------------------------------------------

    #[test]
    fn reported_residual_matches_function() {
        let r = newton(quad, quad_prime, 1.0, Settings::default()).unwrap();
        assert!((r.residual - quad(r.root).abs()).abs() < 1e-15);
    }

    // --- a different test function: a cubic ----------------------------

    #[test]
    fn all_methods_agree_on_cubic_root() {
        // x^3 - x - 2 has a single real root near 1.5213797...
        let cubic = |x: f64| x * x * x - x - 2.0;
        let cubic_prime = |x: f64| 3.0 * x * x - 1.0;
        let truth = 1.521_379_706_804_567_6;

        let b = bisection(cubic, 1.0, 2.0, Settings::default()).unwrap();
        let n = newton(cubic, cubic_prime, 1.5, Settings::default()).unwrap();
        let s = secant(cubic, 1.0, 2.0, Settings::default()).unwrap();

        assert!((b.root - truth).abs() < EPS, "bisection {}", b.root);
        assert!((n.root - truth).abs() < EPS, "newton {}", n.root);
        assert!((s.root - truth).abs() < EPS, "secant {}", s.root);
    }

    // --- method identifiers -------------------------------------------

    #[test]
    fn method_identifiers() {
        assert_eq!(Method::Bisection.as_str(), "bisection");
        assert_eq!(Method::Newton.as_str(), "newton");
        assert_eq!(Method::Secant.as_str(), "secant");
    }

    // --- serde round-trip ---------------------------------------------

    #[test]
    fn root_serde_round_trip() {
        let r = secant(quad, 1.0, 2.0, Settings::default()).unwrap();
        let json = serde_json::to_string(&r).unwrap();
        let back: Root = serde_json::from_str(&json).unwrap();
        assert_eq!(r, back);
    }
}
