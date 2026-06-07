//! The human-readable wind-tunnel run report.
//!
//! An [`AeroResult`] is a big struct; an [`AeroReport`] is the
//! one-page summary a human (or an LLM relaying to a human) actually
//! reads — the headline coefficients, the case conditions, the
//! convergence verdict, and the honest caveats that qualify the
//! numbers.

use crate::api::AeroResult;
use crate::compressible::FlowRegime;
use crate::forces::turbulence_note;

/// A structured, human-readable summary of a wind-tunnel run.
#[derive(Clone, Debug)]
pub struct AeroReport {
    /// The drag coefficient.
    pub cd: f64,
    /// The lift coefficient.
    pub cl: f64,
    /// The side-force coefficient.
    pub cs: f64,
    /// The pitch-moment coefficient.
    pub cm: f64,
    /// The pressure-drag fraction of the total drag (0–1).
    pub pressure_drag_fraction: f64,
    /// The drag area `Cd·A` (m²).
    pub drag_area: f64,
    /// The free-stream dynamic pressure `q∞ = ½·ρ·U∞²` (Pa) of the run.
    pub dynamic_pressure: f64,
    /// The reference area `A` (m²) the coefficients normalise against — the
    /// body's frontal silhouette.
    pub reference_area: f64,
    /// The free-stream Reynolds number.
    pub reynolds_number: f64,
    /// The free-stream Mach number.
    pub mach_number: f64,
    /// Whether the solve converged.
    pub converged: bool,
    /// The number of solver iterations performed.
    pub iterations: usize,
    /// The final residual.
    pub residual: f64,
    /// The mean dimensionless wall distance `y+` over the body.
    pub y_plus_mean: f64,
    /// The peak velocity deficit in the wake.
    pub wake_peak_deficit: f64,
    /// Honest caveat lines qualifying the result.
    pub caveats: Vec<String>,
}

/// The Prandtl–Glauert compressibility factor `1/√(1 − M²)` at Mach `mach` — the
/// subsonic correction that scales a thin-body's incompressible aerodynamic
/// coefficients toward their compressible values (`Cl ≈ Cl_incompressible ·
/// factor`). It is `1` at `M = 0` and grows toward `M = 1`, where the linearised
/// theory diverges; valid only for `0 ≤ M < 1`. Returns `0` outside that range
/// (sonic / supersonic / non-finite), where the correction does not apply.
pub fn prandtl_glauert_factor(mach: f64) -> f64 {
    if mach.is_finite() && (0.0..1.0).contains(&mach) {
        1.0 / (1.0 - mach * mach).sqrt()
    } else {
        0.0
    }
}

/// The finite-wing (3-D) lift-curve slope `a = a₀ / (1 + a₀/(π·e·AR))` (per
/// radian) from Prandtl's lifting-line theory — how a wing of finite aspect
/// ratio `aspect_ratio` and span efficiency `span_efficiency` develops a
/// *gentler* lift slope than its 2-D airfoil section `section_slope_per_rad`
/// (typically ≈ `2π`/rad). The downwash induced by the trailing vortices tilts
/// the local flow and cuts the effective angle of attack, so a given incidence
/// makes less lift. This is the lift-side companion to the induced drag
/// (`crate::sweep::PolarCurve::induced_drag_factor`) — both scale with the same
/// `1/(π·e·AR)` finite-span group. As `AR → ∞` the downwash vanishes and the
/// 2-D slope `a₀` is recovered. Returns `0` for non-physical inputs
/// (`a₀ < 0`, `AR ≤ 0`, `e ≤ 0`, or any non-finite), where the relation does
/// not apply.
pub fn finite_wing_lift_slope(
    section_slope_per_rad: f64,
    aspect_ratio: f64,
    span_efficiency: f64,
) -> f64 {
    if !section_slope_per_rad.is_finite()
        || !aspect_ratio.is_finite()
        || !span_efficiency.is_finite()
        || section_slope_per_rad < 0.0
        || aspect_ratio <= 0.0
        || span_efficiency <= 0.0
    {
        return 0.0;
    }
    section_slope_per_rad
        / (1.0 + section_slope_per_rad / (std::f64::consts::PI * span_efficiency * aspect_ratio))
}

/// The level-flight **stall speed** `V_stall = √(2W / (ρ·S·C_Lmax))` (m/s) — the
/// slowest speed at which a wing of area `wing_area_m2` (m²) can still carry the
/// weight `weight_n` (N) in level flight at air density `air_density` (kg/m³),
/// flying at its maximum lift coefficient `cl_max`. It is exactly the speed at
/// which the available lift `L = ½ρV²S·C_Lmax` equals the weight; any slower and
/// the wing cannot make enough lift and the aircraft stalls. This is the
/// reference speed behind the approach and landing speeds (typically flown at
/// `1.2–1.3·V_stall`). Returns `0` for any non-physical input (non-finite or
/// non-positive).
pub fn stall_speed(weight_n: f64, wing_area_m2: f64, air_density: f64, cl_max: f64) -> f64 {
    if !weight_n.is_finite()
        || weight_n <= 0.0
        || !wing_area_m2.is_finite()
        || wing_area_m2 <= 0.0
        || !air_density.is_finite()
        || air_density <= 0.0
        || !cl_max.is_finite()
        || cl_max <= 0.0
    {
        return 0.0;
    }
    (2.0 * weight_n / (air_density * wing_area_m2 * cl_max)).sqrt()
}

/// The **Breguet range** `R = (V/c)·(L/D)·ln(W₀/W₁)` (m) of a jet in cruise —
/// the classic flight-performance result for how far an aircraft flies burning
/// its fuel. `velocity` `V` (m/s) is the cruise speed, `sfc` `c` (s⁻¹) the
/// thrust-specific fuel consumption (fuel weight per unit thrust per second),
/// `lift_to_drag` `L/D` the aerodynamic efficiency, and `weight_ratio` `W₀/W₁`
/// the start-to-end (takeoff-to-landing) weight ratio — the fuel fraction in
/// log form. It is maximised by cruising at the best-`L/D` point
/// ([`crate::sweep::PolarCurve::best_lift_to_drag_point`]) at the highest
/// `V`/`c`. Returns `0` when no fuel is burned (`W₀/W₁ = 1`) and for any
/// non-physical input (`V`, `c`, or `L/D` non-positive, `W₀/W₁ < 1`, non-finite).
pub fn breguet_range(velocity: f64, sfc: f64, lift_to_drag: f64, weight_ratio: f64) -> f64 {
    if !velocity.is_finite()
        || velocity <= 0.0
        || !sfc.is_finite()
        || sfc <= 0.0
        || !lift_to_drag.is_finite()
        || lift_to_drag <= 0.0
        || !weight_ratio.is_finite()
        || weight_ratio < 1.0
    {
        return 0.0;
    }
    velocity / sfc * lift_to_drag * weight_ratio.ln()
}

/// The **Mach angle** `μ = arcsin(1/M)` (radians) — the half-angle of the Mach
/// cone trailing a body moving at Mach number `mach`, the line along which the
/// pressure disturbances pile into a shock. It is the *supersonic* companion to
/// the subsonic [`prandtl_glauert_factor`]: at `M = 1` the cone is a flat disc
/// (`μ = π/2`), and as the speed rises the cone narrows toward `μ → 0`, sweeping
/// ever more sharply back. There is no Mach cone below the speed of sound, so it
/// returns `0` for subsonic `M < 1` (or non-finite input), exactly as the
/// Prandtl–Glauert factor returns `0` outside its valid range.
pub fn mach_angle(mach: f64) -> f64 {
    if mach.is_finite() && mach >= 1.0 {
        (1.0 / mach).asin()
    } else {
        0.0
    }
}

/// The **dynamic-pressure ratio** `q/p = ½·γ·M²` — the compressible dynamic pressure
/// `q = ½ρV²` normalised by the static pressure `p`, at Mach number `mach` `M` and
/// heat-capacity ratio `gamma` `γ`. Because `ρV² = γ·p·M²` (from `a² = γp/ρ`), the
/// dynamic head is `q = ½γpM²`, so `q/p` depends only on `M` and `γ`.
///
/// In the low-Mach (incompressible) limit it recovers Bernoulli's `p₀ = p + q`: the
/// stagnation excess [`isentropic_stagnation_pressure_ratio`] `p₀/p − 1 → ½γM²`. At
/// higher Mach the compressible stagnation excess exceeds `q/p`. It is `0` at rest and
/// grows with the square of the Mach number. Returns `0` for non-physical input
/// (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn dynamic_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    0.5 * gamma * mach * mach
}

/// The **Prandtl–Meyer function** `ν(M)` (radians) — the angle through which a
/// supersonic stream turns, in a centred expansion fan, as it accelerates
/// isentropically from `M = 1` to Mach number `mach` `M` (heat-capacity ratio
/// `gamma` `γ`):
///
/// ```text
///   ν = √((γ+1)/(γ−1))·atan√((γ−1)/(γ+1)·(M²−1)) − atan√(M²−1)
/// ```
///
/// It is the workhorse of supersonic *expansion* (and the method of
/// characteristics): a wall turning away from the flow expands it, raising the Mach
/// number by the `Δν` swept across the fan. `ν(1) = 0` (no turning at the sonic
/// line) and it rises monotonically, asymptoting to the finite maximum
/// `ν_max = (π/2)·(√((γ+1)/(γ−1)) − 1)` (≈ 130.5° for air, `γ = 1.4`) as `M → ∞` —
/// the most a stream can ever turn while expanding into a vacuum. It is the
/// expansion-side companion to the compression-side [`mach_angle`]; at `M = 2`
/// (`γ = 1.4`) `ν ≈ 26.4°`, at `M = 3` `≈ 49.8°`. Returns `0` for subsonic `M < 1`
/// (no expansion fan) or non-physical input (non-finite `M` or `γ`, or `γ ≤ 1`).
pub fn prandtl_meyer_angle(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 1.0 {
        return 0.0;
    }
    let s = (mach * mach - 1.0).sqrt();
    let k = ((gamma + 1.0) / (gamma - 1.0)).sqrt();
    k * (s / k).atan() - s.atan()
}

/// The **isentropic stagnation temperature ratio** `T₀/T = 1 + ((γ−1)/2)·M²` at
/// Mach number `mach` `M` and heat-capacity ratio `gamma` `γ` — the total-to-
/// static temperature relation that follows directly from adiabatic energy
/// conservation (`cₚ·T₀ = cₚ·T + ½V²`). `T₀` is the temperature the flow reaches
/// when it is brought to rest: the **recovery / stagnation temperature** that
/// drives aerodynamic heating and sizes a high-speed vehicle's thermal
/// protection, the thermal counterpart to the Pitot-airspeed role of the
/// [`isentropic_stagnation_pressure_ratio`].
///
/// It is the more fundamental member of the isentropic-stagnation pair — the
/// pressure ratio is exactly `(T₀/T)^(γ/(γ−1))` — and, like it, is finite and
/// well-behaved across the **whole** range `M ≥ 0`, subsonic and supersonic
/// alike: `1` at `M = 0` (no kinetic energy to recover) and rising linearly in
/// `M²` (a Mach-2 stream stagnates ~80 % hotter in absolute terms; Mach 5, six-
/// fold). Returns `1.0` (the no-rise identity) for non-physical input (non-finite
/// `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn isentropic_stagnation_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 1.0;
    }
    1.0 + 0.5 * (gamma - 1.0) * mach * mach
}

/// The **isentropic stagnation pressure ratio** `p₀/p = (1 + ((γ−1)/2)·M²)^(γ/(γ−1))`
/// at Mach number `mach` `M` and heat-capacity ratio `gamma` `γ` — the exact
/// compressible total-to-static pressure relation for an adiabatic, reversible
/// (isentropic) deceleration of the flow to rest. This is the compressible
/// Pitot law: it is what converts a measured total-minus-static pressure into
/// airspeed once the flow is fast enough that the incompressible `½ρV²` no
/// longer holds.
///
/// It is the *exact* compressibility relation that complements the two
/// linearised/limiting Mach functions here — [`prandtl_glauert_factor`], the
/// thin-body subsonic correction that diverges at `M = 1`, and [`mach_angle`],
/// the supersonic cone half-angle — whereas this stagnation ratio is finite and
/// well-behaved across the **whole** range `M ≥ 0`, subsonic and supersonic
/// alike. It is `1` at `M = 0` (a body at rest compresses nothing) and rises
/// monotonically with Mach; at low speed it reduces to `1 + (γ/2)·M²`, the
/// leading compressible form of the dynamic-pressure rise. Returns `1.0` (the
/// no-correction identity) for non-physical input (non-finite `M` or `γ`,
/// `M < 0`, or `γ ≤ 1`).
pub fn isentropic_stagnation_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 1.0;
    }
    let temperature_ratio = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    temperature_ratio.powf(gamma / (gamma - 1.0))
}

/// The **isentropic stagnation density ratio** `ρ₀/ρ = (1 + ((γ−1)/2)·M²)^(1/(γ−1))`
/// at Mach number `mach` `M` and heat-capacity ratio `gamma` `γ` — the total-to-
/// static *density* relation for an adiabatic, reversible (isentropic) deceleration
/// of the flow to rest. `ρ₀` is the density the gas reaches when brought to rest;
/// this is the third member of the isentropic-stagnation trio, completing the
/// [`isentropic_stagnation_temperature_ratio`] and the
/// [`isentropic_stagnation_pressure_ratio`].
///
/// The three are locked together by the perfect-gas isentrope: writing
/// `τ = T₀/T = 1 + ((γ−1)/2)·M²`, the density ratio is `τ^(1/(γ−1))` and the
/// pressure ratio is `τ^(γ/(γ−1))`, so `p₀/p = (ρ₀/ρ)^γ = (T₀/T)^(γ/(γ−1))` — the
/// polytropic chain `p ∝ ρ^γ ∝ T^(γ/(γ−1))`. Like its companions it is finite and
/// well-behaved across the **whole** range `M ≥ 0`, subsonic and supersonic alike:
/// `1` at `M = 0` (a body at rest compresses nothing) and rising monotonically with
/// Mach (a sonic air stream, `γ = 1.4`, stagnates ~58 % denser — `ρ₀/ρ = 1.2^2.5 ≈
/// 1.577`). It is the quantity behind the density rise in a Pitot/total-condition
/// reduction and in sizing the mass flux of a high-speed intake. Returns `1.0` (the
/// no-compression identity) for non-physical input (non-finite `M` or `γ`, `M < 0`,
/// or `γ ≤ 1`).
pub fn isentropic_stagnation_density_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 1.0;
    }
    let temperature_ratio = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    temperature_ratio.powf(1.0 / (gamma - 1.0))
}

/// The **isentropic area ratio**
/// `A/A* = (1/M)·[(2/(γ+1))·(1 + ((γ−1)/2)·M²)]^((γ+1)/(2(γ−1)))` at Mach number
/// `mach` `M` and heat-capacity ratio `gamma` `γ` — the ratio of the local duct
/// area to the sonic-throat area `A*` in a 1-D isentropic (de Laval) nozzle: the
/// area a streamtube must have to carry the flow from the choked throat (`M = 1`)
/// to Mach `M`, the foundational relation of converging–diverging nozzle design.
///
/// Unlike the monotonic stagnation ratios, `A/A*` has a **minimum of exactly 1 at
/// the sonic throat** `M = 1` and rises on *both* sides: a converging duct
/// accelerates a subsonic flow toward the throat, and only a *diverging* duct
/// downstream accelerates it supersonically — so every area ratio `> 1` is shared
/// by one subsonic and one supersonic solution. It diverges as `M → 0` (an
/// infinite reservoir feeds the throat); at `M = 2` (`γ = 1.4`) the area is
/// `1.6875·A*`, at `M = 3` `4.2346·A*`. Returns `f64::INFINITY` — the `M → 0`
/// limit — for a zero/negative or non-physical Mach (non-finite `M` or `γ`,
/// `M < 0`, or `γ ≤ 1`), distinguishing it from the finite `≥ 1` values of a real
/// `M > 0` flow.
pub fn isentropic_area_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return f64::INFINITY;
    }
    let temperature_ratio = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    let exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    (1.0 / mach) * (2.0 / (gamma + 1.0) * temperature_ratio).powf(exponent)
}

/// The **compressible mass-flow function**
/// `FF = √γ · M · (1 + ((γ−1)/2)·M²)^(−(γ+1)/(2(γ−1)))` at Mach number `mach` `M`
/// and heat-capacity ratio `gamma` `γ` — the non-dimensional mass flow
/// `ṁ·√(R·T₀)/(A·p₀)` an isentropic stream of total temperature `T₀` and total
/// pressure `p₀` carries through area `A`. It is the relation behind nozzle and
/// turbomachinery sizing: for fixed `A`, `p₀`, `T₀` the mass flow is set by `M`.
///
/// Unlike the monotonic stagnation ratios it **peaks at exactly `M = 1`** — the
/// *choking* condition: a converging duct accelerates the flow only up to the sonic
/// throat, where the mass flux is maximal, so no subsonic duct can pass more than
/// `FF(1)·A·p₀/√(R·T₀)`. It is `0` at rest (`M = 0`) and falls again on the
/// supersonic branch, so each value `< FF(1)` is shared by one subsonic and one
/// supersonic solution — the reciprocal face of [`isentropic_area_ratio`], to which
/// it is tied by mass conservation `A/A* = FF(1)/FF(M)`. For `γ = 1.4` the choked
/// peak is `FF(1) = √γ·(2/(γ+1))^((γ+1)/(2(γ−1))) ≈ 0.6847`. Returns `0` for a
/// non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`) — the `M = 0`
/// no-flow limit.
pub fn mass_flow_function(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 0.0;
    }
    let temperature_ratio = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    let exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    gamma.sqrt() * mach * temperature_ratio.powf(-exponent)
}

/// The **characteristic (sonic-referenced) Mach number** `M* = V/a* =
/// M·√((γ+1) / (2 + (γ−1)·M²))` at Mach number `mach` `M` and heat-capacity ratio
/// `gamma` `γ` — the flow speed `V` measured against the *critical* (sonic) speed
/// `a*` rather than the local sound speed `a` (which is what the ordinary Mach
/// number `M = V/a` uses). Because `a*` stays constant along an adiabatic flow while
/// the local `a` varies, `M*` is the more convenient speed variable in several
/// compressible-flow relations.
///
/// It is `0` at rest, crosses `1` exactly at the sonic point (`M = 1`, where
/// `V = a = a*`) so it labels subsonic/supersonic the same way `M` does, and —
/// unlike `M`, which is unbounded — **saturates** at the finite limit
/// `√((γ+1)/(γ−1))` (`≈ 2.449` for air) as `M → ∞`. Its signature use is the
/// **Prandtl relation** for a normal shock, `M₁*·M₂* = 1`: the up- and downstream
/// characteristic Mach numbers are reciprocals (with `M₂` from
/// [`normal_shock_downstream_mach`]). Returns `0` for non-physical input
/// (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn characteristic_mach(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 0.0;
    }
    mach * ((gamma + 1.0) / (2.0 + (gamma - 1.0) * mach * mach)).sqrt()
}

/// The **downstream Mach number** `M₂` behind a stationary **normal shock** with
/// upstream Mach `mach` `M₁` and heat-capacity ratio `gamma` `γ`, from the
/// Rankine–Hugoniot jump conditions:
///
/// ```text
///   M₂² = (1 + ((γ−1)/2)·M₁²) / (γ·M₁² − (γ−1)/2)
/// ```
///
/// Unlike the isentropic stagnation relations above — which describe a *smooth,
/// reversible* deceleration — a normal shock is an abrupt, **irreversible**
/// (entropy-increasing) compression that forms only in **supersonic** flow, and
/// it always leaves the flow **subsonic** (`M₂ < 1` for every `M₁ > 1`). The
/// downstream Mach falls as the shock strengthens, approaching the finite
/// strong-shock limit `√((γ−1)/2γ)` (≈ `0.378` for air) as `M₁ → ∞`. This is the
/// foundational shock relation the static pressure, density and temperature jumps
/// are all built on.
///
/// `M₁ = 1` is the infinitesimal (no-shock) limit, `M₂ = 1`. For **subsonic or
/// sonic** upstream (`M₁ ≤ 1`) no shock forms and the flow passes through
/// unchanged, so the input `mach` is returned. Returns `1.0` (the sonic identity)
/// for non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn normal_shock_downstream_mach(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 1.0;
    }
    if mach <= 1.0 {
        return mach; // subsonic/sonic: no shock forms, the flow is unchanged
    }
    let numerator = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    let denominator = gamma * mach * mach - 0.5 * (gamma - 1.0);
    (numerator / denominator).sqrt()
}

/// The static **pressure ratio** `p₂/p₁ = 1 + (2γ/(γ+1))·(M₁²−1)` across a
/// stationary **normal shock** with upstream Mach `mach` `M₁` and heat-capacity
/// ratio `gamma` `γ` — the Rankine–Hugoniot static-pressure jump, the companion
/// to [`normal_shock_downstream_mach`] in the shock-relations family. (Distinct
/// from the *stagnation* ratio [`isentropic_stagnation_pressure_ratio`], which is
/// the reversible total-to-static relation; this is the irreversible jump across
/// the shock itself.)
///
/// A shock always **compresses**, so `p₂/p₁ > 1` for any `M₁ > 1` and it rises
/// without bound as the shock strengthens (`∝ M₁²` as `M₁ → ∞`) — unlike the
/// density jump, which saturates at `(γ+1)/(γ−1)`. `M₁ = 1` is the infinitesimal
/// (no-shock) limit, `p₂/p₁ = 1`; a Mach-2 shock in air (`γ = 1.4`) raises the
/// static pressure 4.5-fold, a Mach-3 shock 10.3-fold. For **subsonic or sonic**
/// upstream (`M₁ ≤ 1`) no shock forms and the pressure is unchanged (`1.0`).
/// Returns `1.0` (the no-jump identity) for non-physical input (non-finite `M`
/// or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn normal_shock_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 1.0;
    }
    if mach <= 1.0 {
        return 1.0; // subsonic/sonic: no shock forms, the static pressure is unchanged
    }
    1.0 + 2.0 * gamma / (gamma + 1.0) * (mach * mach - 1.0)
}

/// The static **density ratio** `ρ₂/ρ₁ = (γ+1)·M₁² / ((γ−1)·M₁² + 2)` across a
/// stationary **normal shock** with upstream Mach `mach` `M₁` and heat-capacity
/// ratio `gamma` `γ` — the Rankine–Hugoniot density jump, the third member of the
/// shock-relations family with [`normal_shock_downstream_mach`] and
/// [`normal_shock_pressure_ratio`]. (Distinct from the *stagnation* ratio
/// [`isentropic_stagnation_density_ratio`], the reversible total-to-static
/// relation; this is the irreversible jump across the shock itself.)
///
/// Unlike the pressure jump — which grows without bound (`∝ M₁²`) — the density
/// jump **saturates** at the finite strong-shock limit `(γ+1)/(γ−1)` (`= 6` for
/// air, `γ = 1.4`) as `M₁ → ∞`: a gas can be compressed only so far, because the
/// post-shock temperature (and the thermal motion resisting further packing)
/// climbs without limit. `M₁ = 1` is the no-shock limit (`ρ₂/ρ₁ = 1`); a Mach-2
/// shock in air compresses the gas ~2.67-fold, a Mach-3 shock ~3.86-fold. For
/// **subsonic or sonic** upstream (`M₁ ≤ 1`) no shock forms and the density is
/// unchanged (`1.0`). Returns `1.0` (the no-jump identity) for non-physical input
/// (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn normal_shock_density_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 1.0;
    }
    if mach <= 1.0 {
        return 1.0; // subsonic/sonic: no shock forms, the density is unchanged
    }
    (gamma + 1.0) * mach * mach / ((gamma - 1.0) * mach * mach + 2.0)
}

/// The static **temperature ratio** `T₂/T₁ = [2γM₁² − (γ−1)]·[(γ−1)M₁² + 2] /
/// ((γ+1)²·M₁²)` across a stationary **normal shock** with upstream Mach `mach`
/// `M₁` and heat-capacity ratio `gamma` `γ` — the Rankine–Hugoniot temperature
/// jump, completing the static-property trio with [`normal_shock_pressure_ratio`]
/// and [`normal_shock_density_ratio`]. By the ideal-gas law `T = p/(ρR)` it is
/// exactly their quotient, `T₂/T₁ = (p₂/p₁)/(ρ₂/ρ₁)`. (Distinct from the
/// *stagnation* ratio [`isentropic_stagnation_temperature_ratio`], the reversible
/// total-to-static relation; this is the irreversible jump across the shock
/// itself.)
///
/// Unlike the density jump — which saturates at `(γ+1)/(γ−1)` — the temperature
/// jump **grows without bound** (`∝ M₁²` as `M₁ → ∞`): the kinetic energy of the
/// supersonic stream is dumped irreversibly into heat. `M₁ = 1` is the no-shock
/// limit (`T₂/T₁ = 1`); a Mach-2 shock in air (`γ = 1.4`) raises the static
/// temperature 1.6875-fold, a Mach-3 shock ~2.68-fold. For **subsonic or sonic**
/// upstream (`M₁ ≤ 1`) no shock forms and the temperature is unchanged (`1.0`).
/// Returns `1.0` (the no-jump identity) for non-physical input (non-finite `M` or
/// `γ`, `M < 0`, or `γ ≤ 1`).
pub fn normal_shock_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 1.0;
    }
    if mach <= 1.0 {
        return 1.0; // subsonic/sonic: no shock forms, the static temperature is unchanged
    }
    let m2 = mach * mach;
    (2.0 * gamma * m2 - (gamma - 1.0)) * ((gamma - 1.0) * m2 + 2.0)
        / ((gamma + 1.0) * (gamma + 1.0) * m2)
}

/// The **stagnation-pressure ratio** `p₀₂/p₀₁` across a stationary **normal shock**
/// with upstream Mach `mach` `M₁` and heat-capacity ratio `gamma` `γ` — the
/// total-pressure *recovery*, the canonical measure of the shock's irreversible
/// loss:
///
/// ```text
///   p₀₂/p₀₁ = [ (γ+1)M₁² / ((γ−1)M₁² + 2) ]^(γ/(γ−1))
///           · [ (γ+1) / (2γM₁² − (γ−1)) ]^(1/(γ−1))
/// ```
///
/// Unlike the *static* pressure jump [`normal_shock_pressure_ratio`] (which rises
/// without bound), the total pressure is always **lost** across a shock:
/// `p₀₂/p₀₁ < 1` for any `M₁ > 1`, falling monotonically as the shock strengthens
/// (a Mach-2 shock in air, `γ = 1.4`, recovers ~72%, a Mach-3 shock only ~33%).
/// The loss is the thermodynamic signature of the entropy the shock generates,
/// `p₀₂/p₀₁ = e^(−Δs/R)`, and is the headline figure of merit for a supersonic
/// inlet/diffuser. (Distinct from the *isentropic*
/// [`isentropic_stagnation_pressure_ratio`], the reversible total-to-static
/// `p₀/p` of a single stream; this is total-to-total *across* the irreversible
/// jump.) For **subsonic or sonic** upstream (`M₁ ≤ 1`) no shock forms and the
/// total pressure is conserved (`1.0`). Returns `1.0` (the no-loss identity) for
/// non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn normal_shock_stagnation_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 1.0;
    }
    if mach <= 1.0 {
        return 1.0; // subsonic/sonic: no shock forms, the total pressure is conserved
    }
    let m2 = mach * mach;
    let density_term =
        ((gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0)).powf(gamma / (gamma - 1.0));
    let pressure_term =
        ((gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0))).powf(1.0 / (gamma - 1.0));
    density_term * pressure_term
}

/// The **specific entropy rise** `Δs/R` (dimensionless, in units of the gas
/// constant `R`) across a stationary **normal shock** with upstream Mach `mach`
/// `M₁` and heat-capacity ratio `gamma` `γ` — the second-law signature of the
/// shock's irreversibility. A shock is adiabatic but *not* isentropic: it turns
/// ordered supersonic kinetic energy into heat, generating entropy
///
/// ```text
///   Δs/R = (γ/(γ−1))·ln(T₂/T₁) − ln(p₂/p₁)
/// ```
///
/// from the Rankine–Hugoniot static [`normal_shock_temperature_ratio`] and
/// [`normal_shock_pressure_ratio`] jumps. This entropy production is exactly what
/// drives the total-pressure loss [`normal_shock_stagnation_pressure_ratio`]:
/// `p₀₂/p₀₁ = e^(−Δs/R)`, so equivalently `Δs/R = −ln(p₀₂/p₀₁)`. It is `0` at
/// `M₁ = 1` (a vanishingly weak shock is reversible) and **grows monotonically**
/// with shock strength — the thermodynamic reason a stronger shock recovers less
/// total pressure. For **subsonic or sonic** upstream (`M₁ ≤ 1`) no shock forms and
/// no entropy is generated (`0`). Returns `0` for non-physical input (non-finite
/// `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn normal_shock_entropy_rise(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 0.0;
    }
    if mach <= 1.0 {
        return 0.0; // subsonic/sonic: no shock forms, no entropy is generated
    }
    gamma / (gamma - 1.0) * normal_shock_temperature_ratio(mach, gamma).ln()
        - normal_shock_pressure_ratio(mach, gamma).ln()
}

/// The **Rayleigh supersonic pitot ratio** `p₀₂/p₁` — the total pressure a pitot
/// probe reads relative to the freestream *static* pressure at upstream Mach `mach`
/// `M₁` (heat-capacity ratio `gamma` `γ`). It is the working formula of a
/// supersonic pitot tube: above `M = 1` the probe sits behind its own detached bow
/// shock, so it senses the *post-shock* total pressure, not the freestream total —
/// the Rayleigh pitot formula
///
/// ```text
///   p₀₂/p₁ = [ (γ+1)M₁²/2 ]^(γ/(γ−1)) · [ (γ+1)/(2γM₁² − (γ−1)) ]^(1/(γ−1))   (M₁ > 1)
/// ```
///
/// — what inverts a measured pitot-to-static ratio back into a supersonic Mach
/// number. By construction it is the product of the across-shock total-pressure
/// recovery [`normal_shock_stagnation_pressure_ratio`] and the post-shock isentropic
/// rise [`isentropic_stagnation_pressure_ratio`]: `p₀₂/p₁ = (p₀₂/p₀₁)·(p₀₁/p₁)`. For
/// **subsonic** flow (`M ≤ 1`) no shock forms and it reduces to the ordinary
/// isentropic total-to-static ratio [`isentropic_stagnation_pressure_ratio`]; at
/// `M = 2` (`γ = 1.4`) a pitot reads ~5.64× the static pressure. Returns `1.0` for
/// non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn rayleigh_pitot_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || gamma <= 1.0 || mach < 0.0 {
        return 1.0;
    }
    if mach <= 1.0 {
        // Subsonic: no bow shock — the pitot reads the isentropic total pressure.
        return isentropic_stagnation_pressure_ratio(mach, gamma);
    }
    let m2 = mach * mach;
    let total_term = ((gamma + 1.0) * m2 / 2.0).powf(gamma / (gamma - 1.0));
    let shock_term =
        ((gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0))).powf(1.0 / (gamma - 1.0));
    total_term * shock_term
}

/// The **Rayleigh-flow stagnation-temperature ratio**
/// `T₀/T₀* = 2(γ+1)M²(1 + ((γ−1)/2)·M²) / (1 + γM²)²` at Mach number `mach` `M` and
/// heat-capacity ratio `gamma` `γ` — the driving relation of **Rayleigh flow**
/// (frictionless, constant-area duct flow with heat addition), the heat-addition
/// counterpart to the isentropic and normal-shock toolkits.
///
/// Adding heat raises the stagnation temperature `T₀` and drives the Mach number
/// toward `1` from *both* sides, so `T₀/T₀*` against the sonic reference `T₀*`
/// **peaks at exactly `1` at `M = 1`** — the *thermal-choking* limit: a given duct
/// can accept only enough heat to bring the flow to sonic, no more. It is `0` at
/// rest (`M = 0`) and rises on both the subsonic and supersonic branches to that
/// shared maximum (a Mach-2 stream sits at `0.793·T₀*`, Mach 0.5 at `0.691`).
/// Returns `0` for non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn rayleigh_flow_total_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    let m2 = mach * mach;
    let denom = 1.0 + gamma * m2;
    2.0 * (gamma + 1.0) * m2 * (1.0 + 0.5 * (gamma - 1.0) * m2) / (denom * denom)
}

/// The **Rayleigh-flow static temperature ratio** `T/T* = (M(γ+1)/(1+γM²))²` at
/// Mach number `mach` `M` and heat-capacity ratio `gamma` `γ` — the static
/// temperature of a frictionless heat-addition (Rayleigh) flow relative to its sonic
/// value `T*`, the static-property companion to the stagnation
/// [`rayleigh_flow_total_temperature_ratio`].
///
/// Its hallmark is a **maximum at `M = 1/√γ`** (≈ 0.845 for air), *below* the sonic
/// point: as heat is added to a subsonic stream the static temperature rises only
/// until `M = 1/√γ`, then *falls* even as more heat keeps driving the total
/// temperature `T₀` toward the thermal-choking limit — past that Mach the flow
/// accelerates faster than it heats, converting the added enthalpy into kinetic
/// energy. `T/T* = 1` at the sonic point `M = 1`, and `0` at rest (`M = 0`). Returns
/// `0` for non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn rayleigh_flow_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    let ratio = mach * (gamma + 1.0) / (1.0 + gamma * mach * mach);
    ratio * ratio
}

/// The **Rayleigh-flow static-pressure ratio** `p/p* = (1+γ)/(1 + γ·M²)` at Mach
/// number `mach` `M` and heat-capacity ratio `gamma` `γ` — the static pressure
/// referenced to the sonic (thermally choked) state `p*` for **Rayleigh flow**:
/// steady, frictionless, constant-area flow with **heat addition**.
///
/// It completes the Rayleigh static-state set alongside
/// [`rayleigh_flow_temperature_ratio`] (which is exactly `M²·(p/p*)²`) and
/// [`rayleigh_flow_total_temperature_ratio`]. Adding heat drives the Mach number
/// toward `1`: a subsonic stream sees its static pressure fall (`p/p* > 1`, down to
/// `1` at the choke) while a supersonic stream sees it rise (`p/p* < 1`). It is `1`
/// at the sonic point `M = 1` and tends to its maximum `1+γ` as `M → 0`. Returns `0`
/// for non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn rayleigh_flow_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    (1.0 + gamma) / (1.0 + gamma * mach * mach)
}

/// The **Rayleigh-flow velocity ratio** `V/V* = (1+γ)·M²/(1 + γ·M²)` at Mach number
/// `mach` `M` and heat-capacity ratio `gamma` `γ` — the flow speed referenced to the
/// sonic (thermally choked) state `V*` for **Rayleigh flow** (steady, frictionless,
/// constant-area flow with **heat addition**). By continuity `ρV = ρ*V*` it is also
/// the inverse density ratio `ρ*/ρ`.
///
/// It completes the Rayleigh static-state ratios alongside
/// [`rayleigh_flow_temperature_ratio`] and [`rayleigh_flow_pressure_ratio`], and is
/// exactly `M²·(p/p*)`. Adding heat drives the Mach number toward `1`, so the speed
/// rises toward `V*` from below in a subsonic stream (`V/V* < 1`) and falls toward it
/// from above in a supersonic stream (`V/V* > 1`); it is `1` at the sonic point and
/// `0` at rest. Returns `0` for non-physical input (non-finite `M` or `γ`, `M < 0`,
/// or `γ ≤ 1`).
pub fn rayleigh_flow_velocity_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    (1.0 + gamma) * mach * mach / (1.0 + gamma * mach * mach)
}

/// The **Rayleigh-flow stagnation (total) pressure ratio** `p₀/p₀* =
/// ((1+γ)/(1+γ·M²))·((2 + (γ−1)·M²)/(γ+1))^(γ/(γ−1))` at Mach number `mach` `M` and
/// heat-capacity ratio `gamma` `γ` — the stagnation pressure referenced to the sonic
/// (thermally choked) state `p₀*` for **Rayleigh flow** (steady, frictionless,
/// constant-area flow with **heat addition**).
///
/// It completes the Rayleigh family alongside [`rayleigh_flow_temperature_ratio`],
/// [`rayleigh_flow_pressure_ratio`] and [`rayleigh_flow_velocity_ratio`], and equals
/// `(p/p*)·((2+(γ−1)M²)/(γ+1))^(γ/(γ−1))`. It is `1` at the sonic point and `> 1`
/// everywhere else (a minimum at `M = 1`): adding heat always erodes the stagnation
/// pressure toward the choked value, the irreversible Rayleigh loss. Returns `0` for
/// non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn rayleigh_flow_total_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    let exponent = gamma / (gamma - 1.0);
    (gamma + 1.0) / (1.0 + gamma * mach * mach)
        * ((2.0 + (gamma - 1.0) * mach * mach) / (gamma + 1.0)).powf(exponent)
}

/// The **Fanno-flow static temperature ratio** `T/T* = (γ+1)/(2 + (γ−1)·M²)` at
/// Mach number `mach` `M` and heat-capacity ratio `gamma` `γ` — the static
/// temperature referenced to the sonic (choked) state `T*` for **Fanno flow**:
/// steady, adiabatic, constant-area flow with **wall friction**. It is the
/// friction-driven dual of the heat-addition
/// [`rayleigh_flow_temperature_ratio`]: friction, like heat addition, drives the
/// Mach number monotonically toward `1`, where the duct thermally/frictionally
/// chokes.
///
/// Because Fanno flow is adiabatic the **stagnation** temperature `T₀` is
/// conserved, so `T/T*` is fixed entirely by the isentropic static-to-total ratios
/// at `M` and at the sonic reference — `T/T* = (T₀/T*)/(T₀/T)` with
/// `T₀/T* = (γ+1)/2`. It is `(γ+1)/2` at rest (`M = 0`, the hottest static
/// temperature on a Fanno line), falls monotonically through `1` at the sonic
/// point (`M = 1`), and `→ 0` as `M → ∞`. Returns `0` for non-physical input
/// (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn fanno_flow_temperature_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    (gamma + 1.0) / (2.0 + (gamma - 1.0) * mach * mach)
}

/// The **Fanno-flow static pressure ratio** `p/p* = (1/M)·√((γ+1)/(2 + (γ−1)·M²))`
/// at Mach number `mach` `M` and heat-capacity ratio `gamma` `γ` — the static
/// pressure referenced to the sonic (choked) state for **Fanno flow** (steady,
/// adiabatic, constant-area flow with wall friction), the pressure companion to the
/// [`fanno_flow_temperature_ratio`] `T/T*`. It is exactly `√(T/T*)/M`, so it is `1`
/// at the sonic point (`M = 1`), exceeds `1` for subsonic flow, falls below `1` for
/// supersonic flow, and decreases monotonically as friction drives the Mach number
/// toward `1`. It diverges (`→ +∞`) as `M → 0` (the static pressure far from the
/// choked state). Returns `0` for non-physical input (non-finite `M` or `γ`,
/// `M < 0`, or `γ ≤ 1`).
pub fn fanno_flow_pressure_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    (1.0 / mach) * ((gamma + 1.0) / (2.0 + (gamma - 1.0) * mach * mach)).sqrt()
}

/// The **Fanno-flow velocity ratio** `V/V* = M·√((γ+1)/(2 + (γ−1)·M²))` at Mach
/// number `mach` `M` and heat-capacity ratio `gamma` `γ` — the flow speed referenced
/// to the sonic (choked) state for **Fanno flow** (steady, adiabatic, constant-area
/// flow with wall friction), the velocity companion to the
/// [`fanno_flow_temperature_ratio`] `T/T*` and [`fanno_flow_pressure_ratio`] `p/p*`.
/// It is exactly `M·√(T/T*) = M²·(p/p*)`, so it is `1` at the sonic point (`M = 1`),
/// approaches it from below (`V/V* < 1` for subsonic flow) and exceeds it from above
/// (`V/V* > 1` for supersonic) — the opposite trend to the temperature and pressure
/// ratios, since friction accelerates a subsonic stream and decelerates a supersonic
/// one toward `M = 1`. Unlike them it stays bounded as `M → ∞`, approaching the
/// maximum `√((γ+1)/(γ−1))`. Returns `0` at rest (`M = 0`, no flow) and for
/// non-physical input (non-finite `M` or `γ`, `M < 0`, or `γ ≤ 1`).
pub fn fanno_flow_velocity_ratio(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach < 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    mach * ((gamma + 1.0) / (2.0 + (gamma - 1.0) * mach * mach)).sqrt()
}

/// The **Fanno friction parameter** `4fL*/D` (dimensionless) at Mach number `mach`
/// `M` and heat-capacity ratio `gamma` `γ` — the friction factor times the duct
/// length-to-diameter ratio needed to drive **Fanno flow** (steady, adiabatic,
/// constant-area flow with wall friction) from `M` to the sonic choke at `M = 1`:
///
/// ```text
///   4fL*/D = (1 − M²)/(γ·M²) + ((γ+1)/(2γ))·ln[ (γ+1)·M² / (2 + (γ−1)·M²) ]
/// ```
///
/// It is the master Fanno design variable — the *maximum* run of duct a given inlet
/// Mach can take before it chokes — completing the Fanno family with
/// [`fanno_flow_temperature_ratio`], [`fanno_flow_pressure_ratio`] and
/// [`fanno_flow_velocity_ratio`]. It is `0` at the choke (`M = 1`, no length left)
/// and strictly positive on both sides — friction drives a subsonic stream up to
/// `M = 1` and a supersonic stream down to it — diverging as `M → 0`. Returns `0`
/// for non-physical input (non-finite `M` or `γ`, `M ≤ 0`, or `γ ≤ 1`).
pub fn fanno_friction_parameter(mach: f64, gamma: f64) -> f64 {
    if !mach.is_finite() || !gamma.is_finite() || mach <= 0.0 || gamma <= 1.0 {
        return 0.0;
    }
    (1.0 - mach * mach) / (gamma * mach * mach)
        + (gamma + 1.0) / (2.0 * gamma)
            * ((gamma + 1.0) * mach * mach / (2.0 + (gamma - 1.0) * mach * mach)).ln()
}

/// The **induced-drag coefficient** `C_Di = C_L² / (π·e·AR)` of a finite wing
/// (Prandtl lifting-line theory) — the unavoidable "drag-due-to-lift" that comes
/// with making lift at all. A wing of finite aspect ratio `aspect_ratio` `AR`
/// and span efficiency `span_efficiency` `e` sheds trailing vortices whose
/// downwash tilts the lift vector slightly aft, and that backward tilt *is* the
/// induced drag. It is the drag-side companion to [`finite_wing_lift_slope`]:
/// both carry the same `1/(π·e·AR)` finite-span group, here multiplied by the
/// square of the lift coefficient `lift_coefficient` `C_L` (so it quadruples
/// when `C_L` doubles, and is sign-blind — negative-lift downforce induces drag
/// just the same). Elliptical loading (`e = 1`) gives the theoretical minimum
/// `C_Di = C_L²/(π·AR)`, and the term vanishes for an infinite-span wing
/// (`AR → ∞`). Returns `0` at zero lift and for any non-physical input
/// (`AR ≤ 0`, `e ≤ 0`, or non-finite), where the relation does not apply.
pub fn induced_drag_coefficient(
    lift_coefficient: f64,
    aspect_ratio: f64,
    span_efficiency: f64,
) -> f64 {
    if !lift_coefficient.is_finite()
        || !aspect_ratio.is_finite()
        || !span_efficiency.is_finite()
        || aspect_ratio <= 0.0
        || span_efficiency <= 0.0
    {
        return 0.0;
    }
    lift_coefficient * lift_coefficient
        / (std::f64::consts::PI * span_efficiency * aspect_ratio)
}

/// The **maximum lift-to-drag ratio** `(L/D)_max = ½·√(π·e·AR / C_D0)` of a wing
/// with parabolic drag polar `C_D = C_D0 + C_L²/(π·e·AR)` — the single headline
/// aircraft-performance number, the best-glide and best-range/endurance
/// operating point. `cd0` is the zero-lift (parasite) drag coefficient `C_D0`,
/// `aspect_ratio` the wing `AR`, and `span_efficiency` the Oswald factor `e`.
///
/// It is reached at the lift coefficient `C_L* = √(π·e·AR·C_D0)`, the point where
/// the lift-induced drag (see [`induced_drag_coefficient`]) exactly equals the
/// parasite drag `C_D0`; the total drag is then `2·C_D0`, so
/// `(L/D)_max = C_L*/(2·C_D0)`. It improves with span efficiency and aspect ratio
/// (a long, clean wing glides far, `∝ √AR`) and degrades as the parasite drag
/// grows (`∝ 1/√C_D0`). This is the closed-form value from the polar
/// *parameters* — the analytic counterpart to
/// [`crate::sweep::PolarCurve::max_lift_to_drag`], which instead reads the peak
/// off a set of *measured* polar points. Returns `0` for any non-physical input
/// (`C_D0 ≤ 0`, `AR ≤ 0`, `e ≤ 0`, or non-finite).
pub fn max_lift_to_drag_ratio(cd0: f64, aspect_ratio: f64, span_efficiency: f64) -> f64 {
    if !cd0.is_finite()
        || cd0 <= 0.0
        || !aspect_ratio.is_finite()
        || aspect_ratio <= 0.0
        || !span_efficiency.is_finite()
        || span_efficiency <= 0.0
    {
        return 0.0;
    }
    0.5 * (std::f64::consts::PI * span_efficiency * aspect_ratio / cd0).sqrt()
}

impl AeroReport {
    /// Build a report from a completed [`AeroResult`].
    pub fn from_result(result: &AeroResult) -> AeroReport {
        let coeff = &result.coefficients;
        let total_cd = coeff.cd.abs().max(1e-12);
        let pressure_drag_fraction =
            (coeff.cd_pressure / total_cd).clamp(0.0, 1.0);

        let mut caveats = Vec::new();
        // The standing v1 caveat.
        caveats.push(
            "immersed-boundary Cartesian CFD — a real v1, not ANSYS Fluent / \
             STAR-CCM+ parity: no body-fitted meshing, no DES/LES, coefficient \
             accuracy improves with grid resolution"
                .to_string(),
        );
        // Turbulence model note.
        caveats.push(format!(
            "turbulence: {}",
            turbulence_note(result.flow.turbulence.model)
        ));
        // Convergence note.
        if !result.converged {
            caveats.push(format!(
                "the solve did NOT reach the convergence tolerance \
                 (residual {:.2e} after {} iterations) — treat the \
                 coefficients as provisional",
                result.flow.residual, result.flow.iterations
            ));
        }
        // y+ note — a wall function wants y+ in roughly 30–300.
        let yp = result.surface.y_plus_mean;
        if yp > 0.0 && !(1.0..=1000.0).contains(&yp) {
            caveats.push(format!(
                "mean y+ is {yp:.1} — outside the wall-function-friendly band; \
                 the near-wall resolution does not match the high-Re wall \
                 treatment well"
            ));
        }
        // Compressibility note.
        let regime = FlowRegime::classify(result.mach_number);
        if regime != FlowRegime::Incompressible {
            caveats.push(format!("Mach {:.2}: {}", result.mach_number, regime.caveat()));
        }

        AeroReport {
            cd: coeff.cd,
            cl: coeff.cl,
            cs: coeff.cs,
            cm: coeff.cmy,
            pressure_drag_fraction,
            drag_area: result.drag_area(),
            dynamic_pressure: result.tunnel.dynamic_pressure(),
            reference_area: result.tunnel.reference_area,
            reynolds_number: result.reynolds_number,
            mach_number: result.mach_number,
            converged: result.converged,
            iterations: result.flow.iterations,
            residual: result.flow.residual,
            y_plus_mean: result.surface.y_plus_mean,
            wake_peak_deficit: result.wake.peak_deficit(),
            caveats,
        }
    }

    /// The dimensional drag **force** in newtons — `F_D = (Cd·A)·q∞`, the
    /// drag area scaled by the free-stream dynamic pressure. This is the load
    /// the body actually feels, the number behind the dimensionless `cd`.
    pub fn drag_force(&self) -> f64 {
        self.drag_area * self.dynamic_pressure
    }

    /// The aerodynamic lift **force** in newtons — `F_L = Cl·A·q∞` (negative
    /// when the body makes downforce). The dimensional companion to `cl`.
    pub fn lift_force(&self) -> f64 {
        self.cl * self.reference_area * self.dynamic_pressure
    }

    /// The magnitude of the **total** aerodynamic force on the body (N) — the
    /// vector sum of lift, drag, and side force, `√(L² + D² + S²)`. This is the
    /// resultant load the structure must actually react; it always meets or
    /// exceeds the largest single component because the others add in
    /// quadrature.
    pub fn resultant_force(&self) -> f64 {
        let lift = self.lift_force();
        let drag = self.drag_force();
        let side = self.cs * self.reference_area * self.dynamic_pressure;
        (lift * lift + drag * drag + side * side).sqrt()
    }

    /// The lift-to-drag ratio `L/D = Cl / Cd` at this operating point — the
    /// single headline aerodynamic-efficiency number (in unpowered flight, the
    /// glide ratio: horizontal distance travelled per unit height lost). Returns
    /// `0` when the drag is non-positive.
    pub fn lift_to_drag(&self) -> f64 {
        if self.cd > 1e-9 {
            self.cl / self.cd
        } else {
            0.0
        }
    }

    /// The glide angle `γ = atan2(Cd, Cl)` (radians) — the descent slope in
    /// unpowered flight, where `tan γ = D/L = 1/(L/D)`. Defined for any sign of
    /// lift via `atan2`: a lifting body glides shallowly (small `γ`), a
    /// non-lifting or draggy body descends steeply (`γ → π/2` and beyond).
    pub fn glide_angle_rad(&self) -> f64 {
        self.cd.atan2(self.cl)
    }

    /// The Prandtl–Glauert compressibility factor at this run's Mach number —
    /// the subsonic correction that scales incompressible coefficients toward
    /// their compressible values (`Cl ≈ Cl_incompressible · factor`). See the
    /// free [`prandtl_glauert_factor`]; `1` at low speed, growing toward `M = 1`,
    /// `0` once sonic/supersonic (the linearised correction breaks down).
    pub fn prandtl_glauert_factor(&self) -> f64 {
        prandtl_glauert_factor(self.mach_number)
    }

    /// Render the report as a plain-text block — the form a CLI prints
    /// or an LLM relays.
    pub fn to_text(&self) -> String {
        let mut s = String::new();
        s.push_str("=== valenx-aero wind-tunnel report ===\n");
        s.push_str(&format!(
            "  Reynolds number : {:.3e}\n",
            self.reynolds_number
        ));
        s.push_str(&format!("  Mach number     : {:.3}\n", self.mach_number));
        s.push_str(&format!(
            "  P-G factor      : {:.3}  (1/\u{221A}(1\u{2212}M\u{00B2}))\n",
            self.prandtl_glauert_factor()
        ));
        s.push_str(&format!(
            "  converged       : {} ({} iterations, residual {:.2e})\n",
            self.converged, self.iterations, self.residual
        ));
        s.push_str("  --- coefficients ---\n");
        s.push_str(&format!("  drag      Cd : {:+.4}\n", self.cd));
        s.push_str(&format!("  lift      Cl : {:+.4}\n", self.cl));
        s.push_str(&format!("  side      Cs : {:+.4}\n", self.cs));
        s.push_str(&format!("  pitch Cm     : {:+.4}\n", self.cm));
        s.push_str(&format!("  lift/drag L/D: {:+.3}\n", self.lift_to_drag()));
        s.push_str(&format!(
            "  glide angle  : {:.2} deg\n",
            self.glide_angle_rad().to_degrees()
        ));
        s.push_str(&format!(
            "  drag area CdA: {:.4} m^2\n",
            self.drag_area
        ));
        s.push_str(&format!(
            "  dynamic press: {:.1} Pa  (q_inf)\n",
            self.dynamic_pressure
        ));
        s.push_str(&format!(
            "  drag force   : {:.1} N\n",
            self.drag_force()
        ));
        s.push_str(&format!(
            "  ref area A   : {:.4} m^2\n",
            self.reference_area
        ));
        s.push_str(&format!(
            "  lift force   : {:.1} N\n",
            self.lift_force()
        ));
        s.push_str(&format!(
            "  resultant F  : {:.1} N\n",
            self.resultant_force()
        ));
        s.push_str(&format!(
            "  pressure drag: {:.0}% of total ({}% friction)\n",
            100.0 * self.pressure_drag_fraction,
            (100.0 * (1.0 - self.pressure_drag_fraction)).round() as i64
        ));
        s.push_str("  --- flow ---\n");
        s.push_str(&format!("  mean y+         : {:.1}\n", self.y_plus_mean));
        s.push_str(&format!(
            "  wake peak deficit: {:.0}%\n",
            100.0 * self.wake_peak_deficit
        ));
        s.push_str("  --- caveats ---\n");
        for c in &self.caveats {
            s.push_str(&format!("  * {c}\n"));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::{run_windtunnel, AeroRequest};
    use crate::domain::TunnelSizing;
    use crate::geometry::box_body;
    use crate::turbulence::TurbulenceModel;
    use nalgebra::Vector3;

    /// A coarse but real grid for the report-formatting tests — they
    /// assert report *structure* (the headline numbers are copied
    /// through, the caveats are present, the text contains the labels),
    /// none of which needs a fine mesh. A coarse grid keeps the real
    /// end-to-end solve fast.
    fn coarse() -> TunnelSizing {
        TunnelSizing {
            cells_across_body: 4,
            max_cells: 40_000,
            ..TunnelSizing::default()
        }
    }

    #[test]
    fn report_summarises_a_run() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // The headline numbers are carried over verbatim.
        assert_eq!(report.cd, result.coefficients.cd);
        assert_eq!(report.cl, result.coefficients.cl);
        // The pressure-drag fraction is a fraction.
        assert!((0.0..=1.0).contains(&report.pressure_drag_fraction));
        // There is always at least the standing v1 caveat + the
        // turbulence note.
        assert!(report.caveats.len() >= 2);
        assert!(report
            .caveats
            .iter()
            .any(|c| c.contains("not ANSYS Fluent")));
    }

    #[test]
    fn report_text_contains_the_coefficients() {
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let req = AeroRequest::new(15.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(20);
        let result = run_windtunnel(&body, &req).unwrap();
        let text = AeroReport::from_result(&result).to_text();
        assert!(text.contains("wind-tunnel report"));
        assert!(text.contains("drag      Cd"));
        assert!(text.contains("Reynolds number"));
        assert!(text.contains("caveats"));
    }

    #[test]
    fn non_converged_run_gets_a_caveat() {
        // A 1-iteration run will not converge — the report must say so.
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_max_iterations(1);
        let result = run_windtunnel(&body, &req).unwrap();
        if !result.converged {
            let report = AeroReport::from_result(&result);
            assert!(
                report.caveats.iter().any(|c| c.contains("did NOT reach")),
                "a non-converged run should be flagged"
            );
        }
    }

    #[test]
    fn report_carries_dynamic_pressure_and_drag_force() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // q∞ is copied straight from the tunnel's free-stream conditions.
        assert_eq!(report.dynamic_pressure, result.tunnel.dynamic_pressure());
        assert!(report.dynamic_pressure > 0.0, "moving air has positive q");
        // Drag force is the drag area scaled by q: F_D = (Cd·A)·q∞.
        let expected = report.drag_area * report.dynamic_pressure;
        assert!((report.drag_force() - expected).abs() < 1e-9);
        // It also surfaces in the text dump.
        assert!(report.to_text().contains("drag force"));
    }

    #[test]
    fn report_carries_reference_area_and_lift_force() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // Reference area is the tunnel's, and drag_area = Cd·A stays consistent.
        assert_eq!(report.reference_area, result.tunnel.reference_area);
        assert!(report.reference_area > 0.0, "a real body has frontal area");
        assert!((report.drag_area - report.cd * report.reference_area).abs() < 1e-9);
        // Lift force is Cl·A·q (definitional), and it surfaces in the text dump.
        let expected = report.cl * report.reference_area * report.dynamic_pressure;
        assert!((report.lift_force() - expected).abs() < 1e-9);
        assert!(report.to_text().contains("lift force"));
    }

    #[test]
    fn resultant_force_is_the_quadrature_sum_of_lift_drag_side() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // Resultant = √(L² + D² + S²), the quadrature sum of the three forces.
        let side = report.cs * report.reference_area * report.dynamic_pressure;
        let expected =
            (report.lift_force().powi(2) + report.drag_force().powi(2) + side * side).sqrt();
        assert!((report.resultant_force() - expected).abs() < 1e-9);
        // It meets or exceeds every single component (quadrature, never smaller).
        assert!(report.resultant_force() >= report.drag_force().abs() - 1e-9);
        assert!(report.resultant_force() >= report.lift_force().abs() - 1e-9);
        // And it surfaces in the text dump.
        assert!(report.to_text().contains("resultant F"));
    }

    #[test]
    fn lift_to_drag_is_the_coefficient_ratio() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // L/D · Cd = Cl identically (the definitional ratio), and it is finite.
        assert!((report.lift_to_drag() * report.cd - report.cl).abs() < 1e-9);
        assert!(report.lift_to_drag().is_finite());
        // It surfaces in the text dump.
        assert!(report.to_text().contains("L/D"));
    }

    #[test]
    fn glide_angle_is_the_arctangent_of_drag_over_lift() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // γ = atan2(Cd, Cl): a body with positive drag descends at an angle
        // strictly between 0 and π.
        let g = report.glide_angle_rad();
        assert!(g > 0.0 && g < std::f64::consts::PI && g.is_finite(), "γ {g}");
        // For a lifting body, tan γ = D/L = 1/(L/D).
        if report.cl > 0.0 && report.lift_to_drag() > 0.0 {
            assert!((g.tan() - 1.0 / report.lift_to_drag()).abs() < 1e-9, "tan γ = 1/(L/D)");
        }
        // It surfaces in the text dump.
        assert!(report.to_text().contains("glide angle"));
    }

    #[test]
    fn prandtl_glauert_factor_matches_textbook_values() {
        // M = 0 → no compressibility correction.
        assert!((prandtl_glauert_factor(0.0) - 1.0).abs() < 1e-12);
        // M = 0.6 → 1/√(1−0.36) = 1/0.8 = 1.25.
        assert!((prandtl_glauert_factor(0.6) - 1.25).abs() < 1e-12);
        // M = 0.8 → 1/√(1−0.64) = 1/0.6 ≈ 1.6667.
        assert!((prandtl_glauert_factor(0.8) - 1.0 / 0.6).abs() < 1e-12);
        // It rises monotonically through the subsonic range.
        assert!(prandtl_glauert_factor(0.7) > prandtl_glauert_factor(0.3));
        // Sonic / supersonic / out-of-range → 0 (the correction breaks down).
        assert_eq!(prandtl_glauert_factor(1.0), 0.0);
        assert_eq!(prandtl_glauert_factor(1.5), 0.0);
        assert_eq!(prandtl_glauert_factor(-0.1), 0.0);
    }

    #[test]
    fn finite_wing_lift_slope_reduces_below_the_section_value() {
        use std::f64::consts::PI;
        let a0 = 2.0 * PI; // thin-airfoil section slope, per radian
        // AR = 6, e = 1: a = 2π / (1 + 2π/(π·6)) = 2π / (1 + 1/3) = 2π·0.75.
        let a = finite_wing_lift_slope(a0, 6.0, 1.0);
        assert!((a - 2.0 * PI * 0.75).abs() < 1e-9, "AR=6 slope {a}");
        // A finite wing is always gentler than its 2-D section.
        assert!(a < a0, "finite-wing slope must be < section slope");

        // As AR → ∞ the downwash vanishes and the 2-D slope is recovered.
        let a_inf = finite_wing_lift_slope(a0, 1.0e6, 1.0);
        assert!((a_inf - a0).abs() < 1e-3, "AR→∞ should recover a0, got {a_inf}");

        // Monotonic: a higher-aspect-ratio wing has a steeper slope.
        assert!(finite_wing_lift_slope(a0, 12.0, 1.0) > finite_wing_lift_slope(a0, 6.0, 1.0));

        // Non-physical inputs → 0 (the relation does not apply).
        assert_eq!(finite_wing_lift_slope(a0, 0.0, 1.0), 0.0);
        assert_eq!(finite_wing_lift_slope(a0, 6.0, 0.0), 0.0);
        assert_eq!(finite_wing_lift_slope(-1.0, 6.0, 1.0), 0.0);
        assert_eq!(finite_wing_lift_slope(f64::NAN, 6.0, 1.0), 0.0);
    }

    #[test]
    fn stall_speed_is_where_lift_equals_weight_at_max_lift() {
        // A light aircraft: 11 kN weight, 16 m² wing, sea-level air, C_Lmax = 1.5.
        let (w, s, rho, cl_max) = (11_000.0, 16.0, 1.225, 1.5);
        let v = stall_speed(w, s, rho, cl_max);
        assert!((v - 27.36).abs() < 0.1, "stall speed {v} m/s");
        // By construction the lift at V_stall and C_Lmax exactly balances the weight.
        let lift = 0.5 * rho * v * v * s * cl_max;
        assert!((lift - w).abs() < 1e-6, "L={lift} must equal W={w} at the stall");
        // V_stall ∝ 1/√C_Lmax: a higher max lift coefficient lowers the stall speed.
        let v_flapped = stall_speed(w, s, rho, 2.0 * cl_max);
        assert!((v_flapped - v / 2.0_f64.sqrt()).abs() < 1e-9, "∝ 1/√C_Lmax");
        // V_stall ∝ √W: four times the weight doubles the stall speed.
        let v_heavy = stall_speed(4.0 * w, s, rho, cl_max);
        assert!((v_heavy - 2.0 * v).abs() < 1e-9, "∝ √W");
        // Non-physical inputs → 0.
        assert_eq!(stall_speed(0.0, s, rho, cl_max), 0.0);
        assert_eq!(stall_speed(w, s, rho, 0.0), 0.0);
        assert_eq!(stall_speed(w, -1.0, rho, cl_max), 0.0);
        assert_eq!(stall_speed(f64::NAN, s, rho, cl_max), 0.0);
    }

    #[test]
    fn breguet_range_scales_with_speed_efficiency_and_fuel_fraction() {
        // R = (V/c)·(L/D)·ln(W₀/W₁). A long-haul jet: V = 250 m/s, c = 0.6/hr,
        // L/D = 17, W₀/W₁ = 1.5 → ≈ 1.03e7 m (~10,000 km).
        let v = 250.0_f64;
        let c = 0.6 / 3600.0; // 0.6 per hour → per second
        let ld = 17.0_f64;
        let ratio = 1.5_f64;
        let r = breguet_range(v, c, ld, ratio);
        let expected = v / c * ld * ratio.ln();
        assert!((r - expected).abs() < 1e-6, "range {r} vs {expected}");
        // A sensible long-haul figure (within a few thousand km of 10,000).
        assert!((6.0e6..1.4e7).contains(&r), "airliner range {:.0} km", r / 1e3);

        // Burning no fuel (W₀ = W₁) gives zero range.
        assert_eq!(breguet_range(v, c, ld, 1.0), 0.0);
        // R ∝ L/D, ∝ V, and ∝ ln(W₀/W₁) (ratio → ratio² doubles the log).
        assert!((breguet_range(v, c, 2.0 * ld, ratio) - 2.0 * r).abs() < 1e-6, "∝ L/D");
        assert!((breguet_range(2.0 * v, c, ld, ratio) - 2.0 * r).abs() < 1e-6, "∝ V");
        assert!((breguet_range(v, c, ld, ratio * ratio) - 2.0 * r).abs() < 1e-6, "∝ ln(W₀/W₁)");

        // Non-physical inputs → 0 (including an unphysical W₁ > W₀).
        assert_eq!(breguet_range(v, 0.0, ld, ratio), 0.0);
        assert_eq!(breguet_range(-1.0, c, ld, ratio), 0.0);
        assert_eq!(breguet_range(v, c, ld, 0.5), 0.0);
        assert_eq!(breguet_range(v, c, ld, f64::NAN), 0.0);
    }

    #[test]
    fn dynamic_pressure_ratio_is_half_gamma_mach_squared() {
        let g = 1.4;
        // q/p = ½γM²: worked values.
        assert!((dynamic_pressure_ratio(1.0, g) - 0.7).abs() < 1e-12, "M=1 → γ/2 = 0.7");
        assert!((dynamic_pressure_ratio(2.0, g) - 2.8).abs() < 1e-12, "M=2 → 2.8");
        assert!(dynamic_pressure_ratio(0.0, g).abs() < 1e-12, "M=0 → 0");

        // Quadratic in M (2× M → 4×) and linear in γ.
        assert!(
            (dynamic_pressure_ratio(2.0, g) - 4.0 * dynamic_pressure_ratio(1.0, g)).abs() < 1e-12,
            "quadratic in M"
        );
        assert!(
            (dynamic_pressure_ratio(2.0, 2.8) - 2.0 * dynamic_pressure_ratio(2.0, g)).abs() < 1e-12,
            "linear in γ"
        );

        // Incompressible Bernoulli limit: as M → 0, q/p ≈ p₀/p − 1 (threads
        // isentropic_stagnation_pressure_ratio).
        let m = 0.01;
        let q_over_p = dynamic_pressure_ratio(m, g);
        let stag_excess = isentropic_stagnation_pressure_ratio(m, g) - 1.0;
        assert!((q_over_p / stag_excess - 1.0).abs() < 1e-3, "q/p ≈ p₀/p − 1 as M→0");

        // Non-physical input → 0.
        assert_eq!(dynamic_pressure_ratio(-1.0, g), 0.0);
        assert_eq!(dynamic_pressure_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(dynamic_pressure_ratio(f64::NAN, g), 0.0);
    }

    #[test]
    fn mach_angle_is_the_mach_cone_half_angle() {
        use std::f64::consts::PI;
        // At M = 1 the Mach cone is a flat disc: μ = 90°.
        assert!((mach_angle(1.0) - PI / 2.0).abs() < 1e-12);
        // M = 2 → arcsin(0.5) = 30°; M = √2 → arcsin(1/√2) = 45°.
        assert!((mach_angle(2.0) - PI / 6.0).abs() < 1e-12, "M=2 → 30°");
        assert!((mach_angle(2.0_f64.sqrt()) - PI / 4.0).abs() < 1e-12, "M=√2 → 45°");
        // The cone narrows as Mach rises, and → 0 at hypersonic speed.
        assert!(mach_angle(5.0) < mach_angle(2.0), "narrows with Mach");
        assert!(mach_angle(1.0e6) < 1e-3, "→ 0 as M → ∞");
        // Subsonic (no Mach cone) and non-finite → 0.
        assert_eq!(mach_angle(0.8), 0.0);
        assert_eq!(mach_angle(0.0), 0.0);
        assert_eq!(mach_angle(f64::NAN), 0.0);
        assert_eq!(mach_angle(f64::INFINITY), 0.0);
    }

    #[test]
    fn prandtl_meyer_angle_matches_expansion_tables() {
        use std::f64::consts::PI;
        let g = 1.4;
        // ν(1) = 0: no turning at the sonic line.
        assert!(prandtl_meyer_angle(1.0, g).abs() < 1e-12, "ν(1) = 0");
        // Standard expansion-table points (γ = 1.4), in degrees.
        assert!((prandtl_meyer_angle(2.0, g).to_degrees() - 26.380).abs() < 1e-2, "ν(2) ≈ 26.38°");
        assert!((prandtl_meyer_angle(3.0, g).to_degrees() - 49.757).abs() < 1e-2, "ν(3) ≈ 49.76°");
        assert!((prandtl_meyer_angle(5.0, g).to_degrees() - 76.920).abs() < 1e-2, "ν(5) ≈ 76.92°");
        // Monotonically increasing in M for M > 1.
        let (a, b, c) = (
            prandtl_meyer_angle(1.5, g),
            prandtl_meyer_angle(2.5, g),
            prandtl_meyer_angle(4.0, g),
        );
        assert!(a < b && b < c, "ν rises with M: {a} {b} {c}");
        // Bounded above by the vacuum limit ν_max = (π/2)·(√((γ+1)/(γ−1)) − 1), which a
        // very large Mach number approaches from below (≈ 130.45° for air). The limit
        // is an independent closed form; the impl is the two-atan expression.
        let nu_max = 0.5 * PI * (((g + 1.0) / (g - 1.0)).sqrt() - 1.0);
        assert!((nu_max.to_degrees() - 130.454).abs() < 1e-2, "ν_max ≈ 130.45°");
        let nu_big = prandtl_meyer_angle(1.0e6, g);
        assert!(nu_big < nu_max && nu_big > 0.999 * nu_max, "ν(1e6) → ν_max⁻: {nu_big} vs {nu_max}");
        // Subsonic / non-physical → 0 (no expansion fan).
        assert_eq!(prandtl_meyer_angle(0.5, g), 0.0);
        assert_eq!(prandtl_meyer_angle(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(prandtl_meyer_angle(f64::NAN, g), 0.0);
        assert_eq!(prandtl_meyer_angle(2.0, f64::INFINITY), 0.0);
    }

    #[test]
    fn induced_drag_coefficient_is_the_lifting_line_drag_due_to_lift() {
        use std::f64::consts::PI;
        // No lift, no induced drag.
        assert_eq!(induced_drag_coefficient(0.0, 8.0, 0.8), 0.0);
        // Worked point: C_L=0.5, AR=8, e=0.8 → 0.25/(π·6.4) ≈ 0.0124339.
        let cdi = induced_drag_coefficient(0.5, 8.0, 0.8);
        assert!((cdi - 0.25 / (PI * 6.4)).abs() < 1e-12, "C_Di {cdi}");
        assert!((cdi - 0.012_433_9).abs() < 1e-6, "≈0.0124339, got {cdi}");
        // Elliptical loading (e=1) is the theoretical minimum C_Di = C_L²/(π·AR),
        // and it is strictly less than the e=0.8 case.
        let cdi_ell = induced_drag_coefficient(0.5, 8.0, 1.0);
        assert!((cdi_ell - 0.25 / (PI * 8.0)).abs() < 1e-12, "elliptical {cdi_ell}");
        assert!(cdi_ell < cdi, "elliptical loading has the least induced drag");
        // Quadratic in C_L: doubling the lift coefficient quadruples induced drag.
        let base = induced_drag_coefficient(0.4, 8.0, 0.8);
        assert!(
            (induced_drag_coefficient(0.8, 8.0, 0.8) - 4.0 * base).abs() < 1e-12,
            "∝ C_L²"
        );
        // Inverse in AR: doubling the aspect ratio halves the induced drag.
        assert!(
            (induced_drag_coefficient(0.5, 16.0, 0.8) - cdi / 2.0).abs() < 1e-12,
            "∝ 1/AR"
        );
        // Sign-blind: negative-lift downforce induces drag just the same (∝ C_L²).
        assert!(
            (induced_drag_coefficient(-0.5, 8.0, 0.8) - cdi).abs() < 1e-12,
            "C_L² is sign-blind"
        );
        // As AR → ∞ the induced drag vanishes (the downwash disappears).
        assert!(induced_drag_coefficient(0.5, 1.0e6, 1.0) < 1e-6, "AR→∞ → 0");
        // Non-physical inputs → 0 (the relation does not apply).
        assert_eq!(induced_drag_coefficient(0.5, 0.0, 0.8), 0.0);
        assert_eq!(induced_drag_coefficient(0.5, 8.0, 0.0), 0.0);
        assert_eq!(induced_drag_coefficient(f64::NAN, 8.0, 0.8), 0.0);
        assert_eq!(induced_drag_coefficient(0.5, f64::INFINITY, 0.8), 0.0);
    }

    #[test]
    fn max_lift_to_drag_ratio_is_the_best_glide_optimum() {
        use std::f64::consts::PI;
        let (cd0, ar, e) = (0.02, 8.0, 0.8);
        let ld_max = max_lift_to_drag_ratio(cd0, ar, e);
        // Closed form ½·√(π·e·AR/C_D0); worked point ≈ 15.85.
        assert!(
            (ld_max - 0.5 * (PI * e * ar / cd0).sqrt()).abs() < 1e-12,
            "closed form"
        );
        assert!((ld_max - 15.85).abs() < 0.05, "≈15.85, got {ld_max}");
        // The optimum is where induced drag equals parasite drag: at
        // C_L* = √(π·e·AR·C_D0), induced_drag_coefficient(C_L*) = C_D0, the total
        // drag is 2·C_D0 and (L/D)_max = C_L*/(2·C_D0).
        let cl_star = (PI * e * ar * cd0).sqrt();
        let cdi = induced_drag_coefficient(cl_star, ar, e);
        assert!((cdi - cd0).abs() < 1e-12, "induced = parasite at the optimum: {cdi}");
        assert!(
            (ld_max - cl_star / (2.0 * cd0)).abs() < 1e-9,
            "(L/D)_max = C_L*/(2·C_D0)"
        );
        // Scaling: ∝ √AR (4× the aspect ratio doubles it) and ∝ 1/√C_D0
        // (4× the parasite drag halves it).
        assert!(
            (max_lift_to_drag_ratio(cd0, 4.0 * ar, e) - 2.0 * ld_max).abs() < 1e-9,
            "∝ √AR"
        );
        assert!(
            (max_lift_to_drag_ratio(4.0 * cd0, ar, e) - ld_max / 2.0).abs() < 1e-9,
            "∝ 1/√C_D0"
        );
        // Non-physical inputs → 0.
        assert_eq!(max_lift_to_drag_ratio(0.0, ar, e), 0.0);
        assert_eq!(max_lift_to_drag_ratio(cd0, 0.0, e), 0.0);
        assert_eq!(max_lift_to_drag_ratio(cd0, ar, 0.0), 0.0);
        assert_eq!(max_lift_to_drag_ratio(f64::NAN, ar, e), 0.0);
    }

    #[test]
    fn isentropic_stagnation_pressure_ratio_matches_compressible_flow_tables() {
        // M = 0 → no compression, p0/p = 1.
        assert!((isentropic_stagnation_pressure_ratio(0.0, 1.4) - 1.0).abs() < 1e-12);
        // M = 1, γ = 1.4 → 1.2^3.5 ≈ 1.8929 (the sonic stagnation ratio for air).
        assert!((isentropic_stagnation_pressure_ratio(1.0, 1.4) - 1.2_f64.powf(3.5)).abs() < 1e-12);
        assert!((isentropic_stagnation_pressure_ratio(1.0, 1.4) - 1.8929).abs() < 1e-3, "sonic ≈ 1.893");
        // M = 2, γ = 1.4 → 1.8^3.5 ≈ 7.824.
        assert!((isentropic_stagnation_pressure_ratio(2.0, 1.4) - 1.8_f64.powf(3.5)).abs() < 1e-12);
        assert!((isentropic_stagnation_pressure_ratio(2.0, 1.4) - 7.824).abs() < 1e-2, "M=2 ≈ 7.82");
        // Monotonic increasing in Mach, and always ≥ 1.
        let (a, b, c) = (
            isentropic_stagnation_pressure_ratio(0.5, 1.4),
            isentropic_stagnation_pressure_ratio(1.5, 1.4),
            isentropic_stagnation_pressure_ratio(3.0, 1.4),
        );
        assert!(a >= 1.0 && a < b && b < c, "monotone ≥ 1: {a} {b} {c}");
        // Low-Mach limit reduces to 1 + (γ/2)·M² (the incompressible dynamic-pressure form).
        let m = 0.05;
        let exact = isentropic_stagnation_pressure_ratio(m, 1.4);
        assert!((exact - (1.0 + 0.5 * 1.4 * m * m)).abs() < 1e-4, "low-M ≈ 1+(γ/2)M²");
        // Non-physical input → the no-correction identity 1.0.
        assert_eq!(isentropic_stagnation_pressure_ratio(-0.5, 1.4), 1.0);
        assert_eq!(isentropic_stagnation_pressure_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(isentropic_stagnation_pressure_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(isentropic_stagnation_pressure_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn isentropic_stagnation_temperature_ratio_matches_compressible_flow_tables() {
        // M = 0 → no heating, T0/T = 1.
        assert!((isentropic_stagnation_temperature_ratio(0.0, 1.4) - 1.0).abs() < 1e-12);
        // M = 1, γ = 1.4 → 1.2; M = 2 → 1.8; M = 5 → 6.0 (hypersonic).
        assert!((isentropic_stagnation_temperature_ratio(1.0, 1.4) - 1.2).abs() < 1e-12);
        assert!((isentropic_stagnation_temperature_ratio(2.0, 1.4) - 1.8).abs() < 1e-12);
        assert!((isentropic_stagnation_temperature_ratio(5.0, 1.4) - 6.0).abs() < 1e-12);
        // Monotone increasing in Mach, and always ≥ 1.
        let (a, b, c) = (
            isentropic_stagnation_temperature_ratio(0.5, 1.4),
            isentropic_stagnation_temperature_ratio(1.5, 1.4),
            isentropic_stagnation_temperature_ratio(3.0, 1.4),
        );
        assert!(a >= 1.0 && a < b && b < c, "monotone ≥ 1: {a} {b} {c}");
        // Cross-check the isentropic identity p0/p = (T0/T)^(γ/(γ−1)) against #163.
        for m in [0.3_f64, 0.8, 1.0, 2.5, 4.0] {
            let t_ratio = isentropic_stagnation_temperature_ratio(m, 1.4);
            let p_from_t = t_ratio.powf(1.4 / 0.4);
            let p_ratio = isentropic_stagnation_pressure_ratio(m, 1.4);
            assert!(
                (p_from_t - p_ratio).abs() / p_ratio < 1e-12,
                "p0/p = (T0/T)^(γ/(γ−1)) at M={m}"
            );
        }
        // Non-physical input → the no-rise identity 1.0.
        assert_eq!(isentropic_stagnation_temperature_ratio(-0.5, 1.4), 1.0);
        assert_eq!(isentropic_stagnation_temperature_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(isentropic_stagnation_temperature_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(isentropic_stagnation_temperature_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn isentropic_stagnation_density_ratio_matches_compressible_flow_tables() {
        // M = 0 → no compression, ρ0/ρ = 1.
        assert!((isentropic_stagnation_density_ratio(0.0, 1.4) - 1.0).abs() < 1e-12);
        // M = 1, γ = 1.4 → 1.2^2.5 ≈ 1.5774 (the sonic stagnation density ratio for air).
        assert!((isentropic_stagnation_density_ratio(1.0, 1.4) - 1.2_f64.powf(2.5)).abs() < 1e-12);
        assert!((isentropic_stagnation_density_ratio(1.0, 1.4) - 1.5774).abs() < 1e-3, "sonic ≈ 1.577");
        // M = 2, γ = 1.4 → 1.8^2.5 ≈ 4.347.
        assert!((isentropic_stagnation_density_ratio(2.0, 1.4) - 1.8_f64.powf(2.5)).abs() < 1e-12);
        assert!((isentropic_stagnation_density_ratio(2.0, 1.4) - 4.347).abs() < 1e-2, "M=2 ≈ 4.35");
        // Monotone increasing in Mach, and always ≥ 1.
        let (a, b, c) = (
            isentropic_stagnation_density_ratio(0.5, 1.4),
            isentropic_stagnation_density_ratio(1.5, 1.4),
            isentropic_stagnation_density_ratio(3.0, 1.4),
        );
        assert!(a >= 1.0 && a < b && b < c, "monotone ≥ 1: {a} {b} {c}");
        // Cross-checks completing the isentropic trio (ties #163 + #169, non-tautological):
        //   ρ0/ρ = (T0/T)^(1/(γ−1)) = (p0/p)^(1/γ),  and  p0/p = (ρ0/ρ)^γ.
        for m in [0.3_f64, 0.8, 1.0, 2.5, 4.0] {
            let rho = isentropic_stagnation_density_ratio(m, 1.4);
            let t_ratio = isentropic_stagnation_temperature_ratio(m, 1.4);
            let p_ratio = isentropic_stagnation_pressure_ratio(m, 1.4);
            assert!((rho - t_ratio.powf(1.0 / 0.4)).abs() / rho < 1e-12, "ρ0/ρ=(T0/T)^(1/(γ−1)) at M={m}");
            assert!((rho - p_ratio.powf(1.0 / 1.4)).abs() / rho < 1e-12, "ρ0/ρ=(p0/p)^(1/γ) at M={m}");
            assert!((p_ratio - rho.powf(1.4)).abs() / p_ratio < 1e-12, "p0/p=(ρ0/ρ)^γ at M={m}");
        }
        // Non-physical input → the no-compression identity 1.0.
        assert_eq!(isentropic_stagnation_density_ratio(-0.5, 1.4), 1.0);
        assert_eq!(isentropic_stagnation_density_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(isentropic_stagnation_density_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(isentropic_stagnation_density_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn isentropic_area_ratio_matches_compressible_flow_tables() {
        // M = 1 is the sonic throat: A/A* = 1 exactly.
        assert!((isentropic_area_ratio(1.0, 1.4) - 1.0).abs() < 1e-12, "throat A/A* = 1");
        // Standard compressible-flow table points (γ = 1.4).
        assert!((isentropic_area_ratio(2.0, 1.4) - 1.6875).abs() < 1e-3, "M=2 → 1.6875");
        assert!((isentropic_area_ratio(3.0, 1.4) - 4.2346).abs() < 1e-3, "M=3 → 4.2346");
        assert!((isentropic_area_ratio(0.5, 1.4) - 1.3398).abs() < 1e-3, "M=0.5 → 1.3398");
        // A/A* has its MINIMUM (= 1) at the throat and rises on BOTH sides — the same
        // area ratio serves one subsonic and one supersonic solution.
        let throat = isentropic_area_ratio(1.0, 1.4);
        assert!(isentropic_area_ratio(0.9, 1.4) > throat, "subsonic side rises above the throat");
        assert!(isentropic_area_ratio(1.1, 1.4) > throat, "supersonic side rises above the throat");
        // Diverges as M → 0 (an infinite reservoir); non-physical → ∞.
        assert!(isentropic_area_ratio(0.0, 1.4).is_infinite(), "M→0 → ∞");
        assert!(isentropic_area_ratio(-1.0, 1.4).is_infinite(), "M<0 → ∞");
        assert!(isentropic_area_ratio(2.0, 1.0).is_infinite(), "γ≤1 → ∞");
        assert!(isentropic_area_ratio(f64::NAN, 1.4).is_infinite(), "non-finite M → ∞");
        // STRONG non-tautological cross-check via mass conservation ρAV = ρ*A*V*:
        //   A/A* = (ρ*/ρ)·(1/M)·√(T*/T) = [ρ0/ρ(M) / ρ0/ρ*(1)]·(1/M)·√((2/(γ+1))·T0/T),
        // composing the isentropic density ratio #175 and temperature ratio #169 — an
        // independent derivation (those use powf(1/(γ−1)) and a separate sqrt; the impl
        // is the single closed form in M).
        for m in [0.3_f64, 0.8, 1.0, 2.0, 3.5] {
            let g = 1.4;
            let rho = isentropic_stagnation_density_ratio(m, g);
            let rho_star = isentropic_stagnation_density_ratio(1.0, g);
            let t0_over_t = isentropic_stagnation_temperature_ratio(m, g);
            let expected = (rho / rho_star) * (1.0 / m) * (2.0 / (g + 1.0) * t0_over_t).sqrt();
            assert!(
                (isentropic_area_ratio(m, g) - expected).abs() / expected < 1e-9,
                "A/A* from mass conservation at M={m}"
            );
        }
    }

    #[test]
    fn mass_flow_function_peaks_at_choking_and_inverts_the_area_ratio() {
        let g = 1.4;
        // No flow at rest.
        assert!(mass_flow_function(0.0, g).abs() < 1e-12, "FF(0) = 0");
        // The choked peak FF(1) = √γ·(2/(γ+1))^((γ+1)/(2(γ−1))) ≈ 0.6847 for air,
        // re-derived here via the (2/(γ+1))^b form (the impl uses √γ·M·(…)^(−b)).
        let choke = g.sqrt() * (2.0 / (g + 1.0)).powf((g + 1.0) / (2.0 * (g - 1.0)));
        assert!((mass_flow_function(1.0, g) - choke).abs() < 1e-12, "FF(1) = choke constant");
        assert!((mass_flow_function(1.0, g) - 0.684731).abs() < 1e-5, "FF(1) ≈ 0.6847");
        // FF is MAXIMISED at M = 1 — the choking condition: a converging duct cannot
        // pass more than the sonic-throat flux. Both branches sit below the peak.
        let peak = mass_flow_function(1.0, g);
        assert!(mass_flow_function(0.5, g) < peak, "subsonic below the choke peak");
        assert!(mass_flow_function(2.0, g) < peak, "supersonic below the choke peak");
        assert!(
            mass_flow_function(0.9, g) < peak && mass_flow_function(1.1, g) < peak,
            "peak is at M = 1"
        );
        // STRONG non-tautological cross-check threading isentropic_area_ratio: by mass
        // conservation A/A* = FF(1)/FF(M). The area ratio is its own independent closed
        // form ((2/(γ+1)·T0/T)^b / M), so this ties two separately-derived relations.
        let ff_star = mass_flow_function(1.0, g);
        for m in [0.3_f64, 0.5, 2.0, 3.0] {
            let area_from_ff = ff_star / mass_flow_function(m, g);
            assert!(
                (isentropic_area_ratio(m, g) - area_from_ff).abs() / area_from_ff < 1e-9,
                "A/A* = FF(1)/FF(M) at M={m}"
            );
        }
        // Non-physical input → 0 (the no-flow sentinel).
        assert_eq!(mass_flow_function(-1.0, g), 0.0); // M < 0
        assert_eq!(mass_flow_function(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(mass_flow_function(f64::NAN, g), 0.0); // non-finite M
        assert_eq!(mass_flow_function(2.0, f64::INFINITY), 0.0); // non-finite γ
    }

    #[test]
    fn characteristic_mach_matches_the_sonic_referenced_speed() {
        let g = 1.4;
        // M* = 0 at rest; M* = 1 exactly at the sonic point (V = a = a*).
        assert!(characteristic_mach(0.0, g).abs() < 1e-12, "M*(0) = 0");
        assert!((characteristic_mach(1.0, g) - 1.0).abs() < 1e-12, "M*(1) = 1");
        // Worked point: M = 2 → M* = 2·√(2.4/3.6) ≈ 1.633.
        assert!((characteristic_mach(2.0, g) - 1.63299).abs() < 1e-4, "M*(2) ≈ 1.633");
        // M* labels subsonic/supersonic the same way M does (crosses 1 with M).
        assert!(
            characteristic_mach(0.5, g) < 1.0 && characteristic_mach(3.0, g) > 1.0,
            "subsonic < 1 < supersonic"
        );
        // Monotone in M, and SATURATES at the finite limit √((γ+1)/(γ−1)) (≈ 2.449 for
        // air) as M → ∞ — unlike the unbounded ordinary Mach number.
        assert!(characteristic_mach(1.5, g) < characteristic_mach(2.5, g), "monotone in M");
        let limit = ((g + 1.0) / (g - 1.0)).sqrt();
        assert!((limit - 6.0_f64.sqrt()).abs() < 1e-12, "limit = √6 for air");
        let m_big = characteristic_mach(1.0e6, g);
        assert!(m_big < limit && m_big > 0.999 * limit, "M*(1e6) → √6⁻: {m_big} vs {limit}");
        // STRONG cross-check — the PRANDTL relation across a normal shock M₁*·M₂* = 1,
        // with M₂ = normal_shock_downstream_mach(M₁) (#181). Ties #223 to #181: M* via
        // its closed form, M₂ via the Rankine–Hugoniot relation — different derivations.
        for &m1 in &[1.2_f64, 1.5, 2.0, 3.0, 5.0] {
            let m2 = normal_shock_downstream_mach(m1, g);
            let product = characteristic_mach(m1, g) * characteristic_mach(m2, g);
            assert!((product - 1.0).abs() < 1e-9, "M1*·M2* = 1 at M1={m1}: got {product}");
        }
        // Non-physical input → 0.
        assert_eq!(characteristic_mach(-1.0, g), 0.0);
        assert_eq!(characteristic_mach(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(characteristic_mach(f64::NAN, g), 0.0);
        assert_eq!(characteristic_mach(2.0, f64::INFINITY), 0.0);
    }

    #[test]
    fn normal_shock_downstream_mach_matches_compressible_flow_tables() {
        // M1 = 1 is the no-shock (infinitesimal) limit: M2 = 1 exactly.
        assert!((normal_shock_downstream_mach(1.0, 1.4) - 1.0).abs() < 1e-12);
        // Worked points against the normal-shock tables (γ = 1.4):
        // M1 = 2 → M2² = 1.8/5.4 = 1/3 → M2 ≈ 0.5774.
        assert!((normal_shock_downstream_mach(2.0, 1.4) - (1.0_f64 / 3.0).sqrt()).abs() < 1e-12);
        assert!((normal_shock_downstream_mach(2.0, 1.4) - 0.5774).abs() < 1e-3, "M1=2 → 0.5774");
        // M1 = 3 → M2 ≈ 0.4752.
        assert!((normal_shock_downstream_mach(3.0, 1.4) - 0.4752).abs() < 1e-3, "M1=3 → 0.4752");
        // A normal shock is always supersonic → subsonic: M2 < 1 for M1 > 1, and
        // the downstream Mach falls as the shock strengthens.
        let m_15 = normal_shock_downstream_mach(1.5, 1.4);
        let m_5 = normal_shock_downstream_mach(5.0, 1.4);
        assert!(m_15 < 1.0 && m_5 < 1.0 && m_5 < m_15, "M2 < 1 and falls with M1: {m_15} {m_5}");
        // Strong-shock limit M1 → ∞: M2 → √((γ−1)/2γ) ≈ 0.378, approached from above.
        let limit = (0.4_f64 / 2.8).sqrt();
        let m_big = normal_shock_downstream_mach(1.0e6, 1.4);
        assert!(
            (m_big - limit).abs() < 1e-3 && m_big > limit,
            "strong-shock limit ≈ {limit:.4}, got {m_big}"
        );
        // Subsonic/sonic upstream: no shock forms, the flow passes through unchanged.
        assert_eq!(normal_shock_downstream_mach(0.5, 1.4), 0.5);
        // Non-physical input → the sonic identity 1.0.
        assert_eq!(normal_shock_downstream_mach(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(normal_shock_downstream_mach(f64::NAN, 1.4), 1.0);
        assert_eq!(normal_shock_downstream_mach(-1.0, 1.4), 1.0);
        assert_eq!(normal_shock_downstream_mach(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn normal_shock_pressure_ratio_matches_compressible_flow_tables() {
        // M1 = 1 is the no-shock limit: no static-pressure jump.
        assert!((normal_shock_pressure_ratio(1.0, 1.4) - 1.0).abs() < 1e-12);
        // Worked table points (γ = 1.4): M1 = 2 → 4.5 (= 1 + (2.8/2.4)·3), M1 = 3 → 10.33.
        assert!((normal_shock_pressure_ratio(2.0, 1.4) - 4.5).abs() < 1e-12, "M1=2 → 4.5");
        assert!((normal_shock_pressure_ratio(3.0, 1.4) - 10.333).abs() < 1e-3, "M1=3 → 10.33");
        // A shock always compresses: p2/p1 > 1 and strictly increasing for M1 > 1.
        let (a, b, c) = (
            normal_shock_pressure_ratio(1.5, 1.4),
            normal_shock_pressure_ratio(2.5, 1.4),
            normal_shock_pressure_ratio(5.0, 1.4),
        );
        assert!(a > 1.0 && a < b && b < c, "compresses & rises with M1: {a} {b} {c}");
        // Grows without bound (∝ M1²) — unlike the density jump's finite ceiling.
        assert!(normal_shock_pressure_ratio(10.0, 1.4) > 100.0, "unbounded growth");
        // Subsonic/sonic upstream: no shock forms, the static pressure is unchanged.
        assert_eq!(normal_shock_pressure_ratio(0.5, 1.4), 1.0);
        // Non-physical input → the no-jump identity 1.0.
        assert_eq!(normal_shock_pressure_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(normal_shock_pressure_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(normal_shock_pressure_ratio(-1.0, 1.4), 1.0);
        assert_eq!(normal_shock_pressure_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn normal_shock_density_ratio_matches_compressible_flow_tables() {
        // M1 = 1 is the no-shock limit: no density jump.
        assert!((normal_shock_density_ratio(1.0, 1.4) - 1.0).abs() < 1e-12);
        // Worked table points (γ = 1.4): M1 = 2 → 9.6/3.6 = 8/3 ≈ 2.667, M1 = 3 → 3.857.
        assert!((normal_shock_density_ratio(2.0, 1.4) - 8.0 / 3.0).abs() < 1e-12, "M1=2 → 8/3");
        assert!((normal_shock_density_ratio(3.0, 1.4) - 3.857).abs() < 1e-3, "M1=3 → 3.857");
        // A shock always compresses: ρ2/ρ1 > 1 and strictly increasing for M1 > 1.
        let (a, b, c) = (
            normal_shock_density_ratio(1.5, 1.4),
            normal_shock_density_ratio(2.5, 1.4),
            normal_shock_density_ratio(5.0, 1.4),
        );
        assert!(a > 1.0 && a < b && b < c, "compresses & rises with M1: {a} {b} {c}");
        // Strong-shock ceiling M1 → ∞: ρ2/ρ1 → (γ+1)/(γ−1) = 6, approached from BELOW
        // (the finite limit that distinguishes it from the unbounded pressure jump).
        let limit = 2.4 / 0.4; // (γ+1)/(γ−1) = 6 for air
        let r_big = normal_shock_density_ratio(1.0e4, 1.4);
        assert!((r_big - limit).abs() < 1e-3 && r_big < limit, "strong-shock ceiling 6, got {r_big}");
        // Subsonic/sonic upstream: no shock forms, the density is unchanged.
        assert_eq!(normal_shock_density_ratio(0.5, 1.4), 1.0);
        // Non-physical input → the no-jump identity 1.0.
        assert_eq!(normal_shock_density_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(normal_shock_density_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(normal_shock_density_ratio(-1.0, 1.4), 1.0);
        assert_eq!(normal_shock_density_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn normal_shock_temperature_ratio_matches_rankine_hugoniot() {
        // M1 = 1 is the no-shock limit: no temperature jump.
        assert!((normal_shock_temperature_ratio(1.0, 1.4) - 1.0).abs() < 1e-12);
        // Worked table points (γ = 1.4): M1 = 2 → 1.6875 (= 4.5 / (8/3)), M1 = 3 → 2.679.
        assert!(
            (normal_shock_temperature_ratio(2.0, 1.4) - 1.6875).abs() < 1e-12,
            "M1=2 → 1.6875"
        );
        assert!(
            (normal_shock_temperature_ratio(3.0, 1.4) - 2.679).abs() < 1e-3,
            "M1=3 → 2.679"
        );
        // A shock always heats: T2/T1 > 1 and strictly increasing for M1 > 1.
        let (a, b, c) = (
            normal_shock_temperature_ratio(1.5, 1.4),
            normal_shock_temperature_ratio(2.5, 1.4),
            normal_shock_temperature_ratio(5.0, 1.4),
        );
        assert!(a > 1.0 && a < b && b < c, "heats & rises with M1: {a} {b} {c}");
        // Grows without bound (∝ M1²) — exceeds the density jump's finite ceiling (6).
        assert!(normal_shock_temperature_ratio(10.0, 1.4) > 8.0, "unbounded growth");
        // STRONG non-tautological cross-check: by the ideal-gas law T = p/(ρR), the
        // temperature jump is the quotient of the pressure and density jumps. This
        // impl uses the single closed form in M; the check composes the OTHER two
        // Rankine–Hugoniot relations (#187, #193) — different code paths.
        for &m in &[1.2_f64, 1.5, 2.0, 3.0, 5.0, 8.0] {
            let expected =
                normal_shock_pressure_ratio(m, 1.4) / normal_shock_density_ratio(m, 1.4);
            assert!(
                (normal_shock_temperature_ratio(m, 1.4) - expected).abs() < 1e-12,
                "T2/T1 = (p2/p1)/(ρ2/ρ1) at M={m}"
            );
        }
        // Holds for a different γ too (monatomic, γ = 5/3) — not air-specific.
        let g = 5.0 / 3.0;
        let expected = normal_shock_pressure_ratio(2.5, g) / normal_shock_density_ratio(2.5, g);
        assert!(
            (normal_shock_temperature_ratio(2.5, g) - expected).abs() < 1e-12,
            "γ=5/3 cross-check"
        );
        // Subsonic/sonic upstream: no shock forms, the static temperature is unchanged.
        assert_eq!(normal_shock_temperature_ratio(0.5, 1.4), 1.0);
        // Non-physical input → the no-jump identity 1.0.
        assert_eq!(normal_shock_temperature_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(normal_shock_temperature_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(normal_shock_temperature_ratio(-1.0, 1.4), 1.0);
        assert_eq!(normal_shock_temperature_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn normal_shock_stagnation_pressure_ratio_matches_compressible_flow_tables() {
        // M1 = 1 is the isentropic no-shock limit: total pressure conserved.
        assert!((normal_shock_stagnation_pressure_ratio(1.0, 1.4) - 1.0).abs() < 1e-12);
        // Worked table points (γ = 1.4): M1 = 2 → 0.7209, M1 = 3 → 0.3283.
        assert!(
            (normal_shock_stagnation_pressure_ratio(2.0, 1.4) - 0.7209).abs() < 1e-3,
            "M1=2 → 0.7209"
        );
        assert!(
            (normal_shock_stagnation_pressure_ratio(3.0, 1.4) - 0.3283).abs() < 1e-3,
            "M1=3 → 0.3283"
        );
        // The total pressure is always LOST and the loss deepens with M (strictly
        // DECREASING for M1 > 1) — the opposite sense to the rising static ratios.
        let (a, b, c) = (
            normal_shock_stagnation_pressure_ratio(1.5, 1.4),
            normal_shock_stagnation_pressure_ratio(2.5, 1.4),
            normal_shock_stagnation_pressure_ratio(5.0, 1.4),
        );
        assert!(a > b && b > c, "loss deepens with M1: {a} {b} {c}");
        // Strictly between 0 and 1 for any shock.
        assert!(a < 1.0 && c > 0.0, "0 < p02/p01 < 1");
        // STRONG non-tautological cross-check via the ENTROPY relation: the total-
        // pressure loss IS the entropy the shock generates, p02/p01 = exp(−Δs/R) with
        // Δs/R = (γ/(γ−1))·ln(T2/T1) − ln(p2/p1). The impl is the gas-dynamic closed
        // form; this is the thermodynamic identity composing #199 (T2/T1) and #187
        // (p2/p1) — an entirely different derivation.
        for &m in &[1.2_f64, 1.5, 2.0, 3.0, 5.0, 8.0] {
            let g = 1.4;
            let ds_over_r = g / (g - 1.0) * normal_shock_temperature_ratio(m, g).ln()
                - normal_shock_pressure_ratio(m, g).ln();
            assert!(
                (normal_shock_stagnation_pressure_ratio(m, g) - (-ds_over_r).exp()).abs() < 1e-9,
                "p02/p01 = exp(−Δs/R) at M={m}"
            );
        }
        // Holds for a different γ too (monatomic, γ = 5/3) — not air-specific.
        let g = 5.0 / 3.0;
        let ds_over_r = g / (g - 1.0) * normal_shock_temperature_ratio(2.5, g).ln()
            - normal_shock_pressure_ratio(2.5, g).ln();
        assert!(
            (normal_shock_stagnation_pressure_ratio(2.5, g) - (-ds_over_r).exp()).abs() < 1e-9,
            "γ=5/3 entropy cross-check"
        );
        // Subsonic/sonic upstream: no shock forms, total pressure conserved.
        assert_eq!(normal_shock_stagnation_pressure_ratio(0.5, 1.4), 1.0);
        // Non-physical input → the no-loss identity 1.0.
        assert_eq!(normal_shock_stagnation_pressure_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(normal_shock_stagnation_pressure_ratio(f64::NAN, 1.4), 1.0);
        assert_eq!(normal_shock_stagnation_pressure_ratio(-1.0, 1.4), 1.0);
        assert_eq!(normal_shock_stagnation_pressure_ratio(2.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn normal_shock_entropy_rise_matches_the_second_law_total_pressure_loss() {
        let g = 1.4;
        // M = 1: a vanishingly weak shock is reversible → no entropy generated.
        assert!(normal_shock_entropy_rise(1.0, g).abs() < 1e-12, "M=1 → Δs/R = 0");
        // Worked point M=2, γ=1.4: Δs/R = (γ/(γ−1))·ln(T2/T1) − ln(p2/p1) ≈ 0.3273.
        assert!((normal_shock_entropy_rise(2.0, g) - 0.3273).abs() < 1e-3, "M=2 → Δs/R ≈ 0.3273");
        // Monotone increasing with shock strength (a stronger shock generates more
        // entropy — the thermodynamic reason it recovers less total pressure).
        let (s15, s2, s3, s5) = (
            normal_shock_entropy_rise(1.5, g),
            normal_shock_entropy_rise(2.0, g),
            normal_shock_entropy_rise(3.0, g),
            normal_shock_entropy_rise(5.0, g),
        );
        assert!(0.0 < s15 && s15 < s2 && s2 < s3 && s3 < s5, "monotone: {s15} {s2} {s3} {s5}");
        // STRONG non-tautological cross-check via the second-law identity
        // Δs/R = −ln(p02/p01): the impl uses the static T/p jumps, this uses the
        // independently-derived stagnation-pressure recovery — two separate closed
        // forms (and a different γ too, monatomic 5/3, to show it is not air-specific).
        for &(m, gg) in &[(1.5_f64, 1.4_f64), (2.0, 1.4), (3.0, 1.4), (5.0, 1.4), (2.5, 5.0 / 3.0)] {
            let from_static = normal_shock_entropy_rise(m, gg);
            let from_stagnation = -normal_shock_stagnation_pressure_ratio(m, gg).ln();
            assert!(
                (from_static - from_stagnation).abs() / from_static < 1e-9,
                "Δs/R = −ln(p02/p01) at M={m}, γ={gg}: {from_static} vs {from_stagnation}"
            );
        }
        // Subsonic/sonic and non-physical → 0 (no shock, no entropy).
        assert_eq!(normal_shock_entropy_rise(0.5, g), 0.0); // subsonic
        assert_eq!(normal_shock_entropy_rise(-1.0, g), 0.0); // M < 0
        assert_eq!(normal_shock_entropy_rise(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(normal_shock_entropy_rise(f64::NAN, g), 0.0); // non-finite M
        assert_eq!(normal_shock_entropy_rise(2.0, f64::INFINITY), 0.0); // non-finite γ
    }

    #[test]
    fn rayleigh_pitot_ratio_matches_supersonic_pitot_tables() {
        let g = 1.4;
        // Subsonic: no bow shock → the pitot reads the isentropic total-to-static
        // ratio (continuous with the supersonic branch at M = 1).
        for m in [0.3_f64, 0.7, 1.0] {
            assert!(
                (rayleigh_pitot_ratio(m, g) - isentropic_stagnation_pressure_ratio(m, g)).abs()
                    < 1e-12,
                "subsonic pitot = isentropic at M={m}"
            );
        }
        // Standard supersonic-pitot table points (γ = 1.4): M=2 → 5.640, M=3 → 12.06.
        assert!((rayleigh_pitot_ratio(2.0, g) - 5.640).abs() < 1e-2, "M=2 → 5.640");
        assert!((rayleigh_pitot_ratio(3.0, g) - 12.061).abs() < 1e-2, "M=3 → 12.06");
        // Monotonically increasing with M (a pitot reads ever-higher overpressure).
        assert!(rayleigh_pitot_ratio(1.5, g) < rayleigh_pitot_ratio(2.5, g), "rises with M");
        // STRONG cross-check: p₀₂/p₁ = (p₀₂/p₀₁)·(p₀₁/p₁) — the across-shock total
        // recovery (#205) times the post-shock isentropic rise (#163), for several
        // M > 1. The impl uses the combined [(γ+1)M²/2]^… closed form; the check
        // composes the two independent stagnation-ratio fns (different expansions).
        for m in [1.2_f64, 1.5, 2.0, 3.0, 5.0] {
            let expected = normal_shock_stagnation_pressure_ratio(m, g)
                * isentropic_stagnation_pressure_ratio(m, g);
            assert!(
                (rayleigh_pitot_ratio(m, g) - expected).abs() / expected < 1e-9,
                "p02/p1 = (p02/p01)·(p01/p1) at M={m}"
            );
        }
        // Non-physical input → 1.0.
        assert_eq!(rayleigh_pitot_ratio(2.0, 1.0), 1.0); // γ ≤ 1
        assert_eq!(rayleigh_pitot_ratio(f64::NAN, g), 1.0);
        assert_eq!(rayleigh_pitot_ratio(-1.0, g), 1.0);
    }

    #[test]
    fn rayleigh_flow_total_temperature_ratio_matches_heat_addition_tables() {
        let g = 1.4;
        // Thermal-choking limit: T0/T0* = 1 exactly at M = 1 (the maximum).
        assert!((rayleigh_flow_total_temperature_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // No flow at rest.
        assert!(rayleigh_flow_total_temperature_ratio(0.0, g).abs() < 1e-12, "M=0 → 0");
        // Standard Rayleigh-flow table points (γ = 1.4).
        assert!((rayleigh_flow_total_temperature_ratio(2.0, g) - 0.7934).abs() < 1e-3, "M=2 → 0.7934");
        assert!((rayleigh_flow_total_temperature_ratio(0.5, g) - 0.6914).abs() < 1e-3, "M=0.5 → 0.6914");
        assert!((rayleigh_flow_total_temperature_ratio(3.0, g) - 0.6540).abs() < 1e-3, "M=3 → 0.6540");
        // The peak is at M = 1: both subsonic and supersonic branches sit below 1.
        let peak = rayleigh_flow_total_temperature_ratio(1.0, g);
        for m in [0.3_f64, 0.5, 0.9, 1.1, 2.0, 4.0] {
            assert!(rayleigh_flow_total_temperature_ratio(m, g) < peak, "T0/T0* < 1 at M={m}");
        }
        // STRONG non-tautological cross-check: the (1+(γ−1)/2·M²) factor IS the
        // isentropic T0/T, so T0/T0* = 2(γ+1)M²·[T0/T](M)/(1+γM²)² composing the
        // independently-implemented isentropic_stagnation_temperature_ratio.
        for m in [0.3_f64, 0.7, 1.0, 1.5, 2.0, 3.5] {
            let denom = 1.0 + g * m * m;
            let expected = 2.0 * (g + 1.0) * m * m
                * isentropic_stagnation_temperature_ratio(m, g)
                / (denom * denom);
            assert!(
                (rayleigh_flow_total_temperature_ratio(m, g) - expected).abs() / expected < 1e-12,
                "T0/T0* via isentropic T0/T at M={m}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(rayleigh_flow_total_temperature_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(rayleigh_flow_total_temperature_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(rayleigh_flow_total_temperature_ratio(f64::NAN, g), 0.0); // non-finite M
        assert_eq!(rayleigh_flow_total_temperature_ratio(2.0, f64::INFINITY), 0.0); // non-finite γ
    }

    #[test]
    fn rayleigh_flow_total_pressure_ratio_threads_the_pressure_ratio() {
        let g = 1.4;
        // Sonic reference: p₀/p₀* = 1 exactly at M = 1.
        assert!((rayleigh_flow_total_pressure_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // STRONG cross-check threading rayleigh_flow_pressure_ratio:
        // p₀/p₀* = (p/p*)·((2+(γ−1)M²)/(γ+1))^(γ/(γ−1)).
        for &(m, gam) in &[(0.3_f64, 1.4_f64), (0.5, 1.4), (2.0, 1.4), (4.0, 1.3), (0.8, 1.667)] {
            let expected = rayleigh_flow_pressure_ratio(m, gam)
                * ((2.0 + (gam - 1.0) * m * m) / (gam + 1.0)).powf(gam / (gam - 1.0));
            assert!(
                (rayleigh_flow_total_pressure_ratio(m, gam) - expected).abs() / expected < 1e-12,
                "p₀/p₀* = (p/p*)·bracket^(γ/(γ−1)) at M={m}, γ={gam}"
            );
        }
        // Standard Rayleigh-table values (γ = 1.4): M=0.5 → 1.1141, M=2 → 1.5031.
        assert!((rayleigh_flow_total_pressure_ratio(0.5, g) - 1.1141).abs() < 1e-3, "M=0.5 table");
        assert!((rayleigh_flow_total_pressure_ratio(2.0, g) - 1.5031).abs() < 1e-3, "M=2 table");
        // p₀/p₀* ≥ 1 with the minimum at the sonic point (heat addition erodes p₀).
        assert!(
            rayleigh_flow_total_pressure_ratio(0.5, g) > 1.0
                && rayleigh_flow_total_pressure_ratio(2.0, g) > 1.0,
            "p₀/p₀* > 1 away from the choke"
        );
        // Non-physical input → 0.
        assert_eq!(rayleigh_flow_total_pressure_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(rayleigh_flow_total_pressure_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(rayleigh_flow_total_pressure_ratio(f64::NAN, g), 0.0); // non-finite M
    }

    #[test]
    fn rayleigh_flow_velocity_ratio_threads_the_pressure_ratio() {
        let g = 1.4;
        // Sonic reference: V/V* = 1 exactly at M = 1.
        assert!((rayleigh_flow_velocity_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // STRONG cross-check threading rayleigh_flow_pressure_ratio: V/V* = M²·(p/p*).
        for &(m, gam) in &[(0.3_f64, 1.4_f64), (0.5, 1.4), (2.0, 1.4), (4.0, 1.3), (0.8, 1.667)] {
            let expected = m * m * rayleigh_flow_pressure_ratio(m, gam);
            assert!(
                (rayleigh_flow_velocity_ratio(m, gam) - expected).abs() / expected < 1e-12,
                "V/V* = M²·(p/p*) at M={m}, γ={gam}"
            );
        }
        // Standard Rayleigh-table values (γ = 1.4): M=0.5 → 0.44444, M=2 → 1.45454.
        assert!((rayleigh_flow_velocity_ratio(0.5, g) - 0.44444).abs() < 1e-3, "M=0.5 table value");
        assert!((rayleigh_flow_velocity_ratio(2.0, g) - 1.45454).abs() < 1e-3, "M=2 table value");
        // Monotone through the sonic point: subsonic < 1 < supersonic.
        assert!(
            rayleigh_flow_velocity_ratio(0.5, g) < 1.0 && rayleigh_flow_velocity_ratio(2.0, g) > 1.0,
            "subsonic V/V* < 1 < supersonic"
        );
        // Non-physical input → 0.
        assert_eq!(rayleigh_flow_velocity_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(rayleigh_flow_velocity_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(rayleigh_flow_velocity_ratio(f64::NAN, g), 0.0); // non-finite M
    }

    #[test]
    fn rayleigh_flow_pressure_ratio_threads_the_temperature_ratio() {
        let g = 1.4;
        // Sonic reference: p/p* = 1 exactly at M = 1.
        assert!((rayleigh_flow_pressure_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // STRONG cross-check threading rayleigh_flow_temperature_ratio: T/T* = M²·(p/p*)².
        for &(m, gam) in &[(0.3_f64, 1.4_f64), (0.5, 1.4), (2.0, 1.4), (4.0, 1.3), (0.8, 1.667)] {
            let expected = m * m * rayleigh_flow_pressure_ratio(m, gam).powi(2);
            assert!(
                (rayleigh_flow_temperature_ratio(m, gam) - expected).abs() / expected < 1e-12,
                "T/T* = M²·(p/p*)² at M={m}, γ={gam}"
            );
        }
        // Standard Rayleigh-table values (γ = 1.4): M=0.5 → 1.7778, M=2 → 0.36364.
        assert!((rayleigh_flow_pressure_ratio(0.5, g) - 1.7778).abs() < 1e-3, "M=0.5 table value");
        assert!((rayleigh_flow_pressure_ratio(2.0, g) - 0.36364).abs() < 1e-3, "M=2 table value");
        // Subsonic p/p* > 1, supersonic p/p* < 1; tends to the max (1+γ) as M → 0.
        assert!(
            rayleigh_flow_pressure_ratio(0.5, g) > 1.0 && rayleigh_flow_pressure_ratio(2.0, g) < 1.0,
            "subsonic > 1 > supersonic"
        );
        assert!((rayleigh_flow_pressure_ratio(1.0e-6, g) - 2.4).abs() < 1e-3, "M→0 limit 1+γ");
        // Non-physical input → 0.
        assert_eq!(rayleigh_flow_pressure_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(rayleigh_flow_pressure_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(rayleigh_flow_pressure_ratio(f64::NAN, g), 0.0); // non-finite M
    }

    #[test]
    fn rayleigh_flow_temperature_ratio_peaks_below_sonic() {
        let g = 1.4;
        // Sonic reference: T/T* = 1 at M = 1; no flow at rest.
        assert!((rayleigh_flow_temperature_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        assert!(rayleigh_flow_temperature_ratio(0.0, g).abs() < 1e-12, "M=0 → 0");
        // Standard Rayleigh-flow table points (γ = 1.4).
        assert!((rayleigh_flow_temperature_ratio(2.0, g) - 0.5289).abs() < 1e-3, "M=2 → 0.5289");
        assert!((rayleigh_flow_temperature_ratio(0.5, g) - 0.7901).abs() < 1e-3, "M=0.5 → 0.7901");
        // The static-temperature MAX is at M = 1/√γ (≈ 0.845), BELOW the sonic point,
        // and there T/T* > 1 (hotter than at sonic) — the classic Rayleigh feature.
        let m_peak = 1.0 / g.sqrt();
        let t_peak = rayleigh_flow_temperature_ratio(m_peak, g);
        assert!(t_peak > 1.0, "T/T* at M=1/√γ exceeds 1, got {t_peak}");
        for &m in &[0.6_f64, 0.75, 0.95, 1.0, 1.2] {
            assert!(t_peak >= rayleigh_flow_temperature_ratio(m, g), "peak at M=1/√γ, beaten at M={m}");
        }
        // STRONG cross-check threading rayleigh_flow_total_temperature_ratio #247 AND
        // isentropic_stagnation_temperature_ratio: since T0/T0* = (T0/T)·(T/T*)/(T0*/T*)
        // with T0*/T* = (γ+1)/2, we have T/T* = (T0/T0*)·((γ+1)/2)/(T0/T).
        for &m in &[0.3_f64, 0.7, 1.0, 1.5, 2.0, 3.5] {
            let expected = rayleigh_flow_total_temperature_ratio(m, g) * ((g + 1.0) / 2.0)
                / isentropic_stagnation_temperature_ratio(m, g);
            assert!(
                (rayleigh_flow_temperature_ratio(m, g) - expected).abs() / expected < 1e-9,
                "T/T* via T0/T0* and isentropic at M={m}"
            );
        }
        // Non-physical input → 0.
        assert_eq!(rayleigh_flow_temperature_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(rayleigh_flow_temperature_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(rayleigh_flow_temperature_ratio(f64::NAN, g), 0.0); // non-finite M
        assert_eq!(rayleigh_flow_temperature_ratio(2.0, f64::INFINITY), 0.0); // non-finite γ
    }

    #[test]
    fn fanno_friction_parameter_is_the_choking_length() {
        let g = 1.4;
        // At the sonic choke (M = 1) no further friction length remains.
        assert!(fanno_friction_parameter(1.0, g).abs() < 1e-12, "M=1 → 0");
        // STRONG cross-check threading fanno_flow_temperature_ratio (#259): the ln
        // argument (γ+1)M²/(2+(γ−1)M²) is exactly M²·(T/T*).
        for &(m, gam) in &[(0.3_f64, 1.4_f64), (0.5, 1.4), (2.0, 1.4), (4.0, 1.3), (0.8, 1.667)] {
            let expected = (1.0 - m * m) / (gam * m * m)
                + (gam + 1.0) / (2.0 * gam) * (m * m * fanno_flow_temperature_ratio(m, gam)).ln();
            assert!(
                (fanno_friction_parameter(m, gam) - expected).abs() / expected.abs() < 1e-12,
                "4fL*/D via M²·(T/T*) at M={m}, γ={gam}"
            );
        }
        // Standard Fanno-table values (γ = 1.4): M=0.5 → 1.0691, M=2 → 0.3050.
        assert!((fanno_friction_parameter(0.5, g) - 1.0691).abs() < 1e-3, "M=0.5 table value");
        assert!((fanno_friction_parameter(2.0, g) - 0.3050).abs() < 1e-3, "M=2 table value");
        // Positive on both sides of the choke (the minimum is 0 at M = 1).
        assert!(
            fanno_friction_parameter(0.5, g) > 0.0 && fanno_friction_parameter(2.0, g) > 0.0,
            "4fL*/D ≥ 0, positive away from the choke"
        );
        // Non-physical input → 0.
        assert_eq!(fanno_friction_parameter(0.0, g), 0.0); // M = 0
        assert_eq!(fanno_friction_parameter(-1.0, g), 0.0); // M < 0
        assert_eq!(fanno_friction_parameter(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(fanno_friction_parameter(f64::NAN, g), 0.0); // non-finite M
    }

    #[test]
    fn fanno_flow_velocity_ratio_threads_the_fanno_family() {
        let g = 1.4;
        // Sonic reference: V/V* = 1 exactly at M = 1.
        assert!((fanno_flow_velocity_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // STRONG dual cross-check threading the Fanno temperature (#259) and pressure
        // (#265) ratios: V/V* = M·√(T/T*) = M²·(p/p*) — two distinct closed forms.
        for &(m, gam) in &[(0.3_f64, 1.4_f64), (0.5, 1.4), (2.0, 1.4), (4.0, 1.3), (0.8, 1.667)] {
            let v = fanno_flow_velocity_ratio(m, gam);
            let via_t = m * fanno_flow_temperature_ratio(m, gam).sqrt();
            let via_p = m * m * fanno_flow_pressure_ratio(m, gam);
            assert!((v - via_t).abs() / via_t < 1e-12, "V/V* = M·√(T/T*) at M={m}, γ={gam}");
            assert!((v - via_p).abs() / via_p < 1e-12, "V/V* = M²·(p/p*) at M={m}, γ={gam}");
        }
        // Subsonic V/V* < 1, supersonic > 1 — the opposite trend to T/T* and p/p*.
        assert!(fanno_flow_velocity_ratio(0.5, g) < 1.0, "subsonic V/V* < 1");
        assert!(fanno_flow_velocity_ratio(2.0, g) > 1.0, "supersonic V/V* > 1");
        // Monotonic increasing in M.
        assert!(
            fanno_flow_velocity_ratio(1.5, g) > fanno_flow_velocity_ratio(0.5, g),
            "V/V* increases with M"
        );
        // Bounded high-Mach limit: V/V* → √((γ+1)/(γ−1)) (= √6 for γ = 1.4).
        let v_inf = ((g + 1.0) / (g - 1.0)).sqrt();
        assert!(
            (fanno_flow_velocity_ratio(1.0e4, g) - v_inf).abs() / v_inf < 1e-6,
            "V/V* → √((γ+1)/(γ−1))"
        );
        // No flow at rest; non-physical input → 0.
        assert_eq!(fanno_flow_velocity_ratio(0.0, g), 0.0); // M = 0
        assert_eq!(fanno_flow_velocity_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(fanno_flow_velocity_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(fanno_flow_velocity_ratio(f64::NAN, g), 0.0); // non-finite M
    }

    #[test]
    fn fanno_flow_pressure_ratio_threads_the_temperature_ratio() {
        let g = 1.4;
        // Sonic reference: p/p* = 1 exactly at M = 1.
        assert!((fanno_flow_pressure_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // STRONG cross-check threading fanno_flow_temperature_ratio (#259): the Fanno
        // pressure ratio is √(T/T*)/M — a distinct closed form.
        for &(m, gam) in &[(0.3_f64, 1.4_f64), (0.5, 1.4), (2.0, 1.4), (4.0, 1.3), (0.8, 1.667)] {
            let expected = fanno_flow_temperature_ratio(m, gam).sqrt() / m;
            assert!(
                (fanno_flow_pressure_ratio(m, gam) - expected).abs() / expected < 1e-12,
                "p/p* = √(T/T*)/M at M={m}, γ={gam}"
            );
        }
        // Subsonic pressure exceeds sonic; supersonic falls below it.
        assert!(fanno_flow_pressure_ratio(0.5, g) > 1.0, "subsonic p/p* > 1");
        assert!(fanno_flow_pressure_ratio(2.0, g) < 1.0, "supersonic p/p* < 1");
        // Monotonic decreasing in M (friction drives the flow toward sonic).
        assert!(
            fanno_flow_pressure_ratio(0.5, g) > fanno_flow_pressure_ratio(1.5, g),
            "p/p* decreases with M"
        );
        // Diverges at rest: p/p* → ∞ as M → 0.
        assert!(fanno_flow_pressure_ratio(0.0, g).is_infinite(), "M=0 → ∞");
        // Non-physical input → 0.
        assert_eq!(fanno_flow_pressure_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(fanno_flow_pressure_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(fanno_flow_pressure_ratio(f64::NAN, g), 0.0); // non-finite M
    }

    #[test]
    fn fanno_flow_temperature_ratio_is_the_friction_dual_of_rayleigh() {
        let g = 1.4;
        // Sonic reference: T/T* = 1 exactly at M = 1.
        assert!((fanno_flow_temperature_ratio(1.0, g) - 1.0).abs() < 1e-12, "M=1 → 1");
        // Hottest static temperature is at rest: T/T* = (γ+1)/2 at M = 0.
        assert!(
            (fanno_flow_temperature_ratio(0.0, g) - (g + 1.0) / 2.0).abs() < 1e-12,
            "M=0 → (γ+1)/2"
        );
        // STRONG cross-check threading isentropic_stagnation_temperature_ratio: Fanno
        // flow is ADIABATIC ⇒ T0 constant ⇒ T/T* = (T0/T*)/(T0/T) =
        // isentropic(1, γ)/isentropic(M, γ). Two different closed forms (a single
        // rational vs a ratio of two affine-in-M² terms) — non-tautological.
        for &gam in &[1.4_f64, 1.3, 1.667] {
            for &m in &[0.0_f64, 0.3, 0.5, 0.845, 1.0, 2.0, 4.0] {
                let expected = isentropic_stagnation_temperature_ratio(1.0, gam)
                    / isentropic_stagnation_temperature_ratio(m, gam);
                assert!(
                    (fanno_flow_temperature_ratio(m, gam) - expected).abs() < 1e-12,
                    "T/T* = T0T(1)/T0T(M) at M={m}, γ={gam}"
                );
            }
        }
        // Monotonic decreasing in M — friction always drives the flow toward sonic.
        let mut prev = f64::INFINITY;
        for &m in &[0.0_f64, 0.5, 1.0, 2.0, 5.0, 20.0] {
            let t = fanno_flow_temperature_ratio(m, g);
            assert!(t < prev, "T/T* decreases with M at M={m} ({t} !< {prev})");
            prev = t;
        }
        // Vanishes in the high-Mach limit.
        assert!(fanno_flow_temperature_ratio(1.0e4, g) < 1e-6, "large M → 0");
        // Non-physical input → 0.
        assert_eq!(fanno_flow_temperature_ratio(2.0, 1.0), 0.0); // γ ≤ 1
        assert_eq!(fanno_flow_temperature_ratio(-1.0, g), 0.0); // M < 0
        assert_eq!(fanno_flow_temperature_ratio(f64::NAN, g), 0.0); // non-finite M
        assert_eq!(fanno_flow_temperature_ratio(2.0, f64::INFINITY), 0.0); // non-finite γ
    }
}
