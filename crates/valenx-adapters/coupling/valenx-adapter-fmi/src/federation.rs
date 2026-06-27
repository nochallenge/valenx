//! HELICS-style co-simulation **federation** with distributed time
//! coordination — a pure-Rust, in-process reimplementation.
//!
//! ## What this adds over [`crate::cosim`]
//!
//! [`crate::cosim::CoSimMaster`] advances a graph of [`crate::cosim::Subsystem`]s
//! on a single, globally fixed macro-step. That is exactly right for a tightly
//! coupled, equal-rate split. It is *not* right when the parts of a system run
//! at genuinely different rates, or when one part may only influence another
//! after a known minimum delay: forcing them onto one global step either wastes
//! work (stepping the slow part too often) or loses accuracy (stepping the fast
//! part too coarsely).
//!
//! A **federation** generalizes the fixed-step master into *dependency-aware
//! coordinated advancement*. Each [`Federate`] keeps its own current time and
//! its own [`TimePolicy`] (period, lookahead, time-delta). A federate *requests*
//! the next time it would like to be granted; the [`Broker`] *grants* the
//! largest time that is safe for everyone, given each federate's stated minimum
//! influence delay (its **lookahead**). Repeatedly requesting and granting walks
//! the whole federation forward in dependency order, each federate landing only
//! on the times it actually needs to act on.
//!
//! ## Algorithm source / attribution
//!
//! The time-coordination algorithm here — the **time request / grant** model
//! where each federate advertises a *lookahead* (a.k.a. output/impact delay)
//! and a grant is bounded by the minimum over the federation of
//! `granted_time + lookahead` — is reimplemented in pure Rust from the
//! *published* design of **HELICS** (the Hierarchical Engine for Large-scale
//! Infrastructure Co-Simulation), which is distributed under the **BSD-3-Clause**
//! license by the U.S. Department of Energy / LLNL / PNNL / NREL. See the HELICS
//! documentation on "Timing" / "Time Coordination" and the published papers
//! (Palmintier et al., *Design of the HELICS High-Performance Transmission–
//! Distribution–Communication–Market Co-Simulation Framework*, 2017).
//!
//! **This is a clean-room reimplementation of the published algorithm, not a
//! port.** No HELICS C++ source was copied; there is no ZeroMQ, no networking,
//! and no shared C ABI — everything below is in-process safe Rust over the
//! existing [`crate::cosim`] types. The names mirror HELICS concepts
//! (federate / broker / publication / subscription / endpoint / lookahead) so
//! that someone who knows HELICS can read this, but the code is original.
//!
//! ## Honest scope
//!
//! * **In-process only.** "The federation" is a set of [`Federate`]s owned by
//!   one [`Broker`] in one process. There is no transport, no discovery, no
//!   distributed clock — "distributed" here means *dependency-distributed time
//!   coordination*, the algorithm, not multi-machine networking.
//! * **Value exchange is by name** (`publish` / `subscribe`),
//!   decoupled from stepping, exactly like HELICS value federates.
//! * **Message endpoints** carry discrete, timestamped messages between named
//!   endpoints (the HELICS message-federate idea), again in-process.
//! * A co-simulation **FMU** (anything implementing [`crate::cosim::Subsystem`])
//!   can act as a federate via [`SubsystemFederate`], so the FMI importer and
//!   the native subsystems plug straight in. Where a full FMU wire would be
//!   heavy, the coordination is proven with a [`mock`] closure-driven federate.
//!
//! ## Fail-loud
//!
//! Unknown federate / publication / subscription / endpoint names, a negative
//! period or lookahead, a NaN time, or a non-monotone time request are all an
//! [`Err`] with a precise [`FederationError`] — never a silent default.

use std::collections::HashMap;

use crate::cosim::Subsystem;
use crate::error::{FmiError, Result};

// ---------------------------------------------------------------------------
// Time
// ---------------------------------------------------------------------------

/// Simulation time, as a plain `f64` in the federation's time units
/// (seconds, conventionally). Kept as a transparent newtype so the time
/// arithmetic (rounding a request up to the next period, bounding a grant by
/// a dependency's `time + lookahead`) lives in one place and the API reads in
/// HELICS terms.
///
/// Ordinary `f64` comparison is used; the benchmark schedules are chosen so
/// the periods and lookaheads are exactly representable (halves, quarters,
/// integers), making the grant traces bit-exact and deterministic.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Time(pub f64);

impl Time {
    /// Time zero (the federation start).
    pub const ZERO: Time = Time(0.0);

    /// The maximum finite time — the value an unbounded grant collapses to
    /// before the federation horizon clamps it. (Used internally as the
    /// identity for a `min` over dependency bounds.)
    pub const MAX: Time = Time(f64::MAX);

    /// The underlying `f64`.
    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }

    /// `self + dt`, used to form a dependency bound `granted + lookahead`.
    #[inline]
    fn add(self, dt: f64) -> Time {
        Time(self.0 + dt)
    }

    /// The smaller of two times.
    #[inline]
    fn min(self, other: Time) -> Time {
        if other.0 < self.0 {
            other
        } else {
            self
        }
    }

    /// The larger of two times.
    #[inline]
    fn max(self, other: Time) -> Time {
        if other.0 > self.0 {
            other
        } else {
            self
        }
    }
}

// ---------------------------------------------------------------------------
// Time policy
// ---------------------------------------------------------------------------

/// A federate's timing contract — its half of the HELICS time-coordination
/// bargain.
///
/// * `period` — if non-zero, every granted time is snapped **up** to the next
///   integer multiple of `period` (offset by `period_offset`). A periodic
///   federate that requests `0.3` with `period = 0.5` is granted `0.5`. A zero
///   period means "grant exactly what is safe" (an event-driven federate).
/// * `lookahead` — the federate's promise that *nothing it does at time `t`
///   can influence any other federate before `t + lookahead`*. This is the
///   single most important quantity for coordination: a federate with a
///   positive lookahead lets its dependents be granted further ahead (up to
///   `granted + lookahead`) without risking a missed interaction. HELICS calls
///   this the impact/output delay; it must be `>= 0`.
/// * `time_delta` — the minimum amount of time a single grant must advance the
///   federate by (its smallest sensible step). A grant is never issued closer
///   than `time_delta` past the federate's current granted time. Zero means no
///   minimum.
///
/// Construct with [`TimePolicy::periodic`], [`TimePolicy::event_driven`], or
/// the builder setters, all of which validate (a negative `period`,
/// `lookahead`, or `time_delta`, or a non-finite value, is a fail-loud error).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TimePolicy {
    period: f64,
    period_offset: f64,
    lookahead: f64,
    time_delta: f64,
}

impl Default for TimePolicy {
    /// An event-driven federate: no period, no lookahead, no minimum step.
    /// It is granted exactly the time the federation deems safe.
    fn default() -> Self {
        TimePolicy {
            period: 0.0,
            period_offset: 0.0,
            lookahead: 0.0,
            time_delta: 0.0,
        }
    }
}

impl TimePolicy {
    /// A purely periodic federate with the given `period` (and no lookahead
    /// or minimum step). Every grant is snapped up to a multiple of `period`.
    ///
    /// Fail-loud: a non-positive or non-finite `period` is rejected.
    pub fn periodic(period: f64) -> Result<Self> {
        let p = TimePolicy {
            period,
            ..Default::default()
        };
        p.validate()?;
        if period <= 0.0 {
            return Err(err_policy(format!(
                "periodic() needs a strictly positive period, got {period}"
            )));
        }
        Ok(p)
    }

    /// An event-driven federate (the [`Default`]): granted exactly the safe
    /// time, with no period snapping.
    pub fn event_driven() -> Self {
        TimePolicy::default()
    }

    /// Set the period (0 disables period snapping). Returns `self` for
    /// chaining; validates fail-loud.
    pub fn with_period(mut self, period: f64) -> Result<Self> {
        self.period = period;
        self.validate()?;
        Ok(self)
    }

    /// Set the offset applied before period snapping (grants land on
    /// `offset + k*period`). Validates fail-loud (must be finite, `>= 0`).
    pub fn with_period_offset(mut self, offset: f64) -> Result<Self> {
        self.period_offset = offset;
        self.validate()?;
        Ok(self)
    }

    /// Set the lookahead / output delay (must be `>= 0`). Returns `self` for
    /// chaining; validates fail-loud.
    pub fn with_lookahead(mut self, lookahead: f64) -> Result<Self> {
        self.lookahead = lookahead;
        self.validate()?;
        Ok(self)
    }

    /// Set the minimum grant advance (`time_delta`, must be `>= 0`). Returns
    /// `self` for chaining; validates fail-loud.
    pub fn with_time_delta(mut self, time_delta: f64) -> Result<Self> {
        self.time_delta = time_delta;
        self.validate()?;
        Ok(self)
    }

    /// The configured period (0 = no snapping).
    pub fn period(&self) -> f64 {
        self.period
    }

    /// The configured lookahead / output delay.
    pub fn lookahead(&self) -> f64 {
        self.lookahead
    }

    /// The configured minimum grant advance.
    pub fn time_delta(&self) -> f64 {
        self.time_delta
    }

    /// Whether this federate is **interruptible** — i.e. it may be granted a
    /// time *earlier* than the one it requested when the safety frontier sits
    /// below its target (HELICS "may be interrupted by an event"). A pure
    /// event-driven federate (no period) is interruptible; a periodic federate
    /// is *not* — it is only ever granted on its period grid, and simply waits
    /// until the frontier reaches the next grid point.
    fn is_interruptible(&self) -> bool {
        self.period <= 0.0
    }

    /// Reject negative or non-finite parameters.
    fn validate(&self) -> Result<()> {
        for (val, what) in [
            (self.period, "period"),
            (self.period_offset, "period_offset"),
            (self.lookahead, "lookahead"),
            (self.time_delta, "time_delta"),
        ] {
            if !val.is_finite() {
                return Err(err_policy(format!("{what} must be finite, got {val}")));
            }
            if val < 0.0 {
                return Err(err_policy(format!(
                    "{what} must be non-negative, got {val}"
                )));
            }
        }
        Ok(())
    }

    /// Snap `t` **up** to this policy's period grid (`offset + k*period`),
    /// `k` an integer, returning the smallest grid point `>= t`. With no
    /// period this is `t` unchanged.
    ///
    /// This is the HELICS "next valid time" rule: a periodic federate is only
    /// granted times that lie on its period. A tiny epsilon guards against a
    /// `t` that is a floating-point hair below an exact grid point being
    /// pushed a whole period forward.
    fn snap_up(&self, t: Time) -> Time {
        if self.period <= 0.0 {
            return t;
        }
        let rel = t.0 - self.period_offset;
        if rel <= 0.0 {
            return Time(self.period_offset);
        }
        // Number of whole periods at or below `rel`.
        let k = (rel / self.period).floor();
        let grid = self.period_offset + k * self.period;
        // Relative tolerance: treat `t` within 1e-9*period of a grid point as
        // already on the grid (don't bump it a full period forward).
        let tol = self.period * 1e-9;
        if t.0 <= grid + tol {
            Time(grid)
        } else {
            Time(self.period_offset + (k + 1.0) * self.period)
        }
    }
}

// ---------------------------------------------------------------------------
// Value exchange: publications and subscriptions (pub/sub by name)
// ---------------------------------------------------------------------------

/// A typed value carried over a named publication. The federation's value bus is
/// untyped at the wire level (everything is one of these); a subscriber asks
/// for the value it expects. This mirrors HELICS value types without pulling
/// in a serialization dependency.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A real scalar (the common case for physical coupling signals).
    Double(f64),
    /// An integer.
    Int(i64),
    /// A boolean flag.
    Bool(bool),
    /// A free-form string payload.
    Text(String),
    /// A vector of reals (e.g. a state slice).
    Vector(Vec<f64>),
}

impl Value {
    /// Read this value as a `f64`, or a fail-loud type-mismatch error.
    pub fn as_double(&self) -> Result<f64> {
        match self {
            Value::Double(v) => Ok(*v),
            other => Err(err_type("Double", other)),
        }
    }

    /// Read this value as an `i64`, or a fail-loud type-mismatch error.
    pub fn as_int(&self) -> Result<i64> {
        match self {
            Value::Int(v) => Ok(*v),
            other => Err(err_type("Int", other)),
        }
    }

    /// Read this value as a `bool`, or a fail-loud type-mismatch error.
    pub fn as_bool(&self) -> Result<bool> {
        match self {
            Value::Bool(v) => Ok(*v),
            other => Err(err_type("Bool", other)),
        }
    }

    /// Borrow this value as a string slice, or a fail-loud type-mismatch.
    pub fn as_text(&self) -> Result<&str> {
        match self {
            Value::Text(v) => Ok(v.as_str()),
            other => Err(err_type("Text", other)),
        }
    }

    /// Borrow this value as a real vector, or a fail-loud type-mismatch.
    pub fn as_vector(&self) -> Result<&[f64]> {
        match self {
            Value::Vector(v) => Ok(v.as_slice()),
            other => Err(err_type("Vector", other)),
        }
    }

    fn type_name(&self) -> &'static str {
        match self {
            Value::Double(_) => "Double",
            Value::Int(_) => "Int",
            Value::Bool(_) => "Bool",
            Value::Text(_) => "Text",
            Value::Vector(_) => "Vector",
        }
    }
}

/// A discrete, timestamped message between named endpoints.
///
/// Unlike a [`Value`] (which holds *the latest* value and is read at any time),
/// a `Message` is delivered *once*, at the granted time at or after its
/// `time`, modelling HELICS message-federate traffic (packets, events).
#[derive(Clone, Debug, PartialEq)]
pub struct Message {
    /// Source endpoint name.
    pub source: String,
    /// Destination endpoint name.
    pub destination: String,
    /// The time at which this message becomes deliverable.
    pub time: Time,
    /// Opaque payload bytes.
    pub data: Vec<u8>,
}

// ---------------------------------------------------------------------------
// Federate
// ---------------------------------------------------------------------------

/// The behaviour a federate runs when the broker grants it a time.
///
/// On each grant the federate may read its subscribed values, advance its own
/// internal state up to the granted time, and write its publications / send
/// messages. The `ctx` gives it scoped access to the value and message buses
/// (read subscriptions, write publications, send messages); `granted` is the
/// time it has been advanced to.
///
/// Returning `Err` aborts the federation run fail-loud (e.g. an internal
/// model blew up). The default native subsystem driver
/// ([`SubsystemFederate`]) and the [`mock`] closure driver both implement this.
pub trait FederateBehavior {
    /// Advance to `granted`, reading inputs and writing outputs via `ctx`.
    fn on_grant(&mut self, granted: Time, ctx: &mut GrantContext<'_>) -> Result<()>;
}

/// A single member of the federation: a name, its [`TimePolicy`], its current
/// granted time, its pending request, and the behaviour it runs on a grant.
///
/// A federate is registered with the [`Broker`] via [`Broker::add_federate`];
/// its publications, subscriptions, and endpoints are then declared by name on
/// the broker. Direct field access is deliberately not exposed — coordination
/// state is owned and mutated only by the broker's grant loop.
pub struct Federate {
    name: String,
    policy: TimePolicy,
    /// Time this federate has actually been granted (starts at 0; only
    /// meaningful once `ever_granted`).
    granted: Time,
    /// Whether this federate has been granted at least once. Until then its
    /// *target* (next desired grant) is the federation start (time 0), which is
    /// where every federate receives its initial grant.
    ever_granted: bool,
    /// The time this federate last *requested* (its desired next grant).
    /// `None` until it has made a request.
    requested: Option<Time>,
    /// Once a federate has finished (requested at/after the horizon, or
    /// disconnected), it no longer constrains others and is no longer granted.
    finished: bool,
    behavior: Box<dyn FederateBehavior>,
}

impl Federate {
    /// The federate's unique name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The time this federate has been granted so far.
    pub fn granted_time(&self) -> Time {
        self.granted
    }

    /// This federate's time policy.
    pub fn policy(&self) -> &TimePolicy {
        &self.policy
    }

    /// Whether this federate has disconnected (and so no longer constrains
    /// the federation).
    pub fn is_finished(&self) -> bool {
        self.finished
    }

    /// The federate's **target time**: the next time it intends to act,
    /// i.e. its next desired grant.
    ///
    /// * Before its first grant, the target is the federation start (time 0):
    ///   every federate takes an initial grant at 0.
    /// * After a grant, the target is its explicit pending request if it made
    ///   one; otherwise the next natural tick implied by its policy (the next
    ///   period grid point, or `granted + time_delta`, or `granted` for a pure
    ///   event-driven federate with neither). In all cases the result is
    ///   snapped up to the period grid and forced at least `time_delta` past
    ///   the current granted time.
    fn target_time(&self) -> Time {
        let policy = &self.policy;
        if !self.ever_granted {
            // Initial grant lands on the period grid at/after time 0.
            return policy.snap_up(Time::ZERO);
        }
        let base = match self.requested {
            Some(req) => req,
            None => {
                // No explicit request: advance by the smallest natural step.
                if policy.period > 0.0 {
                    // The next grid point strictly after the current grant.
                    policy.snap_up(self.granted.add(policy.period))
                } else if policy.time_delta > 0.0 {
                    self.granted.add(policy.time_delta)
                } else {
                    // Pure event-driven with no minimum step: it has nothing
                    // more to ask for, so its target is where it already is —
                    // it will only be re-granted if the safety frontier moves,
                    // and otherwise contributes no further grants.
                    self.granted
                }
            }
        };
        let min_next = if policy.time_delta > 0.0 {
            self.granted.add(policy.time_delta)
        } else {
            self.granted
        };
        policy.snap_up(base.max(min_next))
    }

    /// The dependency bound this federate imposes on the rest of the
    /// federation: **no other federate may be safely granted past this
    /// federate's next output time plus its lookahead.**
    ///
    /// A federate's outputs are fixed from its last grant until it next acts
    /// (its [`Federate::target_time`]); its effect cannot reach anyone before
    /// `target + lookahead`. This is the HELICS guarantee that lets a fast
    /// federate run ahead of a slow one only as far as the slow one's
    /// lookahead permits. A finished federate imposes no bound ([`Time::MAX`]).
    fn dependency_bound(&self) -> Time {
        if self.finished {
            Time::MAX
        } else {
            self.target_time().add(self.policy.lookahead)
        }
    }
}

// ---------------------------------------------------------------------------
// Broker
// ---------------------------------------------------------------------------

/// The federation coordinator: owns the federates, the named value bus, and
/// the named message bus, and runs the time request/grant algorithm.
///
/// Lifecycle:
/// 1. [`Broker::new`].
/// 2. [`Broker::add_federate`] for each member (returns a [`FederateId`]).
/// 3. [`Broker::register_publication`] / [`Broker::register_subscription`] /
///    [`Broker::register_endpoint`] to declare the named interfaces.
/// 4. [`Broker::run_until`] (or [`Broker::step`] for one grant) to advance the
///    whole federation, optionally collecting a [`GrantRecord`] trace.
///
/// All lookups are by name and fail-loud on an unknown name.
pub struct Broker {
    federates: Vec<Federate>,
    /// Map federate name -> index, for fail-loud name lookups.
    index_of: HashMap<String, FederateId>,

    /// Publication key -> (owning federate, latest value, time it was set).
    publications: HashMap<String, PublicationSlot>,
    /// Subscription key -> the publication key it tracks, per subscriber.
    subscriptions: HashMap<String, SubscriptionSlot>,
    /// Endpoint name -> owning federate.
    endpoints: HashMap<String, FederateId>,
    /// Undelivered messages, kept sorted-enough by delivery time on drain.
    message_queue: Vec<Message>,

    /// Federation horizon: no grant is issued past this time.
    horizon: Time,
}

/// A handle to a registered federate (its index in the broker).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FederateId(pub usize);

/// Internal: a publication's owner and its latest value.
struct PublicationSlot {
    owner: FederateId,
    /// The latest published value, with the time it was published.
    value: Option<(Value, Time)>,
}

/// Internal: a subscription's owner and the publication key it reads.
struct SubscriptionSlot {
    owner: FederateId,
    publication_key: String,
}

/// One row of a federation grant trace: which federate was granted, and the
/// time it was granted to. [`Broker::run_until`] returns these in order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrantRecord {
    /// The federate that was granted.
    pub federate: FederateId,
    /// The time it was granted to.
    pub time: Time,
}

impl Broker {
    /// A new, empty broker with the given federation `horizon` (no grant is
    /// ever issued past this time). Fail-loud on a non-finite or negative
    /// horizon.
    pub fn new(horizon: f64) -> Result<Self> {
        if !horizon.is_finite() || horizon < 0.0 {
            return Err(err_policy(format!(
                "federation horizon must be finite and non-negative, got {horizon}"
            )));
        }
        Ok(Broker {
            federates: Vec::new(),
            index_of: HashMap::new(),
            publications: HashMap::new(),
            subscriptions: HashMap::new(),
            endpoints: HashMap::new(),
            message_queue: Vec::new(),
            horizon: Time(horizon),
        })
    }

    /// Register a federate with a unique `name`, a [`TimePolicy`], and its
    /// [`FederateBehavior`]. Returns its [`FederateId`].
    ///
    /// Fail-loud: a duplicate name is rejected (names address the federate in
    /// every later call).
    pub fn add_federate(
        &mut self,
        name: impl Into<String>,
        policy: TimePolicy,
        behavior: Box<dyn FederateBehavior>,
    ) -> Result<FederateId> {
        policy.validate()?;
        let name = name.into();
        if self.index_of.contains_key(&name) {
            return Err(err_dup_federate(&name));
        }
        let id = FederateId(self.federates.len());
        self.index_of.insert(name.clone(), id);
        self.federates.push(Federate {
            name,
            policy,
            granted: Time::ZERO,
            ever_granted: false,
            requested: None,
            finished: false,
            behavior,
        });
        Ok(id)
    }

    /// Look up a federate id by name, fail-loud on an unknown name.
    pub fn federate_id(&self, name: &str) -> Result<FederateId> {
        self.index_of
            .get(name)
            .copied()
            .ok_or_else(|| err_unknown_federate(name))
    }

    /// Borrow a federate by id, fail-loud on an out-of-range id.
    pub fn federate(&self, id: FederateId) -> Result<&Federate> {
        self.federates
            .get(id.0)
            .ok_or_else(|| err_unknown_federate(&format!("#{}", id.0)))
    }

    /// Number of registered federates.
    pub fn len(&self) -> usize {
        self.federates.len()
    }

    /// Whether no federates are registered.
    pub fn is_empty(&self) -> bool {
        self.federates.is_empty()
    }

    /// The federation horizon (no grant is issued past this time).
    pub fn horizon(&self) -> Time {
        self.horizon
    }

    /// Register a named publication owned by federate `owner`.
    ///
    /// Fail-loud: unknown owner, or a duplicate publication key.
    pub fn register_publication(
        &mut self,
        owner: FederateId,
        key: impl Into<String>,
    ) -> Result<()> {
        self.check_federate(owner)?;
        let key = key.into();
        if self.publications.contains_key(&key) {
            return Err(err_dup_interface("publication", &key));
        }
        self.publications
            .insert(key, PublicationSlot { owner, value: None });
        Ok(())
    }

    /// Register a named subscription owned by `owner` that reads the
    /// publication `publication_key`.
    ///
    /// The publication need not exist yet (HELICS allows declaring a
    /// subscription to a not-yet-registered key), but reading it before any
    /// federate publishes to it yields `None` — fail-loud only when the
    /// *subscription key itself* is duplicated or the owner is unknown.
    pub fn register_subscription(
        &mut self,
        owner: FederateId,
        sub_key: impl Into<String>,
        publication_key: impl Into<String>,
    ) -> Result<()> {
        self.check_federate(owner)?;
        let sub_key = sub_key.into();
        if self.subscriptions.contains_key(&sub_key) {
            return Err(err_dup_interface("subscription", &sub_key));
        }
        self.subscriptions.insert(
            sub_key,
            SubscriptionSlot {
                owner,
                publication_key: publication_key.into(),
            },
        );
        Ok(())
    }

    /// Register a named message endpoint owned by `owner`.
    ///
    /// Fail-loud: unknown owner, or a duplicate endpoint name.
    pub fn register_endpoint(&mut self, owner: FederateId, name: impl Into<String>) -> Result<()> {
        self.check_federate(owner)?;
        let name = name.into();
        if self.endpoints.contains_key(&name) {
            return Err(err_dup_interface("endpoint", &name));
        }
        self.endpoints.insert(name, owner);
        Ok(())
    }

    /// Read the latest value on publication `key`, with the time it was
    /// published — or `None` if nothing has been published yet. Fail-loud on
    /// an unknown publication key.
    pub fn publication_value(&self, key: &str) -> Result<Option<(Value, Time)>> {
        self.publications
            .get(key)
            .map(|slot| slot.value.clone())
            .ok_or_else(|| err_unknown_interface("publication", key))
    }

    // ---- the time-coordination core ------------------------------------

    /// Advance the federation by exactly one grant: pick the federate whose
    /// safe grant time is earliest, advance it, and let it run.
    ///
    /// This is one turn of the HELICS-style loop:
    ///
    /// 1. Every not-yet-finished federate has a *desired* next time — the time
    ///    it last requested, snapped up to its period grid, and at least
    ///    `time_delta` past where it already is. (A federate that has never
    ///    requested defaults to "as soon as possible", i.e. its own minimum
    ///    next step.)
    /// 2. The **granted** time for a federate is its desired time clamped by
    ///    the federation's safety bound: the minimum over *all other*
    ///    federates of `their granted_time + their lookahead`. A federate may
    ///    never be granted past the earliest point another federate could
    ///    still affect it.
    /// 3. The broker grants the federate whose clamped time is smallest
    ///    (ties broken by federate index for determinism), advances it,
    ///    delivers any due messages, and runs its behaviour. The behaviour
    ///    sets the federate's *next* request.
    ///
    /// Returns the [`GrantRecord`] of the federate advanced, or `Ok(None)`
    /// when every federate has finished or reached the horizon (the run is
    /// complete).
    pub fn step(&mut self) -> Result<Option<GrantRecord>> {
        // Compute each active federate's desired next time and the global
        // safety bound, then choose who to grant.
        let chosen = self.select_next()?;
        let Some((id, grant_time)) = chosen else {
            return Ok(None);
        };

        // Advance the chosen federate's clock, then run its behaviour. The
        // behaviour reads subscriptions / writes publications / sends + reads
        // messages through a scoped context, and finally sets its next
        // request (via the context) — or marks itself finished.
        self.federates[id.0].granted = grant_time;
        self.federates[id.0].ever_granted = true;

        // Drain messages now deliverable to this federate's endpoints.
        let inbox = self.drain_messages_for(id, grant_time);

        // Run behaviour with a scoped context borrowing the buses.
        let mut next_request: Option<Time> = None;
        let mut finish = false;
        {
            // Temporarily move the behaviour out so we can borrow `self`'s
            // buses mutably while the behaviour runs. (A behaviour never needs
            // to touch its own `Federate` struct — only the buses.)
            let mut behavior =
                std::mem::replace(&mut self.federates[id.0].behavior, Box::new(NullBehavior));
            let mut ctx = GrantContext {
                broker: self,
                me: id,
                granted: grant_time,
                inbox,
                next_request: &mut next_request,
                finish: &mut finish,
            };
            let result = behavior.on_grant(grant_time, &mut ctx);
            // Restore the behaviour regardless of outcome.
            self.federates[id.0].behavior = behavior;
            result?;
        }

        // Apply the federate's decision about its next time. `requested` is
        // strictly "what this federate asked for during THIS grant": it is
        // reset every grant so that a federate which says nothing falls back to
        // its natural period/delta tick (see `Federate::target_time`), and a
        // federate that wants a specific next time must re-request each grant.
        self.federates[id.0].requested = None;
        if finish {
            self.federates[id.0].finished = true;
        } else if let Some(req) = next_request {
            if !req.0.is_finite() {
                return Err(err_bad_request(self.federates[id.0].name(), req));
            }
            if req < grant_time {
                return Err(err_nonmonotone_request(
                    self.federates[id.0].name(),
                    req,
                    grant_time,
                ));
            }
            if req >= self.horizon {
                // Requesting at/after the horizon means "done".
                self.federates[id.0].finished = true;
            } else {
                self.federates[id.0].requested = Some(req);
            }
        }

        Ok(Some(GrantRecord {
            federate: id,
            time: grant_time,
        }))
    }

    /// Run the federation to completion (or until `max_grants` grants have
    /// been issued, a safety valve against a misconfigured federate that never
    /// advances). Returns the full ordered grant trace.
    ///
    /// "Completion" is when every federate has finished or every federate has
    /// reached the horizon. The trace is exactly the sequence of
    /// [`GrantRecord`]s [`Broker::step`] produced, which the benchmark tests
    /// pin to an exact expected schedule.
    pub fn run_until(&mut self, max_grants: usize) -> Result<Vec<GrantRecord>> {
        let mut trace = Vec::new();
        for _ in 0..max_grants {
            match self.step()? {
                Some(rec) => trace.push(rec),
                None => return Ok(trace),
            }
        }
        Err(err_grant_limit(max_grants))
    }

    /// Choose the next `(federate, granted_time)` to advance, or `None` when
    /// the federation is complete. This is the heart of the coordination.
    ///
    /// For every still-active federate `i`:
    ///
    /// * its **target** `T_i` is [`Federate::target_time`] (its next desired
    ///   grant — time 0 for an ungranted federate, else its request / next
    ///   period tick);
    /// * the **safety frontier** seen by `i` is the minimum over every *other*
    ///   federate `j` of [`Federate::dependency_bound`] = `T_j + lookahead_j`
    ///   (the earliest time `j` could next emit a value that reaches `i`);
    /// * the **grantable time** `G_i` depends on whether `i` is interruptible
    ///   (see [`TimePolicy::is_interruptible`]). A **periodic** federate is
    ///   *gated*: it is granted at exactly its grid `target`, and only once
    ///   `frontier >= target` (otherwise it waits). An **event-driven**
    ///   federate may be woken early, so `G_i = min(T_i, frontier, horizon)`.
    ///   Either way `G_i` is clamped not to fall below `i`'s current grant.
    ///
    /// A federate is **eligible** only if it has never been granted (it still
    /// owes its initial grant at time 0) or its grantable time strictly
    /// exceeds where it already is (`G_i > granted_i`). This is what prevents a
    /// federate that is *done at the current time* from being re-granted there:
    /// when several zero-lookahead federates sit at the same time, each takes
    /// exactly one grant there and then waits for the frontier to move.
    ///
    /// Among eligible federates the broker grants the one with the smallest
    /// `G_i`; ties break by lowest federate index (so at a shared time the
    /// lower-indexed federate publishes before a higher-indexed one reads,
    /// giving Gauss-Seidel-style same-time value visibility). `None` means no
    /// federate is eligible — the federation is complete.
    fn select_next(&self) -> Result<Option<(FederateId, Time)>> {
        let mut best: Option<(FederateId, Time)> = None;
        for (i, f) in self.federates.iter().enumerate() {
            if f.finished || f.granted >= self.horizon {
                continue;
            }
            let id = FederateId(i);

            // Safety frontier: min over OTHER federates of their dependency
            // bound. With no other federates, only the horizon limits us.
            let mut frontier = Time::MAX;
            let mut had_other = false;
            for (k, g) in self.federates.iter().enumerate() {
                if k == i {
                    continue;
                }
                had_other = true;
                frontier = frontier.min(g.dependency_bound());
            }
            if !had_other {
                frontier = self.horizon;
            }

            // Grantable time. A periodic federate is *gated*: it is granted
            // only on its grid `target`, and only once the frontier has
            // reached it (`frontier >= target`); otherwise it must wait. An
            // interruptible (event-driven) federate may instead be woken early
            // at the frontier, so it is grantable at `min(target, frontier)`.
            let target = f.target_time();
            let grant = if f.policy.is_interruptible() {
                target.min(frontier).min(self.horizon).max(f.granted)
            } else {
                // Periodic: must wait until the frontier reaches the target.
                // The initial grant (ungranted, target == 0) is never gated
                // out because frontier >= 0 always.
                if frontier < target {
                    continue;
                }
                target.min(self.horizon).max(f.granted)
            };

            // Eligibility: ungranted federates owe their initial grant; a
            // granted federate is only re-grantable if the frontier has moved
            // it strictly forward.
            let eligible = !f.ever_granted || grant > f.granted;
            if !eligible {
                continue;
            }

            // Earliest grant wins; ties -> lowest federate index.
            match best {
                None => best = Some((id, grant)),
                Some((_, t)) if grant < t => best = Some((id, grant)),
                _ => {}
            }
        }

        Ok(best)
    }

    /// Remove and return the messages now deliverable (time `<= now`) to any
    /// endpoint owned by federate `id`, in delivery-time order.
    fn drain_messages_for(&mut self, id: FederateId, now: Time) -> Vec<Message> {
        // Which endpoint names does this federate own?
        let mine: Vec<&String> = self
            .endpoints
            .iter()
            .filter(|(_, owner)| **owner == id)
            .map(|(name, _)| name)
            .collect();
        if mine.is_empty() {
            return Vec::new();
        }
        let owned: Vec<String> = mine.into_iter().cloned().collect();

        let mut due = Vec::new();
        let mut kept = Vec::with_capacity(self.message_queue.len());
        for m in self.message_queue.drain(..) {
            if m.time <= now && owned.contains(&m.destination) {
                due.push(m);
            } else {
                kept.push(m);
            }
        }
        self.message_queue = kept;
        due.sort_by(|a, b| {
            a.time
                .0
                .partial_cmp(&b.time.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        due
    }

    /// Fail-loud check that a federate id is in range.
    fn check_federate(&self, id: FederateId) -> Result<()> {
        if id.0 >= self.federates.len() {
            return Err(err_unknown_federate(&format!("#{}", id.0)));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// GrantContext — the scoped bus access a behaviour gets on a grant
// ---------------------------------------------------------------------------

/// The scoped capabilities a [`FederateBehavior`] is handed when it is
/// granted a time: read its subscriptions, write its publications, send and
/// receive messages, and declare its next time request.
///
/// All name lookups are checked against ownership — a federate can only write
/// a publication it owns and send from an endpoint it owns — and are fail-loud
/// on an unknown or unauthorized name.
pub struct GrantContext<'b> {
    broker: &'b mut Broker,
    me: FederateId,
    granted: Time,
    inbox: Vec<Message>,
    next_request: &'b mut Option<Time>,
    finish: &'b mut bool,
}

impl GrantContext<'_> {
    /// The time this federate has just been granted.
    pub fn granted_time(&self) -> Time {
        self.granted
    }

    /// This federate's id.
    pub fn me(&self) -> FederateId {
        self.me
    }

    /// Publish `value` on publication `key`. Fail-loud if `key` is unknown or
    /// is owned by a different federate.
    pub fn publish(&mut self, key: &str, value: Value) -> Result<()> {
        let slot = self
            .broker
            .publications
            .get_mut(key)
            .ok_or_else(|| err_unknown_interface("publication", key))?;
        if slot.owner != self.me {
            return Err(err_not_owner("publication", key));
        }
        slot.value = Some((value, self.granted));
        Ok(())
    }

    /// Read the current value of subscription `sub_key` (the latest value on
    /// the publication it tracks), or `None` if nothing has been published
    /// there yet. Fail-loud if the subscription is unknown or owned by a
    /// different federate.
    pub fn value(&self, sub_key: &str) -> Result<Option<Value>> {
        let slot = self
            .broker
            .subscriptions
            .get(sub_key)
            .ok_or_else(|| err_unknown_interface("subscription", sub_key))?;
        if slot.owner != self.me {
            return Err(err_not_owner("subscription", sub_key));
        }
        // The publication may not exist / may be empty -> None.
        match self.broker.publications.get(&slot.publication_key) {
            Some(pubslot) => Ok(pubslot.value.as_ref().map(|(v, _)| v.clone())),
            None => Ok(None),
        }
    }

    /// Read a subscription as a `f64`, treating "no value yet" as the supplied
    /// `default`. Convenience for numeric coupling signals.
    pub fn value_or(&self, sub_key: &str, default: f64) -> Result<f64> {
        match self.value(sub_key)? {
            Some(v) => v.as_double(),
            None => Ok(default),
        }
    }

    /// Send `data` from endpoint `from` (which this federate must own) to
    /// endpoint `to`, deliverable at `deliver_at` (which must be `>=` the
    /// current granted time — a message cannot be sent into the past).
    ///
    /// Fail-loud on an unknown/unauthorized source endpoint, an unknown
    /// destination endpoint, or a past delivery time.
    pub fn send_message(
        &mut self,
        from: &str,
        to: &str,
        deliver_at: Time,
        data: Vec<u8>,
    ) -> Result<()> {
        let owner = self
            .broker
            .endpoints
            .get(from)
            .ok_or_else(|| err_unknown_interface("endpoint", from))?;
        if *owner != self.me {
            return Err(err_not_owner("endpoint", from));
        }
        if !self.broker.endpoints.contains_key(to) {
            return Err(err_unknown_interface("endpoint", to));
        }
        if !deliver_at.0.is_finite() || deliver_at < self.granted {
            return Err(err_past_message(from, to, deliver_at, self.granted));
        }
        self.broker.message_queue.push(Message {
            source: from.to_string(),
            destination: to.to_string(),
            time: deliver_at,
            data,
        });
        Ok(())
    }

    /// The messages delivered to this federate's endpoints at this grant, in
    /// delivery-time order. (Consumed once — call before returning.)
    pub fn messages(&self) -> &[Message] {
        &self.inbox
    }

    /// Request the next time this federate wants to be granted. The broker
    /// will snap it to the period grid and clamp it by the federation's
    /// safety bound. Must be `>=` the current granted time (enforced
    /// fail-loud after the behaviour returns).
    pub fn request_time(&mut self, t: Time) {
        *self.next_request = Some(t);
    }

    /// Declare that this federate is finished and disconnecting: it will not
    /// be granted again and stops constraining the rest of the federation.
    pub fn finish(&mut self) {
        *self.finish = true;
    }
}

// ---------------------------------------------------------------------------
// Subsystem-as-federate bridge
// ---------------------------------------------------------------------------

/// Wrap any [`crate::cosim::Subsystem`] (including a co-simulation FMU, once it
/// implements `Subsystem`) as a [`FederateBehavior`], so the existing native
/// coupling units and the FMI importer plug straight into a federation.
///
/// The bridge maps the subsystem's scalar ports onto named value-bus keys:
///
/// * `input_subs[i]` is the **subscription key** whose latest `Double` feeds
///   the subsystem's input port `i`.
/// * `output_pubs[j]` is the **publication key** the subsystem's output port
///   `j` is published to (as a `Double`) after each step.
///
/// On each grant the bridge reads every input subscription (missing values
/// default to `0.0`, matching the native master's unconnected-input rule),
/// steps the subsystem over `granted - last_granted`, publishes every output,
/// and requests its next periodic grant.
///
/// This is the **full wire** for the native path: a [`crate::cosim::Subsystem`]
/// (the same trait an imported co-sim FMU implements) becomes a first-class
/// time-coordinated federate. Stepping a *binary* `.so/.dll` FMU still sits
/// behind the off-by-default `binary-fmu` feature (see [`crate::fmi`]); that
/// binary, when present, also implements `Subsystem` and so rides this exact
/// bridge with no extra code.
pub struct SubsystemFederate<S: Subsystem> {
    subsystem: S,
    input_subs: Vec<String>,
    output_pubs: Vec<String>,
    /// The granted time of the *previous* grant, so the step size is the
    /// actual coordinated interval rather than a fixed assumption.
    last_granted: Time,
    /// The federate's period, used to schedule the next request.
    period: f64,
}

impl<S: Subsystem> SubsystemFederate<S> {
    /// Wrap `subsystem`, binding its inputs to `input_subs` (subscription
    /// keys) and its outputs to `output_pubs` (publication keys), to be driven
    /// on a `period` cadence.
    ///
    /// Fail-loud: the number of subscription keys must equal the subsystem's
    /// input arity, the number of publication keys must equal its output
    /// arity, and the period must be strictly positive (this driver is
    /// periodic).
    pub fn new(
        subsystem: S,
        input_subs: Vec<String>,
        output_pubs: Vec<String>,
        period: f64,
    ) -> Result<Self> {
        if input_subs.len() != subsystem.n_inputs() {
            return Err(err_arity(
                "input subscription keys",
                input_subs.len(),
                subsystem.n_inputs(),
            ));
        }
        if output_pubs.len() != subsystem.n_outputs() {
            return Err(err_arity(
                "output publication keys",
                output_pubs.len(),
                subsystem.n_outputs(),
            ));
        }
        if !period.is_finite() || period <= 0.0 {
            return Err(err_policy(format!(
                "SubsystemFederate period must be finite and > 0, got {period}"
            )));
        }
        Ok(SubsystemFederate {
            subsystem,
            input_subs,
            output_pubs,
            last_granted: Time::ZERO,
            period,
        })
    }
}

impl<S: Subsystem> FederateBehavior for SubsystemFederate<S> {
    fn on_grant(&mut self, granted: Time, ctx: &mut GrantContext<'_>) -> Result<()> {
        // Gather inputs from the subscriptions (default 0.0, like the native
        // master's unconnected inputs).
        let mut inputs = Vec::with_capacity(self.input_subs.len());
        for key in &self.input_subs {
            inputs.push(ctx.value_or(key, 0.0)?);
        }

        // Step over the actual coordinated interval. A zero-length first step
        // (granted == 0) just samples the initial outputs, matching the native
        // master's priming convention.
        let dt = granted.0 - self.last_granted.0;
        let outputs = self.subsystem.step(self.last_granted.0, dt, &inputs);
        if outputs.len() != self.output_pubs.len() {
            return Err(err_arity(
                "subsystem step outputs",
                outputs.len(),
                self.output_pubs.len(),
            ));
        }

        // Publish outputs by name.
        for (key, val) in self.output_pubs.iter().zip(outputs) {
            ctx.publish(key, Value::Double(val))?;
        }

        self.last_granted = granted;
        // Ask for the next periodic tick.
        ctx.request_time(Time(granted.0 + self.period));
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Null behaviour (placeholder while a real behaviour is borrowed out)
// ---------------------------------------------------------------------------

/// A do-nothing behaviour, used only as a temporary placeholder inside the
/// grant loop while the real behaviour is moved out to be run (so the broker
/// can be borrowed mutably). It is never actually granted.
struct NullBehavior;
impl FederateBehavior for NullBehavior {
    fn on_grant(&mut self, _granted: Time, _ctx: &mut GrantContext<'_>) -> Result<()> {
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Mock federate (closure-driven) — proves coordination without an FMU
// ---------------------------------------------------------------------------

/// A closure-driven federate, for tests and for proving the time-coordination
/// without a full FMU. The closure *is* the federate's behaviour: it advances
/// the federate's own (captured) state on each grant and may read/write the
/// buses through the [`GrantContext`].
///
/// Per the task: where a full FMU wire would be heavy, the coordination is
/// demonstrated with one of these — a closure that advances its own state on
/// grant. [`SubsystemFederate`] is the heavyweight, real-subsystem path.
pub mod mock {
    use super::{FederateBehavior, GrantContext, Result, Time};

    /// A federate whose behaviour is an arbitrary closure.
    pub struct MockFederate<F>
    where
        F: FnMut(Time, &mut GrantContext<'_>) -> Result<()>,
    {
        on_grant: F,
    }

    impl<F> MockFederate<F>
    where
        F: FnMut(Time, &mut GrantContext<'_>) -> Result<()>,
    {
        /// Build a mock federate from a closure run on each grant.
        pub fn new(on_grant: F) -> Self {
            MockFederate { on_grant }
        }

        /// Box it as a [`FederateBehavior`] for [`super::Broker::add_federate`].
        pub fn boxed(self) -> Box<dyn FederateBehavior>
        where
            F: 'static,
        {
            Box::new(self)
        }
    }

    impl<F> FederateBehavior for MockFederate<F>
    where
        F: FnMut(Time, &mut GrantContext<'_>) -> Result<()>,
    {
        fn on_grant(&mut self, granted: Time, ctx: &mut GrantContext<'_>) -> Result<()> {
            (self.on_grant)(granted, ctx)
        }
    }
}

// ---------------------------------------------------------------------------
// Error constructors (mapped onto the crate's existing FmiError taxonomy)
// ---------------------------------------------------------------------------
//
// The federation reuses the crate's single fail-loud error type, `FmiError`.
// A dedicated `FederationError` variant carries a stable machine-readable
// `code` plus a human message, so callers can match on the code while the
// existing FMI/DIS variants stay untouched.

/// Stable, machine-readable error codes for the federation layer. Matching on
/// these (rather than message text) is the supported way to handle a specific
/// failure. Carried inside [`FmiError::Federation`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FederationError {
    /// A time policy parameter was negative or non-finite, or a horizon was
    /// invalid.
    BadPolicy,
    /// A federate name was referenced that is not registered.
    UnknownFederate,
    /// A federate name collided with an already-registered federate.
    DuplicateFederate,
    /// A publication / subscription / endpoint name was referenced that is not
    /// registered.
    UnknownInterface,
    /// A publication / subscription / endpoint name collided with an existing
    /// one.
    DuplicateInterface,
    /// A federate tried to write an interface it does not own.
    NotOwner,
    /// A [`Value`] was read as the wrong type.
    TypeMismatch,
    /// A federate requested a time earlier than the one it was just granted.
    NonMonotoneRequest,
    /// A requested time was non-finite.
    BadRequest,
    /// A message was sent with a delivery time in the past.
    PastMessage,
    /// A port-count / key-count arity did not match.
    Arity,
    /// The grant limit was hit (a federate likely never advances).
    GrantLimitExceeded,
}

fn fed_err(code: FederationError, message: String) -> FmiError {
    FmiError::Federation { code, message }
}

fn err_policy(message: String) -> FmiError {
    fed_err(FederationError::BadPolicy, message)
}

fn err_unknown_federate(name: &str) -> FmiError {
    fed_err(
        FederationError::UnknownFederate,
        format!("unknown federate {name:?}"),
    )
}

fn err_dup_federate(name: &str) -> FmiError {
    fed_err(
        FederationError::DuplicateFederate,
        format!("federate {name:?} is already registered"),
    )
}

fn err_unknown_interface(kind: &str, name: &str) -> FmiError {
    fed_err(
        FederationError::UnknownInterface,
        format!("unknown {kind} {name:?}"),
    )
}

fn err_dup_interface(kind: &str, name: &str) -> FmiError {
    fed_err(
        FederationError::DuplicateInterface,
        format!("{kind} {name:?} is already registered"),
    )
}

fn err_not_owner(kind: &str, name: &str) -> FmiError {
    fed_err(
        FederationError::NotOwner,
        format!("federate does not own {kind} {name:?} and may not write it"),
    )
}

fn err_type(expected: &str, got: &Value) -> FmiError {
    fed_err(
        FederationError::TypeMismatch,
        format!("value is {}, expected {expected}", got.type_name()),
    )
}

fn err_nonmonotone_request(name: &str, req: Time, granted: Time) -> FmiError {
    fed_err(
        FederationError::NonMonotoneRequest,
        format!(
            "federate {name:?} requested time {} earlier than its granted time {}",
            req.0, granted.0
        ),
    )
}

fn err_bad_request(name: &str, req: Time) -> FmiError {
    fed_err(
        FederationError::BadRequest,
        format!("federate {name:?} requested non-finite time {}", req.0),
    )
}

fn err_past_message(from: &str, to: &str, at: Time, now: Time) -> FmiError {
    fed_err(
        FederationError::PastMessage,
        format!(
            "message {from:?} -> {to:?} delivery time {} is before the current grant {}",
            at.0, now.0
        ),
    )
}

fn err_arity(what: &str, got: usize, expected: usize) -> FmiError {
    fed_err(
        FederationError::Arity,
        format!("{what}: got {got}, expected {expected}"),
    )
}

fn err_grant_limit(max: usize) -> FmiError {
    fed_err(
        FederationError::GrantLimitExceeded,
        format!(
            "federation did not complete within {max} grants — a federate may \
             never advance (check periods / lookahead)"
        ),
    )
}

// ---------------------------------------------------------------------------
// Benchmark-pinned tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::mock::MockFederate;
    use super::*;
    use std::cell::RefCell;
    use std::rc::Rc;

    /// Shared log of `(granted_time, value_seen)` a consumer records.
    type SeenLog = Rc<RefCell<Vec<(f64, Option<f64>)>>>;
    /// Shared log of `(granted_time, message_bytes)` a receiver records.
    type MsgLog = Rc<RefCell<Vec<(f64, Vec<u8>)>>>;

    /// Helper: a federate that does nothing but march on its period (its
    /// request is set by nothing here — it relies on the broker's default
    /// "advance by period" when it makes no request).
    fn passive() -> Box<dyn FederateBehavior> {
        MockFederate::new(|_t, _ctx| Ok(())).boxed()
    }

    /// Helper: a federate that on each grant requests the next multiple of
    /// `period` strictly after its granted time. This makes the schedule
    /// explicit and independent of the broker's no-request default.
    fn periodic_requester(period: f64) -> Box<dyn FederateBehavior> {
        MockFederate::new(move |t: Time, ctx: &mut GrantContext<'_>| {
            ctx.request_time(Time(t.0 + period));
            Ok(())
        })
        .boxed()
    }

    // --- Benchmark #1: two periodic federates, periods 1.0 and 0.5, over
    //     [0, 3] -> the EXACT sequence of granted times. ------------------
    #[test]
    fn two_periodic_federates_exact_grant_sequence() {
        let mut b = Broker::new(3.0).unwrap();
        let fast = b
            .add_federate("fast", TimePolicy::periodic(0.5).unwrap(), passive())
            .unwrap();
        let slow = b
            .add_federate("slow", TimePolicy::periodic(1.0).unwrap(), passive())
            .unwrap();

        let trace = b.run_until(1000).unwrap();

        // Expected interleaving of grant times. Both start able to be granted
        // at 0 (initial grant), then march on their periods. Ties at a shared
        // time break by federate index (fast=0 before slow=1). The horizon is
        // *inclusive*: a federate is granted up to and including t = 3.0 (the
        // final time), then stops (its granted time has reached the horizon).
        //
        //  fast grid (period 0.5): 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0
        //  slow grid (period 1.0): 0,      1.0,      2.0,      3.0
        // Sanity: the indices in the trace are exactly the two federate ids.
        assert_eq!(fast, FederateId(0));
        assert_eq!(slow, FederateId(1));

        let times: Vec<(usize, f64)> = trace.iter().map(|r| (r.federate.0, r.time.0)).collect();
        assert_eq!(
            times,
            vec![
                (fast.0, 0.0), // fast @ 0
                (slow.0, 0.0), // slow @ 0
                (fast.0, 0.5), // fast @ 0.5
                (fast.0, 1.0), // fast @ 1.0
                (slow.0, 1.0), // slow @ 1.0
                (fast.0, 1.5), // fast @ 1.5
                (fast.0, 2.0), // fast @ 2.0
                (slow.0, 2.0), // slow @ 2.0
                (fast.0, 2.5), // fast @ 2.5
                (fast.0, 3.0), // fast @ 3.0 (inclusive horizon)
                (slow.0, 3.0), // slow @ 3.0 (inclusive horizon)
            ],
            "exact granted-time interleaving for periods 0.5 and 1.0 over [0,3]"
        );
    }

    // --- Benchmark #2: lookahead on one federate bounds another's grant. --
    #[test]
    fn lookahead_bounds_other_federates_grant() {
        // `producer` is event-driven with a lookahead of 0.4 and a passive
        // behaviour: it takes a single grant at time 0 and then *parks* — it
        // makes no further request, so its target stays at 0 (its last output
        // time). Its dependency bound is therefore `target(0) + lookahead(0.4)
        // = 0.4` forever: nothing it does can reach another federate before
        // 0.4, and (because it has no new output planned) it never advances.
        //
        // `consumer` is event-driven and, once granted, wants to leap straight
        // to 5.0 — but every advance past time 0 is clamped by the producer's
        // 0.4 bound.
        let mut b = Broker::new(10.0).unwrap();

        let producer_pol = TimePolicy::event_driven().with_lookahead(0.4).unwrap();
        let producer = b.add_federate("producer", producer_pol, passive()).unwrap();

        // Consumer wants to jump straight to 5.0 on every grant.
        let consumer = b
            .add_federate(
                "consumer",
                TimePolicy::event_driven(),
                MockFederate::new(|_t: Time, ctx: &mut GrantContext<'_>| {
                    ctx.request_time(Time(5.0));
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();

        // g0: both ungranted with target 0; tie -> producer (index 0) at 0.
        let g0 = b.step().unwrap().unwrap();
        assert_eq!(g0.federate, producer);
        assert_eq!(g0.time, Time(0.0));

        // g1: consumer takes its own initial grant at 0 (its target is 0 until
        // it has been granted once); there it requests 5.0.
        let g1 = b.step().unwrap().unwrap();
        assert_eq!(g1.federate, consumer);
        assert_eq!(g1.time, Time(0.0));

        // g2: consumer now wants 5.0, but the producer's bound
        // (target 0 + lookahead 0.4) clamps it to exactly 0.4.
        let g2 = b.step().unwrap().unwrap();
        assert_eq!(g2.federate, consumer);
        assert_eq!(
            g2.time,
            Time(0.4),
            "consumer's leap to 5.0 is clamped to producer's lookahead bound (0 + 0.4)"
        );
    }

    // --- Benchmark #3: a published value is visible to a subscriber at the
    //     correct granted time. ------------------------------------------
    #[test]
    fn published_value_visible_to_subscriber_at_correct_time() {
        // Producer publishes its granted time (as a Double) onto "src" every
        // 1.0. Consumer subscribes to "src" via "sub" and, on each grant,
        // records (granted_time, value_seen). With Gauss-Seidel-style index
        // ordering, at a shared time the producer (index 0) publishes before
        // the consumer (index 1) reads within the SAME coordinated time.
        let seen: SeenLog = Rc::new(RefCell::new(Vec::new()));

        let mut b = Broker::new(3.0).unwrap();
        let producer = b
            .add_federate(
                "producer",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(|t: Time, ctx: &mut GrantContext<'_>| {
                    // Publish the current time so the subscriber can check it.
                    ctx.publish("src", Value::Double(t.0))?;
                    ctx.request_time(Time(t.0 + 1.0));
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        b.register_publication(producer, "src").unwrap();

        let seen_c = Rc::clone(&seen);
        let consumer = b
            .add_federate(
                "consumer",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(move |t: Time, ctx: &mut GrantContext<'_>| {
                    let v = ctx.value("sub")?.map(|v| v.as_double().unwrap());
                    seen_c.borrow_mut().push((t.0, v));
                    ctx.request_time(Time(t.0 + 1.0));
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        b.register_subscription(consumer, "sub", "src").unwrap();

        b.run_until(1000).unwrap();

        // At each shared time t in {0,1,2}, the producer publishes t and then
        // the consumer (higher index) reads it in the same coordinated time.
        let got = seen.borrow().clone();
        assert_eq!(
            got,
            vec![(0.0, Some(0.0)), (1.0, Some(1.0)), (2.0, Some(2.0))],
            "subscriber sees the producer's value published at the same granted time"
        );

        // And the publication slot holds the last value (2.0) at time 2.0.
        let (val, when) = b.publication_value("src").unwrap().unwrap();
        assert_eq!(val, Value::Double(2.0));
        assert_eq!(when, Time(2.0));
    }

    // --- Benchmark #4: a known small federation -> EXACT grant-time trace.
    //     Three nested-rate periodic federates, the canonical HELICS multi-rate
    //     schedule, pinned over [0, 1]. (Lookahead-bounding is covered
    //     separately and exhaustively by benchmark #2 above.) ----------------
    #[test]
    fn small_mixed_federation_exact_trace() {
        // Nested rates, all zero-lookahead, passive (each marches on its own
        // period):
        //   f0: period 1.0  -> acts at 0, 1
        //   f1: period 0.5  -> acts at 0, 0.5, 1
        //   f2: period 0.25 -> acts at 0, 0.25, 0.5, 0.75, 1
        // Horizon 1.0. At any shared time the lower-indexed federate is granted
        // first (the tie-break), which is exactly the Gauss-Seidel publish-
        // before-read ordering. The schedule is fully determined by the
        // grant/gating rule (a periodic federate is granted only on its grid
        // point, and only once the safety frontier — the minimum of the other
        // federates' next action times — has reached that point).
        let mut b = Broker::new(1.0).unwrap();
        b.add_federate("f0", TimePolicy::periodic(1.0).unwrap(), passive())
            .unwrap();
        b.add_federate("f1", TimePolicy::periodic(0.5).unwrap(), passive())
            .unwrap();
        b.add_federate("f2", TimePolicy::periodic(0.25).unwrap(), passive())
            .unwrap();

        let trace = b.run_until(1000).unwrap();
        let times: Vec<(usize, f64)> = trace.iter().map(|r| (r.federate.0, r.time.0)).collect();

        assert_eq!(
            times,
            vec![
                (0, 0.0),  // f0 @ 0  (tie at 0 -> indices 0,1,2 in order)
                (1, 0.0),  // f1 @ 0
                (2, 0.0),  // f2 @ 0
                (2, 0.25), // f2 @ 0.25 (only f2's grid point here)
                (1, 0.5),  // f1 @ 0.5 (tie at 0.5 -> f1 before f2)
                (2, 0.5),  // f2 @ 0.5
                (2, 0.75), // f2 @ 0.75
                (0, 1.0),  // f0 @ 1   (tie at 1 -> f0, f1, f2 in order)
                (1, 1.0),  // f1 @ 1
                (2, 1.0),  // f2 @ 1
            ],
            "exact nested-rate (1.0 / 0.5 / 0.25) grant interleaving over [0,1]"
        );
    }

    // --- Fail-loud coverage --------------------------------------------
    #[test]
    fn negative_period_is_fail_loud() {
        let err = TimePolicy::periodic(-1.0).unwrap_err();
        assert!(matches!(
            err,
            FmiError::Federation {
                code: FederationError::BadPolicy,
                ..
            }
        ));
    }

    #[test]
    fn negative_lookahead_is_fail_loud() {
        let err = TimePolicy::event_driven().with_lookahead(-0.1).unwrap_err();
        assert!(matches!(
            err,
            FmiError::Federation {
                code: FederationError::BadPolicy,
                ..
            }
        ));
    }

    #[test]
    fn unknown_federate_name_is_fail_loud() {
        let b = Broker::new(1.0).unwrap();
        let err = b.federate_id("nope").unwrap_err();
        assert!(matches!(
            err,
            FmiError::Federation {
                code: FederationError::UnknownFederate,
                ..
            }
        ));
    }

    #[test]
    fn duplicate_federate_name_is_fail_loud() {
        let mut b = Broker::new(1.0).unwrap();
        b.add_federate("dup", TimePolicy::event_driven(), passive())
            .unwrap();
        let err = b
            .add_federate("dup", TimePolicy::event_driven(), passive())
            .unwrap_err();
        assert!(matches!(
            err,
            FmiError::Federation {
                code: FederationError::DuplicateFederate,
                ..
            }
        ));
    }

    #[test]
    fn unknown_publication_is_fail_loud() {
        let b = Broker::new(1.0).unwrap();
        let err = b.publication_value("ghost").unwrap_err();
        assert!(matches!(
            err,
            FmiError::Federation {
                code: FederationError::UnknownInterface,
                ..
            }
        ));
    }

    #[test]
    fn publishing_unowned_interface_is_fail_loud() {
        // f_a owns "p"; f_b tries to publish to it -> NotOwner.
        let mut b = Broker::new(1.0).unwrap();
        let a = b
            .add_federate("a", TimePolicy::event_driven(), passive())
            .unwrap();
        b.register_publication(a, "p").unwrap();
        let _bb = b
            .add_federate(
                "b",
                TimePolicy::event_driven(),
                MockFederate::new(|_t: Time, ctx: &mut GrantContext<'_>| {
                    // b does not own "p".
                    let r = ctx.publish("p", Value::Double(1.0));
                    assert!(matches!(
                        r,
                        Err(FmiError::Federation {
                            code: FederationError::NotOwner,
                            ..
                        })
                    ));
                    ctx.finish();
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        // Run a couple of grants so b's behaviour fires.
        b.run_until(10).unwrap();
    }

    #[test]
    fn nonmonotone_request_is_fail_loud() {
        // A federate granted at >0 that requests a time in its past must error.
        let mut b = Broker::new(5.0).unwrap();
        let _f = b
            .add_federate(
                "rewind",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(|t: Time, ctx: &mut GrantContext<'_>| {
                    if t.0 >= 1.0 {
                        // Ask to go backwards.
                        ctx.request_time(Time(0.5));
                    } else {
                        ctx.request_time(Time(t.0 + 1.0));
                    }
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        let err = b.run_until(100).unwrap_err();
        assert!(matches!(
            err,
            FmiError::Federation {
                code: FederationError::NonMonotoneRequest,
                ..
            }
        ));
    }

    #[test]
    fn value_type_mismatch_is_fail_loud() {
        let v = Value::Double(1.0);
        assert!(matches!(
            v.as_int(),
            Err(FmiError::Federation {
                code: FederationError::TypeMismatch,
                ..
            })
        ));
    }

    // --- SubsystemFederate wire: a native Subsystem as a federate. -------
    #[test]
    fn subsystem_federate_couples_two_integrators() {
        // Two trivial "accumulator" subsystems wired through the value bus:
        //   A: y_A += u_A * dt ; reads B's output, publishes its own.
        //   B: y_B += u_B * dt ; reads A's output, publishes its own.
        // This proves a real `Subsystem` rides the federation as a federate
        // and exchanges values by name across coordinated time.
        use crate::cosim::Subsystem;

        struct Accum {
            y: f64,
        }
        impl Subsystem for Accum {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, dt: f64, inputs: &[f64]) -> Vec<f64> {
                // Integrate the input over the step (forward Euler).
                self.y += inputs[0] * dt;
                vec![self.y]
            }
        }

        let mut b = Broker::new(2.0).unwrap();

        // A reads "b_out" (default 0), publishes "a_out".
        let a = b
            .add_federate(
                "A",
                TimePolicy::periodic(1.0).unwrap(),
                Box::new(
                    SubsystemFederate::new(
                        Accum { y: 1.0 },
                        vec!["a_in".to_string()],
                        vec!["a_out".to_string()],
                        1.0,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        b.register_publication(a, "a_out").unwrap();
        b.register_subscription(a, "a_in", "b_out").unwrap();

        // B reads "a_out", publishes "b_out".
        let bb = b
            .add_federate(
                "B",
                TimePolicy::periodic(1.0).unwrap(),
                Box::new(
                    SubsystemFederate::new(
                        Accum { y: 0.0 },
                        vec!["b_in".to_string()],
                        vec!["b_out".to_string()],
                        1.0,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        b.register_publication(bb, "b_out").unwrap();
        b.register_subscription(bb, "b_in", "a_out").unwrap();

        b.run_until(1000).unwrap();

        // A starts at y=1, B at y=0. Both publish their initial state at t=0
        // (dt=0 step). Over [0,2] with period 1 and Gauss-Seidel index order,
        // the exact values are deterministic; we pin A's final published value.
        let (a_val, _) = b.publication_value("a_out").unwrap().unwrap();
        let (b_val, _) = b.publication_value("b_out").unwrap().unwrap();
        // Hand-trace (Gauss-Seidel order A before B each shared time):
        //   t=0: A.step dt=0 -> y_A=1, pub a_out=1; B reads a_out=1, dt=0 ->
        //        y_B=0, pub b_out=0.
        //   t=1: A reads b_out=0, dt=1 -> y_A=1+0*1=1, pub 1; B reads a_out=1,
        //        dt=1 -> y_B=0+1*1=1, pub 1.
        //   t=2 is the horizon: a request to 2.0 retires the federate (>=
        //   horizon), so the last *grant* is at t=1. Final pubs: a_out=1,
        //   b_out=1.
        assert_eq!(a_val, Value::Double(1.0));
        assert_eq!(b_val, Value::Double(1.0));
    }

    #[test]
    fn subsystem_federate_arity_mismatch_is_fail_loud() {
        use crate::cosim::Subsystem;
        struct OneIn;
        impl Subsystem for OneIn {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, _inputs: &[f64]) -> Vec<f64> {
                vec![0.0]
            }
        }
        // Two subscription keys for a one-input subsystem -> Arity error.
        // (`SubsystemFederate` is not `Debug`, so match on the result rather
        // than calling `unwrap_err`.)
        let result = SubsystemFederate::new(
            OneIn,
            vec!["x".to_string(), "y".to_string()],
            vec!["z".to_string()],
            1.0,
        );
        assert!(matches!(
            result,
            Err(FmiError::Federation {
                code: FederationError::Arity,
                ..
            })
        ));
    }

    // --- Message endpoints: a timed message delivered at the right time. -
    #[test]
    fn endpoint_message_delivered_at_scheduled_time() {
        let received: MsgLog = Rc::new(RefCell::new(Vec::new()));

        let mut b = Broker::new(5.0).unwrap();

        // Sender: on its first grant (t=0) sends a message deliverable at 2.0.
        let sender = b
            .add_federate(
                "sender",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(|t: Time, ctx: &mut GrantContext<'_>| {
                    if t.0 == 0.0 {
                        ctx.send_message("tx", "rx", Time(2.0), vec![0xAB])?;
                    }
                    ctx.request_time(Time(t.0 + 1.0));
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        b.register_endpoint(sender, "tx").unwrap();

        let rec_c = Rc::clone(&received);
        let receiver = b
            .add_federate(
                "receiver",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(move |t: Time, ctx: &mut GrantContext<'_>| {
                    for m in ctx.messages() {
                        rec_c.borrow_mut().push((t.0, m.data.clone()));
                    }
                    ctx.request_time(Time(t.0 + 1.0));
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        b.register_endpoint(receiver, "rx").unwrap();

        b.run_until(1000).unwrap();

        // The message must be delivered exactly at the grant where the
        // receiver's time first reaches 2.0 (its period lands on 2.0).
        let got = received.borrow().clone();
        assert_eq!(got, vec![(2.0, vec![0xAB])], "message delivered at t=2.0");
    }

    #[test]
    fn message_to_unknown_endpoint_is_fail_loud() {
        let mut b = Broker::new(5.0).unwrap();
        let sender = b
            .add_federate(
                "s",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(|_t: Time, ctx: &mut GrantContext<'_>| {
                    let r = ctx.send_message("tx", "nowhere", Time(1.0), vec![1]);
                    assert!(matches!(
                        r,
                        Err(FmiError::Federation {
                            code: FederationError::UnknownInterface,
                            ..
                        })
                    ));
                    ctx.finish();
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        b.register_endpoint(sender, "tx").unwrap();
        b.run_until(10).unwrap();
    }

    #[test]
    fn past_message_delivery_is_fail_loud() {
        let mut b = Broker::new(5.0).unwrap();
        let s = b
            .add_federate(
                "s",
                TimePolicy::periodic(1.0).unwrap(),
                MockFederate::new(|t: Time, ctx: &mut GrantContext<'_>| {
                    if t.0 >= 2.0 {
                        // Try to deliver at 1.0, in the past.
                        let r = ctx.send_message("e", "e", Time(1.0), vec![1]);
                        assert!(matches!(
                            r,
                            Err(FmiError::Federation {
                                code: FederationError::PastMessage,
                                ..
                            })
                        ));
                        ctx.finish();
                    } else {
                        ctx.request_time(Time(t.0 + 1.0));
                    }
                    Ok(())
                })
                .boxed(),
            )
            .unwrap();
        b.register_endpoint(s, "e").unwrap();
        b.run_until(100).unwrap();
    }

    #[test]
    fn lone_federate_marches_on_its_period_to_horizon() {
        // A single periodic federate with no peers is bounded only by the
        // horizon: it should be granted 0, 1, 2 and then retire at 3.
        let mut b = Broker::new(3.0).unwrap();
        let _f = b
            .add_federate(
                "solo",
                TimePolicy::periodic(1.0).unwrap(),
                periodic_requester(1.0),
            )
            .unwrap();
        let trace = b.run_until(100).unwrap();
        let times: Vec<f64> = trace.iter().map(|r| r.time.0).collect();
        assert_eq!(times, vec![0.0, 1.0, 2.0]);
    }
}
