//! A deterministic discrete-event scheduler — a min-heap event queue ordered by
//! simulated time.
//!
//! This is the time engine of the whole framework and is **fully general**: it
//! advances a single monotone simulated clock by popping the earliest pending
//! event, and an event handler may schedule further events into the future.
//! Nothing here knows anything about entities, sensors, or engagements — it is
//! the same discrete-event core a logistics, epidemiology, traffic, or
//! policy-wargaming model would use.
//!
//! ## Determinism (no wall clock)
//!
//! The clock is *simulated time*, seeded and driven entirely by the events you
//! enqueue — there is **no** `Instant::now`, `SystemTime`, or any wall-clock or
//! environment read anywhere in this crate. Two events with the *same* time pop
//! in **insertion order** (a monotonically increasing sequence number breaks the
//! tie), so a given set of events always replays in exactly the same order on
//! every run and machine. This is the property the benchmark test pins: a
//! shuffled batch of inserts pops out in non-decreasing time order.
//!
//! ## The monotonic-time invariant
//!
//! [`Scheduler::schedule`] rejects any event whose time is **strictly before**
//! the current clock (`MissionError::NonMonotonicEvent`) — you cannot schedule
//! into the past. Events *at* the current instant are allowed (zero-delay /
//! same-tick reactions) and run after everything already queued for that instant
//! in insertion order. The clock therefore never moves backward.

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use crate::error::MissionError;

/// A payload carried by a scheduled event.
///
/// The scheduler itself is payload-agnostic; a scenario chooses a concrete
/// payload type `E`. (`valenx-mission-sim`'s own scenario uses [`crate::Event`].)
pub trait EventPayload: Clone {}
impl<T: Clone> EventPayload for T {}

/// An event sitting in the queue: the time it fires, its insertion sequence
/// number (for a stable same-time tiebreak), and its payload.
#[derive(Debug, Clone)]
pub struct ScheduledEvent<E> {
    /// Simulated time at which the event fires (s).
    pub time: f64,
    /// Monotone insertion sequence number — breaks ties between equal times so
    /// ordering is fully deterministic.
    pub seq: u64,
    /// The event payload.
    pub payload: E,
}

impl<E> PartialEq for ScheduledEvent<E> {
    fn eq(&self, other: &Self) -> bool {
        self.time == other.time && self.seq == other.seq
    }
}
impl<E> Eq for ScheduledEvent<E> {}

impl<E> PartialOrd for ScheduledEvent<E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl<E> Ord for ScheduledEvent<E> {
    fn cmp(&self, other: &Self) -> Ordering {
        // Primary key: time. Event times are validated finite on the way in, so
        // `total_cmp` here only ever compares ordinary finite values; it is used
        // purely so the type can be `Ord`. Secondary key: insertion sequence,
        // giving a stable FIFO order among equal times.
        self.time
            .total_cmp(&other.time)
            .then_with(|| self.seq.cmp(&other.seq))
    }
}

/// A deterministic min-heap discrete-event scheduler over payloads of type `E`.
///
/// Wrap each entry in [`Reverse`] inside a [`BinaryHeap`] (a max-heap) so the
/// *earliest* time is what [`Scheduler::pop`] returns.
#[derive(Debug, Clone)]
pub struct Scheduler<E> {
    heap: BinaryHeap<Reverse<ScheduledEvent<E>>>,
    now: f64,
    next_seq: u64,
}

impl<E> Default for Scheduler<E> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E> Scheduler<E> {
    /// A fresh scheduler with the clock at `t = 0` and an empty queue.
    #[must_use]
    pub fn new() -> Self {
        Self {
            heap: BinaryHeap::new(),
            now: 0.0,
            next_seq: 0,
        }
    }

    /// The current simulated clock (s).
    #[must_use]
    pub fn now(&self) -> f64 {
        self.now
    }

    /// How many events are still pending.
    #[must_use]
    pub fn len(&self) -> usize {
        self.heap.len()
    }

    /// Whether the queue is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    /// Peek the time of the next event to fire, if any (does not advance).
    #[must_use]
    pub fn peek_time(&self) -> Option<f64> {
        self.heap.peek().map(|Reverse(ev)| ev.time)
    }

    /// Schedule `payload` to fire at absolute simulated time `time`.
    ///
    /// # Errors
    ///
    /// [`MissionError::NotFinite`] if `time` is not finite;
    /// [`MissionError::NonMonotonicEvent`] if `time` is strictly before the
    /// current clock (scheduling into the past is forbidden — events *at* `now`
    /// are allowed and run in insertion order).
    pub fn schedule(&mut self, time: f64, payload: E) -> Result<(), MissionError> {
        if !time.is_finite() {
            return Err(MissionError::NotFinite {
                quantity: "event time",
                value: time,
            });
        }
        if time < self.now {
            return Err(MissionError::NonMonotonicEvent {
                now: self.now,
                event_time: time,
            });
        }
        let seq = self.next_seq;
        self.next_seq += 1;
        self.heap
            .push(Reverse(ScheduledEvent { time, seq, payload }));
        Ok(())
    }

    /// Schedule `payload` to fire `delay` seconds *after* the current clock.
    ///
    /// # Errors
    ///
    /// [`MissionError::Negative`] if `delay` is negative or not finite.
    pub fn schedule_after(&mut self, delay: f64, payload: E) -> Result<(), MissionError> {
        let delay = crate::error::require_non_negative("event delay", delay)?;
        self.schedule(self.now + delay, payload)
    }

    /// Pop the earliest pending event, advancing the clock to its time.
    ///
    /// Returns `None` when the queue is empty (the clock is left where it was).
    /// Because the clock only ever moves *to* the popped event's time and
    /// [`Scheduler::schedule`] forbids past events, the returned times are
    /// non-decreasing across successive `pop` calls.
    pub fn pop(&mut self) -> Option<ScheduledEvent<E>> {
        let Reverse(ev) = self.heap.pop()?;
        // The clock advances to the event's time. It can only move forward
        // because past events were rejected at schedule time; `max` is a belt-
        // and-braces guard against any float edge so it never slips backward.
        self.now = self.now.max(ev.time);
        Some(ev)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- BENCHMARK PIN: a shuffled insert sequence pops out time-sorted -----
    #[test]
    fn pops_in_nondecreasing_time_order() {
        let mut s: Scheduler<&'static str> = Scheduler::new();
        // Deliberately scramble the insertion order.
        let times = [5.0, 1.0, 3.0, 1.0, 9.0, 2.0, 7.0, 0.0, 4.0, 1.0];
        for &t in &times {
            s.schedule(t, "e").unwrap();
        }
        let mut popped = Vec::new();
        let mut last = f64::NEG_INFINITY;
        while let Some(ev) = s.pop() {
            assert!(ev.time >= last, "{} popped after {}", ev.time, last);
            last = ev.time;
            popped.push(ev.time);
        }
        let mut sorted = times.to_vec();
        sorted.sort_by(f64::total_cmp);
        assert_eq!(popped, sorted);
    }

    #[test]
    fn equal_times_pop_in_insertion_order() {
        let mut s: Scheduler<u32> = Scheduler::new();
        // Five events all at t = 2.0, tagged 0..5 by insertion order.
        for tag in 0..5 {
            s.schedule(2.0, tag).unwrap();
        }
        let order: Vec<u32> = std::iter::from_fn(|| s.pop().map(|e| e.payload)).collect();
        assert_eq!(order, vec![0, 1, 2, 3, 4], "FIFO among equal times");
    }

    #[test]
    fn clock_advances_to_each_popped_event() {
        let mut s: Scheduler<()> = Scheduler::new();
        s.schedule(3.0, ()).unwrap();
        s.schedule(7.5, ()).unwrap();
        assert_eq!(s.now(), 0.0);
        s.pop();
        assert_eq!(s.now(), 3.0);
        s.pop();
        assert_eq!(s.now(), 7.5);
        assert!(s.pop().is_none());
        assert_eq!(s.now(), 7.5, "clock stays put on empty pop");
    }

    #[test]
    fn an_event_may_schedule_a_future_event() {
        let mut s: Scheduler<i32> = Scheduler::new();
        s.schedule(1.0, 1).unwrap();
        let mut fired = Vec::new();
        while let Some(ev) = s.pop() {
            fired.push(ev.payload);
            if ev.payload < 4 {
                // Each event spawns the next one 1 s later.
                s.schedule_after(1.0, ev.payload + 1).unwrap();
            }
        }
        assert_eq!(fired, vec![1, 2, 3, 4]);
        assert_eq!(s.now(), 4.0);
    }

    #[test]
    fn scheduling_into_the_past_is_rejected() {
        let mut s: Scheduler<()> = Scheduler::new();
        s.schedule(5.0, ()).unwrap();
        s.pop(); // clock -> 5.0
                 // Strictly before now -> rejected.
        assert!(matches!(
            s.schedule(4.999, ()),
            Err(MissionError::NonMonotonicEvent { .. })
        ));
        // Exactly at now -> allowed (zero-delay reaction).
        assert!(s.schedule(5.0, ()).is_ok());
        // Non-finite time -> rejected.
        assert!(matches!(
            s.schedule(f64::NAN, ()),
            Err(MissionError::NotFinite { .. })
        ));
    }

    #[test]
    fn empty_scheduler_pops_none() {
        let mut s: Scheduler<()> = Scheduler::new();
        assert!(s.is_empty());
        assert!(s.pop().is_none());
        assert_eq!(s.now(), 0.0);
        assert!(s.peek_time().is_none());
    }
}
