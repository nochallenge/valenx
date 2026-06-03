# RFC 0005: Design Principles

- **Status:** Accepted
- **Author(s):** BDFL
- **Created:** 2026-04-23
- **Discussion PR:** (this commit)
- **Tracking issue:** TBD

---

## Summary

Adopt the five design principles listed in
[DESIGN_PRINCIPLES.md](../DESIGN_PRINCIPLES.md) as binding
commitments for Valenx's user-facing design.

The principles are short enough to carry in your head:

1. Polish before scope
2. Dense, not busy
3. Two speeds, same tool
4. The app is the product
5. Quiet

This RFC exists to give those principles formal standing — any
future change to them goes through the RFC process.

---

## Motivation

A long-lived open-source project needs written values. Without
them, design drifts as contributors join and leave; contentious
calls ("should we ship this modal upsell?" "should the welcome
screen animate?") get relitigated endlessly; and the product loses
coherence year over year.

Five principles, written down, solve this cheaply. They're the
tie-breaker for trade-offs. They let reviewers say *"this violates
principle 4"* and have the violation be specific rather than
subjective.

Specific drivers:

- **We want new contributors to know the stance on telemetry,
  upsells, and cheer** the moment they open the repo. Principle 4
  makes that instant.
- **We want to protect against feature creep** — principle 1
  (polish before scope) is explicit about the order we build.
- **We want simulation users to feel the tool respects them** —
  principles 2 and 5 commit to calm, dense UX over consumer-app
  patterns.
- **We want the app to work for both novices and experts** without
  a "pro mode" toggle. Principle 3 is that commitment in one line.

---

## Guide-level explanation

The principles live in two places:

1. **[DESIGN_PRINCIPLES.md](../DESIGN_PRINCIPLES.md)** — short
   stand-alone file, ~1 page, the canonical source.
2. **[DESIGN.md Section 2](../DESIGN.md#2-design-principles--the-non-negotiables)**
   — the same five principles in the context of the full design
   plan, with background discussion.

Every design PR (anything touching UI) is expected to reference
which principle justifies a contentious choice, if it makes one.
Example review comment: *"This violates principle 4 — the 'free
trial' CTA on the Home screen is an upsell. Remove."*

The principles are **non-negotiable without RFC**. This is the
one place where "we just decided not to" isn't allowed.

---

## Reference-level explanation

The five principles, verbatim:

### 1. Polish before scope
One physics done beautifully beats five done roughly. Ship CFD to
Fusion quality before FEA half-done. This is a deliberate order,
not a rejection of breadth.

### 2. Dense, not busy
Simulation engineers live in numbers. Thirty data points on screen
is fine if they're organized; three data points with wasted space
is not. Fusion 360 is the aesthetic reference — dense, never
crowded.

### 3. Two speeds, same tool
A newcomer runs their first airfoil case by following prompts. An
expert does it in twelve keystrokes through the command palette.
Neither is bolted on. No "pro mode" toggle.

### 4. The app is the product
No modal upsells. No signup walls. No "complete your profile." No
telemetry by default. No marketing copy inside the UI. The app
earns attention by being good.

### 5. Quiet
Animations are subtle. The app does not make sounds. Colors are
restrained. The user's model and results are loud; the app is
quiet around them.

---

### Scope of application

These principles bind:

- **The desktop application** (`valenx-app` and all UI crates)
- **The in-app manual and tutorials** (voice, pacing, cheer level)
- **The public website** when it stands up (it's a sibling effort
  but same voice)
- **Plugin UI contributions** (plugins should follow them; the
  plugin reviewer enforces)

They do **not** bind:

- **Developer-facing output** — compiler errors, rustdoc, CLI
  debug logs — where verbosity serves the developer
- **Marketing materials** produced by third parties (e.g.,
  conference posters) — though the project's own marketing should
  follow principle 4

### How changes happen

Only through a successor RFC. A future maintainer who thinks one
of these is wrong writes a new RFC explaining why, it gets the
standard 10-day discussion window (per [rfcs/README.md](./README.md)),
and passes consensus-or-BDFL. That RFC amends or supersedes this
one.

No hallway-conversation override. No "the new designer says."
Five small rules, protected formally.

---

## Drawbacks

- **Principles this terse invite stretching.** "No exclamation
  points" gets argued as a violation of principle 5 when the real
  question is microcopy voice (DESIGN.md § 11). Principles are a
  starting point for judgment, not a substitute.
- **"Polish before scope" can be misused** to justify indefinitely
  delaying new physics verticals. It's an ordering rule for Year
  1, not a permanent excuse.
- **"Quiet" will get violated first** — it's the easiest to erode
  one "helpful" animation at a time. Requires active review
  enforcement.

---

## Rationale and alternatives

**More principles (e.g., ten like the Rust values).**
Rejected. Five fit in a head; ten don't. Every added principle
dilutes the signal.

**Fewer principles (two or three).**
Rejected. We tried three in early drafts; it left gaps (no
statement on the novice-expert duality, no commitment on
commercial posture). Five is the fewest that covers the space.

**No principles document; rely on code review culture.**
Rejected for a 20-year project. Culture drifts as maintainers
change; written principles don't.

**Make them aspirational ("we aim to...") rather than
non-negotiable.**
Rejected. Soft principles get ignored. These are the hill we
stand on.

---

## Prior art

- **Rust's values** — six short statements on reliability,
  performance, productivity, etc. Influential for Rust's culture
  over a decade. We borrow the format (short, memorable,
  non-negotiable).
- **Basecamp's "Shape Up" principles** — concrete, opinionated,
  short. Similar spirit.
- **Apple's Human Interface Guidelines** — longer but heavily
  curated; too much for our scale, but points toward the idea
  that principles compound over years.
- **Fusion 360's evident priorities** — never explicitly
  published, but we infer consistency, polish, and productivity
  focus; our principles 1, 2, 3 are partly derived from reading
  Fusion's behavior.

---

## Unresolved questions

- **Should principle 4 rule out opt-in telemetry entirely?**
  Current reading: opt-in is fine if user explicitly enables and
  can inspect payloads (DESIGN.md § 25). But principle 4's "no
  telemetry by default" leaves the opt-in question implicit. A
  future RFC may tighten this.
- **How to apply principle 2 (dense) to accessibility?** Dense UI
  can be cognitively harder; principles should not erode WCAG
  compliance. In practice they don't conflict — density comes
  from organized typography and spacing, not visual clutter. But
  worth watching.
- **Does principle 5 preclude completion chimes in any future
  world?** Current answer: yes, by default. Any future sound
  scheme requires a new RFC starting from zero.

---

## Future possibilities

- **Per-vertical sub-principles.** If CFD / FEA / EM end up with
  legitimately different UX needs (they probably do), each
  vertical could publish a short principles addendum that
  specializes these five without contradicting them.
- **Design-review rubric** auto-generated from the principles,
  used as a checklist during PR review.
- **External publication of the principles** — making the
  document public and referenceable by other simulation projects
  (MIT-licensed). "Valenx design principles" as a named artifact.
