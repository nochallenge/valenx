# Valenx Design Principles

Five non-negotiables. When a design trade-off comes up, these break the tie.

The full reasoning, scope, and implementation plan live in
[DESIGN.md](./DESIGN.md). This file is the short version, small enough
to carry in your head.

---

## 1. Polish before scope

One physics done beautifully beats five done roughly. Ship CFD to Fusion
quality before FEA half-done. This is a deliberate order, not a
rejection of breadth.

## 2. Dense, not busy

Simulation engineers live in numbers. Thirty data points on screen is
fine if they're organized; three data points with wasted space is not.
**Fusion 360** is the aesthetic reference — dense, never crowded.

## 3. Two speeds, same tool

A newcomer runs their first airfoil case by following prompts. An
expert does it in twelve keystrokes through the command palette.
Neither is bolted on. No "pro mode" toggle.

## 4. The app is the product

No modal upsells. No signup walls. No "complete your profile." No
telemetry by default. No marketing copy inside the UI. The app earns
attention by being good.

## 5. Quiet

Animations are subtle. The app does not make sounds. Colors are
restrained. The user's model and results are loud; the app is quiet
around them.

---

Changes to this document require an RFC
(see [rfcs/README.md](./rfcs/README.md)). The principles are meant to be
argued with — but only through that process.
