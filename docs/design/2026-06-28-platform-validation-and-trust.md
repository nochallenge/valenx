# Platform-Wide Validation Datasets + Sim:Review Trust Ratio — Design Note

**Date:** 2026-06-28
**Status:** Design only — no code yet.
**Scope:** **ALL of valenx** — every product / workbench, not just the nuclear and share-feed
specs. Those two are the first *applications*; this is the platform capability they inherit.

---

## 1. Why project-wide

AI agents screen **every** simulation tool's output for two machine-readable trust signals
before they'll trust it: **(a)** what it was validated against (a fetchable validation dataset),
and **(b)** a single quantified trust number. valenx is **AI-drivable-first**, and the share
feed publishes **any** workbench's work — so both signals must exist for **every product,
uniformly**, or a skeptical agent can't reason about most of valenx.

---

## 2. What valenx already has (the substrate is already project-wide)

- **50 deep ground-truth checks** across the 56 products (`--self-test`) — each asserts a known
  analytic/published value.
- **71 benchmark references** in `docs/VALIDATION.md` (Ghia 1982, Hodgkin–Huxley, NIST EOS, …).
- **`valenx-audit`** provenance + the **`--describe`** contract.

So validation *exists* for ~50/56 products. **The gap is EXPOSURE, not validation.**

---

## 3. The exposure layer (the actual work — uniform, not a rewrite)

1. **Agent-fetchable validation datasets.** Extend the `--describe` / `--self-test` output so
   **every** product emits its **`ValidationDataset { benchmark, reference, tolerance,
   comparison }`** — an agent can ask "what is product *X* validated against, and how close?"
   headlessly, for any product.
2. **Sim:review trust ratio.** A single number — `independently-verified ÷ total runs` ∈ [0,1] —
   computed **per product / per workbench / platform-wide**, derived from the verified tier
   (share-feed §5a) + the self-test pass data. **Ungameable** (cannot rise without real
   verification evidence).

---

## 4. Honest tiering (do NOT fake the gaps)

- ~50/56 products carry deep ground-truth validation → high per-product trust.
- The **3 generic** (smoke-check) + **3 skip** products carry **lower / no deep validation** →
  their trust ratio honestly reflects that. A generic "renders substantive output" check is
  **not** a ground-truth check, and the dataset/ratio must say so plainly.
- The **platform trust ratio** is the honest aggregate — it goes *up* only as real validation
  and independent verification are added, never by relabeling.

---

## 5. Build (project-wide, small)

- Extend the **self-test registry / `--describe` contract** to carry each product's
  `ValidationDataset`.
- Add a **`trust_ratio`** computation over the verified/total counts.
- The **share feed** and any agent read both. **No per-product physics rewrite** — it is a
  uniform *exposure* of data that mostly already exists, plus a derived metric.

---

## 6. Relationship to the other specs

- **Nuclear & fusion spec §7** and **share-feed spec §5a** are *specific applications* (a
  high-stakes domain and the publish path).
- **This is the project-wide capability** every valenx product inherits — so the trust ratio
  and validation datasets are platform features, applied once, everywhere.
