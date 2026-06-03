# Valenx RFCs

RFC = **Request For Comments**. This is how big decisions on Valenx get
made.

The RFC process is modeled on [Rust's RFC process](https://github.com/rust-lang/rfcs)
— which is, in turn, the gold standard for open-source technical
governance. Shamelessly inspired.

---

## When you need an RFC

You need an RFC for anything substantial:

- New top-level crate in the workspace
- New integrated OSS tool (new adapter)
- Public API changes (Rust crate APIs, plugin ABI, scripting API, CLI)
- `.valenx` file-format changes
- Performance-critical algorithm changes
- Breaking changes of any kind
- Governance, process, or policy changes
- Adding or removing a major dependency

You do **not** need an RFC for:

- Bug fixes
- Documentation improvements
- New tests for existing features
- Refactoring that doesn't change behavior
- Implementing something that already has an accepted RFC
- Minor tweaks to existing features (judgment call — if you're unsure, file one)

---

## The process

```
Idea → Draft RFC → PR → Discussion → Accept / Reject / Postpone → Implement
```

### 1. Before writing, sanity-check

- Search existing RFCs (both merged and open) for prior art
- Open a GitHub Discussion to gauge interest
- Ask in chat — sometimes the RFC isn't needed, sometimes someone is
  already working on it

### 2. Write the draft

- Copy `0000-template.md` to `0000-my-feature.md` (keep `0000-`, PR
  assigns the number)
- Fill in every section honestly — "drawbacks" and "unresolved
  questions" are mandatory, not optional
- Be specific. "Make the CAD kernel faster" is not an RFC; "Cache BRep
  tessellation between viewport updates using a sparse LRU keyed on
  topology hash" is.

### 3. File the PR

- Branch name: `rfc/my-feature`
- PR title: `RFC: My Feature`
- PR description: link to the RFC file, TL;DR in 2-3 sentences
- Tag relevant maintainers

### 4. Discussion

- Minimum **10 calendar days** open before a merge decision
- Substantial changes to the draft reset the 10-day clock
- Comments on the PR are the canonical record — don't have the
  discussion in private channels
- The author drives revisions; maintainers help shape

### 5. Decision

A maintainer proposes a disposition:

- **Accept** — merge the RFC; `0000-` gets replaced with the next number
- **Reject** — close the PR with a comment explaining why (kept in
  `rfcs/rejected/` for reference)
- **Postpone** — close the PR, label "postponed"; reopen when conditions
  change (e.g., blocked on an upstream feature)

Consensus-seeking first. If no consensus after another 10 days, the BDFL
(later TSC) casts a tie-break. Dissenting opinions get captured in the
RFC itself before merge.

### 6. Implementation

- Accepted RFC gets a tracking issue
- PRs reference the RFC number in commit messages
- RFC file may be updated with amendments during implementation — with
  a follow-up PR, not silent edits

---

## RFC lifecycle states

| State | Meaning | Where it lives |
|-------|---------|----------------|
| **Draft** | Being written, not yet filed | Local branch or WIP PR |
| **Open** | PR filed, under discussion | GitHub PR |
| **Accepted** | Merged to `main` | `rfcs/XXXX-name.md` |
| **Implemented** | Code exists, shipping | Unchanged, tracking issue closed |
| **Amended** | Scope or design changed after acceptance | Same file, "Amendments" section at bottom |
| **Superseded** | A later RFC replaced it | File kept, header notes the successor |
| **Withdrawn / rejected** | Not pursued | `rfcs/rejected/XXXX-name.md` |

---

## Numbering

- RFCs are numbered sequentially starting at `0001`
- `0000-` is the placeholder for a draft not yet merged
- Numbers are assigned at merge time, not at PR time, to avoid merge
  conflicts when multiple RFCs land in the same week

---

## What a good RFC looks like

For examples worth copying the style of, browse the [Rust RFC
repository](https://github.com/rust-lang/rfcs/tree/master/text). The
ones that tend to hold up years later share a few traits:

- A clear problem statement that non-experts can follow
- Motivated cost/benefit analysis with honest numbers
- Worked examples of the API or syntax *before and after*
- A real "Drawbacks" section that lists trade-offs the author
  accepted, not just straw men
- "Alternatives considered" with reasons for rejection
- Explicit "Unresolved questions" — the author saying what they
  don't know

Patterns to copy:

- **Motivation before design.** Explain the problem in user terms first
- **Concrete examples** of the API or file before/after
- **Alternatives considered** — even ones you rejected
- **Drawbacks are honest** — every proposal has them
- **Unresolved questions** — say what you don't know

Patterns to avoid:

- "It's obviously better" — prove it with an example
- "We can figure that out later" — if a major question is unresolved,
  either answer it or punt the whole RFC
- Implementation details without motivation
- Quoting statistics without sources

---

## After acceptance: amendments

Accepted RFCs aren't set in stone, but they're not freely editable
either. Once merged:

- **Typos, clarifications, formatting** — direct PR, no process
- **Scope changes, design changes** — file an amendment PR that adds an
  "Amendments" section at the bottom of the RFC, with the old text
  preserved
- **Fundamental redesign** — supersede it with a new RFC; the old one
  stays on record

This matters because implementors reading the RFC in 2028 need to know
what was actually agreed to, not just what the latest rewrite says.

---

## Governance of this process itself

Changes to this RFC process require an RFC (meta, yes). When that
bootstraps, it will land as its own numbered RFC; until then, this
README is the authoritative source.

---

## Index

The numbered RFCs below are the accepted ones. Open PRs are the active
proposals.

- [0001 — `.valenx` project file format](./0001-project-file-format.md)
- [0002 — Adapter contract](./0002-adapter-contract.md)
- [0003 — Plugin API (WIT)](./0003-plugin-api.md)
- [0004 — Results and fields data model](./0004-results-and-fields.md)
- [0005 — Design principles](./0005-design-principles.md)
- [0006 — Design token system and pipeline](./0006-token-system.md)

More will be added as the project grows. See
[rfcs/rejected/](./rejected/) for RFCs that were filed and not
accepted (archived for reference).
