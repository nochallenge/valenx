# Valenx — governance

How decisions get made, who can make them, and how the project moves
from "useful pre-alpha tool by one or two people" to "shippable
infrastructure with a real community."

> **Status: pre-alpha.** No formal foundation, no elected board, no
> commercial entity backing the project. Today the rules below describe
> what we're aiming at, not what has fully landed. As contributors land,
> the rules harden and we move from "BDFL with optional review" toward
> the steering committee model below. Every change to this document is a
> commit on `main` like any other — track it via `git log GOVERNANCE.md`.

## Roles

### Contributor
Anyone who opens a pull request, files an issue, or comments
constructively. No formal commitment. Bound by
[CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md).

### Reviewer
Trusted contributor with merge rights on a specific area
(an adapter crate, the workflow loop, a docs section). Reviewers
can:

- Approve and merge PRs touching only their area.
- Request changes that block a merge.
- Triage issues filed in their area.

Reviewers do **not** approve cross-area changes alone — those need
a maintainer.

### Maintainer
A reviewer who's also accepted responsibility for the project as a
whole. Maintainers can:

- Approve and merge cross-area PRs.
- Land RFCs (see [rfcs/README.md](./rfcs/README.md)).
- Cut releases (see "Releases" below).
- Add and remove reviewers and other maintainers.

Current maintainers live in [MAINTAINERS.md](./MAINTAINERS.md).

### Steering committee
Once there are five or more maintainers, they collectively form a
steering committee that:

- Decides on disputes that two maintainers can't resolve in PR review.
- Sets the high-level roadmap (changes to [ROADMAP.md](./ROADMAP.md)
  need committee sign-off).
- Approves new maintainers.
- Approves changes to this governance document and to
  [POLICIES.md](./POLICIES.md).

Until there are five maintainers, the existing maintainer set acts as
the committee. Decisions are by simple majority, with ties broken in
favour of the status quo.

## Decision-making

### Day-to-day code changes
Open a PR. A reviewer in the affected area approves; the contributor
or any maintainer merges. Small fixes (typos, doc spelling, dead
code) can be self-merged by maintainers without separate review.

### Bigger changes
File an RFC under [rfcs/](./rfcs/). The RFC process is documented in
[rfcs/README.md](./rfcs/README.md). Briefly:

1. Author opens a PR adding a new RFC document under `rfcs/`.
2. Discussion happens in PR review. Authors revise based on feedback.
3. After at least two weeks of discussion (or sooner if there's clear
   consensus), a maintainer either merges (accepting the RFC) or
   closes (rejecting it, with the document moved to `rfcs/rejected/`).

Examples of changes that need an RFC:

- Adapter contract changes (new lifecycle method, breaking trait
  signature change).
- New canonical types in `valenx-fields`, `valenx-mesh`, or
  `valenx-geo`.
- Project file format changes that would invalidate existing
  `.valenx` projects.
- Workflow DAG semantics changes.
- New top-level UI surfaces (a new pane, a new modal mode).

Examples that don't need an RFC:

- Adding a new live adapter that follows the existing pattern.
- New solver / analysis variants inside a live adapter (e.g.
  `rhoSimpleFoam` after `simpleFoam` was already shipping).
- UI polish (status badges, hover tooltips, command-palette
  entries).
- Test fixtures, internal helpers, doc updates.

When in doubt, file the RFC. The cost of an RFC that turns out
unnecessary is small; the cost of a non-RFC change that turns out
to break an established contract is large.

### Disputes
If two maintainers disagree on a PR or RFC outcome, any maintainer can
escalate to the steering committee by adding the `governance`
label and pinging the committee in the PR. The committee resolves
within a week, by simple majority. The decision is recorded as a
comment on the PR and (if it sets a precedent) referenced from the
relevant doc in `docs/decisions/`.

## Releases

### Cadence
Pre-alpha (today): release tags are best-effort, cut by any
maintainer when there's a meaningful body of work to mark
(`0.1.0-alpha.1` was the first; subsequent alphas as the workflow
loop evolves).

Once the project hits `0.1.0` (no alpha suffix):

- **Patch** (`0.1.x`): cut on a 2-4 week cadence as bugs are fixed.
  Any maintainer can cut.
- **Minor** (`0.x.0`): cut every 2-3 months, gathering features
  that landed since the last minor. Any maintainer can cut, but
  the changelog needs review by at least one other maintainer.
- **Major** (`x.0.0`): rare; only when the project file format,
  adapter contract, or another canonical interface needs a
  breaking change. Requires steering committee approval and an
  RFC documenting the migration path.

### Release process
1. Land all in-flight PRs to `main`.
2. Update `CHANGELOG.md`: move "Unreleased" to a new versioned
   section with today's date.
3. Tag the commit: `git tag -a v0.x.y -m "Release 0.x.y"`.
4. Push the tag: `git push origin v0.x.y`.
5. CI builds installers for all three platforms (Phase 10+).
6. Open a PR titled `release: 0.x.y` updating any version-bump
   files. Self-merge after CI green.

## Adding a new maintainer

1. An existing maintainer (or the candidate themselves) opens an
   issue titled `Maintainer nomination: <name>`. The issue lists:
   - Areas the candidate has been active in.
   - Concrete contributions (linked PRs / issues / RFCs).
   - Why the project benefits.
2. Existing maintainers vote in the issue thread. Approval needs a
   simple majority of current maintainers, with at least three
   approvals (or unanimous, if there are fewer than three
   maintainers).
3. On approval, the candidate updates [MAINTAINERS.md](./MAINTAINERS.md)
   in a PR, gets the appropriate GitHub permissions, and is welcomed
   in a release-notes mention.

Removing a maintainer follows the same process, opened by any other
maintainer or by the maintainer themselves (voluntary step-down).

## Conflicts of interest

If a maintainer's employer's product overlaps with a Valenx
adapter's tool (e.g. a maintainer working at a company that sells
proprietary CFD), they:

- Disclose the relationship in [MAINTAINERS.md](./MAINTAINERS.md).
- Recuse themselves from RFC votes and adapter design decisions
  for that specific tool.
- Don't recuse from general project decisions (governance, infra,
  unrelated adapters) — the disclosure is enough.

## Code of conduct enforcement

[CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md) enforcement is handled
by maintainers acting collectively. Reports go to the security
contact in [SECURITY.md](./SECURITY.md) (a private email channel).
Outcomes range from a private warning to a permanent ban; the
escalation ladder is in the Code of Conduct itself.

## Funding

Currently zero. Maintainers contribute on personal time. If the
project ever takes funding (donations, grants, commercial sponsorship)
those rules will need their own RFC and a transparency clause added
to this document.

## Amendments

Changes to this document are PRs like any other. Amendments need
steering committee approval (or, until the committee exists, a
majority of current maintainers). Trivial fixes (typos, broken
links) can be self-merged by any maintainer.
