# Continuous Integration — current state

> **Status (updated 2026-06-12):** `ci.yml` (the push/PR build + lint +
> cross-platform test matrix) is **re-enabled** — it fires on every push
> to `master`/`main` and on every PR, and the `CI OK` summary job is the
> single status check to require in branch protection. The repo is now
> **public**, so GitHub Actions standard runners are free and the
> original billing concern (below) no longer applies to `ci.yml`.
> `ci-nightly.yml` and `release.yml` remain **manual-trigger only**
> (`workflow_dispatch`) — re-enable those per the sections below if you
> want them (the nightly E2E suite and the release/installer pipeline).

---

## Why this is locked down

GitHub Actions billing on private repos uses runner-tier multipliers:

| Runner | Multiplier | Effective rate¹ |
| --- | --- | --- |
| `ubuntu-latest` | 1× | $0.008 / min |
| `windows-latest` | 2× | $0.016 / min |
| `macos-latest` | 10× | $0.080 / min |

¹ For private repos, after the free 2000 min/month allowance is used.

`release.yml` runs all three runners in parallel. A typical full run
for this repo (249+ crates, cold cache) consumes:

| Job | Runner | Duration | Billable minutes |
| --- | --- | --- | --- |
| linux .deb | `ubuntu-latest` | ~17 min | ~17 |
| linux .rpm | `ubuntu-latest` | ~18 min | ~18 |
| macOS .dmg | `macos-latest` | ~13 min | ~130 |
| Windows .msi | `windows-latest` | ~25 min | ~50 |
| **Total** | | | **~215 min** |

At the rates above, one full release run is roughly **$10-15 USD**
per tag push. A single hot-cache rebuild is ~30-40% of that.

`ci.yml` (push / PR matrix) hits ubuntu-only but ran every push to
master + every PR — typically 1-2 runs/day × ~10 min × $0.008 = a
few cents per day, but adds up to a few dollars per month over time.

`ci-nightly.yml` runs the E2E adapter integration tests at 4 AM UTC
nightly on ubuntu — another ~30 min/day × $0.008 = ~$7/month.

---

## What's still possible

Every workflow keeps its `workflow_dispatch` trigger, so you can
still kick off a one-off build from the GitHub UI or CLI when you
want — you just pay for that specific run instead of getting
surprised by auto-fires.

```sh
# CI (build + lint check on demand):
gh workflow run ci.yml

# Nightly E2E (install conda-forge tools + run integration suite):
gh workflow run ci-nightly.yml

# Release pipeline for an existing tag:
gh workflow run release.yml -f tag=v0.1.0-alpha.1
```

The release workflow also accepts a `draft` input (defaults to
`true`) so you can produce installer artifacts attached to a draft
GitHub Release without it going public until you explicitly publish.

---

## Building installers without GitHub Actions

For most one-off builds, prefer the local-build scripts — same
commands the CI jobs run, zero billable minutes:

| Format | Where to build | Command |
| --- | --- | --- |
| Windows `.msi` | Any Windows box with WiX Toolset 3.x | `pwsh scripts/build-installer.ps1` |
| macOS `.dmg` | Any Mac with `cargo-bundle` + `brew install create-dmg jq` | `bash scripts/build-installer.sh` |
| Linux `.deb` + `.rpm` | Linux / WSL with `cargo-deb` + `cargo-generate-rpm` + system libs | `bash scripts/build-installer.sh` |

Full prereqs + commands in `docs/INSTALLER.md`.

If you don't have a Mac, the macOS `.dmg` can be cross-built (with
caveats around signing) using `osxcross`, but the cleanest path is
borrowing a Mac for the one-off when you actually need to ship.

---

## Re-enabling automatic triggers

If you decide a particular workflow's cost is worth it, edit the
`on:` block in the workflow file. The original blocks are shown
below; paste them in place of the existing `workflow_dispatch:`-only
block.

### `.github/workflows/ci.yml` — original triggers

```yaml
on:
  push:
    branches: [master, main]
  pull_request:
    branches: [master, main]
  workflow_dispatch:
```

### `.github/workflows/ci-nightly.yml` — original triggers

```yaml
on:
  schedule:
    # 4 AM UTC every day — well after typical PR cycles
    - cron: '0 4 * * *'
  workflow_dispatch:
```

### `.github/workflows/release.yml` — original triggers

```yaml
on:
  push:
    tags: ['v*']
  workflow_dispatch:
    inputs:
      tag:
        description: 'Tag to build (e.g. v0.1.0-alpha.1)'
        required: true
      draft:
        description: 'Create as draft release?'
        required: false
        default: 'true'
        type: boolean
```

---

## If you want zero-risk no-GHA-ever

The most foolproof "no surprise billing" config is to delete the
workflow files entirely:

```sh
rm .github/workflows/{ci,ci-nightly,release}.yml
git add .github/workflows/ && git commit -m "remove all CI workflows"
git push origin master
```

After this, even the manual `gh workflow run` path is gone. To bring
CI back you'd need to recover the files from git history
(`git show 11d66d4:.github/workflows/release.yml > .github/workflows/release.yml`).

The current manual-only config keeps the workflow files alive so
you can re-enable triggers in seconds; the deletion option keeps you
honest about not running CI at all.

---

## Lower-cost runners (if you want CI but cheaper)

* **Self-hosted runners** — your own hardware, free minutes but you
  manage the runner. https://docs.github.com/en/actions/hosting-your-own-runners
* **Larger Linux runners only** — `ubuntu-latest-4-cores` etc. cost
  1× per minute but finish faster (so cheaper net for cold-cache
  Rust builds). Spec via `runs-on: ubuntu-latest-4-cores`.
* **Public repo** — GitHub Actions on public repos is **free** on
  all runners. Making the repo public would eliminate all billing.
* **Drop the macOS runner** — replace the `macos` job with a no-op
  for now; macOS users build the `.dmg` locally on their own Mac.
  Cuts ~60% of release-pipeline cost.

Pick whichever fits your goals.
