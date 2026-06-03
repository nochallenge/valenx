# Icon inventory

The authoritative list of every icon shipping in Valenx, keyed by
role. This file is the source of truth; the icon SVG under
`crates/valenx-icons/assets/` must exactly match this list.

**Status:** initial inventory drafted below; full Year-1 count
lands when the Home + Workspace mockups stabilize (Phase 0
deliverable).

---

## Format

```
| Role                  | Category | SVG path                         | Status |
|-----------------------|----------|----------------------------------|--------|
| run                   | action   | assets/action/run.svg            | TBD    |
| open-project          | action   | assets/action/open-project.svg   | TBD    |
```

Status values:
- **TBD** — listed, not yet drawn
- **Draft** — SVG exists, review pending
- **Ready** — merged, shipping

---

## Inventory

### Actions (the verbs users invoke)

| Role              | SVG                              | Status |
|-------------------|----------------------------------|--------|
| new-project       | `assets/action/new-project.svg`  | TBD    |
| open-project      | `assets/action/open-project.svg` | TBD    |
| save              | `assets/action/save.svg`         | TBD    |
| save-as           | `assets/action/save-as.svg`      | TBD    |
| run               | `assets/action/run.svg`          | TBD    |
| stop              | `assets/action/stop.svg`         | TBD    |
| pause             | `assets/action/pause.svg`        | TBD    |
| undo              | `assets/action/undo.svg`         | TBD    |
| redo              | `assets/action/redo.svg`         | TBD    |
| cut               | `assets/action/cut.svg`          | TBD    |
| copy              | `assets/action/copy.svg`         | TBD    |
| paste             | `assets/action/paste.svg`        | TBD    |
| delete            | `assets/action/delete.svg`       | TBD    |
| duplicate         | `assets/action/duplicate.svg`    | TBD    |
| export            | `assets/action/export.svg`       | TBD    |
| import            | `assets/action/import.svg`       | TBD    |
| search            | `assets/action/search.svg`       | TBD    |
| settings          | `assets/action/settings.svg`     | TBD    |
| help              | `assets/action/help.svg`         | TBD    |
| close             | `assets/action/close.svg`        | TBD    |

### Navigation and layout

| Role              | SVG                              | Status |
|-------------------|----------------------------------|--------|
| home              | `assets/nav/home.svg`            | TBD    |
| projects          | `assets/nav/projects.svg`        | TBD    |
| learn             | `assets/nav/learn.svg`           | TBD    |
| samples           | `assets/nav/samples.svg`         | TBD    |
| plugins           | `assets/nav/plugins.svg`         | TBD    |
| tools             | `assets/nav/tools.svg`           | TBD    |
| community         | `assets/nav/community.svg`       | TBD    |
| whats-new         | `assets/nav/whats-new.svg`       | TBD    |
| account           | `assets/nav/account.svg`         | TBD    |
| chevron-right     | `assets/nav/chevron-right.svg`   | TBD    |
| chevron-left      | `assets/nav/chevron-left.svg`    | TBD    |
| chevron-up        | `assets/nav/chevron-up.svg`      | TBD    |
| chevron-down      | `assets/nav/chevron-down.svg`    | TBD    |

### Physics (primary signal on project cards and ribbons)

| Role              | SVG                                  | Status |
|-------------------|--------------------------------------|--------|
| physics-cfd       | `assets/physics/cfd.svg`             | TBD    |
| physics-cfd-internal | `assets/physics/cfd-internal.svg` | TBD    |
| physics-fea       | `assets/physics/fea.svg`             | TBD    |
| physics-fea-modal | `assets/physics/fea-modal.svg`       | TBD    |
| physics-thermal   | `assets/physics/thermal.svg`         | TBD    |
| physics-em        | `assets/physics/em.svg`              | TBD    |
| physics-chem      | `assets/physics/chem.svg`            | TBD    |
| physics-md        | `assets/physics/md.svg`              | TBD    |
| physics-battery   | `assets/physics/battery.svg`         | TBD    |
| physics-multi     | `assets/physics/multi.svg`           | TBD    |

### Tools (registered solver / adapter indicators)

| Role                 | SVG                                   | Status |
|----------------------|---------------------------------------|--------|
| tool-openfoam        | `assets/tool/openfoam.svg`            | TBD    |
| tool-su2             | `assets/tool/su2.svg`                 | TBD    |
| tool-calculix        | `assets/tool/calculix.svg`            | TBD    |
| tool-code-aster      | `assets/tool/code-aster.svg`          | TBD    |
| tool-elmer           | `assets/tool/elmer.svg`               | TBD    |
| tool-openems         | `assets/tool/openems.svg`             | TBD    |
| tool-cantera         | `assets/tool/cantera.svg`             | TBD    |
| tool-lammps          | `assets/tool/lammps.svg`              | TBD    |
| tool-gmsh            | `assets/tool/gmsh.svg`                | TBD    |
| tool-freecad         | `assets/tool/freecad.svg`             | TBD    |

### Status glyphs

| Role                  | SVG                                     | Status |
|-----------------------|-----------------------------------------|--------|
| status-ready          | `assets/status/ready.svg`               | TBD    |
| status-running        | `assets/status/running.svg`             | TBD    |
| status-converged      | `assets/status/converged.svg`           | TBD    |
| status-diverged       | `assets/status/diverged.svg`            | TBD    |
| status-warning        | `assets/status/warning.svg`             | TBD    |
| status-error          | `assets/status/error.svg`               | TBD    |
| status-info           | `assets/status/info.svg`                | TBD    |
| status-not-run        | `assets/status/not-run.svg`             | TBD    |
| status-updating       | `assets/status/updating.svg`            | TBD    |
| status-unlocked       | `assets/status/unlocked.svg`            | TBD    |

---

## Workflow for adding an icon

1. Pick a role not already in this file; think about whether an
   existing role covers it (reuse > new icon).
2. Draw the 24×24 SVG. Follow the house style: stroke-based,
   consistent stroke weight, same corner radius as the rest of the
   set.
3. Add a row to the appropriate table, set status **Draft**.
4. Commit SVG under `crates/valenx-icons/assets/<category>/<name>.svg`.
5. Open PR; design-steward reviews.
6. On merge, status flips to **Ready**.

## Size estimates

Counting the rows above at the time of this draft:
**actions 20 + navigation 13 + physics 10 + tool 10 + status 10 = 63**
distinct icon roles. That's a working target, not a commitment —
real count will shift as the Year-1 mockups firm up, probably by
±10. Year-5 projection (all physics verticals in all states + every
solver + every ribbon action): roughly **180–220**, pending
verification from that era's mockups.
