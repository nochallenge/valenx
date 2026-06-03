# RFC 0000: <Short descriptive title>

- **Status:** Draft
- **Author(s):** <GitHub handle(s)>
- **Created:** YYYY-MM-DD
- **Discussion PR:** <link to this PR>
- **Tracking issue:** <filled in after acceptance>

---

## Summary

One paragraph. What are you proposing? What changes?

Should be understandable standalone — a reader shouldn't need to read
the motivation to know what you're suggesting.

---

## Motivation

Why are we doing this? What problem does it solve?

- What use cases does it enable?
- What is the current state of the world, and what's wrong with it?
- If we don't do this, what happens?

Be specific. "It would be cleaner" is not motivation. "Users currently
have to write 40 lines of boilerplate to import a STEP file; this RFC
reduces that to 3 lines" is.

---

## Guide-level explanation

Explain the proposal as if teaching it to a developer or user for the
first time. If it's a user-facing feature, write it the way it would
appear in the user docs. If it's an API, show example code using it.

This section must:

- Define any new terminology
- Show concrete examples of the new functionality in use
- Explain how it affects existing users (if at all)
- Cover the common case; the weird edge cases go in the reference-level
  explanation

---

## Reference-level explanation

The details. Precise, complete, and specific enough that an engineer
could implement from this section alone.

Cover, as applicable:

- Exact API signatures (Rust types, function signatures, trait bounds)
- Exact file-format schema (TOML/JSON structure, field names, defaults)
- Algorithms and their complexity bounds
- Data structures
- Interaction with existing features
- Error handling
- Performance characteristics
- Security considerations
- How it integrates with the rest of the codebase

Include diagrams if they clarify. Mermaid or ASCII art both fine.

---

## Drawbacks

Why should we *not* do this?

Every proposal has costs. Be honest. Examples:

- Adds a new dependency
- Increases binary size by X MB
- Requires a migration for existing users
- Locks us into a design decision
- Learning curve for new contributors
- Performance cost

If this section is short, you haven't thought hard enough about
downsides.

---

## Rationale and alternatives

- Why is this design the best among alternatives?
- What other designs did you consider? Why did you reject them?
- What is the impact of not doing this?
- Are there existing solutions in other projects? How does theirs work,
  and why not just adopt that?

"We couldn't think of any alternatives" means either the problem is
poorly scoped or the author hasn't done enough research. Do the
research.

---

## Prior art

Discuss prior art — both the good and the bad — in relation to this
proposal:

- Other simulation suites (ANSYS, COMSOL, Siemens NX, SimScale)
- Other open-source projects (FreeCAD, Salome, SU2, OpenFOAM itself)
- Academic papers
- Industry standards (STEP AP242, CGNS, VTK formats)

Not every RFC has prior art — fundamentally novel designs are fine.
But if there's prior art you're deliberately diverging from, say why.

---

## Unresolved questions

What parts of the design are still open?

- Questions you want resolved before this RFC merges
- Questions you explicitly expect to defer to a future RFC
- Questions that are out of scope

Being up-front here protects both the author and the reviewers. Nothing
worse than an RFC merging and then six people discovering they each
thought a different thing was agreed to.

---

## Future possibilities

Optional. Where does this lead?

- Natural extensions that would build on this proposal
- Related features that would become easier / possible
- Ecosystem implications

Don't use this section to argue the current proposal. Use it to sketch
the road ahead.

---

## Appendices

Optional. Any of:

- Full API reference for proposed new types
- Worked end-to-end example
- Benchmark results
- Transcript of prototype runs
- Links to related issues / prior discussions

---

## Amendments

(Added post-acceptance only. Each amendment gets a dated heading with a
summary of what changed and why. Original text above stays intact.)
