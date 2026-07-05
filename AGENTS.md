# AGENTS.md

Guidance for AI agents (and other contributors) working in this repository.

**Read [ARCHITECTURE.md](ARCHITECTURE.md) before writing code.** It describes
the package structure, core abstractions, design principles, and the
step-by-step pattern for adding a new algorithm. This file is about how to
work: the approach, the checks, and the commands.

## How to approach design and implementation

The recurring temptation these guard against is **adding machinery where none
was needed**. They are defaults with reasons attached, not commandments — when
a real exception comes up, take it, and say why.

### When two pieces don't fit, don't reach for a bridge

Work through these in order and stop at the first that works:

1. **Fix or extend the shared abstraction at its source.** If a generic
   interface is slow or missing a method for your case, improve the generic
   implementation (or grow a small sibling method, e.g. `getprep` →
   `hasprep`) so every consumer benefits.
2. **Adapt the consumer.** Map, offset, sort, or iterate at the call site.
   A consumer-side adaptation is cheaper and more honest than a producer-side
   duplicate structure.
3. **Accept the mismatch.** Slow paths are allowed; `prepare` and friends are
   the escape hatch (see ARCHITECTURE.md design principles).

A new connecting structure (adapter, composed tree, parallel index, bespoke
traversal helper) is the last resort — not forbidden, but expensive: it forks
maintenance and helps only its own call site, where a fix to the shared
abstraction helps every consumer. If you conclude one is genuinely needed,
build it, but say so prominently and state the tradeoff rather than burying it
in a larger change. Beware the near-miss: avoiding a new *type* but building
new *composition machinery* is the same thing.

### Test genericity against the hardest case

Before proposing a design that claims to be generic, check it against the most
hostile instantiation: an opaque C-library tree, a foreign geometry with
expensive accessors, an unclosed ring. If the design only works for backends
already in the codebase, it isn't generic.

### Don't transform user data silently

No silent copying, re-materializing, ring-closing, or number-type conversion.
Anything expensive or semantics-changing is an explicit opt-in keyword. The
default is the most literal interpretation of what was passed in.

### Build only what has a consumer

No unwired "forward-looking" components, no trait hierarchies nothing
exercises, no type parameters for out-of-scope cases. If generality seems
architecturally tempting, note it as an option in your summary ("this could
grow to support X") instead of building it now.

### After generalizing, subtract

When a new mechanism subsumes an old special case, delete the old path in the
same change rather than leaving cleanup for a later pass. Lean against
single-use helpers: keep one when it earns its place (a performance barrier,
real reuse, a genuine readability win), inline it otherwise.

### Audit yourself for self-consistency

Before presenting a design, check it against principles already established —
in ARCHITECTURE.md, earlier in the session, or in the very file being edited
(e.g. an existing abstract supertype convention). It is much cheaper to catch
"this contradicts a convention I already follow elsewhere" yourself than in
review.

### Report honestly

- Don't launder scope decisions as architecture ("X belongs elsewhere" when the
  truth is "I didn't get to X").
- Every deferred item carries a one-line justification of why it's real work.
- No coined jargon in summaries; every claim should survive one round of
  "why?" without a follow-up.

### Comments describe what, not the story of why

Inline comments briefly state the constraint or behavior the code can't show.
Longer exposition has better homes — the literate header at the top of the
file (see ARCHITECTURE.md), docstrings, `docs/plans/`, or the commit message —
so design history and rationale don't end up narrated between lines of code.

## What consistently works

- Read the relevant code before proposing a design; verify claims against
  source rather than memory.
- After every structural change: run the affected test suites and benchmark
  against recorded baselines before committing.
- Surface genuinely open API choices as explicit questions instead of guessing.
- When challenged on a design, respond with a transparent tradeoff breakdown
  (what's irreducible vs. removable), not defensiveness or blind rework.
- Don't push, open PRs, or merge unless asked.

## Development commands

**IMPORTANT**: when running Julia, always run with `julia --project=docs` so
that you have access to utility packages for loading geometry, etc. If you get
an error about a package not being found, try running
`julia --project=docs -e 'using Pkg; Pkg.instantiate()'` first.

### Testing

Run all tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

To run a single test file in isolation, include it while in the test
environment (test files are self-contained — each one does its own
`using`/`import`):

```bash
julia --project=test -e 'include("test/methods/area.jl")'
```

### Testing against multiple geometry implementations

The **GeometryOpsTestHelpers/** subpackage provides `@test_implementations`
and `@testset_implementations`. These macros run a test block once per
geometry implementation (GeoInterface wrappers always; ArchGDAL,
GeometryBasics, and LibGEOS are added to
`GeometryOpsTestHelpers.TEST_MODULES` via package extensions when those
packages are loaded). Variables prefixed with `$` are converted to each
module's geometry type with `GeoInterface.convert`:

```julia
using GeometryOpsTestHelpers  # also import LibGEOS/ArchGDAL/etc. to activate their extensions

poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]])

@test_implementations GO.area($poly) == 0.5

@testset_implementations "Area" begin
    @test GO.area($poly) == 0.5
end
```

Prefer these macros when testing public API functions, so behavior is verified
across all GeoInterface-compatible geometry libraries.

### Git commit style

- **Use imperative mood**: "Fix bug" not "Fixed bug" or "Fixes bug"
- **Start with a capital letter**: "Add feature" not "add feature"
- **Be concise but descriptive**: Explain what changed, not why (unless
  non-obvious)
- **No trailing periods**: Commit messages don't end with a period
- **Use backticks for code**: Reference functions/types like `smooth` or
  `TraitTarget`
- **No conventional commit prefixes**: Don't use "feat:", "fix:", "docs:", etc.

Examples from the project:

```
Fix type constraint in _smooth function
Add a `smooth` function
Refactor tests to be a bit easier to parse
Tree based acceleration for polygon clipping / boolean ops
Bump version from 0.1.30 to 0.1.31
```
