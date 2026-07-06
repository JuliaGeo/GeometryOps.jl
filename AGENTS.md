# AGENTS.md

Guidance for AI agents (and other contributors) working in this repository.

**Read [ARCHITECTURE.md](ARCHITECTURE.md) before writing code.** It describes
the package structure, core abstractions, design principles, and the
step-by-step pattern for adding a new algorithm. This file is about how to
work: the approach, the checks, and the commands.

## Guidelines

Defaults, not commandments — take a real exception when warranted, and say
why. The common failure they guard against is adding machinery where none was
needed.

- When two pieces don't fit, prefer in order: fix the shared abstraction at
  its source → adapt the consumer at the call site → accept the mismatch
  (slow paths are fine; `prepare` is the escape hatch). A new adapter/index/
  composition layer is the last resort — if truly needed, flag it and its
  tradeoff prominently.
- Generic means generic: a design should survive the hardest case (an opaque
  C-library tree, a foreign geometry with expensive accessors).
- Never transform user data silently; expensive or semantics-changing behavior
  is an explicit opt-in keyword.
- Build only what has a consumer in the same change; note future generality
  instead of building it.
- After generalizing, delete what it subsumed in the same change.
- Check new designs against conventions already in the file/codebase before
  presenting them.
- Report plainly: no scope decisions dressed as architecture, a one-line "why"
  on every deferred item, no coined jargon.
- Comments say what, not the story of why — exposition belongs in the literate
  header, docstrings, or commit messages.
- Read the relevant code before designing; run the affected tests (and
  benchmarks for perf-sensitive changes) before committing.
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
