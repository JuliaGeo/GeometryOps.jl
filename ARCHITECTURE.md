# Architecture

How GeometryOps is put together, and the principles its code is expected to
follow. For workflow guidance and development commands, see
[AGENTS.md](AGENTS.md).

## Monorepo structure

GeometryOps is organized as a Julia workspace (see `[workspace]` in
`Project.toml`; subpackages are resolved by path via `[sources]`):

- **GeometryOpsCore/**: Core abstractions, types, and primitive functions
  (`apply`, `applyreduce`, `flatten`, manifolds, algorithm types)
- **GeometryOpsTestHelpers/**: Test utilities (`@test_implementations`,
  `@testset_implementations`) for running tests across multiple geometry
  implementations
- **src/**: Main package implementation
- **ext/**: Package extensions for optional dependencies (LibGEOS, Proj,
  TGGeometry, DataFrames, FlexiJoins, Makie)

## Core abstractions

**GeoInterface integration**: All functions work with any
GeoInterface-compatible geometry. Dispatch is based on GeoInterface traits
(`PointTrait`, `LineStringTrait`, `PolygonTrait`, ...), not concrete types.

**Manifold system**: Operations can be performed on different manifolds, passed
as the first argument (`area(Spherical(), polygon)`):

- `Planar()`: Euclidean/Cartesian coordinates (default)
- `Spherical()`: Spherical coordinates on a unit sphere
- `Geodesic()`: Geodesic calculations on Earth (requires Proj extension)
- `AutoManifold()`: Automatically select an appropriate manifold

**Apply framework**: `apply` and `applyreduce` from GeometryOpsCore are the
workhorses for geometry operations:

- `apply`: Applies a function to geometries matching a target trait, then
  reconstructs the geometry
- `applyreduce`: Applies a function and reduces the results (e.g. sum, min)
- `TraitTarget`: Specifies which geometry traits to target

```julia
applyreduce(WithTrait((trait, g) -> _area(T, trait, g)), +, _AREA_TARGETS, geom; threaded, init=zero(T))
```

**Algorithm types**: Optional backends dispatch through `Algorithm` subtypes
(defined in `GeometryOpsCore/src/types/algorithm.jl`), always as the first
argument: `buffer(GEOS(), geom, 10.0)` vs the native `buffer(geom, 10.0)`.
Extensions provide `GEOS()` (LibGEOS), `PROJ()` (reprojection/geodesic), and
`TG()` (fast C predicates).

## Prepared geometry and spatial indexing

- `prepare` / `Prepared` (`src/prepared.jl`): materializes a geometry into the
  native layout (preserving coordinate number types) and attaches preparations
  such as `EdgeTree`. Preparations survive GeoInterface decomposition: a ring
  pulled out of a prepared polygon still carries its prep. Query the object you
  hold (`getprep`, `hasprep`) — never a parent-held parallel array.
- `SpatialTreeInterface` (`src/utils/SpatialTreeInterface/`): generic
  single/dual tree traversal. All tree consumers program against this
  interface; none may assume a concrete backend or its traversal order.
- Tree backends: `NaturalIndexing.NaturalIndex` (default),
  `FlexibleRTrees` (STR/HPR/unsorted bulk loading), `STRtree`. Any
  SpatialTreeInterface-compatible tree works, including opaque external ones.
- Clipping accelerators (`src/methods/clipping/clipping_processor.jl`): the
  `IntersectionAccelerator` family — `NestedLoop`, `TreeAccelerator` with
  per-side `TreePolicy` (`IterateEdges`, `BuildTree`), and `AutoAccelerator`.
  Tree sides reuse a `Prepared` input's trees and index anything else as it is;
  `prepare = true` opts into ephemeral preparation.

## Design principles

Deliberate choices the codebase tries to stay consistent with. Deviating is
sometimes right — but do it knowingly, with the reason stated.

1. **Genericity lives in the shared interface, not in special cases.** If a
   generic path is slow for one backend, fix the generic implementation so
   every consumer benefits. A type-specialized bypass forks maintenance and
   must survive the hardest test: "what if the backend is an opaque C-library
   tree?"
2. **Explicit over implicit.** Never silently transform user-supplied data
   (copying, re-materializing, closing rings, converting number types).
   Expensive or semantics-changing conveniences are opt-in keyword arguments.
   APIs read literally at the call site — e.g. curried constructors like
   `EdgeTree(STRtree)` rather than implicit pairing syntax.
3. **Slow paths are acceptable; escape hatches already exist.** Arbitrary
   GeoInterface geometry with a poor memory layout is allowed to be slow — the
   answer is `prepare`, not an adapter layer in the hot loop.
4. **Normalize once at the boundary.** Shape, layout, and closedness concerns
   are resolved at ingestion (`prepare`, entry points); downstream code assumes
   the canonical form and does not re-check defensively per call site.
5. **Dispatch families get an abstract supertype** (`Algorithm`, `Manifold`,
   `IntersectionAccelerator`, `TreePolicy`), so signatures can constrain their
   parameters and the family is discoverable.
6. **No speculative machinery.** Every abstraction, trait, and type parameter
   must have a consumer in the same change. "Forward-looking" unwired
   components don't ship.

## Directory structure

- **src/methods/**: Geometric operations and predicates
  - Basic: `area.jl`, `centroid.jl`, `distance.jl`, `perimeter.jl`, `angles.jl`
  - Spatial relations: `geom_relations/` (contains, intersects, within, ...)
  - Clipping: `clipping/` (intersection, union, difference, cut, coverage)
  - Other: `barycentric.jl`, `buffer.jl`, `convex_hull.jl`, `orientation.jl`,
    `polygonize.jl`
- **src/transformations/**: Geometry transformations
  - `simplify.jl`, `segmentize.jl`, `smooth.jl`, `flip.jl`, `transform.jl`
  - `reproject.jl`: Coordinate reference system transformations
  - `correction/`: Geometry correction utilities
- **src/utils/**: Utility modules
  - `SpatialTreeInterface/`: Generic spatial tree traversal
  - `NaturalIndexing.jl`: Natural indexing (default edge-tree backend)
  - `FlexibleRTrees/`: Bulk-loaded R-trees (STR, HPR, unsorted)
  - `LoopStateMachine/`: State machine for polygon processing
  - `UnitSpherical/`: Spherical geometry utilities
- **test/**: Mirrors the structure of src/ with corresponding test files

## Code organization principles

1. **Literate programming**: Source files put documentation and examples at the
   top, followed by implementation. Examples include visual plots using
   Makie/CairoMakie when appropriate.
2. **One file, one job**: Each file handles one semantic concept (`area.jl`,
   `distance.jl`, ...). Common utilities are extracted to separate files.
3. **Public vs internal**:
   - Public functions: exported, documented, promise API stability
   - Internal functions: prefixed with `_`, not exported, may have comments
     instead of docstrings
4. **File structure pattern** (see `src/methods/area.jl` or
   `src/methods/distance.jl`):

   ```julia
   # # Title
   export function_name

   #=
   ## What is [concept]?
   [Explanation with examples, plots]

   ## Implementation
   [Implementation notes]
   =#

   # Public API with docstring
   function_name(geom, args...) = ...

   # Internal implementation functions
   _function_name(...) = ...
   ```

## Writing a new algorithm

1. **Choose the location**: methods go in `src/methods/`, transformations in
   `src/transformations/`.

2. **Create the file with literate documentation**:

```julia
# # Algorithm Name
export algorithm_name

#=
## What is [algorithm]?
[Clear explanation with visual examples using Makie]

```@example demo
import GeometryOps as GO
import GeoInterface as GI
using CairoMakie

polygon = GI.Polygon([[(0,0), (1,0), (1,1), (0,1), (0,0)]])
result = GO.algorithm_name(polygon)
```

## Implementation
[Notes about the approach, algorithm complexity, special cases]
=#
```

3. **Define target traits** (if using apply/applyreduce):

```julia
const _ALGORITHM_TARGETS = TraitTarget{Union{GI.PolygonTrait,GI.MultiPolygonTrait}}()
```

4. **Implement the public API** with manifold support:

```julia
"""
    algorithm_name(geom, [T = Float64]; kwargs...)

[Docstring with description, parameters, return value]
"""
function algorithm_name(geom, ::Type{T} = Float64; threaded=false, kwargs...) where T
    algorithm_name(Planar(), geom, T; threaded, kwargs...)
end

function algorithm_name(m::Planar, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T
    applyreduce(WithTrait((trait, g) -> _algorithm(T, trait, g)), +, _ALGORITHM_TARGETS, geom; threaded, init=zero(T), kwargs...)
end

function algorithm_name(m::Spherical, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T
    # Spherical implementation if applicable
end
```

5. **Implement internal functions**, dispatched on geometry traits:

```julia
_algorithm(::Type{T}, ::GI.PointTrait, point) where T = ...
_algorithm(::Type{T}, ::GI.LineStringTrait, linestring) where T = ...
_algorithm(::Type{T}, ::GI.PolygonTrait, polygon) where T = ...
```

6. **Add to the main module** in `src/GeometryOps.jl`:

```julia
include("methods/algorithm_name.jl")
```

7. **Write tests** in `test/methods/algorithm_name.jl`, mirroring the source
   structure. Make the file self-contained (it runs in a `@safetestset` and
   should also work via a bare `include`), and use the GeometryOpsTestHelpers
   macros to cover all geometry implementations:

```julia
using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG  # loading these activates the corresponding TestHelpers extensions
import ArchGDAL as AG
using GeometryOpsTestHelpers

poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]])

@testset_implementations "Polygon" begin
    @test GO.algorithm_name($poly) == expected
end

@testset "Edge cases" begin
    # Empty geometries, degenerate cases, etc.
end
```

8. **Register the test** in `test/runtests.jl`:

```julia
@safetestset "Algorithm Name" begin include("methods/algorithm_name.jl") end
```

## Interfacing with input

### GeoInterface pattern

All functions accept any GeoInterface-compatible geometry. Always use
GeoInterface accessors:

```julia
GI.trait(geom)           # Get geometry trait
GI.npoint(geom)          # Number of points
GI.getpoint(geom, i)     # Get point at index i
GI.getpoint(geom)        # Iterator over points
GI.x(point)              # X coordinate
GI.y(point)              # Y coordinate
GI.z(point)              # Z coordinate (if exists)
GI.getexterior(poly)     # Exterior ring
GI.gethole(poly)         # Iterator over holes
GI.isempty(geom)         # Check if empty
```

### Type flexibility

Functions accept a type parameter `T` (defaulting to `Float64`) for numeric
calculations:

```julia
function my_function(geom, ::Type{T} = Float64) where T <: AbstractFloat
    result = zero(T)
    # ... computation using type T
    return result
end
```

### Threading support

Most operations support threading via the `threaded=false` keyword. Threading
happens at the highest level (over arrays, feature collections, or
multi-geometries):

```julia
applyreduce(f, op, targets, geom; threaded=true)
```

### Working with tables

`apply` and `applyreduce` work with any Tables.jl-compatible table:

```julia
result = apply(PointTrait(), df.geometry) do point
    (GI.y(point), GI.x(point))  # Flip coordinates
end
```

## Common patterns

### Error handling for missing extensions

Functions requiring extensions provide error hints:

```julia
function __init__()
    Base.Experimental.register_error_hint(_reproject_error_hinter, MethodError)
    Base.Experimental.register_error_hint(_geodesic_segments_error_hinter, MethodError)
    Base.Experimental.register_error_hint(_buffer_error_hinter, MethodError)
end
```

### TraitTarget usage

```julia
const _TARGETS = TraitTarget{GI.PointTrait}()
const _TARGETS = TraitTarget{Union{GI.PolygonTrait,GI.MultiPolygonTrait}}()
const _TARGETS = TraitTarget{Union{GI.PointTrait,GI.LineStringTrait,GI.LinearRingTrait}}()
```

### WithTrait pattern

Use the `WithTrait` wrapper when you need both the trait and geometry:

```julia
applyreduce(WithTrait((trait, geom) -> _my_function(trait, geom)), +, targets, input)
```

### Handling empty geometries

```julia
GI.isempty(geom) && return zero(T)  # Or appropriate default value
```
