# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Testing
Run all tests:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Note: There is currently no way to run a single test file in isolation.

### Git Commit Style
Commit messages in this repository follow a simple, descriptive style:

- **Use imperative mood**: "Fix bug" not "Fixed bug" or "Fixes bug"
- **Start with a capital letter**: "Add feature" not "add feature"
- **Be concise but descriptive**: Explain what changed, not why (unless non-obvious)
- **No trailing periods**: Commit messages don't end with a period
- **Use backticks for code**: Reference functions/types with backticks like `smooth` or `TraitTarget`
- **No conventional commit prefixes**: Don't use "feat:", "fix:", "docs:", etc.

Examples from the project:
```
Fix type constraint in _smooth function
Add a `smooth` function 
Refactor tests to be a bit easier to parse
Tree based acceleration for polygon clipping / boolean ops
Bump version from 0.1.30 to 0.1.31
```

## High-Level Architecture

### Monorepo Structure
GeometryOps uses a monorepo structure with GeometryOpsCore as a subpackage:
- **GeometryOpsCore/**: Core abstractions, types, and primitive functions (`apply`, `applyreduce`, `flatten`, etc.)
- **src/**: Main package implementation
- **ext/**: Package extensions for optional dependencies (LibGEOS, Proj, TGGeometry, DataFrames, FlexiJoins)

### Core Abstractions

**GeoInterface Integration**: All functions work with any GeoInterface-compatible geometry. Dispatch is based on GeoInterface traits (PointTrait, LineStringTrait, PolygonTrait, etc.), not concrete types.

**Manifold System**: Operations can be performed on different manifolds:
- `Planar()`: Euclidean/Cartesian coordinates (default)
- `Spherical()`: Spherical coordinates on a unit sphere
- `Geodesic()`: Geodesic calculations on Earth (requires Proj extension)
- `AutoManifold()`: Automatically select appropriate manifold

Functions typically accept a manifold as the first argument:
```julia
area(Planar(), polygon)
area(Spherical(), polygon)
```

**Apply Framework**: The `apply` and `applyreduce` functions from GeometryOpsCore are the workhorses for geometry operations:
- `apply`: Applies a function to geometries matching a target trait, then reconstructs the geometry
- `applyreduce`: Applies a function and reduces the results (e.g., sum, min, max)
- `TraitTarget`: Specifies which geometry traits to target (e.g., `TraitTarget{GI.PointTrait}()`)

Example pattern:
```julia
applyreduce(WithTrait((trait, g) -> _area(T, trait, g)), +, _AREA_TARGETS, geom; threaded, init=zero(T))
```

### Code Organization Principles

1. **Literate Programming**: Source files use literate programming with documentation and examples at the top, followed by implementation. Examples should include visual plots using Makie/CairoMakie when appropriate.

2. **One File, One Job**: Each file should handle one semantic concept (e.g., `area.jl`, `distance.jl`, `centroid.jl`). Common utilities can be extracted to separate files.

3. **Public vs Internal**:
   - Public functions: Exported, documented, promise API stability
   - Internal functions: Prefixed with `_`, not exported, may have comments instead of docstrings

4. **File Structure Pattern** (see `src/methods/area.jl` or `src/methods/distance.jl`):
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

### Directory Structure

- **src/methods/**: Geometric operations and predicates
  - Basic: `area.jl`, `centroid.jl`, `distance.jl`, `perimeter.jl`, `angles.jl`
  - Spatial relations: `geom_relations/` (contains, intersects, within, etc.)
  - Clipping: `clipping/` (intersection, union, difference, cut, coverage)
  - Other: `barycentric.jl`, `buffer.jl`, `convex_hull.jl`, `orientation.jl`, `polygonize.jl`

- **src/transformations/**: Geometry transformations
  - `simplify.jl`, `segmentize.jl`, `smooth.jl`, `flip.jl`, `transform.jl`
  - `reproject.jl`: Coordinate reference system transformations
  - `correction/`: Geometry correction utilities

- **src/utils/**: Utility modules
  - `LoopStateMachine/`: State machine for polygon processing
  - `SpatialTreeInterface/`: Spatial indexing (STRtree)
  - `UnitSpherical/`: Spherical geometry utilities
  - `NaturalIndexing.jl`: Natural indexing utilities

- **test/**: Mirror structure of src/ with corresponding test files

## Writing a New Algorithm

### Step-by-Step Pattern

1. **Choose the right location**:
   - Methods (area, distance, etc.) go in `src/methods/`
   - Transformations (simplify, flip, etc.) go in `src/transformations/`

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

# Create example geometry and demonstrate usage
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

4. **Implement public API** with manifold support:
```julia
"""
    algorithm_name(geom, [T = Float64]; kwargs...)

[Docstring with description, parameters, return value]
"""
function algorithm_name(geom, ::Type{T} = Float64; threaded=false, kwargs...) where T
    algorithm_name(Planar(), geom, T; threaded, kwargs...)
end

# Manifold-specific implementations
function algorithm_name(m::Planar, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T
    # Use apply or applyreduce pattern
    applyreduce(WithTrait((trait, g) -> _algorithm(T, trait, g)), +, _ALGORITHM_TARGETS, geom; threaded, init=zero(T), kwargs...)
end

function algorithm_name(m::Spherical, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T
    # Spherical implementation if applicable
end
```

5. **Implement internal functions** (trait-dispatched):
```julia
# Dispatch on different geometry traits
_algorithm(::Type{T}, ::GI.PointTrait, point) where T = # Point implementation
_algorithm(::Type{T}, ::GI.LineStringTrait, linestring) where T = # LineString implementation
_algorithm(::Type{T}, ::GI.PolygonTrait, polygon) where T = # Polygon implementation
```

6. **Add to main module**: Include the file in `src/GeometryOps.jl`:
```julia
include("methods/algorithm_name.jl")
```

7. **Write tests**: Create `test/methods/algorithm_name.jl` mirroring the source structure:
```julia
using Test
using GeometryOps
import GeoInterface as GI

@testset "Algorithm Name" begin
    @testset "Point" begin
        # Test with points
    end

    @testset "Polygon" begin
        # Test with polygons
    end

    @testset "Edge cases" begin
        # Empty geometries, degenerate cases, etc.
    end
end
```

8. **Register test**: Add to `test/runtests.jl`:
```julia
@safetestset "Algorithm Name" begin include("methods/algorithm_name.jl") end
```

## Interfacing with Input

### GeoInterface Pattern
All functions accept any GeoInterface-compatible geometry. Always use GeoInterface accessors:
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

### Type Flexibility
Functions should accept a type parameter `T` (defaulting to `Float64`) for numeric calculations:
```julia
function my_function(geom, ::Type{T} = Float64) where T <: AbstractFloat
    result = zero(T)
    # ... computation using type T
    return result
end
```

### Threading Support
Most operations support threading via `threaded=false` keyword. Threading happens at the highest level (over arrays, feature collections, or multi-geometries):
```julia
applyreduce(f, op, targets, geom; threaded=true)  # Enable threading
```

### Extension-Based Optional Features
Optional dependencies are loaded via extensions:
- **GEOS**: `GEOS()` algorithm for calling LibGEOS functions
- **Proj**: `PROJ()` algorithm for reprojection and geodesic operations
- **TGGeometry**: `TG()` algorithm for fast C-based predicates
- **DataFrames**: `apply` works on DataFrame columns
- **FlexiJoins**: Spatial join operations

Use algorithm types to dispatch to these implementations. **The algorithm is always the first argument:**
```julia
# Use GEOS implementation (algorithm comes first)
buffer(GEOS(), geom, 10.0)

# Use native Julia implementation (default)
buffer(geom, 10.0)
```

The `Algorithm` type and its subtypes are defined in `GeometryOpsCore/src/types/algorithm.jl`.

### Working with Tables
`apply` and `applyreduce` work with any Tables.jl-compatible table (DataFrames, etc.):
```julia
# Apply operation to geometry column in a table
result = apply(PointTrait(), df.geometry) do point
    (GI.y(point), GI.x(point))  # Flip coordinates
end
```

## Common Patterns

### Error Handling for Missing Extensions
Functions requiring extensions should provide helpful error hints:
```julia
function __init__()
    Base.Experimental.register_error_hint(_reproject_error_hinter, MethodError)
    Base.Experimental.register_error_hint(_geodesic_segments_error_hinter, MethodError)
    Base.Experimental.register_error_hint(_buffer_error_hinter, MethodError)
end
```

### TraitTarget Usage
Use `TraitTarget` to specify which geometry types your function operates on:
```julia
# Single trait
const _TARGETS = TraitTarget{GI.PointTrait}()

# Union of traits
const _TARGETS = TraitTarget{Union{GI.PolygonTrait,GI.MultiPolygonTrait}}()

# More complex unions
const _TARGETS = TraitTarget{Union{GI.PointTrait,GI.LineStringTrait,GI.LinearRingTrait}}()
```

### WithTrait Pattern
Use `WithTrait` wrapper when you need both the trait and geometry in your function:
```julia
applyreduce(WithTrait((trait, geom) -> _my_function(trait, geom)), +, targets, input)
```

### Handling Empty Geometries
Always check for empty geometries when appropriate:
```julia
GI.isempty(geom) && return zero(T)  # Or appropriate default value
```
