# Sutherland-Hodgman Convex-Convex Clipping Design

## Overview

Implement the Sutherland-Hodgman polygon clipping algorithm as a new algorithm type for `GO.intersection`. This provides a simpler, faster alternative to Foster-Hormann for convex-convex polygon intersection.

## Decisions

| Decision | Choice |
|----------|--------|
| Algorithm name | `ConvexConvexSutherlandHodgman` |
| Input validation | Trust user, document requirement |
| Return type | Always `GI.Polygon` (empty if disjoint) |
| Integration | Dispatch on `intersection` via algorithm type |
| Manifold | Field included, only `Planar` implemented |
| Reused utilities | `eachedge`, `Predicates.orient`, `_intersection_point` |

## Algorithm Type & File Structure

**New file:** `src/methods/clipping/sutherland_hodgman.jl`

**Algorithm struct:**

```julia
"""
    ConvexConvexSutherlandHodgman{M <: Manifold} <: GeometryOpsCore.Algorithm{M}

Sutherland-Hodgman polygon clipping algorithm optimized for convex-convex intersection.

Both input polygons MUST be convex. If either polygon is non-convex, results are undefined.

This is simpler and faster than Foster-Hormann for small convex polygons, with O(n*m)
complexity where n and m are vertex counts.
"""
struct ConvexConvexSutherlandHodgman{M <: Manifold} <: GeometryOpsCore.Algorithm{M}
    manifold::M
end

# Default constructor uses Planar
ConvexConvexSutherlandHodgman() = ConvexConvexSutherlandHodgman(Planar())
```

**Module integration:** Include after `clipping_processor.jl` in `GeometryOps.jl`:

```julia
include("methods/clipping/sutherland_hodgman.jl")
```

**Export:** Add `ConvexConvexSutherlandHodgman` to the exports.

## Public API & Dispatch

**Method signatures on `intersection`:**

```julia
# Main entry point - algorithm dispatch
function intersection(
    alg::ConvexConvexSutherlandHodgman,
    geom_a,
    geom_b,
    ::Type{T}=Float64;
    kwargs...
) where {T<:AbstractFloat}
    return _intersection_sutherland_hodgman(
        alg, T,
        GI.trait(geom_a), geom_a,
        GI.trait(geom_b), geom_b
    )
end
```

**Trait dispatch (polygon-polygon only):**

```julia
# Polygon-Polygon intersection
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman{Planar},
    ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b
) where {T}
    # Implementation here
end

# Fallback for unsupported geometry combinations
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman,
    ::Type{T},
    trait_a, geom_a,
    trait_b, geom_b
) where {T}
    throw(ArgumentError(
        "ConvexConvexSutherlandHodgman only supports Polygon-Polygon intersection, " *
        "got $(typeof(trait_a)) and $(typeof(trait_b))"
    ))
end
```

No `target` kwarg needed since we always return a polygon.

## Core Algorithm Implementation

The Sutherland-Hodgman algorithm clips a subject polygon against each edge of the clip polygon sequentially. For each clip edge, it processes all vertices of the current output, producing a new vertex list.

**Core logic:**

```julia
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman{Planar},
    ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b
) where {T}
    # Get exterior rings (ignore holes - convex polygons don't have them)
    ring_a = GI.getexterior(poly_a)
    ring_b = GI.getexterior(poly_b)

    # Start with all vertices of poly_a as the output list
    output = collect_points(T, ring_a)

    # Clip against each edge of poly_b
    for (edge_start, edge_end) in eachedge(ring_b)
        isempty(output) && break
        output = _clip_to_edge(output, edge_start, edge_end, T)
    end

    # Return polygon (empty if no intersection)
    return GI.Polygon([output])
end
```

**Edge clipping helper** - for each edge of the clip polygon, iterate through output vertices and apply the four Sutherland-Hodgman cases:

1. Both points inside -> keep end point
2. Start inside, end outside -> keep intersection
3. Both outside -> keep nothing
4. Start outside, end inside -> keep intersection and end point

## Helper Functions

**Edge clipping using existing primitives:**

```julia
function _clip_to_edge(polygon_points, edge_start, edge_end, ::Type{T}) where T
    output = Tuple{T,T}[]

    for i in eachindex(polygon_points)
        current = polygon_points[i]
        next_idx = mod1(i + 1, length(polygon_points))
        next = polygon_points[next_idx]

        # orient returns positive if counter-clockwise (left), negative if clockwise (right)
        current_inside = Predicates.orient(edge_start, edge_end, current) >= 0
        next_inside = Predicates.orient(edge_start, edge_end, next) >= 0

        if current_inside
            push!(output, current)
            if !next_inside
                push!(output, _intersection_point(T, (current, next), (edge_start, edge_end)))
            end
        elseif next_inside
            push!(output, _intersection_point(T, (current, next), (edge_start, edge_end)))
        end
    end

    return output
end
```

No new geometric primitives needed - reuses `Predicates.orient` and `_intersection_point` from existing clipping code.

## Testing

**New file:** `test/methods/clipping/sutherland_hodgman.jl`

**Test cases:**

```julia
using Test
using GeometryOps
import GeoInterface as GI

@testset "ConvexConvexSutherlandHodgman" begin
    @testset "Basic intersection" begin
        # Two overlapping squares
        square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

        result = GO.intersection(ConvexConvexSutherlandHodgman(), square1, square2)
        @test GO.area(result) ≈ 1.0  # 1x1 overlap
    end

    @testset "No intersection" begin
        # Disjoint squares
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0), (5.0, 5.0)]])

        result = GO.intersection(ConvexConvexSutherlandHodgman(), square1, square2)
        @test GO.area(result) ≈ 0.0
    end

    @testset "One contains other" begin
        # Large square contains small square
        large = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
        small = GI.Polygon([[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]])

        result = GO.intersection(ConvexConvexSutherlandHodgman(), large, small)
        @test GO.area(result) ≈ 4.0
    end

    @testset "Triangles" begin
        # Two overlapping triangles
    end

    @testset "Edge cases" begin
        # Shared edge, shared vertex, identical polygons
    end
end
```

**Register in `test/runtests.jl`:**

```julia
@safetestset "Sutherland-Hodgman" begin include("methods/clipping/sutherland_hodgman.jl") end
```

## File Changes Summary

**Files to create:**

- `src/methods/clipping/sutherland_hodgman.jl` - Algorithm struct + implementation
- `test/methods/clipping/sutherland_hodgman.jl` - Tests

**Files to modify:**

- `src/GeometryOps.jl` - Add include and export
- `test/runtests.jl` - Add `@safetestset` for new tests

## Usage

```julia
result = GO.intersection(ConvexConvexSutherlandHodgman(), convex_poly_a, convex_poly_b)
```
