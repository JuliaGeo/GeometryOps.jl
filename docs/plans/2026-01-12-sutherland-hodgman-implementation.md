# Sutherland-Hodgman Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement convex-convex polygon intersection using the Sutherland-Hodgman algorithm.

**Architecture:** New algorithm type `ConvexConvexSutherlandHodgman` that dispatches through the existing `intersection` function. Clips subject polygon against each edge of the clip polygon sequentially.

**Tech Stack:** Julia, GeoInterface, existing GeometryOps utilities (`eachedge`, `Predicates.orient`, `_intersection_point`)

---

### Task 1: Create Algorithm Struct and Export

**Files:**
- Create: `src/methods/clipping/sutherland_hodgman.jl`
- Modify: `src/GeometryOps.jl:76` (add include after `union.jl`)

**Step 1: Create the new file with algorithm struct**

Create `src/methods/clipping/sutherland_hodgman.jl`:

```julia
# # Sutherland-Hodgman Convex-Convex Clipping
export ConvexConvexSutherlandHodgman

"""
    ConvexConvexSutherlandHodgman{M <: Manifold} <: GeometryOpsCore.Algorithm{M}

Sutherland-Hodgman polygon clipping algorithm optimized for convex-convex intersection.

Both input polygons MUST be convex. If either polygon is non-convex, results are undefined.

This is simpler and faster than Foster-Hormann for small convex polygons, with O(n*m)
complexity where n and m are vertex counts.

# Example

```julia
import GeometryOps as GO, GeoInterface as GI

square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
```
"""
struct ConvexConvexSutherlandHodgman{M <: Manifold} <: GeometryOpsCore.Algorithm{M}
    manifold::M
end

# Default constructor uses Planar
ConvexConvexSutherlandHodgman() = ConvexConvexSutherlandHodgman(Planar())
```

**Step 2: Add include to GeometryOps.jl**

In `src/GeometryOps.jl`, after line 76 (`include("methods/clipping/union.jl")`), add:

```julia
include("methods/clipping/sutherland_hodgman.jl")
```

**Step 3: Verify the struct loads**

Run:
```bash
julia --project=. -e 'using GeometryOps; println(GeometryOps.ConvexConvexSutherlandHodgman())'
```

Expected output: `ConvexConvexSutherlandHodgman{Planar}(Planar())`

**Step 4: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl src/GeometryOps.jl
git commit -m "Add ConvexConvexSutherlandHodgman algorithm struct"
```

---

### Task 2: Write Failing Test for Basic Intersection

**Files:**
- Create: `test/methods/clipping/sutherland_hodgman.jl`
- Modify: `test/runtests.jl:40` (add safetestset after Polygon Clipping)

**Step 1: Create test file with basic intersection test**

Create `test/methods/clipping/sutherland_hodgman.jl`:

```julia
using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "ConvexConvexSutherlandHodgman" begin
    @testset "Basic intersection" begin
        # Two overlapping squares - intersection is 1x1 square
        square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 1.0 atol=1e-10
    end
end
```

**Step 2: Register test in runtests.jl**

In `test/runtests.jl`, after line 40 (`@safetestset "Polygon Clipping" ...`), add:

```julia
@safetestset "Sutherland-Hodgman" begin include("methods/clipping/sutherland_hodgman.jl") end
```

**Step 3: Run test to verify it fails**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | grep -A 5 "Sutherland-Hodgman"
```

Expected: FAIL with `MethodError: no method matching intersection(::ConvexConvexSutherlandHodgman{Planar}, ...)`

**Step 4: Commit**

```bash
git add test/methods/clipping/sutherland_hodgman.jl test/runtests.jl
git commit -m "Add failing test for Sutherland-Hodgman intersection"
```

---

### Task 3: Implement Public API Dispatch

**Files:**
- Modify: `src/methods/clipping/sutherland_hodgman.jl`

**Step 1: Add intersection method dispatch**

Append to `src/methods/clipping/sutherland_hodgman.jl`:

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

**Step 2: Run test to verify dispatch works but implementation missing**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | grep -A 5 "Sutherland-Hodgman"
```

Expected: FAIL with `ArgumentError: ConvexConvexSutherlandHodgman only supports Polygon-Polygon intersection`

**Step 3: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl
git commit -m "Add intersection dispatch for ConvexConvexSutherlandHodgman"
```

---

### Task 4: Implement Core Algorithm

**Files:**
- Modify: `src/methods/clipping/sutherland_hodgman.jl`

**Step 1: Add the Polygon-Polygon implementation**

Append to `src/methods/clipping/sutherland_hodgman.jl` (before the fallback method):

```julia
# Polygon-Polygon intersection using Sutherland-Hodgman
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman{Planar},
    ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b
) where {T}
    # Get exterior rings (convex polygons have no holes)
    ring_a = GI.getexterior(poly_a)
    ring_b = GI.getexterior(poly_b)

    # Start with vertices of poly_a as the output list (excluding closing point)
    output = Tuple{T,T}[]
    for point in GI.getpoint(ring_a)
        pt = _tuple_point(point, T)
        # Skip the closing point (same as first)
        if !isempty(output) && pt == output[1]
            continue
        end
        push!(output, pt)
    end

    # Clip against each edge of poly_b
    for (edge_start, edge_end) in eachedge(ring_b, T)
        isempty(output) && break
        output = _sh_clip_to_edge(output, edge_start, edge_end, T)
    end

    # Close the ring if we have points
    if !isempty(output)
        push!(output, output[1])
    end

    # Return polygon
    return GI.Polygon([output])
end

# Clip polygon against a single edge using Sutherland-Hodgman rules
function _sh_clip_to_edge(polygon_points::Vector{Tuple{T,T}}, edge_start, edge_end, ::Type{T}) where T
    output = Tuple{T,T}[]
    n = length(polygon_points)
    n == 0 && return output

    for i in 1:n
        current = polygon_points[i]
        next_pt = polygon_points[mod1(i + 1, n)]

        # Determine if points are inside (left of or on the edge)
        # orient > 0 means left (inside for CCW polygon), == 0 means on edge, < 0 means right (outside)
        current_inside = Predicates.orient(edge_start, edge_end, current; exact=False()) >= 0
        next_inside = Predicates.orient(edge_start, edge_end, next_pt; exact=False()) >= 0

        if current_inside
            push!(output, current)
            if !next_inside
                # Exiting: add intersection point
                intr_pt = _sh_line_intersection(current, next_pt, edge_start, edge_end, T)
                push!(output, intr_pt)
            end
        elseif next_inside
            # Entering: add intersection point
            intr_pt = _sh_line_intersection(current, next_pt, edge_start, edge_end, T)
            push!(output, intr_pt)
        end
        # Both outside: add nothing
    end

    return output
end

# Compute intersection point of line segment (p1, p2) with line through (p3, p4)
function _sh_line_intersection(p1::Tuple{T,T}, p2::Tuple{T,T}, p3::Tuple{T,T}, p4::Tuple{T,T}, ::Type{T}) where T
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    # Lines are parallel - shouldn't happen in valid Sutherland-Hodgman usage
    if abs(denom) < eps(T)
        return p1  # Fallback
    end

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom

    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)

    return (T(x), T(y))
end
```

**Step 2: Run test to verify it passes**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | grep -A 10 "Sutherland-Hodgman"
```

Expected: PASS

**Step 3: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl
git commit -m "Implement Sutherland-Hodgman core algorithm"
```

---

### Task 5: Add More Test Cases

**Files:**
- Modify: `test/methods/clipping/sutherland_hodgman.jl`

**Step 1: Add comprehensive test cases**

Replace the test file with:

```julia
using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "ConvexConvexSutherlandHodgman" begin
    @testset "Basic intersection" begin
        # Two overlapping squares - intersection is 1x1 square
        square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 1.0 atol=1e-10
    end

    @testset "No intersection" begin
        # Disjoint squares
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0), (5.0, 5.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 0.0 atol=1e-10
    end

    @testset "One contains other" begin
        # Large square contains small square
        large = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
        small = GI.Polygon([[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), large, small)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 4.0 atol=1e-10

        # Reverse order should give same result
        result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), small, large)
        @test GO.area(result2) ≈ 4.0 atol=1e-10
    end

    @testset "Triangles" begin
        # Two overlapping triangles
        tri1 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (2.0, 4.0), (0.0, 0.0)]])
        tri2 = GI.Polygon([[(0.0, 2.0), (4.0, 2.0), (2.0, -2.0), (0.0, 2.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), tri1, tri2)
        @test result isa GI.Polygon
        @test GO.area(result) > 0
    end

    @testset "Identical polygons" begin
        # Same polygon should return itself
        square = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square, square)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 4.0 atol=1e-10
    end

    @testset "Shared edge" begin
        # Two squares sharing an edge
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        # Shared edge only - area should be 0 or near 0
        @test GO.area(result) ≈ 0.0 atol=1e-10
    end

    @testset "Unsupported geometry types" begin
        square = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        point = GI.Point(0.5, 0.5)

        @test_throws ArgumentError GO.intersection(GO.ConvexConvexSutherlandHodgman(), square, point)
        @test_throws ArgumentError GO.intersection(GO.ConvexConvexSutherlandHodgman(), point, square)
    end
end
```

**Step 2: Run all tests**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | grep -A 20 "Sutherland-Hodgman"
```

Expected: All tests PASS

**Step 3: Commit**

```bash
git add test/methods/clipping/sutherland_hodgman.jl
git commit -m "Add comprehensive tests for Sutherland-Hodgman"
```

---

### Task 6: Run Full Test Suite

**Files:** None (verification only)

**Step 1: Run full test suite**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: All tests pass

**Step 2: Commit if any fixes needed**

If fixes were needed, commit them:
```bash
git add -A
git commit -m "Fix issues found in full test suite"
```

---

## File Changes Summary

**Files created:**
- `src/methods/clipping/sutherland_hodgman.jl`
- `test/methods/clipping/sutherland_hodgman.jl`

**Files modified:**
- `src/GeometryOps.jl` - Add include
- `test/runtests.jl` - Add safetestset

## Final Usage

```julia
import GeometryOps as GO, GeoInterface as GI

square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
GO.area(result)  # ≈ 1.0
```
