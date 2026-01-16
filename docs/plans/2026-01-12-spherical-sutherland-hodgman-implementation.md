# Spherical Sutherland-Hodgman Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend `ConvexConvexSutherlandHodgman` to support `Spherical()` manifold for convex polygon intersection on the unit sphere.

**Architecture:** Add a new method dispatch on `ConvexConvexSutherlandHodgman{Spherical}` that uses spherical primitives (`spherical_orient`, `spherical_arc_intersection`) instead of planar ones. Key insight: on a sphere, inside/outside requires checking against ALL edges of the clip polygon, not just one.

**Tech Stack:** Julia, GeoInterface, UnitSpherical module (already in codebase)

---

## Task 1: Add `_tuple_point` for UnitSphericalPoint

**Files:**
- Modify: `src/utils/utils.jl:111-112` (after existing `_tuple_point` methods)

**Step 1: Add UnitSphericalPoint methods**

Add after line 112 in `src/utils/utils.jl`:

```julia
_tuple_point(p::UnitSpherical.UnitSphericalPoint{T}, ::Type{T}) where T = p
_tuple_point(p::UnitSpherical.UnitSphericalPoint, ::Type{T}) where T = UnitSpherical.UnitSphericalPoint{T}(p)
```

**Step 2: Verify it compiles**

Run: `julia --project=. -e 'using GeometryOps'`
Expected: No errors

**Step 3: Commit**

```bash
git add src/utils/utils.jl
git commit -m "Add _tuple_point methods for UnitSphericalPoint"
```

---

## Task 2: Add `_point_in_convex_spherical_polygon` with Test

**Files:**
- Modify: `src/methods/clipping/sutherland_hodgman.jl` (add after line 154)
- Modify: `test/methods/clipping/sutherland_hodgman.jl` (add new testset)

**Step 1: Write the failing test**

Add to `test/methods/clipping/sutherland_hodgman.jl` before the final `end`:

```julia
@testset "Spherical helpers" begin
    using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

    @testset "_point_in_convex_spherical_polygon" begin
        transform = UnitSphereFromGeographic()

        # CCW square near equator
        square_pts = UnitSphericalPoint{Float64}[
            transform((0.0, 0.0)),
            transform((2.0, 0.0)),
            transform((2.0, 2.0)),
            transform((0.0, 2.0))
        ]

        inside_pt = transform((1.0, 1.0))
        outside_pt = transform((5.0, 5.0))
        edge_pt = transform((1.0, 0.0))  # On edge

        @test GO._point_in_convex_spherical_polygon(inside_pt, square_pts) == true
        @test GO._point_in_convex_spherical_polygon(outside_pt, square_pts) == false
        # Edge point should be considered inside (>= 0 check)
        @test GO._point_in_convex_spherical_polygon(edge_pt, square_pts) == true
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: Error - `_point_in_convex_spherical_polygon` not defined

**Step 3: Write implementation**

Add to `src/methods/clipping/sutherland_hodgman.jl` before the final fallback function:

```julia
# Point in convex spherical polygon - true if point is on the left of all edges
function _point_in_convex_spherical_polygon(
    point::UnitSpherical.UnitSphericalPoint,
    polygon_points::Vector{<:UnitSpherical.UnitSphericalPoint}
)
    n = length(polygon_points)
    for i in 1:n
        edge_start = polygon_points[i]
        edge_end = polygon_points[mod1(i + 1, n)]
        if UnitSpherical.spherical_orient(edge_start, edge_end, point) < 0
            return false
        end
    end
    return true
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

**Step 5: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl test/methods/clipping/sutherland_hodgman.jl
git commit -m "Add _point_in_convex_spherical_polygon helper"
```

---

## Task 3: Add `_sh_spherical_intersection` Helper

**Files:**
- Modify: `src/methods/clipping/sutherland_hodgman.jl`

**Step 1: Add implementation**

Add after `_point_in_convex_spherical_polygon` in `src/methods/clipping/sutherland_hodgman.jl`:

```julia
# Compute intersection point of two great circle arcs
function _sh_spherical_intersection(
    p1::UnitSpherical.UnitSphericalPoint,
    p2::UnitSpherical.UnitSphericalPoint,
    p3::UnitSpherical.UnitSphericalPoint,
    p4::UnitSpherical.UnitSphericalPoint,
    ::Type{T}
) where T
    result = UnitSpherical.spherical_arc_intersection(p1, p2, p3, p4)

    # Simple extraction (option A)
    # TODO: May need to handle arc_overlap/arc_hinge cases if edge cases arise
    if !isempty(result.points)
        return UnitSpherical.UnitSphericalPoint{T}(result.points[1])
    end

    return UnitSpherical.UnitSphericalPoint{T}(p1)  # Fallback
end
```

**Step 2: Verify it compiles**

Run: `julia --project=. -e 'using GeometryOps'`
Expected: No errors

**Step 3: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl
git commit -m "Add _sh_spherical_intersection helper"
```

---

## Task 4: Add `_sh_clip_to_edge_spherical` Helper

**Files:**
- Modify: `src/methods/clipping/sutherland_hodgman.jl`

**Step 1: Add implementation**

Add after `_sh_spherical_intersection` in `src/methods/clipping/sutherland_hodgman.jl`:

```julia
# Clip polygon against a single edge using Sutherland-Hodgman rules (spherical version)
function _sh_clip_to_edge_spherical(
    polygon_points::Vector{UnitSpherical.UnitSphericalPoint{T}},
    edge_start::UnitSpherical.UnitSphericalPoint,
    edge_end::UnitSpherical.UnitSphericalPoint,
    clip_polygon_points::Vector{UnitSpherical.UnitSphericalPoint{T}},
    ::Type{T}
) where T
    output = UnitSpherical.UnitSphericalPoint{T}[]
    n = length(polygon_points)
    n == 0 && return output

    # Compute for first point, then carry forward (avoid duplicate computations)
    current_inside = _point_in_convex_spherical_polygon(polygon_points[1], clip_polygon_points)

    for i in 1:n
        current = polygon_points[i]
        next_idx = mod1(i + 1, n)
        next_pt = polygon_points[next_idx]

        next_inside = _point_in_convex_spherical_polygon(next_pt, clip_polygon_points)

        if current_inside
            push!(output, current)
            if !next_inside
                intr_pt = _sh_spherical_intersection(current, next_pt, edge_start, edge_end, T)
                push!(output, intr_pt)
            end
        elseif next_inside
            intr_pt = _sh_spherical_intersection(current, next_pt, edge_start, edge_end, T)
            push!(output, intr_pt)
        end

        current_inside = next_inside
    end

    return output
end
```

**Step 2: Verify it compiles**

Run: `julia --project=. -e 'using GeometryOps'`
Expected: No errors

**Step 3: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl
git commit -m "Add _sh_clip_to_edge_spherical helper"
```

---

## Task 5: Add Spherical Entry Point with Basic Test

**Files:**
- Modify: `src/methods/clipping/sutherland_hodgman.jl`
- Modify: `test/methods/clipping/sutherland_hodgman.jl`

**Step 1: Write the failing test**

Add to `test/methods/clipping/sutherland_hodgman.jl`, inside the "Spherical helpers" testset or as a new testset:

```julia
@testset "Spherical ConvexConvexSutherlandHodgman" begin
    using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

    # Helper to create spherical polygon from lon/lat coordinates
    function spherical_polygon(coords)
        transform = UnitSphereFromGeographic()
        points = [transform((lon, lat)) for (lon, lat) in coords]
        push!(points, points[1])  # close ring
        return GI.Polygon([points])
    end

    @testset "Basic intersection - small region" begin
        # Two overlapping squares near equator
        square1 = spherical_polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        square2 = spherical_polygon([(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)])

        result = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            square1, square2
        )
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(GO.Spherical(), result) > 0
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: Error - no method matching for `Spherical` manifold

**Step 3: Write implementation**

Add after `_sh_clip_to_edge_spherical` in `src/methods/clipping/sutherland_hodgman.jl`:

```julia
# Spherical Polygon-Polygon intersection using Sutherland-Hodgman
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman{Spherical},
    ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b
) where {T}
    ring_a = GI.getexterior(poly_a)
    ring_b = GI.getexterior(poly_b)

    # Validate input is UnitSphericalPoint
    first_pt = GI.getpoint(ring_a, 1)
    if !(first_pt isa UnitSpherical.UnitSphericalPoint)
        throw(ArgumentError(
            "Spherical ConvexConvexSutherlandHodgman requires UnitSphericalPoint coordinates, " *
            "got $(typeof(first_pt))"
        ))
    end

    # Collect clip polygon points (excluding closing point)
    clip_points = UnitSpherical.UnitSphericalPoint{T}[]
    for point in GI.getpoint(ring_b)
        if !isempty(clip_points) && point ≈ clip_points[1]
            continue
        end
        push!(clip_points, UnitSpherical.UnitSphericalPoint{T}(point))
    end

    # Build initial output list from poly_a (excluding closing point)
    output = UnitSpherical.UnitSphericalPoint{T}[]
    for point in GI.getpoint(ring_a)
        if !isempty(output) && point ≈ output[1]
            continue
        end
        push!(output, UnitSpherical.UnitSphericalPoint{T}(point))
    end

    # Clip against each edge of poly_b
    n_clip = length(clip_points)
    for i in 1:n_clip
        isempty(output) && break
        edge_start = clip_points[i]
        edge_end = clip_points[mod1(i + 1, n_clip)]
        output = _sh_clip_to_edge_spherical(output, edge_start, edge_end, clip_points, T)
    end

    # Handle empty result - degenerate polygon at north pole
    if isempty(output)
        north_pole = UnitSpherical.UnitSphericalPoint{T}(0, 0, 1)
        return GI.Polygon([[north_pole, north_pole, north_pole]])
    end

    # Close the ring and return
    push!(output, output[1])
    return GI.Polygon([output])
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

**Step 5: Commit**

```bash
git add src/methods/clipping/sutherland_hodgman.jl test/methods/clipping/sutherland_hodgman.jl
git commit -m "Add Spherical Sutherland-Hodgman intersection"
```

---

## Task 6: Add Remaining Spherical Tests

**Files:**
- Modify: `test/methods/clipping/sutherland_hodgman.jl`

**Step 1: Add comprehensive tests**

Expand the "Spherical ConvexConvexSutherlandHodgman" testset:

```julia
@testset "Spherical ConvexConvexSutherlandHodgman" begin
    using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

    function spherical_polygon(coords)
        transform = UnitSphereFromGeographic()
        points = [transform((lon, lat)) for (lon, lat) in coords]
        push!(points, points[1])
        return GI.Polygon([points])
    end

    spherical_area(poly) = GO.area(GO.Spherical(), poly)

    @testset "Basic intersection - small region" begin
        square1 = spherical_polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        square2 = spherical_polygon([(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)])

        result = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            square1, square2
        )
        @test GI.trait(result) isa GI.PolygonTrait
        @test spherical_area(result) > 0
    end

    @testset "No intersection" begin
        square1 = spherical_polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
        square2 = spherical_polygon([(10.0, 10.0), (11.0, 10.0), (11.0, 11.0), (10.0, 11.0)])

        result = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            square1, square2
        )
        @test spherical_area(result) ≈ 0.0 atol=1e-10
    end

    @testset "One contains other" begin
        large = spherical_polygon([(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)])
        small = spherical_polygon([(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)])

        result = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            large, small
        )
        @test spherical_area(result) ≈ spherical_area(small) rtol=1e-3
    end

    @testset "Triangles" begin
        tri1 = spherical_polygon([(0.0, 0.0), (4.0, 0.0), (2.0, 4.0)])
        tri2 = spherical_polygon([(1.0, 1.0), (3.0, 1.0), (2.0, 3.0)])

        result = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            tri1, tri2
        )
        @test GI.trait(result) isa GI.PolygonTrait
        @test spherical_area(result) > 0
    end

    @testset "Near pole" begin
        tri1 = spherical_polygon([(0.0, 85.0), (120.0, 85.0), (240.0, 85.0)])
        tri2 = spherical_polygon([(60.0, 85.0), (180.0, 85.0), (300.0, 85.0)])

        result = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            tri1, tri2
        )
        @test GI.trait(result) isa GI.PolygonTrait
        @test spherical_area(result) > 0
    end

    @testset "Input validation" begin
        planar_poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])

        @test_throws ArgumentError GO.intersection(
            GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
            planar_poly, planar_poly
        )
    end
end
```

**Step 2: Run tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

**Step 3: Commit**

```bash
git add test/methods/clipping/sutherland_hodgman.jl
git commit -m "Add comprehensive Spherical Sutherland-Hodgman tests"
```

---

## Task 7: Final Verification

**Step 1: Run full test suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

**Step 2: Verify git status is clean**

Run: `git status`
Expected: Clean working directory

**Step 3: Review commits**

Run: `git log --oneline -10`
Expected: See all implementation commits
