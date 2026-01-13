# Spherical Sutherland-Hodgman Extension Design

## Overview

Extend `ConvexConvexSutherlandHodgman` to support `Spherical()` manifold for convex polygon intersection on the unit sphere.

## Decisions

| Decision | Choice |
|----------|--------|
| Input format | `UnitSphericalPoint` (validated at runtime) |
| Output format | `UnitSphericalPoint` (same as input) |
| Empty result | Degenerate polygon at north pole |
| Arc intersection | Simple extraction from `ArcIntersectionResult` (option A) |
| Inside/outside test | Full convex point-in-polygon (left of all edges) |
| Type parameter | Explicit `T` parameter (same as Planar) |
| Edge iteration | Refactor `_tuple_point` for `UnitSphericalPoint` |

## Key Insight: Spherical Inside/Outside

On a sphere, orientation relative to a single edge doesn't determine inside/outside. For convex spherical polygons, a point is inside if and only if it's on the left side of ALL edges. This requires passing the full clip polygon to the edge clipping function.

## Implementation

### Entry Point

```julia
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

    # Handle empty result
    if isempty(output)
        north_pole = UnitSpherical.UnitSphericalPoint{T}(0, 0, 1)
        return GI.Polygon([[north_pole, north_pole, north_pole]])
    end

    # Close the ring and return
    push!(output, output[1])
    return GI.Polygon([output])
end
```

### Point in Convex Spherical Polygon

```julia
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

### Edge Clipping (Optimized)

```julia
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

### Spherical Intersection Helper

```julia
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

    return p1  # Fallback
end
```

### `_tuple_point` Refactor

In `src/utils/utils.jl`:

```julia
_tuple_point(p::UnitSpherical.UnitSphericalPoint{T}, ::Type{T}) where T = p
_tuple_point(p::UnitSpherical.UnitSphericalPoint, ::Type{T}) where T = UnitSpherical.UnitSphericalPoint{T}(p)
```

## File Changes

**Modify:**
- `src/methods/clipping/sutherland_hodgman.jl` - Add Spherical implementation
- `src/utils/utils.jl` - Add `_tuple_point` for `UnitSphericalPoint`
- `test/methods/clipping/sutherland_hodgman.jl` - Add Spherical tests

## Testing

### Mirrored Planar Tests
- Basic intersection (overlapping squares near equator)
- No intersection (disjoint squares)
- One contains other
- Triangles

### Spherical-Specific Tests
- Point in convex spherical polygon helper
- Near pole polygons
- Crossing antimeridian
- Input validation (non-UnitSphericalPoint should error)

## Usage

```julia
using GeometryOps.UnitSpherical: UnitSphereFromGeographic

transform = UnitSphereFromGeographic()
poly_a = GI.Polygon([[transform((lon, lat)) for (lon, lat) in coords_a]])
poly_b = GI.Polygon([[transform((lon, lat)) for (lon, lat) in coords_b]])

result = GO.intersection(
    ConvexConvexSutherlandHodgman(Spherical()),
    poly_a,
    poly_b
)
```
