# # Sutherland-Hodgman Convex-Convex Clipping
export ConvexConvexSutherlandHodgman

"""
    ConvexConvexSutherlandHodgman{M <: Manifold} <: GeometryOpsCore.Algorithm{M}

Sutherland-Hodgman polygon clipping algorithm optimized for convex-convex intersection.

Both input polygons MUST be convex. If either polygon is non-convex, results are undefined.

This is simpler and faster than Foster-Hormann for small convex polygons, with O(n*m)
complexity where n and m are vertex counts.

## Spherical manifold

For `Spherical()` manifold, input polygons must have **counter-clockwise winding** when
viewed from outside the sphere (i.e., the interior is on the left when traversing edges).
Polygons with clockwise winding will produce incorrect results (typically a degenerate
polygon). Use `GO.fix(geom; corrections=[GO.ClosedRing(), GO.GeometryCorrection()])` or
manually reverse the coordinates if your input has the wrong winding order.

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

    # Handle empty result (no intersection) - return degenerate polygon with zero area
    if isempty(output)
        zero_pt = (zero(T), zero(T))
        return GI.Polygon([[zero_pt, zero_pt, zero_pt]])
    end

    # Close the ring
    push!(output, output[1])

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

# Compute intersection of subject arc (p1, p2) with the GREAT CIRCLE through (p3, p4)
#
# Note: Sutherland-Hodgman clips against half-planes defined by clip edges. We need to find
# where the subject arc crosses the great circle (infinite line) containing the clip edge,
# NOT where two finite arcs intersect. This is critical because the subject arc may cross
# the clip edge's great circle extension without intersecting the finite clip edge itself.
function _sh_spherical_intersection(
    p1::UnitSpherical.UnitSphericalPoint,
    p2::UnitSpherical.UnitSphericalPoint,
    p3::UnitSpherical.UnitSphericalPoint,
    p4::UnitSpherical.UnitSphericalPoint,
    ::Type{T}
) where T
    tol = eps(T) * 16

    # Get great circle normals
    n_subject = UnitSpherical.robust_cross_product(p1, p2)
    n_clip = UnitSpherical.robust_cross_product(p3, p4)

    n_subject_norm = norm(n_subject)
    n_clip_norm = norm(n_clip)

    # Handle degenerate cases
    if n_subject_norm < tol || n_clip_norm < tol
        return UnitSpherical.UnitSphericalPoint{T}(p1)
    end

    n_subject = n_subject / n_subject_norm
    n_clip = n_clip / n_clip_norm

    # Intersection direction is the cross product of the normals
    intersection_dir = n_subject × n_clip
    dir_norm = norm(intersection_dir)

    # If normals are parallel, great circles are the same (collinear arcs)
    if dir_norm < tol
        # Return midpoint of subject arc as a reasonable fallback
        return UnitSpherical.UnitSphericalPoint{T}(p1)
    end

    intersection_dir = intersection_dir / dir_norm

    # Two candidate intersection points (antipodal)
    cand1 = UnitSpherical.UnitSphericalPoint{T}(intersection_dir)
    cand2 = UnitSpherical.UnitSphericalPoint{T}(-intersection_dir)

    # Return the candidate that lies on the subject arc
    if UnitSpherical.point_on_spherical_arc(cand1, p1, p2)
        return cand1
    elseif UnitSpherical.point_on_spherical_arc(cand2, p1, p2)
        return cand2
    end

    # Fallback: neither candidate is on the arc (shouldn't happen if orient detected a crossing)
    return UnitSpherical.UnitSphericalPoint{T}(p1)
end

# Clip polygon against a single edge using Sutherland-Hodgman rules (spherical version)
function _sh_clip_to_edge_spherical(
    polygon_points::Vector{UnitSpherical.UnitSphericalPoint{T}},
    edge_start::UnitSpherical.UnitSphericalPoint,
    edge_end::UnitSpherical.UnitSphericalPoint,
    ::Type{T}
) where T
    output = UnitSpherical.UnitSphericalPoint{T}[]
    n = length(polygon_points)
    n == 0 && return output

    # Track actual orient values to handle edge cases (orient=0 means exactly on the edge)
    current_orient = UnitSpherical.spherical_orient(edge_start, edge_end, polygon_points[1])

    for i in 1:n
        current = polygon_points[i]
        next_idx = mod1(i + 1, n)
        next_pt = polygon_points[next_idx]

        next_orient = UnitSpherical.spherical_orient(edge_start, edge_end, next_pt)
        current_inside = current_orient >= 0
        next_inside = next_orient >= 0

        if current_inside
            push!(output, current)
            if !next_inside
                # Exiting: add intersection point
                # If current is exactly on the edge (orient=0), it IS the intersection,
                # and we've already added it above - so don't add again
                if current_orient != 0
                    intr_pt = _sh_spherical_intersection(current, next_pt, edge_start, edge_end, T)
                    push!(output, intr_pt)
                end
            end
        elseif next_inside
            # Entering: add intersection point
            # If next is exactly on the edge (orient=0), it IS the intersection,
            # and it will be added in the next iteration - so don't add here
            if next_orient != 0
                intr_pt = _sh_spherical_intersection(current, next_pt, edge_start, edge_end, T)
                push!(output, intr_pt)
            end
        end

        current_orient = next_orient
    end

    return output
end

# Spherical Polygon-Polygon intersection using Sutherland-Hodgman
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman{Spherical{F}},
    ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b
) where {F, T}
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
        if !isempty(output) && point == output[1]
            continue
        end
        push!(output, UnitSpherical.UnitSphericalPoint{T}(point))
    end

    # Save original subject for containment check
    original_subject = copy(output)

    # Clip against each edge of poly_b
    n_clip = length(clip_points)
    for i in 1:n_clip
        isempty(output) && break
        edge_start = clip_points[i]
        edge_end = clip_points[mod1(i + 1, n_clip)]
        # Skip degenerate edges (duplicate vertices)
        edge_start == edge_end && continue
        output = _sh_clip_to_edge_spherical(output, edge_start, edge_end, T)
    end

    # Handle empty result - check if clip polygon is fully inside the original subject
    if isempty(output)
        if !isempty(clip_points) && _point_in_convex_spherical_polygon(clip_points[1], original_subject)
            # Subject contains clip - return clip polygon
            result = copy(clip_points)
            push!(result, result[1])
            return GI.Polygon([result])
        end
        # Truly disjoint - return degenerate polygon with zero area
        north_pole = UnitSpherical.UnitSphericalPoint{T}(0, 0, 1)
        return GI.Polygon([[north_pole, north_pole, north_pole]])
    end

    # Handle degenerate result (1-2 points can't form a valid ring)
    # This happens for adjacent polygons that share only an edge or corner
    if length(output) < 3
        north_pole = UnitSpherical.UnitSphericalPoint{T}(0, 0, 1)
        return GI.Polygon([[north_pole, north_pole, north_pole]])
    end

    # Close the ring and return
    push!(output, output[1])
    return GI.Polygon([output])
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

# =============================================================================
# Spherical Implementation
# =============================================================================

using ..UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic, GeographicFromUnitSphere
using LinearAlgebra: cross, dot, normalize, norm

"""
    spherical_orient(a, b, c)

Compute the orientation of point `c` with respect to the great circle through `a` and `b`.

Returns:
- Positive if `c` is to the left of the directed great circle arc from `a` to `b`
- Negative if `c` is to the right
- Zero if `c` is on the great circle

The orientation is computed as the triple scalar product: (a × b) · c
"""
function spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    # Triple scalar product: (a × b) · c
    return dot(cross(a, b), c)
end

"""
    _point_in_convex_spherical_polygon(point, polygon_points)

Check if a point is STRICTLY inside a convex spherical polygon.

Returns `true` only if the point is strictly inside (all orientations positive).
Returns `false` if the point is on the boundary (any orientation is zero) or outside.

This is critical for the fallback containment check in Sutherland-Hodgman:
points exactly ON the boundary should NOT be considered "inside", as this would
cause adjacent (non-overlapping) polygons to incorrectly return the clip polygon.
"""
function _point_in_convex_spherical_polygon(
    point::UnitSphericalPoint,
    polygon_points::Vector{<:UnitSphericalPoint}
)
    n = length(polygon_points)
    n < 3 && return false

    for i in 1:n
        edge_start = polygon_points[i]
        edge_end = polygon_points[mod1(i + 1, n)]
        orient = spherical_orient(edge_start, edge_end, point)
        if orient <= 0  # strictly outside OR on the boundary
            return false
        end
    end
    return true  # strictly inside (all orient > 0)
end

# Helper to convert a GeoInterface point to UnitSphericalPoint
function _to_unit_spherical(point)
    return UnitSphericalPoint(GI.PointTrait(), point)
end

# Spherical Polygon-Polygon intersection using Sutherland-Hodgman
function _intersection_sutherland_hodgman(
    alg::ConvexConvexSutherlandHodgman{<:Spherical},
    ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b
) where {T}
    # Get exterior rings (convex polygons have no holes)
    ring_a = GI.getexterior(poly_a)
    ring_b = GI.getexterior(poly_b)

    # Convert poly_a vertices to UnitSphericalPoints (excluding closing point)
    output = UnitSphericalPoint{T}[]
    for point in GI.getpoint(ring_a)
        pt = _to_unit_spherical(point)
        pt_typed = UnitSphericalPoint{T}(pt.x, pt.y, pt.z)
        # Skip the closing point (same as first)
        if !isempty(output) && isapprox(pt_typed, output[1]; atol=eps(T)*10)
            continue
        end
        push!(output, pt_typed)
    end

    # Store original subject polygon for containment check
    original_subject = copy(output)

    # Convert clip polygon vertices (excluding closing point)
    clip_points = UnitSphericalPoint{T}[]
    for point in GI.getpoint(ring_b)
        pt = _to_unit_spherical(point)
        pt_typed = UnitSphericalPoint{T}(pt.x, pt.y, pt.z)
        if !isempty(clip_points) && isapprox(pt_typed, clip_points[1]; atol=eps(T)*10)
            continue
        end
        push!(clip_points, pt_typed)
    end

    # Clip against each edge of poly_b
    n_clip = length(clip_points)
    for i in 1:n_clip
        isempty(output) && break
        edge_start = clip_points[i]
        edge_end = clip_points[mod1(i + 1, n_clip)]
        output = _sh_spherical_clip_to_edge(output, edge_start, edge_end, T)
    end

    # Handle empty result - check if clip polygon is fully inside the original subject
    if isempty(output)
        # Use strict interior check: point ON boundary returns false
        if !isempty(clip_points) && _point_in_convex_spherical_polygon(clip_points[1], original_subject)
            # Subject strictly contains clip - return clip polygon
            result = copy(clip_points)
            push!(result, result[1])
            return _spherical_polygon_to_geo(result, T)
        end
        # No intersection - return degenerate polygon with zero area
        zero_pt = (zero(T), zero(T))
        return GI.Polygon([[zero_pt, zero_pt, zero_pt]])
    end

    # Check for degenerate result (collinear points = line segment, not polygon)
    # This happens when adjacent polygons share an edge
    if _is_degenerate_spherical_polygon(output, T)
        zero_pt = (zero(T), zero(T))
        return GI.Polygon([[zero_pt, zero_pt, zero_pt]])
    end

    # Close the ring and convert back to geographic coordinates
    push!(output, output[1])
    return _spherical_polygon_to_geo(output, T)
end

"""
    _is_degenerate_spherical_polygon(points, T)

Check if a set of spherical points forms a degenerate polygon (zero area).

A polygon is degenerate if:
1. It has fewer than 3 distinct points
2. All points are collinear (lie on a single great circle)

This is critical for detecting edge-sharing cases where Sutherland-Hodgman
produces a "polygon" that is actually just a line segment.
"""
function _is_degenerate_spherical_polygon(
    points::Vector{UnitSphericalPoint{T}},
    ::Type{T}
) where T
    n = length(points)
    n < 3 && return true

    # Remove duplicate/very close points
    unique_points = UnitSphericalPoint{T}[]
    for p in points
        is_dup = false
        for existing in unique_points
            if isapprox(p, existing; atol=eps(T)*100)
                is_dup = true
                break
            end
        end
        if !is_dup
            push!(unique_points, p)
        end
    end

    length(unique_points) < 3 && return true

    # Check if all points are collinear on the sphere
    # Points are collinear if they all lie on the same great circle
    # This is true if for any three points A, B, C: (A × B) · C ≈ 0
    p1 = unique_points[1]
    p2 = unique_points[2]

    # Normal to the great circle through p1 and p2
    normal = cross(p1, p2)
    normal_len_sq = dot(normal, normal)

    # If p1 and p2 are nearly identical or antipodal, can't determine great circle
    if normal_len_sq < eps(T)
        return true
    end

    # Check if all other points lie on the same great circle
    for i in 3:length(unique_points)
        p = unique_points[i]
        # Distance from point to great circle plane
        dist = abs(dot(normal, p)) / sqrt(normal_len_sq)
        if dist > sqrt(eps(T))  # Point is not on the great circle
            return false  # Non-degenerate polygon
        end
    end

    return true  # All points are collinear
end

# Convert UnitSphericalPoints back to geographic polygon
function _spherical_polygon_to_geo(points::Vector{UnitSphericalPoint{T}}, ::Type{T}) where T
    geo_transform = GeographicFromUnitSphere()
    geo_points = [geo_transform(p) for p in points]
    return GI.Polygon([geo_points])
end

# Clip polygon against a single spherical edge using Sutherland-Hodgman rules
function _sh_spherical_clip_to_edge(
    polygon_points::Vector{UnitSphericalPoint{T}},
    edge_start::UnitSphericalPoint{T},
    edge_end::UnitSphericalPoint{T},
    ::Type{T}
) where T
    output = UnitSphericalPoint{T}[]
    n = length(polygon_points)
    n == 0 && return output

    for i in 1:n
        current = polygon_points[i]
        next_pt = polygon_points[mod1(i + 1, n)]

        # Determine if points are inside (left of or on the edge)
        # orient > 0 means left (inside for CCW polygon), == 0 means on edge, < 0 means right (outside)
        current_orient = spherical_orient(edge_start, edge_end, current)
        next_orient = spherical_orient(edge_start, edge_end, next_pt)

        current_inside = current_orient >= 0
        next_inside = next_orient >= 0

        if current_inside
            push!(output, current)
            if !next_inside
                # Exiting: add intersection point
                intr_pt = _sh_spherical_line_intersection(current, next_pt, edge_start, edge_end, T)
                if !isnothing(intr_pt)
                    push!(output, intr_pt)
                end
            end
        elseif next_inside
            # Entering: add intersection point
            intr_pt = _sh_spherical_line_intersection(current, next_pt, edge_start, edge_end, T)
            if !isnothing(intr_pt)
                push!(output, intr_pt)
            end
        end
        # Both outside: add nothing
    end

    return output
end

# Compute intersection point of two great circle arcs
function _sh_spherical_line_intersection(
    p1::UnitSphericalPoint{T},
    p2::UnitSphericalPoint{T},
    p3::UnitSphericalPoint{T},
    p4::UnitSphericalPoint{T},
    ::Type{T}
) where T
    # The intersection of two great circles is found by:
    # 1. Compute normal vectors to each great circle plane
    # 2. The intersection is the cross product of these normals (normalized)

    # Normal to great circle through p1, p2
    n1 = cross(p1, p2)
    # Normal to great circle through p3, p4
    n2 = cross(p3, p4)

    # Cross product gives intersection direction
    intersection_dir = cross(n1, n2)

    # Check if great circles are parallel (no intersection or same circle)
    len_sq = dot(intersection_dir, intersection_dir)
    if len_sq < eps(T)^2
        # Great circles are parallel - return midpoint as fallback
        return nothing
    end

    # Normalize to get point on unit sphere
    intersection_pt = normalize(intersection_dir)

    # The cross product gives us one of two antipodal points
    # We need to pick the one that lies on the arcs (or closer to them)
    # Check if this point is on the correct side of both arcs
    # Use dot product to check if point is in the "positive" direction

    # Check which of the two antipodal points is closer to the arc midpoints
    mid1 = normalize(p1 + p2)
    mid2 = normalize(p3 + p4)

    if dot(intersection_pt, mid1) < 0 || dot(intersection_pt, mid2) < 0
        intersection_pt = -intersection_pt
    end

    return UnitSphericalPoint{T}(intersection_pt.x, intersection_pt.y, intersection_pt.z)
end
