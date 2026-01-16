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
