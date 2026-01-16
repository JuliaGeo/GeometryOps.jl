# # Great Circle Arc Intersection
#=
This file implements the intersection of two great circle arcs on the unit sphere.

The fundamental problem is to find where two arcs (AB and CD) intersect. The great circles
containing these arcs have normals N1 = A × B and N2 = C × D. These great circles intersect
at two antipodal points ±(N1 × N2), and we need to determine if either of these points lies
on both arcs.
=#

"""
    ArcIntersectionType

Enumeration of the types of arc intersections:
- `arc_cross`: The arcs cross at a single point in their interiors
- `arc_hinge`: The arcs share exactly one endpoint
- `arc_overlap`: The arcs are collinear and overlap in a segment
- `arc_disjoint`: The arcs do not intersect
"""
@enum ArcIntersectionType arc_cross=1 arc_hinge=2 arc_overlap=3 arc_disjoint=4

"""
    ArcIntersectionResult{T}

Result of computing the intersection between two great circle arcs.

Fields:
- `type::ArcIntersectionType`: The type of intersection
- `points::Vector{UnitSphericalPoint{T}}`: The intersection point(s)
- `fracs::Vector{Tuple{T,T}}`: For each intersection point, the fractional positions (α, β)
  where α is the position along the first arc and β is the position along the second arc.
  Both are in [0, 1], where 0 is the start point and 1 is the end point.
"""
struct ArcIntersectionResult{T}
    type::ArcIntersectionType
    points::Vector{UnitSphericalPoint{T}}
    fracs::Vector{Tuple{T,T}}
end

"""
    spherical_arc_intersection(a1, b1, a2, b2) -> ArcIntersectionResult

Compute the intersection between two great circle arcs on the unit sphere.

The first arc goes from `a1` to `b1`, and the second arc goes from `a2` to `b2`.
All points should be `UnitSphericalPoint` instances.

Returns an `ArcIntersectionResult` containing:
- The type of intersection (cross, hinge, overlap, or disjoint)
- The intersection point(s) if they exist
- The fractional positions along each arc for each intersection point

# Examples

```julia
# Two arcs that cross
a1 = UnitSphereFromGeographic()((-45.0, 0.0))
b1 = UnitSphereFromGeographic()((45.0, 0.0))
a2 = UnitSphereFromGeographic()((0.0, -45.0))
b2 = UnitSphereFromGeographic()((0.0, 45.0))
result = spherical_arc_intersection(a1, b1, a2, b2)
# result.type == arc_cross, with one intersection point at (0°, 0°)
```
"""
function spherical_arc_intersection(
    a1::UnitSphericalPoint{T1},
    b1::UnitSphericalPoint{T2},
    a2::UnitSphericalPoint{T3},
    b2::UnitSphericalPoint{T4}
) where {T1, T2, T3, T4}
    # Promote to common type
    T = promote_type(T1, T2, T3, T4)
    a1 = UnitSphericalPoint{T}(a1)
    b1 = UnitSphericalPoint{T}(b1)
    a2 = UnitSphericalPoint{T}(a2)
    b2 = UnitSphericalPoint{T}(b2)

    tol = eps(T) * 16

    # Check for endpoint equality (hinge cases)
    # There are 4 possible hinge configurations
    if _points_equal(a1, a2, tol)
        α = zero(T)
        β = zero(T)
        return _make_hinge_result(T, a1, α, β)
    elseif _points_equal(a1, b2, tol)
        α = zero(T)
        β = one(T)
        return _make_hinge_result(T, a1, α, β)
    elseif _points_equal(b1, a2, tol)
        α = one(T)
        β = zero(T)
        return _make_hinge_result(T, b1, α, β)
    elseif _points_equal(b1, b2, tol)
        α = one(T)
        β = one(T)
        return _make_hinge_result(T, b1, α, β)
    end

    # Compute great circle normals using robust cross product
    n1 = robust_cross_product(a1, b1)
    n2 = robust_cross_product(a2, b2)

    # Normalize the normals
    n1_norm = norm(n1)
    n2_norm = norm(n2)

    # Check if either arc is degenerate (endpoints are identical or antipodal)
    if n1_norm < tol || n2_norm < tol
        # At least one arc is degenerate
        return ArcIntersectionResult{T}(arc_disjoint, UnitSphericalPoint{T}[], Tuple{T,T}[])
    end

    n1 = n1 / n1_norm
    n2 = n2 / n2_norm

    # Check if the arcs are collinear (same great circle)
    # This happens when n1 and n2 are parallel (same direction or opposite)
    cross_n1_n2 = n1 × n2
    if norm(cross_n1_n2) < tol
        # Arcs are collinear - need special handling
        return _find_collinear_arc_intersection(a1, b1, a2, b2, T)
    end

    # Compute the intersection direction as n1 × n2
    # This gives us two antipodal intersection points
    intersection_dir = cross_n1_n2 / norm(cross_n1_n2)

    # The two candidate intersection points are ±intersection_dir
    candidates = [
        UnitSphericalPoint{T}(intersection_dir),
        UnitSphericalPoint{T}(-intersection_dir)
    ]

    # Check which candidate(s) lie on both arcs
    valid_intersections = UnitSphericalPoint{T}[]
    valid_fracs = Tuple{T,T}[]

    for candidate in candidates
        on_arc1 = point_on_spherical_arc(candidate, a1, b1)
        on_arc2 = point_on_spherical_arc(candidate, a2, b2)

        if on_arc1 && on_arc2
            # This candidate is a valid intersection
            α = _arc_fraction(candidate, a1, b1)
            β = _arc_fraction(candidate, a2, b2)
            push!(valid_intersections, candidate)
            push!(valid_fracs, (α, β))
        end
    end

    # Determine the type of intersection
    if isempty(valid_intersections)
        return ArcIntersectionResult{T}(arc_disjoint, valid_intersections, valid_fracs)
    elseif length(valid_intersections) == 1
        return ArcIntersectionResult{T}(arc_cross, valid_intersections, valid_fracs)
    else
        # Two intersection points - this should be very rare (essentially antipodal arcs)
        # For now, treat as cross and return both points
        return ArcIntersectionResult{T}(arc_cross, valid_intersections, valid_fracs)
    end
end

# Helper functions

"""
    _points_equal(a::UnitSphericalPoint, b::UnitSphericalPoint, tol) -> Bool

Check if two points are equal within a tolerance.
"""
function _points_equal(a::UnitSphericalPoint, b::UnitSphericalPoint, tol)
    return norm(a - b) < tol
end

"""
    _make_hinge_result(T, point, α, β) -> ArcIntersectionResult{T}

Create an `ArcIntersectionResult` for a hinge intersection.
"""
function _make_hinge_result(T, point::UnitSphericalPoint, α, β)
    return ArcIntersectionResult{T}(
        arc_hinge,
        [UnitSphericalPoint{T}(point)],
        [(α, β)]
    )
end

"""
    _arc_fraction(p, a, b) -> T

Compute the fractional position of point `p` along the arc from `a` to `b`.

Returns a value in [0, 1] where 0 corresponds to `a` and 1 corresponds to `b`.
Uses spherical distance to compute the fraction.
"""
function _arc_fraction(p::UnitSphericalPoint{T}, a::UnitSphericalPoint, b::UnitSphericalPoint) where T
    # Compute distances
    dist_ab = spherical_distance(a, b)
    dist_ap = spherical_distance(a, p)

    # Handle edge cases
    if dist_ab < eps(T) * 16
        # Arc is degenerate
        return zero(T)
    end

    # The fraction is dist_ap / dist_ab
    frac = dist_ap / dist_ab

    # Clamp to [0, 1] to handle numerical errors
    return clamp(frac, zero(T), one(T))
end

"""
    _find_collinear_arc_intersection(a1, b1, a2, b2, T) -> ArcIntersectionResult{T}

Find the intersection of two collinear arcs (arcs on the same great circle).

Returns an `ArcIntersectionResult` with type `arc_overlap` if the arcs overlap,
or `arc_disjoint` if they don't.
"""
function _find_collinear_arc_intersection(
    a1::UnitSphericalPoint{T},
    b1::UnitSphericalPoint{T},
    a2::UnitSphericalPoint{T},
    b2::UnitSphericalPoint{T},
    ::Type{T}
) where T
    # For collinear arcs, we need to check if they overlap
    # We do this by checking if any endpoint of one arc lies on the other arc
    # or vice versa

    tol = eps(T) * 16

    # Check all possible overlaps
    a2_on_arc1 = point_on_spherical_arc(a2, a1, b1)
    b2_on_arc1 = point_on_spherical_arc(b2, a1, b1)
    a1_on_arc2 = point_on_spherical_arc(a1, a2, b2)
    b1_on_arc2 = point_on_spherical_arc(b1, a2, b2)

    # Collect all endpoints that define the overlap
    overlap_points = UnitSphericalPoint{T}[]
    overlap_fracs = Tuple{T,T}[]

    # If arcs overlap, the overlap is defined by the "inner" endpoints
    # We need to identify which endpoints bound the overlapping region

    if a2_on_arc1 && !_points_equal(a2, a1, tol) && !_points_equal(a2, b1, tol)
        α = _arc_fraction(a2, a1, b1)
        β = zero(T)
        push!(overlap_points, a2)
        push!(overlap_fracs, (α, β))
    end

    if b2_on_arc1 && !_points_equal(b2, a1, tol) && !_points_equal(b2, b1, tol)
        α = _arc_fraction(b2, a1, b1)
        β = one(T)
        push!(overlap_points, b2)
        push!(overlap_fracs, (α, β))
    end

    if a1_on_arc2 && !_points_equal(a1, a2, tol) && !_points_equal(a1, b2, tol)
        α = zero(T)
        β = _arc_fraction(a1, a2, b2)
        push!(overlap_points, a1)
        push!(overlap_fracs, (α, β))
    end

    if b1_on_arc2 && !_points_equal(b1, a2, tol) && !_points_equal(b1, b2, tol)
        α = one(T)
        β = _arc_fraction(b1, a2, b2)
        push!(overlap_points, b1)
        push!(overlap_fracs, (α, β))
    end

    # If we found overlap points, return them
    if !isempty(overlap_points)
        # Sort by α (fraction along first arc) to get consistent ordering
        perm = sortperm([f[1] for f in overlap_fracs])
        return ArcIntersectionResult{T}(
            arc_overlap,
            overlap_points[perm],
            overlap_fracs[perm]
        )
    end

    # Check if one arc completely contains the other
    if a2_on_arc1 && b2_on_arc1
        # Arc2 is completely inside Arc1
        α1 = _arc_fraction(a2, a1, b1)
        α2 = _arc_fraction(b2, a1, b1)
        return ArcIntersectionResult{T}(
            arc_overlap,
            [a2, b2],
            [(α1, zero(T)), (α2, one(T))]
        )
    elseif a1_on_arc2 && b1_on_arc2
        # Arc1 is completely inside Arc2
        β1 = _arc_fraction(a1, a2, b2)
        β2 = _arc_fraction(b1, a2, b2)
        return ArcIntersectionResult{T}(
            arc_overlap,
            [a1, b1],
            [(zero(T), β1), (one(T), β2)]
        )
    end

    # No overlap
    return ArcIntersectionResult{T}(arc_disjoint, UnitSphericalPoint{T}[], Tuple{T,T}[])
end
