# # Spherical Predicates
#=
This file contains geometric predicates for spherical geometry on the unit sphere.
These predicates determine spatial relationships between points and arcs on the sphere.
=#

"""
    spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint) -> Int

Determine the orientation of point `c` with respect to the great circle arc from `a` to `b`.

Returns:
- `1` if `c` is to the left of the arc (counter-clockwise)
- `-1` if `c` is to the right of the arc (clockwise)
- `0` if `c` is on the great circle (collinear)

Uses [`robust_cross_product`](@ref) for numerical stability with nearly identical or antipodal points.

# Examples
```jldoctest
using GeometryOps.UnitSpherical: UnitSphericalPoint, spherical_orient
a = UnitSphericalPoint(1.0, 0.0, 0.0)
b = UnitSphericalPoint(0.0, 1.0, 0.0)
c = UnitSphericalPoint(0.0, 0.0, 1.0)
spherical_orient(a, b, c)
# output
1
```
"""
function spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    # The orientation is determined by sign((a × b) · c)
    # Use robust_cross_product for numerical stability
    n = robust_cross_product(a, b)
    dot_product = n ⋅ c

    # Use a tolerance for near-zero values
    tol = eps(Float64) * 16  # Same tolerance as S2 geometry
    if abs(dot_product) < tol
        return 0
    end
    return dot_product > 0 ? 1 : -1
    # return ExactPredicates.orient(a, b, UnitSphericalPoint((0., 0., 0.)), c)
end

# Convenience method for raw vectors
function spherical_orient(a::AbstractVector, b::AbstractVector, c::AbstractVector)
    return spherical_orient(
        UnitSphericalPoint(a),
        UnitSphericalPoint(b),
        UnitSphericalPoint(c)
    )
end

"""
    point_on_spherical_arc(p::UnitSphericalPoint, a::UnitSphericalPoint, b::UnitSphericalPoint) -> Bool

Check if point `p` lies on the great circle arc from `a` to `b`.

The arc is the shorter path along the great circle connecting `a` and `b`.
Returns `true` if `p` is on the arc (including endpoints), `false` otherwise.

# Examples
```jldoctest
using GeometryOps.UnitSpherical: UnitSphericalPoint, point_on_spherical_arc
a = UnitSphericalPoint(1.0, 0.0, 0.0)
b = UnitSphericalPoint(0.0, 1.0, 0.0)
mid = UnitSphericalPoint(1/√2, 1/√2, 0.0)
point_on_spherical_arc(mid, a, b)
# output
true
```
"""
function point_on_spherical_arc(p::UnitSphericalPoint, a::UnitSphericalPoint, b::UnitSphericalPoint)
    # First check: is p on the great circle through a and b?
    if spherical_orient(a, b, p) != 0
        return false
    end

    # Second check: is p between a and b on the arc?
    # For the shorter arc, p is between a and b if:
    # (a · p) ≥ (a · b) and (b · p) ≥ (a · b)
    # This works because dot product on unit sphere = cos(angle)
    # If p is between a and b, the angles a-p and b-p are both ≤ angle a-b

    ab = a ⋅ b  # cos(angle between a and b)
    ap = a ⋅ p  # cos(angle between a and p)
    bp = b ⋅ p  # cos(angle between b and p)

    tol = eps(Float64) * 16

    # p is on arc if it's "closer" to both endpoints than they are to each other
    # (in terms of angle, so larger dot product)
    return (ap ≥ ab - tol) && (bp ≥ ab - tol)
end
