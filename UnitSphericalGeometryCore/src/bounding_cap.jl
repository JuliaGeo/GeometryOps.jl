# # Caps and points from geometries

#=
```@docs; canonical=false
spherical_bounding_cap
to_unit_spherical_points
```

These are the bridges from an arbitrary GeoInterface geometry to the unit sphere:
[`to_unit_spherical_points`](@ref) converts the vertices, and
[`spherical_bounding_cap`](@ref) computes a [`SphericalCap`](@ref) that contains the
geometry.  The cap is cheap (two passes, no allocation) rather than minimal; for the
smallest enclosing cap use a minimum bounding circle algorithm.
=#

"""
    to_unit_spherical_points(geom) -> Vector{<:UnitSphericalPoint}

Convert the vertices of `geom` (any GeoInterface geometry with points, such as a ring,
linestring or polygon) to `UnitSphericalPoint`s.  Two-dimensional input is interpreted
geographically, as (longitude, latitude) in degrees; three-dimensional input is taken to
already be unit-sphere `xyz`, and existing `UnitSphericalPoint`s pass through unchanged.
"""
function to_unit_spherical_points(geom)
    return [_unit_spherical_point(p) for p in GI.getpoint(geom)]
end

# Convert one GeoInterface point: lon/lat to unit xyz, or take 3D input as unit xyz.
# Unlike the `UnitSphericalPoint` constructor, this respects `GI.is3d`.
_unit_spherical_point(p::UnitSphericalPoint) = p
function _unit_spherical_point(p)
    if GI.is3d(p)
        return UnitSphericalPoint(GI.x(p), GI.y(p), GI.z(p))
    else
        return UnitSphereFromGeographic()(p)
    end
end

"""
    spherical_bounding_cap(geom) -> SphericalCap

A [`SphericalCap`](@ref) that contains every point of `geom`, including the great-circle
arcs between consecutive vertices.  Points are interpreted as in
[`to_unit_spherical_points`](@ref).

The cap is centered on the normalized sum of the vertices and reaches the furthest vertex,
padded slightly so that comparisons against it stay conservative.  A cap smaller than a
quarter turn is geodesically convex, so it contains the minor arcs between its vertices;
when the vertices spread further than that, or their sum is degenerate (spread over a great
circle), the whole sphere (a cap of radius `π`) is returned instead.

This is a bounding cap, not the minimum one: it is designed to be cheap (two passes over
the vertices, no allocation) for use as a rejection test.

## Example

```jldoctest
julia> using UnitSphericalGeometryCore; import GeoInterface as GI

julia> cap = spherical_bounding_cap(GI.LineString([(0.0, 0.0), (10.0, 0.0)]));

julia> cap.point ≈ UnitSphereFromGeographic()((5.0, 0.0))
true

julia> isapprox(cap.radius, deg2rad(5); atol = 1e-6)
true
```
"""
spherical_bounding_cap(geom) = spherical_bounding_cap(GI.trait(geom), geom)

# `GI.getpoint` has no `PointTrait` method; a point is its own (padded) zero-radius cap.
spherical_bounding_cap(::GI.PointTrait, geom) =
    SphericalCap(_unit_spherical_point(geom), _BOUNDING_CAP_PAD)

function spherical_bounding_cap(::GI.AbstractGeometryTrait, geom)
    sx = sy = sz = 0.0
    n = 0
    for p in GI.getpoint(geom)
        u = _unit_spherical_point(p)
        sx += u[1]; sy += u[2]; sz += u[3]
        n += 1
    end
    n == 0 && return _whole_sphere_cap()
    nrm = sqrt(sx * sx + sy * sy + sz * sz)
    nrm < 1e-9 && return _whole_sphere_cap()  # vertices spread over a great circle
    c = UnitSphericalPoint(sx / nrm, sy / nrm, sz / nrm)
    maxchord2 = 0.0
    for p in GI.getpoint(geom)
        u = _unit_spherical_point(p)
        dx = c[1] - u[1]; dy = c[2] - u[2]; dz = c[3] - u[3]
        ch2 = dx * dx + dy * dy + dz * dz
        ch2 > maxchord2 && (maxchord2 = ch2)
    end
    maxchord2 >= 2.0 && return _whole_sphere_cap()  # radius past a quarter turn
    # Compute the radius from chord length for small-angle precision, then pad it.
    r = 2 * asin(sqrt(maxchord2) / 2)
    return SphericalCap(c, r + _BOUNDING_CAP_PAD + 8 * eps(r))
end

const _BOUNDING_CAP_PAD = 1e-7
_whole_sphere_cap() = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), Float64(π))

# A cap built directly from a non-point geometry is its bounding cap.
SphericalCap(::GI.AbstractGeometryTrait, geom) = spherical_bounding_cap(geom)
