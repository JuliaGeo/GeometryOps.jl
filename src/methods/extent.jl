# # Extent

#=
```@docs; canonical=false
Extents.extent(::Manifold, ::Any)
```

`Extents.extent(m::Manifold, geom)` computes the extent of a geometry *on a
manifold*.  On `Planar()` it delegates to `GI.extent`.  On `Spherical()` it
returns 3D Cartesian `Extent{(:X, :Y, :Z)}`s on the unit sphere — the boxes
that spatial indices over spherical geometry prune with.

On the sphere, an extremum of a coordinate over a *region* is attained
either on the boundary or at one of the six on-sphere critical points of
that coordinate: `(±1,0,0)`, `(0,±1,0)`, `(0,0,±1)`.  A polygon that
strictly encloses one of them (a cell over a pole, say) therefore extends
past every boundary edge's extent, so a region's box is the union of its
edges' [`UnitSpherical.spherical_arc_extent`](@ref)s plus an enclosure
check per axis point.

Enclosure follows S2's loop convention (`s2loop.h`): "All loops are defined
to have a CCW orientation, i.e. the interior of the loop is on the left
side of the edges.  This implies that a clockwise loop enclosing a small
area is interpreted to be a CCW loop enclosing a very large area."  A ring
encloses exactly one of `±eᵢ` iff its winding number about that axis is
`±1`, and the sign picks which; both are enclosed only when the interior is
the larger side of the ring (negative signed area, or area above `2π`), in
which case the axis is extended to `[-1, 1]`.

Two documented approximations, both conservative-safe for meshes: the
winding accumulates wrapped angle deltas, so a single edge must not sweep
more than a half turn about an axis (an edge passing closer to `±eᵢ` than
roughly its own length) — such an edge's own arc extent already reaches
within `(distance)²/2` of `±1`.  And a region whose interior strictly
contains an antipodal pair away from the axis points (a thin tube pole to
pole) is beyond the winding test; its boundary extents again come within
`(distance)²/2` of the truth.
=#

"""
    extent(m::Manifold, geom, [::Type{T} = Float64])::Extents.Extent

The extent of `geom` on the manifold `m` — this method lives on, and
returns an, `Extents.Extent`.

On `Planar()`, `GI.extent(geom)`.  On `Spherical()`, the 3D Cartesian
extent of the geometry on the unit sphere, with geographic (longitude,
latitude) input converted like `UnitSphericalPoint`: curves are covered by
the union of their edges' great-circle arc extents, and rings and polygons
are treated as regions — wound CCW with the interior on the left, per S2's
loop convention — whose extent also covers any enclosed pole or other
on-sphere axis extreme.

## Example

```jldoctest
julia> import GeometryOps as GO, GeoInterface as GI

julia> cap = GI.Polygon([[(lon, 60.0) for lon in 0.0:30.0:360.0]]);  # around the north pole

julia> GO.extent(GO.Spherical(), cap).Z[2]
1.0
```
"""
function Extents.extent(m::Manifold, geom, ::Type{T} = Float64) where T
    return _extent(m, GI.trait(geom), geom, T)
end

_extent(::Planar, trait, geom, ::Type{T}) where T = GI.extent(geom)

_extent(::Spherical, ::GI.PointTrait, geom, ::Type{T}) where T =
    GI.extent(UnitSpherical.UnitSphericalPoint(geom))
_extent(m::Spherical, ::Union{GI.LineTrait, GI.LineStringTrait}, geom, ::Type{T}) where T =
    mapreduce(GI.extent, Extents.union, lazy_edgelist(m, geom, T))
_extent(m::Spherical, ::GI.LinearRingTrait, geom, ::Type{T}) where T =
    _spherical_region_extent(UnitSpherical.to_unit_spherical_points(geom))
_extent(m::Spherical, ::GI.PolygonTrait, geom, ::Type{T}) where T =
    _extent(m, GI.LinearRingTrait(), GI.getexterior(geom), T)
# multi-geometries and collections; holes never extend a polygon's extent,
# so only the exterior ring above matters
_extent(m::Spherical, ::GI.AbstractGeometryTrait, geom, ::Type{T}) where T =
    mapreduce(g -> Extents.extent(m, g, T), Extents.union, GI.getgeom(geom))

function _spherical_region_extent(pts::Vector{<:UnitSpherical.UnitSphericalPoint})
    n = length(pts)
    n > 1 && pts[end] == pts[1] && (n -= 1)
    ext = mapreduce(Extents.union, 1:n) do i
        UnitSpherical.spherical_arc_extent(pts[i], pts[mod1(i + 1, n)])
    end
    n < 3 && return ext

    # winding about each axis from wrapped angle deltas in the plane
    # perpendicular to it; a vertex on an axis makes that winding
    # meaningless, but also puts ±1 into the edge extents above, so skip it
    winding = zeros(MVector{3, Float64})
    onaxis = MVector(false, false, false)
    angles(p) = (atan(p.z, p.y), atan(p.x, p.z), atan(p.y, p.x))
    prev = angles(pts[n])
    for i in 1:n
        p = pts[i]
        onaxis[1] |= p.y == 0 && p.z == 0
        onaxis[2] |= p.z == 0 && p.x == 0
        onaxis[3] |= p.x == 0 && p.y == 0
        cur = angles(p)
        winding .+= rem.(cur .- prev, 2π, RoundNearest)
        prev = cur
    end

    ringarea = sum(i -> _spherical_triangle_area(Girard(), pts[1], pts[i], pts[i + 1]), 2:(n - 1); init = 0.0)
    bigregion = ringarea < 0 || ringarea > 2π

    axis_bounds = values(ext)
    bounds = ntuple(3) do i
        lo, hi = axis_bounds[i]
        if !onaxis[i]
            winding[i] > π && (hi = one(hi))
            winding[i] < -π && (lo = -one(lo))
            if bigregion && abs(winding[i]) <= π
                lo, hi = -one(lo), one(hi)
            end
        end
        (lo, hi)
    end
    return Extents.Extent(X = bounds[1], Y = bounds[2], Z = bounds[3])
end
