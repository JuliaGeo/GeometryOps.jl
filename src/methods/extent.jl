# # Extent

#=
```@docs; canonical=false
Extents.extent(::Manifold, ::Any)
```

`Extents.extent(m::Manifold, geom)` computes the extent of a geometry *on a
manifold*.  On `Planar()` it delegates to `GI.extent`.  On `Spherical()` it
returns 3D Cartesian `Extent{(:X, :Y, :Z)}`s on the unit sphere — the boxes
that spatial indices over spherical geometry prune with.

An extremum of a coordinate over a *region* on the sphere lies either on
the boundary or at an on-sphere critical point of that coordinate — across
the three coordinates, the six axis points `(±1,0,0)`, `(0,±1,0)`,
`(0,0,±1)`.  A region's box is therefore the union of its edges'
[`UnitSpherical.spherical_arc_extent`](@ref)s plus an enclosure check per
axis point, decided by crossing parity as in `S2Loop::InitBound`
(`s2loop.cc`).

Rings follow S2's loop convention (`s2loop.h`): CCW, interior on the left,
so a clockwise ring encloses the complement.  Configurations too close to
degenerate for [`UnitSpherical.spherical_orient`](@ref) to call are retried
with the next edge as anchor; if every anchor fails, the axis is extended
to `±1`, so the box can come out loose but never under-covers.
=#

"""
    extent(m::Manifold, geom, [::Type{T} = Float64])::Extents.Extent

The extent of `geom` on the manifold `m`, as an `Extents.Extent`.  The
method extends `Extents.extent` (GeometryOps does not export `extent`), so
call it as `GO.extent(m, geom)`.

On `Planar()`, `GI.extent(geom)`.  On `Spherical()`, the 3D Cartesian
extent of the geometry on the unit sphere, with geographic (longitude,
latitude) input converted like `UnitSphericalPoint`: curves are covered by
the union of their edges' great-circle arc extents; rings and polygons are
regions — wound CCW with the interior on the left, per S2's loop
convention — whose extent also covers any enclosed axis point (a pole,
say).

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

    lo = MVector(ext.X[1], ext.Y[1], ext.Z[1])
    hi = MVector(ext.X[2], ext.Y[2], ext.Z[2])
    for i in 1:3, s in (1.0, -1.0)
        q = UnitSpherical.UnitSphericalPoint(ntuple(j -> j == i ? s : 0.0, 3))
        inside = UnitSpherical.spherical_ring_contains(pts, n, q)
        # nothing (undecidable) extends too; the box must never under-cover
        if inside === nothing || inside
            s > 0 ? (hi[i] = one(hi[i])) : (lo[i] = -one(lo[i]))
        end
    end
    return Extents.Extent(X = (lo[1], hi[1]), Y = (lo[2], hi[2]), Z = (lo[3], hi[3]))
end
