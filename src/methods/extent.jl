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
area is interpreted to be a CCW loop enclosing a very large area."

The enclosure test is crossing parity, the way `S2Loop::InitBound` decides
pole containment (`s2loop.cc`).  For an anchor edge whose great circle
misses the query point `q`, the side of that edge `q` falls on is the side
the arc from the edge's midpoint to `q` departs into — left is the
interior — and each transversal boundary crossing along the arc flips it.
The departure side equals `q`'s side because the arc can meet the anchor's
great circle again only at the midpoint's antipode, which an arc shorter
than a half turn never reaches.

Where S2 resolves degenerate configurations with exact predicates and
symbolic perturbation, this test detects them — a vertex within
[`UnitSpherical.spherical_orient`](@ref)'s tolerance of a test arc's great
circle, a crossing too close to an arc endpoint to call — and retries with
the next edge as anchor.  If every anchor is degenerate the axis is
extended to `±1`, so the box can come out loose but never under-covers.
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

    lo = MVector(ext.X[1], ext.Y[1], ext.Z[1])
    hi = MVector(ext.X[2], ext.Y[2], ext.Z[2])
    for i in 1:3, s in (1.0, -1.0)
        q = UnitSpherical.UnitSphericalPoint(ntuple(j -> j == i ? s : 0.0, 3))
        inside = _spherical_ring_contains(pts, n, q)
        # nothing (undecidable) extends too; the box must never under-cover
        if inside === nothing || inside
            s > 0 ? (hi[i] = one(hi[i])) : (lo[i] = -one(lo[i]))
        end
    end
    return Extents.Extent(X = (lo[1], hi[1]), Y = (lo[2], hi[2]), Z = (lo[3], hi[3]))
end

# Crossing-parity containment of `q` in the closed region left of the ring,
# after S2Loop::Contains/InitBound.  Returns `nothing` when every anchor edge
# is degenerate with respect to `q`.
function _spherical_ring_contains(pts, n, q)
    for j in 1:n
        UnitSpherical.point_on_spherical_arc(q, pts[j], pts[mod1(j + 1, n)]) && return true
    end
    for j in 1:n
        a, b = pts[j], pts[mod1(j + 1, n)]
        a == b && continue
        side = UnitSpherical.spherical_orient(a, b, q)
        side == 0 && continue
        mid = a + b
        norm(mid) < 1e-9 && continue        # near-antipodal edge, midpoint unstable
        m = UnitSpherical.UnitSphericalPoint(normalize(mid))
        dot(q, m) < -1 + 1e-9 && continue   # test arc q → m would span a half turn
        crossings = 0
        ok = true
        for k in 1:n
            k == j && continue
            c = _arc_crossing_parity(q, m, pts[k], pts[mod1(k + 1, n)])
            if c == -1
                ok = false
                break
            end
            crossings += c
        end
        ok || continue
        # walking from `m` toward `q` departs onto `q`'s side of the anchor
        # edge (the arc meets that great circle again only at `-m`); positive
        # side is the interior, and each crossing flips it
        return isodd(crossings) ? side < 0 : side > 0
    end
    return nothing
end

# Crossing parity of the test arc q → m against ring edge a → b: 1 for a
# transversal crossing, 0 for none, -1 for too close to degenerate to call.
function _arc_crossing_parity(q, m, a, b)
    # a vertex exactly antipodal to `q` lies on every great circle through
    # `q`; its edges can reach the test arc only at `q` itself, excluded by
    # the on-boundary check
    (a == -q || b == -q) && return 0
    sa = UnitSpherical.spherical_orient(q, m, a)
    sb = UnitSpherical.spherical_orient(q, m, b)
    (sa == 0 || sb == 0) && return -1
    sa == sb && return 0
    # `q` on this edge's great circle but off the edge (checked upfront):
    # the circles meet only at `±q`, both out of the test arc's reach — no
    # crossing.  This degeneracy is anchor-independent (a lonlat grid's
    # meridian edges hold `±eₓ`/`±e_y` exactly), so it resolves instead of
    # returning -1.
    sq = UnitSpherical.spherical_orient(a, b, q)
    sq == 0 && return 0
    sm = UnitSpherical.spherical_orient(a, b, m)
    sm == 0 && return -1
    sq == sm && return 0
    # each arc now crosses the other's great circle exactly once, at one of
    # the two antipodal circle intersections; the arcs cross iff those are
    # the same point, i.e. iff the intersection direction `x` points into
    # both arcs' hemispheres
    x = cross(normalize(UnitSpherical.robust_cross_product(q, m)),
              normalize(UnitSpherical.robust_cross_product(a, b)))
    d1 = dot(x, q + m)
    d2 = dot(x, a + b)
    tol = 16 * eps(Float64) * norm(x)
    (abs(d1) <= tol || abs(d2) <= tol) && return -1
    return (d1 > 0) == (d2 > 0) ? 1 : 0
end
