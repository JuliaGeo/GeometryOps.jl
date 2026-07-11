# # Antipodal Edge Split

export AntipodalEdgeSplit

#=
On the sphere an edge between two *antipodal* vertices (a point and its
antipode) has no unique great-circle arc — infinitely many great circles pass
through both — so `relate` with a `Spherical` manifold refuses such an edge
at ingest (the kernel's edge validation). This correction
is the documented remedy: it splits every antipodal edge by inserting the
lon/lat midpoint of its endpoints, replacing one ambiguous edge with two
well-defined (roughly quarter-circle) arcs.

## Example
=#
# ```@example antipodal
# import GeometryOps as GO, GeoInterface as GI
# # the edge (0,0)→(180,0) maps to the antipodal unit vectors (1,0,0) and (-1,0,0)
# polygon = GI.Polygon([GI.LinearRing([(0., 0.), (180., 0.), (90., 80.), (0., 0.)])])
# GO.fix(polygon; corrections = [GO.AntipodalEdgeSplit()])
# ```
#=
The corrected ring carries the inserted midpoint `(90, 0)`, after which
`relate(GO.RelateNG(; manifold = GO.Spherical()), …)` runs without error.

## Implementation
=#

"""
    AntipodalEdgeSplit() <: GeometryCorrection

Split every edge whose endpoints map to exactly-antipodal unit vectors by
inserting the lon/lat midpoint of the edge, so each edge has a well-defined
great-circle arc. This is the remedy for the antipodal-edge `ArgumentError`
thrown by `relate` on the `Spherical` manifold.

It can be called on any geometry as usual (`AntipodalEdgeSplit()(geom)`), or
passed to [`fix`](@ref).

See also [`GeometryCorrection`](@ref).
"""
struct AntipodalEdgeSplit <: GeometryCorrection end

application_level(::AntipodalEdgeSplit) = GI.PolygonTrait

function (::AntipodalEdgeSplit)(::GI.PolygonTrait, polygon)
    exterior = _split_antipodal_curve(GI.getexterior(polygon))
    holes = map(_split_antipodal_curve, GI.gethole(polygon))
    return GI.Wrappers.Polygon([exterior, holes...])
end

(::AntipodalEdgeSplit)(::GI.AbstractCurveTrait, curve) = _split_antipodal_curve(curve)

# Whether the lon/lat points `p`, `q` map to exactly-antipodal unit vectors —
# the same condition the kernel's edge validation throws on (vanishing cross,
# negative dot).
function _is_antipodal_lonlat(p, q)
    up = _spherical_kernel_point(p)
    uq = _spherical_kernel_point(q)
    n = cross(up, uq)
    return iszero(n[1]) && iszero(n[2]) && iszero(n[3]) && (up ⋅ uq) < 0
end

# Insert the lon/lat midpoint into every antipodal edge of a ring/line,
# returning the input unchanged (no copy) when there is nothing to split.
function _split_antipodal_curve(curve)
    pts = [tuples(p) for p in GI.getpoint(curve)]
    split = false
    out = similar(pts, 0)
    push!(out, pts[1])
    for i in 1:(length(pts) - 1)
        p, q = pts[i], pts[i + 1]
        if _is_antipodal_lonlat(p, q)
            push!(out, ((p[1] + q[1]) / 2, (p[2] + q[2]) / 2))
            split = true
        end
        push!(out, q)
    end
    split || return curve
    return GI.geointerface_geomtype(GI.trait(curve))(out)
end
