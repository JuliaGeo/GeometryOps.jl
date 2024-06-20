# # Closed Rings

export ClosedRing

# A closed ring is a ring that has the same start and end point. This is a requirement for a valid polygon (technically, for a valid LinearRing). 
# This correction is used to ensure that the polygon is valid.

# The reason this operates on the polygon level is that several packages are loose about whether they return LinearRings (which is correct) or LineStrings (which is incorrect) for the contents of a polygon.
# Therefore, we decompose manually to ensure correctness.

# ## Example
#=

Many polygon providers do not close their polygons, which makes them invalid 
according to the specification.  Quite a few geometry algorithms assume that 
polygons are closed, and leaving them open can lead to incorrect results!

For example, the following polygon is not valid:

```@example closed-ring
import GeoInterface as GI
polygon = GI.Polygon([[(0, 0), (1, 0), (1, 1), (0, 1)]])
```

even though it will look correct when visualized, and indeed appears correct.

```@example closed-ring
import GeometryOps as GO
GO.fix(polygon, corrections = [GO.ClosedRing()])
```

You can see that the last point of the ring here is equal to the first point. For a polygon with ``n`` sides, there should be ``n+1`` vertices.

=#

# ## Implementation

"""
    ClosedRing() <: GeometryCorrection

This correction ensures that a polygon's exterior and interior rings are closed.

It can be called on any geometry correction as usual.

See also [`GeometryCorrection`](@ref).
"""
struct ClosedRing <: GeometryCorrection end

application_level(::ClosedRing) = GI.PolygonTrait

function (::ClosedRing)(::Type{T}, ::GI.PolygonTrait, polygon) where T
    exterior = _close_linear_ring(T, GI.getexterior(polygon))
    
    holes = map(GI.gethole(polygon)) do hole
        _close_linear_ring(T, hole) # TODO: make this more efficient, or use tuples!
    end

    return GI.Wrappers.Polygon([exterior, holes...])
end

function _close_linear_ring(::Type{T}, ring) where T
    ring = svpoints(ring, T)
    if !equals(GI.getpoint(ring, 1), GI.getpoint(ring, GI.npoint(ring)))
        # Close the ring
        push!(ring.geom, ring.geom[1])
    end
    return ring
end