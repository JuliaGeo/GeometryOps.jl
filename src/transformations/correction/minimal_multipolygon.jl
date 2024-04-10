# # Minimal MultiPolygons

export MinimalMultiPolygon

#=
A minimal multipolygon is a multipolygon where each individual sub-polygon making up the
multipolygon is not intersecting. This is a requirment for a valid multipolygon. However,
different libraries may achieve this in different ways. 

The reason this operates on a multipolygon level is that it is easy for users to mistakenly
create multipolygon's that overlap, which can then be detrimental to polygon clipping
performance and even create wrong answers.
=#

# ## Example
#=

Multipolygon providers may not check that the polygons making up their multipolygons do not
intersect, which makes them invalid according to the specification.

For example, the following multipolygon is not valid:

```@example minimal-multipoly
import GeoInterface as GI
polygon = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
multipolygon = GI.MultiPolygon([polygon, polygon])
```

given that the two sub-polygons are the exact same shape.

```@example minimal-multipoly
import GeometryOps as GO
GO.fix(multipolygon, corrections = [GO.MinimalMultiPolygon()])
```

You can see that the the multipolygon now only contains one sub-polygon, rather than the two
identical ones provided.

=#

# ## Implementation

"""
    MinimalMultiPolygon() <: GeometryCorrection

This correction ensures that the polygon's included in a multipolygon aren't intersecting.
If any polygon's are intersecting, they will be combined to create a unique set of
non-intersecting polygons covering the same area.

It can be called on any geometry correction as usual.

See also [`GeometryCorrection`](@ref).
"""

struct MinimalMultiPolygon <: GeometryCorrection end

application_level(::MinimalMultiPolygon) = GI.MultiPolygonTrait

function (::MinimalMultiPolygon)(::GI.MultiPolygonTrait, multipoly)
    minimal_multipoly = if GI.npolygon(multipoly) > 1
        # Combine any sub-polygons that intersect
        first_poly = GI.getpolygon(multipoly, 1)
        exclude_first_poly = GI.MultiPolygon(collect(Iterators.drop(GI.getpolygon(multipoly), 1)))
        GI.MultiPolygon(union(first_poly, exclude_first_poly; target = GI.PolygonTrait(), fix_multipoly = false))
    else
        multipoly
    end
    return minimal_multipoly
end