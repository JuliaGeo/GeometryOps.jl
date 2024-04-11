# # Intersecting Polygons

export UnionIntersectingPolygons

#=
If the sub-polygons of a multipolygon are intersecting, this makes them invalid according to
specification. Each sub-polygon of a multipolygon being disjoint (other than by a single
point) is a requirment for a valid multipolygon. However, different libraries may achieve
this in different ways. 

For example, taking the union of all sub-polygons of a multipolygon will create a new
multipolygon where each sub-polygon is disjoint. This can be done with the
`UnionIntersectingPolygons` correction.

The reason this operates on a multipolygon level is that it is easy for users to mistakenly
create multipolygon's that overlap, which can then be detrimental to polygon clipping
performance and even create wrong answers.
=#

# ## Example
#=

Multipolygon providers may not check that the polygons making up their multipolygons do not
intersect, which makes them invalid according to the specification.

For example, the following multipolygon is not valid:

```@example union-multipoly
import GeoInterface as GI
polygon = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
multipolygon = GI.MultiPolygon([polygon, polygon])
```

given that the two sub-polygons are the exact same shape.

```@example union-multipoly
import GeometryOps as GO
GO.fix(multipolygon, corrections = [GO.UnionIntersectingPolygons()])
```

You can see that the the multipolygon now only contains one sub-polygon, rather than the two
identical ones provided.
=#

# ## Implementation

"""
    UnionIntersectingPolygons() <: GeometryCorrection

This correction ensures that the polygon's included in a multipolygon aren't intersecting.
If any polygon's are intersecting, they will be combined through the union operation to
create a unique set of disjoint (other than potentially connections by a single point)
polygons covering the same area.

See also [`GeometryCorrection`](@ref).
"""
struct UnionIntersectingPolygons <: GeometryCorrection end

application_level(::UnionIntersectingPolygons) = GI.MultiPolygonTrait

function (::UnionIntersectingPolygons)(::GI.MultiPolygonTrait, multipoly)
    union_multipoly = if GI.npolygon(multipoly) > 1
        # Combine any sub-polygons that intersect
        first_poly = GI.getpolygon(multipoly, 1)
        exclude_first_poly = GI.MultiPolygon(collect(Iterators.drop(GI.getpolygon(multipoly), 1)))
        GI.MultiPolygon(union(first_poly, exclude_first_poly; target = GI.PolygonTrait(), fix_multipoly = nothing))
    else
        tuples(multipoly)
    end
    return union_multipoly
end