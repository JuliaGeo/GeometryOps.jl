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
    union_multipoly = tuples(multipoly)
    n_polys = GI.npolygon(multipoly)
    if n_polys > 1
        keep_idx = trues(n_polys)  # keep track of sub-polygons to remove
        # Combine any sub-polygons that intersect
        for (curr_idx, _) in Iterators.filter(last, Iterators.enumerate(keep_idx))
            curr_poly = union_multipoly.geom[curr_idx]
            poly_disjoint = false
            while !poly_disjoint
                poly_disjoint = true  # assume current polygon is disjoint from others
                for (next_idx, _) in Iterators.filter(last, Iterators.drop(Iterators.enumerate(keep_idx), curr_idx))
                    next_poly = union_multipoly.geom[next_idx]
                    if intersects(curr_poly, next_poly)  # if two polygons intersect
                        new_polys = union(curr_poly, next_poly; target = GI.PolygonTrait())
                        n_new_polys = length(new_polys)
                        if n_new_polys == 1  # if polygons combined
                            poly_disjoint = false
                            union_multipoly.geom[curr_idx] = new_polys[1]
                            curr_poly = union_multipoly.geom[curr_idx]
                            keep_idx[next_idx] = false
                        end
                    end
                end
            end
        end
        keepat!(union_multipoly.geom, keep_idx)
    end
    return union_multipoly
end

struct DiffIntersectingPolygons <: GeometryCorrection end

application_level(::DiffIntersectingPolygons) = GI.MultiPolygonTrait

function (::DiffIntersectingPolygons)(::GI.MultiPolygonTrait, multipoly)
    diff_multipoly = tuples(multipoly)
    n_starting_polys = GI.npolygon(multipoly)
    n_polys = n_starting_polys
    if n_polys > 1
        keep_idx = trues(n_polys)  # keep track of sub-polygons to remove
        # Break apart any sub-polygons that intersect
        for curr_idx in 1:n_starting_polys
            !keep_idx[curr_idx] && continue
            for next_idx in (curr_idx + 1):n_starting_polys
                !keep_idx[next_idx] && continue
                next_poly = diff_multipoly.geom[next_idx]
                n_new_polys = 0
                curr_pieces_added = (n_polys + 1):(n_polys + n_new_polys)
                for curr_piece_idx in Iterators.flatten((curr_idx:curr_idx, curr_pieces_added))
                    !keep_idx[curr_piece_idx] && continue
                    curr_poly = diff_multipoly.geom[curr_piece_idx]
                    if intersects(curr_poly, next_poly)  # if two polygons intersect
                        new_polys = difference(curr_poly, next_poly; target = GI.PolygonTrait())
                        n_new_pieces = length(new_polys) - 1
                        if n_new_pieces < 0  # current polygon is covered by next_polygon
                            keep_idx[curr_piece_idx] = false
                            break
                        elseif n_new_pieces ≥ 0
                            diff_multipoly.geom[curr_piece_idx] = new_polys[1]
                            curr_poly = diff_multipoly.geom[curr_piece_idx]
                            if n_new_pieces > 0 # current polygon breaks into several pieces
                                append!(diff_multipoly.geom, @view new_polys[2:end])
                                append!(keep_idx, trues(n_new_pieces))
                                n_new_polys += n_new_pieces
                            end
                        end
                    end
                end
                n_polys += n_new_polys
            end
        end
        keepat!(diff_multipoly.geom, keep_idx)
    end
    return diff_multipoly
end