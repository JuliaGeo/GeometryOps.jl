# # Overlaps

export overlaps

#=
## What is overlaps?

The overlaps function checks if two geometries overlap. Two geometries overlap
if they have the same dimension, and if they overlap then their interiors
interact, but they both also need interior points exterior to the other
geometry. 

Note that this means it is impossible for a single point to overlap with a
single point and a line only overlaps with another line if only a section of
each line is colinear (crosses don't count for interior points interacting). 

To provide an example, consider these two lines:
```@example overlaps
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (0.0, 10.0)])
l2 = GI.LineString([(0.0, -10.0), (0.0, 3.0)])
f, a, p = lines(GI.getpoint(l1), color = :blue)
scatter!(GI.getpoint(l1), color = :blue)
lines!(GI.getpoint(l2), color = :orange)
scatter!(GI.getpoint(l2), color = :orange)
```
We can see that the two lines overlap in the plot:
```@example overlaps
overlap(l1, l2)
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`overlaps` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are allowed to be in the interior of g2
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are allowed to be in the exterior of g2
    - at least one point of g1 is required to be in the interior of g2
    - at least one point of g2 is required to be in the interior of g1
    - no points of g1 is required to be on the boundary of g2
    - at least one point of g1 is required to be in the exterior of g2
    - at least one point of g2 is required to be in the exterior of g1

The code for the specific implementations is in the geom_geom_processors file.
=#

"""
    overlaps(geom1, geom2)::Bool

Compare two Geometries of the same dimension and return true if their interiors
interact, but they both also have interior points exterior to the other
geometry. 

## Examples
```jldoctest
import GeometryOps as GO, GeoInterface as GI
poly1 = GI.Polygon([[(0,0), (0,5), (5,5), (5,0), (0,0)]])
poly2 = GI.Polygon([[(1,1), (1,6), (6,6), (6,1), (1,1)]])

GO.overlaps(poly1, poly2)
# output
true
```
"""
overlaps(g1, g2)::Bool = _overlaps(GI.trait(g1), g1, GI.trait(g2), g2)


# # Convert features to geometries
_overlaps(::GI.FeatureTrait, g1, ::Any, g2) = overlaps(GI.geometry(g1), g2)
_overlaps(::Any, g1, t2::GI.FeatureTrait, g2) = overlaps(g1, GI.geometry(g2))


# # Non-specified geometries 

# Geometries of different dimensions and points cannot overlap and return false
_overlaps(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g2,
) = false


# # Point disjoint geometries

# Point is disjoint from another point if the points are not equal.
_disjoint(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = !equals(g1, g2)


function overlaps(
    ::GI.MultiPointTrait, points1,
    ::GI.MultiPointTrait, points2,
)
    one_diff = false  # assume that all the points are the same
    one_same = false  # assume that all points are different
    for p1 in GI.getpoint(points1)
        match_point = false
        for p2 in GI.getpoint(points2)
            if equals(p1, p2)  # Point is shared
                one_same = true
                match_point = true
                break
            end
        end
        one_diff |= !match_point  # Point isn't shared
        one_same && one_diff && return true
    end
    return false
end

"""
    overlaps(::GI.LineTrait, line1, ::GI.LineTrait, line)::Bool

If the lines overlap, meaning that they are colinear but each have one endpoint
outside of the other line, return true. Else false.
"""
# overlaps(::GI.LineTrait, line1, ::GI.LineTrait, line) =
#     _overlaps((a1, a2), (b1, b2))

"""
    overlaps(
        ::Union{GI.LineStringTrait, GI.LinearRing}, line1,
        ::Union{GI.LineStringTrait, GI.LinearRing}, line2,
    )::Bool

If the curves overlap, meaning that at least one edge of each curve overlaps,
return true. Else false.
"""
# function overlaps(
#     ::Union{GI.LineStringTrait, GI.LinearRing}, line1,
#     ::Union{GI.LineStringTrait, GI.LinearRing}, line2,
# )
#     edges_a, edges_b = map(sort! ∘ to_edges, (line1, line2))
#     for edge_a in edges_a
#         for edge_b in edges_b
#             _overlaps(edge_a, edge_b) && return true
#         end
#     end
#     return false
# end
# function overlaps(
#     ::GI.LineStringTrait, g1,
#     ::GI.LineStringTrait, g2,
# )
#     cross, overlap = _line_curve_crosses_overlap_interactions(
#         g1, g2;
#         closed_line = false, closed_curve = false,
#     )
#     return !cross && overlap
# end

# function overlaps(
#     ::GI.LineStringTrait, g1,
#     ::GI.LinearRingTrait, g2,
# )
#     cross, overlap = _line_curve_crosses_overlap_interactions(
#         g1, g2;
#         closed_line = true, closed_curve = false,
#     )
#     return !cross && overlap
# end
# function overlaps(
#     ::GI.LinearRingTrait, g1,
#     ::GI.LineStringTrait, g2,
# )
#     cross, overlap = _line_curve_crosses_overlap_interactions(
#         g1, g2;
#         closed_line = false, closed_curve = true,
#     )
#     return !cross && overlap
# end
# function overlaps(
#     ::GI.LinearRingTrait, g1,
#     ::GI.LinearRingTrait, g2,
# )
#     cross, overlap = _line_curve_crosses_overlap_interactions(
#         g1, g2;
#         closed_line = true, closed_curve = true,
#     )
#     return !cross && overlap
# end

overlaps(
    ::GI.LineStringTrait, g1,
    ::GI.LineStringTrait, g2,
) = _line_curve_process(
        g1, g2;
        over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
        in_require = true, on_require = false, out_require = true,
        closed_line = false,
        closed_curve = false,
    ) && _line_curve_process(
            g2, g1;
            over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
            in_require = true, on_require = false, out_require = true,
            closed_line = false,
            closed_curve = false,
        )

overlaps(
    ::GI.LineStringTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
        g1, g2;
        over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
        in_require = true, on_require = false, out_require = true,
        closed_line = false,
        closed_curve = true,
    ) && _line_curve_process(
            g2, g1;
            over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
            in_require = true, on_require = false, out_require = true,
            closed_line = true,
            closed_curve = false,
        )

overlaps(
    ::GI.LinearRingTrait, g1,
    ::GI.LineStringTrait, g2,
) = _line_curve_process(
        g1, g2;
        over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
        in_require = true, on_require = false, out_require = true,
        closed_line = true,
        closed_curve = false,
    ) && _line_curve_process(
            g2, g1;
            over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
            in_require = true, on_require = false, out_require = true,
            closed_line = false,
            closed_curve = true,
        )

overlaps(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
        g1, g2;
        over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
        in_require = true, on_require = false, out_require = true,
        closed_line = true,
        closed_curve = true,
    ) && _line_curve_process(
            g2, g1;
            over_allow = true, cross_allow = false, on_allow = true, out_allow = true,
            in_require = true, on_require = false, out_require = true,
            closed_line = true,
            closed_curve = true,
        )

"""
    overlaps(
        trait_a::GI.PolygonTrait, poly_a,
        trait_b::GI.PolygonTrait, poly_b,
    )::Bool

If the two polygons intersect with one another, but are not equal, return true.
Else false.
"""
# function overlaps(
#     trait_a::GI.PolygonTrait, poly_a,
#     trait_b::GI.PolygonTrait, poly_b,
# )
#     edges_a, edges_b = map(sort! ∘ to_edges, (poly_a, poly_b))
#     return _line_intersects(edges_a, edges_b) &&
#         !equals(trait_a, poly_a, trait_b, poly_b)
# end

"""
    overlaps(
        ::GI.PolygonTrait, poly1,
        ::GI.MultiPolygonTrait, polys2,
    )::Bool

Return true if polygon overlaps with at least one of the polygons within the
multipolygon. Else false.
"""
function overlaps(
    ::GI.PolygonTrait, poly1,
    ::GI.MultiPolygonTrait, polys2,
)
    for poly2 in GI.getgeom(polys2)
        overlaps(poly1, poly2) && return true
    end
    return false
end

"""
    overlaps(
        ::GI.MultiPolygonTrait, polys1,
        ::GI.PolygonTrait, poly2,
    )::Bool

Return true if polygon overlaps with at least one of the polygons within the
multipolygon. Else false.
"""
overlaps(trait1::GI.MultiPolygonTrait, polys1, trait2::GI.PolygonTrait, poly2) = 
    overlaps(trait2, poly2, trait1, polys1)

"""
    overlaps(
        ::GI.MultiPolygonTrait, polys1,
        ::GI.MultiPolygonTrait, polys2,
    )::Bool

Return true if at least one pair of polygons from multipolygons overlap. Else
false.
"""
function overlaps(
    ::GI.MultiPolygonTrait, polys1,
    ::GI.MultiPolygonTrait, polys2,
)
    for poly1 in GI.getgeom(polys1)
        overlaps(poly1, polys2) && return true
    end
    return false
end

"""
    _overlaps(
        (a1, a2)::Edge,
        (b1, b2)::Edge
    )::Bool

If the edges overlap, meaning that they are colinear but each have one endpoint
outside of the other edge, return true. Else false. 
"""
function _overlaps(
    (a1, a2)::Edge,
    (b1, b2)::Edge
)
    # meets in more than one point or at endpoints
    on_top = ExactPredicates.meet(a1, a2, b1, b2) == 0
    on_top || return false
    # check that one endpoint of each edge is within other edge
    a1_in = point_segment_orientation(a1, b1, b2) == point_in
    a2_in = point_segment_orientation(a2, b1, b2) == point_in
    b1_in = point_segment_orientation(b1, a1, a2) == point_in
    b2_in = point_segment_orientation(b2, a1, a2) == point_in
    return (a1_in ⊻ a2_in) && (b1_in ⊻ b2_in)
end


function overlaps(
    ::GI.MultiPointTrait, points1,
    ::GI.MultiPointTrait, points2,
)
    one_diff = false  # assume that all the points are the same
    one_same = false  # assume that all points are different
    for p1 in GI.getpoint(points1)
        match_point = false
        for p2 in GI.getpoint(points2)
            if equals(p1, p2)  # Point is shared
                one_same = true
                match_point = true
                break
            end
        end
        one_diff |= !match_point  # Point isn't shared
        one_same && one_diff && return true
    end
    return false
end