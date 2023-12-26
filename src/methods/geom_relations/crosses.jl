# # Crosses

export crosses

#=
## What is crosses?

The crosses function checks if one geometry is crosses another geometry.
A geometry can only cross another geometry if they are either two lines, or if
one of the geometries has a smaller dimensionality than the other geometry.
If checking two lines, they must meet in one point. If checking two geometries
of different dimensions, the interiors must meet in at least one point and at
least one of the geometries must have a point outside of the other geometry.

Note that points can't cross any geometries, despite different dimension, due to
their inability to be both crosses and exterior to any other shape.

To provide an example, consider these two lines:
```@example cshape
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie


```

```@example cshape

```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

...

The code for the specific implementations is in the geom_geom_processors file.

=#

"""
    crosses(geom1, geom2)::Bool

Return `true` if the first geometry crosses the second geometry. If they are
both lines, they must meet in one point. Otherwise, they must be of different
dimensions, the interiors must intersect, and the interior of the first geometry
must intersect the exterior of the secondary geometry.

## Examples
```jldoctest setup=:(using GeometryOps, GeometryBasics)
import GeometryOps as GO, GeoInterface as GI



# output

```
"""
crosses(g1, g2) = crosses(trait(g1), g1, trait(g2), g2)
crosses(::GI.FeatureTrait, g1, ::Any, g2) = crosses(GI.geometry(g1), g2)
crosses(::Any, g1, t2::GI.FeatureTrait, g2) = crosses(g1, GI.geometry(g2))

"""

"""
crosses(::GI.AbstractTrait, g1, ::GI.AbstractTrait, g2) = false

# Lines crosses geometries
"""
    crosses(::GI.LineStringTrait, g1, ::GI.LineStringTrait, g2)::Bool

A line string is crosses another linestring if the vertices and edges of the
first linestring are crosses the second linestring, including the first and last
vertex. Return true if those conditions are met, else false.
"""
crosses(
    ::GI.LineStringTrait, g1,
    ::GI.LineStringTrait, g2,
) = _line_curve_crosses_overlap_process(
        g1, g2;
        orientation = line_cross,
        closed_line = false, closed_curve = false,
    )

crosses(
    ::GI.LineStringTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_crosses_overlap_process(
        g1, g2;
        orientation = line_cross,
        closed_line = false, closed_curve = true,
    )

crosses(
    ::GI.LinearRingTrait, g1,
    ::GI.LineStringTrait, g2,
) = _line_curve_crosses_overlap_process(
        g1, g2;
        orientation = line_cross,
        closed_line = true, closed_curve = false,
    )

crosses(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_crosses_overlap_process(
    g1, g2;
    orientation = line_cross,
    closed_line = true, closed_curve = true,
)

function _line_curve_crosses_overlap_process(
    line, curve;
    orientation = line_cross,
    closed_line = false, closed_curve = false,
)
    nl = GI.npoint(line)
    nc = GI.npoint(curve)
    first_last_equal_line = equals(GI.getpoint(line, 1), GI.getpoint(line, nl))
    first_last_equal_curve = equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc))
    nl -= first_last_equal_line ? 1 : 0
    nc -= first_last_equal_curve ? 1 : 0
    closed_line |= first_last_equal_line
    closed_curve |= first_last_equal_curve

    # Loop over each line segment
    orientation_req_met = false
    l_start = GI.getpoint(line, closed_line ? nl : 1)
    for i in (closed_line ? 1 : 2):nl
        l_end = GI.getpoint(line, i)
        c_start = GI.getpoint(curve, closed_curve ? nc : 1)
        for j in (closed_curve ? 1 : 2):nc
            c_end = GI.getpoint(curve, j)
            seg_val = _segment_segment_orientation(
                (l_start, l_end),
                (c_start, c_end),
            )
            @show seg_val
            if seg_val == line_over
                return false
            elseif seg_val == line_cross
                orientation_req_met = true
            elseif seg_val == line_hinge && !orientation_req_met
                _, fracs = _intersection_point(
                    (_tuple_point(l_start), _tuple_point(l_end)),
                    (_tuple_point(c_start), _tuple_point(c_end))
                )
                if isnothing(fracs)  # line and curve segments are parallel
                    
                else
                    (α, β) = fracs  # 0 ≤ α ≤ 1 and 0 ≤ β ≤ 1 since hinges
                    β == 0 && break  # if crosses, found on previous segment
                    # curve intersects through line endpoint
                    if !closed_line && (
                        (i == 2 && α == 0) || (i == nl && α == 1)
                    )
                        break  # doesn't cross
                    # curve segment intersects through line vertex or edge
                    elseif 0 < β < 1
                        orientation_req_met = true  # crosses
                    # curve meets line at curve segment endpoint
                    else  # β == 1
                        # no curve segment connecting at intersection point
                        !closed_curve && j == nc && break
                        # see if next curve segment is on the other side of line
                        c_next = GI.getpoint(curve, j < nc ? j + 1 : 1)
                        Δx = GI.x(l_end) - GI.x(l_start)
                        Δy = GI.y(l_end) - GI.y(l_start)
                        if Δx == 0
                            x = GI.x(l_start)
                            x_next = GI.x(c_next)
                            x_next == x && break # next curve segment is on line
                            x_next < x && GI.x(c_start) < x && break
                            x_next > x && GI.x(c_start) > x && break
                            orientation_req_met = true
                        elseif Δy == 0
                            y = GI.y(l_start)
                            y_next = GI.y(c_next)
                            y_next == y && break # next curve segment is on line
                            y_next < y && GI.y(c_start) < y && break
                            y_next > y && GI.y(c_start) > y && break
                            orientation_req_met = true
                        else
                            m = Δy / Δx
                            b = GI.y(c_start) - m * GI.x(c_start)
                            Δy_start = (m * GI.x(c_start) + b) - GI.y(c_start)
                            Δy_next = (m * GI.x(c_next) + b) - GI.y(c_next)
                            Δy_start * Δy_next > 0 && break
                            orientation_req_met = true
                        end

                    end
            end
        end
    end
    return orientation_req_met
end

"""
    crosses(::GI.LineStringTrait, g1, ::GI.LinearRingTrait, g2)::Bool

A line string is crosses a linear ring if the vertices and edges of the
linestring are crosses the linear ring. Return true if those conditions are met,
else false.
"""
# crosses(
#     ::GI.LineStringTrait, g1,
#     ::GI.LinearRingTrait, g2,
# ) = _line_curve_process(
#     g1, g2;
#     in_allow = false, on_allow = true, out_allow = true,
#     in_require = false, on_require = true, out_require = true,
#     closed_line = false,
#     closed_curve = true,
# )

"""
    crosses(::GI.LineStringTrait, g1, ::GI.PolygonTrait, g2)::Bool

A line string is crosses a polygon if the vertices and edges of the
linestring are crosses the polygon. Points of the linestring can be on the
polygon edges, but at least one point must be in the polygon interior. The
linestring also cannot cross through a hole. Return true if those conditions are
met, else false.
"""
crosses(
    ::GI.LineStringTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    in_allow =  false, on_allow = true, out_allow = true,
    in_require = false, on_require = true, out_require = true,
    closed_line = false,
)



"""
     crosses(geom1, geom2)::Bool

Return `true` if the intersection results in a geometry whose dimension is one less than
the maximum dimension of the two source geometries and the intersection set is interior to
both source geometries.

TODO: broken

## Examples 
```julia
import GeoInterface as GI, GeometryOps as GO

line1 = GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])
line2 = GI.LineString([(-2, 2), (4, 2)])

GO.crosses(line1, line2)
# output
true
```
"""
# crosses(g1, g2)::Bool = crosses(trait(g1), g1, trait(g2), g2)::Bool
# crosses(t1::FeatureTrait, g1, t2, g2)::Bool = crosses(GI.geometry(g1), g2)
# crosses(t1, g1, t2::FeatureTrait, g2)::Bool = crosses(g1, geometry(g2))
# crosses(::MultiPointTrait, g1, ::LineStringTrait, g2)::Bool = multipoint_crosses_line(g1, g2)
# crosses(::MultiPointTrait, g1, ::PolygonTrait, g2)::Bool = multipoint_crosses_poly(g1, g2)
# crosses(::LineStringTrait, g1, ::MultiPointTrait, g2)::Bool = multipoint_crosses_lines(g2, g1)
# crosses(::LineStringTrait, g1, ::PolygonTrait, g2)::Bool = line_crosses_poly(g1, g2)
# crosses(::LineStringTrait, g1, ::LineStringTrait, g2)::Bool = line_crosses_line(g1, g2)
# crosses(::PolygonTrait, g1, ::MultiPointTrait, g2)::Bool = multipoint_crosses_poly(g2, g1)
# crosses(::PolygonTrait, g1, ::LineStringTrait, g2)::Bool = line_crosses_poly(g2, g1)

# # function multipoint_crosses_line(geom1, geom2)
# #     int_point = false
# #     ext_point = false
# #     i = 1
# #     np2 = GI.npoint(geom2)

# #     while i < GI.npoint(geom1) && !int_point && !ext_point
# #         for j in 1:GI.npoint(geom2) - 1
# #             exclude_boundary = (j === 1 || j === np2 - 2) ? :none : :both
# #             if point_on_segment(GI.getpoint(geom1, i), (GI.getpoint(geom2, j), GI.getpoint(geom2, j + 1)); exclude_boundary)
# #                 int_point = true
# #             else
# #                 ext_point = true
# #             end
# #         end
# #         i += 1
# #     end

# #     return int_point && ext_point
# # end

# function line_crosses_line(line1, line2)
#     np2 = GI.npoint(line2)
#     if intersects(line1, line2)
#         for i in 1:GI.npoint(line1) - 1
#             for j in 1:GI.npoint(line2) - 1
#                 exclude_boundary = (j === 1 || j === np2 - 2) ? :none : :both
#                 pa = GI.getpoint(line1, i)
#                 pb = GI.getpoint(line1, i + 1)
#                 p = GI.getpoint(line2, j)
#                 te(p, (pa, pb); exclude_boundary) && return true
#             end
#         end
#     end
#     return false
# end

# function line_crosses_poly(line, poly)
#     for l in flatten(AbstractCurveTrait, poly)
#         intersects(line, l) && return true
#     end
#     return false
# end

# function multipoint_crosses_poly(mp, poly)
#     int_point = false
#     ext_point = false

#     for p in GI.getpoint(mp)
#         if point_in_polygon(p, poly)
#             int_point = true
#         else
#             ext_point = true
#         end
#         int_point && ext_point && return true
#     end
#     return false
# end
