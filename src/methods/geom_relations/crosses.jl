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
            if seg_val == line_over
                return false
            elseif seg_val == line_cross
                orientation_req_met = true
            elseif seg_val == line_hinge && !orientation_req_met
                _, fracs = _intersection_point(
                    (_tuple_point(l_start), _tuple_point(l_end)),
                    (_tuple_point(c_start), _tuple_point(c_end))
                )
                if !isnothing(fracs)
                    (α, β) = fracs  # 0 ≤ α ≤ 1 and 0 ≤ β ≤ 1 since hinge
                    if β == 0 # already checked on previous segment
                        c_start = c_end
                        continue
                    elseif α == 0
                        if !closed_line && i == 2
                            c_start = c_end
                            continue
                        end
                    elseif α == 1
                        if !closed_line && i == nl
                            c_start = c_end
                            continue
                        end
                    else # 0 < α < 1, β = 1 (if 0 < β < 1 then seg_val = cross)
                        if !closed_curve && j == nc
                            c_start = c_end
                            continue
                        end
                        c_next = GI.getpoint(curve, j < nc ? j + 1 : 1)
                        x_start, y_start = GI.x(c_start), GI.y(c_start)
                        x_next, y_next = GI.x(c_next), GI.y(c_next)
                        Δx = GI.x(l_end) - GI.x(l_start)
                        Δy = GI.y(l_end) - GI.y(l_start)
                        if Δx == 0
                            x = GI.x(l_start)
                            if (x_next - x) * (x_start - x) ≥ 0
                                c_start = c_end
                                continue
                            end
                        elseif Δy == 0
                            y = GI.y(l_start)
                            if (y_next - y) * (y_start - y) ≥ 0
                                c_start = c_end
                                continue
                            end
                        else
                            m = Δy / Δx
                            b = GI.y(l_start) - m * GI.x(l_start)
                            Δy_start = (m * x_start + b) - y_start
                            Δy_next = (m * x_next + b) - y_next
                            if Δy_start * Δy_next ≥ 0
                                c_start = c_end
                                continue
                            end
                        end
                        orientation_req_met = true
                        continue
                    end
                end
                T = typeof(GI.x(l_start))
                (α, β) =  # α = 0 or α = 1
                    if equals(l_start, c_start)
                        (zero(T), zero(T))
                    elseif equals(l_start, c_end)
                        (zero(T), one(T))
                    elseif equals(l_end, c_start)
                        (one(T), zero(T))
                    elseif equals(l_end, c_end)
                        (one(T), one(T))
                    else
                        fracs
                    end
                if β == 0
                    c_start = c_end 
                    continue
                end
                l1, l2, l3 = α == 0 ?
                    (GI.getpoint(line, i > 2 ? (i - 2) : nl), l_start, l_end) :
                    (l_start, l_end, GI.getpoint(line, i < nl ? (i + 1) : 1))
                θ1 = atan(GI.y(l1) - GI.y(l2), GI.x(l1) - GI.x(l2))
                θ2 = atan(GI.y(l3) - GI.y(l2), GI.x(l3) - GI.x(l2))
                θ1, θ2 = θ1 < θ2 ? (θ1, θ2) : (θ2, θ1)

                c_next = β == 1 ?
                    GI.getpoint(curve, j < nc ? j + 1 : 1) :
                    c_end
                ϕ1 = atan(GI.y(c_start) - GI.y(l2), GI.x(c_start) - GI.x(l2))
                ϕ2 = atan(GI.y(c_next) - GI.y(l2), GI.x(c_next) - GI.x(l2))
                orientation_req_met = ((θ1 < ϕ1 < θ2) ⊻ (θ1 < ϕ2 < θ2))
            end
            c_start = c_end
        end
        l_start = l_end
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
