# export point_in_geom, point_in_polygon

# """
#     point_in_geom(
#         point, geom;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the point is within the geometry (excluding edges and vertices).
# `On` means the point is on an edge or a vertex of the geometry.
# `Out` means the point is outside of the geometry.
# """
# point_in_geom(
#     point, geom;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T} = point_in_geom(
#     GI.trait(point), point,
#     GI.trait(geom), geom;
#     in = in, on = on, out = out,
# )

# """
#     line_in_geom(
#         line, geom;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the line is within geometry (no segments on edges and vertices).
# `On` means the line has at least one segment on a geometry edge or a vertex.
# `Out` means the line has at least one segment outside of the geometry.
# """
# line_in_geom(
#     line, geom;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T} = line_in_geom(
#     GI.trait(line), line,
#     GI.trait(geom), geom;
#     in = in, on = on, out = out,
# )

# # ring_in_geom(ring, geom) = ring_in_geom(
# #     GI.trait(ring), ring,
# #     GI.trait(geom), geom,
# # )

# """
#     point_in_geom(
#         ::GI.PointTrait, point,
#         ::GI.LineStringTrait, line;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# Note that a point can only be within a linestring if the linestring is closed,
# by having an explicilty repeated last point. Even then, this means the point is
# within the ring created by the linestring, not on the linestring itself.
# `In` means the point is within the linestring (excluding edges and vertices).
# `On` means the point is on an edge or a vertex of the linestring.
# `Out` means the point is outside of the linestring.
# """
# function point_in_geom(
#     ::GI.PointTrait, point,
#     ::GI.LineStringTrait, line;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     results = if equals(
#         GI.getpoint(line, 1),
#         GI.getpoint(line, GI.npoint(line)),
#     )
#         _point_in_extent(point, GI.extent(line)) || return out
#         _point_in_closed_curve(point, line; in = in, on = on, out = out)
#     else
#         @warn "Linestring isn't closed. Point cannot be 'in' linestring."
#         out
#     end
#     return results
# end

# """
#     line_in_geom(
#         ::GI.LineStringTrait, line1,
#         ::GI.LineStringTrait, line2;
#         in::T = 1, on::T = -1, out::T = 0,
#     )

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# Note that a linestring can only be within a linestring if that linestring is
# closed by having an explicilty repeated last point. Even then, this means the
# point is within the ring created by the linestring, not on the linestring
# itself.
# `In` means the line is within geometry (no segments on edges and vertices).
# `On` means the line has at least one segment on a geometry edge or a vertex.
# `Out` means the line has at least one segment outside of the geometry.
# """
# function line_in_geom(
#     ::GI.LineStringTrait, line1,
#     ::GI.LineStringTrait, line2;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     results = if equals(  # if line2 is closed by a repeated last point
#         GI.getpoint(line2, 1),
#         GI.getpoint(line2, GI.npoint(line2)),
#     )
#         Extents.intersects(
#             GI.extent(line1),
#             GI.extent(line2),
#         ) || return out
#         _line_in_closed_curve(
#             line1, line2;
#             close = false,
#         )
#     else
#         @warn "Linestring isn't closed. Point cannot be 'in' linestring."
#         out
#     end
#     return results
# end

# # function ring_in_geom(::GI.LinearRingTrait, ring, ::GI.LineStringTrait, line)
# #     results = if equals(
# #         GI.getpoint(line, 1),
# #         GI.getpoint(line, GI.npoint(line)),
# #     )
# #         Extents.intersects(
# #             GI.extent(ring),
# #             GI.extent(line),
# #         ) || return (false, false)
# #         _line_in_closed_curve(ring, line; close = true)
# #     else
# #         @warn "Linestring isn't closed. Point cannot be 'in' linestring."
# #         (false, false)
# #     end
# #     return results
# # end
# """
#     point_in_geom(
#         ::GI.PointTrait, point,
#         ::GI.LinearRingTrait, ring;
#         in::T = 1, on::T = -1, out::T = 0
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the point is within the linear ring (excluding edges and vertices).
# `On` means the point is on an edge or a vertex of the linear ring.
# `Out` means the point is outside of the linear ring.
# """
# function point_in_geom(
#     ::GI.PointTrait, point,
#     ::GI.LinearRingTrait, ring;
#     in::T = 1, on::T = -1, out::T = 0
# ) where {T}
#     _point_in_extent(point, GI.extent(ring)) || return out
#     return _point_in_closed_curve(point, ring; in = in, on = on, out = out)
# end

# """
#     line_in_geom(
#         ::GI.LineStringTrait, line,
#         ::GI.LinearRingTrait, ring;
#         in::T = 1, on::T = -1, out::T = 0,
#     )

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the line is within the ring (no segments on edges and vertices).
# `On` means the line has at least one segment on a ring edge or a vertex.
# `Out` means the line has at least one segment outside of the ring.
# """
# function line_in_geom(
#     ::GI.LineStringTrait, line,
#     ::GI.LinearRingTrait, ring;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     Extents.intersects(GI.extent(line), GI.extent(ring)) || return out
#     return _line_in_closed_curve(
#         line, ring;
#         close = false,
#     )
# end

# # function ring_in_geom(::GI.LinearRingTrait, ring1, ::GI.LinearRingTrait, ring2)
# #     Extents.intersects(GI.extent(ring1), GI.extent(ring2)) || return (false, false)
# #     _line_in_closed_curve(ring1, ring2; close = true)
# # end

# # function polygon_in_geom(::GI.PolygonTrait, poly, ::GI.LinearRingTrait, ring)
# #     Extents.intersects(GI.extent(poly), GI.extent(ring)) || return (false, false)
# #     return _line_in_closed_curve(GI.getexterior(poly), ring; close = true)
# # end

# """
#     point_in_geom(
#         ::GI.PointTrait, point,
#         ::GI.PolygonTrait, poly;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the point is within polygon (excluding edges, vertices, and holes).
# `On` means the point is on an edge or a vertex of the polygon.
# `Out` means the point is outside of the polygon, including within holes.
# """
# function point_in_geom(
#     ::GI.PointTrait, point,
#     ::GI.PolygonTrait, poly;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     _point_in_extent(point, GI.extent(poly)) || return out
#     ext_val = _point_in_closed_curve(
#         point, GI.getexterior(poly);
#         in = in, on = on, out = out,
#     )
#     ext_val == on && return ext_val
#     in_out_counter = (ext_val == in) ? 1 : 0
#     for ring in GI.gethole(poly)
#         hole_val = _point_in_closed_curve(
#             point, ring;
#             in = in, on = on, out = out,
#         )
#         hole_val == on && return hole_val
#         in_out_counter += (hole_val == in) ? 1 : 0
#     end
#     return iseven(in_out_counter) ? out : in
# end

# """
#     point_in_polygon(
#         point, polygon;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the point is within polygon (excluding edges, vertices, and holes).
# `On` means the point is on an edge or a vertex of the polygon.
# `Out` means the point is outside of the polygon, including within holes.
# """
# point_in_polygon(
#     point, polygon;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T} = point_in_polygon(
#     GI.trait(point), point,
#     GI.trait(polygon), polygon;
#     in = in, on = on, out = out,
# )

# """
#     point_in_polygon(
#         trait1::GI.PointTrait, point,
#         trait2::GI.PolygonTrait, poly;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the point is within polygon (excluding edges, vertices, and holes).
# `On` means the point is on an edge or a vertex of the polygon.
# `Out` means the point is outside of the polygon, including within holes.
    
# Note that this is the same as point_in_geom dispatched on a polygon. 
# """
# point_in_polygon(
#     trait1::GI.PointTrait, point,
#     trait2::GI.PolygonTrait, poly;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}= point_in_geom(
#     trait1, point,
#     trait2, poly;
#     in = in, on = on, out = out,
# )

# """
#     point_in_geom(
#         ::GI.LineStringTrait, line,
#         ::GI.PolygonTrait, poly;
#         in::T = 1, on::T = -1, out::T = 0,
#     )

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the line is within the polygon (no segments on edges, vertices, or
#     holes).
# `On` means the line has at least one segment on a polygon edge or a vertex.
# `Out` means the line has at least one segment outside of the polygon (including
#     within a hole).
# """
# function point_in_geom(
#     ::GI.LineStringTrait, line,
#     ::GI.PolygonTrait, poly;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     Extents.intersects(GI.extent(line), GI.extent(ring)) || return out
#     ext_val = _line_in_closed_curve(
#         line, GI.getexterior(poly);
#         close = false, 
#     )
    
#     for ring in GI.gethole(poly)
#         hole_val = _point_in_closed_curve(
#             point, ring;
#             in = in, on = on, out = out,
#         )
#         hole_val == on && return hole_val
#         in_out_counter += (hole_val == in) ? 1 : 0
#     end
#     return iseven(in_out_counter) ? out : in
# end

# # line_in_polygon(
# #     line, poly;
# #     in::T = 1, on::T = -1, out::T = 0,
# # ) where {T} = line_in_geom(
# #         line, GI.trait(line),
# #         poly, GI.trait(poly);
# #         in = in, on = on, out = out,
# #     )
 
# # ring_in_geom(::GI.LinearRingTrait, ring, ::GI.PolygonTrait, poly) = 
# #     _geom_in_polygon(ring, poly; close = true)
    
# # function polygon_in_geom(::GI.PolygonTrait, poly1, ::GI.PolygonTrait, poly2)
# #     # Cheaply check that the point is inside the polygon extent
# #     Extents.intersects(GI.extent(poly1), GI.extent(poly2)) || return (false, false)
# #     # Make sure exterior of poly1 is within exterior of poly2
# #     in_ext, some_on_ext = _line_in_closed_curve(
# #         GI.getexterior(poly1), GI.getexterior(poly2);
# #         close = true,
# #     )
# #     # poly1 not within poly2's external ring
# #     (in_ext || some_on_ext) || return (false, false)
# #     # Check if the geom is in any of the holes
# #     outside_hole, some_on_hole = true, false
# #     for hole in GI.gethole(poly)
# #         outside_hole, some_on_hole = _line_in_closed_curve(
# #             geom, hole;
# #             close = close, in = false,
# #         )
# #         # geom is in a hole -> not in polygon
# #         !(outside_hole || some_on_hole) && return (false, false)
# #     end
# #     return (in_ext && outside_hole, some_on_hole || some_on_ext)  # geom is inside of polygon
# # end

# # function polygon_in_polygon(poly1, poly2)
# #     # edges1, edges2 = to_edges(poly1), to_edges(poly2)
# #     # extent1, extent2 = to_extent(edges1), to_extent(edges2)
# #     # Check the extents intersect
# #     Extents.intersects(GI.extent(poly1), GI.extent(poly2)) || return false

# #     # Check all points in poly1 are in poly2
# #     for point in GI.getpoint(poly1)
# #         point_in_polygon(point, poly2) || return false
# #     end

# #     # Check the line of poly1 does not intersect the line of poly2
# #     intersects(poly1, poly2) && return false

# #     # poly1 must be in poly2
# #     return true
# #  end

# """
#     _point_in_closed_curve(
#         point, curve;
#         in::T = 1, on::T = -1, out::T = 0,
#     )::T

# Determine if point is in, on, or out of a closed curve. Point should be an
# object of Point trait and curve should be a linearstring or ring, that is
# assumed to be closed, regardless of repeated last point.

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means the point is within the closed curve (excluding edges and vertices).
# `On` means the point is on an edge or a vertex of the closed curve.
# `Out` means the point is outside of the closed curve.

# Note that this uses the Algorithm by Hao and Sun (2018):
# https://doi.org/10.3390/sym10100477
# Paper seperates orientation of point and edge into 26 cases. For each case, it
# is either a case where the point is on the edge (returns on), where a ray from
# the point (x, y) to infinity along the line y = y cut through the edge (k += 1),
# or the ray does not pass through the edge (do nothing and continue). If the ray
# passes through an odd number of edges, it is within the curve, else outside of
# of the curve if it didn't return 'on'.
# See paper for more information on cases denoted in comments.
# """
# function _point_in_closed_curve(
#     point, curve;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     x, y = GI.x(point), GI.y(point)
#     n = GI.npoint(curve)
#     n -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, n)) ? 1 : 0
#     k = 0  # counter for ray crossings
#     p_start = GI.getpoint(curve, n)
#     @inbounds for i in 1:n
#         p_end = GI.getpoint(curve, i)
#         v1 = GI.y(p_start) - y
#         v2 = GI.y(p_end) - y
#         if !((v1 < 0 && v2 < 0) || (v1 > 0 && v2 > 0)) # if not cases 11 or 26
#             u1 = GI.x(p_start) - x
#             u2 = GI.x(p_end) - x
#             f = u1 * v2 - u2 * v1
#             if v2 > 0 && v1 ≤ 0                # Case 3, 9, 16, 21, 13, or 24
#                 f == 0 && return on            # Case 16 or 21
#                 f > 0 && (k += 1)              # Case 3 or 9
#             elseif v1 > 0 && v2 ≤ 0            # Case 4, 10, 19, 20, 12, or 25
#                 f == 0 && return on            # Case 19 or 20
#                 f < 0 && (k += 1)              # Case 4 or 10
#             elseif v2 == 0 && v1 < 0           # Case 7, 14, or 17
#                 f == 0 && return on            # Case 17
#             elseif v1 == 0 && v2 < 0           # Case 8, 15, or 18
#                 f == 0 && return on            # Case 18
#             elseif v1 == 0 && v2 == 0          # Case 1, 2, 5, 6, 22, or 23
#                 u2 ≤ 0 && u1 ≥ 0 && return on  # Case 1
#                 u1 ≤ 0 && u2 ≥ 0 && return on  # Case 2
#             end
#         end
#         p_start = p_end
#     end
#     return iseven(k) ? out : in
# end

# """
#     line_in_closed_curve(
#         line, curve;
#         in::T = 1, on::T = -1, out::T = 0,
#         close = false,
#     )

# Determine if line is in, on, or out of a closed curve. Both the line and curve
# should be an object of linestring or linearring trait. The curve is assumed to
# be closed, regardless of repeated last point.

# Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
# `In` means line is within the curve (no segments on edges, vertices, or holes).
# `On` means line has at least one segment on a curve edge or vertex.
# `Out` means the line has at least one segment outside of the curve.

# This algorithm functions by checking if the first point of the line is within
# the curve. If not, then the line is not within the curve, if so, we check for
# intersections between the line and curve, as this would mean a part of the line
# is outside of the curve. We take special care of intersections through vertices
# as it isn't clearcut if those neccesitate a segment of the line being outside
# of the curve.
# """
# _line_in_closed_curve(line, curve;
#         exclude_boundaries = false,
#         close = false,
# ) = _line_in_out_closed_curve(
#     line, curve;
#     disjoint = false,
#     exclude_boundaries = exclude_boundaries,
#     close = close,
# )

# # function _line_in_closed_curve(
# #     line, curve;
# #     in::T = 1, on::T = -1, out::T = 0,
# #     close = false,
# # ) where {T}
# #     # Determine number of points in curve and line
# #     nc = GI.npoint(curve)
# #     nc -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc)) ? 1 : 0
# #     nl = GI.npoint(line)
# #     nl -= (close && equals(GI.getpoint(line, 1), GI.getpoint(line, nl))) ? 1 : 0
# #     # Check to see if first point in line is within curve
# #     point_val = _point_in_closed_curve(
# #         GI.getpoint(line, 1), curve;
# #         in = in, on = on, out = out,
# #     )
# #     # point is outside curve, line can't be within curve
# #     point_val == out && return out
# #     # Check for any intersections between line and curve
# #     line_on_curve = point_val == on  # record if line is "on" part of curve
# #     l_start = _tuple_point(GI.getpoint(line, close ? nl : 1))
# #     for i in (close ? 1 : 2):nl
# #         l_end = _tuple_point(GI.getpoint(line, i))
# #         c_start = _tuple_point(GI.getpoint(curve, nc))
# #         for j in 1:nc
# #             c_end = _tuple_point(GI.getpoint(curve, j))
# #             # Check if edges intersect --> line is not within curve
# #             meet_type = ExactPredicates.meet(l_start, l_end, c_start, c_end)
# #             # open line segments meet in a single point
# #             meet_type == 1 && return out
# #             #=
# #             closed line segments meet in one or several points -> meet at a
# #             vertex or on the edge itself (parallel)
# #             =#
# #             if meet_type == 0
# #                 line_on_curve = true
# #                 # See if segment is parallel and within curve edge
# #                 p1_on_seg = point_on_segment(l_start, c_start, c_end)
# #                 p2_on_seg = point_on_segment(l_end, c_start, c_end)
# #                 # if segment isn't contained within curve edge
# #                 if !p1_on_seg || !p2_on_seg 
# #                     # Make sure l_start is in or on the segment
# #                     p1_in_curve =
# #                         p1_on_seg ||
# #                         _point_in_closed_curve(
# #                             l_start, curve;
# #                             in = in, on = on, out = out,
# #                         ) != out
# #                     !p1_in_curve && return out
# #                     # Make sure l_end is in or on the segment
# #                     p2_in_curve =
# #                         p2_on_seg ||
# #                         _point_in_closed_curve(
# #                             l_end, curve;
# #                             in = in, on = on, out = out,
# #                         ) != out
# #                     !p2_in_curve && return out
# #                     #=
# #                     If both endpoints are within or on the curve, but not
# #                     parallel to the edge, make sure that midpoints between the
# #                     intersections along the segment are within curve
# #                     =# 
# #                     !_segment_mids_in_curve(
# #                         l_start, l_end, curve;
# #                         in = in, on = on, out = out,
# #                     ) && return out  # point of segment is outside of curve
# #                     # line segment is fully within or on curve 
# #                     break 
# #                 end
# #             end
# #             c_start = c_end
# #         end
# #         l_start = l_end
# #     end
# #     # check if line is on any curve edges or vertcies
# #     return line_on_curve ? on : in
# # end

# """
#     _segment_mids_in_curve(
#         l_start, l_end, curve;
#         in::T = 1, on::T = -1, out::T = 0,
#     )

#     Given two points defining a line segment (both with point traits) and a
#     curve (with a linestring or linearring trait), find the intersection points
#     between them and sort them along the segment. Then, make sure that the
#     midpoint between pairs of points along the line is within the curve. Returns
#     true if all of the midpoints are within or on the curve and false otherwise.

#     Note: This function assumes that both of the endpoints of the line segment
#     are on the curve!
# """
# function _segment_mids_in_curve(
#     l_start, l_end, curve;
#     in::T = 1, on::T = -1, out::T = 0,
# ) where {T}
#     # Find intersection points
#     ipoints = intersection_points(
#         GI.Line([l_start, l_end]),
#         curve
#     )
#     npoints = length(ipoints)
#     if npoints < 3  # only intersection points are the endpoints
#         mid_val = _point_in_closed_curve(
#             (l_start .+ l_end) ./ 2, curve;
#             in = in, on = on, out = out,
#         )
#         mid_val == out && return false
#     else  # more intersection points than the endpoints
#         # sort intersection points along the line
#         sort!(ipoints, by = p -> euclid_distance(p, l_start))
#         p_start = ipoints[1]
#         for i in 2:npoints
#             p_end = ipoints[i]
#             # check if midpoint of intersection points is within the curve
#             mid_val = _point_in_closed_curve(
#                 (p_start .+ p_end) ./ 2, curve;
#                 in = in, on = on, out = out,
#             )
#             # if it is out, return false
#             mid_val == out && return false
#             p_start = p_end
#         end
#     end
#     return true  # all intersection point midpoints were in or on the curve
# end

# function _geom_in_polygon(geom, poly; close = false)
#     # Cheaply check that the geom extent is inside the polygon extent
#     Extents.intersects(GI.extent(geom), GI.extent(poly)) || return (false, false)
#     # Check if geom is inside or on the exterior ring
#     in_ext, on_ext = _line_in_closed_curve(
#         geom,
#         GI.getexterior(poly);
#         close = close,
#     )
#     (in_ext || on_ext) || return (false, false)  # geom isn't in external ring
#     # Check if the geom is in any of the holes
#     for hole in GI.gethole(poly)
#         # out_of_hole, some_on_hole = _line_in_closed_curve(
#         #     geom, hole;
#         #     close = close, in = false,
#         # )
#         # geom is in a hole -> not in polygon
#         !(out_of_hole || some_on_hole) && return (false, false)
#     end
#     return (in_ext, some_on_hole)  # geom is inside of polygon
# end

# """
#     _point_in_extent(p, extent::Extents.Extent)::Bool

# Returns true if the point is the bounding box of the extent and false otherwise. 
# """
# function _point_in_extent(p, extent::Extents.Extent)
#     (x1, x2), (y1, y2) = extent.X, extent.Y
#     return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
# end

# function _line_in_out_closed_curve(
#     line, curve;
#     disjoint = false,
#     exclude_boundaries = false,
#     close = false,
# )
#     #=
#     Set variables based off if we are determining within or disjoint.
#     If `_point_in_closed_curve` returns `true_orientation` it is on the right
#     side of the curve for the check. If it returns `false_orientation`, it is
#     on the wrong side of the curve for the check.
#     =#
#     false_orientation = disjoint ? 1 : 0 # if checking within, want points in
#     on = -1  # as used for point in closed curve

#     # Determine number of points in curve and line
#     nc = GI.npoint(curve)
#     nc -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc)) ? 1 : 0
#     nl = GI.npoint(line)
#     nl -= (close && equals(GI.getpoint(line, 1), GI.getpoint(line, nl))) ? 1 : 0

#     # Check to see if first point in line is within curve
#     point_val = _point_in_closed_curve(GI.getpoint(line, 1), curve)
#     # point is out (for within) or in curve (for disjoint) -> wrong orientation
#     point_val == false_orientation && return false
#     # point is on boundary and don't want boundary points -> wrong orientation
#     exclude_boundaries && point_val == on && return false

#     # Check for any intersections between line and curve
#     l_start = _tuple_point(GI.getpoint(line, close ? nl : 1))
#     for i in (close ? 1 : 2):nl
#         l_end = _tuple_point(GI.getpoint(line, i))
#         c_start = _tuple_point(GI.getpoint(curve, nc))
#         for j in 1:nc
#             c_end = _tuple_point(GI.getpoint(curve, j))
#             # Check if edges intersect --> line crosses --> wrong orientation
#             meet_type = ExactPredicates.meet(l_start, l_end, c_start, c_end)
#             # open line segments meet in a single point
#             meet_type == 1 && return false
#             #=
#             closed line segments meet in one or several points -> meet at a
#             vertex or on the edge itself (parallel)
#             =#
#             if meet_type == 0
#                 # See if segment is parallel and within curve edge
#                 p1_on_seg = point_on_segment(l_start, c_start, c_end)
#                 exclude_boundaries && p1_on_seg && return false
#                 p2_on_seg = point_on_segment(l_end, c_start, c_end)
#                 exclude_boundaries && p2_on_seg && return false
#                 # if segment isn't contained within curve edge
#                 if !p1_on_seg || !p2_on_seg 
#                     # Make sure l_start is in corrent orientation
#                     p1_val = p1_on_seg ?
#                         on :
#                         _point_in_closed_curve(l_start, curve)
#                     p1_val == false_orientation && return false
#                     exclude_boundaries && p1_val == on && return false
#                     # Make sure l_end is in is in corrent orientation
#                     p2_val = p2_on_seg ?
#                         on :
#                         _point_in_closed_curve(l_end, curve)
#                     p2_val == false_orientation && return false
#                     exclude_boundaries && p2_val == on && return false
#                     #=
#                     If both endpoints are in the correct orientation, but not
#                     parallel to the edge, make sure that midpoints between the
#                     intersections along the segment are also in the correct
#                     orientation
#                     =# 
#                     !_segment_mids_in_out_curve(
#                         l_start, l_end, curve;
#                         disjoint = disjoint,
#                         exclude_boundaries = exclude_boundaries,
#                     ) && return false  # midpoint on the wrong side of the curve
#                     # line segment is fully within or on curve 
#                     break 
#                 end
#             end
#             c_start = c_end
#         end
#         l_start = l_end
#     end
#     # check if line is on any curve edges or vertcies
#     return true
# end

# function _segment_mids_in_out_curve(
#     l_start, l_end, curve;
#     disjoint = false,
#     exclude_boundaries = false,
# )
#     false_orientation = disjoint ? 1 : 0 # if checking within, want points in
#     on = -1  # as used for point in closed curve
#     # Find intersection points
#     ipoints = intersection_points(
#         GI.Line([l_start, l_end]),
#         curve
#     )
#     npoints = length(ipoints)
#     if npoints < 3  # only intersection points are the endpoints
#         mid_val = _point_in_closed_curve(
#             (l_start .+ l_end) ./ 2, curve;
#             in = in, on = on, out = out,
#         )
#         mid_val == false_orientation && return false
#         exclude_boundaries && mid_val == on && return false
#     else  # more intersection points than the endpoints
#         # sort intersection points along the line
#         sort!(ipoints, by = p -> euclid_distance(p, l_start))
#         p_start = ipoints[1]
#         for i in 2:npoints
#             p_end = ipoints[i]
#             # check if midpoint of intersection points is within the curve
#             mid_val = _point_in_closed_curve(
#                 (p_start .+ p_end) ./ 2, curve;
#                 in = in, on = on, out = out,
#             )
#             # if it is out, return false
#             mid_val == false_orientation && return false
#             exclude_boundaries && mid_val == on && return false
#         end
#     end
#     return true  # all intersection point midpoints were in or on the curve
# end
