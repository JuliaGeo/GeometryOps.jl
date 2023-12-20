@enum PointOrientation point_in=1 point_on=2 point_out=3

@enum LineOrientation line_cross=1 line_hinge=2 line_over=3 line_out=4

"""
    _point_curve_process(
        point, curve;
        process::ProcessType = within_process,
        exclude_boundaries = false,
        repeated_last_coord = false,
    )::Bool

Determines if a point meets the given process checks with respect to a curve.
This curve includes just the line segments that make up the curve. Even if the
curve has a repeated last point that "closes" the curve, this function does not
include the space within the closed curve as a part of the geometry. Point
should be an object of Point trait and curve should be an object with a line
string or a linear ring trait.

If checking within, then the point must be on a segment of the line and if
checking disjoint the point must not be on any segment of the line.

Beyond specifying the process type, user can also specify if the geometry
boundaries should be included in the checks and if the curve should be closed
with repeated a repeated last coordinate matching the first coordinate.
"""
function _point_curve_process(
    point, curve;
    in_allow, on_allow, out_allow,
    repeated_last_coord = false,
)
    n = GI.npoint(curve)
    first_last_equal = equals(GI.getpoint(curve, 1), GI.getpoint(curve, n))
    repeated_last_coord |= first_last_equal
    n -= first_last_equal ? 1 : 0
    # Loop through all curve segments
    p_start = GI.getpoint(curve, repeated_last_coord ? n : 1)
    @inbounds for i in (repeated_last_coord ? 1 : 2):n
        p_end = GI.getpoint(curve, i)
        seg_val = point_segment_orientation(point, p_start, p_end)
        seg_val == point_in && return in_allow
        if seg_val == point_on
            if !repeated_last_coord
                i == 2 && equals(point, p_start) && return on_allow
                i == n && equals(point, p_end) && return on_allow
            end
            return in_allow
        end
        p_start = p_end
    end
    return out_allow
end

"""
    point_segment_orientation(
        point::Point, start::Point, stop::Point;
        in::T = 1, on::T = -1, out::T = 0,
    )::T where {T}

Determines if a point is in, on, or out of a segment. If the point is 'on' the
segment it is on one of the segments endpoints. If it is 'in', it is on any
other point of the segment. If the point is not on any part of the segment, it
is 'out' of the segment.
"""
function point_segment_orientation(
    point, start, stop;
    in::T = point_in, on::T = point_on, out::T = point_out,
) where {T}
    # Parse out points
    x, y = GI.x(point), GI.y(point)
    x1, y1 = GI.x(start), GI.y(start)
    x2, y2 = GI.x(stop), GI.y(stop)
    Δx_seg = x2 - x1
    Δy_seg = y2 - y1
    Δx_pt = x - x1
    Δy_pt = y - y1
    if (Δx_pt == 0 && Δy_pt == 0) || (Δx_pt == Δx_seg && Δy_pt == Δy_seg)
        # If point is equal to the segment start or end points
        return on
    else
        #=
        Determine if the point is on the segment -> see if vector from segment
        start to point is parallel to segment and if point is between the
        segment endpoints
        =#
        on_line = _isparallel(Δx_seg, Δy_seg, Δx_pt, Δy_pt)
        !on_line && return out
        between_endpoints =
            (x2 > x1 ? x1 <= x <= x2 : x2 <= x <= x1) &&
            (y2 > y1 ? y1 <= y <= y2 : y2 <= y <= y1)
        !between_endpoints && return out
    end
    return in
end

"""
    point_filled_curve_orientation(
        point, curve;
        in::T = 1, on::T = -1, out::T = 0,
    )::T where {T}

Determine if point is in, on, or out of a closed curve, which includes the space
enclosed by the closed curve. Point should be an object of Point trait and
curve should be an object with a line string or linear ring trait, that is
assumed to be closed, regardless of repeated last point.

Returns a given in, on, or out value (defaults are 1, -1, and 0) of type T.
`In` means the point is within the closed curve (excluding edges and vertices).
`On` means the point is on an edge or a vertex of the closed curve.
`Out` means the point is outside of the closed curve.

Note that this uses the Algorithm by Hao and Sun (2018):
https://doi.org/10.3390/sym10100477
Paper seperates orientation of point and edge into 26 cases. For each case, it
is either a case where the point is on the edge (returns on), where a ray from
the point (x, y) to infinity along the line y = y cut through the edge (k += 1),
or the ray does not pass through the edge (do nothing and continue). If the ray
passes through an odd number of edges, it is within the curve, else outside of
of the curve if it didn't return 'on'.
See paper for more information on cases denoted in comments.
"""
function point_filled_curve_orientation(
    point, curve;
    in::T = point_in, on::T = point_on, out::T = point_out,
) where {T}
    x, y = GI.x(point), GI.y(point)
    n = GI.npoint(curve)
    n -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, n)) ? 1 : 0
    k = 0  # counter for ray crossings
    p_start = GI.getpoint(curve, n)
    @inbounds for i in 1:n
        p_end = GI.getpoint(curve, i)
        v1 = GI.y(p_start) - y
        v2 = GI.y(p_end) - y
        if !((v1 < 0 && v2 < 0) || (v1 > 0 && v2 > 0)) # if not cases 11 or 26
            u1 = GI.x(p_start) - x
            u2 = GI.x(p_end) - x
            c1 = u1 * v2  # first element of cross product summation
            c2 = u2 * v1  # second element of cross product summation
            f = c1 - c2
            if v2 > 0 && v1 ≤ 0                # Case 3, 9, 16, 21, 13, or 24
                (c1 ≈ c2) && return on         # Case 16 or 21
                f > 0 && (k += 1)              # Case 3 or 9
            elseif v1 > 0 && v2 ≤ 0            # Case 4, 10, 19, 20, 12, or 25
                (c1 ≈ c2) && return on         # Case 19 or 20
                f < 0 && (k += 1)              # Case 4 or 10
            elseif v2 == 0 && v1 < 0           # Case 7, 14, or 17
                (c1 ≈ c2) && return on         # Case 17
            elseif v1 == 0 && v2 < 0           # Case 8, 15, or 18
                (c1 ≈ c2) && return on         # Case 18
            elseif v1 == 0 && v2 == 0          # Case 1, 2, 5, 6, 22, or 23
                u2 ≤ 0 && u1 ≥ 0 && return on  # Case 1
                u1 ≤ 0 && u2 ≥ 0 && return on  # Case 2
            end
        end
        p_start = p_end
    end
    return iseven(k) ? out : in
end

"""
    _point_polygon_process(
        point, polygon;
        process::ProcessType = within_process,
        exclude_boundaries = false,
    )::Bool

Determines if a point meets the given process checks with respect to a polygon,
which excludes any holes specified by the polygon. Point should be an
object of Point trait and polygon should an object with a Polygon trait.
"""
function _point_polygon_process(
    point, polygon;
    in_allow, on_allow, out_allow,
)
    # Check interaction of geom with polygon's exterior boundary
    ext_val = point_filled_curve_orientation(point, GI.getexterior(polygon))
    # If a point is outside, it isn't interacting with any holes
    ext_val == point_out && return out_allow
    # if a point is on an external boundary, it isn't interacting with any holes
    ext_val == point_on && return on_allow
    
    # If geom is within the polygon, need to check interactions with holes
    for hole in GI.gethole(polygon)
        hole_val = point_filled_curve_orientation(point, hole)
        # If a point in in a hole, it is outside of the polygon
        hole_val == point_in && return out_allow
        # If a point in on a hole edge, it is on the edge of the polygon
        hole_val == point_on && return on_allow
    end
    
    # Point is within external boundary and on in/on any holes
    return in_allow
end

function _segment_segment_orientation(
    (a_point, b_point), (c_point, d_point);
    cross::T = line_cross, hinge::T = line_hinge,
    over::T = line_over, out::T = line_out,
) where T
    (ax, ay) = _tuple_point(a_point)
    (bx, by) = _tuple_point(b_point)
    (cx, cy) = _tuple_point(c_point)
    (dx, dy) = _tuple_point(d_point)
    meet_type = ExactPredicates.meet((ax, ay), (bx, by), (cx, cy), (dx, dy))
    # Lines meet at one point within open segments 
    meet_type == 1 && return cross
    # Lines don't meet at any points
    meet_type == -1 && return out
    # Lines meet at one or more points within closed segments
    if _isparallel(((ax, ay), (bx, by)), ((cx, cy), (dx, dy)))
        min_x, max_x = cx < dx ? (cx, dx) : (dx, cx)
        min_y, max_y = cy < dy ? (cy, dy) : (dy, cy)
        if (
            ((ax ≤ min_x && bx ≤ min_x) || (ax ≥ max_x && bx ≥ max_x)) &&
            ((ay ≤ min_y && by ≤ min_y) || (ay ≥ max_y && by ≥ max_y))
        )
            # a_point and b_point are on the same side of segment, don't overlap
            return hinge
        else
            return over
        end
    end
    # if lines aren't parallel then they must hinge
    return hinge
end

function _line_curve_process(
    line, curve;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    closed_line = false,
    closed_curve = false,
)
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # Determine curve endpoints
    nl = GI.npoint(line)
    nc = GI.npoint(curve)
    first_last_equal_line = equals(GI.getpoint(line, 1), GI.getpoint(line, nl))
    first_last_equal_curve = equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc))
    nl -= first_last_equal_line ? 1 : 0
    nc -= first_last_equal_curve ? 1 : 0
    closed_line |= first_last_equal_line
    closed_curve |= first_last_equal_curve
    
    # Loop over each line segment
    l_start = GI.getpoint(line, closed_line ? nl : 1)
    i = closed_line ? 1 : 2
    while i ≤ nl
        l_end = GI.getpoint(line, i)
        c_start = GI.getpoint(curve, closed_curve ? nc : 1)
        # Loop over each curve segment
        for j in (closed_curve ? 1 : 2):nc
            c_end = GI.getpoint(curve, j)
            # Check if line and curve segments meet
            seg_val = _segment_segment_orientation(
                (l_start, l_end),
                (c_start, c_end),
            )
            # if segments are touching
            if seg_val == line_over
                !in_allow && return false
                # at least one point in, meets requirments
                in_req_met = true
                if seg_val == line_over
                    point_val = point_segment_orientation(
                        l_start,
                        c_start, c_end,
                    )
                    if point_val != point_out
                        if point_segment_orientation(
                            l_end,
                            c_start, c_end,
                        ) != point_out
                            l_start = l_end
                            i += 1
                            break
                        elseif point_segment_orientation(
                            c_start,
                            l_start, l_end,
                        ) != point_out
                            l_start = c_start
                            break
                        elseif point_segment_orientation(
                            c_end,
                            l_start, l_end,
                        ) != point_out
                            l_start = c_end
                            break
                        end
                    end
                end
            else
                if seg_val == line_hinge
                    !on_allow && return false
                    # at least one point on, meets requirments
                    on_req_met = true
                elseif seg_val == line_cross
                    !in_allow && return false
                     # at least one point in, meets requirments
                     in_req_met = true
                end
                # no overlap for a give segment
                if j == nc
                    !out_allow && return false
                    out_req_met = true
                end
            end
            c_start = c_end
            j == nc && (i += 1)
        end
    end
    return in_req_met && on_req_met && out_req_met
end

function _line_filled_curve_interactions(
    line, curve;
    closed_line = false,
    filled_line = false,
)
    in_curve = false
    on_curve = false
    out_curve = false

    # Determine number of points in curve and line
    nl = GI.npoint(line)
    nc = GI.npoint(curve)
    first_last_equal_line = equals(GI.getpoint(line, 1), GI.getpoint(line, nl))
    first_last_equal_curve = equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc))
    nl -= first_last_equal_line ? 1 : 0
    nc -= first_last_equal_curve ? 1 : 0
    closed_line |= first_last_equal_line
    filled_line &= closed_line

    # See if first point is in an acceptable orientation
    l_start = GI.getpoint(line, closed_line ? nl : 1)
    point_val = point_filled_curve_orientation(l_start, curve)
    if point_val == point_in
        in_curve = true
    elseif point_val == point_on
        on_curve = true
    else  # point_val == point_out
        out_curve = true
    end

    # Check for any intersections between line and curve
    for i in (closed_line ? 1 : 2):nl
        l_end = GI.getpoint(line, i)
        c_start = GI.getpoint(curve, nc)
        # If already interacted with all regions of curve, can stop
        in_curve && on_curve && out_curve && break
        # Check next segment of line against curve
        for j in 1:nc
            c_end = GI.getpoint(curve, j)
            # Check if two line and curve segments meet
            seg_val = _segment_segment_orientation(
                (l_start, l_end),
                (c_start, c_end),
            )
            if seg_val != line_out
                # If line and curve meet, then at least one point is on boundary
                on_curve = true
                if seg_val == line_cross
                    # When crossing boundary, line is both in and out of curve
                    in_curve = true
                    out_curve = true
                else
                    if seg_val == line_over
                        sp = point_segment_orientation(l_start, c_start, c_end)
                        lp = point_segment_orientation(l_end, c_start, c_end)
                        if sp != point_in || lp != point_in
                            #=
                            Line crosses over segment endpoint, creating a hinge
                            with another segment.
                            =#
                            seg_val = line_hinge
                        end
                    end
                    if seg_val == line_hinge
                        #=
                        Can't determine all types of interactions (in, out) with
                        hinge as it could pass through multiple other segments
                        so calculate if segment endpoints and intersections are
                        in/out of filled curve
                        =#
                        ipoints = intersection_points(
                            GI.Line([l_start, l_end]),
                            curve
                        )
                        npoints = length(ipoints)  # since hinge, at least one
                        sort!(ipoints, by = p -> euclid_distance(p, l_start))
                        p_start = _tuple_point(l_start)
                        for i in 1:(npoints + 1)
                            p_end = i ≤ npoints ?
                                ipoints[i] :
                                _tuple_point(l_end)
                            mid_val = point_filled_curve_orientation(
                                (p_start .+ p_end) ./ 2,
                                curve,
                            )
                            if mid_val == point_in
                                in_curve = true
                            elseif mid_val == point_out
                                out_curve = true
                            end
                        end
                        # already checked segment against whole filled curve
                        l_start = l_end
                        break
                    end
                end
            end
            c_start = c_end
        end
        l_start = l_end
    end
    if filled_line && !in_curve
        if !out_curve  # line overlaps entire curve boundary
            in_curve = true  # line interior overlaps boundary filled interior
        else
            cent = centroid(line)
            if within(cent, curve)
                in_curve = true
                on_curve = true
            end
        end
    end

    return in_curve, on_curve, out_curve
end

function _line_polygon_process(
    line, polygon;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    closed_line = false,
    filled_line = false,
)
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # Check interaction of line with polygon's exterior boundary
    in_curve, on_curve, out_curve = _line_filled_curve_interactions(
        line, GI.getexterior(polygon);
        closed_line = closed_line,
        filled_line = filled_line,
    )
    if on_curve
        !on_allow && return false
        on_req_met = true
    end
    if out_curve
        !out_allow && return false
        out_req_met = true
    end
    !in_curve && return in_req_met && on_req_met && out_req_met

    # Loop over polygon holes
    for hole in GI.gethole(polygon)
        in_hole, on_hole, out_hole =_line_filled_curve_interactions(
            line, hole;
            closed_line = closed_line,
        )
        if in_hole
            !out_allow && return false
            out_req_met = true
        end
        if on_hole
            !on_allow && return false
            on_req_met = true
        end
        if !out_hole  # entire line is in/on hole, can't be in/on other holes
            in_curve = false
            break
        end
    end
    if in_curve
        !in_allow && return false
        in_req_met = true
    end
    return in_req_met && on_req_met && out_req_met
end

function _polygon_polygon_process(
    poly1, poly2;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
)
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require

    ext1 = GI.getexterior(poly1)
    e1_in_p2, e1_on_p2, e1_out_p2 = _line_polygon_process(
        ext1, poly2;
        in_allow = in_allow, in_require = in_require,
        on_allow = on_allow, on_require = on_require,
        out_allow = out_allow, out_require = out_require,
        closed_line = true,
        filled_line = true,
    )
    if e1_on_p2
        !on_allow && return false
        on_req_met = true
    end
    if e1_out_p2
        !out_allow && return false
        out_req_met = true
    end
    !e1_in_p2 && return in_req_met && on_req_met && out_req_met

    # is the part if p1 that is in p2 actually in p2?

    # does p1 touch any other edges or exit p2 at any point?

    # for h2 in GI.gethole(g2)
    #     # check if h2 is inside of e1
    #     e1_in_h2, e1_on_h2, e1_out_h2 = _line_filled_curve_interactions(
    #         ext1, h2;
    #         closed_line = true,
    #         filled_line = true,
    #     )
    #     # skip if poly1 doesn't interact with the hole at all
    #     !e1_in_h2 && !e1_on_h2 && break
    #     # if hole interacts with an edge of poly1
    #     if e1_on_h2
    #         !on_allow && return false
    #         on_req_met = true
    #         #=
    #         we know that h2 touches edge of p1 so:
    #         (1) no hole of p1 can touch the edge of p1 and
    #         (2) no other hole of p2 can line up with current h2
    #         This means there is at least a small border of p1 that is either
    #         inside of p2 (e1_out_h2) or outside of p2 (e1_in_h2)
    #         =#
    #         if e1_out_h2
    #             !in_allow && return false
    #             in_req_met = true
    #         end
    #         if e1_in_h2
    #             !out_allow && return false
    #             out_req_met = true
    #             # entirety of poly1 is within/on h2
    #             !e1_out_h2 && return in_req_met && on_req_met && out_req_met
    #         end
    #     else  # if hole is completly within poly1
    #         !in_allow && return false
    #         in_req_met = true
    #         # Check to see if h2 is within a hole of poly1
    #         for h1 in GI.gethole(poly1)
    #             h2_in_h1, h2_on_h1, h2_out_h1 = _line_filled_curve_interactions(
    #                 h2, h1;
    #                 closed_line = true,
    #                 filled_line = true,
    #             )
    #             if !h2_out_h1
    #                 !out_allow && return false
    #                 out_req_met = true
    #             else

    #             end
    #             # h2 is outside of h1 and cannot be excluded by another hole since it touches the boundary
    #             h2_on_h1 && h2_out_h1 && return false
    #             if !h2_out_h1  #h2 is within bounds of h1, so not in e1
    #                 h2_in_e1 = false
    #                 break
    #             end
    #         end

    #     end

    #     for h1 in GI.gethole(g1)
    #         _, h2_on_h1, h2_out_h1 = _line_filled_curve_interactions(
    #             h2, h1;
    #             closed_line = true,
    #         )
    #         # h2 is outside of h1 and cannot be excluded by another hole since it touches the boundary
    #         h2_on_h1 && h2_out_h1 && return false
    #         if !h2_out_h1  #h2 is within bounds of h1, so not in e1
    #             h2_in_e1 = false
    #             break
    #         end
    #     end
    #     h2_in_e1 && return false
    # end
    # return true
end


function _point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y2) = extent.X, extent.Y
    return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
end
