# # Line-curve interaction

#= Code is based off of DE-9IM Standards (https://en.wikipedia.org/wiki/DE-9IM)
and attempts a standardized solution for most of the functions.
=#

"""
    Enum PointOrientation

Enum for the orientation of a point with respect to a curve. A point can be
`point_in` the curve, `point_on` the curve, or `point_out` of the curve.
"""
@enum PointOrientation point_in=1 point_on=2 point_out=3

#=
Determines if a point meets the given checks with respect to a curve.

If in_allow is true, the point can be on the curve interior.
If on_allow is true, the point can be on the curve boundary.
If out_allow is true, the point can be disjoint from the curve.

If the point is in an "allowed" location, return true. Else, return false.

If closed_curve is true, curve is treated as a closed curve where the first and
last point are connected by a segment.
=#
function _point_curve_process(
    point, curve;
    in_allow, on_allow, out_allow,
    closed_curve = false,
)
    # Determine if curve is closed
    n = GI.npoint(curve)
    first_last_equal = equals(GI.getpoint(curve, 1), GI.getpoint(curve, n))
    closed_curve |= first_last_equal
    n -= first_last_equal ? 1 : 0
    # Loop through all curve segments
    p_start = GI.getpoint(curve, closed_curve ? n : 1)
    @inbounds for i in (closed_curve ? 1 : 2):n
        p_end = GI.getpoint(curve, i)
        seg_val = _point_segment_orientation(point, p_start, p_end)
        seg_val == point_in && return in_allow
        if seg_val == point_on
            if !closed_curve  # if point is on curve endpoints, it is "on"
                i == 2 && equals(point, p_start) && return on_allow
                i == n && equals(point, p_end) && return on_allow
            end
            return in_allow
        end
        p_start = p_end
    end
    return out_allow
end

#=
Determines if a point meets the given checks with respect to a polygon.

If in_allow is true, the point can be within the polygon interior
If on_allow is true, the point can be on the polygon boundary.
If out_allow is true, the point can be disjoint from the polygon.

If the point is in an "allowed" location, return true. Else, return false.
=#
function _point_polygon_process(
    point, polygon;
    in_allow, on_allow, out_allow,
)
    # Check interaction of geom with polygon's exterior boundary
    ext_val = _point_filled_curve_orientation(point, GI.getexterior(polygon))
    # If a point is outside, it isn't interacting with any holes
    ext_val == point_out && return out_allow
    # if a point is on an external boundary, it isn't interacting with any holes
    ext_val == point_on && return on_allow
    
    # If geom is within the polygon, need to check interactions with holes
    for hole in GI.gethole(polygon)
        hole_val = _point_filled_curve_orientation(point, hole)
        # If a point in in a hole, it is outside of the polygon
        hole_val == point_in && return out_allow
        # If a point in on a hole edge, it is on the edge of the polygon
        hole_val == point_on && return on_allow
    end
    
    # Point is within external boundary and on in/on any holes
    return in_allow
end

#=
Determines if a line meets the given checks with respect to a curve.

If over_allow is true, segments of the line and curve can be co-linear.
If cross_allow is true, segments of the line and curve can cross.
If on_allow is true, endpoints of either the line or curve can intersect a 
    segment of the other geometry.
If cross_allow is true, segments of the line and curve can be disjoint.

If in_require is true, the interiors of the line and curve must meet in at least
    one point.
If on_require is true, the bounday of one of the two geometries can meet the
    interior or boundary of the other geometry in at least one point.
If out_require is true, there must be at least one point of the given line that
    is exterior of the curve.

If the point is in an "allowed" location and meets all requirments, return true.
Else, return false.

If closed_line is true, line is treated as a closed line where the first and
last point are connected by a segment. Same with closed_curve.
=#
@inline function _line_curve_process(line, curve; 
    over_allow, cross_allow, kw...
)
    skip, returnval = _maybe_skip_disjoint_extents(line, curve;
        in_allow=(over_allow | cross_allow), kw...
    )
    if skip 
        return returnval
    else
        return _inner_line_curve_process(line, curve; over_allow, cross_allow, kw...)
    end
end

function _inner_line_curve_process(
    line, curve;
    over_allow, cross_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    closed_line = false,
    closed_curve = false,
)
    # Set up requirments
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
    l_start = _tuple_point(GI.getpoint(line, closed_line ? nl : 1))
    i = closed_line ? 1 : 2
    while i ≤ nl
        l_end = _tuple_point(GI.getpoint(line, i))
        c_start = _tuple_point(GI.getpoint(curve, closed_curve ? nc : 1))
        # Loop over each curve segment
        for j in (closed_curve ? 1 : 2):nc
            c_end = _tuple_point(GI.getpoint(curve, j))
            # Check if line and curve segments meet
            seg_val, intr1, _ = _intersection_point(Float64, (l_start, l_end), (c_start, c_end))
            # If segments are co-linear
            if seg_val == line_over
                !over_allow && return false
                # at least one point in, meets requirments
                in_req_met = true
                point_val = _point_segment_orientation(l_start, c_start, c_end)
                # If entire segment isn't covered, consider remaining section
                if point_val != point_out
                    i, l_start, break_off = _find_new_seg(i, l_start, l_end, c_start, c_end)
                    break_off && break
                end
            else
                if seg_val == line_cross
                    !cross_allow && return false
                    in_req_met = true
                elseif seg_val == line_hinge  # could cross or overlap
                    # Determine location of intersection point on each segment
                    (_, (α, β)) = intr1
                    if ( # Don't consider edges of curves as they can't cross
                        (!closed_line && ((α == 0 && i == 2) || (α == 1 && i == nl))) ||
                        (!closed_curve && ((β == 0 && j == 2) || (β == 1 && j == nc)))
                    )
                        !on_allow && return false
                        on_req_met = true
                    else
                        in_req_met = true
                        # If needed, determine if hinge actually crosses
                        if (!cross_allow || !over_allow) && α != 0 && β != 0
                            # Find next pieces of hinge to see if line and curve cross
                            l, c = _find_hinge_next_segments(
                                α, β, l_start, l_end, c_start, c_end,
                                i, line, j, curve,
                            )
                            next_val, _, _ = _intersection_point(Float64, l, c)
                            if next_val == line_hinge
                                !cross_allow && return false
                            else
                                !over_allow && return false
                            end
                        end
                    end
                end
                # no overlap for a give segment, some of segment must be out of curve
                if j == nc
                    !out_allow && return false
                    out_req_met = true
                end
            end
            c_start = c_end  # consider next segment of curve
            if j == nc  # move on to next line segment
                i += 1
                l_start = l_end
            end
        end
    end
    return in_req_met && on_req_met && out_req_met
end

#= If entire segment (le to ls) isn't covered by segment (cs to ce), find remaining section
part of section outside of cs to ce. If completly covered, increase segment index i. =#
function _find_new_seg(i, ls, le, cs, ce)
    break_off = true
    if _point_segment_orientation(le, cs, ce) != point_out
        ls = le
        i += 1
    elseif !equals(ls, cs) && _point_segment_orientation(cs, ls, le) != point_out
        ls = cs
    elseif !equals(ls, ce) && _point_segment_orientation(ce, ls, le) != point_out
        ls = ce
    else
        break_off = false
    end
    return i, ls, break_off
end

#= Find next set of segments needed to determine if given hinge segments cross or not.=#
function _find_hinge_next_segments(α, β, ls, le, cs, ce, i, line, j, curve) 
    next_seg = if β == 1
        if α == 1  # hinge at endpoints, so next segment of both is needed
            ((le, _tuple_point(GI.getpoint(line, i + 1))), (ce, _tuple_point(GI.getpoint(curve, j + 1))))
        else  # hinge at curve endpoint and line interior point, curve next segment needed 
            ((ls, le), (ce, _tuple_point(GI.getpoint(curve, j + 1))))
        end
    else  # hinge at curve interior point and line endpoint, line next segment needed
        ((le, _tuple_point(GI.getpoint(line, i + 1))), (cs, ce))
    end
    return next_seg
end
#=
Determines if a line meets the given checks with respect to a polygon.

If in_allow is true, segments of the line can be in the polygon interior.
If on_allow is true, segments of the line can be on the polygon's boundary.
If out_allow is true, segments of the line can be outside of the polygon.

If in_require is true, the interiors of the line and polygon must meet in at
    least one point.
If on_require is true, the line must have at least one point on the polygon'same
    boundary.
If out_require is true, the line must have at least one point outside of the
    polygon.

If the point is in an "allowed" location and meets all requirments, return true.
Else, return false.

If closed_line is true, line is treated as a closed line where the first and
last point are connected by a segment.
=#
@inline function _line_polygon_process(line, polygon; kw...)
    skip, returnval = _maybe_skip_disjoint_extents(line, polygon; kw...)
    if skip 
        return returnval
    else
        return _inner_line_polygon_process(line, polygon; kw...)
    end
end

function _inner_line_polygon_process(
    line, polygon;
    closed_line=false,
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
)

    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # Check interaction of line with polygon's exterior boundary
    in_curve, on_curve, out_curve = _line_filled_curve_interactions(
        line, GI.getexterior(polygon);
        closed_line = closed_line,
    )
    if on_curve
        !on_allow && return false
        on_req_met = true
    end
    if out_curve
        !out_allow && return false
        out_req_met = true
    end
    # If no points within the polygon, the line is disjoint and we are done
    !in_curve && return in_req_met && on_req_met && out_req_met

    # Loop over polygon holes
    for hole in GI.gethole(polygon)
        in_hole, on_hole, out_hole =_line_filled_curve_interactions(
            line, hole;
            closed_line = closed_line,
        )
        if in_hole  # line in hole is equivalent to being out of polygon
            !out_allow && return false
            out_req_met = true
        end
        if on_hole  # hole bounday is polygon boundary
            !on_allow && return false
            on_req_met = true
        end
        if !out_hole  # entire line is in/on hole, can't be in/on other holes
            in_curve = false
            break
        end
    end
    if in_curve  # entirely of curve isn't within a hole
        !in_allow && return false
        in_req_met = true
    end
    return in_req_met && on_req_met && out_req_met
end

#=
Determines if a polygon meets the given checks with respect to a polygon.

If in_allow is true, the polygon's interiors must intersect.
If on_allow is true, the one of the polygon's boundaries must either interact
    with the other polygon's boundary or interior.
If out_allow is true, the first polygon must have interior regions outside of
    the second polygon.

If in_require is true, the polygon interiors must meet in at least one point.
If on_require is true, one of the polygon's must have at least one boundary
    point in or on the other polygon.
If out_require is true, the first polygon must have at least one interior point
    outside of the second polygon.

If the point is in an "allowed" location and meets all requirments, return true.
Else, return false.
=#
@inline function _polygon_polygon_process(poly1, poly2; kw...)
    skip, returnval = _maybe_skip_disjoint_extents(poly1, poly2; kw...)
    if skip 
        return returnval
    else
        return _polygon_polygon_process(poly1, poly2; kw...)
    end
end

function _inner_polygon_polygon_process(
    poly1, poly2;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
)
    skip, returnval = _maybe_skip_disjoint_extents(poly1, poly2;
        in_allow, on_allow, out_allow, 
        in_require, on_require, out_require,
    )
    skip && return returnval
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # Check if exterior of poly1 is within poly2
    ext1 = GI.getexterior(poly1)
    ext2 = GI.getexterior(poly2)
    # Check if exterior of poly1 is in polygon 2
    e1_in_p2, e1_on_p2, e1_out_p2 = _line_polygon_interactions(
        ext1, poly2;
        closed_line = true,
    )
    if e1_on_p2
        !on_allow && return false
        on_req_met = true
    end
    if e1_out_p2
        !out_allow && return false
        out_req_met = true
    end

    if !e1_in_p2
        # if exterior ring isn't in poly2, check if it surrounds poly2
        _, _, e2_out_e1 = _line_filled_curve_interactions(
            ext2, ext1;
            closed_line = true,
        )  # if they really are disjoint, we are done
        e2_out_e1 && return in_req_met && on_req_met && out_req_met
    end
    # If interiors interact, check if poly2 interacts with any of poly1's holes
    for h1 in GI.gethole(poly1)
        h1_in_p2, h1_on_p2, h1_out_p2 = _line_polygon_interactions(
            h1, poly2;
            closed_line = true,
        )
        if h1_on_p2
            !on_allow && return false
            on_req_met = true
        end
        if h1_out_p2
            !out_allow && return false
            out_req_met = true
        end
        if !h1_in_p2
            # If hole isn't in poly2, see if poly2 is in hole
            _, _, e2_out_h1 = _line_filled_curve_interactions(
                ext2, h1;
                closed_line = true,
            )
            # hole encompasses all of poly2
            !e2_out_h1 && return in_req_met && on_req_met && out_req_met
            break
        end
    end
    #=
    poly2 isn't outside of poly1 and isn't in a hole, poly1 interior must
    interact with poly2 interior
    =#
    !in_allow && return false
    in_req_met = true

    # If any of poly2 holes are within poly1, part of poly1 is exterior to poly2
    for h2 in GI.gethole(poly2)
        h2_in_p1, h2_on_p1, _ = _line_polygon_interactions(
            h2, poly1;
            closed_line = true,
        )
        if h2_on_p1
            !on_allow && return false
            on_req_met = true
        end
        if h2_in_p1
            !out_allow && return false
            out_req_met = true
        end
    end
    return in_req_met && on_req_met && out_req_met 
end

#=
Determines if a point is in, on, or out of a segment. If the point is `on` the
segment it is on one of the segments endpoints. If it is `in`, it is on any
other point of the segment. If the point is not on any part of the segment, it
is `out` of the segment.

Point should be an object of point trait and curve should be an object with a
linestring or linearring trait.

Can provide values of in, on, and out keywords, which determines return values
for each scenario. 
=#
function _point_segment_orientation(
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

#=
Determine if point is in, on, or out of a closed curve, which includes the space
enclosed by the closed curve.

`In` means the point is within the closed curve (excluding edges and vertices).
`On` means the point is on an edge or a vertex of the closed curve.
`Out` means the point is outside of the closed curve.

Point should be an object of point trait and curve should be an object with a
linestring or linearring trait, that is assumed to be closed, regardless of
repeated last point.

Can provide values of in, on, and out keywords, which determines return values
for each scenario. 

Note that this uses the Algorithm by Hao and Sun (2018):
https://doi.org/10.3390/sym10100477
Paper seperates orientation of point and edge into 26 cases. For each case, it
is either a case where the point is on the edge (returns on), where a ray from
the point (x, y) to infinity along the line y = y cut through the edge (k += 1),
or the ray does not pass through the edge (do nothing and continue). If the ray
passes through an odd number of edges, it is within the curve, else outside of
of the curve if it didn't return 'on'.
See paper for more information on cases denoted in comments.
=#
function _point_filled_curve_orientation(
    point, curve;
    in::T = point_in, on::T = point_on, out::T = point_out,
) where {T}
    x, y = GI.x(point), GI.y(point)
    n = GI.npoint(curve)
    n -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, n)) ? 1 : 0
    k = 0  # counter for ray crossings
    p_start = GI.getpoint(curve, n)
    for (i, p_end) in enumerate(GI.getpoint(curve))
        i > n && break
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

#=
Determines the types of interactions of a line with a filled-in curve. By
filled-in curve, I am referring to the exterior ring of a poylgon, for example.

Returns a tuple of booleans: (in_curve, on_curve, out_curve).

If in_curve is true, some of the lines interior points interact with the curve's
    interior points.
If on_curve is true, endpoints of either the line intersect with the curve or
    the line interacts with the polygon boundary.
If out_curve is true, at least one segments of the line is outside the curve.

If closed_line is true, line is treated as a closed line where the first and
last point are connected by a segment.
=#
function _line_filled_curve_interactions(
    line, curve;
    closed_line = false,
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

    # See if first point is in an acceptable orientation
    l_start = _tuple_point(GI.getpoint(line, closed_line ? nl : 1))
    point_val = _point_filled_curve_orientation(l_start, curve)
    if point_val == point_in
        in_curve = true
    elseif point_val == point_on
        on_curve = true
    else  # point_val == point_out
        out_curve = true
    end

    # Check for any intersections between line and curve
    for i in (closed_line ? 1 : 2):nl
        l_end = _tuple_point(GI.getpoint(line, i))
        c_start = _tuple_point(GI.getpoint(curve, nc))
        # If already interacted with all regions of curve, can stop
        in_curve && on_curve && out_curve && break
        # Check next segment of line against curve
        for j in 1:nc
            c_end = _tuple_point(GI.getpoint(curve, j))
            # Check if two line and curve segments meet
            seg_val, _, _ = _intersection_point(Float64, (l_start, l_end), (c_start, c_end))
            if seg_val != line_out
                # If line and curve meet, then at least one point is on boundary
                on_curve = true
                if seg_val == line_cross
                    # When crossing boundary, line is both in and out of curve
                    in_curve = true
                    out_curve = true
                else
                    if seg_val == line_over
                        sp = _point_segment_orientation(l_start, c_start, c_end)
                        lp = _point_segment_orientation(l_end, c_start, c_end)
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
                        dist_from_lstart = let l_start = l_start
                            x -> _euclid_distance(Float64, x, l_start)
                        end
                        sort!(ipoints, by = dist_from_lstart)
                        p_start = _tuple_point(l_start)
                        for i in 1:(npoints + 1)
                            p_end = i ≤ npoints ?
                                _tuple_point(ipoints[i]) :
                                l_end
                            mid_val = _point_filled_curve_orientation(
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
    return in_curve, on_curve, out_curve
end

#=
Determines the types of interactions of a line with a polygon. 

Returns a tuple of booleans: (in_poly, on_poly, out_poly).

If in_poly is true, some of the lines interior points interact with the polygon
    interior points.
If in_poly is true, endpoints of either the line intersect with the polygon or
    the line interacts with the polygon boundary, including hole bounaries.
If out_curve is true, at least one segments of the line is outside the polygon,
    including inside of holes.

If closed_line is true, line is treated as a closed line where the first and
last point are connected by a segment.
=#
function _line_polygon_interactions(
    line, polygon;
    closed_line = false,
)

    in_poly, on_poly, out_poly = _line_filled_curve_interactions(
        line, GI.getexterior(polygon);
        closed_line = closed_line,
    )
    !in_poly && return (in_poly, on_poly, out_poly)
    # Loop over polygon holes
    for hole in GI.gethole(polygon)
        in_hole, on_hole, out_hole =_line_filled_curve_interactions(
            line, hole;
            closed_line = closed_line,
        )
        if in_hole
            out_poly = true
        end
        if on_hole
            on_poly = true
        end
        if !out_hole  # entire line is in/on hole, can't be in/on other holes
            in_poly = false
            return (in_poly, on_poly, out_poly)
        end
    end
    return in_poly, on_poly, out_poly
end

function _point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y2) = extent.X, extent.Y
    return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
end

# Disjoint extent optimisation: skip work based on geom extent intersection
# returns Tuple{Bool, Bool} for (skip, returnval)
@inline function _maybe_skip_disjoint_extents(a, b;
    in_allow, on_allow, out_allow, 
    in_require, on_require, out_require,
    kw...
)
    if (in_allow || in_require || on_allow || on_require)
        # If we need line or interior and no exterior
        if !(out_require || out_allow) && Extents.disjoint(GI.extent(a), GI.extent(b))
            # Return false for disjoint geometries
            return true, false
        end
    else
        # If we need no line or interior, but need exterior
        if (out_require || out_allow) && Extents.disjoint(GI.extent(a), GI.extent(b))
            # Return true for disjoint geometries
            return true, true
        end
    end
    return false, false
end
