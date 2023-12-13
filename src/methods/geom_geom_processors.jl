@enum ProcessType within_process=1 disjoint_process=2 touch_process=3 coverby_process=4

@enum PointOrientation point_in=1 point_on=2 point_off=3

"""
    get_process_return_vals(
        process::ProcessType,
        exclude_boundaries = false,
    )::(Bool, Bool, Bool)

Returns a tuple of booleans which represent the boolean return value for it a
point is in, out, or on a given geometry. This is determined by the process
type as well as by the exclude_boundaries.

For within_process:
    if a point is in a geometry, we should return true
    if a point is out of a geomertry, we should return false
    if a point is on the boundary of a geometry, we should return false if we
        want to exclude boundaries, else we should return true
For disjoint_process:
    if a point is in a geometry, we should return false
    if a point is out of a geometry, we should return true
    if a point is on the boundary of a geometry, we should return true if we
        want to exclude boundaries, else we should return false
For touch_process:
    ???
"""
get_process_return_vals(process::ProcessType, exclude_boundaries = false) =
    (
        process == within_process,         # in value
        process == disjoint_process,       # out value
        (                                  # on value
            (process == touch_process) ||
            (process == within_process ?
                !exclude_boundaries : exclude_boundaries
            )
        )
    )

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
    process::ProcessType = within_process,
    exclude_boundaries = false,
    repeated_last_coord = false,
)
    n = GI.npoint(curve)
    first_last_equal = equals(GI.getpoint(curve, 1), GI.getpoint(curve, n))
    repeated_last_coord |= first_last_equal
    exclude_boundaries |= first_last_equal
    n -= first_last_equal ? 1 : 0
    in_val, out_val, on_val = get_process_return_vals(
        process,
        exclude_boundaries,
    )
    # Loop through all curve segments
    p_start = GI.getpoint(curve, repeated_last_coord ? n : 1)
    @inbounds for i in (repeated_last_coord ? 1 : 2):n
        p_end = GI.getpoint(curve, i)
        seg_val = _point_in_on_out_segment(point, p_start, p_end)
        seg_val == 1 && return in_val
        if seg_val == -1
            i == 2 && equals(point, p_start) && return on_val
            i == n && equals(point, p_end) && return on_val
            return in_val
        end
        p_start = p_end
    end
    return out_val
end

"""
    _point_in_on_out_segment(
        point::Point, start::Point, stop::Point;
        in::T = 1, on::T = -1, out::T = 0,
    )::T where {T}

Determines if a point is in, on, or out of a segment. If the point is 'on' the
segment it is on one of the segments endpoints. If it is 'in', it is on any
other point of the segment. If the point is not on any part of the segment, it
is 'out' of the segment.
"""
function _point_in_on_out_segment(
    point, start, stop;
    in::T = 1, on::T = -1, out::T = 0,
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


point_return_val(point_val, in_val, out_val, on_val) =
    point_val == 1 ?
        in_val :            # point is inside of polygon
        (point_val == 0 ?
            out_val :       # point is outside of polygon
            on_val          # point is on the edge of polygon
        )

"""
    _point_closed_curve_process(
        point, curve;
        process::ProcessType = within_process,
        exclude_boundary = false,
    )::Bool

Determines if a point meets the given process checks with respect to a closed
curve, which includes the space enclosed by the closed curve. Point should be an
object of Point trait and curve should be an object with a line string or linear
ring trait, that is assumed to be closed, regardless of repeated last point.

If checking within, then the point must be within the space enclosed by the
curve and if checking disjoint the point must not be outside of the curve.

Beyond specifying the process type, user can also specify if the geometry
boundaries should be included in the checks and if the curve should be closed
with repeated a repeated last coordinate matching the first coordinate.
"""
function _point_closed_curve_process(
    point, curve;
    process::ProcessType = within_process,
    exclude_boundaries = false,
)
    in_val, out_val, on_val = get_process_return_vals(
        process,
        exclude_boundaries,
    )
    return _point_in_on_out_closed_curve(
        point, curve;
        in = in_val, out = out_val, on = on_val
    )
end

"""
    _point_in_on_out_closed_curve(
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
function _point_in_on_out_closed_curve(
    point, curve;
    in::T = 1, on::T = -1, out::T = 0,
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
            f = u1 * v2 - u2 * v1
            if v2 > 0 && v1 ≤ 0                # Case 3, 9, 16, 21, 13, or 24
                f == 0 && return on            # Case 16 or 21
                f > 0 && (k += 1)              # Case 3 or 9
            elseif v1 > 0 && v2 ≤ 0            # Case 4, 10, 19, 20, 12, or 25
                f == 0 && return on            # Case 19 or 20
                f < 0 && (k += 1)              # Case 4 or 10
            elseif v2 == 0 && v1 < 0           # Case 7, 14, or 17
                f == 0 && return on            # Case 17
            elseif v1 == 0 && v2 < 0           # Case 8, 15, or 18
                f == 0 && return on            # Case 18
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
_point_polygon_process(
    point, polygon;
    process::ProcessType = within_process,
    exclude_boundaries = false,
) = _geom_polygon_process(
    point, polygon,
    _point_closed_curve_process;
    process = process,
    ext_exclude_boundaries = exclude_boundaries,
    hole_exclude_boundaries = !exclude_boundaries
)


function _line_curve_process(
    line, curve;
    process::ProcessType = within_process,
    exclude_boundaries = false,
    closed_line = false,
    closed_curve = false,
)
    nl = GI.npoint(line)
    nc = GI.npoint(curve)
    explicit_closed_line = equals(GI.getpoint(line, 1), GI.getpoint(line, nl))
    explicit_closed_curve = equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc))
    nl -= explicit_closed_line ? 1 : 0
    nc -= explicit_closed_curve ? 1 : 0
    closed_line |= explicit_closed_line
    closed_curve |= explicit_closed_curve
    exclude_boundaries |= explicit_closed_curve

    check_func = 
        if process == within_process
            (args...) -> _line_curve_within_checks(
                args...;
                nc = nc, closed_curve = closed_curve, exclude_boundaries,
            )
        else
            (args...) -> _line_curve_disjoint_checks(
                args...;
                nc = nc, closed_curve = closed_curve, exclude_boundaries,
            )
        end

    l_start = _tuple_point(GI.getpoint(line, closed_line ? nl : 1))
    i = closed_line ? 1 : 2
    while i ≤ nl
        l_end = _tuple_point(GI.getpoint(line, i))
        c_start = _tuple_point(GI.getpoint(curve, closed_curve ? nc : 1))
        for j in (closed_curve ? 1 : 2):nc
            c_end = _tuple_point(GI.getpoint(curve, j))
            
            meet_type = ExactPredicates.meet(l_start, l_end, c_start, c_end)
            passes_checks, break_loop, i, l_start = check_func(
                meet_type,
                l_start, l_end,
                c_start, c_end,
                i, j,
            )
            break_loop && break
            !passes_checks && return false
            c_start = c_end
        end
    end
    return true
end

function _line_curve_within_checks(
    meet_type,
    l_start, l_end,
    c_start, c_end,
    i, j;
    nc,
    closed_curve,
    exclude_boundaries,
)
    is_within = true
    break_loop = false
    #=
        if l_start is in/on curve and curve and line meet either at
        endpoints or are parallel and meet in multiple points
    =#
    if (
        meet_type == 0 &&
        _point_in_on_out_segment(l_start, c_start, c_end) != 0
    )
        #=
        if excluding first and last point of curve, make sure those points
        aren't within the line segment
        =#
        if exclude_boundaries && !closed_curve && ((
            j == 2 && _point_in_on_out_segment(c_start, l_start, l_end) != 0
        ) || (
            j == nc && _point_in_on_out_segment(c_end, l_start, l_end) != 0
        ))
            is_within = false      
        else  # if end points aren't being excluded 
            # if l_end is within curve, whole line is contained in curve
            if _point_in_on_out_segment(l_end, c_start, c_end) != 0
                i += 1
                l_start = l_end
                break_loop = true
            #=
            if c_start is in line, then need to find overlap for c_start
            to l_end as l_start to c_start is overlapping with curve
            =#
            elseif _point_in_on_out_segment(
                c_start,
                l_start, l_end,
            ) == 1
                l_start = c_start
                break_loop = true
            #=
            if c_end is in line, then need to find overlap for c_end to
            l_end as l_start to c_end is overlapping with curve
            =#
            elseif _point_in_on_out_segment(
                c_end,
                l_start, l_end,
            ) == 1
                l_start = c_end
                break_loop = true
            end
        end
    end
    #=
    if line segment has been checked against all curve segments and it isn't
    within any of them, line isn't within curve
    =#
    if j == nc
        is_within = false
    end
    return is_within, break_loop, i, l_start
end

function _line_curve_disjoint_checks(
    meet_type,
    l_start, l_end,
    c_start, c_end,
    i, j;
    nc,
    closed_curve,
    exclude_boundaries,
)
    is_disjoint = true
    break_loop = false
    #=
    if excluding first and last point of curve, line can still cross those
    points and be disjoint
    =#
    if (
        exclude_boundaries && meet_type == 0 &&
        !closed_curve && (j == 2 || j == nc)
    )
        # if line and curve are parallel, they cannot overlap and be disjoint
        if _isparallel(l_start, l_end, c_start, c_end)
            (p1, p2) =
                if j == 2 && equals(c_start, l_start)
                    (l_end, c_end)
                elseif j == 2 && equals(c_start, l_end)
                    (l_start, c_end)
                elseif j == nc && equals(c_end, l_start)
                    (l_end, c_start)
                elseif j == nc &&equals(c_end, l_end)
                    (l_start, c_start)
                else
                    is_disjoint = false
                end
            if is_disjoint && (
                _point_in_on_out_segment(p1, c_start, c_end) ||
                _point_in_on_out_segment(p2, l_start, l_end)
            )
                is_disjoint = false
            end
        #=
        if line and curve aren't parallel, they intersection must be either the
        start or end point of the curve to be disjoint
        =#
        else
            _, (_, c_frac) = _intersection_point(
                (l_start, l_end),
                (c_start, c_end),
            )
            if (
                j == 2 && c_frac != 0 ||
                j == nc && c_frac != 1
            )
                is_disjoint = false
            end
        end
    #=
    if not excluding first and last point of curve, line cannot intersect with
    any points of the curve
    =#
    elseif meet_type != -1
        is_disjoint = false
    end
    #=
    if line segment has been checked against all curve segments and is disjoint
    from all of them, we can now check the next line segment
    =#
    if j == nc
        i += 1
        l_start = l_end
    end
    return is_disjoint, break_loop, i, l_start
end

function _line_closed_curve_process(
    line, curve;
    process::ProcessType = within_process,
    exclude_boundaries = false,
    close = false,
    line_is_poly_ring = false,
)
    #=
    if line isn't the external ring of a polygon, see if at least one point is
    within the closed curve - else, ring is "filled in" and has points within
    closed curve
    =# 
    point_in = line_is_poly_ring
    in_val, out_val, on_val = get_process_return_vals(process, exclude_boundaries)
    # Determine number of points in curve and line
    nc = GI.npoint(curve)
    nc -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc)) ? 1 : 0
    nl = GI.npoint(line)
    nl -= equals(GI.getpoint(line, 1), GI.getpoint(line, nl)) ? 1 : 0

    # Check to see if first point in line is within curve
    l_start = _tuple_point(GI.getpoint(line, close ? nl : 1))
    point_val = _point_in_on_out_closed_curve(l_start, curve)
    point_in |= point_val == 1  # check if point is within closed curve
    point_return = point_return_val(point_val, in_val, out_val, on_val)

    # point is not in correct orientation to curve given process and boundary
    !point_return && return point_return

    # Check for any intersections between line and curve
    for i in (close ? 1 : 2):nl
        l_end = _tuple_point(GI.getpoint(line, i))
        c_start = _tuple_point(GI.getpoint(curve, nc))
        for j in 1:nc
            c_end = _tuple_point(GI.getpoint(curve, j))
            # Check if edges intersect --> line crosses --> wrong orientation
            meet_type = ExactPredicates.meet(l_start, l_end, c_start, c_end)
            # open line segments meet in a single point
            meet_type == 1 && return false
            #=
            closed line segments meet in one or several points -> meet at a
            vertex or on the edge itself (parallel)
            =#
            if meet_type == 0
                # See if segment is parallel and within curve edge
                p1_on_seg = _point_in_on_out_segment(
                    l_start,
                    c_start, c_end,
                ) != 0
                exclude_boundaries && p1_on_seg && return false
                p2_on_seg = _point_in_on_out_segment(
                    l_end,
                    c_start, c_end,
                ) != 0
                exclude_boundaries && p2_on_seg && return false
                # if segment isn't contained within curve edge
                if !p1_on_seg || !p2_on_seg 
                    # Make sure l_start is in corrent orientation
                    if !p1_on_seg
                        p1_val = _point_in_on_out_closed_curve(l_start, curve)
                        # check if point is within closed curve
                        point_in |= p1_on_seg == 1
                        !point_return_val(p1_val, in_val, out_val, on_val) &&
                            return false
                    end
                    # Make sure l_end is in is in corrent orientation
                    if !p2_on_seg
                        p2_val = _point_in_on_out_closed_curve(l_end, curve)
                        # check if point is within closed curve
                        point_in |= p2_val == 1
                        !point_return_val(p2_val, in_val, out_val, on_val) &&
                            return false
                    end
                    #=
                    If both endpoints are in the correct orientation, but not
                    parallel to the edge, make sure that midpoints between the
                    intersections along the segment are also in the correct
                    orientation
                    =# 
                    mid_vals, mid_in = _segment_mids_closed_curve_process(
                        l_start, l_end, curve;
                        process = process,
                        exclude_boundaries = exclude_boundaries,
                    )
                    point_in |= mid_in
                    # midpoint on the wrong side of the curve
                    !mid_vals && return false
                    # line segment is fully within or on curve 
                    break 
                end
            end
            c_start = c_end
        end
        l_start = l_end
    end
    # check if line is on any curve edges or vertcies
    return process == within_process ? point_in : true
end

function _segment_mids_closed_curve_process(
    l_start, l_end, curve;
    process::ProcessType = within_process,
    exclude_boundaries = false,
)
    point_in = false
    in_val, out_val, on_val = get_process_return_vals(process, exclude_boundaries)
    # Find intersection points
    ipoints = intersection_points(
        GI.Line([l_start, l_end]),
        curve
    )
    npoints = length(ipoints)
    if npoints < 3  # only intersection points are the endpoints
        mid_val = _point_in_on_out_closed_curve((l_start .+ l_end) ./ 2, curve)
        point_in |= mid_val == 1
        mid_return = point_return_val(mid_val, in_val, out_val, on_val)
        !mid_return && return (false, point_in)
    else  # more intersection points than the endpoints
        # sort intersection points along the line
        sort!(ipoints, by = p -> euclid_distance(p, l_start))
        p_start = ipoints[1]
        for i in 2:npoints
            p_end = ipoints[i]
            # check if midpoint of intersection points is within the curve
            mid_val = _point_in_on_out_closed_curve(
                (p_start .+ p_end) ./ 2,
                curve,
            )
            point_in |= mid_val == 1
            mid_return = point_return_val(mid_val, in_val, out_val, on_val)
            !mid_return && return (false, point_in)
        end
    end
    # all intersection point midpoints were in or on the curve
    return true, point_in
end

_line_polygon_process(
    line, polygon;
    process::ProcessType = within_process,
    exclude_boundaries = true,
    close = false,
    line_is_poly_ring = false,
) = _geom_polygon_process(
    line, polygon,
    (args...; kwargs...) -> _line_closed_curve_process(
        args...;
        kwargs...,
        close = close,
        line_is_poly_ring = line_is_poly_ring,
    );
    process = process,
    ext_exclude_boundaries = exclude_boundaries,
    hole_exclude_boundaries =
        process == within_process ?
            exclude_boundaries :
            !exclude_boundaries,
)


function _geom_polygon_process(
    geom, polygon, geom_closed_curve_func;
    process::ProcessType = within_process,
    ext_exclude_boundaries = true,
    hole_exclude_boundaries = false
)
    # Check interaction of geom with polygon's exterior boundary
    ext_val = geom_closed_curve_func(
        geom, GI.getexterior(polygon);
        process = process, exclude_boundaries = ext_exclude_boundaries,
    )
    
    #=
    If checking within and geom is outside of exterior ring, return false or
    if checking disjoint and geom is outside of exterior ring, return true.
    =#
    ((process == within_process && !ext_val) ||
        (process == disjoint_process && ext_val)
    ) && return ext_val
    
    # If geom is within the polygon, need to check interactions with holes
    for hole in GI.gethole(polygon)
        hole_val = geom_closed_curve_func(
            geom, hole,
            process = (
                process == within_process ?
                    disjoint_process :
                    within_process
            ),
            exclude_boundaries = hole_exclude_boundaries
        )
        #=
        If checking within and geom is not disjoint from hole, return false or
        if checking disjoint and geom is within hole, return true.
        =#
        process == within_process && !hole_val && return false
        process == disjoint_process && hole_val && return true
    end
    return ext_val
end

function _point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y2) = extent.X, extent.Y
    return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
end
