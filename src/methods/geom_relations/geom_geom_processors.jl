# Preserve spherical Cartesian points in shared relation processors. Converting them
# to (x,y) tuples silently reinterprets the coordinates as longitude and latitude.
_processor_point(::Planar, p) = _tuple_point(p)
_processor_point(::Spherical, p) = _spherical_kernel_point(p)

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
Return whether the point occupies a permitted curve location: `in_allow` for the interior,
`on_allow` for the boundary, or `out_allow` for the exterior.

`closed_curve` connects the first and last vertices.
=#
function _point_curve_process(
    m::Manifold, point, curve;
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
        seg_val = _point_segment_orientation(m, point, p_start, p_end)
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
Return whether the point occupies a permitted polygon location: `in_allow` for the interior,
`on_allow` for the boundary, or `out_allow` for the exterior.
=#
function _point_polygon_process(
    m::Manifold, point, polygon;
    in_allow, on_allow, out_allow, exact,
)
    skip, returnval = _maybe_skip_disjoint_extents(m, point, polygon; in_allow, on_allow, out_allow, on_require = false, out_require = false, in_require = false)
    skip && return returnval
    # Check interaction of geom with polygon's exterior boundary
    ext_val = _point_filled_curve_orientation(m, point, GI.getexterior(polygon); exact)
    # If a point is outside, it isn't interacting with any holes
    ext_val == point_out && return out_allow
    # if a point is on an external boundary, it isn't interacting with any holes
    ext_val == point_on && return on_allow
    
    # If geom is within the polygon, need to check interactions with holes
    for hole in GI.gethole(polygon)
        hole_val = _point_filled_curve_orientation(m, point, hole; exact)
        # If a point in in a hole, it is outside of the polygon
        hole_val == point_in && return out_allow
        # If a point in on a hole edge, it is on the edge of the polygon
        hole_val == point_on && return on_allow
    end
    
    # Point is within external boundary and on in/on any holes
    return in_allow
end

#=
Return whether line/curve interactions satisfy all allowed and required flags.

`over_allow` permits collinear segments; `cross_allow` permits crossings. `on_allow` permits
endpoint contacts, and `out_allow` permits disjoint segments.

`in_require` requires interior contact; `on_require` requires boundary contact with the other
geometry. `out_require` requires a line point outside the curve.

`closed_line` and `closed_curve` connect the corresponding first and last vertices.
=#
@inline function _line_curve_process(m::Manifold, line, curve; 
    over_allow, cross_allow, kw...
)
    skip, returnval = _maybe_skip_disjoint_extents(m, line, curve;
        in_allow=(over_allow | cross_allow), kw...
    )
    skip && return returnval

    return _inner_line_curve_process(m, line, curve; over_allow, cross_allow, kw...)
end

function _inner_line_curve_process(
    m::Manifold, line, curve;
    over_allow, cross_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    closed_line = false, closed_curve = false,
    exact,
)
    # Set up requirements
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
    l_start = _processor_point(m, GI.getpoint(line, closed_line ? nl : 1))
    i = closed_line ? 1 : 2
    while i ≤ nl
        l_end = _processor_point(m, GI.getpoint(line, i))
        c_start = _processor_point(m, GI.getpoint(curve, closed_curve ? nc : 1))
        # Loop over each curve segment
        for j in (closed_curve ? 1 : 2):nc
            c_end = _processor_point(m, GI.getpoint(curve, j))
            # Check if line and curve segments meet
            seg_val, α, β = _seg_seg_orientation(m, l_start, l_end, c_start, c_end; exact)
            # If segments are co-linear
            if seg_val == line_over
                !over_allow && return false
                # at least one point in, meets requirements
                in_req_met = true
                point_val = _point_segment_orientation(m, l_start, c_start, c_end)
                # If entire segment isn't covered, consider remaining section
                if point_val != point_out
                    i, l_start, break_off = _find_new_seg(m, i, l_start, l_end, c_start, c_end)
                    break_off && break
                end
            else
                if seg_val == line_cross
                    !cross_allow && return false
                    in_req_met = true
                elseif seg_val == line_hinge  # could cross or overlap
                    # `α`/`β` locate the intersection along each segment
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
                            l, c = _find_hinge_next_segments(m,
                                α, β, l_start, l_end, c_start, c_end,
                                i, line, j, curve,
                            )
                            next_val, _, _ = _seg_seg_orientation(m, l[1], l[2], c[1], c[2]; exact)
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
part of section outside of cs to ce. If completely covered, increase segment index i. =#
function _find_new_seg(m::Manifold, i, ls, le, cs, ce)
    break_off = true
    if _point_segment_orientation(m, le, cs, ce) != point_out
        ls = le
        i += 1
    elseif !equals(ls, cs) && _point_segment_orientation(m, cs, ls, le) != point_out
        ls = cs
    elseif !equals(ls, ce) && _point_segment_orientation(m, ce, ls, le) != point_out
        ls = ce
    else
        break_off = false
    end
    return i, ls, break_off
end

#= Find next set of segments needed to determine if given hinge segments cross or not.=#
function _find_hinge_next_segments(m::Manifold, α, β, ls, le, cs, ce, i, line, j, curve)
    next_seg = if β == 1
        if α == 1  # hinge at endpoints, so next segment of both is needed
            ((le, _processor_point(m, GI.getpoint(line, i + 1))), (ce, _processor_point(m, GI.getpoint(curve, j + 1))))
        else  # hinge at curve endpoint and line interior point, curve next segment needed 
            ((ls, le), (ce, _processor_point(m, GI.getpoint(curve, j + 1))))
        end
    else  # hinge at curve interior point and line endpoint, line next segment needed
        ((le, _processor_point(m, GI.getpoint(line, i + 1))), (cs, ce))
    end
    return next_seg
end
#=
Return whether line/polygon interactions satisfy all allowed and required flags.

`in_allow`, `on_allow`, and `out_allow` permit line portions in the polygon interior,
boundary, and exterior. The corresponding `_require` flags require at least one point there.

`closed_line` connects the first and last line vertices.
=#
@inline function _line_polygon_process(m::Manifold, line, polygon; kw...)
    skip, returnval = _maybe_skip_disjoint_extents(m, line, polygon; kw...)
    skip && return returnval
    return _inner_line_polygon_process(m, line, polygon; kw...)
end

function _inner_line_polygon_process(
    m::Manifold, line, polygon;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    exact, closed_line = false,
)
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # Check interaction of line with polygon's exterior boundary
    in_curve, on_curve, out_curve = _line_filled_curve_interactions(
        m, line, GI.getexterior(polygon);
        exact, closed_line = closed_line,
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
            m, line, hole;
            exact, closed_line = closed_line,
        )
        if in_hole  # line in hole is equivalent to being out of polygon
            !out_allow && return false
            out_req_met = true
        end
        if on_hole  # hole boundary is polygon boundary
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
Return whether polygon interactions satisfy all allowed and required flags.

`in_allow` permits interior overlap; `on_allow` permits boundary contact with the other
polygon. `out_allow` permits first-polygon interior outside the second.

Each corresponding `_require` flag requires at least one point with that relation.
=#
@inline function _polygon_polygon_process(m::Manifold, poly1, poly2; kw...)
    skip, returnval = _maybe_skip_disjoint_extents(m, poly1, poly2; kw...)
    skip && return returnval
    return _inner_polygon_polygon_process(m, poly1, poly2; kw...)
end

function _inner_polygon_polygon_process(
    m::Manifold, poly1, poly2;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    exact,
)
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # Check if exterior of poly1 is within poly2
    ext1 = GI.getexterior(poly1)
    ext2 = GI.getexterior(poly2)
    # Check if exterior of poly1 is in polygon 2
    e1_in_p2, e1_on_p2, e1_out_p2 = _line_polygon_interactions(
        m, ext1, poly2;
        exact, closed_line = true,
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
            m, ext2, ext1;
            exact, closed_line = true,
        )  # if they really are disjoint, we are done
        e2_out_e1 && return in_req_met && on_req_met && out_req_met
    end
    # If interiors interact, check if poly2 interacts with any of poly1's holes
    for h1 in GI.gethole(poly1)
        h1_in_p2, h1_on_p2, h1_out_p2 = _line_polygon_interactions(
            m, h1, poly2;
            exact, closed_line = true,
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
                m, ext2, h1;
                exact, closed_line = true,
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
            m, h2, poly1;
            exact, closed_line = true,
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
Classify a point against a segment: `on` at an endpoint, `in` in the segment interior, or
`out` elsewhere. Keyword values set the result for each case.

The inputs must have point and line string or linear ring traits.
=#
function _point_segment_orientation(
    ::Planar, point, start, stop;
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
Classify a point against a filled curve: `in` for the interior, `on` for edges or vertices,
and `out` for the exterior. Keywords set the return values.

The curve must be a line string or linear ring. Treat it as closed regardless of a repeated
final vertex.

Uses Hao and Sun (2018), https://doi.org/10.3390/sym10100477. Boundary cases return `on`;
otherwise, an odd horizontal-ray crossing count means inside. Case labels follow the paper.
=#

function _point_filled_curve_orientation(
    ::Planar, point, curve;
    in::T = point_in, on::T = point_on, out::T = point_out, exact,
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
            u1, u2 = GI.x(p_start) - x, GI.x(p_end) - x
            f = Predicates.orient(p_start, p_end, (x, y); exact)
            if v2 > 0 && v1 ≤ 0                # Case 3, 9, 16, 21, 13, or 24
                f == 0 && return on         # Case 16 or 21
                f > 0 && (k += 1)              # Case 3 or 9
            elseif v1 > 0 && v2 ≤ 0            # Case 4, 10, 19, 20, 12, or 25
                f == 0 && return on         # Case 19 or 20
                f < 0 && (k += 1)              # Case 4 or 10
            elseif v2 == 0 && v1 < 0           # Case 7, 14, or 17
                f == 0 && return on         # Case 17
            elseif v1 == 0 && v2 < 0           # Case 8, 15, or 18
                f == 0 && return on         # Case 18
            elseif v1 == 0 && v2 == 0          # Case 1, 2, 5, 6, 22, or 23
                u2 ≤ 0 && u1 ≥ 0 && return on  # Case 1
                u1 ≤ 0 && u2 ≥ 0 && return on  # Case 2
            end
        end
        p_start = p_end
    end
    return iseven(k) ? out : in
end

# Specialized implementation for NaturallyIndexedRing
# This relies on multidispatch.
# TODO: remove?
function _point_filled_curve_orientation(
    ::Planar, point, curve::NaturalIndexing.NaturallyIndexedRing;
    in::T = point_in, on::T = point_on, out::T = point_out, exact,
) where {T}
    x, y = GI.x(GI.PointTrait(), point), GI.y(GI.PointTrait(), point)
    k::Int = 0  # counter for ray crossings

    tree = curve.index

    function per_edge_function(i)
        p_start = _tuple_point(GI.getpoint(curve, i))
        p_end = _tuple_point(GI.getpoint(curve, i + 1))
        v1 = GI.y(p_start) - y
        v2 = GI.y(p_end) - y
        if !((v1 < 0 && v2 < 0) || (v1 > 0 && v2 > 0)) # if not cases 11 or 26
            u1, u2 = GI.x(p_start) - x, GI.x(p_end) - x
            f = Predicates.orient(p_start, p_end, (x, y); exact)
            if v2 > 0 && v1 ≤ 0                # Case 3, 9, 16, 21, 13, or 24
                f == 0 && return LSM.Action{T}(:full_return, on)         # Case 16 or 21
                f > 0 && (k += 1)              # Case 3 or 9
            elseif v1 > 0 && v2 ≤ 0            # Case 4, 10, 19, 20, 12, or 25
                f == 0 && return LSM.Action{T}(:full_return, on)         # Case 19 or 20
                f < 0 && (k += 1)              # Case 4 or 10
            elseif v2 == 0 && v1 < 0           # Case 7, 14, or 17
                f == 0 && return LSM.Action{T}(:full_return, on)         # Case 17
            elseif v1 == 0 && v2 < 0           # Case 8, 15, or 18
                f == 0 && return LSM.Action{T}(:full_return, on)         # Case 18
            elseif v1 == 0 && v2 == 0          # Case 1, 2, 5, 6, 22, or 23
                u2 ≤ 0 && u1 ≥ 0 && return LSM.Action{T}(:full_return, on)  # Case 1
                u1 ≤ 0 && u2 ≥ 0 && return LSM.Action{T}(:full_return, on)  # Case 2
            end
            return LSM.Action(:continue, on)
        end
        p_start = p_end
    end

    result = SpatialTreeInterface.depth_first_search(per_edge_function,extent -> extent.Y[1] <= y <= extent.Y[2], tree)

    if result isa LoopStateMachine.Action
        return result.x
    else
        return iseven(k) ? out : in
    end
end
_point_filled_curve_orientation(
    point, curve;
    in::T = point_in, on::T = point_on, out::T = point_out, exact,
) where {T} = _point_filled_curve_orientation(Planar(), point, curve; in, on, out, exact)

#=
Return `(in_curve, on_curve, out_curve)` for a line against a filled curve.

`in_curve` marks interior contact; `on_curve` marks endpoint or boundary contact. `out_curve`
marks any segment outside the curve.

`closed_line` connects the first and last line vertices.
=#
function _line_filled_curve_interactions(
    m::Manifold, line, curve;
    exact, closed_line = false,
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
    l_start = _processor_point(m, GI.getpoint(line, closed_line ? nl : 1))
    point_val = _point_filled_curve_orientation(m, l_start, curve; exact)
    if point_val == point_in
        in_curve = true
    elseif point_val == point_on
        on_curve = true
    else  # point_val == point_out
        out_curve = true
    end

    # Check for any intersections between line and curve
    for i in (closed_line ? 1 : 2):nl
        l_end = _processor_point(m, GI.getpoint(line, i))
        c_start = _processor_point(m, GI.getpoint(curve, nc))
        # If already interacted with all regions of curve, can stop
        in_curve && on_curve && out_curve && break
        # Check next segment of line against curve
        for j in 1:nc
            c_end = _processor_point(m, GI.getpoint(curve, j))
            # Check if two line and curve segments meet
            seg_val, _, _ = _seg_seg_orientation(m, l_start, l_end, c_start, c_end; exact)
            if seg_val != line_out
                # If line and curve meet, then at least one point is on boundary
                on_curve = true
                if seg_val == line_cross
                    # When crossing boundary, line is both in and out of curve
                    in_curve = true
                    out_curve = true
                else
                    if seg_val == line_over
                        sp = _point_segment_orientation(m, l_start, c_start, c_end)
                        lp = _point_segment_orientation(m, l_end, c_start, c_end)
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
                        in_curve, out_curve = _split_segment_interactions(
                            m, l_start, l_end, curve, in_curve, out_curve; exact,
                        )
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
Return `(in_poly, on_poly, out_poly)` for a line against a polygon.

`in_poly` marks interior contact; `on_poly` marks endpoint or boundary contact, including
holes. `out_poly` marks line portions outside the polygon, including inside holes.

`closed_line` connects the first and last line vertices.
=#
function _line_polygon_interactions(
    m::Manifold, line, polygon;
    exact, closed_line = false,
)

    in_poly, on_poly, out_poly = _line_filled_curve_interactions(
        m, line, GI.getexterior(polygon);
        exact, closed_line = closed_line,
    )
    !in_poly && return (in_poly, on_poly, out_poly)
    # Loop over polygon holes
    for hole in GI.gethole(polygon)
        in_hole, on_hole, out_hole =_line_filled_curve_interactions(
            m, line, hole;
            exact, closed_line = closed_line,
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

# Disjoint extent optimisation: skip work based on geom extent intersection
# returns Tuple{Bool, Bool} for (skip, returnval)
@inline function _maybe_skip_disjoint_extents(::Planar, a, b;
    in_allow, on_allow, out_allow, 
    in_require, on_require, out_require,
    kw...
)
    ext_disjoint = Extents.disjoint(GI.extent(a), GI.extent(b))
    skip, returnval = if !ext_disjoint
        # can't tell anything about this case
        false, false
    elseif out_allow # && ext_disjoint
        if in_require || on_require
            true, false
        else
            true, true
        end
    else  # !out_allow && ext_disjoint
        # points not allowed in exterior, but geoms are disjoint
        true, false
    end
    return skip, returnval
end

#= Planar-defaulting forwarders for callers that do not specify a manifold,
matching the one `_point_filled_curve_orientation` already carries. =#
_line_filled_curve_interactions(line, curve; exact, closed_line = false) =
    _line_filled_curve_interactions(Planar(), line, curve; exact, closed_line)
_line_polygon_interactions(line, polygon; exact, closed_line = false) =
    _line_polygon_interactions(Planar(), line, polygon; exact, closed_line)

#=
Return `line_out`, `line_cross`, `line_hinge`, or `line_over`, plus intersection fractions `α`
along `(a1, a2)` and `β` along `(b1, b2)`.

Callers test fractions only against 0 and 1. Symbolic classifiers need only distinguish
endpoints from interior intersections.
=#
@inline function _seg_seg_orientation(m::Planar, a1, a2, b1, b2; exact)
    seg_val, intr1, _ = _intersection_point(m, Float64, (a1, a2), (b1, b2); exact)
    (_, (α, β)) = intr1
    return seg_val, α, β
end

#=
Split `l_start → l_end` at contacts with `curve`. OR each piece's interior/exterior status
into `in_curve` and `out_curve`.

This resolves hinges spanning multiple curve segments that edge-by-edge classification cannot
settle.
=#
function _split_segment_interactions(m::Planar, l_start, l_end, curve, in_curve, out_curve; exact)
    ipoints = intersection_points(GI.Line(StaticArrays.SVector(l_start, l_end)), curve)
    npoints = length(ipoints)  # since hinge, at least one
    dist_from_lstart = let l_start = l_start
        x -> _euclid_distance(Float64, x, l_start)
    end
    sort!(ipoints, by = dist_from_lstart)
    p_start = _tuple_point(l_start)
    for i in 1:(npoints + 1)
        p_end = i ≤ npoints ? _tuple_point(ipoints[i]) : l_end
        mid_val = _point_filled_curve_orientation(m, (p_start .+ p_end) ./ 2, curve; exact)
        if mid_val == point_in
            in_curve = true
        elseif mid_val == point_out
            out_curve = true
        end
        p_start = p_end
    end
    return in_curve, out_curve
end
