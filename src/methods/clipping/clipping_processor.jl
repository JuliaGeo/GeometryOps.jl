# # Polygon clipping helpers
# This file contains the shared helper functions for the polygon clipping functionalities.

# This enum defines which side of an edge a point is on
@enum PointEdgeSide left=1 right=2 unknown=3

# Constants assigned for readability
const enter, exit = true, false
const crossing, bouncing = true, false

#= A point can either be the start or end of an overlapping chain of points between two
polygons, or not an endpoint of a chain. =#
@enum EndPointType start_chain=1 end_chain=2 not_endpoint=3

#= This is the struct that makes up a_list and b_list. Many values are only used if point is
an intersection point (ipt). =#
@kwdef struct PolyNode{T <: AbstractFloat}
    point::Tuple{T,T}          # (x, y) values of given point
    inter::Bool = false        # If ipt, true, else 0
    neighbor::Int = 0          # If ipt, index of equivalent point in a_list or b_list, else 0
    idx::Int = 0               # If crossing point, index within sorted a_idx_list
    ent_exit::Bool = false     # If ipt, true if enter and false if exit, else false
    crossing::Bool = false     # If ipt, true if intersection crosses from out/in polygon, else false
    endpoint::EndPointType = not_endpoint # If ipt, denotes if point is the start or end of an overlapping chain
    fracs::Tuple{T,T} = (0., 0.) # If ipt, fractions along edges to ipt (a_frac, b_frac), else (0, 0)
end

#= Create a new node with all of the same field values unless alternative values are
provided, in which case those should be used. =#
_update_node(node::PolyNode{T};
    point = node.point, inter = node.inter, neighbor = node.neighbor, idx = node.idx,
    ent_exit = node.ent_exit, crossing = node.crossing, endpoint = node.endpoint,
    fracs = node.fracs,
) where T = PolyNode{T}(;
    point = point, inter = inter, neighbor = neighbor, idx = idx, ent_exit = ent_exit,
    crossing = crossing, endpoint = endpoint, fracs = fracs)

#=
    _build_ab_list(::Type{T}, poly_a, poly_b) -> (a_list, b_list, a_idx_list)

This function takes in two polygon rings and calls '_build_a_list', '_build_b_list', and
'_flag_ent_exit' in order to fully form a_list and b_list. The 'a_list' and 'b_list' that it
returns are the fully updated vectors of PolyNodes that represent the rings 'poly_a' and
'poly_b', respectively. This function also returns 'a_idx_list', which at its "ith" index
stores the index in 'a_list' at which the "ith" intersection point lies.
=#
function _build_ab_list(::Type{T}, poly_a, poly_b, delay_cross_f, delay_bounce_f) where T
    # Make a list for nodes of each polygon
    a_list, a_idx_list, n_b_intrs = _build_a_list(T, poly_a, poly_b)
    b_list = _build_b_list(T, a_idx_list, a_list, n_b_intrs, poly_b)

    # Flag crossings
    _classify_crossing!(T, a_list, b_list)

    # Flag the entry and exits
    _flag_ent_exit!(GI.LinearRingTrait(), poly_b, a_list, (x) -> delay_cross_f(x), (x) -> delay_bounce_f(x, true))
    _flag_ent_exit!(GI.LinearRingTrait(), poly_a, b_list, (x) -> delay_cross_f(x), (x) -> delay_bounce_f(x, false))

    # Set node indices and filter a_idx_list to just crossing points
    _index_crossing_intrs!(a_list, b_list, a_idx_list)

    return a_list, b_list, a_idx_list
end

#=
    _build_a_list(::Type{T}, poly_a, poly_b) -> (a_list, a_idx_list)

This function take in two polygon rings and creates a vector of PolyNodes to represent
poly_a, including its intersection points with poly_b. The information stored in each
PolyNode is needed for clipping using the Greiner-Hormann clipping algorithm.
    
Note: After calling this function, a_list is not fully formed because the neighboring
indicies of the intersection points in b_list still need to be updated. Also we still have
not update the entry and exit flags for a_list.
    
The a_idx_list is a list of the indicies of intersection points in a_list. The value at
index i of a_idx_list is the location in a_list where the ith intersection point lies.
=#
function _build_a_list(::Type{T}, poly_a, poly_b) where T
    n_a_edges = _nedge(poly_a)
    a_list = PolyNode{T}[]  # list of points in poly_a
    sizehint!(a_list, n_a_edges)
    a_idx_list = Vector{Int}()  # finds indices of intersection points in a_list
    a_count = 0  # number of points added to a_list
    n_b_intrs = 0
    # Loop through points of poly_a
    local a_pt1
    for (i, a_p2) in enumerate(GI.getpoint(poly_a))
        a_pt2 = (T(GI.x(a_p2)), T(GI.y(a_p2)))
        if i <= 1 || (a_pt1 == a_pt2)  # don't repeat points
            a_pt1 = a_pt2
            continue
        end
        # Add the first point of the edge to the list of points in a_list
        new_point = PolyNode{T}(;point = a_pt1)
        a_count += 1
        push!(a_list, new_point)
        # Find intersections with edges of poly_b
        local b_pt1
        prev_counter = a_count
        for (j, b_p2) in enumerate(GI.getpoint(poly_b))
            b_pt2 = _tuple_point(b_p2, T)
            if j <= 1 || (b_pt1 == b_pt2)  # don't repeat points
                b_pt1 = b_pt2
                continue
            end
            # Determine if edges intersect and how they intersect
            line_orient, intr1, intr2 = _intersection_point(T, (a_pt1, a_pt2), (b_pt1, b_pt2))
            if line_orient != line_out  # edges intersect
                if line_orient == line_cross  # Intersection point that isn't a vertex
                    int_pt, fracs = intr1
                    new_intr = PolyNode{T}(;
                        point = int_pt, inter = true, neighbor = j - 1,
                        crossing = true, fracs = fracs,
                    )
                    a_count += 1
                    n_b_intrs += 1
                    push!(a_list, new_intr)
                    push!(a_idx_list, a_count)
                else
                    (_, (α1, β1)) = intr1
                    # Determine if a1 or b1 should be added to a_list
                    add_a1 = α1 == 0 && 0 ≤ β1 < 1
                    a1_β = add_a1 ? β1 : zero(T)
                    add_b1 = β1 == 0 && 0 < α1 < 1
                    b1_α = add_b1 ? α1 : zero(T)
                    # If lines are collinear and overlapping, a second intersection exists
                    if line_orient == line_over
                        (_, (α2, β2)) = intr2
                        if α2 == 0 && 0 ≤ β2 < 1
                            add_a1, a1_β = true, β2
                        end
                        if β2 == 0 && 0 < α2 < 1
                            add_b1, b1_α = true, α2
                        end
                    end
                    # Add intersection points determined above
                    if add_a1
                        n_b_intrs += a1_β == 0 ? 0 : 1
                        a_list[prev_counter] = PolyNode{T}(;
                            point = a_pt1, inter = true, neighbor = j - 1,
                            fracs = (zero(T), a1_β),
                        )
                        push!(a_idx_list, prev_counter)
                    end
                    if add_b1
                        new_intr = PolyNode{T}(;
                            point = b_pt1, inter = true, neighbor = j - 1,
                            fracs = (b1_α, zero(T)),
                        )
                        a_count += 1
                        push!(a_list, new_intr)
                        push!(a_idx_list, a_count)
                    end
                end
            end
            b_pt1 = b_pt2
        end
        # Order intersection points by placement along edge using fracs value
        if prev_counter < a_count
            Δintrs = a_count - prev_counter
            inter_points = @view a_list[(a_count - Δintrs + 1):a_count]
            sort!(inter_points, by = x -> x.fracs[1])
        end
        a_pt1 = a_pt2
    end
    return a_list, a_idx_list, n_b_intrs
end

#=
    _build_b_list(::Type{T}, a_idx_list, a_list, poly_b) -> b_list

This function takes in the a_list and a_idx_list build in _build_a_list and poly_b and
creates a vector of PolyNodes to represent poly_b. The information stored in each PolyNode
is needed for clipping using the Greiner-Hormann clipping algorithm.
    
Note: after calling this function, b_list is not fully updated. The entry/exit flags still
need to be updated. However, the neightbor value in a_list is now updated.
=#
function _build_b_list(::Type{T}, a_idx_list, a_list, n_b_intrs, poly_b) where T
    # Sort intersection points by insertion order in b_list
    sort!(a_idx_list, by = x-> a_list[x].neighbor + a_list[x].fracs[2])
    # Initialize needed values and lists
    n_b_edges = _nedge(poly_b)
    n_intr_pts = length(a_idx_list)
    b_list = PolyNode{T}[]
    sizehint!(b_list, n_b_edges + n_b_intrs)
    intr_curr = 1
    b_count = 0
    # Loop over points in poly_b and add each point and intersection point
    local b_pt1
    for (i, b_p2) in enumerate(GI.getpoint(poly_b))
        b_pt2 = _tuple_point(b_p2, T)
        if i ≤ 1 || (b_pt1 == b_pt2)  # don't repeat points
            b_pt1 = b_pt2
            continue
        end
        b_count += 1
        push!(b_list, PolyNode{T}(; point = b_pt1))
        if intr_curr ≤ n_intr_pts
            curr_idx = a_idx_list[intr_curr]
            curr_node = a_list[curr_idx]
            prev_counter = b_count
            while curr_node.neighbor == i - 1  # Add all intersection points on current edge
                b_idx = 0
                new_intr = _update_node(curr_node; neighbor = curr_idx)
                if equals(curr_node.point, b_list[prev_counter].point)
                    # intersection point is vertex of b
                    b_idx = prev_counter
                    b_list[b_idx] = new_intr
                else
                    b_count += 1
                    b_idx = b_count
                    push!(b_list, new_intr)
                end
                a_list[curr_idx] = _update_node(curr_node; neighbor = b_idx)
                intr_curr += 1
                intr_curr > n_intr_pts && break
                curr_idx = a_idx_list[intr_curr]
                curr_node = a_list[curr_idx]
            end
        end
        b_pt1 = b_pt2
    end
    sort!(a_idx_list)  # return a_idx_list to order of points in a_list
    return b_list
end

#=
    _classify_crossing!(T, poly_b, a_list)

This function marks all intersection points as either bouncing or crossing points. "Delayed"
crossing or bouncing intersections (a chain of edges where the central edges overlap and
thus only the first and last edge of the chain determine if the chain is bounding or
crossing) are marked as follows: the first and the last points are marked as crossing if the
chain is crossing and delayed otherwise and all middle points are marked as bouncing.
Additionally, the start and end points of the chain are marked as endpoints using the
endpoints field. 
=#
function _classify_crossing!(::Type{T}, a_list, b_list) where T
    napts = length(a_list)
    nbpts = length(b_list)
    # start centered on last point
    a_prev = a_list[end - 1]
    curr_pt = a_list[end]
    i = napts
    # keep track of unmatched bouncing chains
    start_chain_edge, start_chain_idx = unknown, 0
    unmatched_end_chain_edge, unmatched_end_chain_idx = unknown, 0
    same_winding = true
    # loop over list points
    for next_idx in 1:napts
        a_next = a_list[next_idx]
        if curr_pt.inter && !curr_pt.crossing
            j = curr_pt.neighbor
            b_prev = j == 1 ? b_list[end] : b_list[j-1]
            b_next = j == nbpts ? b_list[1] : b_list[j+1]
            # determine if any segments are on top of one another
            a_prev_is_b_prev = a_prev.inter && a_prev.point == b_prev.point
            a_prev_is_b_next = a_prev.inter && a_prev.point == b_next.point
            a_next_is_b_prev = a_next.inter && a_next.point == b_prev.point
            a_next_is_b_next = a_next.inter && a_next.point == b_next.point
            # determine which side of a segments the p points are on
            b_prev_side = _get_side(b_prev.point, a_prev.point, curr_pt.point, a_next.point)
            b_next_side = _get_side(b_next.point, a_prev.point, curr_pt.point, a_next.point)
            # no sides overlap
            if !a_prev_is_b_prev && !a_prev_is_b_next && !a_next_is_b_prev && !a_next_is_b_next
                if b_prev_side != b_next_side  # lines cross 
                    a_list[i] = _update_node(curr_pt; crossing = true)
                    b_list[j] = _update_node(b_list[j]; crossing = true)
                end
            # end of overlapping chain
            elseif !a_next_is_b_prev && !a_next_is_b_next 
                b_side = a_prev_is_b_prev ? b_next_side : b_prev_side
                if start_chain_edge == unknown  # start loop on overlapping chain
                    unmatched_end_chain_edge = b_side
                    unmatched_end_chain_idx = i
                    same_winding = a_prev_is_b_prev
                else  # close overlapping chain
                    # update end of chain with endpoint and crossing / bouncing tags
                    crossing = b_side != start_chain_edge
                    a_list[i] = _update_node(curr_pt;
                        crossing = crossing,
                        endpoint = end_chain,
                    )
                    b_list[j] = _update_node(b_list[j];
                        crossing = crossing,
                        endpoint = same_winding ? end_chain : start_chain,
                    )
                    # update start of chain with endpoint and crossing / bouncing tags
                    start_pt = a_list[start_chain_idx]
                    a_list[start_chain_idx] = _update_node(start_pt;
                        crossing = crossing,
                        endpoint = start_chain,
                    )
                    b_list[start_pt.neighbor] = _update_node(b_list[start_pt.neighbor];
                        crossing = crossing,
                        endpoint = same_winding ? start_chain : end_chain,
                    )
                end
            # start of overlapping chain
            elseif !a_prev_is_b_prev && !a_prev_is_b_next
                b_side = a_next_is_b_prev ? b_next_side : b_prev_side
                start_chain_edge = b_side
                start_chain_idx = i
                same_winding = a_next_is_b_next
            end
        end
        a_prev = curr_pt
        curr_pt = a_next
        i = next_idx
    end
    # if we started in the middle of overlapping chain, close chain
    if unmatched_end_chain_edge != unknown
        crossing = unmatched_end_chain_edge != start_chain_edge
        # update end of chain with endpoint and crossing / bouncing tags
        end_chain_pt = a_list[unmatched_end_chain_idx]
        a_list[unmatched_end_chain_idx] = _update_node(end_chain_pt;
            crossing = crossing,
            endpoint = end_chain,
        )
        b_list[end_chain_pt.neighbor] = _update_node(b_list[end_chain_pt.neighbor];
            crossing = crossing,
            endpoint = same_winding ? end_chain : start_chain,
        )
        # update start of chain with endpoint and crossing / bouncing tags
        start_pt = a_list[start_chain_idx]
        a_list[start_chain_idx] = _update_node(start_pt;
            crossing = crossing,
            endpoint = start_chain,
        )
        b_list[start_pt.neighbor] = _update_node(start_pt;
            crossing = crossing,
            endpoint = same_winding ? start_chain : end_chain,
        )
    end
end

# Determines if Q lies to the left or right of the line formed by P1-P2-P3
function _get_side(Q, P1, P2, P3)
    s1 = _signed_area_triangle(Q, P1, P2)
    s2 = _signed_area_triangle(Q, P2, P3)
    s3 = _signed_area_triangle(P1, P2, P3)

    side = if s3 ≥ 0
        (s1 < 0) || (s2 < 0) ? right : left
    else #  s3 < 0
        (s1 > 0) || (s2 > 0) ? left : right
    end
    return side
end

# Returns the signed area formed by vertices P, Q, and R
function _signed_area_triangle(P, Q, R)
    return (GI.x(Q)-GI.x(P))*(GI.y(R)-GI.y(P))-(GI.y(Q)-GI.y(P))*(GI.x(R)-GI.x(P))
end

# True if the edge with pt as the starting endpoint is not shared between polygons
_next_edge_off(pt) = !pt.inter || (pt.endpoint == end_chain) || (pt.crossing && pt.endpoint == not_endpoint)

#=
    _flag_ent_exit!(::GI.LinearRingTrait, poly, pt_list, delay_cross_f, delay_bounce_f)

This function flags all the intersection points as either an 'entry' or 'exit' point in
relation to the given polygon. For non-delayed crossings we simply alternate the enter/exit
status. This also holds true for the first and last points of a delayed bouncing, where they
both have an opposite entry/exit flag. Conversely, the first and last point of a delayed
crossing have the same entry/exit status. Furthermore, the crossing/bouncing flag of delayed
crossings and bouncings may be updated. This depends on function specific rules that
determine which of the start or end points (if any) should be marked as crossing for used
during polygon tracing. A consistent rule is that the start and end points of a delayed
crossing will have different crossing/bouncing flags, while a the endpoints of a delayed
bounce will be the same.

Used for clipping polygons by other polygons.
=#
function _flag_ent_exit!(::GI.LinearRingTrait, poly, pt_list, delay_cross_f, delay_bounce_f)
    # Find starting index if there is one
    start_idx = findfirst(_next_edge_off, pt_list)
    isnothing(start_idx) && return
    # Determine if non-overlapping line midpoint is inside or outside of polygon
    npts = length(pt_list)
    next_idx = start_idx < npts ? (start_idx + 1) : 1
    start_val = (pt_list[start_idx].point .+ pt_list[next_idx].point) ./ 2
    status = !_point_filled_curve_orientation(start_val, poly; in = true, on = false, out = false)
    # Loop over points and mark entry and exit status
    start_chain_idx = 0
    for ii in Iterators.flatten((next_idx:npts, 1:start_idx))
        curr_pt = pt_list[ii]
        if curr_pt.endpoint == start_chain
            start_chain_idx = ii
        elseif curr_pt.crossing || curr_pt.endpoint == end_chain
            start_crossing, end_crossing = curr_pt.crossing, curr_pt.crossing
            if curr_pt.endpoint == end_chain  # ending overlapping chain
                start_pt = pt_list[start_chain_idx]
                if curr_pt.crossing  # delayed crossing
                    #= start and end crossing status are different and depend on current
                    entry/exit status =#
                    start_crossing, end_crossing = delay_cross_f(status)
                else  # delayed bouncing
                    next_idx = ii < npts ? (ii + 1) : 1
                    next_val = (curr_pt.point .+ pt_list[next_idx].point) ./ 2
                    pt_in_poly = _point_filled_curve_orientation(next_val, poly; in = true, on = false, out = false)
                    #= start and end crossing status are the same and depend on if adjacent
                    edges of pt_list are within poly =#
                    start_crossing = delay_bounce_f(pt_in_poly)
                    end_crossing = start_crossing
                end
                # update start of chain point
                pt_list[start_chain_idx] = _update_node(start_pt; ent_exit = status, crossing = start_crossing)
                if !curr_pt.crossing
                    status = !status
                end
            end
            pt_list[ii] = _update_node(curr_pt; ent_exit = status, crossing = end_crossing)
            status = !status
        end
    end
    return
end

#=
    _flag_ent_exit!(::GI.LineTrait, line, pt_list)

This function flags all the intersection points as either an 'entry' or 'exit' point in
relation to the given line. Returns true if there are crossing points to classify, else
returns false. Used for cutting polygons by lines.

Assumes that the first point is outside of the polygon and not on an edge.
=#
function _flag_ent_exit!(::GI.LineTrait, poly, pt_list)
    status = !_point_filled_curve_orientation(pt_list[1].point, poly; in = true, on = false, out = false)
    # Loop over points and mark entry and exit status
    for (ii, curr_pt) in enumerate(pt_list)
        if curr_pt.crossing
            pt_list[ii] = _update_node(curr_pt; ent_exit = status)
            status = !status
        end
    end
    return
end

#= Filters a_idx_list to just include crossing points and sets the index of all crossing
points (which element they correspond to within a_idx_list). =#
function _index_crossing_intrs!(a_list, b_list, a_idx_list)
    filter!(x -> a_list[x].crossing, a_idx_list)
    for (i, a_idx) in enumerate(a_idx_list)
        curr_node = a_list[a_idx]
        neighbor_node = b_list[curr_node.neighbor]
        a_list[a_idx] = _update_node(curr_node; idx = i)
        b_list[curr_node.neighbor] = _update_node(neighbor_node; idx = i)
    end
    return
end

#=
    _trace_polynodes(::Type{T}, a_list, b_list, a_idx_list, f_step)::Vector{GI.Polygon}

This function takes the outputs of _build_ab_list and traces the lists to determine which
polygons are formed as described in Greiner and Hormann. The function f_step determines in
which direction the lists are traced.  This function is different for intersection,
difference, and union. f_step must take in two arguments: the most recent intersection
node's entry/exit status and a boolean that is true if we are currently tracing a_list and
false if we are tracing b_list. The functions used for each clipping operation are follows:
    - Intersection: (x, y) -> x ? 1 : (-1)
    - Difference: (x, y) -> (x ⊻ y) ? 1 : (-1)
    - Union: (x, y) -> x ? (-1) : 1

A list of GeoInterface polygons is returned from this function. 
=#
function _trace_polynodes(::Type{T}, a_list, b_list, a_idx_list, f_step) where T
    n_a_pts, n_b_pts = length(a_list), length(b_list)
    n_cross_pts = length(a_idx_list)
    return_polys = Vector{_get_poly_type(T)}(undef, 0)
    # Keep track of number of processed intersection points
    processed_pts = 0
    first_idx = 1
    while processed_pts < n_cross_pts
        curr_list, curr_npoints = a_list, n_a_pts
        on_a_list = true
        # Find first unprocessed intersecting point in subject polygon
        processed_pts += 1
        first_idx = findnext(x -> x != 0, a_idx_list, first_idx)
        idx = a_idx_list[first_idx]
        a_idx_list[first_idx] = 0
        start_pt = a_list[idx]

        # Set first point in polygon
        curr = curr_list[idx]
        pt_list = [curr.point]

        curr_not_start = true
        while curr_not_start
            step = f_step(curr.ent_exit, on_a_list)
            # changed curr_not_intr to curr_not_same_ent_flag
            same_status, prev_status = true, curr.ent_exit
            while same_status
                # Traverse polygon either forwards or backwards
                idx += step
                idx = (idx > curr_npoints) ? mod(idx, curr_npoints) : idx
                idx = (idx == 0) ? curr_npoints : idx

                # Get current node and add to pt_list
                curr = curr_list[idx]
                push!(pt_list, curr.point)
                if (curr.crossing || curr.endpoint != not_endpoint)
                    # Keep track of processed intersection points
                    same_status = curr.ent_exit == prev_status
                    curr_not_start = curr != start_pt && curr != b_list[start_pt.neighbor]
                    !curr_not_start && break
                    if (on_a_list && curr.crossing) || (!on_a_list && a_list[curr.neighbor].crossing)
                        processed_pts += 1
                        a_idx_list[curr.idx] = 0
                    end
                end
            end
            # Switch to next list and next point
            curr_list, curr_npoints = on_a_list ? (b_list, n_b_pts) : (a_list, n_a_pts)
            on_a_list = !on_a_list
            idx = curr.neighbor
            curr = curr_list[idx]
        end
        push!(return_polys, GI.Polygon([pt_list]))
    end
    return return_polys
end

# Get type of polygons that will be made
# TODO: Increase type options
_get_poly_type(::Type{T}) where T =
    GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{Tuple{T, T}}, Nothing, Nothing}}, Nothing, Nothing}

#=
    _find_non_cross_orientation(a_list, b_list, a_poly, b_poly)

For polygns with no crossing intersection points, either one polygon is inside of another,
or they are seperate polygons with no intersection (other than an edge or point).

Return two booleans that represent if a is inside b (potentially with shared edges / points)
and visa versa if b is inside of a.
=#
function _find_non_cross_orientation(a_list, b_list, a_poly, b_poly)
    non_intr_a_idx = findfirst(x -> !x.inter, a_list)
    non_intr_b_idx = findfirst(x -> !x.inter, b_list)
    #= Determine if non-intersection point is in or outside of polygon - if there isn't A
    non-intersection point, then all points are on the polygon edge =#
    a_pt_orient = isnothing(non_intr_a_idx) ? point_on :
        _point_filled_curve_orientation(a_list[non_intr_a_idx].point, b_poly)
    b_pt_orient = isnothing(non_intr_b_idx) ? point_on :
        _point_filled_curve_orientation(b_list[non_intr_b_idx].point, a_poly)
    a_in_b = a_pt_orient != point_out && b_pt_orient != point_in
    b_in_a = b_pt_orient != point_out && a_pt_orient != point_in
    return a_in_b, b_in_a
end

#=
    _add_holes_to_polys!(::Type{T}, return_polys, hole_iterator)

The holes specified by the hole iterator are added to the polygons in the return_polys list.
If this creates more polygon, they are added to the end of the list. If this removes
polygons, they are removed from the list
=#
function _add_holes_to_polys!(::Type{T}, return_polys, hole_iterator) where T
    n_polys = length(return_polys)
    remove_poly_idx = fill(false, n_polys)
    remove_hole_idx = Int[]
    # Remove set of holes from all polygons
    for i in 1:n_polys
        n_new_per_poly = 0
        for curr_hole in hole_iterator # loop through all holes
            # loop through all pieces of original polygon (new pieces added to end of list)
            for j in Iterators.flatten((i:i, (n_polys + 1):(n_polys + n_new_per_poly)))
                curr_poly = return_polys[j]
                remove_poly_idx[j] && continue
                n_existing_holes = GI.nhole(curr_poly)
                curr_poly_ext = n_existing_holes > 0 ? GI.Polygon(StaticArrays.SVector(GI.getexterior(curr_poly))) : curr_poly
                in_ext, on_ext, out_ext = _line_polygon_interactions(curr_hole, curr_poly_ext; closed_line = true)
                if in_ext  # hole is at least partially within the polygon's exterior
                    new_hole, new_hole_poly, n_new_pieces = _combine_holes!(T, curr_hole, curr_poly, return_polys, remove_hole_idx)
                    if n_new_pieces > 0
                        append!(remove_poly_idx, fill(false, n_new_pieces))
                        n_new_per_poly += n_new_pieces
                    end
                    if !on_ext && !out_ext  # hole is completly within exterior
                        push!(curr_poly.geom, new_hole)
                    else  # hole is partially within and outside of polygon's exterior
                        new_polys = difference(curr_poly_ext, new_hole_poly, T; target=GI.PolygonTrait())
                        n_new_polys = length(new_polys) - 1
                        # replace original -> can't have a hole
                        curr_poly.geom[1] = GI.getexterior(new_polys[1])
                        if n_new_polys > 0  # add any extra pieces
                            append!(return_polys, @view new_polys[2:end])
                            append!(remove_poly_idx, fill(false, n_new_polys))
                            n_new_per_poly += n_new_polys
                        end
                    end
                # polygon is completly within hole
                elseif coveredby(curr_poly_ext, GI.Polygon(StaticArrays.SVector(curr_hole)))
                    remove_poly_idx[j] = true
                end
            end
        end
        n_polys += n_new_per_poly
    end
    # Remove all polygon that were marked for removal
    deleteat!(return_polys, remove_poly_idx)
    return
end

#=
    _combine_holes!(::Type{T}, new_hole, curr_poly, return_polys)

The new hole is combined with any existing holes in curr_poly. The holes can be combined
into a larger hole if they are intersecting. If this happens, then the new, combined hole is
returned with the orignal holes making up the new hole removed from curr_poly. Additionally,
if the combined holes form a ring, the interior is added to the return_polys as a new
polygon piece. Additionally, holes leftover after combination will be checked for it they
are in the "main" polygon or in one of these new pieces and moved accordingly. 

If the holes don't touch or curr_poly has no holes, then new_hole is returned without any
changes.
=#
function _combine_holes!(::Type{T}, new_hole, curr_poly, return_polys, remove_hole_idx) where T
    n_new_polys = 0
    empty!(remove_hole_idx)
    new_hole_poly = GI.Polygon(StaticArrays.SVector(new_hole))
    # Combine any existing holes in curr_poly with new hole
    for (k, old_hole) in enumerate(GI.gethole(curr_poly))
        old_hole_poly = GI.Polygon(StaticArrays.SVector(old_hole))
        if intersects(new_hole_poly, old_hole_poly)
            # If the holes intersect, combine them into a bigger hole
            hole_union = union(new_hole_poly, old_hole_poly, T; target = GI.PolygonTrait())[1]
            push!(remove_hole_idx, k + 1)
            new_hole = GI.getexterior(hole_union)
            new_hole_poly = GI.Polygon(StaticArrays.SVector(new_hole))
            n_pieces = GI.nhole(hole_union)
            if n_pieces > 0  # if the hole has a hole, then this is a new polygon piece! 
                append!(return_polys, [GI.Polygon([h]) for h in GI.gethole(hole_union)])
                n_new_polys += n_pieces
            end
        end
    end
    # Remove redundant holes
    deleteat!(curr_poly.geom, remove_hole_idx)
    empty!(remove_hole_idx)
    # If new polygon pieces created, make sure remaining holes are in the correct piece
    @views for piece in return_polys[end - n_new_polys + 1:end]
        for (k, old_hole) in enumerate(GI.gethole(curr_poly))
            if !(k in remove_hole_idx) && within(old_hole, piece)
                push!(remove_hole_idx, k + 1)
                push!(piece.geom, old_hole)
            end
        end
    end
    deleteat!(curr_poly.geom, remove_hole_idx)
    return new_hole, new_hole_poly, n_new_polys
end
