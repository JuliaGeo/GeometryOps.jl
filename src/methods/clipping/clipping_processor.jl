# # Polygon clipping helpers
# This file contains the shared helper functions for the polygon clipping functionalities.

# This enum defines which side of an edge a point is on
@enum PointEdgeSide left=1 right=2 unknown=3

const enter, exit = true, false
const crossing, bouncing = true, false
@enum EndPointType start_chain=1 end_chain=2 not_endpoint=3

#= This is the struct that makes up a_list and b_list. Many values are only used if point is
an intersection point (ipt). =#
@kwdef struct PolyNode{T <: AbstractFloat}
    point::Tuple{T,T}          # (x, y) values of given point
    inter::Bool = false        # If ipt, true, else 0
    neighbor::Int = 0          # If ipt, index of equivalent point in a_list or b_list, else 0
    ent_exit::Bool = false     # If ipt, true if enter and false if exit, else false
    crossing::Bool = false     # If ipt, true if intersection crosses from out/in polygon, else false
    endpoint::EndPointType = not_endpoint # If ipt, true if point is the start of end of an overlapping chain
    fracs::Tuple{T,T} = (0., 0.) # If ipt, fractions along edges to ipt (a_frac, b_frac), else (0, 0)
end

#=
    _build_ab_list(::Type{T}, poly_a, poly_b) -> (a_list, b_list, a_idx_list)

This function takes in two polygon rings and calls '_build_a_list', '_build_b_list', and
'_flag_ent_exit' in order to fully form a_list and b_list. The 'a_list' and 'b_list' that it
returns are the fully updated vectors of PolyNodes that represent the rings 'poly_a' and
'poly_b', respectively. This function also returns 'a_idx_list', which at its "ith" index
stores the index in 'a_list' at which the "ith" intersection point lies.
=#
function _build_ab_list(::Type{T}, poly_a, poly_b) where T
    # Make a list for nodes of each polygon
    a_list, a_idx_list, n_b_intrs = _build_a_list(T, poly_a, poly_b)
    b_list = _build_b_list(T, a_idx_list, a_list, n_b_intrs, poly_b)
    # Flag crossings
    _classify_crossing!(T, a_list, b_list)

    # Flag the entry and exits
    _flag_ent_exit!(GI.LinearRingTrait(), poly_b, a_list)
    _flag_ent_exit!(GI.LinearRingTrait(), poly_a, b_list)

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
    a_list = Vector{PolyNode{T}}(undef, n_a_edges)  # list of points in poly_a
    a_idx_list = Vector{Int}()  # finds indices of intersection points in a_list
    a_count = 0  # number of points added to a_list
    n_b_intrs = 0
    # Loop through points of poly_a
    local a_pt1
    for (i, a_p2) in enumerate(GI.getpoint(poly_a))
        a_pt2 = (T(GI.x(a_p2)), T(GI.y(a_p2)))
        if i <= 1
            a_pt1 = a_pt2
            continue
        end
        # Add the first point of the edge to the list of points in a_list
        new_point = PolyNode{T}(;point = a_pt1)
        a_count += 1
        _add!(a_list, a_count, new_point, n_a_edges)
        # Find intersections with edges of poly_b
        local b_pt1
        prev_counter = a_count
        for (j, b_p2) in enumerate(GI.getpoint(poly_b))
            b_pt2 = _tuple_point(b_p2)
            if j <=1
                b_pt1 = b_pt2
                continue
            end
            int_pt, fracs = _intersection_point(T, (a_pt1, a_pt2), (b_pt1, b_pt2))
            if !isnothing(fracs)
                α, β = fracs
                collinear = isnothing(int_pt)
                # if no intersection point, skip this edge
                if !collinear && 0 < α < 1 && 0 < β < 1
                    # Intersection point that isn't a vertex
                    new_intr = PolyNode{T}(;
                        point = int_pt, inter = true, neighbor = j - 1,
                        crossing = true, fracs = fracs,
                    )
                    a_count += 1
                    n_b_intrs += 1
                    _add!(a_list, a_count, new_intr, n_a_edges)
                    push!(a_idx_list, a_count)
                else
                    if (0 < β < 1 && (collinear || α == 0)) || (α == β == 0)
                        # a_pt1 is an intersection point
                        n_b_intrs += β == 0 ? 0 : 1
                        a_list[prev_counter] = PolyNode{T}(;
                            point = a_pt1, inter = true, neighbor = j - 1,
                            fracs = fracs,
                        )
                        push!(a_idx_list, prev_counter)
                    end
                    if (0 < α < 1 && (collinear || β == 0))
                        # b_pt1 is an intersection point
                        new_intr = PolyNode{T}(;
                            point = b_pt1, inter = true, neighbor = j - 1,
                            fracs = fracs,
                        )
                        a_count += 1
                        _add!(a_list, a_count, new_intr, n_a_edges)
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

# Add value x at index i to given array - if list isn't long enough, push value to array
function _add!(arr::T, i, x, l = length(arr)) where {T <: Vector{<:PolyNode}}
    if i <= l
        arr[i] = x
    else
        push!(arr, x)
    end
    return
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
    b_list = Vector{PolyNode{T}}(undef, n_b_edges + n_b_intrs)
    intr_curr = 1
    b_count = 0
    # Loop over points in poly_b and add each point and intersection point
    for (i, p) in enumerate(GI.getpoint(poly_b))
        (i == n_b_edges + 1) && break
        b_count += 1
        pt = (T(GI.x(p)), T(GI.y(p)))
        b_list[b_count] = PolyNode(;point = pt)
        if intr_curr ≤ n_intr_pts
            curr_idx = a_idx_list[intr_curr]
            curr_node = a_list[curr_idx]
            prev_counter = b_count
            while curr_node.neighbor == i  # Add all intersection points in current edge
                b_idx = if equals(curr_node.point, b_list[prev_counter].point)
                    # intersection point is vertex of b
                    prev_counter
                else
                    b_count += 1
                    b_count
                end
                b_list[b_idx] = PolyNode{T}(;
                    point = curr_node.point, inter = true, neighbor = curr_idx,
                    crossing = curr_node.crossing, fracs = curr_node.fracs,
                )
                a_list[curr_idx] = PolyNode{T}(;
                    point = curr_node.point, inter = true, neighbor = b_idx,
                    crossing = curr_node.crossing, fracs = curr_node.fracs,
                )
                intr_curr += 1
                intr_curr > n_intr_pts && break
                curr_idx = a_idx_list[intr_curr]
                curr_node = a_list[curr_idx]
            end
        end
    end
    sort!(a_idx_list)  # return a_idx_list to order of points in a_list
    return b_list
end

#=
    _classify_crossing!(T, poly_b, a_list)

This function marks all intersection points as either bouncing or crossing points.
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
                    a_list[i] = PolyNode{T}(;
                        point = curr_pt.point, inter = true, neighbor = j,
                        crossing = true, fracs = curr_pt.fracs,
                    )
                    b_list[j] = PolyNode{T}(;
                        point = curr_pt.point, inter = true, neighbor = i,
                        crossing = true, fracs = curr_pt.fracs,
                    )
                end
            # end of overlapping chain
            elseif !a_next_is_b_prev && !a_next_is_b_next 
                b_side = a_prev_is_b_prev ? b_next_side : b_prev_side
                if start_chain_edge == unknown  # start loop on overlapping chain
                    unmatched_end_chain_edge = b_side
                    unmatched_end_chain_idx = i
                else  # close overlapping chain
                    # update end of chain with endpoint and crossing / bouncing tags
                    crossing = b_side != start_chain_edge
                    a_list[i] = PolyNode{T}(;
                        point = curr_pt.point, inter = true, neighbor = j,
                        crossing = crossing, endpoint = end_chain, fracs = curr_pt.fracs,
                    )
                    b_list[j] = PolyNode{T}(;
                        point = curr_pt.point, inter = true, neighbor = i,
                        crossing = crossing, endpoint = end_chain, fracs = curr_pt.fracs,
                    )
                    # update start of chain with endpoint and crossing / bouncing tags
                    start_pt = a_list[start_chain_idx]
                    a_list[start_chain_idx] = PolyNode{T}(;
                        point = start_pt.point, inter = true, neighbor = start_pt.neighbor,
                        crossing = crossing, endpoint = start_chain, fracs = start_pt.fracs,
                    )
                    b_list[start_pt.neighbor] = PolyNode{T}(;
                        point = start_pt.point, inter = true, neighbor = start_chain_idx,
                        crossing = crossing, endpoint = start_chain, fracs = start_pt.fracs,
                    )
                end
            # start of overlapping chain
            elseif !a_prev_is_b_prev && !a_prev_is_b_next
                b_side = a_next_is_b_prev ? b_next_side : b_prev_side
                start_chain_edge = b_side
                start_chain_idx = i
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
        a_list[unmatched_end_chain_idx] = PolyNode{T}(;
            point = end_chain_pt.point, inter = true, neighbor = end_chain_pt.neighbor,
            crossing = crossing, endpoint = end_chain, fracs = end_chain_pt.fracs,
        )
        b_list[end_chain_pt.neighbor] = PolyNode{T}(;
            point = end_chain_pt.point, inter = true, neighbor = unmatched_end_chain_idx,
            crossing = crossing, endpoint = end_chain, fracs = end_chain_pt.fracs,
        )
        # update start of chain with endpoint and crossing / bouncing tags
        start_pt = a_list[start_chain_idx]
        a_list[start_chain_idx] = PolyNode{T}(;
            point = start_pt.point, inter = true, neighbor = start_pt.neighbor,
            crossing = crossing, endpoint = start_chain, fracs = start_pt.fracs,
        )
        b_list[start_pt.neighbor] = PolyNode{T}(;
            point = start_pt.point, inter = true, neighbor = start_chain_idx,
            crossing = crossing, endpoint = start_chain, fracs = start_pt.fracs,
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
    _flag_ent_exit!(::GI.LinearRingTrait, poly_b, a_list)

This function flags all the intersection points as either an 'entry' or 'exit' point in
relation to the given polygon. Returns true if there are crossing points to classify, else
returns false. Used for clipping polygons by other polygons.
=#
function _flag_ent_exit!(::GI.LinearRingTrait, poly, pt_list)
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
            if curr_pt.endpoint == end_chain
                start_pt = pt_list[start_chain_idx]
                if curr_pt.crossing
                    start_crossing, end_crossing = !status, status
                else
                    next_idx = ii < npts ? (ii + 1) : 1
                    next_val = (curr_pt.point .+ pt_list[next_idx].point) ./ 2
                    start_crossing = _point_filled_curve_orientation(next_val, poly; in = true, on = false, out = false)
                    end_crossing = start_crossing
                end
                pt_list[start_chain_idx] = PolyNode(;
                    point = start_pt.point, inter = start_pt.inter, neighbor = start_pt.neighbor,
                    ent_exit = status, crossing = start_crossing, endpoint = start_pt.endpoint,
                    fracs = start_pt.fracs,
                )
                if !curr_pt.crossing
                    status = !status
                end
            end
            pt_list[ii] = PolyNode(;
                point = curr_pt.point, inter = curr_pt.inter, neighbor = curr_pt.neighbor,
                ent_exit = status, crossing = end_crossing, endpoint = curr_pt.endpoint,
                fracs = curr_pt.fracs,
            )
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
            pt_list[ii] = PolyNode(;
                point = curr_pt.point, inter = curr_pt.inter, neighbor = curr_pt.neighbor,
                ent_exit = status, crossing = curr_pt.crossing, fracs = curr_pt.fracs)
            status = !status
        end
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
    # Determine number of crossing intersection points
    n_cross_pts = 0
    for i in eachindex(a_idx_list)
        if a_list[a_idx_list[i]].crossing
            n_cross_pts += 1
        else
            a_idx_list[i] = 0
        end
    end

    return_polys = Vector{_get_poly_type(T)}(undef, 0)
    # Keep track of number of processed intersection points
    processed_pts = 0
    while processed_pts < n_cross_pts
        curr_list, curr_npoints = a_list, n_a_pts
        on_a_list = true
        # Find first unprocessed intersecting point in subject polygon
        processed_pts += 1
        first_idx = findnext(x -> x != 0, a_idx_list, processed_pts)
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
                    if curr.crossing && curr_not_start
                        processed_pts += 1
                        for (i, a_idx) in enumerate(a_idx_list)
                            if a_idx != 0 && equals(a_list[a_idx].point, curr.point)
                                a_idx_list[i] = 0
                                break
                            end
                        end
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

#= Determines if polygons share an edge (in the case where polygons are inside or outside
of one another and only commected by single points or edges) - if they share an edge,
print error message. =#
function share_edge_warn(list, warn_str)
    shared_edge = false
    prev_pt_inter = false
    for pt in list
        shared_edge = prev_pt_inter && pt.inter
        shared_edge && break
        prev_pt_inter = pt.inter
    end
    shared_edge && @warn warn_str
end

#=
    _add_holes_to_polys!(::Type{T}, return_polys, hole_iterator)

The holes specified by the hole iterator are added to the polygons in the return_polys list.
If this creates more polygon, they are added to the end of the list. If this removes
polygons, they are removed from the list
=#
function _add_holes_to_polys!(::Type{T}, return_polys, hole_iterator) where T
    n_polys = length(return_polys)
    # Remove set of holes from all polygons
    for i in 1:n_polys
        n_new_per_poly = 0
        for curr_hole in hole_iterator # loop through all holes
            # loop through all pieces of original polygon (new pieces added to end of list)
            for j in Iterators.flatten((i:i, (n_polys + 1):(n_polys + n_new_per_poly)))
                curr_poly = return_polys[j]
                isnothing(curr_poly) && continue
                n_existing_holes = GI.nhole(curr_poly)
                curr_poly_ext = n_existing_holes > 0 ? GI.Polygon([GI.getexterior(curr_poly)]) : curr_poly
                in_ext, on_ext, out_ext = _line_polygon_interactions(curr_hole, curr_poly_ext; closed_line = true)
                if in_ext  # hole is at least partially within the polygon's exterior
                    new_hole, new_hole_poly, n_new_pieces = _combine_holes!(T, curr_hole, curr_poly, return_polys)
                    n_new_per_poly += n_new_pieces
                    if !on_ext && !out_ext  # hole is completly within exterior
                        push!(curr_poly.geom, new_hole)
                    else  # hole is partially within and outside of polygon's exterior
                        new_polys = difference(curr_poly_ext, new_hole_poly, T; target = GI.PolygonTrait)
                        n_new_polys = length(new_polys) - 1
                        # replace original -> can't have a hole
                        curr_poly.geom[1] = GI.getexterior(new_polys[1])
                        if n_new_polys > 0  # add any extra pieces
                            append!(return_polys, @view new_polys[2:end])
                            n_new_per_poly += n_new_polys
                        end
                    end
                # polygon is completly within hole
                elseif coveredby(curr_poly_ext, GI.Polygon([curr_hole]))
                    return_polys[j] = nothing
                end
            end
        end
        n_polys += n_new_per_poly
    end
    # Remove all polygon that were marked for removal
    filter!(!isnothing, return_polys)
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
function _combine_holes!(::Type{T}, new_hole, curr_poly, return_polys) where T
    n_new_polys = 0
    remove_idx = Int[]
    new_hole_poly = GI.Polygon([new_hole])
    # Combine any existing holes in curr_poly with new hole
    for (k, old_hole) in enumerate(GI.gethole(curr_poly))
        old_hole_poly = GI.Polygon([old_hole])
        if intersects(new_hole_poly, old_hole_poly)
            # If the holes intersect, combine them into a bigger hole
            hole_union = union(new_hole_poly, old_hole_poly, T; target = GI.PolygonTrait)[1]
            push!(remove_idx, k + 1)
            new_hole = GI.getexterior(hole_union)
            new_hole_poly = GI.Polygon([new_hole])
            n_pieces = GI.nhole(hole_union)
            if n_pieces > 0  # if the hole has a hole, then this is a new polygon piece! 
                append!(return_polys, [GI.Polygon([h]) for h in GI.gethole(hole_union)])
                n_new_polys += n_pieces
            end
        end
    end
    # Remove redundant holes
    deleteat!(curr_poly.geom, remove_idx)
    empty!(remove_idx)
    # If new polygon pieces created, make sure remaining holes are in the correct piece
    @views for piece in return_polys[end - n_new_polys + 1:end]
        for (k, old_hole) in enumerate(GI.gethole(curr_poly))
            if !(k in remove_idx) && within(old_hole, piece)
                push!(remove_idx, k + 1)
                push!(piece.geom, old_hole)
            end
        end
    end
    deleteat!(curr_poly.geom, remove_idx)
    return new_hole, new_hole_poly, n_new_polys
end
