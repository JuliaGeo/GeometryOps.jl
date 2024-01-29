# # This file contains the shared helper functions for the polygon clipping functionalities.

#= This is the struct that makes up a_list and b_list. Many values are only used if point is
an intersection point (ipt). =#
mutable struct PolyNode{T <: AbstractFloat}
    idx::Int           # If ipt, index of point in a_idx_list, else 0
    point::Tuple{T,T}  # (x, y) values of given point
    inter::Bool        # If ipt, true, else 0
    neighbor::Int      # If ipt, index of equivalent point in a_list or b_list, else 0
    ent_exit::Bool     # If ipt, true if enter and false if exit, else false
    fracs::Tuple{T,T}  # If ipt, fractions along edges to ipt (a_frac, b_frac), else (0, 0)
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
    a_list, a_idx_list = _build_a_list(T, poly_a, poly_b)
    b_list = _build_b_list(T, a_idx_list, a_list, poly_b)

    # Flag the entry and exits
    _flag_ent_exit!(poly_b, a_list)
    _flag_ent_exit!(poly_a, b_list)

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
    intr_count = 0  # number of intersection points found
    a_count = 0  # number of points added to a_list
    # Loop through points of poly_a
    local a_pt1
    for (i, a_p2) in enumerate(GI.getpoint(poly_a))
        a_pt2 = (T(GI.x(a_p2)), T(GI.y(a_p2)))
        if i <= 1
            continue
        end
        # Add the first point of the edge to the list of points in a_list
        new_point = PolyNode(0, a_pt1, false, 0, false, (zero(T), zero(T)))
        a_count += 1
        _add!(a_list, a_count, new_point, n_a_edges)
        # Find intersections with edges of poly_b
        local b_pt1
        prev_counter = intr_count
        for (j, b_p2) in enumerate(GI.getpoint(poly_b))
            b_pt2 = _tuple_point(b_p2)
            if j <=1
                continue
            end
            int_pt, fracs = _intersection_point(T, (a_pt1, a_pt2), (b_pt1, b_pt2))
            # if no intersection point, skip this edge
            if !isnothing(int_pt) && all(0 .≤ fracs .≤ 1)
                # Set neighbor field to b edge (j-1) to keep track of intersection
                new_intr = PolyNode(intr_count, int_pt, true, j - 1, false, fracs)
                a_count += 1
                intr_count += 1
                _add!(a_list, a_count, new_intr, n_a_edges)
                push!(a_idx_list, a_count)
            end
            b_pt1 = b_pt2
        end

        # Order intersection points by placement along edge using fracs value
        if prev_counter < intr_count
            Δintrs = intr_count - prev_counter
            inter_points = @view a_list[(a_count - Δintrs + 1):a_count]
            sort!(inter_points, by = x -> x.fracs[1])
            for (i, p) in enumerate(inter_points)
                p.idx = prev_counter + i
            end
        end
    
        a_pt1 = a_pt2
    end
    return a_list, a_idx_list
end

# Add value x at index i to given array - if list isn't long enough, push value to array
function _add!(arr, i, x, l = length(arr))
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
function _build_b_list(::Type{T}, a_idx_list, a_list, poly_b) where T
    # Sort intersection points by insertion order in b_list
    sort!(a_idx_list, by = x-> a_list[x].neighbor + a_list[x].fracs[2])
    # Initialize needed values and lists
    n_b_edges = _nedge(poly_b)
    n_intr_pts = length(a_idx_list)
    b_list = Vector{PolyNode{T}}(undef, n_b_edges + n_intr_pts)
    intr_curr = 1
    b_count = 0
    # Loop over points in poly_b and add each point and intersection point
    for (i, p) in enumerate(GI.getpoint(poly_b))
        (i == n_b_edges + 1) && break
        b_count += 1
        pt = (T(GI.x(p)), T(GI.y(p)))
        b_list[b_count] = PolyNode(0, pt, false, 0, false, (zero(T), zero(T)))
        if intr_curr ≤ n_intr_pts
            curr_idx = a_idx_list[intr_curr]
            curr_node = a_list[curr_idx]
            while curr_node.neighbor == i  # Add all intersection points in current edge
                b_count += 1
                b_list[b_count] = PolyNode(curr_node.idx, curr_node.point, true, curr_idx, false, curr_node.fracs)
                curr_node.neighbor = b_count
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
    _flag_ent_exit(poly_b, a_list)

This function flags all the intersection points as either an 'entry' or 'exit' point in
relation to the given polygon.
=#
function _flag_ent_exit!(poly, pt_list)
    local status
    for ii in eachindex(pt_list)
        if ii == 1
            status = !_point_filled_curve_orientation(
                pt_list[ii].point, poly;
                in = true, on = false, out = false
            )
        elseif pt_list[ii].inter
            pt_list[ii].ent_exit = status
            status = !status
        end
    end
    return
end

#=
    _trace_polynodes(a_list, b_list, a_idx_list, f_step)::Vector{GI.Polygon}

This function takes the outputs of _build_ab_list and traces the lists to determine which
polygons are formed as described in Greiner and Hormann. The function f_step determines in
which direction the lists are traced.  This function is different for intersection,
difference, and union. f_step must take in two arguments: the most recent intersection
node's entry/exit status and a boolean that is true if we are currently tracing a_list and
false if we are tracing b_list. The functions used for each clipping operation are follows:
    - Intersection: (x, y) -> x ? 1 : (-1)
    - Difference: (x, y) -> (x ⊻ y) ? 1 : (-1)
    - Union: (x, y) -> (x ⊻ y) ? 1 : (-1)

A list of GeoInterface polygons is returned from this function. 
=#
function _trace_polynodes(a_list, b_list, a_idx_list, f_step)
    n_a_pts, n_b_pts = length(a_list), length(b_list)
    n_intr_pts = length(a_idx_list)
    return_polys = Vector{GI.Polygon}(undef, 0)
    # Keep track of number of processed intersection points
    processed_pts = 0

    while processed_pts < n_intr_pts
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
            curr_not_intr = true
            while curr_not_intr
                # Traverse polygon either forwards or backwards
                idx += step
                idx = (idx > curr_npoints) ? mod(idx, curr_npoints) : idx
                idx = (idx == 0) ? curr_npoints : idx

                # Get current node and add to pt_list
                curr = curr_list[idx]
                push!(pt_list, curr.point)
                if curr.inter 
                    # Keep track of processed intersection points
                    curr_not_start = curr != start_pt && curr != b_list[start_pt.neighbor]
                    if curr_not_start
                        processed_pts += 1
                        a_idx_list[curr.idx] = 0
                    end
                    curr_not_intr = false
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
    return return_polys::Vector{GI.Polygon}
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
        for hole in hole_iterator  # loop through all holes
            hole_poly = GI.Polygon([hole])
            # loop through all pieces of original polygon (new pieces added to end of list)
            for j in Iterators.flatten((i:i, (n_polys + 1):(n_polys + n_new_per_poly)))
                if !isnothing(return_polys[j])
                    new_polys = difference(return_polys[j], hole_poly, T; target = GI.PolygonTrait)
                    n_new_polys = length(new_polys)
                    if n_new_polys == 0  # hole covered whole polygon
                        return_polys[j] = nothing
                    else
                        return_polys[j] = new_polys[1]  # replace original
                        if n_new_polys > 1  # add any extra pieces
                            append!(return_polys, @view new_polys[2:end])
                            n_new_per_poly += n_new_polys - 1
                        end
                    end
                end
            end
        end
        n_polys += n_new_per_poly
    end
    # Remove all polygon that were marked for removal
    filter!(!isnothing, return_polys)::Vector{GI.Polygon}
    return
end