# # This file contains the shared helper functions for the other polygon clipping
# # functionalities.

"""
    _build_a_list(poly_a, poly_b)::a_list, a_idx_list

    This function take in two polygons, the lists of their edges, and some lists containing
    information about the nature of their intersection points. From that, it creates an 
    array to represent polygon_a, and this array contains intersection points. However, after
    calling this function, a_list is not fully formed because the neighboring indicies of the
    intersection points in b_list still need to be updated. We will have fully formed a_list
    after calling _build_ab_list. Also at this point we still have not update the entry and
    exit flags for a_list. The variables poly_a and poly_b are the original polygons that the
    user passed in. The variables edges_a and edges_b are lists of the edges of each polygon, 
    in the form of their cartesian points. Then a_idx_list is a list of the indicies of where in a_list
    an intersection point lies. The value at index i of a_idx_list is the location in a_list where
    the ith intersection point lies. Then alpha_a_list and alpha_b_list are vectors of scalar values which help order
    the intersection points in the case that they lie on the same edge of a polygon. Finally,
    a_list is an empty array of PolyNode objects that gets filled out when this function is run.
    So is b_list. 
"""
function _build_a_list(poly_a, poly_b)
    a_list = Vector{PolyNode}(undef, _nedge(poly_a))
    # Find intersection points and adds them to a_list
    a_idx_list = Vector{Int}()
    # "counter" is used to index all inter-related lists
    counter = 0
    # "acount" is used to index a_list
    acount = 1

    start = true
    local p1
    ii = 0
    for (i, p2) in enumerate(GI.getpoint(poly_a))
        if start
            start = false
        else
            new_point =  PolyNode(i - 1, _tuple_point(p1), false, 0, false, (0.0, 0.0))
            # Add the first point of the edge to the list of points in a_list
            if acount <= length(a_list)
                a_list[acount] = new_point
            else
                push!(a_list, new_point)
            end
            acount = acount + 1

            jj = 0
            start2 = true
            local g1
            prev_counter = counter
            for (j, g2) in enumerate(GI.getpoint(poly_b))
                if start2
                    start2 = false
                else
                    # Check if edge jj of poly_b intersects with edge ii of poly_a
                    if _line_intersects([(_tuple_point(p1), _tuple_point(p2))], [(_tuple_point(g1), _tuple_point(g2))]);
                        int_pt, fracs = _intersection_point((_tuple_point(p1), _tuple_point(p2)), (_tuple_point(g1), _tuple_point(g2)))
                        # if not intersection point, skip this edge (the if statement above should 
                        # catch this but i guess not)
                        isnothing(int_pt) && continue
                        new_intr = PolyNode(counter, int_pt, true, j - 1, false, fracs)
                        if acount <= length(a_list)
                            a_list[acount] = new_intr
                        else
                            push!(a_list, new_intr)
                        end
                        push!(a_idx_list, acount)
                        counter += 1
                        acount += 1
                    end
                end
                jj = jj + 1
                g1 = g2
            end

            # After iterating through all edges of poly_b for edge ii of poly_a,
            # add the intersection points to a_list in CORRECT ORDER if we found any
            if prev_counter < counter
                inter_points = @view a_list[(acount - counter + prev_counter):(acount - 1)]
                sort!(inter_points, by = x -> x.fracs[1])
                for (i, p) in enumerate(inter_points)
                    p.idx = prev_counter + i
                end
            end
        end
        p1 = p2
        ii = ii + 1
    end
    return a_list, a_idx_list
end

"""
    _build_b_list(
        b_list, 
        a_idx_list, 
        alpha_b_list
    )::b_neighbors, b_list

    This function builds b_list. Note that after calling this function, b_list
    is not fully updated. The entry/exit flags still need to be updated. The variables poly_a and poly_b are the original polygons that the
    user passed in. The variables edges_a and edges_b are lists of the edges of each polygon, 
    in the form of their cartesian points. Then a_idx_list is a list of the indicies of where in a_list
    an intersection point lies. The value at index i of a_idx_list is the location in a_list where
    the ith intersection point lies. Then alpha_b_list is a veector of scalar values which help order
    the intersection points in the case that they lie on the same edge of a polygon. Finally,
    b_list starts out as an array of PolyNodes only containing the original points of poly_b
    but after this function is run it include intersection points to.
"""
function _build_b_list(a_idx_list, a_list, poly_b)
    sort_a_idx_list = sort(a_idx_list, by = x-> a_list[x].neighbor+a_list[x].fracs[2])
    b_list = Array{PolyNode, 1}(undef, _nedge(poly_b)+length(a_idx_list))
    counter = 1
    bcounter = 1
    for (i, pi) in enumerate(GI.getpoint(poly_b))
        i == GI.npoint(poly_b) && break

        b_list[bcounter] = PolyNode(i, _tuple_point(pi), false, 0, false, (0.0, 0.0))
        bcounter += 1
        if counter <= length(sort_a_idx_list)
            current_node = a_list[sort_a_idx_list[counter]]
            while current_node.neighbor == i
                b_list[bcounter] = PolyNode(counter, current_node.point, true, sort_a_idx_list[counter], false, current_node.fracs)
                current_node.neighbor = bcounter
                current_node.idx = counter
                bcounter += 1
                counter += 1
                counter>length(sort_a_idx_list) && break
                current_node = a_list[sort_a_idx_list[counter]]
            end
        end
    end
    return b_list, sort_a_idx_list
end

"""
    _flag_ent_exit(poly_b, a_list)::a_list

    This function flags all the intersection points as either an 'entry' or 'exit' point.
"""
function _flag_ent_exit(poly_b, a_list)
    # Put in ent exit flags for poly_a
    status = false
    for ii in eachindex(a_list)
        if ii == 1
            temp = within(a_list[ii].point, poly_b)
            status = !(temp)
            continue
        end
        if a_list[ii].inter
            a_list[ii].ent_exit = status
            status = !status
        end
    end
    return a_list
end

"""
    _build_ab_list(
        poly_a, poly_b
    )::a_list, b_list, a_idx_list

    This function calls '_build_a_list', '_build_b_list', and '_flag_ent_exit'
    in order to fully form a_list and b_list. The 'a_list' and 'b_list'
    that it returns are the fully updated vectors of PolyNodes that represent
    'poly_a' and 'poly_b', respectively. This function also returns
    'a_idx_list', which at its "ith" index stores the index in 'a_list' at
    which the "ith" intersection point lies.
"""
function _build_ab_list(poly_a, poly_b)
    # Make a list for nodes of each polygon
    a_list, a_idx_list = _build_a_list(poly_a, poly_b)
    b_list, sort_a_idx_list = _build_b_list(a_idx_list, a_list, poly_b)

    # Flag the entry and exists
    a_list = _flag_ent_exit(poly_b, a_list)
    b_list = _flag_ent_exit(poly_a, b_list)

    return a_list, b_list, sort_a_idx_list
end

# This is the struct that makes up a_list and b_list
mutable struct PolyNode{T <: AbstractFloat}
    idx::Int
    point::Tuple{T,T}
    inter::Bool
    neighbor::Int
    ent_exit::Bool
    fracs::Tuple{T,T}
end