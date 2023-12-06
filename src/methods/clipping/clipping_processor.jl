# # This file contains the shared helper functions for the other polygon clipping
# # functionalities.

"""
    _lin_ring_to_poly(lin_ring)::GI.Polygon

    This function turns a linear ring into a GeometryInterface polygon.
"""
function _lin_ring_to_poly(lin_ring)
    edges = to_edges(lin_ring)
    points = Vector{Tuple{Float64, Float64}}(undef, length(edges)+1)
    for ii in eachindex(edges)
        points[ii] = edges[ii][1]
    end
    points[end] = edges[1][1]
    return GI.Polygon([points])
end

"""
    _build_a_list(
        ::GI.PolygonTrait, 
        poly_a, ::GI.PolygonTrait, poly_b,
         edges_a, edges_b, intr_list, a_idx_list, 
         b_idx_list, alpha_a_list, alpha_b_list, 
         a_list, b_list
    )::intr_list, a_idx_list, b_idx_list, alpha_a_list, alpha_b_list, counter

    This function take in two polygons, the lists of their edges, and some lists containing
    information about the nature of their intersection points. From that, it creates an 
    array to represent polygon_a, and this array contains intersection points. However, after
    calling this function, a_list is not fully formed because the neighboring indicies of the
    intersection points in b_list still need to be updated. We will have fully formed a_list
    after calling _build_ab_list. Also at this point we still have not update the entry and
    exit flags for a_list.
"""
function _build_a_list(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b, edges_a, edges_b,
                      intr_list, a_idx_list, b_idx_list, alpha_a_list, alpha_b_list, a_list, b_list)
    # Find intersection points and adds them to a_list
    # "counter" is used to index all inter-related lists
    counter = 1
    # "acount" is used to index a_list
    acount = 1
    for ii in eachindex(edges_a)
        # Add the first point of the edge to the list of points in a_list
        if acount <= length(a_list)
            a_list[acount] = PolyNode(ii, false, 0, false, 0)
        else
            push!(a_list, PolyNode(ii, false, 0, false, 0))
        end
        acount = acount + 1

        # Keep track of current position in inter-related lists
        # before finding new intersection points on our new edge
        prev_counter = counter

        for jj in eachindex(edges_b)

            # Add the first point of the edge to b_list
            if ii == 1 
                b_list[jj] = PolyNode(jj, false, 0, false, 0)
            end

            # Check if edge jj of poly_b intersects with edge ii of poly_a
            if _line_intersects([edges_a[ii]], [edges_b[jj]]);
                
                int_pt, alphas = _intersection_point(edges_a[ii], edges_b[jj])
                # if not intersection point, skip this edge (the if statement above should 
                # catch this but i guess not)
                if isnothing(int_pt)
                    continue
                end
                # Store the cartesion coordinates of intersection point
                intr_list[counter] = int_pt
                # Store which edge of poly_b the intersection point lies on
                b_idx_list[counter] = jj
                # Store the alpha values
                alpha_a_list[counter] = alphas[1]
                alpha_b_list[counter] = alphas[2]

                counter = counter + 1
            else
                continue
            end
        end

        # After iterating through all edges of poly_b for edge ii of poly_a,
        # add the intersection points to a_list in CORRECT ORDER if we found any
        if prev_counter < counter
            # Order intersection points based on alpha values
            new_order = sortperm(alpha_a_list[prev_counter:counter-1])
            pts_to_add = Array{PolyNode, 1}(undef, counter - prev_counter)
            for kk in eachindex(new_order)
                # Create PolyNodes of the new intersection points in the correct order
                # and store the correct index in a_idx_list
                pts_to_add[new_order[kk]] = PolyNode(prev_counter+kk-1, true, 0, false, alpha_a_list[prev_counter+kk-1])
                a_idx_list[prev_counter+kk-1] = acount + new_order[kk] - 1;
            end

            # Add the PolyNodes to a_list and update acount
            splice!(a_list, acount:(acount-1), pts_to_add)
            acount = acount + counter - prev_counter
            
        end
        
    end

    # Truncate the inter_related lists
    intr_list = intr_list[1:counter-1] 
    a_idx_list = a_idx_list[1:counter-1]
    b_idx_list = b_idx_list[1:counter-1]
    alpha_a_list = alpha_a_list[1:counter-1]
    alpha_b_list = alpha_b_list[1:counter-1]
    return intr_list, a_idx_list, b_idx_list, alpha_a_list, alpha_b_list, counter
end

"""
    _build_b_list(
        b_list, intr_list, a_idx_list, 
        b_idx_list, alpha_b_list, counter, k
    )::a_list, b_list

    This function builds b_list. Note that after calling this function, b_list
    is not fully updated. The entry/exit flags still need to be updated.
"""

function _build_b_list(b_list, intr_list, a_idx_list, b_idx_list, alpha_b_list, counter, k)
    # Iterate through the b_list and add in intersection points
    # Occasionally I need to skip the new points I added to the array
    skip = false
    num_skips = 0
    b_neighbors = Array{Int, 1}(undef, k)
    for ii in 1:(length(b_list)+length(intr_list))
        # TODO: this skipping scheme could be made nicer
        if skip
            num_skips = num_skips - 1
            if num_skips == 0
                skip = false
            end
            continue
        end
        # Find the index in the intr_list (same as b_idx_list) where the inter point is
        # TODO:find all is inefficient though, so it might be better to make dictionary of indices
        i = findall(x->x==b_list[ii].idx, b_idx_list)
        if !isempty(i)     
            # Order intersection points based on alpha values
            new_order = sortperm(alpha_b_list[i])
            pts_to_add = Array{PolyNode, 1}(undef, length(i))
            for m in eachindex(i)
                pts_to_add[new_order[m]] = PolyNode(i[m], true, a_idx_list[i[m]], false, alpha_b_list[i[m]])
                b_neighbors[i[m]] = ii + new_order[m]
            end   
            # I use splice instead of insert so I can insert array   
            splice!(b_list, ii+1:ii, pts_to_add)
            skip = true
            num_skips = length(i)
        end
    end

    b_neighbors = b_neighbors[1:counter-1]

    return b_neighbors, b_list
end

"""
    _flag_ent_exit(
        ::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, 
        poly_b, a_list, b_list, edges_a, edges_b
    )::a_list, b_list

    This function flags all the intersection points as either an 'entry' or 'exit' point.
"""
function _flag_ent_exit(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b, a_list, b_list, edges_a, edges_b)
    # Put in ent exit flags for poly_a
    status = false
    for ii in eachindex(a_list)
        if ii == 1
            temp = point_in_polygon(edges_a[ii][1], poly_b)
            status = !(temp[1])
            continue
        end
        if a_list[ii].inter
            a_list[ii].ent_exit = status
            status = !status
        end
    end

    # Put in ent exit flags for poly_b
    status = false
    for ii in eachindex(b_list)
        if ii == 1
            temp = point_in_polygon(edges_b[ii][1], poly_a)
            status = !(temp[1])
            continue
        end
        if b_list[ii].inter
            b_list[ii].ent_exit = status
            status = !status
        end
    end

    return a_list, b_list
end

"""
    _build_ab_list(
        ::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b
    )::a_list, b_list, a_idx_list, intr_list, edges_a, edges_b

    This function calls '_build_a_list', '_build_b_list', and '_flag_ent_exit'
    in order to fully form a_list and b_list.
"""
function _build_ab_list(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    # Make a list for nodes of each polygon. Note the definition of PolyNode
    a_list = Array{PolyNode, 1}(undef, _nedge(poly_a))
    b_list = Array{PolyNode, 1}(undef, _nedge(poly_b))

    # Initialize arrays to keep track of the important information 
    # associated with the intersection points of poly_a and poly_b.
    # I initialize these arrays assumed there will be a maximum of 
    # 30 intersection points and then I truncate the arrays later.
    k = 4
    # intr_list stores the cartesian coordintes of the intersection point
    intr_list = Array{Tuple{Real, Real}, 1}(undef, k)
    # At index i, a_idx_list stores the index of the ith intersection point in a_list
    a_idx_list = Array{Int, 1}(undef, k)
    # At index i, b_idx_list stores the edge number of poly_b that the ith intersection
    # point lies on.
    b_idx_list = Array{Int, 1}(undef, k)
    # Alpha values are used to determine the order in which to place
    # intersection points in a polygon list, which is especially
    # useful when many intersection points lie on the same edge of a polygon.
    alpha_a_list = Array{Real, 1}(undef, k)
    alpha_b_list = Array{Real, 1}(undef, k)

    # These lists store the cartesian coordinates of poly_a and poly_b
    edges_a = to_edges(poly_a)
    edges_b = to_edges(poly_b)

    intr_list, a_idx_list, b_idx_list, alpha_a_list, alpha_b_list, counter = _build_a_list(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b, edges_a, edges_b,
                                                                                  intr_list, a_idx_list, b_idx_list, alpha_a_list, alpha_b_list, a_list, b_list)
    b_neighbors, b_list = _build_b_list(b_list, intr_list, a_idx_list, b_idx_list, alpha_b_list, counter, k)

    
    # Iterate through a_list and update the neighbor indices
    for ii in eachindex(a_list)
        if a_list[ii].inter
            a_list[ii].neighbor = b_neighbors[a_list[ii].idx]
        end
    end

    # Flag the entry and exists
    a_list, b_list = _flag_ent_exit(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b, a_list, b_list, edges_a, edges_b)

    return a_list, b_list, a_idx_list, intr_list, edges_a, edges_b
end

# This is the struct that makes up a_list and b_list
mutable struct PolyNode
    idx::Int
    inter::Bool
    neighbor::Int
    ent_exit::Bool
    alpha::Real
end