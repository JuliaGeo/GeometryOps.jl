# #  Union Polygon Clipping
export union

# This 'union' implementation returns the union of two polygons.
# It returns a Vector{Vector{Vector{Tuple{Float}}}. Note that this file only contains
# the functionality of a union of two polygons, and not other geometries.

"""
    _trace_union(::GI.PolygonTrait, poly_a,
     ::GI.PolygonTrait, poly_b, a_list, b_list, a_idx_list,
      intr_list, edges_a, edges_b)::Vector{Vector{Tuple{Float64}}}

Traces the outlines of two polygons in order to find their union.
It returns the outlines of all polygons formed in the union. If
one polygon is completely contained in the other, it returns
the larger one.

## Example

TODO
"""

function _trace_union(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
    # Pre-allocate array for return polygons
    return_polys = Vector{Vector{Tuple{Float64, Float64}}}(undef, 0)
    # Keep track of number of processed intersection points
    processed_pts = 0
    tracker = copy(a_idx_list)

    while processed_pts < length(intr_list)
        # Create variables "list_edges" and "list" so that we can toggle between
        # a_list and b_list
        list_edges = edges_a
        list = a_list

        # Find index of first unprocessed intersecting point in subject polygon
        starting_pt = minimum(tracker)
        idx = starting_pt

        # Get current first unprocessed intersection point PolyNode
        current = a_list[idx]
        # Initialize array to store the intersection polygon cartesian points
        pt_list = Vector{Tuple{Float64, Float64}}(undef, 0)
        # Add the first point to the array
        push!(pt_list, (intr_list[current.idx][1], intr_list[current.idx][2]))
        
        # Mark first intersection point as processed
        processed_pts = processed_pts + 1
        tracker[current.idx] = typemax(Int)

        current_node_not_starting = true
        while current_node_not_starting # While the current node isn't the starting one
            status2 = false
            current_node_not_intersection = true
            while current_node_not_intersection # The current node isn't an intersection
                
                if current.inter
                    status2 = current.ent_exit
                end

                # Depending on status of first intersection point, either
                # traverse polygon forwards or backwards
                if status2
                    idx = idx - 1
                else
                    idx = idx + 1
                end

                # Wrap around the point list
                if idx > length(list)
                    idx = mod(idx, length(list))
                elseif idx == 0
                    idx = length(list)
                end

                # Get current node
                current = list[idx]

                # Add current node to the pt_list
                if current.inter
                    # Add cartesian coordinates from inter_list
                    push!(pt_list, (intr_list[current.idx][1], intr_list[current.idx][2]))
                    
                    # Keep track of processed intersection points
                    if (current != a_list[starting_pt] && current != b_list[a_list[starting_pt].neighbor])
                        processed_pts = processed_pts + 1
                        tracker[current.idx] = typemax(Int)
                    end
                    
                else
                    # Add cartesian coordinates from "list", which should point to either a_list or b_list
                    push!(pt_list, (list_edges[current.idx][1][1], list_edges[current.idx][1][2]))
                end

                current_node_not_intersection = !current.inter
            end
            
            # Break once get back to starting point
            current_node_not_starting = (current != a_list[starting_pt] && current != b_list[a_list[starting_pt].neighbor])

            # Switch to neighbor list
            if list == a_list
                list = b_list
                list_edges = edges_b
            else
                list = a_list
                list_edges = edges_a
            end
            idx = current.neighbor
            current = list[idx]
        end

        push!(return_polys, pt_list)
    end

    # Check if one polygon totally within other and if so, return the larger polygon.
    # TODO: use point list instead of edges once to_points is fixed
    if isempty(return_polys)
        if point_in_polygon(edges_a[1][1], poly_b)[1]
            list = []
            for i in eachindex(edges_b)
                push!(list, edges_b[i][1])
            end
            push!(list, edges_b[1][1])
            push!(return_polys, list)
            return return_polys, false
        elseif point_in_polygon(edges_b[1][1], poly_a)[1]
            list = []
            for i in eachindex(edges_a)
                push!(list, edges_a[i][1])
            end
            push!(list, edges_a[1][1])
            push!(return_polys, list)
            return return_polys, false
        else
            # In the case that the polygons don't intersect and aren't contained in
            # one another, return both polygons.
            push!(return_polys, to_points(poly_a))
            push!(return_polys, to_points(poly_b))
            return return_polys, true
        end
    end

    # It is possible that at this point, we have multiple polygons
    # in 'return_polys'. One of those polygons is the outermost polygon,
    # and the rest are it's hole. In the code below, I find the outermost
    # polygon, make sure it is first in return_polys, and that all it's
    # holes come after it.
    if length(return_polys) > 1
        # Find the index in 'return_polys' of the outermost polygon.
        outer_idx = 1
        for j = 2:(length(return_polys)-1)
            poly1 = GI.Polygon([return_polys[outer_idx]])
            poly2 = GI.Polygon([return_polys[j]])
            if !point_in_polygon(to_edges(poly2)[1][1], poly1)
                outer_idx = j
            end
        end

        # Reorder 'return_polys'
        if outer_idx == 1
            return return_polys, false
        elseif outer_idx == length(return_polys)
            arr = outer_idx
            splice!(arr, 2:1, 1:(length(return_polys)-1))
            return return_polys[arr], false
        else
            arr = 1:(outer_idx-1)
            splice!(arr, (length(arr)+1):length(arr), (outer_idx+1):length(return_polys))
            return return_polys[arr], false
        end
    end

    return return_polys, false
end

"""
    _get_union_holes(return_polys, poly_a, poly_b, ext_poly_b)::Vector{Vector{Vector{Tuple{Float64, Float64}}}}

When the _union_difference function was called, it only took into account the
exteriors of the two polygons when computing the union. The function
'_get_union_holes' takes into account the holes of the original polygons
and adjusts the output of _union_difference (return_polys) accordingly.

## Example

TODO
```
"""

function _get_union_holes(return_polys, poly_a, poly_b, ext_poly_b, diff_polys)
    # First I initiliaze the return object.
    # This will be of length one, but I'm keeping it this type to be consistent
    # with the other polygon clipping functions (intersection and differece).
    final_polys =  Vector{Vector{Vector{Tuple{Float64, Float64}}}}(undef, 0)
    
    if !diff_polys
        # I set mult_poly equal to return_polys because I know that in this case,
        # return_polys only represents one polygon.
        mult_poly = return_polys
        
        for hole in GI.gethole(poly_a)
            # I use ext_poly_b here instead of poly_b in order to not overcount
            # area that is a hole of poly_a and in a hole of poly_b
            new_hole = difference(_lin_ring_to_poly(hole), ext_poly_b)
            for h in new_hole
                if length(h)>0
                    # I claim I can index h at one, because it can't
                    # have a hole because a hole within a hole would
                    # be an invalid polygon
                    append!(mult_poly, [h[1]])
                end
            end        
        end
        
        for hole in GI.gethole(poly_b)
            new_hole = difference(_lin_ring_to_poly(hole), poly_a)
            for h in new_hole
                if length(h)>0
                    append!(mult_poly, [h[1]])
                end
            end 
        end
        push!(final_polys, mult_poly)

        return final_polys
    else
        # If 'return_polys' represents multiple polygons and we are in this case,
        # poly_a and poly_b are completely disjoint. So we need to subtract holes
        # from them separately.
        @assert (length(return_polys)==2) "Since 'diff_polys is true, 'return_polys' should have a length of 2. Instead it has length of $(length(return_polys))."
        for poly in return_polys
            mult_poly = [poly]
            for hole in GI.gethole(poly_a)
                # I use ext_poly_b here instead of poly_b in order to not overcount
                # area that is a hole of poly_a and in a hole of poly_b
                new_hole = difference(_lin_ring_to_poly(hole), ext_poly_b)
                for h in new_hole
                    if length(h)>0
                        # I claim I can index h at one, because it can't
                        # have a hole because a hole within a hole would
                        # be an invalid polygon
                        append!(mult_poly, [h[1]])
                    end
                end        
            end
            
            for hole in GI.gethole(poly_b)
                new_hole = difference(_lin_ring_to_poly(hole), poly_a)
                for h in new_hole
                    if length(h)>0
                        append!(mult_poly, [h[1]])
                    end
                end 
            end

            push!(final_polys, mult_poly)
        end

        return final_polys
    end
        
end


"""
    union(
        ::GI.PolygonTrait, poly_a,
        ::GI.PolygonTrait, poly_b,
    )::Vector{Vector{Vector{Tuple{Float64}}}}    

Calculates the union between two polygons.
## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0), (5.0, -5.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(3.0, 0.0), (8.0, 5.0), (13.0, 0.0), (8.0, -5.0), (3.0, 0.0)]])
GO.union(p1, p2)

# output
1-element Vector{Vector{Vector{Tuple{Float64, Float64}}}}:
[[[(6.5, 3.5), (5.0, 5.0), (0.0, 0.0), (5.0, -5.0), (6.5, -3.5), (8.0, -5.0), (13.0, 0.0), (8.0, 5.0), (6.5, 3.5)]]]
```
"""
function union(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    # First, I get the exteriors of the two polygons
    ext_poly_a = GI.getexterior(poly_a)
    ext_poly_a = _lin_ring_to_poly(ext_poly_a)
    ext_poly_b = GI.getexterior(poly_b)
    ext_poly_b = _lin_ring_to_poly(ext_poly_b)
    # Then, I get the union of the exteriors
    a_list, b_list, a_idx_list, intr_list, edges_a, edges_b = _build_ab_list(GI.trait(ext_poly_a), ext_poly_a, GI.trait(ext_poly_b), ext_poly_b)
    temp = _trace_union(GI.trait(ext_poly_a), ext_poly_a, GI.trait(ext_poly_b), ext_poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
    polys = temp[1]
    diff_polys = temp[2]
    # If the original polygons had holes, we call '_get_union_holes' to take that
    # into account.
    if GI.nhole(poly_a)==0 && GI.nhole(poly_b)==0
        if !diff_polys
            return [polys]
        else
            final_polys =  Vector{Vector{Vector{Tuple{Float64, Float64}}}}(undef, length(polys))
            for i in 1:length(polys)
                final_polys[i] = [polys[i]]
            end
            return final_polys
        end
    else
        return _get_union_holes(polys, poly_a, poly_b, ext_poly_b, diff_polys)
    end 
end

union(geom_a, geom_b) =
    union(GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)