# # Intersection checks

export intersects, intersection, intersection_points, union_test, difference_test

#=
## What is `intersects` vs `intersection` vs `intersection_points`?

The `intersects` methods check whether two geometries intersect with each other.
The `intersection` methods return the geometry intersection between the two
input geometries. The `intersection_points` method returns a list of
intersection points between two geometries.

The `intersects` methods will always return a Boolean. However, note that the
`intersection` methods will not all return the same type. For example, the
intersection of two lines will be a point in most cases, unless the lines are
parallel. On the other hand, the intersection of two polygons will be another
polygon in most cases. Finally, the `intersection_points` method returns a list
of tuple points.

To provide an example, consider these two lines:
```@example intersects_intersection
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie
point1, point2 = Point(124.584961,-12.768946), Point(126.738281,-17.224758)
point3, point4 = Point(123.354492,-15.961329), Point(127.22168,-14.008696)
line1 = Line(point1, point2)
line2 = Line(point3, point4)
f, a, p = lines([point1, point2])
lines!([point3, point4])
```
We can see that they intersect, so we expect intersects to return true, and we
can visualize the intersection point in red.
```@example intersects_intersection
int_bool = GO.intersects(line1, line2)
println(int_bool)
int_point = GO.intersection(line1, line2)
scatter!(int_point, color = :red)
f
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method for intersects, intersection, and
intersection_points that dispatches to the correct implementation based on the
geometry trait. The two underlying helper functions that are widely used in all
geometry dispatches are _line_intersects, which determines if two line segments
intersect and _intersection_point which determines the intersection point
between two line segments.
=#


# tempory poly_in_poly, only works if entirely contained
# Checks if polya in poly b
function poly_in_poly(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    p_a = to_edges(poly_a)[1][1]
    return point_in_polygon(p_a, poly_b)
end

poly_in_poly(geom_a, geom_b) =
    poly_in_poly(GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)

"""
    intersects(geom1, geom2)::Bool

Check if two geometries intersect, returning true if so and false otherwise.

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
GO.intersects(line1, line2)

# output
true
```
"""
intersects(geom1, geom2) = intersects(
    GI.trait(geom1),
    geom1,
    GI.trait(geom2),
    geom2
)

"""
    intersects(::GI.LineTrait, a, ::GI.LineTrait, b)::Bool

Returns true if two line segments intersect and false otherwise.
"""
function intersects(::GI.LineTrait, a, ::GI.LineTrait, b)
    a1 = _tuple_point(GI.getpoint(a, 1))
    a2 = _tuple_point(GI.getpoint(a, 2))
    b1 = _tuple_point(GI.getpoint(b, 1))
    b2 = _tuple_point(GI.getpoint(b, 2))
    meet_type = ExactPredicates.meet(a1, a2, b1, b2)
    return meet_type == 0 || meet_type == 1
end

"""
    intersects(::GI.AbstractTrait, a, ::GI.AbstractTrait, b)::Bool

Returns true if two geometries intersect with one another and false
otherwise. For all geometries but lines, convert the geometry to a list of edges
and cross compare the edges for intersections.
"""
function intersects(
    trait_a::GI.AbstractTrait, a_geom,
    trait_b::GI.AbstractTrait, b_geom,
)   edges_a, edges_b = map(sort! ∘ to_edges, (a_geom, b_geom))
    return _line_intersects(edges_a, edges_b) ||
        within(trait_a, a_geom, trait_b, b_geom) ||
        within(trait_b, b_geom, trait_a, a_geom) 
end

"""
    _line_intersects(
        edges_a::Vector{Edge},
        edges_b::Vector{Edge}
    )::Bool

Returns true if there is at least one intersection between edges within the
two lists of edges.
"""
function _line_intersects(
    edges_a::Vector{Edge},
    edges_b::Vector{Edge}
)
    # Extents.intersects(to_extent(edges_a), to_extent(edges_b)) || return false
    for edge_a in edges_a
        for edge_b in edges_b
            _line_intersects(edge_a, edge_b) && return true 
        end
    end
    return false
end

"""
    _line_intersects(
        edge_a::Edge,
        edge_b::Edge,
    )::Bool

Returns true if there is at least one intersection between two edges.
"""
function _line_intersects(edge_a::Edge, edge_b::Edge)
    meet_type = ExactPredicates.meet(edge_a..., edge_b...)
    return meet_type == 0 || meet_type == 1
end

"""
    intersection(geom_a, geom_b)::Union{Tuple{::Real, ::Real}, ::Nothing}

Return an intersection point between two geometries. Return nothing if none are
found. Else, the return type depends on the input. It will be a union between:
a point, a line, a linear ring, a polygon, or a multipolygon

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
GO.intersection(line1, line2)

# output
(125.58375366067547, -14.83572303404496)
```
"""
intersection(geom_a, geom_b) =
    intersection(GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)

"""
    intersection(
        ::GI.LineTrait, line_a,
        ::GI.LineTrait, line_b,
    )::Union{
        ::Tuple{::Real, ::Real},
        ::Nothing
    }

Calculates the intersection between two line segments. Return nothing if
there isn't one.
"""
function intersection(::GI.LineTrait, line_a, ::GI.LineTrait, line_b)
    # Get start and end points for both lines
    a1 = GI.getpoint(line_a, 1)
    a2 = GI.getpoint(line_a, 2)
    b1 = GI.getpoint(line_b, 1)
    b2 = GI.getpoint(line_b, 2)
    # Determine the intersection point
    point, fracs = _intersection_point((a1, a2), (b1, b2))
    # Determine if intersection point is on line segments
    if !isnothing(point) && 0 <= fracs[1] <= 1 && 0 <= fracs[2] <= 1
        return point
    end
    return nothing
end

intersection(
    trait_a::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom_a,
    trait_b::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom_b,
) = intersection_points(trait_a, geom_a, trait_b, geom_b)


function build_ab_list(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    # Make a list for nodes of each polygon. Note the definition of PolyNode
    a_list = Array{PolyNode, 1}(undef, _nedge(poly_a))
    b_list = Array{PolyNode, 1}(undef, _nedge(poly_b))

    # Initialize arrays to keep track of the important information 
    # associated with the intersection points of poly_a and poly_b.
    # I initialize these arrays assumed there will be a maximum of 
    # 30 intersection points and then I truncate the arrays later.
    k = 30
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

    # Find intersection points and adds them to a_list
    # "counter" is used to index all inter-related lists
    counter = 1
    # "acount" is used to index a_list
    acount = 1
    for ii in eachindex(edges_a)
        # Add the first point of the edge to the list of points in a_list
        if acount <= length(a_list)
            a_list[acount] = PolyNode(ii, false, 0, false, 0, false)
        else
            push!(a_list, PolyNode(ii, false, 0, false, 0, false))
        end
        acount = acount + 1

        # Keep track of current position in inter-related lists
        # before finding new intersection points on our new edge
        prev_counter = counter

        for jj in eachindex(edges_b)

            # Add the first point of the edge to b_list
            if ii == 1 
                b_list[jj] = PolyNode(jj, false, 0, false, 0, false)
            end

            # Check if edge jj of poly_b intersects with edge ii of poly_a
            if _line_intersects([edges_a[ii]], [edges_b[jj]]);
                
                int_pt, alphas = _intersection_point(edges_a[ii], edges_b[jj])
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
                pts_to_add[new_order[kk]] = PolyNode(prev_counter+kk-1, true, 0, false, alpha_a_list[prev_counter+kk-1], false)
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
                pts_to_add[new_order[m]] = PolyNode(i[m], true, a_idx_list[i[m]], false, alpha_b_list[i[m]], false)
                b_neighbors[i[m]] = ii + new_order[m]
            end   
            # I use splice instead of insert so I can insert array   
            splice!(b_list, ii+1:ii, pts_to_add)
            skip = true
            num_skips = length(i)
        end
    end

    b_neighbors = b_neighbors[1:counter-1]
    # Iterate through a_list and update the neighbor indices
    for ii in eachindex(a_list)
        if a_list[ii].inter
            a_list[ii].neighbor = b_neighbors[a_list[ii].idx]
        end
    end

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

    return a_list, b_list, a_idx_list, intr_list, edges_a, edges_b
end

function trace_intersection(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
    # Pre-allocate array for return polygons
    return_polys = Vector{Vector{Tuple{Float64, Float64}}}(undef, 0)
    # TODO: Make sure can remove this counter
    # counter = 0

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
                    idx = idx + 1
                else
                    idx = idx -1
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
        # TODO: Make sure can remove this counter
        # counter = counter + 1
    end

    # Check if one polygon totally within other
    # TODO: use point list instead of edges once to_points is fixed
    if isempty(return_polys)
        if point_in_polygon(edges_a[1][1], poly_b)[1]
            list = []
            for i in eachindex(edges_a)
                push!(list, edges_a[i][1])
            end
            push!(list, edges_a[1][1])
            push!(return_polys, list)
        elseif point_in_polygon(edges_b[1][1], poly_a)[1]
            list = []
            for i in eachindex(edges_b)
                push!(list, edges_b[i][1])
            end
            push!(list, edges_b[1][1])
            push!(return_polys, list)
        end
    end
    return return_polys
end

"""
    intersection(
        ::GI.PolygonTrait, poly_a,
        ::GI.PolygonTrait, poly_b,
    )::Union{
        ::Vector{Vector{Tuple{::Real, ::Real}}}, # is this a good return type?
        ::Nothing
    }

Calculates the intersection between two line segments. Return nothing if
there isn't one.
## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0), (5.0, -5.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(3.0, 0.0), (8.0, 5.0), (13.0, 0.0), (8.0, -5.0), (3.0, 0.0)]])
GO.intersection(p1, p2)

# output
1-element Vector{Vector{Tuple{Float64, Float64}}}:
[(6.5, 3.5), (10.0, 0.0), (6.5, -3.5), (3.0, 0.0), (6.5, 3.5)]
```
"""
function intersection(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    a_list, b_list, a_idx_list, intr_list, edges_a, edges_b = build_ab_list(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b)
    return trace_intersection(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
end

function trace_union(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
    # Pre-allocate array for return polygons
    return_polys = Vector{Vector{Tuple{Float64, Float64}}}(undef, 0)
    # TODO: Make sure can remove this counter
    # counter = 0

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
        # TODO: Make sure can remove this counter
        # counter = counter + 1
    end

    # Check if one polygon totally within other
    # TODO: use point list instead of edges once to_points is fixed
    if isempty(return_polys)
        if point_in_polygon(edges_a[1][1], poly_b)[1]
            list = []
            for i in eachindex(edges_b)
                push!(list, edges_b[i][1])
            end
            push!(list, edges_b[1][1])
            push!(return_polys, list)
        elseif point_in_polygon(edges_b[1][1], poly_a)[1]
            list = []
            for i in eachindex(edges_a)
                push!(list, edges_a[i][1])
            end
            push!(list, edges_a[1][1])
            push!(return_polys, list)
        end
    end

    if length(return_polys) > 1
        outer_idx = 1
        for j = 2:(length(return_polys)-1)
            poly1 = GI.Polygon([return_polys[outer_idx]])
            poly2 = GI.Polygon([return_polys[j]])
            if !poly_in_poly(poly2, poly1)
                outer_idx = j
            end
        end

        if outer_idx == 1
            return return_polys
        elseif outer_idx == length(return_polys)
            arr = outer_idx
            splice!(arr, 2:1, 1:(length(return_polys)-1))
            return return_polys[arr]
        else
            arr = 1:(outer_idx-1)
            splice!(arr, (length(arr)+1):length(arr), (outer_idx+1):length(return_polys))
            return return_polys[arr]
        end
    end

    return return_polys
end

function union_test(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    a_list, b_list, a_idx_list, intr_list, edges_a, edges_b = build_ab_list(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b)
    return trace_union(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
end

union_test(geom_a, geom_b) =
    union_test(GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)

function trace_difference(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
    # Pre-allocate array for return polygons
    return_polys = Vector{Vector{Tuple{Float64, Float64}}}(undef, 0)
    # TODO: Make sure can remove this counter
    # counter = 0

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
                if (!status2 && list == a_list) || (status2 && list == b_list)
                    idx = idx + 1
                else
                    idx = idx - 1
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
        # TODO: Make sure can remove this counter
        # counter = counter + 1
    end

    # # Check if one polygon totally within other
    # # TODO: use point list instead of edges once to_points is fixed
    if isempty(return_polys)
        if point_in_polygon(edges_a[1][1], poly_b)[1]
            return return_polys
        end
        list_b = []
        for i in eachindex(edges_b)
            push!(list_b, edges_b[i][1])
        end
        push!(list_b, edges_b[1][1])

        list_a = []
        for i in eachindex(edges_a)
            push!(list_a, edges_a[i][1])
        end
        push!(list_a, edges_a[1][1])
        
        if point_in_polygon(edges_b[1][1], poly_a)[1]
            push!(return_polys, list_a)
            push!(return_polys, list_b)
            return return_polys
        end
    end
    return return_polys
end
    
function difference_test(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    a_list, b_list, a_idx_list, intr_list, edges_a, edges_b = build_ab_list(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b)
    return trace_difference(GI.trait(poly_a), poly_a, GI.trait(poly_b), poly_b, a_list, b_list, a_idx_list, intr_list, edges_a, edges_b)
end

difference_test(geom_a, geom_b) =
    difference_test(GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)

"""
    intersection(
        ::GI.AbstractTrait, geom_a,
        ::GI.AbstractTrait, geom_b,
    )::Union{
        ::Vector{Vector{Tuple{::Real, ::Real}}}, # is this a good return type?
        ::Nothing
    }

Calculates the intersection between two line segments. Return nothing if
there isn't one.
"""
function intersection(
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b,
)
    @assert(
        false,
        "Intersection between $trait_a and $trait_b isn't implemented yet.",
    )
    return nothing
end

"""
    intersection_points(
        geom_a,
        geom_b,
    )::Union{
        ::Vector{::Tuple{::Real, ::Real}},
        ::Nothing,
    }

Return a list of intersection points between two geometries. If no intersection
point was possible given geometry extents, return nothing. If none are found,
return an empty list.
"""
intersection_points(geom_a, geom_b) =
    intersection_points(GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)

"""
    intersection_points(
        ::GI.AbstractTrait, geom_a,
        ::GI.AbstractTrait, geom_b,
    )::Union{
        ::Vector{::Tuple{::Real, ::Real}},
        ::Nothing,
    }

Calculates the list of intersection points between two geometries, inlcuding
line segments, line strings, linear rings, polygons, and multipolygons. If no
intersection points were possible given geometry extents, return nothing. If
none are found, return an empty list.
"""
function intersection_points(::GI.AbstractTrait, a, ::GI.AbstractTrait, b)
    # Check if the geometries extents even overlap
    Extents.intersects(GI.extent(a), GI.extent(b)) || return nothing
    # Create a list of edges from the two input geometries
    edges_a, edges_b = map(sort! ∘ to_edges, (a, b))
    npoints_a, npoints_b  = length(edges_a), length(edges_b)
    a_closed = npoints_a > 1 && edges_a[1][1] == edges_a[end][1]
    b_closed = npoints_b > 1 && edges_b[1][1] == edges_b[end][1]
    if npoints_a > 0 && npoints_b > 0
        # Initialize an empty list of points
        T = typeof(edges_a[1][1][1]) # x-coordinate of first point in first edge
        result = Tuple{T,T}[]
        # Loop over pairs of edges and add any intersection points to results
        for i in eachindex(edges_a)
            for j in eachindex(edges_b)
                point, fracs = _intersection_point(edges_a[i], edges_b[j])
                if !isnothing(point)
                    #=
                    Determine if point is on edge (all edge endpoints excluded
                    except for the last edge for an open geometry)
                    =#
                    α, β = fracs
                    on_a_edge = (!a_closed && i == npoints_a && 0 <= α <= 1) ||
                        (0 <= α < 1)
                    on_b_edge = (!b_closed && j == npoints_b && 0 <= β <= 1) ||
                        (0 <= β < 1)
                    if on_a_edge && on_b_edge
                        push!(result, point)
                    end
                end
            end
        end
        return result
    end
    return nothing
end

"""
    _intersection_point(
        (a1, a2)::Tuple,
        (b1, b2)::Tuple,
    )

Calculates the intersection point between two lines if it exists, and as if the
line extended to infinity, and the fractional component of each line from the
initial end point to the intersection point.
Inputs:
    (a1, a2)::Tuple{Tuple{::Real, ::Real}, Tuple{::Real, ::Real}} first line
    (b1, b2)::Tuple{Tuple{::Real, ::Real}, Tuple{::Real, ::Real}} second line
Outputs:
    (x, y)::Tuple{::Real, ::Real} intersection point
    (t, u)::Tuple{::Real, ::Real} fractional length of lines to intersection
    Both are ::Nothing if point doesn't exist!

Calculation derivation can be found here:
    https://stackoverflow.com/questions/563198/
"""
function _intersection_point((a1, a2)::Tuple, (b1, b2)::Tuple)
    # First line runs from p to p + r
    px, py = GI.x(a1), GI.y(a1)
    rx, ry = GI.x(a2) - px, GI.y(a2) - py
    # Second line runs from q to q + s 
    qx, qy = GI.x(b1), GI.y(b1)
    sx, sy = GI.x(b2) - qx, GI.y(b2) - qy
    # Intersection will be where p + tr = q + us where 0 < t, u < 1 and
    r_cross_s = rx * sy - ry * sx
    if r_cross_s != 0
        Δqp_x = qx - px
        Δqp_y = qy - py
        t = (Δqp_x * sy - Δqp_y * sx) / r_cross_s
        u = (Δqp_x * ry - Δqp_y * rx) / r_cross_s
        x = px + t * rx
        y = py + t * ry
        return (x, y), (t, u)
    end
    return nothing, nothing
end


mutable struct PolyNode
    idx::Int
    inter::Bool
    neighbor::Int
    ent_exit::Bool
    alpha::Real
    proc::Bool
end