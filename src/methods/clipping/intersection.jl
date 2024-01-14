# #  Intersection Polygon Clipping
export intersection, intersection_points

# This 'intersection' implementation returns the intersection of two polygons.
# It returns a Vector{Vector{Vector{Tuple{Float}}}. The Vector{Vector{Tuple{Float}
# is empty if the two polygons don't intersect. The algorithm to determine the 
# intersection was adapted from "Efficient clipping of efficient polygons," by 
# Greiner and Hormann (1998). DOI: https://doi.org/10.1145/274363.274364

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
        ::GI.PolygonTrait, poly_a,
        ::GI.PolygonTrait, poly_b,
    )::Vector{Vector{Vector{Tuple{Float64}}}}    

Calculates the intersection between two polygons. If the intersection is empty, 
the vector of a vector is empty (note the outermost vector is technically not empty).
## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0), (5.0, -5.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(3.0, 0.0), (8.0, 5.0), (13.0, 0.0), (8.0, -5.0), (3.0, 0.0)]])
GO.intersection(p1, p2)

# output
1-element Vector{Vector{Vector{Tuple{Float64, Float64}}}}:
[[[(6.5, 3.5), (10.0, 0.0), (6.5, -3.5), (3.0, 0.0), (6.5, 3.5)]]]
```
"""

function intersection(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    # First we get the exteriors of 'poly_a' and 'poly_b'
    ext_poly_a = GI.getexterior(poly_a)
    ext_poly_a = GI.Polygon([ext_poly_a])
    ext_poly_b = GI.getexterior(poly_b)
    ext_poly_b = GI.Polygon([ext_poly_b])
    # Then we find the intersection of the exteriors
    a_list, b_list, a_idx_list, intr_list = _build_ab_list(ext_poly_a, ext_poly_b)
    polys = _trace_intersection(ext_poly_a, ext_poly_b, a_list, b_list, a_idx_list, intr_list)
    # If the original polygons had no holes, then we are pretty much done. Otherwise,
    # we call '_get_inter_holes' to take into account the holes.
    if GI.nhole(poly_a)==0 && GI.nhole(poly_b)==0
        final_polys =  Vector{Vector{Vector{Tuple{Float64, Float64}}}}(undef, length(polys))
        for i in 1:length(polys)
            final_polys[i] = [polys[i]]
        end
        return final_polys
    else
        return _get_inter_holes(polys, poly_a, poly_b)
    end    

end



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
    _trace_intersection(poly_a, poly_b, a_list, b_list, a_idx_list,
      intr_list)::Vector{Vector{Tuple{Float64}}}

Traces the outlines of two polygons in order to find their intersection.
It returns the outlines of all polygons formed in the intersection. If
they do not intersect, it returns an empty array.

"""
function _trace_intersection(poly_a, poly_b, a_list, b_list, a_idx_list, intr_list)
    # Pre-allocate array for return polygons
    return_polys = Vector{Vector{Tuple{Float64, Float64}}}(undef, 0)

    # Keep track of number of processed intersection points
    processed_pts = 0
    tracker = copy(a_idx_list)

    while processed_pts < length(intr_list)
        # Create variables "list_edges" and "list" so that we can toggle between
        # a_list and b_list
        list = a_list

        # Find index of first unprocessed intersecting point in subject polygon
        starting_pt = minimum(tracker)
        idx = starting_pt

        # Get current first unprocessed intersection point PolyNode
        current = a_list[idx]
        # Initialize array to store the intersection polygon cartesian points
        pt_list = Vector{Tuple{Float64, Float64}}(undef, 0)
        # Add the first point to the array
        push!(pt_list, current.point)
        
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
                    push!(pt_list, current.point)
                    
                    # Keep track of processed intersection points
                    if (current != a_list[starting_pt] && current != b_list[a_list[starting_pt].neighbor])
                        processed_pts = processed_pts + 1
                        tracker[current.idx] = typemax(Int)
                    end
                    
                else
                    # Add cartesian coordinates from "list", which should point to either a_list or b_list
                    push!(pt_list, current.point)
                end

                current_node_not_intersection = !current.inter
            end
            
            # Break once get back to starting point
            current_node_not_starting = (current != a_list[starting_pt] && current != b_list[a_list[starting_pt].neighbor])

            # Switch to neighbor list
            if list == a_list
                list = b_list
            else
                list = a_list
            end
            idx = current.neighbor
            current = list[idx]
        end

        push!(return_polys, pt_list)
    end

    # Check if one polygon totally within other, and if so
    # return the smaller polygon as the intersection
    if isempty(return_polys)
        if within(a_list[1].point, poly_b)[1]
            list = []
            for point in GI.getpoint(poly_a)
                push!(list, _tuple_point(point))
            end
            push!(return_polys, list)
        elseif within(b_list[1].point, poly_a)[1]
            list = []
            for point in GI.getpoint(poly_b)
                push!(list, _tuple_point(point))
            end
            push!(return_polys, list)
        end
    end

    # If the polygons don't intersect and aren't contained within each
    # other, return_polys will be empty
    return return_polys
end

"""
    _get_inter_holes(return_polys, poly_a, poly_b)::Vector{Vector{Vector{Tuple{Float64, Float64}}}}

When the _trace_difference function was called, it only took into account the
exteriors of the two polygons when computing the difference. The function
'_get_difference_holes' takes into account the holes of the original polygons
and adjust the output of _trace_difference (return_polys) accordingly.

"""

function _get_inter_holes(return_polys, poly_a, poly_b)
    # Initiaze our return object
    final_polys =  Vector{Vector{Vector{Tuple{Float64, Float64}}}}(undef, 0)

    for poly in return_polys
        # Turning polygon into the desired return type I can add more polygons to it
        poly = [[poly]]

        # We subtract the holes of 'poly_a' and 'poly_b' from the output we got
        # from _trace_intersection (return_polys)
        for hole in GI.gethole(poly_a) 
            replacement_p = Vector{Vector{Vector{Tuple{Float64, Float64}}}}(undef, 0)
            for p in poly
                # When we take the difference of our existing intersectio npolygons and 
                # the holes of polygon_a, we might split it up into smaller polygons. 
                new_ps = difference(GI.Polygon(p), GI.Polygon([hole]))
                append!(replacement_p, new_ps)
            end
            poly = replacement_p
        end
        
        for hole in GI.gethole(poly_b)
            replacement_p = Vector{Vector{Vector{Tuple{Float64, Float64}}}}(undef, 0)
            for p in poly
                # When we take the difference of our existing intersectio npolygons and 
                # the holes of polygon_a, we might split it up into smaller polygons. 
                new_ps = difference(GI.Polygon(p), GI.Polygon([hole]))
                append!(replacement_p, new_ps)
            end
            poly = replacement_p
        end
        
        append!(final_polys, poly)
    end

    return final_polys
        
end

"""
>>>>>>> main:src/methods/geom_relations/intersects.jl
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


