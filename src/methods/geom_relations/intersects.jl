# # Intersection checks

export intersects, intersection, intersection_points

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
"""
function intersection(::GI.PolygonTrait, poly_a, ::GI.PolygonTrait, poly_b)
    @assert false "Polygon intersection isn't implemented yet."
    return nothing
end

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
