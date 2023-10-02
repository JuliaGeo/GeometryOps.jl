# # Intersection checks

export intersects, intersection

#=
## What is `intersects` vs `intersection`?

The `intersects` methods check whether two geometries intersect with each other.
The `intersection` methods return the intersection between the two geometries.

The `intersects` methods will always return a Boolean. However, note that the
`intersection` methods will not all return the same type. For example, the
intersection of two lines will be a point in most cases, unless the lines are
parallel. On the other hand, the intersection of two polygons will be another
polygon in most cases.

To provide an example, consider this # TODO update this example:
```@example cshape
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

cshape = Polygon([
    Point(0,0), Point(0,3), Point(3,3), Point(3,2), Point(1,2),
    Point(1,1), Point(3,1), Point(3,0), Point(0,0),
])
f, a, p = poly(cshape; axis = (; aspect = DataAspect()))
```
Let's see what the centroid looks like (plotted in red):
```@example cshape
cent = centroid(cshape)
scatter!(a, GI.x(cent), GI.y(cent), color = :red)
f
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait. This is also used in the
implementation, since it's a lot less work! 

# TODO fill this in!
=#

const MEETS_OPEN = 1
const MEETS_CLOSED = 0

intersects(geom1, geom2) = GO.intersects(
    GI.trait(geom1),
    geom1,
    GI.trait(geom2),
    geom2,
)

GO.intersects(
    trait1::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom1,
    trait2::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom2,
) = line_intersects(trait1, geom1, trait2, geom2)

"""
    line_intersects(line_a, line_b)

Check if `line_a` intersects with `line_b`.

These can be `LineTrait`, `LineStringTrait` or `LinearRingTrait`

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
GO.line_intersects(line1, line2)

# output
true
```
"""
line_intersects(a, b; kw...) = line_intersects(trait(a), a, trait(b), b; kw...)
# Skip to_edges for LineTrait
function line_intersects(::GI.LineTrait, a, ::GI.LineTrait, b; meets=MEETS_OPEN)
    a1 = _tuple_point(GI.getpoint(a, 1))
    b1 = _tuple_point(GI.getpoint(b, 1))
    a2 = _tuple_point(GI.getpoint(a, 2))
    b2 = _tuple_point(GI.getpoint(b, 2))
    return ExactPredicates.meet(a1, a2, b1, b2) == meets
end
function line_intersects(::GI.AbstractTrait, a, ::GI.AbstractTrait, b; kw...)
    edges_a, edges_b = map(sort! ∘ to_edges, (a, b))
    return line_intersects(edges_a, edges_b; kw...)
end
function line_intersects(edges_a::Vector{Edge}, edges_b::Vector{Edge}; meets=MEETS_OPEN)
    # Extents.intersects(to_extent(edges_a), to_extent(edges_b)) || return false
    for edge_a in edges_a
        for edge_b in edges_b
            ExactPredicates.meet(edge_a..., edge_b...) == meets && return true 
        end
    end
    return false
end

"""
    line_intersection(line_a, line_b)

Find a point that intersects LineStrings with two coordinates each.

Returns `nothing` if no point is found.

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
GO.line_intersection(line1, line2)

# output
(125.58375366067547, -14.83572303404496)
```
"""
line_intersection(line_a, line_b) = intersection_points(trait(line_a), line_a, trait(line_b), line_b)

"""
    intersection_points(
        ::GI.AbstractTrait, geom_a,
        ::GI.AbstractTrait, geom_b,
    )::Vector{::Tuple{::Real, ::Real}}

Calculates the list of intersection points between two geometries. 
"""
function intersection_points(::GI.AbstractTrait, a, ::GI.AbstractTrait, b)
    Extents.intersects(GI.extent(a), GI.extent(b)) || return nothing
    result = Tuple{Float64,Float64}[]
    edges_a, edges_b = map(sort! ∘ to_edges, (a, b))
    for edge_a in edges_a
        for edge_b in edges_b
            x = _intersection_point(edge_a, edge_b)
            isnothing(x) || push!(result, x)
        end
    end
    return result
end

"""
    intersection_point(
        ::GI.LineTrait, line_a,
        ::GI.LineTrait, line_b,
    )::Union{
        ::Tuple{::Real, ::Real},
        ::Nothing
    }

Calculates the intersection point between two lines if it exists and return
`nothing` if it doesn't exist.
"""
function intersection_point(::GI.LineTrait, line_a, ::GI.LineTrait, line_b)
    # Get start and end points for both lines
    a1 = GI.getpoint(line_a, 1)
    a2 = GI.getpoint(line_a, 2)
    b1 = GI.getpoint(line_b, 1)
    b2 = GI.getpoint(line_b, 2)
    # Determine the intersection point
    point, _ = _intersection_point((a1, a2), (b1, b2))
    return point
end

"""
    _intersection_point(
        (p11, p12)::Tuple,
        (p21, p22)::Tuple,
    )

Calculates the intersection point between two lines if it exists, and the
fractional component of each line from the initial end point to the
intersection point.
Inputs:
    (p11, p12)::Tuple{Tuple{::Real, ::Real}, Tuple{::Real, ::Real}} first line
    (p21, p22)::Tuple{Tuple{::Real, ::Real}, Tuple{::Real, ::Real}} second line
Outputs:
    (x, y)::Tuple{::Real, ::Real} intersection point
    (t, u)::Tuple{::Real, ::Real} fractional length of lines to intersection
    Both are ::Nothing if point doesn't exist!

Calculation derivation can be found here:
    https://stackoverflow.com/questions/563198/
"""
function _intersection_point((p11, p12)::Tuple, (p21, p22)::Tuple)
    # First line runs from p to p + r
    px, py = GI.x(p11), GI.y(p11)
    rx, ry = GI.x(p12) - px, GI.y(p12) - py
    # Second line runs from q to q + s 
    qx, qy = GI.x(p21), GI.y(p21)
    sx, sy = GI.x(p22) - qx, GI.y(p22) - qy
    # Intersection will be where p + tr = q + us where 0 < t, u < 1 and
    r_cross_s = rx * sy - ry * sx
    if r_cross_s != 0
        Δpq_x = px - qx
        Δpq_y = py - qy
        t = (Δpq_x * sy - Δpq_y * sx) / r_cross_s
        u = (Δpq_x * ry - Δpq_y * rx) / r_cross_s
        if 0 <= t <= 1 && 0 <= u <= 1
            x = px + t * rx
            y = py + t * ry
            return (x, y), (t, u)
        end
    end
    return nothing, nothing
end
