# # Intersection checks

export intersects, intersection

# This code checks whether geometries intersect with each other. 

# !!! note
#     This does not compute intersections, only checks if they exist.

const MEETS_OPEN = 1
const MEETS_CLOSED = 0

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
line_intersection(line_a, line_b) = line_intersection(trait(line_a), line_a, trait(line_b), line_b)
function line_intersection(::GI.AbstractTrait, a, ::GI.AbstractTrait, b)
    Extents.intersects(GI.extent(a), GI.extent(b)) || return nothing
    result = Tuple{Float64,Float64}[]
    edges_a, edges_b = map(sort! ∘ to_edges, (a, b))
    for edge_a in edges_a
        for edge_b in edges_b
            x = _line_intersection(edge_a, edge_b)
            isnothing(x) || push!(result, x)
        end
    end
    return result
end
function line_intersection(::GI.LineTrait, line_a, ::GI.LineTrait, line_b)
    a1 = GI.getpoint(line_a, 1)
    b1 = GI.getpoint(line_b, 1)
    a2 = GI.getpoint(line_a, 2)
    b2 = GI.getpoint(line_b, 2)

    return _line_intersection((a1, a2), (b1, b2))
end
function _line_intersection((p11, p12)::Tuple, (p21, p22)::Tuple)
    # Get points from lines
    x1, y1 = GI.x(p11), GI.y(p11) 
    x2, y2 = GI.x(p12), GI.y(p12)
    x3, y3 = GI.x(p21), GI.y(p21)
    x4, y4 = GI.x(p22), GI.y(p22)

    d = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1))
    a = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3))
    b = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3))

    if d == 0
        if a == 0 && b == 0
            return nothing
        end
        return nothing
    end

    ã  = a / d
    b̃  = b / d

    if ã  >= 0 && ã  <= 1 && b̃  >= 0 && b̃  <= 1
        x = x1 + (ã  * (x2 - x1))
        y = y1 + (ã  * (y2 - y1))
        return (x, y)
    end

    return nothing
end
