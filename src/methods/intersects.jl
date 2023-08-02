# # Intersection checks

export intersects, intersection

# This code checks whether geometries intersect with each other. 

# !!! note
#     This does not compute intersections, only checks if they exist.


"""
    intersects(line_a, line_b)

Check if `line_a` intersects with `line_b`.

These can be `LineTrait`, `LineStringTrait` or `LinearRingTrait`

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
function intersects(a, b)
    Extents.intersects(GI.extent(a), GI.extent(b)) || return false
    return !isnothing(intersection(a, b)) # Probably faster ways to do this
end

"""
    intersection(line_a, line_b)

Find a point that intersects LineStrings with two coordinates each.

Returns `nothing` if no point is found.

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
intersection(line_a, line_b) = _intersection(trait(line_a), line_a, trait(line_b), line_b)
function _intersection(
    ::Union{LineStringTrait,LinearRingTrait}, line_a, 
    ::Union{LineStringTrait,LinearRingTrait}, line_b,
)
    result = Tuple{Float64,Float64}[] # TODO handle 3d, and other Real ?
    Extents.intersects(GI.extent(line_a), GI.extent(line_b)) || return result

    # a1 = GI.getpoint(line_a, 1)
    # b1 = GI.getpoint(line_b, 1)

    for i in 1:GI.npoint(line_a) - 1
        a1 = GI.getpoint(line_a, i)
        a2 = GI.getpoint(line_a, i + 1)
        for j in 1:GI.npoint(line_b) - 1
            b1 = GI.getpoint(line_b, j)
            b2 = GI.getpoint(line_b, j + 1)
            inter = _intersection((a1, a2), (b1, b2))
            @show a1 a2 b1 b2 inter
            isnothing(inter) || push!(result, inter)
            # b1 = b2
        end
        # a1 = a2
    end
    return unique!(result)
end
function _intersection(::LineTrait, line_a, ::LineTrait, line_b)
    a1 = GI.getpoint(line_a, 1)
    b1 = GI.getpoint(line_b, 1)
    a2 = GI.getpoint(line_a, 2)
    b2 = GI.getpoint(line_b, 2)

    return _intersection((a1, a2), (b1, b2))
end
function _intersection((p11, p12)::Tuple, (p21, p22)::Tuple)
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
    else
        @show d a b
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
