"""
    intersects(line_a, line_b)

Check if `line_a` intersects with `line_b`.

These can be `LineTrait`, `LineStringTrait` or `LinearRingTrait`
"""
intersects(a, b) = isnothing(intersection) # Probably faster ways to do this

"""
    intersection(line_a, line_b)

Find a point that intersects LineStrings with two coordinates each.

Returns `nothing` if no point is found.

# Examples

```jldoctest
julia> 
line1 = LineString([[124.584961,-12.768946],[126.738281,-17.224758]])
LineString(Array{Float64,1}[[124.585, -12.7689], [126.738, -17.2248]])

julia> line2 = LineString([[123.354492,-15.961329],[127.22168,-14.008696]])
LineString(Array{Float64,1}[[123.354, -15.9613], [127.222, -14.0087]])

julia> 
intersection(line1, line2)
Point([125.584, -14.8357])
```
"""
intersection(line_a, line_b) = intersection(trait(line_a), line_a, trait(line_b), line_b)
function intersection(
    ::Union{LineStringTrait,LinearRingTrait}, line_a, 
    ::Union{LineStringTrait,LinearRingTrait}, line_b,
)
    result = Tuple{Float64,Float64}[] # TODO handle 3d, and other Real ?
    a1 = GI.getpoint(line_a, 1)
    b1 = GI.getpoint(line_b, 1)

    # TODO we can check all of these against the extent 
    # of line_b and continue the loop if theyre outside
    for i in 1:GI.npoint(line1) - 1
        for j in 1:GI.npoint(line_b) - 1
            a2 = GI.getpoint(line_a, i + 1)
            b2 = GI.getpoint(line_b, j + 1)
            inter = _intersection((a1, a2), (b1, b2))
            isnothing(inter) || push!(result, inter)
            a1 = a2
            b1 = b2
        end
    end
    return unique!(result)
end

function intersection(::LineTrait, line_a, ::LineTrait, line_b)
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
