abstract type SimplifyAlg end

const MIN_POINTS = 3

function checkargs(number, ratio, tol)
    count(isnothing, (number, ratio, tol)) == 2 ||
        error("Must provide one of `number`, `ratio` or `tol` keywords")
    if !isnothing(ratio)
        if ratio <= 0 || ratio > 1
            error("`ratio` must be 0 < ratio <= 1. Got $ratio")
        end
    end
    if !isnothing(number)
        if number < MIN_POINTS
            error("`number` must be $MIN_POINTS or larger. Got $number")
        end
    end
    return nothing
end

"""
    simplify(obj; tol=0.1, prefilter=true)
    simplify(::SimplifyAlg, obj)

Simplify a geometry, feature, feature collection, 
or nested vectors or a table of these.

[`RadialDistance`](@ref), [`DouglasPeucker`](@ref), or 
[`VisvalingamWhyatt`](@ref) algorithms are available, 
listed in order of increasing quality but decreaseing performance.

`PoinTrait` and `MultiPointTrait` are returned unchanged.

The default behaviour is `DouglasPeucker(; tol=0.1, prefilter=true)`.

Pass in constructed `SimplifyAlg`s to use other algorithms.

# Examples

```jldoctest
julia> 
import GeoInterface as GI
import GeometryOps as GO

poly = GI.Polygon([[
    [-70.603637, -33.399918],
    [-70.614624, -33.395332],
    [-70.639343, -33.392466],
    [-70.659942, -33.394759],
    [-70.683975, -33.404504],
    [-70.697021, -33.419406],
    [-70.701141, -33.434306],
    [-70.700454, -33.446339],
    [-70.694274, -33.458369],
    [-70.682601, -33.465816],
    [-70.668869, -33.472117],
    [-70.646209, -33.473835],
    [-70.624923, -33.472117],
    [-70.609817, -33.468107],
    [-70.595397, -33.458369],
    [-70.587158, -33.442901],
    [-70.587158, -33.426283],
    [-70.590591, -33.414248],
    [-70.594711, -33.406224],
    [-70.603637, -33.399918]]])

julia> GO.simplify(DouglasPeucker(; tol=0.01), poly)
Polygon(Array{Array{Float64,1},1}[[[-70.6036, -33.3999], [-70.684, -33.4045], [-70.7011, -33.4343], [-70.6943, -33.4584], [-70.6689, -33.4721], [-70.6098, -33.4681], [-70.5872, -33.4429], [-70.6036, -33.3999]]])
```
"""
simplify(data; tol=0.1, prefilter=true) = _simplify(DouglasPeucker(; tol, prefilter), data)
simplify(alg::SimplifyAlg, data) = _simplify(alg, data)

function _simplify(alg::SimplifyAlg, data)
    # Apply simplication to all curves, multipoints, and points,
    # reconstructing everything else around them.
    simplifier(geom) = _simplify(trait(geom), alg, geom)
    apply(simplifier, Union{PolygonTrait,AbstractCurveTrait,MultiPoint,PointTrait}, data)
end
# For Point and MultiPoint traits we do nothing
_simplify(::PointTrait, alg, geom) = geom
_simplify(::MultiPointTrait, alg, geom) = geom
function _simplify(::PolygonTrait, alg, geom)
    # Force treating children as LinearRing
    rebuilder(g) = rebuild(g, _simplify(LinearRingTrait(), alg, g))
    lrs = map(rebuilder, GI.getgeom(geom))
    return rebuild(geom, lrs)
end
# For curves and rings we simplify
_simplify(::AbstractCurveTrait, alg, geom) = rebuild(geom, simplify(alg, tuple_points(geom)))
function _simplify(::LinearRingTrait, alg, geom)
    # Make a vector of points 
    points = tuple_points(geom)

    # Simplify it once
    simple = _simplify(alg, points)

    return rebuild(geom, simple)
end

"""
    RadialDistance <: SimplifyAlg

Simplifies geometries by removing points less than
`tol` distance from the line between its neighboring points.
"""
struct RadialDistance <: SimplifyAlg 
    number::Union{Int64,Nothing}
    ratio::Union{Float64,Nothing}
    tol::Union{Float64,Nothing}
end
function RadialDistance(; number=nothing, ratio=nothing, tol=nothing)
    checkargs(number, ratio, tol)
    return RadialDistance(number, ratio, tol)
end

settol(alg::RadialDistance, tol) = RadialDistance(alg.number, alg.radius, tol)

function _simplify(alg::RadialDistance, points::Vector)
    previous = first(points)
    distances = Array{Float64}(undef, length(points))
    for i in eachindex(points)
        point = points[i]
        distances[i] = squared_dist(point, previous)
        previous = point
    end
    # This avoids taking the square root of each distance above
    if !isnothing(alg.tol)
        alg = settol(alg, tol^2)
    end
    return _get_points(alg, points, distances)
end

function squared_dist(p1, p2)
    dx = GI.x(p1) - GI.x(p2)
    dy = GI.y(p1) - GI.y(p2)
    return dx^2 + dy^2
end

"""
    DouglasPeucker <: SimplifyAlg

Simplifies geometries by removing points below `tol`
distance from the line between its neighboring points.
"""
struct DouglasPeucker <: SimplifyAlg
    number::Union{Int64,Nothing}
    ratio::Union{Float64,Nothing}
    tol::Union{Float64,Nothing}
    prefilter::Bool
end
function DouglasPeucker(; number=nothing, ratio=nothing, tol=nothing, prefilter=false)
    checkargs(number, ratio, tol)
    return DouglasPeucker(number, ratio, tol, prefilter)
end

settol(alg::DouglasPeucker, tol) = DouglasPeucker(tol, alg.prefilter)

function _simplify(alg::DouglasPeucker, points::Vector)
    length(points) <= MIN_POINTS && return points
    # points = alg.prefilter ? simplify(RadialDistance(alg.tol), points) : points

    distances = _build_tolerances(_squared_segdist, points)
    @assert length(distances) == length(points)
    return _get_points(alg, points, distances)

    # Defined the simplified point vector, starting with the first point
    # new_points = [points[1]]
    # # Iteratively add simplified points
    # _dp_step!(new_points, points, 1, length(points), alg.tol)
    # # Make sure the last point is included
    # push!(new_points, points[end])
end

# function _dp_step!(simplified, points::Vector, first::Integer, last::Integer, tol::Real)
#     max_dist = tol
#     index = 0

#     for i = first+1:last
#         dist = squared_segdist(points[i], points[first], points[last])
#         if dist > max_dist
#             index = i
#             max_dist = dist
#         end
#     end

#     if max_dist > tol
#         if (index - first > 1) 
#             _dp_step!(simplified, points, first, index, tol)
#         end
#         push!(simplified, points[index])
#         if (last - index > 1) 
#             _dp_step!(simplified, points, index, last, tol)
#         end
#     end

#     return nothing
# end

# function _isvalid(ring::Vector)
#     length(ring) < 3 && return false
#     length(ring) == 3 && point_equals_point(ring[3], ring[1]) && return false
#     return true
# end
# point_equals_point(g1, g2) = !(GI.x(g1) == GI.x(g2) && GI.y(g1) == GI.y(g2))

function _squared_segdist(l1, p, l2)
    x, y = GI.x(l1), GI.y(l1)
    dx = GI.x(l2) - x
    dy = GI.y(l2) - y

    if !iszero(dx) || !iszero(dy)
        t = ((GI.x(p) - x) * dx + (GI.y(p) - y) * dy) / (dx * dx + dy * dy)

        if t > 1
            x = GI.x(l2)
            y = GI.y(l2)
        elseif t > 0
            x += dx * t
            y += dy * t
        end
    end

    dx = GI.x(p) - x
    dy = GI.y(p) - y

    return dx^2 + dy^2
end


"""
    VisvalingamWhyatt <: SimplifyAlg

    VisvalingamWhyatt(; kw...)

Simplifies geometries by removing points below `tol`
distance from the line between its neighboring points.

# Keywords

- `number`:
- `ratio`:
- `tol`:
- `prefilter`: wether to use a `RadialDistance()` prefilter - 
    for better performance with a slight loss of quality.
"""
struct VisvalingamWhyatt <: SimplifyAlg 
    number::Union{Int,Nothing}
    ratio::Union{Float64,Nothing}
    tol::Union{Float64,Nothing}
    prefilter::Bool
end
function VisvalingamWhyatt(; number=nothing, ratio=nothing, tol=nothing, prefilter=false)
    checkargs(number, ratio, tol)
    return VisvalingamWhyatt(number, ratio, tol, prefilter)
end

settol(alg::VisvalingamWhyatt, tol) = VisvalingamWhyatt(alg.number, alg.ratio, tol, alg.prefilter)

function _simplify(alg::VisvalingamWhyatt, points::Vector)
    length(points) <= MIN_POINTS && return points
    areas = _build_tolerances(_triangle_double_area, points)

    return _get_points(alg, points, areas)
end

# calculates the area of a triangle given its vertices
_triangle_double_area(p1, p2, p3) =
    abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]))


# Shared utils

function _build_tolerances(f, points)
    nmax = length(points)
    real_tolerances = _flat_tolerances(f, points)

    tolerances = copy(real_tolerances)
    i = collect(1:nmax)

    min_vert = argmin(tolerances)
    this_tolerance = tolerances[min_vert]
    _remove!(tolerances, min_vert)
    deleteat!(i, min_vert)

    while this_tolerance < Inf
        skip = false

        if min_vert < length(i)
            right_tolerance = f(
                points[i[min_vert - 1]],
                points[i[min_vert]],
                points[i[min_vert + 1]],
            )
            if right_tolerance <= this_tolerance
                right_tolerance = this_tolerance
                skip = min_vert == 1
            end

            real_tolerances[i[min_vert]] = right_tolerance
            tolerances[min_vert] = right_tolerance
        end

        if min_vert > 2
            left_tolerance = f(
                points[i[min_vert - 2]],
                points[i[min_vert - 1]],
                points[i[min_vert]],
            )
            if left_tolerance <= this_tolerance
                left_tolerance = this_tolerance
                skip = min_vert == 2
            end
            real_tolerances[i[min_vert - 1]] = left_tolerance
            tolerances[min_vert - 1] = left_tolerance
        end

        if !skip
            min_vert = argmin(tolerances)
        end
        deleteat!(i, min_vert)
        this_tolerance = tolerances[min_vert]
        _remove!(tolerances, min_vert)
    end

    return real_tolerances
end

function tuple_points(geom)
    points = Array{Tuple{Float64,Float64}}(undef, GI.ngeom(geom))
    for (i, p) in enumerate(GI.getpoint(geom))
        points[i] = (GI.x(p), GI.y(p))
    end
    return points
end

function _get_points(alg, points, tolerances)
    (; tol, number, ratio) = alg
    bit_indices = if !isnothing(tol) 
        _tol_indices(alg.tol, points, tolerances)
    elseif !isnothing(number) 
        _number_indices(alg.number, points, tolerances)
    else
        _ratio_indices(alg.ratio, points, tolerances)
    end
    return points[bit_indices]
end

function _tol_indices(tol, points, tolerances)
    tolerances .>= tol
end

function _number_indices(n, points, tolerances)
    tol = partialsort(tolerances, length(points) - n + 1)
    bit_indices = _tol_indices(tol, points, tolerances)
    nselected = sum(bit_indices)
    # If there are multiple values exactly at `tol` we will get 
    # the wrong output length. So we need to remove some.
    while nselected > n
        min_tol = Inf
        min_i = 0
        for i in eachindex(bit_indices)
            bit_indices[i] || continue
            if tolerances[i] < min_tol
                min_tol = tolerances[i]
                min_i = i
            end
        end
        nselected -= 1
        bit_indices[min_i] = false
    end
    return bit_indices 
end

function _ratio_indices(r, points, tolerances)
    n = max(3, round(Int, r * length(points)))
    return _number_indices(n, points, tolerances)
end

function _flat_tolerances(f, points)
    result = Array{Float64}(undef, length(points))
    result[1] = result[end] = Inf

    for i in 2:length(result) - 1
        result[i] = f(points[i-1], points[i], points[i+1])
    end
    return result
end

_remove!(s, i) = s[i:end-1] .= s[i+1:end]
