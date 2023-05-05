abstract type SimplifyAlg end

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
    apply(Union{PolygonTrait,AbstractCurveTrait,MultiPoint,PointTrait}, data) do geom
        _simplify(trait(geom), alg, geom)
    end
end
# For Point and MultiPoint traits we do nothing
_simplify(::PointTrait, alg, geom) = geom
_simplify(::MultiPointTrait, alg, geom) = geom
function _simplify(::PolygonTrait, alg, geom)
    # Force treating children as LinearRing
    lrs = map(GI.getgeom(geom)) do g
        rebuild(g, _simplify(LinearRingTrait(), alg, g))
    end
    return rebuild(geom, lrs)
end
# For curves and rings we simplify
_simplify(::AbstractCurveTrait, alg, geom) = rebuild(geom, simplify(alg, tuple_points(geom)))
function _simplify(::LinearRingTrait, alg, geom)
    GI.npoint(geom) < 4 && throw(ArgumentError("Invalid ring, has less than 4 points."))

    # Make a vector of points 
    points = tuple_points(geom)

    # Simplify it once
    simple = _simplify(alg, points)

    # Reduce the tolerance and simplify until its valid
    while !_isvalid(simple)
        alg = settol(alg, alg.tol * 0.9)
        simple = _simplify(alg, points)
    end

    # Close the ring if its not closed
    point_equals_point(simple[begin], simple[end]) || push!(simple, simple[1])

    return rebuild(geom, simple)
end

"""
    RadialDistance <: SimplifyAlg

Simplifies geometries by removing points less than
`tol` distance from the line between its neighboring points.
"""
@kwdef struct RadialDistance <: SimplifyAlg 
    tol::Float64=0.1
end

settol(::RadialDistance, tol) = RadialDistance(tol)

function _simplify(alg::RadialDistance, points::Vector)
    point = previous = points[1]
    new_points = [previous]

    for i in eachindex(points)
        point = points[i]
        if squared_dist(point, previous) > alg.tol^2
            push!(new_points, point)
            previous = point
        end
    end

    !isequal(previous, point) && push!(new_points, point)

    return new_points
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
@kwdef struct DouglasPeucker <: SimplifyAlg
    tol::Float64=0.1
    prefilter::Bool=true
end

settol(alg::DouglasPeucker, tol) = DouglasPeucker(tol, alg.prefilter)

function _simplify(alg::DouglasPeucker, points::Vector)
    length(points) <= 3 && return points
    points = alg.prefilter ? simplify(RadialDistance(alg.tol), points) : points

    # Defined the simplified point vector, starting with the first point
    new_points = [points[1]]
    # Iteratively add simplified points
    _dp_step!(new_points, points, 1, length(points), alg.tol)
    # Make sure the last point is included
    push!(new_points, points[end])

    return new_points
end

function _dp_step!(simplified, points::Vector, first::Integer, last::Integer, tol::Real)
    max_dist = tol
    index = 0

    for i = first+1:last
        dist = squared_segdist(points[i], points[first], points[last])
        if dist > max_dist
            index = i
            max_dist = dist
        end
    end

    if max_dist > tol
        if (index - first > 1) 
            _dp_step!(simplified, points, first, index, tol)
        end
        push!(simplified, points[index])
        if (last - index > 1) 
            _dp_step!(simplified, points, index, last, tol)
        end
    end

    return nothing
end

function _isvalid(ring::Vector)
    length(ring) < 3 && return false
    length(ring) == 3 && point_equals_point(ring[3], ring[1]) && return false
    return true
end

function squared_segdist(p, l1, l2)
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

point_equals_point(g1, g2) = !(GI.x(g1) == GI.x(g2) && GI.y(g1) == GI.y(g2))
tuple_points(geom) = map(p -> (Float64(GI.x(p)), Float64(GI.y(p))), GI.getpoint(geom))

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
    if count(isnothing, (number, ratio, tol)) == 2
        return VisvalingamWhyatt(number, ratio, tol, prefilter)
    else
        error("Must provide one of `number`, `ratio` or `tol` keywords")
    end
end

settol(alg::VisvalingamWhyatt, tol) = VisvalingamWhyatt(alg.number, alg.ratio, tol, alg.prefilter)

function _simplify(alg::VisvalingamWhyatt, points::Vector)
    length(points) <= 2 && return points
    areas = _build_areas(points)

    (; tol, number, ratio) = alg
    isnothing(tol) || return _by_tol(alg, alg.tol, points, areas)
    isnothing(number) || return _by_number(alg, alg.number, points, areas)
    return _by_ratio(alg, alg.ratio, points, areas)
end

_by_tol(alg, tol, points, areas) = points[areas .>= tol]

function _by_number(alg, n, points, areas)
    tol = partialsort(areas, n)
    return _by_tol(alg, tol, points, areas)[1:n]
end

function _by_ratio(alg, r, points, areas)
    if r <= 0 || r > 1
        error("Ratio must be 0 < r <= 1. Got $r")
    end
    return _by_number(alg, round(Int, r * length(points)), points, areas)
end

function _build_areas(points)
    nmax = length(points)
    real_areas = _triangle_areas(points)

    areas = copy(real_areas)
    i = collect(1:nmax)

    min_vert = argmin(areas)
    this_area = areas[min_vert]
    _remove!(areas, min_vert)
    deleteat!(i, min_vert)

    while this_area < Inf
        skip = false

        if min_vert < length(i)
            right_area = _triangle_area(
                points[i[min_vert - 1]],
                points[i[min_vert]],
                points[i[min_vert + 1]],
            )
            if right_area <= this_area
                right_area = this_area
                skip = min_vert == 1
            end

            real_areas[i[min_vert]] = right_area
            areas[min_vert] = right_area
        end

        if min_vert > 2
            left_area = _triangle_area(
                points[i[min_vert - 2]],
                points[i[min_vert - 1]],
                points[i[min_vert]],
            )
            if left_area <= this_area
                left_area = this_area
                skip = min_vert == 2
            end
            real_areas[i[min_vert - 1]] = left_area
            areas[min_vert - 1] = left_area
        end

        if !skip
            min_vert = argmin(areas)
        end
        deleteat!(i, min_vert)
        this_area = areas[min_vert]
        _remove!(areas, min_vert)
    end

    return real_areas
end


# calculates the area of a triangle given its vertices
_triangle_area(p1, p2, p3) =
    abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2])) / 2.0

function _triangle_areas(points)
    result = Array{Float64}(undef, length(points))
    result[1] = result[end] = Inf

    for i in 2:length(result) - 1
        result[i] = _triangle_area(points[i-1], points[i], points[i+1])
    end
    return result
end

_remove!(s, i) = s[i:end-1] .= s[i+1:end]
