# # Geometry simplification

#=
This file holds implementations for the RadialDistance, Douglas-Peucker, and
Visvalingam-Whyatt algorithms for simplifying geometries (specifically for\
polygons and lines).
=#

export simplify, VisvalingamWhyatt, DouglasPeucker, RadialDistance

const MIN_POINTS = 3
const SIMPLIFY_ALG_KEYWORDS = """
## Keywords

- `ratio`: the fraction of points that should remain after `simplify`. 
    Useful as it will generalise for large collections of objects.
- `number`: the number of points that should remain after `simplify`.
    Less useful for large collections of mixed size objects.
"""
const DOUGLAS_PEUCKER_KEYWORDS = """
$SIMPLIFY_ALG_KEYWORDS
- `tol`: the minimum distance a point will be from the line
    joining its neighboring points.
"""

"""
    abstract type SimplifyAlg

Abstract type for simplification algorithms.

## API

For now, the algorithm must hold the `number`, `ratio` and `tol` properties.  

Simplification algorithm types can hook into the interface by implementing 
the `_simplify(trait, alg, geom)` methods for whichever traits are necessary.
"""
abstract type SimplifyAlg end

"""
    simplify(obj; kw...)
    simplify(::SimplifyAlg, obj; kw...)

Simplify a geometry, feature, feature collection, 
or nested vectors or a table of these.

[`RadialDistance`](@ref), [`DouglasPeucker`](@ref), or 
[`VisvalingamWhyatt`](@ref) algorithms are available, 
listed in order of increasing quality but decreaseing performance.

`PoinTrait` and `MultiPointTrait` are returned unchanged.

The default behaviour is `simplify(DouglasPeucker(; kw...), obj)`.
Pass in other [`SimplifyAlg`](@ref) to use other algorithms.

# Keywords

$APPLY_KEYWORDS

Keywords for DouglasPeucker are allowed when no algorithm is specified:

$DOUGLAS_PEUCKER_KEYWORDS

# Example

Simplify a polygon to have six points:

```jldoctest
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

simple = GO.simplify(poly; number=6)
GI.npoint(simple)

# output
6
```
"""
simplify(alg::SimplifyAlg, data; kw...) = _simplify(alg, data; kw...)
# Default algorithm is DouglasPeucker
simplify(data; prefilter = nothing, calc_extent=false, threaded=false, crs=nothing, kw...) =
    _simplify(DouglasPeucker(; kw...), data; prefilter, calc_extent, threaded, crs)

#= For each algorithm, apply simplication to all curves, multipoints, and
points, reconstructing everything else around them. =#
function _simplify(alg::SimplifyAlg, data; prefilter = nothing, kw...)
    simplifier(geom) = _simplify(GI.trait(geom), alg, geom; prefilter = prefilter)
    return apply(
        simplifier,
        Union{GI.PolygonTrait, GI.AbstractCurveTrait, GI.MultiPointTrait, GI.PointTrait},
        data;
        kw...,
    )
end


## For Point and MultiPoint traits we do nothing
_simplify(::GI.PointTrait, alg, geom; kw...) = geom
_simplify(::GI.MultiPointTrait, alg, geom; kw...) = geom

## For curves, rings, and polygon we simplify
function _simplify(::GI.AbstractCurveTrait, alg, geom; prefilter)
    points = if isnothing(prefilter)
        tuple_points(geom)
    else
        _simplify(prefilter, tuple_points(geom))
    end
    return rebuild(geom, _simplify(alg, points))
end

function _simplify(::GI.PolygonTrait, alg, geom;  kw...)
    ## Force treating children as LinearRing
    simplifier(g) = _simplify(GI.LinearRingTrait(), alg, g; kw...)
    rebuilder(g) = rebuild(g, simplifier(g))
    lrs = map(rebuilder, GI.getgeom(geom))
    return rebuild(geom, lrs)
end


# # Simplify with RadialDistance Algorithm
"""
    RadialDistance <: SimplifyAlg

Simplifies geometries by removing points less than
`tol` distance from the line between its neighboring points.

$SIMPLIFY_ALG_KEYWORDS
- `tol`: the minimum distance between points.

Note: user input `tol` is squared to avoid uneccesary computation in algorithm.
"""
@kwdef struct RadialDistance <: SimplifyAlg 
    number::Union{Int64,Nothing} = nothing
    ratio::Union{Float64,Nothing} = nothing
    tol::Union{Float64,Nothing} = nothing

    function RadialDistance(number, ratio, tol)
        _checkargs(number, ratio, tol)
        # square tolerance for reduced computation
        tol = isnothing(tol) ? tol : tol^2
        new(number, ratio, tol)
    end
end

function _simplify(alg::RadialDistance, points::Vector)
    previous = first(points)
    distances = Array{Float64}(undef, length(points))
    for i in eachindex(points)
        point = points[i]
        distances[i] = _squared_euclid_distance(Float64, point, previous)
        previous = point
    end
    ## Never remove the end points
    distances[begin] = distances[end] = Inf
    return _get_points(alg, points, distances)
end


# # Simplify with DouglasPeucker Algorithm
"""
    DouglasPeucker <: SimplifyAlg

    DouglasPeucker(; number, ratio, tol)

Simplifies geometries by removing points below `tol`
distance from the line between its neighboring points.

$DOUGLAS_PEUCKER_KEYWORDS
Note: user input `tol` is squared to avoid uneccesary computation in algorithm.
"""
@kwdef struct DouglasPeucker <: SimplifyAlg
    number::Union{Int64,Nothing} = nothing
    ratio::Union{Float64,Nothing} = nothing
    tol::Union{Float64,Nothing} = nothing

    function DouglasPeucker(number, ratio, tol)
        _checkargs(number, ratio, tol)
        # square tolerance for reduced computation
        tol = isnothing(tol) ? tol : tol^2
        return new(number, ratio, tol)
    end
end

function _simplify(alg::DouglasPeucker, points::Vector)
    length(points) <= MIN_POINTS && return points
    pts = if !isnothing(alg.tol)
        _simplify_tol(alg, points, 1, length(points))
    else  # TODO: This isn't the correct implementation of DouglasPeucker
        _simplify_interative(alg, points)
    end
    return pts 
end

# Top down DouglasPeucker when given a tolerance
function _simplify_tol(alg, points, start_idx, end_idx)
    max_idx, max_dist = _find_max_dist(points, start_idx, end_idx)
    results = if max_dist < alg.tol
        [points[start_idx], points[end_idx]]
    else
        vcat(
            _simplify_tol(alg, points, start_idx, max_idx)[1:(end - 1)],
            _simplify_tol(alg, points, max_idx, end_idx),
        )
    end
    return results
end

function _find_max_dist(points, start_idx, end_idx)
    max_idx = 0
    max_dist = zero(Float64)
    for i in (start_idx + 1):(end_idx - 1)
        dist = _squared_distance_line(Float64, points[i], points[start_idx], points[end_idx])
        if dist > max_dist
            max_dist = dist
            max_idx = i
        end
    end
    return max_idx, max_dist
end

function _simplify_interative(alg, points)
    npoints = length(points)
    points_to_keep = isnothing(alg.ratio) ?
        alg.number : max(3, round(Int, alg.ratio * npoints))
    results = Vector{Int}(undef, points_to_keep)
    results[1], results[points_to_keep] = 1, npoints
    queue = Vector{Tuple{Int, Int, Int, Float64}}()
    queue_idx, queue_dist = 0, zero(Float64)
    
    i = 2
    start_idx, end_idx = 1, npoints
    max_idx, max_dist = _find_max_dist(points, start_idx, end_idx)
    for i in 2:(points_to_keep - 1)
        results[i] = max_idx
        left_idx, left_dist = _find_max_dist(points, start_idx, max_idx)
        right_idx, right_dist = _find_max_dist(points, max_idx, end_idx)
        # Add values to queue
        if right_dist > 0 && (left_dist > right_dist || queue_dist > right_dist)
            push!(queue, (max_idx, right_idx, end_idx, right_dist))
            if right_dist > queue_dist
                queue_dist = right_dist
                queue_idx = length(queue)
            end
        end
        if left_dist > 0 && (left_dist ≤ right_dist || queue_dist > left_dist)
            push!(queue, (start_idx, left_idx, max_idx, left_dist))
            if left_idx > queue_dist
                queue_dist = left_dist
                queue_idx = length(queue)
            end
        end
        # Determine next value to add to results and remove value from queue if needed
        if queue_dist > left_dist && queue_dist > right_dist
            start_idx, max_idx, end_idx, max_dist = queue[queue_idx]
            deleteat!(queue, queue_idx)
            queue_dist, queue_idx = !isempty(queue) ?
                findmax(x -> x[4], queue) : (zero(Float64), 0)
        elseif left_dist > right_dist
            end_idx = max_idx
            max_idx, max_dist = left_idx, left_dist
        else
            start_idx = max_idx
            max_idx, max_dist = right_idx, right_dist
        end
    end
    return points[sort!(results)]
end

# # Simplify with VisvalingamWhyatt Algorithm
"""
    VisvalingamWhyatt <: SimplifyAlg

    VisvalingamWhyatt(; kw...)

Simplifies geometries by removing points below `tol`
distance from the line between its neighboring points.

$SIMPLIFY_ALG_KEYWORDS
- `tol`: the minimum area of a triangle made with a point and
    its neighboring points.
Note: user input `tol` is doubled to avoid uneccesary computation in algorithm.
"""
@kwdef struct VisvalingamWhyatt <: SimplifyAlg 
    number::Union{Int,Nothing} = nothing
    ratio::Union{Float64,Nothing} = nothing
    tol::Union{Float64,Nothing} = nothing

    function VisvalingamWhyatt(number, ratio, tol)
        _checkargs(number, ratio, tol)
        # double tolerance for reduced computation
        tol = isnothing(tol) ? tol : tol*2
        return new(number, ratio, tol)
    end
end

function _simplify(alg::VisvalingamWhyatt, points::Vector)
    length(points) <= MIN_POINTS && return points
    areas = _build_tolerances(_triangle_double_area, points)
    return _get_points(alg, points, areas)
end

# Calculates double the area of a triangle given its vertices
_triangle_double_area(p2, p1, p3) =
    abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]))


# # Shared utils

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
                points[i[min_vert]],
                points[i[min_vert - 1]],
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
                points[i[min_vert - 1]],
                points[i[min_vert - 2]],
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
    points = Array{Tuple{Float64,Float64}}(undef, GI.npoint(geom))
    for (i, p) in enumerate(GI.getpoint(geom))
        points[i] = (GI.x(p), GI.y(p))
    end
    return points
end

function _get_points(alg, points, tolerances)
    ## This assumes that `alg` has the properties
    ## `tol`, `number`, and `ratio` available...
    tol = alg.tol
    number = alg.number
    ratio = alg.ratio
    bit_indices = if !isnothing(tol) 
        _tol_indices(alg.tol::Float64, points, tolerances)
    elseif !isnothing(number) 
        _number_indices(alg.number::Int64, points, tolerances)
    else
        _ratio_indices(alg.ratio::Float64, points, tolerances)
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
    ## If there are multiple values exactly at `tol` we will get 
    ## the wrong output length. So we need to remove some.
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
        result[i] = f(points[i], points[i-1], points[i+1])
    end
    return result
end

_remove!(s, i) = s[i:end-1] .= s[i+1:end]

# Check SimplifyAlgs inputs to make sure they are valid for below algorithms
function _checkargs(number, ratio, tol)
    count(isnothing, (number, ratio, tol)) == 2 ||
        error("Must provide one of `number`, `ratio` or `tol` keywords")
    if !isnothing(number)
        if number < MIN_POINTS
            error("`number` must be $MIN_POINTS or larger. Got $number")
        end
    elseif !isnothing(ratio)
        if ratio <= 0 || ratio > 1
            error("`ratio` must be 0 < ratio <= 1. Got $ratio")
        end
    else  # !isnothing(tol)
        if tol ≤ 0
            error("`tol` must be a positive number. Got $tol")
        end
    end
    return nothing
end