"""
    abstract type SmoothAlg

Abstract type for smoothing algorithms.
"""
abstract type SmoothAlg end

"""
    Chaikin <: SmoothAlg

    Chaikin(; iterations=1)

Smooths geometries using Chaikin's corner-cutting algorithm.

## Keywords
- `iterations`: the number of times to apply the algorithm.
"""
@kwdef struct Chaikin <: SmoothAlg
    iterations::Int = 1
end

"""
    smooth(alg::SmoothAlg, geom)
    smooth(geom; kw...)

Smooths a geometry using the provided algorithm.

The default algorithm is `Chaikin()`.
"""
smooth(geom; kw...) = smooth(Chaikin(; kw...), geom)

function smooth(alg::SmoothAlg, geom)
    _smooth_function(geom) = _smooth(alg, GI.trait(geom), geom)
    return apply(
        _smooth_function,
        TraitTarget{Union{GI.AbstractCurveTrait,GI.MultiPointTrait,GI.PointTrait}}(),
        geom,
    )
end

_smooth(alg, ::GI.PointTrait, geom) = geom
_smooth(alg, ::GI.MultiPointTrait, geom) = geom

function _smooth(alg::Chaikin, T::Union{GI.LineStringTrait,GI.LinearRingTrait}, geom)
    points = tuple_points(geom)
    needs_ends = T isa GI.LineStringTrait
    smoothed_points = _chaikin_smooth(points, alg.iterations, needs_ends)
    return rebuild(geom, smoothed_points)
end

function _chaikin_smooth(points::Vector{P}, iterations::Int, needs_ends::Bool) where P
    # points is expected to be a vector of points
    num_points = length(points)
    if num_points < 2
        return points
    end

    smoothed_points = points
    for _ in 1:iterations
        if needs_ends 
            n = 2
            # Need to add the first point
            new_points = Vector{P}(undef, num_points * 2)
            new_points[1] = first(points)
        else
            n = 1
            new_points = Vector{P}(undef, num_points * 2 - 2)
        end

        for i in 1:num_points - 1
            p1 = smoothed_points[i]
            p2 = smoothed_points[i+1]
            
            q_x = 0.75 * GI.x(p1) + 0.25 * GI.x(p2)
            q_y = 0.75 * GI.y(p1) + 0.25 * GI.y(p2)
            r_x = 0.25 * GI.x(p1) + 0.75 * GI.x(p2)
            r_y = 0.25 * GI.y(p1) + 0.75 * GI.y(p2)

            new_points[n] = (q_x, q_y)
            new_points[n+1] = (r_x, r_y)
            n += 2
        end

        if needs_ends
            # for open curves, we need to add the first and last points
            # but the loop above does not do that.
            # instead we can just add them here
            new_points[end] = points[end]
        end

        smoothed_points = new_points
    end

    return smoothed_points
end