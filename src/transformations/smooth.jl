"""
    Chaikin <: SmoothAlg

    Chaikin(; iterations=1)

Smooths geometries using Chaikin's corner-cutting algorithm.

## Keywords
- `iterations`: the number of times to apply the algorithm.
"""
@kwdef struct Chaikin{M} <: Algorithm{M}
    manifold::M = Planar()
    iterations::Int = 1
end

"""
    smooth(alg::Algorithm, geom)
    smooth(geom; kw...)

Smooths a geometry using the provided algorithm.

The default algorithm is [`Chaikin()`](@ref), which can be used on the spherical or planar manifolds.
"""
smooth(geom; kw...) = smooth(Chaikin(; kw...), geom)

function smooth(alg::Algorithm, geom)
    _smooth_function(trait, geom) = _smooth(alg, trait, geom)
    return apply(
        WithTrait(_smooth_function),
        TraitTarget{Union{GI.AbstractCurveTrait,GI.MultiPointTrait,GI.PointTrait}}(),
        geom,
    )
end

_smooth(alg, ::GI.PointTrait, geom) = geom
_smooth(alg, ::GI.MultiPointTrait, geom) = geom

function _smooth(alg::Chaikin, T::Union{GI.LineStringTrait,GI.LinearRingTrait}, geom)
    isring = T isa GI.LinearRingTrait
    points = tuple_points(geom)
    if isring && first(points) != last(points)
        push!(points, first(points))
    end
    smoothed_points = _chaikin_smooth(points, alg.iterations, isring)
    return rebuild(geom, smoothed_points)
end

function _chaikin_smooth(points::Vector{P}, iterations::Int, isring::Bool) where P
    # points is expected to be a vector of points
    smoothed_points = points
    for itr in 1:iterations
        num_points = length(smoothed_points)
        if isring 
            n = 1
            new_points = Vector{P}(undef, num_points * 2 - 1)
        else
            n = 2
            # Need to add the first point
            new_points = Vector{P}(undef, num_points * 2)
            new_points[begin] = smoothed_points[begin]
            new_points[end] = smoothed_points[end]
        end
        fill!(new_points, (-9999.0, -9999.0))

        for i in eachindex(smoothed_points)[begin:end-1]
            p1 = smoothed_points[i]
            p2 = smoothed_points[i+1]
            _add_smoothed_points!(new_points, p1, p2, n)
            n += 2
        end
 
        if isring # Close it
            new_points[end] = new_points[begin]
        end

        smoothed_points = new_points
    end

    return smoothed_points
end

function _add_smoothed_points!(new_points, p1, p2, n)
    q_x = 0.75 * GI.x(p1) + 0.25 * GI.x(p2)
    q_y = 0.75 * GI.y(p1) + 0.25 * GI.y(p2)
    r_x = 0.25 * GI.x(p1) + 0.75 * GI.x(p2)
    r_y = 0.25 * GI.y(p1) + 0.75 * GI.y(p2)

    new_points[n] = (q_x, q_y)
    new_points[n+1] = (r_x, r_y)
end
