#=
# Perimeter

The perimeter of a geometry is the length of its boundary.  In many contexts
this is called `length`; to avoid clashing with `Base.length` in Julia, we
call this `perimeter`.

## Examples

=#


function perimeter(geom, ::Type{T} = Float64; threaded=False(), init = zero(T), kwargs...) where T
    perimeter(Planar(), geom, T; threaded, init, kwargs...)
end

#=
The planar implementation is straightforward.  
=#

function perimeter(::Planar, geom, ::Type{T} = Float64; init = zero(T), kwargs...) where T
    function _perimeter_planar_inner(trait, geom)
        @assert GI.npoint(geom) >= 2 "Planar perimeter requires at least 2 points"
        distance = zero(T)
        for (p1, p2) in eachedge(trait, geom, T)
            distance += hypot(GI.x(p2) - GI.x(p1), GI.y(p2) - GI.y(p1))
        end
        return distance
    end
    return applyreduce(
        WithTrait(_perimeter_planar_inner), 
        +, 
        TraitTarget(GI.AbstractCurveTrait), 
        geom; init, kwargs...
    )
end

using .UnitSpherical: UnitSphericalPoint

function perimeter(m::Spherical, geom, ::Type{T} = Float64; init = zero(T), kwargs...) where T
    function _perimeter_spherical_inner(trait, geom)
        @assert GI.npoint(geom) >= 2 "Spherical perimeter requires at least 2 points"
        p1_unknown, rest = Iterators.peel(GI.getpoint(trait, geom))
        p1 = UnitSphericalPoint(GI.PointTrait(), p1_unknown)
        distance = zero(T)
        for p2 in Iterators.map(p -> UnitSphericalPoint(GI.PointTrait(), p), rest)
            distance += spherical_distance(p1, p2)
            p1 = p2
        end
        return distance
    end
    return applyreduce(
        WithTrait(_perimeter_spherical_inner), 
        +, 
        TraitTarget(GI.AbstractCurveTrait), 
        geom; init, kwargs...
    ) * m.radius
end

# The `Geodesic` implementation is in `ext/GeometryOpsProjExt/perimeter.jl`