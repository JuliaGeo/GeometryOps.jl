#=
# Perimeter

The perimeter of a geometry is the length of its boundary.  In many contexts
this is called `length`; to avoid clashing with `Base.length` in Julia, we
call this `perimeter`.

## Examples

=#

"""
    perimeter([m::Manifold = AutoManifold()], geom, [T = Float64])::T

Returns the total length of all curves in `geom`: the length of a line, or the
length of all rings of a polygon, including holes.

## Manifold support

- `AutoManifold()` (default): selects the manifold from the CRS of `geom`, as
  for [`area`](@ref). With the Proj extension loaded, recognized geographic
  CRSs use a geodesic calculation on the CRS ellipsoid, in metres. Without
  Proj, geographic geometries use `Spherical()`. Projected, unknown, and
  CRS-less geometries use `Planar()`.
- `Planar()`: Euclidean length in native coordinate units, regardless of CRS.
- `Spherical()`: great-circle length, with coordinates interpreted as
  (longitude, latitude) in degrees, in units of the sphere's radius.
- `Geodesic()`: geodesic length on the ellipsoid (requires Proj extension).

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
function perimeter(geom, ::Type{T} = Float64; threaded=False(), init = zero(T), kwargs...) where T
    perimeter(AutoManifold(), geom, T; threaded, init, kwargs...)
end

function perimeter(::AutoManifold, geom, ::Type{T} = Float64; kwargs...) where T
    m, scales = _auto_manifold(geom)
    _perimeter_auto(m, geom, T, scales; kwargs...)
end

# The Proj extension adds a `Geodesic` method that applies the CRS axis scales.
_perimeter_auto(m::Manifold, geom, ::Type{T}, scales; kwargs...) where T = perimeter(m, geom, T; kwargs...)

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