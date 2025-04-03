# module UnitSpherical

using CoordinateTransformations
using StaticArrays
import GeoInterface as GI, GeoFormatTypes as GFT
using LinearAlgebra

"""
    UnitSphericalPoint(v)

A unit spherical point, i.e., point living on the 2-sphere (𝕊²),
represented as Cartesian coordinates in ℝ³.

This currently has no support for heights, only going from lat long to spherical
and back again.
"""
struct UnitSphericalPoint{T} <: StaticArrays.StaticVector{3, T}
    data::SVector{3, T}
    UnitSphericalPoint{T}(v::SVector{3, T}) where T = new{T}(v)
    UnitSphericalPoint(x::T, y::T, z::T) where T = new{T}(SVector{3, T}((x, y, z)))
end


function UnitSphericalPoint(v::AbstractVector{T}) where T
    if length(v) == 3
        UnitSphericalPoint{T}(SVector(v[1], v[2], v[3]))
    elseif length(v) == 2
        UnitSphereFromGeographic()(v)
    else
        throw(ArgumentError("Passed a vector of length `$(length(v))` to `UnitSphericalPoint` constructor, which only accepts vectors of length 3 (assumed to be on the unit sphere) or 2 (assumed to be geographic lat/long)."))
    end
end

UnitSphericalPoint(v) = UnitSphericalPoint(GI.trait(v), v)
UnitSphericalPoint(::GI.PointTrait, v) = UnitSphereFromGeographic()(v) # since it works on any GI compatible point

# StaticArraysCore.jl interface implementation
Base.Tuple(p::UnitSphericalPoint) = Base.Tuple(p.data)
Base.@propagate_inbounds @inline Base.getindex(p::UnitSphericalPoint, i::Int64) = p.data[i]
Base.setindex(p::UnitSphericalPoint, args...) = throw(ArgumentError("`setindex!` on a UnitSphericalPoint is not permitted as it is static."))
@generated function StaticArrays.similar_type(::Type{SV}, ::Type{T},
    s::StaticArrays.Size{S}) where {SV <: UnitSphericalPoint,T,S}
    return if length(S) === 1
        UnitSphericalPoint{T}
    else
        StaticArrays.default_similar_type(T, s(), Val{length(S)})
    end
end

# Base math implementation (really just forcing re-wrapping)
Base.:(*)(a::UnitSphericalPoint, b::UnitSphericalPoint) = a .* b
function Base.broadcasted(f, a::AbstractArray{T}, b::UnitSphericalPoint) where {T <: UnitSphericalPoint}
    return Base.broadcasted(f, a, (b,))
end

# GeoInterface implementation: traits
GI.trait(::UnitSphericalPoint) = GI.PointTrait()
GI.geomtrait(::UnitSphericalPoint) = GI.PointTrait()
# GeoInterface implementation: coordinate traits
GI.is3d(::GI.PointTrait, ::UnitSphericalPoint) = true
GI.ismeasured(::GI.PointTrait, ::UnitSphericalPoint) = false
# GeoInterface implementation: accessors
GI.ncoord(::GI.PointTrait, ::UnitSphericalPoint) = 3
GI.getcoord(::GI.PointTrait, p::UnitSphericalPoint) = p[i]
# GeoInterface implementation: metadata (CRS, extent, etc)
GI.crs(::UnitSphericalPoint) = GFT.ProjString("+proj=cart +R=1 +type=crs") # TODO: make this a full WKT definition

# several useful LinearAlgebra functions, forwarded to the static arrays
LinearAlgebra.cross(p1::UnitSphericalPoint, p2::UnitSphericalPoint) = UnitSphericalPoint(LinearAlgebra.cross(p1.data, p2.data))
LinearAlgebra.dot(p1::UnitSphericalPoint, p2::UnitSphericalPoint) = LinearAlgebra.dot(p1.data, p2.data)
LinearAlgebra.normalize(p1::UnitSphericalPoint) = UnitSphericalPoint(LinearAlgebra.normalize(p1.data))

# Spherical cap implementation
struct SphericalCap{T}
    point::UnitSphericalPoint{T}
    radius::T
end

SphericalCap(point::UnitSphericalPoint{T}, radius::Number) where T = SphericalCap{T}(point, convert(T, radius))
SphericalCap(point, radius::Number) = SphericalCap(GI.trait(point), point, radius)
function SphericalCap(::GI.PointTrait, point, radius::Number)
    return SphericalCap(UnitSphereFromGeographic()(point), radius)
end

SphericalCap(geom) = SphericalCap(GI.trait(geom), geom)
SphericalCap(t::GI.PointTrait, geom) = SphericalCap(t, geom, 0)
# TODO: add implementations for line string and polygon traits
# TODO: add implementations to merge two spherical caps
# TODO: add implementations for multitraits based on this

# TODO: this returns an approximately antipodal point...

spherical_distance(x::UnitSphericalPoint, y::UnitSphericalPoint) = acos(clamp(x ⋅ y, -1.0, 1.0))

# TODO: exact-predicate intersection
# This is all inexact and thus subject to floating point error
function _intersects(x::SphericalCap, y::SphericalCap)
    spherical_distance(x.point, y.point) <= max(x.radius, y.radius)
end

_disjoint(x::SphericalCap, y::SphericalCap) = !_intersects(x, y)

function _contains(big::SphericalCap, small::SphericalCap)
    dist = spherical_distance(big.point, small.point)
    return dist < big.radius #=small.point in big=# && dist + small.radius < big.radius
end

function slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01::Number)
    Ω = spherical_distance(a, b)
    sinΩ = sin(Ω)
    return (sin((1-i01)*Ω) / sinΩ) * a + (sin(i01*Ω)/sinΩ) * b
end

function slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01s::AbstractVector{<: Number})
    Ω = spherical_distance(a, b)
    sinΩ = sin(Ω)
    return (sin((1-i01)*Ω) / sinΩ) .* a .+ (sin(i01*Ω)/sinΩ) .* b
end




function circumcenter_on_unit_sphere(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    LinearAlgebra.normalize(a × b + b × c + c × a)
end

"Get the circumcenter of the triangle (a, b, c) on the unit sphere.  Returns a normalized 3-vector."
function SphericalCap(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    circumcenter = circumcenter_on_unit_sphere(a, b, c)
    circumradius = spherical_distance(a, circumcenter)
    return SphericalCap(circumcenter, circumradius)
end

function _is_ccw_unit_sphere(v_0,v_c,v_i)
    # checks if the smaller interior angle for the great circles connecting u-v and v-w is CCW
    return(LinearAlgebra.dot(LinearAlgebra.cross(v_c - v_0,v_i - v_c), v_i) < 0)
end

function angle_between(a, b, c)
    ab = b - a
    bc = c - b
    norm_dot = (ab ⋅ bc) / (LinearAlgebra.norm(ab) * LinearAlgebra.norm(bc))
    angle =  acos(clamp(norm_dot, -1.0, 1.0))
    if _is_ccw_unit_sphere(a, b, c)
        return angle
    else
        return 2π - angle
    end
end


# Coordinate transformations from lat/long to geographic and back
struct UnitSphereFromGeographic <: CoordinateTransformations.Transformation 
end

function (::UnitSphereFromGeographic)(geographic_point)
    # Asssume that geographic_point is GeoInterface compatible
    # Longitude is directly translatable to a spherical coordinate
    # θ (azimuth)
    θ = GI.x(geographic_point)
    # The polar angle is 90 degrees minus the latitude
    # ϕ (polar angle)
    ϕ = 90 - GI.y(geographic_point)
    # Since this is the unit sphere, the radius is assumed to be 1,
    # and we don't need to multiply by it.
    sinϕ, cosϕ = sincosd(ϕ)
    sinθ, cosθ = sincosd(θ)

    return UnitSphericalPoint(
        sinϕ * cosθ,
        sinϕ * sinθ,
        cosϕ
    )
end

struct GeographicFromUnitSphere <: CoordinateTransformations.Transformation 
end

function (::GeographicFromUnitSphere)(xyz::AbstractVector)
    @assert length(xyz) == 3 "GeographicFromUnitCartesian expects a 3D Cartesian vector"
    x, y, z = xyz.data
    return (
        atan(y, x) |> rad2deg,
        90 - (atan(hypot(x, y), z) |> rad2deg),
    )
end

function randsphericalangles(n)
    θ = 2π .* rand(n)
    ϕ = acos.(2 .* rand(n) .- 1)
    return tuple.(θ, ϕ)
end

"""
    randsphere(n)

Return `n` random [`UnitSphericalPoint`](@ref)s spanning the whole sphere 𝕊².
"""
function randsphere(n)
    θϕs = randsphericalangles(n)
    return map(θϕs) do θϕ
        θ, ϕ = θϕ
        sinθ, cosθ = sincos(θ)
        sinϕ, cosϕ = sincos(ϕ)
        UnitSphericalPoint(
            sinϕ * cosθ,
            sinϕ * sinθ,
            cosϕ
        )
    end
end


# end