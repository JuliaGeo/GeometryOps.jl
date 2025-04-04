#=


# UnitSphericalPoint

This file defines the [`UnitSphericalPoint`](@ref) type, which is 
a three-dimensional Cartesian point on the unit 2-sphere (i.e., of radius 1).

This file contains the full implementation of the type as well as a `spherical_distance` function
that computes the great-circle distance between two points on the unit sphere.

```@docs; canonical=false
UnitSphericalPoint
spherical_distance
```
=#

# ## Type definition and constructors
"""
    UnitSphericalPoint(v)

A unit spherical point, i.e., point living on the 2-sphere (𝕊²),
represented as Cartesian coordinates in ℝ³.

This currently has no support for heights, only going from lat long to spherical
and back again.

## Examples

```jldoctest
julia> UnitSphericalPoint(1, 0, 0)
UnitSphericalPoint(1.0, 0.0, 0.0)
```

"""
struct UnitSphericalPoint{T} <: StaticArrays.FieldVector{3, T}
    x::T
    y::T
    z::T
end

UnitSphericalPoint{T}(v::SVector{3, T}) where T = UnitSphericalPoint{T}(v...)
UnitSphericalPoint(v::NTuple{3, T}) where T = UnitSphericalPoint{T}(v...)
UnitSphericalPoint{T}(v::NTuple{3, T}) where T = UnitSphericalPoint{T}(v...)

UnitSphericalPoint(v::SVector{3, T}) where T = UnitSphericalPoint{T}(v...)
## handle the 2-tuple case specifically
UnitSphericalPoint(v::NTuple{2, T}) where T = UnitSphereFromGeographic()(v)
## handle the GeoInterface case, this is the catch-all method
UnitSphericalPoint(v) = UnitSphericalPoint(GI.trait(v), v)
UnitSphericalPoint(::GI.PointTrait, v) = UnitSphereFromGeographic()(v) # since it works on any GI compatible point
## finally, handle the case where a vector is passed in
## we may want it to go to the geographic pipeline _or_ direct materialization
Base.@propagate_inbounds function UnitSphericalPoint(v::AbstractVector{T}) where T
    if length(v) == 3
        UnitSphericalPoint{T}(v[1], v[2], v[3])
    elseif length(v) == 2
        UnitSphereFromGeographic()(v)
    else
        @boundscheck begin
            throw(ArgumentError("""
            Passed a vector of length `$(length(v))` to the `UnitSphericalPoint` constructor, 
            which only accepts vectors of lengths: 
            - **3** (assumed to be on the unit sphere) 
            - **2** (assumed to be geographic lat/long)
            """))
        end
    end
end

Base.show(io::IO, p::UnitSphericalPoint) = print(io, "UnitSphericalPoint($(p.x), $(p.y), $(p.z))")

# ## Interface implementations

# StaticArraysCore.jl interface implementation
Base.setindex(p::UnitSphericalPoint, args...) = throw(ArgumentError("`setindex!` on a UnitSphericalPoint is not permitted as it is static."))
StaticArrays.similar_type(::Type{<: UnitSphericalPoint}, ::Type{Eltype}, ::Size{(3,)}) where Eltype = UnitSphericalPoint{Eltype}
# Base math implementation (really just forcing re-wrapping)
# Base.:(*)(a::UnitSphericalPoint, b::UnitSphericalPoint) = a .* b
function Base.broadcasted(f, a::AbstractArray{T}, b::UnitSphericalPoint) where {T <: UnitSphericalPoint}
    return Base.broadcasted(f, a, (b,))
end
Base.isnan(p::UnitSphericalPoint) = any(isnan, p)
Base.isinf(p::UnitSphericalPoint) = any(isinf, p)
Base.isfinite(p::UnitSphericalPoint) = all(isfinite, p)


# GeoInterface implementation
## Traits:
GI.trait(::UnitSphericalPoint) = GI.PointTrait()
GI.geomtrait(::UnitSphericalPoint) = GI.PointTrait()
## Coordinate traits:
GI.is3d(::GI.PointTrait, ::UnitSphericalPoint) = true
GI.ismeasured(::GI.PointTrait, ::UnitSphericalPoint) = false
## Accessors:
GI.ncoord(::GI.PointTrait, ::UnitSphericalPoint) = 3
GI.getcoord(::GI.PointTrait, p::UnitSphericalPoint) = p[i]
## Metadata (CRS, extent, etc)
GI.crs(::UnitSphericalPoint) = GFT.ProjString("+proj=cart +R=1 +type=crs") # TODO: make this a full WKT definition
# TODO: extent is a little tricky - do we do a spherical cap or an Extents.Extent?

# ## Spherical distance
"""
    spherical_distance(x::UnitSphericalPoint, y::UnitSphericalPoint)

Compute the spherical distance between two points on the unit sphere.  
Returns a `Number`, usually Float64 but that depends on the input type.

# Extended help

## Doctests

```jldoctest
julia> spherical_distance(UnitSphericalPoint(1, 0, 0), UnitSphericalPoint(0, 1, 0))
1.5707963267948966
```

```jldoctest
julia> spherical_distance(UnitSphericalPoint(1, 0, 0), UnitSphericalPoint(1, 0, 0))
0.0
```
"""
spherical_distance(x::UnitSphericalPoint, y::UnitSphericalPoint) = acos(clamp(x ⋅ y, -1.0, 1.0))

# ## Random points
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{UnitSphericalPoint}) = rand(rng, UnitSphericalPoint{Float64})
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{UnitSphericalPoint{T}}) where T <: Number
    ϕ = 2π * rand(rng, T)
    θ = acos(2 * rand(rng, T) - 1)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    return UnitSphericalPoint(
        sinϕ * cosθ,
        sinϕ * sinθ,
        cosϕ
    )
end

# ## Tests

# @testitem "UnitSphericalPoint constructor" begin
#     using GeometryOps.UnitSpherical
#     import GeoInterface as GI

#     northpole = UnitSphericalPoint{Float64}(1, 0, 0)
#     # test that the constructor works for a vector of length 3
#     @test UnitSphericalPoint((1, 0, 0)) == northpole
#     @test UnitSphericalPoint(SVector(1, 0, 0)) == northpole
#     @test UnitSphericalPoint([1, 0, 0]) == northpole
#     # test that the constructor works for a tuple of length 2
#     # and interprets such a thing as a geographic point
#     @test UnitSphericalPoint((90, 0)) == northpole
#     @test UnitSphericalPoint([90, 0]) == northpole
#     @test UnitSphericalPoint(GI.Point((90, 0))) == northpole


# end

#=
```@meta
CollapsedDocStrings = true
```
=#