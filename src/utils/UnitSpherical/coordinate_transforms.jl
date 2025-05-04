#=
# Coordinate transformations
=#
# Coordinate transformations from lat/long to geographic and back
"""
    UnitSphereFromGeographic()

A transformation that converts a geographic point (latitude, longitude) to a 
[`UnitSphericalPoint`] in ℝ³.

Accepts any [GeoInterface-compatible](https://github.com/JuliaGeo/GeoInterface.jl) point.

## Examples

```jldoctest
julia> import GeoInterface as GI; using GeometryOps.UnitSpherical

julia> UnitSphereFromGeographic()(GI.Point(45, 45))
3-element UnitSphericalPoint{Float64} with indices SOneTo(3):
 0.5000000000000001
 0.5000000000000001
 0.7071067811865476
```

```jldoctest
julia> using GeometryOps.UnitSpherical

julia> UnitSphereFromGeographic()((45, 45))
3-element UnitSphericalPoint{Float64} with indices SOneTo(3):
 0.5000000000000001
 0.5000000000000001
 0.7071067811865476
```
"""
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

"""
    GeographicFromUnitSphere()

A transformation that converts a [`UnitSphericalPoint`](@ref) in ℝ³ to a 
2-tuple geographic point (longitude, latitude), in degrees.

Accepts any 3-element vector, but the input is assumed to be on the unit sphere.

## Examples

```jldoctest
julia> using GeometryOps.UnitSpherical

julia> GeographicFromUnitSphere()(UnitSphericalPoint(0.5, 0.5, 1/√(2)))
(45.0, 44.99999999999999)
```
(the inaccuracy is due to the precision of the `atan` function)

"""
struct GeographicFromUnitSphere <: CoordinateTransformations.Transformation 
end

function (::GeographicFromUnitSphere)(xyz::AbstractVector)
    @assert length(xyz) == 3 "GeographicFromUnitCartesian expects a 3D Cartesian vector"
    x, y, z = xyz
    return (
        atand(y, x),
        asind(z),
    )
end