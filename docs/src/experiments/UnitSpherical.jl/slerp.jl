#=

# Slerp (spherical linear interpolation)

```@docs; canonical=false
slerp
```

Slerp is a spherical interpolation method that is used to interpolate between two points on a unit sphere.
It is a generalization of linear interpolation to the sphere.

The algorithm takes two spherical points and a parameter, `i01`, that is a number between 0 and 1.
The algorithm returns a point on the unit sphere that is a linear interpolation between the two points.

The way this works, is that it basically takes the great circle path between the two points and then
interpolates along that path.

=#

"""
    slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01::Number)

Interpolate between `a` and `b`, at a proportion `i01` 
between 0 and 1 along the path from `a` to `b`.

## Examples

```jldoctest
julia> slerp(UnitSphericalPoint(1, 0, 0), UnitSphericalPoint(0, 1, 0), 0.5)
3-element UnitSphericalPoint{Float64} with indices SOneTo(3):
 0.7071067811865475
 0.7071067811865475
 0.0
```
"""
function slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01::Number)
    Ω = spherical_distance(a, b)
    sinΩ = sin(Ω)
    return (sin((1-i01)*Ω) / sinΩ) * a + (sin(i01*Ω)/sinΩ) * b
end

function slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01s::AbstractVector{<: Number})
    Ω = spherical_distance(a, b)
    sinΩ = sin(Ω)
    return @. (sin((1 - i01s) * Ω) / sinΩ) * a + (sin(i01s * Ω) / sinΩ) * b
end

#=
```@meta
CollapsedDocStrings = true
```
=#