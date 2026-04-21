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

Uses the tangent-vector form `cos(r)·a + sin(r)·dir` — where `r = i01 ·
spherical_distance(a, b)` and `dir = normalize(robust_cross_product(a, b) × a)`
is the unit tangent at `a` pointing toward `b`. This avoids the `1/sin(Ω)`
divisor of the classic `sin((1-t)Ω)/sin(Ω) · a + sin(tΩ)/sin(Ω) · b`
formulation, which collapses for near- and exactly-antipodal inputs.
Adapted from Google's S2 geometry library (see [`S2::Interpolate`](https://github.com/google/s2geometry/blob/a4f0cf58a9cfc214585c39de6e3682384fac0917/src/s2/s2edge_distances.cc#L77)
and [`S2::GetPointOnLine`](https://github.com/google/s2geometry/blob/a4f0cf58a9cfc214585c39de6e3682384fac0917/src/s2/s2edge_distances.cc#L47)).

For exactly antipodal `a` and `b` the great circle is mathematically
ambiguous; `robust_cross_product` returns a deterministic perpendicular via
its symbolic-perturbation branch, so the result is still a well-defined unit
vector on *some* great circle through both points.

## Examples

```jldoctest
julia> using GeometryOps.UnitSpherical

julia> slerp(UnitSphericalPoint(1, 0, 0), UnitSphericalPoint(0, 1, 0), 0.5)
3-element UnitSphericalPoint{Float64} with indices SOneTo(3):
 0.7071067811865476
 0.7071067811865475
 0.0
```
"""
function slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01::Number)
    i01 == 0 && return a
    i01 == 1 && return b
    a == b && return a
    Ω = spherical_distance(a, b)
    dir = normalize(cross(robust_cross_product(a, b), a))
    r = i01 * Ω
    return normalize(cos(r) * a + sin(r) * dir)
end

function slerp(a::UnitSphericalPoint, b::UnitSphericalPoint, i01s::AbstractVector{<: Number})
    a == b && return fill(a, size(i01s))
    Ω = spherical_distance(a, b)
    dir = normalize(cross(robust_cross_product(a, b), a))
    return [begin
        t == 0 ? a :
        t == 1 ? b :
        normalize(cos(t * Ω) * a + sin(t * Ω) * dir)
    end for t in i01s]
end

#=
```@meta
CollapsedDocStrings = true
```
=#