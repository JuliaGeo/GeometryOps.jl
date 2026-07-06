# # Spherical arc extents

#=
```@docs; canonical=false
spherical_arc_extent
```

A great-circle arc bulges away from the chord between its endpoints, so the
endpoints' bounding box does not in general contain the arc: two points at
`z = 0.9` on either side of the prime meridian are joined by an arc whose
midpoint has `z = 0.9 / cos(θ/2) > 0.9`.  [`spherical_arc_extent`](@ref)
computes a box that contains the whole arc.

With `t̂ₐ` the unit tangent at `a` pointing along the arc, the arc is
`p(φ) = a cos(φ) + t̂ₐ sin(φ)` for `φ ∈ [0, θ]`, so each Cartesian
coordinate is a sinusoid `pᵢ(φ) = Rᵢ cos(φ - φᵢ)` with amplitude
`Rᵢ = hypot(aᵢ, t̂ₐᵢ)`.  Since `θ ≤ π` there is at most one interior maximum
and one interior minimum per axis: a maximum exists iff `pᵢ` increases at
`a` and decreases at `b` (`t̂ₐᵢ > 0 > t̂ᵦᵢ`), where it attains `Rᵢ`; minima
mirror.  The tangents come from [`robust_cross_product`](@ref), which keeps
nearly-degenerate and nearly-antipodal arcs stable.  Bounds are padded by a
few ulps to absorb floating point error, as S2's `S2LatLngRectBounder` pads
by its maximum error.
=#

"""
    spherical_arc_extent(a, b)::Extents.Extent{(:X, :Y, :Z)}

The 3D Cartesian extent of the shorter great-circle arc between `a` and `b`
on the unit sphere.  Accepts `UnitSphericalPoint`s, or any GeoInterface
point (interpreted geographically, as longitude/latitude, like the
`UnitSphericalPoint` constructor itself).

The extent is exact up to floating point error, padded by a few ulps so it
always contains the arc.  For antipodal endpoints the arc's plane is
ambiguous; the one chosen by [`robust_cross_product`](@ref) is used.

## Example

```jldoctest
julia> using GeometryOps.UnitSpherical

julia> ext = spherical_arc_extent(UnitSphericalPoint(1, 0, 0), UnitSphericalPoint(0, 1, 0));

julia> ext.X[2] ≈ 1 && ext.Y[2] ≈ 1
true
```
"""
spherical_arc_extent(a, b) = spherical_arc_extent(UnitSphericalPoint(a), UnitSphericalPoint(b))
function spherical_arc_extent(a::UnitSphericalPoint{T1}, b::UnitSphericalPoint{T2}) where {T1, T2}
    F = float(promote_type(T1, T2))
    pad = 4 * eps(F)
    if a == b
        bounds = ntuple(i -> (F(a[i]) - pad, F(a[i]) + pad), 3)
    else
        n = robust_cross_product(a, b)
        ta = normalize(cross(n, a))     # unit tangent at `a`, pointing along the arc
        tb = normalize(cross(n, b))     # unit tangent at `b`, pointing along the arc
        bounds = ntuple(3) do i
            lo, hi = minmax(F(a[i]), F(b[i]))
            ta[i] > 0 > tb[i] && (hi = hypot(F(a[i]), F(ta[i])))
            ta[i] < 0 < tb[i] && (lo = -hypot(F(a[i]), F(ta[i])))
            (lo - pad, hi + pad)
        end
    end
    return Extents.Extent(X = bounds[1], Y = bounds[2], Z = bounds[3])
end
