# # Intersection area
export intersection_area

#=
`intersection_area(alg, a, b)` is `area(manifold(alg), intersection(alg, a, b))` with the
result geometry left unbuilt: each engine stops at the rings it would have wrapped and
sums their area in place, through the shared `_ring_area` kernels in `methods/area.jl`.

The per-engine work lives with its engine — `_sh_clip_planar!` / `_sh_clip_spherical!` in
`sutherland_hodgman.jl`, `_overlay_intersection_area` in `overlayng/overlay_ng.jl` — so
this file is only the public surface.
=#

"""
    intersection_area(alg::Algorithm, geom_a, geom_b, [T = Float64]; kwargs...)

Area of the intersection of `geom_a` and `geom_b`, computed without constructing the
intersection geometry.

Equivalent to `area(manifold(alg), intersection(alg, geom_a, geom_b))`, but the result
polygon is never built. On `Spherical()` the result is in square units of the manifold
radius.

The algorithm is required, and the supported list is closed:

| algorithm | manifolds | notes |
|:----------|:----------|:------|
| [`ConvexConvexSutherlandHodgman`](@ref) | `Planar()`, `Spherical()` | takes the same `cache` keyword as `intersection`; with one, the call allocates nothing |
| [`OverlayNG`](@ref) | `Planar()`, `Spherical()` | `T` sets the accumulator and return type only — the arrangement always runs at the algorithm's `point_type` |

## Example

```julia
import GeometryOps as GO, GeoInterface as GI

a = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
b = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

GO.intersection_area(GO.OverlayNG(), a, b)  # 1.0

alg = GO.ConvexConvexSutherlandHodgman()
cache = GO.SutherlandHodgmanCache(alg)
GO.intersection_area(alg, a, b; cache)      # 1.0, allocation-free
```
"""
function intersection_area end

# ## Sutherland-Hodgman
#
# The clip buffers are the result; there is nothing left to do but measure them.

function intersection_area(
    alg::ConvexConvexSutherlandHodgman{Planar}, geom_a, geom_b, ::Type{T}=Float64;
    cache::Union{Nothing, SutherlandHodgmanCache}=nothing
) where {T<:AbstractFloat}
    _sh_check_polygon_traits(GI.trait(geom_a), GI.trait(geom_b))
    cache = isnothing(cache) ? SutherlandHodgmanCache(Planar(), T) : _sh_check_cache(cache, Tuple{T,T})
    return abs(_ring_area(Planar(), _sh_clip_planar!(cache, geom_a, geom_b, T), T))
end

function intersection_area(
    alg::ConvexConvexSutherlandHodgman{Spherical{F}}, geom_a, geom_b, ::Type{T}=Float64;
    cache::Union{Nothing, SutherlandHodgmanCache}=nothing
) where {F, T<:AbstractFloat}
    _sh_check_polygon_traits(GI.trait(geom_a), GI.trait(geom_b))
    m = alg.manifold
    cache = isnothing(cache) ? SutherlandHodgmanCache(m, T) : _sh_check_cache(cache, UnitSphericalPoint{T})
    pts = _sh_clip_spherical!(cache, geom_a, geom_b, T)
    return T(abs(_ring_area(m, pts, T)) * _area_scale(m))
end

# ## OverlayNG

intersection_area(alg::OverlayNG, geom_a, geom_b, ::Type{T}=Float64) where {T<:AbstractFloat} =
    T(_overlay_intersection_area(alg.manifold, T, alg.point_type, geom_a, geom_b, alg.exact))
