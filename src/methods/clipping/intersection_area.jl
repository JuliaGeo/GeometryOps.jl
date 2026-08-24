# # Intersection area
export intersection_area

#=
`intersection_area(alg, a, b)` is `area(manifold(alg), intersection(alg, a, b))` with the
result geometry left unbuilt: each engine stops at the rings it would have wrapped and
sums their area in place, through the shared `_ring_area` kernels in `methods/area.jl`.

The per-engine work lives with its engine — `_sh_clip_planar!` / `_sh_clip_spherical!` in
`sutherland_hodgman.jl`, `_overlay_intersection_area` in `overlayng/overlay_ng.jl`,
`_RingMeasurer` in `clipping_processor.jl` — so this file is only the public surface.
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
| [`FosterHormannClipping`](@ref) | as `intersection` | accumulates during the trace for hole-free polygon pairs; other inputs delegate to the polygon path. Takes a [`FosterHormannCache`](@ref) as `cache` on the traced path, which removes the per-call allocation of its vertex lists |

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
    cache::SutherlandHodgmanCache = SutherlandHodgmanCache(Planar(), T)
) where {T<:AbstractFloat}
    _sh_check_polygon_traits(GI.trait(geom_a), GI.trait(geom_b))
    _sh_check_cache(cache, Tuple{T,T})
    #-- the clip buffers are open: they never repeat the first vertex
    return abs(_ring_area(Planar(), _sh_clip_planar!(cache, geom_a, geom_b, T), T; closed = false))
end

function intersection_area(
    alg::ConvexConvexSutherlandHodgman{Spherical{F}}, geom_a, geom_b, ::Type{T}=Float64;
    cache::SutherlandHodgmanCache = SutherlandHodgmanCache(alg.manifold, T)
) where {F, T<:AbstractFloat}
    _sh_check_polygon_traits(GI.trait(geom_a), GI.trait(geom_b))
    _sh_check_cache(cache, UnitSphericalPoint{T})
    m = alg.manifold
    pts = _sh_clip_spherical!(cache, geom_a, geom_b, T)
    return T(abs(_ring_area(m, pts, T; closed = false)) * _area_scale(m))
end

# ## OverlayNG

intersection_area(alg::OverlayNG, geom_a, geom_b, ::Type{T}=Float64) where {T<:AbstractFloat} =
    T(_overlay_intersection_area(alg.manifold, T, alg.point_type, geom_a, geom_b, alg.exact))

# ## Foster-Hormann
#
# The tracer walks the result ring one vertex at a time, so the area can be accumulated as
# it goes and no ring is ever built (`_RingMeasurer` in clipping_processor.jl).
#
# That covers hole-free polygon pairs. Holes and multipolygons are assembled from the
# traced rings *afterwards* — `_add_holes_to_polys!` cuts holes into result polygons and
# can split, merge and drop them — and there is no ring-level shortcut through that, so
# those inputs delegate to the polygon path. The answer is the same either way; only the
# saving is lost.

#-- `cache` is offered on this path only. The polygon fallback below reaches
#-- `_add_holes_to_polys!`, which clips again from inside the clip it is finishing; handing
#-- that nested call the same buffers would have it overwrite the lists its caller is still
#-- reading. The fast path makes no such nested call, and it is the one a regridder runs.
function intersection_area(
    alg::FosterHormannClipping, geom_a, geom_b, ::Type{T}=Float64;
    cache::Union{Nothing, FosterHormannCache} = nothing, kwargs...
) where {T<:AbstractFloat}
    m = alg.manifold
    if !(GI.trait(geom_a) isa GI.PolygonTrait && GI.trait(geom_b) isa GI.PolygonTrait) ||
       GI.nhole(geom_a) != 0 || GI.nhole(geom_b) != 0
        return _fh_polygon_path_area(alg, m, geom_a, geom_b, T; kwargs...)
    end
    cache === nothing || _fh_check_cache(cache, T, _fh_point_type(m, T))
    ext_a, ext_b = GI.getexterior(geom_a), GI.getexterior(geom_b)
    a_list, b_list, a_idx_list =
        _build_ab_list(alg, T, ext_a, ext_b, _inter_delay_cross_f, _inter_delay_bounce_f; exact = True(), cache)
    sink = _trace_polynodes!(_RingMeasurer(m, T), alg, T, a_list, b_list, a_idx_list,
                             _inter_step, geom_a, geom_b)
    #-- no crossings at all: either one polygon contains the other, or they are disjoint.
    #-- Rare, and the polygon path already decides it in one step.
    sink.nrings == 0 && return _fh_polygon_path_area(alg, m, geom_a, geom_b, T; kwargs...)
    return T(sink.area * _area_scale(m))
end

#-- `target` is spelled out because the areal components are the only ones with an area,
#-- and because `intersection`'s own `target = nothing` default is currently broken
#-- (`TraitTarget(nothing)` has no method).
_fh_polygon_path_area(alg, m, geom_a, geom_b, ::Type{T}; kwargs...) where {T} =
    T(area(m, intersection(alg, geom_a, geom_b, T; target = GI.PolygonTrait(), kwargs...)))
