# # Intersection area
export intersection_area

#=
`intersection_area(alg, a, b)` is `area(manifold(alg), intersection(alg, a, b))` with the
result geometry left unbuilt: each engine stops at the rings it would have wrapped and
accumulates their area in place.

Two engines implement it, on both `Planar()` and `Spherical()`:

- [`ConvexConvexSutherlandHodgman`](@ref) clips into its cache buffers and sums straight
  off them, so with a `cache` the whole call is allocation-free.
- [`OverlayNG`](@ref) runs the arrangement as usual but replaces polygon assembly with a
  ring-area sum, skipping the `Polygon`/`MultiPolygon` wrappers.
=#

"""
    intersection_area(alg::Algorithm, geom_a, geom_b, [T = Float64]; kwargs...)

Area of the intersection of `geom_a` and `geom_b`, computed without constructing the
intersection geometry.

Equivalent to `area(manifold(alg), intersection(alg, geom_a, geom_b))`, but the result
polygon is never built. Supported for [`ConvexConvexSutherlandHodgman`](@ref) (which takes
the same `cache` keyword as `intersection`) and [`OverlayNG`](@ref), on `Planar()` and
`Spherical()` manifolds. On `Spherical()` the result is in square units of the manifold
radius.

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

# ## Ring area
#
# Both engines hand their rings over as a plain vector of points, closed or not.

_area_scale(::Planar) = 1
_area_scale(m::Spherical) = m.radius^2

# Shoelace over the ring, wrapping from the last point to the first. A repeated closing
# point contributes a zero term, so closed and open rings both work.
function _ring_area(::Planar, pts::AbstractVector, ::Type{T}) where {T}
    n = length(pts)
    n < 3 && return zero(T)
    a = zero(T)
    for i in 1:n
        p1, p2 = pts[i], pts[mod1(i + 1, n)]
        a += _area_component(p1, p2)
    end
    return T(a / 2)
end

@inline _unit_spherical(p::UnitSphericalPoint) = p
@inline _unit_spherical(p) = UnitSphericalPoint(GI.PointTrait(), p)

# Signed unit-sphere area by the same fan triangulation `area(::Spherical, geom)` uses.
function _ring_area(::Spherical, pts::AbstractVector, ::Type{T}) where {T}
    n = length(pts)
    n < 3 && return zero(T)
    p1 = _unit_spherical(pts[1])
    # drop the closing point, if the ring carries one
    _unit_spherical(pts[n]) ≈ p1 && (n -= 1)
    n < 3 && return zero(T)
    a = zero(T)
    for i in 2:(n - 1)
        a += _spherical_triangle_area(Eriksson(), p1, _unit_spherical(pts[i]), _unit_spherical(pts[i + 1]))
    end
    return T(a)
end

# ## Sutherland-Hodgman
#
# The clip buffers are the result; there is nothing left to do but measure them.

function intersection_area(
    alg::ConvexConvexSutherlandHodgman{Planar}, geom_a, geom_b, ::Type{T}=Float64;
    cache::Union{Nothing, SutherlandHodgmanCache}=nothing
) where {T<:AbstractFloat}
    cache = isnothing(cache) ? SutherlandHodgmanCache(Planar(), T) : _sh_check_cache(cache, Tuple{T,T})
    return abs(_ring_area(Planar(), _sh_clip_planar!(cache, geom_a, geom_b, T), T))
end

function intersection_area(
    alg::ConvexConvexSutherlandHodgman{Spherical{F}}, geom_a, geom_b, ::Type{T}=Float64;
    cache::Union{Nothing, SutherlandHodgmanCache}=nothing
) where {F, T<:AbstractFloat}
    m = alg.manifold
    cache = isnothing(cache) ? SutherlandHodgmanCache(m, T) : _sh_check_cache(cache, UnitSphericalPoint{T})
    pts = _sh_clip_spherical!(cache, geom_a, geom_b, T)
    return abs(_ring_area(m, pts, T)) * _area_scale(m)
end

# ## OverlayNG
#
# The arrangement runs unchanged; only the extraction is replaced. Anything but two areal
# inputs has an intersection of dimension < 2, hence zero area, and needs no noding at all.

function intersection_area(alg::OverlayNG, geom_a, geom_b, ::Type{T}=Float64) where {T<:AbstractFloat}
    m = alg.manifold
    (_overlay_dimension(geom_a) == 2 && _overlay_dimension(geom_b) == 2) || return zero(T)
    (_ov_isempty(geom_a) || _ov_isempty(geom_b)) && return zero(T)

    input = _OverlayInput(m, geom_a, geom_b, 2, 2, alg.exact, false, false, nothing, nothing)
    ea, eb = _overlay_envelopes(m, OVERLAY_INTERSECTION, input)
    _empty_result_short_circuit(m, OVERLAY_INTERSECTION, input, ea, eb) && return zero(T)

    g = _overlay_marked_graph(m, alg.point_type, OVERLAY_INTERSECTION, geom_a, geom_b,
                              input, alg.exact, nothing, nothing, ea, eb)
    ctx = _build_polygon_ctx(m, g, graph_result_area_edges(g); exact = alg.exact)

    total = zero(T)
    for shell_handle in ctx.shell_list
        shell = ctx.edge_rings[shell_handle]
        a = abs(_ring_area(m, shell.ring_pts, T))
        for hole in shell.holes
            a -= abs(_ring_area(m, ctx.edge_rings[hole].ring_pts, T))
        end
        total += a
    end
    return total * _area_scale(m)
end
