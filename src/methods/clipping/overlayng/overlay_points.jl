# # OverlayPoints — point × point overlay (port of JTS `OverlayPoints`)
#
# Phase 3 of the OverlayNG port (design doc §4). A port of
# `operation/overlayng/OverlayPoints.java`: the overlay of two inputs which are
# *both* point geometries. There is no arrangement to build — the answer is a
# set operation over the two coordinate sets — so this path bypasses the noder,
# graph, labeller and builders entirely, exactly as JTS does.
#
# Semantics (JTS's, minus the precision model):
#
#   * points with identical coordinates are merged to a single point;
#   * the first occurrence of a coordinate supplies the output coordinate;
#   * an empty result is the empty dimension-0 geometry.
#
# Deviations from the Java, each deliberate:
#
#   * **No precision model.** `copyPoint`/`roundCoord` exist only to round to a
#     fixed `PrecisionModel`; there is none here by design (§0), so they are
#     dropped and coordinates pass through bit-exact.
#   * **Point identity is the manifold's kernel point**, not the raw
#     coordinate. On `Planar` that is the same 2-tuple (signed zeros
#     normalized); on `Spherical` it is the normalized unit vector, so
#     longitude ±180° and the pole's arbitrary longitude identify correctly.
#     This is the same identity the noder's vertex nodes use, so point overlay
#     and edge overlay agree on what "the same point" means.
#   * **Output order is input order.** Java iterates a `HashMap`, so its output
#     order is unspecified; keeping first-occurrence order makes results
#     reproducible at no cost.
#
# Nothing here is exported; the public surface is `OverlayNG` (api.jl).

# ## The point map (port of `buildPointMap`)

# JTS keys a `HashMap<Coordinate, Point>` on the coordinate and stores the
# first `Point` seen there. Here `keys` preserves first-occurrence order and
# `coords` maps the kernel point `P` to the output coordinate `T`.
struct _PointMap{P, T}
    keys::Vector{P}
    coords::Dict{P, T}
end

_PointMap(m::Manifold, ::Type{T}) where {T} =
    (P = _kernel_point_type(m); _PointMap{P, T}(P[], Dict{P, T}()))

Base.haskey(pm::_PointMap, k) = haskey(pm.coords, k)
Base.length(pm::_PointMap) = length(pm.keys)

function _point_map(m::Manifold, ::Type{T}, geom) where {T}
    pm = _PointMap(m, T)
    _point_map_add_all!(m, pm, GI.trait(geom), geom)
    return pm
end

# Port of the `GeometryComponentFilter` in `buildPointMap`: visit every
# non-empty `Point` component. Only point geometries reach here (the driver
# dispatches on dimension), so anything else is a bug.
_point_map_add_all!(m, pm, ::GI.PointTrait, geom) =
    (_ov_isempty(geom) || _point_map_add!(m, pm, geom); nothing)

function _point_map_add_all!(m, pm, ::GI.MultiPointTrait, geom)
    for p in GI.getgeom(geom)
        _ov_isempty(p) && continue
        _point_map_add!(m, pm, p)
    end
    return nothing
end

_point_map_add_all!(m, pm, trait, geom) = throw(ArgumentError(
    "OverlayNG point overlay: expected a point geometry, got $(typeof(trait))"))

# Only the first occurrence of a coordinate is kept — this is the merging
# semantics of overlay.
function _point_map_add!(m::Manifold, pm::_PointMap{P, T}, p) where {P, T}
    k = _to_kernel_point(m, p)
    haskey(pm.coords, k) && return nothing
    push!(pm.keys, k)
    pm.coords[k] = _input_output_point(T, p, k)
    return nothing
end

# The output coordinate of an INPUT point — the point paths' counterpart of
# `_emit_node_coord` on a vertex node, and the same pass-through: when the
# output type is the manifold's kernel type the kernel point already IS the
# answer, bit for bit. The lon/lat row takes the input's own coordinates rather
# than converting the kernel point back, which is the same thing one rounding
# earlier.
@inline _input_output_point(::Type{Tuple{Float64, Float64}}, p, k) =
    (Float64(GI.x(p)), Float64(GI.y(p)))
@inline _input_output_point(::Type{<:UnitSphericalPoint}, p, k) = k

# ## The overlay (port of `getResult`)

# Point × point, so the result is dimension 0 under every op. A target above that
# never reaches here — `_overlay_ng` answers it before dispatching — so there is
# nothing to elide, only the result shape to honour.
function _overlay_points(m::Manifold, ::Type{T}, op::_OverlayOpCode, a, b,
        target = nothing) where {T}
    map_a = _point_map(m, T, a)
    map_b = _point_map(m, T, b)

    result = T[]
    if op == OVERLAY_INTERSECTION
        _points_intersection!(result, map_a, map_b)
    elseif op == OVERLAY_UNION
        _points_union!(result, map_a, map_b)
    elseif op == OVERLAY_DIFFERENCE
        _points_difference!(result, map_a, map_b)
    else # OVERLAY_SYMDIFFERENCE
        _points_difference!(result, map_a, map_b)
        _points_difference!(result, map_b, map_a)
    end
    return _dimensional_result(target, 0, _NO_COMPONENTS, _NO_COMPONENTS, result)
end

function _points_intersection!(result, map_a::_PointMap, map_b::_PointMap)
    for k in map_a.keys
        haskey(map_b, k) && push!(result, map_a.coords[k])
    end
    return result
end

function _points_difference!(result, map_a::_PointMap, map_b::_PointMap)
    for k in map_a.keys
        haskey(map_b, k) || push!(result, map_a.coords[k])
    end
    return result
end

function _points_union!(result, map_a::_PointMap, map_b::_PointMap)
    #-- copy all A points, then the B points not in A
    for k in map_a.keys
        push!(result, map_a.coords[k])
    end
    for k in map_b.keys
        haskey(map_a, k) || push!(result, map_b.coords[k])
    end
    return result
end
