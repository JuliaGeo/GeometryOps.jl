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
# `coords` maps the kernel point to the output coordinate.
struct _PointMap{P}
    keys::Vector{P}
    coords::Dict{P, Tuple{Float64, Float64}}
end

_PointMap(m::Manifold) = (P = _kernel_point_type(m);
                          _PointMap{P}(P[], Dict{P, Tuple{Float64, Float64}}()))

Base.haskey(pm::_PointMap, k) = haskey(pm.coords, k)
Base.length(pm::_PointMap) = length(pm.keys)

function _point_map(m::Manifold, geom)
    pm = _PointMap(m)
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
function _point_map_add!(m::Manifold, pm::_PointMap, p)
    k = _to_kernel_point(m, p)
    haskey(pm.coords, k) && return nothing
    push!(pm.keys, k)
    pm.coords[k] = (Float64(GI.x(p)), Float64(GI.y(p)))
    return nothing
end

# ## The overlay (port of `getResult`)

function _overlay_points(m::Manifold, op::_OverlayOpCode, a, b)
    map_a = _point_map(m, a)
    map_b = _point_map(m, b)

    result = Tuple{Float64, Float64}[]
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
    return _create_point_result(result)
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

# The most specific dimension-0 geometry over the result coordinates
# (`GeometryFactory.buildGeometry` restricted to points).
function _create_point_result(coords::Vector{Tuple{Float64, Float64}})
    isempty(coords) && return _empty_geom(0)
    length(coords) == 1 && return GI.Point(coords[1])
    return GI.MultiPoint(coords)
end
