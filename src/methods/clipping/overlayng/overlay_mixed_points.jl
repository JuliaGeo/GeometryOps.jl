# # OverlayMixedPoints — point × (line | area) overlay (port of JTS `OverlayMixedPoints`)
#
# Phase 3 of the OverlayNG port (design doc §4). A port of
# `operation/overlayng/OverlayMixedPoints.java`: the overlay of one point input
# against one higher-dimension input. Points cannot node anything, so the
# result is decided entirely by locating each input point against the non-point
# input and, for the ops that keep it, copying the non-point input through.
# This is also what makes `intersection(OverlayNG(), points, polygon)` an
# efficient "which points are inside" query, as the Javadoc advertises.
#
# Location goes through the kernel locators — `IndexedPointInAreaLocator` for an
# area, `RelatePointLocator` for a line — so both manifolds are exact and
# spherical location follows enclosed-region semantics. `IndexedPointOnLineLocator`
# is not ported: it is a JTS shim whose `locate` delegates straight to
# `PointLocator` (its own `// TODO: optimize this with a segment index`), and
# `RelatePointLocator` is this repo's port of that behaviour.
#
# Deviations from the Java, each deliberate:
#
#   * **`prepareNonPoint` is dropped.** Java re-nodes and rounds the non-point
#     input through `OverlayNG.union(geom, pm)` so the output is valid at the
#     precision model. There is no precision model here (§0), and inputs are
#     contractually valid (§2.2), so a self-union is the identity — the
#     non-point input is copied through as-is. The knock-on JTS caveat that the
#     points are compared against the *un*-rounded geometry when `resultDim ==
#     0` disappears with it: here there is only one geometry to compare against.
#   * **No precision model** anywhere: `OverlayUtil.round` is a no-op and is
#     dropped from `extractCoordinates`.
#   * Duplicate input points are removed in first-occurrence order (Java uses a
#     `HashSet`, whose order is unspecified), reusing `_point_map`.
#
# Nothing here is exported; the public surface is `OverlayNG` (api.jl).

# Port of `getResult` (with the constructor's dimensional naming inlined).
#
# Both halves of the work are targetable and both are worth eliding: locating
# every input point is the whole cost of the point half (this is the "which of my
# points are inside" query), and `_nonpoint_components` deep-copies the entire
# non-point input through `tuples`. An areal target on a point ∪ polygon skips
# the first; a point target on the same skips the second.
function _overlay_mixed_points(m::Manifold, ::Type{T}, op::_OverlayOpCode, a, b, dim_a, dim_b,
        target = nothing; exact) where {T}
    #-- name the dimensional geometries (JTS constructor)
    if dim_a == 0
        geom_point, geom_non_point, non_point_dim, is_point_rhs = a, b, dim_b, false
    else
        geom_point, geom_non_point, non_point_dim, is_point_rhs = b, a, dim_a, true
    end
    result_dim = _result_dimension(op, dim_a, dim_b)

    if op == OVERLAY_INTERSECTION
        points = _mixed_points(m, T, geom_point, geom_non_point, non_point_dim, true, target; exact)
        return _dimensional_result(target, result_dim, _NO_COMPONENTS, _NO_COMPONENTS, points)
    elseif op == OVERLAY_UNION || op == OVERLAY_SYMDIFFERENCE
        #-- UNION and SYMDIFFERENCE have the same output: the non-point input
        #-- plus the points lying outside it
        points = _mixed_points(m, T, geom_point, geom_non_point, non_point_dim, false, target; exact)
        polys, lines = _mixed_components(m, T, geom_non_point, non_point_dim, target)
        return _dimensional_result(target, result_dim, polys, lines, points)
    end
    #-- OVERLAY_DIFFERENCE: removing points from a higher-dimension geometry
    #-- changes nothing, so A survives untouched when the points are on the RHS
    if is_point_rhs
        polys, lines = _mixed_components(m, T, geom_non_point, non_point_dim, target)
        return _dimensional_result(target, result_dim, polys, lines, T[])
    end
    points = _mixed_points(m, T, geom_point, geom_non_point, non_point_dim, false, target; exact)
    return _dimensional_result(target, result_dim, _NO_COMPONENTS, _NO_COMPONENTS, points)
end

# The input points the op keeps: those covered by the non-point input, or those
# outside it. Neither the locator nor the point map is built when the target
# excludes points.
function _mixed_points(m::Manifold, ::Type{T}, geom_point, geom_non_point,
        non_point_dim::Integer, is_covered::Bool, target; exact) where {T}
    _target_needs_dim(target, 0) || return T[]
    locator = _mixed_point_locator(m, geom_non_point, non_point_dim; exact)
    return _find_points(locator, _point_map(m, T, geom_point), is_covered)
end

# The non-point input's components as `(polys, lines)`, one of which is always
# empty — and both are when the target excludes the non-point dimension.
function _mixed_components(m::Manifold, ::Type{T}, geom_non_point, non_point_dim::Integer,
        target) where {T}
    comps = _target_needs_dim(target, non_point_dim) ?
        _nonpoint_components(m, T, geom_non_point) : _NO_COMPONENTS
    return non_point_dim == 2 ? (comps, _NO_COMPONENTS) : (_NO_COMPONENTS, comps)
end

# Port of `createLocator`.
_mixed_point_locator(m::Manifold, geom, dim::Integer; exact) =
    dim == 2 ? IndexedPointInAreaLocator(m, geom; exact) :
               RelatePointLocator(m, geom; exact)

# Port of `findPoints` + `hasLocation` + `createPoints`. The keys of
# `_PointMap` are kernel points, which is what both locators consume.
function _find_points(locator, coords::_PointMap{P, T}, is_covered::Bool) where {P, T}
    out = T[]
    for k in coords.keys
        is_exterior = locate(locator, k) == LOC_EXTERIOR
        (is_covered ? !is_exterior : is_exterior) && push!(out, coords.coords[k])
    end
    return out
end

# Port of `extractPolygons`/`extractLines`: the non-empty atomic components of
# the non-point input, as engine-native (tuple-coordinate) geometries.
#
# Empty components are dropped BEFORE the `tuples` conversion, not after: `tuples`
# rebuilds each ring/linestring from `first(points)` and so cannot represent an
# empty component (`MULTILINESTRING ((10 10, 20 20), EMPTY)`, TestNGOverlayEmpty
# case 9, is exactly that input).
function _nonpoint_components(m::Manifold, ::Type{T}, geom) where {T}
    _ov_isempty(geom) && return _NO_COMPONENTS
    t = GI.trait(geom)
    if t isa GI.MultiPolygonTrait || t isa GI.MultiLineStringTrait
        return [_output_component(m, T, c) for c in GI.getgeom(geom) if !_ov_isempty(c)]
    end
    return [_output_component(m, T, geom)]
end

# One input component rebuilt in the engine's output point type — the copy that
# stands in for JTS's `copyNonPoint`, and the only place in the engine where a
# result geometry is made from input coordinates rather than from arrangement
# nodes. It has to agree with `_Result*` exactly, because `_typed` and the
# `GeometryCollection` element Union both reject anything else.
#
# The lon/lat row goes through `tuples`, unchanged. The xyz row converts each
# vertex with the manifold's own ingest (`_to_kernel_point`), which is precisely
# what the noded path would have done to the same vertex, so a point copied
# through here and the same point emitted from an arrangement are identical.
_output_component(::Manifold, ::Type{Tuple{Float64, Float64}}, geom) = tuples(geom)
_output_component(m::Manifold, ::Type{<:UnitSphericalPoint}, geom) =
    apply(GI.PointTrait(), geom) do p
        _to_kernel_point(m, p)
    end

#-- `copyNonPoint` is not a separate function here: Java copies because its
#-- precision-reduction step may have aliased the input, and `_nonpoint_components`
#-- (whose `tuples` is that copy) already does the job, with `_dimensional_result`
#-- picking the atomic-vs-multi form exactly as `copyNonPoint` did.
