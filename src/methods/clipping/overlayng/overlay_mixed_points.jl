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

# Empty component list for the dimension the non-point input does not have.
# `_create_result_geometry` only ever tests it with `isempty`.
const _NO_COMPONENTS = Any[]

# Port of `getResult` (with the constructor's dimensional naming inlined).
function _overlay_mixed_points(m::Manifold, op::_OverlayOpCode, a, b, dim_a, dim_b; exact)
    #-- name the dimensional geometries (JTS constructor)
    if dim_a == 0
        geom_point, geom_non_point, non_point_dim, is_point_rhs = a, b, dim_b, false
    else
        geom_point, geom_non_point, non_point_dim, is_point_rhs = b, a, dim_a, true
    end

    locator = _mixed_point_locator(m, geom_non_point, non_point_dim; exact)
    coords = _point_map(m, geom_point)

    if op == OVERLAY_INTERSECTION
        return _create_point_result(_find_points(locator, coords, true))
    elseif op == OVERLAY_UNION || op == OVERLAY_SYMDIFFERENCE
        #-- UNION and SYMDIFFERENCE have the same output: the non-point input
        #-- plus the points lying outside it
        points = _find_points(locator, coords, false)
        comps = _nonpoint_components(geom_non_point)
        polys = non_point_dim == 2 ? comps : _NO_COMPONENTS
        lines = non_point_dim == 1 ? comps : _NO_COMPONENTS
        isempty(polys) && isempty(lines) && isempty(points) &&
            return _empty_geom(_result_dimension(op, dim_a, dim_b))
        return _create_result_geometry(polys, lines, points)
    end
    #-- OVERLAY_DIFFERENCE: removing points from a higher-dimension geometry
    #-- changes nothing, so A survives untouched when the points are on the RHS
    is_point_rhs && return _copy_non_point(geom_non_point, non_point_dim)
    return _create_point_result(_find_points(locator, coords, false))
end

# Port of `createLocator`.
_mixed_point_locator(m::Manifold, geom, dim::Integer; exact) =
    dim == 2 ? IndexedPointInAreaLocator(m, geom; exact) :
               RelatePointLocator(m, geom; exact)

# Port of `findPoints` + `hasLocation` + `createPoints`. The keys of
# `_PointMap` are kernel points, which is what both locators consume.
function _find_points(locator, coords::_PointMap, is_covered::Bool)
    out = Tuple{Float64, Float64}[]
    for k in coords.keys
        is_exterior = locate(locator, k) == LOC_EXTERIOR
        (is_covered ? !is_exterior : is_exterior) && push!(out, coords.coords[k])
    end
    return out
end

# Port of `extractPolygons`/`extractLines`: the non-empty atomic components of
# the non-point input, as engine-native (tuple-coordinate) geometries.
function _nonpoint_components(geom)
    _ov_isempty(geom) && return _NO_COMPONENTS
    g = tuples(geom)
    t = GI.trait(g)
    if t isa GI.MultiPolygonTrait || t isa GI.MultiLineStringTrait
        return [c for c in GI.getgeom(g) if !_ov_isempty(c)]
    end
    return [g]
end

# Port of `copyNonPoint`. Java copies because its precision-reduction step may
# have aliased the input; here `tuples` is the copy.
function _copy_non_point(geom, dim::Integer)
    _ov_isempty(geom) && return _empty_geom(dim)
    return tuples(geom)
end
