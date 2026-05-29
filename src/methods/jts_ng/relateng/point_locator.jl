# # RelateNG point-location substrate

const RELATE_EXTERIOR = DimensionLocation(dim_false, loc_exterior)
const _RELATE_ADJACENT_PROBE_DIRECTIONS = (
    (1.0, 0.4142135623730951),
    (0.4142135623730951, 1.0),
    (-0.4142135623730951, 1.0),
    (-1.0, 0.4142135623730951),
    (-1.0, -0.4142135623730951),
    (-0.4142135623730951, -1.0),
    (0.4142135623730951, -1.0),
    (1.0, -0.4142135623730951),
)

"""
    RelatePointLocator(geom; boundary_node_rule, prepared)

Point locator returning both topological location and containing dimension.
"""
struct RelatePointLocator{G,R <: BoundaryNodeRule}
    geom::G
    boundary_node_rule::R
    prepared::Bool
end

function RelatePointLocator(
    geom;
    boundary_node_rule::BoundaryNodeRule = Mod2BoundaryNodeRule(),
    prepared::Bool = false,
)
    return RelatePointLocator(geom, boundary_node_rule, prepared)
end

RelatePointLocator(alg::RelateNG, geom) =
    RelatePointLocator(geom; boundary_node_rule = alg.boundary_node_rule, prepared = alg.prepared)

function _relate_dimension_location(dimension::TopologicalDimension, location::TopologicalLocation)
    location == loc_exterior && return RELATE_EXTERIOR
    return DimensionLocation(dimension, location)
end

_relate_area_location(location::TopologicalLocation) = _relate_dimension_location(dim_area, location)
_relate_line_location(location::TopologicalLocation) = _relate_dimension_location(dim_line, location)
_relate_point_location(location::TopologicalLocation) = _relate_dimension_location(dim_point, location)

"""
    relate_locate(locator, point; exact = True())

Return only the topological location of `point` relative to the locator geometry.
"""
relate_locate(locator::RelatePointLocator, point; kwargs...) =
    relate_locate_with_dim(locator, point; kwargs...).location

"""
    relate_locate_with_dim(locator, point; exact = True())

Return the RelateNG dimension/location for `point`.
"""
function relate_locate_with_dim(
    locator::RelatePointLocator,
    point;
    is_node::Bool = false,
    parent_polygonal = nothing,
    exact = True(),
)
    _relate_is_empty(locator.geom) && return RELATE_EXTERIOR

    if is_node && _relate_is_polygonal(GI.trait(locator.geom), locator.geom)
        return DimensionLocation(dim_area, loc_boundary)
    end
    return _relate_compute_dim_location(locator, point, is_node, parent_polygonal; exact)
end

"""
    relate_locate_node_with_dim(locator, point, parent_polygonal; exact = True())

Locate a point already known to be a node of the target geometry.
"""
relate_locate_node_with_dim(locator::RelatePointLocator, point, parent_polygonal = nothing; kwargs...) =
    relate_locate_with_dim(locator, point; is_node = true, parent_polygonal, kwargs...)

"""
    relate_locate_line_end_with_dim(locator, point; exact = True())

Locate a line endpoint, applying area precedence and the boundary-node rule.
"""
function relate_locate_line_end_with_dim(locator::RelatePointLocator, point; exact = True())
    area_loc = _relate_locate_on_polygonals(locator.geom, point, false, nothing; exact)
    area_loc != loc_exterior && return _relate_area_location(area_loc)

    return _relate_is_line_boundary(locator, point) ?
        DimensionLocation(dim_line, loc_boundary) :
        DimensionLocation(dim_line, loc_interior)
end

function _relate_compute_dim_location(locator::RelatePointLocator, point, is_node::Bool, parent_polygonal; exact)
    area_loc = _relate_locate_on_polygonals(locator.geom, point, is_node, parent_polygonal; exact)
    area_loc != loc_exterior && return _relate_area_location(area_loc)

    line_loc = _relate_locate_on_lines(locator, locator.geom, point, is_node)
    line_loc != loc_exterior && return _relate_line_location(line_loc)

    point_loc = _relate_locate_on_points(locator.geom, point)
    point_loc != loc_exterior && return _relate_point_location(point_loc)

    return RELATE_EXTERIOR
end

_relate_is_polygonal(::GI.PolygonTrait, geom) = true
_relate_is_polygonal(::GI.MultiPolygonTrait, geom) = true
_relate_is_polygonal(trait, geom) = false

function _relate_locate_on_polygonals(geom, point, is_node::Bool, parent_polygonal; exact)
    return _relate_locate_on_polygonals(GI.trait(geom), geom, point, is_node, parent_polygonal; exact)
end

_relate_locate_on_polygonals(::GI.PointTrait, geom, point, is_node, parent_polygonal; exact) = loc_exterior
_relate_locate_on_polygonals(::GI.MultiPointTrait, geom, point, is_node, parent_polygonal; exact) = loc_exterior
_relate_locate_on_polygonals(::GI.AbstractCurveTrait, geom, point, is_node, parent_polygonal; exact) = loc_exterior

function _relate_locate_on_polygonals(::GI.FeatureCollectionTrait, fc, point, is_node, parent_polygonal; exact)
    polygonal_locations = TopologicalLocation[]
    for feature in GI.getfeature(fc)
        loc = _relate_locate_on_polygonals(feature, point, is_node, parent_polygonal; exact)
        loc == loc_interior && return loc_interior
        loc == loc_boundary && push!(polygonal_locations, loc)
    end
    return _relate_resolve_polygonal_boundary_count(length(polygonal_locations), fc, point; exact)
end

function _relate_locate_on_polygonals(::GI.FeatureTrait, feature, point, is_node, parent_polygonal; exact)
    geom = GI.geometry(feature)
    return _relate_locate_on_polygonals(GI.trait(geom), geom, point, is_node, parent_polygonal; exact)
end

function _relate_locate_on_polygonals(::GI.PolygonTrait, polygon, point, is_node, parent_polygonal; exact)
    return _relate_locate_on_polygonal(polygon, point, is_node, parent_polygonal; exact)
end

function _relate_locate_on_polygonals(::GI.MultiPolygonTrait, multipolygon, point, is_node, parent_polygonal; exact)
    return _relate_locate_on_polygonal_collection(
        multipolygon,
        GI.getpolygon(multipolygon),
        point,
        is_node,
        isnothing(parent_polygonal) ? multipolygon : parent_polygonal;
        exact,
    )
end

function _relate_locate_on_polygonals(::GI.AbstractGeometryTrait, geom, point, is_node, parent_polygonal; exact)
    GI.isempty(geom) && return loc_exterior

    polygonal_locations = TopologicalLocation[]
    for child in GI.getgeom(geom)
        loc = _relate_locate_on_polygonals(child, point, is_node, parent_polygonal; exact)
        loc == loc_interior && return loc_interior
        loc == loc_boundary && push!(polygonal_locations, loc)
    end
    return _relate_resolve_polygonal_boundary_count(length(polygonal_locations), geom, point; exact)
end

function _relate_locate_on_polygonals(::Nothing, iterable, point, is_node, parent_polygonal; exact)
    polygonal_locations = TopologicalLocation[]
    for child in iterable
        loc = _relate_locate_on_polygonals(child, point, is_node, parent_polygonal; exact)
        loc == loc_interior && return loc_interior
        loc == loc_boundary && push!(polygonal_locations, loc)
    end
    return length(polygonal_locations) > 0 ? loc_boundary : loc_exterior
end

function _relate_locate_on_polygonal_collection(polygonal_geom, polygons, point, is_node, parent_polygonal; exact)
    boundary_count = 0
    for polygon in polygons
        loc = _relate_locate_on_polygonal(polygon, point, is_node, parent_polygonal; exact)
        loc == loc_interior && return loc_interior
        boundary_count += loc == loc_boundary
    end
    return _relate_resolve_polygonal_boundary_count(boundary_count, polygonal_geom, point; exact)
end

function _relate_resolve_polygonal_boundary_count(boundary_count::Integer, polygonal_geom, point; exact)
    boundary_count == 0 && return loc_exterior
    boundary_count == 1 && return loc_boundary
    return _relate_locate_adjacent_area_boundary(polygonal_geom, point; exact)
end

function _relate_locate_on_polygonal(polygon, point, is_node::Bool, parent_polygonal; exact)
    if is_node && _relate_is_parent_polygonal(parent_polygonal, polygon)
        return loc_boundary
    end
    return _relate_point_polygon_location(point, polygon; exact)
end

function _relate_point_polygon_location(point, polygon; exact)
    exterior_location = _point_filled_curve_orientation(point, GI.getexterior(polygon); exact)
    exterior_location == point_out && return loc_exterior
    exterior_location == point_on && return loc_boundary

    for hole in GI.gethole(polygon)
        hole_location = _point_filled_curve_orientation(point, hole; exact)
        hole_location == point_in && return loc_exterior
        hole_location == point_on && return loc_boundary
    end
    return loc_interior
end

function _relate_locate_adjacent_area_boundary(polygonal_geom, point; exact)
    probe_distance = _relate_probe_distance(polygonal_geom)
    x, y = _tuple_point(point)

    for (dx, dy) in _RELATE_ADJACENT_PROBE_DIRECTIONS
        probe = (x + probe_distance * dx, y + probe_distance * dy)
        loc = _relate_probe_polygonals(GI.trait(polygonal_geom), polygonal_geom, probe; exact)
        loc == loc_exterior && return loc_boundary
    end
    return loc_interior
end

function _relate_probe_distance(geom)
    extent = GI.extent(geom)
    xspan = extent.X[2] - extent.X[1]
    yspan = extent.Y[2] - extent.Y[1]
    span = max(abs(xspan), abs(yspan), 1.0)
    return max(span * sqrt(eps(Float64)), eps(Float64))
end

_relate_probe_polygonals(::GI.PointTrait, geom, point; exact) = loc_exterior
_relate_probe_polygonals(::GI.MultiPointTrait, geom, point; exact) = loc_exterior
_relate_probe_polygonals(::GI.AbstractCurveTrait, geom, point; exact) = loc_exterior

function _relate_probe_polygonals(::GI.FeatureCollectionTrait, fc, point; exact)
    return _relate_probe_polygonals(nothing, GI.getfeature(fc), point; exact)
end

function _relate_probe_polygonals(::GI.FeatureTrait, feature, point; exact)
    geom = GI.geometry(feature)
    return _relate_probe_polygonals(GI.trait(geom), geom, point; exact)
end

function _relate_probe_polygonals(::GI.PolygonTrait, polygon, point; exact)
    return _relate_point_polygon_location(point, polygon; exact)
end

function _relate_probe_polygonals(::GI.MultiPolygonTrait, multipolygon, point; exact)
    boundary_count = 0
    for polygon in GI.getpolygon(multipolygon)
        loc = _relate_point_polygon_location(point, polygon; exact)
        loc == loc_interior && return loc_interior
        boundary_count += loc == loc_boundary
    end
    boundary_count > 0 && return loc_boundary
    return loc_exterior
end

function _relate_probe_polygonals(::GI.AbstractGeometryTrait, geom, point; exact)
    boundary_count = 0
    for child in GI.getgeom(geom)
        loc = _relate_probe_polygonals(GI.trait(child), child, point; exact)
        loc == loc_interior && return loc_interior
        boundary_count += loc == loc_boundary
    end
    boundary_count > 0 && return loc_boundary
    return loc_exterior
end

function _relate_probe_polygonals(::Nothing, iterable, point; exact)
    boundary_count = 0
    for child in iterable
        loc = _relate_probe_polygonals(GI.trait(child), child, point; exact)
        loc == loc_interior && return loc_interior
        boundary_count += loc == loc_boundary
    end
    boundary_count > 0 && return loc_boundary
    return loc_exterior
end

function _relate_is_parent_polygonal(parent_polygonal, polygon)
    isnothing(parent_polygonal) && return false
    parent_polygonal === polygon && return true

    trait = GI.trait(parent_polygonal)
    if trait isa GI.FeatureTrait
        return _relate_is_parent_polygonal(GI.geometry(parent_polygonal), polygon)
    elseif trait isa GI.FeatureCollectionTrait
        return any(feature -> _relate_is_parent_polygonal(feature, polygon), GI.getfeature(parent_polygonal))
    end
    if trait isa Union{GI.MultiPolygonTrait,GI.GeometryCollectionTrait}
        return any(child -> child === polygon, GI.getgeom(parent_polygonal))
    end
    return false
end

function _relate_locate_on_lines(locator::RelatePointLocator, geom, point, is_node::Bool)
    _relate_is_line_boundary(locator, point) && return loc_boundary
    is_node && _relate_has_lines(GI.trait(geom), geom) && return loc_interior
    return _relate_point_on_any_line(GI.trait(geom), geom, point) ? loc_interior : loc_exterior
end

function _relate_is_line_boundary(locator::RelatePointLocator, point)
    boundary_count = _relate_line_boundary_count(GI.trait(locator.geom), locator.geom, point)
    return is_in_boundary(locator.boundary_node_rule, boundary_count)
end

_relate_has_lines(::GI.FeatureTrait, feature) =
    _relate_has_lines(GI.trait(GI.geometry(feature)), GI.geometry(feature))

_relate_has_lines(::GI.FeatureCollectionTrait, fc) =
    _relate_has_lines(nothing, GI.getfeature(fc))

function _relate_has_lines(::Nothing, iterable)
    return any(child -> _relate_has_lines(GI.trait(child), child), iterable)
end

_relate_has_lines(::GI.PointTrait, geom) = false
_relate_has_lines(::GI.MultiPointTrait, geom) = false
_relate_has_lines(::GI.PolygonTrait, geom) = false
_relate_has_lines(::GI.MultiPolygonTrait, geom) = false
_relate_has_lines(::GI.AbstractCurveTrait, geom) = !GI.isempty(geom)

function _relate_has_lines(::GI.AbstractGeometryTrait, geom)
    GI.isempty(geom) && return false
    return any(child -> _relate_has_lines(GI.trait(child), child), GI.getgeom(geom))
end

_relate_line_boundary_count(::GI.PointTrait, geom, point) = 0
_relate_line_boundary_count(::GI.MultiPointTrait, geom, point) = 0
_relate_line_boundary_count(::GI.PolygonTrait, geom, point) = 0
_relate_line_boundary_count(::GI.MultiPolygonTrait, geom, point) = 0

function _relate_line_boundary_count(::GI.FeatureTrait, feature, point)
    geom = GI.geometry(feature)
    return _relate_line_boundary_count(GI.trait(geom), geom, point)
end

function _relate_line_boundary_count(::GI.FeatureCollectionTrait, fc, point)
    return _relate_line_boundary_count(nothing, GI.getfeature(fc), point)
end

function _relate_line_boundary_count(::Nothing, iterable, point)
    count = 0
    for child in iterable
        count += _relate_line_boundary_count(GI.trait(child), child, point)
    end
    return count
end

function _relate_line_boundary_count(::GI.AbstractCurveTrait, curve, point)
    GI.isempty(curve) && return 0

    count = 0
    point = _tuple_point(point)
    count += _tuple_point(GI.getpoint(curve, 1)) == point
    count += _tuple_point(GI.getpoint(curve, GI.npoint(curve))) == point
    return count
end

function _relate_line_boundary_count(::GI.AbstractGeometryTrait, geom, point)
    GI.isempty(geom) && return 0
    count = 0
    for child in GI.getgeom(geom)
        count += _relate_line_boundary_count(GI.trait(child), child, point)
    end
    return count
end

_relate_point_on_any_line(::GI.PointTrait, geom, point) = false
_relate_point_on_any_line(::GI.MultiPointTrait, geom, point) = false
_relate_point_on_any_line(::GI.PolygonTrait, geom, point) = false
_relate_point_on_any_line(::GI.MultiPolygonTrait, geom, point) = false

function _relate_point_on_any_line(::GI.FeatureTrait, feature, point)
    geom = GI.geometry(feature)
    return _relate_point_on_any_line(GI.trait(geom), geom, point)
end

function _relate_point_on_any_line(::GI.FeatureCollectionTrait, fc, point)
    return _relate_point_on_any_line(nothing, GI.getfeature(fc), point)
end

function _relate_point_on_any_line(::Nothing, iterable, point)
    for child in iterable
        _relate_point_on_any_line(GI.trait(child), child, point) && return true
    end
    return false
end

function _relate_point_on_any_line(::GI.AbstractCurveTrait, curve, point)
    GI.isempty(curve) && return false
    _point_in_extent(point, GI.extent(curve)) || return false
    return _relate_point_on_curve(point, curve)
end

function _relate_point_on_any_line(::GI.AbstractGeometryTrait, geom, point)
    GI.isempty(geom) && return false
    for child in GI.getgeom(geom)
        _relate_point_on_any_line(GI.trait(child), child, point) && return true
    end
    return false
end

function _relate_point_on_curve(point, curve)
    npoints = GI.npoint(curve)
    npoints < 2 && return false

    previous = GI.getpoint(curve, 1)
    for i in 2:npoints
        current = GI.getpoint(curve, i)
        if _point_segment_orientation(point, previous, current; in = true, on = true, out = false)
            return true
        end
        previous = current
    end
    return false
end

function _relate_locate_on_points(geom, point)
    points = extract_ng_points(geom)
    point = _tuple_point(point)
    return any(extracted -> extracted.point == point, points) ? loc_interior : loc_exterior
end
