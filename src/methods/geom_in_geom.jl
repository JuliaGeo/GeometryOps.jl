export point_in_geom, point_in_polygon

"""
    point_in_geom(point, geom)::(Bool, Bool)

Returns if a point is within a given geometry. Returns a boolean tuple with two
elements. The first element of the tuple is true if the point is within the
geometry (excluding edges and vertices) and false othereise. The second element
of the tuple is true if the point is on the geometry, the point is on an edge or
is a vertex, and false otherwise.
"""
point_in_geom(point, geom) = point_in_geom(
    GI.trait(point), point,
    GI.trait(geom), geom,
)

line_in_geom(line, geom) = line_in_geom(
    GI.trait(line), line,
    GI.trait(geom), geom,
)

"""
    point_in_geom(
        ::GI.PointTrait, point,
        ::GI.LineStringTrait, linestring,
    )::(Bool, Bool)

Returns a boolean tuple with two elements. The first element is if the point is
within the linestring. Note this is only possible if the linestring is closed.
If the linestring isn't closed (repeated last point), this will throw a warning.
The second element is if the point is on the linestring. 
"""
function point_in_geom(::GI.PointTrait, point, ::GI.LineStringTrait, linestring)
    results = if equals(
        GI.getpoint(linestring, 1),
        GI.getpoint(linestring, GI.npoint(linestring)),
    )
        _point_in_closed_curve(point, linestring)
    else
        @warn "Linestring isn't closed. Point cannot be 'in' linestring."
        (false, false)  # TODO: see if point is actually on linestring!
    end
    return results
end

"""
    point_in_geom(
        ::GI.PointTrait, point,
        ::GI.LinearRingTrait, linearring,
    )::(Bool, Bool)

Returns a boolean tuple with two elements. The first element is if the point is
within the linear ring. The second element is if the point is on the linestring. 
"""
function point_in_geom(::GI.PointTrait, point, ::GI.LinearRingTrait, linearring)
    _point_in_extent(point, GI.extent(linearring)) || return (false, false)
    return _point_in_closed_curve(point, linearring)
end

"""
    point_in_geom(
        ::GI.PointTrait, point,
        ::GI.PolygonTrait, poly,
    )::(Bool, Bool)

Returns a boolean tuple with two elements. The first element is if the point is
within the polygon. This means that it also isn't within any holes. The second
element is if the point is on the polygon, including edges and vertices of the
exterior ring and any holes. 
"""
function point_in_geom(::GI.PointTrait, point, ::GI.PolygonTrait, poly)
    # Cheaply check that the point is inside the polygon extent
    _point_in_extent(point, GI.extent(poly)) || return (false, false)
    # Check if point is inside or on the exterior ring
    in_ext, on_ext = _point_in_closed_curve(point, GI.getexterior(poly))
    on_ext && return (in_ext, on_ext)  # point in on external boundary
    # Check if the point is in any of the holes
    for ring in GI.gethole(poly)
        in_hole, on_hole = _point_in_closed_curve(point, ring)
        in_hole && return (false, false)  # point is in a hole -> not in polygon
        on_hole && return (false, on_hole)  # point is on an edge
    end
    return (in_ext, on_ext)  # point is inside of polygon
end

"""
    point_in_polygon(point, polygon)::(Bool, Bool)

Determines if point is within a polygon, returning a tuple where the first
element is if the point is within the polygon edges, and the second is if the
point is on an edge or vertex.
"""
point_in_polygon(point, polygon) = point_in_polygon(
    GI.trait(point), point,
    GI.trait(polygon), polygon,
)

"""
    point_in_polygon(
        ::GI.PointTrait, point,
        ::GI.PolygonTrait, poly,
    )

Returns a boolean tuple with two elements. The first element is if the point is
within the polygon and the second element is if the point is on the polygon.
    
Note that this is the same as point_in_geom dispatched on a polygon. 
"""
point_in_polygon(trait1::GI.PointTrait, point, trait2::GI.PolygonTrait, poly) =
    point_in_geom(trait1, point, trait2, poly)

line_in_polygon(line, polygon) = any(line_in_geom(line, polygon))
function line_in_geom(::GI.LineStringTrait, line, ::GI.PolygonTrait, poly)
    # Cheaply check that the line extent is inside the polygon extent
    Extents.within(GI.extent(line), GI.extent(poly), ) || return (false, false)
    point_in, point_on = point_in_polygon(GI.getpoint(line, 1), poly)
    (point_in || point_on) || return (false, false)
    vertex_on = point_on
    line_edges, poly_edges = map(sort! ∘ to_edges, (line, poly))
    for l_edge in line_edges
        for p_edge in poly_edges  # need to figure out closed vs not closed 
            _line_intersects(l_edge, p_edge) && return (false, false)
            if !vertex_on
                v1, _ = l_edge
                vertex_on = point_on_segment(v1, p_edge...)
            end
        end
    end
    return (!vertex_on, vertex_on)
end

function polygon_in_polygon(poly1, poly2)
    # edges1, edges2 = to_edges(poly1), to_edges(poly2)
    # extent1, extent2 = to_extent(edges1), to_extent(edges2)
    # Check the extents intersect
    Extents.intersects(GI.extent(poly1), GI.extent(poly2)) || return false

    # Check all points in poly1 are in poly2
    for point in GI.getpoint(poly1)
        point_in_polygon(point, poly2) || return false
    end

    # Check the line of poly1 does not intersect the line of poly2
    intersects(poly1, poly2) && return false

    # poly1 must be in poly2
    return true
 end

"""
    _point_in_closed_curve(point, curve)::(Bool, Bool)

Determine if point is within or on a closed curve. Point should be an object of
Point trait and curve should be a linearstring or ring, that is assumed to be
closed, regardless of repeated last point.

The return object is a boolean tuple (in_bounds, on_bounds). The in_bounds
object means that the point is within the curve, while on_bounds means the point
is on an edge.
"""
function _point_in_closed_curve(point, curve)
    # Determine number of points
    x, y = GI.x(point), GI.y(point)
    n = GI.npoint(curve)
    n -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, n)) ? 1 : 0
    #=
    Check if point is on an edge or if a ray, passing from (x, y) to infinity
    through line y = y intersects with the edge
    =#
    in_bounds = false
    on_bounds = false
    p_start = GI.getpoint(curve, n)
    for i in 1:n
        # Determine endpoints and edge lengths
        p_end = GI.getpoint(curve, i)
        xi, yi = GI.x(p_start), GI.y(p_start)
        xj, yj = GI.x(p_end), GI.y(p_end)
        Δx, Δy = xj - xi, yj - yi
        # Determine if point is on the edge
        on_bounds = point_on_segment(point, p_start, p_end)
        on_bounds && return (false, on_bounds)
        # Edge is vertical, just see if y is between edge endpoints
        if Δx == 0 && x < xi && (yi ≥ yj ? yj ≤ y ≤ yi : yi ≤ y ≤ yj)
            in_bounds = !in_bounds
        #=
        Edge is not vertical, find intersection point on y = y and see if it
        is between edge endpoints.
        =#
        elseif Δx != 0 && Δy != 0
            m = Δy / Δx
            b = yi - m * xi
            x_inter = (y - b) / m
            if (x_inter > x) && (xi ≥ xj ? xj < x_inter ≤ xi : xi ≤ x_inter < xj)
                in_bounds = !in_bounds
            end
        end
        p_start = p_end
    end
    return in_bounds, on_bounds
end

"""
    _point_in_extent(p, extent::Extents.Extent)::Bool

Returns true if the point is the bounding box of the extent and false otherwise. 
"""
function _point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y2) = extent.X, extent.Y
    return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
end
