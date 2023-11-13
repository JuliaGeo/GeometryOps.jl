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

ring_in_geom(ring, geom) = ring_in_geom(
    GI.trait(ring), ring,
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
function point_in_geom(::GI.PointTrait, point, ::GI.LineStringTrait, line)
    results = if equals(
        GI.getpoint(line, 1),
        GI.getpoint(line, GI.npoint(line)),
    )
        _point_in_closed_curve(point, line)
    else
        @warn "Linestring isn't closed. Point cannot be 'in' linestring."
        (false, false)
    end
    return results
end

function line_in_geom(::GI.LineStringTrait, line1, ::GI.LineStringTrait, line2)
    results = if equals(
        GI.getpoint(line2, 1),
        GI.getpoint(line2, GI.npoint(line2)),
    )
        Extents.within(
            GI.extent(line1),
            GI.extent(line2),
        ) || return (false, false)
        _line_in_closed_curve(line1, line2; close = false)
    else
        @warn "Linestring isn't closed. Point cannot be 'in' linestring."
        (false, false)
    end
    return results
end

function ring_in_geom(::GI.LinearRingTrait, ring, ::GI.LineStringTrait, line)
    results = if equals(
        GI.getpoint(line, 1),
        GI.getpoint(line, GI.npoint(line)),
    )
        Extents.within(
            GI.extent(ring),
            GI.extent(line),
        ) || return (false, false)
        _line_in_closed_curve(ring, line; close = true)
    else
        @warn "Linestring isn't closed. Point cannot be 'in' linestring."
        (false, false)
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
function point_in_geom(::GI.PointTrait, point, ::GI.LinearRingTrait, ring)
    _point_in_extent(point, GI.extent(ring)) || return (false, false)
    return _point_in_closed_curve(point, ring)
end

function line_in_geom(::GI.LineStringTrait, line, ::GI.LinearRingTrait, ring)
    Extents.within(GI.extent(line), GI.extent(ring)) || return (false, false)
    return _line_in_closed_curve(line, curve; close = false)
end

function ring_in_geom(::GI.LinearRingTrait, ring1, ::GI.LinearRingTrait, ring2)
    Extents.within(GI.extent(ring1), GI.extent(ring2)) || return (false, false)
    return _line_in_closed_curve(ring1, ring2; close = true)
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
    (in_ext || on_ext) || return (false, false)  # point isn't in external ring
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

line_in_polygon(line, polygon) = line_in_geom(line, polygon)

function line_in_geom(::GI.LineStringTrait, line, ::GI.PolygonTrait, poly)
    # Cheaply check that the line extent is inside the polygon extent
    Extents.within(GI.extent(line), GI.extent(poly)) || return (false, false)
    # Check if point is inside or on the exterior ring
    in_ext, on_ext = _line_in_closed_curve(
        line,
        GI.getexterior(poly);
        close = false,
    )
    (in_ext || on_ext) || return (false, false)  # line isn't in external ring
    # Check if the line is in any of the holes
    for ring in GI.gethole(poly)
        in_hole, on_hole = _line_in_closed_curve(point, ring; close = false)
        # point is in a hole -> not in polygon
        (in_hole || on_hole) && return (false, false)  # TODO: what if all points on the edge of hole?
    end
    return (in_ext, on_ext)  # point is inside of polygon
end

function ring_in_geom(::GI.LinearRingTrait, ring, ::GI.PolygonTrait, poly)
    # Cheaply check that the line extent is inside the polygon extent
    Extents.within(GI.extent(ring), GI.extent(poly)) || return (false, false)
    # Check if point is inside or on the exterior ring
    in_ext, on_ext = _line_in_closed_curve(
        ring,
        GI.getexterior(poly);
        close = false,
    )
    (in_ext || on_ext) || return (false, false)  # line isn't in external ring
    # Check if the line is in any of the holes
    for hole in GI.gethole(poly)
        in_hole, on_hole = _line_in_closed_curve(point, hole; close = true)
        # point is in a hole -> not in polygon
        (in_hole || on_hole) && return (false, false)  # TODO: what if all points on the edge of hole?
    end
    return (in_ext, on_ext)  # point is inside of polygon
end

function polygon_in_geom(::GI.PolygonTrait, poly1, ::GI.PolygonTrait, poly2)
    Extents.intersects(GI.extent(poly1), GI.extent(poly2)) || return false

end

function polygon_in_polygon(poly1, poly2)
    # edges1, edges2 = to_edges(poly1), to_edges(poly2)
    # extent1, extent2 = to_extent(edges1), to_extent(edges2)
    # Check the extents intersect
    Extents.intersects(GI.extent(poly1), GI.extent(poly2)) || return false

    # Check all points in poly1 are in poly2
    for point in GI.getpoint(poly1)
        point_in_polygon(point, poly2)[1] || return false
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

function _line_in_closed_curve(line, curve; close = false)
    # Determine number of points in curve and line
    nc = GI.npoint(curve)
    nc -= equals(GI.getpoint(curve, 1), GI.getpoint(curve, nc)) ? 1 : 0
    nl = GI.npoint(line)
    nl -= (close && equals(GI.getpoint(line, 1), GI.getpoint(line, nl))) ? 1 : 0
    # Check to see if first point in line is within curve
    point_in, point_on = point_in_polygon(GI.getpoint(line, 1), curve)
    (point_in || point_on) || return (false, false)  # point is outside curve
    # Check for any intersections between line and curve
    vertex_on = point_on
    l_start_idx = close ? nl : 1
    l_range = close ? 1 : 2
    c_start = GI.getpoint(curve, nc)
    for i in 1:nc
        c_end = GI.getpoint(curve, i)
        l_start = GI.getpoint(line, l_start_idx)
        for j in l_range:nl
            l_end = GI.getpoint(line, j)
            # Check if edges intersect --> line is not within curve
            _line_intersects(
                (l_start, l_end),
                (c_start, c_end),
            ) && return (false, false)
            # Check if either vertex is on the edge of the curve
            if !vertex_on
                vertex_on = point_on_segment(l_start, c_start, c_end) ||
                    point_on_segment(l_end, c_start, c_end)
            end
            l_start = l_end
        end
        c_start = c_end
    end
    return (!vertex_on, vertex_on)
end

"""
    _point_in_extent(p, extent::Extents.Extent)::Bool

Returns true if the point is the bounding box of the extent and false otherwise. 
"""
function _point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y2) = extent.X, extent.Y
    return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
end
