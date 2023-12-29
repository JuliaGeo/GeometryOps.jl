"""
they have at least one point in common, but their interiors do not intersect.
"""
touches(g1, g2)::Bool = touches(trait(g1), g1, trait(g2), g2)

"""
    touches(::GI.PointTrait, g1, ::GI.PointTrait, g2)::Bool

Two points cannot touch. If they are the same point then their interiors
intersect and if they are different points then they don't share any points.
"""
touches(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = false

"""
    touches(::GI.PointTrait, g1, ::GI.LineStringTrait, g2)::Bool

If a point touches a linestring if it equal to either the first of last point of
the linestring, which make up the linestrings boundaries. If the first and last
point are equal, closing the linestring, then no point can touch the linestring.
"""
function touches(
    ::GI.PointTrait, g1,
    ::GI.LineStringTrait, g2,
)
    n = GI.npoint(g2)
    p1 = GI.getpoint(g2, 1)
    pn = GI.getpoint(g2, n)
    equals(p1, pn) && return false
    return equals(g1, p1) || equals(g1, pn)
end

"""    
    touches(::GI.PointTrait, g1, ::GI.LinearRingTrait, g2)::Bool

If a point cannot 'touch' a linear ring given that the linear ring has no
boundary points. Since the whole ring is "interior", a point cannot touch it.
"""
touches(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = false

"""
    touches(::GI.PointTrait, g1, ::GI.PolygonTrait, g2)::Bool

A point touches a polygon if it is on the boundary of that polygon.
Return true if those conditions are met, else false.
"""
touches(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    in_allow = false, on_allow = true, out_allow = false,
)

"""
    touches(trait1::GI.AbstractTrait, g1, trait2::GI.PointTrait, g2)::Bool

To check if a geometry is touches by a point, switch the order of the
arguments to take advantage of point-geometry touches methods.
"""
touches(
    trait1::GI.AbstractGeometryTrait, g1,
    trait2::GI.PointTrait, g2,
) = touches(trait2, g2, trait1, g1)

# Lines touching geometries
"""
    touches(::GI.LineStringTrait, g1, ::GI.LineStringTrait, g2)::Bool

A line string touches another linestring only if at least one endpoints
(boundary point) of one of the linestrings intersects with the other linestring
"""
touches(
    ::GI.LineStringTrait, g1,
    ::GI.LineStringTrait, g2,
) = _line_curve_process(
    g1, g2;
    over_allow = false, cross_allow = false, on_allow = true, out_allow = true,
    in_require = false, on_require = true, out_require = false,
    closed_line = false,
    closed_curve = false,
)

"""
    touches(::GI.LineStringTrait, g1, ::GI.LinearRingTrait, g2)::Bool

A linestring touches a linear ring if the vertices and edges of the
linestring are touches the linear ring. Return true if those conditions are met,
else false.
"""
touches(
    ::GI.LineStringTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    over_allow = false, cross_allow = false, on_allow = true, out_allow = true,
    in_require = false, on_require = true, out_require = false,
    closed_line = false,
    closed_curve = true,
)
