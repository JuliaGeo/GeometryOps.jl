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

If a point touches a linestring if it equal to 
"""
touches(
    ::GI.PointTrait, g1,
    ::GI.LineStringTrait, g2,
) = _point_curve_process(
    g1, g2;
    process = touch_process,
    repeated_last_coord = false,
)

"""    
    touches(::GI.PointTrait, g1, ::GI.LinearRingTrait, g2)::Bool

If a point is disjoint from a linear ring then it is not on any of the
ring's edges or vertices. If these conditions are met, return true, else false.
"""
touches(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    process = touch_process,
    repeated_last_coord = true,
)