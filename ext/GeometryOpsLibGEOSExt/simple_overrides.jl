#=
# Simple overrides
This file contains simple overrides for GEOS, essentially only those
functions which have direct counterparts in LG and only
require conversion before calling.
=#
# ## Polygon set operations
# ### Difference
function GO.difference(::GEOS, geom_a, geom_b; target=nothing, calc_extent = false)
    return _wrap(LG.difference(GI.convert(LG, geom_a), GI.convert(LG, geom_b)); crs = GI.crs(geom_a), calc_extent)
end
# ### Union
function GO.union(::GEOS, geom_a, geom_b; target=nothing, calc_extent = false)
    return _wrap(LG.union(GI.convert(LG, geom_a), GI.convert(LG, geom_b)); crs = GI.crs(geom_a), calc_extent)
end
# ### Intersection
function GO.intersection(::GEOS, geom_a, geom_b; target=nothing, calc_extent = false)
    return _wrap(LG.intersection(GI.convert(LG, geom_a), GI.convert(LG, geom_b)); crs = GI.crs(geom_a), calc_extent)
end
# ### Symmetric difference
function GO.symdifference(::GEOS, geom_a, geom_b; target=nothing, calc_extent = false)
    return _wrap(LG.symmetric_difference(GI.convert(LG, geom_a), GI.convert(LG, geom_b)); crs = GI.crs(geom_a), calc_extent)
end

# ## DE-9IM boolean methods
# ### Equals
function GO.equals(::GEOS, geom_a, geom_b)
    return LG.equals(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Disjoint
function GO.disjoint(::GEOS, geom_a, geom_b)
    return LG.disjoint(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Touches
function GO.touches(::GEOS, geom_a, geom_b)
    return LG.touches(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Crosses
function GO.crosses(::GEOS, geom_a, geom_b)
    return LG.crosses(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Within
function GO.within(::GEOS, geom_a, geom_b)
    return LG.within(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Contains
function GO.contains(::GEOS, geom_a, geom_b)
    return LG.contains(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Overlaps
function GO.overlaps(::GEOS, geom_a, geom_b)
    return LG.overlaps(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Covers
function GO.covers(::GEOS, geom_a, geom_b)
    return LG.covers(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### CoveredBy
function GO.coveredby(::GEOS, geom_a, geom_b)
    return LG.coveredby(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Intersects
function GO.intersects(::GEOS, geom_a, geom_b)
    return LG.intersects(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end

# ## Convex hull
function GO.convex_hull(::GEOS, geoms)
    return LG.convexhull(
        LG.MultiPoint(
            collect(
                GO.flatten(
                    x -> GI.convert(LG, x), 
                    GI.PointTrait, 
                    geoms
                )
            )
        )
    )
end