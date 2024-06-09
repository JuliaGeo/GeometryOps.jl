#=
# Simple overrides
This file contains simple overrides for GEOS, essentially only those
functions which have direct counterparts in LG and only
require conversion before calling.
=#
# ## Polygon set operations
# ### Difference
function GO.difference(::GEOS, geom_a, geom_b; target=nothing)
    return LG.difference(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Union
function GO.union(::GEOS, geom_a, geom_b; target=nothing)
    return LG.union(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Intersection
function GO.intersection(::GEOS, geom_a, geom_b; target=nothing)
    return LG.intersection(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
end
# ### Symmetric difference
function GO.symdifference(::GEOS, geom_a, geom_b; target=nothing)
    return LG.symmetric_difference(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
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

