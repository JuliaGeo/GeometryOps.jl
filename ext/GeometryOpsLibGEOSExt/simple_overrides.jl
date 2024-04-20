#=
# Simple overrides
This file contains simple overrides for GEOS, essentially only those
functions which have direct counterparts in LibGEOS and only
require conversion before calling.
=#
# ## Polygon set operations
# ### Difference
function GO.difference(::GEOS, geom_a, geom_b; target=nothing)
    return LG.difference(geom_a, geom_b)
end
# ### Union
function GO.union(::GEOS, geom_a, geom_b; target=nothing)
    return LG.union(geom_a, geom_b)
end
# ### Intersection
function GO.intersection(::GEOS, geom_a, geom_b; target=nothing)
    return LG.intersection(geom_a, geom_b)
end

# ## DE-9IM boolean methods
# ### Equals
function GO.equals(::GEOS, geom_a, geom_b)
    return LG.equals(geom_a, geom_b)
end
# ### Disjoint
function GO.disjoint(::GEOS, geom_a, geom_b)
    return LG.disjoint(geom_a, geom_b)
end
# ### Touches
function GO.touches(::GEOS, geom_a, geom_b)
    return LG.touches(geom_a, geom_b)
end
# ### Crosses
function GO.crosses(::GEOS, geom_a, geom_b)
    return LG.crosses(geom_a, geom_b)
end
# ### Within
function GO.within(::GEOS, geom_a, geom_b)
    return LG.within(geom_a, geom_b)
end
# ### Contains
function GO.contains(::GEOS, geom_a, geom_b)
    return LG.contains(geom_a, geom_b)
end
# ### Overlaps
function GO.overlaps(::GEOS, geom_a, geom_b)
    return LG.overlaps(geom_a, geom_b)
end
# ### Covers
function GO.covers(::GEOS, geom_a, geom_b)
    return LG.covers(geom_a, geom_b)
end
# ### CoveredBy
function GO.coveredby(::GEOS, geom_a, geom_b)
    return LG.coveredby(geom_a, geom_b)
end