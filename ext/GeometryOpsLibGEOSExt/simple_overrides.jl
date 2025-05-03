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
# These are all the same so we loop over all names and eval them in
for fn in (:equals, :disjoint, :touches, :crosses, :within, :contains, :overlaps, :covers, :coveredby, :intersects)
    @eval begin
        # The basic method for geometries
        function GO.$fn(::GEOS, geom_a, geom_b)
            return LG.$fn(GI.convert(LG, geom_a), GI.convert(LG, geom_b))
        end
        # Extents and geometries
        function GO.$fn(alg::GEOS, geom_a::GO.Extents.Extent, geom_b)
            return GO.$fn(alg, GO.extent_to_polygon(geom_a), geom_b)
        end
        function GO.$fn(alg::GEOS, geom_a, geom_b::GO.Extents.Extent)
            return GO.$fn(alg, geom_a, GO.extent_to_polygon(geom_b))
        end
        # Pure extents - this should probably be some GEOSRect or something,
        # but for now this works
        function GO.$fn(alg::GEOS, geom_a::GO.Extents.Extent, geom_b::GO.Extents.Extent)
            return GO.$fn(alg, GO.extent_to_polygon(geom_a), GO.extent_to_polygon(geom_b))
        end
    end
end
# ## Convex hull
function GO.convex_hull(::GEOS, geoms)
    chull = LG.convexhull(
        LG.MultiPoint(
            collect(
                GO.flatten(
                    x -> GI.convert(LG.Point, x), 
                    GI.PointTrait, 
                    geoms
                )
            )
        )
    ); 
    return _wrap(
        chull;
        crs = GI.crs(geoms), 
        calc_extent = false
    )
end
