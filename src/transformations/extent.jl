"""
    embed_extent(obj)

Recursively wrap the object with a `GeoInterface.Wrappers` geometry,
calculating and adding an `Extents.Extent` to all objects.

This can improve performance when extents need to be checked multiple times.
"""
embed_extent(x; kw...) = apply(extent_applicator, GI.AbstractTrait, x; kw...)

# We recursively run `embed_extent` so the lowest level
# extent is calculated first and bubbles back up.
# This means we touch each point only once.
extent_applicator(x) = extent_applicator(trait(x), x)
extent_applicator(::Nothing, xs::AbstractArray) = embed_extent.(xs)
function extent_applicator(trait::GI.AbstractGeometryTrait, geom)
    children_with_extents = map(GI.getgeom(geom)) do g
        embed_extent(g)
    end
    wrapper_type = GI.geointerface_geomtype(trait)
    extent = GI.extent(wrapper_type(children_with_extents))
    return wrapper_type(children_with_extents; extent, crs=GI.crs(geom))
end
function extent_applicator(trait::Union{GI.AbstractCurveTrait,GI.MultiPointTrait}, geom)
    wrapper_type = GI.geointerface_geomtype(trait)
    extent = GI.extent(geom)
    return wrapper_type(geom; extent, crs=GI.crs(geom))
end
extent_applicator(::GI.PointTrait, point) = point
