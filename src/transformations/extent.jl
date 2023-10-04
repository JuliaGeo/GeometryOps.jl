"""
    embed_extent(obj)

Recursively wrap the object with a `GeoInterface.Wrappers` geometry,
calculating and adding an `Extents.Extent` to all objects.

This can improve performance when extents need to be checked multiple times.
"""
embed_extent(x) = apply(extent_applicator, AbstractTrait, x)

extent_applicator(x) = extent_applicator(trait(x), x)
extent_applicator(::Nothing, xs::AbstractArray) = embed_extent.(xs)
function extent_applicator(::Union{AbstractCurveTrait,MultiPointTrait}, point) = point
    
function extent_applicator(trait::AbstractGeometryTrait, geom)
    children_with_extents = map(GI.getgeom(geom)) do g
        embed_extent(g)
    end
    wrapper_type = GI.geointerface_geomtype(trait)
    extent = GI.extent(wrapper_type(children_with_extents))
    return wrapper_type(children_with_extents, extent)
end
extent_applicator(::PointTrait, point) = point
extent_applicator(::PointTrait, point) = point
