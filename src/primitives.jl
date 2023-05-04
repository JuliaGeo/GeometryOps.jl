"""
    map(f, target::Type{<:AbstractTrait}, obj; crs)

Reconstruct a geometry or feature using the function `f` on the `target` trait.

`f(target_geom) => x` where `x` also has the `target` trait, or an equivalent.

The result is an functionally similar geometry with values depending on `f`

# Flipped point the order in any feature or geometry, or iterables of either:

```juia
import GeoInterface as GI
import GeometryOps as GO
geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

flipped_geom = GO.map(GI.PointTrait, geom) do p
    (GI.y(p), GI.x(p))
end
"""
function map end
# Add dispatch argument for trait
map(f, target::Type{<:GI.AbstractTrait}, geom; kw...) =
    map(f, target, GI.trait(geom), geom; kw...)
# Try to map over iterables
map(f, target::Type, ::Nothing, iterable; kw...) =
    Base.map(x -> Base.map(f, target, x), iterable; kw...)
# Rewrap feature collections
function map(f, target::Type, ::GI.FeatureCollectionTrait, fc; crs=GI.crs(fc))
    features = Base.map(GI.getfeature(fc)) do feature
        map(f, target, feature)
    end 
    return FeatureCollection(features; crs)
end
# Rewrap features
function map(f, target::Type, ::GI.FeatureTrait, feature; crs=GI.crs(feature))
    properties = GI.properties(feature)
    geometry = map(f, target, geometry(feature); crs)
    return Feature(geometry; properties, crs)
end
# Reconstruct nested geometries
function map(f, target::Type, trait, geom; crs=GI.crs(geom))::(GI.geointerface_geomtype(trait))
    # TODO handle zero length...
    geoms = Base.map(g -> map(f, target, g), GI.getgeom(geom))
    if GI.is3d(geom)
        return GI.geointerface_geomtype(trait){true,false}(geoms; crs)
    else
        return GI.geointerface_geomtype(trait){false,false}(geoms; crs)
    end
end
# Apply f to the target geometry
map(f, ::Type{Target}, ::Trait, geom; crs=nothing) where {Target,Trait<:Target} = f(geom)
# Fail if we hit PointTrait without running `f`
map(f, target::Type, trait::GI.PointTrait, geom; crs=nothing) =
    throw(ArgumentError("target $target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
map(f, target::Type{GI.PointTrait}, trait::GI.PointTrait, geom; crs=nothing) = f(geom)
map(f, target::Type{GI.FeatureTrait}, ::GI.FeatureTrait, feature; crs=nothing) = f(feature)
map(f, target::Type{GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc; crs=nothing) = f(fc)

