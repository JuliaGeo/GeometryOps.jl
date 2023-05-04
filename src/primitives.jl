"""
    apply(f, target::Type{<:AbstractTrait}, obj; crs)

Reconstruct a geometry or feature using the function `f` on the `target` trait.

`f(target_geom) => x` where `x` also has the `target` trait, or an equivalent.

The result is an functionally similar geometry with values depending on `f`

# Flipped point the order in any feature or geometry, or iterables of either:

```juia
import GeoInterface as GI
import GeometryOps as GO
geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

flipped_geom = GO.apply(GI.PointTrait, geom) do p
    (GI.y(p), GI.x(p))
end
"""
function apply end
# Add dispatch argument for trait
apply(f, target::Type{<:GI.AbstractTrait}, geom; kw...) =
    apply(f, target, GI.trait(geom), geom; kw...)
# Try to apply over iterables
apply(f, target::Type, ::Nothing, iterable; kw...) =
    map(x -> apply(f, target, x), iterable; kw...)
# Rewrap feature collections
function apply(f, target::Type, ::GI.FeatureCollectionTrait, fc; crs=GI.crs(fc))
    features = map(GI.getfeature(fc)) do feature
        apply(f, target, feature)
    end 
    return FeatureCollection(features; crs)
end
# Rewrap features
function apply(f, target::Type, ::GI.FeatureTrait, feature; crs=GI.crs(feature))
    properties = GI.properties(feature)
    geometry = apply(f, target, geometry(feature); crs)
    return Feature(geometry; properties, crs)
end
# Reconstruct nested geometries
function apply(f, target::Type, trait, geom; crs=GI.crs(geom))::(GI.geointerface_geomtype(trait))
    # TODO handle zero length...
    geoms = map(g -> apply(f, target, g), GI.getgeom(geom))
    if GI.is3d(geom)
        return GI.geointerface_geomtype(trait){true,false}(geoms; crs)
    els
        return GI.geointerface_geomtype(trait){false,false}(geoms; crs)
    end
end
# Apply f to the target geometry
apply(f, ::Type{Target}, ::Trait, geom; crs=nothing) where {Target,Trait<:Target} = f(geom)
# Fail if we hit PointTrait without running `f`
apply(f, target::Type, trait::GI.PointTrait, geom; crs=nothing) =
    throw(ArgumentError("target $target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
apply(f, target::Type{GI.PointTrait}, trait::GI.PointTrait, geom; crs=nothing) = f(geom)
apply(f, target::Type{GI.FeatureTrait}, ::GI.FeatureTrait, feature; crs=nothing) = f(feature)
apply(f, target::Type{GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc; crs=nothing) = f(fc)
