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
apply(f, ::Type{Target}, geom; kw...)  where Target<:GI.AbstractTrait =
    _apply(f, Target, geom; kw...)

_apply(f, ::Type{Target}, geom; kw...)  where Target =
    _apply(f, Target, GI.trait(geom), geom; kw...)
# Try to _apply over iterables
_apply(f, ::Type{Target}, ::Nothing, iterable; kw...) where Target =
    map(x -> _apply(f, Target, x), iterable; kw...)
# Rewrap feature collections
function _apply(f, ::Type{Target}, ::GI.FeatureCollectionTrait, fc; crs=GI.crs(fc)) where Target
    features = map(GI.getfeature(fc)) do feature
        _apply(f, Target, feature)
    end 
    return FeatureCollection(features; crs)
end
# Rewrap features
function _apply(f, ::Type{Target}, ::GI.FeatureTrait, feature; crs=GI.crs(feature)) where Target
    properties = GI.properties(feature)
    geometry = _apply(f, Target, geometry(feature); crs)
    return Feature(geometry; properties, crs)
end
# Reconstruct nested geometries
function _apply(f, ::Type{Target}, trait, geom; crs=GI.crs(geom))::(GI.geointerface_geomtype(trait)) where Target
    # TODO handle zero length...
    geoms = map(g -> _apply(f, Target, g), GI.getgeom(geom))
    if GI.is3d(geom)
        # The Boolean type parameters here indicate 3d-ness and measure coordinate presence respectively.
        return GI.geointerface_geomtype(trait){true,false}(geoms; crs)
    else
        return GI.geointerface_geomtype(trait){false,false}(geoms; crs)
    end
end
# Apply f to the target geometry
_apply(f, ::Type{Target}, ::Trait, geom; crs=nothing) where {Target,Trait<:Target} = f(geom)
# Fail if we hit PointTrait without running `f`
_apply(f, ::Type{Target}, trait::GI.PointTrait, geom; crs=nothing) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
_apply(f, ::Type{GI.PointTrait}, trait::GI.PointTrait, geom; crs=nothing) = f(geom)
_apply(f, ::Type{GI.FeatureTrait}, ::GI.FeatureTrait, feature; crs=nothing) = f(feature)
_apply(f, ::Type{GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc; crs=nothing) = f(fc)


"""
    flatten(target::Type{<:GI.AbstractTrait}, geom)

Lazily flatten any geometry, feature or iterator of geometries or features
so that objects with the specified trait are returned by the iterator.
"""
flatten(::Type{Target}, geom) where {Target<:GI.AbstractTrait} = _flatten(Target, geom) 

_flatten(::Type{Target}, geom) where Target = _flatten(Target, GI.trait(geom), geom)
# Try to flatten over iterables
_flatten(::Type{Target}, ::Nothing, iterable) where Target = 
    Iterators.flatmap(x -> _flatten(Target, x), iterable)
# Flatten feature collections
function _flatten(::Type{Target}, ::GI.FeatureCollectionTrait, fc) where Target
    Iterators.flatmap(GI.getfeature(fc)) do feature
        _flatten(Target, feature)
    end 
end
_flatten(::Type{Target}, ::GI.FeatureTrait, feature) where Target = 
    _flatten(Target, geometry(feature))
# Apply f to the target geometry
_flatten(::Type{Target}, ::Trait, geom) where {Target,Trait<:Target} = (geom,)
_flatten(::Type{Target}, trait, geom) where Target = 
    Iterators.flatmap(g -> _flatten(Target, g), GI.getgeom(geom))
# Fail if we hit PointTrait without running `f`
_flatten(::Type{Target}, trait::GI.PointTrait, geom) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
_flatten(::Type{<:GI.PointTrait}, ::GI.PointTrait, geom) = (geom,)
_flatten(::Type{<:GI.FeatureTrait}, ::GI.FeatureTrait, feature) = (feature,)
_flatten(::Type{<:GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc) = (fc,)


"""
    reconstruct(geom, components)

Reconstruct `geom` from an iterable of component objects that match its structure.

All objects in `components` must have the same `GeoInterface.trait`.

Ususally used in combination with `flatten`.
"""
reconstruct(geom, components) = _reconstruct(Target, geom) 

_reconstruct(geom, components) = _reconstruct(typeof(GI.trait(first(components))), geom, components) 
_reconstruct(::Type{Target}, geom) where Target = 
    _reconstruct(Target, GI.trait(geom), geom, components)
# Try to reconstruct over iterables
_reconstruct(::Type{Target}, ::Nothing, iterable) where Target = 
    map(x -> _reconstruct(Target, x), iterable)
# Flatten feature collections
function _reconstruct(::Type{Target}, ::GI.FeatureCollectionTrait, fc, components) where Target
    map(GI.getfeature(fc)) do feature
        _reconstruct(Target, feature, components)
    end 
end
_reconstruct(::Type{Target}, ::GI.FeatureTrait, feature) where Target = 
    _reconstruct(Target, geometry(feature))
# Apply f to the target geometry
_reconstruct(::Type{Target}, ::Trait, geom) where {Target,Trait<:Target} = geom
_reconstruct(::Type{Target}, trait, geom) where Target = 
    Iterators.flatmap(g -> _reconstruct(Target, g), GI.getgeom(geom))
# Fail if we hit PointTrait without running `f`
_reconstruct(::Type{Target}, trait::GI.PointTrait, geom) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
_reconstruct(::Type{<:GI.PointTrait}, ::GI.PointTrait, geom) = geom
_reconstruct(::Type{<:GI.FeatureTrait}, ::GI.FeatureTrait, feature) = feature
_reconstruct(::Type{<:GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc) = fc
