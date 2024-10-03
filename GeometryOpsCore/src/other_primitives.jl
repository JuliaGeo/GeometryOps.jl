#=
# Other primitives (unwrap, flatten, etc)

This file defines the following primitives:
```@docs; collapsed=true, canonical=false
unwrap
flatten
reconstruct
rebuild
```
=#

"""
    unwrap(target::Type{<:AbstractTrait}, obj)
    unwrap(f, target::Type{<:AbstractTrait}, obj)

Unwrap the object to vectors, down to the target trait.

If `f` is passed in it will be applied to the target geometries
as they are found.
"""
function unwrap end
unwrap(target::Type, geom) = unwrap(identity, target, geom)
# Add dispatch argument for trait
unwrap(f, target::Type, geom) = unwrap(f, target, GI.trait(geom), geom)
# Try to unwrap over iterables
unwrap(f, target::Type, ::Nothing, iterable) =
    map(x -> unwrap(f, target, x), iterable)
# Rewrap feature collections
unwrap(f, target::Type, ::GI.FeatureCollectionTrait, fc) =
    map(x -> unwrap(f, target, x), GI.getfeature(fc))
unwrap(f, target::Type, ::GI.FeatureTrait, feature) =
    unwrap(f, target, GI.geometry(feature))
unwrap(f, target::Type, trait, geom) = map(g -> unwrap(f, target, g), GI.getgeom(geom))
# Apply f to the target geometry
unwrap(f, ::Type{Target}, ::Trait, geom) where {Target,Trait<:Target} = f(geom)
# Fail if we hit PointTrait
unwrap(f, target::Type, trait::GI.PointTrait, geom) =
    throw(ArgumentError("target $target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
unwrap(f, target::Type{GI.PointTrait}, trait::GI.PointTrait, geom) = f(geom)
unwrap(f, target::Type{GI.FeatureTrait}, ::GI.FeatureTrait, feature) = f(feature)
unwrap(f, target::Type{GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc) = f(fc)

"""
    flatten(target::Type{<:GI.AbstractTrait}, obj)
    flatten(f, target::Type{<:GI.AbstractTrait}, obj)

Lazily flatten any `AbstractArray`, iterator, `FeatureCollectionTrait`,
`FeatureTrait` or `AbstractGeometryTrait` object `obj`, so that objects
with the `target` trait are returned by the iterator.

If `f` is passed in it will be applied to the target geometries.
"""
flatten(::Type{Target}, geom) where {Target<:GI.AbstractTrait} = flatten(identity, Target, geom)
flatten(f, ::Type{Target}, geom) where {Target<:GI.AbstractTrait} = _flatten(f, Target, geom)

_flatten(f, ::Type{Target}, geom) where Target = _flatten(f, Target, GI.trait(geom), geom)
# Try to flatten over iterables
function _flatten(f, ::Type{Target}, ::Nothing, iterable) where Target
    if Tables.istable(iterable)
        column = Tables.getcolumn(iterable, first(GI.geometrycolumns(iterable)))
        Iterators.map(x -> _flatten(f, Target, x), column) |> Iterators.flatten
    else
        Iterators.map(x -> _flatten(f, Target, x), iterable) |> Iterators.flatten
    end
end
# Flatten feature collections
function _flatten(f, ::Type{Target}, ::GI.FeatureCollectionTrait, fc) where Target
    Iterators.map(GI.getfeature(fc)) do feature
        _flatten(f, Target, feature)
    end |> Iterators.flatten
end
_flatten(f, ::Type{Target}, ::GI.FeatureTrait, feature) where Target =
    _flatten(f, Target, GI.geometry(feature))
# Apply f to the target geometry
_flatten(f, ::Type{Target}, ::Trait, geom) where {Target,Trait<:Target} = (f(geom),)
_flatten(f, ::Type{Target}, trait, geom) where Target =
    Iterators.flatten(Iterators.map(g -> _flatten(f, Target, g), GI.getgeom(geom)))
# Fail if we hit PointTrait without running `f`
_flatten(f, ::Type{Target}, trait::GI.PointTrait, geom) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
_flatten(f, ::Type{<:GI.PointTrait}, ::GI.PointTrait, geom) = (f(geom),)
_flatten(f, ::Type{<:GI.FeatureTrait}, ::GI.FeatureTrait, feature) = (f(feature),)
_flatten(f, ::Type{<:GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc) = (f(fc),)


"""
    reconstruct(geom, components)

Reconstruct `geom` from an iterable of component objects that match its structure.

All objects in `components` must have the same `GeoInterface.trait`.

Usually used in combination with `flatten`.
"""
function reconstruct(geom, components)
    obj, iter = _reconstruct(geom, components)
    return obj
end

_reconstruct(geom, components) =
    _reconstruct(typeof(GI.trait(first(components))), geom, components, 1)
_reconstruct(::Type{Target}, geom, components, iter) where Target =
    _reconstruct(Target, GI.trait(geom), geom, components, iter)
# Try to reconstruct over iterables
function _reconstruct(::Type{Target}, ::Nothing, iterable, components, iter) where Target
    vect = map(iterable) do x
        # iter is updated by _reconstruct here
        obj, iter = _reconstruct(Target, x, components, iter)
        obj
    end
    return vect, iter
end
# Reconstruct feature collections
function _reconstruct(::Type{Target}, ::GI.FeatureCollectionTrait, fc, components, iter) where Target
    features = map(GI.getfeature(fc)) do feature
        # iter is updated by _reconstruct here
        newfeature, iter = _reconstruct(Target, feature, components, iter)
        newfeature
    end
    return GI.FeatureCollection(features; crs=GI.crs(fc)), iter
end
function _reconstruct(::Type{Target}, ::GI.FeatureTrait, feature, components, iter) where Target
    geom, iter = _reconstruct(Target, GI.geometry(feature), components, iter)
    return GI.Feature(geom; properties=GI.properties(feature), crs=GI.crs(feature)), iter
end
function _reconstruct(::Type{Target}, trait, geom, components, iter) where Target
    geoms = map(GI.getgeom(geom)) do subgeom
        # iter is updated by _reconstruct here
        subgeom1, iter = _reconstruct(Target, GI.trait(subgeom), subgeom, components, iter)
        subgeom1
    end
    return rebuild(geom, geoms), iter
end
# Apply f to the target geometry
_reconstruct(::Type{Target}, ::Trait, geom, components, iter) where {Target,Trait<:Target} =
    iterate(components, iter)
# Specific cases to avoid method ambiguity
_reconstruct(::Type{<:GI.PointTrait}, ::GI.PointTrait, geom, components, iter) = iterate(components, iter)
_reconstruct(::Type{<:GI.FeatureTrait}, ::GI.FeatureTrait, feature, components, iter) = iterate(feature, iter)
_reconstruct(::Type{<:GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc, components, iter) = iterate(fc, iter)
# Fail if we hit PointTrait without running `f`
_reconstruct(::Type{Target}, trait::GI.PointTrait, geom, components, iter) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))

"""
    rebuild(geom, child_geoms)

Rebuild a geometry from child geometries.

By default geometries will be rebuilt as a `GeoInterface.Wrappers`
geometry, but `rebuild` can have methods added to it to dispatch
on geometries from other packages and specify how to rebuild them.

(Maybe it should go into GeoInterface.jl)
"""
rebuild(geom, child_geoms; kw...) = rebuild(GI.trait(geom), geom, child_geoms; kw...)
function rebuild(trait::GI.AbstractTrait, geom, child_geoms; crs=GI.crs(geom), extent=nothing)
    T = GI.geointerface_geomtype(trait)
    if GI.is3d(geom)
        # The Boolean type parameters here indicate "3d-ness" and "measure" coordinate, respectively.
        return T{true,false}(child_geoms; crs, extent)
    else
        return T{false,false}(child_geoms; crs, extent)
    end
end
