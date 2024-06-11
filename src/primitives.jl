# # Primitive functions

export apply, applyreduce, TraitTarget

#=

This file mainly defines the [`apply`](@ref) and [`applyreduce`](@ref) functions, and some related functionality.

In general, the idea behind the `apply` framework is to take 
as input any geometry, vector of geometries, or feature collection,
deconstruct it to the given trait target (any arbitrary GI.AbstractTrait 
or `TraitTarget` union thereof, like `PointTrait` or `PolygonTrait`) 
and perform some operation on it.  

This allows for a simple and consistent framework within which users can 
define their own operations trivially easily, and removes a lot of the 
complexity involved with handling complex geometry structures.

For example, a simple way to flip the x and y coordinates of a geometry is:

```julia
flipped_geom = GO.apply(GI.PointTrait(), geom) do p
    (GI.y(p), GI.x(p))
end
```

As simple as that.  There's no need to implement your own decomposition because it's done for you.

Functions like [`flip`](@ref), [`reproject`](@ref), [`transform`](@ref), even [`segmentize`](@ref) and [`simplify`](@ref) have been implemented
using the `apply` framework.  Similarly, [`centroid`](@ref), [`area`](@ref) and [`distance`](@ref) have been implemented using the 
[`applyreduce`](@ref) framework.

## Docstrings

### Functions

```@docs
apply
applyreduce
GeometryOps.unwrap
GeometryOps.flatten
GeometryOps.reconstruct
GeometryOps.rebuild
```

## Types

```@docs
TraitTarget
```

## Implementation

=#

const TuplePoint{T} = Tuple{T, T} where T <: AbstractFloat
const TupleEdge{T} = Tuple{TuplePoint{T},TuplePoint{T}} where T

struct SVPoint{N,T,Z,M} <: GeometryBasics.StaticArraysCore.StaticVector{N,T}
    vals::NTuple{N,T}  # TODO: Should Z and M be booleans or BoolsAsTypes?
end

Base.getindex(p::SVPoint, i::Int64) = p.vals[i]

SVPoint_2D(vals::NTuple{2,T}) where T = SVPoint{2, T, _False, _False}(vals)
SVPoint_2D(vals) = SVPoint_2D((GI.x(vals), GI.y(vals)))

SVPoint_3D(vals::NTuple{3,T}) where T = SVPoint{3, T, _True, _False}(vals)
SVPoint_3D(vals) = SVPoint_3D((GI.x(vals), GI.y(vals), GI.z(vals)))

SVPoint_4D(vals::NTuple{4,T}) where T = SVPoint{4, T, _True, _True}(vals)
SVPoint_4D(vals) = SVPoint_4D((GI.x(vals), GI.y(vals), GI.z(vals), GI.m(vals)))

SVPoint(vals::NTuple{2, <:AbstractFloat}) = SVPoint_2D(vals)
SVPoint(vals::NTuple{3, <:AbstractFloat}) = SVPoint_3D(vals)
SVPoint(vals::NTuple{4, <:AbstractFloat}) = SVPoint_4D(vals)

const SVEdge{T} = Tuple{SVPoint{N,T,Z,M}, SVPoint{N,T,Z,M}} where {N,T,Z,M}
const Edge{T} = Union{TupleEdge{T}, SVEdge{T}} where T

const THREADED_KEYWORD = "- `threaded`: `true` or `false`. Whether to use multithreading. Defaults to `false`."
const CRS_KEYWORD = "- `crs`: The CRS to attach to geometries. Defaults to `nothing`."
const CALC_EXTENT_KEYWORD = "- `calc_extent`: `true` or `false`. Whether to calculate the extent. Defaults to `false`."

const APPLY_KEYWORDS = """
$THREADED_KEYWORD
$CRS_KEYWORD
$CALC_EXTENT_KEYWORD
"""

#=
## What is `apply`?

`apply` applies some function to every geometry matching the `Target`
GeoInterface trait, in some arbitrarily nested object made up of:
- `AbstractArray`s (we also try to iterate other non-GeoInteface compatible object)
- `FeatureCollectionTrait` objects
- `FeatureTrait` objects
- `AbstractGeometryTrait` objects

`apply` recursively calls itself through these nested
layers until it reaches objects with the `Target` GeoInterface trait. When found `apply` applies the function `f`, and stops.

The outer recursive functions then progressively rebuild the object
using GeoInterface objects matching the original traits.

If `PointTrait` is found  but it is not the `Target`, an error is thrown.
This likely means the object contains a different geometry trait to
the target, such as `MultiPointTrait` when `LineStringTrait` was specified.

To handle this possibility it may be necessary to make `Target` a
`Union` of traits found at the same level of nesting, and define methods
of `f` to handle all cases.

Be careful making a union across "levels" of nesting, e.g.
`Union{FeatureTrait,PolygonTrait}`, as `_apply` will just never reach
`PolygonTrait` when all the polygons are wrapped in a `FeatureTrait` object.

## Embedding:

`extent` and `crs` can be embedded in all geometries, features, and
feature collections as part of `apply`. Geometries deeper than `Target`
will of course not have new `extent` or `crs` embedded.

- `calc_extent` signals to recalculate an `Extent` and embed it.
- `crs` will be embedded as-is

## Threading

Threading is used at the outermost level possible - over
an array, feature collection, or e.g. a MultiPolygonTrait where
each `PolygonTrait` sub-geometry may be calculated on a different thread.

Currently, threading defaults to `false` for all objects, but can be turned on
by passing the keyword argument `threaded=true` to `apply`.
=#

"""
    apply(f, target::Union{TraitTarget, GI.AbstractTrait}, obj; kw...)

Reconstruct a geometry, feature, feature collection, or nested vectors of
either using the function `f` on the `target` trait.

`f(target_geom) => x` where `x` also has the `target` trait, or a trait that can
be substituted. For example, swapping `PolgonTrait` to `MultiPointTrait` will fail
if the outer object has `MultiPolygonTrait`, but should work if it has `FeatureTrait`.

Objects "shallower" than the target trait are always completely rebuilt, like
a `Vector` of `FeatureCollectionTrait` of `FeatureTrait` when the target
has `PolygonTrait` and is held in the features. These will always be GeoInterface 
geometries/feature/feature collections. But "deeper" objects may remain
unchanged or be whatever GeoInterface compatible objects `f` returns.

The result is a functionally similar geometry with values depending on `f`.

$APPLY_KEYWORDS

## Example

Flipped point the order in any feature or geometry, or iterables of either:

```julia
import GeoInterface as GI
import GeometryOps as GO
geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]),
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

flipped_geom = GO.apply(GI.PointTrait, geom) do p
    (GI.y(p), GI.x(p))
end
```
"""
@inline function apply(
    f::F, target, geom; calc_extent=false, threaded=false, kw...
) where F
    threaded = _booltype(threaded)
    calc_extent = _booltype(calc_extent)
    _apply(f, TraitTarget(target), geom; threaded, calc_extent, kw...)
end

# Call _apply again with the trait of `geom`
@inline _apply(f::F, target, geom; kw...)  where F =
    _apply(f, target, GI.trait(geom), geom; kw...)
# There is no trait and this is an AbstractArray - so just iterate over it calling _apply on the contents
@inline function _apply(f::F, target, ::Nothing, A::AbstractArray; threaded, kw...) where F
    # For an Array there is nothing else to do but map `_apply` over all values
    # _maptasks may run this level threaded if `threaded==true`,
    # but deeper `_apply` called in the closure will not be threaded
    apply_to_array(i) = _apply(f, target, A[i]; threaded=_False(), kw...)
    _maptasks(apply_to_array, eachindex(A), threaded)
end
# There is no trait and this is not an AbstractArray.
# Try to call _apply over it. We can't use threading
# as we don't know if we can can index into it. So just `map`.
@inline function _apply(f::F, target, ::Nothing, iterable::IterableType; threaded, kw...) where {F, IterableType}
    # Try the Tables.jl interface first
    if Tables.istable(iterable)
    _apply_table(f, target, iterable; threaded, kw...)
    else # this is probably some form of iterable...
        if threaded isa _True
            # `collect` first so we can use threads
            _apply(f, target, collect(iterable); threaded, kw...)
        else
            apply_to_iterable(x) = _apply(f, target, x; kw...)
            map(apply_to_iterable, iterable)
        end
    end
end
#= 
Doing this inline in `_apply` is _heavily_ type unstable, so it's best to separate this 
by a function barrier.

This function operates `apply` on the `geometry` column of the table, and returns a new table
with the same schema, but with the new geometry column.

This new table may be of the same type as the old one iff `Tables.materializer` is defined for 
that table.  If not, then a `NamedTuple` is returned.
=#
function _apply_table(f::F, target, iterable::IterableType; threaded, kw...) where {F, IterableType}
    _get_col_pair(colname) = colname => Tables.getcolumn(iterable, colname)
    # We extract the geometry column and run `apply` on it.
    geometry_column = first(GI.geometrycolumns(iterable))
    new_geometry = _apply(f, target, Tables.getcolumn(iterable, geometry_column); threaded, kw...)
    # Then, we obtain the schema of the table,
    old_schema = Tables.schema(iterable)
    # filter the geometry column out,
    new_names = filter(Base.Fix1(!==, geometry_column), old_schema.names)
    # and try to rebuild the same table as the best type - either the original type of `iterable`,
    # or a named tuple which is the default fallback.
    return Tables.materializer(iterable)(
        merge(
            NamedTuple{(geometry_column,), Base.Tuple{typeof(new_geometry)}}((new_geometry,)),
            NamedTuple(Iterators.map(_get_col_pair, new_names))
        )
    )
end

# Rewrap all FeatureCollectionTrait feature collections as GI.FeatureCollection
# Maybe use threads to call _apply on component features
@inline function _apply(f::F, target, ::GI.FeatureCollectionTrait, fc;
    crs=GI.crs(fc), calc_extent=_False(), threaded
) where F

    # Run _apply on all `features` in the feature collection, possibly threaded
    apply_to_feature(i) =
        _apply(f, target, GI.getfeature(fc, i); crs, calc_extent, threaded=_False())::GI.Feature
    features = _maptasks(apply_to_feature, 1:GI.nfeature(fc), threaded)
    if calc_extent isa _True
        # Calculate the extent of the features
        extent = mapreduce(GI.extent, Extents.union, features)
        # Return a FeatureCollection with features, crs and caculated extent
        return GI.FeatureCollection(features; crs, extent)
    else
        # Return a FeatureCollection with features and crs
        return GI.FeatureCollection(features; crs)
    end
end
# Rewrap all FeatureTrait features as GI.Feature, keeping the properties
@inline function _apply(f::F, target, ::GI.FeatureTrait, feature;
    crs=GI.crs(feature), calc_extent=_False(), threaded
) where F
    # Run _apply on the contained geometry
    geometry = _apply(f, target, GI.geometry(feature); crs, calc_extent, threaded)
    # Get the feature properties
    properties = GI.properties(feature)
    if calc_extent isa _True
        # Calculate the extent of the geometry
        extent = GI.extent(geometry)
        # Return a new Feature with the new geometry and calculated extent, but the oroginal properties and crs
        return GI.Feature(geometry; properties, crs, extent)
    else
        # Return a new Feature with the new geometry, but the oroginal properties and crs
        return GI.Feature(geometry; properties, crs)
    end
end
# Reconstruct nested geometries,
# maybe using threads to call _apply on component geoms
@inline function _apply(f::F, target, trait, geom;
    crs=GI.crs(geom), calc_extent=_False(), threaded
)::(GI.geointerface_geomtype(trait)) where F
    # Map `_apply` over all sub geometries of `geom`
    # to create a new vector of geometries
    # TODO handle zero length
    apply_to_geom(i) = _apply(f, target, GI.getgeom(geom, i); crs, calc_extent, threaded=_False())
    geoms = _maptasks(apply_to_geom, 1:GI.ngeom(geom), threaded)
    return _apply_inner(geom, geoms, crs, calc_extent)
end
function _apply_inner(geom, geoms, crs, calc_extent::_True)
    # Calculate the extent of the sub geometries
    extent = mapreduce(GI.extent, Extents.union, geoms)
    # Return a new geometry of the same trait as `geom`,
    # holding tnew `geoms` with `crs` and calcualted extent
    return rebuild(geom, geoms; crs, extent)
end
function _apply_inner(geom, geoms, crs, calc_extent::_False)
    # Return a new geometryof the same trait as `geom`, holding the new `geoms` with `crs`
    return rebuild(geom, geoms; crs)
end
# Fail loudly if we hit PointTrait without running `f`
# (after PointTrait there is no further to dig with `_apply`)
# @inline _apply(f, ::TraitTarget{Target}, trait::GI.PointTrait, geom; crs=nothing, kw...) where Target =
    # throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Finally, these short methods are the main purpose of `apply`.
# The `Trait` is a subtype of the `Target` (or identical to it)
# So the `Target` is found. We apply `f` to geom and return it to previous
# _apply calls to be wrapped with the outer geometries/feature/featurecollection/array.
_apply(f::F, ::TraitTarget{Target}, ::Trait, geom; crs=GI.crs(geom), kw...) where {F,Target,Trait<:Target} = f(geom)
# Define some specific cases of this match to avoid method ambiguity
for T in (
    GI.PointTrait, GI.LinearRing, GI.LineString,
    GI.MultiPoint, GI.FeatureTrait, GI.FeatureCollectionTrait
)
    @eval _apply(f::F, target::TraitTarget{<:$T}, trait::$T, x; kw...) where F = f(x)
end

"""
    applyreduce(f, op, target::Union{TraitTarget, GI.AbstractTrait}, obj; threaded)

Apply function `f` to all objects with the `target` trait,
and reduce the result with an `op` like `+`. 

The order and grouping of application of `op` is not guaranteed.

If `threaded==true` threads will be used over arrays and iterables, 
feature collections and nested geometries.
"""
@inline function applyreduce(
    f::F, op::O, target, geom; threaded=false, init=nothing
) where {F, O}
    threaded = _booltype(threaded)
    _applyreduce(f, op, TraitTarget(target), geom; threaded, init)
end

@inline _applyreduce(f::F, op::O, target, geom; threaded, init) where {F, O} =
    _applyreduce(f, op, target, GI.trait(geom), geom; threaded, init)
# Maybe use threads recucing over arrays
@inline function _applyreduce(f::F, op::O, target, ::Nothing, A::AbstractArray; threaded, init) where {F, O}
    applyreduce_array(i) = _applyreduce(f, op, target, A[i]; threaded=_False(), init)
    _mapreducetasks(applyreduce_array, op, eachindex(A), threaded; init)
end
# Try to applyreduce over iterables
@inline function _applyreduce(f::F, op::O, target, ::Nothing, iterable::IterableType; threaded, init) where {F, O, IterableType}
    if Tables.istable(iterable)
        _applyreduce_table(f, op, target, iterable; threaded, init)
    else
        applyreduce_iterable(i) = _applyreduce(f, op, target, i; threaded=_False(), init)
        if threaded isa _True # Try to `collect` and reduce over the vector with threads
            _applyreduce(f, op, target, collect(iterable); threaded, init)
        else
            # Try to `mapreduce` the iterable as-is
            mapreduce(applyreduce_iterable, op, iterable; init)
        end
    end
end
# In this case, we don't reconstruct the table, but only operate on the geometry column.
function _applyreduce_table(f::F, op::O, target, iterable::IterableType; threaded, init) where {F, O, IterableType}
    # We extract the geometry column and run `applyreduce` on it.
    geometry_column = first(GI.geometrycolumns(iterable))
    return _applyreduce(f, op, target, Tables.getcolumn(iterable, geometry_column); threaded, init)
end
# If `applyreduce` wants features, then applyreduce over the rows as `GI.Feature`s.
function _applyreduce_table(f::F, op::O, target::GI.FeatureTrait, iterable::IterableType; threaded, init) where {F, O, IterableType}
    # We extract the geometry column and run `apply` on it.
    geometry_column = first(GI.geometrycolumns(iterable))
    property_names = Iterators.filter(!=(geometry_column), Tables.schema(iterable).names)
    features = map(Tables.rows(iterable)) do row
        GI.Feature(Tables.getcolumn(row, geometry_column), properties=NamedTuple(Iterators.map(Base.Fix1(_get_col_pair, row), property_names)))
    end
    return _applyreduce(f, op, target, features; threaded, init)
end
# Maybe use threads reducing over features of feature collections
@inline function _applyreduce(f::F, op::O, target, ::GI.FeatureCollectionTrait, fc; threaded, init) where {F, O}
    applyreduce_fc(i) = _applyreduce(f, op, target, GI.getfeature(fc, i); threaded=_False(), init)
    _mapreducetasks(applyreduce_fc, op, 1:GI.nfeature(fc), threaded; init)
end
# Features just applyreduce to their geometry
@inline _applyreduce(f::F, op::O, target, ::GI.FeatureTrait, feature; threaded, init) where {F, O} =
    _applyreduce(f, op, target, GI.geometry(feature); threaded, init)
# Maybe use threads over components of nested geometries
@inline function _applyreduce(f::F, op::O, target, trait, geom; threaded, init) where {F, O}
    applyreduce_geom(i) = _applyreduce(f, op, target, GI.getgeom(geom, i); threaded=_False(), init)
    _mapreducetasks(applyreduce_geom, op, 1:GI.ngeom(geom), threaded; init)
end
# Don't thread over points it won't pay off
@inline function _applyreduce(
    f::F, op::O, target, trait::Union{GI.LinearRing,GI.LineString,GI.MultiPoint}, geom;
    threaded, init
) where {F, O}
    _applyreduce(f, op, target, GI.getgeom(geom); threaded=_False(), init)
end
# Apply f to the target
@inline function _applyreduce(f::F, op::O, ::TraitTarget{Target}, ::Trait, x; kw...) where {F,O,Target,Trait<:Target} 
    f(x)
end
# Fail if we hit PointTrait
# _applyreduce(f, op, target::TraitTarget{Target}, trait::PointTrait, geom; kw...) where Target = 
    # throw(ArgumentError("target $target not found"))
# Specific cases to avoid method ambiguity
for T in (
    GI.PointTrait, GI.LinearRing, GI.LineString, 
    GI.MultiPoint, GI.FeatureTrait, GI.FeatureCollectionTrait
)
    @eval _applyreduce(f::F, op::O, ::TraitTarget{<:$T}, trait::$T, x; kw...) where {F, O} = f(x)
end

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
_flatten(f, ::Type{Target}, ::Nothing, iterable) where Target =
    Iterators.flatten(Iterators.map(x -> _flatten(f, Target, x), iterable))
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

Ususally used in combination with `flatten`.
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


const BasicsGeoms = Union{GB.AbstractGeometry,GB.AbstractFace,GB.AbstractPoint,GB.AbstractMesh,
    GB.AbstractPolygon,GB.LineString,GB.MultiPoint,GB.MultiLineString,GB.MultiPolygon,GB.Mesh}

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

using Base.Threads: nthreads, @threads, @spawn


# Threading utility, modified Mason Protters threading PSA
# run `f` over ntasks, where f recieves an AbstractArray/range
# of linear indices
@inline function _maptasks(f::F, taskrange, threaded::_True)::Vector where F
    ntasks = length(taskrange)
    # Customize this as needed.
    # More tasks have more overhead, but better load balancing
    tasks_per_thread = 2
    chunk_size = max(1, ntasks ÷ (tasks_per_thread * nthreads()))
    # partition the range into chunks
    task_chunks = Iterators.partition(taskrange, chunk_size)
    # Map over the chunks
    tasks = map(task_chunks) do chunk
        # Spawn a task to process this chunk
        @spawn begin
            # Where we map `f` over the chunk indices
            map(f, chunk)
        end
    end

    # Finally we join the results into a new vector
    return mapreduce(fetch, vcat, tasks)
end
#=
Here we use the compiler directive `@assume_effects :foldable` to force the compiler
to lookup through the closure. This alone makes e.g. `flip` 2.5x faster!
=#
Base.@assume_effects :foldable @inline function _maptasks(f::F, taskrange, threaded::_False)::Vector where F
    map(f, taskrange)
end

# Threading utility, modified Mason Protters threading PSA
# run `f` over ntasks, where f recieves an AbstractArray/range
# of linear indices
#
# WARNING: this will not work for mean/median - only ops
# where grouping is possible
@inline function _mapreducetasks(f::F, op, taskrange, threaded::_True; init) where F
    ntasks = length(taskrange)
    # Customize this as needed.
    # More tasks have more overhead, but better load balancing
    tasks_per_thread = 2
    chunk_size = max(1, ntasks ÷ (tasks_per_thread * nthreads()))
    # partition the range into chunks
    task_chunks = Iterators.partition(taskrange, chunk_size)
    # Map over the chunks
    tasks = map(task_chunks) do chunk
        # Spawn a task to process this chunk
        @spawn begin
            # Where we map `f` over the chunk indices
            mapreduce(f, op, chunk; init)
        end
    end

    # Finally we join the results into a new vector
    return mapreduce(fetch, op, tasks; init)
end
Base.@assume_effects :foldable function _mapreducetasks(f::F, op, taskrange, threaded::_False; init) where F
    mapreduce(f, op, taskrange; init)
end
