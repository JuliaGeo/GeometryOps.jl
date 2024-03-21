# # Primitive functions

# This file mainly defines the [`apply`](@ref) function and its relatives.

#=
We pass `threading` and `calc_extent` as types, not simple boolean values.  

This is to help compilation - with a type to hold on to, it's easier for 
the compiler to separate threaded and non-threaded code paths.
=#
abstract type BoolsAsTypes end
struct _True <: BoolsAsTypes end
struct _False <: BoolsAsTypes end

# This struct holds a trait parameter or a union of trait parameters.
struct TraitTarget{T} end
TraitTarget(::Type{T}) where T = TraitTarget{T}()
TraitTarget(::T) where T<:GI.AbstractTrait = TraitTarget{T}()
TraitTarget(::TraitTarget{T}) where T = TraitTarget{T}()
TraitTarget(::Type{<:TraitTarget{T}}) where T = TraitTarget{T}()
TraitTarget(traits::GI.AbstractTrait...) where T = TraitTarget{Union{map(typeof, traits)...}}()

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
    apply(f, target::Type{<:AbstractTrait}, obj; kw...)

Reconstruct a geometry, feature, feature collection, or nested vectors of
either using the function `f` on the `target` trait.

`f(target_geom) => x` where `x` also has the `target` trait, or a trait that can
be substituted. For example, swapping `PolgonTrait` to `MultiPointTrait` will fail
if the outer object has `MultiPolygonTrait`, but should work if it has `FeatureTrait`.

Objects "shallower" than the target trait are always completely rebuilt, like
a `Vector` of `FeatureCollectionTrait` of `FeatureTrait` when the target
has `PolygonTrait` and is held in the features. But "deeper" objects may remain
unchanged - such as points and linear rings if the target is the same `PolygonTrait`.

The result is a functionally similar geometry with values depending on `f`

$APPLY_KEYWORDS

# Example

Flipped point the order in any feature or geometry, or iterables of either:

```julia
import GeoInterface as GI
import GeometryOps as GO
geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]),
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

flipped_geom = GO.apply(GI.PointTrait, geom) do p
    (GI.y(p), GI.x(p))
end
"""
@inline function apply(
    f::F, target, geom; calc_extent=false, threaded=false, kw...
) where F
    threaded = _booltype(threaded)
    calc_extent = _booltype(calc_extent)
    _apply(f, TraitTarget(target), geom; threaded, calc_extent, kw...)
end

@inline _booltype(x::Bool) = x ? _True() : _False()
@inline _booltype(x::Union{_True,_False}) = x

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
@inline function _apply(f::F, target, ::Nothing, iterable; threaded, kw...) where F
    if threaded
        # `collect` first so we can use threads
        _apply(f, target, collect(iterable); threaded, kw...)
    else
        apply_to_iterable(x) = _apply(f, target, x; kw...)
        map(apply_to_iterable, iterable)
    end
end
# Rewrap all FeatureCollectionTrait feature collections as GI.FeatureCollection
# Maybe use threads to call _apply on componenet features
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
    applyreduce(f, op, target::Type{<:AbstractTrait}, obj; threaded)

Apply function `f` to all objects with the `target` trait,
and reduce the result with an `op` like `+`. 

The order and grouping of application of `op` is not guaranteed.

If `threaded==true` threads will be used over arrays and iterables, 
feature collections and nested geometries.
"""
@inline function applyreduce(
    f::F, op, target, geom; threaded=false, init=nothing
) where F
    threaded = _booltype(threaded)
    _applyreduce(f, op, TraitTarget(target), geom; threaded, init)
end

@inline _applyreduce(f::F, op, target, geom; threaded, init) where F =
    _applyreduce(f, op, target, GI.trait(geom), geom; threaded, init)
# Maybe use threads recucing over arrays
@inline function _applyreduce(f::F, op, target, ::Nothing, A::AbstractArray; threaded, init) where F
    applyreduce_array(i) = _applyreduce(f, op, target, A[i]; threaded=_False(), init)
    _mapreducetasks(applyreduce_array, op, eachindex(A), threaded; init)
end
# Try to applyreduce over iterables
@inline function _applyreduce(f::F, op, target, ::Nothing, iterable; threaded, init) where F
    applyreduce_iterable(i) = _applyreduce(f, op, target, x; threaded=_False(), init)
    if threaded # Try to `collect` and reduce over the vector with threads
        _applyreduce(f, op, target, collect(iterable); threaded, init)
    else
        # Try to `mapreduce` the iterable as-is
        mapreduce(applyreduce_iterable, op, iterable; init)
    end
end
# Maybe use threads reducing over features of feature collections
@inline function _applyreduce(f::F, op, target, ::GI.FeatureCollectionTrait, fc; threaded, init) where F
    applyreduce_fc(i) = _applyreduce(f, op, target, GI.getfeature(fc, i); threaded=_False(), init)
    _mapreducetasks(applyreduce_fc, op, 1:GI.nfeature(fc), threaded; init)
end
# Features just applyreduce to their geometry
@inline _applyreduce(f::F, op, target, ::GI.FeatureTrait, feature; threaded, init) where F =
    _applyreduce(f, op, target, GI.geometry(feature); threaded, init)
# Maybe use threads over components of nested geometries
@inline function _applyreduce(f::F, op, target, trait, geom; threaded, init) where F
    applyreduce_geom(i) = _applyreduce(f, op, target, GI.getgeom(geom, i); threaded=_False(), init)
    _mapreducetasks(applyreduce_geom, op, 1:GI.ngeom(geom), threaded; init)
end
# Don't thread over points it won't pay off
@inline function _applyreduce(
    f::F, op, target, trait::Union{GI.LinearRing,GI.LineString,GI.MultiPoint}, geom;
    threaded, init
) where F
    _applyreduce(f, op, target, GI.getgeom(geom); threaded=_False(), init)
end
# Apply f to the target
@inline function _applyreduce(f::F, op, ::TraitTarget{Target}, ::Trait, x; kw...) where {F,Target,Trait<:Target} 
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
    @eval _applyreduce(f::F, op, ::TraitTarget{<:$T}, trait::$T, x; kw...) where F = f(x)
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
        # The Boolean type parameters here indicate 3d-ness and measure coordinate presence respectively.
        return T{true,false}(child_geoms; crs, extent)
    else
        return T{false,false}(child_geoms; crs, extent)
    end
end
# So that GeometryBasics geoms rebuild as themselves
function rebuild(trait::GI.AbstractTrait, geom::BasicsGeoms, child_geoms; crs=nothing)
    GB.geointerface_geomtype(trait)(child_geoms)
end
function rebuild(trait::GI.AbstractTrait, geom::Union{GB.LineString,GB.MultiPoint}, child_geoms; crs=nothing)
    GB.geointerface_geomtype(trait)(GI.convert.(GB.Point, child_geoms))
end
function rebuild(trait::GI.PolygonTrait, geom::GB.Polygon, child_geoms; crs=nothing)
    Polygon(child_geoms[1], child_geoms[2:end])
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
