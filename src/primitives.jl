const THREADED_KEYWORD = "- `threaded`: `true` or `false`. Whether to use multithreading. Defaults to `false`."
const CRS_KEYWORD = "- `crs`: The CRS to attach to geometries. Defaults to `nothing`."
const CALC_EXTENT_KEYWORD = "- `calc_extent`: `true` or `false`. Whether to calculate the extent. Defaults to `false`."

const APPLY_KEYWORDS = """
$THREADED_KEYWORD
$CRS_KEYWORD
$CALC_EXTENT_KEYWORD
"""

# # Primitive functions

# This file mainly defines the [`apply`](@ref) function.

#=
## What is `apply`?

`apply` apples some function to every geometry matching the `Target`
GeoInterface trait, in some abitrarily nested object made up of:
- `AbstractArray`s 
- Some arbitrary iterables may also work here
- `FeatureCollectionTrait` objects
- `FeatureTrait` objects
- `AbstractGeometryTrait` objects

It recussively calls `apply` through these nested
layers until it reaches the `Target`, where it applies `f`, and stops.

The outer recursive functions then progressively rebuild the object
using GeoInterface objects matchching the original traits.

If `PointTrait` is found  but it is not the `Target`, an error is thrown.
This likely means the object contains a different geometry trait to 
the target, such as `MultiPointTrait` when `LineStringTrait` was specified.

To handle this possibility it may be necessary to make `Target` a
`Union` of traits found at the same level of nesting, and define methods
of `f` to handle all cases.

Be careful making a union accross "levels" of nesting, e.g. 
`Union{FeatureTrait,PolygonTrait}`, as `_apply` will just never reach 
`PolygonTrait` when all the polgons are wrapped in a `FeatureTrait` object.

## Embedding:

`extent` and `crs` can be embededd in all geometries, features and
feature collections as part of `apply`. Geometries deeper than `Target`
will of course not hace new `extent` or `crs` embedded.

- `calc_extent` signals to recalculate an `Extent` and embed it. 
- `crs` will be embedded as-is

## Threading

Threading is used at the outermost level possible - over
a array, feature collection or e.g. a MultiPolygonTrait where
each `PolygonTrait` sub geometry may be calculated on a different thread.
=#

"""
    apply(f, target::Type{<:AbstractTrait}, obj; kw...)

Reconstruct a geometry, feature, feature collection or nested vectors of
either using the function `f` on the `target` trait.

`f(target_geom) => x` where `x` also has the `target` trait, or a trait that can 
be substituted. For example, swapping `PolgonTrait` to `MultiPointTrait` will fail
if the outer object has `MultiPolygonTrait`, but should work if it has `FeatureTrait`.

Objects "shallower" than the target trait are always completely rebuilt, like
a `Vector` of `FeatureCollectionTrait` of `FeatureTrait` when the target
has `PolygonTrait` and is held in the features. But "deeper" opjects may remain 
unchanged - such as points and linear rings if the tartet is the same `PolygonTrait`.

The result is an functionally similar geometry with values depending on `f`

$APPLY_KEYWORDS

# Example

Flipped point the order in any feature or geometry, or iterables of either:

```juia
import GeoInterface as GI
import GeometryOps as GO
geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

flipped_geom = GO.apply(GI.PointTrait, geom) do p
    (GI.y(p), GI.x(p))
end
"""
apply(f, ::Type{Target}, geom; kw...) where Target = _apply(f, Target, geom; kw...)

# Call _apply again with the trait of `geom`
_apply(f, ::Type{Target}, geom; kw...)  where Target =
    _apply(f, Target, GI.trait(geom), geom; kw...)
# There is no trait and this is an AbstractArray - so just iterate over it calling _apply on the contents
function _apply(f, ::Type{Target}, ::Nothing, A::AbstractArray; threaded=false, kw...) where Target
    # For an Array there is nothing else to do but map `_apply` over all values
    # _maptasks may run this level threaded if `threaded==true`, 
    # but deeper `_apply` called in the closure will not be threaded 
    _maptasks(eachindex(A); threaded) do i
        _apply(f, Target, A[i]; threaded=false, kw...)
    end
end
# Try to _apply over unknown iterables. We can't use threading on an 
# arbitrary iterable as we maybe can't index into it. So just `map`.
_apply(f, ::Type{Target}, ::Nothing, iterable; kw...) where Target =
    map(x -> _apply(f, Target, x; kw...), iterable)
# Rewrap all FeatureCollectionTrait feature collections as GI.FeatureCollection
function _apply(f, ::Type{Target}, ::GI.FeatureCollectionTrait, fc; 
    crs=GI.crs(fc), calc_extent=false, threaded=false
) where Target
    # Run _apply on all `features` in the feature collection, possibly threaded
    features = _maptasks(1:GI.nfeature(fc); threaded) do i
        feature = GI.getfeature(fc, i)
        _apply(f, Target, feature; crs, calc_extent, threaded=false)::GI.Feature
    end
    if calc_extent
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
function _apply(f, ::Type{Target}, ::GI.FeatureTrait, feature; 
    crs=GI.crs(feature), calc_extent=false, threaded=false
) where Target
    # Run _apply on the contained geometry 
    geometry = _apply(f, Target, GI.geometry(feature); crs, calc_extent, threaded)
    # Get the feature properties
    properties = GI.properties(feature)
    if calc_extent
        # Calculate the extent of the geometry
        extent = GI.extent(geometry)
        # Return a new Feature with the new geometry and calculated extent, but the oroginal properties and crs 
        return GI.Feature(geometry; properties, crs, extent)
    else
        # Return a new Feature with the new geometry, but the oroginal properties and crs 
        return GI.Feature(geometry; properties, crs)
    end
end
# Reconstruct nested geometries
function _apply(f, ::Type{Target}, trait, geom; 
    crs=GI.crs(geom), calc_extent=false, threaded=false
)::(GI.geointerface_geomtype(trait)) where Target
    # Map `_apply` over all sub geometries of `geom`
    # to create a new vector of geometries
    # TODO handle zero length
    geoms = _maptasks(1:GI.ngeom(geom); threaded) do i
        _apply(f, Target, GI.getgeom(geom, i); crs, calc_extent, threaded=false)
    end
    if calc_extent
        # Calculate the extent of the sub geometries
        extent = mapreduce(GI.extent, Extents.union, geoms)
        # Return a new geometry of the same trait as `geom`, 
        # holding tnew `geoms` with `crs` and calcualted extent
        return rebuild(geom, geoms; crs, extent)
    else
        # Return a new geometryof the same trait as `geom`, holding the new `geoms` with `crs`
        return rebuild(geom, geoms; crs)
    end
end
# Fail loudly if we hit PointTrait without running `f`
# (after PointTrait there is no further to dig with `_apply`)
_apply(f, ::Type{Target}, trait::GI.PointTrait, geom; crs=nothing, kw...) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Finally, these short methods are the main purpse of `apply`.
# The Trait is a subtype of the Target (or identical to it)
# So the Target is found. We apply `f` to it and return it to previous 
# _apply calls to be wrapped with the outer geometries/feature/featurecollection/array.
_apply(f, ::Type{Target}, ::Trait, geom; crs=GI.crs(geom), kw...) where {Target,Trait<:Target} = f(geom)
# Define some specific cases of this match to avoid method ambiguity
_apply(f, ::Type{GI.PointTrait}, trait::GI.PointTrait, geom; kw...) = f(geom)
_apply(f, ::Type{GI.FeatureTrait}, ::GI.FeatureTrait, feature; kw...) = f(feature)
_apply(f, ::Type{GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc; kw...) = f(fc)

"""
    unwrap(target::Type{<:AbstractTrait}, obj)
    unwrap(f, target::Type{<:AbstractTrait}, obj)

Unwrap the object newst to vectors, down to the target trait.

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
unwrap(f, target::Type, ::GI.FeatureTrait, feature) = unwrap(f, target, GI.geometry(feature))
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
function _maptasks(f, taskrange; threaded=false)
    if threaded
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
    else
        return map(f, taskrange)
    end
end
