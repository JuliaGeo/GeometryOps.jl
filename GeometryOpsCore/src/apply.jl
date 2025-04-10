# # `apply`

export apply

abstract type Applicator{F,T} end

for T in (:ApplyToGeom, :ApplyToArray, :ApplyToFeatures)
    @eval begin
        struct $T{F,T,O,K} <: Applicator{F,T}
            f::F
            target::T
            obj::O
            kw::K
        end
        $T(f, target; kw...) = $T(f, target, geom, kw)
    end
    # rebuild lets us swap out the function, such as with ThreadFunctors
    rebuild(a::Applicator, f) = $T(f, a.target, a.obj, a.kw) 
end

# Functor definitions
# _maptasks may run this level threaded if `threaded==true`
# but deeper `_apply` calls will not be threaded
# For an Array there is nothing to do but map `_apply` over all values
(a::ApplyToArray)(i::Int) = _apply(a.f, a.target, a.obj[i]; a.kw..., threaded=False())
# For a FeatureCollection or Geometry we need getfeature or getgeom calls
(a::ApplyToFeatures)(i::Int) = _apply(f, target, GI.getfeature(a.obj, i); a.kw..., threaded=False())
(a::ApplyToGeom)(i::Int) = _apply(a.f, a.target, GI.getgeom(a.obj, i); a.kw..., threaded=False())

#=

This file mainly defines the [`apply`](@ref) function.

In general, the idea behind the `apply` framework is to take 
as input any geometry, vector of geometries, or feature collection,
deconstruct it to the given trait target (any arbitrary GI.AbstractTrait 
or `TraitTarget` union thereof, like `PointTrait` or `PolygonTrait`) 
and perform some operation on it.  Then, the geometry or structure is rebuilt.

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

```@docs; collapse=true, canonical=false
apply
```


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

### Embedding

`extent` and `crs` can be embedded in all geometries, features, and
feature collections as part of `apply`. Geometries deeper than `Target`
will of course not have new `extent` or `crs` embedded.

- `calc_extent` signals to recalculate an `Extent` and embed it.
- `crs` will be embedded as-is

### Threading

Threading is used at the outermost level possible - over
an array, feature collection, or e.g. a MultiPolygonTrait where
each `PolygonTrait` sub-geometry may be calculated on a different thread.

Currently, threading defaults to `false` for all objects, but can be turned on
by passing the keyword argument `threaded=true` to `apply`.

Threading uses [StableTasks.jl](https://github.com/JuliaFolds2/StableTasks.jl) to provide
type-stable tasks (base Julia `Threads.@spawn` is not type stable).  This is completely cost-free
and saves some allocations when running multithreaded. 

The current strategy is to launch 2 tasks for each CPU thread, to provide load balancing.  We
assume Julia will manage these tasks efficiently, and we don't want to run too many tasks
since each task does have some overhead when it's created.  This may need revisiting in the future,
but it's a pretty easy heuristic to use.

## Implementation

Literate.jl source code is below.

***

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
    threaded = booltype(threaded)
    calc_extent = booltype(calc_extent)
    _apply(f, TraitTarget(target), geom; threaded, calc_extent, kw...)
end

# Call _apply again with the trait of `geom`
@inline _apply(f::F, target, geom; kw...)  where F =
    _apply(f, target, GI.trait(geom), geom; kw...)
# There is no trait and this is an AbstractArray - so just iterate over it calling _apply on the contents
@inline function _apply(f::F, target, ::Nothing, A::AbstractArray; threaded, kw...) where F
    applicator = ApplyToArray(f, target, A; kw...)
    _maptasks(applicator, eachindex(A), threaded)
end
# There is no trait and this is not an AbstractArray.
# Try to call _apply over it. We can't use threading
# as we don't know if we can can index into it. So just `map`.
@inline function _apply(f::F, target, ::Nothing, iterable::IterableType; threaded, kw...) where {F, IterableType}
    # Try the Tables.jl interface first
    if Tables.istable(iterable)
        _apply_table(f, target, iterable; threaded, kw...)
    else # this is probably some form of iterable...
        if threaded isa True
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
    result = Tables.materializer(iterable)(
        merge(
            NamedTuple{(geometry_column,), Base.Tuple{typeof(new_geometry)}}((new_geometry,)),
            NamedTuple(Iterators.map(_get_col_pair, new_names))
        )
    )
    # Finally, we ensure that metadata is propagated correctly.
    # This can only happen if the original table supports metadata reads,
    # and the result supports metadata writes.
    if DataAPI.metadatasupport(typeof(result)).write
        # Copy over all metadata from the original table to the new table, 
        # if the original table supports metadata reading.
        if DataAPI.metadatasupport(IterableType).read
            for (key, (value, style)) in DataAPI.metadata(iterable; style = true)
                # Default styles are not preserved on data transformation, so we must skip them!
                style == :default && continue
                # We assume that any other style is preserved.
                DataAPI.metadata!(result, key, value; style)
            end
        end
        # We don't usually care about the original table's metadata for GEOINTERFACE namespaced
        # keys, so we should set the crs and geometrycolumns metadata if they are present.
        # Ensure that `GEOINTERFACE:geometrycolumns` and `GEOINTERFACE:crs` are set!
        mdk = DataAPI.metadatakeys(result)
        # If the user has asked for geometry columns to persist, they would be here,
        # so we don't need to set them.
        if !("GEOINTERFACE:geometrycolumns" in mdk)
            # If the geometry columns are not already set, we need to set them.
            DataAPI.metadata!(result, "GEOINTERFACE:geometrycolumns", (geometry_column,); style = :default)
        end
        # Force reset CRS always, since you can pass `crs` to `apply`.
        new_crs = if haskey(kw, :crs)
            kw[:crs]
        else
            GI.crs(iterable) # this will automatically check `GEOINTERFACE:crs` unless the type has a specialized implementation.
        end

        DataAPI.metadata!(result, "GEOINTERFACE:crs", new_crs; style = :default)
    end

    return result
end

# Rewrap all FeatureCollectionTrait feature collections as GI.FeatureCollection
# Maybe use threads to call _apply on component features
@inline function _apply(f::F, target, ::GI.FeatureCollectionTrait, fc;
    crs=GI.crs(fc), calc_extent=False(), threaded
) where F

    # Run _apply on all `features` in the feature collection, possibly threaded
    applicator = ApplyToFeatures(f, target, fc; crs, calc_extent)
    features = _maptasks(applicator, 1:GI.nfeature(fc), threaded)
    if calc_extent isa True
        # Calculate the extent of the features
        extent = mapreduce(GI.extent, Extents.union, features)
        # Return a FeatureCollection with features, crs and calculated extent
        return GI.FeatureCollection(features; crs, extent)
    else
        # Return a FeatureCollection with features and crs
        return GI.FeatureCollection(features; crs)
    end
end
# Rewrap all FeatureTrait features as GI.Feature, keeping the properties
@inline function _apply(f::F, target, ::GI.FeatureTrait, feature;
    crs=GI.crs(feature), calc_extent=False(), threaded
) where F
    # Run _apply on the contained geometry
    geometry = _apply(f, target, GI.geometry(feature); crs, calc_extent, threaded)
    # Get the feature properties
    properties = GI.properties(feature)
    if calc_extent isa True
        # Calculate the extent of the geometry
        extent = GI.extent(geometry)
        # Return a new Feature with the new geometry and calculated extent, but the original properties and crs
        return GI.Feature(geometry; properties, crs, extent)
    else
        # Return a new Feature with the new geometry, but the original properties and crs
        return GI.Feature(geometry; properties, crs)
    end
end
# Reconstruct nested geometries,
# maybe using threads to call _apply on component geoms
@inline function _apply(f::F, target, trait, geom;
    crs=GI.crs(geom), calc_extent=False(), threaded
)::(GI.geointerface_geomtype(trait)) where F
    # Map `_apply` over all sub geometries of `geom`
    # to create a new vector of geometries
    # TODO handle zero length
    applicator = ApplyToGeom(f, target; crs, calc_extent)
    geoms = _maptasks(applicator, 1:GI.ngeom(geom), threaded)
    return _apply_inner(geom, geoms, crs, calc_extent)
end
@inline function _apply(f::F, target::TraitTarget{<:PointTrait}, trait::GI.PolygonTrait, geom;
    crs=GI.crs(geom), calc_extent=False(), threaded
)::(GI.geointerface_geomtype(trait)) where F
    # We need to force rebuilding a LinearRing not a LineString
    geoms = _maptasks(1:GI.ngeom(geom), threaded) do i
        lr = GI.getgeom(geom, i)
        points = map(GI.getgeom(lr)) do p
            _apply(f, target, p; crs, calc_extent, threaded=False())
        end
        _linearring(_apply_inner(lr, points, crs, calc_extent))
    end
    return _apply_inner(geom, geoms, crs, calc_extent)
end
function _apply_inner(geom, geoms, crs, calc_extent::True)
    # Calculate the extent of the sub geometries
    extent = mapreduce(GI.extent, Extents.union, geoms)
    # Return a new geometry of the same trait as `geom`,
    # holding the new `geoms` with `crs` and calculated extent
    return rebuild(geom, geoms; crs, extent)
end
function _apply_inner(geom, geoms, crs, calc_extent::False)
    # Return a new geometry of the same trait as `geom`, holding the new `geoms` with `crs`
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


### `_maptasks` - flexible, threaded `map`

using Base.Threads: nthreads, @threads, @spawn

#=
Here we used to use the compiler directive `@assume_effects :foldable` to force the compiler
to lookup through the closure. This alone makes e.g. `flip` 2.5x faster!

But it caused inference to fail, so we've removed it.  No effect on runtime so far as we can tell, 
at least in Julia 1.11.
=#
@inline function _maptasks(f::F, taskrange, threaded::False)::Vector where F
    map(f, taskrange)
end


# Threading utility, modified Mason Protters threading PSA
# run `f` over ntasks, where f receives an AbstractArray/range
# of linear indices
@inline function _maptasks(f::F, taskrange, threaded::True)::Vector where F
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
        StableTasks.@spawn begin
            # Where we map `f` over the chunk indices
            map(f, chunk)
        end
    end

    # Finally we join the results into a new vector
    return mapreduce(fetch, vcat, tasks)
end
@inline function _maptasks(a::Applicator{<:ThreadFunctors}, taskrange, threaded::True)::Vector
    ntasks = length(taskrange)
    chunk_size = max(1, ntasks ÷ (tf.tasks_per_thread * nthreads()))
    # partition the range into chunks
    task_chunks = Iterators.partition(taskrange, chunk_size)
    # Map over the chunks
    tasks = map(task_chunks, view(a.f.functors, eachindex(task_chunks))) do chunk, ft
        f = rebuild(a, ft)
        # Spawn a task to process this chunk
        StableTasks.@spawn begin
            # Where we map `f` over the chunk indices
            map(f, chunk)
        end
    end

    # Finally we join the results into a new vector
    return mapreduce(fetch, vcat, tasks)
end