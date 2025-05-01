#=
# `applyreduce`
=#

export applyreduce

#=
This file mainly defines the [`applyreduce`](@ref) function.  
    
This performs `apply`, but then reduces the result after flattening instead of rebuilding the geometry.


In general, the idea behind the `apply` framework is to take 
as input any geometry, vector of geometries, or feature collection,
deconstruct it to the given trait target (any arbitrary GI.AbstractTrait 
or `TraitTarget` union thereof, like `PointTrait` or `PolygonTrait`) 
and perform some operation on it.  

[`centroid`](@ref), [`area`](@ref) and [`distance`](@ref) have been implemented using the 
[`applyreduce`](@ref) framework.

```@docs
applyreduce
```


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
    threaded = booltype(threaded)
    _applyreduce(f, op, TraitTarget(target), geom; threaded, init)
end

@inline _applyreduce(f::F, op::O, target, geom; threaded, init) where {F, O} =
    _applyreduce(f, op, target, GI.trait(geom), geom; threaded, init)
# Maybe use threads reducing over arrays
@inline function _applyreduce(f::F, op::O, target, ::Nothing, A::AbstractArray; threaded, init) where {F, O}
    applyreduce_array(i) = _applyreduce(f, op, target, A[i]; threaded=False(), init)
    _mapreducetasks(applyreduce_array, op, eachindex(A), threaded; init)
end
# Try to applyreduce over iterables
@inline function _applyreduce(f::F, op::O, target, ::Nothing, iterable::IterableType; threaded, init) where {F, O, IterableType}
    if Tables.istable(iterable)
        _applyreduce_table(f, op, target, iterable; threaded, init)
    else
        applyreduce_iterable(i) = _applyreduce(f, op, target, i; threaded=False(), init)
        if threaded isa True # Try to `collect` and reduce over the vector with threads
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
    applyreduce_fc(i) = _applyreduce(f, op, target, GI.getfeature(fc, i); threaded=False(), init)
    _mapreducetasks(applyreduce_fc, op, 1:GI.nfeature(fc), threaded; init)
end
# Features just applyreduce to their geometry
@inline _applyreduce(f::F, op::O, target, ::GI.FeatureTrait, feature; threaded, init) where {F, O} =
    _applyreduce(f, op, target, GI.geometry(feature); threaded, init)
# Maybe use threads over components of nested geometries
@inline function _applyreduce(f::F, op::O, target, trait, geom; threaded, init) where {F, O}
    applyreduce_geom(i) = _applyreduce(f, op, target, GI.getgeom(geom, i); threaded=False(), init)
    _mapreducetasks(applyreduce_geom, op, 1:GI.ngeom(geom), threaded; init)
end
# Don't thread over points it won't pay off
@inline function _applyreduce(
    f::F, op::O, target, trait::Union{GI.LinearRing,GI.LineString,GI.MultiPoint}, geom;
    threaded, init
) where {F, O}
    _applyreduce(f, op, target, GI.getgeom(geom); threaded=False(), init)
end
# Apply f to the target
@inline function _applyreduce(f::F, op::O, ::TraitTarget{Target}, ::Trait, x; kw...) where {F,O,Target,Trait<:Target} 
    f(x)
end
@inline function _applyreduce(a::WithTrait{F}, op::O, ::TraitTarget{Target}, trait::Trait, x; kw...) where {F,O,Target,Trait<:Target} 
    a(trait, x; Base.structdiff(values(kw), NamedTuple{(:threaded, :init)})...)
end
# Fail if we hit PointTrait
# _applyreduce(f, op, target::TraitTarget{Target}, trait::PointTrait, geom; kw...) where Target = 
    # throw(ArgumentError("target $target not found"))
# Specific cases to avoid method ambiguity
for T in (
    GI.PointTrait, GI.LinearRing, GI.LineString, 
    GI.MultiPoint, GI.FeatureTrait, GI.FeatureCollectionTrait
)
    @eval begin
        _applyreduce(f::F, op::O, ::TraitTarget{<:$T}, trait::$T, x; kw...) where {F, O} = f(x)
        function _applyreduce(a::WithTrait{F}, op::O, ::TraitTarget{<:$T}, trait::$T, x; kw...) where {F, O}
            a(trait, x; Base.structdiff(values(kw), NamedTuple{(:threaded, :init)})...)
        end
    end
end

### `_mapreducetasks` - flexible, threaded mapreduce

import Base.Threads: nthreads, @threads, @spawn

# Threading utility, modified Mason Protters threading PSA
# run `f` over ntasks, where f receives an AbstractArray/range
# of linear indices
#
# WARNING: this will not work for mean/median - only ops
# where grouping is possible.  That's because the implementation operates
# in chunks, and not globally.
# 
# If you absolutely need a single chunk, then `threaded = false` will always decompose
# to straight `mapreduce` without grouping.
@inline function _mapreducetasks(f::F, op, taskrange, threaded::True; init) where F
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
            mapreduce(f, op, chunk; init)
        end
    end

    # Finally we join the results into a new vector
    return mapreduce(fetch, op, tasks; init)
end

function _mapreducetasks(f::F, op, taskrange, threaded::False; init) where F
    mapreduce(f, op, taskrange; init)
end
