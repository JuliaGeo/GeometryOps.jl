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

"""
    apply(f, target::Type{<:AbstractTrait}, obj; kw...)

Reconstruct a geometry or feature using the function `f` on the `target` trait.

`f(target_geom) => x` where `x` also has the `target` trait, or an equivalent.

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

_apply(f, ::Type{Target}, geom; kw...)  where Target =
    _apply(f, Target, GI.trait(geom), geom; kw...)
function _apply(f, ::Type{Target}, ::Nothing, A::AbstractArray; threaded=false, kw...) where Target
    _maptasks(eachindex(A); threaded) do i
        _apply(f, Target, A[i]; kw...)
    end
end
# Try to _apply over iterables
_apply(f, ::Type{Target}, ::Nothing, iterable; kw...) where Target =
    map(x -> _apply(f, Target, x; kw...), iterable)
# Rewrap feature collections
function _apply(f, ::Type{Target}, ::GI.FeatureCollectionTrait, fc; 
    crs=GI.crs(fc), calc_extent=false, threaded=false
) where Target
    features = _maptasks(1:GI.nfeature(fc); threaded) do i
        feature = GI.getfeature(fc, i)
        _apply(f, Target, feature; crs, calc_extent)::GI.Feature
    end
    if calc_extent
        extent = reduce(features; init=GI.extent(first(features))) do acc, f
            Extents.union(acc, Extents.extent(f))
        end
        return GI.FeatureCollection(features; crs, extent)
    else
        return GI.FeatureCollection(features; crs)
    end
end
# Rewrap features
function _apply(f, ::Type{Target}, ::GI.FeatureTrait, feature; 
    crs=GI.crs(feature), calc_extent=false, threaded=false
) where Target
    properties = GI.properties(feature)
    geometry = _apply(f, Target, GI.geometry(feature); crs, calc_extent)
    if calc_extent
        extent = GI.extent(geometry)
        return GI.Feature(geometry; properties, crs, extent)
    else
        return GI.Feature(geometry; properties, crs)
    end
end
# Reconstruct nested geometries
function _apply(f, ::Type{Target}, trait, geom; 
    crs=GI.crs(geom), calc_extent=false, threaded=false
)::(GI.geointerface_geomtype(trait)) where Target
    # TODO handle zero length...
    geoms = _maptasks(1:GI.ngeom(geom); threaded) do i
        _apply(f, Target, GI.getgeom(geom, i); crs, calc_extent)
    end
    if calc_extent
        extent = GI.extent(first(geoms))
        for g in geoms
            extent = Extents.union(extent, GI.extent(g))
        end
        return rebuild(geom, geoms; crs, extent)
    else
        return rebuild(geom, geoms; crs)
    end
end
# Apply f to the target geometry
_apply(f, ::Type{Target}, ::Trait, geom; crs=GI.crs(geom), kw...) where {Target,Trait<:Target} = f(geom)
# Fail if we hit PointTrait without running `f`
_apply(f, ::Type{Target}, trait::GI.PointTrait, geom; crs=nothing, kw...) where Target =
    throw(ArgumentError("target $Target not found, but reached a `PointTrait` leaf"))
# Specific cases to avoid method ambiguity
_apply(f, ::Type{GI.PointTrait}, trait::GI.PointTrait, geom; kw...) = f(geom)
_apply(f, ::Type{GI.FeatureTrait}, ::GI.FeatureTrait, feature; kw...) = f(feature)
_apply(f, ::Type{GI.FeatureCollectionTrait}, ::GI.FeatureCollectionTrait, fc; kw...) = f(fc)

"""
    unwrap(target::Type{<:AbstractTrait}, obj)
    unwrap(f, target::Type{<:AbstractTrait}, obj)

Unwrap the geometry to vectors, down to the target trait.

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
    flatten(target::Type{<:GI.AbstractTrait}, geom)

Lazily flatten any geometry, feature or iterator of geometries or features
so that objects with the specified trait are returned by the iterator.
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
reconstruct(geom, components) = first(_reconstruct(geom, components))

_reconstruct(geom, components) = 
    _reconstruct(typeof(GI.trait(first(components))), geom, components, 1) 
_reconstruct(::Type{Target}, geom, components, iter) where Target = 
    _reconstruct(Target, GI.trait(geom), geom, components, iter)
# Try to reconstruct over iterables
function _reconstruct(::Type{Target}, ::Nothing, iterable, components, iter) where Target 
    vect = map(iterable) do x
        obj, iter = _reconstruct(Target, x, components, iter)
        obj
    end
    return vect, iter
end
# Reconstruct feature collections
function _reconstruct(::Type{Target}, ::GI.FeatureCollectionTrait, fc, components, iter) where Target
    features = map(GI.getfeature(fc)) do feature
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

By default geometries will be rebuilt as a GeoInterface.Wrappers 
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
        return reduce(vcat, map(fetch, tasks))
    else
        return map(f, taskrange)
    end
end
