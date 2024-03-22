# # Segmentize

export segmentize
export LinearSegments, GeodesicSegments

#=
This function "segmentizes" or "densifies" a geometry by adding 
extra vertices to the geometry so that no segment is longer than 
a given distance.  This is useful for plotting geometries with a 
limited number of vertices, or for ensuring that a geometry is not 
too "coarse" for a given application.

!!! info
    We plan to add interpolated segmentization from DataInterpolations.jl in the future, 
    which will be available to any vector of point-like objects. 

    For now, this function only works on 2D geometries.  We will also support 3D geometries, as well as measure interpolation, in the future.

## Examples

```@example segmentize
import GeometryOps as GO, GeoInterface as GI
rectangle = GI.Wrappers.Polygon([[(0.0, 50.0), (7.071, 57.07), (0, 64.14), (-7.07, 57.07), (0.0, 50.0)]])
linear = GO.segmentize(rectangle; max_distance = 5)
collect(GI.getpoint(linear))
```
You can see that this geometry was segmentized correctly, and now has 8 vertices where it previously had only 4.

Now, we'll also segmentize this using the geodetic method, which is more accurate for lat/lon coordinates.

```@example segmentize
geodesic = GO.segmentize(GO.GeodesicSegments(max_distance = 1000), rectangle)
length(GI.getpoint(geodesic) |> collect)
```
This has a lot of points!  It's important to keep in mind that the `max_distance` is in meters, so this is a very fine-grained segmentation.

Now, let's see what they look like!  To make this fair, we'll use approximately the same number of points for both.

```@example segmentize
using CairoMakie
linear = GO.segmentize(rectangle; max_distance = 0.01)
geodesic = GO.segmentize(GO.GeodesicSegments(; max_distance = 900), rectangle)
f, a, p = poly(collect(GI.getpoint(linear)); label = "Linear", axis = (; aspect = DataAspect()))
p2 = poly!(collect(GI.getpoint(geodesic)); label = "Geodesic")
axislegend(a; position = :lt)
f
```

There are two methods available for segmentizing geometries at the moment: 

```@docs
LinearSegments
GeodesicSegments
```

## Benchmark

We benchmark our method against LibGEOS's `GEOSDensify` method, which is a similar method for densifying geometries.

```@example benchmark
using BenchmarkTools: BenchmarkGroup
using Chairmarks: @be
using Main: plot_trials
using CairoMakie

import GeometryOps as GO, GeoInterface as GI, LibGEOS as LG

segmentize_suite = BenchmarkGroup(["title:Segmentize", "subtitle:Segmentize a rectangle"])

rectangle = GI.Wrappers.Polygon([[(0.0, 50.0), (7.071, 57.07), (0.0, 64.14), (-7.07, 57.07), (0.0, 50.0)]])
lg_rectangle = GI.convert(LG, rectangle)
```

```@example benchmark
# These are initial distances, which yield similar numbers of points 
# in the final geometry.
init_lin = 0.01
init_geo = 900

# LibGEOS.jl doesn't offer this function, so we just wrap it ourselves!
function densify(obj::LG.Geometry, tol::Real, context::LG.GEOSContext = LG.get_context(obj))
    result = LG.GEOSDensify_r(context, obj, tol)
    if result == C_NULL
        error("LibGEOS: Error in GEOSDensify")
    end
    LG.geomFromGEOS(result, context)
end
# now, we get to the actual benchmarking:
for scalefactor in exp10.(LinRange(log10(0.1), log10(10), 5))
    lin_dist = init_lin * scalefactor
    geo_dist = init_geo * scalefactor

    npoints_linear = GI.npoint(GO.segmentize(rectangle; max_distance = lin_dist))
    npoints_geodesic = GO.segmentize(GO.GeodesicSegments(; max_distance = geo_dist), rectangle) |> GI.npoint
    npoints_libgeos = GI.npoint(densify(lg_rectangle, lin_dist))
    
    segmentize_suite["Linear"][npoints_linear] = @be GO.segmentize(GO.LinearSegments(; max_distance = $lin_dist), $rectangle) seconds=1
    segmentize_suite["Geodesic"][npoints_geodesic] = @be GO.segmentize(GO.GeodesicSegments(; max_distance = $geo_dist), $rectangle) seconds=1
    segmentize_suite["LibGEOS"][npoints_libgeos] = @be densify($lg_rectangle, $lin_dist) seconds=1
    
end

plot_trials(segmentize_suite)
```

=#

abstract type SegmentizeMethod end
"""
    LinearSegments(; max_distance::Real)

A method for segmentizing geometries by adding extra vertices to the geometry so that no segment is longer than a given distance.

Here, `max_distance` is a purely nondimensional quantity and will apply in the input space.   This is to say, that if the polygon is
provided in lat/lon coordinates then the `max_distance` will be in degrees of arc.  If the polygon is provided in meters, then the 
`max_distance` will be in meters.
"""
Base.@kwdef struct LinearSegments <: SegmentizeMethod 
    max_distance::Float64
end
"""
    GeodesicSegments(; max_distance::Real, equatorial_radius::Real=6378137, flattening::Real=1/298.257223563)

A method for segmentizing geometries by adding extra vertices to the geometry so that no segment is longer than a given distance.  
This method calculates the distance between points on the geodesic, and assumes input in lat/long coordinates.

!!! warning
    Any input geometries must be in lon/lat coordinates!  If not, the method may fail or error.

## Arguments
- `max_distance::Real`: The maximum distance, **in meters**, between vertices in the geometry.
- `equatorial_radius::Real=6378137`: The equatorial radius of the Earth, in meters.  Passed to `Proj.geod_geodesic`.
- `flattening::Real=1/298.257223563`: The flattening of the Earth, which is the ratio of the difference between the equatorial and polar radii to the equatorial radius.  Passed to `Proj.geod_geodesic`.

One can also omit the `equatorial_radius` and `flattening` keyword arguments, and pass a `geodesic` object directly to the eponymous keyword.

This method uses the Proj/GeographicLib API for geodesic calculations.
"""
struct GeodesicSegments <: SegmentizeMethod 
    geodesic# ::Proj.geod_geodesic
    max_distance::Float64
end

# ## Implementation

"""
    segmentize([method = LinearSegments()], geom; max_distance::Real, threaded)

Segmentize a geometry by adding extra vertices to the geometry so that no segment is longer than a given distance.  
This is useful for plotting geometries with a limited number of vertices, or for ensuring that a geometry is not too "coarse" for a given application.

## Arguments
- `method::SegmentizeMethod = LinearSegments()`: The method to use for segmentizing the geometry.  At the moment, only [`LinearSegments`](@ref) and [`GeodesicSegments`](@ref) are available.
- `geom`: The geometry to segmentize.  Must be a `LineString`, `LinearRing`, or greater in complexity.
- `max_distance::Real`: The maximum distance, **in the input space**, between vertices in the geometry.  Only used if you don't explicitly pass a `method`.

Returns a geometry of similar type to the input geometry, but resampled.
"""
function segmentize(geom; max_distance, threaded::Union{Bool, BoolsAsTypes} = _False())
    return segmentize(LinearSegments(; max_distance), geom; threaded)
end
function segmentize(method::SegmentizeMethod, geom; threaded::Union{Bool, BoolsAsTypes} = _False())
    @assert method.max_distance > 0 "`max_distance` should be positive and nonzero!  Found $(method.max_distance)."
    segmentize_function = Base.Fix1(_segmentize, method)
    return apply(segmentize_function, Union{GI.LinearRingTrait, GI.LineStringTrait}, geom; threaded)
end

_segmentize(method, geom) = _segmentize(method, geom, GI.trait(geom))
#= 
This is a method which performs the common functionality for both linear and geodesic algorithms, 
and calls out to the "kernel" function which we've defined per linesegment.
=#
function _segmentize(method::Union{LinearSegments, GeodesicSegments}, geom, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    first_coord = GI.getpoint(geom, 1)
    x1, y1 = GI.x(first_coord), GI.y(first_coord)
    new_coords = NTuple{2, Float64}[]
    sizehint!(new_coords, GI.npoint(geom))
    push!(new_coords, (x1, y1))
    for coord in Iterators.drop(GI.getpoint(geom), 1)
        x2, y2 = GI.x(coord), GI.y(coord)
        _fill_linear_kernel!(method, new_coords, x1, y1, x2, y2)
        x1, y1 = x2, y2
    end 
    return rebuild(geom, new_coords)
end

function _fill_linear_kernel!(method::LinearSegments, new_coords::Vector, x1, y1, x2, y2)
    dx, dy = x2 - x1, y2 - y1
    distance = hypot(dx, dy) # this is a more stable way to compute the Euclidean distance
    if distance > method.max_distance
        n_segments = ceil(Int, distance / method.max_distance)
        for i in 1:(n_segments - 1)
            t = i / n_segments
            push!(new_coords, (x1 + t * dx, y1 + t * dy))
        end
    end
    # End the line with the original coordinate,
    # to avoid any multiplication errors.
    push!(new_coords, (x2, y2))
    return nothing
end
#=

!!! note
    The `_fill_linear_kernel` definition for `GeodesicSegments` is in the `GeometryOpsProjExt` extension module, in the `segmentize.jl` file.

=#
