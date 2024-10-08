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

Now, we'll also segmentize this using the geodesic method, which is more accurate for lat/lon coordinates.

```@example segmentize
using Proj # required to activate the `GeodesicSegments` method!
geodesic = GO.segmentize(GO.GeodesicSegments(max_distance = 1000), rectangle)
length(GI.getpoint(geodesic) |> collect)
```
This has a lot of points!  It's important to keep in mind that the `max_distance` is in meters, so this is a very fine-grained segmentation.

Now, let's see what they look like!  To make this fair, we'll use approximately the same number of points for both.

```@example segmentize
using CairoMakie
linear = GO.segmentize(rectangle; max_distance = 0.01)
geodesic = GO.segmentize(GO.GeodesicSegments(; max_distance = 1000), rectangle)
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
struct GeodesicSegments{T} <: SegmentizeMethod 
    geodesic::T# ::Proj.geod_geodesic
    max_distance::Float64
end

# Add an error hint for GeodesicSegments if Proj is not loaded!
function _geodesic_segments_error_hinter(io, exc, argtypes, kwargs)
    if isnothing(Base.get_extension(GeometryOps, :GeometryOpsProjExt)) && exc.f == GeodesicSegments
        print(io, "\n\nThe `Geodesic` method requires the Proj.jl package to be explicitly loaded.\n")
        print(io, "You can do this by simply typing ")
        printstyled(io, "using Proj"; color = :cyan, bold = true)
        println(io, " in your REPL, \nor otherwise loading Proj.jl via using or import.")
    end
end

# ## Implementation

"""
    segmentize([method = Planar()], geom; max_distance::Real, threaded)

Segmentize a geometry by adding extra vertices to the geometry so that no segment is longer than a given distance.  
This is useful for plotting geometries with a limited number of vertices, or for ensuring that a geometry is not too "coarse" for a given application.

## Arguments
- `method::Manifold = Planar()`: The method to use for segmentizing the geometry.  At the moment, only [`Planar`](@ref) (assumes a flat plane) and [`Geodesic`](@ref) (assumes geometry on the ellipsoidal Earth and uses Vincenty's formulae) are available.
- `geom`: The geometry to segmentize.  Must be a `LineString`, `LinearRing`, `Polygon`, `MultiPolygon`, or `GeometryCollection`, or some vector or table of those.
- `max_distance::Real`: The maximum distance between vertices in the geometry.  **Beware: for `Planar`, this is in the units of the geometry, but for `Geodesic` and `Spherical` it's in units of the radius of the sphere.**

Returns a geometry of similar type to the input geometry, but resampled.
"""
function segmentize(geom; max_distance, threaded::Union{Bool, BoolsAsTypes} = _False())
    return segmentize(Planar(), geom; max_distance, threaded = _booltype(threaded))
end

# allow three-arg method as well, just in case
segmentize(geom, max_distance::Real; threaded = _False()) = segmentize(Planar(), geom, max_distance; threaded)
segmentize(method::Manifold, geom, max_distance::Real; threaded = _False()) = segmentize(Planar(), geom; max_distance, threaded)

# generic implementation
function segmentize(method::Manifold, geom; max_distance, threaded::Union{Bool, BoolsAsTypes} = _False())
    @assert max_distance > 0 "`max_distance` should be positive and nonzero!  Found $(method.max_distance)."
    _segmentize_function(geom) = _segmentize(method, geom, GI.trait(geom); max_distance)
    return apply(_segmentize_function, TraitTarget(GI.LinearRingTrait(), GI.LineStringTrait()), geom; threaded)
end

function segmentize(method::SegmentizeMethod, geom; threaded::Union{Bool, BoolsAsTypes} = _False())
    @warn "`segmentize(method::$(typeof(method)), geom) is deprecated; use `segmentize($(method isa LinearSegments ? "Planar()" : "Geodesic()"), geom; max_distance, threaded) instead!"  maxlog=3
    @assert method.max_distance > 0 "`max_distance` should be positive and nonzero!  Found $(method.max_distance)."
    new_method = method isa LinearSegments ? Planar() : Geodesic()
    segmentize(new_method, geom; max_distance = method.max_distance, threaded)
end

_segmentize(method, geom) = _segmentize(method, geom, GI.trait(geom))
#= 
This is a method which performs the common functionality for both linear and geodesic algorithms, 
and calls out to the "kernel" function which we've defined per linesegment.
=#
function _segmentize(method::Union{Planar, Spherical}, geom, T::Union{GI.LineStringTrait, GI.LinearRingTrait}; max_distance)
    first_coord = GI.getpoint(geom, 1)
    x1, y1 = GI.x(first_coord), GI.y(first_coord)
    new_coords = NTuple{2, Float64}[]
    sizehint!(new_coords, GI.npoint(geom))
    push!(new_coords, (x1, y1))
    for coord in Iterators.drop(GI.getpoint(geom), 1)
        x2, y2 = GI.x(coord), GI.y(coord)
        _fill_linear_kernel!(method, new_coords, x1, y1, x2, y2; max_distance)
        x1, y1 = x2, y2
    end 
    return rebuild(geom, new_coords)
end

function _fill_linear_kernel!(::Planar, new_coords::Vector, x1, y1, x2, y2; max_distance)
    dx, dy = x2 - x1, y2 - y1
    distance = hypot(dx, dy) # this is a more stable way to compute the Euclidean distance
    if distance > max_distance
        n_segments = ceil(Int, distance / max_distance)
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
