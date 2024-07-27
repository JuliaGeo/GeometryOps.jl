#=
# Convex hull

The _convex hull_ of a set of points is the smallest **convex** polygon that contains all the points.

GeometryOps.jl provides a number of methods for computing the convex hull of a set of points, usually
linked to other Julia packages.  

## Example

### Simple hull
```@example simple
import GeometryOps as GO, GeoInterface as GI
using CairoMakie # to plot

points = randn(GO.Point2f, 100)
f, a, p = plot(points; label = "Points")
hull_poly = GO.convex_hull(points)
lines!(a, hull_poly; label = "Convex hull", color = Makie.wong_colors()[2])
axislegend(a)
f
```

## Convex hull of the USA
```@example usa
import GeometryOps as GO, GeoInterface as GI
using CairoMakie # to plot
using NaturalEarth # for data

all_adm0 = naturalearth("admin_0_countries", 110)
usa = all_adm0.geometry[findfirst(==("USA"), all_adm0.ADM0_A3)]
f, a, p = lines(usa)
lines!(a, GO.convex_hull(usa); color = Makie.wong_colors()[2])
f
```

=#

"""
    convex_hull([method], geometries)

Compute the convex hull of the points in `geometries`.  
Returns a `GI.Polygon` representing the convex hull.

Note that all 

!!! warning
    This interface only computes the 2-dimensional convex hull!

    For higher dimensional hulls, use the relevant package (Qhull.jl or similar).
"""
function convex_hull end

"""
    MonotoneChainMethod()

This is an algorithm for the [`convex_hull`](@ref) function.

Uses [`DelaunayTriangulation.jl`](https://github.com/JuliaGeometry/DelaunayTriangulation.jl) to compute the convex hull.
This is a pure Julia algorithm which provides an optimal Delaunay triangulation.

See also [`convex_hull`](@ref)
"""
struct MonotoneChainMethod end

# GrahamScanMethod, etc. can be implemented in GO as well, if someone wants to.
# If we add an extension on Quickhull.jl, then that would be another algorithm.

convex_hull(geometries) = convex_hull(MonotoneChainMethod(), geometries)

# TODO: have this respect the CRS by pulling it out of `geometries`.
function convex_hull(::MonotoneChainMethod, geometries)
    # Extract all points as tuples.  We have to collect and allocate
    # here, because DelaunayTriangulation only accepts vectors of 
    # point-like geoms.

    # Cleanest would be to use the iterable from GO.flatten directly,
    # but that would require us to implement the convex hull algorithm
    # directly.

    # TODO: create a specialized method that extracts only the information
    # required, GeometryBasics points can be passed through directly.
    points = collect(flatten(tuples, GI.PointTrait, geometries))
    # Compute the convex hull using DelTri (shorthand for DelaunayTriangulation.jl).
    hull = DelaunayTriangulation.convex_hull(points)
    # Convert the result to a `GI.Polygon` and return it.
    # View would be more efficient here, but re-allocating
    # is cleaner.
    return GI.Polygon([GI.LinearRing(DelaunayTriangulation.get_points(hull)[(DelaunayTriangulation.get_vertices(hull))])])
end
