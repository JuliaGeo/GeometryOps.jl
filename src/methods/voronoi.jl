#=
# Voronoi Tessellation

The [_Voronoi tessellation_](https://en.wikipedia.org/wiki/Voronoi_diagram) of a set of points is a partitioning of the plane into regions based on distance to points.
Each region contains all points closer to one generator point than to any other.

GeometryOps.jl provides a method for computing the Voronoi tessellation of a set of points,
using the [DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl) package.

## Example

### Simple tessellation
```@example simple
import GeometryOps as GO, GeoInterface as GI
using CairoMakie # to plot

points = tuple.(randn(20), randn(20))
polygons = GO.voronoi(points)
f, a, p = plot(polygons[1]; label = "Voronoi cell 1")
for (i, poly) in enumerate(polygons[2:end])
    plot!(a, poly; label = "Voronoi cell $(i+1)")
end
scatter!(a, points; color = :black, markersize = 10, label = "Generators")
axislegend(a)
f
```

## Implementation

This implementation mainly just preforms some assertion checks before passing the Arguments
to the DelaunayTriangulation package. It then unpacks the voronoi output into GeoInterface
polygons, whose point match the float type input by the user.
=#

struct __NoCRSProvided end

"""
    voronoi(geometries; clip = true)

Compute the Voronoi tessellation of the points in `geometries`.
Returns a vector of `GI.Polygon` objects representing the Voronoi cells,
in the same order as the input points.

## Arguments
- `geometries`: Any GeoInterface-compatible geometry or collection of geometries
  that can be decomposed into points

## Keyword Arguments  
- `clip`: Whether to clip the Voronoi cells to the convex hull of the input points (default: true)

!!! warning
    This interface only computes the 2-dimensional Voronoi tessellation!
    
!!! note
    The polygons are returned in the same order as the input points after flattening.
    Each polygon corresponds to the Voronoi cell of the point at the same index.
"""
function voronoi(geometries, ::Type{T} = Float64; kwargs...) where T
    return voronoi(Planar(), geometries, T; kwargs...)
end

function voronoi(::Planar, geometries, ::Type{T} = Float64; clip = true, clip_polygon = nothing, crs = __NoCRSProvided()) where T
    # Extract all points as tuples using GO.flatten
    # This handles any GeoInterface-compatible input
    points_iter = collect(flatten(tuples, GI.PointTrait, geometries))
    if crs isa __NoCRSProvided
        crs = GI.crs(geometries)
        # if isnothing(crs)
        #     crs = GI.crs(first(points_iter))
        # end
    end
    # if we have not figured it out yet, we can't do anything
    if crs isa __NoCRSProvided
        error("This code should be unreachable; please file an issue at https://github.com/JuliaGeometry/GeometryOps.jl/issues with the stacktrace and a reproducible example.")
    end
    
    # Handle edge case of too few points
    if length(points_iter) < 3
        throw(ArgumentError("Voronoi tessellation requires at least 3 points, got $(length(points_iter))"))
    end
    
    # Compute Delaunay triangulation
    tri = DelTri.triangulate(points_iter)
    
    # Compute Voronoi tessellation from the triangulation
    _clip_polygon = if isnothing(clip_polygon)
        nothing
    elseif GI.geomtrait(clip_polygon) isa GI.PolygonTrait
        @assert GI.nhole(clip_polygon) == 0
        (collect(flatten(tuples, GI.PointTrait, clip_polygon)), collect(1:GI.npoint(clip_polygon)))
    elseif clip_polygon isa Tuple{Vector{Tuple{T, T}}, Vector{Int}}
        clip_polygon
    elseif clip_polygon isa Tuple{NTuple{<:Any, Tuple{T, T}}, NTuple{<:Any, Int}}
        clip_polygon
    else
        error("Clip polygon must be a polygon or other recognizable form, see the docstring for `DelaunayTriangulation.voronoi` for the recognizable form.  Was neither, got $(typeof(clip_polygon))")
    end
    vorn = DelTri.voronoi(tri; clip = clip, clip_polygon = _clip_polygon)
    
    polygons = GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{Tuple{T, T}}, Nothing, Nothing}}, Nothing, typeof(crs)}[]
    sizehint!(polygons, DelTri.num_polygons(vorn))
    # Implementation below copied from Makie.jl
    # see https://github.com/MakieOrg/Makie.jl/blob/687c4466ce00154714297e36a7f610443c6ad5be/Makie/src/basic_recipes/voronoiplot.jl#L101-L110
    for i in DelTri.each_generator(vorn)
        !DelTri.has_polygon(vorn, i) && continue
        polygon_coords = DelTri.getxy.(DelTri.get_polygon_coordinates(vorn, i))
        push!(polygons, GI.Polygon([GI.LinearRing(polygon_coords)], crs = crs))
        # The code below gets the generator point, but we don't need it
        # gp = DelTri.getxy(DelTri.get_generator(vorn, i))
        # !isempty(polygon_coords) && push!(generators, gp)
    end
    
    return polygons
end
