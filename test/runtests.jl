using GeometryOps
using Test
using SafeTestsets

include("helpers.jl")

@safetestset "Primitives" begin include("primitives.jl") end
@safetestset "Lazy Wrappers" begin include("lazy_wrappers.jl") end
# Methods
@safetestset "Angles" begin include("methods/angles.jl") end
@safetestset "Area" begin include("methods/area.jl") end
# @safetestset "Barycentric coordinate operations" begin include("methods/barycentric.jl") end
@safetestset "Orientation" begin include("methods/orientation.jl") end
@safetestset "Centroid" begin include("methods/centroid.jl") end
@safetestset "Convex Hull" begin include("methods/convex_hull.jl") end
@safetestset "DE-9IM Geom Relations" begin include("methods/geom_relations.jl") end
@safetestset "Distance" begin include("methods/distance.jl") end
@safetestset "Equals" begin include("methods/equals.jl") end
# Clipping
@safetestset "Coverage" begin include("methods/clipping/coverage.jl") end
@safetestset "Cut" begin include("methods/clipping/cut.jl") end
@safetestset "Intersection Point" begin include("methods/clipping/intersection_points.jl") end
@safetestset "Polygon Clipping" begin include("methods/clipping/polygon_clipping.jl") end
# Transformations
@safetestset "Embed Extent" begin include("transformations/extent.jl") end
@safetestset "Reproject" begin include("transformations/reproject.jl") end
@safetestset "Flip" begin include("transformations/flip.jl") end
@safetestset "Simplify" begin include("transformations/simplify.jl") end
@safetestset "Segmentize" begin include("transformations/segmentize.jl") end
@safetestset "Transform" begin include("transformations/transform.jl") end
@safetestset "Geometry Correction" begin include("transformations/correction/geometry_correction.jl") end
@safetestset "Closed Rings" begin include("transformations/correction/closed_ring.jl")  end
@safetestset "Intersecting Polygons" begin include("transformations/correction/intersecting_polygons.jl") end
# Extensions
@safetestset "FlexiJoins" begin include("extensions/flexijoins.jl") end
@safetestset "LibGEOS" begin include("extensions/libgeos.jl") end
