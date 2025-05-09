using GeometryOps
using Test
using SafeTestsets

include("helpers.jl")

@testset "Core" begin
    @safetestset "Algorithm" begin include("core/algorithm.jl") end
    @safetestset "Manifold" begin include("core/manifold.jl") end
    @safetestset "Applicators" begin include("core/applicators.jl") end
end

@safetestset "Types" begin include("types.jl") end
@safetestset "Primitives" begin include("primitives.jl") end

# Utils
@safetestset "Utils" begin include("utils/utils.jl") end
@safetestset "LoopStateMachine" begin include("utils/LoopStateMachine.jl") end
@safetestset "SpatialTreeInterface" begin include("utils/SpatialTreeInterface.jl") end

# Methods
@safetestset "Angles" begin include("methods/angles.jl") end
@safetestset "Area" begin include("methods/area.jl") end
@safetestset "Barycentric coordinate operations" begin include("methods/barycentric.jl") end
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
@safetestset "Force Dimensions" begin include("transformations/forcedims.jl") end
@safetestset "Geometry Correction" begin include("transformations/correction/geometry_correction.jl") end
@safetestset "Closed Rings" begin include("transformations/correction/closed_ring.jl")  end
@safetestset "Intersecting Polygons" begin include("transformations/correction/intersecting_polygons.jl") end
# Extensions
@safetestset "FlexiJoins" begin include("extensions/flexijoins.jl") end
@safetestset "LibGEOS" begin include("extensions/libgeos.jl") end
@safetestset "TGGeometry" begin include("extensions/tggeometry.jl") end
@safetestset "DataFrames" begin include("extensions/dataframes.jl") end