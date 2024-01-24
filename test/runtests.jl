using GeometryOps
using Test

using GeometryOps.GeoInterface
using GeometryOps.GeometryBasics
using ArchGDAL
using LibGEOS
using Random, Distributions

const GI = GeoInterface
const AG = ArchGDAL
const LG = LibGEOS
const GO = GeometryOps

@testset "GeometryOps.jl" begin
    # @testset "Primitives" begin include("primitives.jl") end
    # # Methods
    # @testset "Area" begin include("methods/area.jl") end
    # @testset "Barycentric coordinate operations" begin include("methods/barycentric.jl") end
    # @testset "Bools" begin include("methods/bools.jl") end
    # @testset "Centroid" begin include("methods/centroid.jl") end
    # @testset "DE-9IM Geom Relations" begin include("methods/geom_relations.jl") end
    # @testset "Distance" begin include("methods/distance.jl") end
    # @testset "Equals" begin include("methods/equals.jl") end
    # Clipping
    @testset "Difference" begin include("methods/clipping/difference.jl") end
    @testset "Intersection" begin include("methods/clipping/intersection.jl") end
    @testset "Union" begin include("methods/clipping/union.jl") end
    @testset "Clipping Utils" begin include("methods/clipping/clipping_test_utils.jl") end
    # Transformations
    # @testset "Embed Extent" begin include("transformations/extent.jl") end
    # @testset "Reproject" begin include("transformations/reproject.jl") end
    # @testset "Flip" begin include("transformations/flip.jl") end
    # @testset "Simplify" begin include("transformations/simplify.jl") end
    # @testset "Transform" begin include("transformations/transform.jl") end
end
