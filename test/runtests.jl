using GeometryOps
using Test

using GeometryOps.GeoInterface
using GeometryOps.GeometryBasics
using GeoInterface.Extents: Extents
using ArchGDAL
using LibGEOS
using Random, Distributions
using Proj

const GI = GeoInterface
const AG = ArchGDAL
const LG = LibGEOS
const GO = GeometryOps

@testset "GeometryOps.jl" begin
    @testset "Primitives" begin include("primitives.jl") end
    # # # Methods
    @testset "Angles" begin include("methods/angles.jl") end
    @testset "Area" begin include("methods/area.jl") end
    @testset "Barycentric coordinate operations" begin include("methods/barycentric.jl") end
    @testset "Orientation" begin include("methods/orientation.jl") end
    @testset "Centroid" begin include("methods/centroid.jl") end
    @testset "DE-9IM Geom Relations" begin include("methods/geom_relations.jl") end
    @testset "Distance" begin include("methods/distance.jl") end
    @testset "Equals" begin include("methods/equals.jl") end
    # # # Clipping
    @testset "Coverage" begin include("methods/clipping/coverage.jl") end
    @testset "Cut" begin include("methods/clipping/cut.jl") end
    @testset "Polygon Clipping" begin include("methods/clipping/polygon_clipping.jl") end
    # # Transformations
    @testset "Embed Extent" begin include("transformations/extent.jl") end
    @testset "Reproject" begin include("transformations/reproject.jl") end
    @testset "Flip" begin include("transformations/flip.jl") end
    @testset "Simplify" begin include("transformations/simplify.jl") end
    @testset "Segmentize" begin include("transformations/segmentize.jl") end
    @testset "Transform" begin include("transformations/transform.jl") end
    @testset "Geometry correction" begin 
        include("transformations/correction/geometry_correction.jl")
        include("transformations/correction/closed_ring.jl") 
        include("transformations/correction/minimal_multipolygon.jl")
    end
end
