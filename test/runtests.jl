using GeometryOps
using Test

using GeometryOps.GeoInterface
using GeometryOps.GeometryBasics
using ArchGDAL

const GI = GeoInterface
const AG = ArchGDAL

@testset "GeometryOps.jl" begin
    @testset "Primitives" begin include("primitives.jl") end
    @testset "Bools" begin include("methods/bools.jl") end
    @testset "Signed Area" begin include("methods/signed_area.jl") end
    @testset "Barycentric coordinate operations" begin include("methods/barycentric.jl") end
    @testset "Reproject" begin include("transformations/reproject.jl") end
    @testset "Flip" begin include("transformations/flip.jl") end
    @testset "Simplify" begin include("transformations/simplify.jl") end
end
