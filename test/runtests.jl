using GeometryOps
using Test

using GeometryOps.GeoInterface
using ArchGDAL

const GI = GeoInterface
const AG = ArchGDAL

@testset "GeometryOps.jl" begin
    @testset "Primitives" begin include("primitives.jl") end
    @testset "Signed Area" begin include("methods/signed_area.jl") end
    @testset "reproject" begin include("transformations/reproject.jl") end
    @testset "flip" begin include("transformations/flip.jl") end
    @testset "simplify" begin include("transformations/simplify.jl") end
end
