using GeometryOps
using Test

using GeometryOps.GeoInterface
using ArchGDAL

const GI = GeoInterface
const AG = ArchGDAL

@testset "GeometryOps.jl" begin
    @testset "Primitives" begin include("primitives.jl") end
    @testset "Signed Area" begin include("methods/signed_area.jl") end
end
