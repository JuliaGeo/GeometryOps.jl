using GeometryOps
using Test

using GeometryOps.GeoInterface
using ArchGDAL

const GI = GeoInterface
const AG = ArchGDAL

@testset "GeometryOps.jl" begin
    @testset "Signed Area" include("methods/signed_area.jl")
end
