using Test, UnitSphericalGeometryCore

@testset "UnitSphericalGeometryCore" begin
    @testset "Points, transforms, slerp and caps" begin include("unitspherical.jl") end
    @testset "RobustCrossProduct" begin include("robustcrossproduct.jl") end
end
