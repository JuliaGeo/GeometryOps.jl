using Test
using SpatialTreeInterface

@testset "SpatialTreeInterface" begin
    @testset "LoopStateMachine" begin include("LoopStateMachine.jl") end
    @testset "Trees and traversals" begin include("SpatialTreeInterface.jl") end
end
