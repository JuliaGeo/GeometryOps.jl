using GeometryOpsCore

@testset "Constructing NoAlgorithm" begin
    @test NoAlgorithm() isa NoAlgorithm
end

@testset "Constructing AutoAlgorithm" begin
    @test AutoAlgorithm() isa AutoAlgorithm
    @test AutoAlgorithm(; x = 1) == AutoAlgorithm(AutoManifold(), pairs((; x = 1)))
end

@testset "SingleManifoldAlgorithm" begin
    struct TestSMAlgorithm <: SingleManifoldAlgorithm{Planar}
    end

    @test TestSMAlgorithm() isa TestSMAlgorithm

    @test_throws GeometryOpsCore.WrongManifoldException TestSMAlgorithm(Geodesic())
end

@testset "ManifoldIndependentAlgorithm" begin
    struct TestMIDAlgorithm{M} <: ManifoldIndependentAlgorithm{M}
        m::M
    end

    @test TestMIDAlgorithm(Planar()) isa TestMIDAlgorithm
    @test TestMIDAlgorithm(Spherical()) isa TestMIDAlgorithm
    @test TestMIDAlgorithm(Geodesic()) isa TestMIDAlgorithm
end
