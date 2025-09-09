using Test
import GeometryOps as GO
using GeometryOpsCore
using GeometryOpsCore: MissingKeywordInAlgorithmException, WrongManifoldException

# These tests are more testing the call sites than they are testing the exceptions themselves.
# But this is also important to make sure that the whole pipeline works.
@testset "WrongManifoldException" begin
    @test_throws WrongManifoldException SingleManifoldAlgorithm{Planar}(Spherical())
    @test_throws "Planar" SingleManifoldAlgorithm{Planar}(Spherical())
    @test_throws "called with manifold Spherical" SingleManifoldAlgorithm{Planar}(Spherical())
end

@testset "MissingKeywordInAlgorithmException" begin
    alg = GO.GEOS(; tol = 1.0)
    @test_nowarn GO.enforce(alg, :tol, sum)
    @test_throws MissingKeywordInAlgorithmException GO.enforce(alg, :stat, sum)
    @test_throws "sum" GO.enforce(alg, :stat, sum)
    @test_throws "`stat`" GO.enforce(alg, :stat, sum)
end
