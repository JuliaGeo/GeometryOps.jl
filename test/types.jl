using Test
using GeometryOps
using GeometryOpsCore

using GeometryOps: Planar, Spherical, AutoManifold, rebuild, manifold, best_manifold, enforce, get
using GeometryOps: CLibraryPlanarAlgorithm, GEOS, TG, PROJ

@testset "CLibraryPlanarAlgorithm" begin
    # Test that it's a subtype of SingleManifoldAlgorithm{Planar}
    @test CLibraryPlanarAlgorithm <: GeometryOpsCore.SingleManifoldAlgorithm{Planar}

    # Test that it requires params::NamedTuple
    struct TestAlg <: CLibraryPlanarAlgorithm
        params::NamedTuple
    end

    # Test constructor with keyword arguments
    alg = TestAlg(; a=1, b=2)
    @test alg.params == (; a=1, b=2)

    # Test constructor with NamedTuple
    alg = TestAlg((; a=1, b=2))
    @test alg.params == (; a=1, b=2)

    # Test null constructor
    alg = TestAlg()
    @test alg.params == NamedTuple()

    # Test manifold methods
    @test manifold(alg) == Planar()
    @test best_manifold(alg, nothing) == Planar()

    # Test rebuild methods
    @test rebuild(alg, Planar()) == TestAlg(alg.params)
    @test rebuild(alg, AutoManifold()) == TestAlg(alg.params)
    @test_throws GeometryOpsCore.WrongManifoldException rebuild(alg, Spherical())
    @test rebuild(alg, (; c=3, d=4)) == TestAlg((; c=3, d=4))
    @test rebuild(alg; c=3, d=4) == TestAlg((; c=3, d=4))

    # Test get methods
    @test Base.get(alg, :a, 0) == 0
    @test Base.get(alg, :c, 0) == 0

    # Test enforce method
    @test_throws GeometryOpsCore.MissingKeywordInAlgorithmException enforce(alg, :a, "test")
    @test_throws GeometryOpsCore.MissingKeywordInAlgorithmException enforce(alg, :c, "test")
end

@testset "GEOS" begin
    # Test that it's a subtype of CLibraryPlanarAlgorithm
    @test GEOS() isa CLibraryPlanarAlgorithm

    # Test null constructor
    alg = GEOS()
    @test alg.params == NamedTuple()
    @test manifold(alg) == Planar()

    # Test constructor
    alg = GEOS(; a=1, b=2)
    @test alg.params == (; a=1, b=2)

    # Test manifold methods
    @test manifold(alg) == Planar()
    @test best_manifold(alg, nothing) == Planar()

    # Test rebuild methods
    @test rebuild(alg, Planar()) == GEOS(alg.params)
    @test rebuild(alg, AutoManifold()) == GEOS(alg.params)
    @test_throws GeometryOpsCore.WrongManifoldException rebuild(alg, Spherical())
    @test rebuild(alg, (; c=3, d=4)) == GEOS((; c=3, d=4))
    @test rebuild(alg; c=3, d=4) == GEOS((; c=3, d=4))

    # Test get methods
    @test Base.get(alg, :a, 0) == 1
    @test Base.get(alg, :c, 0) == 0

    # Test enforce method
    @test enforce(alg, :a, "test") == 1
    @test_throws GeometryOpsCore.MissingKeywordInAlgorithmException enforce(alg, :c, "test")
end

@testset "TG" begin
    # Test that it's a subtype of CLibraryPlanarAlgorithm
    @test TG <: CLibraryPlanarAlgorithm

    # Test null constructor
    alg = TG()
    @test alg.params == NamedTuple()
    @test manifold(alg) == Planar()

    # Test constructor
    alg = TG(; a=1, b=2)
    @test alg.params == (; a=1, b=2)

    # Test manifold methods
    @test manifold(alg) == Planar()
    @test best_manifold(alg, nothing) == Planar()

    # Test rebuild methods
    @test rebuild(alg, Planar()) == TG(alg.params)
    @test rebuild(alg, AutoManifold()) == TG(alg.params)
    @test_throws GeometryOpsCore.WrongManifoldException rebuild(alg, Spherical())
    @test rebuild(alg, (; c=3, d=4)) == TG((; c=3, d=4))
    @test rebuild(alg; c=3, d=4) == TG((; c=3, d=4))

    # Test get methods
    @test Base.get(alg, :a, 0) == 1
    @test Base.get(alg, :c, 0) == 0

    # Test enforce method
    @test enforce(alg, :a, "test") == 1
    @test_throws GeometryOpsCore.MissingKeywordInAlgorithmException enforce(alg, :c, "test")
end

@testset "PROJ" begin
    # Test that it's a subtype of Algorithm
    @test PROJ <: GeometryOpsCore.Algorithm

    # Test null constructor
    alg = PROJ()
    @test alg.manifold == Planar()
    @test alg.params == NamedTuple()

    # Test constructors
    alg1 = PROJ(; a=1, b=2)
    @test alg1.manifold == Planar()
    @test alg1.params == (; a=1, b=2)

    alg2 = PROJ(Spherical())
    @test alg2.manifold == Spherical()
    @test alg2.params == NamedTuple()

    # Test manifold method
    @test manifold(alg1) == Planar()
    @test manifold(alg2) == Spherical()

    # Test rebuild methods
    @test rebuild(alg1, Spherical()) == PROJ(Spherical(), alg1.params)
    @test rebuild(alg1, (; c=3, d=4)) == PROJ(Planar(), (; c=3, d=4))

    # Test get methods
    @test Base.get(alg1, :a, 0) == 1
    @test Base.get(alg1, :c, 0) == 0

    # Test enforce method
    @test enforce(alg1, :a, "test") == 1
    @test_throws GeometryOpsCore.MissingKeywordInAlgorithmException enforce(alg1, :c, "test")
end 