using Test
import GeometryOpsCore
import GeometryOpsCore: Prepared, AbstractPreparation, AbstractPreparationTrait,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike,
    preptrait, getprep, Planar, Spherical, manifold
import GeoInterface as GI
import GeoInterface: Extents

# Dummy preparations for interface tests
struct _DummyEdgeTree <: AbstractPreparation
    payload::Int
end
GeometryOpsCore.preptrait(::_DummyEdgeTree) = SpatialEdgeIndexLike()

struct _DummyChildTree <: AbstractPreparation end
GeometryOpsCore.preptrait(::_DummyChildTree) = SpatialIndexLike()

@testset "Preparation retrieval" begin
    ring = GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])
    ext = GI.extent(ring; fallback = true)
    p = Prepared(ring, Planar(), (_DummyEdgeTree(42), _DummyChildTree()), ext)

    @test manifold(p) === Planar()
    @test parent(p) === ring

    # hit: first matching preparation, concretely typed
    et = get(p, SpatialEdgeIndexLike())
    @test et isa _DummyEdgeTree && et.payload == 42
    @test get(p, SpatialIndexLike()) isa _DummyChildTree

    # miss: nothing
    @test get(p, PointInAreaLike()) === nothing

    # plain geometries always miss — uniform call sites
    @test get(ring, SpatialEdgeIndexLike()) === nothing

    # manifold-checked retrieval: mismatch is a miss, not an error
    @test getprep(Planar(), p, SpatialEdgeIndexLike()) === et
    @test getprep(Spherical(), p, SpatialEdgeIndexLike()) === nothing
    @test getprep(Planar(), ring, SpatialEdgeIndexLike()) === nothing
end
