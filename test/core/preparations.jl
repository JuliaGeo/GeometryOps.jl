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

    # first-match-wins tie-break: two preparations with the same capability → the first
    p2 = Prepared(ring, Planar(), (_DummyEdgeTree(1), _DummyEdgeTree(2)), ext)
    @test get(p2, SpatialEdgeIndexLike()).payload == 1

    # constructor validation: a non-geometry is rejected
    @test_throws ArgumentError Prepared(42, Planar(), (), ext)
end

@testset "GeoInterface forwarding" begin
    ring = GI.LinearRing([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 0.0)])
    rext = GI.extent(ring; fallback = true)
    pring = Prepared(ring, Planar(), (_DummyEdgeTree(1),), rext)

    @test GI.isgeometry(pring)
    @test GI.trait(pring) isa GI.LinearRingTrait
    @test GI.npoint(pring) == GI.npoint(ring)
    @test collect(GI.getpoint(pring)) == collect(GI.getpoint(ring))
    @test GI.coordinates(pring) == GI.coordinates(ring)
    @test GI.is3d(pring) == false
    @test GI.extent(pring) == rext          # served from the field
    @test Extents.extent(pring) == rext

    # a Prepared polygon whose ring is itself Prepared: GI hands back the prepared child
    poly = GI.Polygon([pring])
    pext = GI.extent(ring; fallback = true)
    ppoly = Prepared(poly, Planar(), (), pext)
    @test GI.trait(ppoly) isa GI.PolygonTrait
    @test GI.nring(ppoly) == 1
    @test GI.getring(ppoly, 1) === pring
    @test GI.getexterior(ppoly) === pring
    @test get(GI.getring(ppoly, 1), SpatialEdgeIndexLike()) isa _DummyEdgeTree
end
