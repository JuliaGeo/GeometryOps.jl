using Test
import GeometryOps as GO
import GeometryOps: Prepared, prepare, getprep,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike
import GeometryOpsCore
import GeoInterface as GI
import GeoInterface: Extents

@testset "prepare is one generic function" begin
    # the algorithm method and the (future) geometry method share the Core binding
    @test GO.prepare === GeometryOpsCore.prepare
    poly = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])])
    prep = GO.prepare(GO.RelateNG(), poly)
    @test prep isa GO.PreparedRelate
end

# shared fixtures for the remaining testsets
const _PG_RING_OUTER = GI.LinearRing([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)])
const _PG_RING_HOLE  = GI.LinearRing([(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0), (4.0, 4.0)])
const _PG_POLY  = GI.Polygon([_PG_RING_OUTER, _PG_RING_HOLE])
const _PG_POLY2 = GI.Polygon([GI.LinearRing([(20.0, 0.0), (25.0, 0.0), (25.0, 5.0), (20.0, 0.0)])])
const _PG_MP    = GI.MultiPolygon([_PG_POLY, _PG_POLY2])
# mid-latitude spherical quad (lon/lat degrees)
const _PG_SPH_POLY = GI.Polygon([GI.LinearRing([(10.0, 40.0), (20.0, 40.0), (20.0, 50.0), (10.0, 50.0), (10.0, 40.0)])])

@testset "prepare: recursive wrapping, no specs" begin
    p = prepare(_PG_POLY)
    @test p isa Prepared
    @test GI.trait(p) isa GI.PolygonTrait
    @test GeometryOpsCore.manifold(p) === GO.Planar()
    @test p.preps === ()

    # children are themselves Prepared, with per-node extents
    r1 = GI.getring(p, 1)
    @test r1 isa Prepared && GI.trait(r1) isa GI.LinearRingTrait
    @test GI.extent(r1) == GI.extent(_PG_RING_OUTER; fallback = true)
    @test GI.extent(p) == GI.extent(_PG_POLY; fallback = true)

    # geometry content is unchanged
    @test GI.coordinates(p) == GI.coordinates(_PG_POLY)

    # points pass through unwrapped
    pt = GI.Point(1.0, 2.0)
    @test prepare(pt) === pt

    # multipolygon: every level wrapped
    mp = prepare(_PG_MP)
    @test mp isa Prepared && GI.getgeom(mp, 1) isa Prepared
    @test GI.getring(GI.getgeom(mp, 1), 1) isa Prepared

    # spherical: extents are 3D unit-sphere boxes
    ps = prepare(_PG_SPH_POLY; manifold = GO.Spherical())
    @test GI.extent(ps) isa Extents.Extent{(:X, :Y, :Z)}
    @test GI.extent(GI.getring(ps, 1)) isa Extents.Extent{(:X, :Y, :Z)}
end

@testset "RingEdgeIndex" begin
    p = prepare(_PG_POLY; preps = (GO.RingEdgeIndex(),))
    @test p.preps === ()                       # consumed at ring level, not polygon level
    ring = GI.getring(p, 1)
    prep = get(ring, SpatialEdgeIndexLike())
    @test prep isa GO.SpatialEdgeIndex
    tree = prep.tree
    @test tree isa GO.NaturalIndexing.NaturalIndex

    # the tree indexes segments in ring order: query a box covering only segment 1
    # (outer ring segment 1 runs (0,0)->(10,0))
    hits = GO.SpatialTreeInterface.query(tree, Extents.Extent(X = (4.0, 5.0), Y = (-0.1, 0.1)))
    @test hits == [1]

    # spherical: builds 3D arc-extent trees without error
    ps = prepare(_PG_SPH_POLY; manifold = GO.Spherical(), preps = (GO.RingEdgeIndex(),))
    sprep = get(GI.getring(ps, 1), SpatialEdgeIndexLike())
    @test sprep isa GO.SpatialEdgeIndex
    @test GO.SpatialTreeInterface.node_extent(sprep.tree) isa Extents.Extent{(:X, :Y, :Z)}

    # manifold-checked retrieval
    @test getprep(GO.Planar(), GI.getring(p, 1), SpatialEdgeIndexLike()) === prep
    @test getprep(GO.Spherical(), GI.getring(p, 1), SpatialEdgeIndexLike()) === nothing
end
