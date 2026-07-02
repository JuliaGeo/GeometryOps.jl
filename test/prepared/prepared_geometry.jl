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

    # every edge of a square ring is individually isolatable: a thin box centered on
    # each edge midpoint hits exactly that segment
    sq = prepare(GI.LinearRing([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]);
        preps = (GO.RingEdgeIndex(),))
    sqtree = get(sq, SpatialEdgeIndexLike()).tree
    @test GO.SpatialTreeInterface.query(sqtree, Extents.Extent(X = (4.9, 5.1), Y = (-0.1, 0.1))) == [1]  # y ≈ 0 edge
    @test GO.SpatialTreeInterface.query(sqtree, Extents.Extent(X = (9.9, 10.1), Y = (4.9, 5.1))) == [2]  # x ≈ 10 edge
    @test GO.SpatialTreeInterface.query(sqtree, Extents.Extent(X = (4.9, 5.1), Y = (9.9, 10.1))) == [3]  # y ≈ 10 edge
    @test GO.SpatialTreeInterface.query(sqtree, Extents.Extent(X = (-0.1, 0.1), Y = (4.9, 5.1))) == [4]  # x ≈ 0 edge

    # spherical: the tree's stored arc extents round-trip against the kernel's
    # `_segment_extent` — segment 2's own box hits it; a provably disjoint segment's box
    # cannot hit segment 1. (Segment 4 shares vertex (10,40) with segment 1, so adjacent
    # arc extents overlap — use segment 3, the lat-50 edge, and prove disjointness first.)
    spts = collect(GI.getpoint(GI.getring(_PG_SPH_POLY, 1)))
    sseg(i, j) = GO._segment_extent(GO.Spherical(),
        GO._to_kernel_point(GO.Spherical(), spts[i]),
        GO._to_kernel_point(GO.Spherical(), spts[j]))
    @test 2 in GO.SpatialTreeInterface.query(sprep.tree, sseg(2, 3))
    @test !Extents.intersects(sseg(1, 2), sseg(3, 4))
    @test 3 in GO.SpatialTreeInterface.query(sprep.tree, sseg(3, 4))
    @test 1 ∉ GO.SpatialTreeInterface.query(sprep.tree, sseg(3, 4))
end

@testset "ChildTree" begin
    # on a multipolygon: consumed at the MP node (topmost-wins), indexing polygon extents
    mp = prepare(_PG_MP; preps = (GO.ChildTree(),))
    prep = get(mp, SpatialIndexLike())
    @test prep isa GO.SpatialIndex
    # polygon 2 lives at x in (20, 25): a query box there hits only child 2
    @test GO.SpatialTreeInterface.query(prep.tree, Extents.Extent(X = (21.0, 22.0), Y = (0.0, 1.0))) == [2]
    # children did not also get trees
    @test get(GI.getgeom(mp, 1), SpatialIndexLike()) === nothing

    # on a bare polygon: indexes ring extents
    p = prepare(_PG_POLY; preps = (GO.ChildTree(),))
    rprep = get(p, SpatialIndexLike())
    @test rprep isa GO.SpatialIndex
    # the hole (ring 2) occupies (4..6, 4..6)
    @test GO.SpatialTreeInterface.query(rprep.tree, Extents.Extent(X = (4.5, 5.5), Y = (4.5, 5.5))) == [1, 2]
end
