using Test
import GeometryOps as GO
import GeometryOps: Prepared, prepare, getprep,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike
import GeometryOpsCore
import GeoInterface as GI
import GeoInterface: Extents
import LibGEOS as LG  # only for EMPTY geometries — GI wrappers cannot be empty

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

    # `buildprep(::ChildTree, ...)` SKIPS empty children, so tree leaf `i` is the `i`-th
    # NON-EMPTY child in `GI.getgeom` order, not original child `i`. Pin that index shift.
    # (GI wrappers cannot be empty — `GeoInterface.isempty` falls back to `false` for every
    # geometry trait — so the empty middle child is built with LibGEOS, matching the
    # convention in `test/methods/area.jl` and `test/methods/relateng/topology_computer.jl`.)
    mp_gap = prepare(GI.MultiPolygon([_PG_POLY, LG.readgeom("POLYGON EMPTY"), _PG_POLY2]);
        preps = (GO.ChildTree(),))
    prep_gap = get(mp_gap, SpatialIndexLike())
    @test prep_gap isa GO.SpatialIndex
    # the empty middle child is dropped, so the tree has exactly 2 leaves: a box over the
    # whole extent hits both, in tree order [1, 2]
    @test GO.SpatialTreeInterface.query(prep_gap.tree, Extents.Extent(X = (0.0, 25.0), Y = (0.0, 10.0))) == [1, 2]
    # leaf 2 is ORIGINAL child 3 (`_PG_POLY2`, x in (20, 25)): a box over its area returns [2]
    @test GO.SpatialTreeInterface.query(prep_gap.tree, Extents.Extent(X = (21.0, 22.0), Y = (0.0, 1.0))) == [2]
end

@testset "PointInArea" begin
    p = prepare(_PG_POLY; preps = (GO.PointInArea(),))
    prep = get(p, PointInAreaLike())
    @test prep isa GO.PointInAreaIndex
    # ground truth against the unindexed locator: interior, hole, boundary, exterior
    @test GO.locate(prep.locator, (1.0, 1.0)) == GO.LOC_INTERIOR
    @test GO.locate(prep.locator, (5.0, 5.0)) == GO.LOC_EXTERIOR   # inside the hole
    @test GO.locate(prep.locator, (5.0, 0.0)) == GO.LOC_BOUNDARY
    @test GO.locate(prep.locator, (11.0, 5.0)) == GO.LOC_EXTERIOR

    # multipolygon input: consumed at the MP node (relateng's element granularity)
    mp = prepare(_PG_MP; preps = (GO.PointInArea(),))
    @test get(mp, PointInAreaLike()) isa GO.PointInAreaIndex
    @test get(GI.getgeom(mp, 1), PointInAreaLike()) === nothing

    # spherical: a clear, early error (Task 16 of the spherical plan was skipped)
    @test_throws ArgumentError prepare(_PG_SPH_POLY;
        manifold = GO.Spherical(), preps = (GO.PointInArea(),))
end

@testset "relateng seam: Prepared extents are trusted" begin
    m = GO.Planar()
    pp = prepare(_PG_POLY; preps = (GO.RingEdgeIndex(),))
    rg = GO.RelateGeometry(m, pp; exact = GO.True())
    # matched manifold: the Prepared tree is used as-is (no rewrap)
    @test rg.geom === pp

    # mismatched manifold: falls back to the rewrap path (correctness over reuse)
    ps = prepare(_PG_SPH_POLY)   # Planar-prepared (2D extents)
    rgs = GO.RelateGeometry(GO.Spherical(), ps; exact = GO.True())
    @test rgs.geom !== ps
    @test GI.extent(rgs.geom) isa Extents.Extent{(:X, :Y, :Z)}

    # mismatched manifold, reverse direction: spherical-prepared input to a planar
    # relate must be stripped and rebuilt, not silently mis-pruned via cached 3D extents
    sph_prep = prepare(_PG_SPH_POLY; manifold = GO.Spherical())
    overlap_b = GI.Polygon([GI.LinearRing([(15.0, 45.0), (25.0, 45.0), (25.0, 55.0), (15.0, 55.0), (15.0, 45.0)])])
    cross_line = GI.LineString([(5.0, 39.0), (25.0, 51.0)])
    palg = GO.RelateNG()
    @test string(GO.relate(palg, sph_prep, overlap_b)) == string(GO.relate(palg, _PG_SPH_POLY, overlap_b))
    @test string(GO.relate(palg, sph_prep, cross_line)) == string(GO.relate(palg, _PG_SPH_POLY, cross_line))
    @test string(GO.relate(palg, overlap_b, sph_prep)) == string(GO.relate(palg, overlap_b, _PG_SPH_POLY))
    # the rewrap must carry freshly computed 2D extents equal to the plain rebuild's
    rg_mm = GO.RelateGeometry(GO.Planar(), sph_prep; exact = GO.True())
    rg_plain = GO.RelateGeometry(GO.Planar(), _PG_SPH_POLY; exact = GO.True())
    @test GI.extent(rg_mm.geom) == GI.extent(rg_plain.geom)
end

@testset "relateng seam: PointInArea reuse" begin
    m = GO.Planar()
    pp = prepare(_PG_POLY; preps = (GO.PointInArea(),))
    prep = get(pp, PointInAreaLike())

    # `GO.locate(::RelatePointLocator, p)` routes locate_with_dim -> compute_dim_location
    # -> locate_on_polygons -> locate_on_polygonal, which is the seam under test.

    # unprepared locator over a Prepared element uses the preparation immediately
    loc = GO.RelatePointLocator(m, pp; exact = GO.True())
    @test GO.locate(loc, (1.0, 1.0)) == GO.LOC_INTERIOR   # forces the polygonal path
    @test loc.poly_locator[1] === prep.locator             # identity: reused, not rebuilt

    # and in prepared mode too
    ploc = GO.RelatePointLocator(m, pp; exact = GO.True(), is_prepared = true)
    GO.locate(ploc, (1.0, 1.0))
    @test ploc.poly_locator[1] === prep.locator

    # (a) MultiPolygon whole-element pin: `_extract_elements!` keeps a Polygon/MultiPolygon
    # as ONE element, so the prepared MP's own PointInArea node (attached topmost-wins at the
    # MP) is the element at index 1 and is reused there — not a per-polygon rebuild.
    mp_prep = prepare(_PG_MP; preps = (GO.PointInArea(),))
    mp_loc = GO.RelatePointLocator(m, mp_prep; exact = GO.True())
    @test mp_loc.polygons[1] === mp_prep                   # the element IS the whole Prepared MP
    @test GO.locate(mp_loc, (1.0, 1.0)) == GO.LOC_INTERIOR # forces the polygonal path
    @test mp_loc.poly_locator[1] === get(mp_prep, PointInAreaLike()).locator

    # (b) Type-guard fallback: a locator with `exact = False()` cannot reuse a prep whose
    # locator was built with `exact = True()` (the `isa IndexedPointInAreaLocator{M, E}` guard
    # fails on the E mismatch), so a fresh locator is built — but results still match plain.
    guard_loc = GO.RelatePointLocator(m, pp; exact = GO.False())
    @test GO.locate(guard_loc, (1.0, 1.0)) == GO.LOC_INTERIOR
    @test guard_loc.poly_locator[1] !== prep.locator       # NOT reused: `exact` type mismatch
    @test guard_loc.poly_locator[1] !== nothing            # but a locator WAS still built
    plain_loc = GO.RelatePointLocator(m, _PG_POLY; exact = GO.False())
    for q in [(1.0, 1.0), (5.0, 5.0), (5.0, 0.0), (11.0, 5.0)]
        @test GO.locate(guard_loc, q) == GO.locate(plain_loc, q)
    end
end

# ============================================================================
# Task 10 — end-to-end equality: `relate(alg, a, b)` must be string-identical
# whether the inputs are plain or `Prepared`, for ANY preparation combination,
# in BOTH operand positions, on BOTH manifolds. A `relate` string that changes
# with preparation means a seam changed topological semantics — never
# acceptable. The sweep loops over (name, a, b, preps, position) tuples;
# `@test` expressions embed a (name, prep-label, position) tuple so a failure
# names the exact cell.
#
# Axes swept:
#   geometry type : polygon-with-hole, multipolygon, linestring, geometry
#                   collection, point, multipoint
#   preparation   : planar {(), ring, PointInArea, ring+child+PIA};
#                   spherical {(), ring, ring+child}  (no PointInArea on the
#                   sphere — it throws; PointInArea only where polygonal+planar)
#   operand pos.  : A-plain, A-prepared-algorithm, B-plain, both, both+algorithm
#   manifold      : Planar, Spherical
# ============================================================================
@testset "relate: plain == Prepared sweep" begin
    # -- planar prep-combination axis --
    poly_combos = [                                # polygonal a-side: PointInArea allowed
        ("noprep", ()),
        ("ring",   (GO.RingEdgeIndex(),)),
        ("pia",    (GO.PointInArea(),)),
        ("all",    (GO.RingEdgeIndex(), GO.ChildTree(), GO.PointInArea())),
    ]
    nonpoly_combos = [                             # non-polygonal a-side: omit PointInArea
        ("noprep",    ()),
        ("ring",      (GO.RingEdgeIndex(),)),
        ("ringchild", (GO.RingEdgeIndex(), GO.ChildTree())),
    ]

    # -- planar fixtures --
    line    = GI.LineString([(-1.0, -1.0), (5.0, 5.0), (12.0, 5.0)])  # in, through hole, out
    line2   = GI.LineString([(0.0, 5.0), (10.0, 5.0)])               # spans _PG_POLY + its hole
    pt      = GI.Point(5.0, 1.0)                                     # interior to _PG_POLY
    mpt     = GI.MultiPoint([(5.0, 1.0), (30.0, 30.0)])             # one in, one far outside
    overlap = GI.Polygon([GI.LinearRing([(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 5.0)])])
    gc      = GI.GeometryCollection([_PG_POLY2, line2])            # mixed-dimension collection
    pia3    = (GO.RingEdgeIndex(), GO.ChildTree(), GO.PointInArea())

    # each row: (name, a, b, a_combos, b_prep, both_positions)
    planar_pairs = [
        ("poly|poly disjoint", _PG_POLY, _PG_POLY2, poly_combos, pia3, true),
        ("poly|poly overlap",  _PG_POLY, overlap,   poly_combos, pia3, true),
        ("poly|line",          _PG_POLY, line,      poly_combos, (GO.RingEdgeIndex(),), true),
        ("poly|point",         _PG_POLY, pt,        poly_combos, (), false),
        ("mp|poly",            _PG_MP,   _PG_POLY,  poly_combos, pia3, false),
        ("poly|mp",            _PG_POLY, _PG_MP,    poly_combos, (GO.RingEdgeIndex(), GO.ChildTree()), true),
        ("gc|poly",            gc,       _PG_POLY,  nonpoly_combos, pia3, false),
        ("point|multipoint",   pt,       mpt,       nonpoly_combos, (), false),
    ]

    palg = GO.RelateNG()
    for (name, a, b, a_combos, b_prep, both) in planar_pairs
        @testset "$name" begin
            expected = string(GO.relate(palg, a, b))
            # axis: a-side preparation × {plain path, algorithm-prepared path}
            for (plabel, preps) in a_combos
                pa = prepare(a; preps = preps, manifold = GO.Planar())
                @test (name, plabel, :Aplain, string(GO.relate(palg, pa, b))) ==
                      (name, plabel, :Aplain, expected)
                @test (name, plabel, :Aprep, string(GO.relate(GO.prepare(palg, pa), b))) ==
                      (name, plabel, :Aprep, expected)
            end
            # axis: operand position — B-prepared and both-prepared (richest a-combo)
            if both
                pa = prepare(a; preps = last(a_combos)[2], manifold = GO.Planar())
                pb = prepare(b; preps = b_prep, manifold = GO.Planar())
                @test (name, :Bplain, string(GO.relate(palg, a, pb))) ==
                      (name, :Bplain, expected)
                @test (name, :AB, string(GO.relate(palg, pa, pb))) ==
                      (name, :AB, expected)
                @test (name, :ABprep, string(GO.relate(GO.prepare(palg, pa), pb))) ==
                      (name, :ABprep, expected)
            end
        end
    end

    # -- spherical sweep (no PointInArea on the sphere) --
    salg         = GO.RelateNG(; manifold = GO.Spherical())
    sph_overlap  = GI.Polygon([GI.LinearRing([(15.0, 45.0), (25.0, 45.0), (25.0, 55.0), (15.0, 55.0), (15.0, 45.0)])])
    sph_disjoint = GI.Polygon([GI.LinearRing([(40.0, 40.0), (50.0, 40.0), (50.0, 50.0), (40.0, 50.0), (40.0, 40.0)])])
    sph_line     = GI.LineString([(5.0, 39.0), (25.0, 51.0)])   # crosses _PG_SPH_POLY
    sph_combos   = [
        ("noprep",    ()),
        ("ring",      (GO.RingEdgeIndex(),)),
        ("ringchild", (GO.RingEdgeIndex(), GO.ChildTree())),
    ]
    spherical_pairs = [
        ("sph poly|poly overlap",  _PG_SPH_POLY, sph_overlap,  sph_combos, (GO.RingEdgeIndex(), GO.ChildTree()), true),
        ("sph poly|poly disjoint", _PG_SPH_POLY, sph_disjoint, sph_combos, (GO.RingEdgeIndex(),), false),
        ("sph poly|line",          _PG_SPH_POLY, sph_line,     sph_combos, (GO.RingEdgeIndex(),), true),
    ]
    for (name, a, b, a_combos, b_prep, both) in spherical_pairs
        @testset "$name" begin
            expected = string(GO.relate(salg, a, b))
            for (plabel, preps) in a_combos
                pa = prepare(a; preps = preps, manifold = GO.Spherical())
                @test (name, plabel, :Aplain, string(GO.relate(salg, pa, b))) ==
                      (name, plabel, :Aplain, expected)
                @test (name, plabel, :Aprep, string(GO.relate(GO.prepare(salg, pa), b))) ==
                      (name, plabel, :Aprep, expected)
            end
            if both
                pa = prepare(a; preps = last(a_combos)[2], manifold = GO.Spherical())
                pb = prepare(b; preps = b_prep, manifold = GO.Spherical())
                @test (name, :Bplain, string(GO.relate(salg, a, pb))) ==
                      (name, :Bplain, expected)
                @test (name, :AB, string(GO.relate(salg, pa, pb))) ==
                      (name, :AB, expected)
                @test (name, :ABprep, string(GO.relate(GO.prepare(salg, pa), pb))) ==
                      (name, :ABprep, expected)
            end
        end
    end

    # -- REQUIRED: Spherical GeometryCollection with a BARE Point + a Polygon (Task 4 review).
    # The legacy plain-input extent path (`_union_stored_extents`) reads the bare point's
    # extent via `GI.extent(pt; fallback=true)` — a manifold-BLIND 2D lon/lat box — and unions
    # it with the polygon's 3D unit-sphere box; `Extents.union` keeps only shared keys, so the
    # GC WRAPPER caches a dimensionally-collapsed, coordinate-mixed 2D garbage extent.
    # `prepare`'s `_union_prepared_extents` instead uses `rk_interaction_bounds(m, pt)` (a real
    # 3D box), matching `_relate_extent` exactly. The relate STRING is nonetheless IDENTICAL:
    # the garbage is shadowed — `RelateGeometry.extent` (the field driving the envelope
    # short-circuit) is recomputed to the correct 3D box by `_relate_extent` in BOTH paths. So
    # plain == prepared, and both are correct (ground truth `212101212`/overlap, `FF2FF1212`/
    # disjoint verified against the GC's members related separately + planar intuition; the
    # interior point is absorbed by the polygon, so the GC relates to b exactly as the polygon
    # does). This is NOT the latent legacy bug reaching an output — it stays contained.
    @testset "spherical GC{bare point, polygon}" begin
        gc_pt = GI.GeometryCollection([GI.Point(12.0, 42.0), _PG_SPH_POLY])  # point interior to poly
        # Pin that this fixture actually exercises the divergent extent path (else the equality
        # below would be a vacuous no-op): legacy caches a collapsed 2D GC extent, prepare a 3D one.
        legacy_gc_ext = GI.extent(GO._relate_cache_extents(GO.Spherical(), gc_pt))
        prep_gc_ext   = prepare(gc_pt; manifold = GO.Spherical()).extent
        @test keys(legacy_gc_ext) == (:X, :Y)          # legacy: dimensionally collapsed (the bug)
        @test keys(prep_gc_ext)   == (:X, :Y, :Z)      # prepare: correct 3D unit-sphere box

        for (blabel, b) in [(:overlap, sph_overlap), (:disjoint, sph_disjoint)]
            expected = string(GO.relate(salg, gc_pt, b))
            pgc = prepare(gc_pt; manifold = GO.Spherical(), preps = (GO.RingEdgeIndex(),))
            pb  = prepare(b; manifold = GO.Spherical(), preps = (GO.RingEdgeIndex(),))
            @test (blabel, :Aplain, string(GO.relate(salg, pgc, b))) == (blabel, :Aplain, expected)
            @test (blabel, :Aprep, string(GO.relate(GO.prepare(salg, pgc), b))) == (blabel, :Aprep, expected)
            @test (blabel, :Bprep, string(GO.relate(salg, gc_pt, pb))) == (blabel, :Bprep, expected)
        end
    end
end
