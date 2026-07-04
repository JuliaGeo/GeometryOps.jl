using Test

import GeoInterface as GI
import GeometryOps as GO
import GeometryOps: prepare, Prepared, getprep, EdgeTree, AbstractEdgeTree,
    AbstractPreparation, build_edge_tree, edge_tree
import GeometryOps.NaturalIndexing: NaturalIndex
import GeometryOps.SpatialTreeInterface: FlatNoTree
import GeometryOpsCore: True, False
import Extents

import LibGEOS as LG
import ArchGDAL as AG
import GeoJSON
using GeometryOpsTestHelpers

# A user-defined preparation: subtype + any way to build it (here, closures).
struct TagPrep <: AbstractPreparation
    tag::Symbol
end

# --- Fixture geometries ------------------------------------------------------

square_hole = GI.Polygon([
    [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
    [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)],
])
# Same shape, but with unclosed rings (the implicit closing edge must be indexed too).
square_hole_unclosed = GI.Polygon([
    [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)],
])
# Horizontal collinear runs exercise the degenerate ray cases of Hao–Sun.
horiz_degenerate = GI.Polygon([[
    (0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (6.0, 2.0), (8.0, 2.0), (8.0, 4.0),
    (2.0, 4.0), (0.0, 4.0), (0.0, 0.0),
]])
circle = let
    pts = [(5 + 5cos(θ), 5 + 5sin(θ)) for θ in range(0, 2π; length = 257)[1:256]]
    push!(pts, pts[1])
    GI.Polygon([pts])
end
multipoly = GI.MultiPolygon([square_hole, GI.Polygon([[(20.0, 0.0), (25.0, 0.0), (25.0, 5.0), (20.0, 0.0)]])])

# Points that stress every code path: a grid spilling past the extent, every
# vertex, points at vertex height just left/right of it (ray through vertex),
# and every edge midpoint (on-boundary hits).
function probe_points(poly; n = 15)
    ext = GI.extent(poly)
    (x1, x2), (y1, y2) = ext.X, ext.Y
    dx, dy = x2 - x1, y2 - y1
    pts = NTuple{2, Float64}[]
    for x in range(x1 - 0.1dx, x2 + 0.1dx; length = n), y in range(y1 - 0.1dy, y2 + 0.1dy; length = n)
        push!(pts, (x, y))
    end
    for ring in GI.getring(poly)
        np = GI.npoint(ring)
        for i in 1:np
            p = GI.getpoint(ring, i)
            v = (Float64(GI.x(p)), Float64(GI.y(p)))
            push!(pts, v, (v[1] - 0.21dx, v[2]), (v[1] + 0.21dx, v[2]))
            if i < np
                q = GI.getpoint(ring, i + 1)
                push!(pts, ((v[1] + GI.x(q)) / 2, (v[2] + GI.y(q)) / 2))
            end
        end
    end
    return pts
end

# Compare every point-vs-polygon predicate between a plain and a prepared
# polygon; returns a vector of mismatch descriptions (empty = equivalent).
function predicate_mismatches(poly, prep, pts)
    bad = Any[]
    for pt in pts
        for f in (GO.within, GO.coveredby, GO.disjoint, GO.intersects, GO.touches)
            a, b = f(pt, poly), f(pt, prep)
            a == b || push!(bad, (; f, pt, plain = a, prepared = b))
        end
        for f in (GO.contains, GO.covers)
            a, b = f(poly, pt), f(prep, pt)
            a == b || push!(bad, (; f, pt, plain = a, prepared = b))
        end
    end
    return bad
end

@testset "Prepared materializes into a transparent GeoInterface wrapper" begin
    prep = prepare(square_hole)
    @test GI.geomtrait(prep) isa GI.PolygonTrait
    @test GI.isgeometry(typeof(prep))
    @test GI.ngeom(prep) == GI.ngeom(square_hole)
    @test GI.npoint(prep) == GI.npoint(square_hole)
    @test GI.nhole(prep) == 1
    @test !GI.is3d(prep)
    e1, e2 = GI.extent(prep), GI.extent(square_hole)
    @test e1.X == e2.X && e1.Y == e2.Y
    @test Extents.extent(prep) === GI.extent(prep)
    @test GO.equals(prep, square_hole)
    @test GO.area(prep) == GO.area(square_hole)
    @test GO.centroid(prep) == GO.centroid(square_hole)

    # Materialization: `parent` is the converted geometry (tuple storage), not
    # the original object, and the children are themselves `Prepared` nodes.
    @test Base.parent(prep) isa GI.Polygon
    @test Base.parent(prep) !== square_hole
    ring = GI.getexterior(prep)
    @test ring isa Prepared
    @test GI.geomtrait(ring) isa GI.LinearRingTrait
    @test GI.getpoint(ring, 1) === (0.0, 0.0)
    @test collect(GI.gethole(prep))[1] isa Prepared
    # Ring nodes carry their edge tree; the polygon node has no preps of its own.
    @test getprep(ring, AbstractEdgeTree) isa EdgeTree
    @test prep.preps === ()
    @test occursin("Prepared", sprint(show, prep))
    @test occursin("EdgeTree", sprint(show, ring))

    # Preparedness survives decomposition of multi-geometries too.
    prep_mp = prepare(multipoly)
    @test GI.getgeom(prep_mp, 1) isa Prepared
    @test getprep(GI.getexterior(GI.getgeom(prep_mp, 1)), AbstractEdgeTree) isa EdgeTree
    e_mp = GI.extent(prep_mp)
    @test e_mp.X == GI.extent(multipoly).X && e_mp.Y == GI.extent(multipoly).Y

    # Point wrapping hits the point-trait disambiguators.
    ppt = prepare(GI.Point(1.0, 2.0))
    @test GI.geomtrait(ppt) isa GI.PointTrait
    @test GI.ngeom(ppt) == 0
    @test GI.x(ppt) == 1.0 && GI.y(ppt) == 2.0

    # Constructor guardrails.
    @test_throws ArgumentError Prepared(prep, (), nothing)
    @test_throws ArgumentError Prepared("not a geometry", (), nothing)
    @test_throws ArgumentError prepare("not a geometry")
end

@testset "Coordinate number types are preserved" begin
    poly32 = GI.Polygon([[(0.0f0, 0.0f0), (10.0f0, 0.0f0), (10.0f0, 10.0f0), (0.0f0, 10.0f0), (0.0f0, 0.0f0)]])
    prep32 = prepare(poly32)
    @test GI.getpoint(GI.getexterior(prep32), 1) === (0.0f0, 0.0f0)
    polyint = GI.Polygon([[(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]])
    prepint = prepare(polyint)
    @test GI.getpoint(GI.getexterior(prepint), 1) === (0, 0)
    for (plain, prep) in ((poly32, prep32), (polyint, prepint))
        for pt in ((5.0, 5.0), (0.0, 0.0), (11.0, 5.0), (5.0, 0.0), (10.0, 10.0))
            @test GO.within(pt, prep) == GO.within(pt, plain)
        end
    end
end

@testset "Materialization closes rings and preserves UnitSphericalPoints" begin
    # An unclosed input ring is closed during materialization, so every
    # preparation can rely on point 1 == point n.
    unclosed = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]])
    ring = GI.getexterior(prepare(unclosed))
    @test GI.npoint(ring) == 5
    @test GI.getpoint(ring, 1) == GI.getpoint(ring, 5)
    @test GO.area(prepare(unclosed)) == GO.area(unclosed)

    # Points already in unit-spherical representation are stored as-is, not
    # destructured to coordinate tuples; prepare takes the manifold first,
    # like other GeometryOps functions.
    USP = GO.UnitSpherical.UnitSphericalPoint
    pts = [USP(1.0, 0.0, 0.0), USP(0.0, 1.0, 0.0), USP(0.0, 0.0, 1.0)]
    usp_prep = prepare(GO.Spherical(), GI.LinearRing(pts))
    stored = collect(GI.getpoint(parent(usp_prep)))
    @test eltype(stored) <: USP
    @test length(stored) == 4 && stored[end] == stored[1]   # also closed
    # Edge trees are a planar default; the spherical ones land later.
    @test isnothing(getprep(usp_prep, AbstractEdgeTree))
    @test GO.hasprep(GI.getexterior(prepare(GO.Planar(), unclosed)), AbstractEdgeTree)
    # A bare unit-spherical point round-trips as itself.
    @test parent(prepare(GO.Spherical(), USP(0.0, 0.0, 1.0))) isa USP
end

@testset "Polygon rings that report LineStringTrait still get edge trees" begin
    # GeoJSON types polygon rings as line strings; materialization must treat
    # polygon children as rings regardless, or they silently lose their index.
    gj = GeoJSON.read("""{"type": "Polygon", "coordinates": [
        [[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
        [[3, 3], [7, 3], [7, 7], [3, 7], [3, 3]]]}""")
    prep = prepare(gj)
    for ring in GI.getring(prep)
        @test GI.geomtrait(ring) isa GI.LinearRingTrait
        @test getprep(ring, AbstractEdgeTree) isa EdgeTree
    end
    for pt in ((5.0, 5.0), (5.0, 4.9), (0.0, 0.0), (11.0, 5.0), (5.0, 0.0))
        @test GO.within(pt, prep) == GO.within(pt, gj)
    end
end

@testset "getprep: query, get-or-else, precedence" begin
    prep = prepare(square_hole)
    ring = GI.getexterior(prep)
    # Plain geometries always miss, so consumer code needs no special-casing.
    @test getprep(GI.getexterior(square_hole), AbstractEdgeTree) === nothing
    @test getprep(square_hole, AbstractPreparation) === nothing
    # Concrete, abstract-kind, and root-abstract queries all find the prep —
    # and resolve at compile time from the preparation tuple's type, so hits
    # and misses both infer concretely.
    @test @inferred(getprep(ring, EdgeTree)) isa EdgeTree
    @test @inferred(getprep(ring, AbstractEdgeTree)) isa EdgeTree
    @test @inferred(getprep(ring, AbstractPreparation)) isa EdgeTree
    @test @inferred(getprep(ring, TagPrep)) === nothing
    @test @inferred(GO.hasprep(ring, AbstractEdgeTree))
    # Edge trees live on rings, not on the polygon node.
    @test getprep(prep, AbstractEdgeTree) === nothing
    # get-or-else: `f` runs only on a miss.
    @test getprep(() -> TagPrep(:built), square_hole, TagPrep) == TagPrep(:built)
    @test getprep(() -> error("must not build"), ring, EdgeTree) isa EdgeTree

    # Adding preparations to an existing Prepared: storage unchanged, new preps win.
    prep2 = prepare(prep; preps = (g -> TagPrep(:added),))
    @test Base.parent(prep2) === Base.parent(prep)
    @test getprep(prep2, TagPrep) == TagPrep(:added)
    @test getprep(prep2, AbstractPreparation) isa TagPrep  # prepended ⇒ found first
    @test prepare(prep) === prep  # no-op add returns the same object

    # A spec without an `appliesto` declaration applies to the top node only;
    # nodes where no given spec applies still get defaults.
    prep3 = prepare(square_hole; preps = (g -> TagPrep(:spec),))
    @test getprep(prep3, TagPrep) == TagPrep(:spec)
    @test getprep(GI.getexterior(prep3), AbstractEdgeTree) isa EdgeTree

    # A callable selector applies at every node of the recursion.
    prep4 = prepare(square_hole; preps = (t, g) -> ())
    @test prep4.preps === ()
    @test getprep(GI.getexterior(prep4), AbstractEdgeTree) === nothing

    # Non-polygon nodes: no preps of their own, but prepared children.
    prep_mp = prepare(multipoly)
    @test prep_mp.preps === ()
    @test getprep(prep_mp, AbstractEdgeTree) === nothing
end

@testset "EdgeTree construction and backends" begin
    ring = GI.getexterior(square_hole)
    nat = EdgeTree(ring)
    @test edge_tree(nat) isa NaturalIndex
    str = EdgeTree(ring; backend = GO.STRtree)
    @test edge_tree(str) isa GO.STRtree
    flat = EdgeTree(ring; backend = r -> FlatNoTree(GO._edge_extents(r)))
    @test edge_tree(flat) isa FlatNoTree
    @test_throws ArgumentError EdgeTree(square_hole)
    # The curried spec form: `EdgeTree(backend)` applies to every curve of
    # the recursion and nowhere else.
    sel = EdgeTree(GO.STRtree)
    @test GO.appliesto(sel, GI.LinearRingTrait(), false)
    @test !GO.appliesto(sel, GI.PolygonTrait(), true)
    prep = prepare(square_hole; preps = (sel,))
    @test edge_tree(getprep(GI.getexterior(prep), AbstractEdgeTree)) isa GO.STRtree
    @test prep.preps === ()   # the polygon node itself gets no prep from it
    # Unclosed rings index the implicit closing edge.
    @test length(GO._edge_extents(GI.getexterior(square_hole))) == 4
    @test length(GO._edge_extents(GI.getexterior(square_hole_unclosed))) == 4
end

@testset "Indexed point-in-polygon ≡ plain point-in-polygon" begin
    backends = (
        "NaturalIndex" => poly -> prepare(poly),
        "STRtree" => poly -> prepare(poly; preps = (EdgeTree(GO.STRtree),)),
        "FlatNoTree" => poly -> prepare(poly; preps = (EdgeTree(r -> FlatNoTree(GO._edge_extents(r))),)),
    )
    polys = (
        "square with hole" => square_hole,
        "unclosed rings" => square_hole_unclosed,
        "horizontal degeneracies" => horiz_degenerate,
        "circle-256" => circle,
    )
    for (pname, poly) in polys
        pts = probe_points(poly)
        for (bname, prepper) in backends
            @testset "$pname / $bname" begin
                @test predicate_mismatches(poly, prepper(poly), pts) == []
            end
        end
    end

    # Ring-level equivalence, under both exactness settings.
    @testset "ring orientation kernel equivalence" begin
        for poly in (square_hole, horiz_degenerate, circle), exact in (True(), False())
            ring = GI.getexterior(poly)
            tree = build_edge_tree(NaturalIndex, ring)
            bad = [pt for pt in probe_points(poly)
                if GO._point_filled_curve_orientation(GO.Planar(), pt, ring; exact) !==
                   GO._point_filled_curve_orientation(GO.Planar(), pt, ring, tree; exact)]
            @test bad == []
        end
    end

    # Extent-cache-only preparations (no trees) also agree — the `nothing` path.
    prep_extent_only = prepare(square_hole; preps = (t, g) -> ())
    @test predicate_mismatches(square_hole, prep_extent_only, probe_points(square_hole)) == []

    # MultiPolygons are now accelerated through their prepared children.
    mp_prep = prepare(multipoly)
    mp_pts = vcat(probe_points(GI.getgeom(multipoly, 1)), probe_points(GI.getgeom(multipoly, 2)))
    @test predicate_mismatches(multipoly, mp_prep, mp_pts) == []
end

@testset "Clipping reuses prepared ring edge trees" begin
    ngon(c, r, n) = begin
        pts = [(c[1] + r * cos(θ), c[2] + r * sin(θ)) for θ in range(0, 2π; length = n + 1)[1:n]]
        push!(pts, pts[1])
        pts
    end
    donut = GI.Polygon([ngon((5.0, 5.0), 5.0, 64), ngon((5.0, 5.0), 2.0, 16)])
    blob = GI.Polygon([ngon((9.0, 5.0), 5.0, 64)])
    small = GI.Polygon([ngon((7.0, 5.0), 3.0, 8)])

    # A prepared ring's edge tree is reused as-is (indices match `eachedge`),
    # and coordinates read from ring storage match the materialized edge list.
    prep_ring = GI.getexterior(prepare(donut))
    ct, n = GO._curve_trees(prep_ring)
    @test only(ct).tree === edge_tree(getprep(prep_ring, AbstractEdgeTree))   # reused, not rebuilt
    edges = GO.to_edgelist(GI.getexterior(donut), Float64)
    @test n == length(edges)
    for j in (1, 2, n)
        @test GO._edge_coords(ct, j, Float64) == Tuple(edges[j].geom)
    end
    # Any SpatialTreeInterface tree is reused — including one that traverses
    # out of input order, like HPR (a stand-in for e.g. a foreign-library tree).
    hpr_ring = GI.getexterior(prepare(donut; preps = (EdgeTree(GO.FlexibleRTrees.HPR()),)))
    @test only(first(GO._curve_trees(hpr_ring))).tree isa GO.FlexibleRTrees.RTree
    # Unclosed input rings are closed during materialization, so a prepared
    # ring's tree always matches the `eachedge` space and is reused as-is —
    # the closing edge is explicit, never implicit.
    unclosed_ring = GI.getexterior(prepare(GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]])))
    @test GI.npoint(unclosed_ring) == 5   # the first point was repeated at the end
    ct_uc, n_uc = GO._curve_trees(unclosed_ring)
    @test only(ct_uc).tree === edge_tree(getprep(unclosed_ring, AbstractEdgeTree))   # reused
    @test n_uc == 4   # the closing edge is an explicit `eachedge` edge now

    # Prepared clipping ≡ plain nested-loop clipping, for every combination of
    # prepared/plain inputs and backends, including holes and small polygons
    # (prepared inputs override the size heuristic in AutoAccelerator).
    auto = GO.FosterHormannClipping(GO.Planar(), GO.AutoAccelerator())
    plain_alg = GO.FosterHormannClipping()   # NestedLoop ground truth
    prep_nat = prepare
    prep_hpr = g -> prepare(g; preps = (EdgeTree(GO.FlexibleRTrees.HPR()),))
    for (pa, pb) in ((donut, blob), (donut, small), (small, blob))
        for prepper in (prep_nat, prep_hpr)
            for f in (GO.intersection, GO.union, GO.difference)
                expected = f(plain_alg, pa, pb; target = GI.PolygonTrait())
                for (A, B) in ((prepper(pa), pb), (pa, prepper(pb)), (prepper(pa), prepper(pb)))
                    got = f(auto, A, B; target = GI.PolygonTrait())
                    @test length(got) == length(expected)
                    @test all(map(GO.equals, got, expected))
                end
            end
        end
    end
end

@testset "Ring-aware foreach_pair and line-string edge trees" begin
    ngon(c, r, n) = begin
        pts = [(c[1] + r * cos(θ), c[2] + r * sin(θ)) for θ in range(0, 2π; length = n + 1)[1:n]]
        push!(pts, pts[1])
        pts
    end
    donut = GI.Polygon([ngon((5.0, 5.0), 5.0, 64), ngon((5.0, 5.0), 2.0, 16)])
    blob = GI.Polygon([ngon((9.0, 5.0), 5.0, 64)])
    auto = GO.FosterHormannClipping(GO.Planar(), GO.AutoAccelerator())
    plain_alg = GO.FosterHormannClipping()   # NestedLoop ground truth

    # `hasprep` is the node-level boolean companion to `getprep`.
    prep = prepare(donut)
    @test GO.hasprep(GI.getexterior(prep), AbstractEdgeTree)
    @test !GO.hasprep(prep, AbstractEdgeTree)      # edge trees live on the curves
    @test !GO.hasprep(donut, AbstractEdgeTree)     # plain geometries carry no preps

    # Line strings get edge trees indexing exactly their `eachedge` pairs —
    # there is no closing edge to index.
    ls = GI.LineString([(0.0, 3.0), (4.0, 8.0), (8.0, 2.0), (12.0, 7.0)])
    pls = prepare(ls)
    @test GO.hasprep(pls, AbstractEdgeTree)
    ct, n = GO._curve_trees(pls)
    @test only(ct).tree === edge_tree(getprep(pls, AbstractEdgeTree))
    @test n == 3
    ls_edges = GO.to_edgelist(ls, Float64)
    @test all(GO._edge_coords(ct, j, Float64) == Tuple(ls_edges[j].geom) for j in 1:n)
    @test length(GO._edge_extents(ls)) == 3
    @test length(GO._edge_extents(GI.LinearRing(collect(GI.getpoint(ls))))) == 4  # a raw unclosed ring still wraps

    # Point-in-polygon consumers rely on the materialization contract: a
    # polygon prepared whole gets closed rings, so PIP through its (formerly
    # open) boundary matches the plain polygon.
    open_poly = GI.Polygon([collect(GI.getpoint(ls))])
    prep_open_poly = prepare(open_poly)
    for pt in ((6.0, 5.0), (1.0, 6.0), (20.0, 5.0))
        @test GO.within(pt, prep_open_poly) == GO.within(pt, open_poly)
    end

    # Whole geometries decompose into per-curve trees whose offsets match the
    # geometry-global `eachedge` numbering.
    parts, ntot = GO._curve_trees(prep)
    @test length(parts) == 2 && ntot == 64 + 16
    @test parts[1].offset == 0 && parts[2].offset == 64
    donut_edges = GO.to_edgelist(donut, Float64)
    @test length(donut_edges) == ntot
    @test all(GO._edge_coords(parts, j, Float64) == Tuple(donut_edges[j].geom) for j in (1, 64, 65, 80))

    # `intersection_points` on whole geometries: the tree accelerators (parts
    # path) agree with the nested loop, plain and prepared, and the blob
    # crosses the donut's hole so both rings contribute points.
    expected = GO.intersection_points(plain_alg, donut, blob)
    @test !isempty(expected)
    for acc in (GO.SingleNaturalTree(), GO.DoubleNaturalTree())
        @test GO.intersection_points(GO.FosterHormannClipping(GO.Planar(), acc), donut, blob) == expected
    end
    # Every per-side policy shares the `TreePolicy` supertype, and side `b`
    # must carry a tree-building one.
    @test all(p isa GO.TreePolicy for p in (GO.IterateEdges(), GO.BuildTree()))
    @test_throws ArgumentError GO.TreeAccelerator(GO.BuildTree(), GO.IterateEdges())
    @test GO.intersection_points(auto, prepare(donut), prepare(blob)) == expected
    @test GO.intersection_points(auto, prepare(donut), blob) == expected

    # MultiPolygons too — the far-away member is pruned by part extents.
    mp = GI.MultiPolygon([donut, GI.Polygon([ngon((-30.0, 5.0), 2.5, 32)])])
    expected_mp = GO.intersection_points(plain_alg, mp, blob)
    @test GO.intersection_points(auto, prepare(mp), prepare(blob)) == expected_mp
    @test GO.intersection_points(GO.FosterHormannClipping(GO.Planar(), GO.DoubleNaturalTree()), mp, blob) == expected_mp

    # Prepared line strings accelerate on either side.  Ground truth is
    # computed with matching argument order — `_intersection_point` is not
    # bit-symmetric in its operands.
    expected_lp = GO.intersection_points(plain_alg, donut, ls)
    @test !isempty(expected_lp)
    @test GO.intersection_points(auto, donut, pls) == expected_lp                                  # b side: single-tree fast path
    @test GO.intersection_points(auto, pls, donut) == GO.intersection_points(plain_alg, ls, donut) # a side: parts path

    # `cut` goes through the same machinery with a prepared polygon.
    circle = GI.Polygon([ngon((5.0, 5.0), 5.0, 32)])
    cut_line = GI.Line([(-1.0, 5.0), (11.0, 5.0)])
    expected_cut = GO.cut(plain_alg, circle, cut_line)
    got_cut = GO.cut(auto, prepare(circle), cut_line)
    @test length(got_cut) == length(expected_cut)
    @test all(map(GO.equals, got_cut, expected_cut))
end

@testset_implementations "Prepared vs plain across implementations" begin
    # Materialization from each backend's native storage must be exact.
    prep = prepare($square_hole)
    @test GI.getexterior(prep) isa Prepared
    for pt in ((5.0, 5.0), (5.0, 3.0), (0.0, 0.0), (11.0, 5.0), (5.0, 0.0), (5.0, 4.9))
        @test GO.within(pt, prep) == GO.within(pt, $square_hole)
        @test GO.contains(prep, pt) == GO.contains($square_hole, pt)
        @test GO.disjoint(pt, prep) == GO.disjoint(pt, $square_hole)
    end
end
