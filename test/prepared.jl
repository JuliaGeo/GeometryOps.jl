using Test

import GeoInterface as GI
import GeometryOps as GO
import GeometryOps: prepare, Prepared, getprep, RingEdgeTrees, AbstractRingEdgeTrees,
    AbstractPreparation, build_edge_tree
import GeometryOps.NaturalIndexing: NaturalIndex
import GeometryOps.SpatialTreeInterface: FlatNoTree
import GeometryOpsCore: True, False
import Extents

import LibGEOS as LG
import ArchGDAL as AG
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

@testset "Prepared is a transparent GeoInterface wrapper" begin
    prep = prepare(square_hole)
    @test GI.geomtrait(prep) isa GI.PolygonTrait
    @test GI.isgeometry(typeof(prep))
    @test GI.ngeom(prep) == GI.ngeom(square_hole)
    @test GI.npoint(prep) == GI.npoint(square_hole)
    @test GI.nhole(prep) == 1
    @test !GI.is3d(prep)
    @test Base.parent(prep) === square_hole
    e1, e2 = GI.extent(prep), GI.extent(square_hole)
    @test e1.X == e2.X && e1.Y == e2.Y
    @test Extents.extent(prep) === GI.extent(prep)
    @test GO.equals(prep, square_hole)
    @test GO.area(prep) == GO.area(square_hole)
    @test GO.centroid(prep) == GO.centroid(square_hole)
    shown = sprint(show, prep)
    @test occursin("Prepared", shown) && occursin("RingEdgeTrees", shown)

    # Point wrapping hits the point-trait disambiguators.
    ppt = prepare(GI.Point(1.0, 2.0))
    @test GI.geomtrait(ppt) isa GI.PointTrait
    @test GI.ngeom(ppt) == 0
    @test GI.x(ppt) == 1.0 && GI.y(ppt) == 2.0

    # Constructor guardrails.
    @test_throws ArgumentError Prepared(prep, (), nothing)
    @test_throws ArgumentError Prepared("not a geometry", (), nothing)
end

@testset "getprep: query, get-or-else, precedence" begin
    prep = prepare(square_hole)
    # Plain geometries always miss, so consumer code needs no special-casing.
    @test getprep(square_hole, AbstractRingEdgeTrees) === nothing
    @test getprep(square_hole, AbstractPreparation) === nothing
    # Concrete, abstract-kind, and root-abstract queries all find the prep.
    @test getprep(prep, RingEdgeTrees) isa RingEdgeTrees
    @test getprep(prep, AbstractRingEdgeTrees) isa RingEdgeTrees
    @test getprep(prep, AbstractPreparation) isa RingEdgeTrees
    @test getprep(prep, TagPrep) === nothing
    # get-or-else: `f` runs only on a miss.
    @test getprep(() -> TagPrep(:built), square_hole, TagPrep) == TagPrep(:built)
    @test getprep(() -> error("must not build"), prep, RingEdgeTrees) isa RingEdgeTrees

    # Adding preparations to an existing Prepared: parent unchanged, new preps win.
    prep2 = prepare(prep; preps = (g -> TagPrep(:added),))
    @test Base.parent(prep2) === square_hole
    @test getprep(prep2, TagPrep) == TagPrep(:added)
    @test getprep(prep2, RingEdgeTrees) isa RingEdgeTrees
    @test getprep(prep2, AbstractPreparation) isa TagPrep  # prepended ⇒ found first
    @test prepare(prep) === prep  # no-op add returns the same object

    # Custom specs flow through prepare: closures work with no registration.
    prep3 = prepare(square_hole; preps = (g -> TagPrep(:spec), RingEdgeTrees))
    @test getprep(prep3, TagPrep) == TagPrep(:spec)
    @test getprep(prep3, AbstractRingEdgeTrees) isa RingEdgeTrees

    # Non-polygon defaults: extent-only wrapper, predicates still work.
    prep_mp = prepare(multipoly)
    @test prep_mp.preps === ()
    @test getprep(prep_mp, AbstractRingEdgeTrees) === nothing
end

@testset "RingEdgeTrees construction and backends" begin
    nat = RingEdgeTrees(square_hole)
    @test nat.exterior isa NaturalIndex
    @test length(nat.holes) == 1
    str = RingEdgeTrees(square_hole; tree = GO.STRtree)
    @test str.exterior isa GO.STRtree
    flat = RingEdgeTrees(square_hole; tree = r -> FlatNoTree(GO._ring_edge_extents(r)))
    @test flat.exterior isa FlatNoTree
    @test_throws ArgumentError RingEdgeTrees(GI.LineString([(0.0, 0.0), (1.0, 1.0)]))
    # Unclosed rings index the implicit closing edge.
    @test length(GO._ring_edge_extents(GI.getexterior(square_hole))) == 4
    @test length(GO._ring_edge_extents(GI.getexterior(square_hole_unclosed))) == 4
end

@testset "Indexed point-in-polygon ≡ plain point-in-polygon" begin
    backends = (
        "NaturalIndex" => identity,  # default prepare
        "STRtree" => g -> prepare(g; preps = (h -> RingEdgeTrees(h; tree = GO.STRtree),)),
        "FlatNoTree" => g -> prepare(g; preps = (h -> RingEdgeTrees(h; tree = r -> FlatNoTree(GO._ring_edge_extents(r))),)),
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
            prep = bname == "NaturalIndex" ? prepare(poly) : prepper(poly)
            @testset "$pname / $bname" begin
                @test predicate_mismatches(poly, prep, pts) == []
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

    # Extent-only preparations (no trees) also agree — exercises the `nothing` path.
    prep_extent_only = prepare(square_hole; preps = ())
    @test predicate_mismatches(square_hole, prep_extent_only, probe_points(square_hole)) == []

    # MultiPolygon wrappers agree (no acceleration, pure forwarding).
    mp_prep = prepare(multipoly)
    mp_pts = probe_points(GI.getgeom(multipoly, 1))
    @test predicate_mismatches(multipoly, mp_prep, mp_pts) == []
end

@testset_implementations "Prepared vs plain across implementations" begin
    prep = prepare($square_hole)
    for pt in ((5.0, 5.0), (5.0, 3.0), (0.0, 0.0), (11.0, 5.0), (5.0, 0.0), (5.0, 4.9))
        @test GO.within(pt, prep) == GO.within(pt, $square_hole)
        @test GO.contains(prep, pt) == GO.contains($square_hole, pt)
        @test GO.disjoint(pt, prep) == GO.disjoint(pt, $square_hole)
    end
end
