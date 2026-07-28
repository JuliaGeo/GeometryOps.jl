# Tests for the OverlayNG phase-3 point paths (design §4): `_overlay_points`
# (point × point) and `_overlay_mixed_points` (point × line, point × area), plus
# the `OverlayUtil` result-dimension / empty-result rules they exercise through
# the `_overlay_ng` driver.
#
# Structure:
#   1. ported JTS `OverlayNGPointsTest` cases (the floating-safe subset)
#   2. ported JTS `OverlayNGMixedPointsTest` cases (the floating-safe subset)
#   3. result-dimension and empty/degenerate inputs
#   4. spherical equivalents, including cases where a planar even-odd or
#      planar on-segment test gives the *wrong* answer — these pin the
#      requirement that location goes through the kernel locators
#
# Equality: expected WKT compared with GEOS topological `equals` (order- and
# merge-independent); empty results asserted by point count and dimension, since
# GEOS `equals` cannot distinguish empties of different dimension.

using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import GeometryOps: Planar, Spherical, True

const EX = True()

lgc(g) = GI.convert(LG, g)
wkt(s) = LG.readgeom(s)

const OPS = (GO.OVERLAY_INTERSECTION, GO.OVERLAY_UNION,
             GO.OVERLAY_DIFFERENCE, GO.OVERLAY_SYMDIFFERENCE)

ov(m, op, a, b) = GO._overlay_ng(m, op, a, b; exact = EX)

# Result equals the expected WKT (topologically).
eqw(r, ewkt) = LG.equals(lgc(r), wkt(ewkt))

# GeoInterface has no `npoint`/`getpoint` for `PointTrait`; these fill that in.
npts(g) = GI.trait(g) isa GI.PointTrait ? (GI.isempty(g) ? 0 : 1) : GI.npoint(g)
coords(g) = GI.trait(g) isa GI.PointTrait ? [(GI.x(g), GI.y(g))] :
            [(GI.x(p), GI.y(p)) for p in GI.getpoint(g)]

# An empty result of the given dimension (2 → polygonal, 1 → lineal, 0 → puntal).
function is_empty_dim(r, dim::Integer)
    npts(r) == 0 || return false
    t = GI.trait(r)
    dim == 0 && return t isa GI.MultiPointTrait
    dim == 1 && return t isa GI.MultiLineStringTrait
    return t isa GI.MultiPolygonTrait
end

# ---------------------------------------------------------------------------
# 1. Ported JTS OverlayNGPointsTest (point × point)
# ---------------------------------------------------------------------------
#
# SKIPPED, with reasons (fixed-precision cases that cannot apply to a
# floating-only, non-snapping engine — design §0):
#   testDisjointPointsRoundedIntersection — POINT(10.1 10) ∩ POINT(10 10.1) is
#     only nonempty after rounding both to the unit grid.
# testSimpleMergeIntersection is kept, restated at floating precision (the
# extra A points simply do not participate).

@testset "JTS OverlayNGPointsTest (point × point)" begin
    @testset "SimpleIntersection" begin
        a = wkt("MULTIPOINT ((1 1), (2 1))"); b = wkt("POINT (2 1)")
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, a, b), "POINT (2 1)")
    end
    @testset "SimpleMergeIntersection (floating)" begin
        a = wkt("MULTIPOINT ((1 1), (1.5 1.1), (2 1), (2.1 1.1))"); b = wkt("POINT (2 1)")
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, a, b), "POINT (2 1)")
    end
    @testset "SimpleUnion" begin
        a = wkt("MULTIPOINT ((1 1), (2 1))"); b = wkt("POINT (2 1)")
        @test eqw(ov(Planar(), GO.OVERLAY_UNION, a, b), "MULTIPOINT ((1 1), (2 1))")
    end
    @testset "SimpleDifference" begin
        a = wkt("MULTIPOINT ((1 1), (2 1))"); b = wkt("POINT (2 1)")
        @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, a, b), "POINT (1 1)")
    end
    @testset "SimpleSymDifference" begin
        a = wkt("MULTIPOINT ((1 2), (1 1), (2 2), (2 1))")
        b = wkt("MULTIPOINT ((2 2), (2 1), (3 2), (3 1))")
        @test eqw(ov(Planar(), GO.OVERLAY_SYMDIFFERENCE, a, b),
                  "MULTIPOINT ((1 2), (1 1), (3 2), (3 1))")
    end
    @testset "SimpleFloatUnion" begin
        a = wkt("MULTIPOINT ((1 1), (1.5 1.1), (2 1), (2.1 1.1))")
        b = wkt("MULTIPOINT ((1.5 1.1), (2 1), (2 1.2))")
        @test eqw(ov(Planar(), GO.OVERLAY_UNION, a, b),
                  "MULTIPOINT ((1 1), (1.5 1.1), (2 1), (2 1.2), (2.1 1.1))")
    end
    @testset "EmptyIntersection" begin
        a = wkt("MULTIPOINT ((1 1), (3 1))"); b = wkt("POINT (2 1)")
        @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, a, b), 0)
    end
    @testset "EmptyInputIntersection" begin
        a = wkt("MULTIPOINT ((1 1), (3 1))"); b = wkt("POINT EMPTY")
        @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, a, b), 0)
    end
    @testset "EmptyInputUnion" begin
        a = wkt("MULTIPOINT ((1 1), (3 1))"); b = wkt("POINT EMPTY")
        @test eqw(ov(Planar(), GO.OVERLAY_UNION, a, b), "MULTIPOINT ((1 1), (3 1))")
    end
    @testset "EmptyDifference" begin
        a = wkt("MULTIPOINT ((1 1), (3 1))"); b = wkt("MULTIPOINT ((1 1), (2 1), (3 1))")
        @test is_empty_dim(ov(Planar(), GO.OVERLAY_DIFFERENCE, a, b), 0)
    end
end

@testset "point merge and output order" begin
    #-- duplicates merge; the first occurrence supplies the output coordinate,
    #-- and output order is input order (A then the B points not in A)
    a = wkt("MULTIPOINT ((1 1), (2 1), (1 1))")
    b = wkt("MULTIPOINT ((2 1), (3 1))")
    r = ov(Planar(), GO.OVERLAY_UNION, a, b)
    @test coords(r) == [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]
    #-- and the merge is idempotent for the other ops
    @test coords(ov(Planar(), GO.OVERLAY_INTERSECTION, a, b)) == [(2.0, 1.0)]
    @test coords(ov(Planar(), GO.OVERLAY_DIFFERENCE, a, b)) == [(1.0, 1.0)]
    @test coords(ov(Planar(), GO.OVERLAY_SYMDIFFERENCE, a, b)) == [(1.0, 1.0), (3.0, 1.0)]
    #-- signed zeros identify (kernel-point normalization)
    @test npts(ov(Planar(), GO.OVERLAY_INTERSECTION,
                       GI.Point((-0.0, 0.0)), GI.Point((0.0, -0.0)))) == 1
end

# ---------------------------------------------------------------------------
# 2. Ported JTS OverlayNGMixedPointsTest (point × line, point × area)
# ---------------------------------------------------------------------------
#
# SKIPPED, with reasons:
#   testPointLineIntersectionPrec — its expectation ("result is empty because
#     Line is not rounded") is entirely an artifact of the fixed precision
#     model; at floating precision the point is simply not on the line, which
#     the "point off the line" case below covers.

@testset "JTS OverlayNGMixedPointsTest (mixed dimensions)" begin
    L = wkt("LINESTRING (1 1, 9 1)")
    P = wkt("POINT (5 1)")
    @testset "SimpleLineIntersection" begin
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, L, P), "POINT (5 1)")
    end
    @testset "LinePointInOutIntersection" begin
        b = wkt("MULTIPOINT ((5 1), (15 1))")
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, L, b), "POINT (5 1)")
    end
    @testset "SimpleLineUnion" begin
        @test eqw(ov(Planar(), GO.OVERLAY_UNION, L, P), "LINESTRING (1 1, 9 1)")
    end
    @testset "SimpleLineDifference" begin
        @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, L, P), "LINESTRING (1 1, 9 1)")
    end
    @testset "SimpleLineSymDifference" begin
        @test eqw(ov(Planar(), GO.OVERLAY_SYMDIFFERENCE, L, P), "LINESTRING (1 1, 9 1)")
    end
    @testset "LinePointSymDifference" begin
        r = ov(Planar(), GO.OVERLAY_SYMDIFFERENCE, L, wkt("POINT (15 1)"))
        @test GI.trait(r) isa GI.GeometryCollectionTrait
        @test eqw(r, "GEOMETRYCOLLECTION (POINT (15 1), LINESTRING (1 1, 9 1))")
    end
    @testset "PolygonInsideIntersection" begin
        a = wkt("POLYGON ((4 2, 6 2, 6 0, 4 0, 4 2))")
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, a, P), "POINT (5 1)")
    end
    @testset "PolygonDisjointIntersection" begin
        a = wkt("POLYGON ((4 2, 6 2, 6 0, 4 0, 4 2))")
        @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, a, wkt("POINT (15 1)")), 0)
    end
    @testset "PointEmptyLinestringUnion" begin
        r = ov(Planar(), GO.OVERLAY_UNION, wkt("LINESTRING EMPTY"), wkt("POINT (10 10)"))
        @test eqw(r, "POINT (10 10)")
    end
    @testset "LinestringEmptyPointUnion" begin
        r = ov(Planar(), GO.OVERLAY_UNION, wkt("LINESTRING (10 10, 20 20)"), wkt("POINT EMPTY"))
        @test eqw(r, "LINESTRING (10 10, 20 20)")
    end
end

@testset "mixed-dimension symmetry (point on either side)" begin
    L = wkt("LINESTRING (1 1, 9 1)")
    PG = wkt("POLYGON ((4 2, 6 2, 6 0, 4 0, 4 2))")
    Pin = wkt("POINT (5 1)")
    #-- INTERSECTION / UNION / SYMDIFFERENCE are symmetric in the operands
    for (nonpt, ewkt) in ((L, "LINESTRING (1 1, 9 1)"), (PG, "POLYGON ((4 2, 6 2, 6 0, 4 0, 4 2))"))
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, nonpt, Pin), "POINT (5 1)")
        @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, Pin, nonpt), "POINT (5 1)")
        for op in (GO.OVERLAY_UNION, GO.OVERLAY_SYMDIFFERENCE)
            @test eqw(ov(Planar(), op, nonpt, Pin), ewkt)
            @test eqw(ov(Planar(), op, Pin, nonpt), ewkt)
        end
        #-- DIFFERENCE is NOT symmetric: A \ points == A; points \ A drops covered points
        @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, nonpt, Pin), ewkt)
        @test is_empty_dim(ov(Planar(), GO.OVERLAY_DIFFERENCE, Pin, nonpt), 0)
    end
    #-- a point off the non-point input survives the difference
    Pout = wkt("POINT (15 15)")
    @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, Pout, L), "POINT (15 15)")
    @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, Pout, PG), "POINT (15 15)")
end

@testset "multi-component non-point input" begin
    #-- extractPolygons/extractLines: every non-empty component survives, and
    #-- points falling in any of them are absorbed
    MP = wkt("MULTIPOLYGON (((0 0, 2 0, 2 2, 0 2, 0 0)), ((5 5, 7 5, 7 7, 5 7, 5 5)))")
    pts = wkt("MULTIPOINT ((1 1), (6 6), (10 10))")
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, MP, pts), "MULTIPOINT ((1 1), (6 6))")
    r = ov(Planar(), GO.OVERLAY_UNION, MP, pts)
    @test GI.trait(r) isa GI.GeometryCollectionTrait
    @test eqw(r, "GEOMETRYCOLLECTION (POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0)), " *
                 "POLYGON ((5 5, 7 5, 7 7, 5 7, 5 5)), POINT (10 10))")

    ML = wkt("MULTILINESTRING ((0 0, 2 0), (5 0, 7 0))")
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, ML, wkt("MULTIPOINT ((1 0), (6 0), (3 0))")),
              "MULTIPOINT ((1 0), (6 0))")
end

@testset "boundary points count as covered" begin
    #-- JTS `hasLocation` tests only for EXTERIOR, so a point on the boundary of
    #-- an area (or an endpoint of a line) is covered
    PG = wkt("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, PG, wkt("POINT (0 2)")), "POINT (0 2)")
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, PG, wkt("POINT (0 0)")), "POINT (0 0)")
    L = wkt("LINESTRING (0 0, 4 0)")
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, L, wkt("POINT (0 0)")), "POINT (0 0)")
    #-- a hole's interior is exterior to the polygon
    PGH = wkt("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (3 3, 7 3, 7 7, 3 7, 3 3))")
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, PGH, wkt("POINT (5 5)")), 0)
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, PGH, wkt("POINT (3 5)")), "POINT (3 5)")
end

# ---------------------------------------------------------------------------
# 3. Result dimension, empty and degenerate inputs
# ---------------------------------------------------------------------------

@testset "OverlayUtil.resultDimension" begin
    for d0 in 0:2, d1 in 0:2
        @test GO._result_dimension(GO.OVERLAY_INTERSECTION, d0, d1) == min(d0, d1)
        @test GO._result_dimension(GO.OVERLAY_UNION, d0, d1) == max(d0, d1)
        @test GO._result_dimension(GO.OVERLAY_DIFFERENCE, d0, d1) == d0
        @test GO._result_dimension(GO.OVERLAY_SYMDIFFERENCE, d0, d1) == max(d0, d1)
    end
end

@testset "createEmptyResult dimensions" begin
    for (dim, T) in ((0, GI.MultiPointTrait), (1, GI.MultiLineStringTrait),
                     (2, GI.MultiPolygonTrait))
        g = GO._empty_geom(dim)
        @test GI.trait(g) isa T
        @test GI.npoint(g) == 0
        #-- `GI.isempty` falls back to `false` for GeoInterface's own wrappers,
        #-- which is what the engine emits, so the driver uses `_ov_isempty`
        @test GO._ov_isempty(g)
    end
    @test GO._ov_isempty(wkt("POINT EMPTY"))
    @test !GO._ov_isempty(GI.Point((1.0, 2.0)))
    @test !GO._ov_isempty(GI.MultiPoint([(1.0, 2.0)]))
end

@testset "both inputs empty" begin
    EP = wkt("POINT EMPTY")
    EL = wkt("LINESTRING EMPTY")
    EPG = wkt("POLYGON EMPTY")
    for op in OPS
        @test is_empty_dim(ov(Planar(), op, EP, EP), 0)
    end
    #-- empty result dimension follows the op rules
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, EP, EL), 0)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_UNION, EP, EL), 1)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_UNION, EP, EPG), 2)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_DIFFERENCE, EP, EPG), 0)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_DIFFERENCE, EPG, EP), 2)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_SYMDIFFERENCE, EP, EPG), 2)
end

@testset "one empty input" begin
    pts = wkt("MULTIPOINT ((1 1), (2 2))")
    PG = wkt("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
    #-- INTERSECTION with an empty operand is empty at the min dimension
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, pts, wkt("POLYGON EMPTY")), 0)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, PG, wkt("POINT EMPTY")), 0)
    #-- DIFFERENCE with an empty A is empty at dim(A)
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_DIFFERENCE, wkt("POINT EMPTY"), PG), 0)
    #-- DIFFERENCE with an empty B keeps A
    @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, pts, wkt("POLYGON EMPTY")),
              "MULTIPOINT ((1 1), (2 2))")
    @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, PG, wkt("POINT EMPTY")),
              "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
    #-- UNION / SYMDIFFERENCE keep the non-empty operand
    for op in (GO.OVERLAY_UNION, GO.OVERLAY_SYMDIFFERENCE)
        @test eqw(ov(Planar(), op, pts, wkt("POLYGON EMPTY")), "MULTIPOINT ((1 1), (2 2))")
        @test eqw(ov(Planar(), op, PG, wkt("POINT EMPTY")), "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
    end
end

@testset "unsupported inputs still error" begin
    PG = wkt("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
    gc = wkt("GEOMETRYCOLLECTION (POINT (1 1), LINESTRING (0 0, 1 1))")
    @test_throws ArgumentError ov(Planar(), GO.OVERLAY_INTERSECTION, gc, PG)
    @test_throws ArgumentError ov(Planar(), GO.OVERLAY_INTERSECTION, PG, gc)
end

# ---------------------------------------------------------------------------
# 4. Spherical
# ---------------------------------------------------------------------------

@testset "spherical point × point" begin
    a = GI.MultiPoint([(1.0, 1.0), (2.0, 1.0)])
    b = GI.Point((2.0, 1.0))
    @test eqw(ov(Spherical(), GO.OVERLAY_INTERSECTION, a, b), "POINT (2 1)")
    @test eqw(ov(Spherical(), GO.OVERLAY_UNION, a, b), "MULTIPOINT ((1 1), (2 1))")
    @test eqw(ov(Spherical(), GO.OVERLAY_DIFFERENCE, a, b), "POINT (1 1)")
    @test eqw(ov(Spherical(), GO.OVERLAY_SYMDIFFERENCE, a, b), "POINT (1 1)")
    @test is_empty_dim(ov(Spherical(), GO.OVERLAY_INTERSECTION,
                          GI.Point((1.0, 1.0)), GI.Point((3.0, 1.0))), 0)
end

@testset "spherical point identity is geometric, not lon/lat-literal" begin
    #-- ±180° is one point on the sphere (and two distinct planar coordinates);
    #-- the kernel-point key must merge them
    a = GI.MultiPoint([(180.0, 10.0), (0.0, 0.0)])
    b = GI.MultiPoint([(-180.0, 10.0)])
    @test npts(ov(Spherical(), GO.OVERLAY_INTERSECTION, a, b)) == 1
    @test npts(ov(Spherical(), GO.OVERLAY_UNION, a, b)) == 2
    #-- planar treats them as distinct
    @test npts(ov(Planar(), GO.OVERLAY_INTERSECTION, a, b)) == 0
    @test npts(ov(Planar(), GO.OVERLAY_UNION, a, b)) == 3
    #-- the north pole has an arbitrary longitude
    n1 = GI.Point((0.0, 90.0)); n2 = GI.Point((123.0, 90.0))
    @test npts(ov(Spherical(), GO.OVERLAY_INTERSECTION, n1, n2)) == 1
    @test npts(ov(Planar(), GO.OVERLAY_INTERSECTION, n1, n2)) == 0
end

@testset "spherical point × area uses the kernel locator" begin
    #-- the ring (0,0)-(90,0)-(90,60)-(0,60): its northern edge is a great-circle
    #-- arc bowing to ~67.8°N at lon 45, so (45, 62) is INSIDE on the sphere and
    #-- outside the planar lat/lon rectangle. A planar even-odd test on the
    #-- emitted coordinates would get this wrong.
    PG = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (90.0, 60.0), (0.0, 60.0), (0.0, 0.0)]])
    inside = GI.Point((45.0, 62.0))
    @test eqw(ov(Spherical(), GO.OVERLAY_INTERSECTION, PG, inside), "POINT (45 62)")
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, PG, inside), 0)
    #-- and the same point is dropped by the spherical DIFFERENCE
    @test is_empty_dim(ov(Spherical(), GO.OVERLAY_DIFFERENCE, inside, PG), 0)
    @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, inside, PG), "POINT (45 62)")
    #-- a clearly-outside point behaves the same on both manifolds
    outside = GI.Point((45.0, 80.0))
    for m in (Planar(), Spherical())
        @test is_empty_dim(ov(m, GO.OVERLAY_INTERSECTION, PG, outside), 0)
        @test eqw(ov(m, GO.OVERLAY_DIFFERENCE, outside, PG), "POINT (45 80)")
    end
    #-- UNION copies the area through and keeps only exterior points
    r = ov(Spherical(), GO.OVERLAY_UNION, PG, GI.MultiPoint([(45.0, 62.0), (45.0, 80.0)]))
    @test GI.trait(r) isa GI.GeometryCollectionTrait
    @test eqw(r, "GEOMETRYCOLLECTION (POLYGON ((0 0, 90 0, 90 60, 0 60, 0 0)), POINT (45 80))")
end

@testset "spherical point × line uses the kernel locator" begin
    #-- the arc (0,0)→(90,45) passes through ~(45, 35.26), NOT through the
    #-- planar midpoint (45, 22.5)
    L = GI.LineString([(0.0, 0.0), (90.0, 45.0)])
    planar_mid = GI.Point((45.0, 22.5))
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, L, planar_mid), "POINT (45 22.5)")
    @test is_empty_dim(ov(Spherical(), GO.OVERLAY_INTERSECTION, L, planar_mid), 0)
    #-- an endpoint is on the arc on both manifolds
    for m in (Planar(), Spherical())
        @test eqw(ov(m, GO.OVERLAY_INTERSECTION, L, GI.Point((0.0, 0.0))), "POINT (0 0)")
    end
    #-- a point on a constant-latitude planar segment is off the great circle
    L2 = GI.LineString([(0.0, 60.0), (90.0, 60.0)])
    @test is_empty_dim(ov(Spherical(), GO.OVERLAY_INTERSECTION, L2, GI.Point((45.0, 60.0))), 0)
    @test eqw(ov(Planar(), GO.OVERLAY_INTERSECTION, L2, GI.Point((45.0, 60.0))), "POINT (45 60)")
end

@testset "spherical empty and degenerate inputs" begin
    PG = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    EP = wkt("POINT EMPTY")
    @test is_empty_dim(ov(Spherical(), GO.OVERLAY_INTERSECTION, PG, EP), 0)
    @test eqw(ov(Spherical(), GO.OVERLAY_UNION, PG, EP),
              "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
    for op in OPS
        @test is_empty_dim(ov(Spherical(), op, EP, EP), 0)
    end
end
