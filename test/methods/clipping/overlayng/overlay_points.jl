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
include(joinpath(@__DIR__, "common.jl"))
import GeometryOps: Planar, Spherical, True

const OPS = OP_CODES

ov(m, op, a, b; kw...) = GO._overlay_ng(m, op, a, b; exact = EX, kw...)

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
# 1 & 2. Ported JTS OverlayNGPointsTest / OverlayNGMixedPointsTest
# ---------------------------------------------------------------------------
#
# Table-driven: each row is one JTS case, and the name travels as data so a
# failure still says which case it was. `want` is the expected WKT, compared
# topologically, or an `Int` meaning "an empty result of that dimension".
#
# SKIPPED, with reasons (fixed-precision cases that cannot apply to a
# floating-only, non-snapping engine — design §0):
#   testDisjointPointsRoundedIntersection — POINT(10.1 10) ∩ POINT(10 10.1) is
#     only nonempty after rounding both to the unit grid.
#   testPointLineIntersectionPrec — its expectation ("result is empty because
#     Line is not rounded") is entirely an artifact of the fixed precision
#     model; at floating precision the point is simply not on the line, which
#     `mixed-dimension symmetry` covers.
# testSimpleMergeIntersection is kept, restated at floating precision (the
# extra A points simply do not participate).
#
# COVERED ELSEWHERE (ported, then found to be bit-identical calls to a case that
# is still here — deleted rather than duplicated, with the survivor named):
#   SimpleUnion             -> EmptyInputUnion, same op, same merge path
#   EmptyInputIntersection  -> EmptyIntersection (an empty B and a non-matching
#                              B take the same `_ov_isempty` short circuit)
#   SimpleLineIntersection  -> LinePointInOutIntersection, which is the same
#                              call with one extra off-line point
#   SimpleLineUnion / SimpleLineDifference / SimpleLineSymDifference
#                           -> `mixed-dimension symmetry`, which runs all three
#                              ops on the same line/point pair in both orders
#   PolygonInsideIntersection -> `mixed-dimension symmetry` (PG, Pin) and
#                              `boundary points count as covered`
#   PolygonDisjointIntersection -> never reached `_find_points` at all: the
#                              point is envelope-disjoint from the polygon, so
#                              `_empty_result_short_circuit` answers first

const POINT_CASES = (
    #  name                          op                       A                                                  B                                              want
    ("SimpleIntersection",           GO.OVERLAY_INTERSECTION,  "MULTIPOINT ((1 1), (2 1))",                        "POINT (2 1)",                                 "POINT (2 1)"),
    ("SimpleMergeIntersection",      GO.OVERLAY_INTERSECTION,  "MULTIPOINT ((1 1), (1.5 1.1), (2 1), (2.1 1.1))",  "POINT (2 1)",                                 "POINT (2 1)"),
    ("SimpleDifference",             GO.OVERLAY_DIFFERENCE,    "MULTIPOINT ((1 1), (2 1))",                        "POINT (2 1)",                                 "POINT (1 1)"),
    ("SimpleSymDifference",          GO.OVERLAY_SYMDIFFERENCE, "MULTIPOINT ((1 2), (1 1), (2 2), (2 1))",          "MULTIPOINT ((2 2), (2 1), (3 2), (3 1))",     "MULTIPOINT ((1 2), (1 1), (3 2), (3 1))"),
    ("SimpleFloatUnion",             GO.OVERLAY_UNION,         "MULTIPOINT ((1 1), (1.5 1.1), (2 1), (2.1 1.1))",  "MULTIPOINT ((1.5 1.1), (2 1), (2 1.2))",      "MULTIPOINT ((1 1), (1.5 1.1), (2 1), (2 1.2), (2.1 1.1))"),
    ("EmptyIntersection",            GO.OVERLAY_INTERSECTION,  "MULTIPOINT ((1 1), (3 1))",                        "POINT (2 1)",                                 0),
    ("EmptyInputUnion",              GO.OVERLAY_UNION,         "MULTIPOINT ((1 1), (3 1))",                        "POINT EMPTY",                                 "MULTIPOINT ((1 1), (3 1))"),
    ("EmptyDifference",              GO.OVERLAY_DIFFERENCE,    "MULTIPOINT ((1 1), (3 1))",                        "MULTIPOINT ((1 1), (2 1), (3 1))",            0),
    ("LinePointInOutIntersection",   GO.OVERLAY_INTERSECTION,  "LINESTRING (1 1, 9 1)",                            "MULTIPOINT ((5 1), (15 1))",                  "POINT (5 1)"),
    ("PointEmptyLinestringUnion",    GO.OVERLAY_UNION,         "LINESTRING EMPTY",                                 "POINT (10 10)",                               "POINT (10 10)"),
    ("LinestringEmptyPointUnion",    GO.OVERLAY_UNION,         "LINESTRING (10 10, 20 20)",                        "POINT EMPTY",                                 "LINESTRING (10 10, 20 20)"),
)

@testset "JTS OverlayNGPointsTest / OverlayNGMixedPointsTest" begin
    for (name, op, a, b, want) in POINT_CASES
        r = ov(Planar(), op, wkt(a), wkt(b))
        ok = want isa Integer ? is_empty_dim(r, want) : eqw(r, want)
        @test (name, ok) == (name, true)
    end

    #-- kept out of the table: it asserts the result's TYPE as well as its points
    @testset "LinePointSymDifference" begin
        r = ov(Planar(), GO.OVERLAY_SYMDIFFERENCE, wkt("LINESTRING (1 1, 9 1)"), wkt("POINT (15 1)"))
        @test GI.trait(r) isa GI.GeometryCollectionTrait
        @test eqw(r, "GEOMETRYCOLLECTION (POINT (15 1), LINESTRING (1 1, 9 1))")
    end
end

@testset "point merge and output order" begin
    #-- duplicates merge; the first occurrence supplies the output coordinate,
    #-- and output order is input order (A then the B points not in A). The
    #-- planar `(1 1)` duplicate cannot show first-vs-last provenance — both
    #-- occurrences have the same coordinate — so the spherical asserts below,
    #-- where two different lon/lat spellings name one point, are what pin it.
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
    #=
    FIRST occurrence supplies the output coordinate. The fixture is the pole,
    whose longitude is arbitrary: `(0, 90)` and `(45, 90)` are two spellings of
    ONE kernel point, so the merge really does have to choose, and which
    spelling survives is the choice made visible.

    It is only visible in the lon/lat row. The default row emits the kernel
    point itself, so both spellings produce the identical coordinate and the
    provenance question dissolves — which is asserted below rather than skipped,
    because "the output does not depend on which duplicate came first" is a
    stronger property than the rule it replaces, and it is the reason the rule
    stops being observable.
    =#
    ll = Tuple{Float64, Float64}
    north, north2 = (0.0, 90.0), (45.0, 90.0)
    ab = (GI.MultiPoint([north, (2.0, 1.0)]), GI.MultiPoint([north2]))
    ba = (GI.MultiPoint([north2, (2.0, 1.0)]), GI.MultiPoint([north]))
    r_ab = ov(Spherical(), GO.OVERLAY_UNION, ab...; point_type = ll)
    r_ba = ov(Spherical(), GO.OVERLAY_UNION, ba...; point_type = ll)
    #-- fixture premise: the two spellings merged, so one of them was dropped
    @test npts(r_ab) == 2 && npts(r_ba) == 2
    @test coords(r_ab)[1] == north
    @test coords(r_ba)[1] == north2
    #-- and in the default (xyz) row the two orders agree exactly
    x_ab = collect(GI.getpoint(ov(Spherical(), GO.OVERLAY_UNION, ab...)))
    x_ba = collect(GI.getpoint(ov(Spherical(), GO.OVERLAY_UNION, ba...)))
    @test x_ab[1] === x_ba[1] === GO._to_kernel_point(Spherical(), north)
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

    #-- SHAPE, not just point set: one surviving component must come back as the
    #-- atomic geometry and several as the `Multi` form. GEOS `equals` cannot see
    #-- the difference — `POINT (2 1)` and `MULTIPOINT ((2 1))` are equal to it —
    #-- so every assertion above passes against a result that always emits the
    #-- `Multi` form, and `_dimensional_result`'s whole job goes unchecked.
    @test GI.trait(ov(Planar(), GO.OVERLAY_INTERSECTION, MP, wkt("POINT (1 1)"))) isa GI.PointTrait
    @test GI.trait(ov(Planar(), GO.OVERLAY_UNION, MP, wkt("MULTIPOINT ((1 1), (6 6))"))) isa
          GI.MultiPolygonTrait

    #-- an EMPTY member inside a MULTILINESTRING / MULTIPOINT is skipped rather
    #-- than treated as a component at the origin
    ML = wkt("MULTILINESTRING ((0 0, 2 0), EMPTY, (5 0, 7 0))")
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
    #-- `_empty_geom`'s trait/point-count at each dimension is asserted through the
    #-- driver by "both inputs empty" below; what only shows here is that
    #-- `GI.isempty` falls back to `false` for GeoInterface's own wrappers — which
    #-- is what the engine emits — and that is why `_ov_isempty` exists
    for dim in (0, 1, 2), P in (Tuple{Float64, Float64}, USP)
        @test GO._ov_isempty(GO._empty_geom(P, dim))
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

#=
The engine in its lon/lat row. Everything asserted with `eqw` below is a claim
about WHICH components survive, stated against a WKT written in degrees — and
GEOS `equals` is an exact predicate, so it answers "not equal" about two
geometries that differ by one ULP in one vertex.

`OverlayNG(Spherical())` emits unit-sphere xyz by default, and bringing that
back to degrees for GEOS costs 2–4 ULPs of longitude on these fixtures
(`(1, 1)` returns as `(0.9999999999999998, 1)`). Comparing through that
conversion would mean loosening every assertion here from "this point" to "this
point, roughly" — trading an exact oracle for an approximate one to accommodate
a conversion the assertion is not even about.

So the survivorship claims run in the row whose output chart is the WKT's:
there an input vertex no operation cuts passes through bit-for-bit and both
sides of the comparison are exact. `point_type` is spelled for `Planar` too,
where it is the only supported type and therefore a no-op, so the mixed
manifold loops below read the same on both legs.

The default row is pinned elsewhere and not by weaker means: by kernel-point
identity in "point merge and output order" above, by the bit-exact vertex
round-trip in `api.jl`, and — on the areal paths, against s2geography — by
`s2_differential.jl`.
=#
llov(m, op, a, b) = ov(m, op, a, b; point_type = Tuple{Float64, Float64})

@testset "spherical point × point" begin
    a = GI.MultiPoint([(1.0, 1.0), (2.0, 1.0)])
    b = GI.Point((2.0, 1.0))
    @test eqw(llov(Spherical(), GO.OVERLAY_INTERSECTION, a, b), "POINT (2 1)")
    @test eqw(llov(Spherical(), GO.OVERLAY_UNION, a, b), "MULTIPOINT ((1 1), (2 1))")
    @test eqw(llov(Spherical(), GO.OVERLAY_DIFFERENCE, a, b), "POINT (1 1)")
    @test eqw(llov(Spherical(), GO.OVERLAY_SYMDIFFERENCE, a, b), "POINT (1 1)")
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
    #-- planar treats them as distinct. Asserted through UNION, not INTERSECTION:
    #-- these fixtures are envelope-disjoint on the plane, so a planar `∩` never
    #-- reaches the point path at all — `_empty_result_short_circuit` answers
    #-- first and the assertion holds however the locator behaves.
    @test npts(ov(Planar(), GO.OVERLAY_UNION, a, b)) == 3
    #-- the north pole has an arbitrary longitude
    n1 = GI.Point((0.0, 90.0)); n2 = GI.Point((123.0, 90.0))
    @test npts(ov(Spherical(), GO.OVERLAY_INTERSECTION, n1, n2)) == 1
    @test npts(ov(Planar(), GO.OVERLAY_UNION, n1, n2)) == 2
end

@testset "spherical point × area uses the kernel locator" begin
    #-- the ring (0,0)-(90,0)-(90,60)-(0,60): its northern edge is a great-circle
    #-- arc bowing to ~67.8°N at lon 45, so (45, 62) is INSIDE on the sphere and
    #-- outside the planar lat/lon rectangle. The planar contrast is pinned by the
    #-- DIFFERENCE leg below: DIFFERENCE has no disjoint-envelope short circuit,
    #-- so it runs the planar locator; the planar INTERSECTION here is answered by
    #-- that short circuit, since (45, 62) is outside the ring's lat/lon box.
    PG = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (90.0, 60.0), (0.0, 60.0), (0.0, 0.0)]])
    inside = GI.Point((45.0, 62.0))
    @test eqw(llov(Spherical(), GO.OVERLAY_INTERSECTION, PG, inside), "POINT (45 62)")
    @test is_empty_dim(ov(Planar(), GO.OVERLAY_INTERSECTION, PG, inside), 0)
    #-- and the same point is dropped by the spherical DIFFERENCE
    @test is_empty_dim(ov(Spherical(), GO.OVERLAY_DIFFERENCE, inside, PG), 0)
    @test eqw(ov(Planar(), GO.OVERLAY_DIFFERENCE, inside, PG), "POINT (45 62)")
    #-- a clearly-outside point behaves the same on both manifolds
    outside = GI.Point((45.0, 80.0))
    for m in (Planar(), Spherical())
        @test is_empty_dim(ov(m, GO.OVERLAY_INTERSECTION, PG, outside), 0)
        @test eqw(llov(m, GO.OVERLAY_DIFFERENCE, outside, PG), "POINT (45 80)")
    end
    #-- UNION copies the area through and keeps only exterior points
    r = llov(Spherical(), GO.OVERLAY_UNION, PG, GI.MultiPoint([(45.0, 62.0), (45.0, 80.0)]))
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
        @test eqw(llov(m, GO.OVERLAY_INTERSECTION, L, GI.Point((0.0, 0.0))), "POINT (0 0)")
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
    @test eqw(llov(Spherical(), GO.OVERLAY_UNION, PG, EP),
              "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
    for op in OPS
        @test is_empty_dim(ov(Spherical(), op, EP, EP), 0)
    end
end
