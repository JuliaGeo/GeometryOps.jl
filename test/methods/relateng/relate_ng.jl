# Port of JTS RelateNGTest.java (with the helper shape of
# RelateNGTestCase.java). Every test method is ported, in the same order as
# the Java file; the Java's commented-out checks are kept commented for
# parity.
#
# Prepared-mode checks (`check_prepared`/`check_prepared_matches`) follow
# RelateNGTestCase.checkPrepared/checkPreparedMatches: every result through
# the `PreparedRelate` path must equal the unprepared result. In addition
# (GO-side, Task 22) the full-matrix fixture pairs are recorded as they run
# and a wholesale prepared-vs-unprepared loop over a representative sample
# of them runs at the end, plus a cache-reuse smoke test.

using Test
import GeometryOps as GO
import GeoInterface as GI

include(joinpath(@__DIR__, "wkt_util.jl"))

const RUN_PREPARED = true

# =========================================================================
# RelateNGTestCase.java helpers
# =========================================================================

# Every full-matrix fixture pair is recorded as it runs; the wholesale
# prepared-mode loop at the end of the file samples from this list.
const PREPARED_FIXTURES = Tuple{String, String}[]

function check_relate(awkt, bwkt, expected_im::String)
    push!(PREPARED_FIXTURES, (awkt, bwkt))
    a, b = from_wkt(awkt), from_wkt(bwkt)
    @test string(GO.relate(GO.RelateNG(), a, b)) == expected_im
end

function check_relate_matches(awkt, bwkt, pattern::String, expected::Bool)
    check_predicate(() -> GO.pred_matches(pattern), awkt, bwkt, expected)
end

function check_predicate(pred_factory, awkt, bwkt, expected::Bool)
    a, b = from_wkt(awkt), from_wkt(bwkt)
    @test GO.relate_predicate(GO.RelateNG(), pred_factory(), a, b) == expected
end

function check_intersects_disjoint(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_intersects, wkta, wktb, expected)
    check_predicate(GO.pred_intersects, wktb, wkta, expected)
    check_predicate(GO.pred_disjoint, wkta, wktb, !expected)
    check_predicate(GO.pred_disjoint, wktb, wkta, !expected)
end

function check_contains_within(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_contains, wkta, wktb, expected)
    check_predicate(GO.pred_within, wktb, wkta, expected)
end

function check_covers_coveredby(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_covers, wkta, wktb, expected)
    check_predicate(GO.pred_coveredby, wktb, wkta, expected)
end

function check_crosses(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_crosses, wkta, wktb, expected)
    check_predicate(GO.pred_crosses, wktb, wkta, expected)
end

function check_overlaps(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_overlaps, wkta, wktb, expected)
    check_predicate(GO.pred_overlaps, wktb, wkta, expected)
end

function check_touches(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_touches, wkta, wktb, expected)
    check_predicate(GO.pred_touches, wktb, wkta, expected)
end

function check_equals(wkta, wktb, expected::Bool)
    check_predicate(GO.pred_equalstopo, wkta, wktb, expected)
    check_predicate(GO.pred_equalstopo, wktb, wkta, expected)
end

# The named predicates checked by RelateNGTestCase.checkPrepared, in the
# same order as the Java method.
const PREPARED_PREDICATES = [
    ("equalsTopo", GO.pred_equalstopo),
    ("intersects", GO.pred_intersects),
    ("disjoint", GO.pred_disjoint),
    ("covers", GO.pred_covers),
    ("coveredBy", GO.pred_coveredby),
    ("within", GO.pred_within),
    ("contains", GO.pred_contains),
    ("crosses", GO.pred_crosses),
    ("touches", GO.pred_touches),
]

function check_prepared(wkta, wktb)
    if !RUN_PREPARED
        @test_skip "prepared mode (Task 22)"
        return
    end
    a, b = from_wkt(wkta), from_wkt(wktb)
    prep_a = GO.prepare(GO.RelateNG(), a)
    for (name, pred_factory) in PREPARED_PREDICATES
        @test GO.relate_predicate(prep_a, pred_factory(), b) ==
            GO.relate_predicate(GO.RelateNG(), pred_factory(), a, b)
    end
    @test string(GO.relate(prep_a, b)) == string(GO.relate(GO.RelateNG(), a, b))
end

function check_prepared_matches(wkta, wktb, pattern::String)
    if !RUN_PREPARED
        @test_skip "prepared mode (Task 22)"
        return
    end
    a, b = from_wkt(wkta), from_wkt(wktb)
    prep_a = GO.prepare(GO.RelateNG(), a)
    @test GO.relate(prep_a, b, pattern) == GO.relate(GO.RelateNG(), a, b, pattern)
end

# The empty-geometry WKT list (RelateNGTest.java `empties`), hoisted above
# the testset wrapper (`const` is not allowed inside a testset block).
const EMPTIES = [
    "POINT EMPTY",
    "LINESTRING EMPTY",
    "POLYGON EMPTY",
    "MULTIPOINT EMPTY",
    "MULTILINESTRING EMPTY",
    "MULTIPOLYGON EMPTY",
    "GEOMETRYCOLLECTION EMPTY",
]

# =========================================================================
# RelateNGTest.java test methods
# =========================================================================

@testset "RelateNGTest" begin

@testset "testPointsDisjoint" begin
    a = "POINT (0 0)"
    b = "POINT (1 1)"
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
    check_equals(a, b, false)
    check_relate(a, b, "FF0FFF0F2")
end

# ======= P/P  ============

@testset "testPointsContained" begin
    a = "MULTIPOINT (0 0, 1 1, 2 2)"
    b = "MULTIPOINT (1 1, 2 2)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_equals(a, b, false)
    check_relate(a, b, "0F0FFFFF2")
end

@testset "testPointsEqual" begin
    a = "MULTIPOINT (0 0, 1 1, 2 2)"
    b = "MULTIPOINT (0 0, 1 1, 2 2)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_equals(a, b, true)
end

@testset "testValidateRelatePP_13" begin
    a = "MULTIPOINT ((80 70), (140 120), (20 20), (200 170))"
    b = "MULTIPOINT ((80 70), (140 120), (80 170), (200 80))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_contains_within(b, a, false)
    check_covers_coveredby(a, b, false)
    check_overlaps(a, b, true)
    check_touches(a, b, false)
end

# ======= L/P  ============

@testset "testLinePointContains" begin
    a = "LINESTRING (0 0, 1 1, 2 2)"
    b = "MULTIPOINT (0 0, 1 1, 2 2)"
    check_relate(a, b, "0F10FFFF2")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_contains_within(b, a, false)
    check_covers_coveredby(a, b, true)
    check_covers_coveredby(b, a, false)
end

@testset "testLinePointOverlaps" begin
    a = "LINESTRING (0 0, 1 1)"
    b = "MULTIPOINT (0 0, 1 1, 2 2)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_contains_within(b, a, false)
    check_covers_coveredby(a, b, false)
    check_covers_coveredby(b, a, false)
end

@testset "testZeroLengthLinePoint" begin
    a = "LINESTRING (0 0, 0 0)"
    b = "POINT (0 0)"
    check_relate(a, b, "0FFFFFFF2")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_contains_within(b, a, true)
    check_covers_coveredby(a, b, true)
    check_covers_coveredby(b, a, true)
    check_equals(a, b, true)
end

@testset "testZeroLengthLineLine" begin
    a = "LINESTRING (10 10, 10 10, 10 10)"
    b = "LINESTRING (10 10, 10 10)"
    check_relate(a, b, "0FFFFFFF2")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_contains_within(b, a, true)
    check_covers_coveredby(a, b, true)
    check_covers_coveredby(b, a, true)
    check_equals(a, b, true)
end

# tests bug involving checking for non-zero-length lines
@testset "testNonZeroLengthLinePoint" begin
    a = "LINESTRING (0 0, 0 0, 9 9)"
    b = "POINT (1 1)"
    check_relate(a, b, "0F1FF0FF2")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_contains_within(b, a, false)
    check_covers_coveredby(a, b, true)
    check_covers_coveredby(b, a, false)
    check_equals(a, b, false)
end

@testset "testLinePointIntAndExt" begin
    a = "MULTIPOINT((60 60), (100 100))"
    b = "LINESTRING(40 40, 80 80)"
    check_relate(a, b, "0F0FFF102")
end

# ======= L/L  ============

@testset "testLinesCrossProper" begin
    a = "LINESTRING (0 0, 9 9)"
    b = "LINESTRING(0 9, 9 0)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
end

@testset "testLinesOverlap" begin
    a = "LINESTRING (0 0, 5 5)"
    b = "LINESTRING(3 3, 9 9)"
    check_intersects_disjoint(a, b, true)
    check_touches(a, b, false)
    check_overlaps(a, b, true)
end

@testset "testLinesCrossVertex" begin
    a = "LINESTRING (0 0, 8 8)"
    b = "LINESTRING(0 8, 4 4, 8 0)"
    check_intersects_disjoint(a, b, true)
end

@testset "testLinesTouchVertex" begin
    a = "LINESTRING (0 0, 8 0)"
    b = "LINESTRING(0 8, 4 0, 8 8)"
    check_intersects_disjoint(a, b, true)
end

@testset "testLinesDisjointByEnvelope" begin
    a = "LINESTRING (0 0, 9 9)"
    b = "LINESTRING(10 19, 19 10)"
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
end

@testset "testLinesDisjoint" begin
    a = "LINESTRING (0 0, 9 9)"
    b = "LINESTRING (4 2, 8 6)"
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
end

@testset "testLinesClosedEmpty" begin
    a = "MULTILINESTRING ((0 0, 0 1), (0 1, 1 1, 1 0, 0 0))"
    b = "LINESTRING EMPTY"
    check_relate(a, b, "FF1FFFFF2")
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
end

@testset "testLinesRingTouchAtNode" begin
    a = "LINESTRING (5 5, 1 8, 1 1, 5 5)"
    b = "LINESTRING (5 5, 9 5)"
    check_relate(a, b, "F01FFF102")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_touches(a, b, true)
end

@testset "testLinesTouchAtBdy" begin
    a = "LINESTRING (5 5, 1 8)"
    b = "LINESTRING (5 5, 9 5)"
    check_relate(a, b, "FF1F00102")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_touches(a, b, true)
end

@testset "testLinesOverlapWithDisjointLine" begin
    a = "LINESTRING (1 1, 9 9)"
    b = "MULTILINESTRING ((2 2, 8 8), (6 2, 8 4))"
    check_relate(a, b, "101FF0102")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_overlaps(a, b, true)
end

@testset "testLinesDisjointOverlappingEnvelopes" begin
    a = "LINESTRING (60 0, 20 80, 100 80, 80 120, 40 140)"
    b = "LINESTRING (60 40, 140 40, 140 160, 0 160)"
    check_relate(a, b, "FF1FF0102")
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
    check_touches(a, b, false)
end

#=
Case from https://github.com/locationtech/jts/issues/270
Strictly, the lines cross, since their interiors intersect
according to the Orientation predicate.
However, the computation of the intersection point is
non-robust, and reports it as being equal to the endpoint
POINT (-10 0.0000000000000012)
For consistency the relate algorithm uses the intersection node topology.
=#
@testset "testLinesCross_JTS270" begin
    a = "LINESTRING (0 0, -10 0.0000000000000012)"
    b = "LINESTRING (-9.999143275740073 -0.1308959557133398, -10 0.0000000000001054)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_crosses(a, b, false)
    check_overlaps(a, b, false)
    check_touches(a, b, true)
end

@testset "testLinesContained_JTS396" begin
    a = "LINESTRING (1 0, 0 2, 0 0, 2 2)"
    b = "LINESTRING (0 0, 2 2)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_covers_coveredby(a, b, true)
    check_crosses(a, b, false)
    check_overlaps(a, b, false)
    check_touches(a, b, false)
end

#=
This case shows that lines must be self-noded,
so that node topology is constructed correctly
(at least for some predicates).
=#
@testset "testLinesContainedWithSelfIntersection" begin
    a = "LINESTRING (2 0, 0 2, 0 0, 2 2)"
    b = "LINESTRING (0 0, 2 2)"
    #check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_covers_coveredby(a, b, true)
    check_crosses(a, b, false)
    check_overlaps(a, b, false)
    check_touches(a, b, false)
end

@testset "testLineContainedInRing" begin
    a = "LINESTRING(60 60, 100 100, 140 60)"
    b = "LINESTRING(100 100, 180 20, 20 20, 100 100)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(b, a, true)
    check_covers_coveredby(b, a, true)
    check_crosses(a, b, false)
    check_overlaps(a, b, false)
    check_touches(a, b, false)
end

# see https://github.com/libgeos/geos/issues/933
@testset "testLineLineProperIntersection" begin
    a = "MULTILINESTRING ((0 0, 1 1), (0.5 0.5, 1 0.1, -1 0.1))"
    b = "LINESTRING (0 0, 1 1)"
    #check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_covers_coveredby(a, b, true)
    check_crosses(a, b, false)
    check_overlaps(a, b, false)
    check_touches(a, b, false)
end

@testset "testLineSelfIntersectionCollinear" begin
    a = "LINESTRING (9 6, 1 6, 1 0, 5 6, 9 6)"
    b = "LINESTRING (9 9, 3 1)"
    check_relate(a, b, "0F1FFF102")
end

# ======= A/P  ============

@testset "testPolygonPointInside" begin
    a = "POLYGON ((0 10, 10 10, 10 0, 0 0, 0 10))"
    b = "POINT (1 1)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
end

@testset "testPolygonPointOutside" begin
    a = "POLYGON ((10 0, 0 0, 0 10, 10 0))"
    b = "POINT (8 8)"
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
end

@testset "testPolygonPointInBoundary" begin
    a = "POLYGON ((10 0, 0 0, 0 10, 10 0))"
    b = "POINT (1 0)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, true)
end

@testset "testAreaPointInExterior" begin
    a = "POLYGON ((1 5, 5 5, 5 1, 1 1, 1 5))"
    b = "POINT (7 7)"
    check_relate(a, b, "FF2FF10F2")
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_touches(a, b, false)
    check_overlaps(a, b, false)
end

# ======= A/L  ============

@testset "testAreaLineContainedAtLineVertex" begin
    a = "POLYGON ((1 5, 5 5, 5 1, 1 1, 1 5))"
    b = "LINESTRING (2 3, 3 5, 4 3)"
    check_intersects_disjoint(a, b, true)
    #check_contains_within(a, b, true)
    #check_covers_coveredby(a, b, true)
    check_touches(a, b, false)
    check_overlaps(a, b, false)
end

@testset "testAreaLineTouchAtLineVertex" begin
    a = "POLYGON ((1 5, 5 5, 5 1, 1 1, 1 5))"
    b = "LINESTRING (1 8, 3 5, 5 8)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_touches(a, b, true)
    check_overlaps(a, b, false)
end

@testset "testPolygonLineInside" begin
    a = "POLYGON ((0 10, 10 10, 10 0, 0 0, 0 10))"
    b = "LINESTRING (1 8, 3 5, 5 8)"
    check_relate(a, b, "102FF1FF2")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
end

@testset "testPolygonLineOutside" begin
    a = "POLYGON ((10 0, 0 0, 0 10, 10 0))"
    b = "LINESTRING (4 8, 9 3)"
    check_intersects_disjoint(a, b, false)
    check_contains_within(a, b, false)
end

@testset "testPolygonLineInBoundary" begin
    a = "POLYGON ((10 0, 0 0, 0 10, 10 0))"
    b = "LINESTRING (1 0, 9 0)"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, true)
    check_touches(a, b, true)
    check_overlaps(a, b, false)
end

@testset "testPolygonLineCrossingContained" begin
    a = "MULTIPOLYGON (((20 80, 180 80, 100 0, 20 80)), ((20 160, 180 160, 100 80, 20 160)))"
    b = "LINESTRING (100 140, 100 40)"
    check_relate(a, b, "1020F1FF2")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_covers_coveredby(a, b, true)
    check_touches(a, b, false)
    check_overlaps(a, b, false)
end

@testset "testValidateRelateLA_220" begin
    a = "LINESTRING (90 210, 210 90)"
    b = "POLYGON ((150 150, 410 150, 280 20, 20 20, 150 150))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_touches(a, b, false)
    check_overlaps(a, b, false)
end

# See RelateLA.xml (line 585)
@testset "testLineCrossingPolygonAtShellHolePoint" begin
    a = "LINESTRING (60 160, 150 70)"
    b = "POLYGON ((190 190, 360 20, 20 20, 190 190), (110 110, 250 100, 140 30, 110 110))"
    check_relate(a, b, "F01FF0212")
    check_touches(a, b, true)
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_touches(a, b, true)
    check_overlaps(a, b, false)
end

@testset "testLineCrossingPolygonAtNonVertex" begin
    a = "LINESTRING (20 60, 150 60)"
    b = "POLYGON ((150 150, 410 150, 280 20, 20 20, 150 150))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_touches(a, b, false)
    check_overlaps(a, b, false)
end

@testset "testPolygonLinesContainedCollinearEdge" begin
    a = "POLYGON ((110 110, 200 20, 20 20, 110 110))"
    b = "MULTILINESTRING ((110 110, 60 40, 70 20, 150 20, 170 40), (180 30, 40 30, 110 80))"
    check_relate(a, b, "102101FF2")
end

# ======= A/A  ============

@testset "testPolygonsEdgeAdjacent" begin
    a = "POLYGON ((1 3, 3 3, 3 1, 1 1, 1 3))"
    b = "POLYGON ((5 3, 5 1, 3 1, 3 3, 5 3))"
    #check_intersects_disjoint(a, b, true)
    check_overlaps(a, b, false)
    check_touches(a, b, true)
    check_overlaps(a, b, false)
end

@testset "testPolygonsEdgeAdjacent2" begin
    a = "POLYGON ((1 3, 4 3, 3 0, 1 1, 1 3))"
    b = "POLYGON ((5 3, 5 1, 3 0, 4 3, 5 3))"
    #check_intersects_disjoint(a, b, true)
    check_overlaps(a, b, false)
    check_touches(a, b, true)
    check_overlaps(a, b, false)
end

@testset "testPolygonsNested" begin
    a = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))"
    b = "POLYGON ((2 8, 8 8, 8 2, 2 2, 2 8))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_covers_coveredby(a, b, true)
    check_overlaps(a, b, false)
    check_touches(a, b, false)
end

@testset "testPolygonsOverlapProper" begin
    a = "POLYGON ((1 1, 1 7, 7 7, 7 1, 1 1))"
    b = "POLYGON ((2 8, 8 8, 8 2, 2 2, 2 8))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_overlaps(a, b, true)
    check_touches(a, b, false)
end

@testset "testPolygonsOverlapAtNodes" begin
    a = "POLYGON ((1 5, 5 5, 5 1, 1 1, 1 5))"
    b = "POLYGON ((7 3, 5 1, 3 3, 5 5, 7 3))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_overlaps(a, b, true)
    check_touches(a, b, false)
end

@testset "testPolygonsContainedAtNodes" begin
    a = "POLYGON ((1 5, 5 5, 6 2, 1 1, 1 5))"
    b = "POLYGON ((1 1, 5 5, 6 2, 1 1))"
    #check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, true)
    check_covers_coveredby(a, b, true)
    check_overlaps(a, b, false)
    check_touches(a, b, false)
end

@testset "testPolygonsNestedWithHole" begin
    a = "POLYGON ((40 60, 420 60, 420 320, 40 320, 40 60), (200 140, 160 220, 260 200, 200 140))"
    b = "POLYGON ((80 100, 360 100, 360 280, 80 280, 80 100))"
    #check_intersects_disjoint(true, a, b)
    check_contains_within(a, b, false)
    check_contains_within(b, a, false)
    #check_covers_coveredby(false, a, b)
    #check_overlaps(true, a, b)
    check_predicate(GO.pred_contains, a, b, false)
    #check_touches(false, a, b)
end

@testset "testPolygonsOverlappingWithBoundaryInside" begin
    a = "POLYGON ((100 60, 140 100, 100 140, 60 100, 100 60))"
    b = "MULTIPOLYGON (((80 40, 120 40, 120 80, 80 80, 80 40)), ((120 80, 160 80, 160 120, 120 120, 120 80)), ((80 120, 120 120, 120 160, 80 160, 80 120)), ((40 80, 80 80, 80 120, 40 120, 40 80)))"
    check_relate(a, b, "21210F212")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_contains_within(b, a, false)
    check_covers_coveredby(a, b, false)
    check_overlaps(a, b, true)
    check_touches(a, b, false)
end

@testset "testPolygonsOverlapVeryNarrow" begin
    a = "POLYGON ((120 100, 120 200, 200 200, 200 100, 120 100))"
    b = "POLYGON ((100 100, 100000 110, 100000 100, 100 100))"
    check_relate(a, b, "212111212")
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_contains_within(b, a, false)
    #check_covers_coveredby(false, a, b)
    #check_overlaps(true, a, b)
    #check_touches(false, a, b)
end

@testset "testValidateRelateAA_86" begin
    a = "POLYGON ((170 120, 300 120, 250 70, 120 70, 170 120))"
    b = "POLYGON ((150 150, 410 150, 280 20, 20 20, 150 150), (170 120, 330 120, 260 50, 100 50, 170 120))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_overlaps(a, b, false)
    check_predicate(GO.pred_within, a, b, false)
    check_touches(a, b, true)
end

@testset "testValidateRelateAA_97" begin
    a = "POLYGON ((330 150, 200 110, 150 150, 280 190, 330 150))"
    b = "MULTIPOLYGON (((140 110, 260 110, 170 20, 50 20, 140 110)), ((300 270, 420 270, 340 190, 220 190, 300 270)))"
    check_intersects_disjoint(a, b, true)
    check_contains_within(a, b, false)
    check_covers_coveredby(a, b, false)
    check_overlaps(a, b, false)
    check_predicate(GO.pred_within, a, b, false)
    check_touches(a, b, true)
end

@testset "testAdjacentPolygons" begin
    a = "POLYGON ((1 9, 6 9, 6 1, 1 1, 1 9))"
    b = "POLYGON ((9 9, 9 4, 6 4, 6 9, 9 9))"
    check_relate_matches(a, b, GO.IM_PATTERN_ADJACENT, true)
end

@testset "testAdjacentPolygonsTouchingAtPoint" begin
    a = "POLYGON ((1 9, 6 9, 6 1, 1 1, 1 9))"
    b = "POLYGON ((9 9, 9 4, 6 4, 7 9, 9 9))"
    check_relate_matches(a, b, GO.IM_PATTERN_ADJACENT, false)
end

@testset "testAdjacentPolygonsOverlappping" begin
    a = "POLYGON ((1 9, 6 9, 6 1, 1 1, 1 9))"
    b = "POLYGON ((9 9, 9 4, 6 4, 5 9, 9 9))"
    check_relate_matches(a, b, GO.IM_PATTERN_ADJACENT, false)
end

@testset "testContainsProperlyPolygonContained" begin
    a = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))"
    b = "POLYGON ((2 8, 5 8, 5 5, 2 5, 2 8))"
    check_relate_matches(a, b, GO.IM_PATTERN_CONTAINS_PROPERLY, true)
end

@testset "testContainsProperlyPolygonTouching" begin
    a = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))"
    b = "POLYGON ((9 1, 5 1, 5 5, 9 5, 9 1))"
    check_relate_matches(a, b, GO.IM_PATTERN_CONTAINS_PROPERLY, false)
end

@testset "testContainsProperlyPolygonsOverlapping" begin
    a = "GEOMETRYCOLLECTION (POLYGON ((1 9, 6 9, 6 4, 1 4, 1 9)), POLYGON ((2 4, 6 7, 9 1, 2 4)))"
    b = "POLYGON ((5 5, 6 5, 6 4, 5 4, 5 5))"
    check_relate_matches(a, b, GO.IM_PATTERN_CONTAINS_PROPERLY, true)
end

# ================  Repeated Points  =============

@testset "testRepeatedPointLL" begin
    a = "LINESTRING(0 0, 5 5, 5 5, 5 5, 9 9)"
    b = "LINESTRING(0 9, 5 5, 5 5, 5 5, 9 0)"
    check_relate(a, b, "0F1FF0102")
    check_intersects_disjoint(a, b, true)
end

@testset "testRepeatedPointAA" begin
    a = "POLYGON ((1 9, 9 7, 9 1, 1 3, 1 9))"
    b = "POLYGON ((1 3, 1 3, 1 3, 3 7, 9 7, 9 7, 1 3))"
    check_relate(a, b, "212F01FF2")
end

# ================  EMPTY geometries  =============

@testset "testEmptyEmpty" begin
    for a in EMPTIES
        for b in EMPTIES
            check_relate(a, b, "FFFFFFFF2")
            #-- empty geometries are all topologically equal
            check_equals(a, b, true)

            check_intersects_disjoint(a, b, false)
            check_contains_within(a, b, false)
        end
    end
end

@testset "testEmptyNonEmpty" begin
    non_empty_point = "POINT (1 1)"
    non_empty_line = "LINESTRING (1 1, 2 2)"
    non_empty_polygon = "POLYGON ((1 1, 1 2, 2 1, 1 1))"

    for empty in EMPTIES
        check_relate(empty, non_empty_point, "FFFFFF0F2")
        check_relate(non_empty_point, empty, "FF0FFFFF2")

        check_relate(empty, non_empty_line, "FFFFFF102")
        check_relate(non_empty_line, empty, "FF1FF0FF2")

        check_relate(empty, non_empty_polygon, "FFFFFF212")
        check_relate(non_empty_polygon, empty, "FF2FF1FF2")

        check_equals(empty, non_empty_point, false)
        check_equals(empty, non_empty_line, false)
        check_equals(empty, non_empty_polygon, false)

        check_intersects_disjoint(empty, non_empty_point, false)
        check_intersects_disjoint(empty, non_empty_line, false)
        check_intersects_disjoint(empty, non_empty_polygon, false)

        check_contains_within(empty, non_empty_point, false)
        check_contains_within(empty, non_empty_line, false)
        check_contains_within(empty, non_empty_polygon, false)

        check_contains_within(non_empty_point, empty, false)
        check_contains_within(non_empty_line, empty, false)
        check_contains_within(non_empty_polygon, empty, false)
    end
end

# ================  Prepared Relate  =============

@testset "testPreparedAA" begin
    a = "POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))"
    b = "POLYGON((0.5 0.5, 1.5 0.5, 1.5 1.5, 0.5 1.5, 0.5 0.5))"
    check_prepared(a, b)
end

@testset "testPreparedPA" begin
    a = "POINT (5 5)"
    b = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))"
    check_prepared(a, b)
    check_prepared(b, a)

    #-- see https://github.com/libgeos/geos/issues/1275 (not a bug, but a good test to have)
    pattern = "T*****FF*"
    pattern_trans = "T*F**F***"   # IntersectionMatrix.transpose(pattern)
    check_prepared_matches(a, b, pattern)
    check_prepared_matches(b, a, pattern_trans)
end

# ===  Prepared mode, GO-side additions (not in RelateNGTest.java)  ===

# Wholesale prepared-vs-unprepared equality over the full-matrix fixture
# pairs recorded by `check_relate` above (they span P/L/A x P/L/A and the
# empty/GC cases). An evenly spaced sample capped at 40 pairs — plus every
# GeometryCollection pair, and one explicit mixed GC since RelateNGTest's
# only GC matrix fixtures are empty — keeps the added runtime modest
# (`check_prepared` is ~20 engine evaluations per pair).
@testset "testPreparedWholesale" begin
    fixtures = unique(PREPARED_FIXTURES)
    @test length(fixtures) >= 30
    idxs = unique(round.(Int, range(1, length(fixtures); length = min(40, length(fixtures)))))
    sample = fixtures[idxs]
    for fix in fixtures
        is_gc = occursin("GEOMETRYCOLLECTION", fix[1]) || occursin("GEOMETRYCOLLECTION", fix[2])
        is_gc && !(fix in sample) && push!(sample, fix)
    end
    push!(sample, (
        "GEOMETRYCOLLECTION (POINT (1 1), LINESTRING (0 5, 5 5), POLYGON ((0 0, 0 3, 3 3, 3 0, 0 0)))",
        "POLYGON ((2 2, 2 6, 6 6, 6 2, 2 2))",
    ))
    for (awkt, bwkt) in sample
        check_prepared(awkt, bwkt)
    end
end

# Cache-reuse smoke test: one `PreparedRelate` evaluated against several
# different B geometries, each twice, exercising the prebuilt-tree path
# (large A), the below-threshold nested-loop path (small A), and the
# self-noding path that bypasses the cached edges (line A).
@testset "testPreparedCacheReuse" begin
    #-- large A (64 segments >= threshold): the segment tree is prebuilt
    n = 64
    coords = [(5 + 4 * cospi(2k / n), 5 + 4 * sinpi(2k / n)) for k in 0:(n - 1)]
    push!(coords, coords[1])
    a = GI.Polygon([coords])
    prep = GO.prepare(GO.RelateNG(), a)
    @test prep.edge_tree !== nothing
    bs = [
        GI.Polygon([[(4.0, 4.0), (11.0, 4.0), (11.0, 11.0), (4.0, 11.0), (4.0, 4.0)]]),
        GI.LineString([(-1.0, 5.0), (11.0, 5.0)]),
        GI.Point((20.0, 20.0)),
    ]
    for b in bs
        expected = string(GO.relate(GO.RelateNG(), a, b))
        @test string(GO.relate(prep, b)) == expected
        #-- second evaluation against the same prepared instance
        @test string(GO.relate(prep, b)) == expected
        @test GO.relate_predicate(prep, GO.pred_intersects(), b) ==
            GO.relate_predicate(GO.RelateNG(), GO.pred_intersects(), a, b)
    end

    #-- small A: below the threshold no tree is prebuilt (nested-loop reuse)
    a2 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
    prep2 = GO.prepare(GO.RelateNG(), a2)
    @test prep2.edge_tree === nothing
    b21 = GI.Polygon([[(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5)]])
    b22 = GI.LineString([(0.5, -1.0), (0.5, 2.0)])
    @test string(GO.relate(prep2, b21)) == string(GO.relate(GO.RelateNG(), a2, b21))
    @test string(GO.relate(prep2, b22)) == string(GO.relate(GO.RelateNG(), a2, b22))

    #-- self-crossing line A: matrix evaluation requires self-noding, which
    #-- bypasses the cached edges (re-extracted per call, envelope-filtered)
    a3 = GI.LineString([(0.0, 0.0), (5.0, 5.0), (5.0, 0.0), (0.0, 5.0)])
    prep3 = GO.prepare(GO.RelateNG(), a3)
    b3 = GI.LineString([(0.0, 2.0), (6.0, 2.0)])
    expected3 = string(GO.relate(GO.RelateNG(), a3, b3))
    @test string(GO.relate(prep3, b3)) == expected3
    @test string(GO.relate(prep3, b3)) == expected3
end

end # @testset "RelateNGTest"

# =========================================================================
# Public API surface (Task 26): default-algorithm `relate`, exported names,
# and the named-predicate methods `GO.intersects(GO.RelateNG(), a, b)` etc.
# The named predicates are opt-in (design D4): the two-argument defaults
# dispatch to the old engines and are untouched here.
# =========================================================================

@testset "PublicAPI" begin
    pa = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]])
    pb = GI.Polygon([[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0), (2.0, 2.0)]])

    @testset "exports" begin
        for name in (:relate, :DE9IM, :RelateNG, :prepare,
                :Mod2Boundary, :EndpointBoundary,
                :MultivalentEndpointBoundary, :MonovalentEndpointBoundary)
            @test name in names(GO)
        end
    end

    @testset "relate entry points" begin
        @test GO.relate(pa, pb) isa GO.DE9IM
        @test string(GO.relate(pa, pb)) == "212101212"
        @test GO.relate(pa, pb, "T*F**FFF*") isa Bool
        @test GO.relate(pa, pb) == GO.relate(GO.RelateNG(), pa, pb)
        #-- prepared equivalents
        prep = GO.prepare(GO.RelateNG(), pa)
        @test GO.relate(prep, pb) == GO.relate(pa, pb)
        @test GO.relate(prep, pb, "T********") == GO.relate(pa, pb, "T********") == true
    end

    @testset "named predicates (opt-in RelateNG form)" begin
        alg = GO.RelateNG()
        far = GI.Polygon([[(10.0, 10.0), (12.0, 10.0), (12.0, 12.0), (10.0, 12.0), (10.0, 10.0)]])
        inner = GI.Polygon([[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]])
        touching = GI.Polygon([[(4.0, 0.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0), (4.0, 0.0)]])

        #-- each named predicate must agree with the (old-engine) default
        for f in (GO.intersects, GO.disjoint, GO.contains, GO.within,
                GO.covers, GO.coveredby, GO.touches, GO.equals)
            for (x, y) in [(pa, pb), (pa, far), (pa, inner), (inner, pa),
                    (pa, touching), (pa, pa)]
                @test f(alg, x, y) == f(x, y)
            end
        end
        #-- `overlaps`: skip the edge-touching pair, where the old engine
        #-- reports edge intersection (DE-9IM `overlaps` of touching
        #-- polygons is false; the old GO method returns true there)
        for (x, y) in [(pa, pb), (pa, far), (pa, inner), (pa, pa)]
            @test GO.overlaps(alg, x, y) == GO.overlaps(x, y)
        end
        # LibGEOS agrees; old GO.overlaps wrongly returns true (known old-GO gap)
        @test GO.overlaps(GO.RelateNG(), pa, touching) == false
        #-- `crosses`: the old GO method only supports mixed-dimension
        #-- (point/line/polygon) pairs
        crossing_line = GI.LineString([(-1.0, 2.0), (5.0, 2.0)])
        miss_line = GI.LineString([(-2.0, -2.0), (-1.0, -2.0)])
        @test GO.crosses(alg, crossing_line, pa) == GO.crosses(crossing_line, pa) == true
        @test GO.crosses(alg, miss_line, pa) == GO.crosses(miss_line, pa) == false
        #-- swapped argument order (old GO supports poly/line too)
        @test GO.crosses(alg, pa, crossing_line) == GO.crosses(pa, crossing_line) == true
        @test GO.crosses(alg, pa, miss_line) == GO.crosses(pa, miss_line) == false
        #-- poly/poly crosses is always false per DE-9IM dimension rules;
        #-- RelateNG-only capability (old GO has no poly/poly `crosses` method)
        @test GO.crosses(GO.RelateNG(), pa, pb) == false

        #-- spot values
        @test GO.intersects(alg, pa, pb)
        @test GO.disjoint(alg, pa, far)
        @test GO.contains(alg, pa, inner)
        @test GO.within(alg, inner, pa)
        @test GO.covers(alg, pa, inner)
        @test GO.coveredby(alg, inner, pa)
        @test GO.touches(alg, pa, touching)
        @test GO.overlaps(alg, pa, pb)
        #-- `equals(::RelateNG, ...)` is *topological* equality
        #-- (pred_equalstopo): a rotated ring is still equal
        pa_rot = GI.Polygon([[(4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0), (4.0, 0.0)]])
        @test GO.equals(alg, pa, pa_rot) == GO.equals(pa, pa_rot) == true
        @test !GO.equals(alg, pa, pb)
    end
end
