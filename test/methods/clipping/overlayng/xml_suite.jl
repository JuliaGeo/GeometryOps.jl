# JTS overlay XML conformance (phase-3 validation, design §4).
#
# Two runs over the vendored JTS corpora, both driven through the internal
# `_overlay_ng` engine on `Planar()`:
#
# 1. `test/data/jts/overlay/` — the FLOATING-precision general + misc overlay
#    suites. Every `intersection` / `union` / `difference` / `symdifference` op
#    (and its `*NG` spelling) is run and compared to the XML's expected geometry
#    with GEOS topological `equals`. Coordinate equality is deliberately not
#    used: this engine computes an exact arrangement and rounds once at
#    emission, so its vertices need not be bit-identical to JTS's.
# 2. `test/data/jts/overlay_robust/` — JTS's `robust/overlay` corpus of overlay
#    regressions harvested from GEOS / PostGIS / QGIS / shapely / JTS bug
#    reports. Those files mostly carry no expected geometry (they are hard
#    *inputs*), so they are run through JTS's own oracle-free identity check,
#    `overlayAreaTest` / `TestCaseGeometryFunctions.areaDelta`, plus `isValid` on
#    every result.
#
# Skips and known defects are never silent: `overlay_skiplist.jl` carries a
# written justification for every entry, and this file pins the per-file
# pass/skip counts and the exact set of currently-failing robust cases, so both
# a regression and a fix surface as a mismatch.
#
# Runtime: ~15 s warm for both runs.

using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

include(joinpath(@__DIR__, "..", "..", "..", "external", "jts", "overlay_runner.jl"))

const JTS_OVERLAY_DIR = joinpath(@__DIR__, "..", "..", "..", "data", "jts", "overlay")
const JTS_ROBUST_DIR = joinpath(@__DIR__, "..", "..", "..", "data", "jts", "overlay_robust")

const OPCODE = Dict(:intersection => GO.OVERLAY_INTERSECTION, :union => GO.OVERLAY_UNION,
                    :difference => GO.OVERLAY_DIFFERENCE, :symdifference => GO.OVERLAY_SYMDIFFERENCE)

go_overlay(op, a, b) = GO._overlay_ng(GO.Planar(), OPCODE[op], a, b; exact = GO.True())

_xml_files(dir) = sort!(filter!(f -> endswith(f, ".xml"), readdir(dir; join = true)))

# Whether the engine accepts point inputs yet (phase 3, landing on a separate
# branch). Probed rather than assumed, so the point-input suites switch on by
# themselves the moment support lands; the pinned counts below are then relaxed
# to lower bounds instead of going stale.
const POINT_INPUTS_SUPPORTED = try
    GO._overlay_ng(GO.Planar(), GO.OVERLAY_INTERSECTION, GI.Point((0.5, 0.5)),
        GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]); exact = GO.True())
    true
catch
    false
end

# ---------------------------------------------------------------------------
# 1. Expected-result conformance over the general/misc overlay suites
# ---------------------------------------------------------------------------

# Expected (pass, skip) counts per file, pinned as of this branch.
# All 346 skips are structural and recorded by the runner: 310 ops with a point
# operand (`TestOverlayP*`, `TestNGOverlayP`, and the point half of the two
# Empty files) and 36 with a GEOMETRYCOLLECTION operand (`TestNGOverlayGC`, the
# GC half of `TestNGOverlayEmpty`) — both unsupported by the engine on this
# branch — plus the 6 documented `OVERLAY_SKIPLIST` entries.
const OVERLAY_EXPECTED_COUNTS = [
    ("TestNGOverlayA.xml",      88,   0),
    ("TestNGOverlayEmpty.xml",  29,  50),
    ("TestNGOverlayGC.xml",      0,  16),
    ("TestNGOverlayL.xml",      55,   0),
    ("TestNGOverlayP.xml",       0,  60),
    ("TestOverlayAA.xml",       44,   0),
    ("TestOverlayEmpty.xml",   132, 162),
    ("TestOverlayLA.xml",       15,   1),
    ("TestOverlayLL.xml",       25,   0),
    ("TestOverlayMisc.xml",      0,   5),
    ("TestOverlayPA.xml",        0,   9),
    ("TestOverlayPL.xml",        0,  20),
    ("TestOverlayPP.xml",        0,  29),
]

@testset "JTS overlay XML conformance suite" begin
    files = _xml_files(JTS_OVERLAY_DIR)
    @test basename.(files) == first.(OVERLAY_EXPECTED_COUNTS)
    summary = run_overlay_cases(go_overlay, files; point_inputs = POINT_INPUTS_SUPPORTED)
    for f in summary.failures
        println("overlay XML failure: $(f.file) case $(f.case_index) $(f.op) $(f.arg_order) — $(f.detail)")
    end
    @test sum(s -> s.n_fail, summary.per_file) == 0
    for (s, (file, n_pass, n_skip)) in zip(summary.per_file, OVERLAY_EXPECTED_COUNTS)
        @test s.file == file
        if POINT_INPUTS_SUPPORTED
            #-- point support turns skipped ops into executed ones; the pins
            #-- become lower bounds rather than going stale.
            @test (s.file, s.n_pass >= n_pass) == (file, true)
        else
            @test (s.file, s.n_pass) == (file, n_pass)
            @test (s.file, s.n_skip) == (file, n_skip)
        end
    end
    #-- every skip is one of the four accounted-for kinds
    for sk in summary.skipped
        @test sk.reason in ("in skiplist", "point input unsupported",
            "geometry collection input unsupported", "binary op without a B geometry")
    end
end

# ---------------------------------------------------------------------------
# 2. The positive claim behind the three "stale expectation" skiplist entries
# ---------------------------------------------------------------------------
#
# Those three ops are skipped above because GEOS `equals` against the XML's
# expected WKT is false. Skipping alone would lose the finding, so the reason
# they are benign is asserted here: on all three, GEOS 3.14 also disagrees with
# the XML, and our answer is at least as close to it as GEOS's is.

@testset "stale XML expectations (TestOverlayMisc)" begin
    run = load_test_run(joinpath(JTS_OVERLAY_DIR, "TestOverlayMisc.xml"))
    #-- `area_rtol` is per case: cases 1 and 3 are large polygons where all three
    #-- answers agree to 1e-12 relative, case 4 is a 181-unit sliver cut out of
    #-- 8.5e5-magnitude coordinates, where a 2-ulp vertex difference moves the
    #-- area by 2e-11 relative.
    for (case_index, opname, sym_tol, area_rtol) in
            ((1, "union", 1e-9, 1e-12), (3, "intersection", 1e-6, 1e-12), (4, "intersection", 1e-8, 1e-10))
        case = run.cases[case_index]
        a = jts_overlay_operand(case.geom_a)
        b = jts_overlay_operand(case.geom_b)
        la, lb = overlay_to_lg(a), overlay_to_lg(b)
        op = Symbol(opname)
        want = overlay_to_lg(first(i for i in case.items
            if lowercase(i.operation) == opname).expected_result)
        geos = op === :union ? LG.union(la, lb) : LG.intersection(la, lb)
        ours = overlay_to_lg(go_overlay(op, a, b))
        #-- GEOS disagrees with the XML too, so the expectation is the outlier
        @test !LG.equals(geos, want)
        #-- and our answer is no further from the XML than GEOS's is
        d_ours = LG.area(LG.symmetricDifference(ours, want))
        d_geos = LG.area(LG.symmetricDifference(geos, want))
        @test d_ours <= max(d_geos, sym_tol)
        #-- all three engines agree on the area
        @test isapprox(LG.area(ours), LG.area(geos); rtol = area_rtol)
        @test isapprox(LG.area(ours), LG.area(want); rtol = 1e-9)
    end
    #-- case 3's intersection is `equals` to GEOS exactly: our only difference
    #-- from the XML is that its WKT is written to 8 decimal places.
    case3 = run.cases[3]
    a3, b3 = jts_overlay_operand(case3.geom_a), jts_overlay_operand(case3.geom_b)
    @test LG.equals(overlay_to_lg(go_overlay(:intersection, a3, b3)),
                    LG.intersection(overlay_to_lg(a3), overlay_to_lg(b3)))
end

# ---------------------------------------------------------------------------
# 3. Oracle-free identity sweep over the robust corpus
# ---------------------------------------------------------------------------
#
# `overlay_area_delta` is JTS's own `areaDelta`: the five overlay area
# identities must hold, and every result must be `isValid`.
#
# Areas are measured with `LG.area`, not `GO.area`, because this corpus is
# projected data with coordinates up to ~1e6 and `GO.area(Planar())` runs an
# untranslated shoelace: its absolute round-off there is ~2^-11, i.e. ~1e-8
# relative, which would swamp the identity. GEOS (like JTS) subtracts the first
# vertex before summing, which keeps the residual at ~1e-14 relative. That is a
# reported `GO.area` finding, not an overlay one, so the sweep uses the accurate
# oracle and leaves `GO.area` to the Natural Earth sweeps (coordinates ≤ 180,
# where the difference does not arise).

lg_area(g) = LG.area(overlay_to_lg(g))
lg_valid(g) = try LG.isValid(overlay_to_lg(g)) catch; false end

@testset "JTS robust overlay corpus — area identities + validity" begin
    files = _xml_files(JTS_ROBUST_DIR)
    @test length(files) == 44
    summary = run_overlay_identity_cases(go_overlay, lg_area, files;
        valid_fn = lg_valid, result_valid_fn = lg_valid, rtol = 1e-9, record_tests = false)
    observed = Set((f.file, f.case_index) for f in summary.failures)
    #-- Pinned: the exact set of cases the engine currently gets wrong, each one
    #-- an instance of the collapsed-result-ring defect documented in
    #-- `overlay_skiplist.jl`. A fix shrinks this set, a regression grows it;
    #-- either way the assertion below fails and the ledger must be updated.
    newly_failing = sort!(collect(setdiff(observed, ROBUST_KNOWN_DEFECTS)))
    newly_passing = sort!(collect(setdiff(ROBUST_KNOWN_DEFECTS, observed)))
    isempty(newly_failing) || println("robust corpus — NEW failures: ", newly_failing)
    isempty(newly_passing) || println("robust corpus — now passing (update the ledger): ", newly_passing)
    for f in summary.failures
        (f.file, f.case_index) in ROBUST_KNOWN_DEFECTS && continue
        println("robust corpus failure: $(f.file) case $(f.case_index) — $(f.detail)")
    end
    @test isempty(newly_failing)
    @test isempty(newly_passing)
    n_run = sum(s -> s.n_pass + s.n_fail, summary.per_file)
    @test n_run >= 60          # the corpus must not silently shrink
    @test length(summary.skipped) == sum(s -> s.n_skip, summary.per_file)
end

# ---------------------------------------------------------------------------
# 4. Reduced reproducer for the collapsed-result-ring defect
# ---------------------------------------------------------------------------
#
# `TestOverlay-jts-798.xml` case 1, an 8-point pair of VALID polygons (A is a
# very thin sliver triangle whose apex sits just past B's top edge). GEOS emits
# the collapsed part as a result LINESTRING; we keep it inside the area ring as
# an out-and-back spike, which makes the union/difference/symdifference results
# invalid and leaves INTERSECTION — whose result is entirely collapsed — with
# holes and no shell at all.
#
# These are `@test_broken`: they must start passing the moment the defect is
# fixed, at which point the ledger above shrinks too.

@testset "collapsed result ring (jts-798 case 1) — known defect" begin
    A = GO.tuples(LG.readgeom("POLYGON ((66697.40120137333 185279.95469107336, " *
        "66698.375 185273.625, 66697.375 185280.125, 66697.40120137333 185279.95469107336))"))
    B = GO.tuples(LG.readgeom("POLYGON ((66710 185280, 66710 185260, 66690 185260, " *
        "66690 185280, 66710 185280))"))
    #-- both inputs are valid, so the engine's validity contract is satisfied
    @test LG.isValid(overlay_to_lg(A))
    @test LG.isValid(overlay_to_lg(B))
    #-- INTERSECTION: the whole result collapses to a line; we throw instead
    @test_broken (go_overlay(:intersection, A, B); true)
    #-- the other three produce the right area but an invalid geometry
    for op in (:union, :difference, :symdifference)
        r = go_overlay(op, A, B)
        lgr = op === :union ? LG.union(overlay_to_lg(A), overlay_to_lg(B)) :
              op === :difference ? LG.difference(overlay_to_lg(A), overlay_to_lg(B)) :
              LG.symmetricDifference(overlay_to_lg(A), overlay_to_lg(B))
        @test isapprox(lg_area(r), LG.area(lgr); atol = 1e-9)
        @test_broken lg_valid(r)
    end
end
