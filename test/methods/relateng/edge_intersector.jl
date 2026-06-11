# Tests for the RelateNG edge segment intersector (edge_intersector.jl), the
# port of JTS EdgeSegmentIntersector.java over the symbolic segment
# classification `rk_classify_intersection` (no JTS JUnit counterpart; per the
# implementation plan, Task 19).
#
# Strategy: hand-built two-string configurations are fed through
# `process_intersections!` for every segment pair (a nested-loop stand-in for
# the Task 20 enumerator), and the exact set of (NodeKey, section count)
# recorded in a `RelateMatrixPredicate`-backed TopologyComputer is asserted:
#
# - proper crossing (one symbolic crossing node, one section pair),
# - T-touch of a line end on a segment interior,
# - shared endpoint (two incidence flags, ONE geometric point),
# - collinear overlap (up to two distinct points, one section pair each),
# - adjacent-segment shared vertex (sections produced ONCE, not twice),
# - ring-wraparound vertex (closing vertex attributed to the first segment),
# - same-string and same-geometry handling.
#
# The once-only predicate `_is_canonical_incidence` is also unit-tested
# directly (mid-string vertex, string start, string end, ring wraparound).

using Test
import GeometryOps as GO
import GeometryOps: Planar, True, False
import GeoInterface as GI

const M = Planar()
const EX = True()

rgeom(g) = GO.RelateGeometry(M, g; exact = EX)

# A TopologyComputer with a full-matrix predicate attached; returns both.
function im_computer(rga::GO.RelateGeometry, rgb::GO.RelateGeometry)
    pred = GO.RelateMatrixPredicate()
    return GO.TopologyComputer(pred, rga, rgb), pred
end
im_computer(ga, gb) = im_computer(rgeom(ga), rgeom(gb))

imstr(pred) = string(GO.result_im(pred))

# Extract the single segment string of a one-element geometry.
single_ss(rg, is_a) = only(GO.extract_segment_strings(rg, is_a, nothing))

# Nested-loop enumeration of all segment pairs of two strings (the Task 20
# accelerated enumerator replaces this; semantics must match).
function run_pairs!(tc, ss0, ss1)
    for i0 in 1:(length(ss0.pts) - 1), i1 in 1:(length(ss1.pts) - 1)
        GO.process_intersections!(tc, ss0, i0, ss1, i1)
    end
    return tc
end

# The recorded node topology as a Dict of NodeKey => number of sections.
node_counts(tc) = Dict(k => length(v.sections) for (k, v) in tc.node_sections)

# Build a computer over two geometries and run all A x B segment pairs.
function run_case(ga, gb)
    rga, rgb = rgeom(ga), rgeom(gb)
    tc, pred = im_computer(rga, rgb)
    run_pairs!(tc, single_ss(rga, true), single_ss(rgb, false))
    return tc, pred
end

@testset "_is_canonical_incidence" begin
    # open three-point line: segments 1 = (0,0)-(1,0), 2 = (1,0)-(2,0)
    line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])
    ss = single_ss(rgeom(line), true)

    # mid-string vertex: owned by the segment it starts (seg 2), not the
    # segment it ends (seg 1)
    @test !GO._is_canonical_incidence(ss, 1, (1.0, 0.0))
    @test GO._is_canonical_incidence(ss, 2, (1.0, 0.0))
    # string start: owned by the first segment
    @test GO._is_canonical_incidence(ss, 1, (0.0, 0.0))
    # string end: final segment of an open string owns its endpoint
    @test GO._is_canonical_incidence(ss, 2, (2.0, 0.0))
    # segment-interior points are always owned by their segment
    @test GO._is_canonical_incidence(ss, 1, (0.5, 0.0))
    @test GO._is_canonical_incidence(ss, 2, (1.5, 0.0))

    # ring wraparound: the closing vertex is owned by the first segment,
    # NOT by the final segment that ends there
    sq = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    ring = single_ss(rgeom(sq), true)
    @test ring.pts[1] == ring.pts[end]   # closed
    n_seg = length(ring.pts) - 1
    start_pt = ring.pts[1]
    @test GO._is_canonical_incidence(ring, 1, start_pt)
    @test !GO._is_canonical_incidence(ring, n_seg, start_pt)
    # an ordinary ring vertex behaves like a mid-string vertex
    mid_pt = ring.pts[2]
    @test !GO._is_canonical_incidence(ring, 1, mid_pt)
    @test GO._is_canonical_incidence(ring, 2, mid_pt)
end

@testset "proper crossing" begin
    line_a = GI.LineString([(-1.0, 0.0), (1.0, 0.0)])
    line_b = GI.LineString([(0.0, -1.0), (0.0, 1.0)])
    tc, pred = run_case(line_a, line_b)
    key = GO.crossing_node((-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0))
    @test node_counts(tc) == Dict(key => 2)
    # the crossing lies in both line interiors
    @test imstr(pred) == "0FFFFFFF2"

    # argument order is normalized by is_a: feeding (B, A) gives an
    # identical computer state
    rga, rgb = rgeom(line_a), rgeom(line_b)
    tc2, pred2 = im_computer(rga, rgb)
    run_pairs!(tc2, single_ss(rgb, false), single_ss(rga, true))
    @test node_counts(tc2) == Dict(key => 2)
    @test imstr(pred2) == "0FFFFFFF2"
end

@testset "T-touch: line end on segment interior" begin
    line_a = GI.LineString([(-1.0, 0.0), (1.0, 0.0)])
    # touch point is the END of B's final (open) segment
    line_b = GI.LineString([(0.0, 1.0), (0.0, 0.0)])
    tc, _ = run_case(line_a, line_b)
    @test node_counts(tc) == Dict(GO.vertex_node((0.0, 0.0)) => 2)
    # the node lies in the interior of A's segment, so the A-side section is
    # NOT at a vertex; it IS at a vertex of B (B's line end), per the JTS
    # createNodeSection semantics
    sections = tc.node_sections[GO.vertex_node((0.0, 0.0))].sections
    @test !GO.is_node_at_vertex(only(filter(GO.is_a, sections)))
    @test GO.is_node_at_vertex(only(filter(!GO.is_a, sections)))

    # touch point is the START of B
    line_b2 = GI.LineString([(0.0, 0.0), (0.0, 1.0)])
    tc, _ = run_case(line_a, line_b2)
    @test node_counts(tc) == Dict(GO.vertex_node((0.0, 0.0)) => 2)
    sections2 = tc.node_sections[GO.vertex_node((0.0, 0.0))].sections
    @test !GO.is_node_at_vertex(only(filter(GO.is_a, sections2)))
    @test GO.is_node_at_vertex(only(filter(!GO.is_a, sections2)))
end

@testset "shared endpoint: one point, one section pair" begin
    # a0 == b0 sets both a0_on_b and b0_on_a — ONE geometric point
    line_a = GI.LineString([(0.0, 0.0), (1.0, 0.0)])
    line_b = GI.LineString([(0.0, 0.0), (0.0, 1.0)])
    tc, _ = run_case(line_a, line_b)
    @test node_counts(tc) == Dict(GO.vertex_node((0.0, 0.0)) => 2)
end

@testset "collinear overlap" begin
    # partial overlap: the interval ends are a1 = (3,0) and b0 = (1,0)
    line_a = GI.LineString([(0.0, 0.0), (3.0, 0.0)])
    line_b = GI.LineString([(1.0, 0.0), (4.0, 0.0)])
    tc, _ = run_case(line_a, line_b)
    @test node_counts(tc) == Dict(
        GO.vertex_node((1.0, 0.0)) => 2,
        GO.vertex_node((3.0, 0.0)) => 2,
    )

    # containment: B inside one segment of A — both interval ends are B's
    line_a2 = GI.LineString([(0.0, 0.0), (4.0, 0.0)])
    line_b2 = GI.LineString([(1.0, 0.0), (2.0, 0.0)])
    tc, _ = run_case(line_a2, line_b2)
    @test node_counts(tc) == Dict(
        GO.vertex_node((1.0, 0.0)) => 2,
        GO.vertex_node((2.0, 0.0)) => 2,
    )

    # identical segments (reversed orientation): all four flags set, but
    # only two distinct points
    line_b3 = GI.LineString([(3.0, 0.0), (0.0, 0.0)])
    tc, _ = run_case(line_a, line_b3)
    @test node_counts(tc) == Dict(
        GO.vertex_node((0.0, 0.0)) => 2,
        GO.vertex_node((3.0, 0.0)) => 2,
    )
end

@testset "adjacent-segment shared vertex: sections once, not twice" begin
    # A has a mid-string vertex at the origin; B passes through it. Both of
    # A's incident segments touch B there, but only the segment STARTING at
    # the vertex owns the incidence.
    line_a = GI.LineString([(-1.0, 1.0), (0.0, 0.0), (1.0, 1.0)])
    line_b = GI.LineString([(-2.0, 0.0), (2.0, 0.0)])
    tc, _ = run_case(line_a, line_b)
    key = GO.vertex_node((0.0, 0.0))
    @test node_counts(tc) == Dict(key => 2)
    # the A section is built on segment 2, so it sees both incident edges
    sections = tc.node_sections[key].sections
    nsa = only(filter(GO.is_a, sections))
    @test GO.get_vertex(nsa, 0) == (-1.0, 1.0)
    @test GO.get_vertex(nsa, 1) == (1.0, 1.0)
    @test GO.is_node_at_vertex(nsa)
    # the node lies in the interior of B's segment, so the B-side section is
    # not at a vertex
    @test !GO.is_node_at_vertex(only(filter(!GO.is_a, sections)))

    # both strings have a mid-string vertex at the node: 4 touching pairs,
    # still exactly one section pair
    line_b2 = GI.LineString([(-1.0, -1.0), (0.0, 0.0), (1.0, -1.0)])
    tc, _ = run_case(line_a, line_b2)
    @test node_counts(tc) == Dict(key => 2)
end

@testset "ring wraparound vertex" begin
    sq = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    line_b = GI.LineString([(-1.0, -1.0), (0.0, 0.0)])
    rga, rgb = rgeom(sq), rgeom(line_b)
    tc, _ = im_computer(rga, rgb)
    ring = single_ss(rga, true)
    run_pairs!(tc, ring, single_ss(rgb, false))
    key = GO.vertex_node((0.0, 0.0))
    @test node_counts(tc) == Dict(key => 2)
    # the ring section wraps around: its incident edges are the ring
    # vertices either side of the closing point
    nsa = only(filter(GO.is_a, tc.node_sections[key].sections))
    @test GO.get_vertex(nsa, 0) == ring.pts[end - 1]
    @test GO.get_vertex(nsa, 1) == ring.pts[2]
end

@testset "same string and same geometry" begin
    line_pts = [(-1.0, 1.0), (0.0, 0.0), (1.0, 1.0)]
    rga, rgb = rgeom(GI.LineString(line_pts)), rgeom(GI.LineString([(5.0, 5.0), (6.0, 5.0)]))
    tc, pred = im_computer(rga, rgb)
    ssa = single_ss(rga, true)
    # a segment is never intersected with itself
    GO.process_intersections!(tc, ssa, 1, ssa, 1)
    @test isempty(tc.node_sections)
    # adjacent segments of the SAME string share a vertex, but a plain
    # string vertex is not a node: no sections in either pair order
    GO.process_intersections!(tc, ssa, 1, ssa, 2)
    GO.process_intersections!(tc, ssa, 2, ssa, 1)
    @test isempty(tc.node_sections)
    @test imstr(pred) == "FFFFFFFF2"

    # a self-intersection within geometry A records sections (for node
    # analysis) but performs no A/B matrix update
    ml = GI.MultiLineString([
        GI.LineString([(-1.0, -1.0), (1.0, 1.0)]),
        GI.LineString([(-1.0, 1.0), (1.0, -1.0)]),
    ])
    rga2 = rgeom(ml)
    tc, pred = im_computer(rga2, rgb)
    ssa1, ssa2 = GO.extract_segment_strings(rga2, true, nothing)
    GO.process_intersections!(tc, ssa1, 1, ssa2, 1)
    key = GO.crossing_node((-1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (1.0, -1.0))
    @test node_counts(tc) == Dict(key => 2)
    @test imstr(pred) == "FFFFFFFF2"
end

@testset "disjoint segments record nothing" begin
    tc, pred = run_case(
        GI.LineString([(0.0, 0.0), (1.0, 0.0)]),
        GI.LineString([(0.0, 1.0), (1.0, 1.0)]))
    @test isempty(tc.node_sections)
    @test imstr(pred) == "FFFFFFFF2"
end

@testset "collinear overlap across multi-segment strings" begin
    # B overlaps both segments of A; the overlap interval ends (1,0) and
    # (3,0) lie in segment interiors of A, and A's mid-string vertex (2,0)
    # lies in B's interior. Three nodes, one section pair each, regardless
    # of B's direction.
    line_a = GI.LineString([(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)])
    expected = Dict(
        GO.vertex_node((1.0, 0.0)) => 2,
        GO.vertex_node((2.0, 0.0)) => 2,
        GO.vertex_node((3.0, 0.0)) => 2,
    )
    for bpts in ([(1.0, 0.0), (3.0, 0.0)], [(3.0, 0.0), (1.0, 0.0)])
        tc, _ = run_case(line_a, GI.LineString(bpts))
        @test node_counts(tc) == expected
    end
end

@testset "ring x ring sharing the closing corner" begin
    # two squares touching only at the point that closes BOTH rings: the
    # wraparound canonicality rule must fire on both sides at once, so the
    # corner is recorded exactly once with one section pair
    sq_a = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
    sq_b = GI.Polygon([[(0.0, 0.0), (-1.0, 0.0), (-1.0, -1.0), (0.0, -1.0), (0.0, 0.0)]])
    tc, _ = run_case(sq_a, sq_b)
    @test node_counts(tc) == Dict(GO.vertex_node((0.0, 0.0)) => 2)
end

# ---------------------------------------------------------------------------
# Task 20: accelerated edge set enumeration (`process_edge_intersections!`)
# ---------------------------------------------------------------------------

# Run the full A x B enumeration through `process_edge_intersections!` with a
# given accelerator.
function run_enum(ga, gb, accelerator)
    rga, rgb = rgeom(ga), rgeom(gb)
    tc, pred = im_computer(rga, rgb)
    ssa_list = GO.extract_segment_strings(rga, true, nothing)
    ssb_list = GO.extract_segment_strings(rgb, false, nothing)
    GO.process_edge_intersections!(tc, ssa_list, ssb_list, accelerator)
    return tc, pred
end

# Ground truth: the plain all-pairs enumeration used by the Task 19 tests,
# extended over string lists.
function run_all_pairs(ga, gb)
    rga, rgb = rgeom(ga), rgeom(gb)
    tc, pred = im_computer(rga, rgb)
    for ssa in GO.extract_segment_strings(rga, true, nothing),
            ssb in GO.extract_segment_strings(rgb, false, nothing)
        run_pairs!(tc, ssa, ssb)
    end
    return tc, pred
end

# A closed n-gon ring approximating a circle of radius `r` centered at
# `(cx, cy)` (many segments; used to exceed the AutoAccelerator threshold).
ngon(cx, cy, r, n) = GI.Polygon([[
    (cx + r * cos(θ), cy + r * sin(θ)) for θ in range(0, 2π; length = n + 1)
]])

@testset "accelerators produce identical computer state" begin
    fixtures = [
        # proper crossing
        (GI.LineString([(-1.0, 0.0), (1.0, 0.0)]), GI.LineString([(0.0, -1.0), (0.0, 1.0)])),
        # T-touch, both directions of B
        (GI.LineString([(-1.0, 0.0), (1.0, 0.0)]), GI.LineString([(0.0, 1.0), (0.0, 0.0)])),
        (GI.LineString([(-1.0, 0.0), (1.0, 0.0)]), GI.LineString([(0.0, 0.0), (0.0, 1.0)])),
        # shared endpoint
        (GI.LineString([(0.0, 0.0), (1.0, 0.0)]), GI.LineString([(0.0, 0.0), (0.0, 1.0)])),
        # collinear overlaps
        (GI.LineString([(0.0, 0.0), (3.0, 0.0)]), GI.LineString([(1.0, 0.0), (4.0, 0.0)])),
        (GI.LineString([(0.0, 0.0), (4.0, 0.0)]), GI.LineString([(1.0, 0.0), (2.0, 0.0)])),
        (GI.LineString([(0.0, 0.0), (3.0, 0.0)]), GI.LineString([(3.0, 0.0), (0.0, 0.0)])),
        # collinear overlap across multi-segment strings, both directions
        (GI.LineString([(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)]), GI.LineString([(1.0, 0.0), (3.0, 0.0)])),
        (GI.LineString([(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)]), GI.LineString([(3.0, 0.0), (1.0, 0.0)])),
        # adjacent-segment shared vertex
        (GI.LineString([(-1.0, 1.0), (0.0, 0.0), (1.0, 1.0)]), GI.LineString([(-2.0, 0.0), (2.0, 0.0)])),
        # multi-string side: a MultiLineString crossing a line twice
        (GI.MultiLineString([
            GI.LineString([(-1.0, -1.0), (-1.0, 1.0)]),
            GI.LineString([(1.0, -1.0), (1.0, 1.0)]),
        ]), GI.LineString([(-2.0, 0.0), (2.0, 0.0)])),
        # ring x ring at the closing corner; ring x line at the wraparound
        (GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]]),
         GI.Polygon([[(0.0, 0.0), (-1.0, 0.0), (-1.0, -1.0), (0.0, -1.0), (0.0, 0.0)]])),
        (GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]]),
         GI.LineString([(-1.0, -1.0), (0.0, 0.0)])),
        # disjoint
        (GI.LineString([(0.0, 0.0), (1.0, 0.0)]), GI.LineString([(0.0, 1.0), (1.0, 1.0)])),
        # heavily-overlapping many-segment rings: above the AutoAccelerator
        # threshold, so Auto takes the tree path here
        (ngon(0.0, 0.0, 1.0, 64), ngon(0.1, 0.0, 1.0, 64)),
    ]
    for (ga, gb) in fixtures
        tc_truth, pred_truth = run_all_pairs(ga, gb)
        for acc in (GO.NestedLoop(), GO.DoubleSTRtree(), GO.AutoAccelerator())
            tc, pred = run_enum(ga, gb, acc)
            @test node_counts(tc) == node_counts(tc_truth)
            @test imstr(pred) == imstr(pred_truth)
        end
    end
end

# A TopologyPredicate wrapper that counts `is_known` queries. The enumerator
# checks `is_result_known(computer)` (which calls `is_known(predicate)`)
# exactly once after each segment pair it processes, so the count equals the
# number of pairs processed.
mutable struct CountingPredicate{P <: GO.TopologyPredicate} <: GO.TopologyPredicate
    const inner::P
    known_checks::Int
end
CountingPredicate(inner) = CountingPredicate(inner, 0)
GO.predicate_name(p::CountingPredicate) = GO.predicate_name(p.inner)
GO.require_self_noding(p::CountingPredicate) = GO.require_self_noding(p.inner)
GO.require_interaction(p::CountingPredicate) = GO.require_interaction(p.inner)
GO.require_covers(p::CountingPredicate, is_source_a::Bool) = GO.require_covers(p.inner, is_source_a)
GO.require_exterior_check(p::CountingPredicate, is_source_a::Bool) = GO.require_exterior_check(p.inner, is_source_a)
GO.init_dims!(p::CountingPredicate, dimA::Integer, dimB::Integer) = GO.init_dims!(p.inner, dimA, dimB)
GO.init_bounds!(p::CountingPredicate, extA, extB) = GO.init_bounds!(p.inner, extA, extB)
GO.update_dim!(p::CountingPredicate, locA, locB, dim) = GO.update_dim!(p.inner, locA, locB, dim)
GO.finish!(p::CountingPredicate) = GO.finish!(p.inner)
GO.predicate_value(p::CountingPredicate) = GO.predicate_value(p.inner)
function GO.is_known(p::CountingPredicate)
    p.known_checks += 1
    return GO.is_known(p.inner)
end

# Enumerate with a counting predicate; returns the wrapper for inspection.
function run_counted(ga, gb, accelerator, inner_pred)
    rga, rgb = rgeom(ga), rgeom(gb)
    pred = CountingPredicate(inner_pred)
    tc = GO.TopologyComputer(pred, rga, rgb)
    ssa_list = GO.extract_segment_strings(rga, true, nothing)
    ssb_list = GO.extract_segment_strings(rgb, false, nothing)
    GO.process_edge_intersections!(tc, ssa_list, ssb_list, accelerator)
    return pred
end

@testset "early exit when the result is known" begin
    # two heavily-overlapping many-segment rings: an `intersects` predicate
    # becomes known at the first proper edge crossing, so the traversal must
    # stop long before enumerating every extent-overlapping pair
    ga = ngon(0.0, 0.0, 1.0, 64)
    gb = ngon(0.1, 0.0, 1.0, 64)
    for acc in (GO.NestedLoop(), GO.DoubleSTRtree())
        # baseline: the full-matrix predicate is never known early, so its
        # count is the total number of pairs the enumeration would process
        full = run_counted(ga, gb, acc, GO.RelateMatrixPredicate())
        @test !GO.is_known(full.inner)
        @test full.known_checks > 0

        early = run_counted(ga, gb, acc, GO.pred_intersects())
        @test GO.is_known(early.inner)
        @test GO.predicate_value(early.inner)
        @test 0 < early.known_checks < full.known_checks
    end
end

# ---------------------------------------------------------------------------
# Self-pair enumeration (`process_self_intersections!`)
# ---------------------------------------------------------------------------

@testset "self-intersections: accelerators produce identical computer state" begin
    # A self-intersecting zigzag polyline, well above the AutoAccelerator
    # threshold: 41 zigzag vertices oscillating across y = 0, then a return
    # segment along y = 0 that properly crosses every zigzag segment.
    zig_pts = [(Float64(i), iseven(i) ? 1.0 : -1.0) for i in 0:40]
    push!(zig_pts, (40.5, 0.0))
    push!(zig_pts, (-0.5, 0.0))
    zig = GI.LineString(zig_pts)
    @test length(zig_pts) - 1 >= GO.GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS

    # disjoint far-away B side, so the computer state reflects only A's
    # self-noding
    rgb_far = rgeom(GI.LineString([(100.0, 100.0), (101.0, 100.0)]))

    function run_self(acc)
        rga = rgeom(zig)
        tc, pred = im_computer(rga, rgb_far)
        GO.process_self_intersections!(tc,
            GO.extract_segment_strings(rga, true, nothing), acc)
        return node_counts(tc), imstr(pred)
    end

    counts_loop, im_loop = run_self(GO.NestedLoop())
    counts_tree, im_tree = run_self(GO.DoubleSTRtree())
    @test counts_loop == counts_tree
    @test im_loop == im_tree
    # sanity: the return segment properly crosses all 40 zigzag segments
    @test count(k -> k.is_crossing, keys(counts_loop)) == 40
    # self pairs record sections but never update the A/B matrix
    @test im_loop == "FFFFFFFF2"
end

@testset "self-intersections: early exit when the result is known" begin
    # Mirror JTS computeEdgesAll: feed the strings of BOTH inputs through
    # one self-pair enumeration, so the cross (A x B) pairs update the
    # matrix and an `intersects` predicate becomes known at the first proper
    # crossing — the traversal must then stop early.
    ga = ngon(0.0, 0.0, 1.0, 64)
    gb = ngon(0.1, 0.0, 1.0, 64)
    function run_self_counted(inner_pred, acc)
        rga, rgb = rgeom(ga), rgeom(gb)
        pred = CountingPredicate(inner_pred)
        tc = GO.TopologyComputer(pred, rga, rgb)
        ss_list = vcat(GO.extract_segment_strings(rga, true, nothing),
            GO.extract_segment_strings(rgb, false, nothing))
        GO.process_self_intersections!(tc, ss_list, acc)
        return pred
    end
    for acc in (GO.NestedLoop(), GO.DoubleSTRtree())
        # baseline: the full-matrix predicate is never known early, so its
        # count is the total number of pairs the enumeration would process
        full = run_self_counted(GO.RelateMatrixPredicate(), acc)
        @test !GO.is_known(full.inner)
        @test full.known_checks > 0

        early = run_self_counted(GO.pred_intersects(), acc)
        @test GO.is_known(early.inner)
        @test GO.predicate_value(early.inner)
        @test 0 < early.known_checks < full.known_checks
    end
end
