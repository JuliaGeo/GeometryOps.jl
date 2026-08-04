# Tests for the OverlayNG phase-2a engine core (design §3): the half-edge graph
# over the phase-1 `NodedArrangement` — edge sources / depth deltas, the shared
# `OverlayLabel`, the node-pair edge merger, and the CCW-ordered half-edge stars.

using Test
include(joinpath(@__DIR__, "common.jl"))
import GeometryOps: Planar, Spherical, True, False
import GeometryOps: POS_LEFT, POS_RIGHT, POS_ON
import GeometryOps: LOC_INTERIOR, LOC_EXTERIOR, LOC_BOUNDARY, LOC_NONE
import GeometryOps: DIM_A, DIM_L, DIM_COLLAPSE, DIM_NOT_PART

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_crossing(arr) = [i for i in 1:GO.num_nodes(arr) if arr.nodes.keys[i].is_crossing]

# The onext cycle (outgoing half-edge indices, CCW) at node `nid`.
function star_of(g, nid)
    ne = GO.graph_node_edge(g, nid)
    ne == 0 && return Int32[]
    star = Int32[]
    e = ne
    while true
        push!(star, e)
        e = GO.he_onext(g.edges, e)
        e == ne && break
    end
    return star
end

# Every half-edge has a well-formed sym, every node's star is a CCW-ordered
# permutation of its outgoing edges, and degrees match the arrangement incidence.
#
# `fold = true` reports one assertion per property instead of one per element.
# The Natural Earth loop passes it: over 30 country pairs the per-element form
# ran 36,668 assertions to establish the same seven properties the synthetic
# fixtures establish, and a conjunction fails exactly when one of its terms does.
# The small fixtures stay per-element, where a failing index is the diagnosis.
#
# `dest == sym.origin` used to be checked here and is not, because
# `he_dest(edges, i)` is DEFINED as `edges[he_sym(edges, i)].origin`
# (`half_edge.jl:46`) — substituting the definition makes the assertion `x == x`.
function check_graph(m, g; fold = false)
    chk(v) = fold ? v : (@test(v); true)
    sym_ok = true
    for i in eachindex(g.edges)
        s = GO.he_sym(g.edges, i)
        sym_ok &= chk(GO.he_sym(g.edges, s) == i)              # sym is involutive
        sym_ok &= chk(s != i)                                  # never its own sym
    end
    deg_ok = true; origin_ok = true; ccw_ok = true
    for nid in 1:GO.num_nodes(g.arr)
        star = star_of(g, nid)
        isempty(star) && continue
        incidence = count(i -> g.edges[i].origin == nid, eachindex(g.edges))
        deg_ok &= chk(length(star) == incidence)               # degree == incidence
        deg_ok &= chk(GO.he_degree(g.edges, star[1]) == incidence)
        origin_ok &= chk(all(i -> g.edges[i].origin == nid, star))
        #-- strictly CCW-increasing from the representative (lowest) edge
        for t in 1:(length(star) - 1)
            ccw_ok &= chk(GO.he_compare_angular(m, g.edges, g.arr.nodes.keys,
                                                star[t], star[t + 1]; exact = EX) < 0)
        end
    end
    if fold
        @test sym_ok
        @test deg_ok
        @test origin_ok
        @test ccw_ok
    end
    return nothing
end

# ---------------------------------------------------------------------------
# 1. Angular star ordering torture set (§3.1)
# ---------------------------------------------------------------------------

@testset "crossing-node star ordering (both manifolds)" begin
    # row 1: two A lines and one B line through the origin -> one merged node of
    # degree 6 (each line contributes an in- and an out-going half-edge), exercising
    # rk_compare_edge_dir's exact-rational foreign-direction path at a crossing apex.
    # row 2: a near-equatorial line and a meridian crossing at (0,0) -> degree 4,
    # the spherical tangent-ordering case.
    for (A, B, deg) in (
            (GI.MultiLineString([[(-1.0, -1.0), (1.0, 1.0)], [(-1.0, 1.0), (1.0, -1.0)]]),
             GI.LineString([(-1.0, 0.0), (1.0, 0.0)]), 6),
            (GI.LineString([(-10.0, 0.0), (10.0, 0.0)]),
             GI.LineString([(0.0, -10.0), (0.0, 10.0)]), 4),
        ), m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = EX)
        g = GO.OverlayGraph(m, arr; exact = EX)
        cid = findfirst(i -> arr.nodes.keys[i].is_crossing, 1:GO.num_nodes(arr))
        @test cid !== nothing
        @test GO.he_degree(g.edges, GO.graph_node_edge(g, cid)) == deg
        check_graph(m, g)
    end
end

@testset "collinear shared boundary -> one merged edge, both labels" begin
    # edge-adjacent unit squares sharing the vertex-identical edge x = 2
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = EX)
        g = GO.OverlayGraph(m, arr; exact = EX)
        @test length(_crossing(arr)) == 0
        #-- exactly one merged edge is a boundary of BOTH inputs (its sym pair)
        both = [i for i in eachindex(g.edges) if GO.is_boundary_both(g.edges[i].label)]
        @test length(both) == 2
        @test GO.he_sym(g.edges, both[1]) == both[2]
        check_graph(m, g)
    end
    #-- per-side locations on the plane: on the shared x=2 edge the WEST side
    #-- (x<2) is A-interior/B-exterior and the EAST side is B-interior/A-exterior,
    #-- regardless of the canonical ring reorientation applied at ingest — derive
    #-- west/east from the forward half-edge's actual heading.
    m = Planar()
    arr = GO.NodedArrangement(m, A, B; exact = EX)
    g = GO.OverlayGraph(m, arr; exact = EX)
    bf = first(i for i in eachindex(g.edges)
        if GO.is_boundary_both(g.edges[i].label) && g.edges[i].is_forward)
    lbl = g.edges[bf].label
    o = GO.node_point(g.arr, g.edges[bf].origin)
    d = GO.node_point(g.arr, GO.he_dest(g.edges, bf))
    heading_up = d[2] > o[2]                                        # both on x=2
    aL = GO.get_location(lbl, 0, POS_LEFT, true);  aR = GO.get_location(lbl, 0, POS_RIGHT, true)
    bL = GO.get_location(lbl, 1, POS_LEFT, true);  bR = GO.get_location(lbl, 1, POS_RIGHT, true)
    a_west, a_east = heading_up ? (aL, aR) : (aR, aL)
    b_west, b_east = heading_up ? (bL, bR) : (bR, bL)
    @test a_west == LOC_INTERIOR && a_east == LOC_EXTERIOR          # A is west
    @test b_west == LOC_EXTERIOR && b_east == LOC_INTERIOR          # B is east
    @test GO.is_boundary_touch(lbl)                                 # two areas meeting
end

# ---------------------------------------------------------------------------
# 2. Depth-delta algebra and Edge.merge semantics (§2.7, §3.3)
# ---------------------------------------------------------------------------

@testset "depth delta cross-check vs JTS convention" begin
    # JTS: canonical depth delta +1 == Exterior-on-Left (material interior RIGHT),
    # flipped to -1 otherwise. So material_interior_on_left => -1.
    for m in (Planar(), Spherical())
        ccw = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
        kccw = GO._to_kernel_points(m, GI.LinearRing(ccw))
        d = GO._ring_depth_delta(m, kccw, false; exact = EX)
        #-- the expectations are LITERALS, not expressions in `moil`. `_ring_depth_delta`
        #-- IS `moil ? -1 : 1` (edge_source.jl:39-41), so writing them that way made this
        #-- testset an implementation mirror: with `_ring_material_interior_on_left`
        #-- inverted at runtime, `d`, `moil` and both side locations flip together and
        #-- all ten assertions still pass. The convention this testset is named for was
        #-- being computed, never checked.
        @test GO._ring_material_interior_on_left(m, kccw, false; exact = EX)   # CCW shell: interior LEFT
        @test d == -1                                                         # ...so JTS's delta is -1
        @test GO._location_left(d)  == LOC_INTERIOR
        @test GO._location_right(d) == LOC_EXTERIOR
        # reversing the ring flips the delta and the side locations
        kcw = GO._to_kernel_points(m, GI.LinearRing(reverse(ccw)))
        dcw = GO._ring_depth_delta(m, kcw, false; exact = EX)
        @test !GO._ring_material_interior_on_left(m, kcw, false; exact = EX)
        @test dcw == 1
        @test GO._location_left(dcw) == GO._location_right(d)
    end
end

@testset "EdgeSourceInfo: shell/hole deltas and line" begin
    m = Planar()
    shell = GO._to_kernel_points(m, GI.LinearRing([(0.0,0.0),(4.0,0.0),(4.0,4.0),(0.0,4.0),(0.0,0.0)]))
    # a hole is labelled opposite to a shell of the same stored orientation
    d_shell = GO._ring_depth_delta(m, shell, false; exact = EX)
    d_hole  = GO._ring_depth_delta(m, shell, true;  exact = EX)
    @test d_hole == -d_shell
    #-- a line source carries no side labelling. Asserted through `_edge_source_info`,
    #-- which is what actually decides this: constructing an `EdgeSourceInfo` by hand
    #-- and reading its fields back only tests Julia's default constructor.
    line = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
    li = GO._edge_source_info(m, GO._overlay_segstrings(m, line, true)[1]; exact = EX)
    @test li.dim == DIM_L && li.depth_delta == 0
end

@testset "Edge.merge: same-dir sum, opposite-dir collapse, hole role" begin
    srcA(dd, hole) = GO.EdgeSourceInfo(Int8(0), DIM_A, hole, Int8(dd))
    ne(lo, hi) = GO.NodedEdge(Int32(1), Int32(1), Int32(lo), Int32(hi))
    # same direction: deltas add
    base = GO._merge_edge(ne(3, 7), srcA(1, false))
    GO._merge!(base, GO._merge_edge(ne(3, 7), srcA(1, false)))
    @test base.a_depth_delta == 2
    # opposite direction (a-b-a spike): deltas cancel -> DIM_COLLAPSE
    base = GO._merge_edge(ne(3, 7), srcA(-1, false))
    GO._merge!(base, GO._merge_edge(ne(7, 3), srcA(-1, false)))   # reversed
    @test base.a_depth_delta == 0
    lbl = GO._create_label(base)
    @test GO.is_collapse(lbl, 0)
    @test lbl.a_dim == DIM_COLLAPSE
    # hole-role merge: a shell contributor makes the merged edge a shell
    base = GO._merge_edge(ne(3, 7), srcA(-1, true))              # hole
    GO._merge!(base, GO._merge_edge(ne(3, 7), srcA(-1, false)))  # shell
    @test base.a_is_hole == false
    # two holes stay a hole
    base = GO._merge_edge(ne(3, 7), srcA(-1, true))
    GO._merge!(base, GO._merge_edge(ne(3, 7), srcA(-1, true)))
    @test base.a_is_hole == true
end

@testset "geometric a-b-a spike ring produces a collapse edge" begin
    # a self-touching ring with an out-and-back spike (2,4)->(2,6)->(2,4);
    # the noder pairs the two coincident segments and the merger collapses them.
    A = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (2.0, 4.0), (2.0, 6.0),
                     (2.0, 4.0), (0.0, 4.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    m = Planar()
    arr = GO.NodedArrangement(m, A, B; exact = EX)
    g = GO.OverlayGraph(m, arr; exact = EX)
    @test any(i -> GO.is_collapse(g.edges[i].label, 0), eachindex(g.edges))
    check_graph(m, g)
end

# ---------------------------------------------------------------------------
# 3. Label semantics — the four states, predicates, and the isForward swap
# ---------------------------------------------------------------------------

@testset "the four label states from Edge sources" begin
    # Boundary (area, nonzero delta) for A, NotPart for B
    lbl = GO._create_label(GO.MergeEdge(Int32(1), Int32(2), Int32(1), Int32(1),
        DIM_A, Int32(-1), false, DIM_NOT_PART, Int32(0), false))
    @test GO.is_boundary(lbl, 0) && GO.is_not_part(lbl, 1)
    @test GO.is_boundary_either(lbl) && !GO.is_boundary_both(lbl)
    @test GO.is_boundary_singleton(lbl)
    @test GO.has_sides(lbl, 0) && !GO.has_sides(lbl, 1)
    @test GO.get_line_location(lbl, 0) == LOC_INTERIOR   # boundary line loc is INTERIOR
    @test GO.is_known(lbl, 0) && !GO.is_known(lbl, 1)

    # Collapse (area, zero delta)
    lbl = GO._create_label(GO.MergeEdge(Int32(1), Int32(2), Int32(1), Int32(1),
        DIM_A, Int32(0), true, DIM_NOT_PART, Int32(0), false))
    @test GO.is_collapse(lbl, 0) && GO.is_linear(lbl, 0)
    @test GO.is_line_location_unknown(lbl, 0)
    GO.set_location_collapse!(lbl, 0)
    @test GO.get_line_location(lbl, 0) == LOC_INTERIOR   # hole collapse -> INTERIOR

    # Line
    lbl = GO._create_label(GO.MergeEdge(Int32(1), Int32(2), Int32(1), Int32(1),
        DIM_L, Int32(0), false, DIM_NOT_PART, Int32(0), false))
    @test GO.is_line(lbl) && GO.is_line(lbl, 0) && GO.is_linear(lbl, 0)
    @test !GO.is_boundary(lbl, 0)

    # NotPart both
    lbl = GO.OverlayLabel()
    @test GO.is_not_part(lbl, 0) && GO.is_not_part(lbl, 1)
    @test !GO.is_known(lbl, 0) && !GO.is_boundary_either(lbl)
end

@testset "getLocation isForward L/R swap consistency" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    g = GO.OverlayGraph(Planar(), GO.NodedArrangement(Planar(), A, B; exact = EX); exact = EX)
    for i in eachindex(g.edges)
        s = GO.he_sym(g.edges, i)
        lbl = g.edges[i].label
        fi = g.edges[i].is_forward; fs = g.edges[s].is_forward
        @test fi != fs                                  # a pair has opposite directions
        for idx in (0, 1)
            #-- forward.L == reverse.R and forward.R == reverse.L, from ONE shared label
            @test GO.get_location(lbl, idx, POS_LEFT,  fi) == GO.get_location(lbl, idx, POS_RIGHT, fs)
            @test GO.get_location(lbl, idx, POS_RIGHT, fi) == GO.get_location(lbl, idx, POS_LEFT,  fs)
            @test GO.get_location(lbl, idx, POS_ON, fi) == GO.get_location(lbl, idx, POS_ON, fs)
        end
        #-- the OverlayEdge accessor agrees with the label + orientation
        @test GO.oe_get_location(g.edges, i, 0, POS_LEFT) == GO.get_location(lbl, 0, POS_LEFT, fi)
    end
end

@testset "result-marking flags and accessors" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    g = GO.OverlayGraph(Planar(), GO.NodedArrangement(Planar(), A, B; exact = EX); exact = EX)
    i = 1; s = GO.he_sym(g.edges, i)
    @test !GO.oe_in_result_area(g.edges, i)
    GO.oe_mark_in_result_area!(g.edges, i)
    @test GO.oe_in_result_area(g.edges, i) && !GO.oe_in_result_area(g.edges, s)
    GO.oe_mark_in_result_area_both!(g.edges, i)
    @test GO.oe_in_result_area_both(g.edges, i)
    GO.oe_unmark_from_result_area_both!(g.edges, i)
    @test !GO.oe_in_result_area(g.edges, i) && !GO.oe_in_result_area(g.edges, s)
    GO.oe_set_next_result_max!(g.edges, i, s)
    @test GO.oe_next_result_max(g.edges, i) == s && GO.oe_is_result_max_linked(g.edges, i)
    GO.oe_mark_in_result_line!(g.edges, i)
    @test GO.oe_in_result_line(g.edges, i) && GO.oe_in_result_line(g.edges, s)
end

# ---------------------------------------------------------------------------
# 4. Graph invariants on small real inputs (env-gated, phase-1 smoke pattern)
# ---------------------------------------------------------------------------

ne_ok = try
    import NaturalEarth, GeoJSON
    include(joinpath(@__DIR__, "..", "..", "..", "data", "natural_earth_pairs.jl"))
    global ne_names, ne_geoms = load_ne(110)
    length(ne_geoms) > 0
catch err
    @info "Natural Earth subset skipped (data unavailable)" err
    false
end

@testset "Natural Earth country-pair graph invariants" begin
    if !ne_ok
        @test_skip "Natural Earth data unavailable"
    else
        picks = String["Brazil", "France", "Egypt", "India", "Australia"]
        tested = 0
        for nm in picks
            idx = findfirst(==(nm), ne_names)
            idx === nothing && continue
            A = ne_geoms[idx]
            B = shift_geom(A, 0.5, 0.0)
            tested += 1
            for m in (Planar(), Spherical())
                arr = GO.NodedArrangement(m, A, B; exact = EX)
                g = GO.OverlayGraph(m, arr; exact = EX)
                #-- every half-edge has a sym; stars are CCW-consistent; degrees
                #-- match arrangement incidence (checks all invariants at once).
                #-- Folded: this loop's job is breadth over real linework, and one
                #-- assertion per property says the same thing as 36,668 of them.
                #--
                #-- Two assertions used to follow and are gone: `iseven(length(g.edges))`
                #-- (half-edges are pushed in sym pairs, so it is true by construction)
                #-- and a per-node `g.edges[e].origin == nid` loop, which is a strict
                #-- subset of the degree/origin checks `check_graph` already runs.
                check_graph(m, g; fold = true)
            end
        end
        @test tested >= 2
    end
end
