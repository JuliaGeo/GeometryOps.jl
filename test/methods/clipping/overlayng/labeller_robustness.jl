# Labeller robustness on near-coincident inputs (overlay_labeller.jl F0/F1/F2).
#
# Two valid geometries whose boundaries agree to within ±1 ulp still meet in a
# genuine, hair-width sliver, and the exact arrangement resolves it as such. What
# used to break was the LABELLER: it produced an in-result edge set in which some
# node had more result-area boundary arriving than leaving, and the failure only
# surfaced several steps later as
# `_OverlayTopologyError("Ring edge missing at max-ring build")` in
# maximal_edge_ring.jl — a symptom with no pointer to its cause.
#
# This file pins the three pieces that fixed it, in the order they run:
#
#   F0  after result-area marking, in/out degree balance at every node is a hard
#       invariant, checked on every overlay (`_check_result_area_balance`).
#   F1  pass 5 queries the point-in-area locator with the node's KERNEL point,
#       never its emitted output coordinate.
#   F2  a known area location for the other input propagates across nodes that
#       carry no boundary edge of it — where it provably cannot change — so the
#       per-edge point-in-area verdict only ever seeds a whole chain at once.
#
# The synthetic sweep at the bottom is the reduced reproducer: two rectangles
# sharing an edge, one of them with its interior shared-edge vertices moved by
# ±k ulp of latitude. Both inputs are valid and the planar answer is exact and
# clean at every size; spherically the lon/lat -> unit-vector ingest scatters the
# two nominally-identical meridian polylines by ~1 ulp, which is all it takes.
#
# Output VALIDITY is deliberately NOT asserted here. At these separations the
# exact union genuinely contains sliver faces, and how emission should render
# them is a separate defect (`_ring_is_collapsed`, xml_suite.jl's jts-798 case).
# What this file asserts is that the engine does not THROW on valid input and
# that the graph it hands the builder is consistent.

using Test
include(joinpath(@__DIR__, "common.jl"))

const PL = GO.Planar()
const SPH = GO.Spherical()

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# A = rectangle (0,49)-(1,50) whose shared edge lon=1 is subdivided at `lats`.
# B = rectangle (1,49)-(2,50) whose shared edge uses the same latitudes moved by
# ±`k` ulp, alternating, with the two ENDPOINTS left exact so the two rings still
# share their corner vertices bit-for-bit.
function sliver_pair(n::Int, k::Int)
    lats = collect(range(49.0, 50.0; length = n + 2))
    bump(l, up) = (x = l; for _ in 1:k; x = up ? nextfloat(x) : prevfloat(x) end; x)
    pert = [(1.0, (i == 1 || i == length(lats)) ? l : bump(l, isodd(i)))
            for (i, l) in enumerate(lats)]
    a_ring = [(0.0, 49.0); [(1.0, l) for l in lats]; (0.0, 50.0); (0.0, 49.0)]
    b_ring = [pert; (2.0, 50.0); (2.0, 49.0); first(pert)]
    return GI.Polygon([a_ring]), GI.Polygon([b_ring])
end

# A labelled, result-marked graph, stopping just short of ring extraction.
function marked_graph(m, A, B, op)
    arr = GO.NodedArrangement(m, A, B; exact = EX)
    g = GO.OverlayGraph(m, arr; exact = EX)
    input = GO._OverlayInput(m, A, B, 2, 2, EX, false, false, nothing, nothing)
    GO._compute_labelling!(g, input)
    GO._mark_result_area_edges!(g, op)
    GO._unmark_duplicate_edges_from_result_area!(g)     #-- runs F0
    return g
end

nrings(g) = (t = GI.trait(g);
             t isa GI.PolygonTrait ? GI.nring(g) :
             t isa GI.MultiPolygonTrait ? sum(GI.nring(p) for p in GI.getgeom(g); init = 0) : 0)

# `nothing`, or the exception the overlay raised.
function overlay_error(m, op, A, B)
    try
        GO._overlay_ng(m, op, A, B; exact = EX)
        return nothing
    catch e
        return e
    end
end

# ---------------------------------------------------------------------------
# F1 — the point-in-area query runs on kernel points
# ---------------------------------------------------------------------------

@testset "F1: pass 5 locates on kernel points, not emitted coordinates" begin
    A, B = sliver_pair(7, 1)

    #-- PLANAR: the query point is unchanged, node for node. Vertex nodes key by
    #-- their coordinate, and a planar crossing's emitted coordinate is the
    #-- certified correctly-rounded image of the exact rational one — which is
    #-- the best `Tuple{Float64,Float64}` there is, and all the planar locator
    #-- can accept. So F1 is a no-op on the plane, by construction.
    arr_p = GO.NodedArrangement(PL, A, B; exact = EX)
    @test all(GO._node_kernel_point(arr_p, i) === GO.node_point(arr_p, i)
              for i in eachindex(arr_p.nodes.keys))

    #-- SPHERICAL: a vertex node's kernel point is the ingested unit vector,
    #-- bit-for-bit — the same value the locator's own rings are built from.
    arr_s = GO.NodedArrangement(SPH, A, B; exact = EX)
    ring_pts = Set(GO._to_kernel_point(SPH, p) for g in (A, B)
                   for p in GI.getpoint(GI.getexterior(g)))
    vertex_ids = [i for i in eachindex(arr_s.nodes.keys) if !arr_s.nodes.keys[i].is_crossing]
    @test !isempty(vertex_ids)
    @test all(GO._node_kernel_point(arr_s, i) isa GO.UnitSphericalPoint for i in vertex_ids)
    @test all(GO._node_kernel_point(arr_s, i) in ring_pts for i in vertex_ids)

    #-- ...and that is NOT what the emitted coordinate round-trips back to: the
    #-- emission runs `atan`/`asin` and the locator would run `cos`/`sin` on the
    #-- way back in, which is exactly the loss F1 removes.
    round_tripped = [GO._to_kernel_point(SPH, GO.node_point(arr_s, i)) for i in vertex_ids]
    @test any(round_tripped[j] != GO._node_kernel_point(arr_s, vertex_ids[j])
              for j in eachindex(vertex_ids))

    #-- a crossing node's kernel point is the exact on-arc crossing direction,
    #-- normalized — a unit vector, and on the crossing's own great circles
    for i in eachindex(arr_s.nodes.keys)
        arr_s.nodes.keys[i].is_crossing || continue
        p = GO._node_kernel_point(arr_s, i)
        @test p isa GO.UnitSphericalPoint
        @test isapprox(sqrt(sum(abs2, p)), 1.0; atol = 1e-15)
    end
end

# ---------------------------------------------------------------------------
# F2 — one area location per node that cannot separate the other input
# ---------------------------------------------------------------------------

# At a node carrying no boundary edge of input `gi`, `gi`'s area location cannot
# change, so every incident edge that is not part of `gi` must report the same
# `gi` location. This is the property F2 establishes and the property F0's
# balance invariant rests on. (Truncated nodes are excluded: clip pruning thinned
# their star, so a `gi`-boundary edge may simply be absent from it — the same
# exclusion pass 1 makes.)
function transparent_nodes_are_consistent(g, gi)
    edges = g.edges
    truncated = g.arr.truncated
    for nid in eachindex(g.node_edges)
        ne = g.node_edges[nid]
        ne == 0 && continue
        (!isempty(truncated) && truncated[nid]) && continue
        GO._find_propagation_start_edge(edges, ne, gi) == 0 || continue
        loc = nothing
        e = ne
        while true
            l = GO.oe_label(edges, e)
            if GO.is_not_part(l, gi)
                this = GO.get_line_location(l, gi)
                loc === nothing ? (loc = this) : (loc == this || return false)
            end
            e = GO.he_onext(edges, e)
            e == ne && break
        end
    end
    return true
end

@testset "F2: an area location is constant across nodes that cannot change it" begin
    #-- the canonical overlapping squares, both manifolds and all four ops
    for m in (PL, SPH), op in OP_CODES
        g = marked_graph(m, SQ_A, SQ_B, op)
        @test transparent_nodes_are_consistent(g, 0)
        @test transparent_nodes_are_consistent(g, 1)
    end

    #-- and on the sliver pair, where the per-edge point-in-area verdicts that F2
    #-- replaces are coin flips
    for n in (3, 7, 15), k in (1, 2)
        A, B = sliver_pair(n, k)
        for m in (PL, SPH)
            g = marked_graph(m, A, B, GO.OVERLAY_UNION)
            @test transparent_nodes_are_consistent(g, 0)
            @test transparent_nodes_are_consistent(g, 1)
        end
    end
end

# ---------------------------------------------------------------------------
# F0 — the result-area degree-balance invariant
# ---------------------------------------------------------------------------

@testset "F0: result-area in/out degree balance holds on well-behaved overlays" begin
    #-- `_unmark_duplicate_edges_from_result_area!` runs the check itself, so a
    #-- silent return is the assertion; call it explicitly too, for the record
    for m in (PL, SPH), op in OP_CODES
        g = marked_graph(m, SQ_A, SQ_B, op)
        @test GO._check_result_area_balance(g) === nothing
    end

    #-- a hole, a touch, and a multipolygon, on both manifolds
    holed = GI.Polygon([[(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0), (0.0, 0.0)],
                        [(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0), (2.0, 2.0)]])
    touch = GI.Polygon([[(6.0, 0.0), (9.0, 0.0), (9.0, 6.0), (6.0, 6.0), (6.0, 0.0)]])
    multi = GI.MultiPolygon([[[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]],
                             [[(3.0, 3.0), (5.0, 3.0), (5.0, 5.0), (3.0, 5.0), (3.0, 3.0)]]])
    for m in (PL, SPH), op in OP_CODES, B in (touch, multi)
        @test GO._check_result_area_balance(marked_graph(m, holed, B, op)) === nothing
    end
end

@testset "F0: an imbalance is reported at the node, not downstream" begin
    #-- the diagnostic names the node id, its emitted coordinate and its whole
    #-- star, so an imbalance is localized from the message alone
    g = marked_graph(PL, SQ_A, SQ_B, GO.OVERLAY_UNION)
    #-- forge one: unmark a single half-edge that is in the result area
    i = findfirst(j -> GO.oe_in_result_area(g.edges, j), eachindex(g.edges))
    @test i !== nothing
    g.edges[i].in_result_area = false
    err = try
        GO._check_result_area_balance(g)
        nothing
    catch e
        e
    end
    @test err isa GO._OverlayTopologyError
    @test occursin("result-area degree imbalance at node", err.msg)
    @test occursin("star (edge => dest, out/in):", err.msg)
end

# ---------------------------------------------------------------------------
# The reduced sliver sweep
# ---------------------------------------------------------------------------

#=
This sweep needed one fix outside the labeller as well, and the two are
independent — measured, not assumed (throws out of 12):

    JTS pass 5    + shipped comparator   5      F0+F1+F2 + shipped comparator   3
    JTS pass 5    + fixed comparator     4      F0+F1+F2 + fixed comparator     0

The other half was `rk_compare_along_segment(::Spherical, …)`, whose float filter
derived its tolerance from the very float crossing directions whose accuracy was
the problem. On `n = 7, k = 1` it ordered the two crossings of A's segment
[49.25, 49.375] backwards, which made pass 1's own side locations mutually
contradictory across a chain of nodes carrying no B edge at all — a state no
labelling rule can repair, because the location genuinely cannot change there.
That filter is pinned directly in `test/methods/relateng/kernel.jl`.
=#

@testset "sliver sweep: valid near-coincident inputs never reach the ring builder" begin
    println("sliver sweep (spherical union) — n, ulp, status, rings:")
    for n in (1, 3, 7, 15, 31, 63), k in (1, 2)
        A, B = sliver_pair(n, k)
        @test LG.isValid(GI.convert(LG, A))
        @test LG.isValid(GI.convert(LG, B))

        #-- planar: exactly collinear shared edges, one clean rectangle, always
        rp = GO._overlay_ng(PL, GO.OVERLAY_UNION, A, B; exact = EX)
        @test nrings(rp) == 1
        @test isapprox(GO.area(rp), 2.0; rtol = 1e-12)

        #-- spherical: no throw anywhere, and the F0 balance invariant holds
        #-- (`_unmark_duplicate_edges_from_result_area!` checks it inside the run)
        @test overlay_error(SPH, GO.OVERLAY_UNION, A, B) === nothing
        r = GO._overlay_ng(SPH, GO.OVERLAY_UNION, A, B; exact = EX)
        #-- the area is right to machine precision even where the sliver faces
        #-- survive as extra rings
        @test isapprox(GO.area(SPH, r), GO.area(SPH, A) + GO.area(SPH, B); rtol = 1e-12)
        println("  n=$n k=$k  ok  rings=$(nrings(r))")
    end
end
