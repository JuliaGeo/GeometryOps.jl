# # Face enumeration (design doc §5; Layer A of the antimeridian split)
#
# Tests for `_build_faces` / `_face_ring_location` / `_select_faces!` /
# `_build_face_polygons` in `polygon_builder.jl`: the non-dissolving companion to
# the op pipeline that enumerates every boundary cycle of the noded arrangement,
# its one-extraction-per-graph guard, its opt-in dangle / cut-edge hygiene pass,
# and the predicate-dispatch overload of `_is_result_of_op` in
# `overlay_labeller.jl`.

using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

const EX = GO.True()
const PL = GO.Planar()

locname(l) = l == GO.LOC_INTERIOR ? "INT" : l == GO.LOC_EXTERIOR ? "EXT" :
             l == GO.LOC_BOUNDARY ? "BND" : "NONE"

# Signed planar shoelace area: < 0 for clockwise rings, > 0 for counter-clockwise.
signed_ring_area(pts) = 0.5 * sum(pts[i][1] * pts[i + 1][2] - pts[i + 1][1] * pts[i][2]
                                  for i in 1:(length(pts) - 1))

# Build a labelled graph + a face context for A (dim `dim_a`) against B (dim `dim_b`).
function faces_of(m, A, B, dim_a, dim_b; exact = EX)
    arr = GO.NodedArrangement(m, A, B; exact)
    g = GO.OverlayGraph(m, arr; exact)
    input = GO._OverlayInput(m, A, B, dim_a, dim_b, exact, GI.isempty(A), GI.isempty(B),
                             nothing, nothing)
    GO._compute_labelling!(g, input)
    return g, GO._build_faces(m, g; exact)
end

# A freshly labelled graph, for callers that build their own face context.
function labelled_graph(m, A, B, dim_a, dim_b; exact = EX)
    arr = GO.NodedArrangement(m, A, B; exact)
    g = GO.OverlayGraph(m, arr; exact)
    GO._compute_labelling!(g, GO._OverlayInput(m, A, B, dim_a, dim_b, exact,
                                               GI.isempty(A), GI.isempty(B), nothing, nothing))
    return g
end

# ---------------------------------------------------------------------------
@testset "_is_result_of_op: predicate == enum on all 4 ops x 16 location pairs" begin
    preds = Dict(
        GO.OVERLAY_INTERSECTION  => (a, b) -> a == GO.LOC_INTERIOR && b == GO.LOC_INTERIOR,
        GO.OVERLAY_UNION         => (a, b) -> a == GO.LOC_INTERIOR || b == GO.LOC_INTERIOR,
        GO.OVERLAY_DIFFERENCE    => (a, b) -> a == GO.LOC_INTERIOR && b != GO.LOC_INTERIOR,
        GO.OVERLAY_SYMDIFFERENCE => (a, b) -> (a == GO.LOC_INTERIOR) ⊻ (b == GO.LOC_INTERIOR),
    )
    locs = (GO.LOC_INTERIOR, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.LOC_NONE)
    for (op, pred) in preds, l0 in locs, l1 in locs
        @test GO._is_result_of_op(op, l0, l1) == GO._is_result_of_op(pred, l0, l1)
    end

    # ...and the loosened signatures make the dissolving pipeline dispatch a
    # predicate identically to its enum op.
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    for (op, pred) in preds
        ge = labelled_graph(PL, A, B, 2, 2); GO._mark_result_area_edges!(ge, op)
        gp = labelled_graph(PL, A, B, 2, 2); GO._mark_result_area_edges!(gp, pred)
        @test [GO.oe_in_result_area(ge.edges, i) for i in eachindex(ge.edges)] ==
              [GO.oe_in_result_area(gp.edges, i) for i in eachindex(gp.edges)]
    end
end

# ---------------------------------------------------------------------------
@testset "two offset squares -> 4 faces (A∩B, A\\B, B\\A, outer)" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    _, ctx = faces_of(PL, A, B, 2, 2)

    nr = length(ctx.edge_rings)
    @test nr == 4

    labels = Dict{Tuple{Int,Int},Any}()
    outers = 0
    for er in 1:nr
        r = ctx.edge_rings[er]
        la = GO._face_ring_location(ctx, er, 0)
        lb = GO._face_ring_location(ctx, er, 1)
        labels[(Int(la), Int(lb))] = r
        if r.is_hole                       # CCW ring == the unbounded outer face
            outers += 1
            @test la == GO.LOC_EXTERIOR && lb == GO.LOC_EXTERIOR
            @test signed_ring_area(r.ring_pts) > 0          # CCW
        else                               # bounded faces are CW shells
            @test signed_ring_area(r.ring_pts) < 0          # CW
        end
    end
    # exactly one outer (EXT,EXT) face, and it is the only CCW ring
    @test outers == 1

    # the three bounded faces carry the expected per-input locations and areas
    INT = Int(GO.LOC_INTERIOR); EXT = Int(GO.LOC_EXTERIOR)
    @test haskey(labels, (INT, INT))   # A∩B
    @test haskey(labels, (INT, EXT))   # A\B
    @test haskey(labels, (EXT, INT))   # B\A
    @test isapprox(abs(signed_ring_area(labels[(INT, INT)].ring_pts)), 1.0; rtol = 1e-12)
    @test isapprox(abs(signed_ring_area(labels[(INT, EXT)].ring_pts)), 3.0; rtol = 1e-12)
    @test isapprox(abs(signed_ring_area(labels[(EXT, INT)].ring_pts)), 3.0; rtol = 1e-12)

    # face selection via `_build_face_polygons` (fresh graph each time — see the
    # one-extraction-per-graph contract below).
    p_int = GO._build_face_polygons(PL, labelled_graph(PL, A, B, 2, 2),
                                    (a, b) -> a == GO.LOC_INTERIOR && b == GO.LOC_INTERIOR; exact = EX)
    @test length(p_int) == 1
    @test isapprox(GO.area(p_int[1]), 1.0; rtol = 1e-12)

    p_a = GO._build_face_polygons(PL, labelled_graph(PL, A, B, 2, 2),
                                  (a, b) -> a == GO.LOC_INTERIOR; exact = EX)  # A∩B ∪ A\B
    @test length(p_a) == 2
    @test isapprox(sum(GO.area, p_a), 4.0; rtol = 1e-12)
end

# ---------------------------------------------------------------------------
@testset "polygon-with-hole -> 2-ring face polygon, valid" begin
    Ah = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                     [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]])
    Bh = GI.Polygon([[(20.0, 0.0), (21.0, 0.0), (21.0, 1.0), (20.0, 1.0), (20.0, 0.0)]])  # disjoint
    polys = GO._build_face_polygons(PL, labelled_graph(PL, Ah, Bh, 2, 2),
                                    (a, b) -> a == GO.LOC_INTERIOR; exact = EX)
    @test length(polys) == 1
    @test GI.nring(polys[1]) == 2                       # shell + cavity assigned
    @test isapprox(GO.area(polys[1]), 84.0; rtol = 1e-12)   # 100 - 16
    @test LG.isValid(GI.convert(LG, polys[1]))
end

# ---------------------------------------------------------------------------
@testset "dangling line edge doubles in its face ring" begin
    poly = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    line = GI.LineString([(5.0, -5.0), (5.0, 5.0)])   # enters at (5,0), dead-ends at (5,5)
    _, ctx = faces_of(PL, poly, line, 2, 1)

    kept = [er for er in 1:length(ctx.edge_rings)
            if GO._face_ring_location(ctx, er, 0) == GO.LOC_INTERIOR]
    @test length(kept) == 1                            # the square, with a doubled dangle
    r = ctx.edge_rings[kept[1]]
    @test !r.is_hole
    @test isapprox(abs(signed_ring_area(r.ring_pts)), 100.0; rtol = 1e-12)
    # the dangle tip (5,5) is an interior vertex only reachable by walking the
    # dead-end out and back — its presence proves the doubling.
    @test (5.0, 5.0) in r.ring_pts
    @test length(r.ring_pts) == 8                      # 5-pt square + (5,0),(5,5),(5,0) detour
end

# ---------------------------------------------------------------------------
@testset "one-extraction-per-graph contract is enforced, in both orders" begin
    # Ring extraction WRITES the ring-linkage fields (`next_result` /
    # `next_result_max` / `edge_ring` / `max_edge_ring`) of the shared graph and
    # never resets them, so only ONE extraction per graph is valid. Left
    # unchecked this produces silent wrong answers, not errors (measured:
    # op-then-faces gave 3 of 4 rings, faces-then-op gave 0 of 1 polygons), so
    # both entry points detect an already-consumed graph and throw.
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

    op_extract!(g) = GO._build_polygons(PL, g, GO.graph_result_area_edges(g); exact = EX)
    function marked_graph()
        g = labelled_graph(PL, A, B, 2, 2)
        GO._mark_result_area_edges!(g, GO.OVERLAY_INTERSECTION)
        GO._unmark_duplicate_edges_from_result_area!(g)
        return g
    end

    # a fresh graph extracts fine, either way
    @test length(GO._build_faces(PL, labelled_graph(PL, A, B, 2, 2); exact = EX).edge_rings) == 4
    @test length(op_extract!(marked_graph())) == 1

    # faces THEN faces
    g = labelled_graph(PL, A, B, 2, 2)
    GO._build_faces(PL, g; exact = EX)
    @test_throws GO._OverlayTopologyError GO._build_faces(PL, g; exact = EX)

    # faces THEN op (would have silently produced 0 polygons)
    g = marked_graph()
    GO._build_faces(PL, g; exact = EX)
    @test_throws GO._OverlayTopologyError op_extract!(g)

    # op THEN faces (would have silently produced 3 of the 4 rings)
    g = marked_graph()
    op_extract!(g)
    @test_throws GO._OverlayTopologyError GO._build_faces(PL, g; exact = EX)

    # the hygiene pass is not an extraction: it writes no linkage, so it leaves
    # the graph extractable (and the message names the real problem when it is not)
    g = labelled_graph(PL, A, B, 2, 2)
    @test !GO._graph_is_ring_consumed(g)
    GO._remove_dangles!(g)
    @test !GO._graph_is_ring_consumed(g)
    @test length(GO._build_faces(PL, g; exact = EX).edge_rings) == 4
    @test_throws GO._OverlayTopologyError GO._remove_dangles!(g)   # now consumed

    err = try; GO._build_faces(PL, g; exact = EX); catch e; e; end
    @test occursin("already been consumed", sprint(showerror, err))
end

# ---------------------------------------------------------------------------
# ## Face-walk hygiene — dangle and cut-edge removal (opt-in, default off)

# linework with every hygiene case at once: two squares, an interior divider
# splitting the first, a dangling spur off the divider, a bridge (cut edge)
# joining the two squares, and an isolated segment.
const HYG_LINEWORK = GI.MultiLineString([
    [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],    # square 1
    [(5.0, 0.0), (5.0, 10.0)],                                           # divider
    [(5.0, 5.0), (8.0, 5.0)],                                            # dangling spur
    [(20.0, 0.0), (30.0, 0.0), (30.0, 10.0), (20.0, 10.0), (20.0, 0.0)], # square 2
    [(10.0, 5.0), (20.0, 5.0)],                                          # bridge / cut edge
])
const HYG_FAR = GI.LineString([(100.0, 100.0), (101.0, 101.0)])          # isolated segment

hyg_faces(; kw...) = GO._build_faces(PL, labelled_graph(PL, HYG_LINEWORK, HYG_FAR, 1, 1);
                                     exact = EX, kw...)

# the bounded (CW) rings' areas, rounded for comparison
bounded_areas(ctx) = sort([round(abs(signed_ring_area(r.ring_pts)); digits = 9)
                           for r in ctx.edge_rings if !r.is_hole && length(r.ring_pts) >= 4])

@testset "hygiene defaults to off — the raw cycle structure is unchanged" begin
    ctx = hyg_faces()
    @test length(ctx.edge_rings) == 5
    # the spur is traversed out and back inside its face's cycle...
    spur = only(r for r in ctx.edge_rings if (8.0, 5.0) in r.ring_pts)
    @test count(==((5.0, 5.0)), spur.ring_pts) == 2
    # ...and the bridge is doubled in the outer cycle, which therefore spans both
    # squares (a single 15-point cycle of total |area| 200)
    outer = only(r for r in ctx.edge_rings if r.is_hole)
    @test length(outer.ring_pts) == 15
    @test isapprox(abs(signed_ring_area(outer.ring_pts)), 200.0; rtol = 1e-12)
    @test count(==((10.0, 5.0)), outer.ring_pts) == 2
    # the isolated segment's cycle is a degenerate 3-point out-and-back
    @test any(r -> length(r.ring_pts) == 3, ctx.edge_rings)
end

@testset "remove_dangles peels degree-1 chains, leaving cut edges alone" begin
    ctx = hyg_faces(; remove_dangles = true)
    # spur + isolated segment gone; the divided square, its other half, square 2,
    # and the (still bridged) outer cycle remain
    @test length(ctx.edge_rings) == 4
    @test !any(r -> (8.0, 5.0) in r.ring_pts, ctx.edge_rings)
    @test !any(r -> (100.0, 100.0) in r.ring_pts, ctx.edge_rings)
    @test bounded_areas(ctx) == [50.0, 50.0, 100.0]
    # the bridge survives: the outer cycle still spans both squares and doubles it
    outer = only(r for r in ctx.edge_rings if r.is_hole)
    @test count(==((10.0, 5.0)), outer.ring_pts) == 2
end

@testset "remove_cut_edges removes bridges (and subsumes dangles)" begin
    for kw in ((remove_cut_edges = true,), (remove_dangles = true, remove_cut_edges = true))
        ctx = hyg_faces(; kw...)
        # the bridge is gone, so the two squares are separate components: each
        # contributes its own CCW outer cycle
        @test length(ctx.edge_rings) == 5
        @test count(r -> r.is_hole, ctx.edge_rings) == 2
        @test bounded_areas(ctx) == [50.0, 50.0, 100.0]
        # no ring visits any node twice any more — the rings are hygienic
        for r in ctx.edge_rings
            @test length(unique(r.ring_pts)) == length(r.ring_pts) - 1
        end
        @test !any(r -> (8.0, 5.0) in r.ring_pts, ctx.edge_rings)
    end
end

@testset "dangle removal is iterative (whole trees peel off)" begin
    tree = GI.MultiLineString([
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
        [(5.0, 0.0), (5.0, 3.0), (7.0, 4.0)],   # 2-edge chain into the interior
        [(5.0, 3.0), (3.0, 4.0)],               # a branch off the chain
    ])
    g0 = labelled_graph(PL, tree, HYG_FAR, 1, 1)
    ctx0 = GO._build_faces(PL, g0; exact = EX)
    @test length(ctx0.edge_rings) == 3                       # incl. the isolated segment
    @test any(r -> (7.0, 4.0) in r.ring_pts, ctx0.edge_rings)

    g = labelled_graph(PL, tree, HYG_FAR, 1, 1)
    removed = GO._remove_dangles!(g)
    @test length(removed) == 4      # 2 chain edges + branch + the isolated segment
    ctx = GO._build_faces(PL, g; exact = EX)
    @test length(ctx.edge_rings) == 2
    @test bounded_areas(ctx) == [100.0]
    sq = only(r for r in ctx.edge_rings if !r.is_hole)
    @test sq.ring_pts == [(5.0, 0.0), (0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (5.0, 0.0)]
end

@testset "hygiene removes the dangle a face ring would otherwise double" begin
    # the `dangling line edge doubles in its face ring` geometry above, with the
    # dangle removed at the graph level instead of by post-filtering points
    poly = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    line = GI.LineString([(5.0, -5.0), (5.0, 5.0)])
    for kw in ((remove_dangles = true,), (remove_cut_edges = true,))
        g = labelled_graph(PL, poly, line, 2, 1)
        ctx = GO._build_faces(PL, g; exact = EX, kw...)
        @test length(ctx.edge_rings) == 2
        kept = only(er for er in 1:length(ctx.edge_rings)
                    if GO._face_ring_location(ctx, er, 0) == GO.LOC_INTERIOR)
        r = ctx.edge_rings[kept]
        @test !r.is_hole
        @test isapprox(abs(signed_ring_area(r.ring_pts)), 100.0; rtol = 1e-12)
        @test !((5.0, 5.0) in r.ring_pts)            # dangle tip gone
        @test length(r.ring_pts) == 6                # 5-pt square + the (5,0) node
        # `_build_face_polygons` forwards the same opt-in
        gp = labelled_graph(PL, poly, line, 2, 1)
        polys = GO._build_face_polygons(PL, gp, (a, b) -> a == GO.LOC_INTERIOR; exact = EX, kw...)
        @test length(polys) == 1
        @test GI.nring(polys[1]) == 1
        @test isapprox(GO.area(polys[1]), 100.0; rtol = 1e-12)
        @test LG.isValid(GI.convert(LG, polys[1]))
    end
end

@testset "hygiene is a face facility — the op pipeline rejects a filtered graph" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    g = labelled_graph(PL, A, B, 2, 2)
    GO._mark_result_area_edges!(g, GO.OVERLAY_INTERSECTION)
    GO._unmark_duplicate_edges_from_result_area!(g)
    @test isempty(GO._remove_cut_edges!(g))        # two overlapping squares: no bridges
    GO.oe_remove_both!(g.edges, 1)                 # so force a removal
    err = try
        GO._build_polygons(PL, g, GO.graph_result_area_edges(g); exact = EX)
    catch e; e; end
    @test err isa GO._OverlayTopologyError
    @test occursin("does not honour removal", sprint(showerror, err))
end

@testset "hygiene on the sphere" begin
    SPH = GO.Spherical()
    poly = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    line = GI.LineString([(5.0, -5.0), (5.0, 5.0)])
    g = labelled_graph(SPH, poly, line, 2, 1)
    ctx = GO._build_faces(SPH, g; exact = EX, remove_dangles = true)
    @test length(ctx.edge_rings) == 2
    kept = only(er for er in 1:length(ctx.edge_rings)
                if GO._face_ring_location(ctx, er, 0) == GO.LOC_INTERIOR)
    @test !((5.0, 5.0) in ctx.edge_rings[kept].ring_pts)
    @test length(ctx.edge_rings[kept].ring_pts) == 6
    @test isapprox(GO.area(SPH, GI.Polygon([ctx.edge_rings[kept].ring_pts])),
                   GO.area(SPH, poly); rtol = 1e-12)
end

# ---------------------------------------------------------------------------
@testset "collapsed face cycles are dropped by the shared selection path" begin
    # A face cycle with fewer than three distinct emitted vertices bounds no area
    # and is not a legal LinearRing. `_select_faces!` — the one path every face
    # consumer runs through — drops it, exactly as `_assign_shells_and_holes!`
    # drops a collapsed ring from the op result.
    A = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
    B = GI.LineString([(5.0, 5.0), (6.0, 6.0)])
    ctx = GO._build_faces(PL, labelled_graph(PL, A, B, 1, 1); exact = EX)
    @test length(ctx.edge_rings) == 2
    @test all(GO._ring_is_collapsed, ctx.edge_rings)     # both are 3-point out-and-backs
    @test all(r -> !r.is_hole, ctx.edge_rings)           # …and would file as shells

    # `keep` accepts everything, yet nothing is emitted
    GO._select_faces!(ctx, (a, b) -> true)
    @test isempty(ctx.shell_list) && isempty(ctx.free_hole_list)
    @test isempty(GO._build_face_polygons(PL, labelled_graph(PL, A, B, 1, 1),
                                          (a, b) -> true; exact = EX))
end
