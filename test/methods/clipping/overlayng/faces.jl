# # Face enumeration (Layer A of the antimeridian split)
#
# Tests for `_build_faces` / `_face_ring_location` / `_build_face_polygons` in
# `polygon_builder.jl`: the non-dissolving companion to the op pipeline that
# enumerates every minimal ring (face) of the noded arrangement, and the
# predicate-dispatch overload of `_is_result_of_op` in `overlay_labeller.jl`.

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
@testset "one-extraction-per-graph contract" begin
    # Face extraction WRITES the ring-linkage fields (`next_result`/`edge_ring`)
    # of the shared graph. A second `_build_faces` on the same graph therefore
    # finds every edge already ring-tagged and enumerates nothing — the graph is
    # consumed. Callers must build one graph per extraction.
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    g = labelled_graph(PL, A, B, 2, 2)
    ctx1 = GO._build_faces(PL, g; exact = EX)
    @test length(ctx1.edge_rings) == 4
    ctx2 = GO._build_faces(PL, g; exact = EX)           # same graph, already consumed
    @test length(ctx2.edge_rings) == 0
end
