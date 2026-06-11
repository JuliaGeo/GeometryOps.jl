# Allocation discipline and type-stability checks for RelateNG (Task 28).
#
# The hot path (segment-pair classification inside the edge intersector) must
# not allocate per pair; per-call setup (`RelateGeometry`, the topology
# computer's dicts, segment-string extraction) is unavoidable and scales with
# input size, not with the number of segment *pairs*. So: measure `@allocated`
# for a mid-size polygon pair after warmup and assert it stays within 2x of
# the baseline measured when this test was written. A per-segment-pair
# allocation regression on a 256-vertex pair (~65k candidate pairs) blows
# straight through that budget; honest setup-cost drift does not.

using Test
import GeometryOps as GO
import GeoInterface as GI
using Random

include(joinpath(@__DIR__, "..", "..", "data", "polygon_generation.jl"))

const ALG = GO.RelateNG()

# Deterministic mid-size (256-vertex) polygon pairs. Same generator settings
# as the benchmarks (benchmarks/relateng.jl); validity does not matter for an
# allocation measurement, so no oracle check is needed here.
rng = Xoshiro(42)
poly_a = GO.tuples(GI.Polygon(generate_random_poly(0.0, 0.0, 256, 2.0, 0.3, 0.1, rng)))
poly_b = GO.tuples(GI.Polygon(generate_random_poly(2.0, 0.0, 256, 2.0, 0.3, 0.1, rng)))   # overlaps a
poly_c = GO.tuples(GI.Polygon(generate_random_poly(10.0, 0.0, 256, 2.0, 0.3, 0.1, rng)))  # disjoint from a

@testset "allocation budget (256-vertex pairs)" begin
    # Baselines measured 2026-06-11 (Julia 1.12.6, Apple M4 Pro macOS):
    #   intersects, overlapping pair:  355_232 bytes
    #   intersects, disjoint pair:         112 bytes (extent-filter early exit)
    #   full relate, overlapping pair: 587_600 bytes
    # Budget = 2x measured baseline.
    for (name, x, y, budget) in [
        ("intersects overlapping", poly_a, poly_b, 2 * 355_232),
        ("intersects disjoint",    poly_a, poly_c, 2 * 112),
    ]
        @testset "$name" begin
            GO.relate_predicate(ALG, GO.pred_intersects(), x, y)  # warmup
            allocated = @allocated GO.relate_predicate(ALG, GO.pred_intersects(), x, y)
            @test allocated <= budget
        end
    end
    @testset "full relate overlapping" begin
        GO.relate(ALG, poly_a, poly_b)  # warmup
        allocated = @allocated GO.relate(ALG, poly_a, poly_b)
        @test allocated <= 2 * 587_600
    end
end

@testset "type stability" begin
    m = GO.Planar()
    # Kernel classification: concrete SegSegClass out of tuple points.
    @test (@inferred GO.rk_classify_intersection(m,
        (0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0); exact = GO.True())) isa GO.SegSegClass
    # update_dim! on concrete predicate types.
    @inferred GO.update_dim!(GO.pred_intersects(), GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    @inferred GO.update_dim!(GO.RelateMatrixPredicate(), GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    # Full evaluation entry point returns an inferred Bool.
    @test (@inferred GO.relate_predicate(ALG, GO.pred_intersects(), poly_a, poly_b)) isa Bool
    @test GO.relate_predicate(ALG, GO.pred_intersects(), poly_a, poly_b)
    @test !GO.relate_predicate(ALG, GO.pred_intersects(), poly_a, poly_c)
end
