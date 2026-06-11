using Test
import GeometryOps as GO
import GeoInterface as GI

# Two overlapping convex squares
square_a = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
square_b = GI.Polygon([[(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0), (5.0, 5.0)]])
# Square disjoint from all of the above
square_far = GI.Polygon([[(20.0, 20.0), (25.0, 20.0), (25.0, 25.0), (20.0, 25.0), (20.0, 20.0)]])
# Concave pair whose intersection is two separate pieces
u_shape = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (7.0, 10.0), (7.0, 3.0), (3.0, 3.0), (3.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
bar = GI.Polygon([[(-1.0, 6.0), (11.0, 6.0), (11.0, 8.0), (-1.0, 8.0), (-1.0, 6.0)]])
# Square with a hole
holed_a = GI.Polygon([
    [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
    [(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0), (4.0, 4.0)],
])
# Covers all of holed_a, so holed_a's hole is fully within the overlap
big_b = GI.Polygon([[(-1.0, -1.0), (11.0, -1.0), (11.0, 11.0), (-1.0, 11.0), (-1.0, -1.0)]])
#= Left half of holed_a - holed_a's hole straddles the clip boundary at x = 5, which forces
the recursive `difference` call inside `_add_holes_to_polys!` =#
half_b = GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
#= Holed square whose hole overlaps holed_a's hole, which forces the recursive `union` call
inside `_combine_holes!` =#
holed_b = GI.Polygon([
    [(-1.0, -1.0), (9.0, -1.0), (9.0, 9.0), (-1.0, 9.0), (-1.0, -1.0)],
    [(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0), (5.0, 5.0)],
])
#= Holed square contained in the exterior of another holed square, with intersecting holes -
forces the recursive `intersection` call inside `_add_union_holes_contained_polys!` =#
inner_holed = GI.Polygon([
    [(2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0), (2.0, 2.0)],
    [(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0), (4.0, 4.0)],
])
outer_holed = GI.Polygon([
    [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
    [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)],
])

clipping_ops = (GO.intersection, GO.difference, GO.union)

function test_cached_equals_uncached(op, a, b; T = Float64, kwargs...)
    alg = GO.FosterHormannClipping()
    cache = GO.FosterHormannCache(T)
    uncached = op(alg, a, b, T; target = GI.PolygonTrait(), kwargs...)
    cached = op(alg, a, b, T; target = GI.PolygonTrait(), cache, kwargs...)
    @test GI.coordinates.(cached) == GI.coordinates.(uncached)
    @test cache.depth == 0
end

@testset "FosterHormannCache" begin
    @testset "Results identical with and without cache" begin
        fixtures = [
            ("convex pair", square_a, square_b),
            ("disjoint pair", square_a, square_far),
            ("concave multi-piece", u_shape, bar),
            ("hole within overlap", holed_a, big_b),
            ("hole straddling clip boundary", holed_a, half_b),
            ("overlapping holes", holed_a, holed_b),
            ("contained polygons with holes", inner_holed, outer_holed),
        ]
        for (name, a, b) in fixtures
            @testset "$name" begin
                for op in clipping_ops
                    test_cached_equals_uncached(op, a, b)
                    test_cached_equals_uncached(op, b, a)
                end
            end
        end
    end

    @testset "MultiPolygon inputs" begin
        multi_a = GI.MultiPolygon([square_a, square_far])
        for op in clipping_ops
            test_cached_equals_uncached(op, multi_a, square_b)
            test_cached_equals_uncached(op, square_b, multi_a)
            test_cached_equals_uncached(op, multi_a, square_b; fix_multipoly = nothing)
        end
    end

    @testset "Cache reuse across many varied calls" begin
        alg = GO.FosterHormannClipping()
        cache = GO.FosterHormannCache()
        pairs = [(square_a, square_b), (holed_a, holed_b), (u_shape, bar), (holed_a, half_b), (square_a, square_far), (inner_holed, outer_holed)]
        for (a, b) in pairs, op in clipping_ops
            expected = op(alg, a, b; target = GI.PolygonTrait())
            result = op(alg, a, b; target = GI.PolygonTrait(), cache)
            @test GI.coordinates.(result) == GI.coordinates.(expected)
            @test cache.depth == 0
        end
    end

    @testset "Results don't alias the cache" begin
        alg = GO.FosterHormannClipping()
        cache = GO.FosterHormannCache()
        r1 = GO.intersection(alg, square_a, square_b; target = GI.PolygonTrait(), cache)
        snapshot = deepcopy(GI.coordinates.(r1))
        GO.intersection(alg, u_shape, bar; target = GI.PolygonTrait(), cache)
        GO.intersection(alg, holed_a, holed_b; target = GI.PolygonTrait(), cache)
        @test GI.coordinates.(r1) == snapshot
    end

    @testset "Mismatched cache float type throws" begin
        alg = GO.FosterHormannClipping()
        cache32 = GO.FosterHormannCache(Float32)
        @test_throws ArgumentError GO.intersection(alg, square_a, square_b; target = GI.PolygonTrait(), cache = cache32)
        @test_throws ArgumentError GO.difference(alg, square_a, square_b; target = GI.PolygonTrait(), cache = cache32)
        @test_throws ArgumentError GO.union(alg, square_a, square_b; target = GI.PolygonTrait(), cache = cache32)
        # but a matching Float32 cache works
        test_cached_equals_uncached(GO.intersection, square_a, square_b; T = Float32)
    end

    @testset "reset!" begin
        alg = GO.FosterHormannClipping()
        cache = GO.FosterHormannCache()
        #= holed_a's hole straddles half_b's boundary, forcing a recursive difference call,
        so this grows the frame stack beyond the initial single frame =#
        GO.difference(alg, holed_a, half_b; target = GI.PolygonTrait(), cache)
        @test length(cache.frames) > 1
        @test !isempty(cache.frames[1].a_list)
        @test GO.reset!(cache) === cache
        @test cache.depth == 0
        @test length(cache.frames) == 1
        frame = cache.frames[1]
        @test isempty(frame.a_list) && isempty(frame.b_list) && isempty(frame.a_idx_list) &&
            isempty(frame.remove_poly_idx) && isempty(frame.remove_hole_idx)
        # The cache still works correctly after a reset
        expected = GO.difference(alg, holed_a, half_b; target = GI.PolygonTrait())
        result = GO.difference(alg, holed_a, half_b; target = GI.PolygonTrait(), cache)
        @test GI.coordinates.(result) == GI.coordinates.(expected)
        @test cache.depth == 0
    end

    @testset "Cached calls allocate less" begin
        alg = GO.FosterHormannClipping()
        cache = GO.FosterHormannCache()
        run_cached(alg, a, b, cache) = GO.intersection(alg, a, b; target = GI.PolygonTrait(), cache)
        run_uncached(alg, a, b) = GO.intersection(alg, a, b; target = GI.PolygonTrait())
        # warm up both paths, then measure
        run_cached(alg, square_a, square_b, cache)
        run_uncached(alg, square_a, square_b)
        cached_allocs = @allocated run_cached(alg, square_a, square_b, cache)
        uncached_allocs = @allocated run_uncached(alg, square_a, square_b)
        @test cached_allocs < uncached_allocs
    end
end
