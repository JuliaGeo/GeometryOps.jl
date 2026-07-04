using Test

import GeoInterface as GI
import GeometryOps as GO
import GeometryOps.FlexibleRTrees as FRT
import GeometryOps.FlexibleRTrees: RTree, STR, HPR, Unsorted, query, hilbert_key
import GeometryOps.SpatialTreeInterface as STI
import Extents
using Random: Xoshiro

# Random boxes with side lengths ~5% of the unit cube, in N dims.
function random_extents(rng, n, N)
    dims = (:X, :Y, :Z, :M)[1:N]
    return [begin
        lo = ntuple(_ -> rand(rng), N)
        hi = lo .+ 0.05 .* ntuple(_ -> rand(rng), N)
        Extents.Extent(NamedTuple{dims}(tuple.(lo, hi)))
    end for _ in 1:n]
end

brute_force(ext, extents) = findall(e -> Extents.intersects(ext, e), extents)

grow(ext, d) = Extents.buffer(ext, NamedTuple{keys(ext)}(ntuple(_ -> d, length(keys(ext)))))

@testset "query ≡ brute force ($(N)D, $alg, n = $n)" for
        N in (1, 2, 3, 4),
        alg in (STR(), HPR(), Unsorted()),
        n in (1, 5, 16, 17, 100, 256, 1000)
    rng = Xoshiro(hash((N, n)))
    extents = random_extents(rng, n, N)
    tree = RTree(alg, extents; nodecapacity = 8)
    queries = vcat(
        random_extents(rng, 20, N),                          # small probes
        [reduce(Extents.union, extents)],                    # everything
        [grow(reduce(Extents.union, extents), 10.0)],        # superset
        [grow(e, 3.0) for e in random_extents(Xoshiro(0), 5, N)],  # big probes
    )
    for q in queries
        @test query(tree, q) == brute_force(q, extents)
    end
    # A query far outside everything returns nothing.
    faraway = Extents.Extent(NamedTuple{((:X, :Y, :Z, :M)[1:N])}(ntuple(_ -> (99.0, 100.0), N)))
    @test isempty(query(tree, faraway))
end

@testset "construction and type stability" begin
    rng = Xoshiro(7)
    extents = random_extents(rng, 300, 2)
    tree = @inferred RTree(STR(), extents)
    @test tree isa RTree{STR, eltype(extents)}
    @inferred RTree(HPR(), extents)
    @inferred RTree(Unsorted(), extents)
    # The query path is inferrable too (depth_first_search returns Vector{Int}).
    q = Extents.Extent(X = (0.2, 0.4), Y = (0.2, 0.4))
    @inferred query(tree, q)
    # Deep and shallow trees have the SAME concrete type — the point of the flat layout.
    tiny = RTree(STR(), extents[1:3])
    @test typeof(tiny) === typeof(tree)

    # Construction and query are just as inferrable in 3D.
    extents3 = random_extents(rng, 100, 3)
    tree3 = @inferred RTree(STR(), extents3)
    @inferred RTree(HPR(), extents3)
    @inferred query(tree3, Extents.Extent(X = (0.2, 0.4), Y = (0.2, 0.4), Z = (0.2, 0.4)))

    # Geometries as input work through GI.extent.
    lines = GO.to_edgelist(GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]))
    gtree = RTree(HPR(), lines)
    # Edge 2 is the right side; edge 3 (the closing diagonal) has a bbox
    # covering the whole square, so extent-intersects finds it too.
    @test query(gtree, Extents.Extent(X = (0.9, 1.1), Y = (0.4, 0.6))) == [2, 3]

    @test_throws ArgumentError RTree(STR(), Extents.Extent{(:X, :Y)}[])
    @test_throws ArgumentError RTree(STR(), random_extents(rng, 5, 2); nodecapacity = 1)
    @test occursin("RTree{HPR}", sprint(show, gtree))
end

@testset "Hilbert curve properties" begin
    # Order-1 2D curve: the classic U through the four quadrants.
    keys1 = [hilbert_key((UInt32(x), UInt32(y)), 1) for (x, y) in ((0, 0), (0, 1), (1, 1), (1, 0))]
    @test keys1 == [0, 1, 2, 3]
    # In any dimension: the curve visits every grid cell exactly once
    # (bijectivity), and consecutive cells are adjacent (unit Manhattan step).
    for (N, bits) in ((1, 6), (2, 4), (3, 2), (4, 2))
        side = 2^bits
        cells = vec(collect(Iterators.product(ntuple(_ -> 0:(side - 1), N)...)))
        keys = [hilbert_key(UInt32.(c), bits) for c in cells]
        @test allunique(keys)
        path = cells[sortperm(keys)]
        @test all(sum(abs.(path[i + 1] .- path[i])) == 1 for i in 1:(length(path) - 1))
    end
end

@testset "NaturalIndex is dimension-agnostic too" begin
    rng = Xoshiro(11)
    extents = random_extents(rng, 300, 3)
    idx = GO.NaturalIndexing.NaturalIndex(extents)
    for q in random_extents(rng, 20, 3)
        @test STI.query(idx, q) == brute_force(q, extents)
    end
end

@testset "prepared geometry integration" begin
    poly = GI.Polygon([
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
        [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)],
    ])
    for alg in (STR(), HPR(), Unsorted())
        et = GO.EdgeTree(GI.getexterior(poly); backend = alg)   # via the build_edge_tree hook
        @test GO.edge_tree(et) isa RTree
        prep = GO.prepare(poly; preps = GO.EdgeTrees(alg))
        @test GO.edge_tree(GO.getprep(GI.getexterior(prep), GO.AbstractEdgeTree)) isa RTree
        for x in 0.0:0.5:11.0, y in (0.0, 3.0, 5.0, 6.5, 10.0)
            pt = (x, y)
            @test GO.within(pt, prep) == GO.within(pt, poly)
            @test GO.intersects(pt, prep) == GO.intersects(pt, poly)
        end
    end
end
