using Test

import GeoInterface as GI
import GeometryOps as GO
import GeometryOps.FlexibleRTrees as FRT
import GeometryOps.FlexibleRTrees: RTree, STR, HPR, Unsorted, query, hilbert_key
import GeometryOps.SpatialTreeInterface as STI
import Extents
using GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap
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
        N in (2, 3),
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
    for (N, bits) in ((2, 4), (3, 2))
        side = 2^bits
        cells = vec(collect(Iterators.product(ntuple(_ -> 0:(side - 1), N)...)))
        keys = [hilbert_key(UInt32.(c), bits) for c in cells]
        @test allunique(keys)
        path = cells[sortperm(keys)]
        @test all(sum(abs.(path[i + 1] .- path[i])) == 1 for i in 1:(length(path) - 1))
    end
end

@testset "Manifold-aware construction" begin
    band = [GI.Polygon([[(lon, 60.0), (lon + 30.0, 60.0), (lon + 30.0, 80.0), (lon, 80.0), (lon, 60.0)]])
            for lon in 0.0:30.0:330.0]
    cap = GI.Polygon([[(lon, 80.0) for lon in 0.0:30.0:360.0]])  # around the north pole
    geoms = vcat(band, [cap])

    tree = RTree(GO.Spherical(), HPR(), geoms)
    @test Extents.extent(tree) isa Extents.Extent{(:X, :Y, :Z)}

    # only the cap's region reaches the pole; every vertex stops at lat 80
    pole_box = Extents.Extent(X = (-0.01, 0.01), Y = (-0.01, 0.01), Z = (0.99, 1.0))
    @test query(tree, pole_box) == [13]

    # a SphericalCap query against the 3D leaf boxes, end to end
    polecap = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), 0.05)
    @test STI.query(tree, Base.Fix1(Extents.intersects, polecap)) == [13]

    ni = GO.NaturalIndexing.NaturalIndex(GO.Spherical(), geoms)
    @test Extents.extent(ni) isa Extents.Extent{(:X, :Y, :Z)}
    @test STI.query(ni, pole_box) == [13]

    # a ring of UnitSphericalPoints needs no geographic conversion
    z = 0.9; s = sqrt(1 - z^2)
    usp_ring = GI.LinearRing([UnitSphericalPoint(s * cos(t), s * sin(t), z)
                              for t in range(0, 2π; length = 9)[1:8]])
    @test query(RTree(GO.Spherical(), STR(), [usp_ring]), pole_box) == [1]

    # Planar() matches the manifold-less constructors
    pt = RTree(GO.Planar(), HPR(), band)
    t = RTree(HPR(), band)
    @test pt.levels == t.levels && pt.indices == t.indices
end
