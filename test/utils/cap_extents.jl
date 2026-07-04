using Test

import GeometryOps as GO
import GeometryOps.UnitSpherical as US
import GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap
import GeometryOps.SpatialTreeInterface as STI
import GeometryOps.FlexibleRTrees: RTree, STR, HPR, query
import Extents
using Random: Xoshiro
using LinearAlgebra: norm, normalize, cross, dot

randsphere(rng) = normalize(UnitSphericalPoint(randn(rng), randn(rng), randn(rng)))

# A point at angle θ from `c`, in a random direction tangent at `c`.
function point_at_angle(rng, c, θ)
    t = cross(c, randsphere(rng))
    while norm(t) < 1e-6
        t = cross(c, randsphere(rng))
    end
    return cos(θ) * c + sin(θ) * normalize(t)
end

in_box(p, ext) =
    ext.X[1] <= p[1] <= ext.X[2] &&
    ext.Y[1] <= p[2] <= ext.Y[2] &&
    ext.Z[1] <= p[3] <= ext.Z[2]
# With a safety margin, so float classifications imply real-arithmetic ones.
in_box_margin(p, ext, m) =
    ext.X[1] + m <= p[1] <= ext.X[2] - m &&
    ext.Y[1] + m <= p[2] <= ext.Y[2] - m &&
    ext.Z[1] + m <= p[3] <= ext.Z[2] - m
in_cap_margin(cap, p, m) = dot(cap.point, p) >= cap.radiuslike + m

box_around(p, h) = Extents.Extent(
    X = (p[1] - h, p[1] + h),
    Y = (p[2] - h, p[2] + h),
    Z = (p[3] - h, p[3] + h),
)

# Deterministic near-uniform sphere covering, for disjointness spot checks.
function fibonacci_sphere(n)
    ga = π * (3 - sqrt(5))
    return [begin
        z = 1 - (2i - 1) / n
        r = sqrt(max(1 - z * z, 0.0))
        UnitSphericalPoint(r * cos(ga * i), r * sin(ga * i), z)
    end for i in 1:n]
end

@testset "Extents.extent(cap) contains the cap" begin
    rng = Xoshiro(1)
    centers = (
        UnitSphericalPoint(0.0, 0.0, 1.0), UnitSphericalPoint(0.0, -1.0, 0.0),
        randsphere(rng), randsphere(rng), randsphere(rng),
    )
    radii = (0.0, 1e-8, 0.05, 1.0, π / 2, 2.5, Float64(π))
    for c in centers, r in radii
        cap = SphericalCap(c, r)
        box = Extents.extent(cap)
        @test keys(box) == (:X, :Y, :Z)
        @test in_box(c, box)
        # Interior and rim samples all fall inside.
        @test all(in_box(point_at_angle(rng, c, rand(rng) * r), box) for _ in 1:200)
        @test all(in_box(point_at_angle(rng, c, r), box) for _ in 1:50)
        # Never wider than the unit cube.
        @test all(b -> -1 <= b[1] <= b[2] <= 1, (box.X, box.Y, box.Z))
    end
    # The full sphere's box is the whole unit cube.
    full = Extents.extent(SphericalCap(randsphere(rng), Float64(π)))
    @test all(b -> b == (-1.0, 1.0), (full.X, full.Y, full.Z))
end

@testset "sign kernels match exact rational arithmetic" begin
    rng = Xoshiro(2)
    exact_dot3mk(s, c, k) =
        Int(sign(sum(Rational{BigInt}.(s) .* Rational{BigInt}.(c)) - Rational{BigInt}(k)))
    exact_sqnorm3m1(v) = Int(sign(sum(abs2, Rational{BigInt}.(v)) - 1))
    for _ in 1:2000
        c = Tuple(randsphere(rng))
        s = Tuple(2 * rand(rng) * randsphere(rng))
        @test US._sign_dot3mk(s..., c..., 2 * rand(rng) - 1) isa Int
        # Adversarial: k within a few ulps of the float-evaluated dot product,
        # where only the exact fallback can get the sign right.
        kk = (s[1] * c[1] + s[2] * c[2]) + s[3] * c[3]
        for k in (kk, nextfloat(kk), prevfloat(kk), nextfloat(kk, 3), prevfloat(kk, 7))
            @test US._sign_dot3mk(s..., c..., k) == exact_dot3mk(s, c, k)
        end
        # A normalized vector's ‖v‖² − 1 is within rounding of zero.
        v = Tuple(randsphere(rng))
        @test US._sign_sqnorm3m1(v...) == exact_sqnorm3m1(v)
    end
    # Exactly representable zeros.
    @test US._sign_dot3mk(0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5) == 0
    @test US._sign_sqnorm3m1(1.0, 0.0, 0.0) == 0
    @test US._sign_sqnorm3m1(0.0, -1.0, 0.0) == 0
end

@testset "intersects is never false when a witness exists" begin
    rng = Xoshiro(3)
    for _ in 1:2000
        # Radii and interior angles with enough margin that the constructed
        # float witness certifies a real-arithmetic one.
        r = 1e-3 + rand(rng) * (π - 2e-3)
        cap = SphericalCap(randsphere(rng), r)
        p = point_at_angle(rng, cap.point, 0.9 * r * rand(rng))
        h = 10.0^(-7 + 6 * rand(rng))   # box half-sides from 1e-7 to 1e-1
        box = box_around(p, h)
        @test Extents.intersects(cap, box)
        @test Extents.intersects(box, cap)
    end
    # Degenerate exact witnesses: a point cap in a point box.
    cap = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), 0.0)
    @test cap.radiuslike == 1.0
    @test Extents.intersects(cap, box_around(UnitSphericalPoint(0.0, 0.0, 1.0), 0.0))
end

@testset "agreement with a dense spherical lattice" begin
    rng = Xoshiro(4)
    lattice = fibonacci_sphere(200_000)
    nwitnessed = 0
    ndisjoint = 0
    for _ in 1:300
        cap = SphericalCap(randsphere(rng), rand(rng) * π)
        h = 10.0^(-2 + 2 * rand(rng))   # box half-sides from 0.01 to 1
        box = box_around((0.8 + 0.4 * rand(rng)) * randsphere(rng), h)
        pred = Extents.intersects(cap, box)
        # A margin-certified float witness implies a real point in cap ∩ box,
        # so the predicate must say true.  The converse direction may not
        # hold (the filter is conservative), so it is not asserted.
        if any(p -> in_box_margin(p, box, 1e-9) && in_cap_margin(cap, p, 1e-9), lattice)
            @test pred
            nwitnessed += 1
        elseif !pred
            ndisjoint += 1
        end
    end
    # The sweep exercised both regimes.
    @test nwitnessed > 50
    @test ndisjoint > 20
end

@testset "exact boundary discrimination" begin
    # Point cap at the north pole: the cap is exactly {(0, 0, 1)}.
    cap = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), 0.0)
    # A box whose corner touches the pole exactly: intersecting (closed sets).
    @test Extents.intersects(cap, Extents.Extent(X = (0.0, 1.0), Y = (0.0, 1.0), Z = (1.0, 2.0)))
    # One ulp above the sphere: the box misses the shell.
    above = Extents.Extent(X = (-0.001, 0.001), Y = (-0.001, 0.001), Z = (nextfloat(1.0), 2.0))
    @test !Extents.intersects(cap, above)
    # One ulp below the pole: the box misses the cap's half-space.
    below = Extents.Extent(X = (-0.001, 0.001), Y = (-0.001, 0.001), Z = (0.9, prevfloat(1.0)))
    @test !Extents.intersects(cap, below)
    # 2D extents are a clear error, not a wrong answer.
    @test_throws ArgumentError Extents.intersects(cap, Extents.Extent(X = (0.0, 1.0), Y = (0.0, 1.0)))
end

@testset "cap queries through spatial trees" begin
    rng = Xoshiro(5)
    # Boxes around the edge arcs of a few random great circles, padded by
    # the arc's sagitta so each box really contains its arc.
    E3 = Extents.Extent{(:X, :Y, :Z), NTuple{3, NTuple{2, Float64}}}
    boxes = E3[]
    nsegs = 64
    bulge = 1 - cos(π / nsegs) + 1e-9
    for _ in 1:6
        u = randsphere(rng)
        v = normalize(cross(u, randsphere(rng)))
        pts = [cos(θ) * u + sin(θ) * v for θ in range(0, 2π; length = nsegs + 1)]
        for i in 1:nsegs
            a, b = pts[i], pts[i + 1]
            push!(boxes, Extents.Extent(
                X = (min(a[1], b[1]) - bulge, max(a[1], b[1]) + bulge),
                Y = (min(a[2], b[2]) - bulge, max(a[2], b[2]) + bulge),
                Z = (min(a[3], b[3]) - bulge, max(a[3], b[3]) + bulge),
            ))
        end
    end
    trees = (RTree(STR(), boxes), RTree(HPR(), boxes), GO.NaturalIndexing.NaturalIndex(boxes))
    for _ in 1:25
        cap = SphericalCap(randsphere(rng), 0.02 + 2 * rand(rng))
        truth = sort!(findall(b -> Extents.intersects(cap, b), boxes))
        for tree in trees
            @test STI.query(tree, cap) == truth
        end
        # FlexibleRTrees re-exports the same query.
        @test query(trees[1], cap) == truth
    end
    # Plain 3D box queries agree with brute force through the same path.
    for _ in 1:10
        box = box_around(randsphere(rng), 0.2)
        truth = sort!(findall(b -> Extents.intersects(box, b), boxes))
        for tree in trees
            @test STI.query(tree, box) == truth
        end
    end
end
