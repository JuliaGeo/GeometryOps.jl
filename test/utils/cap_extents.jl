using Test

import GeometryOps as GO
import GeometryOps.UnitSpherical as US
import GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap, slerp, arc_extent
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

# A cap straight from its `radiuslike`, the float the predicates treat as
# authoritative (`radius` is kept consistent but unused by them).
cap_k(c, k) = SphericalCap{Float64}(c, acos(clamp(k, -1.0, 1.0)), k)

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

# S2-style adversarial points: with probability 1/3 per coordinate, crush it
# toward a coordinate plane by a log-uniform factor, where degeneracies are
# representable (cf. s2predicates_test.cc ChoosePoint).
function adversarial_point(rng)
    p = randsphere(rng)
    q = ntuple(i -> rand(rng) < 1 / 3 ? p[i] * 10.0^(-50 * rand(rng)) : p[i], 3)
    return normalize(UnitSphericalPoint(q...))
end

# Independent angle-space ground truth at 512-bit precision:
# (intersects, x ⊇ y).  Inputs must not sit on an exact real-arithmetic tie.
function cap_truth(x, y)
    setprecision(BigFloat, 512) do
        cx, cy = BigFloat.(Tuple(x.point)), BigFloat.(Tuple(y.point))
        kx, ky = BigFloat(Float64(x.radiuslike)), BigFloat(Float64(y.radiuslike))
        nx, ny = sqrt(sum(abs2, cx)), sqrt(sum(abs2, cy))
        ky > ny && return (false, true)     # y empty
        kx > nx && return (false, false)    # x empty (y is not)
        kx <= -nx && return (true, true)    # x is the whole sphere
        rx = acos(clamp(kx / nx, big(-1.0), big(1.0)))
        ry = acos(clamp(ky / ny, big(-1.0), big(1.0)))
        d = acos(clamp(sum(cx .* cy) / (nx * ny), big(-1.0), big(1.0)))
        return (d <= rx + ry, d + ry <= rx)
    end
end

@testset "arc_extent contains its arc" begin
    rng = Xoshiro(6)
    for _ in 1:500
        a = randsphere(rng)
        θ = min(10.0^(-8 + 8.49 * rand(rng)), π - 1e-9)   # ~1e-8 up to nearly π
        b = point_at_angle(rng, a, θ)
        box = arc_extent(a, b)
        @test in_box(a, box) && in_box(b, box)
        @test in_box(slerp(a, b, 0.5), box)               # max-sagitta bulge
        @test all(in_box(slerp(a, b, t), box) for t in 0.05:0.05:0.95)
    end
    # Degenerate inputs: coincident, exactly-proportional (non-unit), and
    # antipodal endpoints — finite boxes containing the inputs, never NaN.
    a = randsphere(rng)
    @test in_box(a, arc_extent(a, a))
    prop = arc_extent(a, UnitSphericalPoint((a * prevfloat(1.0))...))
    @test in_box(a, prop) && all(isfinite, Iterators.flatten((prop.X, prop.Y, prop.Z)))
    anti = arc_extent(a, -a)   # ambiguous minor arc: must cover every route
    @test all(in_box(randsphere(rng), anti) for _ in 1:100)
    # The geographic convenience method goes through UnitSphericalPoint.
    geo = arc_extent((0.0, 0.0), (90.0, 0.0))
    @test in_box(GO.UnitSpherical.UnitSphereFromGeographic()((45.0, 0.0)), geo)
end

@testset "Extents.extent(cap) is two-sided tight per axis" begin
    rng = Xoshiro(12)
    pad = 2 * sqrt(eps(Float64))
    for _ in 1:500
        c = randsphere(rng)
        k = 1.98 * rand(rng) - 0.99   # clearly nonempty, non-full caps
        box = Extents.extent(cap_k(c, k))
        for (i, bounds) in enumerate((box.X, box.Y, box.Z))
            hi_true = c[i] >= k ? 1.0 : c[i] * k + sqrt(max((1 - c[i]^2) * (1 - k^2), 0.0))
            lo_true = -c[i] >= k ? -1.0 : c[i] * k - sqrt(max((1 - c[i]^2) * (1 - k^2), 0.0))
            @test lo_true - 3 * pad <= bounds[1] <= lo_true
            @test hi_true <= bounds[2] <= hi_true + 3 * pad
        end
    end
end

@testset "cap–cap predicates match high-precision ground truth" begin
    rng = Xoshiro(7)
    for i in 1:400
        # Half plain random, half adversarial axis-crushed centers with
        # log-uniform radii, where the float screen earns its keep.
        cx = i <= 200 ? randsphere(rng) : adversarial_point(rng)
        cy = i <= 200 ? randsphere(rng) : adversarial_point(rng)
        kx = i <= 200 ? 2 * rand(rng) - 1 : cos(π * 10.0^(-10 * rand(rng)))
        ky = 2 * rand(rng) - 1
        x, y = cap_k(cx, kx), cap_k(cy, ky)
        ti, tc = cap_truth(x, y)
        @test US._intersects(x, y) == ti
        @test US._intersects(y, x) == ti
        @test US._disjoint(x, y) == !ti
        @test US._contains(x, y) == tc
    end
end

@testset "cap–cap tangency to the ulp" begin
    rng = Xoshiro(10)
    for _ in 1:150
        cx, cy = randsphere(rng), randsphere(rng)
        kx = 2 * rand(rng) - 1
        # Solve the external-tangency radiuslike for y in high precision,
        # round to Float64, and sweep ±2 ulps across the boundary.
        roots = setprecision(BigFloat, 512) do
            bcx, bcy, bkx = BigFloat.(Tuple(cx)), BigFloat.(Tuple(cy)), BigFloat(kx)
            Sx, Sy = sum(abs2, bcx), sum(abs2, bcy)
            D = sum(bcx .* bcy)
            disc = (Sx - bkx^2) * (Sx * Sy - D^2)
            disc < 0 ? BigFloat[] : [(bkx * D + sqrt(disc)) / Sx, (bkx * D - sqrt(disc)) / Sx]
        end
        for ky0 in roots
            (-1 <= ky0 <= 1 && BigFloat(Float64(ky0)) != ky0) || continue
            kyf = Float64(ky0)
            for ky in (prevfloat(kyf, 2), prevfloat(kyf), kyf, nextfloat(kyf), nextfloat(kyf, 2))
                x, y = cap_k(cx, kx), cap_k(cy, ky)
                @test US._intersects(x, y) == cap_truth(x, y)[1]
            end
        end
        # Same for internal tangency and containment.
        kb = 2 * rand(rng) - 1
        roots2 = setprecision(BigFloat, 512) do
            bcb, bcs, bkb = BigFloat.(Tuple(cx)), BigFloat.(Tuple(cy)), BigFloat(kb)
            Sb, Ss = sum(abs2, bcb), sum(abs2, bcs)
            D = sum(bcb .* bcs)
            disc = (Sb * Ss - D^2) * (Sb - bkb^2)
            disc < 0 ? BigFloat[] : [(D * bkb + sqrt(disc)) / Sb, (D * bkb - sqrt(disc)) / Sb]
        end
        for ks0 in roots2
            (-1 <= ks0 <= 1 && BigFloat(Float64(ks0)) != ks0) || continue
            ksf = Float64(ks0)
            for ks in (prevfloat(ksf, 2), prevfloat(ksf), ksf, nextfloat(ksf), nextfloat(ksf, 2))
                b, s = cap_k(cx, kb), cap_k(cy, ks)
                @test US._contains(b, s) == cap_truth(b, s)[2]
            end
        end
    end
end

@testset "cap–point containment is exact" begin
    rng = Xoshiro(8)
    # dot(p, c) == k exactly by construction, then one-ulp discrimination.
    for k in (0.75, -0.25, 0.0, 1.0)
        cap = cap_k(UnitSphericalPoint(0.0, 0.0, 1.0), k)
        s = sqrt(max(1 - k^2, 0.0))
        @test US._contains(cap, UnitSphericalPoint(s, 0.0, k))               # exact tie
        @test US._contains(cap, UnitSphericalPoint(s, 0.0, nextfloat(k)))
        @test !US._contains(cap, UnitSphericalPoint(s, 0.0, prevfloat(k)))
    end
    # Random and adversarial agreement with rational arithmetic.
    for i in 1:1000
        cap = cap_k(i <= 500 ? randsphere(rng) : adversarial_point(rng), 2 * rand(rng) - 1)
        p = i <= 500 ? randsphere(rng) : adversarial_point(rng)
        truth = sum(Rational{BigInt}.(Float64.(Tuple(p))) .* Rational{BigInt}.(Float64.(Tuple(cap.point)))) >=
            Rational{BigInt}(Float64(cap.radiuslike))
        @test US._contains(cap, p) == truth
    end
end

@testset "cap predicate consistency laws" begin
    rng = Xoshiro(9)
    # `radiuslike ≤ -‖center‖` is the whole sphere; -1.5 makes that robust
    # to the center's float norm.  With k = -1 exactly and a center of norm
    # just over 1, the cap misses a ~√eps disk at the antipode — that sharp
    # edge of the raw half-space semantics is asserted, not papered over.
    full = cap_k(randsphere(rng), -1.5)
    notquitefull = cap_k(UnitSphericalPoint(1.0, 1e-8, 0.0), -1.0)   # ‖c‖² = 1 + 1e-16
    antipodal = cap_k(UnitSphericalPoint(-1.0, 0.0, 0.0), 0.9)
    @test US._contains(full, antipodal)
    @test !US._contains(notquitefull, antipodal)
    void = SphericalCap{Float64}(UnitSphericalPoint(0.5, 0.0, 0.0), 0.0, 1.0)   # k > ‖c‖
    for _ in 1:300
        x = cap_k(randsphere(rng), 1.98 * rand(rng) - 0.99)
        y = cap_k(randsphere(rng), 1.98 * rand(rng) - 0.99)
        @test US._intersects(x, y) == US._intersects(y, x)
        @test !US._contains(x, y) || US._intersects(x, y)   # containment ⇒ intersection
        @test US._contains(x, x) && US._intersects(x, x)    # closed ⇒ reflexive
        @test US._contains(full, y) && US._intersects(full, y)
        @test !US._contains(void, y) && !US._intersects(void, y)
        @test US._contains(y, void)
    end
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
