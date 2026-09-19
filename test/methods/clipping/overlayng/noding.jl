# Tests for the OverlayNG phase-1 noding substrate (design §2.9): the
# `NodedArrangement` invariants, along-segment ordering vs the exact authority,
# the certified emission fast paths, and the rounded-arrangement / classification
# censuses on a small Natural Earth subset.

using Test
include(joinpath(@__DIR__, "common.jl"))
import GeometryOps: Planar, Spherical, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic
using LinearAlgebra: cross, dot, norm
import Random
using Random: MersenneTwister

@testset "validated UnitSphericalPoint arithmetic" begin
    CF = GO.ValidatedFloats.CrossingFloats
    a = UnitSphericalPoint(CF.CrossingFloat.((0.8, 0.6, 0.0)))
    b = UnitSphericalPoint(CF.CrossingFloat.((0.8, 0.0, 0.6)))
    c = @inferred cross(a, b)
    q = @inferred dot(a, b)
    @test Tuple(c) === cross(Tuple(a), Tuple(b))
    @test q === dot(Tuple(a), Tuple(b))
    @test which(cross, (typeof(a), typeof(b))).module === CF
    @test which(dot, (typeof(a), typeof(b))).module === CF
end

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_crossing_ids(arr) = [Int32(i) for i in 1:GO.num_nodes(arr) if arr.nodes.keys[i].is_crossing]

# Whether an emitted node coordinate is the lon/lat point `ll`. The arrangement
# emits in the manifold's own chart by default — `(x, y)` on the plane,
# unit-sphere xyz on the sphere — so the comparison is against `ll`'s kernel
# image, not against `ll` itself.
_near_kernel(m, p, ll; atol = 1e-12) =
    (q = GO._to_kernel_point(m, ll); all(i -> isapprox(p[i], q[i]; atol), eachindex(q)))

# every proper crossing appears as one shared node id on exactly one A-segment and
# one B-segment interior list (invariant 1); node ids are unique (invariant 2);
# no NodedEdge is zero-length (invariant 3).
function check_invariants(arr, na_strings)
    @test length(unique(arr.nodes.keys)) == GO.num_nodes(arr)          # 2
    @test all(e -> e.node_lo != e.node_hi, arr.edges)                  # 3
    @test all(_crossing_ids(arr)) do cid                               # 1
        a_hits = b_hits = 0
        for (si, _, nid) in arr.seg_nodes
            nid == cid && (si <= na_strings ? (a_hits += 1) : (b_hits += 1))
        end
        a_hits >= 1 && b_hits >= 1
    end
end

# exact-always along-parameter order of a segment's interior node ids
function exact_order(arr, s0, s1, ids)
    R = Rational{BigInt}
    dxr = R(GI.x(s1)) - R(GI.x(s0)); dyr = R(GI.y(s1)) - R(GI.y(s0))
    param(id) = let p = GO._exact_node_point(arr.nodes.keys[id])
        (p[1] - R(GI.x(s0))) * dxr + (p[2] - R(GI.y(s0))) * dyr
    end
    return sort(ids; by = param)
end

# exact-always along-parameter order on the sphere
function exact_order_sph(arr, s0, s1, ids)
    Ne = GO._cross3(GO._vec3(True(), s0), GO._vec3(True(), s1))
    dir(id) = GO._exact_node_dir(True(), arr.nodes.keys[id])
    return sort(ids; lt = (i, j) -> GO._dot3(GO._cross3(dir(i), dir(j)), Ne) > 0)
end

# ---------------------------------------------------------------------------
# 1. Invariants on constructed cases (§2.1)
# ---------------------------------------------------------------------------

@testset "two crossing quads" begin
    A = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0), (2.0, 2.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        na = count(ss -> ss.is_a, arr.segstrings)
        @test length(_crossing_ids(arr)) == 2
        check_invariants(arr, na)
        if m isa Planar
            cpts = sort([GO.node_point(arr, i) for i in _crossing_ids(arr)])
            @test cpts == [(2.0, 4.0), (4.0, 2.0)]     # bit-exact on the integer grid
        end
    end
end

@testset "degree-6 node (tier-2 merge of distinct keys)" begin
    # two A lines and one B line all through the origin -> THREE distinct crossing
    # keys coincident there: crossing_node(lineA1,B), crossing_node(lineA2,B), and
    # crossing_node(lineA1,lineA2) from A's own self-noding. Tier 2 must merge all
    # three into one node — which is what makes the node degree 6.
    A = GI.MultiLineString([[(-1.0, -1.0), (1.0, 1.0)], [(-1.0, 1.0), (1.0, -1.0)]])
    B = GI.LineString([(-1.0, 0.0), (1.0, 0.0)])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        cids = _crossing_ids(arr)
        @test length(cids) == 1                          # merged into one node
        @test _near_kernel(m, GO.node_point(arr, cids[1]), (0.0, 0.0))
        # the merged node is incident to both A lines and B (three parent strings)
        na = count(ss -> ss.is_a, arr.segstrings)
        check_invariants(arr, na)
    end
end

@testset "crossing exactly on a third string's vertex" begin
    # A horizontal line crosses B-line-1 (vertical) at the origin, which is also
    # B-line-2's endpoint vertex: the crossing key and the vertex key coincide
    # and tier 2 merges them.
    A = GI.LineString([(-2.0, 0.0), (2.0, 0.0)])
    B = GI.MultiLineString([[(0.0, -2.0), (0.0, 2.0)], [(0.0, 0.0), (1.0, 1.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        # the origin is a single node shared by A, B-line-1 and B-line-2's vertex
        origin_ids = [i for i in 1:GO.num_nodes(arr)
                      if _near_kernel(m, GO.node_point(arr, i), (0.0, 0.0))]
        @test length(origin_ids) == 1
    end
end

@testset "collinear shared boundary — zero phantom crossings" begin
    # edge-adjacent squares sharing the vertex-identical edge x = 2
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        #-- zero phantom crossings: a crossing node is interned only alongside
        #-- its two interior records, so an empty `seg_nodes` implies none exists
        @test isempty(arr.seg_nodes)
    end
end

@testset "a-b-a spike input" begin
    # B retraces (0,0)->(1,1)->(0,0); A crosses the retraced segment at one point,
    # reported by both candidate pairs but the same canonical crossing key.
    A = GI.LineString([(-1.0, 0.5), (2.0, 0.5)])
    B = GI.LineString([(0.0, 0.0), (1.0, 1.0), (0.0, 0.0)])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())      # must not throw
        cids = _crossing_ids(arr)
        @test length(cids) == 1                                 # one merged node
        m isa Planar && @test GO.node_point(arr, cids[1]) == (0.5, 0.5)
    end
end

# ---------------------------------------------------------------------------
# 2. Ordering cross-check: float-filtered order == exact-always order (§2.5)
# ---------------------------------------------------------------------------

@testset "dense comb ordering matches exact (planar)" begin
    A = GI.LineString([(0.0, 0.0), (201.0, 0.0)])
    B = GI.MultiLineString([[(Float64(i) + 0.3, -1.0), (Float64(i) + 0.3, 1.0)] for i in 1:200])
    arr = GO.NodedArrangement(Planar(), A, B; exact = True())
    # A is string 1, its single segment carries all 200 interior crossings
    ids = [nid for (si, k, nid) in arr.seg_nodes if (si, k) == (Int32(1), Int32(1))]
    @test length(ids) == 200
    s0 = arr.segstrings[1].pts[1]; s1 = arr.segstrings[1].pts[2]
    @test ids == exact_order(arr, s0, s1, ids)          # elementwise
    #-- strictly increasing along the segment. The line above already checks the
    #-- order against an independent `Rational{BigInt}` oracle; this adds only
    #-- comparator antisymmetry, so it does not need 199 assertions to say it.
    @test all(c -> GO.rk_compare_along_segment(Planar(), s0, s1,
                  arr.nodes.keys[ids[c-1]], arr.nodes.keys[ids[c]]; exact = True()) < 0,
              2:length(ids))
end

@testset "dense comb ordering matches exact (spherical)" begin
    A = GI.LineString([(0.0, 0.0), (60.0, 0.0)])
    B = GI.MultiLineString([[(Float64(i) * 0.25 + 0.1, -1.0), (Float64(i) * 0.25 + 0.1, 1.0)] for i in 1:200])
    arr = GO.NodedArrangement(Spherical(), A, B; exact = True())
    ids = [nid for (si, k, nid) in arr.seg_nodes if (si, k) == (Int32(1), Int32(1))]
    @test length(ids) == 200
    s0 = arr.segstrings[1].pts[1]; s1 = arr.segstrings[1].pts[2]
    @test ids == exact_order_sph(arr, s0, s1, ids)
end

# ---------------------------------------------------------------------------
# 3. Emission certificate audit (§2.6)
# ---------------------------------------------------------------------------

@testset "planar emission: certified == rational, every node" begin
    # dense generic-slope grid + shifted-self coastline-like crossings
    Ag = GI.MultiLineString([[(Float64(k) * 4.0, 0.0), (Float64(k) * 4.0 + 0.31, 1000.0)] for k in 1:60])
    Bg = GI.MultiLineString([[(0.0, Float64(j) * 4.0), (1000.0, Float64(j) * 4.0 + 0.29)] for j in 1:60])
    arr = GO.NodedArrangement(Planar(), Ag, Bg; exact = True())
    #-- accumulated, the way the spherical sibling below already does it: the
    #-- per-node form asserted 7 202 times to establish three properties
    ncert = 0; ntot = 0; cert_ok = true; emit_ok = true
    for i in _crossing_ids(arr)
        k = arr.nodes.keys[i]
        (x, y, cert) = GO._certified_crossing(k.pt, k.a1, k.b0, k.b1)
        rx, ry = GO._exact_crossing_point(k)
        rat = (Float64(rx), Float64(ry))
        ntot += 1
        cert && (ncert += 1; ((x, y) == rat) || (cert_ok = false))
        (GO.node_point(arr, i) == rat) || (emit_ok = false)
    end
    @test ntot > 1000
    @test ncert == ntot                               # 100% certified on clean data (S3)
    @test cert_ok                                     # certified => equal to the rational answer
    @test emit_ok                                     # node_point is the rational answer either way
end

function _rational_planar_crossing(a0, a1, b0, b1)
    R = Rational{BigInt}
    ax0, ay0 = R(a0[1]), R(a0[2]); ax1, ay1 = R(a1[1]), R(a1[2])
    bx0, by0 = R(b0[1]), R(b0[2]); bx1, by1 = R(b1[1]), R(b1[2])
    dax, day = ax1 - ax0, ay1 - ay0
    dbx, dby = bx1 - bx0, by1 - by0
    c0x, c0y = bx0 - ax0, by0 - ay0
    t = (c0x * dby - c0y * dbx) / (dax * dby - day * dbx)
    return Float64(ax0 + t * dax), Float64(ay0 + t * day)
end

@testset "planar emission: validated arithmetic domains" begin
    large = 0x1p200
    large_step = eps(large)
    outside = 0x1p500
    outside_step = eps(outside)
    cases = [
        ((0.0, 0.0), (10.0, 10.0), (0.0, 10.0), (10.0, 0.0)),
        ((large, large), (large + 16large_step, large + 12large_step),
         (large, large + 12large_step), (large + 16large_step, large)),
        ((0.0, 0.0), (0x1p-200, 0x1p-200),
         (0.0, 0x1p-200), (0x1p-200, 0.0)),
        ((0.0, 0.0), (1.0, 1.0),
         (0.0, -0x1p-41), (1.0, 1.0 + 0x1p-41)),
        ((outside, outside), (outside + 16outside_step, outside + 12outside_step),
         (outside, outside + 12outside_step), (outside + 16outside_step, outside)),
    ]
    ncert = 0
    for (a0, a1, b0, b1) in cases
        x, y, cert = GO._certified_crossing(a0, a1, b0, b1)
        want = _rational_planar_crossing(a0, a1, b0, b1)
        k = GO.crossing_node(a0, a1, b0, b1)
        @test GO._emit_node_coord(k, Tuple{Float64, Float64}) == want
        if cert
            @test (x, y) == want
            ncert += 1
        end
    end
    @test 0 < ncert < length(cases)
    @test !GO._certified_crossing(cases[end]...)[3]
end

#=
Both spherical output rows, against an oracle independent of the emitter. The
contract is correct rounding: the xyz row emits `(RN(x₁), RN(x₂), RN(x₃))` of
the exact unit crossing direction `x = d/‖d‖`, and the lon/lat row is that
point sent through the vertex row's trigonometry.

The oracle normalizes the exact rational `d` in 4096-bit BigFloat and rounds
each component to Float64 once. That decides the same rounding as exact
arithmetic unless `xᵢ` lies within ~2⁻⁴⁰⁹⁰ of a Float64 midpoint without being
one. `xᵢ` is algebraic of degree ≤ 2 (a rational over a square root) with
coefficient height below 2⁷⁰⁰ for Float64 inputs, and a midpoint is a rational
with denominator ≤ 2¹⁰⁷⁵, so Liouville's bound keeps any such distance above
2⁻³⁰⁰⁰. An exact tie is a rational `xᵢ`, which the oracle would still round
correctly only by luck; none of the populations below produces one (it needs
`‖d‖` rational), and the tie rule is tested on constructed inputs in the
`_round_div_sqrt` testset instead.

Populations are the ones the emitter is sensitive to in different ways:
long arcs (nothing cancels), metre/centimetre/micro arcs (`a0×a1` cancels by
the arc length), near-tangent arcs (`na×nb` cancels by the plane angle, into
the exact fallback), and lon/lat grid data, where an equator or meridian arc
puts the crossing on a coordinate plane and one component is exactly zero.
=#
const _SPH = Spherical()
_canon(v) = GO._spherical_kernel_point(UnitSphericalPoint(v[1], v[2], v[3]))
_logu(rng, lo, hi) = exp(log(lo) + rand(rng) * (log(hi) - log(lo)))
function _frame(p)
    e = zeros(3); e[argmin(abs.(p))] = 1.0
    u = cross(p, e) ./ norm(cross(p, e))
    return (u, cross(p, u))
end
_walk(p, t, s) = cos(s) .* p .+ sin(s) .* t

# a proper crossing through a random point, arc lengths La, Lb, plane angle φ
function _crossing_case(rng, La, Lb, φ)
    X = randn(rng, 3); X ./= norm(X); u, w = _frame(X)
    α = 2π * rand(rng)
    ta = cos(α) .* u .+ sin(α) .* w
    tb = cos(α + φ) .* u .+ sin(α + φ) .* w
    fa = 0.1 + 0.8rand(rng); fb = 0.1 + 0.8rand(rng)
    return (_canon(_walk(X, -ta, fa * La)), _canon(_walk(X, ta, (1 - fa) * La)),
            _canon(_walk(X, -tb, fb * Lb)), _canon(_walk(X, tb, (1 - fb) * Lb)))
end

# grid-ish lon/lat data: integer / half / quarter degrees, equator, meridians
function _lonlat_case(rng)
    g() = rand(rng, (1.0, 0.5, 0.25, 0.1))
    λ0 = round(rand(rng) * 360 - 180); φ0 = round(rand(rng) * 170 - 85)
    kind = rand(rng, 1:3)
    a0, a1, b0, b1 = if kind == 1          # meridian arc vs parallel-endpoint arc
        (λ0, φ0 - g()), (λ0, φ0 + g()), (λ0 - g(), φ0), (λ0 + g(), φ0)
    elseif kind == 2                       # equator arc vs slanted arc
        (λ0 - g(), 0.0), (λ0 + g(), 0.0), (λ0 - g() / 3, -g()), (λ0 + g() / 7, g())
    else                                   # prime meridian vs diagonal
        (0.0, φ0 - g()), (0.0, φ0 + g()), (-g(), φ0 - g() / 2), (g(), φ0 + g() / 3)
    end
    f = GO._spherical_kernel_point
    return (f(a0), f(a1), f(b0), f(b1))
end

_is_proper(c) = GO.rk_classify_intersection(_SPH, c...; exact = True()).kind == GO.SS_PROPER

function _sample_crossings(gen, n)
    ks = GO.NodeKey{USP}[]
    tries = 0
    while length(ks) < n && tries < 50n
        tries += 1
        c = gen()
        _is_proper(c) && push!(ks, GO.crossing_node(c...))
    end
    return ks
end

# the independent oracle: exact rational direction, normalized in BigFloat
function _oracle_usp(k)
    d = GO._sph_crossing_dir(True(), k)
    return setprecision(BigFloat, 4096) do
        x = BigFloat.(collect(d)); x ./= sqrt(sum(x .^ 2))
        UnitSphericalPoint(Float64(x[1]) + 0.0, Float64(x[2]) + 0.0, Float64(x[3]) + 0.0)
    end
end

@testset "spherical emission: certified == correctly rounded, per population" begin
    rng = MersenneTwister(11)
    N = 30
    pops = [
        ("long arcs",          () -> _crossing_case(rng, _logu(rng, 0.01, 1.0), _logu(rng, 0.01, 1.0), (0.05 + 0.9rand(rng)) * π), true),
        ("metre arcs",         () -> _crossing_case(rng, 1.6e-7, 1.6e-7 * (0.5 + rand(rng)), (0.05 + 0.9rand(rng)) * π), true),
        ("centimetre arcs",    () -> _crossing_case(rng, 1.6e-9, 1.6e-9, (0.05 + 0.9rand(rng)) * π), true),
        ("micro arcs",         () -> _crossing_case(rng, _logu(rng, 1e-13, 1e-11), _logu(rng, 1e-13, 1e-11), (0.05 + 0.9rand(rng)) * π), false),
        ("near-tangent 100 km", () -> _crossing_case(rng, 0.0157, 0.0157, _logu(rng, 1e-15, 1e-3)), false),
        ("near-tangent 1 km",  () -> _crossing_case(rng, 1.6e-4, 1.6e-4, _logu(rng, 1e-12, 1e-3)), false),
        ("lon/lat grid",       () -> _lonlat_case(rng), true),
    ]
    nfallback = 0
    for (name, gen, all_fast) in pops
        ks = _sample_crossings(gen, N)
        @testset "$name" begin
            @test length(ks) == N
            nfast = 0
            for k in ks
                want = _oracle_usp(k)
                got = GO._emit_node_coord(k, USP)
                @test got === want                                # bit-for-bit, both signs of zero
                #-- the two paths agree with each other, so the output is a
                #-- function of the node's exact position alone
                (ok, x, y, z) = GO._certified_sph_crossing_fast(k)
                ok && (nfast += 1)
                @test (x + 0.0, y + 0.0, z + 0.0) == GO._certified_sph_crossing_exact(k) || !ok
                #-- the lon/lat row is the vertex row's trigonometry applied to it
                @test GO._emit_node_coord(k, Tuple{Float64, Float64}) == GO._usp_to_lonlat(want)
                #-- re-ingest identity: a cascade level's output re-enters the
                #-- next level bit-for-bit
                @test GO.rk_normalize_usp(got) === got
                @test GO._spherical_kernel_point(got) === got
            end
            all_fast && @test nfast == N
            nfallback += N - nfast
        end
    end
    @test nfallback > 0                     # the near-tangent populations reach the exact fallback
end

# An exactly-zero component certifies on the fast path. The neighbour gap at
# 0.0 is 5e-324, so a bound lumped over the expression would send every
# crossing on a coordinate plane to the exact fallback; the per-operation bounds
# propagate an exact zero with a zero bound instead.
@testset "spherical emission: exact zeros on the fast path" begin
    f = GO._spherical_kernel_point
    #-- equator × meridian: z = 0 and the crossing is on the meridian plane
    k = GO.crossing_node(f((0.0, 0.0)), f((10.0, 0.0)), f((5.0, -1.0)), f((5.0, 1.0)))
    (ok, x, y, z) = GO._certified_sph_crossing_fast(k)
    @test ok
    @test z == 0.0
    @test GO._emit_node_coord(k, USP) === _oracle_usp(k)
    @test GO._certified_sph_crossing_fast(k; safety = 4.0) === (false, 0.0, 0.0, 0.0)
    #-- prime meridian × slanted arc: y = 0
    k = GO.crossing_node(f((0.0, 40.0)), f((0.0, 50.0)), f((-1.0, 44.0)), f((1.0, 46.0)))
    (ok, x, y, z) = GO._certified_sph_crossing_fast(k)
    @test ok
    @test y == 0.0
    @test GO._emit_node_coord(k, USP) === _oracle_usp(k)
end

# The exact fallback's rounding: `_round_div_sqrt(a, S)` is `a/√S` rounded to
# nearest, ties to even, decided by exact midpoint comparison.
@testset "_round_div_sqrt: correct rounding incl. exact ties" begin
    R = Rational{BigInt}
    #-- exactly representable quotients come back exactly
    @test GO._round_div_sqrt(R(3), R(4)) == 1.5
    @test GO._round_div_sqrt(R(-3), R(4)) == -1.5
    @test GO._round_div_sqrt(R(0), R(7)) == 0.0
    @test GO._round_div_sqrt(R(1, 3), R(1, 9)) == 1.0
    #-- exact ties go to the even neighbour, from either side of it
    tie_hi = R(1) + R(1, 2)^53                                   # midpoint of 1.0 and nextfloat(1.0)
    @test GO._round_div_sqrt(tie_hi, R(1)) == 1.0
    @test GO._round_div_sqrt(tie_hi * 3, R(9)) == 1.0
    odd = nextfloat(1.0)                                         # odd significand
    tie_odd = (R(odd) + R(nextfloat(odd))) / 2
    @test GO._round_div_sqrt(tie_odd, R(1)) == nextfloat(odd)    # even neighbour is above
    tie_lo = (R(1.0) + R(prevfloat(1.0))) / 2                    # power-of-two boundary
    @test GO._round_div_sqrt(tie_lo, R(1)) == 1.0
    #-- generic quotients agree with a high-precision evaluation
    rng = MersenneTwister(5)
    for _ in 1:200
        a = R(randn(rng) * 2.0^rand(rng, -40:40)); S = R(abs(randn(rng)) * 2.0^rand(rng, -80:80)) + R(1, 10^9)
        want = setprecision(BigFloat, 4096) do
            Float64(BigFloat(a) / sqrt(BigFloat(S)))
        end
        @test GO._round_div_sqrt(a, S) == want
    end
end

# The same contract through the arrangement: crossing keys as noding produces
# them, realized and cached by `node_point`, on both output rows selected via
# `point_type`. One assertion per property, accumulated over the ~3600
# crossings; the 4096-bit oracle costs ~1 ms per crossing, so this is the
# file's slowest block by design and needs no sampling.
@testset "spherical emission through the arrangement: every node_point == oracle" begin
    Ag = GI.MultiLineString([[(Float64(k) * 0.09 + 0.05, 0.0), (Float64(k) * 0.09 + 0.05 + 0.031, 20.0)] for k in 1:60])
    Bg = GI.MultiLineString([[(0.0, Float64(j) * 0.09 + 0.05), (20.0, Float64(j) * 0.09 + 0.05 + 0.029)] for j in 1:60])

    #-- the default row: unit-sphere xyz, bit-for-bit the correctly rounded direction
    arr_x = GO.NodedArrangement(Spherical(), Ag, Bg; exact = True())
    ids = _crossing_ids(arr_x)
    oracle = Dict(i => _oracle_usp(arr_x.nodes.keys[i]) for i in ids)
    all_usp = true; all_equal = true; all_cached = true
    for i in ids
        emitted = GO.node_point(arr_x, i)
        emitted isa UnitSphericalPoint{Float64} || (all_usp = false)
        emitted === oracle[i] || (all_equal = false)
        GO.node_point(arr_x, i) === emitted || (all_cached = false)   # realized once, read back
    end
    @test length(ids) > 1000
    @test all_usp
    @test all_equal
    @test all_cached

    #-- the lon/lat row: the vertex trigonometry applied to that same point
    arr = GO.NodedArrangement(Spherical(), Ag, Bg; exact = True(),
                              point_type = Tuple{Float64, Float64})
    all_ll = true
    for i in _crossing_ids(arr)
        GO.node_point(arr, i) == GO._usp_to_lonlat(_oracle_usp(arr.nodes.keys[i])) || (all_ll = false)
    end
    @test length(_crossing_ids(arr)) == length(ids)
    @test all_ll
end

# ---------------------------------------------------------------------------
# 4 & 5. Natural Earth subset: rounded-arrangement audit + classification census
# ---------------------------------------------------------------------------

# Rounded-arrangement audit — the one assertion in this file about the property
# the whole substrate exists to guarantee: rounding at emission does not
# introduce topology. No edge incident to a crossing node may properly cross an
# OPPOSITE-SIDE edge once both are realized as Float64.
#
# The scope has to reach PAST the shared node. This function used to compare only
# edges incident to the *same* crossing node, and every such pair carries that
# node's single emitted coordinate as an endpoint — two segments sharing an
# endpoint are `SS_TOUCH` by construction and can never be `SS_PROPER`. Measured
# on the five fixtures below, the old audit examined 184 / 88 / 36 / 136 / 200
# pairs and every one of them shared an endpoint, so it returned 0 by
# construction and would have returned 0 for any arrangement whatsoever.
#
# The corrected audit is discriminating: with emission snapped to 0.1° it finds
# 7 / 0 / 2 / 4 / 8 introduced crossings on the same fixtures, where the old one
# still reports 0.
#
# `m` is the manifold the emitted coordinates are interpreted in, and it is the
# arrangement's own — a spherical arrangement's emitted points are read by
# downstream stages as spherical, so a lon/lat-planar audit of them would be
# asking a different question than the one the guarantee is about. Both default
# arrangements now emit in their manifold's kernel chart, so this lift is the
# identity on both; it stays spelled out because the audit is also the natural
# place to run a lon/lat arrangement through, where it is not.
emitted_pt(::Planar, p) = p
emitted_pt(::Spherical, p) = UnitSphericalPoint(GI.PointTrait(), p)

function rounded_crossings(m, arr)
    ends = [(emitted_pt(m, GO.node_point(arr, e.node_lo)),
             emitted_pt(m, GO.node_point(arr, e.node_hi))) for e in arr.edges]
    is_a = [arr.segstrings[e.string_idx].is_a for e in arr.edges]
    incident = [n for (n, e) in enumerate(arr.edges)
                if arr.nodes.keys[e.node_lo].is_crossing || arr.nodes.keys[e.node_hi].is_crossing]
    introduced = 0
    for n in incident
        (pa0, pa1) = ends[n]
        for mi in eachindex(arr.edges)
            mi == n && continue
            is_a[mi] == is_a[n] && continue          # opposite side only
            (pb0, pb1) = ends[mi]
            GO.rk_classify_intersection(m, pa0, pa1, pb0, pb1; exact = True()).kind ==
                GO.SS_PROPER && (introduced += 1)
        end
    end
    return introduced
end

ne_ok = try
    import NaturalEarth, GeoJSON
    include(joinpath(@__DIR__, "..", "..", "..", "data", "natural_earth_pairs.jl"))
    global ne_names, ne_geoms = load_ne(110)
    length(ne_geoms) > 0
catch err
    @info "Natural Earth subset skipped (data unavailable)" err
    false
end

# The flag check this block used to run — "every touch/collinear classification
# carries a vertex-incidence flag", via a `classify_census` helper that
# re-enumerated the A×B candidate pairs with the same index, the same traversal
# and the same kernel call the noder uses — is a live `@assert` inside
# `collect.jl` three lines from the classification itself. A test-side
# reimplementation of noding stage 1 is not a second opinion; it is the same
# opinion, computed twice.
@testset "Natural Earth subset (rounded-arrangement audit)" begin
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
            arr_p = GO.NodedArrangement(Planar(), A, B; exact = True())
            arr_s = GO.NodedArrangement(Spherical(), A, B; exact = True())
            #-- the fixture must actually cross, or the audit below is vacuous
            @test !isempty(_crossing_ids(arr_p)) && !isempty(_crossing_ids(arr_s))
            @test rounded_crossings(Planar(), arr_p) == 0
            @test rounded_crossings(Spherical(), arr_s) == 0
        end
        @test tested >= 2
    end
end

# ---------------------------------------------------------------------------
# 6. Coincidence sweep positions ill-conditioned crossings exactly
# ---------------------------------------------------------------------------

@testset "coincidence sweep merges an ill-conditioned exact coincidence" begin
    #-- the measured counterexample: a proper crossing whose exact point is (1,1)
    #-- but whose float solve lands ~1.2e-7 off with a determinant-only error
    #-- claim of ~1e-14; the vertex node at (1,1) coincides with it exactly
    a, b = (4.096955890625e8, 1.25542458125e8), (-8.19391175125e8, -2.5108491325e8)
    c, d = (-6.4578703125e6, 2.44491420625e8), (1.93736149375e7, -7.33474257875e8)
    k = GO.crossing_node(a, b, c, d)
    vk = GO.vertex_node((1.0, 1.0))
    @test GO._exact_node_point(k) == (1 // 1, 1 // 1)
    x, y, err = GO._approx_node_point(k)
    @test err >= abs(x - 1.0) + abs(y - 1.0)          # the radius covers the miss
    P = Tuple{Float64, Float64}
    t = GO.NodeTable{P, P}()
    GO._intern_node!(t, vk); GO._intern_node!(t, k)
    parent = Int32[1, 2]
    @test GO._coincidence_sweep!(Planar(), t, parent; exact = True()) == 1
    @test GO._uf_find(parent, Int32(1)) == GO._uf_find(parent, Int32(2))
end

@testset "near-coincident crossing pairs node without throwing" begin
    #-- seeded batch of the same family: a crossing exactly on a third string's
    #-- vertex, with endpoints large enough that the float solve is off by
    #-- far more than the 1e-8 proximity gate
    rng = Random.Xoshiro(54)
    v = (1.0, 1.0)
    dyadic() = rand(rng, -10_000_000_000:10_000_000_000) / 16
    nbatch = 0
    while nbatch < 300
        u = (dyadic(), dyadic()); w = (dyadic(), dyadic())
        u[1] * w[2] == u[2] * w[1] && continue
        a, b, c, d = v .- u, v .+ 2 .* u, v .- w, v .+ 3 .* w
        A = GI.LineString([a, b])
        B = GI.MultiLineString([GI.LineString([c, d]), GI.LineString([(0.0, -5.0), v, (7.0, 0.0)])])
        arr = GO.NodedArrangement(Planar(), A, B; exact = True())
        #-- the crossing and the vertex are one node
        @test count(i -> GO.node_point(arr, i) == v, 1:GO.num_nodes(arr)) == 1
        nbatch += 1
    end
end
