# Tests for the public opt-in `OverlayNG{M}` algorithm (design §4): the four ops
# on both manifolds through the exported surface, the new `symdifference`, the
# algorithm-type plumbing, the error paths, and — the load-bearing one — that
# the DEFAULT engine for `intersection`/`union`/`difference` is still
# Foster–Hormann and behaves exactly as before.
#
# Equality: planar results are compared against LibGEOS's own overlay with GEOS
# topological `equals` plus `isValid`; spherical results are gated on the area
# conservation identities, which are machine-precision gates.

using Test
import ArchGDAL as AG
import GeometryBasics as GB
using GeometryOpsTestHelpers
include(joinpath(@__DIR__, "common.jl"))

const P = GO.Planar()
const S = GO.Spherical()

A, B = SQ_A, SQ_B
#-- a polygon with a hole, and a multipolygon, for the harder shapes
AH = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                 [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]])
BH = GI.Polygon([[(5.0, 5.0), (12.0, 5.0), (12.0, 12.0), (5.0, 12.0), (5.0, 5.0)]])

# ---------------------------------------------------------------------------
# 1. Every op through the public surface, planar, vs LibGEOS
# ---------------------------------------------------------------------------

const PUBLIC_OPS = (
    ("intersection", GO.intersection, LG.intersection),
    ("union",        GO.union,        LG.union),
    ("difference",   GO.difference,   LG.difference),
    ("symdifference", GO.symdifference, LG.symmetricDifference),
)

@testset "planar ops through OverlayNG()" begin
    for (name, gof, lgf) in PUBLIC_OPS, (a, b) in ((A, B), (AH, BH))
        @testset "$name" begin
            r = gof(GO.OverlayNG(), a, b)
            @test LG.isValid(lgc(r))
            @test LG.equals(lgc(r), lgf(lgc(a), lgc(b)))
        end
    end
end

@testset "mixed-dimension and point inputs through the public surface" begin
    L = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    #-- a line clipped by an area
    @test LG.equals(lgc(GO.intersection(GO.OverlayNG(), A, L)),
                    LG.intersection(lgc(A), lgc(L)))
    #-- point inputs
    pts = GI.MultiPoint([(0.5, 0.5), (5.0, 5.0)])
    @test LG.equals(lgc(GO.intersection(GO.OverlayNG(), pts, A)), LG.readgeom("POINT (0.5 0.5)"))
    @test LG.equals(lgc(GO.difference(GO.OverlayNG(), pts, A)), LG.readgeom("POINT (5 5)"))
    #-- union of a point outside an area is a mixed collection
    r = GO.union(GO.OverlayNG(), A, GI.Point((5.0, 5.0)))
    @test GI.trait(r) isa GI.GeometryCollectionTrait
    #-- point × point
    @test LG.equals(lgc(GO.union(GO.OverlayNG(), GI.Point((1.0, 1.0)), GI.Point((2.0, 2.0)))),
                    LG.readgeom("MULTIPOINT ((1 1), (2 2))"))
end

# ---------------------------------------------------------------------------
# 2. Every op through the public surface, spherical
# ---------------------------------------------------------------------------

@testset "spherical ops through OverlayNG(Spherical())" begin
    alg = GO.OverlayNG(S)
    SA = GI.Polygon([[(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0), (0.0, 0.0)]])
    SB = GI.Polygon([[(10.0, 10.0), (30.0, 10.0), (30.0, 30.0), (10.0, 30.0), (10.0, 10.0)]])
    aA = GO.area(S, SA); aB = GO.area(S, SB)
    ai = GO.area(S, GO.intersection(alg, SA, SB))
    au = GO.area(S, GO.union(alg, SA, SB))
    ad = GO.area(S, GO.difference(alg, SA, SB))
    as = GO.area(S, GO.symdifference(alg, SA, SB))
    @test isapprox(au + ai, aA + aB; rtol = 1e-12)
    @test isapprox(ad + ai, aA; rtol = 1e-12)
    @test isapprox(as, au - ai; rtol = 1e-12)
    #-- all four are strictly positive here, i.e. none collapsed to empty
    @test all(>(0), (ai, au, ad, as))
    #-- and the spherical answers differ from the planar ones (great-circle edges)
    @test !isapprox(ai, GO.area(S, GO.intersection(GO.OverlayNG(), SA, SB)); rtol = 1e-9)
end

@testset "spherical point and mixed inputs through the public surface" begin
    alg = GO.OverlayNG(S)
    #-- the ring's northern edge bows to ~67.8°N, so (45, 62) is inside
    PG = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (90.0, 60.0), (0.0, 60.0), (0.0, 0.0)]])
    pt = GI.Point((45.0, 62.0))
    @test GI.trait(GO.intersection(alg, PG, pt)) isa GI.PointTrait
    #-- planar disagrees, which is the point of carrying the manifold: the planar
    #-- result is the EMPTY Multi form, not the atomic point
    plan = GO.intersection(GO.OverlayNG(), PG, pt)
    @test GI.trait(plan) isa GI.MultiPointTrait && GI.npoint(plan) == 0
end

# ---------------------------------------------------------------------------
# 3. `symdifference` — the new public function and its engine decision
# ---------------------------------------------------------------------------

@testset "symdifference entry points" begin
    expected = LG.symmetricDifference(lgc(A), lgc(B))
    #-- manifold-only and bare forms agree (the explicit-algorithm form is the
    #-- symdifference row of the public-ops table above)
    @test LG.equals(lgc(GO.symdifference(P, A, B)), expected)
    @test LG.equals(lgc(GO.symdifference(A, B)), expected)
    #-- the manifold form routes to OverlayNG with that manifold. Compared by
    #-- COORDINATES, not by summed area: the crossing node here is the shared
    #-- corner BETWEEN the symmetric difference's two components, so moving it
    #-- transfers area from one to the other and leaves the sum invariant — an
    #-- area identity passes even when the manifold is dropped entirely.
    @test GI.coordinates(GO.symdifference(S, A, B)) ==
          GI.coordinates(GO.symdifference(GO.OverlayNG(S), A, B))
    @test GI.coordinates(GO.symdifference(S, A, B)) !=
          GI.coordinates(GO.symdifference(P, A, B))
    #-- symdifference(A, B) == union(diff(A,B), diff(B,A)) by area
    @test isapprox(GO.area(GO.symdifference(A, B)),
                   GO.area(GO.difference(GO.OverlayNG(), A, B)) +
                   GO.area(GO.difference(GO.OverlayNG(), B, A)); rtol = 1e-12)
    #-- symmetric in its arguments
    @test LG.equals(lgc(GO.symdifference(A, B)), lgc(GO.symdifference(B, A)))
    #-- disjoint inputs give both operands back
    Far = GI.Polygon([[(50.0, 50.0), (52.0, 50.0), (52.0, 52.0), (50.0, 52.0), (50.0, 50.0)]])
    @test isapprox(GO.area(GO.symdifference(A, Far)), GO.area(A) + GO.area(Far); rtol = 1e-12)
    #-- identical inputs give the empty geometry
    e = GO.symdifference(A, A)
    @test GI.npoint(e) == 0
end

# ---------------------------------------------------------------------------
# 4. Implementation-generic coverage
# ---------------------------------------------------------------------------

@testset_implementations "OverlayNG ops across geometry implementations" begin
    @test isapprox(GO.area(GO.intersection(GO.OverlayNG(), $A, $B)), 1.0; rtol = 1e-12)
    @test isapprox(GO.area(GO.union(GO.OverlayNG(), $A, $B)), 7.0; rtol = 1e-12)
    @test isapprox(GO.area(GO.difference(GO.OverlayNG(), $A, $B)), 3.0; rtol = 1e-12)
    @test isapprox(GO.area(GO.symdifference($A, $B)), 6.0; rtol = 1e-12)
    #-- spherical conservation holds for every implementation too
    @test isapprox(GO.area(S, GO.union(GO.OverlayNG(S), $A, $B)) +
                   GO.area(S, GO.intersection(GO.OverlayNG(S), $A, $B)),
                   GO.area(S, $A) + GO.area(S, $B); rtol = 1e-12)
end

@testset_implementations "OverlayNG with a hole across geometry implementations" begin
    #-- A has a hole that B reaches into (the §2.7 material-interior case)
    @test isapprox(GO.area(GO.intersection(GO.OverlayNG(), $AH, $BH)), 21.0; rtol = 1e-12)
    @test isapprox(GO.area(GO.difference(GO.OverlayNG(), $AH, $BH)),
                   GO.area($AH) - 21.0; rtol = 1e-12)
end

# ---------------------------------------------------------------------------
# 5. Algorithm-type plumbing
# ---------------------------------------------------------------------------

@testset "OverlayNG algorithm type" begin
    @test GO.OverlayNG() isa GO.GeometryOpsCore.Algorithm{<:GO.Planar}
    @test GO.OverlayNG(S) isa GO.GeometryOpsCore.Algorithm{<:GO.Spherical}
    @test GO.GeometryOpsCore.manifold(GO.OverlayNG()) === P
    @test GO.GeometryOpsCore.manifold(GO.OverlayNG(S)) === S
    #-- rebuild swaps the manifold and keeps the exactness setting
    reb = GO.GeometryOpsCore.rebuild(GO.OverlayNG(; exact = GO.False()), S)
    @test GO.GeometryOpsCore.manifold(reb) === S
    @test reb.exact === GO.False()
    #-- the positional-manifold constructor is the keyword one, so comparing the
    #-- two says nothing; what has content is that the DEFAULT manifold is planar
    @test GO.OverlayNG(P) === GO.OverlayNG()
    #-- `exact = False()` still runs (Float64 filters only) and agrees here...
    let alg = GO.OverlayNG(; exact = GO.False())
        @test LG.equals(lgc(GO.intersection(alg, A, B)), LG.intersection(lgc(A), lgc(B)))
        #-- ...and the field really reaches the driver. Without this the wrappers
        #-- could hard-code `exact = True()` and every assertion in this file
        #-- would still pass — only `rebuild` above was pinning `exact` at all.
        @test alg.exact === GO.False()
        @test GI.coordinates(GO.intersection(alg, A, B)) ==
              GI.coordinates(GO._overlay_ng(P, GO.OVERLAY_INTERSECTION, A, B; exact = alg.exact))
    end
end

#=
## `point_type` — the chart the one rounding lands in

The arrangement is exact and symbolic; `node_point` is the only place a Float64
coordinate is ever constructed. `point_type` says which chart that construction
targets, and on the sphere the choice is not neutral, because the engine's own
chart is unit-sphere xyz: emitting xyz hands a vertex node BACK the coordinate
it was ingested as, while emitting (lon, lat) sends it out through `atan`/`asin`
and gets something else.

The tests below are ordered by what they are protecting. The bit-exactness one
is the motivating defect and the reason the default changed — a clipped grid
cell whose vertices the clip did not touch used to come back an ULP off, and
sums over such cells then failed exact-equality identities downstream.
=#
@testset "point_type: the spherical default is the engine's own chart" begin
    #-- default output types, per manifold
    @test GO.OverlayNG().point_type === Tuple{Float64, Float64}
    @test GO.OverlayNG(S).point_type === USP
    #-- and they are exactly `_kernel_point_type`, which is the whole rule
    @test GO.OverlayNG(S).point_type === GO._kernel_point_type(S)
    @test GO.OverlayNG(P).point_type === GO._kernel_point_type(P)

    r = GO.intersection(GO.OverlayNG(S), A, B)
    @test GI.is3d(r)
    @test all(p -> p isa USP, GI.getpoint(r))
    #-- planar is untouched: same default, same 2D result
    rp = GO.intersection(GO.OverlayNG(), A, B)
    @test !GI.is3d(rp)
    @test all(p -> p isa Tuple{Float64, Float64}, GI.getpoint(rp))
end

@testset "point_type: an uncut vertex survives a spherical overlay bit-for-bit" begin
    #-- a cell strictly inside a bigger one: the intersection IS the small cell,
    #-- and the overlay cuts none of its vertices. Every emitted coordinate must
    #-- therefore be the ingested one, bit for bit — not `isapprox` to it.
    inner = GI.Polygon([[(2.0, 49.0), (3.0, 49.0), (3.0, 50.0), (2.0, 50.0), (2.0, 49.0)]])
    outer = GI.Polygon([[(0.0, 40.0), (10.0, 40.0), (10.0, 60.0), (0.0, 60.0), (0.0, 40.0)]])
    r = GO.intersection(GO.OverlayNG(S), inner, outer)
    ingested = Set(GO._to_kernel_point(S, p) for p in GI.getpoint(inner))
    @test Set(GI.getpoint(r)) == ingested

    #-- the lon/lat row cannot do this, and that is not a defect in it: `atan`
    #-- and `asin` of a coordinate that already had an exact image in the output
    #-- format are a second rounding, and there is nowhere for it to go but away
    rll = GO.intersection(GO.OverlayNG(S; point_type = Tuple{Float64, Float64}), inner, outer)
    @test Set(GI.getpoint(rll)) != Set(GI.getpoint(inner))
    @test all(p -> p isa Tuple{Float64, Float64}, GI.getpoint(rll))

    #-- what that buys downstream, measured on the quantity the defect was
    #-- reported through. The bar is NOT exact area equality: the result ring
    #-- starts at a different vertex than the input ring, and a spherical area is
    #-- an ordered sum, so a few ULPs of reassociation survive even when every
    #-- coordinate is identical. What the bar can be is the CONTRAST — with the
    #-- coordinate rounding gone, only that reassociation is left (measured: 15
    #-- ULPs), while the lon/lat row carries both (79 ULPs, 5.4x).
    da = abs(GO.area(S, r) - GO.area(S, inner))
    dll = abs(GO.area(S, rll) - GO.area(S, inner))
    @test da < dll
    @test isapprox(GO.area(S, r), GO.area(S, inner); rtol = 1e-14)
end

@testset "point_type: lon/lat output restores the pre-change spherical behaviour" begin
    alg = GO.OverlayNG(S; point_type = Tuple{Float64, Float64})
    for (name, gof, _) in PUBLIC_OPS
        r = gof(alg, AH, BH)
        x = gof(GO.OverlayNG(S), AH, BH)
        @test !GI.is3d(r)
        @test all(p -> p isa Tuple{Float64, Float64}, GI.getpoint(r))
        #-- the two rows describe the same region; only the chart differs
        @test isapprox(GO.area(S, r), GO.area(S, x); rtol = 1e-12)
    end
    #-- and the lon/lat row's result type is the PLANAR result type: the return
    #-- type is a function of `point_type`, not of the manifold
    @test typeof(GO.intersection(alg, A, B)) ===
          typeof(GO.intersection(GO.OverlayNG(), A, B))
end

@testset "point_type: unsupported manifold/point-type pairs are rejected" begin
    #-- the plane has one chart and no xyz reading of it
    @test_throws ArgumentError GO.OverlayNG(P; point_type = USP)
    #-- neither manifold emits anything else at all
    for pt in (Float64, Tuple{Float32, Float32}, Tuple{Float64, Float64, Float64},
               GI.Point{false, false, Tuple{Float64, Float64}, Nothing})
        @test_throws ArgumentError GO.OverlayNG(P; point_type = pt)
        @test_throws ArgumentError GO.OverlayNG(S; point_type = pt)
    end
    #-- the message names the manifold and lists what it does emit
    msg = try; GO.OverlayNG(P; point_type = USP); catch e; sprint(showerror, e); end
    @test occursin("Planar", msg) && occursin("Tuple{Float64,Float64}", msg)
end

@testset "point_type: rebuild re-derives the default and carries a choice" begin
    Rb = GO.GeometryOpsCore.rebuild
    #-- a defaulted algorithm re-derives on the new manifold
    @test Rb(GO.OverlayNG(P), S).point_type === USP
    @test Rb(GO.OverlayNG(S), P).point_type === Tuple{Float64, Float64}
    #-- an explicit non-default choice survives a rebuild that can honour it
    @test Rb(GO.OverlayNG(S; point_type = Tuple{Float64, Float64}), S).point_type ===
          Tuple{Float64, Float64}
    #-- and `exact` still rides along
    @test Rb(GO.OverlayNG(S; exact = GO.False()), P).exact === GO.False()
end

@testset "error paths" begin
    #-- unsupported manifolds are rejected at construction with a clear message
    @test_throws ArgumentError GO.OverlayNG(GO.AutoManifold())
    #-- geometry collections are rejected by the driver
    gc = GI.GeometryCollection([GI.Point((1.0, 1.0))])
    @test_throws ArgumentError GO.intersection(GO.OverlayNG(), gc, A)
    @test_throws ArgumentError GO.intersection(GO.OverlayNG(), A, gc)
    #-- of the Foster–Hormann keywords, `target` IS accepted (see the `target`
    #-- testsets in overlay_ng.jl); `T` and `fix_multipoly` are not — the result
    #-- is always emitted at Float64, and there is nothing to fix
    @test GO.intersection(GO.OverlayNG(), A, B; target = GI.PolygonTrait()) isa AbstractVector
    @test GO.symdifference(GO.OverlayNG(), A, B; target = GI.PolygonTrait()) isa AbstractVector
    @test_throws MethodError GO.intersection(GO.OverlayNG(), A, B, Float32)
    @test_throws MethodError GO.intersection(GO.OverlayNG(), A, B; fix_multipoly = nothing)
end

# ---------------------------------------------------------------------------
# 6. The DEFAULT algorithm is unchanged
# ---------------------------------------------------------------------------

@testset "defaults still dispatch to Foster-Hormann" begin
    #-- This used to compare TARGETED results, which cannot tell the engines
    #-- apart: `intersection(A, B; target = PolygonTrait())` returns the same
    #-- concrete `Vector{GI.Polygon{...}}` from both, so rewiring the default to
    #-- OverlayNG left 98 of 99 assertions in this file passing. The untargeted
    #-- form is the discriminator — Foster–Hormann returns a Vector where
    #-- OverlayNG returns one geometry.
    #-- The discriminator is `GI.Line`: OverlayNG rejects the trait outright, so
    #-- a default that routed there would throw where Foster–Hormann answers.
    ln1 = GI.Line([(0.0, 0.0), (1.0, 1.0)]); ln2 = GI.Line([(0.0, 1.0), (1.0, 0.0)])
    @test_throws ArgumentError GO.intersection(GO.OverlayNG(), ln1, ln2)
    @test length(GO.intersection(ln1, ln2; target = GI.PointTrait())) == 1

    #-- KNOWN DEFECT in the DEFAULT clipping path, pinned here so that it is
    #-- visible rather than merely unexercised. The *untargeted* algorithm-free
    #-- form is the cleanest engine discriminator there is — Foster–Hormann
    #-- returns a Vector where OverlayNG returns one geometry — but it cannot be
    #-- used, because it throws:
    #--
    #--   MethodError: no method matching GeometryOpsCore.TraitTarget(::Nothing)
    #--
    #-- from the target plumbing in `union.jl` / `intersection.jl`, which passes
    #-- `target = nothing` straight into `TraitTarget`. Nothing to do with
    #-- OverlayNG. These flip to "unexpectedly passing" when it is fixed.
    for gof in (GO.intersection, GO.union, GO.difference)
        @test_broken gof(A, B) isa Vector
        @test_broken gof(P, A, B) isa Vector
    end
    #-- the documented `intersection` example from its docstring still holds
    line1 = GI.Line([(124.584961, -12.768946), (126.738281, -17.224758)])
    line2 = GI.Line([(123.354492, -15.961329), (127.22168, -14.008696)])
    ipts = GO.intersection(line1, line2; target = GI.PointTrait())
    @test length(ipts) == 1
    @test all(isapprox.(GI.coordinates(ipts[1]), [125.58375366067548, -14.83572303404496];
                        rtol = 1e-12))
    #-- areas agree between the default engine and OverlayNG on this simple case
    @test isapprox(sum(GO.area, GO.intersection(A, B; target = GI.PolygonTrait())),
                   GO.area(GO.intersection(GO.OverlayNG(), A, B)); rtol = 1e-12)
    @test isapprox(sum(GO.area, GO.union(A, B; target = GI.PolygonTrait())),
                   GO.area(GO.union(GO.OverlayNG(), A, B)); rtol = 1e-12)
    @test isapprox(sum(GO.area, GO.difference(A, B; target = GI.PolygonTrait())),
                   GO.area(GO.difference(GO.OverlayNG(), A, B)); rtol = 1e-12)
end
