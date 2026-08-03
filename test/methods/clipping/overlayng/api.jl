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
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import ArchGDAL as AG
import GeometryBasics as GB
using GeometryOpsTestHelpers

const P = GO.Planar()
const S = GO.Spherical()

lgc(g) = GI.convert(LG, g)

A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
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

@testset "planar ops through OverlayNG(Planar())" begin
    for (name, gof, lgf) in PUBLIC_OPS
        @test LG.equals(lgc(gof(GO.OverlayNG(P), A, B)), lgf(lgc(A), lgc(B)))
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
    #-- planar disagrees, which is the point of carrying the manifold
    @test GI.trait(GO.intersection(GO.OverlayNG(), PG, pt)) isa GI.MultiPointTrait
end

# ---------------------------------------------------------------------------
# 3. `symdifference` — the new public function and its engine decision
# ---------------------------------------------------------------------------

@testset "symdifference entry points" begin
    expected = LG.symmetricDifference(lgc(A), lgc(B))
    #-- explicit algorithm, manifold-only, and bare forms all agree
    @test LG.equals(lgc(GO.symdifference(GO.OverlayNG(), A, B)), expected)
    @test LG.equals(lgc(GO.symdifference(P, A, B)), expected)
    @test LG.equals(lgc(GO.symdifference(A, B)), expected)
    #-- the manifold form routes to OverlayNG with that manifold
    @test isapprox(GO.area(S, GO.symdifference(S, A, B)),
                   GO.area(S, GO.symdifference(GO.OverlayNG(S), A, B)); rtol = 1e-14)
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
    #-- the keyword and positional-manifold constructors agree
    @test GO.OverlayNG(S) == GO.OverlayNG(; manifold = S)
    #-- `exact = False()` still runs (Float64 filters only) and agrees here
    @test LG.equals(lgc(GO.intersection(GO.OverlayNG(; exact = GO.False()), A, B)),
                    LG.intersection(lgc(A), lgc(B)))
end

@testset "error paths" begin
    #-- unsupported manifolds are rejected at construction with a clear message
    @test_throws ArgumentError GO.OverlayNG(GO.AutoManifold())
    #-- geometry collections are rejected by the driver
    gc = GI.GeometryCollection([GI.Point((1.0, 1.0))])
    for (_, gof, _) in PUBLIC_OPS
        @test_throws ArgumentError gof(GO.OverlayNG(), gc, A)
        @test_throws ArgumentError gof(GO.OverlayNG(), A, gc)
    end
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
    #-- the algorithm-free forms keep returning a Vector of target geometries,
    #-- which is the Foster–Hormann contract (OverlayNG returns one geometry)
    for gof in (GO.intersection, GO.union, GO.difference)
        v = gof(A, B; target = GI.PolygonTrait())
        @test v isa Vector
        @test all(g -> GI.trait(g) isa GI.PolygonTrait, v)
        #-- and so does the manifold-only form
        @test gof(P, A, B; target = GI.PolygonTrait()) isa Vector
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
