using Test
using LinearAlgebra
import GeometryOps as GO, GeoInterface as GI
using GeometryOps.UnitSpherical
using GeometryOpsCore: Planar, Spherical

#=
The lightweight predicates on `Spherical`.

The oracle is RelateNG on the same manifold. Where the two are known to differ
the test says so and pins the reason, rather than asserting the weaker of the
two answers.
=#

# regular n-gon of angular radius `r` (degrees) about (lon0, lat0)
function _ngon(lon0, lat0, r, n; phase = 0.0)
    c = UnitSphereFromGeographic()((lon0, lat0))
    ref = abs(c[3]) < 0.9 ? UnitSphericalPoint(0.0, 0.0, 1.0) : UnitSphericalPoint(1.0, 0.0, 0.0)
    e1 = normalize(cross(c, ref)); e2 = cross(c, e1)
    ρ = deg2rad(r)
    pts = Tuple{Float64,Float64}[]
    for k in 0:(n - 1)
        θ = phase + 2π * k / n
        v = cos(ρ) .* c .+ sin(ρ) .* (cos(θ) .* e1 .+ sin(θ) .* e2)
        p = GeographicFromUnitSphere()(UnitSphericalPoint(normalize(v)))
        push!(pts, (p[1], p[2]))
    end
    push!(pts, pts[1])
    return GI.Polygon([GI.LinearRing(pts)])
end

# alternating radii: genuinely non-convex on the sphere
function _star(lon0, lat0, ro, ri, n)
    c = UnitSphereFromGeographic()((lon0, lat0))
    ref = abs(c[3]) < 0.9 ? UnitSphericalPoint(0.0, 0.0, 1.0) : UnitSphericalPoint(1.0, 0.0, 0.0)
    e1 = normalize(cross(c, ref)); e2 = cross(c, e1)
    pts = Tuple{Float64,Float64}[]
    for k in 0:(2n - 1)
        θ = π * k / n
        ρ = deg2rad(iseven(k) ? ro : ri)
        v = cos(ρ) .* c .+ sin(ρ) .* (cos(θ) .* e1 .+ sin(θ) .* e2)
        p = GeographicFromUnitSphere()(UnitSphericalPoint(normalize(v)))
        push!(pts, (p[1], p[2]))
    end
    push!(pts, pts[1])
    return GI.Polygon([GI.LinearRing(pts)])
end

const TARGET = _ngon(0.0, 0.0, 5.0, 5)
const TVERTS = [(GI.x(p), GI.y(p)) for p in GI.getpoint(GI.getexterior(TARGET))]

@testset "known answers" begin
    inside = _ngon(0.0, 0.0, 2.0, 5)
    around = _ngon(0.0, 0.0, 9.0, 5)
    far    = _ngon(40.0, 0.0, 2.0, 5)

    @test GO.intersects(Spherical(), inside, TARGET)
    @test GO.within(Spherical(), inside, TARGET)
    @test GO.contains(Spherical(), TARGET, inside)
    @test !GO.within(Spherical(), around, TARGET)
    @test GO.contains(Spherical(), around, TARGET)

    @test !GO.intersects(Spherical(), far, TARGET)
    @test GO.disjoint(Spherical(), far, TARGET)
    @test !GO.within(Spherical(), far, TARGET)

    # a polygon is within and contains itself, in either winding
    @test GO.within(Spherical(), TARGET, TARGET)
    @test GO.contains(Spherical(), TARGET, TARGET)
    reversed = GI.Polygon([GI.LinearRing(reverse(TVERTS))])
    @test GO.within(Spherical(), reversed, TARGET)
    @test GO.contains(Spherical(), TARGET, reversed)
end

@testset "a ring's own centre is inside it" begin
    #= The default exterior anchor is the antipode of the vertex mass, which
    puts the ring's centre exactly antipodal to it and the parity walk's test
    arc undefined. Regression for the anchor nudge. =#
    for (lon, lat) in [(0.0, 0.0), (12.3, 41.7), (0.0, 85.0), (179.0, -30.0)]
        ring = _ngon(lon, lat, 3.0, 5)
        @test GO.within(Spherical(), GI.Point((lon, lat)), ring)
        @test GO.intersects(Spherical(), GI.Point((lon, lat)), ring)
    end
end

@testset "points outside a ring are outside it" begin
    #= Regression for the cosine span test, which admitted points on the far
    side of the query and inverted the parity. The pentagon's boundary crosses
    the equator east of the centre at ~4.05 deg, so 4.9 and 5.0 are outside
    even though they are nearer the centre than the vertices are. =#
    for lon in (4.9, 5.0, 6.0, 10.0)
        @test !GO.within(Spherical(), GI.Point((lon, 0.0)), TARGET)
        @test GO.disjoint(Spherical(), GI.Point((lon, 0.0)), TARGET)
    end
    for lon in (0.0, 2.0, 3.5)
        @test GO.within(Spherical(), GI.Point((lon, 0.0)), TARGET)
    end
end

@testset "ring vertices are on the boundary" begin
    for v in TVERTS
        @test GO.intersects(Spherical(), GI.Point(v), TARGET)
        @test !GO.disjoint(Spherical(), GI.Point(v), TARGET)
        @test GO.coveredby(Spherical(), GI.Point(v), TARGET)
        @test !GO.within(Spherical(), GI.Point(v), TARGET)
    end
end

#= The oracle for `within`/`contains`/`covers`/`coveredby` is the DE-9IM matrix,
not `relate_predicate`.

`relate_predicate(RelateNG(Spherical()), pred_within(), a, b)` contradicts the
matrix `relate` computes for the same pair: `relate` returns e.g. `2FF1FF212`,
`relate(alg, a, b, "T*F**F***")` agrees that it matches the within pattern, and
`pred_within` still answers `false`. The planar predicate answers `true` on the
identical matrix, so the defect is in the spherical predicate short-circuit.
Pattern-matching the matrix is the same question without that layer. =#
_within_pattern(a, b) = GO.relate(GO.RelateNG(Spherical()), a, b, "T*F**F***")

@testset "upstream: spherical `pred_within` contradicts its own matrix" begin
    #= Pin the defect so this test starts failing when RelateNG is fixed, at
    which point the workaround above can go. =#
    inner = _ngon(0.0, 0.0, 1.0, 5)
    outer = _star(0.0, 0.0, 8.0, 2.0, 5)
    alg = GO.RelateNG(Spherical())
    @test GO.relate(alg, inner, outer, "T*F**F***")            # matrix says within
    @test GO.within(Spherical(), inner, outer)                 # this path agrees
    @test_broken GO.relate_predicate(alg, GO.pred_within(), inner, outer)
end

@testset "non-convex rings" begin
    #= No separating-axis or convexity shortcut is valid on the sphere: a
    tiling's cells are genuinely non-convex. =#
    star = _star(0.0, 0.0, 8.0, 2.0, 5)
    oracle(p, a, b) = GO.relate_predicate(GO.RelateNG(Spherical()), p(), a, b)
    for other in (_ngon(0.0, 0.0, 1.0, 5), _ngon(6.0, 0.0, 2.0, 5), _star(3.0, 1.0, 6.0, 1.0, 6))
        @test GO.intersects(Spherical(), other, star) == oracle(GO.pred_intersects, other, star)
        @test GO.within(Spherical(), other, star) == _within_pattern(other, star)
    end
end

@testset "degenerate duplicate vertices" begin
    #= Real rings carry coincident vertices; HEALPix rings do near polar
    corners. A zero-length edge must not divide by zero, produce a NaN, or
    change the answer. =#
    base = _ngon(3.0, 1.0, 4.0, 5)
    pts = [(GI.x(p), GI.y(p)) for p in GI.getpoint(GI.getexterior(base))]
    for k in 1:3
        dup = copy(pts); insert!(dup, k + 1, dup[k])
        duped = GI.Polygon([GI.LinearRing(dup)])
        @test GO.intersects(Spherical(), duped, TARGET) == GO.intersects(Spherical(), base, TARGET)
        @test GO.within(Spherical(), duped, TARGET) == GO.within(Spherical(), base, TARGET)
        @test GO.disjoint(Spherical(), duped, TARGET) == GO.disjoint(Spherical(), base, TARGET)
    end
end

@testset "shared edges and vertices" begin
    e1, e2 = TVERTS[1], TVERTS[2]
    m1 = UnitSphereFromGeographic()(e1); m2 = UnitSphereFromGeographic()(e2)
    outward = normalize(m1 .+ m2 .- 0.6 .* UnitSphereFromGeographic()((0.0, 0.0)))
    fg = GeographicFromUnitSphere()(UnitSphericalPoint(outward))
    glued = GI.Polygon([GI.LinearRing([e1, e2, (fg[1], fg[2]), e1])])
    @test GO.intersects(Spherical(), glued, TARGET)
    @test !GO.within(Spherical(), glued, TARGET)
    @test !GO.disjoint(Spherical(), glued, TARGET)
end

@testset "cell-scale geometry" begin
    #= HEALPix L18 cells are ~4e-6 rad across. Three nearly-collinear points on
    a sphere is where an orientation filter is most likely to fail, so this is
    the conditioning that matters, not degree-scale synthetics. =#
    R = 4e-6
    tgt = _ngon(12.3, 41.7, rad2deg(R), 5)
    oracle(p, a, b) = GO.relate_predicate(GO.RelateNG(Spherical()), p(), a, b)
    cases = [
        _ngon(12.3 + rad2deg(R) * 0.6, 41.7, rad2deg(R), 5),
        _ngon(12.3 + rad2deg(R) * 0.6, 41.7, rad2deg(R), 33),
        _ngon(12.3, 41.7, rad2deg(R) * 0.4, 5),
        _ngon(12.3 + rad2deg(R) * 8, 41.7, rad2deg(R), 5),
        _star(12.3, 41.7, rad2deg(R), rad2deg(R) * 0.3, 5),
        tgt,
    ]
    for a in cases
        @test GO.intersects(Spherical(), a, tgt) == oracle(GO.pred_intersects, a, tgt)
        @test GO.within(Spherical(), a, tgt) == GO.relate(GO.RelateNG(Spherical()), a, tgt, "T*F**F***")
        @test GO.contains(Spherical(), a, tgt) == GO.relate(GO.RelateNG(Spherical()), tgt, a, "T*F**F***")
    end
end

@testset "planar behaviour is unchanged" begin
    # the two-argument form must stay planar, CRS or no CRS
    for a in (_ngon(3.0, 1.0, 4.0, 5), _star(0.0, 0.0, 8.0, 2.0, 5), GI.Point((1.0, 1.0)),
              GI.LineString([(-10.0, 0.0), (10.0, 0.0)]))
        for f in (GO.intersects, GO.disjoint, GO.within, GO.contains, GO.covers,
                  GO.coveredby, GO.touches)
            @test f(a, TARGET) === f(Planar(), a, TARGET)
            @test f(a, TARGET) === f(GO.AutoManifold(), a, TARGET)
        end
    end
end

@testset "allocation-free" begin
    #= The point of the path: a yes/no answer on small polygons must not touch
    the heap, because the caller runs tens of thousands of them per query. =#
    a5  = _ngon(3.0, 1.0, 4.0, 5)
    a33 = _ngon(3.0, 1.0, 4.0, 33)
    e1, e2 = TVERTS[1], TVERTS[2]
    m1 = UnitSphereFromGeographic()(e1); m2 = UnitSphereFromGeographic()(e2)
    fg = GeographicFromUnitSphere()(UnitSphericalPoint(normalize(m1 .+ m2 .- 0.6 .* UnitSphereFromGeographic()((0.0, 0.0)))))
    glued = GI.Polygon([GI.LinearRing([e1, e2, (fg[1], fg[2]), e1])])
    far = _ngon(40.0, 0.0, 2.0, 5)

    for a in (a5, a33, glued, far, TARGET)
        for f in (GO.intersects, GO.disjoint, GO.within, GO.coveredby)
            f(Spherical(), a, TARGET)          # compile
            @test (@allocated f(Spherical(), a, TARGET)) == 0
        end
    end
end

@testset "the ring primitives do not allocate" begin
    v = [UnitSphereFromGeographic()(p) for p in TVERTS[1:end-1]]
    q = UnitSphereFromGeographic()((1.0, 1.0))
    anchor = UnitSpherical.spherical_exterior_anchor(v, length(v))
    enc(v, n, q, a) = UnitSpherical.spherical_ring_encloses(v, n, q; anchor = a)
    con(v, n, q) = UnitSpherical.spherical_ring_contains(v, n, q)
    enc(v, length(v), q, anchor); con(v, length(v), q)
    @test (@allocated enc(v, length(v), q, anchor)) == 0
    @test (@allocated con(v, length(v), q)) == 0
end

@testset "`crosses` and `overlaps` refuse a manifold they cannot honour" begin
    a = _ngon(3.0, 1.0, 4.0, 5)
    @test_throws ArgumentError GO.crosses(Spherical(), a, TARGET)
    @test_throws ArgumentError GO.overlaps(Spherical(), a, TARGET)
    # planar still answers, unchanged (on the trait pairs those legacy paths support)
    line = GI.LineString([(-10.0, 0.0), (10.0, 0.0)])
    @test GO.crosses(Planar(), line, TARGET) === GO.crosses(line, TARGET)
    @test GO.overlaps(Planar(), a, TARGET) === GO.overlaps(a, TARGET)
end
