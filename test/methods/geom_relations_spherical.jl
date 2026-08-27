using Test
using LinearAlgebra
using Random
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
    #= Fixed upstream by #472. Kept as a regression test: the three answers must
    agree, and `_within_pattern` above is now redundant with `pred_within`. =#
    inner = _ngon(0.0, 0.0, 1.0, 5)
    outer = _star(0.0, 0.0, 8.0, 2.0, 5)
    alg = GO.RelateNG(Spherical())
    @test GO.relate(alg, inner, outer, "T*F**F***")            # matrix says within
    @test GO.within(Spherical(), inner, outer)                 # this path agrees
    @test GO.relate_predicate(alg, GO.pred_within(), inner, outer)
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

#= Irregular rings, in the shape real data comes in: unequal edge lengths and
turn angles, so a vertex's neighbourhood is not the tidy one a regular n-gon
gives. `Random` is only used to lay the vertices out, from a fixed seed. =#
function _irregular(lon0, lat0, r, n; seed)
    rng = Random.MersenneTwister(seed)
    c = UnitSphereFromGeographic()((lon0, lat0))
    ref = abs(c[3]) < 0.9 ? UnitSphericalPoint(0.0, 0.0, 1.0) : UnitSphericalPoint(1.0, 0.0, 0.0)
    e1 = normalize(cross(c, ref)); e2 = cross(c, e1)
    pts = Tuple{Float64,Float64}[]
    for k in 0:(n - 1)
        θ = 2π * (k + 0.4 * (rand(rng) - 0.5)) / n
        ρ = deg2rad(r * (0.55 + 0.9 * rand(rng)))
        v = cos(ρ) .* c .+ sin(ρ) .* (cos(θ) .* e1 .+ sin(θ) .* e2)
        p = GeographicFromUnitSphere()(UnitSphericalPoint(normalize(v)))
        push!(pts, (p[1], p[2]))
    end
    push!(pts, pts[1])
    return GI.Polygon([GI.LinearRing(pts)])
end

@testset "a ring is within itself, and contains itself" begin
    #= Every edge of a ring tested against its own ring is a shared edge, so
    this is the densest shared-edge case there is, and the one real tilings
    hit constantly.

    Regression for two faults in the segment-splitting walk. The walk ordered
    split points by `dot(A, ·)` from an assumed starting ordinate of exactly
    `1`, but `dot(A, A)` is a rounded sum of three squares and lands an ulp
    below it, so the segment's own start vertex sorted as strictly ahead of
    the start and opened a zero-length piece. That piece was then classified
    by normalizing `p + p`, which moves the point by an ulp — and an ulp past
    a shared vertex is off the end of both arcs meeting there, so a ring
    vertex classified as OUTSIDE the ring it belongs to, and the ring stopped
    being within itself.

    A regular n-gon mostly hides this (the perturbed point stays inside a
    neighbouring edge's band); irregular rings do not. =#
    for n in (5, 8, 12, 20, 38), seed in 1:12
        ring = _irregular(12.3, 41.7, 4.0, n; seed)
        @test GO.within(Spherical(), ring, ring)
        @test GO.contains(Spherical(), ring, ring)
        @test GO.coveredby(Spherical(), ring, ring)
        @test !GO.disjoint(Spherical(), ring, ring)
    end
    # and at cell scale, where the vertices are 4e-6 rad apart
    for n in (5, 12, 38), seed in 1:8
        ring = _irregular(12.3, 41.7, rad2deg(4e-6), n; seed)
        @test GO.within(Spherical(), ring, ring)
    end
end

@testset "a multipolygon contains its own components" begin
    parts = [_irregular(10.0 + 30.0 * i, 20.0, 3.0, 12; seed = 7 + i) for i in 0:3]
    mp = GI.MultiPolygon(parts)
    for part in parts
        @test GO.contains(Spherical(), mp, part)
        @test GO.within(Spherical(), part, mp)
        @test GO.intersects(Spherical(), mp, part)
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
                  GO.coveredby, GO.touches, GO.crosses)
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

@testset "`overlaps` refuses a manifold it cannot honour" begin
    a = _ngon(3.0, 1.0, 4.0, 5)
    @test_throws ArgumentError GO.overlaps(Spherical(), a, TARGET)
    # planar still answers, unchanged (on the trait pairs that legacy path supports)
    @test GO.overlaps(Planar(), a, TARGET) === GO.overlaps(a, TARGET)
end

#= `crosses` runs on the shared processors now, so `Spherical()` reaches the
same four geometric leaves every other lightweight predicate does. The oracle
is RelateNG on the same manifold. =#
@testset "`crosses` on the sphere" begin
    tgt = _ngon(12.3, 41.7, 4.0, 5)
    oracle(a, b) = GO.relate_predicate(GO.RelateNG(Spherical()), GO.pred_crosses(), a, b)

    #= Coordinates are deliberately irregular. RelateNG's exact spherical
    crossing point divides by zero when two arcs meet on a coordinate plane —
    a great circle crossing the equator at longitude 0, say — so a tidy
    fixture cannot be cross-checked against it at all. That defect is upstream
    of this path and unrelated to `crosses`. =#
    through   = GI.LineString([(3.7, 39.1), (21.9, 44.3)])
    inside    = GI.LineString([(11.9, 41.4), (12.8, 42.0)])
    outside   = GI.LineString([(40.0, 10.0), (50.0, 12.0)])
    half_in   = GI.LineString([(12.3, 41.7), (31.4, 44.9)])
    bent      = GI.LineString([(3.7, 39.1), (12.86, 41.86), (21.9, 44.3)])
    crossing  = GI.LineString([(9.1, 47.2), (17.3, 36.4)])
    tee       = GI.LineString([(12.4, 47.9), (12.86, 41.86)])
    abutting  = GI.LineString([(21.9, 44.3), (33.1, 46.8)])

    @testset "line against polygon" begin
        for (line, expected) in ((through, true), (inside, false), (outside, false), (half_in, true))
            @test GO.crosses(Spherical(), line, tgt) == expected
            @test GO.crosses(Spherical(), line, tgt) == oracle(line, tgt)
            # and the transpose, `T*****T**`
            @test GO.crosses(Spherical(), tgt, line) == expected
            @test GO.crosses(Spherical(), tgt, line) == oracle(tgt, line)
        end
    end

    @testset "ring against polygon" begin
        for (ring, expected) in ((GI.getexterior(_ngon(15.1, 42.4, 3.3, 5)), true),
                                 (GI.getexterior(_ngon(12.3, 41.7, 1.1, 5)), false),
                                 (GI.getexterior(_ngon(60.1, 12.4, 3.3, 5)), false))
            @test GO.crosses(Spherical(), ring, tgt) == expected
            @test GO.crosses(Spherical(), ring, tgt) == oracle(ring, tgt)
        end
    end

    @testset "curve against curve" begin
        # interiors meeting in one point: `0********`
        @test GO.crosses(Spherical(), through, crossing) == true
        @test GO.crosses(Spherical(), through, crossing) == oracle(through, crossing)
        #= A vee whose apex sits on the arc's interior touches without passing
        through, which is still a zero-dimensional meeting of the interiors.
        The apex is the normalized midpoint of the arc's endpoints, so it is on
        the great circle by construction rather than by luck. =#
        apex = GeographicFromUnitSphere()(UnitSphericalPoint(normalize(
            UnitSphereFromGeographic()((3.7, 39.1)) .+ UnitSphereFromGeographic()((21.9, 44.3)))))
        vee = GI.LineString([(8.0, 55.0), (apex[1], apex[2]), (18.0, 55.0)])
        @test GO.crosses(Spherical(), vee, through) == true
        @test GO.crosses(Spherical(), vee, through) == oracle(vee, through)
        #= The tee's own endpoint is what lands on the other curve, so the
        interiors never meet and `0********` fails. =#
        @test GO.crosses(Spherical(), tee, bent) == false
        @test GO.crosses(Spherical(), tee, bent) == oracle(tee, bent)
        # a shared arc is one-dimensional, so it is not a crossing
        @test GO.crosses(Spherical(), through, through) == false
        @test GO.crosses(Spherical(), through, through) == oracle(through, through)
        # boundary-to-boundary contact leaves the interiors disjoint
        @test GO.crosses(Spherical(), through, abutting) == false
        @test GO.crosses(Spherical(), through, abutting) == oracle(through, abutting)
    end

    @testset "multipoint against curve and polygon" begin
        in_out  = GI.MultiPoint([(12.3, 41.7), (41.9, 3.2)])
        both_in = GI.MultiPoint([(12.3, 41.7), (12.9, 41.9)])
        on_end  = GI.MultiPoint([(3.7, 39.1), (41.9, 3.2)])
        on_vert = GI.MultiPoint([(12.86, 41.86), (41.9, 3.2)])
        @test GO.crosses(Spherical(), in_out, tgt) == true
        @test GO.crosses(Spherical(), in_out, tgt) == oracle(in_out, tgt)
        @test GO.crosses(Spherical(), both_in, tgt) == false
        @test GO.crosses(Spherical(), both_in, tgt) == oracle(both_in, tgt)
        # a point on the curve's boundary is neither interior to it nor exterior
        @test GO.crosses(Spherical(), on_end, bent) == false
        @test GO.crosses(Spherical(), on_end, bent) == oracle(on_end, bent)
        @test GO.crosses(Spherical(), on_vert, bent) == true
        @test GO.crosses(Spherical(), on_vert, bent) == oracle(on_vert, bent)
    end

    @testset "equal dimensions, and single points" begin
        other = _ngon(15.1, 42.4, 3.3, 5)
        @test GO.crosses(Spherical(), other, tgt) == false
        @test GO.crosses(Spherical(), other, tgt) == oracle(other, tgt)
        @test GO.crosses(Spherical(), GI.Point((12.3, 41.7)), tgt) == false
        @test GO.crosses(Spherical(), GI.Point((12.3, 41.7)), tgt) == oracle(GI.Point((12.3, 41.7)), tgt)
    end

    @testset "the manifold changes the answer" begin
        #= A 120-degree arc at latitude 40.2 bulges to about 59.4 degrees at
        its midpoint longitude, so a small cap sitting on the bulge is crossed
        on the sphere and missed entirely by the planar segment. =#
        arc = GI.LineString([(-61.3, 40.2), (58.7, 40.2)])
        for lat in (59.0, 59.4, 60.0)
            cap = _ngon(-1.3, lat, 1.2, 5)
            @test GO.crosses(Spherical(), arc, cap) == true
            @test GO.crosses(Spherical(), arc, cap) == oracle(arc, cap)
            @test GO.crosses(Planar(), arc, cap) == false
        end
        # below the bulge neither manifold finds a crossing
        for lat in (57.0, 58.0)
            cap = _ngon(-1.3, lat, 1.2, 5)
            @test GO.crosses(Spherical(), arc, cap) == false
            @test GO.crosses(Spherical(), arc, cap) == oracle(arc, cap)
            @test GO.crosses(Planar(), arc, cap) == false
        end
    end
end
