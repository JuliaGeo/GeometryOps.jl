# End-to-end `relate(RelateNG(; manifold = Spherical()), …)` smoke tests
# (Task 14). With the spherical kernel and the manifold-derived point type
# in place, the engine runs on the unindexed NestedLoop path. Away from the
# poles and the antimeridian, spherical and planar topology agree (the spike
# measured 2000/2000), so the patch test compares the two DE-9IM matrices
# directly; the pole-containing ring is a case only the spherical kernel
# gets right.

using Test
import GeometryOps as GO
import GeometryOps: Spherical, RelateNG
import GeoInterface as GI

alg = RelateNG(; manifold = Spherical())

@testset "spherical relate agrees with planar on a small mid-latitude patch" begin
    # two overlapping boxes near (10°, 45°): away from poles/antimeridian the
    # spherical and planar topology agree
    A = GI.Polygon([GI.LinearRing([(0., 40.), (10., 40.), (10., 50.), (0., 50.), (0., 40.)])])
    B = GI.Polygon([GI.LinearRing([(5., 45.), (15., 45.), (15., 55.), (5., 55.), (5., 45.)])])
    @test GO.relate(alg, A, B) == GO.relate(A, B)           # same DE-9IM
    @test GO.relate(alg, A, B, "T*T***T**")                  # overlaps
end

@testset "spherical handles a pole-containing ring" begin
    cap = GI.Polygon([GI.LinearRing([(0., 80.), (120., 80.), (240., 80.), (0., 80.)])])
    pt  = GI.Point(0., 89.)                                  # near north pole
    @test GO.relate(alg, cap, pt, "T*****FF*")               # contains
end

# Task 15: the tree accelerator must build 3D great-circle arc extents, not a
# 2D coordinate box. A long near-equatorial arc (lon 0 → 170) bulges to y ≈ 1
# at lon 90 while its endpoint box has y ∈ [0, 0.17]; a polygon straddling the
# equator at lon 90 crosses it there. A 2D endpoint box prunes that pair away
# (wrong DE-9IM); the bulge-aware `spherical_arc_extent` keeps it.
@testset "spherical tree accelerator agrees with NestedLoop (arc bulge)" begin
    A = GI.Polygon([GI.LinearRing([(0., 0.), (170., 0.), (85., 40.), (0., 0.)])])
    B = GI.Polygon([GI.LinearRing([(88., -2.), (92., -2.), (92., 2.), (88., 2.), (88., -2.)])])
    tree = RelateNG(; manifold = Spherical(), accelerator = GO.DoubleNaturalTree())
    loop = RelateNG(; manifold = Spherical(), accelerator = GO.NestedLoop())
    @test GO.relate(tree, A, B) == GO.relate(loop, A, B)
end

# Task 15: above the 32-segment threshold AutoAccelerator picks the tree path;
# it must give the same DE-9IM as the unindexed nested loop on a real ring.
@testset "spherical AutoAccelerator picks the tree above threshold" begin
    n = 48                                                   # 48 segments > threshold
    ringpts = [(10.0 + 8cosd(t), 45.0 + 5sind(t)) for t in range(0, 360; length = n + 1)]
    A = GI.Polygon([GI.LinearRing(ringpts)])
    B = GI.Polygon([GI.LinearRing([(8., 43.), (20., 43.), (20., 52.), (8., 52.), (8., 43.)])])
    loop = RelateNG(; manifold = Spherical(), accelerator = GO.NestedLoop())
    @test GO.relate(alg, A, B) == GO.relate(loop, A, B)
end

# Task 17: prepared spherical relate (A indexed once, in 3D) must agree with
# the unprepared nested-loop relate over several B geometries. The prepared
# edge index is the dimension-generic natural-order `RTree` over 3D arc extents.
@testset "spherical prepared relate agrees with unprepared" begin
    n = 48
    ringpts = [(10.0 + 8cosd(t), 45.0 + 5sind(t)) for t in range(0, 360; length = n + 1)]
    A = GI.Polygon([GI.LinearRing(ringpts)])
    prep = GO.prepare(alg, A)
    loop = RelateNG(; manifold = Spherical(), accelerator = GO.NestedLoop())
    Bs = (
        GI.Polygon([GI.LinearRing([(8., 43.), (20., 43.), (20., 52.), (8., 52.), (8., 43.)])]),
        GI.Polygon([GI.LinearRing([(11., 45.), (13., 45.), (13., 47.), (11., 47.), (11., 45.)])]),
        GI.Point(10., 45.),
        GI.Point(40., 80.),
    )
    for B in Bs
        @test GO.relate(prep, B) == GO.relate(loop, A, B)
    end
end

# Ring winding must not change which region a polygon bounds: real-world data
# arrives in either winding (Natural Earth ships shapefile-convention CW
# shells), and the kernel contract for `rk_point_in_ring` is the area enclosed
# by the ring — the region smaller than a hemisphere — like the planar
# ray-crossing parity. Under the previous S2 left-of-ring reading, a CW
# continent meant "the whole sphere minus the continent" and intersected
# everything.
@testset "polygon interior is winding-independent" begin
    sq(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s), (x0, y0)]
    ccw(pts) = GI.Polygon([GI.LinearRing(pts)])
    cw(pts) = GI.Polygon([GI.LinearRing(reverse(pts))])
    canada_ish = sq(-100.0, 50.0, 10.0)
    australia_ish = sq(140.0, -35.0, 15.0)
    continent = sq(-140.0, 20.0, 60.0)                       # continent-sized
    far_pt = GI.Point(150.0, -25.0)
    @testset "disjoint stays disjoint in all four winding combinations" begin
        for A in (ccw(canada_ish), cw(canada_ish)), B in (ccw(australia_ish), cw(australia_ish))
            @test !GO.intersects(alg, A, B)
            @test GO.disjoint(alg, A, B)
        end
        for A in (ccw(continent), cw(continent))
            @test !GO.intersects(alg, A, far_pt)
            @test !GO.intersects(alg, A, ccw(australia_ish))
        end
    end
    @testset "containment in all four winding combinations" begin
        inner = sq(-97.0, 53.0, 2.0)
        for A in (ccw(canada_ish), cw(canada_ish)), B in (ccw(inner), cw(inner))
            @test GO.contains(alg, A, B)
            @test GO.within(alg, B, A)
            @test !GO.touches(alg, A, B)
        end
    end
    @testset "shared-edge neighbors touch in all four winding combinations" begin
        left = sq(0.0, 40.0, 10.0)
        right = sq(10.0, 40.0, 10.0)
        for A in (ccw(left), cw(left)), B in (ccw(right), cw(right))
            @test GO.touches(alg, A, B)
            @test GO.intersects(alg, A, B)
        end
    end
    @testset "equator-symmetric ring (planar isCCW tie case)" begin
        # two vertices tie for the extreme xyz-y coordinate, which broke the
        # planar extreme-vertex orientation test on the sphere: both windings
        # read CW
        eq = sq(-5.0, -5.0, 10.0)
        inside = GI.Point(0.0, 0.0)
        for A in (ccw(eq), cw(eq))
            @test GO.contains(alg, A, inside)
            @test !GO.intersects(alg, A, far_pt)
        end
    end
end

# A degenerate zero-area sliver ring with a repeated vertex (NE 110m North
# Korea ships an `[A, A, B, A]` first polygon) must not swallow the sphere:
# the zero-length arc `A → A` used to read every query point as on-boundary.
@testset "degenerate sliver polygon in a multipolygon stays local" begin
    a = (130.78, 42.22)
    sliver = GI.LinearRing([a, a, (130.78002, 42.22), a])
    body = GI.LinearRing([(124., 38.), (130., 38.), (130., 43.), (124., 43.), (124., 38.)])
    mp = GI.MultiPolygon([GI.Polygon([sliver]), GI.Polygon([body])])
    far = GI.Polygon([GI.LinearRing([(-60., -25.), (-55., -25.), (-55., -20.), (-60., -20.), (-60., -25.)])])
    @test !GO.intersects(alg, mp, far)
    @test !GO.intersects(alg, mp, GI.Point(-60.0, -20.0))
    @test GO.intersects(alg, mp, GI.Point(127.0, 40.0))
end

# The winding-independent region also drives the kernel interaction bounds: a
# CW polar cap must still bound the cap (its enclosed pole included), not the
# complement — which would under-cover the pole and prune away a contained
# point.
@testset "CW polar cap still means the cap" begin
    cap_cw = GI.Polygon([GI.LinearRing([(0., 80.), (240., 80.), (120., 80.), (0., 80.)])])
    pt = GI.Point(0., 89.)                                   # near north pole
    @test GO.relate(alg, cap_cw, pt, "T*****FF*")            # contains
    e = GO.rk_interaction_bounds(Spherical(), cap_cw)
    @test e.Z[2] >= 1.0                                      # bounds reach the enclosed pole
end

# Prepared point location goes through the indexed spherical locator (the
# longitude-interval edge index with meridian-arc crossing parity against a
# pole anchor); unprepared location is the exact ring scan. End to end the
# two must produce the same DE-9IM for every point, including the awkward
# geometries: an antimeridian-straddling box (split index intervals), polar
# caps over either pole (enclosed pole, polar-axis queries), and a shell
# with the south pole on its boundary (the parity anchor falls back to the
# north pole). The grid deliberately includes both poles, the ±180° seam,
# and points on ring boundaries.
@testset "prepared point location agrees with unprepared" begin
    rings = (
        [(170., -10.), (-170., -10.), (-170., 10.), (170., 10.), (170., -10.)],  # antimeridian box
        [(0., 80.), (120., 80.), (240., 80.), (0., 80.)],                        # north polar cap
        [(0., -80.), (120., -80.), (240., -80.), (0., -80.)],                    # south polar cap
        [(0., -90.), (20., -60.), (-20., -60.), (0., -90.)],                     # south pole on boundary
    )
    qs = [GI.Point(lon, lat) for lon in -180.0:30.0:180.0 for lat in -90.0:15.0:90.0]
    push!(qs, GI.Point(180.0, 5.0), GI.Point(175.0, 0.0), GI.Point(0.0, -70.0))
    for pts in rings, w in (pts, reverse(pts))
        A = GI.Polygon([GI.LinearRing(w)])
        prep = GO.prepare(alg, A)
        n_mismatch = count(q -> GO.relate(prep, q) != GO.relate(alg, A, q), qs)
        @test n_mismatch == 0
    end
end

# A ring may carry antipodal VERTEX pairs even though antipodal edges are
# rejected at ingest — the AntipodalEdgeSplit output is exactly such a ring.
# The orientation sum must not connect non-adjacent vertices (a chord between
# an antipodal pair has no defined geodesic; the one-apex Girard fan
# degenerated on exactly that, flipping the shared orientation bit and with
# it the polygon interior — CI regression on the fixed ring below. The S2
# turning-angle curvature uses adjacent vertices only.)
@testset "antipodal vertex pair does not flip the interior" begin
    # AntipodalEdgeSplit's output for the (0,0) → (180,0) edge, verbatim:
    # the first vertex (0,0) is antipodal to the mid-ring vertex (180,0)
    split_pts = [(0., 0.), (90., 0.), (180., 0.), (90., 80.), (0., 0.)]
    inside = GI.Point(10., 10.)      # under the (0,0)→(90,80) arc (lat 44.6° at lon 10)
    outside = GI.Point(10., -10.)    # southern hemisphere
    for pts in (split_pts, reverse(split_pts))
        poly = GI.Polygon([GI.LinearRing(pts)])
        @test GO.relate(alg, poly, inside, "T*****FF*")      # contains
        @test GO.contains(alg, poly, inside)
        @test !GO.intersects(alg, poly, outside)
    end
    # a different configuration: the first vertex's antipode mid-ring, with
    # the ring elsewhere entirely off the first vertex's great circles
    pts2 = [(20., 30.), (100., 10.), (-160., -30.), (-100., 10.), (20., 30.)]
    inside2 = GI.Point(0., 40.)
    for pts in (pts2, reverse(pts2))
        poly = GI.Polygon([GI.LinearRing(pts)])
        @test GO.contains(alg, poly, inside2)
        @test !GO.intersects(alg, poly, GI.Point(0., -80.))
    end
end

# `Spherical(; oriented = true)`: ring directions are declared correct
# (exterior rings CCW, interior rings CW), so the polygon interior is the
# region on the LEFT of each ring's stored order (S2 InitOriented semantics)
# — no winding is computed, and regions larger than a hemisphere become
# representable. The default manifold keeps the enclosed-region reading
# tested above.
@testset "oriented manifold: winding is authoritative" begin
    oalg = RelateNG(; manifold = Spherical(oriented = true))
    sq(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s), (x0, y0)]
    ccw(pts...) = GI.Polygon([GI.LinearRing(p) for p in pts])

    @testset "correctly wound data matches the unoriented mode" begin
        A = ccw(sq(0., 40., 10.))
        B = ccw(sq(5., 45., 10.))
        @test GO.relate(oalg, A, B) == GO.relate(alg, A, B)
        @test GO.contains(oalg, A, ccw(sq(2., 42., 2.)))
        @test GO.touches(oalg, A, ccw(sq(10., 40., 10.)))
        @test !GO.intersects(oalg, A, GI.Point(-90., -45.))
    end

    @testset "a CW shell denotes its complement" begin
        Q = sq(0., 40., 10.)
        C = ccw(reverse(Q))                       # the sphere minus the square
        far = GI.Point(-90., -45.)
        inside_sq = GI.Point(5., 45.)
        @test GO.contains(oalg, C, far)
        @test !GO.contains(oalg, C, inside_sq)
        @test !GO.intersects(oalg, C, inside_sq)
        @test GO.intersects(oalg, C, GI.Point(0., 45.))   # boundary point
        # the square and its complement share a boundary and have disjoint
        # interiors: they touch, and neither contains the other
        @test GO.touches(oalg, ccw(Q), C)
        @test !GO.contains(oalg, C, ccw(Q))
        # the same CW ring under the default manifold is still the square
        @test GO.contains(alg, C, inside_sq)
        @test !GO.intersects(alg, C, far)
    end

    @testset "region larger than a hemisphere: sphere minus a polar cap" begin
        cap_cw = [(0., 80.), (240., 80.), (120., 80.), (0., 80.)]
        S = ccw(cap_cw)
        @test GO.contains(oalg, S, GI.Point(0., -89.))    # opposite pole is interior
        @test !GO.intersects(oalg, S, GI.Point(0., 89.))  # the cap is excluded
        # containment of a whole mid-latitude polygon in the complement
        eq_box = ccw(sq(0., -5., 10.))
        @test GO.contains(oalg, S, eq_box)
        @test GO.within(oalg, eq_box, S)
        # a polygon straddling the cap boundary overlaps it
        straddle = ccw(sq(-5., 75., 10.))
        @test GO.relate(oalg, S, straddle, "T*T***T**")   # overlaps
        # default manifold: the same ring is the cap
        @test GO.contains(alg, S, GI.Point(0., 89.))
        @test !GO.intersects(alg, S, GI.Point(0., -89.))
    end

    @testset "hole windings in both modes" begin
        shell = sq(0., 30., 20.)
        hole = sq(5., 35., 5.)
        in_ring = GI.Point(2., 32.)    # inside the shell, outside the hole
        in_hole = GI.Point(7., 37.)    # inside the hole cavity
        donut = ccw(shell, reverse(hole))          # CCW shell + CW hole: the convention
        @test GO.contains(oalg, donut, in_ring)
        @test !GO.contains(oalg, donut, in_hole)
        @test !GO.intersects(oalg, donut, in_hole)
        # a CCW-wound "hole" removes the cavity's complement instead,
        # leaving only the disk: (shell ∩ disk) = disk
        weird = ccw(shell, hole)
        @test GO.contains(oalg, weird, in_hole)
        @test !GO.intersects(oalg, weird, in_ring)
        # CW shell with a CW hole elsewhere: (sphere − Q) − R
        Q, R = shell, sq(100., -20., 10.)
        comp = ccw(reverse(Q), reverse(R))
        @test GO.contains(oalg, comp, GI.Point(-90., -45.))
        @test !GO.intersects(oalg, comp, in_ring)               # inside Q
        @test !GO.intersects(oalg, comp, GI.Point(105., -15.))  # inside R
        # the default manifold reads every winding combination as the
        # enclosed regions (shapefile-convention CW shell + CCW hole too)
        for s in (shell, reverse(shell)), h in (hole, reverse(hole))
            shp = ccw(s, h)
            @test GO.contains(alg, shp, in_ring)
            @test !GO.contains(alg, shp, in_hole)
            @test GO.relate(alg, shp, donut, "T*F**FFF*")       # topologically equal
        end
    end

    @testset "oriented prepared point location agrees with unprepared" begin
        polys = (
            ccw([(0., 80.), (120., 80.), (240., 80.), (0., 80.)]),           # polar cap
            ccw([(0., 80.), (240., 80.), (120., 80.), (0., 80.)]),           # sphere minus the cap
            ccw([(170., -10.), (-170., -10.), (-170., 10.), (170., 10.), (170., -10.)]),  # antimeridian box
            ccw(sq(0., 30., 20.), reverse(sq(5., 35., 5.))),                 # donut
            ccw(reverse(sq(0., 30., 20.))),                                  # complement of a box
        )
        qs = [GI.Point(lon, lat) for lon in -180.0:30.0:180.0 for lat in -90.0:15.0:90.0]
        push!(qs, GI.Point(180.0, 5.0), GI.Point(7.0, 37.0), GI.Point(2.0, 32.0))
        for A in polys
            prep = GO.prepare(oalg, A)
            n_mismatch = count(q -> GO.relate(prep, q) != GO.relate(oalg, A, q), qs)
            @test n_mismatch == 0
        end
    end
end

# Task 18: an exactly-antipodal edge has no unique great-circle arc; the kernel
# refuses it at ingest with a message pointing at the AntipodalEdgeSplit remedy.
@testset "spherical antipodal edge throws informatively" begin
    p0   = GO._to_kernel_point(Spherical(), (0., 0.))       # (1, 0, 0)
    p180 = GO._to_kernel_point(Spherical(), (180., 0.))     # (-1, 0, 0): antipodal
    p90  = GO._to_kernel_point(Spherical(), (90., 0.))      # (0, 1, 0)
    err = try; GO._validate_relate_edges(Spherical(), GI.LineString([p0, p180])); nothing; catch e; e; end
    @test err isa ArgumentError
    @test occursin("AntipodalEdgeSplit", err.msg)
    #-- a normal edge and a repeated (zero-length) vertex do NOT throw
    @test (GO._validate_relate_edges(Spherical(), GI.LineString([p0, p90])); true)
    @test (GO._validate_relate_edges(Spherical(), GI.LineString([p0, p0])); true)
    #-- the whole relate rejects a polygon carrying an antipodal edge at ingest
    bad = GI.Polygon([GI.LinearRing([(0., 0.), (180., 0.), (90., 80.), (0., 0.)])])
    @test_throws ArgumentError GO.relate(alg, bad, GI.Point(10., 10.))
end

# Kernel points are 3D, so every coordinate comparison in the segment-string
# machinery must compare all three components. Two vertices at mirror
# latitudes on one meridian — e.g. (0, -5) and (0, 5) — share x and y on the
# unit sphere and differ only in z; the planar 2D coordinate equality
# (`_equals2`, Java's equals2D) read such a pair as a repeated point and
# DELETED the second vertex at extraction, corrupting edge topology for any
# ring with an equator-mirrored meridian edge.
@testset "equator-mirrored meridian edges survive extraction" begin
    left = GI.Polygon([GI.LinearRing(
        [(-10., -5.), (0., -5.), (0., 5.), (-10., 5.), (-10., -5.)])])
    right = GI.Polygon([GI.LinearRing(
        [(0., -5.), (10., -5.), (10., 5.), (0., 5.), (0., -5.)])])
    #-- the extracted shell must keep all its vertices
    rg = GO.RelateGeometry(Spherical(), left; exact = GO.True())
    ss = GO.extract_segment_strings(rg, GO.GEOM_A, nothing)
    @test length(only(ss).pts) == 5
    #-- shared-meridian neighbors: same DE-9IM as the planar engine
    @test GO.touches(alg, left, right)
    @test string(GO.relate(alg, left, right)) == string(GO.relate(left, right))
    inner = GI.Polygon([GI.LinearRing(
        [(-8., -3.), (-2., -3.), (-2., 3.), (-8., 3.), (-8., -3.)])])
    @test GO.contains(alg, left, inner)
end

# Containment parity is anchored at a DEFINITIONALLY exterior point (the
# antipode of the ring's vertex mass), so a ring that self-intersects on the
# sphere — the class `prepare` validation rejects — degrades to even-odd
# semantics instead of inverting globally: both lobes of a figure-eight read
# IN, everything else OUT, matching planar even-odd ray crossing (and what
# S2Builder's undirected repair produces). The previous wedge-plus-winding
# bootstrap answered with whichever lobe hosted the anchor edge, which on
# real data (NE 110m Sudan) inverted every containment answer on the globe.
@testset "self-crossing rings degrade even-odd, not inverted" begin
    palg = RelateNG()
    #-- an explicit self-crossing quadrilateral bowtie (crossing near (0,0))
    #-- and an asymmetric variant (crossing near (-5.8, 0.3), unequal lobes)
    bowtie = [(-10., -10.), (10., 10.), (10., -10.), (-10., 10.), (-10., -10.)]
    skew = [(-20., -5.), (20., 10.), (20., -10.), (-20., 6.), (-20., -5.)]
    cases = (
        (bowtie, (8., 1.), (-8., -1.), (0., 5.), (0., -5.)),
        (skew, (15., 0.), (-17., 0.), (0., 8.), (0., -8.)),
    )
    far = GI.Point(100., 40.)
    for (pts, in_a, in_b, out_a, out_b) in cases, w in (pts, reverse(pts))
        poly = GI.Polygon([GI.LinearRing(w)])
        #-- even-odd: both lobes IN, between/far OUT, in both windings
        @test GO.contains(alg, poly, GI.Point(in_a))
        @test GO.contains(alg, poly, GI.Point(in_b))
        @test !GO.intersects(alg, poly, GI.Point(out_a))
        @test !GO.intersects(alg, poly, GI.Point(out_b))
        @test !GO.intersects(alg, poly, far)
        #-- and exact agreement with planar even-odd on the same probes
        probes = (GI.Point(in_a), GI.Point(in_b), GI.Point(out_a), GI.Point(out_b), far)
        for q in probes
            @test GO.intersects(alg, poly, q) == GO.intersects(palg, poly, q)
        end
        #-- the indexed locator (prepared with validation bypassed — the
        #-- build-time pole-anchor classification composes with the same
        #-- parity) must agree with the unprepared exact scan
        prep = GO.prepare(alg, poly; validate = false)
        for q in probes
            @test GO.relate(prep, q) == GO.relate(alg, poly, q)
        end
    end
end

# The definitional anchor degenerates for near-hemisphere/vertex-symmetric
# rings (vertex mass ~ 0): those queries fall back to the wedge-plus-winding
# bootstrap — for such rings the enclosed/complement distinction is itself
# near-degenerate, and the fallback keeps the pre-existing behavior: an
# exact-equator ring (every edge on ONE great circle) resolves boundary
# queries exactly and refuses interior location with the documented
# degenerate-ring error, unchanged.
@testset "degenerate vertex mass falls back to the wedge bootstrap" begin
    equator = [(0., 0.), (90., 0.), (180., 0.), (-90., 0.)]
    usps = [GO.UnitSpherical.UnitSphereFromGeographic()(p) for p in equator]
    @test GO.UnitSpherical.spherical_exterior_anchor(usps, 4) === nothing
    hemi = GI.Polygon([GI.LinearRing([equator; [equator[1]]])])
    @test GO.intersects(alg, hemi, GI.Point(45., 0.))            # boundary, exact
    @test_throws ArgumentError GO.relate(alg, hemi, GI.Point(10., 45.))
    #-- a ring 0.5° off the equator has a tiny but usable vertex mass: the
    #-- parity path still answers, 89.5° from the anchor
    near_eq = [(Float64(lon), 0.5) for lon in 0:30:360]
    cap = GI.Polygon([GI.LinearRing(near_eq)])
    @test GO.contains(alg, cap, GI.Point(0., 45.))
    @test !GO.intersects(alg, cap, GI.Point(0., -45.))
end
